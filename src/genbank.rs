//! GenBank flat file parser for CDS coordinates and codon_start.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};

use flate2::read::GzDecoder;

use crate::error::Error;

/// CDS record extracted from a GenBank entry.
#[derive(Debug, Clone)]
pub struct CdsRecord {
    pub cdna_start: i32,
    pub cdna_end: i32,
    pub codon_start: u8,
}

/// Parse a gzip-compressed GenBank flat file.
/// Returns a map from transcript ID to CDS record.
pub fn parse_genbank_gz<R: Read>(reader: R) -> Result<HashMap<String, CdsRecord>, Error> {
    let decoder = GzDecoder::new(reader);
    let buf_reader = BufReader::new(decoder);
    parse_genbank(buf_reader)
}

/// Parse GenBank from a buffered reader.
pub fn parse_genbank<R: BufRead>(reader: R) -> Result<HashMap<String, CdsRecord>, Error> {
    let mut results = HashMap::new();
    let mut lines = reader.lines();
    let mut buffer_line: Option<String> = None;

    loop {
        // Read next line (from buffer or stream), skipping blank lines between records
        let first_line = loop {
            let line = match buffer_line.take() {
                Some(line) => line,
                None => match lines.next() {
                    Some(line) => line?,
                    None => return Ok(results), // EOF
                },
            };
            if !line.trim().is_empty() {
                break line;
            }
        };

        // Must start with LOCUS
        if !first_line.starts_with("LOCUS") {
            return Err(Error::Parse(format!(
                "expected LOCUS, got: '{}'",
                &first_line[..first_line.len().min(20)]
            )));
        }

        // Parse header to find VERSION and FEATURES
        let mut transcript_id: Option<String> = None;

        loop {
            let line = match buffer_line.take() {
                Some(l) => l,
                None => match lines.next() {
                    Some(l) => l?,
                    None => {
                        return Err(Error::Parse("unexpected EOF in GenBank header".to_string()));
                    }
                },
            };

            let tag = if line.len() >= 12 {
                line[..12].trim()
            } else {
                line.trim()
            };

            if tag == "VERSION" {
                let rest = if line.len() > 12 {
                    line[12..].trim()
                } else {
                    ""
                };
                transcript_id = Some(rest.to_string());
            } else if tag == "FEATURES" {
                break;
            }
        }

        let transcript_id = transcript_id
            .ok_or_else(|| Error::Parse("VERSION not found before FEATURES".to_string()))?;

        // Check for predicted transcripts (XM, XR) — skip entire record
        let prefix = transcript_id.split('_').next().unwrap_or("");
        let is_predicted = prefix == "XM" || prefix == "XR";

        // Parse features section
        let mut cds_record: Option<CdsRecord> = None;

        loop {
            let line = match buffer_line.take() {
                Some(l) => l,
                None => match lines.next() {
                    Some(l) => l?,
                    None => break,
                },
            };

            // Check for ORIGIN or record terminator
            if line.starts_with("ORIGIN") || line.starts_with("//") {
                if line.starts_with("//") {
                    break;
                }
                // Skip to terminator
                loop {
                    match lines.next() {
                        Some(Ok(l)) if l.starts_with("//") => break,
                        Some(Ok(_)) => continue,
                        _ => break,
                    }
                }
                break;
            }

            let feature_tag = if line.len() >= 21 {
                line[..21].trim()
            } else {
                line.trim()
            };

            if feature_tag == "CDS" && !is_predicted {
                let content = if line.len() > 21 {
                    line[21..].trim().to_string()
                } else {
                    String::new()
                };

                // Parse CDS coordinates
                let (cdna_start, cdna_end) = parse_cds_range(&content)?;
                let mut codon_start: u8 = 0;

                // Read continuation lines for qualifiers
                for l in lines.by_ref() {
                    let cont_line = l?;

                    // Check if we've left the CDS feature
                    if cont_line.starts_with("ORIGIN") || cont_line.starts_with("//") {
                        buffer_line = Some(cont_line);
                        break;
                    }

                    let cont_tag = if cont_line.len() >= 21 {
                        cont_line[..21].trim()
                    } else {
                        cont_line.trim()
                    };

                    if !cont_tag.is_empty() {
                        // New feature tag — buffer and break
                        buffer_line = Some(cont_line);
                        break;
                    }

                    // Qualifier line
                    let qualifier = if cont_line.len() > 21 {
                        cont_line[21..].trim()
                    } else {
                        ""
                    };

                    if let Some(rest) = qualifier.strip_prefix("/codon_start=") {
                        codon_start = rest.parse().map_err(|e| {
                            Error::Parse(format!("invalid codon_start '{rest}': {e}"))
                        })?;
                    }
                }

                if codon_start == 0 {
                    return Err(Error::Parse(format!(
                        "codon_start not found for CDS in transcript {transcript_id}"
                    )));
                }

                cds_record = Some(CdsRecord {
                    cdna_start,
                    cdna_end,
                    codon_start,
                });
            } else if !feature_tag.is_empty() {
                // Non-CDS feature — skip its continuation lines
                for l in lines.by_ref() {
                    let cont_line = l?;

                    if cont_line.starts_with("ORIGIN") || cont_line.starts_with("//") {
                        buffer_line = Some(cont_line);
                        break;
                    }

                    let cont_tag = if cont_line.len() >= 21 {
                        cont_line[..21].trim()
                    } else {
                        cont_line.trim()
                    };

                    if !cont_tag.is_empty() {
                        buffer_line = Some(cont_line);
                        break;
                    }
                }
            }
        }

        if let Some(record) = cds_record {
            if results.contains_key(&transcript_id) {
                return Err(Error::Parse(format!(
                    "duplicate GenBank transcript ID: {transcript_id}"
                )));
            }
            results.insert(transcript_id, record);
        }
    }
}

/// Parse CDS coordinate range from GenBank.
/// Handles simple ranges "142..1621" and join ranges "join(142..516,609..1621)".
fn parse_cds_range(content: &str) -> Result<(i32, i32), Error> {
    let s = content.trim();

    if let Some(inner) = s.strip_prefix("join(").and_then(|s| s.strip_suffix(')')) {
        // Join range: take first start and last end
        let first_range = inner
            .split(',')
            .next()
            .ok_or_else(|| Error::Parse(format!("empty join in CDS: '{content}'")))?;
        let last_range = inner.split(',').next_back().unwrap();

        let start: i32 = first_range
            .split("..")
            .next()
            .unwrap()
            .trim()
            .parse()
            .map_err(|e| Error::Parse(format!("invalid CDS join start: {e}")))?;
        let end: i32 = last_range
            .split("..")
            .last()
            .unwrap()
            .trim()
            .parse()
            .map_err(|e| Error::Parse(format!("invalid CDS join end: {e}")))?;
        Ok((start, end))
    } else {
        // Simple range
        let parts: Vec<&str> = s.split("..").collect();
        if parts.len() != 2 {
            return Err(Error::Parse(format!("invalid CDS range: '{content}'")));
        }
        let start: i32 = parts[0]
            .trim()
            .parse()
            .map_err(|e| Error::Parse(format!("invalid CDS start: {e}")))?;
        let end: i32 = parts[1]
            .trim()
            .parse()
            .map_err(|e| Error::Parse(format!("invalid CDS end: {e}")))?;
        Ok((start, end))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn make_genbank_record(transcript_id: &str, cds_range: &str, codon_start: u8) -> String {
        format!(
            "\
LOCUS       {transcript_id}    2618 bp    mRNA    linear   PRI 01-JAN-2024
DEFINITION  Test transcript.
VERSION     {transcript_id}
FEATURES             Location/Qualifiers
     source          1..2618
     CDS             {cds_range}
                     /codon_start={codon_start}
                     /product=\"test protein\"
     exon            1..2618
ORIGIN
        1 atgcctcaga
//
"
        )
    }

    #[test]
    fn simple_cds_range() {
        let data = make_genbank_record("NM_001005484.2", "142..1621", 1);
        let reader = Cursor::new(data.as_bytes());
        let results = parse_genbank(std::io::BufReader::new(reader)).unwrap();
        assert_eq!(results.len(), 1);
        let cds = &results["NM_001005484.2"];
        assert_eq!(cds.cdna_start, 142);
        assert_eq!(cds.cdna_end, 1621);
        assert_eq!(cds.codon_start, 1);
    }

    #[test]
    fn join_cds_range() {
        let data = make_genbank_record("NM_000001.1", "join(142..516,609..1621)", 1);
        let reader = Cursor::new(data.as_bytes());
        let results = parse_genbank(std::io::BufReader::new(reader)).unwrap();
        let cds = &results["NM_000001.1"];
        assert_eq!(cds.cdna_start, 142);
        assert_eq!(cds.cdna_end, 1621);
    }

    #[test]
    fn codon_start_extraction() {
        let data = make_genbank_record("NM_000002.1", "100..300", 2);
        let reader = Cursor::new(data.as_bytes());
        let results = parse_genbank(std::io::BufReader::new(reader)).unwrap();
        assert_eq!(results["NM_000002.1"].codon_start, 2);
    }

    #[test]
    fn xm_xr_skipped() {
        let data = format!(
            "{}\n{}",
            make_genbank_record("XM_011541469.3", "100..300", 1),
            make_genbank_record("NM_000001.1", "50..200", 1),
        );
        let reader = Cursor::new(data.as_bytes());
        let results = parse_genbank(std::io::BufReader::new(reader)).unwrap();
        assert_eq!(results.len(), 1);
        assert!(results.contains_key("NM_000001.1"));
        assert!(!results.contains_key("XM_011541469.3"));
    }

    #[test]
    fn nr_transcript_no_cds() {
        let data = "\
LOCUS       NR_046018    2618 bp    mRNA    linear   PRI 01-JAN-2024
VERSION     NR_046018.2
FEATURES             Location/Qualifiers
     source          1..2618
     exon            1..354
ORIGIN
        1 atgcctcaga
//
";
        let reader = Cursor::new(data.as_bytes());
        let results = parse_genbank(std::io::BufReader::new(reader)).unwrap();
        assert!(results.is_empty());
    }
}
