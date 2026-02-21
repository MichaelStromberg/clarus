//! Parser for FASTA sequence files.

use std::io::{BufRead, BufReader, Read};

use flate2::read::GzDecoder;

use crate::error::Error;

/// Reads gzip-compressed FASTA and yields (accession, sequence) pairs.
///
/// The accession is extracted from the header line: the first whitespace-delimited
/// token after `>`. Sequence bases are uppercased.
pub fn parse_fasta_gz<R: Read>(reader: R) -> Result<Vec<(String, Vec<u8>)>, Error> {
    let decoder = GzDecoder::new(reader);
    let buf_reader = BufReader::new(decoder);
    parse_fasta(buf_reader)
}

/// Reads FASTA from a buffered reader and yields (accession, sequence) pairs.
fn parse_fasta<R: BufRead>(reader: R) -> Result<Vec<(String, Vec<u8>)>, Error> {
    let mut results: Vec<(String, Vec<u8>)> = Vec::new();
    let mut current_accession: Option<String> = None;
    let mut current_sequence: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            // Finish previous record
            if let Some(acc) = current_accession.take() {
                results.push((acc, current_sequence));
                current_sequence = Vec::new();
            }
            // Extract accession from header using the three patterns from the spec
            let accession = extract_accession(&line)?;
            current_accession = Some(accession);
        } else if current_accession.is_some() {
            // Append sequence data, uppercased
            let trimmed = line.trim();
            let start = current_sequence.len();
            current_sequence.extend_from_slice(trimmed.as_bytes());
            current_sequence[start..].make_ascii_uppercase();
        }
    }

    // Don't forget the last record
    if let Some(acc) = current_accession {
        results.push((acc, current_sequence));
    }

    Ok(results)
}

/// Extracts the accession from a FASTA header line.
///
/// Tries three patterns in order (per spec Section 4.1):
/// 1. `>gi|<digits>|ref|<accession>|` — legacy NCBI GI format
/// 2. `>ref|<accession>|` — modern NCBI RefSeq format
/// 3. `>(\S+)` — generic first token
fn extract_accession(header: &str) -> Result<String, Error> {
    let header = header.trim_start_matches('>');

    // Pattern 1: >gi|<digits>|ref|<accession>|
    if header.starts_with("gi|") {
        let parts: Vec<&str> = header.split('|').collect();
        if parts.len() >= 4 && parts[2] == "ref" {
            let acc = parts[3].trim();
            if !acc.is_empty() {
                return Ok(acc.to_string());
            }
        }
    }

    // Pattern 2: >ref|<accession>|
    if header.starts_with("ref|") {
        let parts: Vec<&str> = header.split('|').collect();
        if parts.len() >= 2 {
            let acc = parts[1].trim();
            if !acc.is_empty() {
                return Ok(acc.to_string());
            }
        }
    }

    // Pattern 3: first whitespace-delimited token
    let first_token = header.split_whitespace().next().unwrap_or("");
    if first_token.is_empty() {
        return Err(Error::Parse(format!("empty FASTA header: >{header}")));
    }
    Ok(first_token.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    fn make_gz(content: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(content).unwrap();
        encoder.finish().unwrap()
    }

    #[test]
    fn parse_single_sequence() {
        let fasta = b">NC_000001.11 Homo sapiens chromosome 1\nACGTacgt\nNNNN\n";
        let gz = make_gz(fasta);
        let results = parse_fasta_gz(std::io::Cursor::new(gz)).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].0, "NC_000001.11");
        assert_eq!(results[0].1, b"ACGTACGTNNNN");
    }

    #[test]
    fn parse_multiple_sequences() {
        let fasta = b">chr1\nACGT\n>chr2\nTTTT\nAAAA\n>chr3\nGGG\n";
        let gz = make_gz(fasta);
        let results = parse_fasta_gz(std::io::Cursor::new(gz)).unwrap();
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].0, "chr1");
        assert_eq!(results[0].1, b"ACGT");
        assert_eq!(results[1].0, "chr2");
        assert_eq!(results[1].1, b"TTTTAAAA");
        assert_eq!(results[2].0, "chr3");
        assert_eq!(results[2].1, b"GGG");
    }

    #[test]
    fn uppercase_bases() {
        let fasta = b">seq1\nacgtACGTnN\n";
        let gz = make_gz(fasta);
        let results = parse_fasta_gz(std::io::Cursor::new(gz)).unwrap();
        assert_eq!(results[0].1, b"ACGTACGTNN");
    }

    #[test]
    fn extract_gi_format() {
        let acc = extract_accession("gi|224589823|ref|NC_000024.9|").unwrap();
        assert_eq!(acc, "NC_000024.9");
    }

    #[test]
    fn extract_ref_format() {
        let acc = extract_accession("ref|NC_000013.11| Homo sapiens chromosome 13").unwrap();
        assert_eq!(acc, "NC_000013.11");
    }

    #[test]
    fn extract_generic_format() {
        let acc = extract_accession(
            "NC_000001.11 Homo sapiens chromosome 1, GRCh38.p14 Primary Assembly",
        )
        .unwrap();
        assert_eq!(acc, "NC_000001.11");
    }
}
