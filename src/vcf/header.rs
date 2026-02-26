//! VCF header parsing per specification 01.
//!
//! Parses `##fileformat`, `##contig` lines, and the `#CHROM` header line.

use std::io::BufRead;

use crate::error::Error;
use crate::vcf::assembly_inference::normalize_contig_name;

/// Parsed VCF header.
#[derive(Debug)]
pub struct VcfHeader {
    /// VCF version string (e.g., "VCFv4.2").
    pub file_format: String,
    /// Contig definitions: (ensembl_name, length).
    pub contigs: Vec<(String, u32)>,
    /// Sample names from the #CHROM header line.
    pub sample_names: Vec<String>,
}

/// Parse VCF header lines (everything before the first data line).
/// Returns the parsed header and the byte offset of the first data line.
pub fn parse_header(data: &[u8]) -> Result<(VcfHeader, usize), Error> {
    let mut file_format = String::new();
    let mut contigs = Vec::new();
    let mut sample_names = Vec::new();
    let mut found_chrom_line = false;
    let mut offset = 0;

    while offset < data.len() {
        // Find end of current line
        let line_end = memchr::memchr(b'\n', &data[offset..])
            .map(|pos| offset + pos)
            .unwrap_or(data.len());

        let line = &data[offset..line_end];
        // Strip trailing \r for Windows line endings
        let line = if line.last() == Some(&b'\r') {
            &line[..line.len() - 1]
        } else {
            line
        };

        let line_str = std::str::from_utf8(line)
            .map_err(|e| Error::VcfFormat(format!("invalid UTF-8 in header: {e}")))?;

        if let Some(stripped) = line_str.strip_prefix("##fileformat=") {
            file_format = stripped.to_string();
        } else if line_str.starts_with("##contig=<") {
            if let Some((name, length)) = parse_contig_line(line_str) {
                let ensembl_name = normalize_contig_name(&name).to_string();
                contigs.push((ensembl_name, length));
            }
        } else if line_str.starts_with("#CHROM") {
            let fields: Vec<&str> = line_str.split('\t').collect();
            // Sample names start at column index 9
            if fields.len() > 9 {
                sample_names = fields[9..].iter().map(|s| (*s).to_string()).collect();
            }
            found_chrom_line = true;
            // Advance past this line
            offset = if line_end < data.len() {
                line_end + 1
            } else {
                data.len()
            };
            break;
        } else if !line_str.starts_with('#') && !line_str.is_empty() {
            // Non-header, non-empty line — we've passed the header without #CHROM
            return Err(Error::VcfFormat("missing #CHROM header line".to_string()));
        }

        offset = if line_end < data.len() {
            line_end + 1
        } else {
            data.len()
        };
    }

    if file_format.is_empty() {
        return Err(Error::VcfFormat("missing ##fileformat line".to_string()));
    }

    if !found_chrom_line {
        return Err(Error::VcfFormat("missing #CHROM header line".to_string()));
    }

    Ok((
        VcfHeader {
            file_format,
            contigs,
            sample_names,
        },
        offset,
    ))
}

/// Parse VCF header from a streaming `BufRead` source.
/// Returns the parsed header. The reader is positioned at the first data line after return.
pub fn parse_header_from_reader<R: BufRead>(reader: &mut R) -> Result<VcfHeader, Error> {
    let mut file_format = String::new();
    let mut contigs = Vec::new();
    let mut sample_names = Vec::new();
    let mut found_chrom_line = false;
    let mut line_buf = String::new();

    loop {
        line_buf.clear();
        let bytes_read = reader
            .read_line(&mut line_buf)
            .map_err(|e| Error::VcfFormat(format!("error reading header line: {e}")))?;

        if bytes_read == 0 {
            break; // EOF
        }

        let line_str = line_buf.trim_end_matches(['\n', '\r']);

        if let Some(stripped) = line_str.strip_prefix("##fileformat=") {
            file_format = stripped.to_string();
        } else if line_str.starts_with("##contig=<") {
            if let Some((name, length)) = parse_contig_line(line_str) {
                let ensembl_name = normalize_contig_name(&name).to_string();
                contigs.push((ensembl_name, length));
            }
        } else if line_str.starts_with("#CHROM") {
            let fields: Vec<&str> = line_str.split('\t').collect();
            if fields.len() > 9 {
                sample_names = fields[9..].iter().map(|s| (*s).to_string()).collect();
            }
            found_chrom_line = true;
            break;
        } else if !line_str.starts_with('#') && !line_str.is_empty() {
            return Err(Error::VcfFormat("missing #CHROM header line".to_string()));
        }
    }

    if file_format.is_empty() {
        return Err(Error::VcfFormat("missing ##fileformat line".to_string()));
    }

    if !found_chrom_line {
        return Err(Error::VcfFormat("missing #CHROM header line".to_string()));
    }

    Ok(VcfHeader {
        file_format,
        contigs,
        sample_names,
    })
}

/// Parse a `##contig=<ID=chr1,length=248956422>` line.
/// Returns (raw_name, length) or None if parsing fails.
fn parse_contig_line(line: &str) -> Option<(String, u32)> {
    // Strip "##contig=<" prefix and ">" suffix
    let inner = line.strip_prefix("##contig=<")?.strip_suffix('>')?;

    let mut id = None;
    let mut length = None;

    for part in inner.split(',') {
        if let Some(val) = part.strip_prefix("ID=") {
            id = Some(val.to_string());
        } else if let Some(val) = part.strip_prefix("length=") {
            length = val.parse().ok();
        }
    }

    Some((id?, length?))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_minimal_header() {
        let header = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let (h, offset) = parse_header(header).unwrap();
        assert_eq!(h.file_format, "VCFv4.2");
        assert!(h.contigs.is_empty());
        assert!(h.sample_names.is_empty());
        assert_eq!(offset, header.len());
    }

    #[test]
    fn parse_header_with_contigs_and_samples() {
        let header = b"##fileformat=VCFv4.2\n\
            ##contig=<ID=chr1,length=248956422>\n\
            ##contig=<ID=chr2,length=242193529>\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2\n";
        let (h, _) = parse_header(header).unwrap();
        assert_eq!(h.contigs.len(), 2);
        assert_eq!(h.contigs[0], ("1".to_string(), 248_956_422));
        assert_eq!(h.contigs[1], ("2".to_string(), 242_193_529));
        assert_eq!(h.sample_names, vec!["SAMPLE1", "SAMPLE2"]);
    }

    #[test]
    fn missing_fileformat_error() {
        let header = b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        assert!(parse_header(header).is_err());
    }

    #[test]
    fn missing_chrom_line_error() {
        let header = b"##fileformat=VCFv4.2\n";
        assert!(parse_header(header).is_err());
    }

    #[test]
    fn parse_contig_line_valid() {
        let (name, length) = parse_contig_line("##contig=<ID=chr1,length=248956422>").unwrap();
        assert_eq!(name, "chr1");
        assert_eq!(length, 248_956_422);
    }

    #[test]
    fn parse_contig_line_with_extra_fields() {
        let (name, length) =
            parse_contig_line("##contig=<ID=chr1,length=248956422,assembly=GRCh38>").unwrap();
        assert_eq!(name, "chr1");
        assert_eq!(length, 248_956_422);
    }
}
