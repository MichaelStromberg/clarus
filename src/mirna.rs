//! Parser for miRNA region annotations from GFF3 files.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};

use flate2::read::GzDecoder;

use crate::error::Error;

/// A miRNA genomic region defined by start and end coordinates (1-based, inclusive).
#[derive(Debug, Clone, PartialEq)]
pub struct MirnaRegion {
    pub begin: u32,
    pub end: u32,
}

/// Parses miRNA regions from a gzip-compressed GFF3 file.
///
/// Filters for lines where the type column (column 2) is `"miRNA"`, extracts the
/// RefSeq accession (column 0) and coordinates (columns 3-4). Returns regions
/// grouped by RefSeq accession, sorted by begin position within each chromosome.
pub fn parse_mirna_gff3<R: Read>(reader: R) -> Result<HashMap<String, Vec<MirnaRegion>>, Error> {
    let decoder = GzDecoder::new(reader);
    let buf_reader = BufReader::new(decoder);
    let mut result: HashMap<String, Vec<MirnaRegion>> = HashMap::new();

    for (line_num, line_result) in buf_reader.lines().enumerate() {
        let line = line_result?;

        // Skip comment and empty lines
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 9 {
            return Err(Error::Parse(format!(
                "line {}: expected at least 9 columns, got {}",
                line_num + 1,
                columns.len()
            )));
        }

        // Filter: column 2 must be "miRNA"
        if columns[2] != "miRNA" {
            continue;
        }

        let accession = columns[0].to_string();
        let begin: u32 = columns[3].parse().map_err(|e| {
            Error::Parse(format!(
                "line {}: invalid start position '{}': {e}",
                line_num + 1,
                columns[3]
            ))
        })?;
        let end: u32 = columns[4].parse().map_err(|e| {
            Error::Parse(format!(
                "line {}: invalid end position '{}': {e}",
                line_num + 1,
                columns[4]
            ))
        })?;

        result
            .entry(accession)
            .or_default()
            .push(MirnaRegion { begin, end });
    }

    // Sort each chromosome's regions by begin position
    for regions in result.values_mut() {
        regions.sort_by_key(|r| r.begin);
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    fn gzip(data: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(data).unwrap();
        encoder.finish().unwrap()
    }

    #[test]
    fn parse_basic_mirna_lines() {
        let gff = b"\
##gff-version 3
NC_000001.11\tBestRefSeq\tmiRNA\t17369\t17436\t.\t-\t.\tID=rna-NR_106918.1
NC_000001.11\tBestRefSeq\tmiRNA\t30366\t30503\t.\t+\t.\tID=rna-NR_036051.1
";
        let compressed = gzip(gff);
        let result = parse_mirna_gff3(compressed.as_slice()).unwrap();

        assert_eq!(result.len(), 1);
        let regions = result.get("NC_000001.11").unwrap();
        assert_eq!(regions.len(), 2);
        assert_eq!(
            regions[0],
            MirnaRegion {
                begin: 17369,
                end: 17436
            }
        );
        assert_eq!(
            regions[1],
            MirnaRegion {
                begin: 30366,
                end: 30503
            }
        );
    }

    #[test]
    fn skips_comments_and_empty_lines() {
        let gff = b"\
##gff-version 3
# This is a comment

NC_000001.11\tBestRefSeq\tmiRNA\t100\t200\t.\t+\t.\tID=test
";
        let compressed = gzip(gff);
        let result = parse_mirna_gff3(compressed.as_slice()).unwrap();

        assert_eq!(result.len(), 1);
        let regions = result.get("NC_000001.11").unwrap();
        assert_eq!(regions.len(), 1);
    }

    #[test]
    fn filters_non_mirna_types() {
        let gff = b"\
##gff-version 3
NC_000001.11\tBestRefSeq\tgene\t100\t200\t.\t+\t.\tID=gene1
NC_000001.11\tBestRefSeq\tmiRNA\t300\t400\t.\t+\t.\tID=mirna1
NC_000001.11\tBestRefSeq\texon\t300\t400\t.\t+\t.\tID=exon1
";
        let compressed = gzip(gff);
        let result = parse_mirna_gff3(compressed.as_slice()).unwrap();

        assert_eq!(result.len(), 1);
        let regions = result.get("NC_000001.11").unwrap();
        assert_eq!(regions.len(), 1);
        assert_eq!(
            regions[0],
            MirnaRegion {
                begin: 300,
                end: 400
            }
        );
    }

    #[test]
    fn too_few_columns_error() {
        let gff = b"\
##gff-version 3
NC_000001.11\tBestRefSeq\tmiRNA\t100
";
        let compressed = gzip(gff);
        let result = parse_mirna_gff3(compressed.as_slice());

        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("expected at least 9 columns")
        );
    }

    #[test]
    fn multiple_chromosomes() {
        let gff = b"\
##gff-version 3
NC_000002.12\tBestRefSeq\tmiRNA\t500\t600\t.\t+\t.\tID=m2
NC_000001.11\tBestRefSeq\tmiRNA\t200\t300\t.\t-\t.\tID=m1a
NC_000001.11\tBestRefSeq\tmiRNA\t100\t150\t.\t+\t.\tID=m1b
";
        let compressed = gzip(gff);
        let result = parse_mirna_gff3(compressed.as_slice()).unwrap();

        assert_eq!(result.len(), 2);

        // chr1 regions should be sorted by begin
        let chr1 = result.get("NC_000001.11").unwrap();
        assert_eq!(chr1.len(), 2);
        assert_eq!(
            chr1[0],
            MirnaRegion {
                begin: 100,
                end: 150
            }
        );
        assert_eq!(
            chr1[1],
            MirnaRegion {
                begin: 200,
                end: 300
            }
        );

        let chr2 = result.get("NC_000002.12").unwrap();
        assert_eq!(chr2.len(), 1);
        assert_eq!(
            chr2[0],
            MirnaRegion {
                begin: 500,
                end: 600
            }
        );
    }
}
