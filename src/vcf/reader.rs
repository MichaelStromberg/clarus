//! VCF reader orchestrator: decompress → parse header → iterate records.

use std::collections::HashSet;
use std::path::Path;

use crate::error::Error;
use crate::genome_assembly::GenomeAssembly;
use crate::vcf::assembly_inference;
use crate::vcf::bgzf;
use crate::vcf::header::{self, VcfHeader};
use crate::vcf::record::{self, VcfRecord};

/// A fully parsed VCF file.
#[derive(Debug)]
pub struct VcfFile {
    pub header: VcfHeader,
    pub inferred_assembly: GenomeAssembly,
    pub records: Vec<VcfRecord>,
}

/// Read a VCF file: decompress, parse header, infer assembly, parse all records.
/// Validates chromosome contiguity (same chromosome cannot reappear after a different one).
pub fn read_vcf(path: &Path, expected_assembly: GenomeAssembly) -> Result<VcfFile, Error> {
    let data = bgzf::read_vcf_bytes(path)?;
    let (header, data_offset) = header::parse_header(&data)?;

    // Infer assembly from contig lines
    let inferred_assembly = assembly_inference::infer_assembly(&header.contigs);

    // Validate assembly agreement
    if inferred_assembly != GenomeAssembly::Unknown && inferred_assembly != expected_assembly {
        return Err(Error::AssemblyMismatch(format!(
            "VCF contigs indicate {inferred_assembly} but expected {expected_assembly}"
        )));
    }

    let num_samples = header.sample_names.len();
    let data_section = &data[data_offset..];

    // Parse records using memchr for fast line splitting
    let mut records = Vec::new();
    let mut offset = 0;
    let mut seen_chromosomes: HashSet<String> = HashSet::new();
    let mut current_chromosome: Option<String> = None;

    while offset < data_section.len() {
        let line_end = memchr::memchr(b'\n', &data_section[offset..])
            .map(|pos| offset + pos)
            .unwrap_or(data_section.len());

        let line_bytes = &data_section[offset..line_end];
        // Strip trailing \r
        let line_bytes = if line_bytes.last() == Some(&b'\r') {
            &line_bytes[..line_bytes.len() - 1]
        } else {
            line_bytes
        };

        offset = if line_end < data_section.len() {
            line_end + 1
        } else {
            data_section.len()
        };

        // Skip empty lines
        if line_bytes.is_empty() {
            continue;
        }

        let line = std::str::from_utf8(line_bytes)
            .map_err(|e| Error::VcfFormat(format!("invalid UTF-8 in data line: {e}")))?;

        let record = record::parse_record(line, num_samples)?;

        // Skip reference records
        if record.is_reference {
            continue;
        }

        // Skip records with no informative ALT alleles
        if record.alt_alleles.is_empty() {
            continue;
        }

        // Enforce chromosome contiguity
        match &current_chromosome {
            Some(current) if *current == record.chromosome => {
                // Same chromosome, good
            }
            Some(_) | None => {
                // New chromosome — check it hasn't been seen before
                if seen_chromosomes.contains(&record.chromosome) {
                    return Err(Error::VcfFormat(format!(
                        "chromosome {} reappears after other chromosomes (VCF must have contiguous chromosome blocks)",
                        record.chromosome
                    )));
                }
                if let Some(prev) = current_chromosome.take() {
                    seen_chromosomes.insert(prev);
                }
                current_chromosome = Some(record.chromosome.clone());
            }
        }

        records.push(record);
    }

    Ok(VcfFile {
        header,
        inferred_assembly,
        records,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_vcf_bytes(body: &str) -> Vec<u8> {
        let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        format!("{header}{body}").into_bytes()
    }

    #[test]
    fn parse_simple_vcf() {
        let data =
            make_vcf_bytes("chr1\t100\t.\tA\tG\t30\tPASS\t.\nchr1\t200\t.\tC\tT\t40\tPASS\t.\n");
        let (header, offset) = header::parse_header(&data).unwrap();
        let data_section = &data[offset..];
        assert!(!data_section.is_empty());
        assert_eq!(header.sample_names.len(), 0);
    }

    #[test]
    fn reference_records_skipped() {
        let data =
            make_vcf_bytes("chr1\t100\t.\tA\t.\t30\tPASS\t.\nchr1\t200\t.\tC\tT\t40\tPASS\t.\n");
        // Write a temp file to test through read_vcf
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.vcf");
        std::fs::write(&path, &data).unwrap();
        let vcf = read_vcf(&path, GenomeAssembly::Unknown).unwrap();
        assert_eq!(vcf.records.len(), 1);
        assert_eq!(vcf.records[0].position, 200);
    }

    #[test]
    fn chromosome_contiguity_enforced() {
        let data = make_vcf_bytes(
            "chr1\t100\t.\tA\tG\t.\t.\t.\nchr2\t100\t.\tA\tG\t.\t.\t.\nchr1\t200\t.\tC\tT\t.\t.\t.\n",
        );
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.vcf");
        std::fs::write(&path, &data).unwrap();
        let result = read_vcf(&path, GenomeAssembly::Unknown);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("reappears"));
    }

    #[test]
    fn empty_vcf_is_valid() {
        let data = make_vcf_bytes("");
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.vcf");
        std::fs::write(&path, &data).unwrap();
        let vcf = read_vcf(&path, GenomeAssembly::Unknown).unwrap();
        assert!(vcf.records.is_empty());
    }
}
