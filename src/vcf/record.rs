//! VCF data line parsing per specification 01.

use crate::error::Error;
use crate::vcf::info::{self, InfoFields};
use crate::vcf::sample::{self, VcfSample, delimited_field};

/// Non-informative ALT alleles that are silently skipped.
const NON_INFORMATIVE_ALTS: &[&str] = &["*", "<*>", "<M>", "<NON_REF>"];

/// A parsed VCF data record (one line).
#[derive(Debug)]
pub struct VcfRecord {
    pub chromosome: String,
    pub position: i64,
    pub id: Option<String>,
    pub ref_allele: String,
    pub alt_alleles: Vec<String>,
    pub quality: Option<f64>,
    pub filters: Vec<String>,
    pub info: InfoFields,
    pub samples: Vec<VcfSample>,
    /// True if this is a reference record (ALT is "." or "<NON_REF>").
    pub is_reference: bool,
    /// End position: INFO END if present, otherwise POS + len(REF) - 1.
    pub end_position: i64,
}

/// Parse a single VCF data line into a `VcfRecord`.
///
/// `line` should be a tab-delimited VCF data line (not a header line).
/// `num_expected_samples` is the number of samples declared in the header.
pub fn parse_record(line: &str, num_expected_samples: usize) -> Result<VcfRecord, Error> {
    let fields: Vec<&str> = line.split('\t').collect();

    let min_fields = if num_expected_samples > 0 { 9 } else { 8 };
    if fields.len() < min_fields {
        return Err(Error::VcfFormat(format!(
            "expected at least {min_fields} tab-separated columns, got {}",
            fields.len()
        )));
    }

    let chromosome = fields[0].to_string();
    let position: i64 = fields[1]
        .parse()
        .map_err(|e| Error::VcfFormat(format!("invalid POS '{}': {e}", fields[1])))?;

    let id = if fields[2] == "." || fields[2].is_empty() {
        None
    } else {
        Some(fields[2].to_string())
    };

    let ref_allele = fields[3].to_string();
    let alt_field = fields[4];

    // Check for reference record
    let is_reference = alt_field == "." || alt_field == "<NON_REF>";

    // Parse ALT alleles, filtering non-informative ones
    let all_alts: Vec<&str> = alt_field.split(',').collect();
    let alt_alleles: Vec<String> = all_alts
        .iter()
        .filter(|a| !NON_INFORMATIVE_ALTS.contains(a))
        .map(|a| (*a).to_string())
        .collect();

    let quality = if fields[5] == "." || fields[5].is_empty() {
        None
    } else {
        fields[5].parse().ok()
    };

    let filters = if fields[6] == "." || fields[6].is_empty() {
        Vec::new()
    } else {
        fields[6].split(';').map(|s| s.to_string()).collect()
    };

    // Parse INFO
    let info = info::parse_info(fields[7]);

    // Compute end position
    let end_position = info
        .end
        .unwrap_or_else(|| position + ref_allele.len() as i64 - 1);

    // Parse samples
    let samples = if is_reference || fields.len() <= 8 || num_expected_samples == 0 {
        Vec::new()
    } else {
        let format_str = fields[8];
        let sample_strs: Vec<&str> = fields[9..].to_vec();

        if sample_strs.len() != num_expected_samples {
            return Err(Error::VcfFormat(format!(
                "expected {} sample columns, got {}",
                num_expected_samples,
                sample_strs.len()
            )));
        }

        sample::parse_samples(format_str, &sample_strs, alt_alleles.len())
    };

    Ok(VcfRecord {
        chromosome,
        position,
        id,
        ref_allele,
        alt_alleles,
        quality,
        filters,
        info,
        samples,
        is_reference,
        end_position,
    })
}

/// Reusable parser that holds scratch buffers for reduced-allocation record parsing.
///
/// Instead of collecting `line.split('\t')` into a `Vec<&str>` per record (19.6M allocations),
/// this stores byte positions of tab characters in a reusable `Vec<usize>` and accesses
/// fields by offset. Similarly reuses a colon-offset buffer for FORMAT key parsing.
pub struct RecordParser {
    /// Byte positions of tab characters in the current line.
    tabs: Vec<usize>,
    /// Byte positions of colon characters in the current FORMAT field.
    format_colons: Vec<usize>,
}

impl Default for RecordParser {
    fn default() -> Self {
        Self::new()
    }
}

impl RecordParser {
    pub fn new() -> Self {
        Self {
            tabs: Vec::with_capacity(16),
            format_colons: Vec::with_capacity(16),
        }
    }

    /// Parse a VCF data line into a `VcfRecord`, reusing internal scratch buffers.
    pub fn parse(&mut self, line: &str, num_expected_samples: usize) -> Result<VcfRecord, Error> {
        // Find all tab positions (reuses Vec across calls)
        self.tabs.clear();
        for (i, &b) in line.as_bytes().iter().enumerate() {
            if b == b'\t' {
                self.tabs.push(i);
            }
        }

        let num_fields = self.tabs.len() + 1;
        let min_fields = if num_expected_samples > 0 { 9 } else { 8 };
        if num_fields < min_fields {
            return Err(Error::VcfFormat(format!(
                "expected at least {min_fields} tab-separated columns, got {num_fields}"
            )));
        }

        // Extract fields by offset (no Vec<&str> allocation)
        let chrom_field = delimited_field(line, &self.tabs, 0);
        let pos_field = delimited_field(line, &self.tabs, 1);
        let id_field = delimited_field(line, &self.tabs, 2);
        let ref_field = delimited_field(line, &self.tabs, 3);
        let alt_field = delimited_field(line, &self.tabs, 4);
        let qual_field = delimited_field(line, &self.tabs, 5);
        let filter_field = delimited_field(line, &self.tabs, 6);
        let info_field = delimited_field(line, &self.tabs, 7);

        let chromosome = chrom_field.to_string();
        let position: i64 = pos_field
            .parse()
            .map_err(|e| Error::VcfFormat(format!("invalid POS '{pos_field}': {e}")))?;

        let id = if id_field == "." || id_field.is_empty() {
            None
        } else {
            Some(id_field.to_string())
        };

        let ref_allele = ref_field.to_string();

        // Check for reference record
        let is_reference = alt_field == "." || alt_field == "<NON_REF>";

        // Parse ALT alleles, filtering non-informative ones (no intermediate Vec<&str>)
        let alt_alleles: Vec<String> = alt_field
            .split(',')
            .filter(|a| !NON_INFORMATIVE_ALTS.contains(a))
            .map(|a| a.to_string())
            .collect();

        let quality = if qual_field == "." || qual_field.is_empty() {
            None
        } else {
            qual_field.parse().ok()
        };

        let filters = if filter_field == "." || filter_field.is_empty() {
            Vec::new()
        } else {
            filter_field.split(';').map(|s| s.to_string()).collect()
        };

        // Parse INFO
        let info = info::parse_info(info_field);

        // Compute end position
        let end_position = info
            .end
            .unwrap_or_else(|| position + ref_allele.len() as i64 - 1);

        // Parse samples
        let samples = if is_reference || num_fields <= 9 || num_expected_samples == 0 {
            Vec::new()
        } else {
            let format_field = delimited_field(line, &self.tabs, 8);
            let num_sample_fields = num_fields - 9;

            if num_sample_fields != num_expected_samples {
                return Err(Error::VcfFormat(format!(
                    "expected {num_expected_samples} sample columns, got {num_sample_fields}"
                )));
            }

            if format_field == "." || format_field.is_empty() {
                (0..num_expected_samples)
                    .map(|_| VcfSample::default())
                    .collect()
            } else {
                // Find colon positions in FORMAT string (reusable across samples)
                self.format_colons.clear();
                for (i, &b) in format_field.as_bytes().iter().enumerate() {
                    if b == b':' {
                        self.format_colons.push(i);
                    }
                }

                let num_alts = alt_alleles.len();
                let mut samples = Vec::with_capacity(num_expected_samples);
                for s_idx in 0..num_expected_samples {
                    let sample_str = delimited_field(line, &self.tabs, 9 + s_idx);
                    samples.push(sample::parse_sample_from_fields(
                        format_field,
                        &self.format_colons,
                        sample_str,
                        num_alts,
                    ));
                }
                samples
            }
        };

        Ok(VcfRecord {
            chromosome,
            position,
            id,
            ref_allele,
            alt_alleles,
            quality,
            filters,
            info,
            samples,
            is_reference,
            end_position,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_simple_snv() {
        let line = "chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT:DP\t0/1:30";
        let record = parse_record(line, 1).unwrap();
        assert_eq!(record.chromosome, "chr1");
        assert_eq!(record.position, 100);
        assert_eq!(record.id, None);
        assert_eq!(record.ref_allele, "A");
        assert_eq!(record.alt_alleles, vec!["G"]);
        assert!(!record.is_reference);
        assert_eq!(record.end_position, 100);
        assert_eq!(record.samples.len(), 1);
    }

    #[test]
    fn parse_reference_record() {
        let line = "chr1\t100\t.\tA\t.\t.\t.\t.\tGT\t0/0";
        let record = parse_record(line, 1).unwrap();
        assert!(record.is_reference);
        assert!(record.samples.is_empty());
    }

    #[test]
    fn non_informative_alts_filtered() {
        let line = "chr1\t100\t.\tA\tG,*,<NON_REF>\t.\t.\t.";
        let record = parse_record(line, 0).unwrap();
        assert_eq!(record.alt_alleles, vec!["G"]);
    }

    #[test]
    fn multi_allelic() {
        let line = "chr1\t100\t.\tACG\tA,ACGT\t.\t.\t.";
        let record = parse_record(line, 0).unwrap();
        assert_eq!(record.alt_alleles, vec!["A", "ACGT"]);
    }

    #[test]
    fn info_end_position() {
        let line = "chr1\t100\t.\tA\t<DEL>\t.\t.\tEND=500;SVTYPE=DEL";
        let record = parse_record(line, 0).unwrap();
        assert_eq!(record.end_position, 500);
    }

    #[test]
    fn filters_parsed() {
        let line = "chr1\t100\t.\tA\tG\t30\tq10;LowDP\t.";
        let record = parse_record(line, 0).unwrap();
        assert_eq!(record.filters, vec!["q10", "LowDP"]);
    }

    #[test]
    fn vcf_id_parsed() {
        let line = "chr1\t100\trs123\tA\tG\t30\tPASS\t.";
        let record = parse_record(line, 0).unwrap();
        assert_eq!(record.id.as_deref(), Some("rs123"));
    }
}
