//! VCF FORMAT + sample column parsing per specification 01.

/// A parsed sample from one VCF record.
#[derive(Debug, Clone)]
pub struct VcfSample {
    pub genotype: Option<String>,
    pub variant_frequencies: Option<Vec<f64>>,
    pub total_depth: Option<i32>,
    pub genotype_quality: Option<f64>,
    pub copy_number: Option<i32>,
    pub minor_haplotype_copy_number: Option<i32>,
    pub repeat_unit_counts: Option<Vec<i32>>,
    pub allele_depths: Option<Vec<i32>>,
    pub failed_filter: bool,
    pub split_read_counts: Option<Vec<i32>>,
    pub paired_end_read_counts: Option<Vec<i32>>,
    pub is_de_novo: bool,
    pub loss_of_heterozygosity: bool,
    pub somatic_quality: Option<f64>,
    pub is_empty: bool,
    pub genotype_dosage: Option<f64>,
    pub genotype_posteriors: Option<Vec<f64>>,
    pub bin_count: Option<i32>,
}

impl Default for VcfSample {
    fn default() -> Self {
        Self {
            genotype: None,
            variant_frequencies: None,
            total_depth: None,
            genotype_quality: None,
            copy_number: None,
            minor_haplotype_copy_number: None,
            repeat_unit_counts: None,
            allele_depths: None,
            failed_filter: false,
            split_read_counts: None,
            paired_end_read_counts: None,
            is_de_novo: false,
            loss_of_heterozygosity: false,
            somatic_quality: None,
            is_empty: true,
            genotype_dosage: None,
            genotype_posteriors: None,
            bin_count: None,
        }
    }
}

/// Parse all sample columns for a VCF record.
///
/// `format_str` is the FORMAT column (e.g., "GT:DP:AD:GQ").
/// `sample_strs` is a slice of sample column strings.
/// `num_alts` is the number of informative ALT alleles (for VF computation).
pub fn parse_samples(format_str: &str, sample_strs: &[&str], num_alts: usize) -> Vec<VcfSample> {
    if format_str == "." || format_str.is_empty() {
        return sample_strs.iter().map(|_| VcfSample::default()).collect();
    }

    let format_keys: Vec<&str> = format_str.split(':').collect();

    sample_strs
        .iter()
        .map(|s| parse_one_sample(&format_keys, s, num_alts))
        .collect()
}

fn parse_one_sample(format_keys: &[&str], sample_str: &str, num_alts: usize) -> VcfSample {
    let mut sample = VcfSample::default();

    if sample_str == "." || sample_str.is_empty() {
        return sample;
    }

    let mut vf_value: Option<f64> = None;
    let mut values = sample_str.split(':');

    for key in format_keys {
        let value = values.next().unwrap_or(".");
        if !is_missing(value) {
            apply_field(&mut sample, &mut vf_value, key, value);
        }
    }

    finalize_sample(&mut sample, vf_value, num_alts);
    sample
}

/// Parse a single sample using pre-computed format delimiter positions.
/// Avoids `Vec<&str>` allocation for both format keys and sample values.
pub(crate) fn parse_sample_from_fields(
    format_str: &str,
    format_colons: &[usize],
    sample_str: &str,
    num_alts: usize,
) -> VcfSample {
    let mut sample = VcfSample::default();

    if sample_str == "." || sample_str.is_empty() {
        return sample;
    }

    let num_keys = format_colons.len() + 1;
    let mut vf_value: Option<f64> = None;
    let mut values = sample_str.split(':');

    for k in 0..num_keys {
        let key = delimited_field(format_str, format_colons, k);
        let value = values.next().unwrap_or(".");
        if !is_missing(value) {
            apply_field(&mut sample, &mut vf_value, key, value);
        }
    }

    finalize_sample(&mut sample, vf_value, num_alts);
    sample
}

/// Apply a single FORMAT key-value pair to a sample.
#[inline]
fn apply_field(sample: &mut VcfSample, vf_value: &mut Option<f64>, key: &str, value: &str) {
    match key {
        "GT" => {
            sample.genotype = Some(value.to_string());
            sample.is_empty = false;
        }
        "DP" => sample.total_depth = value.parse().ok(),
        "AD" => sample.allele_depths = parse_int_array(value),
        "GQ" => sample.genotype_quality = value.parse().ok(),
        "VF" => *vf_value = value.parse().ok(),
        "CN" => sample.copy_number = value.parse().ok(),
        "MCN" => sample.minor_haplotype_copy_number = value.parse().ok(),
        "REPCN" => {
            sample.repeat_unit_counts = parse_int_array_slash(value);
        }
        "FT" => {
            sample.failed_filter = value != "PASS" && value != ".";
        }
        "SR" => sample.split_read_counts = parse_int_array(value),
        "PR" => sample.paired_end_read_counts = parse_int_array(value),
        "DN" => sample.is_de_novo = value == "DeNovo",
        "SQ" => sample.somatic_quality = value.parse().ok(),
        "DS" => sample.genotype_dosage = value.parse().ok(),
        "GP" => sample.genotype_posteriors = parse_float_array(value),
        "BC" => sample.bin_count = value.parse().ok(),
        "AQ" | "ORIGIN" => {}
        _ => {}
    }
}

/// Post-process a parsed sample: compute variant frequencies and detect LOH.
fn finalize_sample(sample: &mut VcfSample, vf_value: Option<f64>, num_alts: usize) {
    sample.variant_frequencies = compute_variant_frequencies(
        vf_value,
        sample.allele_depths.as_deref(),
        sample.total_depth,
        num_alts,
    );

    if let (Some(mcn), Some(cn), Some(gt)) = (
        sample.minor_haplotype_copy_number,
        sample.copy_number,
        &sample.genotype,
    ) && mcn == 0
        && cn >= 2
        && (gt == "1/2" || gt == "1|2")
    {
        sample.loss_of_heterozygosity = true;
    }
}

/// Get the k-th field from a delimiter-separated string using pre-computed positions.
#[inline]
pub(crate) fn delimited_field<'a>(s: &'a str, delimiters: &[usize], index: usize) -> &'a str {
    let start = if index == 0 {
        0
    } else {
        delimiters[index - 1] + 1
    };
    let end = if index < delimiters.len() {
        delimiters[index]
    } else {
        s.len()
    };
    &s[start..end]
}

fn compute_variant_frequencies(
    vf_value: Option<f64>,
    allele_depths: Option<&[i32]>,
    total_depth: Option<i32>,
    num_alts: usize,
) -> Option<Vec<f64>> {
    // Priority 1: VF field (only for single ALT)
    if let Some(vf) = vf_value
        && num_alts == 1
    {
        return Some(vec![vf]);
    }

    // Priority 2: AD-based calculation
    if let Some(ad) = allele_depths
        && ad.len() > 1
    {
        let total: i64 = if let Some(dp) = total_depth {
            // For mitochondrial variants, DP replaces sum(AD)
            dp as i64
        } else {
            ad.iter().map(|v| *v as i64).sum()
        };

        if total > 0 {
            let freqs: Vec<f64> = ad[1..]
                .iter()
                .map(|&alt_depth| alt_depth as f64 / total as f64)
                .collect();
            return Some(freqs);
        }
    }

    None
}

fn is_missing(value: &str) -> bool {
    value == "." || value.is_empty()
}

fn parse_int_array(value: &str) -> Option<Vec<i32>> {
    let values: Result<Vec<i32>, _> = value.split(',').map(|v| v.parse()).collect();
    values.ok()
}

fn parse_int_array_slash(value: &str) -> Option<Vec<i32>> {
    let values: Result<Vec<i32>, _> = value.split('/').map(|v| v.parse()).collect();
    values.ok()
}

fn parse_float_array(value: &str) -> Option<Vec<f64>> {
    let values: Result<Vec<f64>, _> = value.split(',').map(|v| v.parse()).collect();
    values.ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_gt_dp_ad() {
        let samples = parse_samples("GT:DP:AD", &["0/1:30:10,20"], 1);
        assert_eq!(samples.len(), 1);
        let s = &samples[0];
        assert_eq!(s.genotype.as_deref(), Some("0/1"));
        assert_eq!(s.total_depth, Some(30));
        assert_eq!(s.allele_depths, Some(vec![10, 20]));
        assert!(!s.is_empty);
    }

    #[test]
    fn variant_frequency_from_vf() {
        let samples = parse_samples("GT:VF", &["0/1:0.5"], 1);
        let vf = samples[0].variant_frequencies.as_ref().unwrap();
        assert_eq!(vf.len(), 1);
        assert!((vf[0] - 0.5).abs() < f64::EPSILON);
    }

    #[test]
    fn variant_frequency_from_ad() {
        let samples = parse_samples("GT:AD", &["0/1:10,20"], 1);
        let vf = samples[0].variant_frequencies.as_ref().unwrap();
        assert_eq!(vf.len(), 1);
        assert!((vf[0] - 20.0 / 30.0).abs() < 1e-6);
    }

    #[test]
    fn loh_detection() {
        let samples = parse_samples("GT:CN:MCN", &["1/2:2:0"], 2);
        assert!(samples[0].loss_of_heterozygosity);
    }

    #[test]
    fn no_loh_when_cn_too_low() {
        let samples = parse_samples("GT:CN:MCN", &["1/2:1:0"], 2);
        assert!(!samples[0].loss_of_heterozygosity);
    }

    #[test]
    fn missing_sample() {
        let samples = parse_samples("GT:DP", &["."], 1);
        assert!(samples[0].is_empty);
        assert_eq!(samples[0].genotype, None);
    }

    #[test]
    fn repcn_parsing() {
        let samples = parse_samples("GT:REPCN", &["0/1:20/25"], 1);
        assert_eq!(samples[0].repeat_unit_counts, Some(vec![20, 25]));
    }
}
