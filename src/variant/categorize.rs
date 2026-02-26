//! Variant categorization per specification 02.
//!
//! Decision tree: breakend → SV, symbolic STR → RepeatExpansion,
//! symbolic + CNV/LOH → CNV, symbolic → SV, else SmallVariant.

use crate::variant::types::{VariantCategory, VariantType};

/// Result of categorizing one ALT allele.
#[derive(Debug, Clone)]
pub struct CategorizedVariant {
    pub category: VariantCategory,
    pub variant_type: VariantType,
    /// For repeat expansions: the parsed repeat count.
    pub repeat_count: Option<i32>,
}

/// Categorize and determine the type of a variant from its ALT allele and INFO fields.
///
/// `alt` is a single ALT allele string (already filtered for non-informative alleles).
/// `svtype` is the INFO SVTYPE value if present.
/// `trimmed_ref_len` and `trimmed_alt_len` are the lengths after bidirectional trimming
/// (only used for SmallVariant type determination).
pub fn categorize(
    alt: &str,
    svtype: Option<&str>,
    trimmed_ref_len: usize,
    trimmed_alt_len: usize,
) -> CategorizedVariant {
    // 1. Breakend test: ALT contains '[' or ']'
    if alt.contains('[') || alt.contains(']') {
        return CategorizedVariant {
            category: VariantCategory::Sv,
            variant_type: VariantType::Translocation,
            repeat_count: None,
        };
    }

    // 2. Symbolic allele test: starts with '<' and ends with '>'
    if alt.starts_with('<') && alt.ends_with('>') {
        let inner = &alt[1..alt.len() - 1];

        // STR prefix → RepeatExpansion
        if inner.starts_with("STR") {
            let repeat_count = inner
                .strip_prefix("STR")
                .and_then(|s| s.parse::<i32>().ok());
            return CategorizedVariant {
                category: VariantCategory::RepeatExpansion,
                variant_type: VariantType::ShortTandemRepeatVariation,
                repeat_count,
            };
        }

        // Check SVTYPE for CNV/LOH
        if let Some(st) = svtype
            && (st == "CNV" || st == "LOH")
        {
            let vtype = cnv_type(alt);
            return CategorizedVariant {
                category: VariantCategory::Cnv,
                variant_type: vtype,
                repeat_count: None,
            };
        }

        // Default symbolic → SV
        let vtype = sv_type(alt, svtype);
        return CategorizedVariant {
            category: VariantCategory::Sv,
            variant_type: vtype,
            repeat_count: None,
        };
    }

    // 3. Default: SmallVariant
    let variant_type = small_variant_type(trimmed_ref_len, trimmed_alt_len);
    CategorizedVariant {
        category: VariantCategory::SmallVariant,
        variant_type,
        repeat_count: None,
    }
}

/// Determine SmallVariant type from trimmed allele lengths.
fn small_variant_type(ref_len: usize, alt_len: usize) -> VariantType {
    match (ref_len, alt_len) {
        (1, 1) => VariantType::Snv,
        (r, a) if r == a && r > 1 => VariantType::Mnv,
        (r, 0) if r > 0 => VariantType::Deletion,
        (0, a) if a > 0 => VariantType::Insertion,
        _ => VariantType::Delins,
    }
}

/// Determine SV type from ALT allele and SVTYPE.
fn sv_type(alt: &str, svtype: Option<&str>) -> VariantType {
    match svtype {
        Some("DEL") => VariantType::Deletion,
        Some("INS") => VariantType::Insertion,
        Some("DUP") => {
            if alt == "<DUP:TANDEM>" {
                VariantType::TandemDuplication
            } else {
                VariantType::Duplication
            }
        }
        Some("TDUP") => VariantType::TandemDuplication,
        Some("INV") => VariantType::Inversion,
        Some("BND") => VariantType::Translocation,
        Some("CNV") => VariantType::CopyNumberVariation,
        Some("STR") => VariantType::ShortTandemRepeatVariation,
        Some("ALU" | "LINE1" | "SVA") => VariantType::MobileElementInsertion,
        Some("LOH") => VariantType::CopyNumberVariation,
        _ => VariantType::Unknown,
    }
}

/// Determine CNV type from ALT allele.
fn cnv_type(alt: &str) -> VariantType {
    match alt {
        "<DEL>" => VariantType::CopyNumberLoss,
        "<DUP>" => VariantType::CopyNumberGain,
        _ => VariantType::CopyNumberVariation,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn breakend_is_sv() {
        let result = categorize("A[chr2:1000[", None, 0, 0);
        assert_eq!(result.category, VariantCategory::Sv);
        assert_eq!(result.variant_type, VariantType::Translocation);
    }

    #[test]
    fn str_is_repeat_expansion() {
        let result = categorize("<STR25>", None, 0, 0);
        assert_eq!(result.category, VariantCategory::RepeatExpansion);
        assert_eq!(result.variant_type, VariantType::ShortTandemRepeatVariation);
        assert_eq!(result.repeat_count, Some(25));
    }

    #[test]
    fn str_no_count() {
        let result = categorize("<STR>", None, 0, 0);
        assert_eq!(result.category, VariantCategory::RepeatExpansion);
        assert_eq!(result.repeat_count, None);
    }

    #[test]
    fn del_sv() {
        let result = categorize("<DEL>", Some("DEL"), 0, 0);
        assert_eq!(result.category, VariantCategory::Sv);
        assert_eq!(result.variant_type, VariantType::Deletion);
    }

    #[test]
    fn cnv_with_svtype_cnv() {
        let result = categorize("<CNV>", Some("CNV"), 0, 0);
        assert_eq!(result.category, VariantCategory::Cnv);
        assert_eq!(result.variant_type, VariantType::CopyNumberVariation);
    }

    #[test]
    fn cnv_del() {
        let result = categorize("<DEL>", Some("CNV"), 0, 0);
        assert_eq!(result.category, VariantCategory::Cnv);
        assert_eq!(result.variant_type, VariantType::CopyNumberLoss);
    }

    #[test]
    fn snv() {
        let result = categorize("G", None, 1, 1);
        assert_eq!(result.category, VariantCategory::SmallVariant);
        assert_eq!(result.variant_type, VariantType::Snv);
    }

    #[test]
    fn deletion() {
        let result = categorize("A", None, 3, 0);
        assert_eq!(result.category, VariantCategory::SmallVariant);
        assert_eq!(result.variant_type, VariantType::Deletion);
    }

    #[test]
    fn insertion() {
        let result = categorize("ACGT", None, 0, 3);
        assert_eq!(result.category, VariantCategory::SmallVariant);
        assert_eq!(result.variant_type, VariantType::Insertion);
    }

    #[test]
    fn mnv() {
        let result = categorize("TG", None, 2, 2);
        assert_eq!(result.category, VariantCategory::SmallVariant);
        assert_eq!(result.variant_type, VariantType::Mnv);
    }

    #[test]
    fn delins() {
        let result = categorize("TGC", None, 2, 3);
        assert_eq!(result.category, VariantCategory::SmallVariant);
        assert_eq!(result.variant_type, VariantType::Delins);
    }

    #[test]
    fn tandem_duplication() {
        let result = categorize("<DUP:TANDEM>", Some("DUP"), 0, 0);
        assert_eq!(result.category, VariantCategory::Sv);
        assert_eq!(result.variant_type, VariantType::TandemDuplication);
    }

    #[test]
    fn loh_is_cnv() {
        let result = categorize("<LOH>", Some("LOH"), 0, 0);
        assert_eq!(result.category, VariantCategory::Cnv);
        assert_eq!(result.variant_type, VariantType::CopyNumberVariation);
    }
}
