//! Variant ID (VID) construction per specification 02.
//!
//! Format: `{ensembl_chrom}-{pos}-{ref}-{alt}` for small variants and BND.
//! Format: `{ensembl_chrom}-{pos}-{end}-{ref}-{alt}-{svtype}` for SVs and CNVs.

use crate::variant::types::{VariantCategory, VariantType};

/// Construct a variant ID string.
///
/// Convenience wrapper around [`construct_vid_into`] that returns an owned `String`.
/// Prefer `construct_vid_into` in hot paths to reuse a buffer across calls.
#[allow(clippy::too_many_arguments)]
pub fn construct_vid(
    ensembl_chrom: &str,
    position: i64,
    end_position: i64,
    ref_allele: &str,
    alt_allele: &str,
    category: VariantCategory,
    variant_type: VariantType,
    ref_base: Option<&str>,
) -> String {
    let mut buf = String::with_capacity(64);
    construct_vid_into(
        &mut buf,
        ensembl_chrom,
        position,
        end_position,
        ref_allele,
        alt_allele,
        category,
        variant_type,
        ref_base,
    );
    buf
}

/// Construct a variant ID string directly into a reusable buffer.
///
/// Same logic as `construct_vid` but eliminates intermediate `format!()` allocations
/// by writing directly into the caller's buffer via `push_str` and `itoa`.
#[allow(clippy::too_many_arguments)]
pub fn construct_vid_into(
    buf: &mut String,
    ensembl_chrom: &str,
    position: i64,
    end_position: i64,
    ref_allele: &str,
    alt_allele: &str,
    category: VariantCategory,
    variant_type: VariantType,
    ref_base: Option<&str>,
) {
    buf.clear();

    // Determine display REF: replace "N" with actual reference base
    let display_ref = if ref_allele == "N" {
        ref_base.unwrap_or("N")
    } else {
        ref_allele
    };

    let mut itoa_buf = itoa::Buffer::new();

    // For empty alleles after normalization, add padding base
    if display_ref.is_empty() || alt_allele.is_empty() {
        if let Some(base) = ref_base {
            let final_pos = position - 1;
            match category {
                VariantCategory::SmallVariant
                | VariantCategory::Reference
                | VariantCategory::Sv
                    if variant_type == VariantType::Translocation =>
                {
                    buf.push_str(ensembl_chrom);
                    buf.push('-');
                    buf.push_str(itoa_buf.format(final_pos));
                    buf.push('-');
                    // padded ref
                    if display_ref.is_empty() {
                        buf.push_str(base);
                    } else {
                        buf.push_str(base);
                        buf.push_str(display_ref);
                    }
                    buf.push('-');
                    // padded alt
                    if alt_allele.is_empty() {
                        buf.push_str(base);
                    } else {
                        buf.push_str(base);
                        buf.push_str(alt_allele);
                    }
                }
                VariantCategory::Sv | VariantCategory::Cnv => {
                    let svtype_str = sv_type_for_vid(variant_type);
                    buf.push_str(ensembl_chrom);
                    buf.push('-');
                    buf.push_str(itoa_buf.format(position));
                    buf.push('-');
                    buf.push_str(itoa_buf.format(end_position));
                    buf.push('-');
                    if display_ref.is_empty() {
                        buf.push_str(base);
                    } else {
                        buf.push_str(base);
                        buf.push_str(display_ref);
                    }
                    buf.push('-');
                    if alt_allele.is_empty() {
                        buf.push_str(base);
                    } else {
                        buf.push_str(base);
                        buf.push_str(alt_allele);
                    }
                    buf.push('-');
                    buf.push_str(svtype_str);
                }
                VariantCategory::RepeatExpansion => {
                    buf.push_str(ensembl_chrom);
                    buf.push('-');
                    buf.push_str(itoa_buf.format(position));
                    buf.push('-');
                    buf.push_str(itoa_buf.format(end_position));
                    buf.push('-');
                    if display_ref.is_empty() {
                        buf.push_str(base);
                    } else {
                        buf.push_str(base);
                        buf.push_str(display_ref);
                    }
                    buf.push('-');
                    if alt_allele.is_empty() {
                        buf.push_str(base);
                    } else {
                        buf.push_str(base);
                        buf.push_str(alt_allele);
                    }
                    buf.push_str("-STR");
                }
                _ => {
                    // SmallVariant/Reference with padding
                    buf.push_str(ensembl_chrom);
                    buf.push('-');
                    buf.push_str(itoa_buf.format(final_pos));
                    buf.push('-');
                    if display_ref.is_empty() {
                        buf.push_str(base);
                    } else {
                        buf.push_str(base);
                        buf.push_str(display_ref);
                    }
                    buf.push('-');
                    if alt_allele.is_empty() {
                        buf.push_str(base);
                    } else {
                        buf.push_str(base);
                        buf.push_str(alt_allele);
                    }
                }
            }
        } else {
            write_vid_body(
                buf,
                ensembl_chrom,
                position,
                end_position,
                display_ref,
                alt_allele,
                category,
                variant_type,
                &mut itoa_buf,
            );
        }
    } else {
        write_vid_body(
            buf,
            ensembl_chrom,
            position,
            end_position,
            display_ref,
            alt_allele,
            category,
            variant_type,
            &mut itoa_buf,
        );
    }
}

/// Write the standard (non-padded) VID format into the buffer.
#[allow(clippy::too_many_arguments)]
fn write_vid_body(
    buf: &mut String,
    ensembl_chrom: &str,
    position: i64,
    end_position: i64,
    final_ref: &str,
    final_alt: &str,
    category: VariantCategory,
    variant_type: VariantType,
    itoa_buf: &mut itoa::Buffer,
) {
    match category {
        VariantCategory::SmallVariant | VariantCategory::Reference => {
            buf.push_str(ensembl_chrom);
            buf.push('-');
            buf.push_str(itoa_buf.format(position));
            buf.push('-');
            buf.push_str(final_ref);
            buf.push('-');
            buf.push_str(final_alt);
        }
        VariantCategory::Sv if variant_type == VariantType::Translocation => {
            buf.push_str(ensembl_chrom);
            buf.push('-');
            buf.push_str(itoa_buf.format(position));
            buf.push('-');
            buf.push_str(final_ref);
            buf.push('-');
            buf.push_str(final_alt);
        }
        VariantCategory::Sv | VariantCategory::Cnv => {
            let svtype_str = sv_type_for_vid(variant_type);
            buf.push_str(ensembl_chrom);
            buf.push('-');
            buf.push_str(itoa_buf.format(position));
            buf.push('-');
            buf.push_str(itoa_buf.format(end_position));
            buf.push('-');
            buf.push_str(final_ref);
            buf.push('-');
            buf.push_str(final_alt);
            buf.push('-');
            buf.push_str(svtype_str);
        }
        VariantCategory::RepeatExpansion => {
            buf.push_str(ensembl_chrom);
            buf.push('-');
            buf.push_str(itoa_buf.format(position));
            buf.push('-');
            buf.push_str(itoa_buf.format(end_position));
            buf.push('-');
            buf.push_str(final_ref);
            buf.push('-');
            buf.push_str(final_alt);
            buf.push_str("-STR");
        }
    }
}

fn sv_type_for_vid(variant_type: VariantType) -> &'static str {
    match variant_type {
        VariantType::Deletion => "DEL",
        VariantType::Insertion => "INS",
        VariantType::Duplication | VariantType::TandemDuplication => "DUP",
        VariantType::Inversion => "INV",
        VariantType::CopyNumberVariation
        | VariantType::CopyNumberLoss
        | VariantType::CopyNumberGain => "CNV",
        VariantType::RunOfHomozygosity => "LOH",
        _ => "SV",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snv_vid() {
        let vid = construct_vid(
            "1",
            100,
            100,
            "A",
            "G",
            VariantCategory::SmallVariant,
            VariantType::Snv,
            Some("A"),
        );
        assert_eq!(vid, "1-100-A-G");
    }

    #[test]
    fn insertion_vid_with_padding() {
        let vid = construct_vid(
            "1",
            101,
            103,
            "",
            "CGT",
            VariantCategory::SmallVariant,
            VariantType::Insertion,
            Some("A"),
        );
        // padding base A added before empty ref, position shifted to 100
        assert_eq!(vid, "1-100-A-ACGT");
    }

    #[test]
    fn sv_vid() {
        let vid = construct_vid(
            "1",
            100,
            500,
            "A",
            "<DEL>",
            VariantCategory::Sv,
            VariantType::Deletion,
            Some("A"),
        );
        assert_eq!(vid, "1-100-500-A-<DEL>-DEL");
    }

    #[test]
    fn translocation_vid() {
        let vid = construct_vid(
            "1",
            100,
            100,
            "A",
            "A[chr2:200[",
            VariantCategory::Sv,
            VariantType::Translocation,
            Some("A"),
        );
        assert_eq!(vid, "1-100-A-A[chr2:200[");
    }

    #[test]
    fn repeat_expansion_vid() {
        let vid = construct_vid(
            "1",
            100,
            200,
            "A",
            "<STR25>",
            VariantCategory::RepeatExpansion,
            VariantType::ShortTandemRepeatVariation,
            Some("A"),
        );
        assert_eq!(vid, "1-100-200-A-<STR25>-STR");
    }
}
