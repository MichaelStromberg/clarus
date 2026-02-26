//! Variant category and type enums per specification 02.

use std::fmt;

use serde::Serialize;

/// Top-level variant category determined from the ALT allele.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantCategory {
    Reference,
    SmallVariant,
    Sv,
    Cnv,
    RepeatExpansion,
}

/// Specific variant type within a category.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum VariantType {
    // SmallVariant types
    #[serde(rename = "SNV")]
    Snv,
    #[serde(rename = "insertion")]
    Insertion,
    #[serde(rename = "deletion")]
    Deletion,
    #[serde(rename = "delins")]
    Delins,
    #[serde(rename = "MNV")]
    Mnv,
    // SV types
    #[serde(rename = "duplication")]
    Duplication,
    #[serde(rename = "tandem_duplication")]
    TandemDuplication,
    #[serde(rename = "inversion")]
    Inversion,
    #[serde(rename = "translocation")]
    Translocation,
    #[serde(rename = "mobile_element_insertion")]
    MobileElementInsertion,
    // CNV types
    #[serde(rename = "copy_number_variation")]
    CopyNumberVariation,
    #[serde(rename = "copy_number_loss")]
    CopyNumberLoss,
    #[serde(rename = "copy_number_gain")]
    CopyNumberGain,
    // Shared
    #[serde(rename = "short_tandem_repeat_variation")]
    ShortTandemRepeatVariation,
    #[serde(rename = "run_of_homozygosity")]
    RunOfHomozygosity,
    #[serde(rename = "unknown")]
    Unknown,
}

impl fmt::Display for VariantType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            VariantType::Snv => "SNV",
            VariantType::Insertion => "insertion",
            VariantType::Deletion => "deletion",
            VariantType::Delins => "delins",
            VariantType::Mnv => "MNV",
            VariantType::Duplication => "duplication",
            VariantType::TandemDuplication => "tandem_duplication",
            VariantType::Inversion => "inversion",
            VariantType::Translocation => "translocation",
            VariantType::MobileElementInsertion => "mobile_element_insertion",
            VariantType::CopyNumberVariation => "copy_number_variation",
            VariantType::CopyNumberLoss => "copy_number_loss",
            VariantType::CopyNumberGain => "copy_number_gain",
            VariantType::ShortTandemRepeatVariation => "short_tandem_repeat_variation",
            VariantType::RunOfHomozygosity => "run_of_homozygosity",
            VariantType::Unknown => "unknown",
        };
        write!(f, "{s}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn variant_type_display() {
        assert_eq!(VariantType::Snv.to_string(), "SNV");
        assert_eq!(VariantType::Insertion.to_string(), "insertion");
        assert_eq!(VariantType::Deletion.to_string(), "deletion");
        assert_eq!(VariantType::Delins.to_string(), "delins");
        assert_eq!(VariantType::Mnv.to_string(), "MNV");
        assert_eq!(
            VariantType::TandemDuplication.to_string(),
            "tandem_duplication"
        );
        assert_eq!(
            VariantType::CopyNumberVariation.to_string(),
            "copy_number_variation"
        );
        assert_eq!(
            VariantType::ShortTandemRepeatVariation.to_string(),
            "short_tandem_repeat_variation"
        );
    }
}
