//! Genome assembly identification and reference ID computation.

use std::fmt;

use crate::error::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum GenomeAssembly {
    Unknown = 0,
    GRCh37 = 1,
    GRCh38 = 2,
}

impl GenomeAssembly {
    #[must_use]
    pub fn to_byte(self) -> u8 {
        self as u8
    }
}

impl TryFrom<u8> for GenomeAssembly {
    type Error = Error;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(GenomeAssembly::Unknown),
            1 => Ok(GenomeAssembly::GRCh37),
            2 => Ok(GenomeAssembly::GRCh38),
            _ => Err(Error::Parse(format!(
                "invalid genome assembly byte: {value}"
            ))),
        }
    }
}

impl std::str::FromStr for GenomeAssembly {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "grch37" => Ok(GenomeAssembly::GRCh37),
            "grch38" => Ok(GenomeAssembly::GRCh38),
            "unknown" => Ok(GenomeAssembly::Unknown),
            _ => Err(Error::Parse(format!("unrecognized genome assembly: {s}"))),
        }
    }
}

/// Compute a deterministic reference ID from the assembly and patch level.
///
/// The ID is the first 4 bytes (LE) of the SHA-256 hash of `"<assembly>.p<patch_level>"`.
/// This provides a single uint32 that downstream file types (SA, cache) can embed
/// for quick compatibility checks against the reference file.
#[must_use]
pub fn compute_reference_id(assembly: GenomeAssembly, patch_level: u8) -> u32 {
    use sha2::{Digest, Sha256};
    let input = format!("{assembly}.p{patch_level}");
    let hash = Sha256::digest(input.as_bytes());
    u32::from_le_bytes([hash[0], hash[1], hash[2], hash[3]])
}

impl fmt::Display for GenomeAssembly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GenomeAssembly::Unknown => write!(f, "Unknown"),
            GenomeAssembly::GRCh37 => write!(f, "GRCh37"),
            GenomeAssembly::GRCh38 => write!(f, "GRCh38"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn byte_round_trip() {
        for assembly in [
            GenomeAssembly::Unknown,
            GenomeAssembly::GRCh37,
            GenomeAssembly::GRCh38,
        ] {
            let byte = assembly.to_byte();
            let back = GenomeAssembly::try_from(byte).unwrap();
            assert_eq!(assembly, back);
        }
    }

    #[test]
    fn invalid_byte() {
        assert!(GenomeAssembly::try_from(3).is_err());
    }

    #[test]
    fn compute_reference_id_is_deterministic() {
        let id1 = super::compute_reference_id(GenomeAssembly::GRCh38, 14);
        let id2 = super::compute_reference_id(GenomeAssembly::GRCh38, 14);
        assert_eq!(id1, id2);
        assert_ne!(id1, 0);

        // Different patch level produces different ID
        let id3 = super::compute_reference_id(GenomeAssembly::GRCh38, 13);
        assert_ne!(id1, id3);

        // Different assembly produces different ID
        let id4 = super::compute_reference_id(GenomeAssembly::GRCh37, 14);
        assert_ne!(id1, id4);
    }

    #[test]
    fn parse_from_string() {
        assert_eq!(
            "GRCh37".parse::<GenomeAssembly>().unwrap(),
            GenomeAssembly::GRCh37
        );
        assert_eq!(
            "grch38".parse::<GenomeAssembly>().unwrap(),
            GenomeAssembly::GRCh38
        );
        assert_eq!(
            "GRCH38".parse::<GenomeAssembly>().unwrap(),
            GenomeAssembly::GRCh38
        );
        assert!("hg19".parse::<GenomeAssembly>().is_err());
    }
}
