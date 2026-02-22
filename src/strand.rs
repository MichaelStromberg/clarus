//! Strand orientation for genomic features.

use std::fmt;

use crate::error::Error;

/// Strand orientation of a genomic feature.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum Strand {
    Forward = 0,
    Reverse = 1,
}

impl Strand {
    /// Parse from GFF3 column 7. "-" is reverse; everything else is forward.
    #[must_use]
    pub fn from_gff3(s: &str) -> Self {
        if s == "-" {
            Self::Reverse
        } else {
            Self::Forward
        }
    }

    #[must_use]
    pub fn is_reverse(self) -> bool {
        self == Self::Reverse
    }

    #[must_use]
    pub fn to_byte(self) -> u8 {
        self as u8
    }
}

impl TryFrom<u8> for Strand {
    type Error = Error;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(Self::Forward),
            1 => Ok(Self::Reverse),
            _ => Err(Error::Parse(format!("invalid strand byte: {value}"))),
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Forward => write!(f, "+"),
            Self::Reverse => write!(f, "-"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_gff3() {
        assert_eq!(Strand::from_gff3("+"), Strand::Forward);
        assert_eq!(Strand::from_gff3("-"), Strand::Reverse);
        assert_eq!(Strand::from_gff3("."), Strand::Forward);
    }

    #[test]
    fn byte_round_trip() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let byte = strand.to_byte();
            let back = Strand::try_from(byte).unwrap();
            assert_eq!(strand, back);
        }
    }

    #[test]
    fn invalid_byte() {
        assert!(Strand::try_from(2).is_err());
    }

    #[test]
    fn is_reverse() {
        assert!(!Strand::Forward.is_reverse());
        assert!(Strand::Reverse.is_reverse());
    }
}
