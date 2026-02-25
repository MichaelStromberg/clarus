//! Spatial binning for cache organization.
//!
//! The genome is partitioned into fixed-width bins of 2^20 (1,048,576) base pairs.

use super::BIN_SHIFT;

/// Compute the bin index for a 1-based genomic position.
///
/// The `u8` return type is safe for all supported genomes: the largest human
/// chromosome (chr1, 248 Mbp) uses 238 bins with `BIN_SHIFT = 20`.
#[must_use]
pub fn bin_index(position: i32) -> u8 {
    debug_assert!(
        position >= 1,
        "bin_index requires a 1-based position, got {position}"
    );
    ((position - 1) >> BIN_SHIFT) as u8
}

/// Compute the number of bins needed for a chromosome of the given length.
///
/// The `u8` return type is safe for all supported genomes: the largest human
/// chromosome (chr1, 248 Mbp) requires 238 bins with `BIN_SHIFT = 20`.
#[must_use]
pub fn num_bins(chromosome_length: u32) -> u8 {
    if chromosome_length == 0 {
        return 0;
    }
    (((chromosome_length - 1) >> BIN_SHIFT) + 1) as u8
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bin_index_first_megabase() {
        assert_eq!(bin_index(1), 0);
        assert_eq!(bin_index(1_048_576), 0);
    }

    #[test]
    fn bin_index_second_megabase() {
        assert_eq!(bin_index(1_048_577), 1);
        assert_eq!(bin_index(2_097_152), 1);
    }

    #[test]
    fn num_bins_chr1() {
        // Chr1 GRCh38: 248,956,422 bp
        assert_eq!(num_bins(248_956_422), 238);
    }

    #[test]
    fn num_bins_chrm() {
        // Mitochondria: 16,569 bp → 1 bin
        assert_eq!(num_bins(16_569), 1);
    }

    #[test]
    fn num_bins_chry() {
        // ChrY GRCh38: 57,227,415 bp
        assert_eq!(num_bins(57_227_415), 55);
    }

    #[test]
    fn num_bins_zero_length() {
        assert_eq!(num_bins(0), 0);
    }
}
