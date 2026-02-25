//! Translational slippage detection (Level 4 of the evaluation hierarchy).

use crate::codon::{self, CodonTable};
use crate::transcript::types::TranslationalSlip;

/// Detect translational slippage by checking one-end match.
pub(super) fn detect_slippage(
    cds: &[u8],
    expected: &[u8],
    translated: &[u8],
) -> Option<TranslationalSlip> {
    if expected.len() <= 10 || translated.len() <= 10 {
        return None;
    }

    // Check begin match: positions 2-11 (1-indexed) → 1..11 (0-indexed)
    let begin_match =
        expected.len() >= 11 && translated.len() >= 11 && expected[1..11] == translated[1..11];

    // Check end match: last 11 characters
    let end_match = expected.len() >= 11
        && translated.len() >= 11
        && expected[expected.len() - 11..] == translated[translated.len() - 11..];

    // Exactly one end must match (XOR)
    if !(begin_match ^ end_match) {
        return None;
    }

    // Count 5' matching amino acids (skip first, count it as 1)
    let mut num_5prime_matches: usize = 1;
    let min_len = expected.len().min(translated.len());
    for i in 1..min_len {
        if expected[i] != translated[i] {
            break;
        }
        num_5prime_matches += 1;
    }

    let num_5prime_bases = num_5prime_matches * 3;
    let standard_table = CodonTable::standard();

    // Re-translate with skip+1 and skip+2
    let skip1 = num_5prime_bases + 1;
    let skip2 = num_5prime_bases + 2;

    let retrans1 = if skip1 < cds.len() {
        codon::translate(&cds[skip1..], &standard_table)
    } else {
        Vec::new()
    };
    let retrans2 = if skip2 < cds.len() {
        codon::translate(&cds[skip2..], &standard_table)
    } else {
        Vec::new()
    };

    // Count 3' matching amino acids
    let count_3prime = |retrans: &[u8]| -> usize {
        let mut count = 0;
        let e_len = expected.len();
        let r_len = retrans.len();
        let min = e_len.min(r_len);
        for i in 1..=min {
            if expected[e_len - i] == retrans[r_len - i] {
                count += 1;
            } else {
                break;
            }
        }
        count
    };

    let matches_3prime_1 = count_3prime(&retrans1);
    let matches_3prime_2 = count_3prime(&retrans2);
    let num_3prime_bases = matches_3prime_1.max(matches_3prime_2) * 3;

    let cds_length = cds.len();
    let slip_length = cds_length.saturating_sub(num_5prime_bases + num_3prime_bases);
    let slip_length = slip_length.min(255) as u8;
    let slip_position = (num_5prime_bases + 1) as i32;

    Some(TranslationalSlip {
        position: slip_position,
        length: slip_length,
    })
}
