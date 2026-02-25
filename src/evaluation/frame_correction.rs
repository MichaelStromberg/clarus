//! Frame correction for CDS alignment (Levels 5-6 of the evaluation hierarchy).

use super::alignment;
use super::{EvaluationStats, find_amino_acid_edits, find_amino_acid_edits_reverse};
use crate::codon::{self, CodonTable};
use crate::transcript::types::{AminoAcidEdit, CodingRegion};

/// Apply frame correction results to a coding region.
pub(super) fn accept_frame_correction(
    coding: &mut CodingRegion,
    padding: u8,
    offset: u16,
    edits: Vec<AminoAcidEdit>,
    stats: &mut EvaluationStats,
) {
    coding.cds_padding = padding;
    coding.protein_offset = offset;
    if edits.is_empty() {
        stats.num_contained += 1;
    } else {
        coding.amino_acid_edits = Some(edits);
        stats.num_with_amino_acid_edits += 1;
    }
    stats.num_wrong_frame += 1;
}

/// Attempt frame correction by prepending `padding` bases to CDS.
pub(super) fn try_frame_correction(
    coding: &mut CodingRegion,
    cdna: &[u8],
    table: &CodonTable,
    expected: &[u8],
    is_selenoprotein: bool,
    padding: u8,
    stats: &mut EvaluationStats,
) -> bool {
    let adjusted_start = coding.cdna_start as i64 - 1 - padding as i64;
    let cds_end = coding.cdna_end as usize;

    let adjusted_cds = if adjusted_start < 0 {
        // Pad with N's
        let n_count = (-adjusted_start) as usize;
        let mut padded = vec![b'N'; n_count];
        padded.extend_from_slice(&cdna[..cds_end]);
        padded
    } else {
        cdna[adjusted_start as usize..cds_end].to_vec()
    };

    let mut translated = codon::translate(&adjusted_cds, table);

    codon::strip_trailing_x(&mut translated);

    // For non-selenoproteins, truncate after the first stop codon.
    // Characters after * are translation artifacts (codon::translate continues
    // past stop codons). Selenoproteins are excluded because * represents
    // selenocysteine (U) in those. This is safe in try_frame_correction because
    // we're past slippage detection (Level 4).
    if !is_selenoprotein && let Some(stop_pos) = translated.iter().position(|&b| b == b'*') {
        translated.truncate(stop_pos + 1);
    }

    if translated.is_empty() {
        return false;
    }

    // Try alignment with full translated protein
    if try_frame_alignment(
        coding,
        expected,
        &translated,
        is_selenoprotein,
        padding,
        stats,
    ) {
        return true;
    }

    // Offset fallback: skip the first amino acid, which comes from the extra codon
    // created by frame padding. This handles cases where the padding-induced codon
    // doesn't match any position in expected, causing a systematic 1-position shift
    // that the diagonal-only aligner cannot accommodate.
    if translated.len() > 1
        && try_frame_alignment(
            coding,
            expected,
            &translated[1..],
            is_selenoprotein,
            padding,
            stats,
        )
    {
        return true;
    }

    false
}

/// Try forward and reverse alignment of translated against expected for frame correction.
fn try_frame_alignment(
    coding: &mut CodingRegion,
    expected: &[u8],
    translated: &[u8],
    is_selenoprotein: bool,
    padding: u8,
    stats: &mut EvaluationStats,
) -> bool {
    let max_edits = if is_selenoprotein { 10 } else { 2 };

    // Forward alignment: translated (query) within expected (reference)
    if let Some(result) = alignment::semi_global_align(expected, translated) {
        let edits = find_amino_acid_edits(expected, translated, result.start);
        if edits.len() <= max_edits {
            // safe: protein offsets bounded by protein length (< u16::MAX)
            accept_frame_correction(coding, padding, result.start as u16, edits, stats);
            return true;
        }
    }

    // Reverse alignment: expected (query) within translated (reference).
    // Handles the case where prepending bases creates an extra codon, making
    // translated longer than expected.
    if translated.len() > expected.len()
        && let Some(result) = alignment::semi_global_align(translated, expected)
    {
        let edits = find_amino_acid_edits_reverse(expected, translated, result.start);
        if edits.len() <= max_edits {
            // safe: protein offsets bounded by protein length (< u16::MAX)
            accept_frame_correction(coding, padding, result.start as u16, edits, stats);
            return true;
        }
    }

    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frame_correction_reverse_alignment_no_edits() {
        let mut cdna = Vec::new();
        // Extra upstream base at position 0
        cdna.push(b'G');
        // CDS: 44 bases (cdna_start=2, 1-based -> index 1)
        cdna.extend_from_slice(b"GCGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCATAA");

        assert_eq!(cdna.len(), 45);
        let cds = &cdna[1..45]; // 44 bytes
        assert_eq!(cds.len(), 44);
        assert_eq!(44 % 3, 2); // Not a multiple of 3

        let table = CodonTable::standard();

        // After prepending 1 base: adjusted = cdna[0..45] = 45 bytes = 15 codons
        let adjusted = &cdna[0..45];
        let translated_adjusted = codon::translate(adjusted, &table);
        assert_eq!(translated_adjusted.len(), 15);

        // Expected = last 14 AAs of translated_adjusted (skip first codon)
        let expected: Vec<u8> = translated_adjusted[1..].to_vec();
        assert_eq!(expected.len(), 14);

        // Verify forward alignment would fail (translated is longer than expected)
        assert!(translated_adjusted.len() > expected.len());
        assert!(alignment::semi_global_align(&expected, &translated_adjusted).is_none());

        let mut coding = CodingRegion {
            genomic_start: 2,
            genomic_end: 45,
            cdna_start: 2, // 1-based, CDS starts at index 1
            cdna_end: 45,
            protein_id: "NP_TEST.1".to_string(),
            protein_seq: expected,
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let mut stats = EvaluationStats::default();
        let expected_protein = coding.protein_seq.clone();
        let result = try_frame_correction(
            &mut coding,
            &cdna,
            &table,
            &expected_protein,
            false,
            1,
            &mut stats,
        );

        assert!(result);
        assert_eq!(coding.cds_padding, 1);
        assert_eq!(stats.num_wrong_frame, 1);
        assert_eq!(stats.num_contained, 1); // No edits -> counted as contained
        assert_eq!(stats.num_with_amino_acid_edits, 0);
        assert!(coding.amino_acid_edits.is_none());
    }

    #[test]
    fn frame_correction_reverse_alignment_with_edits() {
        let table = CodonTable::standard();

        let mut cdna = Vec::new();
        // 2 extra upstream bases
        cdna.extend_from_slice(b"GG");
        // CDS: 43 bytes
        cdna.extend_from_slice(b"CGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCATAA");

        assert_eq!(cdna.len(), 45);
        assert_eq!(cdna[2..].len(), 43);
        assert_eq!(43 % 3, 1);

        // adjusted = cdna[0..45] = 45 bytes = 15 codons
        let translated_adjusted = codon::translate(&cdna[0..45], &table);
        assert_eq!(translated_adjusted.len(), 15);

        // Expected = last 14 AAs but with 1 edit (change position 5 to 'W')
        let mut expected: Vec<u8> = translated_adjusted[1..].to_vec();
        assert_eq!(expected.len(), 14);
        expected[4] = b'W'; // Introduce 1 mismatch at position 5 (1-based)

        assert!(alignment::semi_global_align(&expected, &translated_adjusted).is_none());

        let mut coding = CodingRegion {
            genomic_start: 3,
            genomic_end: 45,
            cdna_start: 3,
            cdna_end: 45,
            protein_id: "NP_TEST.1".to_string(),
            protein_seq: expected.clone(),
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let mut stats = EvaluationStats::default();
        let result =
            try_frame_correction(&mut coding, &cdna, &table, &expected, false, 2, &mut stats);

        assert!(result);
        assert_eq!(coding.cds_padding, 2);
        assert_eq!(stats.num_wrong_frame, 1);
        assert_eq!(stats.num_with_amino_acid_edits, 1);
        let edits = coding.amino_acid_edits.as_ref().unwrap();
        assert_eq!(edits.len(), 1);
        assert_eq!(edits[0].position, 5);
        assert_eq!(edits[0].amino_acid, b'W');
    }

    #[test]
    fn frame_correction_reverse_seleno_edits() {
        let table = CodonTable::standard();

        // Selenoprotein through reverse path with 5 edits (>2 but <=10)
        let mut adjusted = Vec::new();
        adjusted.extend_from_slice(b"GGC"); // G
        for _ in 0..30 {
            adjusted.extend_from_slice(b"GCA"); // A
        }
        adjusted.extend_from_slice(b"TAA"); // *
        assert_eq!(adjusted.len(), 96);

        let cdna = adjusted.clone();
        assert_eq!(cdna[1..].len(), 95);
        assert_eq!(95 % 3, 2);

        let translated_adjusted = codon::translate(&adjusted, &table);
        assert_eq!(translated_adjusted.len(), 32);

        // Expected = last 31 AAs with 5 edits
        let mut expected: Vec<u8> = translated_adjusted[1..].to_vec();
        assert_eq!(expected.len(), 31);
        for i in 0..5 {
            expected[i * 3] = b'W'; // 5 mismatches spread out
        }

        let mut coding = CodingRegion {
            genomic_start: 2,
            genomic_end: 96,
            cdna_start: 2,
            cdna_end: 96,
            protein_id: "NP_TEST.1".to_string(),
            protein_seq: expected.clone(),
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let mut stats = EvaluationStats::default();
        let result = try_frame_correction(
            &mut coding,
            &cdna,
            &table,
            &expected,
            true, // selenoprotein
            1,
            &mut stats,
        );

        assert!(result);
        assert_eq!(coding.cds_padding, 1);
        assert_eq!(stats.num_wrong_frame, 1);
        assert_eq!(stats.num_with_amino_acid_edits, 1);
        let edits = coding.amino_acid_edits.as_ref().unwrap();
        assert_eq!(edits.len(), 5);
    }

    #[test]
    fn frame_correction_reverse_non_seleno_edit_limit_falls_through() {
        let table = CodonTable::standard();

        // Non-selenoprotein with 3 edits through reverse path -> falls through (returns Ok(false))
        let mut adjusted = Vec::new();
        adjusted.extend_from_slice(b"GGC"); // G
        for _ in 0..30 {
            adjusted.extend_from_slice(b"GCA"); // A
        }
        adjusted.extend_from_slice(b"TAA"); // *
        assert_eq!(adjusted.len(), 96);

        let mut cdna = Vec::new();
        cdna.push(adjusted[0]);
        cdna.extend_from_slice(&adjusted[1..]);
        assert_eq!(cdna.len(), 96);

        let translated_adjusted = codon::translate(&adjusted, &table);

        let mut expected: Vec<u8> = translated_adjusted[1..].to_vec();
        // Introduce 3 edits (exceeds non-seleno limit of 2)
        expected[0] = b'W';
        expected[5] = b'W';
        expected[10] = b'W';

        let mut coding = CodingRegion {
            genomic_start: 2,
            genomic_end: 96,
            cdna_start: 2,
            cdna_end: 96,
            protein_id: "NP_TEST.1".to_string(),
            protein_seq: expected.clone(),
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let mut stats = EvaluationStats::default();
        let result = try_frame_correction(
            &mut coding,
            &cdna,
            &table,
            &expected,
            false, // not selenoprotein
            1,
            &mut stats,
        );

        // Falls through instead of erroring — returns Ok(false)
        assert!(!result);
    }

    #[test]
    fn frame_correction_reverse_alignment_trailing_x() {
        let table = CodonTable::standard();

        let mut cdna = Vec::new();
        cdna.push(b'A'); // upstream base
        cdna.extend_from_slice(b"TG"); // CDS[0..2]: combined with 'A' → "ATG"
        for _ in 0..20 {
            cdna.extend_from_slice(b"GCA"); // A
        }
        cdna.extend_from_slice(b"TAA"); // *
        cdna.extend_from_slice(b"GC"); // incomplete trailing
        assert_eq!(cdna.len(), 68);
        assert_eq!(cdna[1..].len(), 67); // CDS = 67 bytes
        assert_eq!(67 % 3, 1);

        // Verify adjusted translation
        let translated_raw = codon::translate(&cdna[0..68], &table);
        assert_eq!(translated_raw.len(), 23);
        assert_eq!(translated_raw[0], b'M');
        assert_eq!(translated_raw[21], b'*');
        assert_eq!(translated_raw[22], b'X');

        let expected: Vec<u8> = translated_raw[..22].to_vec(); // M + A*20 + *

        let mut coding = CodingRegion {
            genomic_start: 2,
            genomic_end: 68,
            cdna_start: 2, // 1-based → CDS starts at index 1
            cdna_end: 68,
            protein_id: "NP_TEST.1".to_string(),
            protein_seq: expected.clone(),
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let mut stats = EvaluationStats::default();
        let result = try_frame_correction(
            &mut coding,
            &cdna,
            &table,
            &expected,
            false,
            1, // padding=1
            &mut stats,
        );

        assert!(result);
        assert_eq!(coding.cds_padding, 1);
        assert_eq!(stats.num_wrong_frame, 1);
        assert_eq!(stats.num_contained, 1);
    }

    #[test]
    fn frame_correction_trailing_x_stripped_enables_forward() {
        let table = CodonTable::standard();

        let mut expected: Vec<u8> = [b'M']
            .into_iter()
            .chain(std::iter::repeat_n(b'A', 73))
            .chain([b'*'])
            .collect();
        assert_eq!(expected.len(), 75);

        let mut adjusted = Vec::new();
        adjusted.extend_from_slice(b"GTC"); // V
        for _ in 0..42 {
            adjusted.extend_from_slice(b"GCA"); // A
        }
        adjusted.extend_from_slice(b"TAA"); // *
        adjusted.extend_from_slice(b"GC"); // leftover → X
        assert_eq!(adjusted.len(), 134);

        let raw = codon::translate(&adjusted, &table);
        assert_eq!(raw.len(), 45);
        assert_eq!(raw[0], b'V');
        assert_eq!(raw[43], b'*');
        assert_eq!(raw[44], b'X');

        expected[31] = b'V';

        let mut cdna = vec![adjusted[0]];
        cdna.extend_from_slice(&adjusted[1..]);
        assert_eq!(cdna.len(), 134);

        let mut coding = CodingRegion {
            genomic_start: 2,
            genomic_end: 134,
            cdna_start: 2,
            cdna_end: 134,
            protein_id: "NP_TEST.1".to_string(),
            protein_seq: expected.clone(),
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let mut stats = EvaluationStats::default();
        let result =
            try_frame_correction(&mut coding, &cdna, &table, &expected, false, 1, &mut stats);

        assert!(result);
        assert_eq!(coding.cds_padding, 1);
        assert_eq!(stats.num_wrong_frame, 1);
        assert_eq!(stats.num_contained, 1);
        assert_eq!(coding.protein_offset, 31);
    }

    #[test]
    fn frame_correction_post_stop_truncation() {
        let table = CodonTable::standard();

        let expected: Vec<u8> = [b'M']
            .into_iter()
            .chain(std::iter::repeat_n(b'A', 58))
            .chain([b'*'])
            .collect();
        assert_eq!(expected.len(), 60);

        let mut adjusted = Vec::new();
        adjusted.extend_from_slice(b"GGC"); // G
        for _ in 0..19 {
            adjusted.extend_from_slice(b"GCA"); // A
        }
        adjusted.extend_from_slice(b"TAA"); // *
        adjusted.extend_from_slice(b"GGC"); // G (post-stop artifact)
        assert_eq!(adjusted.len(), 66);

        let raw = codon::translate(&adjusted, &table);
        assert_eq!(raw.len(), 22);
        assert_eq!(raw[0], b'G');
        assert_eq!(raw[20], b'*');
        assert_eq!(raw[21], b'G');

        let cdna = adjusted.clone();
        assert_eq!(cdna.len(), 66);

        let mut coding = CodingRegion {
            genomic_start: 2,
            genomic_end: 66,
            cdna_start: 2,
            cdna_end: 66,
            protein_id: "NP_TEST.1".to_string(),
            protein_seq: expected.clone(),
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let mut stats = EvaluationStats::default();
        let result =
            try_frame_correction(&mut coding, &cdna, &table, &expected, false, 1, &mut stats);

        assert!(result);
        assert_eq!(coding.cds_padding, 1);
        assert_eq!(stats.num_wrong_frame, 1);
        assert_eq!(stats.num_with_amino_acid_edits, 1);
        assert_eq!(coding.protein_offset, 39);
    }

    #[test]
    fn frame_correction_offset_fallback() {
        let table = CodonTable::standard();
        let expected = b"VMA*".to_vec();

        let adjusted = b"TGGATGGCATAA".to_vec();
        assert_eq!(adjusted.len(), 12);

        let raw = codon::translate(&adjusted, &table);
        assert_eq!(&raw, b"WMA*");

        let cdna = adjusted;

        let mut coding = CodingRegion {
            genomic_start: 2,
            genomic_end: 12,
            cdna_start: 2,
            cdna_end: 12,
            protein_id: "NP_TEST.1".to_string(),
            protein_seq: expected.clone(),
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let mut stats = EvaluationStats::default();
        let result =
            try_frame_correction(&mut coding, &cdna, &table, &expected, false, 1, &mut stats);

        assert!(result);
        assert_eq!(coding.cds_padding, 1);
        assert_eq!(stats.num_wrong_frame, 1);
        assert_eq!(stats.num_contained, 1);
        assert_eq!(coding.protein_offset, 1);
    }
}
