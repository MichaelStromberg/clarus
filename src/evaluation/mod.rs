//! Transcript evaluation: validates CDS translation against expected protein.
//!
//! Implements the 7-level resolution hierarchy from spec 07.

pub mod alignment;
pub mod selenoprotein;

use crate::cli;
use crate::codon::{self, CodonTable};
use crate::error::Error;
use crate::transcript::types::{
    AminoAcidEdit, CodingRegion, IntermediateTranscript, TranslationalSlip,
};

/// Statistics from transcript evaluation.
#[derive(Debug, Default)]
pub struct EvaluationStats {
    pub num_perfect_protein: u32,
    pub num_non_coding: u32,
    pub num_contained: u32,
    pub num_with_amino_acid_edits: u32,
    pub num_translational_slippage: u32,
    pub num_wrong_frame: u32,
    pub num_unresolvable: u32,
}

/// Evaluate all coding transcripts, setting correction fields on their coding regions.
pub fn evaluate_transcripts(
    transcripts: &mut [IntermediateTranscript],
    chromosomes: &[crate::chromosome::Chromosome],
) -> Result<EvaluationStats, Error> {
    let mut stats = EvaluationStats::default();

    for transcript in transcripts.iter_mut() {
        let coding = match transcript.coding_region.as_mut() {
            Some(cr) => cr,
            None => {
                stats.num_non_coding += 1;
                continue;
            }
        };

        // Determine codon table
        let is_mito = chromosomes
            .get(transcript.chromosome_index)
            .is_some_and(|chr| chr.ensembl_name == "MT");
        let table = if is_mito {
            CodonTable::mitochondrial()
        } else {
            CodonTable::standard()
        };

        let is_seleno = selenoprotein::is_selenoprotein(transcript.gene.hgnc_id);

        match evaluate_single(
            coding,
            &transcript.cdna_seq,
            &table,
            is_seleno,
            is_mito,
            &transcript.id,
            &mut stats,
        ) {
            Ok(()) => {}
            Err(Error::UnresolvableTranscript(_)) => {
                // Print chromosome context
                if let Some(chr) = chromosomes.get(transcript.chromosome_index) {
                    cli::kv(
                        "Chromosome",
                        &format!(
                            "{} (ensembl: {}, refseq: {}, index: {})",
                            chr.ucsc_name,
                            chr.ensembl_name,
                            chr.refseq_accession,
                            transcript.chromosome_index,
                        ),
                    );
                }
                stats.num_unresolvable += 1;
                transcript.coding_region = None;
            }
            Err(other) => return Err(other),
        }
    }

    Ok(stats)
}

/// Evaluate a single coding transcript.
fn evaluate_single(
    coding: &mut CodingRegion,
    cdna: &[u8],
    table: &CodonTable,
    is_selenoprotein: bool,
    _is_mito: bool,
    transcript_id: &str,
    stats: &mut EvaluationStats,
) -> Result<(), Error> {
    // Clone expected protein to avoid borrowing coding while mutating
    let expected = coding.protein_seq.clone();

    // Extract CDS from cDNA
    let cds_start = (coding.cdna_start - 1) as usize;
    let cds_end = coding.cdna_end as usize;
    if cds_start >= cdna.len() || cds_end > cdna.len() {
        return Err(Error::UnresolvableTranscript(format!(
            "{transcript_id}: CDS coordinates [{}, {}] out of bounds for cDNA length {}",
            coding.cdna_start,
            coding.cdna_end,
            cdna.len()
        )));
    }
    let cds = &cdna[cds_start..cds_end];
    let translated = codon::translate(cds, table);

    // Level 1: Perfect match
    if translated == expected {
        stats.num_perfect_protein += 1;
        return Ok(());
    }

    // Level 2-3: Semi-global alignment
    if let Some(result) = alignment::semi_global_align(&expected, &translated) {
        let edits = find_amino_acid_edits(&expected, &translated, result.start);
        let max_edits = if is_selenoprotein { 10 } else { 2 };

        if edits.is_empty() {
            // Level 2: Contained (no edits)
            coding.protein_offset = result.start as u16;
            stats.num_contained += 1;
            return Ok(());
        }

        if edits.len() <= max_edits {
            // Level 3: AA edits
            coding.protein_offset = result.start as u16;
            coding.amino_acid_edits = Some(edits);
            stats.num_with_amino_acid_edits += 1;
            return Ok(());
        }
        // Too many edits — fall through to try lower levels
    }

    // Level 2-3 reverse: translated longer than expected (e.g., extra AA at start)
    if translated.len() > expected.len()
        && let Some(result) = alignment::semi_global_align(&translated, &expected)
    {
        let edits = find_amino_acid_edits_reverse(&expected, &translated, result.start);
        let max_edits = if is_selenoprotein { 10 } else { 2 };

        if edits.is_empty() {
            coding.protein_offset = 0;
            stats.num_contained += 1;
            return Ok(());
        }

        if edits.len() <= max_edits {
            coding.protein_offset = 0;
            coding.amino_acid_edits = Some(edits);
            stats.num_with_amino_acid_edits += 1;
            return Ok(());
        }
        // Too many edits — fall through to try lower levels
    }

    // Level 4: Translational slippage
    if let Some(slip) = detect_slippage(cds, &expected, &translated)
        && slip.length <= 2
    {
        coding.slip = Some(slip);
        stats.num_translational_slippage += 1;
        return Ok(());
    }

    // Level 5: Frame correction +1
    if try_frame_correction(coding, cdna, table, &expected, is_selenoprotein, 1, stats)? {
        return Ok(());
    }

    // Level 6: Frame correction +2
    if try_frame_correction(coding, cdna, table, &expected, is_selenoprotein, 2, stats)? {
        return Ok(());
    }

    // Level 7: Unresolvable — dump diagnostics and return a non-fatal error
    dump_unresolvable_diagnostics(transcript_id, coding, cdna, table);
    Err(Error::UnresolvableTranscript(transcript_id.to_string()))
}

/// Find amino acid differences between expected and translated proteins at given offset.
/// Used when translated (query) is aligned within expected (reference).
fn find_amino_acid_edits(expected: &[u8], translated: &[u8], offset: usize) -> Vec<AminoAcidEdit> {
    let mut edits = Vec::new();
    for (i, &t) in translated.iter().enumerate() {
        let e_pos = offset + i;
        if e_pos < expected.len() && t != expected[e_pos] {
            edits.push(AminoAcidEdit {
                position: (i + 1) as u32,
                amino_acid: expected[e_pos],
            });
        }
    }
    edits
}

/// Find amino acid edits when expected is aligned within translated (reverse alignment).
/// `offset` is the position in `translated` where `expected` starts.
fn find_amino_acid_edits_reverse(
    expected: &[u8],
    translated: &[u8],
    offset: usize,
) -> Vec<AminoAcidEdit> {
    let mut edits = Vec::new();
    for (i, &e) in expected.iter().enumerate() {
        let t_pos = offset + i;
        if t_pos < translated.len() && e != translated[t_pos] {
            edits.push(AminoAcidEdit {
                position: (i + 1) as u32, // 1-based in expected protein
                amino_acid: e,            // expected amino acid
            });
        }
    }
    edits
}

/// Detect translational slippage by checking one-end match.
fn detect_slippage(cds: &[u8], expected: &[u8], translated: &[u8]) -> Option<TranslationalSlip> {
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

/// Dump diagnostic information for an unresolvable transcript (spec Section 8.5).
fn dump_unresolvable_diagnostics(
    transcript_id: &str,
    coding: &CodingRegion,
    cdna: &[u8],
    table: &CodonTable,
) {
    cli::warning(&format!(
        "unresolvable protein mismatch for transcript {transcript_id}"
    ));

    // Coding region coordinates
    cli::kv(
        "Coding region",
        &format!(
            "genomic [{}, {}], cDNA [{}, {}]",
            coding.genomic_start, coding.genomic_end, coding.cdna_start, coding.cdna_end,
        ),
    );

    // Expected protein sequence
    let expected_str = String::from_utf8_lossy(&coding.protein_seq);
    cli::kv("Expected protein", &expected_str);

    // CDS and translations
    let cds_start = (coding.cdna_start - 1) as usize;
    let cds_end = (coding.cdna_end as usize).min(cdna.len());
    if cds_start < cdna.len() {
        let cds = &cdna[cds_start..cds_end];
        let cds_str = String::from_utf8_lossy(cds);
        cli::kv("CDS sequence", &cds_str);

        let translated = codon::translate(cds, table);
        let translated_str = String::from_utf8_lossy(&translated);
        cli::kv("Translated", &translated_str);

        // Frame +1
        if cds_start + 1 < cdna.len() {
            let frame1 = codon::translate(&cdna[cds_start + 1..cds_end], table);
            cli::kv("Frame +1", &String::from_utf8_lossy(&frame1));
        }
        // Frame +2
        if cds_start + 2 < cdna.len() {
            let frame2 = codon::translate(&cdna[cds_start + 2..cds_end], table);
            cli::kv("Frame +2", &String::from_utf8_lossy(&frame2));
        }
    }
}

/// Apply frame correction results to a coding region.
fn accept_frame_correction(
    coding: &mut CodingRegion,
    padding: u8,
    offset: u16,
    edits: Vec<AminoAcidEdit>,
    stats: &mut EvaluationStats,
) {
    coding.cds_padding = padding;
    coding.protein_offset = offset;
    if !edits.is_empty() {
        coding.amino_acid_edits = Some(edits);
        stats.num_with_amino_acid_edits += 1;
    } else {
        stats.num_contained += 1;
    }
    stats.num_wrong_frame += 1;
}

/// Attempt frame correction by prepending `padding` bases to CDS.
fn try_frame_correction(
    coding: &mut CodingRegion,
    cdna: &[u8],
    table: &CodonTable,
    expected: &[u8],
    is_selenoprotein: bool,
    padding: u8,
    stats: &mut EvaluationStats,
) -> Result<bool, Error> {
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

    // Strip trailing X (incomplete codon artifacts from leftover bases)
    while translated.last() == Some(&b'X') {
        translated.pop();
    }

    // For non-selenoproteins, truncate after the first stop codon.
    // Characters after * are translation artifacts (codon::translate continues
    // past stop codons). Selenoproteins are excluded because * represents
    // selenocysteine (U) in those. This is safe in try_frame_correction because
    // we're past slippage detection (Level 4).
    if !is_selenoprotein && let Some(stop_pos) = translated.iter().position(|&b| b == b'*') {
        translated.truncate(stop_pos + 1);
    }

    if translated.is_empty() {
        return Ok(false);
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
        return Ok(true);
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
        return Ok(true);
    }

    Ok(false)
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
    fn perfect_match_evaluation() {
        // ATG GCA TGC TAA = M A C *
        let cdna = b"ATGGCATGCTAA";
        let expected = b"MAC*".to_vec();

        let mut coding = CodingRegion {
            genomic_start: 1,
            genomic_end: 12,
            cdna_start: 1,
            cdna_end: 12,
            protein_id: "NP_001.1".to_string(),
            protein_seq: expected,
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let table = CodonTable::standard();
        let mut stats = EvaluationStats::default();
        evaluate_single(
            &mut coding,
            cdna,
            &table,
            false,
            false,
            "NM_001.1",
            &mut stats,
        )
        .unwrap();
        assert_eq!(stats.num_perfect_protein, 1);
        assert_eq!(coding.cds_padding, 0);
        assert_eq!(coding.protein_offset, 0);
    }

    #[test]
    fn contained_alignment() {
        // Translated protein is shorter — contained in expected
        let cdna = b"ATGGCATGCTAA";
        // Expected protein is longer but contains MAC*
        let expected = b"XXXMAC*XXX".to_vec();

        let mut coding = CodingRegion {
            genomic_start: 1,
            genomic_end: 12,
            cdna_start: 1,
            cdna_end: 12,
            protein_id: "NP_001.1".to_string(),
            protein_seq: expected,
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let table = CodonTable::standard();
        let mut stats = EvaluationStats::default();
        evaluate_single(
            &mut coding,
            cdna,
            &table,
            false,
            false,
            "NM_001.1",
            &mut stats,
        )
        .unwrap();
        assert_eq!(stats.num_contained, 1);
        assert_eq!(coding.protein_offset, 3); // "MAC*" starts at position 3 in "XXXMAC*XXX"
    }

    #[test]
    fn amino_acid_edits_selenoprotein() {
        // Selenoprotein: UGA (TGA) translates to * but expected has U
        // Build a longer CDS so alignment passes threshold (need ~78% identity)
        // Pattern: ATG + 20 matching codons + TGA(=*) + TGA(=*) + 10 matching codons + TAA
        // GCA=A repeated to pad
        let mut cdna_vec = Vec::new();
        cdna_vec.extend_from_slice(b"ATG"); // M
        for _ in 0..20 {
            cdna_vec.extend_from_slice(b"GCA"); // A
        }
        cdna_vec.extend_from_slice(b"TGA"); // * (should be U)
        cdna_vec.extend_from_slice(b"TGA"); // * (should be U)
        for _ in 0..10 {
            cdna_vec.extend_from_slice(b"GCA"); // A
        }
        cdna_vec.extend_from_slice(b"TAA"); // *
        let cdna = cdna_vec.as_slice();

        // Expected protein: M + 20xA + U + U + 10xA + *
        let mut expected = Vec::new();
        expected.push(b'M');
        for _ in 0..20 {
            expected.push(b'A');
        }
        expected.push(b'U');
        expected.push(b'U');
        for _ in 0..10 {
            expected.push(b'A');
        }
        expected.push(b'*');

        let cds_len = cdna.len() as i32;
        let mut coding = CodingRegion {
            genomic_start: 1,
            genomic_end: cds_len,
            cdna_start: 1,
            cdna_end: cds_len,
            protein_id: "NP_001.1".to_string(),
            protein_seq: expected,
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let table = CodonTable::standard();
        let mut stats = EvaluationStats::default();
        evaluate_single(
            &mut coding,
            cdna,
            &table,
            true,
            false,
            "NM_001.1",
            &mut stats,
        )
        .unwrap();
        assert_eq!(stats.num_with_amino_acid_edits, 1);

        let edits = coding.amino_acid_edits.as_ref().unwrap();
        assert_eq!(edits.len(), 2);
        assert_eq!(edits[0].amino_acid, b'U');
        assert_eq!(edits[1].amino_acid, b'U');
    }

    #[test]
    fn non_selenoprotein_too_many_edits_falls_through() {
        // Non-selenoprotein with 3 edits exceeds limit of 2.
        // Falls through alignment to slippage/frame correction, then unresolvable.
        let mut cdna_vec = Vec::new();
        cdna_vec.extend_from_slice(b"ATG"); // M
        for _ in 0..20 {
            cdna_vec.extend_from_slice(b"GCA"); // A
        }
        cdna_vec.extend_from_slice(b"GGC"); // G (should be W)
        cdna_vec.extend_from_slice(b"GGC"); // G (should be W)
        cdna_vec.extend_from_slice(b"GGC"); // G (should be W)
        for _ in 0..5 {
            cdna_vec.extend_from_slice(b"GCA"); // A
        }
        cdna_vec.extend_from_slice(b"TAA"); // *
        let cdna = cdna_vec.as_slice();

        let mut expected = Vec::new();
        expected.push(b'M');
        for _ in 0..20 {
            expected.push(b'A');
        }
        expected.push(b'W');
        expected.push(b'W');
        expected.push(b'W');
        for _ in 0..5 {
            expected.push(b'A');
        }
        expected.push(b'*');

        let cds_len = cdna.len() as i32;
        let mut coding = CodingRegion {
            genomic_start: 1,
            genomic_end: cds_len,
            cdna_start: 1,
            cdna_end: cds_len,
            protein_id: "NP_001.1".to_string(),
            protein_seq: expected,
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let table = CodonTable::standard();
        let mut stats = EvaluationStats::default();
        let result = evaluate_single(
            &mut coding,
            cdna,
            &table,
            false,
            false,
            "NM_001.1",
            &mut stats,
        );
        // Falls through all levels to unresolvable
        assert!(matches!(result, Err(Error::UnresolvableTranscript(_))));
    }

    #[test]
    fn find_edits_basic() {
        let expected = b"MACK*";
        let translated = b"MXCK*";
        let edits = find_amino_acid_edits(expected, translated, 0);
        assert_eq!(edits.len(), 1);
        assert_eq!(edits[0].position, 2);
        assert_eq!(edits[0].amino_acid, b'A');
    }

    #[test]
    fn unresolvable_transcript_is_skipped() {
        use crate::biotype::BioType;
        use crate::chromosome::Chromosome;
        use crate::strand::Strand;
        use crate::transcript::types::{Gene, Source};

        // CDS translates to "MAC*" but expected protein is completely different
        let cdna = b"ATGGCATGCTAA".to_vec();
        let completely_wrong_protein = b"WXYZQRST*".to_vec();

        let coding = CodingRegion {
            genomic_start: 1,
            genomic_end: 12,
            cdna_start: 1,
            cdna_end: 12,
            protein_id: "NP_999.1".to_string(),
            protein_seq: completely_wrong_protein,
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let mut transcripts = vec![IntermediateTranscript {
            chromosome_index: 0,
            start: 1,
            end: 12,
            id: "NM_FAKE.1".to_string(),
            biotype: BioType::MRna,
            source: Source::RefSeq,
            strand: Strand::Forward,
            gene: Gene {
                chromosome_index: 0,
                symbol: "FAKE".to_string(),
                ncbi_gene_id: None,
                ensembl_id: None,
                hgnc_id: None,
                on_reverse_strand: false,
            },
            transcript_regions: Vec::new(),
            coding_region: Some(coding),
            cdna_seq: cdna,
            is_canonical: false,
            is_mane_select: false,
        }];

        let chromosomes = vec![Chromosome {
            ucsc_name: "chr1".to_string(),
            ensembl_name: "1".to_string(),
            refseq_accession: "NC_000001.11".to_string(),
            genbank_accession: "CM000663.2".to_string(),
            length: 248956422,
            ref_index: 0,
        }];

        let stats = evaluate_transcripts(&mut transcripts, &chromosomes).unwrap();
        assert_eq!(stats.num_unresolvable, 1);
        assert!(transcripts[0].coding_region.is_none());
    }

    #[test]
    fn find_amino_acid_edits_reverse_basic() {
        // expected = "ABCDE", translated = "XABCDE" (expected starts at offset 1)
        let expected = b"ABCDE";
        let translated = b"XABCDE";
        let edits = find_amino_acid_edits_reverse(expected, translated, 1);
        assert!(edits.is_empty());

        // With a mismatch: translated = "XABXDE" (C->X at position 3 in expected)
        let translated2 = b"XABXDE";
        let edits2 = find_amino_acid_edits_reverse(expected, translated2, 1);
        assert_eq!(edits2.len(), 1);
        assert_eq!(edits2[0].position, 3); // 1-based position in expected
        assert_eq!(edits2[0].amino_acid, b'C'); // expected amino acid
    }

    #[test]
    fn frame_correction_reverse_alignment_no_edits() {
        // Build a scenario where CDS length is not a multiple of 3.
        // After prepending 1 base (padding=1), the translated protein is 1 AA longer
        // than expected, but expected is contained within it.
        //
        // cDNA: position 0 is the extra upstream base, positions 1..end is the CDS.
        // CDS = 14 bases (not multiple of 3). After prepending 1 base: 15 bases = 5 codons.
        // Translate to 5 AAs. Expected protein is 4 AAs (the last 4 of the 5).
        //
        // Codons from adjusted CDS: GCA GCA ATG GCA TGC = A A M A C
        // Expected protein: AMAC (positions 2-5 of translated)

        // cDNA layout: [G] [CA ATG GCA TGC TAA GCA GCA GCA GCA GCA GCA GCA TAA]
        //               ^extra   ^CDS starts at position 2 (1-based)
        // We need CDS length not multiple of 3, and after prepending 1 base it IS multiple of 3.
        // CDS length = 14 bases. 14 % 3 = 2. Prepend 1 -> 15 bases = 5 codons.
        // But we need the forward alignment to fail (translated longer than expected).

        // Simpler approach: use try_frame_correction directly.
        // expected protein: 30 AAs of 'A' followed by '*'
        // cDNA: [extra_base] [CDS of 92 bases, not multiple of 3]
        // After prepending 1 base: 93 bases = 31 codons -> 31 AAs
        // translated = 31 AAs, expected = 31 AAs... no that's equal.

        // We want translated to be LONGER than expected.
        // CDS = 92 bases (92 % 3 = 2). Translate 92 bases: 30 codons + 2 leftover = 30 AAs.
        // Prepend 1 base: 93 bases = 31 codons = 31 AAs.
        // If expected is 30 AAs, translated (31 AAs) is longer -> forward alignment fails.
        // But expected must match a substring of translated.

        // Let's construct: expected = 30 x 'A' (no stop, like partial proteins)
        // adjusted CDS (93 bases) translates to: X + 30 x 'A' where X is from the
        // first codon formed by [extra_base + first 2 CDS bases].

        // cDNA layout (1-indexed):
        // pos 1: 'C' (extra upstream base)
        // pos 2..93: CDS (92 bases)
        // CDS content: "GC" + "GCA" * 29 + "TAA" = 2 + 87 + 3 = 92 bases
        // Wait, let me just use actual code sequences.

        // cDNA (0-indexed): C GC [GCA x 29] TAA
        // Positions: extra=C, CDS starts at index 1 (1-based: cdna_start=2)
        // CDS = cdna[1..93] = "GC" + "GCA"*29 + "TAA" = 2 + 87 + 3 = 92 bases
        // Direct translation of CDS (92 bases): codons from "GCG CAG CAG CA..." - wait that's wrong.
        //
        // Let me be more careful. CDS = positions 1..92 (0-indexed).
        // After prepending 1 base (position 0 = 'C'):
        // adjusted_cds = cdna[0..93] = "C" + CDS = "CGC" + "GCA"*29 + "TAA"
        // Wait, I need to think about what CDS content gives the right result.

        // Let me just use a concrete small example with enough length for alignment.
        // I'll make a cDNA where:
        // - cdna_start = 2 (1-based), so CDS starts at index 1
        // - CDS = 44 bases (44 % 3 = 2): translates to 14 AAs (44/3 = 14 remainder 2)
        // - Prepend 1 base: 45 bases = 15 codons = 15 AAs
        // - Expected = last 14 AAs of those 15 (i.e., the substring starting at position 1)
        // Build:
        // cdna[0] = 'G' (extra base)
        // cdna[1..45] = CDS of 44 bases
        // adjusted_cds = cdna[0..45] = 45 bases
        // We want adjusted_cds to translate to: X, A, A, A, A, A, A, A, A, A, A, A, A, A, *
        // where X is anything, and expected = "AAAAAAAAAAAAA*" (13 A's + *)
        //
        // 45 bases as codons: GGC GCA GCA GCA GCA GCA GCA GCA GCA GCA GCA GCA GCA GCA TAA
        // That's 15 codons = G, A, A, A, A, A, A, A, A, A, A, A, A, A, *
        // So adjusted_cds = "GGCGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCATAA"
        // CDS = adjusted_cds[1..] = "GCGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCATAA" (44 bytes)
        // Direct translate of CDS: GCG CAG CAG ... not what we want for expected.
        // expected should match the protein from the adjusted (padded) translation minus first AA.

        let mut cdna = Vec::new();
        // Extra upstream base at position 0
        cdna.push(b'G');
        // CDS: 44 bases (cdna_start=2, 1-based -> index 1)
        // Content: "GCGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCATAA"
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

        assert!(result.unwrap());
        assert_eq!(coding.cds_padding, 1);
        assert_eq!(stats.num_wrong_frame, 1);
        assert_eq!(stats.num_contained, 1); // No edits -> counted as contained
        assert_eq!(stats.num_with_amino_acid_edits, 0);
        assert!(coding.amino_acid_edits.is_none());
    }

    #[test]
    fn frame_correction_reverse_alignment_with_edits() {
        let table = CodonTable::standard();

        // Build cDNA where padding=2 creates a translated protein 1 AA longer than expected,
        // but with 1 mismatch.
        // adjusted_cds (after prepend 2): 45 bytes = 15 codons
        // CDS: 43 bytes (43 % 3 = 1), cdna_start=3 (1-based, index 2)
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

        assert!(result.unwrap());
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
        // Build adjusted CDS of 96 bytes = 32 codons: G + 30*A + *
        let mut adjusted = Vec::new();
        adjusted.extend_from_slice(b"GGC"); // G
        for _ in 0..30 {
            adjusted.extend_from_slice(b"GCA"); // A
        }
        adjusted.extend_from_slice(b"TAA"); // *
        assert_eq!(adjusted.len(), 96);

        // cdna = adjusted (extra base at [0], CDS at [1..96])
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

        assert!(result.unwrap());
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
        assert_eq!(result.unwrap(), false);
    }

    #[test]
    fn post_stop_reverse_alignment_contained() {
        // CDS translates to expected + extra residues after stop.
        // translated = "MAC*Y" (5 AAs), expected = "MAC*" (4 AAs)
        // Reverse alignment finds expected within translated → contained.
        let cdna = b"ATGGCATGCTAATATGCA";
        let expected = b"MAC*".to_vec();

        let mut coding = CodingRegion {
            genomic_start: 1,
            genomic_end: 15,
            cdna_start: 1,
            cdna_end: 15, // includes the extra TAT codon
            protein_id: "NP_001.1".to_string(),
            protein_seq: expected,
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let table = CodonTable::standard();
        let mut stats = EvaluationStats::default();

        evaluate_single(
            &mut coding,
            cdna,
            &table,
            false,
            false,
            "NM_TEST.1",
            &mut stats,
        )
        .unwrap();
        assert_eq!(stats.num_contained, 1);
    }

    #[test]
    fn reverse_alignment_level2_extra_at_start() {
        // Translated = "P" + expected (extra P at start from upstream CDS start)
        // This tests the reverse alignment at Level 2-3
        // Build CDS: CCA(=P) + ATG(=M) + 20*GCA(=A) + TAA(=*)
        // = 3 + 3 + 60 + 3 = 69 bytes
        let mut cdna = Vec::new();
        cdna.extend_from_slice(b"CCA"); // P
        cdna.extend_from_slice(b"ATG"); // M
        for _ in 0..20 {
            cdna.extend_from_slice(b"GCA"); // A
        }
        cdna.extend_from_slice(b"TAA"); // *
        assert_eq!(cdna.len(), 69);

        // Translated = P M A*20 *
        // Expected = M A*20 * (without the leading P)
        let table = CodonTable::standard();
        let translated = codon::translate(&cdna, &table);
        assert_eq!(translated.len(), 23); // P + M + 20*A + *
        let expected: Vec<u8> = translated[1..].to_vec(); // M + 20*A + *

        let mut coding = CodingRegion {
            genomic_start: 1,
            genomic_end: 69,
            cdna_start: 1,
            cdna_end: 69,
            protein_id: "NP_TEST.1".to_string(),
            protein_seq: expected,
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let mut stats = EvaluationStats::default();
        evaluate_single(
            &mut coding,
            &cdna,
            &table,
            false,
            false,
            "NM_TEST.1",
            &mut stats,
        )
        .unwrap();
        assert_eq!(stats.num_contained, 1);
        assert_eq!(coding.protein_offset, 0);
        assert!(coding.amino_acid_edits.is_none());
    }

    #[test]
    fn frame_correction_reverse_alignment_trailing_x() {
        // Frame correction produces translation with trailing X after stop.
        // Without trimming, translated (23 AAs) > expected (22 AAs).
        // Reverse alignment finds expected within translated → contained + wrong_frame.
        //
        // Layout: cdna[0] = 'A' (upstream), cdna[1..68] = CDS (67 bytes, 67%3=1)
        // adjusted_cds = cdna[0..68] = 68 bytes: ATG + GCA*20 + TAA + GC
        // translate(adjusted) = M + A*20 + * + X (23 chars)
        // expected = M + A*20 + * (22 chars)

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

        assert!(result.unwrap());
        assert_eq!(coding.cds_padding, 1);
        assert_eq!(stats.num_wrong_frame, 1);
        assert_eq!(stats.num_contained, 1);
    }

    #[test]
    fn frame_correction_trailing_x_stripped_enables_forward() {
        // When the adjusted CDS has leftover bytes, translate produces trailing X.
        // With X, the query extends 1 past the end of reference, failing the alignment.
        // Stripping X allows the forward alignment to succeed.
        //
        // expected = 30-char protein: V + 28*A + *
        // CDS = 44 bytes (not multiple of 3). padding=1 → 45 bytes = 15 codons.
        // translated(45 bytes) = G + V + A*12 + * (15 AAs, no X since 45%3=0)
        //
        // But we want X at the end: CDS = 44 bytes, padding=1 → 45 bytes = 15 codons, no leftover.
        // Try: CDS = 43 bytes (43%3=1), padding=1 → 44 bytes = 14 codons + 2 leftover → 14 AA + X.
        //
        // Better approach: construct CDS where partial CDS translates correctly in frame+2,
        // but with leftover bytes creating trailing X.
        //
        // cDNA: [G] [TC AAG GCA GCA GCA GCA GCA GCA GCA GCA GCA GCA GCA TAA GC]
        //        ^extra     ^CDS starts at 1-based position 2
        // CDS = 44 bytes (44%3=2). padding=1 → 45 bytes = 15 codons → 15 AA.
        // But we want trailing X from leftover. Let's use padding=2:
        // CDS = 43 bytes (43%3=1), cdna_start=3 → padding=1: 44 bytes → 14 codons + 2 leftover → 14 AA + X.
        //
        // Simplest: build expected as a long protein, CDS covers the C-terminal part
        // in shifted frame, with trailing X from leftover bytes.

        let table = CodonTable::standard();

        // Build a cDNA where:
        // - expected is 76 chars (like NM_001098722.2)
        // - CDS on alt contig covers only C-terminal portion
        // - Frame correction (padding=1) produces correct match + trailing X

        // Expected protein: M + 73*A + * + (extra to pad) = long
        let mut expected = Vec::new();
        expected.push(b'M');
        for _ in 0..73 {
            expected.push(b'A');
        }
        expected.push(b'*');
        assert_eq!(expected.len(), 75);

        // Build adjusted CDS (after prepending 1 base) that translates to:
        // V + last 43 AAs of expected + X
        // = V + A*42 + * + X = 45 chars
        // For this we need 44 codons + 2 leftover bytes = 134 bytes
        // First codon: GTC (=V), then 42 codons of GCA (=A), then TAA (=*), then 2 leftover bytes
        let mut adjusted = Vec::new();
        adjusted.extend_from_slice(b"GTC"); // V
        for _ in 0..42 {
            adjusted.extend_from_slice(b"GCA"); // A
        }
        adjusted.extend_from_slice(b"TAA"); // *
        adjusted.extend_from_slice(b"GC"); // leftover → X
        assert_eq!(adjusted.len(), 134);

        // Raw translate: V + A*42 + * + X = 45 chars
        let raw = codon::translate(&adjusted, &table);
        assert_eq!(raw.len(), 45);
        assert_eq!(raw[0], b'V');
        assert_eq!(raw[43], b'*');
        assert_eq!(raw[44], b'X');

        // After X stripping: 44 chars (V + A*42 + *)
        // Forward alignment: semi_global_align(expected=75, query=44)
        // At d=31: expected[31]=A, query[0]=V → mismatch. At d=32: expected[32]=A, query[0]=V → mismatch.
        // Hmm, V doesn't appear in expected[31..]. Let me make expected[31] = V.
        expected[31] = b'V';

        // Now at d=31: expected[31]=V, query[0]=V → match. expected[32..74]=A, query[1..43]=A → match.
        // expected[74]=*, query[43]=* → match. 44 matches. Score = 44*15 = 660. Per base = 15.0 ✓

        // Build cDNA: [upstream_base] + [CDS of 133 bytes]
        // adjusted = cdna[0..134], so CDS = cdna[1..134] = 133 bytes, cdna_start=2 (1-based)
        let mut cdna = vec![adjusted[0]]; // First byte is part of extra base
        cdna.extend_from_slice(&adjusted[1..]);
        assert_eq!(cdna.len(), 134);
        // CDS = cdna[1..134] = 133 bytes. 133 % 3 = 1.
        // padding=1: adjusted = cdna[0..134] = 134 bytes = 44 codons + 2 leftover → 44 AA + X

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

        assert!(result.unwrap());
        assert_eq!(coding.cds_padding, 1);
        assert_eq!(stats.num_wrong_frame, 1);
        assert_eq!(stats.num_contained, 1);
        assert_eq!(coding.protein_offset, 31);
    }

    #[test]
    fn frame_correction_post_stop_truncation() {
        // Frame correction produces translation with a residue after stop codon
        // (e.g., *G). Without truncation, the extra G causes 3 mismatches at
        // every alignment position (exceeding max_edits=2). Truncating at *
        // removes the trailing G, leaving 21 chars with only 1 mismatch (G≠A)
        // at d=39 in expected.
        //
        // Expected: M + A*58 + * = 60 chars
        // Raw translation: G + A*19 + * + G = 22 chars (post-stop G)
        // After truncation: G + A*19 + * = 21 chars
        // Forward at d=39: 1 mismatch (G≠A) → accepted with 1 AA edit

        let table = CodonTable::standard();

        // Expected protein: M + A*58 + * = 60 chars
        let mut expected = Vec::new();
        expected.push(b'M');
        for _ in 0..58 {
            expected.push(b'A');
        }
        expected.push(b'*');
        assert_eq!(expected.len(), 60);

        // Build adjusted CDS: G + A*19 + * + G = 22 codons = 66 bytes
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

        // cDNA = adjusted. cdna_start=2, padding=1 → adjusted_start=0
        // CDS = cdna[1..66] = 65 bytes (65%3=2). With padding: 66 bytes = 22 codons.
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

        // Truncation at * → G + A*19 + * = 21 chars
        // Forward alignment at d=39: 1 mismatch (G≠A) → accepted as AA edit
        assert!(result.unwrap());
        assert_eq!(coding.cds_padding, 1);
        assert_eq!(stats.num_wrong_frame, 1);
        assert_eq!(stats.num_with_amino_acid_edits, 1);
        assert_eq!(coding.protein_offset, 39);
    }

    #[test]
    fn frame_correction_offset_fallback() {
        // The padding-induced extra codon creates a junk first AA (W) that
        // doesn't appear in expected. For short proteins (≤ 4 AAs), 1 mismatch
        // produces per_base score 9.5 < threshold 10.0, so forward alignment
        // of the full translated fails. Reverse is skipped (same length).
        // The offset fallback (translated[1:]) finds a perfect match.
        //
        // Expected: V M A * = 4 chars
        // Full translated: W M A * = 4 chars
        //   Forward d=0: V≠W → 1 mm, score 38, per_base 9.5 < 10.0 → fails
        //   Reverse: same length → skipped
        // Offset translated[1:] = M A * = 3 chars
        //   Forward d=1: perfect match → accepted as contained

        let table = CodonTable::standard();
        let expected = b"VMA*".to_vec();

        // Build adjusted CDS: TGG(W) ATG(M) GCA(A) TAA(*) = 12 bytes
        let adjusted = b"TGGATGGCATAA".to_vec();
        assert_eq!(adjusted.len(), 12);

        let raw = codon::translate(&adjusted, &table);
        assert_eq!(&raw, b"WMA*");

        // cDNA = adjusted. cdna_start=2, padding=1 → adjusted_start=0
        // CDS = cdna[1..12] = 11 bytes (11%3=2). With padding: 12 bytes = 4 codons.
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

        // translated[1:] = M A * matches expected at d=1 → contained
        assert!(result.unwrap());
        assert_eq!(coding.cds_padding, 1);
        assert_eq!(stats.num_wrong_frame, 1);
        assert_eq!(stats.num_contained, 1);
        assert_eq!(coding.protein_offset, 1);
    }
}
