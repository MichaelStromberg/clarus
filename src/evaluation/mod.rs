//! Transcript evaluation: validates CDS translation against expected protein.
//!
//! Implements the 7-level resolution hierarchy from spec 07.

pub mod alignment;
mod diagnostics;
mod frame_correction;
pub mod selenoprotein;
mod slippage;

use crate::cli;
use crate::codon::{self, CodonTable};
use crate::error::Error;
use crate::transcript::types::{AminoAcidEdit, CodingRegion, IntermediateTranscript};

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
        let Some(coding) = transcript.coding_region.as_mut() else {
            stats.num_non_coding += 1;
            continue;
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
            &transcript.id,
            &mut stats,
        ) {
            Ok(()) => {}
            Err(Error::UnresolvableTranscript(ref msg)) => {
                cli::warning(&format!("unresolvable transcript: {msg}"));
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
                // Dump detailed diagnostics if CDS coordinates are valid
                if let Some(ref cr) = transcript.coding_region {
                    let cds_start = (cr.cdna_start - 1) as usize;
                    let cds_end = cr.cdna_end as usize;
                    if cds_start < transcript.cdna_seq.len() && cds_end <= transcript.cdna_seq.len()
                    {
                        diagnostics::dump_unresolvable_diagnostics(
                            cr,
                            &transcript.cdna_seq,
                            &table,
                        );
                    }
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
    let mut translated = codon::translate(cds, table);

    codon::strip_trailing_x(&mut translated);

    // Level 1: Perfect match
    if translated == expected {
        stats.num_perfect_protein += 1;
        return Ok(());
    }

    let max_edits = if is_selenoprotein { 10 } else { 2 };

    // Level 2-3: Semi-global alignment (forward)
    if let Some(result) = alignment::semi_global_align(&expected, &translated) {
        let edits = find_amino_acid_edits(&expected, &translated, result.start);
        // safe: protein offsets bounded by protein length (< u16::MAX)
        if try_accept_alignment(coding, edits, result.start as u16, max_edits, stats) {
            return Ok(());
        }
    }

    // Level 2-3 reverse: translated longer than expected (e.g., extra AA at start)
    if translated.len() > expected.len()
        && let Some(result) = alignment::semi_global_align(&translated, &expected)
    {
        let edits = find_amino_acid_edits_reverse(&expected, &translated, result.start);
        if try_accept_alignment(coding, edits, 0, max_edits, stats) {
            return Ok(());
        }
    }

    // Level 4: Translational slippage
    if let Some(slip) = slippage::detect_slippage(cds, &expected, &translated)
        && slip.length <= 2
    {
        coding.slip = Some(slip);
        stats.num_translational_slippage += 1;
        return Ok(());
    }

    // Level 5: Frame correction +1
    if frame_correction::try_frame_correction(
        coding,
        cdna,
        table,
        &expected,
        is_selenoprotein,
        1,
        stats,
    ) {
        return Ok(());
    }

    // Level 6: Frame correction +2
    if frame_correction::try_frame_correction(
        coding,
        cdna,
        table,
        &expected,
        is_selenoprotein,
        2,
        stats,
    ) {
        return Ok(());
    }

    // Level 7: Unresolvable
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

/// Try to accept an alignment result at Level 2-3.
/// Returns `true` if accepted (contained or AA edits within limit).
fn try_accept_alignment(
    coding: &mut CodingRegion,
    edits: Vec<AminoAcidEdit>,
    protein_offset: u16,
    max_edits: usize,
    stats: &mut EvaluationStats,
) -> bool {
    if edits.len() > max_edits {
        return false;
    }
    coding.protein_offset = protein_offset;
    if edits.is_empty() {
        stats.num_contained += 1;
    } else {
        coding.amino_acid_edits = Some(edits);
        stats.num_with_amino_acid_edits += 1;
    }
    true
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
        evaluate_single(&mut coding, cdna, &table, false, "NM_001.1", &mut stats).unwrap();
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
        evaluate_single(&mut coding, cdna, &table, false, "NM_001.1", &mut stats).unwrap();
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
        let expected: Vec<u8> = [b'M']
            .into_iter()
            .chain(std::iter::repeat_n(b'A', 20))
            .chain([b'U', b'U'])
            .chain(std::iter::repeat_n(b'A', 10))
            .chain([b'*'])
            .collect();

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
        evaluate_single(&mut coding, cdna, &table, true, "NM_001.1", &mut stats).unwrap();
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

        let expected: Vec<u8> = [b'M']
            .into_iter()
            .chain(std::iter::repeat_n(b'A', 20))
            .chain([b'W', b'W', b'W'])
            .chain(std::iter::repeat_n(b'A', 5))
            .chain([b'*'])
            .collect();

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
        let result = evaluate_single(&mut coding, cdna, &table, false, "NM_001.1", &mut stats);
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
        use std::sync::Arc;

        use crate::biotype::BioType;
        use crate::chromosome::Chromosome;
        use crate::strand::Strand;
        use crate::transcript::types::{Designations, Gene, Source};

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
            gene: Arc::new(Gene {
                chromosome_index: 0,
                symbol: "FAKE".to_string(),
                ncbi_gene_id: None,
                ensembl_id: None,
                hgnc_id: None,
                on_reverse_strand: false,
            }),
            transcript_regions: Vec::new(),
            coding_region: Some(coding),
            cdna_seq: cdna,
            designations: Designations::default(),
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

        evaluate_single(&mut coding, cdna, &table, false, "NM_TEST.1", &mut stats).unwrap();
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
        evaluate_single(&mut coding, &cdna, &table, false, "NM_TEST.1", &mut stats).unwrap();
        assert_eq!(stats.num_contained, 1);
        assert_eq!(coding.protein_offset, 0);
        assert!(coding.amino_acid_edits.is_none());
    }

    #[test]
    fn trailing_x_stripped_resolves_as_contained() {
        // CDS = 8 bytes (not a multiple of 3). Translates to MGX (trailing X from
        // incomplete codon). After stripping X → MG. Expected = MG*. MG is shorter
        // than MG* → semi_global_align finds MG contained within MG* at offset 0.
        //
        // ATG GGT GG → M G X (8 bytes, 2 complete codons + 2 leftover)
        let cdna = b"ATGGGTGG";
        let expected = b"MG*".to_vec();

        let mut coding = CodingRegion {
            genomic_start: 1,
            genomic_end: 8,
            cdna_start: 1,
            cdna_end: 8,
            protein_id: "ENSTEST.1".to_string(),
            protein_seq: expected,
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let table = CodonTable::standard();
        let mut stats = EvaluationStats::default();
        evaluate_single(&mut coding, cdna, &table, false, "ENST_TEST.1", &mut stats).unwrap();
        assert_eq!(stats.num_contained, 1);
    }

    #[test]
    fn short_trailing_x_single_aa() {
        // CDS = 5 bytes. Translates to MX. After stripping → M.
        // Expected = M*. M is contained within M*.
        //
        // ATG GG → M X (5 bytes, 1 complete codon + 2 leftover)
        let cdna = b"ATGGG";
        let expected = b"M*".to_vec();

        let mut coding = CodingRegion {
            genomic_start: 1,
            genomic_end: 5,
            cdna_start: 1,
            cdna_end: 5,
            protein_id: "ENSTEST.1".to_string(),
            protein_seq: expected,
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        };

        let table = CodonTable::standard();
        let mut stats = EvaluationStats::default();
        evaluate_single(&mut coding, cdna, &table, false, "ENST_TEST.1", &mut stats).unwrap();
        assert_eq!(stats.num_contained, 1);
    }
}
