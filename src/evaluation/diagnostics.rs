//! Diagnostic output for unresolvable transcripts (spec Section 8.5).

use crate::cli;
use crate::codon::{self, CodonTable};
use crate::transcript::types::CodingRegion;

/// Dump diagnostic information for an unresolvable transcript.
pub(super) fn dump_unresolvable_diagnostics(
    coding: &CodingRegion,
    cdna: &[u8],
    table: &CodonTable,
) {
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
