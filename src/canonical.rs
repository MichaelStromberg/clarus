//! Transcript designation propagation for PAR region duplicates.
//!
//! Propagates source-specific designation flags (MANE Select, MANE Plus Clinical,
//! RefSeq Select, Ensembl Canonical) across PAR region transcript duplicates.

use std::collections::HashMap;

use crate::transcript::types::{Designations, IntermediateTranscript};

/// Reapply designation flags for PAR region transcripts.
///
/// In GFF3, a transcript appearing on both chrX and chrY may be marked with
/// designation flags on only one chromosome. This step collects the union of
/// all designation flags per transcript ID, then applies to all instances.
///
/// Returns the number of individual flag applications made.
pub fn reapply_par_designations(transcripts: &mut [IntermediateTranscript]) -> usize {
    // Collect the union of flags per transcript ID.
    let mut flags_by_id: HashMap<String, Designations> = HashMap::new();

    for tx in transcripts.iter() {
        if let Some(flags) = flags_by_id.get_mut(tx.id.as_str()) {
            *flags |= tx.designations;
        } else {
            flags_by_id.insert(tx.id.clone(), tx.designations);
        }
    }

    // Apply the union back to all instances, counting newly set flags
    let mut count = 0;
    for tx in transcripts.iter_mut() {
        if let Some(&flags) = flags_by_id.get(&tx.id) {
            let before = tx.designations;
            tx.designations |= flags;

            count += usize::from(tx.designations.is_mane_select && !before.is_mane_select);
            count +=
                usize::from(tx.designations.is_mane_plus_clinical && !before.is_mane_plus_clinical);
            count += usize::from(tx.designations.is_refseq_select && !before.is_refseq_select);
            count +=
                usize::from(tx.designations.is_ensembl_canonical && !before.is_ensembl_canonical);
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::biotype::BioType;
    use crate::strand::Strand;
    use crate::transcript::types::{Gene, Source};
    use std::sync::Arc;

    fn make_tx(id: &str, chromosome_index: usize) -> IntermediateTranscript {
        IntermediateTranscript {
            chromosome_index,
            start: 1,
            end: 100,
            id: id.to_string(),
            biotype: BioType::MRna,
            source: Source::RefSeq,
            strand: Strand::Forward,
            gene: Arc::new(Gene {
                chromosome_index,
                symbol: "SHOX".to_string(),
                ncbi_gene_id: Some("6473".to_string()),
                ensembl_id: None,
                hgnc_id: None,
                on_reverse_strand: false,
            }),
            transcript_regions: vec![],
            coding_region: None,
            cdna_seq: vec![b'A'; 100],
            designations: Designations::default(),
        }
    }

    #[test]
    fn par_reapplies_all_four_flags() {
        let mut tx_x = make_tx("NM_018390.4", 22);
        tx_x.designations.is_mane_select = true;
        tx_x.designations.is_mane_plus_clinical = true;
        tx_x.designations.is_refseq_select = true;
        tx_x.designations.is_ensembl_canonical = true;

        let tx_y = make_tx("NM_018390.4", 23);

        let mut txs = vec![tx_x, tx_y];
        let count = reapply_par_designations(&mut txs);
        assert_eq!(count, 4);
        assert!(txs[1].designations.is_mane_select);
        assert!(txs[1].designations.is_mane_plus_clinical);
        assert!(txs[1].designations.is_refseq_select);
        assert!(txs[1].designations.is_ensembl_canonical);
    }

    #[test]
    fn par_no_reapplication_when_already_set() {
        let mut tx_x = make_tx("NM_018390.4", 22);
        tx_x.designations.is_mane_select = true;
        tx_x.designations.is_refseq_select = true;

        let mut tx_y = make_tx("NM_018390.4", 23);
        tx_y.designations.is_mane_select = true;
        tx_y.designations.is_refseq_select = true;

        let mut txs = vec![tx_x, tx_y];
        let count = reapply_par_designations(&mut txs);
        assert_eq!(count, 0);
    }

    #[test]
    fn par_merges_flags_from_both_copies() {
        // chrX has mane_select, chrY has ensembl_canonical — both should end up with both
        let mut tx_x = make_tx("ENST00000123.1", 22);
        tx_x.designations.is_mane_select = true;

        let mut tx_y = make_tx("ENST00000123.1", 23);
        tx_y.designations.is_ensembl_canonical = true;

        let mut txs = vec![tx_x, tx_y];
        let count = reapply_par_designations(&mut txs);
        assert_eq!(count, 2); // mane_select → chrY, ensembl_canonical → chrX
        assert!(txs[0].designations.is_mane_select);
        assert!(txs[0].designations.is_ensembl_canonical);
        assert!(txs[1].designations.is_mane_select);
        assert!(txs[1].designations.is_ensembl_canonical);
    }
}
