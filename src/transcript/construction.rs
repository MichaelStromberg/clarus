//! Transcript construction: exon normalization, intron insertion, coding region detection.

use crate::error::Error;
use crate::genbank::CdsRecord;
use crate::gff3::entry::{CigarOpType, Gff3Entry, MatchRecord};
use crate::strand::Strand;

use super::types::{CodingRegion, TranscriptRegion, TranscriptRegionType};

/// Normalize exons and insert introns for a transcript.
///
/// Chooses between forward/reverse simple path and CIGAR-aware path
/// based on presence of cDNA_match records.
pub fn build_transcript_regions(
    exons: &[Gff3Entry],
    matches: &[MatchRecord],
    strand: Strand,
) -> Result<Vec<TranscriptRegion>, Error> {
    let mut regions = if matches.is_empty() {
        build_exon_regions(exons, strand)
    } else {
        build_cigar_aware(matches, strand)?
    };

    // Insert introns between consecutive exons
    insert_introns(&mut regions, strand);
    regions.sort_by_key(|r| r.genomic_start);
    Ok(regions)
}

/// Build exon regions with cDNA coordinates assigned in transcript order.
///
/// Forward strand: exons sorted ascending by genomic start (5'→3' matches genome).
/// Reverse strand: exons sorted descending by genomic start (5'→3' is opposite to genome),
/// then re-sorted ascending for output.
fn build_exon_regions(exons: &[Gff3Entry], strand: Strand) -> Vec<TranscriptRegion> {
    let mut sorted: Vec<&Gff3Entry> = exons.iter().collect();
    if strand.is_reverse() {
        sorted.sort_by(|a, b| b.start.cmp(&a.start));
    } else {
        sorted.sort_by_key(|e| e.start);
    }

    let mut regions = Vec::with_capacity(sorted.len());
    let mut cdna_start: i32 = 1;
    let mut id: u16 = 1;

    for exon in &sorted {
        let genomic_length = exon.end - exon.start;
        let cdna_end = cdna_start + genomic_length;

        regions.push(TranscriptRegion {
            region_type: TranscriptRegionType::Exon,
            id,
            genomic_start: exon.start,
            genomic_end: exon.end,
            cdna_start,
            cdna_end,
            cigar_ops: None,
        });

        cdna_start = cdna_end + 1;
        id += 1;
    }

    if strand.is_reverse() {
        regions.sort_by_key(|r| r.genomic_start);
    }
    regions
}

/// CIGAR-aware exon construction from cDNA_match records.
///
/// Each MatchRecord becomes a single exon region. Genomic coordinates come from
/// the entry's start/end, cDNA coordinates from Target start/end, and CIGAR ops
/// are preserved on the region for use by `cigar_walk` during coordinate mapping.
/// D/I gaps within a match are NOT split into separate regions — they are
/// alignment artifacts, not real splice sites.
fn build_cigar_aware(
    matches: &[MatchRecord],
    strand: Strand,
) -> Result<Vec<TranscriptRegion>, Error> {
    let mut sorted_matches: Vec<&MatchRecord> = matches.iter().collect();
    sorted_matches.sort_by_key(|m| m.entry.start);

    let mut regions = Vec::with_capacity(sorted_matches.len());

    for m in &sorted_matches {
        let target_start = m
            .entry
            .attributes
            .target_start
            .ok_or_else(|| Error::Parse("cDNA_match missing TargetStart".into()))?;
        let target_end = m
            .entry
            .attributes
            .target_end
            .ok_or_else(|| Error::Parse("cDNA_match missing TargetEnd".into()))?;

        regions.push(TranscriptRegion {
            region_type: TranscriptRegionType::Exon,
            id: 0, // assigned below
            genomic_start: m.entry.start,
            genomic_end: m.entry.end,
            cdna_start: target_start,
            cdna_end: target_end,
            cigar_ops: m.entry.attributes.cigar_ops.clone(),
        });
    }

    // Assign exon IDs sorted by genomic position.
    // Exon counts per transcript are well under u16::MAX for all known genomes.
    regions.sort_by_key(|r| r.genomic_start);
    let count = regions.len() as u16;
    if strand.is_reverse() {
        for (i, region) in regions.iter_mut().enumerate() {
            region.id = count - (i as u16);
        }
    } else {
        for (i, region) in regions.iter_mut().enumerate() {
            region.id = (i as u16) + 1;
        }
    }

    Ok(regions)
}

/// Insert introns between consecutive exon regions.
fn insert_introns(regions: &mut Vec<TranscriptRegion>, strand: Strand) {
    // Sort by genomic start to find gaps
    regions.sort_by_key(|r| r.genomic_start);

    let mut introns = Vec::new();

    for window in regions.windows(2) {
        let prev = &window[0];
        let curr = &window[1];
        let gap = curr.genomic_start - prev.genomic_end - 1;
        if gap <= 0 {
            continue;
        }

        let (cdna_start, cdna_end, id) = if strand.is_reverse() {
            (curr.cdna_end, prev.cdna_start, curr.id)
        } else {
            (prev.cdna_end, curr.cdna_start, prev.id)
        };

        introns.push(TranscriptRegion {
            region_type: TranscriptRegionType::Intron,
            id,
            genomic_start: prev.genomic_end + 1,
            genomic_end: curr.genomic_start - 1,
            cdna_start,
            cdna_end,
            cigar_ops: None,
        });
    }

    regions.extend(introns);
}

/// Detect coding region boundaries from CDS entries and map to cDNA coordinates.
pub fn detect_coding_region(
    cds_entries: &[Gff3Entry],
    transcript_regions: &[TranscriptRegion],
    strand: Strand,
    protein_id: &str,
    protein_seq: Vec<u8>,
    genbank_cds: Option<&CdsRecord>,
) -> Result<CodingRegion, Error> {
    if cds_entries.is_empty() {
        return Err(Error::Parse("no CDS entries for coding region".to_string()));
    }

    // Find overall CDS genomic boundaries
    let genomic_start = cds_entries
        .iter()
        .map(|c| c.start)
        .min()
        .expect("invariant: slice non-empty after is_empty guard");
    let genomic_end = cds_entries
        .iter()
        .map(|c| c.end)
        .max()
        .expect("invariant: slice non-empty after is_empty guard");

    // Map to cDNA coordinates (strand-aware)
    let (map_start, map_end) = if strand.is_reverse() {
        (genomic_end, genomic_start)
    } else {
        (genomic_start, genomic_end)
    };

    let cdna_start = map_genomic_to_cdna(map_start, transcript_regions, strand)?;
    let cdna_end = map_genomic_to_cdna(map_end, transcript_regions, strand)?;

    // CDS offset from GenBank comparison
    let cds_offset = match genbank_cds {
        Some(gb) if cdna_start > gb.cdna_start => (cdna_start - gb.cdna_start) as u16,
        _ => 0,
    };

    Ok(CodingRegion {
        genomic_start,
        genomic_end,
        cdna_start,
        cdna_end,
        protein_id: protein_id.to_string(),
        protein_seq,
        cds_padding: 0,
        cds_offset,
        protein_offset: 0,
        amino_acid_edits: None,
        slip: None,
    })
}

/// Walk CIGAR operations to map a genomic position to cDNA within a CIGAR-bearing region.
///
/// Implements spec `05_coordinate_mapping.md` Section 3.2.
fn cigar_walk(
    variant_pos: i32,
    region: &TranscriptRegion,
    ops: &[crate::gff3::entry::CigarOp],
    strand: Strand,
) -> i32 {
    // Step 1: Initialize genomic_pos just before exon start
    let mut genomic_pos = if strand.is_reverse() {
        region.genomic_end + 1
    } else {
        region.genomic_start - 1
    };

    // Step 2: Initialize cdna_pos
    let mut cdna_pos = region.cdna_start - 1;

    // Step 3: Walk each CIGAR op
    for op in ops {
        // CIGAR op lengths are small (genomic alignment segments within exons);
        // safe to convert to i32 for signed coordinate arithmetic on reverse strand.
        debug_assert!(
            i32::try_from(op.length).is_ok(),
            "CIGAR op length exceeds i32::MAX"
        );
        let len = op.length as i32;
        let (delta_genome, delta_cdna) = match op.op_type {
            CigarOpType::Match => {
                let dg = if strand.is_reverse() { -len } else { len };
                (dg, len)
            }
            CigarOpType::Insertion => (0, len),
            CigarOpType::Deletion => {
                let dg = if strand.is_reverse() { -len } else { len };
                (dg, 0)
            }
        };

        // Step 4: Next genomic position after this op
        let next_genomic_pos = genomic_pos + delta_genome;

        // Step 5: Check whether target position is reached within this op
        let reached = if strand.is_reverse() {
            next_genomic_pos <= variant_pos
        } else {
            next_genomic_pos >= variant_pos
        };

        if reached {
            // Step 6: Compute remaining distance and constrained deltas
            let remaining = (variant_pos - genomic_pos).abs();
            let constrained_cdna = match op.op_type {
                CigarOpType::Match | CigarOpType::Insertion => remaining,
                CigarOpType::Deletion => 0,
            };
            cdna_pos += constrained_cdna;
            return cdna_pos;
        }

        // Step 7: Accumulate full deltas
        genomic_pos += delta_genome;
        cdna_pos += delta_cdna;
    }

    cdna_pos
}

/// Map a genomic position to a cDNA position through transcript regions.
///
/// Dispatches to `cigar_walk` when the exon has CIGAR ops, otherwise uses
/// the linear offset formula.
fn map_genomic_to_cdna(
    genomic_pos: i32,
    regions: &[TranscriptRegion],
    strand: Strand,
) -> Result<i32, Error> {
    for region in regions {
        if region.region_type != TranscriptRegionType::Exon {
            continue;
        }
        if genomic_pos >= region.genomic_start && genomic_pos <= region.genomic_end {
            if let Some(ops) = &region.cigar_ops {
                return Ok(cigar_walk(genomic_pos, region, ops, strand));
            }
            let cdna_pos = if strand.is_reverse() {
                region.cdna_start + (region.genomic_end - genomic_pos)
            } else {
                region.cdna_start + (genomic_pos - region.genomic_start)
            };
            return Ok(cdna_pos);
        }
    }
    Err(Error::Parse(format!(
        "genomic position {genomic_pos} not found in any exon region"
    )))
}

#[cfg(test)]
mod tests {
    use crate::biotype::BioType;
    use crate::gff3::entry::{CigarOp, CigarOpType, Gff3Attributes};

    use super::*;

    fn make_exon(start: i32, end: i32) -> Gff3Entry {
        Gff3Entry {
            chromosome_index: 0,
            start,
            end,
            biotype: BioType::Exon,
            strand: Strand::Forward,
            attributes: Gff3Attributes::default(),
        }
    }

    #[test]
    fn forward_strand_normalization() {
        // Spec Section 10 worked example
        let exons = vec![
            make_exon(1000, 1200),
            make_exon(1500, 1700),
            make_exon(2000, 2300),
        ];
        let regions = build_transcript_regions(&exons, &[], Strand::Forward).unwrap();

        // 3 exons + 2 introns = 5 regions
        assert_eq!(regions.len(), 5);

        // Exon 1
        assert_eq!(regions[0].region_type, TranscriptRegionType::Exon);
        assert_eq!(regions[0].genomic_start, 1000);
        assert_eq!(regions[0].genomic_end, 1200);
        assert_eq!(regions[0].cdna_start, 1);
        assert_eq!(regions[0].cdna_end, 201);
        assert_eq!(regions[0].id, 1);

        // Intron 1
        assert_eq!(regions[1].region_type, TranscriptRegionType::Intron);
        assert_eq!(regions[1].genomic_start, 1201);
        assert_eq!(regions[1].genomic_end, 1499);
        assert_eq!(regions[1].cdna_start, 201);
        assert_eq!(regions[1].cdna_end, 202);
        assert_eq!(regions[1].id, 1);

        // Exon 2
        assert_eq!(regions[2].region_type, TranscriptRegionType::Exon);
        assert_eq!(regions[2].genomic_start, 1500);
        assert_eq!(regions[2].cdna_start, 202);
        assert_eq!(regions[2].cdna_end, 402);
        assert_eq!(regions[2].id, 2);

        // Intron 2
        assert_eq!(regions[3].region_type, TranscriptRegionType::Intron);

        // Exon 3
        assert_eq!(regions[4].region_type, TranscriptRegionType::Exon);
        assert_eq!(regions[4].cdna_start, 403);
        assert_eq!(regions[4].cdna_end, 703);
        assert_eq!(regions[4].id, 3);
    }

    #[test]
    fn reverse_strand_normalization() {
        // Spec Section 11 worked example
        let exons = vec![
            make_exon(5000, 5100),
            make_exon(5400, 5600),
            make_exon(6000, 6200),
        ];
        let regions = build_transcript_regions(&exons, &[], Strand::Reverse).unwrap();

        assert_eq!(regions.len(), 5);

        // After sorting by genomic start:
        // Exon X (5000-5100) should have highest cDNA
        assert_eq!(regions[0].genomic_start, 5000);
        assert_eq!(regions[0].cdna_start, 403);
        assert_eq!(regions[0].cdna_end, 503);
        assert_eq!(regions[0].id, 3);

        // Intron between X and Y
        assert_eq!(regions[1].region_type, TranscriptRegionType::Intron);
        assert_eq!(regions[1].cdna_start, 402);
        assert_eq!(regions[1].cdna_end, 403);
        assert_eq!(regions[1].id, 2);

        // Exon Y (5400-5600)
        assert_eq!(regions[2].genomic_start, 5400);
        assert_eq!(regions[2].cdna_start, 202);
        assert_eq!(regions[2].cdna_end, 402);
        assert_eq!(regions[2].id, 2);

        // Intron between Y and Z
        assert_eq!(regions[3].region_type, TranscriptRegionType::Intron);
        assert_eq!(regions[3].cdna_start, 201);
        assert_eq!(regions[3].cdna_end, 202);
        assert_eq!(regions[3].id, 1);

        // Exon Z (6000-6200) should have lowest cDNA
        assert_eq!(regions[4].genomic_start, 6000);
        assert_eq!(regions[4].cdna_start, 1);
        assert_eq!(regions[4].cdna_end, 201);
        assert_eq!(regions[4].id, 1);
    }

    #[test]
    fn coding_region_forward() {
        let exons = vec![
            make_exon(1000, 1200),
            make_exon(1500, 1700),
            make_exon(2000, 2300),
        ];
        let regions = build_transcript_regions(&exons, &[], Strand::Forward).unwrap();

        // CDS spanning 1050..2250
        let cds = vec![make_exon(1050, 2250)]; // Using make_exon for CDS geometry
        let coding = detect_coding_region(
            &cds,
            &regions,
            Strand::Forward,
            "NP_001.1",
            b"MACK*".to_vec(),
            None,
        )
        .unwrap();

        assert_eq!(coding.genomic_start, 1050);
        assert_eq!(coding.genomic_end, 2250);
        assert_eq!(coding.cdna_start, 51); // 1 + (1050 - 1000) = 51
        assert_eq!(coding.cdna_end, 653); // 403 + (2250 - 2000) = 653
    }

    #[test]
    fn coding_region_reverse() {
        let exons = vec![
            make_exon(5000, 5100),
            make_exon(5400, 5600),
            make_exon(6000, 6200),
        ];
        let regions = build_transcript_regions(&exons, &[], Strand::Reverse).unwrap();

        let cds = vec![make_exon(5050, 6150)];
        let coding = detect_coding_region(
            &cds,
            &regions,
            Strand::Reverse,
            "NP_001.1",
            b"MACK*".to_vec(),
            None,
        )
        .unwrap();

        // Reverse strand: map_start = 6150, map_end = 5050
        assert_eq!(coding.cdna_start, 51); // 1 + (6200 - 6150) = 51
        assert_eq!(coding.cdna_end, 453); // 403 + (5100 - 5050) = 453
    }

    #[test]
    fn cds_offset_computation() {
        let exons = vec![make_exon(1000, 1200)];
        let regions = build_transcript_regions(&exons, &[], Strand::Forward).unwrap();

        let cds = vec![make_exon(1074, 1200)];
        let gb_cds = CdsRecord {
            cdna_start: 50,
            cdna_end: 200,
            codon_start: 1,
        };

        let coding = detect_coding_region(
            &cds,
            &regions,
            Strand::Forward,
            "NP_001.1",
            b"MACK*".to_vec(),
            Some(&gb_cds),
        )
        .unwrap();

        // cdna_start = 1 + (1074 - 1000) = 75
        // cds_offset = 75 - 50 = 25
        assert_eq!(coding.cdna_start, 75);
        assert_eq!(coding.cds_offset, 25);
    }

    fn make_match(
        start: i32,
        end: i32,
        target_id: &str,
        target_start: i32,
        target_end: i32,
        cigar: Option<Vec<CigarOp>>,
    ) -> MatchRecord {
        MatchRecord {
            entry: Gff3Entry {
                chromosome_index: 0,
                start,
                end,
                biotype: BioType::CdnaMatch,
                strand: Strand::Forward,
                attributes: Gff3Attributes {
                    target_id: Some(target_id.to_string()),
                    target_start: Some(target_start),
                    target_end: Some(target_end),
                    cigar_ops: cigar,
                    ..Default::default()
                },
            },
            exons: Vec::new(),
        }
    }

    #[test]
    fn cigar_aware_single_region_with_ops() {
        // Models NM_001396027.1: Gap=M241 D1 M35 D1 M420
        // Genomic span = 698 bp, cDNA span = 696 bp
        // Single exon region with CIGAR ops preserved for cigar_walk mapping.
        let cigar = vec![
            CigarOp {
                op_type: CigarOpType::Match,
                length: 241,
            },
            CigarOp {
                op_type: CigarOpType::Deletion,
                length: 1,
            },
            CigarOp {
                op_type: CigarOpType::Match,
                length: 35,
            },
            CigarOp {
                op_type: CigarOpType::Deletion,
                length: 1,
            },
            CigarOp {
                op_type: CigarOpType::Match,
                length: 420,
            },
        ];
        let cigar_clone = cigar.clone();
        let matches = vec![make_match(1000, 1697, "NM_TEST.1", 1, 696, Some(cigar))];
        let regions = build_transcript_regions(&[], &matches, Strand::Forward).unwrap();

        // Single exon region (no fake introns from D gaps)
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].region_type, TranscriptRegionType::Exon);
        assert_eq!(regions[0].genomic_start, 1000);
        assert_eq!(regions[0].genomic_end, 1697);
        assert_eq!(regions[0].cdna_start, 1);
        assert_eq!(regions[0].cdna_end, 696);
        assert_eq!(regions[0].id, 1);
        assert_eq!(regions[0].cigar_ops, Some(cigar_clone));
    }

    #[test]
    fn cigar_aware_no_gap_single_region() {
        // No CIGAR ops — single region (existing behavior)
        let matches = vec![make_match(1000, 1200, "NM_TEST.1", 1, 201, None)];
        let regions = build_transcript_regions(&[], &matches, Strand::Forward).unwrap();

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].genomic_start, 1000);
        assert_eq!(regions[0].genomic_end, 1200);
        assert_eq!(regions[0].cdna_start, 1);
        assert_eq!(regions[0].cdna_end, 201);
    }

    #[test]
    fn cigar_aware_single_region_with_insertion() {
        // Gap=M100 I2 M100 — insertion of 2 bases in cDNA (not in genome)
        // Single exon region with CIGAR ops; cigar_walk handles the I gap.
        let cigar = vec![
            CigarOp {
                op_type: CigarOpType::Match,
                length: 100,
            },
            CigarOp {
                op_type: CigarOpType::Insertion,
                length: 2,
            },
            CigarOp {
                op_type: CigarOpType::Match,
                length: 100,
            },
        ];
        let cigar_clone = cigar.clone();
        let matches = vec![make_match(1000, 1199, "NM_TEST.1", 1, 202, Some(cigar))];
        let regions = build_transcript_regions(&[], &matches, Strand::Forward).unwrap();

        // Single exon region (I ops don't split into sub-regions)
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].region_type, TranscriptRegionType::Exon);
        assert_eq!(regions[0].genomic_start, 1000);
        assert_eq!(regions[0].genomic_end, 1199);
        assert_eq!(regions[0].cdna_start, 1);
        assert_eq!(regions[0].cdna_end, 202);
        assert_eq!(regions[0].id, 1);
        assert_eq!(regions[0].cigar_ops, Some(cigar_clone));
    }

    #[test]
    fn cds_mapping_with_cigar_deletions() {
        // Build regions from CIGAR-aware match with D ops (NM_001396027.1 model).
        // CDS spans the full region. cigar_walk produces cdna_end=696 (not 698)
        // because D ops are walked, not computed by linear formula.
        let cigar = vec![
            CigarOp {
                op_type: CigarOpType::Match,
                length: 241,
            },
            CigarOp {
                op_type: CigarOpType::Deletion,
                length: 1,
            },
            CigarOp {
                op_type: CigarOpType::Match,
                length: 35,
            },
            CigarOp {
                op_type: CigarOpType::Deletion,
                length: 1,
            },
            CigarOp {
                op_type: CigarOpType::Match,
                length: 420,
            },
        ];
        let matches = vec![make_match(1000, 1697, "NM_TEST.1", 1, 696, Some(cigar))];
        let regions = build_transcript_regions(&[], &matches, Strand::Forward).unwrap();

        // CDS spanning the full region
        let cds = vec![make_exon(1000, 1697)];
        let coding = detect_coding_region(
            &cds,
            &regions,
            Strand::Forward,
            "NP_TEST.1",
            b"MACK*".to_vec(),
            None,
        )
        .unwrap();

        assert_eq!(coding.cdna_start, 1);
        assert_eq!(coding.cdna_end, 696); // NOT 698
    }

    #[test]
    fn cigar_aware_reverse_strand_cds_mapping() {
        // Reverse-strand cDNA_match with D ops.
        // CIGAR: M241 D1 M35 D1 M420 (same as forward test)
        // Single region with CIGAR ops; cigar_walk handles reverse strand mapping.
        let cigar = vec![
            CigarOp {
                op_type: CigarOpType::Match,
                length: 241,
            },
            CigarOp {
                op_type: CigarOpType::Deletion,
                length: 1,
            },
            CigarOp {
                op_type: CigarOpType::Match,
                length: 35,
            },
            CigarOp {
                op_type: CigarOpType::Deletion,
                length: 1,
            },
            CigarOp {
                op_type: CigarOpType::Match,
                length: 420,
            },
        ];
        let matches = vec![make_match(1000, 1697, "NM_TEST.1", 1, 696, Some(cigar))];
        let regions = build_transcript_regions(&[], &matches, Strand::Reverse).unwrap();

        // Single exon region
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].region_type, TranscriptRegionType::Exon);
        assert_eq!(regions[0].genomic_start, 1000);
        assert_eq!(regions[0].genomic_end, 1697);
        assert_eq!(regions[0].cdna_start, 1);
        assert_eq!(regions[0].cdna_end, 696);
        assert_eq!(regions[0].id, 1);

        // CDS spanning full range maps correctly via cigar_walk
        let cds = vec![make_exon(1000, 1697)];
        let coding = detect_coding_region(
            &cds,
            &regions,
            Strand::Reverse,
            "NP_TEST.1",
            b"M".to_vec(),
            None,
        )
        .unwrap();
        assert_eq!(coding.cdna_start, 1);
        assert_eq!(coding.cdna_end, 696);
    }

    #[test]
    fn cigar_walk_maps_through_deletions() {
        // Test map_genomic_to_cdna at specific positions through M241/D1/M35/D1/M420.
        let cigar = vec![
            CigarOp {
                op_type: CigarOpType::Match,
                length: 241,
            },
            CigarOp {
                op_type: CigarOpType::Deletion,
                length: 1,
            },
            CigarOp {
                op_type: CigarOpType::Match,
                length: 35,
            },
            CigarOp {
                op_type: CigarOpType::Deletion,
                length: 1,
            },
            CigarOp {
                op_type: CigarOpType::Match,
                length: 420,
            },
        ];
        let matches = vec![make_match(1000, 1697, "NM_TEST.1", 1, 696, Some(cigar))];
        let regions = build_transcript_regions(&[], &matches, Strand::Forward).unwrap();

        // First base of M241
        assert_eq!(
            map_genomic_to_cdna(1000, &regions, Strand::Forward).unwrap(),
            1
        );
        // Last base of M241
        assert_eq!(
            map_genomic_to_cdna(1240, &regions, Strand::Forward).unwrap(),
            241
        );
        // D1 gap — maps to last cDNA before deletion
        assert_eq!(
            map_genomic_to_cdna(1241, &regions, Strand::Forward).unwrap(),
            241
        );
        // First base of M35
        assert_eq!(
            map_genomic_to_cdna(1242, &regions, Strand::Forward).unwrap(),
            242
        );
        // Last base of M35
        assert_eq!(
            map_genomic_to_cdna(1276, &regions, Strand::Forward).unwrap(),
            276
        );
        // D1 gap — maps to last cDNA before deletion
        assert_eq!(
            map_genomic_to_cdna(1277, &regions, Strand::Forward).unwrap(),
            276
        );
        // First base of M420
        assert_eq!(
            map_genomic_to_cdna(1278, &regions, Strand::Forward).unwrap(),
            277
        );
        // Last base of M420
        assert_eq!(
            map_genomic_to_cdna(1697, &regions, Strand::Forward).unwrap(),
            696
        );
    }

    #[test]
    fn cigar_walk_maps_through_insertion() {
        // Test map_genomic_to_cdna through M100/I2/M100.
        // cDNA positions 101-102 are insertion bases with no genomic counterpart.
        let cigar = vec![
            CigarOp {
                op_type: CigarOpType::Match,
                length: 100,
            },
            CigarOp {
                op_type: CigarOpType::Insertion,
                length: 2,
            },
            CigarOp {
                op_type: CigarOpType::Match,
                length: 100,
            },
        ];
        let matches = vec![make_match(1000, 1199, "NM_TEST.1", 1, 202, Some(cigar))];
        let regions = build_transcript_regions(&[], &matches, Strand::Forward).unwrap();

        // First base
        assert_eq!(
            map_genomic_to_cdna(1000, &regions, Strand::Forward).unwrap(),
            1
        );
        // Last base of first M block
        assert_eq!(
            map_genomic_to_cdna(1099, &regions, Strand::Forward).unwrap(),
            100
        );
        // First base after insertion — cDNA skips I2 bases (101, 102)
        assert_eq!(
            map_genomic_to_cdna(1100, &regions, Strand::Forward).unwrap(),
            103
        );
        // Last base
        assert_eq!(
            map_genomic_to_cdna(1199, &regions, Strand::Forward).unwrap(),
            202
        );
    }
}
