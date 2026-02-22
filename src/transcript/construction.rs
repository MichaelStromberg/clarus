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
    let mut regions = if !matches.is_empty() {
        build_cigar_aware(matches, strand)?
    } else if strand.is_reverse() {
        build_reverse(exons)
    } else {
        build_forward(exons)
    };

    // Insert introns between consecutive exons
    insert_introns(&mut regions, strand);
    regions.sort_by_key(|r| r.genomic_start);
    Ok(regions)
}

/// Forward strand exon normalization.
fn build_forward(exons: &[Gff3Entry]) -> Vec<TranscriptRegion> {
    let mut sorted: Vec<&Gff3Entry> = exons.iter().collect();
    sorted.sort_by_key(|e| e.start);

    let mut regions = Vec::with_capacity(sorted.len());
    let mut cdna_start: i32 = 1;
    let mut id: u16 = 1;

    for exon in sorted {
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

    regions
}

/// Reverse strand exon normalization.
fn build_reverse(exons: &[Gff3Entry]) -> Vec<TranscriptRegion> {
    let mut sorted: Vec<&Gff3Entry> = exons.iter().collect();
    // Sort descending by start for cDNA assignment (5'→3' in transcript order)
    sorted.sort_by(|a, b| b.start.cmp(&a.start));

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

    // Re-sort by ascending genomic start for output
    regions.sort_by_key(|r| r.genomic_start);
    regions
}

/// Count the number of M (match) blocks in a MatchRecord's CIGAR.
fn count_m_blocks(m: &MatchRecord) -> u16 {
    match &m.entry.attributes.cigar_ops {
        Some(ops) => ops
            .iter()
            .filter(|o| o.op_type == CigarOpType::Match)
            .count() as u16,
        None => 1,
    }
}

/// CIGAR-aware exon normalization from cDNA_match records.
///
/// When CIGAR ops contain D (deletion) or I (insertion) ops, the match is split
/// into separate sub-regions at each D/I boundary. Each M (match) block becomes
/// its own exon region with a 1:1 genomic-to-cDNA mapping, allowing the linear
/// formula in `map_genomic_to_cdna` to work correctly.
fn build_cigar_aware(
    matches: &[MatchRecord],
    strand: Strand,
) -> Result<Vec<TranscriptRegion>, Error> {
    let mut sorted_matches: Vec<&MatchRecord> = matches.iter().collect();
    sorted_matches.sort_by_key(|m| m.entry.start);

    let mut regions = Vec::new();
    let total_blocks: u16 = if strand.is_reverse() {
        sorted_matches.iter().map(|m| count_m_blocks(m)).sum()
    } else {
        0
    };
    let mut id: u16 = if strand.is_reverse() { total_blocks } else { 1 };

    for m in &sorted_matches {
        let target_start = m
            .entry
            .attributes
            .target_start
            .ok_or_else(|| Error::Parse("cDNA_match missing TargetStart".to_string()))?;
        let target_end = m
            .entry
            .attributes
            .target_end
            .ok_or_else(|| Error::Parse("cDNA_match missing TargetEnd".to_string()))?;

        if let Some(ref cigar_ops) = m.entry.attributes.cigar_ops {
            let mut genomic_pos = m.entry.start;
            let mut cdna_pos = if strand.is_reverse() {
                target_end
            } else {
                target_start
            };

            for op in cigar_ops {
                let len = op.length as i32;
                match op.op_type {
                    CigarOpType::Match => {
                        let (cs, ce) = if strand.is_reverse() {
                            (cdna_pos - len + 1, cdna_pos)
                        } else {
                            (cdna_pos, cdna_pos + len - 1)
                        };
                        regions.push(TranscriptRegion {
                            region_type: TranscriptRegionType::Exon,
                            id,
                            genomic_start: genomic_pos,
                            genomic_end: genomic_pos + len - 1,
                            cdna_start: cs,
                            cdna_end: ce,
                            cigar_ops: None,
                        });
                        genomic_pos += len;
                        if strand.is_reverse() {
                            cdna_pos -= len;
                            id = id.saturating_sub(1);
                        } else {
                            cdna_pos += len;
                            id += 1;
                        }
                    }
                    CigarOpType::Deletion => {
                        genomic_pos += len;
                    }
                    CigarOpType::Insertion => {
                        if strand.is_reverse() {
                            cdna_pos -= len;
                        } else {
                            cdna_pos += len;
                        }
                    }
                }
            }
        } else {
            // No CIGAR — single region (preserves existing behavior)
            regions.push(TranscriptRegion {
                region_type: TranscriptRegionType::Exon,
                id,
                genomic_start: m.entry.start,
                genomic_end: m.entry.end,
                cdna_start: target_start,
                cdna_end: target_end,
                cigar_ops: None,
            });
            if strand.is_reverse() {
                id = id.saturating_sub(1);
            } else {
                id += 1;
            }
        }
    }

    Ok(regions)
}

/// Insert introns between consecutive exon regions.
fn insert_introns(regions: &mut Vec<TranscriptRegion>, strand: Strand) {
    // Sort by genomic start to find gaps
    regions.sort_by_key(|r| r.genomic_start);

    let exons: Vec<TranscriptRegion> = regions.clone();
    let mut introns = Vec::new();

    for window in exons.windows(2) {
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
    let genomic_start = cds_entries.iter().map(|c| c.start).min().unwrap();
    let genomic_end = cds_entries.iter().map(|c| c.end).max().unwrap();

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

/// Map a genomic position to a cDNA position through transcript regions.
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
    use crate::gff3::entry::{CigarOp, Gff3Attributes};

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
    fn cigar_aware_splits_at_deletions() {
        // Models NM_001396027.1: Gap=M241 D1 M35 D1 M420
        // Genomic span = 698 bp, cDNA span = 696 bp
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

        // 3 exon sub-regions + 2 introns = 5 regions
        assert_eq!(regions.len(), 5);

        // Exon 1: M241
        assert_eq!(regions[0].region_type, TranscriptRegionType::Exon);
        assert_eq!(regions[0].genomic_start, 1000);
        assert_eq!(regions[0].genomic_end, 1240);
        assert_eq!(regions[0].cdna_start, 1);
        assert_eq!(regions[0].cdna_end, 241);
        assert_eq!(regions[0].id, 1);

        // Intron 1 (D1 gap)
        assert_eq!(regions[1].region_type, TranscriptRegionType::Intron);
        assert_eq!(regions[1].genomic_start, 1241);
        assert_eq!(regions[1].genomic_end, 1241);

        // Exon 2: M35
        assert_eq!(regions[2].region_type, TranscriptRegionType::Exon);
        assert_eq!(regions[2].genomic_start, 1242);
        assert_eq!(regions[2].genomic_end, 1276);
        assert_eq!(regions[2].cdna_start, 242);
        assert_eq!(regions[2].cdna_end, 276);
        assert_eq!(regions[2].id, 2);

        // Intron 2 (D1 gap)
        assert_eq!(regions[3].region_type, TranscriptRegionType::Intron);
        assert_eq!(regions[3].genomic_start, 1277);
        assert_eq!(regions[3].genomic_end, 1277);

        // Exon 3: M420
        assert_eq!(regions[4].region_type, TranscriptRegionType::Exon);
        assert_eq!(regions[4].genomic_start, 1278);
        assert_eq!(regions[4].genomic_end, 1697);
        assert_eq!(regions[4].cdna_start, 277);
        assert_eq!(regions[4].cdna_end, 696);
        assert_eq!(regions[4].id, 3);

        // CDS spanning the full range maps correctly
        let cds = vec![make_exon(1000, 1697)];
        let coding = detect_coding_region(
            &cds,
            &regions,
            Strand::Forward,
            "NP_TEST.1",
            b"M".to_vec(),
            None,
        )
        .unwrap();
        assert_eq!(coding.cdna_start, 1);
        assert_eq!(coding.cdna_end, 696); // NOT 698
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
    fn cigar_aware_insertion_op() {
        // Gap=M100 I2 M100 — insertion of 2 bases in cDNA (not in genome)
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

        // 2 sub-regions, NO intron (genomically adjacent, gap=0)
        assert_eq!(regions.len(), 2);

        // Sub-region 1: M100
        assert_eq!(regions[0].genomic_start, 1000);
        assert_eq!(regions[0].genomic_end, 1099);
        assert_eq!(regions[0].cdna_start, 1);
        assert_eq!(regions[0].cdna_end, 100);

        // Sub-region 2: M100 (cDNA jumps by 2 for the I op)
        assert_eq!(regions[1].genomic_start, 1100);
        assert_eq!(regions[1].genomic_end, 1199);
        assert_eq!(regions[1].cdna_start, 103);
        assert_eq!(regions[1].cdna_end, 202);
    }
}
