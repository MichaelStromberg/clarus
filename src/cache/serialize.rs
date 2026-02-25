//! Per-type binary serialization for cache records.
//!
//! Each type is serialized using fixed-width encoding: u8, u16, u32
//! for integers and u8-prefixed ASCII for strings.

use std::collections::{HashMap, HashSet};
use std::io::Write;

use crate::error::Error;
use crate::gff3::entry::{CigarOpType, Gff3Entry};
use crate::reference::binary_io::BinaryWrite;
use crate::transcript::types::{
    CodingRegion, Gene, IntermediateTranscript, Source, TranscriptRegion, TranscriptRegionType,
};

// ── Coordinate Helpers ──────────────────────────────────────────

/// Convert a signed coordinate to u32, validating it is non-negative.
fn coord_to_u32(value: i32, field: &str) -> Result<u32, Error> {
    u32::try_from(value).map_err(|_| Error::Validation(format!("{field} is negative: {value}")))
}

// ── Gene ID Helpers ─────────────────────────────────────────────

/// Parse an NCBI gene ID string (e.g. "7157") to a u32.
fn parse_ncbi_gene_id(s: &str) -> Result<u32, Error> {
    s.parse::<u32>()
        .map_err(|_| Error::Parse(format!("invalid NCBI gene ID: '{s}'")))
}

/// Extract the numeric suffix from an Ensembl gene ID (e.g. "ENSG00000141510" → 141510).
fn parse_ensembl_gene_suffix(s: &str) -> Result<u32, Error> {
    let digits = s.trim_start_matches(|c: char| !c.is_ascii_digit());
    let value = digits
        .parse::<u64>()
        .map_err(|_| Error::Parse(format!("invalid Ensembl gene ID: '{s}'")))?;
    u32::try_from(value)
        .map_err(|_| Error::Parse(format!("Ensembl gene ID suffix exceeds u32::MAX: '{s}'")))
}

// ── Gene Serialization ──────────────────────────────────────────

/// Serialize a gene record.
///
/// Format: Symbol (prefixed string) + Flags (u8) + conditional fields.
///
/// Flags: bit 0 = reverse strand, bit 1 = has HGNC ID,
///        bit 2 = has NCBI ID, bit 3 = has Ensembl ID.
pub fn write_gene<W: Write>(writer: &mut W, gene: &Gene) -> Result<(), Error> {
    writer.write_prefixed_string(&gene.symbol)?;

    let mut flags: u8 = 0;
    if gene.on_reverse_strand {
        flags |= 0x01;
    }
    if gene.hgnc_id.is_some() {
        flags |= 0x02;
    }
    if gene.ncbi_gene_id.is_some() {
        flags |= 0x04;
    }
    if gene.ensembl_id.is_some() {
        flags |= 0x08;
    }
    writer.write_u8(flags)?;

    if let Some(ref id) = gene.ncbi_gene_id {
        writer.write_u32(parse_ncbi_gene_id(id)?)?;
    }
    if let Some(ref id) = gene.ensembl_id {
        writer.write_u32(parse_ensembl_gene_suffix(id)?)?;
    }
    if let Some(hgnc_id) = gene.hgnc_id {
        writer.write_u32(hgnc_id)?;
    }

    Ok(())
}

// ── Transcript Region Serialization ──────────────────────────────

/// Serialize a transcript region.
///
/// Format: Start(u32), End(u32), CdnaStart(u32), CdnaEnd(u32), Id(u16),
///         Flags(u8), conditional CIGAR ops.
///
/// Flags: bit 0 = type (0=Exon, 1=Intron), bit 1 = has CIGAR ops.
pub fn write_transcript_region<W: Write>(
    writer: &mut W,
    region: &TranscriptRegion,
) -> Result<(), Error> {
    writer.write_u32(coord_to_u32(region.genomic_start, "region.genomic_start")?)?;
    writer.write_u32(coord_to_u32(region.genomic_end, "region.genomic_end")?)?;
    writer.write_u32(coord_to_u32(region.cdna_start, "region.cdna_start")?)?;
    writer.write_u32(coord_to_u32(region.cdna_end, "region.cdna_end")?)?;
    writer.write_u16(region.id)?;

    let type_bits: u8 = match region.region_type {
        TranscriptRegionType::Exon => 0,
        TranscriptRegionType::Intron => 1,
    };
    let cigar_ops = region.cigar_ops.as_deref().unwrap_or_default();
    let flags = type_bits | if cigar_ops.is_empty() { 0 } else { 0x02 };
    writer.write_u8(flags)?;

    if !cigar_ops.is_empty() {
        writer.write_u8(
            u8::try_from(cigar_ops.len())
                .map_err(|_| Error::Validation("CIGAR ops count exceeds u8::MAX".into()))?,
        )?;
        for op in cigar_ops {
            let op_byte: u8 = match op.op_type {
                CigarOpType::Match => 0,
                CigarOpType::Insertion => 1,
                CigarOpType::Deletion => 2,
            };
            writer.write_u8(op_byte)?;
            writer.write_u16(
                u16::try_from(op.length)
                    .map_err(|_| Error::Validation("CIGAR op length exceeds u16::MAX".into()))?,
            )?;
        }
    }

    Ok(())
}

// ── Coding Region Serialization ──────────────────────────────────

/// Serialize a coding region.
pub fn write_coding_region<W: Write>(
    writer: &mut W,
    cr: &CodingRegion,
    protein_seq_index: u16,
) -> Result<(), Error> {
    writer.write_u32(coord_to_u32(cr.genomic_start, "coding.genomic_start")?)?;
    writer.write_u32(coord_to_u32(cr.genomic_end, "coding.genomic_end")?)?;
    writer.write_u32(coord_to_u32(cr.cdna_start, "coding.cdna_start")?)?;
    writer.write_u32(coord_to_u32(cr.cdna_end, "coding.cdna_end")?)?;
    writer.write_prefixed_string(&cr.protein_id)?;
    writer.write_u16(protein_seq_index)?;
    writer.write_u8(cr.cds_padding)?;
    writer.write_u16(cr.cds_offset)?;
    writer.write_u16(cr.protein_offset)?;

    let mut flags: u8 = 0;
    if cr.amino_acid_edits.is_some() {
        flags |= 0x01;
    }
    if cr.slip.is_some() {
        flags |= 0x02;
    }
    writer.write_u8(flags)?;

    if let Some(ref edits) = cr.amino_acid_edits {
        writer
            .write_u8(u8::try_from(edits.len()).map_err(|_| {
                Error::Validation("amino acid edits count exceeds u8::MAX".into())
            })?)?;
        for edit in edits {
            writer.write_u16(u16::try_from(edit.position).map_err(|_| {
                Error::Validation(format!(
                    "AA edit position exceeds u16::MAX: {}",
                    edit.position
                ))
            })?)?;
            writer.write_u8(edit.amino_acid)?;
        }
    }

    if let Some(ref slip) = cr.slip {
        writer.write_u16(u16::try_from(slip.position).map_err(|_| {
            Error::Validation(format!("slip position exceeds u16::MAX: {}", slip.position))
        })?)?;
        writer.write_u8(slip.length)?;
    }

    Ok(())
}

// ── Transcript Serialization ─────────────────────────────────────

/// Build the 8-bit transcript flags.
///
/// Bit 0: Source (0=Ensembl, 1=RefSeq), bit 1: HasCodingRegion,
/// bit 2: MANE Plus Clinical, bit 3: MANE Select,
/// bit 4: RefSeq Select, bit 5: Ensembl Canonical.
fn build_transcript_flags(tx: &IntermediateTranscript) -> u8 {
    let source_bit: u8 = match tx.source {
        Source::Ensembl => 0,
        Source::RefSeq => 1,
    };
    let has_coding = u8::from(tx.coding_region.is_some());
    let mane_plus_clinical = u8::from(tx.designations.is_mane_plus_clinical);
    let mane_select = u8::from(tx.designations.is_mane_select);
    let refseq_select = u8::from(tx.designations.is_refseq_select);
    let ensembl_canonical = u8::from(tx.designations.is_ensembl_canonical);

    source_bit
        | (has_coding << 1)
        | (mane_plus_clinical << 2)
        | (mane_select << 3)
        | (refseq_select << 4)
        | (ensembl_canonical << 5)
}

/// Serialize a complete transcript record.
pub fn write_transcript<W: Write>(
    writer: &mut W,
    tx: &IntermediateTranscript,
    gene_index: u16,
    region_indices: &[u32],
    cdna_seq_index: u16,
    protein_seq_index: u16,
) -> Result<(), Error> {
    writer.write_u32(coord_to_u32(tx.start, "tx.start")?)?;
    writer.write_u32(coord_to_u32(tx.end, "tx.end")?)?;
    writer.write_prefixed_string(&tx.id)?;

    // Biotype as separate u8
    writer.write_u8(tx.biotype as u8)?;

    writer.write_u16(gene_index)?;

    let flags = build_transcript_flags(tx);
    writer.write_u8(flags)?;

    // Transcript region indices
    writer.write_u16(
        u16::try_from(region_indices.len())
            .map_err(|_| Error::Validation("transcript region count exceeds u16::MAX".into()))?,
    )?;
    for &idx in region_indices {
        writer.write_u32(idx)?;
    }

    // Coding region (conditional)
    if let Some(ref cr) = tx.coding_region {
        write_coding_region(writer, cr, protein_seq_index)?;
    }

    // cDNA sequence index
    writer.write_u16(cdna_seq_index)?;

    Ok(())
}

// ── Regulatory Region Serialization ──────────────────────────────

/// Extract the stable ID from a regulatory region's GFF3 ID attribute.
///
/// Strips any prefix before a colon (e.g. `regulatory_region:ENSR1_958` → `ENSR1_958`).
fn extract_regulatory_id(entry: &Gff3Entry) -> &str {
    match entry.attributes.id.find(':') {
        Some(pos) => &entry.attributes.id[pos + 1..],
        None => &entry.attributes.id,
    }
}

/// Serialize a regulatory region.
///
/// Format: Start(u32), End(u32), BioType(u8), Id(prefixed string),
///         Flags(u8), conditional Note/EcoId/PubMedIds.
///
/// Flags: bit 0 = has ECO ID, bit 1 = has PubMed IDs, bit 2 = has note.
pub fn write_regulatory_region<W: Write>(writer: &mut W, entry: &Gff3Entry) -> Result<(), Error> {
    writer.write_u32(coord_to_u32(entry.start, "regulatory.start")?)?;
    writer.write_u32(coord_to_u32(entry.end, "regulatory.end")?)?;
    writer.write_u8(entry.biotype as u8)?;

    let id = extract_regulatory_id(entry);
    writer.write_prefixed_string(id)?;

    let mut flags: u8 = 0;
    if entry.attributes.eco_id.is_some() {
        flags |= 0x01;
    }
    if entry.attributes.pubmed_ids.is_some() {
        flags |= 0x02;
    }
    if entry.attributes.note.is_some() {
        flags |= 0x04;
    }
    writer.write_u8(flags)?;

    if let Some(ref note) = entry.attributes.note {
        writer.write_u16_prefixed_string(note)?;
    }

    if let Some(eco_id) = entry.attributes.eco_id {
        writer.write_u32(eco_id)?;
    }

    if let Some(ref pubmed_ids) = entry.attributes.pubmed_ids {
        writer.write_u8(
            u8::try_from(pubmed_ids.len())
                .map_err(|_| Error::Validation("PubMed IDs count exceeds u8::MAX".into()))?,
        )?;
        for &id in pubmed_ids {
            writer.write_u32(id)?;
        }
    }

    Ok(())
}

// ── Deduplication Helpers ────────────────────────────────────────

/// Build a deduplicated, sorted array of genes from a set of transcripts,
/// and return the index mapping.
pub fn dedup_genes<'a>(
    transcripts: &[&'a IntermediateTranscript],
) -> (Vec<&'a Gene>, HashMap<Gene, u16>) {
    let mut seen: HashSet<&Gene> = HashSet::new();
    let mut genes: Vec<&Gene> = Vec::new();

    for tx in transcripts {
        let gene = tx.gene.as_ref();
        if seen.insert(gene) {
            genes.push(gene);
        }
    }

    genes.sort_by(|a, b| a.symbol.cmp(&b.symbol));

    let index_map = genes
        .iter()
        .enumerate()
        .map(|(i, gene)| ((*gene).clone(), i as u16)) // safe: count checked by write_u16(try_from) at write site
        .collect();

    (genes, index_map)
}

/// Build a deduplicated, sorted array of transcript regions,
/// and return the index mapping.
pub fn dedup_regions<'a>(
    transcripts: &[&'a IntermediateTranscript],
) -> (Vec<&'a TranscriptRegion>, HashMap<TranscriptRegion, u32>) {
    let mut seen: HashSet<&TranscriptRegion> = HashSet::new();
    let mut regions: Vec<&TranscriptRegion> = Vec::new();

    for tx in transcripts {
        for region in &tx.transcript_regions {
            if seen.insert(region) {
                regions.push(region);
            }
        }
    }

    regions.sort_by(|a, b| {
        a.genomic_start
            .cmp(&b.genomic_start)
            .then(a.genomic_end.cmp(&b.genomic_end))
    });

    let index_map = regions
        .iter()
        .enumerate()
        .map(|(i, region)| ((*region).clone(), i as u32)) // safe: regions per chromosome well under u32::MAX
        .collect();

    (regions, index_map)
}

/// Build a deduplicated, sorted array of sequences (cDNA or protein),
/// and return the index mapping.
pub fn dedup_sequences(seqs: Vec<&[u8]>) -> (Vec<&[u8]>, HashMap<Vec<u8>, u16>) {
    let mut seen: HashSet<&[u8]> = HashSet::new();
    let mut unique: Vec<&[u8]> = Vec::new();

    for seq in seqs {
        if seen.insert(seq) {
            unique.push(seq);
        }
    }

    unique.sort();

    let index_map = unique
        .iter()
        .enumerate()
        .map(|(i, seq)| (seq.to_vec(), i as u16)) // safe: count checked by write_u16(try_from) at write site
        .collect();

    (unique, index_map)
}

/// Write a sequence as a u32-length-prefixed byte array.
pub fn write_sequence<W: Write>(writer: &mut W, seq: &[u8]) -> Result<(), Error> {
    writer.write_u32(
        u32::try_from(seq.len())
            .map_err(|_| Error::Validation("sequence length exceeds u32::MAX".into()))?,
    )?;
    writer.write_all(seq)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::biotype::BioType;
    use crate::strand::Strand;
    use crate::transcript::types::Designations;
    use std::sync::Arc;

    fn make_tx(source: Source) -> IntermediateTranscript {
        IntermediateTranscript {
            chromosome_index: 0,
            start: 1,
            end: 100,
            id: "NM_001.1".to_string(),
            biotype: BioType::MRna,
            source,
            strand: Strand::Forward,
            gene: Arc::new(Gene {
                chromosome_index: 0,
                symbol: "TP53".to_string(),
                ncbi_gene_id: None,
                ensembl_id: None,
                hgnc_id: None,
                on_reverse_strand: false,
            }),
            transcript_regions: vec![],
            coding_region: None,
            cdna_seq: vec![],
            designations: Designations::default(),
        }
    }

    #[test]
    fn flags_source_bit() {
        let ensembl = make_tx(Source::Ensembl);
        assert_eq!(build_transcript_flags(&ensembl) & 0x01, 0);

        let refseq = make_tx(Source::RefSeq);
        assert_eq!(build_transcript_flags(&refseq) & 0x01, 1);
    }

    #[test]
    fn flags_has_coding_region() {
        let mut tx = make_tx(Source::RefSeq);
        assert_eq!(build_transcript_flags(&tx) & 0x02, 0);

        tx.coding_region = Some(crate::transcript::types::CodingRegion {
            genomic_start: 10,
            genomic_end: 90,
            cdna_start: 1,
            cdna_end: 81,
            protein_id: "NP_001.1".to_string(),
            protein_seq: vec![],
            cds_padding: 0,
            cds_offset: 0,
            protein_offset: 0,
            amino_acid_edits: None,
            slip: None,
        });
        assert_eq!(build_transcript_flags(&tx) & 0x02, 0x02);
    }

    #[test]
    fn flags_all_designations() {
        let mut tx = make_tx(Source::RefSeq);
        tx.designations.is_mane_plus_clinical = true;
        tx.designations.is_mane_select = true;
        tx.designations.is_refseq_select = true;
        tx.designations.is_ensembl_canonical = true;

        let flags = build_transcript_flags(&tx);
        assert_eq!(flags & 0x04, 0x04, "MANE Plus Clinical at bit 2");
        assert_eq!(flags & 0x08, 0x08, "MANE Select at bit 3");
        assert_eq!(flags & 0x10, 0x10, "RefSeq Select at bit 4");
        assert_eq!(flags & 0x20, 0x20, "Ensembl Canonical at bit 5");
    }

    #[test]
    fn flags_no_designations() {
        let tx = make_tx(Source::Ensembl);
        // Ensembl source, no coding region, no designations → all zeros
        assert_eq!(build_transcript_flags(&tx), 0x00);
    }
}
