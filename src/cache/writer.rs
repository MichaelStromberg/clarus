//! Cache file writer: produces the transcript cache binary file.
//!
//! Layout: Common header → cache header → per-reference compressed blocks.
//! Each reference block contains shared dedup arrays + sparse bin records.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::{Seek, Write};

use crate::cache::bin::bin_index;
use crate::cache::serialize::{
    dedup_genes, dedup_regions, dedup_sequences, write_gene, write_regulatory_region,
    write_sequence, write_transcript, write_transcript_region,
};
use crate::cache::{CACHE_FILE_TYPE, CACHE_FORMAT_VERSION, IndexEntry};
use crate::chromosome::Chromosome;
use crate::error::Error;
use crate::genome_assembly::GenomeAssembly;
use crate::gff3::entry::Gff3Entry;
use crate::reference::binary_io::BinaryWrite;
use crate::reference::common_header::write_common_header;
use crate::transcript::types::IntermediateTranscript;

/// Data source version information for the cache header.
pub struct DataSourceVersion<'a> {
    pub name: &'a str,
    pub version: &'a str,
    pub release_year: u16,
    pub release_month: u8,
    pub release_day: u8,
}

/// Compute a deterministic cache ID from the data source version string.
///
/// The ID is the first 4 bytes (LE) of the SHA-256 hash of the version string.
pub fn compute_cache_id(version: &str) -> u32 {
    use sha2::{Digest, Sha256};
    let hash = Sha256::digest(version.as_bytes());
    u32::from_le_bytes([hash[0], hash[1], hash[2], hash[3]])
}

/// Write a complete transcript cache file.
///
/// Returns a list of `(ref_index, byte_offset)` entries for the index file,
/// along with the computed cache_id.
pub fn write_cache<W: Write + Seek>(
    writer: &mut W,
    assembly: GenomeAssembly,
    reference_id: u32,
    data_source: &DataSourceVersion<'_>,
    chromosomes: &[Chromosome],
    transcripts: &[IntermediateTranscript],
    regulatory_regions: &[Gff3Entry],
) -> Result<(Vec<IndexEntry>, u32), Error> {
    let cache_id = compute_cache_id(data_source.version);

    // ── Pre-scan: find chromosomes with data ────────────
    let chroms_with_data: BTreeSet<usize> = transcripts
        .iter()
        .map(|tx| tx.chromosome_index)
        .chain(regulatory_regions.iter().map(|r| r.chromosome_index))
        .collect();
    let chromosome_count = u16::try_from(chroms_with_data.len())
        .map_err(|_| Error::Validation("chromosome count exceeds u16::MAX".into()))?;

    // ── Header (spec 10 §3.1–3.2) ──────────────────────
    write_common_header(writer, CACHE_FILE_TYPE, CACHE_FORMAT_VERSION, 0)?; // 0 = placeholder
    writer.write_u8(assembly.to_byte())?;
    writer.write_u32(reference_id)?;
    writer.write_u32(cache_id)?;
    writer.write_u16(chromosome_count)?;
    writer.write_prefixed_string(data_source.name)?;
    writer.write_prefixed_string(data_source.version)?;
    writer.write_u16(data_source.release_year)?;
    writer.write_u8(data_source.release_month)?;
    writer.write_u8(data_source.release_day)?;

    // ── Group transcripts and regulatory regions by chromosome ──
    let mut tx_by_chrom: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, tx) in transcripts.iter().enumerate() {
        tx_by_chrom.entry(tx.chromosome_index).or_default().push(i);
    }

    let mut reg_by_chrom: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, reg) in regulatory_regions.iter().enumerate() {
        reg_by_chrom
            .entry(reg.chromosome_index)
            .or_default()
            .push(i);
    }

    // ── Body: per-reference compressed blocks ───────────
    let mut index_entries: Vec<IndexEntry> = Vec::new();

    for &chrom_idx in &chroms_with_data {
        let ref_offset = writer.stream_position()?;
        let ref_index = chromosomes[chrom_idx].ref_index;

        // Build uncompressed reference buffer
        let ref_data = serialize_reference(
            chrom_idx,
            transcripts,
            regulatory_regions,
            &tx_by_chrom,
            &reg_by_chrom,
        )?;

        // Compress with zstd level 21
        let compressed = zstd::encode_all(ref_data.as_slice(), 21)?;

        // Write: ref_index(u16) + decompressed_size(u32) + compressed_size(u32) + data
        writer.write_u16(ref_index)?;
        writer.write_u32(ref_data.len() as u32)?;
        writer.write_u32(compressed.len() as u32)?;
        writer.write_all(&compressed)?;

        index_entries.push((ref_index, ref_offset));
    }

    // ── Patch file_length at offset 12 ──────────────────
    let file_length = writer.stream_position()?;
    writer.seek(std::io::SeekFrom::Start(12))?;
    writer.write_u64(file_length)?;
    writer.seek(std::io::SeekFrom::End(0))?;

    Ok((index_entries, cache_id))
}

/// Serialize a single reference's uncompressed data.
///
/// Contains shared dedup arrays followed by sparse bin records.
fn serialize_reference(
    chrom_idx: usize,
    transcripts: &[IntermediateTranscript],
    regulatory_regions: &[Gff3Entry],
    tx_by_chrom: &HashMap<usize, Vec<usize>>,
    reg_by_chrom: &HashMap<usize, Vec<usize>>,
) -> Result<Vec<u8>, Error> {
    let mut buf = Vec::new();

    // Collect all transcripts for this chromosome, sorted by (start, end)
    let tx_indices = tx_by_chrom.get(&chrom_idx).cloned().unwrap_or_default();
    let mut chrom_txs: Vec<&IntermediateTranscript> =
        tx_indices.iter().map(|&i| &transcripts[i]).collect();
    chrom_txs.sort_by(|a, b| a.start.cmp(&b.start).then(a.end.cmp(&b.end)));

    // ── Shared dedup arrays (spec 10 §4.1) ──────────────
    let (genes, gene_index) = dedup_genes(&chrom_txs);
    let (regions, region_index) = dedup_regions(&chrom_txs);

    let cdna_seqs: Vec<&[u8]> = chrom_txs.iter().map(|tx| tx.cdna_seq.as_slice()).collect();
    let (unique_cdna, cdna_index) = dedup_sequences(cdna_seqs);

    let protein_seqs: Vec<&[u8]> = chrom_txs
        .iter()
        .filter_map(|tx| {
            tx.coding_region
                .as_ref()
                .map(|cr| cr.protein_seq.as_slice())
        })
        .collect();
    let (unique_proteins, protein_index) = dedup_sequences(protein_seqs);

    // Write gene count + genes
    buf.write_u16(
        u16::try_from(genes.len())
            .map_err(|_| Error::Validation("gene count exceeds u16::MAX".into()))?,
    )?;
    for gene in &genes {
        write_gene(&mut buf, gene)?;
    }

    // Write transcript region count + regions
    buf.write_u32(regions.len() as u32)?;
    for region in &regions {
        write_transcript_region(&mut buf, region)?;
    }

    // Write cDNA sequence count + sequences
    buf.write_u16(
        u16::try_from(unique_cdna.len())
            .map_err(|_| Error::Validation("cDNA sequence count exceeds u16::MAX".into()))?,
    )?;
    for seq in &unique_cdna {
        write_sequence(&mut buf, seq)?;
    }

    // Write protein sequence count + sequences
    buf.write_u16(
        u16::try_from(unique_proteins.len())
            .map_err(|_| Error::Validation("protein sequence count exceeds u16::MAX".into()))?,
    )?;
    for seq in &unique_proteins {
        write_sequence(&mut buf, seq)?;
    }

    // ── Group into bins (spec 10 §1.4: store in every bin the feature overlaps) ──
    // Transcripts by bin
    let mut tx_by_bin: BTreeMap<u8, Vec<&IntermediateTranscript>> = BTreeMap::new();
    for tx in &chrom_txs {
        let start_bin = bin_index(tx.start);
        let end_bin = bin_index(tx.end);
        for bin in start_bin..=end_bin {
            tx_by_bin.entry(bin).or_default().push(tx);
        }
    }

    // Regulatory regions by bin
    let reg_indices = reg_by_chrom.get(&chrom_idx).cloned().unwrap_or_default();
    let mut reg_by_bin: BTreeMap<u8, Vec<&Gff3Entry>> = BTreeMap::new();
    for &i in &reg_indices {
        let reg = &regulatory_regions[i];
        let start_bin = bin_index(reg.start);
        let end_bin = bin_index(reg.end);
        for bin in start_bin..=end_bin {
            reg_by_bin.entry(bin).or_default().push(reg);
        }
    }

    // Collect all bins with data (sparse)
    let bins_with_data: BTreeSet<u8> = tx_by_bin.keys().chain(reg_by_bin.keys()).copied().collect();

    // Write num bins with data
    buf.write_u8(
        u8::try_from(bins_with_data.len())
            .map_err(|_| Error::Validation("bin count exceeds u8::MAX".into()))?,
    )?;

    // ── Sparse bin records (spec 10 §4.2) ───────────────
    for &bin_num in &bins_with_data {
        buf.write_u8(bin_num)?;

        // Transcripts in this bin
        let bin_txs: &[_] = tx_by_bin.get(&bin_num).map_or(&[], Vec::as_slice);
        buf.write_u16(
            u16::try_from(bin_txs.len())
                .map_err(|_| Error::Validation("bin transcript count exceeds u16::MAX".into()))?,
        )?;

        for tx in bin_txs {
            let gene_idx = gene_index[tx.gene.as_ref()];
            let region_idxs: Vec<u32> = tx
                .transcript_regions
                .iter()
                .map(|r| region_index[r])
                .collect();
            let cdna_idx = cdna_index[&tx.cdna_seq];
            let protein_idx = tx
                .coding_region
                .as_ref()
                .map_or(0, |cr| protein_index[&cr.protein_seq]);

            write_transcript(&mut buf, tx, gene_idx, &region_idxs, cdna_idx, protein_idx)?;
        }

        // Regulatory regions in this bin
        let bin_regs: &[_] = reg_by_bin.get(&bin_num).map_or(&[], Vec::as_slice);
        buf.write_u16(u16::try_from(bin_regs.len()).map_err(|_| {
            Error::Validation("bin regulatory region count exceeds u16::MAX".into())
        })?)?;

        for reg in bin_regs {
            write_regulatory_region(&mut buf, reg)?;
        }
    }

    Ok(buf)
}
