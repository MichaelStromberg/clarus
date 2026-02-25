//! Cache creation pipeline: orchestrates transcript construction and serialization.

pub mod ensembl;
pub mod refseq;

use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result, bail};

use crate::cache::index::write_index;
use crate::cache::writer::{DataSourceVersion, write_cache};
use crate::canonical;
use crate::chromosome::Chromosome;
use crate::cli;
use crate::evaluation;
use crate::genome_assembly::GenomeAssembly;
use crate::gff3::entry::Gff3Entry;
use crate::hgnc::HgncDatabase;
use crate::reference::reader::ReferenceReader;
use crate::sequence::{RnaSequences, reverse_complement};
use crate::strand::Strand;
use crate::transcript::types::{Gene, IntermediateTranscript};

/// Shared context for pipeline execution.
pub struct PipelineContext<'a> {
    pub reference: &'a ReferenceReader,
    pub ref_reader: &'a mut BufReader<File>,
    pub hgnc_db: &'a HgncDatabase,
    pub base_dir: &'a Path,
    pub out_dir: &'a Path,
    pub assembly: GenomeAssembly,
}

/// Open a file with a contextual error message.
pub fn open_file(path: &Path) -> Result<File> {
    File::open(path).with_context(|| format!("failed to open {}", path.display()))
}

/// Create a file for writing with a contextual error message.
fn create_file(path: &Path) -> Result<File> {
    File::create(path).with_context(|| format!("failed to create {}", path.display()))
}

/// Resolve local path for a downloaded file.
pub fn local_path(base_dir: &Path, section: &str, url: &str) -> Result<std::path::PathBuf> {
    Ok(base_dir
        .join(section)
        .join(crate::download::filename_from_url(url)?))
}

/// Inject mitochondrial cDNA sequences from the reference genome into the RNA map.
///
/// The `id_fn` closure extracts the transcript ID from each entry, allowing
/// RefSeq (name-based) and Ensembl (versioned stable ID) to share this logic.
///
/// Returns the number of sequences injected.
pub fn inject_mt_cdna(
    rna_sequences: &mut RnaSequences,
    genes: &[crate::gff3::entry::GeneRecord],
    reference: &ReferenceReader,
    ref_reader: &mut BufReader<File>,
    id_fn: impl Fn(&crate::gff3::entry::Gff3Entry) -> String,
) -> Result<usize> {
    let Some(mt_index) = reference
        .chromosomes
        .iter()
        .position(|c| c.ensembl_name == "MT")
    else {
        return Ok(0);
    };

    let mt_data = reference.load_chromosome(ref_reader, mt_index)?;
    let mut count = 0;

    for gene_record in genes {
        if gene_record.entry.chromosome_index != mt_index {
            continue;
        }

        for tx_record in &gene_record.transcripts {
            let tx_entry = &tx_record.entry;
            let tx_id = id_fn(tx_entry);

            let start = tx_entry.start;
            let end = tx_entry.end;
            let length = (end - start + 1) as usize;

            if let Some(slice) = mt_data.substring((start - 1) as i64, length) {
                let cdna = if tx_entry.strand == Strand::Reverse {
                    reverse_complement(slice)
                } else {
                    slice.to_vec()
                };

                rna_sequences.insert(tx_id, cdna);
                count += 1;
            }
        }
    }

    Ok(count)
}

/// Print designation flag counts and apply PAR reapplication.
pub fn print_designation_summary(label: &str, transcripts: &mut [IntermediateTranscript]) {
    cli::section(&format!("{label} Designation Summary"));

    let par_count = canonical::reapply_par_designations(transcripts);
    if par_count > 0 {
        cli::kv("PAR reapplied", &cli::num(par_count));
    }

    let counts = [
        (
            "MANE Select",
            transcripts
                .iter()
                .filter(|t| t.designations.is_mane_select)
                .count(),
        ),
        (
            "MANE Plus Clinical",
            transcripts
                .iter()
                .filter(|t| t.designations.is_mane_plus_clinical)
                .count(),
        ),
        (
            "RefSeq Select",
            transcripts
                .iter()
                .filter(|t| t.designations.is_refseq_select)
                .count(),
        ),
        (
            "Ensembl Canonical",
            transcripts
                .iter()
                .filter(|t| t.designations.is_ensembl_canonical)
                .count(),
        ),
    ];

    for (name, count) in counts {
        if count > 0 {
            cli::kv(name, &cli::num(count));
        }
    }
}

/// Print transcript construction statistics.
pub fn print_transcript_stats(transcripts: &[IntermediateTranscript], skipped_no_rna: u32) {
    let num_coding = transcripts
        .iter()
        .filter(|t| t.coding_region.is_some())
        .count();
    let num_non_coding = transcripts.len() - num_coding;

    cli::kv("Transcripts", &cli::num(transcripts.len()));
    cli::kv("Coding", &cli::num(num_coding));
    cli::kv("Non-coding", &cli::num(num_non_coding));
    if skipped_no_rna > 0 {
        cli::kv("Skipped (no RNA)", &cli::num(skipped_no_rna));
    }
}

/// Print evaluation results and bail if there are unresolvable transcripts.
pub fn print_eval_stats(label: &str, stats: &evaluation::EvaluationStats) -> Result<()> {
    cli::kv("Perfect protein", &cli::num(stats.num_perfect_protein));
    cli::kv("Non-coding", &cli::num(stats.num_non_coding));
    cli::kv("Contained", &cli::num(stats.num_contained));
    cli::kv("AA edits", &cli::num(stats.num_with_amino_acid_edits));
    cli::kv("Slippage", &cli::num(stats.num_translational_slippage));
    cli::kv("Wrong frame", &cli::num(stats.num_wrong_frame));
    cli::kv("Unresolvable", &cli::num(stats.num_unresolvable));

    if stats.num_unresolvable > 0 {
        bail!(
            "{} unresolvable transcript(s) in {label} pipeline — cannot write cache",
            stats.num_unresolvable
        );
    }

    Ok(())
}

/// Parse an ISO 8601 date string "YYYY-MM-DD" into (year, month, day).
pub fn parse_release_date(date_str: &str) -> Result<(u16, u8, u8)> {
    let mut parts = date_str.split('-');
    let err = || format!("invalid release date format (expected YYYY-MM-DD): '{date_str}'");

    let year: u16 = parts.next().with_context(err)?.parse().with_context(err)?;
    let month: u8 = parts.next().with_context(err)?.parse().with_context(err)?;
    let day: u8 = parts.next().with_context(err)?.parse().with_context(err)?;
    Ok((year, month, day))
}

/// Build a `Gene` from resolved fields.
pub fn build_gene(
    chromosome_index: usize,
    symbol: String,
    ncbi_gene_id: Option<String>,
    ensembl_id: Option<String>,
    hgnc_id: Option<u32>,
    on_reverse_strand: bool,
) -> Arc<Gene> {
    Arc::new(Gene {
        chromosome_index,
        symbol,
        ncbi_gene_id,
        ensembl_id,
        hgnc_id,
        on_reverse_strand,
    })
}

/// Write cache and index files for a set of transcripts.
pub fn write_cache_files(
    data_source: &DataSourceVersion<'_>,
    assembly: GenomeAssembly,
    reference_id: u32,
    chromosomes: &[Chromosome],
    transcripts: &[IntermediateTranscript],
    regulatory_regions: &[Gff3Entry],
    out_dir: &Path,
) -> Result<()> {
    cli::section("Cache Output");

    let file_prefix = format!(
        "{assembly}_{name}_{year:04}{month:02}{day:02}",
        assembly = assembly,
        name = data_source.name,
        year = data_source.release_year,
        month = data_source.release_month,
        day = data_source.release_day,
    );

    // Write transcript cache file (.cache)
    let cache_path = out_dir.join(format!("{file_prefix}.cache"));
    let mut cache_file = BufWriter::new(create_file(&cache_path)?);
    let (index_entries, cache_id) = write_cache(
        &mut cache_file,
        assembly,
        reference_id,
        data_source,
        chromosomes,
        transcripts,
        regulatory_regions,
    )?;
    drop(cache_file);

    // Write cache index file (.cache.idx)
    let idx_path = out_dir.join(format!("{file_prefix}.cache.idx"));
    let mut idx_file = BufWriter::new(create_file(&idx_path)?);
    write_index(&mut idx_file, cache_id, &index_entries)?;
    drop(idx_file);

    cli::kv("Cache ID", &format!("{cache_id:#010X}"));
    cli::kv("Transcripts", &cli::num(transcripts.len()));
    cli::kv("Regulatory", &cli::num(regulatory_regions.len()));
    eprintln!();
    cli::success(&format!("{}", cache_path.display()));
    cli::success(&format!("{}", idx_path.display()));

    eprintln!();

    Ok(())
}
