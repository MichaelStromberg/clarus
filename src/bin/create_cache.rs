use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Parser;
use colored::Colorize;

use clarus::cli;
use clarus::config::{CacheConfig, EnsemblConfig, RefSeqCacheConfig};
use clarus::download::{self, DownloadTask, filename_from_url};
use clarus::evaluation;
use clarus::genbank;
use clarus::genome_assembly::GenomeAssembly;
use clarus::gff3;
use clarus::reference::reader::ReferenceReader;
use clarus::sequence::{ProteinSequences, RnaSequences};
use clarus::strand::Strand;
use clarus::transcript::construction::{build_transcript_regions, detect_coding_region};
use clarus::transcript::types::{Gene, IntermediateTranscript, Source};

#[derive(Parser)]
#[command(name = "create_cache", about = "Create a Clarus transcript cache file")]
struct Cli {
    /// Path to the JSON configuration file
    #[arg(short = 'c', long = "config")]
    config: PathBuf,

    /// Output data directory
    #[arg(short = 'o', long = "out")]
    out: PathBuf,

    /// Path to the reference sequence file (created by create_ref)
    #[arg(short = 'r', long = "ref")]
    reference: PathBuf,
}

/// Open a file with a contextual error message.
fn open_file(path: &Path) -> Result<File> {
    File::open(path).with_context(|| format!("failed to open {}", path.display()))
}

/// Resolve local path for a downloaded file.
fn local_path(base_dir: &Path, section: &str, url: &str) -> Result<PathBuf> {
    Ok(base_dir.join(section).join(filename_from_url(url)?))
}

fn main() -> Result<()> {
    let start = Instant::now();
    let cli_args = Cli::parse();

    cli::banner("Create Cache");

    // ── Configuration ────────────────────────────────────
    cli::section("Configuration");

    let config = CacheConfig::from_file(&cli_args.config)?;
    let assembly: GenomeAssembly = config.genome_assembly.parse()?;

    cli::kv("Config", &cli_args.config.display().to_string());
    cli::kv("Assembly", &assembly.to_string());
    cli::kv("Reference", &cli_args.reference.display().to_string());

    if let Some(ref ensembl) = config.ensembl {
        cli::kv(
            "Ensembl",
            &format!("v{} ({})", ensembl.version, ensembl.release_date),
        );
    }
    if let Some(ref refseq) = config.refseq {
        cli::kv(
            "RefSeq",
            &format!("{} ({})", refseq.version, refseq.release_date),
        );
    }

    eprintln!();

    // ── Load Reference ───────────────────────────────────
    cli::section("Reference");

    let ref_file = open_file(&cli_args.reference)?;
    let mut ref_reader = BufReader::new(ref_file);
    let reference = ReferenceReader::from_reader(&mut ref_reader)?;

    cli::kv("Chromosomes", &reference.chromosomes.len().to_string());
    cli::kv("Assembly", &reference.assembly.to_string());

    eprintln!();

    // ── File Resolution ──────────────────────────────────
    cli::section("File Resolution");

    let base_dir = PathBuf::from("/tmp/create_cache").join(assembly.to_string());

    // Build download tasks from file_entries (ensembl + refseq)
    let mut all_tasks: Vec<DownloadTask> = Vec::new();

    for (section_name, field_name, entry) in config.file_entries() {
        let filename = filename_from_url(&entry.url)?;
        let dest = base_dir.join(section_name).join(&filename);
        all_tasks.push(DownloadTask {
            name: format!("{section_name}/{field_name}"),
            url: entry.url.clone(),
            dest,
            expected_md5: Some(entry.md5.clone()),
        });
    }

    // Add HGNC separately (no MD5)
    {
        let filename = filename_from_url(&config.hgnc.url)?;
        let dest = base_dir.join("hgnc").join(&filename);
        all_tasks.push(DownloadTask {
            name: "hgnc".to_string(),
            url: config.hgnc.url.clone(),
            dest,
            expected_md5: None,
        });
    }

    // Display cached/missing status
    for task in &all_tasks {
        let filename = task.dest.file_name().unwrap_or_default().to_string_lossy();
        if task.dest.exists() {
            cli::kv(&task.name, &format!("{filename} {}", "(cached)".dimmed()));
        } else {
            cli::kv(&task.name, &format!("{filename} {}", "(missing)".yellow()));
        }
    }

    let to_download = download::resolve_tasks(&all_tasks);

    eprintln!();

    // ── Downloading ──────────────────────────────────────
    if !to_download.is_empty() {
        cli::section("Downloading");
        download::download_and_verify(&to_download)?;
        eprintln!();
    }

    // ── Pipelines ────────────────────────────────────────
    if let Some(ref refseq) = config.refseq {
        run_refseq_pipeline(refseq, &reference, &base_dir)?;
    }
    if let Some(ref ensembl_config) = config.ensembl {
        run_ensembl_pipeline(ensembl_config, &config.hgnc.url, &reference, &base_dir)?;
    }

    // ── Summary ──────────────────────────────────────────
    cli::print_summary(start);
    Ok(())
}

fn run_refseq_pipeline(
    refseq: &RefSeqCacheConfig,
    reference: &ReferenceReader,
    base_dir: &Path,
) -> Result<()> {
    cli::section("RefSeq Parsing");

    let gff3_path = local_path(base_dir, "refseq", &refseq.gff3.url)?;
    let genbank_path = local_path(base_dir, "refseq", &refseq.genbank.url)?;
    let rna_path = local_path(base_dir, "refseq", &refseq.rna_fasta.url)?;
    let protein_path = local_path(base_dir, "refseq", &refseq.protein_fasta.url)?;

    // 1. Parse GFF3
    let gff3_result = gff3::parse_refseq_gff3_gz(open_file(&gff3_path)?, &reference.name_to_index)?;
    cli::kv(
        "GFF3",
        &format!(
            "{} genes, {} regulatory regions",
            gff3_result.genes.len(),
            gff3_result.regulatory_regions.len()
        ),
    );

    // 2. Parse GenBank
    let genbank_records = genbank::parse_genbank_gz(open_file(&genbank_path)?)?;
    cli::kv("GenBank", &format!("{} CDS records", genbank_records.len()));

    // 3. Parse RNA FASTA
    let mut rna_sequences = RnaSequences::from_gz(open_file(&rna_path)?)?;
    cli::kv("RNA FASTA", &format!("{} sequences", rna_sequences.len()));

    // 4. Parse Protein FASTA
    let protein_sequences = ProteinSequences::from_gz(open_file(&protein_path)?)?;
    cli::kv(
        "Protein FASTA",
        &format!("{} sequences", protein_sequences.len()),
    );

    eprintln!();

    // ── Transcript Construction ──────────────────────────
    cli::section("Transcript Construction");

    let mut transcripts: Vec<IntermediateTranscript> = Vec::new();
    let mut skipped_no_rna: u32 = 0;

    for gene_record in &gff3_result.genes {
        let gene_entry = &gene_record.entry;
        let on_reverse_strand = gene_entry.strand == Strand::Reverse;

        let gene = Arc::new(Gene {
            chromosome_index: gene_entry.chromosome_index,
            symbol: gene_entry
                .attributes
                .gene_symbol
                .clone()
                .unwrap_or_default(),
            ncbi_gene_id: gene_entry.attributes.gene_id.clone(),
            ensembl_id: None,
            hgnc_id: gene_entry.attributes.hgnc_id,
            on_reverse_strand,
        });

        for tx_record in &gene_record.transcripts {
            let tx_entry = &tx_record.entry;
            let tx_name = tx_entry
                .attributes
                .name
                .as_deref()
                .unwrap_or(&tx_entry.attributes.id);

            // Take ownership of cDNA sequence to avoid copying.
            // Each tx_name appears at most once across all GFF3 gene→transcript entries.
            let cdna_seq = match rna_sequences.remove(tx_name) {
                Some(seq) => seq,
                None => {
                    skipped_no_rna += 1;
                    continue;
                }
            };

            // Build transcript regions (exon normalization + intron insertion)
            let regions =
                build_transcript_regions(&tx_record.exons, &tx_record.matches, tx_entry.strand)?;

            // Detect coding region if CDS entries exist
            let is_coding = tx_entry.biotype.is_coding() && !tx_record.cds_entries.is_empty();
            let coding_region = if is_coding {
                // Extract protein ID from first CDS entry's Name attribute (RefSeq),
                // falling back to FASTA cross-reference (Ensembl)
                let protein_id = tx_record
                    .cds_entries
                    .first()
                    .and_then(|cds| cds.attributes.name.as_deref())
                    .or_else(|| protein_sequences.protein_id_for_transcript(tx_name))
                    .unwrap_or("");

                if protein_id.is_empty() {
                    None
                } else {
                    let protein_seq = protein_sequences
                        .get_by_protein_id(protein_id)
                        .map(|s| s.to_vec())
                        .unwrap_or_default();

                    let genbank_cds = genbank_records.get(tx_name);

                    match detect_coding_region(
                        &tx_record.cds_entries,
                        &regions,
                        tx_entry.strand,
                        protein_id,
                        protein_seq,
                        genbank_cds,
                    ) {
                        Ok(cr) => Some(cr),
                        Err(e) => {
                            cli::warning(&format!(
                                "coding region detection failed for {tx_name}: {e}"
                            ));
                            None
                        }
                    }
                }
            } else {
                None
            };

            transcripts.push(IntermediateTranscript {
                chromosome_index: tx_entry.chromosome_index,
                start: tx_entry.start,
                end: tx_entry.end,
                id: tx_name.to_string(),
                biotype: tx_entry.biotype,
                source: Source::RefSeq,
                strand: tx_entry.strand,
                gene: Arc::clone(&gene),
                transcript_regions: regions,
                coding_region,
                cdna_seq,
                is_canonical: tx_entry.attributes.is_refseq_select,
                is_mane_select: tx_entry.attributes.is_mane_select,
            });
        }
    }

    let num_coding = transcripts
        .iter()
        .filter(|t| t.coding_region.is_some())
        .count();
    let num_non_coding = transcripts.len() - num_coding;

    cli::kv("Transcripts", &transcripts.len().to_string());
    cli::kv("Coding", &num_coding.to_string());
    cli::kv("Non-coding", &num_non_coding.to_string());
    if skipped_no_rna > 0 {
        cli::kv("Skipped (no RNA)", &skipped_no_rna.to_string());
    }

    eprintln!();

    // ── Evaluation ───────────────────────────────────────
    cli::section("Evaluation");

    let stats = evaluation::evaluate_transcripts(&mut transcripts, &reference.chromosomes)?;

    cli::kv("Perfect protein", &stats.num_perfect_protein.to_string());
    cli::kv("Non-coding", &stats.num_non_coding.to_string());
    cli::kv("Contained", &stats.num_contained.to_string());
    cli::kv("AA edits", &stats.num_with_amino_acid_edits.to_string());
    cli::kv("Slippage", &stats.num_translational_slippage.to_string());
    cli::kv("Wrong frame", &stats.num_wrong_frame.to_string());
    cli::kv("Unresolvable", &stats.num_unresolvable.to_string());

    eprintln!();
    Ok(())
}

fn run_ensembl_pipeline(
    ensembl_config: &EnsemblConfig,
    hgnc_url: &str,
    reference: &ReferenceReader,
    base_dir: &Path,
) -> Result<()> {
    cli::section("Ensembl Parsing");

    // Parse HGNC for Ensembl gene ID → HGNC ID lookup (selenoprotein detection)
    let hgnc_path = local_path(base_dir, "hgnc", hgnc_url)?;
    let ensembl_to_hgnc = clarus::hgnc::parse_hgnc_tsv(open_file(&hgnc_path)?)?;

    let gff3_path = local_path(base_dir, "ensembl", &ensembl_config.gff3.url)?;
    let reg_gff_path = local_path(base_dir, "ensembl", &ensembl_config.regulatory_gff.url)?;
    let cdna_path = local_path(base_dir, "ensembl", &ensembl_config.cdna_fasta.url)?;
    let ncrna_path = local_path(base_dir, "ensembl", &ensembl_config.ncrna_fasta.url)?;
    let peptide_path = local_path(base_dir, "ensembl", &ensembl_config.peptide_fasta.url)?;

    // 1. Parse main GFF3
    let gff3_result =
        gff3::parse_ensembl_gff3_gz(open_file(&gff3_path)?, &reference.name_to_index)?;
    cli::kv(
        "GFF3",
        &format!(
            "{} genes, {} regulatory regions",
            gff3_result.genes.len(),
            gff3_result.regulatory_regions.len()
        ),
    );

    // 2. Parse regulatory GFF3
    let regulatory_regions =
        gff3::parse_regulatory_gff3_gz(open_file(&reg_gff_path)?, &reference.name_to_index)?;
    cli::kv(
        "Regulatory GFF",
        &format!("{} regions", regulatory_regions.len()),
    );

    // 3. Parse RNA FASTA (merge cdna + ncrna)
    let mut rna_sequences = RnaSequences::from_gz(open_file(&cdna_path)?)?;
    let cdna_count = rna_sequences.len();

    let ncrna_sequences = RnaSequences::from_gz(open_file(&ncrna_path)?)?;
    let ncrna_count = ncrna_sequences.len();
    rna_sequences.merge(ncrna_sequences)?;
    cli::kv(
        "RNA FASTA",
        &format!(
            "{} sequences (cdna: {}, ncrna: {})",
            rna_sequences.len(),
            cdna_count,
            ncrna_count
        ),
    );

    // 4. Parse Protein FASTA
    let protein_sequences = ProteinSequences::from_gz(open_file(&peptide_path)?)?;
    cli::kv(
        "Protein FASTA",
        &format!("{} sequences", protein_sequences.len()),
    );

    eprintln!();

    // ── Transcript Construction ──────────────────────────
    cli::section("Ensembl Transcript Construction");

    let mut transcripts: Vec<IntermediateTranscript> = Vec::new();
    let mut skipped_no_rna: u32 = 0;

    for gene_record in &gff3_result.genes {
        let gene_entry = &gene_record.entry;
        let on_reverse_strand = gene_entry.strand == Strand::Reverse;

        let gene = Arc::new(Gene {
            chromosome_index: gene_entry.chromosome_index,
            symbol: gene_entry
                .attributes
                .name
                .clone()
                .or_else(|| gene_entry.attributes.gene_id.clone())
                .unwrap_or_default(),
            ncbi_gene_id: None,
            ensembl_id: gene_entry.attributes.gene_id.clone(),
            hgnc_id: gene_entry
                .attributes
                .gene_id
                .as_deref()
                .and_then(|eid| ensembl_to_hgnc.get(eid).copied()),
            on_reverse_strand,
        });

        for tx_record in &gene_record.transcripts {
            let tx_entry = &tx_record.entry;

            // Build versioned transcript ID: stable_id.version
            let tx_stable_id = tx_entry
                .attributes
                .transcript_id
                .as_deref()
                .unwrap_or_else(|| {
                    // Fallback: extract from ID attr after colon
                    tx_entry
                        .attributes
                        .id
                        .find(':')
                        .map(|pos| &tx_entry.attributes.id[pos + 1..])
                        .unwrap_or(&tx_entry.attributes.id)
                });
            let version = tx_entry.attributes.version.unwrap_or(1);
            let versioned_id = format!("{tx_stable_id}.{version}");

            // Look up cDNA sequence by versioned ID
            let cdna_seq = match rna_sequences.remove(&versioned_id) {
                Some(seq) => seq,
                None => {
                    skipped_no_rna += 1;
                    continue;
                }
            };

            // Build transcript regions (simple path — no cDNA_match in Ensembl)
            let regions =
                build_transcript_regions(&tx_record.exons, &tx_record.matches, tx_entry.strand)?;

            // Detect coding region if CDS entries exist
            let is_coding = tx_entry.biotype.is_coding() && !tx_record.cds_entries.is_empty();
            let coding_region = if is_coding {
                // Protein ID from FASTA cross-reference, falling back to CDS protein_id attr
                let protein_id = protein_sequences
                    .protein_id_for_transcript(&versioned_id)
                    .or_else(|| {
                        tx_record
                            .cds_entries
                            .first()
                            .and_then(|cds| cds.attributes.protein_id.as_deref())
                    })
                    .unwrap_or("");

                if protein_id.is_empty() {
                    None
                } else {
                    let protein_seq = protein_sequences
                        .get_by_protein_id(protein_id)
                        .map(|s| s.to_vec())
                        .unwrap_or_default();

                    // No GenBank CDS data for Ensembl (CDS offset always 0)
                    match detect_coding_region(
                        &tx_record.cds_entries,
                        &regions,
                        tx_entry.strand,
                        protein_id,
                        protein_seq,
                        None,
                    ) {
                        Ok(cr) => Some(cr),
                        Err(e) => {
                            cli::warning(&format!(
                                "coding region detection failed for {versioned_id}: {e}"
                            ));
                            None
                        }
                    }
                }
            } else {
                None
            };

            transcripts.push(IntermediateTranscript {
                chromosome_index: tx_entry.chromosome_index,
                start: tx_entry.start,
                end: tx_entry.end,
                id: versioned_id,
                biotype: tx_entry.biotype,
                source: Source::Ensembl,
                strand: tx_entry.strand,
                gene: Arc::clone(&gene),
                transcript_regions: regions,
                coding_region,
                cdna_seq,
                is_canonical: tx_entry.attributes.is_ensembl_canonical,
                is_mane_select: false,
            });
        }
    }

    let num_coding = transcripts
        .iter()
        .filter(|t| t.coding_region.is_some())
        .count();
    let num_non_coding = transcripts.len() - num_coding;

    cli::kv("Transcripts", &transcripts.len().to_string());
    cli::kv("Coding", &num_coding.to_string());
    cli::kv("Non-coding", &num_non_coding.to_string());
    if skipped_no_rna > 0 {
        cli::kv("Skipped (no RNA)", &skipped_no_rna.to_string());
    }

    eprintln!();

    // ── Evaluation ───────────────────────────────────────
    cli::section("Ensembl Evaluation");

    let stats = evaluation::evaluate_transcripts(&mut transcripts, &reference.chromosomes)?;

    cli::kv("Perfect protein", &stats.num_perfect_protein.to_string());
    cli::kv("Non-coding", &stats.num_non_coding.to_string());
    cli::kv("Contained", &stats.num_contained.to_string());
    cli::kv("AA edits", &stats.num_with_amino_acid_edits.to_string());
    cli::kv("Slippage", &stats.num_translational_slippage.to_string());
    cli::kv("Wrong frame", &stats.num_wrong_frame.to_string());
    cli::kv("Unresolvable", &stats.num_unresolvable.to_string());

    eprintln!();
    Ok(())
}
