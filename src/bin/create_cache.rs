use std::io::BufReader;
use std::path::PathBuf;
use std::time::Instant;

use anyhow::Result;
use clap::Parser;
use colored::Colorize;

use clarus::cli;
use clarus::config::CacheConfig;
use clarus::download::{self, DownloadTask, filename_from_url};
use clarus::genome_assembly::GenomeAssembly;
use clarus::hgnc::HgncDatabase;
use clarus::pipeline::{self, PipelineContext};
use clarus::reference::reader::ReferenceReader;

#[derive(Parser)]
#[command(name = "create_cache", about = "Create a Clarus transcript cache file")]
struct Cli {
    /// Path to the JSON configuration file
    #[arg(short = 'c', long = "config")]
    config: PathBuf,

    /// Output data directory
    #[arg(short = 'o', long = "out")]
    out: PathBuf,

    /// Path to the reference sequence file (created by `create_ref`)
    #[arg(short = 'r', long = "ref")]
    reference: PathBuf,
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
    cli::kv("Output", &cli_args.out.display().to_string());

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

    let ref_file = pipeline::open_file(&cli_args.reference)?;
    let mut ref_reader = BufReader::new(ref_file);
    let reference = ReferenceReader::from_reader(&mut ref_reader)?;

    cli::kv("Chromosomes", &cli::num(reference.chromosomes.len()));
    cli::kv("Assembly", &reference.assembly.to_string());

    eprintln!();

    // ── File Resolution ──────────────────────────────────
    cli::section("File Resolution");

    let base_dir = std::env::temp_dir()
        .join("create_cache")
        .join(assembly.to_string());

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

    // ── Load HGNC Database ───────────────────────────────
    cli::section("HGNC Database");

    let hgnc_path = pipeline::local_path(&base_dir, "hgnc", &config.hgnc.url)?;
    let hgnc_db = HgncDatabase::from_tsv(pipeline::open_file(&hgnc_path)?)?;
    cli::kv("Entries", &cli::num(hgnc_db.len()));

    eprintln!();

    // Ensure output directory exists
    let cache_dir = cli_args.out.join("cache");
    std::fs::create_dir_all(&cache_dir)?;

    // ── Pipelines ────────────────────────────────────────
    if let Some(ref refseq) = config.refseq {
        let mut ctx = PipelineContext {
            reference: &reference,
            ref_reader: &mut ref_reader,
            hgnc_db: &hgnc_db,
            base_dir: &base_dir,
            out_dir: &cache_dir,
            assembly,
        };
        pipeline::refseq::run(&mut ctx, refseq)?;
    }
    if let Some(ref ensembl_config) = config.ensembl {
        let mut ctx = PipelineContext {
            reference: &reference,
            ref_reader: &mut ref_reader,
            hgnc_db: &hgnc_db,
            base_dir: &base_dir,
            out_dir: &cache_dir,
            assembly,
        };
        pipeline::ensembl::run(&mut ctx, ensembl_config)?;
    }

    // ── Summary ──────────────────────────────────────────
    cli::print_summary(start);
    Ok(())
}
