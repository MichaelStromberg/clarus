use std::path::PathBuf;
use std::time::Instant;

use anyhow::Result;
use clap::Parser;
use colored::Colorize;

use clarus::cli;
use clarus::config::CacheConfig;
use clarus::download::{self, DownloadTask, filename_from_url};
use clarus::genome_assembly::GenomeAssembly;

#[derive(Parser)]
#[command(name = "create_cache", about = "Create a Clarus transcript cache file")]
struct Cli {
    /// Path to the JSON configuration file
    #[arg(short = 'c', long = "config")]
    config: PathBuf,

    /// Output data directory
    #[arg(short = 'o', long = "out")]
    out: PathBuf,
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

    // ── Summary ──────────────────────────────────────────
    cli::print_summary(start);
    Ok(())
}
