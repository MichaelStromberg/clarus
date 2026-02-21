use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result, bail};
use clap::Parser;
use colored::Colorize;

use clarus::assembly_report::parse_assembly_report;
use clarus::bands::{Band, parse_ideogram};
use clarus::chromosome::Chromosome;
use clarus::fasta::parse_fasta_gz;
use clarus::genome_assembly::{GenomeAssembly, compute_reference_id};
use clarus::mirna::{MirnaRegion, parse_mirna_gff3};
use clarus::perf;
use clarus::reference::reader::ReferenceReader;
use clarus::reference::writer::ReferenceWriter;

/// Maximum number of N bases allowed in chrY before flagging PAR masking.
const CHRY_N_MASKING_THRESHOLD: usize = 33_720_001;

#[derive(Parser)]
#[command(name = "create_ref", about = "Create a Clarus reference sequence file")]
struct Cli {
    /// Path to the NCBI genome assembly report file
    #[arg(short = 'a', long = "assembly-report")]
    assembly_report: PathBuf,

    /// Path to the gzip-compressed FASTA file
    #[arg(short = 'f', long = "fasta")]
    fasta: PathBuf,

    /// Path to the NCBI GDP ideogram file (cytogenetic bands)
    #[arg(short = 'b', long = "bands")]
    bands: Option<PathBuf>,

    /// Path to a gzip-compressed GFF3 file (for miRNA regions)
    #[arg(long = "gff")]
    gff: Option<PathBuf>,

    /// Output reference file path
    #[arg(short = 'o', long = "out")]
    out: PathBuf,

    /// Genome assembly name (GRCh37 or GRCh38)
    #[arg(short = 'g', long = "assembly", default_value = "GRCh38")]
    assembly: String,

    /// Assembly patch level
    #[arg(short = 'p', long = "patch-level", default_value_t = 14)]
    patch_level: u8,
}

fn banner() {
    eprintln!();
    eprintln!("{} {}", "Clarus".bold().cyan(), "Create Reference".dimmed());
    eprintln!("{}", "(c) 2026 Michael Stromberg".dimmed());
    eprintln!();
}

fn section(title: &str) {
    let bar = "─".repeat(50);
    eprintln!("{} {}", title.bold().blue(), bar.dimmed());
}

fn kv(key: &str, value: &str) {
    eprintln!("  {:<20} {}", key.dimmed(), value);
}

fn success(msg: &str) {
    eprintln!("  {} {}", "✓".green().bold(), msg);
}

fn warning(msg: &str) {
    eprintln!("  {} {}", "⚠".yellow(), msg.yellow());
}

fn main() -> Result<()> {
    let start = Instant::now();
    let cli = Cli::parse();

    let assembly: GenomeAssembly = cli.assembly.parse()?;

    banner();

    // ── Parsing ──────────────────────────────────────────
    section("Parsing");

    let report_file = File::open(&cli.assembly_report).context("failed to open assembly report")?;
    let report_reader = BufReader::new(report_file);
    let (chromosomes, name_to_index) = parse_assembly_report(report_reader)?;

    kv(
        "Assembly report",
        &format!(
            "{} ({} chromosomes)",
            cli.assembly_report.display(),
            chromosomes.len()
        ),
    );

    // Validate ref_index contiguity
    for (i, chr) in chromosomes.iter().enumerate() {
        if chr.ref_index as usize != i {
            bail!(
                "ref_index gap: expected {i}, got {} for chromosome {}",
                chr.ref_index,
                chr.ucsc_name
            );
        }
    }

    let fasta_file = File::open(&cli.fasta).context("failed to open FASTA file")?;
    let fasta_records = parse_fasta_gz(fasta_file)?;

    kv(
        "FASTA",
        &format!(
            "{} ({} sequences)",
            cli.fasta.display(),
            fasta_records.len()
        ),
    );

    // Parse bands if provided
    let bands_map: HashMap<String, Vec<Band>> = if let Some(ref bands_path) = cli.bands {
        let bands_file = File::open(bands_path).context("failed to open bands file")?;
        let bands_reader = BufReader::new(bands_file);
        let parsed = parse_ideogram(bands_reader)?;
        let total_bands: usize = parsed.values().map(|v| v.len()).sum();
        kv(
            "Bands",
            &format!(
                "{} ({} bands across {} chromosomes)",
                bands_path.display(),
                total_bands,
                parsed.len()
            ),
        );
        parsed
    } else {
        HashMap::new()
    };

    // Parse miRNA regions if provided
    let mirna_map: HashMap<String, Vec<MirnaRegion>> = if let Some(ref gff_path) = cli.gff {
        let gff_file = File::open(gff_path).context("failed to open GFF3 file")?;
        let parsed = parse_mirna_gff3(gff_file)?;
        let mut by_ensembl: HashMap<String, Vec<MirnaRegion>> = HashMap::new();
        for (accession, regions) in parsed {
            if let Some(&idx) = name_to_index.get(&accession) {
                by_ensembl.insert(chromosomes[idx].ensembl_name.clone(), regions);
            }
        }
        let total: usize = by_ensembl.values().map(|v| v.len()).sum();
        kv(
            "miRNA",
            &format!(
                "{} ({} regions across {} chromosomes)",
                gff_path.display(),
                total,
                by_ensembl.len()
            ),
        );
        by_ensembl
    } else {
        HashMap::new()
    };

    // Match FASTA sequences to assembly report chromosomes
    let mut sequence_map: HashMap<usize, Vec<u8>> = HashMap::new();
    let mut matched = 0usize;
    let mut skipped = 0usize;

    for (accession, seq) in fasta_records {
        if let Some(&idx) = name_to_index.get(&accession) {
            let expected_len = chromosomes[idx].length as usize;
            if seq.len() != expected_len {
                bail!(
                    "sequence length mismatch for {}: FASTA has {} bases, assembly report says {}",
                    accession,
                    seq.len(),
                    expected_len
                );
            }
            sequence_map.insert(idx, seq);
            matched += 1;
        } else {
            warning(&format!("'{accession}' not in assembly report, skipping"));
            skipped += 1;
        }
    }

    kv("Matched", &format!("{matched} sequences"));
    if skipped > 0 {
        kv("Skipped", &format!("{skipped} sequences"));
    }

    eprintln!();

    // ── Validation ───────────────────────────────────────
    section("Validation");

    if let Some(&chry_idx) = name_to_index.get("chrY")
        && let Some(chry_seq) = sequence_map.get(&chry_idx)
    {
        let n_count = chry_seq.iter().filter(|&&b| b == b'N').count();
        if n_count > CHRY_N_MASKING_THRESHOLD {
            bail!(
                "chrY has {n_count} N bases (threshold: {CHRY_N_MASKING_THRESHOLD}). \
                 The FASTA appears to have PAR masking applied. \
                 Use an unmasked reference genome."
            );
        }
        success(&format!(
            "chrY N-masking ({} N bases, threshold: {})",
            n_count.to_string().bold(),
            CHRY_N_MASKING_THRESHOLD
        ));
    }

    // Build ordered chromosome and sequence lists (only those with sequences)
    let mut output_chroms: Vec<Chromosome> = Vec::new();
    let mut output_sequences: Vec<Vec<u8>> = Vec::new();
    let mut missing = 0usize;

    for (i, mut chr) in chromosomes.into_iter().enumerate() {
        if let Some(seq) = sequence_map.remove(&i) {
            chr.ref_index = u16::try_from(output_chroms.len())
                .context("too many chromosomes for u16 ref_index")?;
            output_chroms.push(chr);
            output_sequences.push(seq);
        } else {
            missing += 1;
        }
    }

    if missing > 0 {
        warning(&format!(
            "{missing} assembly report entries have no FASTA sequence"
        ));
    }

    eprintln!();

    // ── Writing ──────────────────────────────────────────
    section("Writing");
    let reference_id = compute_reference_id(assembly, cli.patch_level);
    kv("Output", &cli.out.display().to_string());
    kv("Assembly", &format!("{assembly}.p{}", cli.patch_level));
    kv("Reference ID", &format!("0x{reference_id:08X}"));
    kv("Chromosomes", &output_chroms.len().to_string());
    if !mirna_map.is_empty() {
        let total_mirna: usize = mirna_map.values().map(|v| v.len()).sum();
        kv("miRNA regions", &total_mirna.to_string());
    }

    let mut out_file = File::create(&cli.out).context("failed to create output file")?;
    ReferenceWriter::write(
        &mut out_file,
        assembly,
        cli.patch_level,
        &output_chroms,
        &output_sequences,
        &bands_map,
        &mirna_map,
    )?;

    let file_size = std::fs::metadata(&cli.out)
        .map(|m| perf::format_bytes(m.len()))
        .unwrap_or_default();
    kv("File size", &file_size);

    eprintln!();

    // ── Verification ─────────────────────────────────────
    section("Verification");
    verify_output(
        &cli.out,
        &output_chroms,
        reference_id,
        &bands_map,
        &mirna_map,
    )?;

    // ── Summary ──────────────────────────────────────────
    let elapsed = start.elapsed();
    eprintln!();
    eprintln!(
        "{}  {}\n{}  {}",
        "Time".dimmed(),
        perf::format_elapsed(elapsed).bold(),
        "Peak memory".dimmed(),
        perf::peak_memory_bytes()
            .map(perf::format_bytes)
            .unwrap_or_else(|| "N/A".to_string())
            .bold(),
    );
    eprintln!();
    Ok(())
}

fn verify_output(
    path: &Path,
    expected_chroms: &[Chromosome],
    expected_reference_id: u32,
    bands_map: &HashMap<String, Vec<Band>>,
    mirna_map: &HashMap<String, Vec<MirnaRegion>>,
) -> Result<()> {
    let verify_file = File::open(path).context("failed to open output for verification")?;
    let mut verify_reader = BufReader::new(verify_file);
    let ref_reader = ReferenceReader::from_reader(&mut verify_reader)?;

    if ref_reader.chromosomes.len() != expected_chroms.len() {
        bail!(
            "verification failed: expected {} chromosomes, got {}",
            expected_chroms.len(),
            ref_reader.chromosomes.len()
        );
    }

    if ref_reader.reference_id != expected_reference_id {
        bail!(
            "verification failed: expected reference ID 0x{expected_reference_id:08X}, got 0x{:08X}",
            ref_reader.reference_id
        );
    }

    for (i, chr) in ref_reader.chromosomes.iter().enumerate() {
        if chr.length != expected_chroms[i].length {
            bail!(
                "verification failed: chromosome {} length mismatch ({} vs {})",
                chr.ucsc_name,
                chr.length,
                expected_chroms[i].length
            );
        }
    }

    // Verify band and miRNA counts per chromosome
    let has_bands = !bands_map.is_empty();
    let has_mirna = !mirna_map.is_empty();
    if has_bands || has_mirna {
        let mut total_bands: usize = 0;
        let mut total_mirna: usize = 0;
        for (i, chr) in ref_reader.chromosomes.iter().enumerate() {
            let chr_data = ref_reader.load_chromosome(&mut verify_reader, i)?;

            if has_bands {
                let expected = bands_map
                    .get(&chr.ensembl_name)
                    .map(|b| b.len())
                    .unwrap_or(0);
                if chr_data.bands.len() != expected {
                    bail!(
                        "verification failed: chromosome {} band count mismatch ({} vs {})",
                        chr.ucsc_name,
                        chr_data.bands.len(),
                        expected
                    );
                }
                total_bands += chr_data.bands.len();
            }

            if has_mirna {
                let expected = mirna_map
                    .get(&chr.ensembl_name)
                    .map(|m| m.len())
                    .unwrap_or(0);
                if chr_data.mirna.len() != expected {
                    bail!(
                        "verification failed: chromosome {} miRNA count mismatch ({} vs {})",
                        chr.ucsc_name,
                        chr_data.mirna.len(),
                        expected
                    );
                }
                total_mirna += chr_data.mirna.len();
            }
        }

        let mut parts = vec![format!("{} chromosomes", ref_reader.chromosomes.len())];
        if has_bands {
            parts.push(format!("{total_bands} bands"));
        }
        if has_mirna {
            parts.push(format!("{total_mirna} miRNA regions"));
        }
        parts.push(format!("assembly={}", ref_reader.assembly));
        parts.push(format!("patch={}", ref_reader.patch_level));
        success(&parts.join(", "));
    } else {
        success(&format!(
            "{} chromosomes, assembly={}, patch={}",
            ref_reader.chromosomes.len(),
            ref_reader.assembly,
            ref_reader.patch_level
        ));
    }

    Ok(())
}
