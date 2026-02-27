//! Clarus annotator: VCF → JSON annotation pipeline.

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result, bail};
use clap::Parser;

use clarus::chromosome::Chromosome;
use clarus::cli;
use clarus::genome_assembly::GenomeAssembly;
use clarus::hgvsg;
use clarus::json::types::{
    AnnotationStats, DataSource, InputInfo, JsonSample, JsonVariant, Metadata, Position,
};
use clarus::json::writer::OutputWriter;
use clarus::perf;
use clarus::reference::reader::{ChromosomeData, ReferenceReader};
use clarus::variant::categorize;
use clarus::variant::normalize;
use clarus::variant::types::{VariantCategory, VariantType};
use clarus::variant::vid;
use clarus::vcf::assembly_inference::normalize_contig_name;
use clarus::vcf::streaming_reader::VcfStream;

#[derive(Parser)]
#[command(
    name = "clarus",
    about = "High-performance genomic variant annotation engine"
)]
struct Cli {
    /// Input VCF file (bgzf-compressed or plain gzip)
    #[arg(short = 'i', long = "input")]
    input: PathBuf,

    /// Output prefix (produces <prefix>_metadata.json.gz, <prefix>_variants.jsonl.gz, <prefix>_genes.jsonl.gz)
    #[arg(short = 'o', long = "output")]
    output: PathBuf,

    /// Data root directory
    #[arg(short = 'd', long = "data")]
    data: PathBuf,

    /// Genome assembly
    #[arg(long = "ga")]
    genome_assembly: String,

    /// Annotation source (refseq or ensembl)
    #[arg(long = "source")]
    source: String,
}

fn main() -> Result<()> {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let start = Instant::now();
    let cli_args = Cli::parse();

    cli::banner("Annotator");

    // ── Configuration ────────────────────────────────────
    cli::section("Configuration");

    let assembly: GenomeAssembly = cli_args.genome_assembly.parse()?;
    if assembly == GenomeAssembly::Unknown {
        bail!("genome assembly must be GRCh37 or GRCh38");
    }

    let source = cli_args.source.to_lowercase();
    if source != "refseq" && source != "ensembl" {
        bail!("source must be 'refseq' or 'ensembl'");
    }

    cli::kv("Input", &cli_args.input.display().to_string());
    cli::kv("Output prefix", &cli_args.output.display().to_string());
    cli::kv("Data root", &cli_args.data.display().to_string());
    cli::kv("Assembly", &assembly.to_string());
    cli::kv("Source", &source);

    eprintln!();

    // ── File Discovery ───────────────────────────────────
    cli::section("File Discovery");

    let ref_path = discover_reference_file(&cli_args.data, assembly)?;
    cli::kv("Reference", &ref_path.display().to_string());

    let cache_path = discover_cache_file(&cli_args.data, assembly, &source)?;
    cli::kv("Cache", &cache_path.display().to_string());

    eprintln!();

    // ── Load Reference ───────────────────────────────────
    cli::section("Reference");

    let ref_file = std::fs::File::open(&ref_path)
        .with_context(|| format!("cannot open reference: {}", ref_path.display()))?;
    let mut ref_reader_file = BufReader::new(ref_file);
    let ref_reader = ReferenceReader::from_reader(&mut ref_reader_file)?;

    if ref_reader.assembly != assembly {
        bail!(
            "reference assembly mismatch: file is {} but expected {}",
            ref_reader.assembly,
            assembly
        );
    }

    cli::kv("Chromosomes", &cli::num(ref_reader.chromosomes.len()));
    cli::kv("Assembly", &ref_reader.assembly.to_string());
    cli::kv("Patch level", &format!("p{}", ref_reader.patch_level));
    eprintln!();

    // ── Read VCF ─────────────────────────────────────────
    cli::section("VCF Input");

    let vcf_start = Instant::now();
    let mut vcf = VcfStream::open(&cli_args.input, assembly)?;

    cli::kv("File format", &vcf.header.file_format);
    cli::kv("Contigs", &cli::num(vcf.header.contigs.len()));
    cli::kv("Samples", &cli::num(vcf.header.sample_names.len()));
    cli::kv("Inferred assembly", &vcf.inferred_assembly.to_string());
    cli::kv("VCF open time", &perf::format_elapsed(vcf_start.elapsed()));
    eprintln!();

    // ── Annotate & Write JSON ────────────────────────────
    cli::section("Annotation");

    let mut writer = OutputWriter::new(&cli_args.output)?;

    let input_info = InputInfo {
        file_name: cli_args
            .input
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default(),
        vcf_version: vcf.header.file_format.clone(),
        contigs: vcf.header.contigs.len(),
        samples: vcf.header.sample_names.len(),
        inferred_assembly: vcf.inferred_assembly.to_string(),
    };

    let mut current_chrom_index: Option<usize> = None;
    let mut current_chrom_data: Option<ChromosomeData> = None;
    let mut positions_written = 0u64;
    let mut variants_written = 0u64;
    let mut skipped_chroms = 0u64;
    let mut vid_buf = String::with_capacity(64);
    let mut hgvsg_buf = String::with_capacity(128);

    while let Some(record) = vcf.next_record()? {
        let ensembl_name = normalize_contig_name(&record.chromosome);

        let chrom_index = match ref_reader.get_index(&record.chromosome) {
            Some(idx) => idx,
            None => {
                skipped_chroms += 1;
                continue;
            }
        };

        // Load chromosome data if different from current
        if current_chrom_index != Some(chrom_index) {
            // Drop the old chromosome before allocating the new one to avoid
            // holding two ~250 MB sequences simultaneously during zstd decode.
            drop(current_chrom_data.take());
            current_chrom_data =
                Some(ref_reader.load_chromosome(&mut ref_reader_file, chrom_index)?);
            current_chrom_index = Some(chrom_index);
        }

        let chrom_data = current_chrom_data
            .as_ref()
            .expect("chrom_data must be loaded after successful chrom_index lookup");
        let chromosome = &ref_reader.chromosomes[chrom_index];

        let ref_accessor =
            |offset: i64, len: usize| -> Option<&[u8]> { chrom_data.substring(offset, len) };

        // Build variants for each ALT allele (borrows record fields)
        let svtype = record.info.svtype.as_deref();
        let mut json_variants: Vec<JsonVariant> = Vec::with_capacity(record.alt_alleles.len());
        for alt in &record.alt_alleles {
            json_variants.push(annotate_alt(
                ensembl_name,
                chromosome,
                record.position,
                &record.ref_allele,
                alt,
                svtype,
                record.end_position,
                &ref_accessor,
                &mut vid_buf,
                &mut hgvsg_buf,
            ));
        }

        // Convert samples by consuming (moves owned fields, zero clones)
        let json_samples: Vec<JsonSample> =
            record.samples.into_iter().map(JsonSample::from).collect();

        let cyto_band = find_cytogenetic_band(chrom_data, record.position).map(str::to_string);

        // Build Position, moving remaining owned fields (zero clones)
        let position = Position {
            chromosome: ensembl_name.to_string(),
            position: record.position,
            ref_allele: record.ref_allele,
            alt_alleles: record.alt_alleles,
            quality: record.quality,
            filters: record.filters,
            strand_bias: record.info.strand_bias,
            fisher_strand_bias: record.info.fisher_strand_bias,
            mapping_quality: record.info.mapping_quality,
            impute_score: record.info.impute_score,
            allele_frequency: record.info.af,
            ref_panel_allele_frequency: record.info.raf,
            cytogenetic_band: cyto_band,
            samples: json_samples,
            variants: json_variants,
        };

        variants_written += position.variants.len() as u64;
        writer.write_position(&position)?;
        positions_written += 1;
    }

    let metadata = Metadata {
        annotator: "Clarus 0.1.0".to_string(),
        creation_time: chrono::Utc::now().format("%Y-%m-%d %H:%M:%S").to_string(),
        genome_assembly: assembly.to_string(),
        schema_version: 1,
        data_sources: vec![DataSource {
            name: "Reference".to_string(),
            version: format!("{}.p{}", assembly, ref_reader.patch_level),
            description: None,
            release_date: None,
        }],
        samples: vcf.header.sample_names.clone(),
        input: input_info,
        stats: AnnotationStats {
            positions: positions_written,
            variants: variants_written,
            skipped_records: skipped_chroms,
        },
    };
    writer.finish(&metadata)?;

    cli::kv("Positions", &cli::num(positions_written));
    cli::kv("Variants", &cli::num(variants_written));
    if skipped_chroms > 0 {
        cli::warning(&format!(
            "Skipped {} records on unrecognized chromosomes",
            cli::num(skipped_chroms)
        ));
    }
    let prefix_str = cli_args.output.display().to_string();
    cli::success(&format!("Wrote {prefix_str}_metadata.json.gz"));
    cli::success(&format!("Wrote {prefix_str}_variants.jsonl.gz"));
    cli::success(&format!("Wrote {prefix_str}_genes.jsonl.gz"));
    eprintln!();

    // ── Summary ──────────────────────────────────────────
    cli::print_summary(start);
    Ok(())
}

/// Discover the reference file with the highest patch number.
/// Pattern: `{data_root}/reference/{assembly}.p*.dat`
fn discover_reference_file(data_root: &Path, assembly: GenomeAssembly) -> Result<PathBuf> {
    let ref_dir = data_root.join("reference");
    let prefix = format!("{}.", assembly);

    let mut best: Option<(u8, PathBuf)> = None;

    if ref_dir.is_dir() {
        for entry in std::fs::read_dir(&ref_dir)? {
            let entry = entry?;
            let name = entry.file_name();
            let name_str = name.to_string_lossy();

            if name_str.starts_with(&prefix) && name_str.ends_with(".dat") {
                // Extract patch number from "GRCh38.p14.dat"
                let middle = &name_str[prefix.len()..name_str.len() - 4];
                if let Some(patch_str) = middle.strip_prefix('p')
                    && let Ok(patch) = patch_str.parse::<u8>()
                    && (best.is_none() || patch > best.as_ref().unwrap().0)
                {
                    best = Some((patch, entry.path()));
                }
            }
        }
    }

    best.map(|(_, path)| path).ok_or_else(|| {
        anyhow::anyhow!(
            "no reference file found matching {prefix}p*.dat in {}",
            ref_dir.display()
        )
    })
}

/// Discover the cache file with the latest date stamp.
/// Pattern: `{data_root}/cache/{assembly}_{source}_*.cache` (case-insensitive source match).
fn discover_cache_file(
    data_root: &Path,
    assembly: GenomeAssembly,
    source: &str,
) -> Result<PathBuf> {
    let cache_dir = data_root.join("cache");
    let assembly_str = assembly.to_string();
    let source_lower = source.to_lowercase();

    let mut best: Option<(String, PathBuf)> = None;

    if cache_dir.is_dir() {
        for entry in std::fs::read_dir(&cache_dir)? {
            let entry = entry?;
            let name = entry.file_name();
            let name_str = name.to_string_lossy();

            // Match: {assembly}_{source}_*.cache with case-insensitive source
            if !name_str.ends_with(".cache") || name_str.ends_with(".cache.idx") {
                continue;
            }

            let name_lower = name_str.to_lowercase();
            let prefix_lower = format!("{}_{}_", assembly_str.to_lowercase(), source_lower);

            if name_lower.starts_with(&prefix_lower) {
                let date_part = &name_str[prefix_lower.len()..name_str.len() - 6];
                let date_str = date_part.to_string();
                if best.is_none() || *date_str > *best.as_ref().unwrap().0 {
                    best = Some((date_str, entry.path()));
                }
            }
        }
    }

    best.map(|(_, path)| path).ok_or_else(|| {
        anyhow::anyhow!(
            "no cache file found matching {assembly_str}_{source}_*.cache in {}",
            cache_dir.display()
        )
    })
}

/// Annotate a single ALT allele: normalize, categorize, compute VID and HGVSg.
#[allow(clippy::too_many_arguments)]
fn annotate_alt<'a>(
    ensembl_name: &str,
    chromosome: &Chromosome,
    record_position: i64,
    ref_allele: &str,
    alt: &str,
    svtype: Option<&str>,
    end_position: i64,
    ref_accessor: &dyn Fn(i64, usize) -> Option<&'a [u8]>,
    vid_buf: &mut String,
    hgvsg_buf: &mut String,
) -> JsonVariant {
    let trimmed = normalize::trim(record_position, ref_allele, alt);

    let categorized = categorize::categorize(
        alt,
        svtype,
        trimmed.ref_allele.len(),
        trimmed.alt_allele.len(),
    );

    let normalized = if categorized.category == VariantCategory::SmallVariant {
        normalize::normalize(record_position, ref_allele, alt, ref_accessor)
    } else {
        trimmed
    };

    let (begin, end) = match categorized.category {
        VariantCategory::SmallVariant => {
            let b = normalized.position;
            let e = if normalized.ref_allele.is_empty() {
                b
            } else {
                b + normalized.ref_allele.len() as i64 - 1
            };
            (b, e)
        }
        VariantCategory::Sv if categorized.variant_type == VariantType::Translocation => {
            (record_position, end_position)
        }
        _ => (record_position + 1, end_position),
    };

    let ref_base =
        ref_accessor(normalized.position - 1, 1).map(|b| String::from_utf8_lossy(b).into_owned());

    vid::construct_vid_into(
        vid_buf,
        ensembl_name,
        begin,
        end,
        &normalized.ref_allele,
        &normalized.alt_allele,
        categorized.category,
        categorized.variant_type,
        ref_base.as_deref(),
    );
    let variant_id = vid_buf.clone();

    let hgvsg_notation = if categorized.category == VariantCategory::SmallVariant
        && (!normalized.ref_allele.is_empty() || !normalized.alt_allele.is_empty())
    {
        let refseq_acc = &chromosome.refseq_accession;
        if !refseq_acc.is_empty() {
            hgvsg::compute_hgvsg_into(
                hgvsg_buf,
                refseq_acc,
                normalized.position,
                &normalized.ref_allele,
                &normalized.alt_allele,
                ref_accessor,
            );
            Some(hgvsg_buf.clone())
        } else {
            None
        }
    } else {
        None
    };

    // Convert Cow alleles to owned display strings
    let display_ref = if normalized.ref_allele.is_empty() {
        "-".to_string()
    } else {
        normalized.ref_allele.into_owned()
    };
    let display_alt = if normalized.alt_allele.is_empty() {
        "-".to_string()
    } else {
        normalized.alt_allele.into_owned()
    };

    JsonVariant {
        vid: variant_id,
        chromosome: ensembl_name.to_string(),
        begin,
        end,
        ref_allele: display_ref,
        alt_allele: display_alt,
        variant_type: categorized.variant_type,
        hgvsg: hgvsg_notation,
    }
}

/// Find the cytogenetic band containing the given position (binary search).
fn find_cytogenetic_band(chrom_data: &ChromosomeData, position: i64) -> Option<&str> {
    let pos = position as u32;
    let bands = &chrom_data.bands;

    // Bands are sorted by begin. Find the rightmost band whose begin <= pos.
    let idx = bands.partition_point(|band| band.begin <= pos);
    if idx == 0 {
        return None;
    }
    let band = &bands[idx - 1];
    if pos <= band.end {
        Some(&band.name)
    } else {
        None
    }
}
