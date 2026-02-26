#![allow(unused_assignments, dead_code)]

use std::io::Write;
use std::time::Instant;

use anyhow::{Result, bail};
use clap::Parser;
use colored::Colorize;
use comfy_table::{Attribute, Cell, CellAlignment, ContentArrangement, Table, presets};

use clarus::cli;
use clarus::json::types::{JsonSample, JsonVariant, Position};
use clarus::perf;
use clarus::variant::types::VariantType;

// ─── Constants ───────────────────────────────────────────

const DEFAULT_ITERATIONS: usize = 100_000;
const DEFAULT_BATCH_SIZE: usize = 10_000;
const DEFAULT_RUNS: usize = 5;

// ─── CLI ─────────────────────────────────────────────────

#[derive(Parser)]
#[command(
    name = "json_serialization",
    about = "Benchmark JSON serialization approaches for Position structs"
)]
struct Cli {
    /// Number of iterations for the single-position benchmark
    #[arg(short = 'n', long = "iterations")]
    iterations: Option<usize>,

    /// Number of positions in the batch benchmark
    #[arg(short = 'b', long = "batch-size")]
    batch_size: Option<usize>,

    /// Number of runs (median is reported)
    #[arg(short = 'r', long = "runs")]
    runs: Option<usize>,

    /// Hidden flag: cargo bench passes --bench automatically
    #[arg(long = "bench", hide = true, action = clap::ArgAction::Count)]
    _bench: u8,
}

// ─── Result types ────────────────────────────────────────

struct RunResult {
    elapsed: std::time::Duration,
    total_bytes: u64,
    position_count: u64,
}

struct BenchmarkResult {
    name: String,
    scenario: String,
    median_secs: f64,
    mbps: f64,
    ns_per_position: f64,
    position_count: u64,
}

// ─── Data generators ─────────────────────────────────────

fn make_snv_position(idx: usize) -> Position {
    Position {
        chromosome: "chr1".to_string(),
        position: 100_000 + idx as i64,
        id: None,
        repeat_unit: None,
        ref_repeat_count: None,
        sv_end: None,
        ref_allele: "A".to_string(),
        alt_alleles: vec!["G".to_string()],
        quality: Some(250.0),
        filters: vec!["PASS".to_string()],
        ci_pos: None,
        ci_end: None,
        sv_length: None,
        breakend_event_id: None,
        strand_bias: None,
        fisher_strand_bias: Some(2.435),
        mapping_quality: Some(60.0),
        impute_score: None,
        allele_frequency: None,
        ref_panel_allele_frequency: None,
        cytogenetic_band: Some("1p36.33".to_string()),
        samples: vec![
            JsonSample {
                is_empty: None,
                genotype: Some("0/1".to_string()),
                variant_frequencies: Some(vec![0.45]),
                total_depth: Some(42),
                genotype_quality: Some(99.0),
                copy_number: None,
                minor_haplotype_copy_number: None,
                repeat_unit_counts: None,
                allele_depths: Some(vec![23, 19]),
                failed_filter: false,
                split_read_counts: None,
                paired_end_read_counts: None,
                is_de_novo: false,
                loss_of_heterozygosity: false,
                somatic_quality: None,
                bin_count: None,
                is_imputed_genotype: None,
                genotype_dosage: None,
                genotype_posteriors: None,
            },
            JsonSample {
                is_empty: None,
                genotype: Some("0/0".to_string()),
                variant_frequencies: Some(vec![0.0]),
                total_depth: Some(38),
                genotype_quality: Some(99.0),
                copy_number: None,
                minor_haplotype_copy_number: None,
                repeat_unit_counts: None,
                allele_depths: Some(vec![38, 0]),
                failed_filter: false,
                split_read_counts: None,
                paired_end_read_counts: None,
                is_de_novo: false,
                loss_of_heterozygosity: false,
                somatic_quality: None,
                bin_count: None,
                is_imputed_genotype: None,
                genotype_dosage: None,
                genotype_posteriors: None,
            },
            JsonSample {
                is_empty: None,
                genotype: Some("0/1".to_string()),
                variant_frequencies: Some(vec![0.52]),
                total_depth: Some(31),
                genotype_quality: Some(89.0),
                copy_number: None,
                minor_haplotype_copy_number: None,
                repeat_unit_counts: None,
                allele_depths: Some(vec![15, 16]),
                failed_filter: false,
                split_read_counts: None,
                paired_end_read_counts: None,
                is_de_novo: false,
                loss_of_heterozygosity: false,
                somatic_quality: None,
                bin_count: None,
                is_imputed_genotype: None,
                genotype_dosage: None,
                genotype_posteriors: None,
            },
        ],
        variants: vec![JsonVariant {
            vid: format!("1-{}-A-G", 100_000 + idx),
            chromosome: "chr1".to_string(),
            begin: 100_000 + idx as i64,
            end: 100_001 + idx as i64,
            is_structural_variant: false,
            ref_allele: "A".to_string(),
            alt_allele: "G".to_string(),
            variant_type: VariantType::Snv,
            hgvsg: Some(format!("NC_000001.11:g.{}A>G", 100_000 + idx)),
        }],
    }
}

fn make_indel_position(idx: usize) -> Position {
    Position {
        chromosome: "chr2".to_string(),
        position: 200_000 + idx as i64,
        id: Some(format!("rs{}", 10_000 + idx)),
        repeat_unit: None,
        ref_repeat_count: None,
        sv_end: None,
        ref_allele: "ATCG".to_string(),
        alt_alleles: vec!["A".to_string()],
        quality: Some(180.5),
        filters: vec!["PASS".to_string()],
        ci_pos: None,
        ci_end: None,
        sv_length: None,
        breakend_event_id: None,
        strand_bias: Some(1.23),
        fisher_strand_bias: Some(3.12),
        mapping_quality: Some(55.0),
        impute_score: None,
        allele_frequency: None,
        ref_panel_allele_frequency: None,
        cytogenetic_band: Some("2p25.3".to_string()),
        samples: vec![
            JsonSample {
                is_empty: None,
                genotype: Some("0/1".to_string()),
                variant_frequencies: Some(vec![0.33]),
                total_depth: Some(55),
                genotype_quality: Some(95.0),
                copy_number: None,
                minor_haplotype_copy_number: None,
                repeat_unit_counts: None,
                allele_depths: Some(vec![37, 18]),
                failed_filter: false,
                split_read_counts: None,
                paired_end_read_counts: None,
                is_de_novo: false,
                loss_of_heterozygosity: false,
                somatic_quality: None,
                bin_count: None,
                is_imputed_genotype: None,
                genotype_dosage: None,
                genotype_posteriors: None,
            },
            JsonSample {
                is_empty: None,
                genotype: Some("0/0".to_string()),
                variant_frequencies: Some(vec![0.0]),
                total_depth: Some(40),
                genotype_quality: Some(99.0),
                copy_number: None,
                minor_haplotype_copy_number: None,
                repeat_unit_counts: None,
                allele_depths: Some(vec![40, 0]),
                failed_filter: false,
                split_read_counts: None,
                paired_end_read_counts: None,
                is_de_novo: false,
                loss_of_heterozygosity: false,
                somatic_quality: None,
                bin_count: None,
                is_imputed_genotype: None,
                genotype_dosage: None,
                genotype_posteriors: None,
            },
            JsonSample {
                is_empty: None,
                genotype: Some("0/1".to_string()),
                variant_frequencies: Some(vec![0.48]),
                total_depth: Some(29),
                genotype_quality: Some(78.0),
                copy_number: None,
                minor_haplotype_copy_number: None,
                repeat_unit_counts: None,
                allele_depths: Some(vec![15, 14]),
                failed_filter: false,
                split_read_counts: None,
                paired_end_read_counts: None,
                is_de_novo: false,
                loss_of_heterozygosity: false,
                somatic_quality: None,
                bin_count: None,
                is_imputed_genotype: None,
                genotype_dosage: None,
                genotype_posteriors: None,
            },
        ],
        variants: vec![JsonVariant {
            vid: format!("2-{}-ATCG-A", 200_000 + idx),
            chromosome: "chr2".to_string(),
            begin: 200_000 + idx as i64,
            end: 200_004 + idx as i64,
            is_structural_variant: false,
            ref_allele: "ATCG".to_string(),
            alt_allele: "A".to_string(),
            variant_type: VariantType::Deletion,
            hgvsg: Some(format!(
                "NC_000002.12:g.{}_{}del",
                200_001 + idx,
                200_003 + idx
            )),
        }],
    }
}

fn make_batch(size: usize) -> Vec<Position> {
    (0..size)
        .map(|i| {
            if i % 100 < 85 {
                make_snv_position(i)
            } else {
                make_indel_position(i)
            }
        })
        .collect()
}

// ─── Manual JSON writer ──────────────────────────────────

fn write_json_string(buf: &mut Vec<u8>, s: &str) {
    buf.push(b'"');
    for byte in s.bytes() {
        match byte {
            b'"' => buf.extend_from_slice(b"\\\""),
            b'\\' => buf.extend_from_slice(b"\\\\"),
            b'\n' => buf.extend_from_slice(b"\\n"),
            b'\r' => buf.extend_from_slice(b"\\r"),
            b'\t' => buf.extend_from_slice(b"\\t"),
            0x00..=0x1f => {
                // Control characters: \u00XX
                let _ = write!(buf, "\\u{byte:04x}");
            }
            _ => buf.push(byte),
        }
    }
    buf.push(b'"');
}

fn write_i64(buf: &mut Vec<u8>, v: i64) {
    let mut tmp = itoa::Buffer::new();
    buf.extend_from_slice(tmp.format(v).as_bytes());
}

fn write_i32(buf: &mut Vec<u8>, v: i32) {
    let mut tmp = itoa::Buffer::new();
    buf.extend_from_slice(tmp.format(v).as_bytes());
}

fn write_f64(buf: &mut Vec<u8>, v: f64) {
    let mut tmp = ryu::Buffer::new();
    buf.extend_from_slice(tmp.format(v).as_bytes());
}

fn write_bool(buf: &mut Vec<u8>, v: bool) {
    if v {
        buf.extend_from_slice(b"true");
    } else {
        buf.extend_from_slice(b"false");
    }
}

fn write_string_array(buf: &mut Vec<u8>, arr: &[String]) {
    buf.push(b'[');
    for (i, s) in arr.iter().enumerate() {
        if i > 0 {
            buf.push(b',');
        }
        write_json_string(buf, s);
    }
    buf.push(b']');
}

fn write_f64_array(buf: &mut Vec<u8>, arr: &[f64]) {
    buf.push(b'[');
    for (i, &v) in arr.iter().enumerate() {
        if i > 0 {
            buf.push(b',');
        }
        write_f64(buf, v);
    }
    buf.push(b']');
}

fn write_i64_array(buf: &mut Vec<u8>, arr: &[i64]) {
    buf.push(b'[');
    for (i, &v) in arr.iter().enumerate() {
        if i > 0 {
            buf.push(b',');
        }
        write_i64(buf, v);
    }
    buf.push(b']');
}

fn write_i32_array(buf: &mut Vec<u8>, arr: &[i32]) {
    buf.push(b'[');
    for (i, &v) in arr.iter().enumerate() {
        if i > 0 {
            buf.push(b',');
        }
        write_i32(buf, v);
    }
    buf.push(b']');
}

fn write_variant_type(buf: &mut Vec<u8>, vt: VariantType) {
    let s = match vt {
        VariantType::Snv => "SNV",
        VariantType::Insertion => "insertion",
        VariantType::Deletion => "deletion",
        VariantType::Delins => "delins",
        VariantType::Mnv => "MNV",
        VariantType::Duplication => "duplication",
        VariantType::TandemDuplication => "tandem_duplication",
        VariantType::Inversion => "inversion",
        VariantType::Translocation => "translocation",
        VariantType::MobileElementInsertion => "mobile_element_insertion",
        VariantType::CopyNumberVariation => "copy_number_variation",
        VariantType::CopyNumberLoss => "copy_number_loss",
        VariantType::CopyNumberGain => "copy_number_gain",
        VariantType::ShortTandemRepeatVariation => "short_tandem_repeat_variation",
        VariantType::RunOfHomozygosity => "run_of_homozygosity",
        VariantType::Unknown => "unknown",
    };
    write_json_string(buf, s);
}

fn write_sample_manual(buf: &mut Vec<u8>, sample: &JsonSample) {
    buf.push(b'{');
    let mut first = true;

    macro_rules! field {
        ($name:expr, $writer:expr) => {
            if !first {
                buf.push(b',');
            }
            first = false;
            write_json_string(buf, $name);
            buf.push(b':');
            $writer;
        };
    }

    if let Some(v) = sample.is_empty {
        field!("isEmpty", write_bool(buf, v));
    }
    if let Some(ref v) = sample.genotype {
        field!("genotype", write_json_string(buf, v));
    }
    if let Some(ref v) = sample.variant_frequencies {
        field!("variantFrequencies", write_f64_array(buf, v));
    }
    if let Some(v) = sample.total_depth {
        field!("totalDepth", write_i32(buf, v));
    }
    if let Some(v) = sample.genotype_quality {
        field!("genotypeQuality", write_f64(buf, v));
    }
    if let Some(v) = sample.copy_number {
        field!("copyNumber", write_i32(buf, v));
    }
    if let Some(v) = sample.minor_haplotype_copy_number {
        field!("minorHaplotypeCopyNumber", write_i32(buf, v));
    }
    if let Some(ref v) = sample.repeat_unit_counts {
        field!("repeatUnitCounts", write_i32_array(buf, v));
    }
    if let Some(ref v) = sample.allele_depths {
        field!("alleleDepths", write_i32_array(buf, v));
    }
    if sample.failed_filter {
        field!("failedFilter", write_bool(buf, true));
    }
    if let Some(ref v) = sample.split_read_counts {
        field!("splitReadCounts", write_i32_array(buf, v));
    }
    if let Some(ref v) = sample.paired_end_read_counts {
        field!("pairedEndReadCounts", write_i32_array(buf, v));
    }
    if sample.is_de_novo {
        field!("isDeNovo", write_bool(buf, true));
    }
    if sample.loss_of_heterozygosity {
        field!("lossOfHeterozygosity", write_bool(buf, true));
    }
    if let Some(v) = sample.somatic_quality {
        field!("somaticQuality", write_f64(buf, v));
    }
    if let Some(v) = sample.bin_count {
        field!("binCount", write_i32(buf, v));
    }
    if let Some(v) = sample.is_imputed_genotype {
        field!("isImputedGenotype", write_bool(buf, v));
    }
    if let Some(v) = sample.genotype_dosage {
        field!("genotypeDosage", write_f64(buf, v));
    }
    if let Some(ref v) = sample.genotype_posteriors {
        field!("genotypePosteriors", write_f64_array(buf, v));
    }

    buf.push(b'}');
}

fn write_variant_manual(buf: &mut Vec<u8>, variant: &JsonVariant) {
    buf.push(b'{');

    write_json_string(buf, "vid");
    buf.push(b':');
    write_json_string(buf, &variant.vid);

    buf.push(b',');
    write_json_string(buf, "chromosome");
    buf.push(b':');
    write_json_string(buf, &variant.chromosome);

    buf.push(b',');
    write_json_string(buf, "begin");
    buf.push(b':');
    write_i64(buf, variant.begin);

    buf.push(b',');
    write_json_string(buf, "end");
    buf.push(b':');
    write_i64(buf, variant.end);

    if variant.is_structural_variant {
        buf.push(b',');
        write_json_string(buf, "isStructuralVariant");
        buf.push(b':');
        write_bool(buf, true);
    }

    buf.push(b',');
    write_json_string(buf, "refAllele");
    buf.push(b':');
    write_json_string(buf, &variant.ref_allele);

    buf.push(b',');
    write_json_string(buf, "altAllele");
    buf.push(b':');
    write_json_string(buf, &variant.alt_allele);

    buf.push(b',');
    write_json_string(buf, "variantType");
    buf.push(b':');
    write_variant_type(buf, variant.variant_type);

    if let Some(ref v) = variant.hgvsg {
        buf.push(b',');
        write_json_string(buf, "hgvsg");
        buf.push(b':');
        write_json_string(buf, v);
    }

    buf.push(b'}');
}

fn write_position_manual(buf: &mut Vec<u8>, pos: &Position) {
    buf.push(b'{');

    // chromosome
    write_json_string(buf, "chromosome");
    buf.push(b':');
    write_json_string(buf, &pos.chromosome);

    // position
    buf.push(b',');
    write_json_string(buf, "position");
    buf.push(b':');
    write_i64(buf, pos.position);

    // id
    if let Some(ref v) = pos.id {
        buf.push(b',');
        write_json_string(buf, "id");
        buf.push(b':');
        write_json_string(buf, v);
    }

    // repeatUnit
    if let Some(ref v) = pos.repeat_unit {
        buf.push(b',');
        write_json_string(buf, "repeatUnit");
        buf.push(b':');
        write_json_string(buf, v);
    }

    // refRepeatCount
    if let Some(v) = pos.ref_repeat_count {
        buf.push(b',');
        write_json_string(buf, "refRepeatCount");
        buf.push(b':');
        write_i32(buf, v);
    }

    // svEnd
    if let Some(v) = pos.sv_end {
        buf.push(b',');
        write_json_string(buf, "svEnd");
        buf.push(b':');
        write_i64(buf, v);
    }

    // refAllele
    buf.push(b',');
    write_json_string(buf, "refAllele");
    buf.push(b':');
    write_json_string(buf, &pos.ref_allele);

    // altAlleles
    buf.push(b',');
    write_json_string(buf, "altAlleles");
    buf.push(b':');
    write_string_array(buf, &pos.alt_alleles);

    // quality
    if let Some(v) = pos.quality {
        buf.push(b',');
        write_json_string(buf, "quality");
        buf.push(b':');
        write_f64(buf, v);
    }

    // filters
    if !pos.filters.is_empty() {
        buf.push(b',');
        write_json_string(buf, "filters");
        buf.push(b':');
        write_string_array(buf, &pos.filters);
    }

    // ciPos
    if let Some(ref v) = pos.ci_pos {
        buf.push(b',');
        write_json_string(buf, "ciPos");
        buf.push(b':');
        write_i64_array(buf, v);
    }

    // ciEnd
    if let Some(ref v) = pos.ci_end {
        buf.push(b',');
        write_json_string(buf, "ciEnd");
        buf.push(b':');
        write_i64_array(buf, v);
    }

    // svLength
    if let Some(v) = pos.sv_length {
        buf.push(b',');
        write_json_string(buf, "svLength");
        buf.push(b':');
        write_i64(buf, v);
    }

    // breakendEventId
    if let Some(ref v) = pos.breakend_event_id {
        buf.push(b',');
        write_json_string(buf, "breakendEventId");
        buf.push(b':');
        write_json_string(buf, v);
    }

    // strandBias
    if let Some(v) = pos.strand_bias {
        buf.push(b',');
        write_json_string(buf, "strandBias");
        buf.push(b':');
        write_f64(buf, v);
    }

    // fisherStrandBias
    if let Some(v) = pos.fisher_strand_bias {
        buf.push(b',');
        write_json_string(buf, "fisherStrandBias");
        buf.push(b':');
        write_f64(buf, v);
    }

    // mappingQuality
    if let Some(v) = pos.mapping_quality {
        buf.push(b',');
        write_json_string(buf, "mappingQuality");
        buf.push(b':');
        write_f64(buf, v);
    }

    // imputeScore
    if let Some(v) = pos.impute_score {
        buf.push(b',');
        write_json_string(buf, "imputeScore");
        buf.push(b':');
        write_f64(buf, v);
    }

    // alleleFrequency
    if let Some(ref v) = pos.allele_frequency {
        buf.push(b',');
        write_json_string(buf, "alleleFrequency");
        buf.push(b':');
        write_f64_array(buf, v);
    }

    // refPanelAlleleFrequency
    if let Some(ref v) = pos.ref_panel_allele_frequency {
        buf.push(b',');
        write_json_string(buf, "refPanelAlleleFrequency");
        buf.push(b':');
        write_f64_array(buf, v);
    }

    // cytogeneticBand
    if let Some(ref v) = pos.cytogenetic_band {
        buf.push(b',');
        write_json_string(buf, "cytogeneticBand");
        buf.push(b':');
        write_json_string(buf, v);
    }

    // samples
    buf.push(b',');
    write_json_string(buf, "samples");
    buf.push(b':');
    buf.push(b'[');
    for (i, sample) in pos.samples.iter().enumerate() {
        if i > 0 {
            buf.push(b',');
        }
        write_sample_manual(buf, sample);
    }
    buf.push(b']');

    // variants
    buf.push(b',');
    write_json_string(buf, "variants");
    buf.push(b':');
    buf.push(b'[');
    for (i, variant) in pos.variants.iter().enumerate() {
        if i > 0 {
            buf.push(b',');
        }
        write_variant_manual(buf, variant);
    }
    buf.push(b']');

    buf.push(b'}');
}

// ─── Benchmark functions ─────────────────────────────────

fn bench_serde_json_single(pos: &Position, iterations: usize) -> RunResult {
    let mut buf = Vec::with_capacity(1024);
    let start = Instant::now();
    let mut total_bytes: u64 = 0;

    for _ in 0..iterations {
        buf.clear();
        serde_json::to_writer(&mut buf, pos).unwrap();
        total_bytes += buf.len() as u64;
    }

    RunResult {
        elapsed: start.elapsed(),
        total_bytes,
        position_count: iterations as u64,
    }
}

fn bench_sonic_rs_single(pos: &Position, iterations: usize) -> RunResult {
    let mut buf = Vec::with_capacity(1024);
    let start = Instant::now();
    let mut total_bytes: u64 = 0;

    for _ in 0..iterations {
        buf.clear();
        sonic_rs::to_writer(&mut buf, pos).unwrap();
        total_bytes += buf.len() as u64;
    }

    RunResult {
        elapsed: start.elapsed(),
        total_bytes,
        position_count: iterations as u64,
    }
}

fn bench_simd_json_single(pos: &Position, iterations: usize) -> RunResult {
    let mut buf = Vec::with_capacity(1024);
    let start = Instant::now();
    let mut total_bytes: u64 = 0;

    for _ in 0..iterations {
        buf.clear();
        simd_json::to_writer(&mut buf, pos).unwrap();
        total_bytes += buf.len() as u64;
    }

    RunResult {
        elapsed: start.elapsed(),
        total_bytes,
        position_count: iterations as u64,
    }
}

fn bench_manual_single(pos: &Position, iterations: usize) -> RunResult {
    let mut buf = Vec::with_capacity(1024);
    let start = Instant::now();
    let mut total_bytes: u64 = 0;

    for _ in 0..iterations {
        buf.clear();
        write_position_manual(&mut buf, pos);
        total_bytes += buf.len() as u64;
    }

    RunResult {
        elapsed: start.elapsed(),
        total_bytes,
        position_count: iterations as u64,
    }
}

fn bench_serde_json_batch(batch: &[Position]) -> RunResult {
    let mut buf = Vec::with_capacity(1024);
    let start = Instant::now();
    let mut total_bytes: u64 = 0;

    for pos in batch {
        buf.clear();
        serde_json::to_writer(&mut buf, pos).unwrap();
        total_bytes += buf.len() as u64;
    }

    RunResult {
        elapsed: start.elapsed(),
        total_bytes,
        position_count: batch.len() as u64,
    }
}

fn bench_sonic_rs_batch(batch: &[Position]) -> RunResult {
    let mut buf = Vec::with_capacity(1024);
    let start = Instant::now();
    let mut total_bytes: u64 = 0;

    for pos in batch {
        buf.clear();
        sonic_rs::to_writer(&mut buf, pos).unwrap();
        total_bytes += buf.len() as u64;
    }

    RunResult {
        elapsed: start.elapsed(),
        total_bytes,
        position_count: batch.len() as u64,
    }
}

fn bench_simd_json_batch(batch: &[Position]) -> RunResult {
    let mut buf = Vec::with_capacity(1024);
    let start = Instant::now();
    let mut total_bytes: u64 = 0;

    for pos in batch {
        buf.clear();
        simd_json::to_writer(&mut buf, pos).unwrap();
        total_bytes += buf.len() as u64;
    }

    RunResult {
        elapsed: start.elapsed(),
        total_bytes,
        position_count: batch.len() as u64,
    }
}

fn bench_manual_batch(batch: &[Position]) -> RunResult {
    let mut buf = Vec::with_capacity(1024);
    let start = Instant::now();
    let mut total_bytes: u64 = 0;

    for pos in batch {
        buf.clear();
        write_position_manual(&mut buf, pos);
        total_bytes += buf.len() as u64;
    }

    RunResult {
        elapsed: start.elapsed(),
        total_bytes,
        position_count: batch.len() as u64,
    }
}

// ─── Runner ──────────────────────────────────────────────

fn run_benchmark<F>(name: &str, scenario: &str, runs: usize, f: F) -> BenchmarkResult
where
    F: Fn() -> RunResult,
{
    eprint!("  {:<40} ", format!("{name} ({scenario})"));

    let mut results: Vec<RunResult> = Vec::with_capacity(runs);
    for _ in 0..runs {
        results.push(f());
    }

    // Sort by elapsed to find median
    results.sort_by(|a, b| a.elapsed.cmp(&b.elapsed));
    let median = &results[runs / 2];
    let median_secs = median.elapsed.as_secs_f64();
    let mbps = (median.total_bytes as f64 / (1024.0 * 1024.0)) / median_secs;
    let ns_per_position = if median.position_count > 0 {
        median_secs * 1e9 / median.position_count as f64
    } else {
        0.0
    };

    eprintln!(
        "{:.3}s  ({:.1} MB/s, {:.0} ns/pos)",
        median_secs, mbps, ns_per_position
    );

    BenchmarkResult {
        name: name.to_string(),
        scenario: scenario.to_string(),
        median_secs,
        mbps,
        ns_per_position,
        position_count: median.position_count,
    }
}

// ─── Results table ───────────────────────────────────────

fn print_results_table(title: &str, results: &[BenchmarkResult]) {
    if results.is_empty() {
        return;
    }

    eprintln!("  {}", title.bold());
    eprintln!();

    let baseline_secs = results[0].median_secs;

    let mut table = Table::new();
    table
        .load_preset(presets::UTF8_FULL_CONDENSED)
        .set_content_arrangement(ContentArrangement::Dynamic);

    table.set_header(vec![
        Cell::new("Approach").add_attribute(Attribute::Bold),
        Cell::new("Median (s)")
            .add_attribute(Attribute::Bold)
            .set_alignment(CellAlignment::Right),
        Cell::new("MB/s")
            .add_attribute(Attribute::Bold)
            .set_alignment(CellAlignment::Right),
        Cell::new("ns/position")
            .add_attribute(Attribute::Bold)
            .set_alignment(CellAlignment::Right),
        Cell::new("Speedup")
            .add_attribute(Attribute::Bold)
            .set_alignment(CellAlignment::Right),
    ]);

    for r in results {
        let speedup = baseline_secs / r.median_secs;
        table.add_row(vec![
            Cell::new(&r.name),
            Cell::new(format!("{:.3}", r.median_secs)).set_alignment(CellAlignment::Right),
            Cell::new(format!("{:.1}", r.mbps)).set_alignment(CellAlignment::Right),
            Cell::new(format!("{:.0}", r.ns_per_position)).set_alignment(CellAlignment::Right),
            Cell::new(format!("{speedup:.2}x")).set_alignment(CellAlignment::Right),
        ]);
    }

    for line in table.to_string().lines() {
        eprintln!("  {line}");
    }
    eprintln!();
}

// ─── Markdown output ─────────────────────────────────────

fn write_markdown(
    single_results: &[BenchmarkResult],
    batch_results: &[BenchmarkResult],
    iterations: usize,
    batch_size: usize,
    runs: usize,
    snv_bytes: usize,
    indel_bytes: usize,
) -> Result<()> {
    let path = std::path::Path::new("docs/clarus");
    std::fs::create_dir_all(path)?;

    let mut out = String::new();

    out.push_str("# JSON Serialization Benchmark Results\n\n");
    out.push_str(&format!(
        "**Date:** {}\n\n",
        chrono::Local::now().format("%Y-%m-%d %H:%M")
    ));

    out.push_str("## Test Setup\n\n");
    out.push_str("| Parameter | Value |\n");
    out.push_str("|-----------|-------|\n");
    out.push_str(&format!(
        "| Single iterations | {} |\n",
        cli::num(iterations)
    ));
    out.push_str(&format!("| Batch size | {} |\n", cli::num(batch_size)));
    out.push_str(&format!("| Runs (median) | {} |\n", runs));
    out.push_str(&format!("| SNV position size | ~{} bytes |\n", snv_bytes));
    out.push_str(&format!(
        "| Indel position size | ~{} bytes |\n",
        indel_bytes
    ));
    out.push_str(&format!("| Batch mix | 85% SNV / 15% indel |\n"));
    out.push_str("\n### Crates\n\n");
    out.push_str("| Crate | Key Feature |\n");
    out.push_str("|-------|-------------|\n");
    out.push_str("| **serde_json** | Baseline — the standard |\n");
    out.push_str("| **sonic-rs** | SIMD string escaping (AVX2/NEON) |\n");
    out.push_str("| **simd-json** | SIMD-accelerated (mainly parsing, also serialization) |\n");
    out.push_str(
        "| **manual writer** | Direct `Vec<u8>` writes with `itoa`/`ryu`, zero serde overhead |\n",
    );
    out.push_str("\n");

    // Single results
    out.push_str("## Single Position (hot cache)\n\n");
    if !single_results.is_empty() {
        let baseline = single_results[0].median_secs;
        out.push_str("| Approach | Median (s) | MB/s | ns/position | Speedup |\n");
        out.push_str("|----------|----------:|-----:|------------:|--------:|\n");
        for r in single_results {
            let speedup = baseline / r.median_secs;
            out.push_str(&format!(
                "| {} | {:.3} | {:.1} | {:.0} | {:.2}x |\n",
                r.name, r.median_secs, r.mbps, r.ns_per_position, speedup
            ));
        }
    }
    out.push_str("\n");

    // Batch results
    out.push_str("## Batch (mixed positions)\n\n");
    if !batch_results.is_empty() {
        let baseline = batch_results[0].median_secs;
        out.push_str("| Approach | Median (s) | MB/s | ns/position | Speedup |\n");
        out.push_str("|----------|----------:|-----:|------------:|--------:|\n");
        for r in batch_results {
            let speedup = baseline / r.median_secs;
            out.push_str(&format!(
                "| {} | {:.3} | {:.1} | {:.0} | {:.2}x |\n",
                r.name, r.median_secs, r.mbps, r.ns_per_position, speedup
            ));
        }
    }
    out.push_str("\n");

    out.push_str("## Analysis\n\n");
    out.push_str("*(To be filled in after reviewing results)*\n");

    std::fs::write("docs/clarus/fast_json_serialization.md", &out)?;
    Ok(())
}

// ─── Main ────────────────────────────────────────────────

fn main() -> Result<()> {
    let start = Instant::now();
    let cli_args = Cli::parse();

    let iterations = cli_args.iterations.unwrap_or(DEFAULT_ITERATIONS);
    let batch_size = cli_args.batch_size.unwrap_or(DEFAULT_BATCH_SIZE);
    let runs = cli_args.runs.unwrap_or(DEFAULT_RUNS);

    cli::banner("JSON Serialization Benchmark");

    // ── Configuration ────────────────────────────────────
    cli::section("Configuration");
    cli::kv("Iterations (single)", &cli::num(iterations));
    cli::kv("Batch size", &cli::num(batch_size));
    cli::kv("Runs (median)", &cli::num(runs));
    eprintln!();

    // ── Data setup ───────────────────────────────────────
    cli::section("Data Setup");

    let snv = make_snv_position(0);
    let snv_bytes = serde_json::to_vec(&snv)?.len();
    cli::kv("SNV position size", &format!("{} bytes", snv_bytes));

    let indel = make_indel_position(0);
    let indel_bytes = serde_json::to_vec(&indel)?.len();
    cli::kv("Indel position size", &format!("{} bytes", indel_bytes));

    let batch = make_batch(batch_size);
    let batch_total: usize = batch
        .iter()
        .map(|p| serde_json::to_vec(p).unwrap().len())
        .sum();
    cli::kv("Batch total JSON", &perf::format_bytes(batch_total as u64));
    eprintln!();

    // ── Validation ───────────────────────────────────────
    cli::section("Validation");

    // Compare manual writer output to serde_json for both SNV and indel
    validate_manual_writer(&snv, "SNV")?;
    validate_manual_writer(&indel, "Indel")?;
    eprintln!();

    // ── Warmup ───────────────────────────────────────────
    cli::section("Warmup");
    eprint!("  {}", "Warming up serializers...".dimmed());
    // Run each once to warm up
    let mut warm_buf = Vec::with_capacity(1024);
    serde_json::to_writer(&mut warm_buf, &snv).unwrap();
    warm_buf.clear();
    sonic_rs::to_writer(&mut warm_buf, &snv).unwrap();
    warm_buf.clear();
    simd_json::to_writer(&mut warm_buf, &snv).unwrap();
    warm_buf.clear();
    write_position_manual(&mut warm_buf, &snv);
    eprintln!(" done");
    eprintln!();

    // ── Single Position Benchmark ────────────────────────
    cli::section("Single Position Benchmark");

    let single_results = vec![
        run_benchmark("serde_json", "single", runs, || {
            bench_serde_json_single(&snv, iterations)
        }),
        run_benchmark("sonic-rs", "single", runs, || {
            bench_sonic_rs_single(&snv, iterations)
        }),
        run_benchmark("simd-json", "single", runs, || {
            bench_simd_json_single(&snv, iterations)
        }),
        run_benchmark("manual writer", "single", runs, || {
            bench_manual_single(&snv, iterations)
        }),
    ];

    eprintln!();

    // ── Batch Benchmark ──────────────────────────────────
    cli::section("Batch Benchmark");

    let batch_results = vec![
        run_benchmark("serde_json", "batch", runs, || {
            bench_serde_json_batch(&batch)
        }),
        run_benchmark("sonic-rs", "batch", runs, || bench_sonic_rs_batch(&batch)),
        run_benchmark("simd-json", "batch", runs, || bench_simd_json_batch(&batch)),
        run_benchmark("manual writer", "batch", runs, || {
            bench_manual_batch(&batch)
        }),
    ];

    eprintln!();

    // ── Results ──────────────────────────────────────────
    cli::section("Results");
    print_results_table("Single Position (hot cache)", &single_results);
    print_results_table("Batch (mixed positions)", &batch_results);

    // ── Write markdown ───────────────────────────────────
    write_markdown(
        &single_results,
        &batch_results,
        iterations,
        batch_size,
        runs,
        snv_bytes,
        indel_bytes,
    )?;
    cli::success("Results written to docs/clarus/fast_json_serialization.md");

    // ── Summary ──────────────────────────────────────────
    cli::print_summary(start);
    Ok(())
}

// ─── Validation helper ───────────────────────────────────

fn validate_manual_writer(pos: &Position, label: &str) -> Result<()> {
    let serde_output = serde_json::to_vec(pos)?;
    let mut manual_output = Vec::with_capacity(1024);
    write_position_manual(&mut manual_output, pos);

    // Parse both into serde_json::Value for semantic comparison
    // (ryu and serde_json may format floats slightly differently)
    let serde_value: serde_json::Value = serde_json::from_slice(&serde_output)?;
    let manual_value: serde_json::Value = serde_json::from_slice(&manual_output).map_err(|e| {
        // On parse failure, show what the manual writer produced for debugging
        let preview = String::from_utf8_lossy(&manual_output);
        anyhow::anyhow!("Manual writer produced invalid JSON for {label}: {e}\nOutput: {preview}")
    })?;

    if serde_value != manual_value {
        let serde_pretty = serde_json::to_string_pretty(&serde_value)?;
        let manual_pretty = serde_json::to_string_pretty(&manual_value)?;
        bail!(
            "Manual writer mismatch for {label}!\n\n--- serde_json ---\n{serde_pretty}\n\n--- manual ---\n{manual_pretty}"
        );
    }

    cli::success(&format!(
        "{label}: manual writer matches serde_json ({} bytes)",
        serde_output.len()
    ));
    Ok(())
}
