use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result, bail};
use clap::Parser;
use colored::Colorize;
use comfy_table::{Attribute, Cell, CellAlignment, ContentArrangement, Table, presets};

use clarus::cli;
use clarus::perf;

/// BGZF magic bytes: gzip magic (0x1f, 0x8b) + deflate method (0x08) + FEXTRA flag (0x04)
const BGZF_MAGIC: [u8; 4] = [0x1f, 0x8b, 0x08, 0x04];

/// Number of benchmark iterations per approach.
const ITERATIONS: usize = 5;

#[derive(Parser)]
#[command(
    name = "bgzf_benchmark",
    about = "Benchmark BGZF decompression approaches for VCF parsing"
)]
struct Cli {
    /// Path to a BGZF-compressed file (e.g., a .vcf.gz)
    #[arg(short = 'i', long = "input")]
    input: PathBuf,

    /// Only run the flate2 benchmark (for testing different feature backends)
    #[arg(long = "flate2-only")]
    flate2_only: bool,

    /// Number of iterations per approach (default: 5)
    #[arg(short = 'n', long = "iterations")]
    iterations: Option<usize>,

    /// Hidden flag: cargo bench passes --bench automatically
    #[arg(long = "bench", hide = true)]
    _bench: bool,
}

/// Result from a single benchmark run.
struct RunResult {
    elapsed: std::time::Duration,
    decompressed_bytes: u64,
    line_count: u64,
}

/// Summary for one benchmark approach.
struct BenchmarkResult {
    name: String,
    median_secs: f64,
    compressed_mbps: f64,
    decompressed_mbps: f64,
    decompressed_bytes: u64,
    line_count: u64,
}

fn main() -> Result<()> {
    let start = Instant::now();
    let cli_args = Cli::parse();
    let iterations = cli_args.iterations.unwrap_or(ITERATIONS);

    cli::banner("BGZF Benchmark");

    // ── Configuration ────────────────────────────────────
    cli::section("Configuration");

    let input_path = &cli_args.input;
    let file_size = std::fs::metadata(input_path)
        .with_context(|| format!("cannot stat {}", input_path.display()))?
        .len();

    cli::kv("Input", &input_path.display().to_string());
    cli::kv("File size", &perf::format_bytes(file_size));
    cli::kv("Iterations", &cli::num(iterations));

    if cli_args.flate2_only {
        cli::kv("Mode", "flate2-only");
    } else {
        cli::kv("Mode", "full suite");
    }

    // Detect flate2 backend
    let flate2_backend = detect_flate2_backend();
    cli::kv("flate2 backend", flate2_backend);

    eprintln!();

    // ── Warmup ───────────────────────────────────────────
    cli::section("Warmup");
    eprintln!("  {}", "Priming OS page cache...".dimmed());
    warmup_page_cache(input_path)?;
    cli::success("Page cache primed");
    eprintln!();

    // ── Benchmarks ───────────────────────────────────────
    cli::section("Benchmarks");

    let mut results: Vec<BenchmarkResult> = Vec::new();

    if cli_args.flate2_only {
        let name = format!("flate2 ({flate2_backend})");
        results.push(run_benchmark(&name, iterations, file_size, || {
            bench_flate2(input_path)
        })?);
    } else {
        // 1. flate2 MultiGzDecoder (baseline)
        let name = format!("flate2 ({flate2_backend})");
        results.push(run_benchmark(&name, iterations, file_size, || {
            bench_flate2(input_path)
        })?);

        // 2. noodles-bgzf with libdeflate
        results.push(run_benchmark(
            "noodles-bgzf (libdeflate, 1T)",
            iterations,
            file_size,
            || bench_noodles(input_path),
        )?);

        // 3. bgzip with libdeflater
        results.push(run_benchmark(
            "bgzip (libdeflater, 1T)",
            iterations,
            file_size,
            || bench_bgzip(input_path),
        )?);

        // 4. gzp BgzfSyncReader
        results.push(run_benchmark(
            "gzp (BgzfSyncReader)",
            iterations,
            file_size,
            || bench_gzp(input_path),
        )?);

        // 5. Manual BGZF + libdeflater (single-threaded)
        results.push(run_benchmark(
            "manual + libdeflater (1T)",
            iterations,
            file_size,
            || bench_manual_single(input_path),
        )?);

        // 6. noodles-bgzf MultithreadedReader
        results.push(run_benchmark(
            "noodles-bgzf (libdeflate, MT)",
            iterations,
            file_size,
            || bench_noodles_mt(input_path),
        )?);

        // 7. Manual BGZF + libdeflater + rayon
        results.push(run_benchmark(
            "manual + libdeflater + rayon",
            iterations,
            file_size,
            || bench_manual_rayon(input_path),
        )?);
    }

    eprintln!();

    // ── Validation ───────────────────────────────────────
    if results.len() > 1 {
        cli::section("Validation");
        let ref_bytes = results[0].decompressed_bytes;
        let ref_lines = results[0].line_count;
        let mut all_match = true;

        for r in &results[1..] {
            if r.decompressed_bytes != ref_bytes || r.line_count != ref_lines {
                cli::warning(&format!(
                    "{}: bytes={} lines={} (expected bytes={} lines={})",
                    r.name,
                    cli::num(r.decompressed_bytes),
                    cli::num(r.line_count),
                    cli::num(ref_bytes),
                    cli::num(ref_lines),
                ));
                all_match = false;
            }
        }

        if all_match {
            cli::success(&format!(
                "All approaches: {} decompressed, {} lines",
                perf::format_bytes(ref_bytes),
                cli::num(ref_lines),
            ));
        }

        eprintln!();
    }

    // ── Results ──────────────────────────────────────────
    cli::section("Results");
    print_results_table(&results);

    // ── Summary ──────────────────────────────────────────
    cli::print_summary(start);
    Ok(())
}

// ─── Benchmark runner ────────────────────────────────────

fn run_benchmark<F>(
    name: &str,
    iterations: usize,
    compressed_size: u64,
    f: F,
) -> Result<BenchmarkResult>
where
    F: Fn() -> Result<RunResult>,
{
    eprint!("  {:<36} ", name);

    let mut runs: Vec<RunResult> = Vec::with_capacity(iterations);
    for _ in 0..iterations {
        runs.push(f()?);
    }

    // Sort by elapsed to find median
    runs.sort_by(|a, b| a.elapsed.cmp(&b.elapsed));
    let median = &runs[iterations / 2];
    let median_secs = median.elapsed.as_secs_f64();

    let compressed_mbps = (compressed_size as f64 / (1024.0 * 1024.0)) / median_secs;
    let decompressed_mbps = (median.decompressed_bytes as f64 / (1024.0 * 1024.0)) / median_secs;

    eprintln!(
        "{:.2}s  ({:.0} MB/s compressed, {:.0} MB/s decompressed)",
        median_secs, compressed_mbps, decompressed_mbps
    );

    Ok(BenchmarkResult {
        name: name.to_string(),
        median_secs,
        compressed_mbps,
        decompressed_mbps,
        decompressed_bytes: median.decompressed_bytes,
        line_count: median.line_count,
    })
}

// ─── Benchmark implementations ───────────────────────────

/// 1. flate2 MultiGzDecoder — baseline
fn bench_flate2(path: &PathBuf) -> Result<RunResult> {
    let start = Instant::now();
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let decoder = flate2::bufread::MultiGzDecoder::new(reader);
    let (bytes, lines) = count_bytes_and_lines(decoder)?;
    Ok(RunResult {
        elapsed: start.elapsed(),
        decompressed_bytes: bytes,
        line_count: lines,
    })
}

/// 2. noodles-bgzf Reader with libdeflate
fn bench_noodles(path: &PathBuf) -> Result<RunResult> {
    use noodles_bgzf as bgzf;

    let start = Instant::now();
    let file = File::open(path)?;
    let reader = bgzf::io::Reader::new(file);
    let (bytes, lines) = count_bytes_and_lines(reader)?;
    Ok(RunResult {
        elapsed: start.elapsed(),
        decompressed_bytes: bytes,
        line_count: lines,
    })
}

/// 3. bgzip BGZFReader with libdeflater
fn bench_bgzip(path: &PathBuf) -> Result<RunResult> {
    let start = Instant::now();
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let decoder = bgzip::BGZFReader::new(reader)?;
    let (bytes, lines) = count_bytes_and_lines(decoder)?;
    Ok(RunResult {
        elapsed: start.elapsed(),
        decompressed_bytes: bytes,
        line_count: lines,
    })
}

/// 4. gzp BgzfSyncReader
fn bench_gzp(path: &PathBuf) -> Result<RunResult> {
    use gzp::bgzf::BgzfSyncReader;

    let start = Instant::now();
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let decoder = BgzfSyncReader::new(reader);
    let (bytes, lines) = count_bytes_and_lines(decoder)?;
    Ok(RunResult {
        elapsed: start.elapsed(),
        decompressed_bytes: bytes,
        line_count: lines,
    })
}

/// 5. Manual BGZF block parsing + libdeflater (single-threaded)
fn bench_manual_single(path: &PathBuf) -> Result<RunResult> {
    let start = Instant::now();
    let raw = std::fs::read(path)?;
    let blocks = parse_bgzf_blocks(&raw)?;

    let mut total_bytes: u64 = 0;
    let mut total_lines: u64 = 0;
    let mut decompressor = libdeflater::Decompressor::new();

    for block in &blocks {
        let payload = &raw[block.payload_offset..block.payload_offset + block.payload_len];
        let mut output = vec![0u8; block.uncompressed_size];
        decompressor
            .deflate_decompress(payload, &mut output)
            .map_err(|e| anyhow::anyhow!("deflate decompress failed: {e}"))?;
        total_bytes += output.len() as u64;
        total_lines += bytecount::count(&output, b'\n') as u64;
    }

    Ok(RunResult {
        elapsed: start.elapsed(),
        decompressed_bytes: total_bytes,
        line_count: total_lines,
    })
}

/// 6. noodles-bgzf MultithreadedReader
fn bench_noodles_mt(path: &PathBuf) -> Result<RunResult> {
    use noodles_bgzf as bgzf;

    let start = Instant::now();
    let file = File::open(path)?;
    let reader = bgzf::io::MultithreadedReader::new(file);
    let (bytes, lines) = count_bytes_and_lines(reader)?;
    Ok(RunResult {
        elapsed: start.elapsed(),
        decompressed_bytes: bytes,
        line_count: lines,
    })
}

/// 7. Manual BGZF + libdeflater + rayon parallel decompression
fn bench_manual_rayon(path: &PathBuf) -> Result<RunResult> {
    use rayon::prelude::*;

    let start = Instant::now();
    let raw = std::fs::read(path)?;
    let blocks = parse_bgzf_blocks(&raw)?;

    let results: Vec<(u64, u64)> = blocks
        .par_iter()
        .map(|block| {
            let mut decompressor = libdeflater::Decompressor::new();
            let payload = &raw[block.payload_offset..block.payload_offset + block.payload_len];
            let mut output = vec![0u8; block.uncompressed_size];
            decompressor
                .deflate_decompress(payload, &mut output)
                .expect("deflate decompress failed");
            let bytes = output.len() as u64;
            let lines = bytecount::count(&output, b'\n') as u64;
            (bytes, lines)
        })
        .collect();

    let (total_bytes, total_lines) = results
        .iter()
        .fold((0u64, 0u64), |(b, l), (rb, rl)| (b + rb, l + rl));

    Ok(RunResult {
        elapsed: start.elapsed(),
        decompressed_bytes: total_bytes,
        line_count: total_lines,
    })
}

// ─── BGZF block parser ──────────────────────────────────

struct BgzfBlock {
    payload_offset: usize,
    payload_len: usize,
    uncompressed_size: usize,
}

fn parse_bgzf_blocks(data: &[u8]) -> Result<Vec<BgzfBlock>> {
    let mut blocks = Vec::new();
    let mut offset = 0;

    while offset < data.len() {
        // An empty EOF block is exactly 28 bytes with ISIZE=0
        if offset + 18 > data.len() {
            break;
        }

        // Validate magic bytes
        if data[offset..offset + 4] != BGZF_MAGIC {
            bail!(
                "invalid BGZF magic at offset {offset}: {:02x} {:02x} {:02x} {:02x}",
                data[offset],
                data[offset + 1],
                data[offset + 2],
                data[offset + 3],
            );
        }

        // BSIZE is at bytes 16-17 (little-endian), block size - 1
        let bsize = u16::from_le_bytes([data[offset + 16], data[offset + 17]]) as usize + 1;

        if offset + bsize > data.len() {
            bail!(
                "BGZF block at offset {offset} extends beyond file (bsize={bsize}, remaining={})",
                data.len() - offset
            );
        }

        // ISIZE is the last 4 bytes of the block (little-endian)
        let isize_offset = offset + bsize - 4;
        let uncompressed_size = u32::from_le_bytes([
            data[isize_offset],
            data[isize_offset + 1],
            data[isize_offset + 2],
            data[isize_offset + 3],
        ]) as usize;

        // The raw DEFLATE payload starts after the gzip header.
        // Standard BGZF header is 18 bytes, payload ends 8 bytes before block end (CRC32 + ISIZE).
        let payload_offset = offset + 18;
        let payload_len = bsize - 18 - 8;

        // Skip empty EOF block
        if uncompressed_size > 0 {
            blocks.push(BgzfBlock {
                payload_offset,
                payload_len,
                uncompressed_size,
            });
        }

        offset += bsize;
    }

    Ok(blocks)
}

// ─── Helpers ─────────────────────────────────────────────

/// Read all bytes from a reader, counting decompressed bytes and newlines.
fn count_bytes_and_lines<R: Read>(reader: R) -> Result<(u64, u64)> {
    let mut buf_reader = BufReader::with_capacity(256 * 1024, reader);
    let mut total_bytes: u64 = 0;
    let mut total_lines: u64 = 0;

    loop {
        let buf = buf_reader.fill_buf()?;
        if buf.is_empty() {
            break;
        }
        let len = buf.len();
        total_bytes += len as u64;
        total_lines += bytecount::count(buf, b'\n') as u64;
        buf_reader.consume(len);
    }

    Ok((total_bytes, total_lines))
}

/// Read the file once to prime OS page cache.
fn warmup_page_cache(path: &PathBuf) -> Result<()> {
    let mut file = File::open(path)?;
    let mut buf = vec![0u8; 1024 * 1024];
    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
    }
    Ok(())
}

/// Detect which flate2 backend is compiled in.
fn detect_flate2_backend() -> &'static str {
    #[cfg(feature = "bench-zlib-ng")]
    {
        return "zlib-ng";
    }
    #[cfg(feature = "bench-zlib-rs")]
    {
        return "zlib-rs";
    }
    #[allow(unreachable_code)]
    "miniz_oxide"
}

/// Print formatted results table using comfy-table.
fn print_results_table(results: &[BenchmarkResult]) {
    if results.is_empty() {
        return;
    }

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
        Cell::new("Comp MB/s")
            .add_attribute(Attribute::Bold)
            .set_alignment(CellAlignment::Right),
        Cell::new("Decomp MB/s")
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
            Cell::new(format!("{:.2}", r.median_secs)).set_alignment(CellAlignment::Right),
            Cell::new(format!("{:.1}", r.compressed_mbps)).set_alignment(CellAlignment::Right),
            Cell::new(format!("{:.1}", r.decompressed_mbps)).set_alignment(CellAlignment::Right),
            Cell::new(format!("{speedup:.2}x")).set_alignment(CellAlignment::Right),
        ]);
    }

    // Print to stderr to match project conventions
    for line in table.to_string().lines() {
        eprintln!("  {line}");
    }
    eprintln!();
}
