# Clarus

High-performance genomic variant annotation and VUS resolution engine written in Rust.

## Building

Requires Rust 2024 edition (1.85+).

```bash
cargo build              # Debug build
cargo build --release    # Release build
```

## Testing

```bash
cargo test               # Run all tests
cargo test <test_name>   # Run a single test
```

## Linting and Formatting

```bash
cargo clippy             # Lint
cargo fmt                # Format code
cargo fmt -- --check     # Check formatting without modifying
```

## Architecture

### Library Modules (`src/`)

| Module | Description |
|--------|-------------|
| `error` | Central error type (`thiserror`-based) with `Io`, `Parse`, `Validation`, `Format`, and `UnresolvableTranscript` variants |
| `chromosome` | Chromosome representation with UCSC, Ensembl, RefSeq, and GenBank naming |
| `genome_assembly` | Genome assembly enum (`GRCh37`, `GRCh38`) and deterministic reference ID computation |
| `config` | JSON configuration file deserialization and validation for pipeline inputs |
| `download` | HTTP file downloading with MD5 checksum verification and parallel execution |
| `assembly_report` | Parser for NCBI genome assembly report files |
| `fasta` | Parser for gzip-compressed FASTA sequence files with transcript ID extraction |
| `bands` | Parser for NCBI cytogenetic band (ideogram) files |
| `biotype` | Enum of 102 Sequence Ontology biotypes with categorization and single-byte serialization |
| `codon` | Codon translation tables (standard and vertebrate mitochondrial genetic codes) |
| `evaluation` | 7-level transcript evaluation hierarchy validating CDS translation against expected protein |
| `genbank` | Parser for gzip-compressed GenBank flat files extracting CDS records |
| `gff3` | GFF3 annotation parser with gene-transcript-exon hierarchy builder and CIGAR support |
| `mirna` | Parser for miRNA region annotations from gzip-compressed GFF3 files |
| `perf` | Performance measurement utilities (elapsed time, peak memory) |
| `reference` | Binary reference file format: reader, writer, and shared I/O traits |
| `sequence` | RNA and protein sequence dictionaries indexed by transcript/protein ID |
| `strand` | Forward/reverse strand type with GFF3 parsing and byte serialization |
| `transcript` | Transcript region construction, coordinate mapping, and intermediate data types |

### Reference File Format

Clarus uses a custom binary format for reference sequences:

- 8-byte magic signature followed by file type and format version
- Per-chromosome data blocks containing Zstandard-compressed sequences (level 21), cytogenetic bands, and miRNA regions
- Header carries genome assembly, patch level, and a SHA-256-derived 32-bit reference ID for compatibility checks
- All coordinates are 1-based, closed intervals

### Binaries (`src/bin/`)

**`create_ref`** — Creates a Clarus reference sequence file from NCBI data sources.

```
Usage: create_ref -c <CONFIG> -o <OUT_DIR>

Options:
  -c, --config <PATH>   JSON configuration file (URLs, MD5 checksums, assembly info)
  -o, --out <PATH>      Output data directory (writes to {out}/reference/{assembly}.p{patch}.dat)
```

Reads a JSON configuration file that declares source URLs and expected MD5 checksums, automatically downloads missing files (with parallel downloads and checksum verification), then parses, validates, and round-trip verifies the reference file.

Example configuration files are in `configurations/reference/`.

**`create_cache`** — Builds transcript cache files from RefSeq and Ensembl data sources.

```
Usage: create_cache -c <CONFIG> -r <REFERENCE> -o <OUT_DIR>

Options:
  -c, --config <PATH>   JSON configuration file (Ensembl/RefSeq versions, URLs, checksums)
  -r, --ref <PATH>      Reference sequence file (from create_ref)
  -o, --out <PATH>      Output data directory
```

Pipeline stages:
1. **File resolution** — resolves 13 data sources (8 Ensembl + 4 RefSeq + HGNC), downloading missing files with MD5 verification
2. **RefSeq parsing** — GFF3 annotations (genes + regulatory regions), GenBank CDS records, RNA FASTA sequences, and protein FASTA sequences
3. **Transcript construction** — builds transcript regions with exon normalization, intron insertion, CIGAR-aware coordinate mapping, and coding region detection
4. **Evaluation** — validates CDS translation against expected protein using a 7-level resolution hierarchy (perfect match, contained, AA edits, translational slippage, frame correction +1/+2, unresolvable)

For GRCh38: processes 105,006 transcripts (78,422 coding, 26,584 non-coding) with all coding transcripts resolved successfully.

Example configuration files are in `configurations/cache/`.

### Benchmarks (`benches/`)

**`bgzf_benchmark`** — Benchmarks BGZF decompression approaches to find the fastest strategy for VCF parsing.

```
Usage: cargo bench --bench bgzf_benchmark -- -i <INPUT> [--flate2-only] [-n <ITERATIONS>]

Options:
  -i, --input <PATH>       BGZF-compressed file (e.g., .vcf.gz)
  --flate2-only            Only run the flate2 benchmark (for testing different backends)
  -n, --iterations <N>     Iterations per approach (default: 5)
```

Compares 7 approaches: flate2 (with swappable miniz_oxide/zlib-rs/zlib-ng backends), noodles-bgzf, bgzip, gzp, and manual BGZF parsing with libdeflater — both single-threaded and parallel (rayon). Feature flags switch the flate2 backend:

```bash
cargo bench --bench bgzf_benchmark -- -i file.vcf.gz                          # full suite
cargo bench --bench bgzf_benchmark --features bench-zlib-rs -- -i file.vcf.gz --flate2-only
cargo bench --bench bgzf_benchmark --features bench-zlib-ng -- -i file.vcf.gz --flate2-only
```

Results on a 951 MB 3-sample WGS VCF (Apple M3 Max, 14 cores): manual BGZF + libdeflater + rayon achieves **11.4x** speedup over flate2 baseline (0.42s vs 4.76s). Full results in [`docs/clarus/fast_bgzf_parsing.md`](docs/clarus/fast_bgzf_parsing.md).

## Current Status

**Active development.** The build-time pipeline is functional: `create_ref` produces reference files and `create_cache` constructs transcript caches with full protein validation. The runtime annotation engine is not yet started.

### Implemented

- Config-driven pipelines with JSON configuration files
- Automatic downloading of source files with MD5 verification
- Binary reference file creation with round-trip verification (`create_ref`)
- GFF3 annotation parsing with gene-transcript-exon hierarchy
- GenBank flat file parsing for CDS records
- Transcript construction with CIGAR-aware coordinate mapping for frameshift transcripts
- 7-level protein evaluation resolving all coding transcripts for GRCh38
- 102 Sequence Ontology biotype classifications
- Standard and vertebrate mitochondrial codon translation
- Selenoprotein detection for relaxed evaluation thresholds
- BGZF decompression benchmarking (7 approaches, 3 flate2 backends) with parallel rayon strategy selected

## License

GPL-3.0
