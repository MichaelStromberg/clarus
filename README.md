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
| `error` | Central error type (`thiserror`-based) with `Io`, `Parse`, `Validation`, and `Format` variants |
| `chromosome` | Chromosome representation with UCSC, Ensembl, RefSeq, and GenBank naming |
| `genome_assembly` | Genome assembly enum (`GRCh37`, `GRCh38`) and deterministic reference ID computation |
| `assembly_report` | Parser for NCBI genome assembly report files |
| `fasta` | Parser for gzip-compressed FASTA sequence files |
| `bands` | Parser for NCBI cytogenetic band (ideogram) files |
| `mirna` | Parser for miRNA region annotations from gzip-compressed GFF3 files |
| `perf` | Performance measurement utilities (elapsed time, peak memory) |
| `reference` | Binary reference file format: reader, writer, and shared I/O traits |

### Reference File Format

Clarus uses a custom binary format for reference sequences:

- 8-byte magic signature followed by file type and format version
- Per-chromosome data blocks containing Zstandard-compressed sequences (level 21), cytogenetic bands, and miRNA regions
- Header carries genome assembly, patch level, and a SHA-256-derived 32-bit reference ID for compatibility checks
- All coordinates are 1-based, closed intervals

### Binaries (`src/bin/`)

**`create_ref`** — Creates a Clarus reference sequence file from NCBI data sources.

```
Usage: create_ref [OPTIONS] -a <ASSEMBLY_REPORT> -f <FASTA> -o <OUT>

Options:
  -a, --assembly-report <PATH>   NCBI genome assembly report file
  -f, --fasta <PATH>             Gzip-compressed FASTA file
  -b, --bands <PATH>             NCBI GDP ideogram file (cytogenetic bands)
      --gff <PATH>               Gzip-compressed GFF3 file (miRNA regions)
  -o, --out <PATH>               Output reference file path
  -g, --assembly <NAME>          Genome assembly name [default: GRCh38]
  -p, --patch-level <N>          Assembly patch level [default: 14]
```

The tool validates input data (sequence lengths, chrY PAR masking, ref_index contiguity) and round-trip verifies the written file by reading it back.

## Current Status

**Early development.** The reference file pipeline (`create_ref`) is implemented and functional. The runtime annotation engine is not yet started.

### Implemented

- NCBI assembly report parsing
- FASTA sequence parsing and compression
- Cytogenetic band parsing (NCBI ideogram format)
- miRNA region parsing (GFF3 format)
- Binary reference file writing and reading with round-trip verification
- Benchmarks for reference file header reading and chromosome decompression

## License

GPL-3.0
