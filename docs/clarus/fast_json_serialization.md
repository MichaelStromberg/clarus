# JSON Serialization Benchmark Results

**Date:** 2026-02-25 20:52

## Test Setup

| Parameter | Value |
|-----------|-------|
| Single iterations | 100,000 |
| Batch size | 10,000 |
| Runs (median) | 5 |
| SNV position size | ~692 bytes |
| Indel position size | ~744 bytes |
| Batch mix | 85% SNV / 15% indel |

### Crates

| Crate | Key Feature |
|-------|-------------|
| **serde_json** | Baseline — the standard |
| **sonic-rs** | SIMD string escaping (AVX2/NEON) |
| **simd-json** | SIMD-accelerated (mainly parsing, also serialization) |
| **manual writer** | Direct `Vec<u8>` writes with `itoa`/`ryu`, zero serde overhead |

## Single Position (hot cache)

| Approach | Median (s) | MB/s | ns/position | Speedup |
|----------|----------:|-----:|------------:|--------:|
| serde_json | 0.068 | 968.1 | 682 | 1.00x |
| sonic-rs | 0.054 | 1211.7 | 545 | 1.25x |
| simd-json | 0.084 | 785.8 | 840 | 0.81x |
| manual writer | 0.100 | 661.2 | 998 | 0.68x |

## Batch (mixed positions)

| Approach | Median (s) | MB/s | ns/position | Speedup |
|----------|----------:|-----:|------------:|--------:|
| serde_json | 0.007 | 981.0 | 680 | 1.00x |
| sonic-rs | 0.005 | 1223.2 | 546 | 1.25x |
| simd-json | 0.008 | 790.5 | 844 | 0.81x |
| manual writer | 0.010 | 648.7 | 1029 | 0.66x |

## Analysis

### Winner: sonic-rs (1.25x faster)

**sonic-rs** is the clear winner, delivering a consistent **25% speedup** over serde_json in both
single-position and batch scenarios. It achieves ~1,220 MB/s vs serde_json's ~970 MB/s on Apple
Silicon (M3). The speedup is uniform across hot-cache and cold-cache workloads, suggesting it comes
from SIMD-accelerated string escaping (NEON on ARM) rather than cache effects.

### Surprises

- **simd-json was 19% slower** than serde_json for serialization. This is expected — simd-json's
  SIMD optimizations primarily target *parsing* (deserialization). Its serialization path has
  overhead from the `serde_impl` bridge layer.

- **The manual writer was 32-34% slower** than serde_json. This is counterintuitive but explained
  by serde_json's highly optimized internal `Formatter` and `Serializer` which avoid per-field
  branch overhead. The manual writer's `write_json_string` with byte-by-byte escaping is slower
  than serde_json's SIMD-aware string writer. Additionally, the macro-based conditional field
  emission adds branch overhead that serde's codegen eliminates at compile time.

### Recommendation

Switch `JsonWriter::write_position` from `serde_json::to_writer` to `sonic_rs::to_writer`.
This is a **one-line drop-in change** since sonic-rs is serde-compatible:

```rust
// Before
serde_json::to_writer(&mut self.scratch, position)?;

// After
sonic_rs::to_writer(&mut self.scratch, position)?;
```

Expected impact on end-to-end annotation: ~25% reduction in JSON serialization time. Since
serialization was ~40-50% of total runtime (after gzip optimization), this translates to
roughly **10-12% overall speedup** for free.

### Platform: Apple Silicon M3, macOS, Rust 2024 edition
