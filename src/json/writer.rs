//! Three-file output writer: metadata, variants JSONL, and genes JSONL.
//!
//! - `{prefix}_variants.jsonl.gz` — BGZF-compressed, one JSON position per line
//! - `{prefix}_metadata.json.gz` — gzip-compressed metadata object
//! - `{prefix}_genes.jsonl.gz` — BGZF-compressed, one JSON gene per line (empty for now)
//!
//! Callers must invoke [`OutputWriter::finish`] to finalize all output streams.
//! If dropped without calling `finish`, the `Drop` impl will attempt a best-effort
//! cleanup (writing BGZF EOF markers and flushing buffers), but I/O errors during
//! drop are silently ignored.

use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use flate2::Compression;
use flate2::write::GzEncoder;
use noodles_bgzf::io::writer::{Builder, CompressionLevel};

use crate::error::Error;
use crate::json::types::{Metadata, Position};

/// Output BufWriter capacity (256 KB for fewer syscalls).
const BUF_CAPACITY: usize = 256 * 1024;

/// Manages the three output streams.
///
/// The variants file is opened immediately on construction. The metadata and genes
/// files are written when [`finish`](Self::finish) is called, since metadata depends
/// on post-annotation statistics.
pub struct OutputWriter {
    metadata_path: PathBuf,
    /// `None` after `finish()` has been called (taken by value).
    variants_writer: Option<noodles_bgzf::io::Writer<BufWriter<std::fs::File>>>,
    genes_path: PathBuf,
    /// Reusable scratch buffer for serializing each position before passing to BGZF.
    scratch: Vec<u8>,
}

impl OutputWriter {
    /// Create a new output writer. Opens `{prefix}_variants.jsonl.gz` immediately.
    /// Stores paths for metadata and genes files (written on `finish()`).
    pub fn new(prefix: &Path) -> Result<Self, Error> {
        let prefix_str = prefix.display().to_string();

        let variants_path = PathBuf::from(format!("{prefix_str}_variants.jsonl.gz"));
        let metadata_path = PathBuf::from(format!("{prefix_str}_metadata.json.gz"));
        let genes_path = PathBuf::from(format!("{prefix_str}_genes.jsonl.gz"));

        let file = std::fs::File::create(&variants_path)?;
        let buf_writer = BufWriter::with_capacity(BUF_CAPACITY, file);
        let variants_writer = Builder::default()
            .set_compression_level(CompressionLevel::FAST)
            .build_from_writer(buf_writer);

        Ok(Self {
            metadata_path,
            variants_writer: Some(variants_writer),
            genes_path,
            scratch: Vec::with_capacity(4096),
        })
    }

    /// Write a single position as one JSONL line to the variants file.
    pub fn write_position(&mut self, position: &Position) -> Result<(), Error> {
        self.scratch.clear();
        sonic_rs::to_writer(&mut self.scratch, position)
            .map_err(|e| Error::Format(format!("JSON serialization error: {e}")))?;
        self.scratch.push(b'\n');
        self.variants_writer
            .as_mut()
            .expect("write_position called after finish")
            .write_all(&self.scratch)?;
        Ok(())
    }

    /// Finish all output streams:
    /// 1. Finish the BGZF variants stream (writes EOF marker) and flush.
    /// 2. Write metadata as gzip-compressed JSON.
    /// 3. Write an empty genes BGZF file (just the EOF marker).
    pub fn finish(mut self, metadata: &Metadata) -> Result<(), Error> {
        // 1. Finish variants BGZF stream and flush the underlying BufWriter
        let variants_bgzf = self
            .variants_writer
            .take()
            .expect("finish called more than once");
        flush_bgzf(variants_bgzf)?;

        // 2. Write metadata (gzip, default compression) — serialize to scratch then compress
        let meta_file = std::fs::File::create(&self.metadata_path)?;
        let buf_writer = BufWriter::with_capacity(BUF_CAPACITY, meta_file);
        let mut encoder = GzEncoder::new(buf_writer, Compression::default());
        let mut meta_buf = std::mem::take(&mut self.scratch);
        meta_buf.clear();
        sonic_rs::to_writer(&mut meta_buf, metadata)
            .map_err(|e| Error::Format(format!("JSON serialization error: {e}")))?;
        meta_buf.push(b'\n');
        encoder.write_all(&meta_buf)?;
        // GzEncoder::finish() returns the inner BufWriter — flush it explicitly
        encoder.finish()?.flush()?;

        // 3. Write empty genes BGZF file (just EOF marker for now)
        let genes_file = std::fs::File::create(&self.genes_path)?;
        let genes_buf = BufWriter::with_capacity(BUF_CAPACITY, genes_file);
        let genes_writer = Builder::default()
            .set_compression_level(CompressionLevel::FAST)
            .build_from_writer(genes_buf);
        flush_bgzf(genes_writer)?;

        Ok(())
    }
}

impl Drop for OutputWriter {
    fn drop(&mut self) {
        if let Some(writer) = self.variants_writer.take() {
            // Best-effort: finish the BGZF stream so the file has a valid EOF marker.
            // Errors are silently ignored since Drop cannot propagate them.
            let _ = flush_bgzf(writer);
        }
    }
}

/// Finish a BGZF writer (writes EOF marker) and flush the underlying BufWriter.
fn flush_bgzf(writer: noodles_bgzf::io::Writer<BufWriter<std::fs::File>>) -> Result<(), Error> {
    let mut buf_writer = writer.finish()?;
    buf_writer.flush()?;
    Ok(())
}
