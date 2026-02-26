//! Streaming gzip JSON writer.
//!
//! Writes the JSON output incrementally: header first, then one position at a time,
//! so memory usage is proportional to a single position, not the entire output.
//!
//! Each position is serialized to a reusable scratch buffer before writing to gzip,
//! so the compressor receives large chunks instead of tiny per-field writes.

use std::io::{BufWriter, Write};
use std::path::Path;

use flate2::Compression;
use flate2::write::GzEncoder;

use crate::error::Error;
use crate::json::types::{Header, Position};

/// Output BufWriter capacity (256 KB for fewer syscalls).
const BUF_CAPACITY: usize = 256 * 1024;

/// A streaming JSON writer that produces gzip-compressed output.
pub struct JsonWriter {
    encoder: GzEncoder<BufWriter<std::fs::File>>,
    position_count: usize,
    /// Reusable scratch buffer for serializing each position before passing to gzip.
    scratch: Vec<u8>,
}

impl JsonWriter {
    /// Create a new JSON writer that writes to a gzip-compressed file.
    pub fn new(path: &Path) -> Result<Self, Error> {
        let file = std::fs::File::create(path)?;
        let buf_writer = BufWriter::with_capacity(BUF_CAPACITY, file);
        let encoder = GzEncoder::new(buf_writer, Compression::fast());

        Ok(Self {
            encoder,
            position_count: 0,
            scratch: Vec::with_capacity(4096),
        })
    }

    /// Write the JSON header and open the positions array.
    pub fn write_header(&mut self, header: &Header) -> Result<(), Error> {
        let header_json = sonic_rs::to_string(header)
            .map_err(|e| Error::Format(format!("JSON serialization error: {e}")))?;

        write!(self.encoder, "{{\"header\":{header_json},\"positions\":[")?;
        Ok(())
    }

    /// Write a single position to the output.
    pub fn write_position(&mut self, position: &Position) -> Result<(), Error> {
        // Serialize to scratch buffer first, then write the whole chunk to gzip.
        // This avoids thousands of tiny writes per position through the compressor.
        self.scratch.clear();
        if self.position_count > 0 {
            self.scratch.push(b',');
        }
        sonic_rs::to_writer(&mut self.scratch, position)
            .map_err(|e| Error::Format(format!("JSON serialization error: {e}")))?;
        self.encoder.write_all(&self.scratch)?;

        self.position_count += 1;
        Ok(())
    }

    /// Close the positions array, close the top-level object, and flush.
    pub fn finish(mut self) -> Result<(), Error> {
        writeln!(self.encoder, "]}}")?;
        self.encoder.finish()?;
        Ok(())
    }
}
