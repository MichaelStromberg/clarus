//! Bounded parallel BGZF streaming reader.
//!
//! Instead of decompressing an entire BGZF file upfront (which requires holding
//! both compressed + decompressed data in memory), this reader:
//!
//! 1. Scans BGZF block headers from disk (~450 KB metadata for a 951 MB file)
//! 2. Decompresses blocks in batches: reads compressed data from disk, decompresses
//!    with rayon + libdeflater, then drops the compressed buffer
//! 3. Exposes a `BufRead` interface that yields decompressed data
//!
//! Memory at any time: one batch of compressed blocks (temporary, ~9 MB) + one batch
//! of decompressed blocks (~8 MB). The compressed file is never held in memory.

use std::fs::File;
use std::io::{BufRead, Read, Seek, SeekFrom};
use std::path::Path;

use rayon::prelude::*;

use crate::error::Error;

/// Default number of BGZF blocks to decompress per batch.
/// 128 blocks x 64 KB max = 8 MB of decompressed data per batch.
const DEFAULT_BATCH_SIZE: usize = 128;

/// BGZF magic bytes: gzip magic (0x1f, 0x8b) + deflate method (0x08) + FEXTRA flag (0x04).
const BGZF_MAGIC: [u8; 4] = [0x1f, 0x8b, 0x08, 0x04];

/// Gzip magic bytes (first two bytes only).
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Metadata for a single BGZF block (parsed from headers, no decompressed data).
struct BgzfBlockMeta {
    /// Absolute offset of the raw DEFLATE payload within the compressed file.
    payload_file_offset: u64,
    /// Length of the raw DEFLATE payload.
    payload_len: usize,
    /// Uncompressed size of this block (from gzip ISIZE trailer).
    uncompressed_size: usize,
}

/// A streaming BGZF reader that decompresses blocks in bounded parallel batches.
///
/// Implements `BufRead` so it can be used with `read_line()` and other standard
/// I/O patterns without loading the entire decompressed file into memory.
pub struct BgzfStreamReader {
    /// File handle for reading compressed blocks on demand (None for non-BGZF files).
    file: Option<File>,
    /// All block metadata parsed upfront from block headers.
    blocks: Vec<BgzfBlockMeta>,
    /// Index into `blocks` where the next batch starts.
    batch_start: usize,
    /// Decompressed block buffers for the current batch.
    batch_buffers: Vec<Vec<u8>>,
    /// Which block within the current batch we're currently reading from.
    batch_index: usize,
    /// Byte offset within the current block being read.
    byte_offset: usize,
    /// Number of blocks per batch.
    batch_size: usize,
}

impl BgzfStreamReader {
    /// Open a BGZF-compressed file and prepare for streaming decompression.
    ///
    /// Scans block headers from disk (never loads the entire compressed file),
    /// then decompresses the first batch. Also handles plain gzip and uncompressed files.
    pub fn new(path: &Path) -> Result<Self, Error> {
        let mut file = File::open(path)?;

        // Read first 4 bytes to detect format
        let mut magic = [0u8; 4];
        let bytes_read = file.read(&mut magic)?;
        file.seek(SeekFrom::Start(0))?;

        if bytes_read < 2 {
            return Err(Error::VcfFormat("file too small to be BGZF or gzip".into()));
        }

        if bytes_read >= 4 && magic == BGZF_MAGIC {
            // BGZF: scan metadata from disk, prepare for batch decompression
            let blocks = scan_block_metadata(&mut file)?;

            let mut reader = Self {
                file: Some(file),
                blocks,
                batch_start: 0,
                batch_buffers: Vec::new(),
                batch_index: 0,
                byte_offset: 0,
                batch_size: DEFAULT_BATCH_SIZE,
            };

            if !reader.blocks.is_empty() {
                reader.decompress_next_batch()?;
            }

            Ok(reader)
        } else if magic[..2] == GZIP_MAGIC {
            Self::from_plain_gzip(file)
        } else {
            Self::from_uncompressed(file)
        }
    }

    /// Decompress plain gzip from file stream (no compressed data held in memory).
    fn from_plain_gzip(file: File) -> Result<Self, Error> {
        let reader = std::io::BufReader::new(file);
        let mut decoder = flate2::bufread::MultiGzDecoder::new(reader);
        let mut output = Vec::new();
        decoder.read_to_end(&mut output)?;

        Ok(Self {
            file: None,
            blocks: Vec::new(),
            batch_start: 0,
            batch_buffers: vec![output],
            batch_index: 0,
            byte_offset: 0,
            batch_size: DEFAULT_BATCH_SIZE,
        })
    }

    /// Read uncompressed file into a single buffer.
    fn from_uncompressed(mut file: File) -> Result<Self, Error> {
        let mut data = Vec::new();
        file.read_to_end(&mut data)?;

        Ok(Self {
            file: None,
            blocks: Vec::new(),
            batch_start: 0,
            batch_buffers: vec![data],
            batch_index: 0,
            byte_offset: 0,
            batch_size: DEFAULT_BATCH_SIZE,
        })
    }

    /// Decompress the next batch of blocks using rayon + libdeflater.
    ///
    /// Reads compressed data for the batch from disk into a temporary buffer,
    /// decompresses in parallel, then drops the compressed buffer.
    fn decompress_next_batch(&mut self) -> Result<(), Error> {
        let end = (self.batch_start + self.batch_size).min(self.blocks.len());
        let batch_blocks = &self.blocks[self.batch_start..end];

        if batch_blocks.is_empty() {
            self.batch_buffers.clear();
            self.batch_index = 0;
            self.byte_offset = 0;
            return Ok(());
        }

        // Read the contiguous compressed range for this batch from disk
        let file = self.file.as_mut().expect("file handle required for BGZF");
        let range_start = batch_blocks[0].payload_file_offset;
        let last = batch_blocks.last().unwrap();
        let range_end = last.payload_file_offset + last.payload_len as u64;
        let range_len = (range_end - range_start) as usize;

        file.seek(SeekFrom::Start(range_start))?;
        let mut compressed = vec![0u8; range_len];
        file.read_exact(&mut compressed)?;

        // Pre-allocate output buffers
        let mut buffers: Vec<Vec<u8>> = batch_blocks
            .iter()
            .map(|b| vec![0u8; b.uncompressed_size])
            .collect();

        // Decompress in parallel
        buffers
            .par_iter_mut()
            .zip(batch_blocks.par_iter())
            .try_for_each(|(buf, block)| -> Result<(), Error> {
                let offset = (block.payload_file_offset - range_start) as usize;
                let payload = &compressed[offset..offset + block.payload_len];
                let mut decompressor = libdeflater::Decompressor::new();
                decompressor.deflate_decompress(payload, buf).map_err(|e| {
                    Error::VcfFormat(format!("BGZF deflate decompress failed: {e}"))
                })?;
                Ok(())
            })?;

        // compressed is dropped here — only decompressed buffers remain
        self.batch_buffers = buffers;
        self.batch_index = 0;
        self.byte_offset = 0;
        self.batch_start = end;

        Ok(())
    }
}

impl Read for BgzfStreamReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let available = self.fill_buf()?;
        if available.is_empty() {
            return Ok(0);
        }
        let n = buf.len().min(available.len());
        buf[..n].copy_from_slice(&available[..n]);
        self.consume(n);
        Ok(n)
    }
}

impl BufRead for BgzfStreamReader {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        // Phase 1: Advance past any exhausted blocks within the current batch
        while self.batch_index < self.batch_buffers.len()
            && self.byte_offset >= self.batch_buffers[self.batch_index].len()
        {
            self.batch_index += 1;
            self.byte_offset = 0;
        }

        // Phase 2: If current batch is exhausted, load the next batch
        if self.batch_index >= self.batch_buffers.len() && self.batch_start < self.blocks.len() {
            self.decompress_next_batch()
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))?;
        }

        // Phase 3: Return a reference to the current position (no more mutation)
        if self.batch_index < self.batch_buffers.len() {
            Ok(&self.batch_buffers[self.batch_index][self.byte_offset..])
        } else {
            Ok(&[])
        }
    }

    fn consume(&mut self, amt: usize) {
        self.byte_offset += amt;
    }
}

/// Scan BGZF block headers from a file to extract metadata without reading payloads.
///
/// Reads only the 18-byte header and 4-byte ISIZE trailer per block, seeking over
/// the payload data. Uses ~22 bytes of stack buffers regardless of file size.
fn scan_block_metadata(file: &mut File) -> Result<Vec<BgzfBlockMeta>, Error> {
    let file_len = file.metadata()?.len();
    file.seek(SeekFrom::Start(0))?;

    let mut blocks = Vec::new();
    let mut file_offset: u64 = 0;
    let mut header_buf = [0u8; 18];
    let mut isize_buf = [0u8; 4];

    while file_offset < file_len {
        if file_offset + 18 > file_len {
            break;
        }

        // Read block header (18 bytes)
        file.read_exact(&mut header_buf)?;

        if header_buf[..4] != BGZF_MAGIC {
            return Err(Error::VcfFormat(format!(
                "invalid BGZF magic at offset {file_offset}: {:02x} {:02x} {:02x} {:02x}",
                header_buf[0], header_buf[1], header_buf[2], header_buf[3],
            )));
        }

        // BSIZE is at bytes 16-17 (little-endian), block size - 1
        let bsize = u16::from_le_bytes([header_buf[16], header_buf[17]]) as u64 + 1;

        if file_offset + bsize > file_len {
            return Err(Error::VcfFormat(format!(
                "BGZF block at offset {file_offset} extends beyond file (bsize={bsize}, remaining={})",
                file_len - file_offset
            )));
        }

        // Seek past payload + CRC32 to ISIZE (last 4 bytes of block)
        // Current position: file_offset + 18. Target: file_offset + bsize - 4.
        // Skip: bsize - 18 - 4 = bsize - 22.
        file.seek(SeekFrom::Current(bsize as i64 - 22))?;
        file.read_exact(&mut isize_buf)?;
        // Now positioned at the start of the next block

        let uncompressed_size = u32::from_le_bytes(isize_buf) as usize;

        // Raw DEFLATE payload: after 18-byte gzip header, before 8-byte trailer
        let payload_file_offset = file_offset + 18;
        let payload_len = bsize as usize - 18 - 8;

        // Skip empty EOF block
        if uncompressed_size > 0 {
            blocks.push(BgzfBlockMeta {
                payload_file_offset,
                payload_len,
                uncompressed_size,
            });
        }

        file_offset += bsize;
    }

    Ok(blocks)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_bgzf_produces_empty_reader() {
        // A valid BGZF EOF block (28 bytes)
        let eof_block: [u8; 28] = [
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.bgzf");
        std::fs::write(&path, &eof_block).unwrap();

        let mut reader = BgzfStreamReader::new(&path).unwrap();
        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).unwrap();
        assert!(buf.is_empty());
    }

    #[test]
    fn invalid_magic_returns_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("bad.bgzf");
        // Valid EOF block followed by invalid data
        let mut data = Vec::new();
        data.extend_from_slice(&[
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ]);
        // Bad magic (20 bytes so it passes the length check)
        data.extend_from_slice(&[0x00u8; 20]);

        std::fs::write(&path, &data).unwrap();
        let result = BgzfStreamReader::new(&path);
        assert!(result.is_err());
    }

    #[test]
    fn scan_metadata_skips_eof_blocks() {
        let eof_block: [u8; 28] = [
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("eof.bgzf");
        std::fs::write(&path, &eof_block).unwrap();

        let mut file = File::open(&path).unwrap();
        let blocks = scan_block_metadata(&mut file).unwrap();
        assert!(blocks.is_empty());
    }
}
