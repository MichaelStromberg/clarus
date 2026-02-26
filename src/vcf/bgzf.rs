//! Parallel BGZF decompression for VCF files.
//!
//! Reads a BGZF-compressed file into memory, parses block headers to find
//! DEFLATE payloads and their uncompressed sizes, then decompresses all blocks
//! in parallel via rayon with per-thread `libdeflater::Decompressor`.
//! Falls back to `flate2::MultiGzDecoder` for plain gzip files.

use std::io::{BufReader, Read};
use std::path::Path;

use rayon::prelude::*;

use crate::error::Error;

/// Wrapper around a raw pointer to allow sending across threads.
/// SAFETY: The caller must guarantee non-overlapping access to the pointed-to memory.
struct SendPtr(*mut u8);
unsafe impl Send for SendPtr {}
unsafe impl Sync for SendPtr {}

/// BGZF magic bytes: gzip magic (0x1f, 0x8b) + deflate method (0x08) + FEXTRA flag (0x04).
const BGZF_MAGIC: [u8; 4] = [0x1f, 0x8b, 0x08, 0x04];

/// Gzip magic bytes (first two bytes only).
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

struct BgzfBlock {
    payload_offset: usize,
    payload_len: usize,
    uncompressed_size: usize,
    /// Cumulative output offset for this block (prefix sum of uncompressed sizes).
    output_offset: usize,
}

/// Read a BGZF or plain gzip file and return the fully decompressed bytes.
pub fn read_vcf_bytes(path: &Path) -> Result<Vec<u8>, Error> {
    let raw = std::fs::read(path)?;

    if raw.len() < 4 {
        return Err(Error::VcfFormat("file too small to be BGZF or gzip".into()));
    }

    // Check for BGZF magic (gzip with FEXTRA flag)
    if raw[..4] == BGZF_MAGIC {
        decompress_bgzf(&raw)
    } else if raw[..2] == GZIP_MAGIC {
        decompress_plain_gzip(&raw)
    } else {
        // Assume uncompressed
        Ok(raw)
    }
}

/// Parse BGZF blocks and decompress all in parallel using rayon + libdeflater.
fn decompress_bgzf(data: &[u8]) -> Result<Vec<u8>, Error> {
    let blocks = parse_bgzf_blocks(data)?;

    if blocks.is_empty() {
        return Ok(Vec::new());
    }

    // Total output size is the last block's offset + its uncompressed size
    let last = blocks.last().unwrap();
    let total_size = last.output_offset + last.uncompressed_size;
    let mut output = vec![0u8; total_size];

    // Decompress blocks in parallel, writing directly to correct output offsets.
    // SAFETY: Each block writes to output[output_offset..output_offset+uncompressed_size].
    // These ranges are non-overlapping because output_offset is a prefix sum of
    // uncompressed sizes computed during block parsing.
    let output_ptr = SendPtr(output.as_mut_ptr());
    let ptr = &output_ptr; // force closure to capture &SendPtr, not the inner *mut u8

    blocks
        .par_iter()
        .try_for_each(|block| -> Result<(), Error> {
            let payload = &data[block.payload_offset..block.payload_offset + block.payload_len];
            let dest = unsafe {
                std::slice::from_raw_parts_mut(
                    ptr.0.add(block.output_offset),
                    block.uncompressed_size,
                )
            };
            let mut decompressor = libdeflater::Decompressor::new();
            decompressor
                .deflate_decompress(payload, dest)
                .map_err(|e| Error::VcfFormat(format!("BGZF deflate decompress failed: {e}")))?;
            Ok(())
        })?;

    Ok(output)
}

/// Fallback: decompress plain gzip using flate2.
fn decompress_plain_gzip(data: &[u8]) -> Result<Vec<u8>, Error> {
    let reader = BufReader::new(data);
    let mut decoder = flate2::bufread::MultiGzDecoder::new(reader);
    let mut output = Vec::new();
    decoder.read_to_end(&mut output)?;
    Ok(output)
}

/// Parse BGZF block headers to extract payload locations and uncompressed sizes.
fn parse_bgzf_blocks(data: &[u8]) -> Result<Vec<BgzfBlock>, Error> {
    let mut blocks = Vec::new();
    let mut offset = 0;
    let mut output_offset = 0;

    while offset < data.len() {
        if offset + 18 > data.len() {
            break;
        }

        if data[offset..offset + 4] != BGZF_MAGIC {
            return Err(Error::VcfFormat(format!(
                "invalid BGZF magic at offset {offset}: {:02x} {:02x} {:02x} {:02x}",
                data[offset],
                data[offset + 1],
                data[offset + 2],
                data[offset + 3],
            )));
        }

        // BSIZE is at bytes 16-17 (little-endian), block size - 1
        let bsize = u16::from_le_bytes([data[offset + 16], data[offset + 17]]) as usize + 1;

        if offset + bsize > data.len() {
            return Err(Error::VcfFormat(format!(
                "BGZF block at offset {offset} extends beyond file (bsize={bsize}, remaining={})",
                data.len() - offset
            )));
        }

        // ISIZE is the last 4 bytes of the block (little-endian)
        let isize_offset = offset + bsize - 4;
        let uncompressed_size = u32::from_le_bytes([
            data[isize_offset],
            data[isize_offset + 1],
            data[isize_offset + 2],
            data[isize_offset + 3],
        ]) as usize;

        // Raw DEFLATE payload: after 18-byte gzip header, before 8-byte trailer (CRC32 + ISIZE)
        let payload_offset = offset + 18;
        let payload_len = bsize - 18 - 8;

        // Skip empty EOF block
        if uncompressed_size > 0 {
            blocks.push(BgzfBlock {
                payload_offset,
                payload_len,
                uncompressed_size,
                output_offset,
            });
            output_offset += uncompressed_size;
        }

        offset += bsize;
    }

    Ok(blocks)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_bgzf_produces_empty_output() {
        // A valid BGZF EOF block (28 bytes)
        let eof_block: [u8; 28] = [
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];
        let result = decompress_bgzf(&eof_block).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn invalid_magic_returns_error() {
        // 20 bytes so we get past the length check and actually test magic validation
        let data = [0x00u8; 20];
        let result = decompress_bgzf(&data);
        assert!(result.is_err());
    }
}
