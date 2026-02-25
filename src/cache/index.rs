//! Cache index file writer.
//!
//! The index file provides random access into the transcript cache file,
//! mapping each chromosome to a byte offset.

use std::io::{Seek, Write};

use crate::cache::{INDEX_FILE_TYPE, INDEX_FORMAT_VERSION, IndexEntry};
use crate::error::Error;
use crate::reference::binary_io::BinaryWrite;
use crate::reference::common_header::write_common_header;

/// Write a cache index file.
///
/// Format:
///   Common header: signature + file_type(3) + format_version(1) + file_length(u64)
///   Index header: cache_id(u32) + chromosome_count(u16)
///   Per reference: ref_index(u16) + file_offset(u64)
pub fn write_index<W: Write + Seek>(
    writer: &mut W,
    cache_id: u32,
    index_entries: &[IndexEntry],
) -> Result<(), Error> {
    // Header (file_length placeholder = 0, will be patched)
    write_common_header(writer, INDEX_FILE_TYPE, INDEX_FORMAT_VERSION, 0)?;

    // Index header
    writer.write_u32(cache_id)?;
    writer.write_u16(
        u16::try_from(index_entries.len())
            .map_err(|_| Error::Validation("index entry count exceeds u16::MAX".into()))?,
    )?;

    // Per-reference entries
    for &(ref_index, file_offset) in index_entries {
        writer.write_u16(ref_index)?;
        writer.write_u64(file_offset)?;
    }

    // Patch file_length at offset 12
    let file_length = writer.stream_position()?;
    writer.seek(std::io::SeekFrom::Start(12))?;
    writer.write_u64(file_length)?;
    writer.seek(std::io::SeekFrom::End(0))?;

    Ok(())
}
