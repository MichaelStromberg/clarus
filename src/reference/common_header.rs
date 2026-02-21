//! Common file header shared across Clarus file types.

use std::io::{Read, Write};

use crate::error::Error;
use crate::reference::binary_io::{BinaryRead, BinaryWrite};

/// Clarus file signature: 0x89 0x43 0x4C 0x41 0x0D 0x0A 0x1A 0x0A
pub const CLARUS_SIGNATURE: u64 = 727_905_341_820_126_089;

/// File type identifier for reference sequence files.
pub const REFERENCE_FILE_TYPE: u16 = 1;

/// Current format version for reference sequence files.
pub const REFERENCE_FORMAT_VERSION: u16 = 1;

/// Writes the common header (signature + file type + format version) to a writer.
pub fn write_common_header<W: Write>(
    writer: &mut W,
    file_type: u16,
    format_version: u16,
) -> Result<(), Error> {
    writer.write_u64(CLARUS_SIGNATURE)?;
    writer.write_u16(file_type)?;
    writer.write_u16(format_version)?;
    Ok(())
}

/// Reads and validates the common header from a reader.
/// Returns (file_type, format_version).
pub fn read_common_header<R: Read>(reader: &mut R) -> Result<(u16, u16), Error> {
    let signature = reader.read_u64()?;
    if signature != CLARUS_SIGNATURE {
        return Err(Error::Format(format!(
            "invalid Clarus file signature: expected {CLARUS_SIGNATURE}, got {signature}"
        )));
    }

    let file_type = reader.read_u16()?;
    let format_version = reader.read_u16()?;

    Ok((file_type, format_version))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn signature_bytes_match_spec() {
        let bytes = CLARUS_SIGNATURE.to_le_bytes();
        assert_eq!(bytes, [0x89, 0x43, 0x4C, 0x41, 0x0D, 0x0A, 0x1A, 0x0A]);
    }

    #[test]
    fn round_trip() {
        let mut buf = Vec::new();
        write_common_header(&mut buf, REFERENCE_FILE_TYPE, REFERENCE_FORMAT_VERSION).unwrap();

        let mut cursor = Cursor::new(buf);
        let (file_type, format_version) = read_common_header(&mut cursor).unwrap();
        assert_eq!(file_type, REFERENCE_FILE_TYPE);
        assert_eq!(format_version, REFERENCE_FORMAT_VERSION);
    }

    #[test]
    fn invalid_signature() {
        let mut buf = Vec::new();
        buf.extend_from_slice(&0u64.to_le_bytes());
        buf.extend_from_slice(&1u16.to_le_bytes());
        buf.extend_from_slice(&1u16.to_le_bytes());

        let mut cursor = Cursor::new(buf);
        let result = read_common_header(&mut cursor);
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("invalid Clarus file signature")
        );
    }
}
