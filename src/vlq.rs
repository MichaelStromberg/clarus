//! Variable-length quantity (VLQ) integer encoding for cache binary format.
//!
//! Uses 7 data bits per byte with the MSB as a continuation flag.
//! Values 0-127 use 1 byte, 128-16383 use 2 bytes, etc.

use std::io::{Read, Write};

use crate::error::Error;

/// Write a VLQ-encoded unsigned integer.
pub fn write_vlq<W: Write>(writer: &mut W, mut value: u64) -> Result<(), Error> {
    loop {
        let mut byte = (value & 0x7F) as u8;
        value >>= 7;
        if value != 0 {
            byte |= 0x80;
        }
        writer.write_all(&[byte])?;
        if value == 0 {
            break;
        }
    }
    Ok(())
}

/// Read a VLQ-encoded unsigned integer.
pub fn read_vlq<R: Read>(reader: &mut R) -> Result<u64, Error> {
    let mut result: u64 = 0;
    let mut shift: u32 = 0;
    loop {
        let mut buf = [0u8; 1];
        reader.read_exact(&mut buf)?;
        let byte = buf[0];
        result |= ((byte & 0x7F) as u64) << shift;
        if byte & 0x80 == 0 {
            break;
        }
        shift += 7;
        if shift >= 64 {
            return Err(Error::Parse("VLQ integer overflow".to_string()));
        }
    }
    Ok(result)
}

/// Write a VLQ-prefixed UTF-8 string.
pub fn write_vlq_string<W: Write>(writer: &mut W, s: &str) -> Result<(), Error> {
    write_vlq(writer, s.len() as u64)?;
    writer.write_all(s.as_bytes())?;
    Ok(())
}

/// Read a VLQ-prefixed UTF-8 string.
pub fn read_vlq_string<R: Read>(reader: &mut R) -> Result<String, Error> {
    let len = read_vlq(reader)? as usize;
    if len == 0 {
        return Ok(String::new());
    }
    let mut buf = vec![0u8; len];
    reader.read_exact(&mut buf)?;
    String::from_utf8(buf).map_err(|e| Error::Parse(format!("invalid UTF-8 in VLQ string: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn round_trip(value: u64) -> u64 {
        let mut buf = Vec::new();
        write_vlq(&mut buf, value).unwrap();
        let mut cursor = Cursor::new(buf);
        read_vlq(&mut cursor).unwrap()
    }

    #[test]
    fn single_byte_values() {
        for v in 0..=127 {
            assert_eq!(round_trip(v), v);
        }
    }

    #[test]
    fn single_byte_encoding_size() {
        let mut buf = Vec::new();
        write_vlq(&mut buf, 0).unwrap();
        assert_eq!(buf.len(), 1);
        assert_eq!(buf[0], 0);

        buf.clear();
        write_vlq(&mut buf, 127).unwrap();
        assert_eq!(buf.len(), 1);
        assert_eq!(buf[0], 127);
    }

    #[test]
    fn two_byte_values() {
        assert_eq!(round_trip(128), 128);
        assert_eq!(round_trip(16383), 16383);
        assert_eq!(round_trip(300), 300);
    }

    #[test]
    fn two_byte_encoding_size() {
        let mut buf = Vec::new();
        write_vlq(&mut buf, 128).unwrap();
        assert_eq!(buf.len(), 2);
    }

    #[test]
    fn multi_byte_values() {
        assert_eq!(round_trip(16384), 16384);
        assert_eq!(round_trip(2_097_151), 2_097_151);
        assert_eq!(round_trip(u32::MAX as u64), u32::MAX as u64);
        assert_eq!(round_trip(u64::MAX), u64::MAX);
    }

    #[test]
    fn string_round_trip() {
        let mut buf = Vec::new();
        write_vlq_string(&mut buf, "hello").unwrap();
        let mut cursor = Cursor::new(buf);
        assert_eq!(read_vlq_string(&mut cursor).unwrap(), "hello");
    }

    #[test]
    fn empty_string() {
        let mut buf = Vec::new();
        write_vlq_string(&mut buf, "").unwrap();
        assert_eq!(buf.len(), 1); // just the 0 length byte
        let mut cursor = Cursor::new(buf);
        assert_eq!(read_vlq_string(&mut cursor).unwrap(), "");
    }

    #[test]
    fn unicode_string() {
        let s = "café ☕";
        let mut buf = Vec::new();
        write_vlq_string(&mut buf, s).unwrap();
        let mut cursor = Cursor::new(buf);
        assert_eq!(read_vlq_string(&mut cursor).unwrap(), s);
    }

    #[test]
    fn sequential_vlq_reads() {
        let mut buf = Vec::new();
        write_vlq(&mut buf, 42).unwrap();
        write_vlq(&mut buf, 300).unwrap();
        write_vlq(&mut buf, 0).unwrap();

        let mut cursor = Cursor::new(buf);
        assert_eq!(read_vlq(&mut cursor).unwrap(), 42);
        assert_eq!(read_vlq(&mut cursor).unwrap(), 300);
        assert_eq!(read_vlq(&mut cursor).unwrap(), 0);
    }
}
