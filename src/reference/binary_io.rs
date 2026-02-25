//! Binary I/O extension traits for reading and writing little-endian primitive types.

use std::io::{Read, Write};

use crate::error::Error;

/// Extension trait for writing little-endian binary values.
pub(crate) trait BinaryWrite: Write {
    fn write_u8(&mut self, value: u8) -> Result<(), Error> {
        self.write_all(&[value])?;
        Ok(())
    }

    fn write_u16(&mut self, value: u16) -> Result<(), Error> {
        self.write_all(&value.to_le_bytes())?;
        Ok(())
    }

    fn write_u32(&mut self, value: u32) -> Result<(), Error> {
        self.write_all(&value.to_le_bytes())?;
        Ok(())
    }

    fn write_u64(&mut self, value: u64) -> Result<(), Error> {
        self.write_all(&value.to_le_bytes())?;
        Ok(())
    }

    fn write_prefixed_string(&mut self, s: &str) -> Result<(), Error> {
        let len = s.len();
        if len > 255 {
            return Err(Error::Validation(format!(
                "string too long for u8 prefix: {len} bytes"
            )));
        }
        self.write_all(&[len as u8])?;
        self.write_all(s.as_bytes())?;
        Ok(())
    }

    fn write_u16_prefixed_string(&mut self, s: &str) -> Result<(), Error> {
        let len = s.len();
        if len > u16::MAX as usize {
            return Err(Error::Validation(format!(
                "string too long for u16 prefix: {len} bytes"
            )));
        }
        self.write_all(&(len as u16).to_le_bytes())?;
        self.write_all(s.as_bytes())?;
        Ok(())
    }
}

/// Extension trait for reading little-endian binary values.
pub(crate) trait BinaryRead: Read {
    fn read_u8(&mut self) -> Result<u8, Error> {
        let mut buf = [0u8; 1];
        self.read_exact(&mut buf)?;
        Ok(buf[0])
    }

    fn read_u16(&mut self) -> Result<u16, Error> {
        let mut buf = [0u8; 2];
        self.read_exact(&mut buf)?;
        Ok(u16::from_le_bytes(buf))
    }

    fn read_u32(&mut self) -> Result<u32, Error> {
        let mut buf = [0u8; 4];
        self.read_exact(&mut buf)?;
        Ok(u32::from_le_bytes(buf))
    }

    fn read_u64(&mut self) -> Result<u64, Error> {
        let mut buf = [0u8; 8];
        self.read_exact(&mut buf)?;
        Ok(u64::from_le_bytes(buf))
    }

    fn read_prefixed_string(&mut self) -> Result<String, Error> {
        let len = self.read_u8()? as usize;
        if len == 0 {
            return Ok(String::new());
        }
        let mut buf = vec![0u8; len];
        self.read_exact(&mut buf)?;
        String::from_utf8(buf).map_err(|e| Error::Parse(format!("invalid UTF-8: {e}")))
    }

    #[allow(dead_code)] // Reader counterpart of write_u16_prefixed_string, used by cache reader
    fn read_u16_prefixed_string(&mut self) -> Result<String, Error> {
        let len = self.read_u16()? as usize;
        if len == 0 {
            return Ok(String::new());
        }
        let mut buf = vec![0u8; len];
        self.read_exact(&mut buf)?;
        String::from_utf8(buf).map_err(|e| Error::Parse(format!("invalid UTF-8: {e}")))
    }
}

impl<W: Write + ?Sized> BinaryWrite for W {}
impl<R: Read + ?Sized> BinaryRead for R {}
