//! VCF file reading: BGZF decompression, header parsing, record iteration.

pub mod assembly_inference;
pub mod bgzf;
pub mod bgzf_stream;
pub mod header;
pub mod info;
pub mod reader;
pub mod record;
pub mod sample;
pub mod streaming_reader;
