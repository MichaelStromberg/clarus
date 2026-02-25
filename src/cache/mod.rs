//! Binary cache file output: transcript cache and index files.

pub mod bin;
pub mod index;
pub mod serialize;
pub mod writer;

/// File type identifier for transcript cache files.
pub const CACHE_FILE_TYPE: u16 = 2;

/// Format version for transcript cache files.
pub const CACHE_FORMAT_VERSION: u16 = 1;

/// File type identifier for cache index files.
pub const INDEX_FILE_TYPE: u16 = 3;

/// Format version for cache index files.
pub const INDEX_FORMAT_VERSION: u16 = 1;

/// Bit shift for spatial binning (2^20 ≈ 1 megabase bins).
pub const BIN_SHIFT: u32 = 20;

/// Index entry: `(ref_index, byte_offset_in_cache_file)`.
pub type IndexEntry = (u16, u64);
