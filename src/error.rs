//! Error types for the Clarus library.

use thiserror::Error;

/// Errors that can occur during Clarus operations.
#[derive(Debug, Error)]
pub enum Error {
    /// An I/O error occurred.
    #[error(transparent)]
    Io(#[from] std::io::Error),

    /// A parse error occurred while reading input data.
    #[error("{0}")]
    Parse(String),

    /// A validation constraint was violated.
    #[error("{0}")]
    Validation(String),

    /// A file format error was detected.
    #[error("{0}")]
    Format(String),

    /// A transcript that could not be resolved during evaluation.
    /// These are excluded from the cache rather than causing a fatal error.
    #[error("unresolvable transcript: {0}")]
    UnresolvableTranscript(String),
}
