//! Clarus: high-performance genomic variant annotation and VUS resolution engine.

// Pedantic lint configuration: allow noisy lints that don't improve code quality
// in this domain. Genomic file formats (GFF3, GenBank, FASTA) use many domain terms
// that aren't Rust identifiers, and every public function returns Result.
#![allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::doc_markdown,
    clippy::must_use_candidate,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    clippy::cast_lossless,
    clippy::cast_possible_wrap,
    clippy::similar_names,
    clippy::module_name_repetitions,
    clippy::implicit_hasher
)]

pub mod error;

pub mod assembly_report;
pub mod bands;
pub mod biotype;
pub mod cache;
pub mod canonical;
pub mod chromosome;
pub mod cli;
pub mod codon;
pub mod config;
pub mod download;
pub mod evaluation;
pub mod fasta;
pub mod genbank;
pub mod genome_assembly;
pub mod gff3;
pub mod hgnc;
pub mod mirna;
pub mod perf;
pub mod pipeline;
pub mod reference;
pub mod sequence;
pub mod strand;
pub mod transcript;
pub mod vlq;
