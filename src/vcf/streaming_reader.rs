//! Streaming VCF reader: decompress in bounded batches, parse records on demand.
//!
//! Uses `BgzfStreamReader` for bounded parallel decompression instead of
//! loading the entire decompressed VCF into memory.

use std::collections::HashSet;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;

use crate::error::Error;
use crate::genome_assembly::GenomeAssembly;
use crate::vcf::assembly_inference;
use crate::vcf::bgzf_stream::BgzfStreamReader;
use crate::vcf::header::{self, VcfHeader};
use crate::vcf::record::{RecordParser, VcfRecord};

/// A streaming VCF reader: header is parsed upfront, records are yielded on demand.
pub struct VcfStream {
    pub header: VcfHeader,
    pub inferred_assembly: GenomeAssembly,
    reader: BufReader<BgzfStreamReader>,
    line_buf: String,
    parser: RecordParser,
    num_samples: usize,
    current_chromosome: Option<String>,
    seen_chromosomes: HashSet<String>,
}

impl VcfStream {
    /// Open a VCF file: set up streaming decompression and parse the header.
    /// Records are NOT parsed yet — call `.next_record()` or iterate.
    pub fn open(path: &Path, expected_assembly: GenomeAssembly) -> Result<Self, Error> {
        let bgzf_reader = BgzfStreamReader::new(path)?;
        let mut reader = BufReader::new(bgzf_reader);

        let header = header::parse_header_from_reader(&mut reader)?;

        let inferred_assembly = assembly_inference::infer_assembly(&header.contigs);

        if inferred_assembly != GenomeAssembly::Unknown && inferred_assembly != expected_assembly {
            return Err(Error::AssemblyMismatch(format!(
                "VCF contigs indicate {inferred_assembly} but expected {expected_assembly}"
            )));
        }

        let num_samples = header.sample_names.len();

        Ok(Self {
            header,
            inferred_assembly,
            reader,
            line_buf: String::new(),
            parser: RecordParser::new(),
            num_samples,
            current_chromosome: None,
            seen_chromosomes: HashSet::new(),
        })
    }

    /// Yield the next non-reference, informative record, or `None` at EOF.
    /// Enforces chromosome contiguity.
    pub fn next_record(&mut self) -> Result<Option<VcfRecord>, Error> {
        loop {
            self.line_buf.clear();
            let bytes_read = self.reader.read_line(&mut self.line_buf)?;

            if bytes_read == 0 {
                return Ok(None);
            }

            let line = self.line_buf.trim_end_matches(['\n', '\r']);

            if line.is_empty() {
                continue;
            }

            let record = self.parser.parse(line, self.num_samples)?;

            if record.is_reference || record.alt_alleles.is_empty() {
                continue;
            }

            // Enforce chromosome contiguity
            match &self.current_chromosome {
                Some(current) if *current == record.chromosome => {}
                _ => {
                    if self.seen_chromosomes.contains(&record.chromosome) {
                        return Err(Error::VcfFormat(format!(
                            "chromosome {} reappears after other chromosomes (VCF must have contiguous chromosome blocks)",
                            record.chromosome
                        )));
                    }
                    if let Some(prev) = self.current_chromosome.take() {
                        self.seen_chromosomes.insert(prev);
                    }
                    self.current_chromosome = Some(record.chromosome.clone());
                }
            }

            return Ok(Some(record));
        }
    }
}
