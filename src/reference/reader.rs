//! Reader for Clarus reference sequence files.

use std::collections::HashMap;
use std::io::{Read, Seek, SeekFrom};

use crate::bands::Band;
use crate::chromosome::Chromosome;
use crate::error::Error;
use crate::genome_assembly::GenomeAssembly;
use crate::mirna::MirnaRegion;
use crate::reference::binary_io::BinaryRead;
use crate::reference::common_header::{
    REFERENCE_FILE_TYPE, REFERENCE_FORMAT_VERSION, read_common_header,
};

/// Represents a loaded chromosome's data.
pub struct ChromosomeData {
    pub sequence_offset: u32,
    pub sequence: Vec<u8>,
    pub bands: Vec<Band>,
    pub mirna: Vec<MirnaRegion>,
}

impl ChromosomeData {
    /// Extract a substring at 0-based offset with given length.
    /// Returns None if out of bounds.
    #[must_use]
    pub fn substring(&self, offset: i64, length: usize) -> Option<&[u8]> {
        // Safe: sequence_offset is u32, so the cast to i64 is always lossless.
        let adjusted = offset - self.sequence_offset as i64;
        let num_bases = self.sequence.len() as i64;

        if adjusted < 0 || length < 1 || adjusted >= num_bases {
            return None;
        }

        // Safe: the adjusted < 0 guard above ensures this is non-negative.
        let start = adjusted as usize;
        let end = (start + length).min(self.sequence.len());

        Some(&self.sequence[start..end])
    }
}

/// Reader for Clarus reference sequence files.
pub struct ReferenceReader {
    pub assembly: GenomeAssembly,
    pub patch_level: u8,
    pub reference_id: u32,
    pub chromosomes: Vec<Chromosome>,
    pub name_to_index: HashMap<String, usize>,
    offsets: Vec<u64>,
}

impl ReferenceReader {
    /// Read and parse a reference file header.
    pub fn from_reader<R: Read + Seek>(reader: &mut R) -> Result<Self, Error> {
        // Read common header
        let (file_type, format_version) = read_common_header(reader)?;

        if file_type != REFERENCE_FILE_TYPE {
            return Err(Error::Format(format!(
                "unexpected file type: expected {REFERENCE_FILE_TYPE}, got {file_type}"
            )));
        }
        if format_version != REFERENCE_FORMAT_VERSION {
            return Err(Error::Format(format!(
                "unexpected format version: expected {REFERENCE_FORMAT_VERSION}, got {format_version}"
            )));
        }

        // Read file length
        let _file_length = reader.read_u64()?;

        // Read assembly, patch level, and reference ID
        let assembly_byte = reader.read_u8()?;
        let assembly = GenomeAssembly::try_from(assembly_byte)?;
        let patch_level = reader.read_u8()?;
        let reference_id = reader.read_u32()?;

        // Read chromosome count
        let chrom_count = reader.read_u32()?;

        // Read chromosome records
        let mut chromosomes = Vec::with_capacity(chrom_count as usize);
        let mut offsets = Vec::with_capacity(chrom_count as usize);
        let mut name_to_index: HashMap<String, usize> = HashMap::new();

        for i in 0..chrom_count as usize {
            let ucsc_name = reader.read_prefixed_string()?;
            let ensembl_name = reader.read_prefixed_string()?;
            let refseq_accession = reader.read_prefixed_string()?;
            let genbank_accession = reader.read_prefixed_string()?;
            let length = reader.read_u32()?;
            let offset = reader.read_u64()?;

            let ref_index = u16::try_from(i)
                .map_err(|_| Error::Validation(format!("chromosome index {i} exceeds u16::MAX")))?;

            let chr = Chromosome {
                ucsc_name,
                ensembl_name,
                refseq_accession,
                genbank_accession,
                length,
                ref_index,
            };

            for name in chr.names() {
                name_to_index.insert(name.to_string(), i);
            }

            chromosomes.push(chr);
            offsets.push(offset);
        }

        Ok(ReferenceReader {
            assembly,
            patch_level,
            reference_id,
            chromosomes,
            name_to_index,
            offsets,
        })
    }

    /// Load a chromosome's data by index.
    pub fn load_chromosome<R: Read + Seek>(
        &self,
        reader: &mut R,
        index: usize,
    ) -> Result<ChromosomeData, Error> {
        if index >= self.chromosomes.len() {
            return Err(Error::Validation(format!(
                "chromosome index {index} out of range (count: {})",
                self.chromosomes.len()
            )));
        }

        reader.seek(SeekFrom::Start(self.offsets[index]))?;

        let sequence_offset = reader.read_u32()?;
        let compressed_size = reader.read_u32()?;

        let mut compressed = vec![0u8; compressed_size as usize];
        reader.read_exact(&mut compressed)?;

        let sequence = zstd::decode_all(compressed.as_slice())?;

        // Read cytogenetic bands
        let band_count = reader.read_u32()?;
        let mut bands = Vec::with_capacity(band_count as usize);
        for _ in 0..band_count {
            let begin = reader.read_u32()?;
            let end = reader.read_u32()?;
            let name = reader.read_prefixed_string()?;
            bands.push(Band { begin, end, name });
        }

        // Read miRNA regions
        let mirna_count = reader.read_u32()?;
        let mut mirna = Vec::with_capacity(mirna_count as usize);
        for _ in 0..mirna_count {
            let begin = reader.read_u32()?;
            let end = reader.read_u32()?;
            mirna.push(MirnaRegion { begin, end });
        }

        Ok(ChromosomeData {
            sequence_offset,
            sequence,
            bands,
            mirna,
        })
    }

    /// Look up a chromosome index by name.
    #[must_use]
    pub fn get_index(&self, name: &str) -> Option<usize> {
        self.name_to_index.get(name).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::writer::ReferenceWriter;
    use std::io::Cursor;

    fn make_test_file() -> Vec<u8> {
        let chroms = vec![
            Chromosome {
                ucsc_name: "chr1".to_string(),
                ensembl_name: "1".to_string(),
                refseq_accession: "NC_000001.11".to_string(),
                genbank_accession: "CM000663.2".to_string(),
                length: 10,
                ref_index: 0,
            },
            Chromosome {
                ucsc_name: "chr2".to_string(),
                ensembl_name: "2".to_string(),
                refseq_accession: "NC_000002.12".to_string(),
                genbank_accession: "CM000664.2".to_string(),
                length: 8,
                ref_index: 1,
            },
            Chromosome {
                ucsc_name: "chrX".to_string(),
                ensembl_name: "X".to_string(),
                refseq_accession: "NC_000023.11".to_string(),
                genbank_accession: "CM000685.2".to_string(),
                length: 5,
                ref_index: 2,
            },
        ];
        let sequences = vec![
            b"ACGTACGTNN".to_vec(),
            b"TTTTAAAA".to_vec(),
            b"GGGCC".to_vec(),
        ];
        let bands = HashMap::new();
        let mirna = HashMap::new();

        let mut buf = Cursor::new(Vec::new());
        ReferenceWriter::write(
            &mut buf,
            GenomeAssembly::GRCh38,
            14,
            &chroms,
            &sequences,
            &bands,
            &mirna,
        )
        .expect("failed to write test reference file");
        buf.into_inner()
    }

    fn make_test_file_with_bands() -> Vec<u8> {
        let chroms = vec![
            Chromosome {
                ucsc_name: "chr1".to_string(),
                ensembl_name: "1".to_string(),
                refseq_accession: "NC_000001.11".to_string(),
                genbank_accession: "CM000663.2".to_string(),
                length: 10,
                ref_index: 0,
            },
            Chromosome {
                ucsc_name: "chr2".to_string(),
                ensembl_name: "2".to_string(),
                refseq_accession: "NC_000002.12".to_string(),
                genbank_accession: "CM000664.2".to_string(),
                length: 8,
                ref_index: 1,
            },
        ];
        let sequences = vec![b"ACGTACGTNN".to_vec(), b"TTTTAAAA".to_vec()];

        let mut bands = HashMap::new();
        bands.insert(
            "1".to_string(),
            vec![
                Band {
                    begin: 1,
                    end: 2_300_000,
                    name: "p36.33".to_string(),
                },
                Band {
                    begin: 2_300_001,
                    end: 5_300_000,
                    name: "p36.32".to_string(),
                },
                Band {
                    begin: 243_500_001,
                    end: 248_956_422,
                    name: "q44".to_string(),
                },
            ],
        );
        // chr2 has no bands
        let mirna = HashMap::new();

        let mut buf = Cursor::new(Vec::new());
        ReferenceWriter::write(
            &mut buf,
            GenomeAssembly::GRCh38,
            14,
            &chroms,
            &sequences,
            &bands,
            &mirna,
        )
        .expect("failed to write test reference file");
        buf.into_inner()
    }

    #[test]
    fn round_trip_header() {
        let data = make_test_file();
        let mut cursor = Cursor::new(data);
        let reader = ReferenceReader::from_reader(&mut cursor).unwrap();

        assert_eq!(reader.assembly, GenomeAssembly::GRCh38);
        assert_eq!(reader.patch_level, 14);
        assert_eq!(reader.chromosomes.len(), 3);

        assert_eq!(reader.chromosomes[0].ucsc_name, "chr1");
        assert_eq!(reader.chromosomes[0].ensembl_name, "1");
        assert_eq!(reader.chromosomes[0].refseq_accession, "NC_000001.11");
        assert_eq!(reader.chromosomes[0].genbank_accession, "CM000663.2");
        assert_eq!(reader.chromosomes[0].length, 10);

        assert_eq!(reader.chromosomes[1].ucsc_name, "chr2");
        assert_eq!(reader.chromosomes[1].length, 8);

        assert_eq!(reader.chromosomes[2].ucsc_name, "chrX");
        assert_eq!(reader.chromosomes[2].length, 5);
    }

    #[test]
    fn round_trip_sequences() {
        let data = make_test_file();
        let mut cursor = Cursor::new(data);
        let reader = ReferenceReader::from_reader(&mut cursor).unwrap();

        let chr1_data = reader.load_chromosome(&mut cursor, 0).unwrap();
        assert_eq!(chr1_data.sequence, b"ACGTACGTNN");
        assert_eq!(chr1_data.sequence_offset, 0);
        assert!(chr1_data.bands.is_empty());

        let chr2_data = reader.load_chromosome(&mut cursor, 1).unwrap();
        assert_eq!(chr2_data.sequence, b"TTTTAAAA");
        assert!(chr2_data.bands.is_empty());

        let chrx_data = reader.load_chromosome(&mut cursor, 2).unwrap();
        assert_eq!(chrx_data.sequence, b"GGGCC");
        assert!(chrx_data.bands.is_empty());
    }

    #[test]
    fn round_trip_bands() {
        let data = make_test_file_with_bands();
        let mut cursor = Cursor::new(data);
        let reader = ReferenceReader::from_reader(&mut cursor).unwrap();

        // chr1 should have 3 bands
        let chr1_data = reader.load_chromosome(&mut cursor, 0).unwrap();
        assert_eq!(chr1_data.sequence, b"ACGTACGTNN");
        assert_eq!(chr1_data.bands.len(), 3);

        assert_eq!(chr1_data.bands[0].begin, 1);
        assert_eq!(chr1_data.bands[0].end, 2_300_000);
        assert_eq!(chr1_data.bands[0].name, "p36.33");

        assert_eq!(chr1_data.bands[1].begin, 2_300_001);
        assert_eq!(chr1_data.bands[1].end, 5_300_000);
        assert_eq!(chr1_data.bands[1].name, "p36.32");

        assert_eq!(chr1_data.bands[2].begin, 243_500_001);
        assert_eq!(chr1_data.bands[2].end, 248_956_422);
        assert_eq!(chr1_data.bands[2].name, "q44");

        // chr2 should have no bands
        let chr2_data = reader.load_chromosome(&mut cursor, 1).unwrap();
        assert_eq!(chr2_data.sequence, b"TTTTAAAA");
        assert!(chr2_data.bands.is_empty());
    }

    #[test]
    fn name_lookups() {
        let data = make_test_file();
        let mut cursor = Cursor::new(data);
        let reader = ReferenceReader::from_reader(&mut cursor).unwrap();

        assert_eq!(reader.get_index("chr1"), Some(0));
        assert_eq!(reader.get_index("1"), Some(0));
        assert_eq!(reader.get_index("NC_000001.11"), Some(0));
        assert_eq!(reader.get_index("CM000663.2"), Some(0));
        assert_eq!(reader.get_index("chr2"), Some(1));
        assert_eq!(reader.get_index("X"), Some(2));
        assert_eq!(reader.get_index("chrY"), None);
    }

    #[test]
    fn file_length_matches() {
        let data = make_test_file();
        let file_length = u64::from_le_bytes(data[12..20].try_into().unwrap());
        assert_eq!(file_length, data.len() as u64);
    }

    #[test]
    fn file_length_matches_with_bands() {
        let data = make_test_file_with_bands();
        let file_length = u64::from_le_bytes(data[12..20].try_into().unwrap());
        assert_eq!(file_length, data.len() as u64);
    }

    #[test]
    fn substring_extraction() {
        let data = ChromosomeData {
            sequence_offset: 0,
            sequence: b"ACGTACGTNN".to_vec(),
            bands: Vec::new(),
            mirna: Vec::new(),
        };

        // Normal extraction
        assert_eq!(data.substring(0, 4), Some(b"ACGT".as_slice()));
        assert_eq!(data.substring(4, 4), Some(b"ACGT".as_slice()));
        assert_eq!(data.substring(8, 2), Some(b"NN".as_slice()));

        // Clipping at end
        assert_eq!(data.substring(8, 5), Some(b"NN".as_slice()));

        // Out of bounds
        assert_eq!(data.substring(-1, 4), None);
        assert_eq!(data.substring(10, 1), None);
        assert_eq!(data.substring(0, 0), None);
    }
}
