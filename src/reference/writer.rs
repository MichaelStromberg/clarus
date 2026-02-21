//! Writer for Clarus reference sequence files.

use std::collections::HashMap;
use std::io::{Seek, Write};

use crate::bands::Band;
use crate::chromosome::Chromosome;
use crate::error::Error;
use crate::genome_assembly::{GenomeAssembly, compute_reference_id};
use crate::mirna::MirnaRegion;
use crate::reference::binary_io::BinaryWrite;
use crate::reference::common_header::{
    REFERENCE_FILE_TYPE, REFERENCE_FORMAT_VERSION, write_common_header,
};

/// Writes a Clarus reference sequence file.
pub struct ReferenceWriter;

impl ReferenceWriter {
    /// Writes a complete reference file.
    ///
    /// `chromosomes` and `sequences` must be the same length, ordered by ref_index.
    /// Each sequence is zstd-compressed at level 21.
    /// `bands` maps ensembl chromosome name to cytogenetic bands.
    pub fn write<W: Write + Seek>(
        writer: &mut W,
        assembly: GenomeAssembly,
        patch_level: u8,
        chromosomes: &[Chromosome],
        sequences: &[Vec<u8>],
        bands: &HashMap<String, Vec<Band>>,
        mirna: &HashMap<String, Vec<MirnaRegion>>,
    ) -> Result<(), Error> {
        if chromosomes.len() != sequences.len() {
            return Err(Error::Validation(format!(
                "chromosome count ({}) != sequence count ({})",
                chromosomes.len(),
                sequences.len()
            )));
        }

        let compressed = Self::compress_sequences(sequences)?;
        let (offsets, file_length) =
            Self::calculate_offsets(chromosomes, &compressed, bands, mirna);

        Self::write_header(
            writer,
            assembly,
            patch_level,
            file_length,
            chromosomes,
            &offsets,
        )?;
        Self::write_data_blocks(writer, chromosomes, &compressed, bands, mirna)?;
        Self::verify_position(writer, file_length)?;

        Ok(())
    }

    fn compress_sequences(sequences: &[Vec<u8>]) -> Result<Vec<Vec<u8>>, Error> {
        sequences
            .iter()
            .map(|seq| Ok(zstd::encode_all(seq.as_slice(), 21)?))
            .collect()
    }

    fn calculate_offsets(
        chromosomes: &[Chromosome],
        compressed: &[Vec<u8>],
        bands: &HashMap<String, Vec<Band>>,
        mirna: &HashMap<String, Vec<MirnaRegion>>,
    ) -> (Vec<u64>, u64) {
        let header_size = Self::calculate_header_size(chromosomes);
        let mut offsets = Vec::with_capacity(chromosomes.len());
        let mut current_offset = header_size as u64;

        for (i, comp) in compressed.iter().enumerate() {
            offsets.push(current_offset);
            let band_size =
                Self::calculate_band_block_size(bands.get(&chromosomes[i].ensembl_name));
            let mirna_size =
                Self::calculate_mirna_block_size(mirna.get(&chromosomes[i].ensembl_name));
            let block_size = 4 + 4 + comp.len() as u64 + band_size as u64 + mirna_size as u64;
            current_offset += block_size;
        }

        (offsets, current_offset)
    }

    fn write_header<W: Write>(
        writer: &mut W,
        assembly: GenomeAssembly,
        patch_level: u8,
        file_length: u64,
        chromosomes: &[Chromosome],
        offsets: &[u64],
    ) -> Result<(), Error> {
        write_common_header(writer, REFERENCE_FILE_TYPE, REFERENCE_FORMAT_VERSION)?;
        writer.write_u64(file_length)?;
        writer.write_u8(assembly.to_byte())?;
        writer.write_u8(patch_level)?;
        writer.write_u32(compute_reference_id(assembly, patch_level))?;
        writer.write_u32(
            u32::try_from(chromosomes.len())
                .map_err(|_| Error::Validation("chromosome count exceeds u32::MAX".into()))?,
        )?;

        for (i, chr) in chromosomes.iter().enumerate() {
            writer.write_prefixed_string(&chr.ucsc_name)?;
            writer.write_prefixed_string(&chr.ensembl_name)?;
            writer.write_prefixed_string(&chr.refseq_accession)?;
            writer.write_prefixed_string(&chr.genbank_accession)?;
            writer.write_u32(chr.length)?;
            writer.write_u64(offsets[i])?;
        }

        Ok(())
    }

    fn write_data_blocks<W: Write>(
        writer: &mut W,
        chromosomes: &[Chromosome],
        compressed: &[Vec<u8>],
        bands: &HashMap<String, Vec<Band>>,
        mirna: &HashMap<String, Vec<MirnaRegion>>,
    ) -> Result<(), Error> {
        for (i, comp) in compressed.iter().enumerate() {
            // Sequence offset (always 0 for complete references)
            writer.write_u32(0)?;
            // Compressed buffer
            writer.write_u32(u32::try_from(comp.len()).map_err(|_| {
                Error::Validation("compressed buffer size exceeds u32::MAX".into())
            })?)?;
            writer.write_all(comp)?;

            Self::write_bands(writer, bands.get(&chromosomes[i].ensembl_name))?;
            Self::write_mirna_regions(writer, mirna.get(&chromosomes[i].ensembl_name))?;
        }
        Ok(())
    }

    fn write_bands<W: Write>(writer: &mut W, bands: Option<&Vec<Band>>) -> Result<(), Error> {
        if let Some(band_list) = bands {
            writer.write_u32(
                u32::try_from(band_list.len())
                    .map_err(|_| Error::Validation("band count exceeds u32::MAX".into()))?,
            )?;
            for band in band_list {
                writer.write_u32(band.begin)?;
                writer.write_u32(band.end)?;
                writer.write_prefixed_string(&band.name)?;
            }
        } else {
            writer.write_u32(0)?;
        }
        Ok(())
    }

    fn write_mirna_regions<W: Write>(
        writer: &mut W,
        mirna: Option<&Vec<MirnaRegion>>,
    ) -> Result<(), Error> {
        if let Some(list) = mirna {
            writer.write_u32(
                u32::try_from(list.len())
                    .map_err(|_| Error::Validation("miRNA count exceeds u32::MAX".into()))?,
            )?;
            for region in list {
                writer.write_u32(region.begin)?;
                writer.write_u32(region.end)?;
            }
        } else {
            writer.write_u32(0)?;
        }
        Ok(())
    }

    fn verify_position<W: Write + Seek>(writer: &mut W, expected: u64) -> Result<(), Error> {
        let actual = writer.stream_position()?;
        if actual != expected {
            return Err(Error::Validation(format!(
                "file length mismatch: expected {expected}, actual {actual}"
            )));
        }
        Ok(())
    }

    /// Calculate the total header size in bytes.
    fn calculate_header_size(chromosomes: &[Chromosome]) -> usize {
        // Common header: signature(8) + file_type(2) + format_version(2) + file_length(8)
        let mut size: usize = 8 + 2 + 2 + 8;
        // Reference header: assembly(1) + patch_level(1) + reference_id(4) + chrom_count(4)
        size += 1 + 1 + 4 + 4;
        // Chromosome records
        for chr in chromosomes {
            size += prefixed_string_size(&chr.ucsc_name);
            size += prefixed_string_size(&chr.ensembl_name);
            size += prefixed_string_size(&chr.refseq_accession);
            size += prefixed_string_size(&chr.genbank_accession);
            size += 4; // length
            size += 8; // offset
        }
        size
    }

    /// Calculate the on-disk size of a band block (count + band records).
    fn calculate_band_block_size(bands: Option<&Vec<Band>>) -> usize {
        let mut size = 4; // band_count (u32)
        if let Some(band_list) = bands {
            for band in band_list {
                size += 4; // begin (u32)
                size += 4; // end (u32)
                size += prefixed_string_size(&band.name);
            }
        }
        size
    }

    /// Calculate the on-disk size of a miRNA block (count + region records).
    fn calculate_mirna_block_size(regions: Option<&Vec<MirnaRegion>>) -> usize {
        let mut size = 4; // count (u32)
        if let Some(list) = regions {
            size += list.len() * 8; // 4 (begin) + 4 (end) per region
        }
        size
    }
}

/// Returns the on-disk size of a uint8-prefix ASCII string.
fn prefixed_string_size(s: &str) -> usize {
    1 + s.len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn header_size_calculation() {
        let chroms = vec![Chromosome {
            ucsc_name: "chr1".to_string(),
            ensembl_name: "1".to_string(),
            refseq_accession: "NC_000001.11".to_string(),
            genbank_accession: "CM000663.2".to_string(),
            length: 248_956_422,
            ref_index: 0,
        }];

        let size = ReferenceWriter::calculate_header_size(&chroms);
        // common header: 20
        // ref header: 10 (assembly(1) + patch_level(1) + reference_id(4) + chrom_count(4))
        // chr record: (1+4) + (1+1) + (1+12) + (1+10) + 4 + 8 = 43
        assert_eq!(size, 20 + 10 + 43);
    }

    #[test]
    fn write_produces_valid_file() {
        let chroms = vec![Chromosome {
            ucsc_name: "chr1".to_string(),
            ensembl_name: "1".to_string(),
            refseq_accession: "NC_000001.11".to_string(),
            genbank_accession: "CM000663.2".to_string(),
            length: 4,
            ref_index: 0,
        }];
        let sequences = vec![b"ACGT".to_vec()];
        let bands = HashMap::new();

        let mut buf = Cursor::new(Vec::new());
        let mirna = HashMap::new();
        ReferenceWriter::write(
            &mut buf,
            GenomeAssembly::GRCh38,
            14,
            &chroms,
            &sequences,
            &bands,
            &mirna,
        )
        .unwrap();

        let data = buf.into_inner();
        // Verify signature
        assert_eq!(&data[0..8], &727_905_341_820_126_089u64.to_le_bytes());
        // Verify file type
        assert_eq!(&data[8..10], &1u16.to_le_bytes());
        // Verify format version
        assert_eq!(&data[10..12], &1u16.to_le_bytes());
        // Verify file length matches
        let file_length = u64::from_le_bytes(data[12..20].try_into().unwrap());
        assert_eq!(file_length, data.len() as u64);
    }

    #[test]
    fn write_with_bands() {
        let chroms = vec![Chromosome {
            ucsc_name: "chr1".to_string(),
            ensembl_name: "1".to_string(),
            refseq_accession: "NC_000001.11".to_string(),
            genbank_accession: "CM000663.2".to_string(),
            length: 4,
            ref_index: 0,
        }];
        let sequences = vec![b"ACGT".to_vec()];
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
            ],
        );

        let mut buf = Cursor::new(Vec::new());
        let mirna = HashMap::new();
        ReferenceWriter::write(
            &mut buf,
            GenomeAssembly::GRCh38,
            14,
            &chroms,
            &sequences,
            &bands,
            &mirna,
        )
        .unwrap();

        let data = buf.into_inner();
        let file_length = u64::from_le_bytes(data[12..20].try_into().unwrap());
        assert_eq!(file_length, data.len() as u64);
    }

    #[test]
    fn mismatched_counts_error() {
        let chroms = vec![Chromosome {
            ucsc_name: "chr1".to_string(),
            ensembl_name: String::new(),
            refseq_accession: String::new(),
            genbank_accession: String::new(),
            length: 4,
            ref_index: 0,
        }];
        let sequences: Vec<Vec<u8>> = vec![];
        let bands = HashMap::new();

        let mut buf = Cursor::new(Vec::new());
        let mirna = HashMap::new();
        let result = ReferenceWriter::write(
            &mut buf,
            GenomeAssembly::GRCh38,
            14,
            &chroms,
            &sequences,
            &bands,
            &mirna,
        );
        assert!(result.is_err());
    }
}
