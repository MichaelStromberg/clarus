//! RNA and protein sequence dictionaries for transcript construction.

use std::collections::HashMap;
use std::io::Read;

use crate::error::Error;
use crate::fasta;

/// RNA sequence dictionary indexed by transcript ID.
pub struct RnaSequences {
    sequences: HashMap<String, Vec<u8>>,
}

impl RnaSequences {
    /// Build from a gzip-compressed RNA FASTA file.
    pub fn from_gz<R: Read>(reader: R) -> Result<Self, Error> {
        let entries = fasta::parse_fasta_gz(reader)?;
        let mut sequences = HashMap::with_capacity(entries.len());
        for (id, seq) in entries {
            if sequences.contains_key(&id) {
                return Err(Error::Validation(format!(
                    "duplicate transcript ID in RNA FASTA: {id}"
                )));
            }
            sequences.insert(id, seq);
        }
        Ok(Self { sequences })
    }

    /// Get a cDNA sequence by transcript ID.
    #[must_use]
    pub fn get(&self, transcript_id: &str) -> Option<&[u8]> {
        self.sequences.get(transcript_id).map(|v| v.as_slice())
    }

    /// Get a mutable reference to the inner map for mitochondrial sequence insertion.
    pub fn inner_mut(&mut self) -> &mut HashMap<String, Vec<u8>> {
        &mut self.sequences
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }
}

/// Protein sequence dictionary with transcript cross-reference.
pub struct ProteinSequences {
    sequences: HashMap<String, Vec<u8>>,
    transcript_to_protein: HashMap<String, String>,
}

impl ProteinSequences {
    /// Build from a gzip-compressed protein FASTA file.
    ///
    /// Parses headers for ` transcript:` marker to build cross-reference.
    /// Appends `*` stop codon to all protein sequences.
    pub fn from_gz<R: Read>(reader: R) -> Result<Self, Error> {
        let entries = fasta::parse_fasta_gz_with_headers(reader)?;
        let mut sequences = HashMap::with_capacity(entries.len());
        let mut transcript_to_protein = HashMap::new();

        for (primary_id, header, mut seq) in entries {
            // Append stop codon
            seq.push(b'*');

            if sequences.contains_key(&primary_id) {
                return Err(Error::Validation(format!(
                    "duplicate protein ID in protein FASTA: {primary_id}"
                )));
            }

            // Parse transcript cross-reference from header
            if let Some(pos) = header.find(" transcript:") {
                let rest = &header[pos + 12..]; // skip " transcript:"
                let tid = rest.split_whitespace().next().unwrap_or("");
                if !tid.is_empty() {
                    if transcript_to_protein.contains_key(tid) {
                        return Err(Error::Validation(format!(
                            "duplicate transcript ID in protein FASTA cross-reference: {tid}"
                        )));
                    }
                    transcript_to_protein.insert(tid.to_string(), primary_id.clone());
                }
            }

            sequences.insert(primary_id, seq);
        }

        Ok(Self {
            sequences,
            transcript_to_protein,
        })
    }

    /// Get a protein sequence by protein ID.
    #[must_use]
    pub fn get_by_protein_id(&self, protein_id: &str) -> Option<&[u8]> {
        self.sequences.get(protein_id).map(|v| v.as_slice())
    }

    /// Get the protein ID for a transcript ID.
    #[must_use]
    pub fn protein_id_for_transcript(&self, transcript_id: &str) -> Option<&str> {
        self.transcript_to_protein
            .get(transcript_id)
            .map(|s| s.as_str())
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    fn make_gz(content: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(content).unwrap();
        encoder.finish().unwrap()
    }

    #[test]
    fn rna_indexing() {
        let fasta = b">NM_001005484.2 Homo sapiens\nACGT\nTTTT\n>NR_046018.2 something\nAAAA\n";
        let gz = make_gz(fasta);
        let rna = RnaSequences::from_gz(std::io::Cursor::new(gz)).unwrap();
        assert_eq!(rna.len(), 2);
        assert_eq!(rna.get("NM_001005484.2"), Some(b"ACGTTTTT".as_slice()));
        assert_eq!(rna.get("NR_046018.2"), Some(b"AAAA".as_slice()));
        assert!(rna.get("MISSING").is_none());
    }

    #[test]
    fn rna_duplicate_error() {
        let fasta = b">NM_001.1 first\nACGT\n>NM_001.1 second\nTTTT\n";
        let gz = make_gz(fasta);
        assert!(RnaSequences::from_gz(std::io::Cursor::new(gz)).is_err());
    }

    #[test]
    fn protein_indexing_with_cross_reference() {
        let fasta = b">NP_001005484.1 protein [Homo sapiens] transcript:NM_001005484.2\nMPQIQK\n";
        let gz = make_gz(fasta);
        let prot = ProteinSequences::from_gz(std::io::Cursor::new(gz)).unwrap();
        assert_eq!(prot.len(), 1);
        // Stop codon appended
        assert_eq!(
            prot.get_by_protein_id("NP_001005484.1"),
            Some(b"MPQIQK*".as_slice())
        );
        assert_eq!(
            prot.protein_id_for_transcript("NM_001005484.2"),
            Some("NP_001005484.1")
        );
    }

    #[test]
    fn protein_without_transcript_marker() {
        let fasta = b">NP_999.1 some protein\nMACK\n";
        let gz = make_gz(fasta);
        let prot = ProteinSequences::from_gz(std::io::Cursor::new(gz)).unwrap();
        assert_eq!(
            prot.get_by_protein_id("NP_999.1"),
            Some(b"MACK*".as_slice())
        );
        assert!(prot.protein_id_for_transcript("anything").is_none());
    }

    #[test]
    fn protein_duplicate_error() {
        let fasta = b">NP_001.1 first\nMACK\n>NP_001.1 second\nMMMM\n";
        let gz = make_gz(fasta);
        assert!(ProteinSequences::from_gz(std::io::Cursor::new(gz)).is_err());
    }
}
