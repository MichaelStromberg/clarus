//! Transcript data model types for intermediate representation.

use std::sync::Arc;

use crate::biotype::BioType;
use crate::gff3::entry::CigarOp;
use crate::strand::Strand;

/// Type of a transcript region.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TranscriptRegionType {
    Exon,
    Intron,
}

/// A region within a transcript (exon or intron) mapping genomic to cDNA coordinates.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct TranscriptRegion {
    pub region_type: TranscriptRegionType,
    pub id: u16,
    pub genomic_start: i32,
    pub genomic_end: i32,
    pub cdna_start: i32,
    pub cdna_end: i32,
    pub cigar_ops: Option<Vec<CigarOp>>,
}

/// An amino acid edit correction from transcript evaluation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AminoAcidEdit {
    pub position: u32,  // 1-based position in translated protein
    pub amino_acid: u8, // Expected amino acid from reference protein
}

/// Translational slippage record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TranslationalSlip {
    pub position: i32, // 1-based CDS position
    pub length: u8,    // 1 or 2
}

/// Coding region of a transcript.
#[derive(Debug, Clone)]
pub struct CodingRegion {
    pub genomic_start: i32,
    pub genomic_end: i32,
    pub cdna_start: i32,
    pub cdna_end: i32,
    pub protein_id: String,
    pub protein_seq: Vec<u8>,
    pub cds_padding: u8,
    pub cds_offset: u16,
    pub protein_offset: u16,
    pub amino_acid_edits: Option<Vec<AminoAcidEdit>>,
    pub slip: Option<TranslationalSlip>,
}

/// Annotation source.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Source {
    RefSeq,
    Ensembl,
}

/// Designation flags for transcript clinical/reference significance.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[allow(clippy::struct_excessive_bools)]
pub struct Designations {
    pub is_mane_select: bool,
    pub is_mane_plus_clinical: bool,
    pub is_refseq_select: bool,
    pub is_ensembl_canonical: bool,
}

impl std::ops::BitOrAssign for Designations {
    fn bitor_assign(&mut self, rhs: Self) {
        self.is_mane_select |= rhs.is_mane_select;
        self.is_mane_plus_clinical |= rhs.is_mane_plus_clinical;
        self.is_refseq_select |= rhs.is_refseq_select;
        self.is_ensembl_canonical |= rhs.is_ensembl_canonical;
    }
}

/// Gene record in the intermediate representation.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Gene {
    pub chromosome_index: usize,
    pub symbol: String,
    pub ncbi_gene_id: Option<String>,
    pub ensembl_id: Option<String>,
    pub hgnc_id: Option<u32>,
    pub on_reverse_strand: bool,
}

/// Intermediate transcript record.
#[derive(Debug)]
pub struct IntermediateTranscript {
    pub chromosome_index: usize,
    pub start: i32,
    pub end: i32,
    pub id: String,
    pub biotype: BioType,
    pub source: Source,
    pub strand: Strand,
    pub gene: Arc<Gene>,
    pub transcript_regions: Vec<TranscriptRegion>,
    pub coding_region: Option<CodingRegion>,
    pub cdna_seq: Vec<u8>,
    pub designations: Designations,
}
