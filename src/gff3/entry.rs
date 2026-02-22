//! GFF3 data structures for parsed entries and hierarchy records.

use crate::biotype::BioType;
use crate::strand::Strand;

/// CIGAR operation type from Gap attributes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOpType {
    Match,
    Insertion,
    Deletion,
}

/// A single CIGAR operation with type and length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarOp {
    pub op_type: CigarOpType,
    pub length: u32,
}

/// Parsed attributes from GFF3 column 9.
#[derive(Debug, Clone, Default)]
pub struct Gff3Attributes {
    pub id: String,
    pub parent_id: Option<String>,
    pub name: Option<String>,
    pub gene_symbol: Option<String>,
    pub gene_id: Option<String>,
    pub hgnc_id: Option<i32>,
    pub transcript_id: Option<String>,
    pub is_refseq_select: bool,
    pub is_mane_select: bool,
    pub target_id: Option<String>,
    pub target_start: Option<i32>,
    pub target_end: Option<i32>,
    pub cigar_ops: Option<Vec<CigarOp>>,
    pub note: Option<String>,
    pub pubmed_ids: Option<Vec<i32>>,
    pub eco_id: Option<i32>,
}

/// A single parsed GFF3 entry (one line).
#[derive(Debug, Clone)]
pub struct Gff3Entry {
    pub chromosome_index: usize,
    pub start: i32,
    pub end: i32,
    pub biotype: BioType,
    pub strand: Strand,
    pub attributes: Gff3Attributes,
}

/// A gene record with its child transcripts.
#[derive(Debug)]
pub struct GeneRecord {
    pub entry: Gff3Entry,
    pub transcripts: Vec<TranscriptRecord>,
}

/// A transcript record with its child exons, CDS entries, and cDNA_match records.
#[derive(Debug)]
pub struct TranscriptRecord {
    pub entry: Gff3Entry,
    pub exons: Vec<Gff3Entry>,
    pub cds_entries: Vec<Gff3Entry>,
    pub matches: Vec<MatchRecord>,
}

/// A cDNA_match record with its associated exon entries.
#[derive(Debug)]
pub struct MatchRecord {
    pub entry: Gff3Entry,
    pub exons: Vec<Gff3Entry>,
}

/// Result of parsing an entire GFF3 file.
#[derive(Debug)]
pub struct Gff3Result {
    pub genes: Vec<GeneRecord>,
    pub regulatory_regions: Vec<Gff3Entry>,
}
