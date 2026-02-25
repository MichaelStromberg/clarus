//! GFF3 data structures for parsed entries and hierarchy records.

use std::collections::HashMap;

use crate::biotype::{BioType, BioTypeCategory};
use crate::error::Error;
use crate::strand::Strand;
use crate::transcript::types::Designations;

/// CIGAR operation type from Gap attributes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CigarOpType {
    Match,
    Insertion,
    Deletion,
}

/// A single CIGAR operation with type and length.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
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
    pub hgnc_id: Option<u32>,
    pub transcript_id: Option<String>,
    pub designations: Designations,
    pub target_id: Option<String>,
    pub target_start: Option<i32>,
    pub target_end: Option<i32>,
    pub cigar_ops: Option<Vec<CigarOp>>,
    pub note: Option<String>,
    pub pubmed_ids: Option<Vec<u32>>,
    pub eco_id: Option<u32>,
    /// Ensembl CDS protein_id attribute.
    pub protein_id: Option<String>,
    /// Ensembl version attribute (integer version number).
    pub version: Option<u16>,
    /// Ensembl biotype attribute from column 9 (distinct from column 3 SO type).
    pub biotype_attr: Option<String>,
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

/// Shared hierarchy data consumed by [`HierarchyData::build`] to produce a [`Gff3Result`].
///
/// Both [`RefSeqHierarchyBuilder`] and [`EnsemblHierarchyBuilder`] delegate their
/// `build()` methods to this struct, which contains the common gene→transcript→exon
/// assembly logic. For Ensembl, `match_children` is empty.
pub(crate) struct HierarchyData {
    pub entries: Vec<Gff3Entry>,
    pub gene_indices: Vec<usize>,
    pub transcript_children: HashMap<usize, Vec<usize>>,
    pub exon_children: HashMap<usize, Vec<usize>>,
    pub cds_children: HashMap<usize, Vec<usize>>,
    pub match_children: HashMap<usize, Vec<usize>>,
    pub regulatory_regions: Vec<Gff3Entry>,
}

impl HierarchyData {
    /// Consume the hierarchy data and produce the final GFF3 result.
    pub fn build(mut self) -> Result<Gff3Result, Error> {
        let mut genes = Vec::with_capacity(self.gene_indices.len());

        let mut entries: Vec<Option<Gff3Entry>> = self.entries.drain(..).map(Some).collect();

        let take = |entries: &mut [Option<Gff3Entry>], i: usize| -> Result<Gff3Entry, Error> {
            entries
                .get_mut(i)
                .and_then(Option::take)
                .ok_or_else(|| Error::Parse(format!("missing GFF3 entry at index {i}")))
        };

        for &gene_idx in &self.gene_indices {
            let gene_entry = take(&mut entries, gene_idx)?;
            let transcript_indices = self
                .transcript_children
                .remove(&gene_idx)
                .unwrap_or_default();

            let mut transcripts = Vec::with_capacity(transcript_indices.len());
            for tx_idx in transcript_indices {
                let tx_entry = take(&mut entries, tx_idx)?;

                let exon_indices = self.exon_children.remove(&tx_idx).unwrap_or_default();
                let mut exons = Vec::with_capacity(exon_indices.len());
                for i in exon_indices {
                    exons.push(take(&mut entries, i)?);
                }

                let cds_indices = self.cds_children.remove(&tx_idx).unwrap_or_default();
                let mut cds_entries = Vec::with_capacity(cds_indices.len());
                for i in cds_indices {
                    cds_entries.push(take(&mut entries, i)?);
                }

                let match_indices = self.match_children.remove(&tx_idx).unwrap_or_default();
                let mut matches = Vec::with_capacity(match_indices.len());
                for i in match_indices {
                    matches.push(MatchRecord {
                        entry: take(&mut entries, i)?,
                        exons: Vec::new(),
                    });
                }

                transcripts.push(TranscriptRecord {
                    entry: tx_entry,
                    exons,
                    cds_entries,
                    matches,
                });
            }

            genes.push(GeneRecord {
                entry: gene_entry,
                transcripts,
            });
        }

        Ok(Gff3Result {
            genes,
            regulatory_regions: self.regulatory_regions,
        })
    }
}

/// Construct a lookup key from an ID string and chromosome index.
pub(crate) fn make_id_key(id: &str, chromosome_index: usize) -> String {
    format!("{id}|{chromosome_index}")
}

/// Result of parsing a single GFF3 line.
pub enum ParsedLine {
    Entry(Box<Gff3Entry>),
    Discarded,
    Comment,
    EndOfSection,
}

/// Parse a single GFF3 line using the given attribute parser.
///
/// Handles comment/directive detection, column splitting, chromosome lookup,
/// biotype parsing, biotype-category discard, start/end parsing, and strand parsing.
/// The `parse_attrs` closure is called with GFF3 column 9 to produce source-specific attributes.
pub fn parse_gff3_line(
    line: &str,
    name_to_index: &HashMap<String, usize>,
    parse_attrs: impl FnOnce(&str) -> Result<Gff3Attributes, Error>,
) -> Result<ParsedLine, Error> {
    if line.starts_with('#') {
        if line == "###" {
            return Ok(ParsedLine::EndOfSection);
        }
        return Ok(ParsedLine::Comment);
    }

    let line = line.trim();
    if line.is_empty() {
        return Ok(ParsedLine::Comment);
    }

    let columns: Vec<&str> = line.split('\t').collect();
    if columns.len() != 9 {
        return Err(Error::Parse(format!(
            "GFF3 line has {} columns, expected 9",
            columns.len()
        )));
    }

    let chr_name = columns[0];
    let Some(&chromosome_index) = name_to_index.get(chr_name) else {
        return Ok(ParsedLine::Discarded);
    };

    let biotype: BioType = columns[2].parse()?;
    let category = biotype.category();
    if matches!(
        category,
        BioTypeCategory::NotUseful | BioTypeCategory::Unsupported
    ) {
        return Ok(ParsedLine::Discarded);
    }

    let start: i32 = columns[3]
        .parse()
        .map_err(|e| Error::Parse(format!("invalid start '{}': {e}", columns[3])))?;
    let end: i32 = columns[4]
        .parse()
        .map_err(|e| Error::Parse(format!("invalid end '{}': {e}", columns[4])))?;

    let strand = Strand::from_gff3(columns[6]);
    let attributes = parse_attrs(columns[8])?;

    if attributes.id.is_empty() {
        return Err(Error::Parse(format!(
            "GFF3 entry missing required ID attribute: {line}"
        )));
    }

    Ok(ParsedLine::Entry(Box::new(Gff3Entry {
        chromosome_index,
        start,
        end,
        biotype,
        strand,
        attributes,
    })))
}
