//! Ensembl GFF3 parent-child hierarchy builder.
//!
//! Ensembl uses colon-delimited IDs (`gene:ENSG...`, `transcript:ENST...`),
//! registers `gene` and `transcript` prefixes, filters transcripts by biotype
//! attribute (rejecting `artifact` and `unconfirmed_transcript`), and never
//! has `cDNA_match` entries.

use std::collections::HashMap;

use crate::biotype::{BioType, BioTypeCategory};
use crate::error::Error;

use super::entry::{Gff3Entry, Gff3Result, HierarchyData, make_id_key};

/// Builds a gene-transcript-exon hierarchy from flat Ensembl GFF3 entries.
pub struct EnsemblHierarchyBuilder {
    /// All entries stored by index.
    entries: Vec<Gff3Entry>,
    /// ID key → entry index in `entries`.
    id_to_index: HashMap<String, usize>,
    /// Gene indices (in order of appearance).
    gene_indices: Vec<usize>,
    /// Transcript children: parent_index → vec of transcript indices.
    transcript_children: HashMap<usize, Vec<usize>>,
    /// Exon children: parent_index → vec of exon indices.
    exon_children: HashMap<usize, Vec<usize>>,
    /// CDS children: parent_index → vec of CDS indices.
    cds_children: HashMap<usize, Vec<usize>>,
    /// Regulatory region entries (collected separately).
    regulatory_regions: Vec<Gff3Entry>,
}

impl Default for EnsemblHierarchyBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl EnsemblHierarchyBuilder {
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
            id_to_index: HashMap::new(),
            gene_indices: Vec::new(),
            transcript_children: HashMap::new(),
            exon_children: HashMap::new(),
            cds_children: HashMap::new(),
            regulatory_regions: Vec::new(),
        }
    }

    /// Add a parsed GFF3 entry to the hierarchy.
    pub fn add(&mut self, entry: Gff3Entry) -> Result<(), Error> {
        let category = entry.biotype.category();

        // Regulatory regions go to a separate list
        if category == BioTypeCategory::Regulatory {
            self.regulatory_regions.push(entry);
            return Ok(());
        }

        let entry_index = self.entries.len();

        // Extract the stable ID from colon-delimited format (e.g., "gene:ENSG..." → "ENSG...")
        let stable_id = ensembl_stable_id(&entry.attributes.id);
        let id_key = make_id_key(stable_id, entry.chromosome_index);

        // Determine entity type from ID prefix
        let prefix = ensembl_id_prefix(&entry.attributes.id);

        match prefix {
            "gene" => {
                // Register gene in lookup
                if self.id_to_index.contains_key(&id_key) {
                    return Err(Error::Parse(format!("duplicate GFF3 ID key: '{id_key}'")));
                }
                self.id_to_index.insert(id_key, entry_index);

                if category == BioTypeCategory::Gene {
                    self.entries.push(entry);
                    self.gene_indices.push(entry_index);
                } else {
                    return Err(Error::Parse(format!(
                        "gene-prefixed entry has unexpected biotype: {:?}",
                        entry.biotype
                    )));
                }
            }
            "transcript" => {
                // Filter out artifact and unconfirmed transcripts
                if let Some(ref bt) = entry.attributes.biotype_attr
                    && (bt == "artifact" || bt == "unconfirmed_transcript")
                {
                    return Ok(());
                }

                // Register transcript in lookup
                if self.id_to_index.contains_key(&id_key) {
                    return Err(Error::Parse(format!("duplicate GFF3 ID key: '{id_key}'")));
                }
                self.id_to_index.insert(id_key, entry_index);

                // Resolve parent gene
                let parent_id = entry.attributes.parent_id.as_ref().ok_or_else(|| {
                    Error::Parse(format!(
                        "transcript entry '{}' has no Parent attribute",
                        entry.attributes.id
                    ))
                })?;
                let parent_stable_id = ensembl_stable_id(parent_id);
                let parent_key = make_id_key(parent_stable_id, entry.chromosome_index);
                let &parent_index = self.id_to_index.get(&parent_key).ok_or_else(|| {
                    Error::Parse(format!(
                        "GFF3 parent not found: '{parent_key}' (child: '{}')",
                        entry.attributes.id
                    ))
                })?;

                self.entries.push(entry);
                self.transcript_children
                    .entry(parent_index)
                    .or_default()
                    .push(entry_index);
            }
            _ => {
                // Exon or CDS — needs a parent
                let parent_id = entry.attributes.parent_id.as_ref().ok_or_else(|| {
                    Error::Parse(format!(
                        "GFF3 entry '{}' has no Parent attribute",
                        entry.attributes.id
                    ))
                })?;
                let parent_stable_id = ensembl_stable_id(parent_id);
                let parent_key = make_id_key(parent_stable_id, entry.chromosome_index);
                let Some(&parent_index) = self.id_to_index.get(&parent_key) else {
                    // Parent transcript may have been filtered (artifact/unconfirmed).
                    // Silently skip children of filtered transcripts.
                    return Ok(());
                };

                match entry.biotype {
                    BioType::Cds => {
                        self.entries.push(entry);
                        self.cds_children
                            .entry(parent_index)
                            .or_default()
                            .push(entry_index);
                    }
                    BioType::Exon => {
                        self.entries.push(entry);
                        self.exon_children
                            .entry(parent_index)
                            .or_default()
                            .push(entry_index);
                    }
                    _ => {
                        return Err(Error::Parse(format!(
                            "unexpected Ensembl GFF3 entry type: {:?} (ID: '{}')",
                            entry.biotype, entry.attributes.id
                        )));
                    }
                }
            }
        }

        Ok(())
    }

    /// Consume the builder and produce the final GFF3 result.
    pub fn build(self) -> Result<Gff3Result, Error> {
        HierarchyData {
            entries: self.entries,
            gene_indices: self.gene_indices,
            transcript_children: self.transcript_children,
            exon_children: self.exon_children,
            cds_children: self.cds_children,
            match_children: HashMap::new(),
            regulatory_regions: self.regulatory_regions,
        }
        .build()
    }
}

/// Extract the entity type prefix from an Ensembl ID.
/// E.g., "gene:ENSG00000141510" → "gene", "ENSE00001146308" → "" (no prefix).
fn ensembl_id_prefix(id: &str) -> &str {
    match id.find(':') {
        Some(pos) => &id[..pos],
        None => "",
    }
}

/// Extract the stable ID from an Ensembl colon-delimited ID.
/// E.g., "gene:ENSG00000141510" → "ENSG00000141510", "ENSE00001146308" → "ENSE00001146308".
fn ensembl_stable_id(id: &str) -> &str {
    match id.find(':') {
        Some(pos) => &id[pos + 1..],
        None => id,
    }
}

#[cfg(test)]
mod tests {
    use crate::biotype::BioType;
    use crate::strand::Strand;

    use super::super::entry::Gff3Attributes;
    use super::*;

    fn make_ensembl_entry(
        chr: usize,
        biotype: BioType,
        start: i32,
        end: i32,
        strand: Strand,
        id: &str,
        parent_id: Option<&str>,
        biotype_attr: Option<&str>,
    ) -> Gff3Entry {
        Gff3Entry {
            chromosome_index: chr,
            start,
            end,
            biotype,
            strand,
            attributes: Gff3Attributes {
                id: id.to_string(),
                parent_id: parent_id.map(String::from),
                gene_id: Some("ENSG00000000001".to_string()),
                biotype_attr: biotype_attr.map(String::from),
                ..Default::default()
            },
        }
    }

    #[test]
    fn simple_gene_transcript_exon_cds() {
        let mut builder = EnsemblHierarchyBuilder::new();

        // Gene
        builder
            .add(make_ensembl_entry(
                0,
                BioType::Gene,
                69091,
                70008,
                Strand::Forward,
                "gene:ENSG00000186092",
                None,
                Some("protein_coding"),
            ))
            .unwrap();

        // Transcript
        builder
            .add(make_ensembl_entry(
                0,
                BioType::MRna,
                69091,
                70008,
                Strand::Forward,
                "transcript:ENST00000335137",
                Some("gene:ENSG00000186092"),
                Some("protein_coding"),
            ))
            .unwrap();

        // Exon
        builder
            .add(make_ensembl_entry(
                0,
                BioType::Exon,
                69091,
                70008,
                Strand::Forward,
                "exon:ENSE00002319515",
                Some("transcript:ENST00000335137"),
                None,
            ))
            .unwrap();

        // CDS
        builder
            .add(make_ensembl_entry(
                0,
                BioType::Cds,
                69091,
                69984,
                Strand::Forward,
                "CDS:ENSP00000334393",
                Some("transcript:ENST00000335137"),
                None,
            ))
            .unwrap();

        let result = builder.build().unwrap();
        assert_eq!(result.genes.len(), 1);
        assert_eq!(result.genes[0].transcripts.len(), 1);
        assert_eq!(result.genes[0].transcripts[0].exons.len(), 1);
        assert_eq!(result.genes[0].transcripts[0].cds_entries.len(), 1);
        assert!(result.genes[0].transcripts[0].matches.is_empty());
    }

    #[test]
    fn ncrna_gene_accepted() {
        let mut builder = EnsemblHierarchyBuilder::new();

        // ncRNA_gene — should be accepted as a gene-category biotype
        builder
            .add(make_ensembl_entry(
                0,
                BioType::NcRnaGene,
                11869,
                14409,
                Strand::Forward,
                "gene:ENSG00000223972",
                None,
                Some("lncRNA"),
            ))
            .unwrap();

        // Transcript under it
        builder
            .add(make_ensembl_entry(
                0,
                BioType::LncRna,
                11869,
                14409,
                Strand::Forward,
                "transcript:ENST00000456328",
                Some("gene:ENSG00000223972"),
                Some("lncRNA"),
            ))
            .unwrap();

        let result = builder.build().unwrap();
        assert_eq!(result.genes.len(), 1);
        assert_eq!(result.genes[0].entry.biotype, BioType::NcRnaGene);
        assert_eq!(result.genes[0].transcripts.len(), 1);
    }

    #[test]
    fn artifact_transcript_skipped() {
        let mut builder = EnsemblHierarchyBuilder::new();

        builder
            .add(make_ensembl_entry(
                0,
                BioType::Gene,
                100,
                500,
                Strand::Forward,
                "gene:ENSG00000000001",
                None,
                Some("protein_coding"),
            ))
            .unwrap();

        // Artifact transcript — should be silently skipped
        builder
            .add(make_ensembl_entry(
                0,
                BioType::MRna,
                100,
                500,
                Strand::Forward,
                "transcript:ENST00000000001",
                Some("gene:ENSG00000000001"),
                Some("artifact"),
            ))
            .unwrap();

        let result = builder.build().unwrap();
        assert_eq!(result.genes[0].transcripts.len(), 0);
    }

    #[test]
    fn unconfirmed_transcript_skipped() {
        let mut builder = EnsemblHierarchyBuilder::new();

        builder
            .add(make_ensembl_entry(
                0,
                BioType::Gene,
                100,
                500,
                Strand::Forward,
                "gene:ENSG00000000001",
                None,
                Some("protein_coding"),
            ))
            .unwrap();

        builder
            .add(make_ensembl_entry(
                0,
                BioType::MRna,
                100,
                500,
                Strand::Forward,
                "transcript:ENST00000000001",
                Some("gene:ENSG00000000001"),
                Some("unconfirmed_transcript"),
            ))
            .unwrap();

        let result = builder.build().unwrap();
        assert_eq!(result.genes[0].transcripts.len(), 0);
    }

    #[test]
    fn children_of_filtered_transcript_skipped() {
        let mut builder = EnsemblHierarchyBuilder::new();

        builder
            .add(make_ensembl_entry(
                0,
                BioType::Gene,
                100,
                500,
                Strand::Forward,
                "gene:ENSG00000000001",
                None,
                None,
            ))
            .unwrap();

        // Artifact transcript — filtered
        builder
            .add(make_ensembl_entry(
                0,
                BioType::MRna,
                100,
                500,
                Strand::Forward,
                "transcript:ENST00000000001",
                Some("gene:ENSG00000000001"),
                Some("artifact"),
            ))
            .unwrap();

        // Exon of filtered transcript — should be silently skipped (parent not in lookup)
        builder
            .add(make_ensembl_entry(
                0,
                BioType::Exon,
                100,
                300,
                Strand::Forward,
                "exon:ENSE00000000001",
                Some("transcript:ENST00000000001"),
                None,
            ))
            .unwrap();

        let result = builder.build().unwrap();
        assert_eq!(result.genes[0].transcripts.len(), 0);
    }

    #[test]
    fn regulatory_separation() {
        let mut builder = EnsemblHierarchyBuilder::new();

        builder
            .add(make_ensembl_entry(
                0,
                BioType::Gene,
                100,
                200,
                Strand::Forward,
                "gene:ENSG00000000001",
                None,
                None,
            ))
            .unwrap();

        let reg = Gff3Entry {
            chromosome_index: 0,
            start: 300,
            end: 400,
            biotype: BioType::Enhancer,
            strand: Strand::Forward,
            attributes: Gff3Attributes {
                id: "enhancer:1".to_string(),
                ..Default::default()
            },
        };
        builder.add(reg).unwrap();

        let result = builder.build().unwrap();
        assert_eq!(result.genes.len(), 1);
        assert_eq!(result.regulatory_regions.len(), 1);
        assert_eq!(result.regulatory_regions[0].biotype, BioType::Enhancer);
    }

    #[test]
    fn duplicate_gene_id_error() {
        let mut builder = EnsemblHierarchyBuilder::new();
        builder
            .add(make_ensembl_entry(
                0,
                BioType::Gene,
                100,
                200,
                Strand::Forward,
                "gene:ENSG00000000001",
                None,
                None,
            ))
            .unwrap();
        let result = builder.add(make_ensembl_entry(
            0,
            BioType::Gene,
            300,
            400,
            Strand::Forward,
            "gene:ENSG00000000001",
            None,
            None,
        ));
        assert!(result.is_err());
    }

    #[test]
    fn missing_parent_error() {
        let mut builder = EnsemblHierarchyBuilder::new();
        let result = builder.add(make_ensembl_entry(
            0,
            BioType::MRna,
            100,
            200,
            Strand::Forward,
            "transcript:ENST00000000001",
            Some("gene:ENSG_MISSING"),
            None,
        ));
        assert!(result.is_err());
    }

    #[test]
    fn id_prefix_extraction() {
        assert_eq!(ensembl_id_prefix("gene:ENSG00000141510"), "gene");
        assert_eq!(
            ensembl_id_prefix("transcript:ENST00000269305"),
            "transcript"
        );
        assert_eq!(ensembl_id_prefix("CDS:ENSP00000269305"), "CDS");
        assert_eq!(ensembl_id_prefix("ENSE00001146308"), "");
    }

    #[test]
    fn stable_id_extraction() {
        assert_eq!(ensembl_stable_id("gene:ENSG00000141510"), "ENSG00000141510");
        assert_eq!(
            ensembl_stable_id("transcript:ENST00000269305"),
            "ENST00000269305"
        );
        assert_eq!(ensembl_stable_id("ENSE00001146308"), "ENSE00001146308");
    }
}
