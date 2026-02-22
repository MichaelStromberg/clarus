//! GFF3 parent-child hierarchy builder.

use std::collections::HashMap;

use crate::biotype::{BioType, BioTypeCategory};
use crate::error::Error;

use super::entry::{GeneRecord, Gff3Entry, Gff3Result, MatchRecord, TranscriptRecord};

/// Builds a gene-transcript-exon hierarchy from flat GFF3 entries.
pub struct HierarchyBuilder {
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
    /// Match entries: parent_index → vec of match indices.
    match_children: HashMap<usize, Vec<usize>>,
    /// Exon key → parent ID string (for cDNA_match association).
    exon_keys: HashMap<String, String>,
    /// Regulatory region entries (collected separately).
    regulatory_regions: Vec<Gff3Entry>,
}

impl Default for HierarchyBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl HierarchyBuilder {
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
            id_to_index: HashMap::new(),
            gene_indices: Vec::new(),
            transcript_children: HashMap::new(),
            exon_children: HashMap::new(),
            cds_children: HashMap::new(),
            match_children: HashMap::new(),
            exon_keys: HashMap::new(),
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
        let id_key = make_id_key(&entry.attributes.id, entry.chromosome_index);

        // Register in ID dictionary if prefix is "gene", "id", or "rna"
        let id_prefix = id_prefix(&entry.attributes.id);
        if matches!(id_prefix, "gene" | "id" | "rna") {
            if self.id_to_index.contains_key(&id_key) {
                return Err(Error::Parse(format!("duplicate GFF3 ID key: '{id_key}'")));
            }
            self.id_to_index.insert(id_key.clone(), entry_index);
        }

        // Gene handling
        if entry.biotype == BioType::Gene || entry.biotype == BioType::Pseudogene {
            self.entries.push(entry);
            self.gene_indices.push(entry_index);
            return Ok(());
        }

        // cDNA_match handling — find parent via exon key
        if entry.biotype == BioType::CdnaMatch {
            let target_id = entry
                .attributes
                .target_id
                .as_deref()
                .ok_or_else(|| Error::Parse("cDNA_match missing Target ID".to_string()))?;
            let exon_key = make_exon_key(entry.chromosome_index, target_id, entry.start);
            let parent_id = self.exon_keys.get(&exon_key).ok_or_else(|| {
                Error::Parse(format!("cDNA_match: no exon key match for '{exon_key}'"))
            })?;
            let parent_id_key = make_id_key(parent_id, entry.chromosome_index);
            let &parent_index = self.id_to_index.get(&parent_id_key).ok_or_else(|| {
                Error::Parse(format!(
                    "cDNA_match: parent ID not found in dictionary: '{parent_id_key}'"
                ))
            })?;
            self.entries.push(entry);
            self.match_children
                .entry(parent_index)
                .or_default()
                .push(entry_index);
            return Ok(());
        }

        // Everything else needs a parent
        let parent_id = entry.attributes.parent_id.as_ref().ok_or_else(|| {
            Error::Parse(format!(
                "GFF3 entry '{}' has no Parent attribute",
                entry.attributes.id
            ))
        })?;

        let parent_id_key = make_id_key(parent_id, entry.chromosome_index);
        let &parent_index = self.id_to_index.get(&parent_id_key).ok_or_else(|| {
            Error::Parse(format!(
                "GFF3 parent not found: '{parent_id_key}' (child: '{}')",
                entry.attributes.id
            ))
        })?;

        match entry.biotype {
            BioType::Cds => {
                self.entries.push(entry);
                self.cds_children
                    .entry(parent_index)
                    .or_default()
                    .push(entry_index);
            }
            BioType::Exon => {
                // Register exon key for cDNA_match association
                if let Some(ref tid) = entry.attributes.transcript_id {
                    let exon_key = make_exon_key(entry.chromosome_index, tid, entry.start);
                    self.exon_keys.insert(exon_key, parent_id.to_owned());
                }
                self.entries.push(entry);
                self.exon_children
                    .entry(parent_index)
                    .or_default()
                    .push(entry_index);
            }
            _ => {
                // Transcript entry — filter by Name prefix (NM/NR only for RefSeq)
                if let Some(ref name) = entry.attributes.name {
                    let name_prefix = name_prefix(name);
                    if name_prefix == "NM" || name_prefix == "NR" {
                        self.entries.push(entry);
                        self.transcript_children
                            .entry(parent_index)
                            .or_default()
                            .push(entry_index);
                    }
                    // Otherwise silently skip (XM, XR, etc.)
                }
                // Entries with no Name are silently skipped
            }
        }

        Ok(())
    }

    /// Consume the builder and produce the final GFF3 result.
    pub fn build(mut self) -> Gff3Result {
        let mut genes = Vec::with_capacity(self.gene_indices.len());

        // We need to take entries out by index. Build a vec of Option to allow taking.
        let mut entries: Vec<Option<Gff3Entry>> = self.entries.drain(..).map(Some).collect();

        for &gene_idx in &self.gene_indices {
            let gene_entry = entries[gene_idx].take().unwrap();
            let transcript_indices = self
                .transcript_children
                .remove(&gene_idx)
                .unwrap_or_default();

            let mut transcripts = Vec::with_capacity(transcript_indices.len());
            for tx_idx in transcript_indices {
                let tx_entry = entries[tx_idx].take().unwrap();

                let exon_indices = self.exon_children.remove(&tx_idx).unwrap_or_default();
                let exons: Vec<Gff3Entry> = exon_indices
                    .into_iter()
                    .map(|i| entries[i].take().unwrap())
                    .collect();

                let cds_indices = self.cds_children.remove(&tx_idx).unwrap_or_default();
                let cds_entries: Vec<Gff3Entry> = cds_indices
                    .into_iter()
                    .map(|i| entries[i].take().unwrap())
                    .collect();

                let match_indices = self.match_children.remove(&tx_idx).unwrap_or_default();
                let matches: Vec<MatchRecord> = match_indices
                    .into_iter()
                    .map(|i| {
                        let match_entry = entries[i].take().unwrap();
                        MatchRecord {
                            entry: match_entry,
                            exons: Vec::new(), // Populated during construction
                        }
                    })
                    .collect();

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

        Gff3Result {
            genes,
            regulatory_regions: self.regulatory_regions,
        }
    }
}

/// Construct an ID key: "id|chromosome_index"
fn make_id_key(id: &str, chromosome_index: usize) -> String {
    format!("{id}|{chromosome_index}")
}

/// Construct an exon key: "chr_index|transcript_id|start"
fn make_exon_key(chromosome_index: usize, transcript_id: &str, start: i32) -> String {
    format!("{chromosome_index}|{transcript_id}|{start}")
}

/// Extract prefix before first '-' from an ID string.
fn id_prefix(id: &str) -> &str {
    id.split('-').next().unwrap_or(id)
}

/// Extract prefix before first '_' from a Name string.
fn name_prefix(name: &str) -> &str {
    name.split('_').next().unwrap_or(name)
}

#[cfg(test)]
mod tests {
    use crate::biotype::BioType;
    use crate::strand::Strand;

    use super::super::entry::Gff3Attributes;
    use super::*;

    fn make_entry(
        chr: usize,
        biotype: BioType,
        start: i32,
        end: i32,
        strand: Strand,
        id: &str,
        parent_id: Option<&str>,
        name: Option<&str>,
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
                name: name.map(String::from),
                gene_symbol: Some("TEST".to_string()),
                gene_id: Some("12345".to_string()),
                ..Default::default()
            },
        }
    }

    #[test]
    fn simple_gene_transcript_exon() {
        let mut builder = HierarchyBuilder::new();

        // Gene
        builder
            .add(make_entry(
                0,
                BioType::Gene,
                11874,
                14409,
                Strand::Forward,
                "gene-DDX11L1",
                None,
                None,
            ))
            .unwrap();

        // Transcript (NR prefix accepted)
        builder
            .add(make_entry(
                0,
                BioType::MRna,
                11874,
                14409,
                Strand::Forward,
                "rna-NR_046018.2",
                Some("gene-DDX11L1"),
                Some("NR_046018.2"),
            ))
            .unwrap();

        // Exon 1
        let mut exon1 = make_entry(
            0,
            BioType::Exon,
            11874,
            12227,
            Strand::Forward,
            "exon-NR_046018.2-1",
            Some("rna-NR_046018.2"),
            None,
        );
        exon1.attributes.transcript_id = Some("NR_046018.2".to_string());
        builder.add(exon1).unwrap();

        // Exon 2
        let mut exon2 = make_entry(
            0,
            BioType::Exon,
            12613,
            12721,
            Strand::Forward,
            "exon-NR_046018.2-2",
            Some("rna-NR_046018.2"),
            None,
        );
        exon2.attributes.transcript_id = Some("NR_046018.2".to_string());
        builder.add(exon2).unwrap();

        let result = builder.build();
        assert_eq!(result.genes.len(), 1);
        assert_eq!(result.genes[0].transcripts.len(), 1);
        assert_eq!(result.genes[0].transcripts[0].exons.len(), 2);
    }

    #[test]
    fn regulatory_separation() {
        let mut builder = HierarchyBuilder::new();
        builder
            .add(make_entry(
                0,
                BioType::Gene,
                100,
                200,
                Strand::Forward,
                "gene-X",
                None,
                None,
            ))
            .unwrap();

        // Regulatory region
        let reg = Gff3Entry {
            chromosome_index: 0,
            start: 300,
            end: 400,
            biotype: BioType::Enhancer,
            strand: Strand::Forward,
            attributes: Gff3Attributes {
                id: "enhancer-1".to_string(),
                gene_id: Some("999".to_string()),
                ..Default::default()
            },
        };
        builder.add(reg).unwrap();

        let result = builder.build();
        assert_eq!(result.genes.len(), 1);
        assert_eq!(result.regulatory_regions.len(), 1);
        assert_eq!(result.regulatory_regions[0].biotype, BioType::Enhancer);
    }

    #[test]
    fn xm_transcript_skipped() {
        let mut builder = HierarchyBuilder::new();
        builder
            .add(make_entry(
                0,
                BioType::Gene,
                100,
                200,
                Strand::Forward,
                "gene-X",
                None,
                None,
            ))
            .unwrap();
        // XM prefix — should be silently skipped
        builder
            .add(make_entry(
                0,
                BioType::MRna,
                100,
                200,
                Strand::Forward,
                "rna-XM_001.1",
                Some("gene-X"),
                Some("XM_001.1"),
            ))
            .unwrap();

        let result = builder.build();
        assert_eq!(result.genes[0].transcripts.len(), 0);
    }

    #[test]
    fn duplicate_id_error() {
        let mut builder = HierarchyBuilder::new();
        builder
            .add(make_entry(
                0,
                BioType::Gene,
                100,
                200,
                Strand::Forward,
                "gene-X",
                None,
                None,
            ))
            .unwrap();
        let result = builder.add(make_entry(
            0,
            BioType::Gene,
            300,
            400,
            Strand::Forward,
            "gene-X",
            None,
            None,
        ));
        assert!(result.is_err());
    }

    #[test]
    fn missing_parent_error() {
        let mut builder = HierarchyBuilder::new();
        let result = builder.add(make_entry(
            0,
            BioType::Exon,
            100,
            200,
            Strand::Forward,
            "exon-1",
            Some("rna-MISSING"),
            None,
        ));
        assert!(result.is_err());
    }
}
