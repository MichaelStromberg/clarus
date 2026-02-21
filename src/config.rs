use std::path::Path;

use anyhow::{Context, Result, bail};
use serde::Deserialize;

use crate::genome_assembly::GenomeAssembly;

#[derive(Debug, Clone, Deserialize)]
pub struct FileEntry {
    pub url: String,
    pub md5: String,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ReferenceConfig {
    pub genome_assembly: String,
    pub patch_level: u8,
    pub assembly_report: FileEntry,
    pub fasta: FileEntry,
    pub gff: Option<FileEntry>,
    pub ideogram: Option<FileEntry>,
}

impl ReferenceConfig {
    pub fn from_file(path: &Path) -> Result<Self> {
        let content = std::fs::read_to_string(path)
            .with_context(|| format!("failed to read config file: {}", path.display()))?;
        let config: Self = serde_json::from_str(&content)
            .with_context(|| format!("failed to parse config file: {}", path.display()))?;
        config.validate()?;
        Ok(config)
    }

    fn validate(&self) -> Result<()> {
        // Verify assembly name parses
        self.genome_assembly
            .parse::<GenomeAssembly>()
            .with_context(|| format!("invalid genome assembly: '{}'", self.genome_assembly))?;

        // Verify MD5 strings are 32 hex characters
        for (name, entry) in self.file_entries() {
            validate_md5(name, &entry.md5)?;
        }

        Ok(())
    }

    /// Iterate all file entries uniformly as (name, entry) pairs.
    pub fn file_entries(&self) -> impl Iterator<Item = (&str, &FileEntry)> {
        [
            ("assembly_report", &self.assembly_report),
            ("fasta", &self.fasta),
        ]
        .into_iter()
        .chain(self.gff.as_ref().map(|e| ("gff", e)))
        .chain(self.ideogram.as_ref().map(|e| ("ideogram", e)))
    }
}

fn validate_md5(name: &str, md5: &str) -> Result<()> {
    if md5.len() != 32 || !md5.chars().all(|c| c.is_ascii_hexdigit()) {
        bail!("invalid MD5 for '{name}': expected 32 hex characters, got '{md5}'");
    }
    Ok(())
}

fn validate_md5_optional(name: &str, md5: &str) -> Result<()> {
    if md5.is_empty() {
        return Ok(());
    }
    validate_md5(name, md5)
}

// ── Cache Config Types ──────────────────────────────────────────

#[derive(Debug, Deserialize)]
pub struct HgncConfig {
    pub url: String,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct EnsemblConfig {
    pub name: String,
    pub version: String,
    pub description: String,
    pub release_date: String,
    pub gff3: FileEntry,
    pub regulatory_gff: FileEntry,
    pub gene_mysql: FileEntry,
    pub transcript_mysql: FileEntry,
    pub ncrna_fasta: FileEntry,
    pub cdna_fasta: FileEntry,
    pub cds_fasta: FileEntry,
    pub peptide_fasta: FileEntry,
}

impl EnsemblConfig {
    pub fn file_entries(&self) -> impl Iterator<Item = (&str, &FileEntry)> {
        [
            ("gff3", &self.gff3),
            ("regulatory_gff", &self.regulatory_gff),
            ("gene_mysql", &self.gene_mysql),
            ("transcript_mysql", &self.transcript_mysql),
            ("ncrna_fasta", &self.ncrna_fasta),
            ("cdna_fasta", &self.cdna_fasta),
            ("cds_fasta", &self.cds_fasta),
            ("peptide_fasta", &self.peptide_fasta),
        ]
        .into_iter()
    }
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RefSeqCacheConfig {
    pub name: String,
    pub version: String,
    pub description: String,
    pub release_date: String,
    pub genbank: FileEntry,
    pub gff3: FileEntry,
    pub rna_fasta: FileEntry,
    pub protein_fasta: FileEntry,
}

impl RefSeqCacheConfig {
    pub fn file_entries(&self) -> impl Iterator<Item = (&str, &FileEntry)> {
        [
            ("genbank", &self.genbank),
            ("gff3", &self.gff3),
            ("rna_fasta", &self.rna_fasta),
            ("protein_fasta", &self.protein_fasta),
        ]
        .into_iter()
    }
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct CacheConfig {
    pub genome_assembly: String,
    pub ensembl: Option<EnsemblConfig>,
    pub hgnc: HgncConfig,
    pub refseq: Option<RefSeqCacheConfig>,
}

impl CacheConfig {
    pub fn from_file(path: &Path) -> Result<Self> {
        let content = std::fs::read_to_string(path)
            .with_context(|| format!("failed to read config file: {}", path.display()))?;
        let config: Self = serde_json::from_str(&content)
            .with_context(|| format!("failed to parse config file: {}", path.display()))?;
        config.validate()?;
        Ok(config)
    }

    fn validate(&self) -> Result<()> {
        self.genome_assembly
            .parse::<GenomeAssembly>()
            .with_context(|| format!("invalid genome assembly: '{}'", self.genome_assembly))?;

        if self.ensembl.is_none() && self.refseq.is_none() {
            bail!("at least one of 'ensembl' or 'refseq' must be present");
        }

        for (section, name, entry) in self.file_entries() {
            validate_md5_optional(&format!("{section}/{name}"), &entry.md5)?;
        }

        Ok(())
    }

    /// Returns `(section, field_name, &FileEntry)` triples for all ensembl + refseq entries.
    /// HGNC is excluded — access it separately via `self.hgnc`.
    pub fn file_entries(&self) -> impl Iterator<Item = (&str, &str, &FileEntry)> {
        let ensembl_iter = self.ensembl.iter().flat_map(|e| {
            e.file_entries()
                .map(|(name, entry)| ("ensembl", name, entry))
        });
        let refseq_iter = self.refseq.iter().flat_map(|r| {
            r.file_entries()
                .map(|(name, entry)| ("refseq", name, entry))
        });
        ensembl_iter.chain(refseq_iter)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_config(json: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(json.as_bytes()).unwrap();
        f
    }

    #[test]
    fn valid_config_all_fields() {
        let json = r#"{
            "genomeAssembly": "GRCh38",
            "patchLevel": 14,
            "assemblyReport": { "url": "https://example.com/report.txt", "md5": "21f3ac4aa8245a99eb874082051b9dde" },
            "fasta": { "url": "https://example.com/genome.fna.gz", "md5": "c30471567037b2b2389d43c908c653e1" },
            "gff": { "url": "https://example.com/genes.gff.gz", "md5": "24b731562b9d4cae9e37b23404b5be16" },
            "ideogram": { "url": "https://example.com/ideogram", "md5": "043761ad1c7ba3c101d2423090603f60" }
        }"#;
        let f = write_config(json);
        let config = ReferenceConfig::from_file(f.path()).unwrap();
        assert_eq!(config.genome_assembly, "GRCh38");
        assert_eq!(config.patch_level, 14);
        assert_eq!(config.file_entries().count(), 4);
    }

    #[test]
    fn valid_config_optional_fields_omitted() {
        let json = r#"{
            "genomeAssembly": "GRCh37",
            "patchLevel": 13,
            "assemblyReport": { "url": "https://example.com/report.txt", "md5": "21f3ac4aa8245a99eb874082051b9dde" },
            "fasta": { "url": "https://example.com/genome.fna.gz", "md5": "c30471567037b2b2389d43c908c653e1" }
        }"#;
        let f = write_config(json);
        let config = ReferenceConfig::from_file(f.path()).unwrap();
        assert_eq!(config.genome_assembly, "GRCh37");
        assert!(config.gff.is_none());
        assert!(config.ideogram.is_none());
        assert_eq!(config.file_entries().count(), 2);
    }

    #[test]
    fn invalid_md5() {
        let json = r#"{
            "genomeAssembly": "GRCh38",
            "patchLevel": 14,
            "assemblyReport": { "url": "https://example.com/report.txt", "md5": "not_a_valid_md5" },
            "fasta": { "url": "https://example.com/genome.fna.gz", "md5": "c30471567037b2b2389d43c908c653e1" }
        }"#;
        let f = write_config(json);
        let err = ReferenceConfig::from_file(f.path()).unwrap_err();
        assert!(err.to_string().contains("invalid MD5"));
    }

    #[test]
    fn invalid_assembly_name() {
        let json = r#"{
            "genomeAssembly": "hg19",
            "patchLevel": 14,
            "assemblyReport": { "url": "https://example.com/report.txt", "md5": "21f3ac4aa8245a99eb874082051b9dde" },
            "fasta": { "url": "https://example.com/genome.fna.gz", "md5": "c30471567037b2b2389d43c908c653e1" }
        }"#;
        let f = write_config(json);
        let err = ReferenceConfig::from_file(f.path()).unwrap_err();
        assert!(err.to_string().contains("invalid genome assembly"));
    }

    // ── CacheConfig tests ───────────────────────────────────

    fn minimal_ensembl_json() -> &'static str {
        r#"{
            "name": "Ensembl", "version": "115", "description": "test", "releaseDate": "2025-09-02",
            "gff3": { "url": "https://example.com/e1.gz", "md5": "" },
            "regulatoryGff": { "url": "https://example.com/e2.gz", "md5": "" },
            "geneMysql": { "url": "https://example.com/e3.gz", "md5": "" },
            "transcriptMysql": { "url": "https://example.com/e4.gz", "md5": "" },
            "ncrnaFasta": { "url": "https://example.com/e5.gz", "md5": "" },
            "cdnaFasta": { "url": "https://example.com/e6.gz", "md5": "" },
            "cdsFasta": { "url": "https://example.com/e7.gz", "md5": "" },
            "peptideFasta": { "url": "https://example.com/e8.gz", "md5": "" }
        }"#
    }

    fn minimal_refseq_json() -> &'static str {
        r#"{
            "name": "RefSeq", "version": "v1", "description": "test", "releaseDate": "2025-08-06",
            "genbank": { "url": "https://example.com/r1.gz", "md5": "dba0915f2560e6b9d4943fac274816b9" },
            "gff3": { "url": "https://example.com/r2.gz", "md5": "24b731562b9d4cae9e37b23404b5be16" },
            "rnaFasta": { "url": "https://example.com/r3.gz", "md5": "b4a2ce202c90c0f24f22850c6bc7d774" },
            "proteinFasta": { "url": "https://example.com/r4.gz", "md5": "25a880236608ad046d7f9b65e895eb4b" }
        }"#
    }

    #[test]
    fn cache_config_both_sources() {
        let json = format!(
            r#"{{
                "genomeAssembly": "GRCh38",
                "ensembl": {ensembl},
                "hgnc": {{ "url": "https://example.com/hgnc.json" }},
                "refseq": {refseq}
            }}"#,
            ensembl = minimal_ensembl_json(),
            refseq = minimal_refseq_json()
        );
        let f = write_config(&json);
        let config = CacheConfig::from_file(f.path()).unwrap();
        assert_eq!(config.genome_assembly, "GRCh38");
        assert!(config.ensembl.is_some());
        assert!(config.refseq.is_some());
        // 8 ensembl + 4 refseq = 12
        assert_eq!(config.file_entries().count(), 12);
    }

    #[test]
    fn cache_config_refseq_only() {
        let json = format!(
            r#"{{
                "genomeAssembly": "GRCh37",
                "hgnc": {{ "url": "https://example.com/hgnc.json" }},
                "refseq": {refseq}
            }}"#,
            refseq = minimal_refseq_json()
        );
        let f = write_config(&json);
        let config = CacheConfig::from_file(f.path()).unwrap();
        assert!(config.ensembl.is_none());
        assert!(config.refseq.is_some());
        assert_eq!(config.file_entries().count(), 4);
    }

    #[test]
    fn cache_config_missing_both_sources() {
        let json = r#"{
            "genomeAssembly": "GRCh38",
            "hgnc": { "url": "https://example.com/hgnc.json" }
        }"#;
        let f = write_config(json);
        let err = CacheConfig::from_file(f.path()).unwrap_err();
        assert!(err.to_string().contains("at least one"));
    }

    #[test]
    fn cache_config_invalid_md5() {
        let json = r#"{
            "genomeAssembly": "GRCh38",
            "hgnc": { "url": "https://example.com/hgnc.json" },
            "refseq": {
                "name": "RefSeq", "version": "v1", "description": "test", "releaseDate": "2025-08-06",
                "genbank": { "url": "https://example.com/r1.gz", "md5": "not_valid" },
                "gff3": { "url": "https://example.com/r2.gz", "md5": "24b731562b9d4cae9e37b23404b5be16" },
                "rnaFasta": { "url": "https://example.com/r3.gz", "md5": "b4a2ce202c90c0f24f22850c6bc7d774" },
                "proteinFasta": { "url": "https://example.com/r4.gz", "md5": "25a880236608ad046d7f9b65e895eb4b" }
            }
        }"#;
        let f = write_config(json);
        let err = CacheConfig::from_file(f.path()).unwrap_err();
        assert!(err.to_string().contains("invalid MD5"));
    }

    #[test]
    fn cache_config_invalid_assembly() {
        let json = format!(
            r#"{{
                "genomeAssembly": "hg19",
                "hgnc": {{ "url": "https://example.com/hgnc.json" }},
                "refseq": {refseq}
            }}"#,
            refseq = minimal_refseq_json()
        );
        let f = write_config(&json);
        let err = CacheConfig::from_file(f.path()).unwrap_err();
        assert!(err.to_string().contains("invalid genome assembly"));
    }

    #[test]
    fn cache_config_empty_md5_passes_validation() {
        let json = format!(
            r#"{{
                "genomeAssembly": "GRCh38",
                "ensembl": {ensembl},
                "hgnc": {{ "url": "https://example.com/hgnc.json" }}
            }}"#,
            ensembl = minimal_ensembl_json()
        );
        let f = write_config(&json);
        let config = CacheConfig::from_file(f.path()).unwrap();
        assert!(config.ensembl.is_some());
        // All 8 ensembl entries have empty MD5 — validation should pass
        assert_eq!(config.file_entries().count(), 8);
    }
}
