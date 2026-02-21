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
}
