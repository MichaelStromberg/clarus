//! JSON output serde types per specification 16.

use serde::Serialize;

use crate::variant::types::VariantType;

/// Metadata written to `{prefix}_metadata.json.gz`.
#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Metadata {
    pub annotator: String,
    pub creation_time: String,
    pub genome_assembly: String,
    pub schema_version: u32,
    pub data_sources: Vec<DataSource>,
    pub samples: Vec<String>,
    pub input: InputInfo,
    pub stats: AnnotationStats,
}

/// VCF input provenance.
#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct InputInfo {
    pub file_name: String,
    pub vcf_version: String,
    pub contigs: usize,
    pub samples: usize,
    pub inferred_assembly: String,
}

/// Annotation statistics, filled after the annotation loop completes.
#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AnnotationStats {
    pub positions: u64,
    pub variants: u64,
    pub skipped_records: u64,
}

/// A data source referenced in the metadata.
#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct DataSource {
    pub name: String,
    pub version: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub release_date: Option<String>,
}

/// One per line in the genes JSONL file.
#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Gene {
    pub name: String,
    pub hgnc_id: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ncbi_gene_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ensembl_gene_id: Option<String>,
}

/// A single genomic position (one VCF record).
#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Position {
    pub chromosome: String,
    pub position: i64,
    pub ref_allele: String,
    pub alt_alleles: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub quality: Option<f64>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub filters: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub strand_bias: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fisher_strand_bias: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mapping_quality: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub impute_score: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub allele_frequency: Option<Vec<f64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ref_panel_allele_frequency: Option<Vec<f64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub locus: Option<String>,
    pub samples: Vec<JsonSample>,
    pub variants: Vec<JsonVariant>,
}

/// A sample in the JSON output.
#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct JsonSample {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genotype: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub variant_frequencies: Option<Vec<f64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub total_depth: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genotype_quality: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub allele_depths: Option<Vec<i32>>,
    #[serde(skip_serializing_if = "is_false")]
    pub failed_filter: bool,
    #[serde(skip_serializing_if = "is_false")]
    pub is_de_novo: bool,
    #[serde(skip_serializing_if = "is_false")]
    pub loss_of_heterozygosity: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bin_count: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub is_imputed_genotype: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genotype_dosage: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genotype_posteriors: Option<Vec<f64>>,
}

/// A variant (one per ALT allele) in the JSON output.
#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct JsonVariant {
    pub vid: String,
    pub chromosome: String,
    pub begin: i64,
    pub end: i64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub variant_type: VariantType,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hgvsg: Option<String>,
}

fn is_false(v: &bool) -> bool {
    !v
}

impl From<crate::vcf::sample::VcfSample> for JsonSample {
    fn from(s: crate::vcf::sample::VcfSample) -> Self {
        Self {
            genotype: s.genotype,
            variant_frequencies: s.variant_frequencies,
            total_depth: s.total_depth,
            genotype_quality: s.genotype_quality,
            allele_depths: s.allele_depths,
            failed_filter: s.failed_filter,
            is_de_novo: s.is_de_novo,
            loss_of_heterozygosity: s.loss_of_heterozygosity,
            bin_count: s.bin_count,
            is_imputed_genotype: None,
            genotype_dosage: s.genotype_dosage,
            genotype_posteriors: s.genotype_posteriors,
        }
    }
}
