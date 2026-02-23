//! Ensembl GFF3 line and attribute parser.

use std::collections::HashMap;
use std::sync::LazyLock;

use crate::biotype::{BioType, BioTypeCategory};
use crate::error::Error;
use crate::strand::Strand;

use super::entry::{Gff3Attributes, Gff3Entry};
use super::refseq_parser::ParsedLine;

/// Attribute keys that are recognized but not needed for transcript construction or evaluation.
/// These are metadata fields from Ensembl GFF3 that don't affect the cache output. Unknown keys
/// are treated as errors to catch schema changes, so all expected keys must be listed here or
/// parsed in the match arms below.
static SKIPPED_KEYS: LazyLock<std::collections::HashSet<&'static str>> = LazyLock::new(|| {
    [
        "Alias",
        "ccdsid",
        "constitutive",
        "description",
        "ensembl_end_phase",
        "external_name",
        "havana_transcript",
        "havana_version",
        "logic_name",
        "rank",
        "bound_end",
        "bound_start",
        "color",
        "extended_end",
        "extended_start",
        "feature_type",
        "gene_biotype",
        "gene_name",
        "transcript_support_level",
    ]
    .into_iter()
    .collect()
});

/// Parse a single Ensembl GFF3 line into a structured entry.
pub fn parse_line(line: &str, name_to_index: &HashMap<String, usize>) -> Result<ParsedLine, Error> {
    // Comments and directives
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

    // Split into 9 tab columns
    let columns: Vec<&str> = line.split('\t').collect();
    if columns.len() != 9 {
        return Err(Error::Parse(format!(
            "GFF3 line has {} columns, expected 9",
            columns.len()
        )));
    }

    // Column 1: chromosome
    let chr_name = columns[0];
    let chromosome_index = match name_to_index.get(chr_name) {
        Some(&idx) => idx,
        None => return Ok(ParsedLine::Discarded),
    };

    // Column 3: biotype (strict — unknown type is an error)
    let biotype: BioType = columns[2].parse()?;

    // Discard not-useful and unsupported
    let category = biotype.category();
    if matches!(
        category,
        BioTypeCategory::NotUseful | BioTypeCategory::Unsupported
    ) {
        return Ok(ParsedLine::Discarded);
    }

    // Column 4 & 5: start and end
    let start: i32 = columns[3]
        .parse()
        .map_err(|e| Error::Parse(format!("invalid start '{}': {e}", columns[3])))?;
    let end: i32 = columns[4]
        .parse()
        .map_err(|e| Error::Parse(format!("invalid end '{}': {e}", columns[4])))?;

    // Column 7: strand
    let strand = Strand::from_gff3(columns[6]);

    // Column 9: attributes
    let attributes = parse_attributes(columns[8])?;

    // ID is required for non-discarded entries
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

/// Parse Ensembl GFF3 column 9 attributes.
fn parse_attributes(attrs_str: &str) -> Result<Gff3Attributes, Error> {
    let mut attrs = Gff3Attributes::default();
    let mut exon_id: Option<&str> = None;

    for pair in attrs_str.split(';') {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }

        let eq_pos = pair
            .find('=')
            .ok_or_else(|| Error::Parse(format!("attribute missing '=': '{pair}'")))?;
        let key = &pair[..eq_pos];
        let value = &pair[eq_pos + 1..];

        match key {
            "ID" => attrs.id = value.to_string(),
            "Parent" => attrs.parent_id = Some(value.to_string()),
            "Name" => attrs.name = Some(value.to_string()),
            "gene_id" => attrs.gene_id = Some(value.to_string()),
            "transcript_id" => attrs.transcript_id = Some(value.to_string()),
            "protein_id" => attrs.protein_id = Some(value.to_string()),
            "version" => {
                attrs.version = Some(
                    value
                        .parse::<u16>()
                        .map_err(|e| Error::Parse(format!("invalid version '{value}': {e}")))?,
                );
            }
            "biotype" => attrs.biotype_attr = Some(value.to_string()),
            "tag" => {
                attrs.is_ensembl_canonical = value.split(',').any(|t| t == "Ensembl_canonical");
            }
            "ensembl_phase" => {} // recognized but not used
            "exon_id" => exon_id = Some(value),
            "Note" | "note" => attrs.note = Some(value.to_string()),
            _ => {
                if !SKIPPED_KEYS.contains(key) {
                    return Err(Error::Parse(format!(
                        "unrecognized GFF3 attribute key: '{key}'"
                    )));
                }
            }
        }
    }

    // Ensembl exons sometimes lack ID= but always have exon_id=
    if attrs.id.is_empty()
        && let Some(eid) = exon_id
    {
        attrs.id = format!("exon:{eid}");
    }

    Ok(attrs)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_name_map() -> HashMap<String, usize> {
        let mut m = HashMap::new();
        m.insert("1".to_string(), 0);
        m.insert("17".to_string(), 16);
        m
    }

    #[test]
    fn parse_ensembl_gene_line() {
        let line = "1\tensembl\tgene\t69091\t70008\t.\t+\t.\tID=gene:ENSG00000186092;Name=OR4F5;gene_id=ENSG00000186092;biotype=protein_coding;version=14";
        let result = parse_line(line, &test_name_map()).unwrap();
        match result {
            ParsedLine::Entry(e) => {
                assert_eq!(e.chromosome_index, 0);
                assert_eq!(e.start, 69091);
                assert_eq!(e.end, 70008);
                assert_eq!(e.biotype, BioType::Gene);
                assert_eq!(e.strand, Strand::Forward);
                assert_eq!(e.attributes.id, "gene:ENSG00000186092");
                assert_eq!(e.attributes.name.as_deref(), Some("OR4F5"));
                assert_eq!(e.attributes.gene_id.as_deref(), Some("ENSG00000186092"));
                assert_eq!(e.attributes.biotype_attr.as_deref(), Some("protein_coding"));
                assert_eq!(e.attributes.version, Some(14));
            }
            _ => panic!("expected Entry"),
        }
    }

    #[test]
    fn parse_ensembl_transcript_line() {
        let line = "1\tensembl\tmRNA\t69091\t70008\t.\t+\t.\tID=transcript:ENST00000335137;Parent=gene:ENSG00000186092;Name=OR4F5-201;transcript_id=ENST00000335137;biotype=protein_coding;version=4";
        let result = parse_line(line, &test_name_map()).unwrap();
        match result {
            ParsedLine::Entry(e) => {
                assert_eq!(e.attributes.id, "transcript:ENST00000335137");
                assert_eq!(
                    e.attributes.parent_id.as_deref(),
                    Some("gene:ENSG00000186092")
                );
                assert_eq!(
                    e.attributes.transcript_id.as_deref(),
                    Some("ENST00000335137")
                );
                assert_eq!(e.attributes.version, Some(4));
            }
            _ => panic!("expected Entry"),
        }
    }

    #[test]
    fn parse_ensembl_cds_line() {
        let line = "1\tensembl\tCDS\t69091\t69984\t.\t+\t0\tID=CDS:ENSP00000334393;Parent=transcript:ENST00000335137;protein_id=ENSP00000334393";
        let result = parse_line(line, &test_name_map()).unwrap();
        match result {
            ParsedLine::Entry(e) => {
                assert_eq!(e.biotype, BioType::Cds);
                assert_eq!(e.attributes.protein_id.as_deref(), Some("ENSP00000334393"));
                assert_eq!(
                    e.attributes.parent_id.as_deref(),
                    Some("transcript:ENST00000335137")
                );
            }
            _ => panic!("expected Entry"),
        }
    }

    #[test]
    fn parse_ensembl_exon_line() {
        let line = "1\tensembl\texon\t69091\t70008\t.\t+\t.\tParent=transcript:ENST00000335137;exon_id=ENSE00002319515;version=1;ID=exon:ENSE00002319515";
        let result = parse_line(line, &test_name_map()).unwrap();
        match result {
            ParsedLine::Entry(e) => {
                assert_eq!(e.biotype, BioType::Exon);
                assert_eq!(e.attributes.id, "exon:ENSE00002319515");
                assert_eq!(
                    e.attributes.parent_id.as_deref(),
                    Some("transcript:ENST00000335137")
                );
            }
            _ => panic!("expected Entry"),
        }
    }

    #[test]
    fn exon_without_id_uses_exon_id_fallback() {
        // Real Ensembl exon line that lacks ID= but has exon_id=
        let line = "1\thavana_tagene\texon\t11121\t11211\t.\t+\t.\tParent=transcript:ENST00000832824;Name=ENSE00004248723;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00004248723;rank=1;version=1";
        let result = parse_line(line, &test_name_map()).unwrap();
        match result {
            ParsedLine::Entry(e) => {
                assert_eq!(e.biotype, BioType::Exon);
                assert_eq!(e.attributes.id, "exon:ENSE00004248723");
                assert_eq!(
                    e.attributes.parent_id.as_deref(),
                    Some("transcript:ENST00000832824")
                );
            }
            _ => panic!("expected Entry"),
        }
    }

    #[test]
    fn unknown_chromosome_discarded() {
        let line = "chrUn\tensembl\tgene\t100\t200\t.\t+\t.\tID=gene:ENSG00000000001;gene_id=ENSG00000000001;biotype=protein_coding;version=1";
        let result = parse_line(line, &test_name_map()).unwrap();
        assert!(matches!(result, ParsedLine::Discarded));
    }

    #[test]
    fn not_useful_biotype_discarded() {
        let line = "1\tensembl\tbiological_region\t100\t200\t.\t+\t.\tID=biological_region:1";
        let result = parse_line(line, &test_name_map()).unwrap();
        assert!(matches!(result, ParsedLine::Discarded));
    }

    #[test]
    fn unknown_attribute_key_error() {
        let line = "1\tensembl\tgene\t100\t200\t.\t+\t.\tID=gene:ENSG00000000001;totally_unknown_key=value";
        assert!(parse_line(line, &test_name_map()).is_err());
    }

    #[test]
    fn skipped_attribute_key() {
        let line = "1\tensembl\tgene\t100\t200\t.\t+\t.\tID=gene:ENSG00000000001;description=some protein;logic_name=ensembl";
        let result = parse_line(line, &test_name_map()).unwrap();
        assert!(matches!(result, ParsedLine::Entry(_)));
    }

    #[test]
    fn ensembl_canonical_tag_parsed() {
        let line = "1\tensembl\tmRNA\t69091\t70008\t.\t+\t.\tID=transcript:ENST00000335137;Parent=gene:ENSG00000186092;transcript_id=ENST00000335137;biotype=protein_coding;tag=gencode_basic,gencode_primary,Ensembl_canonical;version=4";
        let result = parse_line(line, &test_name_map()).unwrap();
        match result {
            ParsedLine::Entry(e) => {
                assert!(e.attributes.is_ensembl_canonical);
            }
            _ => panic!("expected Entry"),
        }
    }

    #[test]
    fn non_canonical_tag_not_set() {
        let line = "1\tensembl\tmRNA\t69091\t70008\t.\t+\t.\tID=transcript:ENST00000335137;Parent=gene:ENSG00000186092;transcript_id=ENST00000335137;biotype=protein_coding;tag=gencode_basic;version=4";
        let result = parse_line(line, &test_name_map()).unwrap();
        match result {
            ParsedLine::Entry(e) => {
                assert!(!e.attributes.is_ensembl_canonical);
            }
            _ => panic!("expected Entry"),
        }
    }

    #[test]
    fn comment_and_section_end() {
        assert!(matches!(
            parse_line("# comment", &test_name_map()).unwrap(),
            ParsedLine::Comment
        ));
        assert!(matches!(
            parse_line("###", &test_name_map()).unwrap(),
            ParsedLine::EndOfSection
        ));
    }
}
