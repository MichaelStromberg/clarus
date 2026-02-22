//! GFF3 line and attribute parser.

use std::collections::HashMap;
use std::sync::LazyLock;

use crate::biotype::{BioType, BioTypeCategory};
use crate::error::Error;
use crate::strand::Strand;

use super::entry::{CigarOp, CigarOpType, Gff3Attributes, Gff3Entry};

/// Result of parsing a single GFF3 line.
pub enum ParsedLine {
    Entry(Box<Gff3Entry>),
    Discarded,
    Comment,
    EndOfSection,
}

/// Attribute keys that are recognized but silently ignored.
static SKIPPED_KEYS: LazyLock<std::collections::HashSet<&'static str>> = LazyLock::new(|| {
    [
        "align_id",
        "allele",
        "anticodon",
        "assembly_bases_aln",
        "assembly_bases_seq",
        "batch_id",
        "bio-material",
        "bit_score",
        "blast_aligner",
        "blast_score",
        "bound_moiety",
        "cell-line",
        "chromosome",
        "codons",
        "common_component",
        "consensus_splices",
        "country",
        "crc32",
        "curated_alignment",
        "description",
        "direction",
        "end_range",
        "exception",
        "exon_identity",
        "exon_number",
        "e_value",
        "feat_class",
        "filter_score",
        "for_remapping",
        "function",
        "gap_count",
        "gbkey",
        "gene_biotype",
        "gene_synonym",
        "genome",
        "hsp_percent_coverage",
        "identity",
        "idty",
        "inference",
        "inversion_merge_aligner",
        "isolate",
        "isolation-source",
        "Is_circular",
        "lxr_locAcc_currStat_120",
        "lxr_locAcc_currStat_35",
        "map",
        "matchable_bases",
        "matched_bases",
        "matches",
        "merge_aligner",
        "mobile_element_type",
        "model_evidence",
        "mol_type",
        "not_for_annotation",
        "number",
        "num_ident",
        "num_mismatch",
        "part",
        "partial",
        "pct_coverage",
        "pct_coverage_hiqual",
        "pct_identity_gap",
        "pct_identity_gapopen_only",
        "pct_identity_ungap",
        "product",
        "product_coverage",
        "promoted_rank",
        "protein_id",
        "pseudo",
        "qtaxid",
        "rank",
        "recombination_class",
        "regulatory_class",
        "rpt_family",
        "rpt_type",
        "rpt_unit_range",
        "rpt_unit_seq",
        "satellite",
        "sex",
        "splices",
        "standard_name",
        "start_range",
        "tissue-type",
        "transl_except",
        "transl_table",
        "weighted_identity",
    ]
    .into_iter()
    .collect()
});

/// Parse a single GFF3 line into a structured entry.
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

/// Parse GFF3 column 9 attributes.
fn parse_attributes(attrs_str: &str) -> Result<Gff3Attributes, Error> {
    let mut attrs = Gff3Attributes::default();

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
            "gene" => attrs.gene_symbol = Some(value.to_string()),
            "Dbxref" => parse_dbxref(value, &mut attrs)?,
            "tag" => {
                if value == "RefSeq Select" {
                    attrs.is_refseq_select = true;
                } else if value == "MANE Select" {
                    attrs.is_mane_select = true;
                }
            }
            "Target" => parse_target(value, &mut attrs)?,
            "Gap" => attrs.cigar_ops = Some(parse_cigar(value)?),
            "transcript_id" => attrs.transcript_id = Some(value.to_string()),
            "Note" | "note" => attrs.note = Some(url_decode(value)),
            "experiment" => parse_experiment(value, &mut attrs)?,
            _ => {
                if !SKIPPED_KEYS.contains(key) {
                    return Err(Error::Parse(format!(
                        "unrecognized GFF3 attribute key: '{key}'"
                    )));
                }
            }
        }
    }

    Ok(attrs)
}

/// Parse Dbxref value: comma-separated database:identifier pairs.
/// Uses last-colon rule for splitting.
fn parse_dbxref(value: &str, attrs: &mut Gff3Attributes) -> Result<(), Error> {
    for xref in value.split(',') {
        let colon_pos = xref
            .rfind(':')
            .ok_or_else(|| Error::Parse(format!("Dbxref entry missing colon: '{xref}'")))?;
        let db_key = &xref[..colon_pos];
        let db_value = &xref[colon_pos + 1..];

        match db_key {
            "GeneID" => attrs.gene_id = Some(db_value.to_string()),
            "HGNC:HGNC" => {
                let id: i32 = db_value
                    .parse()
                    .map_err(|e| Error::Parse(format!("invalid HGNC ID '{db_value}': {e}")))?;
                attrs.hgnc_id = Some(id);
            }
            _ => {} // Silently ignore other databases
        }
    }
    Ok(())
}

/// Parse Target attribute: "transcript_id start end [strand]"
fn parse_target(value: &str, attrs: &mut Gff3Attributes) -> Result<(), Error> {
    let parts: Vec<&str> = value.split(' ').collect();
    if parts.len() < 3 {
        return Err(Error::Parse(format!(
            "Target attribute has {} components, expected at least 3: '{value}'",
            parts.len()
        )));
    }
    attrs.target_id = Some(parts[0].to_string());
    attrs.target_start = Some(
        parts[1]
            .parse()
            .map_err(|e| Error::Parse(format!("invalid Target start '{}': {e}", parts[1])))?,
    );
    attrs.target_end = Some(
        parts[2]
            .parse()
            .map_err(|e| Error::Parse(format!("invalid Target end '{}': {e}", parts[2])))?,
    );
    Ok(())
}

/// Parse CIGAR operations from Gap attribute value.
pub fn parse_cigar(value: &str) -> Result<Vec<CigarOp>, Error> {
    let mut ops = Vec::new();
    for token in value.split(' ') {
        if token.is_empty() {
            continue;
        }
        let op_char = token.as_bytes()[0];
        let op_type = match op_char {
            b'M' => CigarOpType::Match,
            b'I' => CigarOpType::Insertion,
            b'D' => CigarOpType::Deletion,
            _ => {
                return Err(Error::Parse(format!(
                    "invalid CIGAR operation type: '{}'",
                    op_char as char
                )));
            }
        };
        let length: u32 = token[1..]
            .parse()
            .map_err(|e| Error::Parse(format!("invalid CIGAR length in '{token}': {e}")))?;
        ops.push(CigarOp { op_type, length });
    }
    Ok(ops)
}

/// URL-decode Note/note values.
fn url_decode(value: &str) -> String {
    value
        .replace("%25", "%")
        .replace("%2C", ",")
        .replace("%3B", ";")
        .replace("%3D", "=")
}

/// Parse experiment attribute for ECO ID and PubMed IDs.
fn parse_experiment(value: &str, attrs: &mut Gff3Attributes) -> Result<(), Error> {
    let bytes = value.as_bytes();
    let len = bytes.len();
    let mut i = 0;
    let mut eco_id: Option<i32> = None;
    let mut pubmed_ids: Vec<i32> = Vec::new();

    while i < len {
        // Find next '['
        match bytes[i..].iter().position(|&b| b == b'[') {
            Some(offset) => i += offset + 1,
            None => break,
        }

        // Find closing ']'
        let close = bytes[i..]
            .iter()
            .position(|&b| b == b']')
            .ok_or_else(|| Error::Parse("experiment attribute: unclosed bracket".to_string()))?;
        let content = &value[i..i + close];
        i += close + 1;

        if content.starts_with("PMID") {
            let decoded = content.replace("%2C", ",");
            for part in decoded.split(", ") {
                if part.len() > 5 {
                    let id: i32 = part[5..]
                        .parse()
                        .map_err(|e| Error::Parse(format!("invalid PubMed ID in '{part}': {e}")))?;
                    pubmed_ids.push(id);
                }
            }
        } else if content.starts_with("ECO") {
            let id: i32 = content[4..]
                .parse()
                .map_err(|e| Error::Parse(format!("invalid ECO ID in '{content}': {e}")))?;
            eco_id = Some(id);
        } else {
            return Err(Error::Parse(format!(
                "experiment attribute: unrecognized bracket content: '{content}'"
            )));
        }
    }

    match eco_id {
        Some(id) => attrs.eco_id = Some(id),
        None => {
            return Err(Error::Parse(
                "experiment attribute: no ECO ID found".to_string(),
            ));
        }
    }

    if !pubmed_ids.is_empty() {
        attrs.pubmed_ids = Some(pubmed_ids);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_name_map() -> HashMap<String, usize> {
        let mut m = HashMap::new();
        m.insert("NC_000001.11".to_string(), 0);
        m.insert("NC_000002.12".to_string(), 1);
        m
    }

    #[test]
    fn parse_gene_line() {
        let line = "NC_000001.11\tBestRefSeq\tgene\t11874\t14409\t.\t+\t.\tID=gene-DDX11L1;Dbxref=GeneID:100287102,HGNC:HGNC:37102;gene=DDX11L1";
        let result = parse_line(line, &test_name_map()).unwrap();
        match result {
            ParsedLine::Entry(e) => {
                assert_eq!(e.chromosome_index, 0);
                assert_eq!(e.start, 11874);
                assert_eq!(e.end, 14409);
                assert_eq!(e.biotype, BioType::Gene);
                assert_eq!(e.strand, Strand::Forward);
                assert_eq!(e.attributes.id, "gene-DDX11L1");
                assert_eq!(e.attributes.gene_id.as_deref(), Some("100287102"));
                assert_eq!(e.attributes.hgnc_id, Some(37102));
                assert_eq!(e.attributes.gene_symbol.as_deref(), Some("DDX11L1"));
            }
            _ => panic!("expected Entry"),
        }
    }

    #[test]
    fn parse_comment() {
        let result = parse_line("# comment", &test_name_map()).unwrap();
        assert!(matches!(result, ParsedLine::Comment));
    }

    #[test]
    fn parse_end_section() {
        let result = parse_line("###", &test_name_map()).unwrap();
        assert!(matches!(result, ParsedLine::EndOfSection));
    }

    #[test]
    fn unknown_chromosome_discarded() {
        let line = "chrUn\tBestRefSeq\tgene\t100\t200\t.\t+\t.\tID=gene-X;gene=X";
        let result = parse_line(line, &test_name_map()).unwrap();
        assert!(matches!(result, ParsedLine::Discarded));
    }

    #[test]
    fn not_useful_biotype_discarded() {
        let line = "NC_000001.11\tBestRefSeq\tbiological_region\t100\t200\t.\t+\t.\tID=id-BR1";
        let result = parse_line(line, &test_name_map()).unwrap();
        assert!(matches!(result, ParsedLine::Discarded));
    }

    #[test]
    fn unknown_biotype_error() {
        let line = "NC_000001.11\tBestRefSeq\tunknown_type\t100\t200\t.\t+\t.\tID=id-X";
        assert!(parse_line(line, &test_name_map()).is_err());
    }

    #[test]
    fn url_decode_note() {
        let decoded = url_decode("100%25 of items%2C separated%3B equal%3D sign");
        assert_eq!(decoded, "100% of items, separated; equal= sign");
    }

    #[test]
    fn parse_target_attribute() {
        let mut attrs = Gff3Attributes::default();
        parse_target("NM_001005484.2 1 2618 +", &mut attrs).unwrap();
        assert_eq!(attrs.target_id.as_deref(), Some("NM_001005484.2"));
        assert_eq!(attrs.target_start, Some(1));
        assert_eq!(attrs.target_end, Some(2618));
    }

    #[test]
    fn parse_cigar_ops() {
        let ops = parse_cigar("M345 I2 M100 D1 M50").unwrap();
        assert_eq!(ops.len(), 5);
        assert_eq!(
            ops[0],
            CigarOp {
                op_type: CigarOpType::Match,
                length: 345
            }
        );
        assert_eq!(
            ops[1],
            CigarOp {
                op_type: CigarOpType::Insertion,
                length: 2
            }
        );
        assert_eq!(
            ops[2],
            CigarOp {
                op_type: CigarOpType::Match,
                length: 100
            }
        );
        assert_eq!(
            ops[3],
            CigarOp {
                op_type: CigarOpType::Deletion,
                length: 1
            }
        );
        assert_eq!(
            ops[4],
            CigarOp {
                op_type: CigarOpType::Match,
                length: 50
            }
        );
    }

    #[test]
    fn parse_dbxref() {
        let mut attrs = Gff3Attributes::default();
        super::parse_dbxref("GeneID:672,HGNC:HGNC:1100,MIM:113705", &mut attrs).unwrap();
        assert_eq!(attrs.gene_id.as_deref(), Some("672"));
        assert_eq!(attrs.hgnc_id, Some(1100));
    }

    #[test]
    fn parse_experiment_attr() {
        let mut attrs = Gff3Attributes::default();
        parse_experiment(
            "COORDINATES: polyA evidence [ECO:0006100][PMID:29474683, PMID:32163032]",
            &mut attrs,
        )
        .unwrap();
        assert_eq!(attrs.eco_id, Some(6100));
        assert_eq!(attrs.pubmed_ids, Some(vec![29474683, 32163032]));
    }

    #[test]
    fn parse_experiment_eco_only() {
        let mut attrs = Gff3Attributes::default();
        parse_experiment("COORDINATES: evidence [ECO:0006100]", &mut attrs).unwrap();
        assert_eq!(attrs.eco_id, Some(6100));
        assert!(attrs.pubmed_ids.is_none());
    }

    #[test]
    fn parse_experiment_no_eco_error() {
        let mut attrs = Gff3Attributes::default();
        assert!(parse_experiment("[PMID:12345]", &mut attrs).is_err());
    }

    #[test]
    fn tag_refseq_select() {
        let line = "NC_000001.11\tBestRefSeq\tmRNA\t100\t200\t.\t+\t.\tID=rna-NM_001.1;Parent=gene-X;Name=NM_001.1;tag=RefSeq Select";
        let result = parse_line(line, &test_name_map()).unwrap();
        match result {
            ParsedLine::Entry(e) => {
                assert!(e.attributes.is_refseq_select);
                assert!(!e.attributes.is_mane_select);
            }
            _ => panic!("expected Entry"),
        }
    }

    #[test]
    fn tag_mane_select() {
        let line = "NC_000001.11\tBestRefSeq\tmRNA\t100\t200\t.\t+\t.\tID=rna-NM_001.1;Parent=gene-X;Name=NM_001.1;tag=MANE Select";
        let result = parse_line(line, &test_name_map()).unwrap();
        match result {
            ParsedLine::Entry(e) => {
                assert!(e.attributes.is_mane_select);
            }
            _ => panic!("expected Entry"),
        }
    }

    #[test]
    fn unknown_attribute_key_error() {
        let line = "NC_000001.11\tBestRefSeq\tgene\t100\t200\t.\t+\t.\tID=gene-X;totally_unknown_key=value";
        assert!(parse_line(line, &test_name_map()).is_err());
    }

    #[test]
    fn skipped_attribute_key() {
        let line = "NC_000001.11\tBestRefSeq\tgene\t100\t200\t.\t+\t.\tID=gene-X;gene_biotype=protein_coding;gbkey=Gene";
        let result = parse_line(line, &test_name_map()).unwrap();
        assert!(matches!(result, ParsedLine::Entry(_)));
    }
}
