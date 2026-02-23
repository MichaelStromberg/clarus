//! HGNC complete set parser: builds Ensembl gene ID → HGNC numeric ID lookup.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};

use crate::error::Error;

/// Parse the HGNC complete set TSV file and build a mapping from Ensembl gene ID
/// to HGNC numeric ID.
///
/// The file is plain TSV (not gzipped) with a header row. Relevant columns (0-indexed):
/// - Column 0: `hgnc_id` — format `HGNC:nnnnn`, strip prefix, parse as i32
/// - Column 19: `ensembl_gene_id` — e.g. `ENSG00000250722` (may be empty)
///
/// Rows with empty `ensembl_gene_id` are skipped.
pub fn parse_hgnc_tsv<R: Read>(reader: R) -> Result<HashMap<String, i32>, Error> {
    let buf = BufReader::new(reader);
    let mut map = HashMap::new();
    let mut first = true;

    for line in buf.lines() {
        let line = line?;

        // Skip header row
        if first {
            first = false;
            continue;
        }

        let mut cols = line.split('\t');
        let hgnc_raw = match cols.next() {
            Some(s) => s,
            None => continue,
        };

        // Advance to column 19
        let ensembl_gene_id = match cols.nth(18) {
            Some(s) if !s.is_empty() => s,
            _ => continue,
        };

        // Strip "HGNC:" prefix and parse numeric ID
        let hgnc_id: i32 = hgnc_raw
            .strip_prefix("HGNC:")
            .ok_or_else(|| Error::Parse(format!("expected HGNC: prefix in '{hgnc_raw}'")))?
            .parse()
            .map_err(|e| Error::Parse(format!("invalid HGNC ID '{hgnc_raw}': {e}")))?;

        map.insert(ensembl_gene_id.to_string(), hgnc_id);
    }

    Ok(map)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_header() -> &'static str {
        "hgnc_id\tsymbol\tname\tlocus_group\tlocus_type\tstatus\tlocation\tlocation_sortable\talias_symbol\talias_name\tprev_symbol\tprev_name\tgene_family\tgene_family_id\tdate_approved_reserved\tdate_symbol_changed\tdate_name_changed\tdate_modified\tentrez_id\tensembl_gene_id\tvega_id"
    }

    #[test]
    fn parse_basic() {
        let data = format!(
            "{}\n\
             HGNC:10751\tSELENOP\tSelenoprotein P\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tENSG00000250722\t\n\
             HGNC:1234\tFOO\tFoo gene\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tENSG00000000001\t\n\
             HGNC:5678\tBAR\tBar gene\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tENSG00000000002\t\n",
            make_header()
        );
        let map = parse_hgnc_tsv(data.as_bytes()).unwrap();
        assert_eq!(map.len(), 3);
        assert_eq!(map["ENSG00000250722"], 10751);
        assert_eq!(map["ENSG00000000001"], 1234);
        assert_eq!(map["ENSG00000000002"], 5678);
    }

    #[test]
    fn empty_ensembl_field_skipped() {
        let data = format!(
            "{}\n\
             HGNC:10751\tSELENOP\tSelenoprotein P\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tENSG00000250722\t\n\
             HGNC:9999\tNOENS\tNo ensembl\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n",
            make_header()
        );
        let map = parse_hgnc_tsv(data.as_bytes()).unwrap();
        assert_eq!(map.len(), 1);
        assert!(map.contains_key("ENSG00000250722"));
        assert!(!map.contains_key(""));
    }

    #[test]
    fn strips_hgnc_prefix() {
        let data = format!(
            "{}\n\
             HGNC:10751\tSELENOP\tSelenoprotein P\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tENSG00000250722\t\n",
            make_header()
        );
        let map = parse_hgnc_tsv(data.as_bytes()).unwrap();
        assert_eq!(*map.get("ENSG00000250722").unwrap(), 10751);
    }
}
