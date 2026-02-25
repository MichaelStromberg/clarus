//! HGNC complete set parser: builds a multi-indexed gene database.
//!
//! Parses the HGNC complete set TSV file with dynamic column detection from the
//! header row. Provides lookups by HGNC ID, Ensembl gene ID, and NCBI gene ID.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};

use crate::error::Error;

/// A single HGNC gene entry.
pub struct HgncEntry {
    pub hgnc_id: u32,
    pub symbol: String,
    pub entrez_id: Option<String>,
    pub ensembl_gene_id: Option<String>,
    pub refseq_accession: Option<Vec<String>>,
}

/// Multi-indexed HGNC gene database.
pub struct HgncDatabase {
    entries: Vec<HgncEntry>,
    by_hgnc_id: HashMap<u32, usize>,
    by_ensembl_id: HashMap<String, usize>,
    by_ncbi_id: HashMap<String, usize>,
}

impl HgncDatabase {
    /// Parse from the HGNC complete set TSV file.
    ///
    /// Column positions are determined dynamically from the header row by
    /// searching for column names: `hgnc_id`, `symbol`, `entrez_id`,
    /// `ensembl_gene_id`, `refseq_accession`.
    ///
    /// The `refseq_accession` column is pipe-delimited (e.g. `NM_007294|NM_007295`).
    ///
    /// Index behavior:
    /// - `by_hgnc_id`: duplicate HGNC IDs produce an error
    /// - `by_ensembl_id`: last-write-wins (silently overwrites)
    /// - `by_ncbi_id`: last-write-wins (silently overwrites)
    pub fn from_tsv<R: Read>(reader: R) -> Result<Self, Error> {
        let buf = BufReader::new(reader);
        let mut lines = buf.lines();

        // Parse header to find column positions
        let header = lines
            .next()
            .ok_or_else(|| Error::Parse("HGNC TSV: empty file".to_string()))??;

        let columns: Vec<&str> = header.split('\t').collect();
        let col_hgnc_id = find_column(&columns, "hgnc_id")?;
        let col_symbol = find_column(&columns, "symbol")?;
        let col_entrez_id = find_column(&columns, "entrez_id")?;
        let col_ensembl_gene_id = find_column(&columns, "ensembl_gene_id")?;
        let col_refseq_accession = find_column(&columns, "refseq_accession")?;

        let max_col = [
            col_hgnc_id,
            col_symbol,
            col_entrez_id,
            col_ensembl_gene_id,
            col_refseq_accession,
        ]
        .into_iter()
        .max()
        .expect("non-empty array");

        let mut entries = Vec::new();
        let mut by_hgnc_id = HashMap::new();
        let mut by_ensembl_id = HashMap::new();
        let mut by_ncbi_id = HashMap::new();

        for line in lines {
            let line = line?;
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() <= max_col {
                continue;
            }

            // Parse hgnc_id: strip optional "HGNC:" prefix
            let hgnc_raw = cols[col_hgnc_id];
            let hgnc_id: u32 = hgnc_raw
                .strip_prefix("HGNC:")
                .unwrap_or(hgnc_raw)
                .parse()
                .map_err(|e| Error::Parse(format!("invalid HGNC ID '{hgnc_raw}': {e}")))?;

            let symbol = cols[col_symbol].to_string();

            let entrez_id = non_empty(cols[col_entrez_id]);
            let ensembl_gene_id = non_empty(cols[col_ensembl_gene_id]);

            let refseq_accession = {
                let s = cols[col_refseq_accession];
                if s.is_empty() {
                    None
                } else {
                    Some(s.split('|').map(String::from).collect())
                }
            };

            let idx = entries.len();

            // by_hgnc_id: duplicate is an error
            if by_hgnc_id.contains_key(&hgnc_id) {
                return Err(Error::Validation(format!("duplicate HGNC ID: {hgnc_id}")));
            }
            by_hgnc_id.insert(hgnc_id, idx);

            // by_ensembl_id: last-write-wins
            if let Some(ref eid) = ensembl_gene_id {
                by_ensembl_id.insert(eid.clone(), idx);
            }

            // by_ncbi_id: last-write-wins
            if let Some(ref nid) = entrez_id {
                by_ncbi_id.insert(nid.clone(), idx);
            }

            entries.push(HgncEntry {
                hgnc_id,
                symbol,
                entrez_id,
                ensembl_gene_id,
                refseq_accession,
            });
        }

        Ok(HgncDatabase {
            entries,
            by_hgnc_id,
            by_ensembl_id,
            by_ncbi_id,
        })
    }

    /// Look up an entry by numeric HGNC ID.
    #[must_use]
    pub fn get_by_hgnc_id(&self, hgnc_id: u32) -> Option<&HgncEntry> {
        self.by_hgnc_id.get(&hgnc_id).map(|&idx| &self.entries[idx])
    }

    /// Look up an entry by Ensembl gene ID (e.g. `ENSG00000012048`).
    #[must_use]
    pub fn get_by_ensembl_id(&self, ensembl_id: &str) -> Option<&HgncEntry> {
        self.by_ensembl_id
            .get(ensembl_id)
            .map(|&idx| &self.entries[idx])
    }

    /// Look up an entry by NCBI Entrez gene ID (e.g. `672`).
    #[must_use]
    pub fn get_by_ncbi_id(&self, ncbi_id: &str) -> Option<&HgncEntry> {
        self.by_ncbi_id.get(ncbi_id).map(|&idx| &self.entries[idx])
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Return all entries sorted by HGNC ID ascending.
    #[must_use]
    pub fn sorted_entries(&self) -> Vec<&HgncEntry> {
        let mut sorted: Vec<&HgncEntry> = self.entries.iter().collect();
        sorted.sort_by_key(|e| e.hgnc_id);
        sorted
    }
}

fn find_column(columns: &[&str], name: &str) -> Result<usize, Error> {
    columns
        .iter()
        .position(|&c| c == name)
        .ok_or_else(|| Error::Parse(format!("HGNC TSV: column '{name}' not found in header")))
}

fn non_empty(s: &str) -> Option<String> {
    if s.is_empty() {
        None
    } else {
        Some(s.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_tsv(rows: &[&str]) -> String {
        let header = "hgnc_id\tsymbol\tname\tlocus_group\tlocus_type\tstatus\tlocation\tlocation_sortable\talias_symbol\talias_name\tprev_symbol\tprev_name\tgene_family\tgene_family_id\tdate_approved_reserved\tdate_symbol_changed\tdate_name_changed\tdate_modified\tentrez_id\tensembl_gene_id\tvega_id\tucsc_id\tena\trefseq_accession";
        let mut result = header.to_string();
        for row in rows {
            result.push('\n');
            result.push_str(row);
        }
        result
    }

    // Build a row with 24 columns (0-23), filling unused cols with empty strings.
    fn make_row(
        hgnc_id: &str,
        symbol: &str,
        entrez_id: &str,
        ensembl_id: &str,
        refseq: &str,
    ) -> String {
        let mut cols = vec![""; 24];
        cols[0] = hgnc_id;
        cols[1] = symbol;
        cols[18] = entrez_id;
        cols[19] = ensembl_id;
        cols[23] = refseq;
        cols.join("\t")
    }

    #[test]
    fn parse_basic() {
        let row1 = make_row("HGNC:1100", "BRCA1", "672", "ENSG00000012048", "NM_007294");
        let row2 = make_row("HGNC:5678", "TP53", "7157", "ENSG00000141510", "");
        let data = make_tsv(&[&row1, &row2]);

        let db = HgncDatabase::from_tsv(data.as_bytes()).unwrap();
        assert_eq!(db.len(), 2);

        let brca1 = db.get_by_hgnc_id(1100).unwrap();
        assert_eq!(brca1.symbol, "BRCA1");
        assert_eq!(brca1.entrez_id.as_deref(), Some("672"));
        assert_eq!(brca1.ensembl_gene_id.as_deref(), Some("ENSG00000012048"));
        assert_eq!(
            brca1.refseq_accession.as_deref(),
            Some(vec!["NM_007294".to_string()].as_slice())
        );

        let tp53 = db.get_by_ensembl_id("ENSG00000141510").unwrap();
        assert_eq!(tp53.hgnc_id, 5678);
        assert!(tp53.refseq_accession.is_none());

        let by_ncbi = db.get_by_ncbi_id("672").unwrap();
        assert_eq!(by_ncbi.hgnc_id, 1100);
    }

    #[test]
    fn pipe_delimited_refseq() {
        let row = make_row("HGNC:100", "FOO", "", "", "NM_001|NM_002|NM_003");
        let data = make_tsv(&[&row]);
        let db = HgncDatabase::from_tsv(data.as_bytes()).unwrap();
        let entry = db.get_by_hgnc_id(100).unwrap();
        assert_eq!(
            entry.refseq_accession.as_ref().unwrap(),
            &["NM_001", "NM_002", "NM_003"]
        );
    }

    #[test]
    fn empty_optional_fields_skipped() {
        let row = make_row("HGNC:200", "BAR", "", "", "");
        let data = make_tsv(&[&row]);
        let db = HgncDatabase::from_tsv(data.as_bytes()).unwrap();
        let entry = db.get_by_hgnc_id(200).unwrap();
        assert!(entry.entrez_id.is_none());
        assert!(entry.ensembl_gene_id.is_none());
        assert!(entry.refseq_accession.is_none());
    }

    #[test]
    fn duplicate_hgnc_id_errors() {
        let row1 = make_row("HGNC:100", "FOO", "", "", "");
        let row2 = make_row("HGNC:100", "BAR", "", "", "");
        let data = make_tsv(&[&row1, &row2]);
        assert!(HgncDatabase::from_tsv(data.as_bytes()).is_err());
    }

    #[test]
    fn last_write_wins_for_ensembl_id() {
        let row1 = make_row("HGNC:100", "FIRST", "", "ENSG00000000001", "");
        let row2 = make_row("HGNC:200", "SECOND", "", "ENSG00000000001", "");
        let data = make_tsv(&[&row1, &row2]);
        let db = HgncDatabase::from_tsv(data.as_bytes()).unwrap();
        // Last write wins: HGNC:200 should be indexed
        let entry = db.get_by_ensembl_id("ENSG00000000001").unwrap();
        assert_eq!(entry.hgnc_id, 200);
    }

    #[test]
    fn sorted_entries() {
        let row1 = make_row("HGNC:300", "GAMMA", "", "", "");
        let row2 = make_row("HGNC:100", "ALPHA", "", "", "");
        let row3 = make_row("HGNC:200", "BETA", "", "", "");
        let data = make_tsv(&[&row1, &row2, &row3]);
        let db = HgncDatabase::from_tsv(data.as_bytes()).unwrap();
        let sorted = db.sorted_entries();
        assert_eq!(sorted[0].hgnc_id, 100);
        assert_eq!(sorted[1].hgnc_id, 200);
        assert_eq!(sorted[2].hgnc_id, 300);
    }

    #[test]
    fn strips_hgnc_prefix() {
        let row = make_row("HGNC:10751", "SELENOP", "", "ENSG00000250722", "");
        let data = make_tsv(&[&row]);
        let db = HgncDatabase::from_tsv(data.as_bytes()).unwrap();
        assert_eq!(db.get_by_hgnc_id(10751).unwrap().symbol, "SELENOP");
    }

    #[test]
    fn numeric_hgnc_id_without_prefix() {
        let row = make_row("10751", "SELENOP", "", "ENSG00000250722", "");
        let data = make_tsv(&[&row]);
        let db = HgncDatabase::from_tsv(data.as_bytes()).unwrap();
        assert_eq!(db.get_by_hgnc_id(10751).unwrap().symbol, "SELENOP");
    }
}
