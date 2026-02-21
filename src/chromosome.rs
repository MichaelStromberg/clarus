//! Chromosome representation and validation.

use crate::error::Error;

#[derive(Debug, Clone)]
pub struct Chromosome {
    pub ucsc_name: String,
    pub ensembl_name: String,
    pub refseq_accession: String,
    pub genbank_accession: String,
    pub length: u32,
    pub ref_index: u16,
}

impl Chromosome {
    pub fn validate(&self) -> Result<(), Error> {
        if self.ucsc_name.is_empty()
            && self.ensembl_name.is_empty()
            && self.refseq_accession.is_empty()
            && self.genbank_accession.is_empty()
        {
            return Err(Error::Validation(format!(
                "chromosome at ref_index {} has no name fields set",
                self.ref_index
            )));
        }
        Ok(())
    }

    /// Returns all non-empty name fields for this chromosome.
    #[must_use]
    pub fn names(&self) -> Vec<&str> {
        let mut names = Vec::with_capacity(4);
        if !self.ucsc_name.is_empty() {
            names.push(self.ucsc_name.as_str());
        }
        if !self.ensembl_name.is_empty() {
            names.push(self.ensembl_name.as_str());
        }
        if !self.refseq_accession.is_empty() {
            names.push(self.refseq_accession.as_str());
        }
        if !self.genbank_accession.is_empty() {
            names.push(self.genbank_accession.as_str());
        }
        names
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_chromosome() {
        let chr = Chromosome {
            ucsc_name: "chr1".to_string(),
            ensembl_name: "1".to_string(),
            refseq_accession: "NC_000001.11".to_string(),
            genbank_accession: "CM000663.2".to_string(),
            length: 248_956_422,
            ref_index: 0,
        };
        assert!(chr.validate().is_ok());
        assert_eq!(chr.names().len(), 4);
    }

    #[test]
    fn valid_with_partial_names() {
        let chr = Chromosome {
            ucsc_name: String::new(),
            ensembl_name: String::new(),
            refseq_accession: "NC_000001.11".to_string(),
            genbank_accession: String::new(),
            length: 100,
            ref_index: 0,
        };
        assert!(chr.validate().is_ok());
        assert_eq!(chr.names(), vec!["NC_000001.11"]);
    }

    #[test]
    fn invalid_no_names() {
        let chr = Chromosome {
            ucsc_name: String::new(),
            ensembl_name: String::new(),
            refseq_accession: String::new(),
            genbank_accession: String::new(),
            length: 100,
            ref_index: 0,
        };
        assert!(chr.validate().is_err());
    }
}
