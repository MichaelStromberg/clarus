//! Parser for NCBI genome assembly report files.

use std::collections::HashMap;
use std::io::BufRead;

use crate::chromosome::Chromosome;
use crate::error::Error;

/// Parses an NCBI genome assembly report TSV file.
///
/// Returns the chromosomes (in order of appearance) and a name-to-index lookup map.
pub fn parse_assembly_report<R: BufRead>(
    reader: R,
) -> Result<(Vec<Chromosome>, HashMap<String, usize>), Error> {
    let mut chromosomes = Vec::new();
    let mut name_to_index: HashMap<String, usize> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            return Err(Error::Parse(format!(
                "assembly report line has {} columns, expected at least 10: {line}",
                fields.len()
            )));
        }

        let ensembl_name = normalize_na(fields[0]);
        let genbank_accession = normalize_na(fields[4]);
        let refseq_accession = normalize_na(fields[6]);
        let length: u32 = fields[8]
            .parse()
            .map_err(|e| Error::Parse(format!("invalid length '{}': {e}", fields[8])))?;
        let ucsc_name = normalize_na(fields[9]);

        let ref_index = u16::try_from(chromosomes.len()).map_err(|_| {
            Error::Validation(format!(
                "too many chromosomes ({}) for u16 ref_index",
                chromosomes.len()
            ))
        })?;
        let chr = Chromosome {
            ucsc_name,
            ensembl_name,
            refseq_accession,
            genbank_accession,
            length,
            ref_index,
        };
        chr.validate()?;

        let idx = chromosomes.len();
        for name in chr.names() {
            name_to_index.insert(name.to_string(), idx);
        }
        chromosomes.push(chr);
    }

    Ok((chromosomes, name_to_index))
}

fn normalize_na(s: &str) -> String {
    if s.eq_ignore_ascii_case("na") {
        String::new()
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const SAMPLE_REPORT: &str = "\
# Assembly name:  GRCh38.p14
# comment line
1\tassembled-molecule\t1\tChromosome\tCM000663.2\t=\tNC_000001.11\tPrimary Assembly\t248956422\tchr1
2\tassembled-molecule\t2\tChromosome\tCM000664.2\t=\tNC_000002.12\tPrimary Assembly\t242193529\tchr2
HSCHR1_CTG1\tunlocalized-scaffold\t1\tChromosome\tKI270706.1\t=\tNT_187361.1\tPrimary Assembly\t175055\tna
";

    #[test]
    fn parse_basic() {
        let reader = Cursor::new(SAMPLE_REPORT);
        let (chroms, name_map) = parse_assembly_report(std::io::BufReader::new(reader)).unwrap();

        assert_eq!(chroms.len(), 3);

        // First chromosome
        assert_eq!(chroms[0].ucsc_name, "chr1");
        assert_eq!(chroms[0].ensembl_name, "1");
        assert_eq!(chroms[0].refseq_accession, "NC_000001.11");
        assert_eq!(chroms[0].genbank_accession, "CM000663.2");
        assert_eq!(chroms[0].length, 248_956_422);
        assert_eq!(chroms[0].ref_index, 0);

        // Second chromosome
        assert_eq!(chroms[1].ref_index, 1);
        assert_eq!(chroms[1].ucsc_name, "chr2");

        // Third - "na" becomes empty
        assert_eq!(chroms[2].ucsc_name, "");
        assert_eq!(chroms[2].ensembl_name, "HSCHR1_CTG1");
        assert_eq!(chroms[2].ref_index, 2);

        // Name lookups
        assert_eq!(name_map["chr1"], 0);
        assert_eq!(name_map["NC_000001.11"], 0);
        assert_eq!(name_map["1"], 0);
        assert_eq!(name_map["CM000663.2"], 0);
        assert_eq!(name_map["NC_000002.12"], 1);
        assert_eq!(name_map["NT_187361.1"], 2);
        // "na" should not be in the map
        assert!(!name_map.contains_key("na"));
    }

    #[test]
    fn skips_comments_and_empty_lines() {
        let input =
            "# comment\n\n# another\n1\tx\t1\tChromosome\tGB1\t=\tRS1\tPrimary\t100\tucsc1\n";
        let reader = Cursor::new(input);
        let (chroms, _) = parse_assembly_report(std::io::BufReader::new(reader)).unwrap();
        assert_eq!(chroms.len(), 1);
    }

    #[test]
    fn sequential_ref_index() {
        let input = "\
A\tx\t1\tC\tGB1\t=\tRS1\tP\t100\tU1
B\tx\t2\tC\tGB2\t=\tRS2\tP\t200\tU2
C\tx\t3\tC\tGB3\t=\tRS3\tP\t300\tU3
";
        let reader = Cursor::new(input);
        let (chroms, _) = parse_assembly_report(std::io::BufReader::new(reader)).unwrap();
        for (i, chr) in chroms.iter().enumerate() {
            assert_eq!(chr.ref_index, i as u16);
        }
    }
}
