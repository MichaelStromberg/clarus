//! Parser for NCBI cytogenetic band (ideogram) files.

use std::collections::HashMap;
use std::io::BufRead;

use crate::error::Error;

/// A cytogenetic band on a chromosome.
#[derive(Debug, Clone, PartialEq)]
pub struct Band {
    /// 1-based inclusive start position.
    pub begin: u32,
    /// 1-based inclusive end position.
    pub end: u32,
    /// Band name, e.g. "p36.33".
    pub name: String,
}

/// Parses an NCBI GDP ideogram file into a map of ensembl chromosome name to bands.
///
/// Expected format: tab-separated with columns:
///   chromosome, arm, band, iscn_start, iscn_stop, bp_start, bp_stop, stain[, density]
///
/// Lines starting with `#` are skipped.
pub fn parse_ideogram<R: BufRead>(reader: R) -> Result<HashMap<String, Vec<Band>>, Error> {
    let mut bands: HashMap<String, Vec<Band>> = HashMap::new();

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result?;

        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 7 {
            return Err(Error::Parse(format!(
                "line {}: expected at least 7 tab-separated columns, got {}",
                line_num + 1,
                fields.len()
            )));
        }

        let chromosome = fields[0].to_string();
        let arm = fields[1];
        let band_num = fields[2];

        let bp_start: u32 = fields[5].trim().parse().map_err(|e| {
            Error::Parse(format!(
                "line {}: invalid bp_start '{}': {e}",
                line_num + 1,
                fields[5]
            ))
        })?;

        let bp_stop: u32 = fields[6].trim().parse().map_err(|e| {
            Error::Parse(format!(
                "line {}: invalid bp_stop '{}': {e}",
                line_num + 1,
                fields[6]
            ))
        })?;

        let name = format!("{arm}{band_num}");

        bands.entry(chromosome).or_default().push(Band {
            begin: bp_start,
            end: bp_stop,
            name,
        });
    }

    Ok(bands)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const SAMPLE_IDEOGRAM: &str = "\
#chromosome\tarm\tband\tiscn_start\tiscn_stop\tbp_start\tbp_stop\tstain\tdensity
1\tp\t36.33\t0\t100\t1\t2300000\tgneg\t
1\tp\t36.32\t100\t244\t2300001\t5300000\tgpos\t25
1\tq\t44\t5000\t5100\t243500001\t248956422\tacen\t
X\tp\t22.33\t0\t50\t1\t2800000\tgneg\t
";

    #[test]
    fn parse_band_fields() {
        let reader = Cursor::new(SAMPLE_IDEOGRAM);
        let result = parse_ideogram(reader).unwrap();

        let chr1_bands = &result["1"];
        assert_eq!(chr1_bands[0].begin, 1);
        assert_eq!(chr1_bands[0].end, 2_300_000);
        assert_eq!(chr1_bands[0].name, "p36.33");

        assert_eq!(chr1_bands[1].begin, 2_300_001);
        assert_eq!(chr1_bands[1].end, 5_300_000);
        assert_eq!(chr1_bands[1].name, "p36.32");
    }

    #[test]
    fn band_name_construction() {
        let reader = Cursor::new(SAMPLE_IDEOGRAM);
        let result = parse_ideogram(reader).unwrap();

        let chr1_bands = &result["1"];
        assert_eq!(chr1_bands[2].name, "q44");

        let chrx_bands = &result["X"];
        assert_eq!(chrx_bands[0].name, "p22.33");
    }

    #[test]
    fn grouping_by_chromosome() {
        let reader = Cursor::new(SAMPLE_IDEOGRAM);
        let result = parse_ideogram(reader).unwrap();

        assert_eq!(result.len(), 2); // "1" and "X"
        assert_eq!(result["1"].len(), 3);
        assert_eq!(result["X"].len(), 1);
    }

    #[test]
    fn comment_lines_skipped() {
        let input =
            "#this is a comment\n#another comment\n1\tp\t36.33\t0\t100\t1\t2300000\tgneg\t\n";
        let reader = Cursor::new(input);
        let result = parse_ideogram(reader).unwrap();

        assert_eq!(result.len(), 1);
        assert_eq!(result["1"].len(), 1);
    }

    #[test]
    fn empty_lines_skipped() {
        let input = "\n\n1\tp\t36.33\t0\t100\t1\t2300000\tgneg\t\n\n";
        let reader = Cursor::new(input);
        let result = parse_ideogram(reader).unwrap();

        assert_eq!(result.len(), 1);
    }

    #[test]
    fn too_few_columns_error() {
        let input = "1\tp\t36.33\t0\t100\t1\n";
        let reader = Cursor::new(input);
        let result = parse_ideogram(reader);

        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("expected at least 7")
        );
    }
}
