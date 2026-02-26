//! VCF INFO field parser per specification 01.

/// Parsed INFO field values.
#[derive(Debug, Default)]
pub struct InfoFields {
    pub end: Option<i64>,
    pub svtype: Option<String>,
    pub svlen: Option<i64>,
    pub cipos: Option<(i64, i64)>,
    pub ciend: Option<(i64, i64)>,
    pub af: Option<Vec<f64>>,
    pub raf: Option<Vec<f64>>,
    pub event: Option<String>,
    pub fisher_strand_bias: Option<f64>,
    pub mapping_quality: Option<f64>,
    pub strand_bias: Option<f64>,
    pub imprecise: bool,
    pub impute_score: Option<f64>,
    pub repeat_unit: Option<String>,
    pub ref_repeat_count: Option<i32>,
    /// INV3/INV5 flags (always false from VCF, included for data model completeness).
    pub inv3: bool,
    pub inv5: bool,
}

/// Parse a VCF INFO field string (the 8th column).
pub fn parse_info(info: &str) -> InfoFields {
    let mut fields = InfoFields::default();

    if info == "." || info.is_empty() {
        return fields;
    }

    for token in info.split(';') {
        if let Some((key, value)) = token.split_once('=') {
            match key {
                "END" => fields.end = value.parse().ok(),
                "SVTYPE" => fields.svtype = Some(value.to_string()),
                "SVLEN" => {
                    // Take absolute value per spec
                    fields.svlen = value.parse::<i64>().ok().map(|v| v.abs());
                }
                "CIPOS" => fields.cipos = parse_ci_pair(value),
                "CIEND" => fields.ciend = parse_ci_pair(value),
                "AF" => fields.af = parse_float_array(value),
                "RAF" => fields.raf = parse_float_array(value),
                "EVENT" => fields.event = Some(value.to_string()),
                "FS" => fields.fisher_strand_bias = value.parse().ok(),
                "MQ" => fields.mapping_quality = value.parse().ok(),
                "SB" => fields.strand_bias = value.parse().ok(),
                "INFO" => fields.impute_score = value.parse().ok(),
                "RU" => fields.repeat_unit = Some(value.to_string()),
                "REF" => fields.ref_repeat_count = value.parse().ok(),
                _ => {}
            }
        } else {
            // Flag field (no '=')
            if token == "IMPRECISE" {
                fields.imprecise = true;
            }
        }
    }

    fields
}

fn parse_ci_pair(value: &str) -> Option<(i64, i64)> {
    let parts: Vec<&str> = value.split(',').collect();
    if parts.len() == 2 {
        let a = parts[0].parse().ok()?;
        let b = parts[1].parse().ok()?;
        Some((a, b))
    } else {
        None
    }
}

fn parse_float_array(value: &str) -> Option<Vec<f64>> {
    let values: Result<Vec<f64>, _> = value.split(',').map(|v| v.parse()).collect();
    values.ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_simple_info() {
        let info = parse_info("END=1000;SVTYPE=DEL;SVLEN=-500;IMPRECISE");
        assert_eq!(info.end, Some(1000));
        assert_eq!(info.svtype.as_deref(), Some("DEL"));
        assert_eq!(info.svlen, Some(500)); // absolute value
        assert!(info.imprecise);
    }

    #[test]
    fn parse_cipos() {
        let info = parse_info("CIPOS=-10,10;CIEND=-5,5");
        assert_eq!(info.cipos, Some((-10, 10)));
        assert_eq!(info.ciend, Some((-5, 5)));
    }

    #[test]
    fn parse_af() {
        let info = parse_info("AF=0.5,0.3");
        let af = info.af.unwrap();
        assert_eq!(af.len(), 2);
        assert!((af[0] - 0.5).abs() < f64::EPSILON);
    }

    #[test]
    fn dot_info() {
        let info = parse_info(".");
        assert_eq!(info.end, None);
        assert_eq!(info.svtype, None);
        assert!(!info.imprecise);
    }

    #[test]
    fn parse_repeat_fields() {
        let info = parse_info("RU=CAG;REF=20");
        assert_eq!(info.repeat_unit.as_deref(), Some("CAG"));
        assert_eq!(info.ref_repeat_count, Some(20));
    }
}
