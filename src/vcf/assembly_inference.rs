//! Assembly inference from VCF contig lengths per specification 01.
//!
//! Matches (ensembl_chrom_name, length) pairs against known GRCh37/GRCh38 tables.
//! Mitochondrial contigs are excluded from assembly inference.

use crate::genome_assembly::GenomeAssembly;

/// Contig lengths for GRCh37 (Ensembl name → length).
const GRCH37_CONTIGS: &[(&str, u32)] = &[
    ("1", 249_250_621),
    ("2", 243_199_373),
    ("3", 198_022_430),
    ("4", 191_154_276),
    ("5", 180_915_260),
    ("6", 171_115_067),
    ("7", 159_138_663),
    ("8", 146_364_022),
    ("9", 141_213_431),
    ("10", 135_534_747),
    ("11", 135_006_516),
    ("12", 133_851_895),
    ("13", 115_169_878),
    ("14", 107_349_540),
    ("15", 102_531_392),
    ("16", 90_354_753),
    ("17", 81_195_210),
    ("18", 78_077_248),
    ("19", 59_128_983),
    ("20", 63_025_520),
    ("21", 48_129_895),
    ("22", 51_304_566),
    ("X", 155_270_560),
    ("Y", 59_373_566),
];

/// Contig lengths for GRCh38 (Ensembl name → length).
const GRCH38_CONTIGS: &[(&str, u32)] = &[
    ("1", 248_956_422),
    ("2", 242_193_529),
    ("3", 198_295_559),
    ("4", 190_214_555),
    ("5", 181_538_259),
    ("6", 170_805_979),
    ("7", 159_345_973),
    ("8", 145_138_636),
    ("9", 138_394_717),
    ("10", 133_797_422),
    ("11", 135_086_622),
    ("12", 133_275_309),
    ("13", 114_364_328),
    ("14", 107_043_718),
    ("15", 101_991_189),
    ("16", 90_338_345),
    ("17", 83_257_441),
    ("18", 80_373_285),
    ("19", 58_617_616),
    ("20", 64_444_167),
    ("21", 46_709_983),
    ("22", 50_818_468),
    ("X", 156_040_895),
    ("Y", 57_227_415),
];

/// UCSC-to-Ensembl name mapping for standard chromosomes.
fn ucsc_to_ensembl(name: &str) -> Option<&str> {
    if let Some(num) = name.strip_prefix("chr") {
        if num == "M" || num == "MT" {
            Some("MT")
        } else {
            Some(num)
        }
    } else {
        None
    }
}

/// Normalize a contig name to Ensembl format.
/// Accepts: "chr1" → "1", "1" → "1", "chrM" → "MT", "MT" → "MT".
pub fn normalize_contig_name(name: &str) -> &str {
    ucsc_to_ensembl(name).unwrap_or(name)
}

/// Infer genome assembly from a single (ensembl_name, length) pair.
/// Returns `Unknown` for mitochondrial or unrecognized contigs.
pub fn infer_assembly_from_contig(ensembl_name: &str, length: u32) -> GenomeAssembly {
    // Mitochondrial contigs don't contribute to assembly inference
    if ensembl_name == "MT" {
        return GenomeAssembly::Unknown;
    }

    let matches_37 = GRCH37_CONTIGS
        .iter()
        .any(|(name, len)| *name == ensembl_name && *len == length);
    let matches_38 = GRCH38_CONTIGS
        .iter()
        .any(|(name, len)| *name == ensembl_name && *len == length);

    match (matches_37, matches_38) {
        (true, false) => GenomeAssembly::GRCh37,
        (false, true) => GenomeAssembly::GRCh38,
        (true, true) => GenomeAssembly::Unknown, // ambiguous (shouldn't happen for standard chroms)
        (false, false) => GenomeAssembly::Unknown, // unrecognized contig
    }
}

/// Infer assembly from a collection of (ensembl_name, length) pairs.
/// All recognized contigs must agree. Returns `Unknown` if no contigs matched
/// or there is a conflict.
pub fn infer_assembly(contigs: &[(String, u32)]) -> GenomeAssembly {
    let mut result = GenomeAssembly::Unknown;

    for (name, length) in contigs {
        let assembly = infer_assembly_from_contig(name, *length);
        if assembly == GenomeAssembly::Unknown {
            continue;
        }
        if result == GenomeAssembly::Unknown {
            result = assembly;
        } else if result != assembly {
            // Conflict: contigs disagree
            return GenomeAssembly::Unknown;
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grch38_chr1() {
        assert_eq!(
            infer_assembly_from_contig("1", 248_956_422),
            GenomeAssembly::GRCh38
        );
    }

    #[test]
    fn grch37_chr1() {
        assert_eq!(
            infer_assembly_from_contig("1", 249_250_621),
            GenomeAssembly::GRCh37
        );
    }

    #[test]
    fn mitochondrial_excluded() {
        assert_eq!(
            infer_assembly_from_contig("MT", 16569),
            GenomeAssembly::Unknown
        );
    }

    #[test]
    fn unrecognized_contig() {
        assert_eq!(
            infer_assembly_from_contig("1", 999),
            GenomeAssembly::Unknown
        );
    }

    #[test]
    fn infer_assembly_multiple_contigs() {
        let contigs = vec![
            ("1".to_string(), 248_956_422),
            ("2".to_string(), 242_193_529),
            ("MT".to_string(), 16569),
        ];
        assert_eq!(infer_assembly(&contigs), GenomeAssembly::GRCh38);
    }

    #[test]
    fn infer_assembly_conflict() {
        let contigs = vec![
            ("1".to_string(), 248_956_422), // GRCh38
            ("2".to_string(), 243_199_373), // GRCh37
        ];
        assert_eq!(infer_assembly(&contigs), GenomeAssembly::Unknown);
    }

    #[test]
    fn normalize_ucsc_names() {
        assert_eq!(normalize_contig_name("chr1"), "1");
        assert_eq!(normalize_contig_name("chrX"), "X");
        assert_eq!(normalize_contig_name("chrM"), "MT");
        assert_eq!(normalize_contig_name("chrMT"), "MT");
        assert_eq!(normalize_contig_name("1"), "1");
        assert_eq!(normalize_contig_name("MT"), "MT");
    }
}
