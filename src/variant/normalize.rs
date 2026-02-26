//! Variant normalization per specification 03.
//!
//! Two stages:
//! 1. Bidirectional trimming + iterative left-alignment (VCF parsing time)
//! 2. Right-shifting for HGVS (HGVS generation time, in hgvsg.rs)

use std::borrow::Cow;

/// Result of normalization.
///
/// Uses `Cow<str>` to avoid allocation for the common case: `trim()` returns
/// borrowed substrings of the original alleles (~85% SNVs need no trim), and
/// `normalize()` only allocates when left-alignment actually shifts the variant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalizedVariant<'a> {
    pub position: i64,
    pub ref_allele: Cow<'a, str>,
    pub alt_allele: Cow<'a, str>,
}

/// Maximum upstream bases to fetch per left-alignment iteration.
const LEFT_ALIGN_WINDOW: usize = 50;

/// Bidirectional trim: remove common prefix and suffix from REF/ALT alleles.
/// Returns borrowed substrings — zero allocation.
pub fn trim<'a>(position: i64, ref_allele: &'a str, alt_allele: &'a str) -> NormalizedVariant<'a> {
    let ref_bytes = ref_allele.as_bytes();
    let alt_bytes = alt_allele.as_bytes();

    // Fast path: SNVs and identical alleles need no trimming (~85% of variants)
    if ref_bytes.len() == 1 && alt_bytes.len() == 1 || ref_bytes == alt_bytes {
        return NormalizedVariant {
            position,
            ref_allele: Cow::Borrowed(ref_allele),
            alt_allele: Cow::Borrowed(alt_allele),
        };
    }

    // Step 1: Find longest common prefix
    let prefix_len = ref_bytes
        .iter()
        .zip(alt_bytes.iter())
        .take_while(|(a, b)| a == b)
        .count();

    // Fast path: no common prefix or suffix possible (first bases differ)
    if prefix_len == 0 {
        // Check suffix only
        let suffix_len = ref_bytes
            .iter()
            .rev()
            .zip(alt_bytes.iter().rev())
            .take_while(|(a, b)| a == b)
            .count();

        if suffix_len == 0 {
            return NormalizedVariant {
                position,
                ref_allele: Cow::Borrowed(ref_allele),
                alt_allele: Cow::Borrowed(alt_allele),
            };
        }

        // Input is valid UTF-8 &str, so substrings are safe
        return NormalizedVariant {
            position,
            ref_allele: Cow::Borrowed(&ref_allele[..ref_bytes.len() - suffix_len]),
            alt_allele: Cow::Borrowed(&alt_allele[..alt_bytes.len() - suffix_len]),
        };
    }

    let new_pos = position + prefix_len as i64;

    // Step 2: Find longest common suffix (after prefix removal)
    let suffix_len = ref_bytes[prefix_len..]
        .iter()
        .rev()
        .zip(alt_bytes[prefix_len..].iter().rev())
        .take_while(|(a, b)| a == b)
        .count();

    // Input is valid UTF-8 &str, so substring indexing is safe
    NormalizedVariant {
        position: new_pos,
        ref_allele: Cow::Borrowed(&ref_allele[prefix_len..ref_bytes.len() - suffix_len]),
        alt_allele: Cow::Borrowed(&alt_allele[prefix_len..alt_bytes.len() - suffix_len]),
    }
}

/// Full normalization: iterative trim + left-alignment.
///
/// `ref_accessor` fetches reference bases at (0-based offset, length).
/// Returns the fully normalized variant. For SNVs (85% of variants), returns
/// borrowed slices with zero allocation.
pub fn normalize<'a, 'r>(
    position: i64,
    ref_allele: &'a str,
    alt_allele: &'a str,
    ref_accessor: &dyn Fn(i64, usize) -> Option<&'r [u8]>,
) -> NormalizedVariant<'a> {
    let initial = trim(position, ref_allele, alt_allele);

    // Fast path: no left-alignment needed (SNVs, MNVs, delins, symbolic)
    if !can_left_align(&initial) {
        return initial;
    }

    // Slow path: left-alignment loop requires owned data
    let mut pos = initial.position;
    let mut ref_owned = initial.ref_allele.into_owned();
    let mut alt_owned = initial.alt_allele.into_owned();

    loop {
        let prev_pos = pos;

        // Left-align
        let (new_pos, rotated) = left_align_raw(
            pos,
            if ref_owned.is_empty() {
                &alt_owned
            } else {
                &ref_owned
            },
            ref_accessor,
        );
        pos = new_pos;
        if let Some(rotated_str) = rotated {
            if ref_owned.is_empty() {
                alt_owned = rotated_str;
            } else {
                ref_owned = rotated_str;
            }
        }

        if pos == prev_pos {
            break;
        }

        // Re-trim after left-alignment (scoped to release borrows before reassignment)
        let (new_pos, new_ref, new_alt) = {
            let retrimmed = trim(pos, &ref_owned, &alt_owned);
            (
                retrimmed.position,
                retrimmed.ref_allele.into_owned(),
                retrimmed.alt_allele.into_owned(),
            )
        };
        pos = new_pos;
        ref_owned = new_ref;
        alt_owned = new_alt;
    }

    NormalizedVariant {
        position: pos,
        ref_allele: Cow::Owned(ref_owned),
        alt_allele: Cow::Owned(alt_owned),
    }
}

/// Check if a variant can be left-aligned: must be a pure insertion or deletion.
fn can_left_align(variant: &NormalizedVariant<'_>) -> bool {
    // Symbolic alleles cannot be left-aligned
    if variant.alt_allele.starts_with('<')
        || variant.alt_allele.contains('[')
        || variant.alt_allele.contains(']')
    {
        return false;
    }

    // Only pure insertions (ref empty) or pure deletions (alt empty)
    variant.ref_allele.is_empty() || variant.alt_allele.is_empty()
}

/// Left-align a pure insertion or deletion, returning (new_position, rotated_allele).
/// Returns `None` for the rotated allele when shift == 0 (no change).
fn left_align_raw<'r>(
    position: i64,
    allele: &str,
    ref_accessor: &dyn Fn(i64, usize) -> Option<&'r [u8]>,
) -> (i64, Option<String>) {
    let allele_bytes = allele.as_bytes();
    let allele_len = allele_bytes.len();

    if allele_len == 0 {
        return (position, None);
    }

    // Fetch upstream reference window. Position is 1-based, accessor is 0-based.
    let pos_0based = position - 1;
    let window_size = LEFT_ALIGN_WINDOW.min(pos_0based as usize);
    if window_size == 0 {
        return (position, None);
    }

    let upstream_start = pos_0based - window_size as i64;
    let upstream = match ref_accessor(upstream_start, window_size) {
        Some(seq) => seq,
        None => return (position, None),
    };

    // Virtual context = upstream ++ allele_bytes. Walk backward comparing
    // context[i] to context[i - allele_len] without allocating.
    let total_len = upstream.len() + allele_len;
    if total_len <= allele_len {
        return (position, None);
    }

    let ctx = |idx: usize| -> u8 {
        if idx < upstream.len() {
            upstream[idx]
        } else {
            allele_bytes[idx - upstream.len()]
        }
    };

    let mut shift = 0usize;
    let mut i = total_len - 1;

    while i >= allele_len {
        if ctx(i) == ctx(i - allele_len) {
            shift += 1;
            i -= 1;
        } else {
            break;
        }
    }

    if shift == 0 {
        return (position, None);
    }

    let new_pos = position - shift as i64;

    // Extract rotated allele at new position
    let rotated_start = upstream.len() - shift;
    let mut rotated_buf = Vec::with_capacity(allele_len);
    for idx in rotated_start..rotated_start + allele_len {
        rotated_buf.push(ctx(idx));
    }
    let rotated_str = String::from_utf8(rotated_buf)
        .unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned());

    (new_pos, Some(rotated_str))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trim_snv_no_change() {
        let result = trim(100, "A", "G");
        assert_eq!(result.position, 100);
        assert_eq!(result.ref_allele, "A");
        assert_eq!(result.alt_allele, "G");
    }

    #[test]
    fn trim_common_prefix() {
        // REF=ACG, ALT=AT → trim A prefix → pos+1, REF=CG, ALT=T
        let result = trim(100, "ACG", "AT");
        assert_eq!(result.position, 101);
        assert_eq!(result.ref_allele, "CG");
        assert_eq!(result.alt_allele, "T");
    }

    #[test]
    fn trim_common_suffix() {
        // REF=ACG, ALT=TG → trim G suffix → REF=AC, ALT=T, then trim nothing
        let result = trim(100, "ACG", "TG");
        assert_eq!(result.position, 100);
        assert_eq!(result.ref_allele, "AC");
        assert_eq!(result.alt_allele, "T");
    }

    #[test]
    fn trim_insertion_padding() {
        // REF=A, ALT=ACGT → trim A prefix → pos+1, REF="", ALT="CGT"
        let result = trim(100, "A", "ACGT");
        assert_eq!(result.position, 101);
        assert_eq!(result.ref_allele, "");
        assert_eq!(result.alt_allele, "CGT");
    }

    #[test]
    fn trim_deletion_padding() {
        // REF=ACGT, ALT=A → trim A prefix → pos+1, REF="CGT", ALT=""
        let result = trim(100, "ACGT", "A");
        assert_eq!(result.position, 101);
        assert_eq!(result.ref_allele, "CGT");
        assert_eq!(result.alt_allele, "");
    }

    #[test]
    fn trim_identical() {
        let result = trim(100, "ACG", "ACG");
        assert_eq!(result.position, 100);
        assert_eq!(result.ref_allele, "ACG");
        assert_eq!(result.alt_allele, "ACG");
    }

    #[test]
    fn left_align_deletion() {
        // Reference (0-based): A(0) A(1) A(2) C(3) G(4) T(5)
        // VCF 1-based positions: A(1) A(2) A(3) C(4) G(5) T(6)
        // Deletion: position=2, REF="AA", ALT="A" (bases at positions 2-3 are AA)
        // After trim: position=3, REF="A", ALT=""
        // Left-align: position 3 (1-based) has upstream positions 1-2 = "AA"
        // Can shift left through the run of A's to position 1
        let reference = b"AAACGT";
        let accessor = |offset: i64, len: usize| -> Option<&[u8]> {
            if offset < 0 {
                return None;
            }
            let start = offset as usize;
            let end = (start + len).min(reference.len());
            if start < reference.len() {
                Some(&reference[start..end])
            } else {
                None
            }
        };

        let result = normalize(2, "AA", "A", &accessor);
        assert_eq!(result.ref_allele, "A");
        assert_eq!(result.alt_allele, "");
        assert_eq!(result.position, 1); // Left-aligned to first A
    }

    #[test]
    fn normalize_snv_no_change() {
        let accessor = |_offset: i64, _len: usize| -> Option<&[u8]> { None };
        let result = normalize(100, "A", "G", &accessor);
        assert_eq!(result.position, 100);
        assert_eq!(result.ref_allele, "A");
        assert_eq!(result.alt_allele, "G");
    }
}
