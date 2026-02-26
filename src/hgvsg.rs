//! HGVSg (genomic HGVS) notation per specification 09.
//!
//! Format: `{refseq_accession}:g.{change}`
//! Applies right-rotation (forward strand) before formatting.

/// Maximum downstream bases to fetch for right-rotation.
const MAX_DOWNSTREAM: usize = 500;

/// Result of HGVSg computation.
#[derive(Debug, Clone)]
pub struct HgvsgResult {
    /// The full HGVSg string, e.g., "NC_000001.11:g.1000A>G".
    pub notation: String,
    /// Position after right-rotation (for JSON output).
    pub rotated_start: i64,
    /// End position after right-rotation.
    pub rotated_end: i64,
}

/// Positions returned by the buffer-based HGVSg computation.
#[derive(Debug, Clone)]
pub struct HgvsgPositions {
    pub rotated_start: i64,
    pub rotated_end: i64,
}

/// Compute HGVSg notation directly into a reusable buffer.
///
/// Eliminates intermediate `format!()` allocations by writing directly via `push_str`
/// and `itoa`. Returns only the rotated positions; the notation is in `buf`.
pub fn compute_hgvsg_into<'a>(
    buf: &mut String,
    refseq_accession: &str,
    position: i64,
    ref_allele: &str,
    alt_allele: &str,
    ref_accessor: &dyn Fn(i64, usize) -> Option<&'a [u8]>,
) -> HgvsgPositions {
    buf.clear();
    let mut itoa_buf = itoa::Buffer::new();
    let ref_len = ref_allele.len();
    let alt_len = alt_allele.len();

    // Identity: ref == alt
    if ref_allele == alt_allele {
        buf.push_str(refseq_accession);
        buf.push_str(":g.");
        buf.push_str(itoa_buf.format(position));
        buf.push('=');
        return HgvsgPositions {
            rotated_start: position,
            rotated_end: position + ref_len as i64 - 1,
        };
    }

    // Substitution: single base change
    if ref_len == 1 && alt_len == 1 {
        buf.push_str(refseq_accession);
        buf.push_str(":g.");
        buf.push_str(itoa_buf.format(position));
        buf.push_str(ref_allele);
        buf.push('>');
        buf.push_str(alt_allele);
        return HgvsgPositions {
            rotated_start: position,
            rotated_end: position,
        };
    }

    // Pure deletion (alt empty)
    if alt_len == 0 && ref_len > 0 {
        return compute_deletion_hgvsg_into(
            buf,
            refseq_accession,
            position,
            ref_allele,
            ref_accessor,
        );
    }

    // Pure insertion (ref empty)
    if ref_len == 0 && alt_len > 0 {
        return compute_insertion_hgvsg_into(
            buf,
            refseq_accession,
            position,
            alt_allele,
            ref_accessor,
        );
    }

    // Inversion check
    if ref_len == alt_len && ref_len > 1 && is_reverse_complement(ref_allele, alt_allele) {
        let end = position + ref_len as i64 - 1;
        buf.push_str(refseq_accession);
        buf.push_str(":g.");
        buf.push_str(itoa_buf.format(position));
        buf.push('_');
        buf.push_str(itoa_buf.format(end));
        buf.push_str("inv");
        return HgvsgPositions {
            rotated_start: position,
            rotated_end: end,
        };
    }

    // DelIns
    let end = position + ref_len as i64 - 1;
    buf.push_str(refseq_accession);
    buf.push_str(":g.");
    if ref_len == 1 {
        buf.push_str(itoa_buf.format(position));
        buf.push_str("delins");
        buf.push_str(alt_allele);
    } else {
        buf.push_str(itoa_buf.format(position));
        buf.push('_');
        buf.push_str(itoa_buf.format(end));
        buf.push_str("delins");
        buf.push_str(alt_allele);
    }
    HgvsgPositions {
        rotated_start: position,
        rotated_end: end,
    }
}

/// Compute HGVSg notation for a variant.
///
/// Convenience wrapper around [`compute_hgvsg_into`] that returns an owned `HgvsgResult`.
/// Prefer `compute_hgvsg_into` in hot paths to reuse a buffer across calls.
pub fn compute_hgvsg<'a>(
    refseq_accession: &str,
    position: i64,
    ref_allele: &str,
    alt_allele: &str,
    ref_accessor: &dyn Fn(i64, usize) -> Option<&'a [u8]>,
) -> HgvsgResult {
    let mut buf = String::with_capacity(128);
    let positions = compute_hgvsg_into(
        &mut buf,
        refseq_accession,
        position,
        ref_allele,
        alt_allele,
        ref_accessor,
    );
    HgvsgResult {
        notation: buf,
        rotated_start: positions.rotated_start,
        rotated_end: positions.rotated_end,
    }
}

/// Compute HGVSg for a deletion into a buffer, applying right-rotation.
fn compute_deletion_hgvsg_into<'a>(
    buf: &mut String,
    refseq_accession: &str,
    position: i64,
    ref_allele: &str,
    ref_accessor: &dyn Fn(i64, usize) -> Option<&'a [u8]>,
) -> HgvsgPositions {
    let ref_len = ref_allele.len();
    let downstream_start_0based = position + ref_len as i64 - 1;
    let (rotated_pos, _rotated_allele) = right_rotate(
        position,
        downstream_start_0based,
        ref_allele.as_bytes(),
        ref_accessor,
    );

    let start = rotated_pos;
    let end = rotated_pos + ref_len as i64 - 1;
    let mut itoa_buf = itoa::Buffer::new();

    buf.push_str(refseq_accession);
    buf.push_str(":g.");
    if ref_len == 1 {
        buf.push_str(itoa_buf.format(start));
        buf.push_str("del");
    } else {
        buf.push_str(itoa_buf.format(start));
        buf.push('_');
        buf.push_str(itoa_buf.format(end));
        buf.push_str("del");
    }

    HgvsgPositions {
        rotated_start: start,
        rotated_end: end,
    }
}

/// Compute HGVSg for an insertion into a buffer, applying right-rotation and duplication detection.
fn compute_insertion_hgvsg_into<'a>(
    buf: &mut String,
    refseq_accession: &str,
    position: i64,
    alt_allele: &str,
    ref_accessor: &dyn Fn(i64, usize) -> Option<&'a [u8]>,
) -> HgvsgPositions {
    let alt_bytes = alt_allele.as_bytes();
    let alt_len = alt_bytes.len();

    let downstream_start_0based = position - 1;
    let (rotated_pos, rotated_allele) =
        right_rotate(position, downstream_start_0based, alt_bytes, ref_accessor);

    let dup_check_start_0based = rotated_pos - 1 - alt_len as i64;
    let is_dup = if dup_check_start_0based >= 0 {
        ref_accessor(dup_check_start_0based, alt_len)
            .is_some_and(|upstream| upstream == rotated_allele)
    } else {
        false
    };

    let mut itoa_buf = itoa::Buffer::new();

    if is_dup {
        let dup_start = rotated_pos - alt_len as i64;
        let dup_end = rotated_pos - 1;
        buf.push_str(refseq_accession);
        buf.push_str(":g.");
        if alt_len == 1 {
            buf.push_str(itoa_buf.format(dup_start));
            buf.push_str("dup");
        } else {
            buf.push_str(itoa_buf.format(dup_start));
            buf.push('_');
            buf.push_str(itoa_buf.format(dup_end));
            buf.push_str("dup");
        }
        HgvsgPositions {
            rotated_start: dup_start,
            rotated_end: dup_end,
        }
    } else {
        let ins_start = rotated_pos - 1;
        let ins_end = rotated_pos;
        let rotated_str = String::from_utf8_lossy(&rotated_allele);
        buf.push_str(refseq_accession);
        buf.push_str(":g.");
        buf.push_str(itoa_buf.format(ins_start));
        buf.push('_');
        buf.push_str(itoa_buf.format(ins_end));
        buf.push_str("ins");
        buf.push_str(&rotated_str);
        HgvsgPositions {
            rotated_start: ins_start,
            rotated_end: ins_end,
        }
    }
}

/// Right-rotate (shift downstream) an insertion or deletion allele.
///
/// `position_1based` is the 1-based variant position (returned position will also be 1-based).
/// `downstream_start_0based` is the 0-based offset where downstream reference begins.
/// Returns (new_position_1based, rotated_allele_bytes).
fn right_rotate<'a>(
    position_1based: i64,
    downstream_start_0based: i64,
    allele: &[u8],
    ref_accessor: &dyn Fn(i64, usize) -> Option<&'a [u8]>,
) -> (i64, Vec<u8>) {
    let allele_len = allele.len();
    if allele_len == 0 {
        return (position_1based, allele.to_vec());
    }

    let downstream_len = MAX_DOWNSTREAM.max(allele_len);
    let downstream = match ref_accessor(downstream_start_0based, downstream_len) {
        Some(seq) => seq,
        None => return (position_1based, allele.to_vec()),
    };

    // Scan forward: compare virtual context[i] to context[i + allele_len]
    // where context = allele ++ downstream (without actually concatenating).
    let mut shift = 0usize;
    let max_shift = downstream.len();

    while shift < max_shift {
        // context[shift] is allele[shift] when shift < allele_len, else downstream[shift - allele_len]
        let left = if shift < allele_len {
            allele[shift]
        } else {
            downstream[shift - allele_len]
        };
        // context[shift + allele_len] is always downstream[shift]
        if left == downstream[shift] {
            shift += 1;
        } else {
            break;
        }
    }

    if shift == 0 {
        return (position_1based, allele.to_vec());
    }

    let new_pos = position_1based + shift as i64;

    // Extract rotated allele from virtual context[shift..shift+allele_len]
    let mut rotated = Vec::with_capacity(allele_len);
    for i in shift..shift + allele_len {
        if i < allele_len {
            rotated.push(allele[i]);
        } else {
            rotated.push(downstream[i - allele_len]);
        }
    }
    (new_pos, rotated)
}

/// Check if `alt` is the reverse complement of `ref_allele`.
fn is_reverse_complement(ref_allele: &str, alt_allele: &str) -> bool {
    if ref_allele.len() != alt_allele.len() {
        return false;
    }
    let ref_bytes = ref_allele.as_bytes();
    let alt_bytes = alt_allele.as_bytes();

    ref_bytes
        .iter()
        .zip(alt_bytes.iter().rev())
        .all(|(r, a)| complement(*r) == Some(*a))
}

fn complement(base: u8) -> Option<u8> {
    match base {
        b'A' => Some(b'T'),
        b'T' => Some(b'A'),
        b'C' => Some(b'G'),
        b'G' => Some(b'C'),
        b'a' => Some(b't'),
        b't' => Some(b'a'),
        b'c' => Some(b'g'),
        b'g' => Some(b'c'),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_accessor<'a>(sequence: &'a [u8]) -> impl Fn(i64, usize) -> Option<&'a [u8]> + 'a {
        move |offset: i64, len: usize| {
            if offset < 0 {
                return None;
            }
            let start = offset as usize;
            if start >= sequence.len() {
                return None;
            }
            let end = (start + len).min(sequence.len());
            Some(&sequence[start..end])
        }
    }

    #[test]
    fn substitution() {
        let accessor = make_accessor(b"ACGTACGT");
        let result = compute_hgvsg("NC_000001.11", 5, "A", "G", &accessor);
        assert_eq!(result.notation, "NC_000001.11:g.5A>G");
    }

    #[test]
    fn identity() {
        let accessor = make_accessor(b"ACGTACGT");
        let result = compute_hgvsg("NC_000001.11", 5, "A", "A", &accessor);
        assert_eq!(result.notation, "NC_000001.11:g.5=");
    }

    #[test]
    fn single_base_deletion() {
        // Deletion of A at position 5 in reference ACGTACGT
        let accessor = make_accessor(b"ACGTACGTNN");
        let result = compute_hgvsg("NC_000001.11", 5, "A", "", &accessor);
        // Should right-rotate: A at pos 5 matches downstream
        assert!(result.notation.contains("del"));
    }

    #[test]
    fn multi_base_deletion() {
        let accessor = make_accessor(b"ACGTXYZACGT");
        let result = compute_hgvsg("NC_000001.11", 5, "XYZ", "", &accessor);
        assert!(result.notation.contains("del"));
        assert!(result.notation.contains('_'));
    }

    #[test]
    fn insertion_becomes_duplication() {
        // Insert "A" at position 2 in reference "AACGT"
        // Upstream base at position 1 is "A" — matches inserted "A" → duplication
        let accessor = make_accessor(b"AACGTNN");
        let result = compute_hgvsg("NC_000001.11", 2, "", "A", &accessor);
        assert!(
            result.notation.contains("dup"),
            "expected dup, got: {}",
            result.notation
        );
    }

    #[test]
    fn insertion_not_duplication() {
        // Insert "T" at position 2 in reference "AACGT"
        // Upstream base at position 1 is "A" — doesn't match "T"
        let accessor = make_accessor(b"AACGTNN");
        let result = compute_hgvsg("NC_000001.11", 2, "", "T", &accessor);
        assert!(
            result.notation.contains("ins"),
            "expected ins, got: {}",
            result.notation
        );
    }

    #[test]
    fn inversion() {
        let accessor = make_accessor(b"ACGTACGT");
        let _result = compute_hgvsg("NC_000001.11", 1, "ACGT", "ACGT", &accessor);
        // ACGT reversed complement is ACGT — palindromic. Use a non-palindromic sequence.
        let result = compute_hgvsg("NC_000001.11", 1, "AACG", "CGTT", &accessor);
        assert!(
            result.notation.contains("inv"),
            "expected inv, got: {}",
            result.notation
        );
    }

    #[test]
    fn delins_single() {
        let accessor = make_accessor(b"ACGTACGT");
        let result = compute_hgvsg("NC_000001.11", 5, "A", "TG", &accessor);
        assert_eq!(result.notation, "NC_000001.11:g.5delinsTG");
    }

    #[test]
    fn delins_range() {
        let accessor = make_accessor(b"ACGTACGT");
        let result = compute_hgvsg("NC_000001.11", 5, "ACG", "TG", &accessor);
        assert_eq!(result.notation, "NC_000001.11:g.5_7delinsTG");
    }

    #[test]
    fn reverse_complement_check() {
        assert!(is_reverse_complement("AACG", "CGTT"));
        assert!(is_reverse_complement("AT", "AT"));
        assert!(!is_reverse_complement("AA", "TG"));
        assert!(!is_reverse_complement("A", "TT"));
    }
}
