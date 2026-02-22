//! Semi-global diagonal-only DP alignment for protein sequences.

const MATCH_SCORE: i32 = 15;
const MISMATCH_SCORE: i32 = -7;
const SCORE_THRESHOLD: f64 = 10.0;

/// Result of a successful semi-global alignment.
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    pub start: usize,
    pub score: i32,
}

/// Perform semi-global diagonal-only alignment of query against reference.
///
/// Returns the best alignment position if the score meets the threshold
/// (score / query_len >= 10.0).
pub fn semi_global_align(reference: &[u8], query: &[u8]) -> Option<AlignmentResult> {
    if reference.is_empty() || query.is_empty() {
        return None;
    }

    let ref_len = reference.len();
    let query_len = query.len();

    // DP matrix: F[i][j] = F[i-1][j-1] + score(ref[i-1], query[j-1])
    // Only diagonal transitions — we can use O(query_len) space with rolling diagonals.

    // prev_diag[j] stores F[i-1][j-1] effectively
    // We iterate through reference positions and maintain the diagonal scores.

    // For each starting position in reference, compute the diagonal score.
    // F[i][j] = F[i-1][j-1] + score(ref[i-1], query[j-1])
    // F[i][0] = 0, F[0][j] = 0

    // We need to find max over all i of F[i][query_len].
    // F[i][j] = sum of scores along the diagonal from (i-j+1, 1) to (i, j)
    //         = sum_{k=0}^{j-1} score(ref[i-j+k], query[k])

    // For each diagonal d = i - j (the starting ref position when j=0),
    // compute the cumulative score for positions j=1..query_len.

    let mut max_score = i32::MIN;
    let mut max_end: usize = 0;

    // Iterate over all possible alignment start positions in reference
    for d in 0..ref_len {
        let mut score: i32 = 0;
        let available = (ref_len - d).min(query_len);

        for j in 0..available {
            score += if reference[d + j] == query[j] {
                MATCH_SCORE
            } else {
                MISMATCH_SCORE
            };
        }

        // Only valid if we covered the entire query
        if available == query_len && score > max_score {
            max_score = score;
            max_end = d + query_len;
        }
    }

    if max_score == i32::MIN {
        return None;
    }

    let score_per_base = max_score as f64 / query_len as f64;
    if score_per_base < SCORE_THRESHOLD {
        return None;
    }

    Some(AlignmentResult {
        start: max_end - query_len,
        score: max_score,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_match() {
        let reference = b"ABCDEFGHIJ";
        let query = b"DEFGH";
        let result = semi_global_align(reference, query).unwrap();
        assert_eq!(result.start, 3);
        assert_eq!(result.score, 5 * MATCH_SCORE);
    }

    #[test]
    fn perfect_match_at_start() {
        let reference = b"ABCDE";
        let query = b"ABCDE";
        let result = semi_global_align(reference, query).unwrap();
        assert_eq!(result.start, 0);
    }

    #[test]
    fn partial_alignment_with_offset() {
        let reference = b"XXXXXMPQIQKVTPVXXX";
        let query = b"MPQIQKVTPV";
        let result = semi_global_align(reference, query).unwrap();
        assert_eq!(result.start, 5);
    }

    #[test]
    fn below_threshold_returns_none() {
        // All mismatches — score per base = -7 < 10
        let reference = b"AAAAAAAAAA";
        let query = b"ZZZZZZZZZZ";
        assert!(semi_global_align(reference, query).is_none());
    }

    #[test]
    fn empty_input_returns_none() {
        assert!(semi_global_align(b"", b"ABC").is_none());
        assert!(semi_global_align(b"ABC", b"").is_none());
        assert!(semi_global_align(b"", b"").is_none());
    }

    #[test]
    fn query_longer_than_reference() {
        assert!(semi_global_align(b"AB", b"ABCDE").is_none());
    }
}
