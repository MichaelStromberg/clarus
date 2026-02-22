//! Selenoprotein gene identification by HGNC ID.

/// The 16 HGNC IDs for known selenoprotein genes.
static SELENOPROTEIN_HGNC_IDS: [i32; 16] = [
    2884,  // DIO1
    4555,  // GPX1
    4556,  // GPX2
    10751, // SELENON
    10752, // SELENOP
    14133, // SEPHS2
    15999, // SELENOW
    17705, // SELENOT
    18136, // SELENOO
    18251, // SELENOK
    29361, // SELENOF
    30394, // SELENOI
    30395, // SELENOS
    30396, // SELENOV
    30397, // SELENOH
    30399, // SELENOM
];

/// Returns true if the given HGNC ID identifies a selenoprotein gene.
pub fn is_selenoprotein(hgnc_id: Option<i32>) -> bool {
    match hgnc_id {
        Some(id) => SELENOPROTEIN_HGNC_IDS.contains(&id),
        None => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_selenoproteins() {
        assert!(is_selenoprotein(Some(2884))); // DIO1
        assert!(is_selenoprotein(Some(10752))); // SELENOP
        assert!(is_selenoprotein(Some(30399))); // SELENOM
    }

    #[test]
    fn non_selenoprotein() {
        assert!(!is_selenoprotein(Some(12345)));
        assert!(!is_selenoprotein(Some(0)));
    }

    #[test]
    fn null_hgnc_id() {
        assert!(!is_selenoprotein(None));
    }
}
