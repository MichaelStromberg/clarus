//! Codon translation tables (standard and mitochondrial).

/// Lookup table for translating codons to amino acids.
///
/// Indexed by 6-bit codon encoding: A=0, C=1, G=2, T/U=3.
/// Index = first*16 + second*4 + third.
pub struct CodonTable {
    table: [u8; 64],
}

fn base_to_index(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

impl CodonTable {
    /// Standard genetic code (NCBI translation table 1).
    #[must_use]
    pub fn standard() -> Self {
        // Amino acids indexed by codon encoding (A=0, C=1, G=2, T=3)
        // Order: AAA, AAC, AAG, AAT, ACA, ACC, ACG, ACT, AGA, AGC, AGG, AGT,
        //        ATA, ATC, ATG, ATT, CAA, ...
        #[rustfmt::skip]
        let table: [u8; 64] = [
            b'K', b'N', b'K', b'N',  // AA*: AAA=K, AAC=N, AAG=K, AAT=N
            b'T', b'T', b'T', b'T',  // AC*: ACA=T, ACC=T, ACG=T, ACT=T
            b'R', b'S', b'R', b'S',  // AG*: AGA=R, AGC=S, AGG=R, AGT=S
            b'I', b'I', b'M', b'I',  // AT*: ATA=I, ATC=I, ATG=M, ATT=I
            b'Q', b'H', b'Q', b'H',  // CA*: CAA=Q, CAC=H, CAG=Q, CAT=H
            b'P', b'P', b'P', b'P',  // CC*: CCA=P, CCC=P, CCG=P, CCT=P
            b'R', b'R', b'R', b'R',  // CG*: CGA=R, CGC=R, CGG=R, CGT=R
            b'L', b'L', b'L', b'L',  // CT*: CTA=L, CTC=L, CTG=L, CTT=L
            b'E', b'D', b'E', b'D',  // GA*: GAA=E, GAC=D, GAG=E, GAT=D
            b'A', b'A', b'A', b'A',  // GC*: GCA=A, GCC=A, GCG=A, GCT=A
            b'G', b'G', b'G', b'G',  // GG*: GGA=G, GGC=G, GGG=G, GGT=G
            b'V', b'V', b'V', b'V',  // GT*: GTA=V, GTC=V, GTG=V, GTT=V
            b'*', b'Y', b'*', b'Y',  // TA*: TAA=*, TAC=Y, TAG=*, TAT=Y
            b'S', b'S', b'S', b'S',  // TC*: TCA=S, TCC=S, TCG=S, TCT=S
            b'*', b'C', b'W', b'C',  // TG*: TGA=*, TGC=C, TGG=W, TGT=C
            b'L', b'F', b'L', b'F',  // TT*: TTA=L, TTC=F, TTG=L, TTT=F
        ];
        Self { table }
    }

    /// Vertebrate mitochondrial genetic code (NCBI translation table 2).
    /// Differences from standard: UGA→W, AGA→*, AGG→*, AUA→M.
    #[must_use]
    pub fn mitochondrial() -> Self {
        let mut table = Self::standard();
        // TGA: index = 3*16 + 2*4 + 0 = 56 → W (was *)
        table.table[56] = b'W';
        // AGA: index = 0*16 + 2*4 + 0 = 8 → * (was R)
        table.table[8] = b'*';
        // AGG: index = 0*16 + 2*4 + 2 = 10 → * (was R)
        table.table[10] = b'*';
        // ATA: index = 0*16 + 3*4 + 0 = 12 → M (was I)
        table.table[12] = b'M';
        Self { table: table.table }
    }

    /// Translate a single codon (3 bytes) to an amino acid.
    #[must_use]
    pub fn translate_codon(&self, codon: &[u8]) -> u8 {
        if codon.len() < 3 {
            return b'X';
        }
        let i0 = base_to_index(codon[0]);
        let i1 = base_to_index(codon[1]);
        let i2 = base_to_index(codon[2]);
        match (i0, i1, i2) {
            (Some(a), Some(b), Some(c)) => self.table[a * 16 + b * 4 + c],
            _ => b'X',
        }
    }
}

/// Translate a CDS nucleotide sequence to a protein sequence.
pub fn translate(cds: &[u8], table: &CodonTable) -> Vec<u8> {
    let mut protein = Vec::with_capacity(cds.len() / 3 + 1);
    let mut i = 0;
    while i + 2 < cds.len() {
        protein.push(table.translate_codon(&cds[i..i + 3]));
        i += 3;
    }
    // Incomplete codon at end
    if i < cds.len() {
        protein.push(b'X');
    }
    protein
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_start_codon() {
        let table = CodonTable::standard();
        assert_eq!(table.translate_codon(b"ATG"), b'M');
    }

    #[test]
    fn standard_stop_codons() {
        let table = CodonTable::standard();
        assert_eq!(table.translate_codon(b"TAA"), b'*');
        assert_eq!(table.translate_codon(b"TAG"), b'*');
        assert_eq!(table.translate_codon(b"TGA"), b'*');
    }

    #[test]
    fn mitochondrial_differences() {
        let table = CodonTable::mitochondrial();
        assert_eq!(table.translate_codon(b"TGA"), b'W'); // was *
        assert_eq!(table.translate_codon(b"AGA"), b'*'); // was R
        assert_eq!(table.translate_codon(b"AGG"), b'*'); // was R
        assert_eq!(table.translate_codon(b"ATA"), b'M'); // was I
    }

    #[test]
    fn translate_short_orf() {
        let table = CodonTable::standard();
        // ATG GCA TGC TAA = M A C *
        let cds = b"ATGGCATGCTAA";
        let protein = translate(cds, &table);
        assert_eq!(protein, b"MAC*");
    }

    #[test]
    fn translate_incomplete_codon() {
        let table = CodonTable::standard();
        // ATG GC (incomplete) → M X
        let cds = b"ATGGC";
        let protein = translate(cds, &table);
        assert_eq!(protein, b"MX");
    }

    #[test]
    fn translate_empty() {
        let table = CodonTable::standard();
        let protein = translate(b"", &table);
        assert!(protein.is_empty());
    }

    #[test]
    fn ambiguous_base() {
        let table = CodonTable::standard();
        assert_eq!(table.translate_codon(b"NNN"), b'X');
        assert_eq!(table.translate_codon(b"ATN"), b'X');
    }
}
