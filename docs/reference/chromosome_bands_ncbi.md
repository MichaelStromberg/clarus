# Chromosome Bands: NCBI Data Source Research

## Problem

The reference file format (Section 3.3 of the spec) defines a cytogenetic band block per chromosome, but the writer was hardcoding `band_count = 0`. We need an authoritative data source for chromosome band coordinates.

## Approaches Investigated

### 1. GFF3 `map=` Attribute (Unsuitable)

NCBI's GFF3 files include `map=` attributes on region features, e.g. `map=GRCh38.p14;1;p36.33`. However, these map scaffold accessions to band names — they don't provide the genomic coordinate ranges for each band. A single `map=` entry tells you which band a scaffold belongs to, not where the band starts and ends on the chromosome.

### 2. NCBI `assembly_regions.txt` (Insufficient)

The assembly regions file (`*_assembly_regions.txt`) provides centromere and PAR boundaries, but does not include cytogenetic band data.

### 3. UCSC `cytoBand.txt` (Valid but 0-based)

UCSC hosts `cytoBand.txt` files with band coordinates in **0-based, half-open** intervals. This would require coordinate conversion to match Clarus's 1-based closed convention. The band boundaries are identical to NCBI's data.

### 4. NCBI GDP Ideogram Files (Recommended)

NCBI hosts authoritative ideogram files at:
```
ftp.ncbi.nlm.nih.gov/pub/gdp/
```

The file `ideogram_9606_GCF_000001305.16_850_V1` contains 862 bands for GRCh38 primary chromosomes (1-22, X, Y) with **1-based, closed intervals** — matching Clarus's coordinate convention exactly. No coordinate conversion needed.

**URL:** `https://ftp.ncbi.nlm.nih.gov/pub/gdp/ideogram_9606_GCF_000001305.16_850_V1`

**Note on accession:** `GCF_000001305.16` is the **Primary Assembly unit** accession within the full assembly `GCF_000001405.40` (GRCh38.p14). Bands only apply to primary chromosomes, so NCBI keys ideogram files to this unit accession (confirmed from the assembly report header: `GCF_000001305.16 → Primary Assembly`).

## File Format

Tab-separated, 9 columns:

| Column | Name       | Description                          |
|--------|------------|--------------------------------------|
| 0      | chromosome | Ensembl-style name (1, 2, ..., X, Y) |
| 1      | arm        | p or q                               |
| 2      | band       | Band number (e.g. 36.33)             |
| 3      | iscn_start | ISCN start position                  |
| 4      | iscn_stop  | ISCN stop position                   |
| 5      | bp_start   | Genomic start (1-based inclusive)     |
| 6      | bp_stop    | Genomic stop (1-based inclusive)      |
| 7      | stain      | Stain type (gneg, gpos, acen, etc.)  |
| 8      | density    | Optional stain density               |

Example:
```
#chromosome  arm  band   iscn_start  iscn_stop  bp_start   bp_stop    stain  density
1            p    36.33  0           100        1          2300000    gneg
1            p    36.32  100         244        2300001    5300000    gpos   25
```

Band names are constructed as `arm + band` (e.g. `p36.33`). Chromosome names use the Ensembl convention, matching the `ensembl_name` field in our assembly report parser.

## Columns Used

We only need: chromosome (0), arm (1), band (2), bp_start (5), bp_stop (6).

## Recommendation

Use the NCBI GDP ideogram file `ideogram_9606_GCF_000001305.16_850_V1` for GRCh38.p14. It provides:
- Authoritative NCBI data
- 1-based closed intervals (no conversion needed)
- Ensembl-style chromosome names (direct match to our data model)
- 862 bands across 24 primary chromosomes
