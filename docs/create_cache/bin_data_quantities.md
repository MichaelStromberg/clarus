# Per-Bin Data Quantity Maximums

Captured from GRCh38 pipeline run on 2026-02-23 using RefSeq GCF_000001405.40-RS_2025_08 and Ensembl v115.

## RefSeq

| Field               |   Max | Chromosome | Bin |
|---------------------|------:|------------|----:|
| transcripts         |   462 | chr17      |  41 |
| genes               |   108 | chr14      |  96 |
| transcript_regions  |  8803 | chr17      |  41 |
| cdna_sequences      |   462 | chr17      |  41 |
| protein_sequences   |   204 | chr17      |  41 |
| regulatory_regions  |   354 | chr19      |   1 |

## Ensembl

| Field               |   Max | Chromosome | Bin |
|---------------------|------:|------------|----:|
| transcripts         |  1347 | chr1       |   0 |
| genes               |   156 | chr14      | 101 |
| transcript_regions  | 22190 | chr6       |  30 |
| cdna_sequences      |  1344 | chr1       |   0 |
| protein_sequences   |   590 | chr6       |  30 |
| regulatory_regions  |   470 | chr9       | 130 |

## Combined Maximums (across both sources)

| Field               |   Max | Source  | Fits in |
|---------------------|------:|---------|---------|
| transcripts         |  1347 | Ensembl | u16     |
| genes               |   156 | Ensembl | u8      |
| transcript_regions  | 22190 | Ensembl | u16     |
| cdna_sequences      |  1344 | Ensembl | u16     |
| protein_sequences   |   590 | Ensembl | u16     |
| regulatory_regions  |   470 | Ensembl | u16     |

---

# Per-Reference Deduplication Maximums

Captured from the same GRCh38 pipeline run. With per-reference deduplication (all bins within a chromosome share a single deduplicated array), the unique counts are substantially larger than per-bin.

## RefSeq

| Field               |     Max | Chromosome |
|---------------------|--------:|------------|
| genes               |   2,975 | chr1       |
| transcript_regions  | 128,425 | chr1       |
| cdna_sequences      |   8,699 | chr1       |
| protein_sequences   |   5,072 | chr1       |
| transcripts         |   8,739 | chr1       |
| regulatory_regions  |  13,790 | chr1       |

## Ensembl

| Field               |     Max | Chromosome |
|---------------------|--------:|------------|
| genes               |   7,058 | chr1       |
| transcript_regions  | 601,535 | chr1       |
| cdna_sequences      |  47,906 | chr1       |
| protein_sequences   |  15,452 | chr1       |
| transcripts         |  47,973 | chr1       |
| regulatory_regions  |  35,922 | chr1       |

## Combined Maximums (across both sources)

| Field               |     Max | Source  | Fits in |
|---------------------|--------:|---------|---------|
| genes               |   7,058 | Ensembl | u16     |
| transcript_regions  | 601,535 | Ensembl | u32     |
| cdna_sequences      |  47,906 | Ensembl | u16     |
| protein_sequences   |  15,452 | Ensembl | u16     |
| transcripts         |  47,973 | Ensembl | u16     |
| regulatory_regions  |  35,922 | Ensembl | u16     |

---

# Per-Field Integer Maximums

Captured from the same GRCh38 pipeline run. These are the maximum values observed for each integer field across all transcripts, used to determine whether fixed-width types (u8, u16, u32) can replace VLQ encoding.

## RefSeq

| Struct / Field                    |    Max | Type |
|-----------------------------------|-------:|------|
| **CodingRegion**                  |        |      |
| cdna_start                        |   4706 | i32  |
| cdna_end                          | 108201 | i32  |
| cds_padding                       |      2 | u8   |
| cds_offset                        |    209 | u16  |
| protein_offset                    |     69 | u16  |
| amino_acid_edits (count)          |     10 |      |
| **AminoAcidEdit**                 |        |      |
| position                          |    858 | u32  |
| **TranslationalSlip**             |        |      |
| position                          |   1186 | i32  |
| length                            |      2 | u8   |
| **TranscriptRegion**              |        |      |
| cdna_start                        | 107906 | i32  |
| cdna_end                          | 109224 | i32  |
| id                                |    363 | u16  |
| cigar_ops (count)                 |     39 |      |
| **CigarOp**                       |        |      |
| length                            |   6352 | u32  |

## Ensembl

| Struct / Field                    |    Max | Type |
|-----------------------------------|-------:|------|
| **CodingRegion**                  |        |      |
| cdna_start                        |  11620 | i32  |
| cdna_end                          | 108465 | i32  |
| cds_padding                       |      2 | u8   |
| cds_offset                        |      0 | u16  |
| protein_offset                    |   1218 | u16  |
| amino_acid_edits (count)          |     10 |      |
| **AminoAcidEdit**                 |        |      |
| position                          |    858 | u32  |
| **TranslationalSlip**             |        |      |
| position                          |      0 | i32  |
| length                            |      0 | u8   |
| **TranscriptRegion**              |        |      |
| cdna_start                        | 108170 | i32  |
| cdna_end                          | 347561 | i32  |
| id                                |    365 | u16  |
| cigar_ops (count)                 |      0 |      |
| **CigarOp**                       |        |      |
| length                            |      0 | u32  |

## Combined Maximums (across both sources)

| Struct / Field                    |    Max | Source  | Fits in |
|-----------------------------------|-------:|---------|---------|
| **CodingRegion**                  |        |         |         |
| cdna_start                        |  11620 | Ensembl | u16     |
| cdna_end                          | 108465 | Ensembl | u32     |
| cds_padding                       |      2 | Both    | u8      |
| cds_offset                        |    209 | RefSeq  | u16     |
| protein_offset                    |   1218 | Ensembl | u16     |
| amino_acid_edits (count)          |     10 | Both    | u8      |
| **AminoAcidEdit**                 |        |         |         |
| position                          |    858 | Both    | u16     |
| **TranslationalSlip**             |        |         |         |
| position                          |   1186 | RefSeq  | u16     |
| length                            |      2 | RefSeq  | u8      |
| **TranscriptRegion**              |        |         |         |
| cdna_start                        | 108170 | Ensembl | u32     |
| cdna_end                          | 347561 | Ensembl | u32     |
| id                                |    365 | Ensembl | u16     |
| cigar_ops (count)                 |     39 | RefSeq  | u8      |
| **CigarOp**                       |        |         |         |
| length                            |   6352 | RefSeq  | u16     |
