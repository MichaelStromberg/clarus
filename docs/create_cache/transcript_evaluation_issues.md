# Transcript Evaluation Issues

Analysis of unresolvable transcripts encountered during `create_cache` for GRCh38.

## Evaluation Hierarchy

The transcript evaluation module (`src/evaluation/mod.rs`) implements a 7-level resolution hierarchy from specification 07:

1. **Perfect match** — translated CDS exactly equals expected protein
2. **Contained** — translated protein aligns within expected with no edits
3. **AA edits** — aligned with acceptable amino acid edits (2 for normal, 10 for selenoproteins)
4. **Translational slippage** — ribosomal frameshift with slip length <= 2
5. **Frame correction +1** — prepend 1 base to CDS and re-evaluate
6. **Frame correction +2** — prepend 2 bases to CDS and re-evaluate
7. **Unresolvable** — transcript cannot be resolved; coding region is removed

## Root Cause: Alignment Limitation with Longer Translated Proteins

### The Problem

`semi_global_align(reference, query)` requires `query.len() <= reference.len()`. In the frame correction path (`try_frame_correction`), the call is `semi_global_align(expected, &translated)`, treating `expected` as reference and `translated` as query.

When the original CDS length is not a multiple of 3, prepending 1-2 bases can make the adjusted CDS a multiple of 3, producing one extra complete codon. This makes `translated` 1 AA longer than `expected`, causing the alignment to unconditionally return `None`.

### Why Prepending Creates Longer Proteins

Example with padding = 1:
- Original CDS: 92 bases (92 % 3 = 2) -> translates to 30 complete codons + 2 leftover = 30 AAs
- After prepend 1 base: 93 bases (93 % 3 = 0) -> translates to 31 codons = 31 AAs
- `translated` (31 AAs) > `expected` (30 AAs) -> forward alignment fails

### The Fix

Added a reverse alignment fallback in `try_frame_correction`. When `translated.len() > expected.len()`, the function calls `semi_global_align(&translated, expected)` (swapping reference and query), then uses `find_amino_acid_edits_reverse` to compute edits from the reversed perspective.

## Affected Transcripts (GRCh38)

From `create_cache` output showing 37 unresolvable transcripts:

### Frame +1 Matches (7 unique transcripts, ~19 entries)

These transcripts have CDS length % 3 == 2. Prepending 2 bases (Frame +1 diagnostic = skip-1 = prepend-2 reading frame) produces a matching protein:

- **NM_144686.4** (7 entries across alt contigs)
- **NM_001145303.3** (7 entries across alt contigs)
- NM_000314.8
- NM_001004431.3
- NM_001331018.1
- NM_001331021.1
- NM_001331022.1

### Frame +2 Matches (2 unique transcripts, ~3 entries)

CDS length % 3 == 1. Prepending 1 base produces a matching protein:

- **NM_001006607.3** (2 entries)
- NM_001025200.4

### No Frame Match (~13 unique transcripts, ~14 entries)

These transcripts do not match any frame shift and may be truly unresolvable or require other mechanisms:

- NM_001098721.2, NM_001098722.2, NM_001170820.4, NM_001289160.2
- NM_001301773.2, NM_001301774.2, NM_001304388.2, NM_001348294.2
- NM_004485.4, NM_030975.2, NM_152348.4, NM_173683.4
- NM_199340.5 (2 entries)

## Diagnostic Output Interpretation

The `dump_unresolvable_diagnostics` function outputs Frame +1 and Frame +2 translations by skipping 1 or 2 bases from CDS start. This is the **reverse** of what `try_frame_correction` does (prepending bases):

| Diagnostic label | Operation | `try_frame_correction` equivalent |
|---|---|---|
| Frame +1 (skip 1 base) | `translate(cds[1..])` | Level 6: padding=2 (prepend 2) |
| Frame +2 (skip 2 bases) | `translate(cds[2..])` | Level 5: padding=1 (prepend 1) |

## Phase 1 Outcome

With the reverse alignment fallback in `try_frame_correction`:
- Frame +1 and Frame +2 match transcripts (~22 entries, ~9 unique) resolved via levels 5-6
- `num_unresolvable` dropped from 37 to 12

---

## Phase 2: Post-Stop Translation Artifacts and Reverse Alignment

After Phase 1, 12 transcripts remained unresolvable. Analysis revealed two fixable root causes and two genuinely unresolvable categories.

### Root Cause: `codon::translate` continues past stop codons

`codon::translate` does NOT stop at `*` — it translates all remaining codons and appends `X` for incomplete trailing codons. This causes alignment failures because the translated protein has extra residues after the stop codon that don't appear in the expected protein.

### Category A: Post-stop translation artifacts (8 transcripts, fixed)

**A1 — Extra residues after stop at Level 2-3 (3 transcripts):**
- **NM_001304388.2**: translated = expected + `"Y"` (1 extra AA after stop)
- **NM_030975.2**: translated = expected + `"SSPKR"` (5 extra AAs after stop)
- **NM_199340.5** (×2 entries): translated = expected + `"EE"` (2 extra AAs after stop)

These failed because `semi_global_align(expected, translated)` requires `query.len() <= reference.len()`, and translated was longer than expected.

**A2 — Extra residues after stop prevent frame correction alignment (4 transcripts):**
- **NM_001098722.2**, **NM_004485.4**, **NM_001098721.2**: Frame +2 translation matches expected but ends with `*X` (incomplete trailing codon)
- **NM_001170820.4**: Frame +2 translation matches expected but ends with `*G` (extra codon after stop)

These failed because trailing junk after `*` prevented the diagonal aligner from finding a match.

**Fix:** Added `trim_after_stop()` helper that truncates translated protein after the first `*`. Applied in both `evaluate_single` and `try_frame_correction` after each `codon::translate` call. Selenoproteins are excluded since `*` represents selenocysteine (U) in those.

### Category B: Extra residue at START (1 transcript, fixed)

- **NM_152348.4**: translated = `"P"` + expected (extra `P` at start from upstream CDS)

**Fix:** Added reverse alignment at Level 2-3 in `evaluate_single`. When `translated.len() > expected.len()`, calls `semi_global_align(translated, expected)` with swapped arguments.

### Category C: Genuinely unresolvable (3 transcripts, 4 entries)

**C1 — Partial CDS (2 transcripts):**
- **NM_001301774.2**, **NM_001301773.2**: CDS = 18 bp encoding only 6 AAs for a 96 AA expected protein. These are genuinely incomplete CDS records.

**C2 — Same-length frame shift (1 transcript, 2 entries):**
- **NM_001348294.2** (×2): Frame +1 translation matches expected exactly, but CDS length % 3 == 0, so no padding change occurs. The diagonal-only aligner cannot detect shifts when translated and expected are the same length but offset.

### Phase 2 Outcome

- `num_unresolvable` drops from 12 to 4 (the Category C transcripts)
- Category A transcripts resolve via perfect match or Level 2-3 alignment after trimming
- Category B transcript resolves via reverse alignment at Level 2-3

---

## Phase 3: Regression Fix — Remove `trim_after_stop`, Replace `Err(Validation)` with Fall-Through

Phase 2 resolved 8 of 12 unresolvable transcripts, but introduced regressions that were discovered by comparing `create_cache_output2.txt` (pre-Phase 2) vs `create_cache_output3.txt` (post-Phase 2).

### Regressions Discovered

| Stat | Before Phase 2 | After Phase 2 | Delta |
|---|---|---|---|
| Perfect protein | 77318 | 77328 | +10 |
| Contained | 592 | 600 | +8 |
| AA edits | 483 | 490 | +7 |
| **Slippage** | **17** | **0** | **-17** |
| Wrong frame | 294 | 290 | -4 |
| Unresolvable | 12 | 4 | -8 |

All 17 slippage cases were broken, and 3 new unresolvable transcripts appeared on chr7.

### Root Cause 1: `trim_after_stop` breaks slippage detection

`trim_after_stop` cuts the translated protein at the first `*` (stop codon). But some CDS sequences have internal TGA codons that produce `*` in the translation. Cutting there shortens the protein before slippage detection (Level 4), completely breaking all 17 slippage cases.

### Root Cause 2: `Err(Validation)` blocks fall-through to lower levels

When alignment found too many edits (>2 for non-selenoprotein), `Err(Validation)` was returned. This was caught by `evaluate_transcripts` the same way as `Err(UnresolvableTranscript)` — the transcript was counted as unresolvable and its coding region removed. This prevented Levels 4-6 from being tried, causing 3 new chr7 unresolvable entries.

### Fix

**Removed `trim_after_stop` entirely.** The Phase 2 reverse alignment blocks already handle all post-stop and extra-AA cases without trimming:
- Extra AAs after stop → reverse alignment finds expected within translated
- Extra AA at start → reverse alignment finds expected within translated at offset

**Replaced `Err(Validation)` returns with fall-through.** In all four alignment blocks (evaluate_single forward/reverse, try_frame_correction forward/reverse), transcripts with too many edits now fall through to try lower levels (slippage, frame correction) instead of erroring out.

### Phase 3 Outcome

- Slippage restored to ~17
- chr7 regressions resolved (transcripts now reach slippage/frame correction)
- Unresolvable count: ~2-4 (only genuinely irrecoverable: partial CDS + same-length frame shift)

---

## Phase 4: CIGAR-Aware Coordinate Mapping for Frameshift Transcripts

After Phase 3, `create_cache` crashed on NM_001396027.1 with a fatal error:
```
Error: CDS coordinates [1, 698] out of bounds for cDNA length 696 in NM_001396027.1
```

### The Transcript

**NM_001396027.1** (FAM246C, chr22) is a polymorphic pseudogene with 2 single-nucleotide deletions in the GRCh38 reference genome relative to the actual mRNA. The GFF3 models this with:
- 3 exons (241 + 35 + 420 = 696 bp) separated by 1-bp genomic gaps (frameshifts)
- 1 cDNA_match spanning the full region: `Gap=M241 D1 M35 D1 M420`
- Genomic span = 698 bp, cDNA span = 696 bp

### Root Cause: `build_cigar_aware` + linear `map_genomic_to_cdna`

`build_cigar_aware` created a **single** TranscriptRegion for the entire cDNA_match with genomic span 698 bp and cDNA span 696 bp. When `detect_coding_region` called `map_genomic_to_cdna` to map the CDS end position, the linear formula `cdna_start + (genomic_pos - genomic_start)` computed `1 + 697 = 698` instead of the correct `696`, because the 2 D (deletion) ops make the genomic span 2 bp longer than the cDNA span.

The CIGAR ops were stored in the TranscriptRegion but never consulted during coordinate mapping.

### Secondary Issue: Fatal Error Propagation

Phase 3 removed `Error::Validation` from the `evaluate_transcripts` catch clause (correctly, to enable fall-through for alignment edit limits). But the CDS bounds check in `evaluate_single` still returned `Err(Validation)`, which was no longer caught — propagating as a fatal pipeline crash instead of being treated as an unresolvable transcript.

### Fix

**Split CIGAR-aware regions at D/I ops** (`src/transcript/construction.rs`): Refactored `build_cigar_aware` to walk each MatchRecord's CIGAR operations and create a separate sub-region for each M (match) block. D ops create genomic gaps (that `insert_introns` fills). I ops advance cDNA without advancing genomic position. Each sub-region has a 1:1 genomic-to-cDNA mapping, so the existing linear formula in `map_genomic_to_cdna` works correctly.

For `Gap=M241 D1 M35 D1 M420`, the single region is split into 3 sub-regions with 1-bp introns between them — exactly matching the GFF3 exon structure.

**Changed CDS bounds check to UnresolvableTranscript** (`src/evaluation/mod.rs`): As a defensive measure, changed the error type so CDS bounds violations are caught and treated as unresolvable rather than crashing the pipeline.

### Phase 4 Outcome

- Pipeline no longer crashes on NM_001396027.1
- NM_001396027.1 evaluates successfully (CDS 1..696 = 231 aa + stop matches GenBank protein)

---

## Phase 5: Frame Correction Preprocessing and Offset Fallback

After Phase 4, `create_cache` completed successfully with 11 unresolvable transcripts (all on alt contigs/fix patches). Investigation of these transcripts revealed three fixable failure mechanisms in the frame correction path (Levels 5-6).

### Root Cause Analysis

All 11 unresolvable transcripts are on alt contigs where the genomic sequence diverges from the primary assembly. The frame correction path (`try_frame_correction`) was blocked by three artifacts:

**1. Trailing X from incomplete codons**

When CDS + padding is not a multiple of 3, `codon::translate` produces `X` for the incomplete trailing codon. This extends the translated protein past what the reference allows, causing `semi_global_align` to reject the alignment (query extends past the end of the reference at every diagonal position).

Example: CDS = 92 bases, padding = 2 → 94 bases. 94 % 3 = 1 → 31 complete codons + 1 leftover base = 31 AAs + X. The trailing X prevents alignment.

**2. Post-stop translation artifacts**

`codon::translate` continues past stop codons (`*`), translating all remaining codons. In frame correction, this produces characters after the stop codon (e.g., `*G`, `*X`) that cause extra mismatches, pushing the edit count above the max_edits threshold (2 for non-selenoproteins).

**3. Padding-induced extra codon at start**

Prepending 1-2 bases creates a new first codon that translates to a junk amino acid. For very short proteins, the single mismatch from this junk AA can push the per-base alignment score below the threshold (10.0), causing the diagonal-only aligner to reject the alignment entirely.

### Affected Transcripts

| Group | Transcripts | Contig | Mechanism |
|-------|-------------|--------|-----------|
| A | NM_001098722.2, NM_004485.4, NM_001098721.2 | chr1_KQ458384v1_alt | Frame +2 matches but trailing X prevents alignment |
| B | NM_000296.4, NM_001009944.3 | chr16_KI270853v1_alt | Partial CDS (45 bp for 370+ AA protein) — genuinely unresolvable |
| C | NM_001170820.4 | chr11_KI270903v1_alt | Frame +2 matches but post-stop `*G` prevents alignment |
| D | NM_001080978.4, NM_005874.5 | chr19_GL949752v1_alt | Translated offset by 8 positions — genuinely unresolvable |
| E | NM_001278403.3, NM_001278404.3 | chr19_GL949752v1_alt | No frame matches — genuinely unresolvable |
| F | NM_001012709.1 | chr11_MU273369v1_fix | Frame +2 matches but padding-induced extra codon prevents alignment |

### Fix

Three preprocessing steps added to `try_frame_correction` (`src/evaluation/mod.rs`), applied after `codon::translate` and before alignment:

1. **Strip trailing X**: Remove all trailing `X` characters (incomplete codon artifacts from leftover bases after padding).

2. **Truncate after first `*`** (non-selenoproteins only): Characters after the first stop codon are translation artifacts. Truncating removes them. Selenoproteins are excluded because `*` represents selenocysteine (U) in those proteins. This is safe in `try_frame_correction` because it's past slippage detection (Level 4).

3. **Offset fallback**: After trying alignment with the full translated protein, try again with `translated[1..]` (skipping the padding-induced junk first AA). This handles the case where the junk codon doesn't match any position in the expected protein at sufficient alignment quality.

A new `try_frame_alignment` helper was extracted to handle both forward and reverse alignment with edit counting, called by both the main alignment attempt and the offset fallback.

### Phase 5 Outcome

All 11 previously unresolvable transcripts now resolve via frame correction:

| Stat | Before Phase 5 | After Phase 5 | Delta |
|---|---|---|---|
| Perfect protein | 77271 | 77271 | 0 |
| Contained | 584 | 610 | +26 |
| AA edits | 546 | 531 | -15 |
| Slippage | 10 | 10 | 0 |
| Wrong frame | 272 | 283 | +11 |
| **Unresolvable** | **11** | **0** | **-11** |

The +11 wrong_frame matches the 11 previously unresolvable transcripts, all now resolved. The shift from AA edits to contained (+26/-15) reflects improved alignment quality from the preprocessing (fewer spurious post-stop mismatches). Zero regressions in slippage or perfect protein counts.
