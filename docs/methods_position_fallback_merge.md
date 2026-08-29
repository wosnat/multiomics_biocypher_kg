# Methods: Position-Based Fallback Merge for Gene Mapping

## Problem

The gene mapping pipeline (`build_gene_mapping.py`) constructs a unified gene table per strain by merging NCBI GCF genome annotations with Cyanorak database annotations. The primary merge strategy explodes the NCBI `old_locus_tag` attribute (a URL-encoded comma-separated list of historical identifiers) and joins on the Cyanorak `locus_tag`. When the NCBI `old_locus_tag` does not include the identifier form used by Cyanorak, the merge fails and the same gene appears as two separate rows: one NCBI-only (with genomic coordinates but no Cyanorak functional annotations) and one Cyanorak-only (with functional annotations but no NCBI coordinates or cross-references). These split entries propagate downstream as duplicate gene nodes in the knowledge graph.

## Affected Strain

The problem was identified in *Prochlorococcus marinus* MIT9313, whose genome has undergone three rounds of reannotation:

| Generation | Locus tag format | Example | Source |
|---|---|---|---|
| Original (2003) | `PMT` + 4-digit (no underscore) | `PMT0107` | Cyanorak GBK |
| NCBI re-annotation 1 | `PMT_` + 4-digit (with underscore) | `PMT_0107` | NCBI `old_locus_tag` |
| NCBI re-annotation 2 | `RG24_RS` + 5-digit | `RG24_RS00545` | NCBI `old_locus_tag` |
| Current NCBI | `AKG35_RS` + 5-digit | `AKG35_RS00545` | NCBI `locus_tag` |

For approximately 2,060 MIT9313 genes, the NCBI `old_locus_tag` contains all three historical forms (`PMT_####`, `RG24_RS#####`, and the earlier PMT0###-equivalent), and the merge succeeds. However, for approximately 360 genes, the `old_locus_tag` contains only the `PMT_####` and `RG24_RS#####` forms, omitting the `PMT0###` form that Cyanorak uses. This mismatch is an artifact of incomplete carryover during NCBI genome reannotation and is not correctable at the data source.

## Investigation

### Quantifying the split

Analysis of the MIT9313 `gene_mapping.csv` before the fix showed:

| Category | Count |
|---|---|
| Matched (both sources) | 2,159 |
| Cyanorak-only | 674 |
| NCBI-only | 230 |
| **Total rows** | **3,063** |

Of the 674 Cyanorak-only entries, 266 had PMT0### locus tags (without underscore) -- these were the split candidates. The remaining 408 had PMT_#### locus tags and were genuine Cyanorak-only genes with no NCBI counterpart.

### External translation tables

No available data source directly bridges the PMT0### and PMT_#### identifier forms:

- **NCBI GCF GFF** (`old_locus_tag`): Contains `PMT_####` and `RG24_RS#####` but not `PMT0###` for the affected genes.
- **NCBI GCA GFF**: Uses `PMT_####` as the primary `locus_tag`; no `old_locus_tag` attribute at all.
- **Cyanorak GBK**: Uses `PMT0###` as `locus_tag`; does not cross-reference to `PMT_####`.
- **MIT9313_genbank.tsv** (manually curated translation table): Maps `AKG35_RS` to `PMT_` format; does not contain a `PMT0###` column.

String-based normalization (removing the underscore to convert `PMT_0107` to `PMT0107`) was considered but rejected because the identifier correspondence is not always a simple underscore removal and cannot be reliably inferred without additional evidence.

### Position-based matching analysis

Since the NCBI and Cyanorak annotations describe the same physical genome, genes encoding the same protein occupy the same (or very nearly the same) genomic coordinates. A comprehensive coordinate comparison was performed between all 230 NCBI-only entries and all 266 PMT0### Cyanorak-only entries:

**Reciprocal overlap** was defined as `overlap_length / max(ncbi_gene_length, cyanorak_gene_length)`, where `overlap_length = max(0, min(ncbi_end, cyan_end) - max(ncbi_start, cyan_start))`.

Results at various thresholds:

| Threshold | Matches | Conflicts | Clean 1:1 |
|---|---|---|---|
| Overlap >= 50% | 164 | 2 | 162 |
| Overlap >= 90% | 124 | 2 | 122 |
| Overlap >= 90% + start <= 10bp | 113 | 1 | 112 |
| Overlap >= 90% + start <= 10bp + end <= 10bp | 107 | 1 | **105** |

Key observations:
- **End coordinates are highly conserved**: All matches with >= 90% overlap had an end coordinate difference of exactly 0 bp (identical stop codon position). Start coordinates varied due to different CDS start-site annotation between NCBI and Cyanorak.
- **One conflict**: `AKG35_RS05630` (NCBI) matched two Cyanorak entries (`PMT2281` and `PMT2283`) at identical coordinates -- likely a gene fusion/split discrepancy between annotations.
- The 11 matches excluded by the start <= 10bp filter had start differences of 12-108 bp, suggesting potential annotation differences rather than simple coordinate shifts. These were conservatively excluded.

## Implementation

### Identity criterion (revised 2026-08-29)

The first implementation (2026-03) accepted a Cyanorak/NCBI pair only when the
two intervals had ≥ 90 % reciprocal overlap, start coordinates within 50 bp and
end coordinates within 3 bp, strand-blind. That encoded the *wrong* invariant.
A bacterial CDS is identified by its **stop codon and reading frame**, not by
its start: re-annotation routinely moves the start codon (PGAP vs. Cyanorak
differ by up to ~250 bp on MIT9313, with reciprocal overlaps down to 0.18), and
on the − strand the start codon is the genomic *end*, so the ±3 bp end gate
rejected every − strand start re-call. Cyanorak also sometimes reports a CDS
without its stop codon (3′ end 3 bp short) or truncated a few codons early,
still in frame. The thresholds therefore missed 336 same-gene pairs corpus-wide
(42 on MIT9313 — the `PMT0040` / `PMT_0040` "two locus-tag families" backlog
item), while accepting, in principle, two genes whose stop codons differ by
≤ 3 bp (a different frame).

`_position_fallback_merge(ncbi_sourced, cyan_only)` now accepts a Cyanorak
call **C** for an unmatched NCBI call **N** when:

1. **Same contig** — NCBI coordinates are shifted by a per-contig offset (see
   below) before comparison.
2. **Same strand.**
3. **In frame** — the 3′ ends (genomic `end` on +, genomic `start` on −)
   differ by a multiple of 3.
4. **C's 3′ end lies inside N's interval.** In frame and inside means no stop
   codon can separate the two 3′ ends, so C is the same ORF as N up to the
   start-codon choice or a missing stop codon — by construction, not by a
   tolerance.
5. **1 : 1 in both directions** — an NCBI gene matching several Cyanorak calls
   (a fusion/split disagreement, e.g. MIT9313 `AKG35_RS05630` ↔ `PMT2281` +
   `PMT2283`) or a Cyanorak call landing in frame inside two nested NCBI genes
   is skipped with a warning.

No overlap ratio, start tolerance or end tolerance remains, and no distance
bound inside N's interval either: an in-frame Cyanorak ORF ending inside an
NCBI span is the same locus even when N is a nonsense pseudogene or
frameshifted CDS that Cyanorak called as its intact upstream ORF.

### Contig awareness (new 2026-08-29)

Cyanorak stores a draft genome as **one concatenated record** (NCBI contigs in
order), whereas NCBI coordinates are contig-relative. The first implementation
ignored this and produced cross-contig false merges on the two multi-contig
Cyanorak strains (PAC1, 20 contigs; SB, 4). `_contig_offsets()` derives, per
NCBI `seqid`, the constant that maps NCBI to Cyanorak coordinates from the
genes the locus_tag merge already paired, measured at the 3′ end; the offset is
used when a strict majority of a contig's matched genes agree on it (contigs backed by fewer than three matched genes are used with a warning) (they agree 100 % in
practice: PAC1 contig 5 = +229 662, SB contig 2 = +420 116, every closed genome
= 0). A contig with no matched gene is usable at offset 0 only in a
single-sequence assembly; otherwise the fallback never compares coordinates on
it and says so in the step-0 log.

### Traceability and conflict policy

Unchanged: every merged row carries `position_merge_note`
(`position_merge:<cyanorak_tag>→<ncbi_locus_tag>`) and the consumed Cyanorak
locus tag is appended to `old_locus_tags`, so papers keyed on either form
resolve to one gene. Conflicts (rule 5) are logged and left unmerged.

## Verification

The criterion was validated against sequence, not against thresholds: for
every candidate pair the Cyanorak interval was translated from the Cyanorak
GenBank record (table 11) and compared with the NCBI protein (`protein.faa`).
Over the 21 Cyanorak strains, all **923** position-merged pairs in the rebuilt
`gene_mapping.csv` files encode the NCBI protein (one is an in-frame suffix of
the other, no internal stop) — **0 false merges**. The same check on the
pre-revision output found the 11 stop-codon-excluded pairs the old gates had
caught correctly and the 2 cross-contig pairs (PAC1 `EV03_0060→EV03_0395`, SB
`EV02_0773→EV02_1244`) they had merged wrongly.

## Results

Rebuilt 2026-08-29 (sub-step 5, all strains):

| Strain | Rows before → after | Position-merged rows |
|---|---|---|
| MIT9313 | 2,948 → 2,906 (−42) | 157 (was 115) |
| WH8102 | 2,881 → 2,830 (−51) | 174 |
| PAC1 | 2,370 → 2,317 (−53) | 53 (was 1, and wrong) |
| MIT9202 | 2,027 → 2,005 (−22) | 80 |
| BL107 | 2,595 → 2,567 (−28) | 61 |
| WH7803 | 2,621 → 2,590 (−31) | 55 |
| 15 other Cyanorak strains | −196 in total | — |
| **All 21** | **−423** | **923** |

On MIT9313 the 42 collapsed rows were exactly the `PMTnnnn` (Cyanorak-only) /
`PMT_nnnn` (NCBI-only) twins: NCBI's `old_locus_tag` lists `PMT_0040` but not
`PMT0040` for `AKG35_RS00210`, and Cyanorak's `PMT0040` starts 93 bp upstream
with the same stop codon. They are one gene. The remaining `PMT_`-only rows
(≈ 400) are genuine Cyanorak-only calls with no NCBI counterpart — Cyanorak
itself names those `PMT_nnnn` — not a stale second annotation build.

## Test Coverage

`tests/test_build_gene_mapping.py`:

| Test | Scenario |
|---|---|
| `test_position_fallback_merges_unmatched_pair` | identical coordinates, same strand → merged |
| `test_position_fallback_merges_start_codon_recall_plus_strand` | MIT9313 `PMT0040`: + strand, start re-called 93 bp downstream, shared stop (overlap 0.76) → merged |
| `test_position_fallback_merges_start_codon_recall_minus_strand` | MIT9313 `PMT0110`: − strand, genomic `end` moved 24 bp, shared stop → merged |
| `test_position_fallback_merges_cyanorak_stop_codon_excluded_plus` / `_minus` | MED4 `PMM1858a` / `PMM1719`: Cyanorak 3′ end 3 bp short, in frame → merged |
| `test_position_fallback_is_contig_aware` | multi-contig: a Cyanorak gene at the *raw* coordinates of a contig-2 NCBI gene is not merged; one at the offset coordinates is |
| `test_position_fallback_skips_different_strand` | same coordinates, opposite strands → not merged |
| `test_position_fallback_skips_shared_start_different_stop` | stops differ by 3 bp (passed the old ±3 bp gate) → not merged |
| `test_position_fallback_skips_out_of_frame_3prime_end` | stops differ by 4 bp → not merged |
| `test_position_fallback_skips_different_stop_codon_far` / `_near` | stop codons 1000 bp / 10 bp apart → not merged |
| `test_position_fallback_merges_nested_in_frame_orf` | Cyanorak ORF ending in frame 600 bp inside the NCBI span → merged (pins the absence of a distance bound) |
| `test_position_fallback_note_names_ncbi_locus_tag_without_old_locus_tag` | NCBI gene with no `old_locus_tag`: note reads `→<locus_tag_ncbi>`, never `→None` |
| `test_position_fallback_skips_reverse_conflict` | one Cyanorak call in frame inside two nested NCBI genes → neither merged |
| `test_position_fallback_skips_conflict` | one NCBI gene, two Cyanorak calls at the same stop → neither merged |

## Applicability to Other Strains

The fallback runs unconditionally. It fires on every strain with a Cyanorak
annotation (21 as of 2026-08-29; MIT1327's Cyanorak layer is derived from
NCBI's own GBFF, so nothing is left for it to do) and on none of the
NCBI-only strains (there is no Cyanorak side). Its output should be re-checked
with the translation test above whenever a new Cyanorak strain is onboarded.
