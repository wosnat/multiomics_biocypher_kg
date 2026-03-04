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

### Algorithm

A position-based fallback merge was added to `load_gff_from_ncbi_and_cyanorak()` in `build_gene_mapping.py`. The fallback runs after the primary locus_tag-based merge and before the final concatenation, operating only on entries that failed the primary merge.

**Function**: `_position_fallback_merge(ncbi_sourced, cyan_only, min_overlap=0.9, max_bp_diff=10)`

For each unmatched Cyanorak entry (PMT0### locus tag, no NCBI counterpart):

1. Select all unmatched NCBI entries on the same strand.
2. Compute reciprocal overlap, start difference, and end difference against each candidate using vectorized NumPy operations.
3. Accept candidates meeting all three criteria: reciprocal overlap >= 90%, |start difference| <= 10 bp, |end difference| <= 10 bp.
4. Collect all candidate pairs keyed by NCBI locus tag.

After scanning all Cyanorak entries:

5. For each NCBI locus tag with exactly one Cyanorak candidate: merge by copying all Cyanorak annotation columns into the NCBI row in-place.
6. For NCBI locus tags with multiple Cyanorak candidates (conflicts): skip the merge entirely and emit a warning.
7. Add a `position_merge_note` column to each merged row recording the Cyanorak-to-NCBI locus tag mapping (e.g., `position_merge:PMT0107->PMT_0107`).

The consumed Cyanorak locus tags are returned so the caller can remove them from the Cyanorak-only set before concatenation, preventing duplicate rows.

### Traceability

Every position-merged entry is annotated with a `position_merge_note` value in `gene_mapping.csv`. This column is null for all entries merged by the standard locus_tag-based method, making it straightforward to audit which genes were merged by position and to revert or adjust the merge if needed.

### Conflict policy

When one NCBI gene matches multiple Cyanorak entries by position, neither merge is performed. This conservative policy avoids incorrectly collapsing genes that may represent annotation disagreements (e.g., gene fusions in one annotation, separate genes in the other). Conflicts are logged as warnings for manual review.

## Results

MIT9313 `gene_mapping.csv` after the fix:

| Category | Before | After | Change |
|---|---|---|---|
| Both sources (locus_tag merge) | 2,159 | 2,159 | -- |
| Position-merged | 0 | 105 | +105 |
| Both sources (total) | 2,159 | 2,264 | +105 |
| Cyanorak-only | 674 | 569 | -105 |
| NCBI-only | 230 | 125 | -105 |
| Conflicts skipped | -- | 1 | -- |
| **Total rows** | **3,063** | **2,958** | **-105** |

The 105 position-merged entries correspond to genes where the Cyanorak `PMT0###` locus tag had no matching `old_locus_tag` in the NCBI GFF. After the merge, these genes have both NCBI genomic coordinates/cross-references and Cyanorak functional annotations (cluster assignments, ontology terms, pathway data) in a single row, eliminating the downstream duplicate gene node problem.

The remaining 569 Cyanorak-only and 125 NCBI-only entries represent genuine differences between the two annotation sources (genes annotated in one but not the other) and are retained as single-source rows.

## Test Coverage

Five unit tests validate the position fallback behavior:

| Test | Scenario |
|---|---|
| `test_position_fallback_merges_unmatched_pair` | NCBI and Cyanorak entries at overlapping coordinates on the same strand are merged into one row |
| `test_position_fallback_skips_different_strand` | Entries at identical coordinates but on opposite strands are not merged |
| `test_position_fallback_skips_low_overlap` | Entries with < 90% reciprocal overlap are not merged |
| `test_position_fallback_skips_large_coord_diff` | Entries with start difference > 10 bp are not merged |
| `test_position_fallback_skips_conflict` | When one NCBI entry matches two Cyanorak entries, neither is merged |

## Applicability to Other Strains

This fix is specific to strains where the NCBI `old_locus_tag` is incomplete relative to the identifiers used by Cyanorak. Among the 13 strains in the knowledge graph, MIT9313 is the only one exhibiting this pattern. For all other strains, the locus_tag-based merge succeeds for the full set of shared genes, and the position fallback produces zero merges (the function exits early when no unmatched NCBI entries exist). The fallback is therefore safe to run unconditionally for all strains without affecting their existing merge results.
