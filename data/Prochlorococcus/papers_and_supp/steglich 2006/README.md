# Steglich 2006

**Citation:** Steglich C, Futschik M, Rector T, Steen R, Chisholm SW. "Genome-Wide Analysis of Light Sensing in *Prochlorococcus*." *Journal of Bacteriology* 188(22), 7796–7806 (2006).
**DOI:** 10.1128/JB.01097-06
**Organism(s):** *Prochlorococcus* MED4 (axenic)
**Topic:** Affymetrix microarray screen of MED4 responses to 45-min exposure (after 5 h dark adaptation) to six light conditions — high white (55 µE), white/blue/red/green at 13 µE, and DCMU + white light — versus a darkness control. Asks which genes are responsive to light quality vs intensity vs PSII electron flow.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `steglich-et-al-2006-genome-wide-analysis-of-light-sensing-in-prochlorococcus.pdf` | PDF | Main paper | reference | — |
| `suppltable1.pdf` | PDF | PDF version of Table S1 (same content as csv) | skip | Skip (duplicate of csv already in). |
| `supp table legend.txt` | text | Legend for Table S1 | reference | — |
| `supp table s1.csv` | CSV | Genes significant in ≥1 of 6 conditions; log2FC columns HL/dark, WL/dark, BL/dark, RL/dark, GL/dark, DCMU/WL with asterisk-encoded p-values | already in | — |
| `supp table s1_resolved.csv` | CSV | Auto-generated resolution output | — | — (auto-generated; skip) |
| `supp table s1_resolved_report.txt` | text | Resolution report | — | — (auto-generated; skip) |
| `paperconfig.yaml` | YAML | Current integration config | reference | — |
| `paperconfig_orig.yaml` | YAML | Superseded draft | — | — (skip / historical) |

## Current paperconfig summary

- Experiments defined: 6 (one per light condition: high white, white-13, blue-13, red-13, green-13, DCMU+white-13), all MED4 MICROARRAY.
- Statistical analyses (DE edges): 6 — one per experiment, each reading a different log2FC column from the same Table S1 CSV. `pvalue_asterisk_in_logfc: true` throughout (significance encoded as `*`/`**`/`***` on the fold-change value).
- Supplementary materials entry types: 1× `csv` (wide table powering 6 analyses).
- Organisms covered: *Prochlorococcus* MED4.
- Table scope: `filtered_subset` — genes significant in at least one of the six contrasts; all contrasts reported per row.
- Non-DE evidence: none.
- ID resolution: `gene ID` column → `old_locus_tag` (PMM#### style); direct MED4 mapping.

## Recommended actions

1. **No action** on DE — six-condition integration is a clean factorial fit to the asterisk-encoded data.
2. **Consider** — the paper's hierarchical clustering (grouping light responses) produces a qualitative grouping (e.g. "HL-and-BL-responsive", "DCMU-specific") that could be encoded as a `derived_metrics_table` (categorical, small vocabulary) if the clustering table is recoverable from the paper's supplementary figure. Currently PDF-only.
3. **Skip** — `suppltable1.pdf` is a format duplicate of the csv.
4. **Verify** — with `pvalue_asterisk_in_logfc: true` and `logfc_threshold: 1.0`, confirm that genes below both asterisk threshold and fold threshold are classified `not_significant` not dropped.

## Notes

- Gene IDs: MED4 `PMM####` as `old_locus_tag`. Native tier-1 lookup.
- No continuous adjusted_p_value is available — significance comes from asterisk count. `expression_status` will effectively be `significant_*` when any asterisk is present AND |log2FC| ≥ 1.0.
- All six conditions share a single CSV because the paper reported the full matrix on genes significant anywhere; this is the correct `csv` structure (one row → many analyses).
