# Martiny 2006

**Citation:** Martiny AC, Coleman ML, Chisholm SW. "Phosphate acquisition genes in *Prochlorococcus* ecotypes: Evidence for genome-wide adaptation." *PNAS* 103(33), 12552–12557 (2006).
**DOI:** 10.1073/pnas.0601301103
**Organism(s):** *Prochlorococcus* MED4 (HL) and MIT9313 (LL)
**Topic:** Microarray time course (0, 4, 12, 24, 48 h) comparing P-starved vs P-replete cultures for the two ecotypes. Identifies the phoB regulon and ecotype-specific differences in the P-starvation response. Significance: q < 0.05 at t=48 h; full temporal profiles reported for all detected genes passing that filter.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `martiny-et-al-2006-phosphate-acquisition-genes-in-prochlorococcus-ecotypes-evidence-for-genome-wide-adaptation.pdf` | PDF | Main paper | reference | — |
| `01301table_1.csv` | CSV | MED4 DE table: 30 up / 4 down genes, FC + q-value at 5 timepoints | already in | — |
| `01301table_1.xls` | XLS | Excel original of table 1 | skip | Skip (duplicate of CSV already in). |
| `01301table_1_resolved.csv` | CSV | Auto-generated resolution output | — | — (auto-generated; skip) |
| `01301table_1_resolved_report.txt` | text | Resolution report | — | — (auto-generated; skip) |
| `01301table_2.csv` | CSV | MIT9313 DE table: 176 genes, FC + q-value at 4 timepoints (0, 4, 12, 24 h; no 48h column) | already in | — |
| `01301table_2.xls` | XLS | Excel original of table 2 | skip | Skip (duplicate of CSV already in). |
| `01301table_2_resolved.csv` | CSV | Auto-generated resolution output | — | — (auto-generated; skip) |
| `01301table_2_resolved_report.txt` | text | Resolution report | — | — (auto-generated; skip) |
| `01301table3.pdf` | PDF | Table 3: phoB gene-cluster organization across 11 strains (not per-gene DE) | skip | Skip — pathway/locus-context only; no per-gene evidence usable as DerivedMetric or edge. |
| `paperconfig.yaml` | YAML | Current integration config | reference | — |
| `paperconfig_orig.yaml` | YAML | Superseded draft | — | — (skip / historical) |

## Current paperconfig summary

- Experiments defined: 2 — `phosphorus_phosphate_starvation_med4_microarray`, `phosphorus_phosphate_starvation_mit9313_microarray`.
- Statistical analyses (DE edges): 9 total — 5 timepoints for MED4 (0/4/12/24/48 h) × 34 genes + 4 timepoints for MIT9313 (0/4/12/24 h) × 176 genes.
- Supplementary materials entry types: 2× `csv`.
- Organisms covered: *Prochlorococcus* MED4, MIT9313.
- Table scope: `filtered_subset` — genes were pre-filtered for q < 0.05 at the final timepoint, but all timepoint values are reported.
- Non-DE evidence: none integrated. (phoB-cluster content in Table 3 is locus-context, not per-gene DE.)
- ID resolution: `ORF` column → `old_locus_tag` for both strains (PMM#### / PMT####). Both CSVs use `skip_rows: 1` to jump over a title row.

## Recommended actions

1. **No action** — DE integration is complete for both strains at all reported timepoints.
2. **Skip** — `.xls` files are byte-duplicates of the CSVs already loaded.
3. **Consider add** — if Table 3's phoB-cluster organization is ever relevant as a gene attribute (e.g. "in_phoB_cluster" boolean), it could be encoded as a `derived_metrics_table` with `value_kind: boolean`. Currently lives as PDF-only prose.
4. **Verify** — asymmetric timepoint coverage: MIT9313 has no 48 h column; paperconfig correctly omits it.

## Notes

- Gene IDs: MED4 `PMM####` and MIT9313 `PMT####` old locus tags; both resolve via native `old_locus_tag` tier-1 lookup.
- No adjusted-p-value on the MED4 t=48 column where "0.0000" is reported — adapter will treat as literal 0.
- The downregulated genes in MED4 (n=4) and the up/down pattern in MIT9313 (mostly ribosomal downregulation) are all flattened into the single `Changes_expression_of` stream per experiment.
