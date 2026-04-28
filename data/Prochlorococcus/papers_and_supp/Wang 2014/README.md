# Wang 2014

**Citation:** Wang J, Chen L, Huang S, Liu J, Zhang X, Qin S (2014). A transcriptomic study of the *Prochlorococcus* cells' differential response to various nutrients. *BMC Microbiology* 14:11.
**DOI:** 10.1186/1471-2180-14-11
**Organism(s):** *Prochlorococcus marinus* MED4 (axenic)
**Topic:** RNA-seq-based profiling of MED4 grown in Pro99 and AMP media at exponential and log-phase (10 samples total), 21C, continuous light at 28 umol quanta m-2 s-1. All MED4 CDS genes were classified by RPKM quartiles into expression-level categories (VEG/HEG/MEG/LEG/NEG) and Core/Flexible designation; operon predictions are also reported.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `1471-2180-14-11.pdf` | PDF | Main paper | reference | — |
| `GeneClassification_clean.csv` | CSV | MED4 CDS genes classified by expression level (VEG/HEG/MEG/LEG/NEG/--) + Core/Flexible + COG. Cleaned version (title row removed, footnote suffixes stripped) | already in | — |
| `GeneClassification 12866_2013_2178_MOESM3_ESM.csv` | CSV | Original CSV of gene classification (pre-clean) | skip | Superseded by `GeneClassification_clean.csv` |
| `12866_2013_2178_MOESM3_ESM.xlsx` | XLSX | Original Excel of gene classification table | skip | Duplicate |
| `operons 12866_2013_2178_MOESM1_ESM.csv` | CSV | MED4 operon predictions (operon ID, coordinates, member genes) | skip | Operon structure; no node type for operons in KG today |
| `12866_2013_2178_MOESM1_ESM.xlsx` | XLSX | XLSX of operon predictions | skip | Duplicate |
| `12866_2013_2178_MOESM2_ESM.xlsx` | XLSX | Supplementary Table 2 | skip | Pending inspection — if per-gene evidence, re-evaluate |
| `12866_2013_2178_MOESM4_ESM.xlsx` | XLSX | Supplementary Table 4 | skip | Pending inspection — if per-gene evidence, re-evaluate |
| `cluster_extractions/` | dir | Per-cluster description JSONs | reference | — |
| `paperconfig.yaml` | YAML | Active paperconfig | reference | — |

## Classification

**Bucket B — new metrics / DE / resolution (want to add)**

The paper has no pairwise DE — output is per-gene RPKM-quartile categorical labels (VEG/HEG/MEG/LEG/NEG/--) plus a Core/Flexible flag. Currently encoded as a `gene_clusters` entry but the labels are a small fixed vocabulary (no membership scores, no per-cluster narratives) and should be migrated to `derived_metrics_table` with `value_kind: categorical` (allowed_categories listed) + a second boolean DerivedMetric for Core/Flexible. MOESM2 and MOESM4 XLSX files have not been inspected and may carry per-gene evidence. Operon predictions (MOESM1) defer until an Operon node type exists. MED4 is deployed; no new strain needed. Action is a paperconfig edit.

## Current paperconfig summary

- Experiments defined: 1 — `expression_profiling_med4_rnaseq`
- Statistical analyses (DE edges): 0 — no pairwise DE; quartile classification only
- Supplementary materials entry types: `derived_metrics_table` (single entry `med4_expression_level`) with two categorical metrics:
  - `expression_level_class` — `ExpressionLevel` column, `allowed_categories: [VEG, HEG, MEG, LEG, NEG, "--"]`. `field_description` documents the meaning of every value (RPKM-quartile rank + COG / DEG / Ka enrichment notes from Wang 2014 Figures 3-4).
  - `pangenome_membership` — `Core/Flexible` column, `allowed_categories: [Core, Flexible]`. `field_description` documents Core vs Flexible membership and the paper's NEG/VEG core-vs-flexible composition statistics.
- Organisms covered: Prochlorococcus MED4
- Table scope(s): n/a
- Non-DE evidence: 2 DerivedMetrics — ExpressionLevel (VEG 1081 / MEG 445 / HEG 291 / LEG 82 / NEG 66 / `--` 89; 5 NaN skipped) and Core/Flexible (Core 1251 / Flexible 714; 94 NaN skipped).
- ID resolution: `sysName` (PMM####) as `locus_tag`; `locus_tag` column (PMED4_xxxxx) as `alternative_locus_tag`; `name` as `gene_name`. `desc` declared in `product_columns`.

## Recommended actions

1. **Done (2026-04-28)** — migrated `ExpressionLevel` and `Core/Flexible` from `gene_clusters` to `derived_metrics_table` with two categorical metrics. Per-value meanings drawn from `cluster_extractions/med4_expression_level.json` are inlined in each metric's `field_description`. `cluster_extractions/` files are no longer consumed by the adapter (the categorical labels are not soft-clustering memberships) but kept for provenance. `validate_paperconfig.py` passes.
2. **Skip** — operon predictions (MOESM1) do not correspond to any node type currently in the schema; defer until an Operon node type is added.
3. **Add (evaluate)** — inspect MOESM2 and MOESM4 XLSX content; if either contains per-gene DE, periodicity flags, or numeric scores, add as `csv` or `derived_metrics_table`. Not inspected in this pass.
4. **Reference** — main PDF is reference-only.
5. **Optional follow-up** — register `expression_level_class` and `pangenome_membership` in `multiomics_kg/vocab/non_de_evidence.py` if other papers begin emitting equivalent metrics (currently advisory warnings, not errors).

## Notes

- Gene ID format: `sysName` = `PMM####` (native MED4 Cyanorak/NCBI locus tag, direct match); `locus_tag` column = `PMED4_xxxxx` (JGI IMG draft) bridged via `gene_id_mapping.json`.
- `GeneClassification_clean.csv` differs from the raw CSV by having the title row removed and footnote suffixes stripped from column names; original CSV kept for provenance.
- `_resolved.csv` / `_resolved_report.txt` files auto-generated by `prepare_data.sh` step 4.
- Strain coverage: MED4 is deployed in the KG.
- The citation given in the previous README referenced the wrong paper (Wang 2014 Nannochloropsis); corrected here to the BMC Microbiology 14:11 article (Prochlorococcus MED4 transcriptomics) matching the actual data files.
