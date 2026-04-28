# Alonso-Saez 2023

**Citation:** Alonso-Saez L, Palacio AS, Cabello AM, Robaina-Estevez S, Gonzalez JM, Garczarek L, Lopez-Urrutia A (2023). Transcriptional Mechanisms of Thermal Acclimation in *Prochlorococcus*. *mBio* 14(3).
**DOI:** 10.1128/mbio.03425-22
**Organism(s):** *Prochlorococcus marinus* MIT9301 (HLII clade representative)
**Topic:** Quantitative transcriptome (RNA-Seq) analysis of MIT9301 acclimated long-term to 5 temperatures (17, 20, 22, 25, 30C) under a 12h:12h L/D cycle. RNA samples collected 3h after subjective sunrise (day) and sunset (night). Five soft clusters (A-E) identified by fuzzy c-means clustering of day/night expression patterns across the thermal gradient; clusters correspond to core daytime genes, photosystems, stress response, respiration/DNA replication, and night-biased expression.

## Classification

**Bucket B — new metrics / DE / resolution (want to add)**

Cluster integration is already in place. The remaining gap is per-gene transcript abundance (S3 raw counts, S4 normalised counts) which fits the numeric `derived_metrics_table` pattern shipped 2026-04-21 (see `docs/kg-changes/non-de-evidence-extension.md`). No organism blocker — strain is deployed. Action: add one or more `derived_metrics_table` entries with `value_kind: numeric` referencing the BV-BRC IDs already mapped via the cluster wiring.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `alonso-sáez-et-al-2023-...pdf` | PDF | Main paper | reference | — |
| `TABLE S5 softcluster membership.csv` | CSV | Soft cluster assignments for 1831 MIT9301 genes + 5 per-cluster probability columns (A-E) | already in | — |
| `soft_clusters_mbio.03425-22-s0008.xlsx` | XLSX | Original Excel of Table S5 | skip | Duplicate of CSV already in |
| `mbio.03425-22-s0005.xlsx` | XLSX | Table S3: mRNA transcripts per cell (raw abundance, not DE) | skip | Sample x gene expression / abundance matrix |
| `mbio.03425-22-s0006.xlsx` | XLSX | Table S4: mRNA transcripts per biovolume (raw abundance, not DE) | skip | Sample x gene expression / abundance matrix |
| `supp legends.txt` | TXT | Supplementary legends | reference | — |
| `cluster_extractions/` | dir | Per-cluster description JSONs | reference | — |
| `paperconfig.yaml` | YAML | Active paperconfig | reference | — |

## Current paperconfig summary

- Experiments defined: 1 — `thermal_acclimation_mit9301_rnaseq`
- Statistical analyses (DE edges): 0 — no pairwise DE fold-changes/p-values (absolute transcript counts across thermal gradient)
- Supplementary materials entry types: `gene_clusters` (single entry `mit9301_softclusters_thermal_acclimation`, fuzzy c-means K=5)
- Organisms covered: Prochlorococcus MIT9301
- Table scope(s): n/a
- Non-DE evidence: 1 ClusteringAnalysis — 5-cluster fuzzy c-means assignment (clusters A, B, C, D, E)
- ID resolution: `RefSeq Locus Tag` (P9301_NNNNN) as `old_locus_tag` — matches 1760/1831 (96.1%); `Gene  symbol` (two spaces) as `gene_name`. 71 genes have only BV-BRC `fig|167546.4.peg.NNN` IDs and will not resolve without an added id_translation.

## Recommended actions

1. **No action** — Table S5's fuzzy c-means clustering with per-gene soft-membership probabilities is correctly represented as a `gene_clusters` entry per the decision rule (soft-clustering, one cluster per gene, per-cluster functional descriptions available from the paper text and `cluster_extractions/`).
2. **Add (optional)** — The 5 per-cluster probability columns (`Probability score Cluster A..E`) could be added as 5 numeric DerivedMetric entries (`value_kind: numeric`, `rankable: true`) to preserve the full soft-membership data, but this is low-value since membership is already captured via the hard assignment + `score_col` on the cluster edge. Defer unless a downstream query needs per-cluster probabilities.
3. **Add** — Tables S3 (transcripts/cell) and S4 (transcripts/biovolume) are single-condition per-gene quantitative snapshots. Wire as `derived_metrics_table` entries with `value_kind: numeric` per the pattern shipped 2026-04-21 (see `docs/kg-changes/non-de-evidence-extension.md`).
4. **Add id_translation (optional)** — For the 71 genes with BV-BRC `fig|` IDs only, an `id_translation` entry (diamond protein match against MIT9301 protein FASTA) would recover them; low impact (<4% of genes).
5. **Reference** — main PDF is reference-only.

## Notes

- Gene ID format: `RefSeq Locus Tag` = `P9301_NNNNN` (NCBI old locus tag); `BV-BRC_ID` = `fig|167546.4.peg.NNN` (PATRIC/BV-BRC IDs); `Gene  symbol` (note: two spaces in column name) as gene_name fallback.
- No single score column; the CSV has 5 per-cluster probability columns (A-E) so `score_col` is not set in the paperconfig. A preprocessing step could derive the max probability as a single score column if needed.
- The CSV has 2 extra header rows (title + blank line) before actual column names; `skip_rows: 2` handles this.
- `_resolved.csv` / `_resolved_report.txt` files auto-generated by `prepare_data.sh` step 4.
- Strain coverage: MIT9301 is deployed in the KG.
