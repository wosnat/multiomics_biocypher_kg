# Huang 2020 (published 2021)

**Citation:** Huang S, Sun Y, Zhang S, Long L (2021). Temporal transcriptomes of a marine cyanopodovirus and its *Synechococcus* host during infection. *MicrobiologyOpen* 10(1):e1150.
**DOI:** 10.1002/mbo3.1150
**GEO:** GSE150732
**Organism(s):** *Synechococcus* sp. WH7803 (host); cyanopodovirus S-SBP1 (phage, attached via `Tests_coculture_with` as the generic Phage taxon 10239)
**Topic:** RNA-seq time course of S-SBP1 infection at MOI 2. Samples at 15 min, 1, 3, 5 and 7 h post-infection plus two uninfected controls (30 min, 8 h). Latent period 8 h, burst size ~30.

Directory is named for the GEO submission year; the paper published in 2021.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `MicrobiologyOpen - 2021 - Huang - ....pdf` | PDF | Main paper | reference | — |
| `GSE150732_WH7803.Differential_Expression.xls` | XLS | Host DE, 2528 rows x 5 timepoints (GEO original) | source | kept untouched |
| `GSE150732_WH7803_DE.csv` | CSV | Sheet1 exported for the pipeline | already in | — |
| `GSE150732_WH7803_DE_resolved.csv` | CSV | Pre-resolved locus tags (step 4) | generated | — |
| `GSE150732_WH7803_response_groups.csv` | CSV | The 32 classified genes + short `response_group` key | already in | — |
| `cluster_extractions/huang2021_response_groups.json` | JSON | Per-group descriptions (hand-authored from the Results text) | already in | — |
| `GSE150732_WH7803.rpkm.xls` | XLS | Per-sample host RPKM | reference | no DE contrast; see D1 in `plans/geo_paperconfig_updates.md` |
| `GSE150732_S-SBP1.rpkm.xls` | XLS | Per-sample **phage** RPKM | skip | phage genes are not KG Gene nodes |
| `geo_source.txt` | TXT | GEO provenance + column notes | reference | — |
| `data sources.txt` | TXT | Original download note | reference | — |

## Current paperconfig summary

- Experiments: 1 — WH7803 S-SBP1 infection (`treatment_type: [viral]`, `growth_phase: infected`)
- Statistical analyses: 5 — one per timepoint (0.25, 1, 3, 5, 7 h), all vs the 30 min uninfected control
- Clustering analyses: 1 — `huang2021_response_groups`, 3 GeneClusters (1a n=9, 1b n=20, 2 n=3), 32 membership edges
- Table scope: `all_detected_genes` (2528 host genes at every timepoint, unfiltered)
- ID resolution: **97.2%** (2456/2528), tier-1 on native `SYNWH7803_RS#####`, plus 68 rows recovered via the gene-name column. **Zero locus-tag collisions.**
- Expected edges: 2456 x 5 = **12,280** `Changes_expression_of`

## Notes

- **Thresholds are the authors' own**, not the adapter defaults: Methods state "Host genes with log2FC >= 2 as well as p < 0.05 were considered to be upregulated", so `logfc_threshold: 2.0`.
- **p-values are UNADJUSTED.** The source reports `Pvalue-*`, and the authors' own cut is on the unadjusted p. The paperconfig schema has no `pvalue_col`, so the raw p rides in `adjusted_p_value_col` and `Changes_expression_of.adjusted_p_value` holds an unadjusted p for this publication. See the comment block in `paperconfig.yaml`.
- Both controls (30 min, 8 h) were statistically indistinguishable, so the authors compared everything to the 30 min control.
- The 72 unresolved rows are RNA features (`rrf`) and locus tags absent from the current RefSeq annotation.
## Response groups (gene_clusters)

The paper clusters normalised transcript profiles across the five timepoints and classifies the **32 significantly upregulated** genes into three groups (Figure 6 / Table 3):

| key | n | what it is |
|---|---|---|
| `1a` | 9 | Carbon fixation / carboxysome — `rbcL`, RuBisCO small subunit, carboxysome shell, CcmK. Immediate upregulation |
| `1b` | 20 | Ribosomal proteins (all 20 members). Upregulated 1–3 h |
| `2` | 3 | SOS / stress — `umuC`, `umuD` and one hypothetical. Upregulated late, 5–7 h |

- **"Group 3" is deliberately excluded.** It is the residual label the source table gives the other 2,496 rows ("other genes were downregulated or remained unaltered"), not a discovered cluster. Including it would create a 2,424-member GeneCluster spanning ~98% of the genome and report `cluster_count` as 4. `GSE150732_WH7803_DE.csv` still holds every row.
- The entry keys on a `response_group` column (`1a`/`1b`/`2`) added alongside the original `Group` labels, because the key becomes part of the GeneCluster node id and every other `gene_clusters` entry in the repo uses space-free keys.
- Cluster descriptions in `cluster_extractions/` were **hand-authored from the Results paragraph**, not produced by an LLM extraction run — the paper describes each group explicitly, so no inference was needed. `metadata.model` records this.
- All 32 genes resolve (100%).
