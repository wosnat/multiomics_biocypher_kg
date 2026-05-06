# Moreno 2023

**Citation:** Moreno-Cabezuelo JÁ, Gómez-Baena G, Díez J, García-Fernández JM. "Integrated Proteomic and Metabolomic Analyses Show Differential Effects of Glucose Availability in Marine *Synechococcus* and *Prochlorococcus*." *Microbiology Spectrum* 11(2), e03275-22 (2023).
**DOI:** 10.1128/spectrum.03275-22
**Organism(s):** *Prochlorococcus* MED4, SS120 (and MIT9313 in paper, not in DE tables here); *Synechococcus* WH7803, WH8102, BL107; heterotroph fractions *Alteromonas* (MarRef v6) and *Marinobacter* (MarRef v6) as coculture partners.
**Topic:** LC-MS/MS quantitative proteomics (MarRef v6 database match) of 5 non-axenic cyano cultures × 4 glucose × light conditions (light, light+100 nM glc, light+5 mM glc, dark, dark+100 nM glc, dark+5 mM glc; comparisons are glucose-vs-no-glucose within light or dark). Heterotroph fractions (Alteromonas, Marinobacter) were proteomically profiled in parallel from the same cultures. Reference-proteome-match organism pattern (no cultured strain genomes for the heterotroph side).

## Classification

**Bucket A — DE integration complete** (secondary: bucket B for the metabolomics side, now actionable since Phase 2 metabolomics support landed).

All three sides of the proteome DE — cyano S2, Alteromonas S3, Marinobacter S4 — are wired up and resolve natively. Cyano S2 (5 strains × 4 comparisons) resolves against the deployed cyano genomes via the parsed `KEGG` / `extracted_gn` columns (97.1%–100%). Alteromonas S3 resolves natively against the deployed `Alteromonas (MarRef v6)` reference proteome (GCA_003513035.1, locus prefix DEH24_*) at 100% across all 5 cocultures (the originally planned DEH24 → MADE_RS bridge was superseded when GCA_003513035.1 was deployed directly as a reference-proteome-match organism). Marinobacter S4 resolves against `Marinobacter (MarRef v6)` (HP15 / GCF_000166295.1) at 97.2%–97.8%. The paper's metabolomic flux side (Calvin vs OPP, fermentation evidence) is metabolite-level — Metabolite/MetaboliteAssay support has now landed, but the published data is qualitative flux inference (not per-metabolite quantitation per condition) so the natural fit is sparse: a small number of `metabolite_assays_table` boolean flags per Experiment (`Assay_flags_metabolite`) for the metabolites the authors call out. Inspect supplementary material before deciding scope.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `moreno-cabezuelo-et-al-2023-integrated-proteomic-and-metabolomic-analyses-show-differential-effects-of-glucose.pdf` | PDF | Main paper | reference | — |
| `spectrum.03275-22-s0005.pdf` | PDF | Supp PDF (figures/legends) | reference | — |
| `legend.txt` | text | Column-abbreviation legend for S1–S4 | reference | — |
| `spectrum.03275-22-s0001.xlsx` | XLSX | Supp Table S1 — compiled proteome quantitation across all cyanos (all strains in one sheet) | skip | Skip — superseded by the per-strain CSVs in S2 (entry `supp_table_s1` is not defined; S2 is authoritative per-strain). |
| `spectrum.03275-22-s0002.xlsx` | XLSX | Supp Table S2 — per-strain cyano proteomes | skip | Skip (Excel source; per-strain `_modified.csv` copies are loaded). |
| `spectrum.03275-22-s0003.xlsx` | XLSX | Supp Table S3 — per-culture Alteromonas proteomes | skip | Skip (Excel source; per-culture `_modified.csv` copies are loaded). |
| `spectrum.03275-22-s0004.xlsx` | XLSX | Supp Table S4 — per-culture Marinobacter proteomes | skip | Skip (Excel source; per-culture `_modified.csv` copies are loaded). |
| `table s2 MED4 ...s0002.csv`, `table s2 SS120 ...s0002.csv`, `table s2 WH7803 ...s0002.csv`, `table s2 WH8102 ...s0002.csv`, `table s2 BL107 ...s0002.csv` (originals) | CSV | Unmodified per-strain cyano tables | skip | Skip (superseded by `_modified.csv` siblings — `scripts/build_modified_csv/...` adds the parsed locus_tag / log2FC / t-test columns that the paperconfig expects). |
| `table s2 <strain> ...s0002_modified.csv` (×5: MED4, SS120, WH7803, WH8102, BL107) | CSV | Cyano proteome DE — 4 glucose/light comparisons per strain | already in | — |
| `table s3 Alteromonas in <strain> cultures ...s0003.csv` (×5: MED4, SS120, WH7803, WH8102, BL107) (originals) | CSV | Unmodified per-culture Alteromonas tables | skip | Skip (superseded by `_modified.csv`). |
| `table s3 Alteromonas in <strain> cultures ...s0003_modified.csv` (×5) | CSV | Alteromonas (MarRef v6) proteome DE from each cyano coculture; DEH24_* locus tags resolve natively against GCA_003513035.1 | already in | — |
| `table s3 Alteromonas in WH8102 cultures ...s0003.xlsx` | XLSX | Excel original of one S3 sheet | skip | Skip (duplicate format). |
| `table s4 Marinobacter in <strain> cultures ...s0004.csv` (×5: MED4, SS120, WH7803, WH8102, BL107) (originals) | CSV | Unmodified per-culture Marinobacter tables | skip | Skip (superseded by `_modified.csv`). |
| `table s4 Marinobacter in <strain> cultures ...s0004_modified.csv` (×5) | CSV | Marinobacter (MarRef v6) proteome DE from each cyano coculture | already in | — |
| `*_modified_resolved.csv` / `*_modified_resolved_report.txt` (many) | CSV/text | Auto-generated resolution outputs (prepare_data step 4) | — | — (auto-generated; skip) |
| `paperconfig.yaml` | YAML | Current integration config | reference | — |

## Current paperconfig summary

- Experiments defined: 60 across 7 organisms — 20 cyano experiments (5 strains × 4 glucose/light conditions) + 20 Alteromonas (5 cocultures × 4 conditions) + 20 Marinobacter (5 cocultures × 4 conditions). All 60 are active and have statistical_analyses entries.
- Statistical analyses (DE edges): 60 entries across 15 supplementary_materials `csv` blocks (5 cyano S2 + 5 Alteromonas S3 + 5 Marinobacter S4).
- Supplementary materials entry types: 15× `csv` (all `_modified.csv`); no `gene_clusters` or `derived_metrics_table` entries.
- Organisms covered: *Prochlorococcus* MED4, SS120; *Synechococcus* WH7803, WH8102, BL107; heterotrophs *Alteromonas (MarRef v6)* and *Marinobacter (MarRef v6)* with `organism_type: reference_proteome_match`.
- Table scope: `all_detected_genes` — full quantitative proteomes.
- Non-DE evidence: none integrated. Metabolomics side (Calvin/OPP flux inferences, fermentation evidence) is not per-gene and is not modeled.
- ID resolution (all three sides resolve natively against deployed assemblies):
  - **cyano S2** — `KEGG` / `extracted_gn` columns parsed by `build_moreno2023_s3_modified_csv.py` upstream, resolve against native locus_tag per strain. Resolution: MED4 100% (327/327), SS120 100% (930/930), WH7803 100% (359/359), WH8102 100% (895/895), BL107 97.1% (168/173).
  - **Alteromonas S3** — MarRef v6 prefix `DEH24_*` via `locus_tag_deh24` column resolves natively against `Alteromonas (MarRef v6)` (GCA_003513035.1). Resolution: 100% across all 5 tables (MED4 16/16, SS120 54/54, WH7803 4/4, WH8102 78/78, BL107 11/11).
  - **Marinobacter S4** — MarRef v6 locus tags resolved against `Marinobacter (MarRef v6)` (HP15 / GCF_000166295.1). Resolution: 97.2%–97.8% across all 5 tables.

## Recommended actions

1. **No action** — all 60 DE analyses are wired up and resolving (97%–100% across all 15 source tables). The DEH24 → MADE_RS bridge originally planned for the Alteromonas side was superseded by deploying GCA_003513035.1 directly as `Alteromonas (MarRef v6)`.
2. **Skip** — all `.xlsx` originals and non-`_modified` CSVs; they are superseded by the per-strain `_modified.csv` files referenced in paperconfig.
3. **Consider add (Phase 2 metabolomics, now unblocked)** — the paper reports a metabolomic flux inference (glucose → OPP vs Calvin) that is qualitative and experiment-level. Phase 2 `metabolite_assays_table` could carry the called-out metabolites as boolean flags per Experiment (`Assay_flags_metabolite`); inspect supplementary content before scoping.

## Notes

- Reference-proteome-match organism pattern: Alteromonas (MarRef v6) and Marinobacter (MarRef v6) are NOT cultured strains in this study — they are whatever heterotrophs in the non-axenic cultures matched the MarRef v6 reference database. They carry `organism_type: reference_proteome_match` and `reference_database: "MarRef v6"` on their OrganismTaxon nodes. See `docs/kg-changes/reference-proteome-match-organisms.md`.
- `_modified.csv` files are generated by a committed script (`scripts/build_modified_csv/build_moreno2023_s3_modified_csv.py`) that parses the per-row Description field to extract locus_tag, log2FC, and t-test columns under stable names — the paperconfig references the modified versions, not the originals.
- Significance: t-test p-values (`t_test_*` columns) with `pvalue_threshold: 0.05`. Directionality is encoded in log2FC sign of glucose-vs-control; `growth_phase: stationary` across all analyses.
- SS120 Marinobacter `light_high_glucose` analysis was recovered after an initial column-detection issue (commit `a410d91`); all 5 SS120 Marinobacter comparisons are now active.
