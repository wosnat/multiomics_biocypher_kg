# Moreno 2023

**Citation:** Moreno-Cabezuelo JÁ, Gómez-Baena G, Díez J, García-Fernández JM. "Integrated Proteomic and Metabolomic Analyses Show Differential Effects of Glucose Availability in Marine *Synechococcus* and *Prochlorococcus*." *Microbiology Spectrum* 11(2), e03275-22 (2023).
**DOI:** 10.1128/spectrum.03275-22
**Organism(s):** *Prochlorococcus* MED4, SS120 (and MIT9313 in paper, not in DE tables here); *Synechococcus* WH7803, WH8102, BL107; heterotroph fractions *Alteromonas* (MarRef v6) and *Marinobacter* (MarRef v6) as coculture partners.
**Topic:** LC-MS/MS quantitative proteomics (MarRef v6 database match) of 5 non-axenic cyano cultures × 4 glucose × light conditions (light, light+100 nM glc, light+5 mM glc, dark, dark+100 nM glc, dark+5 mM glc; comparisons are glucose-vs-no-glucose within light or dark). Heterotroph fractions (Alteromonas, Marinobacter) were proteomically profiled in parallel from the same cultures. Reference-proteome-match organism pattern (no cultured strain genomes for the heterotroph side).

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
| `table s3 Alteromonas in <strain> cultures ...s0003_modified.csv` (×5) | CSV | Alteromonas (MarRef v6) proteome DE from each cyano coculture; DEH24_* locus tags | already in | — **but DE match rate is ~0 until DEH24 → MADE_RS (or direct DEH24) bridge is built** (flagged at the top of paperconfig.yaml). |
| `table s3 Alteromonas in WH8102 cultures ...s0003.xlsx` | XLSX | Excel original of one S3 sheet | skip | Skip (duplicate format). |
| `table s4 Marinobacter in <strain> cultures ...s0004.csv` (×5: MED4, SS120, WH7803, WH8102, BL107) (originals) | CSV | Unmodified per-culture Marinobacter tables | skip | Skip (superseded by `_modified.csv`). |
| `table s4 Marinobacter in <strain> cultures ...s0004_modified.csv` (×5) | CSV | Marinobacter (MarRef v6) proteome DE from each cyano coculture | already in | — |
| `*_modified_resolved.csv` / `*_modified_resolved_report.txt` (many) | CSV/text | Auto-generated resolution outputs (prepare_data step 4) | — | — (auto-generated; skip) |
| `paperconfig.yaml` | YAML | Current integration config | reference | — |

## Current paperconfig summary

- Experiments defined: 60 across 7 organisms — 20 cyano experiments (5 strains × 4 glucose/light conditions) + 20 Alteromonas (5 cocultures × 4 conditions) + 20 Marinobacter (5 cocultures × 4 conditions). Note the Marinobacter in SS120 is missing one column (the FC G mM/L column), leaving 39 "active" experiments plus 20 Alteromonas experiments that are presently deferred-loaded (see id-resolution note below).
- Statistical analyses (DE edges): ~60 entries across 15 supplementary_materials `csv` blocks (5 cyano S2 + 5 Alteromonas S3 + 5 Marinobacter S4).
- Supplementary materials entry types: 15× `csv` (all `_modified.csv`); no `gene_clusters` or `derived_metrics_table` entries.
- Organisms covered: *Prochlorococcus* MED4, SS120; *Synechococcus* WH7803, WH8102, BL107; heterotrophs *Alteromonas (MarRef v6)* and *Marinobacter (MarRef v6)* with `organism_type: reference_proteome_match`.
- Table scope: `all_detected_genes` — full quantitative proteomes.
- Non-DE evidence: none integrated. Metabolomics side (Calvin/OPP flux inferences, fermentation evidence) is not per-gene and is not modeled.
- ID resolution: **cyano S2** — `locus_tag_*` columns parsed by `build_moreno2023_s3_modified_csv.py` upstream, then resolved via native `locus_tag` per strain. **Alteromonas S3** — custom MarRef v6 prefix `DEH24_*`; paperconfig uses `locus_tag_deh24` as `locus_tag`-typed name_col. Comments inside the paperconfig flag that "DEH24 → MADE_RS bridge is needed" and match rates are expected to be low until the bridge or direct GCA_003513035.1 deployment is complete. **Marinobacter S4** — similar MarRef v6 locus tags, resolved against HP15/Marinobacter MarRef reference proteome (GCF_000166295.1).

## Recommended actions

1. **No action** on cyano S2 integration — 20 DE analyses across 5 strains are wired up and should resolve cleanly against native locus tags.
2. **Build DEH24 bridge** — the top-of-file TODO in paperconfig.yaml is the blocker for ~20 Alteromonas expression edges per coculture. Options: (a) deploy GCA_003513035.1 directly so DEH24_* resolves natively, OR (b) build an `id_translation` entry with `method: diamond_protein_match` mapping DEH24 → MADE_RS against the MADE_RS reference. Once built, rerun prepare_data steps 3+4 and rebuild KG.
3. **Verify Marinobacter match rate** — Marinobacter MarRef reference proteome is deployed as `Marinobacter (MarRef v6)` (see CLAUDE.md). Run `/check-gene-ids` on one `_modified.csv` sample to confirm DEH24-style identifiers also resolve there (the MarRef v6 Marinobacter assembly is GCF_000166295.1/HP15).
4. **Skip** — all `.xlsx` originals and non-`_modified` CSVs; they are superseded.
5. **Consider add** — the paper reports a metabolomic flux inference (glucose → OPP vs Calvin) that is qualitative and experiment-level, not per-gene; skip for now.

## Notes

- Reference-proteome-match organism pattern: Alteromonas (MarRef v6) and Marinobacter (MarRef v6) are NOT cultured strains in this study — they are whatever heterotrophs in the non-axenic cultures matched the MarRef v6 reference database. They carry `organism_type: reference_proteome_match` and `reference_database: "MarRef v6"` on their OrganismTaxon nodes. See `docs/kg-changes/reference-proteome-match-organisms.md`.
- `_modified.csv` files are generated by a committed script (`scripts/build_modified_csv/build_moreno2023_s3_modified_csv.py`) that parses the per-row Description field to extract locus_tag, log2FC, and t-test columns under stable names — the paperconfig references the modified versions, not the originals.
- Significance: t-test p-values (`t_test_*` columns) with `pvalue_threshold: 0.05`. Directionality is encoded in log2FC sign of glucose-vs-control; `growth_phase: stationary` across all analyses.
- SS120 Marinobacter FC column missing for `G mM/L` comparison → correctly omitted from paperconfig.
