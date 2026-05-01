# Zimmerman 2023

**Citation:** Zimmerman AE, Podowski JC, Gallagher GE, Coleman ML, Waldbauer JR. "Tracking nitrogen allocation to proteome biosynthesis in a marine microbial community." *Nature Microbiology* 8:498–509 (March 2023). Published online 2023-01-12.
**DOI:** 10.1038/s41564-022-01303-9
**Organism(s):** Surface mixed-layer microbial community at Station ALOHA (Hawaii Ocean Time-series, North Pacific Subtropical Gyre, 25 m). Most peptides assigned to *Prochlorococcus* (genus 1218), *Pelagibacter* (SAR11), *Synechococcus*, Pelagibacteraceae, Rhodobacteraceae, Rhodospirillaceae, Flavobacteriia, Isochrysidales, Pelagomonadales, Dinophyceae, Archaea — community proteomics, not a single cultured isolate.
**Topic:** ¹⁵N-tracking metaproteomics. On-deck bottle incubations of ALOHA surface seawater amended (2.5 µM) with five ¹⁵N-labelled N substrates — ammonium, glutamate, leucine, nitrate, urea — sampled at 2/4/8/11 h across day/dusk and night/dawn incubations. Combines two LC-MS readouts per peptide: (i) **cPIE** — peptide-level ¹⁵N atom% incorporation rate (atom%/h), measuring active biosynthesis; (ii) **diDO-IPTL** — peptide abundance change (log2/h), measuring net protein expression. Together these resolve N substrate preference and proteome-N allocation per taxon. Headline findings: *Prochlorococcus* dominates ¹⁵N-ammonium and ¹⁵N-urea uptake; allocates daytime N to chaperones / cell division / amino-acid metabolism and nighttime N to translation / ATP synthase / ribosomes / carbon fixation; chaperones (notably GroEL/S) are the single largest biosynthetic N sink across both day and night.

## Classification

**Not yet integrated.** No `paperconfig.yaml` in this directory. This paper sits at the boundary of what the KG schema currently supports — see "Integration considerations" below.

## Available data inventory

| File | Type | Content |
|---|---|---|
| `s41564-022-01303-9.pdf` | PDF | Main paper (12 pages, includes methods) |
| `41564_2022_1303_MOESM1_ESM.pdf` | PDF | Supplementary text + Extended Data figures |
| `41564_2022_1303_MOESM3_ESM.txt` | TSV | **Supplementary Data 1** — 10,016 rows × 25 columns. One row per (peptide × experiment × incubation × N amendment). Key columns: `peptide` (sequence), `expt` (Expt1/Expt2), `incubation` (day_dusk/night_dawn), `amendment` (ammo/glu/leuc/nitr/urea), `IPTL.log2.slope` (diDO-IPTL log2/h), `corrected.15N.slope` (cPIE atom%/h, NA where unmeasured), `n_possible_proteins`, `n_distinct_taxID`, `top1_taxID` / `top1_name` / `top1_level` (LCA taxonomy of best parent protein), `top1_species`/`genus`/`family`/`order`/`class`/`phylum`/`kingdom`/`superkingdom`, `function_category` / `function_subcategory` / `function_name` / `function_ID` (KEGG KO) / `function_description` |

Per-peptide breakdown by genus in MOESM3 (long-format rows, so each peptide appears up to ~10× across amendments × incubations × experiments):

| Genus | Rows |
|---|---:|
| *Prochlorococcus* | 3,724 |
| *Candidatus* Pelagibacter | 1,912 |
| (no genus assigned) | 2,322 + 485 NA |
| *Emiliania* | 275 |
| *Aureococcus* | 254 |
| *Candidatus* Thiodiazotropha | 106 |
| *Symbiodinium* | 75 |
| *Candidatus* Puniceispirillum | 65 |

Within *Prochlorococcus* most peptides are assigned at genus (taxid 1218, 2,105 rows) or species (taxid 1219 *P. marinus*, 1,333 rows) level — none reach a specific cultured strain in our KG.

## Integration considerations

Three structural problems prevent direct paperconfig.yaml integration under the current schema:

1. **Peptide-level, not gene-level.** All values are peptide measurements. The KG's `Changes_expression_of` and `Derived_metric_*` edges target Gene nodes via locus_tag. There is no peptide → gene mapping in the supplementary data; the authors provide LCA taxonomy + KEGG KO per peptide, but no `protein_id` / `WP_*` / locus_tag column. Going from peptide → gene would require re-running the peptide search against a strain-specific FASTA (MED4, MIT9301, etc.) and accepting only single-genome-unambiguous hits — substantial reanalysis, not data import.

2. **Community reference DB, not a cultured strain.** Peptides were searched against the **ALOHA Gene Catalogue v1** (a non-redundant gene set built from ALOHA metagenomes/metatranscriptomes; ~4,565 peptides matched), supplemented by metatranscriptome-derived translations. This is analogous to the Marinobacter / Alt_MarRef situation (see `docs/community_proteomics_marref_saga.md`) but worse — ALOHA Gene Catalogue spans hundreds of taxa, not a single reference proteome we could register as a `reference_proteome_match` OrganismTaxon.

3. **No DE comparison.** The paper's analytical unit is *rate per peptide* (¹⁵N atom%/h, IPTL log2/h), not treatment-vs-control fold change. The ¹⁵N atom% rate is naturally a numeric `Derived_metric_quantifies_gene` (Bucket B). The IPTL log2-slope could in principle map to fold-change-like edges, but there is no contrasted control condition — every amendment is its own track. `compartment` would be `whole_cell`, `treatment_type` likely `nitrogen` + `diel`.

If/when integrated, the natural shape would be: **6 Experiments per amendment** (5 amendments × 2 incubations = 10, less the few with no signal) emitting numeric DerivedMetrics for `15n_incorporation_rate_atom_pct_per_h` and `iptl_log2_slope_per_h`, scoped to *Prochlorococcus* (and possibly *Pelagibacter*, *Synechococcus*) at the genus aggregate via a `reference_proteome_match` organism — pending a peptide → strain mapping pipeline.

## Notes

- Mass spectra: ProteomeXchange MSV000089118 / PRIDE PXD038614. FASTAs: Mendeley Data `data.mendeley.com/datasets/7226j6cpp3/1`.
- ALOHA Gene Catalogue: `https://datacommons.cyverse.org/browse/iplant/home/shared/imicrobe/projects/263/`
- cPIE pipeline code: `github.com/waldbauerlab`
- Two follow-up readouts not in MOESM3: protein-level rollups (Fig. 3–5 derived) and proteome-N-allocation per functional category (Fig. 4b/5). Reproducing those from MOESM3 requires the abundance-weighting in Methods §"Calculation of N allocation to proteome functional categories" (multiplies median ¹⁵N rate per category × median IPTL spectral counts).
- Peer-reviewed Nature Microbiology paper (not preprint).
