# Kaur 2018

**Title:** 100 Days of marine *Synechococcus*-*Ruegeria pomeroyi* interaction: A detailed analysis of the exoproteome

## Classification

**Bucket D -- defer / nothing to do (already integrated)**

100-day exoproteome time course of Synechococcus WH7803 + Ruegeria pomeroyi DSS-3 (both organisms now deployed in the KG, both as primary genome strains). Tables S4a/S4b (Syn) and S5a/S5b (Ruegeria) provide per-protein signed-ratio fold change vs day 1 plus q-values across 7 time points; all four are wired to time-course experiments in paperconfig.yaml and have `*_resolved.csv` siblings, so the per-time-point `Changes_expression_of` edges are produced on every build. The remaining supplementary tables (S1-S3) are sample/protein-id manifests with no DE; nothing actionable left.

**Authors:** Amandeep Kaur, Juan R. Hernandez-Fernaud, Maria del Mar Aguilo-Ferretjans, Elizabeth M. Wellington, Joseph A. Christie-Oleza

**Journal:** Environmental Microbiology (2018) 20(2), 785-799

**DOI:** 10.1111/1462-2920.14012

## Summary

Detailed exoproteomic time-course analysis of a 100-day co-culture between *Synechococcus* sp. WH7803 and *Ruegeria pomeroyi* DSS-3 in both natural oligotrophic seawater (SW) and nutrient-enriched artificial seawater (ASW). Samples collected at 8 timepoints (days 1, 3, 7, 14, 21, 32, 60, 100) from triplicate flasks. The study reveals a transition between initial growth phase and stable state, with the heterotroph switching motility based on organic matter availability and the phototroph adapting to seawater oligotrophy by reducing selective leakiness and increasing nutrient acquisition.

## Organisms

- **Synechococcus sp. WH7803** (NCBI taxid: 32051) -- NOT in KG yet
- **Ruegeria pomeroyi DSS-3** (NCBI taxid: 89184) -- NOT in KG yet

## Experimental Design

- **Omics type:** EXOPROTEOMICS (shotgun proteomics of cell-free exoproteome)
- **Instrument:** nanoLC-ESI-MS/MS, Orbitrap Fusion
- **Software:** MaxQuant 1.5.5.1 (identification), Perseus 1.5.5.3 (quantification/statistics)
- **Quantification:** Label-free, relative protein abundance normalized to protein size, log2 transformed
- **Statistical test:** Two-sample Student's t-test (P=0.05), permutation-based FDR (q=0.05)
- **Growth conditions:** 22C, 10 umol photons m-2 s-1 continuous light, 140 rpm shaking
- **Media:** Natural seawater (SW, Gulf Stream) and enriched artificial seawater (ASW)
- **Timepoints:** Days 1, 3, 7, 14, 21, 32, 60, 100
- **Replicates:** 3 biological replicates per condition per timepoint

## Supplementary Tables Used

| Table | Organism | Medium | CSV File |
|-------|----------|--------|----------|
| S4a | Syn WH7803 | SW | `table s4a Syn SW emi14012-sup-0004-suppinfo4.csv` |
| S4b | Syn WH7803 | ASW | `table s4b Syn ASW emi14012-sup-0004-suppinfo4.csv` |
| S5a | R. pomeroyi DSS-3 | SW | `table s5a Pom SW emi14012-sup-0005-suppinfo5.csv` |
| S5b | R. pomeroyi DSS-3 | ASW | `table s5b Pom ASW emi14012-sup-0005-suppinfo5.csv` |

## Fold Change Convention

Fold change values are **signed ratios** comparing each timepoint to timepoint 1 (day 1 baseline):
- Positive values (e.g., 3.8) = higher abundance at later timepoint vs day 1
- Negative values (e.g., -1.6) = lower abundance at later timepoint vs day 1
- Value of 1 or -1 = no change or not detected at one timepoint

This is the same signed-ratio convention used in Oleza 2017 and Ma 2022. No `fold_change_type` conversion is needed -- values pass through as-is.

## Gene ID Types

- **S4 tables (Synechococcus):** "NCBI reference" column contains RefSeq protein accessions with version (e.g., `YP_001223734.1`). "Other reference" column has SynWH7803_NNNN locus tags.
- **S5 tables (R. pomeroyi):** "NCBI reference" column contains NCBI protein accessions without version (e.g., `AAV93452`). No "Other reference" column.

## Phase 4 Requirements

Neither organism is currently in the KG. Integration will require:
1. Adding Synechococcus WH7803 and R. pomeroyi DSS-3 genomes to the genome registry
2. Building gene mappings for both organisms
3. Protein-level ID resolution (RefSeq YP_ accessions for Syn, GenBank AAV accessions for R. pomeroyi)

## CSV Structure Notes

- Row 1: Table title
- Row 2: Blank
- Row 3: Column headers (data header)
- Row 4+: Data
- `skip_rows: 2` in paperconfig to skip title and blank rows
- S4 tables have 25 columns (including "Other reference"); S5 tables have 24 columns (no "Other reference")
- Column names for TP60 have a double space ("  Time point 60 vs 1") in all 4 CSVs -- preserved exactly in paperconfig
