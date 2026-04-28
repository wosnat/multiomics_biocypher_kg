# Fadeev 2022

**Citation:** Fadeev E, et al. Characterization of membrane vesicles in *Alteromonas macleodii* indicates potential roles in their copiotrophic lifestyle. microLife (Oxford / FEMS), 2022/2023, uqac025.
**DOI:** 10.1093/femsml/uqac025 (filename `uqac025.pdf`)
**Organism(s):** *Alteromonas macleodii* strains AD45, ATCC27126, BGP6, BS11, HOT1A3, MIT1002 (six strains compared)
**Topic:** Comparative proteomic survey of outer membrane vesicles (MVs) vs whole-cell proteomes across six *Alteromonas macleodii* strains. For each strain, a Top-N most-abundant-MV-protein table reports `Prop. Abund. MVs` and `Prop. Abund. Cells` (NSAF-style proportional abundance percentages) with NCBI_PGAP function + accession and COG20 category. Single growth condition per strain; no treatment-vs-control fold-change.

## Classification

**Bucket B — new metrics / DE / resolution (want to add) — partial**

Per-strain "Prop. Abund. MVs" / "Prop. Abund. Cells" tables are single-condition NSAF-style abundance (no fold change). They fit the numeric `derived_metrics_table` pattern shipped 2026-04-21 (see `docs/kg-changes/non-de-evidence-extension.md`) — analogous to the biller 2022 vesicle integration. **HOT1A3 and MIT1002 are deployed**, so those two strains can be wired now. The other four strains (AD45, ATCC27126, BGP6, BS11) are not in the KG; either deploy them as new strains or skip those four CSVs. Action: paperconfig with two `derived_metrics_table` entries (HOT1A3, MIT1002) keyed on `WP_*` accessions, optional follow-up for the other strains.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `uqac025.pdf` | PDF | Main paper | reference | -- |
| `Figure_S1_Supplementary_Data.pdf` | PDF | Supplementary figure S1 | reference | -- |
| `Table_S1-most_abundant_MV_proteins_Supplementary_Data.xlsx` | XLSX | Combined six-sheet workbook (one sheet per strain) | reference | Source for the per-strain CSVs below |
| `Table_S1-...AD45.csv` | CSV | AD45 top MV proteins: WP_ accession + function + COG + `Prop. Abund. MVs %` + `Prop. Abund. Cells %` | blocked (strain not deployed) | Deploy AD45 first, then wire as `derived_metrics_table` (numeric) |
| `Table_S1-...ATCC27126.csv` | CSV | Same shape, ATCC27126 | blocked (strain not deployed) | Deploy ATCC27126 first, then wire as `derived_metrics_table` (numeric) |
| `Table_S1-...BGP6.csv` | CSV | Same shape, BGP6 | blocked (strain not deployed) | Deploy BGP6 first, then wire as `derived_metrics_table` (numeric) |
| `Table_S1-...BS11.csv` | CSV | Same shape, BS11 | blocked (strain not deployed) | Deploy BS11 first, then wire as `derived_metrics_table` (numeric) |
| `Table_S1-...HOT1A3.csv` | CSV | Same shape, HOT1A3 | ready to wire | Add as `derived_metrics_table` with `value_kind: numeric` (prop_abund_MVs, prop_abund_cells) |
| `Table_S1-...MIT1002.csv` | CSV | Same shape, MIT1002 | ready to wire | Add as `derived_metrics_table` with `value_kind: numeric` (prop_abund_MVs, prop_abund_cells) |

## Current paperconfig summary

No `paperconfig.yaml` -- paper is not currently wired into the KG (also not present in `paperconfig_files.txt`).

## Notes

- Gene IDs are NCBI WP_ protein accessions (`NCBI_PGAP_accession` column). For MIT1002/HOT1A3 the standard `protein_id_refseq` Tier-2 path on the existing gene_id_mapping should resolve them via the `derived_metrics_table` loader (shipped 2026-04-21, see `docs/kg-changes/non-de-evidence-extension.md`).
- `Prop. Abund.` columns are formatted as percent strings (e.g. `13.3%`) -- need stripping/normalisation before any numeric ingestion.
- For HOT1A3 and MIT1002 the COG category is already present in the source CSV; the KG already carries CogFunctionalCategory edges per gene from eggNOG, so the COG column is informational only.
- If the MV-vs-cell proportion ratio is the analytically interesting signal, an `MV_enrichment_ratio = log2(MVs / Cells)` could be precomputed and stored as a single per-gene numeric DerivedMetric instead of two.
- Adding AD45, ATCC27126, BGP6, BS11 to the genome registry would be a non-trivial four-strain prepare_data run before any of those tables could be ingested.
