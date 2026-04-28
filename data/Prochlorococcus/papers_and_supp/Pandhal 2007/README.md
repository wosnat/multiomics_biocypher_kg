# Pandhal 2007

**Citation:** Pandhal J, Wright PC, Biggs CA. A quantitative proteomic analysis of light adaptation in a globally significant marine cyanobacterium *Prochlorococcus marinus* MED4. Journal of Proteome Research, 2007, 6(3): 996–1005.
**DOI:** 10.1021/pr060460c
**Organism(s):** *Prochlorococcus marinus* MED4
**Topic:** Quantitative iTRAQ proteomics of MED4 grown under three light intensities (low 20, medium 60, high 100 μeinstein m⁻² s⁻¹). Two iTRAQ experiments (biological + technical replicates): iTRAQ1 = HL/LL/ML/HL; iTRAQ2 = LL/ML/HL/HL. ~184 proteins identified (~11% of predicted coding genes); 15 proteins reported as statistically significantly differentially expressed between light intensities (down-regulation of photosystem proteins, up-regulation of GroEL). Proteins annotated with COG categories.

## Classification

**Bucket B — new metrics / DE / resolution (want to add)**

Per-gene proteomic DE across light intensities is exactly the bucket-B shape: pairwise iTRAQ ratios (HL vs LL, ML vs LL) with weighted error-factor probabilistic p-values, on ~15 significantly DE proteins out of ~184 quantified. Action: convert the two `pr060460csi*.xls` files to CSV, build a `paperconfig.yaml` with `omics_type: PROTEOMICS`, `treatment_type: ["light"]`, one Experiment per pairwise light comparison, and feed the smaller XLS (the 15-protein DE table) as `csv` with `statistical_analyses`. The larger XLS (full ~184-protein quantitation with emPAI / pI / MW) is bucket-C-shaped (single-condition descriptors); skip or hold for DerivedMetric. MED4 is already deployed and has full ID coverage so resolution should be straightforward. Strain status: no new deployment required.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `a-quantitative-proteomic-analysis-of-light-adaptation-in-a-globally-significant-marine-cyanobacterium-prochlorococcus.pdf` | PDF | Main article | reference | — |
| `pr060460csi20061210_124449.xls` | XLS | Larger supplementary table — most likely full list of identified proteins across both iTRAQ experiments (~184 proteins with emPAI, iTRAQ ratios, COG assignments, pI, MW, peptide counts) | add | Export to CSV; candidate for `derived_metrics_table` (numeric `emPAI`, `pI`, `MW`, `peptide_count`, boolean `detected_in_iTRAQ1`/`detected_in_iTRAQ2`) + `csv` with `statistical_analyses` if the iTRAQ ratios + weighted p-values are per-gene DE-shaped |
| `pr060460csi20061210_124513.xls` | XLS | Smaller supplementary table — likely the 15 significantly differentially expressed proteins with fold-change and p-value across light conditions | add | Export to CSV and add as `csv` with `statistical_analyses` (HL vs LL, ML vs LL — 2–3 pairwise comparisons, `omics_type: PROTEOMICS`, `test_type: "iTRAQ weighted error-factor probabilistic p-value"`) |

## Current paperconfig summary

No paperconfig.yaml — paper is not integrated.

## Recommended actions

1. **Add** — Export both XLS files to CSV and inspect column structure. The smaller table (`124513.xls`) is the priority integration — it likely contains per-gene fold-change + p-value across light conditions, matching the `statistical_analyses` schema cleanly.
2. **Add** — The larger table (`124449.xls`) can seed a `derived_metrics_table` with numeric metrics (`emPAI_hl`, `emPAI_ml`, `emPAI_ll`, `MW`, `pI`, `unique_peptides`) + booleans (`detected_in_iTRAQ1`, `detected_in_iTRAQ2`, `cytoplasmic_membrane_predicted`). COG assignments are already covered by eggNOG so don't duplicate.
3. **Create paperconfig.yaml** — Experiment(s): one per pairwise comparison (HL vs LL, ML vs LL, HL vs ML). Organism: `Prochlorococcus MED4`. `omics_type: PROTEOMICS`. `treatment_type: ["light"]`. `light_condition: "high light (100 μE)" / "medium light (60 μE)" / "low light (20 μE)"`.
4. **Skip** — No need to duplicate the XLS files as CSV-XLSX pairs; keep only the CSV derivative for integration.
5. **Check gene IDs** — Column in the XLS is likely protein locus tag (`PMMxxxx`) or UniProt accession. Use `/check-gene-ids` after CSV export; MED4 has excellent ID coverage in `gene_id_mapping.json` so match rates should be high.

## Notes

- MED4 is the KG's primary focus strain — DE edges from this paper will complement the existing Tolonen 2006, Thompson 2011, Coe 2024, Lindell 2005 MED4 proteomics / transcriptomics layers.
- Paper Table 3 (in-PDF) lists COG distribution — not useful for integration since eggNOG already provides this.
- Paper identifies 15 DE proteins out of 184 quantified — table scope for the DE analyses would be `filtered_subset` (the smaller XLS) or `all_detected_genes` (if pulling from the full-184 table with significance flags). Clarify by reading the XLS headers.
- 2007-vintage protein IDs may include obsolete GI numbers; Tier-2 `protein_id` + Tier-3 `gene_name` lookups will handle them.
