# Batch D — control-timing curation

## Summary
- Papers: 6 · DE experiments: 14 · baseline: 1 · mixed: 0 · time_matched: 13 · shared-model ids: 0
- Needs human review:
  - Al-Hosani 2015 `salt_saltacclimated_5_salt_as9601_rnaseq`: arms matched by growth state (exponential) but sampled on different days (control d2+d4, salt d4+d6); set `time_matched` per the "same growth state" rule, with a design note. Confirm that reading.
  - Aharonovich 2016 HOT1A3 experiments: `logfc_col` holds linear ratios, not log2 (see Other findings). Not a timing issue, but it corrupts the edges.
  - Fang 2019: `treatment_organism: Phage` on a vDOM (cell-free lysate) addition (see Other findings).

## Aharonovich 2016
- Sources read: main PDF pp. 1–5 (Methods, Fig. 1, Fig. 2 legend, Results); README; headers of all 5 CSVs. Supplementary information (sampling details) is not in the directory.
- Paper default: `time_matched`
- Sampling: triplicate cocultures and axenic controls set up in parallel. All 12 libraries are from one time point: "Twelve mRNA libraries from one timepoint (20 h)" (p. 2893). Fig. 2 legend: "The samples for t = 20 h were selected for transcriptomic analysis."
- DE fitting: Rockhopper, one pairwise comparison per table (coculture vs axenic for each strain, plus two coculture-vs-coculture comparisons on the HOT1A3 side). No shared model is described, so no id.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `coculture_alteromonas_hot1a3_med4_rnaseq` | 20h | time_matched (default) | — | high | p. 2893 "mRNA libraries from one timepoint (20 h)"; Fig. 2a "Each of these co-cultures was compared with a control of the relevant Prochlorococcus strain growing axenically" |
| `coculture_alteromonas_hot1a3_mit9313_rnaseq_inoc_5e6` | 20h | time_matched (default) | — | high | same |
| `coculture_alteromonas_hot1a3_mit9313_rnaseq_inoc_05e6` | 20h | time_matched (default) | — | high | same |
| `coculture_prochlorococcus_mit9313_high_alteromonas_hot1a3_rnaseq` | 20h | time_matched (default) | — | high | 9313+1A3 vs 9313+diluted 1A3, both 20 h cocultures (CSV header "9313/9313DIL") |
| `coculture_prochlorococcus_mit9313_hot1a3_rnaseq` | 20h | time_matched (default) | — | high | 9313+1A3 vs MED4+1A3, both 20 h (CSV header "9313+1A3/med4+1A3") |

Opens / uncertainties:
- None for timing. Rockhopper can process several conditions in one run, but the paper does not say it did, and each table is a single pairwise comparison, so no id.

Other findings (not fixed):
- **HOT1A3 fold changes are linear, not log2.** `fold change 9313/9313DIL` (min 0.002, median 3.7, max 282,301, no negatives) and `fold change 9313+1A3/med4+1A3` (min 0.0002, median 1.2, max 4,801, no negatives) are plain ratios. Neither analysis sets `fold_change_type: linear`, so the adapter loads them as `log2_fold_change`. Every down-regulated HOT1A3 gene shows up as a positive value.
- **The HOT1A3 control conditions are wrong.** `coculture_prochlorococcus_mit9313_high_alteromonas_hot1a3_rnaseq` and `coculture_prochlorococcus_mit9313_hot1a3_rnaseq` say `control_condition: Pro99 medium growth conditions`. The actual denominators are "HOT1A3 + MIT9313 diluted (0.5e6)" and "HOT1A3 + MED4" cocultures; there is no axenic HOT1A3 arm. The same goes for their `treatment_organism`/role framing: the variable is the partner strain or density, not coculture vs axenic.
- In `de_genes_med4_t20h.csv` the two expression column headers ("med4 t=20" and "med4+1a3 t=20") look swapped. The log2FC sign disagrees with them in 96% of rows, but it matches the paper (thioredoxin TxlA up, recR/recO down in MED4+1A3). The `logfc_col` is correct; only the expression headers are mislabelled.
- The paper's DE criterion is fold change > 2, q < 0.05 and expression > 10. The config uses `pvalue_threshold: 0.01` with no logFC threshold.
- Methods give 20 °C and 27 µE for maintenance, but the experiments use 24C and 10 µE. These may be the experiment conditions from the Supplementary information, which was not available to check.

## Al-Hosani 2015
- Sources read: main PDF pp. 2–3 (Methods in full, Fig. 1 and legend, Results p. 23); README; CSV header.
- Paper default: `time_matched`
- Sampling: cultures were acclimated for 3 transfers, then grown in 300 ml batches with samples drawn every 48 h. RNA-seq used 3 replicates "sampled on 2 consecutive days of exponential growth". Fig. 1 marks the control sampling days (*) as days 2 and 4 and the treatment days (X) as days 4 and 6.
- DE fitting: DESeq, a single two-group comparison with the control (3.8% salt) as reference.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `salt_saltacclimated_5_salt_as9601_rnaseq` | acclimated | time_matched (default) + design_notes | — | medium | Fig. 1 legend: "The * and X indicate the sampling days used for RNAseq analysis for control and treatment respectively"; p. 23 "sampled on 2 consecutive days of exponential growth" |

Opens / uncertainties:
- `salt_saltacclimated_5_salt_as9601_rnaseq`: the two arms were harvested at the same growth state (exponential, verified by the paper's within-group QC) but 2 days apart in calendar time. I chose `time_matched` under the "acclimated cultures harvested at the same growth state" rule and recorded the day offset in `design_notes`. If calendar offset should count, the alternative is `baseline`.

Other findings (not fixed):
- The paper calls DE at raw p < 0.05 ("using the p-value as criterion"), not adjusted p. The config uses `padj` with `prefiltered: true`, which is fine for edges, but significance in the source table is p-value-based.

## Anjur-Dietrich 2025
- Sources read: main PDF (Results "Adherent and planktonic cells are transcriptionally distinct", Methods "RNA extraction and analysis", lines 563–580); supp docx (media-1) legends; media-2 Supplementary Table 5 header; CSV header.
- Paper default: `time_matched`
- Sampling: "Each population was grown to mid-exponential phase in an undisturbed state before collection". Planktonic cells were fixed in suspension and biofilm cells were fixed while adhered. Triplicates.
- DE fitting: edgeR, one planktonic vs biofilm comparison (|logFC| > 0.8, q < 0.05).

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `growth_state_freeliving_planktonic_enrichment_mit9301_rnaseq` | mid-exponential | time_matched (default) | — | medium | Results: "Each population was grown to mid-exponential phase in an undisturbed state before collection" — matched growth state; same calendar day not stated |

Opens / uncertainties:
- The paper does not say whether both enrichments were harvested on the same day; both were harvested at mid-exponential, so `time_matched` by growth state.

Other findings (not fixed):
- `treatment_type: [growth_phase]` is a poor fit. The contrast is lifestyle (planktonic vs biofilm enrichment, both mid-exponential), not growth phase.
- media-2 Supplementary Table 5 is the same 240-row significant-only table as the integrated CSV, so the README's "inspect" action can be closed.

## Bagby & Chisholm 2015
- Sources read: main PDF Methods ("Transcriptional profiling", microarray analysis), Fig. 3/4/6 legends, Results; Table S1 xls notes (moesm33); supp doc (moesm29) figure legends; CSV header.
- Paper default: `time_matched`
- Sampling: three mid-log cultures were pelleted, pooled and split evenly into 4 bottles pre-equilibrated with air or one of 3 shock gases. Sampling at t = 0 (immediately after the split) and then in parallel at the same times. Arrays were run at t = 0, 70 min, 6 h, 12 h and 24 h, with 2 replicates.
- DE fitting: betr Bayesian time-series. Fig. legends: "the baseline used for significance testing (for the air control, expression at t = 0; for experimental conditions, the parallel time course in air)". Table S1 note 3: "For 0.036% CO2 / 21% O2 (air; control) maximal fold change is with respect to t = 0 transcript abundance. For other (experimental) conditions, maximal change is the extremum of log2(fold change(t - t0), gas) - log2(fold change(t - t0), air)." Each condition is tested separately (gas vs air, air vs t0), so no shared model id.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `carbon_air_0036_co2_21_med4_microarray` | max across 3-24h | baseline (override) | — | high | Table S1 note 3: air "maximal fold change is with respect to t = 0" |
| `carbon_ambient_co2_no_o2_med4_microarray` | max across 3-24h | time_matched (default) + design_notes | — | high | Table S1 note 3; legends "for experimental conditions, the parallel time course in air" |
| `carbon_low_co2_0004_co2_med4_microarray` | max across 3-24h | time_matched (default) + design_notes | — | high | same |
| `carbon_low_co2_and_no_med4_microarray` | max across 3-24h | time_matched (default) + design_notes | — | high | same |

This agrees with the plan's prior. The three shock experiments got an optional `design_notes` because their value is a difference of t0-normalised changes (gas minus air), not a direct gas/air ratio.

Opens / uncertainties:
- None.

Other findings (not fixed):
- The `experimental_context` on the 3 shock experiments lists timepoints "3h, 6h, 9h, 12h, 18h, 24h", but samples were taken at 0, 70 min, 3, 6, 9, 12, 18 and 24 h and only 0, 70 min, 6, 12 and 24 h were hybridised. The label "max across 3-24h" may also be off: the extremum could include 70 min.
- `growth_phase: acute_stress` on `gas_shock_air` (the air control) is debatable. It is the handling-stress control, which arguably fits, but it is not a gas shock.

## Fang 2019
- Sources read: main PDF (Results "Prochlorococcus growth is not affected by viral lysis products" and "An overview", Fig. 1 legend, Methods DESeq2 section); `legend of supp tables.txt`; Table S2 header (skip_rows preamble); supp Table S1.
- Paper default: `time_matched`
- Sampling: vDOM was added to mid-log axenic MIT9313 at 1:4 v/v, and controls were amended with Pro99. Transcriptomes were taken "at 0.5, 1, 2, 4, 8, 12, 24, 48 and 72 h after amendments with vDOM or the Pro99 growth medium (control)", with 2 replicates per time point.
- DE fitting: DESeq2, comparing vDOM-amended with control samples at each time point. Table S2 is "log2 ratio (vDOM/Control)". All 9 contrasts sit in one experiment, so no id is needed.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `viral_viral_dom_vdom_addition_mit9313_rnaseq` | 0.5h … 72h (9) | time_matched (default) | — | high | Results: "cells taken at 0.5, 1, 2, 4, 8, 12, 24, 48 and 72 h after amendments with vDOM or the Pro99 growth medium (control)"; supp legend "normalized to those of the control cultures" |

Opens / uncertainties:
- None.

Other findings (not fixed):
- `treatment_organism: Phage` (taxid 10239) on a cell-free viral-lysate (vDOM) addition creates a `Tests_coculture_with` treatment edge to Phage, but no live phage was present (S-ShM2 lysate of WH8102, filtered). `background_factors: [light]` also lacks `axenic`, though the cultures were axenic.
- `growth_phase: exponential` at 48 h and 72 h is consistent with the paper ("continued exponential growth for up to 3 days").

## He 2022
- Sources read: main PDF Methods in full ("Low salinity acclimation", "RNAseq analysis"); Data Sheet 1 docx (DE tables only); `geo_source.txt`; GEO DE and FPKM table headers.
- Paper default: `time_matched`
- Sampling: after 5 acclimation transfers, acclimated cultures (28 psu and 34 psu control, 3 replicates each) were inoculated into fresh medium. "During the exponential growth phase, 100 ml cultures were filtered". The paper does not give the harvest day for each arm.
- DE fitting: DESeq per strain according to Methods. The GEO tables carry the edgeR logCPM signature. There is one two-group comparison per strain, and the strains are different organisms, so no shared id.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `salt_low_salinity_acclimation_28_natl1a_rnaseq` | acclimated | time_matched (default) | — | medium | Methods: "with salinity 34 psu being the control ... During the exponential growth phase, 100 ml cultures were filtered" |
| `salt_low_salinity_acclimation_28_med4_rnaseq` | acclimated | time_matched (default) | — | medium | same |

Opens / uncertainties:
- Both strains: control and treatment were harvested in exponential phase, but the day is not stated, so this is `time_matched` by matched growth state. No design note was added; one could be added if growth-state matching should always be recorded.

Other findings (not fixed):
- The test-tool mismatch (Methods say DESeq, GEO files carry edgeR output) is already recorded in the README.
