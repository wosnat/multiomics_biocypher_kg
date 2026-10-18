# Batch H — control-timing curation

## Summary
- Papers: 7 · DE experiments: 23 · baseline: 5 · mixed: 0 · time_matched: 18 · shared-model ids: 0
- Needs human review:
  - Beliaev 2014 (all 6): steady states were sampled one after another in the same reactor (axenic, then lactate coculture, then HCO3- coculture), not in parallel. Marked `time_matched` because the arms were harvested in the same growth state (carbon-limited steady state). Confirm that this is how sequential chemostat steady states should be classed.
  - Huang 2020 `viral_s_sbp1_infection_wh7803_rnaseq`: marked `baseline`. Every infected timepoint (15 min to 7h) is compared with one uninfected control sampled at 30 min. The 15 min and 1h contrasts are close to 30 min, but none is at the same time.
  - Doron 2016 (all 3): marked `time_matched` (the supplementary methods state that uninfected controls were taken at every timepoint). The 2-5 min control (column `C0`) may have been sampled just before infection. This is recorded in `design_notes`.
  - Ma 2022 `coculture_ecoli_vs_axenic_utex2973_rnaseq`: the RNA sampling day is not stated. `time_matched` is inferred from the word "meanwhile" (low confidence).
  - Oleza 2017: `omics_type: EXOPROTEOMICS` / `compartment: exoproteome` contradicts the methods, which describe cell-pellet (cellular) proteomes. See Other findings.

## Beliaev 2014
- Sources read: main PDF (Methods: Cultivation, Chemostat cultivation, RNA isolation and sequencing; Table 1; Results text introducing the comparisons); SI docs moesm41/42 (supplementary figures and stoichiometry, nothing on sampling); README; CSV headers of S3 (RPKM columns for axenic, lactate coculture and HCO3- coculture, plus 3 fold-change columns).
- Paper default: `time_matched`
- Sampling: steady-state chemostats. The axenic control was a carbon-limited chemostat. The coculture was set up by adding Shewanella to the Synechococcus chemostat and switching it gradually to lactate, then back to HCO3-. "Samples for all the analyses were only taken after at least five volume changes under steady-state conditions." RNA came from duplicate chemostat cultures.
- DE fitting: no statistical test. RPKM ratios with a 2-fold cut-off, one ratio per pair of conditions. No shared model.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `coculture_lactate_vs_axenic_syn7002_rnaseq` | steady state | time_matched (default) | — | medium | p. 3 Methods: "Samples for all the analyses were only taken after at least five volume changes under steady-state conditions" |
| `coculture_lactate_vs_axenic_shew_w3181_rnaseq` | steady state | time_matched (default) | — | medium | same |
| `coculture_hco3_vs_axenic_syn7002_rnaseq` | steady state | time_matched (default) | — | medium | same |
| `coculture_hco3_vs_axenic_shew_w3181_rnaseq` | steady state | time_matched (default) | — | medium | same |
| `coculture_hco3_vs_lactate_syn7002_rnaseq` | steady state | time_matched (default) | — | medium | same |
| `coculture_hco3_vs_lactate_shew_w3181_rnaseq` | steady state | time_matched (default) | — | medium | same |

Opens / uncertainties:
- All 6: the steady states were sequential, not parallel, so they differ in calendar time but are matched in physiological state. I chose `time_matched` under the "acclimated cultures harvested at the same growth state" rule. I added no `design_notes` because the field is optional for `time_matched`. If you want this caveat to be visible in the KG, add a one-line note.

Other findings (not fixed):
- None of note. The test_type `RPKM_fold_change` (no p-values) is consistent with the paper.

## Kratzl 2024
- Sources read: main PDF (Results paragraph on the reference experiment, Fig. 3 legend, Materials and methods "OMICs analysis" and "Multi OMICs methods"); SI PDF (grep for statistics); paperconfig.
- Paper default: `time_matched`
- Sampling: a 9-fold parallel photobioreactor with three set-ups (axenic S. elongatus, axenic P. putida with a sucrose feed, and the coculture), each in biological triplicate. "After ~60 h samples for multi-OMICs analysis were taken". Transcriptomics and proteomics both came from experiment EII: "all samples were taken from the same experiment conducted in the 9-fold" reactor.
- DE fitting: the RNA-seq analysis was outsourced to Eurofins and the methods give no fitting detail. Each experiment is a different organism × omics pair, so a shared id is impossible.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `coculture_vs_pputida_axenic_rnaseq` | 60h | time_matched (default) | — | high | p. 12 Methods: "After ~60 h samples for multi-OMICs analysis were taken"; p. 4: "all samples were taken from the same experiment" |
| `coculture_vs_selongatus_axenic_rnaseq` | 60h | time_matched (default) | — | high | same |
| `coculture_vs_pputida_axenic_proteomics` | 60h | time_matched (default) | — | high | same |
| `coculture_vs_selongatus_axenic_proteomics` | 60h | time_matched (default) | — | high | same |

Opens / uncertainties: none.

Other findings (not fixed):
- `temperature` reads "25-32C" on the RNA experiments and "27.0-27.4C" on the proteomics experiments, but both omics came from the same run (EII, 27.0 °C per p. 4).
- `light_intensity` units are written "umol photons m-2 s-2", copying a typo in the paper (it should be s-1).

## Ma 2022
- Sources read: main PDF (Methods: transcriptomic analysis, quantitative proteomics, metabolomics; Results; supplementary figure legend about the proteome); paperconfig.
- Paper default: `time_matched`
- Sampling: proteomics used "the same weight of 4-day cocultured strains … Axenic Synechococcus 2973 cells with the same incubation time were used as a control." For RNA-seq, coculture cells were collected "meanwhile, the same amount of Synechococcus cscB+ cultivated under axenic conditions was used as a control". No day is given for the RNA-seq.
- DE fitting: RNA-seq was run by GENEWIZ (FC > 1.5, padj ≤ 0.05; tool not stated). Proteomics was run by BGI (iTRAQ, t-test of nine comparison groups). These are different omics, so no shared id.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `coculture_ecoli_vs_axenic_utex2973_rnaseq` | 4 days | time_matched (default) + design_notes | — | low | p. 4 Methods: "meanwhile, the same amount of Synechococcus cscB+ cultivated under axenic conditions was used as a control" |
| `coculture_ecoli_vs_axenic_utex2973_proteomics` | 4 days | time_matched (default) | — | high | p. 4 Methods: "Axenic Synechococcus 2973 cells with the same incubation time were used as a control" |

Opens / uncertainties:
- RNA-seq: the sampling day is not stated, and the "4 days" in the paperconfig was carried over from the proteomics. Checking the GEO record or the SI docx would settle it.

Other findings (not fixed):
- `growth_phase: acclimated_steady_state` on all three analyses. This was a 4-day batch coculture in a dialysis bag, not a chemostat or turbidostat, so `exponential` (or `unknown`) is probably more accurate.
- The RNA-seq methods also mention "ascorbic acid-treated axenic" cells collected alongside. The supplementary tables are coculture vs plain axenic, so nothing is mis-wired, but the extra arm is not represented.

## Oleza 2017 (Christie-Oleza 2017)
- Sources read: main PDF (Methods: Bacterial growth and culture conditions, Preparation of cellular proteomes for nanoLC–MS/MS; Results on ASW and seawater proteomes; Fig. 2 legend; Table 1); SI PDF converted to text (grep only); paperconfig.
- Paper default: `time_matched`
- Sampling: proteomes of Synechococcus "grown in ASW and seawater medium in the presence and absence of R. pomeroyi" and of R. pomeroyi with and without the phototroph in seawater. "At the time point of collection (35 days for ASW and 10 days for seawater)", with triplicate cultures.
- DE fitting: ASW used a t-test on NSAF values, seawater Synechococcus used ANOVA (Progenesis), and R. pomeroyi used a Perseus t-test. These are separate per-comparison tests, so no shared id.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `early_coculture_syn_wh7803_proteomics` | 35 days | time_matched (default) | — | high | p. 8 Methods: "At the time point of collection (35 days for ASW and 10 days for seawater)" |
| `late_coculture_syn_wh7803_proteomics` | 10 days | time_matched (default) | — | high | same |
| `longterm_coculture_rpom_dss3_proteomics` | long-term | time_matched (default) | — | medium | same sentence (seawater = 10 days); the R. pomeroyi monoculture is not separately dated |

Opens / uncertainties: none on timing.

Other findings (not fixed):
- **Compartment is wrong.** All 3 experiments carry `omics_type: EXOPROTEOMICS` / `compartment: exoproteome`, but the methods describe cellular proteomes ("Shotgun proteomic analysis of cellular extracts … Cell pellets were collected by centrifugation"). They should be `PROTEOMICS` / `whole_cell`. (The exoproteome study is Kaur 2018.)
- `longterm_coculture_rpom_dss3_proteomics` uses timepoint "long-term". The methods give 10 days for every seawater sample, so it should be "10 days" / 240 h.
- The experiment names "early coculture" and "late coculture" do not match the paper, where they are ASW vs natural seawater.
- `growth_phase: acclimated_steady_state` on the 35-day ASW samples is doubtful. The paper says axenic ASW cultures die after 4–6 weeks, so the axenic arm at 35 days may already be declining.

## Kaur 2018
- Sources read: main PDF (Methods: culture set-up, exoproteome preparation, Data analysis; Results on growth phases and PCA; Fig. 1 and 2 legends); legends.txt; CSV column names ("Fold change comparison of time point 1 vs different time points. Time point N vs 1").
- Paper default: `baseline`
- Sampling: cocultures of Synechococcus and R. pomeroyi in ASW and SW. "triplicate flasks for each one of the eight time points analyzed, i.e., days 1, 3, 7, 14, 21, 32, 60 and 100." Each later day is compared with day 1. There is no axenic arm.
- DE fitting: Perseus "two-sample Student's t test (P = 0.05) using a permutation-based false discovery rate", one test per pair. All contrasts sit within one experiment per organism × medium, so no shared id.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `syn_wh7803_exoproteome_sw` | day 3, 7, 14, 21, 32, 60, 100 | baseline (default) + design_notes | — | high | CSV header "Time point N vs 1"; p. 796 Methods: "days 1, 3, 7, 14, 21, 32, 60 and 100" |
| `syn_wh7803_exoproteome_asw` | same | baseline (default) + design_notes | — | high | same |
| `rpom_dss3_exoproteome_sw` | same | baseline (default) + design_notes | — | high | same |
| `rpom_dss3_exoproteome_asw` | same | baseline (default) + design_notes | — | high | same |

This matches the known case (4 → baseline vs day 1, and day 1 itself is not reported).

Opens / uncertainties:
- The day 1 control is labelled `exponential` in design_notes, inferred from the growth curves (Fig. 1) and from the text describing an exponential-to-stable-state shift at day 21. The paper does not label day 1 explicitly (medium confidence).

Other findings (not fixed):
- `growth_phase: acclimated_steady_state` is set on every timepoint. The paper says SW cocultures reach a stable state at day 21, and ASW cocultures only enter it after about 100 days (day 21 in ASW is a transient transition). Days 3–14 are exponential.
- The column header reads "time point 1 vs different time points". Check that the sign convention is later/day1 (the paperconfig assumes positive = up vs day 1).

## Huang 2020 (Huang 2021, MicrobiologyOpen)
- Sources read: main PDF (Methods: infection and sampling, transcriptome analysis / edgeR paragraph; Results on host DE); data sources.txt (no SI; the tables are in the text); GEO DE column names (`logFC-<tp>_vs_CT30m`).
- Paper default: `baseline`
- Sampling: "At 15 min, 1 h, 3 h, 5 h, and 7 h after phage lysate addition, 100 ml of culture was sampled … A parallel control group without phage addition was also sampled at the 30 min and 8 h time points." Because the two controls did not differ, "phage-treated samples were only compared to the control sample collected at 30 min."
- DE fitting: edgeR. All contrasts are within the single experiment, so no id.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `viral_s_sbp1_infection_wh7803_rnaseq` | 15 min, 1h, 3h, 5h, 7h | baseline (default) + design_notes | — | high | Methods: "phage-treated samples were only compared to the control sample collected at 30 min" |

Opens / uncertainties:
- The control is a parallel uninfected culture, not the same culture at t0. It is still sampled at one fixed time that matches none of the infected timepoints, so I chose `baseline`. The 30 min control is close to the 15 min and 1h points. The authors showed that 30 min and 8 h controls do not differ, which partly reassures, but the rule says `baseline`.

Other findings (not fixed): none beyond what the paperconfig comments already flag (unadjusted p carried in `adjusted_p_value_col`).

## Doron 2016
- Sources read: main PDF (Methods: Cultures and infection, RNA extraction, Microarray experimentation; Results with the Fig. 7 legend); supplementary PDF moesm36 (Microarray data analyses, Supplementary Fig. S1/S2 legends); geo_source.txt; GEO column names (`S30.C30`, `SYN9.1.C1`, …).
- Paper default: `time_matched`
- Sampling: SI: "three biological replicates for each time point after phage addition (2-5 min and 0.5, 1, 1.5, 2, 3, 4 h) of infected and uninfected control cultures." Each log-ratio column pairs infected and control at the same time (S30.C30).
- DE fitting: limma, one linear model per host strain comparing infected with uninfected over time. Each strain is its own experiment and the strains are different organisms, so no shared id is possible. All contrasts for one strain sit in one experiment.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `syn9_infection_wh7803` | 2-5 min, 30–240 min | time_matched (default) + design_notes | — | high | SI p. 6: "three biological replicates for each time point … of infected and uninfected control cultures"; columns S{t}.C{t} |
| `syn9_infection_wh8102` | 2-5 min, 0.5–4 h | time_matched (default) + design_notes | — | high | same |
| `syn9_infection_wh8109` | 2-5 min, 0.5–4 h | time_matched (default) + design_notes | — | high | same |

Opens / uncertainties:
- 2-5 min point: the control column is `C0`, and the Fig. 7 / Fig. S2 legends mention "the uninfected control cultures prior to infection relative to the infected cultures at 5 min". The earliest control may have been taken a few minutes before infection rather than at 2-5 min. The difference is minutes, and every other timepoint is clearly matched, so I kept `time_matched` and recorded the caveat in design_notes. If the earliest control was in fact a t0 sample, `mixed` would be strictly correct.

Other findings (not fixed): none.
