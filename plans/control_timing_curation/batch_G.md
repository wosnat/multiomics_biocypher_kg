# Batch G — control-timing curation

## Summary
- Papers: 6 · DE experiments: 12 · baseline: 0 · mixed: 0 · time_matched: 12 · shared-model ids: 1 (`mit0604_rnaseq_deseq2`, Hackl 2023, 3 experiments)
- Needs human review:
  - Hackl 2023 (all 3 experiments): methods do not say when the untreated controls were harvested relative to the 2h mitomycin C and 1h post-UV arms. `time_matched` is assumed from the parallel nine-culture design (low confidence).
  - Johnson 2026b (both experiments): no `statistical_model_id`. The methods do not say whether DESeq2 was fitted once on all 44 samples or per contrast, and the counts cannot tell the two apart. One Methods sentence says "ZT 0 ... as baseline", but the published fold changes are the same-ZT contrast (checked against the GEO counts).
  - Fuszard 2012: the replete and deplete arms were harvested at the same growth state (stationary onset, or mid-exponential for NATL2A), not on a stated common date. `time_matched` is by growth state (medium confidence).

## Dominguez-Martin 2017 (`Domínguez 2017/`)
- Sources read: main PDF Methods ("Cell collection", LC-MS/MS, Progenesis), Results section "Changes in the proteome ... provoked by azaserine addition", Fig. 5 legend; legends.txt; paperconfig.
- Paper default: `time_matched`
- Sampling: one 10 L exponential culture was split in two. One half got 100 uM azaserine and the other was the control, and both were kept under the same conditions. Proteomics samples were taken from 3 biological replicates "after 8 h of growth in the absence of azaserine (control cultures) and in the presence of 100 uM azaserine".
- DE fitting: Progenesis QI label-free quantitation, ANOVA p and q < 0.05. There is a single comparison.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `azaserine_ss120_proteomics` | 24h | time_matched (default) | — | high | Methods: "Cultures ... were split into two aliquots: one was used as the control culture, while 100 uM azaserine ... was added to the other"; Results: "collected samples ... after 8 h of growth in the absence of azaserine (control cultures) and in the presence of 100 uM azaserine" |

Opens / uncertainties: none.

Other findings (not fixed):
- Timepoint: the paper samples proteomics after **8 h**. The paperconfig has `timepoint: "24h"` / `timepoint_hours: 24`, and `experimental_context` says "harvested 24 h after azaserine addition". The 24 h figure in the paper is the GS-activity time course, not the proteome.
- Concentration: the paper uses **100 uM** azaserine. `treatment_condition` and the header comment say 10 uM.
- Light: Methods give continuous blue light at **4 uE m-2 s-1**. The paperconfig `light_intensity` says 25 umol photons m-2 s-1 (this may come from another section; worth checking).

## Fuszard 2012
- Sources read: main PDF (Results "Overview", Methods summary on p. 6–7); Additional file 3 Materials and Methods (MOESM3 .DOC); Table S1 CSV header (`Control MEAN, deplete MEAN, deplete SD, Ratio of means`); paperconfig.
- Paper default: `time_matched`
- Sampling: for each strain, triplicate replete (50 uM) and deplete (10 uM Pi) cultures were grown in parallel. Supplementary methods: "MIT 9312 and SS120 cultures were harvested upon the initiation of stationary phase, and NATL2A cultures were harvested in mid-exponential phase". Main text: harvested "once measured optical densities reached 0.4". Both arms were harvested at the same growth state. No common calendar date is stated, and MIT9312 replete grew faster than deplete.
- DE fitting: iTRAQ, one multiplex per strain with replete and deplete channels. Ratio of means with fold-change cutoffs only; no per-protein test. Each strain is a separate run, so there is no shared model.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `pi_limitation_mit9312_itraq` | stationary | time_matched (default) | — | medium | Supp methods: "harvested upon the initiation of stationary phase"; same growth state, not a stated common time |
| `pi_limitation_natl2a_itraq` | mid-exponential | time_matched (default) | — | medium | Supp methods: "NATL2A cultures were harvested in mid-exponential phase" |
| `pi_limitation_ss120_itraq` | stationary | time_matched (default) | — | medium | as MIT9312 |

Opens / uncertainties:
- All 3: `time_matched` here means acclimated arms harvested at the same growth state (OD 0.4 or stationary onset), not the same clock time. I did not add `design_notes`. Add one if the reviewer wants that nuance recorded.

Other findings (not fixed):
- **Axenic status:** supplementary methods open with "Non-axenic Prochlorococcus strains MIT 9312, NATL2A and SS120". All three experiments carry `background_factors: ["light", "axenic"]`.
- **Light intensities swapped:** both the main text ("30, 10 and 20 uE ... respectively" for MIT9312, NATL2A, SS120) and the supplementary methods (MIT9312 30, SS120 20, NATL2A 10) give NATL2A 10 and SS120 20. The paperconfig has NATL2A 20 and SS120 10.
- The main text says 8-plex iTRAQ reagents; the paperconfig and README say 4-plex. This is minor.

## Pandhal 2007
- Sources read: main PDF Methods (culture acclimation, protein extraction, iTRAQ, Table 1 label layout, ProQuant analysis) and Results on differential regulation; paperconfig.
- Paper default: `time_matched`
- Sampling: MED4 was acclimated to each light level (LL 20, ML 60, HL 100 uE) through 3 sub-culturings. For each light level, "50 mL of cells was harvested at mid-exponential phase". Growth was measured at the same time of day. The three light levels were run as parallel acclimated cultures harvested at the same growth state.
- DE fitting: ProQuant / ProGroup per-pair iTRAQ ratios with a weighted-error-factor p-value. HL/LL and ML/LL are ratios against the same LL reporter channel in the same two iTRAQ runs. They are not one fitted statistical model, so no `statistical_model_id` is set. The KG aggregation (mean log2 ratio, max p) is our own processing.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `light_HL_vs_LL_med4_proteomics` | acclimated mid-exponential | time_matched (default) | — | medium | Methods: "subcultured three times in late-exponential phase ... to ensure adaptation"; "50 mL of cells was harvested at mid-exponential phase" |
| `light_ML_vs_LL_med4_proteomics` | acclimated mid-exponential | time_matched (default) | — | medium | same |

Opens / uncertainties:
- Both experiments share the LL reference channel inside the same iTRAQ multiplexes, so their errors are correlated. Under the plan's definition (one fitted model) I did not set an id. Set one if the reviewer counts a shared iTRAQ reference channel as a shared model.

Other findings (not fixed): none.

## Hackl 2023
- Sources read: main PDF STAR Methods "RNA-sequencing experiments" and "RNA-sequencing analyses"; geo_source.txt; legend.txt; headers and first rows of the three GEO sense DESeq2 tables (`baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, gene`).
- Paper default: `time_matched`
- Sampling: "Nine 30 mL cultures ... grown to exponential phase. 3 cultures were treated with mitomycin C ... for 2 hours; 3 cultures were applied a UV shock ... then placed back ... for 1h; the 3 remaining cultures were kept as control with no treatment. Cells were harvested by centrifugation". The harvest time of the controls is not stated.
- DE fitting: DESeq2 v1.24.0 Wald test with BH correction, padj < 0.1. **One model:** `baseMean` is identical across the three GEO tables (e.g. MIT0604_00001 = 2433.1033 in all three). The contrasts are also exactly linear: MMC vs UV = MMC vs control − UV vs control (−0.6846 − (−0.3496) = −0.3350, as reported).

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `mmc_dna_damage_mit0604` | 2h | time_matched (default) | mit0604_rnaseq_deseq2 | low (timing) / high (model) | STAR Methods "the 3 remaining cultures were kept as control with no treatment"; identical baseMean across the GEO tables |
| `uv_dna_damage_mit0604` | 1h post-shock | time_matched (default) | mit0604_rnaseq_deseq2 | low / high | same |
| `mmc_vs_uv_mit0604` | endpoint | time_matched (default) | mit0604_rnaseq_deseq2 | low / high | treatment-vs-treatment contrast from the same fit; exact linear relation verified |

`design_notes` was added to all three. Each note says the control harvest time is unstated and `time_matched` is assumed from the parallel design, and records the shared DESeq2 fit and the exact linear relation.

Opens / uncertainties:
- All 3: the treated arms have different exposure lengths (2h mitomycin C vs 30 s UV plus 1h recovery), so the control cannot be clock-matched to both arms unless the treatments were staggered. Asking the authors, or the GEO sample metadata, would settle it.

Other findings (not fixed):
- The Methods give mitomycin C as "15 mg mL-1" for RNA-seq and "20 mM" in the qPCR survey; the paperconfig uses 15 ug/mL. The paper's units look like a typo. No action needed.

## Johnson 2026b
- Sources read: main PDF Results "light-decoupling" section, Fig. 4 legend, Methods 5.x sampling and 5.8 statistics; media-1 `metadata` sheet (44 Circadian samples); media-2 DESeq2 sheet headers ("2H Dark vs Light", "14H Light vs Dark"); GEO GSE314951 count matrix.
- Paper default: `time_matched`
- Sampling: a diel-entrained MED4 culture was sampled every 2 h over 24 h (ZT 0–22, triplicates). Extra perturbation samples (duplicates) were taken at ZT 2 and 4 under extended darkness and at ZT 14 and 16 under extended light. The Methods say the subset used "the ZT 0 circadian time-course sample as baseline". However, recomputing log2 ratios from the GEO counts (median-of-ratios normalisation) reproduces the published log2FC only against the **same-ZT** normal-diel sample: r = 1.000 and median absolute difference ≤ 0.035 for all four timepoints. The alternatives fit much worse: vs ZT 0, r = 0.41 and 0.54; vs ZT 12, r = 0.81 and 0.57.
- DE fitting: DESeq2 v1.46, padj < 0.05. The design or formula is not stated. It is not stated whether one model was fitted on all samples. Pairwise and all-sample size factors both reproduce the values to about 0.003, so the data cannot settle it.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `extended_darkness_med4` | 2h, 4h | time_matched (default) | — | high (verified from counts) | media-2 header "log2FC 2H Dark vs Light"; GEO counts reproduce 2Hr-dk vs 2Hr (r = 1.000) |
| `extended_light_med4` | 14h, 16h | time_matched (default) | — | high (verified from counts) | "log2FC 14H Light vs Dark"; 14Hr-lt vs 14Hr (r = 1.000) |

`design_notes` was added to both. They name the same-ZT circadian control and record that the Methods "ZT 0 baseline" sentence does not match the published numbers.

Opens / uncertainties:
- No `statistical_model_id`, because the fitting is unknown (see Summary). If the reviewer can confirm one DESeq2 fit over all 44 samples, both experiments should get e.g. `med4_circadian_rnaseq_deseq2`.

Other findings (not fixed):
- The control arm is the triplicate circadian time-course culture and the perturbation arm is a separately collected duplicate subset. They are time-matched in ZT, but it is not clear they are the same flask series.

## Tai 2009 (`Synechococcus/.../Tal 2009/`)
- Sources read: main PDF Fig. 1 legend, Results growth paragraph, Experimental procedures "RNA extraction" and microarray analysis; SI PDF skimmed for sampling (nothing more on timing); paperconfig.
- Paper default: `time_matched`
- Sampling: 1 L batch mono- and co-cultures were inoculated at 1e6 WH8102 cells/mL, "grown to mid-log phase (approximately 1e8 cells/mL ...) and harvested". Growth rates did not differ (0.417 vs 0.416 per day). Both arms were harvested at the same growth state.
- DE fitting: SAM on log2 signal ratios of replicate spots. There is a single comparison.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `vibrio_coculture_wh8102_microarray` | mid-log phase | time_matched (default) | — | medium | Experimental procedures: "From 1 l batch mono- or co-cultures ... grown to mid-log phase (approximately 1e8 cells ml-1, see Fig. 1) ... and harvested"; equal growth rates |

Opens / uncertainties: none beyond medium confidence (same growth state, same inoculum and growth rate, so effectively the same day).

Other findings (not fixed):
- The directory is named "Tal 2009" but the author is Vera **Tai** (the paperconfig `papername` "Tai 2009" is correct). There is no `doi` field in the paperconfig.
- The validator warns that `Gene Name` is typed `locus_tag` but looks like free text. This warning was already there.

## Process note
- While checking validation I briefly ran `git stash -- <path>` on my six paperconfigs. It stashed only those six files. I popped all six stashes straight away, `git stash list` is empty, and the diff is intact: 39 inserted lines over the 6 files. No other paper was touched.
- `validate_paperconfig.py` reports VALIDATION PASSED on all 6 edited files. Every warning it shows existed before (Hackl island `experiments` fields, Tai `Gene Name`).
