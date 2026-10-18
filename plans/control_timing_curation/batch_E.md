# Batch E — control-timing curation

## Summary
- Papers: 5 · DE experiments: 15 · baseline: 2 · mixed: 2 · time_matched: 11 · shared-model ids: 0
- Needs human review:
  - Lin 2015 `phosphorus_plimited_natl2a_rnaseq_uninfected`: set `mixed` because Table S1 has no nutrient-replete library at 50 h, so the "50h (P added)" contrast cannot be time-matched. The paper does not say which replete sample (46 h or 59 h) served as its control.
  - Tetu 2019 `plastic_pvc_plastic_leachate_2_mit9312_rnaseq` / `plastic_pvc_plastic_leachate_2_natl2a_rnaseq`: set `baseline` (override). The evidence is only the GEO GSE118155 protocol: PVC was harvested at 90 min, but the control and HDPE cultures at 180 min, and one control set serves both. The paper does not say when the control was harvested.
  - Tetu 2019 HDPE experiments: the paper says HDPE was exposed for 120 min, GEO says 180 min. The paperconfig timepoint is `120 min`.
  - Tetu 2019 and Lin 2015: DE fitting is not described well enough to decide whether a shared model was used, so no `statistical_model_id` was set (see Opens).

## Lin 2015
- Sources read: main PDF (Results p. 451–453 incl. Fig. 2 legend; Experimental procedures p. 456–457, incl. "Identification of differentially expressed genes"); Table S1.docx (library list); `legend of supp tables.txt`; headers of `table S4A.csv` / `table S4B.csv`; README.md.
- Paper default: `time_matched` (1 experiment `time_matched`, 1 `mixed`; tie broken toward the infected experiment, which is fully matched)
- Sampling:
  - Log-phase NATL2A was washed into P-replete or P-depleted Pro99. Uninfected cultures were sampled at 4, 24, 46 and 59 h (both arms).
  - Aliquots of the P-limited cultures had P re-added at 46 h and were sampled at 50 and 59 h.
  - At 47 h, aliquots of both arms were infected and sampled at 47, 48, 49, 51, 55 and 59 h. P was re-added to infected P-limited aliquots at 51 h and those were sampled at 55 and 59 h.
  - Table S1 lists the replete libraries at exactly those hours. There is no replete-uninfected library at 50 h.
- DE fitting: DESeq (v1) "with default parameters"; "Two biological replicates were conducted for each time point". There is no statement of one model across time points. Each column pair in S4 reads as a separate -P/control or P-added/control contrast.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `phosphorus_plimited_natl2a_rnaseq_uninfected` | 4h, 24h, 46h, 59h, 50h (P added), 59h (P added) | mixed (override) | — | medium | Table S1: uninfected nutrient-replete libraries only at 4/24/46/59 h; S4 legend "Log2 fold change (P-limited/nutrient-replete)"; Fig. 2A legend "P-limited (-P) normalized to nutrient-replete (control) cultures" |
| `phosphorus_plimited_natl2a_rnaseq_infected` | 47h, 48h, 49h, 51h, 55h, 59h, 55h (P added), 59h (P added) | time_matched (default) | — | high | Table S1: infected replete libraries at 47/48/49/51/55/59 h, matching every infected -P and P-reintroduced sample; Fig. 2B legend "aliquots from both P-limited and nutrient-replete cultures were infected" |

Opens / uncertainties:
- `phosphorus_plimited_natl2a_rnaseq_uninfected`: the control for the 50h (P added) column is not stated. It must be a replete library from another hour (46 h or 59 h). I chose `mixed` and recorded this in `design_notes`. The control growth phase is written as exponential, which is inferred. The supplementary figures doc or the authors could resolve it.
- No `statistical_model_id`: the paper does not describe the DESeq design. Its "default parameters" and per-time-point replicates suggest pairwise fits.

Other findings (not fixed):
- The paperconfig has no `doi`. The README gives 10.1111/1462-2920.13104.
- Culture conditions in the paperconfig (24C, 14:10 light:dark, 40 umol, `background_factors` including `diel`) do not match the Methods, which say "maintained in the Pro99 medium at 21°C under constant cool white light (30 μmol quanta m–2 s–1)" (p. 456). The RNA was reused from Zeng & Chisholm 2012, so check that paper's conditions. `diel` looks wrong.
- Both experiments have the same `name` ("NATL2A P-limited vs P-replete control (RNASEQ)"). The infected one should say infected.
- The uninfected experiment's `background_factors` has `axenic`, but the infected one does not. Both cultures are axenic hosts.

## Lindell 2007
- Sources read: main PDF (Methods, Fig. 3 legend); Supplementary information PDF (Methods "Phage infection", "Array Data Analyses"; Suppl. Table 3 legend); `supp table legend.txt`; `supp table 3.csv` header; README.md.
- Paper default: `mixed` (single DE experiment)
- Sampling: triplicate MED4 cultures (10^8 cells/ml) were infected with P-SSP7 (MOI 3). The paired control cultures were amended with filter-sterilized spent medium. Samples were collected each hour for 0–8 h.
- DE fitting: Bayesian t-test (GoldenSpike / hdarray) between infected and control at each time point, with q-values. The tests are per time point, and there is a single experiment, so there is no shared-model id.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `viral_phage_infected_med4_microarray` | 1h, 2h, 3h, 4h, 5h, 6h, 7h, 8h | mixed (default) | — | high (for 5–7 h); medium (for 1–3 h being matched) | Supp. Methods: "Control arrays at 4 and 8 h after infection gave the same expression profiles (data not shown) and were used for comparison with infected cells at 5, 6, 7 h after infection." Main Methods: "paired control cultures were amended with filter-sterilized spent medium" |

Opens / uncertainties:
- 1h, 2h and 3h are recorded as time-matched because the paired controls were sampled alongside the infected cultures, and the 5–7 h exception implies controls existed at other hours. The paper does not list the control arrays explicitly. Even if 1–3 h were not matched, the value stays `mixed`, because 4h and 8h are matched.

Other findings (not fixed):
- `control_condition` is "Pro99 medium phage infection baseline". The actual control is a paired uninfected culture amended with spent medium, so "baseline" is misleading.
- `temperature: 24C` and `light_intensity: 10 umol` disagree with the Supp. Methods: "21 °C under continuous white light at 10-25 µmol photon".
- `treatment_organism: Phage` (taxid 10239) is generic; the phage is P-SSP7.
- `growth_phase: infected` on every analysis. Check that this is in the validator vocabulary.

## Martiny 2006
- Sources read: main PDF (Methods "P-Starvation Time Series", "Array Analysis"; Fig. 1 legend); headers of `01301table_1.csv` / `01301table_2.csv`; README.md.
- Paper default: `time_matched`
- Sampling: triplicate 4 L cultures were split in two, then washed and resuspended in P-replete or P-depleted Pro99. "Samples were taken … at 0, 4, 12, 24, and 48 h after resuspension" from both arms. Tables give FC (-P/+P) per time point.
- DE fitting: Cyber-T Bayesian t-test (separate window sizes per strain) plus QVALUE, one test per time point. The two experiments are different strains, so there is no shared-model id.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `phosphorus_phosphate_starvation_med4_microarray` | 0h, 4h, 12h, 24h, 48h | time_matched (default) | — | high | Methods: cultures "split in two, and washed twice in either P-replete … or -depleted … Samples were taken … at 0, 4, 12, 24, and 48 h"; Table 1 "FC, fold change (-P/+P)" per t |
| `phosphorus_phosphate_starvation_mit9313_microarray` | 0h, 4h, 12h, 24h | time_matched (default) | — | high | same; Table 2 "fold-change (-P/+P)" |

Opens / uncertainties: none.

Other findings (not fixed):
- MIT9313 `table_scope_detail` says "Genes with q < 0.05 at t=48h", but the Table 2 caption says "differentially expressed (q<0.05 at t=24)". MIT9313 has no 48 h column.
- `light_intensity` on both experiments lists both strains' values ("30 uE … (MED4), 12 uE … (MIT9313)"). Each experiment should carry only its own strain's value.

## Steglich 2006
- Sources read: main PDF (Materials and Methods: culture/experimental conditions, "Normalization and LPE test"; Fig. 3 legend); `supp table legend.txt`; `supp table s1.csv` header; README.md.
- Paper default: `time_matched`
- Sampling: triplicate cultures were concentrated, split into seven subcultures, and dark-adapted for 5 h. They were then shifted for 45 min to HL, WL, WL+DCMU, BL, RL or GL; "Control cultures were kept in darkness". All were harvested by rapid filtration.
- DE fitting: one array per sample. Fold change = condition / reference sample (dark; WL for DCMU). Significance comes from a separate LPE z-test plus BH for each comparison. The five light-vs-dark experiments share the dark reference arrays, but each was tested separately, so there is no shared-model id.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `light_stress_high_white_light_55_med4_microarray` | 45 min | time_matched (default) | — | medium | Methods: "divided into seven 20-ml subcultures … shifted to the respective experimental conditions for 45 min … Control cultures were kept in darkness" |
| `light_stress_white_light_13_umol_med4_microarray` | 45 min | time_matched (default) | — | medium | same |
| `light_stress_blue_light_13_umol_med4_microarray` | 45 min | time_matched (default) | — | medium | same |
| `light_stress_red_light_13_umol_med4_microarray` | 45 min | time_matched (default) | — | medium | same |
| `light_stress_green_light_13_umol_med4_microarray` | 45 min | time_matched (default) | — | medium | same |
| `light_stress_dcmu_white_light_13_med4_microarray` | 45 min | time_matched (default) | — | high | DCMU and WL subcultures both received the 45 min WL shift in parallel |

Opens / uncertainties:
- The dark control harvest time is not stated explicitly. The parallel subcultures imply it was harvested alongside the 45 min samples, hence `medium`.

Other findings (not fixed): none.

## Tetu 2019
- Sources read: main PDF (Methods "RNA isolation, sequencing and analysis", "Statistics and reproducibility"; Results on transcriptomics); Supplementary PDFs MOESM1/MOESM2 (dataset legends); headers of both DE CSVs and both GEO count files; `geo_source.txt`; GEO GSE118155 sample records (treatment protocol, fetched).
- Paper default: `time_matched` (2 `time_matched`, 2 `baseline`; tie broken toward HDPE)
- Sampling:
  - The paper says RNA was extracted from triplicate 50 mL control and leachate cultures "after short-term exposure (90 min for PVC, 120 min for HDPE)".
  - GEO GSE118155 (identical text on all 18 samples) says: "PVC spiked cultures were harvested after 90 mins exposure, control and HDPE exposed cultures were harvested after 180 mins".
  - The count files contain a single control set per strain (C4–C6) next to H4–H6 and P1–P3.
- DE fitting: DESeq2 "to identify differential expression between control and leachate exposed cultures" (padj < 0.01). The paper does not say whether one design (~condition over all 9 samples) or two pairwise fits were used. HDPE and PVC share the control libraries.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `plastic_hdpe_plastic_leachate_50_mit9312_rnaseq` | 120 min | time_matched (default) | — | medium | GEO: "control and HDPE exposed cultures were harvested after 180 mins" |
| `plastic_pvc_plastic_leachate_2_mit9312_rnaseq` | 90 min | baseline (override) | — | medium | GEO: "PVC spiked cultures were harvested after 90 mins … control … after 180 mins"; one control set C4–C6 in GSE118155 counts |
| `plastic_hdpe_plastic_leachate_50_natl2a_rnaseq` | 120 min | time_matched (default) | — | medium | as above |
| `plastic_pvc_plastic_leachate_2_natl2a_rnaseq` | 90 min | baseline (override) | — | medium | as above |

Opens / uncertainties:
- PVC `baseline`: this rests on the GEO protocol only; the paper is silent on when the control was harvested. If the authors in fact harvested controls at both 90 and 180 min, the GEO deposit would show two control sets, and it shows one. That makes `baseline` the value the deposit implies. The direction is also unusual: the control was sampled later (180 min) than the treatment (90 min).
- No `statistical_model_id`: the HDPE and PVC experiments of each strain are plausibly two contrasts of one DESeq2 `~condition` fit over C, H and P (one count file per strain, shared control). The paper does not say so, so no id was set. This could be checked with the authors, or by re-running DESeq2 on GSE118155 counts and comparing against Supplementary Datasets 3/4.

Other findings (not fixed):
- HDPE exposure time: the paper (and the paperconfig `timepoint: 120 min`, `timepoint_hours: 2`, and `experimental_context` "120 min") disagree with GEO, which says 180 min. One of the two sources is wrong.

## Validation
- All five paperconfigs parse. `scripts/validate_paperconfig.py` reports VALIDATION PASSED on each after the edit.
