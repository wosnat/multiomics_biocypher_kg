# Batch C — control-timing curation

## Summary
- Papers: 5 · DE experiments: 22 · baseline: 2 · mixed: 0 · time_matched: 20 · shared-model ids: 3 (`ez55_rnaseq_edger_glm` ×8 in Barreto 2022, `mit1002_rnaseq_deseq2` ×2 in Biller 2016, `natl2a_rnaseq_deseq2` ×2 in Biller 2018)
- Needs human review:
  - Barreto 2022, the 6 `coculture_*_at_{400,800}_ez55_rnaseq` experiments: the tables are **not** partner-vs-axenic. They are partner deviations from the partner-averaged coculture response (verified numerically, see Other findings). The experiment names and `control_condition` say otherwise. Recorded in `design_notes` but not fixed.
  - Barreto 2022 `carbon_stress_elevated_pco2_800_ppm_ez55_rnaseq` and `interaction_pco2_coculture_ez55_rnaseq`: the author script defines these coefficients differently from what the paperconfig says (axenic-only pCO2 effect; non-standard "interaction"). Recorded in `design_notes`, not fixed.
  - Barreto 2022 and Hennon 2017 share sequencing libraries: Barreto re-analyses Hennon's MIT9312 cocultures and axenic EZ55. Both papers' DE experiments therefore restate the same samples. `statistical_model_id` is paper-scoped, so it cannot express this.
  - Biller 2016 `mit1002_rnaseq_deseq2` and Biller 2018 `natl2a_rnaseq_deseq2` (medium confidence): the methods describe one DESeq2 dataset per organism with pairwise Wald tests. That is a single fit, but the methods never show the code. In Biller 2016 the two MIT1002 experiments also share the same 12 h reference samples.
  - Hennon 2017: one EZ55 GLM or two fits. The evidence says two separate fits, so no id was set. See "Decisions for the curator".

## Barreto 2022
- Sources read: main PDF pp. 1–4 (Methods: Strains, Culture conditions, Growth experiments, RNA library prep, Gene transcription analysis) and Results pp. 3–5 (Fig. 2 legend). Supp PDF table legends S1–S14 and Fig. S1–S9 legends. Author R scripts `Alteromonas_analysis/dge_analysis_Alteromonas_.txt`, `deg_analysis_MIT9312.txt`, `dge_analysis_CC9311.txt` and `dge_analysis_WH8102.txt`, plus `targets.xlsx`. CSV headers of S1_*, S2, S3, S4, S5 and S9–S14. README.md and mail.md.
- Paper default: `time_matched`
- Sampling: every culture was acclimated to 400 or 800 ppm for three transfer cycles (about 14 generations) and harvested in the final transfer cycle. Synechococcus cultures were diluted whenever they reached 2.6e5 cells/mL, so they stayed in exponential growth (p. 2). The axenic EZ55 and MIT9312 coculture libraries come from the earlier Hennon 2017 experiment: "all reads from our previous experiment with MIT9312 and EZ55 [7] were analyzed together here" (p. 2, ref 7 = Hennon et al. 2018). So all arms were harvested at the same acclimated growth state, but not in parallel on one date.
- DE fitting:
  - Cyanobacteria: one edgeR fit per strain, `glmFit(y, ~Treatment)` then `glmLRT(coef=2)`, i.e. separate fits.
  - EZ55: one edgeR GLM over all 48 EZ55 libraries, `design_1_interactcynoc <- model.matrix(~ 0+ cyano + treatment+cyano:treatment)`, with 9 contrasts via `glmLRT` / `makeContrasts` (script lines 68–599). Methods p. 3: "Statistical analysis was performed based on our fitted model and the following custom contrasts were applied using the function makeContrasts".

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `carbon_stress_elevated_pco2_800_ppm_mit9312_rnaseq` | mid-exponential | time_matched (default) | — | medium | p. 2 "each culture was acclimated to experimental conditions for three transfer cycles"; MIT9312 script: separate `glmFit(y,design_2)` |
| `carbon_stress_elevated_pco2_800_ppm_cc9311_rnaseq` | mid-exponential | time_matched (default) | — | medium | same; CC9311 script: separate fit `~Treatment` |
| `carbon_stress_elevated_pco2_800_ppm_wh8102_rnaseq` | mid-exponential | time_matched (default) | — | medium | same; WH8102 script: separate fit |
| `carbon_stress_elevated_pco2_800_ppm_ez55_rnaseq` | mid-exponential | time_matched (default) | ez55_rnaseq_edger_glm | high (model) / medium (timing) | script l.102 `contrast1<-glmLRT(fit_1_interactcynoc, coef=5)` |
| `coculture_prochlorococcus_mit9312_at_400_ez55_rnaseq` | mid-exponential | time_matched (default) | ez55_rnaseq_edger_glm | high / medium | script contrasts 3–5; S5 arithmetic (below) |
| `coculture_synechococcus_cc9311_at_400_ez55_rnaseq` | mid-exponential | time_matched (default) | ez55_rnaseq_edger_glm | high / medium | same |
| `coculture_synechococcus_wh8102_at_400_ez55_rnaseq` | mid-exponential | time_matched (default) | ez55_rnaseq_edger_glm | high / medium | same |
| `coculture_prochlorococcus_mit9312_at_800_ez55_rnaseq` | mid-exponential | time_matched (default) | ez55_rnaseq_edger_glm | high / medium | script contrasts 7–9; S5 arithmetic |
| `coculture_synechococcus_cc9311_at_800_ez55_rnaseq` | mid-exponential | time_matched (default) | ez55_rnaseq_edger_glm | high / medium | same |
| `coculture_synechococcus_wh8102_at_800_ez55_rnaseq` | mid-exponential | time_matched (default) | ez55_rnaseq_edger_glm | high / medium | same |
| `interaction_pco2_coculture_ez55_rnaseq` | mid-exponential | time_matched (default) | ez55_rnaseq_edger_glm | high / medium | script l.413 contrast6 |

`design_notes` was added on all 8 EZ55 experiments. Each note names the shared model, says what the coefficient or contrast actually is, and records that the axenic and MIT9312 arms are Hennon 2017 libraries.

Opens / uncertainties:
- Timing for all 11 experiments is `time_matched` under the rule "acclimated cultures harvested at the same growth state". The arms were not harvested in parallel, though: the axenic EZ55 and MIT9312 cocultures were grown and sequenced for Hennon 2017, and the Synechococcus cocultures were grown later in a different medium (SEv rather than PEv). This confound between partner and batch is written into `design_notes`. A curator may prefer to treat it as weaker than "time_matched".
- The harvest density in the final transfer cycle is not stated for the Synechococcus cultures (medium confidence on growth state).
- S9–S14 do not appear in the author script, which only goes up to contrast 9. Their identity as contrasts 3–5 and 7–9 is inferred from the S5 arithmetic below. That identity is exact at both pCO2 levels, so I am confident they come from the shared fit.

Other findings (not fixed):
- **S9–S14 are mislabelled as "coculture with X vs axenic".**
  - Table S5 carries both `logFC` and `logFC_rel_axenic`. Its `logFC` equals the S9–S14 logFC for the same gene, and `logFC_rel_axenic − logFC` equals the `General_coculture` logFC at that pCO2 (max deviation 1e-9 over all 25 genes × 3 partners × 2 levels).
  - The three partner logFCs also sum to zero per gene at both 400 and 800 ppm.
  - So each table is the partner's deviation from the general coculture response. That matches the methods' item (iii), "the effect of specific cyanobacteria in addition to the general co-culture response", and the Fig. 2 legend ("black bars in columns 2-4 indicate significant difference between the specific cyanobacterial response and the general coculture response").
  - The experiment `name`, `treatment_condition`, `control_condition` and `experimental_context` of all 6 should be revised. The real partner-vs-axenic effect needs S3 plus S9–S11, which is not available per gene for non-significant genes.
- **`carbon_stress_elevated_pco2_800_ppm_ez55_rnaseq` (S2) is the pCO2 effect in axenic EZ55 only.** It is coefficient 5 (`treatmentpositive`) of `~0+cyano+treatment+cyano:treatment`, and with the cyano:treatment interaction columns that coefficient applies to the reference cyano level (`cyanoNone`). The paperconfig `experimental_context` "across axenic and coculture conditions" and its `background_factors: [axenic, coculture]` / partner `background_organisms` are therefore wrong. The table effectively re-analyses Hennon 2017's axenic experiment.
- **S4 / `interaction_pco2_coculture_ez55_rnaseq` is ambiguous.** The script's contrast6 is `(Pro:trt + CC:trt + W:trt)/3 - treatmentpositive`, which is not a clean interaction term. Separately, Table S5's `General_coculture` 800 ppm logFC equals the S4 logFC exactly (6 overlapping genes), i.e. the authors also present S4 as "general coculture response at 800 ppm". The paper's own description is internally inconsistent; needs author clarification.
- **Culture conditions in the paperconfig contradict the methods.** The paperconfig says `medium: Pro99 medium`, `temperature: 22C`, `light_condition: continuous light`, `light_intensity: 30 umol photons m-2 s-1` on all 11 experiments. Methods p. 2: artificial seawater with SEv (Synechococcus) or PEv (MIT9312) nutrients, "21 °C under 150 μmol photons m−2 s−1 on a 14:10 light:dark cycle". The coculture experiments' `background_factors: [light]` arguably should be `diel`.
- **Cross-paper sample overlap:**
  - Barreto `carbon_stress_elevated_pco2_800_ppm_mit9312_rnaseq` re-analyses the same libraries as Hennon `carbon_stress_elevated_co2_800_ppm_mit9312_rnaseq`.
  - Barreto `carbon_stress_elevated_pco2_800_ppm_ez55_rnaseq` (axenic coefficient) re-analyses the same samples as Hennon `carbon_stress_elevated_co2_800_ppm_ez55_rnaseq_axenic`.
  - These are not independent replications. There is no field to express this.
- Significance is unadjusted p < 0.05 with |logFC| > 1 (methods p. 3: "unadjusted p values < 0.05"). The paperconfig maps `PValue` to `adjusted_p_value_col`, so the edges carry unadjusted p-values labelled as adjusted.

## Hennon 2017
- Sources read: main PDF pp. 520–523 (Materials and methods in full: Cultures, Growth experiments, RNA library prep, Differential gene expression) and Table 2 notes. Supp docx moesm27 (Supp Table 1: carbonate chemistry, RNA co-cultures t = 0 vs RNA axenic t = final) and moesm32. `legends of supp tables.txt`. CSV headers of S2/S3/S4 (columns logFC, logCPM, LR, PValue). README.md.
- Paper default: `time_matched`
- Sampling: cocultures and axenic EZ55 were transferred twice at ambient or elevated CO2 (about 10 generations).
  - Cocultures were split into five tubes and harvested "when at least three of these replicates crossed the cutoff value of about 5 × 10^6 cells ml−1" (p. 522).
  - Axenic EZ55 was grown separately in 1 L bottles and harvested at about 3 × 10^8 cells per 200 ml.
  - Ambient and elevated arms are therefore harvested at the same growth state (density-triggered), not on a common date. The time of day within the 14:10 L:D cycle is not stated.
- DE fitting: edgeR, "a general linear model to estimate the tag-wise common dispersion with a paired-sample design to minimize the impact of strain–strain variability" (p. 523). The tables carry `LR`, so this is a glmLRT. The text's "exact test" is inconsistent with that.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `carbon_stress_elevated_co2_800_ppm_mit9312_rnaseq` | acclimated | time_matched (default) | — | high | p. 522 RNA prep: acclimation then density-triggered harvest of both CO2 arms |
| `carbon_stress_elevated_co2_800_ppm_ez55_rnaseq_coculture` | acclimated | time_matched (default) | — | high | same |
| `carbon_stress_elevated_co2_800_ppm_ez55_rnaseq_axenic` | acclimated | time_matched (default) | — | high | p. 522 "axenic Alteromonas cultures were grown in 1 liter polycarbonate bottles", harvested at about 3e8 cells |

### Decisions for the curator — Hennon 2017: one EZ55 GLM or two fits?

**Recommendation: two separate fits; no `statistical_model_id`. Confidence: high.**

Evidence:
1. **Decisive: `logCPM` differs between the tables for the same gene.** In edgeR, the `logCPM` column of a `glmLRT` table is `aveLogCPM`, computed from the DGEList passed to `glmFit`. Every contrast of one fit therefore reports identical `logCPM` per gene. For the 48 genes present in both S3 (coculture) and S4 (axenic), `logCPM` differs by −2.5 to +3.6 (e.g. AEZ55_2834: 4.02 vs 6.57; AEZ55_4628: 9.30 vs 10.25). The two tables were therefore produced from two different DGEList objects, i.e. two fits.
2. The methods describe a single CO2 contrast "in elevated CO2 compared with ambient" with a paired (clone) design. They do not mention a culture-type factor or a CO2 × culture interaction.
3. Coculture and axenic cultures were grown in different vessels and harvested by different criteria. Their carbonate chemistry was measured at different times (Supp Table 1: co-cultures t = 0, axenic t = final). The libraries were also deposited separately (coculture SRR5320558–69; axenic BioSamples SAMN07374746–57). Together these show two separate sub-experiments.
4. The MIT9312 table comes from the coculture libraries but is a different organism, so it cannot be linked to the EZ55 experiments.

What would overturn this: the authors' code, or an exactly matching `logCPM`, which it does not.

Opens / uncertainties:
- None on control_timing. Minor: the harvest time of day under 14:10 L:D is not stated for either arm.

Other findings (not fixed):
- The methods say significance is from an "exact test", but the tables carry an `LR` column, which comes from a GLM likelihood-ratio test. The paperconfig `test_type: edgeR` is unaffected.
- `PValue` is mapped to `adjusted_p_value_col`. The paper reports unadjusted p < 0.05.
- Shares libraries with Barreto 2022; see the Barreto Other findings.

## Biller 2016
- Sources read: main PDF pp. 2831–2840 (Materials and methods: Culturing and sampling, RNAseq library construction, DE analysis), Fig. 1 legend, Results "Timing and general patterns" and the Alteromonas section with Table 3. CSV headers of S2 (NATL2A, 7 timepoint columns) and S3 (MIT1002, "24 vs 12 hrs", "48 vs 12 hrs"). `geo_source.txt` and the GSE73511 count header (T{n}Ax / T{n}Alt × 3). README.md.
- Paper default: `baseline` (2 of 3 DE experiments). NATL2A override: `time_matched`.
- Sampling: triplicate axenic NATL2A starter cultures were split into 6 bottles. Washed Alteromonas was added to three of them at 0 h, and sterile medium was added to the other three. Both sets were sampled at −12, 2, 4, 6, 12, 24 and 48 h (p. 2832 and Results). Alteromonas transcriptomes exist only in the cocultures and were compared with the 12 h point (Table 3; S3 header). The Alteromonas population declined for the first day after addition, "were stable the following day, and then began to grow" (p. 2835); the 12 h point is "when the Alteromonas population stopped declining".
- DE fitting: DESeq2 1.6.3. "Raw counts of reads mapping to each ORF ... across all time points and replicates were first compiled; reads mapping to the Prochlorococcus and A. macleodii MIT1002 genomes were stored in separate tables and analyzed independently ... Results from the biological replicates were used to estimate the dispersion of each gene. Tests for differential expression were carried out between each pairwise condition of interest ... with the Wald test" (p. 2835).

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `coculture_alteromonas_macleodii_mit1002_natl2a_rnaseq` | -12h, 2h, 4h, 6h, 12h, 24h, 48h | time_matched (override) | — | high | Fig. 1 legend "Co-cultured and axenic samples were compared pairwise for each indicated time point" |
| `coculture_24_hours_after_coculturing_mit1002_rnaseq` | 24h vs 12h | baseline (default) | mit1002_rnaseq_deseq2 | high (timing) / medium (model) | S3 header "24 vs 12 hrs after addition"; Table 3 "vs 12 h time point" |
| `coculture_48_hours_after_coculturing_mit1002_rnaseq` | 48h vs 12h | baseline (default) | mit1002_rnaseq_deseq2 | high / medium | S3 header "48 vs 12 hrs after addition" |

Opens / uncertainties:
- `mit1002_rnaseq_deseq2` rests on the methods text: one compiled per-organism count table with a joint dispersion estimate, then Wald tests per pairwise comparison, which describes one DESeqDataSet. The text never says "one design" explicitly. Independent of how DE was fitted, the two experiments share the same 12 h reference libraries, so they are not independent. If the curator rejects the id, the `design_notes` still record the shared reference.
- The growth phase of the 12 h control is written as "unknown (lag-like)", because the vocabulary has no lag value. A curator may prefer `other:lag`.
- The NATL2A override carries an optional `design_notes` flagging that the `-12h` contrast precedes Alteromonas addition, so both arms are axenic at that point.

Other findings (not fixed):
- `growth_phase: exponential` on the MIT1002 analyses: at 24 h the Alteromonas population was "stable" (not growing) and only began to grow afterwards (p. 2835, Fig. 1a). `exponential` describes NATL2A, not the measured MIT1002 arm.
- MIT1002 `timepoint: "24h vs 12h"` puts the contrast into the label. Under the new contract the label should probably be `24h` / `48h`, with the 12 h reference carried in `design_notes`.

## Biller 2018
- Sources read: main PDF (Materials and methods pp. 13–15 in full; Fig. 1 legend; Results "Culture dynamics" p. 3). `legend for supp tables.txt` (S3, S4, S5, S6). CSV headers of S3 (Axenic/Co-culture × 36/40/48 hours) and S6B. GSE93197 count header (A/H × time × L/D × 3). `geo_source.txt`. The supp docx moesm47 turned out to be a Biller 2016-era Rockhopper legend, not relevant.
- Paper default: `time_matched`
- Sampling: three 4 L bottles each of axenic and coculture NATL2A were grown on 13:11 L:D and sampled from mid-exponential at 0–24 h. "After the first 24 h of sampling, the cultures ... were split in half ... Half of these cultures remained in the 13:11 L:D incubator, and the other half were moved into a dark incubator". Both halves were sampled at 28, 32, 36, 40 and 48 h, and the dark half additionally at 72–144 h (p. 14). The count matrix confirms paired L and D samples at 28–48 h.
- DE fitting: DESeq2 1.10.0. "raw counts ... across all time points and replicates were first compiled; reads mapping to the Prochlorococcus and A. macleodii MIT1002 genomes were stored in separate tables and analyzed independently. The standard DESeq2 pipeline was used to normalize ... and estimate gene-level dispersions ... Tests for differential expression were performed between each pairwise condition of interest (i.e., the cocultured and axenic Prochlorococcus samples from a given time point) with the Wald test" (p. 15).

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `darkness_extended_darkness_natl2a_rnaseq_axenic` | 1h/5h/13h extended darkness (36h/40h/48h) | time_matched (default) | natl2a_rnaseq_deseq2 | high (timing) / medium (model) | p. 14 split at 24 h, both halves sampled at 36, 40, 48 h; S3 legend "extended darkness versus 13:11 light:dark controls" |
| `darkness_extended_darkness_natl2a_rnaseq_coculture` | same | time_matched (default) | natl2a_rnaseq_deseq2 | high / medium | same |
| `darkness_extended_darkness_mit1002_rnaseq` | 1h extended darkness, 5h extended darkness | time_matched (default) | — | high | S6 legend "extended darkness versus diel controls"; only one MIT1002 experiment, so no id |

Opens / uncertainties:
- `natl2a_rnaseq_deseq2`: the methods' own example contrast ("the cocultured and axenic Prochlorococcus samples from a given time point") crosses the axenic and coculture sample sets. That implies both are in one per-organism DESeq2 dataset, hence one fit linking the two NATL2A experiments. Medium, because the code is not available. Note also that the axenic libraries are single-end and the coculture libraries paired-end, which would argue for keeping them apart if the authors had been cautious. The plan's prior ("probably separate fits") is therefore contradicted by the text but not disproven.

Other findings (not fixed):
- `timepoint_hours` is inconsistent between experiments of the same design. The NATL2A analyses use the experiment clock (36/40/48), while the MIT1002 analyses use hours in darkness (1.0/5.0) for the same sampling times (36 h / 40 h).

## Coe 2024
- Sources read: main PDF pp. 1–3 (Materials and methods: Cultures, Transcriptomic sampling, RNA-seq statistical analysis) and figure legends (asterisks P < .01). Supp PDF (methods notes, table legends). `supplementary table legend.txt`. CSV headers of S2/S4 (0–24 hours columns). `geo_source.txt` (42 unlabelled samples). README.md.
- Paper default: `time_matched`
- Sampling: dark-tolerant and parental co-cultures were sampled "from exponentially growing cells in standard 13:11 light:dark conditions, 2 days prior to placement into extended darkness on the seventh transfer ... across one complete diel cycle at 0, 4, 8, 16, 20, and 24 h" (p. 2). The figures and tables compare the two phenotypes "at the indicated time points".
- DE fitting: DESeq2 1.24.0. NATL2A and MIT1002 reads were "analyzed separately"; "Tests for differential expression were then performed on each pairwise comparison of interest ... with the Wald test". The two DE experiments are different organisms, so no id is possible.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `darkness_darktolerant_coculture_under_1311_natl2a_rnaseq` | 0h, 4h, 8h, 13h, 16h, 20h, 24h | time_matched (default) | — | high | p. 2 Transcriptomic sampling; legends "between the dark-tolerant and parental phenotype at the indicated time point" |
| `darkness_darktolerant_coculture_under_1311_mit1002_rnaseq` | same | time_matched (default) | — | high | same |

Opens / uncertainties:
- The methods do not say explicitly that parental cultures were sampled in the same session as dark-tolerant ones. A per-timepoint pairwise design over one diel cycle implies it.

Other findings (not fixed):
- The methods list sampling at 0, 4, 8, 16, 20 and 24 h (no 13 h), while S2/S4 and the paperconfig have a `13h` timepoint. The 13 h point (sunset) is probably real and missing from the methods text, but worth checking against GEO GSE264347 sample labels.
- Significance: paper and legends say adjusted P < .01 ("***"). The paperconfig has `pvalue_threshold: 0.1` with `pvalue_asterisk_in_logfc: true`. If the asterisk is what the adapter uses, the threshold value is cosmetic but misleading.
