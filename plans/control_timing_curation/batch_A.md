# Batch A — control-timing curation

## Summary
- Papers: 1 · DE experiments: 10 · baseline: 9 · mixed: 0 · time_matched: 1 · shared-model ids: 2 (`med4_rnaseq_deseq2` on 3 experiments, `hot1a3_rnaseq_deseq2` on 3 experiments)
- Needs human review:
  - `coculture_alteromonas_hot1a3_med4_rnaseq`: curated **`baseline`**, not the planned `mixed`. The analysis code shows the "day 11" contrast is coculture day 11 vs axenic **day 7**, so neither timepoint is time-matched. This breaks the plan acceptance criterion "Weissberg ... reads `mixed`" (Decision D0).
  - `coculture_prochlorococcus_med4_hot1a3_rnaseq`: `time_matched` (Decision D1, high confidence).
  - Proteomics: no `statistical_model_id`. The proteomics contrasts are DEqMS/limma fits separate from the RNA models, and axenic and coculture are fitted separately (Decision D2). `test_type: DESeq2` on the 4 proteomics experiments is wrong (Other findings).
  - `DE_med4_rnaseq_axenic_day_14_vs_exp_day_7` / `DE_med4_proteome_axenic_day_14_vs_exp_day_7`: `growth_phase` should probably read `death` (Decision D3).
  - **Data errors in the shipped DE tables (not fixed, high impact):** 2 "axenic long-term" RNA tables are actually coculture days 60+89 vs axenic. The HOT1A3 axenic proteome "day 31" table actually holds day 89 vs day 11. See Other findings 1 and 2.

## Weissberg 2025
- Sources read:
  - Main PDF (`cc1a3_bioax_2025.pdf`, 22 pp.), all of it: Results pp. 2–12 including the Fig 1/2/3/6 legends and the "Decomposing Prochlorococcus cells" section; Methods pp. 14–15 (Experiment setup, RNA SEQ, Proteomics, DE analysis); Code availability.
  - `README.md`, `MED4_legend.txt`, `HOT1A3_legend.txt`.
  - The `contrast1` column of every shipped DE CSV.
  - **The authors' analysis repo, cloned locally at `~/github/CC1A3`** (the paper cites github.com/wosnat/CC1A3):
    - `RNASEQ/data/CC1A3_RNASEQ_metadata.xlsx` (sample → group, timepoint, day)
    - `RNASEQ/CC1A3 RNASEQ MED4 DESeq.ipynb`, `... MED4 DESeq_longterm2.ipynb`, `... 1A3 DESeq.ipynb`, `... 1A3 DESeq_longterm2.ipynb`
    - `PROTEOME/CC1A3 MED4 proteome DE.ipynb`, `PROTEOME/CC1A3 HOT1A3 proteome DE.ipynb`
    - `Enrichment/CC1A3 MED4 Enrichment.ipynb` and `... HOT1A3 Enrichment path.ipynb` (the contrast-name lookup)
    - `Enrichment/MED4_export_to_kg.ipynb` and `Enrichment/HOT1A3_export_to_kg.ipynb` (the contrast → CSV filename map)
- Paper default: `baseline`
- Sampling:
  - The main experiment ran coculture, axenic MED4 and axenic HOT1A3 in parallel, with omics sampled on days 11, 18, 31, 60 and 89 from separate triplicate bottles per timepoint (Methods p. 14).
  - The axenic MED4 control was **repeated in a second experiment** with omics sampled on days 7 and 14, "because the sampling on days 11, 18 did not match the timing of log exponential growth and decline ... Instead, Omics measurements were taken on days 7, 14 to match med4 axenic growth dynamics" (Methods p. 14).
  - The RNA metadata confirms the mapping: P1 = day 7, P2 = day 14 (axenic MED4); A1/A2/A3 = days 11/18/31 (axenic HOT1A3); C1–C5 = days 11/18/31/60/89 (coculture).
- DE fitting:
  - **RNA:** DESeq2, one fit per organism over all non-technical-replicate samples, `design = ~groupandtime`. MED4 levels are P1, P2, C1–C5; HOT1A3 levels are A1–A3, C1–C5. The pooled "days 60+89" contrasts come from a **second** DESeq2 fit (`~treatment`, with `longterm` = coculture C4+C5) on the same samples.
  - **Proteomics:** DEqMS (limma `lmFit` + `makeContrasts` + `spectraCounteBayes`; Methods p. 15 "Differential Protein abundance was analysed by DEqMS"). Each organism × culture group (axenic, coculture) is fitted separately, and the pooled "days 60+89" contrast is another separate fit.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `growth_state_pro99lown_nutrient_starvation_hot1a3_rnaseq_axenic` | day 18, day 31, days 60+89 | baseline (default) | hot1a3_rnaseq_deseq2 | high | Exported contrasts `A2vsA1_rna`, `A3vsA1_rna` (both in `~groupandtime` fit) and `LATEvsA1_rna` (longterm2 fit); A1 = day 11 (metadata `Day` column) |
| `growth_state_pro99lown_nutrient_starvation_hot1a3_rnaseq_coculture` | day 18, day 31, day 60, day 89, days 60+89 | baseline (default) | hot1a3_rnaseq_deseq2 | high | `C2..C5vsC1_rna` (`~groupandtime`), `LATEvsC1_rna` (longterm2); C1 = day 11 |
| `coculture_prochlorococcus_med4_hot1a3_rnaseq` | day 11, day 18, day 31 | **time_matched (override)** | hot1a3_rnaseq_deseq2 | high | `C1vsA1_rna`, `C2vsA2_rna`, `C3vsA3_rna`. Metadata: A1/A2/A3 Day = 11/18/31, same run and sequencing batch as C1–C3 |
| `growth_state_pro99lown_nutrient_starvation_hot1a3_proteomics_axenic` | day 18, day 31 | baseline (default) | — | high | DEqMS axenic-only fit: `A2 - A1`, `A5 - A1` (see Other finding 2) |
| `growth_state_pro99lown_nutrient_starvation_hot1a3_proteomics_coculture` | day 18, day 31, day 60, day 89, days 60+89 | baseline (default) | — | high | DEqMS coculture-only fit `C2..C5 - C1`; separate fit `Clongterm - C1` |
| `growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic` | day 14, days 60+89 | baseline (default) | med4_rnaseq_deseq2 | high | `P2vsP1_rna` (`~groupandtime`) and `LATEvsP1_rna` (longterm2); P1 = day 7 |
| `growth_state_pro99lown_nutrient_starvation_med4_rnaseq_coculture` | day 18, day 31, day 60, day 89, days 60+89 | baseline (default) | med4_rnaseq_deseq2 | high | `C2..C5vsC1_rna`, `LATEvsC1_rna`; C1 = day 11, "log exponential" |
| `coculture_alteromonas_hot1a3_med4_rnaseq` | day 11, day 18 | **baseline (default; plan expected mixed)** | med4_rnaseq_deseq2 | high on facts / medium on vocab fit | Export map: day-11 file = `C1vsP1_rna` (coculture d11 vs axenic **d7**), day-18 file = `C2vsP2_rna` (coculture d18 vs axenic d14) |
| `growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic` | day 14, day 31, day 89 | baseline (default) | — | high | DEqMS axenic-only fit: `P2 - P1`, `P3 - P1`, `P5 - P1`; P1 = day 7 |
| `growth_state_pro99lown_nutrient_starvation_med4_proteomics_coculture` | day 18, day 31, day 60, day 89, days 60+89 | baseline (default) | — | high | DEqMS coculture-only fit `C2..C5 - C1`; separate fit `Clongterm - C1` |

Shared-model checks:
- **MED4, exact identity:** `C2vsP2 = C2vsC1 + C1vsP1 − P2vsP1` holds on log2FC to 1e-15 over all 1,855 genes (MLE, no shrinkage).
- **HOT1A3, near-exact:** `C2vsA2 = C2vsC1 + C1vsA1 − A2vsA1` (and the same for day 31) holds on 4,012 of 4,015 genes. The 3 exceptions are genes where DESeq2 sets a contrast LFC to 0 because all counts are zero.
- So the coculture-vs-axenic experiments are not independent of the growth-state experiments in the same model.

Opens / uncertainties:
- `coculture_alteromonas_hot1a3_med4_rnaseq`: the vocabulary has no value for "different calendar day and separate culture run, chosen to match growth state". I used `baseline` because by the literal rule no contrast is time-matched, and `mixed` requires at least one that is. See D0.
- `growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic`:
  - P1/P2 (days 7/14) are from the axenic repeat experiment. Nothing states which run the day 31 and day 89 axenic proteomes (P3/P5) came from.
  - The paper says "Axenic Prochlorococcus was sequenced on 4 timepoints: growth (days 7), decline (day 11) and dead cells (days 30 and 90)" (p. 9). The "day 11" there conflicts with the methods and file names (day 14).
  - So these contrasts may span two culture runs. The `control_timing` value is unaffected.
- The HOT1A3 proteome A5 = day 89 mapping is inferred from the index convention, not stated anywhere:
  - RNA metadata: timepoint 5 = day 89.
  - The MED4 export maps P5 → day 89.
  - The proteome notebooks carry no day column.
  - Medium-high confidence.
- Design_notes on the 4 HOT1A3 experiments call the day-11 HOT1A3 control `nutrient_limited`, not `exponential`. The paper says Alteromonas "There is no exponential growth stage. Instead, the cell numbers begin to decline as the experiment starts, likely due to C starvation" (p. 10).

Other findings (not fixed):
1. **Two "axenic long-term" RNA tables contain coculture data.** In both cases `longterm` in the longterm2 notebooks = coculture C4+C5, and no axenic RNA exists after day 14 (MED4) or day 31 (HOT1A3). The paper says "we could not obtain RNA from these samples", p. 6.
   - `MED4__rnaseq__axenic_MED4__longterm_starvation_days_60_and_89__vs__log_expnential_day_7.csv` (analysis `DE_med4_rnaseq_axenic_longterm_vs_exp_day_7`) is DESeq2 `longterm vs P1` = **coculture** days 60+89 vs axenic day 7.
     - Numerical confirmation: `LATEvsP1 − LATEvsC1 ≈ C1vsP1` (max abs diff 0.013).
   - `HOT1A3__rnaseq__axenic_HOT1A3__longterm_starvation_days_60_and_89__vs__day_11.csv` (`DE_hot1a3_rnaseq_axenic_longterm_vs_day_11`) is `longterm vs A1` = **coculture** days 60+89 vs axenic day 11.
   - Both are currently attached to the `_axenic` experiments with `background_factors: [axenic, ...]`. They should probably be dropped, or moved to a coculture-vs-axenic experiment.
2. **The HOT1A3 axenic proteome "day 31" table is day 89.**
   - In `HOT1A3_export_to_kg.ipynb`, `contrast_map` maps both `A3vsA1_prot` and `A5vsA1_prot` to `HOT1A3__proteome__axenic_HOT1A3__day_31__vs__day_11`. The later write wins.
   - The shipped CSV's `contrast1` column is `A5vsA1_prot` for all 2,229 rows.
   - So `DE_hot1a3_proteome_axenic_day_31_vs_day_11` (timepoint day 31) carries day-89-vs-day-11 data, and the real day 31 contrast was never exported.
   - Fix upstream in CC1A3 and re-export.
3. **The MED4 "day 11 coculture vs axenic" table is not time-matched.**
   - `MED4__rnaseq__log_expnential_day_11__coculture_with_HOT1A3__vs__axenic.csv` = `C1vsP1_rna`, i.e. coculture day 11 vs axenic day 7 (separate experiment run).
   - The filename, legend ("Both samples are using log expnential day 11") and the analysis id `..._exp_day_11` all imply day 11 for both arms.
   - Captured in design_notes. The filename and legend themselves are misleading.
4. **`test_type: DESeq2` is wrong on the 4 proteomics experiments.**
   - Affected: `growth_state_..._hot1a3_proteomics_axenic`, `..._hot1a3_proteomics_coculture`, `..._med4_proteomics_axenic`, `..._med4_proteomics_coculture`.
   - Methods p. 15 and the notebooks use **DEqMS** (limma + `spectraCounteBayes`; the p-values are `sca.P.Value` / `sca.adj.pval`).
   - README.md also says "proteomics via DESeq2", which is wrong.
5. **The HOT1A3 control condition is mislabelled "exponential growth".**
   - Affected: `control_condition: PRO99-lowN exponential growth` on the 4 HOT1A3 growth_state experiments, and `growth_phase: exponential` on `DE_hot1a3_rnaseq_coculture_vs_axenic_day_11`.
   - The paper says HOT1A3 has no exponential phase in this medium and declines from the start (p. 10, Fig 6).
   - Better: `nutrient_limited`, or `other:carbon_starved_decline`. For the MED4 day-11 analysis, `exponential` is correct.
6. **The growth_phase on axenic MED4 day 14 should probably be `death`.**
   - Affected: `DE_med4_rnaseq_axenic_day_14_vs_exp_day_7` and `DE_med4_proteome_axenic_day_14_vs_exp_day_7`, currently `nutrient_limited` (D3).
7. **The paper text uses "day 11" for the axenic MED4 decline sample in two places** ("we compared dying cultures (day 11) to exponential growth (day 7)", p. 5; "decline (day 11)", p. 9). The methods, metadata and files all say day 14. This is a manuscript erratum, not a KG issue.
8. Minor upstream code bug, not affecting the KG: in `CC1A3 HOT1A3 proteome DE.ipynb` cell 65, `C3 - C2`, `C4 - C3` and `C5 - C4` all use `coef_col = 4` (the `C5 - C1` coefficient). None of these are exported to the KG.
9. The paperconfig has no `doi:`. The DOI comes from PDF extraction (10.1101/2025.11.24.690089). This is fine and unchanged.

## Decisions for the curator

**D0 (new, blocks the plan acceptance criterion) — `coculture_alteromonas_hot1a3_med4_rnaseq`: `baseline` or `mixed`?**
- Evidence: both MED4 coculture-vs-axenic contrasts use the axenic **repeat** experiment as control.
  - day 11 = `C1vsP1` = coculture day 11 vs axenic day 7 (both exponential).
  - day 18 = `C2vsP2` = coculture day 18 vs axenic day 14 (axenic in death).
  - Sources: `Enrichment/MED4_export_to_kg.ipynb` `contrast_map`, the shipped CSV `contrast1` column, and the metadata Day column (P1 = 7, P2 = 14).
- The design intent was growth-state matching ("to match med4 axenic growth dynamics", Methods p. 14). That holds for day 11, where both arms are log-exponential. It fails for day 18: the coculture is in slow decline, while the axenic culture lost viability around day 11.
- My curation: `baseline`, the literal vocabulary rule ("every contrast compares against a control sampled at another time"). Design_notes spell out both controls.
- Options:
  - (a) Keep `baseline` and change the plan acceptance criterion to "reads `baseline`".
  - (b) Call day 11 time-matched in the "same growth state" sense of the `time_matched` definition ("acclimated cultures ... harvested at the same growth state") and use `mixed`. This stretches that clause: these are batch cultures from two different experiment runs, not acclimated steady states.
  - (c) Extend the vocabulary (e.g. `phase_matched`).
- **Recommendation: (a).** Calendar day and culture run both differ, which is exactly what the field should warn about. Note also that the phrase "day 11 matched" in the plan's Known-cases table is wrong.

**D1 — HOT1A3-side `coculture_prochlorococcus_med4_hot1a3_rnaseq` control days: `time_matched`.**
- Evidence:
  - `RNASEQ/data/CC1A3_RNASEQ_metadata.xlsx` lists axenic HOT1A3 samples A1A–A3C with Day = 11, 18, 31. These are in the same sequencing batch (TGC #7–15) as the coculture C1–C5 (TGC #16–34), and come from the main experiment, not the repeat.
  - Methods p. 14: the main experiment sampled all cultures on days 11, 18, 31, 60, 89. Only the axenic *MED4* control was repeated.
  - `CC1A3 RNASEQ 1A3 DESeq.ipynb` contrasts `C1 vs A1`, `C2 vs A2`, `C3 vs A3` → exported as the day 11/18/31 coculture-vs-axenic files (`HOT1A3_export_to_kg.ipynb`).
- **Recommendation: `time_matched` (applied as an override), high confidence.**

**D2 — Are the proteomics contrasts separate fits or shared with the RNA models? Separate. No `statistical_model_id` on any proteomics experiment.**
- Evidence:
  - Methods p. 15: "Differential gene expressions were analysed ... by DESeq2 ... Differential Protein abundance was analysed by DEqMS".
  - `PROTEOME/CC1A3 MED4 proteome DE.ipynb` and `... HOT1A3 proteome DE.ipynb` fit limma/DEqMS separately on the axenic samples (`~0+P` / `~0+A`) and on the coculture samples (`~0+C`), plus separate pooled-longterm fits.
  - So there is no model shared between the axenic and coculture proteomics experiments of one organism, and none with RNA.
  - Within each proteomics experiment, the per-day contrasts come from one fit. All of those contrasts sit inside that single experiment, so there is nothing to link.
- **Recommendation: leave unset (applied).** Separately, fix `test_type` to `DEqMS` (Other finding 4).

**D3 — Should the `nutrient_limited` growth_phase on `DE_med4_rnaseq_axenic_day_14_vs_exp_day_7` (and the proteome twin `DE_med4_proteome_axenic_day_14_vs_exp_day_7`) read `death`?**
- Evidence for `death`:
  - Results p. 5: "The axenic culture died around day 11, therefore we compared dying cultures ... to exponential growth (day 7)".
  - Fig 3 legend (p. 6): panels C/D are "axenic death stage", and the only axenic MED4 RNA contrast is day 14 vs day 7.
  - Fig 4 legend: "Pathways enriched when comparing death ... to log exponential".
  - Viability (p. 3): axenic MED4 "does not revive when transferred into fresh media on day 11".
- Evidence against:
  - Methods p. 14 calls the day 14 sample "decline".
  - The notebooks code timepoint 2 as `phase = 'decline'`.
  - It is the onset of nitrogen starvation, so `nutrient_limited` is not wrong biologically.
- The vocabulary has no `decline`. `death` is in it, and the paper's own framing of these contrasts is death.
- The day-18 design_notes already describe the axenic day 14 control as death.
- **Recommendation: `death`** for both analyses. That is consistent with the axenic proteome day 31/89 analyses, which are already `death`. Not changed (out of scope).
