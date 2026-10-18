# Batch B — control-timing curation

## Summary
- Papers: 1 · DE experiments: 60 · baseline: 0 · mixed: 0 · time_matched: 60 · shared-model ids: 0
- Needs human review: none for control_timing / statistical_model_id (both stated explicitly in Methods). Two "Other findings" below (growth_phase/timepoint "stationary", light_intensity) are worth a look but were not changed.

## Moreno 2023
- Sources read: main PDF pp. 1-24 (Results pp. 2-13 incl. Fig. 1 legend p. 4, Table 1 footnote p. 7; Materials and Methods pp. 22-24 in full); Supp PDF s0005 (figures S1-S3 only, no methods); `legend.txt` (Supp Table S1-S4 column legends); `README.md`; CSV headers of `table s2 MED4 ..._modified.csv` and `table s4 Marinobacter in MED4 ..._modified.csv`; the paperconfig (60 experiments, one analysis each).
- Paper default: `time_matched` (added as `publication.control_timing`, right after `papermainpdf`; no `doi` key in this file). No overrides.
- Sampling: p. 22 Methods, "When they reached 0.05 units of absorbance ... cultures were divided into 6 aliquots and subjected to the different experimental conditions: light, light + 100 nM glucose, light + 5 mM glucose, darkness, darkness + 100 nM glucose, and darkness + 5 mM glucose ... After 24 h, cells were harvested by centrifugation". Every arm (glucose-treated and its no-glucose control, in light and in dark) is an aliquot of one parent culture, harvested together 24 h after the split. The Alteromonas (S3) and Marinobacter (S4) proteomes come from the same 18 samples (the heterotroph fraction of the same non-axenic cultures; p. 23 "Proteomic data analysis", MarRef v6 search), so the same timing holds for them. Fig. 1 legend p. 4: "Three biological replicates were prepared per condition."
- DE fitting: label-free quantification in Progenesis QI; per-comparison Student t-test between two groups (p. 23: "For a protein to be considered significantly differentially expressed, it had to be identified and quantified using at least two unique peptides and a P value of <=0.05"; `legend.txt`: "TTEST.L.LGnM, Student's t-test of results from cultures in the light vs cultures in the light with 100 nM glucose", one TTEST column per pair). No shared fitted model -> no `statistical_model_id`.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| Cyano S2 (20): `light_low_glucose_{med4,ss120,wh7803,wh8102,bl107}_proteomics`, `light_high_glucose_{med4,ss120,wh7803,wh8102,bl107}_proteomics` | stationary | time_matched (default) | — | high | p. 22: cultures "divided into 6 aliquots ... After 24 h, cells were harvested"; control = light aliquot, no glucose, same harvest |
| Cyano S2: `dark_low_glucose_{med4,ss120,wh7803,wh8102,bl107}_proteomics`, `dark_high_glucose_{med4,ss120,wh7803,wh8102,bl107}_proteomics` | darkness | time_matched (default) | — | high | same; control = darkness aliquot, no glucose, same harvest |
| Marinobacter S4 (20): `light_low_glucose_marino_in_{med4,ss120,wh7803,wh8102,bl107}_proteomics`, `light_high_glucose_marino_in_{med4,ss120,wh7803,wh8102,bl107}_proteomics` | stationary | time_matched (default) | — | high | same samples as S2 (heterotroph fraction, p. 23 MarRef v6 search); per-pair t-test (`legend.txt`, Table S4) |
| Marinobacter S4: `dark_low_glucose_marino_in_{med4,ss120,wh7803,wh8102,bl107}_proteomics`, `dark_high_glucose_marino_in_{med4,ss120,wh7803,wh8102,bl107}_proteomics` | darkness | time_matched (default) | — | high | same |
| Alteromonas S3 (20): `light_low_glucose_alt_in_{med4,ss120,wh7803,wh8102,bl107}_proteomics`, `light_high_glucose_alt_in_{med4,ss120,wh7803,wh8102,bl107}_proteomics` | stationary | time_matched (default) | — | high | same samples as S2; per-pair t-test (`legend.txt`, Table S3) |
| Alteromonas S3: `dark_low_glucose_alt_in_{med4,ss120,wh7803,wh8102,bl107}_proteomics`, `dark_high_glucose_alt_in_{med4,ss120,wh7803,wh8102,bl107}_proteomics` | darkness | time_matched (default) | — | high | same |

All 60 DE experiment keys are covered by the brace expansions above (3 organism sides x 5 strains/cultures x {light,dark} x {low,high} glucose). No `design_notes` added: none is required (`time_matched`, no model id), and the default needs no qualification.

Opens / uncertainties:
- None on the two curated fields. Note only: the t-tests are pairwise, but all four comparisons per strain share the same control replicates (the light control is the denominator for both light_low and light_high; likewise dark). That is shared samples, not a shared fitted model, so per the rules no `statistical_model_id` was set.

Other findings (not fixed):
- `growth_phase: stationary` / `timepoint: "stationary"` on the 30 light analyses, and `experimental_context` "harvested at stationary growth" (MED4 light_low): the paper does not say stationary. Cultures were split at A674 0.05 (Pro) / A550 0.1 (Syn) and harvested 24 h later (p. 22) — that reads as early/mid exponential, not stationary. Likely should be `exponential` (or `unknown`) with a timepoint like "24 h" / `timepoint_hours: 24`.
- `growth_phase: darkness` / `timepoint: "darkness"` on the 30 dark analyses: timepoint is a condition label, not a time; the sampling time is 24 h after transfer to darkness (`timepoint_hours: 24` would be accurate; currently null).
- `light_intensity: "15 umol photons m-2 s-1"` on all light experiments disagrees with p. 22: Prochlorococcus 4 uE m-2 s-1 for SS120 (low-light ecotype) and 40 uE m-2 s-1 for MED4; Synechococcus 40 uE m-2 s-1. The heterotroph-side experiments inherit their cyano culture value.
- `adjusted_p_value_col` points at the raw Student t-test P value columns (`TTEST.*`, `t.test ...`); the paper applies no multiple-testing correction to the proteomics (only to DAVID enrichment and metabolomics q values). Values carried as `adjusted_p_value` are therefore unadjusted.
- Cyano S2 CSV header labels the dark t-tests `TTEST.D.LGnM` / `TTEST.D.LGmM` (with "L"), while `legend.txt` calls them `TTEST.D.DGnM` / `TTEST.D.DGmM`. The paperconfig maps them to the dark experiments; worth a one-time check that these columns really are the dark comparisons (the per-sample columns D1-3 / DGnM1-3 / DGmM1-3 exist, so a recomputation would settle it).
