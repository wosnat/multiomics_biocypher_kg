# Batch F — control-timing curation

## Summary
- Papers: 5 · DE experiments: 18 · baseline: 3 · mixed: 1 · time_matched: 14 · shared-model ids: 0
- Needs human review:
  - thompson 2011 — paper default is `baseline` (3 of 5 DE experiments), which is unusual for a -Fe vs +Fe design; check you agree with "majority wins" here rather than defaulting to `time_matched` and overriding 3.
  - tolonen 2006 `growth_medium_growth_on_*` (4 experiments) — `time_matched` is medium confidence: the methods do not say both arms were harvested on the same day (noted in `design_notes`).
  - tolonen 2006 `nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray` — `mixed`; the paper says only that *one* of the +NH4 0h arrays was replaced by a 3h array, so 0h is partly mismatched. I counted it as mismatched in the notes; the label would be `mixed` either way.
  - Capovilla 2023 — both DE analyses are labelled `day 1 and day 3`, but the paper's Fig. 5B and Results suggest the sd02 tables are the day-1 (24 h) contrast. This doesn't change `control_timing` (time-matched either way). Listed under Other findings.

## thompson 2011
- Sources read: main PDF pp. 1582–1583 (Methods: microarray experiment, sampling, array normalization and analysis, Cyber-T); `supp table legend.txt`; README.md; headers of `supp table 1.csv` / `supp table 2.csv` (`Gene,Description,Cluster,0,12,24,48,R*` and `…,0,16,28,53,R*`).
- Paper default: `baseline` (3 of 5 DE experiments; the 2 time-matched ones are overrides)
- Sampling: mid-log +Fe cultures were washed and split into +Fe and no-added-iron (-Fe) media. RNA was sampled at 0/12/24/48/70 h (MED4) and 0/16/28/53/72 h (MIT9313). Iron was added back to the -Fe cultures at 49 h (MED4) and 54 h (MIT9313).
- DE fitting: Cyber-T (Bayesian t-test), run per pairwise comparison (five comparisons per strain), with q-values from Storey's Q-VALUE. There is no shared model, so no `statistical_model_id`.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `iron_iron_starvation_med4_med4_microarray_light` | 12h, 24h | time_matched (override) | — | high | p. 1582: "We compared no added-iron (-Fe) to iron-replete (+Fe) expression at each of the four iron-stress time points" |
| `iron_iron_starvation_med4_med4_microarray_light_2` | 48h | baseline (default) | — | high | p. 1582: "Owing to poor hybridization for a MED4 48-h iron-replete sample, 48-h no added-iron samples were compared with 24-h iron-replete samples" |
| `growth_medium_ironreplete_control_med4_med4_microarray` | R (rescue: 70h vs 48h) | baseline (default) | — | high | p. 1582–3: "A fifth comparison (R) was made between -Fe cultures before (48-h MED4 …) and after (70-h MED4 …) iron rescue" |
| `iron_iron_starvation_mit9313_mit9313_microarray` | 16h, 28h, 53h | time_matched (override) | — | high | same "at each of the four iron-stress time points" sentence; no substitution is reported for MIT9313 |
| `growth_medium_ironreplete_control_mit9313_mit9313_microarray` | R (rescue: 72h vs 53h) | baseline (default) | — | high | p. 1583: before (53-h) vs after (72-h) rescue |

Opens / uncertainties:
- None on timing. The only judgement call is the paper default (see the review list above).

Other findings (not fixed):
- The two rescue experiments (`growth_medium_ironreplete_control_*`) have `treatment_condition: Iron-replete control` and `control_condition: Iron starvation`. In fact both arms are the *same* -Fe cultures, after vs before iron addition. The names read as if the comparison were the +Fe control culture against the -Fe culture, which it is not.
- The 0h columns (`0`) in both supp tables are not integrated. They are -Fe vs +Fe right after resuspension, time-matched. This is fine, just noting it.

## tolonen 2006
- Sources read: main PDF Methods (Cell culture, Data analysis, pp. 12–13), Fig. 1 legend, Results "Physiological response"; `legend of supp tables.txt`; `medExpression.1.csv` header (`GENE,t=0h,t=3h,t=6h,t=12h,t=24h,t=48h,Cyanate,Urea,ANNOTATION`).
- Paper default: `time_matched`
- Sampling:
  - N deprivation: mid-log cultures were centrifuged and resuspended in +NH4 or -N Pro99 at t=0, then sampled at 0, 3, 6, 12, 24 and 48 h, in duplicate.
  - Alternative N sources: log-phase cultures were established on urea, cyanate or nitrite, alongside cultures on ammonium. This is a steady-state comparison.
- DE fitting: Goldenspike (Choe et al. 2005) q-values for each -N/+NH4 or alternate-N/NH4 comparison. This is a per-comparison test, with no shared model described.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `nitrogen_nitrogen_deprivation_med4_med4_microarray` | 0h–48h | time_matched (default) | — | high | Methods: resuspended in +NH4 or -N; "Samples were taken at 0, 3, 6, 12, 24, and 48 h"; the hybridization substitutions are MIT9313-only |
| `growth_medium_growth_on_cyanate_as_med4_microarray` | steady state | time_matched (default) | — | medium | Methods: "log-phase cultures … were established in Pro99 medium containing one of the following nitrogen sources" |
| `growth_medium_growth_on_urea_as_med4_microarray` | steady state | time_matched (default) | — | medium | same |
| `nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray` | 0h–48h | mixed (override) | — | high | Data analysis: "Because of microarray hybridization problems with the MIT9313 control (+NH4) samples at t=24 h, the -N expression values at this time point were compared to the +N at t=12 h. Similarly, a MIT9313 +NH4 array at t=3 h was substituted for one of the +NH4 arrays at t=0 h" |
| `growth_medium_growth_on_nitrite_as_mit9313_microarray` | steady state | time_matched (default) | — | medium | as for the cyanate experiment |
| `growth_medium_growth_on_urea_as_mit9313_microarray` | steady state | time_matched (default) | — | medium | as for the cyanate experiment |

Opens / uncertainties:
- `growth_medium_growth_on_*` (4 experiments): harvest timing of the ammonium reference is not stated. I chose `time_matched` because both arms are acclimated log-phase cultures, and recorded this in `design_notes`. It is also unstated whether the ammonium reference arrays are the same arrays as the +NH4 arm of the deprivation series. ArrayExpress E-TABM-91 sample metadata would resolve both points.
- `nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray`: 0h is only partly mismatched (1 of the 2 +NH4 replicates came from 3h). The authors argue the impact is small because correlations between +NH4 time points are as high as between replicates.

Other findings (not fixed):
- The `gene_clusters` entries' `experimental_context` says "27 µmol photons m-2 s-1, 24°C". The paper's Methods say 22 °C, with 50 (MED4) or 10 (MIT9313) µmol quanta m-2 s-1, which matches the experiment blocks. The cluster text looks copied from another paper.
- Every N-deprivation analysis has `growth_phase: acute_stress`, including 24h and 48h. In Fig. 1, fluorescence declines from 12–24 h, so `nutrient_limited` may fit 24h and 48h better.

## Read 2017
- Sources read: main PDF pp. 2268–2269 (Methods: cell cultures, RNA sequencing overview, Rockhopper/DESeq); README.md; header of `DE_genes_3h_post_starvation.csv`.
- Paper default: `time_matched`
- Sampling: triplicate mid-exponential cultures were washed and split into N-replete or N-depleted Pro99, then sampled at 0, 3, 12 and 24 h after resuspension.
- DE fitting: Rockhopper (DESeq negative binomial), run per time point. Only one experiment, so no `statistical_model_id`.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq` | 3h, 12h, 24h | time_matched (default) | — | high | p. 2269: "Differential expression was calculated … between the experimental (N-deprived) cells and the control (N-replete cells) at each time point (0, 3, 12, 24 h post N removal)" |

Opens / uncertainties: none.

Other findings (not fixed):
- `adjusted_p_value_col: p_value`. The Methods say significance was a BH-adjusted P < 0.05 (Rockhopper reports q-values), so the mapping looks right. The README's "verify" item can be closed.

## Thompson 2016
- Sources read: main PDF pp. 3–4 (Methods: growth and infection in light and dark; detection of differential gene expression); Results p. 7 (four-condition design); `geo_source.txt` and the GSE79359 sample list (Ctr/Inf × Dk/Lt × 0,1,2,4,8 h × 2 replicates); README.md.
- Paper default: `time_matched`
- Sampling:
  - Light and dark experiments were run on separate days. In each, a mid-log MED4 mother culture was split into 4 bottles: 2 got phage P-HM2 and 2 got spent medium.
  - All bottles were sampled together at 0, 1, 2, 4 and 8 h, reported as 0.5–8.5 h post-inoculation.
- DE fitting: DESeq2 plus NOISeq, run per time point and per pairwise contrast. The paper says each time point was analysed separately, so there is no shared model.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `light_stress_constant_dark_med4_rnaseq_infected` | 0.5h–8.5h post-inoculation | time_matched (default) | — | high | p. 4: "Transcript abundances were analyzed at each time point separately … comparing dark to light treatments with infected or uninfected constant"; p. 3: "Light and dark experiments were conducted on separate days" |
| `light_stress_constant_dark_med4_rnaseq_dark` | 1.5h–8.5h post-inoculation | time_matched (default) | — | high | same |
| `viral_phage_phm2_lysate_med4_rnaseq_light` | 0.5h–8.5h post-inoculation | time_matched (default) | — | high | p. 4: "comparing infected to uninfected treatments with light or dark constant"; p. 3: spent-medium bottles split from the same mother culture |
| `viral_phage_phm2_lysate_med4_rnaseq_dark` | 2.5h–8.5h post-inoculation | time_matched (default) | — | high | same |

Opens / uncertainties:
- The dark-vs-light experiments compare arms from two different days and mother cultures, matched on hours post-inoculation. I treated this as `time_matched` and recorded it in `design_notes`. The paper also says it chose not to use a t=0 reference ("It was not appropriate to compare each time point to a t = 0 reference").

Other findings (not fixed):
- `light_stress_constant_dark_med4_rnaseq_dark` (the uninfected spent-medium experiment) has `background_factors: [axenic]` with no light or dark factor. Its key suffix `_dark` is also confusing, because the experiment is "uninfected", not "dark".

## Capovilla 2023 (2 DE experiments of 4)
- Sources read: main PDF p. 9 Methods ("RNA-Seq and Metabolites Analysis", "RNA-Seq Analysis"); Fig. 5 legend and Results paragraph (p. 6–7); sd02 xlsx sheet headers; sample/condition columns of `read counts … sd04.csv` / `sd05.csv` (AMP1 vs AMP1_chitosan × day_1 / day_3, triplicates).
- Paper default: `time_matched`
- Sampling: replicate AMP1 cultures were split, and half received 56 µg/mL chitosan. RNA was collected on day 1 and day 3 from both conditions, in triplicate.
- DE fitting: DESeq2 Wald test, with MIT9303 and MIT9313 "analyzed separately". The contrasts are chitosan vs control on day 1 and on day 3. The two DE experiments are different organisms, so no `statistical_model_id` is possible. The metabolomics experiments are not DE and got no field.

| Experiment | timepoints | control_timing | statistical_model_id | confidence | Evidence |
|---|---|---|---|---|---|
| `chitosan_addition_mit9313_rnaseq` | day 1 and day 3 | time_matched (default) | — | high | p. 9: "Differential expression tests were performed on comparisons of cultures with chitosan vs. control cultures 1 and 3 d after chitosan addition"; "collected … in triplicates at each time point and each condition" |
| `chitosan_addition_mit9303_rnaseq` | day 1 and day 3 | time_matched (default) | — | high | same |

Opens / uncertainties: none on timing.

Other findings (not fixed):
- Each strain's single analysis is labelled `timepoint: day 1 and day 3`, but the sd02 sheets have no day column. The Fig. 5B legend says the tables show expression "24 h after addition of chitosan", and Results say significant changes "disappeared after 48 h". The sd02 table is therefore probably the day-1 contrast only (timepoint `day 1`).
- The `experimental_context` contradicts itself ("RNA-Seq 24h after chitosan addition" and "comparisons … on days 1 and 3").
- The Methods say day 1/day 3 "after inoculation" in one sentence and "after chitosan addition" in another.
