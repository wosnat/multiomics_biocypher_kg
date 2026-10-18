# Fix batch E

Validation: `scripts/validate_paperconfig.py` ends in `VALIDATION PASSED` on all 9 paperconfigs. No new
warnings: Capovilla's 4 metabolomics warnings and Tai's 2 `Gene Name` warnings were already there. Nothing was
committed, and `prepare_data.sh` / Docker were not run.

## Aharonovich 2016
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| E1 linear FC | The two HOT1A3 `logfc_col`s hold linear ratios | Yes. From the CSVs: `fold change 9313/9313DIL` has min 0.002, median 3.71, max 282,301, 0 negatives, 422 of 1,966 values < 1. `fold change 9313+1A3/med4+1A3` has min 0.0002, median 1.20, max 4,801, 0 negatives, 874 of 2,381 values < 1. The values are the ratio of the `+.01` expression columns (row 1: 20.01/0.01 = 2001; row 2: 34.01/71.01 = 0.479). | Added `fold_change_type: linear` to `DE_hot1a3_coculture_9313_vs_9313dil_t20h` and `DE_hot1a3_coculture_9313_vs_med4_t20h` |
| HOT1A3 conditions | Both control_conditions say "Pro99 medium growth conditions", but the denominators are cocultures | Yes. CSV headers are "Expression 9313+1A3 / 9313+1A3 dil" and "9313+1A3 / med 4". Fig. 2a legend (p. 2895): high and low HOT1A3 densities are "5 × 10^6 and 0.5 × 10^6 cells per ml"; MED4+1A3 used 5 × 10^6 (p. 2896). | `coculture_prochlorococcus_mit9313_high_alteromonas_hot1a3_rnaseq`: control "Pro99 medium growth conditions" → "Coculture with Prochlorococcus MIT9313 (low Alteromonas inoculum, 0.5x10^6 cells/ml)". Name "... vs Pro99 medium growth conditions" → "... vs coculture with MIT9313 at low Alteromonas inoculum (0.5x10^6 cells/ml) (RNASEQ)". experimental_context gained "; both arms are 20 h cocultures with MIT9313 that differ only in the initial HOT1A3 inoculum (no axenic HOT1A3 arm)". `coculture_prochlorococcus_mit9313_hot1a3_rnaseq`: control → "Coculture with Prochlorococcus MED4". Name → "HOT1A3 Coculture with Prochlorococcus MIT9313 vs coculture with Prochlorococcus MED4 (RNASEQ)". experimental_context gained the same-inoculum / partner-strain note. |
| DE criterion | Paper uses FC > 2, q < 0.05, expression > 10; config has pvalue_threshold 0.01 | Yes, p. 2894: "fold change in expression of co-culture compared with the axenic culture was >2, the q-value was <0.05 and the expression level in either axenic culture or co-culture was >10". This criterion covers the three Prochlorococcus coculture-vs-axenic tables. | `pvalue_threshold` 0.01 → 0.05 on the 3 Prochlorococcus analyses (MED4, MIT9313 5e6, MIT9313 0.5e6). FC > 2 was already enforced by the adapter default `logfc_threshold` 1.0, so no field was added. Expression > 10 cannot be expressed in a paperconfig. The two HOT1A3 tables keep 0.01, because the paper gives no criterion for them. |
| MED4 headers | The expression column headers look swapped | Yes. Row 1 of PMM0735: log2(56.01/0.01) = 12.45 = the log2FC, so the numerator is the column headed "med4 t=20" (the axenic one). log2FC agrees with the paper. | Report only; no change |

Not applied / needs the curator:
- **Adapter bug, which affects E1.** `OMICSAdapter._check_significance` receives the raw, unconverted value (`fc_float`) and tests `abs(fc) >= logfc_threshold` (default 1.0) on it. With `fold_change_type: linear`, a down-regulated gene (raw ratio < 1) can therefore never be called significant, and a raw ratio of 1–2 passes the FC filter. The edges are still better than before: before the fix every value was read as log2, so all HOT1A3 edges were "up"; now down-regulated genes get a negative log2FC but read `not_significant`. The fix is to pass the converted log2 value in `omics_adapter.py` around line 851. The same applies to Beliaev 2014, which also sets `fold_change_type: linear`.
- The two HOT1A3 experiments keep `treatment_type: [coculture]` and `treatment_organism: Prochlorococcus MIT9313`, which produces a `treatment` partner edge. In the 9313-vs-MED4 experiment, MED4 is a partner in the control arm and is not represented. Changing the partner/role framing was out of scope.
- Methods give 20 °C and 27 µE for maintenance, while the config uses 24C and 10 µE. These may be the experiment conditions from the Supplementary information, which is not in the directory. Left unchanged.

Graph-visible changes: the 2 HOT1A3 experiments get new `name`, `control`, `experimental_context`. Their HOT1A3 edges get `log2_fold_change` = log2(ratio) instead of the raw ratio, with direction now signed (`expression_direction` and `expression_status` change; the significant set is unchanged, see the bug above). On the 3 Prochlorococcus tables, `significant` / `expression_status` change for rows with 0.01 < q ≤ 0.05.
Snapshot entries affected: none. The only snapshot entry is the Publication node (`experiment_count` 5, organisms, treatment_types), which these changes leave unchanged.
CHANGELOG line: Aharonovich 2016: the two HOT1A3 tables were loaded as log2 when they hold linear ratios, so every HOT1A3 gene read as up-regulated and the control was "Pro99 medium". The KG now stores log2(ratio) and names the real comparator cocultures (MIT9313 at low HOT1A3 inoculum; MED4). The Prochlorococcus significance cutoff is now the paper's q < 0.05, not 0.01.

## Tetu 2019
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| PVC control_timing | Curator decision: set the 2 PVC experiments to matched | The decision is the curator's. The note text was checked: the paper says "after short-term exposure (90 min for PVC, 120 min for HDPE)" (Methods, RNA isolation) and never gives the control harvest time. GEO GSE118155 says PVC was harvested at 90 min and control and HDPE at 180 min. | Removed `control_timing: baseline` from `plastic_pvc_plastic_leachate_2_mit9312_rnaseq` and `plastic_pvc_plastic_leachate_2_natl2a_rnaseq`, so they inherit the paper default `matched`. Rewrote both `design_notes`: "90 min: PVC-exposed cultures were harvested after 90 min and the shared AMP1 no-leachate control after 180 min (GEO GSE118155 protocol; the paper does not state the control harvest time). Treated as matched: treatment and control are parallel exponential cultures and the exposure length is part of the treatment. The same three control libraries serve as control for the HDPE experiment of this strain." |
| HDPE 120 vs 180 min | Paper and GEO disagree | The paper is internally consistent. The single statement in the paper is the Methods quote above (120 min for HDPE), and it does not mention 180 min anywhere. Only GEO says 180. | None. The paper text supports 120 min, so `timepoint: 120 min` is kept. |

Not applied / needs the curator: HDPE 120 (paper) vs 180 min (GEO) remains an unresolved source conflict. The config follows the paper.
Graph-visible changes: none, except that the `control_timing` value derived for the 2 PVC experiments becomes `matched`.
Snapshot entries affected: none (no Tetu entries).
CHANGELOG line: Tetu 2019: the PVC-leachate experiments were marked as compared against a later-harvested (baseline) control. They are now recorded as matched, with the 90 vs 180 min GEO harvest offset kept in the design note.

## Capovilla 2023
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| sd02 timepoint | The sd02 tables are the day-1 contrast, not "day 1 and day 3" | Yes. Fig. 5B legend: "relative expression of genes ... 24 h after addition of chitosan". Results: changes "disappeared after 48 h". Decisive check: chitosan/control log2 ratios were recomputed from the sd04/sd05 normalized counts. Sign agreement with sd02 is 100% for day 1 in both strains (MIT9313 81 genes r = 0.88; MIT9303 36 genes r = 0.44). For day 3 it is 37% and 50% (r = -0.25 and 0.02). | Both analyses: `timepoint: day 1 and day 3` → `day 1`; `timepoint_hours: null` → `24` |
| experimental_context | Contradicts itself | Yes | "RNA-Seq 24h after chitosan addition. DESeq2 comparisons of chitosan vs control cultures on days 1 and 3 after inoculation." → "RNA-Seq 24h (day 1) after chitosan addition. The DE table is the DESeq2 chitosan vs control contrast on day 1; day 3 was also sequenced, but significant changes disappeared after 48 h and are not in the table." (both RNA-seq experiments) |

Not applied / needs the curator: none.
Graph-visible changes: `time_point` "day 1 and day 3" → "day 1" and `time_point_hours` → 24.0 on all Capovilla DE edges. This flows into the Experiment `time_point_labels` / `time_point_hours` arrays and the `experimental_context` of the 2 RNA-seq experiments.
Snapshot entries affected: none. The 8 Capovilla snapshot entries are all metabolomics/assay binding edges with no properties.
CHANGELOG line: Capovilla 2023: the chitosan DE tables were labelled "day 1 and day 3" with no hours. The KG now labels them day 1 (24 h after chitosan), which a recomputation from the published counts confirms.

## Coe 2024
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| 13 h timepoint | Absent from the methods sampling list | Partly confirmed. The Methods (p. 2) list "0, 4, 8, 16, 20, and 24 h"; the S2–S4 tables have a "13 hours" column. GEO sample labels are not available offline (`GSE264347_raw_counts.tsv.gz` has unlabelled `sample_1..42`). 42 samples equals 7 timepoints × 2 phenotypes × 3 replicates, which fits 13 h being real (6 timepoints would give 36). | None. The table column is real data; the methods omission looks like a text error. |
| pvalue_threshold 0.1 | The paper uses adjusted P < .01 | Yes. Methods: "genes with an adjusted P-value of <.01 were considered to be significantly differentially abundant". The supp legend says "***" means p < 0.01. Adapter: with `pvalue_asterisk_in_logfc: true`, significance comes from the asterisk only, so 0.1 is **not** the effective cutoff. However, the adapter writes `pvalue_threshold` as the placeholder `adjusted_p_value` on every significant edge (`omics_adapter.py` ~l.838). Significant Coe edges therefore carry `adjusted_p_value = 0.1`. | Not applied, per the "fix only if it is the effective cutoff" rule |

Not applied / needs the curator: `pvalue_threshold` 0.1 → 0.01 on all 14 Coe analyses is recommended. It would make the placeholder `adjusted_p_value` on significant edges read 0.01 (the paper's cutoff) instead of 0.1. Significance calls would be unchanged. The 13 h timepoint could be confirmed against the GEO GSE264347 sample table online.
Graph-visible changes: none.
Snapshot entries affected: none.
CHANGELOG line: none (no change).

## Thompson 2016
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| `_dark` experiment name / background_factors | The uninfected experiment is mislabelled and has no light/dark factor | Name: yes. Both light-stress experiments had the identical name "MED4 Constant dark vs Constant light (90 µE m-2 s-1) (RNASEQ)". Neither said infected vs uninfected (paper p. 4: "comparing dark to light treatments with infected or uninfected constant"; spent-medium control bottles, p. 3). background_factors: **not confirmed**. Dark vs light is the *treatment* (`treatment_type: [light]`), so a light or darkness background factor would be wrong. The infected sibling likewise has none. | Names: `light_stress_constant_dark_med4_rnaseq_infected` → "MED4 Constant dark vs Constant light (90 µE m-2 s-1) in phage P-HM2 infected cultures (RNASEQ)". `light_stress_constant_dark_med4_rnaseq_dark` → "MED4 Constant dark vs Constant light (90 µE m-2 s-1) in uninfected cultures (spent medium, no phage) (RNASEQ)". Keys unchanged. background_factors unchanged. |

Not applied / needs the curator: background_factors left as-is (see above). The infected experiment has `[viral]` without `axenic`. The paper defers culture details to ref. [4], so this was not checked.
Graph-visible changes: Experiment `name` on the 2 light-stress experiments.
Snapshot entries affected: none.
CHANGELOG line: Thompson 2016: the two dark-vs-light experiments shared one name. Their names now say which is phage-infected and which is the uninfected spent-medium control.

## Anjur 2025
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| treatment_type | `[growth_phase]` is a poor fit for planktonic vs biofilm | Confirmed as a poor fit: both arms were grown "to mid-exponential phase". However, the closed `Experiment.treatment_type` vocabulary (nitrogen … growth_phase, compartment, chemical, oxygen, rna_decay, tss_mapping, genomic_analysis) has no lifestyle or attachment value. `growth_phase` is documented in the validator as "Growth state / multi-condition comparison", which is the nearest fit. | None |

Not applied / needs the curator: minting a new treatment_type (e.g. `lifestyle` or `attachment`) needs a vocabulary change in `config/controlled_vocabularies.yaml` plus the validator.
Graph-visible changes: none.
Snapshot entries affected: none.
CHANGELOG line: none.

## Bagby 2015
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| experimental_context timepoints | The text lists 3–24 h, but arrays were 0, 70 min, 6, 12, 24 h | Yes. Methods: "Samples were taken at t = 0, 70 min, 3 h, 6 h, 9 h, 12 h, 18 h and 24 h" and "hybridization were performed for two biological replicates at timepoints t = 0, 70 min, 6 h, 12 h and 24 h". | 3 gas-shock experiments: "Highest FC across timepoints (3h, 6h, 9h, 12h, 18h, 24h)." → "Highest FC across the hybridized timepoints (70 min, 6h, 12h, 24h; t = 0 is the reference). Samples were also taken at 3h, 9h and 18h but were not hybridized." |

Not applied / needs the curator: the analysis label `timepoint: max across 3-24h` (4 analyses) may also be wrong (the extremum can include 70 min). It was not on the list, so it was left alone.
Graph-visible changes: `experimental_context` on 3 Experiments.
Snapshot entries affected: none.
CHANGELOG line: Bagby 2015: the gas-shock experiments listed six array timepoints from 3 to 24 h. They now list the four that were hybridized (70 min, 6, 12, 24 h vs t = 0).

## Kratzl 2024
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| temperature | RNA "25-32C" vs proteomics "27.0-27.4C", both from EII | Yes, p. 4: "Samples for transcriptomics and proteomics were derived from EII and samples for metabolomics from EI at 27.0 °C or 27.4 °C, respectively". 27.4 °C belongs to EI (metabolomics, not in this config). 25–32 °C is the generic range in the Fig. 3 legend. | All 4 experiments: `temperature` → `27.0C` (was 25-32C ×2, 27.0-27.4C ×2) |
| light units | "m-2 s-2" should be "m-2 s-1" | Yes. The figure legends carry the s−2 typo, but Methods state "120 µmol photons s−1 m−2" (CellDEG section). | All 8 occurrences (4 `light_intensity` + 4 `experimental_context`) changed from "umol photons m-2 s-2" to "umol photons m-2 s-1" |

Not applied / needs the curator: none.
Graph-visible changes: Experiment `temperature`, `light_intensity` and `experimental_context` on 4 experiments.
Snapshot entries affected: none.
CHANGELOG line: Kratzl 2024: the omics experiments read 25–32 °C (RNA-seq) or 27.0–27.4 °C (proteomics) and "m-2 s-2". All four now read 27.0 °C (run EII) and µmol photons m-2 s-1.

## Tai 2009 (directory "Tal 2009")
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| publication.doi | Add the DOI if it can be confirmed | Yes. The main PDF prints "doi:10.1111/j.1462-2920.2009.01997.x", which matches the SI file name `emi_1997` and the README. | Added `doi: "10.1111/j.1462-2920.2009.01997.x"`. It is identical to the DOI already in `cache/pdf_extraction_cache.json`, so the Publication node id does not change. Directory not renamed. |

Not applied / needs the curator: none.
Graph-visible changes: none (same DOI as before).
Snapshot entries affected: none.
CHANGELOG line: none needed (no graph change).
