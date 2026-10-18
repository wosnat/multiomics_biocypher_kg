# Fix batch F

## Domínguez 2017
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| 1. `axenic` in background_factors | Paper may not say the SS120 cultures were axenic | The paper says nothing either way. Methods p. 12 ("Prochlorococcus strains and culture conditions") says only "Prochlorococcus marinus strain SS120 … was routinely cultured in polycarbonate flasks (Nalgene) using PCR-S11 medium". "axenic" does not appear in the main PDF, in the three supplementary PDFs or in legends.txt. Nothing says non-axenic. | **kept** `axenic`, as the rule requires when the paper is silent. Added `design_notes: Axenic status of the SS120 cultures is not stated in the paper.` |

Not applied / needs the curator: none.

## Christie-Oleza 2017
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| 2. growth_phase on `syn_wh7803_early_coculture_asw` (35 days, ASW) | `acclimated_steady_state` may be wrong; axenic ASW cultures die after 4–6 weeks | Yes. The treatment arm is the coculture with R. pomeroyi. p. 2: "axenic cultures died after 4–6 weeks". Fig. 1c legend, p. 2: "The arrow indicates the time point at which cells were collected for proteomic analysis (35 days)". In Fig. 1c, at the arrow both the coculture and the axenic curves sit at the maximal density plateau (~10^9 cells/ml). Exponential growth ended around day 20, and the axenic decline starts after the arrow. The coculture is at its post-exponential maximum, not in a long-term acclimated state. Its stable long-term state only follows later (weeks 10–40). | **fixed** `growth_phase: acclimated_steady_state` → `stationary`. Added a `design_notes` sentence on the 35-day harvest point (both arms at maximal density, before the axenic die-off). |

Not applied / needs the curator: none.

## Fuszard 2012
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| 3a. iTRAQ plex | The main text says 8-plex; the paperconfig and README say 4-plex | Yes, 8-plex. Main PDF Methods p. 6: "labelled with 8-plex iTRAQ reagents". Supplementary methods (MOESM3) do not name a plex but say "moieties 113,114 and 115 were added to the phosphate replete control replicates, and the moieties 116, 117 and 118" to the deplete ones. Tags 113 and 118 exist only in the 8-plex kit (the 4-plex kit is 114–117). This holds for all three strains, so there is no per-strain difference. | **fixed** on 3 experiments: `experimental_context` `iTRAQ 4-plex labeling` → `iTRAQ 8-plex labeling` (MIT9312: `(6 of 8 tags used: 113-115 replete, 116-118 deplete)`). Also the header comment `iTRAQ 4-plex` → `8-plex`, and README Topic `iTRAQ 4-plex` → `iTRAQ 8-plex`. |
| 3b. light regime | 13:11 L:D cycle → `diel` | Yes. Main PDF Methods p. 6: "in a 13:11 h light:dark regime at 23°C". Supp. methods: "grown under a 13:11 hour light:dark regime" | **fixed** on all 3 experiments: `background_factors: ["light"]` → `["diel"]` |

Not applied / needs the curator: `experimental_context` on all 3 experiments says "2D-LC MALDI-TOF/TOF MS/MS". The paper says SCX fractionation and then "a QStar XL Hybrid ESI Quadrupole time-of-flight" (main p. 6; supp. "ESI-qQ-TOF MS/MS"). It is not on this batch list, so I left it unchanged.

## Moreno 2023
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| 4a. tail and pairing of S3/S4 `t.test` p-values | S2 TTEST values are one-tailed (unpaired for MED4/SS120/BL107, paired for WH7803/WH8102; fixes_C). S3/S4 not yet checked | Recomputed with scipy from the normalized per-sample block (second set of 18 replicate columns) of all 10 S3/S4 tables, for every protein and all 4 t.test columns. I tried Student, Welch and paired tests, each one- and two-tailed. Two-tailed reproduces at most 0–10% of rows in any column, so every column is **one-tailed**. **Paired one-tailed** reproduces 87–100% of rows for Alteromonas in BL107/MED4/SS120/WH8102. It reproduces 53–60% for Alteromonas in WH7803 (only 15 proteins, p rounded to 2 decimals), and 80–100% once the tolerance allows that rounding. It also reproduces 82–99% for Marinobacter in BL107/MED4/WH7803/WH8102. **Marinobacter in SS120 is unpaired Student one-tailed**: 78–97% at rounding tolerance, against ≤11% for paired. The exception is **Marinobacter in WH8102 `t.test G mM/O`** (dark + 5 mM vs dark): no test or group pair reproduces it (≤5%). | **fixed** on all 60 experiments: the E3 `design_notes` sentence is extended. 43 get "The p-values are one-tailed (half the two-tailed value) from a paired t-test." (S2 WH7803/WH8102, all 20 S3, and S4 except SS120). 16 get "… from an unpaired Student t-test." (S2 MED4/SS120/BL107 and S4 SS120). 1 (`dark_high_glucose_marino_in_wh8102_proteomics`) gets "The test behind these p-values could not be reproduced from the per-sample columns; the other three comparisons in this table are one-tailed paired t-tests (half the two-tailed value)." |
| 4b. medium | Methods name PCR-S11 (Prochlorococcus) and ASW (Synechococcus), not Pro99/SN | Yes. Methods p. 22: "Prochlorococcus and marine Synechococcus cells were grown in PCR-S11 (60) and ASW media (61), respectively, in a culture room at 24°C." | **fixed** programmatically on all 60 experiments: `medium` `Pro99` → `PCR-S11` (24: MED4/SS120 and their Alteromonas/Marinobacter fractions) and `SN` → `ASW` (36: WH7803/WH8102/BL107 and their heterotroph fractions). The same wording is fixed inside `experimental_context` and the header comment. Checked by script: every experiment's `medium` equals its host culture's medium, and its `experimental_context` contains it. No `Pro99`/`SN` remains in the file. |

Not applied / needs the curator: none. One p-value column (Marinobacter in WH8102, dark + 5 mM) is unexplained, and its note says so.

## Bagby 2015
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| 5. timepoint label of the 4 gas-shock analyses | "max across 3-24h" is wrong; the maximum is over the hybridized 70 min–24 h | Yes. Methods p. 2234: "Samples were taken at t = 0, 70 min, 3 h, 6 h, 9 h, 12 h, 18 h and 24 h" … "hybridization … for two biological replicates at timepoints t = 0, 70 min, 6 h, 12 h and 24 h". Table S1 note 3 (moesm33): maximal fold change "with respect to t = 0" (air) / "the extremum of log2(fold change(t - t0), gas) - log2(fold change(t - t0), air)". So the extremum runs over the arrays at 70 min, 6, 12 and 24 h. 3 h was never hybridized. The experiments' own `experimental_context` already said so. | **fixed** on 4 analyses: `timepoint: max across 3-24h` → `max across 70min-24h`. `timepoint_hours: null` is kept (correct for an extremum). The same label at the start of the 4 `design_notes` and both README mentions are updated to match. |

Not applied / needs the curator: none.

## Graph-visible changes
- Christie-Oleza 2017 `syn_wh7803_early_coculture_asw`: edge `growth_phase` `acclimated_steady_state` → `stationary`, so `Experiment.growth_phases` for `early_coculture_syn_wh7803_proteomics` changes too. A new `design_notes` value is added.
- Fuszard 2012 (3 experiments): `background_factors` `[light]` → `[diel]` (also rolls up to Publication/OrganismTaxon `background_factors`), and `experimental_context` text changes.
- Moreno 2023 (60 experiments): `medium` and `experimental_context` change, and `design_notes` is extended.
- Bagby 2015 (4 analyses): edge `time_point` label → `max across 70min-24h` (and `Experiment.time_point_labels`). `design_notes` text changes.
- Domínguez 2017: new `design_notes` only.

## Snapshot entries affected
None. `tests/kg_validity/snapshot_data.json` has 0 hits for any experiment id, DOI or analysis id of these five papers.

## CHANGELOG line
Fuszard 2012 now records the 13:11 light:dark regime (`background_factors` `diel`, was `light`) and 8-plex iTRAQ (was 4-plex). The Moreno 2023 cultures are in PCR-S11 / ASW (were Pro99 / SN), and their notes now say the reported p-values are one-tailed (paired or unpaired, per table). The Bagby 2015 gas-shock fold changes are labelled `max across 70min-24h` (were `3-24h`, though 3 h was never hybridized). Christie-Oleza 2017's 35-day ASW coculture proteome is `stationary` (was `acclimated_steady_state`).
