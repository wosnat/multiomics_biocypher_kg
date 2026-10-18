# Fix batch D

All six paperconfigs still validate (`VALIDATION PASSED`). Biller 2018's two warnings (`prefiltered` makes `pvalue_threshold` redundant) were there before this pass. `control_timing` was not touched.

## Biller 2016
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| MIT1002 growth_phase | `exponential` is wrong for the measured MIT1002 arm | Yes. p. 2833: "Heterotroph cell abundances declined for the first day following their introduction to the culture, were stable the following day, and then began to grow (Figure 1a)". p. 2837: "their cell numbers decreased for about 24 h ... they began to grow exponentially again after 48 h". p. 2839: "declined ... leveled off, and then began to increase after 2 days". So at 24 h MIT1002 is stable, and at 48 h it is only starting to grow. Exponential growth here is NATL2A, not MIT1002. No canonical value fits a stable, recovering population (no `lag` value exists), so I used `unknown` and not an `other:` slug. | `DE_alteromonas_MIT1002_coculture_24v12h` and `_48v12h`: `growth_phase: exponential` → `unknown`. A sentence with the page citation was appended to the `design_notes` of both MIT1002 experiments. |
| MIT1002 timepoint labels | "24h vs 12h" / "48h vs 12h" put the contrast in the label | Yes. The S3 header is "24 vs 12 hrs after addition". The 12 h reference is already stated in `control_condition` and `design_notes`. | `timepoint: 24h vs 12h` → `24h`; `48h vs 12h` → `48h`. `timepoint_hours` was already 24.0 / 48.0 and is unchanged. |

## Biller 2018
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| MIT1002 timepoint_hours on the hours-in-darkness clock | MIT1002 1 h / 5 h are the same samples as NATL2A 36 h / 40 h | Yes. The Table 1/2 footnotes say "first 5 h of extended darkness (36- and 40-h time points)". The S3 legend reads "between 36 and 48 h (1 to 13 h in extended darkness)". The S6B legend reads "following 1 to 5 h of extended darkness". MIT1002 and NATL2A were sequenced from the same coculture bottles (one RNA, reads split per genome, p. 15). | `DE_extended_darkness_vs_diel_MIT1002_1h`: `timepoint: 1h extended darkness` / `timepoint_hours: 1.0` → `1h extended darkness (36h)` / `36.0`. `_5h`: `5h extended darkness` / `5.0` → `5h extended darkness (40h)` / `40.0`. The labels now match the NATL2A analyses. |

## Weissberg 2025 (HOT1A3 day-11 control only)
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| HOT1A3 control_condition / name | HOT1A3 has no exponential phase in PRO99-lowN | Yes. p. 10: "There is no exponential growth stage. Instead, the cell numbers begin to decline as the experiment starts, likely due to C starvation" (Fig 6A/B). Results: "Alteromonas go into decline at the beginning of both axenic growth and co-culture". All 4 experiments use the same series sampled on day 11 as the reference (the filenames are `..._vs__day_11`). | 4 HOT1A3 growth_state experiments (RNASEQ/PROTEOMICS × axenic/coculture): `control_condition: PRO99-lowN exponential growth` → `PRO99-lowN day 11`. `name: HOT1A3 PRO99-lowN nutrient starvation vs PRO99-lowN exponential growth (...)` → `... vs PRO99-lowN day 11 (...)`. The MED4 experiments were not touched: their day 7 / day 11 references really are exponential. |
| growth_phase of `DE_hot1a3_rnaseq_coculture_vs_axenic_day_11` | `exponential` is wrong | Yes (same quote). The other HOT1A3 day-11 references are described as nutrient_limited in design_notes. | `growth_phase: exponential` → `nutrient_limited` |

## Lin 2015
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| DOI | README gives 10.1111/1462-2920.13104 | Yes. The PDF header reads "doi:10.1111/1462-2920.13104". It matches the DOI in `cache/pdf_extraction_cache.json`, so the Publication node id does not change. | Added `doi: "10.1111/1462-2920.13104"` under `papername` |
| Culture conditions | The paperconfig says 24C, 14:10 L:D, 40 umol, `diel`. The methods say 21 C, constant light, 30 umol | Partly verified (medium confidence). p. 456, Strains and growth conditions: "Axenic Prochlorococcus NATL2A cultures were maintained in the Pro99 medium at 21°C under constant cool white light (30 μmol quanta m–2 s–1)". The RNA-seq procedure gives no different growth conditions (cells were spun "at 21°C" before resuspension). Fig. S1/S2 (the same experiment, "redrawn from Zeng and Chisholm, 2012") mention no light:dark cycle. The 24C / 14:10 / 40 umol values have no source anywhere in the paper or its supplement. I could not read Zeng & Chisholm 2012 itself. | Both experiments: `temperature: 24C` → `21C`; `light_condition: 14:10 light:dark cycle` → `continuous light`; `light_intensity: 40 umol photons m-2 s-1` → `30 umol photons m-2 s-1`; `diel` → `light` in `background_factors` |
| Duplicate name | Both experiments share one name | Yes | Infected: `NATL2A P-limited vs P-replete control (RNASEQ)` → `NATL2A P-limited vs P-replete control, phage-infected (RNASEQ)`. The uninfected name is unchanged. |
| axenic background on infected | Both hosts are axenic | Yes. p. 456: "Axenic Prochlorococcus NATL2A cultures"; aliquots of the same cultures were infected (Fig. S1) | Infected `background_factors: [viral, diel]` → `[axenic, viral, light]`. Uninfected `[axenic, diel]` → `[axenic, light]`, from the diel fix above. |

## Lindell 2007
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| control_condition | The control is a paired uninfected culture amended with spent medium | Yes. Methods Summary (main PDF): "the paired control cultures were amended with filter-sterilized spent medium". Supp. Methods: "control treatments were amended with filter-sterilized spent medium". | `control_condition: Pro99 medium phage infection baseline` → `Uninfected control amended with filter-sterilized spent medium`. `name: MED4 Phage infected vs Pro99 medium phage infection baseline (MICROARRAY)` → `MED4 Phage infected vs uninfected spent-medium control (MICROARRAY)`, so the name matches the control. |
| temperature / light | 21 C, continuous white light at 10-25 umol | Yes. Supp. Methods: "at 21 °C under continuous white light at 10-25 µmol photon·m-1·s-1"; Methods Summary says the same. | `temperature: 24C` → `21C`; `light_intensity: 10 umol photons m-2 s-1` → `10-25 umol photons m-2 s-1`. `light_condition: continuous light` was already correct. |
| treatment_organism `Phage` → P-SSP7 | The phage is P-SSP7 | Yes, the paper says P-SSP7 | **Not applied.** `data/Prochlorococcus/treatment_organisms.csv` has no P-SSP7 entry (only `10239,Phage`), and P-SSP7 is not in `cyanobacteria_genomes.csv`. `treatment_organism: Phage` / `10239` left as is. |

## Martiny 2006
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|
| MIT9313 table_scope_detail | The Table 2 caption says q<0.05 at t=24 | Yes. The `01301table_2.csv` caption reads "MIT9313 genes identified as differentially expressed (q<0.05 at t=24)", and the columns stop at t=24. Table 1 (MED4) says t=48, so the MED4 detail is correct. | MIT9313: `Genes with q < 0.05 at t=48h, all timepoints reported` → `... at t=24h, all timepoints reported` |
| light_intensity per strain | Each experiment lists both strains | Yes. Methods, Culture Conditions: "12 µE m-2·s-1 ... for MIT9313 ... and 30 µE·m-2·s-1 for MED4" | MED4: `30 uE m-2 s-1 (MED4), 12 uE m-2 s-1 (MIT9313)` → `30 uE m-2 s-1`. MIT9313: same string → `12 uE m-2 s-1`. |

## Not applied / needs the curator
- **Lindell 2007 `treatment_organism`:** stays the generic `Phage` (taxid 10239). There is no canonical P-SSP7 organism in `treatment_organisms.csv` or the genome registry. Adding one would be a registry change outside this pass.
- **Lin 2015 culture conditions (applied, medium confidence):** the paper gives only the maintenance conditions (21 C, constant light, 30 umol). It says the RNA comes from Zeng & Chisholm 2012, which I could not read. The values I replaced had no source in this paper at all. If Zeng & Chisholm 2012 used different experimental light, it should be checked against that paper.
- **Biller 2016 `unknown`:** no canonical value describes a stable, recovering heterotroph population. The curator may prefer `other:lag` or `recovery`. The page citations are in `design_notes`.
- **Seen but not on my list (unchanged):**
  - Lin 2015: the infected analyses carry `growth_phase: nutrient_limited`, not `infected`.
  - Weissberg 2025: the axenic and coculture HOT1A3 growth_state experiments still share a `name` within each omics type (this was so before).

## Graph-visible changes
- Biller 2016, MIT1002 `Changes_expression_of` edges: `time_point` changes `24h vs 12h` → `24h` and `48h vs 12h` → `48h`. Analysis `growth_phase` changes `exponential` → `unknown`, which changes `Experiment.growth_phases` on the 2 MIT1002 experiments.
- Biller 2018, `darkness_extended_darkness_mit1002_rnaseq` edges: `time_point` gains ` (36h)` / ` (40h)`. `time_point_hours` changes 1.0 / 5.0 → 36.0 / 40.0, and with it the Experiment's `time_point_labels` / `time_point_hours` arrays.
- Weissberg 2025: on the 4 HOT1A3 growth_state Experiments, `control` and `name` change from "exponential growth" to "day 11". On `coculture_prochlorococcus_med4_hot1a3_rnaseq`, the growth phase of the day-11 edges changes `exponential` → `nutrient_limited`, so `Experiment.growth_phases` loses `exponential`.
- Lin 2015, both Experiments: `temperature`, `light_condition`, `light_intensity` and `background_factors` change (diel → light; the infected one gains axenic). The infected Experiment gets a new `name`. Publication `doi` is now from the paperconfig and has the same value as before.
- Lindell 2007: Experiment `name`, `control`, `temperature` and `light_intensity` change.
- Martiny 2006: `light_intensity` changes on both Experiments; MIT9313 `table_scope_detail` changes.

## Snapshot entries affected
None. `tests/kg_validity/snapshot_data.json` contains none of these experiment ids, the Biller 2016 / Lindell 2007 / Martiny 2006 / Lin 2015 / Weissberg 2025 DOIs, or the changed strings. The only Biller 2018 entries are the NATL2A DerivedMetric nodes (`s4a_natl2a_axenic`, `s5_natl2a_survival`), which this pass does not change.

## CHANGELOG line
Corrected metadata in six papers:
- Biller 2016: MIT1002 timepoints were labelled "24h vs 12h" / "48h vs 12h" and said to be in exponential growth; they are now "24h" / "48h" with growth phase unknown, because Alteromonas was stable, not growing.
- Biller 2018: MIT1002 extended-darkness points were at 1 h / 5 h on the darkness clock; they are now at 36 h / 40 h on the experiment clock, matching NATL2A.
- Weissberg 2025: the HOT1A3 day-11 reference was called "exponential growth"; it is now "day 11", nutrient_limited, since HOT1A3 declines from the start.
- Lin 2015: conditions were 24 C, 14:10 L:D, 40 umol (diel); they are now 21 C, continuous light, 30 umol. The infected experiment has its own name and an axenic background, and the DOI is recorded.
- Lindell 2007: the control is now the spent-medium uninfected culture, at 21 C and 10-25 umol.
- Martiny 2006: each strain now carries only its own light intensity, and MIT9313 significance is stated at t=24h.
