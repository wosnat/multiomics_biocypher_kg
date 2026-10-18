# Control-timing curation log (C1–C3)

**Status:** first pass done 2026-10-15 by 8 parallel agents; **awaiting human review** before the
validator commit (plan `plans/per_arm_sampling_and_shared_model.md`). Nothing committed.
Agent instructions: `plans/control_timing_curation/INSTRUCTIONS.md`. Per-paper evidence tables
(experiment · value · confidence · quote + page) are in `plans/control_timing_curation/batch_{A..H}.md`.

## Review decisions (2026-10-18, partial)

Applied: **R1** keep growth-state matching, renamed `time_matched` → **`matched`**, matching notes
added (Beliaev ×6, Fuszard ×3, Barreto cyanobacteria ×3). **D0** Weissberg MED4 coculture vs axenic →
`matched` with note. **D1** accepted. **D2** confirmed (author: the four proteomics experiments are
one model; notes updated). **D3** confirmed. **D4** accepted (two fits). **D5 / D5b / D7 / D8 / D9
moot**: `statistical_model_id` dropped entirely (mostly unreported, so a guess); model / fit sentences
removed from every `design_notes` (2026-10-18), keeping what each table compares. **D6** (round 2): Thompson 2011
default flipped to `matched`; the 48h MED4 and both rescue experiments carry explicit `baseline`
overrides with notes. **Round 3 (opens page):** spot checks S1, S2, S4–S10 accepted as
curated (S2 Hackl: the notes already say the control harvest time is unstated); **S3** Tetu PVC →
`matched`. Data errors E1, E2, E4–E7 → fix after re-checking the paper; **E3** raw p-values → keep, with
a note on each affected experiment; the 21 smaller fixes → re-check and apply. Fix pass run by 5
agents per `plans/control_timing_curation/FIX_INSTRUCTIONS.md`; results in
`plans/control_timing_curation/fixes_{A..E}.md`. Counts now: 150 `matched`, 21 `baseline`, 3 `mixed`. In the text
below, `time_matched` reads `matched` and the model-id columns are historical.

## Totals

36 papers with DE contrasts (both `paperconfig_files.txt` lists) · 174 DE experiments ·
`time_matched` 149 · `baseline` 22 · `mixed` 3 · 6 `statistical_model_id` values over 21 experiments.
Mechanical check (2026-10-15): every DE experiment resolves a value; every `baseline`/`mixed` has
`design_notes`; no `|`/`'` in notes; every model id is on ≥ 2 experiments; no field on a non-DE
experiment. `validate_paperconfig.py` unchanged on every edited file.

## 1. Rule question — decide before review (affects the vocabulary text)

Agents met arms that were **harvested in the same growth state but not at the same time** and
classed them `time_matched` (medium confidence). The vocab text says "sampled at the same time".

| Case | Batch |
|---|---|
| Beliaev 2014 — sequential steady states in one chemostat (axenic → lactate coculture → HCO3- coculture) | H |
| Al-Hosani 2015 — both exponential, control days 2/4 vs salt days 4/6 | D |
| Fuszard 2012 — both at stationary onset (NATL2A mid-exponential), no common date | G |

Options: (a) `time_matched` means "same time **or** same physiological state (steady state /
matched growth phase)" — reword the vocab + skill text; (b) strict — these become `baseline`
(misleading: nothing is compared against an earlier reference). Recommendation: (a), with a
`design_notes` line on each such experiment saying how the arms were matched.

## 2. Decisions (step 0 of the plan, now with evidence)

| # | Decision | Finding | Recommendation |
|---|---|---|---|
| D0 | Weissberg `coculture_alteromonas_hot1a3_med4_rnaseq` | "day 11" = coculture d11 vs axenic **d7** (axenic MED4 from a separate repeat experiment); "day 18" = coculture d18 vs axenic d14. Source: author repo `~/github/CC1A3` sample metadata + export notebooks | `baseline` (curated), **not** `mixed`; fix the plan acceptance criterion |
| D1 | Weissberg HOT1A3-side `coculture_prochlorococcus_med4_hot1a3_rnaseq` | axenic HOT1A3 A1/A2/A3 = days 11/18/31, same run and batch; contrasts C1vsA1, C2vsA2, C3vsA3 | `time_matched` (high) |
| D2 | Weissberg proteomics model | DEqMS (limma), axenic and coculture fitted separately per organism; long-term pooled contrast another fit | **APPLIED 2026-10-18** — no `statistical_model_id`; `test_type` should be DEqMS, not DESeq2 |
| D3 | Weissberg axenic MED4 day-14 `growth_phase` | paper: "dying cultures", Fig 3C/D "axenic death stage" | **APPLIED 2026-10-18** — `death` (RNA + proteome) |
| D4 | Hennon 2017 EZ55: one GLM or two | edgeR `logCPM` differs between the two tables for the same gene (−2.5 to +3.6 on 48 shared genes) — impossible in one fit | two fits, no id (high) |
| D5 | Biller 2018 NATL2A | methods: one DESeq2 dataset per organism, pairwise Wald tests (no code) | `natl2a_rnaseq_deseq2` (medium) — contradicts plan prior |
| D6 | Thompson 2011 paper default | 3 of 5 not time-matched → default `baseline`, 2 overrides | accept majority rule, or flip to `time_matched` + 3 overrides (cosmetic) |
| D7 | Pandhal 2007 | HL/LL and ML/LL share the LL iTRAQ channel in the same runs; per-pair ratios | no id (shared reference channel ≠ fitted model) |

Shared models verified numerically: Weissberg MED4 `C2vsP2 = C2vsC1 + C1vsP1 − P2vsP1` to 1e-15
on 1,855 genes (HOT1A3: 4,012/4,015, 3 DESeq2-zeroed); Hackl 2023 identical `baseMean` across 3 GEO
tables and MMC−UV = (MMC−ctrl) − (UV−ctrl) exactly.

## 3. Data bugs found — the live graph is wrong today (outside this curation's scope)

| Sev | Paper | Problem | Fix |
|---|---|---|---|
| **high** | Aharonovich 2016 | 2 HOT1A3 tables (`fold change 9313/9313DIL`, `fold change 9313+1A3/med4+1A3`) are linear ratios (0 negatives, max 282,301 / 4,801; 422 / 874 values < 1) read as log2 — down-regulated genes appear up | `fold_change_type: linear` on both analyses (verified 2026-10-15) |
| **high** | Weissberg 2025 | `DE_med4_rnaseq_axenic_longterm_vs_exp_day_7` and `DE_hot1a3_rnaseq_axenic_longterm_vs_day_11` are **coculture** d60+d89 vs axenic — no axenic RNA exists at those days | **FIXED 2026-10-18** — author re-export (tables dropped / re-exported, paperconfig updated) |
| **high** | Weissberg 2025 | `DE_hot1a3_proteome_axenic_day_31_vs_day_11` holds **day 89 vs day 11** (export notebook overwrote the file); real d31 contrast never exported | **FIXED 2026-10-18** — author re-export (tables dropped / re-exported, paperconfig updated) |
| **high** | Barreto 2022 | 6 `coculture_*_at_{400,800}_ez55_rnaseq` tables are partner deviation from the partner-averaged coculture response, not "vs axenic" (3 partner tables sum to 0 per gene) | rename experiments + `control_condition`; S2 = axenic-only pCO2 effect; S4 ambiguous (author) |
| med | Moreno 2023, Barreto 2022, Hennon 2017, Al-Hosani 2015 | raw p-values mapped to `adjusted_p_value_col` | decide policy (null adjusted p, or keep with note) |
| med | Fang 2019 | `treatment_organism: Phage` on a cell-free lysate → spurious coculture edge; `axenic` missing | remove partner, add `axenic` |
| med | Domínguez 2017 | proteomics at 8 h (cfg 24h); azaserine 100 uM (cfg 10); 4 uE blue light (cfg 25) | fix fields |
| med | Oleza 2017 | `EXOPROTEOMICS` / `exoproteome`, but methods describe cell-pellet proteomes; "long-term" = 10 days | fix omics_type/compartment/timepoint |
| med | Fuszard 2012 | non-axenic per supp methods, cfg says axenic; NATL2A/SS120 light swapped | fix fields |

## 4. Smaller paperconfig findings (not fixed)

- **growth_phase**:
  - Moreno light analyses read `stationary` but look exponential; dark timepoint_hours null (24).
  - Kaur 2018 and Ma 2022 use `acclimated_steady_state` on batch/exponential days.
  - Tolonen N-deprivation 24/48 h are `acute_stress`.
  - Biller 2016 MIT1002 `exponential` is wrong for declining Alteromonas.
  - Weissberg HOT1A3 day-11 is labelled "exponential" (no exponential phase in this medium).
- **culture conditions**:
  - Moreno light (4 / 40 uE, cfg 15).
  - Barreto (ASW + SEv/PEv, 21 C, 150 uE, 14:10; cfg Pro99, 22 C, continuous, 30).
  - Lin 2015 (21 C, constant light; cfg diel).
  - Lindell 2007 (control text, temp, light; generic Phage).
  - Tolonen clusters (temp/light).
  - Kratzl (temperature mismatch, units typo).
  - Martiny (both strains' light on each experiment).
- **labels / scope**:
  - Thompson 2011 rescue condition names.
  - Capovilla "day 1 and day 3" is day 1 only.
  - Martiny MIT9313 `table_scope_detail` t=48h vs caption t=24.
  - Aharonovich HOT1A3 `control_condition`, MED4 columns swapped, threshold differs.
  - Anjur `growth_phase` treatment type.
  - Bagby `experimental_context` timepoints.
  - Thompson 2016 uninfected `_dark` key / background_factors.
  - Biller 2018 MIT1002 vs NATL2A clocks.
  - Coe 13h missing from methods, threshold 0.1 vs 0.01.
  - Tetu HDPE 120 vs 180 min (GEO).
- **metadata**:
  - Lin 2015 and Tal 2009 (author is Tai) have no DOI.
  - Lin 2015 both experiments share one `name`.
- **cross-paper**: Barreto 2022 re-analyses Hennon 2017 libraries — not independent replications;
  a paper-scoped `statistical_model_id` cannot say so (candidate future field).

## 5. Low / medium-confidence values to spot-check

| Paper | Experiment(s) | Value | Why uncertain |
|---|---|---|---|
| Ma 2022 | `coculture_ecoli_vs_axenic_utex2973_rnaseq` | time_matched | RNA sampling day never stated (low) |
| Hackl 2023 | all 3 | time_matched | control harvest time vs 2 h MMC / 1 h UV not stated (low) |
| Tetu 2019 | 2 PVC experiments | baseline | only GEO says PVC 90 min vs control 180 min |
| Lin 2015 | `phosphorus_plimited_natl2a_rnaseq_uninfected` | mixed | no replete library at 50 h; 46 h or 59 h control? |
| Doron 2016 | all 3 | time_matched | 2–5 min control may predate infection (→ mixed) |
| Huang 2020 | infection | baseline | single 30-min uninfected control for all times |
| Tolonen 2006 | 4 `growth_medium_growth_on_*` | time_matched | control harvest time not stated |
| Steglich 2006, Anjur 2025, He 2022 | all | time_matched | matched harvest implied, day not stated |
| Biller 2016 | `mit1002_rnaseq_deseq2` | model id | no code; shared 12 h reference |
| Johnson 2026b | both | no id | one fit vs per-contrast undeterminable |
| Tetu 2019 | per strain HDPE+PVC | no id | possibly one `~condition` fit; paper silent |

## Per-paper result

| Paper | Batch | Default | Overrides | Model ids |
|---|---|---|---|---|
| Weissberg 2025 | A | baseline | HOT1A3-side coculture → time_matched | med4_rnaseq_deseq2 (3), hot1a3_rnaseq_deseq2 (3) |
| Moreno 2023 | B | time_matched | — | — |
| Barreto 2022 | C | time_matched | — | ez55_rnaseq_edger_glm (8) |
| Hennon 2017 | C | time_matched | — | — |
| Biller 2016 | C | baseline | NATL2A → time_matched | mit1002_rnaseq_deseq2 (2) |
| Biller 2018 | C | time_matched | — | natl2a_rnaseq_deseq2 (2) |
| Coe 2024 | C | time_matched | — | — |
| Aharonovich 2016, Al-Hosani 2015, Anjur 2025, Fang 2019, He 2022 | D | time_matched | — | — |
| Bagby 2015 | D | time_matched | air experiment → baseline | — |
| Lin 2015 | E | time_matched | uninfected NATL2A → mixed | — |
| Lindell 2007 | E | mixed | — | — |
| Martiny 2006, Steglich 2006 | E | time_matched | — | — |
| Tetu 2019 | E | time_matched | 2 PVC → baseline | — |
| Thompson 2011 | F | baseline | 2 starvation → time_matched | — |
| Tolonen 2006 | F | time_matched | MIT9313 N-deprivation → mixed | — |
| Read 2017, Thompson 2016, Capovilla 2023 | F | time_matched | — | — |
| Domínguez 2017, Fuszard 2012, Pandhal 2007, Johnson 2026b, Tal 2009 | G | time_matched | — | — |
| Hackl 2023 | G | time_matched | — | mit0604_rnaseq_deseq2 (3) |
| Beliaev 2014, Kratzl 2024, Ma 2022, Oleza 2017, Doron 2016 | H | time_matched | — | — |
| Kaur 2018, Huang 2020 | H | baseline | — | — |
