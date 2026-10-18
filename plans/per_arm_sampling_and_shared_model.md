# Control timing + shared-model identifier on Experiments

**Status:** DONE — implemented 2026-10-15, revised and curation closed 2026-10-18 (on `main`, rebuild
pending). See **§ Revision 2026-10-18** below; the design text after it is the 2026-10-15 version,
kept for history.
**Dialog:** `multiomics_explorer/docs/kg-specs/2026-10-13-per-arm-sampling-shared-model-explorer-asks.md`
**Origin:** `plans/backlog.md` → "Analysis-driven asks — pro_n_import_coculture", P2 item 3.

## Revision 2026-10-18 (curator review, `plans/control_timing_curation.md`)

- **`statistical_model_id` dropped** — field, adapter emission, schema, validator rule (old rule 4)
  and tests. Most papers do not report how DE was fitted, so the id was mostly a guess. Model / fit
  sentences were then removed from every `design_notes` too (Biller 2016/2018, Hackl 2023,
  Weissberg 2025, Barreto 2022), keeping what each table compares.
- **`time_matched` renamed `matched`** so the value does not imply the same clock time: it covers
  arms matched by sampling time **or** by growth state (R1). `design_notes` on growth-state-matched
  experiments says how (added: Beliaev 2014 ×6, Fuszard 2012 ×3, Barreto 2022 cyanobacteria ×3).
- **Weissberg 2025 `coculture_alteromonas_hot1a3_med4_rnaseq` → `matched`** (D0, override; the
  axenic repeat was run to match growth dynamics), with a note naming the day 7 / day 14 controls.
  D1 accepted (`matched`), D2 confirmed (DEqMS; per the author, the four proteomics experiments are
  one model — notes updated), D3 confirmed (`death`), D4 accepted (Hennon: two fits).
- Validator rules are now: 1 required + in vocab · 2 `baseline`/`mixed` ⇒ `design_notes` ·
  3 `mixed` ⇒ ≥ 2 timepoints · 4 no `|`/`'` in notes · 5 warning on a non-DE experiment.
- Counts: 174 DE experiments — 150 `matched`, 21 `baseline`, 3 `mixed`.
- **Explorer contract change:** the accepted contract (T1–T5, SIM-001..004) named `time_matched`
  and `statistical_model_id`; both changed — the explorer dialog needs a note before release.

## Problem

1. **A contrast shows one sampling time, even when its arms were sampled at different times.**
   `timepoint` / `growth_phase` describe the treatment arm (by an unstated convention). Weissberg
   2025 `DE_med4_rnaseq_coculture_vs_axenic_day_18` compares coculture day 18 with axenic day 14,
   which was in death phase; the graph says "day 18".
2. **Nothing says which experiments are one fitted model.** Weissberg's MED4 RNA contrasts come from
   one DESeq2 model (one is an exact linear combination of others), so they are not independent
   replications.

## Design — three curated Experiment fields, nothing computed post-import

| paperconfig | Experiment property | Type | Rule |
|---|---|---|---|
| `publication.control_timing` (paper default) + experiment `control_timing` (override) | `control_timing` | str, closed `{time_matched, baseline, mixed}` | **Required** on every experiment with ≥ 1 `statistical_analyses` entry; absent on experiments with no DE contrast (like `table_scope`). Every value is a curated claim — the adapter never defaults it (explorer T1). |
| experiment `statistical_model_id` | `statistical_model_id` | str, sparse, `"{doi}:{key}"` (adapter adds the DOI prefix, T3) | Only on experiments that are contrasts from one fitted model shared with ≥ 1 other experiment. |
| experiment `design_notes` | `design_notes` | str, sparse, free text | Required when `control_timing ∈ {baseline, mixed}`. The only per-timepoint record of the control arm. Added to `experimentFullText`. |

Linked experiments are found by matching equal `statistical_model_id` (no stored list or count —
R3; explorer SIM-002). `time_point` / `growth_phase` keep names, types and values; their docs say
**treatment arm**.

`design_notes` writing rules (T4; in the `/paperconfig` skill): exact `timepoint` labels; other
experiments by paperconfig key verbatim; control sampling time + growth phase in the `growth_phase`
vocabulary; for a shared model, the model and any exact linear relation; no `|` / `'`.

### Validator (`scripts/validate_paperconfig.py`) — errors unless noted

1. Every experiment with `statistical_analyses` resolves a `control_timing` (override or paper
   default) in the vocab (values loaded from `controlled_vocabularies.yaml`).
2. `baseline` / `mixed` ⇒ `design_notes` present.
3. `mixed` ⇒ the experiment's analyses span ≥ 2 distinct `timepoint` labels (T5).
4. Each `statistical_model_id` value is on ≥ 2 experiments of the paper, all with the same
   `organism` and `omics_type`.
5. `design_notes` contains no `|` or `'`.
6. Warning: `control_timing` / `statistical_model_id` on an experiment with no `statistical_analyses`.

Because `tests/test_paperconfig_validation.py` validates every wired paperconfig, rule 1 can only
land **in the same commit as the curation pass** (otherwise all 49 papers fail).

### Code

- `omics_adapter`: emit the three properties on Experiment nodes (`_clean_str` on `design_notes`);
  `control_timing` only for experiments with analyses; prefix `statistical_model_id` with the DOI.
- `config/schema_config.yaml`: three Experiment properties.
- `scripts/post-import.sh` + `.cypher`: drop + recreate `experimentFullText` with `design_notes`
  (explorer: expect `search_text` goldens to move).
- Vocab: `Experiment.control_timing` (sparse, closed, value descriptions) — text in § Vocabulary
  and documentation.

## Curation pass (required — the fields mean nothing until a curator fills them)

Scope (paperconfig counts, 2026-10-15): **49 papers, 209 experiments, 174 with DE contrasts.** Every
paper with DE contrasts gets `publication.control_timing`; that is the whole point of making it
required.

**C1 — per paper, read the methods** (sampling schedule; which control was sampled when; how DE was
fitted — one model over all samples or one fit per comparison). The `control_condition` text is
not enough: `"Axenic"` hid the Weissberg day-14 control.

**C2 — write the paperconfig**: paper default = the value most of its DE experiments take; override
the rest; `design_notes` on every `baseline` / `mixed` experiment; `statistical_model_id` on shared
fits. Methods that do not say when the control was sampled: pick the value the design implies and
say so in `design_notes` — do not guess silently.

**C3 — curation log** `plans/control_timing_curation.md`: one row per paper — default, overrides,
model ids, and the evidence (methods section / page, one short quote). Reviewed by the user via the
log + `git diff` before the validator commit.

**C4 — validate**: `validate_paperconfig.py` over all papers with the new rules on; fix and re-run.

Papers can be curated in batches by parallel subagents (each returns its log rows + diff); the log
is the review surface.

Known cases going in (still checked in C1):

| Paper | `control_timing` | `statistical_model_id` |
|---|---|---|
| Weissberg 2025 | 8 `growth_state_*` → `baseline` (vs exponential day 7 / day 11); `coculture_alteromonas_hot1a3_med4_rnaseq` → `mixed` (d11 matched; d18 vs axenic d14, death); `coculture_prochlorococcus_med4_hot1a3_rnaseq` (d11/d18/d31) → **curator** | one per organism × omics: MED4 RNA (3 experiments), HOT1A3 RNA (3); proteomics → **curator** |
| Barreto 2022 | expected `time_matched` | one id over the 8 EZ55 experiments (`~0+cyano+treatment+cyano:treatment`); the 3 cyanobacteria-side experiments are separate fits |
| Hennon 2017 | expected `time_matched` | the 2 EZ55 experiments if one GLM → **curator** |
| Kaur 2018 | 4 → `baseline` (vs day 1, which is not itself reported) | — |
| Bagby & Chisholm 2015 | `carbon_air_0036_co2_21_med4_microarray` → `baseline` (vs time 0 pre-shock); other 3 → `time_matched` | — |
| Biller 2018, Coe 2024, Thompson 2016, Tetu 2019 | C1 | probably separate fits — C1 |

## Not changing

- `timepoint`, `timepoint_hours`, `growth_phase`, `time_point_*` arrays, `growth_phases`: names,
  types, values. Only their documentation sharpens (treatment arm).
- No edge properties, no Experiment arrays, no post-import rollups, no Arm/Sample node.
  Withdrawn from the 2026-10-14 contract: `control_time_point` / `_hours` / `control_growth_phase`
  edge props, `time_point_control_*` arrays, `control_growth_phases`, `linked_experiment_ids` /
  `_count`, post-import singleton clearing, the GLM-text heuristic. Structured per-timepoint fields
  can return later as an additive change if `design_notes` proves insufficient.
- No explorer contract break; no `mcp_min_version` bump.

## Explorer dialog (2026-10-15)

SIM-001 accepted; SIM-002 no stored count (explorer computes linked ids in verbose); SIM-003
`design_notes` sparse-compact in `list_experiments`; SIM-004 name `design_notes`. Explorer surfaces
`control_timing` compact + filter + `by_control_timing` breakdown, and on DE rows when not
`time_matched`; `design_notes` in the DE per-experiment envelope. KG adopted T1 (required, paper
default + override), T2 (graph test below), T3 (DOI prefix), T4 (writing rules), T5 (`mixed` ⇒
≥ 2 timepoints; counts are graph counts — 209 Experiment nodes).

## Decisions needed (curator)

1. Weissberg: HOT1A3-side `coculture_*` control days; whether the proteomics layers are separate fits
   from the RNA models; whether the `nutrient_limited` label on `axenic_day_14_vs_exp_day_7` should
   read `death`.
2. Hennon 2017: one EZ55 GLM or two fits.

## Acceptance criteria

- Unit: adapter emits `control_timing` only on experiments with analyses, `statistical_model_id`
  DOI-prefixed, `design_notes` cleaned; validator rules 1–6 each covered by a test;
  `test_paperconfig_validation.py` green over all 49 papers.
- Curation log covers every paper with DE contrasts.
- KG validity: `control_timing` present on every Experiment with `Changes_expression_of` edges and
  in the vocab; `mixed` experiments have `time_point_count ≥ 2`; `baseline` / `mixed` carry
  `design_notes`; every `statistical_model_id` is shared by ≥ 2 Experiment nodes with the same
  `organism_name` + `omics_type` (T2); Weissberg `coculture_alteromonas_hot1a3_med4_rnaseq` reads
  `mixed`, Kaur's four read `baseline`, the 3 MED4 RNA experiments share one id, Barreto's 8 EZ55
  experiments share one id.
- `experimentFullText` includes `design_notes`.
- `pytest -m "not slow and not kg"` + `pytest -m kg` green; CHANGELOG `### Added` entry; CLAUDE.md
  Experiment key-facts bullet updated on implementation.

## Vocabulary and documentation

Applied (2026-10-15, revised 2026-10-18): `config/controlled_vocabularies.yaml`
(`Experiment.control_timing`), `.claude/skills/paperconfig/SKILL.md` ("Control timing and design
notes"), `CLAUDE.md` (template + Experiment key-facts bullet), `CHANGELOG.md`.
