# Control timing + shared-model identifier on Experiments

**Status:** simplified 2026-10-15 (supersedes the 2026-10-13/14 per-arm edge-property design);
explorer accepted the simplified contract 2026-10-15 with T1–T5, all adopted (see § Explorer dialog).
Not started: curation pass, validator, adapter, vocabulary and docs (the exact vocab / skill /
CLAUDE.md text is drafted in § Vocabulary and documentation, to apply with the implementation).
**Dialog:** `multiomics_explorer/docs/kg-specs/2026-10-13-per-arm-sampling-shared-model-explorer-asks.md`
**Origin:** `plans/backlog.md` → "Analysis-driven asks — pro_n_import_coculture", P2 item 3.

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

## Vocabulary and documentation (drafted 2026-10-15, apply with the implementation)

### `config/controlled_vocabularies.yaml` — after `Experiment.table_scope`

```yaml
Experiment.control_timing:
  # provenance: curated per paperconfig (plans/per_arm_sampling_and_shared_model.md, 2026-10-15).
  #   Required on every experiment that has >=1 statistical_analyses entry: a paper-level default
  #   (publication.control_timing) that each experiment may override, so every value is a
  #   curator's claim, never an adapter default (explorer T1). Absent on experiments with no DE
  #   contrast (metabolomics / derived-metric-only), like table_scope. baseline and mixed require
  #   design_notes; mixed requires >=2 timepoints (validator).
  applies_to: Experiment
  applies_to_kind: node
  property: control_timing
  value_type: string
  closed: true
  sparse: true
  values: [time_matched, baseline, mixed]
  description: >
    Whether the control arm of each contrast was sampled at the same time as
    the treatment arm. time_point and growth_phase always describe the
    treatment arm, so on a baseline or mixed experiment they say nothing about
    when or in what state the control was sampled; design_notes does. Absent
    on experiments with no differential-expression contrast.
  value_descriptions:
    time_matched: every contrast compares treatment and control arms sampled at the same time
    baseline: every contrast compares against a control sampled at another time, typically a single early reference point; read design_notes
    mixed: some timepoints are time-matched and others compare against a control sampled at another time; read design_notes for which
```

(Checked 2026-10-15: passes `tests/test_controlled_vocab.py` — no apostrophe / pipe / build
provenance in descriptions. Adding it makes `pytest -m kg::test_every_vocabulary_has_a_node` fail
until the next rebuild, as with any vocab addition.)

### `.claude/skills/paperconfig/SKILL.md`

1. § 1 Publication Metadata — add to the YAML block, after `papermainpdf`:
   `control_timing: time_matched  # paper-level default; required when the paper has DE contrasts (see "Control timing" below)`
   and a bullet: "`control_timing`: *(required when any experiment has `statistical_analyses`)*
   Paper-level default for `Experiment.control_timing`; an experiment may override it. See
   **Control timing, shared models and design notes** below."
2. § 2 Experiments Block — in the example, after `table_scope_detail`:
   ```yaml
      # control_timing: mixed           # overrides publication.control_timing (time_matched | baseline | mixed)
      # statistical_model_id: med4_rnaseq_deseq2   # only when >=2 experiments come from one fitted model
      # design_notes: >-                # required when control_timing is baseline or mixed
      #   Day 18: coculture sampled on day 18 vs axenic sampled on day 14 ...
   ```
3. New subsection before "Canonical Vocabulary (`treatment_type` and `background_factors`)":

   > #### Control timing, shared models and design notes
   >
   > `timepoint` and `growth_phase` on an analysis describe the **treatment arm**. Three experiment
   > fields record what they cannot: when the control arm was sampled, and whether several
   > experiments are contrasts from one fitted model. **Read the paper's methods for both** — the
   > `control_condition` text is not enough (Weissberg 2025's `"Axenic"` hid a day-18 vs day-14
   > comparison).
   >
   > | Field | Where | Rule |
   > |---|---|---|
   > | `control_timing` | `publication` (default) + experiment (override) | **Required** for every experiment with `statistical_analyses`, via the paper default or an override. Absent on experiments with no DE contrast. |
   > | `statistical_model_id` | experiment | Only when ≥ 2 experiments of this paper are contrasts from **one fitted model** (one DESeq2 / edgeR / limma fit over shared samples). Free key, unique within the paper (`med4_rnaseq_deseq2`); the KG stores it as `{doi}:{key}`. All experiments sharing it must have the same `organism` and `omics_type`. Separate fits get no id. |
   > | `design_notes` | experiment | Free text. **Required** when `control_timing` is `baseline` or `mixed`; recommended with `statistical_model_id`. |
   >
   > | `control_timing` | When |
   > |---|---|
   > | `time_matched` | Every contrast compares arms sampled at the same time (treated vs untreated at 1 h, coculture vs axenic both on day 11). Most papers. |
   > | `baseline` | Every contrast compares against a control sampled at another time — usually one early reference ("vs exponential day 7", "vs time 0 pre-shock", "vs day 1"). |
   > | `mixed` | Some timepoints are time-matched, others are not. Needs ≥ 2 timepoints. |
   >
   > Paper default = the value most of the paper's DE experiments take; override the rest.
   > `time_matched` is a claim checked in the methods, not a fallback — if the methods do not say
   > when the control was sampled, pick the value the design implies and say so in `design_notes`.
   >
   > `design_notes` rules: exact `timepoint` labels (`"day 18"`); other experiments by paperconfig
   > key verbatim (graph id `{doi}_{key}`); the control sampling time and its growth phase in the
   > `growth_phase` vocabulary; for a shared model, what the model is and any exact linear
   > relation between contrasts; no `|` or `'`.
   >
   > ```yaml
   >     coculture_alteromonas_hot1a3_med4_rnaseq:
   >       control_timing: mixed
   >       statistical_model_id: med4_rnaseq_deseq2
   >       design_notes: >-
   >         day 11: coculture vs axenic, both sampled on day 11. day 18: coculture sampled on day 18 vs
   >         axenic sampled on day 14, when the axenic culture was in death phase - not time-matched.
   >         One DESeq2 model with growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic and
   >         growth_state_pro99lown_nutrient_starvation_med4_rnaseq_coculture; the day 18 contrast is an
   >         exact linear combination of contrasts in those two experiments.
   > ```
4. Validation Checklist — add:
   - `publication.control_timing` is set (paper with DE contrasts), checked against the methods; experiments that differ override it
   - Every `baseline` / `mixed` experiment has `design_notes` naming exact `timepoint` labels and the control arm sampling time + growth phase; `mixed` has ≥ 2 timepoints
   - `statistical_model_id` is set on every experiment from a shared fitted model (≥ 2 per id, same `organism` + `omics_type`), and on no separately fitted experiment
5. Validation Script list — add validator rules 1–6 once implemented.

### `CLAUDE.md`

- "Adding Omics Data from Publications" template: after `papermainpdf` add
  `control_timing: time_matched  # REQUIRED with DE contrasts: time_matched | baseline | mixed — paper default, checked in the methods`;
  after `light_intensity` add commented `control_timing` override, `design_notes`, and
  `statistical_model_id` lines pointing to the skill section.
- Key graph facts, Experiment bullet: `control_timing` / `statistical_model_id` / `design_notes`
  (after the rebuild).
