# Timepoint & Growth-Phase Backfill — Design

**Date:** 2026-04-12
**Status:** Design approved, awaiting implementation plan
**Scope:** Populate timepoint and growth-phase metadata across all 30 paperconfigs (23 Prochlorococcus + 7 Synechococcus) — both time-series and non-time-series experiments.

## Problem

Every publication reports *when* its cells were sampled and *what physiological state* they were in, but this information is only partially structured in our paperconfigs today:

- `timepoint` / `timepoint_hours` exist but are set to `null` for non-time-series analyses, even though the paper clearly states a sampling time (e.g., Tetu 2019 exposes cells to HDPE leachate for 120 min).
- Physiological state at sampling (log phase, stationary, nutrient-limited, acclimated steady-state) has no structured field at all. It lives in free-text `experimental_context` when present.

This blocks biologically meaningful queries like "all acute (<6 h) nutrient-stress responses in exponential-phase cells" or "compare stationary vs exponential responses to the same treatment."

## Goals

1. Every `statistical_analyses[]` entry in every paperconfig has populated `timepoint`, `timepoint_hours`, and a new `growth_phase` field.
2. The field is propagated through the pipeline: paperconfig → adapter → `Changes_expression_of` edge → aggregated onto the `Experiment` node.
3. Extraction is LLM-assisted with human review (per-paper git-diff review), not manual re-reading of 30 PDFs.

## Non-goals

- No retroactive inference of growth phase from expression data itself.
- No changes to `is_time_course` logic (stays derived from `time_point_count > 1`).
- No stage-3 validation loop like cluster extraction — data is too shallow (4 fields per row) to justify that overhead; direct git-diff review is sufficient.
- No re-extraction of cluster descriptions or other paperconfig fields.

## Design

### 1. Schema additions

**Paperconfig — analysis-level fields (per `statistical_analyses[]` row):**

| Field | Type | Required | Notes |
|---|---|---|---|
| `timepoint` | str | yes | Human-readable label, e.g. `"120 min"`, `"48h"`, `"acclimated"`, `"unknown"` |
| `timepoint_hours` | number \| null | yes | Numeric conversion in hours; `null` only when paper truly doesn't report |
| `growth_phase` | enum str | yes | One of the canonical values below |

**Canonical `growth_phase` enum** (extendable):

| Value | When to use |
|---|---|
| `exponential` | Cells dividing actively (early/mid/late-log collapsed). Default for most short-exposure stress studies. |
| `stationary` | Nutrient-replete stationary phase. Rare in these papers. |
| `nutrient_limited` | Cells arrested due to N/P/Fe/C limitation. Distinct from classic stationary. |
| `acclimated_steady_state` | Grown ≥5–10 generations under the treatment condition (common in pCO2 / temperature / light-intensity acclimation studies). |
| `infected` | Host cells after viral/phage infection; distinct physiology from exponential or stationary. |
| `recovery` | Cells recovering after stress relief (e.g., nutrient re-addition, light restoration post-darkness). |
| `diel` | Cells on a normal light-dark cycle; implies periodic physiology. Use when the paper's sampling is tied to the diel cycle rather than a single snapshot. Not split into dark/light sub-phases. |
| `darkness` | Cells in **extended/prolonged darkness** — distinct from the dark half of a diel cycle. Used in extended-darkness survival experiments. |
| `death` | Late decline / cell-death phase after stationary. Used for long time-courses that sample past viable-cell peak. |
| `acute_stress` | Cells still morphologically exponential (dividing) but transcriptionally in shock response — short-exposure oxidative / osmotic / heat / chemical stress. Narrow scope: **do not** use when a more specific phase (`nutrient_limited`, `infected`, `darkness`) fits. |
| `unknown` | Paper truly does not state the phase at sampling. |

**Escape hatch for novelty:** If the paper clearly describes a phase that doesn't map to any canonical value, the LLM outputs `other:<short_slug>` (e.g., `other:heat_acclimated`, `other:late_decline`). The validator accepts the `other:` prefix; the extraction report aggregates all `other:*` values across the corpus so recurring ones can be promoted into the enum.

**`unknown` vs `other:*`:** A hard semantic distinction.

- `unknown` → the paper gives no information about physiological state at sampling. Nothing to promote.
- `other:<slug>` → the paper *does* describe a state, but it doesn't fit the current enum. These are candidates for enum extension.

**Rationale for per-analysis placement:** `growth_phase` genuinely varies across timepoints of a single experiment. In Weissberg 2025 N-starvation, cells are `exponential` at t=4h and `nutrient_limited` at t=24h and beyond. Placing `growth_phase` on the analysis row (alongside `timepoint_hours`) allows this variation; per-experiment placement would not.

### 2. Propagation to the knowledge graph

**`config/schema_config.yaml`:**

```yaml
changes_expression_of:
  properties:
    ...existing...
    growth_phase: str

experiment:
  properties:
    ...existing...
    growth_phases: str[]   # distinct values across child analyses
```

**`multiomics_kg/adapters/omics_adapter.py`:**

Extend per-edge property dict — add `growth_phase` alongside `time_point`, `time_point_order`, `time_point_hours`. One-line change per edge-construction site. Sanitize via existing `_clean_str()`.

**`scripts/post-import.sh` and `scripts/post-import.cypher`:**

New aggregation step (runs after existing experiment computations):

```cypher
MATCH (e:Experiment)-[r:Changes_expression_of]->(:Gene)
WITH e, collect(DISTINCT r.growth_phase) AS phases
SET e.growth_phases = phases;
```

### 3. Validator changes (`.claude/skills/paperconfig/validate_paperconfig.py`)

New checks on every `statistical_analyses[]` entry:

- `timepoint` required — warn if missing; `"unknown"` accepted.
- `timepoint_hours` required — warn if `null`; allowed when paper truly doesn't report.
- `growth_phase` required — **error** if missing; **error** if not in `VALID_GROWTH_PHASES`.

New constant:

```python
VALID_GROWTH_PHASES = {
    "exponential",
    "stationary",
    "nutrient_limited",
    "acclimated_steady_state",
    "infected",
    "recovery",
    "diel",
    "darkness",
    "death",
    "acute_stress",
    "unknown",
}

# Accepts either a canonical enum value or an `other:<slug>` escape value
def is_valid_growth_phase(value: str) -> bool:
    return value in VALID_GROWTH_PHASES or (
        value.startswith("other:") and len(value) > len("other:")
    )
```

Sanity warning on experiments: `is_time_course=true` but only one analysis references it (or vice versa — 2+ analyses but `is_time_course=false`).

### 4. Extraction pipeline

New module: `multiomics_kg/extraction/timepoint/`

```
extraction/timepoint/
├── extract.py              # CLI; reads paperconfig + PDF, emits JSON proposal next to paperconfig
├── merge.py                # applies approved JSON → paperconfig.yaml via ruamel.yaml round-trip
├── extraction_utils.py     # JSON I/O, paperconfig parsing, analysis lookup, safety checks
└── prompts.py              # SHARED_RULES + per-experiment-type hint templates
```

**Why the intermediate JSON (reversing the earlier decision):** The extraction produces more than just the three paperconfig fields — it also captures evidence quotes, confidence bands, and `other:*` escape values for enum-evolution triage. Committing the JSON gives a durable audit trail tied to each paperconfig's history; the review loop inspects the JSON (not just a git diff on the yaml) so reviewers see the evidence without re-reading PDFs.

**Input package to the LLM (per paper):**

1. **Background context** (text):
   - `publication.papername`, `publication.doi`.
   - Full `experiments` block from the paperconfig (treatment/control conditions, treatment_type, background_factors, experimental_context, medium, temperature, light — all the metadata the LLM needs to reason about phase).
   - Publication metadata from `cache/pdf_extraction_cache.json` if an entry exists for `papermainpdf`: `title`, `abstract`, `description`, `study_type`. Used for orientation and to reduce token waste on basics the LLM would otherwise re-parse from the PDF. Optional — if cache miss, skipped.

2. **Extraction targets** (structured list, per analysis):
   - `id`, `experiment_key`.
   - Existing values for `timepoint`, `timepoint_hours`, `growth_phase` (null/absent when not set).
   - `fields_to_fill: [...]` — the subset of the three fields that need extraction for this analysis (derived from the preprocess; field-level, not row-level).
   - Analyses where all three fields are already set are **omitted entirely** from the prompt (unless `--validate`).

3. **PDFs** (multimodal attachments):
   - `papermainpdf` — always included.
   - `extraction.additional_pdfs` from the paperconfig — same convention as the cluster-extraction pipeline (e.g. Coe 2024 supplement). Included when declared.

**Extraction artifact — written to `<paper_dir>/extractions/timepoint.json`:**

Per-paper extractions live in a new `extractions/` subdirectory next to the `paperconfig.yaml`. This is a new convention going forward — we anticipate additional extraction types (beyond cluster and timepoint) so each paper gets one home directory for all LLM-derived metadata. Existing `cluster_extractions/` directories stay where they are for now; a future refactor may consolidate them under `extractions/cluster/`.

Example layout:

```
data/Prochlorococcus/papers_and_supp/Tetu 2019/
├── paperconfig.yaml
├── s42003-019-0410-x.pdf
├── extractions/
│   └── timepoint.json      ← new
└── ... (supplementary CSVs)
```

**What the LLM actually returns** (per analysis — the payload, not the whole JSON file):

```json
{
  "analysis_id": "DE_hdpe_leachate_vs_control_MIT9312",
  "timepoint": "120 min",
  "timepoint_hours": 2.0,
  "growth_phase": "exponential",
  "self_assessment": "high",
  "assessment_notes": "",
  "supporting_quotes": [
    {"quote": "RNA was extracted 120 min post-exposure from mid-exponential cultures at OD 0.1", "location": "Methods §2.3"}
  ],
  "source_figures": []
}
```

The prompt asks for exactly these fields; the LLM does not emit metadata, lookups, or scope echoes.

**Rules the LLM follows for its output fields:**

- Only emit the three data fields (`timepoint`, `timepoint_hours`, `growth_phase`) that are in `fields_requested` for this analysis. Fields not requested must not appear in the response.
- `self_assessment` ∈ {`high`, `medium`, `low`}.
- `supporting_quotes` — verbatim from the paper; at least one quote required unless `self_assessment: low` with `assessment_notes` explaining why no direct quote exists.
- `source_figures` — figure/table numbers used, empty array if purely text-derived.
- `growth_phase` — either canonical enum value, `unknown`, or `other:<slug>`.

**What `extract.py` wraps around the LLM output (post-call, before writing the JSON):**

- `metadata` block: `paper`, `doi`, `model`, `extracted_at` (ISO timestamp), `input_tokens`, `output_tokens`. Bookkeeping/audit trail, not requested from the LLM.
- `experiment_key` per analysis: deterministic lookup from the paperconfig by `analysis_id`. Included in the JSON for human readability; no reason to round-trip through the LLM.
- `fields_requested` per analysis: set by the preprocess. Echoed into the JSON so reviewers see scope at a glance.

**Final on-disk shape:**

```json
{
  "metadata": {
    "paper": "Tetu 2019",
    "doi": "10.1038/s42003-019-0410-x",
    "model": "gpt-4.1-mini",
    "extracted_at": "2026-04-12T14:32:11",
    "input_tokens": 12450,
    "output_tokens": 890
  },
  "analyses": [
    {
      "analysis_id": "DE_hdpe_leachate_vs_control_MIT9312",
      "experiment_key": "plastic_hdpe_plastic_leachate_50_mit9312_rnaseq",
      "fields_requested": ["timepoint", "timepoint_hours", "growth_phase"],
      "timepoint": "120 min",
      "timepoint_hours": 2.0,
      "growth_phase": "exponential",
      "self_assessment": "high",
      "assessment_notes": "",
      "supporting_quotes": [
        {"quote": "RNA was extracted 120 min post-exposure from mid-exponential cultures at OD 0.1", "location": "Methods §2.3"}
      ],
      "source_figures": []
    }
  ]
}
```

This schema borrows the `metadata`/`supporting_quotes`/`source_figures`/`self_assessment` conventions from the cluster-extraction pipeline (`data/.../zinser 2009/cluster_extractions/med4_diel_clusters.json`) so both extractions share vocabulary for reviewers.

**Review workflow:**

0. (Optional) `extract.py --paper "Tetu 2019" --dry-run` — no LLM call, just prints per-analysis `fields_to_fill` lists. Used to confirm scope before a real run.
1. `extract.py --paper "Tetu 2019"` runs preprocess → calls LLM with the package above → writes `data/.../Tetu 2019/extractions/timepoint.json` (creates `extractions/` subdir if needed).
2. Aggregate report `data/timepoint_extraction_report.md` updates with: confidence breakdown, `unknown` list, sorted `other:*` frequencies, evidence quotes for each flagged item.
3. Reviewer inspects JSON per paper — evidence quotes are inline, no need to re-read PDFs. Edits in place if extraction got it wrong.
4. `merge.py --paper "Tetu 2019"` writes the three paperconfig fields (`timepoint`, `timepoint_hours`, `growth_phase`) into `paperconfig.yaml` via `ruamel.yaml` round-trip.
5. Validator runs against updated paperconfig; reviewer commits both the JSON and the paperconfig diff.

**CLI flags:**

| Flag | Purpose |
|---|---|
| `--paper "<name>"` / `--all` | Scope selection. |
| `--dry-run` | Preprocess only — list per-analysis `fields_to_fill`, no LLM call. |
| `--validate` | Ask LLM about **every** field regardless of existing values. Emits existing vs proposed in JSON. `merge.py` surfaces mismatches as warnings but does not overwrite without `--force`. |
| `--fields timepoint,timepoint_hours,growth_phase` | Restrict extraction to specific fields. Useful when re-running after an enum change: `--fields growth_phase` only. |
| `--force` (on `merge.py`) | Allow overwrite of existing paperconfig values. Never default. |

**Preprocess (runs inside `extract.py` before the LLM call):**

For every analysis in the paperconfig, compute `fields_to_fill`:

| Existing state | Behavior (default) | Behavior (`--validate`) |
|---|---|---|
| field absent | include in `fields_to_fill` | include |
| `timepoint_hours: null` | include (null is a "needs filling" signal) | include |
| `timepoint` is empty string | include | include |
| field has a non-null, non-empty value | skip | include (LLM emits existing+proposed for comparison) |

An analysis with empty `fields_to_fill` (all three set) is omitted from the prompt entirely in default mode. In `--validate`, every analysis is included and every field re-examined.

This keeps re-runs cheap (only the missing subset goes to the LLM) and prevents the LLM from "helpfully" revisiting values a human already curated.

**Matching key:** `analysis_id` alone is sufficient. The paperconfig validator already guarantees analysis IDs are unique within a paperconfig. `merge.py` walks `publication.supplementary_materials.<table>.statistical_analyses[]` looking for matches by `id`.

**Safety checks:**

*In `extract.py`:*
1. Fail if an LLM-proposed `analysis_id` does not exist in the paperconfig (hallucinated ID).
2. Fail if any paperconfig analysis has no proposal (LLM missed a row).

*In `merge.py`:*
3. Warn (do not overwrite) if an analysis already has `growth_phase` set and the JSON value differs — re-runs never silently clobber human edits. Reviewer must pass `--force` to overwrite.
4. Reject any `growth_phase` value that's neither a canonical enum value nor an `other:<slug>` — catches malformed JSON before it enters the paperconfig.

**Prompt design** (key rules):
- `SHARED_RULES`: enum values; "unknown is always better than wrong" (carries forward the no-guessing memory from the cluster-extraction project); must quote evidence from the methods section; confidence band (`high`/`medium`/`low`) required; low-confidence → `growth_phase: unknown`.
- **Escape hatch rule:** If the paper describes a state not in the enum, emit `other:<slug>` rather than forcing a bad fit. `other:*` is preferred over `unknown` whenever the paper gives *any* positional information. Slugs are short snake_case (e.g., `other:heat_acclimated`, `other:late_decline`).
- Per-experiment-type hints via branching on `omics_type` + `treatment_type`:
  - Diel studies (`treatment_type: [diel]`) → `diel` for cycling-phase samples.
  - Extended darkness (`treatment_type: [darkness]`, prolonged dark exposure) → `darkness`. Do not confuse with the dark half of a normal diel cycle (which stays `diel`).
  - Acclimation / chronic exposure (`treatment_type: [carbon]` with ≥5 generations language) → `acclimated_steady_state`.
  - Nutrient starvation (`treatment_type: [nitrogen|phosphorus|iron]`) → `exponential` at early timepoints, `nutrient_limited` once depletion is confirmed, `death` if sampling extends past viable-cell peak.
  - Phage / viral infection → `exponential` at t=0 (pre-infection), `infected` post-infection timepoints.
  - Rescue / re-addition experiments → `recovery` for post-intervention timepoints.
  - Short-exposure stress (≤6 h) at still-dividing cells → `acute_stress` (only when no more-specific phase fits; `nutrient_limited`, `infected`, `darkness` take precedence if applicable).

### 5. Rollout

1. **Schema + validator + adapter + post-import changes first.** Land as an inert PR — new field accepts missing values with warnings, no KG data produced yet. Run `pytest -m "not slow and not kg"` and omics-edge-snapshot before/after to confirm no regressions. Paperconfigs continue to validate (warnings only).
2. **Build the extraction tool** in `multiomics_kg/extraction/timepoint/`.
3. **Run extraction in batches**, grouped by organism for review efficiency:
   - Batch 1: Prochlorococcus coculture papers (Biller 2018, Hennon 2017, Aharonovich 2016, …).
   - Batch 2: Prochlorococcus nutrient-stress papers (Weissberg, Read, Martiny, Thompson 2011).
   - Batch 3: Prochlorococcus diel + short-exposure (Zinser, Thompson 2016, Tetu, Fang).
   - Batch 4: Prochlorococcus infection + remaining (Lindell, Lin, Anjur, …).
   - Batch 5: Synechococcus papers (Tal, Beliaev, Kratzl, Ma, Oleza, Kaur, Bernstein).
4. **Per-batch review loop:**
   - Run `uv run python -m multiomics_kg.extraction.timepoint.extract --paper "<name>"`.
   - Inspect `git diff <paperconfig>`.
   - Edit fields in place if the LLM got something wrong; commit.
5. **Rebuild the KG** once all paperconfigs updated.
6. **New KG validity test** (`tests/kg_validity/test_expression.py`): every `Changes_expression_of` edge has non-null `growth_phase`; every `Experiment` has `growth_phases` array of length ≥1.
7. **Docs updates:**
   - `CLAUDE.md` — add `growth_phase` to the edge-property list; add `growth_phases` to the Experiment node property list.
   - `.claude/skills/paperconfig/SKILL.md` — update the "Required Fields" table, add a Growth Phase section with the enum and decision guidance.
   - `MEMORY.md` — add a project memory summarizing the rollout.

## Review UI

Two-surface review:
1. **`extractions/timepoint.json`** — the primary review artifact. Contains evidence quotes, confidence, and the proposed fields. Reviewer edits in place if extraction got something wrong, then runs `merge.py`.
2. **`paperconfig.yaml` git diff** — confirmation that only the three intended fields (`timepoint`, `timepoint_hours`, `growth_phase`) changed, with expected values.

Both artifacts are committed together per paper.

## Data flow diagram

```
┌────────────────────┐        ┌──────────────────────┐
│ paperconfig.yaml   │───────▶│ extract.py           │
│ (experiments +     │        │  • preprocess:       │
│  analyses)         │        │    compute per-      │
│                    │        │    analysis          │
│ pdf_extraction_    │───────▶│    fields_to_fill    │
│   cache.json       │        │  • render prompt     │
│ (optional, for     │        │    (background +     │
│  orientation)      │        │     targets + PDFs)  │
│                    │        │  • call LLM          │
│ papermainpdf +     │───────▶│  • validate IDs      │
│ additional_pdfs    │        │  • write JSON        │
└────────────────────┘        └──────────┬───────────┘
                                         │
                                         ▼
                              ┌──────────────────────┐
                              │ extractions/         │  ◀── human reviews,
                              │   timepoint.json     │      edits in place
                              │ (evidence quotes,    │      if needed
                              │  confidence, other:* │
                              │  escapes)            │
                              └──────────┬───────────┘
                                         │
                                         ▼
                              ┌──────────────────────┐
                              │ merge.py             │
                              │  • ruamel yaml       │
                              │    round-trip        │
                              │  • overwrite guard   │
                              └──────────┬───────────┘
                                         │
                                         ▼
                              ┌──────────────────────┐
                              │ paperconfig.yaml     │  ◀── human reviews git diff
                              │ (timepoint +         │      commits both JSON
                              │  growth_phase        │      and yaml
                              │  populated)          │
                              └──────────┬───────────┘
                                         │
                    (KG rebuild)         ▼
                              ┌──────────────────────┐
                              │ omics_adapter        │
                              │  → edge.growth_phase │
                              └──────────┬───────────┘
                                         ▼
                              ┌──────────────────────┐
                              │ post-import.cypher   │
                              │  → experiment.       │
                              │     growth_phases[]  │
                              └──────────────────────┘
```

## Enum iteration plan

The `VALID_GROWTH_PHASES` enum is an initial set (7 values including `unknown`). We expect it to grow as the `other:*` escape hatch surfaces recurring novel states. Vocabulary refinement is a planned first-class part of the rollout, not a one-off.

**Iteration loop (cheap — no LLM re-extraction):**

1. **Extract a batch** (organism-grouped, per rollout section). Writes per-paper JSON proposals + updates the corpus-wide `data/timepoint_extraction_report.md`.
2. **Review the report for `other:*` clusters and `unknown` cases:**
   - `other:<slug>` frequencies sorted desc → any slug appearing in ≥2 papers is a strong candidate for promotion.
   - `unknown` count trending high → may signal a missing enum value the LLM couldn't even name; inspect evidence quotes for patterns.
   - Low-confidence canonical extractions (`exponential` confidence=low with ambiguous evidence) → may signal shoehorning.
3. **Propose enum extensions.** Before promoting a slug, sanity-check:
   - Biologically meaningful (implies distinct transcriptional state)?
   - Appears in ≥2 papers?
   - No overlap with an existing value (if yes, split or merge consciously)?
4. **Update the enum** in one place:
   - `VALID_GROWTH_PHASES` in the validator.
   - Enum table in this spec (mark added values with the batch number that prompted them).
   - `.claude/skills/paperconfig/SKILL.md` growth-phase section.
5. **Remap paperconfigs in place** — no LLM re-run. A small `remap.py` script (or manual edit) replaces `other:<promoted_slug>` with `<promoted_slug>` across all paperconfigs + their extraction JSONs. Evidence quotes stay intact; the reviewer doesn't need to re-read methods.
6. **Commit the enum change** alongside the remap diffs. The promoted slug's first appearance (batch number, paper, evidence quote) is noted in the spec's enum table.

**Initial enum covers the predictable cases** (`exponential`, `stationary`, `nutrient_limited`, `acclimated_steady_state`, `infected`, `recovery`, `diel`, `darkness`, `death`, `acute_stress`). Not pre-added, to be surfaced via `other:*` only if the data demands them:

- Splitting `exponential` into `early_exponential` / `late_exponential` — not worth the granularity unless queries actually need it.
- Splitting `diel` into `diel_light` / `diel_dark` — not worth it for our use case; `diel` stays unified.

**Freeze the enum when:** two consecutive batches complete without any new `other:*` slug appearing in ≥2 papers. Document the final enum in CLAUDE.md and MEMORY.md.

## Open questions (none blocking)

- Do we want to publish the extraction report as a committed artifact (`data/timepoint_extraction_report.md`) or keep it local? (Proposal: commit it, matches the cluster-extraction precedent.)
- Should `growth_phase_description` be indexed in the full-text `experimentFullText` index? (Proposal: no — it's per-edge free text, not per-experiment; skipping avoids index bloat.)
