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
| `growth_phase_description` | str | no | Free-text nuance (e.g. "mid-exponential at OD ~0.1 when leachate added") |

**Canonical `growth_phase` enum** (extendable):

| Value | When to use |
|---|---|
| `exponential` | Cells dividing actively (early/mid/late-log collapsed). Default for most short-exposure stress studies. |
| `stationary` | Nutrient-replete stationary phase. Rare in these papers. |
| `nutrient_limited` | Cells arrested due to N/P/Fe/C limitation. Distinct from classic stationary. |
| `acclimated_steady_state` | Grown ≥5–10 generations under the treatment condition (common in pCO2 / temperature / light-intensity acclimation studies). |
| `unknown` | Paper does not state the phase at sampling. |

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
- `growth_phase_description` optional — no validation.

New constant:

```python
VALID_GROWTH_PHASES = {
    "exponential",
    "stationary",
    "nutrient_limited",
    "acclimated_steady_state",
    "unknown",
}
```

Sanity warning on experiments: `is_time_course=true` but only one analysis references it (or vice versa — 2+ analyses but `is_time_course=false`).

### 4. Extraction pipeline

New module: `multiomics_kg/extraction/timepoint/`

```
extraction/timepoint/
├── extract.py              # CLI; reads paperconfig + PDF, writes proposed fields into paperconfig.yaml
├── extraction_utils.py     # ruamel.yaml round-trip, analysis lookup, safety checks
├── prompts.py              # SHARED_RULES + per-experiment-type hint templates
└── (no merge.py — extract writes direct to paperconfig; git diff is the review UI)
```

**Input to the LLM per paper:**
- `papermainpdf` + any `additional_pdfs`.
- A rendered snapshot of the paperconfig's `experiments` block and list of `(analysis_id, experiment_key, existing_timepoint_fields)` — so the LLM sees what analyses exist and their scientific context.

**LLM output** (internal; not persisted to disk):

```json
{
  "paper": "Tetu 2019",
  "doi": "10.1038/s42003-019-0410-x",
  "analyses": [
    {
      "analysis_id": "DE_hdpe_leachate_vs_control_MIT9312",
      "experiment_key": "plastic_hdpe_plastic_leachate_50_mit9312_rnaseq",
      "timepoint": "120 min",
      "timepoint_hours": 2.0,
      "growth_phase": "exponential",
      "growth_phase_description": "mid-exponential cells harvested 120 min after HDPE leachate addition",
      "confidence": "high",
      "evidence": "Methods p.3: 'RNA was extracted 120 min post-exposure from mid-exponential cultures at OD ...'",
      "notes": ""
    }
  ]
}
```

**Output to disk:** `extract.py` walks the paperconfig with `ruamel.yaml` (preserves comments and ordering) and writes the four fields directly into each `statistical_analyses[]` entry by `id`. A summary report `data/timepoint_extraction_report.md` lists confidence breakdown, flagged-for-review items, and per-paper diff stats.

**Matching key:** `analysis_id` alone is sufficient. The paperconfig validator already guarantees analysis IDs are unique within a paperconfig. `extract.py` walks `publication.supplementary_materials.<table>.statistical_analyses[]` looking for matches by `id`.

**Safety checks in `extract.py`:**
1. Fail if an LLM-proposed `analysis_id` does not exist in the paperconfig (hallucinated ID).
2. Fail if any paperconfig analysis has no proposal (LLM missed a row).
3. Warn (do not overwrite) if an analysis already has `growth_phase` set and the new value differs — re-runs never silently clobber human edits.

**Prompt design** (key rules):
- `SHARED_RULES`: enum values; "unknown is always better than wrong" (carries forward the no-guessing memory from the cluster-extraction project); must quote evidence from the methods section; confidence band (`high`/`medium`/`low`) required; low-confidence → `growth_phase: unknown`.
- Per-experiment-type hints via branching on `omics_type` + `treatment_type`:
  - Diel → default `growth_phase: exponential` unless methods explicitly say otherwise.
  - Acclimation / chronic exposure (`treatment_type: [carbon]` with ≥5 generations language) → `acclimated_steady_state`.
  - Nutrient starvation (`treatment_type: [nitrogen|phosphorus|iron]`) → `exponential` at early timepoints, `nutrient_limited` once depletion is confirmed.
  - Phage / viral infection → `exponential` host at infection.
  - Short-exposure stress (≤6 h) → default `exponential`.

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

Git diff on each paperconfig. No sidecar JSON.

## Data flow diagram

```
┌────────────────────┐        ┌──────────────────────┐
│ paperconfig.yaml   │───────▶│ extract.py           │
│ (experiments +     │        │  • read paperconfig  │
│  analyses)         │        │  • render LLM prompt │
│                    │        │  • call LLM          │
│ paper PDF(s)       │───────▶│  • validate IDs      │
└────────────────────┘        │  • ruamel write-back │
                              └──────────┬───────────┘
                                         │
                                         ▼
                              ┌──────────────────────┐
                              │ paperconfig.yaml     │  ◀── human reviews git diff
                              │ (timepoint +         │      commits or reverts
                              │  growth_phase        │
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

## Open questions (none blocking)

- Do we want to publish the extraction report as a committed artifact (`data/timepoint_extraction_report.md`) or keep it local? (Proposal: commit it, matches the cluster-extraction precedent.)
- Should `growth_phase_description` be indexed in the full-text `experimentFullText` index? (Proposal: no — it's per-edge free text, not per-experiment; skipping avoids index bloat.)
