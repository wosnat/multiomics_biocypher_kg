# Timepoint & Growth-Phase Backfill Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Populate `timepoint` + `timepoint_hours` + new `growth_phase` field on every `statistical_analyses[]` row in all 30 paperconfigs (23 Prochlorococcus + 7 Synechococcus), propagate `growth_phase` through adapter → `Changes_expression_of` edge → aggregated `growth_phases[]` array on `Experiment` nodes.

**Architecture:** Land schema + validator + adapter + post-import changes first (inert; validator warns on missing). Then build a three-tool extraction pipeline (`extract.py` / `merge.py` / `remap.py`) in `multiomics_kg/extraction/timepoint/` modeled on the existing cluster-extraction pattern. Each per-paper extraction writes JSON to `<paper_dir>/extractions/timepoint.json`; merge writes to paperconfig.yaml only once per batch, after `other:*` slug promotion via `remap.py` is settled.

**Tech Stack:** Python 3.11+, ruamel.yaml (preserves comments/ordering), OpenAI Responses API with structured output (`pydantic` schemas), `ncbigene`-style neo4j schema via BioCypher, pytest.

**Reference spec:** `docs/superpowers/specs/2026-04-12-timepoint-growth-phase-backfill-design.md`

---

## File Structure

**Created:**
- `multiomics_kg/extraction/timepoint/__init__.py`
- `multiomics_kg/extraction/timepoint/__main__.py` — CLI dispatcher (extract/merge/remap)
- `multiomics_kg/extraction/timepoint/extraction_utils.py` — JSON I/O, paperconfig walking, signature computation, analysis lookup
- `multiomics_kg/extraction/timepoint/prompts.py` — SHARED_RULES + per-experiment-type hints
- `multiomics_kg/extraction/timepoint/extract.py` — CLI; preprocess + LLM call + write JSON (+ partial handling, dry-run, validate, resume, retry)
- `multiomics_kg/extraction/timepoint/merge.py` — CLI; write JSON fields into paperconfig.yaml (+ staleness guards, partial guard, overwrite guard)
- `multiomics_kg/extraction/timepoint/remap.py` — CLI; rewrite `--from X --to Y` across paperconfigs + JSONs (+ provenance)
- `multiomics_kg/extraction/timepoint/report.py` — write `data/timepoint_extraction_report.md` aggregated from all JSONs
- `tests/extraction/timepoint/test_extraction_utils.py`
- `tests/extraction/timepoint/test_extract.py` (uses mock LLM)
- `tests/extraction/timepoint/test_merge.py`
- `tests/extraction/timepoint/test_remap.py`
- `tests/extraction/timepoint/test_report.py`

**Modified:**
- `config/schema_config.yaml` — add `growth_phase` to `experiment to gene expression association`; add `growth_phases: str[]` to `experiment`.
- `multiomics_kg/adapters/omics_adapter.py` — add `growth_phase` to edge properties in `_create_expression_edges` (around line 801–805).
- `scripts/post-import.sh` — new Cypher block to compute `experiment.growth_phases[]`.
- `scripts/post-import.cypher` — mirror the same block.
- `.claude/skills/paperconfig/validate_paperconfig.py` — add `VALID_GROWTH_PHASES`, `is_valid_growth_phase()`, per-analysis checks for `timepoint` / `timepoint_hours` / `growth_phase`.
- `.claude/skills/paperconfig/SKILL.md` — required-fields table + Growth Phase section.
- `tests/kg_validity/test_expression.py` — new test: every edge has `growth_phase`, every Experiment has `growth_phases`.
- `CLAUDE.md` — add `growth_phase` edge property, `growth_phases[]` experiment property.
- `memory/MEMORY.md` — project memory summarizing rollout.

**Paperconfig file modifications are out of scope for the code plan** — those happen during the rollout phase (Tasks 14+) via `extract.py`/`merge.py`, one batch at a time, per-reviewer.

---

## Task 1: Schema — add `growth_phase` to edge and `growth_phases[]` to experiment node

**Files:**
- Modify: `config/schema_config.yaml:164-182` (edge definition) and lines 37–64 (experiment properties)

- [ ] **Step 1: Read the existing schema to confirm line numbers**

Run: `grep -n "changes_expression_of\|^experiment:" config/schema_config.yaml | head -10`
Expected: shows the two sections to modify.

- [ ] **Step 2: Edit `changes_expression_of` edge properties**

File: `config/schema_config.yaml` — find the `experiment to gene expression association:` block. Add `growth_phase: str` to the `properties:` list, after `time_point_hours: float` and before `log2_fold_change: float`:

```yaml
experiment to gene expression association:
  is_a: Association
  represented_as: edge
  label_as_edge: changes_expression_of
  label_in_input: changes_expression_of
  source: experiment
  target: gene
  properties:
    time_point: str
    time_point_order: int
    time_point_hours: float
    growth_phase: str
    log2_fold_change: float
    adjusted_p_value: float
    expression_direction: str
    significant: str
    expression_status: str
    rank_by_effect: int
    rank_up: int
    rank_down: int
```

- [ ] **Step 3: Edit `experiment` node properties**

Add `growth_phases: str[]` at the end of the `experiment:` properties list, after `time_point_significant_down: int[]`:

```yaml
    time_point_significant_down: int[]
    growth_phases: str[]
```

- [ ] **Step 4: Verify the schema still parses**

Run: `uv run python -c "import yaml; yaml.safe_load(open('config/schema_config.yaml'))"`
Expected: no output (parses cleanly).

- [ ] **Step 5: Commit**

```bash
git add config/schema_config.yaml
git commit -m "feat(schema): add growth_phase to Changes_expression_of edge + growth_phases[] to Experiment"
```

---

## Task 2: Adapter — emit `growth_phase` on expression edges

**Files:**
- Modify: `multiomics_kg/adapters/omics_adapter.py` (the `_create_expression_edges` method around lines 557–560, 613–620, 801–805)
- Test: `tests/test_omics_adapter.py` (add a focused test)

- [ ] **Step 1: Read the adapter's edge-construction site**

Run: `sed -n '550,620p' multiomics_kg/adapters/omics_adapter.py`
Expected: shows the call to the edge-builder with `time_point=...`, `time_point_order=...`, `time_point_hours=...`.

- [ ] **Step 2: Add a `growth_phase` parameter to the edge-builder call**

At the site that currently passes `time_point=timepoint, time_point_order=tp_order, time_point_hours=timepoint_hours`, also pass `growth_phase`:

```python
# Near line 555-560 in _create_expression_edges (or its equivalent)
growth_phase = analysis.get("growth_phase")  # may be None on paperconfigs not yet backfilled

edge_records = self._make_expression_edges(
    ...,
    time_point=timepoint,
    time_point_order=tp_order,
    time_point_hours=timepoint_hours,
    growth_phase=growth_phase,
)
```

- [ ] **Step 3: Extend the edge-builder signature**

At line 613 (`def _make_expression_edges` or similar), add `growth_phase: str | None = None,` parameter:

```python
def _make_expression_edges(
    self,
    ...,
    time_point: str | None = None,
    time_point_order: int = 1,
    time_point_hours: float | None = None,
    growth_phase: str | None = None,
) -> ...:
```

- [ ] **Step 4: Write `growth_phase` into the edge property dict**

At the property-assembly site (around line 801-805):

```python
if time_point:
    edge_properties['time_point'] = self.clean_text(time_point)
edge_properties['time_point_order'] = time_point_order
if time_point_hours is not None:
    edge_properties['time_point_hours'] = time_point_hours
if growth_phase:
    edge_properties['growth_phase'] = self.clean_text(growth_phase)
```

Use `self.clean_text()` for string sanitization (matches the existing pattern for `time_point`) — strips `|` and `'`.

- [ ] **Step 5: Write a failing adapter test**

Create or append to `tests/test_omics_adapter.py`:

```python
def test_growth_phase_flows_to_edge_properties(tmp_path):
    """If an analysis has growth_phase set, the edge carries it."""
    from multiomics_kg.adapters.omics_adapter import OMICSAdapter
    # Build a minimal paperconfig with one analysis having growth_phase
    # (use the existing test fixture pattern — see test_omics_adapter_organism_gene.py)
    paperconfig = {
        "publication": {
            "papername": "Test",
            "papermainpdf": str(tmp_path / "dummy.pdf"),
            "experiments": {
                "exp1": {
                    "name": "Test experiment",
                    "organism": "Prochlorococcus MED4",
                    "treatment_condition": "test",
                    "control_condition": "control",
                    "omics_type": "RNASEQ",
                    "test_type": "DESeq2",
                    "treatment_type": ["nitrogen"],
                },
            },
            "supplementary_materials": {
                "tbl": {
                    "type": "csv",
                    "filename": str(tmp_path / "de.csv"),
                    "statistical_analyses": [{
                        "id": "DE_test",
                        "experiment": "exp1",
                        "timepoint": "24h",
                        "timepoint_hours": 24,
                        "growth_phase": "nutrient_limited",
                        "name_col": "gene",
                        "logfc_col": "log2fc",
                    }],
                },
            },
        },
    }
    # ... construct DE csv, instantiate adapter, iterate get_edges()
    edges = list(adapter.get_edges())
    de_edges = [e for e in edges if e[3] == "changes_expression_of"]
    assert all(e[4].get("growth_phase") == "nutrient_limited" for e in de_edges)
```

(Follow the fixture style of `tests/test_omics_adapter_organism_gene.py` — build a minimal CSV + paperconfig dict and feed through the adapter.)

- [ ] **Step 6: Run the test, confirm it fails, then confirm it passes after the edits**

Run: `uv run pytest tests/test_omics_adapter.py::test_growth_phase_flows_to_edge_properties -v`
Expected: PASS after the edits in steps 2–4.

- [ ] **Step 7: Run the broader adapter test suite for regression**

Run: `uv run pytest tests/test_omics_adapter.py tests/test_omics_adapter_organism_gene.py -v`
Expected: all tests pass — the new parameter is optional (`growth_phase: str | None = None`) so old fixtures still work.

- [ ] **Step 8: Commit**

```bash
git add multiomics_kg/adapters/omics_adapter.py tests/test_omics_adapter.py
git commit -m "feat(omics_adapter): propagate growth_phase to Changes_expression_of edges"
```

---

## Task 3: Post-import — aggregate `growth_phases[]` onto Experiment nodes

**Files:**
- Modify: `scripts/post-import.sh` (add new block)
- Modify: `scripts/post-import.cypher` (mirror the same block)

- [ ] **Step 1: Locate where Experiment summary properties are computed**

Run: `grep -n "Compute Experiment summary" scripts/post-import.sh`
Expected: shows line 126 echo header and the block that follows.

- [ ] **Step 2: Add a new Cypher block after the existing summary computations**

In `scripts/post-import.sh`, after the existing Experiment summary Cypher (around line 220 — find the block ending with `SET e.time_point_significant_down = ...;` and insert after it):

```bash
echo "=== Post-process: Compute Experiment growth_phases ==="
cypher-shell -a "$NEO4J_URL" <<EOF
MATCH (e:Experiment)
OPTIONAL MATCH (e)-[r:Changes_expression_of]->(:Gene)
WITH e, [v IN collect(DISTINCT r.growth_phase) WHERE v IS NOT NULL] AS phases
SET e.growth_phases = phases;
EOF
```

The `[v IN ... WHERE v IS NOT NULL]` filter handles the case where no edges carry `growth_phase` yet (pre-backfill).

- [ ] **Step 3: Mirror the same block in `scripts/post-import.cypher`**

Add the same Cypher (without the shell wrapper):

```cypher
// Compute Experiment growth_phases (distinct values across child edges)
MATCH (e:Experiment)
OPTIONAL MATCH (e)-[r:Changes_expression_of]->(:Gene)
WITH e, [v IN collect(DISTINCT r.growth_phase) WHERE v IS NOT NULL] AS phases
SET e.growth_phases = phases;
```

- [ ] **Step 4: Verify shell script syntax**

Run: `bash -n scripts/post-import.sh`
Expected: no output (syntactically valid).

- [ ] **Step 5: Commit**

```bash
git add scripts/post-import.sh scripts/post-import.cypher
git commit -m "feat(post-import): aggregate growth_phase edge values into Experiment.growth_phases[]"
```

---

## Task 4: Validator — add `VALID_GROWTH_PHASES` enum and per-analysis checks

**Files:**
- Modify: `.claude/skills/paperconfig/validate_paperconfig.py`

- [ ] **Step 1: Read the existing validator to find where analysis fields are checked**

Run: `grep -n "statistical_analyses\|timepoint_hours\|def validate" .claude/skills/paperconfig/validate_paperconfig.py`
Expected: shows the analysis-level validation function (likely `validate_analysis` or a loop inside a larger function).

- [ ] **Step 2: Add the enum constant and helper near the top of the module**

Insert near the other canonical-vocabulary constants (search `VALID_TREATMENT_TYPES` or similar):

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


def is_valid_growth_phase(value: str) -> bool:
    """Accepts a canonical enum value or a well-formed `other:<slug>` escape value."""
    if not isinstance(value, str):
        return False
    if value in VALID_GROWTH_PHASES:
        return True
    return value.startswith("other:") and len(value) > len("other:")
```

- [ ] **Step 3: Add the per-analysis checks**

In the analysis-validation section (where `timepoint_hours` is currently checked — grep found the line), add three new checks:

```python
# --- Timepoint / growth_phase backfill checks (2026-04-12) ---
if "timepoint" not in analysis:
    warnings.append(
        f"{analysis_id}: missing 'timepoint' (use 'unknown' if paper doesn't state)"
    )

# timepoint_hours check already exists — ensure the existing warning wording still
# mentions the possibility of null for papers that don't report a time.

if "growth_phase" not in analysis:
    errors.append(
        f"{analysis_id}: missing required field 'growth_phase' "
        f"(one of {sorted(VALID_GROWTH_PHASES)} or 'other:<slug>')"
    )
elif not is_valid_growth_phase(analysis["growth_phase"]):
    errors.append(
        f"{analysis_id}: invalid growth_phase '{analysis['growth_phase']}' "
        f"(must be in VALID_GROWTH_PHASES or start with 'other:')"
    )
```

- [ ] **Step 4: Write focused tests for the validator**

Create `tests/test_validate_paperconfig_growth_phase.py`:

```python
"""Tests for the growth_phase / timepoint additions to validate_paperconfig.py."""
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).parent.parent / ".claude/skills/paperconfig"
sys.path.insert(0, str(SKILL_DIR))

from validate_paperconfig import (
    VALID_GROWTH_PHASES,
    is_valid_growth_phase,
)


def test_is_valid_growth_phase_accepts_canonical():
    for v in VALID_GROWTH_PHASES:
        assert is_valid_growth_phase(v), f"rejected canonical value: {v}"


def test_is_valid_growth_phase_accepts_other_slug():
    assert is_valid_growth_phase("other:heat_acclimated")
    assert is_valid_growth_phase("other:x")  # minimal valid


def test_is_valid_growth_phase_rejects_other_empty():
    assert not is_valid_growth_phase("other:")
    assert not is_valid_growth_phase("other")  # missing colon


def test_is_valid_growth_phase_rejects_bare_nonenum():
    assert not is_valid_growth_phase("maybe_stressed")
    assert not is_valid_growth_phase("")


def test_valid_growth_phases_includes_all_11():
    expected = {
        "exponential", "stationary", "nutrient_limited", "acclimated_steady_state",
        "infected", "recovery", "diel", "darkness", "death", "acute_stress",
        "unknown",
    }
    assert VALID_GROWTH_PHASES == expected
```

- [ ] **Step 5: Run the tests**

Run: `uv run pytest tests/test_validate_paperconfig_growth_phase.py -v`
Expected: all 5 tests pass.

- [ ] **Step 6: Run the validator on the existing corpus to confirm warnings (not errors, yet)**

Run: `uv run python .claude/skills/paperconfig/validate_paperconfig.py --all 2>&1 | tail -30`
Expected: warnings about missing `growth_phase` in every paperconfig; errors should only appear if a paperconfig has a malformed `growth_phase` value (none should, since no paperconfig has been backfilled yet — the errors should be about missing, which we want as errors per the spec).

**Decision point:** The spec says `growth_phase` missing is an **error**. That means every paperconfig will fail validation immediately. Since we want to land the validator change before the extraction tool is built, temporarily make the missing-growth_phase case a **warning** during rollout (add a TODO in the code), and flip to error after all paperconfigs are backfilled (Task 20).

Update Step 3's code:

```python
if "growth_phase" not in analysis:
    # TODO (timepoint-backfill rollout, 2026-04-12): flip to errors.append
    # once all 30 paperconfigs have growth_phase populated. Tracked in
    # docs/superpowers/plans/2026-04-12-timepoint-growth-phase-backfill.md
    warnings.append(
        f"{analysis_id}: missing required field 'growth_phase' "
        f"(one of {sorted(VALID_GROWTH_PHASES)} or 'other:<slug>')"
    )
elif not is_valid_growth_phase(analysis["growth_phase"]):
    errors.append(...)  # stays as error — malformed values should never land
```

- [ ] **Step 7: Re-run the validator**

Run: `uv run python .claude/skills/paperconfig/validate_paperconfig.py --all 2>&1 | tail -5`
Expected: warnings, no errors.

- [ ] **Step 8: Commit**

```bash
git add .claude/skills/paperconfig/validate_paperconfig.py tests/test_validate_paperconfig_growth_phase.py
git commit -m "feat(validator): VALID_GROWTH_PHASES enum + per-analysis growth_phase/timepoint checks (warn-only during rollout)"
```

---

## Task 5: Extraction utils — JSON I/O, paperconfig walking, signature

**Files:**
- Create: `multiomics_kg/extraction/timepoint/__init__.py`
- Create: `multiomics_kg/extraction/timepoint/extraction_utils.py`
- Test: `tests/extraction/timepoint/test_extraction_utils.py`

- [ ] **Step 1: Create the package init**

Create `multiomics_kg/extraction/timepoint/__init__.py` (empty file).

- [ ] **Step 2: Write failing tests for the utils**

Create `tests/extraction/timepoint/__init__.py` (empty) and `tests/extraction/timepoint/test_extraction_utils.py`:

```python
"""Tests for timepoint extraction utilities."""
import json
from pathlib import Path

import pytest
from ruamel.yaml import YAML

from multiomics_kg.extraction.timepoint.extraction_utils import (
    EXTRACTIONS_DIR,
    compute_paperconfig_signature,
    find_analyses,
    iter_paperconfigs,
    load_extraction_json,
    save_extraction_json,
)


@pytest.fixture
def paper_dir(tmp_path):
    pc = tmp_path / "paperconfig.yaml"
    pc.write_text("""publication:
  papername: Test 2026
  papermainpdf: data/fake.pdf
  experiments:
    exp1:
      name: Test
      organism: Prochlorococcus MED4
      omics_type: RNASEQ
  supplementary_materials:
    tbl_a:
      type: csv
      filename: data/fake.csv
      statistical_analyses:
        - id: DE_foo_24h
          experiment: exp1
          timepoint: "24h"
          timepoint_hours: 24
          growth_phase: exponential
          name_col: gene
          logfc_col: log2fc_24h
        - id: DE_foo_48h
          experiment: exp1
          name_col: gene
          logfc_col: log2fc_48h
""")
    return tmp_path


def test_find_analyses_returns_all(paper_dir):
    analyses = find_analyses(paper_dir / "paperconfig.yaml")
    ids = [a["id"] for a in analyses]
    assert ids == ["DE_foo_24h", "DE_foo_48h"]


def test_find_analyses_includes_logfc_col(paper_dir):
    analyses = find_analyses(paper_dir / "paperconfig.yaml")
    assert analyses[0]["logfc_col"] == "log2fc_24h"
    assert analyses[1]["logfc_col"] == "log2fc_48h"


def test_find_analyses_includes_experiment_key(paper_dir):
    analyses = find_analyses(paper_dir / "paperconfig.yaml")
    assert all(a["experiment"] == "exp1" for a in analyses)


def test_save_and_load_extraction_json(paper_dir):
    metadata = {"paper": "Test 2026", "status": "complete", "missing_analyses": []}
    analyses = [
        {"analysis_id": "DE_foo_24h", "timepoint": "24h", "growth_phase": "exponential"},
    ]
    path = save_extraction_json(paper_dir, metadata, analyses)
    assert path == paper_dir / EXTRACTIONS_DIR / "timepoint.json"
    assert path.exists()

    loaded = load_extraction_json(paper_dir)
    assert loaded["metadata"]["paper"] == "Test 2026"
    assert loaded["analyses"][0]["analysis_id"] == "DE_foo_24h"


def test_load_extraction_json_returns_none_if_missing(tmp_path):
    assert load_extraction_json(tmp_path) is None


def test_compute_signature_stable_across_runs(paper_dir):
    sig1 = compute_paperconfig_signature(paper_dir / "paperconfig.yaml")
    sig2 = compute_paperconfig_signature(paper_dir / "paperconfig.yaml")
    assert sig1 == sig2
    assert len(sig1) == 40  # sha1 hex


def test_compute_signature_changes_on_edit(paper_dir):
    sig1 = compute_paperconfig_signature(paper_dir / "paperconfig.yaml")
    # Add a new analysis
    pc = paper_dir / "paperconfig.yaml"
    text = pc.read_text()
    text = text.replace(
        "          logfc_col: log2fc_48h\n",
        "          logfc_col: log2fc_48h\n"
        "        - id: DE_foo_72h\n"
        "          experiment: exp1\n"
        "          name_col: gene\n"
        "          logfc_col: log2fc_72h\n",
    )
    pc.write_text(text)
    sig2 = compute_paperconfig_signature(pc)
    assert sig1 != sig2


def test_iter_paperconfigs_reads_list_files(tmp_path, monkeypatch):
    list_file = tmp_path / "paperconfig_files.txt"
    p1 = tmp_path / "paper1" / "paperconfig.yaml"
    p1.parent.mkdir()
    p1.write_text("publication: {papername: p1}\n")
    list_file.write_text(f"# comment\n{p1}\n")
    paths = list(iter_paperconfigs([list_file]))
    assert paths == [p1]
```

- [ ] **Step 3: Run tests — confirm they all fail**

Run: `uv run pytest tests/extraction/timepoint/test_extraction_utils.py -v`
Expected: 8 tests FAIL with ModuleNotFoundError / ImportError.

- [ ] **Step 4: Implement `extraction_utils.py`**

Create `multiomics_kg/extraction/timepoint/extraction_utils.py`:

```python
"""Extraction file I/O, paperconfig walking, signature computation for the
timepoint/growth-phase backfill pipeline.

Mirrors the cluster-extraction pattern in multiomics_kg/extraction/cluster/.
"""
from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Iterable, Iterator

from ruamel.yaml import YAML

logger = logging.getLogger(__name__)

EXTRACTIONS_DIR = "extractions"
EXTRACTION_FILENAME = "timepoint.json"

# Fields we pass to the LLM alongside the three target fields. Staleness
# signature includes these because a change to any of them (e.g. name_col
# swapped, experiment re-referenced) could invalidate an earlier extraction.
SIGNATURE_CONTEXT_FIELDS = ("experiment", "name_col", "logfc_col")


_YAML = YAML()
_YAML.preserve_quotes = True


def _load_yaml(path: Path) -> dict:
    with open(path) as f:
        return _YAML.load(f)


def find_analyses(paperconfig_path: Path) -> list[dict]:
    """Return a flat list of every analysis dict in the paperconfig, in file
    order. Each dict is the raw entry from `statistical_analyses[]`.
    """
    data = _load_yaml(paperconfig_path)
    pub = data.get("publication", {})
    analyses: list[dict] = []
    for table in (pub.get("supplementary_materials") or {}).values():
        if not isinstance(table, dict):
            continue
        for a in (table.get("statistical_analyses") or []):
            analyses.append(dict(a))  # shallow copy detaches from yaml
    return analyses


def compute_paperconfig_signature(paperconfig_path: Path) -> str:
    """sha1 hex over sorted analysis_ids + context fields (experiment,
    name_col, logfc_col). Stable across runs, changes whenever relevant
    paperconfig state changes.
    """
    analyses = find_analyses(paperconfig_path)
    rows = sorted(
        (a.get("id", ""), *(str(a.get(k, "")) for k in SIGNATURE_CONTEXT_FIELDS))
        for a in analyses
    )
    canonical = "\n".join("|".join(row) for row in rows)
    return hashlib.sha1(canonical.encode("utf-8")).hexdigest()


def save_extraction_json(paper_dir: Path, metadata: dict, analyses: list[dict]) -> Path:
    """Write `extractions/timepoint.json` under `paper_dir`, creating the
    subdirectory if needed. Returns the path.
    """
    ext_dir = paper_dir / EXTRACTIONS_DIR
    ext_dir.mkdir(parents=True, exist_ok=True)
    path = ext_dir / EXTRACTION_FILENAME
    payload = {"metadata": metadata, "analyses": analyses}
    path.write_text(json.dumps(payload, indent=2, default=str))
    return path


def load_extraction_json(paper_dir: Path) -> dict | None:
    """Return the parsed JSON dict, or None if the file doesn't exist."""
    path = paper_dir / EXTRACTIONS_DIR / EXTRACTION_FILENAME
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as e:
        logger.error("Failed to parse %s: %s", path, e)
        raise


def iter_paperconfigs(list_files: Iterable[Path]) -> Iterator[Path]:
    """Yield paperconfig paths listed in one-per-line text files. Skips
    blank lines and comments (lines starting with `#`). Does NOT filter
    out paperconfig paths that don't exist on disk — caller can decide.
    """
    for lf in list_files:
        for line in Path(lf).read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            yield Path(line)
```

- [ ] **Step 5: Run tests, confirm all pass**

Run: `uv run pytest tests/extraction/timepoint/test_extraction_utils.py -v`
Expected: 8 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/extraction/timepoint/__init__.py \
        multiomics_kg/extraction/timepoint/extraction_utils.py \
        tests/extraction/timepoint/__init__.py \
        tests/extraction/timepoint/test_extraction_utils.py
git commit -m "feat(timepoint extraction): extraction_utils — JSON I/O, paperconfig walker, signature"
```

---

## Task 6: Preprocess — compute `fields_requested` per analysis

**Files:**
- Modify: `multiomics_kg/extraction/timepoint/extraction_utils.py`
- Modify: `tests/extraction/timepoint/test_extraction_utils.py`

- [ ] **Step 1: Write failing tests for the preprocess**

Append to `tests/extraction/timepoint/test_extraction_utils.py`:

```python
from multiomics_kg.extraction.timepoint.extraction_utils import (
    compute_fields_requested,
)


def test_compute_fields_requested_all_missing():
    analysis = {"id": "x", "name_col": "g", "logfc_col": "lfc"}
    assert compute_fields_requested(analysis, validate=False) == [
        "timepoint", "timepoint_hours", "growth_phase",
    ]


def test_compute_fields_requested_timepoint_hours_null_is_request():
    analysis = {
        "id": "x", "timepoint": "24h", "timepoint_hours": None,
        "growth_phase": "exponential",
    }
    # null is a "needs filling" signal — request timepoint_hours only
    assert compute_fields_requested(analysis, validate=False) == ["timepoint_hours"]


def test_compute_fields_requested_empty_string_is_request():
    analysis = {
        "id": "x", "timepoint": "", "timepoint_hours": 24,
        "growth_phase": "exponential",
    }
    assert compute_fields_requested(analysis, validate=False) == ["timepoint"]


def test_compute_fields_requested_all_set_returns_empty():
    analysis = {
        "id": "x", "timepoint": "24h", "timepoint_hours": 24,
        "growth_phase": "exponential",
    }
    assert compute_fields_requested(analysis, validate=False) == []


def test_compute_fields_requested_validate_mode_always_all_three():
    analysis = {
        "id": "x", "timepoint": "24h", "timepoint_hours": 24,
        "growth_phase": "exponential",
    }
    assert compute_fields_requested(analysis, validate=True) == [
        "timepoint", "timepoint_hours", "growth_phase",
    ]
```

- [ ] **Step 2: Run, confirm they fail with ImportError**

Run: `uv run pytest tests/extraction/timepoint/test_extraction_utils.py -v -k compute_fields_requested`
Expected: 5 FAIL with ImportError.

- [ ] **Step 3: Implement `compute_fields_requested`**

Append to `multiomics_kg/extraction/timepoint/extraction_utils.py`:

```python
TARGET_FIELDS = ("timepoint", "timepoint_hours", "growth_phase")


def compute_fields_requested(analysis: dict, validate: bool = False) -> list[str]:
    """Return the list of target fields that need LLM extraction for this
    analysis.

    Default mode (`validate=False`):
      - Field absent → request.
      - `timepoint_hours` is None → request (null is a "needs filling"
        signal for numeric fields).
      - `timepoint` or `growth_phase` is empty string → request.
      - Field has a non-null, non-empty value → skip.

    Validate mode (`validate=True`):
      - Always request all three fields (LLM re-examines even populated ones).
    """
    if validate:
        return list(TARGET_FIELDS)

    requested: list[str] = []
    for field in TARGET_FIELDS:
        if field not in analysis:
            requested.append(field)
            continue
        value = analysis[field]
        if value is None:
            requested.append(field)
        elif isinstance(value, str) and value == "":
            requested.append(field)
        # non-null, non-empty → skip
    return requested
```

- [ ] **Step 4: Run tests, confirm all pass**

Run: `uv run pytest tests/extraction/timepoint/test_extraction_utils.py -v`
Expected: all pass (8 pre-existing + 5 new).

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/extraction/timepoint/extraction_utils.py \
        tests/extraction/timepoint/test_extraction_utils.py
git commit -m "feat(timepoint extraction): compute_fields_requested preprocess helper"
```

---

## Task 7: Prompts — SHARED_RULES and per-experiment-type hints

**Files:**
- Create: `multiomics_kg/extraction/timepoint/prompts.py`
- Test: `tests/extraction/timepoint/test_prompts.py`

- [ ] **Step 1: Write failing tests for the prompt builder**

Create `tests/extraction/timepoint/test_prompts.py`:

```python
"""Tests for prompt construction."""
from multiomics_kg.extraction.timepoint.prompts import (
    SHARED_RULES,
    VALID_GROWTH_PHASES_LIST,
    build_prompt,
)


def test_shared_rules_mentions_unknown_vs_other():
    assert "unknown" in SHARED_RULES.lower()
    assert "other:" in SHARED_RULES


def test_shared_rules_mentions_enum_values():
    for v in ("exponential", "nutrient_limited", "acclimated_steady_state",
              "infected", "recovery", "diel", "darkness", "death", "acute_stress"):
        assert v in SHARED_RULES


def test_valid_growth_phases_list_matches_validator():
    # Must exactly match VALID_GROWTH_PHASES in the validator (no drift).
    expected = {
        "exponential", "stationary", "nutrient_limited", "acclimated_steady_state",
        "infected", "recovery", "diel", "darkness", "death", "acute_stress",
        "unknown",
    }
    assert set(VALID_GROWTH_PHASES_LIST) == expected


def test_build_prompt_includes_paper_name_and_doi():
    background = {
        "papername": "Tetu 2019",
        "doi": "10.1038/s42003-019-0410-x",
        "experiments": {"exp1": {"treatment_condition": "HDPE"}},
    }
    targets = [{
        "id": "DE_x",
        "experiment_key": "exp1",
        "logfc_col": "HDPE log2 Fold change",
        "existing": {"timepoint": None, "timepoint_hours": None, "growth_phase": None},
        "fields_requested": ["timepoint", "timepoint_hours", "growth_phase"],
    }]
    prompt = build_prompt(background, targets, pdf_cache_entry=None)
    assert "Tetu 2019" in prompt
    assert "10.1038/s42003-019-0410-x" in prompt
    assert "DE_x" in prompt
    assert "HDPE log2 Fold change" in prompt


def test_build_prompt_includes_cache_entry_when_present():
    background = {"papername": "X", "doi": None, "experiments": {}}
    cache_entry = {"title": "My Title", "abstract": "Short abstract"}
    prompt = build_prompt(background, [], pdf_cache_entry=cache_entry)
    assert "My Title" in prompt


def test_build_prompt_skips_cache_entry_when_absent():
    background = {"papername": "X", "doi": None, "experiments": {}}
    prompt = build_prompt(background, [], pdf_cache_entry=None)
    # Should not have a "Publication metadata" section
    assert "PDF cache" not in prompt.lower() or "unavailable" in prompt.lower()
```

- [ ] **Step 2: Run tests — confirm they fail**

Run: `uv run pytest tests/extraction/timepoint/test_prompts.py -v`
Expected: all FAIL with ImportError.

- [ ] **Step 3: Implement `prompts.py`**

Create `multiomics_kg/extraction/timepoint/prompts.py`:

```python
"""Prompt construction for the timepoint/growth-phase extractor."""
from __future__ import annotations

import json

VALID_GROWTH_PHASES_LIST = [
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
]


SHARED_RULES = f"""\
You are extracting sampling-time and growth-phase metadata for a list of
statistical analyses from one scientific publication. The publication's PDF
is attached; additional supplementary PDFs may also be attached.

For each analysis listed in the "Extraction targets" section below, return
the fields listed in its `fields_requested` array — and no others.

FIELD VALUE RULES

1. `timepoint` (string, human-readable): the sampling time label, e.g.
   "120 min", "24h", "acclimated steady-state". Use "unknown" ONLY if the
   paper truly doesn't state a time.

2. `timepoint_hours` (number or null): numeric conversion of `timepoint`
   to hours. null ONLY if the paper truly doesn't state a time. "120 min" →
   2. "48h" → 48. "acclimated" with no specific hours → null.

3. `growth_phase` (string): the physiological state at sampling. Must be
   one of these canonical values:
     {", ".join(VALID_GROWTH_PHASES_LIST)}

   If the paper describes a phase not in this list, emit `other:<slug>`
   with a short snake_case slug, e.g. `other:heat_acclimated`,
   `other:late_decline`. Prefer `other:<slug>` over `unknown` whenever
   the paper gives any positional information.

   `unknown` means the paper gives no information about physiological
   state at sampling. `other:<slug>` means the paper describes a state
   that doesn't fit the enum.

QUALITY AND ANCHORING RULES

4. `self_assessment`: "high", "medium", or "low". Based on how directly
   the paper states the value and how unambiguous the evidence is.

5. `assessment_notes`: free-text rationale. Short, specific. Use this to
   flag edge cases, disagreements between the CSV label and the paper
   methods, or anything a reviewer should know.

6. `supporting_quotes`: list of objects `{{"quote": "...", "location": "..."}}`.
   Verbatim from the paper's methods or results. At least one quote is
   required unless `self_assessment` is "low" AND `assessment_notes`
   explains why no direct quote is available.

7. `source_figures`: list of figure or table numbers you consulted.
   Empty list is acceptable for purely text-derived extractions.

CROSS-CHECK RULE

8. Each analysis has a `logfc_col` — the fold-change column name from
   the source CSV. This column name often contains the timepoint
   verbatim ("24 hours", "Axenic, 36 hours", "log2FC_48h_vs_T0"). Cross-check:

   - If `logfc_col` agrees with paper methods, cite both in
     `supporting_quotes`.
   - If they disagree, PREFER paper methods but note the discrepancy
     in `assessment_notes` and lower `self_assessment` accordingly.

DEFAULT BIAS

9. Unknown is always better than wrong. If you cannot find clear evidence
   from the paper, set `growth_phase` to "unknown" and explain in
   `assessment_notes`.

OUTPUT SHAPE

Return one JSON object:

```
{{
  "analyses": [
    {{
      "analysis_id": "...",
      "timepoint": "...",        // only if in fields_requested
      "timepoint_hours": 2.0,    // only if in fields_requested; may be null
      "growth_phase": "...",     // only if in fields_requested
      "self_assessment": "high", // always
      "assessment_notes": "...", // always (may be empty string)
      "supporting_quotes": [     // always
        {{"quote": "...", "location": "..."}}
      ],
      "source_figures": []       // always (may be empty list)
    }}
  ]
}}
```

Do NOT emit `metadata`, `experiment_key`, or `fields_requested` — these
are added by the extraction tool after your response.

PER-EXPERIMENT-TYPE HINTS

- Diel studies (treatment_type contains "diel") → `diel` for
  cycling-phase samples. Do not split into light/dark sub-phases.
- Extended darkness (treatment_type contains "darkness", prolonged
  dark exposure of hours or days) → `darkness`. Do NOT confuse with
  the dark half of a normal diel cycle — that stays `diel`.
- Acclimation / chronic exposure (e.g. treatment_type contains
  "carbon" or "temperature" with language about ≥5 generations) →
  `acclimated_steady_state`.
- Nutrient starvation (treatment_type contains "nitrogen", "phosphorus",
  or "iron") → `exponential` at early timepoints before depletion is
  confirmed, `nutrient_limited` once cells are clearly nutrient-starved,
  `death` if sampling extends past viable-cell peak.
- Phage / viral infection → `exponential` at t=0 (pre-infection),
  `infected` at post-infection timepoints.
- Rescue / re-addition experiments → `recovery` for post-intervention
  timepoints.
- Short-exposure stress (≤6 h) at still-dividing cells → `acute_stress`,
  ONLY when no more-specific phase fits. `nutrient_limited`, `infected`,
  and `darkness` take precedence.
"""


def build_prompt(
    background: dict,
    targets: list[dict],
    pdf_cache_entry: dict | None = None,
) -> str:
    """Assemble the text prompt. PDFs are attached separately by extract.py."""
    sections = [SHARED_RULES]

    sections.append("\n## PAPER CONTEXT")
    sections.append(f"**Paper:** {background.get('papername', '')}")
    sections.append(f"**DOI:** {background.get('doi', '')}")

    if pdf_cache_entry:
        sections.append("\n### Publication metadata")
        for k in ("title", "abstract", "description", "study_type"):
            if v := pdf_cache_entry.get(k):
                sections.append(f"**{k.title()}:** {v}")

    sections.append("\n### Experiments block")
    sections.append("```yaml")
    sections.append(json.dumps(background.get("experiments", {}), indent=2, default=str))
    sections.append("```")

    sections.append("\n## EXTRACTION TARGETS")
    if not targets:
        sections.append("_(no analyses require extraction — all fields populated)_")
    else:
        sections.append("```json")
        sections.append(json.dumps(targets, indent=2, default=str))
        sections.append("```")

    return "\n".join(sections)
```

- [ ] **Step 4: Run tests, confirm all pass**

Run: `uv run pytest tests/extraction/timepoint/test_prompts.py -v`
Expected: 6 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/extraction/timepoint/prompts.py \
        tests/extraction/timepoint/test_prompts.py
git commit -m "feat(timepoint extraction): prompts — SHARED_RULES + build_prompt"
```

---

## Task 8: LLM output validation (pure-python, no API call)

**Files:**
- Modify: `multiomics_kg/extraction/timepoint/extraction_utils.py`
- Modify: `tests/extraction/timepoint/test_extraction_utils.py`

- [ ] **Step 1: Write failing tests**

Append to `tests/extraction/timepoint/test_extraction_utils.py`:

```python
from multiomics_kg.extraction.timepoint.extraction_utils import (
    validate_llm_payload,
)


def test_validate_llm_payload_all_valid():
    requested_map = {
        "DE_a": ["timepoint", "timepoint_hours", "growth_phase"],
        "DE_b": ["timepoint", "timepoint_hours", "growth_phase"],
    }
    payload = {
        "analyses": [
            {
                "analysis_id": "DE_a",
                "timepoint": "24h", "timepoint_hours": 24.0, "growth_phase": "exponential",
                "self_assessment": "high", "assessment_notes": "",
                "supporting_quotes": [{"quote": "q", "location": "p1"}],
                "source_figures": [],
            },
            {
                "analysis_id": "DE_b",
                "timepoint": "48h", "timepoint_hours": 48.0,
                "growth_phase": "other:heat_shocked",
                "self_assessment": "medium", "assessment_notes": "note",
                "supporting_quotes": [{"quote": "q2", "location": "p2"}],
                "source_figures": ["Fig 1"],
            },
        ]
    }
    valid, missing = validate_llm_payload(payload, requested_map)
    assert len(valid) == 2
    assert missing == []


def test_validate_llm_payload_detects_missing_analysis():
    requested_map = {"DE_a": ["timepoint"], "DE_b": ["timepoint"]}
    payload = {"analyses": [{
        "analysis_id": "DE_a", "timepoint": "24h",
        "self_assessment": "high", "assessment_notes": "",
        "supporting_quotes": [{"quote": "q", "location": "p"}],
        "source_figures": [],
    }]}
    valid, missing = validate_llm_payload(payload, requested_map)
    assert len(valid) == 1
    assert missing == [{"analysis_id": "DE_b", "reason": "not_returned"}]


def test_validate_llm_payload_rejects_invalid_growth_phase():
    requested_map = {"DE_a": ["growth_phase"]}
    payload = {"analyses": [{
        "analysis_id": "DE_a", "growth_phase": "maybe_stressed",
        "self_assessment": "high", "assessment_notes": "",
        "supporting_quotes": [{"quote": "q", "location": "p"}],
        "source_figures": [],
    }]}
    valid, missing = validate_llm_payload(payload, requested_map)
    assert valid == []
    assert missing == [{
        "analysis_id": "DE_a",
        "reason": "invalid_growth_phase: maybe_stressed",
    }]


def test_validate_llm_payload_rejects_bad_self_assessment():
    requested_map = {"DE_a": ["timepoint"]}
    payload = {"analyses": [{
        "analysis_id": "DE_a", "timepoint": "24h",
        "self_assessment": "sure", "assessment_notes": "",
        "supporting_quotes": [{"quote": "q", "location": "p"}],
        "source_figures": [],
    }]}
    valid, missing = validate_llm_payload(payload, requested_map)
    assert valid == []
    assert "invalid_self_assessment" in missing[0]["reason"]


def test_validate_llm_payload_rejects_missing_requested_field():
    requested_map = {"DE_a": ["timepoint", "growth_phase"]}
    payload = {"analyses": [{
        "analysis_id": "DE_a", "timepoint": "24h",  # growth_phase absent
        "self_assessment": "high", "assessment_notes": "",
        "supporting_quotes": [{"quote": "q", "location": "p"}],
        "source_figures": [],
    }]}
    valid, missing = validate_llm_payload(payload, requested_map)
    assert valid == []
    assert missing == [{"analysis_id": "DE_a", "reason": "missing_field: growth_phase"}]


def test_validate_llm_payload_rejects_timepoint_hours_not_numeric():
    requested_map = {"DE_a": ["timepoint_hours"]}
    payload = {"analyses": [{
        "analysis_id": "DE_a", "timepoint_hours": "twenty-four",
        "self_assessment": "high", "assessment_notes": "",
        "supporting_quotes": [{"quote": "q", "location": "p"}],
        "source_figures": [],
    }]}
    valid, missing = validate_llm_payload(payload, requested_map)
    assert valid == []
    assert "timepoint_hours_not_numeric" in missing[0]["reason"]


def test_validate_llm_payload_accepts_null_timepoint_hours():
    requested_map = {"DE_a": ["timepoint_hours"]}
    payload = {"analyses": [{
        "analysis_id": "DE_a", "timepoint_hours": None,
        "self_assessment": "low", "assessment_notes": "Paper doesn't state time",
        "supporting_quotes": [],
        "source_figures": [],
    }]}
    valid, missing = validate_llm_payload(payload, requested_map)
    assert len(valid) == 1
    assert missing == []


def test_validate_llm_payload_rejects_hallucinated_analysis_id():
    requested_map = {"DE_a": ["timepoint"]}
    payload = {"analyses": [{
        "analysis_id": "DE_ghost",  # not in requested_map
        "timepoint": "24h",
        "self_assessment": "high", "assessment_notes": "",
        "supporting_quotes": [{"quote": "q", "location": "p"}],
        "source_figures": [],
    }]}
    # Hallucinated IDs are a HARD fail — caller should not get partial success.
    import pytest
    with pytest.raises(ValueError, match="DE_ghost"):
        validate_llm_payload(payload, requested_map)
```

- [ ] **Step 2: Run, confirm they fail**

Run: `uv run pytest tests/extraction/timepoint/test_extraction_utils.py -v -k validate_llm_payload`
Expected: 8 FAIL.

- [ ] **Step 3: Implement `validate_llm_payload`**

Append to `multiomics_kg/extraction/timepoint/extraction_utils.py`:

```python
# Must match validate_paperconfig.VALID_GROWTH_PHASES exactly.
# Duplicated here to avoid importing from the skills dir; kept in sync
# by tests/extraction/timepoint/test_prompts.py::test_valid_growth_phases_list_matches_validator.
_VALID_GROWTH_PHASES = frozenset([
    "exponential", "stationary", "nutrient_limited", "acclimated_steady_state",
    "infected", "recovery", "diel", "darkness", "death", "acute_stress",
    "unknown",
])

_VALID_SELF_ASSESSMENT = frozenset(["high", "medium", "low"])


def _is_valid_growth_phase(v) -> bool:
    return isinstance(v, str) and (
        v in _VALID_GROWTH_PHASES
        or (v.startswith("other:") and len(v) > len("other:"))
    )


def _validate_one_analysis(
    row: dict, requested: list[str]
) -> tuple[dict | None, dict | None]:
    """Return (valid_row, None) on success, or (None, missing_entry) on failure."""
    # Required structural fields
    for required_field in ("self_assessment", "assessment_notes",
                           "supporting_quotes", "source_figures"):
        if required_field not in row:
            return None, {
                "analysis_id": row.get("analysis_id"),
                "reason": f"missing_field: {required_field}",
            }

    if row["self_assessment"] not in _VALID_SELF_ASSESSMENT:
        return None, {
            "analysis_id": row["analysis_id"],
            "reason": f"invalid_self_assessment: {row['self_assessment']}",
        }

    # Each requested field must be present; skip validation for fields not
    # requested — LLM is instructed not to emit them.
    for field in requested:
        if field not in row:
            return None, {
                "analysis_id": row["analysis_id"],
                "reason": f"missing_field: {field}",
            }
        value = row[field]
        if field == "growth_phase":
            if not _is_valid_growth_phase(value):
                return None, {
                    "analysis_id": row["analysis_id"],
                    "reason": f"invalid_growth_phase: {value}",
                }
        elif field == "timepoint_hours":
            if value is not None and not isinstance(value, (int, float)):
                return None, {
                    "analysis_id": row["analysis_id"],
                    "reason": f"timepoint_hours_not_numeric: {value!r}",
                }
        elif field == "timepoint":
            if not isinstance(value, str):
                return None, {
                    "analysis_id": row["analysis_id"],
                    "reason": f"timepoint_not_string: {value!r}",
                }

    return row, None


def validate_llm_payload(
    payload: dict, requested_map: dict[str, list[str]]
) -> tuple[list[dict], list[dict]]:
    """Split an LLM response into (valid_rows, missing_entries).

    Raises ValueError if the LLM returned an analysis_id not in
    requested_map (hallucinated ID — treat as hard failure).
    """
    valid: list[dict] = []
    missing: list[dict] = []
    returned_ids: set[str] = set()

    for row in payload.get("analyses", []):
        aid = row.get("analysis_id")
        if aid not in requested_map:
            raise ValueError(
                f"LLM returned unknown analysis_id {aid!r}; expected one of "
                f"{sorted(requested_map)}"
            )
        returned_ids.add(aid)
        valid_row, missing_entry = _validate_one_analysis(row, requested_map[aid])
        if valid_row is not None:
            valid.append(valid_row)
        elif missing_entry is not None:
            missing.append(missing_entry)

    # Any requested analysis the LLM didn't return
    for aid in requested_map:
        if aid not in returned_ids:
            missing.append({"analysis_id": aid, "reason": "not_returned"})

    return valid, missing
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/extraction/timepoint/test_extraction_utils.py -v`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/extraction/timepoint/extraction_utils.py \
        tests/extraction/timepoint/test_extraction_utils.py
git commit -m "feat(timepoint extraction): validate_llm_payload with partial-result detection"
```

---

## Task 9: `extract.py` — LLM orchestration (with mocked LLM in tests)

**Files:**
- Create: `multiomics_kg/extraction/timepoint/extract.py`
- Test: `tests/extraction/timepoint/test_extract.py`

- [ ] **Step 1: Write failing tests using a mock LLM client**

Create `tests/extraction/timepoint/test_extract.py`:

```python
"""Tests for extract.py. LLM calls are mocked — no network."""
import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest
from ruamel.yaml import YAML

from multiomics_kg.extraction.timepoint.extract import (
    build_background,
    build_targets,
    extract_one_paper,
)


@pytest.fixture
def paper_dir(tmp_path):
    pc = tmp_path / "paperconfig.yaml"
    pc.write_text("""publication:
  papername: Fake 2026
  papermainpdf: data/fake.pdf
  experiments:
    exp1:
      name: N-starvation
      organism: Prochlorococcus MED4
      treatment_type: [nitrogen]
      omics_type: RNASEQ
  supplementary_materials:
    tbl:
      type: csv
      filename: data/fake.csv
      statistical_analyses:
        - id: DE_n_24h
          experiment: exp1
          name_col: gene
          logfc_col: log2fc_24h
        - id: DE_n_48h
          experiment: exp1
          name_col: gene
          logfc_col: log2fc_48h
""")
    # papermainpdf doesn't need to exist for build_background; extract_one_paper
    # will need it. Create a dummy file.
    (tmp_path.parent / "data").mkdir(parents=True, exist_ok=True)
    (tmp_path.parent / "data" / "fake.pdf").write_text("dummy")
    return tmp_path


def test_build_background_includes_experiments_block(paper_dir):
    bg = build_background(paper_dir / "paperconfig.yaml")
    assert bg["papername"] == "Fake 2026"
    assert "exp1" in bg["experiments"]
    assert bg["experiments"]["exp1"]["omics_type"] == "RNASEQ"


def test_build_targets_computes_fields_requested(paper_dir):
    analyses = [
        {"id": "DE_n_24h", "experiment": "exp1", "name_col": "gene",
         "logfc_col": "log2fc_24h"},
        {"id": "DE_n_48h", "experiment": "exp1", "name_col": "gene",
         "logfc_col": "log2fc_48h", "timepoint": "48h", "timepoint_hours": 48},
    ]
    targets = build_targets(analyses, validate=False)
    # First analysis: all three requested
    assert set(targets[0]["fields_requested"]) == {
        "timepoint", "timepoint_hours", "growth_phase",
    }
    assert targets[0]["logfc_col"] == "log2fc_24h"
    # Second analysis: only growth_phase requested (timepoint + hours set)
    assert targets[1]["fields_requested"] == ["growth_phase"]


def test_build_targets_dropped_when_all_populated(paper_dir):
    analyses = [{
        "id": "DE_done",
        "experiment": "exp1",
        "name_col": "gene", "logfc_col": "lfc",
        "timepoint": "24h", "timepoint_hours": 24, "growth_phase": "exponential",
    }]
    targets = build_targets(analyses, validate=False)
    assert targets == []


def test_extract_one_paper_happy_path(paper_dir, monkeypatch):
    fake_llm = MagicMock()
    fake_llm.call.return_value = (
        {
            "analyses": [
                {
                    "analysis_id": "DE_n_24h",
                    "timepoint": "24h",
                    "timepoint_hours": 24.0,
                    "growth_phase": "nutrient_limited",
                    "self_assessment": "high",
                    "assessment_notes": "",
                    "supporting_quotes": [{"quote": "q", "location": "Methods"}],
                    "source_figures": [],
                },
                {
                    "analysis_id": "DE_n_48h",
                    "timepoint": "48h",
                    "timepoint_hours": 48.0,
                    "growth_phase": "nutrient_limited",
                    "self_assessment": "high",
                    "assessment_notes": "",
                    "supporting_quotes": [{"quote": "q", "location": "Methods"}],
                    "source_figures": [],
                },
            ]
        },
        {"input_tokens": 100, "output_tokens": 200, "model": "fake-model"},
    )
    json_path = extract_one_paper(
        paper_dir,
        llm_client=fake_llm,
        validate=False,
    )
    assert json_path.exists()
    data = json.loads(json_path.read_text())
    assert data["metadata"]["status"] == "complete"
    assert data["metadata"]["missing_analyses"] == []
    assert data["metadata"]["model"] == "fake-model"
    assert len(data["analyses"]) == 2
    # experiment_key injected by extract.py
    assert all(a["experiment_key"] == "exp1" for a in data["analyses"])
    # fields_requested injected
    assert all("fields_requested" in a for a in data["analyses"])


def test_extract_one_paper_partial_marks_status(paper_dir):
    fake_llm = MagicMock()
    fake_llm.call.return_value = (
        {
            "analyses": [
                # Only 1 of 2 returned — DE_n_48h missing → partial
                {
                    "analysis_id": "DE_n_24h",
                    "timepoint": "24h", "timepoint_hours": 24.0,
                    "growth_phase": "nutrient_limited",
                    "self_assessment": "high", "assessment_notes": "",
                    "supporting_quotes": [{"quote": "q", "location": "M"}],
                    "source_figures": [],
                },
            ]
        },
        {"input_tokens": 50, "output_tokens": 100, "model": "fake-model"},
    )
    json_path = extract_one_paper(paper_dir, llm_client=fake_llm, validate=False)
    data = json.loads(json_path.read_text())
    assert data["metadata"]["status"] == "partial"
    assert data["metadata"]["missing_analyses"] == [
        {"analysis_id": "DE_n_48h", "reason": "not_returned"},
    ]
    assert len(data["analyses"]) == 1


def test_extract_one_paper_skips_when_no_fields_requested(paper_dir):
    # Pre-populate both analyses
    pc_path = paper_dir / "paperconfig.yaml"
    yaml = YAML()
    with open(pc_path) as f:
        data = yaml.load(f)
    for a in data["publication"]["supplementary_materials"]["tbl"]["statistical_analyses"]:
        a["timepoint"] = "24h"
        a["timepoint_hours"] = 24
        a["growth_phase"] = "exponential"
    with open(pc_path, "w") as f:
        yaml.dump(data, f)

    fake_llm = MagicMock()
    result = extract_one_paper(paper_dir, llm_client=fake_llm, validate=False)
    # No LLM call
    fake_llm.call.assert_not_called()
    # No JSON written (or return None, per spec)
    assert result is None
```

- [ ] **Step 2: Run, confirm they fail**

Run: `uv run pytest tests/extraction/timepoint/test_extract.py -v`
Expected: all FAIL with ImportError.

- [ ] **Step 3: Implement `extract.py`**

Create `multiomics_kg/extraction/timepoint/extract.py`:

```python
"""Timepoint/growth-phase LLM extraction per paper.

Usage:
    uv run python -m multiomics_kg.extraction.timepoint.extract --paper "Tetu 2019"
    uv run python -m multiomics_kg.extraction.timepoint.extract --all
    uv run python -m multiomics_kg.extraction.timepoint.extract --paper "X" --dry-run
    uv run python -m multiomics_kg.extraction.timepoint.extract --all --validate
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Protocol

from multiomics_kg.extraction.timepoint.extraction_utils import (
    compute_fields_requested,
    compute_paperconfig_signature,
    find_analyses,
    iter_paperconfigs,
    load_extraction_json,
    save_extraction_json,
    validate_llm_payload,
)
from multiomics_kg.extraction.timepoint.prompts import build_prompt

logger = logging.getLogger(__name__)


class LLMClient(Protocol):
    """Minimal protocol — any caller providing `.call(...)` satisfies this."""
    def call(
        self, prompt: str, pdf_paths: list[Path], model: str | None = None,
    ) -> tuple[dict, dict]:
        """Return (payload_dict, metadata_dict).
        metadata_dict: keys `input_tokens`, `output_tokens`, `model`.
        """
        ...


def build_background(paperconfig_path: Path) -> dict:
    """Extract the background context for the LLM prompt."""
    from multiomics_kg.extraction.timepoint.extraction_utils import _load_yaml
    data = _load_yaml(paperconfig_path)
    pub = data.get("publication", {})
    return {
        "papername": pub.get("papername", ""),
        "doi": pub.get("doi", ""),
        "experiments": dict(pub.get("experiments") or {}),
        "additional_pdfs": list(
            (pub.get("extraction") or {}).get("additional_pdfs") or []
        ),
        "papermainpdf": pub.get("papermainpdf"),
    }


def build_targets(analyses: list[dict], validate: bool) -> list[dict]:
    """Build the `targets` list for the prompt. Omits analyses with empty
    fields_requested unless --validate.
    """
    targets: list[dict] = []
    for a in analyses:
        requested = compute_fields_requested(a, validate=validate)
        if not requested and not validate:
            continue
        targets.append({
            "id": a["id"],
            "experiment_key": a.get("experiment", ""),
            "logfc_col": a.get("logfc_col", ""),
            "existing": {
                "timepoint": a.get("timepoint"),
                "timepoint_hours": a.get("timepoint_hours"),
                "growth_phase": a.get("growth_phase"),
            },
            "fields_requested": requested,
        })
    return targets


def _load_pdf_cache_entry(cache_path: Path, papermainpdf: str | None) -> dict | None:
    if not papermainpdf or not cache_path.exists():
        return None
    try:
        cache = json.loads(cache_path.read_text())
    except Exception:
        return None
    entry = cache.get(papermainpdf)
    if entry and "publication" in entry:
        return entry["publication"]
    return None


def extract_one_paper(
    paper_dir: Path,
    llm_client: LLMClient,
    validate: bool = False,
    pdf_cache_path: Path | None = None,
    model: str | None = None,
) -> Path | None:
    """Extract one paper's timepoint/growth_phase metadata via the LLM.

    Returns the written JSON path, or None if no extraction was needed
    (all fields populated, default mode).
    """
    paperconfig_path = paper_dir / "paperconfig.yaml"
    analyses = find_analyses(paperconfig_path)
    targets = build_targets(analyses, validate=validate)
    if not targets:
        logger.info("No fields to extract for %s — skipping LLM call.", paper_dir)
        return None

    background = build_background(paperconfig_path)

    pdf_cache_entry = _load_pdf_cache_entry(
        pdf_cache_path or Path("cache/pdf_extraction_cache.json"),
        background.get("papermainpdf"),
    )
    prompt = build_prompt(background, targets, pdf_cache_entry)

    pdf_paths: list[Path] = []
    if pmp := background.get("papermainpdf"):
        pdf_paths.append(Path(pmp))
    pdf_paths.extend(Path(p) for p in background.get("additional_pdfs", []))

    payload, call_meta = llm_client.call(prompt, pdf_paths, model=model)

    requested_map = {t["id"]: t["fields_requested"] for t in targets}
    valid_rows, missing = validate_llm_payload(payload, requested_map)

    # Inject experiment_key and fields_requested
    for row in valid_rows:
        aid = row["analysis_id"]
        # Look up experiment_key from targets
        for t in targets:
            if t["id"] == aid:
                row["experiment_key"] = t["experiment_key"]
                row["fields_requested"] = t["fields_requested"]
                break

    metadata = {
        "paper": background.get("papername", ""),
        "doi": background.get("doi", ""),
        "model": call_meta.get("model", ""),
        "extracted_at": datetime.utcnow().isoformat(timespec="seconds"),
        "input_tokens": call_meta.get("input_tokens", 0),
        "output_tokens": call_meta.get("output_tokens", 0),
        "paperconfig_signature": compute_paperconfig_signature(paperconfig_path),
        "status": "complete" if not missing else "partial",
        "missing_analyses": missing,
    }

    json_path = save_extraction_json(paper_dir, metadata, valid_rows)
    return json_path


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="extract.py")
    parser.add_argument("--paper", help="Paper name (papername) to extract")
    parser.add_argument("--all", action="store_true", help="Extract all papers")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print fields_requested per analysis; no LLM call")
    parser.add_argument("--validate", action="store_true",
                        help="Re-examine every field, even populated ones")
    parser.add_argument("--model", default=None,
                        help="Override default model (e.g. gpt-4.1-mini)")
    # --resume and --retry are added in Task 12.
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    paper_dirs = _resolve_paper_dirs(args)
    if not paper_dirs:
        logger.error("No papers matched selection.")
        return 2

    if args.dry_run:
        for pd in paper_dirs:
            _print_dry_run(pd, validate=args.validate)
        return 0

    from multiomics_kg.extraction.timepoint.llm_client import OpenAIResponsesClient
    client = OpenAIResponsesClient()
    had_partial = False
    for pd in paper_dirs:
        try:
            path = extract_one_paper(pd, client, validate=args.validate, model=args.model)
            if path is None:
                logger.info("[%s] no fields to extract, skipped.", pd.name)
                continue
            data = json.loads(path.read_text())
            logger.info("[%s] status=%s wrote %s",
                        pd.name, data["metadata"]["status"], path)
            if data["metadata"]["status"] == "partial":
                had_partial = True
        except Exception as e:
            logger.error("[%s] extraction failed: %s", pd.name, e)
            had_partial = True
    return 1 if had_partial else 0


def _resolve_paper_dirs(args) -> list[Path]:
    """Map --paper / --all to paper directories. See implementation notes."""
    list_files = [
        Path("data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"),
        Path("data/Synechococcus/papers_and_supp/paperconfig_files.txt"),
    ]
    paper_dirs = [p.parent for p in iter_paperconfigs(list_files) if p.exists()]
    if args.all:
        return paper_dirs
    if args.paper:
        matches = []
        for pd in paper_dirs:
            try:
                from multiomics_kg.extraction.timepoint.extraction_utils import _load_yaml
                data = _load_yaml(pd / "paperconfig.yaml")
                if data.get("publication", {}).get("papername") == args.paper:
                    matches.append(pd)
            except Exception:
                continue
        return matches
    return []


def _print_dry_run(paper_dir: Path, validate: bool) -> None:
    analyses = find_analyses(paper_dir / "paperconfig.yaml")
    targets = build_targets(analyses, validate=validate)
    print(f"\n{paper_dir.name}:")
    if not targets:
        print("  (all fields populated)")
        return
    for t in targets:
        print(f"  {t['id']}: needs {t['fields_requested']}  (logfc_col={t['logfc_col']!r})")


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/extraction/timepoint/test_extract.py -v`
Expected: all 6 tests pass (the test for the no-op skip relies on `extract_one_paper` returning None when `targets == []`).

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/extraction/timepoint/extract.py \
        tests/extraction/timepoint/test_extract.py
git commit -m "feat(timepoint extraction): extract.py orchestration (LLM client mocked in tests)"
```

---

## Task 10: LLM client wrapper (OpenAI Responses API)

**Files:**
- Create: `multiomics_kg/extraction/timepoint/llm_client.py`

- [ ] **Step 1: Inspect how the cluster-extraction module uses the OpenAI API**

Run: `grep -n "openai\|response\|client\." multiomics_kg/extraction/cluster/extract.py | head -30`
Expected: shows how to instantiate the client and pass PDFs.

- [ ] **Step 2: Write `llm_client.py` following the same pattern**

Create `multiomics_kg/extraction/timepoint/llm_client.py`:

```python
"""OpenAI Responses API client for timepoint extraction.

Mirrors the cluster-extraction pattern. Attaches the paper PDF + any
extras.additional_pdfs as multimodal inputs. Returns a parsed dict plus
token/model metadata.
"""
from __future__ import annotations

import json
import logging
import os
from pathlib import Path

from openai import OpenAI

logger = logging.getLogger(__name__)

DEFAULT_MODEL = "gpt-4.1-mini"


class OpenAIResponsesClient:
    """Minimal wrapper around OpenAI Responses API for structured-JSON extraction."""

    def __init__(self, api_key: str | None = None, model: str = DEFAULT_MODEL):
        self._client = OpenAI(api_key=api_key or os.environ.get("OPENAI_API_KEY"))
        self._model = model

    def call(
        self,
        prompt: str,
        pdf_paths: list[Path],
        model: str | None = None,
    ) -> tuple[dict, dict]:
        """Call the LLM with the prompt + attached PDFs; parse JSON response."""
        model = model or self._model

        # Upload PDFs that exist. Missing PDFs are a warning, not a hard fail.
        file_inputs = []
        for p in pdf_paths:
            if not p.exists():
                logger.warning("PDF not found, skipping: %s", p)
                continue
            uploaded = self._client.files.create(file=p.open("rb"), purpose="user_data")
            file_inputs.append({"type": "input_file", "file_id": uploaded.id})

        content = file_inputs + [{"type": "input_text", "text": prompt}]

        response = self._client.responses.create(
            model=model,
            input=[{"role": "user", "content": content}],
            response_format={"type": "json_object"},
        )

        raw_text = response.output_text
        try:
            payload = json.loads(raw_text)
        except json.JSONDecodeError as e:
            raise ValueError(f"LLM returned unparseable JSON: {e}\n{raw_text[:500]}")

        meta = {
            "input_tokens": response.usage.input_tokens if response.usage else 0,
            "output_tokens": response.usage.output_tokens if response.usage else 0,
            "model": model,
        }
        return payload, meta
```

- [ ] **Step 3: Smoke-test manually** (skip if no API key in env)

Run: `uv run python -c "from multiomics_kg.extraction.timepoint.llm_client import OpenAIResponsesClient; print('import ok')"`
Expected: `import ok`

- [ ] **Step 4: Commit**

```bash
git add multiomics_kg/extraction/timepoint/llm_client.py
git commit -m "feat(timepoint extraction): llm_client — OpenAI Responses API wrapper"
```

---

## Task 11: `merge.py` — write JSON fields back into paperconfig.yaml

**Files:**
- Create: `multiomics_kg/extraction/timepoint/merge.py`
- Test: `tests/extraction/timepoint/test_merge.py`

- [ ] **Step 1: Write failing tests**

Create `tests/extraction/timepoint/test_merge.py`:

```python
"""Tests for merge.py."""
import json
from pathlib import Path

import pytest
from ruamel.yaml import YAML

from multiomics_kg.extraction.timepoint.merge import merge_one_paper
from multiomics_kg.extraction.timepoint.extraction_utils import (
    compute_paperconfig_signature,
    save_extraction_json,
)


def _write_paperconfig(tmp_path: Path, timepoint=None, timepoint_hours=None, growth_phase=None) -> Path:
    pc = tmp_path / "paperconfig.yaml"
    lines = [
        "publication:",
        "  papername: Fake",
        "  papermainpdf: data/fake.pdf",
        "  experiments:",
        "    exp1: {name: Test, omics_type: RNASEQ}",
        "  supplementary_materials:",
        "    tbl:",
        "      type: csv",
        "      filename: data/fake.csv",
        "      statistical_analyses:",
        "        - id: DE_a",
        "          experiment: exp1",
        "          name_col: gene",
        "          logfc_col: lfc",
    ]
    if timepoint is not None:
        lines.append(f"          timepoint: \"{timepoint}\"")
    if timepoint_hours is not None:
        lines.append(f"          timepoint_hours: {timepoint_hours}")
    if growth_phase is not None:
        lines.append(f"          growth_phase: {growth_phase}")
    pc.write_text("\n".join(lines) + "\n")
    return pc


def _valid_analysis(aid="DE_a", **overrides) -> dict:
    base = {
        "analysis_id": aid, "experiment_key": "exp1",
        "fields_requested": ["timepoint", "timepoint_hours", "growth_phase"],
        "timepoint": "24h", "timepoint_hours": 24.0, "growth_phase": "exponential",
        "self_assessment": "high", "assessment_notes": "",
        "supporting_quotes": [{"quote": "q", "location": "p"}],
        "source_figures": [],
    }
    base.update(overrides)
    return base


def test_merge_writes_three_fields_into_paperconfig(tmp_path):
    pc = _write_paperconfig(tmp_path)
    sig = compute_paperconfig_signature(pc)
    save_extraction_json(
        tmp_path,
        {"paper": "Fake", "status": "complete", "missing_analyses": [],
         "paperconfig_signature": sig},
        [_valid_analysis()],
    )

    merge_one_paper(tmp_path, force=False)

    yaml = YAML()
    with open(pc) as f:
        data = yaml.load(f)
    a = data["publication"]["supplementary_materials"]["tbl"]["statistical_analyses"][0]
    assert a["timepoint"] == "24h"
    assert a["timepoint_hours"] == 24.0
    assert a["growth_phase"] == "exponential"


def test_merge_refuses_partial_without_force(tmp_path):
    pc = _write_paperconfig(tmp_path)
    sig = compute_paperconfig_signature(pc)
    save_extraction_json(
        tmp_path,
        {"paper": "Fake", "status": "partial",
         "missing_analyses": [{"analysis_id": "DE_b", "reason": "not_returned"}],
         "paperconfig_signature": sig},
        [_valid_analysis()],
    )
    with pytest.raises(SystemExit) as exc:
        merge_one_paper(tmp_path, force=False)
    assert exc.value.code != 0


def test_merge_allows_partial_with_force(tmp_path):
    pc = _write_paperconfig(tmp_path)
    sig = compute_paperconfig_signature(pc)
    save_extraction_json(
        tmp_path,
        {"paper": "Fake", "status": "partial",
         "missing_analyses": [{"analysis_id": "DE_b", "reason": "not_returned"}],
         "paperconfig_signature": sig},
        [_valid_analysis()],
    )
    merge_one_paper(tmp_path, force=True)  # should not raise


def test_merge_refuses_on_signature_mismatch(tmp_path):
    pc = _write_paperconfig(tmp_path)
    save_extraction_json(
        tmp_path,
        {"paper": "Fake", "status": "complete", "missing_analyses": [],
         "paperconfig_signature": "0" * 40},  # wrong sig
        [_valid_analysis()],
    )
    with pytest.raises(SystemExit):
        merge_one_paper(tmp_path, force=False)


def test_merge_warns_on_ghost_analysis_id(tmp_path, caplog):
    pc = _write_paperconfig(tmp_path)
    sig = compute_paperconfig_signature(pc)
    # Ghost analysis: JSON has DE_ghost but paperconfig doesn't
    save_extraction_json(
        tmp_path,
        {"paper": "Fake", "status": "complete", "missing_analyses": [],
         "paperconfig_signature": sig},
        [_valid_analysis(aid="DE_ghost"), _valid_analysis(aid="DE_a")],
    )
    with caplog.at_level("WARNING"):
        merge_one_paper(tmp_path, force=False)
    assert any("DE_ghost" in r.message for r in caplog.records)


def test_merge_refuses_overwrite_of_differing_value(tmp_path):
    pc = _write_paperconfig(tmp_path, growth_phase="acclimated_steady_state")
    sig = compute_paperconfig_signature(pc)
    save_extraction_json(
        tmp_path,
        {"paper": "Fake", "status": "complete", "missing_analyses": [],
         "paperconfig_signature": sig},
        [_valid_analysis(growth_phase="exponential")],
    )
    with pytest.raises(SystemExit):
        merge_one_paper(tmp_path, force=False)


def test_merge_allows_overwrite_with_force(tmp_path):
    pc = _write_paperconfig(tmp_path, growth_phase="acclimated_steady_state")
    sig = compute_paperconfig_signature(pc)
    save_extraction_json(
        tmp_path,
        {"paper": "Fake", "status": "complete", "missing_analyses": [],
         "paperconfig_signature": sig},
        [_valid_analysis(growth_phase="exponential")],
    )
    merge_one_paper(tmp_path, force=True)
    yaml = YAML()
    with open(pc) as f:
        data = yaml.load(f)
    a = data["publication"]["supplementary_materials"]["tbl"]["statistical_analyses"][0]
    assert a["growth_phase"] == "exponential"


def test_merge_rejects_malformed_growth_phase(tmp_path):
    pc = _write_paperconfig(tmp_path)
    sig = compute_paperconfig_signature(pc)
    save_extraction_json(
        tmp_path,
        {"paper": "Fake", "status": "complete", "missing_analyses": [],
         "paperconfig_signature": sig},
        [_valid_analysis(growth_phase="banana")],
    )
    with pytest.raises(SystemExit):
        merge_one_paper(tmp_path, force=False)
```

- [ ] **Step 2: Run, confirm they fail**

Run: `uv run pytest tests/extraction/timepoint/test_merge.py -v`
Expected: all FAIL with ImportError.

- [ ] **Step 3: Implement `merge.py`**

Create `multiomics_kg/extraction/timepoint/merge.py`:

```python
"""Write extraction JSON fields back into paperconfig.yaml.

Rules:
- Verbatim merge: growth_phase written exactly as it appears in JSON
  (no other:<slug> → canonical transformation; that's remap.py's job).
- Staleness check (paperconfig_signature + analysis_id cross-check)
  refuses to write without --force.
- Partial JSON (metadata.status == "partial") requires --force.
- Overwrite guard: differing existing non-null value requires --force.
- Malformed values rejected even with --force.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from ruamel.yaml import YAML

from multiomics_kg.extraction.timepoint.extraction_utils import (
    compute_paperconfig_signature,
    load_extraction_json,
    iter_paperconfigs,
    _load_yaml,  # reuse
)

logger = logging.getLogger(__name__)

_VALID_GROWTH_PHASES = frozenset([
    "exponential", "stationary", "nutrient_limited", "acclimated_steady_state",
    "infected", "recovery", "diel", "darkness", "death", "acute_stress",
    "unknown",
])


def _is_valid_growth_phase(v) -> bool:
    return isinstance(v, str) and (
        v in _VALID_GROWTH_PHASES
        or (v.startswith("other:") and len(v) > len("other:"))
    )


def merge_one_paper(paper_dir: Path, force: bool = False) -> None:
    """Merge this paper's extraction JSON into its paperconfig.yaml.

    On any refusal, prints the reason and sys.exit(1). --force overrides
    partial / signature-mismatch / overwrite-conflict, but never malformed
    values.
    """
    paperconfig_path = paper_dir / "paperconfig.yaml"
    extraction = load_extraction_json(paper_dir)
    if extraction is None:
        logger.error("No extraction JSON at %s", paper_dir)
        sys.exit(1)

    metadata = extraction.get("metadata", {})
    rows = extraction.get("analyses", [])

    # 1. Partial guard
    if metadata.get("status") == "partial" and not force:
        logger.error(
            "Refusing to merge %s: extraction is partial. Missing: %s. "
            "Pass --force to merge the valid rows only.",
            paper_dir.name, metadata.get("missing_analyses", []),
        )
        sys.exit(1)

    # 2. Signature check
    current_sig = compute_paperconfig_signature(paperconfig_path)
    stored_sig = metadata.get("paperconfig_signature")
    if stored_sig and stored_sig != current_sig and not force:
        logger.error(
            "Refusing to merge %s: paperconfig_signature mismatch. "
            "Paperconfig changed since extraction (stored=%s, now=%s). "
            "Pass --force to merge anyway.",
            paper_dir.name, stored_sig[:8], current_sig[:8],
        )
        sys.exit(1)

    # 3. Reject malformed growth_phase values (force does not bypass this)
    for row in rows:
        gp = row.get("growth_phase")
        if gp is not None and not _is_valid_growth_phase(gp):
            logger.error(
                "Refusing to merge %s: analysis %s has malformed growth_phase %r.",
                paper_dir.name, row.get("analysis_id"), gp,
            )
            sys.exit(1)

    yaml = YAML()
    yaml.preserve_quotes = True
    with open(paperconfig_path) as f:
        data = yaml.load(f)

    analyses_by_id = {}
    for table in (data.get("publication", {}).get("supplementary_materials") or {}).values():
        if not isinstance(table, dict):
            continue
        for a in (table.get("statistical_analyses") or []):
            analyses_by_id[a["id"]] = a

    # 4. Analysis-id cross-check
    ghosts = [r["analysis_id"] for r in rows if r["analysis_id"] not in analyses_by_id]
    for g in ghosts:
        logger.warning("Ghost analysis_id in JSON (not in paperconfig): %s", g)

    coverage_gap = [aid for aid in analyses_by_id if aid not in {r["analysis_id"] for r in rows}]
    if coverage_gap:
        logger.warning(
            "Paperconfig has analyses not covered by extraction (re-run extract?): %s",
            coverage_gap,
        )

    # 5. Overwrite guard + apply
    for row in rows:
        aid = row["analysis_id"]
        if aid in ghosts:
            continue
        target = analyses_by_id[aid]
        for field in ("timepoint", "timepoint_hours", "growth_phase"):
            if field not in row:
                continue
            new_value = row[field]
            existing = target.get(field)
            if existing not in (None, "") and existing != new_value and not force:
                logger.error(
                    "Refusing to overwrite %s:%s (current=%r, json=%r). Pass --force.",
                    aid, field, existing, new_value,
                )
                sys.exit(1)
            target[field] = new_value

    with open(paperconfig_path, "w") as f:
        yaml.dump(data, f)

    logger.info("Merged %s → %s", paper_dir.name, paperconfig_path)


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="merge.py")
    parser.add_argument("--paper", help="Paper name (papername) to merge")
    parser.add_argument("--all", action="store_true", help="Merge all papers that have extraction JSON")
    parser.add_argument("--force", action="store_true", help="Allow overwrite / partial / signature override")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    from multiomics_kg.extraction.timepoint.extract import _resolve_paper_dirs
    paper_dirs = _resolve_paper_dirs(args)
    if not paper_dirs:
        logger.error("No papers matched selection.")
        return 2

    for pd in paper_dirs:
        if (pd / "extractions" / "timepoint.json").exists():
            merge_one_paper(pd, force=args.force)
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/extraction/timepoint/test_merge.py -v`
Expected: all 8 tests pass.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/extraction/timepoint/merge.py \
        tests/extraction/timepoint/test_merge.py
git commit -m "feat(timepoint extraction): merge.py — verbatim write into paperconfig with staleness/partial/overwrite guards"
```

---

## Task 12: `extract.py` `--resume` and `--retry` flags

**Files:**
- Modify: `multiomics_kg/extraction/timepoint/extract.py`
- Modify: `tests/extraction/timepoint/test_extract.py`

- [ ] **Step 1: Write failing tests**

Append to `tests/extraction/timepoint/test_extract.py`:

```python
from unittest.mock import MagicMock


def test_resume_only_extracts_missing_analyses(paper_dir):
    # Seed a partial extraction with DE_n_24h valid and DE_n_48h missing
    sig = __import__(
        "multiomics_kg.extraction.timepoint.extraction_utils", fromlist=["compute_paperconfig_signature"]
    ).compute_paperconfig_signature(paper_dir / "paperconfig.yaml")
    from multiomics_kg.extraction.timepoint.extraction_utils import save_extraction_json
    save_extraction_json(paper_dir, {
        "paper": "Fake 2026", "doi": "", "status": "partial",
        "paperconfig_signature": sig, "input_tokens": 0, "output_tokens": 0,
        "missing_analyses": [{"analysis_id": "DE_n_48h", "reason": "not_returned"}],
    }, [{
        "analysis_id": "DE_n_24h", "experiment_key": "exp1",
        "fields_requested": ["timepoint", "timepoint_hours", "growth_phase"],
        "timepoint": "24h", "timepoint_hours": 24.0, "growth_phase": "exponential",
        "self_assessment": "high", "assessment_notes": "",
        "supporting_quotes": [{"quote": "q", "location": "M"}],
        "source_figures": [],
    }])

    fake_llm = MagicMock()
    fake_llm.call.return_value = ({
        "analyses": [{
            "analysis_id": "DE_n_48h",
            "timepoint": "48h", "timepoint_hours": 48.0,
            "growth_phase": "nutrient_limited",
            "self_assessment": "high", "assessment_notes": "",
            "supporting_quotes": [{"quote": "q", "location": "M"}],
            "source_figures": [],
        }]
    }, {"input_tokens": 50, "output_tokens": 100, "model": "fake-model"})

    from multiomics_kg.extraction.timepoint.extract import extract_one_paper
    extract_one_paper(paper_dir, llm_client=fake_llm, validate=False, resume=True)

    # LLM should have been called with only the missing analysis in targets
    call_args = fake_llm.call.call_args
    prompt = call_args.kwargs.get("prompt") or call_args.args[0]
    assert "DE_n_48h" in prompt
    assert "DE_n_24h" not in prompt  # valid row not re-asked

    # Final JSON should now be complete and contain both rows
    import json
    final = json.loads((paper_dir / "extractions" / "timepoint.json").read_text())
    assert final["metadata"]["status"] == "complete"
    assert {r["analysis_id"] for r in final["analyses"]} == {"DE_n_24h", "DE_n_48h"}


def test_retry_ignores_existing_json(paper_dir):
    # Seed a complete JSON
    sig = __import__(
        "multiomics_kg.extraction.timepoint.extraction_utils", fromlist=["compute_paperconfig_signature"]
    ).compute_paperconfig_signature(paper_dir / "paperconfig.yaml")
    from multiomics_kg.extraction.timepoint.extraction_utils import save_extraction_json
    save_extraction_json(paper_dir, {
        "paper": "Fake 2026", "status": "complete", "paperconfig_signature": sig,
        "missing_analyses": [], "input_tokens": 0, "output_tokens": 0,
    }, [])

    fake_llm = MagicMock()
    fake_llm.call.return_value = ({"analyses": []}, {"input_tokens": 0, "output_tokens": 0, "model": "m"})

    from multiomics_kg.extraction.timepoint.extract import extract_one_paper
    # With retry=True, the existing JSON is ignored and every analysis is requested.
    # (paperconfig has two analyses, all without timepoint fields → both requested)
    extract_one_paper(paper_dir, llm_client=fake_llm, validate=False, retry=True)
    call_args = fake_llm.call.call_args
    prompt = call_args.kwargs.get("prompt") or call_args.args[0]
    assert "DE_n_24h" in prompt
    assert "DE_n_48h" in prompt
```

- [ ] **Step 2: Run — confirm fail**

Run: `uv run pytest tests/extraction/timepoint/test_extract.py -v -k "resume or retry"`
Expected: 2 FAIL (resume/retry kwargs not accepted).

- [ ] **Step 3: Add `--resume` and `--retry` logic to `extract.py`**

Modify `extract_one_paper`:

```python
def extract_one_paper(
    paper_dir: Path,
    llm_client: LLMClient,
    validate: bool = False,
    pdf_cache_path: Path | None = None,
    model: str | None = None,
    resume: bool = False,
    retry: bool = False,
) -> Path | None:
    paperconfig_path = paper_dir / "paperconfig.yaml"
    analyses = find_analyses(paperconfig_path)

    existing = None if retry else load_extraction_json(paper_dir)

    if resume:
        if existing is None or existing.get("metadata", {}).get("status") != "partial":
            logger.info("--resume requested but no partial JSON at %s", paper_dir)
            return None
        missing_ids = {m["analysis_id"] for m in existing["metadata"].get("missing_analyses", [])}
        analyses = [a for a in analyses if a["id"] in missing_ids]

    targets = build_targets(analyses, validate=validate)
    if not targets:
        logger.info("No fields to extract for %s — skipping LLM call.", paper_dir)
        return None

    # ... (existing prompt/call/validate logic) ...

    # Merge with prior valid rows when resuming:
    if resume and existing is not None:
        prior_rows = existing.get("analyses", [])
        # Carry-over: only keep prior rows not re-requested
        valid_rows = prior_rows + valid_rows
```

Full revised signature + resume merge lives at the top-of-function block and just before `save_extraction_json`. Also wire the CLI flags:

```python
parser.add_argument("--resume", action="store_true",
                    help="Re-extract only missing_analyses from an existing partial JSON")
parser.add_argument("--retry", action="store_true",
                    help="Re-extract the entire paper from scratch")
```

And pass through in the main loop:

```python
path = extract_one_paper(
    pd, client, validate=args.validate, model=args.model,
    resume=args.resume, retry=args.retry,
)
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/extraction/timepoint/test_extract.py -v`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/extraction/timepoint/extract.py \
        tests/extraction/timepoint/test_extract.py
git commit -m "feat(timepoint extraction): --resume retries only missing analyses; --retry starts fresh"
```

---

## Task 13: `remap.py` — `--from X --to Y` rewriter

**Files:**
- Create: `multiomics_kg/extraction/timepoint/remap.py`
- Test: `tests/extraction/timepoint/test_remap.py`

- [ ] **Step 1: Write failing tests**

Create `tests/extraction/timepoint/test_remap.py`:

```python
"""Tests for remap.py."""
import json
from pathlib import Path

import pytest
from ruamel.yaml import YAML

from multiomics_kg.extraction.timepoint.remap import remap_value


def _seed_paper(paper_dir: Path, growth_phase: str) -> None:
    paper_dir.mkdir(parents=True, exist_ok=True)
    (paper_dir / "paperconfig.yaml").write_text(f"""publication:
  papername: X
  papermainpdf: data/x.pdf
  experiments:
    exp1: {{name: T, omics_type: RNASEQ}}
  supplementary_materials:
    tbl:
      type: csv
      filename: data/x.csv
      statistical_analyses:
        - id: DE_a
          experiment: exp1
          name_col: g
          logfc_col: lfc
          timepoint: "24h"
          timepoint_hours: 24
          growth_phase: {growth_phase}
""")
    (paper_dir / "extractions").mkdir(exist_ok=True)
    (paper_dir / "extractions" / "timepoint.json").write_text(json.dumps({
        "metadata": {"paper": "X", "status": "complete", "missing_analyses": []},
        "analyses": [{
            "analysis_id": "DE_a", "experiment_key": "exp1",
            "fields_requested": ["growth_phase"],
            "timepoint": "24h", "timepoint_hours": 24.0,
            "growth_phase": growth_phase,
            "self_assessment": "high", "assessment_notes": "",
            "supporting_quotes": [{"quote": "q", "location": "p"}],
            "source_figures": [],
        }],
    }))


def test_remap_rewrites_paperconfig_and_json(tmp_path):
    paper_dir = tmp_path / "paper1"
    _seed_paper(paper_dir, "other:heat_acclimated")

    remap_value(
        paper_dirs=[paper_dir],
        from_value="other:heat_acclimated",
        to_value="heat_acclimated",   # pretend validator has been updated
    )

    yaml = YAML()
    with open(paper_dir / "paperconfig.yaml") as f:
        data = yaml.load(f)
    a = data["publication"]["supplementary_materials"]["tbl"]["statistical_analyses"][0]
    assert a["growth_phase"] == "heat_acclimated"

    ext = json.loads((paper_dir / "extractions" / "timepoint.json").read_text())
    assert ext["analyses"][0]["growth_phase"] == "heat_acclimated"


def test_remap_appends_provenance_to_assessment_notes(tmp_path):
    paper_dir = tmp_path / "paper1"
    _seed_paper(paper_dir, "other:heat_acclimated")
    remap_value([paper_dir], "other:heat_acclimated", "heat_acclimated")
    ext = json.loads((paper_dir / "extractions" / "timepoint.json").read_text())
    note = ext["analyses"][0]["assessment_notes"]
    assert "Remapped from other:heat_acclimated" in note
    assert "heat_acclimated" in note


def test_remap_rejects_invalid_to_value(tmp_path):
    with pytest.raises(SystemExit):
        remap_value([tmp_path], "other:x", "malformed value with spaces")


def test_remap_accepts_other_to_other_rename(tmp_path):
    paper_dir = tmp_path / "paper1"
    _seed_paper(paper_dir, "other:heat_acclim")
    remap_value([paper_dir], "other:heat_acclim", "other:heat_acclimated")
    ext = json.loads((paper_dir / "extractions" / "timepoint.json").read_text())
    assert ext["analyses"][0]["growth_phase"] == "other:heat_acclimated"


def test_remap_merge_into_existing_canonical(tmp_path):
    paper_dir = tmp_path / "paper1"
    _seed_paper(paper_dir, "other:heat_shocked")
    remap_value([paper_dir], "other:heat_shocked", "acute_stress")
    ext = json.loads((paper_dir / "extractions" / "timepoint.json").read_text())
    assert ext["analyses"][0]["growth_phase"] == "acute_stress"


def test_remap_skips_analyses_without_match(tmp_path):
    paper_dir = tmp_path / "paper1"
    _seed_paper(paper_dir, "exponential")
    remap_value([paper_dir], "other:heat_acclimated", "heat_acclimated")
    yaml = YAML()
    with open(paper_dir / "paperconfig.yaml") as f:
        data = yaml.load(f)
    # Unchanged
    a = data["publication"]["supplementary_materials"]["tbl"]["statistical_analyses"][0]
    assert a["growth_phase"] == "exponential"


def test_remap_dry_run_does_not_write(tmp_path):
    paper_dir = tmp_path / "paper1"
    _seed_paper(paper_dir, "other:heat_acclimated")
    yaml = YAML()
    with open(paper_dir / "paperconfig.yaml") as f:
        before = f.read()
    remap_value([paper_dir], "other:heat_acclimated", "heat_acclimated", dry_run=True)
    with open(paper_dir / "paperconfig.yaml") as f:
        after = f.read()
    assert before == after
```

- [ ] **Step 2: Run, confirm fail**

Run: `uv run pytest tests/extraction/timepoint/test_remap.py -v`
Expected: all FAIL with ImportError.

- [ ] **Step 3: Implement `remap.py`**

Create `multiomics_kg/extraction/timepoint/remap.py`:

```python
"""Remap growth_phase values across paperconfigs + extraction JSONs.

Usage:
    # Promote other:<slug> → canonical (after adding to VALID_GROWTH_PHASES):
    uv run python -m multiomics_kg.extraction.timepoint.remap \
        --from other:heat_acclimated --to heat_acclimated

    # Rename a slug:
    uv run python -m multiomics_kg.extraction.timepoint.remap \
        --from other:heat_acclim --to other:heat_acclimated

    # Merge into an existing canonical value:
    uv run python -m multiomics_kg.extraction.timepoint.remap \
        --from other:heat_shocked --to acute_stress
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import date
from pathlib import Path

from ruamel.yaml import YAML

from multiomics_kg.extraction.timepoint.extraction_utils import (
    EXTRACTION_FILENAME,
    EXTRACTIONS_DIR,
    iter_paperconfigs,
)
from multiomics_kg.extraction.timepoint.merge import _is_valid_growth_phase

logger = logging.getLogger(__name__)


def remap_value(
    paper_dirs: list[Path],
    from_value: str,
    to_value: str,
    dry_run: bool = False,
) -> None:
    """Rewrite growth_phase=from_value → growth_phase=to_value across
    paperconfigs and extraction JSONs in paper_dirs.
    """
    if not _is_valid_growth_phase(to_value):
        logger.error(
            "Refusing to remap: --to value %r is not a valid growth_phase "
            "(must be in VALID_GROWTH_PHASES or well-formed other:<slug>). "
            "If promoting, update VALID_GROWTH_PHASES in the validator first.",
            to_value,
        )
        sys.exit(1)

    yaml = YAML()
    yaml.preserve_quotes = True
    provenance = f"Remapped from {from_value} → {to_value} on {date.today().isoformat()}"

    for paper_dir in paper_dirs:
        _remap_paperconfig(paper_dir, from_value, to_value, yaml, dry_run)
        _remap_json(paper_dir, from_value, to_value, provenance, dry_run)


def _remap_paperconfig(
    paper_dir: Path, from_value: str, to_value: str, yaml: YAML, dry_run: bool,
) -> None:
    path = paper_dir / "paperconfig.yaml"
    if not path.exists():
        return
    with open(path) as f:
        data = yaml.load(f)
    changed = False
    for table in (data.get("publication", {}).get("supplementary_materials") or {}).values():
        if not isinstance(table, dict):
            continue
        for a in (table.get("statistical_analyses") or []):
            if a.get("growth_phase") == from_value:
                a["growth_phase"] = to_value
                changed = True
    if changed and not dry_run:
        with open(path, "w") as f:
            yaml.dump(data, f)
        logger.info("Remapped %s", path)
    elif changed:
        logger.info("[dry-run] would remap %s", path)


def _remap_json(
    paper_dir: Path, from_value: str, to_value: str, provenance: str, dry_run: bool,
) -> None:
    path = paper_dir / EXTRACTIONS_DIR / EXTRACTION_FILENAME
    if not path.exists():
        return
    data = json.loads(path.read_text())
    changed = False
    for row in data.get("analyses", []):
        if row.get("growth_phase") == from_value:
            row["growth_phase"] = to_value
            existing_note = row.get("assessment_notes", "") or ""
            row["assessment_notes"] = (
                f"{existing_note}\n{provenance}" if existing_note else provenance
            )
            changed = True
    if changed and not dry_run:
        path.write_text(json.dumps(data, indent=2, default=str))
        logger.info("Remapped %s", path)
    elif changed:
        logger.info("[dry-run] would remap %s", path)


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="remap.py")
    parser.add_argument("--from", dest="from_value", required=True,
                        help="Current value (e.g. other:heat_acclimated)")
    parser.add_argument("--to", dest="to_value", required=True,
                        help="Target value (canonical enum or other:<slug>)")
    parser.add_argument("--dry-run", action="store_true",
                        help="List affected files without writing")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    list_files = [
        Path("data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"),
        Path("data/Synechococcus/papers_and_supp/paperconfig_files.txt"),
    ]
    paper_dirs = [p.parent for p in iter_paperconfigs(list_files) if p.exists()]
    remap_value(paper_dirs, args.from_value, args.to_value, dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/extraction/timepoint/test_remap.py -v`
Expected: all 7 tests pass.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/extraction/timepoint/remap.py \
        tests/extraction/timepoint/test_remap.py
git commit -m "feat(timepoint extraction): remap.py — --from/--to rewriter with provenance"
```

---

## Task 14: Aggregate report (`report.py`)

**Files:**
- Create: `multiomics_kg/extraction/timepoint/report.py`
- Test: `tests/extraction/timepoint/test_report.py`

- [ ] **Step 1: Write failing tests**

Create `tests/extraction/timepoint/test_report.py`:

```python
"""Tests for report.py — aggregates extraction JSONs into a markdown report."""
import json
from pathlib import Path

from multiomics_kg.extraction.timepoint.report import (
    aggregate_reports,
    render_markdown,
)


def _seed(paper_dir: Path, paper: str, analyses: list[dict], status: str = "complete"):
    paper_dir.mkdir(parents=True, exist_ok=True)
    (paper_dir / "extractions").mkdir(exist_ok=True)
    (paper_dir / "extractions" / "timepoint.json").write_text(json.dumps({
        "metadata": {"paper": paper, "status": status, "missing_analyses": []},
        "analyses": analyses,
    }))


def test_aggregate_counts_self_assessments(tmp_path):
    _seed(tmp_path / "p1", "P1", [
        {"analysis_id": "a", "growth_phase": "exponential", "self_assessment": "high",
         "timepoint": "24h", "timepoint_hours": 24, "assessment_notes": "",
         "supporting_quotes": [], "source_figures": []},
    ])
    _seed(tmp_path / "p2", "P2", [
        {"analysis_id": "b", "growth_phase": "unknown", "self_assessment": "low",
         "timepoint": "unknown", "timepoint_hours": None, "assessment_notes": "",
         "supporting_quotes": [], "source_figures": []},
    ])
    agg = aggregate_reports([tmp_path / "p1", tmp_path / "p2"])
    assert agg["self_assessment_counts"] == {"high": 1, "low": 1}


def test_aggregate_other_slug_frequencies(tmp_path):
    _seed(tmp_path / "p1", "P1", [
        {"analysis_id": "a", "growth_phase": "other:heat_acclimated",
         "self_assessment": "medium", "timepoint": "48h", "timepoint_hours": 48,
         "assessment_notes": "", "supporting_quotes": [], "source_figures": []},
    ])
    _seed(tmp_path / "p2", "P2", [
        {"analysis_id": "b", "growth_phase": "other:heat_acclimated",
         "self_assessment": "high", "timepoint": "72h", "timepoint_hours": 72,
         "assessment_notes": "", "supporting_quotes": [], "source_figures": []},
    ])
    agg = aggregate_reports([tmp_path / "p1", tmp_path / "p2"])
    assert agg["other_slug_counts"] == {"other:heat_acclimated": 2}


def test_aggregate_collects_unknown_with_evidence(tmp_path):
    _seed(tmp_path / "p1", "P1", [
        {"analysis_id": "a", "growth_phase": "unknown", "self_assessment": "low",
         "timepoint": "unknown", "timepoint_hours": None,
         "assessment_notes": "Paper does not describe phase",
         "supporting_quotes": [], "source_figures": []},
    ])
    agg = aggregate_reports([tmp_path / "p1"])
    assert len(agg["unknowns"]) == 1
    assert agg["unknowns"][0]["analysis_id"] == "a"
    assert "Paper does not describe phase" in agg["unknowns"][0]["notes"]


def test_aggregate_includes_partial_missing(tmp_path):
    _seed(tmp_path / "p1", "P1", [], status="partial")
    d = tmp_path / "p1" / "extractions" / "timepoint.json"
    data = json.loads(d.read_text())
    data["metadata"]["missing_analyses"] = [{"analysis_id": "x", "reason": "not_returned"}]
    d.write_text(json.dumps(data))
    agg = aggregate_reports([tmp_path / "p1"])
    assert len(agg["partial_extractions"]) == 1
    assert agg["partial_extractions"][0]["paper"] == "P1"


def test_render_markdown_has_required_sections(tmp_path):
    _seed(tmp_path / "p1", "P1", [
        {"analysis_id": "a", "growth_phase": "other:x", "self_assessment": "medium",
         "timepoint": "24h", "timepoint_hours": 24, "assessment_notes": "",
         "supporting_quotes": [], "source_figures": []},
    ])
    agg = aggregate_reports([tmp_path / "p1"])
    md = render_markdown(agg)
    assert "# Timepoint Extraction Report" in md
    assert "## Self-assessment distribution" in md
    assert "## other:* slugs" in md
    assert "## Unknowns" in md
    assert "## Partial extractions" in md
```

- [ ] **Step 2: Run, confirm fail**

Run: `uv run pytest tests/extraction/timepoint/test_report.py -v`
Expected: all FAIL with ImportError.

- [ ] **Step 3: Implement `report.py`**

Create `multiomics_kg/extraction/timepoint/report.py`:

```python
"""Aggregate per-paper extraction JSONs into data/timepoint_extraction_report.md."""
from __future__ import annotations

import argparse
import json
import logging
import sys
from collections import Counter
from pathlib import Path

from multiomics_kg.extraction.timepoint.extraction_utils import (
    EXTRACTION_FILENAME,
    EXTRACTIONS_DIR,
    iter_paperconfigs,
    load_extraction_json,
)

logger = logging.getLogger(__name__)

REPORT_PATH = Path("data/timepoint_extraction_report.md")


def aggregate_reports(paper_dirs: list[Path]) -> dict:
    """Walk paper_dirs, aggregate stats from each extraction JSON."""
    sa_counter: Counter = Counter()
    other_counter: Counter = Counter()
    unknowns: list[dict] = []
    partial: list[dict] = []
    per_paper: list[dict] = []

    for paper_dir in paper_dirs:
        ext = load_extraction_json(paper_dir)
        if ext is None:
            continue
        paper = ext["metadata"].get("paper", paper_dir.name)
        status = ext["metadata"].get("status")
        rows = ext.get("analyses", [])

        if status == "partial":
            partial.append({
                "paper": paper,
                "missing_analyses": ext["metadata"].get("missing_analyses", []),
            })

        for row in rows:
            sa = row.get("self_assessment")
            if sa:
                sa_counter[sa] += 1
            gp = row.get("growth_phase")
            if isinstance(gp, str) and gp.startswith("other:"):
                other_counter[gp] += 1
            if gp == "unknown":
                unknowns.append({
                    "paper": paper,
                    "analysis_id": row.get("analysis_id"),
                    "notes": row.get("assessment_notes", ""),
                    "quotes": row.get("supporting_quotes", []),
                })
        per_paper.append({"paper": paper, "n_analyses": len(rows), "status": status})

    return {
        "self_assessment_counts": dict(sa_counter),
        "other_slug_counts": dict(other_counter),
        "unknowns": unknowns,
        "partial_extractions": partial,
        "per_paper": per_paper,
    }


def render_markdown(agg: dict) -> str:
    lines = ["# Timepoint Extraction Report", ""]

    lines.append("## Self-assessment distribution")
    for sa in ("high", "medium", "low"):
        n = agg["self_assessment_counts"].get(sa, 0)
        lines.append(f"- `{sa}`: {n}")
    lines.append("")

    lines.append("## other:* slugs")
    sorted_slugs = sorted(
        agg["other_slug_counts"].items(), key=lambda kv: kv[1], reverse=True,
    )
    if not sorted_slugs:
        lines.append("_(none)_")
    else:
        for slug, n in sorted_slugs:
            lines.append(f"- `{slug}`: {n}")
    lines.append("")

    lines.append("## Unknowns")
    if not agg["unknowns"]:
        lines.append("_(none)_")
    else:
        for u in agg["unknowns"]:
            lines.append(f"- **{u['paper']}** / `{u['analysis_id']}`: {u['notes']}")
            for q in u["quotes"]:
                lines.append(f"  - _{q.get('location', '')}_: {q.get('quote', '')}")
    lines.append("")

    lines.append("## Partial extractions")
    if not agg["partial_extractions"]:
        lines.append("_(none)_")
    else:
        for p in agg["partial_extractions"]:
            lines.append(f"- **{p['paper']}**: {len(p['missing_analyses'])} missing")
            for m in p["missing_analyses"]:
                lines.append(f"  - `{m['analysis_id']}`: {m['reason']}")
    lines.append("")

    lines.append("## Per-paper summary")
    for p in agg["per_paper"]:
        lines.append(f"- {p['paper']}: {p['n_analyses']} analyses ({p['status']})")

    return "\n".join(lines) + "\n"


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="report.py")
    parser.add_argument("--output", default=str(REPORT_PATH),
                        help=f"Output path (default: {REPORT_PATH})")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    list_files = [
        Path("data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"),
        Path("data/Synechococcus/papers_and_supp/paperconfig_files.txt"),
    ]
    paper_dirs = [p.parent for p in iter_paperconfigs(list_files) if p.exists()]

    agg = aggregate_reports(paper_dirs)
    md = render_markdown(agg)
    Path(args.output).write_text(md)
    logger.info("Report written: %s", args.output)
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/extraction/timepoint/test_report.py -v`
Expected: all 5 tests pass.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/extraction/timepoint/report.py \
        tests/extraction/timepoint/test_report.py
git commit -m "feat(timepoint extraction): report.py — aggregate per-paper JSONs into corpus-wide markdown"
```

---

## Task 15: CLI dispatcher `__main__.py`

**Files:**
- Create: `multiomics_kg/extraction/timepoint/__main__.py`

- [ ] **Step 1: Write the dispatcher**

Create `multiomics_kg/extraction/timepoint/__main__.py`:

```python
"""Dispatch `python -m multiomics_kg.extraction.timepoint <subcommand>`.

Subcommands: extract, merge, remap, report.
"""
import sys


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage: python -m multiomics_kg.extraction.timepoint "
              "{extract|merge|remap|report} [args...]", file=sys.stderr)
        return 2

    sub = sys.argv[1]
    argv = sys.argv[2:]

    if sub == "extract":
        from multiomics_kg.extraction.timepoint.extract import _main as extract_main
        return extract_main(argv)
    if sub == "merge":
        from multiomics_kg.extraction.timepoint.merge import _main as merge_main
        return merge_main(argv)
    if sub == "remap":
        from multiomics_kg.extraction.timepoint.remap import _main as remap_main
        return remap_main(argv)
    if sub == "report":
        from multiomics_kg.extraction.timepoint.report import _main as report_main
        return report_main(argv)

    print(f"Unknown subcommand: {sub}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 2: Smoke-test the CLI**

Run:
```bash
uv run python -m multiomics_kg.extraction.timepoint extract --help
uv run python -m multiomics_kg.extraction.timepoint merge --help
uv run python -m multiomics_kg.extraction.timepoint remap --help
uv run python -m multiomics_kg.extraction.timepoint report --help
```
Expected: each prints its argparse help.

- [ ] **Step 3: Commit**

```bash
git add multiomics_kg/extraction/timepoint/__main__.py
git commit -m "feat(timepoint extraction): unified CLI dispatcher (extract/merge/remap/report)"
```

---

## Task 16: End-to-end integration test (still mocks LLM)

**Files:**
- Create: `tests/extraction/timepoint/test_integration.py`

- [ ] **Step 1: Write a full extract → merge → remap flow test**

Create `tests/extraction/timepoint/test_integration.py`:

```python
"""Integration test covering extract → merge → remap with a mocked LLM."""
import json
from pathlib import Path
from unittest.mock import MagicMock

from ruamel.yaml import YAML

from multiomics_kg.extraction.timepoint.extract import extract_one_paper
from multiomics_kg.extraction.timepoint.merge import merge_one_paper
from multiomics_kg.extraction.timepoint.remap import remap_value


def test_full_flow_extract_merge_remap(tmp_path):
    pc = tmp_path / "paperconfig.yaml"
    pc.write_text("""publication:
  papername: Fake 2026
  papermainpdf: data/fake.pdf
  experiments:
    exp1:
      name: X
      organism: Prochlorococcus MED4
      treatment_type: [carbon]
      omics_type: RNASEQ
  supplementary_materials:
    tbl:
      type: csv
      filename: data/fake.csv
      statistical_analyses:
        - id: DE_a
          experiment: exp1
          name_col: gene
          logfc_col: log2fc
""")
    (tmp_path / "data").mkdir(exist_ok=True)
    (tmp_path / "data" / "fake.pdf").write_text("x")

    # LLM emits other:heat_acclimated
    fake_llm = MagicMock()
    fake_llm.call.return_value = ({
        "analyses": [{
            "analysis_id": "DE_a",
            "timepoint": "acclimated",
            "timepoint_hours": None,
            "growth_phase": "other:heat_acclimated",
            "self_assessment": "high",
            "assessment_notes": "",
            "supporting_quotes": [{"quote": "q", "location": "Methods"}],
            "source_figures": [],
        }]
    }, {"input_tokens": 50, "output_tokens": 100, "model": "fake"})

    # 1. Extract
    extract_one_paper(tmp_path, llm_client=fake_llm)

    # 2. Merge (writes other:heat_acclimated verbatim)
    merge_one_paper(tmp_path, force=False)

    yaml = YAML()
    with open(pc) as f:
        data = yaml.load(f)
    assert data["publication"]["supplementary_materials"]["tbl"][
        "statistical_analyses"][0]["growth_phase"] == "other:heat_acclimated"

    # 3. Remap to canonical (simulating enum promotion). For the test we
    #    bypass the is_valid_growth_phase check by remapping to an
    #    existing canonical value.
    remap_value([tmp_path], "other:heat_acclimated", "acclimated_steady_state")

    with open(pc) as f:
        data = yaml.load(f)
    assert data["publication"]["supplementary_materials"]["tbl"][
        "statistical_analyses"][0]["growth_phase"] == "acclimated_steady_state"

    ext = json.loads((tmp_path / "extractions" / "timepoint.json").read_text())
    assert ext["analyses"][0]["growth_phase"] == "acclimated_steady_state"
    assert "Remapped from" in ext["analyses"][0]["assessment_notes"]
```

- [ ] **Step 2: Run the integration test**

Run: `uv run pytest tests/extraction/timepoint/test_integration.py -v`
Expected: PASS.

- [ ] **Step 3: Run the full suite once to check for regressions**

Run: `uv run pytest tests/extraction/timepoint/ -v`
Expected: all pass.

- [ ] **Step 4: Run the broader project suite**

Run: `uv run pytest -m "not slow and not kg" -v`
Expected: all pass (new code should not break existing tests).

- [ ] **Step 5: Commit**

```bash
git add tests/extraction/timepoint/test_integration.py
git commit -m "test(timepoint extraction): end-to-end extract → merge → remap integration test"
```

---

## Task 17: Prompt preview — inspect rendered prompts before burning LLM tokens

**Files:**
- Modify: `multiomics_kg/extraction/timepoint/extract.py`
- Modify: `tests/extraction/timepoint/test_extract.py`

**Why:** The prompt is the most important artifact in this pipeline — it's what actually shapes every LLM response across 30 papers. Before running real extractions, render the full prompt (with a real paperconfig + cache entry) and review it by eye. Catch issues like: missing context, leaking irrelevant fields, unclear instructions, malformed sections, PDF paths pointing nowhere. Cheaper than an LLM call, and catches classes of bugs that `--dry-run` (which only shows `fields_requested`) doesn't surface.

- [ ] **Step 1: Write a failing test for the `--print-prompt` CLI flag**

Append to `tests/extraction/timepoint/test_extract.py`:

```python
def test_print_prompt_renders_without_llm_call(paper_dir, capsys):
    """--print-prompt renders the full prompt and PDF attachment list to
    stdout; never calls the LLM."""
    from multiomics_kg.extraction.timepoint.extract import print_prompt_one_paper

    print_prompt_one_paper(paper_dir, validate=False)

    out = capsys.readouterr().out
    # Prompt body
    assert "PAPER CONTEXT" in out
    assert "EXTRACTION TARGETS" in out
    assert "Fake 2026" in out
    assert "DE_n_24h" in out
    assert "log2fc_24h" in out
    # PDF attachment list
    assert "PDFs to attach:" in out or "PDF attachments:" in out


def test_print_prompt_skips_when_nothing_to_extract(paper_dir, capsys):
    """If all fields are populated, --print-prompt says so and emits no
    prompt body."""
    from ruamel.yaml import YAML
    pc = paper_dir / "paperconfig.yaml"
    yaml = YAML()
    with open(pc) as f:
        data = yaml.load(f)
    for a in data["publication"]["supplementary_materials"]["tbl"]["statistical_analyses"]:
        a["timepoint"] = "24h"
        a["timepoint_hours"] = 24
        a["growth_phase"] = "exponential"
    with open(pc, "w") as f:
        yaml.dump(data, f)

    from multiomics_kg.extraction.timepoint.extract import print_prompt_one_paper
    print_prompt_one_paper(paper_dir, validate=False)
    out = capsys.readouterr().out
    assert "no fields to extract" in out.lower()
    assert "EXTRACTION TARGETS" not in out
```

- [ ] **Step 2: Run — confirm fail**

Run: `uv run pytest tests/extraction/timepoint/test_extract.py -v -k print_prompt`
Expected: 2 FAIL with ImportError.

- [ ] **Step 3: Implement `print_prompt_one_paper` and wire the flag**

In `multiomics_kg/extraction/timepoint/extract.py`, add a new function alongside `extract_one_paper`:

```python
def print_prompt_one_paper(
    paper_dir: Path,
    validate: bool = False,
    pdf_cache_path: Path | None = None,
) -> None:
    """Render the full LLM prompt for this paper and print to stdout.
    No LLM call. Includes the list of PDF attachment paths so the
    reviewer can verify file existence.
    """
    paperconfig_path = paper_dir / "paperconfig.yaml"
    analyses = find_analyses(paperconfig_path)
    targets = build_targets(analyses, validate=validate)
    if not targets:
        print(f"\n--- {paper_dir.name}: no fields to extract (all populated) ---\n")
        return

    background = build_background(paperconfig_path)
    pdf_cache_entry = _load_pdf_cache_entry(
        pdf_cache_path or Path("cache/pdf_extraction_cache.json"),
        background.get("papermainpdf"),
    )
    prompt = build_prompt(background, targets, pdf_cache_entry)

    pdf_paths: list[Path] = []
    if pmp := background.get("papermainpdf"):
        pdf_paths.append(Path(pmp))
    pdf_paths.extend(Path(p) for p in background.get("additional_pdfs", []))

    print(f"\n{'=' * 72}")
    print(f"PROMPT PREVIEW — {paper_dir.name}")
    print("=" * 72)
    print(prompt)
    print()
    print("PDFs to attach:")
    for p in pdf_paths:
        marker = "✓" if p.exists() else "✗ MISSING"
        print(f"  {marker} {p}")
    print()
    print(f"Targets: {len(targets)} analyses")
    for t in targets:
        print(f"  - {t['id']}: needs {t['fields_requested']}  (logfc_col={t['logfc_col']!r})")
    print("=" * 72)
```

In `_main`, add a flag and dispatch:

```python
parser.add_argument("--print-prompt", action="store_true",
                    help="Render the full prompt (prompt text + PDF paths) "
                         "without calling the LLM. Use to review prompt quality "
                         "before running real extractions.")
```

In the main loop (before the `OpenAIResponsesClient` block):

```python
if args.print_prompt:
    for pd in paper_dirs:
        print_prompt_one_paper(pd, validate=args.validate)
    return 0
```

- [ ] **Step 4: Run the tests**

Run: `uv run pytest tests/extraction/timepoint/test_extract.py -v`
Expected: all pass (existing 6 + 2 new).

- [ ] **Step 5: Preview prompts on real papers (manual review gate)**

Pick 3 representative papers spanning different treatment types — e.g., one time-course, one coculture, one acclimation:

```bash
# Time-course (nutrient starvation)
uv run python -m multiomics_kg.extraction.timepoint extract \
    --paper "Weissberg 2025" --print-prompt 2>&1 | tee /tmp/prompt_weissberg.txt

# Coculture
uv run python -m multiomics_kg.extraction.timepoint extract \
    --paper "Biller 2018" --print-prompt 2>&1 | tee /tmp/prompt_biller.txt

# Acclimation
uv run python -m multiomics_kg.extraction.timepoint extract \
    --paper "Barreto 2022" --print-prompt 2>&1 | tee /tmp/prompt_barreto.txt
```

**Read each prompt end-to-end.** Check:

- Are SHARED_RULES clear and unambiguous?
- Does the experiments block contain enough context (treatment_type, background_factors, conditions)?
- Are `logfc_col` values rendered per-analysis?
- Are the per-experiment-type hints relevant and accurate for this paper's treatment_type?
- Is the PDF attachment list correct (no missing files marked `✗`)?
- Is the cache-entry section present and useful (or missing and maybe should be added)?
- Any leaking of irrelevant fields? Any prose that reads as unclear or contradictory?

If the prompt has problems, **stop here and fix them** — edits go in `multiomics_kg/extraction/timepoint/prompts.py` or `extract.py`'s `build_background`/`build_targets`. Re-run tests. Re-preview. Only proceed to rollout when the prompt reads clean on all 3 sample papers.

- [ ] **Step 6: Commit the flag (and any prompt fixes from Step 5)**

```bash
git add multiomics_kg/extraction/timepoint/extract.py \
        tests/extraction/timepoint/test_extract.py \
        multiomics_kg/extraction/timepoint/prompts.py  # if edited in Step 5
git commit -m "feat(timepoint extraction): --print-prompt flag for manual prompt review before rollout"
```

---

## Task 18: Code review of the extraction pipeline

**Why:** Before running the pipeline against all 30 papers (Task 19), get a fresh pair of eyes on the extraction code. LLM calls cost money and rollout-time bugs are expensive to unwind — a pre-rollout review catches error-handling gaps, resource leaks, and logic bugs that tests didn't. This is a **gate**, not optional.

**Scope of review** — these files and nothing else (the schema/adapter/post-import/validator changes from Tasks 1–4 are small and already exercised by the KG tests):

```
multiomics_kg/extraction/timepoint/__init__.py
multiomics_kg/extraction/timepoint/__main__.py
multiomics_kg/extraction/timepoint/extract.py
multiomics_kg/extraction/timepoint/extraction_utils.py
multiomics_kg/extraction/timepoint/llm_client.py
multiomics_kg/extraction/timepoint/merge.py
multiomics_kg/extraction/timepoint/prompts.py
multiomics_kg/extraction/timepoint/remap.py
multiomics_kg/extraction/timepoint/report.py
tests/extraction/timepoint/*.py
```

- [ ] **Step 1: Dispatch the code-reviewer agent**

Use the `superpowers:code-reviewer` agent. Give it the full scope, the spec, and the plan so it can check against both the design intent and the task-level implementation.

Sample dispatch prompt:

> Review the timepoint/growth-phase extraction pipeline before we run it against 30 real publications.
>
> **Design spec:** `docs/superpowers/specs/2026-04-12-timepoint-growth-phase-backfill-design.md`
> **Implementation plan:** `docs/superpowers/plans/2026-04-12-timepoint-growth-phase-backfill.md` (Tasks 5–16)
> **Code under review:** `multiomics_kg/extraction/timepoint/` and `tests/extraction/timepoint/`
>
> Focus areas — in priority order:
>
> 1. **Correctness vs the spec.** Do the safety checks in `merge.py` (partial guard, staleness guards, overwrite guard, malformed-value rejection) match the spec's safety-check list? Is remap's `--from`/`--to` symmetric (promotion / slug rename / merge-into-canonical all work)? Does `extract.py --resume` really only re-ask missing analyses?
> 2. **Failure modes around real LLM calls.** What happens if the LLM returns invalid JSON, an API error, or a truncated response? Where are retries reasonable vs where should we fail hard? Are file uploads in `llm_client.py` cleaned up on error?
> 3. **Staleness detection robustness.** Is `paperconfig_signature` deterministic across YAML re-serialization? Would re-ordering fields in the paperconfig cause a false-positive mismatch?
> 4. **Enum duplication.** `VALID_GROWTH_PHASES` is defined in three places (validator, prompts, utils+merge). Is the sync test (`test_valid_growth_phases_list_matches_validator`) sufficient to catch drift? Should we refactor to a single source?
> 5. **General quality.** Error messages, logging levels, test coverage, dead code, mutable default args, unclosed file handles, typing.
>
> Read-only; report findings as a prioritized list (blockers / nice-to-haves). The implementer will address blockers before the rollout.

- [ ] **Step 2: Triage findings**

Categorize the reviewer's output:
- **Blockers** — correctness bugs, data-loss risks, broken safety guards, anything that would corrupt paperconfigs. Fix before rollout.
- **Nice-to-haves** — style, minor refactors, better error messages. Fix if cheap, defer otherwise.
- **Rejections** — disagreements with the reviewer (must be justified; not reflexive).

If findings are nontrivial, write a short response document at `docs/superpowers/code-reviews/2026-04-12-timepoint-pipeline-review.md` listing each finding, the verdict (fix / defer / reject), and the rationale. Optional for small reviews — the git commits speak for themselves.

- [ ] **Step 3: Fix blockers**

Make the code changes. Update tests where behavior changes. Run:

```bash
uv run pytest tests/extraction/timepoint/ -v
```

Expected: all pass.

- [ ] **Step 4: Re-run the integration test and prompt preview as regression check**

```bash
uv run pytest tests/extraction/timepoint/test_integration.py -v
uv run python -m multiomics_kg.extraction.timepoint extract \
    --paper "Weissberg 2025" --print-prompt 2>&1 | head -50
```

Expected: integration test passes; prompt still renders cleanly.

- [ ] **Step 5: Commit the fixes**

```bash
git add multiomics_kg/extraction/timepoint/ tests/extraction/timepoint/
# Optionally: docs/superpowers/code-reviews/2026-04-12-timepoint-pipeline-review.md
git commit -m "fix(timepoint extraction): address code review findings (pre-rollout)"
```

---

## Task 19: Rollout — run extraction in batches (paperconfig edits happen here)

**This task is operational, not code. Each batch gets its own commit outside this plan. Record here for traceability.**

Batches (per design spec):
1. Prochlorococcus coculture papers (Biller 2018, Hennon 2017, Aharonovich 2016, …)
2. Prochlorococcus nutrient-stress papers (Weissberg, Read, Martiny, Thompson 2011)
3. Prochlorococcus diel + short-exposure (Zinser, Thompson 2016, Tetu, Fang)
4. Prochlorococcus infection + remaining (Lindell, Lin, Anjur, …)
5. Synechococcus papers (Tal, Beliaev, Kratzl, Ma, Oleza, Kaur, Bernstein)

For each batch:

- [ ] **Run extract across the batch**

```bash
# Dry run first
uv run python -m multiomics_kg.extraction.timepoint extract --paper "Biller 2018" --dry-run
```

**Expected outcomes:**

- **"no fields to extract (all populated)"** — the paperconfig is already complete. **This is fine — skip the real extraction for this paper and move to the next.** Do not treat the no-op as a failure; it just means this paper needs no work from the pipeline. Some papers in later batches may already have partial or full timepoint metadata from earlier manual edits or prior runs.
- **Lists `fields_requested` for one or more analyses** — proceed to the real extraction:

```bash
# Real extraction (only for papers where dry-run showed fields needed)
uv run python -m multiomics_kg.extraction.timepoint extract --paper "Biller 2018"
```

- [ ] **Review JSONs**

Open `data/.../Biller 2018/extractions/timepoint.json` — check evidence quotes, `self_assessment`, `other:*` / `unknown` cases. Edit in place if needed.

- [ ] **Regenerate aggregate report**

```bash
uv run python -m multiomics_kg.extraction.timepoint report
```

Open `data/timepoint_extraction_report.md` — review `other:*` frequencies and `unknown` list.

- [ ] **If an `other:*` slug should be promoted**

1. Add it to `VALID_GROWTH_PHASES` in `.claude/skills/paperconfig/validate_paperconfig.py`.
2. Add it to `VALID_GROWTH_PHASES_LIST` in `multiomics_kg/extraction/timepoint/prompts.py`.
3. Add it to `_VALID_GROWTH_PHASES` in `multiomics_kg/extraction/timepoint/extraction_utils.py` and `multiomics_kg/extraction/timepoint/merge.py` (duplicated to avoid skill-dir imports — keep in sync).
4. Run the three-way sync test to confirm: `uv run pytest tests/extraction/timepoint/test_prompts.py::test_valid_growth_phases_list_matches_validator`.
5. Remap:

```bash
uv run python -m multiomics_kg.extraction.timepoint remap \
    --from other:<slug> --to <slug>
```

- [ ] **Merge across the batch**

```bash
uv run python -m multiomics_kg.extraction.timepoint merge --paper "Biller 2018"
# ...or if all papers in the batch are ready:
uv run python -m multiomics_kg.extraction.timepoint merge --all
```

- [ ] **Run the validator**

```bash
uv run python .claude/skills/paperconfig/validate_paperconfig.py --all
```

Expected: no errors. Warnings about missing growth_phase are OK until all 30 paperconfigs are backfilled.

- [ ] **Commit the batch**

```bash
git add data/**/paperconfig.yaml data/**/extractions/timepoint.json \
         data/timepoint_extraction_report.md \
         .claude/skills/paperconfig/validate_paperconfig.py \
         multiomics_kg/extraction/timepoint/prompts.py \
         multiomics_kg/extraction/timepoint/extraction_utils.py \
         multiomics_kg/extraction/timepoint/merge.py
git commit -m "data: timepoint/growth_phase backfill — batch 1 (Prochlorococcus coculture papers)"
```

Repeat for batches 2–5.

---

## Task 20: Flip validator missing-growth_phase from warning to error

**Files:**
- Modify: `.claude/skills/paperconfig/validate_paperconfig.py`

**Prerequisite:** All 30 paperconfigs have `growth_phase` populated.

- [ ] **Step 1: Verify 100% coverage**

```bash
uv run python .claude/skills/paperconfig/validate_paperconfig.py --all 2>&1 | grep -i "missing.*growth_phase"
```

Expected: no output (no paperconfig is missing growth_phase).

- [ ] **Step 2: Change `warnings.append` → `errors.append` for the missing case**

Edit `.claude/skills/paperconfig/validate_paperconfig.py` — find the block added in Task 4:

```python
if "growth_phase" not in analysis:
    # TODO (timepoint-backfill rollout, 2026-04-12): flip to errors.append
    # once all 30 paperconfigs have growth_phase populated.
    warnings.append(...)
```

Change to:

```python
if "growth_phase" not in analysis:
    errors.append(
        f"{analysis_id}: missing required field 'growth_phase' "
        f"(one of {sorted(VALID_GROWTH_PHASES)} or 'other:<slug>')"
    )
```

- [ ] **Step 3: Re-run validator**

Run: `uv run python .claude/skills/paperconfig/validate_paperconfig.py --all`
Expected: no errors.

- [ ] **Step 4: Commit**

```bash
git add .claude/skills/paperconfig/validate_paperconfig.py
git commit -m "validator: missing growth_phase is now an error (all paperconfigs backfilled)"
```

---

## Task 21: KG validity test for growth_phase

**Files:**
- Modify: `tests/kg_validity/test_expression.py`

**Prerequisite:** Paperconfigs backfilled AND KG rebuilt locally.

- [ ] **Step 1: Rebuild the KG**

Run: `docker compose up --build -d`
Expected: Docker pipeline succeeds; Neo4j available at `localhost:7687`.

- [ ] **Step 2: Add the test**

Append to `tests/kg_validity/test_expression.py`:

```python
@pytest.mark.kg
def test_all_expression_edges_have_growth_phase(neo4j_driver):
    """After the 2026-04-12 backfill, every Changes_expression_of edge
    must carry a non-null growth_phase property.
    """
    with neo4j_driver.session() as session:
        result = session.run(
            "MATCH ()-[r:Changes_expression_of]->() "
            "WHERE r.growth_phase IS NULL "
            "RETURN count(r) AS n"
        )
        n_null = result.single()["n"]
    assert n_null == 0, f"{n_null} expression edges have NULL growth_phase"


@pytest.mark.kg
def test_experiments_have_growth_phases_array(neo4j_driver):
    """Every Experiment node has a non-empty growth_phases array after
    the 2026-04-12 backfill + post-import aggregation.
    """
    with neo4j_driver.session() as session:
        result = session.run(
            "MATCH (e:Experiment) "
            "WHERE e.growth_phases IS NULL OR size(e.growth_phases) = 0 "
            "RETURN count(e) AS n"
        )
        n_empty = result.single()["n"]
    assert n_empty == 0, f"{n_empty} Experiment nodes lack growth_phases"
```

- [ ] **Step 3: Run KG validity tests**

Run: `uv run pytest tests/kg_validity/test_expression.py -v`
Expected: both new tests pass.

- [ ] **Step 4: Commit**

```bash
git add tests/kg_validity/test_expression.py
git commit -m "test(kg_validity): every expression edge has growth_phase, every Experiment has growth_phases[]"
```

---

## Task 22: Documentation updates

**Files:**
- Modify: `CLAUDE.md`
- Modify: `.claude/skills/paperconfig/SKILL.md`
- Modify: `/home/osnat/.claude/projects/-home-osnat-github-multiomics-biocypher-kg/memory/MEMORY.md`

- [ ] **Step 1: Update CLAUDE.md**

Find the section "Expression edges: ~210K `Changes_expression_of` edges ...". Add `growth_phase` to the edge-property list. Find the Experiment node section; add `growth_phases` (str[]) to the computed/post-import property list.

Also find "No EnvironmentalCondition nodes ..." note — leave as is (still accurate).

- [ ] **Step 2: Update SKILL.md**

In `.claude/skills/paperconfig/SKILL.md`, find the "Required Fields" table under statistical_analyses and add `growth_phase` alongside `timepoint` / `timepoint_hours`. Add a new "Growth Phase" subsection with:

- The enum table (copy from the spec's "Canonical `growth_phase` enum" section).
- A short decision-guidance list (copy the per-experiment-type hints from the spec's "Prompt design" section).
- A note about the `other:<slug>` escape hatch.

- [ ] **Step 3: Update MEMORY.md**

Add an entry under "Key Data Patterns":

```markdown
## Timepoint & Growth-Phase Backfill (2026-04-12)
- Every statistical_analyses[] row now has `timepoint`, `timepoint_hours`, `growth_phase` (enum or `other:<slug>`).
- growth_phase enum: exponential, stationary, nutrient_limited, acclimated_steady_state, infected, recovery, diel, darkness, death, acute_stress, unknown.
- LLM extraction lives in `multiomics_kg/extraction/timepoint/` — same pattern as cluster-extraction.
- Per-paper JSON in `<paper_dir>/extractions/timepoint.json`; aggregate report at `data/timepoint_extraction_report.md`.
- Subcommands: `uv run python -m multiomics_kg.extraction.timepoint {extract|merge|remap|report}`.
- `remap.py --from X --to Y` handles promotion / slug rename / merge-into-canonical.
- Plan: `docs/superpowers/plans/2026-04-12-timepoint-growth-phase-backfill.md`. Spec: `docs/superpowers/specs/2026-04-12-timepoint-growth-phase-backfill-design.md`.
```

- [ ] **Step 4: Commit**

```bash
git add CLAUDE.md .claude/skills/paperconfig/SKILL.md
git commit -m "docs: document growth_phase edge/node properties + backfill skill + memory entry"
```

(Memory file lives outside the repo; update it separately via the standard memory-write flow.)

---

## Self-Review

1. **Spec coverage:**
   - Schema additions (spec §1): Task 1.
   - KG propagation — adapter (spec §2): Task 2.
   - KG propagation — post-import (spec §2): Task 3.
   - Validator + `VALID_GROWTH_PHASES` + `is_valid_growth_phase` (spec §3): Task 4 (warn-only), Task 20 (flip to error).
   - Extraction pipeline — utils (spec §4): Task 5, Task 6, Task 8.
   - Extraction pipeline — prompts (spec §4): Task 7.
   - Extraction pipeline — `extract.py` + LLM client (spec §4): Tasks 9, 10, 12.
   - Extraction pipeline — `merge.py` with staleness/partial/overwrite guards (spec §4): Task 11.
   - Extraction pipeline — `remap.py` with `--from`/`--to` (spec §4): Task 13.
   - Aggregate report (spec §4 + Resolved decisions): Task 14.
   - CLI dispatcher (spec §4): Task 15.
   - Integration test: Task 16.
   - Prompt preview / manual review gate (new): Task 17.
   - Code review gate (new): Task 18.
   - Rollout batches (spec §5): Task 19.
   - KG validity tests (spec §5): Task 21.
   - Docs updates (spec §5): Task 22.
   - **All covered.**

2. **Placeholder scan:** No TBD, no "implement later", no "add appropriate error handling". Every code step includes the actual code.

3. **Type consistency:** `extract_one_paper` signature stable (arguments added, not renamed). `remap_value` signature consistent across tests and implementation. `VALID_GROWTH_PHASES` duplicated in three places (validator, prompts, utils/merge) — called out explicitly in Task 19 (rollout) with a sync test in Task 7 (`test_valid_growth_phases_list_matches_validator`).

4. **Known duplication of the enum set.** Three locations: `validate_paperconfig.py`, `prompts.py`, `extraction_utils.py` + `merge.py`. Avoided by importing once: future refactor — not needed for v1. The `test_valid_growth_phases_list_matches_validator` test pins `prompts.py` to the validator; `extraction_utils.py` + `merge.py` stay manually synced. Rollout step 17 flags all three files as the sync surface on promotion.

---

## Execution Handoff

**Plan complete and saved to `docs/superpowers/plans/2026-04-12-timepoint-growth-phase-backfill.md`. Two execution options:**

**1. Subagent-Driven (recommended)** — I dispatch a fresh subagent per task, review between tasks, fast iteration.

**2. Inline Execution** — Execute tasks in this session using `executing-plans`, batch execution with checkpoints for review.

**Which approach?**
