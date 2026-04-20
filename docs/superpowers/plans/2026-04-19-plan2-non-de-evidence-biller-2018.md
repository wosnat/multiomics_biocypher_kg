# Non-DE Evidence — Biller 2018 Slice, Plan 2 (Adapter + End-to-End Build)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the `ObservationsAdapter` that turns Plan 1's `derived_metrics_table` paperconfig entries into `DerivedMetric` nodes + 3 binding edges + `derived_metric_{flags,classifies,quantifies}_gene` measurement edges, emit `Experiment.compartment` from the omics adapter, wire the new adapter into the build pipeline, and add a synthetic numeric-DM fixture that exercises the numeric code path (which Biller 2018 itself does not hit).

**Architecture:** `observations_adapter.py` is a sibling of `cluster_adapter.py` and follows its exact structural pattern: a per-paperconfig `ObservationsAdapter` class with `get_nodes()` / `get_edges()` / `download_data()`, and a `MultiObservationsAdapter` wrapper that reads `paperconfig_files.txt`, builds `_organism_lookup` from `cyanobacteria_genomes.csv`, and delegates. Denormalized fields on every emitted `DerivedMetric` are computed from the **parent Experiment dict** (looked up via `paperconfig_utils.get_experiments(config)[entry['experiment']]`) — never re-read from the paperconfig entry. The registry from `multiomics_kg/vocab/non_de_evidence.py` is read **only indirectly** (the paperconfig already passed validation that matched registry value_kind); the adapter reads `rankable` / `has_p_value` / `allowed_categories` / `p_value_threshold` / `unit` / `field_description` from the paperconfig metric entry for numeric, and sets `rankable="false"` / `has_p_value="false"` / `allowed_categories=[]` / `unit=""` deterministically for boolean + categorical. Post-import rollups (`total_gene_count`, `growth_phases`, rank/percentile/bucket, significance) are Plan 3 — Plan 2 leaves those properties unset.

**Tech Stack:** Python 3.11+, BioCypher, pandas (CSV reading), pytest (unit tests), Docker Compose (import smoke test).

---

## File Structure

**Create:**
- `multiomics_kg/adapters/observations_adapter.py` — `ObservationsAdapter` + `MultiObservationsAdapter` + helpers (`_clean_str`, `_make_derived_metric_id`, `_resolve_csv_path`, `_parse_boolean_cell`, `_load_pdf_cache`).
- `tests/test_observations_adapter.py` — unit tests for the adapter (all 3 value_kinds, token parser, denormalization, edge skips, organism lookup).
- `tests/fixtures/non_de/__init__.py` — empty marker.
- `tests/fixtures/non_de/synthetic_numeric_metrics.csv` — 100-row numeric-DM test fixture (three metrics, real MED4 locus tags, p-values).
- `tests/fixtures/non_de/synthetic_paperconfig.yaml` — paperconfig referencing the fixture CSV.
- `tests/fixtures/non_de/paperconfig_files.txt` — one-line list pointing at the synthetic paperconfig (test-only, not wired into the production build).

**Modify:**
- `multiomics_kg/adapters/omics_adapter.py` — add `compartment` to the Experiment node properties (default `"whole_cell"`).
- `create_knowledge_graph.py` — instantiate `MultiObservationsAdapter` after `MultiClusterAdapter` and call `bc.write_nodes()` / `bc.write_edges()`.

**No changes:** schema_config.yaml (Plan 1), paperconfig_utils.py (Plan 1), validate_paperconfig.py (Plan 1), build_gene_id_mapping.py (Plan 1), resolve_paper_ids.py (Plan 1), post-import scripts (Plan 3), KG validity tests (Plan 3), docs (Plan 3).

---

## Conventions reused verbatim from `cluster_adapter.py`

These are not negotiable — the code review agent will reject divergence without justification:

- Pub node ID: `doi:{doi}`
- Experiment node ID: `{doi}_{experiment_key}` — **no** `doi:` prefix on experiments
- DerivedMetric node ID: `derived_metric:{doi_short}:{entry_key}:{metric_type}` where `doi_short = doi.rstrip("/").rsplit("/", 1)[-1]` (or paper-name slug fallback)
- OrganismTaxon ID lookup: `_organism_lookup[preferred_name]` → `insdc.gcf:{accession}`, populated from `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv`
- CSV resolution: probe for `<stem>_resolved.csv` via `resolve_paper_ids.get_resolved_path`; prefer column `resolved_locus_tag` over `name_col` when the resolved file exists
- String sanitization: `_clean_str(value) → value.replace("'", "^").replace("|", ",")` — **local helper** per CLAUDE.md, not an import
- PDF cache path: `cache/pdf_extraction_cache.json` (relative to project root)
- Edge ID convention for binding edges: `pub_dm__{dm_id}`, `exp_dm__{dm_id}__{exp_key}`, `dm_org__{dm_id}`; measurement edges: `{dm_id}__{gene_locus}` (each DM emits only one edge type, so collisions are impossible)

---

# Tasks

### Task 0: Polish Biller 2018 paperconfig as the template for future derived-metric paperconfigs

The Biller 2018 paperconfig (authored in Plan 1) is the first real-paper use of `derived_metrics_table`, so its author-chosen strings become the template for future paperconfigs (zinser 2009, Waldbauer 2012, biller 2022). Two improvements land in one commit:

**A. Shorten `entry_key` names.**
Current keys like `natl2a_periodicity_axenic_derived_metrics` are verbose — the `_derived_metrics` suffix duplicates `type: derived_metrics_table`, and the qualifiers pack redundantly. They flow into DerivedMetric node IDs as-is (`derived_metric:{doi_short}:{entry_key}:{metric_type}`), producing ~94-char IDs that are longer than necessary for logs/queries/debugging. Shorter keys produce ~70-char IDs in line with cluster_adapter's existing output. **Convention this establishes for future paperconfigs:** `{supp_table_ref}_{strain}_{qualifier?}` — compact, self-documenting, uses the paper's own table labels (S4A/S4B/S5) as the leading disambiguator.

**B. Add `name:` to every metric entry.**
ObservationsAdapter's `name` fallback is bare `metric_type` — valid but unhelpful in UI and indistinguishable from `metric_type` itself. Explicit `name:` makes DerivedMetric.name human-readable and disambiguates the two DMs that share `metric_type=periodic_in_coculture_LD` (NATL2A coculture entry + MIT1002 entry).

Both changes are local to this one paperconfig file. Entry keys are only referenced as map keys in `supplementary_materials` and as a component of DM node IDs (only produced, never consumed) — no cross-file dependencies.

**Files:**
- Modify: `data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml`

- [ ] **Step 1: Rename the 4 entry keys**

In `data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml`, rename each `supplementary_materials` key:

| Current | New |
|---|---|
| `natl2a_periodicity_axenic_derived_metrics:` | `s4a_natl2a_axenic:` |
| `natl2a_periodicity_coculture_derived_metrics:` | `s4a_natl2a_coculture:` |
| `mit1002_periodicity_derived_metrics:` | `s4b_mit1002:` |
| `natl2a_darkness_survival_derived_metrics:` | `s5_natl2a_survival:` |

Only the key name changes — all nested fields (filename, organism, experiment, metrics, ...) stay identical.

- [ ] **Step 2: Add `name:` to every metric entry**

Insert a `name:` line **immediately after** each metric's `metric_type:` line (preserving canonical ordering: metric_type → name → value_kind → value_col → ...).

**Entry `s4a_natl2a_axenic` (2 metrics):**

```yaml
      - metric_type: periodic_in_axenic_LD
        name: "Periodic in NATL2A axenic L:D (Table S4A)"
        value_kind: boolean
        value_col: "Periodic in axenic, L:D cultures"
        # ... rest unchanged
      - metric_type: periodic_in_axenic_extended_darkness
        name: "Periodic in NATL2A axenic extended darkness (Table S4A)"
        value_kind: boolean
        # ... rest unchanged
```

**Entry `s4a_natl2a_coculture` (2 metrics):**

```yaml
      - metric_type: periodic_in_coculture_LD
        name: "Periodic in NATL2A coculture L:D (Table S4A)"
        value_kind: boolean
        # ... rest unchanged
      - metric_type: periodic_in_coculture_extended_darkness
        name: "Periodic in NATL2A coculture extended darkness (Table S4A)"
        value_kind: boolean
        # ... rest unchanged
```

**Entry `s4b_mit1002` (2 metrics):**

```yaml
      - metric_type: periodic_in_coculture_LD
        name: "Periodic in MIT1002 coculture L:D (Table S4B)"
        value_kind: boolean
        # ... rest unchanged
      - metric_type: periodic_in_coculture_extended_darkness
        name: "Periodic in MIT1002 coculture extended darkness (Table S4B)"
        value_kind: boolean
        # ... rest unchanged
```

**Entry `s5_natl2a_survival` (1 metric):**

```yaml
      - metric_type: darkness_survival_class
        name: "NATL2A darkness survival class (Table S5)"
        value_kind: categorical
        # ... rest unchanged
```

- [ ] **Step 3: Validate the paperconfig**

Run: `uv run python scripts/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml"`
Expected: exit 0. The `name` field is an optional free-text property on metric entries — validator should accept silently, same as it does for `field_description`.

> If the validator rejects `name:` as unknown, Plan 1's validator extension missed it. Add `name` to the per-metric allowed-fields list in `scripts/validate_paperconfig.py` (inside the `derived_metrics_table` dispatch), alongside `field_description`. That change belongs in this same commit.

- [ ] **Step 4: Sanity-check the resulting DM IDs by tracing the format**

Expected DM node IDs after this task (what Tasks 5, 8, 15 will produce):

```
derived_metric:mSystems.00040-18:s4a_natl2a_axenic:periodic_in_axenic_LD
derived_metric:mSystems.00040-18:s4a_natl2a_axenic:periodic_in_axenic_extended_darkness
derived_metric:mSystems.00040-18:s4a_natl2a_coculture:periodic_in_coculture_LD
derived_metric:mSystems.00040-18:s4a_natl2a_coculture:periodic_in_coculture_extended_darkness
derived_metric:mSystems.00040-18:s4b_mit1002:periodic_in_coculture_LD
derived_metric:mSystems.00040-18:s4b_mit1002:periodic_in_coculture_extended_darkness
derived_metric:mSystems.00040-18:s5_natl2a_survival:darkness_survival_class
```

All 7 IDs are unique. The two `periodic_in_coculture_LD` DMs are disambiguated by `s4a_natl2a_coculture` vs `s4b_mit1002`. Length range: 68–85 chars (vs ~94 chars with original keys).

- [ ] **Step 5: Commit**

```bash
git add "data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml"
# If the validator needed an update, also:
# git add scripts/validate_paperconfig.py
git commit -m "Biller 2018 paperconfig: shorten entry_keys + add name: on every metric (template)"
```

---

### Task 1: Emit `Experiment.compartment` from omics_adapter (default `"whole_cell"`)

Plan 1 added `compartment: str` to the Experiment schema block in `config/schema_config.yaml` but no adapter yet emits it. ObservationsAdapter reads `compartment` back off the parent Experiment **dict** (not the node), so this task is independent — but the field must still reach the graph so queries like `MATCH (e:Experiment {compartment: 'vesicle'})` work. Biller 2018's experiments have no explicit `compartment` in the paperconfig, so all three will default to `"whole_cell"`.

**Files:**
- Modify: `multiomics_kg/adapters/omics_adapter.py` — add one line to the `exp_props` dict inside `get_nodes()`.
- Test: `tests/test_omics_adapter_experiment.py` (create).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_omics_adapter_experiment.py
"""Tests for Experiment.compartment emission (Plan 2 Task 1)."""
import os
import tempfile

import pytest
import yaml

from multiomics_kg.adapters.omics_adapter import OMICSAdapter


@pytest.fixture
def _stub_pdf_extractor(monkeypatch):
    """PDFPublicationExtractor hits disk/network; stub it out."""
    def _fake_extract(self, pdf_path):
        return {"publication": {"title": "stub", "doi": None}}
    from multiomics_kg.adapters.pdf_publication_extraction import PDFPublicationExtractor
    monkeypatch.setattr(PDFPublicationExtractor, "extract_from_pdf", _fake_extract)


def _write_paperconfig(tmp_path, compartment=None):
    exp = {
        "name": "Test exp",
        "organism": "Prochlorococcus MED4",
        "treatment_condition": "N-limit",
        "control_condition": "replete",
        "omics_type": "RNASEQ",
        "test_type": "DESeq2",
        "treatment_type": ["nitrogen"],
        "background_factors": [],
    }
    if compartment is not None:
        exp["compartment"] = compartment
    config = {
        "publication": {
            "papername": "Test 2024",
            "doi": "10.1234/test.2024",
            "papermainpdf": str(tmp_path / "fake.pdf"),
            "experiments": {"exp_a": exp},
            "supplementary_materials": {},
        }
    }
    config_path = tmp_path / "paperconfig.yaml"
    config_path.write_text(yaml.dump(config))
    (tmp_path / "fake.pdf").write_bytes(b"")
    return str(config_path)


def test_experiment_compartment_defaults_to_whole_cell(tmp_path, _stub_pdf_extractor):
    config_path = _write_paperconfig(tmp_path)  # no compartment in paperconfig
    adapter = OMICSAdapter(config_file=config_path)
    adapter.download_data()
    nodes = adapter.get_nodes()
    exp_nodes = [(nid, lbl, props) for nid, lbl, props in nodes if lbl == "experiment"]
    assert len(exp_nodes) == 1
    _, _, props = exp_nodes[0]
    assert props["compartment"] == "whole_cell"


def test_experiment_compartment_honours_paperconfig(tmp_path, _stub_pdf_extractor):
    config_path = _write_paperconfig(tmp_path, compartment="vesicle")
    adapter = OMICSAdapter(config_file=config_path)
    adapter.download_data()
    nodes = adapter.get_nodes()
    exp_nodes = [(nid, lbl, props) for nid, lbl, props in nodes if lbl == "experiment"]
    _, _, props = exp_nodes[0]
    assert props["compartment"] == "vesicle"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_omics_adapter_experiment.py -v`
Expected: FAIL with `KeyError: 'compartment'` (the Experiment node dict doesn't contain the key).

- [ ] **Step 3: Implement**

In `multiomics_kg/adapters/omics_adapter.py`, inside `get_nodes()` where `exp_props` is built (look for the existing dict that sets `"name"`, `"organism_name"`, etc.), add **one** line. Locate this block:

```python
            exp_props = {
                    "name": self.clean_text(exp.get("name", "")),
                    "organism_name": self.clean_text(exp.get("organism", "")),
                    "treatment_type": self._normalize_list_field(exp, "treatment_type"),
                    "treatment": self.clean_text(exp.get("treatment_condition", "")),
                    "control": self.clean_text(exp.get("control_condition", "")),
                    "experimental_context": self.clean_text(exp.get("experimental_context", "")),
                    "omics_type": self.clean_text(exp.get("omics_type", "")),
                    "statistical_test": self.clean_text(exp.get("test_type", "")),
                    "is_time_course": "true" if len(timepoints) > 1 else "false",
                    "medium": self.clean_text(exp.get("medium", "")),
                    "temperature": self.clean_text(exp.get("temperature", "")),
                    "light_condition": self.clean_text(exp.get("light_condition", "")),
                    "light_intensity": self.clean_text(exp.get("light_intensity", "")),
                    "table_scope": self.clean_text(exp.get("table_scope", "")),
                    "table_scope_detail": self.clean_text(exp.get("table_scope_detail", "")),
                    "background_factors": self._normalize_list_field(exp, "background_factors"),
            }
```

Add this line **between `"organism_name"` and `"treatment_type"`** so compartment sits near the other denormalized metadata:

```python
                    "compartment": self.clean_text(exp.get("compartment", "whole_cell")),
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_omics_adapter_experiment.py -v`
Expected: PASS (2/2).

- [ ] **Step 5: Run the full adapter test suite to catch regressions**

Run: `uv run pytest tests/test_omics_adapter_organism_gene.py tests/test_omics_adapter_experiment.py -v`
Expected: all pass. The added property is additive — no existing test should break.

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/adapters/omics_adapter.py tests/test_omics_adapter_experiment.py
git commit -m "omics_adapter: emit Experiment.compartment (default whole_cell)"
```

---

### Task 2: Scaffold `observations_adapter.py` — module skeleton + node-id helpers

Stand up the module with the structural pieces all subsequent tasks depend on: the `_clean_str` helper, `_make_derived_metric_id`, `_resolve_csv_path`, `_load_pdf_cache`. No `get_nodes` / `get_edges` logic yet — those land in Tasks 5–11.

**Files:**
- Create: `multiomics_kg/adapters/observations_adapter.py`
- Test: `tests/test_observations_adapter.py` (create)

- [ ] **Step 1: Write the failing test**

```python
# tests/test_observations_adapter.py
"""Tests for ObservationsAdapter (Plan 2)."""
import os
from pathlib import Path

import pandas as pd
import pytest
import yaml

from multiomics_kg.adapters.observations_adapter import (
    ObservationsAdapter,
    MultiObservationsAdapter,
    _clean_str,
    _make_derived_metric_id,
    _resolve_csv_path,
    _parse_boolean_cell,
)


def test_clean_str_replaces_single_quote_and_pipe():
    assert _clean_str("it's | great") == "it^s , great"


def test_clean_str_handles_none():
    assert _clean_str(None) == ""


def test_clean_str_passes_non_string_through():
    # ints stringify but don't go through replacement
    assert _clean_str(42) == "42"


def test_make_derived_metric_id_uses_doi_short():
    dm_id = _make_derived_metric_id(
        doi="10.1128/mSystems.00040-18",
        paper_name="Biller 2018",
        entry_key="s4a_natl2a_ld",
        metric_type="periodic_in_axenic_LD",
    )
    assert dm_id == "derived_metric:mSystems.00040-18:s4a_natl2a_ld:periodic_in_axenic_LD"


def test_make_derived_metric_id_falls_back_to_paper_slug():
    dm_id = _make_derived_metric_id(
        doi="",
        paper_name="Biller 2018",
        entry_key="s4a_natl2a_ld",
        metric_type="periodic_in_axenic_LD",
    )
    assert dm_id == "derived_metric:biller_2018:s4a_natl2a_ld:periodic_in_axenic_LD"


def test_resolve_csv_path_prefers_resolved(tmp_path):
    src = tmp_path / "table.csv"
    resolved = tmp_path / "table_resolved.csv"
    src.write_text("a,b\n1,2\n")
    resolved.write_text("a,b,resolved_locus_tag\n1,2,PMM0001\n")
    path, used_resolved = _resolve_csv_path(str(src))
    assert path == resolved
    assert used_resolved is True


def test_resolve_csv_path_falls_back_to_original(tmp_path):
    src = tmp_path / "table.csv"
    src.write_text("a,b\n1,2\n")
    path, used_resolved = _resolve_csv_path(str(src))
    assert path == src
    assert used_resolved is False
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_observations_adapter.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'multiomics_kg.adapters.observations_adapter'`.

- [ ] **Step 3: Write the minimal implementation**

```python
# multiomics_kg/adapters/observations_adapter.py
"""
ObservationsAdapter — reads derived_metrics_table entries from paperconfig.yaml
files and emits:
  - DerivedMetric nodes (one per Experiment × metric_type)
  - publication_has_derived_metric edges (Publication → DerivedMetric)
  - experiment_has_derived_metric edges (Experiment → DerivedMetric)
  - derived_metric_belongs_to_organism edges (DerivedMetric → OrganismTaxon)
  - derived_metric_quantifies_gene edges (DerivedMetric → Gene, value_kind=numeric)
  - derived_metric_flags_gene edges (DerivedMetric → Gene, value_kind=boolean)
  - derived_metric_classifies_gene edges (DerivedMetric → Gene, value_kind=categorical)

Sibling of cluster_adapter.py — follows the same structural pattern.
"""
import json
import logging
import re
from pathlib import Path

import pandas as pd

from multiomics_kg.download.resolve_paper_ids import get_resolved_path
from multiomics_kg.utils.paperconfig_utils import (
    load_paperconfig,
    load_all_paperconfigs,
    get_paper_name,
    get_experiments,
    iter_derived_metrics_tables,
)
from multiomics_kg.vocab.non_de_evidence import DEFAULT_SKIP_TOKENS, VALID_BLANK_POLICIES

logger = logging.getLogger(__name__)

_DEFAULT_PDF_CACHE = Path(__file__).parent.parent.parent / "cache" / "pdf_extraction_cache.json"


def _clean_str(value) -> str:
    """Sanitize string for BioCypher CSV output (CLAUDE.md convention)."""
    if value is None:
        return ""
    if not isinstance(value, str):
        return str(value)
    return value.replace("'", "^").replace("|", ",")


def _make_derived_metric_id(
    doi: str, paper_name: str, entry_key: str, metric_type: str
) -> str:
    """DerivedMetric node ID: derived_metric:{doi_short}:{entry_key}:{metric_type}."""
    if doi:
        doi_short = doi.rstrip("/").rsplit("/", 1)[-1]
    else:
        doi_short = re.sub(r"[^a-z0-9]+", "_", paper_name.lower()).strip("_")
    return f"derived_metric:{doi_short}:{entry_key}:{metric_type}"


def _resolve_csv_path(csv_path: str) -> tuple[Path, bool]:
    """Probe for pre-resolved CSV (written by resolve_paper_ids.py).

    Returns (path_to_use, is_resolved).
    """
    p = Path(csv_path)
    resolved = get_resolved_path(p)
    if resolved.exists():
        return resolved, True
    return p, False


def _load_pdf_cache(cache_path: Path = _DEFAULT_PDF_CACHE) -> dict:
    if cache_path.exists():
        try:
            with open(cache_path) as f:
                return json.load(f)
        except Exception:
            pass
    return {}


def _parse_boolean_cell(
    value,
    true_tokens: list[str],
    false_tokens: list[str],
    skip_tokens: list[str],
    blank_policy: str,
) -> str | None:
    """Map a single CSV cell value to "true" / "false" / None (=skip).

    Hard-errors on unexpected tokens — per parent spec, the adapter must not
    silently coerce ambiguous values.
    """
    # NaN / None / empty string → apply blank_policy
    if value is None or pd.isna(value):
        return _apply_blank_policy(blank_policy)
    s = str(value).strip()
    if s == "":
        return _apply_blank_policy(blank_policy)
    if s in true_tokens:
        return "true"
    if s in false_tokens:
        return "false"
    if s in skip_tokens:
        return None
    raise ValueError(
        f"Unexpected boolean token {s!r}: not in true_tokens={true_tokens}, "
        f"false_tokens={false_tokens}, or skip_tokens={skip_tokens}. "
        f"Paperconfig author must classify explicitly."
    )


def _apply_blank_policy(blank_policy: str) -> str | None:
    if blank_policy == "skip":
        return None
    if blank_policy == "true":
        return "true"
    if blank_policy == "false":
        return "false"
    raise ValueError(
        f"Invalid blank_policy {blank_policy!r}; must be one of {VALID_BLANK_POLICIES}"
    )


# ─── Adapter classes (get_nodes / get_edges land in Tasks 5–11) ─────────────


class ObservationsAdapter:
    """Adapter for one paperconfig's derived_metrics_table entries."""

    def __init__(self, config_file: str, test_mode: bool = False):
        self.config_file = config_file
        self.test_mode = test_mode
        self.config = load_paperconfig(Path(config_file))
        self.paper_name = get_paper_name(self.config, fallback_path=Path(config_file))
        self.doi = self._extract_doi()
        self._dm_entries = list(iter_derived_metrics_tables(self.config))
        self._organism_lookup: dict[str, str] = {}

    def _extract_doi(self) -> str:
        pub = self.config.get("publication", {})
        doi = pub.get("doi", "")
        if doi:
            return doi
        pdf_path = pub.get("papermainpdf", "")
        if pdf_path:
            cache = _load_pdf_cache()
            cached = cache.get(pdf_path, {})
            doi = cached.get("publication", {}).get("doi", "")
        return doi or ""

    def get_nodes(self) -> list[tuple]:
        return []  # filled in Tasks 5–7

    def get_edges(self) -> list[tuple]:
        return []  # filled in Tasks 8–11

    def download_data(self, **kwargs):
        pass


class MultiObservationsAdapter:
    """Wrapper that reads paperconfig_files.txt and delegates."""

    def __init__(
        self,
        config_list_file: str | list[str],
        genome_config_file: str = None,
        test_mode: bool = False,
        **kwargs,
    ):
        self._organism_lookup: dict[str, str] = {}
        if genome_config_file:
            self._organism_lookup = self._build_organism_lookup(genome_config_file)
        self.adapters: list[ObservationsAdapter] = []
        list_files = config_list_file if isinstance(config_list_file, list) else [config_list_file]
        paperconfigs = load_all_paperconfigs([Path(lf) for lf in list_files])
        for pc_path, config in paperconfigs:
            supp = config.get("publication", {}).get("supplementary_materials", {})
            has_dm = any(
                isinstance(v, dict) and v.get("type") == "derived_metrics_table"
                for v in supp.values()
            )
            if not has_dm:
                continue
            adapter = ObservationsAdapter(config_file=str(pc_path), test_mode=test_mode)
            adapter._organism_lookup = self._organism_lookup
            self.adapters.append(adapter)

    def _build_organism_lookup(self, genome_config_file: str) -> dict[str, str]:
        lookup: dict[str, str] = {}
        try:
            df = pd.read_csv(genome_config_file)
            for _, row in df.iterrows():
                name = row.get("preferred_name", "")
                accession = row.get("ncbi_accession", "")
                if name and accession:
                    lookup[str(name)] = f"insdc.gcf:{accession}"
        except Exception as e:
            logger.warning(f"Could not load genome config '{genome_config_file}': {e}")
        return lookup

    def download_data(self, **kwargs):
        pass

    def get_nodes(self) -> list[tuple]:
        nodes = []
        for adapter in self.adapters:
            nodes.extend(adapter.get_nodes())
        return nodes

    def get_edges(self) -> list[tuple]:
        edges = []
        for adapter in self.adapters:
            edges.extend(adapter.get_edges())
        return edges
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_observations_adapter.py -v`
Expected: 7/7 pass (the six `test_*` functions above all defined in Step 1).

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/observations_adapter.py tests/test_observations_adapter.py
git commit -m "observations_adapter: scaffold module + id helpers + boolean parser"
```

---

### Task 3: Boolean token parser — exhaustive coverage

`_parse_boolean_cell` was defined in Task 2 (so Task 5 onward can call it). Now add the full test matrix: all 4 legal paths + 3 blank-policy variants + unknown-token error + skip-token path.

**Files:**
- Modify: `tests/test_observations_adapter.py` (append test cases)

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_observations_adapter.py`:

```python
def _bp(value, blank_policy="skip"):
    return _parse_boolean_cell(
        value,
        true_tokens=["Y", "yes"],
        false_tokens=["N", "no"],
        skip_tokens=["NA", "n/a"],
        blank_policy=blank_policy,
    )


def test_boolean_true_token():
    assert _bp("Y") == "true"
    assert _bp("yes") == "true"


def test_boolean_false_token():
    assert _bp("N") == "false"
    assert _bp("no") == "false"


def test_boolean_skip_token_returns_none():
    assert _bp("NA") is None
    assert _bp("n/a") is None


def test_boolean_blank_with_skip_policy():
    assert _bp("", blank_policy="skip") is None
    assert _bp(None, blank_policy="skip") is None
    assert _bp(float("nan"), blank_policy="skip") is None


def test_boolean_blank_with_true_policy():
    assert _bp("", blank_policy="true") == "true"
    assert _bp(None, blank_policy="true") == "true"


def test_boolean_blank_with_false_policy():
    assert _bp("", blank_policy="false") == "false"
    assert _bp(float("nan"), blank_policy="false") == "false"


def test_boolean_unknown_token_raises():
    with pytest.raises(ValueError, match="Unexpected boolean token"):
        _bp("maybe")


def test_boolean_invalid_blank_policy_raises():
    with pytest.raises(ValueError, match="Invalid blank_policy"):
        _parse_boolean_cell("", true_tokens=["Y"], false_tokens=[], skip_tokens=[], blank_policy="whatever")


def test_boolean_whitespace_treated_as_blank():
    # Empty-after-strip → blank → skip_policy default returns None
    assert _bp("   ") is None
    # With blank_policy="false", whitespace is blank → "false"
    assert _bp("   ", blank_policy="false") == "false"
```

- [ ] **Step 2: Run tests to verify they all pass**

Run: `uv run pytest tests/test_observations_adapter.py -v`
Expected: all new tests PASS (parser implementation landed in Task 2 already).

- [ ] **Step 3: Commit**

```bash
git add tests/test_observations_adapter.py
git commit -m "observations_adapter: exhaustive boolean parser tests"
```

---

### Task 4: DOI extraction + denormalization helper

Extract the parent-Experiment denormalization logic into a dedicated helper method `_denormalized_fields(experiment_dict) → dict` so the three `get_nodes`-branch tasks (5–7) share one codepath and can't drift. This is the enforcement point for spec invariant #9 ("denormalized fields match parent").

**Files:**
- Modify: `multiomics_kg/adapters/observations_adapter.py` — add `_denormalized_fields` + tests

- [ ] **Step 1: Write the failing test**

Append to `tests/test_observations_adapter.py`:

```python
def test_denormalized_fields_from_experiment():
    """_denormalized_fields produces the exact 9-field dict DerivedMetric.* needs."""
    # Arrange: a minimal paperconfig + ObservationsAdapter instance
    # Easiest: build the adapter from an in-memory config dict
    from multiomics_kg.adapters.observations_adapter import ObservationsAdapter
    exp = {
        "name": "NATL2A extended darkness",
        "organism": "Prochlorococcus NATL2A",
        "omics_type": "RNASEQ",
        "treatment_type": ["darkness"],
        "background_factors": ["axenic", "diel"],
        "treatment_condition": "Extended darkness",
        "light_condition": "continuous darkness",
        "experimental_context": "Axenic NATL2A in Pro99 at 24C",
        "compartment": "whole_cell",
    }
    # Bypass __init__'s file I/O by constructing manually
    adapter = ObservationsAdapter.__new__(ObservationsAdapter)
    adapter.doi = "10.1128/mSystems.00040-18"
    fields = adapter._denormalized_fields(exp)

    assert fields == {
        "organism_name": "Prochlorococcus NATL2A",
        "publication_doi": "10.1128/mSystems.00040-18",
        "compartment": "whole_cell",
        "omics_type": "RNASEQ",
        "treatment_type": ["darkness"],
        "background_factors": ["axenic", "diel"],
        "treatment": "Extended darkness",
        "light_condition": "continuous darkness",
        "experimental_context": "Axenic NATL2A in Pro99 at 24C",
    }


def test_denormalized_fields_compartment_defaults():
    from multiomics_kg.adapters.observations_adapter import ObservationsAdapter
    exp = {
        "name": "x", "organism": "Prochlorococcus MED4",
        "omics_type": "RNASEQ", "treatment_type": [],
        "background_factors": [], "treatment_condition": "",
        "light_condition": "", "experimental_context": "",
        # no compartment key
    }
    adapter = ObservationsAdapter.__new__(ObservationsAdapter)
    adapter.doi = "10.1234/x"
    fields = adapter._denormalized_fields(exp)
    assert fields["compartment"] == "whole_cell"


def test_denormalized_fields_normalizes_scalar_list_fields():
    """paperconfig may accidentally provide a string where list is expected."""
    from multiomics_kg.adapters.observations_adapter import ObservationsAdapter
    exp = {
        "organism": "P. MED4", "omics_type": "RNASEQ",
        "treatment_type": "nitrogen",  # scalar, should become ["nitrogen"]
        "background_factors": "axenic",  # scalar → list
        "treatment_condition": "", "light_condition": "",
        "experimental_context": "", "compartment": "whole_cell",
    }
    adapter = ObservationsAdapter.__new__(ObservationsAdapter)
    adapter.doi = "10.1234/x"
    fields = adapter._denormalized_fields(exp)
    assert fields["treatment_type"] == ["nitrogen"]
    assert fields["background_factors"] == ["axenic"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_observations_adapter.py -v -k denormalized`
Expected: FAIL with `AttributeError: 'ObservationsAdapter' object has no attribute '_denormalized_fields'`.

- [ ] **Step 3: Implement `_denormalized_fields`**

In `multiomics_kg/adapters/observations_adapter.py`, inside class `ObservationsAdapter`, add (place after `_extract_doi`, before `get_nodes`):

```python
    def _denormalized_fields(self, experiment: dict) -> dict:
        """Derive the 9 parent-Experiment fields every DerivedMetric copies.

        Enforces spec invariant #9: DerivedMetric denormalized fields equal
        parent Experiment's values. Never re-read from the paperconfig entry.
        """
        def _as_list(v):
            if v is None:
                return []
            if isinstance(v, str):
                return [_clean_str(v)] if v else []
            return [_clean_str(x) for x in v]

        return {
            "organism_name": _clean_str(experiment.get("organism", "")),
            "publication_doi": _clean_str(self.doi),
            "compartment": _clean_str(experiment.get("compartment", "whole_cell")),
            "omics_type": _clean_str(experiment.get("omics_type", "")),
            "treatment_type": _as_list(experiment.get("treatment_type", [])),
            "background_factors": _as_list(experiment.get("background_factors", [])),
            "treatment": _clean_str(experiment.get("treatment_condition", "")),
            "light_condition": _clean_str(experiment.get("light_condition", "")),
            "experimental_context": _clean_str(experiment.get("experimental_context", "")),
        }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_observations_adapter.py -v -k denormalized`
Expected: 3/3 PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/observations_adapter.py tests/test_observations_adapter.py
git commit -m "observations_adapter: _denormalized_fields helper (spec invariant #9)"
```

---

### Task 5: `get_nodes` — boolean DerivedMetric emission

Implement the boolean branch of `get_nodes()`. Each boolean metric produces one `DerivedMetric` node with `value_kind="boolean"`, `rankable="false"`, `has_p_value="false"`, empty `allowed_categories`, empty `unit`, null `p_value_threshold` (emitted as `0.0` sentinel or schema allows null — we emit `0.0` and document; alternatively leave the key out so it's null in Neo4j). **Decision:** leave `p_value_threshold` out of the props dict for boolean/categorical — BioCypher will emit it as null. This matches the spec ("null when no p-value or no threshold reported").

**Files:**
- Modify: `multiomics_kg/adapters/observations_adapter.py` — fill `get_nodes` (boolean branch only for now)
- Modify: `tests/test_observations_adapter.py` — add boolean-node test

- [ ] **Step 1: Write the failing test**

Append to `tests/test_observations_adapter.py`:

```python
def _write_s4a_like_paperconfig(tmp_path):
    """Minimal in-memory paperconfig mimicking Biller 2018 S4A-axenic shape."""
    csv_path = tmp_path / "s4a_small.csv"
    # 3 rows: 2 Y + 1 blank for axenic_LD; 1 Y + 2 blank for extended_darkness
    csv_path.write_text(
        'NCBI ID_2,"Periodic in axenic, L:D cultures","Periodic in axenic, extended darkness cultures"\n'
        "PMN2A_RS00015,Y,\n"
        "PMN2A_RS00020,Y,Y\n"
        "PMN2A_RS00025,,\n"
    )
    resolved_path = tmp_path / "s4a_small_resolved.csv"
    resolved_path.write_text(
        'NCBI ID_2,"Periodic in axenic, L:D cultures","Periodic in axenic, extended darkness cultures",resolved_locus_tag,resolution_method\n'
        "PMN2A_RS00015,Y,,PMN2A_1328,tier1:NCBI ID_2\n"
        "PMN2A_RS00020,Y,Y,PMN2A_1329,tier1:NCBI ID_2\n"
        "PMN2A_RS00025,,,,unresolved\n"
    )
    config = {
        "publication": {
            "papername": "Biller 2018",
            "doi": "10.1128/mSystems.00040-18",
            "papermainpdf": str(tmp_path / "fake.pdf"),
            "experiments": {
                "axenic_rnaseq": {
                    "name": "NATL2A axenic",
                    "organism": "Prochlorococcus NATL2A",
                    "omics_type": "RNASEQ",
                    "treatment_type": ["darkness"],
                    "background_factors": ["axenic", "diel"],
                    "treatment_condition": "Extended darkness",
                    "light_condition": "continuous darkness",
                    "experimental_context": "Axenic NATL2A in Pro99",
                },
            },
            "supplementary_materials": {
                "s4a_axenic": {
                    "type": "derived_metrics_table",
                    "filename": str(csv_path),
                    "organism": "Prochlorococcus NATL2A",
                    "experiment": "axenic_rnaseq",
                    "name_col": "NCBI ID_2",
                    "metrics": [
                        {
                            "metric_type": "periodic_in_axenic_LD",
                            "value_kind": "boolean",
                            "value_col": "Periodic in axenic, L:D cultures",
                            "true_tokens": ["Y"],
                            "false_tokens": [],
                            "skip_tokens": ["NA", "N/A"],
                            "blank_policy": "skip",
                            "field_description": "RAIN FDR<0.05 axenic L:D",
                        },
                        {
                            "metric_type": "periodic_in_axenic_extended_darkness",
                            "value_kind": "boolean",
                            "value_col": "Periodic in axenic, extended darkness cultures",
                            "true_tokens": ["Y"],
                            "false_tokens": [],
                            "skip_tokens": ["NA", "N/A"],
                            "blank_policy": "skip",
                            "field_description": "RAIN FDR<0.05 axenic extended darkness",
                        },
                    ],
                },
            },
        },
    }
    pc_path = tmp_path / "paperconfig.yaml"
    pc_path.write_text(yaml.dump(config))
    return str(pc_path)


def test_get_nodes_emits_one_dm_per_boolean_metric(tmp_path):
    pc_path = _write_s4a_like_paperconfig(tmp_path)
    adapter = ObservationsAdapter(config_file=pc_path)
    nodes = adapter.get_nodes()

    dm_nodes = [(nid, lbl, props) for nid, lbl, props in nodes if lbl == "derived_metric"]
    assert len(dm_nodes) == 2

    ids = {nid for nid, _, _ in dm_nodes}
    assert "derived_metric:mSystems.00040-18:s4a_axenic:periodic_in_axenic_LD" in ids
    assert "derived_metric:mSystems.00040-18:s4a_axenic:periodic_in_axenic_extended_darkness" in ids


def test_boolean_dm_node_has_expected_props(tmp_path):
    pc_path = _write_s4a_like_paperconfig(tmp_path)
    adapter = ObservationsAdapter(config_file=pc_path)
    nodes = adapter.get_nodes()
    dm = next(
        props for nid, lbl, props in nodes
        if lbl == "derived_metric" and nid.endswith("periodic_in_axenic_LD")
    )
    assert dm["metric_type"] == "periodic_in_axenic_LD"
    assert dm["value_kind"] == "boolean"
    assert dm["rankable"] == "false"
    assert dm["has_p_value"] == "false"
    assert dm["allowed_categories"] == []
    assert dm["unit"] == ""
    assert dm["field_description"] == "RAIN FDR<0.05 axenic L:D"
    # Denormalized from parent Experiment
    assert dm["organism_name"] == "Prochlorococcus NATL2A"
    assert dm["compartment"] == "whole_cell"  # default
    assert dm["omics_type"] == "RNASEQ"
    assert dm["treatment_type"] == ["darkness"]
    assert dm["background_factors"] == ["axenic", "diel"]
    assert dm["treatment"] == "Extended darkness"
    assert dm["light_condition"] == "continuous darkness"
    assert dm["experimental_context"] == "Axenic NATL2A in Pro99"
    assert dm["publication_doi"] == "10.1128/mSystems.00040-18"
    # experiment_id is the raw-doi + exp_key concatenation
    assert dm["experiment_id"] == "10.1128/mSystems.00040-18_axenic_rnaseq"
    assert dm["name"].startswith("periodic_in_axenic_LD")  # default name
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_observations_adapter.py -v -k get_nodes_emits_one_dm`
Expected: FAIL — `get_nodes` currently returns `[]`.

- [ ] **Step 3: Implement `get_nodes` (boolean branch + skeleton for numeric/categorical)**

Replace the placeholder `get_nodes` in `multiomics_kg/adapters/observations_adapter.py` with the full implementation. This covers boolean now; Task 6 and 7 add categorical and numeric (trivial additional branches — structure lands here):

```python
    def get_nodes(self) -> list[tuple]:
        """Emit DerivedMetric nodes — one per (entry × metric_type).

        Matches cluster_adapter's behaviour: if the source CSV is missing,
        skip the whole entry so we never emit dangling DM nodes with no
        measurement edges.
        """
        nodes = []
        experiments = get_experiments(self.config)
        for entry_key, entry in self._dm_entries:
            exp_key = entry.get("experiment")
            if not exp_key or exp_key not in experiments:
                logger.warning(
                    f"derived_metrics_table '{entry_key}' references unknown "
                    f"experiment '{exp_key}' — skipping"
                )
                continue
            csv_path, _ = _resolve_csv_path(entry["filename"])
            if not csv_path.exists():
                logger.warning(
                    f"derived_metrics_table CSV not found: {csv_path} — skipping"
                )
                continue
            exp = experiments[exp_key]
            denorm = self._denormalized_fields(exp)

            for metric in entry.get("metrics", []):
                metric_type = metric.get("metric_type", "")
                value_kind = metric.get("value_kind", "")
                if not metric_type or not value_kind:
                    logger.warning(
                        f"Metric in '{entry_key}' missing metric_type or value_kind — skipping"
                    )
                    continue

                dm_id = _make_derived_metric_id(
                    self.doi, self.paper_name, entry_key, metric_type
                )
                # experiment_id mirrors the omics_adapter's format: raw-doi + "_" + exp_key
                experiment_id = f"{self.doi}_{exp_key}" if self.doi else f"{self.paper_name}_{exp_key}"

                props = {
                    "name": _clean_str(metric.get("name", metric_type)),
                    "experiment_id": _clean_str(experiment_id),
                    "metric_type": _clean_str(metric_type),
                    "value_kind": _clean_str(value_kind),
                    "field_description": _clean_str(metric.get("field_description", "")),
                    **denorm,
                }

                if value_kind == "boolean":
                    props["rankable"] = "false"
                    props["has_p_value"] = "false"
                    props["unit"] = ""
                    props["allowed_categories"] = []
                elif value_kind == "categorical":
                    # Task 6 completes the categorical branch
                    props["rankable"] = "false"
                    props["has_p_value"] = "false"
                    props["unit"] = ""
                    props["allowed_categories"] = [
                        _clean_str(c) for c in metric.get("allowed_categories", [])
                    ]
                elif value_kind == "numeric":
                    # Task 7 completes the numeric branch
                    props["rankable"] = _clean_str(metric.get("rankable", "false"))
                    props["has_p_value"] = _clean_str(metric.get("has_p_value", "false"))
                    props["unit"] = _clean_str(metric.get("unit", ""))
                    props["allowed_categories"] = []
                    pvt = metric.get("p_value_threshold")
                    if pvt is not None:
                        props["p_value_threshold"] = float(pvt)
                else:
                    logger.warning(
                        f"Unknown value_kind '{value_kind}' in '{entry_key}/{metric_type}' — skipping"
                    )
                    continue

                nodes.append((dm_id, "derived_metric", props))
        return nodes
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_observations_adapter.py -v -k get_nodes_emits_one_dm`
Then: `uv run pytest tests/test_observations_adapter.py -v -k boolean_dm_node`
Expected: 2/2 PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/observations_adapter.py tests/test_observations_adapter.py
git commit -m "observations_adapter: get_nodes emits boolean DerivedMetric"
```

---

### Task 6: `get_nodes` — categorical DerivedMetric coverage test

Task 5's implementation already covers categorical emission (the categorical branch is in the same switch). Add a dedicated unit test that exercises a categorical metric with `allowed_categories` populated from paperconfig (option C: inline, not from vocab).

**Files:**
- Modify: `tests/test_observations_adapter.py`

- [ ] **Step 1: Write the test**

Append to `tests/test_observations_adapter.py`:

```python
def _write_s5_like_paperconfig(tmp_path):
    """Minimal in-memory paperconfig mimicking Biller 2018 S5 categorical shape."""
    csv_path = tmp_path / "s5_small.csv"
    csv_path.write_text(
        "NCBI ID_2,darkness_cluster\n"
        "PMN2A_RS00030,darkness_axenic+darkness_coculture\n"
        "PMN2A_RS00040,darkness_coculture+unique_coculture\n"
        "PMN2A_RS00050,darkness_axenic+unique_axenic\n"
    )
    resolved_path = tmp_path / "s5_small_resolved.csv"
    resolved_path.write_text(
        "NCBI ID_2,darkness_cluster,resolved_locus_tag,resolution_method\n"
        "PMN2A_RS00030,darkness_axenic+darkness_coculture,PMN2A_1330,tier1:NCBI ID_2\n"
        "PMN2A_RS00040,darkness_coculture+unique_coculture,PMN2A_1331,tier1:NCBI ID_2\n"
        "PMN2A_RS00050,darkness_axenic+unique_axenic,PMN2A_1332,tier1:NCBI ID_2\n"
    )
    config = {
        "publication": {
            "papername": "Biller 2018",
            "doi": "10.1128/mSystems.00040-18",
            "papermainpdf": str(tmp_path / "fake.pdf"),
            "experiments": {
                "axenic_rnaseq": {
                    "name": "NATL2A axenic",
                    "organism": "Prochlorococcus NATL2A",
                    "omics_type": "RNASEQ",
                    "treatment_type": ["darkness"],
                    "background_factors": ["axenic", "diel"],
                    "treatment_condition": "Extended darkness",
                    "light_condition": "continuous darkness",
                    "experimental_context": "Axenic NATL2A in Pro99",
                },
            },
            "supplementary_materials": {
                "s5_survival": {
                    "type": "derived_metrics_table",
                    "filename": str(csv_path),
                    "organism": "Prochlorococcus NATL2A",
                    "experiment": "axenic_rnaseq",
                    "name_col": "NCBI ID_2",
                    "metrics": [
                        {
                            "metric_type": "darkness_survival_class",
                            "value_kind": "categorical",
                            "value_col": "darkness_cluster",
                            "allowed_categories": [
                                "darkness_axenic+darkness_coculture",
                                "darkness_coculture+unique_coculture",
                                "darkness_axenic+unique_axenic",
                            ],
                            "field_description": "Transcript presence at 72-144h",
                        },
                    ],
                },
            },
        },
    }
    pc_path = tmp_path / "paperconfig.yaml"
    pc_path.write_text(yaml.dump(config))
    return str(pc_path)


def test_categorical_dm_node_has_allowed_categories_from_paperconfig(tmp_path):
    pc_path = _write_s5_like_paperconfig(tmp_path)
    adapter = ObservationsAdapter(config_file=pc_path)
    nodes = adapter.get_nodes()
    dm_nodes = [(nid, lbl, props) for nid, lbl, props in nodes if lbl == "derived_metric"]
    assert len(dm_nodes) == 1
    _, _, props = dm_nodes[0]
    assert props["value_kind"] == "categorical"
    assert props["rankable"] == "false"
    assert props["has_p_value"] == "false"
    assert props["allowed_categories"] == [
        "darkness_axenic+darkness_coculture",
        "darkness_coculture+unique_coculture",
        "darkness_axenic+unique_axenic",
    ]
    assert props["unit"] == ""
    assert props["metric_type"] == "darkness_survival_class"
```

- [ ] **Step 2: Run test to verify it passes**

Run: `uv run pytest tests/test_observations_adapter.py -v -k categorical_dm_node`
Expected: PASS (the categorical branch of `get_nodes` was already written in Task 5).

- [ ] **Step 3: Commit**

```bash
git add tests/test_observations_adapter.py
git commit -m "observations_adapter: test categorical DerivedMetric node emission"
```

---

### Task 7: `get_nodes` — numeric DerivedMetric coverage test

Numeric metrics need their own test because they carry 3 additional props (`rankable`, `has_p_value`, `p_value_threshold`) that are paperconfig-declared, not adapter-forced.

**Files:**
- Modify: `tests/test_observations_adapter.py`

- [ ] **Step 1: Write the test**

Append to `tests/test_observations_adapter.py`:

```python
def _write_numeric_dm_paperconfig(tmp_path):
    csv_path = tmp_path / "numeric_small.csv"
    csv_path.write_text(
        "locus_tag,fourier,fourier_p\n"
        "PMM0001,0.87,0.001\n"
        "PMM0002,0.42,0.12\n"
    )
    # No _resolved.csv — name_col == locus_tag triggers the skip-resolve path in resolve_paper_ids
    config = {
        "publication": {
            "papername": "Synthetic Numeric",
            "doi": "10.9999/synthetic",
            "papermainpdf": str(tmp_path / "fake.pdf"),
            "experiments": {
                "diel_rnaseq": {
                    "name": "MED4 diel",
                    "organism": "Prochlorococcus MED4",
                    "omics_type": "RNASEQ",
                    "treatment_type": ["diel"],
                    "background_factors": [],
                    "treatment_condition": "diel cycle",
                    "light_condition": "14:10 LD",
                    "experimental_context": "MED4 in Pro99",
                },
            },
            "supplementary_materials": {
                "fourier_metrics": {
                    "type": "derived_metrics_table",
                    "filename": str(csv_path),
                    "organism": "Prochlorococcus MED4",
                    "experiment": "diel_rnaseq",
                    "name_col": "locus_tag",
                    "metrics": [
                        {
                            "metric_type": "fourier_score",
                            "value_kind": "numeric",
                            "value_col": "fourier",
                            "unit": "",
                            "rankable": "true",
                            "has_p_value": "true",
                            "p_value_col": "fourier_p",
                            "p_value_threshold": 0.05,
                            "field_description": "Fourier periodicity score",
                        },
                    ],
                },
            },
        },
    }
    pc_path = tmp_path / "paperconfig.yaml"
    pc_path.write_text(yaml.dump(config))
    return str(pc_path)


def test_numeric_dm_node_props_from_paperconfig(tmp_path):
    pc_path = _write_numeric_dm_paperconfig(tmp_path)
    adapter = ObservationsAdapter(config_file=pc_path)
    nodes = adapter.get_nodes()
    dm_nodes = [(nid, lbl, props) for nid, lbl, props in nodes if lbl == "derived_metric"]
    assert len(dm_nodes) == 1
    _, _, props = dm_nodes[0]
    assert props["metric_type"] == "fourier_score"
    assert props["value_kind"] == "numeric"
    assert props["rankable"] == "true"
    assert props["has_p_value"] == "true"
    assert props["p_value_threshold"] == 0.05
    assert props["unit"] == ""
    assert props["allowed_categories"] == []
```

- [ ] **Step 2: Run test to verify it passes**

Run: `uv run pytest tests/test_observations_adapter.py -v -k numeric_dm_node`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add tests/test_observations_adapter.py
git commit -m "observations_adapter: test numeric DerivedMetric node emission"
```

---

### Task 8: `get_edges` — binding edges (pub → DM, exp → DM, DM → org)

For each DM node, emit exactly 3 binding edges (Publication, Experiment, OrganismTaxon). Organism target is looked up via `self._organism_lookup[organism_name]`. No iteration over CSV rows in this task — binding edges are one-time-per-DM.

**Files:**
- Modify: `multiomics_kg/adapters/observations_adapter.py` — implement `get_edges` binding portion
- Modify: `tests/test_observations_adapter.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_observations_adapter.py`:

```python
def test_get_edges_emits_three_binding_edges_per_dm(tmp_path):
    pc_path = _write_s4a_like_paperconfig(tmp_path)
    adapter = ObservationsAdapter(config_file=pc_path)
    adapter._organism_lookup = {"Prochlorococcus NATL2A": "insdc.gcf:GCF_000012465.1"}
    edges = adapter.get_edges()

    pub_edges = [e for e in edges if e[3] == "publication_has_derived_metric"]
    exp_edges = [e for e in edges if e[3] == "experiment_has_derived_metric"]
    org_edges = [e for e in edges if e[3] == "derived_metric_belongs_to_organism"]

    # 2 DMs (axenic_LD + axenic_extended_darkness) × 3 binding types
    assert len(pub_edges) == 2
    assert len(exp_edges) == 2
    assert len(org_edges) == 2

    # Check shape: pub_edges are Pub → DM
    for eid, src, tgt, label, props in pub_edges:
        assert src == "doi:10.1128/mSystems.00040-18"
        assert tgt.startswith("derived_metric:mSystems.00040-18:s4a_axenic:")

    for eid, src, tgt, label, props in exp_edges:
        assert src == "10.1128/mSystems.00040-18_axenic_rnaseq"
        assert tgt.startswith("derived_metric:mSystems.00040-18:s4a_axenic:")

    for eid, src, tgt, label, props in org_edges:
        assert src.startswith("derived_metric:mSystems.00040-18:s4a_axenic:")
        assert tgt == "insdc.gcf:GCF_000012465.1"


def test_get_edges_skips_org_binding_when_lookup_misses(tmp_path):
    """If organism not in the lookup dict, no org edge (edge would dangle)."""
    pc_path = _write_s4a_like_paperconfig(tmp_path)
    adapter = ObservationsAdapter(config_file=pc_path)
    # lookup intentionally empty
    adapter._organism_lookup = {}
    edges = adapter.get_edges()
    org_edges = [e for e in edges if e[3] == "derived_metric_belongs_to_organism"]
    assert len(org_edges) == 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_observations_adapter.py -v -k get_edges_emits_three`
Expected: FAIL — `get_edges` still returns `[]`.

- [ ] **Step 3: Implement binding edges in `get_edges`**

Replace `get_edges` in `multiomics_kg/adapters/observations_adapter.py` with the full implementation (binding + all three measurement branches — stubbed in this task for brevity; measurement branches filled in Tasks 9–11):

```python
    def get_edges(self) -> list[tuple]:
        """Emit binding + measurement edges for each DerivedMetric."""
        edges = []
        experiments = get_experiments(self.config)
        pub_id = f"doi:{self.doi}" if self.doi else None

        for entry_key, entry in self._dm_entries:
            exp_key = entry.get("experiment")
            if not exp_key or exp_key not in experiments:
                continue
            exp = experiments[exp_key]
            organism = exp.get("organism", "")

            csv_path, use_resolved = _resolve_csv_path(entry["filename"])
            if not csv_path.exists():
                logger.warning(f"derived_metrics_table CSV not found: {csv_path}")
                continue
            df = pd.read_csv(csv_path)

            # Column used to reach the Gene node
            if use_resolved and "resolved_locus_tag" in df.columns:
                gene_col = "resolved_locus_tag"
                logger.info(
                    f"Using pre-resolved CSV: {csv_path.name} ({len(df)} rows)"
                )
            else:
                gene_col = entry.get("name_col", "")

            experiment_id = (
                f"{self.doi}_{exp_key}" if self.doi else f"{self.paper_name}_{exp_key}"
            )

            for metric in entry.get("metrics", []):
                metric_type = metric.get("metric_type", "")
                value_kind = metric.get("value_kind", "")
                if not metric_type or not value_kind:
                    continue

                dm_id = _make_derived_metric_id(
                    self.doi, self.paper_name, entry_key, metric_type
                )

                # --- Binding edges ---
                if pub_id:
                    edges.append((
                        f"pub_dm__{dm_id}",
                        pub_id, dm_id, "publication_has_derived_metric", {},
                    ))
                edges.append((
                    f"exp_dm__{dm_id}__{exp_key}",
                    experiment_id, dm_id, "experiment_has_derived_metric", {},
                ))
                if organism and organism in self._organism_lookup:
                    edges.append((
                        f"dm_org__{dm_id}",
                        dm_id, self._organism_lookup[organism],
                        "derived_metric_belongs_to_organism", {},
                    ))

                # --- Measurement edges (Tasks 9–11) ---
                value_col = metric.get("value_col", "")
                if not value_col or value_col not in df.columns:
                    logger.warning(
                        f"value_col '{value_col}' not in {csv_path.name} "
                        f"for '{entry_key}/{metric_type}' — skipping measurement edges"
                    )
                    continue

                if value_kind == "boolean":
                    edges.extend(self._emit_boolean_edges(
                        df, gene_col, value_col, metric, dm_id, metric_type,
                    ))
                elif value_kind == "categorical":
                    edges.extend(self._emit_categorical_edges(
                        df, gene_col, value_col, metric, dm_id, metric_type,
                    ))
                elif value_kind == "numeric":
                    edges.extend(self._emit_numeric_edges(
                        df, gene_col, value_col, metric, dm_id, metric_type,
                    ))

        return edges

    def _emit_boolean_edges(self, df, gene_col, value_col, metric, dm_id, metric_type):
        """Task 9 implements this."""
        return []

    def _emit_categorical_edges(self, df, gene_col, value_col, metric, dm_id, metric_type):
        """Task 10 implements this."""
        return []

    def _emit_numeric_edges(self, df, gene_col, value_col, metric, dm_id, metric_type):
        """Task 11 implements this."""
        return []
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_observations_adapter.py -v -k get_edges`
Expected: binding-edge tests PASS; measurement-edge tests (Tasks 9–11) not yet defined.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/observations_adapter.py tests/test_observations_adapter.py
git commit -m "observations_adapter: get_edges emits 3 binding edges per DerivedMetric"
```

---

### Task 9: `_emit_boolean_edges` — `derived_metric_flags_gene`

Implement the boolean measurement branch using `_parse_boolean_cell`. Rows with `resolved_locus_tag` NaN are skipped. Rows whose cell resolves to `None` (skip_token or blank-skip policy) are skipped. Unknown tokens raise `ValueError`.

**Files:**
- Modify: `multiomics_kg/adapters/observations_adapter.py` — fill `_emit_boolean_edges`
- Modify: `tests/test_observations_adapter.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_observations_adapter.py`:

```python
def test_boolean_edges_emit_true_flag_only_for_Y_rows(tmp_path):
    pc_path = _write_s4a_like_paperconfig(tmp_path)
    adapter = ObservationsAdapter(config_file=pc_path)
    adapter._organism_lookup = {"Prochlorococcus NATL2A": "insdc.gcf:GCF_000012465.1"}
    edges = adapter.get_edges()

    flag_edges = [e for e in edges if e[3] == "derived_metric_flags_gene"]
    # Fixture row expectations:
    #   periodic_in_axenic_LD: 2 Y rows (PMN2A_1328, PMN2A_1329), 1 blank row skipped,
    #                          1 unresolved row also blank (skipped) → 2 edges
    #   periodic_in_axenic_extended_darkness: 1 Y row (PMN2A_1329), 2 blank → 1 edge
    # Total = 3
    assert len(flag_edges) == 3

    for eid, src, tgt, label, props in flag_edges:
        assert src.startswith("derived_metric:mSystems.00040-18:s4a_axenic:")
        assert tgt.startswith("ncbigene:")
        assert props["value_flag"] == "true"
        assert props["metric_type"] in {
            "periodic_in_axenic_LD", "periodic_in_axenic_extended_darkness",
        }


def test_boolean_edges_skip_unresolved_genes(tmp_path):
    """Rows with resolved_locus_tag NaN must not produce edges."""
    pc_path = _write_s4a_like_paperconfig(tmp_path)
    adapter = ObservationsAdapter(config_file=pc_path)
    adapter._organism_lookup = {"Prochlorococcus NATL2A": "insdc.gcf:GCF_000012465.1"}
    edges = adapter.get_edges()
    flag_edges = [e for e in edges if e[3] == "derived_metric_flags_gene"]
    # Row 3 of the fixture has empty resolved_locus_tag — must not appear
    for eid, src, tgt, label, props in flag_edges:
        assert tgt != "ncbigene:"
        assert tgt != "ncbigene:nan"


def test_boolean_edges_raise_on_unknown_token(tmp_path):
    """Unexpected boolean tokens must hard-error."""
    csv_path = tmp_path / "bad.csv"
    csv_path.write_text("locus_tag,periodic\nPMM0001,maybe\n")
    # No resolved csv needed — name_col is locus_tag
    config = {
        "publication": {
            "papername": "Bad Tokens", "doi": "10.9999/bad",
            "papermainpdf": str(tmp_path / "fake.pdf"),
            "experiments": {"e": {
                "organism": "Prochlorococcus MED4", "omics_type": "RNASEQ",
                "treatment_type": [], "background_factors": [],
                "treatment_condition": "", "light_condition": "",
                "experimental_context": "",
            }},
            "supplementary_materials": {
                "entry": {
                    "type": "derived_metrics_table",
                    "filename": str(csv_path),
                    "organism": "Prochlorococcus MED4",
                    "experiment": "e",
                    "name_col": "locus_tag",
                    "metrics": [{
                        "metric_type": "periodic_in_axenic_LD",
                        "value_kind": "boolean",
                        "value_col": "periodic",
                        "true_tokens": ["Y"],
                        "false_tokens": [],
                        "skip_tokens": ["NA"],
                        "blank_policy": "skip",
                    }],
                },
            },
        },
    }
    pc_path = tmp_path / "paperconfig.yaml"
    pc_path.write_text(yaml.dump(config))
    adapter = ObservationsAdapter(config_file=pc_path)
    adapter._organism_lookup = {"Prochlorococcus MED4": "insdc.gcf:GCF_000011465.1"}
    with pytest.raises(ValueError, match="Unexpected boolean token"):
        adapter.get_edges()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_observations_adapter.py -v -k boolean_edges`
Expected: all 3 FAIL — `_emit_boolean_edges` returns `[]`.

- [ ] **Step 3: Implement `_emit_boolean_edges`**

Replace the stub in `multiomics_kg/adapters/observations_adapter.py` with:

```python
    def _emit_boolean_edges(self, df, gene_col, value_col, metric, dm_id, metric_type):
        """Emit derived_metric_flags_gene edges.

        Hard-errors on unexpected tokens (spec invariant: no silent coercion).
        Skips rows with unresolved locus_tag or skip_token cells.

        test_mode scope (intentional divergence from cluster_adapter): the
        100-edge cap is PER METRIC, not per entry. For a 2-metric entry in
        --test mode this caps at 200 edges total (vs cluster_adapter's
        shared-counter approach that caps at 100 per entry). Cleaner for
        DM observability — each DM gets a representative sample.
        """
        true_tokens = metric.get("true_tokens") or []
        false_tokens = metric.get("false_tokens") or []
        skip_tokens = metric.get("skip_tokens") or list(DEFAULT_SKIP_TOKENS)
        blank_policy = metric.get("blank_policy", "skip")

        edges = []
        count = 0
        for _, row in df.iterrows():
            if self.test_mode and count >= 100:
                break
            gene_locus = row.get(gene_col, "")
            if pd.isna(gene_locus):
                continue
            gene_locus = str(gene_locus).strip()
            if not gene_locus:
                continue

            flag = _parse_boolean_cell(
                row.get(value_col),
                true_tokens, false_tokens, skip_tokens, blank_policy,
            )
            if flag is None:
                continue  # skip / blank-skip

            edge_id = f"{dm_id}__{gene_locus}"
            props = {
                "metric_type": _clean_str(metric_type),
                "value_flag": flag,
            }
            edges.append((
                edge_id, dm_id, f"ncbigene:{gene_locus}",
                "derived_metric_flags_gene", props,
            ))
            count += 1
        return edges
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_observations_adapter.py -v -k boolean_edges`
Expected: 3/3 PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/observations_adapter.py tests/test_observations_adapter.py
git commit -m "observations_adapter: _emit_boolean_edges + hard-error on unknown tokens"
```

---

### Task 10: `_emit_categorical_edges` — `derived_metric_classifies_gene`

Categorical measurement: every row with a non-blank cell value whose cell ∈ `allowed_categories` gets one edge. Out-of-set values raise `ValueError`. NaN / blank rows are skipped.

**Files:**
- Modify: `multiomics_kg/adapters/observations_adapter.py`
- Modify: `tests/test_observations_adapter.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_observations_adapter.py`:

```python
def test_categorical_edges_emit_one_per_in_set_row(tmp_path):
    pc_path = _write_s5_like_paperconfig(tmp_path)
    adapter = ObservationsAdapter(config_file=pc_path)
    adapter._organism_lookup = {"Prochlorococcus NATL2A": "insdc.gcf:GCF_000012465.1"}
    edges = adapter.get_edges()
    cls_edges = [e for e in edges if e[3] == "derived_metric_classifies_gene"]
    assert len(cls_edges) == 3

    for eid, src, tgt, label, props in cls_edges:
        assert src.startswith("derived_metric:mSystems.00040-18:s5_survival:")
        assert tgt.startswith("ncbigene:")
        assert props["metric_type"] == "darkness_survival_class"
        assert props["value_text"] in {
            "darkness_axenic+darkness_coculture",
            "darkness_coculture+unique_coculture",
            "darkness_axenic+unique_axenic",
        }


def test_categorical_edges_raise_on_out_of_set(tmp_path):
    csv_path = tmp_path / "bad_cat.csv"
    csv_path.write_text("locus_tag,cluster\nPMM0001,unknown_category\n")
    config = {
        "publication": {
            "papername": "Bad Cat", "doi": "10.9999/badcat",
            "papermainpdf": str(tmp_path / "fake.pdf"),
            "experiments": {"e": {
                "organism": "Prochlorococcus MED4", "omics_type": "RNASEQ",
                "treatment_type": [], "background_factors": [],
                "treatment_condition": "", "light_condition": "",
                "experimental_context": "",
            }},
            "supplementary_materials": {
                "entry": {
                    "type": "derived_metrics_table",
                    "filename": str(csv_path),
                    "organism": "Prochlorococcus MED4",
                    "experiment": "e",
                    "name_col": "locus_tag",
                    "metrics": [{
                        "metric_type": "darkness_survival_class",
                        "value_kind": "categorical",
                        "value_col": "cluster",
                        "allowed_categories": ["a", "b"],
                    }],
                },
            },
        },
    }
    pc_path = tmp_path / "paperconfig.yaml"
    pc_path.write_text(yaml.dump(config))
    adapter = ObservationsAdapter(config_file=pc_path)
    adapter._organism_lookup = {"Prochlorococcus MED4": "insdc.gcf:GCF_000011465.1"}
    with pytest.raises(ValueError, match="out of allowed_categories"):
        adapter.get_edges()


def test_categorical_edges_skip_blank_rows(tmp_path):
    csv_path = tmp_path / "blank_cat.csv"
    csv_path.write_text("locus_tag,cluster\nPMM0001,a\nPMM0002,\nPMM0003,b\n")
    config = {
        "publication": {
            "papername": "Blank Cat", "doi": "10.9999/blankcat",
            "papermainpdf": str(tmp_path / "fake.pdf"),
            "experiments": {"e": {
                "organism": "Prochlorococcus MED4", "omics_type": "RNASEQ",
                "treatment_type": [], "background_factors": [],
                "treatment_condition": "", "light_condition": "",
                "experimental_context": "",
            }},
            "supplementary_materials": {
                "entry": {
                    "type": "derived_metrics_table",
                    "filename": str(csv_path),
                    "organism": "Prochlorococcus MED4",
                    "experiment": "e",
                    "name_col": "locus_tag",
                    "metrics": [{
                        "metric_type": "darkness_survival_class",
                        "value_kind": "categorical",
                        "value_col": "cluster",
                        "allowed_categories": ["a", "b"],
                    }],
                },
            },
        },
    }
    pc_path = tmp_path / "paperconfig.yaml"
    pc_path.write_text(yaml.dump(config))
    adapter = ObservationsAdapter(config_file=pc_path)
    adapter._organism_lookup = {"Prochlorococcus MED4": "insdc.gcf:GCF_000011465.1"}
    edges = adapter.get_edges()
    cls_edges = [e for e in edges if e[3] == "derived_metric_classifies_gene"]
    assert len(cls_edges) == 2  # PMM0002's blank row skipped
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_observations_adapter.py -v -k categorical_edges`
Expected: all 3 FAIL.

- [ ] **Step 3: Implement `_emit_categorical_edges`**

Replace the stub:

```python
    def _emit_categorical_edges(self, df, gene_col, value_col, metric, dm_id, metric_type):
        """Emit derived_metric_classifies_gene edges.

        Hard-errors on values not in allowed_categories (spec invariant #11).
        """
        allowed = set(metric.get("allowed_categories") or [])

        edges = []
        count = 0
        for _, row in df.iterrows():
            if self.test_mode and count >= 100:
                break
            gene_locus = row.get(gene_col, "")
            if pd.isna(gene_locus):
                continue
            gene_locus = str(gene_locus).strip()
            if not gene_locus:
                continue
            val = row.get(value_col)
            if pd.isna(val):
                continue
            s = str(val).strip()
            if s == "":
                continue
            if s not in allowed:
                raise ValueError(
                    f"Categorical value {s!r} out of allowed_categories "
                    f"{sorted(allowed)} for metric {metric_type!r}"
                )
            edge_id = f"{dm_id}__{gene_locus}"
            edges.append((
                edge_id, dm_id, f"ncbigene:{gene_locus}",
                "derived_metric_classifies_gene",
                {"metric_type": _clean_str(metric_type), "value_text": _clean_str(s)},
            ))
            count += 1
        return edges
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_observations_adapter.py -v -k categorical_edges`
Expected: 3/3 PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/observations_adapter.py tests/test_observations_adapter.py
git commit -m "observations_adapter: _emit_categorical_edges + allowed_categories check"
```

---

### Task 11: `_emit_numeric_edges` — `derived_metric_quantifies_gene`

Numeric measurement: `value` is required and must be numeric; `p_value` + `adjusted_p_value` are optional. Rows with NaN `value` are skipped. `p_value` / `adjusted_p_value` columns may or may not be in the paperconfig — when absent, don't set the key on the edge.

**Files:**
- Modify: `multiomics_kg/adapters/observations_adapter.py`
- Modify: `tests/test_observations_adapter.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_observations_adapter.py`:

```python
def test_numeric_edges_emit_one_per_row_with_value(tmp_path):
    pc_path = _write_numeric_dm_paperconfig(tmp_path)
    adapter = ObservationsAdapter(config_file=pc_path)
    adapter._organism_lookup = {"Prochlorococcus MED4": "insdc.gcf:GCF_000011465.1"}
    edges = adapter.get_edges()
    q_edges = [e for e in edges if e[3] == "derived_metric_quantifies_gene"]
    assert len(q_edges) == 2

    for eid, src, tgt, label, props in q_edges:
        assert src.startswith("derived_metric:synthetic:fourier_metrics:")
        assert tgt.startswith("ncbigene:PMM")
        assert props["metric_type"] == "fourier_score"
        assert "value" in props
        assert isinstance(props["value"], float)
        # p_value column declared → key present
        assert "p_value" in props
        assert isinstance(props["p_value"], float)
        # adjusted_p_value not declared in this fixture → key absent
        assert "adjusted_p_value" not in props


def test_numeric_edges_skip_rows_with_nan_value(tmp_path):
    csv_path = tmp_path / "nan_val.csv"
    csv_path.write_text("locus_tag,v\nPMM0001,0.5\nPMM0002,\nPMM0003,0.8\n")
    config = {
        "publication": {
            "papername": "NaN Val", "doi": "10.9999/nan",
            "papermainpdf": str(tmp_path / "fake.pdf"),
            "experiments": {"e": {
                "organism": "Prochlorococcus MED4", "omics_type": "RNASEQ",
                "treatment_type": [], "background_factors": [],
                "treatment_condition": "", "light_condition": "",
                "experimental_context": "",
            }},
            "supplementary_materials": {
                "entry": {
                    "type": "derived_metrics_table",
                    "filename": str(csv_path),
                    "organism": "Prochlorococcus MED4",
                    "experiment": "e",
                    "name_col": "locus_tag",
                    "metrics": [{
                        "metric_type": "fourier_score",
                        "value_kind": "numeric",
                        "value_col": "v",
                        "rankable": "true",
                        "has_p_value": "false",
                    }],
                },
            },
        },
    }
    pc_path = tmp_path / "paperconfig.yaml"
    pc_path.write_text(yaml.dump(config))
    adapter = ObservationsAdapter(config_file=pc_path)
    adapter._organism_lookup = {"Prochlorococcus MED4": "insdc.gcf:GCF_000011465.1"}
    edges = adapter.get_edges()
    q_edges = [e for e in edges if e[3] == "derived_metric_quantifies_gene"]
    assert len(q_edges) == 2  # PMM0002 row with blank `v` skipped
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_observations_adapter.py -v -k numeric_edges`
Expected: FAIL — `_emit_numeric_edges` stub returns `[]`.

- [ ] **Step 3: Implement `_emit_numeric_edges`**

Replace the stub:

```python
    def _emit_numeric_edges(self, df, gene_col, value_col, metric, dm_id, metric_type):
        """Emit derived_metric_quantifies_gene edges.

        Skips rows where `value` (the `value_col` cell) is NaN.
        Carries optional `p_value` / `adjusted_p_value` only when the paperconfig
        declares the corresponding column AND the cell is non-null.
        """
        p_value_col = metric.get("p_value_col")
        adj_p_value_col = metric.get("adjusted_p_value_col")

        edges = []
        count = 0
        for _, row in df.iterrows():
            if self.test_mode and count >= 100:
                break
            gene_locus = row.get(gene_col, "")
            if pd.isna(gene_locus):
                continue
            gene_locus = str(gene_locus).strip()
            if not gene_locus:
                continue
            raw_val = row.get(value_col)
            if pd.isna(raw_val):
                continue
            try:
                value = float(raw_val)
            except (TypeError, ValueError):
                continue  # non-numeric cell → skip (not a hard error; paper may flag as blank)

            props = {
                "metric_type": _clean_str(metric_type),
                "value": value,
            }
            if p_value_col and p_value_col in df.columns:
                pv = row.get(p_value_col)
                if not pd.isna(pv):
                    try:
                        props["p_value"] = float(pv)
                    except (TypeError, ValueError):
                        pass
            if adj_p_value_col and adj_p_value_col in df.columns:
                apv = row.get(adj_p_value_col)
                if not pd.isna(apv):
                    try:
                        props["adjusted_p_value"] = float(apv)
                    except (TypeError, ValueError):
                        pass

            edge_id = f"{dm_id}__{gene_locus}"
            edges.append((
                edge_id, dm_id, f"ncbigene:{gene_locus}",
                "derived_metric_quantifies_gene", props,
            ))
            count += 1
        return edges
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_observations_adapter.py -v -k numeric_edges`
Expected: 2/2 PASS.

- [ ] **Step 5: Full module regression**

Run: `uv run pytest tests/test_observations_adapter.py -v`
Expected: all tests PASS (all prior tests still green).

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/adapters/observations_adapter.py tests/test_observations_adapter.py
git commit -m "observations_adapter: _emit_numeric_edges with optional p-value plumbing"
```

---

### Task 12: `MultiObservationsAdapter` — organism lookup + paperconfig filtering

`MultiObservationsAdapter` was scaffolded in Task 2. Test it now that the real `get_nodes` / `get_edges` are wired.

**Files:**
- Modify: `tests/test_observations_adapter.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_observations_adapter.py`:

```python
def test_multi_adapter_skips_paperconfigs_without_derived_metrics(tmp_path):
    """Paperconfigs with no derived_metrics_table entries must be skipped."""
    # paperconfig_a: has a derived_metrics_table entry
    pc_a = _write_s4a_like_paperconfig(tmp_path / "a")
    # paperconfig_b: a csv-type entry, no derived_metrics_table
    (tmp_path / "b").mkdir()
    pc_b_path = tmp_path / "b" / "paperconfig.yaml"
    pc_b_path.write_text(yaml.dump({
        "publication": {
            "papername": "Other 2020", "doi": "10.9999/other",
            "papermainpdf": str(tmp_path / "b" / "fake.pdf"),
            "experiments": {"e": {
                "organism": "Prochlorococcus MED4", "omics_type": "RNASEQ",
                "treatment_type": [], "background_factors": [],
                "treatment_condition": "", "light_condition": "",
                "experimental_context": "",
            }},
            "supplementary_materials": {
                "non_dm": {"type": "csv", "filename": "x.csv"},
            },
        },
    }))

    list_path = tmp_path / "paperconfig_files.txt"
    list_path.write_text(f"{pc_a}\n{pc_b_path}\n")

    multi = MultiObservationsAdapter(config_list_file=str(list_path))
    assert len(multi.adapters) == 1
    assert multi.adapters[0].config_file == pc_a


def test_multi_adapter_builds_organism_lookup(tmp_path):
    """_organism_lookup keys match the genome CSV preferred_name column."""
    genome_csv = tmp_path / "genomes.csv"
    genome_csv.write_text(
        "ncbi_accession,,taxid,strain,data_dir,clade,preferred_name,organism_type,,\n"
        "GCF_000012465.1,,59920,NATL2A,dir/,LLII,Prochlorococcus NATL2A,genome_strain,,\n"
    )
    multi = MultiObservationsAdapter(
        config_list_file=str(tmp_path / "empty.txt"),
        genome_config_file=str(genome_csv),
    )
    # Even with no paperconfigs, the lookup should be built
    assert multi._organism_lookup.get("Prochlorococcus NATL2A") == "insdc.gcf:GCF_000012465.1"


def test_multi_adapter_propagates_organism_lookup_to_children(tmp_path):
    pc_a = _write_s4a_like_paperconfig(tmp_path / "a")
    genome_csv = tmp_path / "genomes.csv"
    genome_csv.write_text(
        "ncbi_accession,,taxid,strain,data_dir,clade,preferred_name,organism_type,,\n"
        "GCF_000012465.1,,59920,NATL2A,dir/,LLII,Prochlorococcus NATL2A,genome_strain,,\n"
    )
    list_path = tmp_path / "paperconfig_files.txt"
    list_path.write_text(f"{pc_a}\n")

    multi = MultiObservationsAdapter(
        config_list_file=str(list_path),
        genome_config_file=str(genome_csv),
    )
    assert len(multi.adapters) == 1
    assert multi.adapters[0]._organism_lookup["Prochlorococcus NATL2A"] == "insdc.gcf:GCF_000012465.1"

    edges = multi.get_edges()
    org_edges = [e for e in edges if e[3] == "derived_metric_belongs_to_organism"]
    assert len(org_edges) == 2  # one per DM
    for _, _, tgt, _, _ in org_edges:
        assert tgt == "insdc.gcf:GCF_000012465.1"
```

- [ ] **Step 2: Run tests**

Run: `uv run pytest tests/test_observations_adapter.py -v -k multi_adapter`
Expected: 3/3 PASS (MultiObservationsAdapter was wired fully in Task 2).

> **Note:** The `empty.txt` list file in `test_multi_adapter_builds_organism_lookup` doesn't exist; `load_all_paperconfigs` currently opens it via `open(lf)`. If this raises `FileNotFoundError`, touch an empty file in the test before instantiating:
>
> ```python
> (tmp_path / "empty.txt").write_text("")
> ```
> Add this line to the test if needed during Step 2.

- [ ] **Step 3: Commit**

```bash
git add tests/test_observations_adapter.py
git commit -m "observations_adapter: MultiObservationsAdapter organism lookup tests"
```

---

### Task 13: Synthetic numeric-DM paperconfig fixture (100 rows)

The Biller 2018 paper has no numeric metrics, so the numeric branch isn't exercised end-to-end by production data until a future slice (zinser 2009 / Waldbauer 2012). Plan 2 ships a synthetic 100-row fixture to ensure the numeric path survives refactors. Put the fixture under `tests/fixtures/non_de/` and verify via a unit test that `MultiObservationsAdapter` ingests it cleanly. The fixture is **test-only** — it is NOT appended to the production paperconfig_files.txt.

Design decisions (per the user's Plan 2 brief):
- Organism: `Prochlorococcus MED4` (already in cyanobacteria_genomes.csv → organism lookup works)
- Locus tags: 100 real MED4 tags sequential from `PMM0001` through `PMM0100` — all resolve 1:1 to Gene nodes
- Metrics: 3 × numeric (`fourier_score`, `peak_time_h`, `peak_fit_r_squared`) — covers `rankable` + `has_p_value` + `p_value_threshold` variations
- p-values distribution: ~20 rows below 0.05, ~30 rows 0.05–0.2, ~50 rows 0.5+ — exercises bucket splits in Plan 3
- Value distribution: sorted numeric — ensures rank_by_metric / metric_percentile / metric_bucket have deterministic expected values

**Files:**
- Create: `tests/fixtures/non_de/__init__.py` (empty)
- Create: `tests/fixtures/non_de/synthetic_numeric_metrics.csv`
- Create: `tests/fixtures/non_de/synthetic_paperconfig.yaml`
- Create: `tests/fixtures/non_de/paperconfig_files.txt`
- Modify: `tests/test_observations_adapter.py` (add fixture-ingest test)

- [ ] **Step 1: Create the fixture CSV**

Create `tests/fixtures/non_de/synthetic_numeric_metrics.csv` with this Python generator (run once, commit the output):

```bash
uv run python -c "
import csv, math, random
random.seed(42)
rows = []
for i in range(1, 101):
    lt = f'PMM{i:04d}'
    # fourier: descending from 0.98 to 0.01
    fourier = 0.98 - (i - 1) * (0.97 / 99)
    # peak_time: spread across 0-24 hours
    peak_time = round((i * 0.25) % 24, 2)
    # peak_fit: descending from 1.0 to 0.2, noisy
    peak_fit = max(0.2, 1.0 - (i - 1) * (0.8 / 99) + random.uniform(-0.05, 0.05))
    # p-values: first 20 significant, next 30 marginal, last 50 weak
    if i <= 20:
        p_val = 0.001 * (i / 20)
        adj_p = 0.005 * (i / 20)
    elif i <= 50:
        p_val = 0.05 + 0.15 * ((i - 20) / 30)
        adj_p = 0.1 + 0.2 * ((i - 20) / 30)
    else:
        p_val = 0.4 + 0.5 * ((i - 50) / 50)
        adj_p = 0.6 + 0.35 * ((i - 50) / 50)
    rows.append([lt, f'{fourier:.4f}', f'{peak_time:.2f}', f'{peak_fit:.4f}', f'{p_val:.4f}', f'{adj_p:.4f}'])

with open('tests/fixtures/non_de/synthetic_numeric_metrics.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['locus_tag', 'fourier', 'peak_time', 'peak_fit', 'p_value', 'adj_p_value'])
    w.writerows(rows)
print('wrote 100 rows')
"
```

- [ ] **Step 2: Create the synthetic paperconfig YAML**

Create `tests/fixtures/non_de/synthetic_paperconfig.yaml`:

```yaml
publication:
  papername: "Synthetic Numeric DM Fixture"
  doi: "10.9999/synthetic-numeric-dm"
  papermainpdf: "tests/fixtures/non_de/NONEXISTENT.pdf"  # placeholder; no PDF extractor run
  experiments:
    synthetic_diel_rnaseq:
      name: "Synthetic MED4 diel RNA-seq"
      organism: "Prochlorococcus MED4"
      omics_type: RNASEQ
      treatment_type: [diel]
      background_factors: [axenic]
      treatment_condition: "Diel cycle"
      control_condition: ""
      experimental_context: "Synthetic MED4 in Pro99 at 24C"
      light_condition: "14:10 LD"
  supplementary_materials:
    synthetic_fourier_metrics:
      type: derived_metrics_table
      filename: "tests/fixtures/non_de/synthetic_numeric_metrics.csv"
      organism: "Prochlorococcus MED4"
      experiment: synthetic_diel_rnaseq
      name_col: locus_tag
      id_columns:
        - column: locus_tag
          id_type: locus_tag
      metrics:
        - metric_type: fourier_score
          value_kind: numeric
          value_col: fourier
          unit: ""
          rankable: "true"
          has_p_value: "true"
          p_value_col: p_value
          adjusted_p_value_col: adj_p_value
          p_value_threshold: 0.05
          field_description: "Fourier periodicity score (synthetic)"
        - metric_type: peak_time_h
          value_kind: numeric
          value_col: peak_time
          unit: "h"
          rankable: "false"
          has_p_value: "false"
          field_description: "Peak time within 24h diel cycle (synthetic)"
        - metric_type: peak_fit_r_squared
          value_kind: numeric
          value_col: peak_fit
          unit: ""
          rankable: "true"
          has_p_value: "false"
          field_description: "Peak-fit R^2 (synthetic)"
```

- [ ] **Step 3: Create the list file**

Create `tests/fixtures/non_de/paperconfig_files.txt`:

```
tests/fixtures/non_de/synthetic_paperconfig.yaml
```

- [ ] **Step 4: Create the empty `__init__.py`**

Create `tests/fixtures/non_de/__init__.py` (empty file — ensures pytest discovers the directory cleanly).

- [ ] **Step 5: Write the failing integration test**

Append to `tests/test_observations_adapter.py`:

```python
def test_synthetic_numeric_fixture_emits_three_dms_and_300_edges():
    """Regression guard for the numeric code path.

    Fixture: tests/fixtures/non_de/synthetic_paperconfig.yaml references a
    100-row CSV with 3 numeric metrics. Expected:
      - 3 DerivedMetric nodes
      - 3 binding edges × 3 DMs = 9 binding edges (pub/exp/org × 3)
      - 100 × 3 = 300 derived_metric_quantifies_gene edges
    """
    list_file = "tests/fixtures/non_de/paperconfig_files.txt"
    # Genome lookup for organism edge target
    genome_csv = "data/Prochlorococcus/genomes/cyanobacteria_genomes.csv"
    multi = MultiObservationsAdapter(
        config_list_file=list_file,
        genome_config_file=genome_csv,
    )
    assert len(multi.adapters) == 1

    nodes = multi.get_nodes()
    dm_nodes = [n for n in nodes if n[1] == "derived_metric"]
    assert len(dm_nodes) == 3

    edges = multi.get_edges()
    by_type = {}
    for e in edges:
        by_type.setdefault(e[3], 0)
        by_type[e[3]] += 1

    assert by_type.get("publication_has_derived_metric", 0) == 3
    assert by_type.get("experiment_has_derived_metric", 0) == 3
    assert by_type.get("derived_metric_belongs_to_organism", 0) == 3
    assert by_type.get("derived_metric_quantifies_gene", 0) == 300
    # No boolean or categorical edges
    assert "derived_metric_flags_gene" not in by_type
    assert "derived_metric_classifies_gene" not in by_type
```

- [ ] **Step 6: Run the test**

Run: `uv run pytest tests/test_observations_adapter.py -v -k synthetic_numeric_fixture`
Expected: PASS.

> If it fails with a PDF extraction error, the fixture paperconfig's `papermainpdf` must point at a non-existent path (as written above). The ObservationsAdapter only reads the PDF-extraction cache and never invokes the extractor itself — so missing PDF is fine. (The `OMICSAdapter`, which DOES invoke the extractor, is not exercised in this test.)

- [ ] **Step 7: Run the validator against the synthetic paperconfig**

Run: `uv run python scripts/validate_paperconfig.py tests/fixtures/non_de/synthetic_paperconfig.yaml`
Expected: exit 0 (possibly a `novel metric_type` warning if Plan 1 registered only Biller-2018 metrics — that's acceptable; all three fixture metric_types are already in `KNOWN_METRIC_TYPES`).

- [ ] **Step 8: Commit**

```bash
git add tests/fixtures/non_de/ tests/test_observations_adapter.py
git commit -m "observations_adapter: synthetic 100-row numeric-DM fixture + ingest test"
```

---

### Task 14: Wire `MultiObservationsAdapter` into `create_knowledge_graph.py`

After `MultiClusterAdapter`, instantiate `MultiObservationsAdapter` with the same `config_list_file` + `genome_config_file` args and call `bc.write_nodes()` / `bc.write_edges()`.

**Files:**
- Modify: `create_knowledge_graph.py`

- [ ] **Step 1: Add the import at the top of the file**

In `create_knowledge_graph.py`, add this line alongside the existing adapter imports (e.g., right after `from multiomics_kg.adapters.cluster_adapter import MultiClusterAdapter`):

```python
from multiomics_kg.adapters.observations_adapter import MultiObservationsAdapter
```

- [ ] **Step 2: Wire the adapter into `main()`**

Find the existing `MultiClusterAdapter` block (around lines 92–106). Immediately AFTER the existing `if cluster_edges: bc.write_edges(cluster_edges)` block, add:

```python
    # Derived-metric observations (non-DE evidence; Biller 2018 boolean + categorical,
    # future zinser 2009 / Waldbauer 2012 numeric). One DerivedMetric node per
    # (Experiment × metric_type) + 3 binding edges + 1 measurement edge per metric.
    observations_adapter = MultiObservationsAdapter(
        config_list_file=[
            'data/Prochlorococcus/papers_and_supp/paperconfig_files.txt',
            'data/Synechococcus/papers_and_supp/paperconfig_files.txt',
        ],
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        test_mode=TEST_MODE,
    )
    obs_nodes = observations_adapter.get_nodes()
    obs_edges = observations_adapter.get_edges()
    if obs_nodes:
        bc.write_nodes(obs_nodes)
    if obs_edges:
        bc.write_edges(obs_edges)
```

- [ ] **Step 3: Smoke test — unit suite still passes**

Run: `uv run pytest -m "not slow and not kg" -x`
Expected: all tests pass. The new adapter is additive; no existing tests should regress.

- [ ] **Step 4: Smoke test — `create_knowledge_graph.py --test` runs to completion**

Run: `uv run python create_knowledge_graph.py --test`
Expected: the script runs to the end. BioCypher log output should include lines about writing DerivedMetric nodes and new edge types. The final line prints the import-call script path.

- [ ] **Step 5: Verify new CSV outputs exist**

Run: `ls biocypher-log/example_knowledge_graph/ | grep -i derivedmetric`
Expected: at least one `DerivedMetric*.csv` (nodes) and one `DerivedMetricQuantifiesGene*.csv` / `DerivedMetricFlagsGene*.csv` / `DerivedMetricClassifiesGene*.csv` (edges, depending on which entries were picked up in test mode).

> **Debug note:** in `--test` mode each adapter is limited to 100 items. For ObservationsAdapter this means at most 100 measurement edges per metric. That's plenty to verify CSV emission works but too few to validate Biller 2018's full edge counts — which the next task does against a full run.

- [ ] **Step 6: Commit**

```bash
git add create_knowledge_graph.py
git commit -m "create_knowledge_graph: wire MultiObservationsAdapter after cluster adapter"
```

---

### Task 15: End-to-end Biller 2018 verification

The real paperconfig (authored in Plan 1) ships with 4 `derived_metrics_table` entries declaring 7 DerivedMetric specs total:
- `s4a_natl2a_axenic`: 2 boolean metrics (`periodic_in_axenic_LD`, `periodic_in_axenic_extended_darkness`) → DMs bound to `darkness_extended_darkness_natl2a_rnaseq_axenic`
- `s4a_natl2a_coculture`: 2 boolean metrics (`periodic_in_coculture_LD`, `periodic_in_coculture_extended_darkness`) → DMs bound to `darkness_extended_darkness_natl2a_rnaseq_coculture`
- `s4b_mit1002`: 2 boolean metrics (`periodic_in_coculture_LD`, `periodic_in_coculture_extended_darkness`) → DMs bound to `darkness_extended_darkness_mit1002_rnaseq`
- `s5_natl2a_survival`: 1 categorical metric (`darkness_survival_class`) → DM bound to `darkness_extended_darkness_natl2a_rnaseq_axenic`

Expected counts (verified against the _resolved.csv files under `data/Prochlorococcus/papers_and_supp/Biller 2018/`):

| Edge type | Expected count (approx) |
|---|---|
| `derived_metric_flags_gene` | ≈ 4,100 (S4A axenic: 1515 + 127; S4A coculture: 1877 + 555; S4B: 530 + 2 — minus ~14% unresolved filter on NATL2A rows) |
| `derived_metric_classifies_gene` | ≈ 258 (S5 resolved rows; 100 + 90 + 79 − 11 unresolved) |
| `derived_metric_quantifies_gene` | 0 (Biller 2018 has no numeric metrics) |
| `publication_has_derived_metric` | 7 |
| `experiment_has_derived_metric` | 7 |
| `derived_metric_belongs_to_organism` | 7 |
| `DerivedMetric` nodes | 7 |

**Files:**
- Create: `tests/test_observations_adapter_biller2018.py` — full-paperconfig integration test

- [ ] **Step 1: Write the test**

Create `tests/test_observations_adapter_biller2018.py`:

```python
"""End-to-end Biller 2018 integration test for ObservationsAdapter.

Runs against the real paperconfig at
`data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml`.
Asserts expected DerivedMetric + binding-edge + measurement-edge counts
match the Plan 2 budget (see docs/superpowers/plans/2026-04-19-plan2-...).

Requires Plan 1's _resolved.csv files to be present (run prepare_data.sh
step 4 if they're missing).
"""
from pathlib import Path

import pytest

from multiomics_kg.adapters.observations_adapter import ObservationsAdapter


BILLER_PC = Path("data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml")
BILLER_DOI = "10.1128/mSystems.00040-18"


@pytest.mark.skipif(
    not BILLER_PC.exists(),
    reason="Biller 2018 paperconfig not present (Plan 1 author step)",
)
def test_biller_2018_dm_nodes():
    adapter = ObservationsAdapter(config_file=str(BILLER_PC))
    adapter._organism_lookup = {
        "Prochlorococcus NATL2A": "insdc.gcf:GCF_000012465.1",
        "Alteromonas macleodii MIT1002": "insdc.gcf:GCF_901457835.2",
    }
    nodes = adapter.get_nodes()
    dm_nodes = [n for n in nodes if n[1] == "derived_metric"]
    # 4 entries × 2/2/2/1 metrics = 7 DerivedMetric nodes
    assert len(dm_nodes) == 7

    metric_types = sorted({props["metric_type"] for _, _, props in dm_nodes})
    assert metric_types == [
        "darkness_survival_class",
        "periodic_in_axenic_LD",
        "periodic_in_axenic_extended_darkness",
        "periodic_in_coculture_LD",
        "periodic_in_coculture_extended_darkness",
        # the two MIT1002 entries reuse periodic_in_coculture_LD /
        # periodic_in_coculture_extended_darkness metric_types but with
        # different DM node IDs (different entry_key in the ID)
    ]
    # Sanity: 4 distinct DMs carry value_kind=boolean from NATL2A axenic + coculture;
    # 2 boolean from MIT1002; 1 categorical from S5
    by_kind = {}
    for _, _, props in dm_nodes:
        by_kind.setdefault(props["value_kind"], 0)
        by_kind[props["value_kind"]] += 1
    assert by_kind == {"boolean": 6, "categorical": 1}


@pytest.mark.skipif(
    not BILLER_PC.exists(),
    reason="Biller 2018 paperconfig not present",
)
def test_biller_2018_binding_edges():
    adapter = ObservationsAdapter(config_file=str(BILLER_PC))
    adapter._organism_lookup = {
        "Prochlorococcus NATL2A": "insdc.gcf:GCF_000012465.1",
        "Alteromonas macleodii MIT1002": "insdc.gcf:GCF_901457835.2",
    }
    edges = adapter.get_edges()
    by_type = {}
    for e in edges:
        by_type.setdefault(e[3], 0)
        by_type[e[3]] += 1

    # 7 DMs → 7 pub + 7 exp + 7 org binding edges
    assert by_type.get("publication_has_derived_metric", 0) == 7
    assert by_type.get("experiment_has_derived_metric", 0) == 7
    assert by_type.get("derived_metric_belongs_to_organism", 0) == 7


@pytest.mark.skipif(
    not BILLER_PC.exists(),
    reason="Biller 2018 paperconfig not present",
)
def test_biller_2018_flag_edge_count_in_band():
    """Boolean measurement edges — expected band 3,500 - 4,700 (accounts for
    resolved-row filter; exact count depends on prepare_data.sh step 4 version).
    """
    adapter = ObservationsAdapter(config_file=str(BILLER_PC))
    adapter._organism_lookup = {
        "Prochlorococcus NATL2A": "insdc.gcf:GCF_000012465.1",
        "Alteromonas macleodii MIT1002": "insdc.gcf:GCF_901457835.2",
    }
    edges = adapter.get_edges()
    flag_edges = [e for e in edges if e[3] == "derived_metric_flags_gene"]
    assert 3500 <= len(flag_edges) <= 4700, (
        f"derived_metric_flags_gene count {len(flag_edges)} outside expected band"
    )


@pytest.mark.skipif(
    not BILLER_PC.exists(),
    reason="Biller 2018 paperconfig not present",
)
def test_biller_2018_classifies_edge_count():
    """Categorical measurement edges — S5 has 258 resolved rows (verified)."""
    adapter = ObservationsAdapter(config_file=str(BILLER_PC))
    adapter._organism_lookup = {
        "Prochlorococcus NATL2A": "insdc.gcf:GCF_000012465.1",
    }
    edges = adapter.get_edges()
    cls_edges = [e for e in edges if e[3] == "derived_metric_classifies_gene"]
    # Tolerance: 255–270 rows depending on resolver edge cases
    assert 250 <= len(cls_edges) <= 270, (
        f"derived_metric_classifies_gene count {len(cls_edges)} outside expected band"
    )


@pytest.mark.skipif(
    not BILLER_PC.exists(),
    reason="Biller 2018 paperconfig not present",
)
def test_biller_2018_no_numeric_edges():
    """Biller 2018 has no numeric metrics — DM quantifies_gene count must be 0."""
    adapter = ObservationsAdapter(config_file=str(BILLER_PC))
    adapter._organism_lookup = {
        "Prochlorococcus NATL2A": "insdc.gcf:GCF_000012465.1",
        "Alteromonas macleodii MIT1002": "insdc.gcf:GCF_901457835.2",
    }
    edges = adapter.get_edges()
    q_edges = [e for e in edges if e[3] == "derived_metric_quantifies_gene"]
    assert len(q_edges) == 0
```

- [ ] **Step 2: Run the test**

Run: `uv run pytest tests/test_observations_adapter_biller2018.py -v`
Expected: all 5 PASS.

> If `flag_edge_count` or `classifies_edge_count` falls outside the asserted band, either the Plan 1 `_resolved.csv` files are stale (run `uv run python -m multiomics_kg.download.resolve_paper_ids --papers "Biller 2018" --force`) OR the resolved-row counts differ from the numbers captured in the plan. The plan's bands have been deliberately widened to absorb reasonable resolver differences; if the count is wildly off, debug before adjusting the band.

- [ ] **Step 3: Commit**

```bash
git add tests/test_observations_adapter_biller2018.py
git commit -m "observations_adapter: end-to-end Biller 2018 count assertions"
```

---

### Task 16: Full unit-test suite + Docker import smoke + Changes_expression_of parity

Final gate before handing off to Plan 3.

- [ ] **Step 1: Full unit suite**

Run: `uv run pytest -m "not slow and not kg" -v`
Expected: all pass. Take particular note of:
- `tests/test_observations_adapter.py` — all ObservationsAdapter unit tests
- `tests/test_observations_adapter_biller2018.py` — integration assertions
- `tests/test_omics_adapter_experiment.py` — compartment emission (Task 1)
- `tests/test_paperconfig_utils.py` — ensure `iter_derived_metrics_tables` still green (Plan 1 tests)

- [ ] **Step 2: Rebuild the Docker stack end-to-end**

```bash
docker compose down
docker compose up -d --build
```

Wait for all stages (`build`, `import`, `post-process`, `deploy`, `app`) to reach "exited" or "running" state. Watch the logs:

```bash
docker compose logs -f import
```

Expected:
- `import` exits 0
- `output/import.status` contains a `0` exit code line (via `scripts/import.sh`)
- `output/import.report` shows zero skipped relationships in the `non-de` categories:
  ```
  grep -iE "derived_metric|DerivedMetric" output/import.report | grep -iE "skipped|skipping" || echo "no skipped derived-metric relationships"
  ```
  Expected: `no skipped derived-metric relationships`

- [ ] **Step 3: Verify `Changes_expression_of` count is unchanged**

```bash
docker compose exec deploy cypher-shell -u neo4j -p neo4j \
  "MATCH ()-[r:Changes_expression_of]->() RETURN count(r) AS n"
```

Expected: `227361` (matches step-0 baseline in `docs/kg-changes/biller-2018-retrofit-baseline/delta.md`).

- [ ] **Step 4: Verify new DerivedMetric node + edge counts against expected bands**

```bash
docker compose exec deploy cypher-shell -u neo4j -p neo4j "
MATCH (d:DerivedMetric) RETURN count(d) AS dm_nodes;
MATCH ()-[r:Derived_metric_flags_gene]->() RETURN count(r) AS flag_edges;
MATCH ()-[r:Derived_metric_classifies_gene]->() RETURN count(r) AS cls_edges;
MATCH ()-[r:Derived_metric_quantifies_gene]->() RETURN count(r) AS q_edges;
MATCH ()-[r:Publication_has_derived_metric]->() RETURN count(r) AS pub_dm_edges;
MATCH ()-[r:Experiment_has_derived_metric]->() RETURN count(r) AS exp_dm_edges;
MATCH ()-[r:Derived_metric_belongs_to_organism]->() RETURN count(r) AS dm_org_edges;
"
```

Expected:
- `dm_nodes`: 7 (Biller 2018 only — synthetic fixture is test-only, not wired into production)
- `flag_edges`: in band [3500, 4700]
- `cls_edges`: in band [250, 270]
- `q_edges`: 0
- `pub_dm_edges`: 7
- `exp_dm_edges`: 7
- `dm_org_edges`: 7

- [ ] **Step 5: Run `/omics-edge-snapshot` to confirm no DE regression**

The skill counts `Changes_expression_of` edges across papers against the live graph. DerivedMetric edge counts aren't surfaced yet — that's Plan 3. Save a baseline and diff against step-0:

```bash
# Capture after-Plan-2 snapshot from live graph
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --save after_plan2

# Compare step-0 baseline to after_plan2
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py \
    --compare after_biller2018_retrofit_removal --against after_plan2
```

Expected: `Total expression edges: 227,361 → 227,361  (+0)` and `No regressions — no publication lost edges`.

- [ ] **Step 6: Commit the snapshot**

```bash
git add .claude/skills/omics-edge-snapshot/snapshots/after_plan2.json
git commit -m "omics-edge-snapshot: after_plan2 — no DE regression vs step-0 baseline"
```

- [ ] **Step 7: Final verification summary**

If all prior steps passed:
- ✅ `pytest -m "not slow and not kg"` green
- ✅ `docker compose up -d --build` exits clean, zero skipped relationships
- ✅ `Changes_expression_of` count unchanged vs step-0 baseline (227,361)
- ✅ 7 DerivedMetric nodes + expected binding/measurement edge counts
- ✅ `/omics-edge-snapshot` reports no DE regression

Plan 2 DoD is met. Handoff artifacts for Plan 3:
- `docs/superpowers/plans/2026-04-19-plan2-non-de-evidence-biller-2018.md` — this file (committed history)
- `.claude/skills/omics-edge-snapshot/snapshots/after_plan2.json`
- Current KG with Biller 2018's 7 DerivedMetric nodes + binding + measurement edges, ready for Plan 3's post-import Cypher (rank/percentile/bucket/significance, rollups) + KG validity tests.

---

## Self-review

**Spec coverage vs Plan 2 requirements (from slice spec §Plan 2):**

| Requirement | Task(s) |
|---|---|
| `multiomics_kg/adapters/observations_adapter.py` with `ObservationsAdapter` + `MultiObservationsAdapter` | Tasks 2, 5–12 |
| Wired into `create_knowledge_graph.py` after `MultiClusterAdapter` | Task 14 |
| Synthetic numeric-DM paperconfig fixture (~100 rows) under `tests/fixtures/non_de/` | Task 13 |
| Fixture wired via optional test-only `paperconfig_files.txt` | Task 13 (the fixture's own `paperconfig_files.txt`, read only by the unit test, not the production build) |
| One DerivedMetric node per (Experiment × metric_type) | Task 5 (core loop), Tasks 6–7 (value_kind branches) |
| Node ID = `derived_metric:{doi_short}:{entry_key}:{metric_type}` | Task 2 helper, Task 5 emission |
| 3 binding edges per DM: Publication → DM, Experiment → DM, DM → OrganismTaxon | Task 8 |
| Measurement edge type matches value_kind | Tasks 9 (boolean), 10 (categorical), 11 (numeric) |
| Denormalized fields computed from parent Experiment (spec invariant #9) | Task 4 `_denormalized_fields` helper |
| `rankable` / `has_p_value` paperconfig-declared for numeric; adapter forces "false" for boolean/categorical | Tasks 5, 6, 7 |
| `allowed_categories` inline from paperconfig for categorical | Task 6 |
| Boolean token parser hard-errors on unexpected tokens | Tasks 2, 3, 9 |
| Categorical hard-errors on out-of-set values | Task 10 |
| String sanitization via local `_clean_str()` | Task 2 |
| Uses `_resolved.csv` files via same probe pattern as omics_adapter | Task 2 (`_resolve_csv_path`), Task 8 (use in get_edges) |
| Tested against Biller 2018 real paperconfig | Task 15 |
| Tested against synthetic numeric-DM fixture | Task 13 |
| `uv run python create_knowledge_graph.py --test` succeeds | Task 14 Step 3 |
| `docker compose up -d` reaches import clean | Task 16 Step 2 |
| `Changes_expression_of` unchanged | Task 16 Step 3 |
| Adapter unit tests green | Task 16 Step 1 |

All requirements covered.

**Out of scope (confirmed — Plan 3 territory):**
- Post-import Cypher (rank/percentile/bucket, significance, rollups)
- KG validity tests
- `snapshot_data.json` regeneration for DerivedMetric nodes
- `docs/kg-changes/non-de-evidence-extension.md`
- Skill doc updates (paperconfig / cypher-queries / omics-edge-snapshot)

**Placeholder scan:** none found.

**Type consistency:** `_make_derived_metric_id`, `_parse_boolean_cell`, `_denormalized_fields`, `_emit_boolean_edges`, `_emit_categorical_edges`, `_emit_numeric_edges`, `_organism_lookup`, `_dm_entries` — consistent across all tasks. Edge labels exactly match schema (`publication_has_derived_metric`, `experiment_has_derived_metric`, `derived_metric_belongs_to_organism`, `derived_metric_flags_gene`, `derived_metric_classifies_gene`, `derived_metric_quantifies_gene`) and spec.

---

## Confirmed consistency with existing adapters

Prior to execution, these design decisions were verified against the live KG and existing adapters:

**Node ID format — raw string concatenation (no `normalize_curie`):**
- cluster_adapter's raw-string IDs like `clustering_analysis:ycae131:supp_table_3_darktolerant_clusters` land in the live Neo4j graph unchanged (verified via `MATCH (c:ClusteringAnalysis) RETURN c.id`).
- BioCypher + `neo4j-admin import` accept colon-separated multi-part IDs on project-specific prefixes (`cluster:`, `clustering_analysis:`, and now `derived_metric:`) with no prefix registration in bioregistry required.
- `schema_config.yaml`'s `preferred_id: derived_metric_id` is a schema-documentation label, not a validation gate — matches the `cluster_id` / `clustering_analysis_id` pattern.
- **Conclusion:** `normalize_curie` is NOT needed for DerivedMetric IDs. Plan follows cluster_adapter's raw-string style.

**Gene ID target format — `ncbigene:<locus_tag>` raw:**
- cluster_adapter emits raw `f"ncbigene:{gene_locus}"`; omics_adapter runs `normalize_curie("ncbigene:" + locus_tag)`. For registered prefixes (`ncbigene`) `normalize_curie` is a no-op — both paths produce identical strings, so they interoperate with CyanorakNcbi's Gene nodes (also normalize_curie'd). Plan follows cluster_adapter.

**Text sanitization — local `_clean_str(value) → value.replace("'", "^").replace("|", ",")`:**
- Matches cluster_adapter's local helper. Matches CLAUDE.md's stated convention ("Define this helper locally in each adapter"). `'` → `^` prevents single-quote-terminated strings from breaking CSV output; `|` → `,` prevents pipe-as-array-delimiter collisions.

---

## Follow-up refactors (out of scope for Plan 2)

Three DRY opportunities surfaced during plan review. None are blockers; each belongs in a separate post-Plan-2 commit (or a future refactor spec):

1. **DRY `_extract_doi` / `_load_pdf_cache` / `_resolve_csv_path` into `paperconfig_utils`.**
   These three helpers are now duplicated verbatim across `cluster_adapter.py` and (via Plan 2) `observations_adapter.py`. Drift risk is low today (both copies were written from the same cluster_adapter source), but any future change must be mirrored.
   **Proposed shape:**
   ```python
   # multiomics_kg/utils/paperconfig_utils.py
   def get_doi(config: dict, pdf_cache_path: Path | None = None) -> str: ...
   def load_pdf_cache(cache_path: Path | None = None) -> dict: ...
   def probe_resolved_csv(csv_path: str | Path) -> tuple[Path, bool]: ...
   ```
   Both adapters become thin wrappers. Touches 3 files + 3 test files.

2. **Factor test fixture builders (`_write_s4a_like_paperconfig`, `_write_s5_like_paperconfig`, `_write_numeric_dm_paperconfig`) into a shared module.**
   Plan 2 defines 3 similar test helpers. Cluster adapter tests have their own similar helpers (`_write_cluster_csv`, `_write_paperconfig`). Future observation-adapter additions will add more.
   **Proposed shape:** `tests/helpers/paperconfig_factory.py` exposing `build_derived_metrics_paperconfig(kind, metrics, ...)` plus `build_gene_clusters_paperconfig(...)`. Both test suites converge on one helper.

3. **Cosmetic: richer `name:` defaults on DerivedMetric emission.**
   If a paperconfig omits `name:` on a metric entry, adapter currently falls back to bare `metric_type` (Task 0 authors `name:` on all Biller 2018 entries, so in the production graph this fallback is never hit). Cluster_adapter's analogous fallback is `f"Clustering {entry_key}"` — disambiguated by entry_key. If a future paperconfig forgets `name:`, consider `f"{metric_type} ({entry_key})"` as a more informative default. Low priority — depends on whether validator ever makes `name:` required.

Tracking: add these as GitHub issues post-merge, or bundle into a "plan-2-followups" section in `plans/` if the repo pattern is to keep follow-ups as plan files.

---

## Completion status (2026-04-20)

All 15 coding tasks executed via `superpowers:subagent-driven-development`. 15 commits on `dev` from base `99595ec` through HEAD `d6bb17b`:

| Task | Commit(s) | Summary |
|---|---|---|
| 0 | `6075323` | Biller 2018 paperconfig polish (short entry_keys + `name:` per metric) |
| 1 | `b1aed37`, `ac0204c` | `Experiment.compartment` emission from omics_adapter |
| 2 | `b68ed10` | observations_adapter.py scaffold + helpers + boolean parser |
| 3 | `eafdd6d`, `28e1918` | Exhaustive boolean parser tests (incl. pd.NA, numeric, case) |
| 4 | `fd9b64e` | `_denormalized_fields` helper (spec invariant #9) |
| 5 | `ea89927`, `097db22` | `get_nodes` 3-branch impl + defensive paperconfig guards |
| 6 | `3f0e222` | Categorical DM node coverage test |
| 7 | `e1fd4a5` | Numeric DM node coverage test |
| 8 | `2738447`, `00d44c1` | `get_edges` binding edges + malformed-CSV guard |
| 9 | `52869c5` | `_emit_boolean_edges` → `derived_metric_flags_gene` |
| 10 | `8acf183` | `_emit_categorical_edges` → `derived_metric_classifies_gene` |
| 11 | `a5f656e` | `_emit_numeric_edges` → `derived_metric_quantifies_gene` |
| 12 | `eca9e55` | MultiObservationsAdapter organism-lookup tests |
| 13 | `c9ba43b` | Synthetic 100-row numeric-DM fixture + ingest test |
| 14 | `476f7c4` | Wire MultiObservationsAdapter into `create_knowledge_graph.py` |
| 15 | `d6bb17b` | Biller 2018 end-to-end count assertions (5 tests) |

**Test status:** 1,588 unit tests pass (pytest -m "not slow and not kg"). 44 tests in `tests/test_observations_adapter.py`, 5 in `tests/test_observations_adapter_biller2018.py`.

**Verified output** from `uv run python create_knowledge_graph.py --test` (Task 14):
- 7 DerivedMetric nodes
- 7 each of `publication_has_derived_metric` / `experiment_has_derived_metric` / `derived_metric_belongs_to_organism`
- 502 `derived_metric_flags_gene` (test-mode capped; **real count 4,160** per Task 15)
- 100 `derived_metric_classifies_gene` (test-mode capped; **real count 258** per Task 15)
- 0 `derived_metric_quantifies_gene` (Biller 2018 has no numeric)

### Deferred (Plan 2 DoD, user-deferred)

- `docker compose down && docker compose up -d --build` — full KG rebuild + import smoke
- Capture `.claude/skills/omics-edge-snapshot/snapshots/after_plan2.json` via the skill and compare against `after_biller2018_retrofit_removal.json` for DE parity check
- Query live graph to confirm `Changes_expression_of` count equals step-0 baseline (227,361)

These are low-risk (CSVs already validated; compartment emission is additive; DE adapter untouched). Safe to run when KG downtime is acceptable.

### Next: Plan 3

Scope (from slice spec §Plan 3):
- Post-import Cypher additions (`scripts/post-import.sh` + `scripts/post-import.cypher` byte-identical pair): rank/percentile/bucket on numeric DMs; significance derivation; analysis-node rollups (total_gene_count, growth_phases); Experiment/Publication/OrganismTaxon rollups; Gene routing counts; new scalar + full-text indexes
- KG validity tests under `tests/kg_validity/`
- Regenerated `tests/kg_validity/snapshot_data.json`
- `docs/kg-changes/non-de-evidence-extension.md` (downstream comms doc)
- Skill doc updates: `.claude/skills/paperconfig/SKILL.md`, `.claude/skills/cypher-queries/SKILL.md`, `.claude/skills/omics-edge-snapshot/SKILL.md`
- `CLAUDE.md` key-facts update

**Plan 3 expected-values reference (from Task 15 real counts):**
- `Derived_metric_flags_gene` edge count (Biller 2018): 4,160
- `Derived_metric_classifies_gene` edge count (Biller 2018): 258
- Synthetic fixture (Task 13): 100 × 3 metrics = 300 `Derived_metric_quantifies_gene` edges (if fixture wired into production for Plan 3 — currently test-only, NOT in production paperconfig_files.txt)

Plan 3's synthetic fixture decision: the fixture is currently test-only. Plan 3 needs to decide whether to wire it into the production build so numeric-DM post-import Cypher (rank/percentile/bucket/significance) has data to run on in the live graph — OR keep it test-only and only verify via Cypher unit tests. The spec says Plan 3 "relies on the Plan-2 fixture being in Neo4j" (line 240 of slice spec) — so likely wire it in for Plan 3, then remove post-deployment.
