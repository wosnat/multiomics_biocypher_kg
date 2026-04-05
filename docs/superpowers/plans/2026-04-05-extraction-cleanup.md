# Extraction Pipeline Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Clean up the hacky extraction pipeline into a maintainable module with quality verification.

**Architecture:** `extraction_utils.py` owns all file I/O and data logic. `extract.py` owns LLM calls, prompts, and CLI. Adapter delegates to `extraction_utils`. Old pipeline code (visual/semantic/merge/synthesis/validation/RunManager/RAG/review UI) deleted entirely.

**Tech Stack:** Python, OpenAI Responses API, Pydantic, pandas, pytest

**Spec:** `docs/superpowers/specs/2026-04-05-extraction-cleanup-design.md`

---

### Task 1: Create `extraction_utils.py` with file I/O functions

**Files:**
- Create: `multiomics_kg/extraction/cluster/extraction_utils.py`
- Create: `tests/test_extraction_utils.py`

- [ ] **Step 1: Write failing tests for load/save/list**

```python
# tests/test_extraction_utils.py
import json
import pytest
from pathlib import Path


def test_load_extraction(tmp_path):
    """Load clean-format JSON, get clusters dict."""
    from multiomics_kg.extraction.cluster.extraction_utils import load_extraction

    ext_dir = tmp_path / "cluster_extractions"
    ext_dir.mkdir()
    data = {
        "metadata": {"paper": "Test"},
        "clusters": {
            "1": {"id": "test_up", "name": "Cluster 1", "functional_description": "Transport genes"},
            "2": {"id": "test_down", "name": "Cluster 2", "functional_description": "Photosynthesis"},
        },
    }
    (ext_dir / "my_entry.json").write_text(json.dumps(data))

    result = load_extraction(tmp_path, "my_entry")
    assert len(result) == 2
    assert result["1"]["functional_description"] == "Transport genes"


def test_load_wrong_format(tmp_path):
    """Old format with stage2_results returns empty dict."""
    from multiomics_kg.extraction.cluster.extraction_utils import load_extraction

    ext_dir = tmp_path / "cluster_extractions"
    ext_dir.mkdir()
    old = {"stage2_results": {"1": {"name": "old"}}, "stage3_validation": {}}
    (ext_dir / "old_entry.json").write_text(json.dumps(old))

    result = load_extraction(tmp_path, "old_entry")
    assert result == {}


def test_load_missing_file(tmp_path):
    """Missing file returns empty dict, no crash."""
    from multiomics_kg.extraction.cluster.extraction_utils import load_extraction

    result = load_extraction(tmp_path, "nonexistent")
    assert result == {}


def test_save_load_roundtrip(tmp_path):
    """Save then load produces identical data."""
    from multiomics_kg.extraction.cluster.extraction_utils import load_extraction, save_extraction

    metadata = {"paper": "Test 2006", "model": "gpt-4.1-mini"}
    clusters = {
        "1": {
            "id": "test_up",
            "name": "Test cluster 1",
            "functional_description": "Genes involved in transport",
            "behavioral_description": "Upregulated early",
            "peak_time_hours": 6.0,
            "period_hours": None,
        },
    }
    save_extraction(tmp_path, "test_entry", metadata, clusters)
    loaded = load_extraction(tmp_path, "test_entry")
    assert loaded == clusters

    # Also check markdown was created
    md_path = tmp_path / "cluster_extractions" / "test_entry.md"
    assert md_path.exists()
    assert "Test cluster 1" in md_path.read_text()


def test_list_extraction_files(tmp_path):
    """List entry keys from cluster_extractions/ dir."""
    from multiomics_kg.extraction.cluster.extraction_utils import list_extraction_files

    ext_dir = tmp_path / "cluster_extractions"
    ext_dir.mkdir()
    (ext_dir / "entry_a.json").write_text("{}")
    (ext_dir / "entry_b.json").write_text("{}")
    (ext_dir / "entry_c.md").write_text("")  # md file should be ignored

    result = list_extraction_files(tmp_path)
    assert sorted(result) == ["entry_a", "entry_b"]


def test_get_cluster_data():
    """Simple dict lookup by key."""
    from multiomics_kg.extraction.cluster.extraction_utils import get_cluster_data

    clusters = {"1": {"name": "Cluster 1"}, "2": {"name": "Cluster 2"}}
    assert get_cluster_data(clusters, "1") == {"name": "Cluster 1"}
    assert get_cluster_data(clusters, 1) == {"name": "Cluster 1"}  # int key
    assert get_cluster_data(clusters, "99") == {}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_extraction_utils.py -v`
Expected: FAIL with `ModuleNotFoundError` (module doesn't exist yet)

- [ ] **Step 3: Implement `extraction_utils.py`**

```python
# multiomics_kg/extraction/cluster/extraction_utils.py
"""Extraction file I/O and data utilities.

Single source of truth for extraction file locations and JSON structure.
Used by the adapter, extract.py, and report/verify tooling.
"""
import json
import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

EXTRACTIONS_DIR = "cluster_extractions"


# ── File I/O ──


def load_extraction(paper_dir: Path, entry_key: str) -> dict[str, dict]:
    """Load cluster extraction data from JSON.

    Returns the ``clusters`` dict: {cluster_key: {id, name, functional_description, ...}}.
    Returns {} if file is missing or has wrong format.
    """
    json_path = Path(paper_dir) / EXTRACTIONS_DIR / f"{entry_key}.json"
    if not json_path.exists():
        return {}
    try:
        data = json.loads(json_path.read_text())
    except Exception:
        logger.warning("Failed to parse extraction JSON: %s", json_path)
        return {}
    if "clusters" not in data:
        logger.warning("Extraction JSON missing 'clusters' key: %s", json_path)
        return {}
    return data["clusters"]


def save_extraction(
    paper_dir: Path,
    entry_key: str,
    metadata: dict,
    clusters: dict[str, dict],
) -> Path:
    """Write extraction results as JSON + markdown summary.

    Creates ``cluster_extractions/`` dir if needed. Returns path to JSON file.
    """
    ext_dir = Path(paper_dir) / EXTRACTIONS_DIR
    ext_dir.mkdir(parents=True, exist_ok=True)

    # JSON
    output = {"metadata": metadata, "clusters": clusters}
    json_path = ext_dir / f"{entry_key}.json"
    json_path.write_text(json.dumps(output, indent=2, default=str))

    # Markdown summary
    md_lines = [f"# {metadata.get('paper', '')} — {entry_key}\n"]
    for key in sorted(clusters, key=lambda x: (not x.isdigit(), x)):
        c = clusters[key]
        direction = c.get("direction", "")
        assessment = c.get("self_assessment", "")
        md_lines.append(f"## Cluster {key} | {direction} | {assessment}\n")
        md_lines.append(f"**Name:** {c.get('name', '')}")
        md_lines.append(f"**Enrichment:** {c.get('enrichment_category', '')} "
                        f"(p={c.get('enrichment_pvalue', 'N/A')})")
        md_lines.append(f"**Functional:** {c.get('functional_description', '')}\n")
        md_lines.append(f"**Behavioral:** {c.get('behavioral_description', '')}\n")
        notes = c.get("confidence_notes", "")
        if notes:
            md_lines.append(f"**Notes:** {notes}\n")
        quotes = c.get("supporting_quotes", [])
        if quotes:
            md_lines.append("**Quotes:**")
            for q in quotes:
                md_lines.append(f"- [{q.get('location', '')}] {q.get('quote', '')}")
            md_lines.append("")
    md_path = ext_dir / f"{entry_key}.md"
    md_path.write_text("\n".join(md_lines))

    return json_path


def get_cluster_data(clusters: dict, key) -> dict:
    """Get data for a single cluster by key. Returns {} if not found."""
    return clusters.get(str(key), {})


def list_extraction_files(paper_dir: Path) -> list[str]:
    """List entry keys that have extraction JSONs in paper_dir/cluster_extractions/."""
    ext_dir = Path(paper_dir) / EXTRACTIONS_DIR
    if not ext_dir.is_dir():
        return []
    return sorted(p.stem for p in ext_dir.glob("*.json"))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_extraction_utils.py -v`
Expected: all 6 tests PASS

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/extraction/cluster/extraction_utils.py tests/test_extraction_utils.py
git commit -m "feat: add extraction_utils with load/save/list for clean JSON format"
```

---

### Task 2: Add cluster key matching and CSV loading to `extraction_utils.py`

**Files:**
- Modify: `multiomics_kg/extraction/cluster/extraction_utils.py`
- Modify: `tests/test_extraction_utils.py`

- [ ] **Step 1: Write failing tests for match_cluster_keys and load_cluster_summaries**

Append to `tests/test_extraction_utils.py`:

```python
def test_match_cluster_keys_numeric():
    """Numeric keys matched by 'cluster N' in name."""
    from multiomics_kg.extraction.cluster.extraction_utils import match_cluster_keys

    parsed = [
        {"id": "mit_up_1", "name": "MIT9313 cluster 1 (up, transport)"},
        {"id": "mit_down_6", "name": "MIT9313 cluster 6 (down, photosynthesis)"},
        {"id": "mit_down_7", "name": "MIT9313 cluster 7 (down, translation)"},
    ]
    expected_keys = {"1", "2", "3", "4", "5", "6", "7"}
    matched, unmatched = match_cluster_keys(parsed, expected_keys)
    assert "1" in matched
    assert "6" in matched
    assert "7" in matched
    assert matched["1"]["id"] == "mit_up_1"
    assert len(unmatched) == 0


def test_match_cluster_keys_alpha():
    """Alpha keys like HEG, LEG matched case-insensitively."""
    from multiomics_kg.extraction.cluster.extraction_utils import match_cluster_keys

    parsed = [
        {"id": "med4_heg", "name": "MED4 cluster HEG (up, expression)"},
        {"id": "med4_leg", "name": "MED4 cluster LEG (down, expression)"},
    ]
    expected_keys = {"HEG", "LEG", "MEG"}
    matched, unmatched = match_cluster_keys(parsed, expected_keys)
    assert "HEG" in matched
    assert "LEG" in matched
    assert len(unmatched) == 0


def test_match_cluster_keys_composite_positional():
    """Composite keys fall back to positional matching when names don't match."""
    from multiomics_kg.extraction.cluster.extraction_utils import match_cluster_keys

    parsed = [
        {"id": "x_1", "name": "NATL2A cluster 1 (up, periodic)"},
        {"id": "x_2", "name": "NATL2A cluster 2 (down, periodic)"},
    ]
    expected_keys = {"coculture_LD", "coculture_darkness"}
    matched, unmatched = match_cluster_keys(parsed, expected_keys)
    # Positional fallback: 2 parsed, 2 expected → match by position
    assert len(matched) == 2
    assert len(unmatched) == 0


def test_match_cluster_keys_unmatched():
    """Unmatched extractions returned in unmatched list when counts differ."""
    from multiomics_kg.extraction.cluster.extraction_utils import match_cluster_keys

    parsed = [
        {"id": "x_1", "name": "cluster 1"},
        {"id": "x_99", "name": "cluster 99"},
    ]
    expected_keys = {"1", "2", "3"}
    matched, unmatched = match_cluster_keys(parsed, expected_keys)
    assert "1" in matched
    assert len(unmatched) == 1
    assert unmatched[0]["id"] == "x_99"


def test_load_cluster_summaries(tmp_path):
    """Load cluster summaries from CSV."""
    from multiomics_kg.extraction.cluster.extraction_utils import load_cluster_summaries

    csv_path = tmp_path / "clusters.csv"
    csv_path.write_text("gene_id,cluster\nPMT001,1\nPMT002,1\nPMT003,2\n")

    table_config = {
        "filename": str(csv_path),
        "gene_id_col": "gene_id",
        "cluster_col": "cluster",
    }
    result = load_cluster_summaries(table_config)
    assert len(result) == 2
    assert result["1"]["gene_count"] == 2
    assert result["2"]["gene_count"] == 1
    assert "PMT001" in result["1"]["sample_genes"]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_extraction_utils.py::test_match_cluster_keys_numeric tests/test_extraction_utils.py::test_load_cluster_summaries -v`
Expected: FAIL with `ImportError`

- [ ] **Step 3: Implement match_cluster_keys, load_cluster_summaries, find_all_entries**

Append to `multiomics_kg/extraction/cluster/extraction_utils.py`:

```python
import re

from multiomics_kg.utils.paperconfig_utils import (
    load_all_paperconfigs,
    iter_cluster_tables,
)


# ── Data loading ──


def find_all_entries() -> list[tuple[Path, str, dict, dict]]:
    """Discover all gene_clusters entries across all paperconfigs.

    Returns list of (paper_dir, entry_key, table_config, pub_config).
    """
    entries = []
    for pc_path, pc in load_all_paperconfigs():
        paper_dir = pc_path.parent
        pub = pc.get("publication", {})
        for entry_key, table in iter_cluster_tables(pc):
            entries.append((paper_dir, entry_key, table, pub))
    return entries


def load_cluster_summaries(table_config: dict) -> dict[str, dict]:
    """Load cluster gene counts and sample genes from CSV.

    Returns {cluster_key: {gene_count: int, sample_genes: list[str]}}.
    """
    csv_path = Path(table_config["filename"])
    gene_id_col = table_config["gene_id_col"]
    cluster_col = table_config["cluster_col"]
    skip_rows = table_config.get("skip_rows", 0)
    df = pd.read_csv(csv_path, skiprows=skip_rows if skip_rows else None)

    clusters = {}
    for cid, group in df.groupby(cluster_col):
        ckey = _cluster_val_to_str(cid)
        if not ckey:
            continue
        genes = group[gene_id_col].dropna().tolist()
        clusters[ckey] = {
            "gene_count": len(genes),
            "sample_genes": [str(g) for g in genes[:5]],
        }
    return clusters


def _cluster_val_to_str(val) -> str:
    """Convert cluster column value to clean string key."""
    if pd.isna(val):
        return ""
    if isinstance(val, float) and val == int(val):
        return str(int(val))
    return str(val).strip()


# ── Cluster key matching ──


def match_cluster_keys(
    parsed_clusters: list[dict],
    expected_keys: set[str],
) -> tuple[dict[str, dict], list[dict]]:
    """Map model output to actual cluster keys.

    Tries: name regex → id suffix → case-insensitive → positional fallback.
    Returns (matched: {key: extraction_dict}, unmatched: [extraction_dict]).
    """
    remaining_keys = set(expected_keys)
    matched = {}
    unmatched_pass1 = []

    for c in parsed_clusters:
        key = _try_match_key(c, remaining_keys)
        if key:
            matched[key] = c
            remaining_keys.discard(key)
        else:
            unmatched_pass1.append(c)

    # Positional fallback: if unmatched count == remaining keys count
    if unmatched_pass1 and len(unmatched_pass1) == len(remaining_keys):
        remaining_sorted = sorted(remaining_keys, key=lambda x: (not x.isdigit(), x))
        for c, key in zip(unmatched_pass1, remaining_sorted):
            matched[key] = c
        return matched, []

    return matched, unmatched_pass1


def _try_match_key(extraction: dict, expected_keys: set[str]) -> str | None:
    """Try to match a single extraction to a cluster key."""
    name = extraction.get("name", "")
    ext_id = extraction.get("id", "")

    # Try "cluster KEY" in name
    m = re.search(r"cluster\s+(\S+)", name, re.IGNORECASE)
    if m:
        candidate = m.group(1).rstrip(",.):")
        if candidate in expected_keys:
            return candidate
        for k in expected_keys:
            if k.lower() == candidate.lower():
                return k

    # Try suffix of id
    m = re.search(r"_([^_]+)$", ext_id)
    if m:
        candidate = m.group(1)
        if candidate in expected_keys:
            return candidate
        for k in expected_keys:
            if k.lower() == candidate.lower():
                return k

    return None
```

- [ ] **Step 4: Run all tests**

Run: `uv run pytest tests/test_extraction_utils.py -v`
Expected: all 11 tests PASS

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/extraction/cluster/extraction_utils.py tests/test_extraction_utils.py
git commit -m "feat: add cluster key matching and CSV loading to extraction_utils"
```

---

### Task 3: Update adapter to use `extraction_utils`

**Files:**
- Modify: `multiomics_kg/adapters/cluster_adapter.py:62-118` (replace `_load_extraction_json` and `_get_extraction_cluster_data`)
- Modify: `multiomics_kg/adapters/cluster_adapter.py:231,271` (call sites)
- Modify: `tests/test_cluster_adapter.py` (update/remove RunManager tests)

- [ ] **Step 1: Write failing test for adapter with new format**

Add to `tests/test_cluster_adapter.py` (replace the old RunManager-based tests):

```python
def test_adapter_loads_new_extraction_format(tmp_path):
    """Adapter reads from cluster_extractions/{entry_key}.json."""
    import json

    from multiomics_kg.extraction.cluster.extraction_utils import save_extraction

    # Write extraction in new format
    metadata = {"paper": "Test"}
    clusters = {
        "1": {
            "id": "test_up_transport",
            "name": "Test cluster 1 (up, transport)",
            "functional_description": "Transport genes upregulated",
            "behavioral_description": "Rapid upregulation",
            "peak_time_hours": 6.0,
            "period_hours": None,
        },
    }
    save_extraction(tmp_path, "test_entry", metadata, clusters)

    from multiomics_kg.extraction.cluster.extraction_utils import load_extraction, get_cluster_data

    loaded = load_extraction(tmp_path, "test_entry")
    ext = get_cluster_data(loaded, "1")
    assert ext["functional_description"] == "Transport genes upregulated"
    assert ext["name"] == "Test cluster 1 (up, transport)"
    assert ext["peak_time_hours"] == 6.0
```

- [ ] **Step 2: Run test to verify it passes** (this tests extraction_utils, should already pass)

Run: `uv run pytest tests/test_cluster_adapter.py::test_adapter_loads_new_extraction_format -v`
Expected: PASS

- [ ] **Step 3: Replace adapter loading code**

In `multiomics_kg/adapters/cluster_adapter.py`, replace lines 62-118 (both `_load_extraction_json` and `_get_extraction_cluster_data`) with:

```python
from multiomics_kg.extraction.cluster.extraction_utils import load_extraction, get_cluster_data
```

Then update the call sites:
- Line 231: change `extraction = _load_extraction_json(self._paperconfig_dir, entry_key)` → `extraction = load_extraction(self._paperconfig_dir, entry_key)`
- Line 271: change `ext_data = _get_extraction_cluster_data(extraction, cluster_key)` → `ext_data = get_cluster_data(extraction, cluster_key)`

- [ ] **Step 4: Remove old tests that reference RunManager and deleted functions**

In `tests/test_cluster_adapter.py`:
- Remove `_load_extraction_json` and `_get_extraction_cluster_data` from the import block (line 27-28)
- Remove `test_load_extraction_from_cache_dir` (the RunManager-based test, ~line 648-675)
- Remove `test_extraction_review_reject_overrides_verdict` (~line 685-692)
- Remove `test_extraction_review_edit_overrides_fields` (~line 700-719)

- [ ] **Step 5: Run adapter tests**

Run: `uv run pytest tests/test_cluster_adapter.py -v`
Expected: PASS (remaining tests still work)

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/adapters/cluster_adapter.py tests/test_cluster_adapter.py
git commit -m "refactor: adapter uses extraction_utils, remove verdict gating"
```

---

### Task 4: Move extraction JSONs to new location + format

**Files:**
- Create: `scripts/migrate_extractions.py` (one-time migration script)

- [ ] **Step 1: Write migration script**

```python
#!/usr/bin/env python3
"""One-time: migrate extraction JSONs to new format + location.

Reads: data/<paper>/cluster_extraction_{entry}.json (old format with stage2_results)
Writes: data/<paper>/cluster_extractions/{entry}.json (new format with clusters)
Deletes: old files after successful migration
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from multiomics_kg.extraction.cluster.extraction_utils import save_extraction


def migrate():
    old_files = sorted(Path("data").rglob("cluster_extraction_*.json"))
    old_files = [f for f in old_files if "pilot" not in f.name and "cluster_extractions" not in str(f)]

    print(f"Found {len(old_files)} files to migrate\n")

    for old_path in old_files:
        paper_dir = old_path.parent
        # Extract entry_key from filename: cluster_extraction_{entry_key}.json
        entry_key = old_path.stem.replace("cluster_extraction_", "")

        data = json.loads(old_path.read_text())

        # Old format: stage2_results; new format: clusters
        if "clusters" in data:
            clusters = data["clusters"]
        elif "stage2_results" in data:
            clusters = data["stage2_results"]
        else:
            print(f"  SKIP {old_path.name} — no clusters or stage2_results key")
            continue

        metadata = data.get("metadata", {})

        new_path = save_extraction(paper_dir, entry_key, metadata, clusters)
        print(f"  {old_path.name} -> {new_path.relative_to(paper_dir)}")

        # Delete old files
        old_path.unlink()
        old_md = old_path.with_suffix(".md")
        if old_md.exists():
            old_md.unlink()
            print(f"    Deleted {old_md.name}")

    print("\nDone.")


if __name__ == "__main__":
    migrate()
```

- [ ] **Step 2: Run migration**

Run: `uv run python scripts/migrate_extractions.py`
Expected: 15 files migrated, old files deleted

- [ ] **Step 3: Verify new files exist and adapter can read them**

Run:
```bash
ls data/Prochlorococcus/papers_and_supp/tolonen\ 2006/cluster_extractions/
uv run python -c "
import sys; sys.path.insert(0, '.')
from pathlib import Path
from multiomics_kg.extraction.cluster.extraction_utils import load_extraction
clusters = load_extraction(Path('data/Prochlorococcus/papers_and_supp/tolonen 2006'), 'mit9313_kmeans_nstarvation')
print(f'Loaded {len(clusters)} clusters')
print(f'Cluster 6: {clusters[\"6\"][\"name\"][:50]}')
"
```
Expected: `mit9313_kmeans_nstarvation.json` and `.md` in `cluster_extractions/`. 7 clusters loaded. Cluster 6 has "photosynthesis" in name.

- [ ] **Step 4: Verify no old files remain**

Run: `find data -name "cluster_extraction_*.json" -not -path "*/cluster_extractions/*" -not -name "pilot*"`
Expected: no output (all migrated)

- [ ] **Step 5: Commit**

```bash
git add data/*/papers_and_supp/*/cluster_extractions/ scripts/migrate_extractions.py
git rm data/*/papers_and_supp/*/cluster_extraction_*.json data/*/papers_and_supp/*/cluster_extraction_*.md 2>/dev/null || true
git commit -m "refactor: migrate extraction JSONs to cluster_extractions/ with clean format"
```

---

### Task 5: Create clean `extract.py` module

**Files:**
- Create: `multiomics_kg/extraction/cluster/extract.py`
- Delete: `scripts/run_extraction.py`
- Create: `tests/test_extract.py`

- [ ] **Step 1: Write failing tests for prompt and report functions**

```python
# tests/test_extract.py
import json
import pytest
from pathlib import Path


def test_build_context_block():
    """Context block includes organism, method, treatment."""
    from multiomics_kg.extraction.cluster.extract import build_context_block

    table = {
        "name": "MIT9313 N-starvation",
        "organism": "Prochlorococcus MIT9313",
        "cluster_method": "K-means (K=7)",
        "cluster_type": "response_pattern",
        "treatment": "N-starvation time course",
        "omics_type": "MICROARRAY",
    }
    result = build_context_block(table)
    assert "Prochlorococcus MIT9313" in result
    assert "K-means (K=7)" in result
    assert "N-starvation time course" in result
    assert "MICROARRAY" in result


def test_format_cluster_summaries():
    """Summary text has all cluster keys."""
    from multiomics_kg.extraction.cluster.extract import format_cluster_summaries

    clusters = {
        "1": {"gene_count": 7, "sample_genes": ["PMT001", "PMT002"]},
        "6": {"gene_count": 81, "sample_genes": ["PMT100"]},
    }
    result = format_cluster_summaries(clusters)
    assert "Cluster 1:" in result
    assert "Cluster 6:" in result
    assert "7 genes" in result
    assert "81 genes" in result


def test_generate_report_stable_order(tmp_path):
    """Report generated twice is identical (no ordering jitter)."""
    from multiomics_kg.extraction.cluster.extract import generate_report
    from multiomics_kg.extraction.cluster.extraction_utils import save_extraction

    metadata = {"paper": "Test 2006"}
    clusters = {
        "2": {"id": "t_2", "name": "Cluster 2", "direction": "down",
              "self_assessment": "medium", "functional_description": "Photo",
              "behavioral_description": "Down", "confidence_notes": ""},
        "1": {"id": "t_1", "name": "Cluster 1", "direction": "up",
              "self_assessment": "high", "functional_description": "Transport",
              "behavioral_description": "Up", "confidence_notes": ""},
    }
    save_extraction(tmp_path, "test_entry", metadata, clusters)

    entries = [(tmp_path, "test_entry", {}, {"papername": "Test 2006"})]
    report1 = generate_report(entries)
    report2 = generate_report(entries)
    assert report1 == report2
    assert report1.index("Cluster 1") < report1.index("Cluster 2")
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_extract.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Implement `extract.py`**

Create `multiomics_kg/extraction/cluster/extract.py` — move and clean up from `scripts/run_extraction.py`:

```python
# multiomics_kg/extraction/cluster/extract.py
"""Cluster description extraction via OpenAI Responses API.

Usage:
    uv run python -m multiomics_kg.extraction.cluster.extract             # extract all
    uv run python -m multiomics_kg.extraction.cluster.extract --report    # generate report
    uv run python -m multiomics_kg.extraction.cluster.extract --verify    # report + checks
    uv run python -m multiomics_kg.extraction.cluster.extract --dry-run   # show plan
    uv run python -m multiomics_kg.extraction.cluster.extract --paper "Tolonen 2006"
    uv run python -m multiomics_kg.extraction.cluster.extract --entry mit9313_kmeans_nstarvation --force
"""
import argparse
import json
import logging
import os
import re
import time
from datetime import datetime
from pathlib import Path
from typing import Literal, Optional

from pydantic import BaseModel

from multiomics_kg.extraction.cluster.extraction_utils import (
    find_all_entries,
    get_cluster_data,
    list_extraction_files,
    load_cluster_summaries,
    load_extraction,
    match_cluster_keys,
    save_extraction,
)

logger = logging.getLogger(__name__)


# ── Pydantic schemas ──


class SupportingQuote(BaseModel):
    quote: str
    location: str


class ClusterExtraction(BaseModel):
    id: str
    name: str
    functional_description: str
    behavioral_description: str
    peak_time_hours: Optional[float]
    period_hours: Optional[float]
    direction: Literal["up", "down", "mixed", "not_described"]
    enrichment_category: str
    enrichment_pvalue: Optional[float]
    enrichment_significant: bool
    confidence_notes: str
    supporting_quotes: list[SupportingQuote]
    self_assessment: Literal["high", "medium", "low"]
    assessment_notes: str


class AnalysisExtraction(BaseModel):
    clusters: list[ClusterExtraction]


# ── Prompt layer ──


DEVELOPER_MSG_TEMPLATE = """\
You are extracting structured descriptions of gene expression clusters from \
a scientific paper.

{context_block}

For each cluster, extract all fields in the output schema. Key rules:
- Read FIGURES (heatmaps, time-course plots) as PRIMARY source for direction/timing.
- "not described in paper" is always better than guessing.
- Only cite genes mentioned by name in the paper (not locus tags like PMM*, PMT*, etc).
- Do NOT include treatment conditions in descriptions — those live on the analysis node.
- functional_description: 2-3 sentences. behavioral_description: 1-2 sentences.
- self_assessment: your confidence. assessment_notes: what you're uncertain about.
- Each cluster must have a unique id in snake_case: {{organism_short}}_{{direction}}_{{theme}}.
- name format: "{{Organism}} cluster {{KEY}} ({{direction}}, {{theme}})" — under 60 chars.
  Use the EXACT cluster key from the list below.
- For enrichment: use p-values from the paper/figures, max 3 decimal places or scientific notation.
- Max 3-5 named genes per cluster description.

CRITICAL: You MUST extract EXACTLY {n_clusters} clusters, one for each cluster key \
listed below. Use the EXACT cluster keys as they appear — do NOT renumber, skip, \
or merge clusters. Every key must appear exactly once in your output.

{cluster_summaries}
"""


def build_context_block(table: dict) -> str:
    """Build context block from paperconfig entry."""
    parts = [
        f"Analysis: {table.get('name', '')}",
        f"Organism: {table.get('organism', '')}",
        f"Clustering: {table.get('cluster_method', '')}",
        f"Type: {table.get('cluster_type', '')}",
        f"Treatment: {table.get('treatment', '')}",
    ]
    if table.get("experimental_context"):
        parts.append(f"Context: {table['experimental_context']}")
    if table.get("omics_type"):
        parts.append(f"Omics: {table['omics_type']}")
    if table.get("figure_hint"):
        parts.append(f"Key figures: {table['figure_hint']}")
    if table.get("time_points"):
        parts.append(f"Time points (hours): {table['time_points']}")
    return "\n".join(parts)


def format_cluster_summaries(clusters: dict[str, dict]) -> str:
    """Format cluster summaries for the prompt."""
    lines = []
    for key in sorted(clusters, key=lambda x: (not x.isdigit(), x)):
        info = clusters[key]
        sample = ", ".join(info["sample_genes"][:3])
        lines.append(f"Cluster {key}: {info['gene_count']} genes (sample: {sample})")
    return "\n".join(lines)


# ── LLM layer ──


def upload_pdf(client, pdf_path: Path) -> str:
    """Upload PDF via Files API, return file_id."""
    f = client.files.create(file=open(pdf_path, "rb"), purpose="user_data")
    return f.id


def extract_analysis(
    client,
    file_id: str,
    table_config: dict,
    cluster_summaries: dict[str, dict],
    model: str = "gpt-4.1-mini",
    flex: bool = False,
) -> tuple[AnalysisExtraction, dict]:
    """Run extraction for one analysis. Returns (parsed, usage_dict)."""
    ctx = build_context_block(table_config)
    summaries = format_cluster_summaries(cluster_summaries)

    dev_msg = DEVELOPER_MSG_TEMPLATE.format(
        context_block=ctx,
        n_clusters=len(cluster_summaries),
        cluster_summaries=summaries,
    )

    kwargs = dict(
        model=model,
        temperature=0,
        input=[
            {"role": "developer", "content": dev_msg},
            {"role": "user", "content": [
                {"type": "input_file", "file_id": file_id},
                {"type": "input_text", "text": f"Extract descriptions for all {len(cluster_summaries)} clusters."},
            ]},
        ],
        text_format=AnalysisExtraction,
    )
    if flex:
        kwargs["service_tier"] = "flex"

    t0 = time.time()
    resp = client.responses.parse(**kwargs)
    elapsed = time.time() - t0

    usage = {
        "input_tokens": resp.usage.input_tokens,
        "output_tokens": resp.usage.output_tokens,
        "duration_sec": elapsed,
    }
    return resp.output[0].content[0].parsed, usage


def extract_paper(
    client,
    file_id: str,
    tables_and_summaries: list[tuple[dict, dict[str, dict]]],
    model: str = "gpt-4.1-mini",
    flex: bool = False,
) -> tuple[AnalysisExtraction, dict]:
    """Extract all analyses for one paper in a single call.

    Used for multi-organism papers (e.g. Tolonen) where per-analysis
    extraction confuses the model.
    """
    context_parts = []
    all_summaries = []
    total_clusters = 0
    for table_config, cluster_summaries in tables_and_summaries:
        context_parts.append(build_context_block(table_config))
        all_summaries.append(format_cluster_summaries(cluster_summaries))
        total_clusters += len(cluster_summaries)

    ctx = "\n\n".join(context_parts)
    summaries = "\n\n".join(all_summaries)

    dev_msg = DEVELOPER_MSG_TEMPLATE.format(
        context_block=ctx,
        n_clusters=total_clusters,
        cluster_summaries=summaries,
    )

    kwargs = dict(
        model=model,
        temperature=0,
        input=[
            {"role": "developer", "content": dev_msg},
            {"role": "user", "content": [
                {"type": "input_file", "file_id": file_id},
                {"type": "input_text", "text": f"Extract descriptions for all {total_clusters} clusters."},
            ]},
        ],
        text_format=AnalysisExtraction,
    )
    if flex:
        kwargs["service_tier"] = "flex"

    t0 = time.time()
    resp = client.responses.parse(**kwargs)
    elapsed = time.time() - t0

    usage = {
        "input_tokens": resp.usage.input_tokens,
        "output_tokens": resp.usage.output_tokens,
        "duration_sec": elapsed,
    }
    return resp.output[0].content[0].parsed, usage


# ── Report layer ──


def generate_report(entries: list[tuple[Path, str, dict, dict]]) -> str:
    """Generate diff-friendly markdown report from existing extraction JSONs."""
    lines = ["# Cluster Extraction Report\n"]

    for paper_dir, entry_key, table_config, pub in sorted(entries, key=lambda e: (pub_name(e[3]), e[1])):
        clusters = load_extraction(paper_dir, entry_key)
        if not clusters:
            continue

        paper = pub_name(pub)
        lines.append(f"## {paper} / {entry_key}\n")

        for key in sorted(clusters, key=lambda x: (not x.isdigit(), x)):
            c = clusters[key]
            direction = c.get("direction", "")
            assessment = c.get("self_assessment", "")
            lines.append(f"### Cluster {key} | {direction} | {assessment}")
            lines.append(f"**Name:** {c.get('name', '')}")
            enrich = c.get("enrichment_category", "")
            pval = c.get("enrichment_pvalue", "N/A")
            sig = c.get("enrichment_significant", "")
            lines.append(f"**Enrichment:** {enrich} (p={pval}, sig={sig})")
            lines.append(f"**Functional:** {c.get('functional_description', '')}")
            lines.append(f"**Behavioral:** {c.get('behavioral_description', '')}")
            notes = c.get("confidence_notes", "")
            if notes:
                lines.append(f"**Notes:** {notes}")
            lines.append("")

    return "\n".join(lines)


def verify_quality(entries: list[tuple[Path, str, dict, dict]]) -> list[str]:
    """Run programmatic quality checks. Returns list of warning strings."""
    warnings = []

    for paper_dir, entry_key, table_config, pub in entries:
        clusters = load_extraction(paper_dir, entry_key)
        if not clusters:
            continue
        paper = pub_name(pub)
        prefix = f"[{paper} / {entry_key}"

        # Check 1: duplicate ids within analysis
        ids_seen = {}
        for key, c in clusters.items():
            cid = c.get("id", "")
            if cid in ids_seen:
                warnings.append(f"{prefix} / cluster {key}] duplicate id: {cid} (also used by cluster {ids_seen[cid]})")
            ids_seen[cid] = key

        # Check 2: locus tags in descriptions
        locus_pat = re.compile(r"\b(PMM\d{3,}|PMT\d{3,}|P9301_\d+|tll\d{3,}|SY28_\d+|BSR22_\d+)\b")
        for key, c in clusters.items():
            for field in ("functional_description", "behavioral_description"):
                text = c.get(field, "")
                m = locus_pat.search(text)
                if m:
                    warnings.append(f"{prefix} / cluster {key}] locus tag in {field}: {m.group()}")

        # Check 3: empty direction with non-empty description
        for key, c in clusters.items():
            has_desc = len(c.get("functional_description", "")) > 20
            direction = c.get("direction", "")
            if has_desc and not direction:
                warnings.append(f"{prefix} / cluster {key}] has description but empty direction")

    return warnings


def pub_name(pub: dict) -> str:
    """Extract paper name from publication config."""
    return pub.get("papername", "Unknown")


# ── CLI ──


def main():
    from dotenv import load_dotenv
    load_dotenv()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    parser = argparse.ArgumentParser(description="Cluster extraction pipeline")
    parser.add_argument("--paper", help="Filter by paper name (partial match)")
    parser.add_argument("--entry", help="Filter by entry key (exact)")
    parser.add_argument("--model", default=os.environ.get("CLUSTER_EXTRACTION_MODEL", "gpt-4.1-mini"))
    parser.add_argument("--flex", action="store_true", help="Flex processing (50%% cheaper)")
    parser.add_argument("--force", action="store_true", help="Overwrite existing")
    parser.add_argument("--dry-run", action="store_true", help="Show what would run")
    parser.add_argument("--report", action="store_true", help="Generate report (no API calls)")
    parser.add_argument("--verify", action="store_true", help="Run quality checks (implies --report)")
    args = parser.parse_args()

    if args.verify:
        args.report = True

    entries = find_all_entries()

    # Filter
    if args.paper:
        entries = [(d, k, t, p) for d, k, t, p in entries
                   if args.paper.lower() in p.get("papername", "").lower()]
    if args.entry:
        entries = [(d, k, t, p) for d, k, t, p in entries if k == args.entry]

    if args.report:
        report = generate_report(entries)
        if args.verify:
            warnings = verify_quality(entries)
            if warnings:
                report += "\n## Warnings\n\n"
                report += "\n".join(f"- {w}" for w in warnings) + "\n"
            else:
                report += "\n## Warnings\n\nNo issues found.\n"
        report_path = Path("data/cluster_extraction_report.md")
        report_path.write_text(report)
        logger.info(f"Report written to {report_path}")
        if args.verify:
            logger.info(f"{len(warnings)} warnings")
        return

    if not entries:
        logger.warning("No matching entries found")
        return

    if args.dry_run:
        print(f"\nWould process {len(entries)} entries:\n")
        for paper_dir, entry_key, table, pub in entries:
            summaries = load_cluster_summaries(table)
            exists = entry_key in list_extraction_files(paper_dir)
            status = "EXISTS" if exists else "NEW"
            print(f"  {pub_name(pub)} / {entry_key}: {len(summaries)} clusters [{status}]")
        return

    # Run extraction
    from openai import OpenAI
    client = OpenAI()
    total_in = total_out = total_clusters = 0

    # Group by paper for PDF reuse
    by_paper: dict[Path, list] = {}
    for paper_dir, entry_key, table, pub in entries:
        by_paper.setdefault(paper_dir, []).append((entry_key, table, pub))

    for paper_dir, group in by_paper.items():
        paper = pub_name(group[0][2])

        # Find PDF
        pdf_path_str = group[0][2].get("papermainpdf", "")
        if pdf_path_str:
            pdf_path = Path(pdf_path_str)
            if not pdf_path.is_absolute():
                pdf_path = Path.cwd() / pdf_path
        else:
            pdfs = list(paper_dir.glob("*.pdf"))
            pdf_path = pdfs[0] if pdfs else None

        if not pdf_path or not pdf_path.exists():
            logger.warning(f"No PDF for {paper}, skipping")
            continue

        logger.info(f"Paper: {paper} ({len(group)} entries)")
        file_id = upload_pdf(client, pdf_path)

        for entry_key, table, pub in group:
            if entry_key in list_extraction_files(paper_dir) and not args.force:
                logger.info(f"  {entry_key}: exists, skipping")
                continue

            summaries = load_cluster_summaries(table)
            logger.info(f"  {entry_key}: {len(summaries)} clusters")

            try:
                parsed, usage = extract_analysis(
                    client, file_id, table, summaries,
                    model=args.model, flex=args.flex,
                )

                expected_keys = set(summaries.keys())
                matched, unmatched = match_cluster_keys(
                    [c.model_dump() for c in parsed.clusters], expected_keys,
                )

                for c in unmatched:
                    logger.warning(f"    Unmatched: {c.get('name', '?')}")
                missing = expected_keys - set(matched.keys())
                if missing:
                    logger.warning(f"    Missing: {sorted(missing)}")

                metadata = {
                    "paper": paper,
                    "doi": pub.get("doi"),
                    "organism": table.get("organism", ""),
                    "entry_key": entry_key,
                    "model": args.model,
                    "extracted_at": datetime.now().isoformat(),
                    "input_tokens": usage["input_tokens"],
                    "output_tokens": usage["output_tokens"],
                }
                save_extraction(paper_dir, entry_key, metadata, matched)
                total_in += usage["input_tokens"]
                total_out += usage["output_tokens"]
                total_clusters += len(matched)
                logger.info(f"    {len(matched)}/{len(summaries)} matched")

            except Exception as e:
                logger.error(f"  {entry_key}: FAILED — {e}")

            time.sleep(2)

        try:
            client.files.delete(file_id)
        except Exception:
            pass

    logger.info(f"Done: {total_clusters} clusters, {total_in:,} in + {total_out:,} out tokens")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/test_extract.py -v`
Expected: all 3 tests PASS

- [ ] **Step 5: Delete old script**

```bash
rm scripts/run_extraction.py
```

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/extraction/cluster/extract.py tests/test_extract.py
git rm scripts/run_extraction.py
git commit -m "feat: clean extract.py module with report and verify"
```

---

### Task 6: Add verification tests

**Files:**
- Modify: `tests/test_extract.py`

- [ ] **Step 1: Write verification tests**

Append to `tests/test_extract.py`:

```python
def test_detect_duplicate_ids():
    """Duplicate id within analysis → warning."""
    from multiomics_kg.extraction.cluster.extract import verify_quality
    from multiomics_kg.extraction.cluster.extraction_utils import save_extraction

    import tempfile
    with tempfile.TemporaryDirectory() as td:
        paper_dir = Path(td)
        clusters = {
            "1": {"id": "dup_id", "name": "C1", "functional_description": "x" * 30,
                  "behavioral_description": "", "direction": "up"},
            "2": {"id": "dup_id", "name": "C2", "functional_description": "y" * 30,
                  "behavioral_description": "", "direction": "down"},
        }
        save_extraction(paper_dir, "test", {"paper": "Test"}, clusters)
        entries = [(paper_dir, "test", {}, {"papername": "Test"})]
        warnings = verify_quality(entries)
        assert any("duplicate id" in w for w in warnings)


def test_detect_locus_tags():
    """Locus tag in description → warning."""
    from multiomics_kg.extraction.cluster.extract import verify_quality
    from multiomics_kg.extraction.cluster.extraction_utils import save_extraction

    import tempfile
    with tempfile.TemporaryDirectory() as td:
        paper_dir = Path(td)
        clusters = {
            "1": {"id": "x", "name": "C1",
                  "functional_description": "Includes gene PMM0042 involved in transport",
                  "behavioral_description": "Upregulated", "direction": "up"},
        }
        save_extraction(paper_dir, "test", {"paper": "Test"}, clusters)
        entries = [(paper_dir, "test", {}, {"papername": "Test"})]
        warnings = verify_quality(entries)
        assert any("locus tag" in w and "PMM0042" in w for w in warnings)


def test_detect_empty_direction():
    """Non-empty description with empty direction → warning."""
    from multiomics_kg.extraction.cluster.extract import verify_quality
    from multiomics_kg.extraction.cluster.extraction_utils import save_extraction

    import tempfile
    with tempfile.TemporaryDirectory() as td:
        paper_dir = Path(td)
        clusters = {
            "1": {"id": "x", "name": "C1",
                  "functional_description": "A real description with enough text",
                  "behavioral_description": "", "direction": ""},
        }
        save_extraction(paper_dir, "test", {"paper": "Test"}, clusters)
        entries = [(paper_dir, "test", {}, {"papername": "Test"})]
        warnings = verify_quality(entries)
        assert any("empty direction" in w for w in warnings)
```

- [ ] **Step 2: Run tests**

Run: `uv run pytest tests/test_extract.py -v`
Expected: all 6 tests PASS

- [ ] **Step 3: Commit**

```bash
git add tests/test_extract.py
git commit -m "test: add verification checks for duplicate ids, locus tags, empty direction"
```

---

### Task 7: Delete old pipeline code, junk files, old tests

**Files:**
- Delete: 16 source files, 4 test files, junk directories (see list below)

- [ ] **Step 1: Delete old pipeline code**

```bash
git rm multiomics_kg/extraction/cluster/visual.py
git rm multiomics_kg/extraction/cluster/semantic.py
git rm multiomics_kg/extraction/cluster/merge.py
git rm multiomics_kg/extraction/cluster/synthesis.py
git rm multiomics_kg/extraction/cluster/validation.py
git rm multiomics_kg/extraction/cluster/pipeline.py
git rm multiomics_kg/extraction/cluster/migrate.py
git rm multiomics_kg/extraction/cluster/prompts.py
git rm multiomics_kg/extraction/cluster/run_manager.py
git rm multiomics_kg/extraction/rag.py
git rm multiomics_kg/review/cluster_review_app.py
git rm multiomics_kg/review/review_components.py
git rm multiomics_kg/review/review_data.py
```

- [ ] **Step 2: Delete old scripts**

```bash
git rm scripts/pilot_extraction.py
git rm scripts/dry_run_prompts.py
git rm scripts/rag_experiment.py
```

- [ ] **Step 3: Delete old tests**

```bash
git rm tests/test_cluster_extraction_migration.py
git rm tests/test_review_data.py
git rm tests/test_extraction_merge.py
git rm tests/test_run_manager.py
```

- [ ] **Step 4: Delete junk data directories**

```bash
rm -rf data/Prochlorococcus/papers_and_supp/tolonen\ 2006/.extraction_cache
rm -rf data/Prochlorococcus/papers_and_supp/tolonen\ 2006/pilot_comparison
# Find and remove any other .extraction_cache dirs
find data -name ".extraction_cache" -type d -exec rm -rf {} + 2>/dev/null
find data -name "pilot_comparison" -type d -exec rm -rf {} + 2>/dev/null
```

- [ ] **Step 5: Verify no broken imports**

```bash
uv run python -c "from multiomics_kg.adapters.cluster_adapter import ClusterAdapter; print('adapter OK')"
uv run python -c "from multiomics_kg.extraction.cluster.extract import main; print('extract OK')"
uv run python -c "from multiomics_kg.extraction.cluster.extraction_utils import load_extraction; print('utils OK')"
```
Expected: all print OK

- [ ] **Step 6: Grep for stale references**

```bash
grep -r "RunManager\|stage4_review\|stage2_results\|stage3_validation\|stage1_merged" multiomics_kg/ tests/
```
Expected: no output

- [ ] **Step 7: Run all tests**

Run: `uv run pytest -m "not slow and not kg" -v`
Expected: PASS (no test references deleted modules)

- [ ] **Step 8: Commit**

```bash
git add -A
git commit -m "chore: delete old pipeline code, review UI, RAG, old tests, junk dirs"
```

---

### Task 8: Generate initial quality report + verify

**Files:**
- Creates: `data/cluster_extraction_report.md`

- [ ] **Step 1: Generate report with verification**

```bash
uv run python -m multiomics_kg.extraction.cluster.extract --report --verify
```
Expected: `data/cluster_extraction_report.md` created, warnings printed

- [ ] **Step 2: Review warnings**

Read the warnings section. Expected issues:
- Duplicate ids (e.g., multiple "med4_down_not_described" in Zinser diel clusters)
- Possibly some locus tags that slipped through

- [ ] **Step 3: Commit report**

```bash
git add data/cluster_extraction_report.md
git commit -m "docs: initial cluster extraction quality report"
```

---

### Task 9: Delete migration script, final cleanup

**Files:**
- Delete: `scripts/migrate_extractions.py`

- [ ] **Step 1: Delete migration script** (one-time, no longer needed)

```bash
git rm scripts/migrate_extractions.py
```

- [ ] **Step 2: Final verification**

```bash
# All tests pass
uv run pytest -m "not slow and not kg" -v

# Module is runnable
uv run python -m multiomics_kg.extraction.cluster.extract --dry-run

# Report regenerates identically
uv run python -m multiomics_kg.extraction.cluster.extract --report
git diff data/cluster_extraction_report.md  # should be empty
```

- [ ] **Step 3: Commit**

```bash
git rm scripts/migrate_extractions.py
git commit -m "chore: remove one-time migration script, final cleanup"
```

---

## Code Review Checklist

After all tasks complete:
- [ ] No imports of deleted modules anywhere: `grep -r "RunManager\|stage4_review\|stage2_results\|stage3_validation\|stage1_merged" multiomics_kg/ tests/`
- [ ] All 15 extraction JSONs in `cluster_extractions/` with `clusters` key
- [ ] Report generates without errors
- [ ] Verification checks pass (or warnings documented)
- [ ] `uv run pytest -m "not slow and not kg"` passes
- [ ] KG build produces same 82 populated + 33 empty clusters (verify after Docker rebuild)
