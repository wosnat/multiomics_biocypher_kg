# Cluster Extraction Review System Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor the cluster extraction pipeline to use versioned run directories with pre-processing caches, then build a Streamlit review UI for human-in-the-loop quality improvement.

**Architecture:** The pipeline writes per-stage JSON files into timestamped run directories under `{paper_dir}/.extraction_cache/{entry_key}/runs/`. A shared pre-processing cache avoids redundant PDF/table parsing. A Streamlit app reads these files and provides merge-centric cluster review with issue classification. The cluster adapter is updated to respect review status before using extraction data in the KG.

**Tech Stack:** Python, Streamlit, existing OpenAI-based extraction pipeline, JSON file storage, hashlib for input hashing.

**Spec:** `docs/superpowers/specs/2026-04-04-cluster-extraction-review-system-design.md`

---

## File Structure

| File | Responsibility |
|---|---|
| `multiomics_kg/extraction/cluster/run_manager.py` (create) | Run directory creation, symlink management, stage file I/O, input hashing, review copy-forward |
| `multiomics_kg/extraction/cluster/pipeline.py` (modify) | Refactor to use run_manager instead of monolithic JSON output |
| `multiomics_kg/extraction/cluster/migrate.py` (create) | One-time migration of existing `cluster_extraction_*.json` into new structure |
| `multiomics_kg/adapters/cluster_adapter.py` (modify) | Read from new directory structure + respect review status |
| `multiomics_kg/review/__init__.py` (create) | Package init |
| `multiomics_kg/review/cluster_review_app.py` (create) | Streamlit app main entry point, layout, navigation |
| `multiomics_kg/review/review_data.py` (create) | Data loading: scan papers, load runs, load review state, export issue report |
| `multiomics_kg/review/review_components.py` (create) | Streamlit UI components: merge view, review controls, source viewer, diff view |
| `tests/test_run_manager.py` (create) | Tests for run directory management |
| `tests/test_cluster_extraction_migration.py` (create) | Tests for migration of old JSONs |
| `tests/test_review_data.py` (create) | Tests for review data loading and export |

---

## Task 1: Run Manager — Directory Structure and Stage I/O

**Files:**
- Create: `multiomics_kg/extraction/cluster/run_manager.py`
- Test: `tests/test_run_manager.py`

- [ ] **Step 1: Write tests for run directory creation**

```python
# tests/test_run_manager.py
import json
import time
from pathlib import Path
from multiomics_kg.extraction.cluster.run_manager import RunManager


def test_create_run_returns_timestamped_dir(tmp_path):
    """New run creates a timestamped directory under runs/."""
    rm = RunManager(tmp_path, "med4_kmeans_nstarvation")
    run_dir = rm.create_run()
    assert run_dir.parent.name == "runs"
    assert run_dir.parent.parent.name == "med4_kmeans_nstarvation"
    # Timestamp format: YYYY-MM-DDTHH-MM-SS
    assert len(run_dir.name) == 19
    assert run_dir.exists()


def test_create_run_updates_current_symlink(tmp_path):
    """create_run updates the 'current' symlink to point at the new run."""
    rm = RunManager(tmp_path, "med4_kmeans_nstarvation")
    run_dir = rm.create_run()
    current = tmp_path / "med4_kmeans_nstarvation" / "current"
    assert current.is_symlink()
    assert current.resolve() == run_dir.resolve()


def test_create_second_run_preserves_first(tmp_path):
    """Creating a second run does not delete the first."""
    rm = RunManager(tmp_path, "med4_kmeans_nstarvation")
    run1 = rm.create_run()
    time.sleep(1.1)  # ensure different timestamp
    run2 = rm.create_run()
    assert run1.exists()
    assert run2.exists()
    assert run1 != run2
    current = tmp_path / "med4_kmeans_nstarvation" / "current"
    assert current.resolve() == run2.resolve()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_run_manager.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'multiomics_kg.extraction.cluster.run_manager'`

- [ ] **Step 3: Implement RunManager — directory creation**

```python
# multiomics_kg/extraction/cluster/run_manager.py
"""Manages versioned extraction run directories and stage file I/O."""

import hashlib
import json
import os
from datetime import datetime
from pathlib import Path
from typing import Optional


class RunManager:
    """Manages the .extraction_cache/{entry_key}/ directory for one paperconfig entry.

    Directory layout:
        {cache_root}/{entry_key}/
            runs/
                2026-04-04T10-30-00/
                    stage1_merged.json
                    stage2_results.json
                    stage3_validation.json
                    stage4_review.json
                ...
            current -> runs/2026-04-04T10-30-00  (symlink)

        {cache_root}/shared/
            pdf_text.json
            pages/
            embeddings.npz
            rag_chunks.json
            tables/
    """

    def __init__(self, cache_root: Path, entry_key: str):
        self.cache_root = Path(cache_root)
        self.entry_key = entry_key
        self.entry_dir = self.cache_root / entry_key
        self.runs_dir = self.entry_dir / "runs"
        self.current_link = self.entry_dir / "current"
        self.shared_dir = self.cache_root / "shared"

    def create_run(self) -> Path:
        """Create a new timestamped run directory and update 'current' symlink."""
        timestamp = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
        run_dir = self.runs_dir / timestamp
        run_dir.mkdir(parents=True, exist_ok=True)
        # Update symlink atomically
        tmp_link = self.current_link.with_suffix(".tmp")
        if tmp_link.exists() or tmp_link.is_symlink():
            tmp_link.unlink()
        tmp_link.symlink_to(run_dir)
        os.replace(tmp_link, self.current_link)
        return run_dir

    def get_current_run(self) -> Optional[Path]:
        """Return the current run directory, or None if no runs exist."""
        if self.current_link.is_symlink():
            target = self.current_link.resolve()
            if target.exists():
                return target
        return None

    def list_runs(self) -> list[Path]:
        """Return all run directories, oldest first."""
        if not self.runs_dir.exists():
            return []
        return sorted(
            [d for d in self.runs_dir.iterdir() if d.is_dir()],
            key=lambda p: p.name,
        )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_run_manager.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Write tests for stage file I/O**

Add to `tests/test_run_manager.py`:

```python
def test_write_and_read_stage(tmp_path):
    """write_stage persists JSON, read_stage retrieves it."""
    rm = RunManager(tmp_path, "entry1")
    run_dir = rm.create_run()
    data = {"1": {"enrichment_category": [{"value": "Transport", "source": "table"}]}}
    rm.write_stage(run_dir, 1, data)
    assert (run_dir / "stage1_merged.json").exists()
    loaded = rm.read_stage(run_dir, 1)
    assert loaded == data


def test_read_stage_missing_returns_empty_dict(tmp_path):
    """read_stage returns {} when stage file does not exist."""
    rm = RunManager(tmp_path, "entry1")
    run_dir = rm.create_run()
    assert rm.read_stage(run_dir, 2) == {}


def test_read_current_stage(tmp_path):
    """read_current_stage reads from the current symlinked run."""
    rm = RunManager(tmp_path, "entry1")
    run_dir = rm.create_run()
    data = {"1": {"id": "up_transport"}}
    rm.write_stage(run_dir, 2, data)
    assert rm.read_current_stage(2) == data


STAGE_FILENAMES = {1: "stage1_merged.json", 2: "stage2_results.json",
                   3: "stage3_validation.json", 4: "stage4_review.json"}


def test_stage_filenames(tmp_path):
    """Each stage number maps to the correct filename."""
    rm = RunManager(tmp_path, "entry1")
    run_dir = rm.create_run()
    for stage_num, expected_name in STAGE_FILENAMES.items():
        rm.write_stage(run_dir, stage_num, {"test": stage_num})
        assert (run_dir / expected_name).exists()
```

- [ ] **Step 6: Implement stage file I/O methods**

Add to `RunManager` class in `run_manager.py`:

```python
    STAGE_FILES = {
        1: "stage1_merged.json",
        2: "stage2_results.json",
        3: "stage3_validation.json",
        4: "stage4_review.json",
    }

    def write_stage(self, run_dir: Path, stage: int, data: dict) -> Path:
        """Write stage data to the appropriate JSON file in run_dir."""
        filename = self.STAGE_FILES[stage]
        path = run_dir / filename
        with open(path, "w") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        return path

    def read_stage(self, run_dir: Path, stage: int) -> dict:
        """Read stage data from run_dir. Returns {} if file missing."""
        filename = self.STAGE_FILES[stage]
        path = run_dir / filename
        if not path.exists():
            return {}
        with open(path) as f:
            return json.load(f)

    def read_current_stage(self, stage: int) -> dict:
        """Read stage data from the current run. Returns {} if no current run."""
        run_dir = self.get_current_run()
        if run_dir is None:
            return {}
        return self.read_stage(run_dir, stage)
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `pytest tests/test_run_manager.py -v`
Expected: PASS (7 tests)

- [ ] **Step 8: Commit**

```bash
git add multiomics_kg/extraction/cluster/run_manager.py tests/test_run_manager.py
git commit -m "feat: add RunManager for versioned extraction run directories"
```

---

## Task 2: Run Manager — Input Hashing and Review Copy-Forward

**Files:**
- Modify: `multiomics_kg/extraction/cluster/run_manager.py`
- Test: `tests/test_run_manager.py`

- [ ] **Step 1: Write tests for input hashing**

Add to `tests/test_run_manager.py`:

```python
def test_compute_input_hash_deterministic(tmp_path):
    """Same stage1 data produces same hash."""
    rm = RunManager(tmp_path, "entry1")
    stage1 = {
        "1": {"enrichment_category": [{"value": "Transport", "source": "table", "confidence": "very_high"}]},
        "2": {"enrichment_category": [{"value": "Stress", "source": "visual", "confidence": "high"}]},
    }
    h1 = rm.compute_input_hash("1", stage1)
    h2 = rm.compute_input_hash("1", stage1)
    assert h1 == h2


def test_compute_input_hash_differs_across_clusters(tmp_path):
    """Different clusters produce different hashes."""
    rm = RunManager(tmp_path, "entry1")
    stage1 = {
        "1": {"enrichment_category": [{"value": "Transport"}]},
        "2": {"enrichment_category": [{"value": "Stress"}]},
    }
    assert rm.compute_input_hash("1", stage1) != rm.compute_input_hash("2", stage1)


def test_compute_input_hash_changes_on_data_change(tmp_path):
    """Hash changes when stage1 data for a cluster changes."""
    rm = RunManager(tmp_path, "entry1")
    stage1_v1 = {"1": {"enrichment_category": [{"value": "Transport"}]}}
    stage1_v2 = {"1": {"enrichment_category": [{"value": "Nitrogen transport"}]}}
    assert rm.compute_input_hash("1", stage1_v1) != rm.compute_input_hash("1", stage1_v2)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_run_manager.py::test_compute_input_hash_deterministic -v`
Expected: FAIL — `AttributeError: 'RunManager' object has no attribute 'compute_input_hash'`

- [ ] **Step 3: Implement input hashing**

Add to `RunManager` class:

```python
    @staticmethod
    def compute_input_hash(cluster_key: str, stage1_merged: dict) -> str:
        """Compute a deterministic hash of stage1 inputs for one cluster.

        Used to detect whether re-running synthesis would use different inputs.
        """
        cluster_data = stage1_merged.get(str(cluster_key), {})
        serialized = json.dumps(cluster_data, sort_keys=True, ensure_ascii=False)
        return hashlib.sha256(serialized.encode()).hexdigest()[:16]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_run_manager.py -k "input_hash" -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Write tests for review copy-forward**

Add to `tests/test_run_manager.py`:

```python
def test_copy_forward_reviews_unchanged_input(tmp_path):
    """Reviews carry forward with same status when input hash matches."""
    rm = RunManager(tmp_path, "entry1")
    stage1 = {"1": {"enrichment_category": [{"value": "Transport"}]}}

    # First run with review
    run1 = rm.create_run()
    rm.write_stage(run1, 1, stage1)
    review = {
        "1": {
            "status": "approve",
            "issues": [],
            "failing_stages": [],
            "notes": "",
            "edited_fields": {},
            "reviewed_at": "2026-04-04T12:00:00",
            "reviewed_in_run": run1.name,
            "input_hash": rm.compute_input_hash("1", stage1),
        }
    }
    rm.write_stage(run1, 4, review)

    # Second run with identical stage1
    time.sleep(1.1)
    run2 = rm.create_run()
    rm.write_stage(run2, 1, stage1)
    rm.copy_forward_reviews(run1, run2, stage1)

    new_review = rm.read_stage(run2, 4)
    assert new_review["1"]["status"] == "approve"
    assert new_review["1"]["reviewed_in_run"] == run1.name  # preserved


def test_copy_forward_reviews_marks_stale_on_change(tmp_path):
    """Reviews marked stale when input hash changes."""
    rm = RunManager(tmp_path, "entry1")
    stage1_v1 = {"1": {"enrichment_category": [{"value": "Transport"}]}}
    stage1_v2 = {"1": {"enrichment_category": [{"value": "Nitrogen uptake"}]}}

    run1 = rm.create_run()
    rm.write_stage(run1, 1, stage1_v1)
    review = {
        "1": {
            "status": "approve",
            "issues": [],
            "failing_stages": [],
            "notes": "",
            "edited_fields": {},
            "reviewed_at": "2026-04-04T12:00:00",
            "reviewed_in_run": run1.name,
            "input_hash": rm.compute_input_hash("1", stage1_v1),
        }
    }
    rm.write_stage(run1, 4, review)

    time.sleep(1.1)
    run2 = rm.create_run()
    rm.write_stage(run2, 1, stage1_v2)
    rm.copy_forward_reviews(run1, run2, stage1_v2)

    new_review = rm.read_stage(run2, 4)
    assert new_review["1"]["status"] == "stale"


def test_copy_forward_reviews_skips_clusters_not_in_new_run(tmp_path):
    """Clusters absent from new stage1 are not copied."""
    rm = RunManager(tmp_path, "entry1")
    stage1_old = {"1": {"x": []}, "2": {"x": []}}
    stage1_new = {"1": {"x": []}}  # cluster 2 removed

    run1 = rm.create_run()
    rm.write_stage(run1, 1, stage1_old)
    review = {
        "1": {"status": "approve", "input_hash": rm.compute_input_hash("1", stage1_old),
              "issues": [], "failing_stages": [], "notes": "", "edited_fields": {},
              "reviewed_at": "2026-04-04T12:00:00", "reviewed_in_run": run1.name},
        "2": {"status": "reject", "input_hash": rm.compute_input_hash("2", stage1_old),
              "issues": [], "failing_stages": [], "notes": "", "edited_fields": {},
              "reviewed_at": "2026-04-04T12:00:00", "reviewed_in_run": run1.name},
    }
    rm.write_stage(run1, 4, review)

    time.sleep(1.1)
    run2 = rm.create_run()
    rm.write_stage(run2, 1, stage1_new)
    rm.copy_forward_reviews(run1, run2, stage1_new)

    new_review = rm.read_stage(run2, 4)
    assert "1" in new_review
    assert "2" not in new_review
```

- [ ] **Step 6: Run tests to verify they fail**

Run: `pytest tests/test_run_manager.py -k "copy_forward" -v`
Expected: FAIL — `AttributeError: 'RunManager' object has no attribute 'copy_forward_reviews'`

- [ ] **Step 7: Implement review copy-forward**

Add to `RunManager` class:

```python
    def copy_forward_reviews(
        self, prev_run: Path, new_run: Path, new_stage1: dict
    ) -> None:
        """Copy review decisions from prev_run to new_run.

        - If input hash matches: review carries forward as-is
        - If input hash differs: review status set to 'stale'
        - Clusters not in new_stage1 are dropped
        """
        prev_review = self.read_stage(prev_run, 4)
        if not prev_review:
            return

        new_review = {}
        for cluster_key, review_entry in prev_review.items():
            if cluster_key not in new_stage1:
                continue
            new_hash = self.compute_input_hash(cluster_key, new_stage1)
            old_hash = review_entry.get("input_hash", "")
            entry_copy = dict(review_entry)
            if new_hash != old_hash:
                entry_copy["status"] = "stale"
            entry_copy["input_hash"] = new_hash
            new_review[cluster_key] = entry_copy

        if new_review:
            self.write_stage(new_run, 4, new_review)
```

- [ ] **Step 8: Run all tests**

Run: `pytest tests/test_run_manager.py -v`
Expected: PASS (13 tests)

- [ ] **Step 9: Commit**

```bash
git add multiomics_kg/extraction/cluster/run_manager.py tests/test_run_manager.py
git commit -m "feat: add input hashing and review copy-forward to RunManager"
```

---

## Task 3: Migration — Convert Existing Extraction JSONs

**Files:**
- Create: `multiomics_kg/extraction/cluster/migrate.py`
- Test: `tests/test_cluster_extraction_migration.py`

- [ ] **Step 1: Write migration tests**

```python
# tests/test_cluster_extraction_migration.py
import json
from pathlib import Path
from multiomics_kg.extraction.cluster.migrate import migrate_extraction_json
from multiomics_kg.extraction.cluster.run_manager import RunManager


def _make_legacy_json(paper_dir: Path, entry_key: str) -> Path:
    """Create a minimal legacy-format extraction JSON."""
    data = {
        "metadata": {
            "paper": "Test Paper",
            "doi": "10.1234/test",
            "organism": "Prochlorococcus MED4",
            "extracted_at": "2026-04-01T10:32:54",
            "table_key": entry_key,
        },
        "stage1_merged": {
            "1": {"enrichment_category": [{"value": "Transport", "source": "table", "confidence": "very_high"}],
                  "gene_count": 5},
            "2": {"enrichment_category": [{"value": "Stress", "source": "visual", "confidence": "high"}],
                  "gene_count": 3},
        },
        "stage2_results": {
            "1": {"id": "up_transport", "name": "Cluster 1", "functional_description": "Transport",
                  "behavioral_description": "Upregulated", "peak_time_hours": None, "period_hours": None},
            "2": {"id": "stress_resp", "name": "Cluster 2", "functional_description": "Stress",
                  "behavioral_description": "Late response", "peak_time_hours": 12, "period_hours": None},
        },
        "stage3_validation": {
            "1": {"verdict": "pass", "explanation": "Supported by text"},
            "2": {"verdict": "warn", "explanation": "Inferred"},
        },
    }
    path = paper_dir / f"cluster_extraction_{entry_key}.json"
    with open(path, "w") as f:
        json.dump(data, f)
    return path


def test_migrate_creates_run_directory(tmp_path):
    """Migration creates a run with stage1-3 files."""
    paper_dir = tmp_path / "paper"
    paper_dir.mkdir()
    _make_legacy_json(paper_dir, "my_entry")

    cache_dir = paper_dir / ".extraction_cache"
    migrate_extraction_json(paper_dir, "my_entry", cache_dir)

    rm = RunManager(cache_dir, "my_entry")
    run_dir = rm.get_current_run()
    assert run_dir is not None
    assert (run_dir / "stage1_merged.json").exists()
    assert (run_dir / "stage2_results.json").exists()
    assert (run_dir / "stage3_validation.json").exists()


def test_migrate_preserves_data(tmp_path):
    """Migrated stage files contain the correct data."""
    paper_dir = tmp_path / "paper"
    paper_dir.mkdir()
    _make_legacy_json(paper_dir, "my_entry")

    cache_dir = paper_dir / ".extraction_cache"
    migrate_extraction_json(paper_dir, "my_entry", cache_dir)

    rm = RunManager(cache_dir, "my_entry")
    run_dir = rm.get_current_run()
    stage1 = rm.read_stage(run_dir, 1)
    assert "1" in stage1
    assert stage1["1"]["enrichment_category"][0]["value"] == "Transport"

    stage2 = rm.read_stage(run_dir, 2)
    assert stage2["1"]["id"] == "up_transport"

    stage3 = rm.read_stage(run_dir, 3)
    assert stage3["1"]["verdict"] == "pass"


def test_migrate_skips_if_no_legacy_json(tmp_path):
    """Migration is a no-op when no legacy file exists."""
    paper_dir = tmp_path / "paper"
    paper_dir.mkdir()
    cache_dir = paper_dir / ".extraction_cache"
    result = migrate_extraction_json(paper_dir, "missing_entry", cache_dir)
    assert result is None


def test_migrate_stores_metadata(tmp_path):
    """Migration preserves metadata in a separate file in the run directory."""
    paper_dir = tmp_path / "paper"
    paper_dir.mkdir()
    _make_legacy_json(paper_dir, "my_entry")

    cache_dir = paper_dir / ".extraction_cache"
    migrate_extraction_json(paper_dir, "my_entry", cache_dir)

    rm = RunManager(cache_dir, "my_entry")
    run_dir = rm.get_current_run()
    metadata_path = run_dir / "metadata.json"
    assert metadata_path.exists()
    with open(metadata_path) as f:
        meta = json.load(f)
    assert meta["paper"] == "Test Paper"
    assert meta["migrated_from"] is not None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_cluster_extraction_migration.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'multiomics_kg.extraction.cluster.migrate'`

- [ ] **Step 3: Implement migration**

```python
# multiomics_kg/extraction/cluster/migrate.py
"""One-time migration of legacy cluster_extraction_*.json into versioned run structure."""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

from multiomics_kg.extraction.cluster.run_manager import RunManager

logger = logging.getLogger(__name__)


def migrate_extraction_json(
    paper_dir: Path, entry_key: str, cache_dir: Path
) -> Optional[Path]:
    """Migrate a legacy cluster_extraction_{entry_key}.json into the new run structure.

    Args:
        paper_dir: Directory containing the legacy JSON file.
        entry_key: The gene_clusters entry key (e.g., 'med4_kmeans_nstarvation').
        cache_dir: The .extraction_cache directory to write into.

    Returns:
        Path to the new run directory, or None if no legacy file found.
    """
    legacy_path = paper_dir / f"cluster_extraction_{entry_key}.json"
    if not legacy_path.exists():
        return None

    with open(legacy_path) as f:
        data = json.load(f)

    rm = RunManager(cache_dir, entry_key)
    run_dir = rm.create_run()

    # Split monolithic JSON into per-stage files
    if "stage1_merged" in data:
        rm.write_stage(run_dir, 1, data["stage1_merged"])
    if "stage2_results" in data:
        rm.write_stage(run_dir, 2, data["stage2_results"])
    if "stage3_validation" in data:
        rm.write_stage(run_dir, 3, data["stage3_validation"])

    # Save metadata separately
    metadata = dict(data.get("metadata", {}))
    metadata["migrated_from"] = str(legacy_path)
    metadata["migrated_at"] = datetime.now().isoformat(timespec="seconds")
    with open(run_dir / "metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    logger.info("Migrated %s → %s", legacy_path, run_dir)
    return run_dir


def migrate_all_tolonen(project_root: Path) -> list[Path]:
    """Migrate both Tolonen 2006 extraction JSONs."""
    paper_dir = project_root / "data" / "Prochlorococcus" / "papers_and_supp" / "tolonen 2006"
    cache_dir = paper_dir / ".extraction_cache"
    results = []
    for entry_key in ["med4_kmeans_nstarvation", "mit9313_kmeans_nstarvation"]:
        run_dir = migrate_extraction_json(paper_dir, entry_key, cache_dir)
        if run_dir:
            results.append(run_dir)
    return results
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_cluster_extraction_migration.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/extraction/cluster/migrate.py tests/test_cluster_extraction_migration.py
git commit -m "feat: add migration for legacy extraction JSONs to versioned runs"
```

---

## Task 4: Refactor Pipeline to Use Run Manager

**Files:**
- Modify: `multiomics_kg/extraction/cluster/pipeline.py`
- Test: `tests/test_run_manager.py` (integration test)

- [ ] **Step 1: Write integration test for pipeline output structure**

Add to `tests/test_run_manager.py`:

```python
from multiomics_kg.extraction.cluster.run_manager import RunManager


def test_pipeline_output_structure_matches_spec(tmp_path):
    """Verify the directory structure created by RunManager matches the spec.

    This is a structural test — it doesn't run the actual pipeline,
    just validates the file layout that the refactored pipeline should produce.
    """
    cache_dir = tmp_path / ".extraction_cache"
    rm = RunManager(cache_dir, "med4_kmeans")
    run_dir = rm.create_run()

    # Write all 4 stages
    rm.write_stage(run_dir, 1, {"1": {"data": "stage1"}})
    rm.write_stage(run_dir, 2, {"1": {"id": "test"}})
    rm.write_stage(run_dir, 3, {"1": {"verdict": "pass"}})
    rm.write_stage(run_dir, 4, {"1": {"status": "approve"}})

    # Write metadata
    import json
    with open(run_dir / "metadata.json", "w") as f:
        json.dump({"paper": "Test"}, f)

    # Verify structure
    assert (cache_dir / "med4_kmeans" / "runs").is_dir()
    assert (cache_dir / "med4_kmeans" / "current").is_symlink()
    assert (run_dir / "stage1_merged.json").exists()
    assert (run_dir / "stage2_results.json").exists()
    assert (run_dir / "stage3_validation.json").exists()
    assert (run_dir / "stage4_review.json").exists()
    assert (run_dir / "metadata.json").exists()

    # Verify shared cache directory creation
    shared = cache_dir / "shared"
    shared.mkdir(parents=True, exist_ok=True)
    assert shared.is_dir()
```

- [ ] **Step 2: Run test to verify it passes** (structural test, should pass immediately)

Run: `pytest tests/test_run_manager.py::test_pipeline_output_structure_matches_spec -v`
Expected: PASS

- [ ] **Step 3: Refactor pipeline.py — update `run_pipeline` to use RunManager**

Read the current `pipeline.py` first. Then modify `run_pipeline()` (lines 117-239) to:

1. Create a `RunManager` and `create_run()` instead of building `out_path` directly
2. Write `stage1_merged.json`, `stage2_results.json`, `stage3_validation.json` via `rm.write_stage()`
3. Write `metadata.json` to the run directory
4. Copy forward reviews from previous run if one exists
5. Keep the markdown report generation (write to run directory)
6. Add `--force-llm` and `--force-all` CLI flags to `main()`

Key changes to `run_pipeline()`:

```python
def run_pipeline(
    paperconfig_path: Path,
    table_key: Optional[str] = None,
    stage: Optional[int] = None,
    from_stage: Optional[int] = None,
    path_filter: Optional[str] = None,
    force_llm: bool = False,
    force_all: bool = False,
) -> list[Path]:
    """Run extraction pipeline, writing results to versioned run directories."""
    with open(paperconfig_path) as f:
        config = yaml.safe_load(f)

    paper_dir = paperconfig_path.parent
    cache_dir = paper_dir / ".extraction_cache"
    paper_name = config["publication"]["papername"]
    project_root = _find_project_root(paper_dir)
    doi = _lookup_doi(paper_name, project_root)

    tables = {
        k: v for k, v in config["publication"]["supplementary_materials"].items()
        if v.get("type") == "gene_clusters"
    }
    if table_key:
        tables = {table_key: tables[table_key]}

    output_dirs = []
    for tkey, tconf in tables.items():
        rm = RunManager(cache_dir, tkey)
        prev_run = rm.get_current_run()
        run_dir = rm.create_run()

        organism = tconf["organism"]
        treatment = tconf.get("treatment", "")
        cluster_method = tconf.get("cluster_method", "unknown")

        # Stage 1
        effective_stage = stage or from_stage or 1
        if effective_stage <= 1:
            merged = run_stage1(paperconfig_path, tkey, tconf, path_filter)
        elif prev_run:
            merged = rm.read_stage(prev_run, 1)
        else:
            merged = {}
        rm.write_stage(run_dir, 1, merged)

        cluster_keys = sorted(merged.keys())

        # Stage 2
        if effective_stage <= 2 and merged:
            stage2 = run_synthesis(merged, paper_name, organism, treatment, cluster_method)
        elif prev_run:
            stage2 = rm.read_stage(prev_run, 2)
        else:
            stage2 = {}
        rm.write_stage(run_dir, 2, stage2)

        # Stage 3
        if effective_stage <= 3 and stage2:
            main_pdf = _find_main_pdf(paper_dir, config)
            csv_path = Path(tconf["filename"])
            if main_pdf and csv_path.exists():
                stage3 = run_validation(main_pdf, paper_dir, csv_path, stage2, cluster_keys)
            else:
                stage3 = {}
        elif prev_run:
            stage3 = rm.read_stage(prev_run, 3)
        else:
            stage3 = {}
        rm.write_stage(run_dir, 3, stage3)

        # Copy forward reviews from previous run
        if prev_run:
            rm.copy_forward_reviews(prev_run, run_dir, merged)

        # Write metadata
        metadata = {
            "paper": paper_name,
            "doi": doi,
            "organism": organism,
            "extracted_at": datetime.now().isoformat(timespec="seconds"),
            "table_key": tkey,
        }
        with open(run_dir / "metadata.json", "w") as f:
            json.dump(metadata, f, indent=2)

        # Generate report
        report = _generate_report(
            {"stage2_results": stage2, "stage3_validation": stage3},
            cluster_keys,
        )
        with open(run_dir / "report.md", "w") as f:
            f.write(report)

        output_dirs.append(run_dir)

    return output_dirs
```

Update `main()` to add CLI flags:

```python
def main():
    parser = argparse.ArgumentParser(description="Cluster extraction pipeline")
    parser.add_argument("paperconfig", type=Path, help="Path to paperconfig.yaml")
    parser.add_argument("--table-key", help="Specific gene_clusters entry key")
    parser.add_argument("--stage", type=int, help="Run only this stage")
    parser.add_argument("--from-stage", type=int, help="Run from this stage onwards")
    parser.add_argument("--path-filter", help="Run only this extraction path")
    parser.add_argument("--force-llm", action="store_true", help="Re-run all LLM calls")
    parser.add_argument("--force-all", action="store_true", help="Clear all caches and re-extract")
    args = parser.parse_args()

    results = run_pipeline(
        args.paperconfig,
        table_key=args.table_key,
        stage=args.stage,
        from_stage=args.from_stage,
        path_filter=args.path_filter,
        force_llm=args.force_llm,
        force_all=args.force_all,
    )
    for r in results:
        print(f"Output: {r}")
```

- [ ] **Step 4: Run existing pipeline tests to check for regressions**

Run: `pytest tests/test_extraction_merge.py tests/test_extraction_table.py -v`
Expected: PASS (merge and table tests don't depend on pipeline output format)

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/extraction/cluster/pipeline.py tests/test_run_manager.py
git commit -m "refactor: pipeline.py uses RunManager for versioned run output"
```

---

## Task 5: Update Cluster Adapter to Read New Structure

**Files:**
- Modify: `multiomics_kg/adapters/cluster_adapter.py`
- Modify: `tests/test_cluster_adapter.py`

- [ ] **Step 1: Write test for adapter reading new structure**

Add to `tests/test_cluster_adapter.py`:

```python
def test_load_extraction_from_cache_dir(tmp_path):
    """Adapter reads extraction data from .extraction_cache/{entry}/current/ structure."""
    from multiomics_kg.extraction.cluster.run_manager import RunManager

    # Set up cache directory with a run
    cache_dir = tmp_path / ".extraction_cache"
    rm = RunManager(cache_dir, "test_entry")
    run_dir = rm.create_run()
    rm.write_stage(run_dir, 2, {
        "1": {"id": "up_transport", "name": "Cluster 1",
              "functional_description": "Transport genes",
              "behavioral_description": "Rapid upregulation",
              "peak_time_hours": None, "period_hours": None},
    })
    rm.write_stage(run_dir, 3, {
        "1": {"verdict": "pass", "explanation": "OK"},
    })
    rm.write_stage(run_dir, 4, {
        "1": {"status": "approve", "input_hash": "abc123",
              "issues": [], "failing_stages": [], "notes": "",
              "edited_fields": {}, "reviewed_at": "2026-04-04T12:00:00",
              "reviewed_in_run": run_dir.name},
    })

    from multiomics_kg.adapters.cluster_adapter import _load_extraction_json
    result = _load_extraction_json(tmp_path, "test_entry")
    assert result.get("stage2_results", {}).get("1", {}).get("id") == "up_transport"


def test_adapter_respects_review_status(tmp_path):
    """Adapter skips cluster descriptions when review status is not approve/edit."""
    from multiomics_kg.extraction.cluster.run_manager import RunManager
    from multiomics_kg.adapters.cluster_adapter import _get_extraction_cluster_data

    cache_dir = tmp_path / ".extraction_cache"
    rm = RunManager(cache_dir, "test_entry")
    run_dir = rm.create_run()
    rm.write_stage(run_dir, 2, {
        "1": {"id": "up_transport", "functional_description": "Transport"},
    })
    rm.write_stage(run_dir, 3, {"1": {"verdict": "pass"}})
    rm.write_stage(run_dir, 4, {"1": {"status": "reject"}})

    extraction = {
        "stage2_results": rm.read_stage(run_dir, 2),
        "stage3_validation": rm.read_stage(run_dir, 3),
        "stage4_review": rm.read_stage(run_dir, 4),
    }
    result = _get_extraction_cluster_data(extraction, "1")
    assert result == {}  # rejected — no data used


def test_adapter_uses_edited_fields(tmp_path):
    """When review status is 'edit', adapter uses edited_fields over stage2."""
    from multiomics_kg.extraction.cluster.run_manager import RunManager
    from multiomics_kg.adapters.cluster_adapter import _get_extraction_cluster_data

    cache_dir = tmp_path / ".extraction_cache"
    rm = RunManager(cache_dir, "test_entry")
    run_dir = rm.create_run()
    rm.write_stage(run_dir, 2, {
        "1": {"id": "up_transport", "functional_description": "LLM version"},
    })
    rm.write_stage(run_dir, 3, {"1": {"verdict": "pass"}})
    rm.write_stage(run_dir, 4, {
        "1": {"status": "edit", "input_hash": "abc",
              "edited_fields": {"functional_description": "Human-corrected version"}},
    })

    extraction = {
        "stage2_results": rm.read_stage(run_dir, 2),
        "stage3_validation": rm.read_stage(run_dir, 3),
        "stage4_review": rm.read_stage(run_dir, 4),
    }
    result = _get_extraction_cluster_data(extraction, "1")
    assert result["functional_description"] == "Human-corrected version"
```

- [ ] **Step 2: Run new tests to verify they fail**

Run: `pytest tests/test_cluster_adapter.py::test_load_extraction_from_cache_dir -v`
Expected: FAIL (adapter still reads legacy path)

- [ ] **Step 3: Update `_load_extraction_json` to try new structure first, fall back to legacy**

Modify `multiomics_kg/adapters/cluster_adapter.py` lines 62-73:

```python
def _load_extraction_json(paperconfig_dir: Path, entry_key: str) -> dict:
    """Load extraction data for a gene_clusters entry.

    Tries new versioned structure first (.extraction_cache/{entry}/current/),
    falls back to legacy monolithic JSON.
    """
    # New structure: .extraction_cache/{entry_key}/current/
    cache_dir = paperconfig_dir / ".extraction_cache"
    entry_dir = cache_dir / entry_key
    current_link = entry_dir / "current"
    if current_link.is_symlink() or current_link.is_dir():
        run_dir = current_link.resolve() if current_link.is_symlink() else current_link
        if run_dir.exists():
            result = {}
            for stage, filename in RunManager.STAGE_FILES.items():
                stage_path = run_dir / filename
                if stage_path.exists():
                    with open(stage_path) as f:
                        key = {1: "stage1_merged", 2: "stage2_results",
                               3: "stage3_validation", 4: "stage4_review"}[stage]
                        result[key] = json.load(f)
            return result

    # Legacy fallback
    json_path = paperconfig_dir / f"cluster_extraction_{entry_key}.json"
    if not json_path.exists():
        return {}
    try:
        with open(json_path) as f:
            return json.load(f)
    except Exception:
        logger.warning("Failed to load extraction JSON: %s", json_path)
        return {}
```

Add import at top of file:

```python
from multiomics_kg.extraction.cluster.run_manager import RunManager
```

- [ ] **Step 4: Update `_get_extraction_cluster_data` to respect review status**

Modify `multiomics_kg/adapters/cluster_adapter.py` lines 76-84:

```python
def _get_extraction_cluster_data(extraction: dict, cluster_key: str) -> dict:
    """Get per-cluster data from extraction, respecting validation and review status."""
    stage2 = extraction.get("stage2_results", {})
    stage3 = extraction.get("stage3_validation", {})
    stage4 = extraction.get("stage4_review", {})

    cluster_data = stage2.get(str(cluster_key), {})
    verdict = stage3.get(str(cluster_key), {}).get("verdict", "")
    review = stage4.get(str(cluster_key), {})
    review_status = review.get("status", "")

    # If review exists, it takes precedence over verdict
    if review_status:
        if review_status in ("approve", "edit"):
            # Apply human edits over LLM output
            result = dict(cluster_data)
            for field, value in review.get("edited_fields", {}).items():
                result[field] = value
            return result
        else:
            # reject, flag-issue, stale → skip
            return {}

    # No review — fall back to verdict check (legacy behavior)
    if verdict != "pass":
        return {}
    return cluster_data
```

- [ ] **Step 5: Run all adapter tests**

Run: `pytest tests/test_cluster_adapter.py -v`
Expected: PASS (all existing + 3 new tests)

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/adapters/cluster_adapter.py tests/test_cluster_adapter.py
git commit -m "feat: adapter reads versioned extraction runs and respects review status"
```

---

## Task 6: Migrate Tolonen Data and Re-run Extraction

**Files:**
- Run: `multiomics_kg/extraction/cluster/migrate.py`
- Run: `multiomics_kg/extraction/cluster/pipeline.py`

- [ ] **Step 1: Run migration on existing Tolonen JSONs**

```bash
uv run python -c "
from pathlib import Path
from multiomics_kg.extraction.cluster.migrate import migrate_all_tolonen
results = migrate_all_tolonen(Path('.'))
for r in results:
    print(f'Migrated: {r}')
"
```

Expected: Two run directories created under `data/Prochlorococcus/papers_and_supp/tolonen 2006/.extraction_cache/`

- [ ] **Step 2: Verify migration output**

```bash
ls -la "data/Prochlorococcus/papers_and_supp/tolonen 2006/.extraction_cache/med4_kmeans_nstarvation/current/"
ls -la "data/Prochlorococcus/papers_and_supp/tolonen 2006/.extraction_cache/mit9313_kmeans_nstarvation/current/"
```

Expected: Each directory contains `stage1_merged.json`, `stage2_results.json`, `stage3_validation.json`, `metadata.json`

- [ ] **Step 3: Re-run Tolonen extraction on current pipeline**

```bash
uv run python -m multiomics_kg.extraction.cluster.pipeline \
    "data/Prochlorococcus/papers_and_supp/tolonen 2006/paperconfig.yaml"
```

Expected: Two new timestamped run directories created. `current` symlinks updated. Each run contains all stage files + `report.md`.

- [ ] **Step 4: Verify new run output**

```bash
ls -la "data/Prochlorococcus/papers_and_supp/tolonen 2006/.extraction_cache/med4_kmeans_nstarvation/runs/"
cat "data/Prochlorococcus/papers_and_supp/tolonen 2006/.extraction_cache/med4_kmeans_nstarvation/current/report.md"
```

Expected: Two run directories (migrated + fresh). Report shows verdicts for all 9 MED4 clusters.

- [ ] **Step 5: Add `.extraction_cache` to .gitignore**

Check if `.gitignore` exists and add the pattern:

```
# Cluster extraction cache (versioned runs, pre-processing cache)
.extraction_cache/
```

- [ ] **Step 6: Commit**

```bash
git add .gitignore
git commit -m "chore: migrate Tolonen extraction data, re-run pipeline, gitignore cache"
```

---

## Task 7: Review Data Layer

**Files:**
- Create: `multiomics_kg/review/__init__.py`
- Create: `multiomics_kg/review/review_data.py`
- Test: `tests/test_review_data.py`

- [ ] **Step 1: Create package init**

```python
# multiomics_kg/review/__init__.py
```

- [ ] **Step 2: Write tests for data scanning**

```python
# tests/test_review_data.py
import json
from pathlib import Path
from multiomics_kg.extraction.cluster.run_manager import RunManager
from multiomics_kg.review.review_data import (
    scan_papers_with_clusters,
    load_entry_summary,
    compute_entry_status_color,
    export_issue_report,
)


def _setup_paper_with_extraction(tmp_path, paper_name, entry_key, clusters):
    """Create a minimal paperconfig + extraction cache for testing."""
    paper_dir = tmp_path / paper_name
    paper_dir.mkdir(parents=True)

    # Minimal paperconfig
    import yaml
    config = {
        "publication": {
            "papername": paper_name,
            "papermainpdf": "",
            "supplementary_materials": {
                entry_key: {
                    "type": "gene_clusters",
                    "filename": str(paper_dir / "clusters.csv"),
                    "organism": "Test Organism",
                    "gene_id_col": "gene",
                    "cluster_col": "cluster",
                }
            }
        }
    }
    with open(paper_dir / "paperconfig.yaml", "w") as f:
        yaml.dump(config, f)

    # Extraction cache
    cache_dir = paper_dir / ".extraction_cache"
    rm = RunManager(cache_dir, entry_key)
    run_dir = rm.create_run()

    stage1 = {str(k): {"enrichment_category": [{"value": f"cat_{k}"}]} for k in clusters}
    stage2 = {str(k): {"id": f"cluster_{k}", "functional_description": f"desc_{k}"} for k in clusters}
    stage3 = {str(k): {"verdict": v} for k, v in clusters.items()}

    rm.write_stage(run_dir, 1, stage1)
    rm.write_stage(run_dir, 2, stage2)
    rm.write_stage(run_dir, 3, stage3)

    return paper_dir, rm, run_dir


def test_scan_papers_finds_cluster_entries(tmp_path):
    """scan_papers_with_clusters returns paper/entry pairs."""
    paper_dir, _, _ = _setup_paper_with_extraction(
        tmp_path, "Paper1", "entry_a", {1: "pass", 2: "warn"}
    )
    results = scan_papers_with_clusters([paper_dir / "paperconfig.yaml"])
    assert len(results) == 1
    assert results[0]["paper_name"] == "Paper1"
    assert results[0]["entry_key"] == "entry_a"


def test_load_entry_summary(tmp_path):
    """load_entry_summary returns cluster count and verdict breakdown."""
    paper_dir, rm, run_dir = _setup_paper_with_extraction(
        tmp_path, "Paper1", "entry_a", {1: "pass", 2: "warn", 3: "fail"}
    )
    summary = load_entry_summary(paper_dir, "entry_a")
    assert summary["cluster_count"] == 3
    assert summary["verdicts"]["pass"] == 1
    assert summary["verdicts"]["warn"] == 1
    assert summary["verdicts"]["fail"] == 1


def test_compute_entry_status_color_all_approved(tmp_path):
    """Green when all clusters approved in current run."""
    paper_dir, rm, run_dir = _setup_paper_with_extraction(
        tmp_path, "Paper1", "entry_a", {1: "pass"}
    )
    review = {
        "1": {"status": "approve", "reviewed_in_run": run_dir.name,
              "input_hash": rm.compute_input_hash("1", rm.read_stage(run_dir, 1)),
              "issues": [], "failing_stages": [], "notes": "",
              "edited_fields": {}, "reviewed_at": "2026-04-04T12:00:00"},
    }
    rm.write_stage(run_dir, 4, review)
    color = compute_entry_status_color(paper_dir, "entry_a")
    assert color == "green"


def test_compute_entry_status_color_carried_forward(tmp_path):
    """Light green when all approved but from a previous run."""
    paper_dir, rm, run_dir = _setup_paper_with_extraction(
        tmp_path, "Paper1", "entry_a", {1: "pass"}
    )
    review = {
        "1": {"status": "approve", "reviewed_in_run": "some-old-run",
              "input_hash": rm.compute_input_hash("1", rm.read_stage(run_dir, 1)),
              "issues": [], "failing_stages": [], "notes": "",
              "edited_fields": {}, "reviewed_at": "2026-04-04T12:00:00"},
    }
    rm.write_stage(run_dir, 4, review)
    color = compute_entry_status_color(paper_dir, "entry_a")
    assert color == "light_green"


def test_compute_entry_status_color_unreviewed(tmp_path):
    """Red when no review exists."""
    paper_dir, _, _ = _setup_paper_with_extraction(
        tmp_path, "Paper1", "entry_a", {1: "pass"}
    )
    color = compute_entry_status_color(paper_dir, "entry_a")
    assert color == "red"


def test_export_issue_report(tmp_path):
    """Export produces CSV rows with issue details."""
    paper_dir, rm, run_dir = _setup_paper_with_extraction(
        tmp_path, "Paper1", "entry_a", {1: "warn"}
    )
    review = {
        "1": {"status": "flag-issue", "issues": ["hallucinated", "cross_contamination"],
              "failing_stages": ["stage1_visual"], "notes": "Bad quote",
              "input_hash": "abc", "edited_fields": {},
              "reviewed_at": "2026-04-04T12:00:00", "reviewed_in_run": run_dir.name},
    }
    rm.write_stage(run_dir, 4, review)
    rows = export_issue_report([(paper_dir, "entry_a")])
    assert len(rows) == 1
    assert rows[0]["issues"] == "hallucinated,cross_contamination"
    assert rows[0]["failing_stages"] == "stage1_visual"
    assert rows[0]["notes"] == "Bad quote"
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `pytest tests/test_review_data.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'multiomics_kg.review.review_data'`

- [ ] **Step 4: Implement review_data.py**

```python
# multiomics_kg/review/review_data.py
"""Data loading and scanning for the cluster review UI."""

import csv
import io
import json
import logging
from pathlib import Path
from typing import Optional

import yaml

from multiomics_kg.extraction.cluster.run_manager import RunManager

logger = logging.getLogger(__name__)


def scan_papers_with_clusters(
    paperconfig_paths: list[Path],
) -> list[dict]:
    """Scan paperconfigs and return entries that have gene_clusters.

    Returns list of dicts with keys: paper_name, entry_key, paperconfig_path, paper_dir.
    """
    results = []
    for pc_path in paperconfig_paths:
        pc_path = Path(pc_path)
        if not pc_path.exists():
            continue
        with open(pc_path) as f:
            config = yaml.safe_load(f)
        paper_name = config.get("publication", {}).get("papername", pc_path.parent.name)
        supp = config.get("publication", {}).get("supplementary_materials", {})
        for key, entry in supp.items():
            if entry.get("type") == "gene_clusters":
                results.append({
                    "paper_name": paper_name,
                    "entry_key": key,
                    "paperconfig_path": pc_path,
                    "paper_dir": pc_path.parent,
                    "organism": entry.get("organism", ""),
                })
    return results


def load_entry_summary(paper_dir: Path, entry_key: str) -> dict:
    """Load summary stats for one extraction entry.

    Returns dict with cluster_count, verdicts (pass/warn/fail counts),
    review_statuses (approve/edit/reject/flag-issue/stale/unreviewed counts).
    """
    cache_dir = paper_dir / ".extraction_cache"
    rm = RunManager(cache_dir, entry_key)
    run_dir = rm.get_current_run()

    if run_dir is None:
        return {"cluster_count": 0, "verdicts": {}, "review_statuses": {}, "has_run": False}

    stage2 = rm.read_stage(run_dir, 2)
    stage3 = rm.read_stage(run_dir, 3)
    stage4 = rm.read_stage(run_dir, 4)

    cluster_keys = sorted(set(list(stage2.keys()) + list(stage3.keys())))
    verdicts: dict[str, int] = {}
    review_statuses: dict[str, int] = {}

    for ck in cluster_keys:
        v = stage3.get(ck, {}).get("verdict", "none")
        verdicts[v] = verdicts.get(v, 0) + 1

        r = stage4.get(ck, {}).get("status", "unreviewed")
        review_statuses[r] = review_statuses.get(r, 0) + 1

    return {
        "cluster_count": len(cluster_keys),
        "verdicts": verdicts,
        "review_statuses": review_statuses,
        "has_run": True,
    }


def compute_entry_status_color(paper_dir: Path, entry_key: str) -> str:
    """Compute rollup color for an entry based on review state.

    Returns: 'green', 'light_green', 'yellow', 'red'
    """
    cache_dir = paper_dir / ".extraction_cache"
    rm = RunManager(cache_dir, entry_key)
    run_dir = rm.get_current_run()

    if run_dir is None:
        return "red"

    stage2 = rm.read_stage(run_dir, 2)
    stage4 = rm.read_stage(run_dir, 4)
    cluster_keys = list(stage2.keys())

    if not cluster_keys:
        return "red"

    all_approved = True
    any_current_run_review = False
    any_stale = False

    for ck in cluster_keys:
        review = stage4.get(ck, {})
        status = review.get("status", "")
        if status not in ("approve", "edit"):
            all_approved = False
        if status == "stale":
            any_stale = True
        if review.get("reviewed_in_run") == run_dir.name:
            any_current_run_review = True

    if not all_approved:
        if any_stale:
            return "yellow"
        return "red"

    # All approved — check if any were reviewed in current run
    has_all_current = all(
        stage4.get(ck, {}).get("reviewed_in_run") == run_dir.name
        for ck in cluster_keys
    )
    if has_all_current:
        return "green"
    return "light_green"


def export_issue_report(entries: list[tuple[Path, str]]) -> list[dict]:
    """Export issue report rows across all entries.

    Returns list of dicts suitable for CSV export.
    """
    rows = []
    for paper_dir, entry_key in entries:
        cache_dir = paper_dir / ".extraction_cache"
        rm = RunManager(cache_dir, entry_key)
        run_dir = rm.get_current_run()
        if run_dir is None:
            continue

        stage3 = rm.read_stage(run_dir, 3)
        stage4 = rm.read_stage(run_dir, 4)
        metadata_path = run_dir / "metadata.json"
        paper_name = entry_key
        if metadata_path.exists():
            with open(metadata_path) as f:
                meta = json.load(f)
            paper_name = meta.get("paper", entry_key)

        for ck in sorted(stage4.keys()):
            review = stage4[ck]
            if review.get("status") in ("approve",):
                continue  # skip clean approvals
            verdict = stage3.get(ck, {}).get("verdict", "")
            rows.append({
                "paper": paper_name,
                "entry": entry_key,
                "cluster": ck,
                "verdict": verdict,
                "review_status": review.get("status", ""),
                "issues": ",".join(review.get("issues", [])),
                "failing_stages": ",".join(review.get("failing_stages", [])),
                "notes": review.get("notes", ""),
            })
    return rows
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_review_data.py -v`
Expected: PASS (6 tests)

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/review/__init__.py multiomics_kg/review/review_data.py tests/test_review_data.py
git commit -m "feat: add review data layer for scanning papers and loading extraction state"
```

---

## Task 8: Streamlit Review App — Sidebar and Navigation

**Files:**
- Create: `multiomics_kg/review/cluster_review_app.py`
- Modify: `pyproject.toml` (add streamlit dependency)

- [ ] **Step 1: Add streamlit dependency**

Add to `pyproject.toml` dependencies:

```
"streamlit>=1.30.0",
"openpyxl>=3.1.0",
```

Then run:

```bash
uv sync
```

- [ ] **Step 2: Create the Streamlit app with sidebar navigation**

```python
# multiomics_kg/review/cluster_review_app.py
"""Streamlit app for reviewing cluster extraction results."""

import json
import logging
from pathlib import Path

import streamlit as st
import yaml

from multiomics_kg.extraction.cluster.run_manager import RunManager
from multiomics_kg.review.review_data import (
    compute_entry_status_color,
    export_issue_report,
    load_entry_summary,
    scan_papers_with_clusters,
)

logger = logging.getLogger(__name__)

# Color mapping for status indicators
STATUS_COLORS = {
    "green": "\U0001f7e2",       # green circle
    "light_green": "\U0001f7e1", # yellow circle (closest to light green)
    "yellow": "\U0001f7e0",      # orange circle
    "red": "\U0001f534",         # red circle
}

REVIEW_STATUS_ICONS = {
    "approve": "\u2705",
    "edit": "\u270f\ufe0f",
    "reject": "\u274c",
    "flag-issue": "\u26a0\ufe0f",
    "stale": "\U0001f504",
    "unreviewed": "\u2b1c",
}


def find_paperconfig_paths() -> list[Path]:
    """Find all paperconfig list files and collect paths."""
    project_root = Path(__file__).parent.parent.parent
    paths = []
    for list_file in [
        project_root / "data" / "Prochlorococcus" / "papers_and_supp" / "paperconfig_files.txt",
        project_root / "data" / "Synechococcus" / "papers_and_supp" / "paperconfig_files.txt",
    ]:
        if list_file.exists():
            for line in list_file.read_text().splitlines():
                line = line.strip()
                if line and not line.startswith("#"):
                    pc_path = project_root / line
                    if pc_path.exists():
                        paths.append(pc_path)
    return paths


def render_sidebar():
    """Render sidebar with paper/entry/cluster navigation."""
    st.sidebar.title("Cluster Review")

    # Scan for papers with gene_clusters
    if "papers" not in st.session_state:
        pc_paths = find_paperconfig_paths()
        st.session_state.papers = scan_papers_with_clusters(pc_paths)

    papers = st.session_state.papers
    if not papers:
        st.sidebar.warning("No gene_clusters entries found in paperconfigs.")
        return None

    # Group by paper
    paper_names = sorted(set(p["paper_name"] for p in papers))

    # Paper selector with color coding
    paper_options = []
    for pn in paper_names:
        entries = [p for p in papers if p["paper_name"] == pn]
        worst_color = "green"
        for e in entries:
            c = compute_entry_status_color(e["paper_dir"], e["entry_key"])
            if c == "red":
                worst_color = "red"
                break
            elif c == "yellow" and worst_color != "red":
                worst_color = "yellow"
            elif c == "light_green" and worst_color == "green":
                worst_color = "light_green"
        icon = STATUS_COLORS.get(worst_color, "")
        paper_options.append(f"{icon} {pn}")

    selected_paper_label = st.sidebar.selectbox("Paper", paper_options)
    if not selected_paper_label:
        return None
    selected_paper = selected_paper_label.split(" ", 1)[1]

    # Entry selector
    entries = [p for p in papers if p["paper_name"] == selected_paper]
    entry_options = []
    for e in entries:
        color = compute_entry_status_color(e["paper_dir"], e["entry_key"])
        icon = STATUS_COLORS.get(color, "")
        entry_options.append(f"{icon} {e['entry_key']}")

    selected_entry_label = st.sidebar.selectbox("Analysis Entry", entry_options)
    if not selected_entry_label:
        return None
    selected_entry_key = selected_entry_label.split(" ", 1)[1]

    selected = next(
        e for e in entries if e["entry_key"] == selected_entry_key
    )

    # Summary stats
    summary = load_entry_summary(selected["paper_dir"], selected["entry_key"])
    if summary["has_run"]:
        st.sidebar.markdown("---")
        st.sidebar.markdown(f"**Clusters:** {summary['cluster_count']}")
        verdict_parts = [f"{v}: {c}" for v, c in sorted(summary["verdicts"].items())]
        st.sidebar.markdown(f"**Verdicts:** {', '.join(verdict_parts)}")
        review_parts = [
            f"{REVIEW_STATUS_ICONS.get(s, '')} {s}: {c}"
            for s, c in sorted(summary["review_statuses"].items())
        ]
        st.sidebar.markdown(f"**Reviews:** {', '.join(review_parts)}")
    else:
        st.sidebar.warning("No extraction runs found.")

    # Export button
    st.sidebar.markdown("---")
    if st.sidebar.button("Export Issue Report"):
        all_entries = [(p["paper_dir"], p["entry_key"]) for p in papers]
        rows = export_issue_report(all_entries)
        if rows:
            import csv
            import io
            buf = io.StringIO()
            writer = csv.DictWriter(buf, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
            st.sidebar.download_button(
                "Download CSV", buf.getvalue(), "issue_report.csv", "text/csv"
            )
        else:
            st.sidebar.info("No issues to export.")

    return selected


def main():
    st.set_page_config(page_title="Cluster Extraction Review", layout="wide")
    selected = render_sidebar()

    if selected is None:
        st.title("Cluster Extraction Review")
        st.info("Select a paper and analysis entry from the sidebar.")
        return

    st.title(f"{selected['paper_name']} — {selected['entry_key']}")
    st.caption(f"Organism: {selected['organism']}")

    # Load current run data
    cache_dir = selected["paper_dir"] / ".extraction_cache"
    rm = RunManager(cache_dir, selected["entry_key"])
    run_dir = rm.get_current_run()

    if run_dir is None:
        st.warning("No extraction run found. Run the pipeline first.")
        return

    stage1 = rm.read_stage(run_dir, 1)
    stage2 = rm.read_stage(run_dir, 2)
    stage3 = rm.read_stage(run_dir, 3)
    stage4 = rm.read_stage(run_dir, 4)

    # Store in session state for components to access
    st.session_state.current_run_dir = run_dir
    st.session_state.current_rm = rm
    st.session_state.stage1 = stage1
    st.session_state.stage2 = stage2
    st.session_state.stage3 = stage3
    st.session_state.stage4 = stage4
    st.session_state.paper_dir = selected["paper_dir"]

    # Placeholder for center pane (Task 9 will add components)
    cluster_keys = sorted(stage1.keys(), key=lambda k: (not k.isdigit(), k))
    st.markdown(f"**Run:** `{run_dir.name}` | **Clusters:** {len(cluster_keys)}")

    for ck in cluster_keys:
        review = stage4.get(ck, {})
        review_status = review.get("status", "unreviewed")
        reviewed_in = review.get("reviewed_in_run", "")
        is_current = reviewed_in == run_dir.name
        icon = REVIEW_STATUS_ICONS.get(review_status, REVIEW_STATUS_ICONS["unreviewed"])
        if review_status in ("approve", "edit") and not is_current:
            icon = "\U0001f7e1"  # carried forward indicator

        st.markdown(f"### {icon} Cluster {ck}")
        st.caption(f"Review: {review_status}" + ("" if is_current else " (carried forward)"))

        # Placeholder text — will be replaced by components in Task 9
        s2 = stage2.get(ck, {})
        st.text(f"ID: {s2.get('id', '—')} | Name: {s2.get('name', '—')}")
        verdict = stage3.get(ck, {}).get("verdict", "none")
        st.text(f"Verdict: {verdict}")
        st.markdown("---")


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Test the app launches**

```bash
uv run streamlit run multiomics_kg/review/cluster_review_app.py --server.headless true &
sleep 3
curl -s http://localhost:8501 | head -5
kill %1
```

Expected: HTML response from Streamlit (confirms app loads without crash)

- [ ] **Step 4: Commit**

```bash
git add multiomics_kg/review/cluster_review_app.py pyproject.toml uv.lock
git commit -m "feat: Streamlit review app with sidebar navigation and paper scanning"
```

---

## Task 9: Review Components — Merge View and Review Controls

**Files:**
- Create: `multiomics_kg/review/review_components.py`
- Modify: `multiomics_kg/review/cluster_review_app.py`

- [ ] **Step 1: Create review components module**

```python
# multiomics_kg/review/review_components.py
"""Streamlit UI components for cluster review: merge view, review controls, source viewer."""

import json
import subprocess
import webbrowser
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd
import streamlit as st

from multiomics_kg.extraction.cluster.run_manager import RunManager

# Fields displayed in the merge view
MERGE_FIELDS = [
    "enrichment_category",
    "enrichment_details",
    "direction",
    "temporal_pattern",
    "cluster_description",
    "treatment_conditions",
    "peak_time",
    "period_description",
    "light_phase",
]

CONFIDENCE_COLORS = {
    "very_high": "#2e7d32",  # dark green
    "high": "#1565c0",       # blue
    "medium": "#f57f17",     # amber
    "low": "#c62828",        # red
}

ISSUE_OPTIONS = [
    "wrong_cluster",
    "hallucinated",
    "low_info",
    "cross_contamination",
    "partial",
    "source_missing",
]

FAILING_STAGE_OPTIONS = [
    "stage1_visual",
    "stage1_semantic",
    "stage1_table",
    "stage2_synthesis",
    "stage3_validation",
    "merge",
]

REVIEW_STATUS_OPTIONS = ["approve", "edit", "reject", "flag-issue"]


def render_merge_view(cluster_key: str, stage1_cluster: dict, stage2_cluster: dict):
    """Render the merge-centric view for one cluster.

    Shows per-field candidates from each path with confidence badges,
    then the synthesis result.
    """
    for field in MERGE_FIELDS:
        candidates = stage1_cluster.get(field, [])
        if not candidates:
            continue

        with st.expander(f"**{field.replace('_', ' ').title()}**", expanded=False):
            if isinstance(candidates, list) and candidates and isinstance(candidates[0], dict):
                for c in candidates:
                    source = c.get("source", "?")
                    conf = c.get("confidence", "?")
                    value = c.get("value", "")
                    color = CONFIDENCE_COLORS.get(conf, "#666")
                    st.markdown(
                        f"<span style='color:{color};font-weight:bold'>[{source}|{conf}]</span> {value}",
                        unsafe_allow_html=True,
                    )
            else:
                st.text(str(candidates))

    # Synthesis result
    st.markdown("#### Synthesis Result")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown(f"**ID:** `{stage2_cluster.get('id', '—')}`")
        st.markdown(f"**Name:** {stage2_cluster.get('name', '—')}")
        if stage2_cluster.get("peak_time_hours") is not None:
            st.markdown(f"**Peak time:** {stage2_cluster['peak_time_hours']}h")
        if stage2_cluster.get("period_hours") is not None:
            st.markdown(f"**Period:** {stage2_cluster['period_hours']}h")
    with col2:
        st.markdown("**Functional description:**")
        st.text(stage2_cluster.get("functional_description", "—"))
        st.markdown("**Behavioral description:**")
        st.text(stage2_cluster.get("behavioral_description", "—"))


def render_quotes(cluster_key: str, stage1_cluster: dict):
    """Render supporting quotes with source labels and page references."""
    quotes = stage1_cluster.get("supporting_quotes", [])
    if not quotes:
        st.caption("No supporting quotes.")
        return

    with st.expander(f"Supporting quotes ({len(quotes)})", expanded=False):
        for q in quotes:
            source = q.get("source", "?")
            score = q.get("relevance_score", "?")
            location = q.get("location", "")
            quote_text = q.get("quote", "")
            st.markdown(
                f"**[{source}]** (relevance: {score}) {f'— {location}' if location else ''}"
            )
            st.markdown(f"> {quote_text}")


def render_verdict(cluster_key: str, stage3_cluster: dict):
    """Render stage 3 validation verdict."""
    verdict = stage3_cluster.get("verdict", "none")
    explanation = stage3_cluster.get("explanation", "")
    colors = {"pass": "green", "warn": "orange", "fail": "red", "none": "gray"}
    color = colors.get(verdict, "gray")
    st.markdown(
        f"**Verdict:** <span style='color:{color};font-weight:bold'>{verdict.upper()}</span> — {explanation}",
        unsafe_allow_html=True,
    )


def render_review_controls(
    cluster_key: str,
    stage2_cluster: dict,
    current_review: dict,
    run_dir: Path,
    rm: RunManager,
    stage1: dict,
):
    """Render review controls for one cluster: status, issues, failing stages, notes, editable fields."""
    prefix = f"review_{cluster_key}"

    # Status
    current_status = current_review.get("status", "")
    status_idx = REVIEW_STATUS_OPTIONS.index(current_status) if current_status in REVIEW_STATUS_OPTIONS else None
    status = st.selectbox(
        "Review status",
        REVIEW_STATUS_OPTIONS,
        index=status_idx,
        placeholder="Select...",
        key=f"{prefix}_status",
    )

    # Issue classification
    current_issues = current_review.get("issues", [])
    issues = st.multiselect(
        "Issue classification",
        ISSUE_OPTIONS,
        default=[i for i in current_issues if i in ISSUE_OPTIONS],
        key=f"{prefix}_issues",
    )

    # Failing stages
    current_stages = current_review.get("failing_stages", [])
    failing_stages = st.multiselect(
        "Failing stage(s)",
        FAILING_STAGE_OPTIONS,
        default=[s for s in current_stages if s in FAILING_STAGE_OPTIONS],
        key=f"{prefix}_failing_stages",
    )

    # Notes
    notes = st.text_area(
        "Notes",
        value=current_review.get("notes", ""),
        key=f"{prefix}_notes",
    )

    # Editable fields (when status is 'edit')
    edited_fields = dict(current_review.get("edited_fields", {}))
    if status == "edit":
        st.markdown("**Edit descriptions:**")
        edited_fields["functional_description"] = st.text_area(
            "Functional description",
            value=edited_fields.get(
                "functional_description",
                stage2_cluster.get("functional_description", ""),
            ),
            key=f"{prefix}_func_desc",
        )
        edited_fields["behavioral_description"] = st.text_area(
            "Behavioral description",
            value=edited_fields.get(
                "behavioral_description",
                stage2_cluster.get("behavioral_description", ""),
            ),
            key=f"{prefix}_behav_desc",
        )

    # Save button
    if st.button("Save review", key=f"{prefix}_save"):
        review_entry = {
            "status": status,
            "issues": issues,
            "failing_stages": failing_stages,
            "notes": notes,
            "edited_fields": edited_fields if status == "edit" else {},
            "reviewed_at": datetime.now().isoformat(timespec="seconds"),
            "reviewed_in_run": run_dir.name,
            "input_hash": rm.compute_input_hash(cluster_key, stage1),
        }
        # Read current stage4, update this cluster, write back
        stage4 = rm.read_stage(run_dir, 4)
        stage4[cluster_key] = review_entry
        rm.write_stage(run_dir, 4, stage4)
        st.success(f"Review saved for cluster {cluster_key}")
        st.rerun()


def render_source_access(paper_dir: Path, paperconfig_path: Path):
    """Render source file access: inline tables, PDF links, xdg-open for others."""
    with open(paperconfig_path) as f:
        config = yaml.safe_load(f)

    st.markdown("### Source Files")

    # Main PDF
    pdf_rel = config.get("publication", {}).get("papermainpdf", "")
    if pdf_rel:
        project_root = Path(__file__).parent.parent.parent
        pdf_path = project_root / pdf_rel
        if pdf_path.exists():
            if st.button("Open main PDF", key="open_pdf"):
                webbrowser.open(f"file://{pdf_path.resolve()}")

    # List paper directory files
    for f_path in sorted(paper_dir.iterdir()):
        if f_path.name.startswith(".") or f_path.is_dir():
            continue
        suffix = f_path.suffix.lower()

        if suffix in (".csv", ".tsv"):
            with st.expander(f"Table: {f_path.name}"):
                try:
                    sep = "\t" if suffix == ".tsv" else ","
                    df = pd.read_csv(f_path, sep=sep, nrows=200)
                    st.dataframe(df, use_container_width=True)
                except Exception as e:
                    st.error(f"Could not read {f_path.name}: {e}")

        elif suffix in (".xlsx", ".xls"):
            with st.expander(f"Spreadsheet: {f_path.name}"):
                try:
                    xls = pd.ExcelFile(f_path)
                    sheet = st.selectbox(
                        "Sheet", xls.sheet_names, key=f"sheet_{f_path.name}"
                    )
                    df = pd.read_excel(f_path, sheet_name=sheet, nrows=200)
                    st.dataframe(df, use_container_width=True)
                except Exception as e:
                    st.error(f"Could not read {f_path.name}: {e}")

        elif suffix == ".pdf":
            if st.button(f"Open {f_path.name}", key=f"open_{f_path.name}"):
                webbrowser.open(f"file://{f_path.resolve()}")

        else:
            if st.button(f"Open {f_path.name}", key=f"open_{f_path.name}"):
                subprocess.Popen(["xdg-open", str(f_path)])
```

- [ ] **Step 2: Add missing import to review_components.py**

Add at top of file:
```python
import yaml
```

- [ ] **Step 3: Wire components into the main app**

Replace the placeholder cluster loop in `cluster_review_app.py` `main()` (the section after `st.markdown(f"**Run:** ...")`) with:

```python
    from multiomics_kg.review.review_components import (
        render_merge_view,
        render_quotes,
        render_review_controls,
        render_source_access,
        render_verdict,
    )

    # Source access in a sidebar expander or tab
    tab_clusters, tab_sources = st.tabs(["Clusters", "Sources"])

    with tab_sources:
        render_source_access(selected["paper_dir"], selected["paperconfig_path"])

    with tab_clusters:
        cluster_keys = sorted(stage1.keys(), key=lambda k: (not k.isdigit(), k))
        st.markdown(f"**Run:** `{run_dir.name}` | **Clusters:** {len(cluster_keys)}")

        for ck in cluster_keys:
            review = stage4.get(ck, {})
            review_status = review.get("status", "unreviewed")
            reviewed_in = review.get("reviewed_in_run", "")
            is_current = reviewed_in == run_dir.name

            icon = REVIEW_STATUS_ICONS.get(review_status, REVIEW_STATUS_ICONS["unreviewed"])
            if review_status in ("approve", "edit") and not is_current:
                icon = "\U0001f7e1"

            with st.expander(
                f"{icon} Cluster {ck} — {stage2.get(ck, {}).get('name', '(no name)')}",
                expanded=(review_status not in ("approve",)),
            ):
                st.caption(
                    f"Review: {review_status}"
                    + ("" if is_current or not reviewed_in else " (carried forward)")
                )

                # Merge view
                render_merge_view(ck, stage1.get(ck, {}), stage2.get(ck, {}))

                # Quotes
                render_quotes(ck, stage1.get(ck, {}))

                # Verdict
                render_verdict(ck, stage3.get(ck, {}))

                # Review controls
                st.markdown("---")
                render_review_controls(
                    ck, stage2.get(ck, {}), review, run_dir, rm, stage1
                )
```

- [ ] **Step 4: Test the app loads with Tolonen data**

```bash
uv run streamlit run multiomics_kg/review/cluster_review_app.py
```

Expected: App loads in browser. Sidebar shows Tolonen 2006 with color indicator. Selecting it shows 9 MED4 clusters (or 7 MIT9313 clusters) with merge view, quotes, verdict, and review controls.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/review/review_components.py multiomics_kg/review/cluster_review_app.py
git commit -m "feat: merge-centric cluster view with review controls and source access"
```

---

## Task 10: Diff View for Comparing Runs

**Files:**
- Modify: `multiomics_kg/review/review_components.py`
- Modify: `multiomics_kg/review/cluster_review_app.py`

- [ ] **Step 1: Add diff view component**

Add to `multiomics_kg/review/review_components.py`:

```python
def render_diff_view(rm: RunManager, cluster_key: str):
    """Show side-by-side diff between current and previous run's Stage 2 results."""
    runs = rm.list_runs()
    if len(runs) < 2:
        st.caption("No previous run to compare.")
        return

    current_run = rm.get_current_run()
    # Default to the second-most-recent run
    prev_options = [r for r in runs if r != current_run]
    if not prev_options:
        st.caption("No previous run to compare.")
        return

    prev_run_name = st.selectbox(
        "Compare with run:",
        [r.name for r in reversed(prev_options)],
        key=f"diff_run_{cluster_key}",
    )
    prev_run = rm.runs_dir / prev_run_name

    old_s2 = rm.read_stage(prev_run, 2).get(cluster_key, {})
    new_s2 = rm.read_stage(current_run, 2).get(cluster_key, {})

    fields = ["id", "name", "functional_description", "behavioral_description",
              "peak_time_hours", "period_hours"]

    for field in fields:
        old_val = str(old_s2.get(field, "—"))
        new_val = str(new_s2.get(field, "—"))
        if old_val == new_val:
            color = "#888"  # gray — unchanged
            label = "="
        else:
            color = "#c62828"  # red — changed
            label = "CHANGED"
        st.markdown(
            f"**{field}** <span style='color:{color}'>[{label}]</span>",
            unsafe_allow_html=True,
        )
        if old_val != new_val:
            col_old, col_new = st.columns(2)
            with col_old:
                st.markdown(f"**Old:** {old_val}")
            with col_new:
                st.markdown(f"**New:** {new_val}")
```

- [ ] **Step 2: Add diff tab to cluster expander in main app**

In `cluster_review_app.py`, inside the cluster expander (after `render_review_controls`), add:

```python
                # Diff view
                if len(rm.list_runs()) > 1:
                    with st.expander("Compare with previous run"):
                        render_diff_view(rm, ck)
```

Add `render_diff_view` to the import from `review_components`.

- [ ] **Step 3: Test diff view with Tolonen data**

```bash
uv run streamlit run multiomics_kg/review/cluster_review_app.py
```

Expected: After the Tolonen re-run (Task 6 created a second run), each cluster shows a "Compare with previous run" expander showing field-level diffs between the migrated and fresh extraction.

- [ ] **Step 4: Commit**

```bash
git add multiomics_kg/review/review_components.py multiomics_kg/review/cluster_review_app.py
git commit -m "feat: add diff view for comparing extraction runs"
```

---

## Task 11: Re-run Triggers from UI

**Files:**
- Modify: `multiomics_kg/review/cluster_review_app.py`

- [ ] **Step 1: Add re-run buttons to the app**

Add a re-run section to the main app, after the tabs in `main()`:

```python
    # Re-run controls
    st.sidebar.markdown("---")
    st.sidebar.markdown("### Re-run Extraction")

    force_mode = st.sidebar.selectbox(
        "Force mode",
        ["Default (cached pre-processing)", "Force LLM (re-run all LLM)", "Force all (clear cache)"],
        key="force_mode",
    )

    if st.sidebar.button("Re-run this entry"):
        force_llm = "Force LLM" in force_mode
        force_all = "Force all" in force_mode
        with st.spinner(f"Running extraction for {selected['entry_key']}..."):
            from multiomics_kg.extraction.cluster.pipeline import run_pipeline
            results = run_pipeline(
                selected["paperconfig_path"],
                table_key=selected["entry_key"],
                force_llm=force_llm,
                force_all=force_all,
            )
            st.sidebar.success(f"Done: {results[0].name if results else 'no output'}")
            # Clear cached data to force reload
            for key in ["papers", "stage1", "stage2", "stage3", "stage4"]:
                st.session_state.pop(key, None)
            st.rerun()
```

- [ ] **Step 2: Add per-cluster re-synthesis button**

In `review_components.py`, add to `render_review_controls()` before the save button:

```python
    # Re-run buttons
    col_resyn, col_full = st.columns(2)
    with col_resyn:
        if st.button("Re-synthesize (Stage 2+3)", key=f"{prefix}_resynthesize"):
            with st.spinner("Re-running synthesis..."):
                from multiomics_kg.extraction.cluster.pipeline import run_pipeline
                run_pipeline(
                    st.session_state.get("paperconfig_path"),
                    table_key=st.session_state.get("entry_key"),
                    from_stage=2,
                )
                st.rerun()
    with col_full:
        if st.button("Full re-extract", key=f"{prefix}_full_reextract"):
            with st.spinner("Re-running full extraction..."):
                from multiomics_kg.extraction.cluster.pipeline import run_pipeline
                run_pipeline(
                    st.session_state.get("paperconfig_path"),
                    table_key=st.session_state.get("entry_key"),
                )
                st.rerun()
```

Update `cluster_review_app.py` to store `paperconfig_path` and `entry_key` in session state:

```python
    st.session_state.paperconfig_path = selected["paperconfig_path"]
    st.session_state.entry_key = selected["entry_key"]
```

- [ ] **Step 3: Test re-run from UI**

```bash
uv run streamlit run multiomics_kg/review/cluster_review_app.py
```

Expected: Sidebar shows re-run button with force mode dropdown. Per-cluster re-synthesize and full re-extract buttons appear in each cluster card. Clicking triggers extraction and reloads the page.

- [ ] **Step 4: Commit**

```bash
git add multiomics_kg/review/cluster_review_app.py multiomics_kg/review/review_components.py
git commit -m "feat: add re-run triggers from review UI (entry-level and per-cluster)"
```

---

## Task 12: Run All Tests and Final Verification

- [ ] **Step 1: Run all unit tests**

```bash
pytest tests/test_run_manager.py tests/test_cluster_extraction_migration.py tests/test_review_data.py tests/test_cluster_adapter.py tests/test_extraction_merge.py tests/test_extraction_table.py -v
```

Expected: All tests pass.

- [ ] **Step 2: Run full test suite (excluding slow and kg)**

```bash
pytest -m "not slow and not kg" -v
```

Expected: No regressions.

- [ ] **Step 3: Verify Streamlit app works end-to-end**

```bash
uv run streamlit run multiomics_kg/review/cluster_review_app.py
```

Manual verification:
- Sidebar shows Tolonen 2006 with color-coded entries
- Selecting MED4 entry shows 9 clusters with merge view
- Each cluster shows field candidates from all paths with confidence
- Supporting quotes are visible with source labels
- Review controls work: set status, select issues, select failing stages, add notes, save
- Saved reviews persist (reload page, review still there)
- Sources tab shows inline CSV/XLS tables, PDF open links
- Diff view shows comparison between migrated and fresh extraction runs
- Re-run button triggers extraction and reloads

- [ ] **Step 4: Commit any fixes**

```bash
git add -u
git commit -m "fix: address issues found during final verification"
```

---

## Summary

| Task | What it builds | Key files |
|---|---|---|
| 1 | RunManager — directory creation, stage I/O | `run_manager.py`, `test_run_manager.py` |
| 2 | Input hashing, review copy-forward | `run_manager.py`, `test_run_manager.py` |
| 3 | Migration of legacy JSONs | `migrate.py`, `test_cluster_extraction_migration.py` |
| 4 | Pipeline refactor to use RunManager | `pipeline.py` |
| 5 | Adapter reads new structure + review status | `cluster_adapter.py`, `test_cluster_adapter.py` |
| 6 | Migrate Tolonen data + re-run | (runtime) |
| 7 | Review data layer (scanning, summaries, export) | `review_data.py`, `test_review_data.py` |
| 8 | Streamlit app skeleton + sidebar | `cluster_review_app.py` |
| 9 | Merge view, review controls, source access | `review_components.py` |
| 10 | Diff view between runs | `review_components.py` |
| 11 | Re-run triggers from UI | `cluster_review_app.py` |
| 12 | Final verification | (testing) |
