"""Tests for RunManager — versioned extraction run directories and stage file I/O."""

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


STAGE_FILENAMES = {
    1: "stage1_merged.json",
    2: "stage2_results.json",
    3: "stage3_validation.json",
    4: "stage4_review.json",
}


def test_stage_filenames(tmp_path):
    """Each stage number maps to the correct filename."""
    rm = RunManager(tmp_path, "entry1")
    run_dir = rm.create_run()
    for stage_num, expected_name in STAGE_FILENAMES.items():
        rm.write_stage(run_dir, stage_num, {"test": stage_num})
        assert (run_dir / expected_name).exists()


def test_list_runs_returns_oldest_first(tmp_path):
    """list_runs returns all runs sorted oldest-first by timestamp."""
    rm = RunManager(tmp_path, "entry1")
    run1 = rm.create_run()
    time.sleep(1.1)  # ensure different timestamp
    run2 = rm.create_run()
    time.sleep(1.1)
    run3 = rm.create_run()
    runs = rm.list_runs()
    assert len(runs) == 3
    assert runs == [run1, run2, run3]
    # Verify oldest-first ordering
    assert runs[0].name < runs[1].name < runs[2].name


def test_read_current_stage_no_runs_returns_empty(tmp_path):
    """read_current_stage returns {} when no runs exist."""
    rm = RunManager(tmp_path, "entry1")
    assert rm.read_current_stage(1) == {}
    assert rm.read_current_stage(3) == {}


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


def test_copy_forward_reviews_unchanged_input(tmp_path):
    """Reviews carry forward with same status when input hash matches."""
    rm = RunManager(tmp_path, "entry1")
    stage1 = {"1": {"enrichment_category": [{"value": "Transport"}]}}

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

    time.sleep(1.1)
    run2 = rm.create_run()
    rm.write_stage(run2, 1, stage1)
    rm.copy_forward_reviews(run1, run2, stage1)

    new_review = rm.read_stage(run2, 4)
    assert new_review["1"]["status"] == "approve"
    assert new_review["1"]["reviewed_in_run"] == run1.name


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
    stage1_new = {"1": {"x": []}}

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
