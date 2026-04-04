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
    paper_dir = tmp_path / "paper"
    paper_dir.mkdir()
    cache_dir = paper_dir / ".extraction_cache"
    result = migrate_extraction_json(paper_dir, "missing_entry", cache_dir)
    assert result is None


def test_migrate_stores_metadata(tmp_path):
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
