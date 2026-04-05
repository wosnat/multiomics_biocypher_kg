import json
import pytest
from pathlib import Path


def test_load_extraction(tmp_path):
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
    from multiomics_kg.extraction.cluster.extraction_utils import load_extraction
    ext_dir = tmp_path / "cluster_extractions"
    ext_dir.mkdir()
    old = {"stage2_results": {"1": {"name": "old"}}, "stage3_validation": {}}
    (ext_dir / "old_entry.json").write_text(json.dumps(old))
    result = load_extraction(tmp_path, "old_entry")
    assert result == {}


def test_load_missing_file(tmp_path):
    from multiomics_kg.extraction.cluster.extraction_utils import load_extraction
    result = load_extraction(tmp_path, "nonexistent")
    assert result == {}


def test_save_load_roundtrip(tmp_path):
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
    md_path = tmp_path / "cluster_extractions" / "test_entry.md"
    assert md_path.exists()
    assert "Test cluster 1" in md_path.read_text()


def test_list_extraction_files(tmp_path):
    from multiomics_kg.extraction.cluster.extraction_utils import list_extraction_files
    ext_dir = tmp_path / "cluster_extractions"
    ext_dir.mkdir()
    (ext_dir / "entry_a.json").write_text("{}")
    (ext_dir / "entry_b.json").write_text("{}")
    (ext_dir / "entry_c.md").write_text("")
    result = list_extraction_files(tmp_path)
    assert sorted(result) == ["entry_a", "entry_b"]


def test_get_cluster_data():
    from multiomics_kg.extraction.cluster.extraction_utils import get_cluster_data
    clusters = {"1": {"name": "Cluster 1"}, "2": {"name": "Cluster 2"}}
    assert get_cluster_data(clusters, "1") == {"name": "Cluster 1"}
    assert get_cluster_data(clusters, 1) == {"name": "Cluster 1"}
    assert get_cluster_data(clusters, "99") == {}
