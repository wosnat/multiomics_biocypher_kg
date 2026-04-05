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
            "source_figures": ["Figure 2"],
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
