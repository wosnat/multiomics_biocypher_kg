"""
Unit tests for multiomics_kg/utils/kegg_utils.py and the KEGG section of
multiomics_kg/adapters/functional_annotation_adapter.py.

All tests are offline — no real network calls, no live Neo4j required.
"""

import json
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from multiomics_kg.utils.kegg_utils import (
    _parse_ko_names,
    _parse_pathway_names,
    _parse_ko_to_pathways,
    _parse_brite_hierarchy,
    download_kegg_data,
)
from multiomics_kg.adapters.functional_annotation_adapter import (
    KeggAnnotationAdapter,
    MultiKeggAnnotationAdapter,
    _ko_node_id,
    _pathway_node_id,
    _subcat_node_id,
    _cat_node_id,
)


# ---------------------------------------------------------------------------
# Shared minimal test data
# ---------------------------------------------------------------------------

MINI_KEGG_DATA = {
    "ko_names": {
        "K02338": "DNA-directed DNA polymerase; dnaN",
        "K01952": "phosphoribosylformylglycinamidine synthase; purL",
    },
    "pathway_names": {
        "ko03030": "DNA replication",
        "ko00230": "Purine metabolism",
    },
    "ko_to_pathways": {
        "K02338": ["ko03030"],
        "K01952": ["ko00230"],
    },
    "pathway_to_subcategory": {
        "ko03030": "09124",
        "ko00230": "09102",
    },
    "subcategory_names": {
        "09124": "Replication and repair",
        "09102": "Energy metabolism",
    },
    "subcategory_to_category": {
        "09124": "09100",
        "09102": "09100",
    },
    "category_names": {
        "09100": "Metabolism",
    },
}

MINI_GENE_DATA = {
    "PMM0001": {
        "locus_tag": "PMM0001",
        "kegg_ko": ["K02338"],
    },
    "PMM0002": {
        "locus_tag": "PMM0002",
        "kegg_ko": ["K01952"],
    },
    "PMM0003": {
        "locus_tag": "PMM0003",
        "kegg_ko": [],  # no KO annotations
    },
    "PMM0004": {
        "locus_tag": "PMM0004",
        "kegg_ko": None,  # null field
    },
}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def genome_dir(tmp_path):
    """Write MINI_GENE_DATA as gene_annotations_merged.json, return dir path."""
    d = tmp_path / "MED4"
    d.mkdir()
    (d / "gene_annotations_merged.json").write_text(
        json.dumps(MINI_GENE_DATA), encoding="utf-8"
    )
    return d


@pytest.fixture
def kegg_cache(tmp_path):
    """Write MINI_KEGG_DATA as kegg_data.json, return cache_root."""
    kegg_dir = tmp_path / "kegg"
    kegg_dir.mkdir()
    (kegg_dir / "kegg_data.json").write_text(
        json.dumps(MINI_KEGG_DATA), encoding="utf-8"
    )
    return tmp_path  # cache_root


# ---------------------------------------------------------------------------
# kegg_utils: parser unit tests
# ---------------------------------------------------------------------------


def test_parse_ko_names():
    text = "ko:K02338\tDNA polymerase\nko:K01952\tpurL\n"
    result = _parse_ko_names(text)
    assert result == {"K02338": "DNA polymerase", "K01952": "purL"}


def test_parse_ko_names_skips_non_k():
    text = "ko:K02338\tDNA polymerase\nko:M00001\tsomething\n"
    result = _parse_ko_names(text)
    assert "M00001" not in result
    assert "K02338" in result


def test_parse_pathway_names():
    text = "path:ko03030\tDNA replication\npath:map03030\tDNA replication (map)\n"
    result = _parse_pathway_names(text)
    assert "ko03030" in result
    assert result["ko03030"] == "DNA replication"
    assert "map03030" not in result  # only ko-prefixed


def test_parse_ko_to_pathways():
    text = "ko:K02338\tpath:ko03030\nko:K02338\tpath:ko00230\nko:K01952\tpath:ko00230\n"
    result = _parse_ko_to_pathways(text)
    assert set(result["K02338"]) == {"ko03030", "ko00230"}
    assert result["K01952"] == ["ko00230"]


def test_parse_ko_to_pathways_skips_map():
    text = "ko:K02338\tpath:map03030\nko:K02338\tpath:ko03030\n"
    result = _parse_ko_to_pathways(text)
    # map-prefixed pathways should be skipped (only ko-prefix retained)
    assert all(pw.startswith("ko") for pw in result["K02338"])


def test_parse_brite_hierarchy():
    brite_json = {
        "name": "ko00001",
        "children": [
            {
                "name": "09100 Metabolism",
                "children": [
                    {
                        "name": "09102 Energy metabolism",
                        "children": [
                            {"name": "00230 Purine metabolism [PATH:ko00230]", "children": []},
                        ],
                    }
                ],
            },
            {
                "name": "09120 Genetic Information Processing",
                "children": [
                    {
                        "name": "09124 Replication and repair",
                        "children": [
                            {"name": "03030 DNA replication [PATH:ko03030]", "children": []},
                        ],
                    }
                ],
            },
        ],
    }
    pw_to_sc, sc_names, sc_to_cat, cat_names = _parse_brite_hierarchy(brite_json)

    assert pw_to_sc["ko00230"] == "09102"
    assert pw_to_sc["ko03030"] == "09124"
    assert sc_names["09102"] == "Energy metabolism"
    assert sc_names["09124"] == "Replication and repair"
    assert sc_to_cat["09102"] == "09100"
    assert sc_to_cat["09124"] == "09120"
    assert cat_names["09100"] == "Metabolism"
    assert cat_names["09120"] == "Genetic Information Processing"


def test_download_kegg_data_uses_cache(kegg_cache):
    """If cache exists, no network calls are made."""
    with patch("multiomics_kg.utils.kegg_utils.requests.get") as mock_get:
        result = download_kegg_data(kegg_cache)
        mock_get.assert_not_called()
    assert result["ko_names"]["K02338"] == "DNA-directed DNA polymerase; dnaN"


def test_download_kegg_data_writes_cache(tmp_path):
    """With no cache, data is fetched and written."""
    mock_resp = MagicMock()
    mock_resp.raise_for_status = MagicMock()
    mock_resp.text = "ko:K02338\tDNA polymerase\n"
    mock_resp.json.return_value = {"name": "ko00001", "children": []}

    with patch("multiomics_kg.utils.kegg_utils.requests.get", return_value=mock_resp):
        result = download_kegg_data(tmp_path, force=False)

    cache_file = tmp_path / "kegg" / "kegg_data.json"
    assert cache_file.exists()
    assert "ko_names" in result


def test_download_kegg_data_force_refetch(kegg_cache):
    """force=True re-downloads even if cache exists."""
    mock_resp = MagicMock()
    mock_resp.raise_for_status = MagicMock()
    mock_resp.text = "ko:K99999\tNew KO\n"
    mock_resp.json.return_value = {"name": "ko00001", "children": []}

    with patch("multiomics_kg.utils.kegg_utils.requests.get", return_value=mock_resp):
        result = download_kegg_data(kegg_cache, force=True)

    assert "K99999" in result["ko_names"]


# ---------------------------------------------------------------------------
# Node ID helpers
# ---------------------------------------------------------------------------


def test_ko_node_id():
    nid = _ko_node_id("K02338")
    assert "K02338" in nid


def test_pathway_node_id():
    nid = _pathway_node_id("ko03030")
    assert "ko03030" in nid


def test_subcat_node_id():
    assert _subcat_node_id("09124") == "kegg.subcategory:09124"


def test_cat_node_id():
    assert _cat_node_id("09100") == "kegg.category:09100"


# ---------------------------------------------------------------------------
# KeggAnnotationAdapter (per-strain)
# ---------------------------------------------------------------------------


def test_kegg_adapter_get_all_ko_ids(genome_dir):
    adapter = KeggAnnotationAdapter(genome_dir)
    ids = adapter.get_all_ko_ids()
    assert "K02338" in ids
    assert "K01952" in ids


def test_kegg_adapter_skips_empty_and_none(genome_dir):
    adapter = KeggAnnotationAdapter(genome_dir)
    ids = adapter.get_all_ko_ids()
    assert None not in ids
    assert "" not in ids


def test_kegg_adapter_get_edges(genome_dir):
    adapter = KeggAnnotationAdapter(genome_dir)
    edges = list(adapter.get_edges())
    labels = {e[3] for e in edges}
    assert labels == {"gene_has_kegg_ko"}
    # PMM0001 → K02338
    sources = {e[1] for e in edges}
    assert any("PMM0001" in s for s in sources)
    targets = {e[2] for e in edges}
    assert any("K02338" in t for t in targets)


def test_kegg_adapter_test_mode(genome_dir):
    adapter = KeggAnnotationAdapter(genome_dir, test_mode=True)
    edges = list(adapter.get_edges())
    assert len(edges) <= 100


def test_kegg_adapter_missing_file(tmp_path):
    """Missing JSON file should produce no edges (with warning, no crash)."""
    empty_dir = tmp_path / "empty_strain"
    empty_dir.mkdir()
    adapter = KeggAnnotationAdapter(empty_dir)
    assert list(adapter.get_edges()) == []


# ---------------------------------------------------------------------------
# MultiKeggAnnotationAdapter
# ---------------------------------------------------------------------------


@pytest.fixture
def genome_config_file(tmp_path, genome_dir):
    """Write a minimal cyanobacteria_genomes.csv pointing to genome_dir."""
    csv_path = tmp_path / "genomes.csv"
    csv_path.write_text(f"organism,strain,data_dir\nPro,MED4,{genome_dir}\n", encoding="utf-8")
    return str(csv_path)


@pytest.fixture
def multi_kegg(genome_config_file, kegg_cache):
    """MultiKeggAnnotationAdapter loaded from cache (no network)."""
    return MultiKeggAnnotationAdapter(
        genome_config_file=genome_config_file,
        cache_root=kegg_cache,
        test_mode=False,
        cache=True,
    )


def test_multi_kegg_get_nodes_types(multi_kegg):
    nodes = list(multi_kegg.get_nodes())
    node_labels = {n[1] for n in nodes}
    assert "kegg orthologous group" in node_labels
    assert "kegg pathway" in node_labels
    assert "kegg subcategory" in node_labels
    assert "kegg category" in node_labels


def test_multi_kegg_ko_nodes_have_names(multi_kegg):
    nodes = list(multi_kegg.get_nodes())
    ko_nodes = [n for n in nodes if n[1] == "kegg orthologous group"]
    for node_id, label, props in ko_nodes:
        assert "name" in props


def test_multi_kegg_nodes_deduplicated(tmp_path, genome_dir, kegg_cache):
    """Two strains sharing the same KO should yield only one KO node."""
    # Add a second strain with the same gene data
    strain2 = tmp_path / "MIT9312"
    strain2.mkdir()
    (strain2 / "gene_annotations_merged.json").write_text(
        json.dumps(MINI_GENE_DATA), encoding="utf-8"
    )
    csv_path = tmp_path / "genomes.csv"
    csv_path.write_text(
        f"organism,strain,data_dir\nPro,MED4,{genome_dir}\nPro,MIT9312,{strain2}\n",
        encoding="utf-8",
    )
    adapter = MultiKeggAnnotationAdapter(
        genome_config_file=str(csv_path),
        cache_root=kegg_cache,
        cache=True,
    )
    nodes = list(adapter.get_nodes())
    ko_nodes = [n for n in nodes if n[1] == "kegg orthologous group"]
    ko_ids = [n[0] for n in ko_nodes]
    assert len(ko_ids) == len(set(ko_ids)), "Duplicate KO nodes found"


def test_multi_kegg_get_edges_labels(multi_kegg):
    edges = list(multi_kegg.get_edges())
    labels = {e[3] for e in edges}
    assert "gene_has_kegg_ko" in labels
    assert "ko_in_kegg_pathway" in labels
    assert "kegg_pathway_in_kegg_subcategory" in labels
    assert "kegg_subcategory_in_kegg_category" in labels


def test_multi_kegg_ko_pathway_edges(multi_kegg):
    edges = list(multi_kegg.get_edges())
    ko_pw_edges = [e for e in edges if e[3] == "ko_in_kegg_pathway"]
    # K02338 → ko03030 should exist
    assert any("K02338" in e[1] and "ko03030" in e[2] for e in ko_pw_edges)


def test_multi_kegg_pathway_subcat_edges(multi_kegg):
    edges = list(multi_kegg.get_edges())
    pw_sc_edges = [e for e in edges if e[3] == "kegg_pathway_in_kegg_subcategory"]
    assert len(pw_sc_edges) > 0
    # Targets should be subcategory node IDs
    for e in pw_sc_edges:
        assert "kegg.subcategory:" in e[2]


def test_multi_kegg_subcat_cat_edges(multi_kegg):
    edges = list(multi_kegg.get_edges())
    sc_cat_edges = [e for e in edges if e[3] == "kegg_subcategory_in_kegg_category"]
    assert len(sc_cat_edges) > 0
    for e in sc_cat_edges:
        assert "kegg.category:" in e[2]


def test_multi_kegg_clean_str_in_names(multi_kegg):
    """Node names must not contain raw single-quotes or pipe characters."""
    nodes = list(multi_kegg.get_nodes())
    for node_id, label, props in nodes:
        name = props.get("name", "")
        assert "'" not in name, f"Raw single-quote in node name: {name!r}"
        assert "|" not in name, f"Pipe in node name: {name!r}"
