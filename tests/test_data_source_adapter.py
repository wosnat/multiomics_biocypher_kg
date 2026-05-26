"""Unit tests for DataSourceAdapter.

Reads logical_sources from gene_annotations_config.yaml; emits one
DataSource node per logical source with auto-derived info_types.
"""

import pytest
from multiomics_kg.adapters.data_source_adapter import DataSourceAdapter


@pytest.fixture
def adapter():
    return DataSourceAdapter(config_path="config/gene_annotations_config.yaml")


def test_emits_five_nodes(adapter):
    """Five data sources: ncbi, cyanorak, uniprot, eggnog, psortb."""
    adapter.download_data()
    nodes = list(adapter.get_nodes())
    ids = {props["id"] for _, _, props in nodes}
    assert ids == {"ncbi", "cyanorak", "uniprot", "eggnog", "psortb"}


def test_psortb_is_tool_run_gene_level(adapter):
    adapter.download_data()
    nodes = {props["id"]: props for _, _, props in adapter.get_nodes()}
    psortb = nodes["psortb"]
    assert psortb["provenance"] == "tool_run"
    assert psortb["scope"] == "gene_level"
    assert psortb["name"] == "PSORTb"
    # info_types auto-derived from the two psortb field rules
    assert {"psortb_localization", "psortb_score"} <= set(psortb["info_types"])


def test_ncbi_node_properties(adapter):
    adapter.download_data()
    nodes = {props["id"]: props for _, _, props in adapter.get_nodes()}
    ncbi = nodes["ncbi"]
    assert ncbi["scope"] == "gene_level"
    assert ncbi["provenance"] == "download"
    assert ncbi["applies_to_organisms"] == []


def test_cyanorak_organism_restricted(adapter):
    adapter.download_data()
    nodes = {props["id"]: props for _, _, props in adapter.get_nodes()}
    cyano = nodes["cyanorak"]
    assert cyano["scope"] == "organism_restricted"
    assert cyano["provenance"] == "download"
    assert any("Prochlorococcus" in o for o in cyano["applies_to_organisms"])


def test_eggnog_is_tool_run(adapter):
    adapter.download_data()
    nodes = {props["id"]: props for _, _, props in adapter.get_nodes()}
    assert nodes["eggnog"]["provenance"] == "tool_run"


def test_info_types_auto_derived(adapter):
    """info_types is computed from the fields block; eggnog should include
    cog_category, kegg_ko, go_terms, etc."""
    adapter.download_data()
    nodes = {props["id"]: props for _, _, props in adapter.get_nodes()}
    eggnog_info = set(nodes["eggnog"]["info_types"])
    # spot-check a few expected entries from the fields block
    assert "kegg_ko" in eggnog_info
    assert "eggnog_ogs" in eggnog_info


def test_emits_no_edges(adapter):
    adapter.download_data()
    edges = list(adapter.get_edges())
    assert edges == []
