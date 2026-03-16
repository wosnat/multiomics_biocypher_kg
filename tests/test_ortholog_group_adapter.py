"""Tests for OrthologGroupAdapter and MultiOrthologGroupAdapter."""

import json
import os
import tempfile

import pytest

from multiomics_kg.adapters.ortholog_group_adapter import (
    MultiOrthologGroupAdapter,
    OrthologGroupAdapter,
)


# ---------------------------------------------------------------------------
# Test data
# ---------------------------------------------------------------------------

STRAIN1_DATA = {
    "PMM0001": {
        "locus_tag": "PMM0001",
        "ortholog_groups": [
            {"og_id": "cyanorak:CK_00000364", "source": "cyanorak", "taxonomic_level": "curated", "taxon_id": 0},
            {"og_id": "eggnog:COG0592@2", "source": "eggnog", "taxonomic_level": "Bacteria", "taxon_id": 2},
            {"og_id": "eggnog:1MKTR@1212", "source": "eggnog", "taxonomic_level": "Prochloraceae", "taxon_id": 1212},
        ],
    },
    "PMM0002": {
        "locus_tag": "PMM0002",
        "ortholog_groups": [
            {"og_id": "cyanorak:CK_00000363", "source": "cyanorak", "taxonomic_level": "curated", "taxon_id": 0},
            {"og_id": "eggnog:COG0592@2", "source": "eggnog", "taxonomic_level": "Bacteria", "taxon_id": 2},
        ],
    },
    "PMM0003": {
        "locus_tag": "PMM0003",
        "ortholog_groups": [],
    },
}

STRAIN2_DATA = {
    "ALT831_RS00180": {
        "locus_tag": "ALT831_RS00180",
        "ortholog_groups": [
            {"og_id": "eggnog:COG0592@2", "source": "eggnog", "taxonomic_level": "Bacteria", "taxon_id": 2},
            {"og_id": "eggnog:4648R@72275", "source": "eggnog", "taxonomic_level": "Alteromonadaceae", "taxon_id": 72275},
        ],
    },
}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def temp_strain1_dir():
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "gene_annotations_merged.json"), "w") as f:
            json.dump(STRAIN1_DATA, f)
        yield d


@pytest.fixture
def temp_strain2_dir():
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "gene_annotations_merged.json"), "w") as f:
            json.dump(STRAIN2_DATA, f)
        yield d


@pytest.fixture
def multi_config_csv(temp_strain1_dir, temp_strain2_dir, tmp_path):
    csv_path = tmp_path / "genomes.csv"
    with open(csv_path, "w") as f:
        f.write("ncbi_accession,data_dir,strain_name,ncbi_taxon_id,clade\n")
        f.write(f"GCF_000011465.1,{temp_strain1_dir}/,MED4,59919,HLI\n")
        f.write(f"GCF_901457835.2,{temp_strain2_dir}/,MIT1002,28108,\n")
    return str(csv_path)


# ---------------------------------------------------------------------------
# OrthologGroupAdapter (single strain)
# ---------------------------------------------------------------------------


class TestOrthologGroupAdapter:
    def test_loads_memberships(self, temp_strain1_dir):
        adapter = OrthologGroupAdapter(genome_dir=temp_strain1_dir)
        memberships = adapter.get_og_memberships()
        # PMM0001 has 3, PMM0002 has 2, PMM0003 has 0 = total 5
        assert len(memberships) == 5

    def test_membership_structure(self, temp_strain1_dir):
        adapter = OrthologGroupAdapter(genome_dir=temp_strain1_dir)
        for lt, og in adapter.get_og_memberships():
            assert isinstance(lt, str)
            assert "og_id" in og
            assert "source" in og
            assert "taxonomic_level" in og
            assert "taxon_id" in og

    def test_test_mode_limits_results(self, tmp_path):
        # Create data with >100 genes
        data = {}
        for i in range(200):
            lt = f"GENE_{i:04d}"
            data[lt] = {
                "locus_tag": lt,
                "ortholog_groups": [
                    {"og_id": f"eggnog:COG{i:04d}@2", "source": "eggnog", "taxonomic_level": "Bacteria", "taxon_id": 2},
                ],
            }
        with open(tmp_path / "gene_annotations_merged.json", "w") as f:
            json.dump(data, f)
        adapter = OrthologGroupAdapter(genome_dir=tmp_path, test_mode=True)
        memberships = adapter.get_og_memberships()
        assert len(memberships) == 100

    def test_missing_json_returns_empty(self, tmp_path):
        adapter = OrthologGroupAdapter(genome_dir=tmp_path / "nonexistent")
        assert adapter.get_og_memberships() == []


# ---------------------------------------------------------------------------
# MultiOrthologGroupAdapter
# ---------------------------------------------------------------------------


class TestMultiOrthologGroupAdapter:
    def test_yields_nodes_and_edges(self, multi_config_csv):
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        nodes = adapter.get_nodes()
        edges = adapter.get_edges()
        assert len(nodes) > 0
        assert len(edges) > 0

    def test_deduplicates_nodes(self, multi_config_csv):
        """Same OG (COG0592@2) from two strains → one node."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        nodes = adapter.get_nodes()
        og_ids = [n[0] for n in nodes]
        assert og_ids.count("eggnog:COG0592@2") == 1

    def test_node_properties(self, multi_config_csv):
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        nodes = adapter.get_nodes()
        for og_id, label, props in nodes:
            assert label == "ortholog_group"
            assert "source" in props
            assert props["source"] in ("cyanorak", "eggnog")
            assert "taxonomic_level" in props
            assert isinstance(props["taxon_id"], int)

    def test_edge_ids_unique(self, multi_config_csv):
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        edges = adapter.get_edges()
        edge_ids = [e[0] for e in edges]
        assert len(edge_ids) == len(set(edge_ids))

    def test_edge_structure(self, multi_config_csv):
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        for edge_id, source, target, edge_type, props in adapter.get_edges():
            assert edge_type == "gene_in_ortholog_group"
            assert source.startswith("ncbigene:")
            assert isinstance(props, dict)

    def test_node_count(self, multi_config_csv):
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        nodes = adapter.get_nodes()
        # Unique OGs: CK_00000364, CK_00000363, COG0592@2, 1MKTR@1212, 4648R@72275 = 5
        assert len(nodes) == 5

    def test_edge_count(self, multi_config_csv):
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        edges = adapter.get_edges()
        # PMM0001 3 + PMM0002 2 + PMM0003 0 + ALT831 2 = 7
        assert len(edges) == 7

    def test_download_data_is_noop(self, multi_config_csv):
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        adapter.download_data()  # should not raise
