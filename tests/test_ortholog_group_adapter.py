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
        "product": "DNA polymerase III, beta subunit",
        "gene_name": "dnaN",
        "organism_strain": "Prochlorococcus MED4",
        "ortholog_groups": [
            {"og_id": "cyanorak:CK_00000364", "source": "cyanorak", "taxonomic_level": "curated", "taxon_id": 0, "specificity_rank": 0},
            {"og_id": "eggnog:COG0592@2", "source": "eggnog", "taxonomic_level": "Bacteria", "taxon_id": 2, "specificity_rank": 3},
            {"og_id": "eggnog:1MKTR@1212", "source": "eggnog", "taxonomic_level": "Prochloraceae", "taxon_id": 1212, "specificity_rank": 1},
        ],
    },
    "PMM0002": {
        "locus_tag": "PMM0002",
        "product": "conserved hypothetical protein",
        "gene_name": None,
        "organism_strain": "Prochlorococcus MED4",
        "ortholog_groups": [
            {"og_id": "cyanorak:CK_00000363", "source": "cyanorak", "taxonomic_level": "curated", "taxon_id": 0, "specificity_rank": 0},
            {"og_id": "eggnog:COG0592@2", "source": "eggnog", "taxonomic_level": "Bacteria", "taxon_id": 2, "specificity_rank": 3},
        ],
    },
    "PMM0003": {
        "locus_tag": "PMM0003",
        "product": "hypothetical protein",
        "gene_name": None,
        "organism_strain": "Prochlorococcus MED4",
        "ortholog_groups": [],
    },
}

STRAIN2_DATA = {
    "ALT831_RS00180": {
        "locus_tag": "ALT831_RS00180",
        "product": "DNA polymerase III, beta subunit",
        "gene_name": "dnaN",
        "organism_strain": "Alteromonas macleodii MIT1002",
        "ortholog_groups": [
            {"og_id": "eggnog:COG0592@2", "source": "eggnog", "taxonomic_level": "Bacteria", "taxon_id": 2, "specificity_rank": 3},
            {"og_id": "eggnog:4648R@72275", "source": "eggnog", "taxonomic_level": "Alteromonadaceae", "taxon_id": 72275, "specificity_rank": 1},
        ],
    },
}

STRAIN3_DATA = {
    "A9601_00010": {
        "locus_tag": "A9601_00010",
        "product": "DNA polymerase III, beta subunit",
        "gene_name": "dnaN",
        "organism_strain": "Prochlorococcus AS9601",
        "ortholog_groups": [
            {"og_id": "eggnog:COG0592@2", "source": "eggnog", "taxonomic_level": "Bacteria", "taxon_id": 2, "specificity_rank": 3},
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
def temp_strain3_dir():
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "gene_annotations_merged.json"), "w") as f:
            json.dump(STRAIN3_DATA, f)
        yield d


@pytest.fixture
def multi_config_csv(temp_strain1_dir, temp_strain2_dir, temp_strain3_dir, tmp_path):
    csv_path = tmp_path / "genomes.csv"
    with open(csv_path, "w") as f:
        f.write("ncbi_accession,data_dir,strain_name,ncbi_taxon_id,clade\n")
        f.write(f"GCF_000011465.1,{temp_strain1_dir}/,MED4,59919,HLI\n")
        f.write(f"GCF_901457835.2,{temp_strain2_dir}/,MIT1002,28108,\n")
        f.write(f"GCF_000015849.1,{temp_strain3_dir}/,AS9601,146891,HLII\n")
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
                    {"og_id": f"eggnog:COG{i:04d}@2", "source": "eggnog", "taxonomic_level": "Bacteria", "taxon_id": 2, "specificity_rank": 3},
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
            assert "name" in props
            assert "source" in props
            assert props["source"] in ("cyanorak", "eggnog")
            assert "taxonomic_level" in props
            assert isinstance(props["taxon_id"], int)

    def test_name_is_raw_identifier(self, multi_config_csv):
        """name property should be the raw OG ID without the source prefix."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["cyanorak:CK_00000364"]["name"] == "CK_00000364"
        assert node_map["eggnog:COG0592@2"]["name"] == "COG0592@2"
        assert node_map["eggnog:4648R@72275"]["name"] == "4648R@72275"

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
        # PMM0001 3 + PMM0002 2 + PMM0003 0 + ALT831 2 + A9601 1 = 8
        assert len(edges) == 8

    def test_download_data_is_noop(self, multi_config_csv):
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        adapter.download_data()  # should not raise


# ---------------------------------------------------------------------------
# OrthologGroup enrichment properties
# ---------------------------------------------------------------------------


class TestOrthologGroupEnrichment:
    """Tests for consensus/aggregate properties on OrthologGroup nodes."""

    def test_consensus_product_majority_vote(self, multi_config_csv):
        """COG0592@2 has 3 members with 'DNA polymerase III, beta subunit'
        and 1 with 'conserved hypothetical protein' -> majority wins."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["eggnog:COG0592@2"]["consensus_product"] == "DNA polymerase III, beta subunit"

    def test_consensus_product_excludes_hypothetical(self, multi_config_csv):
        """CK_00000364 has 1 member (PMM0001) with real product -> uses it."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["cyanorak:CK_00000364"]["consensus_product"] == "DNA polymerase III, beta subunit"

    def test_consensus_product_falls_back_to_hypothetical(self, multi_config_csv):
        """CK_00000363 has only PMM0002 ('conserved hypothetical protein') -> falls back."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["cyanorak:CK_00000363"]["consensus_product"] == "conserved hypothetical protein"

    def test_consensus_gene_name(self, multi_config_csv):
        """COG0592@2 has 3 members with gene_name='dnaN', 1 with None -> 'dnaN'."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["eggnog:COG0592@2"]["consensus_gene_name"] == "dnaN"

    def test_consensus_gene_name_none_when_all_null(self, multi_config_csv):
        """CK_00000363 has only PMM0002 (gene_name=None) -> None."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["cyanorak:CK_00000363"]["consensus_gene_name"] is None

    def test_member_count(self, multi_config_csv):
        """COG0592@2 should have 4 members (PMM0001, PMM0002, ALT831_RS00180, A9601_00010)."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["eggnog:COG0592@2"]["member_count"] == 4
        assert node_map["cyanorak:CK_00000364"]["member_count"] == 1

    def test_organism_count(self, multi_config_csv):
        """COG0592@2 spans 3 organisms (MED4, MIT1002, AS9601)."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["eggnog:COG0592@2"]["organism_count"] == 3

    def test_genera_list(self, multi_config_csv):
        """COG0592@2 spans Alteromonas + Prochlorococcus -> sorted list."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["eggnog:COG0592@2"]["genera"] == ["Alteromonas", "Prochlorococcus"]

    def test_has_cross_genus_true(self, multi_config_csv):
        """COG0592@2 has both Prochlorococcus and Alteromonas -> True."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["eggnog:COG0592@2"]["has_cross_genus_members"] == "cross_genus"

    def test_has_cross_genus_false(self, multi_config_csv):
        """CK_00000364 has only PMM0001 (Prochlorococcus) -> False."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["cyanorak:CK_00000364"]["has_cross_genus_members"] == "single_genus"

    def test_string_sanitization(self, tmp_path):
        """Product containing ' or | should be cleaned."""
        data = {
            "GENE001": {
                "locus_tag": "GENE001",
                "product": "5'-nucleotidase | phosphatase",
                "gene_name": "ush'A",
                "organism_strain": "Prochlorococcus MED4",
                "ortholog_groups": [
                    {"og_id": "test:OG1", "source": "eggnog", "taxonomic_level": "Bacteria", "taxon_id": 2, "specificity_rank": 3},
                ],
            },
        }
        strain_dir = tmp_path / "strain_sanitize"
        strain_dir.mkdir()
        with open(strain_dir / "gene_annotations_merged.json", "w") as f:
            json.dump(data, f)
        csv_path = tmp_path / "genomes_sanitize.csv"
        with open(csv_path, "w") as f:
            f.write("ncbi_accession,data_dir,strain_name,ncbi_taxon_id,clade\n")
            f.write(f"GCF_TEST,{strain_dir}/,TEST,12345,\n")
        adapter = MultiOrthologGroupAdapter(genome_config_file=str(csv_path))
        nodes = adapter.get_nodes()
        props = nodes[0][2]
        assert "'" not in props["consensus_product"]
        assert "|" not in props["consensus_product"]
        assert "'" not in props["consensus_gene_name"]
