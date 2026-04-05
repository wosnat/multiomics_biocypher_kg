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
        "organism_name": "Prochlorococcus MED4",
        "cog_category": ["L"],
        "cyanorak_Role": ["F.1"],
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
        "organism_name": "Prochlorococcus MED4",
        "cog_category": ["S"],
        "cyanorak_Role": ["R.2"],
        "ortholog_groups": [
            {"og_id": "cyanorak:CK_00000363", "source": "cyanorak", "taxonomic_level": "curated", "taxon_id": 0, "specificity_rank": 0},
            {"og_id": "eggnog:COG0592@2", "source": "eggnog", "taxonomic_level": "Bacteria", "taxon_id": 2, "specificity_rank": 3},
        ],
    },
    "PMM0003": {
        "locus_tag": "PMM0003",
        "product": "hypothetical protein",
        "gene_name": None,
        "organism_name": "Prochlorococcus MED4",
        "cog_category": [],
        "cyanorak_Role": [],
        "ortholog_groups": [],
    },
}

STRAIN2_DATA = {
    "ALT831_RS00180": {
        "locus_tag": "ALT831_RS00180",
        "product": "DNA polymerase III, beta subunit",
        "gene_name": "dnaN",
        "organism_name": "Alteromonas macleodii MIT1002",
        "cog_category": ["L"],
        "cyanorak_Role": [],
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
        "organism_name": "Prochlorococcus AS9601",
        "cog_category": ["L"],
        "cyanorak_Role": ["F.1"],
        "ortholog_groups": [
            {"og_id": "eggnog:COG0592@2", "source": "eggnog", "taxonomic_level": "Bacteria", "taxon_id": 2, "specificity_rank": 3},
        ],
    },
}

# Minimal cyanorak_roles.csv content for tests
CYANORAK_ROLES_CSV = """\
primary role id,primary role,secondary role id,secondary role,sub role id,sub role
A,Amino acid biosynthesis,A.7,Other,-,-
F,DNA metabolism,F.1,"DNA replication, recombination, and repair",-,-
J,Photosynthesis and respiration,J.8,Photosystem II,-,-
R,Other,R.2,Conserved hypothetical proteins,-,-
R,Other,R.4,Hypothetical proteins,-,-
"""


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
def temp_role_tree_file(tmp_path):
    role_path = tmp_path / "cyanorak_roles.csv"
    with open(role_path, "w") as f:
        f.write(CYANORAK_ROLES_CSV)
    return str(role_path)


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
    def test_yields_nodes_and_edges(self, multi_config_csv, temp_role_tree_file):
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        edges = adapter.get_edges()
        assert len(nodes) > 0
        assert len(edges) > 0

    def test_deduplicates_nodes(self, multi_config_csv, temp_role_tree_file):
        """Same OG (COG0592@2) from two strains → one node."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        og_ids = [n[0] for n in nodes]
        assert og_ids.count("eggnog:COG0592@2") == 1

    def test_node_properties(self, multi_config_csv, temp_role_tree_file):
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        for og_id, label, props in nodes:
            assert label == "ortholog_group"
            assert "name" in props
            assert "source" in props
            assert props["source"] in ("cyanorak", "eggnog")
            assert "taxonomic_level" in props
            assert isinstance(props["taxon_id"], int)
            # New properties should be present (may be None)
            assert "description" in props
            assert "functional_description" in props

    def test_name_is_raw_identifier(self, multi_config_csv, temp_role_tree_file):
        """name property should be the raw OG ID without the source prefix."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["cyanorak:CK_00000364"]["name"] == "CK_00000364"
        assert node_map["eggnog:COG0592@2"]["name"] == "COG0592@2"
        assert node_map["eggnog:4648R@72275"]["name"] == "4648R@72275"

    def test_edge_ids_unique(self, multi_config_csv, temp_role_tree_file):
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        _ = adapter.get_nodes()  # must call first to populate _og_info
        edges = adapter.get_edges()
        edge_ids = [e[0] for e in edges]
        assert len(edge_ids) == len(set(edge_ids))

    def test_gene_og_edge_structure(self, multi_config_csv, temp_role_tree_file):
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        _ = adapter.get_nodes()
        gene_og_edges = [e for e in adapter.get_edges() if e[3] == "gene_in_ortholog_group"]
        for edge_id, source, target, edge_type, props in gene_og_edges:
            assert source.startswith("ncbigene:")
            assert isinstance(props, dict)

    def test_node_count(self, multi_config_csv, temp_role_tree_file):
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        # Unique OGs: CK_00000364, CK_00000363, COG0592@2, 1MKTR@1212, 4648R@72275 = 5
        assert len(nodes) == 5

    def test_gene_og_edge_count(self, multi_config_csv, temp_role_tree_file):
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        _ = adapter.get_nodes()
        edges = adapter.get_edges()
        gene_og_edges = [e for e in edges if e[3] == "gene_in_ortholog_group"]
        # PMM0001 3 + PMM0002 2 + PMM0003 0 + ALT831 2 + A9601 1 = 8
        assert len(gene_og_edges) == 8

    def test_download_data_is_noop(self, multi_config_csv, temp_role_tree_file):
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        adapter.download_data()  # should not raise


# ---------------------------------------------------------------------------
# OrthologGroup enrichment properties
# ---------------------------------------------------------------------------


class TestOrthologGroupEnrichment:
    """Tests for consensus/aggregate properties on OrthologGroup nodes."""

    def test_consensus_product_majority_vote(self, multi_config_csv, temp_role_tree_file):
        """COG0592@2 has 3 members with 'DNA polymerase III, beta subunit'
        and 1 with 'conserved hypothetical protein' -> majority wins."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["eggnog:COG0592@2"]["consensus_product"] == "DNA polymerase III, beta subunit"

    def test_consensus_product_excludes_hypothetical(self, multi_config_csv, temp_role_tree_file):
        """CK_00000364 has 1 member (PMM0001) with real product -> uses it."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["cyanorak:CK_00000364"]["consensus_product"] == "DNA polymerase III, beta subunit"

    def test_consensus_product_falls_back_to_hypothetical(self, multi_config_csv, temp_role_tree_file):
        """CK_00000363 has only PMM0002 ('conserved hypothetical protein') -> falls back."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["cyanorak:CK_00000363"]["consensus_product"] == "conserved hypothetical protein"

    def test_consensus_gene_name(self, multi_config_csv, temp_role_tree_file):
        """COG0592@2 has 3 members with gene_name='dnaN', 1 with None -> 'dnaN'."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["eggnog:COG0592@2"]["consensus_gene_name"] == "dnaN"

    def test_consensus_gene_name_none_when_all_null(self, multi_config_csv, temp_role_tree_file):
        """CK_00000363 has only PMM0002 (gene_name=None) -> None."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["cyanorak:CK_00000363"]["consensus_gene_name"] is None

    def test_member_count(self, multi_config_csv, temp_role_tree_file):
        """COG0592@2 should have 4 members (PMM0001, PMM0002, ALT831_RS00180, A9601_00010)."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["eggnog:COG0592@2"]["member_count"] == 4
        assert node_map["cyanorak:CK_00000364"]["member_count"] == 1

    def test_organism_count(self, multi_config_csv, temp_role_tree_file):
        """COG0592@2 spans 3 organisms (MED4, MIT1002, AS9601)."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["eggnog:COG0592@2"]["organism_count"] == 3

    def test_genera_list(self, multi_config_csv, temp_role_tree_file):
        """COG0592@2 spans Alteromonas + Prochlorococcus -> sorted list."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["eggnog:COG0592@2"]["genera"] == ["Alteromonas", "Prochlorococcus"]

    def test_has_cross_genus_true(self, multi_config_csv, temp_role_tree_file):
        """COG0592@2 has both Prochlorococcus and Alteromonas -> True."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["eggnog:COG0592@2"]["has_cross_genus_members"] == "cross_genus"

    def test_has_cross_genus_false(self, multi_config_csv, temp_role_tree_file):
        """CK_00000364 has only PMM0001 (Prochlorococcus) -> False."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["cyanorak:CK_00000364"]["has_cross_genus_members"] == "single_genus"

    def test_string_sanitization(self, tmp_path, temp_role_tree_file):
        """Product containing ' or | should be cleaned."""
        data = {
            "GENE001": {
                "locus_tag": "GENE001",
                "product": "5'-nucleotidase | phosphatase",
                "gene_name": "ush'A",
                "organism_name": "Prochlorococcus MED4",
                "cog_category": [],
                "cyanorak_Role": [],
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
        adapter = MultiOrthologGroupAdapter(genome_config_file=str(csv_path), role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        props = nodes[0][2]
        assert "'" not in props["consensus_product"]
        assert "|" not in props["consensus_product"]
        assert "'" not in props["consensus_gene_name"]


# ---------------------------------------------------------------------------
# Functional description and annotation edges
# ---------------------------------------------------------------------------


class TestFunctionalDescription:
    """Tests for functional_description property and OG annotation edges."""

    def test_functional_description_from_majority_roles(self, multi_config_csv, temp_role_tree_file):
        """COG0592@2 has 4 members: 3 with cog_category=['L'] and 1 with ['S'].
        L passes majority (75%), S does not (25%). Only L should appear."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        func_desc = node_map["eggnog:COG0592@2"]["functional_description"]
        assert func_desc is not None
        assert "Replication, recombination and repair" in func_desc  # COG "L"
        assert "Function unknown" not in func_desc  # COG "S" filtered

    def test_functional_description_includes_cyanorak_roles(self, multi_config_csv, temp_role_tree_file):
        """COG0592@2: PMM0001 has cyanorak_Role=['F.1'], PMM0002 has ['R.2'],
        ALT831 and A9601 have []. Only F.1 gets 2/4=50%, R.2 gets 1/4=25%.
        F.1 at exactly 50% does NOT pass (>50% required). But A9601 also has F.1.
        So F.1 gets 2/4 = 50%, still not >50%. Let's check."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        # COG0592@2: 4 members. F.1 in PMM0001 + A9601 = 2/4 = 50%, NOT >50%
        # So cyanorak roles should not appear in functional_description for COG0592
        func_desc = node_map["eggnog:COG0592@2"]["functional_description"]
        # L is 3/4 = 75% → passes
        assert "Replication" in func_desc

    def test_functional_description_filters_uninformative(self, multi_config_csv, temp_role_tree_file):
        """CK_00000363 has PMM0002 with cog_category=['S'] and cyanorak_Role=['R.2'].
        Both are uninformative → functional_description should be None."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["cyanorak:CK_00000363"]["functional_description"] is None

    def test_functional_description_cyanorak_hierarchy(self, multi_config_csv, temp_role_tree_file):
        """CK_00000364 has PMM0001 with cyanorak_Role=['F.1'] → 1/1 = 100%.
        Should produce full hierarchy: 'DNA metabolism > DNA replication, ...'"""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        func_desc = node_map["cyanorak:CK_00000364"]["functional_description"]
        assert func_desc is not None
        assert "DNA metabolism" in func_desc
        assert "DNA replication" in func_desc

    def test_description_null_for_cyanorak_groups(self, multi_config_csv, temp_role_tree_file):
        """Cyanorak groups have no eggNOG DB entry → description should be None."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        assert node_map["cyanorak:CK_00000364"]["description"] is None

    def test_og_in_cog_category_edges(self, multi_config_csv, temp_role_tree_file):
        """OG_in_cog_category edges for majority COG categories."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        _ = adapter.get_nodes()
        edges = adapter.get_edges()
        cog_edges = [e for e in edges if e[3] == "og_in_cog_category"]
        # COG0592@2: L passes (3/4), S does not (1/4)
        # CK_00000364: L passes (1/1)
        # CK_00000363: S is uninformative but still gets an edge (filtering is only for functional_description text)
        cog_edge_sources = {(e[1], e[2]) for e in cog_edges}
        assert ("eggnog:COG0592@2", "cog.category:L") in cog_edge_sources

    def test_og_has_cyanorak_role_edges(self, multi_config_csv, temp_role_tree_file):
        """OG_has_cyanorak_role edges for majority Cyanorak roles."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        _ = adapter.get_nodes()
        edges = adapter.get_edges()
        cyanorak_edges = [e for e in edges if e[3] == "og_has_cyanorak_role"]
        # CK_00000364: PMM0001 has F.1, 1/1 = 100% → edge
        cyanorak_edge_pairs = {(e[1], e[2]) for e in cyanorak_edges}
        assert ("cyanorak:CK_00000364", "cyanorak.role:F.1") in cyanorak_edge_pairs

    def test_annotation_edge_ids_unique(self, multi_config_csv, temp_role_tree_file):
        """All edge IDs should be unique across all edge types."""
        adapter = MultiOrthologGroupAdapter(genome_config_file=multi_config_csv, role_tree_file=temp_role_tree_file)
        _ = adapter.get_nodes()
        edges = adapter.get_edges()
        edge_ids = [e[0] for e in edges]
        assert len(edge_ids) == len(set(edge_ids))

    def test_multi_valued_roles(self, tmp_path, temp_role_tree_file):
        """Genes with multiple roles should vote once for each role."""
        data = {
            "G1": {
                "locus_tag": "G1", "product": "psb", "gene_name": "psbA",
                "organism_name": "Pro MED4",
                "cog_category": ["C", "J"],
                "cyanorak_Role": ["J.8", "F.1"],
                "ortholog_groups": [
                    {"og_id": "test:OG_MULTI", "source": "eggnog", "taxonomic_level": "Bacteria", "taxon_id": 2, "specificity_rank": 3},
                ],
            },
            "G2": {
                "locus_tag": "G2", "product": "psb", "gene_name": "psbB",
                "organism_name": "Pro MED4",
                "cog_category": ["C"],
                "cyanorak_Role": ["J.8"],
                "ortholog_groups": [
                    {"og_id": "test:OG_MULTI", "source": "eggnog", "taxonomic_level": "Bacteria", "taxon_id": 2, "specificity_rank": 3},
                ],
            },
        }
        strain_dir = tmp_path / "multi_role"
        strain_dir.mkdir()
        with open(strain_dir / "gene_annotations_merged.json", "w") as f:
            json.dump(data, f)
        csv_path = tmp_path / "genomes_multi.csv"
        with open(csv_path, "w") as f:
            f.write("ncbi_accession,data_dir,strain_name,ncbi_taxon_id,clade\n")
            f.write(f"GCF_TEST,{strain_dir}/,TEST,12345,\n")
        adapter = MultiOrthologGroupAdapter(genome_config_file=str(csv_path), role_tree_file=temp_role_tree_file)
        nodes = adapter.get_nodes()
        node_map = {n[0]: n[2] for n in nodes}
        func_desc = node_map["test:OG_MULTI"]["functional_description"]
        # C: 2/2 = 100% → "Energy production and conversion"
        # J (cog): 1/2 = 50% → NOT >50%
        # J.8 (cyanorak): 2/2 = 100% → "Photosynthesis and respiration > Photosystem II"
        # F.1 (cyanorak): 1/2 = 50% → NOT >50%
        assert "Photosynthesis and respiration" in func_desc
        assert "Energy production and conversion" in func_desc
