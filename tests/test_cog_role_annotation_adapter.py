"""
Unit tests for the COG/role section of
multiomics_kg/adapters/functional_annotation_adapter.py and
multiomics_kg/utils/cyanorak_role_utils.py.

All tests are offline — no downloads, no live Neo4j required.
"""

import json
import pytest
from pathlib import Path

from multiomics_kg.utils.cyanorak_role_utils import parse_cyanorak_role_tree, _derive_parent
from multiomics_kg.adapters.functional_annotation_adapter import (
    COG_FUNCTIONAL_CATEGORIES,
    CogRoleAnnotationAdapter,
    MultiCogRoleAnnotationAdapter,
    _cog_cat_node_id,
    _cyanorak_role_node_id,
    _tigr_role_node_id,
)


# ---------------------------------------------------------------------------
# Minimal test data
# ---------------------------------------------------------------------------

# A small but structurally complete Cyanorak roles tree.
# 6 codes: two roots (0, B) with children, and one B.1.1 leaf.
MINI_CYANORAK_ROLES_TXT = """\
Cyanorak Roles

-

0

Non-coding gene (RNA)

0.1

tRNA

-

0.2

rRNA

-

B

Cellular processes

B.1

Signal transduction

-

B.1.1

Two-component systems

-
"""

# Gene data: exercises all three edge types (COG, CyanorakRole, TigrRole),
# plus empty-list and None variants.
MINI_COG_GENE_DATA: dict = {
    "PMM0001": {
        "locus_tag": "PMM0001",
        "cog_category": ["J", "K"],
        "cyanorak_Role": ["0.1", "B.1"],
        "cyanorak_Role_description": ["tRNA", "Signal transduction"],
        "tIGR_Role": ["12345"],
        "tIGR_Role_description": ["Some tIGR role"],
    },
    "PMM0002": {
        "locus_tag": "PMM0002",
        "cog_category": ["J"],
        "cyanorak_Role": ["0.2"],
        "cyanorak_Role_description": ["rRNA"],
        "tIGR_Role": ["67890"],
        "tIGR_Role_description": ["Another tIGR role"],
    },
    "PMM0003": {
        "locus_tag": "PMM0003",
        "cog_category": [],
        "cyanorak_Role": [],
        "cyanorak_Role_description": [],
        "tIGR_Role": [],
        "tIGR_Role_description": [],
    },
    "PMM0004": {
        "locus_tag": "PMM0004",
        "cog_category": None,
        "cyanorak_Role": None,
        "cyanorak_Role_description": None,
        "tIGR_Role": None,
        "tIGR_Role_description": None,
    },
    # Gene with special characters to validate _clean_str
    "PMM0005": {
        "locus_tag": "PMM0005",
        "cog_category": ["S"],
        "cyanorak_Role": [],
        "cyanorak_Role_description": [],
        "tIGR_Role": ["99999"],
        "tIGR_Role_description": ["Role with 'quotes' and | pipes"],
    },
}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def roles_txt_file(tmp_path):
    """Write MINI_CYANORAK_ROLES_TXT to a temp file and return its path."""
    p = tmp_path / "cyanorak_roles.txt"
    p.write_text(MINI_CYANORAK_ROLES_TXT, encoding="utf-8")
    return p


@pytest.fixture
def role_tree(roles_txt_file):
    """Return the parsed tree dict from the mini roles file."""
    return parse_cyanorak_role_tree(roles_txt_file)


@pytest.fixture
def genome_dir(tmp_path):
    """Write MINI_COG_GENE_DATA as gene_annotations_merged.json; return dir."""
    d = tmp_path / "MED4"
    d.mkdir()
    (d / "gene_annotations_merged.json").write_text(
        json.dumps(MINI_COG_GENE_DATA), encoding="utf-8"
    )
    return d


@pytest.fixture
def cog_adapter(genome_dir):
    return CogRoleAnnotationAdapter(genome_dir=genome_dir)


@pytest.fixture
def genome_config_csv(tmp_path, genome_dir):
    """Write a minimal cyanobacteria_genomes.csv with one strain row."""
    csv_path = tmp_path / "genomes.csv"
    csv_path.write_text(
        "ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir,clade\n"
        f"GCF_000011465.1,Pro_MED4,59919,MED4,{genome_dir},HLI\n"
        "# this comment should be skipped\n",
        encoding="utf-8",
    )
    return str(csv_path)


@pytest.fixture
def multi_adapter(genome_config_csv, roles_txt_file):
    return MultiCogRoleAnnotationAdapter(
        genome_config_file=genome_config_csv,
        role_tree_file=roles_txt_file,
    )


# ===========================================================================
# Tests for parse_cyanorak_role_tree
# ===========================================================================


class TestParseCyanorakRoleTree:
    def test_all_codes_present(self, role_tree):
        assert set(role_tree.keys()) == {"0", "0.1", "0.2", "B", "B.1", "B.1.1"}

    def test_root_codes_have_no_parent(self, role_tree):
        assert role_tree["0"]["parent"] is None
        assert role_tree["B"]["parent"] is None

    def test_child_codes_have_correct_parent(self, role_tree):
        assert role_tree["0.1"]["parent"] == "0"
        assert role_tree["0.2"]["parent"] == "0"
        assert role_tree["B.1"]["parent"] == "B"
        assert role_tree["B.1.1"]["parent"] == "B.1"

    def test_descriptions_parsed(self, role_tree):
        assert role_tree["0"]["description"] == "Non-coding gene (RNA)"
        assert role_tree["0.1"]["description"] == "tRNA"
        assert role_tree["B"]["description"] == "Cellular processes"
        assert role_tree["B.1"]["description"] == "Signal transduction"
        assert role_tree["B.1.1"]["description"] == "Two-component systems"

    def test_returns_dict(self, role_tree):
        assert isinstance(role_tree, dict)

    def test_does_not_include_unclassified_header(self, role_tree):
        # "Cyanorak Roles" or "Unclassified" must not appear as codes
        for code in role_tree:
            assert code not in {"Cyanorak Roles", "Unclassified", "-"}


class TestDeriveParent:
    def test_two_segment(self):
        assert _derive_parent("B.1") == "B"

    def test_three_segment(self):
        assert _derive_parent("B.1.1") == "B.1"

    def test_one_segment_returns_none(self):
        assert _derive_parent("B") is None
        assert _derive_parent("0") is None

    def test_numeric_two_segment(self):
        assert _derive_parent("0.1") == "0"


# ===========================================================================
# Tests for helper ID functions
# ===========================================================================


class TestNodeIdHelpers:
    def test_cog_cat_node_id(self):
        assert _cog_cat_node_id("J") == "cog.category:J"

    def test_cog_cat_node_id_single_letter(self):
        for letter in "JABKLDYVTMNZWUOXCGEFHIPQRS":
            nid = _cog_cat_node_id(letter)
            assert nid == f"cog.category:{letter}"

    def test_cyanorak_role_node_id(self):
        assert _cyanorak_role_node_id("B.1") == "cyanorak.role:B.1"
        assert _cyanorak_role_node_id("0") == "cyanorak.role:0"

    def test_tigr_role_node_id(self):
        assert _tigr_role_node_id("12345") == "tigr.role:12345"


# ===========================================================================
# Tests for CogRoleAnnotationAdapter (per-strain)
# ===========================================================================


class TestCogRoleAdapterEdgeStructure:
    def test_all_edges_are_5_tuples(self, cog_adapter):
        for edge in cog_adapter.get_edges():
            assert len(edge) == 5, f"Edge is not a 5-tuple: {edge}"
            edge_id, src, tgt, label, props = edge
            assert isinstance(edge_id, str)
            assert isinstance(src, str)
            assert isinstance(tgt, str)
            assert isinstance(label, str)
            assert isinstance(props, dict)

    def test_gene_node_id_format(self, cog_adapter):
        for _, src, _, _, _ in cog_adapter.get_edges():
            assert src.startswith("ncbigene:"), f"Unexpected source ID: {src}"

    def test_cog_target_id_format(self, cog_adapter):
        edges = [e for e in cog_adapter.get_edges() if e[3] == "gene_in_cog_category"]
        for _, _, tgt, _, _ in edges:
            assert tgt.startswith("cog.category:"), f"Unexpected COG target ID: {tgt}"

    def test_cyanorak_target_id_format(self, cog_adapter):
        edges = [e for e in cog_adapter.get_edges() if e[3] == "gene_has_cyanorak_role"]
        for _, _, tgt, _, _ in edges:
            assert tgt.startswith("cyanorak.role:"), f"Unexpected CyanorakRole target ID: {tgt}"

    def test_tigr_target_id_format(self, cog_adapter):
        edges = [e for e in cog_adapter.get_edges() if e[3] == "gene_has_tigr_role"]
        for _, _, tgt, _, _ in edges:
            assert tgt.startswith("tigr.role:"), f"Unexpected TigrRole target ID: {tgt}"


class TestCogRoleAdapterEdgeCounts:
    def test_pmm0001_cog_cat_edges(self, cog_adapter):
        edges = [
            e for e in cog_adapter.get_edges()
            if "PMM0001" in e[1] and e[3] == "gene_in_cog_category"
        ]
        assert len(edges) == 2  # "J" and "K"

    def test_pmm0001_cyanorak_role_edges(self, cog_adapter):
        edges = [
            e for e in cog_adapter.get_edges()
            if "PMM0001" in e[1] and e[3] == "gene_has_cyanorak_role"
        ]
        assert len(edges) == 2  # "0.1" and "B.1"

    def test_pmm0001_tigr_role_edges(self, cog_adapter):
        edges = [
            e for e in cog_adapter.get_edges()
            if "PMM0001" in e[1] and e[3] == "gene_has_tigr_role"
        ]
        assert len(edges) == 1  # "12345"

    def test_pmm0002_one_cog_edge(self, cog_adapter):
        edges = [
            e for e in cog_adapter.get_edges()
            if "PMM0002" in e[1] and e[3] == "gene_in_cog_category"
        ]
        assert len(edges) == 1  # "J"

    def test_pmm0003_empty_lists_yields_no_edges(self, cog_adapter):
        edges = [e for e in cog_adapter.get_edges() if "PMM0003" in e[1]]
        assert len(edges) == 0

    def test_pmm0004_null_fields_yields_no_edges(self, cog_adapter):
        edges = [e for e in cog_adapter.get_edges() if "PMM0004" in e[1]]
        assert len(edges) == 0


class TestCogRoleAdapterEdgeIds:
    def test_cog_edge_id_format(self, cog_adapter):
        edges = [e for e in cog_adapter.get_edges() if e[3] == "gene_in_cog_category"]
        for edge_id, _, _, _, _ in edges:
            assert "-cogcat-" in edge_id, f"COG edge ID missing '-cogcat-': {edge_id}"

    def test_cyanorak_edge_id_format(self, cog_adapter):
        edges = [e for e in cog_adapter.get_edges() if e[3] == "gene_has_cyanorak_role"]
        for edge_id, _, _, _, _ in edges:
            assert "-cyrole-" in edge_id, f"CyanorakRole edge ID missing '-cyrole-': {edge_id}"

    def test_tigr_edge_id_format(self, cog_adapter):
        edges = [e for e in cog_adapter.get_edges() if e[3] == "gene_has_tigr_role"]
        for edge_id, _, _, _, _ in edges:
            assert "-tigrrole-" in edge_id, f"TigrRole edge ID missing '-tigrrole-': {edge_id}"

    def test_no_duplicate_edge_ids(self, cog_adapter):
        edges = list(cog_adapter.get_edges())
        ids = [e[0] for e in edges]
        assert len(ids) == len(set(ids)), "Duplicate edge IDs in CogRoleAnnotationAdapter"


class TestCogRoleAdapterGetAllCodes:
    def test_get_all_cyanorak_codes(self, cog_adapter):
        codes = cog_adapter.get_all_cyanorak_codes()
        code_set = {c for c, _ in codes}
        assert "0.1" in code_set
        assert "0.2" in code_set
        assert "B.1" in code_set

    def test_cyanorak_codes_include_descriptions(self, cog_adapter):
        codes = dict(cog_adapter.get_all_cyanorak_codes())
        assert codes.get("0.1") == "tRNA"
        assert codes.get("0.2") == "rRNA"

    def test_get_all_tigr_codes(self, cog_adapter):
        codes = cog_adapter.get_all_tigr_codes()
        code_set = {c for c, _ in codes}
        assert "12345" in code_set
        assert "67890" in code_set

    def test_tigr_codes_include_descriptions(self, cog_adapter):
        codes = dict(cog_adapter.get_all_tigr_codes())
        assert codes.get("12345") == "Some tIGR role"

    def test_null_fields_not_in_cyanorak_codes(self, cog_adapter):
        codes = {c for c, _ in cog_adapter.get_all_cyanorak_codes()}
        assert None not in codes
        assert "" not in codes


class TestCogRoleAdapterMissingFile:
    def test_missing_json_yields_no_edges(self, tmp_path):
        empty_dir = tmp_path / "empty_strain"
        empty_dir.mkdir()
        adapter = CogRoleAnnotationAdapter(genome_dir=empty_dir)
        assert list(adapter.get_edges()) == []

    def test_missing_json_get_all_cyanorak_codes_returns_empty(self, tmp_path):
        empty_dir = tmp_path / "empty_strain"
        empty_dir.mkdir()
        adapter = CogRoleAnnotationAdapter(genome_dir=empty_dir)
        assert adapter.get_all_cyanorak_codes() == set()

    def test_missing_json_get_all_tigr_codes_returns_empty(self, tmp_path):
        empty_dir = tmp_path / "empty_strain"
        empty_dir.mkdir()
        adapter = CogRoleAnnotationAdapter(genome_dir=empty_dir)
        assert adapter.get_all_tigr_codes() == set()


class TestCogRoleAdapterTestMode:
    def test_test_mode_limits_cog_edges(self, tmp_path):
        big_genes = {
            f"PMM{i:04d}": {
                "locus_tag": f"PMM{i:04d}",
                "cog_category": ["J"],
                "cyanorak_Role": [],
                "cyanorak_Role_description": [],
                "tIGR_Role": [],
                "tIGR_Role_description": [],
            }
            for i in range(200)
        }
        gene_file = tmp_path / "gene_annotations_merged.json"
        gene_file.write_text(json.dumps(big_genes), encoding="utf-8")
        adapter = CogRoleAnnotationAdapter(genome_dir=tmp_path, test_mode=True)
        edges = list(adapter.get_edges())
        assert len(edges) <= 100


# ===========================================================================
# Tests for MultiCogRoleAnnotationAdapter.get_nodes()
# ===========================================================================


class TestMultiCogRoleAdapterCogNodes:
    def test_cog_node_count_matches_dict(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        cog_nodes = [n for n in nodes if n[1] == "cog functional category"]
        assert len(cog_nodes) == len(COG_FUNCTIONAL_CATEGORIES)

    def test_all_standard_letters_present(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        cog_ids = {n[0] for n in nodes if n[1] == "cog functional category"}
        for letter in COG_FUNCTIONAL_CATEGORIES:
            assert _cog_cat_node_id(letter) in cog_ids, f"Missing COG node for letter {letter}"

    def test_cog_node_has_code_and_name(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        for node_id, label, props in nodes:
            if label == "cog functional category":
                assert "code" in props, f"COG node {node_id} missing 'code'"
                assert "name" in props, f"COG node {node_id} missing 'name'"
                assert isinstance(props["code"], str)
                assert isinstance(props["name"], str)

    def test_cog_node_code_matches_id(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        for node_id, label, props in nodes:
            if label == "cog functional category":
                assert node_id == _cog_cat_node_id(props["code"])


class TestMultiCogRoleAdapterCyanorakNodes:
    def test_cyanorak_role_nodes_from_tree(self, multi_adapter, role_tree):
        nodes = list(multi_adapter.get_nodes())
        cyr_nodes = [n for n in nodes if n[1] == "cyanorak role"]
        # All 6 codes in MINI_CYANORAK_ROLES_TXT should yield nodes
        assert len(cyr_nodes) == len(role_tree)

    def test_cyanorak_node_ids_match_tree(self, multi_adapter, role_tree):
        nodes = list(multi_adapter.get_nodes())
        cyr_ids = {n[0] for n in nodes if n[1] == "cyanorak role"}
        for code in role_tree:
            assert _cyanorak_role_node_id(code) in cyr_ids, f"Missing CyanorakRole node: {code}"

    def test_cyanorak_node_has_code_and_description(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        for node_id, label, props in nodes:
            if label == "cyanorak role":
                assert "code" in props, f"CyanorakRole node {node_id} missing 'code'"
                assert "description" in props, f"CyanorakRole node {node_id} missing 'description'"

    def test_cyanorak_node_description_correct(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        cyr_nodes = {n[0]: n[2] for n in nodes if n[1] == "cyanorak role"}
        assert cyr_nodes[_cyanorak_role_node_id("0.1")]["description"] == "tRNA"
        assert cyr_nodes[_cyanorak_role_node_id("B")]["description"] == "Cellular processes"


class TestMultiCogRoleAdapterTigrNodes:
    def test_tigr_nodes_only_from_data(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        tigr_nodes = [n for n in nodes if n[1] == "tigr role"]
        # Data has: 12345 (PMM0001), 67890 (PMM0002), 99999 (PMM0005)
        assert len(tigr_nodes) == 3

    def test_tigr_node_ids_present(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        tigr_ids = {n[0] for n in nodes if n[1] == "tigr role"}
        assert _tigr_role_node_id("12345") in tigr_ids
        assert _tigr_role_node_id("67890") in tigr_ids
        assert _tigr_role_node_id("99999") in tigr_ids

    def test_tigr_node_has_code_and_description(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        for node_id, label, props in nodes:
            if label == "tigr role":
                assert "code" in props
                assert "description" in props


class TestMultiCogRoleAdapterNoDuplicateNodes:
    def test_no_duplicate_node_ids(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        node_ids = [n[0] for n in nodes]
        assert len(node_ids) == len(set(node_ids)), "Duplicate node IDs in get_nodes()"

    def test_tigr_deduplication_across_strains(self, tmp_path, roles_txt_file):
        """Two strains sharing the same tIGR code should yield only one TigrRole node."""
        shared_genes = {
            "PMM0001": {
                "locus_tag": "PMM0001",
                "cog_category": ["J"],
                "cyanorak_Role": [],
                "cyanorak_Role_description": [],
                "tIGR_Role": ["12345"],
                "tIGR_Role_description": ["Shared role"],
            }
        }
        strain1 = tmp_path / "MED4"
        strain1.mkdir()
        (strain1 / "gene_annotations_merged.json").write_text(json.dumps(shared_genes), encoding="utf-8")
        strain2 = tmp_path / "MIT9312"
        strain2.mkdir()
        (strain2 / "gene_annotations_merged.json").write_text(json.dumps(shared_genes), encoding="utf-8")

        csv_path = tmp_path / "genomes.csv"
        csv_path.write_text(
            "ncbi_accession,ncbi_taxon_id,strain_name,data_dir,clade\n"
            f"GCF_001,59919,MED4,{strain1},HLI\n"
            f"GCF_002,59919,MIT9312,{strain2},HLI\n",
            encoding="utf-8",
        )
        adapter = MultiCogRoleAnnotationAdapter(
            genome_config_file=str(csv_path),
            role_tree_file=roles_txt_file,
        )
        nodes = list(adapter.get_nodes())
        tigr_nodes = [n for n in nodes if n[1] == "tigr role"]
        tigr_ids = [n[0] for n in tigr_nodes]
        assert len(tigr_ids) == len(set(tigr_ids)), "Duplicate TigrRole nodes across strains"


# ===========================================================================
# Tests for MultiCogRoleAnnotationAdapter.get_edges()
# ===========================================================================


class TestMultiCogRoleAdapterEdges:
    def test_all_expected_labels_present(self, multi_adapter):
        edges = list(multi_adapter.get_edges())
        labels = {e[3] for e in edges}
        assert "gene_in_cog_category" in labels
        assert "gene_has_cyanorak_role" in labels
        assert "gene_has_tigr_role" in labels
        assert "cyanorak_role_is_a_cyanorak_role" in labels

    def test_cog_cat_edge_count(self, multi_adapter):
        edges = list(multi_adapter.get_edges())
        cog = [e for e in edges if e[3] == "gene_in_cog_category"]
        # PMM0001: J,K (2); PMM0002: J (1); PMM0003/4: none; PMM0005: S (1) → 4
        assert len(cog) == 4

    def test_cyanorak_role_edge_count(self, multi_adapter):
        edges = list(multi_adapter.get_edges())
        cyr = [e for e in edges if e[3] == "gene_has_cyanorak_role"]
        # PMM0001: 0.1, B.1 (2); PMM0002: 0.2 (1) → 3
        assert len(cyr) == 3

    def test_tigr_role_edge_count(self, multi_adapter):
        edges = list(multi_adapter.get_edges())
        tigr = [e for e in edges if e[3] == "gene_has_tigr_role"]
        # PMM0001: 12345 (1); PMM0002: 67890 (1); PMM0005: 99999 (1) → 3
        assert len(tigr) == 3

    def test_cyanorak_hierarchy_edge_count(self, multi_adapter, role_tree):
        edges = list(multi_adapter.get_edges())
        hier = [e for e in edges if e[3] == "cyanorak_role_is_a_cyanorak_role"]
        # Edges: 0.1→0, 0.2→0, B.1→B, B.1.1→B.1 → 4 hierarchy edges
        expected = sum(1 for entry in role_tree.values() if entry["parent"] is not None)
        assert len(hier) == expected

    def test_hierarchy_edge_ids(self, multi_adapter):
        edges = list(multi_adapter.get_edges())
        for edge in edges:
            if edge[3] == "cyanorak_role_is_a_cyanorak_role":
                edge_id, child, parent, _, _ = edge
                assert edge_id == f"{edge_id}", "Edge ID must be a string"
                assert child.startswith("cyanorak.role:")
                assert parent.startswith("cyanorak.role:")

    def test_hierarchy_edge_parent_child_correctness(self, multi_adapter):
        edges = list(multi_adapter.get_edges())
        hier = {
            (e[1], e[2]) for e in edges if e[3] == "cyanorak_role_is_a_cyanorak_role"
        }
        # 0.1 → 0
        assert (_cyanorak_role_node_id("0.1"), _cyanorak_role_node_id("0")) in hier
        # B.1.1 → B.1
        assert (_cyanorak_role_node_id("B.1.1"), _cyanorak_role_node_id("B.1")) in hier

    def test_gene_cog_edge_targets_known_nodes(self, multi_adapter):
        """All gene→COG edges must point to nodes returned by get_nodes()."""
        node_ids = {n[0] for n in multi_adapter.get_nodes()}
        edges = list(multi_adapter.get_edges())
        for e in edges:
            if e[3] == "gene_in_cog_category":
                assert e[2] in node_ids, f"COG target {e[2]} not in node set"


# ===========================================================================
# Tests for string sanitization (_clean_str)
# ===========================================================================


class TestStringCleanStr:
    def test_single_quotes_sanitized_in_cog_nodes(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        for _, label, props in nodes:
            for val in props.values():
                assert "'" not in str(val), f"Raw single-quote in {label} property: {val!r}"

    def test_pipe_sanitized_in_cog_nodes(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        for _, label, props in nodes:
            for val in props.values():
                assert "|" not in str(val), f"Pipe character in {label} property: {val!r}"

    def test_tigr_description_with_special_chars_is_cleaned(self, multi_adapter):
        """PMM0005's tIGR description contains quotes and pipes — must be sanitized."""
        nodes = list(multi_adapter.get_nodes())
        tigr_99999 = next(
            (n for n in nodes if n[1] == "tigr role" and "99999" in n[0]), None
        )
        assert tigr_99999 is not None, "TigrRole node for 99999 not found"
        desc = tigr_99999[2]["description"]
        assert "'" not in desc
        assert "|" not in desc


# ===========================================================================
# Tests for CSV parsing (comments + empty data_dir)
# ===========================================================================


class TestMultiCogRoleCsvParsing:
    def test_comment_lines_skipped(self, tmp_path, genome_dir, roles_txt_file):
        csv_path = tmp_path / "genomes.csv"
        csv_path.write_text(
            "ncbi_accession,ncbi_taxon_id,strain_name,data_dir,clade\n"
            "# this is a comment\n"
            f"GCF_000011465.1,59919,MED4,{genome_dir},HLI\n",
            encoding="utf-8",
        )
        adapter = MultiCogRoleAnnotationAdapter(
            genome_config_file=str(csv_path),
            role_tree_file=roles_txt_file,
        )
        assert len(adapter._strain_adapters) == 1

    def test_empty_data_dir_skipped(self, tmp_path, roles_txt_file):
        csv_path = tmp_path / "genomes.csv"
        csv_path.write_text(
            "ncbi_accession,ncbi_taxon_id,strain_name,data_dir,clade\n"
            "GCF_000011465.1,59919,MED4,,HLI\n",  # empty data_dir
            encoding="utf-8",
        )
        adapter = MultiCogRoleAnnotationAdapter(
            genome_config_file=str(csv_path),
            role_tree_file=roles_txt_file,
        )
        assert len(adapter._strain_adapters) == 0
