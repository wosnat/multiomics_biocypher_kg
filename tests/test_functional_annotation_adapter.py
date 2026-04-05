"""
Unit tests for multiomics_kg/utils/go_utils.py and
multiomics_kg/adapters/functional_annotation_adapter.py.

All tests are offline — no downloads, no live Neo4j required.
"""

import json
import pytest
from pathlib import Path
from unittest.mock import patch

from multiomics_kg.utils.go_utils import (
    NAMESPACE_TO_LABEL,
    compute_ancestry_closure,
    load_go_data,
    make_go_go_edge_label,
    parse_obo,
)
from multiomics_kg.adapters.functional_annotation_adapter import (
    EcAnnotationAdapter,
    GoAnnotationAdapter,
    MultiEcAnnotationAdapter,
    MultiGoAnnotationAdapter,
)


# ---------------------------------------------------------------------------
# Shared test data
# ---------------------------------------------------------------------------

MINI_OBO_CONTENT = """\
format-version: 1.2

[Term]
id: GO:0005737
name: cytoplasm
namespace: cellular_component
is_a: GO:0005622 ! intracellular anatomical structure

[Term]
id: GO:0005622
name: intracellular anatomical structure
namespace: cellular_component

[Term]
id: GO:0006260
name: DNA replication
namespace: biological_process
is_a: GO:0071897 ! DNA biosynthetic process

[Term]
id: GO:0071897
name: DNA biosynthetic process
namespace: biological_process

[Term]
id: GO:0003677
name: DNA binding
namespace: molecular_function

[Term]
id: GO:0000001
name: obsolete term
namespace: biological_process
is_obsolete: true

[Typedef]
id: part_of
name: part of
"""

# Equivalent of what parse_obo(MINI_OBO_CONTENT) should return
MINI_GO_DATA: dict = {
    "GO:0005737": {
        "name": "cytoplasm",
        "namespace": "cellular_component",
        "parents": [["GO:0005622", "is_a"]],
    },
    "GO:0005622": {
        "name": "intracellular anatomical structure",
        "namespace": "cellular_component",
        "parents": [],
    },
    "GO:0006260": {
        "name": "DNA replication",
        "namespace": "biological_process",
        "parents": [["GO:0071897", "is_a"]],
    },
    "GO:0071897": {
        "name": "DNA biosynthetic process",
        "namespace": "biological_process",
        "parents": [],
    },
    "GO:0003677": {
        "name": "DNA binding",
        "namespace": "molecular_function",
        "parents": [],
    },
}

MINI_GENE_DATA: dict = {
    "PMM0001": {
        "locus_tag": "PMM0001",
        "go_terms": ["GO:0005737", "GO:0006260", "GO:0003677"],
    },
    "PMM0002": {
        "locus_tag": "PMM0002",
        "go_terms": ["GO:0006260"],
    },
    "PMM0003": {
        "locus_tag": "PMM0003",
        "go_terms": [],  # no GO annotations
    },
    "PMM0004": {
        "locus_tag": "PMM0004",
        "go_terms": ["GO:9999999"],  # unknown / not in go_data
    },
}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def obo_file(tmp_path):
    """Write MINI_OBO_CONTENT to a temp .obo file and return its path."""
    p = tmp_path / "go-basic.obo"
    p.write_text(MINI_OBO_CONTENT, encoding="utf-8")
    return p


@pytest.fixture
def tmp_go_cache(tmp_path):
    """Write MINI_GO_DATA as go_namespace_cache.json, return cache_root."""
    go_dir = tmp_path / "go_terms"
    go_dir.mkdir()
    (go_dir / "go_namespace_cache.json").write_text(
        json.dumps(MINI_GO_DATA), encoding="utf-8"
    )
    return tmp_path


@pytest.fixture
def tmp_genome_dir(tmp_path):
    """Write MINI_GENE_DATA as gene_annotations_merged.json, return genome_dir."""
    (tmp_path / "gene_annotations_merged.json").write_text(
        json.dumps(MINI_GENE_DATA), encoding="utf-8"
    )
    return tmp_path


@pytest.fixture
def go_adapter(tmp_genome_dir):
    return GoAnnotationAdapter(genome_dir=tmp_genome_dir, go_data=MINI_GO_DATA)


@pytest.fixture
def multi_config_csv(tmp_path, tmp_genome_dir):
    """Write a cyanobacteria_genomes.csv with one data row + one comment."""
    csv_path = tmp_path / "genomes.csv"
    csv_path.write_text(
        "ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir,clade\n"
        f"GCF_000011465.1,Pro_MED4,59919,MED4,{tmp_genome_dir},HLI\n"
        "# this comment should be skipped\n",
        encoding="utf-8",
    )
    return str(csv_path)


@pytest.fixture
def multi_adapter(multi_config_csv, tmp_go_cache):
    return MultiGoAnnotationAdapter(
        genome_config_file=multi_config_csv,
        cache_root=tmp_go_cache,
    )


# ---------------------------------------------------------------------------
# TestParseObo
# ---------------------------------------------------------------------------


class TestParseObo:
    def test_parses_non_obsolete_terms(self, obo_file):
        result = parse_obo(obo_file)
        assert len(result) == 5  # 5 non-obsolete terms

    def test_skips_obsolete_term(self, obo_file):
        result = parse_obo(obo_file)
        assert "GO:0000001" not in result

    def test_namespace_correct(self, obo_file):
        result = parse_obo(obo_file)
        assert result["GO:0005737"]["namespace"] == "cellular_component"
        assert result["GO:0006260"]["namespace"] == "biological_process"
        assert result["GO:0003677"]["namespace"] == "molecular_function"

    def test_name_correct(self, obo_file):
        result = parse_obo(obo_file)
        assert result["GO:0005737"]["name"] == "cytoplasm"
        assert result["GO:0006260"]["name"] == "DNA replication"

    def test_parents_extracted(self, obo_file):
        result = parse_obo(obo_file)
        assert result["GO:0005737"]["parents"] == [["GO:0005622", "is_a"]]
        assert result["GO:0006260"]["parents"] == [["GO:0071897", "is_a"]]

    def test_root_has_no_parents(self, obo_file):
        result = parse_obo(obo_file)
        assert result["GO:0005622"]["parents"] == []
        assert result["GO:0003677"]["parents"] == []


# ---------------------------------------------------------------------------
# TestAncestryClosure
# ---------------------------------------------------------------------------


class TestAncestryClosure:
    def test_seed_included(self):
        closure = compute_ancestry_closure({"GO:0005737"}, MINI_GO_DATA)
        assert "GO:0005737" in closure

    def test_direct_ancestor_included(self):
        closure = compute_ancestry_closure({"GO:0005737"}, MINI_GO_DATA)
        assert "GO:0005622" in closure

    def test_unrelated_term_excluded(self):
        closure = compute_ancestry_closure({"GO:0005737"}, MINI_GO_DATA)
        assert "GO:0003677" not in closure
        assert "GO:0006260" not in closure

    def test_multi_seed_union(self):
        closure = compute_ancestry_closure({"GO:0005737", "GO:0006260"}, MINI_GO_DATA)
        assert "GO:0005622" in closure   # ancestor of GO:0005737
        assert "GO:0071897" in closure   # ancestor of GO:0006260

    def test_empty_seed_returns_empty(self):
        closure = compute_ancestry_closure(set(), MINI_GO_DATA)
        assert closure == set()

    def test_unknown_seed_ignored(self):
        # Unknown term not in go_data → silently excluded from closure
        closure = compute_ancestry_closure({"GO:9999999"}, MINI_GO_DATA)
        assert closure == set()


# ---------------------------------------------------------------------------
# TestLoadGoData
# ---------------------------------------------------------------------------


class TestLoadGoData:
    def test_loads_from_existing_cache(self, tmp_go_cache):
        data = load_go_data(tmp_go_cache)
        assert "GO:0005737" in data
        assert data["GO:0005737"]["name"] == "cytoplasm"

    def test_does_not_re_download_when_cached(self, tmp_go_cache):
        """load_go_data must not download OBO when cache exists."""
        with patch("multiomics_kg.utils.go_utils._download_obo") as mock_dl:
            load_go_data(tmp_go_cache, force=False)
        mock_dl.assert_not_called()

    def test_force_rebuilds_cache(self, tmp_go_cache, obo_file):
        """With force=True, cache is rebuilt. We mock download to supply our OBO."""
        def fake_download(path):
            import shutil
            shutil.copy(obo_file, path)

        with patch("multiomics_kg.utils.go_utils._download_obo", side_effect=fake_download):
            data = load_go_data(tmp_go_cache, force=True)

        assert "GO:0005737" in data
        assert "GO:0000001" not in data  # obsolete excluded


# ---------------------------------------------------------------------------
# TestMakeGoGoEdgeLabel
# ---------------------------------------------------------------------------


class TestMakeGoGoEdgeLabel:
    def test_valid_bp_is_a_bp(self):
        label = make_go_go_edge_label("biological_process", "is_a", "biological_process")
        assert label == "biological_process_is_a_biological_process"

    def test_valid_cc_part_of_cc(self):
        label = make_go_go_edge_label("cellular_component", "part_of", "cellular_component")
        assert label == "cellular_component_part_of_cellular_component"

    def test_unknown_relation_returns_none(self):
        label = make_go_go_edge_label("biological_process", "occurs_in", "cellular_component")
        assert label is None

    def test_unknown_namespace_combination_returns_none(self):
        # CC → BP not in schema
        label = make_go_go_edge_label("cellular_component", "is_a", "biological_process")
        assert label is None


# ---------------------------------------------------------------------------
# TestGoAnnotationAdapterEdges
# ---------------------------------------------------------------------------


class TestGoAnnotationAdapterEdges:
    def test_pmm0001_yields_three_edges(self, go_adapter):
        edges = [e for e in go_adapter.get_edges() if e[1].endswith("PMM0001")]
        assert len(edges) == 3

    def test_pmm0002_yields_one_edge(self, go_adapter):
        edges = [e for e in go_adapter.get_edges() if e[1].endswith("PMM0002")]
        assert len(edges) == 1

    def test_pmm0003_yields_no_edges(self, go_adapter):
        edges = [e for e in go_adapter.get_edges() if e[1].endswith("PMM0003")]
        assert len(edges) == 0

    def test_unknown_go_id_skipped(self, go_adapter):
        # PMM0004 has GO:9999999 which is not in MINI_GO_DATA
        edges = [e for e in go_adapter.get_edges() if e[1].endswith("PMM0004")]
        assert len(edges) == 0

    def test_edge_structure(self, go_adapter):
        edges = list(go_adapter.get_edges())
        for edge in edges:
            assert len(edge) == 5, "Edge must be a 5-tuple"
            edge_id, source, target, label, props = edge
            assert isinstance(edge_id, str)
            assert isinstance(source, str)
            assert isinstance(target, str)
            assert isinstance(label, str)
            assert isinstance(props, dict)

    def test_biological_process_label(self, go_adapter):
        edges = list(go_adapter.get_edges())
        bp_edges = [e for e in edges if e[3] == "gene_involved_in_biological_process"]
        assert len(bp_edges) >= 1
        # GO:0006260 is biological_process → should produce a BP edge
        assert any("GO:0006260" in e[0] or "GO:0006260" in e[2] for e in bp_edges)

    def test_cellular_component_label(self, go_adapter):
        edges = list(go_adapter.get_edges())
        cc_edges = [e for e in edges if e[3] == "gene_located_in_cellular_component"]
        assert len(cc_edges) >= 1

    def test_molecular_function_label(self, go_adapter):
        edges = list(go_adapter.get_edges())
        mf_edges = [e for e in edges if e[3] == "gene_enables_molecular_function"]
        assert len(mf_edges) >= 1


# ---------------------------------------------------------------------------
# TestGoAnnotationAdapterIds
# ---------------------------------------------------------------------------


class TestGoAnnotationAdapterIds:
    def test_gene_node_id_format(self, go_adapter):
        edges = list(go_adapter.get_edges())
        gene_ids = {e[1] for e in edges}
        assert all(ids.startswith("ncbigene:") for ids in gene_ids)

    def test_go_node_id_format(self, go_adapter):
        edges = list(go_adapter.get_edges())
        go_ids = {e[2] for e in edges}
        # normalize_curie("go:GO:0005737") → "go:0005737"
        assert all(gid.startswith("go:") for gid in go_ids)

    def test_edge_id_format(self, go_adapter):
        edges = list(go_adapter.get_edges())
        for edge_id, source, target, label, _ in edges:
            assert "-go-GO:" in edge_id


# ---------------------------------------------------------------------------
# TestGoAnnotationAdapterTestMode
# ---------------------------------------------------------------------------


class TestGoAnnotationAdapterTestMode:
    def test_test_mode_limits_edges(self, tmp_genome_dir):
        # Create 200 genes each with 1 GO term → normally 200 edges
        big_genes = {
            f"PMM{i:04d}": {"locus_tag": f"PMM{i:04d}", "go_terms": ["GO:0005737"]}
            for i in range(200)
        }
        (tmp_genome_dir / "gene_annotations_merged.json").write_text(
            json.dumps(big_genes), encoding="utf-8"
        )
        adapter = GoAnnotationAdapter(
            genome_dir=tmp_genome_dir, go_data=MINI_GO_DATA, test_mode=True
        )
        edges = list(adapter.get_edges())
        assert len(edges) <= 100


# ---------------------------------------------------------------------------
# TestMultiGoAnnotationAdapterNodes
# ---------------------------------------------------------------------------


class TestMultiGoAnnotationAdapterNodes:
    def test_go_nodes_emitted(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        node_ids = {n[0] for n in nodes}
        # GO:0005737 and its ancestor GO:0005622 should both be present
        assert any("5737" in nid for nid in node_ids)
        assert any("5622" in nid for nid in node_ids)

    def test_closure_includes_ancestors(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        node_ids = {n[0] for n in nodes}
        # GO:0071897 is ancestor of GO:0006260 (both BP)
        assert any("71897" in nid for nid in node_ids)

    def test_node_label_cellular_component(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        cc_nodes = [n for n in nodes if n[1] == "cellular component"]
        assert len(cc_nodes) >= 1

    def test_node_label_biological_process(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        bp_nodes = [n for n in nodes if n[1] == "biological process"]
        assert len(bp_nodes) >= 1

    def test_node_label_molecular_function(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        mf_nodes = [n for n in nodes if n[1] == "molecular function"]
        assert len(mf_nodes) >= 1

    def test_node_has_name_property(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        for node_id, label, props in nodes:
            assert "name" in props, f"Node {node_id} missing 'name' property"
            assert isinstance(props["name"], str)

    def test_no_duplicate_nodes(self, multi_adapter):
        nodes = list(multi_adapter.get_nodes())
        node_ids = [n[0] for n in nodes]
        assert len(node_ids) == len(set(node_ids)), "Duplicate GO node IDs"


# ---------------------------------------------------------------------------
# TestMultiGoAnnotationAdapterEdges
# ---------------------------------------------------------------------------


class TestMultiGoAnnotationAdapterEdges:
    def test_gene_to_go_edges_present(self, multi_adapter):
        edges = list(multi_adapter.get_edges())
        gene_go = [e for e in edges if e[3].startswith("gene_")]
        assert len(gene_go) >= 1

    def test_go_go_edges_present(self, multi_adapter):
        edges = list(multi_adapter.get_edges())
        go_go = [e for e in edges if not e[3].startswith("gene_")]
        assert len(go_go) >= 1

    def test_go_go_edge_label_is_valid(self, multi_adapter):
        from multiomics_kg.utils.go_utils import _ALLOWED_GO_GO_EDGE_LABELS
        edges = list(multi_adapter.get_edges())
        for e in edges:
            if not e[3].startswith("gene_"):
                assert e[3] in _ALLOWED_GO_GO_EDGE_LABELS, f"Invalid GO-GO label: {e[3]}"

    def test_no_duplicate_go_go_edges(self, multi_adapter):
        edges = list(multi_adapter.get_edges())
        go_go = [(e[1], e[2], e[3]) for e in edges if not e[3].startswith("gene_")]
        assert len(go_go) == len(set(go_go)), "Duplicate GO-GO edges"


# ---------------------------------------------------------------------------
# TestMultiAdapterCsvParsing
# ---------------------------------------------------------------------------


class TestMultiAdapterCsvParsing:
    def test_comment_lines_skipped(self, tmp_path, tmp_genome_dir, tmp_go_cache):
        csv_path = tmp_path / "genomes.csv"
        csv_path.write_text(
            "ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir,clade\n"
            "# comment\n"
            f"GCF_000011465.1,Pro_MED4,59919,MED4,{tmp_genome_dir},HLI\n",
            encoding="utf-8",
        )
        adapter = MultiGoAnnotationAdapter(
            genome_config_file=str(csv_path),
            cache_root=tmp_go_cache,
        )
        assert len(adapter.adapters) == 1

    def test_empty_data_dir_skipped(self, tmp_path, tmp_go_cache):
        csv_path = tmp_path / "genomes.csv"
        csv_path.write_text(
            "ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir,clade\n"
            "GCF_000011465.1,Pro_MED4,59919,MED4,,HLI\n",  # empty data_dir
            encoding="utf-8",
        )
        adapter = MultiGoAnnotationAdapter(
            genome_config_file=str(csv_path),
            cache_root=tmp_go_cache,
        )
        assert len(adapter.adapters) == 0


# ===========================================================================
# EC annotation tests
# ===========================================================================

# ---------------------------------------------------------------------------
# Shared EC test data
# ---------------------------------------------------------------------------

# Minimal ec_dict + enzymes structure that mirrors what ec_adapter builds.
# 2 level-1 categories → 2 L2 → 2 L3 → 2 leaf enzymes (1.1.1.1 and 2.7.7.7).
MINI_EC_CACHE: dict = {
    "ec_dict": {
        "1.-.-.-": {
            "name": "Oxidoreductases",
            "1.1.-.-": {
                "name": "Acting on the CH-OH group of donors",
                "1.1.1.-": {
                    "name": "With NAD+ as acceptor",
                    "entries": ["1.1.1.1"],
                },
            },
        },
        "2.-.-.-": {
            "name": "Transferases",
            "2.7.-.-": {
                "name": "Transferring phosphorus-containing groups",
                "2.7.7.-": {
                    "name": "Nucleotidyltransferases",
                    "entries": ["2.7.7.7"],
                },
            },
        },
    },
    "enzymes": {
        "1.1.1.1": {"de": "Alcohol dehydrogenase", "an": [], "ca": [], "cc": []},
        "2.7.7.7": {"de": "DNA-directed DNA polymerase", "an": [], "ca": [], "cc": []},
    },
}

MINI_EC_GENE_DATA: dict = {
    "PMM0001": {"locus_tag": "PMM0001", "ec_numbers": ["2.7.7.7"]},
    "PMM0002": {"locus_tag": "PMM0002", "ec_numbers": ["1.1.1.1", "2.7.7.7"]},
    "PMM0003": {"locus_tag": "PMM0003", "ec_numbers": []},
    "PMM0004": {"locus_tag": "PMM0004", "ec_numbers": None},
}


# ---------------------------------------------------------------------------
# EC Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def tmp_ec_cache_dir(tmp_path):
    """Write MINI_EC_CACHE as ec_data.json; return the cache directory."""
    ec_dir = tmp_path / "ec"
    ec_dir.mkdir()
    (ec_dir / "ec_data.json").write_text(
        json.dumps(MINI_EC_CACHE), encoding="utf-8"
    )
    return ec_dir


@pytest.fixture
def tmp_ec_genome_dir(tmp_path):
    """Write MINI_EC_GENE_DATA as gene_annotations_merged.json; return genome_dir."""
    genome_dir = tmp_path / "genome"
    genome_dir.mkdir()
    (genome_dir / "gene_annotations_merged.json").write_text(
        json.dumps(MINI_EC_GENE_DATA), encoding="utf-8"
    )
    return genome_dir


@pytest.fixture
def ec_strain_adapter(tmp_ec_genome_dir):
    return EcAnnotationAdapter(genome_dir=tmp_ec_genome_dir)


@pytest.fixture
def multi_ec_config_csv(tmp_path, tmp_ec_genome_dir):
    csv_path = tmp_path / "ec_genomes.csv"
    csv_path.write_text(
        "ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir,clade\n"
        f"GCF_000011465.1,Pro_MED4,59919,MED4,{tmp_ec_genome_dir},HLI\n"
        "# comment should be skipped\n",
        encoding="utf-8",
    )
    return str(csv_path)


@pytest.fixture
def multi_ec_adapter(multi_ec_config_csv, tmp_ec_cache_dir):
    return MultiEcAnnotationAdapter(
        genome_config_file=multi_ec_config_csv,
        cache_dir=tmp_ec_cache_dir,
        cache=True,
    )


# ---------------------------------------------------------------------------
# TestEcAnnotationAdapterEdges
# ---------------------------------------------------------------------------


class TestEcAnnotationAdapterEdges:
    def test_pmm0001_yields_one_edge(self, ec_strain_adapter):
        edges = [e for e in ec_strain_adapter.get_edges() if "PMM0001" in e[1]]
        assert len(edges) == 1

    def test_pmm0002_yields_two_edges(self, ec_strain_adapter):
        edges = [e for e in ec_strain_adapter.get_edges() if "PMM0002" in e[1]]
        assert len(edges) == 2

    def test_pmm0003_yields_no_edges(self, ec_strain_adapter):
        edges = [e for e in ec_strain_adapter.get_edges() if "PMM0003" in e[1]]
        assert len(edges) == 0

    def test_pmm0004_null_ec_yields_no_edges(self, ec_strain_adapter):
        edges = [e for e in ec_strain_adapter.get_edges() if "PMM0004" in e[1]]
        assert len(edges) == 0

    def test_edge_structure(self, ec_strain_adapter):
        for edge in ec_strain_adapter.get_edges():
            assert len(edge) == 5, "Edge must be a 5-tuple"
            edge_id, source, target, label, props = edge
            assert isinstance(edge_id, str)
            assert isinstance(source, str)
            assert isinstance(target, str)
            assert isinstance(label, str)
            assert isinstance(props, dict)

    def test_edge_label(self, ec_strain_adapter):
        for _, _, _, label, _ in ec_strain_adapter.get_edges():
            assert label == "gene_catalyzes_ec_number"

    def test_gene_node_id_format(self, ec_strain_adapter):
        for _, source, _, _, _ in ec_strain_adapter.get_edges():
            assert source.startswith("ncbigene:"), f"Unexpected gene ID: {source}"

    def test_ec_node_id_format(self, ec_strain_adapter):
        for _, _, target, _, _ in ec_strain_adapter.get_edges():
            assert target.startswith("ec:"), f"Unexpected EC node ID: {target}"

    def test_edge_id_format(self, ec_strain_adapter):
        for edge_id, _, _, _, _ in ec_strain_adapter.get_edges():
            assert "-ec-" in edge_id, f"Edge ID missing '-ec-': {edge_id}"

    def test_no_duplicate_edges(self, ec_strain_adapter):
        edges = list(ec_strain_adapter.get_edges())
        ids = [e[0] for e in edges]
        assert len(ids) == len(set(ids)), "Duplicate edge IDs"

    def test_missing_json_file_yields_no_edges(self, tmp_path):
        adapter = EcAnnotationAdapter(genome_dir=tmp_path)
        assert list(adapter.get_edges()) == []


# ---------------------------------------------------------------------------
# TestEcAnnotationAdapterTestMode
# ---------------------------------------------------------------------------


class TestEcAnnotationAdapterTestMode:
    def test_test_mode_limits_edges(self, tmp_ec_genome_dir):
        big_genes = {
            f"PMM{i:04d}": {"locus_tag": f"PMM{i:04d}", "ec_numbers": ["2.7.7.7"]}
            for i in range(200)
        }
        (tmp_ec_genome_dir / "gene_annotations_merged.json").write_text(
            json.dumps(big_genes), encoding="utf-8"
        )
        adapter = EcAnnotationAdapter(genome_dir=tmp_ec_genome_dir, test_mode=True)
        assert len(list(adapter.get_edges())) <= 100


# ---------------------------------------------------------------------------
# TestMultiEcAnnotationAdapterNodes
# ---------------------------------------------------------------------------


class TestMultiEcAnnotationAdapterNodes:
    def test_ec_nodes_emitted(self, multi_ec_adapter):
        nodes = list(multi_ec_adapter.get_nodes())
        assert len(nodes) > 0

    def test_all_eight_nodes_present(self, multi_ec_adapter):
        # MINI_EC_CACHE has 2 L1 + 2 L2 + 2 L3 + 2 L4 = 8 nodes
        nodes = list(multi_ec_adapter.get_nodes())
        assert len(nodes) == 8

    def test_leaf_node_ids(self, multi_ec_adapter):
        node_ids = {n[0] for n in multi_ec_adapter.get_nodes()}
        assert "ec:1.1.1.1" in node_ids
        assert "ec:2.7.7.7" in node_ids

    def test_hierarchy_node_ids(self, multi_ec_adapter):
        node_ids = {n[0] for n in multi_ec_adapter.get_nodes()}
        # Level 1 nodes
        assert any("1.-.-.-" in nid or "1.-.-.--" in nid or nid == "ec:1.-.-.-" for nid in node_ids)

    def test_node_label(self, multi_ec_adapter):
        for node_id, label, props in multi_ec_adapter.get_nodes():
            assert label == "ec_number", f"Unexpected label {label!r} for node {node_id}"

    def test_leaf_node_has_name(self, multi_ec_adapter):
        nodes = {n[0]: n[2] for n in multi_ec_adapter.get_nodes()}
        assert nodes.get("ec:2.7.7.7", {}).get("name") == "DNA-directed DNA polymerase"

    def test_no_duplicate_nodes(self, multi_ec_adapter):
        nodes = list(multi_ec_adapter.get_nodes())
        node_ids = [n[0] for n in nodes]
        assert len(node_ids) == len(set(node_ids)), "Duplicate EC node IDs"


# ---------------------------------------------------------------------------
# TestMultiEcAnnotationAdapterEdges
# ---------------------------------------------------------------------------


class TestMultiEcAnnotationAdapterEdges:
    def test_gene_ec_edges_present(self, multi_ec_adapter):
        edges = list(multi_ec_adapter.get_edges())
        gene_ec = [e for e in edges if e[3] == "gene_catalyzes_ec_number"]
        assert len(gene_ec) == 3  # PMM0001×1 + PMM0002×2

    def test_hierarchy_edges_present(self, multi_ec_adapter):
        edges = list(multi_ec_adapter.get_edges())
        hier = [e for e in edges if e[3] == "ec_number_is_a_ec_number"]
        # MINI_EC_CACHE has 6 hierarchy edges
        assert len(hier) == 6

    def test_hierarchy_edge_structure(self, multi_ec_adapter):
        edges = list(multi_ec_adapter.get_edges())
        for e in edges:
            if e[3] == "ec_number_is_a_ec_number":
                edge_id, child, parent, label, props = e
                assert child.startswith("ec:"), f"Bad child ID: {child}"
                assert parent.startswith("ec:"), f"Bad parent ID: {parent}"

    def test_l4_to_l3_hierarchy_edge(self, multi_ec_adapter):
        edges = list(multi_ec_adapter.get_edges())
        hier = {(e[1], e[2]) for e in edges if e[3] == "ec_number_is_a_ec_number"}
        assert ("ec:2.7.7.7", "ec:2.7.7.-") in hier

    def test_gene_ec_edge_targets_known_node(self, multi_ec_adapter):
        nodes = {n[0] for n in multi_ec_adapter.get_nodes()}
        edges = list(multi_ec_adapter.get_edges())
        for e in edges:
            if e[3] == "gene_catalyzes_ec_number":
                assert e[2] in nodes, f"Target EC node {e[2]} not in node set"


# ---------------------------------------------------------------------------
# TestMultiEcCacheBehavior
# ---------------------------------------------------------------------------


class TestMultiEcCacheBehavior:
    def test_loads_from_project_cache(self, multi_ec_config_csv, tmp_ec_cache_dir):
        """EC data loads from our JSON cache without calling pypath."""
        from unittest.mock import patch

        with patch("pypath.inputs.expasy.expasy_enzymes") as mock_enz, \
             patch("pypath.inputs.expasy.expasy_enzyme_classes") as mock_cls:
            adapter = MultiEcAnnotationAdapter(
                genome_config_file=multi_ec_config_csv,
                cache_dir=tmp_ec_cache_dir,
                cache=True,
            )
            list(adapter.get_nodes())

        mock_enz.assert_not_called()
        mock_cls.assert_not_called()

    def test_cache_file_created_after_download(self, multi_ec_config_csv, tmp_path):
        """After a (mocked) download, ec_data.json is written to cache_dir."""
        from unittest.mock import patch

        cache_dir = tmp_path / "new_cache"
        cache_dir.mkdir()

        def fake_enzymes():
            return MINI_EC_CACHE["enzymes"]

        def fake_classes():
            # (level1, level2, level3, name) tuples that prepare_ec_hierarchy_dict expects
            return [
                ("1", None, None, "Oxidoreductases"),
                ("1", "1", None, "Acting on the CH-OH group of donors"),
                ("1", "1", "1", "With NAD+ as acceptor"),
                ("2", None, None, "Transferases"),
                ("2", "7", None, "Transferring phosphorus-containing groups"),
                ("2", "7", "7", "Nucleotidyltransferases"),
            ]

        with patch("pypath.inputs.expasy.expasy_enzymes", side_effect=fake_enzymes), \
             patch("pypath.inputs.expasy.expasy_enzyme_classes", side_effect=fake_classes):
            adapter = MultiEcAnnotationAdapter(
                genome_config_file=multi_ec_config_csv,
                cache_dir=cache_dir,
                cache=False,
            )

        assert (cache_dir / "ec_data.json").exists(), "ec_data.json not created"

    def test_csv_comment_lines_skipped(self, tmp_path, tmp_ec_cache_dir, tmp_ec_genome_dir):
        csv_path = tmp_path / "genomes.csv"
        csv_path.write_text(
            "ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir,clade\n"
            "# comment\n"
            f"GCF_000011465.1,Pro_MED4,59919,MED4,{tmp_ec_genome_dir},HLI\n",
            encoding="utf-8",
        )
        adapter = MultiEcAnnotationAdapter(
            genome_config_file=str(csv_path),
            cache_dir=tmp_ec_cache_dir,
            cache=True,
        )
        assert len(adapter._strain_adapters) == 1

    def test_empty_data_dir_skipped(self, tmp_path, tmp_ec_cache_dir):
        csv_path = tmp_path / "genomes.csv"
        csv_path.write_text(
            "ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir,clade\n"
            "GCF_000011465.1,Pro_MED4,59919,MED4,,HLI\n",
            encoding="utf-8",
        )
        adapter = MultiEcAnnotationAdapter(
            genome_config_file=str(csv_path),
            cache_dir=tmp_ec_cache_dir,
            cache=True,
        )
        assert len(adapter._strain_adapters) == 0
