"""Unit tests for brite_adapter.py — fixture-driven, no HTTP calls."""
import pytest

from multiomics_kg.adapters.brite_adapter import MultiBriteAdapter


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

# Transporters slice: 2 A-level entries, B and C levels, 4 KO leaves total
TRANSPORTER_TREE = {
    "name": "ko02000",
    "children": [
        {
            "name": "1. ABC Transporters",
            "children": [
                {
                    "name": "Phosphate transport system",
                    "children": [
                        {"name": "K02036  pstB; phosphate transport protein"},
                    ],
                },
                {"name": "K06147  ABC-2.B; ABC-2 type transport"},
            ],
        },
        {
            "name": "2. Non-ABC transporters",
            "children": [
                {"name": "K03455  narK; nitrate/nitrite transporter"},
                {"name": "K03284  mhpT; 3-hydroxyphenylpropionic acid transporter"},
            ],
        },
    ],
}

# Peptidase slice: 1 A-level, 1 B-level, 1 KO leaf
PEPTIDASE_TREE = {
    "name": "ko01002",
    "children": [
        {
            "name": "Serine peptidases",
            "children": [
                {"name": "K01313  F11; coagulation factor XI"},
            ],
        }
    ],
}

# Enzyme slice: 4 non-KO levels (A→B→C→D) + KO at E-level
ENZYME_TREE_DEEP = {
    "name": "ko01000",
    "children": [
        {
            "name": "1. Oxidoreductases",
            "children": [
                {
                    "name": "1.1  Acting on the CH-OH group of donors",
                    "children": [
                        {
                            "name": "1.1.1  With NAD+ or NADP+ as acceptor",
                            "children": [
                                {
                                    "name": "1.1.1.1  alcohol dehydrogenase",
                                    "children": [
                                        {"name": "K00001  E1.1.1.1, adh; alcohol dehydrogenase"},
                                    ],
                                },
                            ],
                        },
                    ],
                },
            ],
        },
    ],
}


def _make_adapter(tree_data: dict) -> MultiBriteAdapter:
    """Build adapter pre-loaded with tree data (no network calls)."""
    adapter = MultiBriteAdapter(cache_root="unused", trees=list(tree_data.keys()))
    adapter._tree_data = tree_data
    return adapter


# ---------------------------------------------------------------------------
# Node tests: basic structure
# ---------------------------------------------------------------------------

def test_get_nodes_yields_brite_category_label():
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    nodes = list(adapter.get_nodes())
    assert all(label == "brite category" for _, label, _ in nodes), (
        "All nodes must have label 'brite category'"
    )


def test_get_nodes_excludes_ko_leaves():
    """KO leaf entries (K#####...) must NOT appear as BriteCategory nodes."""
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    names = {props["name"] for _, _, props in adapter.get_nodes()}
    ko_names = [n for n in names if n.startswith("K") and n[1:6].isdigit()]
    assert not ko_names, f"KO leaves found in nodes: {ko_names}"


def test_get_nodes_a_level_has_level_zero():
    """A-level entries must be level=0 (broadest)."""
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    by_name = {props["name"]: props for _, _, props in adapter.get_nodes()}
    assert by_name["1. ABC Transporters"]["level"] == 0
    assert by_name["2. Non-ABC transporters"]["level"] == 0


def test_get_nodes_b_level_has_level_one():
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    by_name = {props["name"]: props for _, _, props in adapter.get_nodes()}
    assert by_name["Phosphate transport system"]["level"] == 1


def test_get_nodes_deep_levels():
    """Enzyme tree: A=0, B=1, C=2, D=3."""
    adapter = _make_adapter({"ko01000": ENZYME_TREE_DEEP})
    by_name = {props["name"]: props for _, _, props in adapter.get_nodes()}
    assert by_name["1. Oxidoreductases"]["level"] == 0
    assert by_name["1.1  Acting on the CH-OH group of donors"]["level"] == 1
    assert by_name["1.1.1  With NAD+ or NADP+ as acceptor"]["level"] == 2
    assert by_name["1.1.1.1  alcohol dehydrogenase"]["level"] == 3


def test_get_nodes_level_kinds():
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    by_name = {props["name"]: props for _, _, props in adapter.get_nodes()}
    assert by_name["1. ABC Transporters"]["level_kind"] == "brite_class"
    assert by_name["Phosphate transport system"]["level_kind"] == "brite_subclass"


def test_get_nodes_tree_properties():
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    by_name = {props["name"]: props for _, _, props in adapter.get_nodes()}
    abc = by_name["1. ABC Transporters"]
    assert abc["tree"] == "transporters"
    assert abc["tree_code"] == "ko02000"


# ---------------------------------------------------------------------------
# Node tests: ID format
# ---------------------------------------------------------------------------

def test_get_nodes_a_level_id_contains_tree_code():
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    by_name = {props["name"]: nid for nid, _, props in adapter.get_nodes()}
    assert "ko02000" in by_name["1. ABC Transporters"]
    assert "ko02000" in by_name["2. Non-ABC transporters"]


def test_get_nodes_a_level_ids_use_positional_suffix():
    """A-level IDs must contain .A1 and .A2 positional suffixes."""
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    by_name = {props["name"]: nid for nid, _, props in adapter.get_nodes()}
    assert ".A1" in by_name["1. ABC Transporters"]
    assert ".A2" in by_name["2. Non-ABC transporters"]


def test_get_nodes_b_level_id_embeds_parent_position():
    """B-level ID must contain the parent A-level positional index."""
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    by_name = {props["name"]: nid for nid, _, props in adapter.get_nodes()}
    # Phosphate transport system is A1.B1
    assert ".A1.B1" in by_name["Phosphate transport system"]


def test_get_nodes_id_stability():
    """Parsing the same fixture twice gives identical IDs in same order."""
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    ids_first = [nid for nid, _, _ in adapter.get_nodes()]
    ids_second = [nid for nid, _, _ in adapter.get_nodes()]
    assert ids_first == ids_second


def test_get_nodes_no_duplicate_ids():
    """No two nodes in the same tree may share an ID."""
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    ids = [nid for nid, _, _ in adapter.get_nodes()]
    assert len(ids) == len(set(ids)), f"Duplicate IDs: {[x for x in ids if ids.count(x) > 1]}"


def test_get_nodes_different_trees_different_id_namespaces():
    """Node IDs from different trees must not collide."""
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE, "ko01002": PEPTIDASE_TREE})
    ids = [nid for nid, _, _ in adapter.get_nodes()]
    assert len(ids) == len(set(ids))


# ---------------------------------------------------------------------------
# Node tests: string sanitization and HTML decoding
# ---------------------------------------------------------------------------

def test_get_nodes_string_sanitization_single_quote():
    """Single quotes in names are replaced with ^."""
    tree = {
        "ko02000": {
            "name": "ko02000",
            "children": [{"name": "It's a transporter", "children": []}],
        }
    }
    adapter = _make_adapter(tree)
    nodes = list(adapter.get_nodes())
    assert nodes[0][2]["name"] == "It^s a transporter"


def test_get_nodes_string_sanitization_pipe():
    """Pipe characters are stripped from names."""
    tree = {
        "ko02000": {
            "name": "ko02000",
            "children": [{"name": "ABC|DEF transporters", "children": []}],
        }
    }
    adapter = _make_adapter(tree)
    nodes = list(adapter.get_nodes())
    assert nodes[0][2]["name"] == "ABCDEF transporters"


def test_get_nodes_html_entity_decoding():
    """HTML entities in category names are decoded before sanitization."""
    tree = {
        "ko02000": {
            "name": "ko02000",
            "children": [
                {"name": "Transporters &amp; channels", "children": []},
                {"name": "ABC &gt; DEF", "children": []},
            ],
        }
    }
    adapter = _make_adapter(tree)
    by_name_raw = {props["name"] for _, _, props in adapter.get_nodes()}
    assert "Transporters & channels" in by_name_raw
    assert "ABC > DEF" in by_name_raw


def test_get_nodes_missing_name_skipped(caplog):
    """Entries with no 'name' key are skipped with a warning logged."""
    tree = {
        "ko02000": {
            "name": "ko02000",
            "children": [
                {"children": []},                    # no 'name' key
                {"name": "Good entry", "children": []},
            ],
        }
    }
    adapter = _make_adapter(tree)
    with caplog.at_level("WARNING"):
        nodes = list(adapter.get_nodes())
    assert len(nodes) == 1
    assert nodes[0][2]["name"] == "Good entry"
    assert any("no name" in r.message.lower() for r in caplog.records)


# ---------------------------------------------------------------------------
# Node tests: multiple trees
# ---------------------------------------------------------------------------

def test_get_nodes_multiple_trees():
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE, "ko01002": PEPTIDASE_TREE})
    nodes = list(adapter.get_nodes())
    trees_seen = {props["tree"] for _, _, props in nodes}
    assert trees_seen == {"transporters", "peptidases"}


def test_get_nodes_test_mode_caps_per_tree():
    """test_mode=True still yields some nodes (fixture is below 100 per tree)."""
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE, "ko01002": PEPTIDASE_TREE})
    adapter.test_mode = True
    nodes = list(adapter.get_nodes())
    assert len(nodes) > 0


# ---------------------------------------------------------------------------
# Edge tests: parent edges
# ---------------------------------------------------------------------------

def test_get_edges_parent_edge_label():
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    parent_edges = [e for e in adapter.get_edges() if e[3] == "brite_category_is_a_brite_category"]
    assert len(parent_edges) > 0


def test_get_edges_parent_edge_child_to_parent():
    """Brite_category_is_a_brite_category: source=child, target=parent."""
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    by_name = {props["name"]: nid for nid, _, props in adapter.get_nodes()}

    parent_edges = {
        (src, tgt)
        for _, src, tgt, label, _ in adapter.get_edges()
        if label == "brite_category_is_a_brite_category"
    }
    # Phosphate transport system (B-level) → 1. ABC Transporters (A-level)
    assert (by_name["Phosphate transport system"], by_name["1. ABC Transporters"]) in parent_edges


def test_get_edges_a_level_has_no_parent_edge():
    """A-level (level=0) nodes must NOT appear as source in parent edges."""
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    by_name = {props["name"]: nid for nid, _, props in adapter.get_nodes()}
    a_level_ids = {by_name["1. ABC Transporters"], by_name["2. Non-ABC transporters"]}

    parent_srcs = {
        src
        for _, src, tgt, label, _ in adapter.get_edges()
        if label == "brite_category_is_a_brite_category"
    }
    assert not (a_level_ids & parent_srcs), (
        "A-level nodes should not have parent edges"
    )


def test_get_edges_deep_parent_chain():
    """Enzyme tree D-level node points to C-level, not A-level."""
    adapter = _make_adapter({"ko01000": ENZYME_TREE_DEEP})
    by_name = {props["name"]: nid for nid, _, props in adapter.get_nodes()}

    parent_map = {
        src: tgt
        for _, src, tgt, label, _ in adapter.get_edges()
        if label == "brite_category_is_a_brite_category"
    }
    # D-level → C-level
    assert parent_map[by_name["1.1.1.1  alcohol dehydrogenase"]] == by_name["1.1.1  With NAD+ or NADP+ as acceptor"]
    # C-level → B-level
    assert parent_map[by_name["1.1.1  With NAD+ or NADP+ as acceptor"]] == by_name["1.1  Acting on the CH-OH group of donors"]


# ---------------------------------------------------------------------------
# Edge tests: KO edges
# ---------------------------------------------------------------------------

def test_get_edges_ko_edge_label():
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    ko_edges = [e for e in adapter.get_edges() if e[3] == "kegg_term_in_brite_category"]
    assert len(ko_edges) >= 4  # fixture has 4 KO leaves


def test_get_edges_ko_source_is_kegg_orthology():
    """KO edge source IDs must use kegg.orthology: prefix."""
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    ko_srcs = [
        src
        for _, src, tgt, label, _ in adapter.get_edges()
        if label == "kegg_term_in_brite_category"
    ]
    assert all("kegg.orthology:" in src for src in ko_srcs), (
        f"Non-kegg.orthology source found: {[s for s in ko_srcs if 'kegg.orthology:' not in s]}"
    )


def test_get_edges_ko_target_is_immediate_parent():
    """KO leaf at B-level of TRANSPORTER_TREE connects to its B-level parent, not A-level."""
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    by_name = {props["name"]: nid for nid, _, props in adapter.get_nodes()}

    ko_edges = {
        src: tgt
        for _, src, tgt, label, _ in adapter.get_edges()
        if label == "kegg_term_in_brite_category"
    }
    # K02036 (pstB) lives under "Phosphate transport system" (B-level, A1.B1)
    phosphate_id = by_name["Phosphate transport system"]
    assert any(tgt == phosphate_id for tgt in ko_edges.values()), (
        "K02036 should point to 'Phosphate transport system' node"
    )


def test_get_edges_ko_id_extracted_from_name():
    """K-ID must be extracted as the first token of the name field."""
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    ko_srcs = [
        src
        for _, src, tgt, label, _ in adapter.get_edges()
        if label == "kegg_term_in_brite_category"
    ]
    # K02036, K06147, K03455, K03284 must all appear
    expected_kos = {"K02036", "K06147", "K03455", "K03284"}
    found_kos = {src.split(":")[-1] for src in ko_srcs}
    assert expected_kos == found_kos, f"Expected {expected_kos}, got {found_kos}"


def test_get_edges_no_duplicate_ko_edges():
    """Same (KO, BriteCategory) pair must not appear twice."""
    adapter = _make_adapter({"ko02000": TRANSPORTER_TREE})
    ko_pairs = [
        (src, tgt)
        for _, src, tgt, label, _ in adapter.get_edges()
        if label == "kegg_term_in_brite_category"
    ]
    assert len(ko_pairs) == len(set(ko_pairs)), "Duplicate KO→BRITE edges found"


def test_get_edges_ko_in_deep_tree():
    """Enzyme tree: KO leaf at E-level points to D-level category."""
    adapter = _make_adapter({"ko01000": ENZYME_TREE_DEEP})
    by_name = {props["name"]: nid for nid, _, props in adapter.get_nodes()}

    ko_edges = {
        src: tgt
        for _, src, tgt, label, _ in adapter.get_edges()
        if label == "kegg_term_in_brite_category"
    }
    d_level_id = by_name["1.1.1.1  alcohol dehydrogenase"]
    assert "kegg.orthology:K00001" in ko_edges
    assert ko_edges["kegg.orthology:K00001"] == d_level_id
