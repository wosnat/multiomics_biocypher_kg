"""Unit tests for the pure InterPro reference parser (prepare_data step 9 core)."""

from __future__ import annotations

from multiomics_kg.utils.interpro_reference import (
    build_reference,
    normalize_type,
    parse_entry_list,
    parse_parent_child_tree,
)

ENTRY_LIST = "\n".join(
    [
        "ENTRY_AC\tENTRY_TYPE\tENTRY_NAME",
        "IPR000001\tDomain\tKringle",
        "IPR000002\tFamily\tCdc20/Fizzy",
        "IPR000003\tHomologous_superfamily\tRetinoid X receptor",
        "IPR000004\tConserved_site\tSome conserved site",
        "IPR000685\tFamily\tRuBisCO large subunit",
        "malformed line without tabs",
    ]
)

# Depth = leading "--" pairs. IPR000002 is a child of IPR000001; IPR000004 a
# grandchild. IPR000685 is a separate root. IPR000003 is not in the tree.
PARENT_CHILD = "\n".join(
    [
        "IPR000001::Kringle::",
        "--IPR000002::Cdc20/Fizzy::",
        "----IPR000004::Some conserved site::",
        "IPR000685::RuBisCO large subunit::",
    ]
)


def test_normalize_type_matches_callsjson_upper():
    assert normalize_type("Homologous_superfamily") == "HOMOLOGOUS_SUPERFAMILY"
    assert normalize_type("Conserved_site") == "CONSERVED_SITE"
    assert normalize_type("Family") == "FAMILY"
    assert normalize_type(None) is None


def test_parse_entry_list_skips_header_and_malformed():
    entries = parse_entry_list(ENTRY_LIST)
    assert set(entries) == {"IPR000001", "IPR000002", "IPR000003", "IPR000004", "IPR000685"}
    assert entries["IPR000001"] == {"type": "DOMAIN", "name": "Kringle"}
    assert entries["IPR000003"]["type"] == "HOMOLOGOUS_SUPERFAMILY"


def test_parse_parent_child_tree_depth_and_parent():
    tree = parse_parent_child_tree(PARENT_CHILD)
    assert tree["IPR000001"] == {"parent": None, "level": 0}
    assert tree["IPR000002"] == {"parent": "IPR000001", "level": 1}
    assert tree["IPR000004"] == {"parent": "IPR000002", "level": 2}
    assert tree["IPR000685"] == {"parent": None, "level": 0}


def test_build_reference_combines_and_defaults():
    ref = build_reference(ENTRY_LIST, PARENT_CHILD)
    # name/type from entry.list; parent/level from tree
    assert ref["IPR000002"] == {
        "name": "Cdc20/Fizzy",
        "type": "FAMILY",
        "parent": "IPR000001",
        "level": 1,
    }
    # entry not in the tree → parentless root at level 0
    assert ref["IPR000003"]["parent"] is None
    assert ref["IPR000003"]["level"] == 0
    assert ref["IPR000003"]["type"] == "HOMOLOGOUS_SUPERFAMILY"


def test_build_reference_hierarchy_never_dangles():
    """Every entry's parent (when set) is itself a key in the reference."""
    ref = build_reference(ENTRY_LIST, PARENT_CHILD)
    for acc, meta in ref.items():
        if meta["parent"] is not None:
            assert meta["parent"] in ref, f"{acc} parent {meta['parent']} missing"
