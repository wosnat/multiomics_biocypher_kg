"""
Unit tests for compute_go_levels (GO BFS level assignment).

Uses synthetic go_data dicts — does not hit the real OBO file or KG.
"""
from multiomics_kg.utils.go_utils import (
    CANONICAL_GO_ROOTS,
    compute_go_levels,
)


def _bp(parents):
    """Helper: build a biological_process term dict entry."""
    return {"name": "", "namespace": "biological_process", "parents": parents}


def test_root_only_gets_level_zero():
    go_data = {"GO:0008150": _bp([])}
    levels, orphans = compute_go_levels(go_data)
    assert levels["GO:0008150"] == (0, False)
    assert orphans == []


def test_linear_chain_assigns_monotonic_levels():
    go_data = {
        "GO:0008150": _bp([]),                       # root
        "GO:0000001": _bp([["GO:0008150", "is_a"]]),  # depth 1
        "GO:0000002": _bp([["GO:0000001", "is_a"]]),  # depth 2
    }
    levels, orphans = compute_go_levels(go_data)
    assert levels["GO:0008150"] == (0, False)
    assert levels["GO:0000001"] == (1, False)
    assert levels["GO:0000002"] == (2, False)
    assert orphans == []


def test_diamond_equal_arms_not_best_effort():
    # Two parents both at depth 1 → confluence node at depth 2 with min==max.
    go_data = {
        "GO:0008150": _bp([]),
        "GO:0000001": _bp([["GO:0008150", "is_a"]]),
        "GO:0000002": _bp([["GO:0008150", "is_a"]]),
        "GO:0000003": _bp([["GO:0000001", "is_a"], ["GO:0000002", "is_a"]]),
    }
    levels, _ = compute_go_levels(go_data)
    assert levels["GO:0000003"] == (2, False)


def test_diamond_unequal_arms_flags_best_effort():
    # One parent at depth 1, another at depth 2 → min=2, max=3, flagged.
    go_data = {
        "GO:0008150": _bp([]),
        "GO:0000001": _bp([["GO:0008150", "is_a"]]),        # depth 1
        "GO:0000002": _bp([["GO:0000001", "is_a"]]),         # depth 2
        "GO:0000003": _bp([["GO:0008150", "is_a"], ["GO:0000002", "is_a"]]),
    }
    levels, _ = compute_go_levels(go_data)
    assert levels["GO:0000003"] == (1, True)  # min=1 via direct root, max=3 via chain


def test_part_of_is_traversed():
    go_data = {
        "GO:0008150": _bp([]),
        "GO:0000001": _bp([["GO:0008150", "part_of"]]),
    }
    levels, _ = compute_go_levels(go_data)
    assert levels["GO:0000001"] == (1, False)


def test_regulates_is_not_traversed():
    # A term whose only parent is via 'regulates' must not be reachable
    # from a root → orphan.
    go_data = {
        "GO:0008150": _bp([]),
        "GO:0000001": _bp([["GO:0008150", "regulates"]]),
    }
    levels, orphans = compute_go_levels(go_data)
    assert "GO:0000001" not in levels
    assert "GO:0000001" in orphans


def test_orphan_not_canonical_root_is_surfaced():
    # Term with empty parents that isn't one of the three canonical roots.
    go_data = {
        "GO:0008150": _bp([]),
        "GO:0999999": _bp([]),
    }
    levels, orphans = compute_go_levels(go_data)
    assert "GO:0999999" not in levels
    assert orphans == ["GO:0999999"]


def test_all_three_namespaces_handled():
    go_data = {
        "GO:0008150": {"name": "", "namespace": "biological_process", "parents": []},
        "GO:0003674": {"name": "", "namespace": "molecular_function", "parents": []},
        "GO:0005575": {"name": "", "namespace": "cellular_component", "parents": []},
    }
    levels, orphans = compute_go_levels(go_data)
    assert levels["GO:0008150"] == (0, False)
    assert levels["GO:0003674"] == (0, False)
    assert levels["GO:0005575"] == (0, False)
    assert orphans == []


def test_canonical_roots_constant():
    assert CANONICAL_GO_ROOTS == frozenset({
        "GO:0008150",  # biological_process
        "GO:0003674",  # molecular_function
        "GO:0005575",  # cellular_component
    })
