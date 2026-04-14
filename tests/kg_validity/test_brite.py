"""
KG validity tests for KEGG BRITE functional hierarchy.

Validates:
- All 12 configured BRITE trees are present
- Level=0 (A-level) nodes exist and have no parent edges
- Non-root nodes have exactly one parent edge
- KO→BRITE edges exist and are sufficiently populated
- Specific biological spot checks (ftsH in peptidases, pstB in transporters)
- Post-import computed properties are populated
- No duplicate KO→BRITE edges
"""

import pytest


pytestmark = pytest.mark.kg


# ---------------------------------------------------------------------------
# Tree presence
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("tree_code,tree_name", [
    ("ko01000", "enzymes"),
    ("ko02000", "transporters"),
    ("ko01002", "peptidases"),
    ("ko03000", "transcription_factors"),
    ("ko02044", "secretion"),
    ("ko02022", "two_component"),
    ("ko02048", "defense"),
    ("ko03110", "chaperones"),
    ("ko03011", "ribosome"),
    ("ko03012", "translation_factors"),
    ("ko03016", "trna_biogenesis"),
    ("ko03032", "dna_replication"),
])
def test_brite_tree_present(run_query, tree_code, tree_name):
    """Configured BRITE trees should have BriteCategory nodes **iff** at least
    one of their KOs is a KeggTerm in the KG.

    After subhierarchy pruning (see plans/2026-04-14-brite-pruning.md), a tree
    whose KOs are all outside the KG legitimately prunes to zero nodes. We assert
    consistency: either the tree has nodes AND has at least one incident KO
    edge, or it has neither.
    """
    node_result = run_query(
        "MATCH (b:BriteCategory {tree_code: $tc}) RETURN count(b) AS cnt",
        tc=tree_code,
    )
    edge_result = run_query(
        """
        MATCH (:KeggTerm)-[:Kegg_term_in_brite_category]->(b:BriteCategory {tree_code: $tc})
        RETURN count(*) AS cnt
        """,
        tc=tree_code,
    )
    node_count = node_result[0]["cnt"]
    edge_count = edge_result[0]["cnt"]

    if node_count == 0:
        assert edge_count == 0, (
            f"Tree {tree_code} ({tree_name}) has 0 nodes but {edge_count} "
            f"KO edges — pruning inconsistency"
        )
    else:
        assert edge_count > 0, (
            f"Tree {tree_code} ({tree_name}) has {node_count} nodes but 0 "
            f"KO edges — pruned tree should always have KO coverage"
        )


# ---------------------------------------------------------------------------
# Level semantics
# ---------------------------------------------------------------------------

def test_brite_level_zero_nodes_exist(run_query):
    """At least one level-0 node must exist per non-empty tree.

    With subhierarchy pruning, some of the 12 configured trees may prune to
    zero nodes if none of their KOs are in the KG. We require at least 6
    level-0 nodes overall (half the trees must remain populated for the
    config to be meaningful).
    """
    result = run_query(
        """
        MATCH (b:BriteCategory {level: 0})
        RETURN count(b) AS cnt
        """
    )
    assert result[0]["cnt"] >= 6, (
        f"Only {result[0]['cnt']} level=0 BriteCategory nodes; expected ≥ 6 "
        f"(half the 12 configured trees)"
    )


def test_brite_level_zero_has_no_parent_edge(run_query):
    """Level=0 nodes must NOT be the source of a Brite_category_is_a_brite_category edge."""
    result = run_query("""
        MATCH (b:BriteCategory {level: 0})-[:Brite_category_is_a_brite_category]->()
        RETURN count(b) AS cnt
    """)
    assert result[0]["cnt"] == 0, (
        f"{result[0]['cnt']} level=0 BriteCategory nodes have parent edges (should be 0)"
    )


def test_brite_non_root_has_exactly_one_parent(run_query):
    """Every non-root BriteCategory (level > 0) must have exactly one parent edge."""
    result = run_query("""
        MATCH (b:BriteCategory)
        WHERE b.level > 0
        WITH b, size([(b)-[:Brite_category_is_a_brite_category]->() | 1]) AS parent_count
        WHERE parent_count <> 1
        RETURN count(b) AS bad, collect(b.name)[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} non-root BriteCategory nodes have ≠ 1 parent: {result[0]['examples']}"
    )


# ---------------------------------------------------------------------------
# KO edges
# ---------------------------------------------------------------------------

def test_brite_ko_edges_exist(run_query):
    """Kegg_term_in_brite_category edges must exist in significant numbers."""
    result = run_query("""
        MATCH ()-[r:Kegg_term_in_brite_category]->()
        RETURN count(r) AS cnt
    """)
    assert result[0]["cnt"] >= 200, (
        f"Only {result[0]['cnt']} Kegg_term_in_brite_category edges; expected ≥ 200"
    )


def test_brite_transporters_ko_count(run_query):
    """Transporters tree (ko02000) must have at least 50 KO leaf edges."""
    result = run_query("""
        MATCH (ko:KeggTerm)-[:Kegg_term_in_brite_category]->(b:BriteCategory {tree_code: 'ko02000'})
        RETURN count(DISTINCT ko) AS cnt
    """)
    assert result[0]["cnt"] >= 50, (
        f"Transporters tree has only {result[0]['cnt']} KO edges; expected ≥ 50"
    )


def test_brite_peptidases_ko_count(run_query):
    """Peptidases tree (ko01002) must have at least 20 KO leaf edges."""
    result = run_query("""
        MATCH (ko:KeggTerm)-[:Kegg_term_in_brite_category]->(b:BriteCategory {tree_code: 'ko01002'})
        RETURN count(DISTINCT ko) AS cnt
    """)
    assert result[0]["cnt"] >= 20, (
        f"Peptidases tree has only {result[0]['cnt']} KO edges; expected ≥ 20"
    )


def test_brite_no_duplicate_ko_edges(run_query):
    """No (KeggTerm, BriteCategory) pair should appear more than once."""
    result = run_query("""
        MATCH (ko:KeggTerm)-[r:Kegg_term_in_brite_category]->(b:BriteCategory)
        WITH ko, b, count(r) AS cnt
        WHERE cnt > 1
        RETURN count(*) AS duplicates
    """)
    assert result[0]["duplicates"] == 0, (
        f"{result[0]['duplicates']} duplicate Kegg_term_in_brite_category (KO,BriteCategory) pairs"
    )


def test_brite_no_dangling_ko_edges(run_query):
    """No Kegg_term_in_brite_category edge may start from a non-KeggTerm node.

    Regression guard for the pruning fix: before pruning, the BRITE adapter
    emitted ~12K edges referring to KO IDs that had no KeggTerm node, causing
    ``neo4j-admin import`` to fail with exit code 70.
    """
    result = run_query(
        """
        MATCH (n)-[r:Kegg_term_in_brite_category]->()
        WHERE NOT n:KeggTerm
        RETURN count(r) AS cnt
        """
    )
    assert result[0]["cnt"] == 0, (
        f"{result[0]['cnt']} Kegg_term_in_brite_category edges have a "
        f"non-KeggTerm source — pruning regression"
    )


# ---------------------------------------------------------------------------
# Biological spot checks
# ---------------------------------------------------------------------------

def test_brite_ftsh_in_peptidases(run_query):
    """ftsH (K03798) must appear in the peptidases tree (ko01002)."""
    result = run_query("""
        MATCH (ko:KeggTerm)-[:Kegg_term_in_brite_category]->(b:BriteCategory {tree_code: 'ko01002'})
        WHERE ko.name STARTS WITH 'K03798'
           OR ko.name CONTAINS 'ftsH'
        RETURN count(ko) AS cnt
    """)
    assert result[0]["cnt"] > 0, (
        "ftsH (K03798) not found in peptidases tree (ko01002)"
    )


def test_brite_pstb_in_transporters(run_query):
    """pstB (K02036) must appear in the transporters tree (ko02000)."""
    result = run_query("""
        MATCH (ko:KeggTerm {level_kind: 'ko'})-[:Kegg_term_in_brite_category]->(b:BriteCategory {tree_code: 'ko02000'})
        WHERE ko.name STARTS WITH 'K02036'
           OR ko.name CONTAINS 'pstB'
        RETURN count(ko) AS cnt
    """)
    assert result[0]["cnt"] > 0, (
        "pstB (K02036) not found in transporters tree (ko02000)"
    )


# ---------------------------------------------------------------------------
# Post-import computed properties
# ---------------------------------------------------------------------------

def test_brite_member_ko_count_populated(run_query):
    """member_ko_count must be set and > 0 on level=0 BriteCategory nodes."""
    result = run_query("""
        MATCH (b:BriteCategory {level: 0})
        WHERE b.member_ko_count IS NULL OR b.member_ko_count = 0
        RETURN count(b) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} level=0 BriteCategory nodes have member_ko_count NULL or 0"
    )


def test_brite_gene_count_populated(run_query):
    """gene_count must be set on level=0 BriteCategory nodes (may be 0 for rare trees)."""
    result = run_query("""
        MATCH (b:BriteCategory {level: 0})
        WHERE b.gene_count IS NULL
        RETURN count(b) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} level=0 BriteCategory nodes have gene_count IS NULL"
    )
