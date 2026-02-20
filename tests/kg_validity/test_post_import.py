"""
Post-import Cypher validation tests.

The post-import.cypher script runs after neo4j-admin import to:
1. Create bidirectional Gene_is_homolog_of_gene edges between genes sharing
   a Cyanorak_cluster (both A→B and B→A, with distance property)
2. Propagate Affects_expression_of_homolog edges: if X affects gene A,
   and A is homolog of B, then X affects_homolog B

These tests verify that post-import ran correctly and produced valid results.
"""

import pytest


pytestmark = pytest.mark.kg


# ---------------------------------------------------------------------------
# Gene_is_homolog_of_gene edges
# ---------------------------------------------------------------------------

def test_homolog_edges_exist(run_query):
    """Gene_is_homolog_of_gene edges must exist after post-import."""
    result = run_query(
        "MATCH ()-[r:Gene_is_homolog_of_gene]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt > 0, (
        "No Gene_is_homolog_of_gene edges found. "
        "Did post-import.cypher run successfully?"
    )


def test_homolog_edge_count_reasonable(run_query):
    """
    With 12 strains and thousands of shared clusters, expect millions of
    homolog edge pairs. A low count indicates post-import failure.
    """
    result = run_query(
        "MATCH ()-[r:Gene_is_homolog_of_gene]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt > 100_000, (
        f"Only {cnt} Gene_is_homolog_of_gene edges found; "
        f"expected > 100,000 for 12 strains with shared Cyanorak clusters"
    )


def test_homolog_edges_are_bidirectional(run_query):
    """
    post-import.cypher creates both A→B and B→A homolog edges so LLM agents
    can traverse from any direction. For a sample of 200 edges, verify that
    the reverse edge always exists.
    """
    result = run_query("""
        MATCH (g1:Gene)-[:Gene_is_homolog_of_gene]->(g2:Gene)
        WITH g1, g2 LIMIT 200
        WHERE NOT (g2)-[:Gene_is_homolog_of_gene]->(g1)
        RETURN count(*) AS missing_reverse
    """)
    missing = result[0]["missing_reverse"]
    assert missing == 0, (
        f"{missing} forward homolog edge(s) have no corresponding reverse edge. "
        f"post-import.cypher bidirectional MERGE may have failed."
    )


def test_homolog_edges_have_distance_property(run_query):
    """Every Gene_is_homolog_of_gene edge must carry a distance property."""
    result = run_query("""
        MATCH ()-[r:Gene_is_homolog_of_gene]->()
        WHERE r.distance IS NULL
        RETURN count(r) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Gene_is_homolog_of_gene edges are missing the distance property"
    )


def test_homolog_edges_have_cluster_id(run_query):
    """Every homolog edge must reference the Cyanorak cluster that caused it."""
    result = run_query("""
        MATCH ()-[r:Gene_is_homolog_of_gene]->()
        WHERE r.cluster_id IS NULL
        RETURN count(r) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Gene_is_homolog_of_gene edges are missing cluster_id"
    )


def test_homolog_edges_source_is_cyanorak(run_query):
    """All homolog edges created by post-import should have source='cyanorak_cluster'."""
    result = run_query("""
        MATCH ()-[r:Gene_is_homolog_of_gene]->()
        WHERE r.source IS NULL OR r.source <> 'cyanorak_cluster'
        RETURN count(r) AS bad_rows,
               collect(DISTINCT r.source)[..5] AS bad_sources
    """)
    bad = result[0]["bad_rows"]
    assert bad == 0, (
        f"{bad} Gene_is_homolog_of_gene edges have unexpected source: "
        f"{result[0]['bad_sources']}"
    )


def test_homolog_distance_valid_values(run_query):
    """
    distance values are produced by a CASE statement in post-import.cypher.
    All values must start with one of the expected prefixes.
    """
    valid_prefixes = [
        "same strain:",
        "same clade:",
        "same species:",
        "same genus:",
        "same order:",
        "cross order",
    ]
    result = run_query("""
        MATCH ()-[r:Gene_is_homolog_of_gene]->()
        WHERE r.distance IS NOT NULL
        RETURN collect(DISTINCT r.distance) AS all_distances
    """)
    if not result or not result[0]["all_distances"]:
        pytest.skip("No distance values to check")

    all_distances = result[0]["all_distances"]
    bad = [
        d for d in all_distances
        if not any(d.startswith(p) for p in valid_prefixes)
    ]
    assert not bad, (
        f"Unexpected distance values found: {bad}\n"
        f"Valid prefixes: {valid_prefixes}"
    )


def test_no_self_homolog_edges(run_query):
    """A gene should not be its own homolog (same node → same node)."""
    result = run_query("""
        MATCH (g:Gene)-[:Gene_is_homolog_of_gene]->(g)
        RETURN count(g) AS self_loops
    """)
    self_loops = result[0]["self_loops"]
    assert self_loops == 0, (
        f"{self_loops} Gene node(s) have self-referential Gene_is_homolog_of_gene edges"
    )


# ---------------------------------------------------------------------------
# Affects_expression_of_homolog edges
# ---------------------------------------------------------------------------

def test_expression_homolog_edges_exist(run_query):
    """Affects_expression_of_homolog edges must exist after post-import."""
    result = run_query(
        "MATCH ()-[r:Affects_expression_of_homolog]->() RETURN count(r) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt > 0, (
        "No Affects_expression_of_homolog edges found. "
        "Did post-import.cypher propagation step run?"
    )


def test_expression_homolog_has_original_gene(run_query):
    """Every Affects_expression_of_homolog must record the original_gene ID."""
    result = run_query("""
        MATCH ()-[r:Affects_expression_of_homolog]->()
        WHERE r.original_gene IS NULL
        RETURN count(r) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Affects_expression_of_homolog edges are missing original_gene"
    )


def test_expression_homolog_has_homology_metadata(run_query):
    """
    Affects_expression_of_homolog edges must carry homology_source,
    homology_cluster_id, and distance (copied from the homolog edge).
    """
    result = run_query("""
        MATCH ()-[r:Affects_expression_of_homolog]->()
        WHERE r.homology_source IS NULL
           OR r.homology_cluster_id IS NULL
           OR r.distance IS NULL
        RETURN count(r) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Affects_expression_of_homolog edges are missing homology metadata "
        f"(homology_source / homology_cluster_id / distance)"
    )


def test_expression_homolog_inherits_log2fc(run_query):
    """
    Expression homolog edges must carry log2_fold_change (always present on the
    original edge). Note: adjusted_p_value may be null if the original study
    did not report it — that is expected and not tested here.
    """
    result = run_query("""
        MATCH ()-[r:Affects_expression_of_homolog]->()
        WHERE r.log2_fold_change IS NULL
        RETURN count(r) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Affects_expression_of_homolog edges are missing log2_fold_change"
    )


def test_expression_homolog_count_vs_direct(run_query):
    """
    The number of homolog expression edges should be substantially larger
    than the number of direct expression edges (homology amplifies coverage).
    If homolog_count < direct_count, propagation likely failed.
    """
    result = run_query("""
        MATCH ()-[:Affects_expression_of]->()
        WITH count(*) AS direct
        MATCH ()-[:Affects_expression_of_homolog]->()
        WITH direct, count(*) AS homolog
        RETURN direct, homolog
    """)
    row = result[0]
    assert row["homolog"] > row["direct"], (
        f"Affects_expression_of_homolog count ({row['homolog']}) is not greater "
        f"than direct Affects_expression_of count ({row['direct']}). "
        f"Homolog propagation may have failed or graph has very few shared clusters."
    )
