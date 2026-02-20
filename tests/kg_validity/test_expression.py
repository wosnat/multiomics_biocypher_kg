"""
Expression edge data quality tests for the multi-omics knowledge graph.

Validates Affects_expression_of edges, which carry quantitative data
from 19+ differential expression studies. Bad data types or invalid
values would silently break downstream LLM analysis.

Checks:
- log2_fold_change and adjusted_p_value are stored as floats (not strings)
- adjusted_p_value is in [0, 1]
- expression_direction is only 'up' or 'down'
- Required properties (publications, control_condition) are present
- At least one expression edge exists per organism/environmental_condition source
"""

import pytest


pytestmark = pytest.mark.kg


# ---------------------------------------------------------------------------
# Edge existence
# ---------------------------------------------------------------------------

def test_expression_edges_exist(run_query):
    """At least one Affects_expression_of edge must exist."""
    result = run_query(
        "MATCH ()-[e:Affects_expression_of]->() RETURN count(e) AS cnt"
    )
    assert result[0]["cnt"] > 0, "No Affects_expression_of edges found"


def test_expression_edge_count_minimum(run_query):
    """
    With 19+ studies across multiple strains, expect tens of thousands of edges.
    A low count indicates an import failure.
    """
    result = run_query(
        "MATCH ()-[e:Affects_expression_of]->() RETURN count(e) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt > 10_000, (
        f"Only {cnt} Affects_expression_of edges found; expected > 10,000 "
        f"from 19+ differential expression studies"
    )


# ---------------------------------------------------------------------------
# Numeric data types (type correctness)
# ---------------------------------------------------------------------------

def test_log2_fold_change_is_numeric(run_query):
    """
    log2_fold_change must be stored as a float.
    BioCypher imports floats correctly; this catches regressions where
    the adapter emits strings (e.g., '1.23' instead of 1.23).

    In Cypher: toFloat(x) returns null if x is already a float (no-op in Neo4j 4.4),
    so we check that the property type is numeric via size() failure heuristic.
    We match edges where log2_fold_change is not null and verify it's usable
    in arithmetic — a string would raise a type error.
    """
    result = run_query("""
        MATCH ()-[e:Affects_expression_of]->()
        WHERE e.log2_fold_change IS NOT NULL
        WITH e.log2_fold_change AS fc
        WHERE fc <> toFloat(toString(fc))
        RETURN count(*) AS bad_rows
    """)
    assert result[0]["bad_rows"] == 0, (
        f"{result[0]['bad_rows']} expression edges have non-numeric log2_fold_change"
    )


def test_adjusted_p_value_is_numeric(run_query):
    """adjusted_p_value must be stored as a float, not a string."""
    result = run_query("""
        MATCH ()-[e:Affects_expression_of]->()
        WHERE e.adjusted_p_value IS NOT NULL
        WITH e.adjusted_p_value AS padj
        WHERE padj <> toFloat(toString(padj))
        RETURN count(*) AS bad_rows
    """)
    assert result[0]["bad_rows"] == 0, (
        f"{result[0]['bad_rows']} expression edges have non-numeric adjusted_p_value"
    )


# ---------------------------------------------------------------------------
# Value range constraints
# ---------------------------------------------------------------------------

def test_adjusted_p_value_in_valid_range(run_query):
    """adjusted_p_value must be in [0, 1] — it is a probability."""
    result = run_query("""
        MATCH ()-[e:Affects_expression_of]->()
        WHERE e.adjusted_p_value IS NOT NULL
          AND (e.adjusted_p_value < 0 OR e.adjusted_p_value > 1)
        RETURN count(e) AS bad_rows,
               collect(e.adjusted_p_value)[..5] AS examples
    """)
    bad = result[0]["bad_rows"]
    assert bad == 0, (
        f"{bad} expression edges have adjusted_p_value outside [0, 1]. "
        f"Example values: {result[0]['examples']}"
    )


def test_expression_direction_valid_values(run_query):
    """
    expression_direction must be 'up' or 'down'.
    Any other value (e.g., 'UP', 'upregulated', '') indicates a parsing bug.
    """
    result = run_query("""
        MATCH ()-[e:Affects_expression_of]->()
        WHERE e.expression_direction IS NOT NULL
          AND NOT e.expression_direction IN ['up', 'down']
        RETURN count(e) AS bad_rows,
               collect(DISTINCT e.expression_direction)[..10] AS bad_values
    """)
    bad = result[0]["bad_rows"]
    assert bad == 0, (
        f"{bad} expression edges have invalid expression_direction. "
        f"Invalid values seen: {result[0]['bad_values']}"
    )


def test_log2_fold_change_direction_consistency(run_query):
    """
    Sanity check: if expression_direction='up', log2_fold_change should be > 0,
    and 'down' should be < 0. A sign flip indicates a parsing error.
    """
    result = run_query("""
        MATCH ()-[e:Affects_expression_of]->()
        WHERE e.expression_direction IS NOT NULL
          AND e.log2_fold_change IS NOT NULL
          AND (
            (e.expression_direction = 'up'   AND e.log2_fold_change < 0)
            OR
            (e.expression_direction = 'down' AND e.log2_fold_change > 0)
          )
        RETURN count(e) AS inconsistent,
               collect(e.expression_direction + ':' + toString(e.log2_fold_change))[..5]
               AS examples
    """)
    bad = result[0]["inconsistent"]
    assert bad == 0, (
        f"{bad} expression edges have mismatched direction vs log2FC sign. "
        f"Examples: {result[0]['examples']}"
    )


# ---------------------------------------------------------------------------
# Required property presence
# ---------------------------------------------------------------------------

def test_expression_edges_have_publications(run_query):
    """Every expression edge must reference at least one publication."""
    result = run_query("""
        MATCH ()-[e:Affects_expression_of]->()
        WHERE e.publications IS NULL OR size(e.publications) = 0
        RETURN count(e) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Affects_expression_of edges are missing the publications property"
    )


def test_expression_edges_have_control_condition(run_query):
    """Every expression edge must have a control_condition (baseline description)."""
    result = run_query("""
        MATCH ()-[e:Affects_expression_of]->()
        WHERE e.control_condition IS NULL OR e.control_condition = ''
        RETURN count(e) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Affects_expression_of edges are missing control_condition"
    )


def test_expression_edges_have_direction(run_query):
    """Every expression edge must have expression_direction set."""
    result = run_query("""
        MATCH ()-[e:Affects_expression_of]->()
        WHERE e.expression_direction IS NULL OR e.expression_direction = ''
        RETURN count(e) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Affects_expression_of edges are missing expression_direction"
    )


# ---------------------------------------------------------------------------
# Source node coverage
# ---------------------------------------------------------------------------

def test_organisms_have_expression_edges(run_query):
    """
    At least some OrganismTaxon nodes must be expression sources
    (coculture experiments use organisms as treatments).
    """
    result = run_query("""
        MATCH (o:OrganismTaxon)-[:Affects_expression_of]->()
        RETURN count(DISTINCT o) AS cnt
    """)
    assert result[0]["cnt"] > 0, (
        "No OrganismTaxon node is a source of Affects_expression_of; "
        "coculture experiment edges may be missing"
    )


def test_environmental_conditions_have_expression_edges(run_query):
    """
    EnvironmentalCondition nodes must source expression edges
    (stress experiments use conditions as treatments).
    """
    result = run_query("""
        MATCH (e:EnvironmentalCondition)-[:Affects_expression_of]->()
        RETURN count(DISTINCT e) AS cnt
    """)
    assert result[0]["cnt"] > 0, (
        "No EnvironmentalCondition node is a source of Affects_expression_of; "
        "stress experiment edges may be missing"
    )


def test_expression_edges_target_genes(run_query):
    """All Affects_expression_of edges must target a Gene node."""
    result = run_query("""
        MATCH ()-[e:Affects_expression_of]->(target)
        WHERE NOT target:Gene
        RETURN count(e) AS wrong_targets,
               collect(DISTINCT labels(target))[..5] AS target_labels
    """)
    bad = result[0]["wrong_targets"]
    assert bad == 0, (
        f"{bad} Affects_expression_of edges point to non-Gene nodes: "
        f"{result[0]['target_labels']}"
    )
