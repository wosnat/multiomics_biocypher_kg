"""
Expression edge data quality tests for the multi-omics knowledge graph.

Validates Changes_expression_of edges (from Experiment nodes to Gene nodes),
which carry quantitative data from 19+ differential expression studies.
Bad data types or invalid values would silently break downstream LLM analysis.

Checks:
- log2_fold_change and adjusted_p_value are stored as floats (not strings)
- adjusted_p_value is in [0, 1]
- expression_direction is only 'up' or 'down'
- Experiment nodes have required properties (control, name)
- Every Experiment is linked to a Publication via Has_experiment
- All expression edges source from Experiment nodes and target Gene nodes
"""

import pytest


pytestmark = pytest.mark.kg


# ---------------------------------------------------------------------------
# Edge existence
# ---------------------------------------------------------------------------

def test_expression_edges_exist(run_query):
    """At least one Changes_expression_of edge must exist."""
    result = run_query(
        "MATCH ()-[e:Changes_expression_of]->() RETURN count(e) AS cnt"
    )
    assert result[0]["cnt"] > 0, "No Changes_expression_of edges found"


def test_expression_edge_count_minimum(run_query):
    """
    With 19+ studies across multiple strains, expect tens of thousands of edges.
    A low count indicates an import failure.
    """
    result = run_query(
        "MATCH ()-[e:Changes_expression_of]->() RETURN count(e) AS cnt"
    )
    cnt = result[0]["cnt"]
    assert cnt > 10_000, (
        f"Only {cnt} expression edges found; expected > 10,000 "
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
    """
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
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
        MATCH ()-[e:Changes_expression_of]->()
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
        MATCH ()-[e:Changes_expression_of]->()
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
        MATCH ()-[e:Changes_expression_of]->()
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
        MATCH ()-[e:Changes_expression_of]->()
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
# Required property presence (on Experiment nodes, not edges)
# ---------------------------------------------------------------------------

def test_expression_edges_have_publications(run_query):
    """Every Experiment must be linked to a Publication via Has_experiment."""
    result = run_query("""
        MATCH (exp:Experiment)
        WHERE NOT (exp)<-[:Has_experiment]-(:Publication)
        RETURN count(exp) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Experiment nodes have no Has_experiment edge from a Publication"
    )


def test_experiments_have_control(run_query):
    """Every Experiment node must have a control property (baseline description)."""
    result = run_query("""
        MATCH (exp:Experiment)
        WHERE exp.control IS NULL OR exp.control = ''
        RETURN count(exp) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Experiment nodes are missing the control property"
    )


def test_expression_edges_have_direction(run_query):
    """Every expression edge must have expression_direction set."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE e.expression_direction IS NULL OR e.expression_direction = ''
        RETURN count(e) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} expression edges are missing expression_direction"
    )


# ---------------------------------------------------------------------------
# Source/target node coverage
# ---------------------------------------------------------------------------

def test_experiments_source_expression_edges(run_query):
    """
    Experiment nodes must be the source of Changes_expression_of edges.
    """
    result = run_query("""
        MATCH (exp:Experiment)-[:Changes_expression_of]->()
        RETURN count(DISTINCT exp) AS cnt
    """)
    assert result[0]["cnt"] > 0, (
        "No Experiment node is a source of Changes_expression_of; "
        "expression edges may be missing"
    )


def test_expression_edges_target_genes(run_query):
    """All expression edges must target a Gene node."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->(target)
        WHERE NOT target:Gene
        RETURN count(e) AS wrong_targets,
               collect(DISTINCT labels(target))[..5] AS target_labels
    """)
    bad = result[0]["wrong_targets"]
    assert bad == 0, (
        f"{bad} expression edges point to non-Gene nodes: "
        f"{result[0]['target_labels']}"
    )


def test_expression_edges_source_from_experiments(run_query):
    """All Changes_expression_of edges must originate from an Experiment node."""
    result = run_query("""
        MATCH (src)-[e:Changes_expression_of]->()
        WHERE NOT src:Experiment
        RETURN count(e) AS wrong_sources,
               collect(DISTINCT labels(src))[..5] AS source_labels
    """)
    bad = result[0]["wrong_sources"]
    assert bad == 0, (
        f"{bad} expression edges originate from non-Experiment nodes: "
        f"{result[0]['source_labels']}"
    )


# ---------------------------------------------------------------------------
# Experiment node properties
# ---------------------------------------------------------------------------

def test_experiment_name_present_on_majority(run_query):
    """
    Most Experiment nodes should have a name property.
    The analysis 'name' field is present in all paperconfigs with named analyses.
    """
    result = run_query("""
        MATCH (exp:Experiment)
        RETURN
          count(exp) AS total,
          count(exp.name) AS with_name
    """)
    total = result[0]["total"]
    with_name = result[0]["with_name"]
    if total == 0:
        pytest.skip("No Experiment nodes found")
    coverage = with_name / total
    assert coverage >= 0.9, (
        f"Only {with_name}/{total} ({coverage:.1%}) Experiment nodes have name; "
        f"expected >= 90%"
    )


def test_experiment_has_treatment_type(run_query):
    """Every Experiment node must have a treatment_type property."""
    result = run_query("""
        MATCH (exp:Experiment)
        WHERE exp.treatment_type IS NULL OR exp.treatment_type = ''
        RETURN count(exp) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Experiment nodes are missing treatment_type"
    )


def test_experiment_has_organism_strain(run_query):
    """Every Experiment node must have an organism_strain property."""
    result = run_query("""
        MATCH (exp:Experiment)
        WHERE exp.organism_strain IS NULL OR exp.organism_strain = ''
        RETURN count(exp) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Experiment nodes are missing organism_strain"
    )


def test_experiment_has_omics_type(run_query):
    """Every Experiment node must have an omics_type property."""
    result = run_query("""
        MATCH (exp:Experiment)
        WHERE exp.omics_type IS NULL OR exp.omics_type = ''
        RETURN count(exp) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Experiment nodes are missing omics_type"
    )


# ---------------------------------------------------------------------------
# No old node/edge types
# ---------------------------------------------------------------------------

def test_no_environmental_condition_nodes(run_query):
    """EnvironmentalCondition nodes must not exist (replaced by Experiment)."""
    result = run_query(
        "MATCH (e:EnvironmentalCondition) RETURN count(e) AS cnt"
    )
    assert result[0]["cnt"] == 0, (
        f"{result[0]['cnt']} EnvironmentalCondition nodes still exist"
    )


def test_no_old_expression_edge_types(run_query):
    """Old Condition_changes_expression_of / Coculture_changes_expression_of edges must not exist."""
    result = run_query(
        "MATCH ()-[r:Condition_changes_expression_of|Coculture_changes_expression_of]->() "
        "RETURN count(r) AS cnt"
    )
    assert result[0]["cnt"] == 0, (
        f"{result[0]['cnt']} old-style expression edges still exist"
    )
