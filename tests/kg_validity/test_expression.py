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
    """Every Experiment node must have a non-empty treatment_type list."""
    result = run_query("""
        MATCH (exp:Experiment)
        WHERE exp.treatment_type IS NULL OR size(exp.treatment_type) = 0
        RETURN count(exp) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Experiment nodes are missing treatment_type"
    )


def test_experiment_treatment_type_values_canonical(run_query):
    """All treatment_type values should be from the canonical vocabulary."""
    result = run_query("""
        MATCH (e:Experiment)
        UNWIND e.treatment_type AS tt
        WITH DISTINCT tt
        RETURN collect(tt) AS all_types
    """)
    known = {
        "nitrogen_stress", "phosphorus_stress", "iron_stress", "carbon_stress",
        "salt_stress", "oxygen_stress", "temperature_stress", "light_stress",
        "darkness", "plastic_stress", "viral", "coculture", "growth_state",
        "growth_medium", "diel", "axenic", "continuous_light", "diel_cycle",
    }
    actual = set(result[0]["all_types"])
    unknown = actual - known
    assert not unknown, f"Unknown treatment_type values in KG: {unknown}"


def test_experiment_has_organism_name(run_query):
    """Every Experiment node must have an organism_name property."""
    result = run_query("""
        MATCH (exp:Experiment)
        WHERE exp.organism_name IS NULL OR exp.organism_name = ''
        RETURN count(exp) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Experiment nodes are missing organism_name"
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


def test_experiment_table_scope_valid(run_query):
    """Every Experiment with table_scope must use one of the allowed values."""
    valid = {
        "all_detected_genes", "significant_any_timepoint",
        "significant_only", "top_n", "filtered_subset",
    }
    result = run_query("""
        MATCH (exp:Experiment)
        WHERE exp.table_scope IS NOT NULL AND exp.table_scope <> ''
        RETURN exp.table_scope AS scope, count(exp) AS cnt
    """)
    bad = [(r["scope"], r["cnt"]) for r in result if r["scope"] not in valid]
    assert not bad, f"Invalid table_scope values: {bad}"


def test_experiment_table_scope_coverage(run_query):
    """Most Experiment nodes should have table_scope set."""
    result = run_query("""
        MATCH (exp:Experiment)
        RETURN
          count(exp) AS total,
          count(CASE WHEN exp.table_scope IS NOT NULL AND exp.table_scope <> ''
                     THEN 1 END) AS with_scope
    """)
    total = result[0]["total"]
    with_scope = result[0]["with_scope"]
    if total == 0:
        pytest.skip("No Experiment nodes found")
    coverage = with_scope / total
    assert coverage >= 0.9, (
        f"Only {with_scope}/{total} ({coverage:.1%}) Experiment nodes have "
        f"table_scope; expected >= 90%"
    )


# ---------------------------------------------------------------------------
# No old node/edge types
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# expression_status (derived property)
# ---------------------------------------------------------------------------

def test_expression_status_populated(run_query):
    """Every expression edge must have expression_status set."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE e.expression_status IS NULL
        RETURN count(e) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} expression edges are missing expression_status"
    )


def test_expression_status_valid_values(run_query):
    """expression_status must be one of the three allowed values."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE NOT e.expression_status IN ['significant_up', 'significant_down', 'not_significant']
        RETURN count(e) AS bad,
               collect(DISTINCT e.expression_status)[..10] AS bad_values
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} edges have invalid expression_status: {result[0]['bad_values']}"
    )


def test_expression_status_consistent_with_significant(run_query):
    """expression_status must be consistent with significant + expression_direction."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE (e.expression_status = 'significant_up'
               AND NOT (e.significant = 'significant' AND e.expression_direction = 'up'))
           OR (e.expression_status = 'significant_down'
               AND NOT (e.significant = 'significant' AND e.expression_direction = 'down'))
           OR (e.expression_status = 'not_significant'
               AND e.significant = 'significant')
        RETURN count(e) AS inconsistent
    """)
    assert result[0]["inconsistent"] == 0, (
        f"{result[0]['inconsistent']} edges have expression_status inconsistent with significant/direction"
    )


def test_expression_status_counts_match_significant(run_query):
    """Total significant_up + significant_down should equal total where significant = 'significant'."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        RETURN
          count(CASE WHEN e.expression_status IN ['significant_up', 'significant_down'] THEN 1 END) AS status_sig,
          count(CASE WHEN e.significant = 'significant' THEN 1 END) AS old_sig
    """)
    assert result[0]["status_sig"] == result[0]["old_sig"], (
        f"expression_status significant count ({result[0]['status_sig']}) != "
        f"significant property count ({result[0]['old_sig']})"
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


# ---------------------------------------------------------------------------
# Directional ranks (rank_up, rank_down — post-import computed)
# ---------------------------------------------------------------------------

def test_rank_up_only_on_significant_up(run_query):
    """rank_up must be null on edges that are not significant_up."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE e.rank_up IS NOT NULL
          AND e.expression_status <> 'significant_up'
        RETURN count(e) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} edges have rank_up set but are not significant_up"
    )


def test_rank_down_only_on_significant_down(run_query):
    """rank_down must be null on edges that are not significant_down."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE e.rank_down IS NOT NULL
          AND e.expression_status <> 'significant_down'
        RETURN count(e) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} edges have rank_down set but are not significant_down"
    )


def test_rank_up_covers_all_significant_up(run_query):
    """Every significant_up edge must have rank_up set."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE e.expression_status = 'significant_up'
          AND e.rank_up IS NULL
        RETURN count(e) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} significant_up edges are missing rank_up"
    )


def test_rank_down_covers_all_significant_down(run_query):
    """Every significant_down edge must have rank_down set."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE e.expression_status = 'significant_down'
          AND e.rank_down IS NULL
        RETURN count(e) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} significant_down edges are missing rank_down"
    )


def test_rank_up_starts_at_one(run_query):
    """Each experiment x timepoint must have a rank_up = 1 if any significant_up edges exist."""
    result = run_query("""
        MATCH (exp:Experiment)-[r:Changes_expression_of]->()
        WHERE r.expression_status = 'significant_up'
        WITH exp.id AS eid, r.time_point_order AS tp,
             min(r.rank_up) AS min_rank
        WHERE min_rank <> 1
        RETURN count(*) AS bad,
               collect(eid + ':tp' + toString(tp))[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} experiment x timepoint groups have rank_up not starting at 1: "
        f"{result[0]['examples']}"
    )


def test_rank_down_starts_at_one(run_query):
    """Each experiment x timepoint must have a rank_down = 1 if any significant_down edges exist."""
    result = run_query("""
        MATCH (exp:Experiment)-[r:Changes_expression_of]->()
        WHERE r.expression_status = 'significant_down'
        WITH exp.id AS eid, r.time_point_order AS tp,
             min(r.rank_down) AS min_rank
        WHERE min_rank <> 1
        RETURN count(*) AS bad,
               collect(eid + ':tp' + toString(tp))[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} experiment x timepoint groups have rank_down not starting at 1: "
        f"{result[0]['examples']}"
    )


def test_rank_up_contiguous(run_query):
    """rank_up values must be contiguous 1..N within each experiment x timepoint."""
    result = run_query("""
        MATCH (exp:Experiment)-[r:Changes_expression_of]->()
        WHERE r.expression_status = 'significant_up'
        WITH exp.id AS eid, r.time_point_order AS tp,
             count(r) AS n, max(r.rank_up) AS max_rank
        WHERE max_rank <> n
        RETURN count(*) AS bad,
               collect(eid + ':tp' + toString(tp) + ' n=' + toString(n) + ' max=' + toString(max_rank))[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} experiment x timepoint groups have non-contiguous rank_up: "
        f"{result[0]['examples']}"
    )


def test_rank_down_contiguous(run_query):
    """rank_down values must be contiguous 1..N within each experiment x timepoint."""
    result = run_query("""
        MATCH (exp:Experiment)-[r:Changes_expression_of]->()
        WHERE r.expression_status = 'significant_down'
        WITH exp.id AS eid, r.time_point_order AS tp,
             count(r) AS n, max(r.rank_down) AS max_rank
        WHERE max_rank <> n
        RETURN count(*) AS bad,
               collect(eid + ':tp' + toString(tp) + ' n=' + toString(n) + ' max=' + toString(max_rank))[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} experiment x timepoint groups have non-contiguous rank_down: "
        f"{result[0]['examples']}"
    )


def test_rank_up_rank_down_mutually_exclusive(run_query):
    """No edge should have both rank_up and rank_down set."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE e.rank_up IS NOT NULL AND e.rank_down IS NOT NULL
        RETURN count(e) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} edges have both rank_up and rank_down set"
    )


def test_not_significant_edges_have_no_directional_rank(run_query):
    """not_significant edges must have both rank_up and rank_down as null."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE e.expression_status = 'not_significant'
          AND (e.rank_up IS NOT NULL OR e.rank_down IS NOT NULL)
        RETURN count(e) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} not_significant edges have a directional rank set"
    )


# ---------------------------------------------------------------------------
# rank_by_effect (post-import computed — all edges regardless of direction)
# ---------------------------------------------------------------------------

def test_rank_by_effect_populated(run_query):
    """Every expression edge must have rank_by_effect set."""
    result = run_query("""
        MATCH ()-[e:Changes_expression_of]->()
        WHERE e.rank_by_effect IS NULL
        RETURN count(e) AS missing
    """)
    assert result[0]["missing"] == 0, (
        f"{result[0]['missing']} expression edges are missing rank_by_effect"
    )


def test_rank_by_effect_starts_at_one(run_query):
    """Each experiment x timepoint must have rank_by_effect = 1."""
    result = run_query("""
        MATCH (exp:Experiment)-[r:Changes_expression_of]->()
        WITH exp.id AS eid, r.time_point_order AS tp,
             min(r.rank_by_effect) AS min_rank
        WHERE min_rank <> 1
        RETURN count(*) AS bad,
               collect(eid + ':tp' + toString(tp))[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} experiment x timepoint groups have rank_by_effect not starting at 1: "
        f"{result[0]['examples']}"
    )


def test_rank_by_effect_contiguous(run_query):
    """rank_by_effect must be contiguous 1..N within each experiment x timepoint."""
    result = run_query("""
        MATCH (exp:Experiment)-[r:Changes_expression_of]->()
        WITH exp.id AS eid, r.time_point_order AS tp,
             count(r) AS n, max(r.rank_by_effect) AS max_rank
        WHERE max_rank <> n
        RETURN count(*) AS bad,
               collect(eid + ':tp' + toString(tp) + ' n=' + toString(n) + ' max=' + toString(max_rank))[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} experiment x timepoint groups have non-contiguous rank_by_effect: "
        f"{result[0]['examples']}"
    )
