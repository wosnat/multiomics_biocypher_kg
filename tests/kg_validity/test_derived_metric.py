"""
KG validity tests for DerivedMetric nodes and the 6 new edge types
(3 binding + 3 measurement). Parent-spec §Success-criteria assertions
filtered to DerivedMetric scope.

Covers:
- Node presence and minimum counts
- Binding edge cardinality (1:many from parent to DM)
- Denormalized fields match parent Experiment
- Edge-target constraints (all measurement edges target :Gene)
- Per-edge-type value constraints (value_flag enum, value_text in allowed_categories)
- Rollup consistency (Experiment/Publication/OrganismTaxon/Gene) with direct query-time counts
- Empty-state defaults on nodes without DM children
- Index existence
"""

import pytest

pytestmark = pytest.mark.kg


# ---------------------------------------------------------------------------
# Node presence
# ---------------------------------------------------------------------------

def test_derived_metric_nodes_exist(run_query):
    """DerivedMetric nodes must exist.

    Biller 2018 contributes 7 (6 boolean + 1 categorical); the synthetic
    fixture (when wired) adds 3 more. We assert the production floor (>=7)
    so this test passes both in fixture-enabled and fixture-free graphs.
    """
    result = run_query("MATCH (dm:DerivedMetric) RETURN count(dm) AS cnt")
    assert result[0]["cnt"] >= 7, (
        f"Expected at least 7 DerivedMetric nodes (Biller 2018 floor), got {result[0]['cnt']}"
    )


@pytest.mark.parametrize("prop", [
    "name",
    "metric_type",
    "value_kind",
    "experiment_id",
    "organism_name",
    "publication_doi",
    "compartment",
    "omics_type",
    "rankable",
    "has_p_value",
    "total_gene_count",
])
def test_derived_metric_required_properties(run_query, prop):
    """Every DerivedMetric must have these non-null properties."""
    result = run_query(
        f"MATCH (dm:DerivedMetric) WHERE dm.{prop} IS NULL RETURN count(dm) AS cnt"
    )
    assert result[0]["cnt"] == 0, f"Found DerivedMetric nodes with null {prop}"


def test_derived_metric_value_kind_enum(run_query):
    """value_kind must be one of numeric/boolean/categorical."""
    result = run_query("""
        MATCH (dm:DerivedMetric)
        WHERE NOT dm.value_kind IN ['numeric', 'boolean', 'categorical']
        RETURN collect(DISTINCT dm.value_kind) AS bad
    """)
    assert result[0]["bad"] == [], f"Unexpected value_kind values: {result[0]['bad']}"


def test_derived_metric_rankable_enum(run_query):
    """rankable must be string 'true' or 'false' (never bool, never null)."""
    result = run_query("""
        MATCH (dm:DerivedMetric)
        WHERE NOT dm.rankable IN ['true', 'false']
        RETURN count(dm) AS bad
    """)
    assert result[0]["bad"] == 0


def test_derived_metric_has_p_value_enum(run_query):
    """has_p_value must be string 'true' or 'false'."""
    result = run_query("""
        MATCH (dm:DerivedMetric)
        WHERE NOT dm.has_p_value IN ['true', 'false']
        RETURN count(dm) AS bad
    """)
    assert result[0]["bad"] == 0


def test_derived_metric_rankable_only_numeric(run_query):
    """rankable='true' must imply value_kind='numeric'."""
    result = run_query("""
        MATCH (dm:DerivedMetric)
        WHERE dm.rankable = 'true' AND dm.value_kind <> 'numeric'
        RETURN count(dm) AS bad, collect(dm.id)[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} DMs have rankable='true' but value_kind != 'numeric': {result[0]['examples']}"
    )


# ---------------------------------------------------------------------------
# Binding edge cardinality (1:many parent -> DM; exactly 1 parent per DM)
# ---------------------------------------------------------------------------

def test_every_dm_has_exactly_one_publication_parent(run_query):
    result = run_query("""
        MATCH (dm:DerivedMetric)
        OPTIONAL MATCH (p:Publication)-[:PublicationHasDerivedMetric]->(dm)
        WITH dm, count(p) AS pc
        WHERE pc <> 1
        RETURN count(dm) AS bad, collect(dm.id)[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} DMs have != 1 Publication parent: {result[0]['examples']}"
    )


def test_every_dm_has_exactly_one_experiment_parent(run_query):
    result = run_query("""
        MATCH (dm:DerivedMetric)
        OPTIONAL MATCH (e:Experiment)-[:ExperimentHasDerivedMetric]->(dm)
        WITH dm, count(e) AS ec
        WHERE ec <> 1
        RETURN count(dm) AS bad, collect(dm.id)[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} DMs have != 1 Experiment parent: {result[0]['examples']}"
    )


def test_every_dm_has_exactly_one_organism(run_query):
    result = run_query("""
        MATCH (dm:DerivedMetric)
        OPTIONAL MATCH (dm)-[:DerivedMetricBelongsToOrganism]->(o:OrganismTaxon)
        WITH dm, count(o) AS oc
        WHERE oc <> 1
        RETURN count(dm) AS bad, collect(dm.id)[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} DMs have != 1 Organism parent: {result[0]['examples']}"
    )


# ---------------------------------------------------------------------------
# Denormalized field consistency
# ---------------------------------------------------------------------------

def test_dm_denormalized_experiment_id_matches_parent(run_query):
    result = run_query("""
        MATCH (e:Experiment)-[:ExperimentHasDerivedMetric]->(dm:DerivedMetric)
        WHERE dm.experiment_id <> e.id
        RETURN count(dm) AS bad, collect(dm.id)[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} DMs have experiment_id mismatch: {result[0]['examples']}"
    )


@pytest.mark.parametrize("field", [
    "organism_name", "compartment", "omics_type", "treatment",
    "light_condition", "experimental_context",
])
def test_dm_denormalized_scalar_matches_parent(run_query, field):
    result = run_query(f"""
        MATCH (e:Experiment)-[:ExperimentHasDerivedMetric]->(dm:DerivedMetric)
        WHERE coalesce(dm.{field}, '') <> coalesce(e.{field}, '')
        RETURN count(dm) AS bad, collect(dm.id)[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} DMs have {field} != parent Experiment: {result[0]['examples']}"
    )


def test_dm_publication_doi_matches_parent(run_query):
    result = run_query("""
        MATCH (p:Publication)-[:PublicationHasDerivedMetric]->(dm:DerivedMetric)
        WHERE dm.publication_doi <> p.doi
        RETURN count(dm) AS bad
    """)
    assert result[0]["bad"] == 0


# ---------------------------------------------------------------------------
# Edge target constraints: all measurement edges must target :Gene
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("rel_type", [
    "Derived_metric_quantifies_gene",
    "Derived_metric_flags_gene",
    "Derived_metric_classifies_gene",
])
def test_measurement_edges_target_gene(run_query, rel_type):
    result = run_query(f"""
        MATCH (dm:DerivedMetric)-[r:`{rel_type}`]->(t)
        WHERE NOT t:Gene
        RETURN count(*) AS bad
    """)
    assert result[0]["bad"] == 0


# ---------------------------------------------------------------------------
# Per-edge-type value constraints
# ---------------------------------------------------------------------------

def test_flag_edges_value_flag_enum(run_query):
    """Every derived_metric_flags_gene edge has value_flag in {'true','false'}."""
    result = run_query("""
        MATCH ()-[r:Derived_metric_flags_gene]->()
        WHERE r.value_flag IS NULL OR NOT r.value_flag IN ['true', 'false']
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0


def test_classify_edges_value_text_non_null(run_query):
    result = run_query("""
        MATCH ()-[r:Derived_metric_classifies_gene]->()
        WHERE r.value_text IS NULL OR r.value_text = ''
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0


def test_classify_edges_value_text_in_allowed_categories(run_query):
    """value_text on every classifies edge must be in parent DM's allowed_categories."""
    result = run_query("""
        MATCH (dm:DerivedMetric)-[r:Derived_metric_classifies_gene]->()
        WHERE NOT r.value_text IN dm.allowed_categories
        RETURN count(r) AS bad, collect(r.value_text)[..5] AS bad_values
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} classifies edges have value_text outside allowed_categories: {result[0]['bad_values']}"
    )


def test_quantify_edges_value_non_null(run_query):
    """Every quantifies edge must have a non-null numeric value."""
    result = run_query("""
        MATCH ()-[r:Derived_metric_quantifies_gene]->()
        WHERE r.value IS NULL
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0


def test_quantify_edges_metric_type_matches_parent(run_query):
    """Edge.metric_type must equal parent DM.metric_type on quantifies edges."""
    result = run_query("""
        MATCH (dm:DerivedMetric)-[r:Derived_metric_quantifies_gene]->()
        WHERE r.metric_type <> dm.metric_type
        RETURN count(r) AS bad
    """)
    assert result[0]["bad"] == 0


# ---------------------------------------------------------------------------
# DerivedMetric total_gene_count consistency
# ---------------------------------------------------------------------------

def test_dm_total_gene_count_matches_edge_count(run_query):
    """total_gene_count equals actual count of outgoing measurement edges."""
    result = run_query("""
        MATCH (dm:DerivedMetric)
        OPTIONAL MATCH (dm)-[r:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->()
        WITH dm, dm.total_gene_count AS declared, count(r) AS actual
        WHERE declared <> actual
        RETURN count(dm) AS mismatched, collect(dm.id)[..3] AS examples
    """)
    assert result[0]["mismatched"] == 0, (
        f"{result[0]['mismatched']} DMs have total_gene_count != edge count: {result[0]['examples']}"
    )


def test_dm_emits_only_one_edge_type(run_query):
    """Each DerivedMetric emits exactly one of the 3 measurement edge types, matching its value_kind."""
    result = run_query("""
        MATCH (dm:DerivedMetric)
        WITH dm,
             count { (dm)-[:Derived_metric_quantifies_gene]->() } AS n_quant,
             count { (dm)-[:Derived_metric_flags_gene]->() } AS n_flag,
             count { (dm)-[:Derived_metric_classifies_gene]->() } AS n_class
        WITH dm, dm.value_kind AS vk, n_quant, n_flag, n_class
        WHERE (vk = 'numeric' AND (n_flag > 0 OR n_class > 0))
           OR (vk = 'boolean' AND (n_quant > 0 OR n_class > 0))
           OR (vk = 'categorical' AND (n_quant > 0 OR n_flag > 0))
        RETURN count(dm) AS bad, collect(dm.id)[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} DMs emit an edge type that doesn't match their value_kind: {result[0]['examples']}"
    )


# ---------------------------------------------------------------------------
# Experiment DM rollup consistency
# ---------------------------------------------------------------------------

def test_experiment_derived_metric_count_matches_edges(run_query):
    result = run_query("""
        MATCH (e:Experiment)
        OPTIONAL MATCH (e)-[:ExperimentHasDerivedMetric]->(dm:DerivedMetric)
        WITH e, e.derived_metric_count AS declared, count(DISTINCT dm) AS actual
        WHERE declared <> actual
        RETURN count(e) AS mismatched, collect(e.id)[..3] AS examples
    """)
    assert result[0]["mismatched"] == 0, f"{result[0]['mismatched']} Experiments: {result[0]['examples']}"


def test_experiment_derived_metric_gene_count_matches_query(run_query):
    result = run_query("""
        MATCH (e:Experiment)
        OPTIONAL MATCH (e)-[:ExperimentHasDerivedMetric]->(:DerivedMetric)
          -[:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(g:Gene)
        WITH e, e.derived_metric_gene_count AS declared, count(DISTINCT g) AS actual
        WHERE declared <> actual
        RETURN count(e) AS mismatched, collect(e.id)[..3] AS examples
    """)
    assert result[0]["mismatched"] == 0


def test_experiment_reports_fold_change_enum(run_query):
    result = run_query("""
        MATCH (e:Experiment)
        WHERE NOT e.reports_fold_change IN ['true', 'false']
        RETURN count(e) AS bad
    """)
    assert result[0]["bad"] == 0


def test_experiment_reports_fold_change_consistent(run_query):
    """reports_fold_change='true' iff Experiment has any Changes_expression_of edge."""
    result = run_query("""
        MATCH (e:Experiment)
        WITH e, e.reports_fold_change AS declared,
             EXISTS { (e)-[:Changes_expression_of]->() } AS has_de
        WHERE (declared = 'true') <> has_de
        RETURN count(e) AS mismatched, collect(e.id)[..3] AS examples
    """)
    assert result[0]["mismatched"] == 0


# ---------------------------------------------------------------------------
# Publication DM rollup consistency
# ---------------------------------------------------------------------------

def test_publication_derived_metric_count_matches_edges(run_query):
    result = run_query("""
        MATCH (p:Publication)
        OPTIONAL MATCH (p)-[:PublicationHasDerivedMetric]->(dm:DerivedMetric)
        WITH p, p.derived_metric_count AS declared, count(DISTINCT dm) AS actual
        WHERE declared <> actual
        RETURN count(p) AS mismatched, collect(p.doi)[..3] AS examples
    """)
    assert result[0]["mismatched"] == 0, f"{result[0]['mismatched']}: {result[0]['examples']}"


def test_publication_derived_metric_gene_count_matches_query(run_query):
    result = run_query("""
        MATCH (p:Publication)
        OPTIONAL MATCH (p)-[:PublicationHasDerivedMetric]->(:DerivedMetric)
          -[:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(g:Gene)
        WITH p, p.derived_metric_gene_count AS declared, count(DISTINCT g) AS actual
        WHERE declared <> actual
        RETURN count(p) AS mismatched, collect(p.doi)[..3] AS examples
    """)
    assert result[0]["mismatched"] == 0


def test_publication_compartments_matches_child_experiments(run_query):
    result = run_query("""
        MATCH (p:Publication)
        OPTIONAL MATCH (p)-[:Has_experiment]->(e:Experiment)
        WITH p, p.compartments AS declared,
             apoc.coll.sort([x IN collect(DISTINCT e.compartment) WHERE x IS NOT NULL]) AS actual
        WHERE declared <> actual
        RETURN count(p) AS mismatched, collect(p.doi)[..3] AS examples
    """)
    assert result[0]["mismatched"] == 0


# ---------------------------------------------------------------------------
# OrganismTaxon DM rollup consistency
# ---------------------------------------------------------------------------

def test_organism_derived_metric_count_matches_edges(run_query):
    result = run_query("""
        MATCH (o:OrganismTaxon)
        OPTIONAL MATCH (dm:DerivedMetric)-[:DerivedMetricBelongsToOrganism]->(o)
        WITH o, o.derived_metric_count AS declared, count(DISTINCT dm) AS actual
        WHERE declared <> actual
        RETURN count(o) AS mismatched, collect(o.preferred_name)[..3] AS examples
    """)
    assert result[0]["mismatched"] == 0, f"{result[0]['mismatched']}: {result[0]['examples']}"


# ---------------------------------------------------------------------------
# Gene routing count consistency
# ---------------------------------------------------------------------------

def test_gene_classifier_flag_count_matches_edges(run_query):
    """On a sample of genes with flag edges, rollup equals distinct-DM count."""
    result = run_query("""
        MATCH (g:Gene)
        WHERE g.classifier_flag_count > 0
        WITH g LIMIT 20
        OPTIONAL MATCH (dm:DerivedMetric)-[:Derived_metric_flags_gene]->(g)
        WITH g, g.classifier_flag_count AS declared, count(DISTINCT dm) AS actual
        WHERE declared <> actual
        RETURN count(g) AS mismatched, collect(g.locus_tag)[..5] AS examples
    """)
    assert result[0]["mismatched"] == 0


def test_gene_classifier_label_count_matches_edges(run_query):
    result = run_query("""
        MATCH (g:Gene)
        WHERE g.classifier_label_count > 0
        WITH g LIMIT 20
        OPTIONAL MATCH (dm:DerivedMetric)-[:Derived_metric_classifies_gene]->(g)
        WITH g, g.classifier_label_count AS declared, count(DISTINCT dm) AS actual
        WHERE declared <> actual
        RETURN count(g) AS mismatched, collect(g.locus_tag)[..5] AS examples
    """)
    assert result[0]["mismatched"] == 0


def test_gene_compartments_observed_consistent(run_query):
    """compartments_observed equals union over all 3 DM edge types."""
    result = run_query("""
        MATCH (g:Gene)
        WHERE size(g.compartments_observed) > 0
        WITH g LIMIT 20
        OPTIONAL MATCH (dm:DerivedMetric)
          -[:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(g)
        WITH g, g.compartments_observed AS declared,
             apoc.coll.sort([x IN collect(DISTINCT dm.compartment) WHERE x IS NOT NULL]) AS actual
        WHERE declared <> actual
        RETURN count(g) AS mismatched, collect(g.locus_tag)[..5] AS examples
    """)
    assert result[0]["mismatched"] == 0


# ---------------------------------------------------------------------------
# Empty-state defaults on nodes without DM children
# ---------------------------------------------------------------------------

def test_experiment_dm_empty_state_defaults(run_query):
    """Experiments without DM children: all DM rollup fields are 0 / [], never null."""
    result = run_query("""
        MATCH (e:Experiment)
        WHERE NOT EXISTS { (e)-[:ExperimentHasDerivedMetric]->() }
          AND (e.derived_metric_count IS NULL
               OR e.derived_metric_count <> 0
               OR e.derived_metric_value_kinds IS NULL
               OR size(e.derived_metric_value_kinds) <> 0
               OR e.reports_derived_metric_types IS NULL
               OR size(e.reports_derived_metric_types) <> 0
               OR e.derived_metric_gene_count IS NULL
               OR e.derived_metric_gene_count <> 0)
        RETURN count(e) AS bad, collect(e.id)[..3] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} Experiments without DM children have wrong defaults: {result[0]['examples']}"
    )


def test_gene_dm_empty_state_defaults(run_query):
    """All genes have non-null DM routing defaults (never null)."""
    result = run_query("""
        MATCH (g:Gene)
        WHERE g.numeric_metric_count IS NULL
           OR g.classifier_flag_count IS NULL
           OR g.classifier_label_count IS NULL
           OR g.numeric_metric_types_observed IS NULL
           OR g.classifier_flag_types_observed IS NULL
           OR g.classifier_label_types_observed IS NULL
           OR g.compartments_observed IS NULL
        RETURN count(g) AS bad
    """)
    assert result[0]["bad"] == 0


# ---------------------------------------------------------------------------
# Index existence (extends test_post_import.EXPECTED_INDEXES)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("idx_name", [
    "derived_metric_metric_type_idx",
    "derived_metric_value_kind_idx",
    "derived_metric_compartment_idx",
    "derived_metric_omics_type_idx",
    "derived_metric_treatment_type_idx",
    "derived_metric_organism_idx",
    "derived_metric_experiment_idx",
    "experiment_compartment_idx",
    "derivedMetricFullText",
])
def test_plan3_index_exists(run_query, idx_name):
    result = run_query(
        f"SHOW INDEXES YIELD name, state WHERE name = '{idx_name}' "
        f"RETURN state"
    )
    assert len(result) == 1, f"Index {idx_name} not found"
    assert result[0]["state"] == "ONLINE", f"Index {idx_name} state = {result[0]['state']}"
