#!/bin/bash
set -euo pipefail

TIMEFORMAT="  [timing] %Rs"

echo "=== Post-process: Starting Neo4j ==="
neo4j start

# Wait until Neo4j accepts Bolt connections. A fixed sleep is unreliable on
# slower machines — Neo4j may need 30-60s to open the Bolt port on a cold start.
echo "=== Post-process: Waiting for Neo4j to accept connections ==="
for i in $(seq 1 60); do
  if cypher-shell -a bolt://localhost:7687 "RETURN 1;" >/dev/null 2>&1; then
    echo "  Neo4j ready after ${i}s"
    break
  fi
  if [ "$i" = "60" ]; then
    echo "ERROR: Neo4j did not become ready within 60s" >&2
    cat /logs/neo4j.log 2>/dev/null || true
    exit 1
  fi
  sleep 1
done

# ─────────────────────────────────────────────────────────────────────────────
# Group 1: all indexes (scalar + full-text across every node type).
# ─────────────────────────────────────────────────────────────────────────────
echo "=== Post-process: Create all indexes ==="
time cypher-shell <<'CYPHER'
// Gene scalar + full-text
CREATE INDEX gene_locus_tag_idx IF NOT EXISTS FOR (g:Gene) ON (g.locus_tag);
CREATE INDEX gene_name_idx IF NOT EXISTS FOR (g:Gene) ON (g.gene_name);
CREATE INDEX gene_organism_name_idx IF NOT EXISTS FOR (g:Gene) ON (g.organism_name);
CREATE FULLTEXT INDEX geneFullText IF NOT EXISTS FOR (n:Gene) ON EACH [
  n.gene_summary, n.all_identifiers, n.gene_name_synonyms,
  n.alternate_functional_descriptions];

// Ontology / role full-text
CREATE FULLTEXT INDEX biologicalProcessFullText IF NOT EXISTS
  FOR (n:BiologicalProcess) ON EACH [n.name];
CREATE FULLTEXT INDEX molecularFunctionFullText IF NOT EXISTS
  FOR (n:MolecularFunction) ON EACH [n.name];
CREATE FULLTEXT INDEX cellularComponentFullText IF NOT EXISTS
  FOR (n:CellularComponent) ON EACH [n.name];
CREATE FULLTEXT INDEX ecNumberFullText IF NOT EXISTS
  FOR (n:EcNumber) ON EACH [n.name];
CREATE FULLTEXT INDEX keggFullText IF NOT EXISTS
  FOR (n:KeggTerm) ON EACH [n.name];
CREATE FULLTEXT INDEX cogCategoryFullText IF NOT EXISTS
  FOR (n:CogFunctionalCategory) ON EACH [n.name];
CREATE FULLTEXT INDEX cyanorakRoleFullText IF NOT EXISTS
  FOR (n:CyanorakRole) ON EACH [n.name];
CREATE FULLTEXT INDEX tigrRoleFullText IF NOT EXISTS
  FOR (n:TigrRole) ON EACH [n.name];

// OrthologGroup
CREATE INDEX ortholog_group_id_idx IF NOT EXISTS FOR (og:OrthologGroup) ON (og.id);
CREATE INDEX ortholog_group_name_idx IF NOT EXISTS FOR (og:OrthologGroup) ON (og.name);
CREATE INDEX ortholog_group_level_idx IF NOT EXISTS FOR (og:OrthologGroup) ON (og.taxonomic_level);
CREATE INDEX ortholog_group_rank_idx IF NOT EXISTS FOR (og:OrthologGroup) ON (og.specificity_rank);
CREATE FULLTEXT INDEX orthologGroupFullText IF NOT EXISTS
  FOR (og:OrthologGroup) ON EACH [og.consensus_product, og.consensus_gene_name, og.description, og.functional_description];

// Pfam
CREATE INDEX pfam_name_idx IF NOT EXISTS FOR (p:Pfam) ON (p.name);
CREATE INDEX pfam_clan_name_idx IF NOT EXISTS FOR (c:PfamClan) ON (c.name);
CREATE FULLTEXT INDEX pfamFullText IF NOT EXISTS
  FOR (p:Pfam) ON EACH [p.name, p.short_name];
CREATE FULLTEXT INDEX pfamClanFullText IF NOT EXISTS
  FOR (c:PfamClan) ON EACH [c.name];

// BriteCategory
CREATE INDEX brite_category_tree_idx IF NOT EXISTS FOR (b:BriteCategory) ON (b.tree_code);
CREATE INDEX brite_category_level_idx IF NOT EXISTS FOR (b:BriteCategory) ON (b.level);
CREATE INDEX brite_category_name_idx IF NOT EXISTS FOR (b:BriteCategory) ON (b.name);
CREATE FULLTEXT INDEX briteCategoryFullText IF NOT EXISTS
  FOR (b:BriteCategory) ON EACH [b.name];

// TCDB / CAZy scalar indexes
CREATE INDEX tcdb_family_level_idx IF NOT EXISTS FOR (t:TcdbFamily) ON (t.level);
CREATE INDEX tcdb_family_level_kind_idx IF NOT EXISTS FOR (t:TcdbFamily) ON (t.level_kind);
CREATE INDEX tcdb_family_tcdb_id_idx IF NOT EXISTS FOR (t:TcdbFamily) ON (t.tcdb_id);
CREATE INDEX tcdb_family_tc_class_id_idx IF NOT EXISTS FOR (t:TcdbFamily) ON (t.tc_class_id);
CREATE INDEX cazy_family_level_idx IF NOT EXISTS FOR (c:CazyFamily) ON (c.level);
CREATE INDEX cazy_family_level_kind_idx IF NOT EXISTS FOR (c:CazyFamily) ON (c.level_kind);
CREATE INDEX cazy_family_cazy_id_idx IF NOT EXISTS FOR (c:CazyFamily) ON (c.cazy_id);

// TCDB / CAZy full-text indexes
CREATE FULLTEXT INDEX tcdbFamilyFullText IF NOT EXISTS
    FOR (t:TcdbFamily) ON EACH [t.name, t.tcdb_id, t.superfamily];
CREATE FULLTEXT INDEX cazyFamilyFullText IF NOT EXISTS
    FOR (c:CazyFamily) ON EACH [c.name, c.cazy_id];

// Publication
// publicationFullText: drop+recreate so DM-search-text + compartments are picked up
// even on reruns against an existing graph (Neo4j won't add properties to an existing index).
DROP INDEX publicationFullText IF EXISTS;
CREATE FULLTEXT INDEX publicationFullText
  FOR (p:Publication) ON EACH [p.title, p.abstract, p.description, p.compartments, p.derived_metric_search_text]
  OPTIONS {
    indexConfig: {
      `fulltext.analyzer`: 'standard-no-stop-words',
      `fulltext.eventually_consistent`: false
    }
  };

// Experiment
CREATE INDEX experiment_id_idx IF NOT EXISTS FOR (e:Experiment) ON (e.id);
CREATE INDEX experiment_organism_idx IF NOT EXISTS FOR (e:Experiment) ON (e.organism_name);
CREATE INDEX experiment_treatment_type_idx IF NOT EXISTS FOR (e:Experiment) ON (e.treatment_type);
CREATE INDEX experiment_background_factors_idx IF NOT EXISTS FOR (e:Experiment) ON (e.background_factors);
CREATE INDEX experiment_omics_type_idx IF NOT EXISTS FOR (e:Experiment) ON (e.omics_type);
// experimentFullText: same drop+recreate as publicationFullText.
DROP INDEX experimentFullText IF EXISTS;
CREATE FULLTEXT INDEX experimentFullText
  FOR (e:Experiment) ON EACH [e.name, e.treatment, e.control, e.experimental_context, e.light_condition, e.compartment, e.derived_metric_search_text]
  OPTIONS {
    indexConfig: {
      `fulltext.analyzer`: 'standard-no-stop-words',
      `fulltext.eventually_consistent`: false
    }
  };

// OrganismTaxon
CREATE INDEX organism_type_idx IF NOT EXISTS FOR (o:OrganismTaxon) ON (o.organism_type);

// ClusteringAnalysis
CREATE INDEX clustering_analysis_organism_idx IF NOT EXISTS FOR (ca:ClusteringAnalysis) ON (ca.organism_name);
CREATE INDEX clustering_analysis_method_idx IF NOT EXISTS FOR (ca:ClusteringAnalysis) ON (ca.cluster_method);
CREATE INDEX clustering_analysis_type_idx IF NOT EXISTS FOR (ca:ClusteringAnalysis) ON (ca.cluster_type);
CREATE FULLTEXT INDEX clusteringAnalysisFullText IF NOT EXISTS
  FOR (ca:ClusteringAnalysis) ON EACH [ca.name, ca.treatment, ca.experimental_context];

// DerivedMetric scalar + full-text indexes
CREATE INDEX derived_metric_metric_type_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.metric_type);
CREATE INDEX derived_metric_value_kind_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.value_kind);
CREATE INDEX derived_metric_compartment_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.compartment);
CREATE INDEX derived_metric_omics_type_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.omics_type);
CREATE INDEX derived_metric_treatment_type_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.treatment_type);
CREATE INDEX derived_metric_organism_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.organism_name);
CREATE INDEX derived_metric_experiment_idx IF NOT EXISTS FOR (dm:DerivedMetric) ON (dm.experiment_id);
CREATE FULLTEXT INDEX derivedMetricFullText IF NOT EXISTS
  FOR (dm:DerivedMetric) ON EACH [dm.name, dm.field_description];

// Experiment.compartment scalar index (adapter-emitted by Plan 2 Task 1)
CREATE INDEX experiment_compartment_idx IF NOT EXISTS FOR (e:Experiment) ON (e.compartment);

// GeneCluster
CREATE FULLTEXT INDEX geneClusterFullText IF NOT EXISTS
  FOR (gc:GeneCluster) ON EACH [gc.name, gc.functional_description, gc.temporal_pattern, gc.expression_dynamics];

// Reaction
CREATE INDEX reaction_id_idx IF NOT EXISTS FOR (r:Reaction) ON (r.id);
CREATE INDEX reaction_kegg_id_idx IF NOT EXISTS FOR (r:Reaction) ON (r.kegg_reaction_id);
CREATE INDEX reaction_mnxr_idx IF NOT EXISTS FOR (r:Reaction) ON (r.mnxr_id);
CREATE FULLTEXT INDEX reactionFullText IF NOT EXISTS FOR (r:Reaction) ON EACH [r.name];

// Metabolite
CREATE INDEX metabolite_id_idx IF NOT EXISTS FOR (m:Metabolite) ON (m.id);
CREATE INDEX metabolite_kegg_id_idx IF NOT EXISTS FOR (m:Metabolite) ON (m.kegg_compound_id);
CREATE INDEX metabolite_mnxm_idx IF NOT EXISTS FOR (m:Metabolite) ON (m.mnxm_id);
CREATE INDEX metabolite_chebi_idx IF NOT EXISTS FOR (m:Metabolite) ON (m.chebi_id);
CREATE FULLTEXT INDEX metaboliteFullText IF NOT EXISTS FOR (m:Metabolite) ON EACH [m.name];
CYPHER

# ─────────────────────────────────────────────────────────────────────────────
# Group 2: small-table aggregations. Ordering matters:
#   growth_phases (on Experiment) feeds Publication rollup;
#   expression_status on edges feeds Experiment summary and rank_up/rank_down.
# ─────────────────────────────────────────────────────────────────────────────
echo "=== Post-process: Compute small-table aggregations ==="
time cypher-shell <<'CYPHER'
// GeneCluster member_count
MATCH (gc:GeneCluster)
OPTIONAL MATCH (gc)-[r:Gene_in_gene_cluster]->()
WITH gc, count(r) AS actual_count
SET gc.member_count = actual_count;

// Experiment growth_phases (must run before Publication rollup)
MATCH (e:Experiment)
OPTIONAL MATCH (e)-[r:Changes_expression_of]->(:Gene)
WITH e, [v IN collect(DISTINCT r.growth_phase) WHERE v IS NOT NULL] AS phases
SET e.growth_phases = phases;

// Publication summary properties
MATCH (p:Publication)
OPTIONAL MATCH (p)-[:Has_experiment]->(e:Experiment)
WITH p,
     count(e) AS ec,
     [x IN collect(DISTINCT e.omics_type) WHERE x IS NOT NULL] AS ots,
     [x IN collect(DISTINCT e.organism_name) WHERE x IS NOT NULL] AS orgs,
     [x IN collect(DISTINCT e.coculture_partner) WHERE x IS NOT NULL AND x <> ''] AS coculture_orgs,
     apoc.coll.toSet(reduce(s = [], t IN collect(coalesce(e.treatment_type, [])) | s + t)) AS tts,
     apoc.coll.toSet(reduce(s = [], t IN collect(coalesce(e.background_factors, [])) | s + t)) AS bfs,
     apoc.coll.toSet(reduce(s = [], t IN collect(coalesce(e.growth_phases, [])) | s + t)) AS gps
SET p.experiment_count = ec,
    p.treatment_types = apoc.coll.sort(tts),
    p.background_factors = apoc.coll.sort(bfs),
    p.omics_types = apoc.coll.sort(ots),
    p.growth_phases = apoc.coll.sort(gps),
    p.organisms = apoc.coll.sort(apoc.coll.toSet(orgs + coculture_orgs));

// expression_status on every Changes_expression_of edge
MATCH ()-[r:Changes_expression_of]->()
SET r.expression_status = CASE
  WHEN r.significant = 'significant' AND r.expression_direction = 'up'   THEN 'significant_up'
  WHEN r.significant = 'significant' AND r.expression_direction = 'down' THEN 'significant_down'
  ELSE 'not_significant'
END;

// Experiment stats defaults (so experiments with 0 edges still get populated props)
// gene_count:           cumulative edge count (sum across timepoints for time-course)
// distinct_gene_count:  distinct gene count (union across timepoints) — populated below
MATCH (e:Experiment)
SET e.gene_count = 0,
    e.distinct_gene_count = 0,
    e.significant_up_count = 0,
    e.significant_down_count = 0,
    e.time_point_count = 0,
    e.time_point_labels = [],
    e.time_point_orders = [],
    e.time_point_hours = [],
    e.time_point_totals = [],
    e.time_point_significant_up = [],
    e.time_point_significant_down = [],
    e.time_point_growth_phases = [];

// Experiment stats computation (overrides defaults where edges exist)
MATCH (e:Experiment)-[r:Changes_expression_of]->(g:Gene)
WITH e,
     COALESCE(r.time_point, '') AS tp,
     r.time_point_order AS tp_order,
     COALESCE(r.time_point_hours, -1.0) AS tp_hours,
     max(COALESCE(r.growth_phase, '')) AS tp_gp,
     count(r) AS total,
     count(CASE WHEN r.expression_status = 'significant_up' THEN 1 END) AS sig_up,
     count(CASE WHEN r.expression_status = 'significant_down' THEN 1 END) AS sig_down
ORDER BY e.id, tp_order
WITH e,
     sum(total) AS gene_count,
     sum(sig_up) AS significant_up_count,
     sum(sig_down) AS significant_down_count,
     collect(tp) AS tp_labels,
     collect(tp_order) AS tp_orders,
     collect(tp_hours) AS tp_hours_list,
     collect(tp_gp) AS tp_gps,
     collect(total) AS tp_totals,
     collect(sig_up) AS tp_sig_up,
     collect(sig_down) AS tp_sig_down
SET e.gene_count = gene_count,
    e.significant_up_count = significant_up_count,
    e.significant_down_count = significant_down_count,
    e.time_point_count = size(tp_labels),
    e.time_point_labels = tp_labels,
    e.time_point_orders = tp_orders,
    e.time_point_hours = tp_hours_list,
    e.time_point_totals = tp_totals,
    e.time_point_significant_up = tp_sig_up,
    e.time_point_significant_down = tp_sig_down,
    e.time_point_growth_phases = tp_gps;

// Distinct gene count per experiment (union across timepoints).
// gene_count above is cumulative; distinct_gene_count is the unique-gene
// scalar for detection-power / pathway-background reasoning.
MATCH (e:Experiment)-[:Changes_expression_of]->(g:Gene)
WITH e, count(DISTINCT g) AS dgc
SET e.distinct_gene_count = dgc;

// OrganismTaxon gene_count
MATCH (o:OrganismTaxon)
OPTIONAL MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o)
WITH o, count(g) AS gc
SET o.gene_count = gc;

// OrganismTaxon aggregate rollups from Publication
MATCH (o:OrganismTaxon)
OPTIONAL MATCH (p:Publication)
  WHERE ANY(org IN p.organisms WHERE org = o.preferred_name)
WITH o,
     count(DISTINCT p) AS pc,
     CASE WHEN count(p) > 0 THEN sum(p.experiment_count) ELSE 0 END AS ec,
     apoc.coll.toSet(reduce(s = [], t IN collect(p.treatment_types) | s + t)) AS tts,
     apoc.coll.toSet(reduce(s = [], t IN collect(p.omics_types) | s + t)) AS ots,
     apoc.coll.toSet(reduce(s = [], t IN collect(coalesce(p.background_factors, [])) | s + t)) AS bfs,
     apoc.coll.toSet(reduce(s = [], t IN collect(coalesce(p.growth_phases, [])) | s + t)) AS gps
SET o.publication_count = pc,
    o.experiment_count = ec,
    o.treatment_types = tts,
    o.omics_types = ots,
    o.background_factors = bfs,
    o.growth_phases = gps;

// ClusteringAnalysis growth_phases (from linked experiments)
MATCH (ca:ClusteringAnalysis)
OPTIONAL MATCH (e:Experiment)-[:ExperimentHasClusteringAnalysis]->(ca)
WITH ca, apoc.coll.toSet(reduce(s = [], t IN collect(coalesce(e.growth_phases, [])) | s + t)) AS gps
SET ca.growth_phases = apoc.coll.sort(gps);

// DerivedMetric total_gene_count: count of outgoing measurement edges.
// Each DM emits exactly ONE of the 3 edge types based on its value_kind,
// so the union across types is unambiguous.
MATCH (dm:DerivedMetric)
OPTIONAL MATCH (dm)-[r:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(:Gene)
WITH dm, count(r) AS total
SET dm.total_gene_count = total;

// DerivedMetric growth_phases: union from parent Experiment
// (mirrors ClusteringAnalysis growth_phases; reads Experiment.growth_phases
// set earlier in Group 2).
MATCH (dm:DerivedMetric)
OPTIONAL MATCH (e:Experiment)-[:ExperimentHasDerivedMetric]->(dm)
WITH dm, apoc.coll.toSet(reduce(s = [], t IN collect(coalesce(e.growth_phases, [])) | s + t)) AS gps
SET dm.growth_phases = apoc.coll.sort(gps);

// DerivedMetric numeric distribution stats (value_kind='numeric' only).
// Aggregates over Derived_metric_quantifies_gene.value. Cypher percentileCont
// is the reference interpolation method consumed by explorer queries.
// Numeric DMs without quantifies edges (theoretical) leave the props null.
MATCH (dm:DerivedMetric {value_kind: 'numeric'})-[r:Derived_metric_quantifies_gene]->(:Gene)
WITH dm,
     min(r.value)                  AS v_min,
     max(r.value)                  AS v_max,
     percentileCont(r.value, 0.25) AS v_q1,
     percentileCont(r.value, 0.5)  AS v_median,
     percentileCont(r.value, 0.75) AS v_q3
SET dm.value_min    = v_min,
    dm.value_max    = v_max,
    dm.value_q1     = v_q1,
    dm.value_median = v_median,
    dm.value_q3     = v_q3;

// DerivedMetric boolean flag counts (value_kind='boolean' only).
// Counts true vs false values on Derived_metric_flags_gene edges.
// Booleans without edges get 0/0 (defaults handled by COALESCE pattern).
MATCH (dm:DerivedMetric {value_kind: 'boolean'})
OPTIONAL MATCH (dm)-[r:Derived_metric_flags_gene]->(:Gene)
WITH dm,
     count(CASE WHEN r.value = 'true'  THEN 1 END) AS n_true,
     count(CASE WHEN r.value = 'false' THEN 1 END) AS n_false
SET dm.flag_true_count  = n_true,
    dm.flag_false_count = n_false;

// DerivedMetric categorical distribution (value_kind='categorical' only).
// Parallel arrays sorted by label so output is deterministic.
// Categoricals without edges leave label/count arrays empty.
MATCH (dm:DerivedMetric {value_kind: 'categorical'})
SET dm.category_labels = [], dm.category_counts = [];
MATCH (dm:DerivedMetric {value_kind: 'categorical'})-[r:Derived_metric_classifies_gene]->(:Gene)
WITH dm, r.value AS cat, count(r) AS cnt
ORDER BY cat
WITH dm, collect(cat) AS labels, collect(cnt) AS counts
SET dm.category_labels = labels,
    dm.category_counts = counts;

// OrganismTaxon clustering rollup
MATCH (o:OrganismTaxon)
OPTIONAL MATCH (ca:ClusteringAnalysis)-[:ClusteringanalysisBelongsToOrganism]->(o)
WITH o,
     count(ca) AS ca_count,
     collect(DISTINCT ca.cluster_type) AS ctypes,
     sum(coalesce(ca.cluster_count, 0)) AS total_clusters
SET o.clustering_analysis_count = ca_count,
    o.cluster_types = ctypes,
    o.cluster_count = total_clusters;

// Publication clustering rollup
MATCH (p:Publication)
OPTIONAL MATCH (p)-[:PublicationHasClusteringAnalysis]->(ca:ClusteringAnalysis)
WITH p,
     count(ca) AS ca_count,
     collect(DISTINCT ca.cluster_type) AS ctypes,
     sum(coalesce(ca.cluster_count, 0)) AS total_clusters
SET p.clustering_analysis_count = ca_count,
    p.cluster_types = ctypes,
    p.cluster_count = total_clusters;

// Publication DM rollup defaults
MATCH (p:Publication)
SET p.derived_metric_count = 0,
    p.derived_metric_gene_count = 0,
    p.compartments = [],
    p.derived_metric_types = [],
    p.derived_metric_value_kinds = [];

// Publication DM compute
MATCH (p:Publication)
OPTIONAL MATCH (p)-[:PublicationHasDerivedMetric]->(dm:DerivedMetric)
WITH p,
     count(DISTINCT dm) AS dm_count,
     [x IN collect(DISTINCT dm.metric_type) WHERE x IS NOT NULL] AS metric_types,
     [x IN collect(DISTINCT dm.value_kind) WHERE x IS NOT NULL] AS value_kinds
SET p.derived_metric_count = dm_count,
    p.derived_metric_types = apoc.coll.sort(metric_types),
    p.derived_metric_value_kinds = apoc.coll.sort(value_kinds);

// Publication compartments: from child Experiments
MATCH (p:Publication)
OPTIONAL MATCH (p)-[:Has_experiment]->(e:Experiment)
WITH p, [x IN collect(DISTINCT e.compartment) WHERE x IS NOT NULL] AS comps
SET p.compartments = apoc.coll.sort(comps);

// Publication derived_metric_gene_count
MATCH (p:Publication)
OPTIONAL MATCH (p)-[:PublicationHasDerivedMetric]->(:DerivedMetric)
  -[:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(g:Gene)
WITH p, count(DISTINCT g) AS dmg_count
SET p.derived_metric_gene_count = dmg_count;

// Publication derived_metric_search_text: aggregated DM tokens for fulltext discovery.
// Tokens: dm.name + dm.metric_type (underscore -> space) + dm.field_description.
// compartment is indexed separately via p.compartments (already computed above).
// Stored as null when no DMs reachable so the fulltext index skips the node.
MATCH (p:Publication)-[:PublicationHasDerivedMetric]->(dm:DerivedMetric)
WITH p,
     [x IN collect(DISTINCT dm.name) WHERE x IS NOT NULL]              AS names,
     [x IN collect(DISTINCT dm.metric_type) WHERE x IS NOT NULL]       AS metric_types,
     [x IN collect(DISTINCT dm.field_description) WHERE x IS NOT NULL] AS descs
SET p.derived_metric_search_text = trim(
      apoc.text.join(names, ' ') + ' '
    + apoc.text.replace(apoc.text.join(metric_types, ' '), '_', ' ') + ' '
    + apoc.text.join(descs, ' ')
);

// OrganismTaxon DM rollup defaults
MATCH (o:OrganismTaxon)
SET o.derived_metric_count = 0,
    o.derived_metric_gene_count = 0,
    o.compartments = [],
    o.derived_metric_types = [],
    o.derived_metric_value_kinds = [];

// OrganismTaxon DM compute
MATCH (o:OrganismTaxon)
OPTIONAL MATCH (dm:DerivedMetric)-[:DerivedMetricBelongsToOrganism]->(o)
WITH o,
     count(DISTINCT dm) AS dm_count,
     [x IN collect(DISTINCT dm.metric_type) WHERE x IS NOT NULL] AS metric_types,
     [x IN collect(DISTINCT dm.value_kind) WHERE x IS NOT NULL] AS value_kinds,
     [x IN collect(DISTINCT dm.compartment) WHERE x IS NOT NULL] AS comps
SET o.derived_metric_count = dm_count,
    o.derived_metric_types = apoc.coll.sort(metric_types),
    o.derived_metric_value_kinds = apoc.coll.sort(value_kinds),
    o.compartments = apoc.coll.sort(comps);

// OrganismTaxon derived_metric_gene_count
MATCH (o:OrganismTaxon)
OPTIONAL MATCH (dm:DerivedMetric)-[:DerivedMetricBelongsToOrganism]->(o)
OPTIONAL MATCH (dm)-[:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(g:Gene)
WITH o, count(DISTINCT g) AS dmg_count
SET o.derived_metric_gene_count = dmg_count;

// Experiment clustering rollup
MATCH (e:Experiment)
OPTIONAL MATCH (e)-[:ExperimentHasClusteringAnalysis]->(ca:ClusteringAnalysis)
WITH e,
     count(ca) AS ca_count,
     collect(DISTINCT ca.cluster_type) AS ctypes,
     sum(coalesce(ca.cluster_count, 0)) AS total_clusters
SET e.clustering_analysis_count = ca_count,
    e.cluster_types = ctypes,
    e.cluster_count = total_clusters;

// Experiment DM rollup defaults (empty-state; compute below overrides where children exist)
MATCH (e:Experiment)
SET e.reports_fold_change = 'false',
    e.reports_derived_metric_types = [],
    e.derived_metric_count = 0,
    e.derived_metric_value_kinds = [],
    e.derived_metric_gene_count = 0;

// Experiment reports_fold_change: 'true' iff outgoing Changes_expression_of exists
MATCH (e:Experiment)
WHERE EXISTS { (e)-[:Changes_expression_of]->() }
SET e.reports_fold_change = 'true';

// Experiment DM compute (overrides defaults)
MATCH (e:Experiment)
OPTIONAL MATCH (e)-[:ExperimentHasDerivedMetric]->(dm:DerivedMetric)
WITH e,
     count(DISTINCT dm) AS dm_count,
     [x IN collect(DISTINCT dm.metric_type) WHERE x IS NOT NULL] AS metric_types,
     [x IN collect(DISTINCT dm.value_kind) WHERE x IS NOT NULL] AS value_kinds
SET e.derived_metric_count = dm_count,
    e.reports_derived_metric_types = apoc.coll.sort(metric_types),
    e.derived_metric_value_kinds = apoc.coll.sort(value_kinds);

// Experiment derived_metric_gene_count: distinct genes reachable via ANY child DM edge type
MATCH (e:Experiment)
OPTIONAL MATCH (e)-[:ExperimentHasDerivedMetric]->(:DerivedMetric)
  -[:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(g:Gene)
WITH e, count(DISTINCT g) AS dmg_count
SET e.derived_metric_gene_count = dmg_count;

// Experiment derived_metric_search_text: same aggregation shape as Publication.
// compartment is indexed separately via e.compartment (adapter-emitted).
// Stored as null when no DMs reachable.
MATCH (e:Experiment)-[:ExperimentHasDerivedMetric]->(dm:DerivedMetric)
WITH e,
     [x IN collect(DISTINCT dm.name) WHERE x IS NOT NULL]              AS names,
     [x IN collect(DISTINCT dm.metric_type) WHERE x IS NOT NULL]       AS metric_types,
     [x IN collect(DISTINCT dm.field_description) WHERE x IS NOT NULL] AS descs
SET e.derived_metric_search_text = trim(
      apoc.text.join(names, ' ') + ' '
    + apoc.text.replace(apoc.text.join(metric_types, ' '), '_', ' ') + ' '
    + apoc.text.join(descs, ' ')
);
CYPHER

# ─────────────────────────────────────────────────────────────────────────────
# Group 3: heavy Gene/edge/BriteCategory writes using `CALL { } IN TRANSACTIONS`.
# Each statement gets its own implicit transaction (cypher-shell default for
# stdin, and required for IN TRANSACTIONS).
# ─────────────────────────────────────────────────────────────────────────────
echo "=== Post-process: Compute Gene routing signals, ranks, BriteCategory ==="
time cypher-shell <<'CYPHER'
// annotation_types
MATCH (g:Gene)
CALL {
  WITH g
  SET g.annotation_types =
    CASE WHEN EXISTS { (g)-[:Gene_involved_in_biological_process]->() } THEN ['go_bp'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_enables_molecular_function]->() } THEN ['go_mf'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_located_in_cellular_component]->() } THEN ['go_cc'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_pfam]->() } THEN ['pfam'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_in_cog_category]->() } THEN ['cog_category'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_kegg_ko]->() } THEN ['kegg'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_kegg_ko]->()-[:Kegg_term_in_brite_category]->() } THEN ['brite'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_catalyzes_ec_number]->() } THEN ['ec'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_cyanorak_role]->() } THEN ['cyanorak_role'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_tigr_role]->() } THEN ['tigr_role'] ELSE [] END
} IN TRANSACTIONS OF 1000 ROWS;

// expression_edge_count + significant_up/down_count
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)<-[e:Changes_expression_of]-()
  WITH g, count(e) AS total,
       sum(CASE WHEN e.expression_status = 'significant_up' THEN 1 ELSE 0 END) AS sig_up,
       sum(CASE WHEN e.expression_status = 'significant_down' THEN 1 ELSE 0 END) AS sig_down
  SET g.expression_edge_count = total,
      g.significant_up_count = sig_up,
      g.significant_down_count = sig_down
} IN TRANSACTIONS OF 500 ROWS;

// rank_by_effect
MATCH (e:Experiment)
WITH e
CALL {
  WITH e
  MATCH (e)-[r:Changes_expression_of]->(g:Gene)
  WITH r.time_point_order AS tp, r, abs(r.log2_fold_change) AS abs_fc,
       coalesce(r.adjusted_p_value, 2.0) AS padj, g.locus_tag AS lt
  ORDER BY tp, abs_fc DESC, padj ASC, lt ASC
  WITH tp, collect(r) AS edges
  UNWIND range(0, size(edges)-1) AS i
  SET (edges[i]).rank_by_effect = i + 1
} IN TRANSACTIONS OF 10 ROWS;

// rank_up (among significant_up per experiment x timepoint)
MATCH (e:Experiment)
WITH e
CALL {
  WITH e
  MATCH (e)-[r:Changes_expression_of]->(g:Gene)
  WHERE r.expression_status = 'significant_up'
  WITH r.time_point_order AS tp, r, abs(r.log2_fold_change) AS abs_fc,
       coalesce(r.adjusted_p_value, 2.0) AS padj, g.locus_tag AS lt
  ORDER BY tp, abs_fc DESC, padj ASC, lt ASC
  WITH tp, collect(r) AS edges
  UNWIND range(0, size(edges)-1) AS i
  SET (edges[i]).rank_up = i + 1
} IN TRANSACTIONS OF 10 ROWS;

// rank_down (among significant_down per experiment x timepoint)
MATCH (e:Experiment)
WITH e
CALL {
  WITH e
  MATCH (e)-[r:Changes_expression_of]->(g:Gene)
  WHERE r.expression_status = 'significant_down'
  WITH r.time_point_order AS tp, r, abs(r.log2_fold_change) AS abs_fc,
       coalesce(r.adjusted_p_value, 2.0) AS padj, g.locus_tag AS lt
  ORDER BY tp, abs_fc DESC, padj ASC, lt ASC
  WITH tp, collect(r) AS edges
  UNWIND range(0, size(edges)-1) AS i
  SET (edges[i]).rank_down = i + 1
} IN TRANSACTIONS OF 10 ROWS;

// Numeric DM rank/percentile/bucket: ranks derived_metric_quantifies_gene edges
// grouped by DerivedMetric, only when parent DM has rankable='true'. Ties on
// value broken by Gene.locus_tag ascending (reproducibility).
// Percentile: rank 1 (highest value) -> 100.0; rank N (lowest) -> 0.0.
// Buckets pinned per slice spec §Post-import (thresholds must not drift).
MATCH (dm:DerivedMetric {rankable: 'true'})
CALL {
  WITH dm
  MATCH (dm)-[r:Derived_metric_quantifies_gene]->(g:Gene)
  WITH r, r.value AS val, g.locus_tag AS lt
  ORDER BY val DESC, lt ASC
  WITH collect(r) AS edges, count(r) AS n
  UNWIND range(0, size(edges) - 1) AS i
  WITH edges[i] AS r, i, n,
       CASE WHEN n = 1 THEN 100.0
            ELSE 100.0 * toFloat(n - i - 1) / toFloat(n - 1)
       END AS pct
  SET r.rank_by_metric = i + 1,
      r.metric_percentile = pct,
      r.metric_bucket = CASE
        WHEN pct >= 90.0 THEN 'top_decile'
        WHEN pct >= 75.0 THEN 'top_quartile'
        WHEN pct >= 25.0 THEN 'mid'
        ELSE 'low'
      END
} IN TRANSACTIONS OF 10 ROWS;

// Numeric DM significance: on derived_metric_quantifies_gene edges,
// only when parent DM has has_p_value='true' AND p_value_threshold IS NOT NULL
// AND the edge's adjusted_p_value is non-null. Left null otherwise.
MATCH (dm:DerivedMetric {has_p_value: 'true'})
WHERE dm.p_value_threshold IS NOT NULL
CALL {
  WITH dm
  MATCH (dm)-[r:Derived_metric_quantifies_gene]->()
  WHERE r.adjusted_p_value IS NOT NULL
  SET r.significant = CASE
    WHEN r.adjusted_p_value < dm.p_value_threshold THEN 'true'
    ELSE 'false'
  END
} IN TRANSACTIONS OF 1000 ROWS;

// Gene DM routing defaults + compute.
// Defaults run first (every gene gets 0 / []), then 4 OPTIONAL-MATCH passes
// override defaults only on genes reached by DM edges.

// Defaults
MATCH (g:Gene)
CALL {
  WITH g
  SET g.numeric_metric_count = 0,
      g.boolean_metric_count = 0,
      g.categorical_metric_count = 0,
      g.numeric_metric_types_observed = [],
      g.boolean_metric_types_observed = [],
      g.categorical_metric_types_observed = [],
      g.compartments_observed = []
} IN TRANSACTIONS OF 1000 ROWS;

// numeric_metric_count + numeric_metric_types_observed
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (dm:DerivedMetric)-[:Derived_metric_quantifies_gene]->(g)
  WITH g,
       count(DISTINCT dm) AS cnt,
       [x IN collect(DISTINCT dm.metric_type) WHERE x IS NOT NULL] AS types
  SET g.numeric_metric_count = cnt,
      g.numeric_metric_types_observed = apoc.coll.sort(types)
} IN TRANSACTIONS OF 1000 ROWS;

// boolean_metric_count + boolean_metric_types_observed
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (dm:DerivedMetric)-[:Derived_metric_flags_gene]->(g)
  WITH g,
       count(DISTINCT dm) AS cnt,
       [x IN collect(DISTINCT dm.metric_type) WHERE x IS NOT NULL] AS types
  SET g.boolean_metric_count = cnt,
      g.boolean_metric_types_observed = apoc.coll.sort(types)
} IN TRANSACTIONS OF 1000 ROWS;

// categorical_metric_count + categorical_metric_types_observed
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (dm:DerivedMetric)-[:Derived_metric_classifies_gene]->(g)
  WITH g,
       count(DISTINCT dm) AS cnt,
       [x IN collect(DISTINCT dm.metric_type) WHERE x IS NOT NULL] AS types
  SET g.categorical_metric_count = cnt,
      g.categorical_metric_types_observed = apoc.coll.sort(types)
} IN TRANSACTIONS OF 1000 ROWS;

// compartments_observed: union across all 3 DM edge types via parent DerivedMetric
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (dm:DerivedMetric)
    -[:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(g)
  WITH g, [x IN collect(DISTINCT dm.compartment) WHERE x IS NOT NULL] AS comps
  SET g.compartments_observed = apoc.coll.sort(comps)
} IN TRANSACTIONS OF 1000 ROWS;

// closest_ortholog_group_size + closest_ortholog_genera
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[:Gene_in_ortholog_group]->(og:OrthologGroup)
  WITH g, og ORDER BY og.specificity_rank ASC LIMIT 1
  SET g.closest_ortholog_group_size = og.member_count,
      g.closest_ortholog_genera = og.genera
} IN TRANSACTIONS OF 1000 ROWS;

// cluster_membership_count + cluster_types
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (gc:GeneCluster)-[:Gene_in_gene_cluster]->(g)
  OPTIONAL MATCH (ca:ClusteringAnalysis)-[:ClusteringAnalysisHasGeneCluster]->(gc)
  WITH g,
       count(DISTINCT gc) AS membership_count,
       collect(DISTINCT ca.cluster_type) AS ctypes
  SET g.cluster_membership_count = membership_count,
      g.cluster_types = CASE WHEN size(ctypes) = 0 THEN [] ELSE ctypes END
} IN TRANSACTIONS OF 1000 ROWS;

// BriteCategory computed properties: member_ko_count, gene_count, organism_count
MATCH (b:BriteCategory)
CALL {
  WITH b
  // member_ko_count: KO leaves directly or transitively under this category
  MATCH (b)<-[:Brite_category_is_a_brite_category*0..]-(leaf:BriteCategory)
  OPTIONAL MATCH (ko:KeggTerm)-[:Kegg_term_in_brite_category]->(leaf)
  WHERE ko.level_kind = 'ko'
  WITH b, count(DISTINCT ko) AS ko_count
  // gene_count and organism_count via KO→gene edges
  OPTIONAL MATCH (b)<-[:Brite_category_is_a_brite_category*0..]-(leaf2:BriteCategory)
  OPTIONAL MATCH (ko2:KeggTerm)-[:Kegg_term_in_brite_category]->(leaf2)
  WHERE ko2.level_kind = 'ko'
  OPTIONAL MATCH (g:Gene)-[:Gene_has_kegg_ko]->(ko2)
  WITH b, ko_count, count(DISTINCT g) AS g_count, collect(DISTINCT g.organism_name) AS orgs
  SET b.member_ko_count = ko_count,
      b.gene_count = g_count,
      b.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 100 ROWS;

// ── TcdbFamily computed properties ───────────────────────────────────────────

// member_count: direct child count
MATCH (t:TcdbFamily)
CALL {
  WITH t
  OPTIONAL MATCH (child:TcdbFamily)-[:Tcdb_family_is_a_tcdb_family]->(t)
  WITH t, count(child) AS mc
  SET t.member_count = mc
} IN TRANSACTIONS OF 1000 ROWS;

// tc_class_id: pointer to the root tc_class node
// Self-reference for class nodes:
MATCH (t:TcdbFamily {level_kind: 'tc_class'})
SET t.tc_class_id = t.id;

// Non-class nodes: walk up to nearest tc_class ancestor
MATCH (t:TcdbFamily) WHERE t.level_kind <> 'tc_class'
CALL {
  WITH t
  MATCH (t)-[:Tcdb_family_is_a_tcdb_family*1..]->(cls:TcdbFamily {level_kind: 'tc_class'})
  WITH t, cls LIMIT 1
  SET t.tc_class_id = cls.id
} IN TRANSACTIONS OF 1000 ROWS;

// gene_count + organism_count: subtree traversal (descendants ∪ self via *0..)
MATCH (t:TcdbFamily)
CALL {
  WITH t
  OPTIONAL MATCH (t)<-[:Tcdb_family_is_a_tcdb_family*0..]-(desc:TcdbFamily)<-[:Gene_has_tcdb_family]-(g:Gene)
  WITH t, count(DISTINCT g) AS gc, collect(DISTINCT g.organism_name) AS orgs
  SET t.gene_count = gc,
      t.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 1000 ROWS;

// metabolite_count: distinct metabolites reachable via Tcdb_family_transports_metabolite in subtree
MATCH (t:TcdbFamily)
CALL {
  WITH t
  OPTIONAL MATCH (t)<-[:Tcdb_family_is_a_tcdb_family*0..]-(desc:TcdbFamily)-[:Tcdb_family_transports_metabolite]->(m:Metabolite)
  WITH t, count(DISTINCT m) AS mc
  SET t.metabolite_count = mc
} IN TRANSACTIONS OF 1000 ROWS;

// ── CazyFamily computed properties ───────────────────────────────────────────

// gene_count + organism_count: subtree traversal (descendants ∪ self via *0..)
MATCH (c:CazyFamily)
CALL {
  WITH c
  OPTIONAL MATCH (c)<-[:Cazy_family_is_a_cazy_family*0..]-(desc:CazyFamily)<-[:Gene_has_cazy_family]-(g:Gene)
  WITH c, count(DISTINCT g) AS gc, collect(DISTINCT g.organism_name) AS orgs
  SET c.gene_count = gc,
      c.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 1000 ROWS;

// ── Gene routing extensions ──────────────────────────────────────────────────

// Extend annotation_types — add 'tcdb' / 'cazy' if any edge of that type exists.
// (The base annotation_types is computed earlier in this same group; we APPEND
//  here rather than rebuild from scratch so the existing values persist.)
MATCH (g:Gene)
CALL {
  WITH g
  WITH g, coalesce(g.annotation_types, []) AS existing
  SET g.annotation_types = existing +
    CASE WHEN EXISTS { (g)-[:Gene_has_tcdb_family]->() } THEN ['tcdb'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_cazy_family]->() } THEN ['cazy'] ELSE [] END
} IN TRANSACTIONS OF 1000 ROWS;

// tcdb_family_count + cazy_family_count single-hop rollups (per TCDB-S1 / TCDB-S2)
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[r1:Gene_has_tcdb_family]->()
  WITH g, count(r1) AS tc
  OPTIONAL MATCH (g)-[r2:Gene_has_cazy_family]->()
  WITH g, tc, count(r2) AS cz
  SET g.tcdb_family_count = tc,
      g.cazy_family_count = cz
} IN TRANSACTIONS OF 1000 ROWS;

// Gene.metabolite_count = UNION across catalysis + transport paths (per TCDB-S3)
// Catalysis arm: Gene -> Reaction -> Metabolite (existing chemistry).
// Transport arm: Gene -> TcdbFamily -> ... -> tc_specificity leaf -> Metabolite.
// Uses apoc.coll.toSet for distinct UNION across the two paths.
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[:Gene_catalyzes_reaction]->(:Reaction)-[:Reaction_has_metabolite]->(m1:Metabolite)
  WITH g, collect(DISTINCT m1) AS m_cat
  OPTIONAL MATCH (g)-[:Gene_has_tcdb_family]->(:TcdbFamily)
                 <-[:Tcdb_family_is_a_tcdb_family*0..]-(:TcdbFamily {level_kind: 'tc_specificity'})
                 -[:Tcdb_family_transports_metabolite]->(m2:Metabolite)
  WITH g, m_cat + collect(DISTINCT m2) AS all_m
  SET g.metabolite_count = size(apoc.coll.toSet(all_m))
} IN TRANSACTIONS OF 500 ROWS;

// ── Metabolism rollups ────────────────────────────────────────────────────

// Reaction.gene_count, organism_count, organisms[]
CALL {
  MATCH (r:Reaction)<-[:Gene_catalyzes_reaction]-(g:Gene)
  WITH r, count(DISTINCT g) AS gene_count, collect(DISTINCT g.organism_name) AS organisms
  SET r.gene_count = gene_count,
      r.organism_count = size(organisms),
      r.organisms = organisms
} IN TRANSACTIONS OF 1000 ROWS;

// Metabolite.gene_count, organism_count, transporter_count
// (UNION-aware: catalysis + transport paths)
CALL {
  MATCH (m:Metabolite)
  OPTIONAL MATCH (m)<-[:Reaction_has_metabolite]-(:Reaction)<-[:Gene_catalyzes_reaction]-(g_cat:Gene)
  WITH m, collect(DISTINCT g_cat) AS gs_cat
  OPTIONAL MATCH (m)<-[:Tcdb_family_transports_metabolite]-(:TcdbFamily)
                 <-[:Gene_has_tcdb_family]-(g_tr:Gene)
  WITH m, gs_cat + collect(DISTINCT g_tr) AS all_g
  WITH m,
       size(apoc.coll.toSet(all_g)) AS gc,
       size(apoc.coll.toSet([x IN all_g | x.organism_name])) AS oc
  SET m.gene_count = gc, m.organism_count = oc
} IN TRANSACTIONS OF 1000 ROWS;

// Metabolite.transporter_count: distinct TcdbFamily ancestors with substrate edge
CALL {
  MATCH (m:Metabolite)
  OPTIONAL MATCH (t:TcdbFamily)-[:Tcdb_family_transports_metabolite]->(m)
  WITH m, count(DISTINCT t) AS tc
  SET m.transporter_count = tc
} IN TRANSACTIONS OF 1000 ROWS;

// Materialize Organism_has_metabolite (2-hop save)
CALL {
  MATCH (o:OrganismTaxon)<-[:Gene_belongs_to_organism]-(g:Gene)
        -[:Gene_catalyzes_reaction]->(:Reaction)
        -[:Reaction_has_metabolite]->(m:Metabolite)
  WITH DISTINCT o, m
  MERGE (o)-[:Organism_has_metabolite]->(m)
} IN TRANSACTIONS OF 1000 ROWS;

// Materialize Organism_has_metabolite (transport arm)
CALL {
  MATCH (o:OrganismTaxon)<-[:Gene_belongs_to_organism]-(g:Gene)
        -[:Gene_has_tcdb_family]->(:TcdbFamily)
        <-[:Tcdb_family_is_a_tcdb_family*0..]-(:TcdbFamily {level_kind: 'tc_specificity'})
        -[:Tcdb_family_transports_metabolite]->(m:Metabolite)
  WITH DISTINCT o, m
  MERGE (o)-[:Organism_has_metabolite]->(m)
} IN TRANSACTIONS OF 1000 ROWS;

// Organism rollup props
CALL {
  MATCH (o:OrganismTaxon)-[:Organism_has_metabolite]->(m:Metabolite)
  WITH o, count(DISTINCT m) AS metabolite_count
  SET o.metabolite_count = metabolite_count
} IN TRANSACTIONS OF 100 ROWS;

CALL {
  MATCH (o:OrganismTaxon)<-[:Gene_belongs_to_organism]-(:Gene)
        -[:Gene_catalyzes_reaction]->(r:Reaction)
  WITH o, count(DISTINCT r) AS reaction_count
  SET o.reaction_count = reaction_count
} IN TRANSACTIONS OF 100 ROWS;

// ── Chemistry slice-1 rollups (KG-A1, KG-A2, KG-A4) ───────────────────────

// KG-A1: Gene.reaction_count -- single-hop count of catalysis edges.
// count(r) returns 0 cleanly on no-match; no defaults pass needed.
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[r:Gene_catalyzes_reaction]->()
  WITH g, count(r) AS rxn_count
  SET g.reaction_count = rxn_count
} IN TRANSACTIONS OF 1000 ROWS;

// KG-A2: Gene.metabolite_count -- defined as UNION across catalysis + transport
// paths. The catalysis arm + transport arm are both implemented earlier in
// this script (see "Gene.metabolite_count = UNION across catalysis + transport"
// block above, post-TCDB-CAZy landing). Chemistry-slice-1's placeholder
// catalysis-only block has been dropped — our UNION subsumes it.

// KG-A4: KeggTerm pathway-level rollups (sparse on pathways only).
// level_kind = 'pathway' filter; KOs / categories left unset.
MATCH (p:KeggTerm) WHERE p.level_kind = 'pathway'
CALL {
  WITH p
  OPTIONAL MATCH (r:Reaction)-[:Reaction_in_kegg_pathway]->(p)
  WITH p, count(r) AS rxn_count
  SET p.reaction_count = rxn_count
} IN TRANSACTIONS OF 100 ROWS;

MATCH (p:KeggTerm) WHERE p.level_kind = 'pathway'
CALL {
  WITH p
  OPTIONAL MATCH (m:Metabolite)-[:Metabolite_in_pathway]->(p)
  WITH p, count(m) AS met_count
  SET p.metabolite_count = met_count
} IN TRANSACTIONS OF 100 ROWS;
CYPHER

echo "=== Post-process complete ==="
neo4j stop
