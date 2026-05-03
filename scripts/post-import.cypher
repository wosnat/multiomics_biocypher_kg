// Post-import Cypher commands
// These run after the knowledge graph is imported into Neo4j
// Indexes only — homolog edges and expression propagation have been replaced
// by OrthologGroup nodes and Gene_in_ortholog_group edges (query-time joins).

// Scalar indexes for get_gene exact lookup
CREATE INDEX gene_locus_tag_idx IF NOT EXISTS FOR (g:Gene) ON (g.locus_tag);
CREATE INDEX gene_name_idx IF NOT EXISTS FOR (g:Gene) ON (g.gene_name);
CREATE INDEX gene_organism_name_idx IF NOT EXISTS FOR (g:Gene) ON (g.organism_name);

// Full-text index for find_gene free-text search
CREATE FULLTEXT INDEX geneFullText IF NOT EXISTS FOR (n:Gene) ON EACH [
  n.gene_summary, n.all_identifiers, n.gene_name_synonyms,
  n.alternate_functional_descriptions];

// Full-text indexes for ontology term search (search_ontology tool)
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

// OrthologGroup indexes for efficient 2-hop homolog lookups
CREATE INDEX ortholog_group_id_idx IF NOT EXISTS FOR (og:OrthologGroup) ON (og.id);
CREATE INDEX ortholog_group_name_idx IF NOT EXISTS FOR (og:OrthologGroup) ON (og.name);
CREATE INDEX ortholog_group_level_idx IF NOT EXISTS FOR (og:OrthologGroup) ON (og.taxonomic_level);
CREATE INDEX ortholog_group_rank_idx IF NOT EXISTS FOR (og:OrthologGroup) ON (og.specificity_rank);

CREATE FULLTEXT INDEX orthologGroupFullText IF NOT EXISTS
  FOR (og:OrthologGroup) ON EACH [og.consensus_product, og.consensus_gene_name, og.description, og.functional_description];

// Pfam domain and clan indexes
CREATE INDEX pfam_name_idx IF NOT EXISTS FOR (p:Pfam) ON (p.name);
CREATE INDEX pfam_clan_name_idx IF NOT EXISTS FOR (c:PfamClan) ON (c.name);

CREATE FULLTEXT INDEX pfamFullText IF NOT EXISTS
  FOR (p:Pfam) ON EACH [p.name, p.short_name];
CREATE FULLTEXT INDEX pfamClanFullText IF NOT EXISTS
  FOR (c:PfamClan) ON EACH [c.name];

// ── BriteCategory indexes ───────────────────────────────────────────────
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

// Publication fulltext index
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

// Experiment indexes
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

// ── OrganismTaxon indexes ──────────────────────────────────────────────
CREATE INDEX organism_type_idx IF NOT EXISTS FOR (o:OrganismTaxon) ON (o.organism_type);

// ── ClusteringAnalysis indexes ──────────────────────────────────────────
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

// ── GeneCluster indexes ─────────────────────────────────────────────────
CREATE FULLTEXT INDEX geneClusterFullText IF NOT EXISTS
  FOR (gc:GeneCluster) ON EACH [gc.name, gc.functional_description, gc.temporal_pattern, gc.expression_dynamics];

// KeggTerm
CREATE INDEX kegg_term_id_idx IF NOT EXISTS FOR (k:KeggTerm) ON (k.id);

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
CREATE INDEX metabolite_hmdb_idx IF NOT EXISTS FOR (m:Metabolite) ON (m.hmdb_id);
CREATE FULLTEXT INDEX metaboliteFullText IF NOT EXISTS FOR (m:Metabolite) ON EACH [m.name];

// -----------------------------------------------------------------------
// ── GeneCluster member_count verification ──────────────────────────────
MATCH (gc:GeneCluster)
OPTIONAL MATCH (gc)-[r:Gene_in_gene_cluster]->()
WITH gc, count(r) AS actual_count
SET gc.member_count = actual_count;

// Experiment growth_phases (must run before Publication rollup uses it)
// -----------------------------------------------------------------------

MATCH (e:Experiment)
OPTIONAL MATCH (e)-[r:Changes_expression_of]->(:Gene)
WITH e, [v IN collect(DISTINCT r.growth_phase) WHERE v IS NOT NULL] AS phases
SET e.growth_phases = phases;

// -----------------------------------------------------------------------
// Publication summary properties (pre-computed for list_publications)
// -----------------------------------------------------------------------

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

// -----------------------------------------------------------------------
// expression_status on edges (derived from significant + expression_direction)
// -----------------------------------------------------------------------

MATCH ()-[r:Changes_expression_of]->()
SET r.expression_status = CASE
  WHEN r.significant = 'significant' AND r.expression_direction = 'up'   THEN 'significant_up'
  WHEN r.significant = 'significant' AND r.expression_direction = 'down' THEN 'significant_down'
  ELSE 'not_significant'
END;

// -----------------------------------------------------------------------
// Experiment summary properties (pre-computed for list_experiments)
// -----------------------------------------------------------------------

// Pass 1: set defaults for all experiments (avoids NULL from OPTIONAL MATCH)
// gene_count:           cumulative edge count (sum across timepoints for time-course)
// distinct_gene_count:  distinct gene count (union across timepoints) — see Pass 3 below
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

// Pass 2: compute actual stats for experiments with expression edges
// Neo4j cannot store nulls in arrays, so we COALESCE:
//   time_point labels: null → "" (non-time-course experiments)
//   time_point hours:  null → -1.0 (unknown conversion)
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

// Pass 3: distinct gene count per experiment (union across timepoints).
// gene_count above is cumulative (sums per-timepoint counts), so for a
// time-course experiment that measures the same gene set at every TP,
// gene_count = distinct_gene_count * time_point_count. Researchers reasoning
// about detection power / pathway-background size need distinct_gene_count,
// not gene_count.
MATCH (e:Experiment)-[:Changes_expression_of]->(g:Gene)
WITH e, count(DISTINCT g) AS dgc
SET e.distinct_gene_count = dgc;

// -----------------------------------------------------------------------
// OrganismTaxon summary properties (pre-computed for list_organisms)
// -----------------------------------------------------------------------

// gene_count: count of Gene_belongs_to_organism edges
MATCH (o:OrganismTaxon)
OPTIONAL MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o)
WITH o, count(g) AS gc
SET o.gene_count = gc;

// publication_count, experiment_count, treatment_types, omics_types, background_factors
// (depends on Publication.organisms being computed first)
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

// ── ClusteringAnalysis summary properties ─────────────────────────────

// ClusteringAnalysis: growth_phases (union of linked Experiment.growth_phases)
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
// Booleans without edges get 0/0.
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

// OrganismTaxon: clustering_analysis_count, cluster_types, cluster_count
MATCH (o:OrganismTaxon)
OPTIONAL MATCH (ca:ClusteringAnalysis)-[:ClusteringanalysisBelongsToOrganism]->(o)
WITH o,
     count(ca) AS ca_count,
     collect(DISTINCT ca.cluster_type) AS ctypes,
     sum(coalesce(ca.cluster_count, 0)) AS total_clusters
SET o.clustering_analysis_count = ca_count,
    o.cluster_types = ctypes,
    o.cluster_count = total_clusters;

// Publication: clustering_analysis_count, cluster_types, cluster_count
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

// Experiment: clustering_analysis_count, cluster_types, cluster_count
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

// =====================================================================
// F1.1: Term-level is_uninformative flag (sentinel str 'true' / absent).
// Driven by config/uninformative_terms.yaml. Vocabulary is small; keys
// are hard-coded here matching the YAML section names.
//
// Guiding principle: only flag terms that convey no class signal at all.
// Pfam DUF/UPF, COG R, all BRITE entries, all EC numbers stay UN-flagged.
// =====================================================================

// Direct ID flags
MATCH (t:BiologicalProcess) WHERE t.id IN ['go:0008150'] SET t.is_uninformative = 'true';
MATCH (t:MolecularFunction) WHERE t.id IN ['go:0003674'] SET t.is_uninformative = 'true';
MATCH (t:CellularComponent) WHERE t.id IN ['go:0005575'] SET t.is_uninformative = 'true';
MATCH (t:CogFunctionalCategory) WHERE t.id IN ['cog.category:S'] SET t.is_uninformative = 'true';

MATCH (t:CyanorakRole)
WHERE t.id IN ['cyanorak.role:R','cyanorak.role:R.1','cyanorak.role:R.2',
               'cyanorak.role:R.4','cyanorak.role:R.5']
SET t.is_uninformative = 'true';

MATCH (t:TigrRole)
WHERE t.id IN ['tigr.role:156','tigr.role:704','tigr.role:856',
               'tigr.role:185','tigr.role:157']
SET t.is_uninformative = 'true';

// Pattern-based flag for KEGG (uncharacterized protein KOs, ~210 nodes)
MATCH (t:KeggTerm)
WHERE t.name =~ '^K\\d+;\\s+uncharacterized protein\\b.*'
SET t.is_uninformative = 'true';

// =====================================================================
// F1.2 + F1.3: annotation_quality (numeric 0-3) + annotation_state (enum)
// from informative_source_count over 8 source buckets.
//
// SOURCE_BUCKETS:start
//   live (8): go, kegg, pfam, ec, role, reaction, transporter, cazy
// SOURCE_BUCKETS:end
//
// Maintenance: when adding a new functional Gene-edge type, append a
// has_<bucket> EXISTS line, include in informative_count sum, and add
// the edge type(s) to has_any_edge. See the design spec section
// 'Source bucket maintenance'.
// =====================================================================

MATCH (g:Gene)
CALL {
  WITH g
  WITH g,
       EXISTS { (g)-[:Gene_involved_in_biological_process|Gene_enables_molecular_function|Gene_located_in_cellular_component]->(t)
                WHERE t.is_uninformative IS NULL } AS has_go,
       EXISTS { (g)-[:Gene_has_kegg_ko]->(t) WHERE t.is_uninformative IS NULL } AS has_kegg,
       EXISTS { (g)-[:Gene_has_pfam]->() } AS has_pfam,
       EXISTS { (g)-[:Gene_catalyzes_ec_number]->() } AS has_ec,
       (g.gene_category IS NOT NULL AND g.gene_category <> 'Unknown') AS has_role,
       EXISTS { (g)-[:Gene_catalyzes_reaction]->() } AS has_reaction,
       EXISTS { (g)-[:Gene_has_tcdb_family]->() } AS has_transporter,
       EXISTS { (g)-[:Gene_has_cazy_family]->() } AS has_cazy,
       EXISTS { (g)-[:Gene_involved_in_biological_process|Gene_enables_molecular_function|Gene_located_in_cellular_component
                     |Gene_has_kegg_ko|Gene_has_pfam|Gene_catalyzes_ec_number
                     |Gene_in_cog_category|Gene_has_cyanorak_role|Gene_has_tigr_role
                     |Gene_catalyzes_reaction|Gene_has_tcdb_family|Gene_has_cazy_family]->() } AS has_any_edge
  WITH g,
       (CASE WHEN has_go THEN 1 ELSE 0 END
        + CASE WHEN has_kegg THEN 1 ELSE 0 END
        + CASE WHEN has_pfam THEN 1 ELSE 0 END
        + CASE WHEN has_ec THEN 1 ELSE 0 END
        + CASE WHEN has_role THEN 1 ELSE 0 END
        + CASE WHEN has_reaction THEN 1 ELSE 0 END
        + CASE WHEN has_transporter THEN 1 ELSE 0 END
        + CASE WHEN has_cazy THEN 1 ELSE 0 END) AS informative_count,
       has_any_edge
  SET g.annotation_state =
        CASE
          WHEN informative_count >= 2 THEN 'informative_multi'
          WHEN informative_count = 1 THEN 'informative_single'
          WHEN has_any_edge THEN 'catch_all_only'
          ELSE 'no_evidence'
        END,
      g.annotation_quality =
        CASE
          WHEN informative_count >= 2 THEN 3
          WHEN informative_count = 1 THEN 2
          WHEN has_any_edge THEN 1
          ELSE 0
        END
} IN TRANSACTIONS OF 500 ROWS;

// -----------------------------------------------------------------------
// Gene routing signals (pre-computed for fast gene_overview queries)
// -----------------------------------------------------------------------

// annotation_types: which ontology edge types exist for each gene
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

// =====================================================================
// F1.4: Gene.informative_annotation_types — granular per-source list,
// only includes a source if at least one connected term is informative.
// Parallel of Gene.annotation_types (presence-by-source) with the
// informativeness filter applied.
// =====================================================================

MATCH (g:Gene)
CALL {
  WITH g
  SET g.informative_annotation_types =
    CASE WHEN EXISTS { (g)-[:Gene_involved_in_biological_process]->(t)
                       WHERE t.is_uninformative IS NULL }
         THEN ['go_bp'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_enables_molecular_function]->(t)
                       WHERE t.is_uninformative IS NULL }
         THEN ['go_mf'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_located_in_cellular_component]->(t)
                       WHERE t.is_uninformative IS NULL }
         THEN ['go_cc'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_pfam]->() } THEN ['pfam'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_in_cog_category]->(t)
                       WHERE t.is_uninformative IS NULL }
         THEN ['cog_category'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_kegg_ko]->(t)
                       WHERE t.is_uninformative IS NULL }
         THEN ['kegg'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_kegg_ko]->()-[:Kegg_term_in_brite_category]->() }
         THEN ['brite'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_catalyzes_ec_number]->() } THEN ['ec'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_cyanorak_role]->(t)
                       WHERE t.is_uninformative IS NULL }
         THEN ['cyanorak_role'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_tigr_role]->(t)
                       WHERE t.is_uninformative IS NULL }
         THEN ['tigr_role'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_catalyzes_reaction]->() } THEN ['reaction'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_tcdb_family]->() } THEN ['transporter'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_cazy_family]->() } THEN ['cazy'] ELSE [] END
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

// rank_by_effect: within each experiment + timepoint, rank by |log2FC| descending
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

// rank_up: among significant_up edges per experiment + timepoint, rank by |log2FC| descending
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

// rank_down: among significant_down edges per experiment + timepoint, rank by |log2FC| descending
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

// Gene: cluster_membership_count + cluster_types
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

// metabolite_count: distinct metabolites reachable via Tcdb_family_transports_metabolite.
// Single-hop: substrate edges are rolled up to every ancestor in the adapter
// (each TcdbFamily has direct edges to all metabolites in its subtree).
MATCH (t:TcdbFamily)
CALL {
  WITH t
  OPTIONAL MATCH (t)-[:Tcdb_family_transports_metabolite]->(m:Metabolite)
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
// Transport arm: Gene -> TcdbFamily -> Metabolite (single-hop; substrate edges
// are rolled up to every ancestor in the adapter, so no descendants walk).
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[:Gene_catalyzes_reaction]->(:Reaction)-[:Reaction_has_metabolite]->(m1:Metabolite)
  WITH g, collect(DISTINCT m1) AS m_cat
  OPTIONAL MATCH (g)-[:Gene_has_tcdb_family]->(:TcdbFamily)-[:Tcdb_family_transports_metabolite]->(m2:Metabolite)
  WITH g, m_cat, collect(DISTINCT m2) AS m_tr
  SET g.metabolite_count = size(apoc.coll.toSet(m_cat + m_tr))
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

// Metabolite.gene_count: UNION across catalysis + transport paths (single-hop).
// Substrate edges are rolled up to every ancestor in the adapter, so the
// transport arm is a 1-hop traversal at any TcdbFamily level the gene is
// annotated at. organism_count is computed below from the materialized
// Organism_has_metabolite edge so size(organism_names) == organism_count
// is invariant by construction (KG-A8).
CALL {
  MATCH (m:Metabolite)
  OPTIONAL MATCH (m)<-[:Reaction_has_metabolite]-(:Reaction)<-[:Gene_catalyzes_reaction]-(g_cat:Gene)
  WITH m, collect(DISTINCT g_cat) AS gs_cat
  OPTIONAL MATCH (m)<-[:Tcdb_family_transports_metabolite]-(:TcdbFamily)<-[:Gene_has_tcdb_family]-(g_tr:Gene)
  WITH m, gs_cat, collect(DISTINCT g_tr) AS gs_tr
  WITH m, apoc.coll.toSet(gs_cat + gs_tr) AS all_g
  SET m.gene_count = size(all_g)
} IN TRANSACTIONS OF 1000 ROWS;

// Metabolite.transporter_count: distinct tc_specificity LEAVES with substrate edge.
// Filter to leaves so the count reflects "actual transporter systems" rather
// than ancestors-via-rollup. Source-level filter recovers pre-rollup semantics.
CALL {
  MATCH (m:Metabolite)
  OPTIONAL MATCH (t:TcdbFamily {level_kind: 'tc_specificity'})-[:Tcdb_family_transports_metabolite]->(m)
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

// Materialize Organism_has_metabolite (transport arm) — single-hop after rollup.
CALL {
  MATCH (o:OrganismTaxon)<-[:Gene_belongs_to_organism]-(g:Gene)
        -[:Gene_has_tcdb_family]->(:TcdbFamily)
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

// ── Chemistry slice-1 follow-up rollups (KG-A5, A6, A7, A8) ────────────────
// Denormalize Metabolite_in_pathway and Organism_has_metabolite edges onto
// Metabolite for list_metabolites flat filter + per-row use.

// KG-A5/A6/A7: pathway_ids / pathway_names / pathway_count.
// collect(DISTINCT p) keeps id/name index-aligned (single ordered node list);
// ORDER BY p.id sorts pathways alphabetically before collect.
MATCH (m:Metabolite)
CALL {
  WITH m
  OPTIONAL MATCH (m)-[:Metabolite_in_pathway]->(p:KeggTerm)
  WITH m, p ORDER BY p.id
  WITH m, collect(DISTINCT p) AS ps
  SET m.pathway_ids   = [x IN ps | x.id],
      m.pathway_names = [x IN ps | x.name],
      m.pathway_count = size(ps)
} IN TRANSACTIONS OF 1000 ROWS;

// KG-A8: organism_names — distinct sorted OrganismTaxon.preferred_name reachable
// via Organism_has_metabolite (UNION of catalysis + transport, materialized above).
// organism_count is recomputed from the same edge so size(organism_names) ==
// organism_count is invariant by construction. The earlier transport-arm rollup
// missed genes annotated above the tc_specificity leaf; deriving from the
// materialized edge picks up the descendants walk done in the materialization.
MATCH (m:Metabolite)
CALL {
  WITH m
  OPTIONAL MATCH (org:OrganismTaxon)-[:Organism_has_metabolite]->(m)
  WITH m, collect(DISTINCT org) AS orgs
  SET m.organism_names = apoc.coll.sort([o IN orgs | o.preferred_name]),
      m.organism_count = size(orgs)
} IN TRANSACTIONS OF 1000 ROWS;
