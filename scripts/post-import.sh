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
// Composite RANGE index backing the gene_neighbors genomic-window query:
// equality prefix (organism_name, contig) + ordered range suffix start.
CREATE INDEX gene_org_contig_start_idx IF NOT EXISTS FOR (g:Gene) ON (g.organism_name, g.contig, g.start);
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

// PSORTb SubcellularLocalization (flat ontology + scored edge)
CREATE INDEX subcellular_localization_level_idx IF NOT EXISTS FOR (n:SubcellularLocalization) ON (n.level);
CREATE INDEX subcellular_localization_id_idx IF NOT EXISTS FOR (n:SubcellularLocalization) ON (n.psortb_id);
CREATE FULLTEXT INDEX subcellularLocalizationFullText IF NOT EXISTS
    FOR (n:SubcellularLocalization) ON EACH [n.name, n.psortb_id];

// SignalP SignalPeptideType (flat ontology + scored edge)
CREATE INDEX signal_peptide_type_level_idx IF NOT EXISTS FOR (n:SignalPeptideType) ON (n.level);
CREATE INDEX signal_peptide_type_id_idx IF NOT EXISTS FOR (n:SignalPeptideType) ON (n.signalp_id);
CREATE FULLTEXT INDEX signalPeptideTypeFullText IF NOT EXISTS
    FOR (n:SignalPeptideType) ON EACH [n.name, n.signalp_id];

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
CYPHER

# ─────────────────────────────────────────────────────────────────────────────
# Group 2: small-table aggregations. Ordering matters:
#   growth_phases (on Experiment) feeds Publication rollup;
#   expression_status on edges feeds Experiment summary and rank_up/rank_down.
# ─────────────────────────────────────────────────────────────────────────────
echo "=== Post-process: Compute small-table aggregations ==="
time cypher-shell <<'CYPHER'
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
} IN TRANSACTIONS OF 1000 ROWS;

// annotation_types (includes tcdb/cazy — folded in from a former extension pass)
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
    CASE WHEN EXISTS { (g)-[:Gene_has_tigr_role]->() } THEN ['tigr_role'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_tcdb_family]->() } THEN ['tcdb'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_cazy_family]->() } THEN ['cazy'] ELSE [] END
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
} IN TRANSACTIONS OF 1000 ROWS;

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
} IN TRANSACTIONS OF 30 ROWS;

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
} IN TRANSACTIONS OF 30 ROWS;

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
} IN TRANSACTIONS OF 30 ROWS;

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
} IN TRANSACTIONS OF 30 ROWS;

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

// Gene DM routing rollups (combined single Gene scan).
// Each OPTIONAL MATCH is followed by a WITH aggregation that collapses rows
// back to one row per gene before the next OPTIONAL MATCH, so row count never
// multiplies. count(DISTINCT) and collect(DISTINCT) handle the no-edge case
// (count=0, types=[]); no separate defaults pass needed.
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (dm_n:DerivedMetric)-[:Derived_metric_quantifies_gene]->(g)
  WITH g,
       count(DISTINCT dm_n) AS n_cnt,
       [x IN collect(DISTINCT dm_n.metric_type) WHERE x IS NOT NULL] AS n_types
  OPTIONAL MATCH (dm_b:DerivedMetric)-[:Derived_metric_flags_gene]->(g)
  WITH g, n_cnt, n_types,
       count(DISTINCT dm_b) AS b_cnt,
       [x IN collect(DISTINCT dm_b.metric_type) WHERE x IS NOT NULL] AS b_types
  OPTIONAL MATCH (dm_c:DerivedMetric)-[:Derived_metric_classifies_gene]->(g)
  WITH g, n_cnt, n_types, b_cnt, b_types,
       count(DISTINCT dm_c) AS c_cnt,
       [x IN collect(DISTINCT dm_c.metric_type) WHERE x IS NOT NULL] AS c_types
  OPTIONAL MATCH (dm_a:DerivedMetric)
    -[:Derived_metric_quantifies_gene|Derived_metric_flags_gene|Derived_metric_classifies_gene]->(g)
  WITH g, n_cnt, n_types, b_cnt, b_types, c_cnt, c_types,
       [x IN collect(DISTINCT dm_a.compartment) WHERE x IS NOT NULL] AS comps
  SET g.numeric_metric_count = n_cnt,
      g.numeric_metric_types_observed = apoc.coll.sort(n_types),
      g.boolean_metric_count = b_cnt,
      g.boolean_metric_types_observed = apoc.coll.sort(b_types),
      g.categorical_metric_count = c_cnt,
      g.categorical_metric_types_observed = apoc.coll.sort(c_types),
      g.compartments_observed = apoc.coll.sort(comps)
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

// BriteCategory computed properties: member_ko_count, gene_count, organism_count.
// Single-walk variant — the subtree is traversed once and all 3 aggregates
// (KO leaves, distinct genes, distinct organism_names) are derived from the
// joint OPTIONAL MATCH product. count(DISTINCT) handles the row multiplication
// from gene-side fan-out per KO.
MATCH (b:BriteCategory)
CALL {
  WITH b
  OPTIONAL MATCH (b)<-[:Brite_category_is_a_brite_category*0..]-(:BriteCategory)
                   <-[:Kegg_term_in_brite_category]-(ko:KeggTerm {level_kind: 'ko'})
  OPTIONAL MATCH (ko)<-[:Gene_has_kegg_ko]-(g:Gene)
  WITH b,
       count(DISTINCT ko) AS ko_count,
       count(DISTINCT g) AS g_count,
       collect(DISTINCT g.organism_name) AS orgs
  SET b.member_ko_count = ko_count,
      b.gene_count = g_count,
      b.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 100 ROWS;

// ── TcdbFamily computed properties ───────────────────────────────────────────
// Combined single scan: tc_class_id, member_count, gene_count, organism_count,
// metabolite_count.

// tc_class_id: walk up at most to nearest tc_class ancestor. *0..* lets a
// tc_class node match itself; non-class nodes walk up (label filter rejects
// the 0-hop self-match for non-class). Done in its own pass because LIMIT 1
// inside a chained-aggregation block doesn't compose cleanly with later
// OPTIONAL MATCHes.
MATCH (t:TcdbFamily)
CALL {
  WITH t
  MATCH (t)-[:Tcdb_family_is_a_tcdb_family*0..]->(cls:TcdbFamily {level_kind: 'tc_class'})
  WITH t, cls LIMIT 1
  SET t.tc_class_id = cls.id
} IN TRANSACTIONS OF 1000 ROWS;

// member_count + gene_count + organism_count + metabolite_count.
// Substrate edges are rolled up to every ancestor in the adapter (each TcdbFamily
// has direct edges to all metabolites in its subtree), so metabolite_count is
// a single-hop count.
MATCH (t:TcdbFamily)
CALL {
  WITH t
  OPTIONAL MATCH (child:TcdbFamily)-[:Tcdb_family_is_a_tcdb_family]->(t)
  WITH t, count(child) AS mc
  OPTIONAL MATCH (t)<-[:Tcdb_family_is_a_tcdb_family*0..]-(:TcdbFamily)<-[:Gene_has_tcdb_family]-(g:Gene)
  WITH t, mc, count(DISTINCT g) AS gc, collect(DISTINCT g.organism_name) AS orgs
  OPTIONAL MATCH (t)-[:Tcdb_family_transports_metabolite]->(m:Metabolite)
  WITH t, mc, gc, orgs, count(DISTINCT m) AS metc
  SET t.member_count = mc,
      t.gene_count = gc,
      t.organism_count = size([x IN orgs WHERE x IS NOT NULL]),
      t.metabolite_count = metc
} IN TRANSACTIONS OF 1000 ROWS;

// is_promiscuous: families with broad substrate or member coverage. Threshold
// (metabolite_count >= 50 OR member_count >= 100) flags ~30 of 12,883 families
// — clearly in the long tail (p95(metabolite_count) on tc_family ≈ 14).
// Consumed by explorer family_inferred-dominance warnings to distinguish
// curation-effort gaps from biologically-promiscuous transporters (KG-MET-006).
MATCH (t:TcdbFamily)
SET t.is_promiscuous =
  (coalesce(t.metabolite_count, 0) >= 50) OR
  (coalesce(t.member_count, 0) >= 100);

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

// ── SubcellularLocalization computed properties (PSORTb; flat ontology) ───────
// gene_count + organism_count: direct gene->node traversal (no *0.. — flat).
MATCH (n:SubcellularLocalization)
CALL {
  WITH n
  OPTIONAL MATCH (n)<-[:Gene_has_subcellular_localization]-(g:Gene)
  WITH n, count(DISTINCT g) AS gc, collect(DISTINCT g.organism_name) AS orgs
  SET n.gene_count = gc,
      n.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 1000 ROWS;

// Gene.subcellular_localization: denormalized 1:1 routing string (the gene's
// single PSORTb call, absent when no confident localization). STRUCTURAL — so
// deliberately NOT folded into annotation_types / annotation_quality.
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[:Gene_has_subcellular_localization]->(loc:SubcellularLocalization)
  WITH g, loc.psortb_id AS lid
  SET g.subcellular_localization = lid
} IN TRANSACTIONS OF 1000 ROWS;

// rank_by_score on Gene_has_subcellular_localization: within each localization,
// rank genes by descending PSORTb score (1 = strongest). Mirrors rank_by_effect.
MATCH (n:SubcellularLocalization)
CALL {
  WITH n
  MATCH (g:Gene)-[r:Gene_has_subcellular_localization]->(n)
  WITH r, r.score AS s, g.locus_tag AS lt
  ORDER BY s DESC, lt ASC
  WITH collect(r) AS edges
  UNWIND range(0, size(edges) - 1) AS i
  SET (edges[i]).rank_by_score = i + 1
} IN TRANSACTIONS OF 1000 ROWS;

// ── SignalPeptideType computed properties (SignalP; flat ontology) ────────────
// gene_count + organism_count: direct gene->node traversal (no *0.. — flat).
MATCH (n:SignalPeptideType)
CALL {
  WITH n
  OPTIONAL MATCH (n)<-[:Gene_has_signal_peptide_type]-(g:Gene)
  WITH n, count(DISTINCT g) AS gc, collect(DISTINCT g.organism_name) AS orgs
  SET n.gene_count = gc,
      n.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 1000 ROWS;

// Gene.signal_peptide_type: denormalized 1:1 routing string (the gene's single
// SignalP call, absent when no confident signal peptide). STRUCTURAL — so
// deliberately NOT folded into annotation_types / annotation_quality.
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[:Gene_has_signal_peptide_type]->(spt:SignalPeptideType)
  WITH g, spt.signalp_id AS sid
  SET g.signal_peptide_type = sid
} IN TRANSACTIONS OF 1000 ROWS;

// rank_by_probability on Gene_has_signal_peptide_type: within each type, rank
// genes by descending SignalP probability (1 = strongest). Mirrors rank_by_score.
MATCH (n:SignalPeptideType)
CALL {
  WITH n
  MATCH (g:Gene)-[r:Gene_has_signal_peptide_type]->(n)
  WITH r, r.probability AS p, g.locus_tag AS lt
  ORDER BY p DESC, lt ASC
  WITH collect(r) AS edges
  UNWIND range(0, size(edges) - 1) AS i
  SET (edges[i]).rank_by_probability = i + 1
} IN TRANSACTIONS OF 1000 ROWS;

// ── Gene routing extensions ──────────────────────────────────────────────────
// (Note: 'tcdb' / 'cazy' membership in annotation_types is folded into the base
// annotation_types statement above — no separate extension pass.)

// Gene metabolism + ontology counts (combined single Gene scan):
//   - tcdb_family_count  (TCDB-S1)
//   - cazy_family_count  (TCDB-S2)
//   - reaction_count     (KG-A1)
//   - metabolite_count   (TCDB-S3 / KG-A2: UNION of catalysis + transport paths)
//
// Each OPTIONAL MATCH is followed by a WITH aggregation so rows don't multiply.
// count(DISTINCT ...) is required for the count rollups because the chained
// metabolite OPTIONAL MATCHes would otherwise duplicate the parent edge per metabolite.
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[r1:Gene_has_tcdb_family]->(tcd:TcdbFamily)
  OPTIONAL MATCH (tcd)-[:Tcdb_family_transports_metabolite]->(m_tr:Metabolite)
  WITH g, count(DISTINCT r1) AS tc_count, collect(DISTINCT m_tr) AS m_transport
  OPTIONAL MATCH (g)-[r2:Gene_has_cazy_family]->()
  WITH g, tc_count, m_transport, count(r2) AS cz_count
  OPTIONAL MATCH (g)-[r3:Gene_catalyzes_reaction]->(rx:Reaction)
  OPTIONAL MATCH (rx)-[:Reaction_has_metabolite]->(m_cat:Metabolite)
  WITH g, tc_count, m_transport, cz_count,
       count(DISTINCT r3) AS rxn_count,
       collect(DISTINCT m_cat) AS m_catalysis
  SET g.tcdb_family_count = tc_count,
      g.cazy_family_count = cz_count,
      g.reaction_count = rxn_count,
      g.metabolite_count = size(apoc.coll.toSet(m_catalysis + m_transport))
} IN TRANSACTIONS OF 1000 ROWS;

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

// Materialize Organism_has_metabolite (catalysis arm) — also tags evidence_sources
// inline at MERGE time (cheap; avoids a per-edge re-traversal pass after).
// measured_* defaults set at ON CREATE; the measurement-arm MERGE further down
// overrides them with real values when an assay also reaches this (o, m) pair.
CALL {
  MATCH (o:OrganismTaxon)<-[:Gene_belongs_to_organism]-(g:Gene)
        -[:Gene_catalyzes_reaction]->(:Reaction)
        -[:Reaction_has_metabolite]->(m:Metabolite)
  WITH DISTINCT o, m
  MERGE (o)-[r:Organism_has_metabolite]->(m)
  ON CREATE SET r.evidence_sources = ['metabolism'],
                r.measured_assay_count = 0,
                r.measured_compartments = [],
                r.measured_paper_count = 0
  ON MATCH  SET r.evidence_sources =
    CASE WHEN 'metabolism' IN coalesce(r.evidence_sources, [])
         THEN r.evidence_sources
         ELSE coalesce(r.evidence_sources, []) + 'metabolism' END
} IN TRANSACTIONS OF 1000 ROWS;

// Materialize Organism_has_metabolite (transport arm) — single-hop after rollup.
CALL {
  MATCH (o:OrganismTaxon)<-[:Gene_belongs_to_organism]-(g:Gene)
        -[:Gene_has_tcdb_family]->(:TcdbFamily)
        -[:Tcdb_family_transports_metabolite]->(m:Metabolite)
  WITH DISTINCT o, m
  MERGE (o)-[r:Organism_has_metabolite]->(m)
  ON CREATE SET r.evidence_sources = ['transport'],
                r.measured_assay_count = 0,
                r.measured_compartments = [],
                r.measured_paper_count = 0
  ON MATCH  SET r.evidence_sources =
    CASE WHEN 'transport' IN coalesce(r.evidence_sources, [])
         THEN r.evidence_sources
         ELSE coalesce(r.evidence_sources, []) + 'transport' END
} IN TRANSACTIONS OF 1000 ROWS;

// Organism rollup props (~30 organisms — single batch fits comfortably)
CALL {
  MATCH (o:OrganismTaxon)-[:Organism_has_metabolite]->(m:Metabolite)
  WITH o, count(DISTINCT m) AS metabolite_count
  SET o.metabolite_count = metabolite_count
} IN TRANSACTIONS OF 1000 ROWS;

CALL {
  MATCH (o:OrganismTaxon)<-[:Gene_belongs_to_organism]-(:Gene)
        -[:Gene_catalyzes_reaction]->(r:Reaction)
  WITH o, count(DISTINCT r) AS reaction_count
  SET o.reaction_count = reaction_count
} IN TRANSACTIONS OF 1000 ROWS;

// ── Chemistry slice-1 rollups (KG-A4) ───────────────────────
// (Note: KG-A1 Gene.reaction_count and KG-A2 Gene.metabolite_count are folded
// into the combined Gene metabolism + ontology rollup statement above.)

// KG-A4: KeggTerm pathway-level rollups (sparse on pathways only).
// level_kind = 'pathway' filter; KOs / categories left unset.
MATCH (p:KeggTerm) WHERE p.level_kind = 'pathway'
CALL {
  WITH p
  OPTIONAL MATCH (r:Reaction)-[:Reaction_in_kegg_pathway]->(p)
  WITH p, count(r) AS rxn_count
  SET p.reaction_count = rxn_count
} IN TRANSACTIONS OF 1000 ROWS;

MATCH (p:KeggTerm) WHERE p.level_kind = 'pathway'
CALL {
  WITH p
  OPTIONAL MATCH (m:Metabolite)-[:Metabolite_in_pathway]->(p)
  WITH p, count(m) AS met_count
  SET p.metabolite_count = met_count
} IN TRANSACTIONS OF 1000 ROWS;

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

// ── Phase 2 metabolomics: MetaboliteAssay indexes + rollups ───────────────────

CREATE INDEX metabolite_assay_organism_idx     IF NOT EXISTS FOR (a:MetaboliteAssay) ON (a.organism_name);
CREATE INDEX metabolite_assay_compartment_idx  IF NOT EXISTS FOR (a:MetaboliteAssay) ON (a.compartment);
CREATE INDEX metabolite_assay_metric_type_idx  IF NOT EXISTS FOR (a:MetaboliteAssay) ON (a.metric_type);
CREATE INDEX metabolite_assay_value_kind_idx   IF NOT EXISTS FOR (a:MetaboliteAssay) ON (a.value_kind);
CREATE INDEX metabolite_assay_experiment_idx   IF NOT EXISTS FOR (a:MetaboliteAssay) ON (a.experiment_id);
CREATE FULLTEXT INDEX metaboliteAssayFullText  IF NOT EXISTS
  FOR (a:MetaboliteAssay) ON EACH [a.name, a.field_description, a.treatment, a.experimental_context];

// MetaboliteAssay numeric ranks: per-assay, only when rankable='true'.
// Mirrors DerivedMetric pattern (per-assay scope, deterministic Metabolite.id tiebreaker,
// pinned bucket thresholds 90/75/25).
MATCH (a:MetaboliteAssay {rankable: 'true'})
CALL {
  WITH a
  MATCH (a)-[r:Assay_quantifies_metabolite]->(m:Metabolite)
  WITH r, r.value AS val, m.id AS mid
  ORDER BY val DESC, mid ASC
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
} IN TRANSACTIONS OF 30 ROWS;

// MetaboliteAssay total_metabolite_count
MATCH (a:MetaboliteAssay)
CALL {
  WITH a
  OPTIONAL MATCH (a)-[r:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
  WITH a, count(DISTINCT m) AS cnt
  SET a.total_metabolite_count = cnt
} IN TRANSACTIONS OF 1000 ROWS;

// Numeric distribution stats (null for boolean assays)
MATCH (a:MetaboliteAssay {value_kind: 'numeric'})
CALL {
  WITH a
  MATCH (a)-[r:Assay_quantifies_metabolite]->()
  WITH a,
       min(r.value) AS vmin, max(r.value) AS vmax,
       percentileDisc(r.value, 0.25) AS q1,
       percentileDisc(r.value, 0.5)  AS med,
       percentileDisc(r.value, 0.75) AS q3
  SET a.value_min = vmin, a.value_max = vmax,
      a.value_q1 = q1, a.value_median = med, a.value_q3 = q3
} IN TRANSACTIONS OF 1000 ROWS;

// Boolean flag counts (null for numeric assays)
MATCH (a:MetaboliteAssay {value_kind: 'boolean'})
CALL {
  WITH a
  MATCH (a)-[r:Assay_flags_metabolite]->()
  WITH a,
       sum(CASE WHEN r.flag_value='true'  THEN 1 ELSE 0 END) AS t,
       sum(CASE WHEN r.flag_value='false' THEN 1 ELSE 0 END) AS f
  SET a.flag_true_count = t, a.flag_false_count = f
} IN TRANSACTIONS OF 1000 ROWS;

// growth_phases from parent Experiment (mirrors DerivedMetric.growth_phases)
MATCH (a:MetaboliteAssay)
CALL {
  WITH a
  OPTIONAL MATCH (a)<-[:ExperimentHasMetaboliteAssay]-(e:Experiment)
  WITH a, coalesce(e.growth_phases, []) AS phases
  SET a.growth_phases = phases
} IN TRANSACTIONS OF 1000 ROWS;

// ── Metabolite measured_* properties (Phase 2) ────────────────────────────────

MATCH (m:Metabolite)
CALL {
  WITH m
  OPTIONAL MATCH (m)<-[:Assay_quantifies_metabolite|Assay_flags_metabolite]-(a:MetaboliteAssay)
  OPTIONAL MATCH (a)-[:MetaboliteAssayBelongsToOrganism]->(o:OrganismTaxon)
  OPTIONAL MATCH (a)<-[:PublicationHasMetaboliteAssay]-(p:Publication)
  WITH m,
       count(DISTINCT a) AS acnt,
       collect(DISTINCT o.preferred_name) AS orgs,
       collect(DISTINCT a.compartment) AS comps,
       count(DISTINCT p) AS pcnt
  SET m.measured_assay_count = acnt,
      m.measured_organisms = apoc.coll.sort([x IN orgs WHERE x IS NOT NULL]),
      m.measured_compartments = apoc.coll.sort([c IN comps WHERE c IS NOT NULL]),
      m.measured_paper_count = pcnt
} IN TRANSACTIONS OF 1000 ROWS;

// Defaults for Metabolite nodes with no assays
MATCH (m:Metabolite) WHERE m.measured_assay_count IS NULL
SET m.measured_assay_count = 0,
    m.measured_organisms = [],
    m.measured_compartments = [],
    m.measured_paper_count = 0;

// ── Organism_has_metabolite measurement-arm materialization (Phase 2) ─────────
// Adds (organism, metabolite) edges where only the measurement path exists
// (no gene-side catalysis or transport) AND tags 'measured' on edges that
// already exist via catalysis/transport. evidence_sources is set inline at
// MERGE time (mirrors the catalysis + transport blocks above) to avoid a
// per-edge re-traversal pass that scales O(edges × organism-genes).
// measured_assay_count / measured_compartments / measured_paper_count are
// also computed inline here (folded in from a former full-edge augmentation
// pass that scanned every Organism_has_metabolite edge).
CALL {
  MATCH (a:MetaboliteAssay)-[:MetaboliteAssayBelongsToOrganism]->(o:OrganismTaxon)
  MATCH (a)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
  OPTIONAL MATCH (a)<-[:PublicationHasMetaboliteAssay]-(p:Publication)
  WITH o, m,
       count(DISTINCT a) AS acnt,
       collect(DISTINCT a.compartment) AS comps,
       count(DISTINCT p) AS pcnt
  MERGE (o)-[r:Organism_has_metabolite]->(m)
  ON CREATE SET r.evidence_sources = ['measured']
  ON MATCH  SET r.evidence_sources =
    CASE WHEN 'measured' IN coalesce(r.evidence_sources, [])
         THEN r.evidence_sources
         ELSE coalesce(r.evidence_sources, []) + 'measured' END
  SET r.measured_assay_count = acnt,
      r.measured_compartments = [c IN comps WHERE c IS NOT NULL],
      r.measured_paper_count = pcnt
} IN TRANSACTIONS OF 1000 ROWS;

// Recompute Metabolite.organism_count + organism_names AFTER measurement-arm
// materialization (so measured-only pairs contribute to the rollup).
MATCH (m:Metabolite)
CALL {
  WITH m
  OPTIONAL MATCH (org:OrganismTaxon)-[:Organism_has_metabolite]->(m)
  WITH m, collect(DISTINCT org) AS orgs
  SET m.organism_names = apoc.coll.sort([o IN orgs | o.preferred_name]),
      m.organism_count = size(orgs)
} IN TRANSACTIONS OF 1000 ROWS;

// (Note: Organism_has_metabolite measured_* properties are written inline by
// the catalysis-arm + transport-arm + measurement-arm MERGE blocks above —
// no separate full-edge augmentation pass.)

// ── Experiment / Publication / OrganismTaxon rollups (Phase 2) ────────────────

MATCH (e:Experiment)
CALL {
  WITH e
  OPTIONAL MATCH (e)-[:ExperimentHasMetaboliteAssay]->(a:MetaboliteAssay)
  OPTIONAL MATCH (a)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
  WITH e,
       count(DISTINCT a) AS acnt,
       collect(DISTINCT a.compartment) AS comps,
       count(DISTINCT m) AS mcnt
  SET e.metabolite_assay_count = acnt,
      e.metabolite_compartments = [c IN comps WHERE c IS NOT NULL],
      e.metabolite_count = mcnt
} IN TRANSACTIONS OF 1000 ROWS;

MATCH (p:Publication)
CALL {
  WITH p
  OPTIONAL MATCH (p)-[:PublicationHasMetaboliteAssay]->(a:MetaboliteAssay)
  OPTIONAL MATCH (a)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
  WITH p,
       count(DISTINCT a) AS acnt,
       collect(DISTINCT a.compartment) AS comps,
       count(DISTINCT m) AS mcnt
  SET p.metabolite_assay_count = acnt,
      p.metabolite_compartments = [c IN comps WHERE c IS NOT NULL],
      p.metabolite_count = mcnt
} IN TRANSACTIONS OF 1000 ROWS;

MATCH (o:OrganismTaxon)
CALL {
  WITH o
  OPTIONAL MATCH (o)<-[:MetaboliteAssayBelongsToOrganism]-(a:MetaboliteAssay)
  OPTIONAL MATCH (a)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
  WITH o, count(DISTINCT m) AS mcnt
  SET o.measured_metabolite_count = mcnt
} IN TRANSACTIONS OF 1000 ROWS;
CYPHER

# ─────────────────────────────────────────────────────────────────────────────
# Group 4: Schema_info release metadata. Separate cypher-shell invocation because
# its params (version + git identity) are interpolated from the environment — the
# quoted heredocs above deliberately do NOT interpolate. Runs on EVERY build:
# dev leaves KG_* unset → '0.0.0-dev' / 'unknown'; /release-kg sets them.
# Defaults are pushed via bash ${VAR:-default} (the `:-` form fires on unset OR
# empty), so an empty env var still yields the default — coalesce() in the Cypher
# is only a secondary guard (it does not catch empty strings).
# Counts are computed (not hardcoded) so they track data drift.
# Keep the MATCH/SET logic byte-identical to the matching block in post-import.cypher.
# ─────────────────────────────────────────────────────────────────────────────
echo "=== Post-process: Stamp Schema_info release metadata ==="
time cypher-shell \
  -P "version          => '${KG_RELEASE_VERSION:-0.0.0-dev}'" \
  -P "git_sha          => '${KG_GIT_SHA:-unknown}'" \
  -P "git_sha_short    => '${KG_GIT_SHA_SHORT:-unknown}'" \
  -P "git_branch       => '${KG_GIT_BRANCH:-unknown}'" \
  -P "git_dirty        => '${KG_GIT_DIRTY:-unknown}'" \
  -P "mcp_min_version  => '${KG_MCP_MIN_VERSION:-0.1.0}'" \
  -P "release_notes_url => '${KG_RELEASE_NOTES_URL:-}'" \
  <<'CYPHER'
MATCH (s:Schema_info {id: 'schema_info'})
SET s.version           = coalesce($version, '0.0.0-dev'),
    s.built_at          = toString(datetime()),
    s.git_sha           = coalesce($git_sha, 'unknown'),
    s.git_sha_short     = coalesce($git_sha_short, 'unknown'),
    s.git_branch        = coalesce($git_branch, 'unknown'),
    s.git_dirty         = coalesce($git_dirty, 'unknown'),
    s.mcp_min_version   = coalesce($mcp_min_version, '0.1.0'),
    s.release_notes_url = coalesce($release_notes_url, '')
WITH s
SET s.paper_count           = COUNT { (:Publication) },
    s.experiment_count      = COUNT { (:Experiment) },
    s.gene_count            = COUNT { (:Gene) },
    s.organism_count        = COUNT { (:OrganismTaxon) },
    s.expression_edge_count = COUNT { ()-[:Changes_expression_of]->() };
CYPHER

echo "=== Post-process complete ==="
neo4j stop
