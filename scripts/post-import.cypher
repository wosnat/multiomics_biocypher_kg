// Post-import Cypher commands
// These run after the knowledge graph is imported into Neo4j
// Indexes only — homolog edges and expression propagation have been replaced
// by OrthologGroup nodes and Gene_in_ortholog_group edges (query-time joins).

// Scalar indexes for get_gene exact lookup
CREATE INDEX gene_locus_tag_idx IF NOT EXISTS FOR (g:Gene) ON (g.locus_tag);
CREATE INDEX gene_name_idx IF NOT EXISTS FOR (g:Gene) ON (g.gene_name);
CREATE INDEX gene_organism_name_idx IF NOT EXISTS FOR (g:Gene) ON (g.organism_name);
// Composite RANGE index backing the gene_neighbors genomic-window query:
// equality prefix (organism_name, contig) + ordered range suffix start.
CREATE INDEX gene_org_contig_start_idx IF NOT EXISTS FOR (g:Gene) ON (g.organism_name, g.contig, g.start);

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
CREATE INDEX tigr_role_level_idx IF NOT EXISTS FOR (t:TigrRole) ON (t.level);
CREATE INDEX tigr_role_level_kind_idx IF NOT EXISTS FOR (t:TigrRole) ON (t.level_kind);
CREATE INDEX cazy_family_cazy_id_idx IF NOT EXISTS FOR (c:CazyFamily) ON (c.cazy_id);

// TCDB / CAZy full-text indexes
// Full-text defs can't be ALTERed — drop + recreate so description is picked up
// even on reruns against an existing graph (T6 node names).
DROP INDEX tcdbFamilyFullText IF EXISTS;
CREATE FULLTEXT INDEX tcdbFamilyFullText IF NOT EXISTS
    FOR (t:TcdbFamily) ON EACH [t.name, t.tcdb_id, t.superfamily, t.description];
CREATE FULLTEXT INDEX cazyFamilyFullText IF NOT EXISTS
    FOR (c:CazyFamily) ON EACH [c.name, c.cazy_id];

// InterPro entry (hierarchical ontology + scored edge). interpro_type is the
// PRIMARY ORA stratification key (breadth); level is secondary (is-a depth).
CREATE INDEX interpro_entry_level_idx IF NOT EXISTS FOR (e:InterproEntry) ON (e.level);
CREATE INDEX interpro_entry_type_idx IF NOT EXISTS FOR (e:InterproEntry) ON (e.interpro_type);
CREATE INDEX interpro_entry_id_idx IF NOT EXISTS FOR (e:InterproEntry) ON (e.interpro_id);
// Relationship-property indexes (2026-08-28, HO-003): the explorer's trust filters
// (`r.evidence IN $evidence`, `r.evidence_score >= $min`) run on every gene→ontology
// edge type; only the >100K-edge types get an index. `sources` is a list — the
// `any(...)` predicate cannot use a range index, so it is deliberately not indexed.
// Gene_has_interpro_entry carries no evidence_score (constant-source edge).
CREATE INDEX gene_go_bp_evidence_idx        IF NOT EXISTS FOR ()-[r:Gene_involved_in_biological_process]-() ON (r.evidence);
CREATE INDEX gene_go_bp_evidence_score_idx  IF NOT EXISTS FOR ()-[r:Gene_involved_in_biological_process]-() ON (r.evidence_score);
CREATE INDEX gene_go_mf_evidence_idx        IF NOT EXISTS FOR ()-[r:Gene_enables_molecular_function]-() ON (r.evidence);
CREATE INDEX gene_go_mf_evidence_score_idx  IF NOT EXISTS FOR ()-[r:Gene_enables_molecular_function]-() ON (r.evidence_score);
CREATE INDEX gene_go_cc_evidence_idx        IF NOT EXISTS FOR ()-[r:Gene_located_in_cellular_component]-() ON (r.evidence);
CREATE INDEX gene_go_cc_evidence_score_idx  IF NOT EXISTS FOR ()-[r:Gene_located_in_cellular_component]-() ON (r.evidence_score);
CREATE INDEX gene_pfam_evidence_idx         IF NOT EXISTS FOR ()-[r:Gene_has_pfam]-() ON (r.evidence);
CREATE INDEX gene_pfam_evidence_score_idx   IF NOT EXISTS FOR ()-[r:Gene_has_pfam]-() ON (r.evidence_score);
CREATE INDEX gene_interpro_evidence_idx     IF NOT EXISTS FOR ()-[r:Gene_has_interpro_entry]-() ON (r.evidence);
// Full-text defs can't be ALTERed — drop + recreate so new fields are picked up
// even on reruns against an existing graph.
DROP INDEX interproEntryFullText IF EXISTS;
CREATE FULLTEXT INDEX interproEntryFullText IF NOT EXISTS
    FOR (e:InterproEntry) ON EACH [e.name, e.description];

// NCBIfam (flat curated-family ontology; family_type is the stratification key)
CREATE INDEX ncbifam_family_id_idx IF NOT EXISTS FOR (n:NcbifamFamily) ON (n.ncbifam_id);
CREATE INDEX ncbifam_family_type_idx IF NOT EXISTS FOR (n:NcbifamFamily) ON (n.family_type);
CREATE INDEX ncbifam_family_level_idx IF NOT EXISTS FOR (n:NcbifamFamily) ON (n.level);
CREATE FULLTEXT INDEX ncbifamFamilyFullText IF NOT EXISTS
    FOR (n:NcbifamFamily) ON EACH [n.name, n.gene_symbol, n.description];

// MEROPS MeropsFamily (hierarchical clan→family→subfamily ontology + scored edge)
CREATE INDEX merops_family_level_idx IF NOT EXISTS FOR (m:MeropsFamily) ON (m.level);
CREATE INDEX merops_family_level_kind_idx IF NOT EXISTS FOR (m:MeropsFamily) ON (m.level_kind);
CREATE INDEX merops_family_id_idx IF NOT EXISTS FOR (m:MeropsFamily) ON (m.merops_id);
// Full-text defs can't be ALTERed — drop + recreate so cleavage_summary is
// picked up even on reruns against an existing graph.
DROP INDEX meropsFamilyFullText IF EXISTS;
CREATE FULLTEXT INDEX meropsFamilyFullText IF NOT EXISTS
    FOR (m:MeropsFamily) ON EACH [m.name, m.merops_id, m.description, m.cleavage_summary];

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
// Organism name search incl. registry synonyms / taxonomy notes (2026-08-28), so a
// query for the current NCBI name ("Meiothermus taiwanensis") finds MruberA even
// though preferred_name keeps the paper's name.
CREATE FULLTEXT INDEX organismTaxonFullText IF NOT EXISTS
  FOR (o:OrganismTaxon) ON EACH [o.preferred_name, o.organism_name, o.strain_name, o.species, o.name_synonyms, o.taxonomy_note];

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

// Dense list properties on experiment-shaped nodes. The adapters emit [] for
// treatment_type / background_factors, but neo4j-admin import materializes an
// empty string[] cell as NO property. Re-materialize [] so every consumer can
// rely on a list: treatment_type = [] means "characterization experiment, no
// perturbation", never missing data (2026-08-27).
MATCH (e:Experiment)
SET e.treatment_type = coalesce(e.treatment_type, []),
    e.background_factors = coalesce(e.background_factors, []);

MATCH (n)
WHERE n:ClusteringAnalysis OR n:DerivedMetric OR n:MetaboliteAssay
SET n.treatment_type = coalesce(n.treatment_type, []),
    n.background_factors = coalesce(n.background_factors, []);

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
     count(CASE WHEN r.value = 'flagged'     THEN 1 END) AS n_true,
     count(CASE WHEN r.value = 'not_flagged' THEN 1 END) AS n_false
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
SET e.reports_fold_change = 'no_fold_change',
    e.reports_derived_metric_types = [],
    e.derived_metric_count = 0,
    e.derived_metric_value_kinds = [],
    e.derived_metric_gene_count = 0;

// Experiment reports_fold_change: 'fold_change' iff outgoing Changes_expression_of exists
// (R5 two-state string: fold_change | no_fold_change)
MATCH (e:Experiment)
WHERE EXISTS { (e)-[:Changes_expression_of]->() }
SET e.reports_fold_change = 'fold_change';

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
WHERE t.id IN ['tigr.role:156','tigr.role:704','tigr.role:856','tigr.role:270',
               'tigr.role:185','tigr.role:157',
               'tigr.role:hypothetical_proteins','tigr.role:unknown_function',
               'tigr.role:unclassified']
SET t.is_uninformative = 'true';

// Pattern-based flag for KEGG (uncharacterized protein KOs, ~210 nodes)
MATCH (t:KeggTerm)
WHERE t.name =~ '^K\\d+;\\s+uncharacterized protein\\b.*'
SET t.is_uninformative = 'true';

// KEGG global / overview maps (the ko011xx-ko013xx block: ko01100 "Metabolic
// pathways" alone rolls up 27K genes). Unions of other pathways, so they carry
// no pathway-level class signal and are the standard KEGG-ORA exclusion. They
// are parentless level-2 nodes in the KG (no "Global and overview maps"
// subcategory node exists to key a structural rule on), hence an id list —
// mirrored in config/uninformative_terms.yaml (DOC-002, 2026-08-29). ko01310
// Nitrogen cycle / ko01320 Sulfur cycle are parentless too but are narrow
// (16 / 22 KOs) class-bearing subsets, not unions — deliberately NOT flagged.
MATCH (t:KeggTerm)
WHERE t.id IN ['kegg.pathway:ko01100', 'kegg.pathway:ko01110', 'kegg.pathway:ko01120', 'kegg.pathway:ko01200',
               'kegg.pathway:ko01210', 'kegg.pathway:ko01212', 'kegg.pathway:ko01220', 'kegg.pathway:ko01230',
               'kegg.pathway:ko01232', 'kegg.pathway:ko01240', 'kegg.pathway:ko01250']
SET t.is_uninformative = 'true';

// InterPro: unknown-function entries (name-pattern rule; uninformative_terms.yaml)
MATCH (t:InterproEntry)
WHERE t.name =~ '^Protein of unknown function.*'
   OR t.name =~ '^Domain of unknown function.*'
   OR t.name =~ '^Uncharacteri[sz]ed protein family.*'
SET t.is_uninformative = 'true';

// NCBIfam: typed rule (family_type — third rule kind) + name-pattern fallback
MATCH (t:NcbifamFamily)
WHERE t.family_type IN ['hypoth_equivalog', 'hypoth_equivalog_domain']
   OR t.name =~ '(?i).*hypothetical.*'
   OR t.name =~ '(?i).*uncharacterized.*'
   OR t.name =~ '.*\\bDUF\\d.*'
SET t.is_uninformative = 'true';

// =====================================================================
// F1.2 + F1.3: annotation_quality (numeric 0-3) + annotation_state (enum)
// from informative_source_count over 9 source buckets.
//
// SOURCE_BUCKETS:start
//   live (9): go, kegg, pfam, ec, role, reaction, transporter, cazy, ncbifam
// TIER GATE on the TCDB buckets below.
// Gene_has_tcdb_family carries two evidence sources. eggNOG's KEGG_TC is ortholog
// transfer; diamond is direct sequence similarity, and its tier is BOTH a
// confidence band and the depth the call was truncated to (1 = tc_specificity,
// identity>=70; 2 = tc_subfamily, >=40; 3 = tc_family, no identity floor).
//
// Tier 3 is conservative-by-design remote homology, not noise (median identity
// 34%, median e-value 2.3e-30), so those edges ARE emitted and are findable. But
// ~22K of them must not silently upgrade annotation_quality, which was calibrated
// against curated sources and drives genes_by_function routing. So the three
// quality/informativeness buckets count a gene as transporter-annotated only when
// eggNOG called it, or diamond called it at tier<=2.
//
// Gene.tcdb_family_count (most-specific attachments) and Gene.catalyzed_metabolite_count are deliberately NOT tier-gated —
// they are routing counts, not quality signals.

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
       EXISTS { MATCH (g)-[rt:Gene_has_tcdb_family]->() WHERE 'eggnog' IN rt.sources OR rt.tier <= 2 } AS has_transporter,
       EXISTS { (g)-[:Gene_has_cazy_family]->() } AS has_cazy,
       EXISTS { (g)-[:Gene_has_ncbifam_family]->(t) WHERE t.is_uninformative IS NULL } AS has_ncbifam,
       EXISTS { (g)-[:Gene_involved_in_biological_process|Gene_enables_molecular_function|Gene_located_in_cellular_component
                     |Gene_has_kegg_ko|Gene_has_pfam|Gene_catalyzes_ec_number
                     |Gene_in_cog_category|Gene_has_cyanorak_role|Gene_has_tigr_role
                     |Gene_catalyzes_reaction|Gene_has_tcdb_family|Gene_has_cazy_family
                     |Gene_has_interpro_entry|Gene_has_ncbifam_family|Gene_has_merops_family]->() } AS has_any_edge
  WITH g,
       (CASE WHEN has_go THEN 1 ELSE 0 END
        + CASE WHEN has_kegg THEN 1 ELSE 0 END
        + CASE WHEN has_pfam THEN 1 ELSE 0 END
        + CASE WHEN has_ec THEN 1 ELSE 0 END
        + CASE WHEN has_role THEN 1 ELSE 0 END
        + CASE WHEN has_reaction THEN 1 ELSE 0 END
        + CASE WHEN has_transporter THEN 1 ELSE 0 END
        + CASE WHEN has_cazy THEN 1 ELSE 0 END
        + CASE WHEN has_ncbifam THEN 1 ELSE 0 END) AS informative_count,
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

// -----------------------------------------------------------------------
// Gene routing signals (pre-computed for fast gene_overview queries)
// -----------------------------------------------------------------------

// annotation_types: which ontology edge types exist for each gene
// (includes tcdb/cazy — folded in from a former extension pass)
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
    CASE WHEN EXISTS { MATCH (g)-[rt:Gene_has_tcdb_family]->() WHERE 'eggnog' IN rt.sources OR rt.tier <= 2 } THEN ['tcdb'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_cazy_family]->() } THEN ['cazy'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_interpro_entry]->() } THEN ['interpro'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_ncbifam_family]->() } THEN ['ncbifam'] ELSE [] END +
    CASE WHEN EXISTS { MATCH (g)-[rm:Gene_has_merops_family]->() WHERE rm.tier <= 2 } THEN ['merops'] ELSE [] END
} IN TRANSACTIONS OF 1000 ROWS;

// NOTE: 'ncbifam' IS folded into both annotation_types (routing) AND the
// informative_annotation_types / annotation_quality 9-bucket count (F1.4 below) —
// NCBIfam is a curated, functionally specific ontology (unlike Pfam-breadth
// InterPro FAMILY/DOMAIN entries), so it earns a quality bucket the same way
// pfam/go/ec do. 'interpro' is deliberately NOT — it stays a conduit/routing
// signal only (highly redundant with pfam/go/ec — would lift few genes' quality
// for large blast radius). See the design spec §4.2 / §5.
//
// NOTE: 'merops' is TIER-GATED (tier <= 2), exactly like the tcdb tier gate:
// tier-3 calls are conservative remote homology (92% of candidates), findable
// through the edges but not allowed to inflate annotation coverage. merops has
// NO annotation_quality bucket (single evidence source, tier-3-dominated —
// see docs/superpowers/specs/2026-08-17-merops-kg-integration-design.md);
// Gene_has_merops_family does count toward has_any_edge (catch_all_only),
// the interpro/tcdb-tier-3 precedent.

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
    CASE WHEN EXISTS { MATCH (g)-[rt:Gene_has_tcdb_family]->() WHERE 'eggnog' IN rt.sources OR rt.tier <= 2 } THEN ['transporter'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_cazy_family]->() } THEN ['cazy'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_ncbifam_family]->(t)
                       WHERE t.is_uninformative IS NULL }
         THEN ['ncbifam'] ELSE [] END +
    CASE WHEN EXISTS { MATCH (g)-[rm:Gene_has_merops_family]->() WHERE rm.tier <= 2 } THEN ['merops'] ELSE [] END
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
} IN TRANSACTIONS OF 30 ROWS;

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
} IN TRANSACTIONS OF 30 ROWS;

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
} IN TRANSACTIONS OF 30 ROWS;

// Numeric DM rank/percentile/bucket: ranks derived_metric_quantifies_gene edges
// grouped by DerivedMetric, only when parent DM has rankable='rankable'. Ties on
// value broken by Gene.locus_tag ascending (reproducibility).
// Percentile: rank 1 (highest value) -> 100.0; rank N (lowest) -> 0.0.
// Buckets pinned per slice spec §Post-import (thresholds must not drift).
MATCH (dm:DerivedMetric {rankable: 'rankable'})
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
// only when parent DM has has_p_value='p_value' AND p_value_threshold IS NOT NULL
// AND the edge's adjusted_p_value is non-null. Left null otherwise.
MATCH (dm:DerivedMetric {has_p_value: 'p_value'})
WHERE dm.p_value_threshold IS NOT NULL
CALL {
  WITH dm
  MATCH (dm)-[r:Derived_metric_quantifies_gene]->()
  WHERE r.adjusted_p_value IS NOT NULL
  SET r.significant = CASE
    WHEN r.adjusted_p_value < dm.p_value_threshold THEN 'significant'
    ELSE 'not_significant'
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
  OPTIONAL MATCH (t)<-[:Tcdb_family_is_a_tcdb_family*0..]-(desc:TcdbFamily)<-[:Gene_has_tcdb_family]-(g:Gene)
  WITH t, mc, count(DISTINCT g) AS gc,
       count(DISTINCT CASE WHEN desc = t THEN g END) AS dgc,
       collect(DISTINCT g.organism_name) AS orgs
  OPTIONAL MATCH (t)-[:Tcdb_family_transports_metabolite]->(m:Metabolite)
  WITH t, mc, gc, dgc, orgs, count(DISTINCT m) AS metc
  SET t.member_count = mc,
      t.gene_count = gc,
      t.direct_gene_count = dgc,
      t.organism_count = size([x IN orgs WHERE x IS NOT NULL]),
      t.metabolite_count = metc
} IN TRANSACTIONS OF 1000 ROWS;

// The old TcdbFamily promiscuity flag was DELETED (2026-08-16, spec §3 R3 +
// §9.8): it restated a threshold — level >= 2 AND metabolite_count >= 50 —
// over metabolite_count, a count the node already publishes. Consumers
// derive it themselves rather than the KG storing a predicate. History for
// context: it flagged a family as transporting MANY DISTINCT SUBSTRATES (so
// inferring what a member gene moves from family membership is weak),
// consumed by explorer family_inferred-dominance warnings (KG-MET-006).
// SUBSTRATE BREADTH ONLY — a `gene_count >= 500` arm was briefly added
// 2026-08-07 and reverted the same day (answers a different question and
// flagged substrate-poor families like 9.B.34). LEVEL-GATED (level >= 2,
// tc_family and deeper) because substrate counts scale mechanically with
// hierarchy level (median metabolite_count is 153 at tc_class vs 1 at
// tc_family — an unrestricted rule fired on 5 of 7 tc_class and 7 of 34
// tc_subclass nodes vacuously). Threshold >= 50 sat at ~p99 within the levels
// it applied to. Formerly flagged 13 families/subfamilies, all textbook
// multi-substrate transporters: ABC Superfamily 3.A.1 (554 substrates), MFS
// 2.A.1 (476), DMT 2.A.7, RND 2.A.6, MOP flippase 2.A.66, APC 2.A.3, P-type
// ATPase 3.A.3. See the one internal read of this threshold, inlined into
// transport_substrate_resolution below.

// ── CazyFamily computed properties ───────────────────────────────────────────

// gene_count + organism_count: subtree traversal (descendants ∪ self via *0..)
MATCH (c:CazyFamily)
CALL {
  WITH c
  OPTIONAL MATCH (c)<-[:Cazy_family_is_a_cazy_family*0..]-(desc:CazyFamily)<-[:Gene_has_cazy_family]-(g:Gene)
  WITH c, count(DISTINCT g) AS gc,
       count(DISTINCT CASE WHEN desc = c THEN g END) AS dgc,
       collect(DISTINCT g.organism_name) AS orgs
  SET c.gene_count = gc,
      c.direct_gene_count = dgc,
      c.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 1000 ROWS;

// ── InterproEntry computed properties (InterProScan; hierarchical ontology) ───
// gene_count / organism_count are SUBTREE (descendants ∪ self via *0..), the
// same semantics as every other hierarchical ontology (KG-SYNC-005 / ONT-009 —
// was DIRECT until 2026-08-27). direct_gene_count keeps the DIRECT number
// (genes with an edge to this exact entry) — the correct per-term count for
// (type,level)-stratified ORA. The is-a hierarchy is sparse (~1.6K edges over
// ~13K nodes, ~86% level-0), so the two differ on few nodes.
// member_count = direct child entries (structural, like TcdbFamily).
MATCH (e:InterproEntry)
CALL {
  WITH e
  OPTIONAL MATCH (child:InterproEntry)-[:Interpro_entry_is_a_interpro_entry]->(e)
  WITH e, count(child) AS mc
  OPTIONAL MATCH (e)<-[:Interpro_entry_is_a_interpro_entry*0..]-(desc:InterproEntry)<-[:Gene_has_interpro_entry]-(g:Gene)
  WITH e, mc, count(DISTINCT g) AS gc,
       count(DISTINCT CASE WHEN desc = e THEN g END) AS dgc,
       collect(DISTINCT g.organism_name) AS orgs
  SET e.member_count = mc,
      e.gene_count = gc,
      e.direct_gene_count = dgc,
      e.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 1000 ROWS;

// The old InterproEntry promiscuity flag was DELETED (2026-08-16, spec §3
// R3 + §9.8): it restated a threshold — gene_count >= 1000 — over
// gene_count, a count the node already publishes. Formerly flagged
// ultra-common entries (broad domains / superfamilies present in a large
// share of genes) so a (type, level)-stratified ORA could down-weight them
// (design spec §8); consumers apply that cutoff to gene_count directly now.

// ── NcbifamFamily computed properties (flat; direct counts) ───────────────────
MATCH (n:NcbifamFamily)
CALL {
  WITH n
  OPTIONAL MATCH (n)<-[:Gene_has_ncbifam_family]-(g:Gene)
  WITH n, count(DISTINCT g) AS gc, collect(DISTINCT g.organism_name) AS orgs
  SET n.gene_count = gc,
      n.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 1000 ROWS;

// ── MeropsFamily computed properties (hierarchical; subtree counts) ───────────
// gene_count + organism_count: subtree traversal (descendants ∪ self via *0..,
// CAZy pattern — a clan counts every gene under its families/subfamilies).
// peptidase_gene_count is the number to use for "how many proteases": it keeps
// only call_class='peptidase' edges, excluding inhibitor families and
// nonpeptidase_homolog calls (fold evidence, not protease evidence). The gap
// between gene_count and peptidase_gene_count is itself a signal — families
// like C26/C44 are dominated by catalytically dead relatives.
// member_count = direct child nodes (structural, TcdbFamily convention).
MATCH (m:MeropsFamily)
CALL {
  WITH m
  OPTIONAL MATCH (child:MeropsFamily)-[:Merops_family_is_a_merops_family]->(m)
  WITH m, count(child) AS mc
  OPTIONAL MATCH (m)<-[:Merops_family_is_a_merops_family*0..]-(desc:MeropsFamily)<-[r:Gene_has_merops_family]-(g:Gene)
  WITH m, mc, count(DISTINCT g) AS gc,
       count(DISTINCT CASE WHEN desc = m THEN g END) AS dgc,
       count(DISTINCT CASE WHEN r.call_class = 'peptidase' THEN g END) AS pgc,
       collect(DISTINCT g.organism_name) AS orgs,
       collect(DISTINCT CASE WHEN r.call_class = 'peptidase' THEN g.organism_name END) AS porgs
  SET m.member_count = mc,
      m.gene_count = gc,
      m.direct_gene_count = dgc,
      m.peptidase_gene_count = pgc,
      m.organism_count = size([x IN orgs WHERE x IS NOT NULL]),
      m.peptidase_organism_count = size([x IN porgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 1000 ROWS;

// Gene_has_merops_family.pfam_support: whether a Pfam on this gene is curated
// into this MEROPS family via Merops_family_has_pfam_domain (interpro.txt).
// Sound direction only — corroborates a known family; never assigns peptidase
// identity backward from a domain. R5 string pair, tcdb pfam_support pattern.
// The bridge attaches at family level (level 1); subfamily-attached gene edges
// check their parent family (single is-a hop).
CALL {
  MATCH (g:Gene)-[r:Gene_has_merops_family]->(m:MeropsFamily)
  WITH r, g, m,
       EXISTS {
         MATCH (g)-[:Gene_has_pfam]->(:Pfam)<-[:Merops_family_has_pfam_domain]-(fam:MeropsFamily)
         WHERE fam = m OR (m)-[:Merops_family_is_a_merops_family]->(fam)
       } AS pfam_ok
  SET r.pfam_support = CASE WHEN pfam_ok THEN 'corroborated' ELSE 'uncorroborated' END
} IN TRANSACTIONS OF 1000 ROWS;

// Gene_has_merops_family.evidence_score (KG-SYNC-005, ONT-006; R4 — one score
// name per concept, float in [0,1]). TWO placement-confidence signals:
//   +1  tier <= 2                      — identity >= 40% (subfamily/identifier band)
//   +1  pfam_support = 'corroborated'  — a Pfam on this gene is curated into the family
// call_class is deliberately NOT a signal: it is the verdict axis (active
// peptidase?), orthogonal to placement confidence — a confidently placed dead
// homolog must not score lower for being honest. Scale is {0, 0.5, 1};
// round(score * 2) recovers the fired count. Advisory, never a filter, and NOT
// calibrated against the 5-signal TCDB score.
CALL {
  MATCH ()-[r:Gene_has_merops_family]->()
  SET r.evidence_score = round(
      ( (CASE WHEN coalesce(r.tier <= 2, false) THEN 1 ELSE 0 END)
      + (CASE WHEN r.pfam_support = 'corroborated' THEN 1 ELSE 0 END) ) / 2.0, 3)
} IN TRANSACTIONS OF 1000 ROWS;

// Gene.merops_evidence_score_max (ONT-013): mirror of tcdb_evidence_score_max.
// SPARSE — only on genes with a Gene_has_merops_family edge; coalesce(x, -1.0)
// for a total order. Read with Gene.merops_classes: the max says how well the
// best call is PLACED, merops_classes says whether any call is an active peptidase.
CALL {
  MATCH (g:Gene)-[r:Gene_has_merops_family]->()
  WITH g, max(r.evidence_score) AS best
  SET g.merops_evidence_score_max = best
} IN TRANSACTIONS OF 1000 ROWS;

// ── gene_count / organism_count on the eggNOG-era ontologies (KG-SYNC-005, ONT-015)
// One semantics everywhere: SUBTREE on every hierarchical label (plus
// direct_gene_count = genes attached to this exact node), DIRECT on flat labels.
// GO walks is_a ∪ part_of (what `level` and the explorer's genes_by_ontology
// both walk; `regulates` is not in the graph). GO and KEGG are DAGs: a gene
// reachable through two parents is counted once per node, so sibling counts do
// NOT sum to the parent's. Descendants are collected DISTINCT first and genes
// matched second, so DAG path multiplicity never multiplies the gene rows.
// PfamClan (genes reach it only through member Pfams) and BriteCategory (only
// through KO leaves; already counted above) get no direct_gene_count — it would
// be a constant 0.

MATCH (n:BiologicalProcess)
CALL {
  WITH n
  OPTIONAL MATCH (n)<-[:Biological_process_is_a_biological_process|Biological_process_part_of_biological_process*0..]-(desc:BiologicalProcess)
  WITH n, collect(DISTINCT desc) AS descs
  UNWIND descs AS d
  OPTIONAL MATCH (d)<-[:Gene_involved_in_biological_process]-(g:Gene)
  WITH n, count(DISTINCT g) AS gc,
       count(DISTINCT CASE WHEN d = n THEN g END) AS dgc,
       collect(DISTINCT g.organism_name) AS orgs
  SET n.gene_count = gc,
      n.direct_gene_count = dgc,
      n.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 100 ROWS;

MATCH (n:MolecularFunction)
CALL {
  WITH n
  OPTIONAL MATCH (n)<-[:Molecular_function_is_a_molecular_function|Molecular_function_part_of_molecular_function*0..]-(desc:MolecularFunction)
  WITH n, collect(DISTINCT desc) AS descs
  UNWIND descs AS d
  OPTIONAL MATCH (d)<-[:Gene_enables_molecular_function]-(g:Gene)
  WITH n, count(DISTINCT g) AS gc,
       count(DISTINCT CASE WHEN d = n THEN g END) AS dgc,
       collect(DISTINCT g.organism_name) AS orgs
  SET n.gene_count = gc,
      n.direct_gene_count = dgc,
      n.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 100 ROWS;

MATCH (n:CellularComponent)
CALL {
  WITH n
  OPTIONAL MATCH (n)<-[:Cellular_component_is_a_cellular_component|Cellular_component_part_of_cellular_component*0..]-(desc:CellularComponent)
  WITH n, collect(DISTINCT desc) AS descs
  UNWIND descs AS d
  OPTIONAL MATCH (d)<-[:Gene_located_in_cellular_component]-(g:Gene)
  WITH n, count(DISTINCT g) AS gc,
       count(DISTINCT CASE WHEN d = n THEN g END) AS dgc,
       collect(DISTINCT g.organism_name) AS orgs
  SET n.gene_count = gc,
      n.direct_gene_count = dgc,
      n.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 100 ROWS;

MATCH (n:EcNumber)
CALL {
  WITH n
  OPTIONAL MATCH (n)<-[:Ec_number_is_a_ec_number*0..]-(desc:EcNumber)
  WITH n, collect(DISTINCT desc) AS descs
  UNWIND descs AS d
  OPTIONAL MATCH (d)<-[:Gene_catalyzes_ec_number]-(g:Gene)
  WITH n, count(DISTINCT g) AS gc,
       count(DISTINCT CASE WHEN d = n THEN g END) AS dgc,
       collect(DISTINCT g.organism_name) AS orgs
  SET n.gene_count = gc,
      n.direct_gene_count = dgc,
      n.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 100 ROWS;

MATCH (n:KeggTerm)
CALL {
  WITH n
  OPTIONAL MATCH (n)<-[:Kegg_term_is_a_kegg_term*0..]-(desc:KeggTerm)
  WITH n, collect(DISTINCT desc) AS descs
  UNWIND descs AS d
  OPTIONAL MATCH (d)<-[:Gene_has_kegg_ko]-(g:Gene)
  WITH n, count(DISTINCT g) AS gc,
       count(DISTINCT CASE WHEN d = n THEN g END) AS dgc,
       collect(DISTINCT g.organism_name) AS orgs
  // direct_gene_count only on KO leaves: genes attach to KOs alone, so on
  // pathway / subcategory / category nodes it is 0 by construction and is
  // OMITTED (BriteCategory / PfamClan precedent; DOC-006, 2026-08-29).
  SET n.gene_count = gc,
      n.direct_gene_count = CASE WHEN n.level_kind = 'ko' THEN dgc ELSE null END,
      n.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 100 ROWS;

MATCH (n:CyanorakRole)
CALL {
  WITH n
  OPTIONAL MATCH (n)<-[:Cyanorak_role_is_a_cyanorak_role*0..]-(desc:CyanorakRole)
  WITH n, collect(DISTINCT desc) AS descs
  UNWIND descs AS d
  OPTIONAL MATCH (d)<-[:Gene_has_cyanorak_role]-(g:Gene)
  WITH n, count(DISTINCT g) AS gc,
       count(DISTINCT CASE WHEN d = n THEN g END) AS dgc,
       collect(DISTINCT g.organism_name) AS orgs
  SET n.gene_count = gc,
      n.direct_gene_count = dgc,
      n.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 100 ROWS;

MATCH (n:PfamClan)
CALL {
  WITH n
  OPTIONAL MATCH (n)<-[:Pfam_in_pfam_clan]-(:Pfam)<-[:Gene_has_pfam]-(g:Gene)
  WITH n, count(DISTINCT g) AS gc, collect(DISTINCT g.organism_name) AS orgs
  SET n.gene_count = gc,
      n.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 100 ROWS;

MATCH (n:Pfam)
CALL {
  WITH n
  OPTIONAL MATCH (n)<-[:Gene_has_pfam]-(g:Gene)
  WITH n, count(DISTINCT g) AS gc, collect(DISTINCT g.organism_name) AS orgs
  SET n.gene_count = gc,
      n.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 1000 ROWS;

// TigrRole — two-level since 2026-08-29: subtree gene_count over
// Tigr_role_is_a_tigr_role (CyanorakRole pattern) + direct_gene_count +
// ncbifam_family_count (incoming Ncbifam_family_has_tigr_role, subtree).
MATCH (n:TigrRole)
CALL {
  WITH n
  OPTIONAL MATCH (n)<-[:Tigr_role_is_a_tigr_role*0..]-(desc:TigrRole)
  WITH n, collect(DISTINCT desc) AS descs
  UNWIND descs AS d
  OPTIONAL MATCH (d)<-[:Gene_has_tigr_role]-(g:Gene)
  WITH n, descs, count(DISTINCT g) AS gc,
       count(DISTINCT CASE WHEN d = n THEN g END) AS dgc,
       collect(DISTINCT g.organism_name) AS orgs
  UNWIND descs AS d2
  OPTIONAL MATCH (d2)<-[:Ncbifam_family_has_tigr_role]-(f:NcbifamFamily)
  WITH n, gc, dgc, orgs, count(DISTINCT f) AS fc
  SET n.gene_count = gc,
      n.direct_gene_count = dgc,
      n.organism_count = size([x IN orgs WHERE x IS NOT NULL]),
      n.ncbifam_family_count = fc
} IN TRANSACTIONS OF 100 ROWS;

MATCH (n:CogFunctionalCategory)
CALL {
  WITH n
  OPTIONAL MATCH (n)<-[:Gene_in_cog_category]-(g:Gene)
  WITH n, count(DISTINCT g) AS gc, collect(DISTINCT g.organism_name) AS orgs
  SET n.gene_count = gc,
      n.organism_count = size([x IN orgs WHERE x IS NOT NULL])
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
//   - cazy_family_count  (TCDB-S2)
//   (tcdb_family_count moved to the transport-arm statement below, where the
//    most-specific attachment set it now counts is already in hand — explorer
//    ask 2026-08-29-gene-overview-family-counts-asks)
//   - reaction_count     (KG-A1)
//   - catalyzed_metabolite_count  (KG-A2 — the catalysis arm; renamed from
//     metabolite_count for the alpha.7 cut, KG-SYNC-001)
//
// Each OPTIONAL MATCH is followed by a WITH aggregation so rows don't multiply.
// count(DISTINCT ...) is required for the count rollups because the chained
// metabolite OPTIONAL MATCHes would otherwise duplicate the parent edge per metabolite.
//
// BREAKING: metabolite_count was the UNION of the catalysis and transport arms.
// The arms have very different epistemics — catalysis p90 = 11 metabolites,
// transport p90 = 554 — because the step-6 rollup materialises every descendant's
// substrates onto each ancestor, so a gene annotated only at ABC superfamily
// 3.A.1 inherited all 554 substrates ABC transporters have ever been curated for.
// 23,137 genes had transport evidence ONLY, so for the majority of genes the
// stored number was entirely the inflated arm with nothing signalling it. The
// transport arm now lives in transported_metabolite_count (next statement).
// The catalysis arm is RENAMED catalyzed_metabolite_count (KG-SYNC-001): keeping
// the bare name after narrowing its meaning would hand every stale reader the
// narrowed number silently; the rename makes them fail loudly (null) instead,
// and each arm now names itself.
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[r2:Gene_has_cazy_family]->()
  WITH g, count(r2) AS cz_count
  // merops_family_count = ALL merops edges (routing, ungated — tcdb precedent).
  // merops_classes = distinct call_class values, the at-a-glance guard so a
  // gene whose only call is a dead homolog never reads as "1 protease".
  OPTIONAL MATCH (g)-[r2b:Gene_has_merops_family]->()
  WITH g, cz_count, count(r2b) AS mer_count,
       apoc.coll.sort([c IN collect(DISTINCT r2b.call_class) WHERE c IS NOT NULL]) AS mer_classes
  OPTIONAL MATCH (g)-[r3:Gene_catalyzes_reaction]->(rx:Reaction)
  OPTIONAL MATCH (rx)-[:Reaction_has_metabolite]->(m_cat:Metabolite)
  WITH g, cz_count, mer_count, mer_classes,
       count(DISTINCT r3) AS rxn_count,
       count(DISTINCT m_cat) AS cat_met_count
  SET g.cazy_family_count = cz_count,
      g.merops_family_count = mer_count,
      g.merops_classes = mer_classes,
      g.reaction_count = rxn_count,
      g.catalyzed_metabolite_count = cat_met_count
} IN TRANSACTIONS OF 1000 ROWS;

// ── Gene_has_tcdb_family.attachment_depth (KG-SYNC-005, ONT-010) ─────────────
// Materializes the deepest-attachment predicate ONCE. 'most_specific' = no other
// Gene_has_tcdb_family edge of the SAME gene lands on a descendant of this node;
// 'superseded' = one does (e.g. an eggNOG 3.A.1 call next to a diamond 3.A.1.14
// call — a less specific call, NOT a wrong one). A structural fact, not a
// threshold (R3 does not apply); R5 string pair mirroring substrate_depth.
//
// Checking DIRECT ancestry is not enough — a gene may be annotated at 3.A.1 and
// 3.A.1.14.2 with no edge to the intervening 3.A.1.14 — hence *1..4 (TCDB is 5
// levels, so 4 hops is the maximum ancestor distance).
//
// The three transport-arm consumers (Gene.transported_metabolite_count,
// Metabolite.transporter_gene_count, Organism_has_metabolite) read this property
// instead of re-deriving the predicate, so they cannot drift from each other or
// from what a consumer sees on the edge.
CALL {
  MATCH (g:Gene)-[r:Gene_has_tcdb_family]->(t:TcdbFamily)
  WITH r, EXISTS {
    MATCH (g)-[:Gene_has_tcdb_family]->(d:TcdbFamily)
    WHERE (d)-[:Tcdb_family_is_a_tcdb_family*1..4]->(t)
  } AS superseded
  SET r.attachment_depth = CASE WHEN superseded THEN 'superseded' ELSE 'most_specific' END
} IN TRANSACTIONS OF 1000 ROWS;

// Gene transport arm: tcdb_family_count + transported_metabolite_count +
// transport_substrate_resolution.
//
// tcdb_family_count (TCDB-S1) = the gene's MOST-SPECIFIC attachments only
// (2026-08-29, explorer ask gene-overview-family-counts; was every edge).
// A superseded ancestor edge (3.A.1 beside the gene's own 3.A.1.14) is a less
// specific restatement of the same membership, not a second family — counting
// it over-read 7,045 genes (PMM0392: 8 → 7). Same projection the two other
// transport-arm counts already use. Still NOT tier-gated: depth and tier are
// different axes, an uncorroborated most_specific DIAMOND hit still counts.
//
// MOST-SPECIFIC ATTACHMENTS ONLY (attachment_depth above). 6,950 genes are
// annotated at both an ancestor and its own descendant (e.g. both 3.A.1 and
// 3.A.1.14); unioning across all of a gene's attachments pulled in the
// ancestor's full rolled-up substrate set even though a more specific call
// existed. Restricting to the most specific attachments keeps 26,813 of 26,894
// genes (99.7%) and cuts p90 from 554 to 97.
//
// NOT tier-gated, deliberately. Tier and substrate resolution are orthogonal:
// 11,871 genes are 'resolved' yet tier-3-only (narrow 2.A.x secondary carriers
// where remote homology could not justify a subfamily call). A tier gate would
// discard exactly those while keeping eggNOG's equally-lumping tc_family edges,
// since eggNOG carries no tier — an artifact of which tool called it. Tier
// already gates annotation_types/annotation_quality and rolls up as
// Gene.tcdb_evidence_score_max; that is its home.
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[:Gene_has_tcdb_family {attachment_depth: 'most_specific'}]->(t:TcdbFamily)
  OPTIONAL MATCH (t)-[:Tcdb_family_transports_metabolite]->(m_tr:Metabolite)
  WITH g,
       count(DISTINCT m_tr) AS tr_met_count,
       count(DISTINCT t) AS n_deepest,
       // Breadth threshold inlined (spec §3 R3): the deleted TcdbFamily
       // promiscuity flag restated a predicate over metabolite_count, which
       // the node already publishes. Consumers apply their own cutoff to the
       // count; this is the KG's, and it lives in exactly one place.
       collect(DISTINCT (coalesce(t.level, 0) >= 2
                         AND coalesce(t.metabolite_count, 0) >= 50)) AS breadth
  SET g.tcdb_family_count = n_deepest,
      g.transported_metabolite_count = tr_met_count,
      // null REMOVES the property, keeping it sparse: absent means "no TCDB
      // edge at all", which must stay distinguishable from a weak-but-present
      // substrate claim.
      g.transport_substrate_resolution =
        CASE WHEN n_deepest = 0 THEN null
             WHEN any(x IN breadth WHERE x = false) THEN 'resolved'
             ELSE 'family_inferred' END
} IN TRANSACTIONS OF 1000 ROWS;

// ── TCDB evidence score (per Gene_has_tcdb_family edge) ──────────────────────
// An ADVISORY ranking aid over five independent supporting signals. It never
// drops an edge and nothing filters on it silently -- the components are stored
// alongside the total so a consumer can always see WHY, and re-weight without
// re-deriving.
//
// This is deliberately NOT shaped like the deleted `filter_action` chain: it is
// additive rather than first-match-wins, every edge is scored on its own
// evidence (no sibling dependence), and it carries no uncalibrated thresholds --
// `tier` is the only cut and it is already a principled, sibling-independent gate.
//
//   +1  curated                — eggNOG called it (ortholog transfer from a curated DB)
//   +1  agrees_across_sources  — eggNOG and diamond concur
//   +1  tier <= 2              — strong direct sequence evidence (identity >= 40%)
//   +1  pfam_corroborated      — a Pfam on this gene is curated into this TC family
//   +1  go_corroborated        — a GO term on this gene is curated onto this TC family
//
// AGREEMENT IS HIERARCHICAL, not exact-node. eggNOG names subfamilies while
// diamond's tier-3 truncation names the parent family, so the two usually concur
// at DIFFERENT depths: exact same-node agreement covers only 3,641 edges, while
// including ancestor/descendant agreement reaches 21,684 (40.3%). Scoring on
// `size(sources)=2` alone would miss 83% of the real corroboration.
//
// Pfam and GO are kept separate because they are not redundant: 11,565 edges
// carry exactly one of the two (P(GO|Pfam)=77% vs P(GO|not Pfam)=20%).
CALL {
  MATCH (g:Gene)-[r:Gene_has_tcdb_family]->(t:TcdbFamily)
  WITH r, g, t,
       ('eggnog' IN r.sources) AS curated,
       (size(r.sources) = 2 OR EXISTS {
          MATCH (g)-[r2:Gene_has_tcdb_family]->(t2:TcdbFamily)
          WHERE t2 <> t AND any(s IN r2.sources WHERE NOT s IN r.sources)
            AND ( (t)-[:Tcdb_family_is_a_tcdb_family*1..4]->(t2)
               OR (t2)-[:Tcdb_family_is_a_tcdb_family*1..4]->(t) )
       }) AS agree,
       coalesce(r.tier <= 2, false) AS strong_seq,
       EXISTS {
         MATCH (g)-[:Gene_has_pfam]->(:Pfam)<-[:Tcdb_family_has_pfam_domain]-(t)
       } AS pfam_ok,
       EXISTS {
         MATCH (g)-[:Gene_involved_in_biological_process|Gene_enables_molecular_function|Gene_located_in_cellular_component]->(o)
               <-[:Tcdb_family_involved_in_biological_process|Tcdb_family_enables_molecular_function|Tcdb_family_located_in_cellular_component]-(t)
       } AS go_ok
  SET r.source_agreement = CASE WHEN agree   THEN 'both_sources' ELSE 'single_source' END,
      r.pfam_support     = CASE WHEN pfam_ok THEN 'corroborated' ELSE 'uncorroborated' END,
      r.go_support       = CASE WHEN go_ok   THEN 'corroborated' ELSE 'uncorroborated' END,
      r.evidence_score   = round(
          ( (CASE WHEN curated THEN 1 ELSE 0 END)
          + (CASE WHEN agree THEN 1 ELSE 0 END)
          + (CASE WHEN strong_seq THEN 1 ELSE 0 END)
          + (CASE WHEN pfam_ok THEN 1 ELSE 0 END)
          + (CASE WHEN go_ok THEN 1 ELSE 0 END) ) / 5.0, 3)
} IN TRANSACTIONS OF 1000 ROWS;

// Gene.tcdb_evidence_score_max: the strongest TCDB claim this gene has (float 0-1).
// Answers "how confident am I that this gene is a transporter at all", where the
// edge-level score answers "how confident am I in THIS particular assignment".
// Worth materializing: 9,792 of 30,076 TCDB-annotated genes (32.6%) carry several
// calls at DIFFERENT scores, so max() is not just a copy of the single edge.
//
// SPARSE BY DESIGN — set only on genes that actually have a Gene_has_tcdb_family
// edge. A gene with no TCDB evidence must stay distinguishable from one whose
// evidence is weak: writing 0 for both would collapse "we never found a
// transporter signal" into "we found a poor one", which is the exact conflation
// the tier gate and the advisory score exist to avoid. Absent means N/A; use
// coalesce(g.tcdb_evidence_score_max, -1.0) if a total order is needed.
//
// Advisory, like the edge score it aggregates: nothing filters on it.
CALL {
  MATCH (g:Gene)-[r:Gene_has_tcdb_family]->()
  WITH g, max(r.evidence_score) AS best
  SET g.tcdb_evidence_score_max = best
} IN TRANSACTIONS OF 1000 ROWS;

// Gene.interpro_entry_count: distinct InterPro entries per gene (routing signal;
// functional. Folded into annotation_types above but not annotation_quality — §5).
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[r:Gene_has_interpro_entry]->()
  WITH g, count(r) AS ic
  SET g.interpro_entry_count = ic
} IN TRANSACTIONS OF 1000 ROWS;

// Gene.ncbifam_family_count: distinct NCBIfam families per gene (routing signal;
// functional. Folded into annotation_types AND informative_annotation_types /
// annotation_quality above — §4.2).
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[r:Gene_has_ncbifam_family]->()
  WITH g, count(r) AS nc
  SET g.ncbifam_family_count = nc
} IN TRANSACTIONS OF 1000 ROWS;

// OrganismTaxon annotation-capability rollups (KG-SYNC-006 ORG-001): distinct
// genes per organism, scoped by Gene_belongs_to_organism like every sibling
// rollup. Dense 0. peptidase / nonpeptidase_homolog read Gene.merops_classes
// (a gene can be in both -- do not subtract); interpro = edge existence;
// ncbifam = Gene.ncbifam_family_count > 0 (retired families still count).
MATCH (o:OrganismTaxon)
OPTIONAL MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o)
WITH o,
     count(DISTINCT CASE WHEN 'peptidase' IN coalesce(g.merops_classes, []) THEN g END) AS pep,
     count(DISTINCT CASE WHEN 'nonpeptidase_homolog' IN coalesce(g.merops_classes, []) THEN g END) AS nonpep,
     count(DISTINCT CASE WHEN EXISTS { (g)-[:Gene_has_interpro_entry]->() } THEN g END) AS ipr,
     count(DISTINCT CASE WHEN coalesce(g.ncbifam_family_count, 0) > 0 THEN g END) AS nf
SET o.peptidase_gene_count = pep,
    o.nonpeptidase_homolog_gene_count = nonpep,
    o.interpro_gene_count = ipr,
    o.ncbifam_gene_count = nf;

// ── Metabolism rollups ────────────────────────────────────────────────────

// Reaction.gene_count, organism_count, organisms[]
CALL {
  MATCH (r:Reaction)<-[:Gene_catalyzes_reaction]-(g:Gene)
  WITH r, count(DISTINCT g) AS gene_count, collect(DISTINCT g.organism_name) AS organisms
  SET r.gene_count = gene_count,
      r.organism_count = size(organisms),
      r.organisms = organisms
} IN TRANSACTIONS OF 1000 ROWS;

// Metabolite.catalyst_gene_count + transporter_gene_count (per gene-link arm).
// Substrate edges are rolled up to every ancestor in the adapter, so the
// transport arm is a 1-hop traversal at any TcdbFamily level the gene is
// annotated at. organism_count is computed below from the materialized
// Organism_has_metabolite edge so size(organism_names) == organism_count
// is invariant by construction (KG-A8).
// BREAKING, mirroring Gene.catalyzed_metabolite_count: the CATALYSIS arm is
// catalyst_gene_count (renamed from the bare gene_count, KG-SYNC-001 — where a
// node has two gene-link arms, each count names its arm); the transport arm is
// transporter_gene_count. Both ends of the transport relation use the SAME
// predicate — the gene's most-specific TC attachments (attachment_depth) — so
// Gene.transported_metabolite_count and Metabolite.transporter_gene_count are
// two projections of one (gene, metabolite) set and agree by construction.
CALL {
  MATCH (m:Metabolite)
  OPTIONAL MATCH (m)<-[:Reaction_has_metabolite]-(:Reaction)<-[:Gene_catalyzes_reaction]-(g_cat:Gene)
  WITH m, count(DISTINCT g_cat) AS cat_gene_count
  OPTIONAL MATCH (m)<-[:Tcdb_family_transports_metabolite]-(t:TcdbFamily)<-[:Gene_has_tcdb_family {attachment_depth: 'most_specific'}]-(g_tr:Gene)
  WITH m, cat_gene_count, count(DISTINCT g_tr) AS tr_gene_count
  SET m.catalyst_gene_count = cat_gene_count,
      m.transporter_gene_count = tr_gene_count
} IN TRANSACTIONS OF 1000 ROWS;

// Metabolite.transporter_count: distinct transporter systems at MAXIMAL DEPTH.
//
// Was `level_kind = 'tc_specificity'`, which made the count 0 for 1,218 of the
// 1,462 transported metabolites (83%). That filter dated from before the
// ancestor-only prune: only 466 of 11,263 substrate edges now sit on
// tc_specificity nodes, because genes mostly annotate at family/subfamily depth
// and the prune keeps no specificity node below them.
//
// substrate_depth = 'most_specific' (set by tcdb_adapter) means no kept child of this
// node also carries the substrate — so counting DISTINCT sources over those
// edges counts each transporter system once, at the finest resolution the pruned
// graph retains, without double-counting an ancestor together with its own
// descendant. All 1,462 transported metabolites now get a non-zero count.
CALL {
  MATCH (m:Metabolite)
  OPTIONAL MATCH (t:TcdbFamily)-[r:Tcdb_family_transports_metabolite]->(m)
    WHERE r.substrate_depth = 'most_specific'
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

// Materialize Organism_has_metabolite (transport arm).
// Uses the SAME most-specific-attachment property (attachment_depth) as
// Gene.transported_metabolite_count so the organism edge and the gene scalar
// cannot disagree. Without it, one gene
// annotated at ABC superfamily 3.A.1 gave its whole organism all 554 ABC
// substrates, which is how every organism came to "have" 63% of all metabolites.
CALL {
  MATCH (o:OrganismTaxon)<-[:Gene_belongs_to_organism]-(g:Gene)
        -[:Gene_has_tcdb_family {attachment_depth: 'most_specific'}]->(t:TcdbFamily)
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
// BREAKING, mirroring Gene / Metabolite: the CATALYSIS arm is
// catalyzed_metabolite_count (renamed from the bare metabolite_count,
// KG-SYNC-001); the transport arm is transported_metabolite_count. (The
// measurement arm already had its own scalar, measured_metabolite_count, so
// this completes the three-way split.) evidence_sources on the edge is the
// discriminator, so neither arm needs re-traversal here.
CALL {
  MATCH (o:OrganismTaxon)-[r:Organism_has_metabolite]->(m:Metabolite)
  WITH o,
       count(DISTINCT CASE WHEN 'metabolism' IN r.evidence_sources THEN m END) AS cat_count,
       count(DISTINCT CASE WHEN 'transport'  IN r.evidence_sources THEN m END) AS tr_count
  SET o.catalyzed_metabolite_count = cat_count,
      o.transported_metabolite_count = tr_count
} IN TRANSACTIONS OF 1000 ROWS;

CALL {
  MATCH (o:OrganismTaxon)<-[:Gene_belongs_to_organism]-(:Gene)
        -[:Gene_catalyzes_reaction]->(r:Reaction)
  WITH o, count(DISTINCT r) AS reaction_count
  SET o.reaction_count = reaction_count
} IN TRANSACTIONS OF 1000 ROWS;

// ── Chemistry slice-1 rollups (KG-A4) ───────────────────────
// (Note: KG-A1 Gene.reaction_count and KG-A2 Gene.catalyzed_metabolite_count are
// folded into the combined Gene metabolism + ontology rollup statement above.)

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

// MetaboliteAssay numeric ranks: per-assay, only when rankable='rankable'.
// Mirrors DerivedMetric pattern (per-assay scope, deterministic Metabolite.id tiebreaker,
// pinned bucket thresholds 90/75/25).
MATCH (a:MetaboliteAssay {rankable: 'rankable'})
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
       sum(CASE WHEN r.flag_value='detected'     THEN 1 ELSE 0 END) AS t,
       sum(CASE WHEN r.flag_value='not_detected' THEN 1 ELSE 0 END) AS f
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

// ── Publication "discusses" edge rollups (literature-index routing signals) ──

// Publication: distinct discussed genes / pathways
MATCH (p:Publication)
CALL {
  WITH p
  OPTIONAL MATCH (p)-[:Publication_discusses_gene]->(g:Gene)
  WITH p, count(DISTINCT g) AS dgc
  OPTIONAL MATCH (p)-[:Publication_discusses_kegg_pathway]->(k:KeggTerm)
  WITH p, dgc, count(DISTINCT k) AS dpc
  SET p.discussed_gene_count = dgc,
      p.discussed_pathway_count = dpc
} IN TRANSACTIONS OF 1000 ROWS;

// Gene: how many publications discuss this gene (router signal)
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (p:Publication)-[:Publication_discusses_gene]->(g)
  WITH g, count(DISTINCT p) AS dipc
  SET g.discussed_in_publication_count = dipc
} IN TRANSACTIONS OF 1000 ROWS;

// ─────────────────────────────────────────────────────────────────────────────
// Schema_info release metadata.
// In Docker this is run by post-import.sh as a SEPARATE cypher-shell invocation
// with -P params from the environment (KG_RELEASE_VERSION etc. — see "Group 4").
// For standalone `cypher-shell -f` runs, the :param lines below supply dev
// defaults; override with -P / :param to stamp a real release. The MATCH/SET
// logic below is byte-identical to the block in post-import.sh.
// Counts are computed (not hardcoded) so they track data drift.
// ─────────────────────────────────────────────────────────────────────────────
:param version          => '0.0.0-dev'
:param git_sha          => 'unknown'
:param git_sha_short    => 'unknown'
:param git_branch       => 'unknown'
:param git_dirty        => 'unknown'
:param mcp_min_version  => '0.1.0'
:param deployment_role  => 'local-dev'
:param release_notes_url => ''
:param release_highlights => ''
:param release_breaking   => ''
:param vocab_hash         => ''

MATCH (s:Schema_info {id: 'schema_info'})
SET s.version           = coalesce($version, '0.0.0-dev'),
    s.built_at          = toString(datetime()),
    s.git_sha           = coalesce($git_sha, 'unknown'),
    s.git_sha_short     = coalesce($git_sha_short, 'unknown'),
    s.git_branch        = coalesce($git_branch, 'unknown'),
    s.git_dirty         = coalesce($git_dirty, 'unknown'),
    s.mcp_min_version   = coalesce($mcp_min_version, '0.1.0'),
    s.deployment_role   = coalesce($deployment_role, 'local-dev'),
    s.release_notes_url = coalesce($release_notes_url, ''),
    // Empty string → real null property; see matching block in post-import.sh.
    s.release_highlights = CASE WHEN coalesce($release_highlights, '') = '' THEN null ELSE $release_highlights END,
    s.release_breaking   = CASE WHEN coalesce($release_breaking, '')   = '' THEN null ELSE $release_breaking   END,
    s.controlled_vocabularies_hash = CASE WHEN coalesce($vocab_hash, '') = ''
                                          THEN null ELSE $vocab_hash END
WITH s
SET s.paper_count           = COUNT { (:Publication) },
    s.experiment_count      = COUNT { (:Experiment) },
    s.gene_count            = COUNT { (:Gene) },
    s.organism_count        = COUNT { (:OrganismTaxon) },
    s.expression_edge_count = COUNT { ()-[:Changes_expression_of]->() };
