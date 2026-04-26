#!/bin/bash
# Capture a deterministic snapshot of post-import-computed properties.
#
# Usage:
#   scripts/post-import-validate.sh > baseline.txt   # before refactor
#   scripts/post-import-validate.sh > after.txt      # after refactor
#   diff baseline.txt after.txt                       # must be empty
#
# Requires the `deploy` docker container to be running (localhost:7687).

set -euo pipefail

CYPHER() {
  docker exec -i deploy cypher-shell -a bolt://localhost:7687 --format plain
}

section() {
  printf '\n\n======== %s ========\n' "$1"
}

# ----------------------------------------------------------------------------
# SHOW INDEXES — every index name, type, label, properties, state
# ----------------------------------------------------------------------------
section "INDEXES"
CYPHER <<'CYPHER'
SHOW INDEXES
YIELD name, type, entityType, labelsOrTypes, properties, state
WHERE type <> 'LOOKUP'
RETURN name, type, entityType, labelsOrTypes, properties, state
ORDER BY name;
CYPHER

# ----------------------------------------------------------------------------
# Small-table full dumps (deterministic sort)
# ----------------------------------------------------------------------------
section "EXPERIMENT (full dump)"
CYPHER <<'CYPHER'
MATCH (e:Experiment)
RETURN e.id AS id,
       e.gene_count AS gene_count,
       e.distinct_gene_count AS distinct_gene_count,
       e.significant_up_count AS sig_up,
       e.significant_down_count AS sig_down,
       e.time_point_count AS tp_count,
       e.time_point_labels AS tp_labels,
       e.time_point_orders AS tp_orders,
       e.time_point_hours AS tp_hours,
       e.time_point_totals AS tp_totals,
       e.time_point_significant_up AS tp_sig_up,
       e.time_point_significant_down AS tp_sig_down,
       e.time_point_growth_phases AS tp_gps,
       e.growth_phases AS growth_phases,
       e.clustering_analysis_count AS ca_count,
       e.cluster_types AS cluster_types,
       e.cluster_count AS cluster_count
ORDER BY e.id;
CYPHER

section "PUBLICATION (full dump)"
CYPHER <<'CYPHER'
MATCH (p:Publication)
RETURN p.id AS id,
       p.experiment_count AS exp_count,
       apoc.coll.sort(coalesce(p.treatment_types, [])) AS treatment_types,
       apoc.coll.sort(coalesce(p.background_factors, [])) AS background_factors,
       apoc.coll.sort(coalesce(p.omics_types, [])) AS omics_types,
       apoc.coll.sort(coalesce(p.growth_phases, [])) AS growth_phases,
       apoc.coll.sort(coalesce(p.organisms, [])) AS organisms,
       p.clustering_analysis_count AS ca_count,
       apoc.coll.sort(coalesce(p.cluster_types, [])) AS cluster_types,
       p.cluster_count AS cluster_count
ORDER BY p.id;
CYPHER

section "ORGANISMTAXON (full dump)"
CYPHER <<'CYPHER'
MATCH (o:OrganismTaxon)
RETURN o.id AS id,
       o.gene_count AS gene_count,
       o.publication_count AS pub_count,
       o.experiment_count AS exp_count,
       apoc.coll.sort(coalesce(o.treatment_types, [])) AS treatment_types,
       apoc.coll.sort(coalesce(o.omics_types, [])) AS omics_types,
       apoc.coll.sort(coalesce(o.background_factors, [])) AS background_factors,
       apoc.coll.sort(coalesce(o.growth_phases, [])) AS growth_phases,
       o.clustering_analysis_count AS ca_count,
       apoc.coll.sort(coalesce(o.cluster_types, [])) AS cluster_types,
       o.cluster_count AS cluster_count
ORDER BY o.id;
CYPHER

section "CLUSTERINGANALYSIS (full dump)"
CYPHER <<'CYPHER'
MATCH (ca:ClusteringAnalysis)
RETURN ca.id AS id,
       apoc.coll.sort(coalesce(ca.growth_phases, [])) AS growth_phases
ORDER BY ca.id;
CYPHER

section "GENECLUSTER (full dump)"
CYPHER <<'CYPHER'
MATCH (gc:GeneCluster)
RETURN gc.id AS id,
       gc.member_count AS member_count
ORDER BY gc.id;
CYPHER

section "BRITECATEGORY (full dump)"
CYPHER <<'CYPHER'
MATCH (b:BriteCategory)
RETURN b.id AS id,
       b.member_ko_count AS member_ko_count,
       b.gene_count AS gene_count,
       b.organism_count AS organism_count
ORDER BY b.id;
CYPHER

# ----------------------------------------------------------------------------
# Gene aggregates (81K nodes — aggregate only)
# ----------------------------------------------------------------------------
section "GENE AGGREGATES"
CYPHER <<'CYPHER'
MATCH (g:Gene)
RETURN count(g) AS total_genes,
       count(g.annotation_types) AS with_annotation_types,
       count(g.expression_edge_count) AS with_expression_edge_count,
       count(g.significant_up_count) AS with_sig_up,
       count(g.significant_down_count) AS with_sig_down,
       count(g.cluster_membership_count) AS with_cluster_memb,
       count(g.closest_ortholog_group_size) AS with_closest_og_size,
       count(g.closest_ortholog_genera) AS with_closest_og_genera,
       sum(g.expression_edge_count) AS sum_expression_edges,
       sum(g.significant_up_count) AS sum_sig_up,
       sum(g.significant_down_count) AS sum_sig_down,
       sum(g.cluster_membership_count) AS sum_cluster_memb,
       sum(g.closest_ortholog_group_size) AS sum_closest_og_size,
       min(g.expression_edge_count) AS min_expr_edges,
       max(g.expression_edge_count) AS max_expr_edges;
CYPHER

section "GENE annotation_types distribution"
CYPHER <<'CYPHER'
MATCH (g:Gene)
UNWIND coalesce(g.annotation_types, []) AS t
RETURN t AS annotation_type, count(*) AS gene_count
ORDER BY annotation_type;
CYPHER

section "GENE cluster_types distribution"
CYPHER <<'CYPHER'
MATCH (g:Gene)
UNWIND coalesce(g.cluster_types, []) AS t
RETURN t AS cluster_type, count(*) AS gene_count
ORDER BY cluster_type;
CYPHER

section "GENE closest_ortholog_genera distribution"
CYPHER <<'CYPHER'
MATCH (g:Gene)
UNWIND coalesce(g.closest_ortholog_genera, []) AS genus
RETURN genus, count(*) AS gene_count
ORDER BY genus;
CYPHER

# ----------------------------------------------------------------------------
# Edge aggregates
# ----------------------------------------------------------------------------
section "EXPRESSION_STATUS distribution"
CYPHER <<'CYPHER'
MATCH ()-[r:Changes_expression_of]->()
RETURN coalesce(r.expression_status, '<null>') AS expression_status, count(r) AS edge_count
ORDER BY expression_status;
CYPHER

section "RANK aggregates"
CYPHER <<'CYPHER'
MATCH ()-[r:Changes_expression_of]->()
RETURN count(r.rank_by_effect) AS with_rank_by_effect,
       min(r.rank_by_effect) AS min_rbe,
       max(r.rank_by_effect) AS max_rbe,
       sum(r.rank_by_effect) AS sum_rbe,
       count(r.rank_up) AS with_rank_up,
       min(r.rank_up) AS min_ru,
       max(r.rank_up) AS max_ru,
       sum(r.rank_up) AS sum_ru,
       count(r.rank_down) AS with_rank_down,
       min(r.rank_down) AS min_rd,
       max(r.rank_down) AS max_rd,
       sum(r.rank_down) AS sum_rd;
CYPHER

# ----------------------------------------------------------------------------
# Rank subsample — deterministic integer-only dump. Top 10 per (experiment,
# time_point_order) by each rank type. Ties broken by locus_tag so the
# sample is fully reproducible.
# ----------------------------------------------------------------------------
section "RANK SUBSAMPLE: top-10 rank_by_effect per (experiment, time_point_order)"
CYPHER <<'CYPHER'
MATCH (e:Experiment)-[r:Changes_expression_of]->(g:Gene)
WHERE r.rank_by_effect IS NOT NULL AND r.rank_by_effect <= 10
RETURN e.id AS experiment_id,
       r.time_point_order AS tp_order,
       g.locus_tag AS locus_tag,
       r.rank_by_effect AS rank
ORDER BY experiment_id, tp_order, rank, locus_tag;
CYPHER

section "RANK SUBSAMPLE: top-10 rank_up per (experiment, time_point_order)"
CYPHER <<'CYPHER'
MATCH (e:Experiment)-[r:Changes_expression_of]->(g:Gene)
WHERE r.rank_up IS NOT NULL AND r.rank_up <= 10
RETURN e.id AS experiment_id,
       r.time_point_order AS tp_order,
       g.locus_tag AS locus_tag,
       r.rank_up AS rank
ORDER BY experiment_id, tp_order, rank, locus_tag;
CYPHER

section "RANK SUBSAMPLE: top-10 rank_down per (experiment, time_point_order)"
CYPHER <<'CYPHER'
MATCH (e:Experiment)-[r:Changes_expression_of]->(g:Gene)
WHERE r.rank_down IS NOT NULL AND r.rank_down <= 10
RETURN e.id AS experiment_id,
       r.time_point_order AS tp_order,
       g.locus_tag AS locus_tag,
       r.rank_down AS rank
ORDER BY experiment_id, tp_order, rank, locus_tag;
CYPHER

section "DERIVEDMETRIC (full dump)"
CYPHER <<'CYPHER'
MATCH (dm:DerivedMetric)
RETURN dm.id AS id,
       dm.metric_type AS metric_type,
       dm.value_kind AS value_kind,
       dm.rankable AS rankable,
       dm.has_p_value AS has_p_value,
       dm.compartment AS compartment,
       dm.total_gene_count AS total_gene_count,
       apoc.coll.sort(coalesce(dm.growth_phases, [])) AS growth_phases
ORDER BY dm.id;
CYPHER

section "EXPERIMENT DM rollup dump"
CYPHER <<'CYPHER'
MATCH (e:Experiment)
WHERE e.derived_metric_count > 0
RETURN e.id AS id,
       e.reports_fold_change AS rfc,
       apoc.coll.sort(coalesce(e.reports_derived_metric_types, [])) AS metric_types,
       e.derived_metric_count AS dm_count,
       apoc.coll.sort(coalesce(e.derived_metric_value_kinds, [])) AS value_kinds,
       e.derived_metric_gene_count AS dm_gene_count
ORDER BY e.id;
CYPHER

section "PUBLICATION DM rollup dump"
CYPHER <<'CYPHER'
MATCH (p:Publication)
WHERE p.derived_metric_count > 0
RETURN p.id AS id,
       p.derived_metric_count AS dm_count,
       p.derived_metric_gene_count AS dm_gene_count,
       apoc.coll.sort(coalesce(p.compartments, [])) AS compartments,
       apoc.coll.sort(coalesce(p.derived_metric_types, [])) AS dm_types,
       apoc.coll.sort(coalesce(p.derived_metric_value_kinds, [])) AS dm_value_kinds
ORDER BY p.id;
CYPHER

section "ORGANISMTAXON DM rollup dump"
CYPHER <<'CYPHER'
MATCH (o:OrganismTaxon)
WHERE o.derived_metric_count > 0
RETURN o.id AS id,
       o.derived_metric_count AS dm_count,
       o.derived_metric_gene_count AS dm_gene_count,
       apoc.coll.sort(coalesce(o.compartments, [])) AS compartments,
       apoc.coll.sort(coalesce(o.derived_metric_types, [])) AS dm_types,
       apoc.coll.sort(coalesce(o.derived_metric_value_kinds, [])) AS dm_value_kinds
ORDER BY o.id;
CYPHER

section "DM EDGE AGGREGATES"
CYPHER <<'CYPHER'
MATCH ()-[r:Derived_metric_flags_gene]->()
RETURN 'flags' AS edge, count(r) AS total, count(r.value_flag) AS with_value,
       count(DISTINCT r.value_flag) AS distinct_values
UNION ALL
MATCH ()-[r:Derived_metric_classifies_gene]->()
RETURN 'classifies' AS edge, count(r) AS total, count(r.value_text) AS with_value,
       count(DISTINCT r.value_text) AS distinct_values
UNION ALL
MATCH ()-[r:Derived_metric_quantifies_gene]->()
RETURN 'quantifies' AS edge, count(r) AS total, count(r.value) AS with_value,
       count(r.rank_by_metric) AS distinct_values
ORDER BY edge;
CYPHER

section "NUMERIC DM RANK SUBSAMPLE: top-5 per DerivedMetric by rank_by_metric"
CYPHER <<'CYPHER'
MATCH (dm:DerivedMetric {rankable: 'true'})-[r:Derived_metric_quantifies_gene]->(g:Gene)
WHERE r.rank_by_metric <= 5
RETURN dm.id AS dm_id,
       r.rank_by_metric AS rank,
       g.locus_tag AS locus_tag,
       r.value AS value,
       r.metric_percentile AS pct,
       r.metric_bucket AS bucket
ORDER BY dm_id, rank, locus_tag;
CYPHER

section "GENE DM ROUTING AGGREGATES"
CYPHER <<'CYPHER'
MATCH (g:Gene)
RETURN count(g) AS total_genes,
       sum(g.numeric_metric_count) AS sum_numeric_count,
       sum(g.classifier_flag_count) AS sum_flag_count,
       sum(g.classifier_label_count) AS sum_label_count,
       count(CASE WHEN size(g.compartments_observed) > 0 THEN 1 END) AS genes_with_compartment;
CYPHER

printf '\n\n======== END ========\n'
