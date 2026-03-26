#!/bin/bash
set -euo pipefail

echo "=== Post-process: Starting Neo4j ==="
neo4j start
sleep 15

echo "=== Post-process: Create scalar indexes ==="
cypher-shell <<'CYPHER'
CREATE INDEX gene_locus_tag_idx IF NOT EXISTS FOR (g:Gene) ON (g.locus_tag);
CREATE INDEX gene_name_idx IF NOT EXISTS FOR (g:Gene) ON (g.gene_name);
CREATE INDEX gene_organism_strain_idx IF NOT EXISTS FOR (g:Gene) ON (g.organism_strain);
CYPHER

echo "=== Post-process: Create full-text indexes ==="
cypher-shell <<'CYPHER'
CREATE FULLTEXT INDEX geneFullText IF NOT EXISTS FOR (n:Gene) ON EACH [
  n.gene_summary, n.all_identifiers, n.gene_name_synonyms,
  n.alternate_functional_descriptions];
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
CYPHER

echo "=== Post-process: Create OrthologGroup indexes ==="
cypher-shell <<'CYPHER'
CREATE INDEX ortholog_group_id_idx IF NOT EXISTS FOR (og:OrthologGroup) ON (og.id);
CREATE INDEX ortholog_group_name_idx IF NOT EXISTS FOR (og:OrthologGroup) ON (og.name);
CREATE INDEX ortholog_group_level_idx IF NOT EXISTS FOR (og:OrthologGroup) ON (og.taxonomic_level);
CREATE INDEX ortholog_group_rank_idx IF NOT EXISTS FOR (og:OrthologGroup) ON (og.specificity_rank);
CYPHER

echo "=== Post-process: Create Pfam indexes ==="
cypher-shell <<'CYPHER'
CREATE INDEX pfam_name_idx IF NOT EXISTS FOR (p:Pfam) ON (p.name);
CREATE INDEX pfam_clan_name_idx IF NOT EXISTS FOR (c:PfamClan) ON (c.name);

CREATE FULLTEXT INDEX pfamFullText IF NOT EXISTS
  FOR (p:Pfam) ON EACH [p.name, p.short_name];
CREATE FULLTEXT INDEX pfamClanFullText IF NOT EXISTS
  FOR (c:PfamClan) ON EACH [c.name];
CYPHER

echo "=== Post-process: Create Publication indexes ==="
cypher-shell <<'CYPHER'
CREATE FULLTEXT INDEX publicationFullText IF NOT EXISTS
  FOR (p:Publication) ON EACH [p.title, p.abstract, p.description];
CYPHER

echo "=== Post-process: Create Experiment indexes ==="
cypher-shell <<'CYPHER'
CREATE INDEX experiment_id_idx IF NOT EXISTS FOR (e:Experiment) ON (e.id);
CREATE INDEX experiment_organism_idx IF NOT EXISTS FOR (e:Experiment) ON (e.organism_strain);
CREATE INDEX experiment_treatment_type_idx IF NOT EXISTS FOR (e:Experiment) ON (e.treatment_type);
CREATE INDEX experiment_omics_type_idx IF NOT EXISTS FOR (e:Experiment) ON (e.omics_type);

CREATE FULLTEXT INDEX experimentFullText IF NOT EXISTS
  FOR (e:Experiment) ON EACH [e.name, e.treatment, e.control, e.experimental_context, e.light_condition];
CYPHER

echo "=== Post-process: Compute Publication summary properties ==="
cypher-shell <<'CYPHER'
MATCH (p:Publication)
OPTIONAL MATCH (p)-[:Has_experiment]->(e:Experiment)
WITH p,
     count(e) AS ec,
     [x IN collect(DISTINCT e.treatment_type) WHERE x IS NOT NULL] AS tts,
     [x IN collect(DISTINCT e.omics_type) WHERE x IS NOT NULL] AS ots,
     [x IN collect(DISTINCT e.organism_strain) WHERE x IS NOT NULL] AS orgs,
     [x IN collect(DISTINCT e.coculture_partner) WHERE x IS NOT NULL AND x <> ''] AS coculture_orgs
SET p.experiment_count = ec,
    p.treatment_types = apoc.coll.sort(tts),
    p.omics_types = apoc.coll.sort(ots),
    p.organisms = apoc.coll.sort(apoc.coll.toSet(orgs + coculture_orgs))
CYPHER

echo "=== Post-process: Compute expression_status on edges ==="
cypher-shell <<'CYPHER'
MATCH ()-[r:Changes_expression_of]->()
SET r.expression_status = CASE
  WHEN r.significant = 'significant' AND r.expression_direction = 'up'   THEN 'significant_up'
  WHEN r.significant = 'significant' AND r.expression_direction = 'down' THEN 'significant_down'
  ELSE 'not_significant'
END
CYPHER

echo "=== Post-process: Compute Experiment summary properties ==="

echo "--- experiment stats defaults ---"
cypher-shell <<'CYPHER'
MATCH (e:Experiment)
SET e.gene_count = 0,
    e.significant_up_count = 0,
    e.significant_down_count = 0,
    e.time_point_count = 0,
    e.time_point_labels = [],
    e.time_point_orders = [],
    e.time_point_hours = [],
    e.time_point_totals = [],
    e.time_point_significant_up = [],
    e.time_point_significant_down = []
CYPHER

echo "--- experiment stats computation ---"
cypher-shell <<'CYPHER'
MATCH (e:Experiment)-[r:Changes_expression_of]->(g:Gene)
WITH e,
     COALESCE(r.time_point, '') AS tp,
     r.time_point_order AS tp_order,
     COALESCE(r.time_point_hours, -1.0) AS tp_hours,
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
    e.time_point_significant_down = tp_sig_down
CYPHER

echo "=== Post-process: Compute OrganismTaxon summary properties ==="

echo "--- gene_count ---"
cypher-shell <<'CYPHER'
MATCH (o:OrganismTaxon)
OPTIONAL MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o)
WITH o, count(g) AS gc
SET o.gene_count = gc
CYPHER

echo "--- publication_count, experiment_count, treatment_types, omics_types ---"
cypher-shell <<'CYPHER'
MATCH (o:OrganismTaxon)
OPTIONAL MATCH (p:Publication)
  WHERE ANY(org IN p.organisms WHERE org = o.preferred_name)
WITH o,
     count(DISTINCT p) AS pc,
     CASE WHEN count(p) > 0 THEN sum(p.experiment_count) ELSE 0 END AS ec,
     apoc.coll.toSet(reduce(s = [], t IN collect(p.treatment_types) | s + t)) AS tts,
     apoc.coll.toSet(reduce(s = [], t IN collect(p.omics_types) | s + t)) AS ots
SET o.publication_count = pc,
    o.experiment_count = ec,
    o.treatment_types = tts,
    o.omics_types = ots
CYPHER

echo "=== Post-process: Compute Gene routing signals ==="

echo "--- annotation_types ---"
cypher-shell <<'CYPHER'
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
    CASE WHEN EXISTS { (g)-[:Gene_catalyzes_ec_number]->() } THEN ['ec'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_cyanorak_role]->() } THEN ['cyanorak_role'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_tigr_role]->() } THEN ['tigr_role'] ELSE [] END
} IN TRANSACTIONS OF 1000 ROWS;
CYPHER

echo "--- expression_edge_count + significant_up/down_count ---"
cypher-shell <<'CYPHER'
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
CYPHER

echo "--- rank_by_effect ---"
cypher-shell <<'CYPHER'
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
CYPHER

echo "--- closest_ortholog_group_size + closest_ortholog_genera ---"
cypher-shell <<'CYPHER'
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[:Gene_in_ortholog_group]->(og:OrthologGroup)
  WITH g, og ORDER BY og.specificity_rank ASC LIMIT 1
  SET g.closest_ortholog_group_size = og.member_count,
      g.closest_ortholog_genera = og.genera
} IN TRANSACTIONS OF 1000 ROWS;
CYPHER

echo "=== Post-process complete ==="
neo4j stop
