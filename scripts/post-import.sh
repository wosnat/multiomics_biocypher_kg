#!/bin/bash
set -euo pipefail

echo "=== Post-process: Starting Neo4j ==="
neo4j start
sleep 15

echo "=== Post-process: Cyanorak homolog edges ==="
cypher-shell <<'CYPHER'
MATCH (g1:Gene)-[:Gene_in_cyanorak_cluster]->(c:Cyanorak_cluster)<-[:Gene_in_cyanorak_cluster]-(g2:Gene)
WHERE id(g1) <> id(g2)
MATCH (g1)-[:Gene_belongs_to_organism]->(o1:OrganismTaxon)
MATCH (g2)-[:Gene_belongs_to_organism]->(o2:OrganismTaxon)
WITH g1, g2, c, o1, o2,
  CASE
    WHEN o1.id = o2.id
      THEN "same strain: " + o1.strain_name
    WHEN o1.clade IS NOT NULL AND o1.clade <> '' AND o1.clade = o2.clade
      THEN "same clade: " + o1.clade
    WHEN o1.species IS NOT NULL AND o1.species <> '' AND o1.species = o2.species
      THEN "same species: " + o1.species
    WHEN o1.genus IS NOT NULL AND o1.genus <> '' AND o1.genus = o2.genus
      THEN "same genus: " + o1.genus
    WHEN o1.order IS NOT NULL AND o1.order <> '' AND o1.order = o2.order
      THEN "same order: " + o1.order
    ELSE "cross order"
  END AS distance
MERGE (g1)-[:Gene_is_homolog_of_gene {source: "cyanorak_cluster", cluster_id: c.cluster_number, distance: distance}]->(g2);
CYPHER

echo "=== Post-process: Alteromonadaceae eggNOG homolog edges ==="
cypher-shell <<'CYPHER'
MATCH (g1:Gene)-[:Gene_belongs_to_organism]->(o1:OrganismTaxon)
MATCH (g2:Gene)-[:Gene_belongs_to_organism]->(o2:OrganismTaxon)
WHERE g1.alteromonadaceae_og IS NOT NULL
  AND g1.alteromonadaceae_og = g2.alteromonadaceae_og
  AND id(g1) <> id(g2)
  AND o1.genus = "Alteromonas"
  AND o2.genus = "Alteromonas"
WITH g1, g2, g1.alteromonadaceae_og AS og_id, o1, o2,
  CASE
    WHEN o1.id = o2.id THEN "same strain: " + o1.strain_name
    ELSE "same family: Alteromonadaceae"
  END AS distance
MERGE (g1)-[:Gene_is_homolog_of_gene {source: "eggnog_alteromonadaceae_og", cluster_id: og_id, distance: distance}]->(g2);
CYPHER

echo "=== Post-process: Cross-phylum COG homolog edges ==="
cypher-shell <<'CYPHER'
MATCH (g1:Gene)-[:Gene_belongs_to_organism]->(o1:OrganismTaxon)
MATCH (g2:Gene)-[:Gene_belongs_to_organism]->(o2:OrganismTaxon)
WHERE g1.bacteria_cog_og IS NOT NULL
  AND g1.bacteria_cog_og = g2.bacteria_cog_og
  AND id(g1) <> id(g2)
  AND (
    (o1.genus = "Alteromonas" AND o2.genus <> "Alteromonas")
    OR (o2.genus = "Alteromonas" AND o1.genus <> "Alteromonas")
  )
MERGE (g1)-[:Gene_is_homolog_of_gene {source: "eggnog_bacteria_cog_og", cluster_id: g1.bacteria_cog_og, distance: "cross phylum"}]->(g2);
CYPHER

echo "=== Post-process: Homolog edge counts ==="
cypher-shell "MATCH ()-[r:Gene_is_homolog_of_gene]->() RETURN r.source AS source, count(r) AS edges ORDER BY edges DESC;"

echo "=== Post-process: Propagate condition expression to homologs (no distance filter, batched) ==="
cypher-shell <<'CYPHER'
MATCH (source)-[e:Condition_changes_expression_of]->(geneA:Gene)-[h:Gene_is_homolog_of_gene]->(geneB:Gene)
WHERE NOT (source)-[:Condition_changes_expression_of_ortholog]->(geneB)
CALL {
  WITH source, e, geneA, h, geneB
  CREATE (source)-[:Condition_changes_expression_of_ortholog {
    expression_direction: e.expression_direction,
    control_condition: e.control_condition,
    experimental_context: e.experimental_context,
    time_point: e.time_point,
    log2_fold_change: e.log2_fold_change,
    adjusted_p_value: e.adjusted_p_value,
    significant: e.significant,
    publications: e.publications,
    omics_type: e.omics_type,
    organism_strain: e.organism_strain,
    treatment_condition: e.treatment_condition,
    statistical_test: e.statistical_test,
    original_gene: geneA.locus_tag,
    homology_source: 'cyanorak_cluster',
    homology_cluster_id: h.cluster_id,
    distance: h.distance
  }]->(geneB)
} IN TRANSACTIONS OF 5000 ROWS;
CYPHER

echo "=== Post-process: Propagate coculture expression to homologs (cross-phylum filter, batched) ==="
cypher-shell <<'CYPHER'
MATCH (source)-[e:Coculture_changes_expression_of]->(geneA:Gene)-[h:Gene_is_homolog_of_gene]->(geneB:Gene)
WHERE h.distance <> "cross phylum"
  AND NOT (source)-[:Coculture_changes_expression_of_ortholog]->(geneB)
CALL {
  WITH source, e, geneA, h, geneB
  CREATE (source)-[:Coculture_changes_expression_of_ortholog {
    expression_direction: e.expression_direction,
    control_condition: e.control_condition,
    experimental_context: e.experimental_context,
    time_point: e.time_point,
    log2_fold_change: e.log2_fold_change,
    adjusted_p_value: e.adjusted_p_value,
    significant: e.significant,
    publications: e.publications,
    omics_type: e.omics_type,
    organism_strain: e.organism_strain,
    treatment_condition: e.treatment_condition,
    statistical_test: e.statistical_test,
    original_gene: geneA.locus_tag,
    homology_source: 'cyanorak_cluster',
    homology_cluster_id: h.cluster_id,
    distance: h.distance
  }]->(geneB)
} IN TRANSACTIONS OF 5000 ROWS;
CYPHER

echo "=== Post-process: Final edge counts ==="
cypher-shell "MATCH ()-[r:Gene_is_homolog_of_gene]->() RETURN count(r) AS homolog_edges;"
cypher-shell "MATCH ()-[r:Condition_changes_expression_of_ortholog|Coculture_changes_expression_of_ortholog]->() RETURN count(r) AS total_ortholog_edges;"

echo "=== Post-process complete ==="
neo4j stop
