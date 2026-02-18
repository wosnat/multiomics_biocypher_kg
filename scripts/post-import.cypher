// Post-import Cypher commands
// These run after the knowledge graph is imported into Neo4j

// Example: Create indexes for better query performance
// CREATE INDEX gene_name_idx IF NOT EXISTS FOR (g:Gene) ON (g.name);
// CREATE INDEX protein_name_idx IF NOT EXISTS FOR (p:Protein) ON (p.name);

// Example: Add computed properties or relationships
// MATCH (g:Gene)-[:ENCODES]->(p:Protein)
// SET g.has_protein = true;

// Add your custom Cypher commands below:

// Create gene_is_homolog_of_gene edges between genes sharing a Cyanorak cluster.
// Edges are created in BOTH directions (A→B and B→A) so LLM agents can
// traverse from any gene to all its homologs without directionless patterns.
MATCH (g1:Gene)-[:Gene_in_cyanorak_cluster]->(c:Cyanorak_cluster)<-[:Gene_in_cyanorak_cluster]-(g2:Gene)
WHERE id(g1) <> id(g2)
MERGE (g1)-[:Gene_is_homolog_of_gene {source: "cyanorak_cluster", cluster_id: c.cluster_number}]->(g2);

// Create affects_expression_of_homolog edges.
// If source X affects_expression_of gene A, and gene A is homolog of gene B,
// then X affects_expression_of_homolog gene B.
// Uses CREATE (not MERGE) so each expression edge propagates independently —
// multiple publications/conditions on the same gene each produce a separate homolog edge.
MATCH (source)-[e:Affects_expression_of]->(geneA:Gene)-[h:Gene_is_homolog_of_gene]->(geneB:Gene)
CREATE (source)-[:Affects_expression_of_homolog {
  expression_direction: e.expression_direction,
  control_condition: e.control_condition,
  experimental_context: e.experimental_context,
  time_point: e.time_point,
  log2_fold_change: e.log2_fold_change,
  adjusted_p_value: e.adjusted_p_value,
  significant: e.significant,
  publications: e.publications,
  original_gene: geneA.id,
  homology_source: h.source,
  homology_cluster_id: h.cluster_id
}]->(geneB);
