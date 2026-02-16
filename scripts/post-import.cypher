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
MATCH (g1:gene)-[:gene_in_cyanorak_cluster]->(c:cyanorak_cluster)<-[:gene_in_cyanorak_cluster]-(g2:gene)
WHERE elementId(g1) <> elementId(g2)
MERGE (g1)-[:gene_is_homolog_of_gene {source: "cyanorak_cluster", cluster_id: c.cluster_number}]->(g2);
