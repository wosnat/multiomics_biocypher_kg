// Post-import Cypher commands
// These run after the knowledge graph is imported into Neo4j

// Example: Create indexes for better query performance
// CREATE INDEX gene_name_idx IF NOT EXISTS FOR (g:Gene) ON (g.name);
// CREATE INDEX protein_name_idx IF NOT EXISTS FOR (p:Protein) ON (p.name);

// Example: Add computed properties or relationships
// MATCH (g:Gene)-[:ENCODES]->(p:Protein)
// SET g.has_protein = true;

// Add your custom Cypher commands below:

