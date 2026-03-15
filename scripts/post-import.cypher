// Post-import Cypher commands
// These run after the knowledge graph is imported into Neo4j

// Example: Create indexes for better query performance
// CREATE INDEX gene_name_idx IF NOT EXISTS FOR (g:Gene) ON (g.name);
// CREATE INDEX protein_name_idx IF NOT EXISTS FOR (p:Protein) ON (p.name);

// Example: Add computed properties or relationships
// MATCH (g:Gene)-[:ENCODES]->(p:Protein)
// SET g.has_protein = true;

// Add your custom Cypher commands below:

// Scalar indexes for get_gene exact lookup
CREATE INDEX gene_locus_tag_idx IF NOT EXISTS FOR (g:Gene) ON (g.locus_tag);
CREATE INDEX gene_name_idx IF NOT EXISTS FOR (g:Gene) ON (g.gene_name);
CREATE INDEX gene_organism_strain_idx IF NOT EXISTS FOR (g:Gene) ON (g.organism_strain);

// Full-text index for find_gene free-text search
CREATE FULLTEXT INDEX geneFullText IF NOT EXISTS FOR (n:Gene) ON EACH [
  n.gene_summary, n.gene_synonyms,
  n.alternate_functional_descriptions, n.pfam_names];

// Create gene_is_homolog_of_gene edges between genes sharing a Cyanorak cluster.
// Edges are created in BOTH directions (A→B and B→A) so LLM agents can
// traverse from any gene to all its homologs without directionless patterns.
// The `distance` property captures phylogenetic proximity between the two organisms.
MATCH (g1:Gene)-[:Gene_in_cyanorak_cluster]->(c:Cyanorak_cluster)<-[:Gene_in_cyanorak_cluster]-(g2:Gene)
WHERE elementId(g1) <> elementId(g2)
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

// Create Gene_is_homolog_of_gene edges between Alteromonas genes sharing an
// Alteromonadaceae-level eggNOG OG (within-Alteromonas orthologs).
// distance = "same family: Alteromonadaceae" for cross-strain pairs.
MATCH (g1:Gene)-[:Gene_belongs_to_organism]->(o1:OrganismTaxon)
MATCH (g2:Gene)-[:Gene_belongs_to_organism]->(o2:OrganismTaxon)
WHERE g1.alteromonadaceae_og IS NOT NULL
  AND g1.alteromonadaceae_og = g2.alteromonadaceae_og
  AND elementId(g1) <> elementId(g2)
  AND o1.genus = "Alteromonas"
  AND o2.genus = "Alteromonas"
WITH g1, g2, g1.alteromonadaceae_og AS og_id, o1, o2,
  CASE
    WHEN o1.id = o2.id THEN "same strain: " + o1.strain_name
    ELSE "same family: Alteromonadaceae"
  END AS distance
MERGE (g1)-[:Gene_is_homolog_of_gene {source: "eggnog_alteromonadaceae_og", cluster_id: og_id, distance: distance}]->(g2);

// Create Gene_is_homolog_of_gene edges between Alteromonas and Pro/Syn genes
// sharing a Bacteria-level COG OG (cross-phylum orthologs).
// Only creates edges where exactly one gene is Alteromonas to avoid duplicating
// the Cyanorak-based Pro/Syn<->Pro/Syn homolog edges above.
MATCH (g1:Gene)-[:Gene_belongs_to_organism]->(o1:OrganismTaxon)
MATCH (g2:Gene)-[:Gene_belongs_to_organism]->(o2:OrganismTaxon)
WHERE g1.bacteria_cog_og IS NOT NULL
  AND g1.bacteria_cog_og = g2.bacteria_cog_og
  AND elementId(g1) <> elementId(g2)
  AND (
    (o1.genus = "Alteromonas" AND o2.genus <> "Alteromonas")
    OR (o2.genus = "Alteromonas" AND o1.genus <> "Alteromonas")
  )
MERGE (g1)-[:Gene_is_homolog_of_gene {source: "eggnog_bacteria_cog_og", cluster_id: g1.bacteria_cog_og, distance: "cross phylum"}]->(g2);

// Create Condition_changes_expression_of_ortholog edges.
// If an EnvironmentalCondition X changes expression of gene A, and gene A is homolog of gene B,
// then X also changes_expression_of_ortholog gene B.
// No distance filter — condition experiments propagate to all homologs including cross-phylum.
// Uses CREATE (not MERGE) so each expression edge propagates independently —
// multiple publications/conditions on the same gene each produce a separate ortholog edge.
MATCH (source)-[e:Condition_changes_expression_of]->(geneA:Gene)-[h:Gene_is_homolog_of_gene]->(geneB:Gene)
WHERE NOT (source)-[:Condition_changes_expression_of_ortholog]->(geneB)
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
  analysis_name: e.analysis_name,
  original_gene: geneA.locus_tag,
  homology_source: 'cyanorak_cluster',
  homology_cluster_id: h.cluster_id,
  distance: h.distance
}]->(geneB);

// Create Coculture_changes_expression_of_ortholog edges.
// If an OrganismTaxon X (coculture partner) changes expression of gene A, and gene A is homolog of gene B,
// then X also changes_expression_of_ortholog gene B.
// Cross-phylum homologs are excluded (WHERE h.distance <> 'cross phylum') because
// coculture effects are expected to be organism-group-specific.
MATCH (source)-[e:Coculture_changes_expression_of]->(geneA:Gene)-[h:Gene_is_homolog_of_gene]->(geneB:Gene)
WHERE h.distance <> 'cross phylum'
  AND NOT (source)-[:Coculture_changes_expression_of_ortholog]->(geneB)
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
  analysis_name: e.analysis_name,
  original_gene: geneA.locus_tag,
  homology_source: 'cyanorak_cluster',
  homology_cluster_id: h.cluster_id,
  distance: h.distance
}]->(geneB);
