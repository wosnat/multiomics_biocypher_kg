// Post-import Cypher commands
// These run after the knowledge graph is imported into Neo4j
// Indexes only — homolog edges and expression propagation have been replaced
// by OrthologGroup nodes and Gene_in_ortholog_group edges (query-time joins).

// Scalar indexes for get_gene exact lookup
CREATE INDEX gene_locus_tag_idx IF NOT EXISTS FOR (g:Gene) ON (g.locus_tag);
CREATE INDEX gene_name_idx IF NOT EXISTS FOR (g:Gene) ON (g.gene_name);
CREATE INDEX gene_organism_strain_idx IF NOT EXISTS FOR (g:Gene) ON (g.organism_strain);

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

// Pfam domain and clan indexes
CREATE INDEX pfam_name_idx IF NOT EXISTS FOR (p:Pfam) ON (p.name);
CREATE INDEX pfam_clan_name_idx IF NOT EXISTS FOR (c:PfamClan) ON (c.name);

CREATE FULLTEXT INDEX pfamFullText IF NOT EXISTS
  FOR (p:Pfam) ON EACH [p.name, p.short_name];
CREATE FULLTEXT INDEX pfamClanFullText IF NOT EXISTS
  FOR (c:PfamClan) ON EACH [c.name];

// Experiment indexes
CREATE INDEX experiment_id_idx IF NOT EXISTS FOR (e:Experiment) ON (e.id);
CREATE INDEX experiment_organism_idx IF NOT EXISTS FOR (e:Experiment) ON (e.organism_strain);
CREATE INDEX experiment_treatment_type_idx IF NOT EXISTS FOR (e:Experiment) ON (e.treatment_type);
CREATE INDEX experiment_omics_type_idx IF NOT EXISTS FOR (e:Experiment) ON (e.omics_type);

CREATE FULLTEXT INDEX experimentFullText IF NOT EXISTS
  FOR (e:Experiment) ON EACH [e.name, e.treatment, e.control, e.experimental_context, e.light_condition];

// Publication fulltext index
CREATE FULLTEXT INDEX publicationFullText IF NOT EXISTS
  FOR (p:Publication) ON EACH [p.title, p.abstract, p.description];

// -----------------------------------------------------------------------
// Publication summary properties (pre-computed for list_publications)
// -----------------------------------------------------------------------

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
    p.organisms = apoc.coll.sort(apoc.coll.toSet(orgs + coculture_orgs));

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
    CASE WHEN EXISTS { (g)-[:Gene_catalyzes_ec_number]->() } THEN ['ec'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_cyanorak_role]->() } THEN ['cyanorak_role'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_tigr_role]->() } THEN ['tigr_role'] ELSE [] END
} IN TRANSACTIONS OF 1000 ROWS;

// expression_edge_count + significant_expression_count
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)<-[e:Changes_expression_of]-()
  WITH g, count(e) AS total,
       sum(CASE WHEN e.significant = 'significant' THEN 1 ELSE 0 END) AS sig
  SET g.expression_edge_count = total,
      g.significant_expression_count = sig
} IN TRANSACTIONS OF 500 ROWS;

// rank_by_effect: within each experiment + timepoint, rank by |log2FC| descending
MATCH (e:Experiment)
WITH e
CALL {
  WITH e
  MATCH (e)-[r:Changes_expression_of]->(g:Gene)
  WITH r.time_point_order AS tp, r, abs(r.log2_fold_change) AS abs_fc
  ORDER BY tp, abs_fc DESC
  WITH tp, collect(r) AS edges
  UNWIND range(0, size(edges)-1) AS i
  SET (edges[i]).rank_by_effect = i + 1
} IN TRANSACTIONS OF 10 ROWS;

// closest_ortholog_group_size + closest_ortholog_genera
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[:Gene_in_ortholog_group]->(og:OrthologGroup)
  WITH g, og ORDER BY og.specificity_rank ASC LIMIT 1
  SET g.closest_ortholog_group_size = og.member_count,
      g.closest_ortholog_genera = og.genera
} IN TRANSACTIONS OF 1000 ROWS;
