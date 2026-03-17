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
  n.gene_summary, n.gene_synonyms,
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

echo "=== Post-process complete ==="
neo4j stop
