# Ortholog Group Migration Guide

This document describes the migration from the old homolog/ortholog edge model to the new OrthologGroup node model. It is intended for anyone writing Cypher queries against the knowledge graph (MCP tools, Biochatter prompts, manual exploration).

## Summary of Removed Elements

| Element | Type | Approximate count |
|---|---|---|
| `Cyanorak_cluster` | Node | 5,619 |
| `Gene_in_cyanorak_cluster` | Edge | 20,657 |
| `Gene_is_homolog_of_gene` | Edge | 365,840 |
| `Condition_changes_expression_of_ortholog` | Edge | 2,005,209 |
| `Coculture_changes_expression_of_ortholog` | Edge | 82,782 |
| Gene property `cluster_number` | Property | -- |
| Gene property `bacteria_cog_og` | Property | -- |
| Gene property `alteromonadaceae_og` | Property | -- |

## Summary of Added Elements

| Element | Type | Approximate count | Properties |
|---|---|---|---|
| `OrthologGroup` | Node | ~21,000 | `name` (raw OG ID, e.g. "CK_00000364", "COG0592@2"), `source` ("cyanorak" or "eggnog"), `taxonomic_level` ("curated", "Prochloraceae", "Synechococcus", "Alteromonadaceae", "Bacteria", "Cyanobacteria", "Gammaproteobacteria"), `taxon_id` (integer), `specificity_rank` (int: 0=curated, 1=family, 2=order, 3=domain/Bacteria) |
| `Gene_in_ortholog_group` | Edge | ~84,500 | (no properties) |

A gene may belong to 1-3 ortholog groups simultaneously:
- One Cyanorak curated cluster (Prochlorococcus and Synechococcus strains only)
- One eggNOG lowest-level OG (organism-group-specific, e.g. Prochloraceae or Alteromonadaceae)
- One eggNOG Bacteria-level COG

## Cypher Migration Table

| Use case | Old Cypher | New Cypher |
|---|---|---|
| Find homologs of a gene | `MATCH (g:Gene)-[:Gene_is_homolog_of_gene]-(h:Gene) WHERE g.locus_tag = $lt RETURN h` | `MATCH (g:Gene)-[:Gene_in_ortholog_group]->(og:OrthologGroup)<-[:Gene_in_ortholog_group]-(h:Gene) WHERE g.locus_tag = $lt AND g <> h RETURN h, og.source, og.taxonomic_level` |
| Find homologs at a specific level | `MATCH (g)-[:Gene_is_homolog_of_gene {source: "cyanorak"}]-(h)` | `MATCH (g)-[:Gene_in_ortholog_group]->(og:OrthologGroup {source: "cyanorak"})<-[:Gene_in_ortholog_group]-(h) WHERE g <> h RETURN h` |
| Ortholog expression propagation (condition) | `MATCH (c)-[e:Condition_changes_expression_of_ortholog]->(g:Gene) WHERE g.locus_tag = $lt RETURN e` | `MATCH (g:Gene)-[:Gene_in_ortholog_group]->(og:OrthologGroup)<-[:Gene_in_ortholog_group]-(h:Gene)<-[e:Condition_changes_expression_of]-(c) WHERE g.locus_tag = $lt AND g <> h RETURN e, h.locus_tag, og.source, og.taxonomic_level` |
| Coculture ortholog propagation | `MATCH (o)-[e:Coculture_changes_expression_of_ortholog]->(g:Gene) WHERE g.locus_tag = $lt RETURN e` | `MATCH (g:Gene)-[:Gene_in_ortholog_group]->(og:OrthologGroup)<-[:Gene_in_ortholog_group]-(h:Gene)<-[e:Coculture_changes_expression_of]-(o) WHERE g.locus_tag = $lt AND g <> h RETURN e, h.locus_tag, og.source, og.taxonomic_level` |
| Cyanorak cluster members | `MATCH (g:Gene)-[:Gene_in_cyanorak_cluster]->(c:Cyanorak_cluster) WHERE c.cluster_number = $ck RETURN g` | `MATCH (g:Gene)-[:Gene_in_ortholog_group]->(og:OrthologGroup {source: "cyanorak"}) WHERE og.name = $ck RETURN g` |
| Cross-phylum COG bridging | `MATCH (g)-[:Gene_is_homolog_of_gene {source: "bacteria_cog"}]-(h)` | `MATCH (g)-[:Gene_in_ortholog_group]->(og:OrthologGroup {source: "eggnog", taxonomic_level: "Bacteria"})<-[:Gene_in_ortholog_group]-(h) WHERE g <> h RETURN h` |

## Performance Notes

- The new model uses 2-hop traversals (`gene -> OrthologGroup <- gene`) instead of the old 1-hop direct homolog edges. This adds slightly higher query latency per traversal.
- However, the total stored edge count dropped dramatically: ~84,500 membership edges replace ~2,474,000 materialized edges (homolog + ortholog expression). This reduces graph size and import time significantly.
- Indexes are available on `OrthologGroup` nodes: `ortholog_group_id_idx` (id), `ortholog_group_name_idx` (name), `ortholog_group_level_idx` (taxonomic_level). Use these to filter by source or level efficiently.
- For expression propagation queries, always include `AND g <> h` to avoid self-joins, and filter on `og.taxonomic_level` or `og.source` when you only need a specific ortholog scope.

## Breaking Changes Checklist

The following elements no longer exist in the graph. Any query, tool, or prompt referencing them must be updated:

- **Removed edge types:**
  - `Gene_is_homolog_of_gene`
  - `Gene_in_cyanorak_cluster`
  - `Condition_changes_expression_of_ortholog`
  - `Coculture_changes_expression_of_ortholog`

- **Removed node type:**
  - `Cyanorak_cluster`

- **Removed gene properties:**
  - `cluster_number`
  - `bacteria_cog_og`
  - `alteromonadaceae_og`

- **Removed edge properties** (were on `Gene_is_homolog_of_gene`):
  - `distance`
  - `cluster_id`
  - `source`