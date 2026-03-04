---
name: cypher-queries
description: Run Cypher queries against the Neo4j knowledge graph. Provides ready-to-use query templates for gene expression, coculture effects, environmental stress, time series, and more. Use when the user wants to query the knowledge graph.
argument-hint: [query-description]
user-invocable: true
allowed-tools: Read, Grep, Glob, Bash(cypher-shell *), Bash(docker exec *)
---

# Cypher Query Skill

Run Cypher queries against the Prochlorococcus-Alteromonas Neo4j knowledge graph.

## Connection

The Neo4j database runs via Docker. Connect using:

```bash
# Via cypher-shell inside the container
docker exec -i biocypher-neo4j cypher-shell -u neo4j -p neo4j
```

To run a query non-interactively:

```bash
echo "MATCH (n) RETURN labels(n), count(*) LIMIT 10;" | docker exec -i biocypher-neo4j cypher-shell -u neo4j -p neo4j
```

If the container name differs, check with `docker ps` first.

---

## Graph Schema

### Node Labels (BioCypher PascalCase)

| Label | Key Properties |
|---|---|
| `Gene` | `id` (NCBI GeneID), `locus_tag`, `locus_tag_ncbi`, `product`, `gene_name`, `function_description`, `cog_category[]`, `go_terms[]`, `kegg_ko[]`, `kegg_pathway[]`, `ec_numbers[]`, `pfam_ids[]`, `bacteria_cog_og`, `alteromonadaceae_og`, `annotation_quality` |
| `Protein` | `id` (UniProt accession), `locus_tag`, `gene_symbol`, `function_description`, `ec_numbers[]`, `go_biological_processes[]`, `go_molecular_functions[]`, `go_cellular_components[]`, `is_reviewed` |
| `OrganismTaxon` | `id` (GCF accession), `organism_name`, `strain_name`, `clade`, `genus`, `species`, `ncbi_taxon_id` |
| `EnvironmentalCondition` | `id`, `name`, `condition_type`, `nitrogen_level`, `phosphate_level`, `light_condition`, `light_intensity`, `oxygen_level` |
| `Publication` | `id` (DOI), `title`, `authors[]`, `journal`, `publication_year` |
| `Cyanorak_cluster` | `id`, `cluster_number` |
| `BiologicalProcess` | `id` (GO:XXXXXXX), `name` |
| `CellularComponent` | `id` (GO:XXXXXXX), `name` |
| `MolecularFunction` | `id` (GO:XXXXXXX), `name` |
| `KeggOrthologousGroup` | `id` (K#####), `name` |
| `KeggPathway` | `id` (map#####), `name` |
| `KeggSubcategory` | `id`, `name` |
| `KeggCategory` | `id`, `name` |
| `EcNumber` | `id` (EC x.x.x.x), `name`, `catalytic_activity[]` |
| `CogFunctionalCategory` | `id`, `code`, `name` |
| `CyanorakRole` | `id`, `code`, `description` |
| `TigrRole` | `id`, `code`, `description` |

### Edge Labels

| Relationship | Source → Target | Key Properties |
|---|---|---|
| `Affects_expression_of` | OrganismTaxon/EnvironmentalCondition → Gene | `log2_fold_change`, `adjusted_p_value`, `expression_direction` (up/down), `control_condition`, `experimental_context`, `time_point`, `publications[]`, `significant` |
| `Affects_expression_of_homolog` | OrganismTaxon/EnvironmentalCondition → Gene | same + `original_gene`, `homology_source`, `homology_cluster_id`, `distance` |
| `Gene_is_homolog_of_gene` | Gene → Gene | `source` (cyanorak_cluster / eggnog_alteromonadaceae_og / eggnog_bacteria_cog_og), `cluster_id`, `distance` |
| `Gene_belongs_to_organism` | Gene → OrganismTaxon | — |
| `Protein_belongs_to_organism` | Protein → OrganismTaxon | — |
| `Gene_encodes_protein` | Gene → Protein | — |
| `Gene_in_cyanorak_cluster` | Gene → Cyanorak_cluster | — |
| `Gene_involved_in_biological_process` | Gene → BiologicalProcess | — |
| `Gene_located_in_cellular_component` | Gene → CellularComponent | — |
| `Gene_enables_molecular_function` | Gene → MolecularFunction | — |
| `Gene_has_kegg_ko` | Gene → KeggOrthologousGroup | — |
| `Ko_in_kegg_pathway` | KeggOrthologousGroup → KeggPathway | — |
| `Kegg_pathway_in_kegg_subcategory` | KeggPathway → KeggSubcategory | — |
| `Kegg_subcategory_in_kegg_category` | KeggSubcategory → KeggCategory | — |
| `Gene_catalyzes_ec_number` | Gene → EcNumber | — |
| `Ec_number_is_a_ec_number` | EcNumber → EcNumber | — |
| `Gene_in_cog_category` | Gene → CogFunctionalCategory | — |
| `Gene_has_cyanorak_role` | Gene → CyanorakRole | — |
| `Cyanorak_role_is_a_cyanorak_role` | CyanorakRole → CyanorakRole | — |
| `Gene_has_tigr_role` | Gene → TigrRole | — |
| `Protein_involved_in_biological_process` | Protein → BiologicalProcess | `evidence_code`, `reference` |
| `Protein_enables_molecular_function` | Protein → MolecularFunction | `evidence_code`, `reference` |
| `Protein_located_in_cellular_component` | Protein → CellularComponent | `evidence_code`, `reference` |
| `Biological_process_is_a_biological_process` | BiologicalProcess → BiologicalProcess | — |
| `Molecular_function_is_a_molecular_function` | MolecularFunction → MolecularFunction | — |
| `Cellular_component_is_a_cellular_component` | CellularComponent → CellularComponent | — |

### Homolog Edge Sources

`Gene_is_homolog_of_gene` edges come from three sources (filter by `source` property):

| `source` value | Coverage | `distance` values |
|---|---|---|
| `cyanorak_cluster` | Pro/Syn ↔ Pro/Syn | "same strain: X", "same clade: HLI", "same species: …", "same genus: …", "same order: …", "cross order" |
| `eggnog_alteromonadaceae_og` | Alteromonas ↔ Alteromonas | "same strain: X", "same family: Alteromonadaceae" |
| `eggnog_bacteria_cog_og` | Alteromonas ↔ Pro/Syn | "cross phylum" |

---

## Query Templates

### Gene Lookup

```cypher
// Find gene by locus tag
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
RETURN g.locus_tag, g.gene_name, g.product, g.function_description, o.strain_name

// Find gene by name (gene_name field)
MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
WHERE g.gene_name = 'dnaA'
RETURN g.locus_tag, g.product, o.strain_name

// Get full annotation for a gene (GO, KEGG, COG as properties)
MATCH (g:Gene {locus_tag: 'PMM0001'})
RETURN g.locus_tag, g.gene_name, g.product, g.function_description,
       g.cog_category, g.go_terms, g.kegg_ko, g.kegg_pathway,
       g.ec_numbers, g.pfam_ids, g.annotation_quality
```

### Gene with Linked Nodes (GO / KEGG / EC / COG)

```cypher
// Gene → GO biological process nodes
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_involved_in_biological_process]->(bp:BiologicalProcess)
RETURN g.locus_tag, bp.id, bp.name

// Gene → GO molecular function nodes
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_enables_molecular_function]->(mf:MolecularFunction)
RETURN mf.id, mf.name

// Gene → GO cellular component nodes
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_located_in_cellular_component]->(cc:CellularComponent)
RETURN cc.id, cc.name

// Gene → KEGG KO → Pathway → Subcategory → Category (full hierarchy)
MATCH (g:Gene {locus_tag: 'PMM0001'})
      -[:Gene_has_kegg_ko]->(ko:KeggOrthologousGroup)
      -[:Ko_in_kegg_pathway]->(pw:KeggPathway)
      -[:Kegg_pathway_in_kegg_subcategory]->(sc:KeggSubcategory)
      -[:Kegg_subcategory_in_kegg_category]->(cat:KeggCategory)
RETURN ko.id, ko.name, pw.name, sc.name, cat.name

// Gene → EC number
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_catalyzes_ec_number]->(ec:EcNumber)
RETURN ec.id, ec.name, ec.catalytic_activity

// Gene → COG functional category
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_in_cog_category]->(cog:CogFunctionalCategory)
RETURN cog.code, cog.name

// Gene → Cyanorak functional role (Pro/Syn only)
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_has_cyanorak_role]->(role:CyanorakRole)
RETURN role.code, role.description

// Gene → tIGR role (Pro/Syn only)
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_has_tigr_role]->(role:TigrRole)
RETURN role.code, role.description

// Combined: all annotations for a gene
MATCH (g:Gene {locus_tag: 'PMM0001'})
OPTIONAL MATCH (g)-[:Gene_involved_in_biological_process]->(bp:BiologicalProcess)
OPTIONAL MATCH (g)-[:Gene_enables_molecular_function]->(mf:MolecularFunction)
OPTIONAL MATCH (g)-[:Gene_has_kegg_ko]->(ko:KeggOrthologousGroup)
OPTIONAL MATCH (g)-[:Gene_catalyzes_ec_number]->(ec:EcNumber)
OPTIONAL MATCH (g)-[:Gene_in_cog_category]->(cog:CogFunctionalCategory)
RETURN g.locus_tag, g.gene_name, g.product,
       collect(DISTINCT bp.name) AS biological_processes,
       collect(DISTINCT mf.name) AS molecular_functions,
       collect(DISTINCT ko.name) AS kegg_kos,
       collect(DISTINCT ec.id) AS ec_numbers,
       collect(DISTINCT cog.name) AS cog_categories
```

### Find What Affects a Gene

```cypher
// All factors affecting a gene (coculture + stress, both directions)
MATCH (factor)-[r:Affects_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
RETURN labels(factor)[0] AS factor_type, factor.id, r.expression_direction,
       r.log2_fold_change, r.adjusted_p_value, r.control_condition,
       r.experimental_context, r.publications
ORDER BY abs(r.log2_fold_change) DESC

// Significantly up-regulated
MATCH (factor)-[r:Affects_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
WHERE r.expression_direction = 'up' AND (r.adjusted_p_value IS NULL OR r.adjusted_p_value < 0.05)
RETURN factor.id, r.log2_fold_change, r.adjusted_p_value, r.publications
ORDER BY r.log2_fold_change DESC
```

### Find Genes Affected by Organism Coculture

```cypher
// Genes upregulated in Prochlorococcus by Alteromonas coculture
MATCH (org:OrganismTaxon)-[r:Affects_expression_of {expression_direction: 'up'}]->(g:Gene)
WHERE org.genus = 'Alteromonas'
RETURN g.locus_tag, g.product, g.function_description,
       r.log2_fold_change, r.adjusted_p_value, r.experimental_context
ORDER BY r.log2_fold_change DESC
LIMIT 20

// Genes downregulated by coculture with a specific strain
MATCH (org:OrganismTaxon {strain_name: 'EZ55'})-[r:Affects_expression_of]->(g:Gene)
WHERE r.expression_direction = 'down'
RETURN g.locus_tag, g.product, r.log2_fold_change, r.control_condition
ORDER BY r.log2_fold_change ASC
LIMIT 20

// Count DE genes per coculture partner
MATCH (org:OrganismTaxon)-[r:Affects_expression_of]->(g:Gene)
WHERE r.significant = true
RETURN org.strain_name, r.expression_direction, count(g) AS gene_count
ORDER BY org.strain_name, r.expression_direction
```

### Find Genes Affected by Environmental Stress

```cypher
// Genes affected by nitrogen starvation
MATCH (env:EnvironmentalCondition {condition_type: 'nutrient_stress'})-[r:Affects_expression_of]->(g:Gene)
WHERE env.nitrogen_level = 'starved'
RETURN g.locus_tag, g.product, r.expression_direction, r.log2_fold_change
ORDER BY abs(r.log2_fold_change) DESC

// Genes affected by light/dark shift
MATCH (env:EnvironmentalCondition)-[r:Affects_expression_of]->(g:Gene)
WHERE env.light_condition = 'darkness'
RETURN g.locus_tag, g.product, r.expression_direction, r.log2_fold_change

// List all distinct environmental conditions in the graph
MATCH (env:EnvironmentalCondition)-[:Affects_expression_of]->(:Gene)
RETURN env.name, env.condition_type, env.nitrogen_level,
       env.light_condition, env.phosphate_level
ORDER BY env.condition_type, env.name

// DE genes from all stress experiments
MATCH (env:EnvironmentalCondition)-[r:Affects_expression_of]->(g:Gene)
WHERE r.significant = true
RETURN env.name, r.expression_direction, count(g) AS gene_count
ORDER BY gene_count DESC
```

### Compare Conditions

```cypher
// Genes affected by BOTH coculture AND nitrogen stress
MATCH (org:OrganismTaxon)-[r1:Affects_expression_of]->(g:Gene)
WHERE org.genus = 'Alteromonas'
MATCH (env:EnvironmentalCondition)-[r2:Affects_expression_of]->(g)
WHERE env.nitrogen_level = 'starved'
RETURN g.locus_tag, g.product,
       r1.expression_direction AS coculture_dir, r1.log2_fold_change AS coculture_fc,
       r2.expression_direction AS stress_dir, r2.log2_fold_change AS stress_fc

// Genes with OPPOSITE response to coculture vs stress (discordant)
MATCH (org:OrganismTaxon)-[r1:Affects_expression_of]->(g:Gene)
MATCH (env:EnvironmentalCondition)-[r2:Affects_expression_of]->(g)
WHERE r1.expression_direction <> r2.expression_direction
RETURN g.locus_tag, g.product, r1.expression_direction, r2.expression_direction

// Genes consistently up across multiple conditions
MATCH (factor)-[r:Affects_expression_of]->(g:Gene)
WHERE r.expression_direction = 'up' AND r.significant = true
WITH g, count(DISTINCT factor) AS up_count
WHERE up_count >= 3
RETURN g.locus_tag, g.product, g.function_description, up_count
ORDER BY up_count DESC
```

### Homologs

```cypher
// All homologs of a gene (all sources)
MATCH (g:Gene {locus_tag: 'PMM0001'})-[h:Gene_is_homolog_of_gene]->(hg:Gene)
      -[:Gene_belongs_to_organism]->(o:OrganismTaxon)
RETURN hg.locus_tag, hg.product, o.strain_name, h.source, h.cluster_id, h.distance
ORDER BY h.distance

// Only Cyanorak-cluster homologs (Pro/Syn ↔ Pro/Syn)
MATCH (g:Gene {locus_tag: 'PMM0001'})-[h:Gene_is_homolog_of_gene {source: 'cyanorak_cluster'}]->(hg:Gene)
      -[:Gene_belongs_to_organism]->(o:OrganismTaxon)
RETURN hg.locus_tag, o.strain_name, h.distance

// Cross-phylum homologs (Alteromonas ↔ Pro/Syn via bacteria COG OG)
MATCH (g:Gene {locus_tag: 'PMM0001'})-[h:Gene_is_homolog_of_gene {source: 'eggnog_bacteria_cog_og'}]->(hg:Gene)
      -[:Gene_belongs_to_organism]->(o:OrganismTaxon)
RETURN hg.locus_tag, hg.product, o.strain_name, h.cluster_id

// Within-Alteromonas homologs
MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon {strain_name: 'EZ55'})
MATCH (g)-[h:Gene_is_homolog_of_gene {source: 'eggnog_alteromonadaceae_og'}]->(hg:Gene)
      -[:Gene_belongs_to_organism]->(o2:OrganismTaxon)
RETURN g.locus_tag, hg.locus_tag, o2.strain_name, h.cluster_id
LIMIT 20

// Expression of homologs in all strains (via direct + propagated edges)
MATCH (g:Gene {locus_tag: 'PMM0001'})
OPTIONAL MATCH (factor1)-[r1:Affects_expression_of]->(g)
OPTIONAL MATCH (g)-[:Gene_is_homolog_of_gene]->(hg:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
OPTIONAL MATCH (factor2)-[r2:Affects_expression_of]->(hg)
RETURN g.locus_tag AS source_gene,
       r1.expression_direction AS direct_dir, r1.log2_fold_change AS direct_fc,
       hg.locus_tag AS homolog, o.strain_name,
       r2.expression_direction AS homolog_dir, r2.log2_fold_change AS homolog_fc
```

### Homolog Expression Propagation

```cypher
// Inferred homolog expression effects (post-import propagated edges)
MATCH (factor)-[r:Affects_expression_of_homolog]->(g:Gene)
WHERE r.distance CONTAINS 'same clade'
RETURN factor.id, g.locus_tag, r.expression_direction,
       r.log2_fold_change, r.original_gene, r.distance
LIMIT 20

// Homolog edges for a gene (inferred from its paralogs' expression)
MATCH (factor)-[r:Affects_expression_of_homolog]->(g:Gene {locus_tag: 'PMM0001'})
RETURN factor.id, r.original_gene, r.expression_direction,
       r.log2_fold_change, r.homology_source, r.distance
ORDER BY abs(r.log2_fold_change) DESC
```

### Functional Enrichment / Pathway Analysis

```cypher
// Genes upregulated by Alteromonas in each KEGG category
MATCH (org:OrganismTaxon)-[r:Affects_expression_of {expression_direction: 'up'}]->(g:Gene)
      -[:Gene_has_kegg_ko]->(ko:KeggOrthologousGroup)
      -[:Ko_in_kegg_pathway]->(pw:KeggPathway)
      -[:Kegg_pathway_in_kegg_subcategory]->(sc:KeggSubcategory)
      -[:Kegg_subcategory_in_kegg_category]->(cat:KeggCategory)
WHERE org.genus = 'Alteromonas'
RETURN cat.name, sc.name, pw.name, count(DISTINCT g) AS gene_count
ORDER BY gene_count DESC
LIMIT 20

// Genes downregulated by nitrogen stress in each COG category
MATCH (env:EnvironmentalCondition)-[r:Affects_expression_of {expression_direction: 'down'}]->(g:Gene)
      -[:Gene_in_cog_category]->(cog:CogFunctionalCategory)
WHERE env.nitrogen_level = 'starved'
RETURN cog.code, cog.name, count(DISTINCT g) AS gene_count
ORDER BY gene_count DESC

// Most common COG categories in all DE genes
MATCH (factor)-[r:Affects_expression_of]->(g:Gene)-[:Gene_in_cog_category]->(cog:CogFunctionalCategory)
WHERE r.significant = true
RETURN cog.code, cog.name, r.expression_direction, count(g) AS gene_count
ORDER BY gene_count DESC
LIMIT 20

// Genes in a specific KEGG pathway
MATCH (g:Gene)-[:Gene_has_kegg_ko]->(ko:KeggOrthologousGroup)
      -[:Ko_in_kegg_pathway]->(pw:KeggPathway {name: 'Photosynthesis'})
      -[:Gene_belongs_to_organism]->(o:OrganismTaxon)
RETURN g.locus_tag, g.product, ko.name, o.strain_name
ORDER BY o.strain_name

// Genes with a specific EC number
MATCH (g:Gene)-[:Gene_catalyzes_ec_number]->(ec:EcNumber)
      -[:Gene_belongs_to_organism]->(o:OrganismTaxon)
WHERE ec.id CONTAINS '2.7'  // kinases
RETURN g.locus_tag, g.product, ec.id, ec.name, o.strain_name
ORDER BY ec.id
```

### Gene + Protein Annotations

```cypher
// Gene with its protein (UniProt-reviewed)
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_encodes_protein]->(p:Protein)
RETURN g.locus_tag, p.id, p.function_description, p.is_reviewed,
       p.go_biological_processes, p.go_molecular_functions, p.keywords

// Genes with reviewed UniProt proteins
MATCH (g:Gene)-[:Gene_encodes_protein]->(p:Protein {is_reviewed: true})
      -[:Gene_belongs_to_organism]->(o:OrganismTaxon)
WHERE o.genus = 'Prochlorococcus'
RETURN g.locus_tag, g.product, p.id, p.function_description
LIMIT 20

// Genes whose protein has GO biological process annotation
MATCH (g:Gene)-[:Gene_encodes_protein]->(p:Protein)
      -[:Protein_involved_in_biological_process]->(bp:BiologicalProcess)
WHERE bp.name CONTAINS 'DNA repair'
RETURN g.locus_tag, g.product, p.id, bp.id, bp.name
```

### Strain / Organism Queries

```cypher
// All strains in the graph
MATCH (o:OrganismTaxon)
RETURN o.strain_name, o.organism_name, o.clade, o.genus, o.ncbi_taxon_id
ORDER BY o.genus, o.clade, o.strain_name

// Genes unique to one clade (not in other clades)
MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon {clade: 'HLI'})
WHERE NOT (g)-[:Gene_is_homolog_of_gene]->(:Gene)-[:Gene_belongs_to_organism]->
          (:OrganismTaxon {clade: 'HLII'})
RETURN g.locus_tag, g.product LIMIT 20

// Genes present in all Pro strains (Cyanorak cluster shared across all)
MATCH (c:Cyanorak_cluster)<-[:Gene_in_cyanorak_cluster]-(g:Gene)
      -[:Gene_belongs_to_organism]->(o:OrganismTaxon)
WHERE o.genus = 'Prochlorococcus'
WITH c, count(DISTINCT o) AS strain_count
WHERE strain_count >= 8
MATCH (c)<-[:Gene_in_cyanorak_cluster]-(g:Gene {locus_tag_ncbi: g.locus_tag_ncbi})
RETURN c.cluster_number, collect(g.locus_tag) AS genes, strain_count
ORDER BY strain_count DESC
LIMIT 10
```

### Significance Filtering

```cypher
// Significantly DE genes (p < 0.05, |logFC| > 1)
MATCH (factor)-[r:Affects_expression_of]->(g:Gene)
WHERE r.adjusted_p_value < 0.05 AND abs(r.log2_fold_change) > 1
RETURN labels(factor)[0] AS source_type, g.locus_tag, g.product,
       r.expression_direction, r.log2_fold_change, r.adjusted_p_value
ORDER BY r.adjusted_p_value ASC
LIMIT 20

// Top genes by fold change, regardless of significance
MATCH (factor)-[r:Affects_expression_of]->(g:Gene)
WHERE r.significant = true
RETURN g.locus_tag, g.product, r.expression_direction,
       max(abs(r.log2_fold_change)) AS max_fc
ORDER BY max_fc DESC
LIMIT 20
```

### Time Series Analysis

```cypher
// Gene expression across time points
MATCH (factor)-[r:Affects_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
WHERE r.time_point IS NOT NULL
RETURN factor.id, r.time_point, r.log2_fold_change, r.expression_direction
ORDER BY r.time_point

// Genes with consistent direction across all time points
MATCH (factor)-[r:Affects_expression_of]->(g:Gene)
WHERE r.time_point IS NOT NULL
WITH g, factor, collect(r.expression_direction) AS directions
WHERE size([d IN directions WHERE d = 'up']) = size(directions)
   OR size([d IN directions WHERE d = 'down']) = size(directions)
RETURN g.locus_tag, g.product, factor.id, directions
```

### Publication Traceability

```cypher
// Which studies reported expression changes for a gene
MATCH (factor)-[r:Affects_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
RETURN factor.id, r.log2_fold_change, r.expression_direction, r.publications

// All genes from a specific publication DOI
MATCH (factor)-[r:Affects_expression_of]->(g:Gene)
WHERE 'doi:10.1234/example' IN r.publications
RETURN g.locus_tag, g.product, factor.id, r.expression_direction

// Gene count per publication
MATCH (factor)-[r:Affects_expression_of]->(g:Gene)
UNWIND r.publications AS pub
WITH pub, count(DISTINCT g) AS gene_count
RETURN pub, gene_count
ORDER BY gene_count DESC
```

### Graph Statistics

```cypher
// Node counts by label
MATCH (n) RETURN labels(n)[0] AS label, count(*) AS count ORDER BY count DESC

// Edge counts by type
MATCH ()-[r]->() RETURN type(r) AS rel_type, count(*) AS count ORDER BY count DESC

// Genes with no expression data
MATCH (g:Gene)
WHERE NOT (g)<-[:Affects_expression_of]-()
RETURN count(g) AS genes_without_expression

// Annotation quality distribution
MATCH (g:Gene) RETURN g.annotation_quality, count(g) AS count ORDER BY g.annotation_quality
```

---

## Workflow

When the user invokes this skill (e.g., `/cypher-queries "genes affected by Alteromonas"`):

1. Understand the user's query intent
2. Select or adapt the most appropriate template above
3. Substitute actual values (locus tags, organism names, conditions) as needed
4. Run the query against Neo4j via `docker exec`
5. Present the results in a readable format
6. Offer follow-up queries if appropriate

**Important:** Always use exact BioCypher PascalCase labels from the schema table above. Common mistakes to avoid:
- Use `Gene` not `gene`, `OrganismTaxon` not `organism`
- Use `Affects_expression_of` not `affects_expression_of`
- Gene properties like `go_terms`, `kegg_ko`, `cog_category` are arrays (`str[]`) on the Gene node itself (denormalized for LLM use) AND available via linked nodes for graph traversal
