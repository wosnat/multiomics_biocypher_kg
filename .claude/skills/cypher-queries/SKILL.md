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

The Neo4j database is accessible at `localhost:7687` (Bolt, no auth). Connect using:

```bash
# Via cypher-shell directly (no auth)
cypher-shell -a bolt://localhost:7687
```

To run a query non-interactively:

```bash
echo "MATCH (n) RETURN labels(n), count(*) LIMIT 10;" | cypher-shell -a bolt://localhost:7687
```

If cypher-shell is not installed locally, fall back to Docker:

```bash
docker exec -i deploy cypher-shell -a bolt://localhost:7687
```

---

## Graph Schema

### Node Labels (BioCypher PascalCase)

| Label | Key Properties |
|---|---|
| `Gene` | `id` (ncbigene:locus_tag), `locus_tag`, `locus_tag_ncbi`, `locus_tag_cyanorak`, `product`, `gene_name`, `gene_synonyms[]`, `function_description`, `cog_category[]`, `go_terms[]`, `kegg_ko[]`, `kegg_pathway[]`, `ec_numbers[]`, `pfam_ids[]`, `pfam_names[]`, `bacteria_cog_og`, `alteromonadaceae_og`, `annotation_quality`, `old_locus_tags[]`, `alternative_locus_tags[]`, `protein_id`, `start`, `end`, `strand` |
| `Protein` | `id` (uniprot:accession), `locus_tag`, `gene_names[]`, `is_reviewed`, `annotation_score`, `sequence_length`, `molecular_mass`, `refseq_ids[]`, `proteome_ids[]`, `protein_synonyms[]`, `string_ids[]` |
| `OrganismTaxon` | `id` (insdc.gcf:accession or ncbitaxon:taxid), `organism_name`, `strain_name`, `clade`, `genus`, `species`, `ncbi_taxon_id`, `family`, `order`, `phylum`, `kingdom`, `superkingdom`, `lineage` |
| `EnvironmentalCondition` | `id`, `name`, `condition_type`, `description`, `light_condition`, `light_intensity`, `oxygen_level`, `nitrogen_source`, `co2_level`, `publications[]` |
| `Publication` | `id` (DOI), `title`, `authors[]`, `journal`, `publication_year`, `abstract`, `doi`, `organism`, `study_type`, `description` |
| `Cyanorak_cluster` | `id`, `cluster_number` |
| `BiologicalProcess` | `id` (go:NNNNNNN), `name` |
| `CellularComponent` | `id` (go:NNNNNNN), `name` |
| `MolecularFunction` | `id` (go:NNNNNNN), `name` |
| `KeggOrthologousGroup` | `id` (K#####), `name` |
| `KeggPathway` | `id` (map#####), `name` |
| `KeggSubcategory` | `id`, `name` |
| `KeggCategory` | `id`, `name` |
| `EcNumber` | `id` (EC x.x.x.x), `name`, `catalytic_activity[]` |
| `CogFunctionalCategory` | `id`, `code`, `name` |
| `CyanorakRole` | `id`, `code`, `description` |
| `TigrRole` | `id`, `code`, `description` |

### Additional Gene Properties (denormalized annotations)

Gene nodes carry many denormalized annotation fields beyond the key properties above:

| Property | Source | Description |
|---|---|---|
| `function_description_source` | provenance | Which source provided function_description |
| `gene_name_source` | provenance | Which source provided gene_name |
| `product_source` | provenance | Which source provided product |
| `alternate_functional_descriptions[]` | merged | Descriptions from non-primary sources |
| `product_cyanorak` | Cyanorak | Cyanorak-specific product name |
| `go_term_descriptions[]` | eggNOG | Human-readable GO term names |
| `kegg_ko_descriptions[]` | eggNOG | Human-readable KO names |
| `kegg_brite[]` | eggNOG | KEGG BRITE hierarchy entries |
| `kegg_module[]` | eggNOG | KEGG module IDs |
| `kegg_reaction[]` | eggNOG | KEGG reaction IDs |
| `pfam_descriptions[]` | eggNOG | Human-readable Pfam domain names |
| `eggnog_ogs[]` | eggNOG | Full eggNOG orthologous group IDs |
| `eggnog_og_descriptions[]` | eggNOG | Human-readable OG descriptions |
| `seed_ortholog` | eggNOG | Best seed ortholog hit |
| `seed_ortholog_evalue` | eggNOG | E-value of seed ortholog hit |
| `max_annot_lvl` | eggNOG | Deepest taxonomic level with annotation |
| `cyanorak_Role` | Cyanorak | Functional role code (Pro/Syn only) |
| `cyanorak_Role_description` | Cyanorak | Role description text |
| `tIGR_Role` | Cyanorak | TIGR role code (Pro/Syn only) |
| `tIGR_Role_description` | Cyanorak | TIGR role description text |
| `signal_peptide` | eggNOG | Signal peptide prediction |
| `transmembrane_regions` | eggNOG | Transmembrane topology |
| `protein_family` | Cyanorak | Protein family classification |
| `transporter_classification` | Cyanorak | TC number |
| `bigg_reaction[]` | eggNOG | BiGG reaction IDs |
| `cazy_ids[]` | eggNOG | CAZy family IDs |

### Edge Labels

| Relationship | Source → Target | Key Properties |
|---|---|---|
| `Affects_expression_of` | OrganismTaxon/EnvironmentalCondition → Gene | `log2_fold_change`, `adjusted_p_value`, `expression_direction` (up/down), `control_condition`, `experimental_context`, `time_point`, `publications[]`, `significant` |
| `Affects_expression_of_homolog` | OrganismTaxon/EnvironmentalCondition → Gene | same + `original_gene`, `homology_source`, `homology_cluster_id`, `distance` |
| `Gene_is_homolog_of_gene` | Gene ↔ Gene | `source` (cyanorak_cluster / eggnog_alteromonadaceae_og / eggnog_bacteria_cog_og), `cluster_id`, `distance` |
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
| `Biological_process_is_a_biological_process` | BiologicalProcess → BiologicalProcess | — |
| `Biological_process_part_of_biological_process` | BiologicalProcess → BiologicalProcess | — |
| `Biological_process_positively_regulates_biological_process` | BiologicalProcess → BiologicalProcess | — |
| `Biological_process_negatively_regulates_biological_process` | BiologicalProcess → BiologicalProcess | — |
| `Molecular_function_is_a_molecular_function` | MolecularFunction → MolecularFunction | — |
| `Molecular_function_part_of_molecular_function` | MolecularFunction → MolecularFunction | — |
| `Cellular_component_is_a_cellular_component` | CellularComponent → CellularComponent | — |
| `Cellular_component_part_of_cellular_component` | CellularComponent → CellularComponent | — |

### Homolog Edge Sources

`Gene_is_homolog_of_gene` edges come from three sources (filter by `source` property):

| `source` value | Coverage | `distance` values |
|---|---|---|
| `cyanorak_cluster` | Pro/Syn ↔ Pro/Syn | "same strain: X", "same clade: HLI", "same species: …", "same genus: …", "same order: …", "cross order" |
| `eggnog_alteromonadaceae_og` | Alteromonas ↔ Alteromonas | "same strain: X", "same family: Alteromonadaceae" |
| `eggnog_bacteria_cog_og` | Alteromonas ↔ Pro/Syn | "cross phylum" |

### Environmental Condition Types

| `condition_type` | Count | Description |
|---|---|---|
| `growth_medium` | 20 | Baseline growth conditions (control) |
| `light_stress` | 9 | Light/dark shifts, high-light stress |
| `gas_shock` | 6 | CO₂/O₂ perturbations |
| `nutrient_stress` | 5 | Nitrogen, iron, phosphate limitation |
| `growth_state` | 4 | Stationary vs exponential phase |
| `salt_stress` | 2 | Salinity changes |
| `pco2` | 2 | pCO₂ level experiments |
| `coculture` | 2 | Coculture conditions (source is EnvironmentalCondition) |
| `plastic_leachate_stress` | 2 | Plastic leachate exposure |
| `iron_stress` | 2 | Iron limitation/addition |
| `dark_tolerance` | 1 | Extended darkness survival |
| `viral_lysis_products` | 1 | Viral lysate exposure |

Note: Nutrient-specific properties (`nitrogen_source`, `co2_level`, `oxygen_level`) are only populated on a subset of conditions. Use `condition_type` and `description` for filtering — `description` contains the most detail.

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

// Search gene by product description (partial match)
MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
WHERE g.product CONTAINS 'photosystem'
RETURN g.locus_tag, g.gene_name, g.product, o.strain_name
ORDER BY o.strain_name

// Get full annotation for a gene (denormalized properties)
MATCH (g:Gene {locus_tag: 'PMM0001'})
RETURN g.locus_tag, g.gene_name, g.product, g.function_description,
       g.cog_category, g.go_terms, g.kegg_ko, g.kegg_pathway,
       g.ec_numbers, g.pfam_ids, g.pfam_names, g.annotation_quality,
       g.bacteria_cog_og, g.alteromonadaceae_og
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

// Gene → Cyanorak functional role hierarchy (Pro/Syn only)
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_has_cyanorak_role]->(role:CyanorakRole)
OPTIONAL MATCH (role)-[:Cyanorak_role_is_a_cyanorak_role]->(parent:CyanorakRole)
RETURN role.code, role.description, parent.code AS parent_code, parent.description AS parent_desc

// Gene → TIGR role (Pro/Syn only)
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_has_tigr_role]->(role:TigrRole)
RETURN role.code, role.description

// Combined: all linked annotations for a gene
MATCH (g:Gene {locus_tag: 'PMM0001'})
OPTIONAL MATCH (g)-[:Gene_involved_in_biological_process]->(bp:BiologicalProcess)
OPTIONAL MATCH (g)-[:Gene_enables_molecular_function]->(mf:MolecularFunction)
OPTIONAL MATCH (g)-[:Gene_has_kegg_ko]->(ko:KeggOrthologousGroup)
OPTIONAL MATCH (g)-[:Gene_catalyzes_ec_number]->(ec:EcNumber)
OPTIONAL MATCH (g)-[:Gene_in_cog_category]->(cog:CogFunctionalCategory)
OPTIONAL MATCH (g)-[:Gene_has_cyanorak_role]->(crole:CyanorakRole)
RETURN g.locus_tag, g.gene_name, g.product,
       collect(DISTINCT bp.name) AS biological_processes,
       collect(DISTINCT mf.name) AS molecular_functions,
       collect(DISTINCT ko.name) AS kegg_kos,
       collect(DISTINCT ec.id) AS ec_numbers,
       collect(DISTINCT cog.name) AS cog_categories,
       collect(DISTINCT crole.description) AS cyanorak_roles
```

### Find What Affects a Gene

```cypher
// All factors affecting a gene (coculture + stress, both directions)
MATCH (factor)-[r:Affects_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
RETURN labels(factor)[0] AS factor_type, factor.id, factor.name,
       r.expression_direction, r.log2_fold_change, r.adjusted_p_value,
       r.control_condition, r.experimental_context, r.publications
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
WHERE org.organism_name = 'Alteromonas'
RETURN g.locus_tag, g.product, g.function_description,
       r.log2_fold_change, r.adjusted_p_value, r.experimental_context
ORDER BY r.log2_fold_change DESC
LIMIT 20

// Genes downregulated by coculture with a specific treatment organism
MATCH (org:OrganismTaxon {organism_name: 'Alteromonas'})-[r:Affects_expression_of]->(g:Gene)
WHERE r.expression_direction = 'down'
RETURN g.locus_tag, g.product, r.log2_fold_change, r.control_condition
ORDER BY r.log2_fold_change ASC
LIMIT 20

// Count DE genes per coculture partner
MATCH (org:OrganismTaxon)-[r:Affects_expression_of]->(g:Gene)
WHERE r.significant = true
RETURN org.organism_name, r.expression_direction, count(g) AS gene_count
ORDER BY org.organism_name, r.expression_direction
```

### Find Genes Affected by Environmental Stress

```cypher
// Genes affected by nutrient stress
MATCH (env:EnvironmentalCondition {condition_type: 'nutrient_stress'})-[r:Affects_expression_of]->(g:Gene)
RETURN g.locus_tag, g.product, r.expression_direction, r.log2_fold_change,
       env.name, env.description
ORDER BY abs(r.log2_fold_change) DESC
LIMIT 20

// Genes affected by light/dark shift
MATCH (env:EnvironmentalCondition {condition_type: 'light_stress'})-[r:Affects_expression_of]->(g:Gene)
RETURN g.locus_tag, g.product, r.expression_direction, r.log2_fold_change,
       env.name, env.light_condition
ORDER BY abs(r.log2_fold_change) DESC
LIMIT 20

// List all distinct environmental conditions in the graph
MATCH (env:EnvironmentalCondition)
RETURN env.name, env.condition_type, env.description,
       env.light_condition, env.oxygen_level
ORDER BY env.condition_type, env.name

// DE genes from all stress experiments
MATCH (env:EnvironmentalCondition)-[r:Affects_expression_of]->(g:Gene)
WHERE r.significant = true
RETURN env.name, env.condition_type, r.expression_direction, count(g) AS gene_count
ORDER BY gene_count DESC
```

### Compare Conditions

```cypher
// Genes affected by BOTH coculture AND nutrient stress
MATCH (org:OrganismTaxon)-[r1:Affects_expression_of]->(g:Gene)
WHERE org.organism_name = 'Alteromonas'
MATCH (env:EnvironmentalCondition {condition_type: 'nutrient_stress'})-[r2:Affects_expression_of]->(g)
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
WHERE org.organism_name = 'Alteromonas'
RETURN cat.name, sc.name, pw.name, count(DISTINCT g) AS gene_count
ORDER BY gene_count DESC
LIMIT 20

// Genes downregulated by nutrient stress in each COG category
MATCH (env:EnvironmentalCondition {condition_type: 'nutrient_stress'})-[r:Affects_expression_of {expression_direction: 'down'}]->(g:Gene)
      -[:Gene_in_cog_category]->(cog:CogFunctionalCategory)
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
      -[:Ko_in_kegg_pathway]->(pw:KeggPathway)
WHERE pw.name CONTAINS 'Photosynthesis'
MATCH (g)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
RETURN g.locus_tag, g.product, ko.name, o.strain_name
ORDER BY o.strain_name

// Genes with a specific EC number
MATCH (g:Gene)-[:Gene_catalyzes_ec_number]->(ec:EcNumber)
WHERE ec.id CONTAINS '2.7'  // kinases
MATCH (g)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
RETURN g.locus_tag, g.product, ec.id, ec.name, o.strain_name
ORDER BY ec.id

// Genes by Cyanorak functional role
MATCH (g:Gene)-[:Gene_has_cyanorak_role]->(role:CyanorakRole)
WHERE role.description CONTAINS 'Photosynthesis'
MATCH (g)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
RETURN g.locus_tag, g.product, role.description, o.strain_name
ORDER BY o.strain_name
```

### Gene + Protein Annotations

```cypher
// Gene with its protein
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_encodes_protein]->(p:Protein)
RETURN g.locus_tag, p.id, p.gene_names, p.is_reviewed,
       p.sequence_length, p.molecular_mass, p.refseq_ids

// Genes with reviewed UniProt proteins
MATCH (g:Gene)-[:Gene_encodes_protein]->(p:Protein {is_reviewed: true})
MATCH (g)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
WHERE o.genus = 'Prochlorococcus'
RETURN g.locus_tag, g.product, p.id, p.gene_names
LIMIT 20

// Genes with GO biological process via linked nodes
MATCH (g:Gene)-[:Gene_involved_in_biological_process]->(bp:BiologicalProcess)
WHERE bp.name CONTAINS 'DNA repair'
MATCH (g)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
RETURN g.locus_tag, g.product, bp.id, bp.name, o.strain_name
```

### GO Hierarchy Traversal

```cypher
// Find all child processes of a GO term (is_a + part_of)
MATCH (child:BiologicalProcess)-[:Biological_process_is_a_biological_process|Biological_process_part_of_biological_process]->(parent:BiologicalProcess {name: 'photosynthesis'})
RETURN child.id, child.name

// GO terms that regulate a process
MATCH (reg:BiologicalProcess)-[:Biological_process_positively_regulates_biological_process|Biological_process_negatively_regulates_biological_process]->(target:BiologicalProcess)
WHERE target.name CONTAINS 'cell cycle'
RETURN reg.id, reg.name, type(reg) AS regulation_type

// Genes in a GO term and all its descendants (transitive closure)
MATCH (bp:BiologicalProcess {name: 'photosynthesis'})
MATCH (child:BiologicalProcess)-[:Biological_process_is_a_biological_process*0..5]->(bp)
MATCH (g:Gene)-[:Gene_involved_in_biological_process]->(child)
MATCH (g)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
RETURN DISTINCT g.locus_tag, g.product, child.name, o.strain_name
ORDER BY o.strain_name
```

### Strain / Organism Queries

```cypher
// All strains in the graph (genome organisms)
MATCH (o:OrganismTaxon)
WHERE o.id STARTS WITH 'insdc.gcf'
RETURN o.strain_name, o.organism_name, o.clade, o.genus, o.ncbi_taxon_id
ORDER BY o.genus, o.clade, o.strain_name

// Treatment organisms (non-genome)
MATCH (o:OrganismTaxon)
WHERE o.id STARTS WITH 'ncbitaxon'
RETURN o.organism_name, o.ncbi_taxon_id

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
MATCH (c)<-[:Gene_in_cyanorak_cluster]-(g2:Gene)
RETURN c.cluster_number, collect(g2.locus_tag) AS genes, strain_count
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
WHERE '10.1038/ismej.2016.70' IN r.publications
RETURN g.locus_tag, g.product, factor.id, r.expression_direction

// Gene count per publication
MATCH (factor)-[r:Affects_expression_of]->(g:Gene)
UNWIND r.publications AS pub
WITH pub, count(DISTINCT g) AS gene_count
RETURN pub, gene_count
ORDER BY gene_count DESC

// Publication node details
MATCH (pub:Publication)
RETURN pub.id, pub.title, pub.journal, pub.publication_year, pub.organism, pub.study_type
ORDER BY pub.publication_year DESC
```

### Graph Statistics

```cypher
// Node counts by label (excluding BioCypher ontology parent labels)
MATCH (n)
WITH labels(n) AS lbls
UNWIND lbls AS l
WITH l WHERE NOT l IN ['NamedThing', 'Entity', 'BiologicalEntity', 'InformationContentEntity',
     'GroupingClass', 'OrganismalEntity', 'AnatomicalEntity', 'Polypeptide',
     'GeneOntology', 'EnvironmentalFeature', 'BiologicalProcessOrActivity', 'Schema_info']
RETURN l AS label, count(*) AS count ORDER BY count DESC

// Edge counts by type
MATCH ()-[r]->() RETURN type(r) AS rel_type, count(*) AS count ORDER BY count DESC

// Genes with no expression data
MATCH (g:Gene)
WHERE NOT (g)<-[:Affects_expression_of]-()
RETURN count(g) AS genes_without_expression

// Annotation quality distribution
MATCH (g:Gene) RETURN g.annotation_quality, count(g) AS count ORDER BY g.annotation_quality

// Genes per organism
MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
RETURN o.strain_name, o.genus, count(g) AS gene_count
ORDER BY gene_count DESC
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
- Treatment organisms (Phage, Marinobacter, etc.) use `organism_name` not `strain_name` (strain_name is null for them)
- EnvironmentalCondition has no `nitrogen_level` or `phosphate_level` — use `condition_type` and `description` for filtering
- Publication DOIs in `r.publications` arrays are bare DOIs (e.g., `10.1038/...`), not prefixed with `doi:`
