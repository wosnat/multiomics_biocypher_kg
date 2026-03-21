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
| `Gene` | `id` (ncbigene:locus_tag), `locus_tag`, `locus_tag_ncbi`, `locus_tag_cyanorak`, `product`, `gene_name`, `gene_synonyms[]`, `function_description`, `cog_category[]`, `go_terms[]`, `kegg_ko[]`, `kegg_pathway[]`, `ec_numbers[]`, `annotation_quality`, `old_locus_tags[]`, `alternative_locus_tags[]`, `protein_id`, `start`, `end`, `strand` |
| `Protein` | `id` (uniprot:accession), `locus_tag`, `gene_names[]`, `is_reviewed`, `annotation_score`, `sequence_length`, `molecular_mass`, `refseq_ids[]`, `proteome_ids[]`, `protein_synonyms[]`, `string_ids[]` |
| `OrganismTaxon` | `id` (insdc.gcf:accession or ncbitaxon:taxid), `organism_name`, `strain_name`, `clade`, `genus`, `species`, `ncbi_taxon_id`, `family`, `order`, `phylum`, `kingdom`, `superkingdom`, `lineage` |
| `Experiment` | `id`, `name`, `organism_strain`, `treatment_type`, `treatment`, `control`, `experimental_context`, `coculture_partner`, `omics_type`, `statistical_test`, `is_time_course`, `medium`, `temperature`, `light_condition`, `light_intensity` |
| `Publication` | `id` (DOI), `title`, `authors[]`, `journal`, `publication_year`, `abstract`, `doi`, `organism`, `study_type`, `description` |
| `OrthologGroup` | `id` (cyanorak:CK_* or eggnog:COG*@*), `name` (raw OG ID), `source` (cyanorak/eggnog), `taxonomic_level` (curated/Prochloraceae/Bacteria/etc.), `taxon_id` |
| `BiologicalProcess` | `id` (go:NNNNNNN), `name` |
| `CellularComponent` | `id` (go:NNNNNNN), `name` |
| `MolecularFunction` | `id` (go:NNNNNNN), `name` |
| `KeggOrthologousGroup` | `id` (K#####), `name` |
| `KeggPathway` | `id` (map#####), `name` |
| `KeggSubcategory` | `id`, `name` |
| `KeggCategory` | `id`, `name` |
| `EcNumber` | `id` (EC x.x.x.x), `name`, `catalytic_activity[]` |
| `CogFunctionalCategory` | `id`, `code`, `name` |
| `CyanorakRole` | `id`, `code`, `name` |
| `TigrRole` | `id`, `code`, `name` |
| `Pfam` | `id` (pfam:PF*), `name` (description), `short_name` (Pfam shortname) |
| `PfamClan` | `id` (pfam.clan:CL*), `name` (clan name) |

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
| `Changes_expression_of` | Experiment → Gene | `time_point`, `time_point_order`, `time_point_hours`, `log2_fold_change`, `adjusted_p_value`, `expression_direction` (up/down), `significant` (significant/not significant/unknown) |
| `Has_experiment` | Publication → Experiment | — |
| `Tests_coculture_with` | Experiment → OrganismTaxon | — |
| `Gene_in_ortholog_group` | Gene → OrthologGroup | — |
| `Gene_belongs_to_organism` | Gene → OrganismTaxon | — |
| `Protein_belongs_to_organism` | Protein → OrganismTaxon | — |
| `Gene_encodes_protein` | Gene → Protein | — |
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
| `Gene_has_pfam` | Gene → Pfam | — |
| `Pfam_in_pfam_clan` | Pfam → PfamClan | — |
| `Biological_process_is_a_biological_process` | BiologicalProcess → BiologicalProcess | — |
| `Biological_process_part_of_biological_process` | BiologicalProcess → BiologicalProcess | — |
| `Biological_process_positively_regulates_biological_process` | BiologicalProcess → BiologicalProcess | — |
| `Biological_process_negatively_regulates_biological_process` | BiologicalProcess → BiologicalProcess | — |
| `Molecular_function_is_a_molecular_function` | MolecularFunction → MolecularFunction | — |
| `Molecular_function_part_of_molecular_function` | MolecularFunction → MolecularFunction | — |
| `Cellular_component_is_a_cellular_component` | CellularComponent → CellularComponent | — |
| `Cellular_component_part_of_cellular_component` | CellularComponent → CellularComponent | — |

### OrthologGroup Sources

Ortholog relationships are modeled as `OrthologGroup` nodes with `Gene_in_ortholog_group` membership edges. Filter by `source` and `taxonomic_level` properties:

| `source` | `taxonomic_level` | Coverage | Description |
|---|---|---|---|
| `cyanorak` | `curated` | Pro/Syn only | Curated Cyanorak ortholog clusters |
| `eggnog` | `Prochloraceae` | Prochlorococcus | Lowest-level eggNOG OG for Pro |
| `eggnog` | `Synechococcus` | Synechococcus | Lowest-level eggNOG OG for Syn |
| `eggnog` | `Alteromonadaceae` | Alteromonas | Lowest-level eggNOG OG for Alt |
| `eggnog` | `Bacteria` | All strains | Bacteria-level COG (cross-phylum bridging) |

Fallback levels (`Cyanobacteria`, `Gammaproteobacteria`) are used when a gene lacks the target-level OG.

### Experiment Treatment Types

The `treatment_type` property on `Experiment` nodes categorizes the experiment. Use it for filtering:

| `treatment_type` | Description |
|---|---|
| `light_stress` | Light/dark shifts, high-light stress |
| `gas_shock` | CO₂/O₂ perturbations |
| `nutrient_stress` | Nitrogen, iron, phosphate limitation |
| `growth_state` | Stationary vs exponential phase |
| `salt_stress` | Salinity changes |
| `pco2` | pCO₂ level experiments |
| `coculture` | Coculture experiments (has `Tests_coculture_with` edge to partner organism) |
| `plastic_leachate_stress` | Plastic leachate exposure |
| `iron_stress` | Iron limitation/addition |
| `dark_tolerance` | Extended darkness survival |
| `viral_lysis_products` | Viral lysate exposure |

Note: Experiment nodes carry `treatment`, `control`, and `experimental_context` for detailed filtering. Coculture experiments also have `coculture_partner` and a `Tests_coculture_with` edge to the partner `OrganismTaxon`.

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
       g.ec_numbers, g.annotation_quality
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
RETURN role.code, role.name, parent.code AS parent_code, parent.name AS parent_desc

// Gene → TIGR role (Pro/Syn only)
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_has_tigr_role]->(role:TigrRole)
RETURN role.code, role.name

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
       collect(DISTINCT crole.name) AS cyanorak_roles
```

### Find What Affects a Gene

```cypher
// All expression edges for a gene
MATCH (exp:Experiment)-[r:Changes_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
RETURN exp.name, exp.treatment_type, exp.treatment, exp.control,
       r.expression_direction, r.log2_fold_change, r.adjusted_p_value,
       r.time_point, r.time_point_hours
ORDER BY abs(r.log2_fold_change) DESC

// All coculture edges for a gene
MATCH (exp:Experiment)-[r:Changes_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
MATCH (exp)-[:Tests_coculture_with]->(partner:OrganismTaxon)
RETURN partner.organism_name, exp.name,
       r.expression_direction, r.log2_fold_change, r.adjusted_p_value
ORDER BY abs(r.log2_fold_change) DESC

// Significantly up-regulated
MATCH (exp:Experiment)-[r:Changes_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
WHERE r.expression_direction = 'up' AND (r.adjusted_p_value IS NULL OR r.adjusted_p_value < 0.05)
RETURN exp.name, exp.treatment_type, r.log2_fold_change, r.adjusted_p_value
ORDER BY r.log2_fold_change DESC

// Filter by treatment_type
MATCH (exp:Experiment {treatment_type: 'nutrient_stress'})-[r:Changes_expression_of]->(g:Gene)
MATCH (g)-[:Gene_belongs_to_organism]->(o:OrganismTaxon {strain_name: 'MED4'})
RETURN g.locus_tag, g.product, r.expression_direction, r.log2_fold_change, exp.name
ORDER BY abs(r.log2_fold_change) DESC
LIMIT 20
```

### Find Genes Affected by Organism Coculture

```cypher
// Genes upregulated in Prochlorococcus by Alteromonas coculture
MATCH (exp:Experiment)-[r:Changes_expression_of {expression_direction: 'up'}]->(g:Gene)
MATCH (exp)-[:Tests_coculture_with]->(partner:OrganismTaxon)
WHERE partner.organism_name = 'Alteromonas'
RETURN g.locus_tag, g.product, g.function_description,
       r.log2_fold_change, r.adjusted_p_value, exp.experimental_context, exp.organism_strain
ORDER BY r.log2_fold_change DESC
LIMIT 20

// Genes downregulated by coculture with a specific treatment organism
MATCH (exp:Experiment)-[r:Changes_expression_of]->(g:Gene)
MATCH (exp)-[:Tests_coculture_with]->(partner:OrganismTaxon {organism_name: 'Alteromonas'})
WHERE r.expression_direction = 'down'
RETURN g.locus_tag, g.product, r.log2_fold_change, exp.control, exp.treatment
ORDER BY r.log2_fold_change ASC
LIMIT 20

// Count DE genes per coculture partner
MATCH (exp:Experiment)-[r:Changes_expression_of]->(g:Gene)
MATCH (exp)-[:Tests_coculture_with]->(partner:OrganismTaxon)
WHERE r.significant = 'significant'
RETURN partner.organism_name, r.expression_direction, count(g) AS gene_count
ORDER BY partner.organism_name, r.expression_direction

// Which publications studied a specific coculture organism
MATCH (pub:Publication)-[:Has_experiment]->(exp:Experiment)-[:Tests_coculture_with]->(partner:OrganismTaxon)
WHERE partner.organism_name CONTAINS 'Alteromonas'
RETURN pub.id, pub.title, pub.publication_year, partner.organism_name, partner.strain_name
ORDER BY pub.publication_year DESC
```

### Find Genes Affected by Environmental Stress

```cypher
// Genes affected by nitrogen stress (using treatment_type)
MATCH (exp:Experiment {treatment_type: 'nutrient_stress'})-[r:Changes_expression_of]->(g:Gene)
RETURN g.locus_tag, g.product, r.expression_direction, r.log2_fold_change,
       exp.name, exp.experimental_context
ORDER BY abs(r.log2_fold_change) DESC
LIMIT 20

// Genes affected by light/dark shift
MATCH (exp:Experiment {treatment_type: 'light_stress'})-[r:Changes_expression_of]->(g:Gene)
RETURN g.locus_tag, g.product, r.expression_direction, r.log2_fold_change,
       exp.name, exp.light_condition
ORDER BY abs(r.log2_fold_change) DESC
LIMIT 20

// List all distinct experiments in the graph
MATCH (exp:Experiment)
RETURN exp.name, exp.treatment_type, exp.organism_strain, exp.treatment, exp.control,
       exp.omics_type, exp.is_time_course
ORDER BY exp.treatment_type, exp.name

// DE genes from all stress experiments
MATCH (exp:Experiment)-[r:Changes_expression_of]->(g:Gene)
WHERE r.significant = 'significant'
RETURN exp.name, exp.treatment_type, r.expression_direction, count(g) AS gene_count
ORDER BY gene_count DESC

// Which publications studied a specific treatment type
MATCH (pub:Publication)-[:Has_experiment]->(exp:Experiment)
WHERE exp.treatment_type = 'iron_stress'
RETURN pub.id, pub.title, pub.publication_year, exp.name
ORDER BY pub.publication_year DESC
```

### Compare Conditions

```cypher
// Genes affected by BOTH coculture AND nutrient stress
MATCH (exp1:Experiment)-[r1:Changes_expression_of]->(g:Gene)
MATCH (exp1)-[:Tests_coculture_with]->(partner:OrganismTaxon {organism_name: 'Alteromonas'})
MATCH (exp2:Experiment {treatment_type: 'nutrient_stress'})-[r2:Changes_expression_of]->(g)
RETURN g.locus_tag, g.product,
       r1.expression_direction AS coculture_dir, r1.log2_fold_change AS coculture_fc,
       r2.expression_direction AS stress_dir, r2.log2_fold_change AS stress_fc

// Genes with OPPOSITE response to coculture vs stress (discordant)
MATCH (exp1:Experiment {treatment_type: 'coculture'})-[r1:Changes_expression_of]->(g:Gene)
MATCH (exp2:Experiment)-[r2:Changes_expression_of]->(g)
WHERE exp2.treatment_type <> 'coculture' AND r1.expression_direction <> r2.expression_direction
RETURN g.locus_tag, g.product, r1.expression_direction AS coculture_dir, r2.expression_direction AS stress_dir

// Genes consistently up across multiple experiments
MATCH (exp:Experiment)-[r:Changes_expression_of]->(g:Gene)
WHERE r.expression_direction = 'up' AND r.significant = 'significant'
WITH g, count(DISTINCT exp) AS up_count
WHERE up_count >= 3
RETURN g.locus_tag, g.product, g.function_description, up_count
ORDER BY up_count DESC

// Cross-paper comparison: how do multiple publications report the same gene under the same treatment_type
MATCH (pub:Publication)-[:Has_experiment]->(exp:Experiment)-[r:Changes_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
RETURN exp.treatment_type, exp.name, r.expression_direction,
       r.log2_fold_change, exp.statistical_test, pub.id AS publication
ORDER BY exp.treatment_type, r.log2_fold_change DESC
```

### Homologs (via OrthologGroup)

```cypher
// All homologs of a gene (all sources, 2-hop via OrthologGroup)
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_in_ortholog_group]->(og:OrthologGroup)<-[:Gene_in_ortholog_group]-(h:Gene)
WHERE g <> h
MATCH (h)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
RETURN h.locus_tag, h.product, o.strain_name, og.source, og.name, og.taxonomic_level
ORDER BY og.source, o.strain_name

// Only Cyanorak-cluster homologs (Pro/Syn ↔ Pro/Syn)
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_in_ortholog_group]->(og:OrthologGroup {source: 'cyanorak'})<-[:Gene_in_ortholog_group]-(h:Gene)
WHERE g <> h
MATCH (h)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
RETURN h.locus_tag, o.strain_name, og.name AS cluster_id

// Cross-phylum homologs (Alteromonas ↔ Pro/Syn via Bacteria-level COG)
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_in_ortholog_group]->(og:OrthologGroup {source: 'eggnog', taxonomic_level: 'Bacteria'})<-[:Gene_in_ortholog_group]-(h:Gene)
WHERE g <> h
MATCH (h)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
RETURN h.locus_tag, h.product, o.strain_name, og.name AS cog_id

// Within-Alteromonas homologs (Alteromonadaceae-level OGs)
MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon {strain_name: 'EZ55'})
MATCH (g)-[:Gene_in_ortholog_group]->(og:OrthologGroup {source: 'eggnog', taxonomic_level: 'Alteromonadaceae'})<-[:Gene_in_ortholog_group]-(h:Gene)
WHERE g <> h
MATCH (h)-[:Gene_belongs_to_organism]->(o2:OrganismTaxon)
RETURN g.locus_tag, h.locus_tag, o2.strain_name, og.name AS og_id
LIMIT 20

// Find members of a specific Cyanorak cluster by name
MATCH (g:Gene)-[:Gene_in_ortholog_group]->(og:OrthologGroup {source: 'cyanorak', name: 'CK_00000364'})
MATCH (g)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
RETURN g.locus_tag, g.product, o.strain_name
```

### Ortholog Expression (query-time propagation)

```cypher
// Expression data for orthologs of a gene (3-hop: gene→OG←homolog←expression)
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_in_ortholog_group]->(og:OrthologGroup)<-[:Gene_in_ortholog_group]-(h:Gene)
WHERE g <> h
MATCH (exp:Experiment)-[e:Changes_expression_of]->(h)
RETURN h.locus_tag, og.source, og.taxonomic_level, e.expression_direction,
       e.log2_fold_change, exp.name AS experiment
ORDER BY abs(e.log2_fold_change) DESC
LIMIT 20

// Coculture ortholog expression (filter by taxonomic level)
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_in_ortholog_group]->(og:OrthologGroup {source: 'cyanorak'})<-[:Gene_in_ortholog_group]-(h:Gene)
WHERE g <> h
MATCH (exp:Experiment)-[e:Changes_expression_of]->(h)
MATCH (exp)-[:Tests_coculture_with]->(partner:OrganismTaxon)
RETURN h.locus_tag, partner.organism_name, e.expression_direction,
       e.log2_fold_change, og.name AS cluster_id
ORDER BY abs(e.log2_fold_change) DESC
LIMIT 20

// Compare direct vs ortholog expression for a gene
MATCH (g:Gene {locus_tag: 'PMM0001'})
OPTIONAL MATCH (exp1:Experiment)-[r1:Changes_expression_of]->(g)
WITH g, collect({experiment: exp1.name, dir: r1.expression_direction, fc: r1.log2_fold_change}) AS direct
OPTIONAL MATCH (g)-[:Gene_in_ortholog_group]->(og:OrthologGroup)<-[:Gene_in_ortholog_group]-(h:Gene)
WHERE g <> h
OPTIONAL MATCH (exp2:Experiment)-[r2:Changes_expression_of]->(h)
RETURN g.locus_tag, direct,
       h.locus_tag AS homolog, og.source, r2.expression_direction AS homolog_dir, r2.log2_fold_change AS homolog_fc
```

### Functional Enrichment / Pathway Analysis

```cypher
// Genes upregulated by Alteromonas coculture in each KEGG category
MATCH (exp:Experiment)-[r:Changes_expression_of {expression_direction: 'up'}]->(g:Gene)
      -[:Gene_has_kegg_ko]->(ko:KeggOrthologousGroup)
      -[:Ko_in_kegg_pathway]->(pw:KeggPathway)
      -[:Kegg_pathway_in_kegg_subcategory]->(sc:KeggSubcategory)
      -[:Kegg_subcategory_in_kegg_category]->(cat:KeggCategory)
MATCH (exp)-[:Tests_coculture_with]->(partner:OrganismTaxon {organism_name: 'Alteromonas'})
RETURN cat.name, sc.name, pw.name, count(DISTINCT g) AS gene_count
ORDER BY gene_count DESC
LIMIT 20

// Genes downregulated by nitrogen stress in each COG category
MATCH (exp:Experiment {treatment_type: 'nutrient_stress'})-[r:Changes_expression_of {expression_direction: 'down'}]->(g:Gene)
      -[:Gene_in_cog_category]->(cog:CogFunctionalCategory)
RETURN cog.code, cog.name, count(DISTINCT g) AS gene_count
ORDER BY gene_count DESC

// Most common COG categories in all DE genes
MATCH (exp:Experiment)-[r:Changes_expression_of]->(g:Gene)-[:Gene_in_cog_category]->(cog:CogFunctionalCategory)
WHERE r.significant = 'significant'
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
WHERE role.name CONTAINS 'Photosynthesis'
MATCH (g)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
RETURN g.locus_tag, g.product, role.name, o.strain_name
ORDER BY o.strain_name
```

### Gene + Protein Annotations

```cypher
// Gene with its protein
MATCH (g:Gene {locus_tag: 'PMM0001'})-[:Gene_encodes_protein]->(p:Protein)
RETURN g.locus_tag, p.id, p.gene_names, p.is_reviewed,
       p.sequence_length, p.molecular_mass, p.refseq_ids

// Genes with reviewed UniProt proteins
MATCH (g:Gene)-[:Gene_encodes_protein]->(p:Protein {is_reviewed: 'reviewed'})
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

### Pfam Domain Queries

```cypher
// Find all genes with a specific Pfam domain (by shortname)
MATCH (g:Gene)-[:Gene_has_pfam]->(p:Pfam {short_name: 'DNA_pol3_beta'})
RETURN g.locus_tag, g.organism_strain

// Find all domains in a clan
MATCH (d:Pfam)-[:Pfam_in_pfam_clan]->(c:PfamClan {name: 'DNA_clamp'})
RETURN d.short_name, d.name

// Find genes in a superfamily (2-hop via clan)
MATCH (c:PfamClan {name: 'DNA_clamp'})<-[:Pfam_in_pfam_clan]-(d:Pfam)<-[:Gene_has_pfam]-(g:Gene)
RETURN d.short_name, g.locus_tag, g.organism_strain

// Full-text search for Pfam domains
CALL db.index.fulltext.queryNodes('pfamFullText', 'polymerase')
YIELD node, score
RETURN node.short_name, node.name, score
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

// Genes unique to one clade (not shared via OrthologGroup with other clades)
MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon {clade: 'HLI'})
WHERE NOT EXISTS {
  MATCH (g)-[:Gene_in_ortholog_group]->(og:OrthologGroup)<-[:Gene_in_ortholog_group]-(h:Gene)
        -[:Gene_belongs_to_organism]->(:OrganismTaxon {clade: 'HLII'})
  WHERE g <> h
}
RETURN g.locus_tag, g.product LIMIT 20

// Genes present in all Pro strains (OrthologGroup shared across all)
MATCH (og:OrthologGroup {source: 'cyanorak'})<-[:Gene_in_ortholog_group]-(g:Gene)
      -[:Gene_belongs_to_organism]->(o:OrganismTaxon)
WHERE o.genus = 'Prochlorococcus'
WITH og, count(DISTINCT o) AS strain_count
WHERE strain_count >= 8
MATCH (og)<-[:Gene_in_ortholog_group]-(g2:Gene)
RETURN og.name AS cluster_id, collect(g2.locus_tag) AS genes, strain_count
ORDER BY strain_count DESC
LIMIT 10
```

### Significance Filtering

```cypher
// Significantly DE genes (p < 0.05, |logFC| > 1)
MATCH (exp:Experiment)-[r:Changes_expression_of]->(g:Gene)
WHERE r.adjusted_p_value < 0.05 AND abs(r.log2_fold_change) > 1
RETURN exp.name, exp.treatment_type, g.locus_tag, g.product,
       r.expression_direction, r.log2_fold_change, r.adjusted_p_value
ORDER BY r.adjusted_p_value ASC
LIMIT 20

// Top genes by fold change from stress experiments
MATCH (exp:Experiment)-[r:Changes_expression_of]->(g:Gene)
WHERE r.significant = 'significant' AND exp.treatment_type <> 'coculture'
RETURN g.locus_tag, g.product, r.expression_direction,
       max(abs(r.log2_fold_change)) AS max_fc
ORDER BY max_fc DESC
LIMIT 20

// Top genes by fold change from coculture experiments
MATCH (exp:Experiment {treatment_type: 'coculture'})-[r:Changes_expression_of]->(g:Gene)
WHERE r.significant = 'significant'
RETURN g.locus_tag, g.product, r.expression_direction,
       max(abs(r.log2_fold_change)) AS max_fc
ORDER BY max_fc DESC
LIMIT 20
```

### Time Series Analysis

```cypher
// Gene expression across time points (ordered by time_point_order)
MATCH (exp:Experiment)-[r:Changes_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
WHERE r.time_point IS NOT NULL
RETURN exp.name, r.time_point, r.time_point_hours, r.time_point_order,
       r.log2_fold_change, r.expression_direction
ORDER BY r.time_point_order

// Time-course experiments only
MATCH (exp:Experiment {is_time_course: true})-[r:Changes_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
RETURN exp.name, r.time_point, r.time_point_hours, r.log2_fold_change, r.expression_direction
ORDER BY exp.name, r.time_point_order

// Genes with consistent direction across all time points in a time-course experiment
MATCH (exp:Experiment {is_time_course: true})-[r:Changes_expression_of]->(g:Gene)
WITH g, exp, collect(r.expression_direction) AS directions
WHERE size([d IN directions WHERE d = 'up']) = size(directions)
   OR size([d IN directions WHERE d = 'down']) = size(directions)
RETURN g.locus_tag, g.product, exp.name, directions
```

### Publication Traceability

```cypher
// Which studies reported expression changes for a gene
MATCH (pub:Publication)-[:Has_experiment]->(exp:Experiment)-[r:Changes_expression_of]->(g:Gene {locus_tag: 'PMM0001'})
RETURN pub.id AS publication, exp.name, r.log2_fold_change, r.expression_direction

// All genes from a specific publication DOI
MATCH (pub:Publication {id: '10.1038/ismej.2016.70'})-[:Has_experiment]->(exp:Experiment)-[r:Changes_expression_of]->(g:Gene)
RETURN g.locus_tag, g.product, exp.name, r.expression_direction

// Gene count per publication
MATCH (pub:Publication)-[:Has_experiment]->(exp:Experiment)-[:Changes_expression_of]->(g:Gene)
WITH pub, count(DISTINCT g) AS gene_count
RETURN pub.id, pub.title, gene_count
ORDER BY gene_count DESC

// Publication node linked to its experiments (via Has_experiment)
MATCH (pub:Publication)-[:Has_experiment]->(exp:Experiment)
RETURN pub.id, pub.title, pub.publication_year,
       exp.name, exp.treatment_type, exp.organism_strain, exp.omics_type
ORDER BY pub.publication_year DESC

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
WHERE NOT (g)<-[:Changes_expression_of]-()
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
- Use `Changes_expression_of` (Experiment → Gene) for all expression queries — not the old `Condition_changes_expression_of` or `Coculture_changes_expression_of`
- For coculture queries, add: `MATCH (exp)-[:Tests_coculture_with]->(partner:OrganismTaxon)`
- For publication linkage, use: `(pub:Publication)-[:Has_experiment]->(exp:Experiment)` — not the old `Published_expression_data_about`
- Experiment properties carry the metadata (treatment_type, organism_strain, treatment, control, omics_type, etc.); expression edge properties are limited to per-row data (log2_fold_change, adjusted_p_value, expression_direction, significant, time_point, time_point_order, time_point_hours)
- For ortholog expression, use 3-hop query-time joins via `Gene_in_ortholog_group` → `OrthologGroup` (no materialized ortholog edges exist)
- Gene properties like `go_terms`, `kegg_ko`, `cog_category` are arrays (`str[]`) on the Gene node itself (denormalized for LLM use) AND available via linked nodes for graph traversal
- Pfam domains are represented only as graph nodes/edges (`Gene_has_pfam`, `Pfam_in_pfam_clan`) — no pfam properties on Gene nodes
- Treatment organisms (Phage, Marinobacter, etc.) use `organism_name` not `strain_name` (strain_name is null for them)
- Use `treatment_type` on Experiment nodes for filtering by experiment category (e.g., `nutrient_stress`, `light_stress`, `coculture`)
- Publication DOIs are stored as Publication node IDs (e.g., `10.1038/...`), not prefixed with `doi:`
