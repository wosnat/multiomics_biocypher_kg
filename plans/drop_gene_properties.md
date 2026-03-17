# Drop Gene Node Properties

Strip 30 properties from Gene nodes that are either redundant with graph edges,
build-time artifacts, or empty. Does NOT touch `gene_annotations_merged.json` —
properties are only removed from schema config and adapter output.

Pfam handled separately (not in scope).

---

## Properties to Drop (30)

### Category 1: Pipeline metadata (6)
| Property | Reason |
|---|---|
| `seed_ortholog` | EggNOG internal |
| `seed_ortholog_evalue` | EggNOG internal |
| `max_annot_lvl` | EggNOG internal |
| `product_source` | Build provenance |
| `function_description_source` | Build provenance |
| `gene_name_source` | Build provenance |

### Category 2: Redundant Cyanorak coordinates + IDs (6)
| Property | Redundant with |
|---|---|
| `start_cyanorak` | `start` |
| `end_cyanorak` | `end` |
| `strand_cyanorak` | `strand` |
| `locus_tag_cyanorak` | `all_identifiers` |
| `locus_tag_ncbi` | `alternative_locus_tags` (confirmed present) |
| `product_cyanorak` | `product` |

### Category 3: Ontology ID lists (10)
| Property | Represented by |
|---|---|
| `go_terms` | `Gene_involved_in_biological_process` etc. edges |
| `go_term_descriptions` | GO node `name` property |
| `kegg_ko` | `Gene_has_kegg_ko` edges |
| `kegg_ko_descriptions` | KeggTerm node `name` |
| `kegg_pathway` | KeggTerm hierarchy edges |
| `kegg_module` | KeggTerm hierarchy edges |
| `kegg_reaction` | KeggTerm hierarchy edges |
| `kegg_brite` | KeggTerm hierarchy edges |
| `ec_numbers` | `Gene_catalyzes_ec_number` edges |
| `ontology_terms` | Never populated (0% coverage) |

### Category 4: Pipeline-internal OG data (3)
| Property | Represented by |
|---|---|
| `eggnog_ogs` | `OrthologGroup` nodes + `Gene_in_ortholog_group` edges |
| `eggnog_og_descriptions` | `OrthologGroup.consensus_product` |
| `bacteria_cog_og` | Never populated (0% coverage) |

### Category 5: Redundant with graph edges (5)
| Property | Redundant with |
|---|---|
| `cog_category` | `Gene_in_cog_category` edges |
| `cyanorak_Role` | `Gene_has_cyanorak_role` edges + `CyanorakRole` node |
| `cyanorak_Role_description` | `CyanorakRole` node `description` property |
| `tIGR_Role` | `Gene_has_tigr_role` edges + `TigrRole` node |
| `tIGR_Role_description` | `TigrRole` node `description` property |

### Category 6: Empty properties (0 — already non-existent)
`subcellular_location` and `uniprot_accession` are not in the current schema or
adapter enum — already absent. No action needed.

## Properties to KEEP (unchanged)

Core: `locus_tag`, `start`, `end`, `strand`, `product`, `protein_id`
Naming: `gene_name`, `gene_synonyms`, `gene_name_synonyms`, `alternative_locus_tags`,
  `alternate_functional_descriptions`, `old_locus_tags`, `function_description`
Pfam: `pfam_ids`, `pfam_names`, `pfam_descriptions` (separate project)
Function: `protein_family`, `catalytic_activities`, `transmembrane_regions`,
  `signal_peptide`, `transporter_classification`, `cazy_ids`, `bigg_reaction`
Computed: `organism_strain`, `gene_summary`, `all_identifiers`
Quality: `annotation_quality`, `gene_category`

---

## Implementation Plan

### Phase 1: Schema config
**File**: `config/schema_config.yaml` (lines 150-220)

Remove these 30 property lines from the `gene:` block:
- `locus_tag_ncbi`, `locus_tag_cyanorak`, `start_cyanorak`, `end_cyanorak`, `strand_cyanorak`, `product_cyanorak`
- `cyanorak_Role`, `cyanorak_Role_description`
- `tIGR_Role`, `tIGR_Role_description`
- `ec_numbers`
- `go_terms`, `go_term_descriptions`
- `kegg_ko`, `kegg_ko_descriptions`, `kegg_pathway`, `kegg_module`, `kegg_reaction`, `kegg_brite`
- `cog_category`, `eggnog_ogs`, `eggnog_og_descriptions`, `seed_ortholog`, `max_annot_lvl`, `seed_ortholog_evalue`
- `gene_name_source`, `product_source`, `function_description_source`

Also remove associated comments (e.g., `# EggNOG / COG functional roles`, `# Gene Ontology`, `# KEGG`).

### Phase 2: Adapter enum
**File**: `multiomics_kg/adapters/cyanorak_ncbi_adapter.py` (lines 99-166)

Remove matching `GeneNodeField` enum members:
- `LOCUS_TAG_NCBI`, `LOCUS_TAG_CYANORAK`, `START_CYANORAK`, `END_CYANORAK`, `STRAND_CYANORAK`, `PRODUCT_CYANORAK`
- `CYANORAK_ROLE`, `CYANORAK_ROLE_DESCRIPTION`
- `TIGR_ROLE`, `TIGR_ROLE_DESCRIPTION`
- `EC_NUMBERS`
- `GO_TERMS`, `GO_TERM_DESCRIPTIONS`
- `KEGG_KO`, `KEGG_KO_DESCRIPTIONS`, `KEGG_PATHWAY`, `KEGG_MODULE`, `KEGG_REACTION`, `KEGG_BRITE`
- `COG_CATEGORY`, `EGGNOG_OGS`, `EGGNOG_OG_DESCRIPTIONS`, `SEED_ORTHOLOG`, `MAX_ANNOT_LVL`, `SEED_ORTHOLOG_EVALUE`
- `GENE_NAME_SOURCE`, `PRODUCT_SOURCE`, `FUNCTION_DESCRIPTION_SOURCE`

Also update the `_get_gene_nodes()` method:
- Remove `start_cyanorak`, `end_cyanorak` from `int_fields` set (line 292)
- Remove `seed_ortholog_evalue` from `float_fields` set (line 293)

### Phase 3: Unit tests
**File**: `tests/test_cyanorak_ncbi_adapter.py`

1. **Update `MERGED_JSON_DATA` fixture** — remove dropped fields from the test data dict
2. **Update enum tests** — remove assertions for dropped enum members (e.g., `GeneNodeField.GO_TERMS`, `GeneNodeField.KEGG_KO`, `GeneNodeField.COG_CATEGORY`, `GeneNodeField.SEED_ORTHOLOG_EVALUE`)
3. **Update property-presence tests** — remove assertions that check dropped properties exist on gene nodes (e.g., `go_terms`, `kegg_ko`, `eggnog_ogs`, `cog_category`, `seed_ortholog`)
4. **Update `go_term_descriptions` list test** — remove (lines ~550-556)
5. **Add new test**: `test_dropped_properties_not_in_gene_nodes` — verify none of the 30 dropped property names appear in gene node output
6. **Update integration tests** — remove assertions for `go_terms`/`kegg_ko` list checks (lines ~1026-1038)

### Phase 4: Run tests
```bash
# Unit tests
pytest tests/test_cyanorak_ncbi_adapter.py -v
pytest -m "not slow and not kg" -v

# If all pass, rebuild KG and run validity tests
pytest -m kg -v
```

### Phase 5: Code review checklist
- [ ] Every dropped property is removed from ALL three locations (schema, enum, test fixture)
- [ ] No dropped property name appears anywhere in `schema_config.yaml` gene block
- [ ] No dropped enum member remains in `GeneNodeField`
- [ ] `int_fields` and `float_fields` sets in `_get_gene_nodes()` updated
- [ ] Post-import scripts (`post-import.sh`, `post-import.cypher`) don't reference dropped properties (verified: they don't)
- [ ] Full-text index `geneFullText` fields still valid (`gene_summary`, `gene_synonyms`, `alternate_functional_descriptions`, `pfam_names` — all kept)
- [ ] No other adapter reads `GeneNodeField` members that were removed
- [ ] KG snapshot tests (`test_snapshot.py`) don't assert dropped properties (verified: they don't)
- [ ] `gene_annotations_merged.json` is NOT modified (out of scope)
- [ ] `CLAUDE.md` "Key graph facts" and property lists updated

### Phase 6: Documentation
- Update `CLAUDE.md` — remove dropped properties from schema documentation
- Update memory if needed

---

## Resolved Decisions

1. **`tIGR_Role` / `tIGR_Role_description`** — drop both. `TigrRole` nodes carry a `description` property, so Gene-level duplication is unnecessary.
2. **`locus_tag_ncbi`** — drop. Confirmed present in `alternative_locus_tags` (e.g., PMM1428 has `TX50_RS07695` in both).

No open decisions remain.

---

## Expected Outcome

Gene nodes go from ~50 properties to ~22. Saves ~45% Neo4j storage per gene node.
No data loss — all dropped info is either in `gene_annotations_merged.json` (queryable
if ever needed) or represented as graph edges/nodes.
