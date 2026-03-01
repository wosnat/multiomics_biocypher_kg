# Plan: Convert Gene Annotation Properties to Graph Edges

## Context

Gene nodes currently store functional annotations (GO terms, EC numbers, KEGG KOs, COG categories, Pfam domains, roles, etc.) as flat array properties. These should be proper graph edges to ontology/classification nodes so that:

1. **Graph traversal works** — LLM agents can query "genes involved in GO:0006260" via edge traversal, not property filtering
2. **Orphan-protein gap is fixed** — ~46% of proteins lack `Gene_encodes_protein` edges, so `Gene→Protein→GO` fails for them; direct `Gene→GO` edges bypass this
3. **New node types enable cross-gene queries** — "all genes with KO:K02338", "all genes in COG category L", etc.

Properties are **kept on Gene nodes** too (for fast `WHERE 'GO:0005737' IN g.go_terms` filtering).

---

## Properties to Convert

### Tier 1: Target nodes already exist

| Gene property | New edge type | Target node type | Target ID |
|---|---|---|---|
| `go_terms` (BP subset) | `gene_involved_in_biological_process` | BiologicalProcess | `GO:XXXXXXX` |
| `go_terms` (CC subset) | `gene_located_in_cellular_component` | CellularComponent | `GO:XXXXXXX` |
| `go_terms` (MF subset) | `gene_enables_molecular_function` | MolecularFunction | `GO:XXXXXXX` |
| `ec_numbers` | `gene_catalyzes_ec_number` | EcNumber | `eccode:X.X.X.X` |

**GO namespace lookup**: Use go_adapter's loaded Pronto ontology to determine each GO term's namespace (biological_process / cellular_component / molecular_function) at edge-generation time. Extract as shared utility.

### Tier 2: New node types (straightforward)

| Gene property | New edge type | New node type | New node ID |
|---|---|---|---|
| `kegg_ko` + `kegg_ko_descriptions` | `gene_has_kegg_ko` | KeggOrthologousGroup | `kegg.orthology:K#####` |
| `kegg_pathway` | `gene_in_kegg_pathway` | KeggPathway | `kegg.pathway:ko#####` |
| `cog_category` + letter description | `gene_in_cog_category` | CogFunctionalCategory | `cog.category:<letter>` |
| `pfam_ids` + `pfam_names` + `pfam_descriptions` | `gene_has_pfam_domain` | PfamFamily | `pfam:PF#####` |

COG category nodes (A–Z): 25 hardcoded nodes using the standard COG functional category descriptions (no download needed).

KEGG pathway names are not in `gene_annotations_merged.json` — create stub nodes with ID only initially; names can be enriched later from KEGG API.

### Tier 3: Project-specific classification

| Gene property | New edge type | New node type | New node ID |
|---|---|---|---|
| `eggnog_ogs` (COG#### only) | `gene_in_cog_og` | CogOrthologousGroup | `cog:COG####` |
| `cyanorak_Role` + `cyanorak_Role_description` | `gene_has_cyanorak_role` | CyanorakRole | `cyanorak.role:<code>` |
| `tIGR_Role` + `tIGR_Role_description` | `gene_has_tigr_role` | TigrRole | `tigr.role:<code>` |

eggNOG OGs: filter `eggnog_ogs` to only `COG####`-prefixed entries for now (skip bactNOG/cyaNOG which are less standard).

---

## Implementation Architecture

### New adapter: `multiomics_kg/adapters/functional_annotation_adapter.py`

**`FunctionalAnnotationAdapter(genome_dir: Path, go_ontology)`** — per-strain adapter:
- `__init__`: loads `gene_annotations_merged.json` for the strain
- `get_edges()`: yields edges for all term types (no new nodes — those come from Multi wrapper)
- Edge ID format: `{locus_tag}-{edge_type}-{term_id}` (e.g., `PMM0001-go-GO:0006260`)

**`MultiFunctionalAnnotationAdapter(genomes_df: DataFrame)`** — multi-strain wrapper:
- `__init__`: loads GO ontology once (reuse go_adapter's Pronto object or reload from OBO)
- `get_nodes()`: aggregates all term IDs across all strains → deduplicates → yields unique KO/Pathway/COG/Pfam/Role nodes
- `get_edges()`: iterates per-strain adapters, yields all edges

Follows the same Multi* pattern as `MultiCyanorakNcbi` in `cyanorak_ncbi_adapter.py:563-645`.

### Shared GO utility: `multiomics_kg/utils/go_utils.py`

```python
def load_go_ontology(obo_path: Path) -> dict[str, str]:
    """Returns {go_id: namespace} where namespace is 'biological_process',
    'cellular_component', or 'molecular_function'."""
```

The go_adapter.py already loads an OBO file (find the path in `go_adapter.py`). Extract the OBO path as a config value so both adapters share it.

### Modified files

| File | Change |
|---|---|
| `config/schema_config.yaml` | Add 7 new node types + 11 new edge type definitions |
| `multiomics_kg/adapters/functional_annotation_adapter.py` | **New file** — see above |
| `multiomics_kg/utils/go_utils.py` | **New utility** — GO namespace lookup |
| `create_knowledge_graph.py` | Instantiate `MultiFunctionalAnnotationAdapter` and add to pipeline |

---

## Schema Changes (`config/schema_config.yaml`)

### New node types
```yaml
kegg orthologous group:
  is_a: biological entity
  preferred_id: kegg.orthology
  properties:
    kegg_ko_id: str
    name: str

kegg pathway:
  is_a: pathway
  preferred_id: kegg.pathway
  properties:
    kegg_pathway_id: str
    name: str

cog functional category:
  is_a: biological entity
  preferred_id: cog.category
  properties:
    code: str
    name: str
    description: str

cog orthologous group:
  is_a: biological entity
  preferred_id: cog
  properties:
    cog_id: str
    description: str

pfam family:
  is_a: protein domain
  preferred_id: pfam
  properties:
    pfam_id: str
    name: str
    description: str

cyanorak role:
  is_a: biological entity
  preferred_id: cyanorak.role
  properties:
    code: str
    description: str

tigr role:
  is_a: biological entity
  preferred_id: tigr.role
  properties:
    code: str
    description: str
```

### New edge types
```yaml
gene_involved_in_biological_process:
  is_a: gene to biological process association
  source: gene
  target: biological process
  label_as_edge: GENE_INVOLVED_IN_BIOLOGICAL_PROCESS

gene_located_in_cellular_component:
  is_a: gene to cellular component association
  source: gene
  target: cellular component
  label_as_edge: GENE_LOCATED_IN_CELLULAR_COMPONENT

gene_enables_molecular_function:
  is_a: gene to molecular activity association
  source: gene
  target: molecular activity
  label_as_edge: GENE_ENABLES_MOLECULAR_FUNCTION

gene_catalyzes_ec_number:
  is_a: association
  source: gene
  target: ec number
  label_as_edge: GENE_CATALYZES_EC_NUMBER

gene_has_kegg_ko:
  is_a: association
  source: gene
  target: kegg orthologous group
  label_as_edge: GENE_HAS_KEGG_KO

gene_in_kegg_pathway:
  is_a: association
  source: gene
  target: kegg pathway
  label_as_edge: GENE_IN_KEGG_PATHWAY

gene_in_cog_category:
  is_a: association
  source: gene
  target: cog functional category
  label_as_edge: GENE_IN_COG_CATEGORY

gene_in_cog_og:
  is_a: association
  source: gene
  target: cog orthologous group
  label_as_edge: GENE_IN_COG_OG

gene_has_pfam_domain:
  is_a: association
  source: gene
  target: pfam family
  label_as_edge: GENE_HAS_PFAM_DOMAIN

gene_has_cyanorak_role:
  is_a: association
  source: gene
  target: cyanorak role
  label_as_edge: GENE_HAS_CYANORAK_ROLE

gene_has_tigr_role:
  is_a: association
  source: gene
  target: tigr role
  label_as_edge: GENE_HAS_TIGR_ROLE
```

---

## Step-by-Step Implementation

1. **Explore go_adapter.py** — find the OBO file path and Pronto usage to extract into `go_utils.py`
2. **Create `multiomics_kg/utils/go_utils.py`** — `load_go_ontology()` returns `{go_id: namespace}` dict
3. **Update `config/schema_config.yaml`** — add 7 node types and 11 edge types
4. **Create `multiomics_kg/adapters/functional_annotation_adapter.py`**:
   - `CogFunctionalCategory` constants dict (A-Z codes with standard descriptions)
   - `FunctionalAnnotationAdapter` (per-strain, reads gene_annotations_merged.json)
   - `MultiFunctionalAnnotationAdapter` (deduplicates nodes, aggregates edges)
5. **Update `create_knowledge_graph.py`** — instantiate and add `MultiFunctionalAnnotationAdapter`
6. **Add tests** — `tests/test_functional_annotation_adapter.py` (unit tests, no KG needed)

---

## Verification

```bash
# Quick build (TEST_MODE=True, CACHE=True in create_knowledge_graph.py)
uv run python create_knowledge_graph.py

# Check output CSVs exist and have rows
ls biocypher-out/*/Gene_involved_in_biological_process*.csv
ls biocypher-out/*/KeggOrthologousGroup*.csv

# Unit tests
pytest tests/test_functional_annotation_adapter.py -v
pytest -m "not slow and not kg"

# KG validity (after docker compose up -d)
pytest tests/kg_validity/ -v
```

Expected new edge counts (approximate, all strains):
- Gene→GO: ~10k edges
- Gene→EC: ~2k edges
- Gene→KO: ~3k edges
- Gene→COG category: ~8k edges (most genes have a category)
- Gene→Pfam: ~5k edges
