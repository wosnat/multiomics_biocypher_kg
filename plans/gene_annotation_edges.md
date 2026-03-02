# Plan: Convert Gene Annotation Properties to Graph Edges

## Context

Gene nodes currently store functional annotations (GO terms, EC numbers, KEGG KOs, COG categories, Pfam domains, roles, etc.) as flat array properties. These should be proper graph edges to ontology/classification nodes so that:

1. **Graph traversal works** — LLM agents can query "genes involved in GO:0006260" via edge traversal, not property filtering
2. **Orphan-protein gap is fixed** — ~46% of proteins lack `Gene_encodes_protein` edges, so `Gene→Protein→GO` fails for them; direct `Gene→GO` edges bypass this
3. **New node types enable cross-gene queries** — "all genes with KO:K02338", "all genes in COG category L", etc.

Properties are **kept on Gene nodes** too (for fast `WHERE 'GO:0005737' IN g.go_terms` filtering).

---

## Annotation Source Reference

Functional annotations flow through two upstream pipelines depending on the organism:

### Prochlorococcus / Synechococcus (Cyanorak strains)

Described in Doré et al. (2020), *Front. Microbiology*, doi:10.3389/fmicb.2020.567431.

Two-stage process:
1. **IGS Annotation Engine** (U. Maryland) — automatic structural + functional annotation; assigns **tIGR_Role** numeric codes using the legacy TIGR role classification system
2. **Cyanorak v2.1** — expert manual curation by domain specialists; assigns **cyanorak_Role** codes using a custom 3-level cyanobacteria-specific role hierarchy

COG annotations in `gene_mapping.csv` (`eggNOG` field) also come from Cyanorak and are hand-curated.

### Alteromonas strains (MIT1002, EZ55, HOT1A3)

**Not in Cyanorak.** No `cyanorak_Role`, `tIGR_Role`, or Cyanorak-sourced `eggNOG` data. Functional annotations come entirely from:
- **eggnog-mapper v2** (automated): COG categories/OGs, KEGG KO/pathways, Pfam, GO
- **UniProt**: GO terms, EC numbers, Pfam

This means `CyanorakRole` and `TigrRole` nodes/edges will only exist for Pro/Syn strains. **COG is the only functional classification that covers all 13 strains.**

---

## Source Priority and Coverage (per MED4 as representative)

| Field | Cyanorak coverage | eggnog-mapper coverage | Priority |
|---|---|---|---|
| COG#### OG IDs | 75% (hand-curated) | 98% (automated) | **Cyanorak first**, emapper fills gaps |
| COG category letter | extractable from description `[L]` | 88% direct | Derive from Cyanorak desc; fall back to emapper |
| cyanorak_Role | 76% | — | Cyanorak only |
| tIGR_Role | 89% | — | tIGR fills ~13% gap not covered by cyanorak_Role |
| GO terms | Cyanorak + NCBI | eggnog | Union of all |
| EC numbers | Cyanorak | eggnog + UniProt | Union of all |

**Important**: When COG#### IDs differ between Cyanorak and eggnog-mapper (~14% of cases), prefer the **Cyanorak** ID — it is hand-curated. Use eggnog-mapper only for genes Cyanorak lacks (no protein_id match or Alteromonas).

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

COG category nodes: 25 hardcoded nodes using the standard COG functional category descriptions (no download needed).

KEGG pathway names are not in `gene_annotations_merged.json` — create stub nodes with ID only initially; names can be enriched later from KEGG API.

### Tier 3: Project-specific classification

| Gene property | New edge type | New node type | New node ID | Strains |
|---|---|---|---|---|
| `eggnog_ogs` (COG#### only) | `gene_in_cog_og` | CogOrthologousGroup | `cog:COG####` | All 13 |
| `cyanorak_Role` + `cyanorak_Role_description` | `gene_has_cyanorak_role` | CyanorakRole | `cyanorak.role:<code>` | Pro/Syn only |
| `tIGR_Role` + `tIGR_Role_description` | `gene_has_tigr_role` | TigrRole | `tigr.role:<code>` | Pro/Syn only |

eggNOG OGs: extract `COG####` IDs from two formats:
- Plain `COG####` entries (Pro/Syn genes, sourced from Cyanorak hand-curation)
- `COG####@1|root` or `COG####@2|Bacteria` at-format entries (Alteromonas genes from emapper, which have **no plain COG#### entries at all**)

Both formats refer to the same COG database OG and produce the same node ID `cog:COG####`. Skip `bactNOG`, `cyaNOG`, `NOG####`, and non-COG at-format entries (e.g. `1RMNP@1236|Gammaproteobacteria`).

tIGR Roles are retained despite being automated because: (a) they cover ~13% of Pro/Syn genes that lack a Cyanorak Role, and (b) for Pro/Syn they represent the upstream annotation that Cyanorak experts refined.

---

## Hierarchy Considerations

### COG (2-level, simple)

COG has only two functional levels — no deep hierarchy like GO:
- **CogFunctionalCategory** (25 letter codes, hardcoded)
- **CogOrthologousGroup** (COG####, ~5K groups, each assigned to 1–2 categories)

Add a `cog_og_in_cog_category` edge from each `CogOrthologousGroup` node to its 1–2 `CogFunctionalCategory` nodes. The category letter is extractable from `eggnog_og_descriptions` strings (the `[L]` bracket format in Cyanorak descriptions) or from `cog_category` when a gene has only one COG####.

Optional: add 4 `CogSuperCategory` nodes (INFORMATION STORAGE AND PROCESSING / CELLULAR PROCESSES AND SIGNALING / METABOLISM / POORLY CHARACTERIZED) with `cog_category_in_super_category` edges. Trivially hardcoded, low priority.

Do **not** model bactNOG/cyaNOG as graph nodes — they are taxonomic re-clusterings at narrower scope, not a functional subsumption hierarchy.

### Cyanorak Roles (3-level tree, hardcoded)

Cyanorak Roles form a **clean tree** with up to 3 levels:
```
B  (Biosynthesis of cofactors, prosthetic groups, and carriers)
  B.5  (Pigments)
    B.5.1  (Carotenoids)
    B.5.2  (Chlorophylls and porphyrins)
    B.5.3  (Hemes and phycobilins)
  B.10  (Vitamins)
    B.10.1  (Biotin)
    ...
```

Genes are assigned at the **most specific applicable level** — the `cyanorak_Role` field stores only leaf codes; parent codes are not stored alongside them but are fully derivable: strip the last `.N` component.

MED4 distribution: 58 level-1 assignments, 1,497 level-2, 398 level-3. Multiple roles per gene are common (comma-separated).

**Hierarchy edges**: Add `cyanorak_role_is_a_cyanorak_role` (parent→child) edges. The full ~150-node tree with ~130 parent-child edges should be **hardcoded as a static dict** in the adapter from the complete hierarchy definition. This enables queries like "all genes involved in cofactor biosynthesis" via ancestor traversal.

Add new edge type to schema:
```yaml
cyanorak_role_is_a_cyanorak_role:
  is_a: association
  source: cyanorak role
  target: cyanorak role
  label_as_edge: CYANORAK_ROLE_IS_A_CYANORAK_ROLE
```

Also add `cog_og_in_cog_category` to schema:
```yaml
cog_og_in_cog_category:
  is_a: association
  source: cog orthologous group
  target: cog functional category
  label_as_edge: COG_OG_IN_COG_CATEGORY
```

### tIGR Roles (flat, no hierarchy)

tIGR codes are opaque integers (e.g., 132, 156) with no derivable parent structure in the code itself. The description embeds a 2-level hint via ` / ` separator ("DNA metabolism / DNA replication, recombination, and repair") but building hierarchy from description text is fragile. Model as **flat nodes only** — no hierarchy edges.

### GO and EC

Already have full hierarchy implementations elsewhere (go_adapter.py, ec_adapter.py). Tier 1 edges reuse those existing nodes.

---

## Implementation Architecture

### New adapter: `multiomics_kg/adapters/functional_annotation_adapter.py`

**`FunctionalAnnotationAdapter(genome_dir: Path, go_ontology)`** — per-strain adapter:
- `__init__`: loads `gene_annotations_merged.json` for the strain
- `get_edges()`: yields edges for all term types (no new nodes — those come from Multi wrapper)
- Edge ID format: `{locus_tag}-{edge_type}-{term_id}` (e.g., `PMM0001-go-GO:0006260`)
- Gracefully yields nothing for Cyanorak/tIGR role edges when those fields are absent (Alteromonas)

**`MultiFunctionalAnnotationAdapter(genomes_df: DataFrame)`** — multi-strain wrapper:
- `__init__`: loads GO ontology once (reuse go_adapter's Pronto object or reload from OBO)
- `get_nodes()`: aggregates all term IDs across all strains → deduplicates → yields unique nodes for:
  - KO, Pathway, CogFunctionalCategory (25, hardcoded), CogOrthologousGroup, PfamFamily
  - CyanorakRole (~150 nodes, hardcoded tree including all ancestors, not just leaf codes in data)
  - TigrRole (only codes actually present in data)
- `get_edges()`: iterates per-strain adapters, yields all gene→term edges; also yields:
  - `cog_og_in_cog_category` edges (COG#### → category letter)
  - `cyanorak_role_is_a_cyanorak_role` edges (child → parent, full tree)

Follows the same Multi* pattern as `MultiCyanorakNcbi` in `cyanorak_ncbi_adapter.py:563-645`.

### Cyanorak Role tree (hardcoded)

Define a module-level constant in the adapter:
```python
CYANORAK_ROLE_TREE: dict[str, tuple[str, str]] = {
    # code: (description, parent_code_or_None)
    "0": ("Non-coding gene (RNA)", None),
    "0.1": ("tRNA", "0"),
    "0.2": ("rRNA", "0"),
    ...
    "B": ("Biosynthesis of cofactors, prosthetic groups, and carriers", None),
    "B.5": ("Pigments", "B"),
    "B.5.1": ("Carotenoids", "B.5"),
    "B.5.2": ("Chlorophylls and porphyrins", "B.5"),
    ...
}
```

All ~150 entries derived from the complete Cyanorak role hierarchy. `get_nodes()` yields all nodes regardless of whether any gene uses them (enables consistent graph structure). `get_edges()` yields only the parent→child edges implied by the tree.

### Modified files

| File | Change |
|---|---|
| `config/schema_config.yaml` | Add 7 new node types + 13 new edge type definitions (11 original + `cyanorak_role_is_a_cyanorak_role` + `cog_og_in_cog_category`) |
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

cog_og_in_cog_category:
  is_a: association
  source: cog orthologous group
  target: cog functional category
  label_as_edge: COG_OG_IN_COG_CATEGORY

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

cyanorak_role_is_a_cyanorak_role:
  is_a: association
  source: cyanorak role
  target: cyanorak role
  label_as_edge: CYANORAK_ROLE_IS_A_CYANORAK_ROLE

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
3. **Update `config/schema_config.yaml`** — add 7 node types and 13 edge types
4. **Create `multiomics_kg/adapters/functional_annotation_adapter.py`**:
   - `CYANORAK_ROLE_TREE` constant dict (~150 entries, full hierarchy from Cyanorak website)
   - `COG_FUNCTIONAL_CATEGORIES` constant dict (25 standard COG letter codes + descriptions)
   - `FunctionalAnnotationAdapter` (per-strain, reads gene_annotations_merged.json; skips absent fields gracefully)
   - `MultiFunctionalAnnotationAdapter` (deduplicates nodes, aggregates edges, emits hierarchy edges)
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
ls biocypher-out/*/CyanorakRole*.csv
ls biocypher-out/*/Cyanorak_role_is_a_cyanorak_role*.csv

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
- Gene→COG category: ~8k edges (most genes; all 13 strains)
- Gene→COG OG: ~6k edges (all 13 strains, including Alteromonas via at-format extraction)
- COG OG→COG category: ~1-2 per OG (small, bounded by ~5K unique OGs)
- Gene→Pfam: ~5k edges
- Gene→CyanorakRole: ~6k edges (Pro/Syn only; leaf assignments only)
- CyanorakRole→parent: ~130 edges (full tree, hardcoded)
- Gene→tIGR Role: ~4k edges (Pro/Syn only; covers gap not in CyanorakRole)

---

## Separate Plan: Alteromonas Ortholog Edges

This is a **distinct concern** from functional annotation. The existing homolog system uses Cyanorak CLOGs to generate `Gene_is_homolog_of_gene` edges for Pro/Syn (via `scripts/post-import.cypher`). Alteromonas is outside Cyanorak and needs its own mechanism.

### Two scopes

**Within-Alteromonas (MIT1002 ↔ EZ55 ↔ HOT1A3)**

Use eggNOG OGs at `72275|Alteromonadaceae` level. These are the most specific OGs shared across all three Alteromonas strains. IDs are alphanumeric (e.g. `4648R@72275|Alteromonadaceae`). Coverage: ~3,424 genes per strain.

Gene pairs from different Alteromonas strains sharing the same Alteromonadaceae-level OG ID → `Gene_is_homolog_of_gene` edges, analogous to how Cyanorak clusters work for Pro/Syn.

**Cross-phylum (Alteromonas ↔ Pro/Syn)**

Use COG#### at `2|Bacteria` level. Extract `COG####` from at-format entries; this is the same ID space as Pro/Syn COG OG nodes. Genes from different phyla sharing a COG OG → `Gene_is_homolog_of_gene` (or a new edge type `Gene_is_distant_homolog_of_gene` to distinguish from same-clade homologs).

Observed shared COGs between MIT1002 and MED4: **869** — mostly deeply conserved housekeeping genes (DNA replication, ribosomal proteins, core metabolism), which is the expected and useful set for cross-phylum comparison.

### eggNOG OG levels in Alteromonas data

| Taxonomic level | Taxon ID | Genes covered | Use for |
|---|---|---|---|
| root | 1 | ~4,107 | (too broad) |
| Bacteria | 2 | ~4,116 | Cross-phylum ortholog queries |
| Proteobacteria | 1224 | ~3,834 | (intermediate) |
| Gammaproteobacteria | 1236 | ~3,772 | (intermediate) |
| **Alteromonadaceae** | **72275** | **~3,424** | **Within-Alteromonas orthologs** |

### Implementation notes

- Store the Alteromonadaceae-level OG ID as a property on Gene nodes (`alteromonadaceae_og`) to enable fast lookup without traversing all eggnog_ogs entries
- The `Gene_is_homolog_of_gene` edge generation for Alteromonas can be done in a post-import Cypher script analogous to `scripts/post-import.cypher`, grouping genes by shared Alteromonadaceae OG ID
- Cross-phylum homolog edges can be generated from shared `CogOrthologousGroup` node membership (genes already connected via `gene_in_cog_og` edges), avoiding the need for a separate edge generation step — query pattern: `(g1:Gene)-[:GENE_IN_COG_OG]->(cog)<-[:GENE_IN_COG_OG]-(g2:Gene)` where g1 and g2 are from different organisms
- This also means the `gene_in_cog_og` edges (from the functional annotation plan above) **double as the cross-phylum ortholog index** — no separate edge storage needed for that scope