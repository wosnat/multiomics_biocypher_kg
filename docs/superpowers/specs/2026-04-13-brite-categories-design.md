# KEGG BRITE Categories — Design Spec

**Date:** 2026-04-13
**Status:** Draft, ready for planning

## Motivation

42% of KEGG Orthology (KO) terms in the KG have no `Kegg_term_is_a_kegg_term` parent edge — 1,824 of 4,367 KOs. This is not a build bug. KEGG simply does not assign every KO to a pathway map. Many KOs (transporters, proteases, restriction enzymes, uncharacterized proteins) live only in KEGG's BRITE functional hierarchies, which this KG currently does not ingest.

**Impact on MED4**: 295 of 1,017 KOs attached to MED4 genes (29%) are pathway-orphans, stranding 332 gene-KO assignments with no hierarchical rollup. Classes most affected: transporters (230 KOs, ~2,000 gene edges across the KG) and peptidases (69 KOs, ~626 gene edges) — both central to research questions about nutrient uptake and protein turnover under stress.

Loading KEGG BRITE trees restores hierarchical classification for these KOs and enables functional enrichment in axes that pathway maps do not cover.

## Scope

### In scope

Ingest 12 KEGG BRITE trees as new hierarchical ontology nodes in the KG:

| Tree code | Name | Primary use |
|---|---|---|
| ko01000 | Enzymes | EC-based enzyme classification (overlaps existing EcNumber hierarchy; kept for BRITE-native navigation and consistency) |
| ko02000 | Transporters | substrate transport classification |
| ko01002 | Peptidases & inhibitors | protease family classification |
| ko03000 | Transcription factors | regulatory protein classification |
| ko02044 | Secretion system | protein export machinery |
| ko02022 | Two-component system | signal transduction pairs |
| ko02048 | Prokaryotic defense system | CRISPR / restriction-modification / T-A systems |
| ko03110 | Chaperones & folding catalysts | protein folding machinery |
| ko03011 | Ribosome | ribosomal protein classification |
| ko03012 | Translation factors | translation machinery |
| ko03016 | Transfer RNA biogenesis | tRNA processing and modification |
| ko03032 | DNA replication proteins | replisome classification |

All trees are loaded. The Enzyme tree (ko01000) overlaps with the existing `EcNumber` hierarchy — the two remain independent and users/tools can query either surface. Tree depth varies (observed max ≤ 4 non-KO levels); actual per-tree depth is determined at first build.

### Out of scope

- **Enrichment tool changes**, new MCP tools, Cypher recipe additions, and explorer-side documentation. Logged in Appendix A as suggestions for the explorer project; to be brainstormed and planned separately.
- **Gene property changes**. No existing Gene property is removed or renamed by this work. BRITE classification is reached via `Gene_has_kegg_ko → Kegg_term_in_brite_category`, not via direct Gene properties.

## Data model

### New node label: `BriteCategory`

| Property | Type | Source |
|---|---|---|
| `id` | str | `normalize_curie(f'kegg.brite:{tree_code}')` for tree root; `normalize_curie(f'kegg.brite:{tree_code}.{path_code}')` for internal nodes, where `path_code` is a positional path string (e.g. `A3.B2.C1` = 3rd top-level group, 2nd subgroup, 1st family). Positional indexing is stable across builds as long as KEGG does not reorder entries (they don't in practice). Example root: `kegg.brite:ko02000`; example leaf: `kegg.brite:ko02000.A3.B2`. |
| `name` | str | KEGG BRITE label, `_clean_str` sanitized |
| `tree` | str | One of: `enzymes`, `transporters`, `peptidases`, `transcription_factors`, `secretion`, `two_component`, `defense`, `chaperones`, `ribosome`, `translation_factors`, `trna_biogenesis`, `dna_replication` |
| `tree_code` | str | Raw KEGG tree ID: `ko02000`, `ko01002`, ... |
| `level` | int | Set at ingest from JSON nesting depth. The top-level JSON object is tree metadata (not a node). A-level `children` = `0` (broadest); B-level = `1`; C-level = `2`; D-level = `3` (if non-KO). Mirrors the existing pattern: KEGG `category=0`, GO namespace roots = `0`. No synthetic root node, no post-import traversal. No `level_is_best_effort` flag (strict tree). |
| `level_kind` | str | Set at ingest from depth: `brite_class` (0), `brite_subclass` (1), `brite_family` (2), `brite_subfamily` (3). Matches the `compute_level_kind(depth)` utility. Trees that are shallower emit fewer distinct values. `brite_root` is dropped — A-level entries are already the roots. |
| `member_ko_count` | int | Post-import: recursive count of KO descendants |
| `gene_count` | int | Post-import: distinct genes reachable via KO leaves |
| `organism_count` | int | Post-import: distinct organisms among those genes |

### New edges

**`Kegg_term_in_brite_category`**: `(KeggTerm {level_kind:'ko'}) → (BriteCategory)` where the target is a BRITE leaf. A KO may have multiple outbound edges (one per tree it appears in). Enables polyhierarchy: the same KO participates in pathway + multiple BRITE trees simultaneously, mirroring KEGG's actual structure.

**`Brite_category_is_a_brite_category`**: `(BriteCategory) → (BriteCategory)` child → parent, within a single tree. Strict tree (no cross-links).

### Edge choice rationale

No direct `Gene → BriteCategory` edge. Gene reaches BRITE via `Gene_has_kegg_ko` → `Kegg_term_in_brite_category`. This avoids duplicating ~4,000 gene-facing edges per tree and keeps KO as the single source of truth for gene function assignment.

### Node label rationale

A distinct `BriteCategory` label (not overloading `KeggTerm` with new `level_kind` values) is chosen because:

- Pathway hierarchy (metabolic maps) and BRITE hierarchies (functional classification) have fundamentally different semantics. Unified `level: int` is meaningful within a single tree but not across pathway + BRITE.
- Query clarity: `(g)-[:Gene_has_kegg_ko]->(:KeggTerm)-[:Kegg_term_in_brite_category]->(b:BriteCategory {tree:'transporters'})` is more explicit than filtering on `level_kind`.
- Precedent: `Pfam` and `PfamClan` are already split along the same axis.

## Implementation

### Download utility: `multiomics_kg/utils/brite_utils.py`

Sits alongside `kegg_utils.py`. Reuses `_fetch_json` / `_fetch_text` from `kegg_utils.py` — no duplicate HTTP logic. Cache files written to `<cache_root>/kegg/brite_<tree_code>.json` (same `kegg/` subdirectory as `kegg_data.json`).

Functions:

- `download_brite_tree(tree_code: str, cache_root: Path, cache: bool = True) -> dict` — fetches `https://rest.kegg.jp/get/br:<tree_code>` via KEGG REST, parses flat format (A/B/C/D/E prefix lines where prefix char indicates depth), returns parsed hierarchy dict. Cache path: `<cache_root>/kegg/brite_<tree_code>.json`. HTML entity decoding for category names.
- `load_brite_trees(cache_root: Path, trees: list[str]) -> dict[str, dict]` — returns `{tree_code: parsed_hierarchy}` for all requested trees, fetching missing ones as needed.
- `compute_level_kind(depth: int) -> str` — generic depth → `level_kind` label: `0=brite_class`, `1=brite_subclass`, `2=brite_family`, `3=brite_subfamily`. Depth > 3 raises `ValueError` (not expected for the 12 configured trees).

Parser rules:
- `kegg_utils.py` already fetches `br:ko00001` in **JSON format** (`/get/br:ko00001/json`). `brite_utils.py` uses the same JSON endpoint for all 12 trees: `https://rest.kegg.jp/get/br:<tree_code>/json`. This avoids writing a flat-format parser — the recursive JSON `children` structure maps directly to the hierarchy.
- Depth is the nesting depth of the `children` array. Synthetic root = depth 0; first `children` level = depth 1 (= `brite_class`), etc.
- KO leaves are `children` entries whose `name` matches `K\d{5}`. They are not yielded as `BriteCategory` nodes; instead they trigger a `Kegg_term_in_brite_category` edge to their parent `BriteCategory`.
- Malformed entry (missing `name` key, unrecognised name format) is logged and skipped, not raised — BRITE data occasionally has annotation-only lines.

### Adapter: `multiomics_kg/adapters/brite_adapter.py`

Single class `MultiBriteAdapter(cache_root: Path, trees: list[str], test_mode: bool = False, cache: bool = True)`. Follows existing Multi* adapter pattern.

- `download_data(cache=True)` — calls `load_brite_trees` for all 12 configured trees.
- `get_nodes()` — yields `(node_id, 'brite category', properties)` for every BRITE category across all trees. `_clean_str` applied to `name`. Deduplicates by `id`.
- `get_edges()` — yields two edge streams:
  - `Brite_category_is_a_brite_category`: internal parent edges, child → parent, within each tree.
  - `Kegg_term_in_brite_category`: KO leaves (by `kegg.orthology:K<nnnnn>` ID) → BRITE leaf category.
- `test_mode=True` caps nodes at 100 per tree, for fast iteration.

### Schema config: `config/schema_config.yaml`

New entries:

- Node: `brite category` with all properties from the data model table.
- Edges: `brite category is a brite category`, `kegg term in brite category`.
- Full-text index entry: `briteCategoryFullText` on `name`.

### Pipeline wiring: `create_knowledge_graph.py`

Add `MultiBriteAdapter` instantiation after `MultiKeggAdapter`. Standard download → write_nodes → write_edges flow.

### Post-import: `scripts/post-import.sh` (and keep `.cypher` reference in sync)

Compute per-BriteCategory (recursive via descendant walk through `Brite_category_is_a_brite_category`):

- `member_ko_count`: count of KO descendants.
- `gene_count`: distinct genes reachable via KO leaves.
- `organism_count`: distinct `organism_name` among those genes.

Add scalar indexes:

- `brite_category_tree_idx` on `BriteCategory(tree)`
- `brite_category_level_idx` on `BriteCategory(level)`
- `brite_category_name_idx` on `BriteCategory(name)`

Full-text index:

- `briteCategoryFullText` on `BriteCategory(name)`.

## Testing

### Unit tests for utils: `tests/test_brite_utils.py`

- `download_brite_tree`: HTTP mocked — correct URL, correct cache path, parsed JSON round-trip. Cache hit path skips HTTP. `cache=False` forces re-fetch. Rate-limit sleep mocked and asserted between sequential fetches. Non-200 / empty response raises a clear error.
- Flat-format parser (pure function): handcrafted multi-line inputs cover depth and parent assignment; HTML entity decoding (`&amp;`, `&gt;`); KO leaf detection; malformed lines raise parse error with line number; empty input returns empty tree.
- `load_brite_trees`: returns dict keyed by tree_code; missing cache files trigger fetch; partial cache only fetches missing trees.
- `compute_level_kind`: depths 0–3 return expected labels (`brite_class`, `brite_subclass`, `brite_family`, `brite_subfamily`); depth ≥ 4 raises `ValueError`.

### Unit tests for adapter: `tests/test_brite_adapter.py`

- 3 handcrafted BRITE flat-format fixtures (Transporters, Peptidases, Two-component slices).
- `get_nodes()`: fixture-driven; verifies IDs, `tree`, `level`, `level_kind`.
- `get_edges()`: verifies one `Kegg_term_in_brite_category` edge per (KO, tree) leaf; verifies parent edges point child → parent.
- ID format: root node ID equals `normalize_curie(f'kegg.brite:{tree_code}')` (e.g. `kegg.brite:ko02000`); internal node IDs match `kegg.brite:{tree_code}.A{n}(.B{n}...)*`.
- ID stability: parsing the same fixture twice yields identical IDs; reordering an unrelated sibling node does not affect other nodes' IDs.
- String sanitization: names containing `'` → `^`, `|` stripped.

### KG validity tests: `tests/kg_validity/test_brite.py` (`@pytest.mark.kg`)

- All 12 trees present: `MATCH (b:BriteCategory) RETURN DISTINCT b.tree` returns the 11 expected names.
- Every tree has at least one `level=0` node (A-level entries); none have a parent edge.
- Every `BriteCategory` with `level > 0` has exactly one `Brite_category_is_a_brite_category` parent.
- Sanity floor: at least 200 KOs have outbound `Kegg_term_in_brite_category` edges.
- Orphan reduction spot-check:
  - K03798 (`ftsH`) reaches a peptidase tree root.
  - K06147 (ABC transporter) reaches a transporter tree root.
- Post-import properties populated: every `BriteCategory` has non-null `member_ko_count`, `gene_count`, `organism_count`.
- No duplicate `Kegg_term_in_brite_category` edges per (KO, BRITE leaf).

### Snapshot fixture update

Regenerate `tests/kg_validity/snapshot_data.json` via `tests/kg_validity/generate_snapshot.py` to include 3 sample `BriteCategory` nodes (one transporter, one peptidase, one two-component) and 3 sample `Kegg_term_in_brite_category` edges.

## Documentation

- `CLAUDE.md`:
  - "Actual Neo4j labels" list: add `BriteCategory`, `Brite_category_is_a_brite_category`, `Kegg_term_in_brite_category`.
  - "Key graph facts": BRITE subsection — tree list, approximate node and edge counts (populated after first build), level semantics, KO-polyhierarchy note.
  - "Post-import indexes" list: add the three new scalar indexes and `briteCategoryFullText`.
- `docs/kg-changes/brite-categories.md` (new): rationale (KEGG pathway-orphan KO analysis with numbers from this spec), tree list, data model diagram, example queries. Style mirrors `docs/kg-changes/ontology-level.md`.
- `memory/MEMORY.md`: add a pointer entry to the change doc.

## Appendix A — Suggestions for the explorer project (not implemented here)

These depend on the KG changes above and will be brainstormed / planned separately on the explorer side.

### Explorer framework extension: `ONTOLOGY_CONFIG` bridge_rel

`multiomics_explorer/multiomics_explorer/kg/queries_lib.py` drives all ontology traversal via `ONTOLOGY_CONFIG`. Each entry describes a one-hop gene→term connection:

```python
"kegg": {
    "label": "KeggTerm",
    "gene_rel": "Gene_has_kegg_ko",
    "hierarchy_rels": ["Kegg_term_is_a_kegg_term"],
    "gene_connects_to_level": "ko",  # filter descendants to ko-level before matching genes
}
```

BRITE requires a **two-hop bridge**: Gene → KeggTerm(ko) → BriteCategory(leaf). The current framework assumes `MATCH (g:Gene)-[:gene_rel]->(descendant)` — one hop. A new optional `bridge_rel` field is needed:

```python
"kegg_brite_transporters": {
    "label": "BriteCategory",
    "gene_rel": "Gene_has_kegg_ko",         # Gene → bridge (KeggTerm)
    "bridge_rel": "Kegg_term_in_brite_category",  # bridge → descendant (BriteCategory)
    "hierarchy_rels": ["Brite_category_is_a_brite_category"],
    "fulltext_index": "briteCategoryFullText",
    "tree_filter": "transporters",           # WHERE root.tree = 'transporters'
}
# ... repeat for each of the 12 trees
```

The generated Cypher pattern (in `_genes_by_ontology_cfg`) changes to:

```cypher
MATCH (root:BriteCategory) WHERE root.id = tid
MATCH (root)<-[:Brite_category_is_a_brite_category*0..15]-(descendant)
MATCH (bridge)-[:Kegg_term_in_brite_category]->(descendant)  -- bridge_rel
MATCH (g:Gene)-[:Gene_has_kegg_ko]->(bridge)                 -- gene_rel
```

This is a minimal two-line change to `_genes_by_ontology_cfg`: when `bridge_rel` is present, emit the bridge MATCH between the expansion and the gene MATCH; replace `descendant` with `bridge` in the gene MATCH.

### Other explorer work

1. **`enrich_functional` tool gains BRITE namespaces**. Add one `ONTOLOGY_CONFIG` entry per BRITE tree. `namespace=None` runs all 13 namespaces (pathway + 12 BRITE) independently — no cross-namespace FDR pooling — and returns a stacked result keyed by `(namespace, level)`. Add `min_members=3` to prune tiny BRITE families. Dedup at KO-leaf level: when `level` targets leaves, same gene set returns once tagged with all parent namespaces.
2. **New MCP tool `brite_categories_for_gene(locus_tag)`**. Returns the full BRITE classification path for a gene's KOs across all 12 trees.
3. **Cypher recipe additions for the `/cypher-queries` skill**: BRITE tree census per organism; family-level rollup for a gene; transporter subfamily enrichment in a DE gene set.
4. **Tool docs and about-content notes**: explain that pathway coverage is structurally partial in KEGG and BRITE trees complement pathway enrichment.

## Open questions

None at this time.
