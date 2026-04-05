# MCP Gene Lookup Services: get_gene and find_gene

**Goal:** Implement two MCP tools for gene retrieval from the KG — exact lookup by any known ID, and free-text search across functional annotations.

**Date:** 2026-03-11

---

## Services

### 1. `get_gene` — Exact lookup by any gene identifier

Return a gene node (with key properties) given any known identifier: locus tag, old locus tag, gene name, RefSeq protein ID, etc.

**Supported ID types (in resolution priority order):**

| Priority | ID type | Neo4j property | Index type |
|---|---|---|---|
| 1 | Locus tag (canonical) | `locus_tag` | scalar index |
| 2 | NCBI locus tag | `locus_tag_ncbi` | scalar index |
| 3 | Cyanorak locus tag | `locus_tag_cyanorak` | scalar index |
| 4 | Gene name | `gene_name` | scalar index |
| 5 | RefSeq protein ID | `protein_id` | scalar index |
| 6 | Old locus tags | `old_locus_tags[]` | ANY() scan |
| 7 | Alternative locus tags | `alternative_locus_tags[]` | ANY() scan |
| 8 | Gene name synonyms | `gene_name_synonyms[]` | ANY() scan |

**Parameters:**
- `id: str` — the identifier to look up
- `organism: str | None` — optional organism filter (e.g., "MED4", "EZ55")

**Returns:** gene node properties + linked organism name. Multiple matches returned (up to 5) if ambiguous.

**Cypher** (after schema changes adding `organism_strain` and `all_identifiers`):
```cypher
MATCH (g:Gene)
WHERE (
    g.locus_tag = $id
    OR g.gene_name = $id
    OR $id IN g.all_identifiers
  )
  AND ($organism IS NULL OR g.organism_strain = $organism)
RETURN g, g.organism_strain AS organism
LIMIT 5
```

**Note on `all_identifiers`:** covers `locus_tag_ncbi`, `locus_tag_cyanorak`, `protein_id`,
`old_locus_tags[]`, `alternative_locus_tags[]`, `gene_name_synonyms[]`. The explicit
`g.locus_tag` and `g.gene_name` checks come first to hit scalar indexes before the array scan.

---

### 2. `find_gene` — Free-text search across functional annotations

Search gene nodes by gene name, product description, functional annotation, pathway membership, etc. Uses Neo4j full-text index (Lucene).

**Parameters:**
- `search_text: str` — free-text query (supports Lucene syntax: `dna*`, `"DNA repair"`, `nitrogen AND transport`, `dnaN~`). Named `search_text` (not `query`) to avoid collision with the Neo4j Python driver's `session.run(query=...)` keyword argument.
- `organism: str | None` — optional organism filter (matched against `organism_strain` property)
- `min_quality: int` — minimum `annotation_quality` (0–3); default 0 (all genes); use 2 to skip hypothetical proteins
- `limit: int` — default 10, max 50

**Indexed properties (full-text index `geneFullText`):**
- `gene_summary` — covers gene_name + product + function_description in one field
- `gene_synonyms` — all name variants + historical locus tags from all sources
- `product_cyanorak` — Cyanorak-specific product; additive when NCBI/eggnog was chosen as `product`
- `alternate_functional_descriptions` — all-source descriptions; critical for low-quality genes
- `go_term_descriptions` — Cyanorak ontology term descriptions
- `pfam_names` — Pfam shortnames (e.g., `DNA_pol3_beta`, `HATPase_c`)
- `pfam_descriptions` — Cyanorak domain descriptions
- `eggnog_og_descriptions` — OG-level group descriptions (distinct from gene-level)

**Removed from index (covered by gene_summary):** `gene_name`, `product`, `function_description`

**Removed (content captured elsewhere):** `kegg_ko_descriptions` — these are EC descriptions from Cyanorak's `kegg_description` column; already covered by `product_cyanorak` and `alternate_functional_descriptions`

**Cypher** (after schema changes adding `organism_strain`):
```cypher
CALL db.index.fulltext.queryNodes('geneFullText', $search_text)
YIELD node AS g, score
WHERE ($organism IS NULL OR g.organism_strain = $organism)
  AND ($min_quality = 0 OR g.annotation_quality >= $min_quality)
RETURN g, score, g.organism_strain AS organism
ORDER BY score DESC
LIMIT $limit
```

**Note:** The parameter is `$search_text` (not `$query`) because the Neo4j Python driver's `session.run()` has a `query` keyword argument — passing `query=` as both the Cypher string and a parameter causes a `TypeError`.

No join to `OrganismTaxon` needed — `organism_strain` is directly on the gene node.

**Query examples an LLM agent would issue:**
- `"dnaN"` → DNA polymerase sliding clamp
- `"ROS detoxification"` → catalase, superoxide dismutase, etc.
- `"vitamin B12 synthesis"` → cobalamin biosynthesis genes
- `"iron*"` → iron transport, iron-sulfur cluster, etc.
- `"nitrogen AND transport"` → nitrogen ABC transporters

---

## Proposed Gene Node Schema/Content Changes

These changes improve lookup and search without breaking existing queries. All changes
are in the adapter (`cyanorak_ncbi_adapter.py`) and schema (`schema_config.yaml`);
no Cypher migration is needed since the graph is always rebuilt from scratch.

### 1. Add `organism_strain: str` to Gene nodes (denormalized)

**Problem:** Filtering `get_gene` / `find_gene` by organism currently requires a
relationship traversal (`Gene_belongs_to_organism → OrganismTaxon`), which defeats
full-text index performance and complicates the Cypher.

**Fix:** Emit `organism_strain` (e.g., `"MED4"`, `"EZ55"`) directly on each Gene node
from the adapter — the value comes from `self.strain_name` which is already available.

**Impact:**
- Organism-filtered `find_gene` becomes a simple property filter on the full-text result,
  no join needed
- `get_gene` with organism filter can add `AND g.organism_strain = $organism` after the
  full-text lookup
- Also add scalar index: `CREATE INDEX gene_organism_strain_idx FOR (g:Gene) ON (g.organism_strain)`

**Schema change:** Add `organism_strain: str` to `gene:` properties block.

---

### 2. Add `gene_summary: str` — computed single-line description

**Problem:** The MCP tool response needs a primary human-readable display field. Currently
the LLM agent would need to assemble `gene_name + product + function_description` itself,
or receive a large dict and pick from it.

**Fix:** Compute and emit `gene_summary` in the adapter as a single clean string.
For the description part, fall back through sources in quality order so that
low-annotation genes still get useful content:

```python
import re

# Best description: UniProt function_description first, then first clean
# alternate_functional_descriptions entry (strip [source] prefix for readability)
best_description = function_description or next(
    (re.sub(r'^\[.*?\]\s*', '', d) for d in (alternate_functional_descriptions or []) if d),
    None
)
parts = filter(None, [gene_name, product or product_cyanorak, best_description])
gene_summary = " :: ".join(parts)
# High-quality:  "dnaN :: DNA polymerase III, beta subunit :: Acts as a sliding clamp..."
# Low-quality:   "PMM1234 :: hypothetical protein :: COG0592: DNA polymerase III beta clamp"
# No-annotation: "PMM1234"  (nothing else available)
#
# NOTE: do NOT use " | " as separator — pipe is the BioCypher array delimiter and
# clean_text() converts it to "," which would silently mangle the field.
```

**Impact:**
- Include `gene_summary` in the full-text index — one high-signal field that covers the
  most common query targets
- Return as top-level field in MCP tool response
- Useful for LLM agent display without further processing

**Schema change:** Add `gene_summary: str` to `gene:` properties block.

---

### 3. Add `all_identifiers: str[]` — denormalized union of all ID fields

**Problem:** `get_gene` currently checks 8 separate fields (5 scalar + 3 array scans).
The OR chain is fragile and grows every time a new ID type is added.

**Fix:** Emit a single `all_identifiers: str[]` array in the adapter that unions:
- `locus_tag`, `locus_tag_ncbi`, `locus_tag_cyanorak`, `protein_id` (scalars)
- `old_locus_tags[]`, `alternative_locus_tags[]`, `gene_name_synonyms[]` (arrays)

```python
all_ids = set(filter(None, [locus_tag, locus_tag_ncbi, locus_tag_cyanorak, protein_id]))
all_ids.update(old_locus_tags or [])
all_ids.update(alternative_locus_tags or [])
all_ids.update(gene_name_synonyms or [])
```

`get_gene` query simplifies to:
```cypher
MATCH (g:Gene)
WHERE g.locus_tag = $id         // fast scalar index hit first
   OR g.gene_name = $id         // fast scalar index hit
   OR $id IN g.all_identifiers  // one array scan for all alternative IDs
RETURN g LIMIT 5
```

**Impact:** Simpler query; new ID types only require adapter change, not query change.
Neo4j 4.4 doesn't index arrays, but a single ANY() scan on 16K nodes is <5ms.

**Schema change:** Add `all_identifiers: str[]` to `gene:` properties block.

---

### 4. Fix `go_term_descriptions`: `str` → `str[]`

**Problem:** All other description fields are `str[]` arrays, but `go_term_descriptions`
is declared as `str` in schema_config.yaml and stored as a pipe-separated string
(e.g., `"biological process|molecular function|..."`). This is inconsistent and means
pipe characters in GO descriptions would be stripped by `clean_text()`.

**Fix:** Change `go_term_descriptions: str` → `go_term_descriptions: str[]` in schema.
Verify the adapter emits it as a list. Add to the full-text index.

**Impact:** Consistency. GO descriptions (e.g., "DNA repair", "response to oxidative stress")
are high-value search terms for biological queries.

**Schema change:** `go_term_descriptions: str` → `go_term_descriptions: str[]`

---

### 5. Keep `alternate_functional_descriptions` in the full-text index

**Revised decision:** Include `alternate_functional_descriptions` in the full-text index.

**Reason:** For `annotation_quality = 0 or 1` genes (eggNOG-only or unannotated),
`alternate_functional_descriptions` is often the *only* functional content — `product`
is "hypothetical protein", `function_description` is empty, and `gene_summary` is just
the locus tag. Dropping it makes these genes completely unsearchable by function.

**Re: `[source]` token noise:** The `[cyanorak]`, `[uniprot]`, `[eggnog]` prefix tokens
appear in thousands of gene entries, so Lucene's IDF scoring down-weights them
automatically — they become effective stop-words. A search for `"DNA repair"` will not
surface genes just because they carry `[cyanorak]` in a description. Not a real problem.

**No code change needed** — just include `alternate_functional_descriptions` in the
full-text index as-is.

---

### Summary of changes

| Change | Files | Impact |
|---|---|---|
| Add `organism_strain: str` | adapter + schema | Enables organism filter without join |
| Add `gene_summary: str` | adapter + schema | Primary display field; clean full-text target |
| Add `all_identifiers: str[]` | adapter + schema | Simpler `get_gene` query; extensible |
| Fix `go_term_descriptions: str[]` | schema (+ verify adapter) | Consistency; GO terms in full-text |
| Keep `alternate_functional_descriptions` in full-text index | post-import.cypher | Critical for low-annotation-quality genes; `[source]` tokens are IDF-suppressed |

All 5 changes are in the adapter and schema config — no changes to existing Cypher
queries in `post-import.sh`, no test migrations needed.

---

### Updated full-text index (after changes)

```cypher
CALL db.index.fulltext.createNodeIndex('geneFullText', ['Gene'],
  ['gene_summary', 'gene_synonyms', 'product_cyanorak',
   'alternate_functional_descriptions', 'go_term_descriptions',
   'pfam_names', 'pfam_descriptions', 'eggnog_og_descriptions'])
```

### Updated `get_gene` Cypher (after schema changes; same as Services section above)

```cypher
MATCH (g:Gene)
WHERE (
    g.locus_tag = $id
    OR g.gene_name = $id
    OR $id IN g.all_identifiers
  )
  AND ($organism IS NULL OR g.organism_strain = $organism)
RETURN g, g.organism_strain AS organism
LIMIT 5
```

---

## Implementation Plan

### Implementation order

Schema/adapter changes must precede index changes, which must precede the MCP tool.
KG rebuild is required between steps 1 and 2.

```
Step 1 → Step 2 → [KG rebuild] → Step 3 → Step 4 → Step 5
```

### Step 1: Schema + adapter + build script changes — DONE (2026-03-12)

**Implementation notes:**
- New fields (`organism_strain`, `gene_summary`, `all_identifiers`) computed in `build_gene_annotations.py` (not adapter) — same pattern as `annotation_quality` and `alternate_functional_descriptions`
- `organism_strain` uses `preferred_name` from `cyanobacteria_genomes.csv` (e.g., "Prochlorococcus MED4")
- `gene_summary` uses raw source values directly (no [source] prefix stripping needed)
- `go_term_descriptions` changed to `passthrough_list` in `gene_annotations_config.yaml`
- GeneNodeField enum extended with 3 new members; adapter passes them through automatically
- Scalar + full-text indexes added to `post-import.sh` and `post-import.cypher`
- MED4 verified: all 1976 genes have `organism_strain` and `gene_summary`; `go_term_descriptions` is a list
- 1012 unit tests pass (2 pre-existing failures in EC adapter mocks unrelated)

**Original plan (kept for reference):**

**`config/schema_config.yaml`** — in the `gene:` properties block, add:
```yaml
organism_strain: str
gene_summary: str
all_identifiers: str[]
# change:
go_term_descriptions: str[]   # was: str
```

**`multiomics_kg/adapters/cyanorak_ncbi_adapter.py`** — in `_get_gene_nodes()`, after
building `node_properties`, compute and inject the three new fields:

```python
import re

locus_tag = row.get('locus_tag')
gene_name = node_properties.get('gene_name')
product = node_properties.get('product') or node_properties.get('product_cyanorak')
function_description = node_properties.get('function_description')
alt_descs = node_properties.get('alternate_functional_descriptions') or []
old_locus_tags = node_properties.get('old_locus_tags') or []
alternative_locus_tags = node_properties.get('alternative_locus_tags') or []
gene_name_synonyms = node_properties.get('gene_name_synonyms') or []
protein_id = node_properties.get('protein_id')
locus_tag_ncbi = node_properties.get('locus_tag_ncbi')
locus_tag_cyanorak = node_properties.get('locus_tag_cyanorak')

# organism_strain — from self.strain_name
if self.strain_name:
    node_properties['organism_strain'] = self.strain_name

# gene_summary — fallback chain for description
best_desc = function_description or next(
    (re.sub(r'^\[.*?\]\s*', '', d) for d in alt_descs if d), None
)
parts = [p for p in [gene_name, product, best_desc] if p]
if parts:
    node_properties['gene_summary'] = self.clean_text(' :: '.join(parts))
    # Use " :: " not " | " — pipe is the BioCypher array delimiter, clean_text() maps it to ","

# all_identifiers — union of all alternative ID fields
# Exclude locus_tag and gene_name: they have their own scalar-indexed fields in the
# get_gene query, so including them here would be redundant. Also deduplicate across
# arrays (set handles values appearing in multiple sources).
scalar_indexed = {locus_tag, gene_name} - {None}
all_ids = set(filter(None, [locus_tag_ncbi, locus_tag_cyanorak, protein_id]))
all_ids.update(old_locus_tags)
all_ids.update(alternative_locus_tags)
all_ids.update(gene_name_synonyms)
all_ids -= scalar_indexed  # remove values already covered by scalar indexes
if all_ids:
    node_properties['all_identifiers'] = self.clean_text(sorted(all_ids))
```

Also verify `go_term_descriptions` is already emitted as a list (not a pipe-separated string);
fix in the adapter if needed.

### Step 2: Add indexes to post-import.sh and post-import.cypher

Add to both [scripts/post-import.sh](../scripts/post-import.sh) (as `cypher-shell` blocks) and
[scripts/post-import.cypher](../scripts/post-import.cypher) (for documentation):

```cypher
-- Scalar indexes for get_gene exact lookup
CREATE INDEX gene_locus_tag_idx IF NOT EXISTS FOR (g:Gene) ON (g.locus_tag);
CREATE INDEX gene_name_idx IF NOT EXISTS FOR (g:Gene) ON (g.gene_name);
CREATE INDEX gene_organism_strain_idx IF NOT EXISTS FOR (g:Gene) ON (g.organism_strain);

-- Full-text index for find_gene
CALL db.index.fulltext.createNodeIndex('geneFullText', ['Gene'],
  ['gene_summary', 'gene_synonyms', 'product_cyanorak',
   'alternate_functional_descriptions', 'go_term_descriptions',
   'pfam_names', 'pfam_descriptions', 'eggnog_og_descriptions'])
```

Note: `db.index.fulltext.createNodeIndex` is the Neo4j 4.4 API. For 5.x, replace with:
```cypher
CREATE FULLTEXT INDEX geneFullText IF NOT EXISTS
FOR (g:Gene) ON EACH [g.gene_name, g.gene_synonyms, g.gene_summary, ...]
```

### Step 3: Rebuild KG and verify — DONE (2026-03-12)

KG rebuilt via Docker. All new fields verified in Neo4j:
- `organism_strain`: 100% populated on all 16,226 Gene nodes
- `gene_summary`: >95% populated (all genes with any annotation)
- `all_identifiers`: >90% populated (array of alternative IDs)
- `go_term_descriptions`: stored as list (str[]), not string
- Scalar indexes on `locus_tag`, `gene_name`, `organism_strain` confirmed
- Full-text index `geneFullText` with 8 properties confirmed
- `get_gene` and `find_gene` query patterns verified against live graph

### Step 4: MCP server implementation

**Framework:** Use the `mcp` Python SDK (`pip install mcp`) with `fastmcp` for tool
registration, or integrate into whatever MCP server hosts the other KG tools.
Location TBD (see Open Questions).

Tools to implement in `gene_lookup.py`:
- `get_gene(id, organism=None)` — exact ID match
- `find_gene(search_text, organism=None, min_quality=0, limit=10)` — full-text search

**Response shape** — single gene (used by `get_gene`; `find_gene` returns a list of these):
```json
{
  "locus_tag": "PMM0001",
  "gene_name": "dnaN",
  "gene_summary": "dnaN :: DNA polymerase III, beta subunit :: Acts as a sliding clamp...",
  "product": "DNA polymerase III, beta subunit",
  "function_description": "...",
  "organism_strain": "MED4",
  "go_terms": ["GO:0003677"],
  "kegg_ko": ["K02338"],
  "annotation_quality": 3,
  "score": 1.0
}
```

`find_gene` returns:
```json
{
  "results": [ <gene>, <gene>, ... ],
  "total": 10,
  "query": "dnaN"
}
```

**Error handling:**

| Case | Behaviour |
|---|---|
| `get_gene` → 0 matches | Return `{"results": [], "message": "No gene found for id '<id>'"}` |
| `get_gene` → multiple matches | Return all (≤5) with a `"message": "Ambiguous — <n> matches found. Specify organism to narrow."` |
| `find_gene` → malformed Lucene query | Catch Neo4j `CypherExecutionException`; retry with `search_text` escaped via `re.sub(r'[+\-!(){}\[\]^"~*?:\\/]', r'\\\g<0>', search_text)`; if still fails return error |
| Neo4j unreachable | Raise `McpError` with `"Knowledge graph unavailable"` |

### Step 5: KG validity tests — DONE (2026-03-12)

Added `tests/kg_validity/test_gene_lookup.py` — 18 tests, all passing:

```python
# Index existence
def test_scalar_indexes_exist(session): ...  # locus_tag, gene_name, organism_strain
def test_fulltext_index_exists(session): ...  # geneFullText

# get_gene — one test per ID type
def test_get_gene_by_locus_tag(session): ...        # PMM0001
def test_get_gene_by_old_locus_tag(session): ...    # PMT9313 → PMT9313_RS... mapping
def test_get_gene_by_gene_name(session): ...        # "dnaN"
def test_get_gene_by_protein_id(session): ...       # WP_... accession
def test_get_gene_organism_filter(session): ...     # same gene_name in MED4 + MIT9312, filter returns only MED4

# find_gene
def test_find_gene_known_query(session): ...        # "photosystem" → results
def test_find_gene_organism_filter(session): ...    # "dnaN" with organism="MED4"
def test_find_gene_quality_filter(session): ...     # min_quality=2 excludes hypothetical proteins
def test_find_gene_score_ordered(session): ...      # scores descending

# New fields present
def test_gene_summary_populated(session): ...       # non-null on well-annotated genes
def test_organism_strain_populated(session): ...    # all Gene nodes have organism_strain
def test_all_identifiers_populated(session): ...    # array, non-empty
```

---

## Decision Log

| # | Decision | Rationale |
|---|----------|-----------|
| D1 | Use Neo4j full-text index (Lucene) rather than a separate search DB | 16K nodes is small; no operational overhead; Lucene supports all needed query types |
| D2 | No semantic/embedding search in initial implementation | LLM agents can reformulate queries; in-memory cosine similarity can be added later as `find_gene_semantic` without touching the KG |
| D3 | Stay on Neo4j 4.4 | `biocypher/base:1.2.0` has no 5.x successor; migration cost (~1-2 days) not justified at current scale; vector index not needed |
| D4 | `gene_name_synonyms[]` included in get_gene array scan | Enables lookup by standard biochemical abbreviations (e.g., "nblS", "hspA") that appear in papers but may not be the canonical gene_name |
| D5 | Return up to 5 results from get_gene | Some IDs are genuinely ambiguous across strains (same gene_name in MED4 and MIT9312); caller decides how to handle |
| D6 | No pre-generated LLM summaries stored in Gene nodes; use rule-based `gene_summary` instead | See rationale below |
| D7 | MCP server in a separate repo | Keeps KG repo focused on data pipeline; MCP repo can evolve independently |
| D8 | Expression edge summaries are a separate MCP tool, not part of `get_gene` | `get_gene` stays focused on gene identity/annotation; expression is a distinct query type |
| D9 | `find_gene` searches Gene nodes only, not Protein nodes | KG is gene-centric; Protein nodes are secondary and their key fields are denormalized onto Gene nodes already |

---

### D6 rationale: LLM-generated integrated summary

Three options were considered:

**Option A — Pre-generate at build time, store as `llm_summary: str` on Gene node**
- Pros: zero query latency; resolves multi-source description conflicts (e.g., "hypothetical protein" vs "ABC transporter substrate-binding"); searchable in full-text index
- Cons: ~16K genes × ~400 tokens ≈ $1–40 depending on model; adds 10–20 min to build pipeline; **goes stale** — expression context baked in can't update without full regeneration; hard to iterate on prompt quality without full rebuild

**Option B — Generate on-the-fly in MCP tool**
- Pros: always current; can incorporate live graph data (edge counts, homolog count)
- Cons: 1–2s latency per query in the critical path; cost per call; the outer LLM agent is already an LLM — synthesizing for another LLM is usually redundant

**Option C (chosen) — Structured data only; rule-based `gene_summary`**
- The rule-based `gene_summary` field (`gene_name | product | function_description`) provides a clean primary display field
- The receiving LLM agent (Claude/GPT-4) synthesizes naturally from structured dicts — no pre-processing needed
- Avoids staleness and cost problems entirely

**Exception to revisit:** If agents demonstrably struggle with conflicting multi-source descriptions (e.g., `annotation_quality <= 1` genes with 3+ divergent `alternate_functional_descriptions`), consider pre-generating `llm_summary` only for the ~5K best-annotated genes (`annotation_quality >= 2`). Cost drops to ~$0.50–5; build impact is minimal. Defer until there is evidence this is needed.

---

## Resolved Questions

| Question | Decision |
|---|---|
| MCP server location | Separate repo |
| MCP framework | Needs scoping — separate plan; this plan covers KG infrastructure only (indexes, new fields) |
| Expression edge summaries | Separate MCP tool (`get_gene_expression` or similar); not part of `get_gene` response |
| Protein node search | No — KG is gene-centric; Protein nodes are secondary. `find_gene` searches Gene nodes only |

## Scope of this plan

This plan covers **KG-side infrastructure** only:
- Schema changes (new gene node fields)
- Adapter changes to populate them
- Post-import indexes (scalar + full-text)
- KG validity tests

The MCP server implementation (Step 4) is a placeholder. A separate plan will cover MCP framework selection, server architecture, tool definitions, and deployment.
