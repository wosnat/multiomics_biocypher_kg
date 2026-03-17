# Pfam Normalization + Graph Nodes

**Status: COMPLETE** (2026-03-17)
- Steps 1-8 implemented and verified
- KG rebuilt and deployed with Pfam/PfamClan nodes
- 357 KG validity tests pass
- Schema fix applied: pfam clan uses `input_label` override for BioCypher compatibility

## Background

**Pfam** (Protein Families) is a database of ~27.5K conserved protein domain families maintained at EMBL-EBI (now part of InterPro). Each entry represents a structural/functional unit within proteins — e.g., `PF00712` (DNA_pol3_beta) is the DNA polymerase III beta subunit domain. Pfam domains are the standard way to characterize protein function by structure, complementary to:
- **GO** (function/process/location ontology)
- **EC** (enzymatic reaction classification)
- **KEGG** (metabolic pathway membership)
- **COG** (evolutionary conservation groups)

Pfam adds a **structural domain perspective** that none of the above provide: which modular structural units does a protein contain, and which clan (superfamily) do they belong to?

## Problem (resolved)

Gene nodes had three Pfam properties (`pfam_ids`, `pfam_names`, `pfam_descriptions`) sourced from three independent databases. The lists were misaligned. All three properties have been dropped from Gene nodes and replaced by `Pfam`/`PfamClan` graph nodes with `Gene_has_pfam` and `Pfam_in_pfam_clan` edges.

### State in the KG before this task (2026-03-17)

| Metric | Count |
|--------|-------|
| Total genes | 35,226 |
| Has `pfam_ids` (PF* accessions) | 18,055 |
| Has `pfam_names` (shortnames) | 26,607 |
| Has `pfam_descriptions` | 14,910 |
| **Has names but NO ids** | **9,126** |
| Has ids but NO names | 574 |

### Per-strain pfam_ids coverage

| Strain | Genes with pfam_ids | Notes |
|--------|-------------------|-------|
| Synechococcus CC9311 | 2,015 | |
| Synechococcus WH8102 | 1,997 | |
| Prochlorococcus MIT9313 | 1,814 | |
| Prochlorococcus NATL1A | 1,509 | |
| Prochlorococcus NATL2A | 1,492 | |
| Prochlorococcus MIT9312 | 1,478 | |
| Prochlorococcus MED4 | 1,467 | |
| Prochlorococcus AS9601 | 1,465 | |
| Prochlorococcus MIT9301 | 1,455 | |
| Alteromonas HOT1A3 | 1,262 | No Cyanorak data |
| Alteromonas EZ55 | 1,116 | No Cyanorak data |
| Alteromonas MIT1002 | 985 | No Cyanorak data |
| Prochlorococcus RSP50 | **0** | Missing — needs investigation |

### Data sources and their raw values

From `gene_annotations_wide.json` (PMM0001 = dnaN, DNA polymerase III beta subunit):

| Source | Field | Raw value |
|--------|-------|-----------|
| **Cyanorak** (gene_mapping) | `protein_domains` | `TIGR00663,PF00712,PF02768,PF02767,IPR022634,IPR022635,IPR022637,IPR001001` |
| **Cyanorak** (gene_mapping) | `protein_domains_description` | 8 mixed descriptions (TIGR + Pfam + InterPro) |
| **EggNOG** | `PFAMs` | `DNA_pol3_beta,DNA_pol3_beta_2,DNA_pol3_beta_3` |
| **UniProt** | `pfam_ids` | `['PF00712', 'PF02767', 'PF02768']` |

After merging via `gene_annotations_config.yaml`:
- `pfam_ids`: `['PF00712', 'PF02768', 'PF02767']` — union of Cyanorak (PF* tokens only via `extract_pfam_ids`) + UniProt + eggNOG (PF* tokens only). **EggNOG shortnames (non-PF* tokens) are not captured** in `pfam_ids`.
- `pfam_names`: `['DNA_pol3_beta', 'DNA_pol3_beta_2', 'DNA_pol3_beta_3']` — shortnames from eggNOG only. No way to pair them with IDs.
- `pfam_descriptions`: 8 descriptions — **all** domain descriptions from `protein_domains_description`, unfiltered (includes TIGR, InterPro, etc.). Counts never match pfam_ids.

### The eggNOG-only gap (was)

PMM0002 (photosynthetic reaction center protein): eggNOG said it had Pfam domain `PRC`, but `pfam_ids` was null because there was no Cyanorak or UniProt PF* ID for this gene. The domain info existed but was invisible to any ID-based query. Now resolved by shortname-to-accession lookup via Pfam reference data.

**Alteromonas was worst hit** (no Cyanorak genome data):
- HOT1A3: 3,511 pfam_names vs 1,262 pfam_ids — **2,249 genes with domain info but no accession**
- EZ55: 3,605 pfam_names vs 1,116 pfam_ids — **2,489 genes with domain info but no accession**

### Impact on full-text search (was)

`post-import.sh` indexed `pfam_names` in `geneFullText`. Now replaced by `pfamFullText` on `Pfam` nodes (covers `name` and `short_name`).

## ROI

### What this enables (not currently possible)

1. **Domain-based gene search**: Find all genes with a specific structural domain across all strains
2. **Domain architecture comparison**: Compare domain repertoires between strains/ecotypes
3. **Clan-level grouping**: Find all genes in a superfamily (e.g., "DNA clamp" clan CL0060)
4. **Cross-strain domain conservation**: Which domains are universal vs strain-specific?
5. **Domain co-occurrence**: Which domains tend to appear together in Prochlorococcus proteins?

### What this fixed

- Recovered many of the ~9,126 genes that had domain names but no accessions (eggNOG shortname to PF* ID resolution via Pfam reference data)
- Eliminated `pfam_names`, `pfam_descriptions`, and the broken `extract_pfam_ids`/`extract_pfam_names` transforms -- replaced by a single raw union + post-merge enrichment
- Dropped all pfam properties from Gene nodes -- fully represented as graph edges
- Aligned pfam data into a clean, queryable graph structure matching the GO/EC/KEGG/COG pattern

## Solution

Two-part approach:

### Part 1: Normalize pfam fields in prepare_data pipeline (step 2)

Download `Pfam-A.clans.tsv.gz` (494 KB, ~27.5K entries) as reference:

```
PF00712  CL0060  DNA_clamp  DNA_pol3_beta  DNA polymerase III beta subunit, N-terminal domain
```

Columns: `accession | clan_accession | clan_name | shortname | description`

Use as single source of truth. Resolve eggNOG shortnames → PF* IDs via reverse lookup. Look up names/descriptions by accession. Result: aligned `pfam_ids` in `gene_annotations_merged.json`.

### Part 2: Pfam + PfamClan nodes in functional_annotation_adapter.py

Follow the established pattern (GO, EC, KEGG, COG). **Two node types** with correct bioregistry CURIEs:

| Component | Details |
|---|---|
| **Pfam nodes** (domains) | ~2K entries (only those referenced by genes), Neo4j label `Pfam`, `preferred_id: pfam`, node ID `pfam:PF00712`, properties: `name` (description), `short_name` (Pfam shortname) |
| **PfamClan nodes** (superfamilies) | ~800 clans (only those referenced by emitted domains), Neo4j label `PfamClan`, `preferred_id: pfam.clan`, node ID `pfam.clan:CL0060`, properties: `name` (clan name) |
| **Gene_has_pfam edges** | Gene → Pfam, from `pfam_ids` in gene annotations |
| **Pfam_in_pfam_clan edges** | Pfam → PfamClan (flat, one level) |

Gene node Pfam properties were **dropped** -- `pfam_names` and `pfam_descriptions` were eliminated at the config level (replaced by a single raw `pfam_ids` union), and `pfam_ids` was removed from Gene nodes (replaced by `Gene_has_pfam` edges). The `geneFullText` index was updated to remove `pfam_names` (covered by `pfamFullText`).

## Design

### New utility: `multiomics_kg/utils/pfam_utils.py`

Follow the `go_utils.py` pattern (lives in `multiomics_kg/utils/`, not `download/utils/`):

```python
PFAM_CLANS_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz"

@dataclass
class PfamEntry:
    accession: str        # PF00712
    shortname: str        # DNA_pol3_beta
    description: str      # DNA polymerase III beta subunit, N-terminal domain
    clan_accession: str   # CL0060 or ""
    clan_name: str        # DNA_clamp or ""

@dataclass
class PfamData:
    by_accession: dict[str, PfamEntry]     # PF00712 → PfamEntry
    by_shortname: dict[str, str]           # DNA_pol3_beta → PF00712 (reverse lookup)
    clans: dict[str, str]                  # CL0060 → "DNA_clamp"

def load_pfam_data(cache_root: Path, force: bool = False) -> PfamData:
    """Load or download+cache Pfam reference.

    Cache: cache_root/pfam/pfam_reference.json
    Downloads Pfam-A.clans.tsv.gz, parses, caches as JSON, deletes TSV.
    """
```

### Changes to `gene_annotations_config.yaml`

Replace all three pfam field definitions (`pfam_ids`, `pfam_names`, `pfam_descriptions`) with a single raw union:

```yaml
pfam_ids:
  # Post-merge enriched: raw tokens resolved to clean PF* IDs via Pfam reference.
  # The union collects ALL tokens from all sources (PF* IDs, shortnames, TIGR*, IPR*).
  # enrich_pfam_fields() resolves shortnames and filters to PF* only.
  type: union
  sources:
    - source: gene_mapping
      field: protein_domains
      delimiter: ","
    - source: eggnog
      field: PFAMs
      delimiter: ","
    - source: uniprot
      field: pfam_ids
```

No per-source transforms — all filtering and resolution happens in the post-merge enrichment step.

### Cleanup in `annotation_transforms.py` and `build_gene_annotations.py`

- Remove `_tx_extract_pfam_ids` and `_tx_extract_pfam_names` from `annotation_transforms.py`
- Remove their entries from the `TRANSFORMS` registry dict
- Remove the special-case `extract_pfam_ids` / `extract_pfam_names` branches in `_resolve_union()`

### Changes to `build_gene_annotations.py` (step 2)

Add post-merge enrichment (after `extract_ortholog_groups`, before writing JSON):

```python
def enrich_pfam_fields(gene: dict, pfam_data: PfamData) -> None:
    """Resolve raw pfam_ids tokens to clean PF* accessions.

    pfam_ids at this point is a raw union of all tokens from Cyanorak protein_domains,
    eggNOG PFAMs, and UniProt pfam_ids — a mix of PF* IDs, shortnames, TIGR*, IPR*.

    1. For each token: PF* → keep; in pfam_data.by_shortname → resolve to PF* ID; else drop
    2. Deduplicate
    3. Overwrite pfam_ids with clean PF* list (or delete if empty)
    4. Return unresolved shortnames for logging
    """
```

### Changes to `cyanorak_ncbi_adapter.py`

Remove `pfam_ids`, `pfam_names`, `pfam_descriptions` from the Gene node properties dict emitted by the adapter. (`pfam_names` and `pfam_descriptions` no longer exist in merged JSON; `pfam_ids` is consumed by the adapter, not written to Gene nodes.)

### Full-text index update for Gene nodes

Update `geneFullText` in both `post-import.sh` and `post-import.cypher`:

```cypher
-- BEFORE (current):
CREATE FULLTEXT INDEX geneFullText IF NOT EXISTS FOR (n:Gene) ON EACH [
  n.gene_summary, n.gene_synonyms,
  n.alternate_functional_descriptions, n.pfam_names];

-- AFTER (pfam_names removed, covered by pfamFullText):
CREATE FULLTEXT INDEX geneFullText IF NOT EXISTS FOR (n:Gene) ON EACH [
  n.gene_summary, n.gene_synonyms,
  n.alternate_functional_descriptions];
```

### New adapter: `MultiPfamAnnotationAdapter` in `functional_annotation_adapter.py`

Follows `MultiKeggAnnotationAdapter` pattern. **Two node types** with separate bioregistry prefixes:

```python
class PfamAnnotationAdapter:
    """Per-strain: reads gene_annotations_merged.json, yields gene→Pfam edges."""

    def get_all_pfam_ids(self) -> set[str]:
        """All PF* IDs referenced by any gene in this strain."""

    def get_edges(self):
        """Yield (edge_id, gene_node_id, pfam_node_id, 'gene_has_pfam', {})"""


class MultiPfamAnnotationAdapter:
    """Multi-strain: owns Pfam + PfamClan nodes and all edges."""

    def __init__(self, genome_config_file, cache_root, test_mode, cache):
        self.pfam_data = load_pfam_data(cache_root, force=not cache)
        # ... build per-strain adapters

    def get_nodes(self):
        """Yield two node types:

        1. Pfam (domain) nodes — only ~2K referenced by genes
           - Node ID: normalize_curie("pfam:PF00712")
           - label: "pfam"
           - name: "DNA polymerase III beta subunit, N-terminal domain"
           - short_name: "DNA_pol3_beta"
        2. PfamClan nodes — only clans referenced by emitted domains
           - Node ID: normalize_curie("pfam.clan:CL0060")
           - label: "pfam clan"
           - name: "DNA_clamp"
        """

    def get_edges(self):
        """Yield:
        1. gene → Pfam edges (gene_has_pfam)
        2. Pfam → PfamClan edges (pfam_in_pfam_clan)
        """
```

### Schema additions (`schema_config.yaml`)

Two node types with correct bioregistry CURIEs:

```yaml
# Pfam domain families
pfam:
  is_a: named thing
  represented_as: node
  preferred_id: pfam
  label_in_input: pfam
  properties:
    name: str            # "DNA polymerase III beta subunit, N-terminal domain"
    short_name: str      # "DNA_pol3_beta"

# Pfam clan superfamilies
pfam clan:
  is_a: named thing
  represented_as: node
  preferred_id: pfam.clan
  label_in_input: pfam_clan
  properties:
    name: str            # "DNA_clamp"

# Edges
gene to pfam association:
  is_a: association
  represented_as: edge
  label_as_edge: Gene_has_pfam
  source: gene
  target: pfam
  label_in_input: gene_has_pfam

pfam in pfam clan:
  is_a: association
  represented_as: edge
  label_as_edge: Pfam_in_pfam_clan
  source: pfam
  target: pfam clan
  label_in_input: pfam_in_pfam_clan
```

### Integration into `create_knowledge_graph.py`

Add after the existing functional annotation adapters:

```python
pfam_adapter = MultiPfamAnnotationAdapter(
    genome_config_file=genome_config_file,
    cache_root=Path("cache/data"),
    test_mode=test_mode,
    cache=cache,
)
bc.write_nodes(pfam_adapter.get_nodes())
bc.write_edges(pfam_adapter.get_edges())
```

### Post-import indexes (`scripts/post-import.sh`)

```cypher
-- Pfam domain indexes
CREATE INDEX pfam_name_idx IF NOT EXISTS FOR (p:Pfam) ON (p.name);

-- PfamClan indexes
CREATE INDEX pfam_clan_name_idx IF NOT EXISTS FOR (c:PfamClan) ON (c.name);

-- Full-text indexes for domain/clan name search
CREATE FULLTEXT INDEX pfamFullText IF NOT EXISTS
  FOR (p:Pfam) ON EACH [p.name, p.short_name];
CREATE FULLTEXT INDEX pfamClanFullText IF NOT EXISTS
  FOR (c:PfamClan) ON EACH [c.name];
```

## Edge Cases

### Obsolete/renamed eggNOG shortnames

EggNOG v2.1.13 uses Pfam 33.1 (2020); current `Pfam-A.clans.tsv` is Pfam 37+ (2024). Some shortnames may have been renamed or retired between versions.

**Mitigation**: `enrich_pfam_fields()` must log unresolved shortnames with counts per strain. Output a summary to the build log, e.g.:
```
[pfam] MED4: resolved 1423/1445 eggNOG shortnames; 22 unresolved: PRC (15 genes), FtsH_ext (3 genes), ...
```

Unresolved shortnames are simply not added to `pfam_ids`. This is acceptable — a few percent loss from version skew is better than the current 34% gap.

### Shortname uniqueness

Pfam shortnames are unique within Pfam-A — each PF accession has exactly one shortname and vice versa. The reverse lookup (`by_shortname`) is safe; no collision handling needed.

### Domains without clan membership

Not all Pfam domains belong to a clan (~45% of the 27.5K entries have no clan). Domains without a clan still get Pfam nodes and gene→domain edges; they just have no `Pfam_in_pfam_clan` edge. The adapter handles this by checking `entry.clan_accession` before emitting clan-related entities.

### RSP50 missing pfam_ids

RSP50 shows 0 genes with `pfam_ids` in the live KG. Needs investigation:
- Does RSP50 have eggnog annotations? (check `cache/data/Prochlorococcus/genomes/RSP50/eggnog/`)
- Does RSP50 have Cyanorak protein_domains? (check `gene_mapping.csv`)
- Is RSP50 included in UniProt download? (check taxid coverage)

This is a pre-existing data gap, not caused by this plan — but the enrichment step should recover RSP50 coverage if eggnog annotations exist.

## Implementation Steps

### Step 1: `pfam_utils.py`
- Create `multiomics_kg/utils/pfam_utils.py`
- Download, parse, cache `Pfam-A.clans.tsv.gz`
- Expose `load_pfam_data()` returning `PfamData`
- Cache at `cache/data/pfam/pfam_reference.json`
- Unit test: verify reference loading, shortname→accession lookup, cache round-trip

### Step 2: Simplify `gene_annotations_config.yaml` + remove old transforms
- Replace `pfam_ids`, `pfam_names`, `pfam_descriptions` with single raw `pfam_ids` union (no transforms)
- Remove `_tx_extract_pfam_ids` and `_tx_extract_pfam_names` from `annotation_transforms.py`
- Remove their special-case branches in `_resolve_union()` in `build_gene_annotations.py`

### Step 3: Post-merge enrichment in `build_gene_annotations.py`
- Add `enrich_pfam_fields()` — resolves raw token list to clean PF* IDs via reference
- Log unresolved shortnames with counts per strain
- Unit test: verify enrichment with mock data (PF* passthrough, shortname resolution, TIGR/IPR filtering, dedup)

### CHECKPOINT: Verify gene_annotations_merged.json
- Run `prepare_data.sh --steps 2` for at least MED4 and one Alteromonas strain
- Spot-check `pfam_ids` in merged JSON: should be clean PF* lists only (no TIGR*, IPR*, shortnames)
- Verify PMM0001 has `['PF00712', 'PF02767', 'PF02768']` (or superset if shortnames resolved new IDs)
- Verify genes that previously had `pfam_names` but no `pfam_ids` now have PF* IDs (e.g., PMM0002)
- Verify `pfam_names` and `pfam_descriptions` do NOT exist in merged JSON
- Check unresolved shortname log output — should be a small percentage
- **Stop and review before proceeding to adapter/graph changes**

### Step 4: Drop pfam properties from Gene nodes
- Remove `pfam_ids`, `pfam_names`, `pfam_descriptions` from `cyanorak_ncbi_adapter.py` Gene node properties (`pfam_names`/`pfam_descriptions` no longer in merged JSON; `pfam_ids` consumed by adapter not Gene)
- Remove `pfam_ids`, `pfam_names`, `pfam_descriptions` from Gene properties in `schema_config.yaml`
- Update `geneFullText` index in `post-import.sh` and `post-import.cypher` to remove `pfam_names`
- Update `CLAUDE.md` Gene node property list

### Step 5: Adapter in `functional_annotation_adapter.py`
- `PfamAnnotationAdapter` (per-strain)
- `MultiPfamAnnotationAdapter` (multi-strain, nodes + edges)
- Two node types: `Pfam` (`preferred_id: pfam`) and `PfamClan` (`preferred_id: pfam.clan`)
- Pfam: `name` (description), `short_name` (Pfam shortname)
- PfamClan: `name` (clan name)
- Unit test: verify adapter yields correct nodes/edges for sample data

### Step 6: Schema + pipeline wiring
- Add `pfam` node type, `pfam clan` node type, two edge types to `schema_config.yaml`
- Add adapter to `create_knowledge_graph.py`
- Add indexes to `post-import.sh` and `post-import.cypher`
- Add `pfamFullText` + `pfamClanFullText` full-text indexes

### Step 7: Test + validate
- Run `prepare_data.sh --steps 2` to regenerate `gene_annotations_merged.json`
- Run `create_knowledge_graph.py --test` to verify adapter output
- Full KG build + verify counts (see Verification section)
- Check RSP50 pfam coverage after enrichment
- Review unresolved shortname logs

### Step 8: Documentation
- Update `CLAUDE.md`: add Pfam to "Actual Neo4j labels" (nodes: `Pfam`, `PfamClan`; edges: `Gene_has_pfam`, `Pfam_in_pfam_clan`), post-import indexes, key graph facts, drop pfam properties from Gene node list
- Update `.claude/skills/cypher-queries/SKILL.md` with Pfam query examples
- KG change doc for explorer consumption: `plans/kg_changes_for_pfam.md` (already written)
- Update memory files if needed

## Expected graph output

| Entity | Estimated count |
|---|---|
| Pfam nodes (domains) | ~2,000 (unique PF* IDs across all 13 strains) |
| PfamClan nodes (superfamilies) | ~800 (clans referenced by those domains; 812 total in reference) |
| Gene_has_pfam edges | ~25,000 (genes x avg ~1.5 domains) |
| Pfam_in_pfam_clan edges | ~1,500 (domains that belong to a clan — not all do) |

## Example queries enabled

```cypher
-- Find all genes with a specific Pfam domain (by shortname)
MATCH (g:Gene)-[:Gene_has_pfam]->(p:Pfam {short_name: 'DNA_pol3_beta'})
RETURN g.locus_tag, g.organism_strain

-- Find all genes with a specific Pfam domain (by description/name)
MATCH (g:Gene)-[:Gene_has_pfam]->(p:Pfam)
WHERE p.name CONTAINS 'polymerase'
RETURN g.locus_tag, p.name, g.organism_strain

-- Find genes sharing domains with a gene of interest
MATCH (g1:Gene {locus_tag: 'PMM0001'})-[:Gene_has_pfam]->(p:Pfam)<-[:Gene_has_pfam]-(g2:Gene)
WHERE g1 <> g2
RETURN p.name, p.short_name, g2.locus_tag, g2.organism_strain

-- Find all domains in a clan
MATCH (d:Pfam)-[:Pfam_in_pfam_clan]->(c:PfamClan {name: 'DNA_clamp'})
RETURN d.short_name, d.name

-- Domain architecture comparison across strains
MATCH (g:Gene)-[:Gene_has_pfam]->(p:Pfam)
WHERE g.organism_strain = 'Prochlorococcus MED4'
WITH p, count(g) AS gene_count
RETURN p.short_name, p.name, gene_count
ORDER BY gene_count DESC

-- Full-text search for domains by name or shortname
CALL db.index.fulltext.queryNodes('pfamFullText', 'polymerase')
YIELD node, score
RETURN node.short_name, node.name, score

-- Find genes via domain search (2-hop)
CALL db.index.fulltext.queryNodes('pfamFullText', 'transporter')
YIELD node AS p
MATCH (g:Gene)-[:Gene_has_pfam]->(p)
RETURN p.short_name, p.name, g.locus_tag, g.organism_strain

-- List all clans and their domain count
MATCH (d:Pfam)-[:Pfam_in_pfam_clan]->(c:PfamClan)
RETURN c.name, count(d) AS domain_count
ORDER BY domain_count DESC

-- Find genes in a superfamily (2-hop via clan)
MATCH (c:PfamClan {name: 'DNA_clamp'})<-[:Pfam_in_pfam_clan]-(d:Pfam)<-[:Gene_has_pfam]-(g:Gene)
RETURN d.short_name, g.locus_tag, g.organism_strain
```

## Verification

After full rebuild:
```cypher
-- Pfam domain nodes exist
MATCH (p:Pfam) RETURN count(p)
-- Expected: ~2000

-- PfamClan nodes exist
MATCH (c:PfamClan) RETURN count(c)
-- Expected: ~800

-- All Pfam nodes have required properties
MATCH (p:Pfam) WHERE p.name IS NULL OR p.short_name IS NULL
RETURN count(p)
-- Expected: 0

-- All PfamClan nodes have name
MATCH (c:PfamClan) WHERE c.name IS NULL RETURN count(c)
-- Expected: 0

-- gene→domain edges
MATCH ()-[r:Gene_has_pfam]->() RETURN count(r)
-- Expected: ~25000

-- domain→clan edges
MATCH ()-[r:Pfam_in_pfam_clan]->() RETURN count(r)
-- Expected: ~1500

-- gene edges only point to Pfam, not PfamClan
MATCH (g:Gene)-[:Gene_has_pfam]->(c:PfamClan) RETURN count(c)
-- Expected: 0

-- pfam_ids, pfam_names, pfam_descriptions dropped from Gene nodes
MATCH (g:Gene) WHERE g.pfam_ids IS NOT NULL RETURN count(g)
-- Expected: 0

MATCH (g:Gene) WHERE g.pfam_names IS NOT NULL RETURN count(g)
-- Expected: 0

-- RSP50 should now have pfam coverage
MATCH (g:Gene {organism_strain: 'Prochlorococcus RSP50'})-[:Gene_has_pfam]->()
RETURN count(g)
-- Expected: >0 (was 0 before)

-- Per-strain coverage
MATCH (g:Gene)-[:Gene_has_pfam]->(p:Pfam)
RETURN g.organism_strain, count(DISTINCT g) AS genes_with_domains
ORDER BY genes_with_domains DESC

-- Fulltext indexes work
CALL db.index.fulltext.queryNodes('pfamFullText', 'polymerase')
YIELD node, score RETURN node.name, score LIMIT 5

CALL db.index.fulltext.queryNodes('pfamClanFullText', 'clamp')
YIELD node, score RETURN node.name, score LIMIT 5
```

## Architecture rationale: prepare_data vs adapter split

The split follows the same pattern as GO/EC/KEGG:

| Layer | Responsibility | Why here |
|-------|---------------|----------|
| **prepare_data (step 2)** | Resolve eggNOG shortnames → PF* IDs, deduplicate, write clean `pfam_ids` | Cross-source join needs reference data; result cached in `gene_annotations_merged.json` for all downstream consumers |
| **Adapter** | Create Pfam/PfamClan nodes, gene→domain edges, domain→clan edges | Graph materialization is the adapter's job |

This matches GO exactly: `gene_annotations_merged.json` supplies clean `go_terms` (IDs), then `MultiGoAnnotationAdapter` yields nodes and edges. The enrichment is a post-merge secondary derivation (same pattern as `extract_ortholog_groups`), not a config-driven transform — it requires a cross-source join via reference data.

## Decision log

- **Two node types** (not single label with `entry_type`): Pfam domains (`pfam:PF*`) and clans (`pfam.clan:CL*`) are distinct bioregistry namespaces. Using correct CURIEs is cleaner than a single label with pragmatic but semantically wrong IDs.
- **`Pfam_in_pfam_clan`** (not `Pfam_is_a_pfam`): cross-type edge; `is_a` reserved for same-type hierarchies.
- **Single raw union + post-merge enrichment** (not per-source transforms): config collects all raw tokens into one list; enrichment resolves shortnames and filters via reference data. Eliminates `extract_pfam_ids`/`extract_pfam_names` transforms and the `pfam_names`/`pfam_descriptions` intermediate fields. Follows `extract_ortholog_groups` pattern.
- **Drop pfam Gene properties now** (not staged): fully replaced by graph edges in this task; no transition period needed
- **Lazy download in step 2**: no new prepare_data step; `load_pfam_data()` downloads on first call, same as EC normalization
- **Full-text index on Pfam + PfamClan nodes**: enables domain/clan name search without going through Gene nodes
- **Shortname uniqueness is guaranteed**: Pfam-A shortnames are 1:1 with accessions; reverse lookup is safe
- **Unresolved shortnames are logged, not fatal**: a few percent loss from version skew is acceptable vs the current 34% gap
- **File location `multiomics_kg/utils/pfam_utils.py`**: matches `go_utils.py` pattern; used by both build_gene_annotations (enrichment) and adapter (node emission)
- **Protein node `pfam_ids` intentionally kept**: `uniprot_adapter.py` emits `pfam_ids` on Protein nodes (schema_config.yaml line 230). This is a separate property on a different node type — not affected by the Gene-side cleanup. Can be revisited later if Protein→Pfam edges are added.
