# KG Changes for Pfam Ontology Support

Spec for `multiomics_biocypher_kg` changes needed to support Pfam domain/clan
queries in the explorer's ontology MCP tools.

Companion to `pfam_normalization.md` (implementation plan). This document
describes the KG output that the explorer should expect after the Pfam changes
are deployed.

---

## 1. New node types

Two separate node types with correct bioregistry CURIEs:

### Pfam (domain families)

| Property | Type | Example |
|----------|------|---------|
| `id` | string | `pfam:PF00712` |
| `name` | string | `"DNA polymerase III beta subunit, N-terminal domain"` |
| `short_name` | string | `"DNA_pol3_beta"` |

Expected count: ~2,000 (only domains referenced by genes in the KG).

### PfamClan (superfamilies)

| Property | Type | Example |
|----------|------|---------|
| `id` | string | `pfam.clan:CL0060` |
| `name` | string | `"DNA_clamp"` |

Expected count: ~800 (only clans referenced by emitted domains).

### Node ID formats

- Domain IDs: `pfam:PF00712` (bioregistry prefix `pfam`, pattern `^PF\d{5}$`)
- Clan IDs: `pfam.clan:CL0060` (bioregistry prefix `pfam.clan`, pattern `^CL\d+$`)

These are distinct bioregistry namespaces — the ID prefix distinguishes the
node type unambiguously.

## 2. New edge types

| Edge label | Source | Target | Count | Description |
|-----------|--------|--------|-------|-------------|
| `Gene_has_pfam` | `Gene` | `Pfam` | ~25,000 | Gene has this protein domain |
| `Pfam_in_pfam_clan` | `Pfam` | `PfamClan` | ~1,500 | Domain belongs to superfamily clan |

The hierarchy is flat (2 levels only: domain → clan). ~45% of domains have
no clan membership and therefore no `Pfam_in_pfam_clan` edge.

## 3. Gene node property changes

The following properties are **dropped from Gene nodes** in this change:
- `pfam_ids` — replaced by `Gene_has_pfam` edges
- `pfam_names` — replaced by `Pfam.short_name` property
- `pfam_descriptions` — was broken (mixed TIGR/InterPro/Pfam), replaced by `Pfam.name`

The `geneFullText` index is updated to remove `pfam_names`. Pfam search
moves to the new `pfamFullText` index on Pfam nodes.

## 4. Indexes

```cypher
-- Scalar indexes
CREATE INDEX pfam_name_idx IF NOT EXISTS FOR (p:Pfam) ON (p.name);
CREATE INDEX pfam_clan_name_idx IF NOT EXISTS FOR (c:PfamClan) ON (c.name);

-- Full-text indexes
CREATE FULLTEXT INDEX pfamFullText IF NOT EXISTS
  FOR (p:Pfam) ON EACH [p.name, p.short_name];
CREATE FULLTEXT INDEX pfamClanFullText IF NOT EXISTS
  FOR (c:PfamClan) ON EACH [c.name];
```

## 5. Key differences from other ontologies

| Aspect | Pfam | GO/EC/KEGG |
|--------|------|------------|
| Node types | 2 (`Pfam`, `PfamClan`) | 1-3 per ontology |
| Hierarchy depth | 2 (domain → clan) | Deep DAG |
| Not all terms have parents | Yes (~45% of domains have no clan) | Rare |
| Extra text property | `short_name` (indexed) | None |
| Hierarchy edge name | `Pfam_in_pfam_clan` (cross-type) | `X_is_a_X` (same-type) |

## Verification

After KG rebuild:

```cypher
-- Pfam domain nodes
MATCH (p:Pfam) RETURN count(p)
-- Expected: ~2000

-- PfamClan nodes
MATCH (c:PfamClan) RETURN count(c)
-- Expected: ~800

-- All Pfam nodes have required properties
MATCH (p:Pfam) WHERE p.name IS NULL OR p.short_name IS NULL
RETURN count(p)
-- Expected: 0

-- Gene→domain edges
MATCH ()-[r:Gene_has_pfam]->() RETURN count(r)
-- Expected: ~25000

-- Domain→clan hierarchy edges
MATCH ()-[r:Pfam_in_pfam_clan]->() RETURN count(r)
-- Expected: ~1500

-- Gene edges only point to Pfam, not PfamClan
MATCH (g:Gene)-[:Gene_has_pfam]->(c:PfamClan) RETURN count(c)
-- Expected: 0

-- pfam_ids dropped from Gene nodes
MATCH (g:Gene) WHERE g.pfam_ids IS NOT NULL RETURN count(g)
-- Expected: 0

-- Fulltext search works
CALL db.index.fulltext.queryNodes('pfamFullText', 'polymerase')
YIELD node, score RETURN node.name, score LIMIT 5

-- Hierarchy expansion (find genes in a clan)
MATCH (c:PfamClan {name: 'DNA_clamp'})
  <-[:Pfam_in_pfam_clan]-(d:Pfam)
  <-[:Gene_has_pfam]-(g:Gene)
RETURN d.short_name, g.locus_tag, g.organism_strain
```
