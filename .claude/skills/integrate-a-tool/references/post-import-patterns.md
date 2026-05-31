# Post-import patterns — rollups, indexes & cross-node fold-in (Step 4)

Edit **both** `scripts/post-import.sh` and `scripts/post-import.cypher` and keep the
Cypher byte-identical (the `.sh` groups statements into three `cypher-shell`
invocations; the `.cypher` is the reference copy). Capture
`scripts/post-import-validate.sh > baseline.txt` **before** the rebuild; **after**,
the diff must be **additions-only** (this is an addition, not a pure refactor — so
unlike a refactor it is *not* byte-identical, but nothing existing should change).

Slot by group: indexes → Group 1; node/edge rollups + annotation_types fold-in →
Group 3 (the `CALL { … } IN TRANSACTIONS` writes).

Line anchors verified 2026-05-26 — they drift; re-`grep`. See
`assets/post_import_rollup_snippet.cypher` for the copy-paste template.

## 1. Node `gene_count` / `organism_count` (+ `member_count`)

**FLAT** (psortb, signalp — one level, no hierarchy): direct gene→node traversal.

```cypher
MATCH (n:<NodeLabel>)
CALL {
  WITH n
  OPTIONAL MATCH (n)<-[:Gene_has_<x>]-(g:Gene)
  WITH n, count(DISTINCT g) AS gc, collect(DISTINCT g.organism_name) AS orgs
  SET n.gene_count = gc,
      n.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 1000 ROWS;
```

**HIERARCHICAL** (an ancestor's gene_count includes descendants): the cazy `*0..`
subtree walk — verbatim from the live `CazyFamily` block
([post-import.cypher:~889](../../../../scripts/post-import.cypher#L889)):

```cypher
MATCH (c:CazyFamily)
CALL {
  WITH c
  OPTIONAL MATCH (c)<-[:Cazy_family_is_a_cazy_family*0..]-(desc:CazyFamily)<-[:Gene_has_cazy_family]-(g:Gene)
  WITH c, count(DISTINCT g) AS gc, collect(DISTINCT g.organism_name) AS orgs
  SET c.gene_count = gc,
      c.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 1000 ROWS;
```

`member_count` (hierarchical only) = direct-child count, computed alongside — see the
`TcdbFamily` block ([post-import.cypher:~861](../../../../scripts/post-import.cypher#L861)),
which does `member_count` + `gene_count` + `organism_count` + `metabolite_count` in one pass.

## 2. Indexes (Group 1)

Model the cazy/tcdb index cluster ([post-import.cypher:~61–73](../../../../scripts/post-import.cypher#L61)):

```cypher
CREATE INDEX <node>_level_idx IF NOT EXISTS FOR (n:<NodeLabel>) ON (n.level);
CREATE INDEX <node>_id_idx    IF NOT EXISTS FOR (n:<NodeLabel>) ON (n.<tool>_id);
-- CREATE INDEX <node>_level_kind_idx IF NOT EXISTS FOR (n:<NodeLabel>) ON (n.level_kind);  -- HIERARCHICAL only
CREATE FULLTEXT INDEX <node>FullText IF NOT EXISTS FOR (n:<NodeLabel>) ON EACH [n.name, n.<tool>_id];
```

## 3. Cross-node rollups — CONDITIONAL on functional-vs-structural (decision-tree.md, Decision 2)

**FUNCTIONAL ontology only** — fold the new edge into the two annotation arrays so
`gene_overview` / `genes_by_function` reflect it. Add **one line** to each of the two
existing `SET` statements:
- `g.annotation_types` ([post-import.cypher:~609](../../../../scripts/post-import.cypher#L609)) — already ends with the `tcdb`/`cazy` lines; add:
  ```cypher
      CASE WHEN EXISTS { (g)-[:Gene_has_<x>]->() } THEN ['<source>'] ELSE [] END +
  ```
- `g.informative_annotation_types` ([post-import.cypher:~634](../../../../scripts/post-import.cypher#L634)) — same pattern.

If it is an **informative** functional source, also add it to the
`annotation_state`/`annotation_quality` **8-bucket count**
([post-import.cypher:~544–591](../../../../scripts/post-import.cypher#L544)). The bucket
list is *enumerated*, not auto-discovered — follow the maintenance procedure in
`docs/superpowers/specs/2026-05-01-explorer-frictions-resolution-design.md` (touches
`.cypher` **+** `.sh` **+** CLAUDE.md **+** the bucket-count test) and note that it
**shifts `genes_by_function` `min_quality` results**.

**STRUCTURAL ontology** (psortb localization, signalp signal-peptide): **do NOT** fold
into `annotation_types`/`annotation_quality` — it would inflate functional-annotation
coverage and skew `min_quality`. At most add a dedicated `<tool>_count` Gene signal
(surfaces via `gene_details` `g{.*}` for free; adding it to `gene_overview`'s curated
list is an MCP-side change, out of scope):

```cypher
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[r:Gene_has_<x>]->()
  WITH g, count(DISTINCT r) AS c
  SET g.<tool>_count = c
} IN TRANSACTIONS OF 1000 ROWS;
```

`OrganismTaxon` and `Publication` need **no per-ontology rollup**.

## 4. Per-gene rank on a SCORED edge (OPTIONAL)

Adapted from the expression `rank_by_effect` UNWIND pattern
([post-import.cypher:~678](../../../../scripts/post-import.cypher#L678)). Ranks genes
WITHIN each ontology node by descending score (1 = strongest evidence):

```cypher
MATCH (n:<NodeLabel>)
CALL {
  WITH n
  MATCH (g:Gene)-[r:Gene_has_<x>]->(n)
  WITH r, r.score AS s, g.locus_tag AS lt
  ORDER BY s DESC, lt ASC
  WITH collect(r) AS edges
  UNWIND range(0, size(edges) - 1) AS i
  SET (edges[i]).rank_by_score = i + 1
} IN TRANSACTIONS OF 1000 ROWS;
```
