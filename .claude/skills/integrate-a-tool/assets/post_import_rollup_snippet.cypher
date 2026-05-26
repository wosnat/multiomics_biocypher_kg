// ── <Tool> ontology — post-import snippet (TEMPLATE) ─────────────────────────
// Paste the SAME Cypher into BOTH scripts/post-import.sh and scripts/post-import.cypher
// (they must stay byte-identical). Slot by group: indexes → Group 1; node/edge rollups
// + annotation_types fold-in → Group 3 (the CALL { } IN TRANSACTIONS writes).
// Capture scripts/post-import-validate.sh > baseline.txt BEFORE the rebuild; AFTER, the diff
// must be additions-only. Replace every <…>.

// ═══ GROUP 1 — indexes ═══════════════════════════════════════════════════════
CREATE INDEX <node>_level_idx       IF NOT EXISTS FOR (n:<NodeLabel>) ON (n.level);
CREATE INDEX <node>_id_idx          IF NOT EXISTS FOR (n:<NodeLabel>) ON (n.<tool>_id);
// CREATE INDEX <node>_level_kind_idx IF NOT EXISTS FOR (n:<NodeLabel>) ON (n.level_kind);  // HIERARCHICAL only
CREATE FULLTEXT INDEX <node>FullText IF NOT EXISTS FOR (n:<NodeLabel>) ON EACH [n.name, n.<tool>_id];

// ═══ GROUP 3 — node gene_count / organism_count ══════════════════════════════
// FLAT (no hierarchy): direct gene→node traversal.
MATCH (n:<NodeLabel>)
CALL {
  WITH n
  OPTIONAL MATCH (n)<-[:Gene_has_<x>]-(g:Gene)
  WITH n, count(DISTINCT g) AS gc, collect(DISTINCT g.organism_name) AS orgs
  SET n.gene_count = gc,
      n.organism_count = size([x IN orgs WHERE x IS NOT NULL])
} IN TRANSACTIONS OF 1000 ROWS;

// HIERARCHICAL variant (ancestor gene_count includes descendants — the cazy `*0..` walk):
// MATCH (n:<NodeLabel>)
// CALL {
//   WITH n
//   OPTIONAL MATCH (n)<-[:<X>_is_a_<x>*0..]-(:<NodeLabel>)<-[:Gene_has_<x>]-(g:Gene)
//   WITH n, count(DISTINCT g) AS gc, collect(DISTINCT g.organism_name) AS orgs
//   SET n.gene_count = gc, n.organism_count = size([x IN orgs WHERE x IS NOT NULL])
//       // , n.member_count = ...        // count descendant nodes if hierarchical
// } IN TRANSACTIONS OF 1000 ROWS;

// ═══ GROUP 3 — cross-node Gene rollups (CONDITIONAL — see references/post-import-patterns.md) ═══
// FUNCTIONAL ontology ONLY: fold the new edge into the annotation_types arrays so
// gene_overview / genes_by_function reflect it. Add ONE line to EACH of the two existing
// statements (post-import.cypher ~L609 annotation_types, ~L634 informative_annotation_types):
//     CASE WHEN EXISTS { (g)-[:Gene_has_<x>]->() } THEN ['<source>'] ELSE [] END +
// If it is an INFORMATIVE functional source, also add it to the annotation_state/annotation_quality
// 8-bucket count (~L548-576) — follow the bucket-maintenance procedure in
// docs/superpowers/specs/2026-05-01-explorer-frictions-resolution-design.md (touches .cypher + .sh +
// CLAUDE.md + the bucket-count test; shifts genes_by_function min_quality).
//
// STRUCTURAL ontology (psortb localization, signalp signal-peptide): DO NOT fold into annotation_types.
// At most add an optional routing count (surfaces via gene_details g{.*}; gene_overview's curated list
// is MCP-side / out of scope):
MATCH (g:Gene)
CALL {
  WITH g
  OPTIONAL MATCH (g)-[r:Gene_has_<x>]->()
  WITH g, count(DISTINCT r) AS c
  SET g.<tool>_count = c
} IN TRANSACTIONS OF 1000 ROWS;

// ═══ GROUP 3 — per-gene rank on a SCORED edge (OPTIONAL) ══════════════════════
// Adapted from the expression rank_by_effect pattern (post-import.cypher ~L679). Ranks genes
// WITHIN each ontology node by descending score (1 = strongest evidence for that call).
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
