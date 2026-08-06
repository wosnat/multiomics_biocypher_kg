# InterProScan InterPro entries (`InterproEntry` nodes)

**Date:** 2026-07-26
**Spec:** [`docs/superpowers/specs/2026-07-26-interproscan-kg-integration-design.md`](../superpowers/specs/2026-07-26-interproscan-kg-integration-design.md)
**Track:** 3A ontology · **hierarchical** · **scored** edge · **functional** (what the gene does)
**Skill:** produced with [`/integrate-a-tool`](../../.claude/skills/integrate-a-tool/SKILL.md) (Phase-2 partner to `/add-a-tool`)
**Commits:** `b4a85cc2` (cache/merge) · `00423e7e` (adapter/schema/post-import)

> **For the MCP / explorer:** this doc is the integration contract — node/edge labels,
> full property list, id conventions, indexes, and the ORA guidance needed to surface
> InterPro without re-deriving any of it.

## What's changing

InterProScan 5 predicts protein domains/families by direct HMM/profile matching across ~15 member
DBs (Pfam, NCBIfam/TIGRFAM, PANTHER, Gene3D, PROSITE, SFLD, HAMAP, …), integrating each hit into a
unified **InterPro entry** (`IPRxxxxxx`). Phase 1 (`/interproscan-run`) ran IPS on 42 strains,
writing committed `<strain>.interproscan.calls.json`. This change wires those calls into the KG as a
new **hierarchical ontology**: distinct InterPro entries become `InterproEntry` nodes, each gene gets
one `Gene_has_interpro_entry` edge per entry (carrying domain coordinates + best e-value/score +
member DBs), the InterPro is-a hierarchy becomes `Interpro_entry_is_a_interpro_entry` edges, and a
`Pfam_in_interpro_entry` bridge links the existing eggNOG Pfam layer to the InterPro layer.

Two data inputs: (1) per-strain calls.json for the per-match edge evidence; (2) a committed
**InterPro reference** (`cache/data/interpro/interpro_reference.json`, built by **prepare_data step
9** from InterPro's `entry.list` + `ParentChildTreeFile.txt`) for node names/types + the is-a tree.

**Why it's genuinely additive over eggNOG/UniProt (not a merge):** InterProScan is a
*method-independent* cross-check on eggNOG's orthology-transferred Pfam/GO. On MED4, ~88% of Pfam
accessions and 99% of per-gene Pfam calls agree — but IPS uniquely contributes NCBIfam/TIGRFAM +
structural superfamilies + coordinates + the integrating InterPro-entry layer. See spec §0.

## New node label: `InterproEntry`

- **Count:** ~12,999 nodes (entries observed on any strain **+** their is-a ancestors — pruned like
  TCDB/BRITE so hierarchy edges never dangle).
- **ID prefix `interpro:`** (colon — `interpro` **is** a registered bioregistry prefix, unlike
  `psortb`/`signalp`), e.g. `interpro:IPR000685`.
- **Properties (adapter-emitted):**
  - `name` (str) — InterPro entry description
  - `interpro_id` (str) — `IPR*`
  - `interpro_type` (str) — one of `FAMILY | DOMAIN | HOMOLOGOUS_SUPERFAMILY | REPEAT | CONSERVED_SITE | ACTIVE_SITE | BINDING_SITE | PTM`
  - `level` (int) — is-a depth in ParentChildTree (`0` = parentless root). **Sparse hierarchy:** ~86% of entries are `level 0`.
- **Computed (post-import):** `gene_count` (DIRECT genes via `Gene_has_interpro_entry`),
  `organism_count`, `member_count` (direct child entries), `is_promiscuous` (bool — `gene_count >= 1000`).

## New edge types

| Edge | Source → Target | Count | Properties |
|---|---|---|---|
| `Gene_has_interpro_entry` | Gene → InterproEntry | ~397K | `start:int`, `end:int` (domain envelope: min start / max end across matches); `evalue:float` (best/min non-null — **omitted** when no match reports one); `score:float` (best/max non-null — omitted when N/A); `libraries:str[]` (member DBs, e.g. `PFAM\|NCBIFAM`); `match_count:int` |
| `Interpro_entry_is_a_interpro_entry` | InterproEntry → InterproEntry | ~1,569 | — (child → parent; is-a only, **100% within-type**) |
| `Pfam_in_interpro_entry` | Pfam → InterproEntry | ~5,236 | — (member-of bridge; overlap link, **not** a merge) |

**Scored ontology edge:** `Gene_has_interpro_entry` carries per-match evidence on the edge (modeled
on the `Changes_expression_of` property block, as with PSORTb/SignalP). One edge per (gene, entry);
a gene's multiple matches to one entry are folded into a single edge (envelope + best scores +
member-DB list).

### No e-value cutoff (important for query semantics)

`evalue`/`score` are **evidence/ranking payload, never a filter**. InterProScan pre-filters every
match with each member DB's own curated threshold (Pfam/NCBIfam gathering cutoffs; profile/pattern
score cutoffs), so all matches are already significant. ~42% of matches legitimately carry **no
e-value** (HAMAP/PROSITE/PANTHER/PRINTS profile & pattern methods) — an e-value filter would wrongly
drop them. Do **not** threshold on `evalue` in the MCP; use it only to rank a gene's domains.

## Cross-node rollups (functional, bounded)

- **Folded in:** `Gene.annotation_types` gains `'interpro'`; new `Gene.interpro_entry_count` (int)
  routing signal.
- **Deliberately deferred** (recorded in spec §5, **not** done): `Gene.informative_annotation_types`
  (needs an InterPro term-informativeness filter) and the `annotation_quality` 8-bucket count
  (InterPro is highly redundant with the `pfam`/`go`/`ec` buckets → would lift few genes' quality for
  a large blast radius). **So `annotation_quality` / `informative_annotation_types` are unchanged by
  this extension.**

## ORA guidance (for the explorer) — stratify by `(interpro_type, level)`

InterPro **inverts** the usual level-stratification. Elsewhere (GO/KEGG/TCDB/CAZy) `level` is the
granularity proxy; here the hierarchy is shallow/sparse (most entries `level 0`), so **`interpro_type`
is the PRIMARY ORA stratification axis** (FAMILY vs promiscuous DOMAIN vs broad
HOMOLOGOUS_SUPERFAMILY) and `level` is secondary (separates parent/child generations where the
hierarchy exists). Over-Representation Analysis over InterPro is valid (domain/family enrichment,
complementary to pathway/GO ORA) **only** when stratified by `(interpro_type, level)` — run naively
over all types at once and broad domains/superfamilies dominate. `is_promiscuous` flags ultra-common
entries so a band can down-weight them. Both `interpro_type` and `level` are indexed.

## Indexes

- Scalar: `interpro_entry_level_idx` (level), `interpro_entry_type_idx` (interpro_type — the ORA key),
  `interpro_entry_id_idx` (interpro_id).
- Full-text: `interproEntryFullText` on `name`.

## Merge / provenance

New `interproscan` logical source in `gene_annotations_config.yaml` → merged field `interpro_entries`
(list[str]) joined on RefSeq WP_ `protein_id`. Adds a 7th `DataSource` node
(`data_source:interproscan`, provenance `tool_run`) and `'interproscan'` to each annotated gene's
`contributing_sources`. The rich per-match evidence is **not** in the merged JSON — the adapter reads
the calls.json directly (like `tcdb_adapter` reads `tcdb_pruned.json`).

## Out of scope / follow-ups

- **Member-DB signature nodes** (the ~31% unintegrated hits) — dropped by design (redundant
  structural / non-family predictors; ~10-12 genes/strain lose only a functional signal, all retaining
  eggNOG/UniProt/Pfam).
- **InterPro → GO / pathway cross-links** — empty in current artifacts; needs a Phase-1 re-run with
  `--goterms --pathways`. This is the *enrichment-enabling* follow-up for pathway/GO ORA (most
  valuable for heterotrophs). This pass is deliberately enrichment-neutral.
- **Domain-composition** (`contains`/`found-in`) cross-type relationship — not in ParentChildTreeFile;
  clean future `Interpro_entry_contains_interpro_entry` edge.
- `informative_annotation_types` + `annotation_quality` bucket — spec §5.
- MCP/explorer surfacing (incl. `(type, level)` ORA) — separate task, this doc is its contract.

## Validation results (deployed graph, 2026-08-06)

- **Counts:** 12,999 `InterproEntry` · 397,342 `Gene_has_interpro_entry` · 1,569 hierarchy · 5,236
  Pfam bridge. Node types: FAMILY 6,490 · DOMAIN 4,355 · HOMOLOGOUS_SUPERFAMILY 1,533 ·
  CONSERVED_SITE 390 · ACTIVE_SITE 95 · REPEAT 74 · BINDING_SITE 55 · PTM 7.
- **Coverage:** 102,895 genes (~85%) carry an InterPro annotation (`'interpro'` in `annotation_types`).
- **Additive-only:** `/omics-edge-snapshot` before vs after — expression edges 232,758 → 232,758 (+0),
  no publication lost edges.
- **`is_promiscuous` (22 flagged, threshold ≥1000):** the flagged set is exactly the textbook-ubiquitous
  entries — IPR027417 P-loop NTPase superfamily (6,909 genes), IPR036291 NAD(P)-binding superfamily
  (3,163), IPR036388 Winged-helix DNA-binding (2,942), IPR003593 AAA+ ATPase (2,892), IPR013785
  TIM-barrel (1,846) — mostly HOMOLOGOUS_SUPERFAMILY, confirming the `(type, level)` ORA rationale.
  gene_count distribution: p50=12, p95=94, p99=307.
- **Reference vs calls.json:** 0 observed entries missing, 0 type disagreements, 100% within-type is-a.
- **No dangling edges;** all 23 kg-validity assertions pass (`test_interpro.py` + `test_data_source.py`);
  full `tests/kg_validity/` suite green.
