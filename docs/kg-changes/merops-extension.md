# MEROPS peptidase ontology (MeropsFamily nodes)

**Date:** 2026-08-17
**Spec:** [`docs/superpowers/specs/2026-08-17-merops-kg-integration-design.md`](../superpowers/specs/2026-08-17-merops-kg-integration-design.md)
**Phase 1:** [`/merops-diamond`](../../.claude/skills/merops-diamond/SKILL.md) ([design](../superpowers/specs/2026-08-17-merops-diamond-design.md))
**Track:** 3A ontology · hierarchical (clan → family → subfamily) · scored · functional (tier-gated)

## What's changing

The KG gains the authoritative protease classification. `/merops-diamond`'s
per-strain `<strain>.merops.calls.json` (diamond blastp vs. the 5,009-sequence
MEROPS scan library, release 12.5) now flows through the gene-annotation merge
as the `merops_diamond` source (field `merops_ids: str[]`, join `protein_id`),
and a new `merops_adapter` materializes an observed-only `MeropsFamily`
ontology with scored `Gene_has_merops_family` edges. Single evidence source —
eggNOG does not emit MEROPS ids and InterProScan has no MEROPS member DB.

The central design concern is **making three easy-to-miss distinctions
un-missable to an LLM reader**, at each read-point:

1. **On the edge — `call_class`** (`peptidase` | `inhibitor` |
   `nonpeptidase_homolog`): MEROPS files catalytically dead relatives inside
   peptidase families under `.9xx` identifiers; a gene whose best hits are
   those is fold evidence, NOT protease evidence. Threshold-free (derived from
   `best_hit_kind` + the I-family flag), `expression_status` pattern.
2. **On the Gene — `merops_classes: str[]`**: distinct `call_class` values, so
   `gene_details` never shows a bare count that reads as "N proteases" when
   the calls are dead homologs or inhibitors.
3. **On the family node — `peptidase_gene_count` vs `gene_count`**: the
   metabolite-count-split lesson. C26/C44-style families are dominated by dead
   relatives; the gap between the two counts shows it before anyone does ORA
   or per-organism protease counting with the wrong number.

## New node label: `MeropsFamily`

- **155 nodes** (observed-only + ancestors, CAZy pruning): 41 `merops_clan`
  (level 0) · 97 `merops_family` (level 1) · 17 `merops_subfamily` (level 2).
- IDs: `merops.clan:SC` / `merops.family:S14` / `merops.family:S08A` — both
  prefixes registered in bioregistry → colon CURIEs.
- **No identifier-level (S01.001) nodes** — only 6 tier-1 calls exist in the
  whole corpus; they attach at subfamily/family and the identifier survives as
  edge `best_hit_id` (for tier-1 calls identical to the called code).
- Properties (adapter): `name` (family/subfamily: MEROPS type-example name —
  family.txt where curated, else derived from the scan library's `.001`/lowest
  holotype, organism suffix stripped; 100% of observed families named; clans:
  code), `description` (sparse, clans only — clan.txt fold/mechanism text),
  `merops_id`, `level`, `level_kind`, `family_type`
  (`peptidase`|`inhibitor`), `catalytic_type` (**full words**: serine,
  cysteine, metallo, aspartic, threonine, glutamic, asparagine_lyase, mixed,
  unknown; null on inhibitors — one-letter MEROPS codes are insider jargon).
- Computed (post-import): `gene_count` + `organism_count` (subtree `*0..`,
  CAZy pattern), `peptidase_gene_count` (subtree, `call_class='peptidase'`
  edges only), `member_count` (direct children).
- Node metadata source: **committed `cache/data/merops/merops_reference.json`**
  (594 families / 139 subfamily names / 123 clans), built by prepare_data
  **step 9** (`build_merops_reference.py`) from `family.txt` + `clan.txt` +
  `merops_scan.lib` — reads `$MEROPS_DATA_DIR/DB/` when present (Phase-1
  install), else self-downloads (~2 MB). KG build itself needs zero network.

## New edge types

| Edge | Source → Target | Count | Properties |
|---|---|---|---|
| `Gene_has_merops_family` | Gene → MeropsFamily | ~4.2K (≤4,254 candidates; one per candidate, multi-domain proteins fan out) | `call_class`, `tier` (1–3), `confidence_score` (0–1), `identity`, `qcov`, `evalue`, `consensus_n`, `best_hit_id`, `best_hit_kind` (`holotype`\|`putative`\|`nonpeptidase_homolog`) |
| `Merops_family_is_a_merops_family` | MeropsFamily → MeropsFamily (child → parent) | 108 | none |

- Edges attach at the **called depth** (tier-truncated: tier 1 = identifier
  claim → attaches at subfamily/family; tier 2 = subfamily; tier 3 = family).
  Families MEROPS leaves clan-unassigned are hierarchy roots.
- Expected `call_class` split (42 strains): peptidase 3,425 ·
  nonpeptidase_homolog 777 · inhibitor 52. Tier shape: ~92% tier 3 (remote
  homology, conservative-by-design — not a defect), ~8% tier 2, 6 tier-1.
- Deliberately NOT on the edge (stay in committed calls.json): `scov`,
  `length`, `consensus_agreement`, `homolog_hit_fraction` (redundant once
  `call_class` exists).

## Cross-node rollups

- **Functional, TIER-GATED (tcdb precedent):** `'merops'` folds into
  `Gene.annotation_types` and `Gene.informative_annotation_types` **only for
  genes with a tier ≤ 2 edge** (~330 genes) — tier-3 calls stay findable
  through the edges without inflating annotation coverage.
- **NO `annotation_quality` bucket** (bucket count stays 9): single evidence
  source, tier-3-dominated — the Phase-1 spec's explicit call; BRITE precedent
  for informative-without-bucket. Revisit if a second MEROPS evidence source
  lands (`interpro.txt` Pfam bridge / Foldseek).
- `Gene_has_merops_family` **does** count toward `has_any_edge`
  (`catch_all_only` vs `no_evidence`) — interpro/tcdb-tier-3 precedent.
- Gene routing signals (post-import): `merops_family_count` (int, ALL edges,
  ungated), `merops_classes` (str[], sorted distinct `call_class` values).
- `contributing_sources` += `'merops_diamond'` (build-time, via the merge);
  9th `DataSource` node (`data_source:merops_diamond`, "MEROPS (diamond)").
- `OrganismTaxon` / `Publication`: no per-ontology rollup (unchanged).

## Indexes

Scalar: `merops_family_level_idx`, `merops_family_level_kind_idx`,
`merops_family_id_idx` (on `merops_id`). Full-text: `meropsFamilyFullText` on
`name`, `merops_id`, `description`.

## Validation results

- Unit: `tests/test_merops_adapter.py` (vocab helpers, node/edge shape,
  tier-1 attachment, call_class, reference-missing failure) — 27 tests; full
  fast suite 2,318 passed.
- MED4 real-data smoke: 55 nodes / 33 is-a / 59 gene edges (= 59 candidates);
  ClpP (PMM1313) → `merops.family:S14` `{call_class: peptidase, tier: 3,
  best_hit_id: S14.008}`.
- `/omics-edge-snapshot`: `pre_merops_integration` baseline captured;
  compare after the Docker rebuild (expression edges must be unchanged).
- `pytest -m kg`: `tests/kg_validity/test_merops.py` added (node counts +
  level distribution, CURIE form, S14/C26 spot checks, hierarchy sanity,
  call_class/tier/score ranges, tier-gate consistency, rollup invariants,
  MED4-vs-MIT1002 lifestyle gradient) — runs against the rebuilt graph.
- `snapshot_data.json`: regenerate after the rebuild.

## What this enables

- "Every subtilisin-family gene across Alteromonas, ranked by confidence":
  `MATCH (g:Gene)-[r:Gene_has_merops_family]->(m:MeropsFamily {merops_id:'S08'})
  WHERE r.call_class = 'peptidase' RETURN g, r ORDER BY r.confidence_score DESC`
- Honest per-organism protease counting via `peptidase_gene_count` /
  `call_class`, joined with the SignalP + PSORTb layers for the secreted
  exoprotease story (heterotroph DOM degradation).
- Protease-inhibitor biology via `family_type = 'inhibitor'`.

## What does NOT change

- Phase-1 `/merops-diamond` runner + artifacts untouched.
- No existing node/edge/property modified — additions-only, EXCEPT:
  `annotation_state`/`annotation_quality` can shift `no_evidence` →
  `catch_all_only` for genes whose only functional edge is a merops call
  (deliberate, has_any_edge extension).
- **MCP / explorer integration is out of scope** — registering `merops` in the
  explorer's `ONTOLOGY_CONFIG` + ontology enum is an explorer-repo follow-up.
  Until then, `run_cypher` reaches everything.
- Deferred (recorded): MEROPS `interpro.txt` Pfam/InterPro bridges,
  `GO_annotation.txt` GO bridges, cleavage-specificity (`Substrate_search.txt`)
  node properties.
