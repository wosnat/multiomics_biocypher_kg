# MEROPS KG Integration (Phase 2) — Design

**Date:** 2026-08-17
**Status:** approved (user-reviewed schema), implementing
**Phase-1 partner:** `/merops-diamond` (`docs/superpowers/specs/2026-08-17-merops-diamond-design.md`)
**Workflow:** `/integrate-a-tool`

## Step-0 verdicts

1. **What the tool predicts:** per-protein MEROPS peptidase/inhibitor classification
   (clan → family → subfamily → identifier) from diamond blastp vs. the 5,009-sequence
   scan library, with tcdb-parity tiers. 42 strains, 4,254 candidates: 92.1% tier-3
   (family), 331 tier-2 (subfamily or ragged family), 6 tier-1 (identifier).
   97 distinct families, ~41 clans, ~30 observed subfamilies.
2. **Artifact shape:** `<strain>.merops.calls.json` keyed by WP_ accession (join key =
   `protein_id`, eggNOG-style); `{"calls": [candidate, ...]}` sorted by confidence.
   No nulls in the join key; `subfamily`/`clan`/`catalytic_type` nullable per candidate.
3. **Track: 3A ontology, hierarchical, scored** (decision tree: shared closed vocabulary
   with hierarchy + per-call score). Matches the Phase-1 doc's Phase-2 sketch — no
   supersession needed.
4. **Functional-vs-structural: FUNCTIONAL, tier-gated** (tcdb-diamond precedent):
   `annotation_types`/`informative_annotation_types` += `'merops'` only for genes with
   a tier ≤ 2 edge; **no** `annotation_quality` bucket (single-source, tier-3-dominated —
   the Phase-1 sketch's explicit call; BRITE precedent for informative-without-bucket).

## Node: `MeropsFamily` (~170 nodes, observed-only, CAZy-pattern pruning)

Levels: clan (`level 0`, `merops_clan`) → family (`level 1`, `merops_family`) →
subfamily (`level 2`, `merops_subfamily`). **No identifier-level nodes** — only 6
tier-1 calls exist; the identifier is preserved on the edge as `best_hit_id` (for
tier-1 calls the called code equals `best_hit_id`, so nothing is lost).

Node ids: `merops.clan:<CODE>` and `merops.family:<CODE>` — both are **registered
bioregistry prefixes** (verified: `normalize_curie` accepts clan, family, subfamily,
and inhibitor-family shapes) → colon CURIEs, pfam/pfam.clan-style. Subfamilies use
`merops.family:` (a subfamily code is a family-code extension, matching MEROPS usage).

| Property | Type | Source |
|---|---|---|
| `name` | str | family/subfamily: family.txt type-example name, fallback = code; clan: code |
| `description` | str, sparse | clans only — clan.txt fold/mechanism text |
| `merops_id` | str | bare code (`SC`, `S14`, `S08A`) |
| `level` / `level_kind` | int / str | adapter-emitted |
| `family_type` | str | `peptidase` \| `inhibitor` (renamed from MEROPS "entry_type" jargon) |
| `catalytic_type` | str, nullable | **full words**: serine, cysteine, metallo, aspartic, threonine, glutamic, asparagine_lyase, mixed, unknown; null for inhibitors |
| `gene_count`, `organism_count`, `member_count`, `peptidase_gene_count` | int | post-import |

## Edge: `Gene_has_merops_family` (scored; one edge per candidate; multi-domain fan-out)

Attaches at the called depth (tier-1 candidates attach at their subfamily, else family).

| Property | Type | Notes |
|---|---|---|
| `call_class` | str | **read-first verdict**: `peptidase` \| `inhibitor` \| `nonpeptidase_homolog`. Derived: inhibitor ⇐ I-family; nonpeptidase_homolog ⇐ `best_hit_kind == nonpeptidase_homolog`; else peptidase. Threshold-free. |
| `tier` | int 1–3 | tcdb convention (claim-depth confidence band) |
| `confidence_score` | float 0–1 | `(identity/100)·(qcov/100)·agreement_weight` |
| `identity`, `qcov`, `evalue` | float | best-hit alignment stats |
| `consensus_n` | int | hits backing the call |
| `best_hit_id` | str | e.g. `S14.008` — the matched MEROPS reference |
| `best_hit_kind` | str | `holotype` \| `putative` \| `nonpeptidase_homolog` |

**Deliberately dropped from the edge** (stay in committed calls.json): `scov`, `length`,
`consensus_agreement`, `homolog_hit_fraction` (redundant once `call_class` exists; the
fraction is trimodal — 2,703 at 0, 606 ≥ 0.5 — and `best_hit_kind` agrees with it in
the vast majority of cases).

Hierarchy edge: `Merops_family_is_a_merops_family` (child → parent, `{}`). Families
with MEROPS's unassigned-clan sentinel (`clan` codes ending `-`, normalized to None)
are parentless roots (InterPro precedent).

**Why one edge type, not three:** the peptidase/homolog distinction is per-edge, not
per-family (C26 holds both active peptidases and dead homologs); `expression_status` and
two-source `Gene_has_tcdb_family` are the precedents for categorical properties over
edge-type splitting; splitting would break the uniform `Gene_has_<x>` traversal pattern
and triple every rollup.

## Gene additions

- `merops_family_count` (int, post-import; ALL edges — routing, ungated, tcdb precedent)
- `merops_classes` (str[], post-import; distinct `call_class` values — the at-a-glance
  guard so `gene_details` readers never mistake a dead homolog for a protease)
- `annotation_types` / `informative_annotation_types` += `'merops'` **iff** the gene has
  a tier ≤ 2 edge
- `contributing_sources` += `'merops_diamond'` (build-time, via the merge)

## Merge front door (prepare_data step 2)

New hand-wired source `merops_diamond`, join `protein_id`:
- `load_merops()` reads `<strain>.merops.calls.json` → `{WP_: {"merops_ids": [codes]}}`
- field `merops_ids: str[]` (`union`) lands in `gene_annotations_merged.json`
- threading through `_get_raw` → all six `_resolve_*` → `build_wide`/`build_merged`
  (new arg **optional-defaulted** — the positional-call test landmine) →
  `process_strain` row join → `_compute_contributing_sources`
- DataSource: `data_source_adapter.py` `_name_for`/`_description_for` + node-count bump
  in `tests/test_data_source_adapter.py` AND `tests/kg_validity/test_data_source.py`

The adapter reads calls.json directly for per-call edge evidence (interpro/tcdb
pattern); `merops_ids` provides the merge presence + contributing_sources for free.

## Reference cache (prepare_data step 10)

`multiomics_kg/download/build_merops_reference.py` →
**committed** `cache/data/merops/merops_reference.json`:
`{"families": {code: {name, clan, family_type}}, "clans": {code: {description, family_type}}}`.
Input: `$MEROPS_DATA_DIR/DB/{family,clan}.txt` if present (Phase-1 runner put them
there), else self-downloads to `cache/data/merops/raw/` (gitignored, ~300 KB, EBI FTP).
KG build itself needs **zero** network/`~/tools` access — reference + calls.json are
both committed. Adapter fails loudly if the reference is missing. Effectively run-once
(MEROPS maintenance mode since 12.5 / Sept 2023).

## Post-import (both `post-import.sh` and `post-import.cypher`, byte-identical)

- MeropsFamily rollups: `gene_count`/`organism_count` (cazy `*0..` subtree walk),
  `member_count` (direct children), `peptidase_gene_count` (subtree, only
  `call_class = 'peptidase'` edges)
- Gene: `merops_family_count`, `merops_classes`, tier-gated fold-in to
  `annotation_types` + `informative_annotation_types`
- Indexes: `merops_family_level_idx`, `merops_family_level_kind_idx`,
  `merops_family_merops_id_idx`; full-text `meropsFamilyFullText` on `name`, `description`

## Validation gate

`/omics-edge-snapshot` before/after (expression edges unchanged); new
`tests/kg_validity/test_merops.py` (node count sanity, no `Unknown`-style junk, tier ∈
1–3, `confidence_score` ∈ (0,1], `call_class`/`best_hit_kind` vocab, hierarchy
child-level > parent-level, no dangling edges, rollup sanity, DataSource presence);
`pytest -m kg`; `snapshot_data.json` regeneration. Docker rebuild happens in the
`multiomics_biocypher_kg` clone (two-clone workflow) — hand-off to user.

## Implementation checklist (status 2026-08-17)

- [x] Merge front door: `merops_diamond` source + `merops_ids` field wired through
      `build_gene_annotations.py` (optional-defaulted `mer` arg); all 42 strains re-merged
- [x] DataSource: adapter name/description + count bumped to 9 in both tests
- [x] Reference: `build_merops_reference.py` (step 9 sub-builder) → committed
      `merops_reference.json` (100% of observed families named via scan-lib fallback)
- [x] Utils: `classify_code` / `catalytic_type_word` / `family_type` / `call_class` /
      `edge_target_code` / `parse_clan_txt` / `type_example_names` (unit-tested)
- [x] Schema: `merops family` node + scored `gene_has_merops_family` +
      `merops_family_is_a_merops_family`
- [x] Adapter: `merops_adapter.py` (per-strain + Multi) wired into `create_knowledge_graph.py`;
      MED4 smoke: 55 nodes / 33 is-a / 59 gene edges, ClpP exact
- [x] Post-import (.sh + .cypher identical): indexes + fulltext, subtree rollups +
      `peptidase_gene_count`, `merops_family_count` + `merops_classes`, tier-gated
      annotation_types/informative fold-in, has_any_edge extension
- [x] Unit suite green (2,318 passed incl. 27 new merops tests + bucket gate)
- [x] Docs: `docs/kg-changes/merops-extension.md`, CLAUDE.md, CHANGELOG.md,
      merops-diamond SKILL.md Phase-2, add-a-tool hand-off redirect
- [x] `/omics-edge-snapshot --save pre_merops_integration` captured
- [x] Docker rebuild (multiomics_biocypher_kg clone, user-driven) → post-import
- [x] `/omics-edge-snapshot --compare pre_merops_integration` — exit 0, zero regressions,
      all expression/DM/metabolism counts unchanged
- [x] `pytest -m kg` — 1,122 passed (incl. `tests/kg_validity/test_merops.py`); live figures:
      155 nodes · 108 is-a · 4,257 gene edges (3,427 peptidase / 778 nonpeptidase_homolog /
      52 inhibitor) · 339 genes tier-gated into annotation_types · C26 gap 272 vs 41 visible
- [x] `snapshot_data.json` regenerated (462 snapshot assertions green)

## Out of scope (recorded follow-ups)

- MCP/explorer registration (`ONTOLOGY_CONFIG` + ontology enum) — explorer repo
- `interpro.txt` Pfam/InterPro bridge + `GO_annotation.txt` GO bridges (TCDB-bridge-shaped)
- Cleavage-specificity node properties (`Substrate_search.txt`)
- Any `annotation_quality` bucket for merops (revisit if a second evidence source lands)
