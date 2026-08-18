# MEROPS Pfam Bridge + Cleavage Specificity — Design

**Date:** 2026-08-18
**Status:** approved design (user-reviewed in chat), spec for review
**Parent:** `docs/superpowers/specs/2026-08-17-merops-kg-integration-design.md`
  (Phase-2 MEROPS integration — this picks up two of its recorded deferrals)
**Backlog items closed:** `plans/backlog.md` → "MEROPS cross-ontology bridges" +
  "MEROPS cleavage specificity as node properties"

## Intent

Give the MEROPS layer (a) an independent corroboration signal for its
single-source diamond calls, via MEROPS's own curated family→Pfam map, and
(b) functional depth on the family nodes, via curated cleavage-site
specificity — turning "this gene is an S01 peptidase" into "likely cleaves
after Lys/Arg".

## Measured inputs (2026-08-18, MEROPS release 12.5)

| File | Size | Content | Coverage of our 97 observed families |
|---|---|---|---|
| `interpro.txt` | 182 MB | one row per curated UniProt member of a MEROPS subfamily/identifier, with its Pfam domains (`PF00026/PF14543`) + InterPro entry | **95/97**, collapsing to **183 distinct (family, Pfam) pairs**; 62 families have exactly 1 Pfam, 18 have 2; tail (S09: 18, C26: 5) is shared-fold biology (α/β-hydrolase, GATase), not noise |
| `Substrate_search.txt` | 29 MB | ~100K curated cleavage records (`CLE*` ids): enzyme identifier, substrate, position, P4→P4' residues, physiological/non-physiological/synthetic flag | **88/97**; support ranges 39,567 events (S01) to a handful; profiles are textbook (S01 P1: Lys 36% · Arg 34%; M41: Leu + hydrophobics) |
| `GO_annotation.txt` | 1 MB | per-identifier GO terms from member proteins (all kingdoms) | 81/97 — but median 19 / max 389 terms per family incl. eukaryote-only terms |

Decision-driving live-graph measurements: only **338** of 4,197 merops-called
genes lack GO entirely and only **1** is weakly annotated (`annotation_state <=
catch_all_only`); only **4** peptidase-called genes lack Pfam. So neither
bridge has a completeness win — the Pfam bridge's value is **corroboration**,
and the GO bridge has no compensating value for its noise.

## In scope

### 1. `Merops_family_has_pfam_domain` (MeropsFamily → Pfam, ~183 edges)

- **Verb: `has`** (user-confirmed). Same claim shape as
  `Tcdb_family_has_pfam_domain`: *"MEROPS's curated member proteins of this
  family carry this domain"* — a factual composition statement. `related_to`
  stays reserved for weak recall-biased routers (InterPro Layer A);
  `associated` is used nowhere as an edge verb.
- **Family-level attachment only** (level-1 nodes; 95 of 97 observed
  families). Subfamily attachment adds nothing at median 1 Pfam/family.
- Edge property: `member_id_count` (int) — how many distinct curated MEROPS
  identifiers back the (family, Pfam) pair. **Deliberate deviation** from
  TCDB's `curated_tcids` list: MEROPS backing lists run to hundreds of
  identifiers (every characterized S01 member), so the count carries the
  support signal without a bloated array.
- **Dangling-proof** (TCDB-bridge precedent): emit only pairs whose Pfam is in
  the injected `pfam_node_ids` (`MultiPfamAnnotationAdapter.all_pfam_ids()`)
  AND whose family is a kept MeropsFamily node. `pfam_node_ids=None` → no
  bridge edges.
- **Read-direction semantics (document, don't encode):** forward
  (gene's family known → do its domains agree?) is the sound direction and the
  only one `pfam_support` uses. Backward (has PF00117 → is a C26
  peptidase) is unsafe in exactly the shared-fold families where `call_class`
  already flags dead homologs — the same biology, two coherent guards.

### 2. `pfam_support: corroborated | uncorroborated` on `Gene_has_merops_family` (post-import)

- Definition (sound direction only): `corroborated` ⇔ the gene has a
  `Gene_has_pfam` edge to a Pfam that this edge's family — the target node, or
  its level-1 ancestor when the edge attaches at a subfamily — bridges to via
  `Merops_family_has_pfam_domain`.
- **Name + shape follow the vocabulary contract's R5 conversion** (spec
  2026-08-16): the TCDB analog is `Gene_has_tcdb_family.pfam_support`
  (`corroborated|uncorroborated` string pair) — the old `pfam_corroborated`
  bool no longer exists; native bool is not admissible even for
  post-import-set properties. Set by post-import Cypher, identical wording to
  the TCDB entry.
- **Advisory, never a filter.** No composite `evidence_score`: with two
  signals (`tier`, `pfam_support`) a synthesized score adds nothing over
  reading both. The no-`annotation_quality`-bucket decision stands —
  corroboration is not a second independent evidence *source*.

### 3. Cleavage specificity — three sparse MeropsFamily properties

On family-level (level-1) nodes; **sparse** — absent when MEROPS has no
cleavage data for the family (88/97 covered), never empty:

| Property | Type | Content |
|---|---|---|
| `cleavage_p1_residues` | str[] | top P1 residues (residue N-terminal of the cut — Schechter–Berger nomenclature, the dominant specificity determinant): share ≥ 10%, max 3, order = frequency |
| `cleavage_summary` | str | readable, e.g. `"cleaves after Lys (36%) / Arg (34%) - 39567 known cleavages (25% physiological)"` — no `\|` or `'` (clean-string constraint applies to computed strings too) |
| `known_cleavage_count` | int | number of curated cleavage records backing the profile ("known cleavages" is MEROPS's own phrasing; the `CLE*` record ids literally mean cleavage) |

- Aggregated over **all** record types — synthetic assays are legitimate
  specificity evidence; the physiological share is *reported* in the summary,
  not used as a filter.
- **Honesty caveats carried by the design:** (a) support varies 4 orders of
  magnitude — hence the count; (b) family-level aggregation blends members
  (S01 blends trypsin's Lys/Arg with chymotrypsin's aromatics) — the property
  reads "what this family is known to cut", never "what your gene cuts";
  (c) events are not independent evidence (one degradomics paper contributes
  many rows).
- `cleavage_summary` joins `meropsFamilyFullText`.

## Explicitly NOT in scope (measured rejections)

- **GO bridge (`GO_annotation.txt`)** — skipped (user decision, 2026-08-18).
  All-kingdom member rollup gives median 19 / max 389 terms per family
  including biologically-wrong eukaryote terms; the completeness win is 338
  genes (~8%) vs. the 1,311 that justified TCDB's GO bridge. Revisit only
  with a concrete use case + a filtering design.
- **Merops→InterPro edge** — derivable in 2 hops via the existing
  `Pfam_in_interpro_entry` bridge; materializing it would duplicate.
- **`annotation_quality` bucket for merops** — unchanged; `pfam_support`
  is corroboration, not a second source.
- **MCP/explorer surfacing** — separate follow-up, unchanged in
  `plans/backlog.md`.

## Implementation shape

1. **Parsers** (`multiomics_kg/utils/merops_diamond.py`, pure, unit-tested):
   - `parse_interpro_txt_stream(lines) -> {family: {pfam: member_id_count}}`
     (stream — never hold 182 MB; split `PF..../PF....` tokens; family from
     the identifier column, `_MEROPS_ID_RE` reuse)
   - `aggregate_cleavages(lines) -> {family: {p1_counter, physiological_n, total_n}}`
     + `cleavage_properties(agg) -> {cleavage_p1_residues, cleavage_summary,
     known_cleavage_count}` (latin-1; threshold + top-3 + formatting in one
     tested place)
2. **Builder** (`build_merops_reference.py`, **moved to step 10** — see §
   "prepare_data changes"): download `interpro.txt` (182 MB) +
   `Substrate_search.txt` (29 MB) into `cache/data/merops/raw/` (already
   gitignored) when missing/`--refetch-raw`; add two blocks to the committed
   `merops_reference.json`: `"pfam_bridge": {family: {pfam: count}}` and
   `"cleavage": {family: {...the three properties...}}` (a few KB). QC
   warnings on empty parses (column-drift guards, existing pattern).
3. **Adapter** (`merops_adapter.py`):
   - `MultiMeropsAnnotationAdapter(pfam_node_ids=...)` new injected arg
     (TCDB pattern); `create_knowledge_graph.py` passes
     `pfam_adapter.all_pfam_ids()`.
   - `get_nodes()`: merge the three cleavage properties onto level-1 nodes
     (sparse, `_clean_str` on strings).
   - `get_edges()`: emit `merops_family_has_pfam_domain` edges
     (`{member_id_count}`) after the is-a hierarchy, pruned as above.
4. **Schema** (`schema_config.yaml`): `merops family to pfam association`
   (edge, `label_as_edge: Merops_family_has_pfam_domain`,
   `properties: {member_id_count: int}`); three new properties on
   `merops family`.
5. **Post-import** (`.sh` + `.cypher`, identical): one `CALL IN TRANSACTIONS`
   block setting `pfam_support` on every `Gene_has_merops_family` edge
   (`'uncorroborated'` default; `'corroborated'` per §2 definition). No new
   indexes (bridge is tiny); `meropsFamilyFullText` gains `cleavage_summary`
   (drop+recreate, the interproEntryFullText rerun pattern).
6. **Vocabulary contract**: add the 10 entries of §"Controlled-vocabulary
   registration" to `config/controlled_vocabularies.yaml`; import literals
   from `controlled_vocab.py` in the adapter/post-import where the tcdb
   pattern does; the existing contract kg-tests then verify the live graph
   against them.
7. **Tests**: unit (both parsers, property formatting incl. the 20-residue
   filter, adapter emission + pruning, injection-None → no edges); kg-validity
   additions to `test_merops.py` (bridge count ≈183 ± few, no dangling ends,
   S14→PF00574 spot check, corroborated-⇒-intersection-exists consistency,
   S01/S08 cleavage spot checks, sparse-property discipline, summary contains
   "known cleavages").
8. **prepare_data**: step-10 move + `--rebuild` flag per §"prepare_data
   changes".
9. **Docs**: extend `docs/kg-changes/merops-extension.md` (new section) +
   CLAUDE.md (MeropsFamily bullet, labels list, step 9/10, `--rebuild`)
   + CHANGELOG (`[Unreleased]`) + remove the two `plans/backlog.md` bullets +
   note the GO skip there is now a measured rejection.
10. **Validation gate**: fast suite → test-mode build → commit/push → user
    rebuilds → `/omics-edge-snapshot` compare + `pytest -m kg` + snapshot
    regeneration (unchanged process).

## Controlled-vocabulary registration (user-flagged, 2026-08-18)

`config/controlled_vocabularies.yaml` is THE source of truth (67 entries;
adapters import literals via `multiomics_kg/utils/controlled_vocab.py`, tests
assert the live graph matches). The shipped merops integration registered
**nothing** — this work closes that gap and registers its own additions:

**Gap-fix (shipped 2026-08-17 properties):**

| Entry | Shape |
|---|---|
| `Gene_has_merops_family.call_class` | string, closed: `peptidase`, `inhibitor`, `nonpeptidase_homolog` |
| `Gene_has_merops_family.tier` | int, 1–3 (NOT sparse — unlike tcdb, every merops edge carries tier) |
| `Gene_has_merops_family.best_hit_kind` | string, closed: `holotype`, `putative`, `nonpeptidase_homolog` |
| `MeropsFamily.level_kind` | string, closed: `merops_clan`, `merops_family`, `merops_subfamily` |
| `MeropsFamily.family_type` | string, closed: `peptidase`, `inhibitor` |
| `MeropsFamily.catalytic_type` | string, closed (9 values, full words), sparse (absent on inhibitors) |
| `Gene.merops_classes` | string_array, closed (3 call_class values), default `[]` |

**New in this work:**

| Entry | Shape |
|---|---|
| `Gene_has_merops_family.pfam_support` | string, closed: `corroborated`, `uncorroborated` (R5 pair, tcdb wording) |
| `Merops_family_has_pfam_domain.member_id_count` | int, min 1 |
| `MeropsFamily.cleavage_p1_residues` | string_array, closed: the 20 standard amino-acid three-letter codes, sparse |

Consequence for the cleavage parser: **filter P1 to the 20 standard residues**
(the raw file contains modified-residue codes like `TyI`) so the vocabulary
can be closed. Adapter/post-import literals should be imported from the
contract where `controlled_vocab.py` supports it (tcdb pattern), not
re-hardcoded.

## prepare_data changes (user-requested, 2026-08-18)

1. **MEROPS reference moves step 9 → new step 10.** Step 9's identity is
   "central references consumed by the step-2 merge" — the reason it is
   ordered before step 2 in the default run. The MEROPS reference is consumed
   only at KG-build time (`merops_adapter`), never by the merge, so parking it
   in step 9 muddied that dependency story (and now adds a 182 MB download to
   a step people run for merge freshness). Changes: `10)` case in
   `prepare_data.sh` (own `logs/prepare_data_step10.log`); default
   `STEPS="0 9 1 2 3 4 5 6 7 8 10"`; usage/labels/valid-steps text; the
   adapter's fail-loudly message (`--steps 9` → `--steps 10`); CLAUDE.md
   (step 9 + new step 10 blocks, Data Locations); `docs/kg-changes/merops-extension.md`.
2. **New `--rebuild` flag**: reruns every *derived* step with `--force` —
   equivalent to `--steps 9 1 2 3 4 5 6 7 8 10 --force` (dependency order:
   9 before 2; step 0 downloads deliberately excluded). Mutually exclusive
   with `--steps`; composes with `--strains` / `--skip-cyanorak` /
   `--refetch-raw` as usual.

## Predicted post-rebuild figures (assert ranges in kg-validity)

- `Merops_family_has_pfam_domain`: ~183 (minus any Pfams absent from the
  graph's Pfam node set — measure at build, expect small)
- Families with cleavage properties: ≤ 88 of 97 level-1 nodes
- `pfam_support = 'corroborated'`: measure and record in the kg-changes doc —
  expected majority of `call_class='peptidase'` edges, depressed among
  `nonpeptidase_homolog` edges (shared-fold families) — that contrast is
  itself a soft validation of `call_class`.
