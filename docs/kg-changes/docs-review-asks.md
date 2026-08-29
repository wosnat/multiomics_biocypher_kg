# Explorer docs-review asks (DOC-001…008) — what changed in the KG

**Date:** 2026-08-29 · **Driver:** `multiomics_explorer/docs/kg-specs/2026-08-29-docs-review-kg-asks.md`
(the KG review is appended there as §4). **Scope:** four graph changes, one vocabulary hygiene pass,
three doc confirmations. No adapter emits a new node or edge type; no `prepare_data` rerun is needed.

## DOC-001 — eggNOG transfer is `family_inferred` everywhere (graph, BREAKING)

**Root cause.** `annotation_provenance.annotation_edge_props` (GO×3 / EC / Pfam / CAZy) read the rung
from the sparse `<field>_evidence` map and defaulted a missing entry to `curated`; eggNOG-only tokens
never get an entry (InterPro never touches them), so all ~750K of them read `curated`. Two copies of
`_CURATED_SOURCES = {ncbi, cyanorak, uniprot, eggnog}` (provenance module + step-2 fold) made any
eggNOG-backed token InterPro *did* touch `curated` as well. The Pfam-vs-CAZy asymmetry on
`['eggnog','interproscan']` is a third mechanism: at fold time `pfam_ids_source` is still keyed by
eggNOG *shortnames* while InterPro contributes `PF*` accessions, so the eggNOG source is invisible and
Pfam lands on `signature`; `enrich_pfam_fields` later re-keys the source map (union becomes
`['eggnog','interproscan']`) but the rung is not revisited.

**Fix.** `eggnog` removed from both `_CURATED_SOURCES`. New `derive_evidence(sources, recorded)` resolves
the rung from `sources` at KG-build time, never trusting the merge value verbatim: curated source ⇒
`curated`; recorded `signature`/`homology` kept; `eggnog ∈ sources` ⇒ `family_inferred` (beats
`domain_inferred`); interproscan-only keeps its recorded strength; no provenance at all ⇒ `curated`
(legacy rows). `_fold_interpro_field` writes the same answer into `gene_annotations_merged.json` on
the next step-2 run, but the adapter-side derivation means the fix lands on the next `docker compose`
build without one. `Gene_has_pfam.evidence` vocabulary gains `family_inferred`.

**Expected movement** (from the 2026-08-29 live cross-tab): eggNOG-only edges — GO-BP 434,094,
GO-MF 180,473, GO-CC 100,481, EC 11,676, Pfam 24,700, CAZy 744 — go `curated → family_inferred`;
`['eggnog','interproscan']` pairs on EC (2,464), CAZy (801) and GO go the same way; Pfam's 59,338 stay
`signature`. `evidence_score` on the moved edges loses the "curated or signature" signal (2/3 → 1/3 for
single-source). Live gate: `tests/kg_validity/test_annotation_trust.py::test_eggnog_only_edges_are_family_inferred`
+ `test_pfam_and_cazy_agree_on_the_eggnog_interproscan_pair`.

## DOC-002 — KEGG global / overview maps flagged (graph + vocab)

`config/uninformative_terms.yaml` `kegg_term` gains an `ids:` list of **11 of the 13** parentless pathway
nodes (`ko01100 … ko01250`, enumerated live), mirrored in the F1.1 block of `post-import.sh` / `.cypher`.
`ko01310` Nitrogen cycle (16 KOs — a strict subset of `ko00910` Nitrogen metabolism, 36) and `ko01320`
Sulfur cycle (22 KOs) are parentless overview maps too, but they are narrow, class-bearing subsets rather
than unions (`ko01100` has 1,635 KOs; the median pathway has 4), so they stay **informative** — a nitrogen
researcher should never see "Nitrogen cycle" called uninformative.
Category / subcategory nodes (6 + 46) stay **unflagged** on purpose: "Carbohydrate metabolism" carries a
class signal, which is the file's guiding principle; consumers gate those with `level`. Live gate:
`test_uninformative_terms.py::test_kegg_global_maps_are_flagged` (yaml ids ∪ the two kept maps == live parentless set, flag ⇔ listed) and
`test_kegg_category_levels_stay_unflagged`. `annotation_state` does not move — the `kegg` bucket reads
`Gene_has_kegg_ko` edges to KO nodes only.

## DOC-003 — vocabulary descriptions are researcher-facing (vocab)

17 entries rewritten (`Experiment.compartment`, `MetaboliteAssay.value_kind`, `DerivedMetric.metric_type`,
`BriteCategory.tree` / `.level_kind`, `Metabolite.evidence_sources`, `Organism_has_metabolite.evidence_sources`,
`Changes_expression_of.expression_status`, both `metric_bucket`s, `Gene.annotation_state` / `.gene_category`,
`TcdbFamily.level_kind`, `MeropsFamily.level_kind` / `.catalytic_type`, `Gene_has_merops_family.call_class` /
`.best_hit_kind`) plus the 9 `is_uninformative` presence markers (script path dropped). Provenance now sits
in a `# provenance (moved out of the researcher-facing description, 2026-08-29):` comment above each key
— the `table_scope` precedent. Unit gate: `test_descriptions_carry_no_build_provenance` (regex
`CLAUDE\.md|\.py|\.yaml|\.cypher|\.sh|Phase-|addendum|harvested|controller ruling|spec §|lines ~N`).
The `description` text is **not** hashed, so `controlled_vocabularies_hash` is unchanged by this pass.
Not done: a property-level description for cross-edge props (`evidence`, `sources`) — the explorer
labels the owner instead; file a follow-up if the envelope text is wrong often enough to matter.

## DOC-004 — declared-unused values (vocab)

`spent_medium` / `lysate` pruned from `COMPARTMENTS` (`multiomics_kg/vocab/non_de_evidence.py`), the
vocab entry, the paperconfig skill doc and CLAUDE.md; the validator rejects them until re-added.
`MeropsFamily.catalytic_type` keeps `glutamic`: it mirrors MEROPS's complete code set (an external
vocabulary mirrored whole, like `NcbifamFamily.family_type`), and the description now says a value may
be absent from a given build. No `declared_only` marker was added — two values did not justify a loader
feature; prune when KG-minted, describe when external.

## DOC-006 — `KeggTerm.direct_gene_count` KO-only (graph)

`SET n.direct_gene_count = CASE WHEN n.level_kind = 'ko' THEN dgc ELSE null END` in the KEGG rollup
(both post-import files) — a null SET removes the property. Live gate:
`test_annotation_trust.py::test_kegg_direct_gene_count_only_on_kos`.

## DOC-005 / 007 / 008 — confirmations, no change

- **DOC-005:** `Gene.ec_numbers`, `cog_category`, `kegg_ko`, `kegg_pathway` (etc.) were Gene properties
  until `285d95a5` (2026-03-16, "remove lots of gene properties") — three months before the first
  tagged release (`kg-0.1.0-alpha.3`, 2026-06-06). `ko_terms` / `kegg_ids` / `cog_categories` never
  existed under those names. Nothing to log under `breaking_changes`; the explorer's
  `schema_baseline.yaml` line 743 `ec_numbers` belongs to `Protein` / `Reaction`, both of which still
  carry it. `Gene.catalytic_activities` (8,084 genes, UniProt) is the only Gene-level chemistry scalar.
- **DOC-007:** done in `e91f20ce` — `vocabulary-contract.md` R5 names the two presence markers.
- **DOC-008:** F1.1 flags exactly 9 `TigrRole` nodes: subroles `156`, `704`, `185`, `157`; level-0 roots
  `856` ("Not Found") and `270` ("Disrupted reading frame /"); mainroles `hypothetical_proteins`,
  `unknown_function`, `unclassified`. `141` / `703` deliberately unflagged (class known).

## Verification

Unit: `pytest -m "not slow and not kg"`. Live (after the next Docker rebuild): the three new kg tests
above, `pytest tests/kg_validity -q`, and `capture_annotation_state.py --compare` (expected: 0 bucket
movement — neither change touches the 9 source buckets).
