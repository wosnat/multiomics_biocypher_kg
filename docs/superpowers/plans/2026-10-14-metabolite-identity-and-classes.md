# Metabolite Identity + KEGG Compound Classes — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** One `Metabolite` node per metabolite (canonical KEGG → primary ChEBI → MNX identity, primary ChEBI ids, no role terms, no class nodes carrying a monomer structure), deterministic paper-metabolite resolution with recorded stereo assumptions, a curated `MetaboliteClass` layer from KEGG compound BRITE, a `ChemicalRole` hierarchy for role-typed TCDB substrates, and a pinned TCDB substrate-lumping set — all decided in prepare-data step 6 (step 7 folded in).

**Architecture:** Every identity / stereo / role / class decision is made in step 6 (`build_kegg_metabolism_xrefs.py`) by pure functions in a new package `multiomics_kg/utils/metabolite_identity/` plus a ChEBI loader `multiomics_kg/utils/chebi_utils.py`. Step 6 runs stages A–P in memory and writes every output file at once in stage P, only after validation passes. Adapters only copy fields; post-import only adds rollups and indexes.

**Tech Stack:** Python 3 (uv), pandas, PyYAML, sqlite3 (MNX resolver), `chemparse`, BioCypher, Neo4j 5 Cypher (`scripts/post-import.sh` + mirror `.cypher`), pytest.

**Spec:** `docs/superpowers/specs/2026-10-13-metabolite-identity-and-classes-design.md` (signed off 2026-10-14). **Explorer contract:** the top section (Summary → §2b) of `docs/kg-changes/metabolite-identity-and-classes.md`. Executors read the spec section each task cites; this plan does not restate rationale.

## Global Constraints

- Every string property passes `_clean_str` (`'` → `^`, `|` removed) before an adapter yields it; never use `|` or `'` as a literal separator in computed strings (CLAUDE.md, "String property sanitization").
- R1: values the KG mints are lowercase `snake_case`; external names (ChEBI role names, KEGG BRITE names) are verbatim.
- R2: every `sources` value joins a `DataSource` node as `d.id = 'data_source:' + s`.
- R5: no native `bool` property; two-state facts are categorical strings.
- A closed controlled vocabulary describes all of its values or none (`value_descriptions`).
- `scripts/post-import.sh` and `scripts/post-import.cypher` carry identical Cypher logic.
- `normalize_curie` returns `None` for `chebi.role:` and `mnx:` (bioregistry, measured 2026-10-14). Build `chebi.role:<n>` and `kegg.brite:br08xxx…` ids with plain f-strings, never through `normalize_curie` (the `tigr.role:` precedent in `ncbifam_adapter.py:78`).
- Step 6 is deterministic: sorted keys and lists, `json.dumps(..., indent=2, sort_keys=True)`; the only allowed byte difference between two runs is `frozen_on` on newly frozen entries.
- A failed step-6 run writes **no** file (spec §5.2 stage P).
- Commit from this clone (`second_multiomics`); Docker rebuilds happen in `/home/osnat/github/multiomics_biocypher_kg` after a push over SSH: `git push git@github.com:wosnat/multiomics_biocypher_kg.git HEAD:main`. Never run `docker compose` in this clone.
- Unit-test command: `uv run pytest -m "not slow and not kg" -q`. KG tests: `uv run pytest -m kg -q` (needs the live graph at `bolt://localhost:7687`).

## Review Focus

1. **A paper alias file referenced by only some of a paper's tables.** Kujawinski 2023 references `metabolite_aliases.yaml` from 1 of 7 tables; the other 6 use `id_col: KEGG`. A reasonable person expects a curated override (GSH → `C00051`) to win on every table of that paper. The plan defines a paper's overrides as the union of every `aliases_file` its tables reference (Task 7 test `test_overrides_apply_to_every_table_of_the_paper`).
2. **Element symbols that start with `R`.** `Rb` (rubidium cation `C17061`, `mass` 85.468) must stay `defined`, with its mass, while `R` / `Rn` in a generic formula is a placeholder. The naive regex the live-graph invariant uses misfires on `Rb` (measured). Task 2 tests `C17061`-style `Rb`, `Ru`, `Re` and `C5H8O2Rn`.
3. **A case-significant name that is also an acronym.** `CO` vs `Co`, `NO` vs `No` stay distinct; `NAD` / `nad` fold; `AMP` is formula-like and does not fold (harmless). Task 2 pins all of them; Task 7 pins that the resolver skips its lower-cased name step for them.
4. **Re-running step 6 after a successful run.** A second run with the committed lock and pin must produce byte-identical files and must not re-freeze (`frozen_on` unchanged). Task 12 test `test_rerun_with_committed_lock_and_pin_is_byte_identical`; Task 14 Step 6 repeats it on the real caches.
5. **A dropped measurement-only node that is really the same entity as its target** (SAM `chebi:142094` vs `C00019`). It must become a merge (id in `all_identifiers`), never a drop (id in `recurated_from_ids`). Task 14 step "same-entity review of every drop" plus the KG test `test_recurated_ids_are_not_identifiers`.

---

## Measured baseline (re-measured 2026-10-14 for this plan)

Sources: live KG (`Schema_info.version 0.0.0-dev`, the 2026-10-13 build), committed caches at `54689c7c`, ChEBI release **256** (`chebi.obo.gz`, data-version 256, 219,012 terms), KEGG REST `br:br08001…br08021` fetched 2026-10-14. Scratch scripts were run from the session scratchpad and are not committed. The step-6 identity report is the authoritative number after Task 14; these are the plan's expectations.

| Spec claim | Spec value | Re-measured | Note | After rebuild (2026-10-15) |
|---|---|---|---| --- |
| `Metabolite` nodes | 3,356 (2,717 KEGG + 635 `chebi:` + 4 `mnx:`) | **3,356** (2,717 + 635 + 4) | live | **3,203** (2,776 + 427 + 0) |
| ChEBI ids that are secondary | 1,627 (1,454 on KEGG + 173 node ids) | **1,628** (1,455 + 173) | 0 ids unknown to ChEBI 256 | all node ids primary; 79 remap rows `chebi_primary_id` |
| KEGG nodes whose `chebi_id` disagrees with ChEBI's KEGG xref | 331 | **331** | of 2,427 KEGG nodes ChEBI cross-references | `chebi_id` = ChEBI's KEGG xref entity |
| KEGG nodes with > 1 ChEBI entity term | 375 | **375** | | one `chebi_id` (+ `chebi_ids`) |
| Role `Metabolite` nodes / transport edges / families | 57 / 392 / 116 | **57 / 392 / 116** | 8 of the 57 carry a secondary id | **0** role metabolites; 392 role edges / 116 families on `Tcdb_family_transports_chemical_role` |
| Max roles per TC node | 33 | **33** | | **33** (role edges per TC node) |
| `ChemicalRole` hierarchy (no `role` root) | ≈ 111 | **111** | 3 branches present; 132 `is_a` edges inside the set; levels 0:3 · 1:28 · 2:33 · 3:29 · 4:12 · 5:5 · 6:1 | **111** (levels 0:3 · 1:28 · 2:33 · 3:29 · 4:12 · 5:5 · 6:1), 132 `is_a` |
| Multi-branch TCDB roles | 10 | **10** | `siderophore` = `CHEBI:26672`, `phytosiderophore` = `CHEBI:38155` | `chebi.role:26672` reached by 7 most_specific + 11 inherited families |
| Generic-formula nodes | 420 (live regex) / 422 (step-6 rule) | **420** | with `R` not followed by `[abefghu]`; the naive regex counts `Rb` as generic | `structure_kind = generic` **622** |
| `mass` on a generic formula | 0 | **0** | the 1 naive hit is rubidium `C17061` (`Rb`) | **0** |
| Formula-less non-role nodes | 316 | **316** | | **244** formula-less (all nodes) |
| … that gain a formula from ChEBI | ≈ 49 (44 transport) | **52 (44 transport)** | 8 generic (`R` / `*` / `)n`), 44 concrete (e.g. `chebi:12923` D-fructose); most concrete ones are re-keyed or merged first, so the final figure is lower | 33 filled from ChEBI |
| InChIKey on two nodes | 5 (B2) | **22 pairs** = 17 MNX same-entity pairs + the 5 B2 pairs | the MNX pairs disappear by merge | **0** duplicate InChIKeys |
| MNX same-entity pairs | 18 | **18** | | merged |
| Merge groups / merged away / re-keyed | 69 / 74 / 89 | **71 / 77 / 88** (scratch union-find, approximate) | the scratch run applied conjugate/tautomer links without the secondary-id nuance; Task 14 reports the real numbers | 40 merge groups; 84 old ids collapse into existing nodes; 124 `rekey_chebi_kegg_xref` (49 into existing / 75 new) |
| Guard rejections (ChEBI-only KEGG↔KEGG) | 4 | **8 distinct pairs** | conjugate: water/hydroxide `C00001`/`C01328`, cyanide `C00177`/`C01326`; shared ChEBI term: acetoin `C00466`/`C01769`, quinone `C00472`/`C15602`, hydroquinone `C00530`/`C15603`; tautomer: `C04287`/`C06891`, `C05601`/`C06201`, `C15650`/`C15651` | 10 `rejected_links` in the identity report |
| Secondary `chebi:` node ids that only change to their primary id (no role, merge or re-key) | not stated | **60** | needs a remap basis; see "Spec gaps" | 79 `chebi_primary_id` remap rows (incl. fallback relocations) |
| Stereo `is_a` links / enantiomer pairs | ≈ 60 / ≈ 25 | **73 / 29** (before the B2 scrub and stage K) | scrub + drop reduce them; Task 14 reports the final counts | **68 / 28** |
| Paper metabolites / names / direct ids | 149 / 211 / 64 | **149 / 211 / 64** | methods: `kegg_direct` 64 · `ambiguous_multi_id` 36 · `name_match` 35 · `alias_override` 14 | 134 `metabolomics`-tagged nodes |
| Ambiguous paper names (B5) | 36 | **36** | list in Task 14 | 0 ambiguous; 114 repoint rows |
| Unresolved paper names | 1 (6 rows) | **1** (`s-(1-2-dicarboxyethyl)glutathione pos`) | | **1** (6 rows) |
| Dominant-form moves | 7, smallest margin ≈ 5× | **7**, smallest 520 / 103 = 5.0× (L-/D-Ala) | margins use `catalyst_gene_count` (see "Spec gaps") | **7**, all onto L-forms |
| Same-slot duplicate assay edges | 6 | **6** (all on `C00127`, 4 Biller + 2 Kujawinski) | | split onto GSH `C00051` / GSSG `C00127` |
| BRITE: metabolites / class nodes / memberships | 867 / ≈ 431 / ≈ 1,130 | **867 / 431 / 1,130** | per tree identical to spec §3.7; class `is_a` edges **395** (spec ≈ 420); br08001 level-0 = 9 names; max depth 4 levels | **881** classified metabolites / **445** / **1,151**; class `is_a` **408** |
| Lumping set | 13 ids | **13**, identical list | live `level >= 2 AND metabolite_count >= 50` | **13**, identical (`3.A.3` at 47, held by the pin) |
| `transport_substrate_resolution` | — | `resolved` **28,854** · `family_inferred` **1,693** | the post-rebuild check must equal this | **identical** (28,854 · 1,693) |
| `DataSource` nodes | 9 | **9** | → 11 | **11** |
| Assay edges | — | `Assay_quantifies_metabolite` **1,480** · `Assay_flags_metabolite` **186** | | **1,480 · 186** (unchanged) |
| `Reaction_has_metabolite` / `Organism_has_metabolite` / `Metabolite_in_pathway` | — | **10,149 / 92,081 / 9,269** | | **10,149 / 88,320 / 9,342** |
| `Tcdb_family_transports_metabolite` | — | **11,263** | expect −392 roles, minus merge dedup | **10,790** (−392 roles, −81 net dedup) |
| Pathway-level `KeggTerm` / with `gene_count = 0` | — | **447 / 61** | the post-rebuild doc reports the delta | **452 / 66** (+5: `ko01521`, `ko04145`, `ko04510`, `ko05212`, `ko05332`) |

## Spec gaps resolved in this plan (no design change; flagged for the user)

1. **Dominant-form count.** Spec §3.6 says "the most gene-reachable catalysing reactions", but its measured numbers and the lock example (`alanine`: winner 520, runner-up 103; Glu 1,905 vs 69) are **`catalyst_gene_count`**, i.e. distinct catalysing genes. Reaction counts would be L-Glu 86 vs Glu 2. The plan uses **distinct catalysing genes**, computed in step 6 from `gene_annotations_merged.json` `kegg_reactions`, the same set post-import's `Metabolite.catalyst_gene_count` counts.
2. **Remap basis for plain secondary → primary ChEBI ids.** 60 `chebi:` node ids change only because the id was secondary. Spec §12's closed basis list has no value for them. The plan adds **`chebi_primary_id`** to the `metabolite_id_remap.tsv` `basis` column (a release-file column, not a graph vocabulary, spec §7.3).
3. **`DataSource` registry for `kegg` / `chebi`** (spec §8.5 "decided in the plan"). A new top-level `reference_sources:` list in `config/gene_annotations_config.yaml`, beside `sources:`. The merge code reads only `sources` / `fields` (`build_gene_annotations.py:879,1057`), so it is unaffected. `DataSourceAdapter` and the two `sources`-pinning tests read the union.
4. **ChEBI release into `Schema_info`** (spec §8.6). Same pattern as the vocabulary hash: `create_knowledge_graph.py` writes `chebi_release.txt` beside `controlled_vocabularies.sha256`; post-import Group 4 reads it into `Schema_info.chebi_release`. `data_source:chebi.version` also carries it.
5. **Anomer detection** (spec §3.6 "settled in the plan"). ChEBI 256 has **no** anomer relation (its typedefs are `has part`, `has role`, conjugate base/acid, tautomer, substituent group, functional parent, enantiomer, parent hydride). So two members are anomers iff their names (KEGG name, else ChEBI name) are equal after removing `alpha-` / `beta-` / `α-` / `β-` tokens **and** at least one carries a ChEBI ancestor in {`CHEBI:16646` carbohydrate, `CHEBI:63299` carbohydrate derivative, `CHEBI:24400` glycoside, `CHEBI:78616` carbohydrates and carbohydrate derivatives}. Verified: `C00031` D-Glucose / `C00221` beta-D-Glucose, `C00140` / `C03878` (GlcNAc), MurNAc all sit under `CHEBI:63299` or `CHEBI:16646`.
6. **`structure_kind` split rule** (spec §3.3 "settled in the plan"), in this order: `polymer` if the formula has a repeat marker (`)n`, `)m`, `)x`, a trailing `n`), or the name ends in `)n`. Else `generic` if the formula has `*`, a placeholder `R` (`R` not followed by one of `a b e f g h u`), any other non-element symbol (`A`, `X`), or the node's ChEBI term is a class (spec §3.3: "or the ChEBI term is a class"). Else `defined`. KEGG nodes are classified against the non-class ChEBI term `pick_chebi_id` would choose, so a KEGG compound never turns `generic` through a class cross-reference.
7. **ChEBI "class" term** (spec §6.4 `is_class`): no InChIKey, **and** either no generalized formula or a generic one (rule 6).
8. **`DataSource.scope` / `provenance` for `kegg` / `chebi`.** Spec §8.5 says "scope global, provenance reference_db", copied from CLAUDE.md, but the live values are `scope ∈ {gene_level, organism_restricted}` and `provenance ∈ {download, tool_run}`; neither `global` nor `reference_db` exists. The plan uses **`scope: global`** (new, because neither existing value describes a non-gene source) and **`provenance: download`** (existing; both are downloaded reference data), and corrects CLAUDE.md to the values actually in use.
9. **Bare-id full-text search** (spec §9.2 "if the analyser doesn't produce bare tokens …"). Lucene's standard tokenizer (UAX#29) does not break at `:` or `.` between letters and digits, so `kegg.compound:C01353` is **one** token and a bare `C01353` query would miss it. The plan therefore takes the spec's fallback up front: a post-import, search-only `Metabolite.identifier_tokens: str[]` (bare local parts of `all_identifiers`), and `metaboliteFullText` indexes `[name, formula, identifier_tokens]`. The KG test `test_fulltext_hits` still checks the three bare ids.

---

## File Structure

**Create**

| File | Responsibility |
|---|---|
| `multiomics_kg/utils/chebi_utils.py` | `chebi.obo.gz` download + parse; derived cache build; `ChebiIndex` lookups (primary id, KEGG xrefs, ancestry, role / class / carbohydrate tests) |
| `multiomics_kg/utils/metabolite_identity/__init__.py` | Re-exports the public names below, so `metabolite_identity.name_key` / `metabolite_identity.write_resolved_tables` (spec names) work |
| `…/metabolite_identity/names.py` | `name_key`, periodic table, formula classification, `parse_elements`, stereo-descriptor test, anomer-token strip, HMDB normalisation |
| `…/metabolite_identity/canonical.py` | Merge links, union-find with the KEGG guard and overrides, winner rule, id rewrite (stages E–F) |
| `…/metabolite_identity/structure.py` | `structure_kind`, B2 scrub, formula fill, mass rule (stage G) |
| `…/metabolite_identity/crossrefs.py` | `chebi_id` pick, id arrays, `all_identifiers`, `identity_basis`, `name` (stage H) |
| `…/metabolite_identity/stereo.py` | Stereo links, families, anomer test, final stereo lists (stages I, L) |
| `…/metabolite_identity/paper.py` | Alias-file loading, paper-table harvest, per-row resolution, dominant-form rule + lock (stages D-harvest, J) |
| `…/metabolite_identity/roles.py` | Role partition, pruned role hierarchy, `subtree_roles` (stage C) |
| `…/metabolite_identity/classes.py` | KEGG compound BRITE fetch, parse, canonical prune (stage M) |
| `…/metabolite_identity/lumping.py` | Lumping pin seed / drift / auto-append (stage O) |
| `…/metabolite_identity/outputs.py` | Stage-P validation, atomic multi-file write, `_resolved.csv` text, mapping v2, identity report, remap / repoint rows |
| `…/metabolite_identity/pipeline.py` | `Step6Inputs` → `run_identity_stages` (stages E–M, O, validation) → `Step6Outputs` → `render_outputs` (file texts). Pure: no MNX, no network, no writes |
| `multiomics_kg/adapters/metabolite_class_adapter.py` | `MetaboliteClass` nodes + membership / is-a edges |
| `config/metabolite_identity_overrides.yaml` | Curated `merge:` / `keep_separate:` (starts empty) |
| `scripts/metabolite_identity_tsvs.py` | Writes the two release TSVs from the identity report |
| `tests/metabolite_identity_helpers.py` | `make_chebi(...)` fixture builder shared by the unit tests |
| `tests/test_chebi_utils.py`, `tests/test_mi_names.py`, `tests/test_mi_canonical.py`, `tests/test_mi_structure.py`, `tests/test_mi_crossrefs.py`, `tests/test_mi_stereo.py`, `tests/test_mi_paper.py`, `tests/test_mi_roles.py`, `tests/test_mi_classes.py`, `tests/test_mi_lumping.py`, `tests/test_mi_outputs.py`, `tests/test_step6_identity_pipeline.py`, `tests/test_metabolite_class_adapter.py` | Unit tests (no network) |
| `tests/fixtures/chebi_mini.obo` | 20-term obo fixture for the parser test |
| `tests/kg_validity/test_metabolite_identity.py` | Live-graph invariants (spec §11 KG validity) |
| `tests/kg_validity/fixtures/metabolite_elements_pre_identity.json` | Pre-change `{id: elements}` export for the "no node loses `elements`" test |

**Modify**

| File | Change |
|---|---|
| `multiomics_kg/download/build_kegg_metabolism_xrefs.py` | Stages A–P; `_resolve_substrate_strings` gains ChEBI re-key; harvest returns candidates; atomic write |
| `multiomics_kg/utils/metabolite_utils.py` | `resolve_metabolite_candidates` (all candidates; skips lower-cased name step for case-significant names) |
| `multiomics_kg/download/resolve_paper_metabolites.py` | Retired: prints "folded into step 6", exits 0 |
| `scripts/prepare_data.sh` | Step-7 alias, `STEPS` / `--rebuild` lists, header comments |
| `config/schema_config.yaml` | New nodes / edges / properties (spec §7) |
| `config/controlled_vocabularies.yaml` | New entries (spec §7.3) |
| `config/gene_annotations_config.yaml` | `reference_sources:` (`kegg`, `chebi`) |
| `multiomics_kg/adapters/data_source_adapter.py` | Read `reference_sources`; names / descriptions; `chebi` version |
| `multiomics_kg/adapters/metabolism_adapter.py` | New fields; `_parse_elements` allow-list; stereo + enantiomer edges; `all_metabolite_ids()` |
| `multiomics_kg/adapters/tcdb_adapter.py` | `ChemicalRole` nodes + edges; depth refactor; `substrate_breadth` |
| `multiomics_kg/adapters/metabolite_assay_adapter.py` | `stereo_assignment` pass-through; docstring |
| `create_knowledge_graph.py` | `MetaboliteClassAdapter`; `chebi_release.txt` |
| `scripts/post-import.sh`, `scripts/post-import.cypher` | Indexes, rollups, `transport_substrate_resolution`, `Schema_info.chebi_release` |
| `scripts/validate_paperconfig.py` | Extended alias-file validation |
| `data/Prochlorococcus/papers_and_supp/{biller 2022,Capovilla 2023,Kujawinski 2023}/metabolite_aliases.yaml` | Explicit `{id, stereo_assignment, reason}` form + new curated entries |
| `.claude/skills/{paperconfig,add-a-strain,cypher-queries,omics-edge-snapshot}/…` | Spec §14 |
| `CLAUDE.md`, `CHANGELOG.md`, `plans/backlog.md`, `docs/kg-changes/{vocabulary-contract,tcdb-two-source-upgrade,metabolomics-extension,metabolite-identity-and-classes}.md`, `.gitignore` | Spec §14 |
| Tests that pin old behaviour | `tests/test_data_source_adapter.py`, `tests/kg_validity/test_data_source.py`, `tests/test_tcdb_adapter.py` (`…threshold_is_inlined`), `tests/test_controlled_vocab.py`, `tests/test_annotation_provenance.py`, `tests/test_resolve_paper_metabolites.py` (ported, then deleted), `tests/test_build_kegg_metabolism_xrefs_metabolomics.py`, `tests/kg_validity/snapshot_data.json` (regenerated) |

**Shared in-memory data model (stages D–P).** Step 6 converts `kegg_data["compounds"]` (bare `C…` keys) and `kegg_data["additional_compounds"]` (CURIE keys) into one dict **`table: dict[str, dict]` keyed by full node id** (`kegg.compound:C00244`, `chebi:17632`, `mnx:MNXM…`) right after `_fold_substrates_into_kegg_data`. Each entry has the existing compound fields plus `evidence_sources: list[str]`. Stage P converts back. Every `metabolite_identity` function below takes / returns this keyed-by-node-id shape.

---

## Phase 1 — ChEBI cache + `metabolite_identity` (pure functions, unit-tested)

### Task 1: ChEBI loader and derived cache (`chebi_utils.py`)

**Files:**
- Create: `multiomics_kg/utils/chebi_utils.py`, `tests/test_chebi_utils.py`, `tests/fixtures/chebi_mini.obo`, `tests/metabolite_identity_helpers.py`
- Modify: `.gitignore`

**Interfaces:**
- Consumes: `metabolite_identity.names.is_generic_formula` / `is_polymer_formula` (Task 2), used by `is_class`. **Do Task 2 before this task's Step 5** (or implement both in one session); Task 2 has no dependency on Task 1.
- Produces:
  - `ROLE_ROOT = "CHEBI:50906"`, `ROLE_BRANCHES: frozenset[str]` (`CHEBI:24432` biological role, `CHEBI:51086` chemical role, `CHEBI:33232` application), `CARBOHYDRATE_ANCESTORS: frozenset[str]` (`CHEBI:16646`, `CHEBI:63299`, `CHEBI:24400`, `CHEBI:78616`), `CHEBI_OBO_URL`.
  - `parse_chebi_obo(path: Path) -> tuple[str, dict[str, dict]]` → `(release, terms)`; terms keyed `"CHEBI:<n>"` with `{name, charge, alt_ids, kegg, conjugates, tautomers, enantiomers, is_a, generalized_formula, inchikey}`.
  - `class ChebiIndex(terms: dict[str, dict], release: str)` with `release`, `touched: set[str]`, `missing: set[str]`, and methods `primary(cid) -> str | None` (accepts `CHEBI:n` / `chebi:n` / `n`; maps secondary → primary; records touched / missing), `term(cid) -> dict | None`, `kegg_terms(kegg_id: str) -> list[str]` (sorted primary ids whose term cross-references `kegg_id`), `ancestors(cid) -> frozenset[str]` (transitive `is_a`), `is_role(cid) -> bool`, `is_class(cid) -> bool`, `is_carbohydrate(cid) -> bool`, `secondary_ids(cid) -> list[str]`, `charge(cid) -> int | None`.
  - `build_chebi_cache(index: ChebiIndex) -> dict` → `{"release", "terms": {...}}` restricted to spec §6.4 coverage: every term with a KEGG xref ∪ `index.touched` ∪ their `is_a` ancestors ∪ conjugate / tautomer / enantiomer partners. Cached term adds `is_role` / `is_class` booleans **inside the JSON cache only** (not a graph property, so R5 doesn't apply).
  - `index_from_cache(cache: dict) -> ChebiIndex`.
  - `download_chebi_obo(raw_dir: Path, force: bool = False) -> Path`.

- [ ] **Step 1: Write the obo fixture.** `tests/fixtures/chebi_mini.obo` holds, in obo 1.2 syntax with `data-version: 256` in the header:
  - `CHEBI:17632` nitrate (`alt_id: CHEBI:14654`, no KEGG xref), charge −1, `relationship: RO:0018033 CHEBI:48107`;
  - `CHEBI:48107` nitric acid, charge 0, `xref: kegg.compound:C00244`, `relationship: RO:0018034 CHEBI:17632`;
  - `CHEBI:50906` role; `CHEBI:24432` biological role `is_a CHEBI:50906`; `CHEBI:51086` chemical role `is_a CHEBI:50906`; `CHEBI:33232` application `is_a CHEBI:50906`;
  - `CHEBI:26672` siderophore `is_a CHEBI:24432`; `CHEBI:38155` phytosiderophore `is_a CHEBI:26672`;
  - `CHEBI:33229` vitamin (role) `alt_id: CHEBI:10004`, `is_a CHEBI:24432`, `is_a CHEBI:33232`;
  - `CHEBI:16977` L-alanine (`xref: kegg.compound:C00041`, `relationship: RO:0018039 CHEBI:15570`, `is_a CHEBI:16449`, inchikey `QNAYBMKLOCPYGJ-REOHCLBHSA-N`);
  - `CHEBI:15570` D-alanine (`xref: kegg.compound:C00133`, enantiomer of 16977, `is_a CHEBI:16449`, inchikey `QNAYBMKLOCPYGJ-UWTATZPHSA-N`);
  - `CHEBI:16449` alanine (`xref: kegg.compound:C01401`, inchikey `QNAYBMKLOCPYGJ-UHFFFAOYSA-N`);
  - `CHEBI:16646` carbohydrate; `CHEBI:4167` D-glucopyranose `is_a CHEBI:16646`, `xref: kegg.compound:C00031`; `CHEBI:15903` beta-D-glucose `is_a CHEBI:4167`, `xref: kegg.compound:C00221`;
  - `CHEBI:17754` glycerol (`xref: kegg.compound:C00116`); `CHEBI:17522` alditol (no inchikey, `generalized_empirical_formula "C2H6O2(CH2O)n"`).

  Write every `property_value` line as `property_value: chemrof:<prop> "<value>" xsd:string` (charge as `xsd:integer`), exactly the release-256 syntax.

- [ ] **Step 2: Write the shared test helper.** `tests/metabolite_identity_helpers.py`:

```python
"""Shared builders for metabolite-identity unit tests (no network, no MNX)."""
from multiomics_kg.utils.chebi_utils import ChebiIndex


def term(name, *, charge=None, alt=(), kegg=(), conj=(), taut=(), enan=(), is_a=(),
         gformula=None, inchikey=None):
    return {"name": name, "charge": charge, "alt_ids": list(alt), "kegg": list(kegg),
            "conjugates": list(conj), "tautomers": list(taut), "enantiomers": list(enan),
            "is_a": list(is_a), "generalized_formula": gformula, "inchikey": inchikey}


def make_chebi(**terms) -> ChebiIndex:
    """make_chebi(CHEBI_17632=term(...)) — underscores become colons."""
    return ChebiIndex({k.replace("_", ":"): v for k, v in terms.items()}, release="256")


def cpd(name, *, formula=None, mass=None, inchikey=None, smiles=None, chebi_id=None,
        mnxm_id=None, hmdb_id=None, ev=("metabolism",), pathways=()):
    return {"name": name, "formula": formula, "mass": mass, "inchikey": inchikey,
            "smiles": smiles, "chebi_id": chebi_id, "mnxm_id": mnxm_id, "hmdb_id": hmdb_id,
            "evidence_sources": list(ev), "pathways": list(pathways)}
```

- [ ] **Step 3: Write the failing tests** (`tests/test_chebi_utils.py`):

```python
from pathlib import Path

from multiomics_kg.utils.chebi_utils import (
    build_chebi_cache, index_from_cache, parse_chebi_obo, ChebiIndex)

FIX = Path(__file__).parent / "fixtures" / "chebi_mini.obo"


def _index():
    release, terms = parse_chebi_obo(FIX)
    return ChebiIndex(terms, release)


def test_parse_release_and_fields():
    release, terms = parse_chebi_obo(FIX)
    assert release == "256"
    t = terms["CHEBI:48107"]
    assert t["kegg"] == ["C00244"] and t["charge"] == 0
    assert terms["CHEBI:17632"]["conjugates"] == ["CHEBI:48107"]
    assert terms["CHEBI:16977"]["enantiomers"] == ["CHEBI:15570"]
    assert terms["CHEBI:17522"]["generalized_formula"] == "C2H6O2(CH2O)n"


def test_secondary_maps_to_primary_and_is_recorded():
    ix = _index()
    assert ix.primary("chebi:10004") == "CHEBI:33229"
    assert ix.primary("10004") == "CHEBI:33229"
    assert ix.primary("CHEBI:999999") is None
    assert "CHEBI:33229" in ix.touched and "CHEBI:999999" in ix.missing
    assert ix.secondary_ids("CHEBI:33229") == ["CHEBI:10004"]


def test_kegg_terms_and_roles():
    ix = _index()
    assert ix.kegg_terms("C00244") == ["CHEBI:48107"]
    assert ix.is_role("CHEBI:38155") and ix.is_role("chebi:10004")
    assert not ix.is_role("CHEBI:16977")
    assert ix.ancestors("CHEBI:38155") >= {"CHEBI:26672", "CHEBI:24432", "CHEBI:50906"}


def test_class_and_carbohydrate():
    ix = _index()
    assert ix.is_class("CHEBI:17522")          # alditol: no inchikey, generic formula
    assert not ix.is_class("CHEBI:16977")
    assert ix.is_carbohydrate("CHEBI:15903")


def test_cache_coverage_and_round_trip():
    ix = _index()
    ix.primary("chebi:10004")                   # touched → kept with ancestors
    cache = build_chebi_cache(ix)
    assert cache["release"] == "256"
    kept = set(cache["terms"])
    assert {"CHEBI:48107", "CHEBI:17632"} <= kept      # KEGG xref + its conjugate
    assert {"CHEBI:33229", "CHEBI:24432", "CHEBI:50906"} <= kept
    assert cache["terms"]["CHEBI:33229"]["is_role"] is True
    again = index_from_cache(cache)
    assert again.primary("10004") == "CHEBI:33229" and again.release == "256"
```

- [ ] **Step 4: Run to verify failure.** `uv run pytest tests/test_chebi_utils.py -q` → FAIL (`ModuleNotFoundError: multiomics_kg.utils.chebi_utils`).

- [ ] **Step 5: Implement `multiomics_kg/utils/chebi_utils.py`.**

```python
"""ChEBI ontology loader for step 6 (metabolite identity).

Reads chebi.obo.gz (release-tagged) into a term dict, answers the lookups the
identity pipeline needs, and builds the committed derived cache
cache/data/chebi/chebi_metabolite_xrefs.json (spec §6.4). Read by step 6 only.
ChEBI is CC BY 4.0.
"""
from __future__ import annotations

import gzip
import logging
import re
from pathlib import Path

import requests

log = logging.getLogger(__name__)

CHEBI_OBO_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz"
ROLE_ROOT = "CHEBI:50906"
ROLE_BRANCHES = frozenset({"CHEBI:24432", "CHEBI:51086", "CHEBI:33232"})
CARBOHYDRATE_ANCESTORS = frozenset({"CHEBI:16646", "CHEBI:63299", "CHEBI:24400", "CHEBI:78616"})

_REL = {"RO:0018033": "conjugates", "RO:0018034": "conjugates",
        "RO:0018036": "tautomers", "RO:0018039": "enantiomers"}
_PROP = re.compile(r'chemrof:(\S+) "(.*)" xsd:\w+')


def _new_term() -> dict:
    return {"name": None, "charge": None, "alt_ids": [], "kegg": [], "conjugates": [],
            "tautomers": [], "enantiomers": [], "is_a": [], "generalized_formula": None,
            "inchikey": None}


def parse_chebi_obo(path: Path) -> tuple[str, dict[str, dict]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    release, terms, cur, obsolete = "", {}, None, False
    with opener(path, "rt", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("data-version:"):
                release = line.split(":", 1)[1].strip()
            if line == "[Term]":
                cur, obsolete = _new_term(), False
                continue
            if line.startswith("["):
                cur = None
                continue
            if cur is None or ": " not in line:
                continue
            key, _, val = line.partition(": ")
            if key == "id":
                terms[val] = cur
            elif key == "name":
                cur["name"] = val
            elif key == "alt_id":
                cur["alt_ids"].append(val)
            elif key == "is_a":
                cur["is_a"].append(val.split(" ")[0])
            elif key == "is_obsolete" and val == "true":
                obsolete = True
                cur["obsolete"] = True
            elif key == "xref" and val.startswith("kegg.compound:"):
                cur["kegg"].append(val.split(" ")[0].split(":", 1)[1])
            elif key == "relationship":
                rel, tgt = val.split(" ")[:2]
                if rel in _REL:
                    cur[_REL[rel]].append(tgt)
            elif key == "property_value":
                m = _PROP.match(val)
                if m:
                    prop, pval = m.groups()
                    if prop == "generalized_empirical_formula":
                        cur["generalized_formula"] = pval
                    elif prop == "charge":
                        cur["charge"] = int(pval)
                    elif prop == "inchi_key_string":
                        cur["inchikey"] = pval
    terms = {k: v for k, v in terms.items() if not v.pop("obsolete", False)}
    for t in terms.values():
        for k in ("alt_ids", "kegg", "conjugates", "tautomers", "enantiomers", "is_a"):
            t[k] = sorted(set(t[k]))
    if not release:
        raise ValueError(f"{path}: no data-version header (ChEBI format drift?)")
    if not terms:
        raise ValueError(f"{path}: parsed 0 terms (ChEBI format drift?)")
    return release, terms


def _generic(formula: str | None) -> bool:
    # Local import keeps chebi_utils free of a package-level cycle.
    from multiomics_kg.utils.metabolite_identity.names import is_generic_formula, is_polymer_formula
    return formula is not None and (is_generic_formula(formula) or is_polymer_formula(formula))


class ChebiIndex:
    def __init__(self, terms: dict[str, dict], release: str):
        self.release = release
        self._terms = terms
        self._alt = {a: p for p, t in terms.items() for a in t.get("alt_ids", [])}
        self._kegg: dict[str, list[str]] = {}
        for cid, t in terms.items():
            for k in t.get("kegg", []):
                self._kegg.setdefault(k, []).append(cid)
        for v in self._kegg.values():
            v.sort(key=lambda c: int(c.split(":")[1]))
        self._anc: dict[str, frozenset[str]] = {}
        self.touched: set[str] = set()
        self.missing: set[str] = set()

    @staticmethod
    def _norm(cid: str) -> str:
        cid = cid.strip()
        local = cid.split(":", 1)[1] if ":" in cid else cid
        return f"CHEBI:{int(local)}" if local.isdigit() else f"CHEBI:{local}"

    def primary(self, cid: str) -> str | None:
        c = self._norm(cid)
        p = c if c in self._terms else self._alt.get(c)
        if p is None:
            self.missing.add(c)
        else:
            self.touched.add(p)
        return p

    def term(self, cid: str) -> dict | None:
        p = self.primary(cid)
        return self._terms.get(p) if p else None

    def kegg_terms(self, kegg_id: str) -> list[str]:
        out = self._kegg.get(kegg_id, [])
        self.touched.update(out)
        return list(out)

    def ancestors(self, cid: str) -> frozenset[str]:
        p = self.primary(cid)
        if p is None:
            return frozenset()
        if p not in self._anc:
            seen: set[str] = set()
            stack = list(self._terms[p].get("is_a", []))
            while stack:
                x = stack.pop()
                if x in seen:
                    continue
                seen.add(x)
                stack.extend(self._terms.get(x, {}).get("is_a", []))
            self._anc[p] = frozenset(seen)
        return self._anc[p]

    def is_role(self, cid: str) -> bool:
        p = self.primary(cid)
        if p is None:
            return False
        t = self._terms[p]
        if "is_role" in t:                      # loaded from the derived cache
            return t["is_role"]
        return p == ROLE_ROOT or ROLE_ROOT in self.ancestors(p)

    def is_class(self, cid: str) -> bool:
        p = self.primary(cid)
        if p is None:
            return False
        t = self._terms[p]
        if "is_class" in t:
            return t["is_class"]
        return t.get("inchikey") is None and (
            t.get("generalized_formula") is None or _generic(t.get("generalized_formula")))

    def is_carbohydrate(self, cid: str) -> bool:
        p = self.primary(cid)
        return p is not None and bool((self.ancestors(p) | {p}) & CARBOHYDRATE_ANCESTORS)

    def secondary_ids(self, cid: str) -> list[str]:
        p = self.primary(cid)
        return list(self._terms[p].get("alt_ids", [])) if p else []

    def charge(self, cid: str) -> int | None:
        t = self.term(cid)
        return t.get("charge") if t else None

    def all_terms(self) -> dict[str, dict]:
        return self._terms


def build_chebi_cache(index: ChebiIndex) -> dict:
    terms = index.all_terms()
    keep = {c for c, t in terms.items() if t.get("kegg")} | set(index.touched)
    frontier = list(keep)
    while frontier:
        c = frontier.pop()
        t = terms.get(c)
        if t is None:
            continue
        for nxt in t["is_a"] + t["conjugates"] + t["tautomers"] + t["enantiomers"]:
            if nxt in terms and nxt not in keep:
                keep.add(nxt)
                frontier.append(nxt)
    out = {}
    for c in sorted(keep, key=lambda x: int(x.split(":")[1])):
        t = dict(terms[c])
        t["is_role"] = index.is_role(c)
        t["is_class"] = index.is_class(c)
        out[c] = t
    return {"release": index.release, "terms": out}


def index_from_cache(cache: dict) -> ChebiIndex:
    return ChebiIndex(cache["terms"], cache["release"])


def download_chebi_obo(raw_dir: Path, force: bool = False) -> Path:
    raw_dir.mkdir(parents=True, exist_ok=True)
    dest = raw_dir / "chebi.obo.gz"
    if dest.exists() and not force:
        return dest
    tmp = dest.with_suffix(".gz.part")
    with requests.get(CHEBI_OBO_URL, stream=True, timeout=300) as r:
        r.raise_for_status()
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(1 << 20):
                fh.write(chunk)
    tmp.replace(dest)
    return dest
```

- [ ] **Step 6: Add `cache/data/chebi/raw` to `.gitignore`** next to the other `…/raw` lines (around `.gitignore:180`), with a one-line comment "ChEBI obo (step 6); derived chebi_metabolite_xrefs.json one level up is tracked".

- [ ] **Step 7: Run the tests.** `uv run pytest tests/test_chebi_utils.py -q` → PASS. (`is_class` imports from `metabolite_identity.names`, which Task 2 creates; until then, run Task 2 first or temporarily skip `test_class_and_carbohydrate`. Recommended: implement Task 2 Step 3 before this step.)

- [ ] **Step 8: Commit.**

```bash
git add multiomics_kg/utils/chebi_utils.py tests/test_chebi_utils.py tests/fixtures/chebi_mini.obo tests/metabolite_identity_helpers.py .gitignore
git commit -m "feat(metabolite): ChEBI obo loader + derived xref cache (step-6 input)"
```

### Task 2: Names, formulas and elements (`metabolite_identity/names.py`)

**Files:**
- Create: `multiomics_kg/utils/metabolite_identity/__init__.py`, `multiomics_kg/utils/metabolite_identity/names.py`, `tests/test_mi_names.py`

**Interfaces:**
- Produces: `PERIODIC_TABLE: frozenset[str]` (all 118 symbols), `is_formula_like(name: str) -> bool`, `name_key(name: str) -> str`, `is_stereo_qualified(name: str) -> bool`, `strip_anomer_tokens(name: str) -> str`, `is_polymer_formula(formula: str | None) -> bool`, `is_generic_formula(formula: str | None) -> bool`, `parse_elements(formula: str | None) -> list[str]`, `normalize_hmdb(raw: str | None) -> str | None`.

- [ ] **Step 1: Write the failing tests** (`tests/test_mi_names.py`):

```python
import pytest

from multiomics_kg.utils.metabolite_identity import names as N


@pytest.mark.parametrize("a,b", [("Glutamic acid", "glutamic acid"), ("NAD", "nad"),
                                 ("  L-Alanine ", "l-alanine"), ("GABA", "gaba")])
def test_name_key_folds(a, b):
    assert N.name_key(a) == N.name_key(b)


@pytest.mark.parametrize("a,b", [("Co", "CO"), ("NO", "No"), ("CO2", "Co2"), ("HF", "Hf")])
def test_name_key_keeps_formula_case(a, b):
    assert N.name_key(a) != N.name_key(b)


def test_formula_like():
    assert N.is_formula_like("NH4+") and N.is_formula_like("Fe3+") and N.is_formula_like("CO2")
    assert N.is_formula_like("AMP")          # Am + P: harmless false positive, not folded
    assert not N.is_formula_like("NAD") and not N.is_formula_like("DMSP")


@pytest.mark.parametrize("name", [
    "N-Acetyl-D-glucosamine", "S-Adenosyl-L-methionine", "UDP-D-glucose", "Glycyl-L-leucine",
    "(2R,3S)-isocitrate", "(S)-malate", "L-Alanine", "beta-D-glucose", "(R,S)-lactate",
    "(+)-catechin", "(−)-citronellol", "(-)-citronellol", "meso-2,6-diaminopimelate",
    "DL-Methionine", "alpha-L-fucose"])
def test_stereo_qualified(name):
    assert N.is_stereo_qualified(name)


@pytest.mark.parametrize("name", [
    "beta-Citronellol", "alpha-Ketoglutarate", "S-Adenosylmethionine", "Alanine",
    "Glutamate", "R-lipoate", "trans-Aconitate", "cis-Aconitate", "(E)-cinnamate"])
def test_stereo_unqualified(name):
    assert not N.is_stereo_qualified(name)


def test_strip_anomer_tokens():
    assert N.strip_anomer_tokens("beta-D-Glucose") == N.strip_anomer_tokens("D-Glucose")
    assert (N.strip_anomer_tokens("N-Acetyl-beta-D-glucosamine")
            == N.strip_anomer_tokens("N-Acetyl-D-glucosamine"))


@pytest.mark.parametrize("f,generic,polymer", [
    ("C6H12O6", False, False), ("Rb", False, False), ("Ru", False, False), ("ReO4", False, False),
    ("C5H8O2R", True, False), ("C5H8O2Rn", True, False), ("C3H3O6R", True, False),
    ("H2A", True, False), ("RX", True, False), ("C10H17N3O6S*", True, False),
    ("(C6H10O5)n", False, True), ("C12H22O11(C6H10O5)n", False, True), ("(C4H6O2)n", False, True),
    ("C10H17O8PR2(C5H8O5PR)n", True, True)])
def test_formula_classes(f, generic, polymer):
    assert N.is_generic_formula(f) is generic
    assert N.is_polymer_formula(f) is polymer


@pytest.mark.parametrize("f,els", [
    ("C6H12O6", ["C", "H", "O"]), ("Rb", ["Rb"]), ("C5H8O2Rn", ["C", "H", "O"]),
    ("C3H3O6R", ["C", "H", "O"]), ("(C6H10O5)n", ["C", "H", "O"]), ("H2A", ["H"]),
    ("Fe+3", ["Fe"]), (None, []), ("", [])])
def test_parse_elements_allow_list(f, els):
    assert N.parse_elements(f) == els


def test_normalize_hmdb():
    assert N.normalize_hmdb("HMDB02878") == "HMDB0002878"
    assert N.normalize_hmdb("HMDB0002878") == "HMDB0002878"
    assert N.normalize_hmdb("") is None and N.normalize_hmdb(None) is None
```

- [ ] **Step 2: Run** `uv run pytest tests/test_mi_names.py -q` → FAIL (module missing).

- [ ] **Step 3: Implement `names.py`.**

```python
"""Name keys, stereo descriptors and formula rules for metabolite identity (spec §3.3, §3.6)."""
from __future__ import annotations

import re

import chemparse

PERIODIC_TABLE = frozenset("""
H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se
Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb
Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm
Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og
""".split())

# A formula-like name: element symbols with optional counts, optional trailing charge.
_FORMULA_NAME = re.compile(r"^(?:[A-Z][a-z]?\d*)+(?:\d*[+-]|[+-]\d*)?$")
_SYMBOL = re.compile(r"[A-Z][a-z]?")
_WS = re.compile(r"\s+")
# Placeholder R: an R not starting one of Ra Rb Re Rf Rg Rh Ru. 'Rn' reads as R + n (spec §3.3).
_PLACEHOLDER_R = re.compile(r"R(?![abefghu])")
_POLYMER = re.compile(r"\)[nmx]\b|\)[nmx]$|\)[nmx](?=[.(A-Z])|(?<=[0-9A-Za-z])n$")
_DESCRIPTOR = re.compile(
    r"(?:^|(?<=[-\s(,]))"
    r"(?:L-|D-|DL-|meso-"
    r"|\((?:\d*[RS])(?:,\d*[RS])*\)-"
    r"|\([+\-−]\)-"
    r"|(?:alpha|beta|α|β)-(?=[DL]-))")
_ANOMER = re.compile(r"(?:^|(?<=[-\s(,]))(?:alpha|beta|α|β)-", re.IGNORECASE)


def is_formula_like(name: str) -> bool:
    s = name.strip()
    if not _FORMULA_NAME.match(s):
        return False
    return all(sym in PERIODIC_TABLE for sym in _SYMBOL.findall(s))


def name_key(name: str) -> str:
    s = _WS.sub(" ", name.strip())
    return s if is_formula_like(s) else s.casefold()


def is_stereo_qualified(name: str) -> bool:
    return bool(_DESCRIPTOR.search(name))


def strip_anomer_tokens(name: str) -> str:
    return _ANOMER.sub("", name).casefold().strip()


def is_polymer_formula(formula: str | None) -> bool:
    return bool(formula) and bool(_POLYMER.search(formula))


def is_generic_formula(formula: str | None) -> bool:
    if not formula:
        return False
    if "*" in formula or _PLACEHOLDER_R.search(formula):
        return True
    return any(sym not in PERIODIC_TABLE for sym in _SYMBOL.findall(formula))


def parse_elements(formula: str | None) -> list[str]:
    if not formula:
        return []
    cleaned = _PLACEHOLDER_R.sub("", formula).replace("*", "")
    cleaned = re.sub(r"\)[nmx]", ")", cleaned)
    cleaned = re.sub(r"(?<=[0-9A-Za-z])n$", "", cleaned)
    cleaned = re.sub(r"[+-]\d*$|\d*[+-]$", "", cleaned)
    try:
        found = chemparse.parse_formula(cleaned).keys()
    except Exception:
        found = _SYMBOL.findall(cleaned)
    return sorted({e for e in found if e in PERIODIC_TABLE})


def normalize_hmdb(raw: str | None) -> str | None:
    if not raw:
        return None
    m = re.fullmatch(r"HMDB(\d+)", raw.strip())
    return f"HMDB{int(m.group(1)):07d}" if m else None
```

Then `__init__.py`:

```python
"""Metabolite identity pipeline for prepare-data step 6 (spec 2026-10-13)."""
from multiomics_kg.utils.metabolite_identity.names import (  # noqa: F401
    PERIODIC_TABLE, is_formula_like, is_generic_formula, is_polymer_formula,
    is_stereo_qualified, name_key, normalize_hmdb, parse_elements, strip_anomer_tokens)
```

(Later tasks append their public names to this `__init__`.)

- [ ] **Step 4: Run** `uv run pytest tests/test_mi_names.py tests/test_chebi_utils.py -q` → PASS. If a descriptor or formula case fails, fix the regex, not the test: the parametrized cases are spec §3.6 / §3.3 verbatim.

- [ ] **Step 5: Commit.**

```bash
git add multiomics_kg/utils/metabolite_identity tests/test_mi_names.py
git commit -m "feat(metabolite): name keys, stereo descriptors, formula + element rules"
```

### Task 3: Merge links and canonicalisation (`canonical.py`, stages E–F)

**Files:**
- Create: `multiomics_kg/utils/metabolite_identity/canonical.py`, `config/metabolite_identity_overrides.yaml`, `tests/test_mi_canonical.py`

**Interfaces:**
- Consumes: `ChebiIndex` (Task 1).
- Produces:
  - `MERGE_BASES = ("mnx_same_entity", "chebi_kegg_xref", "chebi_conjugate", "chebi_tautomer", "manual")`.
  - `@dataclass(frozen=True) MergeLink(a: str, b: str, basis: str)` with `a < b`.
  - `@dataclass IdentityOverrides(merges: list[MergeLink], keep_separate: set[frozenset[str]])`; `load_identity_overrides(path: Path) -> IdentityOverrides` (missing file or empty lists → empty; a pair in both lists → `ValueError`).
  - `entity_terms(node_id: str, entry: dict, chebi: ChebiIndex) -> set[str]`: KEGG node → `chebi.kegg_terms(C…)`; `chebi:` node → `{primary}`; `mnx:` node → `{primary(entry["chebi_id"])}` if set; else `set()`.
  - `collect_merge_links(table: dict[str, dict], chebi: ChebiIndex) -> list[MergeLink]` (sorted, deduplicated).
  - `@dataclass IdentityResult(canonical_of: dict[str, str], members: dict[str, list[str]], basis: dict[str, list[str]], member_basis: dict[str, str], rejected: list[MergeLink], refused_collisions: list[dict], stale_overrides: list[dict])`. `canonical_of` covers every input node id (singletons map to themselves); `members` only holds groups of size > 1; `member_basis[m]` is the basis of the link that attached non-canonical member `m`.
  - `canonicalise(table, links, *, overrides: IdentityOverrides, reaction_refs: dict[str, int], pathway_refs: dict[str, int], chebi: ChebiIndex, table_pins: list[tuple[str, str, str]]) -> IdentityResult`. `table_pins` are `(table_key, name_key, node_id)` for rows whose node is fixed before resolution (direct id or override target); used to refuse a merge that would join two distinct names of one table.
  - `rewrite_to_canonical(table, reactions: dict[str, dict], subtree_substrates: dict[str, list[str]], canonical_of: dict[str, str]) -> dict[str, dict]` (stage F): returns the new table keyed by canonical ids with list fields unioned (`pathways`, `evidence_sources`), first-non-null scalars **from the canonical member**, and rewrites `reactions[*]["compounds"]` (bare C ids, deduplicated, sorted) and `subtree_substrates` (deduplicated, sorted) in place.

- [ ] **Step 1: Write the failing tests** (`tests/test_mi_canonical.py`):

```python
from tests.metabolite_identity_helpers import cpd, make_chebi, term
from multiomics_kg.utils.metabolite_identity.canonical import (
    IdentityOverrides, MergeLink, canonicalise, collect_merge_links, load_identity_overrides,
    rewrite_to_canonical)

CHEBI = make_chebi(
    CHEBI_17632=term("nitrate", charge=-1, alt=["CHEBI:14654"], conj=["CHEBI:48107"]),
    CHEBI_48107=term("nitric acid", charge=0, kegg=["C00244"], conj=["CHEBI:17632"]),
    CHEBI_15377=term("water", charge=0, kegg=["C00001"], conj=["CHEBI:16234"]),
    CHEBI_16234=term("hydroxide", charge=-1, kegg=["C01328"], conj=["CHEBI:15377"]),
    CHEBI_17234=term("glucose", kegg=["C00293"]),
    CHEBI_28260=term("galactose", kegg=["C01582"]),
)
NO_OV = IdentityOverrides([], set())


def _canon(table, links, **kw):
    kw.setdefault("overrides", NO_OV)
    kw.setdefault("reaction_refs", {})
    kw.setdefault("pathway_refs", {})
    kw.setdefault("table_pins", [])
    return canonicalise(table, links, chebi=CHEBI, **kw)


def test_nitrate_merges_by_conjugate_onto_kegg():
    table = {"kegg.compound:C00244": cpd("Nitrate"),
             "chebi:14654": cpd("nitrate", chebi_id="14654", ev=("transport",))}
    links = collect_merge_links(table, CHEBI)
    assert MergeLink("chebi:14654", "kegg.compound:C00244", "chebi_conjugate") in links
    res = _canon(table, links)
    assert res.canonical_of["chebi:14654"] == "kegg.compound:C00244"
    assert res.basis["kegg.compound:C00244"] == ["chebi_conjugate"]


def test_mnx_same_entity_winner_by_reaction_count():
    table = {"kegg.compound:C00288": cpd("HCO3-", mnxm_id="MNXM60"),
             "kegg.compound:C01353": cpd("Carbonic acid", mnxm_id="MNXM60")}
    res = _canon(table, collect_merge_links(table, CHEBI),
                 reaction_refs={"kegg.compound:C00288": 40, "kegg.compound:C01353": 3})
    assert res.canonical_of["kegg.compound:C01353"] == "kegg.compound:C00288"
    assert res.basis["kegg.compound:C00288"] == ["mnx_same_entity"]


def test_chebi_only_kegg_link_rejected_water_hydroxide():
    table = {"kegg.compound:C00001": cpd("H2O"), "kegg.compound:C01328": cpd("HO-")}
    res = _canon(table, collect_merge_links(table, CHEBI))
    assert res.canonical_of["kegg.compound:C01328"] == "kegg.compound:C01328"
    assert [(l.a, l.b) for l in res.rejected] == [("kegg.compound:C00001", "kegg.compound:C01328")]


def test_manual_merge_overrides_guard_and_keep_separate_blocks_mnx():
    table = {"kegg.compound:C00001": cpd("H2O", mnxm_id="M1"),
             "kegg.compound:C01328": cpd("HO-"),
             "kegg.compound:C00288": cpd("HCO3-", mnxm_id="MNXM60"),
             "kegg.compound:C01353": cpd("Carbonic acid", mnxm_id="MNXM60")}
    ov = IdentityOverrides(
        merges=[MergeLink("kegg.compound:C00001", "kegg.compound:C01328", "manual")],
        keep_separate={frozenset({"kegg.compound:C00288", "kegg.compound:C01353"})})
    res = _canon(table, collect_merge_links(table, CHEBI), overrides=ov)
    assert res.canonical_of["kegg.compound:C01328"] == "kegg.compound:C00001"
    assert "manual" in res.basis["kegg.compound:C00001"]
    assert res.canonical_of["kegg.compound:C01353"] == "kegg.compound:C01353"


def test_stale_identity_override_reported_and_skipped():
    ov = IdentityOverrides(merges=[MergeLink("kegg.compound:C99998", "kegg.compound:C99999", "manual")],
                           keep_separate=set())
    res = _canon({"kegg.compound:C00001": cpd("H2O")}, [], overrides=ov)
    assert res.stale_overrides and res.stale_overrides[0]["kind"] == "merge"


def test_assay_collision_merge_refused():
    table = {"kegg.compound:C00288": cpd("HCO3-", mnxm_id="MNXM60"),
             "kegg.compound:C01353": cpd("Carbonic acid", mnxm_id="MNXM60")}
    pins = [("p|t", "bicarbonate", "kegg.compound:C00288"),
            ("p|t", "carbonic acid", "kegg.compound:C01353")]
    res = _canon(table, collect_merge_links(table, CHEBI), table_pins=pins)
    assert res.canonical_of["kegg.compound:C01353"] == "kegg.compound:C01353"
    assert res.refused_collisions


def test_glucose_galactose_stay_apart():
    table = {"kegg.compound:C00293": cpd("Glucose", inchikey="WQZGKKKJIJFFOK-UHFFFAOYSA-N"),
             "kegg.compound:C01582": cpd("Galactose", inchikey="WQZGKKKJIJFFOK-UHFFFAOYSA-N")}
    res = _canon(table, collect_merge_links(table, CHEBI))
    assert not res.members


def test_chebi_group_prefers_charge0_then_lowest():
    chebi = make_chebi(CHEBI_100=term("x acid", charge=0, conj=["CHEBI:50"]),
                       CHEBI_50=term("x anion", charge=-1, conj=["CHEBI:100"]))
    table = {"chebi:100": cpd("x acid"), "chebi:50": cpd("x anion")}
    res = canonicalise(table, collect_merge_links(table, chebi), overrides=NO_OV,
                       reaction_refs={}, pathway_refs={}, chebi=chebi, table_pins=[])
    assert res.canonical_of["chebi:50"] == "chebi:100"


def test_rewrite_unions_lists_and_dedups_reactions():
    table = {"kegg.compound:C00288": cpd("HCO3-", pathways=["ko00910"], ev=("metabolism",)),
             "kegg.compound:C01353": cpd("Carbonic acid", pathways=["ko00720"], ev=("transport",))}
    canon = {"kegg.compound:C00288": "kegg.compound:C00288",
             "kegg.compound:C01353": "kegg.compound:C00288"}
    rxns = {"R1": {"compounds": ["C00288", "C01353", "C00001"]}}
    subs = {"2.A.1": ["kegg.compound:C01353", "kegg.compound:C00288"]}
    new = rewrite_to_canonical(table, rxns, subs, canon)
    assert list(new) == ["kegg.compound:C00288"]
    assert new["kegg.compound:C00288"]["pathways"] == ["ko00720", "ko00910"]
    assert new["kegg.compound:C00288"]["evidence_sources"] == ["metabolism", "transport"]
    assert rxns["R1"]["compounds"] == ["C00001", "C00288"]
    assert subs["2.A.1"] == ["kegg.compound:C00288"]


def test_load_overrides_contradiction_fails(tmp_path):
    p = tmp_path / "o.yaml"
    p.write_text("merge:\n  - {a: kegg.compound:C1, b: kegg.compound:C2, reason: x}\n"
                 "keep_separate:\n  - {a: kegg.compound:C2, b: kegg.compound:C1, reason: y}\n")
    import pytest
    with pytest.raises(ValueError):
        load_identity_overrides(p)
```

- [ ] **Step 2: Run** `uv run pytest tests/test_mi_canonical.py -q` → FAIL.

- [ ] **Step 3: Implement `canonical.py`.**

```python
"""Stages E-F: merge links, union-find with the KEGG guard, winner rule, id rewrite (spec §3.1)."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from multiomics_kg.utils.chebi_utils import ChebiIndex

MERGE_BASES = ("mnx_same_entity", "chebi_kegg_xref", "chebi_conjugate", "chebi_tautomer", "manual")
_KEGG = "kegg.compound:"


@dataclass(frozen=True, order=True)
class MergeLink:
    a: str
    b: str
    basis: str


def _link(x: str, y: str, basis: str) -> MergeLink:
    a, b = sorted((x, y))
    return MergeLink(a, b, basis)


@dataclass
class IdentityOverrides:
    merges: list[MergeLink]
    keep_separate: set[frozenset[str]]


def load_identity_overrides(path: Path) -> IdentityOverrides:
    if not path.exists():
        return IdentityOverrides([], set())
    doc = yaml.safe_load(path.read_text()) or {}
    merges = [_link(e["a"], e["b"], "manual") for e in (doc.get("merge") or [])]
    keep = {frozenset((e["a"], e["b"])) for e in (doc.get("keep_separate") or [])}
    clash = {frozenset((m.a, m.b)) for m in merges} & keep
    if clash:
        raise ValueError(f"metabolite_identity_overrides.yaml: pairs in both merge and "
                         f"keep_separate: {sorted(sorted(c) for c in clash)}")
    return IdentityOverrides(merges, keep)


def entity_terms(node_id: str, entry: dict, chebi: ChebiIndex) -> set[str]:
    if node_id.startswith(_KEGG):
        return set(chebi.kegg_terms(node_id[len(_KEGG):]))
    if node_id.startswith("chebi:"):
        p = chebi.primary(node_id)
        return {p} if p else set()
    if entry.get("chebi_id"):
        p = chebi.primary(entry["chebi_id"])
        return {p} if p else set()
    return set()


def collect_merge_links(table: dict[str, dict], chebi: ChebiIndex) -> list[MergeLink]:
    links: set[MergeLink] = set()
    by_mnx: dict[str, list[str]] = defaultdict(list)
    terms_of = {n: entity_terms(n, e, chebi) for n, e in table.items()}
    by_term: dict[str, list[str]] = defaultdict(list)
    for n, e in table.items():
        if e.get("mnxm_id"):
            by_mnx[e["mnxm_id"]].append(n)
        for t in terms_of[n]:
            by_term[t].append(n)
    for nodes in by_mnx.values():
        nodes = sorted(nodes)
        links.update(_link(nodes[0], o, "mnx_same_entity") for o in nodes[1:])
    for nodes in by_term.values():
        nodes = sorted(nodes)
        links.update(_link(nodes[0], o, "chebi_kegg_xref") for o in nodes[1:])
    for n, ts in terms_of.items():
        for t in ts:
            tt = chebi.term(t) or {}
            for rel, basis in (("conjugates", "chebi_conjugate"), ("tautomers", "chebi_tautomer")):
                for partner in tt.get(rel, []):
                    p = chebi.primary(partner)
                    for other in by_term.get(p, []):
                        if other != n:
                            links.add(_link(n, other, basis))
    return sorted(links)


@dataclass
class IdentityResult:
    canonical_of: dict[str, str]
    members: dict[str, list[str]]
    basis: dict[str, list[str]]
    member_basis: dict[str, str]
    rejected: list[MergeLink] = field(default_factory=list)
    refused_collisions: list[dict] = field(default_factory=list)
    stale_overrides: list[dict] = field(default_factory=list)


def _choose(members, reaction_refs, pathway_refs, chebi) -> str:
    kegg = [m for m in members if m.startswith(_KEGG)]
    if kegg:
        return sorted(kegg, key=lambda m: (-reaction_refs.get(m, 0), -pathway_refs.get(m, 0), m))[0]
    ch = [m for m in members if m.startswith("chebi:")]
    if ch:
        return sorted(ch, key=lambda m: (chebi.charge(m) != 0, int(m.split(":")[1])))[0]
    return sorted(members)[0]


def canonicalise(table, links, *, overrides, reaction_refs, pathway_refs, chebi, table_pins):
    parent = {n: n for n in table}
    members = {n: {n} for n in table}
    attach: dict[str, str] = {}
    rejected, refused, stale = [], [], []
    pins_by_node: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for tkey, nkey, node in table_pins:
        pins_by_node[node].add((tkey, nkey))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def kegg_ids(root):
        return {m for m in members[root] if m.startswith(_KEGG)}

    def pins(root):
        return set().union(*(pins_by_node.get(m, set()) for m in members[root]))

    for m in overrides.merges:
        if m.a not in table or m.b not in table:
            stale.append({"kind": "merge", "a": m.a, "b": m.b})
    for pair in overrides.keep_separate:
        if not all(x in table for x in pair):
            stale.append({"kind": "keep_separate", "pair": sorted(pair)})
    manual = [m for m in overrides.merges if m.a in table and m.b in table]

    ordered = manual + sorted((l for l in links if l.basis == "mnx_same_entity")) + \
        sorted((l for l in links if l.basis != "mnx_same_entity"))
    for link in ordered:
        ra, rb = find(link.a), find(link.b)
        if ra == rb:
            continue
        ga, gb = members[ra], members[rb]
        if any(frozenset((x, y)) in overrides.keep_separate for x in ga for y in gb):
            continue
        ka, kb = kegg_ids(ra), kegg_ids(rb)
        if link.basis not in ("mnx_same_entity", "manual") and ka and kb and ka != kb:
            rejected.append(link)
            continue
        pa, pb = pins(ra), pins(rb)
        clash = [(t, n1, n2) for t, n1 in pa for t2, n2 in pb if t == t2 and n1 != n2]
        if clash:
            refused.append({"link": [link.a, link.b, link.basis],
                            "tables": sorted({c[0] for c in clash})})
            continue
        parent[rb] = ra
        members[ra] |= members.pop(rb)
        # Each member is attributed the basis of the first link that joined it to anything.
        for m in ga | gb:
            attach.setdefault(m, link.basis)

    canonical_of, groups, basis, member_basis = {}, {}, {}, {}
    for root, ms in members.items():
        canon = _choose(ms, reaction_refs, pathway_refs, chebi)
        for m in ms:
            canonical_of[m] = canon
        if len(ms) > 1:
            groups[canon] = sorted(ms)
            for m in ms:
                if m != canon:
                    member_basis[m] = attach[m]
            basis[canon] = sorted({member_basis[m] for m in ms if m != canon})
    return IdentityResult(canonical_of, groups, basis, member_basis, sorted(rejected), refused, stale)


def _union(a, b):
    return sorted(set(a or []) | set(b or []))


def rewrite_to_canonical(table, reactions, subtree_substrates, canonical_of):
    new: dict[str, dict] = {}
    for node in sorted(table, key=lambda n: (canonical_of[n] != n, n)):
        canon = canonical_of[node]
        src = table[node]
        if canon not in new:
            new[canon] = dict(table[canon]) if canon in table else dict(src)
            new[canon]["pathways"] = sorted(set(new[canon].get("pathways") or []))
            new[canon]["evidence_sources"] = sorted(set(new[canon].get("evidence_sources") or []))
        tgt = new[canon]
        tgt["pathways"] = _union(tgt.get("pathways"), src.get("pathways"))
        tgt["evidence_sources"] = _union(tgt.get("evidence_sources"), src.get("evidence_sources"))
        for k, v in src.items():
            if tgt.get(k) is None and v is not None:
                tgt[k] = v
    bare = {n[len(_KEGG):]: c[len(_KEGG):] for n, c in canonical_of.items()
            if n.startswith(_KEGG) and c.startswith(_KEGG)}
    for r in reactions.values():
        r["compounds"] = sorted({bare.get(c, c) for c in r.get("compounds", [])})
    for tc, ids in subtree_substrates.items():
        subtree_substrates[tc] = sorted({canonical_of.get(i, i) for i in ids})
    return new
```

Also create `config/metabolite_identity_overrides.yaml`:

```yaml
# Curated metabolite identity decisions (spec §3.1, §3.6a). Read by prepare-data step 6.
# merge:         two node ids that ARE one metabolite but no curated source links them.
#                Bypasses the KEGG↔KEGG guard; recorded as identity_basis [manual].
# keep_separate: two node ids that must never merge, whatever the basis (MNX included).
# Each entry: {a: <node id>, b: <node id>, reason: "<one line>"}. Neither kind adds a node.
# Candidates to review: metabolite_identity_report.json → heuristic_name_candidates.
merge: []
keep_separate: []
```

- [ ] **Step 4: Run** `uv run pytest tests/test_mi_canonical.py -q` → PASS.

- [ ] **Step 5: Commit.**

```bash
git add multiomics_kg/utils/metabolite_identity/canonical.py config/metabolite_identity_overrides.yaml tests/test_mi_canonical.py
git commit -m "feat(metabolite): merge links, guarded union-find, canonical winner, id rewrite"
```

### Task 4: Structure and cross-reference fields (`structure.py`, `crossrefs.py`, stages G–H)

**Files:**
- Create: `multiomics_kg/utils/metabolite_identity/structure.py`, `multiomics_kg/utils/metabolite_identity/crossrefs.py`, `tests/test_mi_structure.py`, `tests/test_mi_crossrefs.py`

**Interfaces:**
- Consumes: `names.is_generic_formula`, `names.is_polymer_formula`, `names.parse_elements`, `names.normalize_hmdb` (Task 2); `canonical.IdentityResult`, `canonical.entity_terms` (Task 3); `ChebiIndex` (Task 1).
- Produces:
  - `STRUCTURE_KINDS = ("defined", "generic", "polymer")`.
  - `classify_structure(entry: dict, chebi_term_id: str | None, chebi: ChebiIndex) -> str` (plan "Spec gaps" rule 6).
  - `apply_structure(table: dict[str, dict], chebi: ChebiIndex, chebi_of: dict[str, str | None]) -> dict` (stage G): sets `structure_kind` and `elements` on every entry; fills a missing formula from the ChEBI term's generalized formula; scrubs `inchikey` / `smiles` / concrete formula / `mass` per spec §3.3 (both conditions); drops `mass` on generic / polymer; keeps `elements`; returns `{"scrubbed": [...], "filled": [...], "counts": {...}}`. Raises `ValueError` listing any InChIKey still on two nodes.
  - `pick_chebi_id(node_id: str, entry: dict, chebi: ChebiIndex) -> str | None` → bare primary number.
  - `build_identity_fields(table, identity: IdentityResult, rekeys: dict[str, str], chebi: ChebiIndex, old_entries: dict[str, dict]) -> dict[str, str]` (stage H): sets `kegg_compound_id`, `kegg_compound_ids`, `chebi_id`, `chebi_ids`, `mnxm_id`, `mnxm_ids` (sparse), `hmdb_id`, `all_identifiers`, `identity_basis` (sparse), `name`; returns `identifier_index: {curie: canonical node id}` and raises `ValueError` if any identifier would sit on two nodes. `old_entries` holds the pre-merge entries of every member (for their ids / MNX / HMDB / InChIKey); `rekeys` maps the pre-re-key node id (`chebi:…`) to its KEGG node id.

- [ ] **Step 1: Write the failing tests.** `tests/test_mi_structure.py`:

```python
import pytest

from tests.metabolite_identity_helpers import cpd, make_chebi, term
from multiomics_kg.utils.metabolite_identity.structure import apply_structure, classify_structure

CHEBI = make_chebi(
    CHEBI_17522=term("alditol", gformula="C2H6O2(CH2O)n"),
    CHEBI_17754=term("glycerol", charge=0, kegg=["C00116"], inchikey="PEDCQBHIVMGVHV-UHFFFAOYSA-N"),
    CHEBI_37684=term("mannose", gformula="C6H12O6"))


def test_classify():
    assert classify_structure(cpd("Rubidium cation", formula="Rb"), None, CHEBI) == "defined"
    assert classify_structure(cpd("x", formula="C5H8O2R"), None, CHEBI) == "generic"
    assert classify_structure(cpd("PHB", formula="(C4H6O2)n"), None, CHEBI) == "polymer"
    assert classify_structure(cpd("(1,4-alpha-D-glucosyl)n"), None, CHEBI) == "polymer"
    assert classify_structure(cpd("alditol"), "CHEBI:17522", CHEBI) == "polymer"   # formula fill decides
    assert classify_structure(cpd("thing"), None, CHEBI) == "defined"


def test_b2_scrub_by_class_and_by_shared_inchikey():
    ik = "PEDCQBHIVMGVHV-UHFFFAOYSA-N"
    table = {"kegg.compound:C00116": cpd("Glycerol", formula="C3H8O3", mass=92.09, inchikey=ik),
             "chebi:17522": cpd("alditol", formula="C3H8O3", mass=92.09, inchikey=ik),
             "kegg.compound:C00734": cpd("Chitosan", formula="C6H13NO5", mass=179.2,
                                         inchikey="MSWZFWKMSRAUBD-IVMDWMLBSA-N"),
             "kegg.compound:C00329": cpd("Glucosamine", formula="C6H13NO5", mass=179.2,
                                         inchikey="MSWZFWKMSRAUBD-IVMDWMLBSA-N")}
    table["kegg.compound:C00734"]["name"] = "Chitosan"
    rep = apply_structure(table, CHEBI, {"chebi:17522": "CHEBI:17522"})
    alditol = table["chebi:17522"]
    assert alditol["structure_kind"] in ("generic", "polymer")
    assert alditol.get("inchikey") is None and alditol.get("mass") is None
    assert alditol["elements"] == ["C", "H", "O"]
    assert table["kegg.compound:C00116"]["inchikey"] == ik            # defined side kept
    assert table["kegg.compound:C00734"].get("inchikey") is None       # defined vs defined: KEGG polymer by name
    assert set(rep["scrubbed"]) >= {"chebi:17522", "kegg.compound:C00734"}


def test_formula_fill_and_no_mass_on_generic():
    table = {"chebi:37684": cpd("mannose"), "kegg.compound:C99999": cpd("R-thing", formula="C2H4R", mass=12.0)}
    apply_structure(table, CHEBI, {"chebi:37684": "CHEBI:37684"})
    assert table["chebi:37684"]["formula"] == "C6H12O6" and table["chebi:37684"]["elements"] == ["C", "H", "O"]
    assert table["kegg.compound:C99999"].get("mass") is None


def test_leftover_inchikey_collision_fails():
    ik = "AAAAAAAAAAAAAA-BBBBBBBBBB-N"
    table = {"kegg.compound:C1": cpd("a", formula="CH4", inchikey=ik),
             "kegg.compound:C2": cpd("b", formula="CH4", inchikey=ik)}
    with pytest.raises(ValueError, match="InChIKey"):
        apply_structure(table, CHEBI, {})
```

  The two-defined-nodes-share-an-InChIKey case (chitosan / glucosamine) is resolved by the name: a KEGG entry whose name is a known polymer (`Chitosan`, `Chitin`, `Poly…`, `…ose)n`) or whose KEGG formula carries `)n` is `polymer`. The executor adds a `_POLYMER_NAMES` regex `^(poly|chitosan|chitin|cellulose|starch|glycogen|amylose|peptidoglycan)` (case-insensitive) to `classify_structure` and checks that the B2 inventory (spec §1.1: alditol / glycerol, (1,4-α-D-glucosyl)n / glucose, chitosan / glucosamine, PHB / 3-hydroxybutanoate, chitin / GlcNAc) all resolve with the generic side scrubbed. A pair of two `defined` nodes left sharing a key fails the step (last test).

  `tests/test_mi_crossrefs.py`:

```python
from tests.metabolite_identity_helpers import cpd, make_chebi, term
from multiomics_kg.utils.metabolite_identity.canonical import IdentityResult
from multiomics_kg.utils.metabolite_identity.crossrefs import build_identity_fields, pick_chebi_id

CHEBI = make_chebi(
    CHEBI_17632=term("nitrate", charge=-1, alt=["CHEBI:14654", "CHEBI:44487"], kegg=["C00244"],
                     inchikey="NHNBFGGVMKEFGY-UHFFFAOYSA-N"),
    CHEBI_48107=term("nitric acid", charge=0, kegg=["C00244"], inchikey="GRYLNZFGIOXLOG-UHFFFAOYSA-N"),
    CHEBI_25545=term("nitrates"),                                   # class: no InChIKey
    CHEBI_31000=term("succinic acid-d4", inchikey="KDYFGRWQOYBRFD-NMQOAUCRSA-N"),
    CHEBI_15741=term("succinic acid", charge=0, kegg=["C00042"], inchikey="KDYFGRWQOYBRFD-UHFFFAOYSA-N"),
    CHEBI_26806=term("succinate(2-)", charge=-2, inchikey="KDYFGRWQOYBRFD-UHFFFAOYSA-L"),
    CHEBI_30031=term("succinate(1-)", charge=-1, inchikey="KDYFGRWQOYBRFD-UHFFFAOYSA-M"))


def test_chebi_id_prefers_kegg_xref_charge0():
    assert pick_chebi_id("kegg.compound:C00244", cpd("Nitrate", chebi_id="14654"), CHEBI) == "48107"


def test_isotopologue_never_chebi_id():
    assert pick_chebi_id("kegg.compound:C00042", cpd("Succinate", chebi_id="31000"), CHEBI) == "15741"


def test_mnx_fallback_ranks_aliases_and_skips_isotopologues_and_classes():
    e = cpd("Succinate-ish")
    e["mnx_chebi_aliases"] = ["31000", "25545", "30031", "26806"]
    # C09999 has no ChEBI cross-reference, so the MNX aliases decide:
    # isotopologue 31000 and class 25545 are excluded; no charge-0 member → lowest number.
    assert pick_chebi_id("kegg.compound:C09999", e, CHEBI) == "26806"


def test_identity_fields_nitrate():
    table = {"kegg.compound:C00244": cpd("Nitrate", chebi_id="14654", mnxm_id="MNXM732399",
                                         hmdb_id="HMDB02878", inchikey="GRYLNZFGIOXLOG-UHFFFAOYSA-N")}
    old = {"kegg.compound:C00244": dict(table["kegg.compound:C00244"]),
           "chebi:14654": cpd("nitrate", chebi_id="14654", mnxm_id="MNXM732398")}
    ident = IdentityResult(
        canonical_of={"kegg.compound:C00244": "kegg.compound:C00244", "chebi:14654": "kegg.compound:C00244"},
        members={"kegg.compound:C00244": ["chebi:14654", "kegg.compound:C00244"]},
        basis={"kegg.compound:C00244": ["chebi_conjugate"]},
        member_basis={"chebi:14654": "chebi_conjugate"})
    index = build_identity_fields(table, ident, {}, CHEBI, old)
    e = table["kegg.compound:C00244"]
    assert e["chebi_id"] == "48107" and e["chebi_ids"][0] == "48107"
    assert set(e["chebi_ids"]) == {"48107", "17632"}
    assert e["mnxm_ids"] == ["MNXM732399", "MNXM732398"]
    assert e["hmdb_id"] == "HMDB0002878" and e["identity_basis"] == ["chebi_conjugate"]
    assert {"kegg.compound:C00244", "chebi:14654", "chebi:17632", "chebi:44487", "chebi:48107",
            "mnx:MNXM732398", "mnx:MNXM732399", "hmdb:HMDB0002878",
            "inchikey:GRYLNZFGIOXLOG-UHFFFAOYSA-N"} <= set(e["all_identifiers"])
    assert all(":" in x for x in e["all_identifiers"])
    assert index["chebi:14654"] == "kegg.compound:C00244"


def test_rekeyed_node_has_basis_and_old_id():
    table = {"kegg.compound:C00244": cpd("Nitrate")}
    ident = IdentityResult({"kegg.compound:C00244": "kegg.compound:C00244"}, {}, {}, {})
    build_identity_fields(table, ident, {"chebi:17632": "kegg.compound:C00244"}, CHEBI,
                          {"kegg.compound:C00244": dict(table["kegg.compound:C00244"])})
    e = table["kegg.compound:C00244"]
    assert e["identity_basis"] == ["chebi_kegg_xref"] and "chebi:17632" in e["all_identifiers"]


def test_chebi_keyed_node_chebi_id_equals_local_part():
    table = {"chebi:17632": cpd("nitrate", chebi_id="44487")}
    ident = IdentityResult({"chebi:17632": "chebi:17632"}, {}, {}, {})
    build_identity_fields(table, ident, {}, CHEBI, {"chebi:17632": dict(table["chebi:17632"])})
    assert table["chebi:17632"]["chebi_id"] == "17632"
    assert table["chebi:17632"].get("kegg_compound_ids") is None    # dropped as empty at stage P
```

- [ ] **Step 2: Run** `uv run pytest tests/test_mi_structure.py tests/test_mi_crossrefs.py -q` → FAIL.

- [ ] **Step 3: Implement.** `structure.py`:

```python
"""Stage G: structure_kind, B2 scrub, generalized-formula fill (spec §3.3)."""
from __future__ import annotations

import re
from collections import defaultdict

from multiomics_kg.utils.chebi_utils import ChebiIndex
from multiomics_kg.utils.metabolite_identity.names import (
    is_generic_formula, is_polymer_formula, parse_elements)

STRUCTURE_KINDS = ("defined", "generic", "polymer")
_POLYMER_NAMES = re.compile(
    r"^(poly|chitosan|chitin|cellulose|starch|glycogen|amylose|peptidoglycan)|\)n$", re.I)


def classify_structure(entry: dict, chebi_term_id: str | None, chebi: ChebiIndex) -> str:
    formula = entry.get("formula")
    if formula is None and chebi_term_id:
        formula = (chebi.term(chebi_term_id) or {}).get("generalized_formula")
    if is_polymer_formula(formula) or _POLYMER_NAMES.search(entry.get("name") or ""):
        return "polymer"
    if is_generic_formula(formula):
        return "generic"
    if chebi_term_id and chebi.is_class(chebi_term_id):
        return "generic"
    return "defined"


def _is_placeholder(formula: str | None) -> bool:
    return is_generic_formula(formula) or is_polymer_formula(formula)


def apply_structure(table, chebi, chebi_of):
    filled, scrubbed = [], []
    for node, e in table.items():
        term_id = chebi_of.get(node)
        gformula = (chebi.term(term_id) or {}).get("generalized_formula") if term_id else None
        kind = classify_structure(e, term_id, chebi)
        e["structure_kind"] = kind
        if kind != "defined":
            concrete = bool(e.get("formula")) and not _is_placeholder(e["formula"])
            if concrete or e.get("inchikey") or e.get("smiles"):
                # B2: drop the representative structure, keep its element set (spec §3.3).
                e["_elements_keep"] = parse_elements(e.get("formula"))
                e["inchikey"] = e["smiles"] = None
                if concrete:
                    e["formula"] = gformula if _is_placeholder(gformula) else None
                scrubbed.append(node)
            e["mass"] = None
        if e.get("formula") is None and gformula:
            e["formula"] = gformula
            filled.append(node)
        if _is_placeholder(e.get("formula")):
            e["mass"] = None
    # Condition 2: a non-defined node sharing an InChIKey with a defined node is already scrubbed;
    # anything left sharing a key is a build error.
    by_key = defaultdict(list)
    for node, e in table.items():
        if e.get("inchikey"):
            by_key[e["inchikey"]].append(node)
    dup = {k: v for k, v in by_key.items() if len(v) > 1}
    if dup:
        raise ValueError(f"InChIKey on more than one Metabolite after the B2 scrub: {dup}")
    counts = defaultdict(int)
    for e in table.values():
        keep = e.pop("_elements_keep", None)
        e["elements"] = sorted(set(parse_elements(e.get("formula"))) | set(keep or []))
        counts[e["structure_kind"]] += 1
    return {"scrubbed": sorted(scrubbed), "filled": sorted(filled), "counts": dict(counts)}
```

  `crossrefs.py`:

```python
"""Stage H: chebi_id pick, id arrays, all_identifiers, identity_basis, name (spec §3.2)."""
from __future__ import annotations

import re

from multiomics_kg.utils.chebi_utils import ChebiIndex
from multiomics_kg.utils.metabolite_identity.canonical import IdentityResult
from multiomics_kg.utils.metabolite_identity.names import normalize_hmdb

_KEGG = "kegg.compound:"


def _bare(cid: str) -> str:
    return cid.split(":", 1)[1]


def _entity_rank(chebi: ChebiIndex):
    return lambda c: (chebi.charge(c) != 0, int(_bare(c)))


# Isotopologue names (spec §1.1 B1b): "succinic acid-d4", "UREA C 14", "calcium-43(2+)", "[13C]…", "(15N)…".
# Deliberately narrow: must NOT match "1H-indole", "Vitamin B12", "(2E)-hexenal", "(3S)-…".
_ISOTOPE = re.compile(
    r"-d\d+\b|isotop|deuter|[\[(]\d{1,3}(?:C|N|O|H|S|P)[\])]|\b(?:C|N|O|H) 1[0-9]\b|-\d{1,3}\(\d", re.I)


def _usable(p: str | None, node_id: str, chebi: ChebiIndex) -> bool:
    """An entity term for this node: not a class, not an isotopologue, not another KEGG compound's term."""
    if p is None or chebi.is_class(p):
        return False
    t = chebi.term(p) or {}
    if _ISOTOPE.search(t.get("name") or ""):
        return False
    xref = t.get("kegg", [])
    return not xref or (node_id.startswith(_KEGG) and node_id[len(_KEGG):] in xref)


def pick_chebi_id(node_id, entry, chebi):
    if node_id.startswith("chebi:"):
        p = chebi.primary(node_id)
        return _bare(p) if p else _bare(node_id)
    cands = []
    if node_id.startswith(_KEGG):
        cands = [c for c in chebi.kegg_terms(node_id[len(_KEGG):]) if not chebi.is_class(c)]
    if not cands:
        aliases = list(entry.get("mnx_chebi_aliases") or []) + ([entry["chebi_id"]] if entry.get("chebi_id") else [])
        cands = sorted({p for p in (chebi.primary(a) for a in aliases) if _usable(p, node_id, chebi)})
    if not cands:
        return None
    return _bare(sorted(cands, key=_entity_rank(chebi))[0])


def build_identity_fields(table, identity: IdentityResult, rekeys, chebi, old_entries):
    rekeyed_into: dict[str, list[str]] = {}
    for old, new in rekeys.items():
        rekeyed_into.setdefault(new, []).append(old)
    index: dict[str, str] = {}
    for node in sorted(table):
        e = table[node]
        group = identity.members.get(node, [node])
        olds = sorted(set(group) | set(rekeyed_into.get(node, [])))
        kegg = [node[len(_KEGG):]] if node.startswith(_KEGG) else []
        kegg += sorted({m[len(_KEGG):] for m in olds if m.startswith(_KEGG)} - set(kegg))
        e["kegg_compound_id"] = kegg[0] if kegg else None
        e["kegg_compound_ids"] = kegg or None
        chebi_id = pick_chebi_id(node, e, chebi)
        e["chebi_id"] = chebi_id
        chebi_terms = set()
        for m in olds:
            src = old_entries.get(m, {})
            if m.startswith("chebi:"):
                p = chebi.primary(m)
                if p:
                    chebi_terms.add(p)
            if m.startswith(_KEGG):
                chebi_terms |= {c for c in chebi.kegg_terms(m[len(_KEGG):]) if not chebi.is_class(c)}
            for a in list(src.get("mnx_chebi_aliases") or []) + ([src["chebi_id"]] if src.get("chebi_id") else []):
                p = chebi.primary(a)
                if _usable(p, node, chebi):
                    chebi_terms.add(p)
        ids = [chebi_id] if chebi_id else []
        ids += sorted({_bare(c) for c in chebi_terms} - set(ids), key=int)
        e["chebi_ids"] = ids or None
        mnx = [e["mnxm_id"]] if e.get("mnxm_id") else []
        mnx += sorted({old_entries.get(m, {}).get("mnxm_id") for m in olds} - {None} - set(mnx))
        e["mnxm_id"] = mnx[0] if mnx else None
        e["mnxm_ids"] = mnx if len(mnx) > 1 else None
        e["hmdb_id"] = normalize_hmdb(e.get("hmdb_id"))
        ident = {node} | set(olds)
        ident |= {f"{_KEGG}{k}" for k in kegg}
        for c in chebi_terms | ({f"CHEBI:{chebi_id}"} if chebi_id else set()):
            ident.add(f"chebi:{_bare(c)}")
            ident |= {f"chebi:{_bare(a)}" for a in chebi.secondary_ids(c)}
        ident |= {f"mnx:{m}" for m in mnx}
        hmdbs = {normalize_hmdb(old_entries.get(m, {}).get("hmdb_id")) for m in olds} | {e["hmdb_id"]}
        ident |= {f"hmdb:{h}" for h in hmdbs if h}
        if e.get("inchikey"):
            ident.add(f"inchikey:{e['inchikey'].upper()}")
        e["all_identifiers"] = sorted(ident)
        basis = set(identity.basis.get(node, []))
        if rekeyed_into.get(node):
            basis.add("chebi_kegg_xref")
        e["identity_basis"] = sorted(basis) or None
        for x in e["all_identifiers"]:
            if index.setdefault(x, node) != node:
                raise ValueError(f"identifier {x} on two Metabolite nodes: {index[x]} and {node}")
    return index
```

  `mnx_chebi_aliases` is a transient field stage D adds (every MNX `chebi` alias of the entry's MNXM id, Task 11); stage P removes it before writing. The `_ISOTOPE` pattern must reject the three B1b examples (spec §1.1: "succinic acid-d4", "UREA C 14", "calcium-43(2+)"); add them to `test_mnx_fallback_ranks_aliases_and_skips_isotopologues_and_classes` as extra excluded aliases, and add a negative test that `"1H-indole-3-acetate"`, `"Vitamin B12"`, `"(2E)-hexenal"` and `"(3S)-citramalate"` are **not** isotopologues.

  Name rule (spec §3.2): KEGG-keyed nodes keep the KEGG name already in the entry; on `chebi:` nodes set `e["name"] = chebi.term(node)["name"]` inside the loop when the node id starts with `chebi:`.

- [ ] **Step 4: Run** `uv run pytest tests/test_mi_structure.py tests/test_mi_crossrefs.py -q` → PASS.

- [ ] **Step 5: Commit.**

```bash
git add multiomics_kg/utils/metabolite_identity/{structure,crossrefs}.py tests/test_mi_structure.py tests/test_mi_crossrefs.py
git commit -m "feat(metabolite): structure_kind + B2 scrub, canonical cross-reference fields"
```

### Task 5: Stereo links, families and the anomer test (`stereo.py`, stages I and L)

**Files:**
- Create: `multiomics_kg/utils/metabolite_identity/stereo.py`, `tests/test_mi_stereo.py`

**Interfaces:**
- Consumes: `ChebiIndex`; `names.strip_anomer_tokens`; entries carrying `structure_kind` (Task 4).
- Produces:
  - `stereo_links(table: dict[str, dict], terms_of: dict[str, set[str]], chebi: ChebiIndex) -> tuple[set[tuple[str, str]], set[tuple[str, str]]]` → `(stereo_pairs, enantiomer_pairs)`. `stereo_pairs` holds `(specific, unspecified)` where the specific node's ChEBI term `is_a` the other's; `enantiomer_pairs` holds `(lower_id, higher_id)` from `is_enantiomer_of`. Both require both nodes to have an InChIKey with an equal first block (the scrubbed classes have none, so a class can't join). Pairs between members of one merge group can't occur: stage F already collapsed them.
  - `stereo_families(nodes: set[str], stereo_pairs, enantiomer_pairs) -> dict[str, frozenset[str]]`: connected components; every node in `nodes` maps to its family (a singleton when unlinked).
  - `unspecified_members(stereo_pairs) -> set[str]` = the set of pair targets.
  - `is_anomer_pair(a: str, b: str, table, terms_of, chebi) -> bool` (plan "Spec gaps" rule 5).
  - `finalize_stereo(table, stereo_pairs, enantiomer_pairs, survivors: set[str]) -> dict[str, int]` (stage L): writes `stereo_form_of` (sorted full ids) and `enantiomer_of` (on the lower id only) onto surviving entries, skipping any pair with a non-survivor; returns `{"stereo_form_of": n, "enantiomer_of": m}`.

- [ ] **Step 1: Write the failing tests** (`tests/test_mi_stereo.py`):

```python
from tests.metabolite_identity_helpers import cpd, make_chebi, term
from multiomics_kg.utils.metabolite_identity.stereo import (
    finalize_stereo, is_anomer_pair, stereo_families, stereo_links, unspecified_members)

CHEBI = make_chebi(
    CHEBI_16449=term("alanine"),
    CHEBI_16977=term("L-alanine", is_a=["CHEBI:16449"], enan=["CHEBI:15570"]),
    CHEBI_15570=term("D-alanine", is_a=["CHEBI:16449"], enan=["CHEBI:16977"]),
    CHEBI_16414=term("L-valine", enan=["CHEBI:27266"]),
    CHEBI_27266=term("D-valine", enan=["CHEBI:16414"]),
    CHEBI_16646=term("carbohydrate"),
    CHEBI_4167=term("D-glucopyranose", is_a=["CHEBI:16646"]),
    CHEBI_15903=term("beta-D-glucose", is_a=["CHEBI:4167"]),
    CHEBI_17522=term("alditol"),
    CHEBI_17754=term("glycerol", is_a=["CHEBI:17522"]),
    CHEBI_100=term("citronellol"),
    CHEBI_101=term("beta-citronellol", enan=["CHEBI:102"]),
    CHEBI_102=term("(-)-citronellol", enan=["CHEBI:101"]))

def _d(name, ik):
    e = cpd(name, inchikey=ik)
    e["structure_kind"] = "defined"
    return e

T = {"kegg.compound:C01401": _d("Alanine", "QNAYBMKLOCPYGJ-UHFFFAOYSA-N"),
     "kegg.compound:C00041": _d("L-Alanine", "QNAYBMKLOCPYGJ-REOHCLBHSA-N"),
     "kegg.compound:C00133": _d("D-Alanine", "QNAYBMKLOCPYGJ-UWTATZPHSA-N"),
     "kegg.compound:C00183": _d("L-Valine", "KZSNJWFQEVHDMF-BYPYZUCNSA-N"),
     "kegg.compound:C06417": _d("D-Valine", "KZSNJWFQEVHDMF-SCSAIBSYSA-N"),
     "kegg.compound:C00031": _d("D-Glucose", "WQZGKKKJIJFFOK-GASJEMHNSA-N"),
     "kegg.compound:C00221": _d("beta-D-Glucose", "WQZGKKKJIJFFOK-VFUOTHLCSA-N"),
     "chebi:17522": dict(cpd("alditol"), structure_kind="generic"),
     "kegg.compound:C00116": _d("Glycerol", "PEDCQBHIVMGVHV-UHFFFAOYSA-N"),
     "chebi:101": _d("beta-Citronellol", "QMVPMAAFGQKVCJ-JTQLQIEISA-N"),
     "chebi:102": _d("(-)-Citronellol", "QMVPMAAFGQKVCJ-SNVBAGLBSA-N")}
TERMS = {"kegg.compound:C01401": {"CHEBI:16449"}, "kegg.compound:C00041": {"CHEBI:16977"},
         "kegg.compound:C00133": {"CHEBI:15570"}, "kegg.compound:C00183": {"CHEBI:16414"},
         "kegg.compound:C06417": {"CHEBI:27266"}, "kegg.compound:C00031": {"CHEBI:4167"},
         "kegg.compound:C00221": {"CHEBI:15903"}, "chebi:17522": {"CHEBI:17522"},
         "kegg.compound:C00116": {"CHEBI:17754"}, "chebi:101": {"CHEBI:101"}, "chebi:102": {"CHEBI:102"}}


def test_links_need_isa_or_enantiomer_and_same_first_block():
    s, e = stereo_links(T, TERMS, CHEBI)
    assert ("kegg.compound:C00041", "kegg.compound:C01401") in s
    assert ("kegg.compound:C00133", "kegg.compound:C01401") in s
    assert ("kegg.compound:C00041", "kegg.compound:C00133") in e
    assert not any("chebi:17522" in p for p in s)          # class scrubbed: no InChIKey


def test_enantiomers_without_unspecified_node_stored_once():
    _, e = stereo_links(T, TERMS, CHEBI)
    assert ("kegg.compound:C00183", "kegg.compound:C06417") in e
    assert ("kegg.compound:C06417", "kegg.compound:C00183") not in e


def test_families_and_unspecified():
    s, e = stereo_links(T, TERMS, CHEBI)
    fam = stereo_families(set(T), s, e)
    assert fam["kegg.compound:C00041"] == fam["kegg.compound:C00133"] == fam["kegg.compound:C01401"]
    assert "kegg.compound:C01401" in unspecified_members(s)


def test_anomer_guard():
    assert is_anomer_pair("kegg.compound:C00031", "kegg.compound:C00221", T, TERMS, CHEBI)
    assert not is_anomer_pair("chebi:101", "chebi:102", T, TERMS, CHEBI)   # beta-Citronellol: not a sugar
    assert not is_anomer_pair("kegg.compound:C00041", "kegg.compound:C00133", T, TERMS, CHEBI)


def test_finalize_drops_pairs_with_non_survivors():
    s, e = stereo_links(T, TERMS, CHEBI)
    table = {k: dict(v) for k, v in T.items()}
    survivors = set(table) - {"kegg.compound:C00133"}
    counts = finalize_stereo(table, s, e, survivors)
    assert "enantiomer_of" not in table["kegg.compound:C00041"]
    assert table["kegg.compound:C00041"]["stereo_form_of"] == ["kegg.compound:C01401"]
    assert table["kegg.compound:C00183"]["enantiomer_of"] == ["kegg.compound:C06417"]
    assert counts["enantiomer_of"] >= 1
```

- [ ] **Step 2: Run** `uv run pytest tests/test_mi_stereo.py -q` → FAIL.

- [ ] **Step 3: Implement `stereo.py`.**

```python
"""Stages I and L: stereo families (linked, not merged) and the anomer guard (spec §3.5, §3.6)."""
from __future__ import annotations

from collections import defaultdict

from multiomics_kg.utils.chebi_utils import ChebiIndex
from multiomics_kg.utils.metabolite_identity.names import strip_anomer_tokens


def _block(entry: dict) -> str | None:
    ik = entry.get("inchikey")
    return ik.split("-")[0] if ik else None


def stereo_links(table, terms_of, chebi: ChebiIndex):
    by_term: dict[str, list[str]] = defaultdict(list)
    for n, ts in terms_of.items():
        if n in table:
            for t in ts:
                by_term[t].append(n)
    stereo, enan = set(), set()
    for n, ts in terms_of.items():
        if n not in table or not _block(table[n]):
            continue
        for t in ts:
            tt = chebi.term(t) or {}
            for rel, out in (("is_a", stereo), ("enantiomers", enan)):
                for partner in tt.get(rel, []):
                    for o in by_term.get(chebi.primary(partner) or "", []):
                        if o == n or _block(table[o]) != _block(table[n]):
                            continue
                        out.add((n, o) if rel == "is_a" else tuple(sorted((n, o))))
    return stereo, enan


def stereo_families(nodes, stereo_pairs, enantiomer_pairs):
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in list(stereo_pairs) + list(enantiomer_pairs):
        if a in parent and b in parent:
            parent[find(a)] = find(b)
    comps: dict[str, set[str]] = defaultdict(set)
    for n in nodes:
        comps[find(n)].add(n)
    return {n: frozenset(comps[find(n)]) for n in nodes}


def unspecified_members(stereo_pairs):
    return {b for _, b in stereo_pairs}


def is_anomer_pair(a, b, table, terms_of, chebi):
    names_a = {table[a].get("name") or ""} | {(chebi.term(t) or {}).get("name") or "" for t in terms_of.get(a, ())}
    names_b = {table[b].get("name") or ""} | {(chebi.term(t) or {}).get("name") or "" for t in terms_of.get(b, ())}
    same = any(x and y and x.casefold() != y.casefold() and strip_anomer_tokens(x) == strip_anomer_tokens(y)
               for x in names_a for y in names_b)
    sugar = any(chebi.is_carbohydrate(t) for t in terms_of.get(a, set()) | terms_of.get(b, set()))
    return same and sugar


def finalize_stereo(table, stereo_pairs, enantiomer_pairs, survivors):
    for e in table.values():
        e.pop("stereo_form_of", None)
        e.pop("enantiomer_of", None)
    n_s = n_e = 0
    for a, b in sorted(stereo_pairs):
        if a in survivors and b in survivors:
            table[a].setdefault("stereo_form_of", []).append(b)
            n_s += 1
    for a, b in sorted(enantiomer_pairs):
        if a in survivors and b in survivors:
            table[a].setdefault("enantiomer_of", []).append(b)
            n_e += 1
    for e in table.values():
        for k in ("stereo_form_of", "enantiomer_of"):
            if k in e:
                e[k] = sorted(set(e[k]))
    return {"stereo_form_of": n_s, "enantiomer_of": n_e}
```

- [ ] **Step 4: Run** `uv run pytest tests/test_mi_stereo.py -q` → PASS.

- [ ] **Step 5: Commit.**

```bash
git add multiomics_kg/utils/metabolite_identity/stereo.py tests/test_mi_stereo.py
git commit -m "feat(metabolite): stereo families, enantiomer pairs, anomer guard"
```

### Task 6: Role partition and the `ChemicalRole` hierarchy (`roles.py`, stage C)

**Files:**
- Create: `multiomics_kg/utils/metabolite_identity/roles.py`, `tests/test_mi_roles.py`

**Interfaces:**
- Consumes: `ChebiIndex`, `ROLE_ROOT`, `ROLE_BRANCHES` (Task 1).
- Produces:
  - `role_node_id(primary: str) -> str` → `chebi.role:<n>` (plain f-string; never `normalize_curie`).
  - `partition_role_strings(strings: set[str], chebi: ChebiIndex) -> tuple[set[str], dict[str, str]]` → `(molecule_strings, {substrate_string: role_node_id})`; a string is a role when its **primary** ChEBI term has `CHEBI:50906` among its ancestors (spec §3.4).
  - `build_role_hierarchy(role_ids: set[str], chebi: ChebiIndex) -> dict[str, dict]` → the `chemical_roles` block: `{role_id: {name, parents: [role_id…], level, level_kind, all_identifiers}}` covering the given roles plus every `is_a` ancestor that is a role, excluding `CHEBI:50906`; `level` = shortest `is_a` distance to a branch; `level_kind` = `chebi_role_branch` at level 0 else `chebi_role`; `all_identifiers` = `[chebi.role:<n>, chebi:<n>, chebi:<each secondary>]` sorted.
  - `build_subtree_roles(subtree_strings: dict[str, set[str]], role_map: dict[str, str]) -> dict[str, list[str]]`.

- [ ] **Step 1: Write the failing tests** (`tests/test_mi_roles.py`):

```python
from tests.metabolite_identity_helpers import make_chebi, term
from multiomics_kg.utils.metabolite_identity.roles import (
    build_role_hierarchy, build_subtree_roles, partition_role_strings)

CHEBI = make_chebi(
    CHEBI_50906=term("role"),
    CHEBI_24432=term("biological role", is_a=["CHEBI:50906"]),
    CHEBI_51086=term("chemical role", is_a=["CHEBI:50906"]),
    CHEBI_33232=term("application", is_a=["CHEBI:50906"]),
    CHEBI_26672=term("siderophore", is_a=["CHEBI:24432"]),
    CHEBI_38155=term("phytosiderophore", is_a=["CHEBI:26672"]),
    CHEBI_33229=term("vitamin (role)", alt=["CHEBI:10004"], is_a=["CHEBI:24432", "CHEBI:33232"]),
    CHEBI_17632=term("nitrate"))


def test_partition_detects_secondary_id_roles():
    mols, roles = partition_role_strings(
        {"CHEBI:10004;vitamin", "CHEBI:38155;phytosiderophore", "CHEBI:17632;nitrate"}, CHEBI)
    assert mols == {"CHEBI:17632;nitrate"}
    assert roles == {"CHEBI:10004;vitamin": "chebi.role:33229",
                     "CHEBI:38155;phytosiderophore": "chebi.role:38155"}


def test_hierarchy_ancestors_no_root_levels_and_identifiers():
    h = build_role_hierarchy({"chebi.role:38155", "chebi.role:33229"}, CHEBI)
    assert "chebi.role:50906" not in h
    assert h["chebi.role:24432"]["level"] == 0 and h["chebi.role:24432"]["level_kind"] == "chebi_role_branch"
    assert h["chebi.role:26672"]["level"] == 1 and h["chebi.role:38155"]["level"] == 2
    assert h["chebi.role:38155"]["parents"] == ["chebi.role:26672"]
    assert sorted(h["chebi.role:33229"]["parents"]) == ["chebi.role:24432", "chebi.role:33232"]
    assert h["chebi.role:33229"]["all_identifiers"] == ["chebi.role:33229", "chebi:10004", "chebi:33229"]
    assert all(p in h for v in h.values() for p in v["parents"])


def test_subtree_roles_sorted_nonempty_only():
    sub = build_subtree_roles({"2.A.1": {"CHEBI:10004;vitamin", "CHEBI:17632;nitrate"}, "3.A.1": {"CHEBI:17632;nitrate"}},
                              {"CHEBI:10004;vitamin": "chebi.role:33229"})
    assert sub == {"2.A.1": ["chebi.role:33229"]}
```

- [ ] **Step 2: Run** `uv run pytest tests/test_mi_roles.py -q` → FAIL.

- [ ] **Step 3: Implement `roles.py`.**

```python
"""Stage C: role-typed TCDB substrates → ChemicalRole hierarchy (spec §3.4)."""
from __future__ import annotations

from collections import deque

from multiomics_kg.utils.chebi_utils import ROLE_BRANCHES, ROLE_ROOT, ChebiIndex


def role_node_id(primary: str) -> str:
    return f"chebi.role:{primary.split(':', 1)[1]}"


def _primary_of_string(s: str, chebi: ChebiIndex) -> str | None:
    head = s.partition(";")[0]
    prefix, _, local = head.partition(":")
    if prefix.strip().upper() != "CHEBI" or not local.strip():
        return None
    return chebi.primary(local.strip())


def partition_role_strings(strings, chebi):
    molecules, roles = set(), {}
    for s in sorted(strings):
        p = _primary_of_string(s, chebi)
        if p is not None and chebi.is_role(p):
            roles[s] = role_node_id(p)
        else:
            molecules.add(s)
    return molecules, roles


def _level(p: str, chebi: ChebiIndex) -> int:
    q, seen = deque([(p, 0)]), {p}
    while q:
        x, d = q.popleft()
        if x in ROLE_BRANCHES:
            return d
        for y in (chebi.term(x) or {}).get("is_a", []):
            if y not in seen:
                seen.add(y)
                q.append((y, d + 1))
    raise ValueError(f"role {p} has no path to a ChEBI role branch")


def build_role_hierarchy(role_ids, chebi):
    primaries = {f"CHEBI:{r.split(':', 1)[1]}" for r in role_ids}
    keep = set()
    for p in primaries:
        keep.add(p)
        keep |= {a for a in chebi.ancestors(p) if a != ROLE_ROOT and chebi.is_role(a)}
    out = {}
    for p in sorted(keep, key=lambda c: int(c.split(":")[1])):
        rid = role_node_id(p)
        t = chebi.term(p)
        out[rid] = {
            "name": t["name"],
            "parents": sorted(role_node_id(x) for x in t.get("is_a", []) if x in keep),
            "level": _level(p, chebi),
            "level_kind": "chebi_role_branch" if p in ROLE_BRANCHES else "chebi_role",
            "all_identifiers": sorted([rid, f"chebi:{p.split(':')[1]}"]
                                      + [f"chebi:{a.split(':')[1]}" for a in chebi.secondary_ids(p)]),
        }
    return out


def build_subtree_roles(subtree_strings, role_map):
    out = {}
    for tc, strings in sorted(subtree_strings.items()):
        ids = sorted({role_map[s] for s in strings if s in role_map})
        if ids:
            out[tc] = ids
    return out
```

- [ ] **Step 4: Run** `uv run pytest tests/test_mi_roles.py -q` → PASS.

- [ ] **Step 5: Commit.**

```bash
git add multiomics_kg/utils/metabolite_identity/roles.py tests/test_mi_roles.py
git commit -m "feat(metabolite): ChEBI role partition + pruned ChemicalRole hierarchy"
```

### Task 7: Paper-metabolite resolution (`paper.py`, stages D-harvest and J)

**Files:**
- Create: `multiomics_kg/utils/metabolite_identity/paper.py`, `tests/test_mi_paper.py`
- Modify: `multiomics_kg/utils/metabolite_utils.py` (add `resolve_metabolite_candidates`)

**Interfaces:**
- Consumes: `names.name_key`, `names.is_stereo_qualified`, `names.is_formula_like` (Task 2); `paperconfig_utils.iter_metabolite_assays_tables` (existing).
- Produces:
  - `STEREO_ASSIGNMENTS = ("as_reported", "dominant_form_assumed")`.
  - `@dataclass(frozen=True) Override(id: str, stereo_assignment: str | None, reason: str, raw_name: str)`.
  - `@dataclass PaperRow(paper: str, table: str, row_index: int, name: str, direct_id: str | None, id_method: str | None)`; `name_key` property.
  - `@dataclass PaperTable(paper: str, table: str, csv_path: Path, skip_rows: int, rows: list[PaperRow])`; `key` property `f"{paper}|{table}"`.
  - `load_paper_overrides(pc_path: Path, cfg: dict) -> dict[str, Override]`: union of every distinct `aliases_file` referenced by any of the paper's `metabolite_assays_table` entries; keyed by `name_key`; accepts the legacy `name: "id"` string form (`stereo_assignment=None`) and `{id, stereo_assignment, reason}`; raises `ValueError` on a same-key / different-id pair or an unknown `stereo_assignment`.
  - `harvest_paper_tables(paperconfigs: list[tuple[Path, dict]], project_root: Path) -> tuple[list[PaperTable], dict[str, dict[str, Override]]]` (paper key = paperconfig path relative to `project_root`, today's `per_paper` key).
  - `@dataclass Decision(metabolite_id: str | None, resolution_method: str, stereo_assignment: str | None)`.
  - `@dataclass StereoContext(families: dict[str, frozenset[str]], unspecified: set[str], is_anomer: Callable[[str, str], bool], catalyst_counts: dict[str, int])` with method `variants(node) -> frozenset[str]` (its family minus members that are anomers of `node`).
  - `dominant_form(candidates: set[str], ctx: StereoContext) -> tuple[str | None, dict]` → `(winner or None, {"winner", "counts": {"winner", "runner_up"}, "runner_up_id"})`; `None` when the winner has 0 catalysts, ties, or has < 2× the runner-up.
  - `@dataclass ResolutionResult(decisions: dict[tuple[str, str, int], Decision], lock: dict, drift: list[dict], ambiguous: list[dict], unresolved: list[dict], collisions: list[dict], errors: list[str], stale_overrides: list[dict], stale_lock: list[dict], lock_remapped: list[dict], name_lookup: dict[str, str])`.
  - `resolve_paper_rows(tables, *, overrides: dict[str, dict[str, Override]], name_candidates: dict[str, set[str]], id_index: dict[str, str], ctx: StereoContext, lock: dict, today: str) -> ResolutionResult` (the spec §3.6 order; full logic below).
  - `metabolite_utils.resolve_metabolite_candidates(value: str, conn: sqlite3.Connection) -> list[str]`: every MNXM id the alias step **and** the normalized-name step return (sorted, deduplicated); for `is_formula_like(value)` names it skips the normalized-name step.

- [ ] **Step 1: Write the failing tests** (`tests/test_mi_paper.py`). The `_ctx()` helper wires the Ala / Glu / GlcNAc / Lys families used across the tests.

```python
import pytest

from multiomics_kg.utils.metabolite_identity.paper import (
    Override, PaperRow, PaperTable, StereoContext, dominant_form, load_paper_overrides,
    resolve_paper_rows)

ALA, DALA, ALA_U = "kegg.compound:C00041", "kegg.compound:C00133", "kegg.compound:C01401"
GLU, DGLU, GLU_U = "kegg.compound:C00025", "kegg.compound:C00217", "kegg.compound:C00302"
GLCNAC, BGLCNAC = "kegg.compound:C00140", "kegg.compound:C03878"
LYS, DLYS = "kegg.compound:C00047", "kegg.compound:C00739"
GSH, GSSG = "kegg.compound:C00051", "kegg.compound:C00127"
FAM = {
    ALA: frozenset({ALA, DALA, ALA_U}), DALA: frozenset({ALA, DALA, ALA_U}), ALA_U: frozenset({ALA, DALA, ALA_U}),
    GLU: frozenset({GLU, DGLU, GLU_U}), DGLU: frozenset({GLU, DGLU, GLU_U}), GLU_U: frozenset({GLU, DGLU, GLU_U}),
    GLCNAC: frozenset({GLCNAC, BGLCNAC}), BGLCNAC: frozenset({GLCNAC, BGLCNAC}),
    LYS: frozenset({LYS, DLYS}), DLYS: frozenset({LYS, DLYS}),
    GSH: frozenset({GSH}), GSSG: frozenset({GSSG})}
COUNTS = {ALA: 520, DALA: 103, ALA_U: 0, GLU: 1905, DGLU: 101, GLU_U: 69, GLCNAC: 65, BGLCNAC: 0,
          LYS: 300, DLYS: 10}


def _ctx(counts=COUNTS):
    return StereoContext(families=FAM, unspecified={ALA_U, GLU_U, GLCNAC},
                         is_anomer=lambda a, b: {a, b} == {GLCNAC, BGLCNAC}, catalyst_counts=counts)


def _table(*rows, paper="p", table="t"):
    return PaperTable(paper, table, None, 0,
                      [PaperRow(paper, table, i, n, d, "kegg_direct" if d else None) for i, (n, d) in enumerate(rows)])


IDX = {x: x for x in FAM} | {"kegg.compound:C99999": ALA}


def _run(tables, *, cands=None, overrides=None, lock=None, ctx=None):
    return resolve_paper_rows(tables, overrides=overrides or {}, name_candidates=cands or {},
                              id_index=IDX, ctx=ctx or _ctx(), lock=lock or {"names": {}, "direct_ids": {}},
                              today="2026-10-15")


def test_unqualified_name_goes_to_dominant_form_and_freezes():
    r = _run([_table(("Alanine", None))], cands={"alanine": {DALA, ALA}})
    d = r.decisions[("p", "t", 0)]
    assert (d.metabolite_id, d.resolution_method, d.stereo_assignment) == (ALA, "dominant_stereo_form", "dominant_form_assumed")
    assert r.lock["names"]["alanine"]["id"] == ALA
    assert r.lock["names"]["alanine"]["counts"] == {"winner": 520, "runner_up": 103}


def test_qualified_name_stays():
    r = _run([_table(("D-Alanine", None))], cands={"d-alanine": {DALA}})
    assert r.decisions[("p", "t", 0)].metabolite_id == DALA
    assert r.decisions[("p", "t", 0)].stereo_assignment == "as_reported"


def test_exact_match_to_unspecified_node_goes_through_rule():
    r = _run([_table(("Glutamate", None))], cands={"glutamate": {GLU_U}})
    assert r.decisions[("p", "t", 0)].metabolite_id == GLU


def test_unspecified_node_can_win():
    counts = dict(COUNTS, **{GLU_U: 9000})
    r = _run([_table(("Glutamate", None))], cands={"glutamate": {GLU_U}}, ctx=_ctx(counts))
    d = r.decisions[("p", "t", 0)]
    assert d.metabolite_id == GLU_U and d.stereo_assignment == "dominant_form_assumed"


def test_margin_below_2x_is_ambiguous():
    counts = dict(COUNTS, **{DALA: 300})
    r = _run([_table(("Alanine", None))], cands={"alanine": {ALA, DALA}}, ctx=_ctx(counts))
    assert r.ambiguous and ("p", "t", 0) not in r.decisions


def test_single_specified_match_is_assumed():
    r = _run([_table(("Lysine", None))], cands={"lysine": {LYS}})
    d = r.decisions[("p", "t", 0)]
    assert (d.metabolite_id, d.resolution_method, d.stereo_assignment) == (LYS, "name_match", "dominant_form_assumed")


def test_anomer_only_family_keeps_target_as_reported():
    r = _run([_table(("N-Acetylglucosamine", GLCNAC))])
    d = r.decisions[("p", "t", 0)]
    assert (d.metabolite_id, d.stereo_assignment) == (GLCNAC, "as_reported")


def test_direct_id_canonicalised():
    r = _run([_table(("Whatever", "kegg.compound:C99999"))])
    assert r.decisions[("p", "t", 0)].metabolite_id == ALA


def test_stereo_unspecified_direct_id_goes_through_rule():
    r = _run([_table(("Alanine", ALA_U))])
    d = r.decisions[("p", "t", 0)]
    assert (d.metabolite_id, d.stereo_assignment) == (ALA, "dominant_form_assumed")
    assert r.lock["direct_ids"][ALA_U]["id"] == ALA


def test_override_beats_direct_id():
    ov = {"p": {"glutathione": Override(GSH, "as_reported", "GSH, reduced", "Glutathione")}}
    r = _run([_table(("Glutathione", GSSG))], overrides=ov)
    d = r.decisions[("p", "t", 0)]
    assert (d.metabolite_id, d.resolution_method) == (GSH, "alias_override")


def test_lock_wins_and_drift_reported_not_applied():
    lock = {"names": {"alanine": {"id": DALA, "frozen_on": "2026-01-01",
                                  "counts": {"winner": 1, "runner_up": 0}, "runner_up_id": ALA}},
            "direct_ids": {}}
    r = _run([_table(("Alanine", None))], cands={"alanine": {ALA, DALA}}, lock=lock)
    assert r.decisions[("p", "t", 0)].metabolite_id == DALA
    assert r.drift and r.drift[0]["key"] == "alanine"


def test_lock_target_remapped_through_aliases():
    lock = {"names": {"alanine": {"id": "kegg.compound:C99999", "frozen_on": "x",
                                  "counts": {"winner": 1, "runner_up": 0}, "runner_up_id": None}},
            "direct_ids": {}}
    r = _run([_table(("Alanine", None))], cands={"alanine": {ALA, DALA}}, lock=lock)
    assert r.decisions[("p", "t", 0)].metabolite_id == ALA and r.lock_remapped


def test_stale_lock_entry_reported_not_deleted():
    lock = {"names": {"sarcosine": {"id": ALA, "frozen_on": "x", "counts": {}, "runner_up_id": None}},
            "direct_ids": {}}
    r = _run([_table(("Lysine", None))], cands={"lysine": {LYS}}, lock=lock)
    assert r.stale_lock and "sarcosine" in r.lock["names"]


def test_unmatched_is_unresolved_not_fatal():
    r = _run([_table(("s-(1-2-dicarboxyethyl)glutathione pos", None))])
    assert r.unresolved and not r.ambiguous and ("p", "t", 0) not in r.decisions


def test_distinct_names_same_target_collide_unless_both_overridden():
    t = _table(("Glutathione", GSSG), ("Glutathione Disulfide", GSSG))
    r = _run([t])
    assert r.collisions
    ov = {"p": {"glutathione": Override(GSSG, "as_reported", "", "Glutathione"),
                "glutathione disulfide": Override(GSSG, "as_reported", "", "Glutathione Disulfide")}}
    assert not _run([t], overrides=ov).collisions


def test_case_variants_are_one_key_not_a_collision():
    r = _run([_table(("Glutamic acid", None), ("glutamic acid", None))], cands={"glutamic acid": {GLU, DGLU, GLU_U}})
    assert not r.collisions and r.decisions[("p", "t", 1)].metabolite_id == GLU


def test_two_papers_map_one_name_differently():
    ov = {"a": {"glutathione": Override(GSH, "as_reported", "", "Glutathione")},
          "b": {"glutathione": Override(GSSG, "as_reported", "", "Glutathione")}}
    r = _run([_table(("Glutathione", None), paper="a"), _table(("Glutathione", None), paper="b")], overrides=ov)
    assert r.decisions[("a", "t", 0)].metabolite_id == GSH and r.decisions[("b", "t", 0)].metabolite_id == GSSG


def test_legacy_override_on_stereo_target_needs_assignment():
    ov = {"p": {"alanine": Override(ALA, None, "", "Alanine")}}
    r = _run([_table(("Alanine", None))], overrides=ov)
    assert r.errors and "stereo_assignment" in r.errors[0]


def test_stale_alias_override_reported():
    ov = {"p": {"cystine": Override(GSSG, "as_reported", "", "Cystine")}}
    r = _run([_table(("Lysine", None))], cands={"lysine": {LYS}}, overrides=ov)
    assert r.stale_overrides == [{"paper": "p", "name_key": "cystine", "id": GSSG}]


def test_unresolvable_override_target_fails():
    ov = {"p": {"lysine": Override("kegg.compound:C00000", "as_reported", "", "Lysine")}}
    r = _run([_table(("Lysine", None))], overrides=ov)
    assert r.errors


def test_dominant_form_tie_refused():
    counts = dict(COUNTS, **{DALA: 520})
    winner, info = dominant_form({ALA, DALA}, _ctx(counts))
    assert winner is None and info["counts"] == {"winner": 520, "runner_up": 520}


def test_load_overrides_union_and_key_collision(tmp_path):
    (tmp_path / "a.yaml").write_text('Glutathione: {id: "kegg.compound:C00051", stereo_assignment: as_reported, reason: GSH}\n'
                                     'CO2: "kegg.compound:C00011"\nCo2: "kegg.compound:C00175"\n')
    cfg = {"publication": {"supplementary_materials": {
        "t1": {"type": "metabolite_assays_table", "aliases_file": "a.yaml"},
        "t2": {"type": "metabolite_assays_table"}}}}
    ov = load_paper_overrides(tmp_path / "paperconfig.yaml", cfg)
    assert ov["glutathione"].id == "kegg.compound:C00051"
    assert ov["CO2"].id != ov["Co2"].id
    (tmp_path / "a.yaml").write_text('Glutathione: "kegg.compound:C00051"\nglutathione: "kegg.compound:C00127"\n')
    with pytest.raises(ValueError):
        load_paper_overrides(tmp_path / "paperconfig.yaml", cfg)


def test_overrides_apply_to_every_table_of_the_paper():
    """Review focus 1: Kujawinski references its alias file from 1 of 7 tables."""
    ov = {"p": {"glutathione": Override(GSH, "as_reported", "", "Glutathione")}}
    r = _run([_table(("Glutathione", GSSG), table="kegg_table"), _table(("Glutathione", None), table="flags")],
             overrides=ov)
    assert r.decisions[("p", "kegg_table", 0)].metabolite_id == GSH
    assert r.decisions[("p", "flags", 0)].metabolite_id == GSH
```

  Add to the same file a resolver test against a stub DB (copy `_make_stub_resolver` from `tests/test_build_kegg_metabolism_xrefs_metabolomics.py` and add `INSERT INTO compound_names VALUES ('MNXM13', 'co')`, `INSERT INTO compound_aliases VALUES ('MNXM13','name','CO')`):

```python
def test_resolver_skips_lowercased_step_for_formula_names(tmp_path):
    from multiomics_kg.utils.metabolite_utils import resolve_metabolite_candidates
    conn = _stub(tmp_path)          # helper defined above in this file
    assert resolve_metabolite_candidates("CO", conn) == ["MNXM13"]      # exact alias hit
    assert resolve_metabolite_candidates("Co", conn) == []              # no alias, name step skipped
    assert resolve_metabolite_candidates("glucose", conn) == ["MNXM41"]
```

- [ ] **Step 2: Run** `uv run pytest tests/test_mi_paper.py -q` → FAIL.

- [ ] **Step 3: Implement `resolve_metabolite_candidates`** in `multiomics_kg/utils/metabolite_utils.py` below `resolve_metabolite`:

```python
def resolve_metabolite_candidates(value: str, conn: sqlite3.Connection) -> list[str]:
    """Every MNXM id a name or id can mean — alias hits plus normalized-name hits.

    Unlike resolve_metabolite (first hit wins), this never picks: the caller decides
    (spec 2026-10-13 §3.6). Names that parse as formulas (CO, Co, NO) are case-significant,
    so the lower-cased normalized-name step is skipped for them.
    """
    from multiomics_kg.utils.metabolite_identity.names import is_formula_like
    v = value.strip()
    if _MNXM_RE.match(v) if "_MNXM_RE" in globals() else v.startswith("MNXM"):
        row = conn.execute("SELECT mnxm_id FROM compounds WHERE mnxm_id = ?", (v,)).fetchone()
        if row:
            return [row[0]]
    hits = {r[0] for r in conn.execute(
        "SELECT DISTINCT mnxm_id FROM compound_aliases WHERE value = ?", (v,))}
    if not is_formula_like(v):
        hits |= {r[0] for r in conn.execute(
            "SELECT DISTINCT mnxm_id FROM compound_names WHERE name_normalized = ?",
            (_normalize_name(v),))}
    return sorted(hits)
```

  Replace the `_MNXM_RE` conditional with whatever MNXM test `resolve_metabolite` already uses (read lines 53–83 of the file; reuse its exact check), so there's one rule.

- [ ] **Step 4: Implement `paper.py`.**

```python
"""Paper-metabolite resolution (spec §3.6, §3.6a, §6.5, §6.6). Stages D (harvest) and J."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import pandas as pd
import yaml

from multiomics_kg.utils.metabolite_identity.names import is_stereo_qualified, name_key
from multiomics_kg.utils.paperconfig_utils import iter_metabolite_assays_tables

STEREO_ASSIGNMENTS = ("as_reported", "dominant_form_assumed")
_NO_ID = {"0", "-", "na", "n/a", "none"}
_ID_METHOD = {"kegg.compound": "kegg_direct", "chebi": "chebi_direct", "mnx": "mnx_direct"}


@dataclass(frozen=True)
class Override:
    id: str
    stereo_assignment: str | None
    reason: str
    raw_name: str


@dataclass
class PaperRow:
    paper: str
    table: str
    row_index: int
    name: str
    direct_id: str | None
    id_method: str | None

    @property
    def name_key(self) -> str:
        return name_key(self.name)


@dataclass
class PaperTable:
    paper: str
    table: str
    csv_path: Path | None
    skip_rows: int
    rows: list[PaperRow]
    name_col: str = "compound"

    @property
    def key(self) -> str:
        return f"{self.paper}|{self.table}"


def load_paper_overrides(pc_path: Path, cfg: dict) -> dict[str, Override]:
    files = sorted({e.get("aliases_file") for _, e in iter_metabolite_assays_tables(cfg)
                    if e.get("aliases_file")})
    out: dict[str, Override] = {}
    for f in files:
        path = pc_path.parent / f
        if not path.exists():
            continue
        doc = yaml.safe_load(path.read_text()) or {}
        for raw, val in doc.items():
            if isinstance(val, dict):
                ov = Override(str(val["id"]), val.get("stereo_assignment"), str(val.get("reason", "")), str(raw))
            else:
                ov = Override(str(val), None, "", str(raw))
            if ov.stereo_assignment not in (None, *STEREO_ASSIGNMENTS):
                raise ValueError(f"{path}: {raw!r} stereo_assignment {ov.stereo_assignment!r} "
                                 f"not in {STEREO_ASSIGNMENTS}")
            k = name_key(str(raw))
            if k in out and out[k].id != ov.id:
                raise ValueError(f"{path}: {raw!r} and {out[k].raw_name!r} share name_key {k!r} "
                                 f"but point to {ov.id} vs {out[k].id}")
            out[k] = ov
    return out


def harvest_paper_tables(paperconfigs, project_root: Path):
    tables, overrides = [], {}
    for pc_path, cfg in paperconfigs:
        try:
            paper = str(Path(pc_path).resolve().relative_to(project_root.resolve()))
        except ValueError:
            paper = str(pc_path)
        overrides[paper] = load_paper_overrides(Path(pc_path), cfg)
        for entry_key, entry in iter_metabolite_assays_tables(cfg):
            csv = entry.get("filename")
            if not csv:
                continue
            csv_full = Path(csv) if Path(csv).is_absolute() else project_root / csv
            skip = int(entry.get("skip_rows", 0) or 0)
            df = pd.read_csv(csv_full, dtype=str, keep_default_na=False, skiprows=skip or None)
            name_col = entry.get("name_col") or "compound"
            id_col, id_type = entry.get("id_col") or "", entry.get("id_type") or ""
            rows = []
            for i, row in df.iterrows():
                name = (row.get(name_col) or "").strip()
                if not name:
                    continue
                direct, method = None, None
                if id_col and id_type:
                    raw = (row.get(id_col) or "").strip()
                    if raw and raw.lower() not in _NO_ID:
                        direct = raw if ":" in raw else f"{id_type}:{raw}"
                        method = _ID_METHOD.get(id_type, "id_direct")
                rows.append(PaperRow(paper, entry_key, int(i), name, direct, method))
            tables.append(PaperTable(paper, entry_key, csv_full, skip, rows, name_col))
    return tables, overrides


@dataclass
class Decision:
    metabolite_id: str | None
    resolution_method: str
    stereo_assignment: str | None


@dataclass
class StereoContext:
    families: dict[str, frozenset[str]]
    unspecified: set[str]
    is_anomer: Callable[[str, str], bool]
    catalyst_counts: dict[str, int]

    def variants(self, node: str) -> frozenset[str]:
        fam = self.families.get(node, frozenset({node}))
        return frozenset(m for m in fam if m == node or not self.is_anomer(m, node))


def dominant_form(candidates, ctx: StereoContext):
    pool = set()
    for c in candidates:
        pool |= ctx.variants(c)
    pool = {m for m in pool if m in candidates or not any(ctx.is_anomer(m, c) for c in candidates)}
    ranked = sorted(pool, key=lambda m: (-ctx.catalyst_counts.get(m, 0), m))
    win = ranked[0]
    run = ranked[1] if len(ranked) > 1 else None
    wc, rc = ctx.catalyst_counts.get(win, 0), ctx.catalyst_counts.get(run, 0) if run else 0
    info = {"winner": win, "counts": {"winner": wc, "runner_up": rc}, "runner_up_id": run}
    ok = wc > 0 and wc >= 2 * rc and (run is None or wc != rc)
    return (win if ok else None), info


@dataclass
class ResolutionResult:
    decisions: dict = field(default_factory=dict)
    lock: dict = field(default_factory=dict)
    drift: list = field(default_factory=list)
    ambiguous: list = field(default_factory=list)
    unresolved: list = field(default_factory=list)
    collisions: list = field(default_factory=list)
    errors: list = field(default_factory=list)
    stale_overrides: list = field(default_factory=list)
    stale_lock: list = field(default_factory=list)
    lock_remapped: list = field(default_factory=list)
    name_lookup: dict = field(default_factory=dict)


def resolve_paper_rows(tables, *, overrides, name_candidates, id_index, ctx, lock, today):
    res = ResolutionResult()
    lock = {"names": dict(lock.get("names", {})), "direct_ids": dict(lock.get("direct_ids", {}))}
    # Re-map lock keys / targets through the current aliases (spec §3.6a).
    for section in ("names", "direct_ids"):
        for key in list(lock[section]):
            entry = dict(lock[section][key])
            new_key = id_index.get(key, key) if section == "direct_ids" else key
            tgt = id_index.get(entry["id"])
            if tgt is None:
                res.stale_lock.append({"section": section, "key": key, "id": entry["id"], "why": "target gone"})
                continue
            if tgt != entry["id"]:
                res.lock_remapped.append({"section": section, "key": key, "old": entry["id"], "new": tgt})
                entry["id"] = tgt
            del lock[section][key]
            lock[section][new_key] = entry
    used_keys = {"names": set(), "direct_ids": set()}
    used_overrides: set[tuple[str, str]] = set()

    def has_variants(node):
        return len(ctx.variants(node)) > 1

    def rule(section, key, candidates):
        if key in lock[section]:
            used_keys[section].add(key)
            entry = lock[section][key]
            win, info = dominant_form(candidates, ctx)
            if win != entry["id"]:
                res.drift.append({"section": section, "key": key, "frozen": entry["id"],
                                  "recomputed": win, "counts": info["counts"]})
            return entry["id"]
        win, info = dominant_form(candidates, ctx)
        if win is None:
            return None
        used_keys[section].add(key)
        lock[section][key] = {"id": win, "frozen_on": today, "counts": info["counts"],
                              "runner_up_id": info["runner_up_id"]}
        return win

    for t in tables:
        for row in t.rows:
            k = row.name_key
            key3 = (row.paper, row.table, row.row_index)
            ov = overrides.get(row.paper, {}).get(k)
            if ov is not None:
                used_overrides.add((row.paper, k))
                tgt = id_index.get(ov.id)
                if tgt is None:
                    res.errors.append(f"{row.paper}: alias {ov.raw_name!r} → {ov.id} resolves to no compound")
                    continue
                sa = ov.stereo_assignment
                if sa is None:
                    if has_variants(tgt):
                        res.errors.append(f"{row.paper}: alias {ov.raw_name!r} → {tgt} has stereo variants; "
                                          f"add an explicit stereo_assignment")
                        continue
                    sa = "as_reported"
                res.decisions[key3] = Decision(tgt, "alias_override", sa)
                continue
            if row.direct_id:
                tgt = id_index.get(row.direct_id)
                if tgt is None:
                    res.unresolved.append({"paper": row.paper, "table": row.table, "name": row.name,
                                           "direct_id": row.direct_id})
                    continue
                if tgt in ctx.unspecified and has_variants(tgt):
                    chosen = rule("direct_ids", tgt, {tgt})
                    if chosen is None:
                        res.ambiguous.append({"paper": row.paper, "table": row.table, "name": row.name,
                                              "candidates": sorted(ctx.variants(tgt))})
                        continue
                    method = "dominant_stereo_form" if chosen != tgt else row.id_method
                    res.decisions[key3] = Decision(chosen, method, "dominant_form_assumed")
                else:
                    res.decisions[key3] = Decision(tgt, row.id_method, "as_reported")
                continue
            cands = {id_index.get(c, c) for c in name_candidates.get(k, set())}
            qualified = is_stereo_qualified(row.name)
            if k in lock["names"]:
                chosen = rule("names", k, cands or {lock["names"][k]["id"]})
                res.decisions[key3] = Decision(chosen, "dominant_stereo_form", "dominant_form_assumed")
                continue
            if not cands:
                res.unresolved.append({"paper": row.paper, "table": row.table, "name": row.name})
                continue
            fams = {ctx.families.get(c, frozenset({c})) for c in cands}
            if len(cands) == 1:
                (c,) = cands
                if not has_variants(c) or qualified:
                    res.decisions[key3] = Decision(c, "name_match", "as_reported")
                elif c in ctx.unspecified:
                    chosen = rule("names", k, {c})
                    if chosen is None:
                        res.ambiguous.append({"paper": row.paper, "table": row.table, "name": row.name,
                                              "candidates": sorted(ctx.variants(c))})
                        continue
                    method = "dominant_stereo_form" if chosen != c else "name_match"
                    res.decisions[key3] = Decision(chosen, method, "dominant_form_assumed")
                else:
                    res.decisions[key3] = Decision(c, "name_match", "dominant_form_assumed")
                continue
            if len(fams) == 1 and not qualified:
                chosen = rule("names", k, cands)
                if chosen is not None:
                    res.decisions[key3] = Decision(chosen, "dominant_stereo_form", "dominant_form_assumed")
                    continue
            res.ambiguous.append({"paper": row.paper, "table": row.table, "name": row.name,
                                  "candidates": sorted(cands)})

    # Same-table collisions (spec §3.6): distinct name_keys on one node need both overridden.
    for t in tables:
        by_node: dict[str, set[str]] = defaultdict(set)
        for row in t.rows:
            d = res.decisions.get((row.paper, row.table, row.row_index))
            if d and d.metabolite_id:
                by_node[d.metabolite_id].add(row.name_key)
        for node, keys in sorted(by_node.items()):
            if len(keys) < 2:
                continue
            ovs = overrides.get(t.paper, {})
            if all(k in ovs and id_index.get(ovs[k].id) == node for k in keys):
                continue
            res.collisions.append({"paper": t.paper, "table": t.table, "node": node, "names": sorted(keys)})

    seen_names = {row.name_key for t in tables for row in t.rows}
    seen_ids = {id_index.get(row.direct_id) for t in tables for row in t.rows if row.direct_id}
    for k in sorted(set(lock["names"]) - seen_names):
        res.stale_lock.append({"section": "names", "key": k, "id": lock["names"][k]["id"], "why": "no row"})
    for k in sorted(set(lock["direct_ids"]) - seen_ids):
        res.stale_lock.append({"section": "direct_ids", "key": k, "id": lock["direct_ids"][k]["id"], "why": "no row"})
    for paper, ovs in sorted(overrides.items()):
        for k, ov in sorted(ovs.items()):
            if (paper, k) not in used_overrides:
                res.stale_overrides.append({"paper": paper, "name_key": k, "id": ov.id})
    for k, cands in sorted(name_candidates.items()):
        canon = {id_index.get(c, c) for c in cands}
        if len(canon) == 1:
            res.name_lookup[k] = next(iter(canon))
    res.lock = lock
    return res
```

  Two behaviours the tests pin and the code must keep: stale lock entries stay in `res.lock` (never deleted; a human deletes them), and stale alias overrides keep no node alive (stage K only keeps final targets).

- [ ] **Step 5: Run** `uv run pytest tests/test_mi_paper.py -q` → PASS.

- [ ] **Step 6: Commit.**

```bash
git add multiomics_kg/utils/metabolite_identity/paper.py multiomics_kg/utils/metabolite_utils.py tests/test_mi_paper.py
git commit -m "feat(metabolite): deterministic paper-name resolution, dominant stereo form, lock"
```

### Task 8: KEGG compound classes (`classes.py`, stage M)

**Files:**
- Create: `multiomics_kg/utils/metabolite_identity/classes.py`, `tests/test_mi_classes.py`

**Interfaces:**
- Produces:
  - `COMPOUND_BRITE_TREES = {"br08001": "compounds_with_biological_roles", "br08002": "lipids", "br08003": "phytochemical_compounds", "br08005": "bioactive_peptides", "br08006": "endocrine_disrupting_compounds", "br08007": "pesticides", "br08009": "natural_toxins", "br08021": "glycosides"}`.
  - `LEVEL_KINDS = ("kegg_compound_class", "kegg_compound_subclass", "kegg_compound_group", "kegg_compound_subgroup")`.
  - `BR08001_LEVEL0 = ("Organic acids", "Lipids", "Carbohydrates", "Nucleic acids", "Peptides", "Vitamins and cofactors", "Steroids", "Hormones and transmitters", "Antibiotics")` (measured 2026-10-14; used by the vocabulary + KG test).
  - `parse_compound_brite(tree_code: str, raw: dict) -> tuple[dict[str, dict], list[tuple[str, str]]]` → `({class_id: {parent, name, level}}, [(kegg_C, class_id)])`; class id `kegg.brite:<tree>.A<i>.B<j>…` (1-based positional, the `BriteCategory` scheme); name = `html.unescape`, then strip `<…>` tags and ` [Fig]` / ` [BR:…]` / ` [PATH:…]` suffixes; raises `ValueError` when `raw` has children but 0 compound leaves (fail-loud, `_require_parsed` analogue) or depth > 4.
  - `prune_compound_brite(tree_code, nodes, members, *, canonical_of: dict[str, str], written: set[str]) -> dict` → spec §6.2 shape `{tree_code, tree, nodes: [{id, parent, name, level}], members: [{metabolite_id, class_id}]}` (sorted, deduplicated, ancestors-only).
  - `download_compound_brite(raw_dir: Path, force: bool) -> dict[str, dict]` (per-tree raw JSON at `raw_dir / f"br_{code}.json"`, via `kegg_utils._fetch_json`).

- [ ] **Step 1: Write the failing tests** (`tests/test_mi_classes.py`):

```python
import pytest

from multiomics_kg.utils.metabolite_identity.classes import parse_compound_brite, prune_compound_brite

RAW = {"name": "br08001", "children": [
    {"name": "Organic acids", "children": [
        {"name": "Carboxylic acids [Fig]", "children": [
            {"name": "Monocarboxylic acids", "children": [
                {"name": "C00058  Formate; Methanoate"}, {"name": "C00033  Acetate"}]}]}]},
    {"name": "Peptides", "children": [
        {"name": "Amino acids", "children": [{"name": "C00041  L-Alanine"}, {"name": "C00033  Acetate"}]}]}]}


def test_parse_ids_names_levels():
    nodes, members = parse_compound_brite("br08001", RAW)
    assert nodes["kegg.brite:br08001.A1"] == {"parent": None, "name": "Organic acids", "level": 0}
    assert nodes["kegg.brite:br08001.A1.B1"]["name"] == "Carboxylic acids"
    assert ("C00041", "kegg.brite:br08001.A2.B1") in members
    assert [m for m in members if m[0] == "C00033"] == [("C00033", "kegg.brite:br08001.A1.B1.C1"),
                                                        ("C00033", "kegg.brite:br08001.A2.B1")]


def test_parse_fails_loud_on_empty():
    with pytest.raises(ValueError):
        parse_compound_brite("br08001", {"name": "br08001", "children": [{"name": "X", "children": []}]})


def test_prune_ancestors_only_and_alias_lands_on_canonical():
    nodes, members = parse_compound_brite("br08001", RAW)
    out = prune_compound_brite("br08001", nodes, members,
                               canonical_of={"kegg.compound:C00058": "kegg.compound:C00033"},
                               written={"kegg.compound:C00033"})
    assert {n["id"] for n in out["nodes"]} == {"kegg.brite:br08001.A1", "kegg.brite:br08001.A1.B1",
                                                "kegg.brite:br08001.A1.B1.C1", "kegg.brite:br08001.A2",
                                                "kegg.brite:br08001.A2.B1"}
    assert {"metabolite_id": "kegg.compound:C00033", "class_id": "kegg.brite:br08001.A1.B1.C1"} in out["members"]
    assert len(out["members"]) == len({(m["metabolite_id"], m["class_id"]) for m in out["members"]})
    assert out["tree"] == "compounds_with_biological_roles"
```

- [ ] **Step 2: Run** `uv run pytest tests/test_mi_classes.py -q` → FAIL.

- [ ] **Step 3: Implement `classes.py`.**

```python
"""Stage M: KEGG compound BRITE → MetaboliteClass trees, pruned to written metabolites (spec §3.7, §6.2)."""
from __future__ import annotations

import html
import json
import re
from pathlib import Path

from multiomics_kg.utils import kegg_utils

COMPOUND_BRITE_TREES = {
    "br08001": "compounds_with_biological_roles", "br08002": "lipids",
    "br08003": "phytochemical_compounds", "br08005": "bioactive_peptides",
    "br08006": "endocrine_disrupting_compounds", "br08007": "pesticides",
    "br08009": "natural_toxins", "br08021": "glycosides"}
LEVEL_KINDS = ("kegg_compound_class", "kegg_compound_subclass", "kegg_compound_group", "kegg_compound_subgroup")
BR08001_LEVEL0 = ("Organic acids", "Lipids", "Carbohydrates", "Nucleic acids", "Peptides",
                  "Vitamins and cofactors", "Steroids", "Hormones and transmitters", "Antibiotics")
_LEAF = re.compile(r"^(C\d{5})\s")
_SUFFIX = re.compile(r"\s*\[(?:Fig|BR:[^\]]*|PATH:[^\]]*)\]")
_TAG = re.compile(r"<[^>]+>")
_URL = "https://rest.kegg.jp/get/br:{code}/json"


def _clean(name: str) -> str:
    return _SUFFIX.sub("", _TAG.sub("", html.unescape(name))).strip()


def parse_compound_brite(tree_code, raw):
    nodes, members = {}, []

    def walk(children, path, parent):
        for i, c in enumerate(children, 1):
            m = _LEAF.match(c.get("name", ""))
            if m and not c.get("children"):
                if parent:
                    members.append((m.group(1), parent))
                continue
            if len(path) >= len(LEVEL_KINDS):
                raise ValueError(f"{tree_code}: class depth > {len(LEVEL_KINDS)} at {c.get('name')!r}")
            seg = path + [f"{'ABCD'[len(path)]}{i}"]
            cid = f"kegg.brite:{tree_code}." + ".".join(seg)
            nodes[cid] = {"parent": parent, "name": _clean(c["name"]), "level": len(path)}
            walk(c.get("children", []), seg, cid)

    walk(raw.get("children", []), [], None)
    if raw.get("children") and not members:
        raise ValueError(f"{tree_code}: parsed 0 compound members (KEGG BRITE format drift?)")
    return nodes, members


def prune_compound_brite(tree_code, nodes, members, *, canonical_of, written):
    kept_members = set()
    for c, cls in members:
        mid = canonical_of.get(f"kegg.compound:{c}", f"kegg.compound:{c}")
        if mid in written:
            kept_members.add((mid, cls))
    keep = set()
    for _, cls in kept_members:
        while cls and cls not in keep:
            keep.add(cls)
            cls = nodes[cls]["parent"]
    return {
        "tree_code": tree_code, "tree": COMPOUND_BRITE_TREES[tree_code],
        "nodes": [{"id": i, "parent": nodes[i]["parent"], "name": nodes[i]["name"], "level": nodes[i]["level"]}
                  for i in sorted(keep)],
        "members": [{"metabolite_id": m, "class_id": c} for m, c in sorted(kept_members)],
    }


def download_compound_brite(raw_dir: Path, force: bool = False) -> dict[str, dict]:
    raw_dir.mkdir(parents=True, exist_ok=True)
    out = {}
    for code in COMPOUND_BRITE_TREES:
        path = raw_dir / f"br_{code}.json"
        if path.exists() and not force:
            out[code] = json.loads(path.read_text())
            continue
        data = kegg_utils._fetch_json(_URL.format(code=code))
        path.write_text(json.dumps(data))
        out[code] = data
    return out
```

- [ ] **Step 4: Run** `uv run pytest tests/test_mi_classes.py -q` → PASS.

- [ ] **Step 5: Commit.**

```bash
git add multiomics_kg/utils/metabolite_identity/classes.py tests/test_mi_classes.py
git commit -m "feat(metabolite): KEGG compound BRITE parse + canonical pruning"
```

### Task 9: Lumping pin (`lumping.py`, stage O)

**Files:**
- Create: `multiomics_kg/utils/metabolite_identity/lumping.py`, `tests/test_mi_lumping.py`

**Interfaces:**
- Produces:
  - `SEED_LUMPING_IDS = ("2.A.1", "2.A.1.1", "2.A.1.2", "2.A.1.3", "2.A.3", "2.A.6", "2.A.6.2", "2.A.7", "2.A.66", "2.A.66.1", "3.A.1", "3.A.1.1", "3.A.3")` (re-confirmed live 2026-10-14).
  - `LUMPING_MIN_LEVEL = 2`, `LUMPING_MIN_METABOLITES = 50`.
  - `threshold_lumping(kept: set[str], levels: dict[str, int], substrate_counts: dict[str, int]) -> set[str]`.
  - `update_lumping_pin(pin: dict | None, kept: set[str], levels, substrate_counts, today: str) -> tuple[dict, dict]` → `(pin, {"lumping_set_drift": [...], "lumping_set_added": [...], "seeded": bool})`.
  - `substrate_breadth(tc_id: str, level: int, pin: dict) -> str` → `lumping` | `specific` | `not_assessed` (shared by the adapter, Task 17).

- [ ] **Step 1: Write the failing tests** (`tests/test_mi_lumping.py`):

```python
from multiomics_kg.utils.metabolite_identity.lumping import (
    SEED_LUMPING_IDS, substrate_breadth, update_lumping_pin)

KEPT = set(SEED_LUMPING_IDS) | {"2.A", "9.B.1"}
LEVELS = {i: i.count(".") for i in KEPT}
COUNTS = {i: 60 for i in SEED_LUMPING_IDS} | {"3.A.3": 47, "9.B.1": 3, "2.A": 900}


def test_first_run_seeds_literal_13_not_threshold():
    pin, rep = update_lumping_pin(None, KEPT, LEVELS, COUNTS, "2026-10-15")
    assert pin["lumping_tcdb_ids"] == sorted(SEED_LUMPING_IDS) and rep["seeded"]
    assert "3.A.3" in pin["lumping_tcdb_ids"]                 # 47 molecules, still pinned
    assert pin["frozen_kept_ids"] == sorted(KEPT)


def test_frozen_id_drift_reported_not_applied():
    pin, _ = update_lumping_pin(None, KEPT, LEVELS, COUNTS, "2026-10-15")
    pin2, rep = update_lumping_pin(pin, KEPT, LEVELS, dict(COUNTS, **{"9.B.1": 99}), "2026-11-01")
    assert "9.B.1" not in pin2["lumping_tcdb_ids"]
    assert [d["tc_id"] for d in rep["lumping_set_drift"]] == ["3.A.3", "9.B.1"]
    assert pin2["frozen_on"] == "2026-10-15"


def test_newly_kept_id_auto_classified_and_appended():
    pin, _ = update_lumping_pin(None, KEPT, LEVELS, COUNTS, "2026-10-15")
    kept = KEPT | {"2.A.99"}
    pin2, rep = update_lumping_pin(pin, kept, dict(LEVELS, **{"2.A.99": 2}), dict(COUNTS, **{"2.A.99": 300}), "x")
    assert "2.A.99" in pin2["lumping_tcdb_ids"] and rep["lumping_set_added"] == ["2.A.99"]
    assert "2.A.99" not in pin2["frozen_kept_ids"]


def test_substrate_breadth():
    pin, _ = update_lumping_pin(None, KEPT, LEVELS, COUNTS, "d")
    assert substrate_breadth("3.A.3", 2, pin) == "lumping"
    assert substrate_breadth("9.B.1", 2, pin) == "specific"
    assert substrate_breadth("2.A", 1, pin) == "not_assessed"
```

- [ ] **Step 2: Run** → FAIL. **Step 3: Implement `lumping.py`.**

```python
"""Stage O: the pinned TCDB substrate-lumping set (spec §3.4a). A documented R3 exception."""
from __future__ import annotations

SEED_LUMPING_IDS = ("2.A.1", "2.A.1.1", "2.A.1.2", "2.A.1.3", "2.A.3", "2.A.6", "2.A.6.2",
                    "2.A.7", "2.A.66", "2.A.66.1", "3.A.1", "3.A.1.1", "3.A.3")
LUMPING_MIN_LEVEL = 2
LUMPING_MIN_METABOLITES = 50


def threshold_lumping(kept, levels, substrate_counts):
    return {t for t in kept if levels.get(t, 0) >= LUMPING_MIN_LEVEL
            and substrate_counts.get(t, 0) >= LUMPING_MIN_METABOLITES}


def update_lumping_pin(pin, kept, levels, substrate_counts, today):
    seeded = pin is None
    if seeded:
        pin = {"frozen_on": today, "frozen_kept_ids": sorted(kept),
               "lumping_tcdb_ids": sorted(SEED_LUMPING_IDS)}
    thr = threshold_lumping(kept, levels, substrate_counts)
    frozen = set(pin["frozen_kept_ids"])
    current = set(pin["lumping_tcdb_ids"])
    judged = frozen | (current - frozen)          # frozen ids + ids appended on earlier builds
    drift = [{"tc_id": t, "pinned": "lumping" if t in current else "specific",
              "threshold": "lumping" if t in thr else "specific",
              "metabolite_count": substrate_counts.get(t, 0)}
             for t in sorted(kept & judged) if (t in current) != (t in thr)]
    added = sorted(t for t in kept - judged if t in thr)
    pin = dict(pin, lumping_tcdb_ids=sorted(current | set(added)))
    return pin, {"lumping_set_drift": drift, "lumping_set_added": added, "seeded": seeded}


def substrate_breadth(tc_id, level, pin):
    if level < LUMPING_MIN_LEVEL:
        return "not_assessed"
    return "lumping" if tc_id in set(pin.get("lumping_tcdb_ids", [])) else "specific"
```

- [ ] **Step 4: Run** `uv run pytest tests/test_mi_lumping.py -q` → PASS. **Step 5: Commit.**

```bash
git add multiomics_kg/utils/metabolite_identity/lumping.py tests/test_mi_lumping.py
git commit -m "feat(tcdb): pinned substrate-lumping set with drift / auto-append"
```

### Task 10: Stage-P validation and outputs (`outputs.py`)

**Files:**
- Create: `multiomics_kg/utils/metabolite_identity/outputs.py`, `tests/test_mi_outputs.py`
- Modify: `multiomics_kg/utils/metabolite_identity/__init__.py` (export every public name from Tasks 3–10)

**Interfaces:**
- Consumes: `PaperTable`, `Decision`, `ResolutionResult` (Task 7); `IdentityResult` (Task 3).
- Produces:
  - `validate_build(*, table: dict[str, dict], roles: dict[str, dict], subtree_substrates: dict[str, list[str]], subtree_roles: dict[str, list[str]], class_trees: dict[str, dict], resolution: ResolutionResult, chebi_missing_needed: set[str]) -> list[str]` (empty = OK). Checks (spec §5.2 stage P): ambiguous / collisions / resolution errors; an `all_identifiers` entry on two nodes across `Metabolite` ∪ `ChemicalRole`; an InChIKey on two nodes; dangling `stereo_form_of` / `enantiomer_of` / class member / class parent / role parent / `subtree_roles` / `subtree_substrates` targets; a resolved decision whose `metabolite_id` isn't in `table` or has no `stereo_assignment`; `chebi_missing_needed`.
  - `write_all_atomic(files: dict[Path, str]) -> None`: writes every file to `<name>.tmp`, then `os.replace` each; on an exception before the replace loop, deletes the `.tmp` files and re-raises.
  - `write_resolved_tables(tables: list[PaperTable], decisions: dict) -> dict[Path, str]` (spec name): per table, re-reads the CSV with step 7's arguments (`dtype=str, keep_default_na=False, skiprows=skip_rows or None`), appends `metabolite_id`, `resolution_method`, `stereo_assignment` (unresolved rows: `""`, `unresolved`, `""`), and returns `{<stem>_resolved.csv: csv_text, <stem>_resolution_report.json: json_text}` (report keys as step 7: `source_csv, resolved_csv, total_rows, resolved, unresolved, match_rate, method_counts`).
  - `build_mapping_v2(*, id_aliases: dict[str, str], resolution: ResolutionResult, overrides: dict, tables, compounds_summary: dict, per_paper_stats: dict) -> dict` (spec §6.5; `schema_version: 2`; no `specific_lookup` / `multi_lookup`).
  - `repoint_reason(*, old_id: str | None, new: Decision, row, identity: IdentityResult, rekeys: dict[str, str], dropped: dict[str, str], case_fold_keys: set[str]) -> str | None` (spec §12 order: `alias_override`, `dominant_stereo_form`, `case_fold_convergence`, `direct_id_canonicalised`, `placeholder_recurated`, `remeasured_elsewhere`, `merged`; `None` when `old_id == new.metabolite_id`).
  - `remap_rows(*, old_ids: set[str], table: dict[str, dict], identity: IdentityResult, rekeys: dict[str, str], role_of: dict[str, str], dropped: dict[str, str], normalize: Callable[[str], str]) -> tuple[list[dict], list[str]]` (`normalize` maps an old `chebi:` id to `chebi:<primary>` and leaves other ids unchanged; stage D keyed every node by its primary id, so `rekeys` / `member_basis` / `dropped` use normalised ids) → `(rows, unaccounted)`; each row `{old_id, new_id, basis, role_id}`; basis ∈ merge basis | `rekey_chebi_kegg_xref` | `removed_role` | `removed_placeholder` | `remeasured_elsewhere` | `chebi_primary_id` (plan "Spec gaps" 2); `unaccounted` lists old ids none of these explain (Task 14 requires it empty).
  - `build_identity_report(**sections) -> dict`: a sorted-keys dict of the sections listed in spec §6.7 plus `remap`, `repoints`, `unresolved_paper_names`, `structure_counts`, `stereo_counts`, `lumping`, `rekeys`, `dropped`, `chebi_release`.

- [ ] **Step 1: Write the failing tests** (`tests/test_mi_outputs.py`):

```python
import json

import pandas as pd

from tests.metabolite_identity_helpers import cpd
from multiomics_kg.utils.metabolite_identity.canonical import IdentityResult
from multiomics_kg.utils.metabolite_identity.outputs import (
    remap_rows, repoint_reason, validate_build, write_all_atomic, write_resolved_tables)
from multiomics_kg.utils.metabolite_identity.paper import Decision, PaperRow, PaperTable, ResolutionResult


def _entry(ids, **kw):
    e = cpd("x", **kw)
    e["all_identifiers"] = ids
    return e


def test_validate_flags_duplicate_identifier_across_roles():
    table = {"chebi:1": _entry(["chebi:1", "chebi:10004"])}
    roles = {"chebi.role:33229": {"parents": [], "all_identifiers": ["chebi.role:33229", "chebi:10004"]}}
    errs = validate_build(table=table, roles=roles, subtree_substrates={}, subtree_roles={},
                          class_trees={}, resolution=ResolutionResult(), chebi_missing_needed=set())
    assert any("chebi:10004" in e for e in errs)


def test_validate_flags_dangling_and_missing_assignment():
    table = {"kegg.compound:C1": dict(_entry(["kegg.compound:C1"]), stereo_form_of=["kegg.compound:C9"])}
    res = ResolutionResult(decisions={("p", "t", 0): Decision("kegg.compound:C1", "name_match", None)})
    errs = validate_build(table=table, roles={}, subtree_substrates={"2.A.1": ["kegg.compound:C7"]},
                          subtree_roles={}, class_trees={}, resolution=res, chebi_missing_needed=set())
    assert any("C9" in e for e in errs) and any("C7" in e for e in errs)
    assert any("stereo_assignment" in e for e in errs)


def test_write_all_atomic(tmp_path):
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    write_all_atomic({a: "1", b: "2"})
    assert a.read_text() == "1" and not list(tmp_path.glob("*.tmp"))


def test_write_resolved_tables(tmp_path):
    csv = tmp_path / "m.csv"
    pd.DataFrame({"compound": ["Alanine", "Mystery"], "v": ["1", "2"]}).to_csv(csv, index=False)
    t = PaperTable("p", "t", csv, 0, [PaperRow("p", "t", 0, "Alanine", None, None),
                                      PaperRow("p", "t", 1, "Mystery", None, None)])
    out = write_resolved_tables([t], {("p", "t", 0): Decision("kegg.compound:C00041", "dominant_stereo_form",
                                                               "dominant_form_assumed")})
    text = out[tmp_path / "m_resolved.csv"]
    df = pd.read_csv(__import__("io").StringIO(text), dtype=str, keep_default_na=False)
    assert list(df.columns)[-3:] == ["metabolite_id", "resolution_method", "stereo_assignment"]
    assert df.iloc[1].tolist()[-3:] == ["", "unresolved", ""]
    rep = json.loads(out[tmp_path / "m_resolution_report.json"])
    assert rep["resolved"] == 1 and rep["unresolved"] == 1


def test_repoint_reason_priority():
    ident = IdentityResult({"kegg.compound:C00133": "kegg.compound:C00133"}, {}, {}, {})
    row = PaperRow("p", "t", 0, "Alanine", None, None)
    r = repoint_reason(old_id="kegg.compound:C00133",
                       new=Decision("kegg.compound:C00041", "dominant_stereo_form", "dominant_form_assumed"),
                       row=row, identity=ident, rekeys={}, dropped={}, case_fold_keys=set())
    assert r == "dominant_stereo_form"
    assert repoint_reason(old_id="kegg.compound:C00041", new=Decision("kegg.compound:C00041", "name_match", "x"),
                          row=row, identity=ident, rekeys={}, dropped={}, case_fold_keys=set()) is None


def test_remap_rows_bases():
    table = {"kegg.compound:C00244": _entry(["kegg.compound:C00244", "chebi:14654"]),
             "chebi:17315": _entry(["chebi:17315", "chebi:12961"])}
    ident = IdentityResult({"chebi:14654": "kegg.compound:C00244", "kegg.compound:C00244": "kegg.compound:C00244"},
                           {"kegg.compound:C00244": ["chebi:14654", "kegg.compound:C00244"]},
                           {"kegg.compound:C00244": ["chebi_conjugate"]}, {"chebi:14654": "chebi_conjugate"})
    rows, unaccounted = remap_rows(
        old_ids={"kegg.compound:C00244", "chebi:14654", "chebi:10004", "mnx:MNXM1565340", "chebi:12961", "chebi:999"},
        table=table, identity=ident, rekeys={}, role_of={"chebi:10004": "chebi.role:33229"},
        dropped={"mnx:MNXM1565340": "kegg.compound:C00062"},
        normalize=lambda x: {"chebi:12961": "chebi:17315"}.get(x, x))
    by = {r["old_id"]: r for r in rows}
    assert by["chebi:14654"]["basis"] == "chebi_conjugate"
    assert by["chebi:10004"] == {"old_id": "chebi:10004", "new_id": "", "basis": "removed_role", "role_id": "chebi.role:33229"}
    assert by["mnx:MNXM1565340"]["basis"] == "removed_placeholder"
    assert by["chebi:12961"]["basis"] == "chebi_primary_id"
    assert unaccounted == ["chebi:999"]
```

- [ ] **Step 2: Run** → FAIL. **Step 3: Implement `outputs.py`.**

```python
"""Stage P: validation, atomic multi-file write, resolved tables, mapping v2, report rows."""
from __future__ import annotations

import io
import json
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable

import pandas as pd

REASON_ORDER = ("alias_override", "dominant_stereo_form", "case_fold_convergence",
                "direct_id_canonicalised", "placeholder_recurated", "remeasured_elsewhere", "merged")


def validate_build(*, table, roles, subtree_substrates, subtree_roles, class_trees, resolution,
                   chebi_missing_needed):
    errs = []
    errs += [f"ambiguous paper name (needs an override): {a}" for a in resolution.ambiguous]
    errs += [f"same-table collision: {c}" for c in resolution.collisions]
    errs += list(resolution.errors)
    owner: dict[str, str] = {}
    for kind, items in (("Metabolite", table), ("ChemicalRole", roles)):
        for node, e in items.items():
            for x in e.get("all_identifiers", []):
                prev = owner.setdefault(x, f"{kind} {node}")
                if prev != f"{kind} {node}":
                    errs.append(f"identifier {x} on two nodes: {prev} and {kind} {node}")
    iks = Counter(e["inchikey"] for e in table.values() if e.get("inchikey"))
    errs += [f"InChIKey {k} on {n} Metabolite nodes" for k, n in iks.items() if n > 1]
    for node, e in table.items():
        for k in ("stereo_form_of", "enantiomer_of"):
            errs += [f"{node}.{k} → {t} is not a written compound" for t in e.get(k, []) if t not in table]
    for tc, ids in subtree_substrates.items():
        errs += [f"subtree_substrates[{tc}] → {i} is not a written compound" for i in ids if i not in table]
    for tc, ids in subtree_roles.items():
        errs += [f"subtree_roles[{tc}] → {i} is not a ChemicalRole" for i in ids if i not in roles]
    for r, e in roles.items():
        errs += [f"ChemicalRole {r} parent {p} missing" for p in e.get("parents", []) if p not in roles]
    for code, tree in class_trees.items():
        ids = {n["id"] for n in tree["nodes"]}
        errs += [f"{code}: parent {n['parent']} missing" for n in tree["nodes"] if n["parent"] and n["parent"] not in ids]
        errs += [f"{code}: member {m} dangles" for m in tree["members"]
                 if m["metabolite_id"] not in table or m["class_id"] not in ids]
    for key, d in resolution.decisions.items():
        if d.metabolite_id and d.metabolite_id not in table:
            errs.append(f"row {key} resolved to {d.metabolite_id}, which step 6 did not write")
        if d.metabolite_id and not d.stereo_assignment:
            errs.append(f"row {key} resolved without a stereo_assignment")
    errs += [f"ChEBI id {c} needed by a canonical node is missing from the committed ChEBI cache"
             for c in sorted(chebi_missing_needed)]
    return errs


def write_all_atomic(files):
    tmps = []
    try:
        for path, text in files.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_name(path.name + ".tmp")
            tmp.write_text(text)
            tmps.append((tmp, path))
    except BaseException:
        for tmp, _ in tmps:
            tmp.unlink(missing_ok=True)
        raise
    for tmp, path in tmps:
        os.replace(tmp, path)


def write_resolved_tables(tables, decisions):
    out = {}
    for t in tables:
        df = pd.read_csv(t.csv_path, dtype=str, keep_default_na=False, skiprows=t.skip_rows or None)
        mid, meth, sa = [""] * len(df), ["unresolved"] * len(df), [""] * len(df)
        for row in t.rows:
            d = decisions.get((row.paper, row.table, row.row_index))
            if d and d.metabolite_id:
                mid[row.row_index], meth[row.row_index], sa[row.row_index] = (
                    d.metabolite_id, d.resolution_method, d.stereo_assignment)
        df["metabolite_id"], df["resolution_method"], df["stereo_assignment"] = mid, meth, sa
        resolved_path = t.csv_path.with_name(t.csv_path.stem + "_resolved.csv")
        buf = io.StringIO()
        df.to_csv(buf, index=False)
        out[resolved_path] = buf.getvalue()
        n_res = sum(1 for m in mid if m)
        report = {"source_csv": str(t.csv_path), "resolved_csv": str(resolved_path),
                  "total_rows": len(df), "resolved": n_res, "unresolved": len(df) - n_res,
                  "match_rate": round(n_res / len(df), 4) if len(df) else 0.0,
                  "method_counts": dict(sorted(Counter(meth).items()))}
        out[t.csv_path.with_name(t.csv_path.stem + "_resolution_report.json")] = json.dumps(report, indent=2, sort_keys=True)
    return out


def build_mapping_v2(*, id_aliases, resolution, overrides, tables, compounds_summary, per_paper_stats):
    decisions = defaultdict(lambda: defaultdict(dict))
    for t in tables:
        for row in t.rows:
            d = resolution.decisions.get((row.paper, row.table, row.row_index))
            if d:
                key = row.direct_id if row.direct_id and d.resolution_method.endswith("_direct") else row.name_key
                decisions[row.paper][row.table][key] = {"metabolite_id": d.metabolite_id,
                                                        "resolution_method": d.resolution_method,
                                                        "stereo_assignment": d.stereo_assignment}
    unresolved = defaultdict(lambda: defaultdict(list))
    for u in resolution.unresolved:
        unresolved[u["paper"]][u["table"]].append(u.get("name"))
    return {
        "schema_version": 2,
        "id_aliases": dict(sorted(id_aliases.items())),
        "name_lookup": dict(sorted(resolution.name_lookup.items())),
        "paper_overrides": {p: {k: {"id": o.id, "stereo_assignment": o.stereo_assignment, "reason": o.reason}
                                for k, o in sorted(ovs.items())} for p, ovs in sorted(overrides.items())},
        "decisions": json.loads(json.dumps(decisions, sort_keys=True)),
        "unresolved": json.loads(json.dumps(unresolved, sort_keys=True)),
        "ambiguous": {},
        "compounds": compounds_summary,
        "conflicts": {},
        "per_paper_stats": per_paper_stats,
    }


def repoint_reason(*, old_id, new, row, identity, rekeys, dropped, case_fold_keys):
    if not old_id or old_id == new.metabolite_id:
        return None
    if new.resolution_method == "alias_override":
        return "alias_override"
    if new.resolution_method == "dominant_stereo_form":
        return "dominant_stereo_form"
    if row.name_key in case_fold_keys and new.resolution_method == "name_match":
        return "case_fold_convergence"
    if new.resolution_method.endswith("_direct"):
        return "direct_id_canonicalised"
    if old_id.startswith("mnx:"):
        return "placeholder_recurated"
    if old_id in dropped:
        return "remeasured_elsewhere"
    if identity.canonical_of.get(rekeys.get(old_id, old_id)) == new.metabolite_id:
        return "merged"
    return "merged"


def remap_rows(*, old_ids, table, identity, rekeys, role_of, dropped, normalize):
    rows, unaccounted = [], []
    for old in sorted(old_ids):
        if old in table:
            continue
        norm = normalize(old)
        if old in role_of:
            rows.append({"old_id": old, "new_id": "", "basis": "removed_role", "role_id": role_of[old]})
        elif norm in dropped:
            basis = "removed_placeholder" if old.startswith("mnx:") else "remeasured_elsewhere"
            rows.append({"old_id": old, "new_id": dropped[norm], "basis": basis, "role_id": ""})
        elif norm in rekeys:
            new = identity.canonical_of.get(rekeys[norm], rekeys[norm])
            rows.append({"old_id": old, "new_id": new, "basis": "rekey_chebi_kegg_xref", "role_id": ""})
        elif norm in identity.member_basis:
            rows.append({"old_id": old, "new_id": identity.canonical_of[norm],
                         "basis": identity.member_basis[norm], "role_id": ""})
        elif norm in table and old.startswith("chebi:"):
            rows.append({"old_id": old, "new_id": norm, "basis": "chebi_primary_id", "role_id": ""})
        else:
            unaccounted.append(old)
    return rows, unaccounted


def build_identity_report(**sections):
    return json.loads(json.dumps(sections, sort_keys=True, default=sorted))
```

- [ ] **Step 4: Update `__init__.py`** to re-export the public names of `canonical`, `structure`, `crossrefs`, `stereo`, `paper`, `roles`, `classes`, `lumping`, `outputs` (one `from … import (…)  # noqa: F401` block per module).

- [ ] **Step 5: Run the whole Phase-1 suite.** `uv run pytest tests/test_chebi_utils.py tests/test_mi_*.py -q` → PASS.

- [ ] **Step 6: Commit.**

```bash
git add multiomics_kg/utils/metabolite_identity tests/test_mi_outputs.py
git commit -m "feat(metabolite): stage-P validation, atomic write, resolved tables, remap/repoint rows"
```

---

## Phase 2 — Step 6 stages A–P (step 7 folded in)

Step 6 keeps its I/O stages in `build_kegg_metabolism_xrefs.py` (A: ChEBI, B: existing TCDB rollup + prune, C: role partition, D: MNX resolution + ChEBI re-key + paper-candidate harvest + KEGG enrichment, N: TCDB bridges). Stages E–M, O and the validation part of P are the pure `run_identity_stages` (Task 12). `main()` writes files only after `run_identity_stages` returns (Task 13).

**Pre-existing writes that stay where they are** (not stage-P outputs; spec §5.2 lists the stage-P files): the gitignored raw downloads (`cache/data/{kegg,tcdb,chebi}/raw/`) and `tcdb_hierarchy.json`, which `_build_tcdb_hierarchy` rebuilds deterministically from the raw TSVs before stage B, as today.

### Task 11: Stages A, C, D in `build_kegg_metabolism_xrefs.py`

**Files:**
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py`
- Test: `tests/test_build_kegg_metabolism_xrefs.py` (extend), `tests/test_build_kegg_metabolism_xrefs_metabolomics.py` (replace: see Step 6)

**Interfaces:**
- Consumes: `chebi_utils` (Task 1), `roles.partition_role_strings` / `build_role_hierarchy` / `build_subtree_roles` (Task 6), `paper.harvest_paper_tables` / `PaperTable` / `Override` (Task 7), `metabolite_utils.resolve_metabolite_candidates` (Task 7).
- Produces (all module-level in `build_kegg_metabolism_xrefs.py`):
  - Path constants: `CHEBI_DIR = PROJECT_ROOT/"cache"/"data"/"chebi"`, `CHEBI_CACHE_FILE = CHEBI_DIR/"chebi_metabolite_xrefs.json"`, `METABOLOMICS_DIR = PROJECT_ROOT/"cache"/"data"/"metabolomics"`, `LOCK_FILE = METABOLOMICS_DIR/"dominant_stereo_forms.json"`, `IDENTITY_REPORT_FILE = METABOLOMICS_DIR/"metabolite_identity_report.json"`, `MAPPING_FILE = METABOLOMICS_DIR/"metabolite_id_mapping.json"`, `IDENTITY_OVERRIDES_FILE = PROJECT_ROOT/"config"/"metabolite_identity_overrides.yaml"`.
  - `_load_chebi(cache_root: Path, refetch_raw: bool) -> tuple[ChebiIndex, bool]` (stage A; `(index, from_raw)`).
  - `_gene_reachable_sets(raw)` additionally returns key `"rxn_genes": dict[str, set[str]]` (reaction → `"{data_dir}:{locus_tag}"` gene keys).
  - `_resolve_substrate_strings(substrate_strings, conn, chebi: ChebiIndex | None = None, rekeys: dict[str, str] | None = None)` — unchanged 2-tuple return; with `chebi`: secondary → primary before MNX lookup, and a ChEBI KEGG xref supplies the KEGG id when MNX has none (recorded in `rekeys` as `chebi:<primary n>` → `kegg.compound:C…`).
  - `_mnxm_primary(mnxm: str, conn, chebi: ChebiIndex, rekeys: dict[str, str]) -> str`.
  - `_normalize_paper_id(curie: str, conn, chebi: ChebiIndex, rekeys: dict[str, str]) -> str`.
  - `@dataclass PaperCandidates(name_candidates: dict[str, set[str]], ids: set[str], non_kegg_props: dict[str, dict])`.
  - `_harvest_paper_candidates(tables: list[PaperTable], overrides: dict[str, dict[str, Override]], lock: dict, conn, chebi: ChebiIndex, rekeys: dict[str, str]) -> PaperCandidates`.
  - `_add_paper_candidates(kegg_data: dict, non_kegg_props: dict[str, dict]) -> None` (inserts missing `additional_compounds` with `evidence_sources: []`; **never** tags `metabolomics` — stage K does).
  - `_attach_mnx_chebi_aliases(table: dict[str, dict], conn) -> None` (transient field `mnx_chebi_aliases`).
  - `_to_table(kegg_data: dict) -> dict[str, dict]` (keyed by full node id).
  - **Removed:** `_harvest_paper_metabolites`, `_fold_paper_metabolites_into_kegg_data`, `_write_metabolite_id_mapping` (their jobs move to `paper.py` / stage K / `outputs.build_mapping_v2`).

- [ ] **Step 1: Write the failing tests** (append to `tests/test_build_kegg_metabolism_xrefs.py`; reuse its `_resolve_subs_db(tmp_path)` helper, adding rows as needed):

```python
from tests.metabolite_identity_helpers import make_chebi, term


def _subs_conn_with_chebi_only(tmp_path):
    """MNX knows CHEBI:17234 (glucose) under MNXM41 but has NO kegg alias for it."""
    import sqlite3
    conn = sqlite3.connect(tmp_path / "r.db")
    conn.executescript("""
      CREATE TABLE compounds (mnxm_id TEXT, name TEXT, formula TEXT, mass REAL, inchikey TEXT, smiles TEXT);
      CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT);
      CREATE TABLE compound_names (mnxm_id TEXT, name_normalized TEXT);
      INSERT INTO compounds VALUES ('MNXM41','glucose','C6H12O6',180.0,'WQZGKKKJIJFFOK-UHFFFAOYSA-N',NULL);
      INSERT INTO compound_aliases VALUES ('chebi','4167','MNXM41');
      INSERT INTO compound_aliases VALUES ('chebi','10004','MNXM9');
      INSERT INTO compounds VALUES ('MNXM9','vitamin (role)',NULL,NULL,NULL,NULL);
    """)
    return conn


def test_resolve_substrate_strings_rekeys_via_chebi_kegg_xref(tmp_path):
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _resolve_substrate_strings
    chebi = make_chebi(CHEBI_4167=term("D-glucopyranose", alt=["CHEBI:12345"], kegg=["C00031"]))
    rekeys = {}
    s2p, props = _resolve_substrate_strings({"CHEBI:12345;glucose"}, _subs_conn_with_chebi_only(tmp_path),
                                            chebi=chebi, rekeys=rekeys)
    assert s2p == {"CHEBI:12345;glucose": "kegg.compound:C00031"}
    assert rekeys == {"chebi:4167": "kegg.compound:C00031"} and props == {}


def test_resolve_substrate_strings_without_chebi_is_unchanged(tmp_path):
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _resolve_substrate_strings
    s2p, _ = _resolve_substrate_strings({"CHEBI:4167;glucose"}, _subs_conn_with_chebi_only(tmp_path))
    assert s2p == {"CHEBI:4167;glucose": "chebi:4167"}


def test_gene_reachable_sets_returns_rxn_genes(tmp_path, monkeypatch, synthetic_kegg_raw):
    import multiomics_kg.download.build_kegg_metabolism_xrefs as bx
    strain = tmp_path / "s"
    _write_strain_annotations(strain, {"G1": {"kegg_reactions": ["R00200"]},
                                       "G2": {"kegg_reactions": ["R00200"]}})
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [{"data_dir": str(strain)}])
    sets = bx._gene_reachable_sets(synthetic_kegg_raw)
    assert sets["rxn_genes"]["R00200"] == {f"{strain}:G1", f"{strain}:G2"}
```

  `synthetic_kegg_raw` and `_write_strain_annotations` already exist in this file; if the fixture's shape is the parsed `raw` dict, pass it; if it is a path, call `bx._parse_raw_into_dict` on it first (read the fixture before writing the test).

  Paper-candidate harvest test (new file section in `tests/test_build_kegg_metabolism_xrefs_metabolomics.py`, replacing the four `_harvest_paper_metabolites` / fold / mapping tests, whose behaviour is now covered by `tests/test_mi_paper.py`):

```python
def test_harvest_paper_candidates_collects_names_ids_overrides_lock(tmp_path):
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _harvest_paper_candidates
    from multiomics_kg.utils.metabolite_identity.paper import Override, PaperRow, PaperTable
    from tests.metabolite_identity_helpers import make_chebi, term
    conn = _make_stub_resolver(tmp_path)
    chebi = make_chebi(CHEBI_17234=term("glucose", kegg=["C00031"]))
    t = PaperTable("p", "t", None, 0, [PaperRow("p", "t", 0, "glucose", None, None),
                                       PaperRow("p", "t", 1, "Tyrosine", "chebi:17234", "chebi_direct")])
    ov = {"p": {"pyruvic": Override("kegg.compound:C00022", "as_reported", "", "Pyruvic")}}
    lock = {"names": {"x": {"id": "kegg.compound:C00041"}}, "direct_ids": {}}
    rekeys = {}
    pc = _harvest_paper_candidates([t], ov, lock, conn, chebi, rekeys)
    assert pc.name_candidates["glucose"] == {"kegg.compound:C00031"}
    assert {"kegg.compound:C00031", "kegg.compound:C00022", "kegg.compound:C00041"} <= pc.ids
    assert rekeys == {"chebi:17234": "kegg.compound:C00031"}   # paper-path ChEBI candidate re-keyed
```

- [ ] **Step 2: Run** `uv run pytest tests/test_build_kegg_metabolism_xrefs.py tests/test_build_kegg_metabolism_xrefs_metabolomics.py -q` → the new tests FAIL.

- [ ] **Step 3: Implement stage A.**

```python
from multiomics_kg.utils.chebi_utils import (
    ChebiIndex, download_chebi_obo, index_from_cache, parse_chebi_obo)


def _load_chebi(cache_root: Path, refetch_raw: bool) -> tuple[ChebiIndex, bool]:
    """Stage A (spec §5.1, §6.4): rebuild from the raw obo whenever it is available; reuse the
    committed derived cache only when the download fails (TCDB outage precedent)."""
    raw_dir = cache_root / "chebi" / "raw"
    try:
        obo = download_chebi_obo(raw_dir, force=refetch_raw)
        release, terms = parse_chebi_obo(obo)
        log.info(f"  ChEBI release {release}: {len(terms)} terms (from raw obo)")
        return ChebiIndex(terms, release), True
    except Exception as exc:
        if not CHEBI_CACHE_FILE.exists():
            raise
        log.warning(f"  ChEBI download/parse failed ({exc}); reusing the committed "
                    f"{CHEBI_CACHE_FILE.name}. Ids it lacks are reported, and step 6 fails "
                    f"if a canonical node needs one.")
        return index_from_cache(json.loads(CHEBI_CACHE_FILE.read_text())), False
```

- [ ] **Step 4: Implement the stage-D changes.**

  4a. `_gene_reachable_sets`: inside the gene loop, `for gid, g in genes.items():` (the dict key is the locus tag), add

```python
            for rxn in g.get("kegg_reactions") or []:
                if isinstance(rxn, str) and rxn.startswith("R"):
                    rxns.add(rxn)
                    rxn_genes.setdefault(rxn, set()).add(f"{row['data_dir']}:{gid}")
```

  and return `"rxn_genes": rxn_genes` alongside the existing keys.

  4b. `_resolve_substrate_strings(substrate_strings, conn, chebi=None, rekeys=None)`. With `chebi` set:
  - Pass 1 records `primary_num = chebi.primary(chebi_value)` (bare number) per parsed string; the MNX lookup in pass 2 queries **both** the original value and the primary number (`IN (...)` over their union), and `chebi_to_mnxm[original]` takes the first hit for either.
  - `_primary_for(mnxm, primary_num)`: (1) MNX KEGG alias → `kegg.compound:…`; (2) else, if `chebi.term(primary_num)["kegg"]` is non-empty → `kegg.compound:<lowest C>` and `rekeys[f"chebi:{primary_num}"] = that`; (3) else `chebi:{primary_num}`; (4) else (no ChEBI term) the old `chebi:{mnx alias}` / `mnx:{mnxm}`.
  - Pass 5's `compound_props[primary]["chebi_id"]` = `primary_num` for `chebi:` primaries.
  With `chebi=None` every pass is byte-for-byte today's behaviour (the existing five `test_resolve_substrate_strings_*` tests must still pass unchanged).

  4c. `_mnxm_primary` and `_normalize_paper_id`:

```python
def _mnxm_primary(mnxm: str, conn: sqlite3.Connection, chebi: ChebiIndex, rekeys: dict[str, str]) -> str:
    row = conn.execute("SELECT value FROM compound_aliases WHERE mnxm_id=? AND source='kegg.compound' "
                       "ORDER BY value LIMIT 1", (mnxm,)).fetchone()
    if row:
        return f"kegg.compound:{row[0]}"
    prims = sorted({p for (v,) in conn.execute(
        "SELECT value FROM compound_aliases WHERE mnxm_id=? AND source='chebi'", (mnxm,))
        if (p := chebi.primary(v))}, key=lambda c: (chebi.charge(c) != 0, int(c.split(":")[1])))
    for p in prims:
        kegg = sorted((chebi.term(p) or {}).get("kegg", []))
        if kegg:
            rekeys[f"chebi:{p.split(':')[1]}"] = f"kegg.compound:{kegg[0]}"
            return f"kegg.compound:{kegg[0]}"
    if prims:
        return f"chebi:{prims[0].split(':')[1]}"
    return f"mnx:{mnxm}"


def _normalize_paper_id(curie: str, conn, chebi: ChebiIndex, rekeys: dict[str, str]) -> str:
    prefix, _, local = curie.partition(":")
    if prefix == "kegg.compound":
        return f"kegg.compound:{local.strip().upper()}"
    if prefix == "chebi":
        p = chebi.primary(local)
        if p is None:
            return f"chebi:{int(local)}" if local.isdigit() else curie
        kegg = sorted((chebi.term(p) or {}).get("kegg", []))
        if kegg:
            rekeys[f"chebi:{p.split(':')[1]}"] = f"kegg.compound:{kegg[0]}"
            return f"kegg.compound:{kegg[0]}"
        return f"chebi:{p.split(':')[1]}"
    if prefix == "mnx":
        return _mnxm_primary(local, conn, chebi, rekeys)
    return curie
```

  4d. `_harvest_paper_candidates`:

```python
@dataclass
class PaperCandidates:
    name_candidates: dict[str, set[str]]
    ids: set[str]
    non_kegg_props: dict[str, dict]


def _harvest_paper_candidates(tables, overrides, lock, conn, chebi, rekeys) -> PaperCandidates:
    """Stage D, paper side (spec §5.2 D): collect CANDIDATES, never decisions. Every paper
    row's MNX name candidates, direct id, the paper's override targets and the lock targets
    are all enriched like any compound; stage J decides, stage K drops the losers."""
    name_candidates: dict[str, set[str]] = defaultdict(set)
    ids: set[str] = set()
    for t in tables:
        for row in t.rows:
            for mnxm in resolve_metabolite_candidates(row.name, conn):
                name_candidates[row.name_key].add(_mnxm_primary(mnxm, conn, chebi, rekeys))
            if row.direct_id:
                ids.add(_normalize_paper_id(row.direct_id, conn, chebi, rekeys))
    for ovs in overrides.values():
        for ov in ovs.values():
            ids.add(_normalize_paper_id(ov.id, conn, chebi, rekeys))
    for section in ("names", "direct_ids"):
        for key, entry in lock.get(section, {}).items():
            ids.add(_normalize_paper_id(entry["id"], conn, chebi, rekeys))
            if section == "direct_ids":
                ids.add(_normalize_paper_id(key, conn, chebi, rekeys))
    for cands in name_candidates.values():
        ids |= cands
    non_kegg = {}
    for pid in sorted(i for i in ids if not i.startswith("kegg.compound:")):
        if pid.startswith("mnx:"):
            mnxm = pid[4:]
        else:
            row = conn.execute("SELECT mnxm_id FROM compound_aliases WHERE source='chebi' AND value=? "
                               "ORDER BY mnxm_id LIMIT 1", (pid.split(":")[1],)).fetchone()
            mnxm = row[0] if row else None
        props = conn.execute("SELECT name, formula, mass, inchikey FROM compounds WHERE mnxm_id=?",
                             (mnxm,)).fetchone() if mnxm else None
        name, formula, mass, inchikey = props if props else (None, None, None, None)
        term = chebi.term(pid) if pid.startswith("chebi:") else None
        non_kegg[pid] = {"name": (term or {}).get("name") or name, "formula": formula or None,
                         "mass": mass, "inchikey": inchikey or None, "mnxm_id": mnxm,
                         "chebi_id": pid.split(":")[1] if pid.startswith("chebi:") else None}
    return PaperCandidates(dict(name_candidates), ids, non_kegg)
```

  4e. `_add_paper_candidates`, `_attach_mnx_chebi_aliases`, `_to_table`:

```python
def _add_paper_candidates(kegg_data: dict, non_kegg_props: dict[str, dict]) -> None:
    additional = kegg_data.setdefault("additional_compounds", {})
    for pid, props in non_kegg_props.items():
        additional.setdefault(pid, {**props, "evidence_sources": []})


def _attach_mnx_chebi_aliases(table: dict[str, dict], conn: sqlite3.Connection) -> None:
    by_mnx: dict[str, list[str]] = defaultdict(list)
    for node, e in table.items():
        if e.get("mnxm_id"):
            by_mnx[e["mnxm_id"]].append(node)
    found: dict[str, set[str]] = defaultdict(set)
    for chunk in _batched(sorted(by_mnx)):
        q = ",".join("?" * len(chunk))
        for mnxm, value in conn.execute(
                f"SELECT mnxm_id, value FROM compound_aliases WHERE source='chebi' AND mnxm_id IN ({q})", chunk):
            found[mnxm].add(value)
    for mnxm, nodes in by_mnx.items():
        for n in nodes:
            table[n]["mnx_chebi_aliases"] = sorted(found.get(mnxm, ()), key=lambda v: (len(v), v))


def _to_table(kegg_data: dict) -> dict[str, dict]:
    table = {f"kegg.compound:{k}": dict(v) for k, v in kegg_data["compounds"].items()}
    table.update({k: dict(v) for k, v in kegg_data.get("additional_compounds", {}).items()})
    return table
```

  Delete `_harvest_paper_metabolites`, `_fold_paper_metabolites_into_kegg_data` and `_write_metabolite_id_mapping`, and the tests that exercised them (`test_harvest_paper_metabolites_*` ×4, `test_fold_paper_metabolites_unions_evidence_sources`, `test_write_metabolite_id_mapping_creates_three_tier_skeleton`). `main()` still calls them until Task 13; until then leave `main()` calling the new helpers only where Task 13 says, and accept that `main()` is mid-refactor. Do **not** run step 6 for real between Task 11 and Task 13.

- [ ] **Step 5: Run** `uv run pytest tests/test_build_kegg_metabolism_xrefs.py tests/test_build_kegg_metabolism_xrefs_metabolomics.py -q` → PASS.

- [ ] **Step 6: Commit.**

```bash
git add multiomics_kg/download/build_kegg_metabolism_xrefs.py tests/test_build_kegg_metabolism_xrefs.py tests/test_build_kegg_metabolism_xrefs_metabolomics.py
git commit -m "feat(step6): ChEBI stage A, secondary→primary + KEGG re-key, paper candidate harvest"
```

### Task 12: The pure identity pipeline (`pipeline.py`, stages E–M, O, P-validation)

**Files:**
- Create: `multiomics_kg/utils/metabolite_identity/pipeline.py`, `tests/test_step6_identity_pipeline.py`

**Interfaces:**
- Consumes: everything from Tasks 1–10.
- Produces:
  - `class Step6ValidationError(RuntimeError)` with `.errors: list[str]`.
  - `@dataclass Step6Inputs` — fields: `table: dict[str, dict]`; `kegg_rest: dict` (`kos`, `pathways`, `subcategories`, `categories`, `reactions`); `subtree_substrates: dict[str, list[str]]`; `chemical_roles: dict[str, dict]`; `subtree_roles: dict[str, list[str]]`; `rekeys: dict[str, str]`; `role_of_old: dict[str, str]`; `tables: list[PaperTable]`; `paper_overrides: dict[str, dict[str, Override]]`; `name_candidates: dict[str, set[str]]`; `lock: dict`; `lumping_pin: dict | None`; `identity_overrides: IdentityOverrides`; `chebi: ChebiIndex`; `chebi_from_raw: bool`; `rxn_genes: dict[str, set[str]]`; `compound_to_pathways: dict[str, list[str]]`; `base_pathways: set[str]` (catalysis ∪ transport); `brite_raw: dict[str, dict]`; `brite_committed: dict[str, dict]`; `tcdb_kept: set[str]`; `tcdb_levels: dict[str, int]`; `tcdb_static: dict` (`kept_tcdb_ids`, `pfam_bridge`, `go_bridge`, `seed_aliases`); `previous_ids: set[str]`; `old_resolved: dict[tuple[str, str, int], str]`; `today: str`.
  - `@dataclass Step6Outputs(kegg_data: dict, tcdb_pruned: dict, class_trees: dict[str, dict], chebi_cache: dict | None, mapping: dict, lock: dict, report: dict, resolved_files: dict[Path, str])`.
  - `run_identity_stages(inp: Step6Inputs) -> Step6Outputs` (raises `Step6ValidationError`).
  - `@dataclass Step6Paths(kegg_data: Path, tcdb_pruned: Path, brite_dir: Path, chebi_cache: Path, mapping: Path, lock: Path, report: Path)`.
  - `render_outputs(out: Step6Outputs, paths: Step6Paths) -> dict[Path, str]` (JSON `indent=2, sort_keys=True` + trailing newline; BRITE files at `brite_dir / f"brite_{code}.json"`; the ChEBI cache only when `out.chebi_cache` is not None).

- [ ] **Step 1: Write the failing tests** (`tests/test_step6_identity_pipeline.py`). Build a synthetic `Step6Inputs` with: KEGG nitrate `C00244` + `chebi:17632` (merge by conjugate); the Ala family (`C00041` 2 genes, `C00133` 1 gene via reactions, `C01401` 0) linked by ChEBI `is_a` / enantiomer with InChIKeys sharing a first block; a paper table (tmp CSV) with rows "Alanine" and "D-Tryptophan"-style measurement-only candidate `kegg.compound:C00525` whose old resolved id was `C00525` and that loses to `C00078`; one role substrate already partitioned (`chemical_roles` with `chebi.role:33229`); a one-tree `brite_raw` for `br08001` containing `C00041`; `tcdb_kept = {"2.A.1", "3.A.3"}`.

```python
import copy
from pathlib import Path

import pytest

from multiomics_kg.utils.metabolite_identity.pipeline import (
    Step6Paths, Step6ValidationError, render_outputs, run_identity_stages)
from tests.step6_fixture import make_inputs     # helper module created in Step 3 below


def _render(tmp_path, inp):
    out = run_identity_stages(inp)
    paths = Step6Paths(tmp_path / "kegg_data.json", tmp_path / "tcdb_pruned.json", tmp_path,
                       tmp_path / "chebi.json", tmp_path / "map.json", tmp_path / "lock.json",
                       tmp_path / "report.json")
    return out, render_outputs(out, paths)


def test_canonical_graph_shape(tmp_path):
    out, _ = _render(tmp_path, make_inputs(tmp_path))
    kd = out.kegg_data
    assert "C00244" in kd["compounds"] and "chebi:17632" not in kd["additional_compounds"]
    assert kd["compounds"]["C00244"]["identity_basis"] == ["chebi_conjugate"]
    assert "metabolomics" in kd["compounds"]["C00041"]["evidence_sources"]
    assert "C00525" not in kd["compounds"]                                    # measurement-only, dropped
    assert kd["compounds"]["C00078"]["recurated_from_ids"] == ["kegg.compound:C00525"]
    assert "kegg.compound:C00525" not in kd["compounds"]["C00078"]["all_identifiers"]
    assert kd["compounds"]["C00041"]["stereo_form_of"] == ["kegg.compound:C01401"]
    assert out.tcdb_pruned["lumping_pin"]["lumping_tcdb_ids"]                 # seeded
    assert out.class_trees["br08001"]["members"] == [
        {"metabolite_id": "kegg.compound:C00041", "class_id": "kegg.brite:br08001.A1.B1"}]


def test_metabolomics_tag_only_on_final_targets(tmp_path):
    out, _ = _render(tmp_path, make_inputs(tmp_path))
    tagged = {k for k, v in out.kegg_data["compounds"].items() if "metabolomics" in v["evidence_sources"]}
    assert tagged == {"C00041", "C00078"}


def test_rerun_with_committed_lock_and_pin_is_byte_identical(tmp_path):
    inp = make_inputs(tmp_path)
    out1, files1 = _render(tmp_path, copy.deepcopy(inp))
    inp2 = make_inputs(tmp_path)
    inp2.lock, inp2.lumping_pin, inp2.today = out1.lock, out1.tcdb_pruned["lumping_pin"], "2099-01-01"
    _, files2 = _render(tmp_path, inp2)
    assert files1 == files2                                                   # frozen_on unchanged


def test_ambiguous_name_fails_and_renders_nothing(tmp_path):
    inp = make_inputs(tmp_path, extra_rows=[("Mystery", None)],
                      extra_candidates={"mystery": {"kegg.compound:C00244", "kegg.compound:C00041"}})
    with pytest.raises(Step6ValidationError) as e:
        run_identity_stages(inp)
    assert any("ambiguous" in m for m in e.value.errors)


def test_resolved_csv_targets_are_written_compounds(tmp_path):
    out, _ = _render(tmp_path, make_inputs(tmp_path))
    import io, pandas as pd
    for path, text in out.resolved_files.items():
        if path.name.endswith("_resolved.csv"):
            df = pd.read_csv(io.StringIO(text), dtype=str, keep_default_na=False)
            for mid in df["metabolite_id"]:
                assert not mid or mid.split(":", 1)[1] in out.kegg_data["compounds"] \
                    or mid in out.kegg_data["additional_compounds"]
```

- [ ] **Step 2: Run** `uv run pytest tests/test_step6_identity_pipeline.py -q` → FAIL.

- [ ] **Step 3: Write the fixture builder** `tests/step6_fixture.py` exporting `make_inputs(tmp_path, extra_rows=(), extra_candidates=None) -> Step6Inputs` for the scenario in Step 1. Use `tests.metabolite_identity_helpers.make_chebi/term/cpd`; reactions `{"R1": {"compounds": ["C00041", "C00244"], "pathways": []}, "R2": {"compounds": ["C00133"], "pathways": []}, "R3": {"compounds": ["C00078"], "pathways": []}}`, `rxn_genes={"R1": {"s:G1", "s:G2"}, "R2": {"s:G3"}, "R3": {"s:G4"}}`; `old_resolved` maps row 1 (`"D-Tryptophan"`-named row deliberately unqualified, name "Tryptophan") to `kegg.compound:C00525`; `name_candidates={"alanine": {C00133, C00041}, "tryptophan": {C00525, C00078}}` with a Trp family (`C00078` L-Trp, `C00525` D-Trp, enantiomers, InChIKey first block shared) and `catalyst` counts L 1 vs D 0; `previous_ids` = the table's ids plus `"chebi:10004"`; `role_of_old={"chebi:10004": "chebi.role:33229"}`; `kegg_rest={"kos": {}, "pathways": {}, "subcategories": {}, "categories": {}, "reactions": reactions}`; `brite_raw={"br08001": {"name": "br08001", "children": [{"name": "Peptides", "children": [{"name": "Amino acids", "children": [{"name": "C00041  L-Alanine"}]}]}]}}` and `brite_committed={}` with every other tree code mapped to `{"name": code, "children": [{"name": "X", "children": [{"name": "C99999  none"}]}]}` (a tree whose only member isn't ours → empty pruned tree); `tcdb_levels={"2.A.1": 2, "3.A.3": 2}`; `tcdb_static={"kept_tcdb_ids": ["2.A.1", "3.A.3"], "pfam_bridge": {}, "go_bridge": {}, "seed_aliases": {}}`; `today="2026-10-15"`. `extra_rows` appends rows to the paper table; `extra_candidates` merges into `name_candidates`.

- [ ] **Step 4: Implement `pipeline.py`.**

```python
"""Stages E-M, O and the stage-P validation of prepare-data step 6 as one pure function."""
from __future__ import annotations

import copy
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from multiomics_kg.utils.chebi_utils import ChebiIndex, build_chebi_cache
from multiomics_kg.utils.metabolite_identity.canonical import (
    IdentityOverrides, canonicalise, collect_merge_links, entity_terms, rewrite_to_canonical)
from multiomics_kg.utils.metabolite_identity.classes import (
    COMPOUND_BRITE_TREES, parse_compound_brite, prune_compound_brite)
from multiomics_kg.utils.metabolite_identity.crossrefs import build_identity_fields, pick_chebi_id
from multiomics_kg.utils.metabolite_identity.lumping import update_lumping_pin
from multiomics_kg.utils.metabolite_identity.outputs import (
    build_identity_report, build_mapping_v2, remap_rows, repoint_reason, validate_build,
    write_resolved_tables)
from multiomics_kg.utils.metabolite_identity.paper import (
    Override, PaperTable, StereoContext, resolve_paper_rows)
from multiomics_kg.utils.metabolite_identity.stereo import (
    finalize_stereo, is_anomer_pair, stereo_families, stereo_links, unspecified_members)
from multiomics_kg.utils.metabolite_identity.structure import apply_structure

_KEGG = "kegg.compound:"
_TRANSIENT = ("mnx_chebi_aliases",)


class Step6ValidationError(RuntimeError):
    def __init__(self, errors: list[str]):
        super().__init__(f"step 6 refused to write ({len(errors)} problem(s)):\n  " + "\n  ".join(errors))
        self.errors = errors


@dataclass
class Step6Inputs:
    table: dict
    kegg_rest: dict
    subtree_substrates: dict
    chemical_roles: dict
    subtree_roles: dict
    rekeys: dict
    role_of_old: dict
    tables: list
    paper_overrides: dict
    name_candidates: dict
    lock: dict
    lumping_pin: dict | None
    identity_overrides: IdentityOverrides
    chebi: ChebiIndex
    chebi_from_raw: bool
    rxn_genes: dict
    compound_to_pathways: dict
    base_pathways: set
    brite_raw: dict
    brite_committed: dict
    tcdb_kept: set
    tcdb_levels: dict
    tcdb_static: dict
    previous_ids: set
    old_resolved: dict
    today: str


@dataclass
class Step6Outputs:
    kegg_data: dict
    tcdb_pruned: dict
    class_trees: dict
    chebi_cache: dict | None
    mapping: dict
    lock: dict
    report: dict
    resolved_files: dict


@dataclass
class Step6Paths:
    kegg_data: Path
    tcdb_pruned: Path
    brite_dir: Path
    chebi_cache: Path
    mapping: Path
    lock: Path
    report: Path


def _catalyst_gene_counts(reactions, rxn_genes):
    genes = defaultdict(set)
    for rid, r in reactions.items():
        for c in r.get("compounds", []):
            genes[f"{_KEGG}{c}"] |= rxn_genes.get(rid, set())
    return {n: len(g) for n, g in genes.items()}


def _normalizer(chebi: ChebiIndex):
    def norm(x: str) -> str:
        if x.startswith("chebi:"):
            p = chebi.primary(x)
            return f"chebi:{p.split(':')[1]}" if p else x
        return x
    return norm


def _measured_set(table, res, old_resolved, norm):
    """Stage K (spec §3.6a): the metabolomics tag only on final targets; drop candidates
    with no catalysis / transport / measurement reason to exist."""
    targets = {d.metabolite_id for d in res.decisions.values() if d.metabolite_id}
    for n, e in table.items():
        ev = set(e.get("evidence_sources") or []) - {"metabolomics"}
        if n in targets:
            ev.add("metabolomics")
        e["evidence_sources"] = sorted(ev)
    moved: dict[str, Counter] = defaultdict(Counter)
    for key, old in old_resolved.items():
        d = res.decisions.get(key)
        if old and d and d.metabolite_id:
            moved[norm(old)][d.metabolite_id] += 1
    dropped = {}
    for n in sorted(n for n, e in table.items() if not e["evidence_sources"]):
        if moved.get(n):
            tgt = sorted(moved[n].items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
            dropped[n] = tgt
            table[tgt]["recurated_from_ids"] = sorted(set(table[tgt].get("recurated_from_ids") or []) | {n})
        del table[n]
    return targets, dropped


def _reprune_pathways(kegg_rest, table, targets, inp):
    allowed = set(inp.base_pathways)
    for t in targets:
        if t.startswith(_KEGG):
            allowed |= set(inp.compound_to_pathways.get(t[len(_KEGG):], []))
    pathways = {p: v for p, v in kegg_rest["pathways"].items() if p in allowed}
    for e in table.values():
        if e.get("pathways"):
            e["pathways"] = sorted(p for p in e["pathways"] if p in pathways)
    subs = {v.get("subcategory") for v in pathways.values() if v.get("subcategory")}
    subcategories = {s: v for s, v in kegg_rest["subcategories"].items() if s in subs}
    cats = {v.get("category") for v in subcategories.values() if v.get("category")}
    categories = {c: v for c, v in kegg_rest["categories"].items() if c in cats}
    return pathways, subcategories, categories


def _class_trees(inp, id_index, survivors):
    out = {}
    for code in COMPOUND_BRITE_TREES:
        raw = inp.brite_raw.get(code)
        if raw:
            nodes, members = parse_compound_brite(code, raw)
        elif code in inp.brite_committed:
            committed = inp.brite_committed[code]
            nodes = {n["id"]: {"parent": n["parent"], "name": n["name"], "level": n["level"]}
                     for n in committed["nodes"]}
            members = [(m["metabolite_id"][len(_KEGG):], m["class_id"])
                       if m["metabolite_id"].startswith(_KEGG) else (m["metabolite_id"], m["class_id"])
                       for m in committed["members"]]
        else:
            raise Step6ValidationError([f"no KEGG BRITE {code}: download failed and no committed file"])
        out[code] = prune_compound_brite(code, nodes, members, canonical_of=id_index, written=survivors)
    return out


def run_identity_stages(inp: Step6Inputs) -> Step6Outputs:
    chebi = inp.chebi
    norm = _normalizer(chebi)
    table = copy.deepcopy(inp.table)
    old_entries = copy.deepcopy(table)
    kegg_rest = copy.deepcopy(inp.kegg_rest)
    reactions = kegg_rest["reactions"]
    subs = {k: list(v) for k, v in inp.subtree_substrates.items()}

    # E — canonicalise
    reaction_refs = Counter(f"{_KEGG}{c}" for r in reactions.values() for c in r.get("compounds", []))
    pathway_refs = {n: len(e.get("pathways") or []) for n, e in table.items()}
    pins = []
    for t in inp.tables:
        for row in t.rows:
            ov = inp.paper_overrides.get(t.paper, {}).get(row.name_key)
            fixed = (ov.id if ov else None) or row.direct_id
            if fixed and fixed in table:
                pins.append((t.key, row.name_key, fixed))
    identity = canonicalise(table, collect_merge_links(table, chebi), overrides=inp.identity_overrides,
                            reaction_refs=reaction_refs, pathway_refs=pathway_refs, chebi=chebi,
                            table_pins=pins)
    # F — rewrite every id to canonical
    table = rewrite_to_canonical(table, reactions, subs, identity.canonical_of)
    # G — structure
    chebi_of = {}
    for n, e in table.items():
        cid = pick_chebi_id(n, e, chebi)
        chebi_of[n] = f"CHEBI:{cid}" if cid else None
    structure = apply_structure(table, chebi, chebi_of)
    # H — cross-reference fields
    identifier_index = build_identity_fields(table, identity, inp.rekeys, chebi, old_entries)
    id_index = dict(identifier_index)
    for old, canon in identity.canonical_of.items():
        id_index.setdefault(old, canon)
    for old, new in inp.rekeys.items():
        id_index.setdefault(old, identity.canonical_of.get(new, new))
    # I — stereo families
    terms_of = {n: entity_terms(n, e, chebi) for n, e in table.items()}
    stereo_pairs, enan_pairs = stereo_links(table, terms_of, chebi)
    families = stereo_families(set(table), stereo_pairs, enan_pairs)
    # J — resolve every paper row
    ctx = StereoContext(families, unspecified_members(stereo_pairs),
                        lambda a, b: is_anomer_pair(a, b, table, terms_of, chebi),
                        _catalyst_gene_counts(reactions, inp.rxn_genes))
    name_candidates = {k: {id_index.get(norm(c), norm(c)) for c in v} for k, v in inp.name_candidates.items()}
    res = resolve_paper_rows(inp.tables, overrides=inp.paper_overrides, name_candidates=name_candidates,
                             id_index=id_index, ctx=ctx, lock=inp.lock, today=inp.today)
    # K — measured set; drop measurement-only losers
    targets, dropped = _measured_set(table, res, inp.old_resolved, norm)
    pathways, subcategories, categories = _reprune_pathways(kegg_rest, table, targets, inp)
    survivors = set(table)
    # L — stereo lists on survivors only
    stereo_counts = finalize_stereo(table, stereo_pairs, enan_pairs, survivors)
    # M — compound classes
    class_trees = _class_trees(inp, id_index, survivors)
    # O — lumping pin
    pin, lumping = update_lumping_pin(inp.lumping_pin, inp.tcdb_kept, inp.tcdb_levels,
                                      {tc: len(ids) for tc, ids in subs.items()}, inp.today)
    # P — validate
    missing = set()
    if not inp.chebi_from_raw:
        missing = {n for n in table if n.startswith("chebi:") and chebi.primary(n) is None}
    errors = validate_build(table=table, roles=inp.chemical_roles, subtree_substrates=subs,
                            subtree_roles=inp.subtree_roles, class_trees=class_trees, resolution=res,
                            chebi_missing_needed=missing)
    if errors:
        raise Step6ValidationError(errors)

    compounds, additional = {}, {}
    for n in sorted(table):
        e = {k: v for k, v in table[n].items() if k not in _TRANSIENT and v is not None and v != []}
        (compounds if n.startswith(_KEGG) else additional)[n[len(_KEGG):] if n.startswith(_KEGG) else n] = e
    kegg_data = dict(kegg_rest, pathways=pathways, subcategories=subcategories, categories=categories,
                     compounds=compounds, additional_compounds=additional)
    tcdb_pruned = dict(inp.tcdb_static, subtree_substrates=dict(sorted(subs.items())),
                       chemical_roles=inp.chemical_roles, subtree_roles=inp.subtree_roles, lumping_pin=pin)

    raw_names = defaultdict(set)
    for t in inp.tables:
        for row in t.rows:
            raw_names[row.name_key].add(row.name)
    case_fold_keys = {k for k, v in raw_names.items() if len(v) > 1}
    repoints = []
    for t in inp.tables:
        for row in t.rows:
            key = (row.paper, row.table, row.row_index)
            d, old = res.decisions.get(key), inp.old_resolved.get(key)
            if d is None:
                continue
            reason = repoint_reason(old_id=old and norm(old), new=d, row=row, identity=identity,
                                    rekeys=inp.rekeys, dropped=dropped, case_fold_keys=case_fold_keys)
            if reason:
                repoints.append({"paper": row.paper, "table": row.table, "paper_name": row.name,
                                 "old_id": old, "new_id": d.metabolite_id,
                                 "stereo_assignment": d.stereo_assignment, "reason": reason})
    remap, unaccounted = remap_rows(old_ids=inp.previous_ids, table=table, identity=identity,
                                    rekeys=inp.rekeys, role_of=inp.role_of_old, dropped=dropped,
                                    normalize=norm)
    compounds_summary = defaultdict(lambda: {"aliases": set(), "methods": set()})
    for t in inp.tables:
        for row in t.rows:
            d = res.decisions.get((row.paper, row.table, row.row_index))
            if d and d.metabolite_id:
                compounds_summary[d.metabolite_id]["aliases"].add(row.name)
                compounds_summary[d.metabolite_id]["methods"].add(d.resolution_method)
    compounds_summary = {k: {"aliases": sorted(v["aliases"]), "methods": sorted(v["methods"])}
                         for k, v in sorted(compounds_summary.items())}
    per_paper = defaultdict(lambda: {"resolved": 0, "unresolved": 0, "total": 0})
    for t in inp.tables:
        for row in t.rows:
            s = per_paper[row.paper]
            s["total"] += 1
            s["resolved" if (row.paper, row.table, row.row_index) in res.decisions else "unresolved"] += 1
    mapping = build_mapping_v2(id_aliases={k: v for k, v in id_index.items() if k != v and v in table},
                               resolution=res,
                               overrides=inp.paper_overrides, tables=inp.tables,
                               compounds_summary=compounds_summary, per_paper_stats=dict(per_paper))
    report = build_identity_report(
        chebi_release=chebi.release,
        merge_groups={c: {"members": m, "basis": identity.basis.get(c, [])} for c, m in identity.members.items()},
        rekeys=inp.rekeys, rejected_links=[[l.a, l.b, l.basis] for l in identity.rejected],
        refused_assay_collisions=identity.refused_collisions,
        stale_identity_override=identity.stale_overrides, stale_alias_override=res.stale_overrides,
        stale_lock_entry=res.stale_lock, lock_target_remapped=res.lock_remapped,
        roles_removed=inp.role_of_old, dropped=dropped,
        dropped_same_entity_review=[{"dropped": d, "target": t,
                                     "dropped_name": old_entries.get(d, {}).get("name"),
                                     "target_name": table[t].get("name")} for d, t in sorted(dropped.items())],
        heuristic_name_candidates=_exact_name_pairs(table, identity),
        ambiguous_paper_names=res.ambiguous, unresolved_paper_names=res.unresolved,
        dominant_form_drift=res.drift, lumping=lumping, structure=structure, stereo_counts=stereo_counts,
        remap=remap, unaccounted_old_ids=unaccounted, repoints=repoints,
        counts={"metabolites": len(table), "kegg": len(compounds), "additional": len(additional),
                "merged_away": sum(len(m) - 1 for m in identity.members.values()),
                "rekeyed": len(inp.rekeys), "roles": len(inp.chemical_roles), "dropped": len(dropped),
                "class_nodes": sum(len(t["nodes"]) for t in class_trees.values()),
                "class_members": sum(len(t["members"]) for t in class_trees.values())})
    return Step6Outputs(kegg_data, tcdb_pruned, class_trees,
                        build_chebi_cache(chebi) if inp.chebi_from_raw else None,
                        mapping, res.lock, report, write_resolved_tables(inp.tables, res.decisions))


def _exact_name_pairs(table, identity):
    """Spec §1.2 'exact-name match not covered above': report only, never applied."""
    by_name = defaultdict(list)
    for n, e in table.items():
        if e.get("name"):
            by_name[e["name"].casefold()].append(n)
    return [sorted(v) for k, v in sorted(by_name.items()) if len(v) > 1]


def _dump(obj) -> str:
    return json.dumps(obj, indent=2, sort_keys=True) + "\n"


def render_outputs(out: Step6Outputs, paths: Step6Paths) -> dict[Path, str]:
    files = {paths.kegg_data: _dump(out.kegg_data), paths.tcdb_pruned: _dump(out.tcdb_pruned),
             paths.mapping: _dump(out.mapping), paths.lock: _dump(out.lock), paths.report: _dump(out.report)}
    for code, tree in out.class_trees.items():
        files[paths.brite_dir / f"brite_{code}.json"] = _dump(tree)
    if out.chebi_cache is not None:
        files[paths.chebi_cache] = _dump(out.chebi_cache)
    files.update(out.resolved_files)
    return files
```

  Two determinism rules the rerun test enforces: every list in every output is sorted (or ordered by a documented rule), and `today` reaches only newly frozen lock entries and a newly seeded pin.

- [ ] **Step 5: Run** `uv run pytest tests/test_step6_identity_pipeline.py tests/test_mi_*.py -q` → PASS.

- [ ] **Step 6: Commit.**

```bash
git add multiomics_kg/utils/metabolite_identity/pipeline.py tests/test_step6_identity_pipeline.py tests/step6_fixture.py
git commit -m "feat(step6): pure identity pipeline (stages E-M, O, P-validation) + render"
```

### Task 13: Wire `main()`, write once at stage P, retire step 7

**Files:**
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py` (`main`), `multiomics_kg/download/resolve_paper_metabolites.py`, `scripts/prepare_data.sh`
- Delete: `tests/test_resolve_paper_metabolites.py` (after porting, Step 1)
- Test: `tests/test_mi_outputs.py` (ported cases), `tests/test_prepare_data_step7_alias.py` (new)

**Interfaces:**
- Consumes: Tasks 11–12.
- Produces: `main(force: bool = False, refetch_raw: bool = False, today: str | None = None) -> None`; step 7 module prints `"step 7 is folded into step 6 (build_kegg_metabolism_xrefs); nothing to do"` and returns 0.

- [ ] **Step 1: Port step 7's tests.** `tests/test_resolve_paper_metabolites.py` covers four behaviours; each maps to `write_resolved_tables` / `harvest_paper_tables`. Add to `tests/test_mi_outputs.py`:
  - `test_write_resolved_tables_honours_skip_rows`: a CSV with 2 preamble lines, `skip_rows=2`; assert the resolved text has no preamble and row 0 is the first data row (step 7's `test_resolve_paper_metabolites_skip_rows`).
  - `test_harvest_id_col_with_prefix`: a row whose `id_col` cell is `kegg.compound:C00031` keeps its prefix; a bare `C00031` gets `id_type` prepended; `"0"` / `"NA"` are treated as no id (step 7's `test_resolve_paper_metabolites_id_col_with_prefix`).
  - `test_resolve_paper_metabolites_ambiguous_picks_first_sorted` is **intentionally not ported**: that behaviour is retired (spec §3.6 "an ambiguous match is never applied"); `test_ambiguous_name_fails_and_renders_nothing` (Task 12) replaces it.
  - `test_resolve_paper_metabolites_writes_resolved_csv` → already covered by `test_write_resolved_tables`.
  Then `git rm tests/test_resolve_paper_metabolites.py`.

- [ ] **Step 2: Write the step-7 alias test** (`tests/test_prepare_data_step7_alias.py`):

```python
import subprocess
import sys


def test_step7_module_is_a_noop():
    r = subprocess.run([sys.executable, "-m", "multiomics_kg.download.resolve_paper_metabolites", "--force"],
                       capture_output=True, text=True)
    assert r.returncode == 0 and "folded into step 6" in r.stdout


def test_prepare_data_default_steps_drop_7():
    text = open("scripts/prepare_data.sh").read()
    assert 'STEPS="0 9 1 2 3 4 5 6 8 10"' in text
    assert 'STEPS="9 1 2 3 4 5 6 8 10"' in text
    assert "folded into step 6" in text
```

- [ ] **Step 3: Rewrite `main()`.** Keep everything up to and including the TCDB prune and the `subtree_strings` / `all_substrate_strings` computation (current lines 1276–1354). Replace from `log.info(f"Resolving {len(all_substrate_strings)} …")` to the end of `main()` with:

```python
        # A — ChEBI
        chebi, chebi_from_raw = _load_chebi(cache_root, refetch_raw)
        today = today or date.today().isoformat()

        # previous outputs: old ids for the remap, the frozen lock and pin, old paper resolutions
        prev_kegg = json.loads(KEGG_DATA_FILE.read_text()) if KEGG_DATA_FILE.exists() else {}
        previous_ids = ({f"kegg.compound:{k}" for k in prev_kegg.get("compounds", {})}
                        | set(prev_kegg.get("additional_compounds", {})))
        prev_pruned_path = cache_root / "tcdb" / "tcdb_pruned.json"
        prev_pruned = json.loads(prev_pruned_path.read_text()) if prev_pruned_path.exists() else {}
        lock = json.loads(LOCK_FILE.read_text()) if LOCK_FILE.exists() else {"names": {}, "direct_ids": {}}

        # C — roles out of the substrate strings, before MNX
        molecule_strings, role_map = partition_role_strings(all_substrate_strings, chebi)
        chemical_roles = build_role_hierarchy(set(role_map.values()), chebi)
        subtree_roles = build_subtree_roles(subtree_strings, role_map)
        role_of_old = {}
        for old in prev_kegg.get("additional_compounds", {}):
            if old.startswith("chebi:") and chebi.is_role(old):
                role_of_old[old] = role_node_id(chebi.primary(old))
        log.info(f"  roles: {len(set(role_map.values()))} TCDB roles → {len(chemical_roles)} ChemicalRole nodes")

        # D — resolve molecules (secondary→primary, ChEBI KEGG re-key), harvest paper candidates
        rekeys: dict[str, str] = {}
        string_to_primary, compound_props = _resolve_substrate_strings(
            molecule_strings, conn, chebi=chebi, rekeys=rekeys)
        subtree_substrates = {}
        for tc, strings in subtree_strings.items():
            primaries = {string_to_primary[s] for s in strings if s in string_to_primary}
            if primaries:
                subtree_substrates[tc] = sorted(primaries)
        distinct_substrate_primaries = {p for ids in subtree_substrates.values() for p in ids}
        substrate_kegg_cpds = {p[len("kegg.compound:"):] for p in distinct_substrate_primaries
                               if p.startswith("kegg.compound:")}
        cpd_to_pw = raw.get("compound_to_pathways", {})
        transport_pws = {p for cpd in substrate_kegg_cpds for p in cpd_to_pw.get(cpd, [])}

        tables, paper_overrides = harvest_paper_tables(load_all_paperconfigs(), PROJECT_ROOT)
        paper = _harvest_paper_candidates(tables, paper_overrides, lock, conn, chebi, rekeys)
        paper_kegg_cpds = {i[len("kegg.compound:"):] for i in paper.ids if i.startswith("kegg.compound:")}
        paper_pws = {p for cpd in paper_kegg_cpds for p in cpd_to_pw.get(cpd, [])}

        extended_sets = dict(catalysis_sets)
        extended_sets["cpds"] = catalysis_cpds | substrate_kegg_cpds | paper_kegg_cpds
        extended_sets["pws"] = catalysis_pws | transport_pws | paper_pws
        kegg_data = build_pruned_kegg_data(raw, conn, sets=extended_sets)
        _fold_substrates_into_kegg_data(kegg_data, catalysis_cpds=catalysis_cpds,
                                        substrate_kegg_cpds=substrate_kegg_cpds,
                                        substrate_primaries=distinct_substrate_primaries,
                                        compound_props=compound_props)
        _add_paper_candidates(kegg_data, paper.non_kegg_props)
        table = _to_table(kegg_data)
        _attach_mnx_chebi_aliases(table, conn)

        # N — TCDB Pfam / GO bridges (existing loop, unchanged; it only fills `tcdb_bridges`)
        tcdb_bridges = _tcdb_bridges(cache_root, hierarchy, kept, prev_pruned, refetch_raw)

        # M inputs — KEGG compound BRITE (raw, else committed on outage)
        brite_raw, brite_committed = _load_compound_brite(cache_root, refetch_raw)

        old_resolved = _read_old_resolved(tables)

    inp = Step6Inputs(
        table=table,
        kegg_rest={k: kegg_data[k] for k in ("kos", "pathways", "subcategories", "categories", "reactions")},
        subtree_substrates=subtree_substrates, chemical_roles=chemical_roles, subtree_roles=subtree_roles,
        rekeys=rekeys, role_of_old=role_of_old, tables=tables, paper_overrides=paper_overrides,
        name_candidates=paper.name_candidates, lock=lock, lumping_pin=prev_pruned.get("lumping_pin"),
        identity_overrides=load_identity_overrides(IDENTITY_OVERRIDES_FILE), chebi=chebi,
        chebi_from_raw=chebi_from_raw, rxn_genes=catalysis_sets["rxn_genes"],
        compound_to_pathways=cpd_to_pw, base_pathways=catalysis_pws | transport_pws,
        brite_raw=brite_raw, brite_committed=brite_committed, tcdb_kept=kept,
        tcdb_levels={t: int(hierarchy.get(t, {}).get("level", 0)) for t in kept},
        tcdb_static={"kept_tcdb_ids": sorted(kept), "pfam_bridge": tcdb_bridges["pfam"],
                     "go_bridge": tcdb_bridges["go"], "seed_aliases": dict(sorted(upgraded.items()))},
        previous_ids=previous_ids, old_resolved=old_resolved, today=today)
    out = run_identity_stages(inp)               # raises Step6ValidationError → nothing written
    files = render_outputs(out, Step6Paths(
        kegg_data=KEGG_DATA_FILE, tcdb_pruned=prev_pruned_path, brite_dir=KEGG_CACHE_DIR,
        chebi_cache=CHEBI_CACHE_FILE, mapping=MAPPING_FILE, lock=LOCK_FILE, report=IDENTITY_REPORT_FILE))
    write_all_atomic(files)
    c = out.report["counts"]
    log.info(f"  step 6 wrote {len(files)} files: {c['metabolites']} metabolites "
             f"({c['merged_away']} merged away, {c['rekeyed']} re-keyed, {c['dropped']} dropped), "
             f"{c['roles']} ChemicalRole, {c['class_nodes']} MetaboliteClass")
    log.info("Done.")
```

  Supporting helpers to add in the same module:
  - `_tcdb_bridges(cache_root, hierarchy, kept, previous, refetch_raw) -> dict[str, dict]`: the existing bridge loop (current lines 1474–1543) moved verbatim into a function that **returns** `tcdb_bridges` instead of writing `tcdb_pruned.json`. Its guards keep raising `RuntimeError`, which happens before any stage-P write.
  - `_load_compound_brite(cache_root, refetch_raw) -> tuple[dict, dict]`: `try: raw = download_compound_brite(cache_root / "kegg" / "raw", force=refetch_raw)` (Task 8); on exception log a warning and use `{}`; committed = `{code: json.loads(p.read_text()) for code in COMPOUND_BRITE_TREES if (p := KEGG_CACHE_DIR / f"brite_{code}.json").exists()}`.
  - `_read_old_resolved(tables) -> dict[tuple[str, str, int], str]`: for each table read the existing `<stem>_resolved.csv` (if present) with the same `read_csv` arguments and map `(paper, table, row_index)` → its `metabolite_id` (empty → skip).
  - Imports: `from datetime import date`; the `metabolite_identity` names used above; `write_all_atomic` from `outputs`.
  - Delete the now-dead tail of the old `main()` (the old `_fold_paper…` call, the mapping write, the `KEGG_DATA_FILE.write_text`, the old `tcdb_pruned.json` write).
  - The `if __name__ == "__main__"` block is unchanged (`today` is not a CLI flag).

- [ ] **Step 4: Retire step 7.** Replace the body of `multiomics_kg/download/resolve_paper_metabolites.py` with:

```python
"""Retired (2026-10): step 7 is folded into step 6.

Step 6 (build_kegg_metabolism_xrefs) decides every paper-metabolite resolution and writes
each <stem>_resolved.csv from the same in-memory decisions, so the graph's metabolite set
and the paper tables can't drift (spec 2026-10-13 §5.3). Re-resolve a paper with
`bash scripts/prepare_data.sh --steps 6 --force` (about a minute, no network).
"""
import sys


def main(argv=None) -> int:
    print("step 7 is folded into step 6 (build_kegg_metabolism_xrefs); nothing to do")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 5: Update `scripts/prepare_data.sh`.**
  - Default `STEPS="0 9 1 2 3 4 5 6 8 10"` (line ~114) and the `--rebuild` list `STEPS="9 1 2 3 4 5 6 8 10"` (lines ~144–151).
  - The `7)` case: `echo "[step 7] folded into step 6 — nothing to do (re-run --steps 6 --force to re-resolve papers)"`; no `run_step`.
  - Header comments: step 6 now "writes kegg_data.json, tcdb_pruned.json (incl. chemical_roles / lumping_pin), the 8 KEGG compound BRITE files, the ChEBI cache, metabolite_id_mapping.json, the dominant-stereo-form lock, the identity report, and every paper `<stem>_resolved.csv`"; step 7 "retired, folded into step 6 (accepted as a no-op)"; fix "3 TCDB reference TSVs" → 4; usage line `--steps 6 7 --force` → `--steps 6 --force`; the L199 echo line.

- [ ] **Step 6: Run** `uv run pytest -m "not slow and not kg" -q` → PASS (any failure in files this plan hasn't reached yet, e.g. vocabulary tests, is fixed in its own task; note it and continue only if it's in a later task's scope).

- [ ] **Step 7: Commit.**

```bash
git add multiomics_kg/download/build_kegg_metabolism_xrefs.py multiomics_kg/download/resolve_paper_metabolites.py scripts/prepare_data.sh tests/test_mi_outputs.py tests/test_prepare_data_step7_alias.py
git rm tests/test_resolve_paper_metabolites.py
git commit -m "feat(step6): stages A-P end to end, single atomic write; retire step 7"
```

### Task 14: Curate overrides, run step 6 for real, re-measure, commit caches

**Files:**
- Modify: `data/Prochlorococcus/papers_and_supp/{biller 2022,Capovilla 2023,Kujawinski 2023}/metabolite_aliases.yaml`, `config/metabolite_identity_overrides.yaml` (only if the review finds a same-entity drop)
- Generated + committed: `cache/data/kegg/kegg_data.json`, `cache/data/kegg/brite_br08{001,002,003,005,006,007,009,021}.json`, `cache/data/tcdb/tcdb_pruned.json`, `cache/data/chebi/chebi_metabolite_xrefs.json`, `cache/data/metabolomics/{metabolite_id_mapping,dominant_stereo_forms,metabolite_identity_report}.json`, the 7 paper `*_resolved.csv` + `*_resolution_report.json`

- [ ] **Step 1: Convert every existing alias entry to the explicit form.** For each of the 26 entries in the three files, write `{id, stereo_assignment, reason}`, keeping the existing comment text as `reason`. `stereo_assignment`: `as_reported` when the name is stereo-qualified (`is_stereo_qualified`) or the target has no stereo variant; otherwise `dominant_form_assumed`. Targets known to have variants: `C00019` (SAM; "Adenosyl Methionine" is unqualified → `dominant_form_assumed`), `C00041` ("alanine (isom. sarcosine)" → `dominant_form_assumed`), `C00188` threonine (`dominant_form_assumed`), `C02713` MurNAc (anomer-only → `as_reported`). Step 6 enforces the rest: an unconverted legacy entry on a stereo target fails with an explicit message.

- [ ] **Step 2: Add the spec §3.6 curated overrides.**
  - Biller 2022 and Kujawinski 2023: `Glutathione` / `glutathione` → `{id: "kegg.compound:C00051", stereo_assignment: as_reported, reason: "GSH, reduced"}`; `Glutathione Disulfide` / `glutathione oxidized` → `{id: "kegg.compound:C00127", stereo_assignment: as_reported, reason: "GSSG"}`.
  - The paper whose table carries `Cystine` and `indole 3-acetic acid` (both landed on `C00954`): `Cystine` → `kegg.compound:C00491` (L-cystine; `dominant_form_assumed`), `indole 3-acetic acid` → `kegg.compound:C00954` (`as_reported`).
  - Placeholders: `Arginine` → `C00062` (L-Arg), `cysteine` → `C00097` (L-Cys), `Serine` → `C00065` (L-Ser), each `dominant_form_assumed`, reason "paper gives no stereo; was an mnx: placeholder"; the compartment-worded names (`Arginine hydrochloride`, `Chlorplastic Serine`, …: read the actual names from `metabolite_id_mapping.json` v1 `compounds["mnx:…"]["aliases"]` before the rerun) likewise; `Cytosolic Ribose 5 Phosphate` → `kegg.compound:C00117` (D-ribose 5-phosphate, `as_reported`).
  - Confirm (spec §3.6 "placeholder collision check") that the assay node's own `compartment` field carries the compartment for the "Chlorplastic" / "Cytosolic" rows: `grep -n compartment "data/Prochlorococcus/papers_and_supp/<paper>/paperconfig.yaml"` on the tables those rows sit in.

- [ ] **Step 3: Run step 6 and iterate until it writes.**

```bash
bash scripts/prepare_data.sh --steps 6 --force; tail -40 logs/prepare_data_step6.log
```

  A `Step6ValidationError` lists each problem. Ambiguous names: add an override that names "the form the paper names" (spec §3.6). The 36 names measured on 2026-10-14 and their MNX candidates are below; after canonicalisation many collapse to one node, and only those still listed in the error need an entry:

| Name(s) | Candidates (2026-10-14) | Expected outcome |
|---|---|---|
| S-adenosyl-L-methionine | `chebi:142094`, `C00019` | merge (same entity) or override → `C00019` |
| chitotriose | `chebi:143141`, `chebi:84373`, 2 `mnx:` | override → the chitotriose term the paper means |
| Tyrosine / tyrosine | `chebi:173245`, `chebi:32783`, `C00082` | rule → L-Tyr `C00082` |
| Vitamin B2 | `chebi:176838`, `C00255` | override → riboflavin `C00255` |
| Vitamin B5 | `chebi:176840`, `C00864`, `chebi:14739` | override → pantothenate `C00864` |
| Cysteic Acid | `chebi:21260`, `C00506` | override → L-cysteate `C00506` |
| Glucosylglycerol | `chebi:91033`, `chebi:24287` | override → the paper's form |
| acetyl coenzyme A / UDP-glucose / Cytidine / UDP-glucosamine | KEGG id + 1 `mnx:` | merge or override → the KEGG id |
| Lysine | `C00047`, `mnx:MNXM1565193` | → `C00047`, `dominant_form_assumed` |
| D-Ribose 5-phosphate | `C00117`, `chebi:13013`, `chebi:78679` | → `C00117` |
| Glutathione (+ variants) | `C00127`, `C00051`, 2 `chebi:` | overrides above |
| Alanine / Glutamic acid / Tryptophan / Glutamine / Methionine / Phenylalanine / Asparagine | family members | rule → L-form (7 moves) |
| Proline / Adenosine / Homoserine / D-glucosamine / D-glucosamine 6-phosphate | KEGG id + 1 alternative | rule or override |
| Aconitic Acid | `C00417`, `C02341`, `chebi:22210` | override → cis-aconitate `C00417` (geometric isomers are separate ids) |
| Chitobiose (+ pos) | `C01674` + 8 | existing override → `C01674` |
| Methionine Sulfoxide | `C02989` + 3 `chebi:` | override → `C02989` |
| Vitamin B3 | `C03150`, `chebi:176839`, `C00153`, `C00253` | override → the paper's form (nicotinate `C00253` vs nicotinamide `C00153`); check the source table's context |
| Vitamin K1 | `C02059`, `C03313` | override → phylloquinone `C02059` |
| Aspartic acid / aspartic acid | `C16433`, `C00402`, `C00049` | rule → L-Asp `C00049` |
| Arginine / cysteine / Serine | `mnx:` placeholder + family | overrides above |

  For each override decision the reason must say why ("paper's Table S7 lists it beside nicotinamide riboside, so nicotinamide"). Re-run until step 6 writes.

- [ ] **Step 4: Review the identity report.**

```bash
uv run python - <<'EOF'
import json
r = json.load(open("cache/data/metabolomics/metabolite_identity_report.json"))
print(json.dumps(r["counts"], indent=1))
for k in ("rejected_links", "refused_assay_collisions", "stale_alias_override", "stale_identity_override",
          "stale_lock_entry", "dominant_form_drift", "unaccounted_old_ids", "dropped_same_entity_review",
          "heuristic_name_candidates", "unresolved_paper_names"):
    print(k, len(r[k])); print("  ", r[k][:12])
print("lumping", r["lumping"]); print("structure", r["structure"]["counts"]); print("stereo", r["stereo_counts"])
EOF
```

  Required before committing:
  - `unaccounted_old_ids` is **empty** (every pre-change id is in the remap with a basis).
  - **Same-entity review of every drop** (spec §3.6, Review Focus 5): for each `dropped_same_entity_review` row, decide whether `dropped` and `target` are one entity (e.g. a ChEBI node naming the same molecule as the KEGG target). If yes, add a `merge:` entry to `config/metabolite_identity_overrides.yaml` (reason required) and re-run: the id moves from `recurated_from_ids` into `all_identifiers`. D-Trp `C00525` → L-Trp `C00078` is a stereoisomer, not the same entity: it stays a drop.
  - `heuristic_name_candidates`: read each pair; record a human decision only where it's clearly the same entity (`merge:`) or clearly not (`keep_separate:`); leave the rest for the report.
  - `stale_*` lists: expected empty except legitimately stale entries; fix typos in alias files.
  - `refused_assay_collisions`: 0 expected (spec §3.1).
  - `rejected_links` includes water / hydroxide (`C00001` / `C01328`) and cyanide (`C00177` / `C01326`).

- [ ] **Step 5: Re-measure against the plan's baseline** and record the numbers in the plan's measured-baseline table under a new column "after step 6" (edit this file):

```bash
uv run python - <<'EOF'
import json, collections
kd = json.load(open("cache/data/kegg/kegg_data.json"))
allc = {**{f"kegg.compound:{k}": v for k, v in kd["compounds"].items()}, **kd["additional_compounds"]}
print("metabolites", len(allc), collections.Counter(k.split(":")[0] for k in allc))
print("structure_kind", collections.Counter(v["structure_kind"] for v in allc.values()))
print("with elements", sum(1 for v in allc.values() if v.get("elements")))
print("stereo_form_of", sum(len(v.get("stereo_form_of", [])) for v in allc.values()),
      "enantiomer_of", sum(len(v.get("enantiomer_of", [])) for v in allc.values()))
print("metabolomics-tagged", sum("metabolomics" in v["evidence_sources"] for v in allc.values()))
tp = json.load(open("cache/data/tcdb/tcdb_pruned.json"))
print("roles", len(tp["chemical_roles"]), "role edges", sum(map(len, tp["subtree_roles"].values())),
      "substrate edges", sum(map(len, tp["subtree_substrates"].values())))
print("lumping", tp["lumping_pin"]["lumping_tcdb_ids"])
for c in ("br08001", "br08002", "br08003", "br08005", "br08006", "br08007", "br08009", "br08021"):
    b = json.load(open(f"cache/data/kegg/brite_{c}.json")); print(c, len(b["nodes"]), len(b["members"]))
EOF
```

  Expectations (stop and investigate, don't commit, if far off): metabolites ≈ 3,214 (−57 roles, −≈ 74 merged, −≈ 11 dropped from 3,356); 0 roles in `additional_compounds`; role edges = 392; substrate edges ≈ 11,263 − 392 − merge dedup; lumping list = the 13 seeds; BRITE ≈ 431 nodes / ≈ 1,130 members (fewer only by merges); `metabolomics`-tagged = number of distinct final targets.

- [ ] **Step 6: Determinism check.** `cp -r cache/data /tmp/claude-1001/-home-osnat-github-second-multiomics/ee41053f-9e3e-4291-be5c-1e845f384bc4/scratchpad/step6_run1 && bash scripts/prepare_data.sh --steps 6 --force && diff -r cache/data/kegg /…/step6_run1/kegg && diff -r cache/data/tcdb /…/step6_run1/tcdb && diff -r cache/data/metabolomics /…/step6_run1/metabolomics` → no output (Review Focus 4). Also `git diff --stat -- 'data/**/*_resolved.csv'` shows only the expected tables.

- [ ] **Step 7: Commit the curated files and the generated caches** (one commit, data only):

```bash
git add "data/Prochlorococcus/papers_and_supp/biller 2022/metabolite_aliases.yaml" \
        "data/Prochlorococcus/papers_and_supp/Capovilla 2023/metabolite_aliases.yaml" \
        "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/metabolite_aliases.yaml" \
        config/metabolite_identity_overrides.yaml cache/data/kegg/kegg_data.json cache/data/kegg/brite_br08*.json \
        cache/data/tcdb/tcdb_pruned.json cache/data/chebi/chebi_metabolite_xrefs.json \
        cache/data/metabolomics/*.json
git add -u data/Prochlorococcus/papers_and_supp      # the regenerated *_resolved.csv / *_resolution_report.json
git status --short | grep -v '^??'                     # nothing unexpected staged
git commit -m "chore(cache): step 6 rerun — canonical metabolites, ChEBI cache, compound BRITE, lock, curated aliases"
```

---

## Phase 3 — Schema, vocabularies, adapters

### Task 15: Schema, controlled vocabularies, `kegg` / `chebi` DataSources

**Files:**
- Modify: `config/schema_config.yaml`, `config/controlled_vocabularies.yaml`, `config/gene_annotations_config.yaml`, `multiomics_kg/adapters/data_source_adapter.py`
- Test: `tests/test_data_source_adapter.py`, `tests/kg_validity/test_data_source.py`, `tests/test_controlled_vocab.py`, `tests/test_annotation_provenance.py`

**Interfaces:**
- Produces: `data_source_adapter.data_source_ids(config: dict) -> set[str]` (union of `sources.*.logical_sources[].id` and `reference_sources[].id`), used by every test that pins `sources` values; DataSource nodes `data_source:kegg`, `data_source:chebi`.

- [ ] **Step 1: Update the tests first.**
  - `tests/test_data_source_adapter.py::test_emits_nine_nodes` → rename `test_emits_eleven_nodes`; add `"kegg", "chebi"` to the expected set; new test `test_reference_sources_have_info_types_and_chebi_version` asserting both nodes have non-empty `info_types`, `scope == "global"`, `provenance == "download"`, and `data_source:chebi` `version == "256"` when `cache/data/chebi/chebi_metabolite_xrefs.json` exists (skip otherwise), and that the chebi `description` contains `"CC BY 4.0"`.
  - `tests/kg_validity/test_data_source.py::test_nine_data_source_nodes` → `test_eleven_data_source_nodes`, expected sorted list plus `chebi`, `kegg`.
  - `tests/test_controlled_vocab.py::test_sources_descriptions_do_not_drift_from_gene_annotations_config`: replace the inline `logical = {…}` with `logical = data_source_ids(cfg)`. Add `"identity_basis", "structure_kind", "substrate_breadth", "stereo_assignment"` to the `trust_props` set two tests above it.
  - `tests/test_annotation_provenance.py::test_source_values_are_all_data_source_ids`: use `data_source_ids(cfg)`.

- [ ] **Step 2: Run** `uv run pytest tests/test_data_source_adapter.py tests/test_controlled_vocab.py tests/test_annotation_provenance.py -q` → FAIL.

- [ ] **Step 3: Registry + adapter.** Append to `config/gene_annotations_config.yaml` after the `sources:` block (before `fields:`):

```yaml
# Non-gene reference sources: DataSource nodes for provenance joins (R2) on edges that
# never touch a gene (metabolite identity, stereo relations, compound classes). Not read
# by the gene-annotation merge (build_gene_annotations reads only `sources` / `fields`).
reference_sources:
  - id: kegg
    scope: global
    provenance: download
    info_types_extra: [compound_classes]
  - id: chebi
    scope: global
    provenance: download
    info_types_extra: [metabolite_identity, stereo_relations, chemical_roles]
```

  In `data_source_adapter.py`: add `data_source_ids(config)`; make `_all_logical_sources` also yield `{"_yaml_source_name": ls["id"], **ls}` for each `reference_sources` entry; add `"kegg": "KEGG"`, `"chebi": "ChEBI"` to `_name_for`; add descriptions:
  - `kegg`: `"KEGG COMPOUND and its compound BRITE classification trees (br08001-br08021), via the KEGG REST API; source of MetaboliteClass memberships."`
  - `chebi`: `"ChEBI (Chemical Entities of Biological Interest, EMBL-EBI): metabolite identity, stereo relations and the role ontology behind ChemicalRole. Licensed CC BY 4.0; attribution: ChEBI, EMBL-EBI."`
  - `version`: `""` except `chebi`, which reads `release` from `cache/data/chebi/chebi_metabolite_xrefs.json` when present (`_version_for(src_id)`).

- [ ] **Step 4: Schema** (`config/schema_config.yaml`). Under `metabolite:` properties add, keeping existing ones:

```yaml
    kegg_compound_ids: str[]     # bare, scalar value first; present iff kegg_compound_id is
    chebi_ids: str[]             # bare primary ChEBI numbers, scalar value first; present iff chebi_id is
    mnxm_ids: str[]              # sparse: only when the node spans > 1 MNX entity
    all_identifiers: str[]       # dense, CURIE-only (kegg.compound: / chebi: / mnx: / hmdb: / inchikey:), unlike Gene.all_identifiers (bare)
    identifier_tokens: str[]     # post-import, search-only: bare local parts of all_identifiers (metaboliteFullText)
    recurated_from_ids: str[]    # sparse, NO identity meaning: dropped nodes whose paper rows now point here
    identity_basis: str[]        # sparse: merged or re-keyed nodes only
    structure_kind: str          # dense: defined | generic | polymer
    kegg_compound_categories: str[]   # post-import, sparse: br08001 level-0 names
    measured_compartments: str[] # post-import (already set; declared here for completeness)
```

  Under `tcdb family:` add `substrate_breadth: str  # dense: lumping | specific | not_assessed (pinned set, spec 2026-10-13 §3.4a; an R3 exception)` and replace the comment that tells consumers to query `level >= 2 AND metabolite_count >= 50` with "filter on substrate_breadth; the threshold query no longer reproduces the pinned set". Add `direct_gene_count: int` (already set by post-import, undeclared).

  New entries (append near `brite category`):

```yaml
metabolite class:
  is_a: named thing
  represented_as: node
  preferred_id: kegg.brite
  label_in_input: metabolite class
  properties:
    name: str
    tree: str                    # snake_case KEGG compound tree name
    tree_code: str               # br08001 … br08021
    level: int                   # 0 = broadest
    level_kind: str              # kegg_compound_class | _subclass | _group | _subgroup
    metabolite_count: int        # post-import: subtree
    direct_metabolite_count: int # post-import
    organism_count: int          # post-import

metabolite class is a metabolite class:
  is_a: association
  represented_as: edge
  label_as_edge: metabolite_class_is_a_metabolite_class
  source: metabolite class
  target: metabolite class
  label_in_input: metabolite_class_is_a_metabolite_class

metabolite in metabolite class:
  is_a: association
  represented_as: edge
  label_as_edge: metabolite_in_metabolite_class
  source: metabolite
  target: metabolite class
  label_in_input: metabolite_in_metabolite_class
  properties:
    sources: str[]
    evidence: str

chemical role:
  is_a: named thing
  represented_as: node
  preferred_id: chebi.role       # house-minted prefix; never routed as a chebi: metabolite
  label_in_input: chemical role
  properties:
    name: str
    level: int                   # shortest is_a distance to a ChEBI role branch (DAG: approximate)
    level_kind: str              # chebi_role_branch | chebi_role
    all_identifiers: str[]       # chebi.role:<n> + chebi:<primary and secondary ids>
    transporter_count: int       # post-import, subtree, most_specific edges, DISTINCT families
    direct_transporter_count: int
    tcdb_family_count: int       # post-import, subtree, any depth, DISTINCT families
    direct_tcdb_family_count: int

chemical role is a chemical role:
  is_a: association
  represented_as: edge
  label_as_edge: chemical_role_is_a_chemical_role
  source: chemical role
  target: chemical role
  label_in_input: chemical_role_is_a_chemical_role

tcdb family transports chemical role:
  is_a: association
  represented_as: edge
  label_as_edge: tcdb_family_transports_chemical_role
  source: tcdb family
  target: chemical role
  label_in_input: tcdb_family_transports_chemical_role
  properties:
    # "This family transports something TCDB calls <role>" — NOT a role of the family or
    # of its metabolites. Same rollup + depth rule as tcdb_family_transports_metabolite.
    substrate_depth: str

metabolite is stereo form of metabolite:
  is_a: association
  represented_as: edge
  label_as_edge: metabolite_is_stereo_form_of
  source: metabolite
  target: metabolite
  label_in_input: metabolite_is_stereo_form_of
  properties:
    sources: str[]
    evidence: str

metabolite is enantiomer of metabolite:
  is_a: association
  represented_as: edge
  label_as_edge: metabolite_is_enantiomer_of     # symmetric: stored once, lower → higher id
  source: metabolite
  target: metabolite
  label_in_input: metabolite_is_enantiomer_of
  properties:
    sources: str[]
    evidence: str
```

  Under `assay quantifies metabolite:` and `assay flags metabolite:` properties add `stereo_assignment: str   # as_reported | dominant_form_assumed (dense on every edge)`.

- [ ] **Step 5: Vocabularies** (`config/controlled_vocabularies.yaml`). Add these entries (each `applies_to` / `applies_to_kind` / `property` / `value_type` / `closed` / `values` / `description`, plus `sparse: true` where marked and `value_descriptions` where listed; descriptions are researcher-facing, provenance goes in a `# provenance:` comment, per DOC-003):

| Key | type | closed | sparse | values | value_descriptions |
|---|---|---|---|---|---|
| `Metabolite.identity_basis` | string_array | yes | yes | `mnx_same_entity, chebi_kegg_xref, chebi_conjugate, chebi_tautomer, manual` | `mnx_same_entity`: MetaNetX reconciles the merged ids as one entity · `chebi_kegg_xref`: a ChEBI term and the KEGG compound it cross-references · `chebi_conjugate`: conjugate acid / base forms of one physiological species · `chebi_tautomer`: tautomers of one compound · `manual`: a curated merge with a recorded reason |
| `Metabolite.structure_kind` | string | yes | no | `defined, generic, polymer` | `defined`: one concrete structure · `generic`: a class or a structure with R groups; no single formula mass · `polymer`: a repeat-unit polymer; formula is the repeat unit |
| `Metabolite.kegg_compound_categories` | string_array | yes | yes | the 9 `BR08001_LEVEL0` names | — |
| `MetaboliteClass.level_kind` | string | yes | no | `kegg_compound_class, kegg_compound_subclass, kegg_compound_group, kegg_compound_subgroup` | — |
| `MetaboliteClass.tree` | string | yes | no | the 8 `COMPOUND_BRITE_TREES` values | — |
| `ChemicalRole.level_kind` | string | yes | no | `chebi_role_branch, chebi_role` | — |
| `Tcdb_family_transports_chemical_role.substrate_depth` | string | yes | no | `most_specific, inherited` | copy the two descriptions from `Tcdb_family_transports_metabolite.substrate_depth`, replacing "substrate" with "role-typed substrate" |
| `TcdbFamily.substrate_breadth` | string | yes | no | `lumping, specific, not_assessed` | `lumping`: in the pinned substrate-lumping set — a substrate count here is reachability, not capability · `specific`: level ≥ 2 and not lumping · `not_assessed`: tc_class / tc_subclass, where substrate counts scale with level by construction |
| `Metabolite_in_metabolite_class.sources` / `.evidence` | string_array / string | yes | no | `[kegg]` / `[curated]` | `kegg`: KEGG compound BRITE trees · `curated`: assigned by a human curator or a curated reference database |
| `Metabolite_is_stereo_form_of.sources` / `.evidence` | string_array / string | yes | no | `[chebi]` / `[curated]` | `chebi`: ChEBI is_a between the two terms, with matching InChIKey connectivity · `curated`: as above |
| `Metabolite_is_enantiomer_of.sources` / `.evidence` | string_array / string | yes | no | `[chebi]` / `[curated]` | `chebi`: ChEBI is_enantiomer_of · `curated`: as above |
| `Assay_quantifies_metabolite.stereo_assignment`, `Assay_flags_metabolite.stereo_assignment` | string | yes | no | `as_reported, dominant_form_assumed` | `as_reported`: the paper named the form, or the metabolite has no stereo variant other than anomers · `dominant_form_assumed`: the paper gave no stereo; the row was placed on one form of a stereo family (usually the one with the most catalysing genes) — treat the form as an assumption |

  The `sources` descriptions use the exact text keyed by the DataSource ids, so `test_sources_descriptions_do_not_drift_from_gene_annotations_config` passes.

- [ ] **Step 6: Run** `uv run pytest tests/test_data_source_adapter.py tests/test_controlled_vocab.py tests/test_controlled_vocabulary_adapter.py tests/test_annotation_provenance.py -q` → PASS.

- [ ] **Step 7: Commit.**

```bash
git add config/schema_config.yaml config/controlled_vocabularies.yaml config/gene_annotations_config.yaml multiomics_kg/adapters/data_source_adapter.py tests/test_data_source_adapter.py tests/kg_validity/test_data_source.py tests/test_controlled_vocab.py tests/test_annotation_provenance.py
git commit -m "feat(schema): MetaboliteClass, ChemicalRole, stereo edges, identity props; kegg/chebi DataSources"
```

### Task 16: `metabolism_adapter.py` — identity fields, elements allow-list, stereo + enantiomer edges

**Files:**
- Modify: `multiomics_kg/adapters/metabolism_adapter.py`
- Test: `tests/test_metabolism_adapter.py`, `tests/test_metabolism_adapter_elements.py`

**Interfaces:**
- Consumes: `metabolite_identity.names.parse_elements` (Task 2); `kegg_data.json` shape from Task 13.
- Produces: `MetabolismAdapter.all_metabolite_ids() -> set[str]` (ids actually emitted by `get_nodes`, or every compound id if `get_nodes` hasn't run); edges `metabolite_is_stereo_form_of` (id `m2s:{src}:{tgt}`) and `metabolite_is_enantiomer_of` (id `m2e:{lower}:{higher}`), both `{sources: ['chebi'], evidence: 'curated'}`.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_metabolism_adapter.py`; build the kegg_data dict inline and write it to `tmp_path`):

```python
def _kd_identity(tmp_path):
    import json
    kd = {"reactions": {}, "kos": {}, "pathways": {}, "subcategories": {}, "categories": {},
          "compounds": {
              "C00041": {"name": "L-Alanine", "formula": "C3H7NO2", "elements": ["C", "H", "N", "O"],
                         "kegg_compound_ids": ["C00041"], "chebi_id": "16977", "chebi_ids": ["16977"],
                         "all_identifiers": ["chebi:16977", "kegg.compound:C00041"], "structure_kind": "defined",
                         "evidence_sources": ["metabolism", "metabolomics"],
                         "stereo_form_of": ["kegg.compound:C01401"], "enantiomer_of": ["kegg.compound:C00133"]},
              "C00133": {"name": "D-Alanine", "structure_kind": "defined", "evidence_sources": ["metabolism"],
                         "all_identifiers": ["kegg.compound:C00133"], "stereo_form_of": ["kegg.compound:C01401"]},
              "C01401": {"name": "Alanine", "structure_kind": "defined", "evidence_sources": ["metabolism"],
                         "all_identifiers": ["kegg.compound:C01401"]},
              "C00078": {"name": "L-Tryptophan", "structure_kind": "defined", "evidence_sources": ["metabolomics"],
                         "all_identifiers": ["kegg.compound:C00078"], "recurated_from_ids": ["kegg.compound:C00525"]}},
          "additional_compounds": {
              "chebi:17522": {"name": "alditol", "formula": "C2H6O2(CH2O)n", "elements": ["C", "H", "O"],
                              "structure_kind": "polymer", "chebi_id": "17522", "chebi_ids": ["17522"],
                              "all_identifiers": ["chebi:17522"], "evidence_sources": ["transport"]}}}
    p = tmp_path / "kegg_data.json"
    p.write_text(json.dumps(kd))
    return p


def test_identity_fields_copied(tmp_path):
    from multiomics_kg.adapters.metabolism_adapter import MetabolismAdapter
    a = MetabolismAdapter(kegg_data_path=_kd_identity(tmp_path))
    nodes = {n[0]: n[2] for n in a.get_nodes() if n[1] == "metabolite"}
    ala = nodes["kegg.compound:C00041"]
    assert ala["all_identifiers"] == ["chebi:16977", "kegg.compound:C00041"]
    assert ala["structure_kind"] == "defined" and ala["chebi_ids"] == ["16977"]
    assert nodes["kegg.compound:C00078"]["recurated_from_ids"] == ["kegg.compound:C00525"]
    assert nodes["chebi:17522"]["elements"] == ["C", "H", "O"] and "mass" not in nodes["chebi:17522"]
    assert "stereo_form_of" not in ala and "enantiomer_of" not in ala      # edges, not properties


def test_stereo_and_enantiomer_edges(tmp_path):
    from multiomics_kg.adapters.metabolism_adapter import MetabolismAdapter
    a = MetabolismAdapter(kegg_data_path=_kd_identity(tmp_path))
    list(a.get_nodes())
    edges = {e[0]: e for e in a.get_edges()}
    s = edges["m2s:kegg.compound:C00041:kegg.compound:C01401"]
    assert s[3] == "metabolite_is_stereo_form_of" and s[4] == {"sources": ["chebi"], "evidence": "curated"}
    e = edges["m2e:kegg.compound:C00041:kegg.compound:C00133"]
    assert e[3] == "metabolite_is_enantiomer_of" and e[1] < e[2]
    assert not any(k.startswith("m2e:kegg.compound:C00133") for k in edges)


def test_stereo_edge_guard_skips_unemitted_targets(tmp_path):
    from multiomics_kg.adapters.metabolism_adapter import MetabolismAdapter
    a = MetabolismAdapter(kegg_data_path=_kd_identity(tmp_path))
    list(a.get_nodes())
    a._emitted_metabolite_ids.discard("kegg.compound:C01401")
    assert not any(e[2] == "kegg.compound:C01401" for e in a.get_edges() if e[3] == "metabolite_is_stereo_form_of")
```

  And in `tests/test_metabolism_adapter_elements.py` add `("Rb", ["Rb"])`, `("C5H8O2Rn", ["C", "H", "O"])`, `("C3H3O6R", ["C", "H", "O"])` to its parametrized cases (or add a new parametrized test calling `_parse_elements`).

- [ ] **Step 2: Run** `uv run pytest tests/test_metabolism_adapter.py tests/test_metabolism_adapter_elements.py -q` → FAIL.

- [ ] **Step 3: Implement.**
  - `_parse_elements(formula)` → `return parse_elements(formula)` (import from `multiomics_kg.utils.metabolite_identity.names`), keeping the function name so the elements tests keep importing it.
  - A single `_metabolite_props(cpd: dict, kegg_id: str | None) -> dict` replacing the two duplicated property blocks:

```python
_MET_STR = ("name", "formula", "inchikey", "smiles", "mnxm_id", "chebi_id", "hmdb_id", "structure_kind")
_MET_LIST = ("kegg_compound_ids", "chebi_ids", "mnxm_ids", "all_identifiers", "recurated_from_ids",
             "identity_basis", "evidence_sources")


def _metabolite_props(cpd: dict, kegg_id: str | None) -> dict:
    props = {"kegg_compound_id": kegg_id, "mass": cpd.get("mass"),
             "elements": cpd.get("elements") or _parse_elements(cpd.get("formula"))}
    props.update({k: _clean_str(cpd.get(k)) for k in _MET_STR})
    props.update({k: [_clean_str(v) for v in (cpd.get(k) or [])] for k in _MET_LIST})
    return _drop_nulls(props)
```

  - `get_nodes`: use `_metabolite_props(cpd, cpd_id)` for `compounds` and `_metabolite_props(cpd, None)` for `additional_compounds`; record every yielded metabolite id in `self._emitted_metabolite_ids: set[str]` (initialised in `__init__`). The test-mode cap stays.
  - `all_metabolite_ids()`: `return set(self._emitted_metabolite_ids) or {f"kegg.compound:{k}" for k in data["compounds"]} | set(data.get("additional_compounds", {}))`.
  - Edge emitters:

```python
    def _metabolite_relation_edges(self):
        data = self._load()
        known = self.all_metabolite_ids()
        entries = [(f"kegg.compound:{k}", v) for k, v in data["compounds"].items()]
        entries += list(data.get("additional_compounds", {}).items())
        skipped = 0
        for src, cpd in entries:
            if src not in known:
                continue
            for tgt in cpd.get("stereo_form_of") or []:
                if tgt not in known:
                    skipped += 1
                    continue
                yield (f"m2s:{src}:{tgt}", src, tgt, "metabolite_is_stereo_form_of",
                       {"sources": ["chebi"], "evidence": "curated"})
            for tgt in cpd.get("enantiomer_of") or []:
                if tgt not in known:
                    skipped += 1
                    continue
                assert src < tgt, f"enantiomer pair stored high→low: {src} {tgt}"
                yield (f"m2e:{src}:{tgt}", src, tgt, "metabolite_is_enantiomer_of",
                       {"sources": ["chebi"], "evidence": "curated"})
        if skipped:
            logger.info(f"metabolism: skipped {skipped} stereo/enantiomer edges to unemitted metabolites")
```

  and `yield from self._metabolite_relation_edges()` at the end of `MetabolismAdapter.get_edges` (the Multi subclass already delegates to it). Use the module's existing logger name.

- [ ] **Step 4: Run** `uv run pytest tests/test_metabolism_adapter.py tests/test_metabolism_adapter_elements.py tests/test_metabolism_smoke.py -q` → PASS.

- [ ] **Step 5: Commit.**

```bash
git add multiomics_kg/adapters/metabolism_adapter.py tests/test_metabolism_adapter.py tests/test_metabolism_adapter_elements.py
git commit -m "feat(metabolism-adapter): identity fields, element allow-list, stereo + enantiomer edges"
```

### Task 17: `tcdb_adapter.py` — `ChemicalRole`, role edges, depth refactor, `substrate_breadth`

**Files:**
- Modify: `multiomics_kg/adapters/tcdb_adapter.py`
- Test: `tests/test_tcdb_adapter.py`

**Interfaces:**
- Consumes: `tcdb_pruned.json` blocks `chemical_roles`, `subtree_roles`, `lumping_pin` (Task 13); `lumping.substrate_breadth` (Task 9).
- Produces: `MultiTcdbAnnotationAdapter._compute_substrate_depth(self, subtree_map: dict[str, list[str]]) -> dict[str, set[str]]`; `chemical role` nodes; `chemical_role_is_a_chemical_role` (id `{child}-isa-{parent}`) and `tcdb_family_transports_chemical_role` (id `t2r:{tc_id}:{role_id}`, `{substrate_depth}`) edges; `TcdbFamily.substrate_breadth`.

- [ ] **Step 1: Write the failing tests.** The file's `cache_root` fixture (`tests/test_tcdb_adapter.py:29`) writes a `tcdb_pruned.json` with kept ids `1 › 1.A › 1.A.1 › 1.A.1.5 › 1.A.1.5.2` and **no `lumping_pin`**. Because the adapter will now refuse a pruned file without a pin, first add these keys to that fixture's `tcdb_pruned.json` dict (every existing test keeps passing):

```python
        "chemical_roles": {
            "chebi.role:24432": {"name": "biological role", "parents": [], "level": 0,
                                 "level_kind": "chebi_role_branch",
                                 "all_identifiers": ["chebi.role:24432", "chebi:24432"]},
            "chebi.role:26672": {"name": "siderophore", "parents": ["chebi.role:24432"], "level": 1,
                                 "level_kind": "chebi_role",
                                 "all_identifiers": ["chebi.role:26672", "chebi:26672"]}},
        "subtree_roles": {"1.A.1.5": ["chebi.role:26672"], "1.A.1.5.2": ["chebi.role:26672"]},
        "lumping_pin": {"frozen_on": "2026-10-15",
                        "frozen_kept_ids": ["1", "1.A", "1.A.1", "1.A.1.5", "1.A.1.5.2"],
                        "lumping_tcdb_ids": ["1.A.1"]},
```

  Then add the tests (they use the file's `_make_orchestrator(cache_root, strain_dir)` helper, which returns a downloaded `MultiTcdbAnnotationAdapter`):

```python
def test_chemical_role_nodes_and_isa(cache_root, strain_dir):
    a = _make_orchestrator(cache_root, strain_dir)
    nodes = {n[0]: n for n in a.get_nodes()}
    r = nodes["chebi.role:26672"]
    assert r[1] == "chemical role" and r[2]["level_kind"] == "chebi_role" and r[2]["level"] == 1
    assert all(i.startswith(("chebi.role:", "chebi:")) for i in r[2]["all_identifiers"])
    assert ("chebi.role:26672-isa-chebi.role:24432", "chebi.role:26672", "chebi.role:24432",
            "chemical_role_is_a_chemical_role", {}) in list(a.get_edges())


def test_role_edges_share_the_substrate_depth_rule(cache_root, strain_dir):
    a = _make_orchestrator(cache_root, strain_dir)
    role = {e[1]: e[4]["substrate_depth"] for e in a.get_edges()
            if e[3] == "tcdb_family_transports_chemical_role"}
    assert role == {"tcdb:1.A.1.5": "inherited", "tcdb:1.A.1.5.2": "most_specific"}


def test_substrate_breadth(cache_root, strain_dir):
    a = _make_orchestrator(cache_root, strain_dir)
    props = {n[2]["tcdb_id"]: n[2] for n in a.get_nodes() if n[1] == "tcdb family"}
    assert props["1.A.1"]["substrate_breadth"] == "lumping"
    assert props["1.A.1.5"]["substrate_breadth"] == "specific"
    assert props["1"]["substrate_breadth"] == props["1.A"]["substrate_breadth"] == "not_assessed"


def test_missing_lumping_pin_fails(cache_root, strain_dir):
    path = cache_root / "tcdb" / "tcdb_pruned.json"
    pruned = json.loads(path.read_text())
    del pruned["lumping_pin"]
    path.write_text(json.dumps(pruned))
    with pytest.raises(RuntimeError, match="lumping_pin"):
        _make_orchestrator(cache_root, strain_dir)
```

  Read `_make_orchestrator` first: if it doesn't call `download_data()`, call it in each test before `get_nodes()`.

- [ ] **Step 2: Run** `uv run pytest tests/test_tcdb_adapter.py -q` → FAIL.

- [ ] **Step 3: Implement.**
  - `download_data`: also read `self._chemical_roles = pruned.get("chemical_roles", {})`, `self._subtree_roles = pruned.get("subtree_roles", {})`, and `self._lumping_pin = pruned.get("lumping_pin")`; if the pin is `None`, raise `RuntimeError("tcdb_pruned.json has no lumping_pin — re-run prepare_data step 6")`.
  - `_compute_substrate_depth(self, subtree_map)`: today's body with `self._subtree_substrates` replaced by `subtree_map`; call it as `self._compute_substrate_depth(self._subtree_substrates)` for substrates and `self._compute_substrate_depth(self._subtree_roles)` for roles.
  - `get_nodes`: add `"substrate_breadth": substrate_breadth(tcdb_id, level, self._lumping_pin)` (and `VOCAB.check("TcdbFamily", "substrate_breadth", …)`) to TcdbFamily props; then yield role nodes:

```python
        for role_id, r in sorted(self._chemical_roles.items()):
            yield (role_id, "chemical role", {
                "name": _clean_str(r["name"]), "level": int(r["level"]),
                "level_kind": VOCAB.check("ChemicalRole", "level_kind", r["level_kind"]),
                "all_identifiers": [_clean_str(i) for i in r["all_identifiers"]]})
```

  - `get_edges` after the substrate edges:

```python
        for role_id, r in sorted(self._chemical_roles.items()):
            for parent in r["parents"]:
                yield (f"{role_id}-isa-{parent}", role_id, parent, "chemical_role_is_a_chemical_role", {})
        role_depth = self._compute_substrate_depth(self._subtree_roles)
        for tcdb_id, roles in sorted(self._subtree_roles.items()):
            for role_id in roles:
                depth = "most_specific" if role_id in role_depth.get(tcdb_id, set()) else "inherited"
                yield (f"t2r:{tcdb_id}:{role_id}", _tcdb_node_id(tcdb_id), role_id,
                       "tcdb_family_transports_chemical_role",
                       {"substrate_depth": VOCAB.check("Tcdb_family_transports_chemical_role",
                                                       "substrate_depth", depth)})
```

  Role ids are used verbatim (never `normalize_curie`, which returns `None` for `chebi.role:`).

- [ ] **Step 4: Run** `uv run pytest tests/test_tcdb_adapter.py -q` → PASS, except `test_transport_substrate_resolution_threshold_is_inlined`, which Task 19 rewrites.

- [ ] **Step 5: Commit.**

```bash
git add multiomics_kg/adapters/tcdb_adapter.py tests/test_tcdb_adapter.py
git commit -m "feat(tcdb-adapter): ChemicalRole hierarchy + role edges, shared depth rule, substrate_breadth"
```

### Task 18: Assay `stereo_assignment`, `MetaboliteClassAdapter`, `create_knowledge_graph.py`

**Files:**
- Create: `multiomics_kg/adapters/metabolite_class_adapter.py`, `tests/test_metabolite_class_adapter.py`
- Modify: `multiomics_kg/adapters/metabolite_assay_adapter.py`, `create_knowledge_graph.py`, `tests/test_metabolite_assay_adapter.py`

**Interfaces:**
- Consumes: `brite_br08xxx.json` (Task 13); `MetabolismAdapter.all_metabolite_ids()` (Task 16); `classes.LEVEL_KINDS` (Task 8).
- Produces: `MetaboliteClassAdapter(cache_root: Path, metabolite_node_ids: set[str], test_mode: bool = False)` with `download_data()`, `get_nodes()`, `get_edges()`; `bc._output_directory / "chebi_release.txt"`.

- [ ] **Step 1: Failing tests.** `tests/test_metabolite_class_adapter.py`:

```python
import json

from multiomics_kg.adapters.metabolite_class_adapter import MetaboliteClassAdapter


def _write(tmp_path):
    kegg = tmp_path / "kegg"
    kegg.mkdir()
    (kegg / "brite_br08001.json").write_text(json.dumps({
        "tree_code": "br08001", "tree": "compounds_with_biological_roles",
        "nodes": [{"id": "kegg.brite:br08001.A1", "parent": None, "name": "Peptides", "level": 0},
                  {"id": "kegg.brite:br08001.A1.B1", "parent": "kegg.brite:br08001.A1",
                   "name": "Amino acids <i>x</i>", "level": 1}],
        "members": [{"metabolite_id": "kegg.compound:C00041", "class_id": "kegg.brite:br08001.A1.B1"},
                    {"metabolite_id": "kegg.compound:C99999", "class_id": "kegg.brite:br08001.A1.B1"}]}))
    return tmp_path


def test_nodes_edges_and_dangling_proofing(tmp_path):
    a = MetaboliteClassAdapter(cache_root=_write(tmp_path), metabolite_node_ids={"kegg.compound:C00041"})
    a.download_data()
    nodes = {n[0]: n for n in a.get_nodes()}
    assert nodes["kegg.brite:br08001.A1.B1"][1] == "metabolite class"
    assert nodes["kegg.brite:br08001.A1.B1"][2]["level_kind"] == "kegg_compound_subclass"
    edges = list(a.get_edges())
    member = [e for e in edges if e[3] == "metabolite_in_metabolite_class"]
    assert [(e[1], e[4]) for e in member] == [("kegg.compound:C00041", {"sources": ["kegg"], "evidence": "curated"})]
    assert any(e[3] == "metabolite_class_is_a_metabolite_class" and e[2] == "kegg.brite:br08001.A1" for e in edges)
```

  In `tests/test_metabolite_assay_adapter.py`, extend the existing resolved-CSV fixture with a `stereo_assignment` column and assert the value reaches both edge types; add `test_resolved_row_without_stereo_assignment_raises` (`pytest.raises(ValueError, match="stereo_assignment")`).

- [ ] **Step 2: Run** → FAIL.

- [ ] **Step 3: Implement `metabolite_class_adapter.py`.**

```python
"""MetaboliteClass nodes from KEGG compound BRITE (pruned in step 6; spec 2026-10-13 §3.7, §8.4).

Global data, so no Multi* wrapper. Membership edges are dangling-proofed against the metabolite
ids the metabolism adapter actually emitted (BRITE known_ko_ids precedent)."""
from __future__ import annotations

import json
import logging
from pathlib import Path

from multiomics_kg.utils.controlled_vocab import VOCAB
from multiomics_kg.utils.metabolite_identity.classes import COMPOUND_BRITE_TREES, LEVEL_KINDS, _clean as _clean_brite

logger = logging.getLogger(__name__)


def _clean_str(value: str) -> str:
    return value.replace("'", "^").replace("|", "") if isinstance(value, str) else value


class MetaboliteClassAdapter:
    def __init__(self, cache_root: Path, metabolite_node_ids: set[str], test_mode: bool = False):
        self.cache_root = Path(cache_root)
        self.metabolite_node_ids = set(metabolite_node_ids)
        self.test_mode = test_mode
        self._trees: dict[str, dict] = {}

    def download_data(self, cache: bool = True) -> None:
        for code in COMPOUND_BRITE_TREES:
            p = self.cache_root / "kegg" / f"brite_{code}.json"
            if p.exists():
                self._trees[code] = json.loads(p.read_text())
        if not self._trees:
            raise FileNotFoundError(f"no brite_br08xxx.json under {self.cache_root / 'kegg'}; run step 6")

    def get_nodes(self):
        for code, tree in sorted(self._trees.items()):
            for n in tree["nodes"]:
                yield (n["id"], "metabolite class", {
                    "name": _clean_str(_clean_brite(n["name"])), "tree": tree["tree"], "tree_code": code,
                    "level": int(n["level"]),
                    "level_kind": VOCAB.check("MetaboliteClass", "level_kind", LEVEL_KINDS[n["level"]])})

    def get_edges(self):
        skipped = 0
        for code, tree in sorted(self._trees.items()):
            for n in tree["nodes"]:
                if n["parent"]:
                    yield (f"{n['id']}-isa-{n['parent']}", n["id"], n["parent"],
                           "metabolite_class_is_a_metabolite_class", {})
            for m in tree["members"]:
                if m["metabolite_id"] not in self.metabolite_node_ids:
                    skipped += 1
                    continue
                yield (f"m2c:{m['metabolite_id']}:{m['class_id']}", m["metabolite_id"], m["class_id"],
                       "metabolite_in_metabolite_class", {"sources": ["kegg"], "evidence": "curated"})
        if skipped:
            logger.info(f"metabolite classes: skipped {skipped} memberships to unemitted metabolites")
```

  (Test mode is not capped here: ≈ 431 nodes; the metabolism cap is handled by dangling-proofing.)

- [ ] **Step 4: Assay adapter.** In `_emit_numeric_edges` and `_emit_boolean_edges`, read `sa = (row.get("stereo_assignment") or "").strip()`; if `primary` is set and `sa` is empty, raise `ValueError(f"{csv_path}: row {row_idx} resolved to {primary} without stereo_assignment — re-run prepare_data step 6")`; add `"stereo_assignment": VOCAB.check("Assay_quantifies_metabolite" | "Assay_flags_metabolite", "stereo_assignment", sa)` to the edge props. Fix the module docstring (lines 13–16) and the missing-file warning: "written by step 6 (build_kegg_metabolism_xrefs)" / "run prepare_data.sh --steps 6".

- [ ] **Step 5: `create_knowledge_graph.py`.** Right after the `MultiMetabolismAdapter` write (lines ~129–134):

```python
    # MetaboliteClass (KEGG compound BRITE) — after the metabolism adapter, so memberships are
    # dangling-proofed against the metabolite ids it actually emitted.
    from multiomics_kg.adapters.metabolite_class_adapter import MetaboliteClassAdapter
    metabolite_class_adapter = MetaboliteClassAdapter(
        cache_root=Path("cache/data"), metabolite_node_ids=metabolism_adapter.all_metabolite_ids(),
        test_mode=TEST_MODE)
    metabolite_class_adapter.download_data()
    bc.write_nodes(metabolite_class_adapter.get_nodes())
    bc.write_edges(metabolite_class_adapter.get_edges())
```

  Beside the vocabulary-hash write (line ~447):

```python
    chebi_cache = Path("cache/data/chebi/chebi_metabolite_xrefs.json")
    if chebi_cache.exists():
        (Path(bc._output_directory) / "chebi_release.txt").write_text(json.loads(chebi_cache.read_text())["release"])
```

  Fix the stale `# 4 nodes` DataSource comment at line ~95 → "one DataSource node per logical / reference source (11)".

- [ ] **Step 6: Run** `uv run pytest tests/test_metabolite_class_adapter.py tests/test_metabolite_assay_adapter.py -q` → PASS; then `uv run python create_knowledge_graph.py --test` → completes; `grep -c "" <output>/MetaboliteClass-part000.csv` (BioCypher's file naming; `ls` the output dir) > 0.

- [ ] **Step 7: Commit.**

```bash
git add multiomics_kg/adapters/metabolite_class_adapter.py multiomics_kg/adapters/metabolite_assay_adapter.py create_knowledge_graph.py tests/test_metabolite_class_adapter.py tests/test_metabolite_assay_adapter.py
git commit -m "feat(adapters): MetaboliteClass adapter, assay stereo_assignment, ChEBI release into the build"
```

---

## Phase 4 — Post-import

### Task 19: Rollups, indexes, `transport_substrate_resolution`, `Schema_info.chebi_release`

**Files:**
- Modify: `scripts/post-import.sh`, `scripts/post-import.cypher` (identical Cypher), `tests/test_tcdb_adapter.py`, `tests/kg_validity/test_schema_info.py`
- Create: `tests/test_post_import_metabolite_identity.py`

- [ ] **Step 1: Failing static tests.** In `tests/test_tcdb_adapter.py` replace `test_transport_substrate_resolution_threshold_is_inlined` with:

```python
def test_transport_substrate_resolution_reads_the_pinned_set():
    import pathlib
    for p in ("scripts/post-import.cypher", "scripts/post-import.sh"):
        text = pathlib.Path(p).read_text()
        assert "t.substrate_breadth = 'lumping'" in text, p
        assert "metabolite_count, 0) >= 50" not in text, f"{p}: threshold must not be re-derived"
        assert "R3 exception" in text, f"{p}: cite the documented R3 exception"
```

  `tests/test_post_import_metabolite_identity.py` asserts each of these literals appears in **both** files: `metaboliteClassFullText`, `chemicalRoleFullText`, `metabolite_structure_kind_idx`, `metabolite_class_id_idx`, `chemical_role_id_idx`, `m.identifier_tokens`, `[m.name, m.formula, m.identifier_tokens]`, `c.direct_metabolite_count`, `r.direct_transporter_count`, `m.kegg_compound_categories`, `s.chebi_release`. In `tests/kg_validity/test_schema_info.py` add `test_chebi_release_stamped` (`Schema_info.chebi_release` equals the committed ChEBI cache's `release`).

- [ ] **Step 2: Run** `uv run pytest tests/test_tcdb_adapter.py tests/test_post_import_metabolite_identity.py -q` → FAIL.

- [ ] **Step 3: Group 1 (indexes),** next to the Metabolite indexes (`post-import.sh` lines ~228–234):

```cypher
CREATE INDEX metabolite_structure_kind_idx IF NOT EXISTS FOR (m:Metabolite) ON (m.structure_kind);
DROP INDEX metaboliteFullText IF EXISTS;
// identifier_tokens = bare local parts of all_identifiers (set in group 3): the standard
// analyser keeps "kegg.compound:C01353" as ONE token (UAX#29: ':' and '.' between letters /
// digits don't break), so a bare "C01353" query would miss the CURIE form.
CREATE FULLTEXT INDEX metaboliteFullText IF NOT EXISTS FOR (m:Metabolite) ON EACH [m.name, m.formula, m.identifier_tokens];
CREATE INDEX metabolite_class_id_idx IF NOT EXISTS FOR (c:MetaboliteClass) ON (c.id);
CREATE INDEX metabolite_class_tree_idx IF NOT EXISTS FOR (c:MetaboliteClass) ON (c.tree);
CREATE INDEX metabolite_class_level_idx IF NOT EXISTS FOR (c:MetaboliteClass) ON (c.level);
CREATE INDEX metabolite_class_level_kind_idx IF NOT EXISTS FOR (c:MetaboliteClass) ON (c.level_kind);
CREATE FULLTEXT INDEX metaboliteClassFullText IF NOT EXISTS FOR (c:MetaboliteClass) ON EACH [c.name];
CREATE INDEX chemical_role_id_idx IF NOT EXISTS FOR (r:ChemicalRole) ON (r.id);
CREATE INDEX chemical_role_level_idx IF NOT EXISTS FOR (r:ChemicalRole) ON (r.level);
CREATE INDEX chemical_role_level_kind_idx IF NOT EXISTS FOR (r:ChemicalRole) ON (r.level_kind);
CREATE FULLTEXT INDEX chemicalRoleFullText IF NOT EXISTS FOR (r:ChemicalRole) ON EACH [r.name];
```

  (`all_identifiers` isn't indexed directly: the bare tokens cover id search, and CURIE look-ups use `IN m.all_identifiers` with an exact match.)

- [ ] **Step 4: Group 3.** Replace the breadth line in the Gene transport arm (`post-import.sh` ~1504–1528) and its comment:

```cypher
// transport_substrate_resolution reads the PINNED lumping set (TcdbFamily.substrate_breadth,
// spec 2026-10-13 §3.4a) — a documented R3 exception (vocabulary-contract.md): the threshold
// level >= 2 AND metabolite_count >= 50 no longer reproduces the set once metabolite identity
// merges substrates (3.A.3: 47 molecules, still lumping). Sunset with the deferred TCDB spec.
         collect(DISTINCT (t.substrate_breadth = 'lumping')) AS breadth
```

  (the `CASE … any(x IN breadth WHERE x = false) …` logic is unchanged.)

  Add, **after** the measured `Organism_has_metabolite` MERGE and the `organism_count` / `organism_names` recompute (`post-import.sh` ~1952–1959), so organism counts see all three arms:

```cypher
// Metabolite.identifier_tokens — search-only bare ids for metaboliteFullText
CALL {
  MATCH (m:Metabolite) WHERE m.all_identifiers IS NOT NULL
  SET m.identifier_tokens = [x IN m.all_identifiers | split(x, ':')[1]]
} IN TRANSACTIONS OF 1000 ROWS;

// MetaboliteClass rollups (subtree / direct / organisms)
MATCH (c:MetaboliteClass)
CALL {
  WITH c
  OPTIONAL MATCH (c)<-[:Metabolite_class_is_a_metabolite_class*0..]-(:MetaboliteClass)
                   <-[:Metabolite_in_metabolite_class]-(m:Metabolite)
  WITH c, collect(DISTINCT m) AS ms
  OPTIONAL MATCH (c)<-[:Metabolite_in_metabolite_class]-(dm:Metabolite)
  WITH c, ms, count(DISTINCT dm) AS direct
  UNWIND (CASE WHEN size(ms) = 0 THEN [null] ELSE ms END) AS m
  OPTIONAL MATCH (o:OrganismTaxon)-[:Organism_has_metabolite]->(m)
  WITH c, size(ms) AS mc, direct, count(DISTINCT o) AS oc
  SET c.metabolite_count = mc, c.direct_metabolite_count = direct, c.organism_count = oc
} IN TRANSACTIONS OF 100 ROWS;

// Metabolite.kegg_compound_categories — br08001 level-0 names, sorted, sparse
MATCH (m:Metabolite)
CALL {
  WITH m
  OPTIONAL MATCH (m)-[:Metabolite_in_metabolite_class]->(:MetaboliteClass)
                 -[:Metabolite_class_is_a_metabolite_class*0..]->(top:MetaboliteClass {tree_code: 'br08001', level: 0})
  WITH m, top.name AS n ORDER BY n
  WITH m, collect(DISTINCT n) AS cats
  SET m.kegg_compound_categories = CASE WHEN size(cats) = 0 THEN null ELSE cats END
} IN TRANSACTIONS OF 1000 ROWS;

// ChemicalRole family counts: subtree over the is_a DAG, DISTINCT families; never propagated
// to Gene / OrganismTaxon / Metabolite rollups (spec §3.4 "not propagated").
MATCH (r:ChemicalRole)
CALL {
  WITH r
  OPTIONAL MATCH (r)<-[:Chemical_role_is_a_chemical_role*0..]-(d:ChemicalRole)
  WITH r, collect(DISTINCT d) AS descs
  UNWIND descs AS d
  OPTIONAL MATCH (t:TcdbFamily)-[e:Tcdb_family_transports_chemical_role]->(d)
  WITH r,
       count(DISTINCT t) AS any_t,
       count(DISTINCT CASE WHEN e.substrate_depth = 'most_specific' THEN t END) AS ms_t,
       count(DISTINCT CASE WHEN d = r THEN t END) AS d_any_t,
       count(DISTINCT CASE WHEN d = r AND e.substrate_depth = 'most_specific' THEN t END) AS d_ms_t
  SET r.tcdb_family_count = any_t, r.transporter_count = ms_t,
      r.direct_tcdb_family_count = d_any_t, r.direct_transporter_count = d_ms_t
} IN TRANSACTIONS OF 100 ROWS;
```

  Check every relationship name against the BioCypher output label casing in the `import` script (`Metabolite_in_metabolite_class`, `Metabolite_class_is_a_metabolite_class`, `Chemical_role_is_a_chemical_role`, `Tcdb_family_transports_chemical_role`) after Task 18's `--test` build: `grep -o "relationships=[A-Za-z_]*" <output>/neo4j-admin-import-call.sh | sort -u`.

- [ ] **Step 5: Group 4.** Beside `KG_VOCAB_HASH`, read `/data/build2neo/chebi_release.txt` (same flat-or-timestamped `ls -t` glob, same warning when absent) into `KG_CHEBI_RELEASE`; add `-P "chebi_release => '${KG_CHEBI_RELEASE}'"` and `s.chebi_release = CASE WHEN coalesce($chebi_release, '') = '' THEN null ELSE $chebi_release END`.

- [ ] **Step 6: Mirror** every change into `scripts/post-import.cypher` (same Cypher; its Group-4 equivalent uses the same `$chebi_release` parameter).

- [ ] **Step 7: Run** `uv run pytest tests/test_tcdb_adapter.py tests/test_post_import_metabolite_identity.py tests/test_annotation_quality_buckets.py -q` → PASS.

- [ ] **Step 8: Commit.**

```bash
git add scripts/post-import.sh scripts/post-import.cypher tests/test_tcdb_adapter.py tests/test_post_import_metabolite_identity.py tests/kg_validity/test_schema_info.py
git commit -m "feat(post-import): MetaboliteClass/ChemicalRole rollups + indexes, pinned lumping, chebi_release"
```

---

## Phase 5 — Skills and docs

### Task 20: Skills and `validate_paperconfig.py`

**Files:**
- Modify: `scripts/validate_paperconfig.py`, `tests/test_validate_paperconfig_metabolite_assays.py`, `.claude/skills/paperconfig/SKILL.md`, `.claude/skills/add-a-strain/SKILL.md`, `.claude/skills/cypher-queries/SKILL.md`, `.claude/skills/omics-edge-snapshot/{SKILL.md,omics_edge_snapshot.py}`

- [ ] **Step 1: Validator tests first** (append to `tests/test_validate_paperconfig_metabolite_assays.py`): an alias yaml with `{id, stereo_assignment, reason}` validates; `stereo_assignment: maybe` errors; two entries with one `name_key` and different ids error (`Glutathione` / `glutathione`); `CO2` / `Co2` with different ids pass. Run → FAIL.

- [ ] **Step 2: Implement** in `validate_paperconfig.py` (after the `aliases_file` existence check, ~line 1306): when the file exists, call `multiomics_kg.utils.metabolite_identity.paper.load_paper_overrides(Path(config_path), cfg)` inside `try` and append the `ValueError` text to `errors`. The "legacy string form on a stereo-variant target" rule needs the KG's stereo families, so it stays in step 6 (which fails with the explicit message); the validator emits a **warning** for every legacy string entry: "prefer {id, stereo_assignment, reason}". Run → PASS.

- [ ] **Step 3: `paperconfig/SKILL.md`.** Rewrite the "ID resolution priority" section to the spec §3.6 order (paper alias override → direct id → dominant-form lock → name match → dominant-form rule; ambiguous names fail step 6; unmatched names are unresolved); document the `{id, stereo_assignment, reason}` alias form (legacy string valid only for a target with no stereo variant), the stereo-qualified name rule (descriptors list from spec §3.6, bare `R-` / `S-` excluded), case-insensitive keys except formula-like names, per-paper overrides (union over the paper's alias files); replace "step 7 just looks up …" with "step 6 writes `<stem>_resolved.csv` (incl. `stereo_assignment`)"; perf loop `--steps 6 7 --force` → `--steps 6 --force`.

- [ ] **Step 4: `add-a-strain/SKILL.md`.** Step-output table: step 6 also writes the paper `_resolved.csv` files; step 7 is retired (`--steps … 7` still accepted as a no-op). Add a post-onboarding check: read `cache/data/metabolomics/metabolite_identity_report.json` → `lumping.lumping_set_added`, every `stale_*` list, and `rekeys` added since the last commit (`git diff` the report).

- [ ] **Step 5: `cypher-queries/SKILL.md`.** Add templates:

```cypher
// Metabolite by any identifier (CURIE form)
MATCH (m:Metabolite) WHERE m.id = $curie OR $curie IN m.all_identifiers RETURN m.id, m.name;
// Families that transport something TCDB calls a siderophore (incl. sub-roles)
MATCH (r:ChemicalRole {name: 'siderophore'})<-[:Chemical_role_is_a_chemical_role*0..]-(sub:ChemicalRole)
      <-[e:Tcdb_family_transports_chemical_role {substrate_depth: 'most_specific'}]-(t:TcdbFamily)
RETURN DISTINCT t.tcdb_id, t.name, collect(DISTINCT sub.name) AS roles;
// Carbohydrates an organism has (any arm)
MATCH (o:OrganismTaxon {preferred_name: $org})-[:Organism_has_metabolite]->(m:Metabolite)
WHERE 'Carbohydrates' IN coalesce(m.kegg_compound_categories, []) RETURN m.name ORDER BY m.name;
// Stereo hop: the unspecified form and the mirror image of L-Ala
MATCH (m:Metabolite {id: 'kegg.compound:C00041'})
OPTIONAL MATCH (m)-[:Metabolite_is_stereo_form_of]->(u) OPTIONAL MATCH (m)-[:Metabolite_is_enantiomer_of]-(en)
RETURN u.name, en.name;
```

  and fix the stale hint that tells users to filter `Tcdb_family_transports_metabolite` by `level_kind = 'tc_specificity'` → `substrate_depth = 'most_specific'`.

- [ ] **Step 6: `omics-edge-snapshot`.** In `omics_edge_snapshot.py` (lines ~210–248): `METABOLISM_NODE_COUNTS += ["MetaboliteClass", "ChemicalRole"]`; `METABOLISM_EDGE_COUNTS += ["Tcdb_family_transports_metabolite", "Tcdb_family_transports_chemical_role", "Metabolite_is_stereo_form_of", "Metabolite_is_enantiomer_of", "Metabolite_in_metabolite_class", "Assay_quantifies_metabolite", "Assay_flags_metabolite"]`. `compare_snapshots` must tolerate keys absent from an old snapshot (treat as 0). In SKILL.md add a "Metabolite identity rebuild (2026-10)" note: the metabolism block shows expected drops (`Metabolite`, `Reaction_has_metabolite`, `Organism_has_metabolite`, `Tcdb_family_transports_metabolite`); compare against the expected deltas in `docs/kg-changes/metabolite-identity-and-classes.md`, not the regression flags.

- [ ] **Step 7: Run** `uv run pytest tests/test_validate_paperconfig_metabolite_assays.py -q && uv run python scripts/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/paperconfig.yaml"` → PASS, no errors.

- [ ] **Step 8: Commit.**

```bash
git add scripts/validate_paperconfig.py tests/test_validate_paperconfig_metabolite_assays.py .claude/skills/paperconfig/SKILL.md .claude/skills/add-a-strain/SKILL.md .claude/skills/cypher-queries/SKILL.md .claude/skills/omics-edge-snapshot/SKILL.md .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py
git commit -m "docs(skills): metabolite identity in paperconfig / add-a-strain / cypher / snapshot; validator"
```

### Task 21: Repository docs

**Files:**
- Modify: `CLAUDE.md`, `CHANGELOG.md`, `plans/backlog.md`, `docs/kg-changes/vocabulary-contract.md`, `docs/kg-changes/tcdb-two-source-upgrade.md`, `docs/kg-changes/metabolomics-extension.md`

- [ ] **Step 1: `CLAUDE.md`** (spec §14):
  - Metabolite bullet: canonical identity (KEGG → primary ChEBI → MNX), the new properties with their sparsity, `all_identifiers` CURIE-only vs `Gene.all_identifiers` bare, `recurated_from_ids` (no identity meaning), `structure_kind`, `name` source.
  - New MetaboliteClass bullet (8 trees, `kegg.brite:br08xxx…` ids, `level_kind` names, partial coverage, separate from BriteCategory).
  - New ChemicalRole bullet ("this family transports something TCDB calls X"; `chebi.role:` prefix; hierarchy; not propagated).
  - TcdbFamily bullet: `substrate_breadth` and the pinned set replace the "query the threshold directly" text; the R3 exception and its sunset.
  - Stereo edges + `stereo_assignment`.
  - Label / relationship lists: add `MetaboliteClass`, `ChemicalRole`, `Metabolite_is_stereo_form_of`, `Metabolite_is_enantiomer_of`, `Metabolite_in_metabolite_class`, `Metabolite_class_is_a_metabolite_class`, `Chemical_role_is_a_chemical_role`, `Tcdb_family_transports_chemical_role`.
  - Step 6 description: now also ChEBI (stage A), roles, identity, classes, lock, report, the paper `_resolved.csv` files; step 7 entry: retired (no-op alias). Default `STEPS` string.
  - Data Locations: `cache/data/chebi/` (`chebi_metabolite_xrefs.json` committed, `raw/` gitignored), `cache/data/kegg/brite_br08xxx.json`, `cache/data/metabolomics/dominant_stereo_forms.json` (lock: edit by hand only), `metabolite_identity_report.json`, `config/metabolite_identity_overrides.yaml`.
  - Post-import indexes list: the new ones; correct `metaboliteFullText` (was documented as name + formula, live was name only; now name + formula + `identifier_tokens`).
  - DataSource count 9 → 11 and the `scope` / `provenance` values actually used (`gene_level` / `organism_restricted` / `global`; `download` / `tool_run`).
  - Trust-vocabulary list: add `identity_basis`, `structure_kind`, `substrate_breadth`, `stereo_assignment`; update the count.

- [ ] **Step 2: `CHANGELOG.md`** under `## [Unreleased]`:
  - `### Breaking`: metabolite id remaps (merged / re-keyed / secondary→primary; see `metabolite_id_remap.tsv`); `chebi_id` values change on most nodes (now primary, entity-level); 57 role terms removed from `Metabolite`; the 4 `mnx:` placeholders and the measurement-only nodes removed; paper metabolite re-resolutions (see `metabolite_assay_repoints.tsv`); `name` values on `chebi:` nodes; `TcdbFamily.metabolite_count` loses role entries; the threshold query no longer reproduces the lumping set; `metaboliteFullText` fields.
  - `### Added`: `MetaboliteClass`; identity properties; stereo-form + enantiomer edges; `stereo_assignment`; `ChemicalRole` + `Tcdb_family_transports_chemical_role`; `substrate_breadth` (pinned, no value change); `data_source:kegg` / `data_source:chebi`; `Schema_info.chebi_release`.
  - `### Changed`: step 7 folded into step 6.
  - `### Highlights`: "Ask for every carbohydrate or peptide an organism makes, takes up or was measured with."

- [ ] **Step 3: `docs/kg-changes/vocabulary-contract.md` and `tcdb-two-source-upgrade.md` §4:** the pinned lumping set, the R3 exception (why, scope, sunset with the deferred TCDB spec), and "filter on `substrate_breadth`, don't re-derive".

- [ ] **Step 4: `docs/kg-changes/metabolomics-extension.md`:** step 7 → step 6 writes `_resolved.csv`; the `stereo_assignment` column.

- [ ] **Step 5: `plans/backlog.md`:** close "Metabolite equivalence across namespaces and stereo/charge variants" (P1) and "Chemical categories for `Metabolite`" (P2) with a pointer to this plan; file the deferred TCDB spec item (plausibility, `substrate_depth` semantics, per-organism reach, TCDB name hygiene, lumping-threshold re-tune + retiring the R3 exception).

- [ ] **Step 6: Commit.**

```bash
git add CLAUDE.md CHANGELOG.md plans/backlog.md docs/kg-changes/vocabulary-contract.md docs/kg-changes/tcdb-two-source-upgrade.md docs/kg-changes/metabolomics-extension.md
git commit -m "docs: metabolite identity + classes — CLAUDE.md, CHANGELOG, contracts, backlog"
```

---

## Phase 6 — Validation

The live graph at `bolt://localhost:7687` is the pre-change build until Task 24. **Capture Task 22's baselines before anyone rebuilds the Docker clone**; if a rebuild for other work is imminent, run Task 22 first, out of order: it doesn't depend on any code in this plan except Task 20's snapshot-script update, and a baseline taken without that update is still valid (the new types read as absent = 0).

### Task 22: Baselines (before the rebuild)

**Files:**
- Create: `tests/kg_validity/fixtures/metabolite_elements_pre_identity.json`, `.claude/skills/omics-edge-snapshot/snapshots/pre_metabolite_identity.json` (written by the script)
- Modify: `tests/kg_validity/annotation_state_baseline.json` (written by the capture script)

- [ ] **Step 1: Edge snapshot.** Invoke the `omics-edge-snapshot` skill to save a baseline: `uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --save pre_metabolite_identity`.
- [ ] **Step 2: Annotation state.** `uv run python tests/kg_validity/capture_annotation_state.py --save`.
- [ ] **Step 3: Pre-change `elements` export** (for the "no node loses `elements`" test):

```bash
uv run python - <<'EOF'
import json
from neo4j import GraphDatabase
d = GraphDatabase.driver("bolt://localhost:7687")
with d.session() as s:
    rows = s.run("MATCH (m:Metabolite) WHERE m.elements IS NOT NULL RETURN m.id AS id, m.elements AS el").data()
json.dump({r["id"]: sorted(r["el"]) for r in rows},
          open("tests/kg_validity/fixtures/metabolite_elements_pre_identity.json", "w"), indent=1, sort_keys=True)
print(len(rows))
EOF
```

  Expected: 2,982.
- [ ] **Step 4: Record the other pre-change counts** in the post-rebuild doc draft (Task 25): `transport_substrate_resolution` (resolved 28,854 / family_inferred 1,693), the metabolism counts and the assay counts from the plan's measured-baseline table. Re-run those queries now: if any differ from the table, the live graph changed since 2026-10-14; use today's numbers.
- [ ] **Step 5: Commit** `git add tests/kg_validity/fixtures/metabolite_elements_pre_identity.json tests/kg_validity/annotation_state_baseline.json .claude/skills/omics-edge-snapshot/snapshots/pre_metabolite_identity.json && git commit -m "test(kg): pre-change baselines for the metabolite identity rebuild"`.

### Task 23: KG validity tests (`tests/kg_validity/test_metabolite_identity.py`)

**Files:**
- Create: `tests/kg_validity/test_metabolite_identity.py`
- Modify: `tests/kg_validity/test_metabolism.py` (`test_metabolite_primary_ids_predominantly_kegg` stays; check its threshold still holds after Task 14's numbers)

- [ ] **Step 1: Write the tests** (`pytestmark = pytest.mark.kg`; use the `run_query` fixture). One test per spec §11 KG-validity bullet:

```python
import json
import pathlib

import pytest

pytestmark = pytest.mark.kg
REMAP = pathlib.Path("docs/kg-changes/metabolite-identity/metabolite_id_remap.tsv")


def _one(run_query, q, **kw):
    return run_query(q, **kw)[0]


def test_no_duplicate_identifiers(run_query):
    for prop in ("kegg_compound_id", "chebi_id", "mnxm_id"):
        dup = run_query(f"MATCH (m:Metabolite) WHERE m.{prop} IS NOT NULL WITH m.{prop} AS v, count(*) AS n "
                        f"WHERE n > 1 RETURN v, n LIMIT 5")
        assert not dup, (prop, dup)
    dup = run_query("""
        CALL { MATCH (m:Metabolite) UNWIND m.all_identifiers AS x RETURN x
               UNION ALL MATCH (r:ChemicalRole) UNWIND r.all_identifiers AS x RETURN x }
        WITH x, count(*) AS n WHERE n > 1 RETURN x, n LIMIT 5""")
    assert not dup


def test_no_inchikey_on_two_nodes(run_query):
    assert not run_query("MATCH (m:Metabolite) WHERE m.inchikey IS NOT NULL WITH m.inchikey AS k, count(*) AS n "
                         "WHERE n > 1 RETURN k LIMIT 5")


def test_chebi_keyed_nodes_chebi_id_matches(run_query):
    assert not run_query("MATCH (m:Metabolite) WHERE m.id STARTS WITH 'chebi:' AND m.chebi_id <> split(m.id, ':')[1] "
                         "RETURN m.id LIMIT 5")


def test_all_identifiers_are_curies(run_query):
    assert not run_query("MATCH (m:Metabolite) UNWIND m.all_identifiers AS x WITH m, x "
                         "WHERE NOT x =~ '(kegg\\\\.compound|chebi|mnx|hmdb|inchikey):.+' RETURN m.id, x LIMIT 5")


def test_no_role_metabolites_and_role_ids(run_query):
    role_ids = [r["id"] for r in run_query("MATCH (r:ChemicalRole) RETURN r.id AS id")]
    assert role_ids and all(i.startswith("chebi.role:") and i.split(":")[1].isdigit() for i in role_ids)
    secondary_hits = run_query("""MATCH (r:ChemicalRole) UNWIND r.all_identifiers AS x
        WITH x WHERE x STARTS WITH 'chebi:' MATCH (m:Metabolite {id: x}) RETURN m.id LIMIT 5""")
    assert not secondary_hits
    assert not run_query("MATCH (r:ChemicalRole) WHERE r.level_is_best_effort IS NOT NULL RETURN r.id LIMIT 1")


def test_old_secondary_role_id_resolves(run_query):
    assert run_query("MATCH (r:ChemicalRole) WHERE 'chebi:10004' IN r.all_identifiers RETURN r.id AS id")[0]["id"] \
        == "chebi.role:33229"


def test_role_subtree_counts_distinct_and_ge_direct(run_query):
    bad = run_query("MATCH (r:ChemicalRole) WHERE r.tcdb_family_count < r.direct_tcdb_family_count "
                    "OR r.transporter_count < r.direct_transporter_count RETURN r.id LIMIT 5")
    assert not bad
    row = _one(run_query, """MATCH (r:ChemicalRole) WHERE r.level = 0
        OPTIONAL MATCH (r)<-[:Chemical_role_is_a_chemical_role*0..]-(d)<-[:Tcdb_family_transports_chemical_role]-(t)
        WITH r, count(DISTINCT t) AS n RETURN collect(r.tcdb_family_count = n) AS ok""")
    assert all(row["ok"])


def test_roles_not_propagated(run_query):
    assert not run_query("MATCH (r:ChemicalRole)<-[]-(x) WHERE x:Gene OR x:OrganismTaxon OR x:Metabolite RETURN x LIMIT 1")


def test_no_isolated_metabolite_and_metabolomics_tag_matches_assays(run_query):
    assert not run_query("MATCH (m:Metabolite) WHERE NOT EXISTS { (m)--() } RETURN m.id LIMIT 5")
    assert not run_query("""MATCH (m:Metabolite) WHERE 'metabolomics' IN m.evidence_sources
        AND NOT EXISTS { (:MetaboliteAssay)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m) } RETURN m.id LIMIT 5""")
    assert not run_query("""MATCH (:MetaboliteAssay)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
        WHERE NOT 'metabolomics' IN m.evidence_sources RETURN DISTINCT m.id LIMIT 5""")


def test_recurated_ids_are_not_identifiers(run_query):
    assert not run_query("MATCH (m:Metabolite) UNWIND coalesce(m.recurated_from_ids, []) AS x "
                         "MATCH (n:Metabolite) WHERE x IN n.all_identifiers OR n.id = x RETURN x LIMIT 5")
    dup = run_query("MATCH (m:Metabolite) UNWIND coalesce(m.recurated_from_ids, []) AS x WITH x, count(*) AS n "
                    "WHERE n > 1 RETURN x")
    assert not dup


def test_structure_and_elements(run_query):
    assert not run_query("MATCH (m:Metabolite) WHERE m.structure_kind IS NULL RETURN m.id LIMIT 5")
    assert not run_query("MATCH (m:Metabolite) WHERE m.structure_kind <> 'defined' AND m.inchikey IS NOT NULL "
                         "RETURN m.id LIMIT 5")
    assert not run_query("MATCH (m:Metabolite) WHERE m.structure_kind <> 'defined' AND m.mass IS NOT NULL "
                         "RETURN m.id LIMIT 5")


def test_no_node_loses_elements(run_query):
    pre = json.loads(pathlib.Path("tests/kg_validity/fixtures/metabolite_elements_pre_identity.json").read_text())
    remap = {}
    for line in REMAP.read_text().splitlines()[1:]:
        old, new, basis, role = (line.split("\t") + ["", "", "", ""])[:4]
        remap[old] = new
    with_el = {r["id"] for r in run_query("MATCH (m:Metabolite) WHERE m.elements IS NOT NULL RETURN m.id AS id")}
    lost = [o for o in pre if (remap.get(o, o) or None) and remap.get(o, o) not in with_el
            and not (o in remap and remap[o] == "")]          # removed roles are not metabolites any more
    assert not lost, lost[:10]


def test_identity_basis_only_on_merged_or_rekeyed(run_query):
    bad = run_query("MATCH (m:Metabolite) WHERE m.identity_basis IS NOT NULL AND size(m.all_identifiers) < 2 "
                    "RETURN m.id LIMIT 5")
    assert not bad


def test_stereo_edges(run_query):
    assert not run_query("MATCH (a)-[e:Metabolite_is_stereo_form_of|Metabolite_is_enantiomer_of]->(b) "
                         "WHERE a = b OR e.sources <> ['chebi'] OR e.evidence <> 'curated' RETURN a.id LIMIT 5")
    assert not run_query("MATCH (a)-[:Metabolite_is_enantiomer_of]->(b) WHERE a.id > b.id RETURN a.id LIMIT 5")
    assert not run_query("MATCH (a)-[:Metabolite_is_enantiomer_of]->(b)-[:Metabolite_is_enantiomer_of]->(a) "
                         "RETURN a.id LIMIT 1")


def test_no_anomer_enantiomer_edges(run_query):
    names = run_query("MATCH (a)-[:Metabolite_is_enantiomer_of]->(b) RETURN a.name AS a, b.name AS b")
    from multiomics_kg.utils.metabolite_identity.names import strip_anomer_tokens
    assert not [r for r in names if strip_anomer_tokens(r["a"]) == strip_anomer_tokens(r["b"])]


def test_assay_edges(run_query):
    assert not run_query("MATCH ()-[e:Assay_quantifies_metabolite|Assay_flags_metabolite]->() "
                         "WHERE e.stereo_assignment IS NULL RETURN e LIMIT 1")
    assert not run_query("""MATCH (a:MetaboliteAssay)-[e:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m)
        WITH a, m, e.condition_label AS c, e.time_point AS t, count(e) AS n WHERE n > 1 RETURN a.id, m.id LIMIT 5""")


def test_metabolite_class(run_query):
    assert not run_query("MATCH (c:MetaboliteClass) WHERE c.metabolite_count < c.direct_metabolite_count RETURN c.id LIMIT 5")
    from multiomics_kg.utils.metabolite_identity.classes import BR08001_LEVEL0
    cats = {x["c"] for x in run_query("MATCH (m:Metabolite) UNWIND m.kegg_compound_categories AS c RETURN DISTINCT c")}
    assert cats and cats <= set(BR08001_LEVEL0)


def test_new_edge_sources_join_data_sources(run_query):
    bad = run_query("""MATCH ()-[e:Metabolite_in_metabolite_class|Metabolite_is_stereo_form_of|Metabolite_is_enantiomer_of]->()
        UNWIND e.sources AS s WITH DISTINCT s WHERE NOT EXISTS { MATCH (d:DataSource {id: 'data_source:' + s}) }
        RETURN s""")
    assert not bad


def test_fulltext_hits(run_query):
    for q in ("C01353", "MNXM732398", "HMDB0002878"):
        hits = run_query("CALL db.index.fulltext.queryNodes('metaboliteFullText', $q) YIELD node RETURN node.id LIMIT 3", q=q)
        assert hits, q
    assert run_query("CALL db.index.fulltext.queryNodes('chemicalRoleFullText', 'siderophore') YIELD node RETURN node.id")


def test_lumping_pin(run_query):
    pin = json.loads(pathlib.Path("cache/data/tcdb/tcdb_pruned.json").read_text())
    expected = set(pin["lumping_pin"]["lumping_tcdb_ids"]) & set(pin["kept_tcdb_ids"])
    live = {r["id"] for r in run_query("MATCH (t:TcdbFamily {substrate_breadth: 'lumping'}) RETURN t.tcdb_id AS id")}
    assert live == expected
    assert not run_query("MATCH (t:TcdbFamily) WHERE t.level < 2 AND t.substrate_breadth <> 'not_assessed' RETURN t.tcdb_id LIMIT 3")


def test_transport_substrate_resolution_unchanged(run_query):
    dist = {r["v"]: r["n"] for r in run_query("MATCH (g:Gene) WHERE g.transport_substrate_resolution IS NOT NULL "
                                              "RETURN g.transport_substrate_resolution AS v, count(*) AS n")}
    assert dist == {"resolved": 28854, "family_inferred": 1693}      # Task 22 baseline; update if it moved


def test_every_rekeyed_node_has_kegg_enrichment(run_query):
    bad = run_query("""MATCH (m:Metabolite) WHERE 'chebi_kegg_xref' IN coalesce(m.identity_basis, [])
        AND m.id STARTS WITH 'kegg.compound:' AND m.structure_kind = 'defined' AND m.formula IS NULL RETURN m.id LIMIT 5""")
    assert not bad
```

  The nitrate (`chebi:14654` → `kegg.compound:C00244`) and L-Ala-measured spot checks belong to the snapshot regeneration (Task 24) and the post-rebuild doc (Task 25).

- [ ] **Step 2: Run against the current (pre-change) graph** `uv run pytest tests/kg_validity/test_metabolite_identity.py -q` → most tests FAIL (expected: the graph is pre-change). This confirms they test something.

- [ ] **Step 3: Commit** `git add tests/kg_validity/test_metabolite_identity.py && git commit -m "test(kg): metabolite identity, roles, classes, stereo, lumping invariants"`.

### Task 24: Push, Docker rebuild in the other clone, KG tests, comparisons

- [ ] **Step 1: Full unit suite** `uv run pytest -m "not slow and not kg" -q` → PASS. Also `uv run python create_knowledge_graph.py --test` → completes.
- [ ] **Step 2: Generate the release TSVs** now (the KG test `test_no_node_loses_elements` reads the remap): `uv run python scripts/metabolite_identity_tsvs.py` (write it in this step: read `cache/data/metabolomics/metabolite_identity_report.json`, write `docs/kg-changes/metabolite-identity/metabolite_id_remap.tsv` with header `old_id\tnew_id\tbasis\trole_id` from `report["remap"]`, and `metabolite_assay_repoints.tsv` with header `paper\ttable\tpaper_name\told_id\tnew_id\tstereo_assignment\treason` from `report["repoints"]`, both sorted). Commit the script and both TSVs.
- [ ] **Step 3: Push** (this clone's remote is HTTPS without credentials): `git push git@github.com:wosnat/multiomics_biocypher_kg.git HEAD:main`. Ask the user before pushing if the branch isn't `main` or the push isn't a fast-forward.
- [ ] **Step 4: Rebuild** — hand off to the user: in `/home/osnat/github/multiomics_biocypher_kg`: `git pull && docker compose up -d --build`. Wait for `post-process` to exit 0 (`docker compose logs post-process | tail`). Then check `./output/import.status` and `./output/import.report` there: no skipped relationships (`pytest tests/kg_validity/test_import_report.py` run from that clone, or read the file).
- [ ] **Step 5: KG tests** (from this clone, against `bolt://localhost:7687`): `uv run pytest -m kg -q`. Expected failures to fix in the same change, not to skip: `test_snapshot.py` (it pins `chebi:10004` / `chebi:10074` as `Metabolite`) → regenerate in Step 7; anything pinning 9 DataSources → already updated (Task 15). Any other failure is a bug: debug it (superpowers:systematic-debugging), fix, re-run step 6 / rebuild as needed.
- [ ] **Step 6: Comparisons.**
  - `uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --compare pre_metabolite_identity` → **expression edges unchanged**; metabolism drops only on `Metabolite`, `Reaction_has_metabolite`, `Organism_has_metabolite`, `Tcdb_family_transports_metabolite`, matching the step-6 report (−roles −merged −dropped nodes; −392 role edges minus merge dedup); assay edge totals change only by the GSH/GSSG split and the documented re-resolutions.
  - `uv run python tests/kg_validity/capture_annotation_state.py --compare` → **0 gene movement**.
  - Pathway-level `KeggTerm` count and the zero-gene subset vs the 447 / 61 baseline (record the delta).
- [ ] **Step 7: Regenerate the snapshot fixture** `uv run python tests/kg_validity/generate_snapshot.py`, review the diff (role metabolites gone, `chebi:` → canonical ids), commit with `git commit -m "test(kg): regenerate snapshot after the metabolite identity rebuild"`.

### Task 25: Post-rebuild doc with the two TSVs and measured counts

**Files:**
- Modify: `docs/kg-changes/metabolite-identity-and-classes.md` (top section), `docs/superpowers/plans/2026-10-14-metabolite-identity-and-classes.md` (record final numbers)
- Create: `docs/kg-changes/metabolite-identity/metabolite_id_remap.tsv`, `metabolite_assay_repoints.tsv` (Task 24 Step 2)

- [ ] **Step 1: Measure on the live graph** and replace every "≈" in the doc's top section (Summary → §2b) with the measured value. Queries:

```cypher
MATCH (m:Metabolite) RETURN split(m.id, ':')[0] AS prefix, count(*) ORDER BY prefix;
MATCH (m:Metabolite) RETURN m.structure_kind, count(*), sum(CASE WHEN EXISTS { (m)<-[:Tcdb_family_transports_metabolite]-() } THEN 1 ELSE 0 END) AS with_transport;
MATCH ()-[e:Metabolite_is_stereo_form_of]->() RETURN count(e);
MATCH ()-[e:Metabolite_is_enantiomer_of]->() RETURN count(e);
MATCH (c:MetaboliteClass) RETURN c.tree, count(*);
MATCH ()-[e:Metabolite_in_metabolite_class]->() RETURN count(e);
MATCH (r:ChemicalRole) RETURN r.level, count(*) ORDER BY r.level;
MATCH ()-[e:Tcdb_family_transports_chemical_role]->() RETURN e.substrate_depth, count(*);
MATCH (k:KeggTerm {level_kind: 'pathway'}) RETURN count(k), sum(CASE WHEN coalesce(k.gene_count, 0) = 0 THEN 1 ELSE 0 END);
```

- [ ] **Step 2: Per-name stereo table** (spec §12): for each of the 36 B5 names plus every `case_fold_convergence` row in `metabolite_assay_repoints.tsv`: old id, new id, `stereo_assignment`, and `catalyst_gene_count` / `transporter_gene_count` on the new node (one Cypher `UNWIND $ids` query).
- [ ] **Step 3: Expected metabolism deltas** table (old → new, with the reason) for `Metabolite`, `Reaction_has_metabolite`, `Organism_has_metabolite`, `Tcdb_family_transports_metabolite`, assay edges, so `/omics-edge-snapshot` drops are not read as regressions.
- [ ] **Step 4: Link the two TSVs** from the doc, state the `chebi_primary_id` remap basis (plan "Spec gaps" 2), and change the doc's Status line to "implemented; measured on the <date> build".
- [ ] **Step 5: Spot checks** with the `cypher-queries` skill: nitrate is one node (`'chebi:14654' IN m.all_identifiers` on `kegg.compound:C00244`); "Alanine" measurements sit on L-Ala (`C00041`, `dominant_form_assumed`); the siderophore route returns families; L-Ala → Ala stereo hop and the L/D-Ala enantiomer edge exist.
- [ ] **Step 6: Commit** `git add docs/kg-changes/metabolite-identity-and-classes.md docs/kg-changes/metabolite-identity docs/superpowers/plans/2026-10-14-metabolite-identity-and-classes.md && git commit -m "docs(kg-changes): metabolite identity + classes — measured post-rebuild numbers and release TSVs"`.
- [ ] **Step 7: Release** stays held (spec §12): the bundled release follows the explorer's E1 resolver; `/release-kg` sets `mcp_min_version` and `breaking_changes` then. Nothing to do here beyond noting it in the doc.
