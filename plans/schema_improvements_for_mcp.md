# Schema Improvements for MCP/LLM Queryability

**Goal:** Make the KG queryable by an LLM agent via MCP, focused on analyzing
Prochlorococcus–Alteromonas positive interactions (ROS detoxification, nutrient
exchange, exoenzymes) from Weissberg 2025 and cross-paper evidence.

**Date:** 2026-03-11

---

## Decision Log

| # | Decision | Rationale |
|---|----------|-----------|
| D1 | Split `Affects_expression_of` into `Condition_changes_expression_of` + `Coculture_changes_expression_of` | Clearer semantics: "condition caused this" vs "coculture caused this" |
| D2 | Split `Affects_expression_of_homolog` into `Condition_changes_expression_of_ortholog` + `Coculture_changes_expression_of_ortholog` | Consistency; "ortholog" is more precise than "homolog" for Cyanorak clusters |
| D3 | Restrict coculture ortholog propagation to same taxonomic group | Picocyanobacteria↔picocyanobacteria, Alteromonas↔Alteromonas only. Filter: `h.distance <> 'cross phylum'` |
| D4 | `significant` bug tracked separately (see `memory/known_bugs.md`) | Different root cause, independent fix |
| D5 | Standardize all paperconfig organism names, condition_type values, and test_type values | Consistent vocabulary enables reliable LLM querying and cross-paper comparison |
| D6 | `pvalue_threshold` missing from Weissberg 2025 tracked separately (see `memory/known_bugs.md`) | Independent fix, unrelated to env conditions |
| D7 | `environmental_control_condition_id` is optional metadata — validator emits a **warning** (not error) when `environmental_treatment_condition_id` is set but `environmental_control_condition_id` is missing | Some controls have no defined env condition node (e.g., Bagby `gas_shock_air` control is "time 0 pre-shock"). Adapter ignores this field entirely; it does not affect edge routing or Publication edge emission. |

## Resolved Questions

1. **WH8102 genus** — keep as `Synechococcus WH8102` in paperconfigs and `preferred_name`. Papers universally use "Synechococcus"; the OrganismTaxon node already has `genus: Parasynechococcus` from NCBI for taxonomic accuracy. LLMs find it either way.
2. **Publication connectivity** — adapter-emitted edges (in `omics_adapter.py`). More explicit than post-import Cypher, and the adapter already knows which publications use which sources.
3. **Aharonovich/Biller 2016 env conditions** — yes, retrofit with `environmental_conditions` blocks for the common experimental setup (growth medium, light, temperature). These describe the baseline conditions, not the treatment (the treatment is the coculture organism). This gives these papers EnvironmentalCondition nodes even though the expression edges use OrganismTaxon as source.

---

## Reference: Design Details

### Edge Type Split

Current → Proposed:
```
(EnvironmentalCondition)-[:Affects_expression_of]->(Gene)        → Condition_changes_expression_of
(OrganismTaxon)-[:Affects_expression_of]->(Gene)                 → Coculture_changes_expression_of
(EnvironmentalCondition)-[:Affects_expression_of_homolog]->(Gene) → Condition_changes_expression_of_ortholog
(OrganismTaxon)-[:Affects_expression_of_homolog]->(Gene)          → Coculture_changes_expression_of_ortholog
```

Naming: "changes" (concrete), "coculture" (clear), "ortholog" (precise). `Gene_is_homolog_of_gene` keeps its name.

### Ortholog propagation restriction

| Group | Organisms | Allowed `distance` values |
|---|---|---|
| Picocyanobacteria (Synechococcales) | MED4, AS9601, MIT9301, MIT9312, MIT9313, NATL1A, NATL2A, RSP50, CC9311, WH8102 | `same species`, `same clade`, `same order` |
| Alteromonadaceae | MIT1002, EZ55, HOT1A3 | `same family` |

Filter: `WHERE h.distance <> 'cross phylum'`. Condition orthologs: no restriction.

### New edge properties (all 4 types)

- `omics_type: str` — RNASEQ | PROTEOMICS | METABOLOMICS | MICROARRAY
- `organism_strain: str` — which organism's genes are measured (e.g., "MED4")
- `treatment_condition: str` — what was applied (e.g., "Coculture with Alteromonas HOT1A3")
- `statistical_test: str` — DESeq2, edgeR, etc.

### `condition_category` derivation

`condition_category` on `EnvironmentalCondition` nodes is derived 1:1 from the canonical `condition_type` value (i.e., `condition_category = condition_type`). The controlled vocabulary table in "Paperconfig standardization" is the authoritative mapping — every canonical `condition_type` in that table is also a valid `condition_category` value.

### `Published_expression_data_about` edge

A `Publication` node emits one edge per distinct source node referenced in its analyses:
- `(Publication)-[:Published_expression_data_about]->(EnvironmentalCondition)` — for analyses using `environmental_treatment_condition_id`
- `(Publication)-[:Published_expression_data_about]->(OrganismTaxon)` — for analyses using `treatment_organism`

A single publication may have both types if it contains both condition-based and coculture analyses.

Note: `environmental_control_condition_id` is **not** used for Publication edges — it is paperconfig metadata only (adapter ignores it). Publication edges target treatment nodes only (what the experiment studied), not baseline conditions.

### Paperconfig standardization tables

**Organism names:**

| Current variants | Canonical form | Files to fix |
|---|---|---|
| `Alteromonas HOT1A3` | `Alteromonas macleodii HOT1A3` | Aharonovich 2016 |
| `Alteromonas EZ55` | `Alteromonas macleodii EZ55` | Hennon 2017 |
| `Prochlorococcus marinus MIT9313` (treatment_organism) | `Prochlorococcus MIT9313` | Aharonovich 2016 |

**`condition_type` controlled vocabulary:**

| Current value(s) | Canonical `condition_type` | Papers |
|---|---|---|
| `growth_medium` | `growth_medium` | Aharonovich, Biller 2016, Coe, Fang, He, Martiny, Read, Thompson 2016, Tolonen, Weissberg |
| `nutrient_stress` (nitrogen) | `nitrogen_stress` | Tolonen 2006, Read 2017 |
| `nutrient_stress` (phosphorus) | `phosphorus_stress` | Lin 2015, Martiny 2006 |
| `iron_stress` | `iron_stress` | Thompson 2011 |
| `salt_stress`, `salt_acclimation` | `salt_stress` | Al-Hosani 2015, He 2022 |
| `gas_shock`, `pco2` | `carbon_stress` | Bagby 2015, Hennon 2017, Barreto 2022 |
| `light_stress` | `light_stress` | Steglich 2006, Thompson 2016 |
| `light_stress` (extended darkness), `dark_tolerance` | `darkness` | Biller 2018, Coe 2024 |
| `plastic_leachate_stress` | `plastic_stress` | Tetu 2019 |
| `viral_lysis_products` | `viral` | Fang 2019 |
| `coculture` | `coculture` | Biller 2016 |
| `growth_state` | `growth_state` | Anjur 2025, Weissberg 2025 |
| `temperature_stress` | `temperature_stress` | Labban 2022 |

**`test_type` canonical forms:**

| Current value(s) | Canonical |
|---|---|
| `DESeq2` | `DESeq2` |
| `DESeq` | `DESeq` |
| `edgeR` | `edgeR` |
| `Rockhopper`, `Rockhopper/DESeq` | `Rockhopper` |
| `DESeq2/NOISeq` | `DESeq2` |
| `Affymetrix microarray` | `microarray` |
| `Affymetrix microarray with Cyber-T Bayesian t-test`, `...with Bayesian t-test (Cyber-T)` | `microarray_Cyber-T` |
| `Affymetrix microarray with LPE test` | `microarray_LPE` |
| `Affymetrix microarray with Goldenspike` | `microarray_Goldenspike` |

---

## Agent Roles

> Executable specs (frontmatter + system prompts) are in [`.claude/agents/`](../.claude/agents/). This table is the authoritative source; the spec files are derived from it.

| Agent | Role | Files owned | Skills used |
|---|---|---|---|
| **A: Paperconfig** | Standardize + add missing fields in all 24 YAML files | All `paperconfig.yaml` files | `/paperconfig` for validation |
| **B: Schema + Adapters** | Core code: schema, omics_adapter, cyanorak_ncbi_adapter, post-import.sh | `schema_config.yaml`, `omics_adapter.py`, `cyanorak_ncbi_adapter.py`, `post-import.sh` | — |
| **C: Tests** | Write + run unit tests AND KG validity tests | `tests/test_paperconfig_validation.py`, `tests/test_omics_adapter_organism_gene.py`, `tests/test_cyanorak_ncbi_adapter.py`, `tests/test_create_knowledge_graph.py`, `tests/kg_validity/` | `/omics-edge-snapshot` for before/after |
| **D: Skills** | Update/create skill definitions | `.claude/skills/` (`validate_paperconfig.py`, `omics-edge-snapshot`, `cypher-queries`, `paperconfig`) | — |
| **E: Docs** | Update documentation | `CLAUDE.md`, `plans/`, `.claude/memory/MEMORY.md` | — |
| **F: Validation runner** | Run skills at phase gates | Read-only | `/omics-edge-snapshot`, `/check-gene-ids`, `/cypher-queries` |
| **G: Code Review** | Review all changes before each phase gate | Read-only | — |
| **H: Plan Manager** | Maintain master plan, create detailed phase sub-plans, track decisions and blockers | `plans/schema_improvements_for_mcp.md`, `plans/phase_*.md` | — |

### Agent H: Plan Manager — Responsibilities

**Before each phase:**
- Create `plans/phase_N_<name>.md` with detailed sub-plan:
  - Exact files to edit with line-number references
  - Specific strings to search/replace (for paperconfig standardization)
  - Expected before/after for each change
  - Concrete test commands and expected output
  - Acceptance criteria checklist for G (code review)
- Assign tasks to specific agents with unambiguous scope

**During each phase:**
- Track blockers and decisions that arise
- Update master plan if scope changes (e.g., new inconsistency discovered in paperconfigs)
- Coordinate agent dependencies (e.g., "B can start T3b after G approves T3a")

**After each phase:**
- Record gate results in master plan (pass/fail, edge counts, any surprises)
- Update phase sub-plan with lessons learned
- Adjust downstream phases if needed (e.g., Phase 1 discovered new paperconfig issues → add to Phase 2)

### Agent G: Code Review — Checklist Per Phase

**Every phase:**
- No references to old edge labels (`Affects_expression_of`, `Affects_expression_of_homolog`) in changed files — except `Gene_is_homolog_of_gene` which is kept
- No old paperconfig vocabulary in changed files (organism names, condition_type, test_type)
- String literals match between schema `label_in_input` ↔ adapter edge emission ↔ test assertions ↔ skill templates

**Phase 1 (paperconfig standardization):**
- All 24 paperconfigs use canonical organism names (no `Alteromonas HOT1A3`, no `Prochlorococcus marinus`)
- All `condition_type` values are in the controlled enum
- All `test_type` values are in the canonical set
- All statistical_analyses have `treatment_condition` and `type`
- `validate_paperconfig.py` enforces ALL the above

**Phase 2 (missing env conditions):**
- New `environmental_conditions` blocks have `condition_type` from controlled enum
- `environmental_control_condition_id` / `environmental_treatment_condition_id` reference valid env condition keys
- Existing analyses in retrofitted papers still have correct source routing (treatment_organism → coculture, env_id → condition)

**Phase 3 (preferred_name):**
- `preferred_name` format: `"Prochlorococcus MED4"`, `"Alteromonas macleodii HOT1A3"`, `"Synechococcus WH8102"` etc. — consistent with paperconfig canonical names
- Treatment organisms (Phage, Marinobacter, etc.) also have `preferred_name`

**Phase 4 (edge split + new properties):**
- `schema_config.yaml`: 4 edge types defined with matching `label_as_edge` and `label_in_input`; all new properties listed on all 4 types; `Published_expression_data_about` edge defined (targets: both EnvironmentalCondition and OrganismTaxon); `condition_category` on EnvironmentalCondition
- `omics_adapter.py`: emits correct `label_in_input` strings matching schema; all new properties populated (no `None` for required fields); `condition_category` = `condition_type` value (must cover every value in controlled vocabulary table); Publication edges emitted for both source types
- No leftover `affects_expression_of` string literals in adapter
- Tests in `test_omics_adapter_organism_gene.py` assert both edge labels and all new properties

**Phase 5 (post-import homolog):**
- `post-import.sh`: two separate heredoc queries for condition and coculture propagation
- Condition query: no distance filter (propagates to all homologs)
- Coculture query: `WHERE h.distance <> 'cross phylum'`
- Both queries copy ALL new properties (`omics_type`, `organism_strain`, `treatment_condition`, `statistical_test`) plus existing ones
- Edge labels in Cypher match schema exactly (PascalCase as BioCypher outputs)

**Phase 6 (final):**
- `grep -r "Affects_expression_of" --include="*.py" --include="*.yaml" --include="*.sh" --include="*.cypher" --include="*.md"` returns ZERO hits (except historical plan files and this document)
- All test files reference new edge labels
- `snapshot_data.json` regenerated with new edge types
- CLAUDE.md "Actual Neo4j labels" section updated
- Skill SKILL.md files reference new edge types

---

## Tests Inventory

### Existing tests that need updating

| Test file | Current coverage | Changes needed |
|---|---|---|
| `test_paperconfig_validation.py` | Calls `validate_paperconfig.validate()` on all 24 configs | Passes after D updates validation skill |
| `test_omics_adapter_organism_gene.py` | Edge creation, properties, env condition nodes | Update for new edge labels (`condition_changes_expression_of`, `coculture_changes_expression_of`), new properties |
| `test_cyanorak_ncbi_adapter.py` | Gene/organism node creation | Add test for `preferred_name` on OrganismTaxon |
| `test_create_knowledge_graph.py` | Full integration (`@slow`) | Update expected edge labels |
| `tests/kg_validity/test_expression.py` | Expression edge properties, direction/sign consistency | Update for both new edge labels |
| `tests/kg_validity/test_post_import.py` | Homolog edges exist, bidirectional, properties | Update for both ortholog edge labels, add cross-phylum filter assertion |
| `tests/kg_validity/test_snapshot.py` | Regression snapshot | Regenerate `snapshot_data.json` after final rebuild |

### New tests to write

| Test | File | What it validates |
|---|---|---|
| **Paperconfig vocabulary** | `test_paperconfig_validation.py` (extend) | All `organism` values use canonical names; all `condition_type` in controlled enum; all `test_type` in canonical set; all analyses have `treatment_condition` and `type` |
| **Edge label routing** | `test_omics_adapter_organism_gene.py` (extend) | Config with `environmental_treatment_condition_id` → `condition_changes_expression_of`; config with `treatment_organism` → `coculture_changes_expression_of` |
| **New edge properties** | `test_omics_adapter_organism_gene.py` (extend) | `omics_type`, `organism_strain`, `treatment_condition`, `statistical_test` populated on edges |
| **Publication edges** | `test_omics_adapter_organism_gene.py` (extend) | `published_expression_data_about` edges emitted from Publication → source nodes |
| **condition_category derivation** | `test_omics_adapter_organism_gene.py` (extend) | `condition_category` on EnvironmentalCondition nodes derived from `condition_type` |
| **Preferred name** | `test_cyanorak_ncbi_adapter.py` (extend) | OrganismTaxon nodes have `preferred_name` = genus + strain |
| **No cross-phylum coculture orthologs** | `tests/kg_validity/test_post_import.py` (extend) | `Coculture_changes_expression_of_ortholog` edges all have `distance <> 'cross phylum'` |
| **Ortholog new properties** | `tests/kg_validity/test_post_import.py` (extend) | `omics_type`, `organism_strain`, `treatment_condition`, `statistical_test` propagated to ortholog edges |

### Skills that need updating (owned by Agent D)

| Skill | What changes | Impact |
|---|---|---|
| `validate_paperconfig.py` | Add canonical vocabulary checks: organism names, condition_type enum, test_type set, required fields | `test_paperconfig_validation.py` inherits these checks automatically |
| `omics-edge-snapshot/SKILL.md` | Count `Condition_changes_expression_of` + `Coculture_changes_expression_of` instead of `Affects_expression_of`; report separately | `/deploy-strain` inherits the fix |
| `cypher-queries/SKILL.md` | Replace `Affects_expression_of` templates with new edge type templates; add condition vs coculture query examples | — |
| `paperconfig/SKILL.md` | Document which fields map to which edge labels; clarify `treatment_organism` → coculture edge, `environmental_treatment_condition_id` → condition edge | — |

---

## Phased Execution Plan

Each phase: **scope → implement → test → gate** before proceeding.

---

### Phase 1: Paperconfig Standardization

**Scope:** Normalize all 24 paperconfig files to use consistent vocabulary. Update validation skill to enforce canonical vocabulary. No adapter/schema code changes.

**Agents:** H (create `plans/phase_1_paperconfig_standardization.md`), A (paperconfig edits), D (update `validate_paperconfig.py`), C (run tests), G (review), F (validation)

**Implement (sequential ordering matters):**
1. **D (first):** Update `validate_paperconfig.py` to check canonical vocabulary (organism names, condition_type enum, test_type set, required fields) — A needs this to verify their edits
2. **A (after D):** Standardize organism names, condition_type, test_type, fill missing treatment_condition across all 24 paperconfigs (see standardization tables above), using updated validator to confirm each file
3. **C (after D, can run in parallel with A):** Extend `test_paperconfig_validation.py` with vocabulary-specific assertions if needed beyond what the skill validates

**Review (G):** Verify all 24 paperconfigs use canonical vocabulary, `validate_paperconfig.py` enforces all constraints

**Test:**
- `pytest tests/test_paperconfig_validation.py -v` — all 24 configs pass with new vocabulary checks
- `uv run python create_knowledge_graph.py --test` — no build regressions
- `pytest -m "not slow and not kg"` — all unit tests pass
- Grep paperconfigs for old values: no `Alteromonas HOT1A3`, no `nutrient_stress`, no `Affymetrix microarray with`

**Gate:** G approves, all tests pass, no old vocabulary remains. Commit.

---

### Phase 2: Add Missing Environmental Conditions

**Scope:** Retrofit Aharonovich 2016, Biller 2016, Lindell 2007 with `environmental_conditions` blocks for experimental baselines. (`pvalue_threshold` for Weissberg 2025 is tracked separately in `memory/known_bugs.md` — not in scope here.)

**Agents:** H (create `plans/phase_2_missing_env_conditions.md`), A (paperconfig edits), C (run tests), G (review), F (validation)

**Implement:**
1. **A:** Add `environmental_conditions` blocks to Aharonovich, Biller 2016, Lindell 2007; add env IDs to their analyses (no code changes needed — paperconfig only)

**Review (G):** Verify env condition IDs are valid references, condition_type uses canonical values, existing analysis source routing unchanged

**Test:**
- `pytest tests/test_paperconfig_validation.py -v` — passes
- `uv run python create_knowledge_graph.py --test` — builds successfully
- Verify new EnvironmentalCondition nodes appear (count > 56)
- Verify expression edge counts unchanged (188K — edge routing not changed yet)
- `pytest -m "not slow and not kg"` — passes

**Gate:** G approves, new env condition nodes created, existing edges intact. Commit.

**Gate result (2026-03-11): PASSED**
- Lindell 2007: added `environmental_conditions.lindell_2007_pro99` (`condition_type: growth_medium`); `environmental_control_condition_id` added to all 8 analyses. Node confirmed created in test build.
- Aharonovich 2016 and Biller 2016: already had `environmental_conditions` blocks; `environmental_control_condition_id` added to all analyses.
- Bagby 2015 (expanded scope): 3 of 4 analyses had matching control env condition (`bagby_2015_air`) — added. 1 (`gas_shock_air`) legitimately has no control env node (time 0 baseline).
- Validator updated: warns (not errors) when `environmental_treatment_condition_id` set but `environmental_control_condition_id` missing. 35/35 tests pass.

---

### Phase 3: OrganismTaxon Preferred Names

**Scope:** Small, isolated change to `cyanorak_ncbi_adapter.py`. Done before the edge split to get an easy win committed independently.

**Agents:** H (create `plans/phase_3_organism_names.md`), B (adapter code), C (tests), G (review)

**Implement (sequential):**
1. **B:** Update `cyanorak_ncbi_adapter.py` — set `preferred_name` from organism config
2. **C:** Add test in `test_cyanorak_ncbi_adapter.py` for `preferred_name` — after B is done

**Review (G):** `preferred_name` format consistent with canonical organism names from paperconfigs; treatment organisms also have `preferred_name`

**Test:**
- `pytest tests/test_cyanorak_ncbi_adapter.py -v` — passes with preferred_name assertion
- `uv run python create_knowledge_graph.py --test` — verify in output CSVs
- `pytest -m "not slow and not kg"` — passes

**Gate:** G approves, all OrganismTaxon nodes have readable names. Commit.

---

### Phase 4: Edge Type Split + New Properties

**Status:** DONE — Gate PASSED 2026-03-11. 188,501 → 188,501 edges (170,904 Condition + 17,597 Coculture). Zero regressions.

**Scope:** The core schema change. Split expression edges into 4 types. Add new properties. Add Publication connectivity edges. Add `condition_category` to EnvironmentalCondition nodes.

**Agents:** H (create `plans/phase_4_edge_split.md`), B (schema + adapter code), C (update + run tests), D (update skills), E (update docs), G (review)

**Pre-condition (F):** Run `/omics-edge-snapshot --save pre_phase4` *before* any code changes to establish baseline. Verify current `Affects_expression_of` total matches expected ~188K. This snapshot is the reference for data-loss detection after the split.

**Implement (sequential ordering matters):**
1. **B (schema first):** Update `schema_config.yaml` — 4 new edge type definitions, new properties, `condition_category` on EnvironmentalCondition, `Published_expression_data_about` edge type
2. **B (adapter second, after schema):** Update `omics_adapter.py` — edge label routing by source type, populate `omics_type`/`organism_strain`/`treatment_condition`/`statistical_test`, derive `condition_category` from `condition_type` (see controlled vocabulary table above — `condition_category = condition_type`), emit Publication→source edges (both `EnvironmentalCondition` and `OrganismTaxon` targets)
3. **C (after B):** Update `test_omics_adapter_organism_gene.py` — test edge label routing, new properties, Publication edges, condition_category. Run `pytest tests/test_omics_adapter_organism_gene.py -v`
4. **D (after B, can run in parallel with C):**
  - Update `omics-edge-snapshot` skill — count new edge types instead of `Affects_expression_of`
  - Update `cypher-queries` skill — new edge type templates
  - Update `paperconfig` skill docs — document which config maps to which edge label
5. **E (after B, can run in parallel with C and D):**
  - Update `CLAUDE.md` — new edge labels in "Actual Neo4j labels" section, update query examples

**Review (G):** Schema `label_in_input` strings match adapter emission; all 4 edge types have identical property lists in schema; `condition_category` derivation covers all canonical `condition_type` values (every value in the controlled vocabulary table); no leftover `affects_expression_of` literals; test assertions cover both edge labels + all new properties; Publication edge targets include both EnvironmentalCondition and OrganismTaxon

**Test:**
- `pytest tests/test_omics_adapter_organism_gene.py -v` — new edge label tests pass
- `uv run python create_knowledge_graph.py` — full build
- `/omics-edge-snapshot --compare pre_phase4` — verify condition + coculture total ≈ 188K (no data loss); expected split: ~171K condition + ~17K coculture
- Verify new properties populated: `grep omics_type biocypher-out/*/Condition*`
- Verify Publication edges: count in output CSVs
- Verify `condition_category` on EnvironmentalCondition nodes
- `pytest -m "not slow and not kg"` — passes

**Gate:** G approves, correct edge split, all new properties populated, Publication edges exist, no data loss vs pre_phase4 snapshot. Commit.

---

### Phase 5: Post-Import Homolog Propagation

**Scope:** Split homolog propagation in `post-import.sh`. Add cross-phylum filter for coculture. Copy new properties to ortholog edges.

**Agents:** H (create `plans/phase_5_post_import.md`), B (post-import.sh), C (tests), G (review), F (validation)

**Implement (sequential ordering matters):**
1. **B:** Update `post-import.sh` heredocs:
  - Split into two queries: `Condition_changes_expression_of` → `Condition_changes_expression_of_ortholog` (all homologs), `Coculture_changes_expression_of` → `Coculture_changes_expression_of_ortholog` (WHERE `h.distance <> 'cross phylum'`)
  - Copy new properties (`omics_type`, `organism_strain`, `treatment_condition`, `statistical_test`)
  - Update `post-import.cypher` documentation
2. **C (after B):** Update `tests/kg_validity/test_post_import.py`:
  - Test both ortholog edge types exist
  - Assert zero cross-phylum coculture ortholog edges
  - Assert new properties propagated

**Review (G):** Cypher edge labels match BioCypher PascalCase output exactly; condition query has no distance filter; coculture query has `h.distance <> 'cross phylum'`; all new + existing properties copied in both queries; no remaining `Affects_expression_of_homolog` references

**Test:**
> **Note:** Full Docker build takes ~1 hour. Plan to run overnight or in background. Check with `docker compose logs -f post-process` to monitor homolog propagation step.
- Full Docker build: `docker compose up -d`
- `/omics-edge-snapshot --save post_phase5` — compare with pre_phase4
- Verify: `Condition_changes_expression_of_ortholog` ≈ 778K
- Verify: `Coculture_changes_expression_of_ortholog` < 83K (cross-phylum removed)
- Cypher: `MATCH ()-[e:Coculture_changes_expression_of_ortholog]->() WHERE e.distance = 'cross phylum' RETURN count(e)` → 0
- `pytest tests/kg_validity/test_post_import.py -v` — passes

**Gate:** G approves, ortholog split correct, cross-phylum filter working, new properties propagated. Commit.

---

### Phase 6: Final Validation + Docs

**Status:** DONE — Gate PASSED 2026-03-11. 376 KG validity tests pass (2 known orphan-protein failures, pre-existing). Global sweep clean. Snapshot regenerated.


**Scope:** Update remaining KG validity tests, regenerate snapshot, run full validation suite, finalize documentation.

**Agents:** H (create `plans/phase_6_final_validation.md`), C (all tests), D (final skill polish), E (docs), G (final review), F (run verification queries)

**Implement:**
- **C:**
  - Update `tests/kg_validity/test_expression.py` — both new edge labels, new properties
  - Regenerate `snapshot_data.json`: `uv run python tests/kg_validity/generate_snapshot.py`
  - Update `test_create_knowledge_graph.py` if needed
- **D:** Final pass on all skill docs
- **E:** Update `CLAUDE.md` sections: "Actual Neo4j labels", "Key graph facts", expression edge documentation. Update `MEMORY.md`.

**Review (G):** Full codebase sweep: `grep -r "Affects_expression_of" --include="*.py" --include="*.yaml" --include="*.sh" --include="*.cypher"` returns zero hits (except plan files); all test files reference new edge labels; CLAUDE.md consistent with actual schema; skill templates match deployed edge types

**Test:**
- `pytest tests/kg_validity/ -v` — all KG validity tests pass
- `pytest -m "not slow and not kg"` — all unit tests pass
- **F:** Run verification queries (see below) via `/cypher-queries`
- **F:** `/check-gene-ids all` — no regressions in gene ID matching

**Gate:** G approves full sweep, all tests green, all verification queries return expected results. Commit + tag release.

---

## Verification Queries (Phase 6, run by Agent F via `/cypher-queries`)

```cypher
-- 1. ROS detoxification genes upregulated in Alteromonas during coculture with Pro
MATCH (src:OrganismTaxon)-[e:Coculture_changes_expression_of]->(g:Gene)
      -[:Gene_belongs_to_organism]->(org:OrganismTaxon)
WHERE e.omics_type = 'RNASEQ'
  AND e.expression_direction = 'up'
  AND e.organism_strain = 'HOT1A3'
  AND (g.product CONTAINS 'catalase' OR g.product CONTAINS 'superoxide'
       OR g.product CONTAINS 'peroxidase' OR g.product CONTAINS 'thioredoxin')
RETURN src.preferred_name AS partner, g.gene_name, g.product,
       e.log2_fold_change, e.time_point
ORDER BY e.log2_fold_change DESC

-- 2. Cross-paper nitrogen stress in Prochlorococcus
MATCH (env:EnvironmentalCondition)-[e:Condition_changes_expression_of]->(g:Gene)
WHERE env.condition_category = 'nitrogen_stress'
  AND e.organism_strain = 'MED4'
  AND e.significant = true
RETURN env.name, e.publications, g.gene_name, g.product, e.log2_fold_change
ORDER BY abs(e.log2_fold_change) DESC LIMIT 20

-- 3. Nutrient exchange — Alteromonas transporters responding to coculture
MATCH (src:OrganismTaxon)-[e:Coculture_changes_expression_of]->(g:Gene)
WHERE e.organism_strain = 'HOT1A3'
  AND e.expression_direction = 'up'
  AND (g.product CONTAINS 'transporter' OR g.product CONTAINS 'permease'
       OR g.product CONTAINS 'porin' OR g.product CONTAINS 'secretion')
RETURN src.preferred_name AS partner, g.gene_name, g.product,
       e.log2_fold_change, e.time_point, e.omics_type
ORDER BY e.log2_fold_change DESC

-- 4. Ortholog inference within Alteromonas
MATCH (src:OrganismTaxon)-[e:Coculture_changes_expression_of_ortholog]->(g:Gene)
      -[:Gene_belongs_to_organism]->(org:OrganismTaxon)
WHERE e.organism_strain IN ['EZ55', 'MIT1002']
  AND e.expression_direction = 'up'
  AND g.product CONTAINS 'catalase'
RETURN src.preferred_name, org.strain_name, g.gene_name, g.product,
       e.log2_fold_change, e.distance, e.original_gene

-- 5. Publication connectivity
MATCH (pub:Publication)-[:Published_expression_data_about]->(src)
RETURN pub.title, labels(src)[0] AS source_type, src.name

-- 6. Verify no cross-phylum coculture ortholog edges
MATCH ()-[e:Coculture_changes_expression_of_ortholog]->()
WHERE e.distance = 'cross phylum'
RETURN count(e)  // should be 0

-- 7. Verify standardized organism_strain values on coculture edges
MATCH ()-[e:Coculture_changes_expression_of]->()
RETURN DISTINCT e.organism_strain ORDER BY e.organism_strain

-- 8. Verify standardized condition_category values
MATCH (env:EnvironmentalCondition)
RETURN DISTINCT env.condition_category ORDER BY env.condition_category

-- 9. Verify new properties on condition edges
MATCH (env:EnvironmentalCondition)-[e:Condition_changes_expression_of]->(g:Gene)
WHERE e.omics_type IS NOT NULL
RETURN e.omics_type, e.statistical_test, e.organism_strain, count(*) AS cnt
ORDER BY cnt DESC LIMIT 20

-- 10. Edge count summary
MATCH ()-[e:Condition_changes_expression_of]->() RETURN 'condition_direct' AS type, count(e) AS cnt
UNION ALL
MATCH ()-[e:Coculture_changes_expression_of]->() RETURN 'coculture_direct' AS type, count(e) AS cnt
UNION ALL
MATCH ()-[e:Condition_changes_expression_of_ortholog]->() RETURN 'condition_ortholog' AS type, count(e) AS cnt
UNION ALL
MATCH ()-[e:Coculture_changes_expression_of_ortholog]->() RETURN 'coculture_ortholog' AS type, count(e) AS cnt
```
