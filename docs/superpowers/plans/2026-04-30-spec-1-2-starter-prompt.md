# Spec 1.2 starter prompt — for the next chat

Copy the block below into the new chat to resume.

---

I'm continuing the metabolite-scaffold work. Phase 1.1B (resolver + accessors + step-2 transforms) is done on branch `metabolite-scaffold` — 15 commits, not merged per my instruction. We had a long design discussion that pivoted Spec 1.2 in three load-bearing ways. Read these before doing anything:

1. **Revised spec:** `docs/superpowers/specs/2026-04-30-metabolite-reactions-scaffold-revised.md` (supersedes the 2026-04-28 spec, which now carries a SUPERSEDED banner)
2. **Design decisions memory:** `memory/spec_1_2_design_decisions.md` — the *why* behind every load-bearing choice
3. **Phase 1.1B completed work:** branch `metabolite-scaffold`, commits `8743e0e..0102115`, 15 commits

**Headline decisions to anchor on (don't relitigate without strong reason):**

- **KG is fully KEGG-based.** Primary node IDs are KEGG (`kegg.reaction:R*`, `kegg.compound:C*`). MNX is the dedup hub for (a) cross-reference enrichment and (b) Phase 2 paper integration only. Mixed-prefix fallback (`chebi:` then `mnx:`) only for compounds without a KEGG entry.
- **Reaction direction is not modelled.** Single direction-agnostic `Reaction_has_metabolite` edge. NO substrate/product split — we verified empirically that all 83,796 MNX reactions use ` = ` and KEGG's `<=>` is editorial convention.
- **No build-time resolution.** The 2.6 GB MNX SQLite resolver opens at *prepare-data* time only. Build-time adapter reads pre-computed JSON files. Pattern: `build_og_descriptions.py` (step 5) — prune to gene-reachable subset, enrich, write small JSON cache.

**Schema (Spec 1.2):**
- New nodes: `Reaction`, `Metabolite`
- New edges: `Gene_catalyzes_reaction`, `Reaction_has_metabolite`, `Reaction_in_kegg_pathway`, `Organism_has_metabolite`
- No `Reaction_in_organism` edge — denormalized `Reaction.organisms[]` property instead.

**Pipeline:**
- New step 6 in `prepare_data.sh`: `build_kegg_metabolism_xrefs.py` — prunes KEGG to gene-reachable R-numbers/C-numbers, enriches with KEGG metadata + MNX xrefs, writes `cache/data/kegg/kegg_metabolism_xrefs.json` (~1 MB).
- Extend `multiomics_kg/utils/kegg_utils.py` with 5 new REST endpoints: `/list/reaction`, `/list/compound`, `/link/compound/reaction`, `/link/pathway/reaction`, `/link/pathway/compound`.
- Drop `transform: resolve_kegg_reaction_to_mnxr` line from `config/gene_annotations_config.yaml` so `gene["kegg_reactions"]` stays as raw KEGG R-numbers.
- Add 2 indexes to MNX resolver DDL: `idx_compound_aliases_mnxm`, `idx_reaction_aliases_mnxr`.
- Add helpers `mnxm_to_primary_id()` and `mnxr_to_primary_id()` to `metabolite_utils.py` (KEGG > ChEBI/Rhea > MNX fallback).

**Task right now:** Write the Spec 1.2 implementation plan (mirror the format of `docs/superpowers/plans/2026-04-28-metabolite-foundation-phase-b-resolver.md` — 13 TDD-discipline tasks, each with files-changed list, test-first steps, and a commit message). The plan should split into:

1. Extend `kegg_utils.py` with 5 new endpoints + tests
2. Drop the YAML transform + verify gene_annotations_merged.json shape unchanged
3. Add resolver DDL indexes + helpers in `metabolite_utils.py` + tests
4. New `build_kegg_metabolism_xrefs.py` (the prune-then-enrich step) + tests + wire into `prepare_data.sh`
5. `schema_config.yaml` additions (2 nodes, 4 edges)
6. New `metabolism_adapter.py` + unit tests
7. Wire into `create_knowledge_graph.py`
8. Post-import Cypher (rollups + indexes + Organism_has_metabolite materialization)
9. KG validity tests
10. Integration smoke test against real cache (slow marker)
11. Validation gate against MED4 (real data — confirm ≥2K Reaction nodes, ≥1K Metabolite nodes, sensible counts)
12. Update `CLAUDE.md` with new node types, edge types, and step 6 docs

Use the `superpowers:writing-plans` skill if available. After the plan is written, we'll execute via `superpowers:subagent-driven-development` (same workflow as Phase 1.1B).

Don't start implementing yet. Just read the spec + memory + 1.1B plan as a template, then propose the task breakdown.

---

# End of starter prompt
