# InterPro multi-ontology redesign — open items / backlog

Follow-ups from the 2026-08-17 redesign (spec
`docs/superpowers/specs/2026-08-17-interpro-multi-ontology-redesign-design.md`,
shipped per `docs/kg-changes/interpro-multi-ontology.md`, merged `17fdfd51`).
All items below were explicitly triaged **OK-to-defer** by the final
whole-branch review — none affects correctness of the shipped graph.

## Small code/test cleanups

All 7 items landed 2026-08-18 (each claim re-verified before changing; the
DUF regex tightening was confirmed a no-op on the live graph — same 9
NcbifamFamily nodes match either pattern). See CHANGELOG `[Unreleased]` →
`### Fixed`. The `acc_to_ipr` last-write-wins watch note now lives as a code
comment in `ncbifam_adapter.get_edges` — re-check it after any InterProScan
version-bump re-scan.

## Process / next rebuild

- [x] **Capture `annotation_state` distribution pre-rebuild** — landed
  2026-08-18: `tests/kg_validity/capture_annotation_state.py`
  (`--save` / `--compare`, omics-edge-snapshot pattern) + the committed
  `annotation_state_baseline.json` captured from the 2026-08-18 build
  (124,751 genes: no_evidence 12,061 / catch_all_only 5,988 /
  informative_single 12,642 / informative_multi 94,060). Run `--save`
  before each rebuild, `--compare` after.
- [x] **Orphan-protein Known Issue** — verified 2026-08-18 on a fresh rebuild:
  **still real** (38% orphaned, 25,441/67,024); the tests pass only because
  their thresholds were loosened to <50%. Not closeable — tracked as a live
  investigation in `plans/backlog.md` → `plans/orphan_proteins.md`.

## Research-driven follow-ups (spec §8, interaction-mechanism framing)

- [ ] **antiSMASH BGC detection** — `/add-a-tool` candidate; siderophore /
  vitamin cross-feeding currency; per-genome region calls (new calls.json
  shape per the add-a-tool table).
- [ ] **`Ncbifam_family_has_tigr_role` bridge** — TIGR\* families → existing
  `TigrRole` nodes via JCVI's frozen `TIGRFAMS_ROLE_LINK` archive. Gives
  non-Cyanorak strains (all Alteromonas) role-category access via 2-hop
  `gene → NcbifamFamily → TigrRole`, plus a Cyanorak-role QC cross-check.
  Partial by nature (NF* families have no roles).
- [ ] **NCBIfam EC/GO xref propagation to genes** — `ncbifam_reference.json`
  already carries them; add as another gated enrichment source in the
  step-2 merge (same pattern as InterPro's).
- [ ] **MCP / explorer surfacing** — NcbifamFamily ontology, InterPro/NCBIfam
  descriptions, source/evidence edge filters, router-mode traversals
  (carried over from the two-layer doc's §MCP follow-up; still not
  surfaced).
- [ ] **Further member-DB elevation** (PANTHER, HAMAP-as-ontology, CDD,
  folds…) — each is one merge-config line + small adapter, per the faceted
  design. Only on demonstrated need.
- [ ] **MetaCyc** — dormant in `interpro_reference.json` (`pathways` key);
  measured 2026-08-18 and found thin (dark-gene rescue 554/271, no pathway
  names without the MetaCyc license) — see the evidence-backed deferral in
  `plans/backlog.md`; elevate only with a concrete use case.

## Operational notes (not tasks)

- InterPro descriptions shipped **observed-only** (full set = 27.1 MB > 25 MB
  gate): adding a strain ⇒ re-run the reference step (prepare_data step 9)
  or new entries fail-soft to missing descriptions. Already documented in
  CLAUDE.md; noted here because it's the one corpus-coupling in the design.
- If several adapters ever need the same pruned InterPro tree, that is the
  moment to materialize a committed pruned artifact (spec §3.2 rationale).
