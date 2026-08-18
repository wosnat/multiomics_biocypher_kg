# InterPro multi-ontology redesign — open items / backlog

Follow-ups from the 2026-08-17 redesign (spec
`docs/superpowers/specs/2026-08-17-interpro-multi-ontology-redesign-design.md`,
shipped per `docs/kg-changes/interpro-multi-ontology.md`, merged `17fdfd51`).
All items below were explicitly triaged **OK-to-defer** by the final
whole-branch review — none affects correctness of the shipped graph.

## Small code/test cleanups

- [ ] **Multi-xref fan-out test** — `tests/test_interproscan.py` has no
  fixture with a result whose `xref` list carries 2+ accessions (the
  sequence-dedup fan-out path). A 3-line fixture change would also pin the
  shallow-copy behavior (`dict(record)` shares nested `libraries` /
  `interpro_entries` structures across fanned-out accessions — read-only
  today; a docstring note or deepcopy if a consumer ever mutates).
- [ ] **DUF regex alignment** — post-import Cypher flags
  `.*DUF\d.*` while `config/uninformative_terms.yaml` says `\bDUF\d`.
  Align (Cypher needs `\\b`, both `.cypher` and `.sh`) on the next
  post-import touch. Benign sentinel over-match meanwhile.
- [ ] **DataSource `interproscan.info_types`** — auto-derived from declared
  YAML fields only, so it reads `["interpro_entries", "ncbifam_ids"]` and
  misses the enrichment contributions (pfam/go/ec/cazy/naming recovery).
  Spec §4.2 wanted the full list. Fix = extend `_derive_info_types` or add
  an explicit override in the YAML. Pre-existing mechanism limit (old
  Layer B was equally invisible).
- [ ] **`normalize_strain` sentinel_rate denominator** — `--normalize` passes
  `input_proteins=len(calls)` (calls_made), while scan mode uses the FASTA
  count; `sentinel_rate` semantics differ slightly between modes. Cosmetic
  QC drift; fix = count FASTA headers or document the difference in the
  runner.
- [ ] **Stale mid-branch comments** — `tests/test_interproscan_consistency.py`
  "Mid-branch note (Task 10)" docstring and `create_knowledge_graph.py`
  "empty until Task 18 lands" comment are now historical; trim on next
  touch.
- [ ] **`enrich_interpro_fields` docstring** — claims the `go_term_donors`
  fallback matches the old behavior "exactly"; evidence-strength labeling
  can differ in mixed-type donor-order cases (new behavior is more
  correct). Reword.
- [ ] **`acc_to_ipr` last-write-wins** — if one NCBIfam accession were ever
  attributed two different IPR entries (cross-release drift), the bridge
  keeps the last. Mirrors the pre-existing `pf_to_ipr` convention; only
  matters after an InterProScan version-bump re-scan — check then.

## Process / next rebuild

- [ ] **Capture `annotation_state` distribution pre-rebuild** — alongside the
  `/omics-edge-snapshot` baseline, so bucket-movement claims are
  independently reproducible (the 2026-08-17 baseline was supplied context,
  not a captured artifact).
- [ ] **Orphan-protein Known Issue** — `test_no_orphan_proteins*` passed on
  the 2026-08-17 rebuild; verify on the next rebuild, then close the
  CLAUDE.md Known Issue + `plans/orphan_proteins.md`.

## Research-driven follow-ups (spec §8, interaction-mechanism framing)

- [ ] **MEROPS peptidase classification** — `/add-a-tool` candidate; the one
  real ontology gap for the exoproteolysis mechanism (secreted-protease
  shortlists = MEROPS family × SignalP × PSORTb × coculture expression).
  BRITE ko01002 + InterPro FAMILY entries approximate it meanwhile.
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
  elevate only if a MetaCyc ontology ever earns its keep vs KEGG.

## Operational notes (not tasks)

- InterPro descriptions shipped **observed-only** (full set = 27.1 MB > 25 MB
  gate): adding a strain ⇒ re-run the reference step (prepare_data step 9)
  or new entries fail-soft to missing descriptions. Already documented in
  CLAUDE.md; noted here because it's the one corpus-coupling in the design.
- If several adapters ever need the same pruned InterPro tree, that is the
  moment to materialize a committed pruned artifact (spec §3.2 rationale).
