# Two-state strings (R5 completion) — 2026-08-28

**Scope.** Convert the 7 grandfathered `"true"`/`"false"` properties + the undeclared
`Assay_flags_metabolite.flag_value` to named pairs; keep property names; keep the
paperconfig author contract (`"true"`/`"false"` tokens) and map in the adapters; drop
`bool_string` from `VALUE_TYPES`. Rode along: Meiothermus `name_synonyms`/`taxonomy_note`
+ `organismTaxonFullText` (backlog #3, `preferred_name` unchanged by decision), the
treatment taxid 1299 → 277 bug found on the way, and the TCDB threshold re-measurement
(kept at 50, documented).

**Not changing.** MCP output shape (explorer coerces to bool at its boundary), any
paperconfig, `Schema_info.git_dirty`, post-import row selection.

**Acceptance.** `pytest -m "not slow and not kg"` green; after rebuild `pytest -m kg`
green incl. the enum tests and `test_controlled_vocabularies`; verification queries in
`docs/kg-changes/two-state-strings.md` return only the pairs; snapshot regenerated.

**Status.** DONE 2026-08-28 — rebuilt, `pytest -m kg` 1197 passed (after fixing the one
missed post-import site: DerivedMetric `flag_true_count` compared `r.value = 'true'`),
snapshot regenerated, hand-off verification queries all clean (0 residual literals,
9 relationship indexes ONLINE, `ncbitaxon:277` targeted). Explorer pick-up pending.
