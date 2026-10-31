# Transport programme — small follow-ups

Branch `worktree-transport-small-followups` (worktree off `main` @ `0b91c60f`).
Source: `plans/backlog.md` → Transport programme, the open sub-bullets under
Project 2 / 3a / the TCDB integration review. No Docker rebuild and no
prepare-data run without explicit approval.

## Status (2026-10-31): DONE — merged to main, step 6 + step 2 run, Docker rebuilt (`ef72f299`), `pytest -m kg` 1,431 passed, omics snapshot + annotation-state unchanged; explorer asks doc updated with live counts.

## Status (2026-10-29)

Code + unit tests + docs DONE in the worktree (unit suite 3,582 passed). Waiting
on approval for the three runs: step 2 `--force` (A), step 6 `--force` (B), Docker
rebuild (A, B, C), then `pytest -m kg` + `/omics-edge-snapshot`.

Decisions taken while implementing:
- **C ships as a new dense `display_name`, not a changed `name`.** `t.name = t.tcdb_id`
  is a deliberate marker (2026-08-12 names spec 5c chose it over a `name_source`
  property; `tests/kg_validity/test_tcdb_cazy.py` and the tcdb-diamond skill use it).
- **A dry run:** 16,039 lines over 117 strains (12,317 rank-1 hits below the gate,
  0 not kept), 8,643 with substrates; MED4 urtA →
  `[tcdb_best_hit] The high-affinity (<1 μM) urea porter [3.A.1.4.4] :: [tcdb_substrates] Urea` (format revised 2026-10-31 after the explorer response: bracketed underscore fixed words, cap 15 → 8).
- **B measured:** 52 nodes satisfy the rule; 4 paper rows flip (name matches only).
  Rows giving a KEGG id stay `as_reported` (the paper named that entry), and
  the 6 curated alias rows on the same nodes (biller SAH / GSSG / Vitamin B5,
  Kujawinski GSSG / N-acetylglutamate / N-acetylmuramate) keep their explicit
  `as_reported`. **Resolved (user, 2026-10-29):** those values were set while fixing
  IDs, not as stereo calls. Aliases now take the name-match rule when they omit
  `stereo_assignment`; an explicit value wins but is listed in the identity report's
  `alias_stereo_review` when it differs. 7 aliases corrected to `dominant_form_assumed`
  (11 paper names flip in all). biller SAM / glucose 1-phosphate / methionine sulfoxide
  keep their curated `dominant_form_assumed` (the rule has no ChEBI evidence for them),
  and they stay in the review list.
- **Possible next step (not done):** a third signal, "the node's own name is
  stereo-qualified (L-, D-, S-, R-) and the paper's is not", would catch those three
  and their name-matched equivalents. It needs its own measurement, because geometric
  prefixes would be false positives.

## Triage (measured 2026-10-29)

| # | Backlog item | Finding | Action |
|---|---|---|---|
| 1 | MEROPS `qcov` declared as a fraction | Already fixed in `9cf8516a` (yaml says 0–100) | Close bullet |
| 2 | 4 unnamed best-hit systems | 3.A.1.129.1 is re-anchored to 3.A.1.106.22 (`tcdb_moved.json`); 1.A.77.3.21 is named "AtpI of Synechococcus sp." from its protein cell (42 genes). Left: 1.B.1.12.1 (0 genes) and 2.A.49.2.14 (1 gene) — not on TCDB's pages, bare TCID by design | Close bullet with that note |
| 3 | `metabolite_identity_report.json` `remap` drifts on rerun | By design: `remap`/`repoints`/`dropped` are diffs against the *previous* `kegg_data.json` (code comment in `pipeline.py`); the release record is `docs/kg-changes/metabolite-identity/metabolite_id_remap.tsv` | Close bullet; make the report say so (`baseline` text) |
| 4 | ChEBI class substrates are generic Metabolites | Correct by design; the explorer has to handle it | Move to Explorer / MCP coordination as a hand-off |
| 5 | Reciprocal-check recall cost | A design question for 3b (family-level reciprocity), not a small fix | Leave in the backlog, re-parented under 3b |
| 6 | Best-hit names / substrates searchable from the gene | 16,039 genes have a rank-1 best hit at ≥ 40% identity; 16,020 of those systems have a scraped name | **Implement (A)** |
| 7 | Stereo forms with no ChEBI enantiomer read `as_reported` | A rule is possible that avoids the rejected `is_a` approach | **Implement (B)**, measure first |
| 8 | TCDB name quality | Systems already fall back to their protein cell; what's left is about 1,186 header-less subfamilies (551 with genes), still named by bare TCID | **Optional (C)**, your call |

## A. `[tcdb_best_hit]` description lines (step 2)

- In `build_gene_annotations.py`, next to `load_tcdb_diamond`: for the protein's
  **rank-1** best-hit system with **identity ≥ 40** (the `direct` gate's best-hit arm:
  rank 1, tier ≤ 2), append one line to `alternate_functional_descriptions`:
  `[tcdb_best_hit] <system name> (TC <id>) :: substrates: <name>, <name>`.
  The substrate clause is omitted when the system has no substrate.
- Names come from `tcdb_names.json` (falls back to the TCID). Substrates come from
  `tcdb_pruned.json` `system_substrates` → `kegg_data.json` compound /
  `additional_compounds` names. No `|` or `'` (handled by `clean_text`).
- Not included: near-tie / tier-3 best hits, eggNOG and diamond-tier calls, family names.
  The line mirrors the `direct` edge, so search never says more than the edge does.
- Reaches about 16K genes. `geneFullText` already indexes the field, so no index change.
  "cyanate" / "urea porter" will find genes.
- Tests: a unit test for the line builder (gate, missing name, missing substrates,
  sanitisation) and a merge test on a fixture calls.json.
- **Needs (your approval):** `prepare_data.sh --steps 2 --force` over all 117 strains.
  It rewrites the committed `gene_annotations_merged.json` files; the diff should be
  description lines only, and I'll check that. Then a Docker rebuild.

## B. Stereo assignment without a ChEBI enantiomer (step 6)

- Rule (in `stereo.chebi_stereo_forms`, or a sibling set read at `paper.py:325`):
  a defined-structure node counts as "one stereo form" when the **same KEGG id** is
  cross-referenced by one ChEBI term with a stereo layer and one term with the same
  connectivity block and **no** stereo layer (`UHFFFAOYSA`). That is, KEGG itself
  files a generic and a specific form under the one id.
  - homocystine C01817 (L,L- + unspecified) → flagged ✓
  - pantothenate C00864 ((R)- + unspecified) → flagged ✓
  - fumarate C00122 (only E forms) → not flagged ✓ — the cis/trans false positive that sank the `is_a` rule
  - xanthine C00385 (tautomers, both unspecified) → not flagged ✓
- The alias-override path (`paper.py` ~line 285) defaults `as_reported` when there are no variants.
  Biller 2022 maps pantothenate by alias, so the override path needs the same check, or
  an explicit `stereo_assignment: dominant_form_assumed` in that alias file. I'd do the
  code check plus an error for overrides that contradict it.
- **Measure before deciding:** a dry run listing every row whose `stereo_assignment`
  would flip. It needs the MNX resolver (`MNX_DATA_DIR`) and the gitignored raw
  `chebi.obo` from the main checkout. I'd run it read-only into the scratchpad.
- Tests: unit cases for the four examples above.
- **Needs (your approval):** `prepare_data.sh --steps 6 --force` (rewrites the committed
  `kegg_data.json`, the identity report and the paper `_resolved.csv` files), then a rebuild.

## C. (optional) Subfamily display names

- For a `tc_subfamily` with no upstream header, use
  `name = "<parent family name> — subfamily <id>"` instead of the bare TCID. This changes
  the documented `t.name = t.tcdb_id` fallback, so it goes under CHANGELOG `### Changed`
  and the CLAUDE.md TcdbFamily bullet.

## Done when

- The backlog bullets are closed or moved as triaged; CHANGELOG `[Unreleased]` has
  entries for A (and B/C).
- `pytest -m "not slow and not kg"` passes in the worktree.
- After approved prepare-data runs, the diffs contain only the intended changes.
  After an approved rebuild, `pytest -m kg` and `/omics-edge-snapshot` show no regressions.
