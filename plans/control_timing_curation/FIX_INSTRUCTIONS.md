# Paperconfig fix pass (2026-10-18) — instructions for fix agents

The control-timing review (`plans/control_timing_curation.md`) found metadata errors in several
paperconfigs. The curator asked to **double check each one against the paper and apply it**. You own
a set of papers (given in your prompt); touch no other paper's files.

## Rules

1. **Verify before changing.** For every item, re-read the relevant part of the paper (main PDF via
   the Read tool with `pages`, supplementary files, table legends, CSV headers). Apply the change only
   if the paper confirms it. If the paper contradicts the claim or is ambiguous, change nothing and
   report it as `not applied` with the evidence.
2. **Minimal edits.** Use the Edit tool with small exact replacements; keep comments, key order,
   quoting and indentation. Never reformat the file.
3. **Do not change**: experiment keys, analysis `id`s, `name_col`, `id_columns`, filenames,
   `control_timing` (unless your prompt says so), or anything not on your list. Do not commit, do
   not run `prepare_data.sh` or Docker.
4. **Vocabularies.** `treatment_type`, `background_factors`, `growth_phase`, `omics_type`,
   `compartment`, `test_type` must stay canonical (see `scripts/validate_paperconfig.py` and
   `config/controlled_vocabularies.yaml`). `treatment_type` and `background_factors` must be
   non-empty lists. A `viral` or `coculture` `treatment_type` requires `treatment_organism`.
5. **Strings**: no `|` and no `'` (apostrophe) in any value you write.
6. **design_notes** (free text on an experiment): when you add a sentence, append it to the existing
   note, or create `design_notes: >-` as the last key of the experiment block. Wrap folded lines
   without breaking a hyphenated word across lines (`nutrient-replete` must stay on one line).
   Re-read the parsed value with `yaml.safe_load` to confirm it is what you meant.
7. **Raw p-values (item E3)** — the curator chose "keep with a note": on each listed experiment, add
   this sentence to `design_notes`: `adjusted_p_value on these edges is the unadjusted p-value
   reported by the paper (no multiple-testing correction).` Confirm in the paper first.
8. **Validate**: `uv run python scripts/validate_paperconfig.py "<path>"` must still end in
   `VALIDATION PASSED`, with no new errors.
9. **Graph-visible changes** (timepoint labels, `timepoint_hours`, `growth_phase`, `omics_type`,
   `compartment`, conditions, `fold_change_type`): list them; also `grep -c` the experiment ids of
   your papers in `tests/kg_validity/snapshot_data.json` and report any snapshot entry whose
   properties your change alters (do not edit the snapshot).

## Report — `plans/control_timing_curation/fixes_<batch>.md`

```markdown
# Fix batch <X>

## <Paper>
| Item | Claim | Verified? (quote, page) | Change applied (before → after) |
|---|---|---|---|

Not applied / needs the curator: <list with reason>
Graph-visible changes: <list>
Snapshot entries affected: <list or none>
CHANGELOG line: <one sentence: what the KG asserted before, what it asserts now>
```

When done, reply with only: the report path, a per-paper one-line summary, and the
"not applied / needs the curator" items.
