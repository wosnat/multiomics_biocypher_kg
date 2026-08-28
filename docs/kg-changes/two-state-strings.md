# Two-state facts are named string pairs — the last `"true"`/`"false"` properties converted (2026-08-28)

**Status:** code landed 2026-08-28; needs a Docker rebuild to go live. **BREAKING** for
every consumer that compared against the literal strings `'true'` / `'false'`.
Explorer hand-off: [`2026-08-28-explorer-handoff.md`](2026-08-28-explorer-handoff.md).

## What changed

House rule **R5** of the vocabulary contract (`vocabulary-contract.md`) forbids native
`bool` and deprecates stringified booleans: a two-state fact is a categorical string
that names *both* states, so a reader never has to remember which side of a bool a
property was (`substrate_depth: most_specific | inherited`,
`has_cross_genus_members: cross_genus | single_genus`). Seven released properties
predated R5 and were grandfathered as `value_type: bool_string`; an eighth
(`Assay_flags_metabolite.flag_value`) had the same shape and no vocabulary entry at
all. All eight are now named pairs. **Property names are unchanged** — only the values.

| Property | Was | Now | Live counts (pre-rebuild, 2026-08-28) |
|---|---|---|---|
| `Experiment.is_time_course` | `"true"`/`"false"` | `time_course` \| `single_time_point` | 41 / 168 |
| `Experiment.reports_fold_change` (post-import) | `"true"`/`"false"` | `fold_change` \| `no_fold_change` | 174 / 35 |
| `DerivedMetric.rankable` | `"true"`/`"false"` | `rankable` \| `not_rankable` | 42 / 41 |
| `DerivedMetric.has_p_value` | `"true"`/`"false"` | `p_value` \| `no_p_value` | 0 / 83 |
| `MetaboliteAssay.rankable` | `"true"`/`"false"` | `rankable` \| `not_rankable` | 12 / 2 |
| `Derived_metric_quantifies_gene.significant` (post-import, sparse) | `"true"`/`"false"` | `significant` \| `not_significant` | none set today (no DM has `p_value`) |
| `Derived_metric_flags_gene.value` | `"true"`/`"false"` | `flagged` \| `not_flagged` | 8,126 / 3,773 |
| `Assay_flags_metabolite.flag_value` | `"true"`/`"false"` (undeclared) | `detected` \| `not_detected` | 58 / 128 |

Why these names: `rankable`/`p_value` read as the fact they assert; `significant |
not_significant` reuses the token `Changes_expression_of.expression_status` already
publishes; `detected | not_detected` is `Assay_quantifies_metabolite.detection_status`
minus its `sporadic` state (a single presence call cannot be sporadic); `flagged |
not_flagged` is deliberately neutral because a DM boolean's polarity is defined by the
paper's `true_tokens` (`periodic_in_axenic_LD`, `exoproteome_detected_ancestor`, …).

## What did NOT change

- **The paperconfig author contract.** `rankable: "true"`, `has_p_value: "true"`,
  `true_tokens` / `false_tokens` / `blank_policy: "true"|"false"` keep their wording
  and `scripts/validate_paperconfig.py` still enforces `{"true","false"}`. The
  adapters map at emit time (`observations_adapter.RANKABLE / HAS_P_VALUE /
  FLAG_VALUE`, `metabolite_assay_adapter.RANKABLE / FLAG_VALUE`) and raise on
  anything else. No paperconfig was edited for this change.
- Every gate that reads these values: post-import ranks only `rankable = 'rankable'`
  parents, derives `significant` only under `has_p_value = 'p_value'`, counts
  `flag_true_count` from `flag_value = 'detected'`. Same rows, same numbers.
- `Schema_info.git_dirty` (build metadata, `"true"`/`"false"`) is out of scope.

## Explorer questions answered (R1, R2 — 2026-08-28)

- **R1 — `Derived_metric_flags_gene.value = 'not_flagged'` DOES occur.** The
  "positive-only" picture is stale: it held for the Biller 2018 / Coe 2016 / Biller 2014
  DMs (`flag_false_count = 0`), but 11 of the 27 boolean DMs store both states — live,
  3,773 `not_flagged` vs 8,126 `flagged` edges. Per DM (false / true): Voigt 2014
  `has_primary_tss` MED4 805/824 · MIT9313 1,110/1,029; Steglich 2010
  `expressed_above_background` 920/1,010; Biller 2022 (aem.00798-26) six
  `*_detected_*` metrics 59–205 false each; Hennon 2015 (ismej.2015.36) two
  `rapid_recovery_*` 38/13 and 46/5. So `genes_by_boolean_metric(flag=False)` returns
  rows for those DMs; `DerivedMetric.flag_false_count` tells you per DM whether the
  paper reported the negative state. Same as before — only the token changed.
- **R2 — all eight are `ControlledVocabulary` entries** (`closed: true`,
  `value_type: string`; `Derived_metric_quantifies_gene.significant` is `sparse`),
  so `list_filter_values` can serve them and `test_observed_values_are_declared`
  covers them. Pin the entries, not just the literals.

## Vocabulary

The eight entries in `config/controlled_vocabularies.yaml` (section "Two-state facts")
are `value_type: string`, `closed: true`; `significant` is `sparse: true`. The
`bool_string` value type is **removed** from `controlled_vocab.VALUE_TYPES` — nothing
may declare one again. `Schema_info.controlled_vocabularies_hash` changes.

## Where the strings live (KG side)

`omics_adapter.py` (`is_time_course`), `observations_adapter.py`,
`metabolite_assay_adapter.py`, `scripts/post-import.{sh,cypher}` (5 sites, identical
in both files), `scripts/post-import-validate.sh`, `config/schema_config.yaml`
comments, `tests/kg_validity/{test_derived_metric,test_numeric_derived_metric}.py`,
`tests/kg_validity/snapshot_data.json` (37 values rewritten in place),
`.claude/skills/{cypher-queries,paperconfig}/SKILL.md`.

## Verification (post-rebuild)

```cypher
MATCH (e:Experiment) RETURN e.is_time_course, e.reports_fold_change, count(*);
MATCH (d:DerivedMetric) RETURN d.rankable, d.has_p_value, count(*);
MATCH ()-[r:Derived_metric_flags_gene]->() RETURN r.value, count(*);
MATCH ()-[r:Assay_flags_metabolite]->() RETURN r.flag_value, count(*);
```
Expect exactly the pairs above, no `'true'` / `'false'`, and the same counts as the
table. Then `pytest -m kg` (the enum tests assert the pairs) and
`tests/kg_validity/test_controlled_vocabularies.py` (observed ⊆ declared).
