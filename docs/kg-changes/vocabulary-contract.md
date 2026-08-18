# Controlled-vocabulary contract (`ControlledVocabulary` nodes)

**Date:** 2026-08-18
**Spec:** [`docs/superpowers/specs/2026-08-16-vocabulary-contract-design.md`](../superpowers/specs/2026-08-16-vocabulary-contract-design.md) (rev 5, R1 rescoped 2026-08-17 — see §3 there)
**Driver:** the MCP/explorer's `multiomics_explorer/docs/kg-specs/2026-08-16-interpro-tcdb-asks.md` + follow-up asks
**Track:** graph-wide contract, not a new ontology

## What's changing

Value sets the KG has always enforced internally — which `evidence` strings
are possible on which edge type, which `libraries` InterProScan reports,
which `sources` a `Gene_has_tcdb_family` edge can carry — were only ever
documented as prose (this file, `CLAUDE.md`) or hard-coded on the MCP/explorer
side. When the KG added or renamed a value, the consumer's hard-coded set
silently drifted out of sync and produced wrong (not erroring) answers.

This change publishes those value sets **as data**: one `ControlledVocabulary`
node per (thing-it-applies-to, property) pair, loaded from the same YAML file
by `multiomics_kg/utils/controlled_vocab.py`. A four-gate test suite (unit,
`--test`-build CSV scan, `slow`-build CSV scan, live-graph) checks a `closed`
vocabulary's declared set against what the graph actually contains — an
undeclared emitted value always fails the gate; the reverse direction
(everything declared was actually observed) is a separate, opt-in coverage
check (`exhaustive: true`) run only where the declared set is known to equal
the emitted set. **This is a detection net, not a build-time guard**: the
loader's `VOCAB.check()` helper is wired into exactly one adapter today
(`tcdb_adapter.py`, for `substrate_depth`) — every other vocabulary is emitted
with no loader involvement, so an adapter drifting from the YAML is caught by
the next test run, not blocked at write time.

Along the way, five vocabularies that had accumulated inconsistent naming
(casing, near-duplicate score names, three different two-state encodings)
were aligned to five house rules (§3 below), and four properties that turned
out to restate a threshold over a count the node already publishes were
deleted rather than renamed.

## The contract: `ControlledVocabulary` nodes

Source of truth: `config/controlled_vocabularies.yaml`, loaded by
`multiomics_kg/utils/controlled_vocab.py` and emitted by
`controlled_vocabulary_adapter.py` (node-only, same pattern as the
`DataSource` adapter). One node per declared (`applies_to`, `property`) pair:

```cypher
(:ControlledVocabulary {
  id:              'Gene_has_tcdb_family.evidence_score',
  applies_to:      'Gene_has_tcdb_family',   // node label or edge type
  applies_to_kind: 'edge',                   // 'node' | 'edge'
  property:        'evidence_score',
  value_type:      'float',                  // string | string_array | float
                                              //   | int | bool_string
                                              // NOTE: 'bool' is not admissible — R5
  closed:          'true',                   // "true"/"false" string, not native bool (R5)
  sparse:          'false',
  expected_empty:  'false',
  values:          [],                       // empty here because value_type is
                                              // numeric — see min_value/max_value
  min_value:       0.0,                      // sparse: numeric vocabularies only
  max_value:       1.0,                      // sparse
  signal_count:    5,                        // sparse: evidence_score-shaped only
  signals:         ['eggnog_called', 'source_agreement', 'tier_le_2',
                     'pfam_support', 'go_support'],
  description:     '...'
})
```

For a `string` / `string_array` vocabulary, `values` carries the enumerable
set instead:

```cypher
MATCH (v:ControlledVocabulary {id: 'InterproEntry.interpro_type'})
RETURN v.values
// -> ['FAMILY', 'DOMAIN', 'HOMOLOGOUS_SUPERFAMILY', 'REPEAT', 'CONSERVED_SITE',
//     'ACTIVE_SITE', 'BINDING_SITE', 'PTM']
```

### How to query it

- **"What values can property X take on label/edge Y?"**
  ```cypher
  MATCH (v:ControlledVocabulary {applies_to: 'Gene_has_tcdb_family', property: 'evidence'})
  RETURN v.values, v.closed, v.description
  ```
- **"Is this vocabulary closed (contract-enumerable) or open (enumerate from the graph)?"**
  Read `v.closed`. `"false"` means the value set is genuinely open and
  data-driven (e.g. `Gene.gene_category`, `Experiment.growth_phases`) — do
  not hard-code it, query the graph instead. A `closed: "false"` entry may
  still have `values: []` by design.
- **"Is an empty value set intentional?"**
  Read `v.expected_empty`. `InterproEntry.level_kind` is `closed: "true"`,
  `values: []`, `expected_empty: "true"` — InterPro depth tiers have no
  natural names, so the property is deliberately always null. That is a
  contract, not an oversight.
- **"Has the vocabulary set changed since I last read it?"**
  `Schema_info.controlled_vocabularies_hash` is a sha256 over the emitted
  vocabulary set, stamped at post-import time from the build's
  `config/controlled_vocabularies.yaml`. Compare it release-to-release
  (surfaced by `kg_release_info` on the MCP side) instead of discovering
  drift through a wrong answer.
- **List everything the contract covers:**
  ```cypher
  MATCH (v:ControlledVocabulary) RETURN v.applies_to, v.property ORDER BY 1, 2
  ```

### `applies_to` is per edge type, not per property name

`evidence` and `sources` each get **one `ControlledVocabulary` node per edge
type**, not one shared node — because the domains genuinely differ. Example:
`domain_inferred` is not a possible value of `evidence` on
`Gene_catalyzes_ec_number` (InterPro contributes EC only from FAMILY entries
with a single EC), but it is possible on the three GO edges. That fact is now
structural (absence of the string from that node's `values`), not a caveat
buried in a description field. Likewise `sources` on `Gene_has_cazy_family`
never contains `ncbi` — only the GO edges have an NCBI-sourced component.

## House rules

Five standing conventions, applied to every vocabulary this change touches
and binding for anything the KG adds from here on:

- **R1 — lowercase `snake_case` for values the KG mints; external database
  terms preserved verbatim.** A value that is a controlled term owned by an
  external database (InterPro's `interpro_type`, InterProScan's `libraries` /
  `evalue_library`, NCBIfam's `family_type` including `PfamEq` /
  `PfamAutoEq`) is kept in that database's own casing, so the graph stays
  directly comparable to the source. Everything the KG itself invents
  (`evidence`, `sources`, `substrate_depth`, `level_kind`,
  `MeropsFamily.catalytic_type`, `call_class`, `best_hit_kind`, …) is
  lowercase `snake_case`.
- **R1b — namespace only when values collide across labels.** A value gets a
  namespace prefix only when the same property name holds values from
  different ontologies on different labels — `level_kind` (five labels:
  `tc_family` vs `cazy_family` vs …) needs it, `interpro_type` (one label)
  does not.
- **R2 — every `sources` value is the `id` of a `DataSource` node.**
  Provenance becomes joinable: `gene→ontology edges.sources: 'interpro' ->
  'interproscan'`, `Gene_has_tcdb_family.sources: 'diamond' ->
  'tcdb_diamond'`. Enforced by a kg-validity test — no `sources` value may
  lack a matching `DataSource`.
- **R3 — do not materialize a threshold over a stored count.** A property
  that is only ever `some_count >= N` restates data the node already
  publishes and goes stale silently when the threshold is recalibrated.
  `TcdbFamily.is_promiscuous` and `InterproEntry.is_promiscuous` are deleted
  under this rule — see the derivation recipes below.
- **R4 — one score name per concept, on one scale.** `evidence_score` is now
  a **float in `[0, 1]`** on every annotation edge — fired signals ÷ total
  signals, rounded to 3 decimals — replacing three near-identical names on
  two different integer scales (`evidence_score` 0–3, `tcdb_evidence_score`
  0–5, `tcdb_best_evidence_score`). `signal_count` + `signals` are published
  alongside so `round(score × signal_count)` recovers the original integer
  count.
- **R5 — no native `bool`; a two-state fact is a meaningful categorical
  string.** BioCypher does not round-trip adapter-emitted `bool` properties —
  the one place this shipped (`Interpro_entry_related_to_ec_number` /
  `_cazy_family` `.ambiguous`) was silently `false` on every edge in every
  build ever deployed. `value_type` in the contract admits `string`,
  `string_array`, `float`, `int`, `bool_string` — **not** `bool`; a property
  declared `bool` in `schema_config.yaml` fails the vocabulary test.
  Sentinel-or-absent stays legal for rare-exception flags
  (`is_uninformative`, `level_is_best_effort`) where absence *is* the
  meaning.

## Rename / deletion table

| Object | Property | Before | After | Released before this branch? |
|---|---|---|---|---|
| gene→ontology edges (GO×3, EC, Pfam, CAZy) | `sources` | `'interpro'` | `'interproscan'` | no |
| `Gene_has_tcdb_family` | `sources` | `'diamond'` | `'tcdb_diamond'` | no |
| gene→ontology edges | `evidence_score` | int 0–3 | float 0–1 | no |
| `Gene_has_tcdb_family` | `tcdb_evidence_score` | int 0–5 | renamed to `evidence_score`, float 0–1 | no |
| `Gene` | `tcdb_best_evidence_score` | int 0–5 | renamed to `tcdb_evidence_score_max`, float 0–1 | no |
| `Gene_has_tcdb_family` | `agrees_across_sources` | native `bool` | renamed to `source_agreement`: `both_sources` \| `single_source` | no |
| `Gene_has_tcdb_family` | `pfam_corroborated` | native `bool` | renamed to `pfam_support`: `corroborated` \| `uncorroborated` | no |
| `Gene_has_tcdb_family` | `go_corroborated` | native `bool` | renamed to `go_support`: `corroborated` \| `uncorroborated` | no |
| `Tcdb_family_transports_metabolite` | `substrate_depth` | `deepest` / `ancestor` | `most_specific` / `inherited` | no |
| `InterproEntry` | `interpro_type` | `FAMILY`, `DOMAIN`, … | **unchanged** — external term, kept verbatim under the R1 rescope | no |
| `Gene_has_interpro_entry` | `libraries` / `evalue_library` | `PFAM`, `SUPERFAMILY`, … | **unchanged** — external term, kept verbatim | no |
| `TcdbFamily` | `is_promiscuous` | native `bool` | **deleted** — see derivation below | **yes** (`kg-0.1.0-alpha.6`) |
| `InterproEntry` | `is_promiscuous` | native `bool` | **deleted** — see derivation below | no |
| `Interpro_entry_related_to_{ec_number,cazy_family}` | `ambiguous` | native `bool` (always `false`) | **deleted** — see derivation below | no |
| `Interpro_entry_related_to_{ec_number,cazy_family}` | `source_db` | constant `"interpro.xml"` | **deleted** — hardcoded literal, not data | no |

`TcdbFamily.is_promiscuous` is the one row that touches a released property;
the design spec's pre-change audit of the explorer's own code found zero
consumers of it. Everything else in this table was unreleased.

## Derivation recipes for the four deleted properties

A consumer computing the old value from the current graph:

**`TcdbFamily.is_promiscuous`** (was `metabolite_count >= 50 AND level >= 2`):
```cypher
MATCH (t:TcdbFamily)
WHERE t.level >= 2 AND t.metabolite_count >= 50
RETURN t
```

**`InterproEntry.is_promiscuous`** (was `gene_count >= 1000`):
```cypher
MATCH (n:InterproEntry)
WHERE n.gene_count >= 1000
RETURN n
```

**`Interpro_entry_related_to_{ec_number,cazy_family}.ambiguous`** (was
`len(ecs) > 1 OR interpro_type != 'FAMILY'`, and was uniformly `false` in
every shipped build because of the R5 native-`bool` defect — so there is no
"old value" to reproduce faithfully). The type arm restates a property of the
edge's own source node; the multiplicity arm is the entry's out-degree on
that edge type, and post-pruning out-degree (below) is *more* correct than
what the old flag computed, because the flag was set pre-pruning against
reference EC lists that include obsolete/invalid numbers:
```cypher
MATCH (n:InterproEntry)-[r:Interpro_entry_related_to_ec_number]->()
WITH n, count(r) AS k
WHERE k > 1 OR n.interpro_type <> 'FAMILY'
RETURN n
```

**`Interpro_entry_related_to_{ec_number,cazy_family}.source_db`** (was
always the literal string `"interpro.xml"`) — no derivation needed, it never
varied. Both Layer-A router edge types now carry no properties at all.

## Migration notes for existing queries

- **Any query reading a native `bool` property fails to parse against the new
  graph if it relied on Cypher boolean semantics on a property that is now a
  string.** Concretely: `WHERE r.ambiguous` and `WHERE t.is_promiscuous`
  patterns must be rewritten using the recipes above; `WHERE
  r.source_agreement` (now a string) must become
  `WHERE r.source_agreement = 'both_sources'`, not a bare truthy check.
- **`size(sources) = 2` is still not the right corroboration test** on
  `Gene_has_tcdb_family` — agreement between eggNOG and diamond is
  hierarchical (they often name a family at different depths). Use
  `source_agreement = 'both_sources'` (set by the post-import Cypher using
  the hierarchical check), not list length.
- **Score comparisons across edge types are no longer meaningful as raw
  numbers on different scales** — they were already not comparable before
  this change (0–3 vs 0–5), but now every `evidence_score` is on the same
  `[0, 1]` scale, so a naive `evidence_score > 0.5` filter is at least
  applying the same fraction everywhere. It is still not a calibrated
  probability — read the paired `signal_count` / `signals` to see what
  actually fired.
- **Substrate-depth filters**: replace `substrate_depth = 'deepest'` with
  `substrate_depth = 'most_specific'`, and `'ancestor'` with `'inherited'`.
- **`sources` membership checks**: replace `'interpro' IN r.sources` with
  `'interproscan' IN r.sources`, and `'diamond' IN r.sources` with
  `'tcdb_diamond' IN r.sources`. Prefer `IN` membership over list equality in
  general — it survives a future third source.
- **Do not hard-code any of the value sets in this document long-term** —
  that is precisely the drift this change exists to prevent. Read them from
  `ControlledVocabulary` nodes at startup or on a schedule, and use
  `Schema_info.controlled_vocabularies_hash` to detect when a re-read is due.

## See also

- [`tcdb-two-source-upgrade.md`](tcdb-two-source-upgrade.md) — the TCDB
  two-source contract this change renames properties on (§7.1 there has the
  five chemistry-count definitions)
- [`interpro-multi-ontology.md`](interpro-multi-ontology.md) — the InterPro /
  NCBIfam redesign whose vocabularies are also declared here
- [`interpro-two-layer.md`](interpro-two-layer.md) — the Layer A/B edge
  provenance this change renamed `sources` / `evidence_score` on
