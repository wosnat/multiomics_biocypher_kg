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
                                              //   | int
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

**Per-value descriptions (explorer B1, 2026-08-29).** The trust vocabularies
the explorer filters on — every `*.evidence` and `*.sources` edge vocabulary,
plus `call_class`, `best_hit_kind`, `attachment_depth`, `substrate_depth`,
`pfam_support`, `go_support`, `source_agreement`, `detection_status`,
`table_scope` and `annotation_state` (39 nodes) — additionally carry a sparse
`value_descriptions: str[]`, one `'<value>: <one line>'` element per declared
value in `values` order (Neo4j has no map property; split on the first `': '`):

```cypher
MATCH (v:ControlledVocabulary {id: 'Gene_has_merops_family.call_class'})
UNWIND v.value_descriptions AS d
RETURN split(d, ': ')[0] AS value, substring(d, size(split(d, ': ')[0]) + 2) AS meaning
```

Declared as an optional `value_descriptions` map in the YAML. The loader
rejects a described value that is not in `values`, and a closed vocabulary that
describes some but not all of its values. The key is **not** in
`controlled_vocabularies_hash` — wording can improve without a re-pin. A unit
test (`test_sources_descriptions_do_not_drift_from_gene_annotations_config`)
requires every described `sources` value to be a `logical_sources` id in
`config/gene_annotations_config.yaml`, so the vocabulary text and the
`DataSource` nodes derive from the same source list.

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

  **Recipe** (`multiomics_kg/utils/controlled_vocab.py::vocabularies_hash`,
  written by `create_knowledge_graph.py` to `controlled_vocabularies.sha256`
  in the build output and read by post-import): for every entry, serialise
  `{id, value_type, closed, values (sorted), sparse, expected_empty,
  exhaustive, min_value, max_value, signal_count, signals (sorted)}` with
  `json.dumps(..., sort_keys=True)`; **sort** the resulting strings; join with
  `\n`; sha256 over the UTF-8 bytes; store as `"sha256:" + <64 hex chars>`.
  **Not hashed:** `description` and `applies_to_kind` (`id` is
  `<applies_to>.<property>`, which already pins the target). **Guarantee:**
  nothing build-specific enters the digest — no timestamp, node id, YAML
  order or emission order — so an unchanged vocabulary *set* yields an
  identical string on every rebuild, and a description-only edit does not
  change it. A consumer may pin the string and treat a mismatch as "re-read
  the vocabulary nodes" (the explorer folds this into its
  `kg_release_info` verdict as `warn`). **Enforced at release:** `/release-kg`
  Phase 5 recomputes the hash from the tag checkout and refuses to release a
  staged graph whose stamp is null or different; `metadata.json` on every
  `kg-*` GitHub Release carries `controlled_vocabularies: {hash, entry_count,
  entry_ids}`, and the release notes' "What changed since" block lists the
  added/removed entries whenever the hash moved. Expect it to move on most
  data releases — closed vocabularies gain values as papers land.
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
- **R2 — every `sources` value corresponds to a `DataSource` node.**
  Provenance becomes joinable: `gene→ontology edges.sources: 'interpro' ->
  'interproscan'`, `Gene_has_tcdb_family.sources: 'diamond' ->
  'tcdb_diamond'`. The join is prefixed — a `sources` value `s` matches the
  `DataSource` node whose `id` property is `'data_source:' + s` (e.g.
  `sources` value `'tcdb_diamond'` joins `DataSource {id:
  'data_source:tcdb_diamond'}`), not `d.id = s`. Enforced by a kg-validity
  test — no `sources` value may lack a matching `DataSource` under that
  prefixed join.
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
  `string_array`, `float`, `int` — **not** `bool` (and, since 2026-08-28, no `bool_string`: every two-state fact is a named pair, see `two-state-strings.md`); a property
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
| `TigrRole` | `level_kind` | (did not exist — the label was flat) | **new**: `tigr_mainrole` \| `tigr_subrole`, closed (R1 KG-minted snake_case; R1b-clean — `level_kind` is the cross-label hierarchy-position property, same role as on TCDB/CAZy/MEROPS/BRITE) | no |
| `Gene_has_tigr_role` | `sources` | `[cyanorak]` | **widened**: `[cyanorak, interproscan]`, closed (R2 — both join a `DataSource`) | no |
| `Gene_has_tigr_role` | `evidence` | `[curated]` | **widened**: `[curated, family_inferred]` — two rungs of the one KG-SYNC-005 ladder | no |
| `Ncbifam_family_has_tigr_role` | — | (new edge type) | **no properties at all** (R3/R5): the archive is one frozen source, so provenance is documented on the edge type, not repeated per edge; there is no two-state fact to name | no |

The last four rows are additions/widenings rather than renames or deletions —
they are listed here because the table is where a consumer looks to find out
whether a value set it hard-coded is still complete. Nothing is renamed or
removed by them; a query reading `Gene_has_tigr_role` still parses, it just
now sees values it did not before (see Migration notes).

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
- **`Gene_has_tigr_role` is the worked example of that last point** (2026-08-29):
  it was a single-source edge (`['cyanorak']`) and now carries a second source,
  so `r.sources = ['cyanorak']` silently drops every edge where the curated
  Cyanorak role and the inferred equivalog-family role agree — exactly the
  edges a curated-only query most wants. Use `'cyanorak' IN r.sources`. Its
  `evidence` likewise gains `family_inferred`, so `evidence = 'curated'` is now
  a *filter* rather than a tautology, and a query that wants everything must
  stop asserting the old single value.
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
- [`tigr-role-bridge.md`](tigr-role-bridge.md) — the TigrRole hierarchy /
  NCBIfam role bridge whose `level_kind` and widened
  `Gene_has_tigr_role.sources` / `.evidence` are declared here

## KG-SYNC-005 addendum (2026-08-27) — one evidence ladder, everywhere

`evidence` and `sources` now exist on all **14** gene→ontology edge types (previously 7 / 6). The ladder
gains one rung — **`curated > signature > homology > family_inferred > domain_inferred`** — where
`homology` is a direct sequence-similarity hit vs. a curated reference DB (diamond: TCDB, MEROPS);
strength *within* `homology` is read from the edge's `tier`, never split into further ladder values.
Per-edge `ControlledVocabulary` entries keep declaring their own possible **subset** (spec §5.2 stands).
New R3-exempt structural pair `Gene_has_tcdb_family.attachment_depth: most_specific | superseded`. Renames:
`Gene_has_ncbifam_family.score` → `bit_score` (R4 spirit — `score` is PSORTb's scale);
`MeropsFamily.family_type` → `family_class` (R1b collision with `NcbifamFamily.family_type`).
`NcbifamFamily.family_type` gains the single KG-minted sentinel `retired` inside an otherwise
external-verbatim set (R1 exception, documented on the entry). Full contract:
`docs/kg-changes/annotation-trust-surface.md`.
