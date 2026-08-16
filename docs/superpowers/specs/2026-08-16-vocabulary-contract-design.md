# Controlled-vocabulary contract + cross-ontology vocabulary alignment

**Date:** 2026-08-16
**Status:** Design — approved in brainstorming, not yet implemented
**Driver:** `multiomics_explorer/docs/kg-specs/2026-08-16-interpro-tcdb-asks.md`
(KG-IPT-001 … 007)
**Verified against:** KG `0.0.0-dev`, `built_at 2026-08-13T12:19:46.858Z`; last
release tag `kg-0.1.0-alpha.6`

---

## 1. Problem

The InterPro + TCDB two-source integration introduced five new controlled
vocabularies. They were each designed in isolation, against their own source
database's conventions, and the result is not uniform:

| # | Inconsistency | Where |
|---|---|---|
| 1 | `interpro_type` values are `UPPERCASE`; every other categorical in the KG is lowercase `snake_case` | `InterproEntry` |
| 2 | One provider, three spellings: edge `sources` says `interpro`, `DataSource` / `Gene.contributing_sources` say `interproscan`; the TCDB edge says `diamond`, `DataSource` says `tcdb_diamond` | gene→ontology edges, `Gene_has_tcdb_family` |
| 3 | `is_promiscuous` means *many substrates* on `TcdbFamily` but *many genes* on `InterproEntry` — one name, two axes | `TcdbFamily`, `InterproEntry` |
| 4 | Three near-identical score names on two different scales: `evidence_score` (0–3), `tcdb_evidence_score` (0–5), `tcdb_best_evidence_score` | annotation edges, `Gene` |
| 5 | `substrate_depth: {deepest, ancestor}` mixes a superlative with a structural relation, and `ancestor` describes the *node* while the property describes a *(node, substrate) fact* | `Tcdb_family_transports_metabolite` |

None of these are defects — the pre-flight in the asks doc ran 37 invariant
checks and found zero. They are naming and contract problems, and they are
cheapest to fix now, before the vocabularies are baked into a published MCP
enum.

Separately, and more importantly: **none of these vocabularies are published as
data.** The explorer hard-codes all of them in `list_filter_values` and in
Pydantic `Literal` types. When the KG adds a value, the MCP does not error — it
silently returns a wrong answer.

---

## 2. Scope

### 2.1 What is frozen

Verified against `kg-0.1.0-alpha.6` and against the explorer's actual reads
(`multiomics_explorer/mcp_server/tools.py`, `api/functions.py`):

- `level_kind` values (`tc_class` … `tc_specificity`, `cazy_*`, `brite_*`, KEGG's
  `category|subcategory|pathway|ko`). The explorer's entire TCDB surface derives
  from these — `transport_confidence` is computed from `tc_family` vs
  `tc_subfamily`/`tc_specificity`.
- `annotation_quality` (0–3) and `annotation_state`.
- `DataSource` node ids and `Gene.contributing_sources`.
- The existence of `Gene_has_tcdb_family`.

### 2.2 What is free

Everything the integration introduced. Confirmed unreleased (absent from
`kg-0.1.0-alpha.6:config/schema_config.yaml` and `:scripts/post-import.cypher`)
**and** unread by the MCP:

`sources` · `evidence` · `evidence_score` on gene→ontology edges ·
`Gene_has_tcdb_family.{sources, tier, tcdb_evidence_score, agrees_across_sources,
pfam_corroborated, go_corroborated}` · `Gene.tcdb_best_evidence_score` ·
`InterproEntry.{interpro_type, is_promiscuous}` · `substrate_depth` ·
`Gene.transport_substrate_resolution` · Layer-A `{ambiguous, source_db}`.

`TcdbFamily.is_promiscuous` shipped in alpha.6 but is read by nothing in the
explorer, so renaming it costs one `### Breaking` bullet and no consumer work.

### 2.3 Non-goals

- Declaring *every* vocabulary in the graph. v1 seeds the contract with the
  integration vocabularies plus the ones the explorer already hard-codes
  (§5.3); the mechanism is extensible and later releases add to it.
- Relationship-property indexes (KG-IPT-008 — explicitly deferred by the
  explorer).
- Release-process items (`mcp_min_version`, `deployment_role`, build pinning).
- Any change to `annotation_quality` / `annotation_state`.

---

## 3. House rules

Four rules, each resolving a numbered inconsistency from §1. These are the
standing conventions for every vocabulary the KG adds from here on.

### R1 — lowercase `snake_case`; namespace only when values collide across labels

A value carries a namespace prefix **only** when the same property name holds
values from different ontologies across labels. `level_kind` lives on five
labels, so `tc_family` must be distinguishable from `cazy_family`.
`interpro_type` lives on one label, so its values are bare.

```
interpro_type:  FAMILY, DOMAIN, HOMOLOGOUS_SUPERFAMILY, REPEAT,
                CONSERVED_SITE, ACTIVE_SITE, BINDING_SITE, PTM
             -> family, domain, homologous_superfamily, repeat,
                conserved_site, active_site, binding_site, ptm
```

### R2 — every `sources` value is the `id` of a `DataSource` node

Provenance becomes joinable, and the three-spellings problem disappears:

```
gene→ontology edges:      'interpro'  -> 'interproscan'
Gene_has_tcdb_family:     'diamond'   -> 'tcdb_diamond'
```

Enforced by a kg-validity test: no `sources` value may lack a matching
`DataSource`.

### R3 — breadth flags are named for *what* is broad

```
TcdbFamily.is_promiscuous    -> is_multi_substrate
  (many distinct substrates => inferring a member gene's cargo is weak)
InterproEntry.is_promiscuous -> is_multi_gene
  (matched by a large share of genes => unstratified ORA is invalid)
```

Both remain advisory. Neither is ever a default filter — this is the confirmed
answer to KG-IPT-006, and the rename makes the intent legible without reading a
doc.

**Why not `is_informative`.** The KG already has a released informativeness
concept, and it is a *different axis*. `is_uninformative`
(`config/uninformative_terms.yaml`) is **curated content-emptiness** — "conveys
no class signal at all": the three GO roots, `cog.category:S`,
`cyanorak.role:R.4`. Its stated principle is that a term naming the broad class
stays informative even when the sub-class is unknown, which is why DUF/UPF Pfams
and "Transport > Unknown substrate" are deliberately not flagged.

`is_multi_gene` is **measured frequency**. Its members are the opposite of
content-free: IPR027417 (P-loop NTPase superfamily) is highly informative about
fold and function, and is flagged only because it appears in 6,909 genes and so
cannot discriminate in an ORA. The two are independent — informative-but-
ubiquitous and uninformative-but-rare are both real quadrants.

Naming the frequency flag `is_informative` would make it read as the negation of
a released flag that means something else, and would assert that P-loop NTPase
carries no functional content. `is_multi_gene` names the axis that is actually
broad and sits unambiguously beside `is_uninformative`.

### R4 — one score name per concept, on one scale

`evidence_score` is a **float in `[0,1]`** on every annotation edge: fired
signals ÷ total signals, rounded to 3 decimals.

| Edge | Signals | Values |
|---|---|---|
| `Gene_has_pfam`, `Gene_catalyzes_ec_number`, `Gene_involved_in_biological_process`, `Gene_enables_molecular_function`, `Gene_located_in_cellular_component`, `Gene_has_cazy_family` | 3 | `0.0, 0.333, 0.667, 1.0` |
| `Gene_has_tcdb_family` | 5 | `0.0, 0.2, 0.4, 0.6, 0.8, 1.0` |

```
Gene_has_tcdb_family.tcdb_evidence_score -> evidence_score   (float, /5)
Gene.tcdb_best_evidence_score            -> tcdb_evidence_score_max (float)
```

**Why normalize rather than publish two integer ranges.** The integer is not
lost: the contract publishes `signal_count`, so `0.6 × 5 = 3` recovers it
exactly. What normalizing buys is the failure mode of a consumer who never reads
the contract — under two scales they compare a Pfam `3` to a TCDB `3` and are
wrong by a factor; normalized, they degrade to "roughly comparable fraction".
Neither scale is calibrated across edge types and the contract says so
explicitly.

`annotation_quality` stays an integer 0–3: it is released, the MCP reads it, and
it is an ordinal *state* (`no_evidence` → `informative_multi`), not a signal
fraction.

### Applied to the remaining two vocabularies

```
substrate_depth: {deepest, ancestor} -> {most_specific, inherited}
```
`ancestor` describes the node; the property describes a *(node, substrate)* fact
— `2.A.1` is `most_specific` for one substrate and `inherited` for another.
`inherited` states what is actually true of the fact.

```
Gene.transport_substrate_resolution: {resolved, family_inferred}   UNCHANGED
```
This answers **KG-IPT-002: yes, `resolved` is final.** The explorer retires
`substrate_confirmed`. `family_inferred` is deliberately a shared house token
meaning "inherited from a broader family rather than pinned", used identically
on the `evidence` axis — the reuse is the uniformity, not a collision.

---

## 4. Rename summary (the explorer's migration list)

| Object | Property | Before | After | Released? |
|---|---|---|---|---|
| `InterproEntry` | `interpro_type` | `FAMILY`, `DOMAIN`, … | `family`, `domain`, … | no |
| gene→ontology edges | `sources` | `interpro` | `interproscan` | no |
| `Gene_has_tcdb_family` | `sources` | `diamond` | `tcdb_diamond` | no |
| `TcdbFamily` | `is_promiscuous` | — | `is_multi_substrate` | **yes** |
| `InterproEntry` | `is_promiscuous` | — | `is_multi_gene` | no |
| gene→ontology edges | `evidence_score` | int 0–3 | float 0–1 | no |
| `Gene_has_tcdb_family` | `tcdb_evidence_score` | int 0–5 | `evidence_score`, float 0–1 | no |
| `Gene` | `tcdb_best_evidence_score` | int 0–5 | `tcdb_evidence_score_max`, float 0–1 | no |
| `Tcdb_family_transports_metabolite` | `substrate_depth` | `deepest` / `ancestor` | `most_specific` / `inherited` | no |

Exactly one row is breaking, and it has no consumer.

---

## 5. The contract: `ControlledVocabulary` nodes

### 5.1 Node shape

Source of truth is `config/controlled_vocabularies.yaml`, loaded by
`multiomics_kg/utils/controlled_vocab.py` and emitted by a node-only adapter —
the same pattern `config/gene_annotations_config.yaml` already uses to drive the
`DataSource` nodes.

```
(:ControlledVocabulary {
   id:              'Gene_catalyzes_ec_number.evidence',
   applies_to:      'Gene_catalyzes_ec_number',
   applies_to_kind: 'edge',              // edge | node
   property:        'evidence',
   value_type:      'string',            // string | string_array | float | int
                                         //   | bool | bool_string
   closed:          true,                // see 5.2
   values:          ['curated', 'family_inferred'],
   sparse:          false,               // property may be absent on some rows
   expected_empty:  false,               // declared, but no values by design
   description:     '...'
})
```

Numeric vocabularies add `min_value` / `max_value`, and scores add
`signal_count` + `signals`:

```
(:ControlledVocabulary {
   id:           'Gene_has_tcdb_family.evidence_score',
   value_type:   'float',
   min_value:    0.0,
   max_value:    1.0,
   signal_count: 5,
   signals:      ['eggnog_called', 'agrees_across_sources', 'tier_le_2',
                  'pfam_corroborated', 'go_corroborated'],
   description:  'Advisory ranking, never a filter. Fraction of independent
                  corroborating signals that fired; multiply by signal_count
                  for the raw count. NOT calibrated against other edge types
                  — compare within an edge type. Float: sort and threshold,
                  never test equality.'
})
```

Publishing `signals` makes machine-readable what is currently prose in
`tcdb-two-source-upgrade.md`. The three component booleans
(`agrees_across_sources`, `pfam_corroborated`, `go_corroborated`) are already
stored beside the score precisely so consumers can re-weight; now they can
discover the component names.

`Schema_info.controlled_vocabularies_hash` carries a sha256 over the emitted
vocabulary set, so `kg_release_info` detects drift instead of the explorer
discovering it through a wrong answer.

### 5.2 Three properties that do real work

**`applies_to` is per edge type, so `evidence` subsetting is structural.** The
asks doc's central point — that `domain_inferred` is *not a possible value* on
`Gene_catalyzes_ec_number` — is not a caveat in a description field; it is the
absence of that string from that node's `values`. Measured on the deployed
build:

| Edge type | `evidence` domain |
|---|---|
| `Gene_has_pfam` | `curated`, `signature` |
| `Gene_catalyzes_ec_number` | `curated`, `family_inferred` |
| `Gene_involved_in_biological_process` (and MF / CC) | `curated`, `family_inferred`, `domain_inferred` |
| `Gene_has_cazy_family` | `curated`, `domain_inferred`, `family_inferred` |

**`closed` separates enumerable-by-contract from enumerable-at-runtime.** Some
of what `list_filter_values` serves is genuinely closed (`omics_type`,
`value_kind`, `compartment`); some is open and data-driven (`gene_category`
derives from role annotations, `growth_phase` from paperconfig free text).
Declaring `closed: false` with an empty `values` tells the explorer to enumerate
from the graph rather than hard-code, and stops drift tests from firing on
values that are *supposed* to grow.

**`expected_empty` answers KG-IPT-003 as data rather than prose.** `InterproEntry`
ships:

```yaml
InterproEntry.level_kind:
  value_type: string
  closed: true
  values: []
  expected_empty: true
  description: >
    InterPro depth tiers have no natural names. Emit null. Stratify ORA by
    (interpro_type, level) — interpro_type primary (breadth), level secondary.
```

`ontology_landscape` then emits null **on contract**, not by accident, and does
not have to invent a value.

### 5.3 v1 seed

Two groups. Everything else is a follow-up, added as later releases touch it.

**Integration vocabularies** (all values final per §3–4): gene→ontology
`sources` / `evidence` / `evidence_score`; `Gene_has_tcdb_family` `sources` /
`tier` / `evidence_score`; `Gene.tcdb_evidence_score_max`;
`InterproEntry.interpro_type` / `level_kind` (empty) / `is_multi_gene`;
`TcdbFamily.is_multi_substrate`; `Tcdb_family_transports_metabolite.substrate_depth`;
`Gene.transport_substrate_resolution`; Layer-A `ambiguous` / `source_db`.

**Vocabularies the explorer already hard-codes**, i.e. the eight
`list_filter_values` filters plus their neighbours: `omics_type`, `value_kind`,
`compartment`, `metric_type`, `brite_tree`, `Metabolite.evidence_sources`,
`expression_status`, `metric_bucket`, `detection_status`, `annotation_state`,
`level_kind` (per label), `organism_type`, `table_scope`, `treatment_type`,
`background_factors`, plus the open `gene_category` and `growth_phase`.

Exact value sets for the second group are harvested from the code and the live
graph during implementation and reviewed against `CLAUDE.md` before landing —
this spec does not restate them, to avoid becoming a second source of truth for
values it did not verify.

---

## 6. Enforcement

Adapters import their literals from the loader, so an undeclared value cannot be
emitted. Post-import Cypher cannot import Python, so its literals are covered by
the graph-level test only. Four gates, each catching something the others
structurally cannot:

| Gate | Marker | Runs | Checks |
|---|---|---|---|
| Adapter units — `VOCAB.check()` raises at emit | none | always | undeclared value, in seconds |
| `--test` build → CSV scan | none | always | observed ⊆ declared, pre-import, no Neo4j |
| Full build → CSV scan | `slow` | opt-in | + declared − observed ⊆ `expected_empty` |
| Live graph | `kg` | every rebuild | both directions · post-import-only properties · R2's `DataSource` join · `controlled_vocabularies_hash` |

`tests/test_create_knowledge_graph.py` is restructured to parameterize one
helper over the two build modes:

```
_run_build(mode) -> output_dir            # module-scoped, one build per mode

test_build_no_errors[test]                # --test, fast, default run
test_vocabularies_in_csvs[test]

test_build_no_errors[full]                # @pytest.mark.slow, behaviour unchanged
test_vocabularies_in_csvs[full]
```

This also moves the existing ERROR-line smoke check into the default
`pytest -m "not slow and not kg"` run, where today it only exists behind an
opt-in hour-long test — a gain independent of this change.

Two honest limitations. Both build gates require a populated `cache/`, exactly
as today's slow test does, so neither is portable to a fresh CI checkout. And
`--test` mode (100 items/adapter) cannot contain rare values — `interpro_type:
ptm` has 7 entries graph-wide — so only the `slow` and `kg` gates can run the
declared − observed direction. If `--test` measures in minutes rather than
seconds it moves behind its own marker; this is measured during implementation,
not assumed here.

---

## 7. Answers to the explorer's asks

| ID | Answer |
|---|---|
| **KG-IPT-001** | Implemented as `ControlledVocabulary` nodes (§5), versioned by `Schema_info.controlled_vocabularies_hash`. `evidence` is per edge type structurally. |
| **KG-IPT-002** | **Yes — `resolved` is final.** Retire `substrate_confirmed`. |
| **KG-IPT-003** | Null-by-design, shipped as a `ControlledVocabulary` node with `expected_empty: true`. Stratify by `(interpro_type, level)`. |
| **KG-IPT-004** | Five one-line definitions added to `docs/kg-changes/tcdb-two-source-upgrade.md` (§7.1). |
| **KG-IPT-005** | **Premise is incorrect — see §7.2.** |
| **KG-IPT-006** | Confirmed: advisory, never a default filter. R3 renames it to say so. |
| **KG-IPT-007** | CHANGELOG prose corrected to `control`. Convention adopted: prose naming a graph property uses the exact property name. |

### 7.1 The five chemistry counts (KG-IPT-004)

To be added verbatim beside the existing property tables in
`tcdb-two-source-upgrade.md`, and quotable in MCP field descriptions:

- `Gene.metabolite_count` — distinct metabolites reachable from this gene by
  **catalysis only** (`Gene_catalyzes_reaction` → `Reaction_has_metabolite`).
- `Gene.transported_metabolite_count` — distinct metabolites reachable by
  **transport only**, counting the gene's **deepest** TCDB attachments so an
  ancestor's rollup is not inherited alongside its own descendant.
- `Metabolite.gene_count` — distinct genes reaching this metabolite by
  **catalysis only**. The catalysis-arm mirror of `Gene.metabolite_count`.
- `Metabolite.transporter_gene_count` — distinct **genes** reaching this
  metabolite by transport, same deepest-attachment predicate. Agrees with
  `Gene.transported_metabolite_count` by construction — two projections of one
  (gene, metabolite) set.
- `Metabolite.transporter_count` — distinct **transporter systems**
  (`TcdbFamily` nodes), not genes, whose substrate edge carries
  `substrate_depth = 'most_specific'`.

The explorer's stated understanding of all five is correct.

### 7.2 KG-IPT-005 — InterPro GO xrefs already landed

The ask asks to confirm that InterPro → GO / pathway xrefs are **not** landing
this release, so the explorer does not design for them. The GO half of that is
wrong, and the correction changes what the explorer must build.

InterPro's entry-level GO xrefs **did land**, as Layer B enrichment on the
**existing** `Gene_involved_in_biological_process` / `_enables_molecular_function`
/ `_located_in_cellular_component` edges — see
`multiomics_kg/download/build_gene_annotations.py:397`
(`_fold_interpro_field(gene, "go_terms", go_new)`), gated to FAMILY + DOMAIN
entries. It contributes ~45K net-new GO terms.

The conclusion came from `db.relationshipTypes()`, which cannot see it: no new
relationship type was created. The evidence is already in the asks doc's own
KG-IPT-001 table — `Gene_involved_in_biological_process` carries
`family_inferred` 7,918 and `domain_inferred` 7,782, and `interpro` appears in
its `sources`. Those rows *are* the InterPro GO contribution.

What **is** deferred: InterPro's MetaCyc **pathway** xrefs (no pathway edges from
InterPro at all; InterPro ships no KEGG xrefs, verified 0 occurrences in the
2026-08-06 release), and any `Interpro_entry_related_to_go` Layer-A router.

Consequence for the explorer: the provenance-surfacing design was about to be
deferred on the belief that GO provenance would change later. It will not — the
`sources` / `evidence` shape on GO edges is final this release and should be
built against now.

---

## 8. Implementation surface

**New**
- `config/controlled_vocabularies.yaml`
- `multiomics_kg/utils/controlled_vocab.py` — loader + `VOCAB.check()`
- `multiomics_kg/adapters/controlled_vocabulary_adapter.py` — node-only,
  modelled on the `DataSource` adapter
- `tests/test_controlled_vocab.py`
- `tests/kg_validity/test_controlled_vocabularies.py`

**Changed**
- `create_knowledge_graph.py` — register the adapter; write the vocabulary hash
  to the output dir for post-import to stamp
- `config/schema_config.yaml` — `ControlledVocabulary` node type; renamed
  properties
- `multiomics_kg/utils/annotation_provenance.py` — normalized `evidence_score`
- `multiomics_kg/adapters/{interpro_adapter,tcdb_adapter}.py`,
  `multiomics_kg/download/build_gene_annotations.py` — renamed literals, import
  from the loader
- `scripts/post-import.sh` + `scripts/post-import.cypher` — renames, normalized
  `tcdb_evidence_score_max`, hash stamp in Group 4. **Kept byte-identical to
  each other.**
- `tests/test_create_knowledge_graph.py` — parameterized over build mode (§6)
- `tests/kg_validity/generate_snapshot.py` output — regenerated
- `docs/kg-changes/{interproscan-extension,interpro-two-layer,tcdb-two-source-upgrade}.md`
- `CLAUDE.md`, `CHANGELOG.md`

**Validation gate.** `scripts/post-import-validate.sh` baseline captured against
the currently deployed graph *before* any change; after the rebuild the diff must
be empty except for the renamed properties and the normalized score columns —
every count and every other value byte-identical. Then
`pytest -m "not slow and not kg"`, Docker rebuild, `pytest -m kg`, and
`/omics-edge-snapshot`.

**CHANGELOG.** One `### Breaking` bullet (`is_promiscuous` → `is_multi_substrate`;
`evidence_score` integer → normalized float), one `### Added` bullet for the
contract, plus the KG-IPT-007 prose correction to the existing biller 2016 entry.

**Reply to the explorer.** A short doc answering all seven asks, pointing at §7
and at the rename table in §4, and flagging §7.2 as a correction that unblocks
work they were about to defer.

---

## 9. Open items

1. `create_knowledge_graph.py --test` wall-clock — measured during
   implementation; if minutes, the fast gate moves behind its own marker.
2. Exact value sets for the §5.3 second group, harvested and reviewed before
   landing.
3. Whether `is_multi_gene` / `is_multi_substrate` thresholds are re-calibrated in
   this change. Recommendation: **no** — the TCDB threshold was calibrated
   against the pre-pruning node set and is flagged in `CLAUDE.md` as worth
   revisiting, but mixing a threshold change into a rename would make the
   `post-import-validate` diff unreadable. Separate change.
4. **InterPro has no `is_uninformative` coverage** — the one ontology missing
   from the informativeness filter, which is why `scripts/post-import.cypher`
   (~line 664) excludes `interpro` from `informative_annotation_types` and from
   the `annotation_quality` 8-bucket count. Closing it means an `interpro:`
   section in `uninformative_terms.yaml`, plausibly a `name_patterns` rule over
   "Domain of unknown function" / "Uncharacterised protein family" entries,
   mirroring the existing KEGG pattern. **Deliberately out of scope here:** it
   shifts `annotation_quality`, a released property the MCP reads and gates
   defaults on, and mixing a semantic shift into a rename would make the
   `post-import-validate` diff unreadable — and that diff is the only evidence
   the renames changed nothing else. Separate change, own baseline, own
   `### Breaking` bullet.
