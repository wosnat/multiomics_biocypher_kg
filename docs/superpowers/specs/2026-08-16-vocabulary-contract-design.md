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
| 6 | Four different encodings for a two-state fact — meaningful pair, `"true"`/`"false"` string, native `bool`, sentinel-or-absent — and the native `bool` arm is silently broken (§9.1) | graph-wide |

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
TcdbFamily.is_promiscuous    -> substrate_breadth: multi_substrate | typical
  (many distinct substrates => inferring a member gene's cargo is weak)
InterproEntry.is_promiscuous -> gene_breadth:      ubiquitous | typical
  (matched by a large share of genes => unstratified ORA is invalid)
```

These are **breadth tiers, not booleans**, per R5. The binary framing was
dishonest: `is_multi_gene` would fire at `gene_count >= 1000`, so its false case
means "not ubiquitous", not "one gene" — there is no truthful name for it. A
tier vocabulary states what each side actually is, and leaves room for a middle
band if the thresholds are ever recalibrated (§10.3).

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

`gene_breadth` is **measured frequency**. Its `ubiquitous` members are the opposite of
content-free: IPR027417 (P-loop NTPase superfamily) is highly informative about
fold and function, and is flagged only because it appears in 6,909 genes and so
cannot discriminate in an ORA. The two are independent — informative-but-
ubiquitous and uninformative-but-rare are both real quadrants.

Naming the frequency flag `is_informative` would make it read as the negation of
a released flag that means something else, and would assert that P-loop NTPase
carries no functional content. `gene_breadth` names the axis that is actually
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
lost: the contract publishes `signal_count`, so `round(0.6 × 5) = 3` recovers
it. The published instruction must say **`round`**, not "multiply" — on the
3-signal scale `0.333 × 3 = 0.999`, and a consumer that truncates recovers 0
instead of 1, silently converting the weakest positive evidence into none
(KG-IPT-012). What normalizing buys is the failure mode of a consumer who never reads
the contract — under two scales they compare a Pfam `3` to a TCDB `3` and are
wrong by a factor; normalized, they degrade to "roughly comparable fraction".
Neither scale is calibrated across edge types and the contract says so
explicitly.

`annotation_quality` stays an integer 0–3: it is released, the MCP reads it, and
it is an ordinal *state* (`no_evidence` → `informative_multi`), not a signal
fraction.

### R5 — no native `bool`; a two-state fact is a meaningful categorical string

A boolean-valued property is stored as a string naming **both** states in terms
that read correctly with the column name stripped off — the existing
`OrthologGroup.has_cross_genus_members: cross_genus | single_genus` is the
precedent to generalize. `"true"` / `"false"` is deprecated for new properties;
native `bool` is forbidden outright.

Two independent reasons, and the rule needs only the first:

1. **Adapter-emitted `bool` is broken.** BioCypher does not round-trip it — the
   documented reason `substrate_depth` and `rankable` are already strings. §9.1
   is the first live casualty.
2. **A bare `true` is unreadable in a result row.** These values are consumed by
   an LLM at query time, one row at a time, often with the property name far from
   the value. `multi_substrate` survives that; `true` does not.

**Sentinel-or-absent stays legal** for rare-exception flags — `is_uninformative`,
`level_is_best_effort` — where absence *is* the meaning and the flagged set is a
small minority. Converting those would mean writing a value onto every ontology
term in order to say nothing.

Enforcement: the contract's `value_type` admits `string`, `string_array`,
`float`, `int` and `bool_string`, but **not** `bool`. A property declared `bool`
in `schema_config.yaml` fails the vocabulary test.

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
| `TcdbFamily` | `is_promiscuous` | `bool` | `substrate_breadth: multi_substrate \| typical` | **yes** |
| `InterproEntry` | `is_promiscuous` | `bool` | `gene_breadth: ubiquitous \| typical` | no |
| `Gene_has_tcdb_family` | `agrees_across_sources` | `bool` | `source_agreement: both_sources \| single_source` | no |
| `Gene_has_tcdb_family` | `pfam_corroborated` | `bool` | `pfam_support: corroborated \| uncorroborated` | no |
| `Gene_has_tcdb_family` | `go_corroborated` | `bool` | `go_support: corroborated \| uncorroborated` | no |
| gene→ontology edges | `evidence_score` | int 0–3 | float 0–1 | no |
| `Gene_has_tcdb_family` | `tcdb_evidence_score` | int 0–5 | `evidence_score`, float 0–1 | no |
| `Gene` | `tcdb_best_evidence_score` | int 0–5 | `tcdb_evidence_score_max`, float 0–1 | no |
| `Tcdb_family_transports_metabolite` | `substrate_depth` | `deepest` / `ancestor` | `most_specific` / `inherited` | no |
| `Gene_has_interpro_entry` | `libraries` | `PFAM`, `SUPERFAMILY`, … | `pfam`, `superfamily`, … (13) | no |
| `Interpro_entry_related_to_{ec_number,cazy_family}` | `ambiguous` | `bool` (broken — always false) | `xref_specificity: one_of_several \| sole_xref` | no |

Two rows touch released properties (`substrate_breadth`, and `is_promiscuous`'s
`TcdbFamily` arm is the same row) — and neither has a consumer: the explorer's
audit and ours independently found zero references. Everything else is
unreleased.

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
                                         //   | bool_string
                                         // 'bool' is NOT admissible — see R5
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
                  corroborating signals that fired; round(score * signal_count)
                  recovers the raw count. NOT calibrated against other edge types
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
`sources` / `evidence` / `evidence_score` — **one node per edge type for
`sources` as well as `evidence`**, since the domains differ (`ncbi` appears only
on the GO edges; `Gene_has_cazy_family` admits only `eggnog` + `interproscan`),
and collapsing them would offer `ncbi` as a CAZy filter that can never match
(KG-IPT-011); `Gene_has_interpro_entry.libraries` (13 values, closed,
`string_array`); `Gene_has_tcdb_family` `sources` / `tier` / `evidence_score`;
`Gene.tcdb_evidence_score_max`;
`InterproEntry.interpro_type` / `level_kind` (empty) / `gene_breadth`;
`TcdbFamily.substrate_breadth`; `Gene_has_tcdb_family.{source_agreement,
pfam_support, go_support}`; `Tcdb_family_transports_metabolite.substrate_depth`;
`Gene.transport_substrate_resolution`; Layer-A `xref_specificity` / `source_db`.

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

**CHANGELOG.** One `### Breaking` bullet (`TcdbFamily.is_promiscuous` →
`substrate_breadth`; `evidence_score` integer → normalized float), one
`### Added` bullet for the contract and the R5 boolean rule, plus the KG-IPT-007
prose correction to the existing biller 2016 entry.

**Reply to the explorer.** A short doc answering all seven asks, pointing at §7
and at the rename table in §4, and flagging §7.2 as a correction that unblocks
work they were about to defer.

---

## 9. Follow-up review (explorer, 2026-08-16)

The explorer approved the design and raised four items against it
(`docs/kg-specs/2026-08-16-interpro-tcdb-followup-asks.md`, KG-IPT-009…012).
All four are accepted and folded in above. One needed its diagnosis corrected.

### 9.1 KG-IPT-009 — `ambiguous` is uniformly false (P1, real defect)

The observation is correct: `ambiguous = false` on all 6,854 EC and all 122 CAZy
router edges, including the 3,863 EC edges from non-FAMILY entries that the
documented rule says must be `true`.

**The proposed fix — "fix the `ambiguous` computation" — is wrong.** The
computation is already correct:

```python
# multiomics_kg/adapters/interpro_adapter.py:388
amb = len(ecs) > 1 or etype != "FAMILY"
```

The defect is the **property type**. `config/schema_config.yaml` declares
`ambiguous: bool`, and this KG has a documented BioCypher boolean problem — it
is precisely why `substrate_depth` is a categorical string and why `rankable` /
`has_p_value` are `"true"` / `"false"` strings. Confirmed by blast radius: the
schema has five `bool` properties, and the split is exactly diagnostic —

| Property | Set by | State |
|---|---|---|
| `agrees_across_sources`, `pfam_corroborated`, `go_corroborated` | post-import Cypher | fine |
| `ambiguous` (×2 router edge types) | **adapter** | broken |

Adapter-emitted booleans are broken; post-import booleans are not. That is the
rule R5 now states outright: native `bool` is forbidden, and `value_type` does
not admit it. This defect is what promoted that from a documented quirk to a
rule.

**Fix:** `ambiguous` becomes `xref_specificity: one_of_several | sole_xref` under
R5 — a meaningful pair rather than a stringified boolean. No logic change; only
the property type, name and the two output literals.

**One R1 interaction to catch in implementation:** the guard compares
`etype != "FAMILY"`. R1 lowercases `interpro_type`, so this must become
`"family"` in the same change, or the flag inverts to true everywhere — failing
loudly rather than silently, but still wrong.

**Consequence for the explorer's entry criteria.** Three of their step-2 checks
read properties as native booleans and must be rewritten — R5 removes native
`bool` from the graph entirely, including the post-import ones:

```cypher
r.ambiguous                  ->  r.xref_specificity = 'one_of_several'
t.is_multi_substrate         ->  t.substrate_breadth = 'multi_substrate'
n.is_multi_gene              ->  n.gene_breadth = 'ubiquitous'
```

Expected counts are unchanged: ≥ 3,863 `one_of_several` on the EC router, 13
`multi_substrate`, 22 `ubiquitous`.

### 9.2 KG-IPT-010 — `libraries` violates R1

Accepted. 13 UPPERCASE values, unreleased, the identical problem R1 fixes for
`interpro_type`. Lowercased and added to the v1 seed as `string_array`,
`closed: true`. Their argument for seeding it now rather than later is sound:
`libraries` is the member-DB granularity of the signature-vs-inferred
distinction, which is how a user actually asks for method-independent
corroboration of an eggNOG-transferred annotation.

### 9.3 KG-IPT-011 — `sources` must be per edge type

Accepted; §5.3 amended. They are right that §5.2's argument was made for
`evidence` and then not applied to `sources`, whose domains are equally
non-uniform. The node-id scheme already handles it — this was only a seed-listing
error.

### 9.4 KG-IPT-012 — integer recovery must round

Accepted; §3 R4 and the published `description` both now say
`round(score × signal_count)`.

### 9.5 Accepted minor items

- `source_db` is single-valued (`interpro.xml`, 6,976 edges). Declared as such;
  not a harvest error.
- `Gene.tcdb_evidence_score_max` sentinel guidance in
  `tcdb-two-source-upgrade.md` §2 becomes `coalesce(..., -1.0)` — under
  normalization `0.0` is a legitimate value and the integer `-1` no longer types.

### 9.6 Sequencing — agreed

1. KG lands the §4 renames + the KG-IPT-009 fix, and redeploys.
2. Explorer implements W1 + W4 against the renamed graph.
3. Coordinated KG + MCP release.

`mcp_min_version` protects old-MCP-against-new-KG but not the reverse, which is
the window the explorer would otherwise be developing in. Their entry-criteria
queries are the acceptance check for step 1; the structural counts they list
(`InterproEntry` 12,999 · `Gene_has_interpro_entry` 397,342 · `TcdbFamily` 1,515
· `Gene_has_tcdb_family` 53,763 · `Tcdb_family_transports_metabolite` 11,263 ·
GO inferred 45,226) must be unchanged by the rename pass, and are covered by the
`post-import-validate` diff in §8.

---

## 10. Open items

1. `create_knowledge_graph.py --test` wall-clock — measured during
   implementation; if minutes, the fast gate moves behind its own marker.
2. Exact value sets for the §5.3 second group, harvested and reviewed before
   landing.
3. Whether the `gene_breadth` / `substrate_breadth` thresholds are re-calibrated
   in this change. Recommendation: **no** — the TCDB threshold was calibrated
   against the pre-pruning node set and is flagged in `CLAUDE.md` as worth
   revisiting, but mixing a threshold change into a rename would make the
   `post-import-validate` diff unreadable. Separate change.
5. **Released `"true"` / `"false"` properties predate R5** — `is_time_course`,
   `reports_fold_change`, `rankable` (×2), `has_p_value`, `significant`, and the
   `DerivedMetric` flag `value`. All are MCP-read, so converting them to
   meaningful pairs is breaking and belongs in its own change with its own
   baseline. R5 governs new properties; these are grandfathered until then.
   → `plans/backlog.md`

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
