# Controlled-vocabulary contract + cross-ontology vocabulary alignment

**Date:** 2026-08-16
**Status:** Design — rev 5, not yet implemented. Rev 1 and rev 2 were reviewed
by the explorer; rev 3 folded in their round-2 review; **revs 4–5 delete four
properties rather than rename them — see §0.**
**Driver:** `multiomics_explorer/docs/kg-specs/2026-08-16-interpro-tcdb-asks.md`
(KG-IPT-001 … 008) + `…-interpro-tcdb-followup-asks.md` (KG-IPT-009 … 013)
**Verified against:** KG `0.0.0-dev`, `built_at 2026-08-13T12:19:46.858Z`; last
release tag `kg-0.1.0-alpha.6`

---

## 0. What changed, by revision

### Rev 5 — the derivability audit, applied to everything

Rev 4 deleted one property by asking *what does this tell a consumer that the
graph does not?* Rev 5 applies that test to every flag this integration touches.
Three more fail it.

| Property | Derivable? | Read? | Verdict |
|---|---|---|---|
| Layer-A `source_db` | **a hardcoded literal** — `"interpro.xml"` at `interpro_adapter.py:399` / `:417`, single-valued by construction | no | **deleted** |
| `TcdbFamily.is_promiscuous` | `metabolite_count >= 50 AND level >= 2`, both stored | one internal read | **deleted**, threshold inlined |
| `InterproEntry.is_promiscuous` | `gene_count >= 1000`, stored | no | **deleted** |
| `agrees_across_sources`, `pfam_corroborated`, `go_corroborated` | only via hierarchical traversal / 2-hop bridge checks | yes | keep |
| `substrate_depth` | only via a per-(node, substrate) hierarchy walk | yes | keep |
| `transport_substrate_resolution` | multi-hop chain | yes | keep |
| `tcdb_evidence_score` | from `sources` / `tier` + the three booleans | yes | keep — cheap advisory sort key |

The line that separates them: **materialize traversals, not predicates.** The
keepers stand in for work that is expensive or awkward in Cypher; the deletions
restate an adjacent node's property, a stored count, or a constant.

Consequence: the Layer-A router edges now carry **no properties at all**, which
is correct for a pure router — the edge type is the entire fact. R3 is rewritten
from a naming rule into this modeling rule.

**Two more of the explorer's §6 entry criteria disappear** (`is_multi_substrate`
→ 13, `is_multi_gene` → 22). Replacements, if those exact sets are wanted:

```cypher
MATCH (t:TcdbFamily)    WHERE t.level >= 2 AND t.metabolite_count >= 50  // 13
MATCH (n:InterproEntry) WHERE n.gene_count >= 1000                       // 22
```

### Rev 4 — the Layer-A flag is deleted, not renamed

Three review rounds went into naming a property that turns out to carry no
information. Deleting it resolves **KG-IPT-009 and KG-IPT-013 together**, and
supersedes both the rev-2 rename and the rev-3 split.

| # | Change | Why |
|---|---|---|
| 8 | **`ambiguous` is removed from both Layer-A router edge types.** No `xref_specificity`, no `xref_multiplicity` | Both arms of `len(ecs) > 1 or etype != "FAMILY"` are recoverable from the graph, and the part that is not is arguably wrong — see §9.8 |

A consumer wanting the old flag writes it from data already present:

```cypher
MATCH (n:InterproEntry)-[r:Interpro_entry_related_to_ec_number]->()
WITH n, count(r) AS k
WHERE k > 1 OR n.interpro_type <> 'family'
```

### Rev 3 — round-2 review folded in

The explorer endorsed R3, R4 and R5 with no further comment, and independently
confirmed R5's diagnostic claim. One item is accepted against rev 2:

| # | Change | Why |
|---|---|---|
| 6 | **`xref_specificity` is withdrawn.** The Layer-A edge carries `xref_multiplicity: one_of_several \| sole_xref` — the multiplicity arm **only** | KG-IPT-013. `ambiguous` fuses two orthogonal facts, and my rev-2 name stated only one of them: 1,922 EC router edges (39.5% of flagged) have exactly one xref and are flagged solely for being non-FAMILY. `one_of_several` would contradict their own data — worse than the uninformative flag it replaced |
| 7 | R5's evidence claim is tempered, and its blast-radius table corrected to 7 `bool` pairs across 5 names | The two `is_promiscuous` pairs were omitted from rev 2's table. They are post-import and *work*, which strengthens the case — but the causal claim still rests on one adapter-emitted property, and rev 2 overstated the evidence base |

Rev 3 also corrects an expected count in the explorer's own §6: the flagged
total under `ambiguous` semantics is **4,865**, not 3,863 (3,863 non-FAMILY +
1,002 multi-xref FAMILY). Independently reconciled — the four arms sum to 6,854.

### Rev 2 — changes that postdated the rev-1 approval

Rev 1 is the version the explorer approved. Everything below was new at rev 2,
and two items **supersede things that were approved**. Nothing here changes the
sequencing agreed in §9.6.

| # | Change | Why | Action for the explorer |
|---|---|---|---|
| 1 | **New house rule R5** — no native `bool`; a two-state fact is a categorical string naming both states meaningfully | KG-IPT-009 proved adapter-emitted `bool` is broken, and a bare `true` is unreadable in a result row an LLM reads one line at a time | Review R5 (§3) |
| 2 | **R3 revised: breadth *tiers*, not flags.** `is_multi_substrate` / `is_multi_gene` **no longer exist** — they are `substrate_breadth: multi_substrate \| typical` and `gene_breadth: ubiquitous \| typical` | R5 exposed the binary as dishonest: `is_multi_gene` fires at `gene_count >= 1000`, so its false case means "not ubiquitous", not "one gene" — the negative had no truthful name | **Supersedes an approved name.** Update §6 entry criteria |
| 3 | **KG-IPT-009 diagnosis corrected, and the fix is larger than asked.** The `ambiguous` computation is already correct; the `bool` *type* is the defect. ~~The property becomes `xref_specificity`~~ — **superseded at rev 3, see #6 above** | See §9.1 — the blast radius across all 7 `bool` (entity, property) pairs is exactly diagnostic: post-import bools work, the one adapter-emitted bool does not | Note the fix is a rename, not just a retype |
| 4 | Three of your §6 step-2 entry-criteria queries read native booleans and no longer parse | R5 removes native `bool` from the graph entirely | **Rewrite required** — mapping in §9.1; expected counts unchanged (≥3,863 · 13 · 22) |
| 5 | Seven released `"true"` / `"false"` properties are grandfathered, not converted | All MCP-read; converting is breaking and would obscure the rename pass's validate diff | None now — backlogged (§10.5) |

Accepted from your followup with no further comment: KG-IPT-010 (`libraries`
lowercased + seeded), KG-IPT-011 (`sources` seeded per edge type),
KG-IPT-012 (`round(score × signal_count)`), and both §5 minor items.

---

## 1. Problem

The InterPro + TCDB two-source integration introduced five new controlled
vocabularies. They were each designed in isolation, against their own source
database's conventions, and the result is not uniform:

| # | Inconsistency | Where |
|---|---|---|
| 1 | `interpro_type` values are `UPPERCASE`; every other categorical in the KG is lowercase `snake_case` | `InterproEntry` |
| 2 | One provider, three spellings: edge `sources` says `interpro`, `DataSource` / `Gene.contributing_sources` say `interproscan`; the TCDB edge says `diamond`, `DataSource` says `tcdb_diamond` | gene→ontology edges, `Gene_has_tcdb_family` |
| 3 | `is_promiscuous` means *many substrates* on `TcdbFamily` but *many genes* on `InterproEntry` — one name, two axes. *Resolved by deletion at rev 5, not by renaming* | `TcdbFamily`, `InterproEntry` |
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
`Gene.transport_substrate_resolution` · Layer-A `{ambiguous, source_db}`
(both deleted at revs 4–5 rather than aligned).

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

### R1 — lowercase `snake_case` for values the KG mints; external terms verbatim

**Scope, rescoped 2026-08-17.** R1 governs values the KG **mints**. A value that is an
**external database's own controlled term** is preserved **verbatim**, because mutating
it makes the graph confusing to compare against the source. Casing uniformity is worth
less than traceability to InterPro / InterProScan / NCBI.

| Preserved verbatim (external) | Lowercased (KG-minted) |
|---|---|
| `InterproEntry.interpro_type` — `FAMILY`, `DOMAIN`, … | `evidence`, `sources`, `substrate_depth` |
| `Gene_has_interpro_entry.libraries` / `evalue_library` — `PFAM`, `GENE3D`, … | `transport_substrate_resolution`, `level_kind` |
| `NcbifamFamily.family_type` — incl. `PfamEq`, `PfamAutoEq` | `annotation_state`, `expression_status`, `metric_bucket` |
| | `MeropsFamily.catalytic_type`, `call_class`, `best_hit_kind` |

`catalytic_type` is KG-minted despite describing MEROPS data: the source ships single
letters (`S`, `C`, `M`, …) and the KG chooses the words (`serine`, `asparagine_lyase`).

The contract records the reason on each external entry's `description`, so a later reader
does not "fix" the casing back.

### R1b — namespace only when values collide across labels

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

### R2 — every `sources` value corresponds to a `DataSource` node

Provenance becomes joinable, and the three-spellings problem disappears. The
`sources` value itself is bare (e.g. `'tcdb_diamond'`); the matching
`DataSource` node's `id` property carries the `data_source:` prefix
(`'data_source:tcdb_diamond'`), so the join is `d.id = 'data_source:' + s`,
not `d.id = s`:

```
gene→ontology edges:      'interpro'  -> 'interproscan'
Gene_has_tcdb_family:     'diamond'   -> 'tcdb_diamond'
```

Enforced by a kg-validity test: no `sources` value may lack a matching
`DataSource` under that prefixed join.

### R3 — do not materialize a threshold over a stored count

```
TcdbFamily.is_promiscuous     -> DELETED   (was metabolite_count >= 50 AND level >= 2)
InterproEntry.is_promiscuous  -> DELETED   (was gene_count >= 1000)
```

Both were booleans computed from a count the node already publishes. A consumer
can apply any cutoff to `metabolite_count` / `gene_count` directly, and sees the
magnitude rather than a bit. Materializing the KG's cutoff instead hides that a
judgement was made, and goes stale silently — both thresholds have been
recalibrated once already, and `CLAUDE.md` flags the TCDB one as due again.

The KG's own cutoffs are **documented, not stored**: `metabolite_count >= 50` at
`level >= 2` flags 13 TCDB families; `gene_count >= 1000` flags 22 InterPro
entries. Anyone wanting exactly those sets writes the predicate.

**The one internal consumer is inlined.** `post-import.cypher:1131` reads
`t.is_promiscuous` to compute `Gene.transport_substrate_resolution`; the
threshold moves into that expression. `transport_substrate_resolution` stays —
unlike a threshold restatement, it materializes a multi-hop traversal (gene →
deepest TC attachments → substrate breadth) that is genuinely awkward at query
time. That is the line R3 draws: **materialize traversals, not predicates.**

This resolves inconsistency #3 by deletion rather than by naming, the same shape
as #6. It does not weaken KG-IPT-006 — the ask was that breadth never becomes a
default filter, and a property that does not exist cannot be one.

**Why this is not the `is_informative` axis.** Worth recording, because the two
look adjacent and are not. `is_uninformative`
(`config/uninformative_terms.yaml`) is **curated content-emptiness** — "conveys
no class signal at all": the three GO roots, `cog.category:S`,
`cyanorak.role:R.4`. Its stated principle is that a term naming the broad class
stays informative even when the sub-class is unknown, which is why DUF/UPF Pfams
and "Transport > Unknown substrate" are deliberately not flagged. It is curated,
not derivable, and therefore *does* earn storage — the exact opposite of the
breadth flags.

Breadth was **measured frequency**. IPR027417 (P-loop NTPase superfamily) is
highly informative about fold and function and would have been flagged only
because it appears in 6,909 genes. The two axes are independent —
informative-but-ubiquitous and uninformative-but-rare are both real quadrants —
which is why the frequency axis was never a candidate for the `is_informative`
name, and is now not a stored property at all.

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
   is the first live casualty. *Stated precisely:* the causal claim rests on a
   single property, because `ambiguous` is the only adapter-emitted `bool` in
   the schema. The prior is strong and the fix is correct either way, since R5
   removes the category rather than repairing it — but the evidence base is one
   property, not a survey.
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
| `TcdbFamily` | `is_promiscuous` | `bool` | **removed** (§3 R3) | **yes** |
| `InterproEntry` | `is_promiscuous` | `bool` | **removed** (§3 R3) | no |
| `Gene_has_tcdb_family` | `agrees_across_sources` | `bool` | `source_agreement: both_sources \| single_source` | no |
| `Gene_has_tcdb_family` | `pfam_corroborated` | `bool` | `pfam_support: corroborated \| uncorroborated` | no |
| `Gene_has_tcdb_family` | `go_corroborated` | `bool` | `go_support: corroborated \| uncorroborated` | no |
| gene→ontology edges | `evidence_score` | int 0–3 | float 0–1 | no |
| `Gene_has_tcdb_family` | `tcdb_evidence_score` | int 0–5 | `evidence_score`, float 0–1 | no |
| `Gene` | `tcdb_best_evidence_score` | int 0–5 | `tcdb_evidence_score_max`, float 0–1 | no |
| `Tcdb_family_transports_metabolite` | `substrate_depth` | `deepest` / `ancestor` | `most_specific` / `inherited` | no |
| `Gene_has_interpro_entry` | `libraries` | `PFAM`, `SUPERFAMILY`, … | `pfam`, `superfamily`, … (13) | no |
| `Interpro_entry_related_to_{ec_number,cazy_family}` | `ambiguous` | `bool` (broken — always false) | **removed** (§9.8) | no |
| `Interpro_entry_related_to_{ec_number,cazy_family}` | `source_db` | `"interpro.xml"` constant | **removed** (§0 rev 5) | no |

One row touches a released property — `TcdbFamily.is_promiscuous`, deleted —
and it has no consumer: the explorer's audit and ours independently found zero
references. Everything else is unreleased.

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
`InterproEntry.interpro_type` / `level_kind` (empty);
`Gene_has_tcdb_family.{source_agreement, pfam_support, go_support}`;
`Tcdb_family_transports_metabolite.substrate_depth`;
`Gene.transport_substrate_resolution`. The Layer-A router edges seed nothing —
they carry no properties.

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

**As implemented, this is a detection net across test runs, not a build-time
guard on every emitter.** The loader's `VOCAB.check()` helper is wired into
exactly one adapter (`tcdb_adapter.py`, for `substrate_depth`); every other
vocabulary is emitted with no loader involvement, so an adapter drifting from
the YAML is caught by the CSV/live-graph scans below, not blocked at emit
time. Post-import Cypher cannot import Python either, so its literals are
covered by the graph-level test only. Four gates, each catching something the
others structurally cannot:

| Gate | Marker | Runs | Checks |
|---|---|---|---|
| Adapter units — `VOCAB.check()` raises at emit | none | always | undeclared value, in seconds, but only where an adapter calls it (currently `tcdb_adapter.py` only) |
| `--test` build → CSV scan | none | always | observed ⊆ declared, pre-import, no Neo4j |
| Full build → CSV scan | `slow` | opt-in | + declared − observed ⊆ `expected_empty`, but only on `exhaustive: true` entries |
| Live graph | `kg` | every rebuild | observed ⊆ declared always; declared − observed only on `exhaustive: true` entries · post-import-only properties · R2's `DataSource` join · `controlled_vocabularies_hash` |

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
  — including two stale items in `tcdb-two-source-upgrade.md` §2: the
  `tcdb_evidence_score` distribution predates the 2026-08-13 rebuild (doc says
  `0 → 17,422 · 5 → 1,081`; graph has `0 → 17,045 · 1 → 8,599 · 2 → 9,461 ·
  3 → 7,541 · 4 → 9,957 · 5 → 1,160`), and it is re-scaled to floats by R4
  anyway; and the `coalesce(..., -1)` sentinel guidance becomes `-1.0`
- `CLAUDE.md`, `CHANGELOG.md`

**Validation gate.** `scripts/post-import-validate.sh` baseline captured against
the currently deployed graph *before* any change; after the rebuild the diff must
be empty except for the renamed properties and the normalized score columns —
every count and every other value byte-identical. Then
`pytest -m "not slow and not kg"`, Docker rebuild, `pytest -m kg`, and
`/omics-edge-snapshot`.

**CHANGELOG.** One `### Breaking` bullet (`TcdbFamily.is_promiscuous` →
deleted; `evidence_score` integer → normalized float), one
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

| Set by | Property | Live distribution | State |
|---|---|---|---|
| post-import | `agrees_across_sources` | 21,684 true / 32,079 false | fine |
| post-import | `pfam_corroborated` | 23,634 / 30,129 | fine |
| post-import | `go_corroborated` | 26,885 / 26,878 | fine |
| post-import | `TcdbFamily.is_promiscuous` | 13 true | fine |
| post-import | `InterproEntry.is_promiscuous` | 22 true | fine |
| **adapter** | `ambiguous` (×2 router edge types) | **6,976 false, 0 true** | **broken** |

Seven `(entity, property)` pairs across five names. Every post-import one works;
the only adapter-emitted one does not.

Adapter-emitted booleans are broken; post-import booleans are not. That is the
rule R5 now states outright: native `bool` is forbidden, and `value_type` does
not admit it. This defect is what promoted that from a documented quirk to a
rule.

**Fix (final, rev 4 — see §9.8): the property is deleted.** Rev 2 renamed it,
rev 3 split it; rev 4 establishes that neither arm is worth storing.

**One R1 interaction to catch in implementation:** the guard compares
`etype != "FAMILY"`. R1 lowercases `interpro_type`, so this must become
`"family"` in the same change, or the flag inverts to true everywhere — failing
loudly rather than silently, but still wrong.

**Consequence for the explorer's entry criteria.** Three of their step-2 checks
read properties as native booleans and must be rewritten — R5 removes native
`bool` from the graph entirely, including the post-import ones:

```cypher
// no replacement properties — all three derive from what is already stored:
t.is_multi_substrate  ->  t.level >= 2 AND t.metabolite_count >= 50    // 13
n.is_multi_gene       ->  n.gene_count >= 1000                         // 22

// `r.ambiguous` (4,865 on the EC router):
MATCH (n:InterproEntry)-[r:Interpro_entry_related_to_ec_number]->()
WITH n, count(r) AS k
WHERE k > 1 OR n.interpro_type <> 'family'
```

An emptiness check replaces the old single-valuedness assertion: the router
edges must carry no properties at all.

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

- `source_db` is single-valued (`interpro.xml`, 6,976 edges) — correctly
  observed, and the reason is stronger than they knew: it is a hardcoded string
  literal in the adapter, not harvested data. **Deleted at rev 5** rather than
  declared.
- `Gene.tcdb_evidence_score_max` sentinel guidance in
  `tcdb-two-source-upgrade.md` §2 becomes `coalesce(..., -1.0)` — under
  normalization `0.0` is a legitimate value and the integer `-1` no longer types.

### 9.7 Round-2 review (rev 2 → rev 3)

R3, R4 and R5 endorsed with no further comment; the `gene_breadth` /
breadth reframing was called an improvement on the flags approved at rev 1 —
since superseded, as rev 5 deletes both flags outright. One ask raised and accepted:

**KG-IPT-013 — `xref_specificity` names only one of two arms.** Correct, and the
error is provable from the code without touching the graph:

```python
amb = len(ecs) > 1 or etype != "FAMILY"     # interpro_adapter.py:388
```

A DOMAIN entry carrying exactly one EC satisfies the second disjunct. Under my
rev-2 name that edge would read `one_of_several` while being, demonstrably, the
sole xref on its entry — the value contradicted by the row's own data. On the EC
router that is **1,922 of 4,865 flagged edges, 39.5%**. Worse than the defect it
replaced: `ambiguous` was uninformative, this would have been actively false, and
published inside a contract asserting the value set is meaningful.

Their four-arm decomposition reconciles exactly and I verified the arithmetic
independently: multi/FAMILY 1,002 + multi/non-FAMILY 1,941 + sole/non-FAMILY
1,922 + sole/FAMILY 1,989 = 6,854, with the non-FAMILY arms summing to the 3,863
from §9.1's type breakdown.

**Accepted as recommended: split the axes rather than rename the fusion.** The
edge carries multiplicity only; type stays on the source node where it already
lives. Strictly more informative than the fused flag, R5-compliant, and every
value is true of every row carrying it — which is the acceptance criterion they
proposed and the right one.

### 9.8 Rev 4 — the flag is deleted (supersedes 9.7)

KG-IPT-013's decomposition prompted the question the three naming rounds had
skipped: *what does this property tell a consumer that the graph does not?*
Nothing, and the residue is arguably wrong.

**The type arm carries zero information.** The edge's source is an
`InterproEntry` and `interpro_type` is on it. `etype != "FAMILY"` restates a
property of the node the edge starts at.

**The multiplicity arm is the entry's out-degree** on that edge type, which any
consumer can count. There is one delta: the adapter computes `len(ecs) > 1` from
the full reference list and only afterwards prunes edges to ECs that have nodes,
so the stored value is *pre-pruning*. The pruning is real — the reference carries
entries with up to **23** EC xrefs (2,596 of 10,849 entries have more than one)
while the graph tops out at 14.

**But that delta is a defect, not information.** The pruned ECs are the obsolete
and invalid ones — the tokens `normalize_ec` cannot remap and Expasy has no node
for, the same class that motivated the EC dangling-proof fix. If an entry lists
23 ECs and 22 are dead, the survivor is the *only* valid claim, and flagging that
link "one of several" is wrong. Post-pruning out-degree is the better semantics,
and it is free.

So the property is derivable, denormalized (an entry-level fact copied onto every
edge leaving that entry), and where it is not derivable it is less correct. It has
also never been read: Layer-A routers are deferred explorer-side under W3.

**Deleted from both router edge types**, along with `source_db` at rev 5, so
these edges carry no properties at all. KG-IPT-009 (the flag is uniformly false) and KG-IPT-013 (the replacement name is
untruthful) are both resolved — there is no flag.

The guardrail KG-IPT-009 was defending is not weakened. It never lived in this
property, which was `false` on 100% of edges in every build ever deployed. It
lives in the deliberately weak `related_to` verb, in the edge type name, and in
the documented ~31% reverse precision — none of which this change touches.

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
3. Whether the breadth thresholds are re-calibrated in this change. Only one
   survives, inlined into `transport_substrate_resolution` (§3 R3).
   Recommendation: **no** — the TCDB threshold was calibrated
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

---

## 11. Implementation status (2026-08-18)

Branch `docs/vocabulary-contract-spec`. **Implemented and reviewed; NOT yet validated
against a live graph.**

### Done

Tasks 1–11 of `plans/vocabulary_contract.md`, each with a task review; one Critical was
caught and fixed mid-flight (Task 6), and a final whole-branch review on the strongest
model produced 1 Critical + 7 Important + ~8 Minor, all addressed in one fix wave with a
scoped re-review. Fast suite green throughout (2,365 passed). The build-level drift gate
passes with no drift.

The branch was **merged with `origin/main` partway through**, absorbing the InterPro
multi-ontology redesign plus the MEROPS and NCBIfam ontologies, and the contract was
extended to cover them (66 declared vocabularies).

### Changed since this spec was written

- **R1 was rescoped** (§3): external database terms are preserved verbatim. `interpro_type`
  stays `FAMILY`, `libraries` stays `PFAM`, `NcbifamFamily.family_type` keeps `PfamEq`.
  The normalizers were deleted outright.
- **`DerivedMetric.metric_type` is `closed: false`**, not closed. `KNOWN_METRIC_TYPES` is
  documented as a *soft allowlist* — a paper may legitimately introduce a new type — so a
  closed declaration would fail every time a paper is added.
- **`exhaustive: bool` was added** to the contract. `closed: true` means "no value outside
  this set may appear" and is asserted unconditionally; the declared-minus-observed
  *coverage* direction is opt-in, because several entries legitimately declare a
  documented superset of what the graph emits. No entry opts in yet.
- Gate 1 as described in §6 **does not exist**: `VOCAB.check` is wired at one call site,
  not into every emitter. The real coverage is the CSV scan and the live-graph suite.

### Not done — Task 12, the live validation gate

Requires a person present: it runs `docker compose down` plus a ~1 h rebuild, and its
acceptance baseline must be captured from the **currently deployed pre-change graph**
before anything is rebuilt. Nothing on this branch has been built end to end.

`pytest -m kg` will open with roughly **14 red tests**, of which **one** is a real
assertion about new behaviour and the rest are stale expectations this branch invalidated.
Fix these as part of Task 12 rather than treating them as regressions:

| File | What is stale |
|---|---|
| `tests/kg_validity/test_interpro.py:112-126` | requires the deleted `InterproEntry.is_promiscuous` |
| `tests/kg_validity/test_tcdb_cazy.py` | ~10 tests: `tcdb_evidence_score` as int 0-5 (237-259); the three booleans (269-297); `is_promiscuous` (314, 322-354, 388-398); `tcdb_best_evidence_score` (370-381); `'diamond'` (41, 59); `substrate_depth` (419-441) |
| `tests/kg_validity/test_ncbifam.py:111` | asserts `"interpro" in sources` |
| `tests/kg_validity/snapshot_data.json:1443,1458` | `"diamond"` — regenerate the snapshot |

Then run the §8 validation gate: `post-import-validate` baseline before the rebuild, diff
after (only renamed properties and re-scaled score columns may differ), `pytest -m kg`,
and `/omics-edge-snapshot` to confirm zero per-paper expression-edge deltas.

### Sequencing with the explorer

Unchanged from §9.6 — the explorer's W1 + W4 work is gated on this landing and
redeploying. Their step-2 entry criteria need three query rewrites (§0 revs 4–5, §9.1),
and note that under the R1 rescope `interpro_type` is **not** lowercased after all.
