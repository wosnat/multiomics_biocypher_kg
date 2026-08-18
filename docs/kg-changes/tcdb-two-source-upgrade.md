# TCDB two-source upgrade — MCP / explorer contract

**Status:** LANDED 2026-08-07 · extended 2026-08-12 with the substrate-rollup
depth fix (§7) — same unreleased change set, not yet consumed by the explorer,
so both land as ONE upgrade.
**Audience:** the MCP / explorer team. This is an **upgrade delta**, not a full
description — it assumes the existing TCDB surface documented in
[`tcdb-cazy-ontologies.md`](tcdb-cazy-ontologies.md) is already implemented in the API.
**Design:** [`2026-08-06-tcdb-diamond-kg-integration-design.md`](../superpowers/specs/2026-08-06-tcdb-diamond-kg-integration-design.md)

---

## TL;DR — what breaks, what's new

| | |
|---|---|
| ⚠️ **Behaviour change** | `'tcdb' IN g.annotation_types` no longer follows from "has a `Gene_has_tcdb_family` edge" — it is tier-gated. Affects **15,558 genes**. |
| ⚠️ **BREAKING (2026-08-18)** | `TcdbFamily.is_promiscuous` is **deleted** (vocabulary-contract R3) — it only ever restated `level >= 2 AND metabolite_count >= 50`, both already stored. Query the threshold directly. Node set went 12,902 → 1,515, so any hard-coded count is stale regardless. |
| ⚠️ **Scale change** | `Gene_has_tcdb_family` 16,806 → **53,763** edges; genes with a TC call 11,103 → **30,076**. Result-size assumptions and default limits need review. |
| ✅ **New** | Edge provenance (`sources`) + diamond evidence + an advisory `evidence_score` (renamed 2026-08-18 from `tcdb_evidence_score`, now a float `[0,1]` — see §2). |
| ✅ **New** | 4 ontology→ontology bridge edges (TcdbFamily → Pfam / GO). |
| ✅ **New** | `Gene.tcdb_evidence_score_max` (renamed 2026-08-18 from `tcdb_best_evidence_score`, now a float `[0,1]`). |
| ⚠️ **BREAKING** | `metabolite_count` no longer unions catalysis with transport. Gene / Metabolite / OrganismTaxon each gain a separate transport-arm property. Gene p90 was 554 (transport-dominated), catalysis-only p90 is 11. See §7.3. |
| ⚠️ **Semantics change** | `Metabolite.transporter_count` was 0 for 83% of transported metabolites; now non-zero for all 1,462. See §7.2. |
| ✅ **New** | `Tcdb_family_transports_metabolite.substrate_depth` (`'most_specific'`/`'inherited'`, renamed 2026-08-18 from `'deepest'`/`'ancestor'`) — tells a real transporter system from a rollup ancestor. See §7.1. |
| ✅ **New** | `Gene.transport_substrate_resolution` (`'resolved'`/`'family_inferred'`) — is this gene's substrate tag usable? See §7.4. |

---

## 1. `Gene_has_tcdb_family` now has two evidence sources

eggNOG's `KEGG_TC` (ortholog **transfer**) and `/tcdb-diamond` (direct **sequence
similarity**) both emit TCDB TC IDs. They feed **one ontology** and **one edge per
(gene, TC id)** — never parallel edges.

```cypher
(Gene)-[:Gene_has_tcdb_family {
    sources:          ['eggnog','tcdb_diamond'],  // always present, SORTED
    tier:             2,                      // ↓ tcdb_diamond only — SPARSE, absent on eggNOG-only
    confidence_score: 0.85,
    identity:         50.5,
    qcov:             78.0,
    evalue:           1.15e-74,
    consensus_n:      4
}]->(TcdbFamily)
```

| `sources` | Edges |
|---|---|
| `['tcdb_diamond']` | 36,957 |
| `['eggnog']` | 13,165 |
| `['eggnog','tcdb_diamond']` | 3,641 |

**Do not assume `tier` exists.** It is absent on all 13,165 eggNOG-only edges. Use
`coalesce(r.tier, …)` or branch on `'tcdb_diamond' IN r.sources`.

**`sources` is sorted** — the canonical two-source value is `['eggnog','tcdb_diamond']`
(`tcdb_diamond` sorts after `eggnog`). Even so, prefer membership (`'eggnog' IN
r.sources`) over list equality; it survives a future third source.

### Corroboration is hierarchical — this is the easy thing to get wrong

`size(r.sources) = 2` is **not** the corroboration test. The two sources usually agree
at *different depths*: eggNOG names subfamilies (`3.A.5.1`), diamond's tier-3
truncation names the parent family (`3.A.5`). Same biology, different TC ids, so
different edges.

| Cross-source agreement | Edges |
|---|---|
| exact same node (`size(sources)=2`) | 3,641 |
| **hierarchical** (one is the other's ancestor) | 18,043 |
| **any** | **21,684** |

Filtering on `size(sources)=2` silently discards **83%** of real agreement. Use the
precomputed `r.source_agreement = 'both_sources'` (below) rather than recomputing.

---

## 2. Advisory scores — `evidence_score` and the gene rollup

**⚠️ Renamed and re-scaled (2026-08-18, `docs/kg-changes/vocabulary-contract.md` R4).**
The edge property `tcdb_evidence_score` is now **`evidence_score`** — the same name
every other annotation edge type uses — and it is a **float in `[0,1]`**, not an
integer 0–5. The two component booleans are also gone as native `bool` (R5); they are
now the categorical strings `source_agreement`, `pfam_support`, `go_support` (see the
rename table). The distribution figures and `coalesce` guidance below **predate both
changes** and are corrected in place rather than reproduced twice.

**Edge property `evidence_score` (float `[0,1]`)** — five independent signals,
`+1` each, then divided by 5 and rounded to 3 decimals; `round(score × 5)` recovers the
raw fired-signal count:

| Component | Stored as |
|---|---|
| eggNOG called it | `'eggnog' IN sources` |
| sources agree (hierarchically) | `source_agreement = 'both_sources'` |
| strong direct sequence evidence | `tier <= 2` |
| a gene Pfam is curated into this family | `pfam_support = 'corroborated'` |
| a gene GO term is curated onto this family | `go_support = 'corroborated'` |

The 2026-08-13 distribution figures this doc previously quoted for the 0–5 integer scale
are stale relative to both the InterPro two-layer merge (which shifted the underlying
counts) and this rename (which changed the scale entirely) — not re-verified here; do
not treat any previously-published distribution as current. Re-measure against a live
build if a distribution is needed.

**Gene property `tcdb_evidence_score_max` (float `[0,1]`, renamed from
`tcdb_best_evidence_score`)** = `max()` over the gene's `evidence_score` values.
Answers *"is this gene a transporter at all"*, where the edge score answers *"is THIS
assignment right"*.

> **SPARSE.** Set only on genes with ≥1 TCDB edge. A gene with **no** transporter
> evidence is deliberately distinct from one with **weak** evidence — do not
> `coalesce(..., 0)`, that collapses the distinction. Use `coalesce(..., -1.0)` — not
> `-1` — if a total order is required: under the float scale, `0.0` is now a
> legitimate value, and an integer sentinel no longer even types against the property.

### Both are ADVISORY — please keep them that way

The components are stored beside the total precisely so the API can show *why* and
re-weight without re-deriving. Two asks:

- **Don't hard-filter on the score by default.** Rank with it, expose it, let the user
  threshold. The predecessor of this score (`filter_action`) was deleted after it was
  measured to exclude ~13K candidates on incoherent grounds; this one is additive and
  sibling-independent specifically so that cannot recur.
- **If you do offer a cutoff, surface the components alongside it** so a filtered-out
  call is explainable.

---

## 3. ⚠️ `annotation_types` is tier-gated — behaviour change

`'tcdb'` enters `Gene.annotation_types`, `informative_annotation_types` and the
`annotation_quality` bucket **only** when the gene has an eggNOG-sourced **or**
`tier <= 2` edge.

| | Genes |
|---|---|
| with ≥1 `Gene_has_tcdb_family` edge | 30,076 |
| counted in `annotation_types` | **14,518** |
| **carry the edge but are NOT counted** | **15,558** |

Rationale: tier 3 is conservative-by-design remote homology (median identity 34%,
median e-value 2.3e-30) — good enough to be *findable*, not good enough to inflate a
quality signal calibrated against curated sources.

**Action:** any explorer query that treats "has a TCDB edge" and "`'tcdb'` in
`annotation_types`" as equivalent is now wrong for 15,558 genes. Pick deliberately:
edge presence for recall, `annotation_types` for quality.

`Gene.tcdb_family_count` and `Gene.metabolite_count` are **not** gated — they count all
edges, being routing signals rather than quality signals.

---

## 4. ⚠️ `is_promiscuous` — DELETED (2026-08-18)

**This property no longer exists on `TcdbFamily`.** Per `docs/kg-changes/vocabulary-contract.md`
house rule R3, a property that only ever restates a threshold over a count the node
already publishes is deleted rather than maintained — `is_promiscuous` was exactly
`level >= 2 AND metabolite_count >= 50`, and both `level` and `metabolite_count` are
already stored. A consumer wanting the same 13 flagged nodes queries the predicate
directly:

```cypher
MATCH (t:TcdbFamily) WHERE t.level >= 2 AND t.metabolite_count >= 50 RETURN t
```

The rest of this section is kept as historical record of the level-gating rationale
(the threshold itself is unchanged by the deletion) — **the property name in the code
blocks below no longer exists on the graph.**

**Meaning is unchanged and is what it always was: this family transports MANY DISTINCT
SUBSTRATES**, so inferring what a member gene moves from family membership is weak.

What changed is that it is now **level-gated** (`level >= 2`, i.e. `tc_family` and
deeper) and it means only that.

```cypher
is_promiscuous = level >= 2 AND metabolite_count >= 50
```

- **Why level-gated:** substrate counts scale mechanically with level, since the step-6
  rollup puts every descendant's substrates on each ancestor (median `metabolite_count`
  153 at `tc_class` vs 1 at `tc_family`). The old rule fired on 5 of 7 `tc_class` and 7
  of 34 `tc_subclass` nodes — vacuous, because "Channels and Pores transports many
  things" is what a class *is*.
- **The old `member_count >= 100` arm is gone.** It was dead where it mattered: max
  `member_count` at `tc_family` is 55, so it could only ever fire at `tc_subclass`.

> **Note for anyone who saw the 2026-08-07 morning build:** a `gene_count >= 500` arm
> was briefly added and reverted the same day. It answered a different question ("is
> this a large bucket of genes?") and flagged substrate-*poor* families — e.g. `9.B.34`
> (KPSH) with **zero** substrates. Overloading one boolean with two axes destroyed the
> term. **For family size, filter `t.gene_count` directly** — it is already on the node,
> so no second boolean was added.

**Nothing you already rely on changes.** At `tc_family` / `tc_subfamily` the old and new
rules agree exactly (7 and 6 flags respectively) — the `member_count` arm never fired
below `tc_subclass`, so those levels were always driven by `metabolite_count >= 50`.
The change removes only the 12 vacuous class/subclass flags: **25 → 13**.

All 13 are textbook multi-substrate transporters: ABC `3.A.1` (554 substrates) ·
MFS `2.A.1` (476) · DMT `2.A.7` · RND `2.A.6` · MOP flippase `2.A.66` · APC `2.A.3` ·
P-type ATPase `3.A.3`.

⚠️ [`metabolomics-extension.md`](metabolomics-extension.md) previously documented the
old rule (`metabolite_count >= 50 OR member_count >= 100`, "~30 of 12,883 families").
That row is now corrected; this doc is authoritative.

---

## 5. New ontology→ontology bridge edges

From TCDB's published Pfam→TC and GO→TC maps. **No gene-touching edges were added.**

| Edge | Source → Target | Count |
|---|---|---|
| `Tcdb_family_has_pfam_domain` | TcdbFamily → Pfam | 1,853 |
| `Tcdb_family_involved_in_biological_process` | TcdbFamily → BiologicalProcess | 2,575 |
| `Tcdb_family_enables_molecular_function` | TcdbFamily → MolecularFunction | 2,263 |
| `Tcdb_family_located_in_cellular_component` | TcdbFamily → CellularComponent | 1,708 |

Each carries `curated_tcids: str[]` — the original published 5-part TCIDs, so roll-up
onto surviving ancestors loses no precision.

### 🚫 Direction is load-bearing — do not traverse these backwards

They assert *"TCDB's curated reference proteins for this family carry this domain / GO
term"* — a statement about the **composition of the family**.

| Reading | Precision |
|---|---|
| **Forward** — family known → does the xref agree? | **85% corroborate**, 2% contradict |
| **Reverse** — gene has the domain → is it a transporter? | **~31%** |

TCDB curates whole transport *systems*, so accessory proteins are included:
response-regulator receiver, histidine kinase, TPR and GGDEF domains all map to TC
families without being transporters.

- ✅ Corroborate a `Gene_has_tcdb_family` call you already have (this is what
  `pfam_support` / `go_support` do).
- ✅ "What domains/functions characterise this family?"
- ✅ `gene → TcdbFamily → GO` — 2 hops, first is direct evidence. Of the 4,997 TCDB-annotated
  genes with no direct GO, **2,914** gain reachable GO terms this way.
- ❌ **Never** assign transporter identity from `Pfam → TcdbFamily`. A
  `gene → InterproEntry → Pfam → TcdbFamily` traversal is a recall-biased router
  (like `Publication_discusses_gene`), not an annotation.

---

## 6. Node-set and count changes

`TcdbFamily`: **12,902 → 1,515**. Pruning is now ancestor-only — the old bidirectional
walk left 94.5% of nodes with `gene_count = 0`. Zero metabolite reachability was lost.

| Level | Nodes |
|---|---|
| `tc_class` | 7 |
| `tc_subclass` | 34 |
| `tc_family` | 592 |
| `tc_subfamily` | 596 |
| `tc_specificity` | 286 |

`Tcdb_family_transports_metabolite`: 22,483 → 11,263.

⚠️ **`level_kind = 'tc_specificity'` shifted meaning.** It now selects the specificity
nodes genes actually annotate (286), not every specificity node under a reachable
family. Queries using that filter to recover "leaf-only" substrate semantics still
work, but now count real transporter systems present in our organisms.

⚠️ **Do NOT use `level_kind = 'tc_specificity'` as a substrate filter any more.**
Only **466 of 11,263** substrate edges now sit on specificity nodes, because genes
mostly annotate at family/subfamily depth and the ancestor-only prune keeps no
specificity node beneath them. Use `substrate_depth` (§7) instead.

---

## 7. Substrate rollup depth — new edge property + rebuilt counts

The step-6 rollup materialises **every descendant's substrates onto every
ancestor**, computed over the FULL TCDB hierarchy before pruning. That is
deliberate (it is why an ancestor stays substrate-reachable after its leaves are
pruned away), but it left no way to tell *"this node is the transporter system"*
from *"this node is an ancestor of one"*. Three scalars silently inherited the
ambiguity. All are fixed together here.

### 7.1 New: `Tcdb_family_transports_metabolite.substrate_depth`

```cypher
(TcdbFamily)-[:Tcdb_family_transports_metabolite {
    substrate_depth: 'most_specific'   // no kept CHILD of this node carries this substrate
                                        // → the most specific surviving system for it
  //substrate_depth: 'inherited'     // carries it only via the subtree rollup
}]->(Metabolite)
```

Property renamed 2026-08-18 (`docs/kg-changes/vocabulary-contract.md` R5, `'deepest'`/
`'ancestor'` → `'most_specific'`/`'inherited'`). Categorical **string**, not boolean
(BioCypher mishandles bool properties — the KG uses string vocabularies throughout).
Counts not re-verified against a live build since the rename — treat any previously-
published `deepest`/`ancestor` split as describing the pre-rename property only.

It is a **(node, substrate) fact, not a node fact**: `2.A.1` can be `most_specific` for
one substrate and `inherited` for another, depending on which of its children survived.

This is *not* "curated vs inherited" — only `tc_specificity` nodes own
`substrate_classes` and they are leaves, so curated-vs-inherited is already exactly
`level_kind = 'tc_specificity'`. Depth is the part that needed materialising.

### 7.2 ⚠️ `Metabolite.transporter_count` — redefined

Was `tc_specificity`-only, which read **0 for 1,218 of the 1,462 transported
metabolites (83%)** — a filter written before the ancestor-only prune. Now counts
DISTINCT sources over `substrate_depth = 'most_specific'` edges: every transported
metabolite gets a non-zero count (median 1, max 123), with no ancestor counted
alongside its own descendant.

### 7.3 ⚠️ BREAKING — metabolite counts split by evidence arm

`metabolite_count` used to union catalysis with transport. The arms are not
comparable — catalysis **p90 = 11**, transport **p90 = 554** — and **23,137 genes**
had transport evidence only, so for most genes the stored number was entirely the
inflated arm with nothing marking it.

| Node | `metabolite_count` (was union) | new transport property |
|---|---|---|
| `Gene` | catalysis only | `transported_metabolite_count` |
| `Metabolite` | `gene_count` = catalysis only | `transporter_gene_count` |
| `OrganismTaxon` | `'metabolism'` arm only | `transported_metabolite_count` |

**The five chemistry counts, defined** (verbatim from
`docs/superpowers/specs/2026-08-16-vocabulary-contract-design.md` §7.1, answering
KG-IPT-004 — quotable in MCP field descriptions):

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
  `substrate_depth = 'most_specific'` (renamed 2026-08-18 from `'deepest'` —
  see §7.1 above and `docs/kg-changes/vocabulary-contract.md`).

The transport arm counts each gene's **DEEPEST TC attachments only**. 6,950 genes
are annotated at both an ancestor and its own descendant (e.g. `3.A.1` *and*
`3.A.1.14`); unioning across all attachments pulled in the superfamily's whole
rollup despite a more specific call existing. Restricting to the deepest keeps
26,813 of 26,894 genes (99.7%) and takes p90 from **554 → 97**.

Both ends use the same predicate, so `Gene.transported_metabolite_count` and
`Metabolite.transporter_gene_count` are two projections of one (gene, metabolite)
set and agree by construction. `Organism_has_metabolite`'s transport arm uses it
too — previously one gene annotated at `3.A.1` gave its whole organism all 554 ABC
substrates, which is how every organism came to "have" 63% of all metabolite nodes.

### 7.4 New: `Gene.transport_substrate_resolution`

`'resolved'` | `'family_inferred'` — categorical string, **sparse** (absent when the
gene has no TCDB edge, so "no transporter evidence" stays distinguishable from a
weak substrate claim).

| value | genes | `transported_metabolite_count` p90 |
|---|---|---|
| `resolved` | 28,405 | 33 |
| `family_inferred` | 1,671 | 554 |

`family_inferred` means the gene's deepest TC attachment is a substrate-lumping
node (`is_promiscuous`, e.g. ABC superfamily `3.A.1` with 554 substrates) — the
count is reachability, not capability; take the substrate from `product` / COG /
`function_description` instead. **This answers the friction reported by the
Alteromonas coculture analysis** ("carry a confident-vs-inferred flag on every
substrate tag"), which previously had to be derived by the caller.

**Deliberately NOT tier-gated.** Tier and substrate resolution are orthogonal:
**11,871 `resolved` genes are tier-3-only** — narrow `2.A.x` secondary carriers
where remote homology could not justify a subfamily call, and precisely the ones
TCDB resolves well. A tier gate would discard them while keeping eggNOG's equally
lumping `tc_family` edges, since eggNOG carries no `tier` at all — an artifact of
which tool made the call. Tier keeps its existing homes: the edge property, the
`annotation_types` / `annotation_quality` gate, and `Gene.tcdb_evidence_score_max`.

---

## 8. Checklist

- [ ] `tier` treated as optional (absent on 13,165 eggNOG-only edges)
- [ ] Corroboration uses `source_agreement = 'both_sources'`, **not** `size(sources) = 2`
- [ ] `tcdb_evidence_score_max` not `coalesce`d to `0.0` (sparse is meaningful; use `coalesce(..., -1.0)`)
- [ ] Scores used for ranking, not silent filtering; components surfaced
- [ ] Queries conflating "has TCDB edge" with `'tcdb' IN annotation_types` audited (15,558 genes)
- [ ] Bridge edges never traversed Pfam→TcdbFamily to assign gene function
- [ ] Hard-coded `TcdbFamily` counts / result-size limits reviewed (12,902 → 1,515 nodes; 16.8K → 53.8K edges)
- [ ] `is_promiscuous` consumers migrated — the property is **deleted** (2026-08-18); query `level >= 2 AND metabolite_count >= 50` directly
- [ ] **`metabolite_count` readers audited** — it is catalysis-only now; add `transported_metabolite_count` where the transport arm was wanted (Gene / Metabolite `transporter_gene_count` / OrganismTaxon)
- [ ] **`Metabolite.transporter_count` thresholds re-checked** — it was 0 for 83% of transported metabolites, so any `> 0` filter silently excluded them
- [ ] **Substrate queries use `substrate_depth = 'most_specific'`** (renamed 2026-08-18 from `'deepest'`), not `level_kind = 'tc_specificity'` (466 of 11,263 edges)
- [ ] **`transport_substrate_resolution` surfaced** wherever a substrate is shown — `family_inferred` (1,671 genes) means "do not read the substrate off TCDB"
- [ ] Neither new property treated as a confidence/tier signal — resolution is substrate BREADTH, orthogonal to `tier` (11,871 `resolved` genes are tier-3-only)

## See also

- [`tcdb-cazy-ontologies.md`](tcdb-cazy-ontologies.md) — the base TCDB/CAZy surface this upgrades
- [design spec §5b/§5c](../superpowers/specs/2026-08-06-tcdb-diamond-kg-integration-design.md) — decisions and measured outcomes
- [`.claude/skills/tcdb-diamond/SKILL.md`](../../.claude/skills/tcdb-diamond/SKILL.md) — regenerating the Phase-1 artifacts
