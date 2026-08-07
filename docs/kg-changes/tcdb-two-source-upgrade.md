# TCDB two-source upgrade — MCP / explorer contract

**Status:** LANDED 2026-08-07 · branch `tcdb-diamond-cleanup`
**Audience:** the MCP / explorer team. This is an **upgrade delta**, not a full
description — it assumes the existing TCDB surface documented in
[`tcdb-cazy-ontologies.md`](tcdb-cazy-ontologies.md) is already implemented in the API.
**Design:** [`2026-08-06-tcdb-diamond-kg-integration-design.md`](../superpowers/specs/2026-08-06-tcdb-diamond-kg-integration-design.md)

---

## TL;DR — what breaks, what's new

| | |
|---|---|
| ⚠️ **Behaviour change** | `'tcdb' IN g.annotation_types` no longer follows from "has a `Gene_has_tcdb_family` edge" — it is tier-gated. Affects **15,558 genes**. |
| ⚠️ **Semantics change** | `TcdbFamily.is_promiscuous` is now level-gated (class/subclass always false). Family-level flags are **unchanged**. Node set went 12,902 → 1,515, so any hard-coded count is stale. |
| ⚠️ **Scale change** | `Gene_has_tcdb_family` 16,806 → **53,763** edges; genes with a TC call 11,103 → **30,076**. Result-size assumptions and default limits need review. |
| ✅ **New** | Edge provenance (`sources`) + diamond evidence + an advisory `tcdb_evidence_score`. |
| ✅ **New** | 4 ontology→ontology bridge edges (TcdbFamily → Pfam / GO). |
| ✅ **New** | `Gene.tcdb_best_evidence_score`. |

---

## 1. `Gene_has_tcdb_family` now has two evidence sources

eggNOG's `KEGG_TC` (ortholog **transfer**) and `/tcdb-diamond` (direct **sequence
similarity**) both emit TCDB TC IDs. They feed **one ontology** and **one edge per
(gene, TC id)** — never parallel edges.

```cypher
(Gene)-[:Gene_has_tcdb_family {
    sources:          ['eggnog','diamond'],  // always present
    tier:             2,                      // ↓ diamond only — SPARSE, absent on eggNOG-only
    confidence_score: 0.85,
    identity:         50.5,
    qcov:             78.0,
    evalue:           1.15e-74,
    consensus_n:      4
}]->(TcdbFamily)
```

| `sources` | Edges |
|---|---|
| `['diamond']` | 36,957 |
| `['eggnog']` | 13,165 |
| `['eggnog','diamond']` | 3,641 |

**Do not assume `tier` exists.** It is absent on all 13,165 eggNOG-only edges. Use
`coalesce(r.tier, …)` or branch on `'diamond' IN r.sources`.

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
precomputed `r.agrees_across_sources` (below) rather than recomputing.

---

## 2. Advisory scores — `tcdb_evidence_score` and the gene rollup

**Edge property `tcdb_evidence_score` (int 0–5)** — five independent signals, `+1` each:

| Component | Stored as |
|---|---|
| eggNOG called it | `'eggnog' IN sources` |
| sources agree (hierarchically) | `agrees_across_sources` (bool) |
| strong direct sequence evidence | `tier <= 2` |
| a gene Pfam is curated into this family | `pfam_corroborated` (bool) |
| a gene GO term is curated onto this family | `go_corroborated` (bool) |

Distribution: 0 → 17,422 · 1 → 9,039 · 2 → 9,265 · 3 → 7,812 · 4 → 9,144 · 5 → 1,081.

**Gene property `tcdb_best_evidence_score` (int 0–5)** = `max()` over the gene's edges.
Answers *"is this gene a transporter at all"*, where the edge score answers *"is THIS
assignment right"*. Not a copy: 32.6% of annotated genes carry calls at differing scores.

> **SPARSE.** Set only on genes with ≥1 TCDB edge. A gene with **no** transporter
> evidence is deliberately distinct from one with **weak** evidence — do not
> `coalesce(..., 0)`, that collapses the distinction. Use `coalesce(..., -1)` if a
> total order is required.

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

## 4. ⚠️ `is_promiscuous` — corrected semantics

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
  `pfam_corroborated` / `go_corroborated` do).
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
`Metabolite.transporter_count` drops correspondingly.

---

## 7. Checklist

- [ ] `tier` treated as optional (absent on 13,165 eggNOG-only edges)
- [ ] Corroboration uses `agrees_across_sources`, **not** `size(sources) = 2`
- [ ] `tcdb_best_evidence_score` not `coalesce`d to 0 (sparse is meaningful)
- [ ] Scores used for ranking, not silent filtering; components surfaced
- [ ] Queries conflating "has TCDB edge" with `'tcdb' IN annotation_types` audited (15,558 genes)
- [ ] Bridge edges never traversed Pfam→TcdbFamily to assign gene function
- [ ] Hard-coded `TcdbFamily` counts / result-size limits reviewed (12,902 → 1,515 nodes; 16.8K → 53.8K edges)
- [ ] `is_promiscuous` consumers unaffected (meaning preserved), but any class/subclass-level use now returns false

## See also

- [`tcdb-cazy-ontologies.md`](tcdb-cazy-ontologies.md) — the base TCDB/CAZy surface this upgrades
- [design spec §5b/§5c](../superpowers/specs/2026-08-06-tcdb-diamond-kg-integration-design.md) — decisions and measured outcomes
- [`.claude/skills/tcdb-diamond/SKILL.md`](../../.claude/skills/tcdb-diamond/SKILL.md) — regenerating the Phase-1 artifacts
