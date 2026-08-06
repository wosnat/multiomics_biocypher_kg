# TCDB-Diamond → KG Integration (Phase 2) — Design

**Date:** 2026-08-06
**Status:** Design agreed; implementation not started
**Supersedes:** the Phase-2 sketch in [`2026-05-10-tcdb-diamond-augmentation-design.md §7`](2026-05-10-tcdb-diamond-augmentation-design.md) — that sketch predates real data and picks a different shape (see §2).
**Related:** `.claude/skills/tcdb-diamond/`, `multiomics_kg/utils/tcdb_diamond.py`, `multiomics_kg/adapters/tcdb_adapter.py`, `multiomics_kg/download/build_kegg_metabolism_xrefs.py` (step 6), `docs/kg-changes/tcdb-cazy-ontologies.md`

---

## 1 — What makes TCDB different from every prior tool integration

`/integrate-a-tool` assumes a tool introduces a **new** ontology (PSORTb → `SubcellularLocalization`, SignalP → `SignalPeptideType`, InterProScan → `InterproEntry`). TCDB is the first case where:

1. The ontology **already exists** in the KG (`TcdbFamily`, 12,902 nodes), built from eggNOG's `KEGG_TC` column.
2. The new tool emits the **same vocabulary** — TCDB TC IDs — not a parallel one.

So this is not "add an ontology". It is **"add a second evidence source to an existing ontology"**, which the skill has no template for. The whole design follows from that.

### 1.1 — Why NOT two ontologies (the Pfam/InterPro analogy fails)

`Pfam` and `InterproEntry` are separate node types bridged by `Pfam_in_interpro_entry` because they are genuinely **different vocabularies** (`PF*` vs `IPR*`) that merely overlap.

eggNOG-TC and diamond-TC both emit `2.A.6.1`. Two node sets would make eggNOG's `2.A.6.1` and diamond's `2.A.6.1` *different nodes*, splitting the is-a hierarchy and duplicating all 22,483 `Tcdb_family_transports_metabolite` edges.

> **Decision D1 — One `TcdbFamily` ontology. Two sources, carried as provenance on the edge.**

---

## 2 — What the Phase-1 data actually showed (42 strains, 29,186 proteins)

The Phase-1 spec justified diamond on three goals. Real data inverts two of them.

| Phase-1 claim | Reality |
|---|---|
| "eggNOG `KEGG_TC` is family-level (3-part) only" | **False.** eggNOG: 240 distinct 5-part, 231 4-part, 57 3-part |
| "diamond unlocks the `tc_specificity` leaves" | **Barely.** diamond: 478 3-part, 414 4-part, **36** 5-part. Only 281 calls are strictly deeper than their eggNOG ancestor |
| "no verification surface" | **Confirmed — and this is the real payoff.** 9,912 proteins carry both sources |

**Diamond is a breadth source, not a depth source.** The headline win is coverage: **15,074 proteins** receive a TC call with no eggNOG TC at all, and **762 TC IDs** eggNOG never mentions.

This is why the §7 sketch is superseded: it was built around a specificity win that does not exist.

---

## 3 — Defects found in the Phase-1 artifact (must be fixed before integration)

### 3.1 — `filter_action` is incoherent

`annotate_candidate_filters` ([`tcdb_diamond.py:369`](../../../multiomics_kg/utils/tcdb_diamond.py#L369)) applies a 5-rule first-match chain. Measured against all 40,520 candidates:

| Defect | Evidence |
|---|---|
| **Verdict depends on unrelated siblings.** Rules 1, 2, 4 fire only when `len(cands) > 1`. Identical evidence → opposite verdict. | 920 lone `contradicts_both` + 239 lone `conflicts` kept purely for lacking a sibling |
| **Disqualified by a ghost.** `any_pfam_confirm`/`any_egn_confirm` are computed over all siblings *before* any are dropped, so the "supported alternative" justifying a drop may itself be dropped. | 1,452 candidates |
| **First-match-wins** makes the recorded reason arbitrary among several applicable rules. | structural |
| **Uncalibrated magic numbers** `0.20` / `0.25` drive 6,221 drops. | structural |
| **Silently vanishes proteins.** | 4,200 proteins have zero kept candidates |

The decisive number: asking the **sibling-independent** question — do *both* evidence sources contradict this candidate? — yields **25 of 40,520**. The filter excludes **12,988**. ~13K exclusions rest on reasoning the first two defects show is incoherent.

> **Decision D2 — Delete `filter_action` and `annotate_candidate_filters` entirely.**
> The tier policy (e-value ≤ 0.001, HSP ≥ 50 aa, identity/coverage bands) is already a principled, sibling-independent quality gate. That stays. Only the post-hoc chain dies.

### 3.2 — The runner has a circular dependency

```
step 2 build_gene_annotations ──> gene_annotations_merged.json
                                            │
                                            ▼  run_tcdb_diamond.py:274
                                  <strain>.tcdb.calls.json
                                            │
        ┌───────────────────────────────────┘
        ▼  Phase 2 wants to merge calls.json back into step 2
step 2 build_gene_annotations
```

The runner reads `gene_annotations_merged.json` (step 2's **output**) for Pfam annotations, and `.emapper.annotations` for eggNOG TCs. Consequences today: re-running step 2 mutates `calls.json` even when the diamond blast is byte-identical. Phase 2 would close the loop into a true cycle.

> **Decision D3 — The runner becomes a pure primary-evidence producer.**
> Inputs: `protein.faa` + `tcdb.dmnd`. Nothing else. No eggNOG, no merged annotations.

### 3.3 — Everything removed is recoverable downstream

| Field | Read from | Replacement |
|---|---|---|
| `egn_tcids` | eggNOG | Redundant — eggNOG TCs are already in the merge |
| `egn_agreement` | eggNOG | Derivable: `confirms`/`refines`/`extends`/`conflicts` are all views of the union + is-a hierarchy |
| `pfam_ids` | **merged json (cycle)** | Already `Gene_has_pfam` edges |
| `pfam_tc_families` | TCDB Pfam→TC map | → `Pfam_in_tcdb_family` bridge edges (D5) |
| `pfam_agreement` | **merged json (cycle)** | Derivable via the bridge |
| `filter_action` | all of the above | Deleted (D2) |

No information is lost — it moves from a frozen artifact tag to a queryable graph relation.

---

## 4 — Target design

### D4 — Edge model: one edge per (gene, TC ID), union of both sources

```cypher
(Gene)-[:Gene_has_tcdb_family {
    sources:           ['eggnog','diamond'],   // provenance; corroboration = size() = 2
    tier:              2,                       // diamond only (null when eggnog-only)
    confidence_score:  0.85,                    // diamond only
    identity:          50.5,                    // diamond only
    qcov:              78.0,
    evalue:            1.15e-74,
    consensus_n:       4,
    incompletely_characterized: false
}]->(TcdbFamily {tcdb_id: '2.A.6.1'})
```

- Both sources calling the same TC → **one** edge, `sources: ['eggnog','diamond']`. Cross-source corroboration for free — the one Phase-1 goal that genuinely pays off.
- Genuine disagreement → two edges to different TC IDs, both visible. Nothing silently reconciled.
- Edge ID stays `{locus_tag}-has_tcdb-{tc}`, already unique per (gene, TC).

### D5 — `Pfam_in_tcdb_family` bridge edge

TCDB's curated Pfam→TC map (~8.3K pairs / ~1.3K Pfams) is the **only** non-derivable piece — real external curation. It becomes a bridge, mirroring `Pfam_in_interpro_entry` exactly:

```cypher
(Pfam)-[:Pfam_in_tcdb_family]->(TcdbFamily)
```

Pfam corroboration becomes a query instead of a frozen tag:

```cypher
MATCH (g:Gene)-[:Gene_has_pfam]->(:Pfam)-[:Pfam_in_tcdb_family]->(t:TcdbFamily)
MATCH (g)-[e:Gene_has_tcdb_family]->(t)
RETURN g, t, e   // Pfam-corroborated calls
```

**The map download moves out of the run script into prepare_data step 6**, which already pulls 4 TCDB TSVs into `cache/data/tcdb/raw/`. This becomes the 5th.

### D6 — All candidates become edges; the tier gate applies only to rollups

With `filter_action` gone, **all 40,520 diamond candidates** become edges (up from the 27,532 the filter endorsed). The tier policy remains the quality gate, and every raw signal rides on the edge so consumers filter transparently.

| Layer | Rule |
|---|---|
| `Gene_has_tcdb_family` edges | all candidates, every tier |
| `Gene.annotation_types` / `informative_annotation_types` / `annotation_quality` | `'tcdb'` folds in from **eggNOG** or **tier ≤ 2** only |
| `Gene.tcdb_family_count` | counts all edges (routing count, not a quality signal) |

Rationale: tier 3 is conservative-by-design remote homology, not noise — median identity 34%, **median e-value 2.3e-30**, zero calls violating the 0.001 gate, and 3.9% are ≥70% identity demoted for *subfamily ambiguity*, not weak similarity. Family-level claims at 30–40% identity are defensible for structurally-conserved membrane transporters. But 22K such calls must not silently upgrade `annotation_quality`, which was calibrated against curated sources and drives `genes_by_function` routing.

---

## 5 — Known consequence: the pruning cascade

Extending `transporter_classification` with diamond's 762 new TC IDs re-seeds step 6's pruning ([`build_kegg_metabolism_xrefs.py:346`](../../../multiomics_kg/download/build_kegg_metabolism_xrefs.py#L346)), which walks **up and down** each seed. All 928 diamond TC IDs are in the curated hierarchy, so **no new `seed_aliases` are needed**. But the closure roughly doubles (~5.6K → ~13.4K on a direct model of the seeding).

Under the *current* pruning model that pushes `TcdbFamily` from 12,902 toward **~21K nodes**, most additions being zero-gene `tc_specificity` leaves. It also cascades: `subtree_substrates` → transport-reachable KEGG pathways → new `Metabolite` nodes.

**This is accepted, not solved here.** The pruning cleanup is a **separate follow-up task** (§7).

---

## 6 — Scope

**In scope:** runner purification (D2, D3); Pfam→TC map moved to step 6 (D5); merge of diamond TCs into `transporter_classification` + provenance sidecar; adapter provenance properties (D4); `Pfam_in_tcdb_family` bridge (D5); post-import tier gate (D6); kg-validity assertions; release-notes doc.

**Explicitly NOT in scope:**
- The TCDB pruning cleanup (§7).
- Reworking the tier policy thresholds.
- MCP / explorer registration.
- Any change to `Tcdb_family_transports_metabolite` semantics or the `subtree_substrates` rollup.

---

## 7 — Deferred follow-up: TCDB pruning cleanup

**94.5% of `TcdbFamily` nodes (12,198 of 12,902) have `gene_count = 0`** — 10,292 of them `tc_specificity`. Cause: the "below" arm of bidirectional pruning keeps every descendant of a gene-annotated family, *and* step 6 separately pre-rolls substrates up to ancestors via `subtree_substrates`. Two overlapping mechanisms solve the same problem; the rollup alone suffices, and it is what forced the `is_promiscuous` warning flag (33 families) onto consumers.

Dropping the "below" arm would take `TcdbFamily` to ~2,500 nodes with substrates still reachable via the ancestor rollup. Tracked as its own task so the diffs stay reviewable.

> **Ordering note:** doing the cleanup *first* would let diamond land on a lean hierarchy instead of inflating a bloated one and then deflating it. Deferral was chosen deliberately for reviewability.

---

## 8 — Acceptance criteria

1. `build_strain_calls(tsv_path)` takes **no** eggNOG / merged-annotations arguments; `grep` finds no read of `gene_annotations_merged.json` or `.emapper.annotations` under `.claude/skills/tcdb-diamond/`.
2. `filter_action`, `annotate_candidate_filters`, `compute_egn_agreement`, `compute_pfam_agreement` are gone; no consumer references them.
3. Regenerated `calls.json` for all 42 strains carries only sequence-derived fields; re-running prepare_data step 2 leaves them byte-identical.
4. `cache/data/tcdb/` carries the Pfam→TC map, downloaded by step 6, committed.
5. `Gene_has_tcdb_family` edges carry `sources`; both-source edges have `size(sources) = 2`.
6. `Pfam_in_tcdb_family` edges exist and connect existing `Pfam` nodes to existing `TcdbFamily` nodes (no dangling).
7. `Gene.annotation_types` contains `'tcdb'` only where eggNOG-sourced or tier ≤ 2.
8. `/omics-edge-snapshot`: expression edges unchanged.
9. `pytest -m "not slow and not kg"` and `pytest -m kg` green; `snapshot_data.json` regenerated.
10. `docs/kg-changes/tcdb-cazy-ontologies.md` updated (or a new `tcdb-diamond-extension.md`), CLAUDE.md "Key graph facts" + "Actual Neo4j labels" updated.

---

## 9 — Decision log

| # | Decision |
|---|---|
| D1 | One `TcdbFamily` ontology; sources as edge provenance — **not** two ontologies |
| D2 | Delete `filter_action` / `annotate_candidate_filters` entirely |
| D3 | Runner = pure sequence evidence (`protein.faa` + `tcdb.dmnd`); breaks the step-2 cycle |
| D4 | One `Gene_has_tcdb_family` edge per (gene, TC), `sources: str[]` + diamond evidence props |
| D5 | Pfam→TC map → `Pfam_in_tcdb_family` bridge edges; download moves to prepare_data step 6 |
| D6 | All candidates → edges; tier ≤ 2 gate applies **only** to `annotation_types` / `annotation_quality` |
| D7 | TCDB pruning cleanup deferred to a separate task |
| D8 | Rollout: rewrite → rerun MED4 → review → batch remaining 41 |
