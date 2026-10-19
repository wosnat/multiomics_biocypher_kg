# TCDB diamond best-hit details on `Gene_has_tcdb_family`

**Date:** 2026-10-19
**Spec:** [`docs/superpowers/specs/2026-10-18-tcdb-diamond-best-hit-design.md`](../superpowers/specs/2026-10-18-tcdb-diamond-best-hit-design.md)
**Plan:** `docs/superpowers/plans/2026-10-19-tcdb-diamond-best-hit.md`
**Track:** transport programme, project 2 (slice 1 = [GO transport substrates](go-transport-substrates.md))
**Rule version:** `best_hit_v1`

## What changed

`/tcdb-diamond` blasts every protein against TCDB, then the tier policy truncates the call to
the depth the evidence justifies — usually the 3-part family. Through the pre-rolled
`Tcdb_family_transports_metabolite` edges a family-level call reads as every substrate of every
system in the family (`3.A.1` lists 554). The protein's **actual best-hit system** was computed
and thrown away. It is now kept, as six sparse properties on the diamond-sourced
`Gene_has_tcdb_family` edge in the top hit's family:

| Property | Type | Content |
|---|---|---|
| `best_hit_tcids` | str[] (1–3) | Best-hit 5-part TCDB systems, best bitscore first. Near-ties within 90% of the top bitscore, from any family. Plain TCDB ids — generally not `TcdbFamily` nodes |
| `best_hit_names` | str[] | TCDB name of each system, aligned with `best_hit_tcids`; the TC id when TCDB gives none |
| `best_hit_identity` | float | Percent identity of the top system |
| `best_hit_substrate_tcids` | str[] | Long-form table, one entry per (system, substrate): the system |
| `best_hit_substrate_ids` | str[] | …the substrate, as a canonical Metabolite / ChemicalRole id |
| `best_hit_substrate_names` | str[] | …its name (node name, else the TCDB name) |

Example, MED4 cynD (PMM0372), tier call `tcdb:3.A.1`:
`best_hit_tcids = ['3.A.1.16.2']`, `best_hit_names = ['Bispecific cyanate/nitrite transporter.']`,
`best_hit_identity = 59.2`, substrates nitrite + cyanate.

Nothing else in the graph changed: no node or edge was added or removed, no existing property
moved. The vocabulary gained six `ControlledVocabulary` entries (so
`Schema_info.controlled_vocabularies_hash` moved) and `data_source:tcdb_diamond.info_types` gained
`tcdb_best_hit_ids`.

One existing property did change, as a bug fix: 31 `TcdbFamily.name` values that the T6 scraper
had cut mid-parenthesis (a bare `<` in TCDB text was read as an HTML tag) were re-derived — for
example `3.A.1.4.4` "The high-affinity (" is now "The high-affinity (<1 μM) urea porter".

## The rule (`best_hit_v1`)

Computed from the committed `<strain>.tcdb.tsv` by `select_best_hits`
(`multiomics_kg/utils/tcdb_diamond.py`):

1. **Gate:** query coverage ≥ 75% **and** reference coverage ≥ 75% (gapseq's reference-coverage
   rule), TC class 1–4. No identity or length floor — identity is published instead.
2. **Reciprocal best hit:** the protein's top hit must also pick this protein as its best gated
   match among the strain's proteins. Removes paralogs borrowing a sibling's substrate and
   non-transporters that hit TCDB's protein-import machinery (clpB, dnaK, ftsH).
3. **Window:** gated hits within 90% of the top bitscore, up to 3 distinct systems, any family.

The record sits on the diamond edge in the top hit's 3-part family (`attach_tcid`). 188 proteins
have a best hit but no tier call (top alignment < 50 aa — small photosystem and complex
subunits); they have no edge, so the best hit lives only in the merged annotations'
`tcdb_best_hit_ids`.

## Caveat paragraph for the explorer

> `best_hit_*` on a `Gene_has_tcdb_family` edge is a raw clue: the curated TCDB system(s) this
> gene's protein most resembles, and their substrates as TCDB lists them. It is not a call — the
> edge's tier and TC id stay the conservative claim. Read `best_hit_identity` as how far the
> substrates can be trusted for this gene (roughly ≥ 50–60% fairly safe, < 40% weak). One id is a
> clean best hit; several are a near-tie — check whether their substrates agree (read the three
> `best_hit_substrate_*` lists together: they are one table, one row per system and substrate).
> Near-ties can come from another TC family. Substrates are unfiltered (H⁺, Na⁺, generic classes,
> protein/DNA appear) and an id need not be a node. No substrate lists means TCDB lists none for
> those systems. Not the same as MEROPS `best_hit_id`, which is a single hit.

## Query examples

```cypher
// the specific system behind a family-level call
MATCH (g:Gene {locus_tag: 'PMM0372'})-[r:Gene_has_tcdb_family]->(t)
WHERE r.best_hit_tcids IS NOT NULL
RETURN t.tcdb_id, r.tier, r.best_hit_tcids, r.best_hit_names, r.best_hit_identity

// (system, substrate) rows for one gene
MATCH (g:Gene {locus_tag: 'PMM0372'})-[r:Gene_has_tcdb_family]->()
WHERE r.best_hit_substrate_ids IS NOT NULL
UNWIND range(0, size(r.best_hit_substrate_ids) - 1) AS i
RETURN r.best_hit_substrate_tcids[i] AS system, r.best_hit_substrate_ids[i] AS id,
       r.best_hit_substrate_names[i] AS substrate

// genes whose best hit names urea, per organism
MATCH (g:Gene)-[r:Gene_has_tcdb_family]->()
WHERE any(s IN r.best_hit_substrate_names WHERE toLower(s) = 'urea')
RETURN g.organism_name, count(g) ORDER BY count(g) DESC
```

## Validation (live graph, 2026-10-19)

`uv run python scripts/tcdb_best_hit_validation.py`. The rule was frozen before this ran; MED4
(where it was measured) is training, the other 42 strains hold-out. The comparators share
homology with the TCDB calls, so agreement is an upper bound.

**Build check.** 12,231 edges carry the record (= the count predicted from calls.json + merged
annotations; 0 dropped by pruning, 0 fallback attachments). All 18,763 (edge, substrate) pairs
point at an existing Metabolite / ChemicalRole node. `Gene_has_tcdb_family` (54,627),
`TcdbFamily` (1,515), `Metabolite` (3,243) and `Tcdb_family_transports_metabolite` (10,790)
are unchanged; `/omics-edge-snapshot` shows no regressions.

**Agreement with existing routes (hold-out, 12,074 edges):**

| Comparator | Result |
|---|---|
| eggNOG TC (7,407 edges with one) | 93.7% consistent: same system 853, best hit inside eggNOG's call 5,405, same family other branch 680; different family 469 (6.3%) |
| Pfam bridge of the best-hit family (8,757 where one exists) | 65% corroborated by the gene's own Pfams |
| KEGG BRITE transporter KO (`ko02000`) | 45% of edges (many genes carry no KO) |
| Cyanorak role naming transport (3,481 Pro/Syn edges with a role) | 58% |
| GO transport substrates (3,137 edges where both name substrates) | 47% share a substrate; disjoint cases are mostly naming granularity and unfiltered H⁺ / Na⁺ |

**Strength profile.** Best-hit identity < 40%: 5,292 edges (44%); 40–70%: 5,438; ≥ 70%: 1,344.
3,476 edges (29%) carry more than one system (near-tie). 3,466 substrate pairs are ids in
slice 1's GO target filter (H⁺, Na⁺, upper-level classes) — kept, unfiltered, by design.

**eggNOG different-family sample (30, seed 0, classified by hand):** 18 both valid — mostly one
protein TCDB curates under two systems in different families (LptB in Lpt export, YidC/ALB3 in
the Sec translocase, MexE/OprN in MexEF-OprN, CzcC, FlhA, PilQ) or eggNOG's class 8/9 placement
the rule excludes; 3 eggNOG id stale or malformed (`2.A.1.40` for NCS2 permeases, `9.A.58`);
**4 best hit wrong** (CcdA ×2 whose true family is class 5, AtpI, a K⁺ channel) — all at
24–41% identity, which is what the identity caveat is for; 5 unclear.

**Reciprocal check.** 2,993 gate-passing proteins fail it (MED4 17). On hold-out strains about
1,300 have transporter-like products (paralogous ABC components in heterotrophs) — a recall
cost, recorded for project 3, not tuned here.

**Nitrogen worked example** (MED4, MIT9313, WH8102, HOT1A3). The family-level calls that read
as 518 substrates now resolve: cynA/B/D → `3.A.1.16.2` cyanate + nitrite (MED4 59–71%, WH8102
61–69%); urtA–E → `3.A.1.4.4` urea (MED4, MIT9313, WH8102; 47–62%); amt1 → `1.A.11.2.5`
ammonia (81–88%); HOT1A3 nrtC / nasD → `3.A.1.16.1` nitrate + nitrite (35–55%); dpp/opp →
peptide systems (27–36%). Weak or misleading rows are visible by identity: MED4/MIT9313/WH8102
`clc2`/`eriC` (ClC chloride antiporter, 28–30%) list nitrate among ClC substrates; `acs`
(26%) hits a carnitine CoA ligase. Gaps that are not the rule's: WH8102's `nrtP` (MFS nitrate
transporter) has no diamond hit at all, so it has no call to refine; the SYNW1415–1417
nitrate/sulfonate/bicarbonate ABC resolves to taurine/sulfonate systems (`3.A.1.17`, 34–44%).

## Pipeline

calls.json `best_hit` (`run_tcdb_diamond.py --from-tsv`) → step 2 merge `tcdb_best_hit_ids` →
TCDB name scraper (scope includes best-hit systems; merge-preserving) → step 6
`tcdb_pruned.json` `best_hit_substrates` (separate resolve; `kegg_data.json` untouched) →
`tcdb_adapter`. Full data-flow diagram in the spec.
