# TCDB diamond best-hit details — design

**Date:** 2026-10-18 · **Status:** draft for review · **Project 2 of the transport programme**
(programme framing: `2026-10-18-go-transport-substrates-design.md` §1; slice 1 merged c52c3580)

## 1. Why

`/tcdb-diamond` blasts every protein against TCDB, then the tier policy plus per-family
consensus collapse truncates the call to the depth the evidence "justifies" — usually
`tc_family` (tier 3). The KG keeps only that truncated id. Through
`Tcdb_family_transports_metabolite` (pre-rolled up over the full subtree) a tier-3 call
reads as every substrate of every system in the family: up to ~500 substrates per gene
(`3.A.1` = 554).

The protein's **actual best-hit system** — the curated 5-part TCDB system it most
resembles — is computed and then thrown away. It is not in the KG, so no consumer can
recover it. Example: MED4 `cynD` (PMM0372) is called `tcdb:3.A.1`; its best hit is
`3.A.1.16.2` (cyanate) at 59% identity, 95% query coverage.

This project keeps that one clue. It does **not** decide what a gene transports; it adds
the best hit as one more route, with its caveat, beside eggNOG TC, KEGG BRITE
transporters, Pfam/InterPro domains, Cyanorak roles and GO transport (slice 1). The
reader (an LLM agent via the explorer) weighs routes; the final call stays in the
analysis.

Motivating use (evaluation only, not solved here): "which transporters does this genome
have for nitrogen sources, and what are their substrates or substrate classes"
(ammonium, urea, cyanate vs nitrate/nitrite, amino acids, peptides…).

## 2. Goals, non-goals, anti-goals

**Goals**
- Publish each diamond call's gated best-hit TCDB system(s) and their own substrates,
  on the existing diamond-sourced `Gene_has_tcdb_family` edge.
- Smallest form that survives the project-3 TCDB redesign: six edge properties (ids +
  aligned names, so the record reads without joins), no new node or edge types, no new
  Metabolite nodes.
- Cheap validation against transport evidence already in the KG.

**Non-goals (deferred)**
- A per-gene `Gene → Metabolite` transport edge (after project 3).
- Substrate-class fallback for ambiguous/low-identity calls (derivable later from what is
  stored here).
- Grouping subunits into systems; transport direction.
- Changing tier / consensus / `tcdb_family` attachment, `substrate_depth`,
  `substrate_breadth`, `transport_substrate_resolution`, or any rollup. Project 3 owns these.
- Answering the nitrogen question.

**Anti-goals**
- Storing raw diamond hits (no hit tables, no runner-up scores, no qcov/scov/bitscore per hit).
  (Within the six properties nothing is dropped or filtered — they are raw material.)
- Hand-picked TC-family exclusion lists derived from the genes we inspected (overfitting).
- A fused "confidence" flag that hides the evidence (programme rule: routes with caveats).
- Tuning thresholds on the evaluation strains.

## 3. The rule (frozen before evaluation)

Input: the committed per-strain diamond table
`cache/data/<org>/genomes/<strain>/tcdb/<strain>.tcdb.tsv` (columns: qid, subject
`<acc>-<TCID>`, identity, qcov, scov, length, evalue, bitscore; `--max-target-seqs 25`,
`--evalue 0.001`, `--more-sensitive`). No new diamond run.

For each protein:

1. **Hit gate.** Keep hits with qcov ≥ 75 **and** scov ≥ 75, and TC class ∈ {1, 2, 3, 4}.
   - qcov ≥ 75: measured on MED4 2026-10-18 (hit covers most of our protein).
   - scov ≥ 75: gapseq's reference-coverage rule (`src/transporter.sh`, jotech/gapseq
     master 2dce3a0: TCDB is the blast *query*, `-qcov_hsp_perc 75`, i.e. coverage of the
     reference). Rejects single-shared-domain hits.
   - Classes 1–4: gapseq's `TC_types` (5, 6, 8, 9 ignored: electron carriers,
     accessory/incompletely characterized). Excludes our earlier "8/9 only" by also
     dropping 5 and 6.
   - gapseq's bitscore 200 is **not** adopted (length-biased: on MED4 it drops 26
     proteins with ≥ 50% identity); no identity floor (identity is published instead).
   - Subject format is `lcl|<acc>-<TCID>`; reuse `parse_tcdb_subject_id`. One TCDB
     reference carries a 6-part id (`lcl|Q39253-1.A.1.9.2.3`, 6 of 261,292 rows
     corpus-wide); truncate any id longer than 5 parts to its first 5 (`1.A.1.9.2`, a
     curated system with substrates). The parser must accept it rather than reject it.
2. **Reciprocal best hit.** The protein's top-bitscore gated hit `X` (a subject row) must
   also be `X`'s top-bitscore hit among this strain's **gated** rows (same table, grouped
   by subject; a partial-coverage paralog is not a competitor). Fails → no best-hit
   details for this protein. Standard ortholog-vs-paralog test, no parameter. Approximate
   because the forward run keeps 25 targets per protein; the error leans toward passing.
   Unit = query id (WP_ accession): identical protein copies share one WP_ and therefore
   one result, applied to each of their genes. Ties (my call): order by bitscore desc,
   then TCID / query id ascending — deterministic, no biological meaning.
3. **Window.** Keep gated hits with bitscore ≥ 0.9 × top bitscore, up to 3, best first,
   distinct 5-part TCIDs (a TCID appearing twice keeps its best row), **from any TC
   family** (user, 2026-10-18). Cross-family near-ties are 244 of 12,379 proteins (2%),
   mostly one protein curated as a component of systems in two families (3.D.3 + 3.E.2
   b6f/reaction centre 70; 1.B.42 + 3.A.1 Lpt export 63; 3.A.11 + 3.A.15 22) — real
   information, not competing guesses.
4. **Attach.** The top hit decides the edge: the record goes on the protein's
   diamond-sourced `Gene_has_tcdb_family` edge in the top hit's 3-part family (if
   several, the deepest ancestor-or-self of the top hit). If the gene has no diamond edge
   in that family (188 of 12,379 proteins), the record goes on its highest-confidence
   diamond edge instead. **Nothing is dropped** (user, 2026-10-18: this is raw material for
   project 3, not a finished feature). A gene can have several diamond edges (one per
   called family; 6,752 of 29,642 proteins have more than one); the record sits on exactly
   one of them. Consequence: window ids — even the top one, in the fallback case — are
   **not** guaranteed to fall under the edge's TC id; consumers must not assume it.

The rule is frozen; the counts below are re-measured by the plan's first task from a
committed script implementing exactly this text, and that output becomes the reference.

Effect measured 2026-10-18 (exploratory, gate as above):
- MED4: 180 proteins pass the gate, 163 pass the reciprocal check, 46 of those carry more
  than one window TCID. The 17–18 failures are
  mostly non-transporters (clpB, dnaK, ftsH, gap3 hitting protein-import/secretion
  systems), photosystem paralogs (pcb, ndhD2) and ABC paralogs borrowing a sibling's
  substrate (sufC, PMM1029/1031, modF). cyn/urt/amt1/dpp/ktrA keep their calls.
- 43 strains: ~12.4K proteins receive details, naming ~1,990 distinct 5-part systems
  (~200 of them already TcdbFamily nodes).
- Builder check found that about a third of gate-passing 5-part systems have no
  substrates in TCDB (1,533 of 4,594 pre-reciprocal), so many edges will carry
  `best_hit_tcids` without `best_hit_substrate_ids`.

### Relation to the tier policy

Tier expresses uncertainty by truncating depth (weak hit → shallower claim); the
best-hit record keeps full depth and publishes identity + near-ties instead. The gate is
stricter on coverage (both sides ≥ 75 vs tier 1/2's qcov ≥ 70/60 and tier 3's qcov **or**
scov ≥ 40) and on paralogs (reciprocal), looser on identity (none). Exploratory overlap
over 43 strains, per diamond call, best hit in the same 3-part family:

| Edge tier | Calls | With details | Reciprocal fail | Coverage/class fail | Best-hit identity median (p10–p90) |
|---|---|---|---|---|---|
| 1 | 192 | 138 | 3 | 51 | 74 (70–96) |
| 2 | 6,693 | 3,271 | 359 | 3,063 | 53 (42–74) |
| 3 | 34,277 | 8,782 | 2,630 | 22,865 | 37 (26–67) |

All 138 tier-1 details name the tier call's own system; all 3,271 tier-2 details fall
inside the called subfamily. The new information is on tier 3 (776 at ≥ 70% identity,
demoted by consensus; 5,292 below 40%). 864 tier-2 and 2,505 tier-3 calls carry more than
one window system. Tier-1 calls failing the gate get no details — no exception (user,
2026-10-18).

## 4. What is stored

On the diamond-sourced `Gene_has_tcdb_family` edge (sparse; absent when the rule yields
nothing):

| Property | Type | Content |
|---|---|---|
| `best_hit_tcids` | str[] (1–3) | Window 5-part TCIDs, best bitscore first. Plain TCDB ids (`3.A.1.16.2`), not node ids — these systems are generally not KG nodes. |
| `best_hit_names` | str[] | Aligned index-for-index with `best_hit_tcids`: TCDB system name from `tcdb_names.json` (e.g. `3.A.1.16.2` → "Bispecific cyanate/nitrite transporter."); falls back to the TCID itself when unnamed (same marker as `TcdbFamily.name = tcdb_id`). |
| `best_hit_identity` | float | Percent identity of the top hit. |
| `best_hit_substrate_tcids` | str[] | Long form, one entry per (window system, substrate) pair: the system the substrate belongs to. Aligned index-for-index with the two lists below; ordered by window rank, then substrate id. A substrate two systems share appears once per system; a system with no substrates has no entries. |
| `best_hit_substrate_ids` | str[] | Aligned with `best_hit_substrate_tcids`: that system's own substrate (`substrates.tsv` row for the TCID), resolved like every TCDB substrate in step 6 (secondary → primary ChEBI, KEGG/ChEBI/MNX canonicalisation, role partition) to a canonical Metabolite / ChemicalRole id. **Unfiltered and kept whether or not a node exists** (user, 2026-10-18: raw material for project 3): H⁺, Na⁺, upper-level classes and polymers (protein, DNA) stay in; an id the resolver cannot canonicalise is stored as its primary `chebi:` CURIE. |
| `best_hit_substrate_names` | str[] | Aligned with the two lists above: the canonical node's name when the node exists, else the name TCDB gives with the ChEBI id in `substrates.tsv`. |

Layout (user, 2026-10-19): the substrate lists are a long table, not a union, so the
reader can tell which system names which substrate — the point of a split call (MED4
`sul1`: one system → sulfate, the other → fumarate). Neo4j has no nested lists, so this is
the flat form; it can change when the design matures (project 3).

Read together: one TCID = clean call; several = near-tie, and the reader checks whether
their substrates agree (e.g. MED4 `sul1`: fumarate vs sulfate → genuinely ambiguous, not
"transports both"). Identity says how far the substrate transfer can be trusted (rough
literature guide: ≥ ~50–60% fairly safe, < ~40% weak).

Absent `best_hit_substrate_ids` means TCDB lists no substrate for those systems (about a
third of gate-passing systems). Because nothing is filtered, consumers must expect H⁺,
Na⁺, generic classes and polymers among the ids; filtering (e.g. slice 1's
`go_transport_target_filter.yaml`) is a reading choice, stated in the kg-changes caveat.

Note: the existing edge `identity` / `qcov` come from the highest-identity hit in the
consensus family group, pre-gate; `best_hit_identity` is the gated top-bitscore hit. They
usually coincide; both stay.

**No new Metabolite nodes, no dangling edges.** The ids live in a property, not in
edges, so an id without a node cannot break import. Most will match emitted nodes anyway
(the top hit descends from a kept TC family whose node carries the full-subtree substrate
rollup; 0 gated TCIDs were missing from the step-6 hierarchy on 2026-10-18). Validation
reports how many stored ids have no node.

## 5. Where it is computed (my call; builder may adjust)

- **Phase-1 artifact:** `run_tcdb_diamond.py` adds a per-protein `best_hit` block to
  `<strain>.tcdb.calls.json` (`tcids`, `identity`, plus the rule version), derived from the
  same hit table it already parses. Diamond runs fast, so a plain `--force` rerun of all
  strains is acceptable (user, 2026-10-18); no normalize mode needed. Diamond parameters
  stay as they are (`--max-target-seqs 25` etc.): raising them would change the existing
  tier/consensus calls, which are out of scope. Rule lives in one function with unit tests.
- **Step 6:** collects all best-hit TCIDs across strains, resolves their own substrates
  through the existing substrate resolver, and writes `best_hit_substrates:
  {tcid5: [canonical_id, …]}` into `tcdb_pruned.json` (unfiltered, §4). Canonical ids exist
  only after the step-6 identity pipeline, so the map is written from there, not where
  `subtree_substrates` is first built.
- **TCDB names:** `scrape_tcdb_names.py` keeps the name of **every** 5-part system on each
  family page it parses (today it keeps only systems that are KG nodes: 203 of ~1,990
  best-hit systems named), and is topped up with `--families` for the 251 families the
  best-hit systems fall in (≈ 200 not yet scraped, ~10 min at the 2.5 s cadence; the raw
  page cache is reused). `tcdb_names.json` stays committed. Step 6 adds the names to
  `best_hit_substrates`' sibling map `best_hit_names: {tcid5: name}` in `tcdb_pruned.json`
  (or the adapter reads `tcdb_names.json` directly — builder's choice) and warns on
  best-hit systems still unnamed.
- **Adapter:** `tcdb_adapter` reads `best_hit` from calls.json (already reads it for
  evidence) and `best_hit_substrates` from `tcdb_pruned.json`, attaches the six properties
  per §3 step 4.
- Merge (`gene_annotations_merged.json`): untouched.

Vocabulary/contract: register the six properties in `config/controlled_vocabularies.yaml`
if numeric/list edge props are registered there (KG-SYNC-005 practice); sanitize strings.

### Data flow

New pieces are marked **(new)**; everything else exists today and is unchanged.

```mermaid
flowchart TD
    FAA["protein.faa<br/>(per strain)"] --> DIA["diamond blastp vs TCDB FASTA<br/>(run_tcdb_diamond.py, params unchanged)"]
    DIA --> TSV["&lt;strain&gt;.tcdb.tsv<br/>raw hits, committed"]
    TSV --> TIER["tier policy + consensus collapse<br/>(unchanged)"]
    TSV --> RULE["best-hit rule best_hit_v1 (new)<br/>gate · reciprocal · window"]
    TIER --> CALLS["&lt;strain&gt;.tcdb.calls.json<br/>calls[] + best_hit {tcids, identity, attach_tcid} (new)"]
    RULE --> CALLS

    CALLS --> S2["step 2 merge<br/>tcdb_diamond_ids (unchanged)"]
    S2 --> MERGED["gene_annotations_merged.json"]

    CALLS -->|best-hit TCIDs| SCRAPE["scrape_tcdb_names.py<br/>scope += best-hit families,<br/>keep every system name (new)"]
    TCDBWEB[("tcdb.org family pages")] --> SCRAPE
    SCRAPE --> NAMES["tcdb_names.json (committed)"]

    CALLS -->|best-hit TCIDs| S6["step 6<br/>own substrate strings per TCID →<br/>separate resolve → canonical ids (new)"]
    SUBS[("substrates.tsv<br/>via tcdb_hierarchy.json")] --> S6
    IDP["identity pipeline<br/>(ChEBI, MNX, canonical_of)"] --> S6
    S6 --> PRUNED["tcdb_pruned.json<br/>+ best_hit_substrates {tcid: [[id, name]]} (new)"]
    MERGED --> S6

    CALLS --> ADP["tcdb_adapter"]
    MERGED --> ADP
    NAMES --> ADP
    PRUNED --> ADP
    ADP --> EDGE["Gene_has_tcdb_family edge (diamond-sourced)<br/>+ 6 best_hit_* properties (new)"]
    EDGE --> KG[("Neo4j KG")]
```

| Stage | Reads | Writes | What it contributes |
|---|---|---|---|
| Diamond runner | `protein.faa`, TCDB FASTA | `.tcdb.tsv`, `.tcdb.calls.json` | Per protein: `best_hit.tcids` (≤ 3), `best_hit.identity`, `best_hit.attach_tcid` (which of its calls carries the record), `best_hit.rule` |
| Name scraper | calls.json best-hit TCIDs, tcdb.org | `tcdb_names.json` | Names for best-hit systems (≥ 95% target; else the TCID is the name) |
| Step 6 | calls.json best-hit TCIDs, hierarchy substrate strings, ChEBI/MNX, identity pipeline | `tcdb_pruned.json` `best_hit_substrates` | Per system: its own substrates as `[canonical id, name]`, unfiltered, nodes not required. Uses a separate resolve so `kegg_data.json` and every existing key stay byte-identical |
| Step 2 merge | calls.json | merged annotations | Unchanged — the edge set (which TC ids a gene gets) comes from here as before |
| `tcdb_adapter` | calls.json, merged annotations, `tcdb_names.json`, `tcdb_pruned.json` | `Gene_has_tcdb_family` edges | Strain adapter puts `best_hit_tcids` + `best_hit_identity` on the edge whose TC id = `attach_tcid`; Multi adapter adds `best_hit_names` and the three long-form substrate lists, and carries `best_hit_*` through the seed-alias edge merge |

Rerun order after a change: diamond runner → name scraper (only if new best-hit
families) → `prepare_data.sh --steps 2 6 --force` → KG rebuild. Step 6 warns when
best-hit systems are unnamed.

## 6. Validation (light; existing sources; no tuning)

Rule frozen as §3 before running. MED4 is training data; report on hold-out strains
(e.g. MIT9313, WH8102, HOT1A3, KT2440, DSS-3) and corpus-wide. One throwaway-quality
script (kept under `scripts/` or the plan's scratch), output summarized in the kg-changes doc.

Per gene with best-hit details, measure agreement / inclusion with:
1. **eggNOG TC** — best-hit TCID equal to, descendant of, or ancestor of eggNOG's TC id;
   or contradicting (different family).
2. **KEGG BRITE transporters** (`ko02000` via `Kegg_term_in_brite_category`) — gene's KO
   is a transporter; where the KO name states a substrate, does it overlap.
3. **Domains** — best-hit 3-part family bridges (`Tcdb_family_has_pfam_domain`) to a Pfam
   the gene carries.
4. **Cyanorak role** — gene in a transport role (Pro/Syn).
5. **GO transport (slice 1)** — gene→GO→Metabolite substrates overlap
   `best_hit_substrate_ids` (same / related / disjoint).
6. **Nitrogen sanity table** (a few strains, e.g. MED4, MIT9313, WH8102, one Alteromonas):
   for N-source genes, best-hit substrates vs today's tier-rollup substrate count.

Report also: genes gaining details, reciprocal-check failures by category, split calls
(window members with disjoint substrates), fallback attachments, stored substrate ids with
no node, and how many would be removed by slice 1's target filter. Disagreements:
random sample (~30), classified rule-wrong / comparator-wrong / unclear. If anything tempts
a threshold change, re-evaluate on fresh strains. Comparators are correlated (shared
homology), so agreement is an upper bound.

**Acceptance:** every best-hit protein's record present on exactly one of its gene's diamond edges in a live rebuild; cynA/cynD →
cyanate system ("Bispecific cyanate/nitrite transporter."), amt1 → ammonium, urt → urea on
MED4; ≥ 95% of stored best-hit TCIDs carry a real name (not the id fallback); validation table recorded;
`pytest -m "not slow and not kg"` and `pytest -m kg` green; no existing edge/property
changes (`/omics-edge-snapshot` and TCDB edge counts unchanged).

## 7. Documentation

kg-changes doc (`docs/kg-changes/tcdb-diamond-best-hit.md`) with the caveat paragraph for
the explorer; CLAUDE.md TCDB bullet; tcdb-diamond SKILL.md (artifact schema); CHANGELOG
`### Added`; backlog transport section.

## 8. Left to the builder

Exact function boundaries, test fixtures, the normalize-mode flag name, whether the
`best_hit` block also records the gated top hit's qcov/scov for audit in the artifact
(artifact only, never the KG), log wording.
