# InterProScan → KG integration — design (Phase 2)

**Date:** 2026-07-26
**Skill:** `/integrate-a-tool interproscan`
**Phase-1 doc:** `docs/superpowers/specs/2026-07-22-interproscan-domains-design.md`
**Status:** design agreed; implementation in progress

## 0. Intent

Wire the committed per-strain `<strain>.interproscan.calls.json` artifacts (42
strains) into the live Neo4j graph as an **InterPro-entry ontology** with a
per-gene, coordinate- and evidence-bearing membership edge.

### What the tool predicts + value space
Per protein, a list of member-DB matches (Pfam, NCBIfam, PANTHER, Gene3D,
PROSITE, SFLD, HAMAP, …), each mapped (when integrated) into a unified
**InterPro entry** `IPRxxxxxx`. Each entry has a `type`
(FAMILY / DOMAIN / HOMOLOGOUS_SUPERFAMILY / CONSERVED_SITE / REPEAT /
BINDING_SITE / ACTIVE_SITE / PTM) and sits in a shallow parent/child hierarchy.
Distinct entries observed across strains: thousands (MED4 alone 3,673;
KT2440 7,189).

### calls.json shape + join key + nulls
- Top-level dict **keyed by RefSeq WP_ accession** (= Gene `protein_id`, exactly
  like eggNOG). Absence of a key = "not processed"; `match_count == 0` = "no
  domain found".
- Per protein: `interpro_entries` (list[str]), `libraries` (list[str]),
  `matches` (one record per match×location). Per match: `library`,
  `signature_accession`, `interpro_accession` (**null when unintegrated**),
  `interpro_type` (null when unintegrated), `start`, `end`, `evalue`
  (**null for profile/pattern hits**), `score` (**null when N/A**),
  `go_terms`, `pathways`.
- **GO terms and pathways are EMPTY in every committed artifact** — the Phase-1
  run did not enable the GO/pathway lookup. Not integrable from current data.

### Routing verdicts
- **Decision 1 (categorical-vs-property): ONTOLOGY (3A).** A large shared,
  hierarchical vocabulary (InterPro entries) with cross-gene grouping.
- **Scored edge:** matches carry `evalue`/`score`/coordinates → these ride as
  **edge properties** on `Gene_has_interpro_entry` (the Changes_expression_of
  precedent, as with psortb/signalp).
- **Multiple calls per gene** → one edge per (gene, entry).
- **Decision 2 (functional-vs-structural): FUNCTIONAL** — InterPro entries are
  domain/family evidence of *what the gene does*. Folds into Gene routing
  (`annotation_types`). See §5 for the bounded fold-in decision.

### Agreed scope decisions (user, 2026-07-26)
1. **Node scope = InterPro entries only.** The ~31% unintegrated member-DB
   signatures are dropped; they are overwhelmingly redundant structural
   assignments (Gene3D/FunFam/SUPERFAMILY duplicate the integrated
   HOMOLOGOUS_SUPERFAMILY entries) or non-family predictors (MobiDB-Lite
   disorder, Coils). Net loss: only ~10–12 genes/strain get no IPS edge, all of
   which retain eggNOG/UniProt/Pfam annotations. Member DBs are preserved as the
   `libraries` list on each edge.
2. **Hierarchy = included now.** Requires an InterPro reference download
   (`entry.list` + `ParentChildTreeFile.txt`) — not in calls.json. Builds
   `Interpro_entry_is_a_interpro_entry` edges + real `level` depths.
3. **GO/pathways = deferred.** Empty in artifacts; GO already comes from
   eggNOG+UniProt. Follow-up: re-run Phase-1 with `--goterms` if wanted.

This design **supersedes** the Phase-1 doc's "Phase 2 sketch" (which left node
granularity, hierarchy, and merge-vs-ontology open).

## 1. Target graph objects

### Node: `InterproEntry` (Neo4j `InterproEntry`)
- **id:** `interpro:IPR003686` (`interpro` is a registered bioregistry prefix →
  colon form via `normalize_curie('interpro:…') or 'interpro_…'`).
- **props:** `name` (entry description), `interpro_id` (`IPR*`),
  `interpro_type` (FAMILY/DOMAIN/…), `level` (int, depth in the parent/child
  tree; 0 = root/no-parent). No `level_kind` (InterPro depth is not a named
  tier system; `interpro_type` carries the categorical kind instead).
- **Pruned to observed** entries (any `IPR*` referenced by a gene) **+ their
  ancestors** (so hierarchy edges never dangle), mirroring TCDB/BRITE pruning.

### Edge: `Gene_has_interpro_entry` (Gene → InterproEntry)
One per (gene, entry). **Scored ontology edge** (Changes_expression_of
property-block precedent). Properties, aggregated across that gene's matches to
the entry:
- `start` (int) = min start, `end` (int) = max end (domain envelope)
- `evalue` (float, nullable — best/min non-null across matches)
- `score` (float, nullable — max non-null across matches)
- `libraries` (str[]) = member DBs supporting this entry for this gene
- `match_count` (int) = number of match×location rows folded into this edge

**No e-value cutoff — by design.** InterProScan is an orchestrator over ~15
member DBs, each applying its **own curated significance threshold** before
reporting (HMM DBs like Pfam/NCBIfam use per-model *gathering (GA) cutoffs*;
profile/pattern DBs use per-profile score cutoffs; patterns are deterministic).
So every match in calls.json is already significant by that DB's principled
standard — there is no sub-threshold noise to filter. Moreover **~42% of matches
carry no e-value at all** (MED4: 5,872/14,025 — HAMAP/PROSITE/PANTHER/PRINTS/
MobiDB/Coils report a score or nothing), and e-values are **not comparable across
member DBs**. An extra global e-value filter would therefore (a) be less
principled than the per-DB GA thresholds already applied and (b) wrongly drop
every non-HMM match. `evalue`/`score` are kept as **nullable evidence/ranking
payload only** — never a filter. Low-information predictors (MobiDB-Lite, Coils)
are excluded *structurally* by the entries-only design (never integrated to IPR).

### Edge: `Interpro_entry_is_a_interpro_entry` (child → parent)
From `ParentChildTreeFile.txt`, pruned to kept entries.

### Edge: `Pfam_in_interpro_entry` (Pfam → InterproEntry)  — overlap bridge
Member-of link connecting the existing eggNOG-sourced Pfam layer to the new
InterPro layer (they overlap on Pfam as a member DB — MED4: ~88% of Pfam
accessions shared, 99% per-gene agreement — but neither replaces the other:
InterProScan is a method-independent cross-check and each source contributes
net-new accessions/genes). **Not a merge** — a navigable link.
- Derived free from the calls.json PFAM matches: `signature_accession` (PF*) →
  `interpro_accession` (IPR*).
- **Dangling-proof pruning:** emit only where the source Pfam node exists (PF ∈
  the KG's global Pfam node set) AND the target InterproEntry is kept. The Pfam
  node set is injected into `MultiInterproAnnotationAdapter` (BRITE-`known_ko_ids`
  precedent) so the bridge never references a non-existent Pfam node (the ~162
  InterProScan-only Pfam accessions per strain have no eggNOG-built node).

## 2. Data flow

```
<strain>.interproscan.calls.json  (WP → matches; committed Phase-1 artifact)
        │
        ├── MERGE (Step 2): interproscan logical source → merged field
        │     interpro_entries: list[str]   (join protein_id = WP_)
        │     → gives contributing_sources='interproscan' + DataSource node
        │       + Gene routing fold-in (annotation_types, interpro_entry_count)
        │
interpro_reference.json (entry.list + ParentChildTree → {id:{name,type,parent,level}})
   committed cache: cache/data/interpro/interpro_reference.json  (built by prepare_data step 9)
        │
        └── ADAPTER (Step 3A): interpro_adapter.py
              • per-strain: reads merged JSON (locus_tag→{protein_id, interpro_entries})
                + the strain calls.json (WP→matches) → Gene_has_interpro_entry edges
                with aggregated coords/evalue/score/libraries
              • multi (orchestrator): owns InterproEntry nodes (observed+ancestors)
                + Interpro_entry_is_a_interpro_entry edges, using the reference
```

Rationale for adapter-reads-calls.json: the rich per-match evidence
(coordinates/e-value/score) is heavy for `gene_annotations_merged.json`; the
merge only needs the light `interpro_entries` list for source bookkeeping, and
the adapter reads the Phase-1 artifact directly — exactly like `tcdb_adapter`
reads `tcdb_pruned.json` and `metabolism_adapter` reads `kegg_data.json`.

## 3. Reference download — wired into `prepare_data.sh` as **step 9**

All data download flows through `prepare_data.sh` (no adapter-time network). Two
pieces:
- **InterProScan calls.json** (the per-strain match data) — Phase-1, produced by
  the `/interproscan-run` skill (a heavy Docker tool, exactly like
  eggnog/psortb/signalp; **not** a prepare_data step). Already committed for 42
  strains.
- **InterPro reference** (names/types/hierarchy) — a lightweight global download
  → **new `prepare_data.sh` step 9**.

**Step 9** (`multiomics_kg/download/build_interpro_reference.py`) downloads from
InterPro `current_release` FTP:
- `entry.list` → `{IPR*: (type, name)}`
- `ParentChildTreeFile.txt` → parent map + is-a depth
and writes `cache/data/interpro/interpro_reference.json`
(`{id: {name, type, parent, level}}`, indented for git-friendly diffs).
**Committed** so the Docker build runs offline (Pfam/KEGG precedent). Pure parse
lives in `multiomics_kg/utils/interpro_reference.py` (unit-testable, no
filesystem); the step-9 module handles download + I/O + caching. `--force`
rebuilds from cached raw; `--refetch-raw` re-pulls the FTP files (only on an
InterPro release). Independent of steps 0–8; add `9` to the default `STEPS` list.

## 4. Post-import rollups + indexes

Both `post-import.sh` and `post-import.cypher` (byte-identical Cypher):
- On `InterproEntry`: `gene_count` (distinct genes via `Gene_has_interpro_entry`),
  `organism_count`, `member_count` (recursive subtree gene count via the cazy
  `*0..` walk — hierarchical). Optional per-entry `rank_by_score` on the edge
  (mirrors expression `rank_by_effect`) — **deferred** unless cheap.
- `is_promiscuous` (bool) — flag ultra-common entries (high `gene_count`,
  modelled on tcdb's `is_promiscuous`) so a `(type, level)`-stratified ORA can
  down-weight/flag them within a band. Cheap (gene_count is computed anyway).
  See §8.
- Indexes: scalar on `level` + `interpro_id` + `interpro_type`; full-text
  `interproEntryFullText` on `name`. `interpro_type` **must** be indexed — it is
  the primary ORA stratification key (§8).

## 5. Cross-node rollups (functional, bounded)

- **Fold in (minimal):** add `interpro` to Gene `annotation_types`; add a
  `Gene.interpro_entry_count` routing signal (distinct `Gene_has_interpro_entry`
  edges).
- **Deferred (recorded, not done in first pass):**
  - `informative_annotation_types` — needs an InterPro term-informativeness
    filter (which entries are catch-all); design separately.
  - `annotation_quality` 9th bucket — InterPro is highly redundant with the
    existing `pfam`/`go`/`ec` buckets, so it would lift few genes' quality while
    perturbing `annotation_state`, the bucket-count test, and `genes_by_function`
    `min_quality`. Not worth the blast radius in pass 1.
- `OrganismTaxon`/`Publication`: nothing.

## 6. Validation gate (Step 5)

`/omics-edge-snapshot` before → `docker compose up -d --build` → post-import →
`/omics-edge-snapshot` after (expression edges unchanged; new edges appear) →
kg-validity assertions (InterproEntry count > 0; no `Unknown`-style sentinel;
edge `evalue`/`score` ranges; no orphan gene→entry edges; hierarchy acyclic;
rollup sanity) → `pytest -m kg` → regenerate `snapshot_data.json`.

## 8. InterPro & ORA enrichment (design note — informs modeling, not built here)

The explorer's Over-Representation Analysis (ORA) stratifies terms **by `level`**
to keep each band at comparable granularity and reduce parent/child dependence.
That works for GO/KEGG/tcdb/cazy because there `level` (is-a depth) is a good
granularity proxy. **InterPro inverts this:**

- The InterPro hierarchy is **shallow and sparse** — most entries have no curated
  parent and sit at `level 0`, so `level 0` mixes specific FAMILY entries with
  promiscuous DOMAIN entries and broad HOMOLOGOUS_SUPERFAMILY entries. `level`
  alone barely partitions InterPro and does **not** separate broad from specific.
- Breadth is carried by **`interpro_type`** (FAMILY vs DOMAIN vs
  HOMOLOGOUS_SUPERFAMILY vs SITE/REPEAT/PTM), orthogonal to depth.

**Verdict:** ORA over InterPro is valid and complementary to pathway/GO ORA
(it asks "which domains/families are over-represented?"), **but only if it
stratifies by `(interpro_type, level)`** — with `interpro_type` as the **primary**
axis (the reverse of every other KG ontology, where `level` dominates) and
`level` secondary (separates parent/child generations where the sparse hierarchy
exists; degenerate-but-harmless elsewhere). Consequences for this build:
- Populate `level` as **true is-a depth** from ParentChildTreeFile (0 = parentless
  root, +1 per descent). Do **not** encode breadth into `level`; keep breadth in
  `interpro_type` (preserves the KG-wide `level` invariant).
- Keep `interpro_type` a **first-class, indexed, stratifiable** property → the
  one-label + type-property modeling (not per-type labels) is what makes this work.
- `is_promiscuous` (§4) lets ORA down-weight ultra-common domains within a band.

The ORA itself lives in the explorer/MCP and is **out of scope** here; this note
exists so the node/index modeling supports well-formed ORA later.

## 7. Out of scope / follow-ups

- Member-DB signature nodes (unintegrated hits) — dropped by design.
- **InterPro → GO / pathway cross-links** — empty in current artifacts; needs a
  Phase-1 re-run with `--goterms --pathways`. This is the *enrichment-enabling*
  follow-up: InterPro's value for **pathway/GO ORA** flows entirely through these
  xrefs (most valuable for heterotrophs, where eggNOG/UniProt GO is thin and
  InterPro is method-independent). This pass is deliberately enrichment-neutral.
- **Domain-composition relationship** — InterPro also has a cross-type
  "contains / found-in" relationship (family *contains* domain) **not** in
  ParentChildTreeFile (is-a only). Clean future add:
  `Interpro_entry_contains_interpro_entry`.
- `informative_annotation_types` + `annotation_quality` bucket — §5.
- MCP/explorer surfacing of InterproEntry (incl. `(type, level)`-stratified ORA) —
  separate, out-of-scope task.
