# InterPro two-layer integration — richer functional annotation + ORA — design

**Date:** 2026-08-10
**Status:** design proposed; not implemented
**Predecessors:**
- `docs/superpowers/specs/2026-07-22-interproscan-domains-design.md` (Phase 1 — calls.json artifacts)
- `docs/superpowers/specs/2026-07-26-interproscan-kg-integration-design.md` (Phase 2 — `InterproEntry` ontology, **live in code**)
- `docs/superpowers/specs/2026-08-07-interpro-merge-source-design.md` (gene-field propagation + union attribution — **this spec supersedes/extends it**)

## 0. Why this spec exists

The 2026-08-07 spec proposed folding InterPro's entry-level cross-references (GO,
EC, CAZy, Pfam) into the **gene-annotation merge fields** so they enrich genes.
That is correct but **incomplete**: it treats cross-references that are unsafe to
push onto genes (multi-EC catalytic *domains*, fold-level CAZy) as *noise to
drop*.

This spec adds the missing half. Every InterPro cross-reference is real curated
data. Whether it belongs on a **gene** depends on the entry *type* (§2); but even
when it doesn't, it still belongs on the **entry**. So we model **two
complementary layers** and route each cross-reference to whichever fits — nothing
is discarded.

It also adopts **edge-level source attribution** (values like `eggnog`,
`uniprot`, `cyanorak`, `interpro`) on the gene→ontology edges, so provenance is
queryable in the graph and ORA can be filtered/weighted by source. That requires
fixing the union-field provenance defect the 2026-08-07 spec identified (§6).

### Goals
- **Richer functional annotation** — more genes with GO/EC/CAZy/Pfam, especially
  heterotrophs where UniProt/CyanoRak are thin.
- **Support ORA** (over-representation analysis) — both direct 1-hop
  `gene → term` and 2-hop `gene → InterproEntry → term`.
- **Provenance** — every functional edge says which source(s) asserted it.

### Explicitly out of scope
- MCP/explorer implementation. This is a **KG-side** spec. §7 records the exact
  MCP changes the KG enables so the MCP work can follow, but builds none of it.
- Residue-level sites / unintegrated signatures / PANTHER subfamily GO — deferred
  (§10.1), mentioned only.
- NCBIfam `hmm_PGAP.tsv` route and the TIGRFAM→TigrRole heterotroph-role bridge —
  deferred (§10.2).
- MetaCyc pathway layer — deferred (§10.3).

## 1. Decisions locked (2026-08-10, with the user)

1. **Two layers, routed per cross-reference** (§2) — this is the core new idea.
2. **Edge-level source attribution** (TCDB `sources: str[]` precedent) on the
   gene→ontology edges (§5).
3. **Propagation filters are type-aware, per-field, and noise-gated** (§3/§5.1):
   Pfam = direct HMM hits; GO = FAMILY+DOMAIN (**fold excluded**); EC =
   **FAMILY + single-EC** (multi-EC entries are noise > data → Layer A only);
   CAZy = FAMILY+DOMAIN, fold excluded.
4. **Uncertainty is explicit, not implied** (§2.1, §5.3): every gene→ontology edge
   carries `evidence` (`curated | family_inferred | domain_inferred`) **and** a
   small advisory `evidence_score`, so an LLM/MCP can tell a curated fact from a
   domain guess without decoding `sources`. The label alone must never be read as
   a confident assertion.
5. **Noise-gate principle** (user, 2026-08-10): *do not create a gene edge when the
   propagating entry is ambiguous enough that the individual edge is more likely
   wrong than right* (the multi-EC catalytic-domain case). Such links live in
   Layer A only, marked (§4).
6. **Layer A trimmed to marked EC/CAZy** (§4): drop Layer A GO (redundant — is-a
   ancestors add only 5 entries, so Layer B already covers those genes). Keep
   EC/CAZy with `interpro_type` + an `ambiguous` flag, **weak verbs** (`related_to`,
   never `has`/`enables`/`catalyzes`), documented as a **recall-biased router**
   (`Publication_discusses` / TCDB-bridge precedent) — not an annotation.
7. **NCBIfam / hmm_PGAP + TigrRole bridge → deferred** (§10.2).
8. **Residue-level sites / PANTHER subfamily → deferred**, mentioned only (§10.1).
9. **Phased delivery with validation checkpoints** (§8).
10. **Node injection over pruning** for the Layer-A edge targets, so ontology-level
    ORA works (TCDB bridge precedent, §4.3).

## 2. The two-layer model

InterPro attaches cross-references (GO/EC/CAZy/pathway) to an **entry**.
"Propagation" = copying an entry's cross-references onto the genes that matched
it. Whether that is biologically valid depends on the entry **type**:

- `FAMILY` — "this protein **is an** X" → licenses a functional assignment.
- `DOMAIN` — "this protein **contains** an X-like region" — one piece of the
  protein.
- `HOMOLOGOUS_SUPERFAMILY` — **fold-level**; shared 3D shape, unrelated function.

Given that, each cross-reference is routed to one or both layers:

| Xref | **Layer A** — ontology→ontology `InterproEntry → X` (compact, marked, router) | **Layer B** — gene→ontology (propagated, filtered, ORA-ready) |
|---|---|---|
| **GO**  | **none** — dropped as redundant (§4.1); genes get GO via Layer B | fold into `go_terms`, **FAMILY+DOMAIN (fold excluded)** |
| **EC**  | `Interpro_entry_related_to_ec_number` — **home for the multi-EC / DOMAIN ECs** Layer B refuses (kept, marked `ambiguous`, router) | fold into `ec_numbers`, **FAMILY + single-EC** (noise-gated) |
| **CAZy**| `Interpro_entry_related_to_cazy_family` — **home for fold-level CAZy** | fold into `cazy_ids`, **FAMILY+DOMAIN (fold excluded)** |
| **Pfam**| `Pfam_in_interpro_entry` — **already exists** ✓ (2026-07-26) | direct PFAM HMM hits → `pfam_ids` |
| **Pathway (MetaCyc)** | natural Layer-A home once a Pathway node exists — **deferred (§10.3)** | (no clean gene target) |

**Reading the table:** Layer B is the subset safe enough to enrich genes directly,
and each such edge is stamped with its `evidence` (§5.3). Layer A keeps the
*un-propagatable* links (multi-EC domains, fold CAZy) — the ones the noise-gate
excludes from genes — as an explicitly-marked, neutrally-named **router**. A gene
with a GH13 *catalytic domain* gets **no** EC on its gene edges; an analyst who
wants the 22 domain-level candidates traverses `gene → entry → EC` through the
`related_to` edges, which announce their looseness in the label itself.

**Where each layer earns its place:**
- Layer A exists **only for EC and CAZy**, and only for the links Layer B must
  refuse (multi-EC / DOMAIN ECs, fold CAZy). It is compact, uses deliberately
  weak verbs (`related_to`), and is documented as recall-biased — like
  `Publication_discusses_*` and the TCDB bridges. Layer A GO was dropped: is-a
  ancestors add only 5 entries, so it duplicated Layer B for no gain.
- Layer B makes gene annotations richer, single-hop, ORA-ready, and folds into
  `Gene.annotation_types` / `informative_annotation_types` / `annotation_quality`
  — with `evidence` preserving the curated-vs-inferred distinction.

### 2.1 Uncertainty is carried in properties + labels, never assumed from the label

An LLM reads the edge label as an assertion. Two mechanisms keep that honest:
- **Layer B (existing labels kept):** `Gene_catalyzes_ec_number` etc. stay (renaming
  live edge types is out of scope — blast radius). Confidence rides in **`evidence`**
  (`curated | family_inferred | domain_inferred`) + **`evidence_score`** (§5.3), which
  the MCP filters on. `sources` carries provenance; `evidence` carries strength.
- **Layer A (new labels):** deliberately weak verbs — **`related_to`**, never
  `has`/`enables`/`catalyzes` — plus an `ambiguous` edge flag and the entry's
  `interpro_type`. The label alone signals "cross-reference, not claim."

## 3. Current state (verified 2026-08-10)

| Component | State |
|---|---|
| `InterproEntry` nodes + `Gene_has_interpro_entry` + `Interpro_entry_is_a_interpro_entry` + `Pfam_in_interpro_entry` | ✅ in code (2026-07-26). **Not in the currently-deployed graph** — the live MCP graph predates it; a fresh Docker rebuild is required for any validation. |
| `interpro_reference.json` (step 9) | has `go_terms` (14,803) + `pathways` (5,091); **NO `ec_numbers`, NO `cazy`** |
| Per-strain `entry_xrefs.json` | GO + pathways populated (re-run `b8172a29`) |
| `_resolve_union` source attribution | **defect present** — never reads `track_source`; `go_terms_source` is a silent no-op |
| `interproscan` merge source | passes through only `interpro_entries` (no GO/EC/CAZy/Pfam contribution) |
| Gene→GO/EC/Pfam/CAZy edges | **empty properties** — no `sources` attribution |

So: Layer-A Pfam exists; **everything else in this spec is to-build.**

## 4. Layer A — ontology→ontology edges

### 4.1 New edge types — EC and CAZy only, deliberately weak verbs

Two new edge types (GO dropped — redundant, §2). Both use `related_to`, not
`has`/`enables`/`catalyzes`, so the label itself reads as a loose cross-reference:
- `Interpro_entry_related_to_ec_number` (InterproEntry → EcNumber)
- `Interpro_entry_related_to_cazy_family` (InterproEntry → CazyFamily)
- `Pfam_in_interpro_entry` — **already exists**, no change.

**Properties:**
- `source_db` (str — `interpro.xml`) — xref origin/provenance.
- `ambiguous` (bool) — `true` when the entry maps to multiple mutually-exclusive
  terms of this kind (the multi-EC catalytic domain, one entry → 22 ECs) OR the
  entry is a non-FAMILY type. This is the *"more noise than data"* marker: an
  `ambiguous=true` edge is one candidate among several, not a specific claim.
- `interpro_type` is read from the target-adjacent `InterproEntry` node, not
  duplicated on the edge.

**Documentation requirement:** the schema comment + CLAUDE.md must state, verbatim
in spirit, what the TCDB bridge says — *these are a recall-biased router; read
outward from a gene's known family they corroborate, read backward (carries EC →
therefore is that enzyme) they are low-precision; never use them to assign a gene
its function.* That is what `Gene_catalyzes_ec_number` (Layer B) is for.

### 4.2 Who emits them

`multiomics_kg/adapters/interpro_adapter.py` (`MultiInterproAnnotationAdapter`),
reading the extended `interpro_reference.json` (§6.1). It already owns
`InterproEntry` nodes + hierarchy + the Pfam bridge, so the new cross-links belong
there. Emit one edge per (entry, EC|CAZy xref) for every kept `InterproEntry` node
that carries one, computing `ambiguous` from the entry's term-cardinality + type.

### 4.3 Dangling-proofing — inject, don't prune

Layer-A targets include EC/CAZy terms **no gene otherwise carries** (DOMAIN-level
EC, fold-level CAZy). For ontology-level ORA those target nodes must exist. Follow
the **injection** precedent (`create_knowledge_graph.py` already injects
`pfam_node_ids` + `go_terms` into the interpro adapter, and the TCDB bridges inject
their targets):

- Feed the InterPro-referenced EC/CAZy ids into the **seed sets** of the
  respective node adapters (`ec_adapter`, `cazy_adapter`) so those adapters build
  the term nodes **and their hierarchy/ancestry** properly (BRITE `known_ko_ids`
  precedent). The interpro adapter then emits only the cross-link edges,
  dangling-proof.
- Overlap note: EC/CAZy that Layer B propagates to genes are already seeded via
  gene annotations; injection is only strictly needed for the **un-propagated**
  xrefs (multi-EC / DOMAIN EC, fold CAZy).

## 5. Layer B — gene→ontology propagation + edge-level source attribution

### 5.1 Post-merge enrichment (mirror `enrich_pfam_fields`)

Add an `enrich_interpro_fields()`-style post-merge function in
`multiomics_kg/download/build_gene_annotations.py`. For each gene, using the
gene's matched entries (`interpro_entries` in the merged JSON) + the extended
reference (`interpro_reference.json`) + per-strain `entry_xrefs.json`:

| Field | Contribution | Filter (noise-gated) | `evidence` stamped |
|---|---|---|---|
| `pfam_ids` | direct PFAM signature hits from `calls.json` | none (raw signature) | `signature` |
| `go_terms` | entry GO (per-protein via `entry_xrefs.json`, fallback reference) | **FAMILY+DOMAIN; exclude fold** | `family_inferred` / `domain_inferred` |
| `ec_numbers` | entry EC (reference only) | **FAMILY + single-EC**; `normalize_ec`; bare 3-level `3.4.21` → `3.4.21.-` | `family_inferred` |
| `cazy_ids` | entry CAZy (reference only) | **FAMILY+DOMAIN; exclude fold** | `family_inferred` / `domain_inferred` |
| `alternate_functional_descriptions` | entry names, `[interpro]`-prefixed (mirrors `[pfam]`) | — | — |

**Noise gate (the "more noise than data" rule).** An entry only donates a term to
a gene when that term is a *specific* claim for the entry:
- **EC → FAMILY + single-EC.** A FAMILY entry carrying one EC is a definitive
  "this *is* enzyme X". A FAMILY (or DOMAIN) entry carrying *N* mutually-exclusive
  ECs would make each gene edge ≈1/N precise — below break-even for N≥2 — so those
  do **not** become gene edges; they live in Layer A as `Interpro_entry_related_to_ec_number`
  with `ambiguous=true`. (2026-08-07 §3.1 measured FAMILY+single-EC at the cleanest
  ratio, 0.5.)
- **GO / CAZy → FAMILY+DOMAIN, fold excluded.** GO terms and CAZy families are not
  mutually exclusive (a protein legitimately has several), so cardinality is not
  noise there — only fold-level (shape-only) entries are excluded.

Entry `type` + term-cardinality come from the reference (`type` + `ec_numbers` /
`cazy_ids` keys). This is why the filter must be a post-merge function, not
declarative config: it is type-aware, cardinality-aware, and per-field.

### 5.2 Union source attribution (blocking prerequisite — 2026-08-07 §4)

Edge-level `sources` (§5.3) is only possible if the merge records, per token,
which source(s) contributed it. Today `_resolve_union` ignores `track_source`.

**Shape:** per-token map (2026-08-07 Option B) —
`go_terms_source: {"GO:0003677": ["uniprot","interpro"], ...}`, likewise
`ec_numbers_source`, `pfam_ids_source`, `cazy_ids_source`. This is the only shape
that survives filtering/reordering and represents multi-source tokens honestly
(the common case: 146,576 InterPro GO pairs already come from another source).

**Work:**
- Implement per-token source capture in `_resolve_union`
  (`multiomics_kg/download/build_gene_annotations.py`).
- Rewrite `_compute_contributing_sources` / `_has_source_label` to read the map
  shape (they assume scalar `*_source` fields today).
- Decide compat for the 3 existing scalar `*_source` fields (`product_source`,
  `gene_name_source`, `function_description_source`) — keep scalar; the map shape
  is additive for union fields only.

### 5.3 Edge-level provenance + confidence on gene→ontology edges

To these edge types — `Gene_involved_in_biological_process` /
`_located_in_cellular_component` / `_enables_molecular_function`,
`Gene_catalyzes_ec_number`, `Gene_has_pfam`, `Gene_has_cazy_family` — add:

- **`sources: str[]`** (sorted; `ncbi|cyanorak|eggnog|uniprot|interpro`) —
  *provenance*: who asserted it. From the per-token `*_source` map (§5.2).
- **`evidence: str`** — *inference strength*, the single most important field for an
  LLM: `curated` (UniProt/CyanoRak/eggNOG), `signature` (direct Pfam HMM hit),
  `family_inferred`, `domain_inferred`. When multiple sources agree, the strongest
  wins (`curated` > `signature` > `family_inferred` > `domain_inferred`).
- **`evidence_score: int` (0–3, advisory — like `tcdb_evidence_score`)**: `+1`
  corroborated by ≥2 independent sources (Pfam-via-eggNOG and Pfam-via-InterPro are
  **not** independent — §9.2); `+1` from a curated source or direct signature;
  `+1` from a FAMILY entry (vs DOMAIN). Advisory only — nothing filters on it; it
  gives the MCP a ready sort key.

The label stays a coarse relation; `evidence` + `evidence_score` carry the
confidence the label overclaims. Adapters: `functional_annotation_adapter.py`
(GO/EC/Pfam), `cazy_adapter.py` — currently emit `{}`. Follow the
`tcdb_adapter.py` `sources` pattern: one edge per (gene, term),
`sources = sorted(set(...))`, then derive `evidence` / `evidence_score`.

`schema_config.yaml`: add `sources: str[]`, `evidence: str`, `evidence_score: int`
to each edge above.

## 6. Reference + data-flow work items

### 6.1 Extend step 9 to parse EC + CAZy
`multiomics_kg/utils/interpro_reference.py` + `build_interpro_reference.py`:
extend the streaming `interpro.xml.gz` db_xref parser (`_XML_DB_XREF_RE`) to also
keep `db="EC"` and `db="CAZY"` (verify exact tokens on refetch — the raw XML is
gitignored, not in-checkout). Emit sparse `ec_numbers` / `cazy_ids` keys per
entry alongside `type`. Spec-§2 measured ~5,100 entries with EC, 113 with CAZy.

### 6.2 Lazy reference load inside step 2
Mirror `load_pfam_data`: `build_gene_annotations.main()` calls
`build_interpro_reference.build()` (idempotent — returns cached dict unless
`--force`). Step 9 remains the explicit refresh entry point. No `prepare_data.sh`
renumbering.

### 6.3 Extend the `interproscan` source in config
`config/gene_annotations_config.yaml`: widen the `interproscan` source to declare
it contributes `go_terms` / `ec_numbers` / `cazy_ids` / `pfam_ids` /
`alternate_functional_descriptions` (actual filter logic lives in the §5.1
post-merge function, not declarative config).

## 7. ORA support + MCP follow-up (spec only — do not build)

The KG changes enable two ORA paths; the MCP engine work is a separate follow-up:

1. **Direct 1-hop** `gene → term` — from Layer-B edges. Existing MCP tools
   (`pathway_enrichment`, `genes_by_ontology`, `ontology_landscape`) already count
   these; they gain InterPro-derived coverage automatically once Layer B lands.
   With `sources` on edges, these tools **can** expose a source filter
   (curated-only vs including domain-inferred).
2. **2-hop** `gene → InterproEntry → term` — from Layer-A edges. Lets ORA use even
   un-propagated links (DOMAIN EC, fold CAZy) when an analyst opts in. No MCP tool
   does this today.

**MCP changes the KG enables (for the follow-up spec):**
- New optional `source_filter` / `evidence_filter` (e.g. `curated` only vs
  include inferred) on ontology-edge tools, reading the edge `sources` / `evidence`.
- InterPro ORA must **stratify by `(interpro_type, level)`** — `interpro_type`
  primary (breadth), `level` secondary (2026-07-26 §8). `InterproEntry.is_promiscuous`
  already flags ubiquitous entries for down-weighting.
- Expose the 2-hop `gene → entry → EC|CAZy` (`related_to`, `ambiguous`) traversal as
  an explicit opt-in "candidate/router" mode — never in default annotation output.
- Surface edge `sources` + `evidence` in `gene_ontology_terms` / `gene_overview`.

## 8. Phased delivery (checkpoints between phases)

**Phase 1 — foundation (independent, parallelizable internally).**
- 1a. Extend step-9 reference: EC + CAZy parsing (§6.1) + `normalize_ec` + bare-EC
  normalization. Unit tests on XML fixtures.
- 1b. Implement union `track_source` + rewrite `_compute_contributing_sources`
  (§5.2). Unit tests.
- **Checkpoint:** `pytest -m "not slow and not kg"` clean; merged-JSON for MED4
  shows per-token `*_source` maps + `contributing_sources` includes `interproscan`;
  reference cache now carries `ec_numbers`/`cazy_ids`.

**Phase 2 — Layer B (gene propagation + edge evidence). Depends on Phase 1.**
- 2a. `enrich_interpro_fields()` post-merge with §5.1 noise-gated filters.
- 2b. `sources` + `evidence` + `evidence_score` on Gene→GO/EC/Pfam/CAZy edges +
  adapters + schema (§5.3).
- 2c. Rebuild step 2 across 42 strains; Docker KG rebuild; post-import; indexes.
- **Checkpoint:** `/omics-edge-snapshot` before/after (expression edges
  unchanged); per-field net-new within tolerance of §11 (`pfam_ids` +13,977;
  `ec_numbers` +10,750 FAMILY+single-EC; `cazy_ids` +648; `go_terms` fold-excluded,
  record the measured value); every new edge carries `evidence`; `pytest -m kg`;
  regenerate `snapshot_data.json`.

**Phase 3 — Layer A (EC/CAZy cross-links). Depends on Phase 1 (reference EC/CAZy).**
- 3a. `Interpro_entry_related_to_{ec_number,cazy_family}` in `interpro_adapter.py`
  (§4.1) with `ambiguous` + `source_db`; node injection (§4.3).
- 3b. `schema_config.yaml` edge definitions + the router documentation (§4.1);
  post-import rollups (optional `InterproEntry.ec_count` / `cazy_count`); indexes.
- **Checkpoint:** rebuild; assert no dangling Layer-A edges; a multi-EC domain
  entry is reachable via `gene → entry → EC` (`related_to`, `ambiguous=true`) but
  its ECs are **absent** from the gene's `Gene_catalyzes_ec_number` edges;
  `pytest -m kg`.

**Phase 4 — docs + handoff.**
- `docs/kg-changes/interpro-two-layer.md` (what-changed / release notes).
- CLAUDE.md: new edges, `sources` property, key facts, `informative_annotation_types`
  / `annotation_quality` decision (§9.2).
- Register deferred items (§10). Record the MCP follow-up (§7).

Phases 2 and 3 are independent given Phase 1 and may run in parallel.

## 9. Validation gates (all phases)

- Unit tests for new parsers + filters (fixtures, not live downloads).
- `pytest -m "not slow and not kg"` — baseline 2,127 passed (fresh worktree
  without `.env` shows unrelated omics failures — missing `OPENAI_API_KEY`).
- `/omics-edge-snapshot` before/after — expression edges must be unchanged.
- `pytest -m kg` after each rebuild; regenerate `snapshot_data.json`.
- Per-field delta assertions vs 2026-08-07 §2 (a large divergence = a filter
  regressed).
- Assert `contributing_sources` includes `interproscan` for genes with
  InterPro-derived tokens (observable proof §5.2 landed).
- Assert Gene→ontology edges carry non-empty `sources`; multi-source tokens show
  ≥2 entries.

## 9.1 `annotation_quality` / `informative_annotation_types` semantics — decide in Phase 4
~7,000 genes gain a **first** GO term from family/domain inference (slightly fewer
after fold exclusion), shifting the `annotation_state` distribution. Is inferred
GO "informative"? Defensible yes, but the 8-bucket source list is hand-maintained
(`docs/superpowers/specs/2026-05-01-explorer-frictions-resolution-design.md`).
**Recommendation:** InterPro-inferred GO/EC/CAZy **do** count toward
`informative_annotation_types`, but the edge `evidence` property preserves the
curated-vs-inferred distinction for consumers that care. Confirm during Phase 4.

## 9.2 Circularity with Pfam (2026-08-07 §9.3)
eggNOG → Pfam; InterPro integrates Pfam. Adding InterPro's Pfam hits to `pfam_ids`
is **corroboration, not independent evidence** (0.10 ratio). The `sources` edge
property makes this explicit — a Pfam edge with `sources:['eggnog','interpro']`
must not be counted as two independent lines of evidence downstream.

## 10. Deferred (mentioned, not built)

### 10.1 Residue-level sites + unintegrated signatures + PANTHER subfamily GO
The only source is the raw `*.interproscan.raw.json` (~17 GB, gitignored, **not in
any clone**). Serves *"what does this mutation mean"* (active/binding sites with
exact residues) and sharper *"what is this hypothetical"* (PANTHER subfamily GO).
Does not fit the ontology pattern (a site is positional, on a gene×signature
pair). **⚠️ Perishability:** extraction is a cheap re-parse *while the raw files
survive*; once deleted, recovery is a ~27h re-scan. The user holds a copy on
another machine — if this is ever picked up, extract to committed artifacts
first. Own spec.

### 10.2 NCBIfam via `hmm_PGAP.tsv` + TigrRole heterotroph-role bridge
Highest-precision source measured (EC ratio 0.06) and the only functional-role
path for heterotrophs (Alteromonas/Shewanella/Pseudomonas/Ruegeria have no roles
today). Gated on the TIGRFAM→TigrRole archive hunt (2026-08-07 §9.1). Own spec.

### 10.3 MetaCyc pathway layer
InterPro carries 79,683 MetaCyc `PWY-*` xrefs; complementary to MNX
(reactions/compounds, not pathways). Needs a new Pathway node type + no consumer
yet. Natural Layer-A home (`Interpro_entry_in_pathway`) when wired.
**No InterPro→KEGG mapping exists** — KEGG pathways stay eggNOG-KO-derived.

## 11. Edge-count estimates

Grounded in the scan artifacts (42 strains, **12,994** observed InterPro entries;
+5 is-a ancestors → **12,999** kept — the hierarchy is nearly flat) plus the
2026-08-07 measured deltas where the cache doesn't yet hold EC/CAZy.

### New edges — Layer A (`InterproEntry → EC|CAZy`), after the trim

| Edge type | Estimate | Basis |
|---|---:|---|
| `Interpro_entry_related_to_ec_number` | **~16,000** | spec §3.2 (~5,100 entries carry EC; multi-EC domains, avg ~3.2/entry); most flagged `ambiguous=true` |
| `Interpro_entry_related_to_cazy_family` | **~120** | spec (113 entries carry CAZy) |
| ~~`Interpro_entry_*_go`~~ | **0 (dropped)** | redundant — ancestors add only 5 entries |
| `Pfam_in_interpro_entry` | 0 new | already exists (~5,236) |
| **Layer A total** | **~16,000** | |

### New edges — Layer B (gene→ontology, net-new only), after the noise gate

| Edge type | Net-new | Filter |
|---|---:|---|
| `Gene_has_pfam` | **+13,977** | direct HMM hits |
| `Gene_*_go` (3 labels) | **~+40,000 (est.)** | FAMILY+DOMAIN, **fold excluded** (down from +50,762 all-types; exact TBD Phase 2) |
| `Gene_catalyzes_ec_number` | **+10,750** | **FAMILY + single-EC** (down from +31,492 FAMILY-only) |
| `Gene_has_cazy_family` | **+648** | exclude fold |
| **Layer B total** | **~65,000** | |

**Grand total ≈ 81,000 new edges** (~16K Layer A + ~65K Layer B) — down from the
~123K first draft; a modest addition on top of the existing ~397K
`Gene_has_interpro_entry` edges. The trims (drop Layer A GO, EC noise-gate to
single-EC, exclude fold GO) removed ~42K lower-confidence edges.

### Edges modified (not new) — `sources` + `evidence`
`sources: str[]`, `evidence: str`, `evidence_score: int` are added to **every**
existing edge of the four gene→ontology types. Existing edges gaining an
`interpro` tag (2026-08-07 "corroborating"): `Gene_has_pfam` ~133K, `Gene_*_go`
~147K, `Gene_catalyzes_ec_number` ~30K, `Gene_has_cazy_family` ~800. All *other*
edges of these types (eggNOG/UniProt/CyanoRak, uncorroborated) also get the
properties — `sources` without `interpro`, `evidence=curated`.

### Why Layer A is cheap on purpose
The ~16K multi-EC/domain links, if propagated to genes, would splatter into
**117,542** gene→EC pairs (spec §3.1 "all entries") — the ~11× blow-up the
single-EC gate + Layer A together avoid.

**Caveats:** EC/CAZy Layer-A counts are from the 2026-08-07 measurement, not the
live cache (which lacks EC/CAZy until Phase 1) — Phase 1 yields exact numbers.
The fold-excluded GO count is an estimate; Phase 2 measures it. Layer-B net-new
predates the two most recently added strains, so expect a slight uptick.

## 12. Key files

| Purpose | Path |
|---|---|
| Reference builder / parser | `multiomics_kg/download/build_interpro_reference.py`, `multiomics_kg/utils/interpro_reference.py` |
| Gene-annotation merge | `multiomics_kg/download/build_gene_annotations.py` (`_resolve_union`, `_compute_contributing_sources`, new `enrich_interpro_fields`) |
| Merge config | `config/gene_annotations_config.yaml` (`interproscan` source) |
| InterPro adapter (Layer A) | `multiomics_kg/adapters/interpro_adapter.py` |
| Gene→ontology adapters (Layer B `sources`) | `multiomics_kg/adapters/functional_annotation_adapter.py`, `cazy_adapter.py` |
| Edge-attribution precedent | `multiomics_kg/adapters/tcdb_adapter.py` (`sources` merge) |
| Schema | `config/schema_config.yaml` |
| Node injection wiring | `create_knowledge_graph.py` |
| Post-import | `scripts/post-import.sh` (+ `.cypher` mirror) |
