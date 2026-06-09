# Plan: Publication "discusses" edges

**Scope:** KG repo only. Add narrative-level `Publication → {Gene, KEGG pathway}`
edges extracted from full paper PDFs.

**Goal:** Capture genes and pathways a paper *discusses in its narrative* — which
the supplementary DE tables do **not** capture (regulators named in the discussion,
genes mentioned but not differentially expressed, pathways referred to by name). The
edges act as a **router**: "which papers should I open about gene/pathway X?" → the
agent then feeds those whole PDFs to the LLM (the retrieval mode that empirically
beat chunked RAG). No chunking, no embeddings, no `Chunk` nodes, no vector index.

**Explicitly out of scope:** `get_pdf` and any new PDF-reference Publication
properties (`pdf_path`/`has_pdf`) — that whole content-delivery side is a separate
downstream task and is **not** touched here.

---

## Design decisions (settled)

- **No chunking / no RAG store.** Retrieval (full PDF to the LLM) is a separate
  downstream concern. The graph only holds `Discusses` edges as a discovery index.
- **Edges are recall-biased, not ground truth.** Since the agent re-reads the actual
  PDF, edge precision matters less than coverage. Unresolved mentions are dropped.
- **Entity *linking*, not open extraction.** Mentions resolve to **existing canonical
  nodes** (Gene via `gene_id_mapping`, KEGG pathway via existing KeggTerm pathway
  nodes). We never mint new entity nodes from text — that would create a shadow
  population to reconcile against the curated graph.
- **Pathways = KEGG pathway nodes only for v1** (`KeggTerm`, `level_kind='pathway'`).
  BRITE / GO deferred.
- **Deterministic build.** LLM extraction is non-deterministic, so it runs once at
  extract time and is cached as JSON (same pattern as `cluster_extraction_{key}.json`).
  The adapter reads the resolved cache deterministically — preserves the
  byte-identical-dump guarantee.

## What is NOT changing

- No changes to expression edges, Experiment nodes, or the supp-table pipeline.
- No open relation extraction (no "gene X regulates gene Y" text triples).
- **No `get_pdf`, no PDF hosting, no new PDF-reference Publication properties** — the
  entire content-delivery side is a separate downstream task.
- No BRITE / GO / Metabolite mention edges (possible v2).

---

## Pipeline: extract → resolve → adapter

Mirrors the gene-side pipeline (cluster-extraction skill → prepare-data step 4 →
omics_adapter). Three artifacts per paper, each stage deterministic except Stage 1.

### Stage 1 — Extraction (LLM, per paper, non-deterministic)
- New module `multiomics_kg/extraction/topics/` (`extract.py` + `extraction_utils.py`)
  **modelled on `extraction/cluster/`** and **sharing** its common utils: top-level
  `extraction/pdf_utils.py` (PDF loading), the OpenAI Responses API + pydantic-schema
  client pattern, and load/save-extraction json helpers (factor any cluster-local helper
  that's actually generic up to a shared location rather than copy it). Invoked
  `python -m multiomics_kg.extraction.topics.extract` with `--paper`, `--force`,
  `--dry-run`, `--report` like cluster.
- Exposed as a standalone Claude skill `.claude/skills/extract-discussed-topics/`
  wrapping that module (interactive per-paper run + review), modelled on the cluster flow.
- Reads `papermainpdf` — **full PDF to the LLM** — plus the paper's **candidate strain
  set** (union of `experiments[].organism` from the paperconfig), passed into the prompt
  so the model attributes each gene to one of the paper's strains.
- Per **gene** item: verbatim **`identifiers`** (any locus tag / ORF / gene id /
  protein id the paper prints — high-value: resolve 1:1 via Tier-1 `specific_lookup`,
  and the locus-tag prefix doubles as the strain fingerprint), **`gene_name`**,
  **`strain`** attribution (one of the paper strains | `all` | `unspecified`), and the
  single **`prominence`** signal (`central` | `peripheral`).
- Per **pathway** item: surface form, **KEGG pathway id when the model knows it**,
  `prominence`.
- (`mention_count` dropped — one signal is enough.) Capturing identifiers is an explicit
  requirement, not best-effort: when the paper prints IDs, they must be extracted.
- **Triage signal** `uncaptured_identifiers` (`none`|`some`|`many`) + note in metadata:
  the model self-reports how many printed identifiers it left out (e.g. abbreviated
  per-strain locus numbers in tables like `ntcA (0246*)`). Surfaced in `--report` and the
  run log so we can later tell a systematic gap (refine the prompt) from a paper-specific
  one. v1 deliberately accepts current capture quality rather than chasing abbreviated
  table numbers.
- Writes raw `publication_topics/topics.json` (per-paper subfolder; doi in metadata;
  legacy flat files auto-migrated). Re-run to regenerate.

**Stage 1 status: built + validated** on Aharonovich 2016 (multi-strain). Strain
attribution, full-tag identifier capture (CCRG-1/2 → MIT9313 `PMT*`/`P9313_*`),
prominence, and the triage signal (`many` on this table-heavy paper) all confirmed.

### Stage 2 — Resolution (deterministic, prepare-data)
New prepare-data step that reuses gene-resolution utils but does **not** overload
step 4's CSV flow (different input/output artifacts). Input: raw topics json + paper
organism. Output: `publication_topics/topics_resolved.json` + a `resolution_report.txt`
(per-paper stats, method breakdown, truncated-id count, and an `unresolved_reasons` tally:
`descriptive_only` | `truncated_id_only` | `lookup_miss`). Best-effort: gaps are reported,
not chased.
- **Gene branch (strain-aware — see "Multi-strain handling" below):** resolve via the
  per-strain `gene_id_mapping`, **identifiers first** (Tier-1 `specific_lookup`, 1:1),
  gene_name second; attach resolved `locus_tag` + `strain` + `resolution_method`.
  Unresolved kept with null `locus_tag` (like step 4) so the adapter drops them and they
  stay auditable.
- **Pathway branch (global — see "KEGG pathway resolution flow" below).**

**Stage 2 status: built + validated.** `multiomics_kg/download/resolve_paper_topics.py`
(reuses `get_genome_dir`/`load_mapping_v2`/`resolve_row`; KEGG lookup from
`kegg_data.json['pathways']`). On Aharonovich 2016: 16/22 gene records + 8/10 pathways
resolved; `ntcA` strain=`all` correctly fanned to PMM0246 (MED4) + PMT1831 (MIT9313) via
per-strain name resolution; `ko04110` (not in pruned set) correctly dropped (no dangling).
27 unit tests in `tests/test_resolve_paper_topics.py`. Unresolved cases are exactly the
abbreviated bare numbers flagged by `uncaptured_identifiers`.

### Stage 3 — Adapter (build time, deterministic)
- New `multiomics_kg/adapters/publication_topics_adapter.py` + `Multi*` wrapper
  (adapter pattern; `get_nodes` no-op, `get_edges`, `download_data` loads the resolved
  json). Iterates paperconfigs (like `MultiClusterAdapter`); no lookups at build time —
  reads resolved ids only.
- **Edge source = publication DOI from the paperconfig** (`publication.doi`), built with
  the same `add_prefix_to_id("doi", …)` + override precedence the omics adapter uses, so
  the source matches the existing Publication node id exactly.
- **Targets:** gene → `ncbigene:{locus_tag}`; pathway → `kegg.pathway:ko*` (`_pathway_node_id`).
- Emit `Publication_discusses_gene` / `Publication_discusses_kegg_pathway` for resolved
  entries only, props `prominence` + `evidence` (resolution info stays in the json, off
  the edge). Dedup by (publication, target), central-wins. `_clean_str` on every string.
- Register in `create_knowledge_graph.py` (Multi* wrapper only).

### KEGG pathway resolution flow (the "different flow")
Pathways have no per-organism mapping — resolution is **global** against the KG's
**pruned** KEGG pathway set, and the lookup IS that set so dangling edges are impossible:
- Build the lookup once from `kegg_data.json` (the same pruned pathway set the kegg
  adapter emits as `KeggTerm` `level_kind='pathway'` nodes): `{kegg_id → name}` and a
  normalized `{lower(name) → kegg_id}`.
- Resolve each mention **id-first** (LLM-supplied KEGG id, normalized across
  `ko`/`map`/bare-number variants), then **normalized-name match** fallback.
- Resolves ⇒ guaranteed real `KeggTerm` node (**no dangling edges by construction**);
  outside the set ⇒ dropped (pathway no marine bacterium has).
- Runs in Stage 2 alongside the gene branch, loading the global lookup once (not per
  organism).

### Multi-strain handling
Many papers study several strains; Gene nodes are strain-specific (locus tags are
per-strain), so each gene mention must resolve to a *specific* strain's Gene node.
- **Paper strain set** = union of `experiments[].organism` from the paperconfig.
  Passed into the extraction prompt (Stage 1) so the LLM attributes mentions.
- **Resolution rules (Stage 2):**
  - identifier present → strain implied by the locus-tag prefix; resolve 1:1 in that
    strain (prefix mismatch with the attributed strain = logged, identifier wins).
  - gene_name + specific strain → resolve in that strain's `gene_id_mapping`.
  - gene_name + `all`/`unspecified` → resolve in **each paper-strain** mapping; emit one
    edge per strain that resolves. **Bounded to the paper's strains, never all 23** — a
    genus-level "ntcA" in a MED4+MIT9313 paper yields 2 edges, not 23.
- Pathways are strain-agnostic (global `ko*`) — multi-strain does not affect them.
- **v2 option (noted, not built):** genus-level concepts could instead link
  `Publication → OrthologGroup` (the cross-strain bridge) rather than per-strain Gene
  edges — cleaner semantics, larger change. v1 stays with per-strain Gene edges.

### Schema
- `config/schema_config.yaml`: two new edge types following the
  `Gene_has_tcdb_family` naming convention:
  - `Publication_discusses_gene` (Publication → Gene), props: `prominence`, `evidence`
  - `Publication_discusses_kegg_pathway` (Publication → KeggTerm), same props
  - Edge carries the extraction `evidence` quote; resolution info (`resolved_from`,
    `resolution_method`) stays in the resolved json for audit only, NOT on the edge.

### Tests
- Unit: adapter resolves a known gene mention to the right locus tag; drops an
  unresolvable mention; resolves a KEGG pathway id and name to an existing node;
  sanitization.
- KG validity (`tests/kg_validity/`): edges exist, endpoints valid (no dangling
  targets), `salience` ∈ {focal, passing}.

### Docs
- CLAUDE.md: new edges + new extraction skill / resolution step / adapter facts.

---

## Acceptance criteria

- [x] `publication_topics.json` (raw) produced per paper (Stage 1 done).
- [x] `publication_topics_resolved.json` produced per paper; resolution deterministic and
      `--force`-rebuildable (Stage 2 done; full corpus 41/41).
- [x] `Publication_discusses_gene` + `Publication_discusses_kegg_pathway` edges built —
      BioCypher writes 1099 gene + 140 pathway edges (offline). **Zero-dangling check** is
      a `-m kg` test (`tests/kg_validity/test_discuss_edges.py`) run against a live Docker
      graph — **deferred to a Docker machine** (none on this box).
- [x] Unit tests pass (44). KG-validity tests written (auto-skip without Neo4j).
      `/omics-edge-snapshot` regression check → run on the Docker rebuild.
- [x] CLAUDE.md updated.

## Final status (this machine, no Docker)

All three stages implemented, tested (44 unit tests), committed; merged with 56 upstream
commits (2 conflicts resolved: CLAUDE.md, post-import.cypher). Corpus: 41/41 papers
extracted (chunking handled large PDFs), 1127/1327 genes + 145/218 pathways resolved.
**Must run on a Docker box:** full `create_knowledge_graph.py` + Docker import + post-import
+ `pytest -m kg` (`test_discuss_edges.py`) to confirm zero dangling endpoints + the
post-import rollups, plus `/omics-edge-snapshot` for no expression regression.

## Open questions for implement phase

- Whether `prominence` should also gate edge emission (e.g. drop pure `peripheral`
  name-drops) or keep all and let the consumer rank — default: keep all, rank downstream.
