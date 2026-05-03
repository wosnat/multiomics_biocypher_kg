# 2026-05-01 KG release — Explorer frictions F1-F4

## F1 — Term informativeness + annotation rollup refinement

**New fields:**
- `Gene.annotation_state: str` ∈ `{no_evidence, catch_all_only, informative_single, informative_multi}`
- `Gene.informative_annotation_types: str[]` (granular)
- `<term node>.is_uninformative: 'true' | absent` on 7 ontology types
- `Gene.annotation_quality: int` — REDEFINED as 0-3 numeric encoding of `annotation_state`

**Semantic shift on `annotation_quality`:**
- *Before*: 0-3 mixing product-name regex with arbitrary 4-source structured-count.
- *After*: pure informative-evidence richness over 8 source buckets, no product-name dependency.
- Existing `WHERE g.annotation_quality <= 1` queries silently shift meaning (today: hypothetical-named genes; refined: low-evidence genes — spirit of the filter, now correctly).

**Worked example.** Gene with `pfam_ids=[PF00712]` + `go_terms=[GO:0003674]` (MF root) + `gene_category='Unknown'`, no other functional edges:
- *Before*: structured_count=2 → score 3.
- *After*: GO MF root flagged uninformative; only Pfam contributes → informative_count=1 → score 2, state `informative_single`.

## F2 — Data source surfacing

**New nodes:** 4 `DataSource` nodes (`ncbi`, `cyanorak`, `uniprot`, `eggnog`) with metadata (scope, provenance, info_types). Node IDs prefixed with `data_source:` (e.g. `data_source:cyanorak`).

**New field:** `Gene.contributing_sources: str[]` — per-gene data-source presence (unprefixed labels: `["ncbi", "cyanorak", ...]`).

**Worked example.** TX50_RS09500 → `contributing_sources=['ncbi']`. Distinguishes "ncbi present, others didn't run/match" from "no data" cleanly.

## F3 — `cluster_type` vocab rename

`classification` → `expression_bin`. No paperconfig migrations needed (none used the old label). Validator vocab + paperconfig SKILL.md updated. Convention: `expression_bin` analyses do not carry per-cluster `functional_description`.

## F4 — Gene properties for sparsely-annotated genes

**New fields:**
- `Gene.contig: str` — chromosome / contig name from GFF seqid (required for genomic-neighbor lookups)
- `Gene.seed_ortholog: str` — eggNOG seed-ortholog identifier (taxid.identifier)
- `Gene.seed_ortholog_evalue: float`

**Honest framing.** These primarily benefit the broader sparsely-annotated population, not the strict floor (~14 TX50_RS-style genes). Strict floor cases gain `contig` + neighbor lookup; broader population gains the seed-ortholog "this resembles X (E=Y)" pointer.

## Bucket inventory

The annotation_quality / annotation_state source-bucket list has 8 live buckets as of 2026-05-02: `go`, `kegg`, `pfam`, `ec`, `role`, `reaction`, `transporter` (via `Gene_has_tcdb_family`), `cazy` (via `Gene_has_cazy_family`). TCDB and CAZy ontologies merged into main on 2026-05-02 from a separate spec ([`2026-05-01-tcdb-cazy-ontologies-design.md`](../specs/2026-05-01-tcdb-cazy-ontologies-design.md)).

## KG-A9 + KG-A10 (chemistry indexes piggybacked on this rebuild)

Two RANGE indexes added to support the explorer team's `list_metabolites` tool build:
- `kegg_term_id_idx` on `KeggTerm.id` (5,058 nodes; backs pathway lookups)
- `metabolite_hmdb_idx` on `Metabolite.hmdb_id` (1,425 / 3,025 metabolites populated; sibling parity with `metabolite_chebi_idx` / `metabolite_kegg_id_idx` / `metabolite_mnxm_idx`)

Both ONLINE after Task 13 rebuild. See `/tmp/kg-team-prompt-2026-05-03-index-asks.md` for the original ask context.
