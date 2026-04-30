# Triage: explorer-side KG asks (2026-04-30)

**Source:** `multiomics_explorer/docs/superpowers/specs/2026-04-30-kg-side-requirements-and-questions.md`
**Triaged against:** live KG (`bolt://localhost:7687`) + this repo
**Triager:** KG side
**Status:** triage + KG-1 and KG-2 implemented (commit pending). Other items still need
explorer-side back-and-forth before implementation.

For each item: **decision** (addressable / partially addressable / not addressable),
**root cause** (where in the code/data), **proposed fix shape**, and any **caveats** the
explorer side needs to know.

## Summary table

| # | Item | Decision | Effort | Notes |
|---|---|---|---|---|
| KG-1 | `Gene.function_description = "Alternative locus ID"` | **Implemented** ✓ | S | Stub filter at `load_eggnog` + `_query_descriptions`; needs step-2 + step-5 rerun + KG rebuild |
| KG-2 | `GeneCluster "N/A"` → NULL | **Implemented** ✓ | XS | Adapter omits stub fields from props dict; needs KG rebuild |
| KG-3 | TX50_RS pipeline rerun | **Mostly NOT addressable** | — | Root cause is upstream (Cyanorak doesn't include these genes; original studies didn't measure them) |
| KG-4 | Term-level informativeness flag | **Addressable** — design choice required | M | Recommend curated YAML + post-import Cypher; needs vocabulary curation |
| KG-5 | Per-gene pipeline-scope tracking | **Addressable** — pure post-import derivation | S–M | Derive from existing properties; no new ingest work |
| KG-6 | Cluster curation coverage | **Partially addressable** | M | (b) `clustering_intent` flag is cheap; (a) re-curate uncurated analyses is expensive |
| KG-7a | AA sequence string | **Implemented** ✓ | S | Sequence on `Gene` (not `Protein`) — UniProt coverage is too low; loaded from per-strain `protein.faa` via `protein_id` join |
| KG-7b | InterPro | **Addressable but deferred** | L | Requires InterProScan run; defer pending demand evidence |
| KG-7c | AlphaFold | **NOT addressable in current architecture** | — | Defer pending demand evidence; structure files don't belong in Neo4j |
| MCP-7 | Surface Polypeptide fields | Out of KG scope | — | Explorer side |

---

## KG-1 — `function_description = "Alternative locus ID"` (HIGH)

**Decision:** **Addressable.** Clean fix at the gene-annotation merge step.

**Root cause (verified):** eggNOG `emapper.annotations` files contain the literal string
`"Alternative locus ID"` in the `Description` column for ~786 Prochlorococcus protein hits
per strain (and similar counts elsewhere). Example from `cache/data/Prochlorococcus/genomes/AS9601/eggnog/AS9601.emapper.annotations`:

```
WP_002805301.1	146891.A9601_13051	5.12e-51	161.0	2B8WT@1|root,…	1117|Cyanobacteria	S	Alternative locus ID	-	-	-	…
```

The string flows into `Gene.function_description` via `config/gene_annotations_config.yaml:234-243`:

```yaml
function_description:
  type: single
  candidates:
    - source: uniprot
      field: function_description
    - source: eggnog
      field: Description
```

When UniProt has no function description, eggNOG's stub wins. Live verification:
**100 % of the 7,240 affected genes carry `[eggnog] Alternative locus ID` in
`alternate_functional_descriptions`** (no other source contributes the string).

**Proposed fix:** Add a new transform `_tx_filter_eggnog_stubs` to
`multiomics_kg/download/utils/annotation_transforms.py` that returns `""` for known
stub strings (`"Alternative locus ID"`, etc.). Apply it to the eggNOG candidate in
`gene_annotations_config.yaml`. Rebuild step 2 (`build_gene_annotations.py`) on all
strains, then rebuild the KG.

**Adjacent stubs found in the same scan (top of `function_description` value frequency):**

| value | count | likely stub? |
|---|---:|---|
| `Alternative locus ID` | 7,240 | yes (eggNOG bug) |
| `protein conserved in bacteria` | 513 | borderline (vague but technically content) |
| `transcriptional regulator` / `Transcriptional regulator` | 495 | content but generic |
| `membrane` | 145 | borderline |
| `transcriptional` | 80 | yes (truncated) |
| `integral membrane protein` | 51 | borderline |

Recommend filtering only `"Alternative locus ID"` and `"transcriptional"` (clearly
truncated/broken) in this pass; the borderline ones are weak content but not stubs.

**Open question for the explorer side:** Want any of the borderline ones filtered too?
Or surfaced as content?

**Effort:** S. ~10 lines of code change + step 2 rebuild + KG rebuild.

---

## KG-2 — `GeneCluster "N/A"` → NULL (HIGH)

**Decision:** **Addressable.** Trivial adapter-side fix.

**Root cause (verified):** The LLM extraction prompt at
`multiomics_kg/extraction/cluster/extract.py:70-89,132-134,193,229` *instructs* the
model to write the literal string `"N/A"` when no info is available. The cluster
adapter at `multiomics_kg/adapters/cluster_adapter.py:221-229` then writes the string
through unchanged via `_clean_str(...)`.

```python
"functional_description": _clean_str(ext_data.get("functional_description", "")),
"temporal_pattern": _clean_str(ext_data.get("temporal_pattern", "")),
"expression_dynamics": _clean_str(ext_data.get("expression_dynamics", "")),
```

**Other GeneCluster fields are NOT affected** (live verification): `id`, `name` carry
no `"N/A"` values; `functional_description` has 0 NULLs and 1 empty-string. The
convention is uniform: every uncurated description is the literal `"N/A"`.

**Proposed fix:** Add a `_normalize_extraction_value(s) -> str | None` helper in
`cluster_adapter.py` that returns `None` for `s in {"N/A", "", "NA", "n/a"}` after
cleanup. Apply to the three description fields. Cluster nodes will then carry
`None` (BioCypher writes nothing), and explorer-side `IS NULL` filters work.

**Decision call-out for the explorer:** if the explorer normalizes the strings to
NULL on its side instead, the KG side fix becomes redundant. We should fix it on
the KG side anyway because all consumers benefit.

**Effort:** XS. ~5 lines in `cluster_adapter.py`. No rebuild beyond KG export.

---

## KG-3 — TX50_RS pipeline rerun (MED)

**Decision:** **Mostly NOT addressable** on the KG side. Documenting the breakdown so
the explorer team can recalibrate expectations.

**Verified state of the 14 TX50_RS genes (live):**

| og_count | gene count |
|---:|---:|
| 0 | 5 |
| 1 | 1 |
| 2 | 1 |
| 3 | 6 |
| 4 | 1 |

**eggNOG already covers 9 of 14** (og_count ≥ 1 means at least one eggNOG-derived OG).
The remaining 5 have no OG of any kind.

**Why each pipeline misses these genes:**

| pipeline | covers TX50_RS genes? | remediation possible? |
|---|---|---|
| Cyanorak ortholog grouping | No | **No.** Cyanorak is an externally curated database; we can't add genes to it. The Cyanorak GBK file (`Pro_MED4.gbk`) doesn't list these RefSeq-only loci. |
| eggNOG | Yes for 9/14 | **Partial.** Re-running eggNOG-mapper on these 5 genes would either produce hits (if the database has been updated since the last run) or confirm true "no homolog" status. Cheap but unlikely to change much. |
| Per-paper clustering analyses | No | **No.** Each clustering analysis was run on the *original study's* gene set. The studies that produced our clusters (Tolonen 2006 microarray, Lin 2023, Anderson 2024 etc.) measured PMM* locus tags only. These TX50_RS-only genes aren't in those studies' candidate sets. Re-running their clustering would require their raw data on the full RefSeq locus set, which doesn't exist. |

**The mismatch:** the spec frames KG-3 as "the pipelines never ran on these locus tags."
That's not quite right. The pipelines *did* run on the full proteome of MED4 (UniProt
proteins → eggNOG; UniProt annotations → KEGG/Pfam/etc.). The 14 TX50_RS-only genes are
genuinely under-annotated because:
1. They're hypothetical proteins with no informative homologs (5 of 14)
2. They weren't part of the original Cyanorak curation (all 14)
3. They weren't measured in any of the studies whose clustering we ingest (all 14)

**Recommendation to the explorer side:** treat these as **floor-case genes** (KG has all
the data the upstream universe has), not as "pipeline coverage gaps." KG-5 (pipeline-scope
tracking) is the right structural answer — it lets the explorer surface
`out_of_pipeline_scope` envelopes that explain "Cyanorak doesn't cover this gene" rather
than implying re-running fixes anything.

**Caveat:** running eggNOG-mapper once more on the 5 zero-OG genes is cheap and would
confirm they're truly under-annotated. We can do that as part of routine maintenance.

**Effort:** ~0 (eggNOG rerun is `bash scripts/prepare_data.sh --steps 0 --force --strains MED4`
followed by step 1, 2; or scoped to the 5 protein FASTA entries).

---

## KG-4 — Term-level informativeness flag (MED)

**Decision:** **Addressable.** Design choice between three patterns.

**Root cause (verified):** No `is_uninformative`, `informativeness_class`, or similar
flag on any ontology node label. Confirmed.

**Recommendation: hybrid pattern (option c from the spec).** Maintain a curated YAML
list in this repo, applied at post-import time:

```yaml
# config/uninformative_terms.yaml
terms:
  - id: "COG:S"   # "Function unknown"
  - id: "TIGR:Hypothetical proteins / Conserved"
  - id: "TIGR:Not Found"
  - id: "GO:0003674"  # "molecular_function" (root)
  - id: "GO:0005575"  # "cellular_component" (root)
  - id: "GO:0008150"  # "biological_process" (root)
  - prefix_match: "DUF"   # all Pfam DUF* domains
  - prefix_match: "Domain of unknown function"
  - prefix_match: "Protein of unknown function"
```

Post-import Cypher walks the list and sets `is_uninformative = "true"` on matching
nodes. Property is sparse (only present where true) — consumers test
`WHERE NOT (term.is_uninformative IS NOT NULL)` to filter.

**Why hybrid over KG-side flag (option a):**
- Single source of truth (the YAML)
- No data-modeling debt (we'd otherwise need to decide whether the flag is per-node
  or per-(term, source) — irrelevant if it's a list)
- Cheap to evolve when new catch-all terms are recognized

**Why not pure explorer-side list (option b):**
- KG-side filtering enables Cypher consumers to use it (e.g. ontology-landscape tools,
  enrichment analyses). Worth the small cost.

**Open questions for the explorer side:**
1. Do you want a binary `is_uninformative` flag, or a finer
   `informativeness_class ∈ {informative, catchall_unknown, catchall_hypothetical, root}`?
2. Should DUF/UPF Pfams really be marked uninformative? They *are* genuine domain
   annotations — they just say "we have a structural class but no function." Different
   from "Function unknown."

**Effort:** M. ~30 lines of Cypher + a curated YAML + post-import wiring. Coordination
required on the term list.

---

## KG-5 — Per-gene pipeline-scope tracking (MED)

**Decision:** **Addressable as pure post-import derivation.** No adapter work.

**Root cause (verified):** No property tracks per-gene pipeline coverage.

**Proposed fix:** Add a post-import Cypher that derives a `Gene.processed_by_pipelines: list[str]`
property from existing data. Each pipeline's "coverage" is computable from already-present
properties or edges:

| pipeline | gene is "covered" iff... |
|---|---|
| `cyanorak_curation` | `gene.locus_tag_cyanorak IS NOT NULL` |
| `eggnog_v5` | `(gene)-[:Gene_in_ortholog_group]->(og {source:'eggnog'})` exists |
| `cyanorak_orthogrouping` | `(gene)-[:Gene_in_ortholog_group]->(og {source:'cyanorak'})` exists |
| `kegg_ko_assignment` | `(gene)-[:Gene_has_kegg_ko]->()` exists |
| `pfam_assignment` | `(gene)-[:Gene_has_pfam]->()` exists |
| `go_assignment` | any of the three Gene_*_*GO edges exist |
| `cog_assignment` | `(gene)-[:Gene_in_cog_category]->()` exists |
| `tigr_role_assignment` | `(gene)-[:Gene_has_tigr_role]->()` exists |
| `cyanorak_role_assignment` | `(gene)-[:Gene_has_cyanorak_role]->()` exists |

This is a derived property; recomputed each rebuild; reflects the current state of the
graph exactly. No adapter changes; no schema modeling debate (Gene property vs Pipeline
node) — pure post-import.

**Open questions for the explorer side:**
1. Do you also need pipeline version / run date? If so, we'd need a `Pipeline` node
   per pipeline; gene → pipeline edges; metadata on the edges. That's heavier and
   probably premature.
2. Is the proposed list complete? (eggNOG / Cyanorak / KEGG / Pfam / GO / COG / TIGR /
   Cyanorak roles — anything missing?)

**Effort:** S–M. ~50 lines of post-import Cypher.

---

## KG-6 — Cluster description curation coverage (MED)

**Decision:** **Partially addressable.** Two of three sub-asks are cheap; one is
expensive.

**(a) Curate the empty analyses (NATL2A diel × 2 = 30 clusters):**
Expensive. Requires re-running the LLM cluster extraction (`multiomics_kg/extraction/cluster/extract.py`)
on the source paper, with manual review. Worth scheduling but not in this triage.

**(b) Add `clustering_intent` flag on `ClusteringAnalysis`:**
Cheap. The `MED4 gene expression level classification` analysis is a per-gene RPKM
quartile classifier — the absence of `functional_description` per cluster is **by design**
(the clusters are VEG/HEG/MEG/LEG/NEG bins, not functional groupings).

Schema addition:
```yaml
ClusteringAnalysis:
  properties:
    clustering_intent: str   # "functional" | "expression_bin" | "condition_response" | "diel_phase" | "other"
```

Adapter sets it from a new optional `clustering_intent` field in the paperconfig
`gene_clusters` entry. Default: `"functional"`.

**(c) Bin-meaning metadata on per-cluster nodes:**
Cheap, mechanism already exists. The paperconfig's per-cluster extraction JSON can
carry text like `"VEG = top RPKM quartile (highly expressed)"` in `expression_dynamics`.
This just needs paperconfig curation, not new code.

**Recommendation:** ship (b), then update existing paperconfigs that have non-functional
analyses (MED4 expression-level classification at minimum) to set `clustering_intent` and
provide bin-meaning text. Defer (a) for NATL2A re-curation.

**Caveat:** the spec's table of curation rates was true on 2026-04-30. The "MED4 phage-upregulated
transcription groups" 100 % rate is from a published paper that explicitly named its clusters;
that's the rule, not the exception. The 0 % cases (NATL2A diel) come from analyses where the
paper presented clusters as numbered groups without functional commentary. Some of those
genuinely have no functional content to extract.

**Open questions for the explorer side:**
1. Confirm `clustering_intent` enum values. Initial proposal:
   `functional | expression_bin | condition_response | diel_phase | other`.
2. Is a per-cluster `intent_role: "VEG" | "HEG" | …` flag worth having for explicit bin
   labeling? Probably overkill; the cluster `name` already carries this.

**Effort:** XS for (b) + (c) once paperconfigs updated. M for (a) per analysis.

---

## KG-7a — AA sequence string (LOW–MED)

**Decision:** **Addressable, cheap.** Per-node string limit is not binding in Neo4j 5.x.
The data is already on disk per strain — no new download or external dependency.

**Verified state:** No `aa_sequence` / `sequence` property on Polypeptide. **NCBI
protein FASTAs are already cached per strain** at
`cache/data/<organism>/genomes/<strain>/protein.faa`. Coverage on 2026-04-30:

- 25 cultured-strain genomes (10 Prochlorococcus, 7 Synechococcus/Parasynechococcus
  /Thermosynechococcus, 3 Alteromonas, 4 heterotrophs)
- **2 reference-proteome match organisms** (`Marinobacter/genomes/HP15/protein.faa`
  for MarRef v6, `Alteromonas/genomes/Alt_MarRef/protein.faa`) — same path
  convention, no special-casing needed
- (4 extra cached strains exist in `cache/data/Alteromonas/genomes/{AD45, ATCC27126,
  BGP6, BS11}` that are not currently in the graph; harmless)

Every protein-carrying OrganismTaxon node in the graph has a `protein.faa` at the
canonical path. Sample header from MED4:

```
>WP_002805124.1 MULTISPECIES: photosystem II reaction center protein I [Prochlorococcus]
MLALKISVYTIVFFFVGIFLFGFLASDPTRTPNRKDLESPQD
```

Headers are keyed by `WP_*` RefSeq protein accessions — the same ID used in
`gene_mapping.csv` to link Gene → Protein. Direct join, no ID translation.

**Implemented:** Sequence is added to **`Gene`** nodes (not `Protein`). Reason:
UniProt coverage on this graph is too low (~46 % of UniProt proteins are
orphaned; many genes have no Protein node), so a Protein-level property would
miss large fractions of each genome. Gene is the universal anchor: every CDS
with a `protein_id` in `gene_annotations_merged.json` gets a sequence when the
matching `WP_*` accession is in its strain's `protein.faa`.

The cyanorak/NCBI adapter ([`cyanorak_ncbi_adapter.py`](../multiomics_kg/adapters/cyanorak_ncbi_adapter.py))
loads `cache/data/<organism>/genomes/<strain>/protein.faa` once per strain in
`download_data()`, then injects `gene["sequence"]` from the dict during
`_get_gene_nodes()`. Schema entry added under `gene:` in
[`schema_config.yaml`](../config/schema_config.yaml).

**Cost (measured 2026-04-30 across 31 cached FAAs, 94,291 proteins):**

| metric | value |
|---|---|
| total sequence chars across cache | 29.2 MB |
| FAA files on disk | 36 MB |
| median protein length | 263 AA |
| mean protein length | 309 AA |
| P90 length | 571 AA |
| max length (one outlier in WH8102) | 10,781 AA |

Only ~55,863 of those proteins become Polypeptide nodes (4 cached strains aren't in
the graph; some proteins filter out during ingest). KG-incorporated subset:

| metric | estimate |
|---|---|
| sequence chars in graph | ~17.3 MB |
| Neo4j string-property overhead (~18 bytes × 55,863 nodes) | ~1 MB |
| **estimated store growth** | **~18–20 MB** |

For context: the existing graph store is on the order of GB. Adding sequences is
**<1 % growth**. The 10,781 AA outlier is well within Neo4j 5.x's 4 GB per-property
string limit.

**Why on Gene, not Protein:**
- UniProt coverage gap: ~46 % of UniProt proteins are orphaned (no Protein node
  edges). Putting `sequence` on Protein would silently miss those Genes' sequences.
- Universal anchor: every CDS with a `protein_id` is a Gene; the FAA is keyed by
  the same `protein_id` (`WP_*`); no traversal needed.
- Reference-proteome organisms (Marinobacter MarRef v6, Alteromonas MarRef v6)
  have their own `protein.faa` at the canonical path, so they're covered uniformly
  with cultured strains — no special-casing.

**Per-node size:** The 4 GB string limit per property in Neo4j 5.x is not
binding; the largest sequence in our cache is 10,781 AA.

**Tests** (in [`tests/test_cyanorak_ncbi_adapter.py`](../tests/test_cyanorak_ncbi_adapter.py)):
- `TestLoadProteinSequences` (4 cases): FASTA parser handles multi-line records,
  whitespace-tokenized headers, missing files, empty files.
- `TestGeneSequence` (3 cases): gene nodes carry the right sequence; missing
  WP_ → property omitted (not empty string); no `protein.faa` → no error, no
  property.

**To activate:** rebuild the KG. No prepare_data rerun needed — the adapter
reads `protein.faa` directly at build time, and the file is already cached for
every strain in the graph.

**Effort:** S. ~50 lines net in [cyanorak_ncbi_adapter.py](../multiomics_kg/adapters/cyanorak_ncbi_adapter.py).

---

## KG-7b — InterPro (LOW)

**Decision:** **Addressable but deferred.** No demand evidence yet.

**Cost:** Running InterProScan on ~50 K proteins is ~100 hours of CPU. Output is ~5 GB.
Pfam (already integrated) is the largest single InterPro contributor.

**Recommendation:** defer. Revisit when a concrete analysis actually hits a wall that
InterPro would solve. The marginal coverage over Pfam is small for *Prochlorococcus*
proteins specifically because Pfam dominates the cyano protein space.

**Effort:** L. Run InterProScan + adapter + schema. Probably 1-2 weeks of work
including download/ingest.

---

## KG-7c — AlphaFold structure summaries (LOW)

**Decision:** **NOT addressable in the current architecture.** Defer pending demand.

**Why:** AlphaFold structures are PDB / mmCIF files, ~50 KB-1 MB each. Storing structure
*files* in Neo4j is the wrong tool. Storing *structure summary metrics* (pTM, pLDDT,
fold class) is possible but speculative — there's no current analysis that would use
these.

**If demand emerges:** the right design is probably an external blob store (S3, local FS)
with `Polypeptide.alphafold_url` and `Polypeptide.fold_class_summary` properties pointing
into it.

**Open question for the explorer side:** Can you point to a single analysis that would
use AlphaFold output? If not, defer indefinitely.

**Effort:** L–XL once a design exists.

---

## MCP-7 — Surface Polypeptide fields via explorer

Out of KG-side scope. Documented for cross-reference. The KG already has the data
(`Polypeptide.signal_peptide`, `transmembrane_regions`, `is_reviewed`, etc.) — the gap
is in `gene_details` not traversing to Polypeptide. Explorer-side fix.

---

## Cross-cutting observations

### Stub-string sweep (cross-cutting note in spec)

The spec's suggested cross-cutting query for stub strings is good. The same eggNOG
`Description` source that contaminates `function_description` likely contaminates
`alternate_functional_descriptions` and `gene_summary`. Sample evidence:

```
TX50_RS09260 → gene_summary = "DUF4278 domain-containing protein :: Alternative locus ID"
TX50_RS09730 → gene_summary = "hli :: possible high light inducible protein :: Alternative locus ID"
```

Fixing KG-1 (filter "Alternative locus ID" at the eggNOG merge) automatically removes
the stub from `gene_summary` (computed from the same `best_desc`) and from
`alternate_functional_descriptions` (the `_add_desc("eggnog", eg_desc)` call).

So **KG-1 + sweep is a single fix**, not two.

### Additional stubs caught by the sweep that aren't in the spec

`gene_summary` ending in `:: -` for some TX50_RS rows (`TX50_RS03815: "hypothetical protein :: -"`).
The dash is leaking from a cyanorak/NCBI product field. Worth a separate one-line strip
in `gene_summary` construction. Not in scope of this triage but the explorer team should
expect it to be a separate small fix.

---

## Implementation log (2026-04-30)

### KG-1 implemented

- New helper `is_eggnog_description_stub()` + module constant
  `EGGNOG_DESCRIPTION_STUBS = frozenset({"Alternative locus ID"})` in
  [`multiomics_kg/download/utils/annotation_transforms.py`](../multiomics_kg/download/utils/annotation_transforms.py).
- [`load_eggnog`](../multiomics_kg/download/build_gene_annotations.py) replaces stub
  Description values with `""` at read time. All downstream consumers
  (`product`, `function_description`, `alternate_functional_descriptions`,
  `gene_summary`) inherit the fix automatically.
- [`build_og_descriptions._query_descriptions`](../multiomics_kg/download/build_og_descriptions.py)
  filters stubs at write time so the OG-description cache stays clean.
- Tests: `TestIsEggnogDescriptionStub` (6 cases) +
  `TestLoadEggnog::test_strips_alternative_locus_id_stub_from_description` in
  [`tests/test_build_gene_annotations.py`](../tests/test_build_gene_annotations.py).

**To activate the fix:**
```bash
bash scripts/prepare_data.sh --steps 2 5 --force   # rebuild gene_annotations + OG cache
docker compose build && docker compose up -d       # rebuild KG
```

### KG-2 implemented

- New helper `_normalize_extraction_value()` + module constant
  `_EXTRACTION_STUB_VALUES = frozenset({"N/A", "n/a", "NA", "n.a."})` in
  [`multiomics_kg/adapters/cluster_adapter.py`](../multiomics_kg/adapters/cluster_adapter.py).
- `GeneCluster` description fields (`functional_description`, `temporal_pattern`,
  `expression_dynamics`) are now omitted from the props dict when the extraction
  returned a stub. BioCypher writes nothing → Neo4j NULL.
- Other GeneCluster fields (`id`, `name`) are unchanged because live verification
  showed 0 stubs there.
- Tests: 9 new cases in [`tests/test_cluster_adapter.py`](../tests/test_cluster_adapter.py)
  (parametrized stub coverage + integration test confirming omission from props).

**To activate the fix:** rebuild the KG (no upstream pipeline rerun needed).

### Verification

`pytest -m "not slow and not kg"` → 1672 passed, 0 failures.
Live-graph verification deferred until after rebuild.

---

## Recommended sequencing

If the KG team takes these on:

1. **Sprint 1 (HIGH, ~1 day):** KG-1 + KG-2 fixes. One gene-annotation rebuild + one
   KG rebuild. Closes the two HIGH-severity stub-pollution gaps. Touches every gene and
   every cluster.

2. **Sprint 2 (~2-3 days):** KG-5 (post-import pipeline-scope tracking). Pure derived
   property; no adapter work; ships with next KG rebuild. Unlocks the
   `out_of_pipeline_scope` envelope on the explorer side.

3. **Sprint 3 (~3-5 days):** KG-4 (term-level informativeness, hybrid pattern) +
   KG-6(b)+(c) (clustering_intent flag + paperconfig updates for the MED4 expression-level
   analysis). Coordination with the explorer team on enum values needed.

4. **Backlog:** KG-7a (AA sequence — small standalone change), KG-3 eggNOG rerun on the
   5 zero-OG TX50_RS genes (cheap maintenance), KG-6(a) NATL2A re-curation (per-analysis
   effort).

5. **Deferred until demand evidence:** KG-7b (InterPro), KG-7c (AlphaFold).

6. **Out of KG scope:** MCP-7, KG-3 Cyanorak/clustering reruns (upstream limitations,
   not pipeline-fix-able).
