# Gene Cluster Nodes — Design Spec

## Summary

Add `GeneCluster` nodes to the knowledge graph to capture published co-expression clusters, diel periodicity groups, and expression-level classifications from paper supplementary data. Enables queries like "what cluster is gene X in?", "give me all genes in cluster Y", and "do genes X and Y co-cluster anywhere?"

## Motivation

From [gaps_and_friction.md](../../../analyses/) (nitrogen stress MED4 analysis):
> No clustering/co-expression — Gene co-expression clusters and dynamic pattern classifications from the time courses. Would enable "which cluster does this gene belong to?" queries.

A survey of all 25 papers in `data/Prochlorococcus/papers_and_supp/` found 9 papers with cluster data, 7 with machine-readable supplementary tables suitable for integration. See the data inventory below.

## Graph Model

### Node: `GeneCluster`

**ID format:** `cluster:{doi_short}:{cluster_id}` (e.g., `cluster:zinser2009:cluster_5`, `cluster:tolonen2006:med4_k3`)

| Property | Type | Required | Description | Example |
|---|---|---|---|---|
| `name` | string | yes | Cluster label from paper | "Cluster 5" |
| `source_paper` | string | yes | Short paper reference | "Zinser 2009" |
| `organism_name` | string | yes | Same vocabulary as Gene/Experiment | "Prochlorococcus MED4" |
| `cluster_method` | string | yes | Algorithm used | "Mfuzz soft clustering" |
| `cluster_type` | string | yes | Category enum | "diel_periodicity" / "stress_response" / "expression_level" |
| `treatment_type` | string[] | yes | Same vocabulary as Experiment (shared enum) | ["carbon_stress", "oxygen_stress"] |
| `treatment` | string | yes | Free-text condition description | "CO2/O2 gas shock vs air control" |
| `omics_type` | string | yes | Platform | "MICROARRAY" / "RNASEQ" / "PROTEOMICS" |
| `light_condition` | string | yes | Light regime (important for cyanobacteria) | "14:10 L:D" / "continuous light" |
| `member_count` | int | yes | Number of genes in cluster | 142 |
| `functional_description` | string | yes | What the genes ARE (enrichment summary) | "Enriched for PSI and PSII genes (FDR 1.5e-9), ATP synthase" |
| `behavioral_description` | string | yes | What the genes DO together (temporal/response pattern) | "Expression peaks at dawn, drops through day, nadir at dusk" |
| `peak_time_hours` | float | no | Peak expression time for diel clusters (null for non-diel) | 2.0 (= 2h after lights on) |
| `period_hours` | float | no | Oscillation period (null for non-periodic) | 24.0 |
| `experimental_context` | string | yes | Experimental setup details | "2h sampling over 2 days, custom Affymetrix array" |

**Design notes:**
- `functional_description` and `behavioral_description` are deliberately free-text. Each paper describes cluster behavior in different experimental frames (diel timing vs N-starvation kinetics vs CO2 shock recovery). Structured fields would be too generic or require paper-specific schemas. Free text works well for LLM-driven queries.
- `treatment_type` is a string array (not single string) because some clusters span multiple conditions (e.g., Bagby 2015 CO2+O2). This differs from Experiment nodes which currently use a single string — a separate task should migrate Experiment.treatment_type to string[] for consistency.
- Clusters belong to Publications (paper-level analytical product), not Experiments. Some clusters are cross-experiment (e.g., Bagby's 4 Mfuzz clusters derived from genes DE in any of 3 gas shock conditions).

### Edges

| Edge | Direction | Properties |
|---|---|---|
| `Publication_has_gene_cluster` | Publication → GeneCluster | (none) |
| `Gene_in_gene_cluster` | GeneCluster → Gene | `membership_score` (float, nullable — Mfuzz/SoftCluster probability), `p_value` (float, nullable — periodicity/assignment p-value) |
| `Genecluster_belongs_to_organism` | GeneCluster → OrganismTaxon | (none) |

### treatment_type Vocabulary

Current values on Experiment nodes (shared enum):

| Value | Count (Experiments) |
|---|---|
| `coculture` | 16 |
| `nitrogen_stress` | 15 |
| `carbon_stress` | 11 |
| `light_stress` | 8 |
| `darkness` | 5 |
| `iron_stress` | 5 |
| `viral` | 4 |
| `phosphorus_stress` | 4 |
| `plastic_stress` | 4 |
| `salt_stress` | 3 |
| `growth_state` | 1 |

New values needed for clusters:

| Value | Paper |
|---|---|
| `diel` | Zinser 2009, Biller 2018 |
| `oxygen_stress` | Bagby 2015 (combined with `carbon_stress`) |

Note: `temperature_stress` already exists in the validator's `CANONICAL_CONDITION_TYPES` (covers Alonso-Saez 2023). Only `diel` and `oxygen_stress` are truly new.

### Post-import

- Indexes (scalar): `gene_cluster_organism_idx`, `gene_cluster_treatment_type_idx`, `gene_cluster_type_idx`
- Indexes (full-text): `geneClusterFullText` on `name`, `functional_description`, `behavioral_description`, `experimental_context`
- Computed: verify `member_count` matches actual `Gene_in_gene_cluster` edge count

## Pipeline Integration

### Paperconfig format

New `type: gene_clusters` entry under `supplementary_materials`:

```yaml
supplementary_materials:
  # existing DE tables remain unchanged...

  cluster_table_1:
    type: gene_clusters
    filename: "data/.../Table_S1_clusters.csv"
    organism: "Prochlorococcus MED4"
    gene_id_col: "ORF"
    cluster_col: "cluster"
    score_col: "membership"         # optional: fuzzy membership score column
    omics_type: MICROARRAY
    light_condition: "14:10 L:D"
    treatment_type: ["diel"]
    treatment: "Diel transcriptome, 2h sampling"
    experimental_context: "Custom Affymetrix array, 14:10 L:D cycle, 2 days"
    clusters:
      cluster_1:
        name: "Cluster 1"
        cluster_type: "diel_periodicity"
        functional_description: "PSI and PSII genes (FDR 1.5e-9)"
        behavioral_description: "Peaks at dawn, drops through day"
      cluster_5:
        name: "Cluster 5"
        cluster_type: "diel_periodicity"
        functional_description: "Respiratory terminal oxidases"
        behavioral_description: "Peaks mid-day, follows light intensity"
```

**Key design decisions:**
- Shared fields (`organism`, `omics_type`, `treatment_type`, `light_condition`, `treatment`, `experimental_context`) at table level — same for all clusters from one analysis
- Per-cluster fields (`name`, `cluster_type`, `functional_description`, `behavioral_description`) under `clusters:` — filled by reading the paper
- `gene_id_col` feeds into step 4 resolution pipeline (same as DE tables)
- For papers with separate clusters per organism (e.g., Tolonen 2006 MED4 + MIT9313), use two `type: gene_clusters` entries

### Adapter

New file: `multiomics_kg/adapters/cluster_adapter.py`

- `ClusterAdapter` — per-paper, reads one `gene_clusters` supplementary entry
- `MultiClusterAdapter` — iterates all paperconfigs, aggregates nodes/edges (same pattern as `MultiOMICSAdapter`)
- Reads resolved CSVs (from step 4) for gene ID mapping
- Emits: `GeneCluster` nodes, `Gene_in_gene_cluster` edges, `Genecluster_belongs_to_organism` edges, `Publication_has_gene_cluster` edges

### Gene ID Resolution

Cluster CSVs go through the same step 4 (`resolve_paper_ids.py`) pipeline as DE tables. The `gene_id_col` in the paperconfig specifies which column contains gene identifiers, and `resolve_paper_ids` produces `<stem>_resolved.csv` with `locus_tag` and `resolution_method` columns.

### Schema config

Add to `config/schema_config.yaml`:
- `gene cluster` node type with all properties
- `publication has gene cluster` edge type
- `gene in gene cluster` edge type with `membership_score` property
- `gene cluster belongs to organism` edge type

### Main entry point

Add to `create_knowledge_graph.py`:
- Instantiate `MultiClusterAdapter` (after omics adapter, since it reads the same paperconfigs)
- Call `write_nodes()` and `write_edges()`

## Data Inventory

### Papers with machine-readable cluster data (integration candidates)

| Paper | Organism | Clusters | Method | Type | Supp Table |
|---|---|---|---|---|---|
| **Zinser 2009** | MED4 | 16 diel + 2 aperiodic | Mfuzz | diel_periodicity | Table S1/S2 |
| **Tolonen 2006** | MED4 (9), MIT9313 (7) | 16 total | K-means | stress_response | Supp Info |
| **Alonso-Saez 2023** | MIT9301 | 5 | Fuzzy c-means (SoftCluster) | stress_response | Table S5 (XLSX) |
| **Bagby 2015** | MED4 | 4 | Mfuzz | stress_response | Table S1 / GEO |
| **Steglich 2006** | MED4 | 3 light-response groups | Manual + hierarchical | stress_response | Table 3 |
| **Ziegler 2025** | MED4 | TBD | Euclidean/pheatmap | stress_response | Supp Tables 5-8 |
| **Wang 2014** | MED4 | 5 expression classes | RPKM quartiles | expression_level | Additional file 3 |

### Papers with figure-only clusters (not integrable without manual extraction)

| Paper | Organism | Details |
|---|---|---|
| Al-Hosani 2015 | AS9601 | Hierarchical clustering heatmap, no gene membership table |
| Biller 2018 | NATL2A | Periodicity groups (Table S4A usable — revisit in Phase 3) |
| Lin 2015 | Phage P-SSM2 | Phage gene clusters only (outside KG scope) |

### Papers requiring new paperconfigs (no DE data in KG yet)

| Paper | Status | Notes |
|---|---|---|
| Zinser 2009 | No paperconfig | Definitive MED4 diel study. Microarray, MED4 PMM locus tags. Priority for Phase 1 |
| Alonso-Saez 2023 | No paperconfig | MIT9301. Expression is mRNA copies/cell, not standard log2FC/padj. Clusters integrable even without DE |
| Wang 2014 | No paperconfig | MED4. Expression subclasses (not DE). Clusters integrable without DE |

## Phasing

### Phase 1: Core model + first two papers
- Schema: add GeneCluster node and 3 edge types to `schema_config.yaml`
- Adapter: `cluster_adapter.py` with `ClusterAdapter` + `MultiClusterAdapter`
- Paperconfig format: `type: gene_clusters` support
- Gene ID resolution: extend `resolve_paper_ids.py` to handle `gene_clusters` entries
- Post-import: indexes + member_count verification
- Data: Tolonen 2006 (~16 clusters, MED4 + MIT9313 N-stress) + Zinser 2009 (~18 clusters, MED4 diel)
- Pipeline: integrate into `create_knowledge_graph.py`
- Estimated: ~34 GeneCluster nodes, ~3000-5000 Gene_in_gene_cluster edges

### Phase 2: More organisms and stress types
- Alonso-Saez 2023 (MIT9301, temperature, 5 clusters)
- Bagby 2015 (MED4, CO2/O2, 4 clusters)
- Steglich 2006 (MED4, light, 3 groups)
- Estimated: ~12 additional clusters

### Phase 3: Coculture, expression levels, periodicity
- Ziegler 2025 (MED4, coculture clusters)
- Wang 2014 (MED4, expression-level classes)
- Biller 2018 (NATL2A, diel periodicity from Table S4A)
- Estimated: ~25+ additional clusters

## Tooling Updates

### Paperconfig skill (`/paperconfig`)

Update `.claude/skills/paperconfig/SKILL.md` to document the `type: gene_clusters` format. The interactive wizard should guide users through:
- Identifying cluster tables in supplementary materials
- Specifying `gene_id_col`, `cluster_col`, `score_col`
- Filling in per-cluster descriptions (or triggering LLM extraction — see below)

### Validation script (`validate_paperconfig.py`)

Update `.claude/skills/paperconfig/validate_paperconfig.py`:
- Recognize `type: gene_clusters` as a valid supplementary entry type
- Validate required fields: `filename`, `organism`, `gene_id_col`, `cluster_col`, `clusters` block
- Validate `organism` against `CANONICAL_GENOMIC_ORGANISMS`
- Validate `omics_type` against `VALID_TYPES`
- Validate `treatment_type` array values against `CANONICAL_CONDITION_TYPES` (already has `temperature_stress`, so the new vocabulary value for Alonso-Saez 2023 is `temperature_stress` not `temperature`)
- Check CSV exists and `gene_id_col`/`cluster_col`/`score_col` columns match CSV headers
- Warn if `clusters` block is empty or has entries without `functional_description`/`behavioral_description`

Note: `CANONICAL_CONDITION_TYPES` needs new entries: `diel`, `oxygen_stress`. The existing `temperature_stress` covers Alonso-Saez 2023.

### LLM-based cluster description extraction

New script: `multiomics_kg/extraction/extract_cluster_descriptions.py`

Follows the same pattern as `multiomics_kg/adapters/pdf_publication_extraction.py` (PDF text → LLM prompt → structured JSON → cache).

**Architecture: Two-stage extraction (retrieve then synthesize)**

Based on retrospective of manually-written Tolonen 2006 descriptions, the biggest risk is cross-contamination — citing genes or stats from cluster N in the description of cluster M. Papers discuss clusters thematically (e.g., all hli genes together), not one-cluster-at-a-time. Splitting into two stages prevents this.

**Stage 1: Retrieval — find cluster references in paper text and supplementary materials**

- Input: paper PDF text + supplementary files from the paper directory + list of cluster IDs/names from the CSV
- Sources to search (in order):
  1. Main paper text (Results, Discussion sections)
  2. Figure and table legends (often contain enrichment p-values and category names)
  3. Supplementary text files (HTML, PDF, TXT in the paper directory — e.g., `clusterpvalues.html`)
  4. Supplementary table headers/metadata (sheet names, column headers in XLS that characterize clusters)
- Task: for each cluster, find all paragraphs/sentences/legends that explicitly mention it
- Output: `{cluster_key: [{quote: "...", location: "Results, paragraph 3", source_file: "main_pdf"}, ...]}`
- Can be partially non-LLM: regex search for "cluster 1", "cluster 2", etc. finds most references. LLM only needed for ambiguous references ("the upregulated clusters", "the first group").
- Validation: check coverage — which clusters have zero references? Flag for manual attention.
- The paper directory listing tells us what supplementary files are available — scan all readable formats (PDF, HTML, TXT, XLS sheet names).

**Stage 2: Synthesis — write descriptions from retrieved quotes only**

- Input: per-cluster quotes from stage 1 + cluster metadata (member count, enrichment if available)
- Task: write `functional_description` and `behavioral_description` from ONLY the provided quotes
- Output: `{cluster_key: {functional_description, behavioral_description, confidence, source}}`
- Constraint: the LLM does NOT see the full paper text — only the quotes from stage 1. This prevents cross-contamination.
- `confidence`: `"explicit"` (quotes directly describe this cluster), `"inferred"` (quotes are about a group that includes this cluster), `"not_available"` (no quotes found).

**Why two stages:**
- Stage 1 output is auditable — you can see exactly what text was found before any interpretation
- Stage 2 input is scoped — the LLM can only synthesize from retrieved quotes, not the whole paper
- Errors are traceable — wrong description? Check whether stage 1 found the wrong text or stage 2 interpreted it wrong

**Stage 1 prompt sketch:**
```
You are finding references to gene expression clusters in a scientific paper.
The paper reports {n} clusters from {method} analysis of {organism} under
{treatment}.

For each cluster listed below, find ALL sentences, paragraphs, or figure/table
legends that reference it. Include:
- Direct mentions: "cluster 1", "cluster 3 and 4", "the first cluster"
- Figure/table references: "Figure 3A shows cluster 1...", enrichment p-values
- Gene lists attributed to specific clusters

Return JSON: {cluster_key: [{quote: "exact quote", location: "section/figure"}]}
If a cluster is not mentioned anywhere, return an empty list.

Clusters: [1, 2, 3, ...]
```

**Stage 2 prompt sketch:**
```
You are writing descriptions for gene expression clusters.
For each cluster below, you are given ONLY the relevant quotes from the paper.
Do NOT add information not present in the quotes.

For each cluster, write:
1. functional_description — what types of genes are in this cluster
2. behavioral_description — the temporal/response pattern
3. confidence — "explicit" if quotes directly describe this cluster,
   "inferred" if quotes are about a broader group, "not_available" if
   no information
4. id — a short snake_case identifier summarizing the cluster (e.g.,
   "up_n_transport", "down_translation")

Cluster {key} ({n} genes):
Quotes: [...]
```

**Per-paper output file:**

Each paper gets a `cluster_extraction.json` in its directory:
`data/Prochlorococcus/papers_and_supp/<Author Year>/cluster_extraction.json`

```json
{
  "metadata": {
    "paper": "Tolonen 2006",
    "doi": "10.1038/msb4100087",
    "extracted_at": "2026-03-31T14:00:00",
    "model": "gpt-5-nano",
    "sources_searched": ["msb4100087.pdf", "Kmeans/clusterpvalues.html", "Kmeans/medSigClusters.xls"]
  },
  "clusters": {
    "1": {
      "organism": "Prochlorococcus MED4",
      "member_count": 5,
      "stage1_quotes": [
        {"quote": "Cluster 1, the most rapidly...", "location": "Results, paragraph 2", "source_file": "main_pdf"},
        {"quote": "1 transport and binding P=0.01", "location": "Figure 3 legend", "source_file": "main_pdf"}
      ],
      "stage2_result": {
        "id": "up_n_transport",
        "functional_description": "N transport genes...",
        "behavioral_description": "Rapid and strong...",
        "confidence": "explicit"
      },
      "reviewed": false,
      "reviewer_edits": null
    }
  }
}
```

This preserves stage 1 quotes for audit, tracks review status, and captures human edits (which can become few-shot examples for future extractions).

**Validation and review workflow:**

1. **Extraction run** — script produces `cluster_extraction.json` with both stages' output
2. **Auto-validation** — script flags issues automatically:
   - Clusters with zero stage 1 quotes → `confidence: "not_available"`, flagged for manual attention
   - Clusters marked `inferred` → highlighted for careful review
   - Description length sanity (too short = likely incomplete, too long = likely hallucinated)
   - Cross-check: does the described member count match the CSV?
3. **Review report** — generates a human-readable markdown summary:
   ```
   ## Cluster Extraction Review: Zinser 2009

   ### cluster_1 (142 genes) — EXPLICIT ✓  [3 quotes found]
   - functional: "PSI and PSII genes (FDR 1.5e-9)"
   - behavioral: "Peaks at dawn, drops through day"
   - sources: Figure 3 legend, Table 1 row 1

   ### cluster_7 (45 genes) — NOT AVAILABLE ⚠  [0 quotes found]
   - ACTION NEEDED: fill in manually or confirm empty

   ### cluster_12 (67 genes) — INFERRED ⚠  [1 quote, indirect]
   - functional: "Possibly ribosomal genes based on Figure 4"
   - ACTION NEEDED: verify against paper
   ```
4. **Human edits** — user edits `cluster_extraction.json` directly, sets `reviewed: true` and optionally adds `reviewer_edits` explaining changes
5. **Finalize** — script reads reviewed JSON → generates YAML `clusters:` block for paperconfig

**Iterative refinement:** When the user corrects a description, the `reviewer_edits` field captures what was changed and why. These corrections can be fed back as few-shot examples for future extractions from similar papers.

**Design decisions:**
- Semi-automated, not fully automated: the script drafts descriptions, the human reviews and edits before committing to paperconfig. This is important because papers vary in how explicitly they describe clusters, and `behavioral_description` requires accurate interpretation.
- `confidence` + `source` fields make review tractable — the human can spot-check by looking at the cited figure/table instead of re-reading the whole paper.
- Reuses the LangChain `init_chat_model` + `JsonOutputParser` pattern from `pdf_publication_extraction.py`
- Cache prevents redundant LLM calls when iterating on descriptions
- Could also be integrated into the `/paperconfig` skill as an optional step: "Would you like to auto-extract cluster descriptions from the PDF?"

## Related Work (separate tasks)

- **Migrate Experiment.treatment_type to string[]** — for consistency with GeneCluster
- **MCP tools** (in multiomics_explorer repo): `list_gene_clusters`, `gene_clusters_by_gene`, `genes_in_cluster`
- **Future: enrichment edges** — promote `functional_description` to `Cluster_enriched_for` edges linking GeneCluster → GO/KEGG/CyanoBase nodes with p-value properties
- **Future: physiological data** — most papers have growth curves / Fv/Fm but figure-only (not machine-readable)

## Not in Scope

- Computing new clusters from raw expression data (we use published clusters only)
- Enrichment as separate graph edges (captured as free-text `functional_description` for now)
- Physiological data nodes
- Diel expression values (raw time-series data from Zinser/Biller — only the cluster assignments)
- MCP tool implementation (separate repo, separate task)
