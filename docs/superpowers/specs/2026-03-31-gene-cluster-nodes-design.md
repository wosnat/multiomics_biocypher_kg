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

**Flow:**
1. Read the paper PDF (reuse `_extract_pdf_text` approach from `PDFPublicationExtractor`)
2. Read the cluster CSV to get cluster IDs and member gene counts
3. Prompt the LLM with paper text + cluster list, asking for `functional_description` and `behavioral_description` per cluster
4. Cache results as JSON (keyed by DOI + cluster table filename)
5. Output as a YAML fragment for pasting into the paperconfig `clusters:` block

**LLM prompt principles:**
- Emphasize: extract ONLY what the paper explicitly states. Do not infer, synthesize, or guess.
- For each cluster, return a `confidence` field: `"explicit"` (paper describes this cluster directly), `"inferred"` (derived from figures/context but not stated for this cluster), or `"not_available"` (paper does not describe this cluster).
- When `not_available`, return empty string for the description — do not fabricate.
- Include `source` field: quote or paraphrase the specific sentence/figure reference the description came from (e.g., "Figure 3B caption", "Results section paragraph 4", "Table 1 row header").

**Prompt sketch:**
```
You are extracting cluster descriptions from a scientific paper.
CRITICAL RULES:
- Only report what the paper EXPLICITLY states about each cluster.
- If the paper does not describe a cluster, set confidence to "not_available"
  and leave the description as "".
- Do NOT infer functional descriptions from gene names or IDs.
- Do NOT synthesize descriptions from general paper context.
- Include the source (figure, table, section) for each description.

Paper: {organism} under {treatment}. Authors identified {n} clusters
via {method}.

For each cluster, extract:
1. functional_description — what types of genes (enrichment, categories)
2. behavioral_description — temporal/response pattern (peak time, shape)
3. confidence — "explicit" | "inferred" | "not_available"
4. source — where in the paper this information comes from

Clusters found in the data:
- cluster_1: 142 genes
- cluster_2: 89 genes
- ...

Return JSON: {cluster_id: {functional_description, behavioral_description,
confidence, source}}
```

**Validation and review workflow:**

The extraction output needs human review before going into paperconfig. A multi-step workflow:

1. **Extraction run** — script produces a JSON draft with descriptions + confidence + source citations
2. **Auto-validation** — script flags issues automatically:
   - Clusters marked `not_available` → human must fill in manually or confirm empty
   - Clusters marked `inferred` → highlighted for careful review
   - Description length sanity (too short = likely incomplete, too long = likely hallucinated)
   - Cross-check: does the described member count match the CSV? (e.g., LLM says "142 genes" but CSV has 89)
3. **Review report** — generates a human-readable markdown report:
   ```
   ## Cluster Extraction Review: Zinser 2009

   ### cluster_1 (142 genes) — EXPLICIT ✓
   - functional: "PSI and PSII genes (FDR 1.5e-9)"
   - behavioral: "Peaks at dawn, drops through day"
   - source: "Table 1, row 1; Figure 2 caption"

   ### cluster_7 (45 genes) — NOT AVAILABLE ⚠
   - functional: ""
   - behavioral: ""
   - ACTION NEEDED: fill in manually or confirm empty

   ### cluster_12 (67 genes) — INFERRED ⚠
   - functional: "Possibly ribosomal genes based on Figure 4"
   - behavioral: "Appears to peak mid-day"
   - source: "Inferred from Figure 4 heatmap position"
   - ACTION NEEDED: verify against paper
   ```
4. **Human edits** — user edits the report or the draft JSON directly
5. **Finalize** — script converts reviewed JSON → YAML `clusters:` block for paperconfig

**Optional: iterative refinement** — if the user corrects a description, the correction can be fed back as a few-shot example for future extractions from the same paper or similar papers. Not essential for Phase 1 but would improve quality over time.

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
