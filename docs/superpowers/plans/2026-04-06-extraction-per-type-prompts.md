# Per-Type Extraction Prompts Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Get cluster extraction to deployment-ready quality with per-type prompts, paperconfig hints, and field renames.

**Architecture:** Refactor `extract.py` prompt into shared rules + per-type templates. Add `extraction:` top-level section to paperconfig for scope/additional PDFs. Rename extraction fields. Update paperconfig validator. Re-extract all 115 clusters.

**Tech Stack:** Python, Pydantic, OpenAI Responses API (gpt-4.1-mini), YAML paperconfigs

**Spec:** `docs/superpowers/specs/2026-04-06-extraction-per-type-prompts-design.md`

---

### Task 1: Rename Pydantic fields in extraction schema

**Files:**
- Modify: `multiomics_kg/extraction/cluster/extract.py:44-57`

- [ ] **Step 1: Rename fields in `ClusterExtraction` class**

In `multiomics_kg/extraction/cluster/extract.py`, replace:

```python
class ClusterExtraction(BaseModel):
    id: str
    name: str
    functional_description: str
    behavioral_description: str
    direction: Literal["up", "down", "mixed", "not_described"]
    enrichment_category: str
    enrichment_pvalue: Optional[float]
    enrichment_significant: bool
    confidence_notes: str
    supporting_quotes: list[SupportingQuote]
    source_figures: list[str]
    self_assessment: Literal["high", "medium", "low"]
    assessment_notes: str
```

with:

```python
class ClusterExtraction(BaseModel):
    id: str
    name: str
    functional_description: str
    temporal_pattern: str
    expression_dynamics: str
    enrichment_category: str
    enrichment_pvalue: Optional[float]
    enrichment_significant: bool
    confidence_notes: str
    supporting_quotes: list[SupportingQuote]
    source_figures: list[str]
    self_assessment: Literal["high", "medium", "low"]
    assessment_notes: str
```

Changes: `behavioral_description` → `temporal_pattern`, `direction` (Literal enum) → `expression_dynamics` (free-text str).

- [ ] **Step 2: Verify import works**

Run: `uv run python -c "from multiomics_kg.extraction.cluster.extract import ClusterExtraction; print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add multiomics_kg/extraction/cluster/extract.py
git commit -m "extract: rename direction→expression_dynamics, behavioral_description→temporal_pattern"
```

---

### Task 2: Update paperconfig validator

**Files:**
- Modify: `.claude/skills/paperconfig/validate_paperconfig.py:160-170` (VALID_CLUSTER_TYPES)
- Modify: `.claude/skills/paperconfig/validate_paperconfig.py:173-175` (REQUIRED_CLUSTER_TABLE_FIELDS area)
- Modify: `.claude/skills/paperconfig/validate_paperconfig.py:733-737` (top-level key validation)

- [ ] **Step 1: Update `VALID_CLUSTER_TYPES`**

In `.claude/skills/paperconfig/validate_paperconfig.py`, replace:

```python
VALID_CLUSTER_TYPES = {
    "diel_cycle",
    "diel_cycling",
    "time_series_dynamics",
    "response_pattern",
    "expression_pattern",
    "expression_level",
    "expression_classification",
    "periodicity_classification",
    "diel_expression_pattern",
}
```

with:

```python
VALID_CLUSTER_TYPES = {
    "time_course",
    "diel",
    "condition_comparison",
    "classification",
}
```

- [ ] **Step 2: Accept `extraction_notes` on gene_clusters entries**

Find the `REQUIRED_CLUSTER_TABLE_FIELDS` list (line ~173). Below it, add an optional fields set that the validator should not warn about:

```python
OPTIONAL_CLUSTER_TABLE_FIELDS = {
    "cluster_method", "omics_type", "light_condition", "treatment_type",
    "background_factors", "treatment", "experimental_context", "experiments",
    "id_columns", "time_points", "skip_rows", "figure_hint", "extraction_notes",
    "score_col", "p_value_col",
}
```

Then in the `_validate_gene_clusters_table()` function (around line 580), after the required fields check, add a check for unknown fields:

```python
    known_fields = set(REQUIRED_CLUSTER_TABLE_FIELDS) | OPTIONAL_CLUSTER_TABLE_FIELDS | {"type"}
    for field in table:
        if field not in known_fields:
            warnings.append(f"  [{key}] Unknown field '{field}' in gene_clusters entry")
```

- [ ] **Step 3: Accept top-level `extraction:` section**

In the `validate()` function (around line 733), after the `publication` key check, add acceptance of the `extraction:` section. Find the section where unknown top-level keys would be flagged. The validator currently only checks for `publication`. Add validation for the `extraction` key:

After the `pub = get_publication(config)` line, add:

```python
    # Validate extraction section (optional)
    extraction = config.get("extraction", {})
    if extraction:
        valid_extraction_keys = {"scope", "additional_pdfs"}
        for ek in extraction:
            if ek not in valid_extraction_keys:
                warnings.append(f"  extraction: unknown key '{ek}'")
        scope = extraction.get("scope", "analysis")
        if scope not in ("paper", "analysis"):
            errors.append(f"  extraction.scope must be 'paper' or 'analysis', got '{scope}'")
        additional_pdfs = extraction.get("additional_pdfs", [])
        if additional_pdfs:
            if not isinstance(additional_pdfs, list):
                errors.append("  extraction.additional_pdfs must be a list")
            else:
                for pdf_path in additional_pdfs:
                    if not Path(pdf_path).exists():
                        warnings.append(f"  extraction.additional_pdfs: file not found: {pdf_path}")
```

- [ ] **Step 4: Run validator smoke test**

Run: `uv run python .claude/skills/paperconfig/validate_paperconfig.py data/Prochlorococcus/papers_and_supp/tolonen\ 2006/paperconfig.yaml`
Expected: PASS (no errors on existing config, before paperconfig changes)

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/paperconfig/validate_paperconfig.py
git commit -m "validate: update cluster types (7→4), accept extraction section and extraction_notes"
```

---

### Task 3: Update paperconfigs — cluster_type consolidation

**Files:**
- Modify: `data/Prochlorococcus/papers_and_supp/tolonen 2006/paperconfig.yaml`
- Modify: `data/Prochlorococcus/papers_and_supp/lindell 2007/paperconfig.yaml`
- Modify: `data/Prochlorococcus/papers_and_supp/zinser 2009/paperconfig.yaml`
- Modify: `data/Prochlorococcus/papers_and_supp/alonso 2023/paperconfig.yaml`
- Modify: `data/Prochlorococcus/papers_and_supp/Wang 2014/paperconfig.yaml`
- Modify: `data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml`
- Modify: `data/Prochlorococcus/papers_and_supp/coe 2024/paperconfig.yaml`
- Modify: `data/Synechococcus/papers_and_supp/Bernstein 2017/paperconfig.yaml`

- [ ] **Step 1: Tolonen 2006 — `response_pattern` → `time_course`**

In `data/Prochlorococcus/papers_and_supp/tolonen 2006/paperconfig.yaml`, replace all `cluster_type: "response_pattern"` with `cluster_type: "time_course"` (2 entries: `med4_kmeans_nstarvation` and `mit9313_kmeans_nstarvation`).

- [ ] **Step 2: Lindell 2007 — `response_pattern` → `time_course`**

In `data/Prochlorococcus/papers_and_supp/lindell 2007/paperconfig.yaml`, replace `cluster_type: "response_pattern"` with `cluster_type: "time_course"` (1 entry).

- [ ] **Step 3: Zinser 2009 — `diel_cycling` → `diel`**

In `data/Prochlorococcus/papers_and_supp/zinser 2009/paperconfig.yaml`, replace `cluster_type: "diel_cycling"` with `cluster_type: "diel"` (1 entry).

- [ ] **Step 4: Coe 2024 — `diel_expression_pattern` → `diel`**

In `data/Prochlorococcus/papers_and_supp/coe 2024/paperconfig.yaml`, replace all `cluster_type: "diel_expression_pattern"` with `cluster_type: "diel"` (2 entries).

- [ ] **Step 5: Alonso-Saez 2023 — `expression_pattern` → `condition_comparison`**

In `data/Prochlorococcus/papers_and_supp/alonso 2023/paperconfig.yaml`, replace `cluster_type: "expression_pattern"` with `cluster_type: "condition_comparison"` (1 entry).

- [ ] **Step 6: Wang 2014 — `expression_level` → `condition_comparison`**

In `data/Prochlorococcus/papers_and_supp/Wang 2014/paperconfig.yaml`, replace `cluster_type: "expression_level"` with `cluster_type: "condition_comparison"` (1 entry).

- [ ] **Step 7: Biller 2018 — `periodicity_classification`/`expression_classification` → `classification`**

In `data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml`, replace `cluster_type: "periodicity_classification"` and `cluster_type: "expression_classification"` with `cluster_type: "classification"` (3 entries).

- [ ] **Step 8: Bernstein 2017 — `response_pattern` → `condition_comparison`**

In `data/Synechococcus/papers_and_supp/Bernstein 2017/paperconfig.yaml`, replace all `cluster_type: "response_pattern"` with `cluster_type: "condition_comparison"` (4 entries).

- [ ] **Step 9: Verify all cluster_type values are new vocabulary**

Run: `grep -r 'cluster_type:' data/ --include='*.yaml' | grep -v '#'`
Expected: Only `time_course`, `diel`, `condition_comparison`, `classification` values.

- [ ] **Step 10: Commit**

```bash
git add data/Prochlorococcus/papers_and_supp/*/paperconfig.yaml
git add data/Synechococcus/papers_and_supp/*/paperconfig.yaml
git commit -m "paperconfig: consolidate cluster_type vocabulary (7→4)"
```

---

### Task 4: Update paperconfigs — Bernstein cluster_col and hints

**Files:**
- Modify: `data/Synechococcus/papers_and_supp/Bernstein 2017/paperconfig.yaml`

- [ ] **Step 1: Change `cluster_col` for light entries**

In `data/Synechococcus/papers_and_supp/Bernstein 2017/paperconfig.yaml`, for `bp1_light_clusters` (line 30) and `mruber_light_clusters` (line 78), change:

```yaml
      cluster_col: "ClustID_light"
```

to:

```yaml
      cluster_col: "Clust name_light"
```

- [ ] **Step 2: Change `cluster_col` for oxygen entries**

For `bp1_oxygen_clusters` (line 54) and `mruber_oxygen_clusters` (line 104), change:

```yaml
      cluster_col: "ClustID_ox"
```

to:

```yaml
      cluster_col: "Clust name_ox"
```

- [ ] **Step 3: Add `figure_hint` to all 4 entries**

Add `figure_hint` to each entry after `experimental_context`:

- `bp1_light_clusters`: `figure_hint: "Figure 6, pages 6-7"`
- `bp1_oxygen_clusters`: `figure_hint: "Figure 6, pages 6-7"`
- `mruber_light_clusters`: `figure_hint: "Figure 6, pages 6-7"`
- `mruber_oxygen_clusters`: `figure_hint: "Figure 6, pages 6-7"`

- [ ] **Step 4: Add `extraction_notes` to all 4 entries**

Add `extraction_notes` to each entry:

For bp1_light_clusters and bp1_oxygen_clusters:
```yaml
      extraction_notes: "Paper discusses clusters A-D (irradiance-responsive) and E-H (pO2-responsive) as joint T. elongatus + M. ruber clusters. This entry covers only the T. elongatus genes. Map cluster letters to this organism's subset."
```

For mruber_light_clusters and mruber_oxygen_clusters:
```yaml
      extraction_notes: "Paper discusses clusters A-D (irradiance-responsive) and E-H (pO2-responsive) as joint T. elongatus + M. ruber clusters. This entry covers only the M. ruber genes. Map cluster letters to this organism's subset."
```

- [ ] **Step 5: Verify CSV column exists**

Run: `head -1 "data/Synechococcus/papers_and_supp/Bernstein 2017/data set s2 clustered_Avg_Cond_Light sys002172092sd2_bp1.csv" | tr ',' '\n' | grep -i clust`
Expected: Shows `ClustID_light`, `Clust name_light`, `ClustID_ox`, `Clust name_ox`

- [ ] **Step 6: Commit**

```bash
git add data/Synechococcus/papers_and_supp/Bernstein\ 2017/paperconfig.yaml
git commit -m "paperconfig: Bernstein 2017 — use cluster name column, add figure hints and extraction notes"
```

---

### Task 5: Update paperconfigs — extraction section and Coe/Tolonen hints

**Files:**
- Modify: `data/Prochlorococcus/papers_and_supp/coe 2024/paperconfig.yaml`
- Modify: `data/Prochlorococcus/papers_and_supp/tolonen 2006/paperconfig.yaml`

- [ ] **Step 1: Coe 2024 — add `extraction:` section with supplementary PDF**

At the end of `data/Prochlorococcus/papers_and_supp/coe 2024/paperconfig.yaml`, add:

```yaml

extraction:
  additional_pdfs:
    - data/Prochlorococcus/papers_and_supp/coe 2024/coe_isme_supp_final_rev_final_ycae131.pdf
```

- [ ] **Step 2: Coe 2024 — add `extraction_notes` to both cluster entries**

Add to both `supp_table_3_darktolerant_clusters` and `supp_table_3_parental_clusters`:

```yaml
      extraction_notes: "Individual clusters are not discussed in the main paper text. The supplementary PDF may contain additional cluster characterization. If clusters are not described, use N/A."
```

- [ ] **Step 3: Tolonen 2006 — add `extraction:` section with paper scope**

At the end of `data/Prochlorococcus/papers_and_supp/tolonen 2006/paperconfig.yaml`, add:

```yaml

extraction:
  scope: paper
```

- [ ] **Step 4: Commit**

```bash
git add data/Prochlorococcus/papers_and_supp/coe\ 2024/paperconfig.yaml
git add data/Prochlorococcus/papers_and_supp/tolonen\ 2006/paperconfig.yaml
git commit -m "paperconfig: add extraction sections — Coe supp PDF, Tolonen paper scope"
```

---

### Task 6: Run paperconfig validation

**Files:** None (verification only)

- [ ] **Step 1: Validate all paperconfigs**

Run: `uv run python .claude/skills/paperconfig/validate_paperconfig.py --all 2>&1 | tail -30`

If `--all` is not supported, run on each paperconfig list file:

```bash
for f in $(cat data/Prochlorococcus/papers_and_supp/paperconfig_files.txt | grep -v '^#'); do
  echo "=== $f ==="
  uv run python .claude/skills/paperconfig/validate_paperconfig.py "$f" 2>&1 | tail -5
done
```

Expected: All PASS. No errors about unknown cluster_type. May have pre-existing warnings (non-canonical organisms, etc.) — those are OK.

Also validate Synechococcus paperconfigs:
```bash
for f in $(cat data/Synechococcus/papers_and_supp/paperconfig_files.txt | grep -v '^#'); do
  echo "=== $f ==="
  uv run python .claude/skills/paperconfig/validate_paperconfig.py "$f" 2>&1 | tail -5
done
```

- [ ] **Step 2: Fix any validation errors**

If any errors found, fix the paperconfig and re-validate.

---

### Task 7: Refactor prompt — shared rules + per-type templates

**Files:**
- Modify: `multiomics_kg/extraction/cluster/extract.py:67-253`

This is the largest task. Replace `DEVELOPER_MSG_TEMPLATE` and `CLUSTER_TYPE_GUIDANCE` with the new prompt architecture.

- [ ] **Step 1: Replace prompt constants**

In `multiomics_kg/extraction/cluster/extract.py`, replace everything from `DEVELOPER_MSG_TEMPLATE = """\` through the end of `CLUSTER_TYPE_GUIDANCE` dict (lines 67-215) with:

```python
SHARED_RULES = """\
## Rules

- functional_description: What genes/pathways are in this cluster?
  Can be: enrichment category with p-value, specific highlighted genes, both,
  or "N/A" when not mentioned. 2-3 sentences max.
  Only cite gene names from the paper (e.g., psbA, rbcLS), never locus tags
  (PMM*, PMT*, P9301_*, NATL2_*, PMN2A_RS*, tll*, SY28_*, BSR22_*,
  ALT831_RS*, MIT1002_*, SYNW*, sync_*, A9601_*, WP_*, cds-*).

- temporal_pattern: Describe ONLY the observed expression pattern — what happens
  to transcript levels over time or across conditions.
  Do NOT add biological interpretation (e.g., "indicating stress response",
  "suggesting metabolic coupling").
  1-2 sentences. Include timing numbers when available from the paper.

- expression_dynamics: A short free-text label summarizing the expression pattern.
  Type-specific — see guidance below.

- If the paper does NOT explicitly describe a cluster's function or behavior, set
  the field to "N/A". Do NOT say what the paper didn't do. Each field is
  independent: a cluster can have a functional_description but temporal_pattern = "N/A".

- ONLY state what the paper EXPLICITLY says about THIS SPECIFIC cluster.
  Do NOT synthesize information from other clusters or your own knowledge.
- Do NOT include treatment conditions in descriptions — those live on the analysis node.
- Do NOT describe cluster membership statistics (gene counts, sample gene IDs)
  or restate the cluster definition.
- For enrichment: only report enrichment the PAPER performed with p-values.
  Max 3 decimal places or scientific notation.
- Only cite genes the paper explicitly highlights for this cluster. Max 3-5 genes.
- self_assessment: your confidence. assessment_notes: what you're uncertain about.
- id format: {organism_short}_{dynamics}_{theme} in snake_case.
- name format: "{Organism} cluster {KEY} ({theme})" — under 60 chars.
  Use the EXACT cluster key from the list below.
- supporting_quotes: direct quotes from the paper.
- source_figures: figure/table references used as evidence.
"""

SELF_VERIFICATION = """\
## Self-Verification

Before outputting, verify for each cluster:
1. Does each supporting_quote explicitly mention THIS cluster (by number/name)?
2. If a quote discusses a broader set of genes, is it correctly attributed?
3. Would the functional_description still be accurate if you removed the quotes?
If a quote doesn't match this specific cluster, remove it and set description to "N/A".
"""

TYPE_RULES = {
    "time_course": """\
## Type: time_course — treatment/stress time-series

- expression_dynamics: response timing label — e.g., "early transient", \
"late sustained", "rapid then declining", "gradual increase"
- temporal_pattern: describe when change begins, how fast, whether it \
persists or reverses. Include timepoints from the paper.

### Example — time-course cluster with enrichment:
{{"id": "mit9313_up_transport_binding", "name": "MIT9313 cluster 1 (up, transport and binding)", \
"functional_description": "Enriched for transport and binding (p=0.04). Contains nitrogen \
transport genes urtA and the nitrite permease, and hli genes hliS and hli7.", \
"temporal_pattern": "Most rapidly and strongly upregulated cluster, with genes responding \
within the first hours of the time course.", \
"expression_dynamics": "early rapid upregulation", \
"enrichment_category": "transport and binding", \
"enrichment_pvalue": 0.04, "enrichment_significant": true, \
"self_assessment": "high", "assessment_notes": "", \
"confidence_notes": "", \
"supporting_quotes": [{{"quote": "Cluster 1, the most rapidly and highly upregulated genes \
in each strain, contains N transport genes such as MED4 and MIT9313 urtA", "location": "Page 3"}}], \
"source_figures": ["Figure 2"]}}

### Example — not discussed:
{{"id": "med4_cluster_9", "name": "MED4 cluster 9 (not discussed)", \
"functional_description": "N/A", "temporal_pattern": "N/A", \
"expression_dynamics": "N/A", "enrichment_category": "", \
"enrichment_pvalue": null, "enrichment_significant": false, \
"self_assessment": "low", "assessment_notes": "Paper does not discuss this cluster.", \
"confidence_notes": "", "supporting_quotes": [], "source_figures": []}}
""",

    "diel": """\
## Type: diel — 24h cycling

- expression_dynamics: peak phase label — e.g., "peaks at dawn", "peaks midday", \
"peaks at dusk", "peaks at night"
- temporal_pattern: describe peak timing (hours), periodicity, phase relative \
to light/dark cycle.

### Example — diel cluster with timing:
{{"id": "med4_up_photosynthesis", "name": "Prochlorococcus cluster 1 (photosynthesis)", \
"functional_description": "Enriched for photosystem I and II components (p=1.5e-9). \
Includes genes psbA, psbD, and psaA.", \
"temporal_pattern": "Genes peak in expression near dawn (8.3 h) with 24h periodicity.", \
"expression_dynamics": "peaks at dawn", \
"enrichment_category": "Photosystem I and II", \
"enrichment_pvalue": 1.5e-09, "enrichment_significant": true, \
"self_assessment": "high", "assessment_notes": "", \
"confidence_notes": "", \
"supporting_quotes": [{{"quote": "Expression of approximately half of photosystem (PS) II genes, \
including reaction center genes psbA and psbD, peak in abundance at mid-day", "location": "Page 6"}}], \
"source_figures": ["Figure 4A", "Table S4"]}}
""",

    "condition_comparison": """\
## Type: condition_comparison — across discrete conditions

- expression_dynamics: condition-response label — e.g., "up with light", \
"down at cold", "variable", "high across all"
- temporal_pattern: describe which conditions drive changes and direction/magnitude.

### Example — condition comparison cluster:
{{"id": "mit9301_up_cold_stress", "name": "MIT9301 cluster C (cold stress response)", \
"functional_description": "Enriched for stress response genes including chaperones \
(groES/groEL, dnaK, clpBCP), fatty acid desaturases (desA, desC), and oxidative \
damage protection (recA, sod).", \
"temporal_pattern": "Strongly upregulated at minimum temperature (17°C) during both \
daytime and nighttime.", \
"expression_dynamics": "up at cold", \
"enrichment_category": "stress response", \
"enrichment_pvalue": null, "enrichment_significant": false, \
"self_assessment": "high", "assessment_notes": "", \
"confidence_notes": "", \
"supporting_quotes": [{{"quote": "Clusters C and D genes included different elements \
of the global stress response, such as cellular chaperones (groES/groES, dnaK, and \
clpBCP) and fatty acid desaturases (desA and desC)", "location": "Page 4"}}], \
"source_figures": ["Figure 2", "Figure 3"]}}
""",

    "classification": """\
## Type: classification — categorizing genes by condition/periodicity

- expression_dynamics: category label — e.g., "periodic in L:D only", \
"not periodic", "present in darkness", "condition-specific"
- temporal_pattern: describe which conditions/categories the genes fall into.
- For periodicity clusters: which conditions show 24h periodicity?
- For expression classification: presence/absence pattern across conditions.

### Example — classification cluster:
{{"id": "natl2a_periodic_coculture_LD", "name": "NATL2A cluster coculture_LD (periodic)", \
"functional_description": "N/A", \
"temporal_pattern": "Genes show 24-h periodic oscillations only in coculture L:D condition, \
not in axenic or darkness.", \
"expression_dynamics": "periodic in coculture L:D only", \
"enrichment_category": "", \
"enrichment_pvalue": null, "enrichment_significant": false, \
"self_assessment": "medium", "assessment_notes": "Paper describes periodicity patterns at \
aggregate level, not per-cluster functional detail.", \
"confidence_notes": "", \
"supporting_quotes": [{{"quote": "A subset of genes maintained 24-h periodicity only \
under the coculture light:dark cycle", "location": "Page 5"}}], \
"source_figures": ["Figure 3"]}}
""",
}
```

- [ ] **Step 2: Replace `build_context_block()` and add `build_prompt()`**

Replace `build_context_block`, `get_cluster_type_guidance`, and `format_cluster_summaries` with:

```python
def build_context_block(table: dict) -> str:
    """Build context block from paperconfig entry."""
    parts = [
        f"Analysis: {table.get('name', '')}",
        f"Organism: {table.get('organism', '')}",
        f"Clustering: {table.get('cluster_method', '')}",
        f"Type: {table.get('cluster_type', '')}",
        f"Treatment: {table.get('treatment', '')}",
    ]
    if table.get("treatment_type"):
        tt = table["treatment_type"]
        if isinstance(tt, list):
            parts.append(f"Treatment categories: {', '.join(tt)}")
        else:
            parts.append(f"Treatment categories: {tt}")
    if table.get("background_factors"):
        bf = table["background_factors"]
        if isinstance(bf, list):
            parts.append(f"Background factors: {', '.join(bf)}")
        else:
            parts.append(f"Background factors: {bf}")
    if table.get("experimental_context"):
        parts.append(f"Context: {table['experimental_context']}")
    if table.get("omics_type"):
        parts.append(f"Omics: {table['omics_type']}")
    if table.get("figure_hint"):
        parts.append(f"Key figures: {table['figure_hint']}")
    if table.get("time_points"):
        parts.append(f"Time points (hours): {table['time_points']}")
    if table.get("extraction_notes"):
        parts.append(f"\nExtraction guidance: {table['extraction_notes']}")
    return "\n".join(parts)


def format_cluster_summaries(clusters: dict[str, dict]) -> str:
    """Format cluster summaries for the prompt."""
    lines = []
    for key in sorted(clusters, key=lambda x: (not x.isdigit(), x)):
        info = clusters[key]
        sample = ", ".join(info["sample_genes"][:3])
        lines.append(f"Cluster {key}: {info['gene_count']} genes (sample: {sample})")
    return "\n".join(lines)


def build_prompt(table_config: dict, cluster_summaries: dict) -> str:
    """Assemble full developer prompt from shared rules + type rules."""
    ctx = build_context_block(table_config)
    summaries = format_cluster_summaries(cluster_summaries)
    n = len(cluster_summaries)

    cluster_type = table_config.get("cluster_type", "")
    type_rules = TYPE_RULES.get(cluster_type, f"## Type: {cluster_type}\n\nNo specific guidance available for this cluster type.")

    parts = [
        f"You are extracting structured descriptions of gene expression clusters from a scientific paper.\n\n{ctx}",
        f"For each cluster, extract all fields in the output schema.\n",
        SHARED_RULES,
        type_rules,
        SELF_VERIFICATION,
        f"CRITICAL: You MUST extract EXACTLY {n} clusters, one for each cluster key "
        f"listed below. Use the EXACT cluster keys as they appear — do NOT renumber, skip, "
        f"or merge clusters. Every key must appear exactly once in your output.\n\n{summaries}",
    ]
    return "\n".join(parts)
```

- [ ] **Step 3: Update `extract_analysis()` to use `build_prompt()`**

Replace the `dev_msg` construction in `extract_analysis()`:

```python
def extract_analysis(client, file_ids, table_config, cluster_summaries, model="gpt-4.1-mini", flex=False):
    """Run extraction for one analysis. Returns (parsed, usage_dict)."""
    dev_msg = build_prompt(table_config, cluster_summaries)

    file_inputs = [{"type": "input_file", "file_id": fid} for fid in file_ids]

    kwargs = dict(
        model=model, temperature=0,
        input=[
            {"role": "developer", "content": dev_msg},
            {"role": "user", "content": file_inputs + [
                {"type": "input_text", "text": f"Extract descriptions for all {len(cluster_summaries)} clusters."},
            ]},
        ],
        text_format=AnalysisExtraction,
    )
    if flex:
        kwargs["service_tier"] = "flex"

    t0 = time.time()
    resp = client.responses.parse(**kwargs)
    elapsed = time.time() - t0

    usage = {"input_tokens": resp.usage.input_tokens, "output_tokens": resp.usage.output_tokens, "duration_sec": elapsed}
    return resp.output[0].content[0].parsed, usage
```

Note: `file_id` parameter changed to `file_ids` (list) to support additional PDFs.

- [ ] **Step 4: Update `extract_paper()` similarly**

```python
def extract_paper(client, file_ids, tables_and_summaries, model="gpt-4.1-mini", flex=False):
    """Extract all analyses for one paper in a single call."""
    context_parts = []
    all_summaries = []
    type_parts = []
    total_clusters = 0
    for table_config, cluster_summaries in tables_and_summaries:
        context_parts.append(build_context_block(table_config))
        all_summaries.append(format_cluster_summaries(cluster_summaries))
        total_clusters += len(cluster_summaries)
        ct = table_config.get("cluster_type", "")
        type_rules = TYPE_RULES.get(ct, "")
        if type_rules and type_rules not in type_parts:
            type_parts.append(type_rules)

    ctx = "\n\n".join(context_parts)
    summaries = "\n\n".join(all_summaries)
    type_section = "\n\n".join(type_parts) if type_parts else ""

    parts = [
        f"You are extracting structured descriptions of gene expression clusters from a scientific paper.\n\n{ctx}",
        f"For each cluster, extract all fields in the output schema.\n",
        SHARED_RULES,
        type_section,
        SELF_VERIFICATION,
        f"CRITICAL: You MUST extract EXACTLY {total_clusters} clusters, one for each cluster key "
        f"listed below. Use the EXACT cluster keys as they appear — do NOT renumber, skip, "
        f"or merge clusters. Every key must appear exactly once in your output.\n\n{summaries}",
    ]
    dev_msg = "\n".join(parts)

    file_inputs = [{"type": "input_file", "file_id": fid} for fid in file_ids]

    kwargs = dict(
        model=model, temperature=0,
        input=[
            {"role": "developer", "content": dev_msg},
            {"role": "user", "content": file_inputs + [
                {"type": "input_text", "text": f"Extract descriptions for all {total_clusters} clusters."},
            ]},
        ],
        text_format=AnalysisExtraction,
    )
    if flex:
        kwargs["service_tier"] = "flex"

    t0 = time.time()
    resp = client.responses.parse(**kwargs)
    elapsed = time.time() - t0

    usage = {"input_tokens": resp.usage.input_tokens, "output_tokens": resp.usage.output_tokens, "duration_sec": elapsed}
    return resp.output[0].content[0].parsed, usage
```

- [ ] **Step 5: Verify prompt renders correctly**

Run:
```bash
uv run python -c "
from multiomics_kg.extraction.cluster.extract import build_prompt
table = {'cluster_type': 'diel', 'name': 'Test', 'organism': 'Test', 'cluster_method': 'K-means', 'treatment': 'light'}
result = build_prompt(table, {'1': {'gene_count': 50, 'sample_genes': ['g1']}})
assert 'N/A' in result
assert 'peaks at dawn' in result
assert 'Self-Verification' in result
assert 'EXACTLY 1 clusters' in result
print(f'Lines: {len(result.splitlines())}')
print('OK')
"
```

Expected: `Lines: <number>` then `OK`

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/extraction/cluster/extract.py
git commit -m "extract: refactor prompt into shared rules + per-type templates with self-verification"
```

---

### Task 8: Update `find_all_entries` and CLI for extraction config

**Files:**
- Modify: `multiomics_kg/extraction/cluster/extraction_utils.py:95-105`
- Modify: `multiomics_kg/extraction/cluster/extract.py` (main function)

- [ ] **Step 1: Update `find_all_entries()` to pass extraction config**

In `multiomics_kg/extraction/cluster/extraction_utils.py`, replace `find_all_entries`:

```python
def find_all_entries() -> list[tuple[Path, str, dict, dict, dict]]:
    """Discover all gene_clusters entries across all paperconfigs.
    Returns list of (paper_dir, entry_key, table_config, pub_config, extraction_config).
    """
    entries = []
    for pc_path, pc in load_all_paperconfigs():
        paper_dir = pc_path.parent
        pub = pc.get("publication", {})
        extraction = pc.get("extraction", {})
        for entry_key, table in iter_cluster_tables(pc):
            entries.append((paper_dir, entry_key, table, pub, extraction))
    return entries
```

- [ ] **Step 2: Update `main()` to handle extraction config, scope, and additional PDFs**

In `extract.py`, update `main()`. The key changes:
1. Unpack the new 5th element (extraction config) from entries
2. Upload additional PDFs alongside main PDF
3. Use `extraction.scope` to choose `extract_paper()` vs `extract_analysis()`

Replace the extraction loop (from `by_paper` grouping through end of loop) with:

```python
    by_paper: dict[Path, list] = {}
    for paper_dir, entry_key, table, pub, extraction in entries:
        by_paper.setdefault(paper_dir, []).append((entry_key, table, pub, extraction))

    for paper_dir, group in by_paper.items():
        paper = pub_name(group[0][2])
        extraction = group[0][3]  # extraction config is per-paper

        pdf_path_str = group[0][2].get("papermainpdf", "")
        if pdf_path_str:
            pdf_path = Path(pdf_path_str)
            if not pdf_path.is_absolute():
                pdf_path = Path.cwd() / pdf_path
        else:
            pdfs = list(paper_dir.glob("*.pdf"))
            pdf_path = pdfs[0] if pdfs else None

        if not pdf_path or not pdf_path.exists():
            logger.warning(f"No PDF for {paper}, skipping")
            continue

        logger.info(f"Paper: {paper} ({len(group)} entries)")

        # Upload main PDF + additional PDFs
        file_ids = [upload_pdf(client, pdf_path)]
        for extra_pdf in extraction.get("additional_pdfs", []):
            extra_path = Path(extra_pdf)
            if not extra_path.is_absolute():
                extra_path = Path.cwd() / extra_path
            if extra_path.exists():
                file_ids.append(upload_pdf(client, extra_path))
                logger.info(f"  Uploaded additional PDF: {extra_path.name}")
            else:
                logger.warning(f"  Additional PDF not found: {extra_pdf}")

        scope = extraction.get("scope", "analysis")

        if scope == "paper":
            # Paper-level extraction: one call for all entries
            tables_and_summaries = []
            all_entry_keys = []
            for entry_key, table, pub, _ in group:
                if entry_key in list_extraction_files(paper_dir) and not args.force:
                    logger.info(f"  {entry_key}: exists, skipping")
                    continue
                summaries = load_cluster_summaries(table)
                tables_and_summaries.append((table, summaries))
                all_entry_keys.append(entry_key)
                logger.info(f"  {entry_key}: {len(summaries)} clusters")

            if tables_and_summaries:
                try:
                    parsed, usage = extract_paper(
                        client, file_ids, tables_and_summaries,
                        model=args.model, flex=args.flex,
                    )
                    # Split parsed clusters back to entries
                    all_clusters = [c.model_dump() for c in parsed.clusters]
                    offset = 0
                    for i, (table, summaries) in enumerate(tables_and_summaries):
                        n = len(summaries)
                        entry_clusters = all_clusters[offset:offset + n]
                        offset += n
                        expected_keys = set(summaries.keys())
                        matched, unmatched = match_cluster_keys(entry_clusters, expected_keys)
                        for c in unmatched:
                            logger.warning(f"    Unmatched: {c.get('name', '?')}")
                        missing = expected_keys - set(matched.keys())
                        if missing:
                            logger.warning(f"    Missing: {sorted(missing)}")
                        metadata = {
                            "paper": paper, "doi": group[i][2].get("doi"),
                            "organism": table.get("organism", ""),
                            "entry_key": all_entry_keys[i], "model": args.model,
                            "extracted_at": datetime.now().isoformat(),
                            "input_tokens": usage["input_tokens"],
                            "output_tokens": usage["output_tokens"],
                        }
                        save_extraction(paper_dir, all_entry_keys[i], metadata, matched)
                        total_clusters += len(matched)
                        logger.info(f"    {all_entry_keys[i]}: {len(matched)}/{n} matched")
                    total_in += usage["input_tokens"]
                    total_out += usage["output_tokens"]
                except Exception as e:
                    logger.error(f"  Paper-level extraction FAILED — {e}")
        else:
            # Analysis-level extraction: one call per entry
            for entry_key, table, pub, _ in group:
                if entry_key in list_extraction_files(paper_dir) and not args.force:
                    logger.info(f"  {entry_key}: exists, skipping")
                    continue

                summaries = load_cluster_summaries(table)
                logger.info(f"  {entry_key}: {len(summaries)} clusters")

                try:
                    parsed, usage = extract_analysis(
                        client, file_ids, table, summaries,
                        model=args.model, flex=args.flex,
                    )

                    expected_keys = set(summaries.keys())
                    matched, unmatched = match_cluster_keys(
                        [c.model_dump() for c in parsed.clusters], expected_keys,
                    )

                    for c in unmatched:
                        logger.warning(f"    Unmatched: {c.get('name', '?')}")
                    missing = expected_keys - set(matched.keys())
                    if missing:
                        logger.warning(f"    Missing: {sorted(missing)}")

                    metadata = {
                        "paper": paper, "doi": pub.get("doi"),
                        "organism": table.get("organism", ""),
                        "entry_key": entry_key, "model": args.model,
                        "extracted_at": datetime.now().isoformat(),
                        "input_tokens": usage["input_tokens"],
                        "output_tokens": usage["output_tokens"],
                    }
                    save_extraction(paper_dir, entry_key, metadata, matched)
                    total_in += usage["input_tokens"]
                    total_out += usage["output_tokens"]
                    total_clusters += len(matched)
                    logger.info(f"    {len(matched)}/{len(summaries)} matched")

                except Exception as e:
                    logger.error(f"  {entry_key}: FAILED — {e}")

                time.sleep(2)

        # Cleanup uploaded files
        for fid in file_ids:
            try:
                client.files.delete(fid)
            except Exception:
                pass
```

- [ ] **Step 3: Update filter/report/verify calls to handle 5-tuple**

In `main()`, update the `args.paper` and `args.entry` filters:

```python
    if args.paper:
        entries = [(d, k, t, p, e) for d, k, t, p, e in entries
                   if args.paper.lower() in p.get("papername", "").lower()]
    if args.entry:
        entries = [(d, k, t, p, e) for d, k, t, p, e in entries if k == args.entry]
```

Update `generate_report` and `verify_quality` calls — these functions take 4-tuples, so strip the 5th element:

```python
    if args.report:
        entries_4 = [(d, k, t, p) for d, k, t, p, e in entries]
        report = generate_report(entries_4)
        if args.verify:
            warnings = verify_quality(entries_4)
```

Update the dry-run section similarly:

```python
    if args.dry_run:
        print(f"\nWould process {len(entries)} entries:\n")
        for paper_dir, entry_key, table, pub, extraction in entries:
            summaries = load_cluster_summaries(table)
            exists = entry_key in list_extraction_files(paper_dir)
            status = "EXISTS" if exists else "NEW"
            scope = extraction.get("scope", "analysis")
            extra = f" [scope={scope}]" if scope != "analysis" else ""
            print(f"  {pub_name(pub)} / {entry_key}: {len(summaries)} clusters [{status}]{extra}")
        return
```

- [ ] **Step 4: Verify dry-run works**

Run: `uv run python -m multiomics_kg.extraction.cluster.extract --dry-run`
Expected: Lists 15 entries with scope info for Tolonen.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/extraction/cluster/extraction_utils.py
git add multiomics_kg/extraction/cluster/extract.py
git commit -m "extract: support extraction config — scope, additional_pdfs, extraction_notes"
```

---

### Task 9: Update verification and report for new field names

**Files:**
- Modify: `multiomics_kg/extraction/cluster/extract.py` (verify_quality and generate_report functions)

- [ ] **Step 1: Update `verify_quality()` for new field names and sentinel**

Replace `verify_quality`:

```python
def verify_quality(entries: list[tuple[Path, str, dict, dict]]) -> list[str]:
    """Run programmatic quality checks. Returns list of warning strings."""
    warnings = []

    locus_pat = re.compile(
        r"\b("
        r"PMM\d{3,}|PMT\d{3,}|P9301_\d+|tll\d{3,}|SY28_\d+|BSR22_\d+"
        r"|NATL[12]_\d+|PMN2A_RS\d+|ALT831_RS\d+|MIT1002_\d+"
        r"|SYNW\d{4}|sync_\d{4}|A9601_\d+"
        r"|WP_\d{6,}"
        r"|cds-[A-Z]{2}_\d+"
        r")\b"
    )

    filler_phrases = [
        "likely", "possibly", "not described but",
        "specific functions are not detailed",
        "not explicitly described",
    ]

    for paper_dir, entry_key, table_config, pub in entries:
        clusters = load_extraction(paper_dir, entry_key)
        if not clusters:
            continue
        paper = pub_name(pub)
        prefix = f"[{paper} / {entry_key}"

        # Check 1: locus tags in descriptions
        for key, c in clusters.items():
            for field in ("functional_description", "temporal_pattern"):
                text = c.get(field, "")
                m = locus_pat.search(text)
                if m:
                    warnings.append(
                        f"{prefix} / cluster {key}] locus tag in {field}: {m.group()}"
                    )

        # Check 2: filler on low-confidence clusters
        for key, c in clusters.items():
            if c.get("self_assessment") == "low":
                for field in ("functional_description", "temporal_pattern"):
                    text = c.get(field, "")
                    if text and text != "N/A":
                        warnings.append(
                            f"{prefix} / cluster {key}] low confidence but {field} "
                            f"is not 'N/A': {text[:60]}..."
                        )

        # Check 3: filler phrases in descriptions
        for key, c in clusters.items():
            for field in ("functional_description", "temporal_pattern"):
                text = c.get(field, "").lower()
                for phrase in filler_phrases:
                    if phrase in text:
                        warnings.append(
                            f"{prefix} / cluster {key}] filler phrase '{phrase}' in {field}"
                        )
                        break

        # Check 4: near-identical descriptions within analysis
        descs = {}
        for key, c in clusters.items():
            fd = c.get("functional_description", "")
            if fd and fd != "N/A" and len(fd) > 50:
                prefix_50 = fd[:50]
                if prefix_50 in descs:
                    warnings.append(
                        f"{prefix} / cluster {key}] near-identical functional_description "
                        f"as cluster {descs[prefix_50]}"
                    )
                else:
                    descs[prefix_50] = key

    return warnings
```

Note: removed the "empty direction" check (Check 5) — `direction` no longer exists. `expression_dynamics` is free-text and can be "N/A".

- [ ] **Step 2: Update `generate_report()` for new field names**

Replace `generate_report`:

```python
def generate_report(entries: list[tuple[Path, str, dict, dict]]) -> str:
    """Generate diff-friendly markdown report from existing extraction JSONs."""
    lines = ["# Cluster Extraction Report\n"]

    for paper_dir, entry_key, table_config, pub in sorted(entries, key=lambda e: (pub_name(e[3]), e[1])):
        clusters = load_extraction(paper_dir, entry_key)
        if not clusters:
            continue

        paper = pub_name(pub)
        lines.append(f"## {paper} / {entry_key}\n")

        for key in sorted(clusters, key=lambda x: (not x.isdigit(), x)):
            c = clusters[key]
            dynamics = c.get("expression_dynamics", "")
            assessment = c.get("self_assessment", "")
            lines.append(f"### Cluster {key} | {dynamics} | {assessment}")
            lines.append(f"**Name:** {c.get('name', '')}")
            enrich = c.get("enrichment_category", "")
            pval = c.get("enrichment_pvalue", "N/A")
            sig = c.get("enrichment_significant", "")
            lines.append(f"**Enrichment:** {enrich} (p={pval}, sig={sig})")
            lines.append(f"**Functional:** {c.get('functional_description', '')}")
            lines.append(f"**Temporal pattern:** {c.get('temporal_pattern', '')}")
            notes = c.get("confidence_notes", "")
            if notes:
                lines.append(f"**Confidence notes:** {notes}")
            assess = c.get("assessment_notes", "")
            if assess:
                lines.append(f"**Assessment notes:** {assess}")
            figs = c.get("source_figures", [])
            if figs:
                lines.append(f"**Sources:** {', '.join(figs)}")
            quotes = c.get("supporting_quotes", [])
            if quotes:
                lines.append("**Quotes:**")
                for q in quotes:
                    lines.append(f"- [{q.get('location', '')}] {q.get('quote', '')}")
            lines.append("")

    return "\n".join(lines)
```

- [ ] **Step 3: Commit**

```bash
git add multiomics_kg/extraction/cluster/extract.py
git commit -m "extract: update verify and report for new field names and N/A sentinel"
```

---

### Task 10: Update tests

**Files:**
- Modify: `tests/test_extract.py`
- Modify: `tests/test_extraction_utils.py`

- [ ] **Step 1: Update `test_extract.py` for new field names**

In `tests/test_extract.py`, make these replacements throughout the file:

- `"behavioral_description"` → `"temporal_pattern"`
- `"direction"` → `"expression_dynamics"` (in cluster dicts only, not in assertion strings)
- `"Not discussed in paper."` → `"N/A"`
- `"empty direction"` → remove `test_detect_empty_direction` entirely (direction no longer exists)
- `build_context_block` test: change `"response_pattern"` → `"time_course"` in the table dict

Replace the full file:

```python
# tests/test_extract.py
import json
import pytest
from pathlib import Path


def test_build_context_block():
    """Context block includes organism, method, treatment, extraction_notes."""
    from multiomics_kg.extraction.cluster.extract import build_context_block

    table = {
        "name": "MIT9313 N-starvation",
        "organism": "Prochlorococcus MIT9313",
        "cluster_method": "K-means (K=7)",
        "cluster_type": "time_course",
        "treatment": "N-starvation time course",
        "omics_type": "MICROARRAY",
        "extraction_notes": "Paper discusses both MED4 and MIT9313 jointly.",
    }
    result = build_context_block(table)
    assert "Prochlorococcus MIT9313" in result
    assert "K-means (K=7)" in result
    assert "N-starvation time course" in result
    assert "MICROARRAY" in result
    assert "Paper discusses both MED4 and MIT9313 jointly." in result


def test_build_prompt():
    """build_prompt assembles shared rules + type rules."""
    from multiomics_kg.extraction.cluster.extract import build_prompt

    table = {
        "cluster_type": "diel",
        "name": "Test",
        "organism": "Test",
        "cluster_method": "K-means",
        "treatment": "light",
    }
    summaries = {"1": {"gene_count": 50, "sample_genes": ["g1"]}}
    result = build_prompt(table, summaries)
    assert "N/A" in result
    assert "peaks at dawn" in result
    assert "Self-Verification" in result
    assert "EXACTLY 1 clusters" in result


def test_format_cluster_summaries():
    """Summary text has all cluster keys."""
    from multiomics_kg.extraction.cluster.extract import format_cluster_summaries

    clusters = {
        "1": {"gene_count": 7, "sample_genes": ["PMT001", "PMT002"]},
        "6": {"gene_count": 81, "sample_genes": ["PMT100"]},
    }
    result = format_cluster_summaries(clusters)
    assert "Cluster 1:" in result
    assert "Cluster 6:" in result
    assert "7 genes" in result
    assert "81 genes" in result


def test_generate_report_stable_order(tmp_path):
    """Report generated twice is identical (no ordering jitter)."""
    from multiomics_kg.extraction.cluster.extract import generate_report
    from multiomics_kg.extraction.cluster.extraction_utils import save_extraction

    metadata = {"paper": "Test 2006"}
    clusters = {
        "2": {"id": "t_2", "name": "Cluster 2", "expression_dynamics": "late sustained",
              "self_assessment": "medium", "functional_description": "Photo",
              "temporal_pattern": "Down", "confidence_notes": ""},
        "1": {"id": "t_1", "name": "Cluster 1", "expression_dynamics": "early transient",
              "self_assessment": "high", "functional_description": "Transport",
              "temporal_pattern": "Up", "confidence_notes": ""},
    }
    save_extraction(tmp_path, "test_entry", metadata, clusters)

    entries = [(tmp_path, "test_entry", {}, {"papername": "Test 2006"})]
    report1 = generate_report(entries)
    report2 = generate_report(entries)
    assert report1 == report2
    assert report1.index("Cluster 1") < report1.index("Cluster 2")


def test_detect_filler_on_low_confidence():
    """Low-confidence cluster with non-N/A description produces warning."""
    from multiomics_kg.extraction.cluster.extract import verify_quality
    from multiomics_kg.extraction.cluster.extraction_utils import save_extraction

    import tempfile
    with tempfile.TemporaryDirectory() as td:
        paper_dir = Path(td)
        clusters = {
            "1": {"id": "test_low", "name": "C1",
                  "functional_description": "Some vague filler description here",
                  "temporal_pattern": "N/A",
                  "expression_dynamics": "N/A", "self_assessment": "low"},
        }
        save_extraction(paper_dir, "test", {"paper": "Test"}, clusters)
        entries = [(paper_dir, "test", {}, {"papername": "Test"})]
        warnings = verify_quality(entries)
        assert any("low confidence" in w for w in warnings)


def test_detect_locus_tags():
    """Locus tag in description produces warning."""
    from multiomics_kg.extraction.cluster.extract import verify_quality
    from multiomics_kg.extraction.cluster.extraction_utils import save_extraction

    import tempfile
    with tempfile.TemporaryDirectory() as td:
        paper_dir = Path(td)
        clusters = {
            "1": {"id": "x", "name": "C1",
                  "functional_description": "Includes gene PMM0042 involved in transport",
                  "temporal_pattern": "Upregulated",
                  "expression_dynamics": "up"},
        }
        save_extraction(paper_dir, "test", {"paper": "Test"}, clusters)
        entries = [(paper_dir, "test", {}, {"papername": "Test"})]
        warnings = verify_quality(entries)
        assert any("locus tag" in w and "PMM0042" in w for w in warnings)
```

- [ ] **Step 2: Update `test_extraction_utils.py` for new field names**

In `tests/test_extraction_utils.py`, update the roundtrip test fixture:

Replace:

```python
    clusters = {
        "1": {
            "id": "test_up",
            "name": "Test cluster 1",
            "functional_description": "Genes involved in transport",
            "behavioral_description": "Upregulated early",
            "source_figures": ["Figure 2"],
        },
    }
```

with:

```python
    clusters = {
        "1": {
            "id": "test_up",
            "name": "Test cluster 1",
            "functional_description": "Genes involved in transport",
            "temporal_pattern": "Upregulated early",
            "expression_dynamics": "early transient",
            "source_figures": ["Figure 2"],
        },
    }
```

- [ ] **Step 3: Run tests**

Run: `uv run pytest tests/test_extract.py tests/test_extraction_utils.py -v`
Expected: All tests pass.

- [ ] **Step 4: Commit**

```bash
git add tests/test_extract.py tests/test_extraction_utils.py
git commit -m "test: update extraction tests for field renames and N/A sentinel"
```

---

### Task 11: Run full test suite

**Files:** None (verification only)

- [ ] **Step 1: Run unit tests**

Run: `uv run pytest -m "not slow and not kg" -v --tb=short 2>&1 | tail -40`
Expected: All tests pass. Watch for failures in `test_cluster_adapter.py` — the adapter still reads `behavioral_description` from extraction JSON. Since we renamed the field in extraction output but the adapter hasn't changed (deferred), new extractions will write `temporal_pattern` but the adapter expects `behavioral_description`. This is OK for now — the adapter `.get("behavioral_description", "")` will just return empty string for new extraction data. Existing extraction JSONs still have `behavioral_description`.

- [ ] **Step 2: Run paperconfig validation tests**

Run: `uv run pytest tests/test_paperconfig_validation.py -v --tb=short`
Expected: All tests pass with the new cluster types.

- [ ] **Step 3: Fix any failures**

If any test fails, fix and re-run.

- [ ] **Step 4: Commit fixes if needed**

Only if step 3 required changes.

---

### Task 12: Re-extract all 15 entries

**Files:** None (runs extraction CLI)

**Prerequisites:** `OPENAI_API_KEY` set in `.env`.

- [ ] **Step 1: Dry run**

Run: `uv run python -m multiomics_kg.extraction.cluster.extract --dry-run`
Expected: 15 entries listed. Tolonen shows `[scope=paper]`.

- [ ] **Step 2: Re-extract all with force**

Run: `uv run python -m multiomics_kg.extraction.cluster.extract --force 2>&1 | tee extraction_log.txt`

Watch for:
- `FAILED` lines (API errors)
- `Unmatched` / `Missing` lines (key matching issues)
- `Uploaded additional PDF` for Coe 2024
- Tolonen should show a single paper-level call
- Final summary: "Done: N clusters, X in + Y out tokens"

Expected: ~115 clusters extracted across 15 entries.

- [ ] **Step 3: Generate report with verification**

Run: `uv run python -m multiomics_kg.extraction.cluster.extract --report --verify`
Expected: `data/cluster_extraction_report.md` written.

- [ ] **Step 4: Review warnings**

Run: `tail -30 data/cluster_extraction_report.md`

Check:
1. Warning count < 10 (was 27)
2. Bernstein clusters have real descriptions (not all "N/A")
3. Coe 2024 clusters say "N/A" (not filler)
4. Biller 2018 uses classification language
5. No locus tags anywhere
6. No filler phrases in described clusters

- [ ] **Step 5: Commit baseline**

```bash
git add data/*/papers_and_supp/*/cluster_extractions/*.json
git add data/*/papers_and_supp/*/cluster_extractions/*.md
git add data/cluster_extraction_report.md
git commit -m "extract: re-extract all 115 clusters with per-type prompts (baseline v2)"
```

- [ ] **Step 6: Clean up**

```bash
rm -f extraction_log.txt
```

---

### Task 13: Iterate on prompts if needed

**Files:**
- Modify: `multiomics_kg/extraction/cluster/extract.py` (prompt tweaks)

This task is open-ended — repeat as needed based on Task 12 Step 4 review.

- [ ] **Step 1: Identify issues from report**

Look at warnings and spot-check descriptions. Common issues:
- Clusters that should say "N/A" but have filler
- Locus tags that slipped through
- Wrong expression_dynamics labels
- temporal_pattern doesn't match cluster type guidance
- Bernstein still mostly "N/A" despite hints

- [ ] **Step 2: Tweak prompt and re-extract specific entries**

```bash
# Re-extract just one paper
uv run python -m multiomics_kg.extraction.cluster.extract --paper "Bernstein 2017" --force

# Regenerate report
uv run python -m multiomics_kg.extraction.cluster.extract --report --verify

# Diff against last commit
git diff data/cluster_extraction_report.md
```

- [ ] **Step 3: Commit improvements**

```bash
git add multiomics_kg/extraction/cluster/extract.py
git add data/*/papers_and_supp/*/cluster_extractions/*.json
git add data/*/papers_and_supp/*/cluster_extractions/*.md
git add data/cluster_extraction_report.md
git commit -m "extract: prompt iteration — <describe what changed>"
```
