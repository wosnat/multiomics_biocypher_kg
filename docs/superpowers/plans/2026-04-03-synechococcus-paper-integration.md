# Synechococcus + New Paper Integration Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Integrate 9 new Synechococcus coculture papers and update Prochlorococcus paperconfigs, adding ~12 new organisms and `fold_change_type` support to the omics adapter.

**Architecture:** Paper-by-paper workflow with user review gates (Phases 1-2), then infrastructure code changes for `fold_change_type` and multi-paperconfig support, followed by genome setup (Phase 3), gene ID resolution (Phase 4), and KG build + validation (Phase 5).

**Tech Stack:** Python, BioCypher, Neo4j, YAML paperconfigs, pytest

**Spec:** `docs/superpowers/specs/2026-04-03-synechococcus-paper-integration-design.md`

---

## Phase 1: Prochlorococcus Paperconfig Updates

For each paper: explore data, create/update paperconfig + README.md, present for user approval.

### Task 1.1: Capovilla 2023 — New paperconfig (DE, MIT9303 + MIT9313)

**Files:**
- Create: `data/Prochlorococcus/papers_and_supp/Capovilla 2023/paperconfig.yaml`
- Create: `data/Prochlorococcus/papers_and_supp/Capovilla 2023/README.md`

- [ ] **Step 1: Explore data files**

Read the CSV files to understand their exact structure:
```bash
head -20 "data/Prochlorococcus/papers_and_supp/Capovilla 2023/DE genes RNA-Seq_MIT9303 pnas.2213271120.sd02.csv"
head -20 "data/Prochlorococcus/papers_and_supp/Capovilla 2023/DE genes RNA-Seq_MIT9313 pnas.2213271120.sd02.csv"
```
Check if CSV has both up- and down-regulated sections, or only up. Identify: column headers, gene ID column, logFC column, whether p-values exist, number of rows.

- [ ] **Step 2: Read the paper PDF (if available) for methods context**

Check for a paper PDF in the directory. If present, read methods section for: organism strains, treatment conditions, statistical test used, what the DE comparison is (treatment vs control).

- [ ] **Step 3: Create paperconfig.yaml**

Use `/paperconfig` skill to create the config. Key decisions:
- `table_scope`: likely `significant_only` (paper says "DE genes")
- `name_col`: use `ncbi_gene_locus_tag` or `ncbi_cds_locus_tag` (whichever maps better)
- `logfc_col`: `log2 fold change`
- `adjusted_p_value_col`: likely absent — set to null or omit
- Two experiments: one for MIT9303, one for MIT9313 (different organisms, same treatment)
- MIT9303 is NOT in the KG yet — note this in README, the experiment will produce dangling edges until Phase 3 adds the genome
- `treatment_type`, `background_factors`, `omics_type`: derive from paper

- [ ] **Step 4: Create README.md**

Write a README.md in the paper directory covering:
- Paper citation and DOI
- Data files present (CSVs, xlsx, PDFs)
- What each file contains
- KG integration decision and rationale
- Known issues (MIT9303 not yet in KG, no p-values)

- [ ] **Step 5: Present for user approval**

Show the paperconfig.yaml and README.md contents to the user. Wait for approval before proceeding.

- [ ] **Step 6: Add to paperconfig_files.txt and commit**

```bash
echo "data/Prochlorococcus/papers_and_supp/Capovilla 2023/paperconfig.yaml" >> "data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"
git add "data/Prochlorococcus/papers_and_supp/Capovilla 2023/paperconfig.yaml" \
        "data/Prochlorococcus/papers_and_supp/Capovilla 2023/README.md" \
        "data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"
git commit -m "feat: add Capovilla 2023 paperconfig (MIT9303 + MIT9313 chitosan DE)"
```

---

### Task 1.2: zinser 2009 — New paperconfig (diel clusters, MED4)

**Files:**
- Create: `data/Prochlorococcus/papers_and_supp/zinser 2009/paperconfig.yaml`
- Create: `data/Prochlorococcus/papers_and_supp/zinser 2009/README.md`

- [ ] **Step 1: Explore data files**

```bash
head -5 "data/Prochlorococcus/papers_and_supp/zinser 2009/Table_S1.csv"
cat "data/Prochlorococcus/papers_and_supp/zinser 2009/supp table legends.txt"
```
Identify: cluster column name, cluster membership score column, gene ID column (`Gene or region` likely has PMM#### IDs). Count unique clusters. Check if any columns could serve as logFC (likely not — raw expression).

- [ ] **Step 2: Create paperconfig.yaml**

This paper has clustering data only (no DE). Create a paperconfig with:
- Publication metadata
- `supplementary_materials` with a `gene_clusters` entry:
  - `type: gene_clusters`
  - `gene_id_col`: the gene ID column from Table_S1.csv
  - `cluster_col`: `Cluster`
  - `score_col`: `Cluster membership score`
  - `cluster_type`: `diel_response`
  - `cluster_method`: identify from paper (likely Fourier analysis + K-means or similar)
  - `omics_type`: MICROARRAY
  - `treatment_type`: ["diel_cycle"]
  - `background_factors`: ["continuous_light"] or ["diel_cycle"] — check paper
  - No `statistical_analyses` block (no DE data)

- [ ] **Step 3: Create README.md**

Cover: diel expression time-series (50 timepoints), Fourier periodicity analysis, cluster assignments. Note that raw expression values are NOT integrated as DE edges — only cluster assignments.

- [ ] **Step 4: Present for user approval**

- [ ] **Step 5: Add to paperconfig_files.txt and commit**

```bash
echo "data/Prochlorococcus/papers_and_supp/zinser 2009/paperconfig.yaml" >> "data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"
git add "data/Prochlorococcus/papers_and_supp/zinser 2009/paperconfig.yaml" \
        "data/Prochlorococcus/papers_and_supp/zinser 2009/README.md" \
        "data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"
git commit -m "feat: add zinser 2009 paperconfig (MED4 diel clusters)"
```

---

### Task 1.3: alonso 2023 — New paperconfig (soft clusters, MED4)

**Files:**
- Create: `data/Prochlorococcus/papers_and_supp/alonso 2023/paperconfig.yaml`
- Create: `data/Prochlorococcus/papers_and_supp/alonso 2023/README.md`

- [ ] **Step 1: Explore data files**

```bash
head -5 "data/Prochlorococcus/papers_and_supp/alonso 2023/TABLE S5 softcluster membership.csv"
cat "data/Prochlorococcus/papers_and_supp/alonso 2023/supp legends.txt"
```
Identify: gene ID column (`RefSeq Locus Tag`), cluster column (`Soft cluster`), probability score columns (one per cluster A-E). Check if RefSeq locus tags are standard PMM#### format.

- [ ] **Step 2: Create paperconfig.yaml**

Clustering-only paperconfig:
- `gene_clusters` entry with `score_col` pointing to the max probability or the assigned cluster's probability
- `gene_id_col`: `RefSeq Locus Tag`
- `cluster_col`: `Soft cluster`
- Note: BV-BRC IDs (`fig|167546.4.peg.NNN`) are also available but RefSeq locus tags should map directly
- `cluster_type`: `expression_pattern` or `soft_cluster`
- Derive `treatment_type`, `omics_type` etc. from paper

- [ ] **Step 3: Create README.md**

- [ ] **Step 4: Present for user approval**

- [ ] **Step 5: Add to paperconfig_files.txt and commit**

```bash
echo "data/Prochlorococcus/papers_and_supp/alonso 2023/paperconfig.yaml" >> "data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"
git add "data/Prochlorococcus/papers_and_supp/alonso 2023/paperconfig.yaml" \
        "data/Prochlorococcus/papers_and_supp/alonso 2023/README.md" \
        "data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"
git commit -m "feat: add alonso 2023 paperconfig (MED4 soft clusters)"
```

---

### Task 1.4: coe 2024 — Add gene_clusters to existing paperconfig

**Files:**
- Modify: `data/Prochlorococcus/papers_and_supp/coe 2024/paperconfig.yaml`
- Create: `data/Prochlorococcus/papers_and_supp/coe 2024/README.md`

- [ ] **Step 1: Explore the new CSV**

```bash
head -10 "data/Prochlorococcus/papers_and_supp/coe 2024/Supplemental Table 3.csv"
wc -l "data/Prochlorococcus/papers_and_supp/coe 2024/Supplemental Table 3.csv"
```
Identify columns: `gene`, `dark-tolerant_cluster`, `parental_cluster`, `same_periodicity_pattern`. Count unique clusters in each column.

- [ ] **Step 2: Read existing paperconfig**

Read `data/Prochlorococcus/papers_and_supp/coe 2024/paperconfig.yaml` to understand the current structure (2 experiments, 2 supp tables with DE analyses).

- [ ] **Step 3: Add gene_clusters entry**

Append a `gene_clusters` entry to the existing `supplementary_materials` section. Two possible clusterings:
- dark-tolerant clusters (`dark-tolerant_cluster` column)
- parental clusters (`parental_cluster` column)

Decide whether to create one or two `gene_clusters` entries (one per clustering). The `gene` column is the gene ID.

- [ ] **Step 4: Create README.md**

Summarize the full paper: existing DE analyses + new clustering data.

- [ ] **Step 5: Present for user approval**

- [ ] **Step 6: Commit**

```bash
git add "data/Prochlorococcus/papers_and_supp/coe 2024/paperconfig.yaml" \
        "data/Prochlorococcus/papers_and_supp/coe 2024/README.md"
git commit -m "feat: add gene_clusters entries to coe 2024 paperconfig"
```

---

### Task 1.5: lindell 2007 — Check if paperconfig needs update

**Files:**
- Possibly modify: `data/Prochlorococcus/papers_and_supp/lindell 2007/paperconfig.yaml`
- Create: `data/Prochlorococcus/papers_and_supp/lindell 2007/README.md`

- [ ] **Step 1: Check what changed in the commit**

```bash
git diff 7652261..8499e3e -- "data/Prochlorococcus/papers_and_supp/lindell 2007/"
```
See if `supp table 3.csv` content changed (rows added/removed/modified) or just formatting.

- [ ] **Step 2: Read existing paperconfig and check if it references supp table 3**

```bash
cat "data/Prochlorococcus/papers_and_supp/lindell 2007/paperconfig.yaml"
```

- [ ] **Step 3: Update paperconfig if needed, create README.md**

- [ ] **Step 4: Present for user approval**

- [ ] **Step 5: Commit if changes were made**

---

### Task 1.6: Biller 2018, Wang 2014, Waldbauer 2012, Hennon 2017, tetu 2019 — READMEs for skipped/unchanged papers

**Files:**
- Create: `data/Prochlorococcus/papers_and_supp/Biller 2018/README.md`
- Create: `data/Prochlorococcus/papers_and_supp/Wang 2014/README.md`
- Create: `data/Prochlorococcus/papers_and_supp/Waldbauer  2012/README.md`

- [ ] **Step 1: Explore each paper's data to write informed READMEs**

For each: read CSV headers, legends, check what data is available and why it's being skipped or not updated.

- [ ] **Step 2: Create README.md for each**

- Biller 2018: existing paperconfig covers DE data; new S4A/B (periodicity) and S5 (darkness categories) skipped — not standard DE or clustering
- Wang 2014: gene classification + operon predictions — not DE or clustering; skip
- Waldbauer 2012: PDF-only supplementary tables — no machine-readable data
- Hennon 2017 and tetu 2019: no CSV changes in commit, paperconfigs already complete — note in README if one doesn't exist yet

- [ ] **Step 3: Present all READMEs for user approval**

- [ ] **Step 4: Commit**

```bash
git add "data/Prochlorococcus/papers_and_supp/Biller 2018/README.md" \
        "data/Prochlorococcus/papers_and_supp/Wang 2014/README.md" \
        "data/Prochlorococcus/papers_and_supp/Waldbauer  2012/README.md"
git commit -m "docs: add READMEs for Biller 2018, Wang 2014, Waldbauer 2012"
```

---

## Phase 2: Infrastructure Code Changes + Synechococcus Paperconfigs

### Task 2.1: Add `fold_change_type` support to omics adapter (TDD)

**Files:**
- Modify: `multiomics_kg/adapters/omics_adapter.py:657-691`
- Modify: `.claude/skills/paperconfig/validate_paperconfig.py`
- Create: `tests/test_omics_adapter_fold_change_type.py`

- [ ] **Step 1: Write failing test for linear fold-change conversion**

Create `tests/test_omics_adapter_fold_change_type.py`:

```python
"""Tests for fold_change_type support in omics adapter."""
import math
import pytest

# We'll test the conversion logic in isolation first,
# then integration with a minimal paperconfig.


def test_linear_fc_converted_to_log2():
    """Linear fold-change 2.0 should become log2(2.0) = 1.0."""
    from multiomics_kg.adapters.omics_adapter import _convert_fold_change
    assert _convert_fold_change(2.0, "linear") == pytest.approx(1.0)


def test_linear_fc_less_than_one():
    """Linear fold-change 0.5 should become log2(0.5) = -1.0."""
    from multiomics_kg.adapters.omics_adapter import _convert_fold_change
    assert _convert_fold_change(0.5, "linear") == pytest.approx(-1.0)


def test_linear_fc_of_one():
    """Linear fold-change 1.0 (no change) should become 0.0."""
    from multiomics_kg.adapters.omics_adapter import _convert_fold_change
    assert _convert_fold_change(1.0, "linear") == pytest.approx(0.0)


def test_linear_fc_zero_returns_none():
    """Linear fold-change 0 is undefined in log2 -- return None."""
    from multiomics_kg.adapters.omics_adapter import _convert_fold_change
    assert _convert_fold_change(0.0, "linear") is None


def test_linear_fc_negative_returns_none():
    """Negative linear fold-change is invalid -- return None."""
    from multiomics_kg.adapters.omics_adapter import _convert_fold_change
    assert _convert_fold_change(-1.5, "linear") is None


def test_log2_fc_passthrough():
    """log2 fold-change should pass through unchanged."""
    from multiomics_kg.adapters.omics_adapter import _convert_fold_change
    assert _convert_fold_change(-2.5, "log2") == pytest.approx(-2.5)


def test_default_fc_type_is_log2():
    """When fold_change_type is None/omitted, treat as log2."""
    from multiomics_kg.adapters.omics_adapter import _convert_fold_change
    assert _convert_fold_change(1.5, None) == pytest.approx(1.5)
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/test_omics_adapter_fold_change_type.py -v
```
Expected: FAIL with `ImportError: cannot import name '_convert_fold_change'`

- [ ] **Step 3: Implement `_convert_fold_change` function**

Add to `multiomics_kg/adapters/omics_adapter.py` before the `OMICSAdapter` class (around line 80):

```python
def _convert_fold_change(fc_float: float, fold_change_type: str | None) -> float | None:
    """Convert fold-change value to log2 scale.

    Args:
        fc_float: The raw fold-change value from the CSV.
        fold_change_type: 'log2' (default), 'linear', or None (treated as log2).

    Returns:
        log2-scale fold-change, or None if conversion is invalid (e.g., linear FC <= 0).
    """
    if fold_change_type is None or fold_change_type == "log2":
        return fc_float
    if fold_change_type == "linear":
        if fc_float <= 0:
            return None
        return math.log2(fc_float)
    return fc_float
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/test_omics_adapter_fold_change_type.py -v
```
Expected: All PASS

- [ ] **Step 5: Integrate into the adapter's edge-generation loop**

In `omics_adapter.py`, find the block at line ~689 where `fc_float` is stored:

```python
                if fc_float is not None:
                    edge_properties['log2_fold_change'] = fc_float
```

Replace with:

```python
                if fc_float is not None:
                    fold_change_type = analysis.get('fold_change_type', None)
                    converted = _convert_fold_change(fc_float, fold_change_type)
                    if converted is None:
                        logger.debug(
                            "Skipping row with invalid linear FC %.4f in %s",
                            fc_float, filename
                        )
                        skipped_count += 1
                        continue
                    edge_properties['log2_fold_change'] = converted
```

- [ ] **Step 6: Add `fold_change_type` to paperconfig validation**

In `.claude/skills/paperconfig/validate_paperconfig.py`, add `"fold_change_type"` to the set of known optional analysis fields (find where other optional fields like `pvalue_asterisk_in_logfc`, `timepoint`, `timepoint_hours` are listed). Add validation that value must be `"log2"` or `"linear"` if present.

- [ ] **Step 7: Run full test suite**

```bash
pytest -m "not slow and not kg" -v
```
Expected: All existing tests pass, new tests pass.

- [ ] **Step 8: Commit**

```bash
git add multiomics_kg/adapters/omics_adapter.py \
        tests/test_omics_adapter_fold_change_type.py \
        .claude/skills/paperconfig/validate_paperconfig.py
git commit -m "feat: add fold_change_type support (linear FC → log2 conversion)"
```

---

### Task 2.2: Add multi-paperconfig-list support to create_knowledge_graph.py

**Files:**
- Modify: `create_knowledge_graph.py:80-86`
- Create: `data/Synechococcus/papers_and_supp/paperconfig_files.txt`

- [ ] **Step 1: Create empty Synechococcus paperconfig_files.txt**

```bash
touch "data/Synechococcus/papers_and_supp/paperconfig_files.txt"
```

- [ ] **Step 2: Update create_knowledge_graph.py to read both list files**

Find the `MultiOMICSAdapter` instantiation at line 80:

```python
    omics_adapter = MultiOMICSAdapter(
        config_list_file='data/Prochlorococcus/papers_and_supp/paperconfig_files.txt',
        test_mode=TEST_MODE,
    )
```

Check if `MultiOMICSAdapter` accepts a list of config_list_files or only a single one. Read the `MultiOMICSAdapter.__init__` to see the parameter type.

If it only accepts a single file, the simplest approach is to change `config_list_file` to accept a list of paths:

Option A — change the parameter to accept a list:
```python
    omics_adapter = MultiOMICSAdapter(
        config_list_file=[
            'data/Prochlorococcus/papers_and_supp/paperconfig_files.txt',
            'data/Synechococcus/papers_and_supp/paperconfig_files.txt',
        ],
        test_mode=TEST_MODE,
    )
```

Option B — if the adapter doesn't support lists, modify `MultiOMICSAdapter.__init__` to accept `str | list[str]` and concatenate the file lists internally.

Read the `MultiOMICSAdapter` class to determine which approach is needed and implement it.

- [ ] **Step 3: Run tests to verify no regressions**

```bash
pytest -m "not slow and not kg" -v
```

- [ ] **Step 4: Commit**

```bash
git add create_knowledge_graph.py \
        "data/Synechococcus/papers_and_supp/paperconfig_files.txt" \
        multiomics_kg/adapters/omics_adapter.py
git commit -m "feat: support multiple paperconfig list files (Prochlorococcus + Synechococcus)"
```

---

### Task 2.3: Tal 2009 — New paperconfig (microarray, WH8102)

This paper is the simplest Synechococcus paper (WH8102 already in KG, log2 FC, SYNW#### locus tags).

**Files:**
- Create: `data/Synechococcus/papers_and_supp/Tal 2009/paperconfig.yaml`
- Create: `data/Synechococcus/papers_and_supp/Tal 2009/README.md`

- [ ] **Step 1: Explore data**

```bash
head -10 "data/Synechococcus/papers_and_supp/Tal 2009/table s1.csv"
head -10 "data/Synechococcus/papers_and_supp/Tal 2009/table s2.csv"
```
Verify: Gene ID column (`Gene ID` with SYNW#### values), `Log2 fold change` column, `Score` column. Tables s1 = upregulated, s2 = downregulated (both contain log2 FC).

- [ ] **Step 2: Read paper PDF for methods context**

Check for PDF, read methods for: organism strain (WH8102), coculture partner (Vibrio parahaemolyticus — which strain?), treatment conditions, statistical approach.

- [ ] **Step 3: Create paperconfig.yaml**

Key fields:
- `omics_type`: MICROARRAY
- `organism`: "Synechococcus WH8102" (already in KG but may need to check preferred_name format)
- `treatment_organism`: Vibrio parahaemolyticus (strain from paper)
- `name_col`: `Gene ID`
- `logfc_col`: `Log2 fold change`
- `adjusted_p_value_col`: omit (no p-values; Score is not a p-value)
- `table_scope`: `significant_only` (tables contain only up/down regulated genes)
- Two supp tables pointing to table s1 and s2, or one merged — depends on whether they can share an analysis ID

- [ ] **Step 4: Create README.md**

- [ ] **Step 5: Present for user approval**

- [ ] **Step 6: Add to Synechococcus paperconfig_files.txt and commit**

```bash
echo "data/Synechococcus/papers_and_supp/Tal 2009/paperconfig.yaml" >> "data/Synechococcus/papers_and_supp/paperconfig_files.txt"
git add "data/Synechococcus/papers_and_supp/Tal 2009/paperconfig.yaml" \
        "data/Synechococcus/papers_and_supp/Tal 2009/README.md" \
        "data/Synechococcus/papers_and_supp/paperconfig_files.txt"
git commit -m "feat: add Tal 2009 paperconfig (WH8102 Vibrio coculture microarray)"
```

---

### Task 2.4: Beliaev 2014 — New paperconfig (Syn PCC 7002 + Shewanella, linear FC)

**Files:**
- Create: `data/Synechococcus/papers_and_supp/Beliaev 2014/paperconfig.yaml`
- Create: `data/Synechococcus/papers_and_supp/Beliaev 2014/README.md`

- [ ] **Step 1: Explore all 8 CSV tables**

Read headers and sample rows from each table. Determine:
- Which tables are full gene lists vs filtered DE genes
- What the fold-change columns represent (which comparisons)
- Confirm linear FC (values > 0, no negative values except possibly as ratios < 1)

- [ ] **Step 2: Read paper PDF for methods**

Identify: exact Synechococcus PCC 7002 strain, Shewanella W3-18-1 details, experimental conditions (lactate vs HCO3- vs axenic), how tables s1-s8 map to comparisons.

- [ ] **Step 3: Create paperconfig.yaml**

This paper has multiple experiments and both organisms:
- Experiments: lactate_vs_axenic, hco3_vs_axenic, hco3_vs_lactate (for each organism)
- `fold_change_type: linear` on all analyses
- Tables alternate: odd = Synechococcus, even = Shewanella
- `organism`: one per CSV table (Syn PCC 7002 or Shewanella W3-18-1)
- Need to determine which tables to use (avoid duplicate edges from overlapping tables)
- New organisms: both need genome entries (Phase 3)

- [ ] **Step 4: Create README.md**

- [ ] **Step 5: Present for user approval**

- [ ] **Step 6: Add to paperconfig_files.txt and commit**

---

### Task 2.5: Kratzl 2024 — New paperconfig (S. elongatus + P. putida, RNA-seq + proteomics)

**Files:**
- Create: `data/Synechococcus/papers_and_supp/Kratzl 2024/paperconfig.yaml`
- Create: `data/Synechococcus/papers_and_supp/Kratzl 2024/README.md`

- [ ] **Step 1: Explore CSV files and legends**

```bash
head -5 "data/Synechococcus/papers_and_supp/Kratzl 2024/data set s1 6_coculture vs Selongatus_all 42003_2024_6098_MOESM3_ESM.csv"
head -5 "data/Synechococcus/papers_and_supp/Kratzl 2024/data set s1 1_coculture vs Pputida_all 42003_2024_6098_MOESM3_ESM.csv"
cat "data/Synechococcus/papers_and_supp/Kratzl 2024/data set s1 legends  42003_2024_6098_MOESM3_ESM.txt"
```
Identify gene ID columns, strip `gene-` prefix strategy, proteomics protein ID columns.

- [ ] **Step 2: Read paper PDF for methods, verify strains and assemblies**

Confirm: S. elongatus PCC 7942 strain, P. putida KT2440 strain, NCBI accessions used, statistical methods.

- [ ] **Step 3: Create paperconfig.yaml**

- 4 analyses: 2 RNA-seq + 2 proteomics (one per organism per omics type)
- RNA-seq: `fold_change_type: log2`, has adjusted p-value
- Proteomics: `fold_change_type: log2`, has adjusted p-value
- Gene IDs: RNA-seq uses `gene-PREFIX_RS#####` (need `gene-` stripped — can use `id_columns` or pre-process); proteomics uses UniProt accessions
- Note: `gene-` prefix stripping may need a custom `name_col` approach or pre-processing the CSV

- [ ] **Step 4: Create README.md**

- [ ] **Step 5: Present for user approval**

- [ ] **Step 6: Commit**

---

### Task 2.6: Ma 2022 — New paperconfig (S. elongatus + E. coli)

**Files:**
- Create: `data/Synechococcus/papers_and_supp/Ma 2022/paperconfig.yaml`
- Create: `data/Synechococcus/papers_and_supp/Ma 2022/README.md`

- [ ] **Step 1: Explore CSV files**

Read all 3 CSVs. Verify:
- S1 and S2 are split by direction (up/down) — both have log2FC
- S4 is proteomics with Mean Ratio (linear)
- Gene ID format: `M744_#####`

- [ ] **Step 2: Read paper PDF, identify exact strain and E. coli strain**

Determine: is this PCC 7942 or UTEX 2973? Which E. coli strain? What is `M744_` prefix — look up on NCBI.

- [ ] **Step 3: Create paperconfig.yaml**

- S1 + S2: can use two analyses (one per file) with same experiment, or merge CSVs
- S4: proteomics with `fold_change_type: linear`
- E. coli: treatment organism only (no DE data for E. coli side)

- [ ] **Step 4: Create README.md**

- [ ] **Step 5: Present for user approval**

- [ ] **Step 6: Commit**

---

### Task 2.7: Oleza 2015 — New paperconfig (WH7803 + R. pomeroyi exoproteome)

**Files:**
- Create: `data/Synechococcus/papers_and_supp/Oleza 2015/paperconfig.yaml`
- Create: `data/Synechococcus/papers_and_supp/Oleza 2015/README.md`

- [ ] **Step 1: Check if CSVs exist**

This paper may only have xlsx files. If no CSVs exist, need to export from xlsx first or skip.

```bash
ls "data/Synechococcus/papers_and_supp/Oleza 2015/"
```

- [ ] **Step 2: If xlsx only, determine if CSV export is needed**

If only xlsx: either export relevant sheets to CSV (as a prep step) or skip this paper until CSVs are prepared. Present decision to user.

- [ ] **Step 3: Create paperconfig.yaml (if CSVs available)**

- Table s2a: Syn WH7803 exoproteome, `fold_change_type: linear`, has p-value
- Table s2b: R. pomeroyi — abundance only, skip for DE
- Gene IDs: RefSeq protein accessions (YP_) — will need protein-level resolution (Phase 4)

- [ ] **Step 4: Create README.md**

- [ ] **Step 5: Present for user approval**

- [ ] **Step 6: Commit**

---

### Task 2.8: Oleza 2017 — New paperconfig (WH7803 + R. pomeroyi long-term)

**Files:**
- Create: `data/Synechococcus/papers_and_supp/Oleza 2017/paperconfig.yaml`
- Create: `data/Synechococcus/papers_and_supp/Oleza 2017/README.md`

- [ ] **Step 1: Explore CSV files**

5 CSVs available. Read headers, identify which have fold-change + p-values vs abundance-only.

- [ ] **Step 2: Read paper PDF for experiment structure (early vs late coculture phases)**

- [ ] **Step 3: Create paperconfig.yaml**

- Tables s2a, s5a: Syn WH7803 proteomics (early + late), `fold_change_type: linear`
- Table s7: R. pomeroyi long-term proteomics, `fold_change_type: linear`, has q-value
- Tables s2b, s5b: R. pomeroyi abundance-only — skip for DE
- Multiple experiments: early_coculture, late_coculture, long_term

- [ ] **Step 4: Create README.md**

- [ ] **Step 5: Present for user approval**

- [ ] **Step 6: Commit**

---

### Task 2.9: kaur 2018 — New paperconfig (WH7803 + R. pomeroyi exoproteome time-series)

**Files:**
- Create: `data/Synechococcus/papers_and_supp/kaur 2018/paperconfig.yaml`
- Create: `data/Synechococcus/papers_and_supp/kaur 2018/README.md`

- [ ] **Step 1: Explore CSV files and legends**

```bash
head -5 "data/Synechococcus/papers_and_supp/kaur 2018/table s4a Syn SW emi14012-sup-0004-suppinfo4.csv"
cat "data/Synechococcus/papers_and_supp/kaur 2018/legends.txt"
```
Key question: are the fold-change values linear or log2? Check for negative values and range.

- [ ] **Step 2: Read paper PDF, verify FC type and experimental design**

Timepoint comparisons (TP3 vs 1, TP7 vs 1, etc.) — this is a time-series proteomics experiment. Determine the time-course structure.

- [ ] **Step 3: Create paperconfig.yaml**

- 4 tables: s4a/b (Syn WH7803 in SW/ASW), s5a/b (R. pomeroyi in SW/ASW)
- Multiple fold-change columns per table (TP3 vs 1, TP7 vs 1) — each becomes a separate timepoint analysis
- Has q-values for significance
- `fold_change_type`: `linear` or `log2` (determined in step 1)

- [ ] **Step 4: Create README.md**

- [ ] **Step 5: Present for user approval**

- [ ] **Step 6: Commit**

---

### Task 2.10: Bernstein 2017 — New paperconfig (T. elongatus clustering only)

**Files:**
- Create: `data/Synechococcus/papers_and_supp/Bernstein 2017/paperconfig.yaml`
- Create: `data/Synechococcus/papers_and_supp/Bernstein 2017/README.md`

- [ ] **Step 1: Explore cluster data**

```bash
head -10 "data/Synechococcus/papers_and_supp/Bernstein 2017/data set s2 clustered_Avg_Cond_Light sys002172092sd2.csv"
head -10 "data/Synechococcus/papers_and_supp/Bernstein 2017/data set s2 clusterd_Avg_Cond_Oxygen sys002172092sd2.csv"
cat "data/Synechococcus/papers_and_supp/Bernstein 2017/legends.txt"
```
Identify: cluster columns (`ClustID_light`, `ClustID_ox`), gene ID column (NCBI Locus Tag = `SY28_RS#####`), unique cluster count.

- [ ] **Step 2: Create paperconfig.yaml**

Clustering-only paperconfig (no DE data):
- Two `gene_clusters` entries: one for light clustering, one for oxygen clustering
- `gene_id_col`: `NCBI Locus Tag`
- `cluster_col`: `ClustID_light` or `ClustID_ox`
- Coculture partner: Meiothermus ruber (treatment organism only)
- New organism needed: T. elongatus BP-1 (verify NCBI name/accession)

- [ ] **Step 3: Create README.md**

Note: data set s1 has raw expression values (not integrated), data set s2 has clustering only.

- [ ] **Step 4: Present for user approval**

- [ ] **Step 5: Commit**

---

### Task 2.11: Zhang 2021 — README only (skip)

**Files:**
- Create: `data/Synechococcus/papers_and_supp/Zhang 2021/README.md`

- [ ] **Step 1: Create README.md**

Note: doc supplement only, no machine-readable expression data. Skip for KG integration.

- [ ] **Step 2: Present for user approval**

- [ ] **Step 3: Commit**

---

## Phase 3: Genome Setup + Data Download

### Task 3.1: Look up NCBI accessions for all new organisms

- [ ] **Step 1: For each new organism, search NCBI for RefSeq assembly accession**

Use the organism names and strain identifiers confirmed during Phase 2 paperconfig creation. For each:
- Search NCBI Assembly database
- Prefer RefSeq (GCF_) accessions; use GenBank (GCA_) if no RefSeq
- Record: accession, taxid, genome size, assembly level
- Cross-reference with Cyanorak tables (`data/Cyanorak Organism Table synechococcus.csv`) for marine Synechococcus strains

Organisms to look up:
- Synechococcus sp. WH7803 (Cyanorak has GenBank CT971583, RefSeq NC_009481.1 — need GCF accession)
- Synechococcus sp. PCC 7002
- Thermosynechococcus elongatus BP-1
- Synechococcus elongatus PCC 7942 (or UTEX 2973 — resolved in Phase 2)
- Prochlorococcus MIT9303
- Shewanella putrefaciens W3-18-1
- Pseudomonas putida KT2440
- Ruegeria pomeroyi DSS-3

- [ ] **Step 2: Present organism inventory with accessions for user approval**

---

### Task 3.2: Add new organisms to cyanobacteria_genomes.csv

**Files:**
- Modify: `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv`

- [ ] **Step 1: Add rows for each new genome organism**

For each organism with DE data, add a row with:
`ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir,clade,preferred_name`

- Marine Synechococcus: populate `cyanorak_organism` from the Cyanorak table (e.g., `Syn_WH7803`)
- Freshwater/other cyanobacteria: leave `cyanorak_organism` blank
- Heterotrophs with DE data: leave `cyanorak_organism` blank
- `data_dir`: follow pattern `cache/data/<OrgGroup>/genomes/<Strain>/`

- [ ] **Step 2: Commit**

```bash
git add data/Prochlorococcus/genomes/cyanobacteria_genomes.csv
git commit -m "feat: add new organism genomes (Synechococcus + heterotroph strains)"
```

---

### Task 3.3: Add treatment organisms

**Files:**
- Modify: `data/Prochlorococcus/treatment_organisms.csv`

- [ ] **Step 1: Add rows for heterotrophs without DE data**

Add: Meiothermus ruber, E. coli (strain), Vibrio parahaemolyticus (strain). Format:
`ncbi_taxon_id,organism_name`

- [ ] **Step 2: Commit**

```bash
git add data/Prochlorococcus/treatment_organisms.csv
git commit -m "feat: add treatment organisms (Meiothermus, E. coli, Vibrio)"
```

---

### Task 3.4: Run data download pipeline

- [ ] **Step 1: Run prepare_data.sh step 0 for new strains**

```bash
bash scripts/prepare_data.sh --steps 0 --strains WH7803 PCC7002 BP1 PCC7942 MIT9303 W3-18-1 KT2440 DSS-3
```
(Adjust strain names to match `strain_name` in cyanobacteria_genomes.csv)

This downloads: NCBI genome files, Cyanorak data (for marine Syn only), UniProt data per taxid.

- [ ] **Step 2: Verify downloads completed**

Check that `cache/data/<OrgGroup>/genomes/<Strain>/` directories contain expected files (GFF, protein FASTA, GBFF).

- [ ] **Step 3: Run steps 1-2 (build annotations without eggNOG)**

```bash
bash scripts/prepare_data.sh --steps 1 2 --strains WH7803 PCC7002 BP1 PCC7942 MIT9303 W3-18-1 KT2440 DSS-3
```

- [ ] **Step 4: Verify gene_annotations_merged.json created for each strain**

---

## Phase 4: Gene ID Resolution

### Task 4.1: Build gene ID mappings and resolve paper IDs

- [ ] **Step 1: Run steps 3-4 for ALL strains**

```bash
bash scripts/prepare_data.sh --steps 3 4
```
This rebuilds `gene_id_mapping.json` for all strains (including existing ones that may pick up new `id_columns` from new paperconfigs) and resolves paper CSVs.

- [ ] **Step 2: Check gene ID match rates**

Run `/check-gene-ids` for every paper with expression data. Focus on:
- New Synechococcus papers
- Capovilla 2023 (MIT9303/MIT9313)
- Any existing papers that may be affected

- [ ] **Step 3: Fix mismatches iteratively**

For papers with low match rates:
- Add `id_translation` entries to paperconfigs for non-standard IDs
- Add `annotation_gff` entries where GFF bridges are needed
- For protein-level IDs (kaur 2018, Oleza, Kratzl proteomics): add `id_columns` with `protein_id` type
- Re-run steps 3-4 after each fix
- Present match rate improvements to user

- [ ] **Step 4: Commit all paperconfig updates**

```bash
git add data/*/papers_and_supp/*/paperconfig.yaml
git commit -m "fix: gene ID resolution fixes for new papers"
```

---

## Phase 5: Build and Validate

### Task 5.1: Build KG and validate

- [ ] **Step 1: Take edge snapshot**

Run `/omics-edge-snapshot` to capture current state before rebuild.

- [ ] **Step 2: Build KG**

```bash
uv run python create_knowledge_graph.py
```

- [ ] **Step 3: Deploy to Docker**

```bash
docker compose up -d
```

- [ ] **Step 4: Compare edge snapshot**

Run `/omics-edge-snapshot` again and compare. Verify:
- New edges from new papers
- No regressions on existing papers (same or higher edge counts)

- [ ] **Step 5: Run KG validity tests**

```bash
pytest -m kg -v
```

- [ ] **Step 6: Update test expectations**

Update `tests/kg_validity/` test files for new expected counts:
- Organism count (was 15)
- Gene count (was ~35K)
- Expression edge count (was ~188K)
- Experiment count (was 76)

- [ ] **Step 7: Regenerate snapshot fixtures**

```bash
uv run python tests/kg_validity/generate_snapshot.py
```

- [ ] **Step 8: Final commit**

```bash
git add tests/kg_validity/
git commit -m "test: update KG validity expectations for new organisms and papers"
```
