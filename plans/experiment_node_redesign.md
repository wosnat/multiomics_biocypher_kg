# Plan: Experiment Node Redesign — KG-side changes

**Parent plan:** `multiomics_explorer/plans/redefine_mcp_tools/expression_tools_redesign.md`
**Scope:** Steps 1–4 of Phase 1 (KG repo only). Paperconfig audit, schema changes,
adapter rewrite, post-import updates, rebuild.

**Goal:** Replace the flat `Factor → Gene` expression edge model with a structured
`Publication → Experiment → Gene` model. Introduce `Experiment` nodes that group
time-series analyses, move metadata off edges and onto nodes, and unify coculture
and environmental expression edges into a single `Changes_expression_of` edge type
(source: Experiment, target: Gene).

---

## Current state

- **~188K expression edges** split across two types:
  - `Condition_changes_expression_of` (~171K): EnvironmentalCondition → Gene
  - `Coculture_changes_expression_of` (~17K): OrganismTaxon → Gene
- **~40 EnvironmentalCondition nodes** (pub-scoped: `{doi}_{env_id}`)
- **~21 Publication nodes** with `Published_expression_data_about` edges to sources
- **26 active paperconfigs** in `paperconfig_files.txt`
- **~180 statistical_analyses entries** across all papers
- Edge properties carry duplicated metadata: organism_strain, treatment_condition,
  control_condition, experimental_context, omics_type, statistical_test, analysis_name,
  publications (all repeated per gene per analysis)
- `timepoint` field exists on many analyses but is only passed through as edge
  `time_point` string property — no grouping, no ordering, no numeric normalization
- Schema already has unused `timeseries_cluster` / `cluster_in_publication` /
  `molecular_in_cluster` definitions (will be replaced in Phase 2)

## Target state

- **New `Experiment` node** (~100–120 nodes) grouping analyses by scientific question
- **Single `Changes_expression_of` edge type** (Experiment → Gene) — unifies the two old
  edge labels, with minimal per-gene properties (fold change, p-value, timepoint)
- **New structural edges:** `HAS_EXPERIMENT`, `TESTS_COCULTURE_WITH`
- **Removed:** `Condition_changes_expression_of`, `Coculture_changes_expression_of`,
  `Published_expression_data_about`, `EnvironmentalCondition` nodes, `TESTS_CONDITION`
- **Publication nodes:** retained as-is
- All paperconfigs restructured: new `experiments` block defines experiment-level
  metadata; analyses reference experiments and carry only per-timepoint/CSV fields

## Why EnvironmentalCondition is absorbed into Experiment

The current graph has separate `EnvironmentalCondition` nodes (~40) with structured
properties (medium, temperature, light_condition, etc.) linked to expression edges.
The redesign **merges these into the Experiment node** rather than keeping them as
separate nodes. Rationale:

1. **Already pub-scoped, no reuse:** EnvironmentalCondition nodes are keyed as
   `{doi}_{env_id}` — each paper creates its own. There's no cross-paper sharing or
   normalization benefit from keeping them as separate nodes.

2. **Structured properties are sparsely populated and rarely queried:** Of ~40 condition
   nodes, most only have `treatment_type` and a text description. The detailed fields
   (nitrogen_level, phosphate_level, light_intensity, etc.) are used by a few papers.
   The MCP tools would almost never traverse `Experiment → TESTS_CONDITION →
   EnvironmentalCondition` to read medium or temperature.

3. **treatment_type already denormalized:** The key filterable property (`treatment_type`)
   is going onto the Experiment node regardless. Medium and temperature (the next most
   useful) can go directly onto Experiment too — they're almost always shared between
   treatment and control (the variable is one factor, not the growth conditions).

4. **Eliminates a node type + edge type:** No `EnvironmentalCondition` nodes, no
   `TESTS_CONDITION` edges, simpler graph, simpler Cypher.

**What moves to Experiment:**
- `treatment_type` — from the condition's `treatment_type` field (or `"coculture"`)
- `medium` — from the environmental_conditions block (shared between treatment/control)
- `temperature` — from the environmental_conditions block (shared)
- `light_condition` — light regime: continuous, diel cycle, darkness (critical for cyanobacteria)
- `light_intensity` — PAR in µmol photons m⁻² s⁻¹

**Population note:** Most current paperconfigs do not have `environmental_conditions`
blocks with `light_condition` or `light_intensity` populated. The migration script should
set these to the sentinel value `"NEEDS_CURATION"` when the source data is missing, so
that incomplete fields can be found later with: `grep -r "NEEDS_CURATION" data/`

**What happens to the paperconfig `environmental_conditions` block:**
Replaced by the new `experiments` block (see Step 1). The `experiments` block defines
experiment-level metadata directly — including `treatment_type`, `medium`, `temperature`
— instead of indirecting through a separate conditions block. The
`environmental_conditions` block and the `environmental_*_condition_id` reference fields
on analyses are removed.

**Graph structure (simplified):**

```
Publication ──HAS_EXPERIMENT──> Experiment ──CHANGES_EXPRESSION_OF──> Gene
                                    │
                                    └──TESTS_COCULTURE_WITH──> OrganismTaxon
                                       (coculture experiments only)
```

---

## Step 1: Paperconfig audit and cleanup

**What:** Restructure all 26 paperconfigs: add a new `experiments` block as a
first-class section, move shared experiment-level fields out of individual analyses,
remove `environmental_conditions` block (absorbed into experiments). Automated by
a one-time migration script.

**Files:**
- All 26 `paperconfig.yaml` files under `data/Prochlorococcus/papers_and_supp/`
- New: `scripts/migrate_paperconfigs.py` (one-time migration script)
- Updated: `.claude/skills/paperconfig/` (new format for future papers)
- Updated: `multiomics_kg/utils/validate_paperconfig.py` (new format validation)

### 1a. New paperconfig format

**Before (current):**
```yaml
publication:
  papername: "Weissberg 2025"
  papermainpdf: "data/.../paper.pdf"

  environmental_conditions:
    weissberg_pro99_exponential:
      condition_type: growth_state
      name: "PRO99-lowN exponential growth"
      medium: "PRO99-lowN"
      temperature: "24C"
    weissberg_pro99_starvation:
      condition_type: growth_state
      name: "PRO99-lowN nutrient starvation"
      medium: "PRO99-lowN"
      temperature: "24C"

  supplementary_materials:
    hot1a3_rnaseq_axenic_day_18:
      type: csv
      filename: "data/.../day_18.csv"
      organism: "Alteromonas macleodii HOT1A3"
      statistical_analyses:
        - type: RNASEQ
          id: DE_hot1a3_axenic_d18
          name: "DE of HOT1A3 axenic day 18 vs day 11"
          test_type: DESeq2
          organism: "Alteromonas macleodii HOT1A3"          # repeated
          control_condition: "Exponential growth (day 11)"  # repeated
          treatment_condition: "Nutrient starvation (day 18)"
          experimental_context: "Axenic in PRO99-lowN at 24C"  # repeated
          environmental_control_condition_id: weissberg_pro99_exponential
          environmental_treatment_condition_id: weissberg_pro99_starvation
          timepoint: "day 18"
          name_col: "gene_id"
          logfc_col: "logFC"
          adjusted_p_value_col: "padj"
          pvalue_threshold: 0.05
```

**After (new):**
```yaml
publication:
  papername: "Weissberg 2025"
  papermainpdf: "data/.../paper.pdf"

  experiments:
    n_starvation_axenic_hot1a3_rnaseq:
      name: "HOT1A3 axenic nutrient starvation time course"
      organism: "Alteromonas macleodii HOT1A3"
      treatment_condition: "Nutrient starvation"
      control_condition: "Exponential growth (day 11)"
      experimental_context: "Axenic in PRO99-lowN at 24C under continuous light"
      omics_type: RNASEQ
      test_type: DESeq2
      treatment_type: growth_state
      medium: "PRO99-lowN"
      temperature: "24C"
      light_condition: "continuous light"
      light_intensity: ""  # not reported in paper
      # coculture fields (only for coculture experiments):
      # treatment_organism: "Alteromonas macleodii HOT1A3"
      # treatment_taxid: 28108
      # treatment_assembly_accession: "GCF_001077695.1"

  supplementary_materials:
    hot1a3_rnaseq_axenic_day_18:
      type: csv
      filename: "data/.../day_18.csv"
      organism: "Alteromonas macleodii HOT1A3"
      statistical_analyses:
        - id: DE_hot1a3_axenic_d18
          experiment: n_starvation_axenic_hot1a3_rnaseq  # reference to experiments block
          timepoint: "day 18"
          timepoint_hours: 432
          name_col: "gene_id"
          logfc_col: "logFC"
          adjusted_p_value_col: "padj"
          pvalue_threshold: 0.05
```

### 1b. What moves where

**From `environmental_conditions` → `experiments` block:**
- `condition_type` → `treatment_type` (renamed)
- `medium`, `temperature` → directly on experiment
- `light_condition`, `light_intensity` → directly on experiment

**From each `statistical_analyses` entry → `experiments` block (shared, defined once):**
- `organism`
- `type` → `omics_type`
- `test_type`
- `treatment_condition` (time-stripped)
- `control_condition`
- `experimental_context`
- `name` (experiment-level name, not per-timepoint)
- `treatment_organism`, `treatment_taxid`, `treatment_assembly_accession` (coculture)

**Stays on each `statistical_analyses` entry (per-timepoint/per-CSV):**
- `id` — unique analysis ID
- `experiment` — reference to experiments block key
- `timepoint` — original string label
- `timepoint_hours` — numeric hours (new)
- `name_col`, `logfc_col`, `adjusted_p_value_col` — CSV column names
- `pvalue_threshold`, `logfc_threshold` — significance thresholds
- `significance_mode` — if present
- `id_columns` — gene ID column config (if present)

**Removed entirely:**
- `environmental_conditions` block
- `environmental_control_condition_id` / `environmental_treatment_condition_id` references
- `type` on analyses (moved to experiment as `omics_type`)
- `name` on analyses (moved to experiment — experiment-level name, not per-timepoint)
- Repeated `organism`, `control_condition`, `treatment_condition`, `experimental_context`,
  `test_type` on each analysis

**Note on `analysis_name` edge property:** The current adapter emits `analysis_name`
(from the per-analysis `name` field) on expression edges. This property is removed in
the redesign — the experiment node's `name` property serves this purpose. MCP tools
that display `analysis_name` must switch to reading from the Experiment node.

### 1c. Experiment groups — what defines a group

Each experiment is **scoped within a single publication**. It represents one
biological question: same {organism, treatment (time-stripped), control,
experimental_context, omics_type}. Analyses within a group differ only by timepoint.

**Naming convention:** `{condition}_{organism_short}_{omics}` as a short slug.
Only needs to be unique within the publication.

Examples:
- `n_starvation_axenic_hot1a3_rnaseq`
- `coculture_hot1a3_med4_rnaseq`
- `phage_infection_med4_microarray`
- `p_deprivation_med4_microarray`

For single-point papers, each analysis is its own experiment (one entry in the
experiments block per analysis).

### 1d. Migration script (`scripts/migrate_paperconfigs.py`)

One-time Python script to automate the conversion. Run once, eyeball output, commit,
delete the script.

**What it does:**
1. Reads each paperconfig in old format
2. For each paper, groups `statistical_analyses` by
   {organism, treatment_condition (time-stripped), control_condition,
   experimental_context, type (omics_type), test_type}
3. Generates experiment ID slugs from the grouping key
4. Pulls `treatment_type`, `medium`, `temperature` from `environmental_conditions` block
   (via `environmental_treatment_condition_id` or `environmental_control_condition_id`)
5. For coculture experiments (has `treatment_organism`): sets `treatment_type: coculture`
6. Strips time from `treatment_condition` (Weissberg pattern: `"... (day N)"` → `"..."`)
7. Computes `timepoint_hours` from `timepoint` using normalization rules (see below)
8. Writes new format paperconfig
9. Prints summary: papers processed, experiments created, analyses migrated, any issues

**Validation checks in the script:**
- All analyses in a group agree on shared fields (organism, control, context, etc.)
- No orphan `environmental_conditions` references (all mapped to experiments)
- `treatment_type` is set for every experiment
- Round-trip check: same number of analyses before and after

**Manual review needed after migration:**
- Experiment names (auto-generated from first analysis's `name` field — may need cleanup)
- Weissberg grouping edge cases (phase transitions producing multiple groups)
- Biller 2016 late comparisons (different control → separate experiment)
- Lin 2015 P starvation vs P re-addition (separate experiments)

### 1e. timepoint_hours normalization rules

Used by the migration script and by `validate_paperconfig.py`.

| Pattern | Conversion | Examples |
|---|---|---|
| `"Nh"` or `"N h"` | N | `"4h"` → 4, `"0.5h"` → 0.5 |
| `"-Nh"` | -N | `"-12h"` → -12 |
| `"day N"` / `"Day N"` | N × 24 | `"day 18"` → 432, `"Day 2"` → 48 |
| `"R (rescue: ...)"` | null | Non-numeric |
| `"Nh extended darkness (Mh)"` | M (absolute) | `"1h extended darkness (36h)"` → 36 |
| `"Nh (P added)"` | N | `"50h (P added)"` → 50 |
| `"days N+M"` | null | Pooled samples, no single time |
| `"Nh post-inoculation"` | N | `"0.5h post-inoculation"` → 0.5 |
| `null` / missing | null | Analyses without timepoints (single-point or steady-state) |

**Mixed time-course and single-point papers:** Some papers (e.g., Tolonen 2006) have
time-course analyses alongside single-point steady-state comparisons within the same
supplementary table. Analyses without a `timepoint` field get `timepoint_hours: null`
and belong to a **separate experiment group** (different treatment_condition → different
grouping key). The migration script handles this automatically because the grouping key
includes `treatment_condition`, and the steady-state comparisons (e.g., "Cyanate vs
ammonium") have different treatment text than the time-course (e.g., "N deprivation").

### 1f. Paper inventory

**Time-course papers** (multiple timepoints per experiment):

| Paper | # analyses | estimated # experiments | timepoint range | Notes |
|---|---|---|---|---|
| Weissberg 2025 | ~35 | ~12 | day 11–89 | NEEDS treatment_condition cleanup |
| Coe 2024 | 14 | 2 | 0h–24h | 2 organisms (MIT1002, NATL2A) |
| Lindell 2007 | 8 | 1 | 1h–8h | Phage infection |
| Fang 2019 | 9 | 1 | 0.5h–72h | Phage infection |
| Thompson 2016 | 15 | 4 | 0.5h–8.5h | 4 phage strains × time |
| Thompson 2011 | 8 | 2 | 12h–R | N and Fe starvation + rescue |
| Tolonen 2006 | 15 | 3 | 0h–48h + single-point | N deprivation + alt N sources |
| Martiny 2006 | 9 | 2 | 0h–48h | P starvation: MED4 + MIT9313 |
| Lin 2015 | 14 | 2+ | 4h–59h | P starvation + P re-addition |
| Read 2017 | 3 | 1 | 3h–24h | N removal |
| Biller 2016 | 9 | 2+ | -12h–48h | Coculture time course + late comparisons |
| Biller 2018 | 9 | 3 | 1h–13h extended darkness | Axenic + coculture × strains |
| Ziegler 2025 | 30 | 10 | Day 2–5 | Multiple organisms × coculture partners |
| Aharonovich 2016 | 5 | 5 | 20h (single) | 5 separate pairwise comparisons |

**Single-point papers** (one experiment per analysis):

| Paper | # analyses | Notes |
|---|---|---|
| Al-Hosani 2015 | 1 | Salt acclimation |
| Anjur 2025 | 1 | High light |
| Bagby 2015 | 2 | Gas shock (O2, CO2) |
| Barreto 2022 | 3 | 3 coculture pairs |
| He 2022 | 1 | Coculture |
| Hennon 2017 | 4 | 2 organisms × 2 tables |
| MIT9313_resources | 2 | Light/dark |
| Steglich 2006 | 1 | Light exposure |
| Tetu 2019 | 2 | Phosphonate + methylphosphonate |

### 1g. Weissberg-specific notes

**Experimental design:** Compares late timepoints back to earliest timepoint within
the same culture (e.g., "day 18 vs day 11"). Fold-changes are cumulative from a single
baseline, not independent treatment-vs-control snapshots. Both designs group correctly
(control_condition is constant within a group).

**treatment_condition cleanup:** Strip parenthesized time: `"Nutrient starvation (day 18)"`
→ `"Nutrient starvation"`. Keep control_condition as-is (e.g., "Exponential growth
(day 11)") — the baseline day is important context.

**Phase transitions:** Treatment text encodes biological phases (Nutrient starvation
vs Long-term nutrient starvation vs Decline vs Death/decomposition). After stripping
time, these become separate experiments. This is acceptable — each represents a
distinct biological phase.

### 1h. Acceptance criteria

- [ ] All 26 paperconfigs converted to new format
- [ ] Each paperconfig has an `experiments` block with all experiment definitions
- [ ] Each `statistical_analyses` entry has `experiment` reference to experiments block
- [ ] Each analysis has `timepoint_hours` (or null with comment for unparseable)
- [ ] No `environmental_conditions` block in any paperconfig
- [ ] No `environmental_*_condition_id` references in any analysis
- [ ] No repeated experiment-level fields on individual analyses
- [ ] Weissberg treatment_conditions are time-stripped
- [ ] Every experiment has `treatment_type`
- [ ] `validate_paperconfig.py` updated to validate new format
- [ ] `/paperconfig` skill updated to generate new format
- [ ] Migration script runs cleanly on all 26 papers

### 1i. Full impact analysis — all paperconfig consumers

Every file that reads paperconfig.yaml must be updated for the new format. The changes
fall into three categories: **must change** (reads removed/restructured fields),
**check only** (reads fields that stay on analyses), and **docs/skills** (reference material).

**Must change — code reads `environmental_conditions` or experiment-level fields from analyses:**

| File | What it reads | Required change |
|---|---|---|
| `multiomics_kg/adapters/omics_adapter.py` | `environmental_conditions`, `statistical_analyses` (all fields), `_determine_source_id()` | **Major rewrite** (Step 3). Read `experiments` block for nodes, analyses for edges. Remove env condition node emission. |
| `.claude/skills/paperconfig/validate_paperconfig.py` | All fields — complete schema enforcement | **Rewrite validation**: check `experiments` block, check `experiment` reference on analyses, remove env_conditions validation, add `treatment_type`/`timepoint_hours` validation |
| `tests/test_paperconfig_validation.py` | Creates mock paperconfigs with old format | **Update all fixtures** to new format. Update canonical vocabulary tests. |
| `tests/test_omics_adapter_organism_gene.py` | `sample_config` fixture with `environmental_conditions` and old analysis fields | **Update fixture** to new format with `experiments` block |

**Must change — prepare_data pipeline scripts (steps 3 and 4):**

These run via `bash scripts/prepare_data.sh --steps 3 4` every time a paper is added,
IDs are remapped, or a new strain is deployed. They are **not one-time scripts** — they
must work with the new paperconfig format from the moment the migration is committed.
Breaking them blocks the entire data pipeline.

| File | What it reads | Required change |
|---|---|---|
| `multiomics_kg/download/build_gene_id_mapping.py` (step 3) | `organism` from table or first analysis, `id_columns`, `id_translation`, `annotation_gff` entries | `organism` now comes from experiment (via `analysis.experiment` → `experiments[key].organism`). `id_translation` and `annotation_gff` entries are unchanged (not inside `statistical_analyses`). Must use `get_organism_for_analysis()` helper. |
| `multiomics_kg/download/resolve_paper_ids.py` (step 4) | `organism`, `name_col` from analyses, `environmental_treatment_condition_id` for source detection | `organism` via experiment lookup. Remove `environmental_treatment_condition_id` usage. `name_col` stays on analysis. Must use helper. |

**These must be updated and tested in Commit 2 alongside the paperconfig migration.**
The migrated paperconfigs, new-format helpers in `paperconfig_utils.py`, and all
consumer updates ship together as one atomic commit. If the paperconfigs change but
the pipeline scripts aren't updated, `prepare_data.sh` breaks.

**Must change — skills (run on demand, less critical but still needed):**

| File | What it reads | Required change |
|---|---|---|
| `.claude/skills/check-gene-ids/check_gene_ids.py` | `organism`, `name_col` from analyses | `organism` via experiment lookup |
| `.claude/skills/fix-gene-ids/fix_gene_ids.py` | `organism`, `name_col`, `id_columns` from analyses | `organism` via experiment lookup |
| `.claude/skills/build-gene-mapping-supp/build_gene_mapping_supp.py` | `organism`, `name_col` from analyses | `organism` via experiment lookup |

**Check only — reads fields that stay on analyses (likely no change needed):**

| File | What it reads | Expected impact |
|---|---|---|
| `scripts/validate_annotations.py` | `filename`, `sep`, `skip_rows`, `name_col` | **Minimal** — these stay on analyses/tables |
| `scripts/update_paperconfig_id_columns.py` | Hard-coded mappings, adds `id_columns` | **None** — one-time script, already run |
| `create_knowledge_graph.py` | Passes config file path to MultiOMICSAdapter | **None** — just the path |

**Documentation and skills (reference material):**

| File | Required change |
|---|---|
| `.claude/skills/paperconfig/SKILL.md` | **Rewrite** — new format reference, examples, field definitions |
| `CLAUDE.md` | **Update** — paperconfig format section, examples |
| `.claude/agents/agent-a-paperconfig.md` | **Update** — if it references paperconfig structure |

### 1j. `paperconfig_utils.py` — shared utility module

**File:** `multiomics_kg/utils/paperconfig_utils.py` (new)

Currently, paperconfig loading and traversal logic is duplicated across 8+ files:
- `load_all_paperconfigs()` lives in `build_gene_id_mapping.py`, imported by
  `resolve_paper_ids.py`; other consumers re-implement the same YAML loading loop
- `get_organism_for_entry()` lives in `build_gene_id_mapping.py`; other files
  duplicate the organism extraction logic
- Every consumer independently does `yaml.safe_load()`, iterates
  `supplementary_materials`, filters by `type == "csv"`, etc.

This module centralizes all common paperconfig operations. Created in two phases:
**Commit 0** (old-format functions) and **Commit 2** (new-format additions).

**Commit 0 — old-format functions (pure refactor):**

```python
"""Shared utilities for reading and traversing paperconfig.yaml files."""
import re
from pathlib import Path
import yaml

PROJECT_ROOT = Path(__file__).parent.parent.parent
PAPERCONFIG_FILES_TXT = (
    PROJECT_ROOT / "data/Prochlorococcus/papers_and_supp/paperconfig_files.txt"
)


# ─── Loading ──────────────────────────────────────────────────────────

def load_paperconfig(path: Path) -> dict:
    """Load a single paperconfig.yaml file."""
    with open(path) as f:
        return yaml.safe_load(f) or {}


def load_all_paperconfigs(
    list_file: Path = PAPERCONFIG_FILES_TXT,
) -> list[tuple[Path, dict]]:
    """Load all paperconfigs from paperconfig_files.txt.

    Returns list of (path, config_dict). Skips comments, blank lines,
    missing files (with warning).
    """
    results = []
    with open(list_file) as f:
        for line in f:
            path_str = line.strip()
            if not path_str or path_str.startswith('#'):
                continue
            path = PROJECT_ROOT / path_str
            if not path.exists():
                print(f"  [warn] paperconfig not found: {path}")
                continue
            results.append((path, load_paperconfig(path)))
    return results


# ─── Traversal helpers ───────────────────────────────────────────────

def get_publication(config: dict) -> dict:
    """Get the publication block."""
    return config.get("publication", {})


def get_paper_name(config: dict, fallback_path: Path | None = None) -> str:
    """Get paper name, with fallback to directory name."""
    name = get_publication(config).get("papername")
    if name:
        return name
    if fallback_path:
        return fallback_path.parent.name
    return "unknown"


def get_supplementary_materials(config: dict) -> dict:
    """Get supplementary_materials dict."""
    pub = get_publication(config)
    return pub.get("supplementary_materials", {})


def iter_csv_tables(config: dict):
    """Yield (table_key, table_config) for csv-type supplementary tables."""
    for key, table in get_supplementary_materials(config).items():
        if table.get("type", "csv") == "csv":
            yield key, table


def iter_analyses(config: dict):
    """Yield (table_key, table_config, analysis) for all analyses across all csv tables."""
    for table_key, table in iter_csv_tables(config):
        for analysis in table.get("statistical_analyses", []):
            yield table_key, table, analysis


def get_organism_for_entry(config: dict, entry: dict) -> str | None:
    """Get organism for a supplementary table entry.

    Old format: reads from entry-level 'organism' or first analysis's 'organism'.
    Falls back to entry-level 'organism' for id_translation/annotation_gff entries.
    """
    # Direct organism on entry (id_translation, annotation_gff, or table-level)
    org = entry.get("organism")
    if org:
        return str(org).strip().strip('"')
    # From first analysis's organism field (old format)
    for analysis in entry.get("statistical_analyses", []):
        org = analysis.get("organism")
        if org:
            return str(org).strip().strip('"')
    return None


# ─── Timepoint normalization ─────────────────────────────────────────

def parse_timepoint_hours(timepoint: str | None) -> float | None:
    """Parse a timepoint string to hours. Returns None for unparseable values.

    Handles: "4h", "0.5h", "-12h", "day 18", "Day 2", "50h (P added)",
    "0.5h post-inoculation", "1h extended darkness (36h)".
    Returns None for: "R (rescue: ...)", "days 60+89", None.
    """
    if not timepoint:
        return None
    tp = timepoint.strip()

    # "days N+M" — pooled, unparseable
    if tp.lower().startswith("days ") and "+" in tp:
        return None
    # "R (rescue: ...)" — non-numeric
    if tp.startswith("R ") or tp == "R":
        return None
    # "Nh extended darkness (Mh)" — use absolute time M
    if "extended darkness" in tp and "(" in tp:
        m = re.search(r'\((\d+(?:\.\d+)?)h\)', tp)
        return float(m.group(1)) if m else None
    # "day N" / "Day N"
    if tp.lower().startswith("day "):
        m = re.match(r'[Dd]ay\s+(\d+(?:\.\d+)?)', tp)
        return float(m.group(1)) * 24 if m else None
    # "Nh", "N h", "-Nh", "0.5h post-inoculation", "50h (P added)"
    m = re.match(r'(-?\d+(?:\.\d+)?)\s*h', tp)
    return float(m.group(1)) if m else None
```

**Commit 2 — new-format additions (after paperconfigs are migrated):**

```python
# ─── Experiment lookup (new format) ──────────────────────────────────

def get_experiments(config: dict) -> dict:
    """Get the experiments block from a paperconfig."""
    return config.get("publication", {}).get("experiments", {})


def get_experiment_for_analysis(config: dict, analysis: dict) -> dict:
    """Look up the experiment definition for a given analysis entry."""
    exp_key = analysis.get("experiment")
    if not exp_key:
        raise ValueError(
            f"Analysis '{analysis.get('id')}' missing 'experiment' reference"
        )
    experiments = get_experiments(config)
    if exp_key not in experiments:
        raise ValueError(
            f"Analysis '{analysis.get('id')}' references unknown "
            f"experiment '{exp_key}'"
        )
    return experiments[exp_key]


def get_organism_for_analysis(config: dict, analysis: dict) -> str:
    """Get organism for an analysis (from its experiment block)."""
    return get_experiment_for_analysis(config, analysis)["organism"]
```

**Commit 2 also updates `get_organism_for_entry()`** to use the new-format path:

```python
def get_organism_for_entry(config: dict, entry: dict) -> str | None:
    """Get organism for a supplementary table entry.

    New format: looks up from first analysis's experiment block.
    Falls back to entry-level 'organism' for id_translation/annotation_gff entries.
    """
    # Direct organism on entry (id_translation, annotation_gff)
    org = entry.get("organism")
    if org:
        return str(org).strip().strip('"')
    # From first analysis's experiment block (new format)
    for analysis in entry.get("statistical_analyses", []):
        try:
            return get_organism_for_analysis(config, analysis)
        except (ValueError, KeyError):
            pass
    return None
```

**Used by:** `omics_adapter.py`, `build_gene_id_mapping.py`, `resolve_paper_ids.py`,
`check_gene_ids.py`, `fix_gene_ids.py`, `build_gene_mapping_supp.py`,
`validate_paperconfig.py`, `migrate_paperconfigs.py`, and test files.

**Migration note:** `build_gene_id_mapping.py` currently owns `load_all_paperconfigs()`
and `get_organism_for_entry()`. These move to `paperconfig_utils.py` in commit 0;
delete the old functions, don't leave wrappers.

---

## Step 2: Schema changes (`config/schema_config.yaml`)

### 2a. Add Experiment node

```yaml
experiment:
  is_a: InvestigativeProcess
  represented_as: node
  preferred_id: id
  label_in_input: experiment
  properties:
    name: str                    # human-readable description
    organism_strain: str         # "Prochlorococcus MED4"
    treatment_type: str          # "nitrogen_stress", "coculture", "light_stress", etc.
    treatment: str               # time-stripped treatment condition
    control: str                 # control condition
    experimental_context: str    # other factors held constant
    coculture_partner: str       # partner organism name (null for environmental)
    omics_type: str              # RNASEQ | PROTEOMICS | METABOLOMICS | MICROARRAY
    statistical_test: str        # DESeq2, edgeR, etc.
    is_time_course: str          # "true" | "false" (Neo4j import uses strings)
    medium: str                  # growth medium (e.g., "Pro99", "PRO99-lowN")
    temperature: str             # growth temperature (e.g., "24C")
    light_condition: str         # light regime (e.g., "continuous_light", "13:11 light:dark cycle")
    light_intensity: str         # PAR intensity (e.g., "30 umol photons m-2 s-1")
```

### 2b. Add new edge types

```yaml
# Publication → Experiment
publication has experiment:
  is_a: Association
  represented_as: edge
  label_as_edge: has_experiment
  label_in_input: has_experiment
  source: publication
  target: experiment

# Experiment → OrganismTaxon (coculture partner)
experiment tests coculture:
  is_a: Association
  represented_as: edge
  label_as_edge: tests_coculture_with
  label_in_input: tests_coculture_with
  source: experiment
  target: organism taxon

# Experiment → Gene (expression results)
experiment to gene expression association:
  is_a: Association
  represented_as: edge
  label_as_edge: changes_expression_of
  label_in_input: changes_expression_of
  source: experiment
  target: gene
  properties:
    time_point: str                # original label (null for single-point)
    time_point_order: int          # ordinal position (1-indexed)
    time_point_hours: float        # normalized to hours (nullable)
    log2_fold_change: float
    adjusted_p_value: float
    expression_direction: str      # "up" | "down"
    significant: str               # "significant" | "not significant" | "unknown"
```

Note: `label_as_edge` and `label_in_input` are both lowercase with underscores, following
the dominant pattern in the existing schema (e.g., `coculture_changes_expression_of`,
`gene_in_ortholog_group`). BioCypher v0.12.5 converts these to `Capital_snake_case` for
Neo4j labels (e.g., `changes_expression_of` → `Changes_expression_of`,
`has_experiment` → `Has_experiment`). Verified via toy KG test — see resolved item #1.

This unifies the two old edge labels (`Condition_changes_expression_of`,
`Coculture_changes_expression_of`) into a single `Changes_expression_of` — the source
node is always `Experiment` now, so the condition/coculture prefix is no longer needed.
Reads naturally: "Experiment changes expression of Gene."

### 2c. Remove old types

Delete from schema:
- `coculture to gene expression association` (coculture_changes_expression_of)
- `condition to gene expression association` (condition_changes_expression_of)
- `publication to expression source association` (published_expression_data_about)
- `environmental condition` node type (absorbed into Experiment)

### 2d. Remove unused timeseries_cluster types

Delete from schema (will be replaced by Phase 2 ExpressionCluster design):
- `timeseries_cluster`
- `cluster_in_publication`
- `molecular_in_cluster`

### 2e. Acceptance criteria

- [ ] Schema has `experiment` node with all properties (including `medium`, `temperature`)
- [ ] Schema has `has_experiment`, `tests_coculture_with`, `changes_expression_of` edges
- [ ] Old expression edge types removed
- [ ] `environmental condition` node type removed
- [ ] Old timeseries_cluster types removed
- [ ] BioCypher can parse the schema without errors (`uv run python create_knowledge_graph.py --test`)

---

## Step 3: Adapter changes (`multiomics_kg/adapters/omics_adapter.py`)

This is the largest change. The adapter needs to:
1. Read experiment definitions from the new `experiments` block in paperconfig
2. Emit Experiment nodes with properties directly from the block
3. Emit structural edges (HAS_EXPERIMENT, TESTS_COCULTURE_WITH)
4. Emit `changes_expression_of` edges (Experiment → Gene) instead of the old source → Gene edges
5. Stop emitting EnvironmentalCondition nodes, Published_expression_data_about edges

### 3a. Reading experiments from paperconfig

In `get_nodes()`, after loading the paperconfig:

1. Read the `experiments` block — each key is an experiment ID, value is experiment metadata
2. For each experiment, scan all `statistical_analyses` entries that reference it
   (via `analysis.experiment == experiment_key`)
3. Collect timepoints from referencing analyses → sorted by `timepoint_hours`
4. Determine `is_time_course` = count of analyses with timepoints > 1
5. Build Experiment node properties directly from the experiment block fields

No grouping logic needed — the paperconfig defines the grouping explicitly.
If an analysis is missing the `experiment` reference, **log a warning and skip it**.

### 3b. Experiment node ID

Format: `{doi}_{experiment_group_id}`

Example: `doi:10.1038/ismej.2016.70_coculture_hot1a3_med4_rnaseq`

The DOI prefix ensures global uniqueness. The `experiment_group_id` slug provides
human readability and stability.

### 3c. Experiment node emission (get_nodes)

For each experiment in the `experiments` block:

```python
# experiment_key = "n_starvation_axenic_hot1a3_rnaseq"
# exp = the experiment dict from paperconfig
yield (
    experiment_id,        # "doi:10.xxx_n_starvation_axenic_hot1a3_rnaseq"
    "experiment",         # label_in_input
    {
        "name": exp["name"],
        "organism_strain": exp["organism"],
        "treatment_type": exp["treatment_type"],
        "treatment": exp["treatment_condition"],
        "control": exp["control_condition"],
        "experimental_context": exp.get("experimental_context"),
        "coculture_partner": exp.get("treatment_organism"),  # None for environmental
        "omics_type": exp["omics_type"],
        "statistical_test": exp["test_type"],
        "is_time_course": "true" if num_timepoints > 1 else "false",
        "medium": exp.get("medium"),
        "temperature": exp.get("temperature"),
        "light_condition": exp.get("light_condition"),
        "light_intensity": exp.get("light_intensity"),
    }
)
```

All properties come directly from the paperconfig experiments block — no lookup tables,
no derivation logic. The migration script already placed everything in the right spot.

**NEEDS_CURATION handling:** Properties set to `"NEEDS_CURATION"` in paperconfigs are
emitted as-is to the graph. The adapter does NOT convert them to null — keeping the
sentinel makes incomplete data grepable in both paperconfigs and the graph. Post-migration
curation pass will fill these in from paper PDFs. Find all with:
`grep -r "NEEDS_CURATION" data/` (paperconfigs) or
`MATCH (e:Experiment) WHERE e.light_condition = 'NEEDS_CURATION' RETURN e` (graph).

### 3d. Structural edge emission (get_edges)

For each experiment group, emit:

1. **HAS_EXPERIMENT**: `(Publication) → (Experiment)`
   - Edge ID: `{pub_id}_has_exp_{experiment_group_id}`

2. **TESTS_COCULTURE_WITH** (coculture experiments only):
   - Source: Experiment node
   - Target: OrganismTaxon node (treatment organism)
   - Edge ID: `{experiment_id}_tests_coculture_{organism_id}`

### 3e. Changes_expression_of edge emission (get_edges)

Replace `_load_and_create_edges()` to emit `changes_expression_of` edges:

For each analysis within an experiment group, for each row in the data file:

```python
yield (
    None,               # relationship_id (auto or manual)
    experiment_id,      # source: Experiment node
    gene_id,            # target: Gene node (ncbigene:locus_tag)
    "changes_expression_of",  # label_in_input
    {
        "time_point": timepoint or None,
        "time_point_order": timepoint_order,   # 1-indexed position in series
        "time_point_hours": timepoint_hours,   # float or None
        "log2_fold_change": log2fc,
        "adjusted_p_value": padj,
        "expression_direction": direction,
        "significant": significance,
    }
)
```

**Edge ID:** `{pub_id}_{analysis_id}_{gene_id}` (same pattern as current — unique
per analysis per gene).

**time_point_order:** Computed from the sorted list of timepoints within the experiment
group. For single-point experiments and analyses without a `timepoint` field,
`time_point_order = 1` (not null).

### 3f. Removed from adapter

- EnvironmentalCondition node emission (absorbed into Experiment)
- `environmental_conditions` block parsing (no longer needed — experiments block has everything)
- `Published_expression_data_about` edge emission (the `seen_sources` tracking loop)
- `_determine_source_id()` method (no longer needed — source is always Experiment)
- Edge properties: `publications`, `omics_type`, `organism_strain`, `treatment_condition`,
  `control_condition`, `experimental_context`, `statistical_test`, `analysis_name`
  (all on Experiment node now)

### 3g. Retained in adapter

- Publication node emission (unchanged)
- PDF metadata extraction (unchanged)
- Pre-resolved CSV probing (unchanged)
- Significance determination logic (unchanged)
- Gene ID resolution logic (unchanged)
- String sanitization (unchanged)

### 3h. Acceptance criteria

- [ ] Experiment nodes are created with correct properties from experiments block
- [ ] Analyses missing `experiment` reference are logged and skipped
- [ ] HAS_EXPERIMENT edges connect publications to experiments
- [ ] TESTS_COCULTURE_WITH edges connect experiments to organism taxa
- [ ] No EnvironmentalCondition nodes are emitted
- [ ] `changes_expression_of` edges have correct time_point, time_point_order, time_point_hours
- [ ] No old edge types are emitted
- [ ] Total `changes_expression_of` edge count matches current total expression edge count (~188K)
- [ ] `--test` mode works correctly

---

## Step 4: Post-import updates (`scripts/post-import.sh`, `scripts/post-import.cypher`)

### 4a. New indexes

```cypher
// Experiment node indexes
CREATE INDEX experiment_id_idx FOR (e:Experiment) ON (e.id);
CREATE INDEX experiment_organism_idx FOR (e:Experiment) ON (e.organism_strain);
CREATE INDEX experiment_treatment_type_idx FOR (e:Experiment) ON (e.treatment_type);
CREATE INDEX experiment_omics_type_idx FOR (e:Experiment) ON (e.omics_type);

// Full-text index on Experiment for keyword search
CREATE FULLTEXT INDEX experimentFullText FOR (e:Experiment)
ON EACH [e.name, e.treatment, e.control, e.experimental_context, e.light_condition];
```

### 4b. Update Gene routing signals

Current Cypher references `Condition_changes_expression_of` and
`Coculture_changes_expression_of`. Update to use `Changes_expression_of`:

```cypher
// expression_edge_count — total Changes_expression_of edges per gene
MATCH (g:Gene)
OPTIONAL MATCH (e:Experiment)-[r:Changes_expression_of]->(g)
WITH g, count(r) AS expr_count
SET g.expression_edge_count = expr_count;

// significant_expression_count
MATCH (g:Gene)
OPTIONAL MATCH (e:Experiment)-[r:Changes_expression_of {significant: 'significant'}]->(g)
WITH g, count(r) AS sig_count
SET g.significant_expression_count = sig_count;
```

### 4c. Compute `rank_by_effect` on Changes_expression_of edges

Rank each gene by absolute fold-change within each experiment + timepoint. Computed
post-import because the adapter processes rows one at a time and doesn't know all genes
in an experiment. This gives the LLM instant context: "PMM0120 is the 15th most
changed gene (out of 1696) in this experiment at this timepoint."

```cypher
// Rank genes by |log2FC| within each experiment + timepoint (1 = largest effect)
// Batched per-experiment to bound memory (~1K–3K edges per group)
MATCH (e:Experiment)
WITH e
CALL {
  WITH e
  MATCH (e)-[r:Changes_expression_of]->(g:Gene)
  WITH r.time_point_order AS tp, r, abs(r.log2_fold_change) AS abs_fc
  ORDER BY tp, abs_fc DESC
  WITH tp, collect(r) AS edges
  UNWIND range(0, size(edges)-1) AS i
  SET (edges[i]).rank_by_effect = i + 1
} IN TRANSACTIONS OF 10 ROWS;
```

**Notes:**
- Rank 1 = largest absolute fold-change in that experiment+timepoint
- Ties broken by Neo4j's internal ordering (stable but arbitrary) — acceptable
- Not ranking by padj because it's often null (prefiltered datasets)
- `IN TRANSACTIONS OF 10 ROWS` processes 10 experiments per transaction, bounding
  memory to ~10K–30K edges at a time. Each experiment's `collect()` holds ~1K–3K edges.

### 4d. Remove old references

Remove any Cypher that references:
- `Condition_changes_expression_of`
- `Coculture_changes_expression_of`
- `Published_expression_data_about`

### 4e. Acceptance criteria

- [ ] New indexes created for Experiment nodes
- [ ] Gene routing signals updated to use `Changes_expression_of`
- [ ] `rank_by_effect` computed on all `Changes_expression_of` edges (spot-check: rank 1 has largest |FC|)
- [ ] No references to old edge types in post-import scripts
- [ ] Both `post-import.sh` and `post-import.cypher` are in sync

---

## Step 5: Test updates

### 5a. Unit tests

**Adapter tests (`tests/test_omics_adapter*.py`):**
- Update existing expression edge tests to expect `changes_expression_of` label
- Add tests for Experiment node creation
- Add tests for experiment grouping logic (single-point, time-course, auto-fallback)
- Add tests for structural edges (HAS_EXPERIMENT, TESTS_COCULTURE_WITH)
- Add tests for time_point_order and time_point_hours computation
- Verify edge properties match new minimal schema (no duplicated metadata)

**paperconfig_utils tests (`tests/test_paperconfig_utils.py`) — new file, created in Commit 0:**
- `parse_timepoint_hours()`: all patterns from 1e table (hours, days, negative, extended
  darkness, rescue, pooled days, post-inoculation, null input)
- `iter_analyses()` / `iter_csv_tables()`: correct traversal of multi-table configs
- `get_experiment_for_analysis()`: valid reference, missing reference, unknown key
- `get_organism_for_analysis()`: organism from experiment block
- `get_organism_for_entry()`: fallback to entry-level organism for id_translation entries
- `load_all_paperconfigs()`: skips comments, blank lines, missing files

**Paperconfig validation tests (`tests/test_paperconfig_validation.py`):**
- Update all fixtures from old format (with `environmental_conditions`) to new format
  (with `experiments` block)
- Update canonical vocabulary tests for `treatment_type` values
- Add tests for `experiment` reference validation on analyses
- Add tests for `timepoint_hours` validation

### 5b. KG validity tests (`tests/kg_validity/`)

- `test_expression.py`: Update to query `Changes_expression_of` edges instead of old types
- `test_structure.py`: Add Experiment node presence checks, minimum counts
- `test_post_import.py`: Verify new indexes, update routing signal checks
- `test_snapshot.py`: Regenerate snapshot after rebuild
- Remove EnvironmentalCondition node checks (nodes no longer exist)
- New: `test_experiment.py` — Experiment node properties (including medium, temperature),
  HAS_EXPERIMENT connectivity, time-course grouping, treatment_type values,
  coculture_partner presence, TESTS_COCULTURE_WITH edges

### 5c. Acceptance criteria

- [ ] `pytest -m "not slow and not kg"` passes
- [ ] `pytest -m kg` passes against rebuilt graph
- [ ] Snapshot regenerated

---

## Step 6: Documentation and cleanup

### 6a. CLAUDE.md updates

- Update "Key graph facts" section: Experiment nodes, DE_IN edges, removed edge types,
  EnvironmentalCondition nodes removed
- Update "Actual Neo4j labels" lists (remove EnvironmentalCondition, add Experiment)
- Update post-import description
- Update edge count estimates

### 6b. Validate paperconfig skill

Update `.claude/skills/paperconfig/` to include `experiment_group_id`,
`timepoint_hours` in the template and validation.

### 6c. Omics edge snapshot skill

Update `.claude/skills/omics-edge-snapshot/` to count `Changes_expression_of` instead of old edge types.

### 6d. Cypher queries skill

Update `.claude/skills/cypher-queries/` templates to use new edge types.

---

## Implementation order

```
── PRE-WORK: baseline ─────────────────────────────────────────────────
│
│  git tag pre-experiment-redesign         ← rollback point
│  Run /omics-edge-snapshot                ← "before" baseline while old edge types exist
│
── COMMIT 0 (mandatory): paperconfig_utils.py ─────────────────────────
│
│  Step 0a: Create paperconfig_utils.py    ← loading, traversal, timepoint parsing
│  Step 0b: Write tests for paperconfig_utils ← test_paperconfig_utils.py
│  Step 0c: Migrate ALL paperconfig consumers to use paperconfig_utils:
│           │  - build_gene_id_mapping.py   ← owns load_all_paperconfigs(), delete original
│           │  - resolve_paper_ids.py       ← imports from build_gene_id_mapping, switch
│           │  - omics_adapter.py           ← YAML loading, supp_materials iteration
│           │  - validate_paperconfig.py    ← YAML loading, field traversal
│           │  - check_gene_ids.py          ← organism/name_col extraction
│           │  - fix_gene_ids.py            ← organism/name_col/id_columns extraction
│           │  - build_gene_mapping_supp.py ← organism/name_col extraction
│           │  - validate_annotations.py    ← filename/name_col extraction (if applicable)
│  Step 0d: Run pytest -m "not slow and not kg" ← verify nothing breaks
│  Step 0e: Run prepare_data.sh            ← verify pipeline still works
│  Step 0f: CODE REVIEW                    ← review diff: no old loading/traversal code
│           │                                 left, all consumers use paperconfig_utils
│  Step 0g: MANUAL COMMIT
│
│  This commit is pure refactor — no format change yet. All paperconfig
│  consumers switch to shared utils. Must land first to reduce the size
│  and risk of commit 1. Delete old duplicated loading/traversal code —
│  don't leave wrappers.
│
── COMMIT 1: dry-run migration + validate grouping ────────────────────
│
│  Step 1a: Write migration script         ← scripts/migrate_paperconfigs.py
│           │                                 (uses paperconfig_utils for parsing)
│           │                                 Supports --dry-run: writes to temp dir,
│           │                                 does NOT overwrite originals
│  Step 1b: Run --dry-run on all 26 papers ← outputs new-format paperconfigs to temp dir
│  Step 1c: Manual inspection              ← review experiment grouping for each paper:
│           │  - Are groups correct? (same biological question)
│           │  - Are experiment names readable?
│           │  - treatment_condition time-stripping correct?
│           │  - timepoint_hours values correct?
│           │  - treatment_type assigned for every experiment?
│           │  - Coculture experiments have treatment_organism?
│  Step 1d: Script flags for review        ← migration script reports:
│           │  - Near-duplicate groups (similar but not identical grouping keys)
│           │  - Analyses with null timepoint_hours (unparseable)
│           │  - Experiments with only 1 timepoint (expected single-point?)
│           │  - Missing treatment_type or treatment_organism
│           │  - Round-trip validation: same # analyses before and after
│  Step 1e: Resolve any flagged issues     ← fix paperconfig text inconsistencies,
│           │                                 adjust grouping rules if needed
│  Step 1f: CODE REVIEW                    ← review migration script logic: grouping,
│           │                                 time-stripping, slug generation, flag reports
│  Step 1g: MANUAL COMMIT
│
│  This commit lands the migration script only. No paperconfig files
│  changed yet. Gives confidence that the conversion is sound before
│  touching any consumers.
│
── ATOMIC COMMIT 2: apply migration + update consumers ────────────────
│
│  Step 2a: Run migration for real         ← overwrites all 26 paperconfigs
│  Step 2b: Add new-format functions to    ← get_experiment_for_analysis(),
│           paperconfig_utils.py             get_organism_for_analysis()
│  Step 2c: Update all consumers to use    ← all already on paperconfig_utils (commit 0),
│           new-format helpers                now switch organism lookup to experiment block
│           │  - build_gene_id_mapping.py
│           │  - resolve_paper_ids.py
│           │  - check_gene_ids.py, fix_gene_ids.py, build_gene_mapping_supp.py
│  Step 2d: Update validate_paperconfig.py ← new format validation (experiments block,
│           │                                 experiment references, timepoint_hours)
│  Step 2e: Update test fixtures           ← test_paperconfig_validation.py,
│           │                                 test_omics_adapter_organism_gene.py
│  Step 2f: Update /paperconfig skill      ← new format template and examples
│  Step 2g: Run pytest -m "not slow and not kg" ← adapter tests may fail (expected —
│           │                                      adapter still emits old edge types,
│           │                                      rewritten in commit 3)
│  Step 2h: Run prepare_data.sh --steps 3 4 ← verify pipeline still works
│  Step 2i: CODE REVIEW                    ← review: migrated paperconfigs correct,
│           │                                 all consumers use new-format helpers,
│           │                                 no old env_conditions references remain
│  Step 2j: MANUAL COMMIT
│
│  ALL of 2a–2j committed together. Pipeline must not break.
│
── COMMIT 3: adapter + schema + post-import ───────────────────────────
│
│  Step 3: Schema changes                  ← small, mechanical
│  Step 4: Adapter rewrite                 ← reads experiments block (depends on 2+3)
│  Step 5: Post-import updates             ← indexes, routing signals, rank_by_effect
│  Step 5b: Update /omics-edge-snapshot    ← count Changes_expression_of instead of old types
│  Step 6: Unit tests for adapter          ← new tests for Experiment nodes/edges
│  Step 6b: Run pytest -m "not slow and not kg" ← verify unit tests pass
│  Step 6c: CODE REVIEW                    ← review: schema + adapter + post-import
│           │                                 consistent, no old edge types remain,
│           │                                 snapshot skill updated
│  Step 6d: MANUAL COMMIT
│
── REBUILD + VALIDATE (manual) ────────────────────────────────────────
│
│  Step 7a: Full rebuild                   ← manual: docker compose up
│  Step 7b: Run pytest -m kg               ← KG validity tests against rebuilt graph
│  Step 7c: Run /omics-edge-snapshot       ← compare against pre-work baseline
│  Step 7d: Manual spot-checks             ← verify Experiment nodes, edge counts, properties
│
── COMMIT 4: docs and cleanup ─────────────────────────────────────────
│
│  Step 8: Docs                            ← CLAUDE.md, SKILL.md, agent docs
│  Step 8b: CODE REVIEW                    ← review: docs match new schema/edges,
│           │                                 no stale references
│  Step 8c: MANUAL COMMIT
```

**Critical constraint:** Commit 0 must land first — all consumers on shared utils
before any format change. Commit 1 is dry-run only — validates grouping without
touching paperconfigs or consumers. Commit 2 is atomic — the migrated paperconfigs
and ALL their consumers (now using new-format helpers) must ship together.

**Critical path:** paperconfig_utils (0) → dry-run migration + review (1) →
apply migration + consumer updates (2) → adapter rewrite (4) → rebuild (7)

**No parallel track** — schema changes are small and included in commit 3.

---

## Implementation concerns (from review)

1. **Migration script grouping key is fragile.** Grouping by exact string match on
   `{organism, treatment_condition (time-stripped), control_condition, experimental_context,
   type, test_type}` will split on minor text differences (e.g., `"24C"` vs `"24°C"`).
   Add a normalization pass or emit warnings for near-duplicate groups.

2. **`paperconfig_utils.py` must replace, not wrap.** When moving `load_all_paperconfigs()`
   and `get_organism_for_entry()` from `build_gene_id_mapping.py`, delete the old functions
   entirely — don't leave thin wrappers that rot.

3. **~~`check-gene-ids` and `fix-gene-ids` skills must be in Commit 2, not parallel.~~**
   — **Added to plan.** Skills updated in Commit 2, Phase 2.2 (Agent D).

4. **~~`/omics-edge-snapshot` must be updated (or "before" snapshot taken) before rebuild.~~**
   — **Added to plan.** "Before" snapshot taken in PRE-WORK; skill updated in Commit 3,
   Phase 3.2 (Agent D); "after" snapshot taken post-rebuild.

5. **~~`significant` property bug~~** — **RESOLVED.** Root cause was a BioCypher bug with
   boolean properties. Already fixed; no action needed in this redesign.

6. **~~MCP tool impact~~** — Handled separately in the MCP repo
   (`multiomics_explorer/plans/redefine_mcp_tools/expression_tools_redesign.md`).

7. **~~Rollback strategy~~** — Git tag `pre-experiment-redesign` before rebuild is sufficient.
   Added to plan.

---

## Risk assessment

| Risk | Mitigation |
|---|---|
| Edge count mismatch after rebuild | Snapshot before/after with `/omics-edge-snapshot` |
| Migration script produces wrong grouping | Manual review of each paper's experiments block after migration; round-trip validation in script |
| Weissberg treatment_condition cleanup loses biological meaning | Keep phase-descriptive treatment text ("Nutrient starvation", "Long-term nutrient starvation") — only strip parenthesized time |
| MCP tools broken during transition | Coordinated cutover: both repos deploy together |
| BioCypher label case mismatch | **Resolved.** `changes_expression_of` → `Changes_expression_of` in Neo4j. Matches existing pattern. |
| `id_translation` / `annotation_gff` entries break | Migration script must preserve non-csv supplementary_materials entries unchanged |

---

## Open items (to resolve during implementation)

1. **~~BioCypher label for expression edges~~** — **RESOLVED.** Edge label is
   `changes_expression_of` (input) → `Changes_expression_of` (Neo4j). Unifies
   `Condition_changes_expression_of` + `Coculture_changes_expression_of` — the
   condition/coculture prefix is no longer needed since the source is always Experiment.
   Verified via toy KG that BioCypher v0.12.5 uses `Capital_snake_case`.

2. **~~Lin 2015 grouping~~** — **RESOLVED.** Two separate experiment groups:
   `p_starvation_med4_microarray` and `p_readd_med4_microarray`. Different
   biological questions (starvation vs rescue).

3. **~~Biller 2016 late comparisons~~** — **RESOLVED.** Separate experiment group.
   Different control (12h coculture instead of -12h axenic) = different experiment.

4. **~~Weissberg "days 60+89" pooled samples~~** — **RESOLVED.** Include in the
   time-course experiment with `time_point_order = max+1` and `timepoint_hours: null`.
   Keeps the relationship to the same biological trajectory.

5. **`name` field on experiments:** Migration script auto-generates using the pattern
   `"{organism_short} {treatment} vs {control} ({omics_type})"`. Manual cleanup after
   migration if any names are awkward.

---

## Agent assignments

### File ownership

| Agent | Owned files |
|-------|-------------|
| **B** (Schema + Adapters) | `multiomics_kg/utils/paperconfig_utils.py` (new), `config/schema_config.yaml`, `multiomics_kg/adapters/omics_adapter.py`, `multiomics_kg/download/build_gene_id_mapping.py`, `multiomics_kg/download/resolve_paper_ids.py`, `scripts/post-import.sh`, `scripts/post-import.cypher`, `scripts/migrate_paperconfigs.py` (new), `scripts/validate_annotations.py` |
| **A** (Paperconfig) | All 26 `paperconfig.yaml` files |
| **C** (Tests) | `tests/test_paperconfig_utils.py` (new), `tests/test_paperconfig_validation.py`, `tests/test_omics_adapter_organism_gene.py`, `tests/kg_validity/test_expression.py`, `tests/kg_validity/test_experiment.py` (new), `tests/kg_validity/test_post_import.py`, `tests/kg_validity/test_snapshot.py`, `tests/kg_validity/snapshot_data.json` |
| **D** (Skills) | `.claude/skills/paperconfig/validate_paperconfig.py`, `.claude/skills/paperconfig/SKILL.md`, `.claude/skills/omics-edge-snapshot/SKILL.md`, `.claude/skills/cypher-queries/SKILL.md`, `.claude/skills/check-gene-ids/check_gene_ids.py`, `.claude/skills/fix-gene-ids/fix_gene_ids.py`, `.claude/skills/build-gene-mapping-supp/build_gene_mapping_supp.py` |
| **E** (Docs) | `CLAUDE.md`, memory files |
| **F** (Validation) | Read-only — runs tests, skills, snapshots |
| **G** (Code review) | Read-only — reviews diffs at each gate |
| **H** (Plan manager) | `plans/experiment_node_redesign.md` |

### Task list

```
═══ PRE-WORK ═══════════════════════════════════════════════════════════

  H: Update .claude/agents/ definitions for experiment redesign
  F: Run /omics-edge-snapshot (save "before" baseline)
  Manual: git tag pre-experiment-redesign

═══ COMMIT 0: paperconfig_utils.py ═════════════════════════════════════

  ── Phase 0.1 (sequential) ──────────────────────────────────────────
  │
  │  B: Create multiomics_kg/utils/paperconfig_utils.py
  │     (loading, traversal, timepoint parsing — old-format only)
  │
  ── Phase 0.2 (parallel — no file conflicts) ────────────────────────
  │
  │  B: Migrate own files to use paperconfig_utils:
  │     - build_gene_id_mapping.py  (delete load_all_paperconfigs, get_organism_for_entry)
  │     - resolve_paper_ids.py      (switch imports)
  │     - omics_adapter.py          (YAML loading, supp_materials iteration)
  │     - validate_annotations.py   (if applicable)
  │
  │  D: Migrate skill files to use paperconfig_utils:       ║ parallel
  │     - validate_paperconfig.py                           ║ with B
  │     - check_gene_ids.py                                 ║
  │     - fix_gene_ids.py                                   ║
  │     - build_gene_mapping_supp.py                        ║
  │
  │  C: Write tests/test_paperconfig_utils.py               ║ parallel
  │     (parse_timepoint_hours, iter_analyses,               ║ with B, D
  │      get_organism_for_entry, load_all_paperconfigs)      ║
  │
  ── Phase 0.3 (sequential — after 0.2) ──────────────────────────────
  │
  │  F: Run pytest -m "not slow and not kg"
  │  F: Run prepare_data.sh
  │  G: Code review — no old loading/traversal code left
  │  Manual: COMMIT

═══ COMMIT 1: dry-run migration ════════════════════════════════════════

  ── Phase 1.1 (sequential) ──────────────────────────────────────────
  │
  │  B: Write scripts/migrate_paperconfigs.py
  │     (grouping, time-stripping, slug generation, --dry-run mode,
  │      flag reports for near-duplicates / missing fields)
  │
  ── Phase 1.2 (sequential) ──────────────────────────────────────────
  │
  │  B: Run --dry-run on all 26 papers
  │  Manual: Inspect output — grouping, names, timepoint_hours, flags
  │  Manual: Resolve any flagged issues
  │  G: Code review — migration script logic
  │  Manual: COMMIT

═══ COMMIT 2: apply migration + update consumers ══════════════════════

  ── Phase 2.1 (sequential) ──────────────────────────────────────────
  │
  │  A: Run migration for real (overwrites all 26 paperconfigs)
  │
  │  B: Add new-format functions to paperconfig_utils.py
  │     (get_experiment_for_analysis, get_organism_for_analysis)
  │
  ── Phase 2.2 (parallel — after 2.1, no file conflicts) ─────────────
  │
  │  B: Update own files for new-format helpers:
  │     - build_gene_id_mapping.py  (organism via experiment block)
  │     - resolve_paper_ids.py      (organism via experiment block)
  │
  │  D: Update skill files for new format:                  ║ parallel
  │     - validate_paperconfig.py   (experiments block,     ║ with B
  │       experiment refs, timepoint_hours validation)      ║
  │     - check_gene_ids.py         (organism lookup)       ║
  │     - fix_gene_ids.py           (organism lookup)       ║
  │     - build_gene_mapping_supp.py (organism lookup)      ║
  │     - SKILL.md (paperconfig)    (new format template)   ║
  │
  │  C: Update test fixtures for new format:                ║ parallel
  │     - test_paperconfig_validation.py                    ║ with B, D
  │     - test_omics_adapter_organism_gene.py               ║
  │
  ── Phase 2.3 (sequential — after 2.2) ──────────────────────────────
  │
  │  F: Run pytest -m "not slow and not kg"
  │     (adapter tests may fail — expected, rewritten in commit 3)
  │  F: Run prepare_data.sh --steps 3 4
  │  G: Code review — migrated paperconfigs correct, no old
  │     env_conditions references remain, all consumers updated
  │  Manual: COMMIT

═══ COMMIT 3: adapter + schema + post-import ═══════════════════════════

  ── Phase 3.1 (sequential) ──────────────────────────────────────────
  │
  │  B: All core changes:
  │     - config/schema_config.yaml  (Experiment node, new edges, remove old)
  │     - omics_adapter.py           (major rewrite: Experiment nodes,
  │       Changes_expression_of edges, Has_experiment, Tests_coculture_with)
  │     - post-import.sh + post-import.cypher  (new indexes,
  │       routing signals, rank_by_effect)
  │
  ── Phase 3.2 (parallel — after 3.1, no file conflicts) ─────────────
  │
  │  D: Update skill files:                                 ║
  │     - omics-edge-snapshot SKILL.md                      ║ parallel
  │       (count Changes_expression_of)                     ║ with C
  │     - cypher-queries SKILL.md                           ║
  │       (new edge type templates)                         ║
  │
  │  C: Write/update tests:                                 ║ parallel
  │     - test_omics_adapter_organism_gene.py               ║ with D
  │       (Experiment nodes, new edges, time_point_order)   ║
  │     - tests/kg_validity/test_experiment.py (new)        ║
  │     - tests/kg_validity/test_expression.py              ║
  │     - tests/kg_validity/test_post_import.py             ║
  │
  ── Phase 3.3 (sequential — after 3.2) ──────────────────────────────
  │
  │  F: Run pytest -m "not slow and not kg"
  │  G: Code review — schema + adapter + post-import consistent,
  │     no old edge types remain, snapshot skill updated
  │  Manual: COMMIT

═══ REBUILD + VALIDATE (manual) ════════════════════════════════════════

  Manual: docker compose up (full rebuild)
  F: Run pytest -m kg
  F: Run /omics-edge-snapshot (compare against pre-work baseline)
  Manual: Spot-checks — Experiment nodes, edge counts, properties

═══ COMMIT 4: docs and cleanup ═════════════════════════════════════════

  ── Phase 4.1 (parallel — no file conflicts) ─────────────────────────
  │
  │  E: Update CLAUDE.md (new schema, edge types, counts, labels)
  │     Update memory files
  │
  │  C: Regenerate tests/kg_validity/snapshot_data.json     ║ parallel
  │     Update test_snapshot.py if needed                   ║ with E
  │
  ── Phase 4.2 (sequential — after 4.1) ──────────────────────────────
  │
  │  G: Code review — docs match new schema, no stale references
  │  Manual: COMMIT
```

### Parallel execution rules

- Within a phase, agents marked "parallel" can run simultaneously
- No agent edits a file owned by another agent
- Agents G and F never edit files — they can run after any implementing agent
- Agent H is the only agent that edits plan files
- If two agents need to **read** the same file, that is fine
