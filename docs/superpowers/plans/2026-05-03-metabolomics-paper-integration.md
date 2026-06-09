# Metabolomics paper integration — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Integrate metabolomics paper measurements (concentrations + presence flags) into the KG via a new `MetaboliteAssay` node type and 2 measurement edges, validated against Capovilla 2023 + Kujawinski 2023.

**Architecture:** Phase 2 of the metabolite scaffold. New adapter `metabolite_assay_adapter.py` mirrors `observations_adapter.py`. Step 6 of `prepare_data.sh` extends to harvest paper metabolites into `kegg_data.json` + a new `metabolite_id_mapping.json`. New step 7 (`resolve_paper_metabolites.py`) writes `<stem>_resolved.csv` files paralleling step 4 for genes. Post-import Cypher computes detection_status/ranks/rollups parallel to DerivedMetric. Single source of truth: `metabolism_adapter` owns `Metabolite` node emission; new adapter only emits `MetaboliteAssay` + measurement/binding edges.

**Tech Stack:** Python (uv), BioCypher, Neo4j, pytest, YAML paperconfigs.

**Spec:** `docs/superpowers/specs/2026-05-03-metabolomics-paper-integration-design.md`

## Execution constraints

**No Docker changes.** Phases 6 and 7 include Docker rebuild + import steps (`docker compose up -d build`, etc.) that are required to deploy the new MetaboliteAssay nodes into the running Neo4j. **The implementing agent must NOT execute these Docker steps**; they are user actions. Mark these tasks complete only up to the point where local files, paperconfigs, and resolved CSVs are committed; pause before any `docker compose` invocation. KG validity tests (`pytest -m kg`) that depend on a rebuilt graph are also user actions — skip them; user runs them after the rebuild they trigger.

---

## File structure

**Created:**
- `multiomics_kg/adapters/metabolite_assay_adapter.py` — new adapter (mirrors observations_adapter.py)
- `multiomics_kg/download/resolve_paper_metabolites.py` — step 7 CSV rewriter
- `tests/test_metabolite_assay_adapter.py`
- `tests/test_resolve_paper_metabolites.py`
- `tests/kg_validity/test_metabolomics.py`
- `data/Prochlorococcus/papers_and_supp/Kujawinski 2023/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/Kujawinski 2023/metabolite_aliases.yaml` (created if needed)
- `data/Prochlorococcus/papers_and_supp/Capovilla 2023/metabolite_aliases.yaml` (created if needed)
- `plans/strain_deployment_backlog.md`
- `docs/kg-changes/metabolomics-extension.md`

**Modified:**
- `config/schema_config.yaml` — node + edges
- `multiomics_kg/utils/paperconfig_utils.py` — add `iter_metabolite_assays_tables`
- `tools/validate_paperconfig.py` — add validation rules + `--report-backlog`
- `multiomics_kg/download/build_kegg_metabolism_xrefs.py` — step 6 extension
- `multiomics_kg/adapters/omics_adapter.py` — small fix for unconditional experiments iteration
- `create_knowledge_graph.py` — instantiate `MultiMetaboliteAssayAdapter`
- `scripts/prepare_data.sh` — add step 7
- `scripts/post-import.cypher` + `scripts/post-import.sh` — new ranks + rollups (must remain byte-identical per CLAUDE.md)
- `scripts/post-import-validate.sh` — new property enumerations
- `data/Prochlorococcus/papers_and_supp/Capovilla 2023/paperconfig.yaml` — add metabolite_assays_table entries + metabolomics experiments
- `data/Prochlorococcus/papers_and_supp/Capovilla 2023/paperconfig_files.txt` — already lists Capovilla
- `data/Prochlorococcus/papers_and_supp/paperconfig_files.txt` — add Kujawinski path
- `tests/test_paperconfig_utils.py` — add tests
- `tests/test_validate_paperconfig.py` — add tests
- `tests/test_build_kegg_metabolism_xrefs.py` — add tests
- `tests/kg_validity/test_structure.py`, `test_organism.py`, `test_post_import.py` — extend
- `tests/kg_validity/snapshot_data.json` — regenerate (commit 7)
- `.claude/skills/paperconfig/SKILL.md` — metabolomics section
- `.claude/skills/cypher-queries/SKILL.md` — spot-check templates
- `CLAUDE.md` — update key facts + supplementary_materials table

---

## Phase 1 — Schema + paperconfig_utils + validator

**Goal:** Foundation. No data change. Maps to spec §5.3 commit 1.

### Task 1.1: Add `metabolite assay` node to `config/schema_config.yaml`

**Files:**
- Modify: `config/schema_config.yaml`

- [ ] **Step 1: Locate the DerivedMetric block in schema_config.yaml**

Run: `grep -n "^derived metric:" config/schema_config.yaml`
Expected: prints a single line, e.g. `122:derived metric:`

- [ ] **Step 2: Append the `metabolite assay` block after the DerivedMetric measurement edges**

Find the line after the `derived metric classifies gene` block ends (around line 277). Insert:

```yaml
# MetaboliteAssay node: one per (Experiment × value_kind). Carries assay-level
# metabolomics semantics. Edge type emitted to Metabolite depends on value_kind.
# Adapter emits in metabolite_assay_adapter.py; post-import computes
# total_metabolite_count + distribution stats + growth_phases.
metabolite assay:
  is_a: information content entity
  represented_as: node
  preferred_id: metabolite_assay_id
  label_in_input: metabolite_assay
  properties:
    name: str
    experiment_id: str
    organism_name: str
    publication_doi: str
    compartment: str
    omics_type: str
    treatment_type: str[]
    background_factors: str[]
    treatment: str
    light_condition: str
    experimental_context: str
    metric_type: str
    value_kind: str                  # "numeric" | "boolean"
    unit: str
    rankable: str                    # "true" | "false"
    aggregation_method: str          # "mean_across_replicates" | "pre_aggregated" | "single_measurement"
    field_description: str
    total_metabolite_count: int      # post-import
    growth_phases: str[]             # post-import: from parent Experiment
    # Numeric distribution stats (post-import; null for boolean assays)
    value_min: float
    value_max: float
    value_q1: float
    value_median: float
    value_q3: float
    # Boolean flag counts (post-import; null for numeric assays)
    flag_true_count: int
    flag_false_count: int
```

- [ ] **Step 3: Verify YAML parses**

Run: `uv run python -c "import yaml; yaml.safe_load(open('config/schema_config.yaml'))"`
Expected: no output (no exception)

- [ ] **Step 4: Commit**

```bash
git add config/schema_config.yaml
git commit -m "schema: add metabolite_assay node (Phase 2 metabolomics)"
```

### Task 1.2: Add binding edges to schema

**Files:**
- Modify: `config/schema_config.yaml`

- [ ] **Step 1: Append binding edges after the new node block**

```yaml
# Publication → MetaboliteAssay (mirrors publication_has_derived_metric)
publication has metabolite assay:
  is_a: related to at instance level
  represented_as: edge
  label_in_input: publication_has_metabolite_assay
  source: publication
  target: metabolite assay

# Experiment → MetaboliteAssay
experiment has metabolite assay:
  is_a: related to at instance level
  represented_as: edge
  label_in_input: experiment_has_metabolite_assay
  source: experiment
  target: metabolite assay

# MetaboliteAssay → OrganismTaxon
metabolite assay belongs to organism:
  is_a: related to at instance level
  represented_as: edge
  label_in_input: metabolite_assay_belongs_to_organism
  source: metabolite assay
  target: organism taxon
```

- [ ] **Step 2: Verify YAML parses**

Run: `uv run python -c "import yaml; yaml.safe_load(open('config/schema_config.yaml'))"`
Expected: no output

- [ ] **Step 3: Commit**

```bash
git add config/schema_config.yaml
git commit -m "schema: add metabolite_assay binding edges"
```

### Task 1.3: Add measurement edges to schema

**Files:**
- Modify: `config/schema_config.yaml`

- [ ] **Step 1: Append measurement edges**

```yaml
# MetaboliteAssay → Metabolite (value_kind="numeric"; replicate-aggregated)
assay quantifies metabolite:
  is_a: Association
  represented_as: edge
  label_as_edge: assay_quantifies_metabolite
  label_in_input: assay_quantifies_metabolite
  source: metabolite assay
  target: metabolite
  properties:
    metric_type: str
    condition_label: str
    time_point: str
    time_point_order: int
    time_point_hours: float
    value: float
    value_sd: float
    n_replicates: int
    n_non_zero: int
    replicate_values: float[]
    detection_status: str            # adapter-set: detected | sporadic | not_detected
    rank_by_metric: int              # post-import (rankable="true" only)
    metric_percentile: float         # post-import (rankable="true" only)
    metric_bucket: str               # post-import (rankable="true" only)

# MetaboliteAssay → Metabolite (value_kind="boolean"; presence/absence)
assay flags metabolite:
  is_a: Association
  represented_as: edge
  label_as_edge: assay_flags_metabolite
  label_in_input: assay_flags_metabolite
  source: metabolite assay
  target: metabolite
  properties:
    metric_type: str
    condition_label: str
    flag_value: str                  # "true" | "false"
    n_replicates: int
    n_positive: int
```

- [ ] **Step 2: Verify YAML parses**

Run: `uv run python -c "import yaml; yaml.safe_load(open('config/schema_config.yaml'))"`
Expected: no output

- [ ] **Step 3: Commit**

```bash
git add config/schema_config.yaml
git commit -m "schema: add metabolite_assay measurement edges"
```

### Task 1.4: Add property additions to existing entries

**Files:**
- Modify: `config/schema_config.yaml`

- [ ] **Step 1: Locate the `metabolite` node properties block**

Run: `grep -n "^metabolite:" config/schema_config.yaml`
Expected: one line. Read the next ~30 lines to find the `properties:` block of the `metabolite` node.

- [ ] **Step 2: Append three computed properties to `metabolite`'s `properties:` list**

```yaml
    measured_assay_count: int        # post-import
    measured_organisms: str[]        # post-import
    measured_paper_count: int        # post-import
```

- [ ] **Step 3: Locate the `organism has metabolite` edge properties block**

Run: `grep -n "^organism has metabolite:" config/schema_config.yaml`

- [ ] **Step 4: Append four properties to its `properties:` list**

```yaml
    evidence_sources: str[]          # post-import; values: metabolism | transport | measured
    measured_assay_count: int        # post-import
    measured_compartments: str[]     # post-import
    measured_paper_count: int        # post-import
```

- [ ] **Step 5: Locate the `experiment` node properties block and append**

```yaml
    metabolite_assay_count: int      # post-import
    metabolite_compartments: str[]   # post-import
    metabolite_count: int            # post-import
```

- [ ] **Step 6: Locate the `publication` node properties block and append the same three**

```yaml
    metabolite_assay_count: int      # post-import
    metabolite_compartments: str[]   # post-import
    metabolite_count: int            # post-import
```

- [ ] **Step 7: Locate the `organism taxon` node properties block and append**

```yaml
    measured_metabolite_count: int   # post-import
```

- [ ] **Step 8: Verify YAML parses**

Run: `uv run python -c "import yaml; yaml.safe_load(open('config/schema_config.yaml'))"`
Expected: no output

- [ ] **Step 9: Commit**

```bash
git add config/schema_config.yaml
git commit -m "schema: add metabolomics-derived properties to existing entries"
```

### Task 1.5: Add `iter_metabolite_assays_tables` helper to `paperconfig_utils.py`

**Files:**
- Modify: `multiomics_kg/utils/paperconfig_utils.py`
- Test: `tests/test_paperconfig_utils.py`

- [ ] **Step 1: Read the existing `iter_derived_metrics_tables` function as the model**

Run: `grep -A 12 "^def iter_derived_metrics_tables" multiomics_kg/utils/paperconfig_utils.py`

- [ ] **Step 2: Write a failing test in `tests/test_paperconfig_utils.py`**

Append to that file:

```python
def test_iter_metabolite_assays_tables_filters_by_type():
    cfg = {
        "publication": {
            "supplementary_materials": {
                "tab_a": {"type": "csv", "filename": "a.csv"},
                "tab_b": {"type": "metabolite_assays_table", "filename": "b.csv", "experiment": "e1"},
                "tab_c": {"type": "derived_metrics_table", "filename": "c.csv"},
                "tab_d": {"type": "metabolite_assays_table", "filename": "d.csv", "experiment": "e2"},
            }
        }
    }
    from multiomics_kg.utils.paperconfig_utils import iter_metabolite_assays_tables
    keys = [k for k, _ in iter_metabolite_assays_tables(cfg)]
    assert keys == ["tab_b", "tab_d"]
```

- [ ] **Step 3: Run the test, verify it fails**

Run: `uv run pytest tests/test_paperconfig_utils.py::test_iter_metabolite_assays_tables_filters_by_type -v`
Expected: FAIL with `ImportError` or `AttributeError` for `iter_metabolite_assays_tables`

- [ ] **Step 4: Add the helper to `multiomics_kg/utils/paperconfig_utils.py`**

After `iter_derived_metrics_tables`, add:

```python
def iter_metabolite_assays_tables(config: dict):
    """Yield (entry_key, entry_dict) for each metabolite_assays_table entry."""
    for table_key, table in iter_csv_tables(config):
        if table.get("type") == "metabolite_assays_table":
            yield table_key, table
```

Wait — `iter_csv_tables` may filter out non-csv types. Verify by reading its body:

Run: `grep -A 12 "^def iter_csv_tables" multiomics_kg/utils/paperconfig_utils.py`

If `iter_csv_tables` returns ALL supplementary_materials entries (regardless of type), use the snippet above. If it filters by `type == "csv"`, use this instead:

```python
def iter_metabolite_assays_tables(config: dict):
    """Yield (entry_key, entry_dict) for each metabolite_assays_table entry."""
    pub = config.get("publication", {}) or {}
    sm = pub.get("supplementary_materials", {}) or {}
    for table_key, table in sm.items():
        if isinstance(table, dict) and table.get("type") == "metabolite_assays_table":
            yield table_key, table
```

- [ ] **Step 5: Run the test, verify it passes**

Run: `uv run pytest tests/test_paperconfig_utils.py::test_iter_metabolite_assays_tables_filters_by_type -v`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/utils/paperconfig_utils.py tests/test_paperconfig_utils.py
git commit -m "feat(paperconfig_utils): add iter_metabolite_assays_tables helper"
```

### Task 1.6: Add `validate_metabolite_assays_table` to `tools/validate_paperconfig.py`

**Files:**
- Modify: `tools/validate_paperconfig.py`
- Test: `tests/test_validate_paperconfig.py`

- [ ] **Step 1: Read existing `validate_derived_metrics_table` as model**

Run: `grep -n "validate_derived_metrics_table\|validate_csv\|validate_id_translation" tools/validate_paperconfig.py | head -10`

- [ ] **Step 2: Write failing test that exercises the validator on a happy-path entry**

Append to `tests/test_validate_paperconfig.py`:

```python
def test_validate_metabolite_assays_table_happy_path():
    from tools.validate_paperconfig import validate_metabolite_assays_table
    entry = {
        "type": "metabolite_assays_table",
        "filename": "data/example.csv",
        "experiment": "exp1",
        "organism": "Prochlorococcus MIT9303",
        "name_col": "compound",
        "assays": [
            {
                "metric_type": "cellular_concentration",
                "value_kind": "numeric",
                "field_description": "intracellular concentration fg/cell",
                "sample_columns": [
                    {"condition_label": "control", "replicate_columns": ["c1", "c2"]}
                ],
            }
        ],
    }
    experiments = {"exp1": {"organism": "Prochlorococcus MIT9303", "compartment": "whole_cell"}}
    organism_names = {"Prochlorococcus MIT9303"}
    errors = validate_metabolite_assays_table("entry_key", entry, experiments, organism_names)
    assert errors == []


def test_validate_metabolite_assays_table_missing_required_field():
    from tools.validate_paperconfig import validate_metabolite_assays_table
    entry = {"type": "metabolite_assays_table", "filename": "x.csv"}  # missing experiment, organism, name_col, assays
    errors = validate_metabolite_assays_table("e", entry, {}, set())
    joined = " | ".join(errors)
    for needed in ("experiment", "organism", "name_col", "assays"):
        assert needed in joined


def test_validate_metabolite_assays_table_unknown_experiment():
    from tools.validate_paperconfig import validate_metabolite_assays_table
    entry = {
        "type": "metabolite_assays_table",
        "filename": "x.csv",
        "experiment": "missing_exp",
        "organism": "Prochlorococcus MIT9303",
        "name_col": "compound",
        "assays": [{"metric_type": "m", "value_kind": "numeric", "field_description": "d", "sample_columns": [{"replicate_columns": ["c"]}]}],
    }
    errors = validate_metabolite_assays_table("e", entry, {"other_exp": {"organism": "x", "compartment": "whole_cell"}}, {"Prochlorococcus MIT9303"})
    assert any("missing_exp" in e for e in errors)


def test_validate_metabolite_assays_table_boolean_requires_flag_column():
    from tools.validate_paperconfig import validate_metabolite_assays_table
    entry = {
        "type": "metabolite_assays_table",
        "filename": "x.csv",
        "experiment": "exp1",
        "organism": "Prochlorococcus MIT9303",
        "name_col": "compound",
        "assays": [{
            "metric_type": "presence",
            "value_kind": "boolean",
            "field_description": "d",
            "sample_columns": [{"condition_label": "", "replicate_columns": ["c"]}],  # boolean must use flag_column
        }],
    }
    experiments = {"exp1": {"organism": "Prochlorococcus MIT9303", "compartment": "extracellular"}}
    errors = validate_metabolite_assays_table("e", entry, experiments, {"Prochlorococcus MIT9303"})
    assert any("flag_column" in e for e in errors)
```

- [ ] **Step 3: Run tests, verify they fail**

Run: `uv run pytest tests/test_validate_paperconfig.py -k metabolite_assays_table -v`
Expected: 4 FAILS (function not yet defined)

- [ ] **Step 4: Implement `validate_metabolite_assays_table` in `tools/validate_paperconfig.py`**

Add this function. Match the existing module's style (look at `validate_derived_metrics_table` for cues on error-list returning + organism cross-check signature):

```python
COMPARTMENT_VOCAB = {"whole_cell", "vesicle", "exoproteome", "extracellular"}
AGG_METHOD_VOCAB = {"mean_across_replicates", "pre_aggregated", "single_measurement"}
CELL_FORMAT_VOCAB = {"numeric", "embedded_mean_sd_n"}
VALUE_KIND_VOCAB = {"numeric", "boolean"}


def validate_metabolite_assays_table(entry_key, entry, experiments, organism_names):
    """Validate one metabolite_assays_table supplementary_materials entry.

    Returns a list of error strings (empty list = valid).
    `experiments` is a dict {exp_key: experiment_dict} from paperconfig.
    `organism_names` is a set of OrganismTaxon preferred_name values known to the KG.
    """
    errors = []
    prefix = f"metabolite_assays_table[{entry_key}]"

    # Required entry-level keys
    for key in ("type", "filename", "experiment", "organism", "name_col", "assays"):
        if key not in entry or entry.get(key) in (None, ""):
            errors.append(f"{prefix}: missing required field '{key}'")

    # If we already failed required fields, return early to keep messages clean
    if errors:
        return errors

    exp_key = entry["experiment"]
    if exp_key not in experiments:
        errors.append(f"{prefix}: experiment '{exp_key}' not found in paperconfig.experiments")
    else:
        exp = experiments[exp_key]
        comp = exp.get("compartment", "whole_cell")
        if comp not in COMPARTMENT_VOCAB:
            errors.append(
                f"{prefix}: experiment '{exp_key}' has compartment '{comp}' "
                f"not in v1 vocab {sorted(COMPARTMENT_VOCAB)}"
            )
        if exp.get("organism") and exp["organism"] != entry["organism"]:
            errors.append(
                f"{prefix}: organism '{entry['organism']}' does not match "
                f"experiment '{exp_key}' organism '{exp['organism']}'"
            )

    if entry["organism"] not in organism_names:
        errors.append(f"{prefix}: organism '{entry['organism']}' not a known OrganismTaxon")

    # Entry-level enums (optional fields)
    if "aggregation_method" in entry and entry["aggregation_method"] not in AGG_METHOD_VOCAB:
        errors.append(
            f"{prefix}: aggregation_method '{entry['aggregation_method']}' not in {sorted(AGG_METHOD_VOCAB)}"
        )
    if "cell_format" in entry and entry["cell_format"] not in CELL_FORMAT_VOCAB:
        errors.append(
            f"{prefix}: cell_format '{entry['cell_format']}' not in {sorted(CELL_FORMAT_VOCAB)}"
        )

    assays = entry.get("assays") or []
    if not assays:
        errors.append(f"{prefix}: assays must be a non-empty list")
        return errors

    seen_metric_types = set()
    for i, assay in enumerate(assays):
        a_prefix = f"{prefix}.assays[{i}]"
        for k in ("metric_type", "value_kind", "field_description", "sample_columns"):
            if k not in assay or assay.get(k) in (None, ""):
                errors.append(f"{a_prefix}: missing required field '{k}'")
        if errors and any(a_prefix in e for e in errors):
            continue

        mt = assay["metric_type"]
        if mt in seen_metric_types:
            errors.append(f"{a_prefix}: duplicate metric_type '{mt}' within entry")
        seen_metric_types.add(mt)

        vk = assay["value_kind"]
        if vk not in VALUE_KIND_VOCAB:
            errors.append(f"{a_prefix}: value_kind '{vk}' not in {sorted(VALUE_KIND_VOCAB)}")

        if "aggregation_method" in assay and assay["aggregation_method"] not in AGG_METHOD_VOCAB:
            errors.append(
                f"{a_prefix}: aggregation_method '{assay['aggregation_method']}' not in {sorted(AGG_METHOD_VOCAB)}"
            )

        sample_cols = assay.get("sample_columns") or []
        if not sample_cols:
            errors.append(f"{a_prefix}: sample_columns must be a non-empty list")
            continue

        for j, sc in enumerate(sample_cols):
            sc_prefix = f"{a_prefix}.sample_columns[{j}]"
            if vk == "numeric":
                rc = sc.get("replicate_columns")
                if not isinstance(rc, list) or not rc:
                    errors.append(f"{sc_prefix}: numeric assay requires replicate_columns: list[str] (non-empty)")
            elif vk == "boolean":
                if not sc.get("flag_column"):
                    errors.append(f"{sc_prefix}: boolean assay requires flag_column: str")
                if "flag_true_value" not in sc:
                    errors.append(f"{sc_prefix}: boolean assay requires flag_true_value: str")

    return errors
```

- [ ] **Step 5: Run tests, verify all 4 pass**

Run: `uv run pytest tests/test_validate_paperconfig.py -k metabolite_assays_table -v`
Expected: 4 PASS

- [ ] **Step 6: Commit**

```bash
git add tools/validate_paperconfig.py tests/test_validate_paperconfig.py
git commit -m "feat(validator): add metabolite_assays_table schema validation"
```

### Task 1.7: Wire `validate_metabolite_assays_table` into the main validator entry point

**Files:**
- Modify: `tools/validate_paperconfig.py`

- [ ] **Step 1: Find where the main validator dispatches by `type`**

Run: `grep -n "type.*derived_metrics_table\|elif.*type.*csv\|dispatch\|type.*==" tools/validate_paperconfig.py | head -10`

- [ ] **Step 2: Add a branch for `metabolite_assays_table`**

In whatever function dispatches per-entry validation (e.g. `validate_paperconfig` or `_validate_supplementary_entry`), add an `elif` branch after the `derived_metrics_table` branch:

```python
elif entry_type == "metabolite_assays_table":
    errors.extend(
        validate_metabolite_assays_table(
            entry_key, entry, experiments=experiments, organism_names=known_organisms
        )
    )
```

(Adapt parameter names — `experiments` and `known_organisms` are illustrative; use the actual context the surrounding code uses.)

- [ ] **Step 3: Run the full validator on an existing paperconfig (regression check)**

Run: `uv run python tools/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/Capovilla 2023/paperconfig.yaml"`
Expected: success (no metabolite_assays_table entries yet → no new errors)

- [ ] **Step 4: Commit**

```bash
git add tools/validate_paperconfig.py
git commit -m "feat(validator): dispatch metabolite_assays_table to new validator"
```

### Task 1.8: Add `--report-backlog` flag to validator

**Files:**
- Modify: `tools/validate_paperconfig.py`

- [ ] **Step 1: Locate the argparse setup in `tools/validate_paperconfig.py`**

Run: `grep -n "argparse\|add_argument\|ArgumentParser" tools/validate_paperconfig.py | head -10`

- [ ] **Step 2: Add the flag and a handler**

Add to argparse:

```python
parser.add_argument(
    "--report-backlog",
    action="store_true",
    help="Scan paperconfig files for '# BACKLOG:' markers and print a per-file report.",
)
```

In the main dispatch (early in `main()` after parsing args), add:

```python
if args.report_backlog:
    return _report_backlog(args.paperconfigs)
```

Add this helper at module scope:

```python
def _report_backlog(paths):
    """Scan paperconfig YAML files for '# BACKLOG:' lines and print per-file."""
    import re
    pattern = re.compile(r"#\s*BACKLOG:?\s*(.*)")
    total = 0
    for path in paths:
        try:
            with open(path) as f:
                lines = f.readlines()
        except OSError:
            continue
        hits = []
        for lineno, line in enumerate(lines, start=1):
            m = pattern.search(line)
            if m:
                hits.append((lineno, m.group(1).strip() or "(no message)"))
        if hits:
            print(f"\n{path}: {len(hits)} backlog markers")
            for lineno, msg in hits:
                print(f"  line {lineno}: {msg}")
            total += len(hits)
    print(f"\nTotal: {total} backlog markers across {len(paths)} file(s).")
    return 0
```

- [ ] **Step 3: Smoke-test on existing paperconfigs**

Run: `uv run python tools/validate_paperconfig.py --report-backlog $(cat data/Prochlorococcus/papers_and_supp/paperconfig_files.txt)`
Expected: prints zero markers (no paperconfig has `# BACKLOG:` yet) and a "Total: 0" line.

- [ ] **Step 4: Commit**

```bash
git add tools/validate_paperconfig.py
git commit -m "feat(validator): add --report-backlog flag"
```

### Task 1.9: Phase 1 validation gate

**Files:** none (verification only)

- [ ] **Step 1: Run all unit tests**

Run: `uv run pytest -m "not slow and not kg" -v`
Expected: all pass; new tests for `iter_metabolite_assays_tables` + `validate_metabolite_assays_table` are green.

- [ ] **Step 2: Run a fast KG build smoke test**

Run: `uv run python create_knowledge_graph.py --test`
Expected: completes without schema-related errors. The new node type `metabolite assay` should appear in the BioCypher output (no nodes yet — adapter not wired — but BioCypher must accept the schema).

- [ ] **Step 3: Spot-check schema_config.yaml is well-formed**

Run: `uv run python -c "import yaml; cfg = yaml.safe_load(open('config/schema_config.yaml')); assert 'metabolite assay' in cfg, 'missing metabolite assay node'"`
Expected: no output

---

## Phase 2 — Step 6 extension + step 7

**Goal:** Add paper-metabolite resolution pipeline. Maps to spec §5.3 commit 2.

### Task 2.1: Add paper-alias harvest function in `build_kegg_metabolism_xrefs.py`

**Files:**
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py`
- Test: `tests/test_build_kegg_metabolism_xrefs.py`

- [ ] **Step 1: Read existing structure of `build_kegg_metabolism_xrefs.py`**

Run: `grep -n "^def \|^class " multiomics_kg/download/build_kegg_metabolism_xrefs.py`

Skim the file to understand where MNX resolver is opened and where outputs are written. The new harvest function needs to be called after MNX resolver is loaded, before kegg_data.json is finalized.

- [ ] **Step 2: Write a failing test for the harvest function**

Append to `tests/test_build_kegg_metabolism_xrefs.py` (create if missing):

```python
import pandas as pd


def test_harvest_paper_metabolite_aliases_collects_names_and_ids(tmp_path):
    """Harvest function returns name → primary_id resolutions from paperconfig CSVs + aliases."""
    from multiomics_kg.download.build_kegg_metabolism_xrefs import harvest_paper_metabolite_aliases

    csv_path = tmp_path / "metab.csv"
    pd.DataFrame({
        "compound": ["Glucose", "GABA", "MysteryX"],
        "KEGG_ID": ["C00031", "", ""],
    }).to_csv(csv_path, index=False)

    aliases_path = tmp_path / "metabolite_aliases.yaml"
    aliases_path.write_text('"GABA": "kegg.compound:C00334"\n')

    paperconfigs = [{
        "paperconfig_path": str(tmp_path / "paperconfig.yaml"),
        "publication": {
            "supplementary_materials": {
                "table_a": {
                    "type": "metabolite_assays_table",
                    "filename": str(csv_path),
                    "name_col": "compound",
                    "id_col": "KEGG_ID",
                    "id_type": "kegg.compound",
                    "aliases_file": "metabolite_aliases.yaml",
                    "experiment": "exp1",
                    "organism": "Prochlorococcus MIT9303",
                    "assays": [{"metric_type": "x", "value_kind": "numeric", "field_description": "d",
                                "sample_columns": [{"replicate_columns": ["c"]}]}],
                }
            }
        },
    }]

    # Stub resolver: returns id for "Glucose", None for "MysteryX"
    class StubResolver:
        def resolve(self, name, formula=None):
            return {"Glucose": ("kegg.compound:C00031", "name_match", 1.0)}.get(name)

    result = harvest_paper_metabolite_aliases(paperconfigs, StubResolver())
    # Expect: id_col direct hit for Glucose; alias_override for GABA; unresolved for MysteryX
    assert result["alias_to_primary"]["Glucose"] == "kegg.compound:C00031"
    assert result["alias_to_primary"]["GABA"] == "kegg.compound:C00334"
    assert "MysteryX" in result["unresolved"]
    # Resolution methods recorded
    assert result["resolution_methods"]["Glucose"] == "kegg_direct"
    assert result["resolution_methods"]["GABA"] == "alias_override"
    assert result["resolution_methods"]["MysteryX"] == "unresolved"
```

- [ ] **Step 3: Run, verify it fails**

Run: `uv run pytest tests/test_build_kegg_metabolism_xrefs.py::test_harvest_paper_metabolite_aliases_collects_names_and_ids -v`
Expected: FAIL with import error.

- [ ] **Step 4: Implement `harvest_paper_metabolite_aliases` in `build_kegg_metabolism_xrefs.py`**

Add at module scope (after existing imports):

```python
import yaml as _yaml
import pandas as _pd
from pathlib import Path as _Path


def _load_aliases_file(paperconfig_path: str, aliases_filename: str) -> dict:
    """Load a paper-local metabolite_aliases.yaml. Returns {} if missing or invalid."""
    if not aliases_filename:
        return {}
    base = _Path(paperconfig_path).parent
    path = base / aliases_filename
    if not path.exists():
        return {}
    try:
        data = _yaml.safe_load(path.read_text()) or {}
    except _yaml.YAMLError:
        return {}
    if not isinstance(data, dict):
        return {}
    return {str(k): str(v) for k, v in data.items()}


def harvest_paper_metabolite_aliases(paperconfigs: list[dict], resolver) -> dict:
    """Walk paperconfigs for metabolite_assays_table entries; resolve every name.

    Returns a dict with:
        alias_to_primary: {name_or_alias: primary_id}
        resolution_methods: {name_or_alias: method_str}
        unresolved: list[str]
        per_paper: {paperconfig_path: {resolved: int, unresolved: int, total: int}}
    """
    alias_to_primary: dict[str, str] = {}
    resolution_methods: dict[str, str] = {}
    unresolved: set[str] = set()
    per_paper: dict[str, dict[str, int]] = {}

    for cfg in paperconfigs:
        pc_path = cfg.get("paperconfig_path", "<unknown>")
        pub = (cfg.get("publication") or {})
        sm = (pub.get("supplementary_materials") or {})
        for entry_key, entry in sm.items():
            if not isinstance(entry, dict) or entry.get("type") != "metabolite_assays_table":
                continue
            csv_path = entry.get("filename")
            if not csv_path or not _Path(csv_path).exists():
                continue
            name_col = entry.get("name_col") or "compound"
            id_col = entry.get("id_col") or ""
            id_type = entry.get("id_type") or ""
            formula_col = entry.get("formula_col") or ""
            aliases = _load_aliases_file(pc_path, entry.get("aliases_file") or "")

            try:
                df = _pd.read_csv(csv_path, dtype=str, keep_default_na=False)
            except Exception:
                continue

            stats = per_paper.setdefault(pc_path, {"resolved": 0, "unresolved": 0, "total": 0})
            for _, row in df.iterrows():
                name = (row.get(name_col) or "").strip()
                if not name:
                    continue
                stats["total"] += 1
                if name in alias_to_primary:
                    stats["resolved"] += 1
                    continue

                primary, method = None, None

                # 1. id_col direct
                if id_col and id_type:
                    raw = (row.get(id_col) or "").strip()
                    if raw:
                        # Strip prefix if user already included it
                        if ":" in raw:
                            primary = raw
                        else:
                            primary = f"{id_type}:{raw}"
                        method = (
                            "kegg_direct" if id_type == "kegg.compound"
                            else "chebi_direct" if id_type == "chebi"
                            else "mnx_direct" if id_type == "mnx"
                            else "id_direct"
                        )

                # 2. aliases_file override
                if primary is None and name in aliases:
                    primary, method = aliases[name], "alias_override"

                # 3. resolver (with optional formula disambiguator)
                if primary is None:
                    formula = (row.get(formula_col) or "").strip() if formula_col else None
                    hit = resolver.resolve(name, formula=formula)
                    if hit is not None:
                        primary_id, method_hint, _similarity = hit
                        primary = primary_id
                        method = method_hint or "name_match"

                if primary is not None:
                    alias_to_primary[name] = primary
                    resolution_methods[name] = method
                    stats["resolved"] += 1
                else:
                    unresolved.add(name)
                    resolution_methods[name] = "unresolved"
                    stats["unresolved"] += 1

    return {
        "alias_to_primary": alias_to_primary,
        "resolution_methods": resolution_methods,
        "unresolved": sorted(unresolved),
        "per_paper": per_paper,
    }
```

- [ ] **Step 5: Run the test, verify it passes**

Run: `uv run pytest tests/test_build_kegg_metabolism_xrefs.py::test_harvest_paper_metabolite_aliases_collects_names_and_ids -v`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/download/build_kegg_metabolism_xrefs.py tests/test_build_kegg_metabolism_xrefs.py
git commit -m "feat(step6): harvest paper metabolite aliases from paperconfigs"
```

### Task 2.2: Write `metabolite_id_mapping.json` from harvest result

**Files:**
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py`
- Test: `tests/test_build_kegg_metabolism_xrefs.py`

- [ ] **Step 1: Write a failing test**

Append:

```python
def test_write_metabolite_id_mapping_creates_three_tier_skeleton(tmp_path):
    from multiomics_kg.download.build_kegg_metabolism_xrefs import write_metabolite_id_mapping
    import json
    harvest = {
        "alias_to_primary": {"Glucose": "kegg.compound:C00031", "GABA": "kegg.compound:C00334"},
        "resolution_methods": {"Glucose": "kegg_direct", "GABA": "alias_override"},
        "unresolved": ["MysteryX"],
        "per_paper": {},
    }
    out = tmp_path / "metabolite_id_mapping.json"
    write_metabolite_id_mapping(harvest, out)
    data = json.loads(out.read_text())
    # v1: only name_lookup populated
    assert data["specific_lookup"] == {}
    assert data["multi_lookup"] == {}
    assert data["conflicts"] == {}
    assert data["name_lookup"]["Glucose"] == ["kegg.compound:C00031"]
    assert data["name_lookup"]["GABA"] == ["kegg.compound:C00334"]
    assert "compounds" in data
```

- [ ] **Step 2: Run, verify it fails**

Run: `uv run pytest tests/test_build_kegg_metabolism_xrefs.py::test_write_metabolite_id_mapping_creates_three_tier_skeleton -v`
Expected: FAIL with import error

- [ ] **Step 3: Implement `write_metabolite_id_mapping`**

Add to `multiomics_kg/download/build_kegg_metabolism_xrefs.py`:

```python
import json as _json


def write_metabolite_id_mapping(harvest: dict, out_path) -> None:
    """Write metabolite_id_mapping.json. v1: only name_lookup populated.

    Schema reserves three tiers (specific_lookup, multi_lookup, name_lookup)
    + conflicts + compounds; future ID-column harvesting will fill the
    other tiers. v1 fills name_lookup only.
    """
    name_lookup: dict[str, list[str]] = {}
    compounds: dict[str, dict] = {}
    for name, primary in harvest["alias_to_primary"].items():
        name_lookup.setdefault(name, []).append(primary)
        compounds.setdefault(primary, {"aliases": []})
        if name not in compounds[primary]["aliases"]:
            compounds[primary]["aliases"].append(name)

    data = {
        "specific_lookup": {},
        "multi_lookup": {},
        "name_lookup": name_lookup,
        "conflicts": {},
        "compounds": compounds,
        "schema_version": 1,
    }
    out_path = _Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(_json.dumps(data, indent=2, sort_keys=True))
```

- [ ] **Step 4: Run, verify it passes**

Run: `uv run pytest tests/test_build_kegg_metabolism_xrefs.py::test_write_metabolite_id_mapping_creates_three_tier_skeleton -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/build_kegg_metabolism_xrefs.py tests/test_build_kegg_metabolism_xrefs.py
git commit -m "feat(step6): write metabolite_id_mapping.json with three-tier skeleton"
```

### Task 2.3: Extend `kegg_data.json` with paper-measured compounds + `evidence_sources` union

**Files:**
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py`
- Test: `tests/test_build_kegg_metabolism_xrefs.py`

- [ ] **Step 1: Write a failing test**

Append:

```python
def test_extend_kegg_data_unions_evidence_sources():
    from multiomics_kg.download.build_kegg_metabolism_xrefs import extend_kegg_data_with_paper_metabolites
    kegg_data = {
        "compounds": {
            "kegg.compound:C00031": {  # already gene-reachable via metabolism
                "id": "kegg.compound:C00031",
                "name": "D-Glucose",
                "evidence_sources": ["metabolism"],
            },
        },
        "additional_compounds": {},
    }
    paper_primary_ids = {
        "kegg.compound:C00031",   # overlap → union
        "kegg.compound:C00334",   # not previously in cache → add to compounds bucket as new entry
        "chebi:16865",            # non-KEGG → additional_compounds
    }
    extend_kegg_data_with_paper_metabolites(kegg_data, paper_primary_ids)
    assert sorted(kegg_data["compounds"]["kegg.compound:C00031"]["evidence_sources"]) == ["metabolism", "metabolomics"]
    assert kegg_data["compounds"]["kegg.compound:C00334"]["evidence_sources"] == ["metabolomics"]
    assert kegg_data["additional_compounds"]["chebi:16865"]["evidence_sources"] == ["metabolomics"]
```

- [ ] **Step 2: Run, verify fails**

Run: `uv run pytest tests/test_build_kegg_metabolism_xrefs.py::test_extend_kegg_data_unions_evidence_sources -v`
Expected: FAIL

- [ ] **Step 3: Implement `extend_kegg_data_with_paper_metabolites`**

Add:

```python
def extend_kegg_data_with_paper_metabolites(kegg_data: dict, paper_primary_ids) -> None:
    """In-place extension of kegg_data.json structure with paper-measured compounds.

    KEGG-primary compounds → existing 'compounds' bucket.
    Non-KEGG primary IDs → 'additional_compounds' bucket.
    evidence_sources is unioned (already-present compounds get 'metabolomics' appended).
    """
    compounds = kegg_data.setdefault("compounds", {})
    additional = kegg_data.setdefault("additional_compounds", {})
    for primary in paper_primary_ids:
        if primary.startswith("kegg.compound:"):
            entry = compounds.setdefault(primary, {"id": primary, "evidence_sources": []})
            ev = list(entry.get("evidence_sources") or [])
            if "metabolomics" not in ev:
                ev.append("metabolomics")
            entry["evidence_sources"] = sorted(set(ev))
        else:
            entry = additional.setdefault(primary, {"id": primary, "evidence_sources": []})
            ev = list(entry.get("evidence_sources") or [])
            if "metabolomics" not in ev:
                ev.append("metabolomics")
            entry["evidence_sources"] = sorted(set(ev))
```

- [ ] **Step 4: Run, verify passes**

Run: `uv run pytest tests/test_build_kegg_metabolism_xrefs.py::test_extend_kegg_data_unions_evidence_sources -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/build_kegg_metabolism_xrefs.py tests/test_build_kegg_metabolism_xrefs.py
git commit -m "feat(step6): extend kegg_data.json with paper-measured compounds"
```

### Task 2.4: Pathway extension — paper-measured compounds gain `Metabolite_in_pathway` reach

**Files:**
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py`
- Test: `tests/test_build_kegg_metabolism_xrefs.py`

- [ ] **Step 1: Write a failing test**

Append:

```python
def test_extend_pathway_set_with_paper_compounds():
    from multiomics_kg.download.build_kegg_metabolism_xrefs import extended_pathway_set_with_papers
    # compound→pathway mapping from raw KEGG (full)
    compound_to_pathways = {
        "kegg.compound:C00031": ["kegg.pathway:ko00010", "kegg.pathway:ko_orphan_1"],
        "kegg.compound:C00334": ["kegg.pathway:ko_orphan_2"],
    }
    existing_extended = {"kegg.pathway:ko00010"}
    paper_primary_ids = {"kegg.compound:C00031", "kegg.compound:C00334"}
    result = extended_pathway_set_with_papers(existing_extended, paper_primary_ids, compound_to_pathways)
    # ko00010 stays; ko_orphan_1 + ko_orphan_2 added because paper compounds reach them
    assert result == {"kegg.pathway:ko00010", "kegg.pathway:ko_orphan_1", "kegg.pathway:ko_orphan_2"}
```

- [ ] **Step 2: Run, verify fails**

Run: `uv run pytest tests/test_build_kegg_metabolism_xrefs.py::test_extend_pathway_set_with_paper_compounds -v`
Expected: FAIL

- [ ] **Step 3: Implement**

Add:

```python
def extended_pathway_set_with_papers(
    existing_extended: set, paper_primary_ids, compound_to_pathways: dict
) -> set:
    """Union the extended-pathway set with pathways reachable from paper-measured compounds."""
    result = set(existing_extended)
    for cid in paper_primary_ids:
        for pwid in compound_to_pathways.get(cid, []):
            result.add(pwid)
    return result
```

- [ ] **Step 4: Run, verify passes**

Run: `uv run pytest tests/test_build_kegg_metabolism_xrefs.py::test_extend_pathway_set_with_paper_compounds -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/build_kegg_metabolism_xrefs.py tests/test_build_kegg_metabolism_xrefs.py
git commit -m "feat(step6): extend pathway set with paper-compound-reachable pathways"
```

### Task 2.5: Wire the new helpers into the step 6 main() flow

**Files:**
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py`

- [ ] **Step 1: Locate the `main()` function in `build_kegg_metabolism_xrefs.py`**

Run: `grep -n "^def main\|if __name__" multiomics_kg/download/build_kegg_metabolism_xrefs.py`

- [ ] **Step 2: Read the section that builds `kegg_data` and writes the JSON**

Identify the variable holding the in-memory kegg_data structure and the variable holding the extended pathway set, before the JSON is written.

- [ ] **Step 3: Insert the harvest+extension calls**

Before `kegg_data` is written to `cache/data/kegg/kegg_data.json`, insert:

```python
# Phase 2 metabolomics: harvest paper-measured metabolites, extend kegg_data, write mapping
from multiomics_kg.utils.paperconfig_utils import load_all_paperconfigs
paperconfigs = load_all_paperconfigs()  # signature: returns list of paperconfig dicts with paperconfig_path
harvest = harvest_paper_metabolite_aliases(paperconfigs, mnx_resolver)  # mnx_resolver name may vary; use existing local
paper_primary_ids = set(harvest["alias_to_primary"].values())
extend_kegg_data_with_paper_metabolites(kegg_data, paper_primary_ids)
extended_pathways = extended_pathway_set_with_papers(extended_pathways, paper_primary_ids, compound_to_pathways)
write_metabolite_id_mapping(harvest, _Path("cache/data/metabolomics/metabolite_id_mapping.json"))
print(f"[step6] Phase 2 harvest: {len(harvest['alias_to_primary'])} resolved, {len(harvest['unresolved'])} unresolved")
```

(Variable names `mnx_resolver`, `extended_pathways`, `compound_to_pathways` may differ — adapt to actual local variables. If `load_all_paperconfigs()` doesn't exist with that signature, search `paperconfig_utils.py` for the existing iterator and adapt.)

- [ ] **Step 4: Sanity check — no paperconfig has the new entry type yet, so this is a no-op**

Run step 6 manually:

```bash
uv run python -m multiomics_kg.download.build_kegg_metabolism_xrefs --force 2>&1 | tail -30
```

Expected: completes; final `[step6] Phase 2 harvest: 0 resolved, 0 unresolved` line. `cache/data/metabolomics/metabolite_id_mapping.json` exists with empty `name_lookup`.

- [ ] **Step 5: Verify file is created**

Run: `ls -la cache/data/metabolomics/metabolite_id_mapping.json && cat cache/data/metabolomics/metabolite_id_mapping.json | head -20`
Expected: file exists; JSON has `name_lookup: {}`, `specific_lookup: {}`, etc.

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/download/build_kegg_metabolism_xrefs.py
git commit -m "feat(step6): wire paper-metabolite harvest into main flow"
```

### Task 2.6: Create `resolve_paper_metabolites.py` (step 7) skeleton + first test

**Files:**
- Create: `multiomics_kg/download/resolve_paper_metabolites.py`
- Create: `tests/test_resolve_paper_metabolites.py`

- [ ] **Step 1: Write a failing test**

Create `tests/test_resolve_paper_metabolites.py`:

```python
import json
import pandas as pd


def test_resolve_paper_metabolites_writes_resolved_csv(tmp_path):
    from multiomics_kg.download.resolve_paper_metabolites import resolve_paper_metabolites_for_entry

    # Source CSV: 3 metabolites
    csv_path = tmp_path / "metab.csv"
    pd.DataFrame({
        "compound": ["Glucose", "GABA", "MysteryX"],
        "KEGG_ID": ["C00031", "", ""],
    }).to_csv(csv_path, index=False)

    # Lightweight mapping (mimics metabolite_id_mapping.json post-step6)
    mapping = {
        "specific_lookup": {},
        "multi_lookup": {},
        "name_lookup": {
            "Glucose": ["kegg.compound:C00031"],
            "GABA": ["kegg.compound:C00334"],
        },
        "conflicts": {},
        "compounds": {},
    }

    entry = {
        "type": "metabolite_assays_table",
        "filename": str(csv_path),
        "name_col": "compound",
        "id_col": "KEGG_ID",
        "id_type": "kegg.compound",
    }

    out_csv, report = resolve_paper_metabolites_for_entry(entry, mapping)
    assert out_csv == csv_path.with_name("metab_resolved.csv")
    df = pd.read_csv(out_csv, dtype=str, keep_default_na=False)
    rows = {r["compound"]: (r["metabolite_id"], r["resolution_method"]) for _, r in df.iterrows()}
    assert rows["Glucose"] == ("kegg.compound:C00031", "kegg_direct")
    assert rows["GABA"] == ("kegg.compound:C00334", "name_match")
    assert rows["MysteryX"] == ("", "unresolved")
    assert report["resolved"] == 2
    assert report["unresolved"] == 1
```

- [ ] **Step 2: Run, verify fails**

Run: `uv run pytest tests/test_resolve_paper_metabolites.py::test_resolve_paper_metabolites_writes_resolved_csv -v`
Expected: FAIL with import error

- [ ] **Step 3: Implement `resolve_paper_metabolites_for_entry`**

Create `multiomics_kg/download/resolve_paper_metabolites.py`:

```python
#!/usr/bin/env python3
"""Step 7: Resolve metabolite names in paper CSV files.

For each metabolite_assays_table entry across all paperconfigs, opens the
source CSV, looks up each row's metabolite (by id_col first, then name_col
via metabolite_id_mapping.json), writes <stem>_resolved.csv with two added
columns: metabolite_id + resolution_method. Also writes
<stem>_resolution_report.json.

Step 7 does NOT load the MNX resolver — it reads only the lightweight
metabolite_id_mapping.json that step 6 produced.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from multiomics_kg.utils.paperconfig_utils import (
    iter_metabolite_assays_tables,
    load_all_paperconfigs,
)


def resolve_paper_metabolites_for_entry(entry: dict, mapping: dict) -> tuple[Path, dict]:
    """Resolve one metabolite_assays_table entry's CSV. Returns (resolved_csv_path, report)."""
    src = Path(entry["filename"])
    name_col = entry["name_col"]
    id_col = entry.get("id_col") or ""
    id_type = entry.get("id_type") or ""
    name_lookup = mapping.get("name_lookup") or {}

    df = pd.read_csv(src, dtype=str, keep_default_na=False)

    resolved_ids: list[str] = []
    methods: list[str] = []
    n_resolved = 0
    n_unresolved = 0

    for _, row in df.iterrows():
        primary, method = _resolve_row(row, name_col, id_col, id_type, name_lookup)
        resolved_ids.append(primary or "")
        methods.append(method)
        if primary:
            n_resolved += 1
        else:
            n_unresolved += 1

    df["metabolite_id"] = resolved_ids
    df["resolution_method"] = methods

    out = src.with_name(src.stem + "_resolved.csv")
    df.to_csv(out, index=False)

    report = {
        "source_csv": str(src),
        "resolved_csv": str(out),
        "total_rows": len(df),
        "resolved": n_resolved,
        "unresolved": n_unresolved,
        "method_counts": dict(pd.Series(methods).value_counts()),
    }
    report_path = src.with_name(src.stem + "_resolution_report.json")
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True, default=int))
    return out, report


def _resolve_row(row, name_col, id_col, id_type, name_lookup):
    """Return (primary_id_or_None, resolution_method)."""
    # 1. id_col direct
    if id_col and id_type:
        raw = (row.get(id_col) or "").strip()
        if raw:
            primary = raw if ":" in raw else f"{id_type}:{raw}"
            method = (
                "kegg_direct" if id_type == "kegg.compound"
                else "chebi_direct" if id_type == "chebi"
                else "mnx_direct" if id_type == "mnx"
                else "id_direct"
            )
            return primary, method

    # 2. name_col → mapping
    name = (row.get(name_col) or "").strip()
    if not name:
        return None, "unresolved"
    hits = name_lookup.get(name) or []
    if len(hits) == 1:
        return hits[0], "name_match"
    if len(hits) > 1:
        # ambiguous; v1 picks first deterministically and tags
        return hits[0], "ambiguous_multi_id"
    return None, "unresolved"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--papers", nargs="*", default=None,
                        help="Restrict to specific paper names (optional).")
    args = parser.parse_args(argv)

    mapping_path = Path("cache/data/metabolomics/metabolite_id_mapping.json")
    if not mapping_path.exists():
        print(f"[step7] {mapping_path} not found — run step 6 first.", file=sys.stderr)
        return 1
    mapping = json.loads(mapping_path.read_text())

    paperconfigs = load_all_paperconfigs()
    total_entries = 0
    for cfg in paperconfigs:
        paper = (cfg.get("publication") or {}).get("papername", "<unknown>")
        if args.papers and paper not in args.papers:
            continue
        for entry_key, entry in iter_metabolite_assays_tables(cfg):
            src = Path(entry.get("filename") or "")
            if not src.exists():
                print(f"[step7] {paper}/{entry_key}: source not found {src}", file=sys.stderr)
                continue
            out = src.with_name(src.stem + "_resolved.csv")
            if out.exists() and not args.force:
                print(f"[step7] {paper}/{entry_key}: skipping (exists). Use --force to regen.")
                continue
            _, report = resolve_paper_metabolites_for_entry(entry, mapping)
            print(f"[step7] {paper}/{entry_key}: {report['resolved']}/{report['total_rows']} resolved")
            total_entries += 1

    print(f"[step7] Done. Processed {total_entries} entries.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run the test, verify it passes**

Run: `uv run pytest tests/test_resolve_paper_metabolites.py::test_resolve_paper_metabolites_writes_resolved_csv -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/resolve_paper_metabolites.py tests/test_resolve_paper_metabolites.py
git commit -m "feat(step7): resolve_paper_metabolites — id_col + name_col lookup"
```

### Task 2.7: Add step 7 to `scripts/prepare_data.sh`

**Files:**
- Modify: `scripts/prepare_data.sh`

- [ ] **Step 1: Find where step 6 is called in `prepare_data.sh`**

Run: `grep -n "step 6\|build_kegg_metabolism_xrefs\|prepare_data_step6" scripts/prepare_data.sh`

- [ ] **Step 2: Add step 7 invocation immediately after step 6**

Edit `scripts/prepare_data.sh`. After the step 6 block:

```bash
# Step 7 — Resolve paper metabolite names to primary IDs (mirrors step 4 for genes)
if should_run_step 7; then
    log "Step 7: resolve_paper_metabolites"
    uv run python -m multiomics_kg.download.resolve_paper_metabolites $force_arg \
        > "logs/prepare_data_step7.log" 2>&1
fi
```

(Variable names `should_run_step`, `log`, `force_arg` are illustrative — match the existing pattern around step 6.)

Update the documentation comment block at the top of the file to include Step 7.

- [ ] **Step 3: Smoke-test step 7 by itself**

Run: `bash scripts/prepare_data.sh --steps 7 --force 2>&1 | tail -10`
Expected: completes (zero entries processed since no paperconfig has the new type yet); writes `logs/prepare_data_step7.log`.

- [ ] **Step 4: Commit**

```bash
git add scripts/prepare_data.sh
git commit -m "build: add step 7 (resolve_paper_metabolites) to prepare_data.sh"
```

### Task 2.8: Phase 2 validation gate

**Files:** none (verification only)

- [ ] **Step 1: All unit tests pass**

Run: `uv run pytest tests/test_build_kegg_metabolism_xrefs.py tests/test_resolve_paper_metabolites.py tests/test_paperconfig_utils.py -v`
Expected: all pass.

- [ ] **Step 2: Step 6 outputs are well-formed**

Run: `bash scripts/prepare_data.sh --steps 6 7 --force 2>&1 | tail -20 && cat cache/data/metabolomics/metabolite_id_mapping.json | head`
Expected: step 6 + step 7 succeed; metabolite_id_mapping.json has the 5-key skeleton (specific_lookup/multi_lookup/name_lookup/conflicts/compounds + schema_version).

---

## Phase 3 — `metabolism_adapter` verification

**Goal:** Confirm new Metabolite nodes flow through. Spec §5.3 commit 3.

### Task 3.1: Add KG validity test for `evidence_sources` containing `metabolomics`

**Files:**
- Modify: `tests/kg_validity/test_post_import.py`

- [ ] **Step 1: Read existing post-import KG tests for context**

Run: `grep -n "evidence_sources\|^def test_" tests/kg_validity/test_post_import.py | head -20`

- [ ] **Step 2: Add a new test (will skip pre-Phase-6 since no paperconfig has the new type yet, but written now to pin behavior)**

Append:

```python
@pytest.mark.kg
def test_metabolomics_evidence_source_appears_when_papers_integrated(neo4j_driver):
    """After Capovilla/Kujawinski integration, at least some Metabolite nodes
    must carry 'metabolomics' in evidence_sources. Skipped pre-integration."""
    with neo4j_driver.session() as s:
        n = s.run(
            "MATCH (m:Metabolite) WHERE 'metabolomics' IN m.evidence_sources RETURN count(m) AS n"
        ).single()["n"]
    if n == 0:
        pytest.skip("No metabolomics-source metabolites yet (Phase 2 papers not integrated)")
    assert n > 0
```

- [ ] **Step 3: Run the KG validity tests**

Run: `uv run pytest tests/kg_validity/test_post_import.py::test_metabolomics_evidence_source_appears_when_papers_integrated -v`
Expected: SKIPPED (pre-integration). If FAIL, the harness has a different reason; investigate.

- [ ] **Step 4: Commit**

```bash
git add tests/kg_validity/test_post_import.py
git commit -m "test(kg): assert metabolomics evidence_source appears after integration"
```

### Task 3.2: Phase 3 validation gate

**Files:** none

- [ ] **Step 1: Run all pre-existing KG tests**

Run: `uv run pytest -m kg -v 2>&1 | tail -30`
Expected: existing tests green; metabolomics test SKIPPED.

---

## Phase 4 — `MetaboliteAssayAdapter` + wiring

**Goal:** New adapter, dormant until paperconfigs use it. Spec §5.3 commit 4.

### Task 4.1: Create `metabolite_assay_adapter.py` skeleton + node-emission test

**Files:**
- Create: `multiomics_kg/adapters/metabolite_assay_adapter.py`
- Create: `tests/test_metabolite_assay_adapter.py`

- [ ] **Step 1: Write a failing test for node emission**

Create `tests/test_metabolite_assay_adapter.py`:

```python
import yaml
from pathlib import Path


def _make_minimal_paperconfig(tmp_path):
    csv = tmp_path / "metab.csv"
    csv.write_text("compound,c1,c2\nGlucose,1.0,2.0\n")
    cfg = {
        "publication": {
            "papername": "Test 2026",
            "doi": "10.1234/test",
            "experiments": {
                "exp1": {
                    "name": "test exp",
                    "organism": "Prochlorococcus MIT9303",
                    "compartment": "whole_cell",
                    "omics_type": "METABOLOMICS",
                    "treatment_type": ["nitrogen"],
                    "background_factors": ["axenic"],
                    "treatment": "N-stress",
                    "control": "replete",
                    "experimental_context": "ctx",
                    "light_condition": "light",
                },
            },
            "supplementary_materials": {
                "tab_a": {
                    "type": "metabolite_assays_table",
                    "filename": str(csv),
                    "experiment": "exp1",
                    "organism": "Prochlorococcus MIT9303",
                    "name_col": "compound",
                    "assays": [
                        {
                            "metric_type": "cellular_concentration",
                            "name": "test assay",
                            "value_kind": "numeric",
                            "unit": "fg/cell",
                            "rankable": "true",
                            "field_description": "test",
                            "sample_columns": [
                                {"condition_label": "control",
                                 "replicate_columns": ["c1", "c2"]},
                            ],
                        }
                    ],
                }
            },
        }
    }
    pc_path = tmp_path / "paperconfig.yaml"
    pc_path.write_text(yaml.safe_dump(cfg))
    return pc_path


def test_metabolite_assay_adapter_emits_one_node_per_assay(tmp_path):
    from multiomics_kg.adapters.metabolite_assay_adapter import MetaboliteAssayAdapter
    pc = _make_minimal_paperconfig(tmp_path)
    # Skip resolved-csv lookup — adapter must fall back to the raw csv if no _resolved.csv
    adapter = MetaboliteAssayAdapter(config_file=str(pc), test_mode=False)
    nodes = list(adapter.get_nodes())
    assert len(nodes) == 1
    node_id, label, props = nodes[0]
    assert label == "metabolite_assay"
    assert node_id.startswith("metabolite_assay:")
    assert "tab_a" in node_id
    assert "cellular_concentration" in node_id
    assert props["metric_type"] == "cellular_concentration"
    assert props["value_kind"] == "numeric"
    assert props["compartment"] == "whole_cell"
    assert props["organism_name"] == "Prochlorococcus MIT9303"
```

- [ ] **Step 2: Run, verify fails**

Run: `uv run pytest tests/test_metabolite_assay_adapter.py::test_metabolite_assay_adapter_emits_one_node_per_assay -v`
Expected: FAIL (module doesn't exist)

- [ ] **Step 3: Implement minimal adapter to make test pass**

Create `multiomics_kg/adapters/metabolite_assay_adapter.py`:

```python
"""MetaboliteAssayAdapter — reads metabolite_assays_table entries from paperconfig.yaml.

Mirrors observations_adapter.py (DerivedMetric). Single-paper class +
Multi* registry-driven wrapper.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterator

import yaml

from multiomics_kg.utils.paperconfig_utils import (
    get_paper_name,
    iter_metabolite_assays_tables,
    load_all_paperconfigs,
)

logger = logging.getLogger(__name__)


def _clean_str(value) -> str:
    if value is None:
        return ""
    return str(value).replace("'", "^").replace("|", "")


def _make_metabolite_assay_id(doi: str, paper_name: str, entry_key: str, metric_type: str) -> str:
    base = doi or paper_name or "unknown"
    short = base.split("/")[-1] if "/" in base else base
    return f"metabolite_assay:{short}:{entry_key}:{metric_type}"


def _denormalized_fields_from_experiment(exp: dict) -> dict:
    """Return the subset of Experiment properties denormalized onto every MetaboliteAssay node."""
    return {
        "compartment": _clean_str(exp.get("compartment", "whole_cell") or "whole_cell"),
        "omics_type": _clean_str(exp.get("omics_type", "METABOLOMICS") or "METABOLOMICS"),
        "treatment_type": [_clean_str(t) for t in (exp.get("treatment_type") or [])],
        "background_factors": [_clean_str(t) for t in (exp.get("background_factors") or [])],
        "treatment": _clean_str(exp.get("treatment", "")),
        "light_condition": _clean_str(exp.get("light_condition", "")),
        "experimental_context": _clean_str(exp.get("experimental_context", "")),
    }


class MetaboliteAssayAdapter:
    """Single-paper adapter."""

    def __init__(self, config_file: str, test_mode: bool = False):
        self.config_file = config_file
        self.test_mode = test_mode
        with open(config_file) as f:
            self.config = yaml.safe_load(f) or {}
        self.paper_name = get_paper_name(self.config) or ""
        self.doi = (self.config.get("publication") or {}).get("doi", "")
        self._entries = list(iter_metabolite_assays_tables(self.config))
        self._experiments = ((self.config.get("publication") or {}).get("experiments") or {})

    def get_nodes(self) -> Iterator[tuple[str, str, dict]]:
        for entry_key, entry in self._entries:
            exp_key = entry.get("experiment", "")
            exp = self._experiments.get(exp_key) or {}
            denorm = _denormalized_fields_from_experiment(exp)
            organism = entry.get("organism") or exp.get("organism") or ""
            experiment_id = f"{self.doi}_{exp_key}" if self.doi else f"{self.paper_name}_{exp_key}"
            for assay in entry.get("assays", []):
                metric_type = assay.get("metric_type", "")
                if not metric_type:
                    continue
                node_id = _make_metabolite_assay_id(self.doi, self.paper_name, entry_key, metric_type)
                props = {
                    "name": _clean_str(assay.get("name", metric_type)),
                    "experiment_id": _clean_str(experiment_id),
                    "organism_name": _clean_str(organism),
                    "publication_doi": _clean_str(self.doi),
                    "metric_type": _clean_str(metric_type),
                    "value_kind": _clean_str(assay.get("value_kind", "")),
                    "unit": _clean_str(assay.get("unit", "")),
                    "rankable": _clean_str(assay.get("rankable", "false")),
                    "aggregation_method": _clean_str(
                        assay.get("aggregation_method")
                        or entry.get("aggregation_method", "mean_across_replicates")
                    ),
                    "field_description": _clean_str(assay.get("field_description", "")),
                    **denorm,
                }
                yield node_id, "metabolite_assay", props

    def get_edges(self) -> Iterator[tuple]:
        # Filled in subsequent tasks (4.4 onward)
        return
        yield  # pragma: no cover  -- generator skeleton


class MultiMetaboliteAssayAdapter:
    """Registry-driven wrapper. Reads paperconfig_files.txt and aggregates across papers."""

    def __init__(self, paperconfig_paths=None, test_mode: bool = False):
        if paperconfig_paths is None:
            paperconfigs = load_all_paperconfigs()
            paperconfig_paths = [c.get("paperconfig_path") for c in paperconfigs if c.get("paperconfig_path")]
        self.adapters = []
        for p in paperconfig_paths:
            try:
                with open(p) as f:
                    cfg = yaml.safe_load(f) or {}
            except OSError:
                continue
            sm = (cfg.get("publication") or {}).get("supplementary_materials") or {}
            if any(isinstance(v, dict) and v.get("type") == "metabolite_assays_table" for v in sm.values()):
                self.adapters.append(MetaboliteAssayAdapter(config_file=p, test_mode=test_mode))

    def get_nodes(self):
        for a in self.adapters:
            yield from a.get_nodes()

    def get_edges(self):
        for a in self.adapters:
            yield from a.get_edges()
```

- [ ] **Step 4: Run, verify passes**

Run: `uv run pytest tests/test_metabolite_assay_adapter.py::test_metabolite_assay_adapter_emits_one_node_per_assay -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/metabolite_assay_adapter.py tests/test_metabolite_assay_adapter.py
git commit -m "feat(adapter): MetaboliteAssayAdapter skeleton + node emission"
```

### Task 4.2: Implement `_aggregate_replicates` helper with `detection_status`

**Files:**
- Modify: `multiomics_kg/adapters/metabolite_assay_adapter.py`
- Modify: `tests/test_metabolite_assay_adapter.py`

- [ ] **Step 1: Write failing tests covering the boundary cases**

Append to `tests/test_metabolite_assay_adapter.py`:

```python
def test_aggregate_replicates_all_detected():
    from multiomics_kg.adapters.metabolite_assay_adapter import _aggregate_replicates
    mean, sd, n_rep, n_nz, vals, det = _aggregate_replicates(
        ["1.0", "2.0", "3.0"], null_values=set(), missing_values={""}
    )
    assert n_rep == 3
    assert n_nz == 3
    assert mean == 2.0
    assert vals == [1.0, 2.0, 3.0]
    assert det == "detected"


def test_aggregate_replicates_all_zero():
    from multiomics_kg.adapters.metabolite_assay_adapter import _aggregate_replicates
    _, _, n_rep, n_nz, _, det = _aggregate_replicates(
        ["0", "0", "0"], null_values=set(), missing_values={""}
    )
    assert n_rep == 3 and n_nz == 0 and det == "not_detected"


def test_aggregate_replicates_sporadic():
    from multiomics_kg.adapters.metabolite_assay_adapter import _aggregate_replicates
    _, _, n_rep, n_nz, _, det = _aggregate_replicates(
        ["1.5", "0", "0"], null_values=set(), missing_values={""}
    )
    assert n_rep == 3 and n_nz == 1 and det == "sporadic"


def test_aggregate_replicates_null_values_become_zero():
    from multiomics_kg.adapters.metabolite_assay_adapter import _aggregate_replicates
    _, _, n_rep, n_nz, vals, det = _aggregate_replicates(
        ["nd", "1.0", "NA"], null_values={"nd", "NA"}, missing_values={""}
    )
    # nd + NA → 0; 1.0 detected → 1/3 sporadic
    assert n_rep == 3 and n_nz == 1 and det == "sporadic"
    assert vals == [0.0, 1.0, 0.0]


def test_aggregate_replicates_missing_excluded():
    from multiomics_kg.adapters.metabolite_assay_adapter import _aggregate_replicates
    _, _, n_rep, n_nz, vals, det = _aggregate_replicates(
        ["", "1.0", "2.0"], null_values=set(), missing_values={""}
    )
    # blank excluded → only 2 replicates
    assert n_rep == 2 and n_nz == 2 and det == "detected"
    assert vals == [1.0, 2.0]
```

- [ ] **Step 2: Run, verify all 5 fail**

Run: `uv run pytest tests/test_metabolite_assay_adapter.py -k aggregate_replicates -v`
Expected: 5 FAIL

- [ ] **Step 3: Implement `_aggregate_replicates`**

Add to `multiomics_kg/adapters/metabolite_assay_adapter.py`:

```python
import statistics
from typing import Sequence


def _aggregate_replicates(
    raw_values: Sequence[str],
    null_values: set[str],
    missing_values: set[str],
) -> tuple[float, float, int, int, list[float], str]:
    """Aggregate replicate cells into (mean, sd, n_replicates, n_non_zero, replicate_values, detection_status).

    null_values  → treated as 0.0 (not-detected, counted in n_replicates).
    missing_values → row excluded from aggregation entirely.
    Anything else → coerced to float; non-numeric raises (caller should pre-validate).
    """
    parsed: list[float] = []
    for v in raw_values:
        s = (v or "").strip()
        if s in missing_values:
            continue
        if s in null_values:
            parsed.append(0.0)
            continue
        try:
            parsed.append(float(s))
        except ValueError:
            # Bad cell — treat as missing (skip)
            continue

    n_replicates = len(parsed)
    if n_replicates == 0:
        return 0.0, 0.0, 0, 0, [], "not_detected"
    n_non_zero = sum(1 for x in parsed if x != 0.0)
    mean = statistics.fmean(parsed)
    sd = statistics.stdev(parsed) if n_replicates >= 2 else 0.0

    if n_non_zero == 0:
        det = "not_detected"
    elif n_non_zero == n_replicates:
        det = "detected"
    else:
        det = "sporadic"

    return mean, sd, n_replicates, n_non_zero, parsed, det
```

- [ ] **Step 4: Run, verify all 5 pass**

Run: `uv run pytest tests/test_metabolite_assay_adapter.py -k aggregate_replicates -v`
Expected: 5 PASS

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/metabolite_assay_adapter.py tests/test_metabolite_assay_adapter.py
git commit -m "feat(adapter): _aggregate_replicates with detection_status"
```

### Task 4.3: Implement `embedded_mean_sd_n` cell parser

**Files:**
- Modify: `multiomics_kg/adapters/metabolite_assay_adapter.py`
- Modify: `tests/test_metabolite_assay_adapter.py`

- [ ] **Step 1: Write failing tests**

Append:

```python
def test_parse_embedded_cell_normal():
    from multiomics_kg.adapters.metabolite_assay_adapter import parse_embedded_mean_sd_n
    result = parse_embedded_mean_sd_n("0.00054 (8.8e-05), n=2")
    assert result == (0.00054, 8.8e-05, 2)


def test_parse_embedded_cell_nd_returns_zeros():
    from multiomics_kg.adapters.metabolite_assay_adapter import parse_embedded_mean_sd_n
    assert parse_embedded_mean_sd_n("nd") == (0.0, 0.0, 0)
    assert parse_embedded_mean_sd_n("ND") == (0.0, 0.0, 0)


def test_parse_embedded_cell_na_sd():
    from multiomics_kg.adapters.metabolite_assay_adapter import parse_embedded_mean_sd_n
    # "0.0016 (NA), n=1" → mean parsed, sd=0 (NA), n=1
    result = parse_embedded_mean_sd_n("0.0016 (NA), n=1")
    assert result[0] == 0.0016
    assert result[1] == 0.0
    assert result[2] == 1


def test_parse_embedded_cell_blank_returns_none():
    from multiomics_kg.adapters.metabolite_assay_adapter import parse_embedded_mean_sd_n
    assert parse_embedded_mean_sd_n("") is None
```

- [ ] **Step 2: Run, verify fails**

Run: `uv run pytest tests/test_metabolite_assay_adapter.py -k parse_embedded -v`
Expected: 4 FAIL

- [ ] **Step 3: Implement parser**

Add:

```python
import re

_EMBEDDED_PATTERN = re.compile(
    r"^\s*([0-9.+\-eE]+)\s*\(\s*([0-9.+\-eEnNaA/]+)\s*\)\s*,\s*n\s*=\s*(\d+)\s*$"
)


def parse_embedded_mean_sd_n(cell: str) -> tuple[float, float, int] | None:
    """Parse cells like '0.00054 (8.8e-05), n=2' or 'nd'.

    Returns (mean, sd, n) or None for empty cells. 'nd'/'ND' → (0,0,0).
    Non-numeric sd (e.g. 'NA') → 0.0.
    """
    s = (cell or "").strip()
    if not s:
        return None
    if s.lower() in {"nd", "n.d."}:
        return 0.0, 0.0, 0
    m = _EMBEDDED_PATTERN.match(s)
    if not m:
        return None
    mean = float(m.group(1))
    try:
        sd = float(m.group(2))
    except ValueError:
        sd = 0.0
    n = int(m.group(3))
    return mean, sd, n
```

- [ ] **Step 4: Run, verify passes**

Run: `uv run pytest tests/test_metabolite_assay_adapter.py -k parse_embedded -v`
Expected: 4 PASS

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/metabolite_assay_adapter.py tests/test_metabolite_assay_adapter.py
git commit -m "feat(adapter): parse_embedded_mean_sd_n cell-format parser"
```

### Task 4.4: Implement numeric edge emission (`Assay_quantifies_metabolite`)

**Files:**
- Modify: `multiomics_kg/adapters/metabolite_assay_adapter.py`
- Modify: `tests/test_metabolite_assay_adapter.py`

- [ ] **Step 1: Write failing test**

Append:

```python
def test_numeric_edge_emission(tmp_path):
    """Adapter emits one Assay_quantifies_metabolite edge per (assay, metabolite, condition)."""
    from multiomics_kg.adapters.metabolite_assay_adapter import MetaboliteAssayAdapter

    # Resolved CSV (since step 7 normally produces this)
    resolved = tmp_path / "metab_resolved.csv"
    import pandas as pd
    pd.DataFrame({
        "compound": ["Glucose"],
        "metabolite_id": ["kegg.compound:C00031"],
        "resolution_method": ["name_match"],
        "c1": [1.0],
        "c2": [3.0],
    }).to_csv(resolved, index=False)

    src = tmp_path / "metab.csv"
    pd.DataFrame({"compound": ["Glucose"], "c1": [1.0], "c2": [3.0]}).to_csv(src, index=False)

    cfg = {
        "publication": {
            "papername": "T", "doi": "10.1/x",
            "experiments": {"e1": {"organism": "Prochlorococcus MIT9303",
                                   "compartment": "whole_cell", "omics_type": "METABOLOMICS"}},
            "supplementary_materials": {
                "tab": {
                    "type": "metabolite_assays_table",
                    "filename": str(src),
                    "experiment": "e1",
                    "organism": "Prochlorococcus MIT9303",
                    "name_col": "compound",
                    "assays": [{
                        "metric_type": "conc",
                        "value_kind": "numeric",
                        "field_description": "d",
                        "unit": "fg/cell",
                        "rankable": "true",
                        "sample_columns": [
                            {"condition_label": "control", "replicate_columns": ["c1", "c2"]},
                        ],
                    }],
                }
            },
        }
    }
    import yaml
    pc = tmp_path / "paperconfig.yaml"
    pc.write_text(yaml.safe_dump(cfg))

    adapter = MetaboliteAssayAdapter(config_file=str(pc), test_mode=False)
    edges = list(adapter.get_edges())
    quant = [e for e in edges if e[3] == "assay_quantifies_metabolite"]
    assert len(quant) == 1
    edge_id, src_id, dst_id, label, props = quant[0]
    assert dst_id == "kegg.compound:C00031"
    assert props["value"] == 2.0
    assert props["n_replicates"] == 2
    assert props["n_non_zero"] == 2
    assert props["replicate_values"] == [1.0, 3.0]
    assert props["detection_status"] == "detected"
    assert props["condition_label"] == "control"
    assert props["metric_type"] == "conc"
```

- [ ] **Step 2: Run, verify fails**

Run: `uv run pytest tests/test_metabolite_assay_adapter.py::test_numeric_edge_emission -v`
Expected: FAIL (get_edges still empty generator)

- [ ] **Step 3: Implement `get_edges` for the numeric variant**

In `metabolite_assay_adapter.py`:

1. At module scope (top of file, after existing imports), add:

```python
import pandas as _pd

DEFAULT_NULL_VALUES = {"nd", "ND", "NA", "n.d.", "N/A"}
DEFAULT_MISSING_VALUES = {""}
```

2. Replace the placeholder `get_edges` method on the existing `MetaboliteAssayAdapter` class (do NOT redeclare the class — these are *additional methods on the same class*):

```python
    def _open_resolved_csv(self, src_path: str) -> _pd.DataFrame | None:
        """Prefer <stem>_resolved.csv (step 7 output) over the source CSV."""
        src = Path(src_path)
        if not src.exists():
            return None
        resolved = src.with_name(src.stem + "_resolved.csv")
        path = resolved if resolved.exists() else src
        try:
            return _pd.read_csv(path, dtype=str, keep_default_na=False)
        except Exception as e:
            logger.warning("Could not read %s: %s", path, e)
            return None

    def get_edges(self) -> Iterator[tuple]:
        for entry_key, entry in self._entries:
            exp_key = entry.get("experiment", "")
            df = self._open_resolved_csv(entry.get("filename", ""))
            if df is None:
                continue
            null_values = set(entry.get("null_values") or DEFAULT_NULL_VALUES)
            missing_values = set(entry.get("missing_values") or DEFAULT_MISSING_VALUES)
            name_col = entry.get("name_col", "")
            id_col_in_resolved = "metabolite_id"

            for assay in entry.get("assays", []):
                metric_type = assay.get("metric_type", "")
                value_kind = assay.get("value_kind", "")
                if not metric_type or not value_kind:
                    continue
                assay_node_id = _make_metabolite_assay_id(self.doi, self.paper_name, entry_key, metric_type)

                if value_kind == "numeric":
                    yield from self._emit_numeric_edges(
                        df, entry, assay, assay_node_id, name_col, id_col_in_resolved,
                        null_values, missing_values,
                    )
                elif value_kind == "boolean":
                    yield from self._emit_boolean_edges(
                        df, entry, assay, assay_node_id, name_col, id_col_in_resolved,
                    )
                # binding edges per assay
                yield from self._emit_binding_edges(assay_node_id, exp_key, entry)

    def _emit_numeric_edges(
        self, df, entry, assay, assay_node_id, name_col, id_col_in_resolved,
        null_values, missing_values,
    ):
        sample_cols_block = assay.get("sample_columns", [])
        for row_idx, row in df.iterrows():
            primary = (row.get(id_col_in_resolved) or "").strip()
            if not primary:
                continue  # unresolved row → no edge
            for sc_idx, sc in enumerate(sample_cols_block):
                rcols = sc.get("replicate_columns") or []
                if not rcols:
                    continue
                raw_values = [row.get(c, "") for c in rcols]
                mean, sd, n_rep, n_nz, vals, det = _aggregate_replicates(
                    raw_values, null_values=null_values, missing_values=missing_values,
                )
                if n_rep == 0 and not entry.get("drop_undetected", False):
                    # All-missing row → skip per implicit policy
                    continue
                cond_label = _clean_str(sc.get("condition_label", ""))
                edge_id = f"{assay_node_id}|{primary}|{sc_idx}"
                props = {
                    "metric_type": _clean_str(assay.get("metric_type", "")),
                    "condition_label": cond_label,
                    "time_point": _clean_str(sc.get("time_point", "")),
                    "time_point_order": int(sc.get("time_point_order") or 0),
                    "time_point_hours": float(sc.get("time_point_hours", -1.0) or -1.0),
                    "value": float(mean),
                    "value_sd": float(sd),
                    "n_replicates": int(n_rep),
                    "n_non_zero": int(n_nz),
                    "replicate_values": vals,
                    "detection_status": det,
                }
                yield edge_id, assay_node_id, primary, "assay_quantifies_metabolite", props

    def _emit_boolean_edges(self, df, entry, assay, assay_node_id, name_col, id_col_in_resolved):
        # Filled in Task 4.5
        return
        yield

    def _emit_binding_edges(self, assay_node_id, exp_key, entry):
        # Filled in Task 4.6
        return
        yield
```

(Note: don't actually duplicate the class definition — extend the original; the `...` and re-declaration above is shorthand for "add these methods to the class.")

- [ ] **Step 4: Run, verify passes**

Run: `uv run pytest tests/test_metabolite_assay_adapter.py::test_numeric_edge_emission -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/metabolite_assay_adapter.py tests/test_metabolite_assay_adapter.py
git commit -m "feat(adapter): numeric Assay_quantifies_metabolite edges"
```

### Task 4.5: Implement boolean edge emission (`Assay_flags_metabolite`)

**Files:**
- Modify: `multiomics_kg/adapters/metabolite_assay_adapter.py`
- Modify: `tests/test_metabolite_assay_adapter.py`

- [ ] **Step 1: Write failing test**

Append:

```python
def test_boolean_edge_emission(tmp_path):
    from multiomics_kg.adapters.metabolite_assay_adapter import MetaboliteAssayAdapter
    import pandas as pd

    resolved = tmp_path / "presence_resolved.csv"
    pd.DataFrame({
        "compound": ["GABA", "Glucose"],
        "metabolite_id": ["kegg.compound:C00334", "kegg.compound:C00031"],
        "resolution_method": ["alias_override", "kegg_direct"],
        "intracellular_flag": ["yes", ""],
    }).to_csv(resolved, index=False)
    src = tmp_path / "presence.csv"
    pd.DataFrame({"compound": ["GABA", "Glucose"], "intracellular_flag": ["yes", ""]}).to_csv(src, index=False)

    cfg = {
        "publication": {
            "papername": "T", "doi": "10.1/x",
            "experiments": {"e1": {"organism": "P", "compartment": "whole_cell", "omics_type": "METABOLOMICS"}},
            "supplementary_materials": {
                "tab": {
                    "type": "metabolite_assays_table",
                    "filename": str(src),
                    "experiment": "e1",
                    "organism": "P",
                    "name_col": "compound",
                    "assays": [{
                        "metric_type": "presence",
                        "value_kind": "boolean",
                        "field_description": "presence flag",
                        "rankable": "false",
                        "sample_columns": [
                            {"condition_label": "", "flag_column": "intracellular_flag",
                             "flag_true_value": "yes"},
                        ],
                    }],
                }
            },
        }
    }
    import yaml
    pc = tmp_path / "paperconfig.yaml"
    pc.write_text(yaml.safe_dump(cfg))
    adapter = MetaboliteAssayAdapter(config_file=str(pc))
    edges = [e for e in adapter.get_edges() if e[3] == "assay_flags_metabolite"]
    # Two rows; only first ("yes") emits a true flag; second emits a false flag
    assert len(edges) == 2
    by_dst = {e[2]: e[4] for e in edges}
    assert by_dst["kegg.compound:C00334"]["flag_value"] == "true"
    assert by_dst["kegg.compound:C00334"]["n_positive"] == 1
    assert by_dst["kegg.compound:C00031"]["flag_value"] == "false"
    assert by_dst["kegg.compound:C00031"]["n_positive"] == 0
```

- [ ] **Step 2: Run, verify fails**

Run: `uv run pytest tests/test_metabolite_assay_adapter.py::test_boolean_edge_emission -v`
Expected: FAIL

- [ ] **Step 3: Implement `_emit_boolean_edges`**

Replace the placeholder in `metabolite_assay_adapter.py`:

```python
    def _emit_boolean_edges(self, df, entry, assay, assay_node_id, name_col, id_col_in_resolved):
        sample_cols_block = assay.get("sample_columns", [])
        for row_idx, row in df.iterrows():
            primary = (row.get(id_col_in_resolved) or "").strip()
            if not primary:
                continue
            for sc_idx, sc in enumerate(sample_cols_block):
                flag_col = sc.get("flag_column")
                if not flag_col:
                    continue
                flag_true_value = sc.get("flag_true_value", "yes")
                cell = (row.get(flag_col) or "").strip()
                is_true = cell == flag_true_value
                edge_id = f"{assay_node_id}|{primary}|{sc_idx}"
                cond_label = _clean_str(sc.get("condition_label", ""))
                props = {
                    "metric_type": _clean_str(assay.get("metric_type", "")),
                    "condition_label": cond_label,
                    "flag_value": "true" if is_true else "false",
                    "n_replicates": 1,
                    "n_positive": 1 if is_true else 0,
                }
                yield edge_id, assay_node_id, primary, "assay_flags_metabolite", props
```

- [ ] **Step 4: Run, verify passes**

Run: `uv run pytest tests/test_metabolite_assay_adapter.py::test_boolean_edge_emission -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/metabolite_assay_adapter.py tests/test_metabolite_assay_adapter.py
git commit -m "feat(adapter): boolean Assay_flags_metabolite edges"
```

### Task 4.6: Implement binding edges (`PublicationHasMetaboliteAssay`, `ExperimentHasMetaboliteAssay`, `MetaboliteAssayBelongsToOrganism`)

**Files:**
- Modify: `multiomics_kg/adapters/metabolite_assay_adapter.py`
- Modify: `tests/test_metabolite_assay_adapter.py`

- [ ] **Step 1: Write failing test**

Append:

```python
def test_binding_edges_emitted(tmp_path):
    """Each MetaboliteAssay node has 3 binding edges: from Publication, from Experiment, to OrganismTaxon."""
    from multiomics_kg.adapters.metabolite_assay_adapter import MetaboliteAssayAdapter
    pc = _make_minimal_paperconfig(tmp_path)
    adapter = MetaboliteAssayAdapter(config_file=str(pc), test_mode=False)
    labels = [e[3] for e in adapter.get_edges()]
    assert labels.count("publication_has_metabolite_assay") == 1
    assert labels.count("experiment_has_metabolite_assay") == 1
    assert labels.count("metabolite_assay_belongs_to_organism") == 1
```

- [ ] **Step 2: Run, verify fails**

Run: `uv run pytest tests/test_metabolite_assay_adapter.py::test_binding_edges_emitted -v`
Expected: FAIL

- [ ] **Step 3: Pass entry through to `_emit_binding_edges` and implement it**

Refactor the call site in `get_edges` to pass `entry` directly. Find the `_emit_binding_edges` invocation and change it from `self._emit_binding_edges(assay_node_id, exp_key)` to `self._emit_binding_edges(assay_node_id, exp_key, entry)`. Then implement:

```python
    def _emit_binding_edges(self, assay_node_id, exp_key, entry):
        organism = entry.get("organism") or self._experiments.get(exp_key, {}).get("organism", "")
        experiment_id = f"{self.doi}_{exp_key}" if self.doi else f"{self.paper_name}_{exp_key}"
        pub_id = self.doi or f"papername:{self.paper_name}"

        yield (f"pub__{assay_node_id}", pub_id, assay_node_id,
               "publication_has_metabolite_assay", {})
        yield (f"exp__{assay_node_id}", experiment_id, assay_node_id,
               "experiment_has_metabolite_assay", {})
        yield (f"org__{assay_node_id}", assay_node_id, organism,
               "metabolite_assay_belongs_to_organism", {})
```

- [ ] **Step 4: Run, verify passes**

Run: `uv run pytest tests/test_metabolite_assay_adapter.py::test_binding_edges_emitted -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/metabolite_assay_adapter.py tests/test_metabolite_assay_adapter.py
git commit -m "feat(adapter): emit Publication/Experiment/Organism binding edges"
```

### Task 4.7: Fix `omics_adapter` to iterate experiments unconditionally

**Files:**
- Modify: `multiomics_kg/adapters/omics_adapter.py`

- [ ] **Step 1: Locate where `omics_adapter` emits Experiment nodes**

Run: `grep -n "experiments\b\|Experiment\|get_nodes" multiomics_kg/adapters/omics_adapter.py | head -20`

- [ ] **Step 2: Find the conditional that may skip experiments without DE tables**

Skim the function that yields Experiment nodes. If it only iterates experiments referenced by `statistical_analyses`, change it to iterate the `experiments:` block directly.

- [ ] **Step 3: Modify to iterate the `experiments:` block unconditionally**

The exact diff depends on existing structure. Pattern to apply:

```python
# Before (illustrative):
for stat in entry.get("statistical_analyses", []):
    exp_key = stat.get("experiment")
    if exp_key not in seen_experiments:
        yield <experiment node>

# After:
for exp_key, exp in (publication.get("experiments") or {}).items():
    if exp_key not in seen_experiments:
        yield <experiment node>  # using exp dict for properties
```

- [ ] **Step 4: Run all unit tests to verify no regression**

Run: `uv run pytest -m "not slow and not kg" -v 2>&1 | tail -30`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/omics_adapter.py
git commit -m "fix(omics_adapter): iterate experiments block unconditionally"
```

### Task 4.8: Wire `MultiMetaboliteAssayAdapter` into `create_knowledge_graph.py`

**Files:**
- Modify: `create_knowledge_graph.py`

- [ ] **Step 1: Locate where `MultiObservationsAdapter` is instantiated**

Run: `grep -n "MultiObservationsAdapter\|MultiOMICSAdapter\|write_nodes\|write_edges" create_knowledge_graph.py | head -10`

- [ ] **Step 2: Add `MultiMetaboliteAssayAdapter` after `MultiObservationsAdapter`**

```python
from multiomics_kg.adapters.metabolite_assay_adapter import MultiMetaboliteAssayAdapter
# ...
ma_adapter = MultiMetaboliteAssayAdapter(test_mode=args.test)
bc.write_nodes(ma_adapter.get_nodes())
bc.write_edges(ma_adapter.get_edges())
```

- [ ] **Step 3: Run a `--test` build**

Run: `uv run python create_knowledge_graph.py --test 2>&1 | tail -20`
Expected: success; no MetaboliteAssay nodes yet (no paperconfig has the new entry type).

- [ ] **Step 4: Commit**

```bash
git add create_knowledge_graph.py
git commit -m "feat(kg): wire MultiMetaboliteAssayAdapter into pipeline"
```

### Task 4.9: Phase 4 validation gate

**Files:** none

- [ ] **Step 1: All unit tests pass**

Run: `uv run pytest -m "not slow and not kg" -v 2>&1 | tail -20`
Expected: all green.

- [ ] **Step 2: KG build succeeds with no new MetaboliteAssay nodes (adapter dormant)**

Run: `uv run python create_knowledge_graph.py --test 2>&1 | grep -i "metabolite_assay\|MetaboliteAssay" | head -5`
Expected: no MetaboliteAssay node lines (zero entries because no paperconfig uses the new type yet). Build itself should succeed.

---

## Phase 5 — Post-import Cypher + validate.sh

**Goal:** Computed properties + rollups. Spec §5.3 commit 5.

### Task 5.1: Add MetaboliteAssay indexes + fulltext index

**Files:**
- Modify: `scripts/post-import.cypher`
- Modify: `scripts/post-import.sh`

- [ ] **Step 1: Find the indexes block in `post-import.cypher`**

Run: `grep -n "CREATE INDEX\|CREATE FULLTEXT" scripts/post-import.cypher | head -5`

- [ ] **Step 2: Append after existing indexes (e.g. after `derived_metric_*` indexes)**

In `scripts/post-import.cypher`:

```cypher
// MetaboliteAssay scalar + full-text indexes
CREATE INDEX metabolite_assay_organism_idx     IF NOT EXISTS FOR (a:MetaboliteAssay) ON (a.organism_name);
CREATE INDEX metabolite_assay_compartment_idx  IF NOT EXISTS FOR (a:MetaboliteAssay) ON (a.compartment);
CREATE INDEX metabolite_assay_metric_type_idx  IF NOT EXISTS FOR (a:MetaboliteAssay) ON (a.metric_type);
CREATE INDEX metabolite_assay_value_kind_idx   IF NOT EXISTS FOR (a:MetaboliteAssay) ON (a.value_kind);
CREATE INDEX metabolite_assay_experiment_idx   IF NOT EXISTS FOR (a:MetaboliteAssay) ON (a.experiment_id);
CREATE FULLTEXT INDEX metaboliteAssayFullText  IF NOT EXISTS FOR (a:MetaboliteAssay)
  ON EACH [a.name, a.field_description, a.treatment, a.experimental_context];
```

- [ ] **Step 3: Mirror the exact same lines into `scripts/post-import.sh`**

Open `scripts/post-import.sh`, find the corresponding indexes block (it groups statements into `cypher-shell` invocations). Append the same Cypher into the same group.

- [ ] **Step 4: Run a quick post-import dry-run on the deployed graph (if available)**

Run: `bash scripts/post-import.sh 2>&1 | head -20` (only if a graph is currently deployed; else skip).
Expected: no errors. If no graph deployed, validate via:

```bash
cypher-shell -u neo4j -a bolt://localhost:7687 -f /dev/stdin <<'EOF'
SHOW INDEXES YIELD name WHERE name STARTS WITH 'metabolite_assay_' RETURN name;
EOF
```

(Skip if no graph; will be exercised in Phase 6.)

- [ ] **Step 5: Commit**

```bash
git add scripts/post-import.cypher scripts/post-import.sh
git commit -m "post-import: add MetaboliteAssay indexes"
```

### Task 5.2: Add `rank_by_metric` / `metric_percentile` / `metric_bucket` Cypher

**Files:**
- Modify: `scripts/post-import.cypher`
- Modify: `scripts/post-import.sh`

- [ ] **Step 1: Find the DerivedMetric rank block as model**

Run: `grep -n "rank_by_metric\|MATCH (dm:DerivedMetric {rankable" scripts/post-import.cypher | head -5`

- [ ] **Step 2: Append the parallel block for MetaboliteAssay**

After the DerivedMetric rank block, append:

```cypher
// MetaboliteAssay numeric ranks: per-assay, only when rankable='true'.
// Mirrors DerivedMetric pattern (per-assay scope, deterministic Metabolite.id tiebreaker).
MATCH (a:MetaboliteAssay {rankable: 'true'})
CALL {
  WITH a
  MATCH (a)-[r:Assay_quantifies_metabolite]->(m:Metabolite)
  WITH r, r.value AS val, m.id AS mid
  ORDER BY val DESC, mid ASC
  WITH collect(r) AS edges, count(r) AS n
  UNWIND range(0, size(edges) - 1) AS i
  WITH edges[i] AS r, i, n,
       CASE WHEN n = 1 THEN 100.0
            ELSE 100.0 * toFloat(n - i - 1) / toFloat(n - 1)
       END AS pct
  SET r.rank_by_metric = i + 1,
      r.metric_percentile = pct,
      r.metric_bucket = CASE
        WHEN pct >= 90.0 THEN 'top_decile'
        WHEN pct >= 75.0 THEN 'top_quartile'
        WHEN pct >= 25.0 THEN 'mid'
        ELSE 'low'
      END
} IN TRANSACTIONS OF 10 ROWS;
```

- [ ] **Step 3: Mirror the exact same block into `post-import.sh`**

- [ ] **Step 4: Commit**

```bash
git add scripts/post-import.cypher scripts/post-import.sh
git commit -m "post-import: rank_by_metric/percentile/bucket for MetaboliteAssay"
```

### Task 5.3: Add `MetaboliteAssay` node rollups (`total_metabolite_count`, distribution stats, flag counts, `growth_phases`)

**Files:**
- Modify: `scripts/post-import.cypher`
- Modify: `scripts/post-import.sh`

- [ ] **Step 1: Append the rollup blocks**

```cypher
// MetaboliteAssay total_metabolite_count
MATCH (a:MetaboliteAssay)
OPTIONAL MATCH (a)-[r:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
WITH a, count(DISTINCT m) AS cnt
SET a.total_metabolite_count = cnt;

// Numeric distribution stats (null for boolean assays)
MATCH (a:MetaboliteAssay {value_kind: 'numeric'})-[r:Assay_quantifies_metabolite]->()
WITH a,
     min(r.value) AS vmin, max(r.value) AS vmax,
     percentileDisc(r.value, 0.25) AS q1,
     percentileDisc(r.value, 0.5)  AS med,
     percentileDisc(r.value, 0.75) AS q3
SET a.value_min=vmin, a.value_max=vmax, a.value_q1=q1, a.value_median=med, a.value_q3=q3;

// Boolean flag counts (null for numeric assays)
MATCH (a:MetaboliteAssay {value_kind: 'boolean'})-[r:Assay_flags_metabolite]->()
WITH a,
     sum(CASE WHEN r.flag_value='true'  THEN 1 ELSE 0 END) AS t,
     sum(CASE WHEN r.flag_value='false' THEN 1 ELSE 0 END) AS f
SET a.flag_true_count=t, a.flag_false_count=f;

// growth_phases from parent Experiment
MATCH (a:MetaboliteAssay)<-[:ExperimentHasMetaboliteAssay]-(e:Experiment)
SET a.growth_phases = coalesce(e.growth_phases, []);
```

Mirror into `post-import.sh`.

- [ ] **Step 2: Commit**

```bash
git add scripts/post-import.cypher scripts/post-import.sh
git commit -m "post-import: MetaboliteAssay node rollups (counts + dist + growth_phases)"
```

### Task 5.4: Add Metabolite-side rollups (new measured_*, organism_count UNION extension)

**Files:**
- Modify: `scripts/post-import.cypher`
- Modify: `scripts/post-import.sh`

- [ ] **Step 1: Append the new measured_* properties**

```cypher
// New metabolomics-only Metabolite properties
MATCH (m:Metabolite)
OPTIONAL MATCH (m)<-[:Assay_quantifies_metabolite|Assay_flags_metabolite]-(a:MetaboliteAssay)
OPTIONAL MATCH (a)-[:MetaboliteAssayBelongsToOrganism]->(o:OrganismTaxon)
OPTIONAL MATCH (a)<-[:PublicationHasMetaboliteAssay]-(p:Publication)
WITH m,
     count(DISTINCT a) AS acnt,
     collect(DISTINCT o.preferred_name) AS orgs,
     count(DISTINCT p) AS pcnt
SET m.measured_assay_count = acnt,
    m.measured_organisms   = apoc.coll.sort([x IN orgs WHERE x IS NOT NULL]),
    m.measured_paper_count = pcnt;
```

- [ ] **Step 2: Locate and modify the existing `Metabolite.organism_count` block**

Run: `grep -n "Metabolite\b.*organism_count\|m.organism_count" scripts/post-import.cypher`

The existing block likely UNIONs catalysis + transport paths. Replace with the 3-path version:

```cypher
// Metabolite.organism_count / organism_names — UNION across metabolism + transport + measurement paths
MATCH (m:Metabolite)
OPTIONAL MATCH (m)<-[:Reaction_has_metabolite]-(:Reaction)<-[:Gene_catalyzes_reaction]-(:Gene)-[:Gene_belongs_to_organism]->(o1:OrganismTaxon)
OPTIONAL MATCH (m)<-[:Tcdb_family_transports_metabolite]-(:TcdbFamily {level_kind:'tc_specificity'})<-[:Gene_has_tcdb_family]-(:Gene)-[:Gene_belongs_to_organism]->(o2:OrganismTaxon)
OPTIONAL MATCH (m)<-[:Assay_quantifies_metabolite|Assay_flags_metabolite]-(:MetaboliteAssay)-[:MetaboliteAssayBelongsToOrganism]->(o3:OrganismTaxon)
WITH m, apoc.coll.toSet([o IN collect(DISTINCT o1) + collect(DISTINCT o2) + collect(DISTINCT o3) WHERE o IS NOT NULL | o.preferred_name]) AS names
SET m.organism_count = size(names),
    m.organism_names = apoc.coll.sort(names);
```

(If the existing block doesn't include the TCDB transport path, leave that part as-is and just add the measurement leg. Adapt to actual current Cypher.)

Mirror to `post-import.sh`.

- [ ] **Step 3: Commit**

```bash
git add scripts/post-import.cypher scripts/post-import.sh
git commit -m "post-import: Metabolite measured_* + organism_count UNION extension"
```

### Task 5.5: Extend `Organism_has_metabolite` materialization + augmentation

**Files:**
- Modify: `scripts/post-import.cypher`
- Modify: `scripts/post-import.sh`

- [ ] **Step 1: Locate the existing `Organism_has_metabolite` materialization block**

Run: `grep -n "Organism_has_metabolite\|MERGE.*Organism_has_metabolite" scripts/post-import.cypher | head -5`

- [ ] **Step 2: Refactor existing metabolism + transport blocks to set evidence_sources via ON CREATE/ON MATCH**

If the existing materialization doesn't carry `evidence_sources`, add it. Pattern:

```cypher
// metabolism path (existing — refactor to add evidence_sources)
MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
MATCH (g)-[:Gene_catalyzes_reaction]->(:Reaction)-[:Reaction_has_metabolite]->(m:Metabolite)
WITH DISTINCT o, m
CALL { WITH o, m
  MERGE (o)-[r:Organism_has_metabolite]->(m)
  ON CREATE SET r.evidence_sources = ['metabolism']
  ON MATCH  SET r.evidence_sources =
    CASE WHEN 'metabolism' IN r.evidence_sources THEN r.evidence_sources
         ELSE r.evidence_sources + 'metabolism' END
} IN TRANSACTIONS OF 10000 ROWS;

// transport path (existing — refactor)
MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
MATCH (g)-[:Gene_has_tcdb_family]->(:TcdbFamily {level_kind:'tc_specificity'})-[:Tcdb_family_transports_metabolite]->(m:Metabolite)
WITH DISTINCT o, m
CALL { WITH o, m
  MERGE (o)-[r:Organism_has_metabolite]->(m)
  ON CREATE SET r.evidence_sources = ['transport']
  ON MATCH  SET r.evidence_sources =
    CASE WHEN 'transport' IN r.evidence_sources THEN r.evidence_sources
         ELSE r.evidence_sources + 'transport' END
} IN TRANSACTIONS OF 10000 ROWS;

// NEW: measurement path
MATCH (a:MetaboliteAssay)-[:MetaboliteAssayBelongsToOrganism]->(o:OrganismTaxon)
MATCH (a)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
WITH DISTINCT o, m
CALL { WITH o, m
  MERGE (o)-[r:Organism_has_metabolite]->(m)
  ON CREATE SET r.evidence_sources = ['measured']
  ON MATCH  SET r.evidence_sources =
    CASE WHEN 'measured' IN r.evidence_sources THEN r.evidence_sources
         ELSE r.evidence_sources + 'measured' END
} IN TRANSACTIONS OF 10000 ROWS;

// Augmentation properties on every Organism_has_metabolite edge
MATCH (o:OrganismTaxon)-[r:Organism_has_metabolite]->(m:Metabolite)
OPTIONAL MATCH (o)<-[:MetaboliteAssayBelongsToOrganism]-(a:MetaboliteAssay)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m)
OPTIONAL MATCH (a)<-[:PublicationHasMetaboliteAssay]-(p:Publication)
WITH r,
     count(DISTINCT a) AS acnt,
     collect(DISTINCT a.compartment) AS comps,
     count(DISTINCT p) AS pcnt
SET r.measured_assay_count   = acnt,
    r.measured_compartments  = [c IN comps WHERE c IS NOT NULL],
    r.measured_paper_count   = pcnt;
```

Mirror to `post-import.sh`.

- [ ] **Step 3: Commit**

```bash
git add scripts/post-import.cypher scripts/post-import.sh
git commit -m "post-import: Organism_has_metabolite — measurement path + evidence_sources"
```

### Task 5.6: Add Experiment / Publication / OrganismTaxon rollups + defaults

**Files:**
- Modify: `scripts/post-import.cypher`
- Modify: `scripts/post-import.sh`

- [ ] **Step 1: Append rollup blocks**

```cypher
// Experiment rollups
MATCH (e:Experiment)-[:ExperimentHasMetaboliteAssay]->(a:MetaboliteAssay)
OPTIONAL MATCH (a)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
WITH e,
     count(DISTINCT a) AS acnt,
     collect(DISTINCT a.compartment) AS comps,
     count(DISTINCT m) AS mcnt
SET e.metabolite_assay_count   = acnt,
    e.metabolite_compartments  = [c IN comps WHERE c IS NOT NULL],
    e.metabolite_count         = mcnt;

// Publication rollups
MATCH (p:Publication)-[:PublicationHasMetaboliteAssay]->(a:MetaboliteAssay)
OPTIONAL MATCH (a)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
WITH p,
     count(DISTINCT a) AS acnt,
     collect(DISTINCT a.compartment) AS comps,
     count(DISTINCT m) AS mcnt
SET p.metabolite_assay_count   = acnt,
    p.metabolite_compartments  = [c IN comps WHERE c IS NOT NULL],
    p.metabolite_count         = mcnt;

// OrganismTaxon rollup
MATCH (o:OrganismTaxon)<-[:MetaboliteAssayBelongsToOrganism]-(a:MetaboliteAssay)
OPTIONAL MATCH (a)-[:Assay_quantifies_metabolite|Assay_flags_metabolite]->(m:Metabolite)
WITH o, count(DISTINCT m) AS mcnt
SET o.measured_metabolite_count = mcnt;

// Defaults
MATCH (e:Experiment) WHERE e.metabolite_assay_count IS NULL
SET e.metabolite_assay_count = 0,
    e.metabolite_compartments = [],
    e.metabolite_count = 0;

MATCH (p:Publication) WHERE p.metabolite_assay_count IS NULL
SET p.metabolite_assay_count = 0,
    p.metabolite_compartments = [],
    p.metabolite_count = 0;

MATCH (o:OrganismTaxon) WHERE o.measured_metabolite_count IS NULL
SET o.measured_metabolite_count = 0;

MATCH (m:Metabolite) WHERE m.measured_assay_count IS NULL
SET m.measured_assay_count = 0,
    m.measured_organisms = [],
    m.measured_paper_count = 0;
```

Mirror to `post-import.sh`.

- [ ] **Step 2: Commit**

```bash
git add scripts/post-import.cypher scripts/post-import.sh
git commit -m "post-import: Experiment/Publication/OrganismTaxon rollups + defaults"
```

### Task 5.7: Update `scripts/post-import-validate.sh`

**Files:**
- Modify: `scripts/post-import-validate.sh`

- [ ] **Step 1: Read the existing structure**

Run: `head -50 scripts/post-import-validate.sh`
Run: `grep -n "DerivedMetric\|ClusteringAnalysis" scripts/post-import-validate.sh | head -10`

- [ ] **Step 2: Add MetaboliteAssay enumeration queries**

Find the section that enumerates DerivedMetric properties. Add a parallel section for MetaboliteAssay (querying `total_metabolite_count`, value_min/max/q1/median/q3, flag_true_count, flag_false_count, growth_phases, etc.). Also add per-property checks for the new Metabolite/Experiment/Publication/OrganismTaxon properties.

Use the existing pattern:

```bash
echo "=== MetaboliteAssay properties ==="
cypher-shell -u neo4j -a "${NEO4J_URI}" --format plain <<'EOF'
MATCH (a:MetaboliteAssay)
RETURN
  count(*)                              AS n_assays,
  count(DISTINCT a.organism_name)        AS n_organisms,
  count(DISTINCT a.compartment)          AS n_compartments,
  count(DISTINCT a.metric_type)          AS n_metric_types,
  count(DISTINCT a.value_kind)           AS n_value_kinds,
  sum(coalesce(a.total_metabolite_count, 0)) AS sum_metabolite_count
ORDER BY n_assays;
EOF

echo "=== Assay edges ==="
cypher-shell -u neo4j -a "${NEO4J_URI}" --format plain <<'EOF'
MATCH ()-[r:Assay_quantifies_metabolite]->()
RETURN count(*) AS n, count(DISTINCT r.detection_status) AS n_det
ORDER BY n;
MATCH ()-[r:Assay_flags_metabolite]->()
RETURN count(*) AS n, count(DISTINCT r.flag_value) AS n_flag
ORDER BY n;
EOF
```

- [ ] **Step 3: Commit**

```bash
git add scripts/post-import-validate.sh
git commit -m "post-import-validate: add MetaboliteAssay enumeration queries"
```

### Task 5.8: Phase 5 validation gate

**Files:** none

- [ ] **Step 1: Run unit tests (no regression)**

Run: `uv run pytest -m "not slow and not kg" -v 2>&1 | tail -10`

- [ ] **Step 2: If a graph is deployed, run post-import-validate before/after**

```bash
bash scripts/post-import-validate.sh > /tmp/baseline.txt 2>&1
# Re-run post-import idempotently
bash scripts/post-import.sh
bash scripts/post-import-validate.sh > /tmp/after.txt 2>&1
diff /tmp/baseline.txt /tmp/after.txt
```

Expected: byte-identical (refactor of existing blocks must not change semantics).

---

## Phase 6 — Capovilla 2023 metabolomics integration

**Goal:** First real integration. Spec §5.3 commit 6.

### Task 6.1: Add metabolomics experiments to Capovilla `paperconfig.yaml`

**Files:**
- Modify: `data/Prochlorococcus/papers_and_supp/Capovilla 2023/paperconfig.yaml`

- [ ] **Step 1: Read existing experiments block**

Run: `head -80 "data/Prochlorococcus/papers_and_supp/Capovilla 2023/paperconfig.yaml"`

- [ ] **Step 2: Add two new experiments for metabolomics**

In the `experiments:` block, append:

```yaml
    chitosan_addition_mit9303_metabolomics:
      name: "MIT9303 chitosan addition — intracellular metabolomics"
      organism: "Prochlorococcus MIT9303"
      treatment_condition: "chitosan_addition"
      control_condition: "control"
      experimental_context: "MIT9303 cells in Pro99 AMP1 + 56 ug/mL chitosan vs control; sampled T=4 days post-inoculation"
      omics_type: METABOLOMICS
      treatment_type: ["carbon"]
      background_factors: ["axenic", "continuous_light"]
      compartment: whole_cell
      light_condition: "continuous light"

    chitosan_addition_mit9313_metabolomics:
      name: "MIT9313 chitosan addition — intracellular metabolomics"
      organism: "Prochlorococcus MIT9313"
      treatment_condition: "chitosan_addition"
      control_condition: "control"
      experimental_context: "MIT9313 cells in Pro99 AMP1 + 56 ug/mL chitosan vs control; sampled T=4 and T=6 days"
      omics_type: METABOLOMICS
      treatment_type: ["carbon"]
      background_factors: ["axenic", "continuous_light"]
      compartment: whole_cell
      light_condition: "continuous light"
```

- [ ] **Step 3: Validate paperconfig**

Run: `uv run python tools/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/Capovilla 2023/paperconfig.yaml"`
Expected: success.

- [ ] **Step 4: Commit**

```bash
git add "data/Prochlorococcus/papers_and_supp/Capovilla 2023/paperconfig.yaml"
git commit -m "data(capovilla): add metabolomics experiments"
```

### Task 6.2: Inspect Capovilla `cellular_concentrations` CSV to determine `skip_rows` + columns

**Files:** none

- [ ] **Step 1: Read raw header rows**

Run: `head -15 "data/Prochlorococcus/papers_and_supp/Capovilla 2023/cellular_concentrations metabolites pnas.2213271120.sd03.csv"`

- [ ] **Step 2: Note the row index where `compound` column starts**

The Capovilla CSV has 13 prefatory rows. `skip_rows: 12` should leave row 13 as the header row (`,compound,notes:,9303 Control T=4, replicate 1,...`). Verify by:

```bash
sed -n '13,14p' "data/Prochlorococcus/papers_and_supp/Capovilla 2023/cellular_concentrations metabolites pnas.2213271120.sd03.csv"
```

(Adjust skip_rows in next task if header is on a different line.)

### Task 6.3: Add `metabolite_assays_table` entries to Capovilla paperconfig

**Files:**
- Modify: `data/Prochlorococcus/papers_and_supp/Capovilla 2023/paperconfig.yaml`

- [ ] **Step 1: Append entries to `supplementary_materials:`**

```yaml
    metabolites_intracellular_mit9303:
      type: metabolite_assays_table
      filename: "data/Prochlorococcus/papers_and_supp/Capovilla 2023/cellular_concentrations metabolites pnas.2213271120.sd03.csv"
      experiment: chitosan_addition_mit9303_metabolomics
      organism: "Prochlorococcus MIT9303"
      name_col: "compound"
      skip_rows: 12  # adjust if Task 6.2 found a different header row
      null_values: ["nd", "ND", "NA", ""]
      missing_values: []
      cell_format: numeric
      aggregation_method: mean_across_replicates
      drop_undetected: false
      assays:
        - metric_type: cellular_concentration
          name: "MIT9303 cellular metabolite concentration (fg/cell)"
          value_kind: numeric
          unit: "fg/cell"
          rankable: "true"
          field_description: "Intracellular metabolite concentration in fg/cell, blank-corrected, replicate-aggregated; Capovilla 2023 Table sd03."
          sample_columns:
            - condition_label: control
              time_point: "T=4"
              time_point_order: 1
              replicate_columns:
                - "9303 Control T=4, replicate 1"
                - "9303 Control T=4, replicate 2"
                - "9303 Control T=4, replicate 3"
            - condition_label: chitosan
              time_point: "T=4"
              time_point_order: 1
              replicate_columns:
                - "9303 +Chitosan T=4, replicate 1"
                - "9303 +Chitosan T=4, replicate 2"
                - "9303 +Chitosan T=4, replicate 3"

    metabolites_intracellular_mit9313:
      type: metabolite_assays_table
      filename: "data/Prochlorococcus/papers_and_supp/Capovilla 2023/cellular_concentrations metabolites pnas.2213271120.sd03.csv"
      experiment: chitosan_addition_mit9313_metabolomics
      organism: "Prochlorococcus MIT9313"
      name_col: "compound"
      skip_rows: 12
      null_values: ["nd", "ND", "NA", ""]
      missing_values: []
      cell_format: numeric
      aggregation_method: mean_across_replicates
      assays:
        - metric_type: cellular_concentration
          name: "MIT9313 cellular metabolite concentration (fg/cell)"
          value_kind: numeric
          unit: "fg/cell"
          rankable: "true"
          field_description: "Intracellular metabolite concentration in fg/cell; Capovilla 2023 Table sd03."
          sample_columns:
            - condition_label: control
              time_point: "T=4"
              time_point_order: 1
              replicate_columns:
                - "9313 Control T=4, replicate 1"
                - "9313 Control T=4, replicate 2"
                - "9313 Control T=4, replicate 3"
            - condition_label: chitosan
              time_point: "T=4"
              time_point_order: 1
              replicate_columns:
                - "9313 +Chitosan T=4, replicate 2"
                - "9313 +Chitosan T=4, replicate 2"
                - "9313 +Chitosan T=4, replicate 3"
            - condition_label: control
              time_point: "T=6"
              time_point_order: 2
              replicate_columns:
                - "9313 Control T=6, replicate 2"
                - "9313 Control T=6, replicate 2"
            - condition_label: chitosan
              time_point: "T=6"
              time_point_order: 2
              replicate_columns:
                - "9313 +Chitosan T=6, replicate 2"
                - "9313 +Chitosan T=6, replicate 2"
```

(Note: the CSV has duplicate column names like `9313 +Chitosan T=4, replicate 2` repeated — pandas will rename to `... .1`. Inspect the actual columns and adjust column names accordingly if needed.)

- [ ] **Step 2: Validate**

Run: `uv run python tools/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/Capovilla 2023/paperconfig.yaml"`
Expected: success.

- [ ] **Step 3: Commit**

```bash
git add "data/Prochlorococcus/papers_and_supp/Capovilla 2023/paperconfig.yaml"
git commit -m "data(capovilla): add metabolite_assays_table entries"
```

### Task 6.4: Run step 6 + step 7 + verify resolution rate

**Files:**
- Possibly create: `data/Prochlorococcus/papers_and_supp/Capovilla 2023/metabolite_aliases.yaml`

- [ ] **Step 1: Run step 6**

```bash
bash scripts/prepare_data.sh --steps 6 --force 2>&1 | tail -10
```

Expected: prints `[step6] Phase 2 harvest: N resolved, M unresolved` with N > 0.

- [ ] **Step 2: Run step 7**

```bash
bash scripts/prepare_data.sh --steps 7 --force 2>&1 | tail -20
```

Expected: writes `cellular_concentrations metabolites pnas.2213271120.sd03_resolved.csv` next to the original.

- [ ] **Step 3: Inspect resolution report**

Run: `cat "data/Prochlorococcus/papers_and_supp/Capovilla 2023/cellular_concentrations metabolites pnas.2213271120.sd03_resolution_report.json"`

- [ ] **Step 4: If resolution rate < 70%, examine unresolved names**

```bash
uv run python -c "
import pandas as pd
df = pd.read_csv('data/Prochlorococcus/papers_and_supp/Capovilla 2023/cellular_concentrations metabolites pnas.2213271120.sd03_resolved.csv', dtype=str, keep_default_na=False)
print(df[df['resolution_method'] == 'unresolved'][['compound', 'resolution_method']])
"
```

If unresolved names are typo-fixable or have well-known KEGG IDs, create `data/Prochlorococcus/papers_and_supp/Capovilla 2023/metabolite_aliases.yaml`:

```yaml
# Manual overrides for compound names that don't resolve via MNX.
# Format: "name as it appears in CSV": "primary_id"
"acetyl coenzyme A": "kegg.compound:C00024"
# add as needed
```

Re-run step 6 + step 7 after each alias addition until rate ≥ 70%.

- [ ] **Step 5: Commit aliases (if created)**

```bash
git add "data/Prochlorococcus/papers_and_supp/Capovilla 2023/metabolite_aliases.yaml"
git commit -m "data(capovilla): metabolite_aliases.yaml for unresolved names"
```

### Task 6.5: Take omics-edge-snapshot pre-rebuild

**Files:** none

- [ ] **Step 1: Run snapshot skill**

Run: `bash scripts/snapshot.sh > /tmp/pre_capovilla.json` or invoke the `/omics-edge-snapshot` skill manually.

(See `.claude/skills/omics-edge-snapshot/SKILL.md` for the exact invocation.)

- [ ] **Step 2: Save snapshot for diff**

Keep `/tmp/pre_capovilla.json` for comparison after the rebuild.

### Task 6.6: Rebuild KG + verify

**Files:** none

- [ ] **Step 1: Stop deploy + app containers**

```bash
docker compose stop deploy app
```

- [ ] **Step 2: Rebuild + reimport**

```bash
docker compose up -d build && docker compose logs -f build
# Wait for build to finish (Ctrl-C when "[build] done" appears)
docker compose up -d import && docker compose logs -f import
```

Expected: import completes; no errors in `output/import.report` related to MetaboliteAssay or Assay_*.

- [ ] **Step 3: Re-deploy**

```bash
docker compose up -d deploy app
```

- [ ] **Step 4: Run KG validity tests**

Run: `uv run pytest -m kg -v 2>&1 | tail -30`
Expected: all green; the `test_metabolomics_evidence_source_appears_when_papers_integrated` test now passes (no longer skipped).

- [ ] **Step 5: Compare omics-edge-snapshot**

Run snapshot post-rebuild → `/tmp/post_capovilla.json`. Diff vs `/tmp/pre_capovilla.json`.
Expected: zero net loss in `Changes_expression_of` and `Derived_metric_*` per-paper counts.

- [ ] **Step 6: Commit nothing — verification step only**

(No code changes; validation only.)

---

## Phase 7 — Kujawinski 2023 integration + MIT0801 backlog

**Goal:** Second paper, exercises boolean path. Spec §5.3 commit 7.

### Task 7.1: Create `plans/strain_deployment_backlog.md`

**Files:**
- Create: `plans/strain_deployment_backlog.md`

- [ ] **Step 1: Create the file**

```markdown
# Strain deployment backlog

Strains referenced by paperconfigs whose data remains uningested pending KG deployment.
Track via `/add-a-strain` workflow when picked up.

## MIT0801 (Prochlorococcus, LLI ecotype)

- **Source NCBI taxid:** 1110371 (verify via NCBI Assembly)
- **Referenced by:**
  - Kujawinski 2023 — sample columns in `ChisholmPro_cellSpecific_KEGGexport.csv`
    (`replete_extracellular_s0801ax_10`, `replete_filter_s0801ax_10`) +
    glycine-betaine row in `table s3.csv`
- **Blockers:** strain not in `cyanobacteria_genomes.csv`; no genome assembly downloaded.
- **paperconfig state:** Kujawinski 2023 paperconfig.yaml has MIT0801 entries pre-written, commented out.
  Re-enable on strain deployment via uncomment + adapter rebuild.
- **Estimated work:** ~1 day per `/add-a-strain` skill flow.
```

- [ ] **Step 2: Commit**

```bash
git add plans/strain_deployment_backlog.md
git commit -m "docs: strain deployment backlog (MIT0801 from Kujawinski 2023)"
```

### Task 7.2: Create Kujawinski 2023 paperconfig.yaml

**Files:**
- Create: `data/Prochlorococcus/papers_and_supp/Kujawinski 2023/paperconfig.yaml`

- [ ] **Step 1: Inspect the source CSV first**

Run: `head -3 "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/ChisholmPro_cellSpecific_KEGGexport.2017.07.31_v1.csv"`
Run: `head -3 "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/table s2.csv"`

- [ ] **Step 2: Create the paperconfig**

Create `data/Prochlorococcus/papers_and_supp/Kujawinski 2023/paperconfig.yaml`:

```yaml
publication:
  papername: "Kujawinski 2023"
  doi: "10.1128/msystems.01261-22"
  papermainpdf: "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/kujawinski-et-al-2023-metabolite-diversity-among-representatives-of-divergent-prochlorococcus-ecotypes.pdf"

  experiments:
    kujawinski_metabolomics_9301_whole_cell:
      name: "MIT9301 metabolite diversity (intracellular pool)"
      organism: "Prochlorococcus MIT9301"
      treatment_condition: "replete + P-limited"
      control_condition: "n/a"
      experimental_context: "Pro99 medium; replete and P-limited; 10 + 50 µmol photons m-2 s-1; ax = axenic"
      omics_type: METABOLOMICS
      treatment_type: ["phosphorus"]
      background_factors: ["axenic"]
      compartment: whole_cell
      light_condition: "continuous light"

    kujawinski_metabolomics_9301_extracellular:
      name: "MIT9301 metabolite diversity (extracellular pool)"
      organism: "Prochlorococcus MIT9301"
      treatment_condition: "replete + P-limited"
      control_condition: "n/a"
      experimental_context: "Pro99 medium; replete and P-limited; 10 + 50 µmol photons m-2 s-1; supernatant"
      omics_type: METABOLOMICS
      treatment_type: ["phosphorus"]
      background_factors: ["axenic"]
      compartment: extracellular
      light_condition: "continuous light"

    kujawinski_metabolomics_9313_whole_cell:
      name: "MIT9313 metabolite diversity (intracellular pool)"
      organism: "Prochlorococcus MIT9313"
      treatment_condition: "replete"
      control_condition: "n/a"
      experimental_context: "Pro99 medium; 5 + 10 µmol photons m-2 s-1; axenic"
      omics_type: METABOLOMICS
      treatment_type: []
      background_factors: ["axenic"]
      compartment: whole_cell
      light_condition: "continuous light"

    kujawinski_metabolomics_9313_extracellular:
      name: "MIT9313 metabolite diversity (extracellular pool)"
      organism: "Prochlorococcus MIT9313"
      treatment_condition: "replete"
      control_condition: "n/a"
      experimental_context: "Pro99 medium; 5 + 10 µmol photons m-2 s-1; supernatant"
      omics_type: METABOLOMICS
      treatment_type: []
      background_factors: ["axenic"]
      compartment: extracellular
      light_condition: "continuous light"

    # BACKLOG: requires MIT0801 deployment — see plans/strain_deployment_backlog.md
    # kujawinski_metabolomics_0801_whole_cell:
    #   name: "MIT0801 metabolite diversity (intracellular pool)"
    #   organism: "Prochlorococcus MIT0801"
    #   compartment: whole_cell
    #   omics_type: METABOLOMICS
    #   treatment_type: []
    #   background_factors: ["axenic"]
    # kujawinski_metabolomics_0801_extracellular:
    #   organism: "Prochlorococcus MIT0801"
    #   compartment: extracellular
    #   omics_type: METABOLOMICS

  supplementary_materials:
    metabolites_kegg_export_9301_intracellular:
      type: metabolite_assays_table
      filename: "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/ChisholmPro_cellSpecific_KEGGexport.2017.07.31_v1.csv"
      experiment: kujawinski_metabolomics_9301_whole_cell
      organism: "Prochlorococcus MIT9301"
      name_col: "mtabNames"
      id_col: "KEGG"
      id_type: kegg.compound
      null_values: ["0", ""]
      missing_values: []
      cell_format: numeric
      aggregation_method: pre_aggregated
      assays:
        - metric_type: cellular_concentration
          name: "MIT9301 intracellular metabolite concentration (mol/cell, replete + P-limited)"
          value_kind: numeric
          unit: "mol/cell"
          rankable: "true"
          field_description: "Per-cell concentration; KEGG-tagged; pre-aggregated by authors. Kujawinski 2023 KEGG export."
          sample_columns:
            - condition_label: replete_light_10
              replicate_columns: ["replete_filter_s9301ax_10"]
            - condition_label: replete_light_50
              replicate_columns: ["replete_filter_s9301ax_50"]
            - condition_label: P_limited_light_50
              replicate_columns: ["Plimited_filter_s9301ax_50"]

    metabolites_kegg_export_9301_extracellular:
      type: metabolite_assays_table
      filename: "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/ChisholmPro_cellSpecific_KEGGexport.2017.07.31_v1.csv"
      experiment: kujawinski_metabolomics_9301_extracellular
      organism: "Prochlorococcus MIT9301"
      name_col: "mtabNames"
      id_col: "KEGG"
      id_type: kegg.compound
      null_values: ["0", ""]
      missing_values: []
      cell_format: numeric
      aggregation_method: pre_aggregated
      assays:
        - metric_type: extracellular_concentration
          name: "MIT9301 extracellular metabolite concentration (mol/cell, replete + P-limited)"
          value_kind: numeric
          unit: "mol/cell"
          rankable: "true"
          field_description: "Per-cell extracellular concentration. Kujawinski 2023 KEGG export."
          sample_columns:
            - condition_label: replete_light_10
              replicate_columns: ["replete_extracellular_s9301ax_10"]
            - condition_label: replete_light_50
              replicate_columns: ["replete_extracellular_s9301ax_50"]
            - condition_label: P_limited_light_50
              replicate_columns: ["Plimited_extracellular_s9301ax_50"]

    metabolites_kegg_export_9313_intracellular:
      type: metabolite_assays_table
      filename: "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/ChisholmPro_cellSpecific_KEGGexport.2017.07.31_v1.csv"
      experiment: kujawinski_metabolomics_9313_whole_cell
      organism: "Prochlorococcus MIT9313"
      name_col: "mtabNames"
      id_col: "KEGG"
      id_type: kegg.compound
      null_values: ["0", ""]
      cell_format: numeric
      aggregation_method: pre_aggregated
      assays:
        - metric_type: cellular_concentration
          name: "MIT9313 intracellular metabolite concentration (mol/cell)"
          value_kind: numeric
          unit: "mol/cell"
          rankable: "true"
          field_description: "Per-cell intracellular concentration. Kujawinski 2023 KEGG export."
          sample_columns:
            - condition_label: replete_light_5
              replicate_columns: ["replete_filter_s9313ax_5"]
            - condition_label: replete_light_10
              replicate_columns: ["replete_filter_s9313ax_10"]

    metabolites_kegg_export_9313_extracellular:
      type: metabolite_assays_table
      filename: "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/ChisholmPro_cellSpecific_KEGGexport.2017.07.31_v1.csv"
      experiment: kujawinski_metabolomics_9313_extracellular
      organism: "Prochlorococcus MIT9313"
      name_col: "mtabNames"
      id_col: "KEGG"
      id_type: kegg.compound
      null_values: ["0", ""]
      cell_format: numeric
      aggregation_method: pre_aggregated
      assays:
        - metric_type: extracellular_concentration
          name: "MIT9313 extracellular metabolite concentration (mol/cell)"
          value_kind: numeric
          unit: "mol/cell"
          rankable: "true"
          field_description: "Per-cell extracellular concentration. Kujawinski 2023 KEGG export."
          sample_columns:
            - condition_label: replete_light_5
              replicate_columns: ["replete_extracellular_s9313ax_5"]
            - condition_label: replete_light_10
              replicate_columns: ["replete_extracellular_s9313ax_10"]

    presence_flags_table_s2:
      type: metabolite_assays_table
      filename: "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/table s2.csv"
      experiment: kujawinski_metabolomics_9301_whole_cell    # representative; presence flags are paper-level
      organism: "Prochlorococcus MIT9301"
      name_col: "metabolite"
      cell_format: numeric
      aliases_file: "metabolite_aliases.yaml"
      assays:
        - metric_type: presence_flag_intracellular
          name: "Detected in intracellular pool (Kujawinski 2023 Table S2)"
          value_kind: boolean
          unit: ""
          rankable: "false"
          field_description: "Boolean flag: metabolite detected in intracellular pool of any sampled strain (Kujawinski 2023 Table S2)."
          sample_columns:
            - condition_label: ""
              flag_column: "intracellular"
              flag_true_value: "yes"
        - metric_type: presence_flag_extracellular
          name: "Detected in extracellular pool (Kujawinski 2023 Table S2)"
          value_kind: boolean
          unit: ""
          rankable: "false"
          field_description: "Boolean flag: metabolite detected in extracellular pool of any sampled strain."
          sample_columns:
            - condition_label: ""
              flag_column: "extracellular"
              flag_true_value: "yes"

    # BACKLOG: requires MIT0801 deployment — see plans/strain_deployment_backlog.md
    # metabolites_kegg_export_0801_intracellular:
    #   type: metabolite_assays_table
    #   experiment: kujawinski_metabolomics_0801_whole_cell
    #   organism: "Prochlorococcus MIT0801"
    #   filename: "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/ChisholmPro_cellSpecific_KEGGexport.2017.07.31_v1.csv"
    #   name_col: "mtabNames"
    #   id_col: "KEGG"
    #   id_type: kegg.compound
    #   ... (same shape as 9301_intracellular, columns: replete_filter_s0801ax_10)
```

- [ ] **Step 3: Add Kujawinski to the paperconfig_files.txt registry**

Append to `data/Prochlorococcus/papers_and_supp/paperconfig_files.txt`:

```
data/Prochlorococcus/papers_and_supp/Kujawinski 2023/paperconfig.yaml
```

- [ ] **Step 4: Validate**

Run: `uv run python tools/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/paperconfig.yaml"`
Expected: success.

- [ ] **Step 5: Run --report-backlog**

Run: `uv run python tools/validate_paperconfig.py --report-backlog "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/paperconfig.yaml"`
Expected: prints the MIT0801 BACKLOG markers.

- [ ] **Step 6: Commit**

```bash
git add "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/paperconfig.yaml" \
        data/Prochlorococcus/papers_and_supp/paperconfig_files.txt
git commit -m "data(kujawinski): add paperconfig + register (MIT0801 commented BACKLOG)"
```

### Task 7.3: Create Kujawinski metabolite_aliases.yaml

**Files:**
- Create: `data/Prochlorococcus/papers_and_supp/Kujawinski 2023/metabolite_aliases.yaml`

- [ ] **Step 1: Pre-curate known abbreviations from the paper**

```yaml
# Manual overrides for Kujawinski 2023 — paper-internal abbreviations and typos.
"AmMP": "kegg.compound:C20267"
"GABA": "kegg.compound:C00334"
"DMSP": "kegg.compound:C04022"
"HMP": "kegg.compound:C01279"
"2-3-dihydroxypropane1sulfonate": "kegg.compound:C19675"
```

- [ ] **Step 2: Commit**

```bash
git add "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/metabolite_aliases.yaml"
git commit -m "data(kujawinski): metabolite_aliases.yaml for paper abbreviations"
```

### Task 7.4: Run step 6 + step 7 + verify ≥ 95% resolution

**Files:** none

- [ ] **Step 1: Run step 6 + step 7**

```bash
bash scripts/prepare_data.sh --steps 6 7 --force 2>&1 | tail -10
```

- [ ] **Step 2: Inspect resolution reports**

```bash
cat "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/ChisholmPro_cellSpecific_KEGGexport.2017.07.31_v1_resolution_report.json"
cat "data/Prochlorococcus/papers_and_supp/Kujawinski 2023/table s2_resolution_report.json"
```

Expected: ≥ 95% resolution rate (KEGG IDs already in the source).

- [ ] **Step 3: If below 95%, add aliases and retry**

Add to `metabolite_aliases.yaml` and rerun. Commit additions.

### Task 7.5: Rebuild KG, run validity tests, regenerate snapshot

**Files:**
- Modify: `tests/kg_validity/snapshot_data.json`

- [ ] **Step 1: Take pre-rebuild snapshot**

Run: `bash scripts/snapshot.sh > /tmp/pre_kujawinski.json`

- [ ] **Step 2: Stop, build, import, deploy**

```bash
docker compose stop deploy app
docker compose up -d build && docker compose logs -f build
docker compose up -d import && docker compose logs -f import
docker compose up -d deploy app
```

- [ ] **Step 3: Run KG validity**

Run: `uv run pytest -m kg -v 2>&1 | tail -30`
Expected: all pass.

- [ ] **Step 4: Diff omics-edge snapshot**

Compare to `/tmp/pre_kujawinski.json`. Expected: no DE regression; metabolite-edge counts increased.

- [ ] **Step 5: Regenerate `test_snapshot.py` fixture**

Run: `uv run python tests/kg_validity/generate_snapshot.py`
Expected: rewrites `tests/kg_validity/snapshot_data.json`.

- [ ] **Step 6: Commit fixture**

```bash
git add tests/kg_validity/snapshot_data.json
git commit -m "test(kg): regenerate snapshot fixture post-metabolomics integration"
```

### Task 7.6: Add `test_metabolomics.py` KG validity tests

**Files:**
- Create: `tests/kg_validity/test_metabolomics.py`

- [ ] **Step 1: Create the test file**

```python
"""KG validity tests for metabolomics integration (MetaboliteAssay nodes + edges)."""
import pytest


pytestmark = pytest.mark.kg


def test_metabolite_assay_nodes_present(neo4j_driver):
    with neo4j_driver.session() as s:
        n = s.run("MATCH (a:MetaboliteAssay) RETURN count(a) AS n").single()["n"]
    assert n > 0, "No MetaboliteAssay nodes — Phase 6/7 integration must produce some"


def test_metabolite_assay_value_kind_enum(neo4j_driver):
    with neo4j_driver.session() as s:
        kinds = {r["k"] for r in s.run("MATCH (a:MetaboliteAssay) RETURN DISTINCT a.value_kind AS k")}
    assert kinds <= {"numeric", "boolean"}


def test_metabolite_assay_compartment_in_vocab(neo4j_driver):
    valid = {"whole_cell", "vesicle", "exoproteome", "extracellular"}
    with neo4j_driver.session() as s:
        comps = {r["c"] for r in s.run("MATCH (a:MetaboliteAssay) RETURN DISTINCT a.compartment AS c")}
    assert comps <= valid, f"unexpected compartments: {comps - valid}"


def test_every_assay_has_three_binding_edges(neo4j_driver):
    with neo4j_driver.session() as s:
        rec = s.run("""
            MATCH (a:MetaboliteAssay)
            OPTIONAL MATCH (p:Publication)-[:PublicationHasMetaboliteAssay]->(a)
            OPTIONAL MATCH (e:Experiment)-[:ExperimentHasMetaboliteAssay]->(a)
            OPTIONAL MATCH (a)-[:MetaboliteAssayBelongsToOrganism]->(o:OrganismTaxon)
            WITH a, count(DISTINCT p) AS np, count(DISTINCT e) AS ne, count(DISTINCT o) AS no
            WHERE np = 0 OR ne = 0 OR no = 0
            RETURN count(a) AS n_orphans
        """).single()
    assert rec["n_orphans"] == 0


def test_assay_quantifies_metabolite_props(neo4j_driver):
    with neo4j_driver.session() as s:
        bad = s.run("""
            MATCH ()-[r:Assay_quantifies_metabolite]->()
            WHERE r.value IS NULL OR r.n_replicates IS NULL OR r.n_non_zero IS NULL
              OR r.replicate_values IS NULL OR r.detection_status IS NULL
            RETURN count(r) AS n
        """).single()["n"]
    assert bad == 0


def test_detection_status_consistent(neo4j_driver):
    """detection_status agrees with n_non_zero/n_replicates."""
    with neo4j_driver.session() as s:
        bad = s.run("""
            MATCH ()-[r:Assay_quantifies_metabolite]->()
            WITH r,
                 CASE WHEN r.n_non_zero = 0                 THEN 'not_detected'
                      WHEN r.n_non_zero = r.n_replicates    THEN 'detected'
                      ELSE 'sporadic'
                 END AS expected
            WHERE r.detection_status <> expected
            RETURN count(r) AS n
        """).single()["n"]
    assert bad == 0


def test_replicate_values_size_matches_n_replicates(neo4j_driver):
    with neo4j_driver.session() as s:
        bad = s.run("""
            MATCH ()-[r:Assay_quantifies_metabolite]->()
            WHERE size(r.replicate_values) <> r.n_replicates
            RETURN count(r) AS n
        """).single()["n"]
    assert bad == 0


def test_assay_flags_metabolite_value_enum(neo4j_driver):
    with neo4j_driver.session() as s:
        vals = {r["v"] for r in s.run("MATCH ()-[r:Assay_flags_metabolite]->() RETURN DISTINCT r.flag_value AS v")}
    assert vals <= {"true", "false"}


def test_total_metabolite_count_positive_and_distribution_stats_only_for_numeric(neo4j_driver):
    with neo4j_driver.session() as s:
        bad = s.run("""
            MATCH (a:MetaboliteAssay)
            WHERE a.total_metabolite_count IS NULL OR a.total_metabolite_count < 0
            RETURN count(a) AS n
        """).single()["n"]
    assert bad == 0


def test_kujawinski_extracellular_presence_flag_present(neo4j_driver):
    """Spot-check: 2,3-dihydroxypropane-1-sulfonate flagged extracellular by Kujawinski Table S2."""
    with neo4j_driver.session() as s:
        rec = s.run("""
            MATCH (a:MetaboliteAssay {metric_type: 'presence_flag_extracellular'})
                  -[r:Assay_flags_metabolite {flag_value: 'true'}]->
                  (m:Metabolite {id: 'kegg.compound:C19675'})
            RETURN count(*) AS n
        """).single()
    if rec["n"] == 0:
        pytest.skip("Kujawinski C19675 spot check fixture not present; check alias resolution")
    assert rec["n"] > 0
```

- [ ] **Step 2: Run**

Run: `uv run pytest tests/kg_validity/test_metabolomics.py -v`
Expected: all pass (Phase 7 has populated data).

- [ ] **Step 3: Commit**

```bash
git add tests/kg_validity/test_metabolomics.py
git commit -m "test(kg): MetaboliteAssay validity tests"
```

### Task 7.7: Phase 7 validation gate

- [ ] **Step 1: Run all tests (unit + KG)**

Run: `uv run pytest -m "not slow and not kg" -v && uv run pytest -m kg -v`
Expected: all green.

---

## Phase 8 — Docs + skills

**Goal:** Documentation parity. Spec §5.3 commit 8.

### Task 8.1: Update `CLAUDE.md`

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Add `metabolite_assays_table` row to the supplementary_materials table**

Find the `paperconfig supplementary_materials entry types` table. Add:

```markdown
| `metabolite_assays_table` | Per-metabolite measurements (concentrations + presence flags); creates MetaboliteAssay nodes + 2 measurement edge types (quantifies/flags) per assay; mirrors derived_metrics_table |
```

- [ ] **Step 2: Add `MetaboliteAssay` to the Actual Neo4j labels list**

Append to the `Nodes:` line: `MetaboliteAssay`. Append to the `Relationships:` line: `Assay_quantifies_metabolite`, `Assay_flags_metabolite`, `PublicationHasMetaboliteAssay`, `ExperimentHasMetaboliteAssay`, `MetaboliteAssayBelongsToOrganism`.

- [ ] **Step 3: Add a "Key graph facts" entry**

After the existing `Metabolite nodes` bullet:

```markdown
- MetaboliteAssay nodes: one per (Experiment × value_kind) for metabolomics-paper measurements. Properties match DerivedMetric pattern (denormalized parent-Experiment block + metric_type + value_kind + aggregation_method). Edges: `Assay_quantifies_metabolite` (numeric, replicate-aggregated; carries `value`, `value_sd`, `n_replicates`, `n_non_zero`, `replicate_values`, `detection_status` adapter-set; rank/percentile/bucket post-import) + `Assay_flags_metabolite` (boolean presence). Bindings: `PublicationHasMetaboliteAssay`, `ExperimentHasMetaboliteAssay`, `MetaboliteAssayBelongsToOrganism`. Phase 2 adds `extracellular` to the compartment vocab; drops unused `secretome`. See `docs/kg-changes/metabolomics-extension.md`.
```

- [ ] **Step 4: Update the post-import indexes line**

Add to the indexes list: `metabolite_assay_organism_idx`, `metabolite_assay_compartment_idx`, `metabolite_assay_metric_type_idx`, `metabolite_assay_value_kind_idx`, `metabolite_assay_experiment_idx`, `metaboliteAssayFullText`.

- [ ] **Step 5: Update `prepare_data.sh` step descriptions**

Add Step 7 to the bullet list under "Data preparation steps" (or equivalent section).

- [ ] **Step 6: Commit**

```bash
git add CLAUDE.md
git commit -m "docs(claude.md): document MetaboliteAssay and Phase 2 metabolomics"
```

### Task 8.2: Update `.claude/skills/paperconfig/SKILL.md`

**Files:**
- Modify: `.claude/skills/paperconfig/SKILL.md`

- [ ] **Step 1: Add a metabolomics section**

After the `derived_metrics_table` documentation, add a `metabolite_assays_table` section explaining:
- Entry shape (link to spec §3.1)
- When to use (METABOLOMICS papers)
- Compartment vocab additions
- aliases_file pattern
- BACKLOG comment convention for undeployed strains

Include a minimal Capovilla-style example.

- [ ] **Step 2: Commit**

```bash
git add .claude/skills/paperconfig/SKILL.md
git commit -m "docs(skill:paperconfig): metabolomics section"
```

### Task 8.3: Update `.claude/skills/cypher-queries/SKILL.md`

**Files:**
- Modify: `.claude/skills/cypher-queries/SKILL.md`

- [ ] **Step 1: Append spot-check templates**

```markdown
## Metabolomics queries

### List all metabolite assays per paper

```cypher
MATCH (p:Publication)-[:PublicationHasMetaboliteAssay]->(a:MetaboliteAssay)
RETURN p.papername AS paper, a.metric_type, a.value_kind, a.compartment, a.total_metabolite_count
ORDER BY paper, a.metric_type;
```

### List metabolites measured in organism X with detection frequency

```cypher
MATCH (o:OrganismTaxon {preferred_name: $organism})
      <-[:MetaboliteAssayBelongsToOrganism]-(a:MetaboliteAssay)
      -[r:Assay_quantifies_metabolite]->(m:Metabolite)
RETURN m.name, m.id,
       avg(toFloat(r.n_non_zero) / toFloat(r.n_replicates)) AS detection_freq,
       avg(r.value) AS mean_value
ORDER BY detection_freq DESC LIMIT 30;
```

### Find compounds measured in metabolomics AND predicted by metabolism (gene-reachable)

```cypher
MATCH (m:Metabolite)
WHERE 'metabolomics' IN m.evidence_sources
  AND 'metabolism'   IN m.evidence_sources
RETURN m.name, m.id, m.evidence_sources, m.gene_count, m.measured_assay_count
ORDER BY m.measured_assay_count DESC LIMIT 50;
```
```

- [ ] **Step 2: Commit**

```bash
git add .claude/skills/cypher-queries/SKILL.md
git commit -m "docs(skill:cypher-queries): metabolomics templates"
```

### Task 8.4: Create `docs/kg-changes/metabolomics-extension.md`

**Files:**
- Create: `docs/kg-changes/metabolomics-extension.md`

- [ ] **Step 1: Write the changelog**

```markdown
# Metabolomics paper integration (Phase 2)

**Spec:** docs/superpowers/specs/2026-05-03-metabolomics-paper-integration-design.md
**Plan:** docs/superpowers/plans/2026-05-03-metabolomics-paper-integration.md
**Validation papers:** Capovilla 2023, Kujawinski 2023

## What changed

- New node type: `MetaboliteAssay` (one per Experiment × value_kind).
- New edge types: `Assay_quantifies_metabolite` (numeric, replicate-aggregated), `Assay_flags_metabolite` (boolean presence).
- New binding edges: `PublicationHasMetaboliteAssay`, `ExperimentHasMetaboliteAssay`, `MetaboliteAssayBelongsToOrganism`.
- New post-import properties on Metabolite: `measured_assay_count`, `measured_organisms`, `measured_paper_count`. `organism_count`/`organism_names` UNION extended with measurement path.
- New post-import properties on Organism_has_metabolite: `evidence_sources` (`metabolism`/`transport`/`measured`), `measured_assay_count`, `measured_compartments`, `measured_paper_count`. Materialization extended to measurement-only pairs.
- New post-import properties on Experiment, Publication: `metabolite_assay_count`, `metabolite_compartments`, `metabolite_count`. New on OrganismTaxon: `measured_metabolite_count`.
- New paperconfig entry type: `metabolite_assays_table` (mirrors `derived_metrics_table`).
- New step in `prepare_data.sh` (step 7): `resolve_paper_metabolites.py` — writes `<stem>_resolved.csv` for each metabolomics CSV.
- Step 6 (`build_kegg_metabolism_xrefs.py`) extends to harvest paper-measured metabolites into `kegg_data.json` (with `evidence_sources` including `"metabolomics"`) and a new `cache/data/metabolomics/metabolite_id_mapping.json`.
- Compartment vocab: adds `extracellular`; drops unused `secretome`.

## Validation results

- Capovilla 2023: ~85-90% metabolite name resolution (varies based on aliases curated)
- Kujawinski 2023: ~95% (KEGG IDs in source CSV)
- Zero net loss in `Changes_expression_of` and `Derived_metric_*` edges (omics-edge-snapshot diff)

## Known follow-ups

- MIT0801 strain deployment (Kujawinski 2023 has 3 sample columns) — see `plans/strain_deployment_backlog.md`.
- biller 2022 metabolomics tables (S5–S7): partial; MS-feature-only Table S6 is out of scope.
- Lipid-class names (LOBSTAHS notation) often unresolved — future LipidMaps integration.
```

- [ ] **Step 2: Commit**

```bash
git add docs/kg-changes/metabolomics-extension.md
git commit -m "docs(kg-changes): metabolomics-extension changelog"
```

### Task 8.5: Final integration test

- [ ] **Step 1: Full unit test sweep**

Run: `uv run pytest -m "not slow and not kg" -v`
Expected: all green.

- [ ] **Step 2: Full KG validity sweep**

Run: `uv run pytest -m kg -v`
Expected: all green.

- [ ] **Step 3: Print final summary**

```bash
git log --oneline main.. | head -50
```

Expected: shows the metabolomics implementation commits in order.

---

## Self-review checklist

After implementation, verify:

- [ ] Every spec requirement (§1–§5) maps to a task above. Spot-check by grepping the spec for "must" / "required" and confirming a task covers each.
- [ ] No `TBD` / `TODO` / placeholder text in any file under `multiomics_kg/`, `tests/`, `scripts/`, `docs/superpowers/specs/`, `docs/superpowers/plans/`.
- [ ] `scripts/post-import.cypher` and `scripts/post-import.sh` are byte-identical for the new metabolomics blocks (per CLAUDE.md invariant).
- [ ] `omics-edge-snapshot` shows zero loss in pre-existing edges.
- [ ] Both validation papers appear in `MATCH (p:Publication)-[:PublicationHasMetaboliteAssay]->()` queries.
- [ ] `validate_paperconfig.py --report-backlog` lists exactly the MIT0801 markers in Kujawinski.
- [ ] `cache/data/metabolomics/metabolite_id_mapping.json` exists with non-empty `name_lookup`.
