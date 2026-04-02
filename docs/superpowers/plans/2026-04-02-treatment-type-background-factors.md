# Treatment Type & Background Factors Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert `treatment_type` from single string to `str[]` on Experiment nodes, add `background_factors: str[]` to Experiment/ClusteringAnalysis/GeneCluster nodes, and update all 23 paperconfigs with proper values.

**Architecture:** Two independent categorical list fields sharing a canonical vocabulary. `treatment_type` captures what the DE comparison tests; `background_factors` captures conditions present but not compared. Both propagate to Publication and OrganismTaxon via post-import aggregation.

**Tech Stack:** Python adapters, YAML schema/paperconfigs, Cypher post-import scripts, pytest

**Spec:** `docs/superpowers/specs/2026-04-02-treatment-type-background-factors-design.md`

---

## File Map

| File | Action | Responsibility |
|---|---|---|
| `config/schema_config.yaml` | Modify | Add `background_factors: str[]` to Experiment, ClusteringAnalysis, GeneCluster, Publication, OrganismTaxon; change Experiment `treatment_type` to `str[]` |
| `.claude/skills/paperconfig/validate_paperconfig.py` | Modify | Extend canonical set, validate both fields as lists |
| `tests/test_paperconfig_validation.py` | Modify | Tests for list treatment_type, background_factors validation |
| `multiomics_kg/adapters/omics_adapter.py` | Modify | Read/normalize treatment_type as list, read background_factors |
| `multiomics_kg/adapters/cluster_adapter.py` | Modify | Read background_factors, pass to nodes |
| `scripts/post-import.sh` | Modify | Aggregate background_factors on Publication/OrganismTaxon, fix treatment_type aggregation for arrays |
| `scripts/post-import.cypher` | Modify | Same Cypher changes (kept in sync) |
| `tests/kg_validity/test_expression.py` | Modify | Update treatment_type assertion for array property |
| `tests/kg_validity/test_post_import.py` | Modify | Update coculture consistency check, add background_factors tests |
| All 23 paperconfig.yaml files | Modify | Pass 1: mechanical conversion. Pass 2: manual background_factors |
| `.claude/skills/paperconfig/SKILL.md` | Modify | Update templates and guidance |
| `CLAUDE.md` | Modify | Update schema docs |

---

### Task 1: Schema — Add background_factors, change Experiment treatment_type to array

**Files:**
- Modify: `config/schema_config.yaml:39` (Experiment treatment_type), `:26` (Publication), `:76` (GeneCluster), `:100` (ClusteringAnalysis), `:336` (OrganismTaxon)

- [ ] **Step 1: Change Experiment treatment_type from str to str[]**

In `config/schema_config.yaml`, line 39, change:
```yaml
    treatment_type: str
```
to:
```yaml
    treatment_type: str[]
```

- [ ] **Step 2: Add background_factors to Experiment node**

After `treatment_type: str[]` (line 39), add:
```yaml
    background_factors: str[]
```

- [ ] **Step 3: Add background_factors to GeneCluster node**

After line 76 (`treatment_type: str[]`), add:
```yaml
    background_factors: str[]
```

- [ ] **Step 4: Add background_factors to ClusteringAnalysis node**

After line 100 (`treatment_type: str[]`), add:
```yaml
    background_factors: str[]
```

- [ ] **Step 5: Add background_factors to Publication node**

After line 26 (`treatment_types: str[]`), add:
```yaml
    background_factors: str[]            # post-import: distinct background_factors from experiments
```

- [ ] **Step 6: Add background_factors to OrganismTaxon node**

After line 336 (`treatment_types: str[]`), add:
```yaml
    background_factors: str[]        # pre-computed: distinct background_factors across matched publications
```

- [ ] **Step 7: Commit**

```bash
git add config/schema_config.yaml
git commit -m "schema: add background_factors str[], change Experiment treatment_type to str[]"
```

---

### Task 2: Validator — Extend canonical vocabulary, validate both fields as lists

**Files:**
- Modify: `.claude/skills/paperconfig/validate_paperconfig.py:115-132` (canonical set), `:304-314` (treatment_type validation), `:363-388` (coculture consistency)
- Modify: `tests/test_paperconfig_validation.py:186-212` (treatment_type tests)

- [ ] **Step 1: Write failing test for list treatment_type acceptance**

In `tests/test_paperconfig_validation.py`, add after the `TestCanonicalTreatmentType` class (after line 212):

```python
class TestTreatmentTypeList:
    """Validator accepts treatment_type as a list."""

    def test_treatment_type_list_is_accepted(self, tmp_path, monkeypatch):
        """A list of canonical treatment_type values must be accepted."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(
            csv,
            experiment_overrides={"treatment_type": ["coculture", "darkness"]},
        )
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is True, "Expected validation to pass for list treatment_type"

    def test_treatment_type_list_with_invalid_value_rejected(self, tmp_path, monkeypatch):
        """A list containing a non-canonical value must fail."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(
            csv,
            experiment_overrides={"treatment_type": ["coculture", "fake_stress"]},
        )
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail for list with non-canonical value"
```

- [ ] **Step 2: Write failing test for background_factors validation**

Add after the `TestTreatmentTypeList` class:

```python
class TestBackgroundFactors:
    """Validator accepts and validates background_factors field."""

    def test_background_factors_list_accepted(self, tmp_path, monkeypatch):
        """A list of canonical background_factors values must be accepted."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(
            csv,
            experiment_overrides={
                "treatment_type": "nitrogen_stress",
                "background_factors": ["axenic", "continuous_light"],
                # Remove coculture-specific fields since treatment_type is not coculture
                "treatment_organism": None,
                "treatment_taxid": None,
            },
        )
        # Remove None-valued keys (treatment_organism, treatment_taxid)
        exp = config["publication"]["experiments"]["exp_coculture"]
        exp.pop("treatment_organism", None)
        exp.pop("treatment_taxid", None)
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is True, "Expected validation to pass for valid background_factors"

    def test_background_factors_with_invalid_value_rejected(self, tmp_path, monkeypatch):
        """A background_factors list with non-canonical value must fail."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(
            csv,
            experiment_overrides={
                "background_factors": ["axenic", "made_up_factor"],
            },
        )
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is False, "Expected validation to fail for non-canonical background_factors"

    def test_background_factors_optional(self, tmp_path, monkeypatch):
        """Missing background_factors should not cause validation failure."""
        monkeypatch.chdir(PROJECT_ROOT)
        csv = _write_minimal_csv(tmp_path)
        config = _make_valid_config(csv)
        cfg_file = _write_config(tmp_path, config)
        result = validate(str(cfg_file))
        assert result is True, "Expected validation to pass without background_factors"
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `pytest tests/test_paperconfig_validation.py::TestTreatmentTypeList -v && pytest tests/test_paperconfig_validation.py::TestBackgroundFactors -v`
Expected: FAIL (list treatment_type not yet handled, background_factors not yet validated)

- [ ] **Step 4: Extend CANONICAL_CONDITION_TYPES with new values**

In `validate_paperconfig.py`, lines 115-132, add the three new values to the set:

```python
CANONICAL_CONDITION_TYPES = {
    "growth_medium",
    "nitrogen_stress",
    "phosphorus_stress",
    "iron_stress",
    "salt_stress",
    "carbon_stress",
    "light_stress",
    "darkness",
    "plastic_stress",
    "viral",
    "coculture",
    "growth_state",
    "temperature_stress",
    "diel",
    "oxygen_stress",
    "axenic",
    "continuous_light",
    "diel_cycle",
}
```

- [ ] **Step 5: Update treatment_type validation to handle lists**

Replace lines 304-314 in `validate_paperconfig.py`:

```python
        # Canonical treatment_type (string or list)
        treatment_type = exp.get("treatment_type", "")
        if isinstance(treatment_type, list):
            for tt in treatment_type:
                if tt not in CANONICAL_CONDITION_TYPES:
                    errors.append(
                        _canonical_field_error(
                            config_path, f"experiments.{exp_key}",
                            "treatment_type", tt, CANONICAL_CONDITION_TYPES,
                        )
                    )
                else:
                    print(f"    treatment_type '{tt}': OK")
        elif treatment_type and treatment_type not in CANONICAL_CONDITION_TYPES:
            errors.append(
                _canonical_field_error(
                    config_path, f"experiments.{exp_key}",
                    "treatment_type", treatment_type, CANONICAL_CONDITION_TYPES,
                )
            )
        elif treatment_type:
            print(f"    treatment_type '{treatment_type}': OK")
```

- [ ] **Step 6: Add background_factors validation**

After the treatment_type validation block (after the code added in Step 5), add:

```python
        # Canonical background_factors (optional list)
        background_factors = exp.get("background_factors", [])
        if isinstance(background_factors, str):
            background_factors = [background_factors]
        for bf in background_factors:
            if bf not in CANONICAL_CONDITION_TYPES:
                errors.append(
                    _canonical_field_error(
                        config_path, f"experiments.{exp_key}",
                        "background_factors", bf, CANONICAL_CONDITION_TYPES,
                    )
                )
            else:
                print(f"    background_factors '{bf}': OK")
```

- [ ] **Step 7: Update coculture/viral consistency checks for list treatment_type**

Replace lines 363-388 in `validate_paperconfig.py`. The `tt` variable (line 364) is currently a string; it needs to handle lists:

```python
        # Coculture/viral consistency checks
        raw_tt = exp.get("treatment_type", "")
        tt_list = raw_tt if isinstance(raw_tt, list) else ([raw_tt] if raw_tt else [])
        t_org = exp.get("treatment_organism", "")
        if ("coculture" in tt_list or "viral" in tt_list) and not t_org:
            errors.append(
                f"{config_path} | experiments.{exp_key} | "
                f"treatment_type includes 'coculture' or 'viral' but missing treatment_organism"
            )
        if t_org and "coculture" not in tt_list and "viral" not in tt_list:
            warnings.append(
                f"{config_path} | experiments.{exp_key} | "
                f"has treatment_organism '{t_org}' but treatment_type {tt_list} "
                f"does not include 'coculture' or 'viral'"
            )
        if t_org and str(t_org).strip().lower() in ("phage", "virus", "bacteriophage") and "viral" not in tt_list:
            warnings.append(
                f"{config_path} | experiments.{exp_key} | "
                f"treatment_organism is '{t_org}' but treatment_type {tt_list} "
                f"does not include 'viral' (should it?)"
            )
        if ("coculture" in tt_list or "viral" in tt_list) and "treatment_taxid" not in exp:
            warnings.append(
                f"{config_path} | experiments.{exp_key} | "
                f"treatment_type includes 'coculture' or 'viral' but missing treatment_taxid"
            )
```

- [ ] **Step 8: Run tests to verify they pass**

Run: `pytest tests/test_paperconfig_validation.py -v`
Expected: ALL PASS (including the old tests — they still pass string values which are still accepted)

- [ ] **Step 9: Commit**

```bash
git add .claude/skills/paperconfig/validate_paperconfig.py tests/test_paperconfig_validation.py
git commit -m "feat: validate treatment_type as list, add background_factors validation"
```

---

### Task 3: Paperconfig skill documentation — Update templates and guidance

**Files:**
- Modify: `.claude/skills/paperconfig/SKILL.md:62,87,91-101,113-130,283-294,300-317`

- [ ] **Step 1: Update experiment template to show list format and background_factors**

In `SKILL.md`, replace line 62:
```yaml
      treatment_type: coculture       # canonical treatment category (see below)
```
with:
```yaml
      treatment_type: [coculture]     # canonical treatment category list (see below)
      background_factors: []          # experimental context factors list (see below, optional)
```

- [ ] **Step 2: Update required fields table**

In the Required Experiment Fields table (line 87), change the `treatment_type` row:
```
| `treatment_type` | Canonical treatment category (list) | See **Treatment Types** below |
```

- [ ] **Step 3: Add background_factors to optional fields table**

In the Optional Experiment Fields table (after line 101), add a row:
```
| `background_factors` | Experimental context factors not being compared in DE (list) | When conditions like coculture/axenic status or light regime are relevant but not the treatment variable. Values from same vocabulary as `treatment_type` plus: `axenic`, `continuous_light`, `diel_cycle` |
```

- [ ] **Step 4: Update Treatment Types table**

After the existing Treatment Types table (line 130), add the three new values:
```
| `diel` | Diel cycle as treatment | Circadian gene expression |
| `oxygen_stress` | Oxygen level changes | Microaerobic/anoxic conditions |
| `axenic` | Organism grown alone (background only) | Used in `background_factors`, not `treatment_type` |
| `continuous_light` | Standard continuous illumination (background only) | Used in `background_factors` |
| `diel_cycle` | Light:dark cycling regime (background only) | Used in `background_factors` |
```

Add a paragraph explaining the two-field model:
```
**`treatment_type` vs `background_factors`:** `treatment_type` (required, list) captures what the DE comparison tests. `background_factors` (optional, list) captures conditions present in the experiment but not compared. Both use the same vocabulary. Example: a darkness experiment run in coculture → `treatment_type: [darkness]`, `background_factors: [coculture]`. Axenic/coculture status always goes in `background_factors` unless coculture IS the DE comparison.
```

- [ ] **Step 5: Update gene_clusters template to include background_factors**

In the optional fields table for gene_clusters (lines 283-294), add a row after `treatment_type`:
```
| `background_factors` | str[] | Array of background condition factors (same vocabulary as `treatment_type`) |
```

In the gene_clusters example (lines 300-317), add after `treatment_type: ["nitrogen_stress"]`:
```yaml
  background_factors: ["axenic", "continuous_light"]
```

- [ ] **Step 6: Add "Extending the Vocabulary" section to SKILL.md**

After the `treatment_type` vs `background_factors` explanation paragraph, add a new section:

```markdown
#### Extending the Canonical Vocabulary

The canonical vocabulary is intentionally minimal. When creating a paperconfig for a new paper:

1. **Check existing values first.** Can this condition map to an existing category?
   - "UV exposure" → `light_stress`
   - "CO2 enrichment" → `carbon_stress`
   - "DCMU treatment" → describe in `treatment_condition`/`experimental_context`; use closest category or leave `background_factors` empty
2. **Prefer general categories over specific ones.** `chemical_stress` is better than `dcmu_inhibitor`. The specifics go in `treatment_condition` and `experimental_context`.
3. **If no existing value fits**, flag it to the user: "This paper studies [X], which doesn't map cleanly to any existing canonical value. Closest match: `Y`. Should I use `Y` or propose a new value?"
4. **User decides.** If a new value is approved, update all three locations in a single commit:
   - `CANONICAL_CONDITION_TYPES` in `validate_paperconfig.py`
   - Treatment Types table in this SKILL.md
   - `CLAUDE.md` key graph facts section
5. **Never invent new canonical values without user approval.**
```

- [ ] **Step 7: Commit**

```bash
git add .claude/skills/paperconfig/SKILL.md
git commit -m "docs: update paperconfig skill for list treatment_type and background_factors"
```

---

### Task 4: Omics adapter — Read treatment_type as list, read background_factors

**Files:**
- Modify: `multiomics_kg/adapters/omics_adapter.py:278` (treatment_type property), `:290` (after table_scope_detail)

- [ ] **Step 1: Update treatment_type to normalize as list**

In `omics_adapter.py`, replace line 278:
```python
                    "treatment_type": self.clean_text(exp.get("treatment_type", "")),
```
with:
```python
                    "treatment_type": self._normalize_list_field(exp, "treatment_type"),
```

- [ ] **Step 2: Add background_factors to experiment properties**

After line 290 (`"table_scope_detail": ...`), add:
```python
                    "background_factors": self._normalize_list_field(exp, "background_factors"),
```

- [ ] **Step 3: Add the _normalize_list_field helper method**

Add this method to the `OMICSAdapter` class (before `get_nodes`):

```python
    def _normalize_list_field(self, exp: dict, field: str) -> list[str]:
        """Normalize a paperconfig field to a list of cleaned strings."""
        val = exp.get(field, [])
        if isinstance(val, str):
            val = [val] if val else []
        return [self.clean_text(v) for v in val]
```

- [ ] **Step 4: Run existing omics adapter tests**

Run: `pytest tests/test_omics_adapter*.py -v -m "not slow and not kg"`
Expected: PASS (backward compatible — existing paperconfigs still have string treatment_type which gets normalized to single-element list)

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/omics_adapter.py
git commit -m "feat: omics adapter reads treatment_type as list, adds background_factors"
```

---

### Task 5: Cluster adapter — Read background_factors

**Files:**
- Modify: `multiomics_kg/adapters/cluster_adapter.py:161-175` (ClusteringAnalysis props), `:195-217` (GeneCluster props)

- [ ] **Step 1: Add background_factors to ClusteringAnalysis node**

In `cluster_adapter.py`, after line 164 (`treatment_type = [treatment_type]`), add:
```python
            # background_factors may be a list or a string
            background_factors = table.get("background_factors", [])
            if isinstance(background_factors, str):
                background_factors = [background_factors]
```

Then in the `analysis_props` dict (after line 175 `"treatment_type": treatment_type,`), add:
```python
                "background_factors": background_factors,
```

- [ ] **Step 2: Add background_factors to GeneCluster nodes**

In the `props` dict for GeneCluster (after line 201 `"treatment_type": treatment_type,`), add:
```python
                    "background_factors": background_factors,
```

- [ ] **Step 3: Run cluster adapter tests**

Run: `pytest tests/test_cluster_adapter.py -v -m "not slow and not kg"`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add multiomics_kg/adapters/cluster_adapter.py
git commit -m "feat: cluster adapter reads background_factors for ClusteringAnalysis and GeneCluster"
```

---

### Task 6: Post-import scripts — Aggregate background_factors, fix treatment_type for arrays

**Files:**
- Modify: `scripts/post-import.sh:70,85,100-114,186-200`
- Modify: `scripts/post-import.cypher` (same Cypher, kept in sync)

- [ ] **Step 1: Fix Publication aggregation for array treatment_type + add background_factors**

In `scripts/post-import.sh`, replace lines 100-114:

```bash
echo "=== Post-process: Compute Publication summary properties ==="
cypher-shell <<'CYPHER'
MATCH (p:Publication)
OPTIONAL MATCH (p)-[:Has_experiment]->(e:Experiment)
WITH p,
     count(e) AS ec,
     [x IN collect(DISTINCT e.omics_type) WHERE x IS NOT NULL] AS ots,
     [x IN collect(DISTINCT e.organism_name) WHERE x IS NOT NULL] AS orgs,
     [x IN collect(DISTINCT e.coculture_partner) WHERE x IS NOT NULL AND x <> ''] AS coculture_orgs,
     apoc.coll.toSet(reduce(s = [], t IN collect(coalesce(e.treatment_type, [])) | s + t)) AS tts,
     apoc.coll.toSet(reduce(s = [], t IN collect(coalesce(e.background_factors, [])) | s + t)) AS bfs
SET p.experiment_count = ec,
    p.treatment_types = apoc.coll.sort(tts),
    p.background_factors = apoc.coll.sort(bfs),
    p.omics_types = apoc.coll.sort(ots),
    p.organisms = apoc.coll.sort(apoc.coll.toSet(orgs + coculture_orgs))
CYPHER
```

Note: Uses `reduce` + `apoc.coll.toSet` to flatten array properties (same pattern as OrganismTaxon aggregation), avoiding the cross-product issue of double UNWIND.

- [ ] **Step 2: Fix OrganismTaxon aggregation to include background_factors**

In `scripts/post-import.sh`, replace lines 186-200:

```bash
echo "--- publication_count, experiment_count, treatment_types, omics_types, background_factors ---"
cypher-shell <<'CYPHER'
MATCH (o:OrganismTaxon)
OPTIONAL MATCH (p:Publication)
  WHERE ANY(org IN p.organisms WHERE org = o.preferred_name)
WITH o,
     count(DISTINCT p) AS pc,
     CASE WHEN count(p) > 0 THEN sum(p.experiment_count) ELSE 0 END AS ec,
     apoc.coll.toSet(reduce(s = [], t IN collect(p.treatment_types) | s + t)) AS tts,
     apoc.coll.toSet(reduce(s = [], t IN collect(p.omics_types) | s + t)) AS ots,
     apoc.coll.toSet(reduce(s = [], t IN collect(coalesce(p.background_factors, [])) | s + t)) AS bfs
SET o.publication_count = pc,
    o.experiment_count = ec,
    o.treatment_types = tts,
    o.omics_types = ots,
    o.background_factors = bfs
CYPHER
```

- [ ] **Step 3: Add background_factors index**

In `scripts/post-import.sh`, after line 70 (experiment_treatment_type_idx), add:
```cypher
CREATE INDEX experiment_background_factors_idx IF NOT EXISTS FOR (e:Experiment) ON (e.background_factors);
```

- [ ] **Step 4: Apply identical changes to post-import.cypher**

Mirror all three changes (Publication aggregation, OrganismTaxon aggregation, index) in `scripts/post-import.cypher`. The Cypher content must be identical between the two files.

- [ ] **Step 5: Commit**

```bash
git add scripts/post-import.sh scripts/post-import.cypher
git commit -m "post-import: aggregate background_factors, fix treatment_type for arrays"
```

---

### Task 7: Paperconfig pass 1 — Mechanical conversion

**Files:**
- Modify: All 23 paperconfig.yaml files listed in `data/Prochlorococcus/papers_and_supp/paperconfig_files.txt`

- [ ] **Step 1: Convert all treatment_type from string to single-element list**

For every experiment in every paperconfig, change:
```yaml
      treatment_type: nitrogen_stress
```
to:
```yaml
      treatment_type: [nitrogen_stress]
```

- [ ] **Step 2: Add obvious background_factors**

Apply these rules mechanically:
- If experiment has `treatment_organism` set AND `treatment_type` is NOT `[coculture]` or `[viral]`: add `background_factors: [coculture]`
- If experiment has NO `treatment_organism` AND `treatment_type` is NOT `[coculture]` or `[viral]`: add `background_factors: [axenic]`
- If experiment `light_condition` contains "continuous" and treatment_type is NOT light-related (`[light_stress]`, `[darkness]`, `[diel]`): add `continuous_light` to background_factors
- If experiment has `treatment_type: [coculture]` or `[viral]`: do NOT add `axenic` or `coculture` to background_factors (it's already in treatment_type)

- [ ] **Step 3: Flag experiments where background conditions don't map to current vocabulary**

For any experiment where a notable background condition has no canonical match, leave `background_factors` empty (or partial) and add a YAML comment:
```yaml
      background_factors: [axenic]  # TODO: DCMU photosynthesis inhibitor present, no canonical match
```

Known cases to flag:
- **Steglich 2006**: DCMU + white light experiment (chemical inhibitor as background)
- Any other experiments where `experimental_context` describes conditions not covered by the current vocabulary

- [ ] **Step 4: Apply known multi-factor experiments from brainstorming analysis**

These experiments need specific `background_factors` based on the paper analysis (coculture/light context is in `experimental_context` free text, not in structured fields):
- **Biller 2018**: `darkness_extended_darkness_natl2a_rnaseq_coculture` → `background_factors: [coculture]`; `darkness_extended_darkness_mit1002_rnaseq` → `background_factors: [coculture]`; axenic one → `background_factors: [axenic]`
- **Coe 2024**: coculture experiments under diel cycle → `background_factors: [coculture, diel_cycle]`
- **Hennon 2017**: carbon stress experiments — MIT9312 and EZ55 coculture → `background_factors: [coculture, diel_cycle]`; EZ55 axenic → `background_factors: [axenic, diel_cycle]`
- **Thompson 2016**: `light_stress_..._infected` → `background_factors: [viral]`; `light_stress_..._dark` → `background_factors: [axenic]`; `viral_..._light` → `background_factors: [continuous_light]`; `viral_..._dark` → `background_factors: [darkness]`
- **Lin 2015**: phosphorus stress under diel cycle — infected experiment → `background_factors: [viral, diel_cycle]`; uninfected → `background_factors: [axenic, diel_cycle]`
- **Weissberg 2025**: N-starvation coculture experiments → `background_factors: [coculture]`; axenic ones → `background_factors: [axenic]`; coculture treatment_type experiments → no axenic/coculture in background

- [ ] **Step 5: Run validator on all paperconfigs**

Run: `uv run python -c "from multiomics_kg.skills.paperconfig.validate_paperconfig import validate; import sys; [validate(line.strip()) for line in open('data/Prochlorococcus/papers_and_supp/paperconfig_files.txt') if line.strip() and not line.strip().startswith('#')]"`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add data/Prochlorococcus/papers_and_supp/
git commit -m "paperconfig: convert treatment_type to lists, add background_factors (pass 1)"
```

---

### Task 8: Paperconfig pass 2 — Review flagged experiments and vocabulary gaps

**Files:**
- Modify: Paperconfig files with TODO comments from Task 7 Step 3
- Possibly modify: `validate_paperconfig.py`, `SKILL.md`, `CLAUDE.md` (if vocabulary is extended)

- [ ] **Step 1: Collect all TODO-flagged experiments**

Grep all paperconfigs for `# TODO:` comments added in pass 1. Present a summary table to the user:

| Paper | Experiment | Condition | Closest canonical | Proposal |
|---|---|---|---|---|
| Steglich 2006 | DCMU + light | DCMU inhibitor | `light_stress`? | Use existing or add new? |
| Barreto 2022 | carbon_stress experiments | DE pools axenic + coculture samples | `axenic`? `coculture`? both? | Ambiguous — user decides |

- [ ] **Step 2: User decides on vocabulary extensions**

For each flagged condition, the user chooses:
- **(a)** Map to existing canonical value — update the paperconfig
- **(b)** Approve a new canonical value — then update `CANONICAL_CONDITION_TYPES` in `validate_paperconfig.py`, Treatment Types table in `SKILL.md`, and `CLAUDE.md`
- **(c)** Leave `background_factors` empty for this experiment — remove the TODO comment, accept that this condition is only in free-text `experimental_context`

- [ ] **Step 3: If vocabulary was extended, commit canonical set changes first**

```bash
git add .claude/skills/paperconfig/validate_paperconfig.py .claude/skills/paperconfig/SKILL.md CLAUDE.md
git commit -m "feat: extend canonical vocabulary with [new_value]"
```

- [ ] **Step 4: Fill in remaining background_factors and remove TODO comments**

Apply user decisions to all flagged paperconfigs.

- [ ] **Step 5: Run validator on all paperconfigs**

Run: `pytest tests/test_paperconfig_validation.py -v`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add data/Prochlorococcus/papers_and_supp/
git commit -m "paperconfig: finalize background_factors after manual review (pass 2)"
```

---

### Task 9: KG validity tests — Update for array treatment_type and background_factors

**Files:**
- Modify: `tests/kg_validity/test_expression.py:267-277`
- Modify: `tests/kg_validity/test_post_import.py:552-561,674-688`

- [ ] **Step 1: Update test_experiment_has_treatment_type for array property**

In `tests/kg_validity/test_expression.py`, replace lines 267-277:

```python
def test_experiment_has_treatment_type(run_query):
    """Every Experiment node must have a non-empty treatment_type list."""
    result = run_query("""
        MATCH (exp:Experiment)
        WHERE exp.treatment_type IS NULL OR size(exp.treatment_type) = 0
        RETURN count(exp) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, (
        f"{missing} Experiment nodes are missing treatment_type"
    )
```

- [ ] **Step 2: Update coculture_partner consistency check for list treatment_type**

In `tests/kg_validity/test_post_import.py`, replace lines 674-688:

```python
def test_experiment_coculture_partner_null_when_no_partner(run_query):
    """Experiments without a treatment organism should have null coculture_partner.

    Both coculture and viral experiments have a treatment organism (partner).
    All other treatment types should have null coculture_partner.
    """
    result = run_query("""
        MATCH (e:Experiment)
        WHERE NOT 'coculture' IN e.treatment_type
          AND NOT 'viral' IN e.treatment_type
          AND e.coculture_partner IS NOT NULL
        RETURN count(e) AS bad, collect(e.id)[..5] AS examples
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} non-partner experiments have coculture_partner set: {result[0]['examples']}"
    )
```

- [ ] **Step 3: Add test for background_factors no nulls on Publication**

In `tests/kg_validity/test_post_import.py`, after `test_publication_treatment_types_no_nulls`:

```python
def test_publication_background_factors_no_nulls(run_query):
    """background_factors list should contain no null entries."""
    result = run_query("""
        MATCH (p:Publication)
        WHERE p.background_factors IS NOT NULL
          AND ANY(x IN p.background_factors WHERE x IS NULL)
        RETURN count(p) AS bad
    """)
    assert result[0]["bad"] == 0, (
        f"{result[0]['bad']} publications have null entries in background_factors"
    )
```

- [ ] **Step 4: Add test for treatment_type values are canonical**

In `tests/kg_validity/test_expression.py`, after `test_experiment_has_treatment_type`:

```python
def test_experiment_treatment_type_values_canonical(run_query):
    """All treatment_type values should be from the canonical vocabulary."""
    result = run_query("""
        MATCH (e:Experiment)
        UNWIND e.treatment_type AS tt
        WITH DISTINCT tt
        RETURN collect(tt) AS all_types
    """)
    known = {
        "nitrogen_stress", "phosphorus_stress", "iron_stress", "carbon_stress",
        "salt_stress", "oxygen_stress", "temperature_stress", "light_stress",
        "darkness", "plastic_stress", "viral", "coculture", "growth_state",
        "growth_medium", "diel", "axenic", "continuous_light", "diel_cycle",
    }
    actual = set(result[0]["all_types"])
    unknown = actual - known
    assert not unknown, f"Unknown treatment_type values in KG: {unknown}"
```

- [ ] **Step 5: Commit**

```bash
git add tests/kg_validity/test_expression.py tests/kg_validity/test_post_import.py
git commit -m "test: update KG validity tests for array treatment_type and background_factors"
```

---

### Task 10: Update CLAUDE.md

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Update Experiment node documentation**

In the paperconfig template section, update to show list format for `treatment_type` and add `background_factors`.

- [ ] **Step 2: Update Key graph facts**

Update the Experiment node properties description to reflect `treatment_type: str[]` and `background_factors: str[]`.

- [ ] **Step 3: Update canonical values list if present**

Ensure any mention of canonical treatment_type values includes the three new ones.

- [ ] **Step 4: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: update CLAUDE.md for treatment_type list and background_factors"
```
