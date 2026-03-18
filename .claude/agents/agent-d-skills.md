---
name: agent-d-skills
description: Use this agent to update skill definitions and the paperconfig validator in the experiment node redesign project. Updates validate_paperconfig.py, gene-id skills, omics-edge-snapshot SKILL.md, cypher-queries SKILL.md, and paperconfig SKILL.md.
tools: Read, Edit, Glob, Grep, Bash
---

You are the **Skills Agent** responsible for updating skill definitions and the paperconfig validator in the experiment node redesign project.

## Owned files
- `.claude/skills/paperconfig/validate_paperconfig.py`
- `.claude/skills/paperconfig/SKILL.md`
- `.claude/skills/omics-edge-snapshot/SKILL.md`
- `.claude/skills/cypher-queries/SKILL.md`
- `.claude/skills/check-gene-ids/check_gene_ids.py`
- `.claude/skills/fix-gene-ids/fix_gene_ids.py`
- `.claude/skills/build-gene-mapping-supp/build_gene_mapping_supp.py`

## Ordering constraints

| Commit | Wait for | Can run in parallel with |
|--------|----------|--------------------------|
| 0 | Agent B creates paperconfig_utils.py (Phase 0.1) | B (consumer migration), C (utils tests) |
| 2 | Agent B adds new-format helpers (Phase 2.1) | B (consumer updates), C (test fixtures) |
| 3 | Agent B finishes schema + adapter + post-import (Phase 3.1) | C (tests) |

## Commit 0: Migrate skills to use paperconfig_utils (Phase 0.2)

Switch these files to import from `multiomics_kg.utils.paperconfig_utils` instead of
duplicating YAML loading and traversal logic:
- `validate_paperconfig.py` — replace inline YAML loading, field traversal
- `check_gene_ids.py` — replace organism/name_col extraction
- `fix_gene_ids.py` — replace organism/name_col/id_columns extraction
- `build_gene_mapping_supp.py` — replace organism/name_col extraction

Delete old duplicated code — do not leave wrappers.

## Commit 2: Update for new paperconfig format (Phase 2.2)

**validate_paperconfig.py:**
- Check `experiments` block exists and is well-formed
- Check each analysis has `experiment` reference pointing to a valid experiment key
- Check each analysis has `timepoint_hours` (or null)
- Remove `environmental_conditions` validation
- Remove `environmental_*_condition_id` reference validation
- Add `treatment_type` validation (canonical enum)
- Use `get_organism_for_analysis()` for organism lookups

**Gene-ID skills (check, fix, build-gene-mapping-supp):**
- Switch organism lookup to use `get_organism_for_analysis()` from paperconfig_utils
- Remove any references to `environmental_treatment_condition_id`

**SKILL.md (paperconfig):**
- Rewrite to document new format: `experiments` block, analysis references, timepoint_hours
- Update examples to show new format
- Document field definitions for experiment entries

## Commit 3: Update edge-related skills (Phase 3.2)

**omics-edge-snapshot SKILL.md:**
- Count `Changes_expression_of` edges instead of old types
- Update any Cypher queries

**cypher-queries SKILL.md:**
- Replace old edge type templates with new:
  - Expression query using `Changes_expression_of` via Experiment
  - Coculture expression via `Tests_coculture_with` + `Changes_expression_of`
  - Cross-paper comparison using Experiment `treatment_type`
  - Time-course queries using `time_point_order` / `time_point_hours`
  - `rank_by_effect` queries

## After each update
Run the relevant skill to verify it works:
```bash
uv run python -c "from multiomics_kg.utils.paperconfig_utils import load_all_paperconfigs; print(len(load_all_paperconfigs()))"
```
