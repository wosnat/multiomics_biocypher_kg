# Schema Improvements Subagent Spec

Agents for executing `plans/schema_improvements_for_mcp.md`.

## Agents

| Agent | File | Role |
|-------|------|------|
| H | `agent-h-plan-manager.md` | Create phase sub-plans, track decisions, record gate results |
| G | `agent-g-code-review.md` | Read-only review at each phase gate |
| F | `agent-f-validation-runner.md` | Run skills/tests at phase gates |
| A | `agent-a-paperconfig.md` | Edit all 24 paperconfig.yaml files |
| B | `agent-b-schema-adapters.md` | Edit schema_config.yaml, omics_adapter.py, cyanorak_ncbi_adapter.py, post-import.sh |
| C | `agent-c-tests.md` | Write and run unit + KG validity tests |
| D | `agent-d-skills.md` | Update validate_paperconfig.py and SKILL.md files |
| E | `agent-e-docs.md` | Update CLAUDE.md, MEMORY.md |

## Phase execution order

```
Phase 1 → Phase 2 → Phase 3 → Phase 4 → Phase 5 → Phase 6
```
Never start phase N+1 until Agent G has approved phase N gate.

## Per-phase agent ordering

### Phase 1: Paperconfig standardization
```
H (phase sub-plan)
  → D (update validate_paperconfig.py)  ← MUST BE FIRST
      → A (standardize paperconfigs)    ← parallel with C
      → C (extend test_paperconfig_validation.py)
          → G (review)
              → F (run tests + grep for old values)
                  → COMMIT
```

### Phase 2: Missing environmental conditions
```
H (phase sub-plan)
  → A (add env conditions to 3 paperconfigs)
      → C (run tests)
          → G (review)
              → F (run tests + check env node count)
                  → COMMIT
```

### Phase 3: OrganismTaxon preferred names
```
H (phase sub-plan)
  → B (update cyanorak_ncbi_adapter.py)
      → C (add preferred_name test + run)
          → G (review)
              → COMMIT
```

### Phase 4: Edge type split + new properties
```
F (save pre_phase4 snapshot — BEFORE any code changes)

H (phase sub-plan)
  → B-schema (update schema_config.yaml)
      → B-adapter (update omics_adapter.py)
          → C (update + run test_omics_adapter_organism_gene.py)  ┐ parallel
          → D (update omics-edge-snapshot, cypher-queries, paperconfig skills)  ┤ parallel
          → E (update CLAUDE.md edge labels section)  ┘ parallel
              → G (review)
                  → full KG build
                      → F (compare snapshot, verify counts)
                          → COMMIT
```

### Phase 5: Post-import homolog propagation
```
H (phase sub-plan)
  → B (update post-import.sh — 2 queries with cross-phylum filter)
      → C (update test_post_import.py)
          → G (review)
              → Docker build (overnight)
                  → F (save post_phase5 snapshot, verify ortholog counts, run KG tests)
                      → COMMIT
```

### Phase 6: Final validation + docs
```
H (phase sub-plan)
  → C (update test_expression.py + regenerate snapshot_data.json)  ┐ parallel
  → D (final skill polish)  ┤ parallel
  → E (update CLAUDE.md counts + MEMORY.md)  ┘ parallel
      → G (full codebase sweep — grep for old labels)
          → F (run all 10 verification queries + full test suite + /check-gene-ids)
              → COMMIT + TAG RELEASE
```

## Parallel execution rules

- Within a phase, agents marked "parallel" can run simultaneously
- No agent should edit a file owned by another agent
- If two agents need to read the same file, that is fine
- Agent G and Agent F never edit files
- Agent H is the only agent that edits plan files
