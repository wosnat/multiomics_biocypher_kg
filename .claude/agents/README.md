# Experiment Node Redesign — Subagent Spec

Agents for executing `plans/experiment_node_redesign.md`.

## Agents

| Agent | File | Role |
|-------|------|------|
| H | `agent-h-plan-manager.md` | Track decisions, record gate results, update plan |
| G | `agent-g-code-review.md` | Read-only review at each commit gate |
| F | `agent-f-validation-runner.md` | Run tests, skills, snapshots at gates |
| A | `agent-a-paperconfig.md` | Run migration script on all 26 paperconfig.yaml files |
| B | `agent-b-schema-adapters.md` | Core code: paperconfig_utils.py, schema, adapter, post-import, migration script |
| C | `agent-c-tests.md` | Write and run unit + KG validity tests |
| D | `agent-d-skills.md` | Update validate_paperconfig.py, skill SKILL.md files, gene-id skills |
| E | `agent-e-docs.md` | Update CLAUDE.md, memory files |

## Commit execution order

```
PRE-WORK → Commit 0 → Commit 1 → Commit 2 → Commit 3 → Rebuild → Commit 4
```
Never start commit N+1 until Agent G has approved commit N gate.

## Per-commit agent ordering

### PRE-WORK
```
H (update agent definitions)
  → F (save "before" snapshot)
      → Manual: git tag
```

### Commit 0: paperconfig_utils.py (pure refactor)
```
B (create paperconfig_utils.py)
  → B + D + C (parallel — migrate consumers, write utils tests)
      → F (run tests + prepare_data.sh)
          → G (review)
              → Manual: COMMIT
```

### Commit 1: dry-run migration
```
B (write migration script)
  → B (run --dry-run)
      → Manual: inspect output
          → G (review)
              → Manual: COMMIT
```

### Commit 2: apply migration + update consumers (ATOMIC)
```
A (run migration) + B (add new-format helpers to utils)
  → B + D + C (parallel — update consumers, skills, test fixtures)
      → F (run tests + prepare_data.sh)
          → G (review)
              → Manual: COMMIT
```

### Commit 3: adapter + schema + post-import
```
B (schema + adapter rewrite + post-import)
  → D + C (parallel — update skills, write tests)
      → F (run tests)
          → G (review)
              → Manual: COMMIT
```

### Rebuild + Validate
```
Manual: docker compose up
  → F (pytest -m kg + /omics-edge-snapshot compare)
      → Manual: spot-checks
```

### Commit 4: docs and cleanup
```
E + C (parallel — CLAUDE.md + snapshot regeneration)
  → G (review)
      → Manual: COMMIT
```

## Parallel execution rules

- Within a commit, agents marked "parallel" can run simultaneously
- No agent should edit a file owned by another agent
- If two agents need to read the same file, that is fine
- Agent G and Agent F never edit files
- Agent H is the only agent that edits plan files
