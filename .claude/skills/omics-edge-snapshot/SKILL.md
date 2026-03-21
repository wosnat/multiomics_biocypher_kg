---
name: omics-edge-snapshot
description: Snapshot and compare expression edge counts (Changes_expression_of from Experiment nodes) in the Neo4j knowledge graph. Use before and after omics adapter changes to verify no edges were lost. Detects per-paper regressions (lost edges) vs improvements (gained edges from better ID resolution). Backward compatible with old snapshots that used condition_edges/coculture_edges.
argument-hint: [--save NAME | --compare NAME | --list]
user-invocable: true
allowed-tools: Read, Bash(uv *), Bash(docker exec *), Bash(docker compose *)
---

# Omics Edge Snapshot Skill

Capture and compare expression edge counts per publication before and after changes to the omics adapter or prepare_data pipeline. The graph uses a single edge type for expression data:

- `Changes_expression_of` -- source is an `Experiment` node, target is a `Gene` node

Publications are linked to experiments via `Has_experiment` edges (`Publication -> Experiment`). Ensures that rebuilds don't silently lose valid edges.

## Quick Start

```bash
# Step 1 -- before the rebuild: capture baseline
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --save baseline

# ... rebuild: docker compose up -d --build, or create_knowledge_graph.py + reimport ...

# Step 2 -- after the rebuild: compare
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --compare baseline

# List saved snapshots
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --list

# Compare two saved snapshots (no live Neo4j required)
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --compare before --against after
```

## What It Captures

For each snapshot:

| Metric | Description |
|--------|-------------|
| `total_edges` | Count of `Changes_expression_of` edges (Experiment -> Gene) |
| `per_publication` | Edge count per DOI (via `Publication -[:Has_experiment]-> Experiment -[:Changes_expression_of]-> Gene`) |
| `per_publication_by_direction` | Breakdown into `up` / `down` per publication |
| `by_organism` | Counts grouped by the target gene's `organism_strain` property |
| `per_publication_genes` | Set of locus_tags with edges per publication (for gene-level diff) |

## Key Cypher Queries

```cypher
-- Total edges
MATCH ()-[r:Changes_expression_of]->() RETURN count(r) AS total

-- Per publication
MATCH (pub:Publication)-[:Has_experiment]->(e:Experiment)-[r:Changes_expression_of]->(g:Gene)
RETURN pub.doi AS publication, count(r) AS edges

-- By organism
MATCH (e:Experiment)-[r:Changes_expression_of]->(g:Gene)
RETURN g.organism_strain AS organism, count(r) AS edges
```

## Interpreting the Report

| Outcome | Meaning | Action |
|---------|---------|--------|
| **No regressions** | All publications kept or gained edges | Proceed |
| **IMPROVED** | Publication gained edges | Better ID resolution -- expected after resolve stage |
| **NEW** | New publication appeared | New data added |
| **LOST** (regression) | Publication lost edges | Investigate -- check `_resolved.csv`, re-run `/check-gene-ids` |
| **GONE** (regression) | Publication disappeared entirely | Investigate -- paperconfig issue or adapter bug |

**Exit code**: 0 if no regressions, 1 if any publication lost edges or disappeared.

### Dangling edge removal is expected

When the `resolve_paper_ids.py` stage (prepare_data step 4) creates `_resolved.csv` files, the omics adapter skips rows where `locus_tag` is NaN instead of creating dangling edges. This means:
- **Total edges may decrease** (dangling edges removed) -- this is correct behaviour
- **Per-publication matched counts should stay the same or increase** (better resolution)

### Backward compatibility

When comparing against old snapshots that stored `condition_edges` and `coculture_edges` separately (from the pre-Experiment-node era), the tool sums them to compute the old total. The comparison works correctly across format versions.

## Snapshot Files

Snapshots are saved to `.claude/skills/omics-edge-snapshot/snapshots/<name>.json`.

Suggested naming convention:
- `baseline` -- the current production graph before any changes
- `pre_experiment_redesign` -- last snapshot before the Experiment node migration
- `post_resolve` -- after running prepare_data steps 3+4
- `post_rebuild` -- after a full Docker rebuild

## Workflow

When invoked (e.g., `/omics-edge-snapshot --save baseline`):

1. If `--save NAME`: run the script to query live Neo4j and write `snapshots/NAME.json`
2. If `--compare NAME`: query live Neo4j for current state, load the saved snapshot, print comparison
3. If `--compare NAME --against OTHER`: load both snapshots and compare (no live Neo4j needed)
4. Check exit code: 0 = OK, 1 = regressions found
5. For any LOST/GONE publication, run `/check-gene-ids "Paper Name"` to diagnose

## Connection

Uses `docker exec deploy cypher-shell -u neo4j -p neo4j` (same as `gene-protein-quality` skill). Requires the Docker stack to be running:

```bash
docker compose up -d
```
