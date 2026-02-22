---
name: gene-protein-quality
description: Report gene and protein annotation quality per organism strain. Shows coverage of proteins, functional descriptions, GO terms, KEGG, eggNOG, and expression data. Tracks changes over time via snapshots. Use to monitor annotation completeness across strains.
argument-hint: [--save | --diff | --strain <name>]
user-invocable: true
allowed-tools: Read, Bash(uv *), Bash(python *), Bash(docker exec *)
---

# Gene/Protein Quality Report Skill

Per-strain report of annotation completeness across all gene nodes in the knowledge graph.

## Quick Start

```bash
# Full report, all strains
uv run python .claude/skills/gene-protein-quality/gene_protein_quality.py

# Single strain
uv run python .claude/skills/gene-protein-quality/gene_protein_quality.py --strain MED4

# Save current state as snapshot (for future diffing)
uv run python .claude/skills/gene-protein-quality/gene_protein_quality.py --save

# Compare current state to saved snapshot
uv run python .claude/skills/gene-protein-quality/gene_protein_quality.py --diff

# Save + diff in one pass (useful after a rebuild)
uv run python .claude/skills/gene-protein-quality/gene_protein_quality.py --diff --save
```

## Metrics Reported (per strain)

| Column | Source | What it measures |
|--------|--------|-----------------|
| `Genes` | NCBI/Cyanorak GFF | Total gene nodes |
| `Protein%` | UniProt (`Gene_encodes_protein` edge) | Genes with a linked UniProt protein |
| `Function%` | UniProt (denorm → `function_description`) | Genes with a non-empty functional description |
| `GO_BP%` | UniProt+Cyanorak (`go_biological_processes`) | Genes with ≥1 GO biological process term |
| `KEGG%` | Cyanorak/NCBI (`kegg` property) | Genes with a KEGG annotation |
| `eggNOG%` | Cyanorak (`eggNOG` property) | Genes with an eggNOG cluster assignment |
| `CyanoRole%` | Cyanorak (`cyanorak_Role`) | Genes with a Cyanorak functional role |
| `Product%` | NCBI/Cyanorak (`product`/`product_cyanorak`) | Genes with any product name |
| `Domains%` | InterPro/Cyanorak (`protein_domains`) | Genes with protein domain annotations |
| `Expr%` | OMICS papers (`Affects_expression_of` edges) | Genes with ≥1 direct expression measurement |

## Snapshot Files

Snapshots are saved to `.claude/skills/gene-protein-quality/quality_snapshot.json`. Commit this file to track annotation quality over time.

## Workflow

When invoked:

1. Connect to Neo4j via `docker exec deploy cypher-shell -u neo4j -p neo4j`
2. Run per-strain aggregation queries
3. Print a formatted table (absolute counts + % coverage)
4. If `--diff`: compare against saved snapshot and highlight changes
5. If `--save`: write current metrics to `quality_snapshot.json`
