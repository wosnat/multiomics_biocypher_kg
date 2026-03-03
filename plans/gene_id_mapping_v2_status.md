# Gene ID Mapping v2 — Implementation Status

## Plan Reference
Full approved plan: `/home/osnat/.claude/plans/sleepy-splashing-kettle.md`

## What Was Implemented

### Core files (all complete)

| File | Status | Notes |
|------|--------|-------|
| `multiomics_kg/download/gene_id_graph.py` | ✅ DONE | GeneIdGraph class, three-tier classification, iterative convergence, normalize_id, build_diagnostic_report, to_json_structure |
| `multiomics_kg/download/build_gene_id_mapping.py` | ✅ DONE | Rewritten to use GeneIdGraph; collects all source rows, process_all_rows(), writes v2 JSON + diagnostic report; NO gene_mapping_supp.csv |
| `multiomics_kg/utils/gene_id_utils.py` | ✅ DONE | Added MappingData, load_mapping_v2(), expand_list(), _heuristic_candidates(), resolve_row() |
| `multiomics_kg/download/resolve_paper_ids.py` | ✅ DONE | Uses load_mapping_v2() + resolve_row(), multi-column fallback, writes _resolved_report.txt |
| `tests/test_gene_id_graph.py` | ✅ DONE | 40 unit tests, all passing (897 total tests pass) |
| `.claude/skills/fix-gene-ids/fix_gene_ids.py` | ✅ DONE | Updated to load_mapping_v2() + resolve_row(); fallback to legacy |
| `.claude/skills/fix-gene-ids/SKILL.md` | ✅ DONE | Updated description |
| `.claude/skills/check-gene-ids/check_gene_ids.py` | ✅ DONE | Added check_gene_id_mapping_v2(), new USE_STEP4 / AMBIGUOUS_IN_MAPPING fix strategies |
| `docs/methods_gene_id_mapping.md` | ✅ DONE | ~500-word scientific methods section |
| `CLAUDE.md` | ✅ DONE | Gene ID Mapping section updated with tier table, v2 schema, method strings |

## What Remains

### Phase 1 deployment: MIT9312 (NOT YET STARTED)

The snapshot step was interrupted. Steps to complete:

```bash
# Step 1 — snapshot before rebuild (Neo4j must be running)
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --save before_mit9312

# Step 2 — rebuild gene_id_mapping.json for MIT9312
uv run python -m multiomics_kg.download.build_gene_id_mapping --strains MIT9312 --force

# Review output:
# - cache/data/Prochlorococcus/genomes/MIT9312/gene_id_mapping.json (v2)
# - cache/data/Prochlorococcus/genomes/MIT9312/gene_id_mapping_report.json

# Step 3 — re-resolve paper CSVs for MIT9312 papers
uv run python -m multiomics_kg.download.resolve_paper_ids --force
# (Only MIT9312 papers will actually re-run; others are already up to date)

# Review _resolved_report.txt files in:
# data/Prochlorococcus/papers_and_supp/barreto 2022/
# data/Prochlorococcus/papers_and_supp/tetu 2019/
# data/Prochlorococcus/papers_and_supp/Hennon 2017/

# Step 4 — verify match rates
uv run python .claude/skills/check-gene-ids/check_gene_ids.py

# Step 5 — rebuild KG + compare snapshot
docker compose up -d --build
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --compare before_mit9312
```

MIT9312 papers with known non-standard IDs:
- **barreto 2022**: uses gene symbols (Tier 3 gene_name) — need Tier 1 from annotations
- **tetu 2019**: uses 2017-format RS locus tags (MIT9312_RS…) — need old_locus_tag mapping from NCBI GFF
- **Hennon 2017**: uses CDS accession IDs (`lcl|CP000111.1_cds_*`) — likely need annotation_gff entry

### Phase 2-N: Remaining strains (NOT STARTED)

After MIT9312 validates, proceed strain by strain:
MED4 → MIT9313 → NATL2A → MIT9301 → AS9601 → NATL1A → RSP50 → CC9311 → WH8102 → MIT1002 → EZ55 → HOT1A3

Each strain: snapshot → `build_gene_id_mapping --strains X --force` → `resolve_paper_ids --force` → `check-gene-ids` → rebuild KG → compare snapshot.

## Key Design Decisions (for reference)

- **Three-tier classification**: Tier 1 (gene-unique) → `specific_lookup` 1:1; Tier 2 (protein-level, paralogs OK) → `multi_lookup`; Tier 3 (generic names) → `multi_lookup`
- **Iterative convergence** (not Union-Find): all sources processed together in passes until stable; typically 2-3 passes; order-independent
- **Singleton rule for Tier 2+3**: `multi_lookup` entries used for resolution only when they map to exactly 1 gene in the organism
- **Never silent skip**: every row gets explicit `resolution_method` string in `_resolved.csv`
- **List expansion**: cells like `"PMM0001, PMM0002"` split on `,` and `;`; each part tried; full value tried first

## Files NOT changed (intentionally)

- `multiomics_kg/adapters/omics_adapter.py` — reads `_resolved.csv` unchanged
- `gene_mapping_supp.csv` — not yet deleted (some strains may still need legacy fallback until v2 rebuilt for all)
- `gene_id_mapping.json` files in cache — not yet regenerated (step 3 pending per strain)
