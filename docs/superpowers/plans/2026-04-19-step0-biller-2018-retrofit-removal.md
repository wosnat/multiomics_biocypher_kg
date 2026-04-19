# Step 0 — Biller 2018 Retrofit Removal + Baseline Capture

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Strip the 3 force-fit `gene_clusters` entries from Biller 2018's paperconfig, rebuild the KG cleanly, and commit a baseline (pre/post `/omics-edge-snapshot` + per-paper Cypher dump) so subsequent plans can assert deltas against a known state.

**Architecture:** No code changes. This is a paperconfig edit + cached-preprocessing cleanup + Neo4j rebuild + state capture. The retrofit creates a transient data-loss window (Tables S4A / S4B / S5 evidence is absent from the KG between this step and Plan 2).

**Tech Stack:** YAML, pandas-free bash, Docker Compose, cypher-shell, pytest, `/omics-edge-snapshot` skill.

**Spec reference:** `docs/superpowers/specs/2026-04-19-non-de-evidence-biller-2018-slice.md` §"Step 0 — Retrofit removal + baseline capture".

---

## File Structure

**Modify:**
- `data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml` — remove lines 180–254 (the 3 `gene_clusters` entries + the preceding comment header)

**Delete:**
- `data/Prochlorococcus/papers_and_supp/Biller 2018/table s4A sys003182233st4_resolved.csv`
- `data/Prochlorococcus/papers_and_supp/Biller 2018/table s4A sys003182233st4_resolved_report.txt`
- `data/Prochlorococcus/papers_and_supp/Biller 2018/table s4B sys003182233st4_resolved.csv`
- `data/Prochlorococcus/papers_and_supp/Biller 2018/table s4B sys003182233st4_resolved_report.txt`
- `data/Prochlorococcus/papers_and_supp/Biller 2018/table S5 sys003182233st5_resolved.csv`
- `data/Prochlorococcus/papers_and_supp/Biller 2018/table S5 sys003182233st5_resolved_report.txt`

**Create:**
- `docs/kg-changes/biller-2018-retrofit-baseline/before.json` — Biller-2018-specific Cypher dump before removal
- `docs/kg-changes/biller-2018-retrofit-baseline/after.json` — Biller-2018-specific Cypher dump after removal
- `docs/kg-changes/biller-2018-retrofit-baseline/delta.md` — human-readable summary of what changed
- `docs/kg-changes/biller-2018-retrofit-baseline/capture.cypher` — the reusable Cypher query script

**Skill outputs** (via `/omics-edge-snapshot`, auto-saved to `.claude/skills/omics-edge-snapshot/snapshots/`):
- `before_biller2018_retrofit_removal.json`
- `after_biller2018_retrofit_removal.json`

**Kept (read in subsequent plans):**
- `data/Prochlorococcus/papers_and_supp/Biller 2018/table s4A sys003182233st4.csv` (source; Plan 1 re-authors a `derived_metrics_table` entry pointing at this file)
- `data/Prochlorococcus/papers_and_supp/Biller 2018/table s4B sys003182233st4.csv`
- `data/Prochlorococcus/papers_and_supp/Biller 2018/table S5 sys003182233st5.csv`

**Not modified:**
- `tests/kg_validity/snapshot_data.json` — verified to contain no `Gene_in_gene_cluster` / `Clustering_analysis_has_gene_cluster` / `Publication_has_clustering_analysis` fixtures (grep'd at plan-writing time); no regeneration needed. If `pytest -m kg` surprises with a failure, Task 8 includes a contingency.

---

## Task Overview

| # | Task | Produces |
|---|---|---|
| 1 | Pre-flight verification | Clean git state + running Neo4j confirmed |
| 2 | Capture pre-removal `/omics-edge-snapshot` | `before_biller2018_retrofit_removal.json` |
| 3 | Capture pre-removal Biller 2018 Cypher dump | `docs/kg-changes/biller-2018-retrofit-baseline/before.json` + `capture.cypher` |
| 4 | Edit Biller 2018 paperconfig | 3 `gene_clusters` entries removed |
| 5 | Delete stale `_resolved.csv` files | 6 files gone |
| 6 | Re-run preprocessing (prepare_data steps 3–4) | Refreshed `gene_id_mapping.json` + `_resolved.csv` for remaining CSVs |
| 7 | Rebuild the KG | `docker compose up -d` reaches `deploy` clean |
| 8 | Run test suites | `pytest -m "not slow and not kg"` + `pytest -m kg` green |
| 9 | Capture post-removal `/omics-edge-snapshot` + Cypher dump | `after_biller2018_retrofit_removal.json` + `after.json` |
| 10 | Compare + write delta summary | `delta.md` documenting the drop |
| 11 | Commit | Paperconfig edit + baselines + delta summary in one commit |

---

## Task 1: Pre-flight verification

**Files:** (read-only)

- [ ] **Step 1: Confirm clean git state**

```bash
git status --short
```

Expected: no output (or only the paperconfig intentionally being edited in this plan). If there's unrelated dirty state, stop and ask the user whether to stash.

- [ ] **Step 2: Confirm Neo4j is running**

```bash
docker compose ps
```

Expected output includes a `deploy` (or `neo4j`) service in state `Up`. If not, run:

```bash
docker compose up -d
```

and wait until `curl -s http://localhost:7474` returns an HTTP 200-like response.

- [ ] **Step 3: Confirm the pre-flight KG has the 3 retrofitted ClusteringAnalysis nodes** (anchors what we're removing)

Run via your Cypher preferred tool (`mcp__multiomics-kg__run_cypher` or cypher-shell). Query:

```cypher
MATCH (ca:ClusteringAnalysis)
WHERE ca.id CONTAINS 'natl2a_periodicity'
   OR ca.id CONTAINS 'mit1002_periodicity'
   OR ca.id CONTAINS 'natl2a_darkness_survival'
RETURN ca.id, ca.name, ca.total_gene_count
ORDER BY ca.id;
```

Expected: exactly 3 rows. If 0 rows, the retrofit has already been stripped from the deployed KG — confirm with the user before proceeding (step 0 may already be partially done). If more than 3, scope has drifted; investigate.

- [ ] **Step 4: Record the Biller 2018 DOI** (for subsequent Cypher queries)

```cypher
MATCH (p:Publication)
WHERE p.title CONTAINS 'heterotroph interactions'
   OR p.title CONTAINS 'Prochlorococcus transcriptome'
RETURN p.id, p.doi, p.title
LIMIT 5;
```

Note the `p.id` value (expected form: `doi:10.1038/s41396-...`). Store it as `$BILLER_PUB_ID` in your notes — every subsequent Cypher query in this plan uses this literal.

If 0 rows are returned, the Publication node may be indexed by `name` instead. Fall back:

```cypher
MATCH (p:Publication)-[:Has_experiment]->(e:Experiment)
WHERE e.id CONTAINS 'darkness_extended_darkness_natl2a_rnaseq'
RETURN DISTINCT p.id, p.doi, p.title
LIMIT 5;
```

Still 0? Ask the user — the Publication node may not exist, indicating a deeper issue outside this step's scope.

---

## Task 2: Capture pre-removal `/omics-edge-snapshot`

**Files:**
- Create (via skill): `.claude/skills/omics-edge-snapshot/snapshots/before_biller2018_retrofit_removal.json`

**Why:** Establishes the reference `Changes_expression_of` edge count per Publication. Since step 0 doesn't touch DE edges, Biller 2018's `Changes_expression_of` count must be identical in the post-removal snapshot.

- [ ] **Step 1: Run the snapshot tool with --save**

```bash
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py \
  --save before_biller2018_retrofit_removal
```

Expected output: a summary table with per-publication `Changes_expression_of` counts, ending with "Saved to …before_biller2018_retrofit_removal.json".

- [ ] **Step 2: Verify file exists and contains Biller 2018 counts**

```bash
ls -l .claude/skills/omics-edge-snapshot/snapshots/before_biller2018_retrofit_removal.json
grep '"doi:10.1038/s41396' .claude/skills/omics-edge-snapshot/snapshots/before_biller2018_retrofit_removal.json || echo "MISSING BILLER 2018 IN SNAPSHOT"
```

Expected: file listed with non-zero size; grep prints one or more matching lines with Biller 2018's doi. If the grep falls back to "MISSING BILLER 2018 IN SNAPSHOT", the snapshot's DOI prefix is different — inspect the JSON manually to find Biller 2018's entry. Biller 2018 DE edges exist (Tables S3 and S6B), so it MUST appear in the snapshot.

---

## Task 3: Capture pre-removal Biller 2018 Cypher dump

**Files:**
- Create: `docs/kg-changes/biller-2018-retrofit-baseline/capture.cypher` (reusable query script)
- Create: `docs/kg-changes/biller-2018-retrofit-baseline/before.json` (output of the query against the current KG)

**Why:** `/omics-edge-snapshot` only captures `Changes_expression_of` counts. The baseline also needs ClusteringAnalysis node counts + `Gene_in_gene_cluster` edge counts + `Publication_has_clustering_analysis` / `Clustering_analysis_has_gene_cluster` / `Experiment_has_clustering_analysis` counts scoped to Biller 2018 — so Plan 2 can assert "the retrofit removed exactly these cluster edges and replaced them with DerivedMetric edges."

- [ ] **Step 1: Create the Cypher capture script**

Create `docs/kg-changes/biller-2018-retrofit-baseline/capture.cypher` (the `:param` approach lets one script serve both pre- and post-runs):

```cypher
// Biller 2018 state dump. Parameter: $pub_id = "doi:10.1038/s41396-..."
// (see Task 1 Step 4 for how to obtain the value).
:param pub_id => "doi:10.1038/s41396-018-0129-6";

// --- Publication presence ---
MATCH (p:Publication {id: $pub_id})
RETURN 'publication_exists' AS metric, count(p) AS value
UNION ALL
// --- Experiments ---
MATCH (p:Publication {id: $pub_id})-[:Has_experiment]->(e:Experiment)
RETURN 'experiment_count' AS metric, count(e) AS value
UNION ALL
// --- Changes_expression_of edges (must be unchanged after step 0) ---
MATCH (p:Publication {id: $pub_id})-[:Has_experiment]->(e:Experiment)-[r:Changes_expression_of]->(:Gene)
RETURN 'changes_expression_of_edges' AS metric, count(r) AS value
UNION ALL
// --- ClusteringAnalysis nodes (the 3 retrofitted + any non-retrofitted) ---
MATCH (p:Publication {id: $pub_id})-[:Publication_has_clustering_analysis]->(ca:ClusteringAnalysis)
RETURN 'clustering_analysis_count' AS metric, count(ca) AS value
UNION ALL
// --- Retrofitted ClusteringAnalysis presence ---
MATCH (ca:ClusteringAnalysis)
WHERE ca.id CONTAINS 'natl2a_periodicity'
   OR ca.id CONTAINS 'mit1002_periodicity'
   OR ca.id CONTAINS 'natl2a_darkness_survival'
RETURN 'retrofitted_clustering_analysis_count' AS metric, count(ca) AS value
UNION ALL
// --- GeneCluster nodes from retrofitted analyses ---
MATCH (ca:ClusteringAnalysis)-[:Clustering_analysis_has_gene_cluster]->(gc:GeneCluster)
WHERE ca.id CONTAINS 'natl2a_periodicity'
   OR ca.id CONTAINS 'mit1002_periodicity'
   OR ca.id CONTAINS 'natl2a_darkness_survival'
RETURN 'retrofitted_gene_cluster_count' AS metric, count(gc) AS value
UNION ALL
// --- Gene_in_gene_cluster edges from retrofitted analyses ---
MATCH (ca:ClusteringAnalysis)-[:Clustering_analysis_has_gene_cluster]->(gc:GeneCluster)-[r:Gene_in_gene_cluster]->(:Gene)
WHERE ca.id CONTAINS 'natl2a_periodicity'
   OR ca.id CONTAINS 'mit1002_periodicity'
   OR ca.id CONTAINS 'natl2a_darkness_survival'
RETURN 'retrofitted_gene_in_gene_cluster_edges' AS metric, count(r) AS value
UNION ALL
// --- Experiment_has_clustering_analysis edges pointing at retrofitted ---
MATCH (:Experiment)-[r:Experiment_has_clustering_analysis]->(ca:ClusteringAnalysis)
WHERE ca.id CONTAINS 'natl2a_periodicity'
   OR ca.id CONTAINS 'mit1002_periodicity'
   OR ca.id CONTAINS 'natl2a_darkness_survival'
RETURN 'retrofitted_experiment_has_ca_edges' AS metric, count(r) AS value
UNION ALL
// --- Per-experiment clustering_analysis_count property (if set by post-import) ---
MATCH (p:Publication {id: $pub_id})-[:Has_experiment]->(e:Experiment)
RETURN 'experiment_clustering_analysis_counts' AS metric, collect([e.id, e.clustering_analysis_count]) AS value;
```

If `:param` isn't supported in your cypher-shell version, substitute `$pub_id` with the literal DOI string (from Task 1 Step 4) inline, save as `capture.cypher`, and proceed.

- [ ] **Step 2: Run the capture against the live KG**

Make the output directory first:

```bash
mkdir -p docs/kg-changes/biller-2018-retrofit-baseline
```

Run (substitute `$BILLER_PUB_ID` from Task 1 Step 4):

```bash
cypher-shell -a bolt://localhost:7687 --format plain \
  -P "pub_id => \"$BILLER_PUB_ID\"" \
  < docs/kg-changes/biller-2018-retrofit-baseline/capture.cypher \
  > docs/kg-changes/biller-2018-retrofit-baseline/before.raw.txt
```

If cypher-shell isn't on your host, use the mcp tool instead:

```
mcp__multiomics-kg__run_cypher with the script text, capture the JSON result array
```

- [ ] **Step 3: Convert to structured JSON**

Write a small Python one-liner — save output as valid JSON for consistent diffing later:

```bash
uv run python - <<'PYEOF' > docs/kg-changes/biller-2018-retrofit-baseline/before.json
import json, sys
from pathlib import Path

raw = Path("docs/kg-changes/biller-2018-retrofit-baseline/before.raw.txt").read_text()
# Parse the plain-format cypher-shell output: header row + data rows, pipe-separated
rows = [line.split("|") for line in raw.strip().splitlines() if "|" in line]
header, *data = rows
records = [dict(zip([h.strip() for h in header], [c.strip() for c in row])) for row in data]
json.dump({"snapshot": "before_biller2018_retrofit_removal", "records": records}, sys.stdout, indent=2)
PYEOF
```

(If you used the mcp tool instead of cypher-shell, just pretty-print the JSON response directly into `before.json`.)

- [ ] **Step 4: Verify `before.json` is well-formed**

```bash
uv run python -c "import json; d = json.load(open('docs/kg-changes/biller-2018-retrofit-baseline/before.json')); print(len(d['records']), 'records')"
```

Expected: `8 records` (or `9 records` depending on how `experiment_clustering_analysis_counts` is serialized — the point is that the file parses as JSON and has the expected metrics).

Sanity assertions on the pre-removal state (inspect manually from `before.json`):
- `publication_exists` = 1
- `retrofitted_clustering_analysis_count` = 3
- `retrofitted_gene_in_gene_cluster_edges` > 0 (expect thousands — roughly 1800 + 530 + 270 per the parent spec's backlog table)
- `changes_expression_of_edges` > 0

If any sanity check fails, stop and ask the user — the KG state differs from what this plan assumes.

---

## Task 4: Remove the 3 `gene_clusters` entries from Biller 2018's paperconfig

**Files:**
- Modify: `data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml:180-254`

**Why:** The spec's retrofit step — these are the force-fit cluster entries that will be re-expressed as `derived_metrics_table` in Plan 1.

- [ ] **Step 1: Confirm the exact lines to remove**

```bash
wc -l "data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml"
sed -n '177,185p;250,254p' "data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml"
```

Expected at plan-writing time: file has 254 lines total.
- Line 177: `        pvalue_threshold: 0.1`
- Line 178: `        prefiltered: true` (last field of the preceding DE entry, `s6b_coculture_mit1002_5hr`)
- Line 179: blank
- Line 180: `    # ── Gene cluster data: periodicity classifications ───────────────────`
- Line 181: `    natl2a_periodicity:` (start of the first retrofit entry)
- Line 254: `        - darkness_extended_darkness_natl2a_rnaseq_coculture` (last field of `natl2a_darkness_survival`, also the last line of the file)

If the line numbers don't match (someone edited the file between plan-writing and execution), recompute the deletion range by finding each entry's key/end.

- [ ] **Step 2: Remove the block**

Use the Edit tool to replace the block with a clean empty span. Old string:

```yaml

    # ── Gene cluster data: periodicity classifications ───────────────────
    natl2a_periodicity:
```

Show the full block starting from the blank line before the comment through the last line of `natl2a_darkness_survival`. Use the Edit tool with a single replacement (old_string → empty).

A safer alternative if Edit struggles with the long block: use `sed` to excise lines 179-254 inclusive:

```bash
sed -i '179,254d' "data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml"
```

- [ ] **Step 3: Verify the block is gone**

```bash
grep -n 'natl2a_periodicity\|mit1002_periodicity\|natl2a_darkness_survival\|periodicity_cluster\|darkness_cluster' "data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml" || echo "CLEAN"
```

Expected: `CLEAN` (no matches).

- [ ] **Step 4: Verify the YAML still parses**

```bash
uv run python -c "import yaml; d = yaml.safe_load(open('data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml')); print('OK, experiments:', list(d['publication']['experiments'].keys()))"
```

Expected: `OK, experiments: [...]` with the expected experiment keys (`darkness_extended_darkness_natl2a_rnaseq_axenic`, `darkness_extended_darkness_natl2a_rnaseq_coculture`, `darkness_extended_darkness_mit1002_rnaseq`, etc.). If this errors, the YAML got corrupted during the edit — restore from git and retry.

- [ ] **Step 5: Verify against the existing validator**

```bash
uv run python scripts/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml"
```

Expected: exits 0 with no errors. If errors appear relating to the removed entries still being referenced, another part of the paperconfig points at them — investigate and widen the removal scope.

---

## Task 5: Delete stale `_resolved.csv` files

**Files:**
- Delete (6 files): all `*_resolved.csv` and `*_resolved_report.txt` under `data/Prochlorococcus/papers_and_supp/Biller 2018/` corresponding to Tables S4A / S4B / S5

**Why:** The resolved files were generated by `resolve_paper_ids.py` against the now-removed paperconfig entries. Leaving them would be preprocessing drift — future `--force` runs would skip them if mtime is fresh, and a human inspecting the directory would be confused.

- [ ] **Step 1: List the files to be deleted**

```bash
ls -l "data/Prochlorococcus/papers_and_supp/Biller 2018/"*_resolved*
```

Expected: 6 matching files (3 `_resolved.csv` + 3 `_resolved_report.txt`) for S4A, S4B, S5. Plus possibly `_resolved.csv` / `_resolved_report.txt` for tables NOT being retrofitted — if more than 6 files match, inspect individually; only delete the 6 for the 3 retrofitted tables.

- [ ] **Step 2: Delete the 6 files**

```bash
rm "data/Prochlorococcus/papers_and_supp/Biller 2018/table s4A sys003182233st4_resolved.csv"
rm "data/Prochlorococcus/papers_and_supp/Biller 2018/table s4A sys003182233st4_resolved_report.txt"
rm "data/Prochlorococcus/papers_and_supp/Biller 2018/table s4B sys003182233st4_resolved.csv"
rm "data/Prochlorococcus/papers_and_supp/Biller 2018/table s4B sys003182233st4_resolved_report.txt"
rm "data/Prochlorococcus/papers_and_supp/Biller 2018/table S5 sys003182233st5_resolved.csv"
rm "data/Prochlorococcus/papers_and_supp/Biller 2018/table S5 sys003182233st5_resolved_report.txt"
```

- [ ] **Step 3: Verify deletion**

```bash
ls "data/Prochlorococcus/papers_and_supp/Biller 2018/"*S4A*_resolved* 2>&1 || echo "S4A resolved gone"
ls "data/Prochlorococcus/papers_and_supp/Biller 2018/"*s4A*_resolved* 2>&1 || echo "s4A resolved gone"
ls "data/Prochlorococcus/papers_and_supp/Biller 2018/"*S5*_resolved* 2>&1 || echo "S5 resolved gone"
```

Expected: each `ls` reports "No such file or directory" and the fallback message fires. (Note the case-sensitive S4A vs s4A — the actual filename uses lowercase `s4A` / `s4B`.)

The source CSVs themselves (`table s4A sys003182233st4.csv`, etc.) must still exist — Plan 1 re-references them from new `derived_metrics_table` entries:

```bash
ls "data/Prochlorococcus/papers_and_supp/Biller 2018/table s4A sys003182233st4.csv"
ls "data/Prochlorococcus/papers_and_supp/Biller 2018/table s4B sys003182233st4.csv"
ls "data/Prochlorococcus/papers_and_supp/Biller 2018/table S5 sys003182233st5.csv"
```

Expected: all 3 listed.

---

## Task 6: Re-run preprocessing for NATL2A + MIT1002

**Files:**
- May modify: `cache/data/Prochlorococcus/genomes/NATL2A/gene_id_mapping.json` (rebuilt with S4A/S4B/S5 id_columns no longer contributing)
- May modify: `cache/data/Prochlorococcus/genomes/NATL2A/gene_mapping_supp.csv` (backward-compat view)
- May modify: `cache/data/Alteromonas/genomes/MIT1002/gene_id_mapping.json`

**Why:** `build_gene_id_mapping.py` harvests `id_columns` from paperconfig entries. The 3 removed entries contributed NATL2A `locus_tag_cyanorak` / `old_locus_tag` / `gene_name` pairs and MIT1002 `locus_tag` pairs. Refreshing with `--force` rebuilds the mapping without those contributions. Any IDs that were *only* sourced from the retrofitted tables (unlikely but possible) drop out.

- [ ] **Step 1: Run prepare_data steps 3–4 with --force for the two strains**

```bash
bash scripts/prepare_data.sh --steps 3 4 --strains NATL2A MIT1002 --force
```

Expected duration: ~1-5 minutes. Output logs go to `logs/prepare_data_step3.log` and `logs/prepare_data_step4.log`.

- [ ] **Step 2: Verify clean completion**

```bash
tail -5 logs/prepare_data_step3.log
tail -5 logs/prepare_data_step4.log
```

Expected: both logs end with a success line (look for "Done", "Wrote", or a summary). No tracebacks. If either log shows a traceback, investigate (likely a referenced CSV path is wrong).

- [ ] **Step 3: Verify the gene_id_mapping.json files were rewritten**

```bash
stat -c '%y %n' cache/data/Prochlorococcus/genomes/NATL2A/gene_id_mapping.json cache/data/Alteromonas/genomes/MIT1002/gene_id_mapping.json
```

Expected: both files have mtime within the last few minutes.

---

## Task 7: Rebuild the KG end-to-end

**Files:** (no edits — a Docker lifecycle operation)

**Why:** The `import` container needs a clean Neo4j data volume to re-ingest the updated CSVs without conflicting with the previous state. We also want the `post-process` step to recompute rollups (e.g., `Publication.clustering_analysis_count`) against the new cluster set.

- [ ] **Step 1: Stop the deployed graph (releases the volume lock)**

```bash
docker compose down
```

Expected: all services stopped. If `deploy` / `app` are shown to have taken several seconds to release, that's normal.

- [ ] **Step 2: Start the full pipeline**

```bash
docker compose up -d
```

This triggers: `build` (runs `create_knowledge_graph.py`) → `import` (runs `neo4j-admin import`) → `post-process` → `deploy` → `app`.

Watch the build progress:

```bash
docker compose logs -f build
```

(Ctrl-C once the `build` service finishes — it will exit cleanly.)

- [ ] **Step 3: Wait for import to complete**

```bash
docker compose logs -f import
```

Expected: final line is something like `Imported X nodes and Y relationships`. Ctrl-C when the container exits.

- [ ] **Step 4: Verify the import report is clean**

```bash
cat output/import.status
grep -c 'Skipped' output/import.report || echo "No Skipped lines — clean"
```

Expected: `import.status` contains `0` (exit code). `import.report` has zero "Skipped" lines — dangling-relationship count is 0. If either check fails, stop: a broken import means the paperconfig edit had unintended fallout.

- [ ] **Step 5: Wait for post-process to complete**

```bash
docker compose logs -f post-process
```

Expected: final line is a `[timing]` summary. Ctrl-C when the container exits.

- [ ] **Step 6: Wait for deploy + app**

```bash
docker compose ps
```

Expected: `deploy` in state `Up`, `app` either `Up` or finished its startup. Neo4j is now live.

Confirm Neo4j accepts connections:

```bash
curl -s http://localhost:7474 > /dev/null && echo "Neo4j UP"
```

Expected: `Neo4j UP`.

---

## Task 8: Run the test suites

**Files:** (no edits unless snapshot regeneration is needed)

**Why:** Spec's definition-of-done requires both suites green. The snapshot test suite does NOT currently depend on retrofitted fixtures (verified at plan-writing time via `grep -c 'natl2a_periodicity\|…' snapshot_data.json` → 0), so it should pass without regeneration.

- [ ] **Step 1: Run the non-KG test suite**

```bash
uv run pytest -m "not slow and not kg" -q
```

Expected: all pass. If any fail, investigate — Task 4's YAML edit or Task 6's preprocessing may have broken something. Fix before proceeding.

- [ ] **Step 2: Run the KG validity suite**

```bash
uv run pytest -m kg -v
```

Expected: all pass. If `tests/kg_validity/test_snapshot.py` fails with missing nodes/edges, the snapshot's sample unexpectedly included a retrofitted cluster node — proceed to Step 3 contingency. If other tests fail (e.g., cluster-count assertions elsewhere), they may encode the old state — investigate per-test.

- [ ] **Step 3 (contingency — only if Step 2 snapshot-test failed): Regenerate the snapshot**

```bash
uv run python tests/kg_validity/generate_snapshot.py
```

Expected: writes a new `tests/kg_validity/snapshot_data.json`. Re-run `pytest -m kg` — should now pass. Add the regenerated file to Task 11's commit.

---

## Task 9: Capture post-removal baselines

**Files:**
- Create (via skill): `.claude/skills/omics-edge-snapshot/snapshots/after_biller2018_retrofit_removal.json`
- Create: `docs/kg-changes/biller-2018-retrofit-baseline/after.json`

- [ ] **Step 1: Run `/omics-edge-snapshot --save after_biller2018_retrofit_removal`**

```bash
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py \
  --save after_biller2018_retrofit_removal
```

Expected: summary table + "Saved to …after_biller2018_retrofit_removal.json".

- [ ] **Step 2: Run the Biller 2018 Cypher dump again**

Use the same `capture.cypher` from Task 3 (with the same `$BILLER_PUB_ID`):

```bash
cypher-shell -a bolt://localhost:7687 --format plain \
  -P "pub_id => \"$BILLER_PUB_ID\"" \
  < docs/kg-changes/biller-2018-retrofit-baseline/capture.cypher \
  > docs/kg-changes/biller-2018-retrofit-baseline/after.raw.txt
```

- [ ] **Step 3: Convert to JSON**

```bash
uv run python - <<'PYEOF' > docs/kg-changes/biller-2018-retrofit-baseline/after.json
import json, sys
from pathlib import Path

raw = Path("docs/kg-changes/biller-2018-retrofit-baseline/after.raw.txt").read_text()
rows = [line.split("|") for line in raw.strip().splitlines() if "|" in line]
header, *data = rows
records = [dict(zip([h.strip() for h in header], [c.strip() for c in row])) for row in data]
json.dump({"snapshot": "after_biller2018_retrofit_removal", "records": records}, sys.stdout, indent=2)
PYEOF
```

- [ ] **Step 4: Verify the post-removal state**

```bash
cat docs/kg-changes/biller-2018-retrofit-baseline/after.json
```

Expected values:
- `publication_exists` = 1 (unchanged)
- `experiment_count` = unchanged vs `before.json`
- `changes_expression_of_edges` = **identical** to `before.json` (DE path untouched)
- `clustering_analysis_count` = `before.json` value **minus 3** (exactly 3 CAs dropped)
- `retrofitted_clustering_analysis_count` = 0
- `retrofitted_gene_cluster_count` = 0
- `retrofitted_gene_in_gene_cluster_edges` = 0
- `retrofitted_experiment_has_ca_edges` = 0

If any of these don't match exactly — stop and investigate. The most likely causes: a stale Neo4j volume (run Task 7 again from `docker compose down`), or the DOI recorded in Task 1 Step 4 was wrong (rerun that step and redo Tasks 3 and 9).

- [ ] **Step 5: Delete the intermediate `.raw.txt` files** (keep only the `.json` artifacts)

```bash
rm docs/kg-changes/biller-2018-retrofit-baseline/before.raw.txt
rm docs/kg-changes/biller-2018-retrofit-baseline/after.raw.txt
```

---

## Task 10: Write the delta summary

**Files:**
- Create: `docs/kg-changes/biller-2018-retrofit-baseline/delta.md`

**Why:** Human-readable summary of what this step changed in the KG. Future plans (especially Plan 2, which restores evidence through a different shape) cite this file to explain edge-count deltas vs step-0 baseline.

- [ ] **Step 1: Run the comparison snapshot tool**

```bash
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py \
  --compare before_biller2018_retrofit_removal \
  --against after_biller2018_retrofit_removal
```

Expected output: a diff table per-publication; Biller 2018's `Changes_expression_of` count should be unchanged, and no per-paper regressions for other publications.

Copy the output (or key lines from it) — it's an input to the summary below.

- [ ] **Step 2: Write the delta summary**

Create `docs/kg-changes/biller-2018-retrofit-baseline/delta.md` with the following structure. Fill in the exact numbers from `before.json` and `after.json`:

```markdown
# Biller 2018 Retrofit Removal — Step 0 Baseline

**Date:** YYYY-MM-DD (today)
**Plan:** `docs/superpowers/plans/2026-04-19-step0-biller-2018-retrofit-removal.md`
**Spec:** `docs/superpowers/specs/2026-04-19-non-de-evidence-biller-2018-slice.md`

## What changed

Removed 3 `gene_clusters` entries from Biller 2018's paperconfig and rebuilt the KG:
- `natl2a_periodicity` — sourced from Table S4A
- `mit1002_periodicity` — sourced from Table S4B
- `natl2a_darkness_survival` — sourced from Table S5

Plan 2 re-expresses this evidence through `derived_metrics_table` entries emitting `derived_metric_flags_gene` (periodicity) + `derived_metric_classifies_gene` (darkness survival) edges. Between step 0 and Plan 2, this evidence is absent from the KG — Cypher queries referencing the retrofitted cluster node IDs will return 0 rows.

## Biller 2018 per-paper state

| Metric | Before | After | Delta |
|---|---|---|---|
| `Changes_expression_of` edges | FILL | FILL | **0 (unchanged)** |
| `ClusteringAnalysis` nodes | FILL | FILL | -3 |
| Retrofitted `ClusteringAnalysis` | 3 | 0 | -3 |
| Retrofitted `GeneCluster` nodes | FILL | 0 | FILL |
| Retrofitted `Gene_in_gene_cluster` edges | FILL | 0 | FILL |
| Retrofitted `Experiment_has_clustering_analysis` edges | FILL | 0 | FILL |

## `/omics-edge-snapshot` delta (cross-paper)

Biller 2018's `Changes_expression_of` count must be unchanged; all other publications' counts must also be unchanged (DE path untouched project-wide). Full diff output:

```
(paste Task 10 Step 1 output here — the --compare table)
```

## Artifacts committed alongside this document

- `before.json` — Biller-2018-specific Cypher dump before removal
- `after.json` — Biller-2018-specific Cypher dump after removal
- `capture.cypher` — the Cypher script used for both dumps (reusable)
- `.claude/skills/omics-edge-snapshot/snapshots/before_biller2018_retrofit_removal.json` — cross-paper DE snapshot
- `.claude/skills/omics-edge-snapshot/snapshots/after_biller2018_retrofit_removal.json` — cross-paper DE snapshot

## Follow-up

Plan 1 (vocab + schema + paperconfig preprocessing + Biller 2018 paperconfig authoring) adds the 3 `derived_metrics_table` entries that re-cover this evidence in the new shape; Plan 2 wires the adapter so the KG emits `derived_metric_flags_gene` + `derived_metric_classifies_gene` edges from those entries; Plan 3 adds post-import rollups and KG validity assertions.
```

Fill the `FILL` cells with the integer values from your `before.json` and `after.json`. `Delta = after - before` (signed).

- [ ] **Step 3: Verify the delta.md is valid markdown**

```bash
head -40 docs/kg-changes/biller-2018-retrofit-baseline/delta.md
```

Expected: rendered clean, table shows correct numbers, the fenced code block has the actual `--compare` output (not a placeholder).

---

## Task 11: Commit

- [ ] **Step 1: Review everything that will be committed**

```bash
git status --short
git diff --stat
```

Expected to see:
- `M data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml`
- `D data/Prochlorococcus/papers_and_supp/Biller 2018/table s4A sys003182233st4_resolved.csv`
- `D data/Prochlorococcus/papers_and_supp/Biller 2018/table s4A sys003182233st4_resolved_report.txt`
- `D data/Prochlorococcus/papers_and_supp/Biller 2018/table s4B sys003182233st4_resolved.csv`
- `D data/Prochlorococcus/papers_and_supp/Biller 2018/table s4B sys003182233st4_resolved_report.txt`
- `D data/Prochlorococcus/papers_and_supp/Biller 2018/table S5 sys003182233st5_resolved.csv`
- `D data/Prochlorococcus/papers_and_supp/Biller 2018/table S5 sys003182233st5_resolved_report.txt`
- `?? .claude/skills/omics-edge-snapshot/snapshots/before_biller2018_retrofit_removal.json`
- `?? .claude/skills/omics-edge-snapshot/snapshots/after_biller2018_retrofit_removal.json`
- `?? docs/kg-changes/biller-2018-retrofit-baseline/before.json`
- `?? docs/kg-changes/biller-2018-retrofit-baseline/after.json`
- `?? docs/kg-changes/biller-2018-retrofit-baseline/capture.cypher`
- `?? docs/kg-changes/biller-2018-retrofit-baseline/delta.md`
- (contingency) `M tests/kg_validity/snapshot_data.json` — only if Task 8 Step 3 fired

Also refreshed (keep uncommitted OR include — decide with the user; they're cache artifacts):
- `M cache/data/Prochlorococcus/genomes/NATL2A/gene_id_mapping.json` (or similar)
- `M cache/data/Alteromonas/genomes/MIT1002/gene_id_mapping.json`

Cache files under `cache/` are typically gitignored — check `git check-ignore cache/data/Prochlorococcus/genomes/NATL2A/gene_id_mapping.json` and only stage them if they're tracked.

- [ ] **Step 2: Stage the intentional changes**

```bash
git add "data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml"
git add -u "data/Prochlorococcus/papers_and_supp/Biller 2018/"    # stages the 6 deletions
git add .claude/skills/omics-edge-snapshot/snapshots/before_biller2018_retrofit_removal.json
git add .claude/skills/omics-edge-snapshot/snapshots/after_biller2018_retrofit_removal.json
git add docs/kg-changes/biller-2018-retrofit-baseline/
# Contingency — only if Task 8 Step 3 fired:
# git add tests/kg_validity/snapshot_data.json
```

Verify:

```bash
git status --short
```

Only the intended additions/deletions/modifications should be staged.

- [ ] **Step 3: Commit**

```bash
git commit -m "$(cat <<'EOF'
step 0: strip Biller 2018 force-fit clusters, capture retrofit baseline

Removes 3 gene_clusters entries from Biller 2018's paperconfig
(natl2a_periodicity, mit1002_periodicity, natl2a_darkness_survival).
Deletes the stale _resolved.csv files for the 3 source tables (S4A/S4B/S5).
Re-runs prepare_data steps 3-4 for NATL2A + MIT1002 and rebuilds the KG.

Captures a per-paper Cypher dump + an omics-edge-snapshot as the baseline
for the Biller 2018 slice. Plan 1 re-authors the 3 tables as
derived_metrics_table entries; Plan 2 adapters consume them; Plan 3
verifies rollups. Between this commit and Plan 2, the 3 tables'
evidence is absent from the KG — intentional transient state.

Spec: docs/superpowers/specs/2026-04-19-non-de-evidence-biller-2018-slice.md

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 4: Verify commit landed**

```bash
git log --oneline -3
git show --stat HEAD
```

Expected: HEAD is the new commit; `--stat` shows the paperconfig edit + 6 deletions + new baseline files.

---

## Definition of done

All of the following must be green before handing off to Plan 1:

- [ ] `git status --short` is clean (or contains only intentional cache-file changes under `cache/` if those are tracked).
- [ ] `uv run pytest -m "not slow and not kg" -q` passes.
- [ ] `uv run pytest -m kg -v` passes (with snapshot regen if it fired in Task 8 Step 3).
- [ ] `docs/kg-changes/biller-2018-retrofit-baseline/after.json` shows `changes_expression_of_edges` identical to `before.json`, `retrofitted_*` counts all 0, and `clustering_analysis_count` reduced by exactly 3.
- [ ] `.claude/skills/omics-edge-snapshot/snapshots/` contains both `before_` and `after_biller2018_retrofit_removal.json`, and their cross-paper comparison shows no unintended per-paper regressions.
- [ ] `docs/kg-changes/biller-2018-retrofit-baseline/delta.md` has no `FILL` placeholders remaining.
- [ ] The commit message references the spec and names the 3 removed entries.

## Handoff to Plan 1

Plan 1 (`docs/superpowers/plans/YYYY-MM-DD-plan1-…`, to be written after step 0 lands) picks up from this committed state and adds:
- `multiomics_kg/vocab/non_de_evidence.py` vocabulary module
- `DerivedMetric` + edges in `config/schema_config.yaml`
- `iter_derived_metrics_tables` + validator + `build_gene_id_mapping` + `resolve_paper_ids` extensions
- **Biller 2018 paperconfig authoring** — adds 3 new `derived_metrics_table` entries to the file this plan emptied
- Adapter + KG build wire-in is deferred to Plan 2.
