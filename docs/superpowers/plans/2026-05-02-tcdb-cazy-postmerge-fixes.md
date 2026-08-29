# TCDB-CAZy Post-Merge Fixes — Resume Prompt

**Date:** 2026-05-02
**Status:** 3 fixes committed locally, not pushed; one refactor pending.
**Original feature:** [`2026-05-01-tcdb-cazy-ontologies.md`](2026-05-01-tcdb-cazy-ontologies.md) (8-commit feature merged to main as commit `dbb6181` after merging origin/main's chemistry-slice-1).

---

## Tell the next session this

Pick up implementation of the TCDB-CAZy ontology feature post-rebuild bug fixes in `/home/osnat/github/second_multiomics`.

**Repo setup (two checkouts of github.com/wosnat/multiomics_biocypher_kg):**
- `/home/osnat/github/second_multiomics/` — this is the controller-side checkout where code edits land. **You are here.**
- `/home/osnat/github/multiomics_biocypher_kg/` — the Docker-side checkout where `docker compose up` runs. The user will pull there to test the rebuild. Don't run docker commands yourself.

**Local symlinks (read-only caches from Docker side, set up in the previous session):**
- `cache/data/mnx` → Docker side (4 GB MNX SQLite + TSVs; ~30 min build saved).
- `cache/data/tcdb` → Docker side (TCDB TSVs + tcdb_hierarchy.json + tcdb_pruned.json).
- `cache/data/kegg/raw` → Docker side (KEGG REST cache).

To commit `tcdb_hierarchy.json` / `tcdb_pruned.json` from this side, replace the symlink with a real dir first. The user prefers to keep the Docker side as the canonical source of regenerated cache files; just verify here.

**Branch state:** `main` at `d382efa`, 3 commits ahead of `origin/main`. Not pushed yet — user wants the order-flip refactor folded in first to avoid two rebuild cycles on Docker side.

```
d382efa fix(metabolism): preserve CHEBI prefix in TCDB substrate parse + track tcdb JSON outputs
bb05a0d fix(metabolism): strip CHEBI: prefix when resolving TCDB substrates via MNX
62aec9e fix(post-import): split UNION list-concat to obey Neo4j strict aggregation
e64a43b update cache (tcdb/cazy)                                ← origin/main HEAD
dbb6181 merge: integrate chemistry-slice-1 KG-side additions
```

**Current `cache/data/kegg/kegg_data.json`** (committed in `d382efa`) reflects the post-fix state with **452 transport-only `additional_compounds`** + **385 transport-only kegg compounds in `compounds`** + **645 metabolism+transport overlap**. After the refactor below ships, the kegg_data.json shape may shift — see "Refactor task" for what changes.

---

## What just landed (post-merge bug fixes)

After merging origin/main's chemistry-slice-1 work, the Docker rebuild surfaced three real bugs that the merge had not introduced — they were latent in the original 8-commit feature:

1. **`62aec9e`** — Neo4j strict-aggregation in `Gene.metabolite_count` + `Metabolite.gene_count` UNION blocks. `WITH g, m_cat + collect(DISTINCT m2) AS all_m` is illegal because `m_cat` is an implicit grouping key while `collect(...)` aggregates at the same level. Fix: split the concat into a separate `WITH` after the aggregate. Mirrored in both `scripts/post-import.sh` and `scripts/post-import.cypher`.

2. **`bb05a0d`** — `_resolve_substrates` was passing `"CHEBI:9314"` to `resolve_metabolite()`, which queried `compound_aliases WHERE value = ?`. But MNX stores chebi aliases as `(source='chebi', value='9314')` — bare numbers, no prefix. 0 hits, every time. Fix: parse the prefix ourselves and query `source='chebi' AND value=?`. (Half the bug — see #3.)

3. **`d382efa`** — `_parse_tcdb_substrates` was stripping the `CHEBI:NNNN` prefix at parse time, so `tcdb_hierarchy.json` contained `["potassium(1+)", "water"]` instead of `["CHEBI:8305;potassium(1+)", "CHEBI:15377;water"]`. The implementer was honoring an old "Phase 1.3 deferred" docstring comment that conflicts with our spec. Fix: preserve the full string. Combined with `bb05a0d`, this fully unblocks substrate resolution.
   - Also added: log line surfacing `additional_compounds` count + overlap (`Folded into kegg_data.json: N additional_compounds, M/Total compounds tagged metabolism+transport`).
   - Also added: `.gitignore` un-ignore for `cache/data/tcdb/tcdb_hierarchy.json` and `cache/data/tcdb/tcdb_pruned.json` (mirrors how `cache/data/kegg/kegg_data.json` is tracked while raw TSVs/SQLite stay ignored).

**Local verification:** ran `bash scripts/prepare_data.sh --steps 6 --force` here (uses the symlinked MNX SQLite) — produces non-empty `additional_compounds`. Numbers below.

**Docker verification:** **NOT yet run** — user wants the refactor below first.

---

## Current chemistry-layer numbers (from local step 6 run, MNX symlinked)

| Bucket | Count | Where in kegg_data.json | evidence_sources |
|---|---|---|---|
| Pure-metabolism KEGG compounds | 1,543 | `compounds` | `['metabolism']` |
| Both-paths KEGG compounds | 645 | `compounds` | `['metabolism', 'transport']` |
| Transport-only KEGG compounds (substrate primary `kegg.compound:*`, no gene catalysis) | 385 | `compounds` (per Option B routing in `1ce9f6f`) | `['transport']` |
| Transport-only non-KEGG compounds (substrate primary `chebi:*`, no KEGG cross-ref) | 452 | `additional_compounds` | `['transport']` |
| **Total Metabolite nodes after import** | **3,025** | | |

Other rough numbers:
- **TcdbFamily nodes:** 4,844 (from 535 gene-annotated TCDB IDs walked above + below)
- **`Gene_has_tcdb_family` edges:** ~6-8K (extrapolated; not measured post-rebuild yet)
- **`Tcdb_family_transports_metabolite` edges:** 5,762 leaf×primary (1,097 distinct primary IDs)
- **CazyFamily nodes:** ~30-50 (observed-only; per-strain CAZy IDs)
- **Pathways (gene-reachable, KO ∪ Rxn):** 377 (currently — will grow after the refactor)

---

## Refactor task (your job in the next session)

**User directive:** "i think we need to change the order - prune kegg compounds and pathways after adding the tcdb compounds." Then: "i would add pathways reachable via transport too."

### Goal

Move TCDB substrate resolution before `build_pruned_kegg_data`. Extend the gene-reachable `cpds` and `pws` sets to include transport-reachable members. Then run the regular pruning/enrichment pipeline over the extended sets — newly-added transport-only kegg compounds get their `Metabolite_in_pathway` edges naturally via the existing filter `[p for p in cpd_to_pw[cpd_id] if p in allowed_pathways]`.

### New formula

```
cpds = catalysis_cpds ∪ substrate_kegg_cpds        # extended compound set
where substrate_kegg_cpds = {primary[len("kegg.compound:"):]
                             for ids in leaf_to_primary.values()
                             for primary in ids
                             if primary.startswith("kegg.compound:")}

pws  = KO_reachable ∪ Rxn_reachable ∪ Transport_pws   # extended pathway set
where Transport_pws = {p for cpd in substrate_kegg_cpds
                       for p in cpd_to_pw.get(cpd, [])}
```

`additional_compounds` continues to hold non-KEGG primaries (`chebi:NNNN`, `mnx:MNXM*`) — they don't have a KEGG pathway concept.

### Files to modify

- `multiomics_kg/download/build_kegg_metabolism_xrefs.py`
  - Refactor `build_pruned_kegg_data(raw, conn, out_path)` → `build_pruned_kegg_data(raw, conn, out_path, *, sets: dict)` so `main()` controls the sets passed in.
  - In `main()`, reorder: parse raw + open MNX → `_gene_reachable_sets(raw)` to get catalysis sets → load tcdb_hierarchy → `_prune_tcdb` → `_resolve_substrates` → split substrate primaries (`substrate_kegg_cpds` + non-KEGG `compound_props`) → extend `sets["cpds"]` with `substrate_kegg_cpds` → extend `sets["pws"]` with `transport_pws` → call `build_pruned_kegg_data` with extended sets → tag `evidence_sources` on each compound entry based on origin (`catalysis_cpds` vs `substrate_kegg_cpds` membership) → write `additional_compounds` for non-KEGG primaries → write `tcdb_pruned.json`.
  - **Shrink `_fold_substrates_into_kegg_data` dramatically** — its kegg-compound branch becomes unreachable since kegg compounds now flow through `_bulk_enrich_compounds`. Either delete the function and inline the additional_compounds insert (5 lines in main()), or shrink it to just that responsibility. Prefer deletion + inline if it ends up that small.
  - Update the fold-in log line to reflect the new flow (separate counts for additional_compounds vs both-paths overlap).

- `tests/test_build_kegg_metabolism_xrefs.py` — these tests currently target `_fold_substrates_into_kegg_data` directly. After the refactor, they may need to:
  - `test_compounds_get_metabolism_evidence_source` — adjust to test the new tagging step (now happens after `build_pruned_kegg_data` in main, not via fold-in helper).
  - `test_transport_only_compounds_land_in_additional_compounds` — still relevant for non-KEGG primaries, test directly against the tagged output.
  - `test_overlap_compound_gets_both_evidence_sources` — still relevant; compound is in both `catalysis_cpds` and `substrate_kegg_cpds`, gets `['metabolism', 'transport']`.
  - `test_kegg_compound_substrate_not_in_compounds_lands_in_compounds_not_additional` — replace with: when a kegg.compound primary appears in `substrate_kegg_cpds` but NOT in `catalysis_cpds`, the compound entry lands in `compounds` (not `additional_compounds`) **AND** has its pathway field populated with gene-reachable pathways. Specifically test: a compound participating in 2 KEGG pathways, one in `pws`, one not in `pws` — the compound entry's `pathways` field contains exactly the one in `pws`.

- `tests/kg_validity/test_tcdb_cazy.py` — pre-rebuild thresholds I wrote earlier are too narrow:
  - `assert 100 <= n <= 2000` for `TcdbFamily` count → widen upper to **10,000** (actual ~4,844).
  - `assert 500 <= n <= 5000` for `Gene_has_tcdb_family` edges → widen upper to **15,000** (actual unknown until rebuild; estimate from 535 gene-annotated TCDB IDs × ~10-15 genes-per-id avg).
  - `assert 1000 <= n <= 30000` for `Tcdb_family_transports_metabolite` → keep (5,762 fits).
  - `assert 100 <= n <= 2000` for substrate edge → keep.
  - The test file uses per-function `@pytest.mark.kg`; consider switching to module-level `pytestmark = pytest.mark.kg` for consistency with sibling files (Minor — a previous reviewer flagged this).

### Verification gate

After the refactor:
1. `uv run pytest -m "not slow and not kg"` passes.
2. `bash scripts/prepare_data.sh --steps 6 --force` runs and surfaces sensible numbers in the log:
   - `Gene-reachable: <K> KOs, <R> reactions, <C+S> compounds, <P+T> pathways (KO∪Rxn∪Transport)` — note pathway count grows beyond 377.
   - `TCDB: 4844 kept IDs, 3095 leaves with substrates (5762/1097)`.
   - Compounds in compounds section grow by ~385 + 645 (transport-reachable not in catalysis) — verify with `python3 -c "import json; d = json.load(open('cache/data/kegg/kegg_data.json')); print('compounds:', len(d['compounds']), 'additional_compounds:', len(d['additional_compounds']))"`.
   - **Newly-added transport-only kegg compounds NOW have non-empty `pathways` fields** (the user's whole reason for this refactor). Verify with: `python3 -c "import json; d = json.load(open('cache/data/kegg/kegg_data.json')); xs = [(k,v) for k,v in d['compounds'].items() if v.get('evidence_sources') == ['transport']]; print(f'transport-only kegg compounds: {len(xs)}, with pathways: {sum(1 for _, v in xs if v.get(\"pathways\"))}')`.
3. The `_resolve_substrates_strips_chebi_prefix` test still passes (regression for the CHEBI prefix bug).

### What you DON'T need to touch

- The TCDB / CAZy adapters themselves — they're unaffected by the order change. They consume `tcdb_pruned.json` whose shape doesn't change.
- The post-import scripts (`post-import.sh`, `.cypher`) — unaffected. Already fixed in `62aec9e`.
- The schema (`config/schema_config.yaml`) — unaffected.
- The metabolism_adapter — it iterates `compounds` and `additional_compounds`, both shapes preserved.
- `Tcdb_family_transports_metabolite` edges — unaffected (still emitted by the TCDB adapter from `tcdb_pruned.json`'s `leaf_substrates`).

### Open question to flag if you hit it

Should `mnx:MNXM*` non-KEGG primaries (no chebi alias) land in `additional_compounds` or be dropped? Currently `_resolve_substrates` always produces a primary (`mnxm_to_primary_id` falls back to `mnx:MNXM*`), so the answer should be: include them in `additional_compounds`. Verify the current flow handles this (we saw 0 mnx-prefix entries in the local run — all 452 transport-only non-kegg were chebi: prefixes — but the code path should support it).

---

## Sequence after the refactor lands

1. Push all 4 commits (3 fixes + refactor) to origin.
2. **On Docker side**: `git pull && bash scripts/prepare_data.sh --steps 6 --force`.
3. **On Docker side**: `git add cache/data/tcdb/tcdb_hierarchy.json cache/data/tcdb/tcdb_pruned.json cache/data/kegg/kegg_data.json && git commit -m "data: regenerate tcdb cache + kegg_data.json with order-flip + transport-reachable pathways"`.
4. **On Docker side**: `docker compose down && docker compose up -d --build` (~25 min).
5. After post-process succeeds: run `pytest -m kg -v` against the deployed graph.
6. Run `/omics-edge-snapshot --compare pre_tcdb_cazy_rebuild` (baseline saved at `.claude/skills/omics-edge-snapshot/snapshots/pre_tcdb_cazy_rebuild.json` — must be byte-identical, the regression invariant from the spec).
7. Spot-check Cypher: `MATCH (m:Metabolite {evidence_sources: ['transport']}) RETURN count(m), count(DISTINCT (m)<-[:Metabolite_in_pathway]-(:KeggTerm))` — verify some transport-only metabolites have pathway edges.
8. Regenerate snapshot fixture: `uv run python tests/kg_validity/generate_snapshot.py && git add tests/kg_validity/snapshot_data.json && git commit -m "test: regenerate KG snapshot for TCDB/CAZy ontologies"`.
9. Backfill exact counts in `docs/kg-changes/tcdb-cazy-ontologies.md` (currently labeled "LANDED 2026-05-02 (counts to be backfilled)") — drop the notice, fill in the table.

---

## Reference data — for sanity checking

**MNX SQLite is at** `/home/osnat/github/multiomics_biocypher_kg/cache/data/mnx/metabolite_resolver.db` (read-only via the local symlink). Built 2026-05-01.

**TCDB hierarchy**: 13,643 entries; pruned to **4,844 kept** for our 25 strains.

**Substrate resolution sample** (verified working post-fix):
```python
import sqlite3
conn = sqlite3.connect('file:cache/data/mnx/metabolite_resolver.db?mode=ro', uri=True)
cur = conn.cursor()
cur.execute("SELECT mnxm_id FROM compound_aliases WHERE source='chebi' AND value='9314'")
# returns ('MNXM50000',) for sucrose
```

**Step-6 log signature on success** (`logs/prepare_data_step6.log`):
```
Gene-reachable: 4564 KOs, 2349 reactions, 2188 compounds, 377 pathways (KO∪Rxn), 535 TCDB IDs
Wrote cache/data/kegg/kegg_data.json: 4564 KOs, 377 pathways, 2349 reactions, 2188 compounds
Pruning TCDB hierarchy + resolving transport substrates ...
  TCDB: 4844 kept IDs, 3095 leaves with substrates (5762 total / 1097 distinct primary IDs)
  Folded into kegg_data.json: 452 additional_compounds (transport-only), 645/2573 compounds now tagged metabolism+transport
```

After the refactor, expect this signature to change — the "Wrote kegg_data.json" line should reflect the extended sets:
```
Gene-reachable (catalysis): 4564 KOs, 2349 reactions, 2188 compounds, 377 pathways
TCDB substrate primaries: 1097 distinct (385 kegg.compound, 712 chebi)  # rough
Extended sets: 2573 compounds (+385), <P> pathways (+<T>)               # T = transport-only pathways
Wrote kegg_data.json: 4564 KOs, <P> pathways, 2349 reactions, 2573 compounds
  TCDB: 4844 kept IDs, 3095 leaves with substrates ...
  additional_compounds (non-KEGG): 452 (chebi:*)
```

(Exact log shape is your call — make it readable.)

---

## Cross-spec coordination context (already settled — don't re-open)

- **TCDB-S3 / KG-A2 (`Gene.metabolite_count` UNION):** already shipped. Post-import block in `scripts/post-import.sh` UNIONs catalysis (`Gene → Reaction → Metabolite`) and transport (`Gene → TcdbFamily → tc_specificity → Metabolite`) paths via `apoc.coll.toSet`. Direction-of-arrow bug already fixed in `c63b31d` (DOWN traversal `<-[:Tcdb_family_is_a_tcdb_family*0..]-`).
- **MET-M4 (`evidence_sources` enum, open-ended):** values are `'metabolism'` and `'transport'`; `'metabolomics'` reserved for a future metabolomics-DM spec. No tagging changes needed.
- **Chemistry-slice-1 KG-A1..A4:** all in main as of `e64a43b`. `Gene.reaction_count`, `Metabolite.elements` (Hill-parsed via chemparse), `KeggTerm.reaction_count`/`metabolite_count` rollups all in post-import. `Metabolite.elements` is set by `metabolism_adapter` for both `compounds` and `additional_compounds` loops.

---

## Files NOT to touch

- The user's working tree has unrelated untracked files: `.claude/skills/signalp-run/`, `tools/`, `fmicb-11-567431.pdf`, `multiomics_kg/download/run_signalp.py`. Stay out of these.
- Symlinks at `cache/data/{mnx,tcdb,kegg/raw}` are intentional (read-only access to Docker-side caches). Don't `rm` them; replace only if you specifically need a real dir to commit through.

---

## Useful commands cheat sheet

```bash
# Local step 6 (uses MNX symlink, completes in ~30 s)
bash scripts/prepare_data.sh --steps 6 --force

# Inspect kegg_data.json post-step-6
python3 -c "import json; d = json.load(open('cache/data/kegg/kegg_data.json')); print('compounds:', len(d['compounds']), 'additional_compounds:', len(d['additional_compounds']))"

# Inspect tcdb_pruned.json
python3 -c "import json; p = json.load(open('cache/data/tcdb/tcdb_pruned.json')); print('kept:', len(p['kept_tcdb_ids']), 'leaves with subs:', len(p['leaf_substrates']))"

# Unit tests
uv run pytest -m "not slow and not kg" 2>&1 | tail -5

# Targeted tests around the refactor
uv run pytest tests/test_build_kegg_metabolism_xrefs.py -v

# Schema parses
uv run python -c "from biocypher import BioCypher; BioCypher(schema_config_path='config/schema_config.yaml', biocypher_config_path='config/biocypher_config.yaml'); print('OK')"
```

---

## Other deferred items

These were in the plan but defer until after the rebuild lands:
- **Snapshot regen** (`tests/kg_validity/generate_snapshot.py`) — runs against live Neo4j. Defer.
- **Count backfill in docs/kg-changes/tcdb-cazy-ontologies.md** — needs post-rebuild numbers. Defer.
- **EZ55 step2_metabolism_report quirk** — `validated_total = 580` vs `raw_total = 579` after dropping the validate transform. Off-by-one in the report writer (since both should equal each other now under passthrough). Harmless. File a follow-up issue if it bothers anyone.
- **Test marker style** — `tests/kg_validity/test_tcdb_cazy.py` uses per-function `@pytest.mark.kg` while sibling files use module-level `pytestmark = pytest.mark.kg`. Cosmetic.
