---
name: refresh-mnx
description: Build or refresh the MetaNetX (MNX) metabolite/reaction cross-reference database — downloads the 4 MNX TSVs (~1.5 GB) and builds the SQLite resolver (~2.5 GB, ~30 min). Use whenever the user wants to create, rebuild, refresh, or troubleshoot the MNX database / resolver / cross-reference cache, or mentions a new MNX release, a stale resolver, or a missing `metabolite_resolver.db`. The KG metabolism pipeline (prepare_data step 6) depends on this — if step 6 errors complaining about a missing resolver or `MNX_DATA_DIR`, this skill is the fix.
argument-hint: [--force]
user-invocable: true
allowed-tools: Read, Bash(bash *), Bash(uv *), Bash(df *), Bash(du *), Bash(ls *), Bash(echo *), Bash(grep *), Bash(printenv *)
---

# Refresh MNX Database Skill

Wrapper + troubleshooting checklist around `scripts/refresh_mnx.sh`. The MNX database is a global cross-reference resolver (compounds + reactions across KEGG, ChEBI, HMDB, BiGG, …) used by `prepare_data.sh --steps 6` to enrich KEGG with cross-refs. **Run this only when MNX releases a new version or the resolver is missing/corrupt** — it's a ~30-minute one-shot, not a routine step.

## Quick Start

```bash
# Download + build (skips if resolver already present)
bash scripts/refresh_mnx.sh

# Force re-download + rebuild (use after a new MNX release)
bash scripts/refresh_mnx.sh --force
```

After completion, **remind the user** (do not run automatically):

```bash
bash scripts/prepare_data.sh --steps 6 --force
```

This regenerates `cache/data/kegg/kegg_data.json` with the refreshed MNX cross-refs. Without it the resolver build is wasted — downstream KG builds keep using the stale enrichment.

## Workflow When Invoked

1. **Pre-flight checks** (catch problems before the 30-minute build):
   - Resolve target dir: `printenv MNX_DATA_DIR` — if unset, default is `cache/data/mnx/` under the project root. Tell the user where files will land.
   - Check free disk space at the target dir: `df -h "$MNX_DATA_DIR"` (or the cache dir). Need ≥ **5 GB** free (1.5 GB raw TSVs + ~2.5 GB SQLite + scratch). Bail with a clear message if not.
   - Show what's already there: `ls -lh "$MNX_DATA_DIR" 2>/dev/null || echo "(empty / does not exist)"`. The 4 expected TSVs are `chem_prop.tsv`, `chem_xref.tsv`, `reac_prop.tsv`, `reac_xref.tsv`; the resolver is `metabolite_resolver.db`.
   - If `metabolite_resolver.db` already exists and the user did NOT pass `--force`, confirm they want to re-run — without `--force` the rebuild is a no-op and they may have meant to pass it.

2. **Run** `bash scripts/refresh_mnx.sh [--force]`. Tail `logs/refresh_mnx.log` if it appears stuck (the SQLite build is silent for stretches; that's normal).

3. **Verify outputs** after completion:
   - All 4 TSVs present and non-empty in `$MNX_DATA_DIR`.
   - `metabolite_resolver.db` present, ~2.5 GB.
   - No traceback at the tail of `logs/refresh_mnx.log`.

4. **Tell the user** to run `bash scripts/prepare_data.sh --steps 6 --force` next. **Do not run it for them** — step 6 has its own runtime cost and the user may want to defer or batch it.

## Troubleshooting

### Pre-build: where will files land?

`MNX_DATA_DIR` is resolved by [`multiomics_kg/utils/metabolite_utils.py`](../../../multiomics_kg/utils/metabolite_utils.py) — env var wins, otherwise falls back to `cache/data/mnx/` under the project root. The env var supports `~` expansion. To share one ~4 GB cache across multiple checkouts, set it in `.env` (e.g. `MNX_DATA_DIR=~/tools/mnx`) — pre-existing files in that dir will be reused, no rebuild needed.

### Download fails (network / 404)

The download step uses `multiomics_kg.download.download_metabolism_reference`. MNX URLs can occasionally 404 during a release transition. Symptoms: `requests` `HTTPError` in the log. Fixes:
- Retry once — transient.
- Check whether MNX shipped a new release that broke the URL pattern: visit https://www.metanetx.org/mnxdoc/mnxref.html and confirm the file names still match `chem_prop.tsv`, `chem_xref.tsv`, `reac_prop.tsv`, `reac_xref.tsv`.
- Partial download? The script writes to the final path directly, so a half-finished file looks "present". Pass `--force` to redownload.

### SQLite build crashes mid-way

The build is a single transaction per table (compounds → compound_aliases → compound_names → reactions → reaction_aliases). If it crashes, `metabolite_resolver.db` is left in an inconsistent state. Recovery:
- Delete the DB and rerun: `rm "$MNX_DATA_DIR/metabolite_resolver.db" && bash scripts/refresh_mnx.sh`.
- Don't pass `--force` for this — TSVs are fine, only the DB is bad.
- If the crash is `disk full`, free space and retry. SQLite scratch can momentarily exceed the final 2.5 GB.

### Step 6 complains about missing resolver after a successful rebuild

```
... or set MNX_DATA_DIR to a checkout that already has the resolver.
```

This means step 6's process saw a different `MNX_DATA_DIR` than the rebuild did. Common causes:
- `MNX_DATA_DIR` set in your shell but not in `.env` (or vice versa) — the build read one, step 6 reads the other. Align them.
- Tilde not expanded somewhere. Hardcode an absolute path to test.
- Two checkouts of the project, each with its own default `cache/data/mnx/`. Confirm `printenv MNX_DATA_DIR` prints the same path in both shells, or set it explicitly.

### "Resolver DB present; skip (use --force to rebuild)"

That's the build_mnx_resolver script doing the right thing. If you actually want a rebuild (new MNX release, or you suspect corruption), pass `--force`:

```bash
bash scripts/refresh_mnx.sh --force
```

### How long should this take?

- Downloads: ~5-10 min depending on bandwidth (1.5 GB total).
- SQLite build: ~25-30 min, mostly compound_aliases + reaction_aliases inserts. If it's been over an hour, check `logs/refresh_mnx.log` for a hung HTTP retry or a stuck UNIQUE-constraint loop.

## Key Paths

- Wrapper: [scripts/refresh_mnx.sh](../../../scripts/refresh_mnx.sh)
- Download module: `multiomics_kg.download.download_metabolism_reference` (`--sources mnx`)
- Resolver build module: `multiomics_kg.download.build_mnx_resolver`
- Path resolution: [multiomics_kg/utils/metabolite_utils.py](../../../multiomics_kg/utils/metabolite_utils.py) (`get_mnx_data_dir`)
- Log: `logs/refresh_mnx.log`
- Outputs: `$MNX_DATA_DIR/{chem_prop,chem_xref,reac_prop,reac_xref}.tsv` + `metabolite_resolver.db`

## Why Not Run Step 6 Automatically?

The wrapper script's final message and this skill both nudge toward `prepare_data.sh --steps 6 --force` because that's what consumes the new MNX data. But step 6 is a separate ~minutes-long run that touches every strain's annotations and rewrites `cache/data/kegg/kegg_data.json` — the user may want to defer it, batch it with other prepare_data steps, or coordinate it with a Docker rebuild. Surfacing it as the next step (rather than chaining it) keeps the operator in control.
