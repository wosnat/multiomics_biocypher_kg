# PSORTb subcellular-localization skill — Phase 1 (artifacts)

**Status:** Draft (2026-05-10) · **Owner:** Osnat · **Pairs with:** [tcdb-diamond](2026-05-10-tcdb-diamond-augmentation-design.md)

## 1. Motivation

PSORTb v3.0.3 is the gold-standard predictor of bacterial subcellular localization (precision ≥ 95% on its training set). Localization is orthogonal to the functional annotations we already carry (eggNOG, UniProt, TCDB, Pfam) and is a strong routing signal for chemistry / transport / community-interaction questions:

- Distinguishing inner vs. outer-membrane transporters (currently TCDB-family-inferred only).
- Identifying secreted exoproteome candidates against the existing exoproteomics papers.
- Filtering "membrane noise" out of cytoplasmic regulatory comparisons.

Phase 1 produces inspectable per-strain artifacts. KG integration is **explicitly deferred** to a future Phase 2 spec (cf. tcdb-diamond's same staging).

## 2. Scope

### In scope
- New skill `/psortb-run` (parallels `/eggnog-run`, `/tcdb-diamond`).
- Runs PSORTb v3.0.3 via the official Docker image `brinkmanlab/psortb_commandline:1.0.2` against `protein.faa` for every row in `cyanobacteria_genomes.csv` that has a `protein.faa` (currently 32 rows — all 30 `genome_strain` rows plus the 2 `reference_proteome_match` rows, both of which carry NCBI RefSeq FASTAs). Rows are loaded via the shared `multiomics_kg.download.utils.cli.load_genome_rows` helper; rows without `protein.faa` are skipped with a `MISSING_INPUT` row in the status table.
- Emits per-strain `<strain>.psortb.calls.json` + `<strain>.psortb.skill_summary.json` under `cache/data/<organism>/genomes/<strain>/psortb/`.
- All organisms hardcoded as Gram-negative (`--negative`). The KG has zero Gram-positive strains.

### Out of scope (deferred)
- Schema changes, post-import Cypher, KG validity tests.
- Adding a `subcellular_localization` field to Gene nodes.
- Gram-positive / Archaea support.
- Re-running on protein-FASTA changes (caller invokes `--force`).

## 3. Inputs / outputs

### Input per strain
- `<data_dir>/protein.faa` (NCBI RefSeq protein FASTA, written by step 0). Identifiers are WP_ accessions matching the `protein_id` column of `gene_mapping.csv`.

### Output per strain (`cache/data/<organism>/genomes/<strain>/psortb/`)
- `<strain>.psortb.terse.tsv` — raw PSORTb terse output (kept for traceability, gitignored via `cache/`).
- `<strain>.psortb.calls.json` — per-protein records keyed by WP_ accession:
  ```json
  {
    "WP_011131900.1": {
      "localization": "OuterMembrane",
      "score": 9.97,
      "secondary_localization": null,
      "secondary_score": null,
      "is_multi_localized": false,
      "is_unknown": false
    }
  }
  ```
  Field semantics:
  - `localization` ∈ {`Cytoplasmic`, `CytoplasmicMembrane`, `Periplasmic`, `OuterMembrane`, `Extracellular`, `Unknown`}.
  - `score` ∈ [0, 10]; PSORTb assigns `Unknown` when no class scores ≥ 7.5 (`is_unknown=true`).
  - Multi-localization: the parser handles slash-separated values (e.g. `Cytoplasmic/CytoplasmicMembrane` with scores `9.83/9.71`) by ordering primary/secondary by descending score and setting `is_multi_localized=true`. **However, `--output terse` empirically emits a single `Final_Localization` per protein even when two classes both score ≥ 7.5**, so `is_multi_localized` is effectively always `false` in Phase 1. To expose true multi-loc, switch to `--output long` and parse the per-class score block (Phase 2 if needed).
- `<strain>.psortb.skill_summary.json` — counts per localization class, Unknown rate, multi-localization rate, runtime, image digest.
- Diamond-style status table to stdout; raw container stdout/stderr captured to `logs/psortb_<strain>.log` (already covered by `*.log` gitignore).

## 4. Runner behaviour

```bash
uv run python .claude/skills/psortb-run/run_psortb.py --prepare-image  # one-time: docker pull (~2 GB)
uv run python .claude/skills/psortb-run/run_psortb.py            # all strains, skip already-done
uv run python .claude/skills/psortb-run/run_psortb.py --strain MED4
uv run python .claude/skills/psortb-run/run_psortb.py --force    # re-run even if calls.json exists
uv run python .claude/skills/psortb-run/run_psortb.py --strain MIT1002 --limit 100   # smoke test
```

### Algorithm per strain
1. Resolve `data_dir` from `cyanobacteria_genomes.csv` via the shared
   `multiomics_kg.download.utils.cli.load_genome_rows` helper.
2. Skip if `<data_dir>/psortb/<strain>.psortb.calls.json` exists and `--force` not given.
3. If `<data_dir>/protein.faa` is missing, mark `MISSING_INPUT` and continue.
4. (Optional `--limit N`) — write a temp truncated FASTA of the first N records alongside the original.
5. Run the container:
   ```
   docker run --rm --user $(id -u):$(id -g) \
     -v <data_dir>:/tmp/results \
     brinkmanlab/psortb_commandline:1.0.2 \
     /usr/local/psortb/bin/psort --negative --output terse -i /tmp/results/protein.faa
   ```
   Mount/invocation rationale (discovered during smoke-test):
   - The brinkmanlab `psort` Perl wrapper hardcodes its output directory to
     `/tmp/results/` inside the container and writes
     `<YYYYMMDDHHMMSS>_psortb_gramneg.txt` there — not to stdout. Mounting the
     strain's `data_dir` at `/tmp/results` lands the file in our host directory.
   - `-i` is mandatory; the wrapper does not accept a positional file argument.
   - `--user` makes produced files host-owned (no sudo needed to clean up).
   - Container stdout/stderr captured to `logs/psortb_<strain>.log`.
6. Rename the timestamped output file to `<data_dir>/psortb/<strain>.psortb.terse.tsv`.
7. Parse the 3-column terse TSV (`SeqID`, `Localization`, `Score`); the parser
   splits slash-separated multi-localization values, though `--output terse`
   in practice emits a single `Final_Localization` per protein (see §3 multi-loc note).
8. Write `calls.json` + `skill_summary.json`.

### Image management
- Docker stores the image in its daemon storage (`/var/lib/docker/overlay2/...`), shared across all projects on the host. No per-skill cache directory.
- `--prepare-image` flag pulls `brinkmanlab/psortb_commandline:1.0.2` (~2 GB) and prints its digest, then exits. Run this once before the first batch.
- A normal run (`--strain` or no args) refuses to start unless the image is already present locally — pointing the user at `--prepare-image`. This avoids accidental 2 GB downloads inside a long batch.
- `--refresh-image` re-pulls (e.g. when Brinkman lab releases a new patch). Mutually exclusive with normal-run flags.
- The skill records `image_digest` (from `docker inspect`) in each `skill_summary.json` so artifacts are reproducibility-traceable.

### Concurrency
- PSORTb is single-threaded per protein and the Docker image bundles every classifier; no `--threads` flag passed. Strains run sequentially; users wanting parallelism can spawn multiple invocations with `--strain`.

## 5. Failure modes

| Mode | Behaviour |
|---|---|
| Docker not on PATH | Print install hint, exit 2. |
| Image not present and `--prepare-image` not yet run | Print "run `--prepare-image` first", exit 2. Never silently pulls during a batch. |
| `--prepare-image` pull fails (offline / Docker Hub down) | Print error, exit 1. |
| `protein.faa` missing | Skip strain, report `MISSING_INPUT` in status table. |
| Container non-zero exit | Capture stderr to log, mark strain FAILED, leave any partial files alone. |
| Truncated protein < 20 aa | PSORTb emits `Unknown` — kept as-is (counts toward Unknown rate). |

## 6. Phase 2 (deferred — not part of this change)

- Add `subcellular_localization`, `subcellular_localization_score`, `subcellular_localization_is_multi` properties to Gene nodes via the `functional_annotation_adapter` (or a dedicated psortb_adapter).
- Decide on confidence-floor behaviour for `Unknown` calls (drop vs. keep with `Unknown` literal).
- Post-import: organism-level rollup (`outer_membrane_gene_count`, `extracellular_gene_count`).
- KG validity tests: distribution sanity (e.g. ≥ 0.5% genes Extracellular per strain — the marine surface lifestyle floor).

Phase 2 will get its own design doc once Phase 1 artifacts have been spot-checked against published localization predictions for MED4 (Thompson 2011) and MIT9301 (Bertilsson 2003).

## 7. Acceptance criteria for Phase 1

- `uv run python .claude/skills/psortb-run/run_psortb.py --strain MED4 --limit 100` produces a syntactically valid `MED4.psortb.calls.json` covering 100 proteins in ≤ 5 minutes.
- Full run on all 32 rows completes in a single invocation; each strain ends with a row in the status table; failures do not abort the batch.
- No changes to `schema_config.yaml`, adapters, or post-import scripts.
- `cache/data/**/psortb/` is already covered by the existing `cache/` gitignore — no new ignore rules needed.

## 8. As-built notes (Phase 1 ship, 2026-05-11)

Three discoveries during the first end-to-end run that aren't obvious from PSORTb docs and shaped the final §4 algorithm:

1. **Wrapper hardcodes `/tmp/results/`** — the brinkmanlab `psort` Perl wrapper writes its output to `/tmp/results/<timestamp>_psortb_gramneg.txt` regardless of stdout redirection, container CWD, or `-w`. The original spec said "mount at `/tmp/psortb`, capture stdout" — both turned out to be no-ops. Fix: mount the strain's `data_dir` at `/tmp/results` and rename the timestamped file after the container exits.
2. **`-i` flag is mandatory** — the wrapper does not consume a positional FASTA argument; it requires `-i <path>` or `--seq <path>`. Without it, the wrapper prints "No such file" and the usage banner.
3. **`--output terse` collapses multi-localizations** — the terse format always emits one `Final_Localization` per protein, even when two classes both score ≥ 7.5. The slash-handling parser is still in place (so it'll correctly decompose if we ever switch to `--output long`), but `is_multi_localized` was 0 across all 32 strains in the Phase 1 batch.

Phase 1 batch (32 strains, 96,418 proteins, 8 h sequential): 41.9% Cytoplasmic, 22.0% CytoplasmicMembrane, 1.7% Periplasmic, 1.8% OuterMembrane, 1.0% Extracellular, 31.6% Unknown. OuterMembrane density tracks lifestyle: Prochlorococcus ~0.7% vs heterotrophs (Pseudomonas, Alteromonas) ~3.3%. Spot-check against canonical KT2440 proteins (OprD-family porins, TolC, GroEL, RNA polymerase subunits) all in expected localizations with scores ≥ 9.27.
