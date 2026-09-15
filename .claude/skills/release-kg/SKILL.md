---
name: release-kg
description: Cut, build, verify, and publish a versioned knowledge-graph release. Use when the user says "release the KG", "cut an alpha", "tag X.Y.Z", "/release-kg <version>", or wants to produce a tagged GitHub Release with a verified staged build. Runs preflight → CHANGELOG cut (pauses for polish) → commit/tag/push → clean clone of the tag → docker build into a staging stack → verify `Schema_info.version` matches the tag → deploy (`staging` default; `local` deploys to the Track A lab box via a blue/green volume flip + shared `explorer` login) → publish GitHub Release with metadata.json. Idempotent re-runs. `--dry-run` exercises every phase without mutating anything.
argument-hint: <version> [--target staging|local] [--allow-dirty] [--draft] [--dry-run] [--resume] [--skip-kg-tests]
user-invocable: true
allowed-tools: Read, Edit, Write, Bash(uv run *), Bash(python *), Bash(python3 *), Bash(git *), Bash(docker *), Bash(docker compose *), Bash(gh *), Bash(cypher-shell *), Bash(rm -rf /tmp/kg-release-*)
---

# Release KG Skill

Release orchestrator for the Track A alpha (lab box at `132.75.249.47`). The build/verify/publish core is host-agnostic by design; the deploy step is the pluggable seam. `--target staging` (default) and `--target local` (the lab-box blue/green deploy) are both fully wired (`local` landed in kg-0.1.0-alpha.5). Aura (was Track B) was archived 2026-06-06; see `plans/alpha_release.md` §7.3.

See `plans/alpha_release.md` §2.3 for the full design and `references/PHASES.md` here for phase-by-phase mechanics.

## When to use

Trigger phrasings: "release the KG", "cut an alpha", "tag a release", "/release-kg X.Y.Z", "produce a release build", "publish kg-X.Y.Z".

Do **not** use for routine `docker compose up -d` rebuilds — those stamp dev `Schema_info` with `0.0.0-dev` automatically via post-import Group 4. This skill is specifically for producing a tagged, GitHub-published release.

## Args

| Arg | Default | Meaning |
|---|---|---|
| `<version>` | — | `X.Y.Z[-(alpha\|beta\|rc).N]`. Becomes git tag `kg-<version>`. |
| `--target` | `staging` | Where to deploy after a verified build. `staging` = leave the staging stack up on `:27687`, no host touched. `local` = Track A lab-box deploy on `132.75.249.47:17474`/`:17687` via a blue/green volume flip, shared `explorer` login provisioning, and a `DOCKER-USER` firewall-presence check (warn-only). |
| `--mcp-min` | `[tool.release-kg].mcp_min_version` in `pyproject.toml` (hard fallback `0.1.0` if missing) | Stamped onto `Schema_info.mcp_min_version`; declares min compatible explorer MCP. Per-release override; for repo-wide bumps edit pyproject. |
| `--allow-dirty` | off | Skip the working-tree-clean and behind-origin checks. |
| `--draft` | off | Publish GitHub Release as draft. |
| `--dry-run` | off | Every phase logs `[dry-run] would <action>`; mutates nothing (no commits, tags, push, docker, gh). Use to exercise the pipeline. |
| `--resume` | off | Skip the post-CHANGELOG-cut pause (use on the second invocation after polishing). |
| `--skip-kg-tests` | off | Emergency override: skip the `pytest tests/kg_validity/ --neo4j-url <staging>` step that gates the release in Phase 5. Default is to block the release if any KG-validity test fails. |

## Flow

The script runs six phases in order; each phase is idempotent. The default invocation **pauses once** — after Phase 2 cuts the CHANGELOG — so the operator can polish prose. Re-run with `--resume` to continue from commit onward.

1. **Preflight** — version regex, tooling (`docker`/`git`/`gh`), `gh auth`, repo `.env` present, git on a branch, working tree clean (unless `--allow-dirty`), in sync with origin. Captures SHA / short / branch / dirty. Also runs a **non-fatal unlogged-data check**: if any `paperconfig.yaml` / `paperconfig_files.txt` / `cyanobacteria_genomes.csv` / `treatment_organisms.csv` changed since the last `kg-*` tag but `[Unreleased]` has no `### Data` entries, it prints a WARNING listing the files and continues — it cannot tell a semantic change from a whitespace fix, so it never blocks. Supplementary CSVs are excluded (they churn on every re-resolve without changing what the KG asserts). A second non-fatal check, **unlogged-vocabulary**, warns when `config/controlled_vocabularies.yaml` changed since the last tag but no `### Added` / `### Changed` / `### Breaking` entry mentions the vocabulary — every such change moves `Schema_info.controlled_vocabularies_hash`, which consumers pin. (Phase 6 uses `docker exec staging-deploy cypher-shell …` — the image ships it — so host `cypher-shell` is **not** required.)
2. **CHANGELOG cut** — rename `## [Unreleased]` → `## [<version>] - YYYY-MM-DD`, open fresh empty `## [Unreleased]` above. Idempotent: if `## [<version>]` already exists, no-op.
3. **Commit + tag + push** — `chore(release): kg-<version>` commit, annotated tag `kg-<version>`, `git push --follow-tags`. Each step idempotent.
4. **Clean clone of the tag** — `git clone --branch kg-<version> --depth 1 origin /tmp/kg-release-<version>`, copy `.env` in, set release-stamp env vars (`KG_RELEASE_VERSION`, `KG_GIT_*`, `KG_MCP_MIN_VERSION`).
5. **Build into staging stack** — `docker compose -p kg-release-staging up --build -d build import post-process deploy` on **temp ports `:27474` / `:27687`** (does **not** touch the dev `:7474` / `:7687`). Post-import Group 4 stamps `Schema_info` from the env vars set above. **Then verifies in two layers:** (L0) query staging `Schema_info`, assert `version == tag`, capture counts, and **enforce the vocabulary contract**: recompute `vocabularies_hash()` from the tag checkout's `config/controlled_vocabularies.yaml` and die if `Schema_info.controlled_vocabularies_hash` is null (the build's `controlled_vocabularies.sha256` never reached post-import — post-import itself only warns) or differs (staged graph not built from this checkout); the `--target local` alpha-color build is checked the same way; (L2) run `uv run pytest tests/kg_validity/ --neo4j-url <staging Bolt> -q` to gate the release on structural and semantic validity (~73 s, 1012 assertions). Use `--skip-kg-tests` to bypass in emergencies. On any failure, the staging stack is left running for inspection. (Explorer smoke test = L4, out of scope here.)
7. **Deploy** — `--target staging`: leave staging up, print the Bolt URI. `--target local`: build+verify the alpha-inactive color via a transient `kg-alpha-build` project, tear it down (volume kept), then flip the live `alpha-deploy` onto the new color bound to `ALPHA_BIND_IP` with `up -d --no-deps deploy` (the `--no-deps` is load-bearing — it mounts the already-built/stamped color volume instead of letting `deploy`'s `depends_on` re-run build→import→post-process unstamped), provision the shared `explorer` login (`CREATE OR REPLACE USER`), check the `DOCKER-USER` firewall allowlist (warn-only — needs root), and write `.alpha_active_color`.
8. **Publish** — compose `metadata.json` (version + sha + counts + timestamps + `per_publication_edges` map + `controlled_vocabularies` manifest `{hash, entry_count, entry_ids}` — the hash consumers pin), extract the `[<version>]` CHANGELOG section to a notes fragment, **find the most recent published `kg-*` release on origin and download its `metadata.json`**, render a "What changed since kg-X.Y.Z" diff block (headline-count deltas from Schema_info + per-publication added/changed/removed + a "Controlled vocabularies" section whenever the hash moved: old → new hash, added/removed entry ids, or "same entry set — value lists changed"), prepend the diff block to the notes fragment if non-empty, `gh release create kg-<version> --notes-file <fragment> --prerelease`, `gh release upload` the manifest. First-ever releases or older prior releases without `metadata.json` skip the diff (non-fatal).

## Bringing up an existing release on a new machine (`--bringup`)

The "just run it" path: no cut, tag, push or publish — stand up an **already-published** `kg-<version>` on this box on the release port with the shared `explorer` login, exactly as `--target local` would have left it.

```bash
git clone git@github.com:wosnat/multiomics_biocypher_kg.git && cd multiomics_biocypher_kg
cp .env.alpha.example .env.alpha   # fill ALPHA_BIND_IP / NEO4J_AUTH / ALPHA_EXPLORER_PASSWORD
uv run python .claude/skills/release-kg/release_kg.py 0.1.0-alpha.6 --bringup --target local
```

What it does: verifies the tag exists on origin → creates an empty `.env` if missing (build.sh copies it; nothing at build time needs a value — every build input is committed at the tag, so **no MNX, eggNOG DB or API keys**) → `gh release download <tag> metadata.json` (the release's published stamp; a missing manifest only warns) → clean clone of the tag → derives the stamp (SHA from the tag, `mcp_min` from the manifest else the clone's `pyproject.toml`, vocabulary manifest from the clone) → the same blue/green alpha deploy as Phase 6 `--target local` (build into the inactive color on a temp port, assert `Schema_info.version` + `controlled_vocabularies_hash`, flip `alpha-deploy` onto `${ALPHA_BIND_IP}:17687`, `CREATE OR REPLACE USER explorer`, firewall check) → **compares the rebuilt `Schema_info` against `metadata.json`** (papers / experiments / genes / organisms / expression edges, vocabulary hash, git SHA) and dies — stack left up for inspection — if the rebuild does not reproduce the release. Needs `docker`, `git`, `gh` (auth'd) and ~1 h for the build. `--dry-run` walks every step. Staging is not supported here (it is a verification stack, not a deployment). A prebuilt Neo4j dump asset would make this minutes instead of an hour — backlog.

## Examples

```bash
# Dry-run end-to-end — exercises every phase, mutates nothing
uv run python .claude/skills/release-kg/release_kg.py 0.1.0-alpha.1 --dry-run

# Real cut, target staging — pauses after CHANGELOG cut
uv run python .claude/skills/release-kg/release_kg.py 0.1.0-alpha.1

# (operator polishes CHANGELOG.md in $EDITOR)

# Resume after polish
uv run python .claude/skills/release-kg/release_kg.py 0.1.0-alpha.1 --resume

# New machine: bring up an existing release on the release port + explorer login
uv run python .claude/skills/release-kg/release_kg.py 0.1.0-alpha.6 --bringup --target local
```

## What this skill does NOT do (yet)

- **Aura (`--target aura`)** — archived 2026-06-06 (leadership chose Track A). See `plans/alpha_release.md` §7.3 for the analysis; not callable from this skill.
- **Explorer smoke test in Phase 6** — **out of scope** for this skill (the MCP compatibility contract is explorer-repo work, handled separately). The KG side — `Schema_info.mcp_min_version` — is already stamped by post-import Group 4. Operators wanting a smoke test today should run it manually against the staging Bolt URI.
- **`--skip-rebuild`** — Schema_info-only path against an already-live alpha stack. Still pending; `--target local` has landed without it, so each release rebuilds the alpha color from scratch (two full builds: staging verify + alpha).
- **Auto-fill changelog stats** — the operator writes the prose; `metadata.json` carries the authoritative numbers.

## Gotchas

- **`--dry-run` only needs `git`** — preflight skips the `docker` / `gh` availability checks and the `gh auth status` check when `--dry-run` is set, since those tools are never invoked in a dry run. Real runs require all three (`docker`, `git`, `gh`); `cypher-shell` comes from the staging container via `docker exec`.
- **Staging container names** — the override `docker-compose.staging.yml` renames `build` → `staging-build`, … `deploy` → `staging-deploy`. Without it, the literal `container_name:` directives in the base file would collide with a running dev stack.
- **Stale `staging-deploy` locks the volume on re-runs (auto-handled).** The staging `deploy` container is a long-lived service, so a prior release leaves it `Up`. On the next release its lock on the Neo4j volume makes Phase 5's `neo4j-admin import` fail with *"The database is in use."* (`staging-import` exits 1). Phase 5 now runs `docker compose -p kg-release-staging … down --remove-orphans` (keeps the named volume) **before** the build to release the lock — no-op on a first run. If you ever hit this on an older release, tear staging down manually then re-`--resume`. (Caused the alpha.6 first-cut Phase 5 failure; fixed.)
- **Run from the repo root, not the `/tmp/kg-release-<version>` clone.** The script derives the git branch from its cwd. If you `cd` into the clone (e.g. to inspect staging logs) and then launch `--resume` from there, preflight reports `branch: HEAD` (detached) and Phase 3's `git push origin HEAD` fails with *"destination … is not a full refname."* Always invoke from `/home/osnat/github/multiomics_biocypher_kg`. (Caused the alpha.6 second-cut Phase 3 failure; operator error, not a code bug.)
- **The stamp SHA is re-read from the tag after Phase 3** (fixed after alpha.7). Preflight's HEAD is the pre-cut commit; the `chore(release)` commit moves HEAD, and `--bringup` compares the tag clone's SHA to `metadata.json.git_sha`. alpha.7's graph carries the pre-cut `6bf6e513`; its `metadata.json` was patched to the tag's `9a2c9610` so bringup works.
- **`--resume` is your friend after the CHANGELOG pause** — re-running without it just re-prints "polished? rerun with --resume" because the script detects the cut is already done.
- **Idempotency cuts both ways** — re-running with the same version when the release already exists on GitHub will skip the commit/tag steps but `gh release create` will fail with "release already exists." That's intentional; bump the version or delete the release first.
- **Branch behind origin** triggers a fatal preflight error. Pull (merge or rebase, your call) before retrying.
- **The live `--target local` flip uses `up -d --no-deps deploy`.** Dropping `--no-deps` lets the `deploy` service's `depends_on` re-run build→import→post-process in the live `kg-alpha` project — and because that up intentionally omits `KG_RELEASE_VERSION`, the re-import re-stamps `Schema_info` to the `0.0.0-dev` fallback. Symptom: lab graph serves correct *data* but `version=0.0.0-dev`. (Caused the alpha.5 first-cut failure; fixed.)
- **Explorer provisioning uses `CREATE OR REPLACE USER`,** not `CREATE IF NOT EXISTS` + `ALTER`. On a fresh blue/green volume the latter's `ALTER` to the just-created (identical) password raises "Old and new password cannot be the same" and aborts the deploy after the graph is already live. (Also caused the alpha.5 first-cut failure; fixed.)
