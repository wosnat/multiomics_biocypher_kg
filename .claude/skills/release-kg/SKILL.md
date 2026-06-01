---
name: release-kg
description: Cut, build, verify, and publish a versioned knowledge-graph release. Use when the user says "release the KG", "cut an alpha", "tag X.Y.Z", "/release-kg <version>", or wants to produce a tagged GitHub Release with a verified staged build. Runs preflight → CHANGELOG cut (pauses for polish) → commit/tag/push → clean clone of the tag → docker build into a staging stack → verify `Schema_info.version` matches the tag → deploy (`staging` default; `local`/`aura` stubbed until the hosting decision lands) → publish GitHub Release with metadata.json. Idempotent re-runs. `--dry-run` exercises every phase without mutating anything.
argument-hint: <version> [--target staging|local|aura] [--allow-dirty] [--draft] [--dry-run] [--resume]
user-invocable: true
allowed-tools: Read, Edit, Write, Bash(uv run *), Bash(python *), Bash(python3 *), Bash(git *), Bash(docker *), Bash(docker compose *), Bash(gh *), Bash(cypher-shell *), Bash(rm -rf /tmp/kg-release-*)
---

# Release KG Skill

Agnostic-core release orchestrator: same workflow whether the eventual host is the local box (Track A) or Aura (Track B). The deploy step is a pluggable seam — only `staging` is wired today; `local` and `aura` raise a clear `NotImplementedError` pointing at the plan section where their design lives.

See `plans/alpha_release.md` §2.3 for the full design and `references/PHASES.md` here for phase-by-phase mechanics.

## When to use

Trigger phrasings: "release the KG", "cut an alpha", "tag a release", "/release-kg X.Y.Z", "produce a release build", "publish kg-X.Y.Z".

Do **not** use for routine `docker compose up -d` rebuilds — those stamp dev `Schema_info` with `0.0.0-dev` automatically via post-import Group 4. This skill is specifically for producing a tagged, GitHub-published release.

## Args

| Arg | Default | Meaning |
|---|---|---|
| `<version>` | — | `X.Y.Z[-(alpha\|beta\|rc).N]`. Becomes git tag `kg-<version>`. |
| `--target` | `staging` | Where to deploy after a verified build. `staging` = leave the staging stack up on `:27687`, no host touched. `local` / `aura` = stubbed. |
| `--mcp-min` | `0.1.0` | Stamped onto `Schema_info.mcp_min_version`; declares min compatible explorer MCP. |
| `--allow-dirty` | off | Skip the working-tree-clean and behind-origin checks. |
| `--draft` | off | Publish GitHub Release as draft. |
| `--dry-run` | off | Every phase logs `[dry-run] would <action>`; mutates nothing (no commits, tags, push, docker, gh). Use to exercise the pipeline. |
| `--resume` | off | Skip the post-CHANGELOG-cut pause (use on the second invocation after polishing). |

## Flow

The script runs six phases in order; each phase is idempotent. The default invocation **pauses once** — after Phase 2 cuts the CHANGELOG — so the operator can polish prose. Re-run with `--resume` to continue from commit onward.

1. **Preflight** — version regex, tooling (`docker`/`gh`/`git`/`cypher-shell`), `gh auth`, repo `.env` present, git on a branch, working tree clean (unless `--allow-dirty`), in sync with origin. Captures SHA / short / branch / dirty.
2. **CHANGELOG cut** — rename `## [Unreleased]` → `## [<version>] - YYYY-MM-DD`, open fresh empty `## [Unreleased]` above. Idempotent: if `## [<version>]` already exists, no-op.
3. **Commit + tag + push** — `chore(release): kg-<version>` commit, annotated tag `kg-<version>`, `git push --follow-tags`. Each step idempotent.
4. **Clean clone of the tag** — `git clone --branch kg-<version> --depth 1 origin /tmp/kg-release-<version>`, copy `.env` in, set release-stamp env vars (`KG_RELEASE_VERSION`, `KG_GIT_*`, `KG_MCP_MIN_VERSION`).
5. **Build into staging stack** — `docker compose -p kg-release-staging up --build -d build import post-process deploy` on **temp ports `:27474` / `:27687`** (does **not** touch the dev `:7474` / `:7687`). Post-import Group 4 stamps `Schema_info` from the env vars set above.
6. **Verify** — query staging `Schema_info`; assert `version == tag`; capture counts. (Explorer smoke test deferred until MCP compatibility contract lands.)
7. **Deploy** — `--target staging`: leave staging up, print the Bolt URI. `local`/`aura`: raise `NotImplementedError` with a pointer to the plan.
8. **Publish** — compose `metadata.json` (version + sha + counts + timestamps), extract the `[<version>]` CHANGELOG section to a notes fragment, `gh release create kg-<version> --notes-file <fragment> --prerelease`, `gh release upload` the manifest.

## Examples

```bash
# Dry-run end-to-end — exercises every phase, mutates nothing
uv run python .claude/skills/release-kg/release_kg.py 0.1.0-alpha.1 --dry-run

# Real cut, target staging — pauses after CHANGELOG cut
uv run python .claude/skills/release-kg/release_kg.py 0.1.0-alpha.1

# (operator polishes CHANGELOG.md in $EDITOR)

# Resume after polish
uv run python .claude/skills/release-kg/release_kg.py 0.1.0-alpha.1 --resume
```

## What this skill does NOT do (yet)

- **Track A deploy (`--target local`)** — blue/green volume flip, firewall allowlist, shared `explorer` login. Gated on hosting decision; see `plans/alpha_release.md` §2.2–2.6.
- **Track B deploy (`--target aura`)** — `neo4j-admin database upload`, real `reader` role, per-user accounts. Gated on hosting decision; see §7.3.
- **Explorer smoke test in Phase 6** — **out of scope** for this skill (the MCP compatibility contract is explorer-repo work, handled separately). The KG side — `Schema_info.mcp_min_version` — is already stamped by post-import Group 4. Operators wanting a smoke test today should run it manually against the staging Bolt URI.
- **`--skip-rebuild`** — Schema_info-only path against an already-live alpha stack. Will land alongside `--target local`.
- **Auto-fill changelog stats** — the operator writes the prose; `metadata.json` carries the authoritative numbers.

## Gotchas

- **`--dry-run` only needs `git`** — preflight skips the `docker` / `gh` / `cypher-shell` availability checks and the `gh auth status` check when `--dry-run` is set, since those tools are never invoked in a dry run. Real runs require all four.
- **First real `--target staging` run** needs `docker-compose.yml` to parameterize the deploy port bindings (currently hardcoded to 7474/7687). That edit is small (`${KG_DEPLOY_HTTP_BIND:-127.0.0.1:7474}:7474` etc.) but lands as a separate follow-up so this skill commit stays self-contained. `--dry-run` runs fine without it.
- **`--resume` is your friend after the CHANGELOG pause** — re-running without it just re-prints "polished? rerun with --resume" because the script detects the cut is already done.
- **Idempotency cuts both ways** — re-running with the same version when the release already exists on GitHub will skip the commit/tag steps but `gh release create` will fail with "release already exists." That's intentional; bump the version or delete the release first.
- **Branch behind origin** triggers a fatal preflight error. Pull (merge or rebase, your call) before retrying.
