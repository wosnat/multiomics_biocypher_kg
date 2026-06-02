# `release-kg` — Phase-by-phase reference

Detailed mechanics. See the parent `SKILL.md` for arg/usage summary and `plans/alpha_release.md` §2.3 for the design rationale.

## Idempotency contract

Every phase is re-runnable. The script captures no state between invocations — the state lives in the repo (CHANGELOG, git log, tags) and on disk (`/tmp/kg-release-<version>`, staging docker volume). Re-running with the same `<version>` either no-ops phases that are already done or reproduces them.

| Phase | Idempotent? | How |
|---|---|---|
| 1 Preflight | yes | Read-only |
| 2 CHANGELOG cut | yes | Skips if `## [<version>]` section already exists |
| 3 Commit/tag/push | yes | Skips commit if HEAD subject matches; skips tag if exists; `git push` is naturally idempotent |
| 4 Clean clone | yes | Removes existing `/tmp/kg-release-<version>` before re-cloning |
| 5 Build + verify | rebuild-heavy | `docker compose up` re-attaches if the stack is already up; the build layer cache makes re-builds fast. Verify is read-only. |
| 6 Deploy | depends on backend | `staging`: no-op deploy. `local`/`aura`: each will document its own idempotency when implemented. |
| 7 Publish | yes | `gh release create` fails fast if release exists (intentional — bump version or delete). |

## Phase 1: Preflight

Fail-fast read-only checks, no mutation.

| Check | Failure mode | Override |
|---|---|---|
| Version matches `^\d+\.\d+\.\d+(-(alpha\|beta\|rc)\.\d+)?$` | exit 1 | — |
| Tools on PATH: `docker`, `git`, `gh` (cypher-shell is invoked via `docker exec staging-deploy …`, no host install needed) | exit 1 | — |
| `gh auth status` succeeds | exit 1 | run `gh auth login` |
| Repo `.env` exists | exit 1 | create the file (build.sh hardcodes `cp /src/.env .`) |
| `git status --porcelain` is empty | exit 1 | `--allow-dirty` |
| `HEAD..origin/<branch>` count is 0 | exit 1 | `--allow-dirty` |

On success, captures `git_sha`, `git_sha_short`, `git_branch`, `git_dirty` into the run context — these flow into the `Schema_info` stamp via env vars in Phase 5.

## Phase 2: CHANGELOG cut

Uses the Keep-a-Changelog "accumulate-then-cut" model.

- Reads `CHANGELOG.md` from repo root.
- Locates `## [Unreleased]` (fatal if missing).
- Captures the section body (everything until the next `## [` or EOF).
- Rewrites the file as:

```markdown
## [Unreleased]

### Added

### Changed

### Fixed

## [<version>] - <YYYY-MM-DD>

<the old Unreleased body, verbatim>
```

- Idempotent: if `## [<version>]` already exists anywhere in the file, no-op.
- If `[Unreleased]` had no `-` bullets, the new version section gets `_No entries — placeholder._` so the section is still well-formed.

After this phase, the default invocation **pauses** (return 0) with the message:

```
=== PAUSE: CHANGELOG cut done ===
Review and polish CHANGELOG.md as needed, then re-run with --resume to continue.
```

`--resume` skips the pause. `--dry-run` skips the pause (and doesn't actually write).

## Phase 3: Commit + tag + push

- Subject: `chore(release): kg-<version>`. Skipped if HEAD already has it.
- Annotated tag: `kg-<version>` with message `KG release <version>`. Skipped if tag exists locally.
- Push: `git push origin <branch> --follow-tags`. Naturally idempotent.

The tag is pushed *now* (not after the build) so Phase 4 can clone it.

## Phase 4: Clean clone of the tag

- `git clone --branch kg-<version> --depth 1 <origin-url> /tmp/kg-release-<version>`.
- Removes any pre-existing `/tmp/kg-release-<version>` first.
- Copies the repo's gitignored `.env` into the clone (build.sh inside Docker does `cp /src/.env .` and will abort under `set -e` otherwise).

The build reads the **tag**, never the live working tree — so dev work can continue during a release and the release is reproducible from `git clone --branch kg-<version>`.

**Determinism caveat:** the build still consumes committed `cache/` and the warm `build_cache` volume. Reproducibility is "clean clone of the tag + warm cache," not fully hermetic.

## Phase 5: Build into staging + verify

Brings up a parallel docker stack named `kg-release-staging` on temp ports `127.0.0.1:27474` (HTTP) / `127.0.0.1:27687` (Bolt). The dev stack on `:7474` / `:7687` is untouched.

Env vars set for `post-process` (consumed by post-import.sh Group 4 to stamp `Schema_info`):

```
KG_RELEASE_VERSION=<version>
KG_GIT_SHA=<full sha>
KG_GIT_SHA_SHORT=<short sha>
KG_GIT_BRANCH=<branch>
KG_GIT_DIRTY=true|false
KG_MCP_MIN_VERSION=<--mcp-min, default 0.1.0>
KG_DEPLOY_HTTP_BIND=127.0.0.1:27474
KG_DEPLOY_BOLT_BIND=127.0.0.1:27687
```

**Isolation from the dev stack:** Phase 5 passes `-f docker-compose.yml -f docker-compose.staging.yml`. The override renames the four `container_name`s with a `staging-` prefix (`staging-build`, `staging-import`, `staging-post-process`, `staging-deploy`) so the staging containers don't collide with a running dev stack — which holds the unprefixed literal names from the base file. Compose project namespacing already isolates the named volumes.

After bringing up `deploy`, the script polls `docker exec staging-deploy cypher-shell -a bolt://localhost:7687 RETURN 1;` for up to 120s, then queries `Schema_info`. Talking to Bolt over the container's loopback (port 7687, not the host-published 27687) removes the host `cypher-shell` dependency. Asserts `s.version == <tag>`. Captures `gene_count` / `experiment_count` / `paper_count` / `organism_count` / `expression_edge_count` / `built_at` into context for Phase 7.

**Explorer smoke test:** **out of scope** for this skill (the MCP compatibility contract is explorer-repo work, handled separately — decided 2026-06-01). The KG side of the contract is `Schema_info.mcp_min_version`, already stamped by post-import Group 4. Phase 5 logs an "out-of-scope" note; operators wanting a smoke test today should point the explorer at the staging Bolt URI manually.

## Phase 6: Deploy

The pluggable seam — the entire skill spine drains into one of three backend functions selected by `--target`:

- **`staging`** (implemented): no-op. Leaves the staging stack up on `:27687` so the operator can poke around. Tear down with `docker compose -p kg-release-staging down`.
- **`local`** (stub): raises `NotImplementedError`. Will implement when the hosting decision lands Track A — blue/green volume flip + provision shared `explorer` login + firewall confirmation. See plan §2.2–2.6.
- **`aura`** (stub): raises `NotImplementedError`. Will implement for Track B — `neo4j-admin database dump` the staging volume, `neo4j-admin database upload --to-uri neo4j+s://…`, verify `Schema_info` round-trip, ensure `reader` role + user. See plan §7.3.

## Phase 7: Publish

- **`metadata.json`** — written to the clone dir. Shape:

```json
{
  "version": "0.1.0-alpha.1",
  "tag": "kg-0.1.0-alpha.1",
  "git_sha": "…",
  "git_sha_short": "…",
  "git_branch": "main",
  "git_dirty": false,
  "mcp_min_version": "0.1.0",
  "target": "staging",
  "built_at": "<from Schema_info.built_at>",
  "counts": {
    "papers": …, "experiments": …, "genes": …,
    "organisms": …, "expression_edges": …
  },
  "stamped_at": "<utc iso now>"
}
```

This is the authoritative numbers carrier — the CHANGELOG holds prose, `metadata.json` holds counts.

- **Release-notes fragment** — `.release-notes-<version>.md` extracted from the CHANGELOG's `[<version>]` section.

- **GitHub Release** — `gh release create kg-<version> --notes-file <fragment> --prerelease` (`--draft` adds `--draft`). The `--prerelease` flag is unconditional because any `X.Y.Z-(alpha|beta|rc).N` tag is a pre-release and bare `X.Y.Z` is rare here.

- **Asset upload** — `gh release upload kg-<version> metadata.json`.

- **Operator checklist** — printed at the end: where to find the release, the staging URI, tear-down commands, clone cleanup.

## Tester announcement template

After a successful release:

> **KG release `kg-<version>` is live.**
>
> - Release notes + manifest: `gh release view kg-<version>` (or the GitHub URL).
> - Connect: Bolt URI / `.mcp.json` snippet — see `docs/kg_mcp_guide.md` (and the connection section once the hosting decision lands).
> - `Schema_info.mcp_min_version = <mcp-min>` — make sure your explorer MCP is at least that version.
> - File issues at `<repo>/issues` with the version tag in the title.

## Common operator flows

```bash
# Sanity-check the pipeline against the current commit
uv run python .claude/skills/release-kg/release_kg.py 0.1.0-alpha.1 --dry-run

# Cut, polish, finish — the typical real release
uv run python .claude/skills/release-kg/release_kg.py 0.1.0-alpha.1
# … review CHANGELOG.md, edit if needed …
uv run python .claude/skills/release-kg/release_kg.py 0.1.0-alpha.1 --resume

# Release a draft (don't publish publicly yet)
uv run python .claude/skills/release-kg/release_kg.py 0.1.0-alpha.1 --draft

# Release with a dirty working tree (rare, e.g. emergency)
uv run python .claude/skills/release-kg/release_kg.py 0.1.0-alpha.1 --allow-dirty
```
