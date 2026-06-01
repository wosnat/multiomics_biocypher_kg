# Plan: Alpha release of the KG — agnostic core now; hosting (local box OR Aura) pending budget

**Status:** Draft — hosting decision **reopened 2026-05-25** (local box vs Aura, decided on budget grounds by leadership). Near-term work proceeds on the **deployment-agnostic core** (see "Near-term work" below); none of it depends on the choice. Local-box design retained as **Track A** (§2.2–2.6); Aura as **Track B** (§7.3). Supersedes the lab-only-final framing of 2026-05-24.
**Goal:** Cut tagged, versioned alpha releases of the KG served from a *second* Neo4j stack on the same Linux dev/release box (`132.75.249.47`), so the lab's alpha testers — **physically in the lab, on the lab subnet** — can drive the graph from their own LLM agents via the explorer MCP + research-skills plugin.

## Scope decision log

The access mechanism narrowed across several decisions; recording them so we don't relitigate:

- **Hosting reopened (2026-05-25).** The "Not Aura" call below is **superseded**: leadership will pick **local box vs Aura on budget grounds**. Until that lands, work proceeds on the **deployment-agnostic core** (see "Near-term work" below) — none of it depends on the choice. Both analyses stand: local-box design in §2.2–2.6, Aura in §7.3.
- **Not Aura (2026-05-24).** **[Superseded 2026-05-25 — see entry above.]** Aura Pro is the smallest tier that fits our graph (250K nodes / 2.2M rels > Aura Free's 200K/400K cap) and runs ~$1.5K+/yr. Local box is $0 and iterates in minutes (`docker compose up -d` on a tagged commit) vs. tens of minutes for dump→upload→restore. Deferred to §7.
- **Not Tailscale (2026-05-24).** The original local-box plan reached testers over Tailscale. Dropped because: (a) all current testers are **lab-local**; (b) the off-site case is already covered by the **university VPN**, which puts a home machine back on the campus network — Tailscale would be a redundant second overlay; (c) Tailscale's free Personal plan is **non-commercial-use + user-capped**, and the clean paid tier (Standard, $8/user/mo ≈ $576/yr for 6 users) is real money for redundant reach. Deferred to §7.
- **Lab-only (this draft).** Simplest, $0, smallest exposure. The box has a campus-routable IP; bind the alpha Neo4j to it, turn on Neo4j auth, and **firewall the alpha ports to the lab subnet**. Trivially expandable later (add the VPN-pool CIDR to the allowlist — one rule) without rearchitecting. **[Now the Track-A choice, conditional on picking local.]**

## Near-term work: deployment-agnostic core (do now, before the hosting decision)

The dividing line: **agnostic** = build the graph → stamp its identity → version/release it → tell testers how to drive the MCP. **Forked** = where the DB runs, how testers reach it, and how read-only is enforced. The following is identical under local **and** Aura and proceeds now; the fork is isolated to a single pluggable deploy step.

| # | Item | Status | Section |
|---|---|---|---|
| 1 | `Schema_info` release metadata (version/git_sha/counts, post-import) — **the foundation** | to do | §2.1 |
| 2 | `CHANGELOG.md` + tag scheme `kg-X.Y.Z-alpha.N` + GitHub Release flow | to do | §2.3 |
| 3 | `/release-kg` skill — agnostic spine (preflight, changelog, commit/tag/push, clean-clone-of-tag, publish); the deploy step is the **one seam** that forks per host | to do | §2.3 |
| 4 | KG-validity tests `test_schema_info_release_properties`, `test_schema_info_counts_match` (run vs dev `localhost`) | to do | §3 |
| 5 | ~~MCP compatibility contract (explorer check + version-pin cadence)~~ | **out of scope (2026-06-01)** — explorer-repo work, handled separately. The KG side of the contract (`Schema_info.mcp_min_version`) is in place from item 1. | §2.1, §6.4(5) |
| 6 | `docs/kg_mcp_guide.md` body — refactor so all connection specifics sit in **one** swappable section | ✅ **done** — new §2 "Connection details" carries host-A-vs-B tables (URI, reach-the-host check, credentials, read-only model); rest of the guide is host-agnostic. | §2.7 |

**Deployment-specific — deferred until the hosting decision lands:**
- **Track A (local):** `docker-compose.alpha.yml`, `.env.alpha`, `alpha_up/down.sh`, `DOCKER-USER` firewall allowlist, shared `explorer` login + MCP-enforced read-only, blue/green volume flip, `docs/kg_alpha_it_approval.md` (§2.2–2.6, §6).
- **Track B (Aura):** provision the Aura instance, dump→upload→restore release path, real `reader` role for read-only, `neo4j+s://` URI + per-user accounts (§7.3).
- **Either way:** the connection section of `docs/kg_mcp_guide.md` (item 6 leaves a hole for it).

## Workstreams

| ID | Title | Scope |
|---|---|---|
| **A** | Two-stack local Docker (dev + alpha side-by-side) | **in scope** |
| **B** | Lab-network access (bind to lab IP + Neo4j auth + host-firewall allowlist) | **in scope** |
| **C** | Alpha-user-facing guide (two-repo install: `multiomics_explorer` + `multiomics_research`) | **in scope** |
| **D** | Remote & multi-user access (university-VPN allowlist, Tailscale/Cloudflare, **remote-hosted MCP + per-user auth**, Enterprise/Aura read-only roles) | **deferred** — see §7 |

For the alpha: each lab tester installs the explorer MCP and the research-skills plugin locally and points the MCP at `bolt://132.75.249.47:17687` with a **single shared `explorer` credential** (read access; read-only is enforced by the MCP, not the DB — see §2.5). Per-user identity is low-ROI here and deferred to the remote-MCP-hosting / Enterprise phase (§7).

## 1. Scope

### What is changing

**A. Two-stack Docker layout.**
- Keep `container_name:` on the base file (dev keeps `build`/`import`/`post-process`/`deploy`/`app`, so the pervasive `docker exec deploy …` / `docker logs import` workflow keeps working). The alpha **override sets distinct names** (`alpha-deploy`, etc.) so the two projects never collide. (See review note B2.)
- New `docker-compose.alpha.yml` override: different ports (`17474`/`17687`), separate import-artifact mount (`./output-alpha`), auth on, deploy bound to the lab IP, no `app` (Biochatter).
- New `.env.alpha.example` with `KG_RELEASE_VERSION`, `KG_GIT_SHA`, `KG_GIT_BRANCH`, `ALPHA_BIND_IP`, `NEO4J_AUTH`.
- New `scripts/alpha_up.sh` / `scripts/alpha_down.sh` wrappers around the `-p kg-alpha -f … -f …` invocations.

**B. Lab-network access + auth.**
- Bind the alpha `deploy` ports to the box's **lab IP only** (`${ALPHA_BIND_IP}:17687:7687`), never `0.0.0.0`.
- **Host-firewall allowlist** restricting `17687`/`17474` to the lab subnet — required, because the bind IP is campus-routable (§2.6).
- Neo4j auth ON in the alpha stack (`NEO4J_AUTH=neo4j/<strong-shared-password>`, `NEO4J_dbms_security_auth__enabled: true`).
- **Shared `explorer` login**: one shared read credential (not the `neo4j` admin), distributed out-of-band. Read-only is enforced by the **explorer MCP**, not the DB (Community has no roles → all logins are full-access anyway). Per-user accounts are low-ROI here and deferred (§2.5, §7).

**C. Existing pieces retained.**
- `Schema_info` release-metadata properties (§2.1) — decouples build-time vs. release-time identity.
- `CHANGELOG.md` at repo root (Keep a Changelog format).
- `/release-kg` skill — reshaped around the local two-stack, lab-only model (§2.3).
- Tag scheme `kg-X.Y.Z[-(alpha|beta|rc).N]`; GitHub Release per tag, pre-release for any non-stable suffix.

**D. New alpha-user-facing guide.**
- `docs/kg_mcp_guide.md` rewritten for the lab-local flow (reach the box on the lab subnet → install MCP + research plugin → point at the lab Bolt URI with shared creds).
- Win11 and Linux instructions in parallel where they diverge (uv install, connectivity check).

### What is NOT changing

- The build / `neo4j-admin import` / post-import core pipeline. We append one block to post-import and run the *same* pipeline on the alpha stack from a tagged commit.
- `schema_config.yaml` — `Schema_info` is a BioCypher built-in (verified: `bc.write_schema_info(as_node=True)` in `create_knowledge_graph.py`; the node exists with `id: "schema_info"`).
- The explorer MCP's tool set or transport (still stdio).
- Dev workflow on `localhost:7687`. Adding the alpha stack doesn't touch the dev stack or its container names.

### Acceptance criteria

1. After a fresh build on **either** stack, `MATCH (s:Schema_info) RETURN s.version, s.built_at, s.git_sha, s.gene_count, s.expression_edge_count` returns populated values; dev reports `0.0.0-dev+<sha>`, alpha reports the tagged version.
2. `docker compose -p kg-dev up -d` and the alpha stack both stay running on the box for ≥1 hour without one OOM-killing the other (RAM budget §6.1).
3. From a **lab Win11 laptop**, `Test-NetConnection -ComputerName 132.75.249.47 -Port 17687` reports `TcpTestSucceeded : True`; from a host **outside** the lab subnet the same test reports `False` (firewall allowlist verified). (Linux operator: `nc -zv 132.75.249.47 17687`.)
4. `/release-kg 0.1.0-alpha.1` produces a tagged commit, a build **from a clean clone of that tag into the inactive color** (the live release keeps serving — zero build downtime), populated `Schema_info.version=0.1.0-alpha.1`, the shared `explorer` login provisioned, a GitHub Release with changelog + lab Bolt URI, and an operator checklist — in **≤60 min with warm caches**. Promotion is a single `alpha-deploy` restart; on failure the previous color keeps serving / is flipped back (no stranded testers).
5. After deploy, the shared `explorer` login can read; the explorer MCP rejects writes; and a connection attempt from outside the lab subnet is dropped by the firewall. (Read-only is MCP-enforced — Community has no roles, so every login is full-access; per-user / read-only-role accounts are deferred to the remote-MCP / Enterprise phase. See §2.5.)
6. A lab tester following `docs/kg_mcp_guide.md` can: reach the box on the lab subnet, install `uv`, point the explorer MCP at the lab URI with the shared credentials, install the `multiomics_research` skills plugin, run the compatibility check, and run one example from each section of the guide.

## 2. Design

### 2.1 Extend `Schema_info` with release metadata

The dev stack runs the post-import block with `KG_RELEASE_VERSION` unset (defaults to `0.0.0-dev`, `git_sha=unknown`); the release path sets `KG_RELEASE_VERSION=0.1.0-alpha.1` + the git vars. Same Cypher, two different stamps.

**Decision (a) (2026-05-25):** dev stamps **plain `0.0.0-dev`** — no `+<short_sha>` suffix — so bare `docker compose up -d` stays zero-friction (computing the SHA would need a wrapper or git inside the neo4j image, which it lacks). The dev stamp's only job is "not a release." Defaults are supplied by the **wrapper** (bash `${VAR:-default}` in `post-import.sh`; `:param` lines in `post-import.cypher`), so an unset *or empty* env var still yields the default — `coalesce()` alone would not, since it does not catch empty strings.

Properties added (all post-import; no build-time changes):

| Property | Type | Source |
|---|---|---|
| `version` | str | env `KG_RELEASE_VERSION`; defaults to `0.0.0-dev` |
| `built_at` | str (ISO 8601 UTC) | `toString(datetime())` in Cypher |
| `git_sha` | str | env `KG_GIT_SHA` (40-char) |
| `git_sha_short` | str | env `KG_GIT_SHA_SHORT` (7-char) |
| `git_branch` | str | env `KG_GIT_BRANCH` |
| `git_dirty` | str (`"true"`/`"false"`) | env `KG_GIT_DIRTY` (project convention: bools as strings) |
| `mcp_min_version` | str | env, default `"0.1.0"` |
| `release_notes_url` | str | env, optional |
| `paper_count` | int | `COUNT { (:Publication) }` |
| `experiment_count` | int | `COUNT { (:Experiment) }` |
| `gene_count` | int | `COUNT { (:Gene) }` |
| `organism_count` | int | `COUNT { (:OrganismTaxon) }` |
| `expression_edge_count` | int | `COUNT { ()-[:Changes_expression_of]->() }` |

(Counts are computed, not hardcoded — robust to data drift. Live graph today: 249,519 nodes, 2,188,444 rels, 99,871 genes, 232,758 expression edges, 197 experiments, 43 publications, 37 organisms.)

Append to `scripts/post-import.cypher` and the matching group in `scripts/post-import.sh` (both files must remain byte-identical in Cypher logic):

```cypher
MATCH (s:Schema_info {id: 'schema_info'})
SET s.version            = coalesce($version, '0.0.0-dev'),
    s.built_at           = toString(datetime()),
    s.git_sha            = coalesce($git_sha, 'unknown'),
    s.git_sha_short      = coalesce($git_sha_short, 'unknown'),
    s.git_branch         = coalesce($git_branch, 'unknown'),
    s.git_dirty          = coalesce($git_dirty, 'unknown'),
    s.mcp_min_version    = coalesce($mcp_min_version, '0.1.0'),
    s.release_notes_url  = coalesce($release_notes_url, '')
WITH s
SET s.paper_count           = COUNT { (:Publication) },
    s.experiment_count      = COUNT { (:Experiment) },
    s.gene_count            = COUNT { (:Gene) },
    s.organism_count        = COUNT { (:OrganismTaxon) },
    s.expression_edge_count = COUNT { ()-[:Changes_expression_of]->() };
```

**Run mode:** the existing post-import groups use quoted heredocs (`cypher-shell <<'CYPHER'`), which do not interpolate. Add this as a **separate, final `cypher-shell` invocation** that passes the env vars as params, reading them from the post-process container environment:

```bash
# Defaults are pushed through bash (${VAR:-default}) — the `:-` form fires on
# unset OR empty, so coalesce() in the Cypher is only a secondary guard.
cypher-shell \
  -P "version          => '${KG_RELEASE_VERSION:-0.0.0-dev}'" \
  -P "git_sha          => '${KG_GIT_SHA:-unknown}'" \
  -P "git_sha_short    => '${KG_GIT_SHA_SHORT:-unknown}'" \
  -P "git_branch       => '${KG_GIT_BRANCH:-unknown}'" \
  -P "git_dirty        => '${KG_GIT_DIRTY:-unknown}'" \
  -P "mcp_min_version  => '${KG_MCP_MIN_VERSION:-0.1.0}'" \
  -P "release_notes_url => '${KG_RELEASE_NOTES_URL:-}'" \
  <<'CYPHER'
... block above ...
CYPHER
```

### 2.2 Two-stack Docker layout

**Base `docker-compose.yml` edits (minimal):**
- **Keep** `container_name:` on all services (dev unaffected; `docker exec deploy …` and the docs that use it keep working).
- Parameterize the deploy port binding: `${KG_DEPLOY_HTTP_BIND:-127.0.0.1:7474}:7474` (same for `7687`), so the alpha override can rebind without forking the line. Dev defaults to loopback.
- `NEO4J_AUTH: ${NEO4J_AUTH:-none}` on `deploy` so the alpha override can flip it on; dev keeps `none`.

> **Review note B2:** the prior draft proposed *stripping* `container_name` so compose could auto-name per project. Rejected — `docker exec deploy cypher-shell …` and `docker logs import` are used throughout the dev workflow and ~7 plan/spec docs; auto-naming would rename `deploy` → `kg-dev-deploy-1` and break all of them. Overriding the names in the alpha file (below) gives two coexisting stacks **without** a repo-wide rename.

**Invocation:**
```bash
# Dev (unchanged)
docker compose -p kg-dev up -d

# Alpha
docker compose -p kg-alpha \
  -f docker-compose.yml \
  -f docker-compose.alpha.yml \
  --env-file .env.alpha \
  up -d
```

**Volumes auto-namespace** by project: `kg-dev_biocypher_neo4j_volume` vs `kg-alpha_biocypher_neo4j_volume` — fully isolated on disk. (The `container_name` override does not change volume namespacing.) **Blue/green note:** the alpha override (§2.4) must make the alpha data volume an **external, color-selected name** (e.g. `${ALPHA_DATA_VOLUME}` → `kg-alpha-blue` / `kg-alpha-green`) rather than the compose-managed default, so the release script can build the inactive color and flip `alpha-deploy` between them (§2.3 Phase 3/4). To wire in when §2.4 is implemented.

**Wrapper scripts** so the operator doesn't memorize the invocation:
```bash
#!/usr/bin/env bash
# scripts/alpha_up.sh
set -euo pipefail
cd "$(dirname "$0")/.."
exec docker compose -p kg-alpha \
  -f docker-compose.yml \
  -f docker-compose.alpha.yml \
  --env-file .env.alpha "$@"
```
Operator runs `./scripts/alpha_up.sh up -d` / `down` / `logs -f alpha-deploy`.

### 2.3 `/release-kg` skill — the release script

One script (`.claude/skills/release-kg/SKILL.md` + a Python helper) that does **all** of: GitHub bookkeeping, a **clean clone of the tag**, and the alpha Docker build + release. Runs on the box.

**Inputs:** `version` (required, e.g. `0.1.0-alpha.1`); `--mcp-min VERSION`; `--allow-dirty`; `--skip-rebuild` (Schema_info-only fixes against the live alpha stack); `--draft` (GitHub Release as draft); `--dry-run` (print every action, mutate nothing).

**Phase 1 — Preflight (fail fast, before any mutation):**
- Version matches `^\d+\.\d+\.\d+(-(alpha|beta|rc)\.\d+)?$`.
- Tooling ready: `docker`, `gh auth status` OK, `cypher-shell`, enough free disk for a clone + a fresh volume + output.
- Secrets present (gitignored, never logged): `.env` exists (the build needs it — C3); `.env.alpha` has `ALPHA_BIND_IP`, `NEO4J_AUTH` (bootstrap admin), `ALPHA_EXPLORER_PASSWORD` (shared read login).
- Git: on the release branch, in sync with origin, working tree clean (abort unless `--allow-dirty`). Capture SHA / short SHA / branch / dirty.

**Phase 2 — GitHub bookkeeping:**
- `CHANGELOG.md` (**Keep a Changelog — accumulate-then-cut**): the file is a living repo artifact. Notable changes are logged under a top `## [Unreleased]` section *as they land* (during each change's Document phase), **not** reconstructed from git log at release time. Phase 2 **cuts** that section: rename `## [Unreleased]` → `## [${version}] - <date>` and open a fresh empty `## [Unreleased]` above it. The operator polishes the human bits; the stats line (gene/edge/experiment counts) is auto-filled in Phase 4 from the **built graph** (C5). The GitHub Release notes (Phase 5) are **extracted from this one version's section** — the changelog is the source of truth, the Release is a rendering of it. (`release-please`-style auto-generation stays out of scope — see §2.3 close.)
- Commit `chore(release): kg-${version}`; annotated tag `kg-${version}`.
- Push commit + tag: `git push origin <branch> --follow-tags`. (Tag pushed now so the Phase-3 clean clone pulls exactly what is released.)

**Phase 3 — Clean clone + build into the *inactive* color (zero-downtime build):**
- **Clean clone of the tag** into a throwaway dir: `git clone --branch kg-${version} --depth 1 <origin> /tmp/kg-release-${version}`. The build reads the *tag*, never the live working tree (the C2 fix) — dev work continues and the release is reproducible from the tag. *Determinism caveat:* the build still consumes committed `cache/` + the warm `build_cache` volume — reproducible-given-cache, not hermetic. Don't `--no-cache` a routine alpha (~1 h+ re-download).
- Copy the gitignored `.env` + `.env.alpha` into the clone; write the release-stamp vars (`KG_RELEASE_VERSION`, `KG_GIT_SHA[_SHORT]`, `KG_GIT_BRANCH`, `KG_GIT_DIRTY`) into the clone's `.env.alpha`.
- **Blue/green volumes.** Two alpha data volumes alternate per release (`…-blue` / `…-green`); the live `alpha-deploy` mounts the *active* color, tracked in a small state marker (e.g. `.alpha_active_color`). The release builds into the **inactive** color as a separate transient project (`kg-alpha-staging`) bound to a **temporary localhost port** → `build` → `import` → `post-process` (Schema_info stamp). **The active color keeps serving throughout** — the long build/import/post-process incurs *zero* downtime.

**Phase 4 — Verify the staged color, then flip (minimal downtime):**
- Against the staged color (temp localhost port): sanity-check Schema_info (`version` matches the tag, counts roughly match dev); fill the changelog stats from these.
- **Smoke test via the real explorer.** Point a throwaway explorer instance (the currently-pinned compatible version) at the staged color and run the compatibility-check query **and one starter query** — this exercises the actual KG↔explorer contract, catching schema/tool mismatches a bare Cypher check would miss.
- **Flip.** Stop the active `alpha-deploy`; bring the staged color up on the real `${ALPHA_BIND_IP}:17687` ports (auth on) and update the active-color marker. **Downtime = one deploy restart (≈ under a minute)**, not the build duration.
- **Provision the shared read login** on the now-live color: `CREATE USER explorer IF NOT EXISTS SET PASSWORD '${ALPHA_EXPLORER_PASSWORD}' CHANGE NOT REQUIRED;` (re-run each release — fresh volume). `neo4j` bootstrap admin stays operator-only. Confirm the `DOCKER-USER` allowlist is present on the live ports.
- **On failure at any point:** the previous color is untouched and still serving (or flip back to it). Never leave testers stranded; the retired color is overwritten only on the *next* release.

**Phase 5 — Publish + handoff:**
- `gh release create kg-${version} --notes-file <changelog-fragment> --prerelease` (+ `--draft` if asked). Notes = changelog section + a "Connecting" block (lab Bolt URI + `docs/kg_mcp_guide.md` link). Attach a small `metadata.json` manifest (version, sha, counts, timestamps) for provenance — **no graph dump** (the box IS the distribution).
- Print the operator checklist + tester-announcement template: lab Bolt URI, `.mcp.json` snippet, MCP guide link, GitHub Release URL, reminder that the shared `explorer` password goes out-of-band, confirmation the firewall allowlist is active. Delete the throwaway clone.

**Cross-repo coordination (decided):** KG and explorer **release on independent cycles** — not every KG change needs an explorer change or vice versa — so `/release-kg` does **not** tag the sibling repos. The contract between them is `Schema_info.mcp_min_version` (the `--mcp-min` input): the KG release **declares and records the known-compatible explorer version** in its notes, and the Phase-4 smoke test runs the actual pinned explorer against the new KG to prove it. `multiomics_research` is a **tester-facing fork** — testers fork it and pull upstream for new versions; it is versioned on its own, not by `/release-kg`. When a change genuinely spans both repos, cut an explorer release first, bump `--mcp-min`, then cut the KG release.

**Out of scope for v1 (defer):** `release-please`-style auto-CHANGELOG; automated cross-repo (explorer/research) tagging. (Blue/green build-then-flip + prev-color rollback is **in** scope — Phase 3/4.)

### 2.4 `docker-compose.alpha.yml` (sketch)

```yaml
# Override: alpha stack (lab testers, read-only via MCP, auth on, bound to the lab IP).
# Invoked as: docker compose -p kg-alpha -f docker-compose.yml -f docker-compose.alpha.yml --env-file .env.alpha up -d
services:

  build:
    container_name: alpha-build
    # NOTE: build.sh hardcodes `cp /src/.env .` and builds from the live working
    # tree — so the alpha build uses the repo's .env (dev secrets, fine) and whatever
    # is checked out. The release flow checks out the tag first (C2). No release env
    # vars here: the Schema_info stamp happens in post-process, not build.

  import:
    container_name: alpha-import
    volumes:
      - ./output-alpha:/output   # isolate import.status / import.report from dev's ./output

  post-process:
    container_name: alpha-post-process
    environment:
      NEO4J_dbms_security_auth__enabled: "false"   # auth enabled only on `deploy`
      KG_RELEASE_VERSION: ${KG_RELEASE_VERSION:?must be set in .env.alpha}
      KG_GIT_SHA:         ${KG_GIT_SHA:?}
      KG_GIT_SHA_SHORT:   ${KG_GIT_SHA_SHORT:?}
      KG_GIT_BRANCH:      ${KG_GIT_BRANCH:?}
      KG_GIT_DIRTY:       ${KG_GIT_DIRTY:-false}
      KG_MCP_MIN_VERSION: ${KG_MCP_MIN_VERSION:-0.1.0}

  deploy:
    container_name: alpha-deploy
    ports:
      # Bind to the box's lab IP only — never 0.0.0.0. ${ALPHA_BIND_IP} from .env.alpha.
      - "${ALPHA_BIND_IP}:17474:7474"
      - "${ALPHA_BIND_IP}:17687:7687"
    environment:
      NEO4J_AUTH: ${NEO4J_AUTH:?e.g. neo4j/<strong-shared-password>}
      NEO4J_dbms_security_auth__enabled: "true"
      # Smaller heap than dev — alpha is read-only, low-concurrency
      NEO4J_server_memory_heap_initial__size: "512M"
      NEO4J_server_memory_heap_max__size: "1G"

  # Drop the Biochatter UI; testers use their own MCP clients
  app:
    profiles: ["never"]
```

`.env.alpha.example` (committed; real `.env.alpha` is gitignored — see C4):
```bash
# Versioning — set by /release-kg, do not edit manually
KG_RELEASE_VERSION=0.1.0-alpha.1
KG_GIT_SHA=
KG_GIT_SHA_SHORT=
KG_GIT_BRANCH=main
KG_GIT_DIRTY=false
KG_MCP_MIN_VERSION=0.1.0

# Operator-set, persistent
ALPHA_BIND_IP=132.75.249.47                       # the box's lab IP — request a DHCP reservation / static IP from IT (see §6.4)
NEO4J_AUTH=neo4j/<strong-admin-password>          # bootstrap admin (neo4j) password, operator-only. Never commit.
ALPHA_EXPLORER_PASSWORD=<strong-shared-password>  # shared `explorer` read login; provisioned each release (§2.3 Phase 4), distributed out-of-band. Never commit.
```

### 2.5 Access control (shared `explorer` login + MCP read-only + firewall)

> **Capability vs. choice (verified empirically on `neo4j:5.15-community`):** Community *can* create multiple named users (`CREATE USER alice` succeeds; `SHOW USERS` lists it) but supports **no roles or privileges** (`SHOW ROLES`, `GRANT ROLE`, `GRANT/DENY … ON GRAPH/DATABASE` all return *"Unsupported administration command"*). So every login is **full read + write + schema** — you cannot make a read-only one. *Capable* of per-user logins ≠ *worth it* here; see the ROI note below.

The alpha access model is three layers:

1. **One shared `explorer` login.** A single read credential (`CREATE USER explorer …`, re-provisioned each release by the release script — §2.3 Phase 4), distributed out-of-band, used by all testers. It is **not** the `neo4j` bootstrap admin (operator-only), so a leak rotates `explorer`, not admin.
2. **Explorer MCP enforces read-only** — its `run_cypher` tool blocks writes. The MCP connects **straight to `NEO4J_URI` over Bolt** (no separate tunnel; stdio is only the local MCP↔client transport) using the shared `explorer` credential in each tester's `.mcp.json`. Read-only is enforced in MCP code, **not** by a DB role.
3. **Host-firewall allowlist** (§2.6) pins the ports to the lab subnet.

**Why shared, not per-user (ROI):** per-user Neo4j accounts buy ~nothing here. *No* security gain (all Community logins are full-access; read-only is the MCP regardless); revocation is already free because the release **wipes the volume every release** (§2.3) so accounts are re-created each time anyway; and attribution is better captured by `multiomics_research`'s MCP-side usage logging than by DB-level "who connected." Net: per-user = recurring re-provisioning toil + N secrets to manage, for marginal benefit. A single shared `explorer` login is the pragmatic and **customary** choice for a local **stdio** MCP — stdio servers have no per-user auth layer; a shared service credential to the backing store is the norm. Per-user *identity* is a remote-hosted-MCP concern (§7.2/§7.4).

**Honest limitation:** any holder of the shared credential who bypasses the MCP (raw driver / cypher-shell) **can write** — Community has no read-only enforcement. Acceptable for a physically-present, trusted lab alpha. **Real read-but-not-write and per-user identity arrive together** at the remote-MCP-hosting (OAuth) / Enterprise / Aura phase — deferred to §7. Acceptance criterion #5 reflects this.

### 2.6 Host-firewall allowlist (required)

`132.75.249.47` is a **campus-routable** IP, so binding the alpha ports to it exposes them to anything on campus that can route there — not just the lab. The allowlist is **not optional**.

**Gotcha:** Docker publishes ports by inserting rules into the `DOCKER`/`DOCKER-USER` iptables chains, which are consulted **before** ufw's `INPUT` rules — so a plain `ufw allow/deny` does **not** filter Docker-published ports. Filter in the `DOCKER-USER` chain instead:

```bash
# Drop traffic to the alpha ports from anything outside the lab subnet.
sudo iptables -I DOCKER-USER -p tcp --dport 17687 ! -s <LAB_SUBNET_CIDR> -j DROP
sudo iptables -I DOCKER-USER -p tcp --dport 17474 ! -s <LAB_SUBNET_CIDR> -j DROP
# Persist across reboots (Debian/Ubuntu):
sudo apt-get install -y iptables-persistent && sudo netfilter-persistent save
```
- Confirm `<LAB_SUBNET_CIDR>` with IT (e.g. `132.75.249.0/24`). **Do not** allowlist the whole campus `/16`.
- The `ufw-docker` helper is an alternative to hand-written `DOCKER-USER` rules.
- **To add home-via-VPN access later:** add the VPN-pool CIDR with another `DOCKER-USER` allow/insert rule. No other change. (This is the cheap, reversible expansion path.)

### 2.7 Alpha-user-facing guide (`docs/kg_mcp_guide.md`)

Section structure (lab-local flow):
1. **What this is.**
2. **Reach the box** — you must be on the lab subnet (in the lab, or later on the university VPN). Connectivity check, Win11 + Linux side-by-side:
   - Win11: `Test-NetConnection -ComputerName 132.75.249.47 -Port 17687` (run in **PowerShell**, not Command Prompt — look for `TcpTestSucceeded : True`).
   - Linux: `nc -zv 132.75.249.47 17687`.
3. **Install the explorer MCP** — `uvx --from git+…`.
4. **Fork & install the research skills plugin** — fork `multiomics_research`, clone your fork (`uv sync` pulls the explorer), register as a Claude Code plugin per its README; pull upstream for new versions.
5. **Configure your MCP client** — `NEO4J_URI: bolt://132.75.249.47:17687`, username/password from the operator. `bolt://` (no `+s`) — trusted lab subnet.
6. **First query — compatibility check** — matches against `Schema_info.version` / `mcp_min_version`.
7. **Find your way around / Starter queries / Conventions and gotchas.**
8. **Limits and support** — lab-local only for now (must be on the lab subnet); how to file issues; dev iteration cadence.
9. **Updating** — releases on GitHub; pull new explorer/research; rerun the compatibility check.

## 3. Test

- Unit-test the version-regex in the release skill helper.
- `pytest -m kg` additions (against `bolt://localhost:7687`, dev stack):
  - `test_schema_info_release_properties` — `Schema_info` has `version`, `built_at`, `git_sha`, `gene_count`, `expression_edge_count` populated and right-typed.
  - `test_schema_info_counts_match` — `s.gene_count == COUNT { (:Gene) }` etc., catches stale stats.
- Manual:
  - Bring up both stacks. `docker ps` shows `deploy` and `alpha-deploy`, ports differ, volumes differ. Probe both: dev → `0.0.0-dev+…`, alpha → tagged version.
  - From a lab laptop: connectivity check (criterion #3), MCP compatibility check, one tool per starter-query category.
  - From an off-subnet host: connectivity check **fails** (firewall verified).
  - Run `/release-kg 0.0.0-dev.0` end-to-end against the alpha checkout; verify the alpha stack comes back up with the new `Schema_info` properties.

## 4. Review

- `post-import.cypher` and `post-import.sh` must remain byte-identical in Cypher logic — run `scripts/post-import-validate.sh` before/after.
- `Schema_info` is BioCypher-owned; do not duplicate with a separate `KGRelease` label.
- Alpha port bindings must use `${ALPHA_BIND_IP}` (the lab IP), never `0.0.0.0`; the `DOCKER-USER` allowlist must be present and persisted — grep `docker-compose.alpha.yml` and confirm `iptables -L DOCKER-USER`.
- `.env.alpha` must NOT be committed — add the explicit pattern to `.gitignore` (the existing `.env` line does not match `.env.alpha`; see C4). Only `.env.alpha.example` is committed.
- Secrets discipline: `NEO4J_AUTH` (admin) and `ALPHA_EXPLORER_PASSWORD` (shared read login) come only from `.env.alpha`; the release script must never echo them to logs or the GitHub Release notes. The shared `explorer` login is re-provisioned every release (the volume wipe drops it).
- The alpha build reads the repo `.env` (build.sh `cp /src/.env .`); ensure a `.env` exists on the box or the build aborts under `set -e` (C3).

## 5. Document

- This plan: mark workstreams complete as they ship.
- `CLAUDE.md`: add a "Two-stack deploy" subsection under Docker Pipeline Stages (dev vs alpha split, container-name convention, lab-IP bind + firewall), and a brief "Releases" pointer to the skill + `Schema_info` properties.
- `docs/kg_mcp_guide.md`: ship rewritten alongside the first release.
- Memory: short reference memory pointing at this plan + the MCP guide.

## 6. Sizing, corrections, and open questions

### 6.1 Box resource budget

Two Neo4j containers (steady state = two `deploy` containers; `build`/`import`/`post-process` are transient during a rebuild):

| | Dev `deploy` | Alpha `alpha-deploy` | Notes |
|---|---:|---:|---|
| Heap (max) | 2 GB | 1 GB | Smaller on alpha — read-only, low-concurrency |
| Page cache | ~1 GB | ~1 GB | Neo4j default heuristic |
| Volume disk (DB + build2neo CSVs) | ~3–5 GB | ~3–5 GB | Same data, separate project volumes |
| Import artifacts (`./output*`) | ~KB | ~KB | **Only** `import.status` + `import.report` live here — *not* the CSVs (those are in the volume; see C1) |

**This box (measured 2026-05-24): 125 GiB RAM, 743 GB free disk (`/dev/sda5`, 16% used); each KG data volume is 3.14 GB.** Keeping three KGs (dev + both alpha colors) ≈ **9.4 GB of data volumes** (~18–22 GB with per-project build caches + shared images + transient build layers), and peak RAM during a release (staged `post-process` 2G/4G heap *while* active `alpha-deploy` 1G + dev `deploy` 2G keep serving) is **<15 GB heap + page caches**. That's **~3% of free disk and ~10% of RAM** — ample headroom; the generic "16 GB-box, tight" caveat does not apply to this machine. Pause dev rebuilds during a release only to avoid two full CPU-bound builds competing, not for memory.

### 6.2 Corrections folded in from review

- **C1** — `./output` holds only the tiny `import.status`/`import.report` (per `scripts/import.sh`); the build CSVs go to the project-namespaced volume at `/data/build2neo`. The `./output-alpha` override is still correct and necessary (otherwise both stacks race on `./output/import.report`, which the `test_import_report` KG-validity test reads).
- **C2** — "the tag IS the release" needs care: the base `build` bind-mounts the live working tree (`.:/src/`) and `build.sh` does `cp -r /src/* .`, plus the build reads `cache/` + `build_cache` (not hermetic). Resolved: the release script builds from a **clean clone of the tag** (§2.3 Phase 3), so the live working tree is never the build input; reproducibility is "clean clone of the tag + warm cache," not fully hermetic.
- **C3** — `build.sh` hardcodes `cp /src/.env .`, so the alpha build consumes the repo `.env`, and release env vars on the `build` service would be dead weight (the stamp is in post-process). Dropped them from the build service; documented the `.env`-must-exist dependency.
- **C4** — `.gitignore` has an exact `.env` line, which does **not** match `.env.alpha`. Add `.env.alpha` explicitly.
- **C5** — changelog stats are pulled from the **alpha** graph post-rebuild (§2.3 step 9), not the dev KG, so the numbers match what ships.
- **B1 / B2** — see §2.5 and §2.2.
- `nc` is not native on Win11 → connectivity check uses `Test-NetConnection` (PowerShell) for testers; `nc` only for the Linux operator.

### 6.3 Release mechanics: blue/green build-then-flip

The alpha rebuilds from a clean clone of the tag (~30–60 min warm) into the **inactive** color volume while the active color keeps serving, then flips (§2.3 Phase 3/4). Downtime is one `alpha-deploy` restart (≈ under a minute), and the previous color is the rollback target (flip back). Costs: two alpha data volumes (~2× alpha disk) and build-window RAM pressure (§6.1). A pure dev→alpha volume-copy path (skip the alpha rebuild entirely: `docker run --rm -v <dev_vol>:/from -v <alpha_inactive_vol>:/to alpine cp -a /from/. /to/`, then re-stamp `Schema_info` and flip) remains a future option if per-release rebuild ever becomes painful.

### 6.4 Open questions

1. **Lab subnet CIDR + box IP stability.** Confirm `<LAB_SUBNET_CIDR>` with IT, and pin the box IP: `132.75.249.47` is likely DHCP — request a **DHCP reservation / static IP** (or a lab DNS name testers can use), otherwise a reboot changes the URI and breaks every tester. *(Pre-flight before inviting testers.)*
2. **Auth init across the auth-disabled→enabled volume handoff.** Import + post-process run with auth disabled; `alpha-deploy` then starts with auth **on** + `NEO4J_AUTH`. This should set the initial `neo4j` password on first auth-enabled start — confirm on the first alpha cut (criterion #5 covers it).
3. **`DOCKER-USER` allowlist actually filters.** Verify from an off-subnet host that `17687` is dropped (Docker-bypasses-ufw gotcha, §2.6).
4. **Win11 uv install ergonomics.** Document `winget install astral-sh.uv` + the PowerShell installer; confirm on a clean Win11 box before inviting testers. Note: connectivity check runs in **PowerShell**, not cmd.
5. **MCP version pinning.** `Schema_info.mcp_min_version` is the contract; needs a real version-bump cadence on `multiomics_explorer` and a guide pointer to the right tag.
6. **Box reliability.** Single point of failure (box uptime, lab power, the dev stack hogging RAM mid-release). Document a "alpha box is down" channel for testers. Not a blocker for the first cut.

## 7. Deferred: remote & multi-user access (workstream D)

Revisit when testers go beyond the lab subnet, when per-user database accounts are needed, or when box reliability becomes an issue.

### 7.1 Remote reach beyond the lab

- **University VPN (cheapest):** testers at home on the university VPN land back on the campus network. Add the VPN-pool CIDR to the `DOCKER-USER` allowlist (§2.6) — one rule, no rearchitecting. *This is the first expansion step and needs no new SaaS.* Verify the VPN actually routes to the box's subnet (some university VPNs are split-tunnel and only route specific resources).
- **Tailscale:** an encrypted overlay reaching the box from anywhere. Outbound-only, so it traverses the university firewall without an inbound exception. **Costs (corrected):** the free Personal plan now allows **6 users / unlimited devices but is non-commercial-use** — a lab serving data to colleagues is plausibly outside that; the clean paid self-serve tier is **Standard at $8/user/mo (~$576/yr for 6 users)** — note "Starter ~$6/user" is a **legacy** tier, not current. Redundant with the university VPN for our tester population; only worth it for wire encryption + tight (tailnet-only) exposure if we ever distrust the campus subnet.
- **Cloudflare Tunnel:** also outbound-only; free tier covers it; testers tunnel Bolt via `cloudflared` (extra per-tester step). Viable VPN-alternative fallback.
- **Public IP + Let's Encrypt:** **not viable behind the university firewall** (no inbound, likely no public routable port). Struck.

### 7.2 Per-user identity & read-only enforcement

Two things the alpha deliberately skips, which arrive together later:
- **Per-user identity / auth** is, for a local **stdio** MCP, not a layer that exists — the server runs under the launching user; a shared backing-store credential is customary. It becomes real when the **MCP is remote-hosted** (§7.4): the MCP spec's OAuth 2.1 authorization gates per-user access at the *MCP* layer, which is its right home — not per-user Neo4j accounts.
- **Read-but-not-write at the DB** requires the `reader` role / privilege grants — **Neo4j Enterprise** or **Aura** only. Until then read-only is MCP-enforced (§2.5).

So per-user identity rides with remote-MCP-hosting + Enterprise/Aura, not with the lab-only alpha.

### 7.3 Aura migration (hosting comparison / PI brief)

**KG sizing (measured 2026-05-24):** 249,519 nodes, 2,188,444 relationships, 3.1 GB on-disk Neo4j volume.

| # | Option | Year-1 cost | Reach | Reliability burden | Iteration |
|---|---|---:|---|---|---|
| **1** | **Local box, lab-only** (this plan) | **$0** | Lab subnet (+ univ-VPN via allowlist) | Operator (box uptime) | Minutes |
| 2 | Aura Free | $0 | Public, anywhere | None | Slow (manual dump) |
| 3 | Aura Pro 2 GB / 4 GB | ~$1,577 always-on (~$690 w/ pause) | Public, anywhere | None | Slow |
| **4** | **Aura Pro 4 GB / 8 GB** | **~$3,154 always-on (~$1,380 w/ pause)** | Public, anywhere | None | Slow |
| 5 | Aura Pro 8 GB / 16 GB | ~$6,307 always-on | Public, anywhere | None | Slow |

- **Option 2 (Aura Free) is not viable:** caps at 200K nodes / 400K rels; we're at 249K / 2.2M.
- **Option 4 is the honest Aura tier** for this KG (comfortable cache; room to ~2×).

**Migrate from local-box (Option 1) to Aura (Option 4) when any of:**
1. Testers need access beyond the campus network / university VPN (external collaborators).
2. Box uptime becomes a recurring complaint.
3. Per-user accounts / audit become a requirement (and Enterprise-on-box isn't pursued).
4. We want a public-launch posture ("here's the URI, no VPN required").
5. University IT forbids serving the DB on the campus network / over the VPN.

**Budget framing:** Option 1 is $0; Option 4 is ~$1–3K/yr. Not worth it for a lab-local alpha; likely worth it at beta. **PI ask:** approve Option 1 for the alpha; heads-up on the eventual Option 4 line item if a migration trigger hits.

### 7.4 Remote-host the explorer MCP

Today the MCP is stdio — each tester runs it locally (clone repos + `uv sync`). Hosting it as an HTTP/SSE service (Fly.io / Modal / a VM) is its own decision: it is the natural home for **per-user auth** (OAuth 2.1 per the MCP spec, §7.2) and removes the per-tester install. Deferred; revisit alongside Aura.
