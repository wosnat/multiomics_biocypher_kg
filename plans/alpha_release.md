# Plan: Alpha release on local Docker via Tailscale

**Status:** Draft (rewritten 2026-05-24, supersedes the Aura-targeted draft of 2026-05-14)
**Goal:** Cut tagged, versioned alpha releases of the KG served from a *second* Neo4j stack on the same Linux dev box, so ~5 alpha testers (each with their own readonly Neo4j account, reaching the box over Tailscale) can drive the graph from their own LLM agents via the explorer MCP + research-skills plugin.

## Why not Aura

The prior draft targeted Neo4j Aura Professional. Decision reversed 2026-05-24 because:
- **Cost.** Aura Pro is the smallest tier that fits our graph (250K nodes / 2.2M relationships > Aura Free's 200K/400K cap) and runs ~$65+/month idle. Local box is $0 marginal.
- **Iteration speed.** Aura release = `neo4j-admin database dump` → manual upload via Aura Console → wait for restore. Local release = `docker compose up -d` on a tagged commit. Minutes vs. tens of minutes.
- **Scope (5 testers, LAN).** All current testers are on the same LAN as the box. Aura's "reach from anywhere" advantage doesn't apply yet. Tailscale gives the same reach (over an encrypted overlay) without the dump/restore overhead, and includes the off-LAN escape hatch if a tester travels.
- **Two-deploy model.** The dev iteration loop already runs `docker compose up -d` on `localhost:7687`. Adding a parallel "alpha" stack on the same box for testers (different ports, different volume, auth on) is cheaper than maintaining a dev↔Aura sync.

Aura remains a viable phase-2 path (workstream D) for a broader public release; deferred.

## Workstreams

| ID | Title | Scope |
|---|---|---|
| **A** | Two-stack local Docker (dev + alpha side-by-side) | **in scope** |
| **B** | Tailscale + per-user Neo4j readonly accounts | **in scope** |
| **C** | Alpha-user-facing guide (two-repo install: `multiomics_explorer` + `multiomics_research`) | **in scope** |
| **D** | Aura migration | **deferred** — revisit when testers go beyond the Tailnet or the box reliability becomes an issue |
| **E** | Remote-host the explorer MCP | **deferred** — same reason as before; FastMCP via Fly.io/Modal/VM is its own decision |

Until D and E land: each alpha user installs the explorer MCP and the research-skills plugin locally, joins the Tailnet, and points the MCP at `bolt://kg-box.<tailnet>.ts.net:17687` with their personal readonly credentials.

## 1. Scope

### What is changing

**A. Two-stack Docker layout.**
- Strip `container_name:` from `docker-compose.yml` so compose can run multiple project-scoped copies of the same service.
- New `docker-compose.alpha.yml` override: different ports (`17474`/`17687`), separate output mount (`./output-alpha`), auth turned on, no `app` (Biochatter) service.
- New `.env.alpha.example` with `KG_RELEASE_VERSION`, `KG_GIT_SHA`, `KG_GIT_BRANCH`, `NEO4J_AUTH` template.
- New `scripts/alpha_up.sh` / `scripts/alpha_down.sh` thin wrappers around the `-p kg-alpha -f … -f …` invocations so the operator doesn't have to remember the full command.

**B. Tailscale + auth.**
- Operator-side runbook (in the alpha guide): install Tailscale on the box, enable MagicDNS, expose the box as `kg-box`, invite 5 testers to the tailnet.
- Bind alpha ports to the Tailscale IP only (`100.x.y.z:17687:7687`), never `0.0.0.0` — defense in depth on top of Neo4j auth.
- Optional Tailscale ACL pinning testers to the box's port `17687` only (recommended; sample ACL in the guide).
- Neo4j auth turned ON in the alpha stack (`NEO4J_AUTH=neo4j/<strong>`, `NEO4J_dbms_security_auth__enabled: true`).
- Post-deploy Cypher template for per-user `reader` accounts (Neo4j built-in role) against the `system` database.

**C. Existing pieces retained.**
- `Schema_info` release-metadata properties (see 2.1) — still useful, decouples build-time vs. release-time identity.
- `CHANGELOG.md` at repo root (Keep a Changelog format).
- `/release-kg` skill — reshaped around the local two-stack model (see 2.3).
- Tag scheme `kg-X.Y.Z[-(alpha|beta|rc).N]`.
- GitHub Release per tag, marked pre-release for any non-stable suffix.

**D. New alpha-user-facing guide.**
- `docs/kg_mcp_guide.md` rewritten for the Tailscale + local-box flow.
- New section: install `multiomics_research` plugin in addition to the explorer MCP.
- Win11 and Linux instructions in parallel where they diverge (uv install, Tailscale install).

### What is NOT changing

- The build / `neo4j-admin import` / post-import core pipeline. We append one block at the end of post-import and we run the *same* pipeline on the alpha stack from a tagged commit.
- `schema_config.yaml` — `Schema_info` is a BioCypher built-in.
- The explorer MCP's tool set or transport (still stdio).
- Dev workflow on `localhost:7687`. Adding the alpha stack doesn't touch the dev stack.

### Acceptance criteria

1. After a fresh build on **either** stack, `MATCH (s:Schema_info) RETURN s.version, s.built_at, s.git_sha, s.gene_count, s.expression_edge_count` returns populated, non-empty values; the dev stack reports `0.0.0-dev+<sha>` and the alpha stack reports the tagged version.
2. `docker compose -p kg-dev up -d` and `docker compose -p kg-alpha -f docker-compose.yml -f docker-compose.alpha.yml up -d` both stay running on the same box for ≥1 hour without one OOM-killing the other (RAM budget verified in §6).
3. From a Win11 tester laptop, joined to the tailnet, `nslookup kg-box.<tailnet>.ts.net` resolves and `nc -zv kg-box.<tailnet>.ts.net 17687` succeeds.
4. `/release-kg 0.1.0-alpha.1` produces a tagged commit, an alpha-stack rebuild against that tag, populated `Schema_info.version=0.1.0-alpha.1`, a GitHub Release with changelog + Tailscale Bolt URI, and a printed operator checklist — in under 60 min including the alpha-stack rebuild.
5. After running the per-user Cypher template, the operator can verify each of 5 users can read but not write (`CALL dbms.security.listUsers()` shows 5 + neo4j; a `CREATE (:Test)` from a tester account is rejected).
6. An alpha tester following `docs/kg_mcp_guide.md` can: install Tailscale and join the tailnet, install `uv`, point the explorer MCP at the alpha URI with their credentials, install the `multiomics_research` skills plugin, run the compatibility check, and run one example from each section of the guide.

## 2. Design

### 2.1 Extend `Schema_info` with release metadata

Unchanged from the prior draft — the property set, the post-import block, and the env-var plumbing all carry over verbatim. The only difference is *who calls it*: the dev stack runs it with `KG_RELEASE_VERSION` unset (defaults to `0.0.0-dev+<short_sha>`), and the alpha stack runs it with `KG_RELEASE_VERSION=0.1.0-alpha.1` set in `.env.alpha`. Same Cypher, two different stamps.

Properties added (all post-import; no build-time changes):

| Property | Type | Source |
|---|---|---|
| `version` | str | env `KG_RELEASE_VERSION`; defaults to `0.0.0-dev+<short_sha>` |
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
MATCH (s)
SET s.paper_count           = COUNT { (:Publication) },
    s.experiment_count      = COUNT { (:Experiment) },
    s.gene_count            = COUNT { (:Gene) },
    s.organism_count        = COUNT { (:OrganismTaxon) },
    s.expression_edge_count = COUNT { ()-[:Changes_expression_of]->() };
```

Run mode: pass the env vars through to `cypher-shell -P "version=>'$KG_RELEASE_VERSION'" …` inside the post-process container.

### 2.2 Two-stack Docker layout

**Base `docker-compose.yml` edits:**
- Remove all `container_name:` lines. Compose will auto-name containers `<project>-<service>-<n>`, which lets two projects coexist. Nothing in our scripts shells out to a container by hard-coded name.
- Switch the deploy port binding from `127.0.0.1:7474:7474` to `${KG_DEPLOY_HTTP_BIND:-127.0.0.1:7474}:7474` (same for `7687`) so the alpha override can rebind to the Tailscale IP without forking the whole port line.
- Move `NEO4J_AUTH` to env-var-with-default: `NEO4J_AUTH: ${NEO4J_AUTH:-none}`. Dev keeps `none`; alpha overrides.

**New `docker-compose.alpha.yml`:** see §2.4 — full file is sketched separately so it can be reviewed as a diff.

**Invocation:**
```bash
# Dev (no change in behavior)
docker compose -p kg-dev up -d

# Alpha
docker compose -p kg-alpha \
  -f docker-compose.yml \
  -f docker-compose.alpha.yml \
  --env-file .env.alpha \
  up -d
```

**Volumes auto-namespace:** `kg-dev_biocypher_neo4j_volume` vs `kg-alpha_biocypher_neo4j_volume`. The two stacks are fully isolated on disk.

**Networks auto-namespace** the same way — `kg-dev_biochatter` vs `kg-alpha_biochatter` (alpha drops the `app` service so the network is unused; can be removed in the override).

**Wrapper scripts** (`scripts/alpha_up.sh`, `scripts/alpha_down.sh`) so the operator doesn't have to remember the full invocation:
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
Operator runs `./scripts/alpha_up.sh up -d` / `down` / `logs -f deploy` / etc.

### 2.3 `/release-kg` skill (rewritten for local two-stack)

Layout: `.claude/skills/release-kg/SKILL.md` + small Python helper.

**Inputs:** `version` (required, e.g. `0.1.0-alpha.1`), `--mcp-min VERSION`, `--allow-dirty`, `--skip-rebuild` (operate against the already-built alpha stack — for `Schema_info`-only fixes), `--draft` (publish GitHub Release as draft).

**Steps:**
1. Validate version matches `^\d+\.\d+\.\d+(-(alpha|beta|rc)\.\d+)?$`.
2. Capture git state: SHA, short SHA, branch, dirty. Abort if dirty unless `--allow-dirty`.
3. Open `CHANGELOG.md` in `$EDITOR` pre-filled with a template stub for the new version, auto-filled with stats queried from the **dev** KG (`localhost:7687`). Operator fills in the human bits and saves.
4. Commit: `chore(release): kg-${version}`. Tag: `git tag -a kg-${version} -m "alpha release ${version}"`.
5. Update `.env.alpha`: set `KG_RELEASE_VERSION=${version}`, `KG_GIT_SHA=$(git rev-parse HEAD)`, etc. (idempotent overwrite of the release-stamp lines).
6. **Bring down the alpha stack first** — `./scripts/alpha_up.sh down` — because alpha's `deploy` holds a lock on its volume, and we're about to re-run import.
7. Wipe the alpha volume (`docker volume rm kg-alpha_biocypher_neo4j_volume`) — clean slate per release. Alpha is read-only by design; nothing of value lives there outside the `Schema_info` stamp.
8. Unless `--skip-rebuild`: `./scripts/alpha_up.sh up -d` → triggers `build` → `import` → `post-process` → `deploy` on the alpha side from the tagged commit. Stream the post-process logs until the `Schema_info` block runs.
9. Sanity-check via cypher-shell against `bolt://localhost:17687`: `MATCH (s:Schema_info) RETURN s.version, s.gene_count, s.expression_edge_count`. Must match the tag and roughly match dev counts.
10. Provision/refresh the per-user readonly accounts (§2.5) — the skill prompts whether to re-run.
11. `git push origin main --follow-tags`.
12. Extract the new changelog section, append a "Connecting" block with the Tailscale URI + a pointer to `docs/kg_mcp_guide.md`, then `gh release create kg-${version} --title "KG ${version}" --notes-file <changelog-fragment> --prerelease`. No dump attached — the alpha stack on the box IS the distribution.
13. Print the operator checklist: tailnet URI, MCP guide path, GitHub Release URL, sample `.mcp.json` snippet pre-filled with the URI for testers, reminder to distribute per-user passwords out-of-band.

**Out of scope for v1 (defer):**
- Build-once + volume-copy release path (dev→alpha) as an alternative to step 7-8. Listed in §6 as a future option if the per-release rebuild ever becomes painful.
- Snapshotting the alpha volume before destroy. We can always rebuild from the tag.
- `release-please`-style auto-CHANGELOG. Overkill for a single maintainer.

### 2.4 `docker-compose.alpha.yml` (sketch)

```yaml
# Override: alpha stack (testers, read-only, auth on, Tailscale-only)
# Invoked as: docker compose -p kg-alpha -f docker-compose.yml -f docker-compose.alpha.yml --env-file .env.alpha up -d
services:

  build:
    # Pin to the tagged commit; same image but reads release env vars
    environment:
      KG_RELEASE_VERSION: ${KG_RELEASE_VERSION:?must be set in .env.alpha}
      KG_GIT_SHA:         ${KG_GIT_SHA:?}
      KG_GIT_SHA_SHORT:   ${KG_GIT_SHA_SHORT:?}
      KG_GIT_BRANCH:      ${KG_GIT_BRANCH:?}
      KG_GIT_DIRTY:       ${KG_GIT_DIRTY:-false}

  import:
    volumes:
      - ./output-alpha:/output  # separate from dev's ./output

  post-process:
    environment:
      # Auth still off during post-process — auth is enabled only on `deploy`
      # (post-process and deploy share the volume, but auth init runs lazily)
      NEO4J_dbms_security_auth__enabled: "false"
      KG_RELEASE_VERSION: ${KG_RELEASE_VERSION}
      KG_GIT_SHA:         ${KG_GIT_SHA}
      KG_GIT_SHA_SHORT:   ${KG_GIT_SHA_SHORT}
      KG_GIT_BRANCH:      ${KG_GIT_BRANCH}
      KG_GIT_DIRTY:       ${KG_GIT_DIRTY}
      KG_MCP_MIN_VERSION: ${KG_MCP_MIN_VERSION:-0.1.0}

  deploy:
    # Bind to the Tailscale IP only — never 0.0.0.0
    # ${TAILSCALE_IP} comes from .env.alpha; operator sets it after `tailscale ip -4` on the box
    ports:
      - "${TAILSCALE_IP}:17474:7474"
      - "${TAILSCALE_IP}:17687:7687"
    environment:
      NEO4J_AUTH: ${NEO4J_AUTH:?must be set, e.g. neo4j/<strong-password>}
      NEO4J_dbms_security_auth__enabled: "true"
      # Smaller heap than dev — alpha is read-only, low-concurrency
      NEO4J_server_memory_heap_initial__size: "512M"
      NEO4J_server_memory_heap_max__size: "1G"

  # Drop the Biochatter UI; testers use their own MCP clients
  app:
    profiles: ["never"]
```

`profiles: ["never"]` is the idiom for "don't start this service under any normal `up`" — testers connect via MCP, not Biochatter.

`.env.alpha.example`:
```bash
# Versioning — set by /release-kg, do not edit manually
KG_RELEASE_VERSION=0.1.0-alpha.1
KG_GIT_SHA=
KG_GIT_SHA_SHORT=
KG_GIT_BRANCH=main
KG_GIT_DIRTY=false
KG_MCP_MIN_VERSION=0.1.0

# Operator-set, persistent
TAILSCALE_IP=100.x.y.z          # run `tailscale ip -4` on the box
NEO4J_AUTH=neo4j/<strong-password>  # the admin password; per-user readonly accounts created via Cypher (§2.5)
```

### 2.5 Per-user Neo4j readonly accounts

Run once per release (or whenever a tester is added/removed), via cypher-shell against `bolt://localhost:17687` from the box itself, authenticated as `neo4j`:

```cypher
:use system
CREATE USER alice IF NOT EXISTS
  SET PASSWORD '<generated-strong-password>'
  CHANGE NOT REQUIRED;
GRANT ROLE reader TO alice;
// Repeat for bob, carol, dave, eve.
```

`reader` is the Neo4j built-in role: read on all graphs, no writes, no schema changes. Available in Neo4j Community Edition.

**Caveat for Community Edition multi-user:** Neo4j 5.x Community supports multiple users and the built-in `reader` role, but role customization is limited to the four built-ins (admin, architect, publisher, reader). That's fine for our needs.

Distribute the per-user passwords out-of-band (Signal, encrypted email, password manager share). Never commit them.

### 2.6 Tailscale setup (operator runbook excerpt)

Lives in the alpha guide under "Operator setup". Summary:

1. On the box: `curl -fsSL https://tailscale.com/install.sh | sh && sudo tailscale up`. Confirm with `tailscale ip -4` → record the IP for `.env.alpha`.
2. In the Tailscale admin console: enable MagicDNS, rename the box machine to `kg-box`. After this, `kg-box.<tailnet>.ts.net` resolves from any device on the tailnet.
3. Invite each tester via the admin console (email). They install Tailscale on Win11 (`https://tailscale.com/download/windows`) or Linux (`https://tailscale.com/install`), sign in, get added.
4. Optional but recommended ACL: restrict testers to only port `17687` on `kg-box` (sample policy in the guide).
5. Verify from a tester laptop: `nslookup kg-box.<tailnet>.ts.net` resolves, `nc -zv kg-box.<tailnet>.ts.net 17687` succeeds.

Tailscale is free for personal use up to 100 devices and 3 users. 5 alpha testers + the operator = within free tier. If we ever exceed 3 *users*, the next tier is ~$6/user/month.

### 2.7 Alpha-user-facing guide (`docs/kg_mcp_guide.md`)

Rewrite the existing draft for the Tailscale + local-box flow. Section structure:

1. **What this is** — unchanged.
2. **Install Tailscale and join the tailnet** — new section. Win11 + Linux side-by-side. Verify connectivity step.
3. **Install the explorer MCP** — `uvx --from git+...` line, same as before.
4. **Install the research skills plugin** — new section. Clone `multiomics_research`, follow its README to register as a Claude Code plugin.
5. **Configure your MCP client** — `NEO4J_URI: bolt://kg-box.<tailnet>.ts.net:17687`, username/password from operator. `bolt://` (no `+s`) — Tailscale handles encryption end-to-end; Neo4j cert hassle skipped.
6. **First query — compatibility check** — unchanged; matches against `Schema_info.version`.
7. **Find your way around** — unchanged.
8. **Starter queries** — unchanged.
9. **Conventions and gotchas** — unchanged.
10. **Limits and support** — adjust: drop Aura timeout note; mention dev iteration cadence; how to file issues.
11. **Updating** — releases on GitHub; pull new explorer/research; rerun the compatibility check.

## 3. Test

- Unit-test the version-regex in the release skill helper.
- `pytest -m kg` additions (run against `bolt://localhost:7687`, the dev stack):
  - `test_schema_info_release_properties` — asserts `Schema_info` has `version`, `built_at`, `git_sha`, `gene_count`, `expression_edge_count` populated and the right type.
  - `test_schema_info_counts_match` — `s.gene_count == COUNT { (:Gene) }` etc., catches stale stats.
- Manual:
  - Bring up both stacks on the box. Verify `docker ps` shows two `deploy` containers, ports differ, volumes differ. Run a probe query against both: dev returns `0.0.0-dev+…`, alpha returns the tagged version.
  - From a tester laptop: tailscale connectivity probe, MCP compatibility check, run one tool per starter-query category.
  - Run `/release-kg 0.0.0-dev.0` against an alpha checkout end-to-end; verify the alpha stack comes back up with the new `Schema_info` properties.

## 4. Review

- `post-import.cypher` and `post-import.sh` must remain byte-identical in Cypher logic — run `scripts/post-import-validate.sh` before/after.
- `Schema_info` is a BioCypher-owned node; do not duplicate by also writing a `KGRelease` label.
- Alpha port bindings must use the Tailscale IP, never `0.0.0.0`. Defense in depth above auth — grep `docker-compose.alpha.yml` and `scripts/alpha_up.sh` to confirm.
- `.env.alpha` must NOT be committed — add to `.gitignore`. Only `.env.alpha.example` is committed.
- Never run `neo4j-admin database dump` against the alpha stack without stopping its `deploy` first (locks the DB files). N/A unless we add the dump-snapshot fallback from §6.

## 5. Document

- This plan: mark workstreams complete as they ship.
- `CLAUDE.md`: add a "Two-stack deploy" subsection under Docker Pipeline Stages explaining the dev vs alpha split, and a brief "Releases" pointer to the skill + `Schema_info` properties.
- `docs/kg_mcp_guide.md`: ship rewritten alongside the first release.
- Memory: short reference memory pointing at this plan + the MCP guide so future conversations find them.

## 6. Sizing, costs, and open questions

### 6.1 Box resource budget

Two Neo4j containers on the same host:

| | Dev `deploy` | Alpha `deploy` | Notes |
|---|---:|---:|---|
| Heap (max) | 2 GB | 1 GB | Smaller on alpha — read-only, low-concurrency |
| Page cache | ~1 GB | ~1 GB | Neo4j default heuristic |
| Volume disk | ~3-5 GB | ~3-5 GB | Same data, separate copies |
| Build output | ~1-2 GB in `./output` | ~1-2 GB in `./output-alpha` | Both checked into bind-mount |

Floor: ~6 GB RAM for the two Neo4js, ~10-15 GB disk. Comfortable on any modern dev box (16 GB+ RAM, ~100 GB free); flag if the box is tighter than that.

Build-time CPU: `create_knowledge_graph.py` is single-threaded-ish; running a dev rebuild while alpha is up is fine. Running the alpha rebuild via `/release-kg` while dev is mid-build will compete — operator should pause dev rebuilds during a release.

### 6.2 Cost

- Tailscale free tier (≤100 devices, ≤3 users). Our 5 testers may push us over the user cap; if so, the next tier is "Starter" at ~$6/user/month → ~$36/month for operator + 5 testers. Verify on `https://tailscale.com/pricing/` before inviting.
- Everything else $0.

If we exceed the Tailscale free-user cap and don't want to pay, fallbacks: WireGuard hand-rolled (more setup, no MagicDNS), Cloudflare Tunnel (one tunnel, public hostname behind Cloudflare Access SSO), or simply public-Internet + Let's Encrypt + Neo4j auth (uglier, more exposed). All worse trade-offs than the Tailscale Starter tier for our scale.

### 6.3 Build-each-side vs. volume-copy release path

Default (this plan): alpha rebuilds from the tagged commit. ~30-60 min per release; high reproducibility (the tag IS the release).

Future alternative: build on dev only, then `docker run --rm -v kg-dev_biocypher_neo4j_volume:/from -v kg-alpha_biocypher_neo4j_volume:/to alpine cp -a /from/. /to/` (with both `deploy` containers stopped). Faster (~30s) but adds a "did the copy succeed?" failure mode and a "remember to bump `Schema_info` after the copy" step. Not worth it until per-release rebuild becomes painful.

### 6.4 Open questions

1. **Neo4j Community per-user accounts.** Built-in `reader` role exists in Community Edition, but verify on the actual Neo4j 5.15-community image during the first alpha cut. Fallback: shared `neo4j` admin password + rely on the MCP's read-only enforcement (acceptable for trusted alpha testers; not for broader release).
2. **Win11 uv install ergonomics.** Document `winget install astral-sh.uv` and the PowerShell installer. Confirm by running through the guide on a clean Win11 install before inviting testers.
3. **Tailscale user-count.** Confirm 5 testers + operator fits the free tier *as users* (devices is not the binding limit). If not, factor the Starter tier cost into the alpha announcement.
4. **Box reliability.** Single point of failure. Document a graceful "alpha box is down" message and where testers should look (operator's Signal? GitHub issue?). Not a blocker for first cut.
5. **Aura migration trigger.** Define ahead of time: when do we move from local-box to Aura? Candidates: (a) testers beyond the tailnet, (b) box uptime issues, (c) more than ~15 concurrent users. Record the decision once we hit one of these.
6. **MCP version pinning.** `Schema_info.mcp_min_version` is the contract. Need a real version bump cadence on `multiomics_explorer` + a way for the guide to point testers at the right tag. Same open question as before; not local-box-specific.
