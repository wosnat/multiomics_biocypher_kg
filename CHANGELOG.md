# Changelog

All notable changes to the multi-omics knowledge graph are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versions use the KG release scheme `X.Y.Z[-(alpha|beta|rc).N]` and are tagged
`kg-X.Y.Z…`.

**Process (accumulate-then-cut):** log notable changes under **[Unreleased]** as
they land. At release time, `/release-kg` *cuts* `[Unreleased]` into a dated
version section, stamps the same version onto `Schema_info.version`, and renders
the GitHub Release notes from that section. The changelog is the source of truth;
the GitHub Release is a rendering of one section. See `plans/alpha_release.md` §2.3.

**Authoring conventions for preflight subsections.** Each version section MAY
include two preflight-facing subsections that are extracted at build time and
stamped onto `Schema_info.release_highlights` / `Schema_info.release_breaking`.
MCP/explorer clients surface them when a user connects, so they answer two
preflight questions:

- `### Highlights` — "what can I now ask that I couldn't before?" One bullet
  per user-facing capability or data layer added in this release. Cap at ~5
  bullets. Plain prose, no internal jargon.
- `### Breaking` — "did anything change meaning under me?" One bullet per
  silently-changed semantics (redefinitions, repurposed fields, default flips).
  Errors/removals are also fine here. **Omit the subsection entirely when there
  are no breaking items** — do not leave an empty `### Breaking` heading
  (renders as a blank bullet on the client).

Both are extracted verbatim (markdown). The rest of the version section
(`### Added` / `### Changed` / `### Fixed`) is unchanged in role.

## [Unreleased]

### Highlights
- **Preflight release summary in your MCP/explorer client.** Each KG release
  can now ship two short markdown bullets — what's new and what changed
  meaning — that your client surfaces the moment it connects. Lets you
  answer "what can I now ask that I couldn't before?" and "did anything
  silently change under me?" before you start querying, instead of
  discovering a redefined field through a wrong answer.
- **Two more Prochlorococcus strains: MIT1327 (LLIV) and MIT1314 (HLII).**
  The last two strains of the Soussan 2025 N/P-starvation panel are now in the
  graph, so cross-strain queries and ortholog comparisons span the full panel.
  MIT1327 includes curated CyanoRak ortholog clusters and roles.

### Breaking

### Added
- **Prochlorococcus MIT1327 + MIT1314** — the two Soussan 2025 N/P-starvation
  panel strains previously deferred for lack of a public assembly, completing
  the 15-strain panel (47 OrganismTaxon nodes / 40 genome strains total).
  Genome-only (no DE paper). MIT1327 (`GCF_001632125.1`, taxid 1801626, clade
  LLIV) carries full CyanoRak annotation (CK ortholog clusters + CyanoRak/TIGR
  roles) on 2308 genes, bridged from a CyanoRak-team gbff via diamond
  protein-sequence matching because the strain has no public CyanoRak server
  export — see the new reusable converter
  `multiomics_kg/download/convert_cyanorak_gbff_to_gff.py`. MIT1314
  (`GCF_034093315.1`, taxid 3096220, clade HLII) is NCBI/eggNOG-only. Both
  carry eggNOG, PSORTb, SignalP, and tcdb-diamond annotation.
- New optional string properties on the `Schema_info` node:
  `release_highlights` and `release_breaking`. Sourced from `### Highlights`
  and `### Breaking` subsections inside each version's CHANGELOG entry,
  extracted at build time by `extract_preflight_subsection()` in
  `release_kg.py`, and stamped by `post-import.sh` (Group 4). Absent or
  empty subsections → real `null` property (not empty-string), so legacy
  KGs and no-subsection releases are indistinguishable on the wire.
  Existing `kg_release_info` callers pick both up automatically (they
  already return `Schema_info { .* }`).
- CHANGELOG preamble: authoring conventions for the two preflight
  subsections (~5 bullets each, omit `### Breaking` entirely when nothing
  to flag rather than leaving an empty heading).

### Changed

### Fixed
- `/release-kg --target local`: two first-cut deploy bugs surfaced while cutting
  kg-0.1.0-alpha.5 to the lab box.
  - **Live flip re-stamped `Schema_info` to `0.0.0-dev`.** `_alpha_flip_live_deploy`
    brought the live `alpha-deploy` up with `docker compose -p kg-alpha up -d deploy`;
    `deploy`'s `depends_on` chain re-ran build→import→post-process in the live
    project, and because that env intentionally omits `KG_RELEASE_VERSION`, the
    re-import re-stamped the release graph with the dev fallback version. Fixed by
    adding `--no-deps` so the live deploy just *mounts* the already-built/verified
    `kg-alpha-<color>` volume. (Also added to the rollback path.)
  - **Explorer provisioning aborted the deploy on a fresh volume.**
    `_alpha_provision_explorer_user` ran `CREATE USER explorer IF NOT EXISTS …` then
    `ALTER USER explorer SET PASSWORD …` to the same value; on a fresh blue/green
    volume the `ALTER`-to-identical-password is rejected by Neo4j ("Old and new
    password cannot be the same"), crashing the deploy after the graph was already
    live. Replaced with a single idempotent `CREATE OR REPLACE USER explorer …`.
- `/release-kg` SKILL.md: un-stubbed the `--target local` documentation (it was
  still described as raising `NotImplementedError`); documented the two flip/provision
  gotchas above.
- `scripts/alpha_firewall.sh`: new operator helper to apply the `DOCKER-USER`
  allowlist restricting the alpha ports (`17474`/`17687`) to a confirmed lab subnet
  (plan §2.6). Parameterized on the CIDR, idempotent (`-C` check before `-I`), and
  refuses to run against the campus-wide `/16`.

## [0.1.0-alpha.5] - 2026-06-09

### Added
- **Publication "discusses" edges** — a narrative literature index linking each
  publication to the genes and KEGG pathways it discusses in prose (regulators,
  model genes, pathways named in text), distinct from the supplementary DE-table
  expression data. Recall-biased *router* ("which papers discuss gene/pathway X?"),
  best-effort, not exhaustive. New relationship types:
  - `Publication_discusses_gene` (Publication → Gene)
  - `Publication_discusses_kegg_pathway` (Publication → KeggTerm, pathway-level)

  Both carry `prominence` (`central` | `peripheral`) and the extraction `evidence`
  quote. ~1,099 gene + ~140 pathway edges across 40 publications.

  Three-stage pipeline (`plans/publication_discusses_edges.md`):
  - **Extract** — `multiomics_kg/extraction/topics/` + `/extract-discussed-topics`
    skill: full-PDF LLM extraction with 15-page chunking and reference-page skipping
    (fits per-request token caps; large PDFs no longer fail). Writes
    `publication_topics/topics.json`. Strain-aware (attributes each gene to one of the
    paper's strains); captures verbatim locus tags; self-reports an
    `uncaptured_identifiers` triage signal.
  - **Resolve** — `prepare_data.sh` **step 8**
    (`multiomics_kg/download/resolve_paper_topics.py`): genes resolved
    per-strain via `gene_id_mapping` (identifiers-first; gene families fan out to one
    edge per member; `all`/`unspecified` mentions resolve in each paper strain);
    pathways via a global `kegg_data.json` lookup (dangling-proof — only resolves to
    KeggTerm nodes that exist). Writes `topics_resolved.json` + a `resolution_report.txt`
    (per-paper stats, method breakdown, truncated-id count, `unresolved_reasons` tally).
  - **Adapter** — `multiomics_kg/adapters/publication_topics_adapter.py`: pure edge
    adapter; source = the paper's Publication node id (DOI from paperconfig or
    PDF-extraction cache); targets `ncbigene:{locus_tag}` / `kegg.pathway:ko*`. Wired
    into `create_knowledge_graph.py` after the omics adapter.
- Post-import rollups: `Publication.discussed_gene_count` / `discussed_pathway_count`,
  `Gene.discussed_in_publication_count` (in both `post-import.cypher` and `.sh`).
- Two `config/schema_config.yaml` edge types for the above.
- Tests: 44 unit tests (resolution, chunking/merge, adapter) + 8 KG-validity tests
  (`tests/kg_validity/test_discuss_edges.py`: endpoint correctness/no-dangling,
  prominence enum, evidence presence, post-import rollups).
- Shared `multiomics_kg/extraction/pdf_utils.py` helpers: `upload_pdf` (lifted from
  cluster extraction), `count_pages`, `find_references_page`, `page_chunks`,
  `write_page_range_pdf`.

### Changed

### Fixed

## [0.1.0-alpha.4] - 2026-06-08

### Added
- Track A infrastructure scaffolding (no behavior yet — these wire
  into `/release-kg --target local` in a follow-up):
  - `docker-compose.alpha.yml` — alpha-stack compose override
    matching the plan §2.4 sketch. Renames containers `alpha-build` /
    `alpha-import` / `alpha-post-process` / `alpha-deploy`, forwards
    `KG_*` env vars to post-process, drops the Biochatter UI, and
    hands the data volume off to a `${ALPHA_DATA_VOLUME}`-selected
    external name (default `kg-alpha-blue`) for the blue/green flip.
  - `.env.alpha.example` — committed template for the operator-only
    secrets (`ALPHA_BIND_IP`, `NEO4J_AUTH`, `ALPHA_EXPLORER_PASSWORD`).
  - `scripts/alpha_up.sh` / `scripts/alpha_down.sh` — operator
    wrappers around the `-p kg-alpha -f … -f …` invocation. Read the
    active color from `.alpha_active_color` (written by the release
    flow) and refuse to run if `.env.alpha` or the marker is missing.
  - `.gitignore` adds `.env.alpha` and `.alpha_active_color`.
- `/release-kg --target local` now implements the Track A lab-box
  deploy (was a stub). Replaces `deploy_local_stub` with the real
  blue/green flip orchestration in `release_kg.py`:
  - Reads `.env.alpha` for `ALPHA_BIND_IP`, `NEO4J_AUTH`,
    `ALPHA_EXPLORER_PASSWORD`; refuses to run on unfilled
    `REPLACE_WITH_…` / `<...>` placeholders.
  - Determines active/inactive colors from `.alpha_active_color`
    (first cut bootstraps into `blue`).
  - Builds the alpha-inactive color via a transient
    `kg-alpha-build` compose project on temp localhost ports;
    verifies Schema_info via `docker exec`. Pytest L2 is skipped
    here because Phase 5 already ran it against the same tagged
    KG on the staging stack (same tag + env = same KG).
  - Flips: stops the live alpha-deploy (if any), brings it up on
    the new color bound to `ALPHA_BIND_IP`. Best-effort rollback
    to the previous color on flip failure.
  - Provisions the shared `explorer` read login idempotently via
    `CREATE USER explorer IF NOT EXISTS … + ALTER … CHANGE NOT
    REQUIRED`.
  - Warns (does not block) if the `DOCKER-USER` chain has no rule
    for `:17474`/`:17687` — the firewall allowlist needs sudo on
    the lab box.
  - Updates `.alpha_active_color`, prints an operator summary
    (Bolt URI, browser URI, distribute-credentials reminder).
  - 18 new unit tests cover the deterministic pieces — color
    rotation, marker round-trip, env-alpha parsing, validation
    (missing keys, unfilled placeholders, malformed NEO4J_AUTH).

### Changed
- Decision (2026-06-06): the alpha runs on **Track A** — the lab box
  at `132.75.249.47` (leadership choice). Aura (was Track B) archived;
  `/release-kg` drops the `--target aura` backend. Plan
  (`plans/alpha_release.md`) and tester guide
  (`docs/kg_mcp_guide.md` §2) reframed accordingly. Remaining
  agnostic-vs-local design split stays intact; `--target local`
  is now implemented (see Added above) and this release is the
  first cut deployed to the lab box via it.

### Fixed

## [0.1.0-alpha.3] - 2026-06-06

### Added
- `/release-kg` Phase 5 now gates the release on the KG validity suite
  passing against the staging stack (`uv run pytest tests/kg_validity/
  --neo4j-url bolt://localhost:27687 -q`, ~73 s, 1012 assertions).
  Catches structural / semantic regressions before Phase 7 publishes a
  GitHub Release. `--skip-kg-tests` flag bypasses for emergencies. On
  failure, the staging stack is left running so the operator can
  inspect.
- `/release-kg` Phase 7 now embeds a "What changed since kg-X.Y.Z"
  diff block in the GitHub Release notes. Compares the prior published
  release's `metadata.json` to the current build's metadata: headline
  Schema_info count deltas (papers / experiments / genes / organisms /
  expression edges) + per-publication expression-edge changes (new /
  changed / removed). The per-publication detail catches the
  net-zero-but-paper-A-lost-paper-B-gained regression class that
  Schema_info totals would hide. `metadata.json` now also carries the
  full `per_publication_edges: {doi: count}` map for the *next*
  release's diff. Soft-fails on first-ever release / older prior
  releases without `metadata.json` — release still publishes.

### Changed
- `/release-kg`'s default `--mcp-min` value now reads from
  `[tool.release-kg].mcp_min_version` in `pyproject.toml` (hard
  fallback `0.1.0` if the file/section is missing). Previously the
  default was a Python constant inside the skill script. Repo-wide
  bumps of the cross-repo explorer-MCP contract now live in
  declarative config; per-release override via `--mcp-min` is
  unchanged.
- Bumped `[tool.release-kg].mcp_min_version` from `0.1.0` to
  `0.1.0a1` (PEP 440 alpha). Matches the current
  `multiomics_explorer` version, so the explorer's
  `kg_release_info()` compat check passes against builds off this
  branch (string equality on `mcp_min_version`).

### Fixed

## [0.1.0-alpha.2] - 2026-06-02

### Added

### Changed

### Fixed
- `docker-compose.staging.yml` now forwards `KG_RELEASE_VERSION` and
  `KG_GIT_*` env vars from the compose process into the `post-process`
  container, so post-import.sh Group 4 stamps `Schema_info.version` with
  the tagged version instead of silently falling back to `0.0.0-dev`.
  The 0.1.0-alpha.1 cut built and deployed successfully but stamped the
  wrong version because of this gap, which is why the tag exists but no
  GitHub Release was published for it; this release is the first cut
  with a correctly-stamped staging verify.

## [0.1.0-alpha.1] - 2026-06-02

### Added
- `Schema_info` release metadata, stamped at post-import (every build): `version`,
  `built_at`, `git_sha`, `git_sha_short`, `git_branch`, `git_dirty`,
  `mcp_min_version`, `release_notes_url`, plus computed counts (`paper_count`,
  `experiment_count`, `gene_count`, `organism_count`, `expression_edge_count`).
  Dev builds stamp `0.0.0-dev`; releases stamp the tagged version. Added in both
  `scripts/post-import.sh` (Group 4) and `scripts/post-import.cypher`.

### Changed

### Fixed
