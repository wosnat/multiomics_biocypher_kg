# InterProScan `/interproscan-run` skill — implementation plan

Tracks `/add-a-tool interproscan`. Spec:
`docs/superpowers/specs/2026-07-22-interproscan-domains-design.md`.

## Step 0 — intent
- [x] 4-bullet intent block (see spec)
- [x] Design spec written
- [x] This plan written

## Step 1 — install script + SKILL.md skeleton
- [x] `.claude/skills/interproscan-run/SKILL.md` from template, frontmatter +
      Quick Start + One-Time Setup filled for the flag surface
- [x] Runner supports `--prepare-image` (docker pull + digest, exit) and
      `--refresh-data` (download + md5-verify + extract data tarball, exit)

## Step 2 — actually install + fill One-Time Setup
- [x] `--prepare-image` run on this host (docker pull 5.78-109.0).
      RepoDigest `sha256:dc58b7c147fbbf00c2dd4f5ced42121fc1e8841fcbc7cc2c484380248ff76d11`
- [x] `--refresh-data` run (6.4 G tarball → `~/tools/InterProScan/`; md5-verified + extracted)
- [x] Extracted footprint **34 G**; 17 member-DB dirs (incl. panther/gene3d/superfamily). Recorded in SKILL.md
- [x] `INTERPROSCAN_DATA_DIR` documented in One-Time Setup

## Step 3 — runner + parser + IO module + tests
- [x] `multiomics_kg/utils/tool_calls_io.py` (copied from template)
- [x] `tests/test_tool_calls_io.py`
- [x] `multiomics_kg/utils/interproscan.py` (pure JSON parser + summarizer)
- [x] `tests/test_interproscan.py`
- [x] `.claude/skills/interproscan-run/run_interproscan.py`
- [x] SKILL.md What It Does / Output Schema / QC field docs / spot-checks table
- [ ] `.gitignore` rule `cache/data/*/genomes/*/interproscan/*.limited_*`

## Step 4 — smoke test + validate  ✅
- [x] `run_interproscan.py --strains MED4 --limit 100` — OK, 100 calls / 1000
      matches / 5 no-match / 226 s; all 18 apps ran; `--user`+`:ro` mount clean
      (host-owned outputs); parse_failures=0; sentinel_rate=0.05
- [x] MED4 spot checks pass: rbcL WP_002805854.1 → PF00016 + IPR000685;
      atpH WP_002805169.1 → PF00137 + IPR000454
- [x] Smoke artifacts confirmed gitignored (limited_* / raw.json / temp/)
- [x] Workflow When Invoked section finalized

## Step 5 — full batch (DEFERRED — user kicks off)
- [ ] `nohup … run_interproscan.py > logs/interproscan/batch.log 2>&1 &`
- [ ] Observed batch results table + cross-strain QC narrative + threshold calibration
- [ ] Commit per-strain calls.json + skill_summary.json

## Step 6 — register on new strains
- [x] `add-a-strain/SKILL.md` step 3 gets a `nohup` invocation
- [x] `CLAUDE.md` skills list note
