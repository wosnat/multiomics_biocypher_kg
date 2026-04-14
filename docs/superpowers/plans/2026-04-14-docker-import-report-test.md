# Docker Import Report Test Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a pytest test in `tests/kg_validity/` that fails when `neo4j-admin import` reports skipped rows or exits non-zero, plus the `scripts/import.sh` changes needed to expose those artifacts to the host.

**Architecture:** `scripts/import.sh` writes two host-visible files via a `trap` on EXIT — `/output/import.status` (contains `exit_code=N`) and `/output/import.report` (neo4j-admin's skipped-rows log). A new file `tests/kg_validity/test_import_report.py` reads both via relative paths, auto-skipping when absent.

**Tech Stack:** Bash (trap), pytest, docker compose (Neo4j 5.15 community).

**Reference spec:** `docs/superpowers/specs/2026-04-14-docker-import-report-test-design.md`

---

## File Structure

- **Create:** `tests/kg_validity/test_import_report.py` — three `@pytest.mark.kg` tests: status file present, exit_code == 0, import.report empty.
- **Modify:** `scripts/import.sh` — replace the linear "run import then copy report" flow with a trap-on-EXIT pattern that always emits both artifacts.
- **No changes:** `docker-compose.yml`, `tests/kg_validity/conftest.py`, any other test files.

---

### Task 1: Write the failing tests first

TDD order: tests are written while the artifacts don't yet exist. They will skip cleanly (no failure) until `import.sh` is updated and docker re-runs. This verifies the skip path works.

**Files:**
- Create: `tests/kg_validity/test_import_report.py`

- [ ] **Step 1: Create the test file**

Write this exact content to `tests/kg_validity/test_import_report.py`:

```python
"""Validate neo4j-admin import artifacts from the docker build.

Skips when the artifacts are absent (local dev without a recent
`docker compose up`). Runs as part of `pytest -m kg`.

See docs/superpowers/specs/2026-04-14-docker-import-report-test-design.md
"""
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = REPO_ROOT / "output"
REPORT_PATH = OUTPUT_DIR / "import.report"
STATUS_PATH = OUTPUT_DIR / "import.status"


pytestmark = pytest.mark.kg


def _parse_status(path: Path) -> int:
    for line in path.read_text().splitlines():
        key, _, value = line.partition("=")
        if key.strip() == "exit_code":
            return int(value.strip())
    raise ValueError(f"no exit_code in {path}")


def test_import_status_file_present():
    if not STATUS_PATH.exists():
        pytest.skip(f"{STATUS_PATH} missing — run `docker compose up` first")
    _parse_status(STATUS_PATH)


def test_import_exit_code_zero():
    if not STATUS_PATH.exists():
        pytest.skip(f"{STATUS_PATH} missing — run `docker compose up` first")
    exit_code = _parse_status(STATUS_PATH)
    assert exit_code == 0, (
        f"neo4j-admin import exited {exit_code}; "
        f"inspect {REPORT_PATH} and `docker logs import` for details"
    )


def test_import_report_empty():
    if not REPORT_PATH.exists():
        pytest.skip(f"{REPORT_PATH} missing — run `docker compose up` first")
    contents = REPORT_PATH.read_text()
    if contents.strip():
        lines = contents.splitlines()
        preview = "\n".join(lines[:10])
        pytest.fail(
            f"{REPORT_PATH} has {len(lines)} skipped-row entries "
            f"(dangling relationships or duplicate nodes). First 10:\n{preview}"
        )
```

- [ ] **Step 2: Verify the tests skip cleanly while artifacts are absent**

First delete any stale artifacts from prior runs (the `./output/` dir currently contains stale files owned by root):

```bash
sudo rm -f ./output/import.report ./output/import.status
```

Then run just the new test file:

```bash
uv run pytest tests/kg_validity/test_import_report.py -v
```

Expected output: 3 SKIPPED, 0 failed, 0 passed. Each skip message should mention `run \`docker compose up\` first`.

- [ ] **Step 3: Commit the test file**

```bash
git add tests/kg_validity/test_import_report.py
git commit -m "$(cat <<'EOF'
test(kg): add docker import.report validity tests

Three @pytest.mark.kg tests that read ./output/import.status and
./output/import.report and fail if neo4j-admin skipped any
relationships or exited non-zero. Auto-skip when artifacts are
absent so the kg suite still runs cleanly without a fresh docker
build.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Update import.sh to always emit artifacts via trap

**Files:**
- Modify: `scripts/import.sh` (complete rewrite — 16 lines → ~35 lines)

- [ ] **Step 1: Replace the script contents**

Overwrite `scripts/import.sh` with this exact content:

```bash
#!/bin/bash
set -euo pipefail

OUTPUT_DIR=/output
REPORT_SRC=/var/lib/neo4j/import.report
VOLUME_REPORT=/data/build2neo/import.report

mkdir -p "$OUTPUT_DIR"
# The host-side ./output/ mount may be root-owned from prior runs; make it
# writable so cp below can always succeed.
chmod 777 "$OUTPUT_DIR" || true

IMPORT_EXIT=1
write_status() {
  local trap_exit=$?
  echo "exit_code=${IMPORT_EXIT}" > "$OUTPUT_DIR/import.status" || true
  if [ -f "$REPORT_SRC" ]; then
    cp "$REPORT_SRC" "$OUTPUT_DIR/import.report" || true
    cp "$REPORT_SRC" "$VOLUME_REPORT" || true
  else
    : > "$OUTPUT_DIR/import.report" || true
  fi
  exit "$trap_exit"
}
trap write_status EXIT

sleep 2
if [ ! -f /data/build2neo/neo4j-admin-import-call.sh ]; then
  echo "ERROR: /data/build2neo/neo4j-admin-import-call.sh not found. Build step may have failed." >&2
  exit 1
fi
chmod +x /data/build2neo/neo4j-admin-import-call.sh

# Run the import without letting `set -e` abort before we capture the exit code.
set +e
/data/build2neo/neo4j-admin-import-call.sh
IMPORT_EXIT=$?
set -e

if [ "$IMPORT_EXIT" -ne 0 ]; then
  echo "neo4j-admin import failed with exit code $IMPORT_EXIT" >&2
  exit "$IMPORT_EXIT"
fi

neo4j start
sleep 10
neo4j stop
```

Key behavior changes vs. the old script:
- `IMPORT_EXIT` defaults to 1 (so any exit before the import runs — e.g., missing CSV files — is recorded as a failure).
- `IMPORT_EXIT` is set to neo4j-admin's actual exit code after the import call.
- The EXIT trap runs on success, explicit `exit`, SIGTERM, etc. — it always writes `import.status` and copies `import.report`.
- Both `/output/import.report` and the in-volume `/data/build2neo/import.report` are populated (belt and suspenders for forensics).

- [ ] **Step 2: Shell-lint the script**

```bash
bash -n scripts/import.sh
```

Expected: no output (syntax OK). If you have shellcheck available:

```bash
shellcheck scripts/import.sh
```

Expected: no errors. A `SC2317` info about unreachable code after `exit` in the trap is acceptable; ignore.

- [ ] **Step 3: Commit the script change**

```bash
git add scripts/import.sh
git commit -m "$(cat <<'EOF'
fix(docker): always emit import.report + exit-code on import container exit

Replace the post-success cp with a trap-on-EXIT that writes
/output/import.status (exit_code=N) and copies
/var/lib/neo4j/import.report to both /output and the volume.
Guarantees the host sees the skipped-row report even when
neo4j-admin hits --bad-tolerance and exits 70. Also chmod 777 on
/output to recover from root-owned bind mounts left by prior runs.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Re-run the docker import stage and verify artifacts appear

**Files:** none modified.

- [ ] **Step 1: Remove the stopped import container so compose re-creates it**

```bash
docker compose rm -f import
```

Expected output: `Going to remove import` then `Removed import`. If the container is already gone, `docker compose rm` exits 0 with no action.

- [ ] **Step 2: Run just the import stage**

```bash
docker compose up import
```

The build stage is already completed successfully and its volume is cached, so this re-runs `neo4j-admin import` against the existing `/data/build2neo/` CSVs (~5s import + ~10s neo4j start/stop). Watch for `IMPORT DONE in Xs` and a clean `Stopping Neo4j` at the end. Exit status of `docker compose up import` should be 0.

- [ ] **Step 3: Verify host artifacts**

```bash
ls -la output/import.status output/import.report
cat output/import.status
wc -c output/import.report
```

Expected:
- Both files exist and are owned by `root` or `neo4j` (not missing).
- `output/import.status` contains exactly `exit_code=0`.
- `output/import.report` is 0 bytes.

- [ ] **Step 4: Run the new tests against the real artifacts**

```bash
uv run pytest tests/kg_validity/test_import_report.py -v
```

Expected: 3 PASSED, 0 failed, 0 skipped.

- [ ] **Step 5: Run the full kg suite to check nothing else regressed**

```bash
uv run pytest -m kg -v
```

Expected: all previously-passing tests still pass; the 3 new ones pass; any pre-existing failures (e.g. the documented orphan-protein failures) remain unchanged. No new failures caused by this change.

No commit this task — only verification.

---

### Task 4: Negative control — inject a bad row and confirm the test fails

This is a one-off sanity check that the test actually catches a regression. Do not commit the injected data; it gets reverted at the end.

**Files:** temporary, no commits.

- [ ] **Step 1: Inspect a small relationship CSV to find its format**

```bash
docker run --rm -v biocypher_neo4j_volume:/data neo4j:5.15-community head -3 /data/build2neo/Has_experiment-header.csv /data/build2neo/Has_experiment-part000.csv
```

Note the two ID columns (`:START_ID` / `:END_ID`) so you can construct a valid-looking bad row.

- [ ] **Step 2: Append a row with a non-existent target ID**

```bash
docker run --rm -v biocypher_neo4j_volume:/data neo4j:5.15-community \
  bash -c 'printf "10.9999/fake-doi_bad_experiment\tncbigene:THIS_GENE_DOES_NOT_EXIST\tHas_experiment\n" >> /data/build2neo/Has_experiment-part000.csv'
```

Adjust column order to match the header printed in Step 1 if different.

- [ ] **Step 3: Re-run the import**

```bash
docker compose rm -f import
docker compose up import
```

The import should still exit 0 (one bad row is below `--bad-tolerance=1000`), but `import.report` should now have one line.

- [ ] **Step 4: Confirm the test now fails**

```bash
uv run pytest tests/kg_validity/test_import_report.py::test_import_report_empty -v
```

Expected: FAIL with a message like `output/import.report has 1 skipped-row entries ...` and a preview line referencing the fake gene ID.

- [ ] **Step 5: Revert the injected row and rebuild**

```bash
docker compose rm -f build import
docker compose up build
docker compose up import
```

This re-runs the build from scratch, overwriting the tampered CSV, then re-imports.

- [ ] **Step 6: Confirm green state**

```bash
uv run pytest tests/kg_validity/test_import_report.py -v
```

Expected: 3 PASSED.

No commit this task — negative-control only.

---

### Task 5: Update CLAUDE.md KG Validity section

**Files:**
- Modify: `CLAUDE.md` — append one row to the KG validity test table.

- [ ] **Step 1: Edit the test-files table**

Find the markdown table under `### Test files` in the "KG Validity Tests" section of `CLAUDE.md` and append this row before the closing lines:

```markdown
| `test_import_report.py` | Docker import artifact validation: reads `./output/import.status` + `./output/import.report`, fails if `neo4j-admin import` exited non-zero or skipped any relationships. Auto-skips when artifacts are missing. |
```

- [ ] **Step 2: Commit the docs update**

```bash
git add CLAUDE.md
git commit -m "$(cat <<'EOF'
docs: note test_import_report.py in KG validity test table

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Acceptance criteria (from the spec)

1. `./output/import.status` contains `exit_code=0` after a clean `docker compose up`.
2. `./output/import.report` exists and is 0 bytes.
3. `pytest tests/kg_validity/test_import_report.py -v` → 3 passed.
4. After injecting one bad row (Task 4): `test_import_report_empty` fails with a useful preview.
5. On a fresh checkout with no `./output/` artifacts: all three tests skip (no failures).

## Self-review notes

- Spec coverage: components 1 (import.sh), 2 (compose — no-op confirmed), 3 (test file), 4 (no conftest changes) all mapped to tasks.
- TDD order: tests written before the script change (Task 1 before Task 2). Tests skip cleanly while artifacts are absent, then turn green once Task 3 regenerates them.
- Type/path consistency: `OUTPUT_DIR` = `/output` in `import.sh` ↔ `REPO_ROOT / "output"` in the test; `import.status` and `import.report` filenames identical in both places.
- Negative control is included as a dedicated task (Task 4) so the engineer actually verifies the test catches the regression it was built to catch.
