---
date: 2026-04-14
topic: Docker import.report test + pipeline changes
status: approved (pending user review of written spec)
---

# Docker Import Report Test

## Goal

Detect silent data loss and import failures in the docker-built knowledge graph by inspecting `neo4j-admin import` artifacts from a pytest test in `tests/kg_validity/`. The test must run after `docker compose up` and fail if any relationship was skipped or if the import container exited non-zero.

## Motivation

`neo4j-admin database import full ... --skip-bad-relationships=true --skip-duplicate-nodes=true` writes skipped rows to `import.report` and, by default, continues until the skipped count exceeds `--bad-tolerance` (1000). Two failure modes are currently invisible to the test suite:

1. **Soft failure** — some relationships skipped but count below tolerance; import exits 0 and the graph is built with silent data loss. This is what BRITE dangling edges would have caused had the count been lower.
2. **Hard failure** — skipped count exceeds tolerance; neo4j-admin exits 70. This is the actual BRITE regression (2026-04-14 pruning work).

Today, `import.sh` copies `import.report` to the host-mounted `/output/` directory only on success (`cp ... || true` after import completes). On hard failure the file is never exposed. Host-side `./output/import.report` ends up empty anyway because the directory is root-owned from prior runs and the copy silently fails.

## Components

### 1. `scripts/import.sh` — capture artifacts on every exit

Replace the current "cp after success" pattern with a `trap` that runs on EXIT regardless of success or failure:

```bash
#!/bin/bash
set -euo pipefail

OUTPUT_DIR=/output
REPORT_SRC=/var/lib/neo4j/import.report

mkdir -p "$OUTPUT_DIR"
chmod 777 "$OUTPUT_DIR" || true  # fix root-owned dir from prior runs

IMPORT_EXIT=1
write_status() {
  local exit_code=$?
  echo "exit_code=${IMPORT_EXIT}" > "$OUTPUT_DIR/import.status"
  cp "$REPORT_SRC" "$OUTPUT_DIR/import.report" 2>/dev/null || echo "" > "$OUTPUT_DIR/import.report"
  # Also keep a copy inside the volume for forensics
  cp "$REPORT_SRC" /data/build2neo/import.report 2>/dev/null || true
  exit "$exit_code"
}
trap write_status EXIT

sleep 2
if [ ! -f /data/build2neo/neo4j-admin-import-call.sh ]; then
  echo "ERROR: /data/build2neo/neo4j-admin-import-call.sh not found. Build step may have failed." >&2
  exit 1
fi
chmod +x /data/build2neo/neo4j-admin-import-call.sh

# Run the import without letting `set -e` abort before we write the status file.
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

Post-conditions after this runs (success or failure):
- `/output/import.status` exists and contains `exit_code=N`.
- `/output/import.report` exists (may be empty).

### 2. `docker-compose.yml` — no structural changes

The existing bind mount `./output:/output` is sufficient. The only change would be documenting that `./output/` must be writable by the neo4j UID in the import container; `import.sh` handles this via `chmod 777`. No compose edits required.

### 3. `tests/kg_validity/test_import_report.py` — new test file

```python
"""Validate neo4j-admin import artifacts from the docker build.

Skips when the artifacts are absent (local dev without a recent docker build).
Runs as part of `pytest -m kg`.
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
    # Sanity: the file must be parseable.
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

Three test functions, each auto-skipping if its artifact is missing, so running `pytest -m kg` without a recent docker build still passes the suite.

### 4. No conftest or fixture changes

Paths are derived from `__file__`. Existing `tests/kg_validity/conftest.py` (neo4j driver fixture) is untouched.

## Data flow

```
docker compose up
  └─ build       → writes CSVs to /data/build2neo/
  └─ import      → runs neo4j-admin import
                   └─ trap on EXIT:
                      • /output/import.status   ← "exit_code=N"
                      • /output/import.report   ← skipped rows (may be empty)
  └─ post-process / deploy / app

pytest -m kg
  └─ test_import_report.py reads ./output/import.{status,report}
```

## Acceptance criteria

1. After `docker compose up` on a clean build:
   - `./output/import.status` contains `exit_code=0`.
   - `./output/import.report` exists and is empty.
   - `pytest tests/kg_validity/test_import_report.py -v` passes 3/3.
2. After a deliberate regression (e.g., re-introduce BRITE dangling edges or corrupt a Gene ID):
   - `./output/import.report` is non-empty OR `./output/import.status` contains a non-zero exit code.
   - `pytest tests/kg_validity/test_import_report.py -v` fails with a message pointing at the report file.
3. Running `pytest -m kg` on a checkout with no `./output/` artifacts: all three new tests skip (no failures).

## Non-goals

- Parsing build-container `WARNING --` lines. Too noisy; will be a separate spec if needed.
- Post-process (cypher-shell) error detection. Out of scope; post-process already exits non-zero on hard failures and `service_completed_successfully` gates the downstream stages.
- Thresholded tolerance for bad rows. Strict zero-tolerance matches the snapshot regression pattern; add a threshold later only if a known-noise source appears.
- Docker-daemon introspection (`docker inspect`). The status-file approach keeps the test pure-filesystem and works in any environment where `./output/` was populated.

## Testing plan

1. Apply the `import.sh` change locally.
2. Run `docker compose up` (volumes already exist; import should complete in <10s given the cached build).
3. Verify `./output/import.status` and `./output/import.report` appear on the host.
4. Run `pytest tests/kg_validity/test_import_report.py -v` — expect 3 passes.
5. Manually inject one bad relationship CSV row, rebuild, confirm `test_import_report_empty` fails with the expected preview message.
6. Restore, rebuild, confirm the suite returns to green.
