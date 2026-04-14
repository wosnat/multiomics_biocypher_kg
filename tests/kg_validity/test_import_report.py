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
