"""
Smoke test for the main knowledge graph build script.

Runs create_knowledge_graph.py as a subprocess and fails if:
- The process exits with a non-zero return code
- Any line in stdout/stderr starts with "ERROR"
"""

import subprocess
import sys
import os
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent


@pytest.mark.slow
def test_create_knowledge_graph_no_errors():
    """Run the full KG build and assert no ERROR lines appear."""
    result = subprocess.run(
        [sys.executable, "create_knowledge_graph.py"],
        cwd=str(PROJECT_ROOT),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=3600,  # 1 hour max
    )

    combined_output = (result.stdout or "") + "\n" + (result.stderr or "")

    # Collect all lines that start with ERROR
    error_lines = [
        line for line in combined_output.splitlines()
        if line.strip().startswith("ERROR")
    ]

    # Fail with details if the process crashed
    assert result.returncode == 0, (
        f"Script exited with code {result.returncode}.\n"
        f"stderr:\n{(result.stderr or '')[-2000:]}"
    )

    # Fail if any ERROR lines were logged
    assert len(error_lines) == 0, (
        f"Found {len(error_lines)} ERROR line(s) in output:\n"
        + "\n".join(error_lines)
    )
