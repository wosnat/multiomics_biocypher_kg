"""
Tests that validate all paperconfig.yaml files listed in paperconfig_files.txt.

Reads the file list from data/Prochlorococcus/papers_and_supp/paperconfig_files.txt
and runs the validate_paperconfig.validate() function on each one, ensuring every
config passes validation.
"""

import os
import sys
import pytest
from pathlib import Path

# Add the validate script to the import path
VALIDATE_SCRIPT_DIR = os.path.join(
    os.path.dirname(__file__), os.pardir, ".claude", "skills", "paperconfig"
)
sys.path.insert(0, os.path.abspath(VALIDATE_SCRIPT_DIR))

from validate_paperconfig import validate

# Project root (one level up from tests/)
PROJECT_ROOT = Path(__file__).resolve().parent.parent

PAPERCONFIG_LIST = PROJECT_ROOT / "data" / "Prochlorococcus" / "papers_and_supp" / "paperconfig_files.txt"


def _load_paperconfig_paths() -> list[str]:
    """Read paperconfig file paths from the listing file."""
    assert PAPERCONFIG_LIST.exists(), f"Paperconfig list file not found: {PAPERCONFIG_LIST}"
    paths = []
    with open(PAPERCONFIG_LIST) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                paths.append(line)
    assert len(paths) > 0, "No paperconfig paths found in listing file"
    return paths


PAPERCONFIG_PATHS = _load_paperconfig_paths()


@pytest.mark.parametrize(
    "config_path",
    PAPERCONFIG_PATHS,
    ids=[Path(p).parent.name for p in PAPERCONFIG_PATHS],
)
def test_paperconfig_validates(config_path: str, monkeypatch):
    """Each paperconfig.yaml listed in paperconfig_files.txt must pass validation."""
    # Run validation from the project root so relative paths in configs resolve
    monkeypatch.chdir(PROJECT_ROOT)

    assert os.path.exists(config_path), f"Config file not found: {config_path}"
    result = validate(config_path)
    assert result is True, f"Validation failed for {config_path}"
