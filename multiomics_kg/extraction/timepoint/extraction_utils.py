"""Extraction file I/O, paperconfig walking, signature computation for the
timepoint/growth-phase backfill pipeline.

Mirrors the cluster-extraction pattern in multiomics_kg/extraction/cluster/.
"""
from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Iterable, Iterator

from ruamel.yaml import YAML

logger = logging.getLogger(__name__)

EXTRACTIONS_DIR = "extractions"
EXTRACTION_FILENAME = "timepoint.json"

# Fields we pass to the LLM alongside the three target fields. Staleness
# signature includes these because a change to any of them (e.g. name_col
# swapped, experiment re-referenced) could invalidate an earlier extraction.
SIGNATURE_CONTEXT_FIELDS = ("experiment", "name_col", "logfc_col")


_YAML = YAML()
_YAML.preserve_quotes = True


def _load_yaml(path: Path) -> dict:
    with open(path) as f:
        return _YAML.load(f)


def find_analyses(paperconfig_path: Path) -> list[dict]:
    """Return a flat list of every analysis dict in the paperconfig, in file
    order. Each dict is the raw entry from `statistical_analyses[]`.
    """
    data = _load_yaml(paperconfig_path)
    pub = data.get("publication", {})
    analyses: list[dict] = []
    for table in (pub.get("supplementary_materials") or {}).values():
        if not isinstance(table, dict):
            continue
        for a in (table.get("statistical_analyses") or []):
            analyses.append(dict(a))  # shallow copy detaches from yaml
    return analyses


def compute_paperconfig_signature(paperconfig_path: Path) -> str:
    """sha1 hex over sorted analysis_ids + context fields (experiment,
    name_col, logfc_col). Stable across runs, changes whenever relevant
    paperconfig state changes.
    """
    analyses = find_analyses(paperconfig_path)
    rows = sorted(
        (a.get("id", ""), *(str(a.get(k, "")) for k in SIGNATURE_CONTEXT_FIELDS))
        for a in analyses
    )
    canonical = "\n".join("|".join(row) for row in rows)
    return hashlib.sha1(canonical.encode("utf-8")).hexdigest()


def save_extraction_json(paper_dir: Path, metadata: dict, analyses: list[dict]) -> Path:
    """Write `extractions/timepoint.json` under `paper_dir`, creating the
    subdirectory if needed. Returns the path.
    """
    ext_dir = paper_dir / EXTRACTIONS_DIR
    ext_dir.mkdir(parents=True, exist_ok=True)
    path = ext_dir / EXTRACTION_FILENAME
    payload = {"metadata": metadata, "analyses": analyses}
    path.write_text(json.dumps(payload, indent=2, default=str))
    return path


def load_extraction_json(paper_dir: Path) -> dict | None:
    """Return the parsed JSON dict, or None if the file doesn't exist."""
    path = paper_dir / EXTRACTIONS_DIR / EXTRACTION_FILENAME
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as e:
        logger.error("Failed to parse %s: %s", path, e)
        raise


def iter_paperconfigs(list_files: Iterable[Path]) -> Iterator[Path]:
    """Yield paperconfig paths listed in one-per-line text files. Skips
    blank lines and comments (lines starting with `#`). Does NOT filter
    out paperconfig paths that don't exist on disk — caller can decide.
    """
    for lf in list_files:
        for line in Path(lf).read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            yield Path(line)


TARGET_FIELDS = ("timepoint", "timepoint_hours", "growth_phase")


def compute_fields_requested(analysis: dict, validate: bool = False) -> list[str]:
    """Return the list of target fields that need LLM extraction for this
    analysis.

    Default mode (`validate=False`):
      - Field absent → request.
      - `timepoint_hours` is None → request (null is a "needs filling"
        signal for numeric fields).
      - `timepoint` or `growth_phase` is empty string → request.
      - Field has a non-null, non-empty value → skip.

    Validate mode (`validate=True`):
      - Always request all three fields (LLM re-examines even populated ones).
    """
    if validate:
        return list(TARGET_FIELDS)

    requested: list[str] = []
    for field in TARGET_FIELDS:
        if field not in analysis:
            requested.append(field)
            continue
        value = analysis[field]
        if value is None:
            requested.append(field)
        elif isinstance(value, str) and value == "":
            requested.append(field)
        # non-null, non-empty → skip
    return requested
