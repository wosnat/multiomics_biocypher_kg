"""Shared I/O for per-strain tool calls.json + skill_summary.json files.

Single source for path resolution, JSON read/write, and cross-strain
iteration. Used by:

- Each tool runner (`.claude/skills/<tool>-run/run_<tool>.py`) — to write
  per-strain outputs.
- Cross-strain QC narratives + future spot-check verifiers — to walk all
  strains and read their calls.
- Phase-2 KG adapters that merge tool output into
  `gene_annotations_merged.json` or the BioCypher graph — to walk all
  strains' calls.json without each adapter reimplementing path logic.

Schema-aware iteration ("yield every per-protein record" vs "yield every
per-region record") stays in the per-tool `multiomics_kg/utils/<tool>.py`
module because it depends on Step 0 Q1 (the prediction shape). This module
only handles paths + JSON load/save + cross-strain walks; it returns raw
`Any` for calls (dict for per-protein tools, list for per-region tools).

ONE-TIME CREATION. The first new tool that follows `/add-a-tool` copies
this template into `multiomics_kg/utils/tool_calls_io.py` and commits it.
Subsequent tools just `from multiomics_kg.utils.tool_calls_io import ...`.
Later, opportunistically migrate the four existing runners (psortb-run,
tcdb-diamond, eggnog-run, signalp-run) off their inline I/O to use this
module — not blocking work.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterator

from multiomics_kg.download.utils.cli import load_genome_rows


# ----------------------------------------------------------------------------
# Path resolution
# ----------------------------------------------------------------------------

def _tool_dir(data_dir: Path, tool: str) -> Path:
    """`<data_dir>/<tool>/` — created on demand by save_*."""
    return Path(data_dir) / tool


def calls_path(
    data_dir: Path,
    tool: str,
    strain: str,
    *,
    limited: int | None = None,
) -> Path:
    """Path to `<data_dir>/<tool>/<strain>.<tool>.calls.json`.

    When `limited` is set (smoke test via `--limit N`), the filename gets a
    `.limited_<N>.` infix so the artifact is auto-gitignored via the
    project's `cache/data/*/genomes/*/<tool>/*.limited_*` rule.
    """
    infix = f".limited_{limited}" if limited is not None else ""
    return _tool_dir(data_dir, tool) / f"{strain}.{tool}{infix}.calls.json"


def skill_summary_path(
    data_dir: Path,
    tool: str,
    strain: str,
    *,
    limited: int | None = None,
) -> Path:
    """Path to `<data_dir>/<tool>/<strain>.<tool>.skill_summary.json`."""
    infix = f".limited_{limited}" if limited is not None else ""
    return _tool_dir(data_dir, tool) / f"{strain}.{tool}{infix}.skill_summary.json"


# ----------------------------------------------------------------------------
# Read / write
# ----------------------------------------------------------------------------

def load_calls(data_dir: Path, tool: str, strain: str) -> Any:
    """Read `<strain>.<tool>.calls.json`.

    Return shape depends on the tool (Step 0 Q1):
    - per-protein tools: `dict[str, dict]` keyed by WP_ accession
    - per-region tools: `list[dict]` of `{contig, start, end, ...}`
    - per-genome scalar tools: `dict` with `value` + metadata

    Use the per-tool `multiomics_kg/utils/<tool>.py` module for
    schema-aware iteration over the returned object.
    """
    with calls_path(data_dir, tool, strain).open() as f:
        return json.load(f)


def save_calls(
    data_dir: Path,
    tool: str,
    strain: str,
    calls: Any,
    *,
    limited: int | None = None,
) -> None:
    """Write calls.json. Indented + sorted-keys for clean PR diffs."""
    path = calls_path(data_dir, tool, strain, limited=limited)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(calls, f, indent=2, sort_keys=True)


def load_skill_summary(
    data_dir: Path,
    tool: str,
    strain: str,
    *,
    limited: int | None = None,
) -> dict:
    """Read `<strain>.<tool>.skill_summary.json`."""
    with skill_summary_path(data_dir, tool, strain, limited=limited).open() as f:
        return json.load(f)


def save_skill_summary(
    data_dir: Path,
    tool: str,
    strain: str,
    summary: dict,
    *,
    limited: int | None = None,
) -> None:
    """Write skill_summary.json. Indented for human inspection."""
    path = skill_summary_path(data_dir, tool, strain, limited=limited)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(summary, f, indent=2)


# ----------------------------------------------------------------------------
# Cross-strain iteration (for QC narratives + Phase-2 adapters)
# ----------------------------------------------------------------------------

def iter_strain_calls(
    tool: str,
    strains: list[str] | None = None,
) -> Iterator[tuple[str, Path, Any]]:
    """Yield `(strain_name, data_dir, calls)` for every strain in
    `cyanobacteria_genomes.csv` that has a calls.json for `tool`.

    Strains without output for this tool are silently skipped — use
    `iter_strain_calls_status` if you need to distinguish "ran but empty"
    from "didn't run yet".
    """
    for row in load_genome_rows(strains=strains):
        strain = row["strain_name"]
        data_dir = Path(row["data_dir"])
        if calls_path(data_dir, tool, strain).exists():
            yield strain, data_dir, load_calls(data_dir, tool, strain)


def iter_strain_calls_status(
    tool: str,
    strains: list[str] | None = None,
) -> Iterator[tuple[str, Path, Any | None]]:
    """Like `iter_strain_calls` but yields `(strain, data_dir, None)` for
    strains where the calls.json is missing. Useful for batch coverage
    reports."""
    for row in load_genome_rows(strains=strains):
        strain = row["strain_name"]
        data_dir = Path(row["data_dir"])
        path = calls_path(data_dir, tool, strain)
        calls = load_calls(data_dir, tool, strain) if path.exists() else None
        yield strain, data_dir, calls
