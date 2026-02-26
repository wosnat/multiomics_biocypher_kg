"""Shared CLI helpers for download build scripts (argparse, config, genome CSV)."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Any

import yaml

from multiomics_kg.download.utils.paths import GENOMES_CSV


def add_common_args(parser: argparse.ArgumentParser, default_config: Path) -> None:
    """Add --strains, --force, and --config arguments to *parser*."""
    parser.add_argument(
        "--strains", nargs="+", metavar="STRAIN",
        help="Only process these strains (default: all)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Overwrite existing output files",
    )
    parser.add_argument(
        "--config", default=str(default_config),
        metavar="PATH",
        help=f"YAML config path (default: {default_config})",
    )


def load_config(config_path: str | Path) -> dict[str, Any]:
    """Load a YAML config file and return its contents as a dict.

    Exits with an error message if the file does not exist.
    """
    if not Path(config_path).exists():
        print(f"Error: config not found: {config_path}", file=sys.stderr)
        sys.exit(1)
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_genome_rows(strains: list[str] | None = None) -> list[dict]:
    """Load genome rows from cyanobacteria_genomes.csv.

    Skips comment lines (starting with '#') and rows without both
    ``strain_name`` and ``data_dir``.  If *strains* is given, only rows whose
    ``strain_name`` is in that list are returned.

    Exits with an error message if no rows remain after filtering.
    """
    rows: list[dict] = []
    with open(GENOMES_CSV, newline="") as f:
        reader = csv.DictReader(
            (line for line in f if not line.strip().startswith("#"))
        )
        for row in reader:
            if row.get("strain_name") and row.get("data_dir"):
                rows.append(row)

    if strains:
        rows = [r for r in rows if r["strain_name"] in strains]
        if not rows:
            print(f"Error: no strains matched {strains}", file=sys.stderr)
            sys.exit(1)

    return rows
