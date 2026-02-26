"""Shared path constants and organism-group inference for download build scripts."""

from __future__ import annotations

from pathlib import Path

# utils/ → download/ → multiomics_kg/ → project root
SCRIPT_DIR = Path(__file__).parent.parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
GENOMES_CSV = PROJECT_ROOT / "data/Prochlorococcus/genomes/cyanobacteria_genomes.csv"


def infer_organism_group(path: str | Path) -> str:
    """Infer organism group from a data or cache path.

    'cache/data/Prochlorococcus/genomes/MED4/' → 'Prochlorococcus'
    """
    parts = Path(path).parts
    for i, part in enumerate(parts):
        if part == "data" and i + 1 < len(parts):
            return parts[i + 1]
    return "Prochlorococcus"
