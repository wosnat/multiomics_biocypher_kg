"""Extract OrthologGroup descriptions from eggnog.db into a lightweight JSON cache.

Reads all ortholog_groups from gene_annotations_merged.json across strains,
queries the local eggnog.db SQLite database for descriptions, and writes
a cache file at cache/data/eggnog/og_descriptions.json.

This avoids needing the 39GB eggnog.db at KG build time (e.g., in Docker).

Usage:
    uv run python -m multiomics_kg.download.build_og_descriptions [--force]
"""

from __future__ import annotations

import argparse
import json
import os
import sqlite3
import sys
from pathlib import Path

from multiomics_kg.download.utils.annotation_transforms import is_eggnog_description_stub
from multiomics_kg.download.utils.cli import load_genome_rows
from multiomics_kg.utils.gene_id_utils import load_gene_annotations

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent

CACHE_FILE = PROJECT_ROOT / "cache" / "data" / "eggnog" / "og_descriptions.json"


def _find_eggnog_db() -> Path | None:
    """Find eggnog.db from EGGNOG_DATA_DIR env var or .env file."""
    eggnog_dir = os.environ.get("EGGNOG_DATA_DIR", "")
    if not eggnog_dir:
        env_path = PROJECT_ROOT / ".env"
        if env_path.exists():
            for line in env_path.read_text().splitlines():
                line = line.strip()
                if line.startswith("EGGNOG_DATA_DIR="):
                    eggnog_dir = line.split("=", 1)[1].strip().strip('"').strip("'")
                    break
    if not eggnog_dir:
        return None
    db_path = Path(eggnog_dir).expanduser() / "eggnog.db"
    return db_path if db_path.exists() else None


def _collect_eggnog_og_names() -> set[str]:
    """Collect all unique eggNOG OG names from gene_annotations_merged.json."""
    names: set[str] = set()

    for row in load_genome_rows():
        genes = load_gene_annotations(row["data_dir"])
        if genes is None:
            continue
        for gene in genes.values():
            for og in gene.get("ortholog_groups", []):
                og_id = og.get("og_id", "")
                if not og_id.startswith("eggnog:"):
                    continue
                raw = og_id[len("eggnog:"):]
                if "@" in raw:
                    names.add(raw.rsplit("@", 1)[0])

    return names


def _query_descriptions(db_path: Path, og_names: set[str]) -> dict[str, str]:
    """Query eggnog.db for descriptions. Returns {og_name@level: description}."""
    if not og_names:
        return {}

    descriptions: dict[str, str] = {}
    batch_size = 500
    og_list = list(og_names)

    with sqlite3.connect(db_path) as conn:
        for i in range(0, len(og_list), batch_size):
            batch = og_list[i:i + batch_size]
            placeholders = ",".join("?" * len(batch))
            cursor = conn.execute(
                f"SELECT og, level, description FROM og WHERE og IN ({placeholders})",
                batch,
            )
            for og, level, desc in cursor:
                if desc and not is_eggnog_description_stub(desc):
                    descriptions[f"{og}@{level}"] = desc

    return descriptions


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract eggNOG OG descriptions to cache JSON")
    parser.add_argument("--force", action="store_true", help="Overwrite existing cache file")
    args = parser.parse_args()

    if CACHE_FILE.exists() and not args.force:
        print(f"Cache file already exists: {CACHE_FILE}")
        print("Use --force to overwrite.")
        return

    db_path = _find_eggnog_db()
    if not db_path:
        print("ERROR: eggnog.db not found. Set EGGNOG_DATA_DIR or add it to .env", file=sys.stderr)
        sys.exit(1)

    print(f"eggnog.db: {db_path}")

    print("Collecting OG names from gene_annotations_merged.json...")
    og_names = _collect_eggnog_og_names()
    print(f"  Found {len(og_names)} unique eggNOG OG names")

    print("Querying eggnog.db...")
    descriptions = _query_descriptions(db_path, og_names)
    print(f"  Got {len(descriptions)} descriptions")

    # Write cache
    CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(CACHE_FILE, "w") as f:
        json.dump(descriptions, f, indent=1, sort_keys=True)
    print(f"Wrote {CACHE_FILE} ({CACHE_FILE.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
