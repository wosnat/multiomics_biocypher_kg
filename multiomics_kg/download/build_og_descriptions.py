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
import csv
import json
import os
import sqlite3
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent

CACHE_FILE = PROJECT_ROOT / "cache" / "data" / "eggnog" / "og_descriptions.json"
GENOMES_CSV = PROJECT_ROOT / "data" / "Prochlorococcus" / "genomes" / "cyanobacteria_genomes.csv"


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


def _collect_eggnog_og_keys(genomes_csv: Path) -> set[tuple[str, str]]:
    """Collect all (og_name, level) pairs from gene_annotations_merged.json."""
    keys: set[tuple[str, str]] = set()

    with open(genomes_csv, newline="", encoding="utf-8") as fh:
        lines = [line for line in fh if not line.lstrip().startswith("#")]
    reader = csv.DictReader(lines)

    for row in reader:
        data_dir = row.get("data_dir", "").strip()
        if not data_dir:
            continue
        merged_path = Path(data_dir) / "gene_annotations_merged.json"
        if not merged_path.exists():
            continue
        with open(merged_path) as f:
            genes = json.load(f)
        for gene in genes.values():
            for og in gene.get("ortholog_groups", []):
                og_id = og.get("og_id", "")
                if not og_id.startswith("eggnog:"):
                    continue
                raw = og_id[len("eggnog:"):]
                if "@" not in raw:
                    continue
                name, level = raw.rsplit("@", 1)
                keys.add((name, level))

    return keys


def _query_descriptions(db_path: Path, og_keys: set[tuple[str, str]]) -> dict[str, str]:
    """Query eggnog.db for descriptions. Returns {og_name@level: description}."""
    if not og_keys:
        return {}

    conn = sqlite3.connect(str(db_path))
    descriptions: dict[str, str] = {}

    unique_names = list({name for name, _ in og_keys})
    batch_size = 500

    for i in range(0, len(unique_names), batch_size):
        batch = unique_names[i:i + batch_size]
        placeholders = ",".join("?" * len(batch))
        cursor = conn.execute(
            f"SELECT og, level, description FROM og WHERE og IN ({placeholders})",
            batch,
        )
        for og, level, desc in cursor:
            if desc:
                key = f"{og}@{level}"
                descriptions[key] = desc

    conn.close()
    return descriptions


def main():
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
    print(f"Genomes CSV: {GENOMES_CSV}")

    # Collect all eggNOG OG keys from gene annotations
    print("Collecting OG keys from gene_annotations_merged.json...")
    og_keys = _collect_eggnog_og_keys(GENOMES_CSV)
    print(f"  Found {len(og_keys)} unique eggNOG OG (name, level) pairs")

    # Query descriptions
    print("Querying eggnog.db...")
    descriptions = _query_descriptions(db_path, og_keys)
    print(f"  Got {len(descriptions)} descriptions")

    # Write cache
    CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(CACHE_FILE, "w") as f:
        json.dump(descriptions, f, indent=1)
    print(f"Wrote {CACHE_FILE} ({CACHE_FILE.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
