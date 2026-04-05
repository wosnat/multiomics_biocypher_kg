#!/usr/bin/env python3
"""One-time: migrate extraction JSONs to new format + location.

Reads: data/<paper>/cluster_extraction_{entry}.json (old format with stage2_results)
Writes: data/<paper>/cluster_extractions/{entry}.json (new format with clusters)
Deletes: old files after successful migration
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from multiomics_kg.extraction.cluster.extraction_utils import save_extraction


def migrate():
    old_files = sorted(Path("data").rglob("cluster_extraction_*.json"))
    old_files = [f for f in old_files if "pilot" not in f.name and "cluster_extractions" not in str(f)]

    print(f"Found {len(old_files)} files to migrate\n")

    for old_path in old_files:
        paper_dir = old_path.parent
        entry_key = old_path.stem.replace("cluster_extraction_", "")

        data = json.loads(old_path.read_text())

        if "clusters" in data:
            clusters = data["clusters"]
        elif "stage2_results" in data:
            clusters = data["stage2_results"]
        else:
            print(f"  SKIP {old_path.name} — no clusters or stage2_results key")
            continue

        metadata = data.get("metadata", {})

        new_path = save_extraction(paper_dir, entry_key, metadata, clusters)
        print(f"  {old_path.name} -> {new_path.relative_to(paper_dir)}")

        # Delete old files
        old_path.unlink()
        old_md = old_path.with_suffix(".md")
        if old_md.exists():
            old_md.unlink()
            print(f"    Deleted {old_md.name}")

    print("\nDone.")


if __name__ == "__main__":
    migrate()
