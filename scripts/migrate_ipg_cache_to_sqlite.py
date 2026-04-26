#!/usr/bin/env python3
"""One-time migration: cache/data/ncbi_ipg/*.tsv  →  cache/data/ncbi_ipg.sqlite.

Per-accession TSV files were the original cache format for NCBI IPG lookups.
At ~5000+ accessions they impose ~4× filesystem block-size overhead, slow
git status, and bloat directory listings. SQLite stores the same data in
one file with indexed lookups and atomic writes.

Idempotent: re-running with both the directory and the DB present is a
no-op (INSERT OR REPLACE just rewrites the same body). After a successful
migration, ``git rm -r cache/data/ncbi_ipg/`` removes the old files.
"""
from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

OLD_DIR = Path("cache/data/ncbi_ipg")
NEW_DB = Path("cache/data/ncbi_ipg.sqlite")


def main() -> None:
    if not OLD_DIR.exists():
        print(f"No old cache directory at {OLD_DIR}; nothing to migrate.")
        return

    NEW_DB.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(NEW_DB)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS ipg_cache (
            accession  TEXT PRIMARY KEY,
            body       TEXT NOT NULL,
            fetched_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)

    tsv_files = sorted(OLD_DIR.glob("*.tsv"))
    if not tsv_files:
        print(f"No *.tsv files in {OLD_DIR}; nothing to migrate.")
        conn.close()
        return

    n = 0
    with conn:  # atomic transaction
        for tsv in tsv_files:
            accession = tsv.stem  # filename minus .tsv == "NP_892211.1"
            body = tsv.read_text()
            conn.execute(
                "INSERT OR REPLACE INTO ipg_cache (accession, body) VALUES (?, ?)",
                (accession, body),
            )
            n += 1
    conn.close()

    db_size_mb = NEW_DB.stat().st_size / 1024 / 1024
    print(f"Migrated {n} cache entries → {NEW_DB} ({db_size_mb:.1f} MB)")
    print(f"Old directory ({OLD_DIR}) can now be removed with `git rm -r`.")


if __name__ == "__main__":
    sys.exit(main())
