"""Step 0 sub-step 7 — Build resolver + hierarchy caches.

Reads the seven sub-step-6 cache files (4 MNX TSVs + 3 TCDB TSVs) plus raw
eggNOG annotation files from configured strains, and produces:

- cache/data/mnx/metabolite_resolver.db    (SQLite)
- cache/data/tcdb/tcdb_hierarchy.json
- cache/data/cazy/cazy_hierarchy.json
- cache/data/mnx/metabolite_id_mapping_report.json
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import sqlite3
from pathlib import Path

log = logging.getLogger(__name__)


# ── compounds ─────────────────────────────────────────────────────────────────

_COMPOUNDS_DDL = """
CREATE TABLE IF NOT EXISTS compounds (
    mnxm_id   TEXT PRIMARY KEY,
    name      TEXT,
    reference TEXT,
    formula   TEXT,
    charge    INTEGER,
    mass      REAL,
    inchi     TEXT,
    inchikey  TEXT,
    smiles    TEXT
);
"""


def _iter_data_rows(path: Path):
    """Yield non-comment, non-blank lines from a tab-delimited MNX file."""
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            yield line.rstrip("\n").split("\t")


def _coerce_int(s: str) -> int | None:
    return int(s) if s and s.lstrip("-").isdigit() else None


def _coerce_float(s: str) -> float | None:
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def build_compounds_table(conn: sqlite3.Connection, chem_prop_path: Path) -> int:
    """Populate the `compounds` table from chem_prop.tsv. Returns row count."""
    cur = conn.cursor()
    cur.executescript(_COMPOUNDS_DDL)
    n = 0
    for fields in _iter_data_rows(chem_prop_path):
        # 9 columns expected; pad with empties if shorter
        fields = (fields + [""] * 9)[:9]
        mnxm_id, name, reference, formula, charge, mass, inchi, inchikey, smiles = fields
        cur.execute(
            "INSERT OR REPLACE INTO compounds "
            "(mnxm_id, name, reference, formula, charge, mass, inchi, inchikey, smiles) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (mnxm_id, name, reference, formula,
             _coerce_int(charge), _coerce_float(mass),
             inchi, inchikey, smiles),
        )
        n += 1
    log.info(f"  compounds: {n} rows")
    return n


def main(force: bool = False) -> None:
    raise NotImplementedError("Phase 1.1B — see plan for follow-up tasks")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Rebuild caches even if they exist.")
    args = parser.parse_args()
    main(force=args.force)
