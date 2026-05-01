"""Build the MNX metabolite resolver SQLite. Invoked by scripts/refresh_mnx.sh.

This module does only one thing: build the heavy MNX SQLite (~30 min). It reads
the 4 MNX TSVs and produces:

- cache/data/mnx/metabolite_resolver.db    (SQLite)
- cache/data/mnx/metabolite_id_mapping_report.json

(TCDB hierarchy parsing has moved to build_kegg_metabolism_xrefs.py — step 6 of
prepare_data.sh.)
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


# ── aliases ───────────────────────────────────────────────────────────────────

# Source-prefix normalization: collapse MNX's dual short/long forms.
# Long bioregistry-style spelling wins. Keys: any prefix MNX uses; values:
# canonical form. Prefixes not in the map are passed through unchanged.
_SOURCE_NORM = {
    "CHEBI":           "chebi",
    "chebi":           "chebi",
    "keggC":           "kegg.compound",
    "kegg.compound":   "kegg.compound",
    "keggD":           "kegg.drug",
    "kegg.drug":       "kegg.drug",
    "keggG":           "kegg.glycan",
    "kegg.glycan":     "kegg.glycan",
    "biggM":           "bigg.metabolite",
    "bigg.metabolite": "bigg.metabolite",
    "metacycM":        "metacyc.compound",
    "metacyc.compound":"metacyc.compound",
    "seedM":           "seed.compound",
    "seed.compound":   "seed.compound",
    "vmhM":            "vmh",
    "vmh":             "vmh",
    "lipidmapsM":      "lipidmaps",
    "lipidmaps":       "lipidmaps",
    "SLM":             "slm",
    "slm":             "slm",
    "reactomeM":       "reactome",
    "reactome":        "reactome",
    "envipathM":       "envipath",
    "envipath":        "envipath",
    "sabiorkM":        "sabiork",
    "sabiork":         "sabiork",
    # Reaction-side equivalents
    "keggR":           "kegg.reaction",
    "kegg.reaction":   "kegg.reaction",
    "biggR":           "bigg.reaction",
    "bigg.reaction":   "bigg.reaction",
    "seedR":           "seed.reaction",
    "seed.reaction":   "seed.reaction",
    "vmhreaction":     "vmhreaction",
    "vmhR":            "vmhreaction",
    "metacycR":        "metacyc.reaction",
    "metacyc.reaction":"metacyc.reaction",
    "sabiork.reaction":"sabiork.reaction",
    "rh":              "rhea",
    "rhea":            "rhea",
}


def _split_source_value(combined: str) -> tuple[str, str] | None:
    """Split MNX xref col-1 ('chebi:17234') into (canonical_source, value).

    Returns None for self-references (bare MNXM/MNXR ids without a colon).
    """
    if ":" not in combined:
        return None
    src, _, value = combined.partition(":")
    canonical = _SOURCE_NORM.get(src, src)
    return canonical, value


_ALIAS_DDL = """
CREATE TABLE IF NOT EXISTS compound_aliases (
    source  TEXT,
    value   TEXT,
    mnxm_id TEXT,
    PRIMARY KEY (source, value, mnxm_id)
);
CREATE INDEX IF NOT EXISTS idx_compound_aliases_value
    ON compound_aliases(value);
CREATE INDEX IF NOT EXISTS idx_compound_aliases_source_value
    ON compound_aliases(source, value);
CREATE INDEX IF NOT EXISTS idx_compound_aliases_mnxm
    ON compound_aliases(mnxm_id);
"""


def build_compound_aliases_table(conn: sqlite3.Connection, chem_xref_path: Path) -> int:
    """Populate compound_aliases from chem_xref.tsv with source normalization."""
    cur = conn.cursor()
    cur.executescript(_ALIAS_DDL)
    n = 0
    for fields in _iter_data_rows(chem_xref_path):
        if len(fields) < 2:
            continue
        combined, mnxm_id = fields[0], fields[1]
        parsed = _split_source_value(combined)
        if parsed is None:
            continue  # skip self-references
        source, value = parsed
        cur.execute(
            "INSERT OR IGNORE INTO compound_aliases (source, value, mnxm_id) "
            "VALUES (?, ?, ?)",
            (source, value, mnxm_id),
        )
        n += 1
    log.info(f"  compound_aliases: {n} rows scanned")
    return n


# ── names ─────────────────────────────────────────────────────────────────────

_NAMES_DDL = """
CREATE TABLE IF NOT EXISTS compound_names (
    name_normalized TEXT,
    mnxm_id         TEXT,
    PRIMARY KEY (name_normalized, mnxm_id)
);
CREATE INDEX IF NOT EXISTS idx_compound_names
    ON compound_names(name_normalized);
"""


_NAME_NORM_RE = re.compile(r"[^\w+\-/ ]")


def _normalize_name(s: str) -> str:
    """Lowercase, collapse whitespace, strip punctuation except + - /."""
    s = _NAME_NORM_RE.sub(" ", s.lower())
    return " ".join(s.split())


def build_compound_names_table(conn: sqlite3.Connection, chem_xref_path: Path) -> int:
    """Populate compound_names from the description column of chem_xref.tsv.

    Description uses '||' as alternate-name separator. Each non-empty alt-name
    contributes one (name_normalized, mnxm_id) row.
    """
    cur = conn.cursor()
    cur.executescript(_NAMES_DDL)
    n = 0
    for fields in _iter_data_rows(chem_xref_path):
        if len(fields) < 3:
            continue
        mnxm_id, description = fields[1], fields[2]
        for raw in description.split("||"):
            normalized = _normalize_name(raw)
            if not normalized:
                continue
            cur.execute(
                "INSERT OR IGNORE INTO compound_names (name_normalized, mnxm_id) "
                "VALUES (?, ?)",
                (normalized, mnxm_id),
            )
            n += 1
    log.info(f"  compound_names: {n} rows scanned")
    return n


# ── reactions ─────────────────────────────────────────────────────────────────

_REACTIONS_DDL = """
CREATE TABLE IF NOT EXISTS reactions (
    mnxr_id      TEXT PRIMARY KEY,
    mnx_equation TEXT,
    reference    TEXT,
    classifs     TEXT,
    is_balanced  TEXT,
    is_transport TEXT
);
"""


def build_reactions_table(conn: sqlite3.Connection, reac_prop_path: Path) -> int:
    """Populate the `reactions` table from reac_prop.tsv (6 columns)."""
    cur = conn.cursor()
    cur.executescript(_REACTIONS_DDL)
    n = 0
    for fields in _iter_data_rows(reac_prop_path):
        fields = (fields + [""] * 6)[:6]
        mnxr_id, mnx_equation, reference, classifs, is_balanced, is_transport = fields
        cur.execute(
            "INSERT OR REPLACE INTO reactions "
            "(mnxr_id, mnx_equation, reference, classifs, is_balanced, is_transport) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (mnxr_id, mnx_equation, reference, classifs,
             is_balanced or None, is_transport or None),
        )
        n += 1
    log.info(f"  reactions: {n} rows")
    return n


_REAC_ALIAS_DDL = """
CREATE TABLE IF NOT EXISTS reaction_aliases (
    source  TEXT,
    value   TEXT,
    mnxr_id TEXT,
    PRIMARY KEY (source, value, mnxr_id)
);
CREATE INDEX IF NOT EXISTS idx_reaction_aliases_value
    ON reaction_aliases(value);
CREATE INDEX IF NOT EXISTS idx_reaction_aliases_source_value
    ON reaction_aliases(source, value);
CREATE INDEX IF NOT EXISTS idx_reaction_aliases_mnxr
    ON reaction_aliases(mnxr_id);
"""


def build_reaction_aliases_table(conn: sqlite3.Connection, reac_xref_path: Path) -> int:
    """Populate reaction_aliases from reac_xref.tsv with source normalization."""
    cur = conn.cursor()
    cur.executescript(_REAC_ALIAS_DDL)
    n = 0
    for fields in _iter_data_rows(reac_xref_path):
        if len(fields) < 2:
            continue
        combined, mnxr_id = fields[0], fields[1]
        parsed = _split_source_value(combined)
        if parsed is None:
            continue
        source, value = parsed
        cur.execute(
            "INSERT OR IGNORE INTO reaction_aliases (source, value, mnxr_id) "
            "VALUES (?, ?, ?)",
            (source, value, mnxr_id),
        )
        n += 1
    log.info(f"  reaction_aliases: {n} rows scanned")
    return n


# ── orchestration ─────────────────────────────────────────────────────────────

CACHE_ROOT = Path("cache/data")


def _resolve_paths(cache_root: Path) -> dict[str, Path]:
    return {
        "chem_prop":     cache_root / "mnx" / "chem_prop.tsv",
        "chem_xref":     cache_root / "mnx" / "chem_xref.tsv",
        "reac_prop":     cache_root / "mnx" / "reac_prop.tsv",
        "reac_xref":     cache_root / "mnx" / "reac_xref.tsv",
        "resolver_db":   cache_root / "mnx" / "metabolite_resolver.db",
        "report":        cache_root / "mnx" / "metabolite_id_mapping_report.json",
    }


def _read_mnx_release(chem_prop_path: Path) -> str:
    """Extract '#VERSION' and '#DATE' from chem_prop.tsv comment header."""
    version = ""
    date = ""
    with open(chem_prop_path, encoding="utf-8") as f:
        for line in f:
            if not line.startswith("#"):
                break
            if "VERSION" in line:
                version = line.split(":", 1)[1].strip().rstrip("\t")
            elif "DATE" in line:
                date = line.split(":", 1)[1].strip().rstrip("\t")
            if version and date:
                break
    if version and date:
        # Normalize MNX's `YYYY/MM/DD` to `YYYY-MM-DD` for stable downstream IDs
        date = date.replace("/", "-")
        return f"MNXref {version} ({date})"
    return "MNXref unknown"


def _alias_count_by_source(conn: sqlite3.Connection) -> dict[str, int]:
    cur = conn.cursor()
    cur.execute("SELECT source, COUNT(*) FROM compound_aliases GROUP BY source ORDER BY 2 DESC")
    return dict(cur.fetchall())


def main(force: bool = False) -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    paths = _resolve_paths(CACHE_ROOT)

    if not force and paths["resolver_db"].exists():
        log.info("Resolver DB present; skip (use --force to rebuild).")
        return

    log.info("Building MNX metabolite resolver SQLite")

    # SQLite resolver
    paths["resolver_db"].parent.mkdir(parents=True, exist_ok=True)
    if paths["resolver_db"].exists():
        paths["resolver_db"].unlink()
    conn = sqlite3.connect(paths["resolver_db"])
    n_compounds = build_compounds_table(conn, paths["chem_prop"])
    build_compound_aliases_table(conn, paths["chem_xref"])
    build_compound_names_table(conn, paths["chem_xref"])
    n_reactions = build_reactions_table(conn, paths["reac_prop"])
    build_reaction_aliases_table(conn, paths["reac_xref"])
    conn.commit()

    alias_counts = _alias_count_by_source(conn)
    conn.close()

    # Diagnostic report
    report = {
        "mnx_release":              _read_mnx_release(paths["chem_prop"]),
        "compound_count":           n_compounds,
        "reaction_count":           n_reactions,
        "alias_counts_by_source":   alias_counts,
    }
    paths["report"].write_text(json.dumps(report, indent=2))
    log.info(f"  metabolite_id_mapping_report.json written ({len(report)} fields)")
    log.info("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Rebuild caches even if they exist.")
    args = parser.parse_args()
    main(force=args.force)
