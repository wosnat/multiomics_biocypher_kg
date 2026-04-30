"""Resolver accessor API.

Reads the SQLite resolver DB built by sub-step 7. Used by step 2 transforms,
the Spec 1.2 scaffold builder, and Phase 2 paper-measurement extraction.
"""
from __future__ import annotations

import re
import sqlite3
from pathlib import Path

DEFAULT_DB_PATH = Path("cache/data/mnx/metabolite_resolver.db")

_NAME_NORM_RE = re.compile(r"[^\w+\-/ ]")


def _normalize_name(s: str) -> str:
    """Lowercase, collapse whitespace, strip punctuation except + - /."""
    s = _NAME_NORM_RE.sub(" ", s.lower())
    return " ".join(s.split())


def open_resolver(path: Path | None = None) -> sqlite3.Connection:
    """Open a read-only connection to the resolver DB."""
    p = Path(path) if path else DEFAULT_DB_PATH
    return sqlite3.connect(f"file:{p}?mode=ro", uri=True)


def resolve_metabolite(value: str, conn: sqlite3.Connection) -> tuple[str | None, str]:
    """Resolve a value to an MNXM ID. See spec for method enum."""
    value = value.strip()
    cur = conn.cursor()

    # 1. Direct MNXM
    if value.startswith("MNXM"):
        cur.execute("SELECT 1 FROM compounds WHERE mnxm_id = ?", (value,))
        if cur.fetchone():
            return value, "xref:exact"

    # 2. Alias match (across all sources)
    cur.execute("SELECT mnxm_id FROM compound_aliases WHERE value = ? ORDER BY mnxm_id", (value,))
    rows = cur.fetchall()
    if len(rows) == 1:
        return rows[0][0], "xref:exact"
    if len(rows) > 1:
        return rows[0][0], "xref:ambiguous"

    # 3. Name match (after normalization)
    normalized = _normalize_name(value)
    if normalized:
        cur.execute("SELECT mnxm_id FROM compound_names WHERE name_normalized = ? ORDER BY mnxm_id", (normalized,))
        rows = cur.fetchall()
        if len(rows) == 1:
            return rows[0][0], "name:normalized"
        if len(rows) > 1:
            return rows[0][0], "ambiguous"

    return None, "unresolved"


def resolve_reaction(value: str, conn: sqlite3.Connection) -> tuple[str | None, str]:
    """Resolve a value to an MNXR ID. Alias-only (no name normalization)."""
    value = value.strip()
    cur = conn.cursor()

    if value.startswith("MNXR"):
        cur.execute("SELECT 1 FROM reactions WHERE mnxr_id = ?", (value,))
        if cur.fetchone():
            return value, "xref:exact"

    cur.execute("SELECT mnxr_id FROM reaction_aliases WHERE value = ? ORDER BY mnxr_id", (value,))
    rows = cur.fetchall()
    if len(rows) == 1:
        return rows[0][0], "xref:exact"
    if len(rows) > 1:
        return rows[0][0], "xref:ambiguous"

    return None, "unresolved"


# ── Primary-ID resolution (KG-side canonical IDs) ────────────────────────────
#
# Per Spec 1.2: KEGG-native primary IDs win, with bioregistry-style fallbacks.
# Compound: kegg.compound > chebi > mnx
# Reaction: kegg.reaction > rhea > mnx

_COMPOUND_PRIMARY_PRIORITY = ("kegg.compound", "chebi")
_REACTION_PRIMARY_PRIORITY = ("kegg.reaction", "rhea")


def mnxm_to_primary_id(mnxm_id: str, conn: sqlite3.Connection) -> str:
    """Map an MNXM* ID to its canonical KG primary ID.

    Priority: kegg.compound > chebi > fallback `mnx:<MNXM>`.
    """
    cur = conn.cursor()
    for source in _COMPOUND_PRIMARY_PRIORITY:
        cur.execute(
            "SELECT value FROM compound_aliases WHERE mnxm_id = ? AND source = ? ORDER BY value LIMIT 1",
            (mnxm_id, source),
        )
        row = cur.fetchone()
        if row:
            return f"{source}:{row[0]}"
    return f"mnx:{mnxm_id}"


def mnxr_to_primary_id(mnxr_id: str, conn: sqlite3.Connection) -> str:
    """Map an MNXR* ID to its canonical KG primary ID.

    Priority: kegg.reaction > rhea > fallback `mnx:<MNXR>`.
    """
    cur = conn.cursor()
    for source in _REACTION_PRIMARY_PRIORITY:
        cur.execute(
            "SELECT value FROM reaction_aliases WHERE mnxr_id = ? AND source = ? ORDER BY value LIMIT 1",
            (mnxr_id, source),
        )
        row = cur.fetchone()
        if row:
            return f"{source}:{row[0]}"
    return f"mnx:{mnxr_id}"
