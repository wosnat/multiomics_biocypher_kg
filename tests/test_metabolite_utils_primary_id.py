"""Unit tests for mnxm_to_primary_id() / mnxr_to_primary_id()."""
from __future__ import annotations

import sqlite3

from multiomics_kg.utils import metabolite_utils as mu


def _make_db(tmp_path):
    """Create a minimal in-memory resolver DB with synthetic aliases."""
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id));
        CREATE INDEX idx_compound_aliases_mnxm ON compound_aliases(mnxm_id);
        CREATE TABLE reaction_aliases (source TEXT, value TEXT, mnxr_id TEXT,
            PRIMARY KEY (source, value, mnxr_id));
        CREATE INDEX idx_reaction_aliases_mnxr ON reaction_aliases(mnxr_id);

        -- MNXM41 has both KEGG and ChEBI: KEGG should win
        INSERT INTO compound_aliases VALUES ('kegg.compound', 'C00031', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('chebi', '17234', 'MNXM41');

        -- MNXM999 has only ChEBI: chebi should win
        INSERT INTO compound_aliases VALUES ('chebi', '12345', 'MNXM999');

        -- MNXM_orphan has nothing: fallback to mnx
        -- (no rows inserted)

        -- MNXR101234 has both KEGG and Rhea: KEGG should win
        INSERT INTO reaction_aliases VALUES ('kegg.reaction', 'R00200', 'MNXR101234');
        INSERT INTO reaction_aliases VALUES ('rhea', '10828', 'MNXR101234');

        -- MNXR_rhea_only has only Rhea
        INSERT INTO reaction_aliases VALUES ('rhea', '99999', 'MNXR_rhea_only');
    """)
    conn.commit()
    return conn


def test_mnxm_kegg_wins(tmp_path):
    conn = _make_db(tmp_path)
    assert mu.mnxm_to_primary_id("MNXM41", conn) == "kegg.compound:C00031"


def test_mnxm_chebi_when_no_kegg(tmp_path):
    conn = _make_db(tmp_path)
    assert mu.mnxm_to_primary_id("MNXM999", conn) == "chebi:12345"


def test_mnxm_orphan_falls_back_to_mnx(tmp_path):
    conn = _make_db(tmp_path)
    assert mu.mnxm_to_primary_id("MNXM_orphan", conn) == "mnx:MNXM_orphan"


def test_mnxr_kegg_wins(tmp_path):
    conn = _make_db(tmp_path)
    assert mu.mnxr_to_primary_id("MNXR101234", conn) == "kegg.reaction:R00200"


def test_mnxr_rhea_when_no_kegg(tmp_path):
    conn = _make_db(tmp_path)
    assert mu.mnxr_to_primary_id("MNXR_rhea_only", conn) == "rhea:99999"


def test_mnxr_orphan_falls_back_to_mnx(tmp_path):
    conn = _make_db(tmp_path)
    assert mu.mnxr_to_primary_id("MNXR_unknown", conn) == "mnx:MNXR_unknown"
