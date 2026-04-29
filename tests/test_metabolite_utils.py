"""Tests for utils/metabolite_utils.py — resolver accessors."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from multiomics_kg.utils import metabolite_utils as mu


@pytest.fixture
def tiny_resolver(tmp_path) -> Path:
    """Build a tiny resolver DB inline: 2 compounds, 4 aliases, 3 names; 2 reactions, 3 aliases."""
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    cur.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, reference TEXT,
                                formula TEXT, charge INTEGER, mass REAL,
                                inchi TEXT, inchikey TEXT, smiles TEXT);
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
                                       PRIMARY KEY (source, value, mnxm_id));
        CREATE TABLE compound_names (name_normalized TEXT, mnxm_id TEXT,
                                     PRIMARY KEY (name_normalized, mnxm_id));
        CREATE TABLE reactions (mnxr_id TEXT PRIMARY KEY, mnx_equation TEXT,
                                reference TEXT, classifs TEXT,
                                is_balanced TEXT, is_transport TEXT);
        CREATE TABLE reaction_aliases (source TEXT, value TEXT, mnxr_id TEXT,
                                       PRIMARY KEY (source, value, mnxr_id));
    """)
    cur.execute("INSERT INTO compounds VALUES ('MNXM41', 'D-glucose', '', '', 0, 0, '', '', '')")
    cur.execute("INSERT INTO compounds VALUES ('MNXM999', 'ambig-name', '', '', 0, 0, '', '', '')")
    cur.executemany("INSERT INTO compound_aliases VALUES (?, ?, ?)", [
        ("chebi", "17234", "MNXM41"),
        ("kegg.compound", "C00031", "MNXM41"),
        ("ambiguous-x", "X", "MNXM41"),
        ("ambiguous-x", "X", "MNXM999"),
    ])
    cur.executemany("INSERT INTO compound_names VALUES (?, ?)", [
        ("d-glucose", "MNXM41"),
        ("dextrose", "MNXM41"),
        ("ambig-name", "MNXM999"),
    ])
    cur.execute("INSERT INTO reactions VALUES ('MNXR101234', '1 MNXM3 = 1 MNXM41', '', '', 'B', NULL)")
    cur.executemany("INSERT INTO reaction_aliases VALUES (?, ?, ?)", [
        ("kegg.reaction", "R00299", "MNXR101234"),
        ("rhea", "16332", "MNXR101234"),
        ("rhea", "DUP", "MNXR101234"),
    ])
    conn.commit()
    conn.close()
    return db


def test_open_resolver_reads_db(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    cur = conn.cursor()
    cur.execute("SELECT COUNT(*) FROM compounds")
    assert cur.fetchone()[0] == 2
    conn.close()


def test_resolve_metabolite_direct_mnxm(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    assert mu.resolve_metabolite("MNXM41", conn) == ("MNXM41", "xref:exact")


def test_resolve_metabolite_alias(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    assert mu.resolve_metabolite("17234", conn) == ("MNXM41", "xref:exact")
    assert mu.resolve_metabolite("C00031", conn) == ("MNXM41", "xref:exact")


def test_resolve_metabolite_alias_ambiguous(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    mnxm, method = mu.resolve_metabolite("X", conn)
    assert mnxm in ("MNXM41", "MNXM999")
    assert method == "xref:ambiguous"


def test_resolve_metabolite_name(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    assert mu.resolve_metabolite("D-Glucose", conn) == ("MNXM41", "name:normalized")
    # Whitespace and case variation also resolves
    assert mu.resolve_metabolite("  Dextrose  ", conn) == ("MNXM41", "name:normalized")


def test_resolve_metabolite_unresolved(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    assert mu.resolve_metabolite("does-not-exist", conn) == (None, "unresolved")


def test_resolve_reaction_direct(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    assert mu.resolve_reaction("MNXR101234", conn) == ("MNXR101234", "xref:exact")


def test_resolve_reaction_alias(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    assert mu.resolve_reaction("R00299", conn) == ("MNXR101234", "xref:exact")
    assert mu.resolve_reaction("16332", conn) == ("MNXR101234", "xref:exact")


def test_resolve_reaction_unresolved(tiny_resolver):
    conn = mu.open_resolver(tiny_resolver)
    assert mu.resolve_reaction("R99999", conn) == (None, "unresolved")
