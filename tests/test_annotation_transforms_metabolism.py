"""Tests for the three Phase 1.1B metabolism transforms."""
from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest

from multiomics_kg.download.utils import annotation_transforms as at


@pytest.fixture(autouse=True)
def patch_caches(monkeypatch, tmp_path):
    """Set up tiny resolver DB + tiny TCDB/CAZy hierarchies, redirect module caches."""
    # Resolver DB
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, reference TEXT,
                                formula TEXT, charge INTEGER, mass REAL,
                                inchi TEXT, inchikey TEXT, smiles TEXT);
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
                                       PRIMARY KEY(source, value, mnxm_id));
        CREATE TABLE compound_names (name_normalized TEXT, mnxm_id TEXT,
                                     PRIMARY KEY(name_normalized, mnxm_id));
        CREATE TABLE reactions (mnxr_id TEXT PRIMARY KEY, mnx_equation TEXT,
                                reference TEXT, classifs TEXT,
                                is_balanced TEXT, is_transport TEXT);
        CREATE TABLE reaction_aliases (source TEXT, value TEXT, mnxr_id TEXT,
                                       PRIMARY KEY(source, value, mnxr_id));
    """)
    conn.execute("INSERT INTO reactions VALUES ('MNXR101234', '', '', '', 'B', NULL)")
    conn.execute("INSERT INTO reaction_aliases VALUES ('kegg.reaction', 'R00299', 'MNXR101234')")
    conn.commit()
    conn.close()

    from multiomics_kg.utils import metabolite_utils as mu
    monkeypatch.setattr(mu, "DEFAULT_DB_PATH", db)
    monkeypatch.setattr(at, "_RESOLVER_CONN", None)

    # TCDB hierarchy
    from multiomics_kg.utils import tcdb_utils as tu
    tcdb_path = tmp_path / "tcdb_hierarchy.json"
    tcdb_path.write_text(json.dumps({"3.A.1": {}, "3.A.1.1": {}, "3.A.1.1.1": {}}))
    monkeypatch.setattr(tu, "DEFAULT_PATH", tcdb_path)
    monkeypatch.setattr(tu, "_CACHE", None)

    # CAZy hierarchy
    from multiomics_kg.utils import cazy_utils as cu
    cazy_path = tmp_path / "cazy_hierarchy.json"
    cazy_path.write_text(json.dumps({"GH": {}, "GH13": {}, "GH13_1": {}, "GT": {}, "GT19": {}}))
    monkeypatch.setattr(cu, "DEFAULT_PATH", cazy_path)
    monkeypatch.setattr(cu, "_CACHE", None)

    yield


def test_resolve_kegg_reaction_resolves_known():
    assert at._tx_resolve_kegg_reaction_to_mnxr("R00299") == "MNXR101234"


def test_resolve_kegg_reaction_drops_unknown():
    assert at._tx_resolve_kegg_reaction_to_mnxr("R99999") is None


def test_validate_tcdb_keeps_known():
    assert at._tx_validate_tcdb("3.A.1.1.1") == "3.A.1.1.1"
    assert at._tx_validate_tcdb("3.A.1") == "3.A.1"


def test_validate_tcdb_drops_unknown():
    assert at._tx_validate_tcdb("99.X.99") is None


def test_validate_cazy_keeps_known():
    assert at._tx_validate_cazy("GH13_1") == "GH13_1"
    assert at._tx_validate_cazy("GT19") == "GT19"


def test_validate_cazy_drops_unknown():
    assert at._tx_validate_cazy("XX99") is None
