"""Tests for Phase 2 metabolomics extension of step 6 (build_kegg_metabolism_xrefs)."""
import sqlite3
from pathlib import Path


def _make_stub_resolver(tmp_path: Path) -> sqlite3.Connection:
    """Build a tiny SQLite resolver that mimics the MNX schema.

    Just enough for resolve_metabolite + mnxm_to_primary_id to work for our cases.
    """
    db_path = tmp_path / "stub.db"
    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()
    cur.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, formula TEXT, mass REAL, inchikey TEXT);
        CREATE TABLE compound_aliases (mnxm_id TEXT, source TEXT, value TEXT);
        CREATE TABLE compound_names (mnxm_id TEXT, name_normalized TEXT);
        INSERT INTO compounds VALUES ('MNXM41', 'D-Glucose', 'C6H12O6', 180.06, 'WQZGKKKJIJFFOK-GASJEMHNSA-N');
        INSERT INTO compounds VALUES ('MNXM58', 'Pyruvate', 'C3H3O3', 87.03, 'LCTONWCANYUPML-UHFFFAOYSA-M');
        INSERT INTO compound_aliases VALUES ('MNXM41', 'kegg.compound', 'C00031');
        INSERT INTO compound_aliases VALUES ('MNXM58', 'kegg.compound', 'C00022');
        INSERT INTO compound_names VALUES ('MNXM41', 'glucose');
        INSERT INTO compound_names VALUES ('MNXM58', 'pyruvate');
    """)
    conn.commit()
    return conn


def test_harvest_paper_candidates_collects_names_ids_overrides_lock(tmp_path):
    from multiomics_kg.download.build_kegg_metabolism_xrefs import _harvest_paper_candidates
    from multiomics_kg.utils.metabolite_identity.paper import Override, PaperRow, PaperTable
    from tests.metabolite_identity_helpers import make_chebi, term
    conn = _make_stub_resolver(tmp_path)
    chebi = make_chebi(CHEBI_17234=term("glucose", kegg=["C00031"]))
    t = PaperTable("p", "t", None, 0, [PaperRow("p", "t", 0, "glucose", None, None),
                                       PaperRow("p", "t", 1, "Tyrosine", "chebi:17234", "chebi_direct")])
    ov = {"p": {"pyruvic": Override("kegg.compound:C00022", "as_reported", "", "Pyruvic")}}
    lock = {"names": {"x": {"id": "kegg.compound:C00041"}}, "direct_ids": {}}
    rekeys = {}
    pc = _harvest_paper_candidates([t], ov, lock, conn, chebi, rekeys)
    assert pc.name_candidates["glucose"] == {"kegg.compound:C00031"}
    assert {"kegg.compound:C00031", "kegg.compound:C00022", "kegg.compound:C00041"} <= pc.ids
    assert rekeys == {"chebi:17234": "kegg.compound:C00031"}   # paper-path ChEBI candidate re-keyed
