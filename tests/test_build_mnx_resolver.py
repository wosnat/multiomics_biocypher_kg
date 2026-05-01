"""Unit tests for the MNX resolver builder (invoked by scripts/refresh_mnx.sh)."""
from __future__ import annotations

import json
import sqlite3
import textwrap
from pathlib import Path

from multiomics_kg.download import build_mnx_resolver as bmr


CHEM_PROP_FIXTURE = textwrap.dedent("""\
    ### MetaNetX/MNXref reconciliation ###
    #VERSION:   4.5
    #DATE:      2025/08/13
    #ID\tname\treference\tformula\tcharge\tmass\tInChI\tInChIKey\tSMILES
    MNXM01\tPMF\tmnx:PMF\tH\t1\t1.00794\tInChI=1S/p+1\tGPRLSGONYQIRFK-UHFFFAOYSA-N\t[H+]
    MNXM41\tD-glucose\tchebi:17234\tC6H12O6\t0\t180.063\tInChI=1S/C6H12O6\tWQZGKKKJIJFFOK-GASJEMHNSA-N\tOC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O
    BIOMASS\tBIOMASS\tmnx:BIOMASS\t\t\t\t\t\t
""")


def test_build_compounds_table(tmp_path):
    chem_prop = tmp_path / "chem_prop.tsv"
    chem_prop.write_text(CHEM_PROP_FIXTURE)
    db_path = tmp_path / "resolver.db"

    conn = sqlite3.connect(db_path)
    bmr.build_compounds_table(conn, chem_prop)
    conn.commit()

    cur = conn.cursor()
    cur.execute("SELECT mnxm_id, name, reference, formula, charge, mass, inchi, inchikey, smiles FROM compounds ORDER BY mnxm_id")
    rows = cur.fetchall()

    # 3 rows expected (BIOMASS, MNXM01, MNXM41), but BIOMASS has empty fields
    assert len(rows) == 3
    biomass = [r for r in rows if r[0] == "BIOMASS"][0]
    assert biomass[1] == "BIOMASS"
    assert biomass[3] == ""  # empty formula

    pmf = [r for r in rows if r[0] == "MNXM01"][0]
    assert pmf[1] == "PMF"
    assert pmf[2] == "mnx:PMF"
    assert pmf[3] == "H"
    assert pmf[4] == 1
    assert abs(pmf[5] - 1.00794) < 1e-5
    assert pmf[7] == "GPRLSGONYQIRFK-UHFFFAOYSA-N"

    glucose = [r for r in rows if r[0] == "MNXM41"][0]
    assert glucose[1] == "D-glucose"
    assert glucose[3] == "C6H12O6"
    assert abs(glucose[5] - 180.063) < 1e-3

    conn.close()


CHEM_XREF_FIXTURE = textwrap.dedent("""\
    #source\tID\tdescription
    MNXM41\tMNXM41\tD-glucose
    chebi:17234\tMNXM41\tD-glucose||Dextrose
    CHEBI:17234\tMNXM41\tD-glucose
    kegg.compound:C00031\tMNXM41\tD-glucose
    keggC:C00031\tMNXM41\tD-glucose
    bigg.metabolite:glc__D\tMNXM41\tD-glucose
    biggM:glc__D\tMNXM41\tD-glucose
    mnx:PMF\tMNXM01\tPMF
""")


def test_build_compound_aliases_normalizes_dual_prefixes(tmp_path):
    """Dual short/long source forms (CHEBI/chebi, keggC/kegg.compound) collapse
    to a single canonical prefix; values dedup per (source, value, mnxm)."""
    # First populate compounds (alias FK target)
    chem_prop = tmp_path / "chem_prop.tsv"
    chem_prop.write_text(CHEM_PROP_FIXTURE)
    chem_xref = tmp_path / "chem_xref.tsv"
    chem_xref.write_text(CHEM_XREF_FIXTURE)
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    bmr.build_compounds_table(conn, chem_prop)

    bmr.build_compound_aliases_table(conn, chem_xref)
    conn.commit()

    cur = conn.cursor()
    cur.execute("SELECT source, value, mnxm_id FROM compound_aliases ORDER BY source, value")
    rows = cur.fetchall()

    # CHEBI:17234 should be normalized to source='chebi', value='17234', and
    # both `chebi:17234` and `CHEBI:17234` rows should collapse to one
    chebi_rows = [r for r in rows if r[0] == "chebi"]
    assert len(chebi_rows) == 1
    assert chebi_rows[0] == ("chebi", "17234", "MNXM41")

    # kegg.compound:C00031 and keggC:C00031 collapse to one
    kegg_rows = [r for r in rows if r[0] == "kegg.compound"]
    assert kegg_rows == [("kegg.compound", "C00031", "MNXM41")]

    # bigg.metabolite/biggM collapse
    bigg_rows = [r for r in rows if r[0] == "bigg.metabolite"]
    assert bigg_rows == [("bigg.metabolite", "glc__D", "MNXM41")]

    # Self-references (bare MNXM41) skipped
    assert not any(r[0] == "" for r in rows)

    conn.close()


def test_build_compound_names_normalizes(tmp_path):
    """Names from `name` and description (||-split) get normalized + indexed."""
    chem_prop = tmp_path / "chem_prop.tsv"
    chem_prop.write_text(CHEM_PROP_FIXTURE)
    chem_xref = tmp_path / "chem_xref.tsv"
    chem_xref.write_text(CHEM_XREF_FIXTURE)
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    bmr.build_compounds_table(conn, chem_prop)
    bmr.build_compound_names_table(conn, chem_xref)
    conn.commit()

    cur = conn.cursor()
    cur.execute("SELECT DISTINCT name_normalized, mnxm_id FROM compound_names "
                "WHERE mnxm_id = 'MNXM41' ORDER BY name_normalized")
    rows = cur.fetchall()
    names = {r[0] for r in rows}
    assert "d-glucose" in names
    assert "dextrose" in names

    conn.close()


REAC_PROP_FIXTURE = textwrap.dedent("""\
    #ID\tmnx_equation\treference\tclassifs\tis_balanced\tis_transport
    MNXR101234\t1 MNXM3@MNXD1 + 1 MNXM41@MNXD1 = 1 MNXM7@MNXD1 + 1 MNXM58@MNXD1\tmnx:MNXR101234\t2.7.1.1\tB\t
    MNXR02\t1 MNXM1@MNXD1 = 1 MNXM1@MNXD2\tmnx:MNXR02\t\tB\tT
    EMPTY\t = \tmnx:EMPTY\t\tB\t
""")


REAC_XREF_FIXTURE = textwrap.dedent("""\
    #source\tID\tdescription
    MNXR101234\tMNXR101234\thexokinase
    kegg.reaction:R00299\tMNXR101234\thexokinase
    keggR:R00299\tMNXR101234\thexokinase
    rhea:16332\tMNXR101234\thexokinase
    rh:16332\tMNXR101234\thexokinase
""")


def test_build_reactions_table(tmp_path):
    """Reactions table mirrors reac_prop.tsv 6 columns."""
    reac_prop = tmp_path / "reac_prop.tsv"
    reac_prop.write_text(REAC_PROP_FIXTURE)
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)

    bmr.build_reactions_table(conn, reac_prop)
    conn.commit()

    cur = conn.cursor()
    cur.execute(
        "SELECT mnxr_id, mnx_equation, reference, classifs, is_balanced, is_transport "
        "FROM reactions ORDER BY mnxr_id"
    )
    rows = cur.fetchall()
    assert len(rows) == 3

    # Look up MNXR101234
    hk = [r for r in rows if r[0] == "MNXR101234"][0]
    assert "MNXM41" in hk[1]
    assert "MNXM7" in hk[1]
    assert hk[3] == "2.7.1.1"
    assert hk[4] == "B"
    assert hk[5] in (None, "")

    # MNXR02 is a transport reaction
    transport = [r for r in rows if r[0] == "MNXR02"][0]
    assert transport[5] == "T"

    conn.close()


def test_build_reaction_aliases_normalizes(tmp_path):
    """Reaction xrefs canonicalize keggR→kegg.reaction and rh→rhea."""
    reac_xref = tmp_path / "reac_xref.tsv"
    reac_xref.write_text(REAC_XREF_FIXTURE)
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)

    bmr.build_reaction_aliases_table(conn, reac_xref)
    conn.commit()

    cur = conn.cursor()
    cur.execute("SELECT source, value, mnxr_id FROM reaction_aliases ORDER BY source, value")
    rows = cur.fetchall()

    kegg = [r for r in rows if r[0] == "kegg.reaction"]
    assert kegg == [("kegg.reaction", "R00299", "MNXR101234")]

    rhea = [r for r in rows if r[0] == "rhea"]
    assert rhea == [("rhea", "16332", "MNXR101234")]

    # Self-reference dropped
    assert not any(r[0] == "" for r in rows)

    conn.close()


def test_resolver_has_mnx_side_indexes(tmp_path):
    """Spec 1.2: resolver DB must have indexes on mnxm_id / mnxr_id columns."""
    chem_xref = tmp_path / "chem_xref.tsv"
    chem_xref.write_text(
        "#xref\tmnxm_id\tdescription\n"
        "chebi:17234\tMNXM41\tD-glucose\n"
    )
    reac_xref = tmp_path / "reac_xref.tsv"
    reac_xref.write_text(
        "#xref\tmnxr_id\n"
        "kegg.reaction:R00200\tMNXR101234\n"
    )

    db_path = tmp_path / "resolver.db"
    conn = sqlite3.connect(db_path)
    bmr.build_compound_aliases_table(conn, chem_xref)
    bmr.build_reaction_aliases_table(conn, reac_xref)
    conn.commit()

    cur = conn.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='index'")
    indexes = {row[0] for row in cur.fetchall()}

    assert "idx_compound_aliases_mnxm" in indexes
    assert "idx_reaction_aliases_mnxr" in indexes
    conn.close()


def test_build_main_end_to_end(tmp_path, monkeypatch):
    """main() wires all builders + writes diagnostic report. Synthetic cache."""
    cache = tmp_path / "cache" / "data"
    (cache / "mnx").mkdir(parents=True)

    # Synthesize the MNX input files
    (cache / "mnx" / "chem_prop.tsv").write_text(CHEM_PROP_FIXTURE)
    (cache / "mnx" / "chem_xref.tsv").write_text(CHEM_XREF_FIXTURE)
    (cache / "mnx" / "reac_prop.tsv").write_text(REAC_PROP_FIXTURE)
    (cache / "mnx" / "reac_xref.tsv").write_text(REAC_XREF_FIXTURE)

    monkeypatch.chdir(tmp_path)
    bmr.main(force=True)

    # MNX outputs exist
    assert (cache / "mnx" / "metabolite_resolver.db").exists()
    assert (cache / "mnx" / "metabolite_id_mapping_report.json").exists()

    # Diagnostic report has expected keys
    report = json.loads((cache / "mnx" / "metabolite_id_mapping_report.json").read_text())
    assert report["mnx_release"] == "MNXref 4.5 (2025-08-13)"
    assert report["compound_count"] == 3
    assert report["reaction_count"] == 3
