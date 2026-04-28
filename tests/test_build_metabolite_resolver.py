"""Unit tests for sub-step 7 resolver builder."""
from __future__ import annotations

import sqlite3
import textwrap
from pathlib import Path

from multiomics_kg.download import build_metabolite_resolver as bmr


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
