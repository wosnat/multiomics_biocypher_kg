"""Unit tests for sub-step 7 resolver builder."""
from __future__ import annotations

import json
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


TCDB_FAMILIES_FIXTURE = textwrap.dedent("""\
    1.A.1\tThe Voltage-gated Ion Channel (VIC) Superfamily
    1.A.10\tThe Glutamate-gated Ion Channel (GIC) Family of Neurotransmitter Receptors
""")


TCDB_SUPERFAMILIES_FIXTURE = textwrap.dedent("""\
    #TCID\tSubfamily\tFamily\tFam_abbreviation\tSuperfamily
    1.A.1.1.1\t1.A.1.1\t1.A.1\tVIC\tVIC Superfamily
    1.A.1.1.2\t1.A.1.1\t1.A.1\tVIC\tVIC Superfamily
    1.A.10.1.1\t1.A.10.1\t1.A.10\tGIC\tGIC Superfamily
""")


TCDB_SUBSTRATES_FIXTURE = textwrap.dedent("""\
    1.A.1.1.1\tCHEBI:29103;potassium(1+)
    1.A.10.1.1\tCHEBI:29987;glutamate(2-)|CHEBI:33709;amino acid
""")


def test_build_tcdb_hierarchy(tmp_path):
    """TCDB hierarchy JSON joins families + superfamilies + substrates."""
    fams = tmp_path / "families.tsv"
    fams.write_text(TCDB_FAMILIES_FIXTURE)
    supers = tmp_path / "superfamilies.tsv"
    supers.write_text(TCDB_SUPERFAMILIES_FIXTURE)
    subs = tmp_path / "substrates.tsv"
    subs.write_text(TCDB_SUBSTRATES_FIXTURE)
    out = tmp_path / "tcdb_hierarchy.json"

    bmr.build_tcdb_hierarchy(out, fams, supers, subs)

    h = json.loads(out.read_text())

    # Class (level 0) synthesized
    assert h["1"]["level"] == 0
    assert h["1"]["level_kind"] == "tc_class"
    assert h["1"]["parent"] is None

    # Subclass (level 1) synthesized
    assert h["1.A"]["level"] == 1
    assert h["1.A"]["parent"] == "1"

    # Family (level 2) — name from families.tsv
    assert h["1.A.1"]["level"] == 2
    assert "Voltage-gated Ion Channel" in h["1.A.1"]["name"]
    assert h["1.A.1"]["parent"] == "1.A"
    assert h["1.A.1"]["abbreviation"] == "VIC"

    # Subfamily (level 3) — derived from col 2 of superfamilies
    assert h["1.A.1.1"]["level"] == 3
    assert h["1.A.1.1"]["parent"] == "1.A.1"

    # Specificity (level 4) — substrate joined
    sp = h["1.A.1.1.1"]
    assert sp["level"] == 4
    assert sp["parent"] == "1.A.1.1"
    assert sp["substrate_classes"] == ["potassium(1+)"]
    assert sp["superfamily"] == "VIC Superfamily"

    # Multi-substrate
    glu = h["1.A.10.1.1"]
    assert set(glu["substrate_classes"]) == {"glutamate(2-)", "amino acid"}


EGGNOG_FIXTURE = textwrap.dedent("""\
    ##   eggnog-mapper
    #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
    WP_001.1\to\t1e-50\t100\tCOG1\t-\tE\tdesc\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tGT19\t-\t-
    WP_002.1\to\t1e-50\t100\tCOG1\t-\tE\tdesc\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tGH13_1,CBM48\t-\t-
    WP_003.1\to\t1e-50\t100\tCOG1\t-\tE\tdesc\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-
    WP_004.1\to\t1e-50\t100\tCOG1\t-\tE\tdesc\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tGH32\t-\t-
""")


def test_build_cazy_hierarchy(tmp_path):
    """CAZy hierarchy bootstrapped from raw eggNOG col 19, all 6 classes hardcoded."""
    eggnog = tmp_path / "MED4.emapper.annotations"
    eggnog.write_text(EGGNOG_FIXTURE)
    out = tmp_path / "cazy_hierarchy.json"

    bmr.build_cazy_hierarchy(out, [eggnog])

    h = json.loads(out.read_text())

    # All 6 classes present even if not all observed
    for cls in ("GH", "GT", "PL", "CE", "AA", "CBM"):
        assert h[cls]["level"] == 0
        assert h[cls]["level_kind"] == "cazy_class"
        assert h[cls]["parent"] is None
        assert h[cls]["class"] == cls

    # Observed families derived
    assert h["GT19"]["level"] == 1
    assert h["GT19"]["parent"] == "GT"
    assert h["GT19"]["class"] == "GT"

    assert h["GH13"]["level"] == 1
    assert h["GH13"]["parent"] == "GH"

    # Observed subfamily derived (and its parent family also present)
    assert h["GH13_1"]["level"] == 2
    assert h["GH13_1"]["level_kind"] == "cazy_subfamily"
    assert h["GH13_1"]["parent"] == "GH13"
    assert h["GH13_1"]["class"] == "GH"

    assert h["CBM48"]["parent"] == "CBM"
    assert h["GH32"]["parent"] == "GH"


def test_build_main_end_to_end(tmp_path, monkeypatch):
    """main() wires all builders + writes diagnostic report. Synthetic cache."""
    cache = tmp_path / "cache" / "data"
    (cache / "mnx").mkdir(parents=True)
    (cache / "tcdb").mkdir(parents=True)

    # Synthesize all the input files
    (cache / "mnx" / "chem_prop.tsv").write_text(CHEM_PROP_FIXTURE)
    (cache / "mnx" / "chem_xref.tsv").write_text(CHEM_XREF_FIXTURE)
    (cache / "mnx" / "reac_prop.tsv").write_text(REAC_PROP_FIXTURE)
    (cache / "mnx" / "reac_xref.tsv").write_text(REAC_XREF_FIXTURE)
    (cache / "tcdb" / "families.tsv").write_text(TCDB_FAMILIES_FIXTURE)
    (cache / "tcdb" / "superfamilies.tsv").write_text(TCDB_SUPERFAMILIES_FIXTURE)
    (cache / "tcdb" / "substrates.tsv").write_text(TCDB_SUBSTRATES_FIXTURE)

    # Synthesize a tiny eggNOG annotation for one strain
    strain_dir = cache / "Prochlorococcus" / "genomes" / "MED4" / "eggnog"
    strain_dir.mkdir(parents=True)
    (strain_dir / "MED4.emapper.annotations").write_text(EGGNOG_FIXTURE)

    monkeypatch.chdir(tmp_path)
    bmr.main(force=True)

    # All four outputs exist
    assert (cache / "mnx" / "metabolite_resolver.db").exists()
    assert (cache / "tcdb" / "tcdb_hierarchy.json").exists()
    assert (cache / "cazy" / "cazy_hierarchy.json").exists()
    assert (cache / "mnx" / "metabolite_id_mapping_report.json").exists()

    # Diagnostic report has expected keys
    report = json.loads((cache / "mnx" / "metabolite_id_mapping_report.json").read_text())
    assert report["mnx_release"] == "MNXref 4.5 (2025-08-13)"
    assert report["compound_count"] == 3
    assert report["reaction_count"] == 3
    assert report["tcdb_hierarchy_entry_count"] >= 7
    assert report["cazy_hierarchy_entry_count"] >= 6  # at least the 6 classes
