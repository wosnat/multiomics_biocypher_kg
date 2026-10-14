import pytest

from tests.metabolite_identity_helpers import cpd, make_chebi, term
from multiomics_kg.utils.metabolite_identity.structure import apply_structure, classify_structure

CHEBI = make_chebi(
    CHEBI_17522=term("alditol", gformula="C2H6O2(CH2O)n"),
    CHEBI_17754=term("glycerol", charge=0, kegg=["C00116"], inchikey="PEDCQBHIVMGVHV-UHFFFAOYSA-N"),
    CHEBI_37684=term("mannose", gformula="C6H12O6"))


def test_classify():
    assert classify_structure(cpd("Rubidium cation", formula="Rb"), None, CHEBI) == "defined"
    assert classify_structure(cpd("x", formula="C5H8O2R"), None, CHEBI) == "generic"
    assert classify_structure(cpd("PHB", formula="(C4H6O2)n"), None, CHEBI) == "polymer"
    assert classify_structure(cpd("(1,4-alpha-D-glucosyl)n"), None, CHEBI) == "polymer"
    assert classify_structure(cpd("alditol"), "CHEBI:17522", CHEBI) == "polymer"
    assert classify_structure(cpd("thing"), None, CHEBI) == "defined"


def test_b2_scrub_by_class_and_by_shared_inchikey():
    ik = "PEDCQBHIVMGVHV-UHFFFAOYSA-N"
    table = {"kegg.compound:C00116": cpd("Glycerol", formula="C3H8O3", mass=92.09, inchikey=ik),
             "chebi:17522": cpd("alditol", formula="C3H8O3", mass=92.09, inchikey=ik),
             "kegg.compound:C00734": cpd("Chitosan", formula="C6H13NO5", mass=179.2,
                                         inchikey="MSWZFWKMSRAUBD-IVMDWMLBSA-N"),
             "kegg.compound:C00329": cpd("Glucosamine", formula="C6H13NO5", mass=179.2,
                                         inchikey="MSWZFWKMSRAUBD-IVMDWMLBSA-N")}
    rep = apply_structure(table, CHEBI, {"chebi:17522": "CHEBI:17522"})
    alditol = table["chebi:17522"]
    assert alditol["structure_kind"] in ("generic", "polymer")
    assert alditol.get("inchikey") is None and alditol.get("mass") is None
    assert alditol["elements"] == ["C", "H", "O"]
    assert table["kegg.compound:C00116"]["inchikey"] == ik
    assert table["kegg.compound:C00734"]["structure_kind"] == "polymer"
    assert table["kegg.compound:C00734"].get("inchikey") is None
    assert table["kegg.compound:C00329"]["structure_kind"] == "defined"
    assert table["kegg.compound:C00329"]["inchikey"] == "MSWZFWKMSRAUBD-IVMDWMLBSA-N"
    assert set(rep["scrubbed"]) >= {"chebi:17522", "kegg.compound:C00734"}


def test_formula_fill_and_no_mass_on_generic():
    table = {"chebi:37684": cpd("mannose"), "kegg.compound:C99999": cpd("R-thing", formula="C2H4R", mass=12.0)}
    apply_structure(table, CHEBI, {"chebi:37684": "CHEBI:37684"})
    assert table["chebi:37684"]["formula"] == "C6H12O6" and table["chebi:37684"]["elements"] == ["C", "H", "O"]
    assert table["kegg.compound:C99999"].get("mass") is None


def test_leftover_inchikey_collision_fails():
    ik = "AAAAAAAAAAAAAA-BBBBBBBBBB-N"
    table = {"kegg.compound:C1": cpd("a", formula="CH4", inchikey=ik),
             "kegg.compound:C2": cpd("b", formula="CH4", inchikey=ik)}
    with pytest.raises(ValueError, match="InChIKey"):
        apply_structure(table, CHEBI, {})
