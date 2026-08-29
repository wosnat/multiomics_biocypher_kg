import json

from multiomics_kg.utils.tigr_roles import (
    EQUIVALOG_TYPES,
    inferred_roles,
    load_tigr_roles,
    mainrole_slug,
    role_name,
    split_role_name,
)

TIGR_ROLES = {
    "release": "TIGRFAMs 15.0 (frozen 2018)",
    "roles": {
        "132": {"mainrole": "DNA metabolism", "sub1role": "DNA replication, recombination, and repair"},
        "120": {"mainrole": "Energy metabolism", "sub1role": "TCA cycle"},
        "156": {"mainrole": "Hypothetical proteins", "sub1role": "Conserved"},
    },
    "family_role": {"TIGR00001": ["132"], "TIGR00002": ["120"], "TIGR00003": ["156"],
                    "TIGR00004": ["120"], "TIGR00005": ["120", "132"]},
}
NCBIFAM_REF = {
    "TIGR00001": {"name": "a", "family_type": "equivalog"},
    "TIGR00002": {"name": "b", "family_type": "subfamily"},
    "TIGR00003": {"name": "c", "family_type": "equivalog"},
    "TIGR00004": {"name": "d", "family_type": "equivalog"},
    "TIGR00005": {"name": "e", "family_type": "equivalog"},
}


def test_equivalog_types_is_exactly_equivalog():
    assert EQUIVALOG_TYPES == frozenset({"equivalog"})


def test_mainrole_slug():
    assert mainrole_slug("Energy metabolism") == "energy_metabolism"
    assert mainrole_slug("Purines, pyrimidines, nucleosides, and nucleotides") == \
        "purines_pyrimidines_nucleosides_and_nucleotides"
    assert mainrole_slug("Biosynthesis of cofactors, prosthetic groups, and carriers") == \
        "biosynthesis_of_cofactors_prosthetic_groups_and_carriers"


def test_role_name_and_split_round_trip():
    assert role_name("120", TIGR_ROLES) == "Energy metabolism / TCA cycle"
    assert split_role_name("Energy metabolism / TCA cycle") == ("Energy metabolism", "TCA cycle")
    assert split_role_name("Not Found") == ("Not Found", None)


def test_inferred_roles_equivalog_gate():
    out = inferred_roles(["TIGR00001", "TIGR00002", "NF000001"], TIGR_ROLES, NCBIFAM_REF)
    assert out == {"132": ["TIGR00001"]}          # TIGR00002 is subfamily → dropped


def test_inferred_roles_groups_supporting_accessions():
    out = inferred_roles(["TIGR00004", "TIGR00002"], TIGR_ROLES, NCBIFAM_REF)
    assert out == {"120": ["TIGR00004"]}


def test_inferred_roles_fans_out_multi_role_family():
    assert inferred_roles(["TIGR00005"], TIGR_ROLES, NCBIFAM_REF) == {"120": ["TIGR00005"], "132": ["TIGR00005"]}


def test_inferred_roles_keeps_junk_roles():
    assert inferred_roles(["TIGR00003"], TIGR_ROLES, NCBIFAM_REF) == {"156": ["TIGR00003"]}


def test_inferred_roles_empty_inputs():
    assert inferred_roles(None, TIGR_ROLES, NCBIFAM_REF) == {}
    assert inferred_roles(["TIGR00001"], TIGR_ROLES, {}) == {}   # no ref → no family_type → no edge


def test_load_tigr_roles(tmp_path):
    assert load_tigr_roles(tmp_path) is None
    (tmp_path / "ncbifam").mkdir()
    (tmp_path / "ncbifam" / "tigr_roles.json").write_text(json.dumps(TIGR_ROLES))
    assert load_tigr_roles(tmp_path) == TIGR_ROLES
