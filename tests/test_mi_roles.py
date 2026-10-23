from tests.metabolite_identity_helpers import make_chebi, term
from multiomics_kg.utils.metabolite_identity.roles import (
    build_role_hierarchy, build_system_roles, partition_role_strings)

CHEBI = make_chebi(
    CHEBI_50906=term("role"),
    CHEBI_24432=term("biological role", is_a=["CHEBI:50906"]),
    CHEBI_51086=term("chemical role", is_a=["CHEBI:50906"]),
    CHEBI_33232=term("application", is_a=["CHEBI:50906"]),
    CHEBI_26672=term("siderophore", is_a=["CHEBI:24432"]),
    CHEBI_38155=term("phytosiderophore", is_a=["CHEBI:26672"]),
    CHEBI_33229=term("vitamin (role)", alt=["CHEBI:10004"], is_a=["CHEBI:24432", "CHEBI:33232"]),
    CHEBI_17632=term("nitrate"))


def test_partition_detects_secondary_id_roles():
    mols, roles = partition_role_strings(
        {"CHEBI:10004;vitamin", "CHEBI:38155;phytosiderophore", "CHEBI:17632;nitrate"}, CHEBI)
    assert mols == {"CHEBI:17632;nitrate"}
    assert roles == {"CHEBI:10004;vitamin": "chebi.role:33229",
                     "CHEBI:38155;phytosiderophore": "chebi.role:38155"}


def test_hierarchy_ancestors_no_root_levels_and_identifiers():
    h = build_role_hierarchy({"chebi.role:38155", "chebi.role:33229"}, CHEBI)
    assert "chebi.role:50906" not in h
    assert h["chebi.role:24432"]["level"] == 0 and h["chebi.role:24432"]["level_kind"] == "chebi_role_branch"
    assert h["chebi.role:26672"]["level"] == 1 and h["chebi.role:38155"]["level"] == 2
    assert h["chebi.role:38155"]["parents"] == ["chebi.role:26672"]
    assert sorted(h["chebi.role:33229"]["parents"]) == ["chebi.role:24432", "chebi.role:33232"]
    assert h["chebi.role:33229"]["all_identifiers"] == ["chebi.role:33229", "chebi:10004", "chebi:33229"]
    assert all(p in h for v in h.values() for p in v["parents"])


def test_system_roles_sorted_nonempty_only():
    sub = build_system_roles({"2.A.1.1.1": {"CHEBI:10004;vitamin", "CHEBI:17632;nitrate"}, "3.A.1.1.1": {"CHEBI:17632;nitrate"}},
                              {"CHEBI:10004;vitamin": "chebi.role:33229"})
    assert sub == {"2.A.1.1.1": ["chebi.role:33229"]}
