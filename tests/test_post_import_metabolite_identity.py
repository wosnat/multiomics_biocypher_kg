"""Static checks: metabolite-identity post-import Cypher lives in BOTH scripts."""
import pathlib

import pytest

LITERALS = [
    "metaboliteClassFullText",
    "chemicalRoleFullText",
    "metabolite_structure_kind_idx",
    "metabolite_class_id_idx",
    "chemical_role_id_idx",
    "m.identifier_tokens",
    "[m.name, m.formula, m.identifier_tokens]",
    "c.direct_metabolite_count",
    "r.direct_transporter_count",
    "m.kegg_functional_categories",
    "m.kegg_brite_classes",
    "s.chebi_release",
]


@pytest.mark.parametrize("path", ["scripts/post-import.sh", "scripts/post-import.cypher"])
@pytest.mark.parametrize("literal", LITERALS)
def test_literal_present(path, literal):
    assert literal in pathlib.Path(path).read_text(), f"{path} missing {literal!r}"
