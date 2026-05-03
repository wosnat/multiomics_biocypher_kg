"""KG validity: PMM0001 has eggNOG seed_ortholog populated.

F4 ship 4.2 surfaces eggNOG seed-ortholog pointers as Gene properties.
Spot-check that the field lands with expected format on a known gene.
"""

import pytest

pytestmark = pytest.mark.kg


def test_med4_pmm0001_has_seed_ortholog(run_query):
    """Spot-check: PMM0001 should have eggNOG seed_ortholog = '59919.PMM0001'
    (self-match — MED4 taxid 59919 is in eggNOG's reference DB)."""
    result = run_query("""
        MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
        WHERE o.preferred_name = 'Prochlorococcus MED4' AND g.locus_tag = 'PMM0001'
        RETURN g.seed_ortholog AS so, g.seed_ortholog_evalue AS ev
    """)
    assert len(result) == 1, f"PMM0001 not found"
    assert result[0]["so"] == "59919.PMM0001", f"unexpected seed_ortholog: {result[0]['so']!r}"
    assert isinstance(result[0]["ev"], float), f"seed_ortholog_evalue should be float, got {type(result[0]['ev'])}"
    assert result[0]["ev"] < 1e-50, f"seed_ortholog_evalue should be very small for self-match, got {result[0]['ev']}"
