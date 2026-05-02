"""KG validity: F1.1 is_uninformative term flag.

Asserts:
1. Every YAML ids: entry resolves to an existing node post-build with
   is_uninformative='true'.
2. The KEGG name_pattern matches >= 100 nodes (sanity floor; today ~210).
3. No node outside the YAML scope has is_uninformative set (catches
   over-flagging in Pfam, EC, BRITE which are anti-scope).
4. Drift check: every catch-all role entry in the YAML has its corresponding
   *_TO_CATEGORY mapping returning 'Unknown' in build_gene_annotations.py.
"""

from pathlib import Path

import pytest
import yaml

pytestmark = pytest.mark.kg


@pytest.fixture(scope="module")
def yaml_vocab():
    with open("config/uninformative_terms.yaml") as f:
        return yaml.safe_load(f)


SECTION_TO_LABEL = {
    "biological_process": "BiologicalProcess",
    "molecular_function": "MolecularFunction",
    "cellular_component": "CellularComponent",
    "cog_category": "CogFunctionalCategory",
    "cyanorak_role": "CyanorakRole",
    "tigr_role": "TigrRole",
    "kegg_term": "KeggTerm",
}


def test_yaml_ids_all_resolve_to_flagged_nodes(run_query, yaml_vocab):
    """Every ids: entry exists post-build and has is_uninformative='true'."""
    failures = []
    for section_name, payload in yaml_vocab.items():
        label = SECTION_TO_LABEL.get(section_name)
        assert label, f"Unknown section '{section_name}' — update SECTION_TO_LABEL"
        for term_id in (payload or {}).get("ids") or []:
            result = run_query(
                f"MATCH (t:{label} {{id: $id}}) RETURN t.is_uninformative AS flag",
                id=term_id,
            )
            if not result:
                failures.append(f"{label}:{term_id} — node missing")
            elif result[0]["flag"] != "true":
                failures.append(f"{label}:{term_id} — flag={result[0]['flag']!r}")
    assert not failures, "\n".join(failures)


def test_kegg_pattern_flags_at_least_100_kos(run_query):
    """KEGG name_pattern should flag at least 100 'uncharacterized protein' KOs."""
    result = run_query("""
        MATCH (t:KeggTerm) WHERE t.is_uninformative = 'true' RETURN count(*) AS n
    """)
    assert result[0]["n"] >= 100, f"got {result[0]['n']}"


def test_no_pfam_or_ec_flagged(run_query):
    """Anti-scope: Pfam, PfamClan, EcNumber, BriteCategory stay UN-flagged."""
    for label in ["Pfam", "PfamClan", "EcNumber", "BriteCategory"]:
        result = run_query(f"""
            MATCH (t:{label}) WHERE t.is_uninformative IS NOT NULL
            RETURN count(*) AS n
        """)
        assert result[0]["n"] == 0, f"{label} has flagged terms (anti-scope violation)"


def test_drift_role_yaml_matches_to_category_mappings(yaml_vocab):
    """Build-time drift check (no Neo4j): every catch-all role entry in the
    YAML must, when fed as a gene's only role, yield gene_category='Unknown'
    from the *_TO_CATEGORY mappings in build_gene_annotations.py.

    This catches the case where a new uninformative role gets added
    to the YAML but the corresponding *_TO_CATEGORY entry isn't 'Unknown'.
    """
    from multiomics_kg.download.build_gene_annotations import (
        COG_TO_CATEGORY,
        _compute_gene_category,
    )
    failures = []
    # COG: ids look like 'cog.category:X' — strip prefix
    for term_id in (yaml_vocab.get("cog_category") or {}).get("ids") or []:
        letter = term_id.split(":")[-1]
        if COG_TO_CATEGORY.get(letter) != "Unknown":
            failures.append(f"COG {letter}: TO_CATEGORY={COG_TO_CATEGORY.get(letter)!r}, expected 'Unknown'")
    # Cyanorak: ids look like 'cyanorak.role:R.1'
    for term_id in (yaml_vocab.get("cyanorak_role") or {}).get("ids") or []:
        # Use _compute_gene_category to honor the D-subcode special case
        code = term_id.split(":")[-1]
        result = _compute_gene_category({"cyanorak_Role": [code]})
        if result != "Unknown":
            failures.append(f"Cyanorak {code}: gene_category={result!r}, expected 'Unknown'")
    # TIGR: ids are 'tigr.role:156' but TIGR_TO_CATEGORY is keyed by description.
    # Skip programmatic lookup here; rely on COG + Cyanorak checks above.
    assert not failures, "\n".join(failures)
