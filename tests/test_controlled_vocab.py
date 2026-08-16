import pytest
from multiomics_kg.utils.controlled_vocab import (
    VOCAB, load_vocabularies, vocabularies_hash,
)


def test_loads_the_shipped_config():
    entries = load_vocabularies()
    assert "Gene_catalyzes_ec_number.evidence" in entries
    e = entries["Gene_catalyzes_ec_number.evidence"]
    assert e.applies_to_kind == "edge"
    assert e.value_type == "string"
    assert set(e.values) == {"curated", "family_inferred"}


def test_evidence_domain_is_per_edge_type():
    """domain_inferred is not a possible value on the EC edge (spec §5.2)."""
    entries = load_vocabularies()
    assert "domain_inferred" not in entries["Gene_catalyzes_ec_number.evidence"].values
    assert "domain_inferred" in entries[
        "Gene_involved_in_biological_process.evidence"].values


def test_check_returns_declared_value():
    assert VOCAB.check("Gene_catalyzes_ec_number", "evidence", "curated") == "curated"


def test_check_raises_on_undeclared_value():
    with pytest.raises(ValueError, match="domain_inferred"):
        VOCAB.check("Gene_catalyzes_ec_number", "evidence", "domain_inferred")


def test_check_passes_through_open_vocabularies():
    """closed: false means enumerate at runtime — check must not reject."""
    assert VOCAB.check("Gene", "gene_category", "anything at all") == "anything at all"


def test_bool_value_type_is_rejected():
    """R5: native bool is not an admissible value_type."""
    import yaml, tempfile, pathlib
    bad = {"Bad.prop": {"applies_to": "Bad", "applies_to_kind": "node",
                        "property": "prop", "value_type": "bool",
                        "closed": True, "values": [], "description": "x"}}
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as f:
        yaml.safe_dump(bad, f)
        path = f.name
    with pytest.raises(ValueError, match="bool"):
        load_vocabularies(path)
    pathlib.Path(path).unlink()


def test_hash_is_stable_and_order_independent():
    entries = list(load_vocabularies().values())
    h1 = vocabularies_hash(entries)
    h2 = vocabularies_hash(list(reversed(entries)))
    assert h1 == h2
    assert h1.startswith("sha256:")
