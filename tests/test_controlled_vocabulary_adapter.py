from multiomics_kg.adapters.controlled_vocabulary_adapter import (
    ControlledVocabularyAdapter,
)


def _nodes():
    a = ControlledVocabularyAdapter()
    a.download_data()
    return {nid: props for nid, _label, props in a.get_nodes()}


def test_emits_one_node_per_vocabulary():
    nodes = _nodes()
    assert "Gene_catalyzes_ec_number.evidence" in nodes
    assert "InterproEntry.interpro_type" in nodes


def test_label_is_controlled_vocabulary():
    a = ControlledVocabularyAdapter()
    a.download_data()
    assert {label for _n, label, _p in a.get_nodes()} == {"controlled_vocabulary"}


def test_node_carries_the_declared_values():
    props = _nodes()["InterproEntry.interpro_type"]
    assert props["applies_to"] == "InterproEntry"
    assert props["applies_to_kind"] == "node"
    assert props["value_type"] == "string"
    assert "HOMOLOGOUS_SUPERFAMILY" in props["values"]
    assert props["closed"] == "true"          # R5: bool_string, not bool


def test_expected_empty_vocabulary_emits_empty_values():
    props = _nodes()["InterproEntry.level_kind"]
    assert props["values"] == []
    assert props["expected_empty"] == "true"


def test_score_vocabulary_publishes_its_signals():
    props = _nodes()["Gene_has_tcdb_family.evidence_score"]
    assert props["signal_count"] == 5
    assert "pfam_corroborated" in props["signals"]
    assert props["min_value"] == 0.0 and props["max_value"] == 1.0


def test_strings_are_sanitized():
    """No property value may contain ' or | (BioCypher CSV + array delimiter)."""
    for props in _nodes().values():
        for v in props.values():
            for s in (v if isinstance(v, list) else [v]):
                if isinstance(s, str):
                    assert "'" not in s and "|" not in s
