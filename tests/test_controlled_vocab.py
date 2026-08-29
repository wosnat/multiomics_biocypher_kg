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


def test_every_list_filter_values_filter_is_declared():
    """The 8 filters the MCP serves must all be declared, closed or open."""
    from multiomics_kg.utils.controlled_vocab import load_vocabularies
    declared = load_vocabularies()
    required = [
        "Gene.gene_category", "BriteCategory.tree", "Experiment.growth_phases",
        "DerivedMetric.metric_type", "DerivedMetric.value_kind",
        "Experiment.compartment", "Experiment.omics_type",
        "Metabolite.evidence_sources",
    ]
    missing = [k for k in required if k not in declared]
    assert not missing, f"list_filter_values filters not declared: {missing}"


# ── KG-SYNC-005: every gene→ontology edge type declares sources + evidence ─────

_GENE_ONTOLOGY_EDGES = [
    "Gene_involved_in_biological_process", "Gene_enables_molecular_function",
    "Gene_located_in_cellular_component", "Gene_catalyzes_ec_number",
    "Gene_has_pfam", "Gene_has_cazy_family", "Gene_has_tcdb_family",
    "Gene_has_kegg_ko", "Gene_in_cog_category", "Gene_has_cyanorak_role",
    "Gene_has_tigr_role", "Gene_has_interpro_entry", "Gene_has_ncbifam_family",
    "Gene_has_merops_family",
]
_LADDER = ["curated", "signature", "homology", "family_inferred", "domain_inferred"]


def test_every_gene_ontology_edge_declares_sources_and_evidence():
    from multiomics_kg.utils.controlled_vocab import load_vocabularies
    vocab = load_vocabularies()
    for edge in _GENE_ONTOLOGY_EDGES:
        assert f"{edge}.sources" in vocab, edge
        ev = vocab.get(f"{edge}.evidence")
        assert ev is not None, edge
        assert ev.closed
        assert set(ev.values) <= set(_LADDER), (edge, ev.values)


def test_evidence_ladder_placements():
    from multiomics_kg.utils.controlled_vocab import load_vocabularies
    vocab = load_vocabularies()
    assert vocab["Gene_has_tcdb_family.evidence"].values == ["family_inferred", "homology"]
    assert vocab["Gene_has_merops_family.evidence"].values == ["homology"]
    assert vocab["Gene_has_interpro_entry.evidence"].values == ["signature"]
    assert vocab["Gene_has_cyanorak_role.evidence"].values == ["curated"]
    assert vocab["Gene_has_kegg_ko.evidence"].values == ["family_inferred"]
    ms = vocab["Gene_has_merops_family.evidence_score"]
    assert ms.signal_count == 2 and sorted(ms.signals) == ["pfam_support", "tier_le_2"]
    assert vocab["Gene_has_tcdb_family.attachment_depth"].values == ["most_specific", "superseded"]
    assert "retired" in vocab["NcbifamFamily.family_type"].values
    assert "MeropsFamily.family_class" in vocab and "MeropsFamily.family_type" not in vocab


# --- value_descriptions (explorer B1, 2026-08-29) ---------------------------

_VD_BASE = """
Edge.evidence:
  applies_to: Edge
  applies_to_kind: edge
  property: evidence
  value_type: string
  closed: true
  values: [curated, signature]
  description: x
"""


def _write(tmp_path, extra):
    p = tmp_path / "v.yaml"
    p.write_text(_VD_BASE + extra)
    return p


def test_value_descriptions_are_loaded_and_whitespace_normalized(tmp_path):
    p = _write(tmp_path, """  value_descriptions:
    curated: >
      hand-curated
      by a person
    signature: direct HMM hit
""")
    e = load_vocabularies(p)["Edge.evidence"]
    assert e.value_descriptions == {"curated": "hand-curated by a person",
                                    "signature": "direct HMM hit"}


def test_value_descriptions_default_to_empty(tmp_path):
    e = load_vocabularies(_write(tmp_path, ""))["Edge.evidence"]
    assert e.value_descriptions == {}


def test_value_descriptions_reject_undeclared_value(tmp_path):
    p = _write(tmp_path, """  value_descriptions:
    curated: a
    signature: b
    homology: not a declared value
""")
    with pytest.raises(ValueError, match="homology"):
        load_vocabularies(p)


def test_value_descriptions_on_closed_vocab_must_cover_every_value(tmp_path):
    p = _write(tmp_path, """  value_descriptions:
    curated: a
""")
    with pytest.raises(ValueError, match="signature"):
        load_vocabularies(p)


def test_value_descriptions_do_not_affect_the_hash(tmp_path):
    without = list(load_vocabularies(_write(tmp_path, "")).values())
    with_ = list(load_vocabularies(_write(tmp_path, """  value_descriptions:
    curated: a
    signature: b
""")).values())
    assert vocabularies_hash(without) == vocabularies_hash(with_)


def test_trust_slice_carries_value_descriptions():
    """Every trust vocabulary the explorer filters on is described per value."""
    entries = load_vocabularies()
    trust_props = {"evidence", "sources", "call_class", "best_hit_kind",
                   "attachment_depth", "substrate_depth", "pfam_support",
                   "go_support", "source_agreement", "detection_status",
                   "table_scope", "annotation_state"}
    missing = [e.id for e in entries.values()
               if e.property in trust_props and e.closed and not e.value_descriptions]
    assert not missing, missing


def test_sources_descriptions_do_not_drift_from_gene_annotations_config():
    """Every described `sources` value is a logical_sources id, and every
    source any `sources` vocabulary declares is described somewhere."""
    import yaml
    cfg = yaml.safe_load(open("config/gene_annotations_config.yaml"))
    logical = {ls["id"] for src in cfg["sources"].values()
               for ls in src["logical_sources"]}
    entries = [e for e in load_vocabularies().values() if e.property == "sources"]
    described = {v for e in entries for v in e.value_descriptions}
    declared = {v for e in entries for v in e.values}
    assert described <= logical, described - logical
    assert declared <= described, declared - described


# ── presence markers (backlog 2026-08-29): the two sentinel-or-absent flags ──

_UNINFORMATIVE_LABELS = [
    "BiologicalProcess", "MolecularFunction", "CellularComponent",
    "CogFunctionalCategory", "CyanorakRole", "TigrRole", "KeggTerm",
    "InterproEntry", "NcbifamFamily",
]
_BEST_EFFORT_LABELS = ["BiologicalProcess", "MolecularFunction", "CellularComponent"]


@pytest.mark.parametrize("entry_id", (
    [f"{lbl}.is_uninformative" for lbl in _UNINFORMATIVE_LABELS]
    + [f"{lbl}.level_is_best_effort" for lbl in _BEST_EFFORT_LABELS]
))
def test_presence_markers_are_declared(entry_id):
    """The two R5 exceptions are declared with the presence-marker shape:
    closed, single value 'true', sparse (absence IS the negative state)."""
    e = load_vocabularies()[entry_id]
    assert e.value_type == "string"
    assert e.closed is True
    assert e.values == ["true"]
    assert e.sparse is True
    assert "absen" in e.value_descriptions["true"].lower()


def test_presence_marker_must_be_sparse(tmp_path):
    """A 'true'-only vocabulary that is not sparse is a dense stringified bool,
    which R5 forbids — the loader refuses it."""
    p = tmp_path / "v.yaml"
    p.write_text("""Node.flag:
  applies_to: Node
  applies_to_kind: node
  property: flag
  value_type: string
  closed: true
  values: ['true']
  description: x
""")
    with pytest.raises(ValueError, match="presence marker"):
        load_vocabularies(p)
