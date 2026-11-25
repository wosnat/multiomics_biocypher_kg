"""
Every property in the live graph is declared in the schema the graph publishes.

`Schema_info.schema_info` is the JSON dump of `config/schema_config.yaml`
(`bc.write_schema_info(as_node=True)`), and the explorer reads it as the KG's
schema. Properties the post-import scripts compute (`gene_count`,
`rank_by_effect`, `is_uninformative`, ...) never pass through an adapter, so
nothing forces them into the YAML; this test does. Declaring one costs an empty
CSV column that post-import later fills.

Fails when a node label or relationship type carries a property that its schema
entry does not declare, or declares with a different type, and when a declared
property is set on no node or edge at all (unless allowlisted as optional below).
"""

import json

import pytest


pytestmark = pytest.mark.kg


# Set by BioCypher itself, not by the schema.
BIOCYPHER_PROPS = {"id", "preferred_id", "source", "licence", "version"}

# Labels outside schema_config.yaml.
UNSCHEMED_LABELS = {"Schema_info"}

# Declared properties a build may legitimately leave unset everywhere. Each
# needs a reason; anything else declared but never set should be removed from
# schema_config.yaml, since consumers read the schema as the graph's contents.
OPTIONAL_PROPS = {
    # Only with an embedding file passed to the --go / --ec ontology adapters.
    ("BiologicalProcess", "anc2vec_embedding"),
    ("CellularComponent", "anc2vec_embedding"),
    ("MolecularFunction", "anc2vec_embedding"),
    ("EcNumber", "rxnfp_embedding"),
    # Only when a paperconfig metric reports p-values.
    ("DerivedMetric", "p_value_threshold"),
    ("Derived_metric_quantifies_gene", "p_value"),
    ("Derived_metric_quantifies_gene", "adjusted_p_value"),
    ("Derived_metric_quantifies_gene", "significant"),
    # Only when a gene_clusters entry sets p_value_col.
    ("Gene_in_gene_cluster", "p_value"),
    # Only when the source record carries one (PDF extraction / UniProt).
    ("Publication", "pmid"),
    ("Protein", "subcellular_location"),
}

# Neo4j property type -> schema_config type spellings.
TYPE_MAP = {
    "String": {"str"},
    "Long": {"int"},
    "Double": {"float"},
    "StringArray": {"str[]", "list[str]"},
    "LongArray": {"int[]"},
    "DoubleArray": {"float[]"},
}


def _label(key, entry):
    """BioCypher's Neo4j label for a schema entry."""
    label_as_edge = entry.get("label_as_edge")
    if label_as_edge:
        return label_as_edge[0].upper() + label_as_edge[1:]
    return "".join(w[:1].upper() + w[1:] for w in key.replace("_", " ").split())


@pytest.fixture(scope="module")
def declared(run_query):
    rows = run_query("MATCH (s:Schema_info) RETURN s.schema_info AS schema")
    assert rows, "no Schema_info node"
    schema = json.loads(rows[0]["schema"])
    out = {}
    for key, entry in schema.items():
        if not isinstance(entry, dict) or "represented_as" not in entry:
            continue
        out.setdefault(_label(key, entry), {}).update(entry.get("properties") or {})
    return out


@pytest.fixture(scope="module")
def live(run_query):
    """(label, property) -> set of Neo4j types, nodes and relationships."""
    out = {}
    for r in run_query(
        "CALL db.schema.nodeTypeProperties() "
        "YIELD nodeLabels, propertyName, propertyTypes "
        "RETURN nodeLabels AS labels, propertyName AS prop, propertyTypes AS types"
    ):
        for label in r["labels"]:
            if r["prop"] is not None:
                out.setdefault((label, r["prop"]), set()).update(r["types"])
    for r in run_query(
        "CALL db.schema.relTypeProperties() "
        "YIELD relType, propertyName, propertyTypes "
        "RETURN relType AS label, propertyName AS prop, propertyTypes AS types"
    ):
        if r["prop"] is not None:
            label = r["label"].strip(":`")
            out.setdefault((label, r["prop"]), set()).update(r["types"])
    return out


def _schemed(live, declared):
    for (label, prop), types in sorted(live.items()):
        if label in UNSCHEMED_LABELS or prop in BIOCYPHER_PROPS:
            continue
        if label not in declared:
            continue  # extra BioCypher ancestor labels (e.g. NamedThing)
        yield label, prop, types


def test_every_label_has_a_schema_entry(live, declared):
    labels = {label for label, _ in live} - UNSCHEMED_LABELS
    # Ancestor labels BioCypher adds carry the same properties as the leaf
    # label, so a label is covered when some co-label is declared. Only
    # relationship types are single-labelled; check those.
    rel_labels = {label for label in labels if "_" in label}
    missing = sorted(rel_labels - set(declared))
    assert not missing, f"relationship types with no schema entry: {missing}"


def test_every_live_property_is_declared(live, declared):
    missing = [
        f"{label}.{prop}"
        for label, prop, _ in _schemed(live, declared)
        if prop not in declared[label]
    ]
    assert not missing, (
        f"{len(missing)} properties in the graph are not declared in "
        f"config/schema_config.yaml (declare post-import properties too): {missing}"
    )


def test_declared_types_match_live_types(live, declared):
    wrong = []
    for label, prop, types in _schemed(live, declared):
        want = declared[label].get(prop)
        if want is None:
            continue
        ok = set().union(*(TYPE_MAP.get(t, {t}) for t in types))
        if str(want) not in ok:
            wrong.append(f"{label}.{prop}: declared {want}, live {sorted(types)}")
    assert not wrong, f"schema types disagree with the graph: {wrong}"


def test_every_declared_property_is_set_somewhere(live, declared):
    labels = {label for label, _ in live}
    unset = [
        f"{label}.{prop}"
        for label, props in sorted(declared.items())
        if label in labels  # entities the build did not populate are out of scope
        for prop in sorted(props)
        if (label, prop) not in live and (label, prop) not in OPTIONAL_PROPS
    ]
    assert not unset, (
        f"{len(unset)} properties declared in config/schema_config.yaml are set on "
        f"no node or edge; remove them or add them to OPTIONAL_PROPS with a reason: {unset}"
    )
