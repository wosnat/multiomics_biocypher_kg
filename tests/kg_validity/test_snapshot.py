"""
Snapshot regression tests for the multi-omics knowledge graph.

Loads a previously generated snapshot (snapshot_data.json) and verifies that
every sampled node and edge still exists in the current graph with the expected
properties. This catches silent data loss from adapter changes.

Regenerate the snapshot after intentional changes:
    uv run python tests/kg_validity/generate_snapshot.py
"""

import json
from pathlib import Path

import pytest


pytestmark = pytest.mark.kg

SNAPSHOT_PATH = Path(__file__).parent / "snapshot_data.json"


def _load_snapshot():
    if not SNAPSHOT_PATH.exists():
        pytest.skip(f"Snapshot file not found: {SNAPSHOT_PATH}")
    with open(SNAPSHOT_PATH) as f:
        return json.load(f)


def _node_ids():
    """Generate (test_id, node_dict) pairs for parametrize."""
    try:
        data = json.loads(SNAPSHOT_PATH.read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return []
    return [
        (f"{n['label']}:{n['id']}", n)
        for n in data.get("nodes", [])
    ]


def _edge_ids():
    """Generate (test_id, edge_dict) pairs for parametrize."""
    try:
        data = json.loads(SNAPSHOT_PATH.read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return []
    return [
        (f"{e['type']}:{e['source']}->{e['target']}", e)
        for e in data.get("edges", [])
    ]


# ---------------------------------------------------------------------------
# Node existence
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("test_id,node", _node_ids(), ids=lambda x: x if isinstance(x, str) else "")
def test_snapshot_node_exists(run_query, test_id, node):
    """Verify that a previously observed node still exists with expected properties."""
    result = run_query(
        f"MATCH (n:{node['label']}) WHERE n.id = $nid RETURN n",
        nid=node["id"],
    )
    assert len(result) >= 1, (
        f"Snapshot node missing from graph: {node['label']} with id={node['id']}"
    )
    actual = dict(result[0]["n"])
    for prop, expected_val in node.get("properties", {}).items():
        actual_val = actual.get(prop)
        assert actual_val == expected_val, (
            f"Property mismatch on {node['label']} id={node['id']}: "
            f"{prop}={actual_val!r}, expected {expected_val!r}"
        )


# ---------------------------------------------------------------------------
# Edge existence
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("test_id,edge", _edge_ids(), ids=lambda x: x if isinstance(x, str) else "")
def test_snapshot_edge_exists(run_query, test_id, edge):
    """Verify that a previously observed edge still exists between the expected nodes."""
    result = run_query(
        f"MATCH (src)-[r:`{edge['type']}`]->(tgt) "
        f"WHERE src.id = $src AND tgt.id = $tgt "
        f"RETURN properties(r) AS props",
        src=edge["source"],
        tgt=edge["target"],
    )
    assert len(result) >= 1, (
        f"Snapshot edge missing from graph: {edge['type']} "
        f"from {edge['source']} -> {edge['target']}"
    )
    actual = result[0]["props"]
    for prop, expected_val in edge.get("properties", {}).items():
        actual_val = actual.get(prop)
        if isinstance(expected_val, float):
            assert actual_val is not None and abs(actual_val - expected_val) < 1e-6, (
                f"Property mismatch on {edge['type']} "
                f"{edge['source']}->{edge['target']}: "
                f"{prop}={actual_val!r}, expected {expected_val!r}"
            )
        else:
            assert actual_val == expected_val, (
                f"Property mismatch on {edge['type']} "
                f"{edge['source']}->{edge['target']}: "
                f"{prop}={actual_val!r}, expected {expected_val!r}"
            )
