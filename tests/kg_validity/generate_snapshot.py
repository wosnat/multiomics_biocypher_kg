"""
Generate a snapshot of sample nodes and edges from the live Neo4j knowledge graph.

This snapshot is used by test_snapshot.py as a regression fixture — if any of the
sampled items go missing after a rebuild, the test fails.

Usage:
    uv run python tests/kg_validity/generate_snapshot.py
    uv run python tests/kg_validity/generate_snapshot.py --neo4j-url bolt://localhost:7688
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from neo4j import GraphDatabase


SNAPSHOT_PATH = Path(__file__).parent / "snapshot_data.json"

# Hand-picked anchor nodes that should always be present
ANCHOR_NODES = {
    "Gene": [
        "ncbigene:PMM1375",   # MED4 psbA
        "ncbigene:PMM0001",   # MED4 first gene
    ],
    "OrganismTaxon": [
        "insdc.gcf:GCF_000011465.1",  # MED4
    ],
    # GO root terms — always present as ancestors of any annotated GO term
    "BiologicalProcess": [
        "go:0008150",  # biological_process (root)
    ],
    "CellularComponent": [
        "go:0005575",  # cellular_component (root)
    ],
    "MolecularFunction": [
        "go:0003674",  # molecular_function (root)
    ],
    # EC top-level class — always present in full Expasy hierarchy
    "EcNumber": [
        "ec:1.-.-.-",  # Oxidoreductases (EC class 1)
    ],
}

# Key properties to capture per node type
NODE_PROPERTIES = {
    "Gene": ["locus_tag", "product", "gene_name", "gene_synonyms", "strand"],
    "Protein": ["protein_synonyms", "gene_names"],
    "OrganismTaxon": ["organism_name", "strain_name", "genus", "clade", "ncbi_taxon_id"],
    "Publication": ["doi", "pmid", "title"],
    "EnvironmentalCondition": ["name", "condition_type"],
    "Cyanorak_cluster": ["cluster_number"],
    # GO term node types (from functional_annotation_adapter.py)
    "BiologicalProcess": ["name"],
    "CellularComponent": ["name"],
    "MolecularFunction": ["name"],
    # EC number nodes (from ec_adapter / functional_annotation_adapter)
    "EcNumber": ["name"],
}

# Key properties to capture per edge type
EDGE_PROPERTIES = {
    "Gene_belongs_to_organism": [],
    "Protein_belongs_to_organism": [],
    "Gene_encodes_protein": [],
    "Gene_in_cyanorak_cluster": [],
    "Gene_is_homolog_of_gene": ["cluster_id", "distance"],
    "Affects_expression_of": [
        "expression_direction", "log2_fold_change", "adjusted_p_value",
        "control_condition",
    ],
    "Affects_expression_of_homolog": [
        "expression_direction", "log2_fold_change", "original_gene",
        "homology_cluster_id",
    ],
    # Gene → GO edges (functional_annotation_adapter.py; label_as_edge values)
    "Gene_involved_in_biological_process": [],
    "Gene_located_in_cellular_component": [],
    "Gene_enables_molecular_function": [],
    # GO-GO hierarchy edges (label_as_edge in schema; BioCypher capitalizes first letter)
    "Biological_process_is_a_biological_process": [],
    "Cellular_component_is_a_cellular_component": [],
    "Molecular_function_is_a_molecular_function": [],
    # EC number edges (functional_annotation_adapter)
    "Ec_number_is_a_ec_number": [],
    "Gene_catalyzes_ec_number": [],
}

SAMPLE_SIZE = 5


def run_query(session, cypher, **params):
    return session.run(cypher, **params).data()


def sample_nodes(session):
    """Sample nodes: anchors + deterministic sample per label."""
    nodes = []

    for label, props in NODE_PROPERTIES.items():
        # Anchors first
        anchor_ids = ANCHOR_NODES.get(label, [])
        for nid in anchor_ids:
            result = run_query(
                session,
                f"MATCH (n:{label}) WHERE n.id = $nid RETURN n",
                nid=nid,
            )
            if result:
                node_props = {k: result[0]["n"].get(k) for k in props if result[0]["n"].get(k) is not None}
                nodes.append({
                    "id": nid,
                    "label": label,
                    "properties": node_props,
                    "anchor": True,
                })

        # Deterministic sample (ORDER BY n.id for reproducibility)
        result = run_query(
            session,
            f"MATCH (n:{label}) RETURN n ORDER BY n.id LIMIT $limit",
            limit=SAMPLE_SIZE,
        )
        for row in result:
            nid = row["n"].get("id")
            if nid and not any(n["id"] == nid for n in nodes):
                node_props = {k: row["n"].get(k) for k in props if row["n"].get(k) is not None}
                nodes.append({
                    "id": nid,
                    "label": label,
                    "properties": node_props,
                    "anchor": False,
                })

    return nodes


def sample_edges(session):
    """Sample edges: deterministic sample per relationship type."""
    edges = []

    for rel_type, props in EDGE_PROPERTIES.items():
        result = run_query(
            session,
            f"MATCH (src)-[r:`{rel_type}`]->(tgt) "
            f"RETURN src.id AS source, tgt.id AS target, properties(r) AS rprops "
            f"ORDER BY src.id, tgt.id LIMIT $limit",
            limit=SAMPLE_SIZE,
        )
        for row in result:
            edge_props = {}
            for k in props:
                val = row["rprops"].get(k)
                if val is not None:
                    # Convert neo4j types to JSON-serializable
                    if isinstance(val, float) and (val != val):  # NaN check
                        continue
                    edge_props[k] = val
            edges.append({
                "type": rel_type,
                "source": row["source"],
                "target": row["target"],
                "properties": edge_props,
            })

    return edges


def generate_snapshot(neo4j_url):
    driver = GraphDatabase.driver(neo4j_url, auth=None)
    driver.verify_connectivity()

    with driver.session() as session:
        nodes = sample_nodes(session)
        edges = sample_edges(session)

    driver.close()

    snapshot = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "description": "Regression snapshot of KG nodes/edges. Regenerate with: uv run python tests/kg_validity/generate_snapshot.py",
        "nodes": nodes,
        "edges": edges,
    }

    SNAPSHOT_PATH.write_text(json.dumps(snapshot, indent=2, default=str) + "\n")
    print(f"Snapshot written to {SNAPSHOT_PATH}")
    print(f"  Nodes: {len(nodes)}")
    print(f"  Edges: {len(edges)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate KG snapshot fixture")
    parser.add_argument(
        "--neo4j-url",
        default="bolt://localhost:7687",
        help="Bolt URL for the Neo4j instance (default: bolt://localhost:7687)",
    )
    args = parser.parse_args()
    generate_snapshot(args.neo4j_url)
