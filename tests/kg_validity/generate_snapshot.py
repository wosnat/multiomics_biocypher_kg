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
    "Publication": [
        "doi:10.1038/msb4100087",       # Tolonen 2006 — nitrogen stress, MED4+MIT9313
        "doi:10.1038/ismej.2016.70",     # Biller 2016 — coculture, MED4+MIT9313+HOT1A3
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
    # KEGG anchor nodes
    "KeggTerm": [
        "kegg.orthology:K02338",  # dnaN — DNA pol III beta subunit, common in prokaryotes
        "kegg.category:09100",  # Metabolism — top-level BRITE category
    ],
    # COG anchor nodes — stable hardcoded categories
    "CogFunctionalCategory": [
        "cog.category:J",  # Translation, ribosomal structure and biogenesis
        "cog.category:C",  # Energy production and conversion
    ],
    # CyanorakRole anchor nodes — stable top-level roles
    "CyanorakRole": [
        "cyanorak.role:0",   # Non-coding gene (RNA) — root of ncRNA subtree
        "cyanorak.role:A",   # Amino acid biosynthesis — always present
    ],
    "DerivedMetric": [
        # Biller 2018 retrofitted boolean + categorical DMs — always present post-Plan-3
        "derived_metric:mSystems.00040-18:s4a_natl2a_axenic:periodic_in_axenic_LD",
        "derived_metric:mSystems.00040-18:s5_natl2a_survival:darkness_survival_class",
        # Waldbauer 2012 numeric DM — stable distribution (peak phase 0–24h)
        "derived_metric:journal.pone.0043432:table_s2_waldbauer_diel_metrics:peak_time_protein_h",
    ],
    # Metabolism scaffold (KEGG-primary)
    "Reaction": [
        "kegg.reaction:R11945",   # NADH:ubiquinone oxidoreductase — universal, top gene_count
    ],
    "Metabolite": [
        "kegg.compound:C00002",   # ATP — universal, multi-pathway
        "kegg.compound:C00031",   # D-Glucose
    ],
    # Phase 2 metabolomics
    "MetaboliteAssay": [
        # Kujawinski 2023 MIT9301 intracellular — stable across rebuilds
        "metabolite_assay:msystems.01261-22:metabolites_kegg_export_9301_intracellular:cellular_concentration",
    ],
    # Pfam normalization
    "Pfam": [
        "pfam:PF00005",   # ABC_tran — ubiquitous transporter domain
    ],
    "PfamClan": [
        "pfam.clan:CL0023",   # P-loop_NTPase — large, stable
    ],
    # TCDB transport ontology
    "TcdbFamily": [
        "tcdb:1",            # Channels and Pores — class root
        "tcdb:3.A.1.4.4",    # Polar amino acid ABC transporter — leaf with many genes
        # Tripwires for the acc2tcid hierarchy-seeding fix (May 2026):
        # these were silently dropped before the build seeded from acc2tcid.tsv.
        "tcdb:1.A.11",       # Ammonium Transporter Channel (Amt) Family — 3-part family from families.tsv
        "tcdb:1.A.12.2.2",   # chloride channel — 5-part specificity exclusive to acc2tcid
    ],
    # CAZy carbohydrate-active enzymes
    "CazyFamily": [
        "cazy:GH",   # Glycoside Hydrolases class root
        "cazy:GT2",  # GT2 family — top gene count
    ],
    # KEGG BRITE functional hierarchies
    "BriteCategory": [
        "kegg.brite:ko01000.A1",  # Enzymes / Oxidoreductases — A-level root
    ],
    # Clustering layer
    "ClusteringAnalysis": [
        "clustering_analysis:mSystems.00181-16:bp1_light_clusters",
    ],
    "GeneCluster": [
        "cluster:mSystems.00181-16:bp1_light_clusters:C",  # 929 members
    ],
    # Provenance nodes (4 total, all stable)
    "DataSource": [
        "data_source:ncbi",
        "data_source:cyanorak",
        "data_source:eggnog",
        "data_source:uniprot",
    ],
}

# Key properties to capture per node type
NODE_PROPERTIES = {
    "Gene": ["locus_tag", "product", "gene_name", "gene_name_synonyms", "strand"],
    "Protein": ["protein_synonyms", "gene_names"],
    "OrganismTaxon": ["organism_name", "strain_name", "genus", "clade", "ncbi_taxon_id"],
    "Publication": ["doi", "pmid", "title", "experiment_count", "treatment_types", "omics_types", "organisms"],
    "Experiment": ["name", "organism_name", "treatment_type", "treatment", "control", "omics_type", "is_time_course"],
    "OrthologGroup": ["name", "source", "taxonomic_level", "taxon_id"],
    # GO term node types (from functional_annotation_adapter.py)
    "BiologicalProcess": ["name"],
    "CellularComponent": ["name"],
    "MolecularFunction": ["name"],
    # EC number nodes (from ec_adapter / functional_annotation_adapter)
    "EcNumber": ["name"],
    # KEGG nodes (from MultiKeggAnnotationAdapter)
    "KeggTerm": ["name", "level"],
    # COG / role nodes (from MultiCogRoleAnnotationAdapter)
    "CogFunctionalCategory": ["code", "name"],
    "CyanorakRole": ["code", "description"],
    "TigrRole": ["code", "description"],
    "DerivedMetric": [
        "name", "metric_type", "value_kind", "experiment_id",
        "organism_name", "compartment", "omics_type", "total_gene_count",
        "rankable", "has_p_value",
        # numeric distribution stats (numeric DMs only)
        "value_min", "value_max", "value_q1", "value_median", "value_q3",
        # boolean flag counts (boolean DMs only)
        "flag_true_count", "flag_false_count",
        # categorical distribution (categorical DMs only)
        "category_labels", "category_counts",
    ],
    # Metabolism scaffold (KEGG-primary IDs)
    "Reaction": [
        "name", "kegg_reaction_id", "reaction_class", "mass_balance",
    ],
    "Metabolite": [
        "name", "kegg_compound_id", "formula", "mass",
    ],
    # Phase 2 metabolomics
    "MetaboliteAssay": [
        "name", "metric_type", "value_kind", "experiment_id",
        "organism_name", "compartment", "omics_type",
        "rankable", "aggregation_method",
    ],
    # Pfam
    "Pfam": ["name", "short_name"],
    "PfamClan": ["name"],
    # TCDB / CAZy / BRITE — share the same level/level_kind unified-ontology shape
    "TcdbFamily": ["name", "tcdb_id", "level", "level_kind"],
    "CazyFamily": ["name", "cazy_id", "level", "level_kind"],
    "BriteCategory": ["name", "tree", "tree_code", "level", "level_kind"],
    # Clustering layer
    "ClusteringAnalysis": [
        "name", "organism_name", "cluster_method", "cluster_type",
        "cluster_count", "omics_type",
    ],
    "GeneCluster": ["name", "organism_name", "member_count"],
    # Provenance — 4 stable nodes
    "DataSource": ["scope", "provenance", "description"],
}

# Key properties to capture per edge type
EDGE_PROPERTIES = {
    "Gene_belongs_to_organism": [],
    "Protein_belongs_to_organism": [],
    "Gene_encodes_protein": [],
    "Gene_in_ortholog_group": [],
    "Changes_expression_of": [
        "expression_direction", "log2_fold_change", "adjusted_p_value",
        "time_point", "time_point_order", "time_point_hours",
    ],
    "Has_experiment": [],
    "Tests_coculture_with": [],
    # Gene → GO edges (functional_annotation_adapter.py; label_as_edge values)
    "Gene_involved_in_biological_process": [],
    "Gene_located_in_cellular_component": [],
    "Gene_enables_molecular_function": [],
    # GO-GO hierarchy edges (label_as_edge in schema; BioCypher capitalizes first letter)
    "Biological_process_is_a_biological_process": [],
    "Biological_process_part_of_biological_process": [],
    "Biological_process_negatively_regulates_biological_process": [],
    "Biological_process_positively_regulates_biological_process": [],
    "Cellular_component_is_a_cellular_component": [],
    "Cellular_component_part_of_cellular_component": [],
    "Molecular_function_is_a_molecular_function": [],
    "Molecular_function_part_of_molecular_function": [],
    # EC number edges (functional_annotation_adapter)
    "Ec_number_is_a_ec_number": [],
    "Gene_catalyzes_ec_number": [],
    # KEGG edges (MultiKeggAnnotationAdapter)
    "Gene_has_kegg_ko": [],
    "Kegg_term_is_a_kegg_term": [],
    # COG / role edges (MultiCogRoleAnnotationAdapter)
    "Gene_in_cog_category": [],
    "Gene_has_cyanorak_role": [],
    "Gene_has_tigr_role": [],
    "Cyanorak_role_is_a_cyanorak_role": [],
    # DerivedMetric binding edges (Plan 3)
    "PublicationHasDerivedMetric": [],
    "ExperimentHasDerivedMetric": [],
    "DerivedMetricBelongsToOrganism": [],
    # DerivedMetric measurement edges (Plan 3)
    "Derived_metric_flags_gene": ["metric_type", "value"],
    "Derived_metric_classifies_gene": ["metric_type", "value"],
    "Derived_metric_quantifies_gene": [
        "metric_type", "value", "adjusted_p_value",
        "rank_by_metric", "metric_percentile", "metric_bucket", "significant",
    ],
    # Metabolism scaffold edges
    "Gene_catalyzes_reaction": [],
    "Reaction_has_metabolite": [],
    "Reaction_in_kegg_pathway": [],
    "Metabolite_in_pathway": [],
    "Organism_has_metabolite": [],
    # Phase 2 metabolomics — measurement + binding edges
    "Assay_quantifies_metabolite": [
        "value", "n_replicates", "metric_bucket",
    ],
    "Assay_flags_metabolite": [
        "flag_value", "n_replicates", "n_positive",
    ],
    "PublicationHasMetaboliteAssay": [],
    "ExperimentHasMetaboliteAssay": [],
    "MetaboliteAssayBelongsToOrganism": [],
    # Pfam edges
    "Gene_has_pfam": [],
    "Pfam_in_pfam_clan": [],
    # TCDB edges (gene-attach is at the most specific eggNOG-annotated level;
    # transports rolled up from tc_specificity leaves)
    "Gene_has_tcdb_family": [],
    "Tcdb_family_is_a_tcdb_family": [],
    "Tcdb_family_transports_metabolite": [],
    # CAZy edges
    "Gene_has_cazy_family": [],
    "Cazy_family_is_a_cazy_family": [],
    # KEGG BRITE edges
    "Brite_category_is_a_brite_category": [],
    "Kegg_term_in_brite_category": [],
    # OrthologGroup rollup edges (majority-vote from member genes)
    "Og_has_cyanorak_role": [],
    "Og_in_cog_category": [],
    # Clustering layer edges
    "PublicationHasClusteringAnalysis": [],
    "ClusteringAnalysisHasGeneCluster": [],
    "ClusteringanalysisBelongsToOrganism": [],
    "ExperimentHasClusteringAnalysis": [],
    "Gene_in_gene_cluster": [],
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


def _collect_edge_rows(result, props):
    """Convert query rows into edge dicts with the requested properties."""
    edges = []
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
            "source": row["source"],
            "target": row["target"],
            "properties": edge_props,
        })
    return edges


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
        for d in _collect_edge_rows(result, props):
            edges.append({"type": rel_type, **d})

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

    SNAPSHOT_PATH.write_text(json.dumps(snapshot, indent=2, default=str, sort_keys=True) + "\n")
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
