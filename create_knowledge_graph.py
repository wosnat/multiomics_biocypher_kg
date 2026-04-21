import argparse
import os
from pathlib import Path

from biocypher import BioCypher
from multiomics_kg.adapters.omics_adapter import MultiOMICSAdapter
from multiomics_kg.adapters.cluster_adapter import MultiClusterAdapter
from multiomics_kg.adapters.observations_adapter import MultiObservationsAdapter
from multiomics_kg.adapters.uniprot_adapter import MultiUniprot
from multiomics_kg.adapters.go_adapter import GO

from multiomics_kg.adapters.cyanorak_ncbi_adapter import MultiCyanorakNcbi
from multiomics_kg.adapters.ortholog_group_adapter import MultiOrthologGroupAdapter
from multiomics_kg.adapters.functional_annotation_adapter import (
    MultiGoAnnotationAdapter,
    MultiEcAnnotationAdapter,
    MultiKeggAnnotationAdapter,
    MultiCogRoleAnnotationAdapter,
    MultiPfamAnnotationAdapter,
)
from multiomics_kg.adapters.brite_adapter import MultiBriteAdapter


def parse_args():
    parser = argparse.ArgumentParser(description="Build the multiomics BioCypher knowledge graph.")
    parser.add_argument("--test", action="store_true", help="Test mode: stop each adapter after 100 items.")
    parser.add_argument("--go", action="store_true", help="Download and write full GO ontology nodes/edges.")
    parser.add_argument("--no-cache", action="store_true", help="Re-fetch data instead of using cached files.")
    parser.add_argument("--output-dir", default="./biocypher-log/example_knowledge_graph/",
                        help="Output directory for CSV exports (default: ./biocypher-log/example_knowledge_graph/).")
    return parser.parse_args()


def main():
    args = parse_args()

    CACHE = not args.no_cache
    export_as_csv = True
    TEST_MODE = args.test
    download_GO_data = args.go
    output_dir_path = args.output_dir

    os.makedirs(output_dir_path, exist_ok=True)

    # Instantiate the BioCypher interface
    # You can use `config/biocypher_config.yaml` to configure the framework or
    # supply settings via parameters below
    bc = BioCypher()


    # CyanorakNcbi adapter MUST run before UniProt: it creates gene_mapping.csv
    # files in each data_dir, which UniProt uses for GENE_TO_PROTEIN edges.
    ncbi_cyanorak_adapter = MultiCyanorakNcbi(
        config_list_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        treatment_organisms_file='data/Prochlorococcus/treatment_organisms.csv',
        test_mode=TEST_MODE,
    )
    ncbi_cyanorak_adapter.download_data(cache=CACHE)
    bc.write_nodes(ncbi_cyanorak_adapter.get_nodes())
    bc.write_edges(ncbi_cyanorak_adapter.get_edges())

    # OrthologGroup adapter: reads ortholog_groups from gene_annotations_merged.json
    og_adapter = MultiOrthologGroupAdapter(
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        test_mode=TEST_MODE,
    )
    bc.write_nodes(og_adapter.get_nodes())
    bc.write_edges(og_adapter.get_edges())

    # MultiUniprot adapter reads pre-built protein_annotations.json files.
    # Requires: prepare_data.sh steps 0 + 2 run beforehand.
    # GENE_TO_PROTEIN edges depend on gene_mapping.csv from CyanorakNcbi above.
    uniprot_adapter = MultiUniprot(
        config_list_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        test_mode=TEST_MODE,
    )
    uniprot_adapter.download_data(cache=CACHE)
    bc.write_nodes(uniprot_adapter.get_nodes())
    bc.write_edges(uniprot_adapter.get_edges())

    # omics data
    omics_adapter = MultiOMICSAdapter(
        config_list_file=[
            'data/Prochlorococcus/papers_and_supp/paperconfig_files.txt',
            'data/Synechococcus/papers_and_supp/paperconfig_files.txt',
        ],
        test_mode=TEST_MODE,
    )
    omics_adapter.download_data(cache=CACHE)
    bc.write_nodes(omics_adapter.get_nodes())
    bc.write_edges(omics_adapter.get_edges())

    # Gene cluster data (co-expression clusters from publications)
    cluster_adapter = MultiClusterAdapter(
        config_list_file=[
            'data/Prochlorococcus/papers_and_supp/paperconfig_files.txt',
            'data/Synechococcus/papers_and_supp/paperconfig_files.txt',
        ],
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        test_mode=TEST_MODE,
    )
    cluster_nodes = cluster_adapter.get_nodes()
    cluster_edges = cluster_adapter.get_edges()
    if cluster_nodes:
        bc.write_nodes(cluster_nodes)
    if cluster_edges:
        bc.write_edges(cluster_edges)

    # Derived-metric observations (non-DE evidence; Biller 2018 boolean + categorical,
    # future zinser 2009 / Waldbauer 2012 numeric). One DerivedMetric node per
    # (Experiment × metric_type) + 3 binding edges + 1 measurement edge per metric.
    observations_adapter = MultiObservationsAdapter(
        config_list_file=[
            'data/Prochlorococcus/papers_and_supp/paperconfig_files.txt',
            'data/Synechococcus/papers_and_supp/paperconfig_files.txt',
            'tests/fixtures/non_de/paperconfig_files.txt',  # TEMP Plan 3 Task 10 — unwire before Task 15
        ],
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        test_mode=TEST_MODE,
    )
    obs_nodes = observations_adapter.get_nodes()
    obs_edges = observations_adapter.get_edges()
    if obs_nodes:
        bc.write_nodes(obs_nodes)
    if obs_edges:
        bc.write_edges(obs_edges)

    # Gene → GO annotation edges + GO hierarchy subset (lightweight, always runs)
    go_anno_adapter = MultiGoAnnotationAdapter(
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        cache_root=Path("cache/data"),
        test_mode=TEST_MODE,
        cache=CACHE,
    )
    bc.write_nodes(go_anno_adapter.get_nodes())
    bc.write_edges(go_anno_adapter.get_edges())

    # EC number nodes + hierarchy edges + gene→EC edges (always runs, cached)
    ec_anno_adapter = MultiEcAnnotationAdapter(
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        cache_dir=Path("cache/data/ec"),
        test_mode=TEST_MODE,
        cache=CACHE,
    )
    bc.write_nodes(ec_anno_adapter.get_nodes())
    bc.write_edges(ec_anno_adapter.get_edges())

    # KEGG 4-level hierarchy: KO → Pathway → Subcategory → Category (always runs, cached)
    kegg_anno_adapter = MultiKeggAnnotationAdapter(
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        cache_root=Path("cache/data"),
        test_mode=TEST_MODE,
        cache=CACHE,
    )
    bc.write_nodes(kegg_anno_adapter.get_nodes())
    bc.write_edges(kegg_anno_adapter.get_edges())

    # KEGG BRITE functional hierarchies: 12 trees → BriteCategory nodes + hierarchy/KO edges.
    # Pruned to the subhierarchy reachable from KeggTerm nodes (gene-referenced KOs).
    brite_adapter = MultiBriteAdapter(
        cache_root=Path("cache/data"),
        known_ko_ids=kegg_anno_adapter.all_ko_ids(),
        test_mode=TEST_MODE,
        cache=CACHE,
    )
    brite_adapter.download_data(cache=CACHE)
    bc.write_nodes(brite_adapter.get_nodes())
    bc.write_edges(brite_adapter.get_edges())

    # COG functional categories + Cyanorak roles + tIGR roles
    cog_role_adapter = MultiCogRoleAnnotationAdapter(
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        role_tree_file=Path("data/cyanorak_roles.csv"),
        test_mode=TEST_MODE,
    )
    bc.write_nodes(cog_role_adapter.get_nodes())
    bc.write_edges(cog_role_adapter.get_edges())

    # Pfam domain families + PfamClan superfamilies + gene→Pfam edges (always runs, cached)
    pfam_adapter = MultiPfamAnnotationAdapter(
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        cache_root=Path("cache/data"),
        test_mode=TEST_MODE,
        cache=CACHE,
    )
    bc.write_nodes(pfam_adapter.get_nodes())
    bc.write_edges(pfam_adapter.get_edges())

    # Full GO ontology (all 30K nodes + GO-GO hierarchy) — optional, slow.
    # NOTE: do not run with --go simultaneously; GO node IDs would conflict.

    if download_GO_data:
        go_adapter = GO(
            test_mode=TEST_MODE,
            cache_root=Path("cache/data"),
        )
        go_adapter.download_go_data(cache=CACHE)
        bc.write_nodes(go_adapter.get_go_nodes())
        bc.write_edges(go_adapter.get_go_edges())
        if export_as_csv:
            go_adapter.export_as_csv(path=output_dir_path)


    # Write import call and other post-processing
    bc.write_schema_info(as_node=True)

    bc.write_import_call()

    # bc.summary()


if __name__ == "__main__":
    main()
