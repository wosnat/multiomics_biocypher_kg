import os, shutil

from biocypher import BioCypher, FileDownload
from multiomics_kg.adapters.ec_adapter import EC
from multiomics_kg.adapters.omics_adapter import MultiOMICSAdapter
from multiomics_kg.adapters.uniprot_adapter import (
    Uniprot,
    MultiUniprot,
    UniprotNodeType,
    UniprotNodeField,
    UniprotEdgeType,
    UniprotIDField,
)
from multiomics_kg.adapters.go_adapter import (
    GO
)

from multiomics_kg.adapters.cyanorak_ncbi_adapter import MultiCyanorakNcbi


def main():
    # Whether to cache data by pypath for future usage
    CACHE = True

    # Flag for exporting node and edge files as csv format
    export_as_csv = True

    # Flag for test mode
    TEST_MODE = False

    # remove these for quick testing. Don't forget to set back to True
    download_GO_data = False
    download_EC_data = False


    # dirs
    output_dir_path = "./biocypher-log/example_knowledge_graph/"
    os.makedirs(output_dir_path, exist_ok=True)

    # Instantiate the BioCypher interface
    # You can use `config/biocypher_config.yaml` to configure the framework or
    # supply settings via parameters below
    bc = BioCypher()


    # uniprot configuration
    uniprot_node_types = [
        UniprotNodeType.PROTEIN,
        #UniprotNodeType.GENE,
        # Organism nodes are created by CyanorakNcbi adapter (single source of truth)
        #UniprotNodeType.ORGANISM,
    ]

    uniprot_node_fields = [
        UniprotNodeField.PRIMARY_GENE_NAME,
        UniprotNodeField.LENGTH,
        UniprotNodeField.MASS,
        UniprotNodeField.ORGANISM,
        UniprotNodeField.ORGANISM_ID,
        UniprotNodeField.PROTEIN_NAMES,
        UniprotNodeField.PROTEIN_GENE_NAMES,
        UniprotNodeField.KEGG_IDS,
        UniprotNodeField.PROTEOME,
        #UniprotNodeField.SEQUENCE, # remove for now to reduce size
        UniprotNodeField.SUBCELLULAR_LOCATION,
        UniprotNodeField.EC,
        #UniprotNodeField.PROTT5_EMBEDDING,
        #UniprotNodeField.ESM2_EMBEDDING,
        #UniprotNodeField.NT_EMBEDDING,
        UniprotNodeField.GENE_ORDERED_LOCUS,
        UniprotNodeField.cc_catalytic_activity,
        UniprotNodeField.cc_cofactor,
        UniprotNodeField.cc_function,
        UniprotNodeField.cc_pathway,
        UniprotNodeField.annotation_score,
        UniprotNodeField.cc_caution,
        UniprotNodeField.keywordid,
        UniprotNodeField.keyword,
        UniprotNodeField.reviewed,
        UniprotNodeField.cc_interaction,
        UniprotNodeField.CELLULAR_COMPONENT,
        UniprotNodeField.BIOLOGICAL_PROCESS,
        UniprotNodeField.MOLECULAR_FUNCTION,
        UniprotNodeField.ft_transmem,
        UniprotNodeField.ft_signal,
        UniprotNodeField.cc_domain,
        UniprotNodeField.ft_motif,
        UniprotNodeField.protein_families,
        UniprotNodeField.xref_refseq,
        UniprotNodeField.xref_string,
        UniprotNodeField.xref_eggnog,
        UniprotNodeField.xref_pfam,

    ]





    uniprot_edge_types = [
        UniprotEdgeType.PROTEIN_TO_ORGANISM,
        UniprotEdgeType.PROTEIN_TO_EC,
        #UniprotEdgeType.GENE_TO_PROTEIN,
        UniprotEdgeType.PROTEIN_TO_CELLULAR_COMPONENT,
        UniprotEdgeType.PROTEIN_TO_BIOLOGICAL_PROCESS,
        UniprotEdgeType.PROTEIN_TO_MOLECULAR_FUNCTION,
    ]

    uniprot_id_type = [
         UniprotIDField.GENE_ENTREZ_ID,
    ]



    # MultiUniprot adapter - uses same config file as MultiCyanorakNcbi
    uniprot_adapter = MultiUniprot(
        config_list_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        rev=False,  # whether to include unreviewed entries
        node_types=uniprot_node_types,
        node_fields=uniprot_node_fields,
        edge_types=uniprot_edge_types,
        id_fields=uniprot_id_type,
        test_mode=TEST_MODE,
    )

    uniprot_adapter.download_uniprot_data(cache=CACHE, retries=6, debug=True)

    bc.write_nodes(uniprot_adapter.get_nodes())
    bc.write_edges(uniprot_adapter.get_edges())

    ncbi_cyanorak_adapter = MultiCyanorakNcbi(
        config_list_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        treatment_organisms_file='data/Prochlorococcus/treatment_organisms.csv',
        test_mode=TEST_MODE,
    )
    ncbi_cyanorak_adapter.download_data(cache=CACHE)
    bc.write_nodes(ncbi_cyanorak_adapter.get_nodes())
    #bc.write_edges(ncbi_cyanorak_adapter.get_edges())
    # if export_as_csv:
    #     ncbi_cyanorak_adapter.export_as_csv(path=output_dir_path)

    # omics data
    omics_adapter = MultiOMICSAdapter(
        config_list_file='data/Prochlorococcus/papers_and_supp/paperconfig_files.txt',
        test_mode=TEST_MODE,
    )
    omics_adapter.download_data(cache=CACHE)
    bc.write_nodes(omics_adapter.get_nodes())
    bc.write_edges(omics_adapter.get_edges())

    # gene ontology
    go_adapter = GO(
        test_mode=TEST_MODE
    )

    if download_GO_data:
        go_adapter.download_go_data(cache=CACHE)
        bc.write_nodes(go_adapter.get_go_nodes())
        bc.write_edges(go_adapter.get_go_edges())
        if export_as_csv:
            go_adapter.export_as_csv(path=output_dir_path)


    if download_EC_data:
        # enzyme commission data
        ec_adapter = EC(
            export_csv=export_as_csv,
            output_dir=output_dir_path,
            test_mode=TEST_MODE
        )
        ec_adapter.download_ec_data(cache=CACHE)
        bc.write_nodes(ec_adapter.get_nodes())
        bc.write_edges(ec_adapter.get_edges())


    # Write import call and other post-processing
    bc.write_schema_info(as_node=True)

    bc.write_import_call()

    # bc.summary()


if __name__ == "__main__":
    main()



