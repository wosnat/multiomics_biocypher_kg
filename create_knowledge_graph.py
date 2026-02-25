import os

from biocypher import BioCypher
from multiomics_kg.adapters.ec_adapter import EC
from multiomics_kg.adapters.omics_adapter import MultiOMICSAdapter
from multiomics_kg.adapters.uniprot_adapter import MultiUniprot
from multiomics_kg.adapters.go_adapter import GO

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
