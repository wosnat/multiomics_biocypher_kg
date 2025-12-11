import os, shutil
from biocypher import BioCypher, FileDownload
from template_package.adapters.uniprot_adapter import (
    Uniprot,
    UniprotNodeType,
    UniprotNodeField,
    UniprotEdgeType,
    UniprotIDField,
)



# Whether to cache data by pypath for future usage
CACHE = True

# Flag for exporting node and edge files as csv format
export_as_csv = True

# Flag for test mode
TEST_MODE = False

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
    UniprotNodeType.GENE,
    UniprotNodeType.ORGANISM,
]

uniprot_node_fields = [
    UniprotNodeField.PRIMARY_GENE_NAME,
    UniprotNodeField.LENGTH,
    UniprotNodeField.MASS,
    UniprotNodeField.ORGANISM,
    UniprotNodeField.ORGANISM_ID,
    UniprotNodeField.PROTEIN_NAMES,
    UniprotNodeField.PROTEIN_GENE_NAMES,
    UniprotNodeField.ENSEMBL_TRANSCRIPT_IDS,
    UniprotNodeField.ENSEMBL_GENE_IDS,
    UniprotNodeField.ENTREZ_GENE_IDS,
    UniprotNodeField.KEGG_IDS,
    UniprotNodeField.PROTEOME,
    UniprotNodeField.SEQUENCE,
    #UniprotNodeField.PROTT5_EMBEDDING,
    #UniprotNodeField.ESM2_EMBEDDING,
    #UniprotNodeField.NT_EMBEDDING,
]

uniprot_edge_types = [
     UniprotEdgeType.PROTEIN_TO_ORGANISM,
     UniprotEdgeType.GENE_TO_PROTEIN,
]

uniprot_id_type = [
     UniprotIDField.GENE_ENTREZ_ID,
]


uniprot_adapter = Uniprot(
        organism="59919", # MED4
        node_types=uniprot_node_types,
        node_fields=uniprot_node_fields,
        edge_types=uniprot_edge_types,
        id_fields=uniprot_id_type,
        test_mode=TEST_MODE,
    )

uniprot_adapter.download_uniprot_data(cache=CACHE, retries=6)

uniprot_nodes = uniprot_adapter.get_nodes()
uniprot_edges = uniprot_adapter.get_edges()


bc.write_nodes(uniprot_nodes)
bc.write_edges(uniprot_edges)


if export_as_csv:
    uniprot_adapter.export_data_to_csv(path=output_dir_path,
                                    node_data=uniprot_nodes,
                                    edge_data=uniprot_edges)


# Write import call and other post-processing
bc.write_import_call()
#bc.summary()



