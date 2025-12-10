from biocypher import BioCypher, FileDownload
from template_package.adapters.example_adapter import (
    ExampleAdapter,
    ExampleAdapterNodeType,
    ExampleAdapterEdgeType,
    ExampleAdapterProteinField,
    ExampleAdapterDiseaseField,
)
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
TEST_MODE = True

# dirs
output_dir_path = "."


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
bc.summary()



# # Download and cache resources (change the directory in the options if needed)
# urls = "https://file-examples.com/wp-content/storage/2017/02/file_example_CSV_5000.csv"
# resource = FileDownload(
#     name="Example resource",  # Name of the resource
#     url_s=urls,  # URL to the resource(s)
#     lifetime=7,  # seven days cache lifetime
# )
# paths = bc.download(resource)  # Downloads to '.cache' by default
# print(paths)
# You can use the list of paths returned to read the resource into your adapter

# Choose node types to include in the knowledge graph.
# These are defined in the adapter (`adapter.py`).
# node_types = [
#     ExampleAdapterNodeType.PROTEIN,
#     ExampleAdapterNodeType.DISEASE,
# ]

# # Choose protein adapter fields to include in the knowledge graph.
# # These are defined in the adapter (`adapter.py`).
# node_fields = [
#     # Proteins
#     ExampleAdapterProteinField.ID,
#     ExampleAdapterProteinField.SEQUENCE,
#     ExampleAdapterProteinField.DESCRIPTION,
#     ExampleAdapterProteinField.TAXON,
#     # Diseases
#     ExampleAdapterDiseaseField.ID,
#     ExampleAdapterDiseaseField.NAME,
#     ExampleAdapterDiseaseField.DESCRIPTION,
# ]

# edge_types = [
#     ExampleAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION,
#     ExampleAdapterEdgeType.PROTEIN_DISEASE_ASSOCIATION,
# ]

# # Create a protein adapter instance
# adapter = ExampleAdapter(
#     node_types=node_types,
#     node_fields=node_fields,
#     edge_types=edge_types,
#     # we can leave edge fields empty, defaulting to all fields in the adapter
# )


# # Create a knowledge graph from the adapter
# bc.write_nodes(adapter.get_nodes())
# bc.write_edges(adapter.get_edges())

# # Write admin import statement
# bc.write_import_call()

# # Print summary
# bc.summary()
