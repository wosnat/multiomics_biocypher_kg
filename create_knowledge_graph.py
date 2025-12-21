import os, shutil

from biocypher import BioCypher, FileDownload
from template_package.adapters.uniprot_adapter import (
    Uniprot,
    UniprotNodeType,
    UniprotNodeField,
    UniprotEdgeType,
    UniprotIDField,
)
from template_package.adapters.go_adapter import (
    GO
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
    #UniprotNodeType.GENE,
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
    UniprotNodeField.go,
    UniprotNodeField.go_id,
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
     #UniprotEdgeType.GENE_TO_PROTEIN,
]

uniprot_id_type = [
     UniprotIDField.GENE_ENTREZ_ID,
]



# organism to use: MED4 (59919)
organism="59919" # MED4

uniprot_adapter = Uniprot(
        organism=organism, # MED4
        rev= False, # whether to include unreviewed entries
        node_types=uniprot_node_types,
        node_fields=uniprot_node_fields,
        edge_types=uniprot_edge_types,
        id_fields=uniprot_id_type,
        test_mode=TEST_MODE,
    )

uniprot_adapter.download_uniprot_data(cache=CACHE, retries=6, debug=True)

uniprot_nodes = uniprot_adapter.get_nodes()
uniprot_edges = uniprot_adapter.get_edges()



bc.write_nodes(uniprot_nodes)
bc.write_edges(uniprot_edges)


if export_as_csv:
    uniprot_adapter.export_data_to_csv(path=output_dir_path,
                                    node_data=uniprot_nodes,
                                    edge_data=uniprot_edges)


# gene ontology

go_adapter = GO(
    organism=organism, 
    test_mode=TEST_MODE
)
go_adapter.download_go_data(cache=CACHE)
bc.write_nodes(go_adapter.get_go_nodes())
bc.write_edges(go_adapter.get_go_edges())
if export_as_csv:
    go_adapter.export_as_csv(path=output_dir_path)


# Write import call and other post-processing
bc.write_import_call()
# bc.summary()



