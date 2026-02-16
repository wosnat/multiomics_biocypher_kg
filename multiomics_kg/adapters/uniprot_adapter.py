from time import time
import collections
import csv
from typing import Optional, Union, Literal
from collections.abc import Generator
from enum import Enum, EnumMeta, auto
from functools import lru_cache
import pandas as pd
import numpy as np

import os
import requests
import h5py

from tqdm import tqdm  # progress bar
from pypath.share import curl, settings
from pypath.utils import mapping
from pypath.inputs import uniprot
from biocypher._logger import logger
from contextlib import ExitStack
from bioregistry import normalize_curie

from pydantic import BaseModel, DirectoryPath, FilePath, HttpUrl, validate_call

logger.debug(f"Loading module {__name__}.")


# Uniprot is a hub of many biological databases, the result of this adapter is:
# - protein nodes with properties from uniprot
# - edge nodes linking to various databases:
#   - organisms (NCBI Taxonomy)
#   - KEGG pathways (KEGG)
#   - EC numbers (ExPASy)
#   - GO terms (Gene Ontology) - handled in GO adapter

# edges that we plan to add in the second phase:
#   - genes (NCBI GeneID)
#   - Proteomes (Uniprot Proteomes)
#   - RefSeq (NCBI RefSeq)
#   - STRING (STRING)
#   - EggNOG (EggNOG)
#   - KO (KEGG Orthology)
#   - PFAM (PFAM)









class UniprotEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class UniprotNodeType(Enum, metaclass=UniprotEnumMeta):
    """
    Node types of the UniProt API represented in this adapter.
    """

    PROTEIN = auto()
    GENE = auto()
    ORGANISM = auto()
    CELLULAR_COMPARTMENT = auto()


class UniprotNodeField(Enum, metaclass=UniprotEnumMeta):
    """
    Fields of nodes the UniProt API represented in this adapter. Overview of
    uniprot fields: https://www.uniprot.org/help/return_fields
    """

    # core attributes
    LENGTH = "length"
    SUBCELLULAR_LOCATION = "subcellular_location"
    MASS = "mass"
    ORGANISM = "organism_name"
    ORGANISM_ID = "organism_id"
    PROTEIN_NAMES = "protein_name"
    EC = "ec"
    PROTEIN_GENE_NAMES = "gene_names"
    PRIMARY_GENE_NAME = "gene_primary"
    SEQUENCE = "sequence"

    # added to check 
    GENE_ORDERED_LOCUS = "gene_oln"
    cc_catalytic_activity = "cc_catalytic_activity"
    cc_cofactor = "cc_cofactor"
    cc_function = "cc_function"
    cc_pathway = "cc_pathway"
    annotation_score = "annotation_score"
    cc_caution = "cc_caution"
    keywordid = "keywordid"
    keyword = "keyword"
    reviewed = "reviewed"
    cc_interaction = "cc_interaction"
    CELLULAR_COMPONENT = "go_c"
    BIOLOGICAL_PROCESS = "go_p"
    MOLECULAR_FUNCTION = "go_f"

    # post processing of go fields to extract go ids
    CELLULAR_COMPONENT_ID = "go_c_id"
    BIOLOGICAL_PROCESS_ID = "go_p_id"
    MOLECULAR_FUNCTION_ID = "go_f_id"

    #go_id = "go_id"
    ft_transmem = "ft_transmem"
    ft_signal = "ft_signal"
    cc_domain = "cc_domain"
    #ft_domain = "ft_domain"
    ft_motif = "ft_motif"
    protein_families = "protein_families"
    #ft_region = "ft_region"
    xref_refseq = "xref_refseq"
    xref_string = "xref_string"
    xref_eggnog = "xref_eggnog"
    xref_ko = "xref_ko"
    xref_pfam   = "xref_pfam"






    # xref attributes
    PROTEOME = "xref_proteomes"
    ENTREZ_GENE_IDS = "xref_geneid"
    KEGG_IDS = "xref_kegg"


    # not from uniprot REST
    # we provide these by downloading the ProtT5 embeddings from uniprot
    PROTT5_EMBEDDING = "prott5_embedding"

    # not from uniprot REST
    # we provide these by getting embeddings from ESM2 650M model
    ESM2_EMBEDDING = "esm2_embedding"

    # not from uniprot REST
    # we provide these by getting embeddings from Nucletide tranformer 2 model
    NT_EMBEDDING = "nt_embedding"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None
    @classmethod
    def get_protein_properties(cls):
        return [
            cls.LENGTH.value,
            cls.MASS.value,
            cls.PROTEIN_NAMES.value,
            cls.PROTEOME.value,
            cls.EC.value,
            cls.ORGANISM.value,
            cls.ORGANISM_ID.value,
            cls.SEQUENCE.value,
            cls.PROTT5_EMBEDDING.value,
            cls.ESM2_EMBEDDING.value,
            cls.PROTEIN_GENE_NAMES.value,
            cls.ENTREZ_GENE_IDS.value,
            cls.KEGG_IDS.value,
            cls.PRIMARY_GENE_NAME.value,
            cls.NT_EMBEDDING.value,
            cls.SUBCELLULAR_LOCATION.value,
            # added to test 
            cls.GENE_ORDERED_LOCUS.value,   
            cls.cc_catalytic_activity.value,
            cls.cc_cofactor.value,
            cls.cc_function.value,
            cls.cc_pathway.value,
            cls.annotation_score.value,
            cls.cc_caution.value,
            cls.keywordid.value,
            cls.keyword.value,
            cls.reviewed.value,
            cls.cc_interaction.value,
            cls.CELLULAR_COMPONENT.value,
            cls.BIOLOGICAL_PROCESS.value,
            cls.MOLECULAR_FUNCTION.value,
            cls.ft_transmem.value,
            cls.ft_signal.value,
            cls.cc_domain.value,
            cls.ft_motif.value,
            cls.protein_families.value,
            cls.xref_refseq.value,
            cls.xref_string.value,
            cls.xref_eggnog.value,
            cls.xref_ko.value,
            cls.xref_pfam.value,
        ]
    @classmethod
    def get_gene_properties(cls) -> list:
        return [
            cls.PROTEIN_GENE_NAMES.value,
            cls.ENTREZ_GENE_IDS.value,
            cls.KEGG_IDS.value,
            cls.PRIMARY_GENE_NAME.value,
            cls.NT_EMBEDDING.value,
        ]
    
    @classmethod
    def get_organism_properties(cls) -> list:
        return [cls.ORGANISM.value]
    
    @classmethod
    def get_split_fields(cls) -> list:
        return [
            cls.PROTEOME.value,
            cls.PROTEIN_GENE_NAMES.value,
            cls.EC.value,
            cls.ENTREZ_GENE_IDS.value,
            #cls.ENSEMBL_TRANSCRIPT_IDS.value,
            cls.KEGG_IDS.value,
            cls.CELLULAR_COMPONENT.value,
            cls.MOLECULAR_FUNCTION.value,
            cls.BIOLOGICAL_PROCESS.value,
        ]
    
    @classmethod
    def get_nonuniprot_api_fields(cls) -> list:
        return [
            cls.PROTT5_EMBEDDING.value,
            cls.ESM2_EMBEDDING.value,
            cls.NT_EMBEDDING.value,
            cls.CELLULAR_COMPONENT_ID.value,
            cls.MOLECULAR_FUNCTION_ID.value,
            cls.BIOLOGICAL_PROCESS_ID.value,
        ]


class UniprotEdgeType(Enum, metaclass=UniprotEnumMeta):
    """
    Edge types of the UniProt API represented in this adapter.
    """

    PROTEIN_TO_ORGANISM = auto()
    GENE_TO_PROTEIN = auto()
    PROTEIN_TO_EC = auto()
    PROTEIN_TO_CELLULAR_COMPONENT = auto()
    PROTEIN_TO_BIOLOGICAL_PROCESS = auto()
    PROTEIN_TO_MOLECULAR_FUNCTION = auto()



# to add an edge you need: uniprot field name, target edge name, prefix for target database, edge label
class UniprotEdgeModel(BaseModel):
    """
        to add an edge you need: uniprot field name, target edge name, prefix for target database, edge label
    """
    uniprot_field: UniprotNodeField
    edge_type: UniprotEdgeType
    target_prefix: str
    edge_label: str

UNIPROT_EDGE_MODELS = [
    UniprotEdgeModel(
        uniprot_field=UniprotNodeField.ORGANISM_ID,
        edge_type=UniprotEdgeType.PROTEIN_TO_ORGANISM,
        target_prefix="ncbitaxon",
        edge_label="Protein_belongs_to_organism",
    ),
    UniprotEdgeModel(
        uniprot_field=UniprotNodeField.ENTREZ_GENE_IDS,
        edge_type=UniprotEdgeType.GENE_TO_PROTEIN,
        target_prefix="ncbigene",
        edge_label="Gene_encodes_protein",
    ),
    UniprotEdgeModel(
        uniprot_field=UniprotNodeField.EC,
        edge_type=UniprotEdgeType.PROTEIN_TO_EC,
        target_prefix="eccode",
        edge_label="protein_catalyzes_ec_number",
    ),
    UniprotEdgeModel(
        uniprot_field=UniprotNodeField.CELLULAR_COMPONENT_ID,
        edge_type=UniprotEdgeType.PROTEIN_TO_CELLULAR_COMPONENT,
        target_prefix="go",
        edge_label="protein_located_in_cellular_component",
    ),
    UniprotEdgeModel(
        uniprot_field=UniprotNodeField.BIOLOGICAL_PROCESS_ID,
        edge_type=UniprotEdgeType.PROTEIN_TO_BIOLOGICAL_PROCESS,
        target_prefix="go",
        edge_label="protein_involved_in_biological_process",
    ),
    UniprotEdgeModel(
        uniprot_field=UniprotNodeField.MOLECULAR_FUNCTION_ID,
        edge_type=UniprotEdgeType.PROTEIN_TO_MOLECULAR_FUNCTION,
        target_prefix="go",
        edge_label="protein_contributes_to_molecular_function",
    ),
]





class UniprotIDField(Enum, metaclass=UniprotEnumMeta):
    """
    Fields of edges of the UniProt API represented in this adapter. Used to
    assign source and target identifiers in `get_edges()`.
    """

    # default
    PROTEIN_UNIPROT_ACCESSION = auto()
    GENE_ENTREZ_ID = auto()
    ORGANISM_NCBI_TAXONOMY_ID = auto()

    # optional
    GENE_ENSEMBL_GENE_ID = auto()


class UniProtModel(BaseModel):
    organism: Literal["*"] | int | None = "*"
    rev: bool = True
    node_types: Optional[Union[list[UniprotNodeType], None]] = None
    node_fields: Optional[Union[list[UniprotNodeField], None]] = None
    edge_types: Optional[Union[list[UniprotEdgeType], None]] = None
    id_fields: Optional[Union[list[UniprotIDField], None]] = None
    add_prefix: bool = True
    test_mode: bool = False


class Uniprot:
    """
    Class that downloads uniprot data using pypath and reformats it to be ready
    for import into a BioCypher database.

    Args:
        organism: organism code in NCBI taxid format, e.g. "9606" for human.
        rev: if True, it downloads swissprot (i.e., reviewed) entries only.
        node_types: `UniprotNodeType` fields that will be included in graph, if it is None, select all fields.
        node_fields: `UniprotNodeField` fields that will be included in graph, if it is None, select all fields.
        edge_types: `UniprotEdgeType` fields that will be included in graph, if it is None, select all fields.
        id_fields: `UniprotIDField` field that will be included in graph as node identifier, if it is None, selects first 3 fields.
        add_prefix: if True, add prefix to database identifiers.
        test_mode: if True, limits amount of output data.
    """

    def __init__(
        self,
        organism: Optional[Literal["*"] | int | None] = "*",
        rev: Optional[bool] = True,
        node_types: Optional[Union[list[UniprotNodeType], None]] = None,
        node_fields: Optional[Union[list[UniprotNodeField], None]] = None,
        edge_types: Optional[Union[list[UniprotEdgeType], None]] = None,
        id_fields: Optional[Union[list[UniprotIDField], None]] = None,
        add_prefix: Optional[bool] = True,
        test_mode: Optional[bool] = False,
        assembly_info: Optional[list[dict]] = None,
    ):
        """
        Args:
            assembly_info: Optional list of dicts with keys 'accession', 'strain_name', 'ncbi_taxon_id'.
                When provided, organism nodes use insdc.gcf:<accession> as ID and
                protein-to-organism edges target assemblies instead of taxids.
        """
        model = UniProtModel(
            organism=organism,
            rev=rev,
            node_types=node_types,
            node_fields=node_fields,
            edge_types=edge_types,
            id_fields=id_fields,
            add_prefix=add_prefix,
            test_mode=test_mode,
        ).model_dump()

        # params
        self.organism = model["organism"]
        self.rev = model["rev"]
        self.add_prefix = model["add_prefix"]
        self.test_mode = model["test_mode"]

        # provenance
        self.data_source = "uniprot"
        self.data_version = "2024_03"
        self.data_licence = "CC BY 4.0"

        self._configure_fields()

        self._set_node_and_edge_fields(
            node_types=model["node_types"],
            node_fields=model["node_fields"],
            edge_types=model["edge_types"],
        )

        self.set_id_fields(id_fields=model["id_fields"])

        # Assembly info for assembly-based organism IDs
        self.assembly_info = assembly_info or []

        # loading of ligands and receptors sets
        self.ligands = self._read_ligands_set()
        self.receptors = self._read_receptors_set()

        # loading of subcellular locations set
        self.locations = set()

        # gene property name mappings that will be used gene node properties in KG
        self.gene_property_name_mappings = {"gene_primary":"gene_symbol",
                                            "xref_kegg":"kegg_ids",
                                            }
        # protein property name mappings that will be used protein node properties in KG
        self.protein_property_name_mappings = {
            "protein_name":"protein_names",
            "xref_kegg":"kegg_ids",
        }

    def _read_ligands_set(self) -> set:
        # check if ligands file exists
        if not os.path.isfile("data/ligands_curated.csv"):
            return set()

        ligand_file = pd.read_csv("data/ligands_curated.csv", header=None)
        return set(ligand_file[0])

    def _read_receptors_set(self) -> set:
        # check if receptors file exists
        if not os.path.isfile("data/receptors_curated.csv"):
            return set()

        receptor_file = pd.read_csv("data/receptors_curated.csv", header=None)
        return set(receptor_file[0])

    @validate_call
    def download_uniprot_data(
        self,
        cache: bool = False,
        debug: bool = False,
        retries: int = 3,
        prott5_embedding_output_path: FilePath | None = None,
        esm2_embedding_path: FilePath | None = None,
        nucleotide_transformer_embedding_path: FilePath | None = None,
    ):
        """
        Wrapper function to download uniprot data using pypath; used to access
        settings.

        Args:
            cache: if True, it uses the cached version of the data, otherwise
            forces download.
            debug: if True, turns on debug mode in pypath.
            retries: number of retries in case of download error.
        """

        # stack pypath context managers
        with ExitStack() as stack:

            #stack.enter_context(settings.setup(retries=retries))

            if debug:
                stack.enter_context(curl.debug_on())

            if not cache:
                stack.enter_context(curl.cache_off())

            self._download_uniprot_data(
                prott5_embedding_output_path=prott5_embedding_output_path,
                esm2_embedding_path=esm2_embedding_path,
                nucleotide_transformer_embedding_path=nucleotide_transformer_embedding_path,
            )

            import json
            with open('uniprot_raw_data.json', 'w') as f:
                json.dump(self.data, f, indent=4, default=str)  

            # preprocess data
            self._preprocess_uniprot_data()

            # preprocess organisms
            self._preprocess_organisms()
            
            import json
            with open('uniprot_preprocess_data.json', 'w') as f:
                json.dump(self.data, f, indent=4, default=str)  

    @validate_call
    def _download_uniprot_data(
        self, 
        prott5_embedding_output_path: FilePath | None = None,
        esm2_embedding_path: FilePath | None = None,
        nucleotide_transformer_embedding_path: FilePath | None = None,
    ):
        """
        Download uniprot data from uniprot.org through pypath.

        Here is an overview of uniprot return fields:
        https://www.uniprot.org/help/return_fields

        TODO make use of multi-field query
        """

        logger.info("Downloading uniprot data...")

        t0 = time()

        # download all swissprot ids
        self.uniprot_ids = set(uniprot._all_uniprots(self.organism, self.rev))
        logger.debug(f"found {len(self.uniprot_ids)} uniprot ids")
        # limit to 100 for testing
        if self.test_mode:
            self.uniprot_ids = set(list(self.uniprot_ids)[:100])

        # download attribute dicts
        self.data = {}
        for query_key in tqdm(self.node_fields, desc="Downloading uniprot fields"):
            if query_key in [
                UniprotNodeField.PROTT5_EMBEDDING.value,
                UniprotNodeField.ESM2_EMBEDDING.value,
                UniprotNodeField.NT_EMBEDDING.value,
            ]:
                # downloaded separately and not from the uniprot REST API
                continue

            elif query_key == UniprotNodeField.SUBCELLULAR_LOCATION.value:
                self.data[query_key] = uniprot.uniprot_locations(
                    self.organism, self.rev
                )
            else:
                self.data[query_key] = uniprot.uniprot_data(
                    fields = query_key,
                    organism = self.organism,
                    reviewed = self.rev,
                )

            logger.debug(f"{query_key} field is downloaded")


        if UniprotNodeField.PROTT5_EMBEDDING.value in self.node_fields:
            self.data[UniprotNodeField.PROTT5_EMBEDDING.value] = {}
            self.download_prott5_embeddings(
                prott5_embedding_output_path=prott5_embedding_output_path
            )
        
        if UniprotNodeField.ESM2_EMBEDDING.value in self.node_fields:
            self.data[UniprotNodeField.ESM2_EMBEDDING.value] = {}
            self.retrieve_esm2_embeddings(esm2_embedding_path)

        if UniprotNodeField.NT_EMBEDDING.value in self.node_fields:
            self.data[UniprotNodeField.NT_EMBEDDING.value] = {}
            self.retrieve_nucleotide_transformer_embeddings(
                nucleotide_transformer_embedding_path
            )

        t1 = time()
        msg = f"Acquired UniProt data in {round((t1-t0) / 60, 2)} mins."
        logger.info(msg)

    @validate_call
    def download_prott5_embeddings(
        self, 
        prott5_embedding_output_path: FilePath | None = None
    ):
        """
        Downloads ProtT5 embedding from uniprot website
        If the files exists in a defined file path, then
        directly read it.

        Args:
            prott5_embedding_output_path (FilePath, optional): Defaults to None.
        """
        url: HttpUrl = (
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/uniprot_sprot/per-protein.h5"
        )
        if prott5_embedding_output_path:
            full_path = prott5_embedding_output_path
        else:
            full_path = os.path.join(os.getcwd(), "prott5_protein_embeddings.h5")

        if not os.path.isfile(full_path):
            logger.info("Downloading ProtT5 embeddings...")

            with requests.get(url, stream=True) as response:
                with open(full_path, "wb") as f:
                    for chunk in response.iter_content(512 * 1024):
                        if chunk:
                            f.write(chunk)
        else:
            logger.info("ProtT5 Embedding file is exists. Reading from the file..")

        df_list = []
        with h5py.File(full_path, "r") as file:

            for uniprot_id, embedding in file.items():
                
                embedding = np.array(embedding).astype(np.float16)
                count = np.count_nonzero(~np.isnan(embedding))
                if (
                    self.organism not in ("*", None)
                    and uniprot_id in self.uniprot_ids
                    and count == 1024
                ):
                    # self.data[UniprotNodeField.PROTT5_EMBEDDING.value][
                    #     uniprot_id
                    # ] = embedding
                    df_list.append((uniprot_id, embedding))
                elif count == 1024:
                    # self.data[UniprotNodeField.PROTT5_EMBEDDING.value][
                    #     uniprot_id
                    # ] = embedding
                    df_list.append((uniprot_id, embedding))

        
        self.prott5_embedding_df = pd.DataFrame(df_list, columns=['uniprot_id', 'embedding'])

        del df_list
                    
    @validate_call
    def retrieve_esm2_embeddings(self, 
                                 esm2_embedding_path: FilePath | None = None) -> None:
        
        logger.info("Retrieving ESM2 embeddings...")

        df_list = []
        with h5py.File(esm2_embedding_path, "r") as file:

            for uniprot_id, embedding in file.items():
                embedding = np.array(embedding).astype(np.float16)
                count = np.count_nonzero(~np.isnan(embedding))
                if (
                    self.organism not in ("*", None)
                    and uniprot_id in self.uniprot_ids
                    and count == 1280
                ):
                    # self.data[UniprotNodeField.ESM2_EMBEDDING.value][
                    #     uniprot_id
                    # ] = embedding
                    df_list.append((uniprot_id, embedding))
                elif count == 1280:
                    # self.data[UniprotNodeField.ESM2_EMBEDDING.value][
                    #     uniprot_id
                    # ] = embedding
                    df_list.append((uniprot_id, embedding))

        self.esm2_embedding_df = pd.DataFrame(df_list, columns=['uniprot_id', 'embedding'])

        del df_list

    def retrieve_nucleotide_transformer_embeddings(self,
                                                   nucleotide_transformer_embedding_path: FilePath | None = None) -> None:
        
        logger.info("Retrieving Nucleotide Transformer embeddings...")

        self.entrez_id_to_nucleotide_transformer_embedding = {}
        with h5py.File(nucleotide_transformer_embedding_path, "r") as file:
            for entrez_id, embedding in file.items():
                embedding = np.array(embedding).astype(np.float16)
                self.entrez_id_to_nucleotide_transformer_embedding[entrez_id] = embedding
                
    def _preprocess_uniprot_data(self):
        """
        Preprocess uniprot data to make it ready for import. First, three types
        of processing are applied:
        - nothing is done (for ensembl gene ids, which come from pypath)
        - simple string replacement
        - replace separators in integers and convert to int
        - field splitting

        Then, special treatment is applied to some fields:
        - ensg ids are extracted from the ensembl transcript ids
        - protein names and virus hosts have dedicated normalisation functions
        """

        logger.info("Preprocessing UniProt data.")

        for arg in tqdm(self.node_fields, desc="Processing uniprot fields"):

            # do not process ensembl gene ids (we will get them from pypath)
            # and prott5 embeddings
            if arg in self.nonuniprot_api_fields:
                pass

            elif arg in [
                UniprotNodeField.LENGTH.value,
                UniprotNodeField.MASS.value,
                UniprotNodeField.ORGANISM_ID.value,
            ]:
                for protein, attribute_value in self.data.get(arg).items():
                    self.data[arg][protein] = int(
                        str(attribute_value).replace(",", "")
                    )

            elif arg not in self.split_fields:
                if arg != UniprotNodeField.SUBCELLULAR_LOCATION.value:
                    for protein, attribute_value in self.data.get(arg).items():

                        self.data[arg][protein] = (
                            attribute_value.replace("|", ",")
                            .replace("'", "^")
                            .strip()
                        )

            else:

                for protein, attribute_value in self.data.get(arg).items():
                    # Field splitting
                    self.data[arg][protein] = self._split_fields(
                        arg, attribute_value
                    )


            # Protein names
            if arg == UniprotNodeField.PROTEIN_NAMES.value:

                pass # skip this split - its buggy for now
                # for protein, attribute_value in self.data.get(arg).items():
                #     self.data[arg][protein] = self._split_protein_names_field(
                #         attribute_value
                #     )

            elif arg == UniprotNodeField.SUBCELLULAR_LOCATION.value:
                for protein, attribute_value in self.data.get(arg).items():
                    individual_protein_locations = []
                    for element in attribute_value:
                        loc = (
                            str(element.location)
                            .replace("'", "")
                            .replace("[", "")
                            .replace("]", "")
                            .strip()
                        )
                        individual_protein_locations.append(loc)
                        self.locations.add(loc)

                    self.data[arg][protein] = individual_protein_locations
        # preprocess go fields to extract go ids
        # assume split_fields already applied
        self._preprocess_go_fields()


    def _preprocess_organisms(self):
        """"
        create self.organism_df dataframe for organism nodes
        """
        organism_df = pd.DataFrame.from_dict(self.data[UniprotNodeField.ORGANISM_ID.value], orient='index', columns=['organism_id'])
        for k in self.organism_properties:
            if k in self.data:
                organism_df[k] = organism_df.index.map(self.data[k])
        organism_df = organism_df.drop_duplicates(subset=['organism_id'])
        self.organism_df = organism_df


    def _extract_go_id(self, go_term: str) -> str:
        """
        Extract GO id from GO term string.
        Example input: "aspartate-semialdehyde dehydrogenase activity [GO:0004073]"
        """

        if "GO:" in go_term:
            go_id = go_term.split("GO:")[1].split("]")[0].strip()
            return go_id
        return None
    
    def _preprocess_go_fields(self):
        """
        Preprocess GO fields to extract GO ids from the GO term strings.
        """

        go_fields = [
            UniprotNodeField.CELLULAR_COMPONENT.value,
            UniprotNodeField.BIOLOGICAL_PROCESS.value,
            UniprotNodeField.MOLECULAR_FUNCTION.value,
        ]
        for go_field in go_fields:
            go_id_field = go_field + "_id"
            self.data[go_id_field] = dict()

            for protein, attribute_value in self.data.get(go_field).items():
                # example attribute_value: "aspartate-semialdehyde dehydrogenase activity [GO:0004073]"
                go_ids = [self._extract_go_id(go_term) for go_term in attribute_value if go_term and "GO:" in go_term]
                self.data[go_id_field][protein] = go_ids

    @validate_call
    def _get_ligand_or_receptor(self, uniprot_id: str):
        """
        Tell if UniProt protein node is a L, R or nothing.
        """

        uniprot_id = uniprot_id[8:]

        if uniprot_id in self.ligands:
            return "ligand"
        return "receptor" if uniprot_id in self.receptors else "protein"
    
    def _get_biocypher_property_name(self, original_property_name: str) -> str:
        """
        Get the property name that will be used in the BioCypher KG.
        """

        propname =  self.protein_property_name_mappings.get(original_property_name, original_property_name)
        propname = propname.replace(" ", "_").replace("-", "_")

        return propname

    def get_protein_node(self, protein: str) -> tuple[str, str, dict]:
        """
        Convert a protein to a node representation.
        """
        protein_label = "protein"
        all_props = {arg: self.data.get(arg).get(protein) for arg in self.node_fields}
        protein_id = self.add_prefix_to_id("uniprot", protein)
        protein_props = self._get_protein_properties(all_props, protein_id.split(":")[1])
        return (protein_id, protein_label, protein_props)


    @validate_call
    def get_nodes(
        self, 
        protein_label: str = "protein",
        gene_label: str = "gene",
        organism_label: str = "organism"
    ) -> Generator[tuple[str, str, dict]]:
        """
        Yield nodes (protein, gene, organism) from UniProt data.
        """
        logger.info(
            "Preparing UniProt nodes of the types "
            f"{[type.name for type in self.node_types]}."
        )
        
        node_list = []
        if UniprotNodeType.PROTEIN in self.node_types:
            node_list.extend([self.get_protein_node(protein) for protein in tqdm(self.uniprot_ids)])


        # append organism node to output if desired
        if UniprotNodeType.ORGANISM in self.node_types:
            if self.assembly_info:
                node_list.extend(self._get_assembly_organism_nodes())
            else:
                node_list.extend([self.get_organism_node(i) for i in self.organism_df['organism_id']])

        return node_list
    
    def get_organism_node(self, organism: str) -> tuple[str, str, dict]:
        """
        Get organism node representation from UniProt data.
        When assembly_info is set, this is not used directly — see
        _get_assembly_organism_nodes instead.
        """
        organism_label = "organism"
        organism_id = self.add_prefix_to_id("ncbitaxon",str(organism))

        row = self.organism_df.loc[self.organism_df['organism_id'].isin([organism])].squeeze()
        organism_props = {
            self._get_biocypher_property_name(k): row[k]
            for k in self.organism_df.columns
        }

        # source, licence, and version fields
        organism_props["source"] = self.data_source
        organism_props["licence"] = self.data_licence
        organism_props["version"] = self.data_version

        return (organism_id, organism_label, organism_props)

    def _get_assembly_organism_nodes(self) -> list[tuple]:
        """Create organism nodes keyed by assembly accession (insdc.gcf).

        One node per entry in self.assembly_info.
        """
        nodes = []
        for info in self.assembly_info:
            accession = info['accession']
            node_id = self.add_prefix_to_id("insdc.gcf", accession)
            props = {
                'strain_name': info.get('strain_name', ''),
                'ncbi_taxon_id': info.get('ncbi_taxon_id'),
                'source': self.data_source,
                'licence': self.data_licence,
                'version': self.data_version,
            }
            # Add organism_name from UniProt data if available
            if hasattr(self, 'organism_df') and not self.organism_df.empty:
                org_row = self.organism_df.loc[
                    self.organism_df['organism_id'] == info.get('ncbi_taxon_id')
                ]
                if not org_row.empty:
                    org_name = org_row.iloc[0].get('organism_name')
                    if org_name:
                        props['organism_name'] = org_name
            nodes.append((node_id, "organism", props))
        return nodes
    
    @validate_call
    def _get_organism_nodes(self) -> Generator[tuple[str, str, dict]]:
        """"
        Get organism node representation from UniProt data.
        """
        return[
            self.get_organism_node(i) 
            for i in self.organism_df['organism_id']
        ]
    

    def _create_edges_of_type(self, edge_model: UniprotEdgeModel):
        """
        Create edges of a specific type from UniProt data.
        For PROTEIN_TO_ORGANISM with assembly_info, creates edges to each
        assembly accession (insdc.gcf) instead of taxid.
        """

        # create lists of edges
        edge_list = []

        # generic properties for all edges for now
        properties = {
            "source": self.data_source,
            "licence": self.data_licence,
            "version": self.data_version,
        }
        if (edge_model.edge_type in self.edge_types) and (edge_model.uniprot_field.value in self.data.keys()):
            # Special handling for PROTEIN_TO_ORGANISM when assembly_info is set
            use_assembly = (
                edge_model.edge_type == UniprotEdgeType.PROTEIN_TO_ORGANISM
                and self.assembly_info
            )

            for protein, targets in tqdm(self.data[edge_model.uniprot_field.value].items(), desc=f"Getting {edge_model.edge_type.name} edges"):

                if targets is None or targets == "":
                    continue
                targets = self._ensure_iterable(targets)
                protein_id = self.add_prefix_to_id("uniprot", protein)

                if use_assembly:
                    # Create one edge per assembly accession
                    for info in self.assembly_info:
                        target_id = self.add_prefix_to_id(
                            "insdc.gcf",
                            info['accession'],
                        )
                        edge_list.append(
                            (
                                None,
                                protein_id,
                                target_id,
                                edge_model.edge_label,
                                properties,
                            )
                        )
                else:
                    for target in targets:
                        if not target:
                            continue

                        target_id = self.add_prefix_to_id(
                            edge_model.target_prefix,
                            str(target), # needed for fields that are integers
                        )
                        edge_list.append(
                            (
                                None,
                                protein_id,
                                target_id,
                                edge_model.edge_label,
                                properties,
                            )
                        )
        return edge_list
    


    @validate_call
    def get_edges(self) -> Generator[tuple[None, str, str, str, dict]]:
        """
        Get nodes and edges from UniProt data.
        """

        logger.info(
            "Preparing UniProt edges of the types "
            f"{[type.name for type in self.edge_types]}."
        )

        # create lists of edges
        edge_list = []

        for edge_model in UNIPROT_EDGE_MODELS:
            edges_of_type = self._create_edges_of_type(edge_model)
            edge_list.extend(edges_of_type)


        return edge_list



    @validate_call
    def _get_gene(self, all_props: dict) -> list:
        """
        Get gene node representation from UniProt data per protein. Since one
        protein can have multiple genes, return a list of tuples.
        """

        # if genes and database(GeneID) fields exist, define gene_properties
        if not (
            UniprotNodeField.PROTEIN_GENE_NAMES.value in all_props.keys()
            and UniprotNodeField.ENTREZ_GENE_IDS.value in all_props.keys()
        ):
            return []

        # Find preferred identifier for gene and check if it exists
        if UniprotIDField.GENE_ENTREZ_ID in self.id_fields:

            id_type = UniprotNodeField.ENTREZ_GENE_IDS.value


        gene_raw = all_props.pop(id_type)

        if not gene_raw:
            return []

        type_dict = {
            UniprotNodeField.ENTREZ_GENE_IDS.value: "ncbigene",
        }

        genes = self._ensure_iterable(gene_raw)
        gene_props = {}

        for k in all_props.keys():

            if k not in self.gene_properties:
                continue

            propname = self._get_biocypher_property_name(k)

            # select parenthesis content in field names and make lowercase
            gene_props[propname] = all_props[k]

            if k == UniprotNodeField.NT_EMBEDDING.value and self.entrez_id_to_nucleotide_transformer_embedding.get(genes[0]) is not None:
                gene_props[propname] = [str(emb) for emb in self.entrez_id_to_nucleotide_transformer_embedding[genes[0]]]
        
        # source, licence, and version fields
        gene_props["source"] = self.data_source
        gene_props["licence"] = self.data_licence
        gene_props["version"] = self.data_version

        gene_list = []

        for gene in genes:

            gene_id = self.add_prefix_to_id(
                type_dict[id_type],
                gene,
            )

            gene_list.append((gene_id, gene_props))

        return gene_list

    @validate_call
    def _get_organism(self, all_props: dict):

        organism_props = {}

        organism_id = self.add_prefix_to_id(
            "ncbitaxon",
            str(all_props.pop(UniprotNodeField.ORGANISM_ID.value)),
        )

        for k in all_props.keys():

            if k in self.organism_properties:
                organism_props[k] = all_props[k]

        # source, licence, and version fields
        organism_props["source"] = self.data_source
        organism_props["licence"] = self.data_licence
        organism_props["version"] = self.data_version

        return organism_id, organism_props

    @validate_call
    def _get_protein_properties(self, all_props: dict, protein_id: str) -> dict:
        protein_props = {
            self._get_biocypher_property_name(k): v 
            for k, v in all_props.items()
            if k in self.protein_properties        
        }

        if UniprotNodeField.PROTEIN_NAMES.value in all_props.keys():
                protein_props["primary_protein_name"] = (
                    self._ensure_iterable(all_props[UniprotNodeField.PROTEIN_NAMES.value])[0]
                    if all_props[UniprotNodeField.PROTEIN_NAMES.value]
                    else None
                )


        if UniprotNodeField.PROTT5_EMBEDDING.value in all_props.keys():
            res = self.prott5_embedding_df[self.prott5_embedding_df["uniprot_id"] == protein_id]["embedding"]
            if not res.empty:
            # embedding = np.array(self.prott5_embedding[k]).astype(np.float16)
                propname = self._get_biocypher_property_name(UniprotNodeField.PROTT5_EMBEDDING.value)
                protein_props[propname] = [str(emb) for emb in res.values[0]]

        if UniprotNodeField.ESM2_EMBEDDING.value in all_props.keys():
            res = self.esm2_embedding_df[self.esm2_embedding_df["uniprot_id"] == protein_id]["embedding"]
            if not res.empty:
            # embedding = np.array(self.esm2_embedding[k]).astype(np.float16)
                propname = self._get_biocypher_property_name(UniprotNodeField.ESM2_EMBEDDING.value)
                protein_props[propname] = [str(emb) for emb in res.values[0]]


        # source, licence, and version fields
        protein_props["source"] = self.data_source
        protein_props["licence"] = self.data_licence
        protein_props["version"] = self.data_version

        return protein_props

    def _split_fields(self, field_key, field_value):
        """
        Split fields with multiple entries in uniprot
        Args:
            field_key: field name
            field_value: entry of the field
        """
        if field_value is None or field_value == "":
            return field_value
        
        # replace sensitive elements for admin-import
        field_value = (
            field_value.replace("|", ",").replace("'", "^").strip()
        )

        # define fields that will not be splitted by semicolon
        split_dict = {
            UniprotNodeField.PROTEOME.value: ",",
            UniprotNodeField.PROTEIN_GENE_NAMES.value: " ",
        }

        split_char = split_dict.get(field_key, ";")

        # if field in split_dict split accordingly
        field_value = [i.strip() for i in field_value.strip(split_char).split(split_char)]

        # split colons (":") in kegg field
        if field_key == UniprotNodeField.KEGG_IDS.value:
            _list = [e.split(":")[1].strip() for e in field_value]
            field_value = _list

        # take first element in database(GeneID) field
        if field_key == UniprotNodeField.ENTREZ_GENE_IDS.value:
            field_value = field_value[0]

        # if field has just one element in the list make it string
        if isinstance(field_value, list) and len(field_value) == 1:
            field_value = field_value[0]

        return field_value

    # TODO fix this function to support strings like:
    #    "Q7TU21": "Aspartate-semialdehyde dehydrogenase (ASA dehydrogenase) (ASADH) (EC 1.2.1.11) (Aspartate-beta-semialdehyde dehydrogenase)",
    # "Q7V0G8": "ADP-dependent (S)-NAD(P)H-hydrate dehydratase (EC 4.2.1.136) (ADP-dependent NAD(P)HX dehydratase)",
    # "Q7V0H7": "Multifunctional fusion protein [Includes: 2-C-methyl-D-erythritol 2,4-cyclodiphosphate synthase (MECDP-synthase) (MECPP-synthase) (MECPS) (EC 4.6.1.12); tRNA (guanine-N(1)-)-methyltransferase (EC 2.1.1.228) (M1G-methyltransferase) (tRNA [GM37] methyltransferase)]",
    # "Q7V0W8": "Protein nucleotidyltransferase YdiU (EC 2.7.7.-) (Protein adenylyltransferase YdiU) (EC 2.7.7.108) (Protein uridylyltransferase YdiU) (EC 2.7.7.-)",
    # "Q7V1H8": "Riboflavin biosynthesis protein RibBA [Includes: 3,4-dihydroxy-2-butanone 4-phosphate synthase (DHBP synthase) (EC 4.1.99.12); GTP cyclohydrolase-2 (EC 3.5.4.25) (GTP cyclohydrolase II)]",
    # "Q7V1L1": "ATP-dependent zinc metalloprotease FtsH (EC 3.4.24.-)",
    # "Q7V1N8": "Carbamoyl phosphate synthase large chain (EC 6.3.4.16) (EC 6.3.5.5) (Carbamoyl phosphate synthetase ammonia chain)",
    # "Q7V359": "Coenzyme A biosynthesis bifunctional protein CoaBC (DNA/pantothenate metabolism flavoprotein) (Phosphopantothenoylcysteine synthetase/decarboxylase) (PPCS-PPCDC) [Includes: Phosphopantothenoylcysteine decarboxylase (PPC decarboxylase) (PPC-DC) (EC 4.1.1.36) (CoaC); Phosphopantothenate--cysteine ligase (EC 6.3.2.5) (CoaB) (Phosphopantothenoylcysteine synthetase) (PPC synthetase) (PPC-S)]",
    # "Q7V362": "ATP-dependent zinc metalloprotease FtsH (EC 3.4.24.-)",
    # "Q7TU15": "Threonine synthase (EC 4.2.3.1)",
    # "Q7TU17": "Polyphosphate kinase (EC 2.7.4.1) (ATP-polyphosphate phosphotransferase) (Polyphosphoric acid kinase)",

    def _split_protein_names_field(self, field_value):
        """
        Split protein names field in uniprot
        Args:
            field_value: entry of the protein names field
        Example:
            "Acetate kinase (EC 2.7.2.1) (Acetokinase)" -> ["Acetate kinase", "Acetokinase"]
        """
        field_value = field_value.replace("|", ",").replace(
            "'", "^"
        )  # replace sensitive elements

        if "[Cleaved" in field_value:
            # discarding part after the "[Cleaved"
            clip_index = field_value.index("[Cleaved")
            protein_names = (
                field_value[:clip_index].replace("(Fragment)", "").strip()
            )

            # handling multiple protein names
            if "(EC" in protein_names[0]:
                splitted = protein_names[0].split(" (")
                protein_names = []

                for name in splitted:
                    if not name.strip().startswith("EC") and not name.strip().startswith("Fragm"):
                        protein_names.append(name.rstrip(")").strip())

            elif " (" in protein_names[0]:
                splitted = protein_names[0].split(" (")
                protein_names = [
                    name.rstrip(")").strip()
                    for name in splitted
                    if not name.strip().startswith("Fragm")
                ]
        elif "[Includes" in field_value:
            # discarding part after the "[Includes"
            clip_index = field_value.index("[Includes")
            protein_names = (
                field_value[:clip_index].replace("(Fragment)", "").strip()
            )
            # handling multiple protein names
            if "(EC" in protein_names[0]:

                splitted = protein_names[0].split(" (")
                protein_names = [
                    name.rstrip(")").strip()
                    for name in splitted
                    if not name.strip().startswith("EC")
                    and not name.strip().startswith("Fragm")
                ]
            elif " (" in protein_names[0]:
                splitted = protein_names[0].split(" (")
                protein_names = [
                    name.rstrip(")").strip()
                    for name in splitted
                    if not name.strip().startswith("Fragm")
                ]
        elif "(EC" in field_value.replace("(Fragment)", ""):
            splitted = field_value.split(" (")
            protein_names = [
                name.rstrip(")").strip()
                for name in splitted
                if not name.strip().startswith("EC")
                and not name.strip().startswith("Fragm")
            ]
        elif " (" in field_value.replace("(Fragment)", ""):
            splitted = field_value.split(" (")
            protein_names = [
                name.rstrip(")").strip()
                for name in splitted
                if not name.strip().startswith("Fragm")
            ]
        else:
            protein_names = field_value.replace("(Fragment)", "").strip()

        return protein_names

    def _find_ensg_from_enst(self, enst_list):
        """
        take ensembl transcript ids, return ensembl gene ids by using pypath mapping tool

        Args:
            field_value: ensembl transcript list

        """
        if enst_list is None:
            return None, None
        
        enst_list = self._ensure_iterable(enst_list)

        enst_list = [enst.split(" [")[0] for enst in enst_list]

        ensg_ids = set()
        for enst_id in enst_list:
            ensg_id = list(
                mapping.map_name(
                    enst_id.split(".")[0], "enst_biomart", "ensg_biomart"
                )
            )
            ensg_id = ensg_id[0] if ensg_id else None
            if ensg_id:
                ensg_ids.add(ensg_id)

        ensg_ids = list(ensg_ids)

        if len(ensg_ids) == 1:
            ensg_ids = ensg_ids[0]

        if len(enst_list) == 1:
            enst_list = enst_list[0]

        return enst_list, ensg_ids

    @lru_cache
    def _normalise_curie_cached(
        self, prefix: str, identifier: str, sep: str = ":"
    ) -> Optional[str]:
        """
        Wrapper to call and cache `normalize_curie()` from Bioregistry.
        """

        if not self.normalise_curies:
            return identifier

        return normalize_curie(f"{prefix}{sep}{identifier}", sep=sep)

    @lru_cache
    @validate_call
    def add_prefix_to_id(
        self, prefix: str = None, identifier: str = None, sep: str = ":"
    ) -> str:
        """
        Adds prefix to database id
        """
        if self.add_prefix and identifier:
            return normalize_curie(prefix + sep + identifier)

        return identifier

    def _configure_fields(self):
        # fields that need splitting
        self.split_fields = UniprotNodeField.get_split_fields()

        # fields that are not downloaded from uniprot REST API
        self.nonuniprot_api_fields = UniprotNodeField.get_nonuniprot_api_fields()

        # properties of nodes
        self.protein_properties = UniprotNodeField.get_protein_properties()

        self.gene_properties = UniprotNodeField.get_gene_properties()

        self.organism_properties = UniprotNodeField.get_organism_properties()

    def _set_node_and_edge_fields(
        self,
        node_types,
        node_fields,
        edge_types,
    ):

        # ensure computation of ENSGs

        # check which node types and fields to include
        self.node_types = node_types or list(UniprotNodeType)

        self.node_fields = [field.value for field in node_fields] if node_fields else [field.value for field in UniprotNodeField]

        # check which edge types and fields to include
        self.edge_types = edge_types or list(UniprotEdgeType)

    def set_id_fields(self, id_fields):
        self.id_fields = id_fields or list(UniprotIDField)[:3]

    def _ensure_iterable(self, value):
        return value if isinstance(value, list) else [value]

    @validate_call
    def export_data_to_csv(
        self,
        node_data: Generator[tuple[str, str, dict]] = None,
        edge_data: Generator[tuple[None, str, str, str, dict]] = None,
        path: DirectoryPath | None = None,
    ) -> None:
        """
        Save node and edge data to csv
            node_data: output of `get_nodes()` function
            edge_data: output of `get_edges()` function
            path: Directory to save the output csv file
        """
        if node_data:
            logger.debug("Saving uniprot node data as csv")
            node_types_dict = collections.defaultdict(list)
            for _id, _type, props in node_data:
                _dict = {"id": _id} | props
                node_types_dict[_type].append(_dict)

            for _type, values in node_types_dict.items():
                df = pd.DataFrame.from_records(values)
                if path:
                    full_path = os.path.join(path, f"{_type.capitalize()}.csv")
                else:
                    full_path = os.path.join(
                        os.getcwd(), f"{_type.capitalize()}.csv"
                    )

                df.to_csv(full_path, index=False)
                logger.info(
                    f"{_type.capitalize()} data is written: {full_path}"
                )

        if edge_data:
            logger.debug("Saving uniprot edge data as csv")
            edge_types_dict = collections.defaultdict(list)
            for _, source_id, target_id, _type, props in edge_data:
                _dict = {"source_id": source_id, "target_id": target_id} | props
                edge_types_dict[_type].append(_dict)

            for _type, values in edge_types_dict.items():
                df = pd.DataFrame.from_records(values)
                if path:
                    full_path = os.path.join(path, f"{_type.capitalize()}.csv")
                else:
                    full_path = os.path.join(
                        os.getcwd(), f"{_type.capitalize()}.csv"
                    )

                df.to_csv(full_path, index=False)
                logger.info(
                    f"{_type.capitalize()} data is written: {full_path}"
                )


class MultiUniprot:
    """Wrapper that reads a CSV file listing organisms and runs the Uniprot
    adapter for each unique taxid.  Multiple assemblies sharing the same
    taxid are grouped so that a single download covers all of them, while
    protein-to-organism edges fan out to every assembly."""

    def __init__(self, config_list_file: str, **kwargs):
        """
        Args:
            config_list_file: Path to a CSV file with columns:
                ncbi_accession, cyanorak_organism, ncbi_taxon_id, strain_name, data_dir
                Lines starting with # are treated as comments and skipped.
                One Uniprot adapter is created per *unique* ncbi_taxon_id.
                Assembly info (accession, strain_name) is passed through so
                organism nodes use insdc.gcf IDs.
            **kwargs: Additional arguments passed to each Uniprot adapter.
        """
        self.adapters = []
        self.organism_ids = []

        # Build taxid -> [assembly_info] mapping, deduplicating by taxid
        from collections import OrderedDict
        taxid_to_assemblies: dict[int, list[dict]] = OrderedDict()

        with open(config_list_file, 'r') as f:
            # Filter out comment lines (starting with #) before parsing CSV
            lines = [line for line in f if not line.strip().startswith('#')]
            reader = csv.DictReader(lines)

            for row in reader:
                organism_id = int(row['ncbi_taxon_id'])
                info = {
                    'accession': row['ncbi_accession'],
                    'strain_name': row.get('strain_name') or '',
                    'ncbi_taxon_id': organism_id,
                }
                taxid_to_assemblies.setdefault(organism_id, []).append(info)

        # Create one adapter per unique taxid, passing all assembly info
        for organism_id, assembly_list in taxid_to_assemblies.items():
            self.organism_ids.append(organism_id)
            adapter = Uniprot(
                organism=organism_id,
                assembly_info=assembly_list,
                **kwargs,
            )
            self.adapters.append(adapter)

        logger.info(
            f"Loaded {len(self.adapters)} unique organism(s) "
            f"({sum(len(v) for v in taxid_to_assemblies.values())} assemblies) "
            f"from {config_list_file}"
        )

    def download_uniprot_data(self, **kwargs):
        """Download UniProt data for all organisms."""
        for adapter in self.adapters:
            adapter.download_uniprot_data(**kwargs)

    def get_nodes(self, **kwargs):
        """Get nodes from all adapters, deduplicating organism nodes."""
        nodes = []
        seen_organism_ids = set()

        for adapter in self.adapters:
            adapter_nodes = adapter.get_nodes(**kwargs)
            for node in adapter_nodes:
                node_id, label, props = node
                # Deduplicate organism nodes (they may be shared across adapters)
                if label == "organism":
                    if node_id in seen_organism_ids:
                        continue
                    seen_organism_ids.add(node_id)
                nodes.append(node)
        return nodes

    def get_edges(self, **kwargs):
        """Get edges from all adapters."""
        edges = []
        for adapter in self.adapters:
            edges.extend(adapter.get_edges(**kwargs))
        return edges
