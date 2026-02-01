from __future__ import annotations
from pypath.share import curl, settings

from contextlib import ExitStack

from bioregistry import normalize_curie
from tqdm import tqdm
from time import time
from biocypher._logger import logger

import collections
import csv
import os
import h5py
import gzip

from pydantic import BaseModel, DirectoryPath, FilePath, EmailStr, validate_call
from typing import Literal, Union

from enum import Enum, EnumMeta, auto

import pandas as pd
import numpy as np
from Bio import SeqIO
import gffpandas.gffpandas as gffpd
from urllib.parse import unquote
import re

logger.debug(f"Loading module {__name__}.")


class GeneEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class GeneNodeField(Enum, metaclass=GeneEnumMeta):

    GENE_NAMES = 'gene_names'
    GENE_NAMES_CYANORAK = 'gene_names_cyanorak'
    LOCUS_TAG = 'locus_tag'
    LOCUS_TAG_NCBI = 'locus_tag_ncbi'
    LOCUS_TAG_CYANORAK = 'locus_tag_cyanorak'
    START = 'start'
    END = 'end'
    START_CYANORAK = 'start_cyanorak'
    END_CYANORAK = 'end_cyanorak'
    STRAND = 'strand'
    STRAND_CYANORAK = 'strand_cyanorak'
    PRODUCT = 'product'
    PRODUCT_CYANORAK = 'product_cyanorak'
    PROTEIN_ID = 'protein_id'
    ONTOLOGY_TERM = 'Ontology_term'
    ONTOLOGY_TERM_DESCRIPTION = 'ontology_term_description'
    EGGNOG = 'eggNOG'
    EGGNOG_DESCRIPTION = 'eggNOG_description'
    KEGG = 'kegg'
    KEGG_DESCRIPTION = 'kegg_description'
    CYANORAK_ROLE = 'cyanorak_Role'
    CYANORAK_ROLE_DESCRIPTION = 'cyanorak_Role_description'
    TIGR_ROLE = 'tIGR_Role'
    TIGR_ROLE_DESCRIPTION = 'tIGR_Role_description'
    CLUSTER_NUMBER = 'cluster_number'
    PROTEIN_DOMAINS = 'protein_domains'
    PROTEIN_DOMAINS_DESCRIPTION = 'protein_domains_description'


    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class GeneEdgeType(Enum, metaclass=GeneEnumMeta):
    #TODO
    GENE_TO_PROTEIN = auto()


class GeneModel(BaseModel):
    gene_node_fields: Union[list[GeneNodeField], None] = None
    edge_types: Union[list[GeneEdgeType], None] = None
    test_mode: bool = False
    export_csv: bool = False
    output_dir: DirectoryPath | None = None
    add_prefix: bool = True
    organism: int | Literal["*"] | None = None


class CyanorakNcbi:
    def __init__(
        self,
        gene_node_fields: Union[list[GeneNodeField], None] = None,
        edge_types: Union[list[GeneEdgeType], None] = None,
        test_mode: bool = False,
        export_csv: bool = False,
        output_dir: DirectoryPath | None = None,
        add_prefix: bool = True,
        organism: int | Literal["*"] | None = None,
        ncbi_gff_file: str = None,
        cyan_gff_file: str = None,
        cyan_gbk_file: str = None,
    ):

        model = GeneModel(
            gene_node_fields=gene_node_fields,
            edge_types=edge_types,
            test_mode=test_mode,
            export_csv=export_csv,
            output_dir=output_dir,
            add_prefix=add_prefix,
            organism=organism,
        ).model_dump()

        self.export_csv = model["export_csv"]
        self.output_dir = model["output_dir"]
        self.add_prefix = model["add_prefix"]

        self.cyan_gff_file = cyan_gff_file
        self.cyan_gbk_file = cyan_gbk_file
        self.ncbi_gff_file = ncbi_gff_file

        # no need becuase we are not creating protein to ec edges here
        # if model["organism"] in ("*", None):
        #     self.swissprots = set(uniprot._all_uniprots("*", True))
        # else:
        #     self.swissprots = set(uniprot._all_uniprots(model["organism"], True))

        # set node fields
        self.set_node_fields(gene_node_fields=model["gene_node_fields"])

        # set edge types
        self.set_edge_types(edge_types=model["edge_types"])

        # set early_stopping, if test_mode true
        self.early_stopping = None
        if model["test_mode"]:
            self.early_stopping = 100

    @validate_call
    def download_data(self, cache: bool = False) -> None:
        ''' Download the cyanorak gbk file and ncbi gff file from the provided URLs.'''
        self.data_df = self.load_gff_from_ncbi_and_cynorak(
            ncbi_gff_file=self.ncbi_gff_file,
            cyan_gff_file=self.cyan_gff_file,
            cyan_gbk_file=self.cyan_gbk_file
        )

    @validate_call
    def get_nodes(self, label: str = "gene") -> list[tuple]:

        logger.info("Started writing gene nodes")

        node_list = []
        for _, row in self.data_df.iterrows():
            node_properties = {}
            for field in self.gene_node_fields:
                value = row.get(field)
                if not pd.isna(value):
                    if isinstance(value, str):
                        value = unquote(value)
                    value = self.clean_text(value)
                    value = self._split_field(field, value)
                    node_properties[field] = value

            node_id = self.add_prefix_to_id(
                prefix="ncbigene",
                identifier=row.get("locus_tag"),
            )

            node_list.append((node_id, label, node_properties))

        logger.info(f"Finished writing {len(node_list)} gene nodes")
        return node_list

    @validate_call
    def get_edges(self) -> list[tuple]:
        # for now, no edges to create
        edge_list = []
        logger.info("Started writing gene edges")
        return edge_list




    def _get_cyanorak_cols_to_keep(self) -> list[str]:
        cyan_cols_to_keep = [
            'start', 'end',  'strand', 
            'ID', 'Name', 'Ontology_term', 'cluster_number',
            'cyanorak_Role', 'cyanorak_Role_description', 'eggNOG',
            'eggNOG_description', 'kegg', 'kegg_description',
            'ontology_term_description', 'product', 'protein_domains',
            'protein_domains_description', 'tIGR_Role', 'tIGR_Role_description',
            'locus_tag']
        return cyan_cols_to_keep

    def _get_ncbi_cols_to_keep_map(self) -> dict[str, str]:
        ncbi_cols_to_keep = [
        'Name_gene',
        'gene_gene',
        'locus_tag_cds',
        'old_locus_tag_gene',
        'source_cds',
        'start_cds',
        'end_cds',
        'strand_cds',
        'Note_cds',
        'exception_cds',
        'inference_cds',
        'product_cds',
        'protein_id_cds']
        col_rename_map = {c: c.replace('_gene', '').replace('_cds', '') for c in ncbi_cols_to_keep}
        col_rename_map['locus_tag_cds'] = 'locus_tag_ncbi'
        col_rename_map['old_locus_tag_gene'] = 'locus_tag'
        return col_rename_map

    def _get_final_merged_columns_map(self) -> dict[str, str]:
        return {
            'Name_ncbi' : 'gene_names', 
            'Name_cyanorak': 'gene_names_cyanorak',
            'start_ncbi': 'start', 'end_ncbi': 'end', 'strand_ncbi': 'strand', 
            'start_cyanorak': 'start_cyanorak', 'end_cyanorak': 'end_cyanorak', 'strand_cyanorak': 'strand_cyanorak',
            'product_ncbi': 'product', 'product_cyanorak': 'product_cyanorak',
            'ID': 'locus_tag_cyanoak',
            }

    def _get_split_character(self, field: str) -> str:
        ''' get split character for multiple values in a text field. Return None if not a multivalue field.'''
        comma_split_cols = [
            'cyanorak_Role', 'cyanorak_Role_description',
            'Ontology_term',  'ontology_term_description',
            'eggNOG', 'eggNOG_description', 
            'kegg', #'kegg_description',  
            'protein_domains', 'protein_domains_description', 
            'tIGR_Role', 'tIGR_Role_description'
        ]
        space_split_cols = [
            'gene_names', 'gene', 
        ]
        semmicommas_split_cols = [
            'kegg', 'kegg_description'
        ]

        if field in comma_split_cols:
            return ','
        if field in space_split_cols:
            return ' '
        if field in semmicommas_split_cols:
            return ';'
        return None
    
    def _split_field(self, field: str, value: str) -> list[str] | str:
        ''' Clean up and split a text field based on its expected split character.'''
        split_char = self._get_split_character(field)
        if split_char is not None:
            repattern = re.escape(split_char) + r'(?! )'  # split_char not followed by space
            value = [v.strip() for v in re.split(repattern, str(value))]
        return value
        
    def _get_cynaorak_ID(self, rec) -> str:
        ''' Extract the Cyanorak ID from a GenBank record of cyanorak gbk file.'''
        note_name = 'cyanorak ORF Id:'
        cyanorak_ID =[i for i in rec.qualifiers['note'] if i.startswith(note_name)]
        if len(cyanorak_ID) > 1:
            print(f"Warning: multiple cyanorak IDs found for record {rec.id} - {cyanorak_ID}")
        return cyanorak_ID[0].replace(note_name, '').strip()


    def _get_cyanorak_id_map_from_gbk(self, gbk_file: str) -> dict[str, str]:
        ''' Create a mapping from locus_tag to Cyanorak ID from a GenBank file of cyanorak gbk file.'''
        seq_records = [rec for rec in SeqIO.read(gbk_file, "genbank").features if rec.type in ["CDS"]]
        seq_records_map = {self._get_cynaorak_ID(rec) : locus_tag for rec in seq_records for locus_tag in rec.qualifiers['locus_tag']}
        return seq_records_map


    def ncbi_merge_cds_and_gene_entries(self, ncbi_gff_df: pd.DataFrame) -> pd.DataFrame:
        ''' Merge gene and CDS entries from NCBI GFF DataFrame based on ID and Parent attributes.'''
        gene_df = ncbi_gff_df.loc[ncbi_gff_df.type.isin(['gene'])]
        cds_df = ncbi_gff_df.loc[ncbi_gff_df.type.isin(['CDS'])]
        df = pd.merge(gene_df, cds_df, left_on='ID', right_on='Parent', suffixes=['_gene', '_cds'], )
        ncbi_col_rename_map = self._get_ncbi_cols_to_keep_map()
        df_filter = df[list(ncbi_col_rename_map.keys())].rename(columns=ncbi_col_rename_map)
        return df_filter

    def load_gff(self, gff_file: str) -> pd.DataFrame:
        ''' Load a GFF3 file and return a DataFrame object.'''
        gff_annotation = gffpd.read_gff3(gff_file)
        gff_df = gff_annotation.attributes_to_columns()
        gff_cds_df = gff_df.loc[gff_df.type.isin(['CDS', 'gene'])]
        return gff_cds_df

    def load_gff_from_ncbi_and_cynorak(self, ncbi_gff_file: str, cyan_gff_file: str, cyan_gbk_file: str) -> pd.DataFrame:
        ''' Load a GFF3 file from NCBI and map the locus_tag to Cyanorak ID using the provided GenBank file.'''
        cyan_df = self.load_gff(cyan_gff_file)
        ncbi_df = self.load_gff(ncbi_gff_file)
        ncbi_df = self.ncbi_merge_cds_and_gene_entries(ncbi_df)
        cyanID2locus_tag_map = self._get_cyanorak_id_map_from_gbk(cyan_gbk_file)
        cyan_df['locus_tag'] = cyan_df['ID'].map(cyanID2locus_tag_map)
        cyan_df = cyan_df[self._get_cyanorak_cols_to_keep()]
        merge_df = pd.merge(ncbi_df, cyan_df, on='locus_tag', suffixes=['_ncbi', '_cyanorak'])  
        merge_df = merge_df.rename(columns=self._get_final_merged_columns_map())
        return merge_df

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

    def set_node_fields(self, gene_node_fields):
        if gene_node_fields:
            self.gene_node_fields = [field.value for field in gene_node_fields]
        else:
            self.gene_node_fields = [field.value for field in GeneNodeField]

    def set_edge_types(self, edge_types):
        self.edge_types = edge_types or list(GeneEdgeType)

    def clean_text(
        self, text: str = None,
    ) -> str:
        """
        remove biocypher special characters from text fields
        """
        special_chars_map = {
            "|": ",", # pipe used to separate multiple values in biocypher
            "'": "^", # single quote used to quote strings in biocypher
        }
        if isinstance(text, str):
            for char, replacement in special_chars_map.items():
                text = text.replace(char, replacement)
            return text
        elif isinstance(text, list):
            return [self.clean_text(t) for t in text]
        else:
            return text
        

class MultiCyanorakNcbi:
    """Wrapper that reads a CSV file listing cyanobacteria genome file paths
    and delegates to CyanorakNcbi instances."""

    def __init__(self, config_list_file: str, **kwargs):
        """
        Args:
            config_list_file: Path to a CSV file with columns:
                genome_dir, ncbi_gff, cyan_gff, cyan_gbk
            **kwargs: Additional arguments passed to each CyanorakNcbi.
        """
        self.adapters = []
        with open(config_list_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                genome_dir = row['genome_dir']
                adapter = CyanorakNcbi(
                    ncbi_gff_file=os.path.join(genome_dir, row['ncbi_gff']),
                    cyan_gff_file=os.path.join(genome_dir, row['cyan_gff']),
                    cyan_gbk_file=os.path.join(genome_dir, row['cyan_gbk']),
                    **kwargs,
                )
                self.adapters.append(adapter)
        logger.info(f"Loaded {len(self.adapters)} genome configs from {config_list_file}")

    def download_data(self, **kwargs):
        for adapter in self.adapters:
            adapter.download_data(**kwargs)

    def get_nodes(self):
        nodes = []
        for adapter in self.adapters:
            nodes.extend(adapter.get_nodes())
        return nodes

    def get_edges(self):
        edges = []
        for adapter in self.adapters:
            edges.extend(adapter.get_edges())
        return edges


if __name__ == "__main__":
    # code for testing
    dpath = 'data/Prochlorococcus/genomes/MED4/'
    cgff_fpath = os.path.join(dpath, 'cyanorak/Pro_MED4.gff')
    ngff_fpath = os.path.join(dpath, 'genomic.gff')
    cgbk_fpath = os.path.join(dpath, 'cyanorak/Pro_MED4.gbk')

    adapter = CyanorakNcbi(
        ncbi_gff_file=ngff_fpath,
        cyan_gff_file=cgff_fpath,
        cyan_gbk_file=cgbk_fpath,
    )
    adapter.download_data()
    nodes = adapter.get_nodes()
    import json
    with open('cyanorak_ncbi_gene_nodes.json', 'w') as f:
        json.dump(nodes, f, indent=4)

    # print(adapter.data_df.columns.tolist())
    # for col in adapter.data_df.columns:
    #     commas = adapter.data_df[col].astype(str).str.contains('%3b')
    #     semmicommas = adapter.data_df[col].astype(str).str.contains('%2c')

    #     spaces = adapter.data_df[col].astype(str).str.contains(' ')
    #     if commas.any():
    #         print(f"Column {col} contains commas")

    #     elif semmicommas.any():
    #         print(f"Column {col} contains semmicommas")

# columns in the final merged dataframe:
# ['Name_ncbi', 'gene', 'locus_tag_ncbi', 'locus_tag', 'source',
#        'start_ncbi', 'end_ncbi', 'strand_ncbi', 'Note', 'exception',
#        'inference', 'product_ncbi', 'protein_id', 'start_cyanorak',
#        'end_cyanorak', 'strand_cyanorak', 'ID', 'Name_cyanorak',
#        'Ontology_term', 'cluster_number', 'cyanorak_Role',
#        'cyanorak_Role_description', 'eggNOG', 'eggNOG_description', 'kegg',
#        'kegg_description', 'ontology_term_description', 'product_cyanorak',
#        'protein_domains', 'protein_domains_description', 'tIGR_Role',
#        'tIGR_Role_description']