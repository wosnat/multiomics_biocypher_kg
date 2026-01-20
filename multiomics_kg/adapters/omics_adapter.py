from __future__ import annotations
from pypath.share import curl, settings

from contextlib import ExitStack

from bioregistry import normalize_curie
from tqdm import tqdm
from time import time
from biocypher._logger import logger

import collections
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


class OmicsEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class PublicationNodeField(Enum, metaclass=OmicsEnumMeta):
    TITLE = 'title'
    AUTHORS = 'authors'
    JOURNAL = 'journal'
    PUBLICATION_DATE = 'publication_date'
    DOI = 'doi'
    PMID = 'pmid'
    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None

class StudyNodeField(Enum, metaclass=OmicsEnumMeta):
    TITLE = 'title'
    DESCRIPTION = 'description'
    ABSTRACT = 'abstract'
    STUDY_TYPE = 'study_type' # e.g., Transcriptomics, Proteomics
    ORGANISM = 'organism'

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class StatisticalTestNodeField(Enum, metaclass=OmicsEnumMeta):
    NAME = 'name'
    TEST_TYPE = 'test_type'
    CONTROL_CONDITION = 'control_condition'
    TREATMENT_CONDITION = 'treatment_condition'
    TIMEPOINT = 'timepoint'
    REFERENCE_TIMEPOINT = 'reference_timepoint'
    METHOD = 'method'

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None
    

class TimeSeriesNodeField(Enum, metaclass=OmicsEnumMeta):
    NAME = 'name'
    METHOD = 'method' # e.g., WGCNA, k-means
    DESCRIPTION = 'description'

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None



class OMICSEdgeType(Enum, metaclass=OmicsEnumMeta):
    study_published_in = auto()
    test_in_study = auto()
    study_has_time_series = auto()
    molecular_result_from_test = auto()
    cluster_in_study = auto()
    molecular_in_cluster = auto()


class OMICSModel(BaseModel):
    publication_node_fields: Union[list[PublicationNodeField], None] = None
    edge_types: Union[list[OMICSEdgeType], None] = None
    test_mode: bool = False
    export_csv: bool = False
    output_dir: DirectoryPath | None = None
    add_prefix: bool = True
    organism: int | Literal["*"] | None = None


class OMICSAdapter:
    def __init__(
        self,
        publication_node_fields: Union[list[PublicationNodeField], None] = None,
        edge_types: Union[list[OMICSEdgeType], None] = None,
        test_mode: bool = False,
        export_csv: bool = False,
        output_dir: DirectoryPath | None = None,
        add_prefix: bool = True,
        organism: int | Literal["*"] | None = None,
        ncbi_gff_file: str = None,
        cyan_gff_file: str = None,
        cyan_gbk_file: str = None,
    ):

        model = OMICSModel(
            publication_node_fields=publication_node_fields,
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



        # set node fields
        self.set_node_fields(publication_node_fields=model["publication_node_fields"])

        # set edge types
        self.set_edge_types(edge_types=model["edge_types"])

        # set early_stopping, if test_mode true
        self.early_stopping = None
        if model["test_mode"]:
            self.early_stopping = 100

    @validate_call
    def download_data(self, cache: bool = False) -> None:
        """Download necessary data files for the adapter.

        Args:
            cache (bool, optional): Whether to cache downloaded files. Defaults to False.
        """
        pass
    
