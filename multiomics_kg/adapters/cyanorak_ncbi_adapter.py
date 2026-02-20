from __future__ import annotations
from pypath.share import curl, settings

from contextlib import ExitStack

from bioregistry import normalize_curie
from tqdm import tqdm
from time import time
from biocypher._logger import logger

import collections
import csv
import json
import os
import h5py
import gzip
import xml.etree.ElementTree as ET

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

# NCBI rank names → our property names
_NCBI_RANK_TO_PROP = {
    'superkingdom': 'superkingdom',
    'domain': 'superkingdom',   # NCBI now uses 'domain' for Bacteria / Eukaryota
    'kingdom': 'kingdom',
    'phylum': 'phylum',
    'class': 'tax_class',       # avoid collision with Python built-in
    'order': 'order',
    'family': 'family',
    'genus': 'genus',
    'species': 'species',
}


def _fetch_ncbi_taxonomy(taxid: int, cache_dir: str) -> dict:
    """Fetch taxonomy lineage from NCBI efetch API with file-based caching.

    Uses the same pypath curl library as the genome downloads.

    Args:
        taxid: NCBI taxonomy ID.
        cache_dir: Directory for caching the result as taxonomy_<taxid>.json.

    Returns:
        Dict with keys: lineage (full semicolon string), superkingdom, kingdom,
        phylum, tax_class, order, family, genus, species — only populated
        ranks are included.
    """
    os.makedirs(cache_dir, exist_ok=True)
    cache_path = os.path.join(cache_dir, f"taxonomy_{taxid}.json")

    if os.path.exists(cache_path):
        logger.info(f"Taxonomy cache hit: {cache_path}")
        with open(cache_path) as f:
            return json.load(f)

    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=taxonomy&id={taxid}&rettype=xml&retmode=xml"
    )
    c = curl.Curl(url, silent=False)
    if c.result is None:
        logger.warning(f"Failed to fetch taxonomy for taxid {taxid} from {url}")
        return {}

    result = {}
    try:
        root = ET.fromstring(c.result)
        taxon = root.find('Taxon')
        if taxon is not None:
            result['lineage'] = taxon.findtext('Lineage') or ''
            lineage_ex = taxon.find('LineageEx')
            if lineage_ex is not None:
                for t in lineage_ex.findall('Taxon'):
                    rank = (t.findtext('Rank') or '').lower()
                    name = t.findtext('ScientificName') or ''
                    prop = _NCBI_RANK_TO_PROP.get(rank)
                    if prop and name:
                        result[prop] = name
    except ET.ParseError as e:
        logger.warning(f"Failed to parse taxonomy XML for taxid {taxid}: {e}")
        return {}

    with open(cache_path, 'w') as f:
        json.dump(result, f, indent=2)
    logger.info(f"Taxonomy saved to {cache_path}")
    return result


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
    OLD_LOCUS_TAGS = 'old_locus_tags'
    PROTEIN_DOMAINS = 'protein_domains'
    PROTEIN_DOMAINS_DESCRIPTION = 'protein_domains_description'


    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class ClusterNodeField(Enum, metaclass=GeneEnumMeta):
    CLUSTER_NUMBER = 'cluster_number'

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
    GENE_IN_CLUSTER = auto()


class GeneModel(BaseModel):
    gene_node_fields: Union[list[GeneNodeField], None] = None
    cluster_node_fields: Union[list[ClusterNodeField], None] = None
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
        cluster_node_fields: Union[list[ClusterNodeField], None] = None,
        edge_types: Union[list[GeneEdgeType], None] = None,
        test_mode: bool = False,
        export_csv: bool = False,
        output_dir: DirectoryPath | None = None,
        add_prefix: bool = True,
        organism: int | Literal["*"] | None = None,
        ncbi_gff_file: str = None,
        cyan_gff_file: str = None,
        cyan_gbk_file: str = None,
        ncbi_accession: str = None,
        cyanorak_organism: str = None,
        data_dir: str = None,
        strain_name: str = None,
        ncbi_taxon_id: int = None,
        clade: str = None,
    ):

        model = GeneModel(
            gene_node_fields=gene_node_fields,
            cluster_node_fields=cluster_node_fields,
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

        self.ncbi_accession = ncbi_accession
        self.cyanorak_organism = cyanorak_organism
        self.data_dir = data_dir
        self.strain_name = strain_name
        self.ncbi_taxon_id = ncbi_taxon_id
        self.clade = clade
        self.taxonomy = {}  # populated by download_data()

        # no need becuase we are not creating protein to ec edges here
        # if model["organism"] in ("*", None):
        #     self.swissprots = set(uniprot._all_uniprots("*", True))
        # else:
        #     self.swissprots = set(uniprot._all_uniprots(model["organism"], True))

        # set node fields
        self.set_node_fields(gene_node_fields=model["gene_node_fields"])

        # set cluster node fields
        self.set_cluster_node_fields(cluster_node_fields=model["cluster_node_fields"])

        # set edge types
        self.set_edge_types(edge_types=model["edge_types"])

        # set early_stopping, if test_mode true
        self.early_stopping = None
        if model["test_mode"]:
            self.early_stopping = 100

    def _download_ncbi_genome(self) -> str:
        """Download NCBI genome zip and extract GFF using pypath Curl."""
        gff_path = os.path.join(self.data_dir, "genomic.gff")
        if os.path.exists(gff_path):
            logger.info(f"NCBI GFF already exists at {gff_path}, skipping download")
            return gff_path
        url = (
            f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
            f"{self.ncbi_accession}/download"
            f"?include_annotation_type=GENOME_GFF"
            f"&include_annotation_type=GENOME_FASTA"
            f"&include_annotation_type=PROT_FASTA"
            f"&include_annotation_type=SEQUENCE_REPORT"
            f"&hydrated=FULLY_HYDRATED"
        )
        c = curl.Curl(
            url,
            silent=False,
            compr='zip',
            large=False,
        )
        if c.result is None:
            raise ConnectionError(
                f"Failed to download NCBI genome for {self.ncbi_accession} from {url}"
            )
        # c.result is a dict of {filename: content} for zip files
        gff_content = None
        for name, content in c.result.items():
            if name.endswith('genomic.gff'):
                gff_content = content
                break
        if gff_content is None:
            raise ValueError(f"No genomic.gff found in NCBI download for {self.ncbi_accession}")

        os.makedirs(self.data_dir, exist_ok=True)
        with open(gff_path, 'w') as f:
            f.write(gff_content)
        logger.info(f"NCBI GFF saved to {gff_path}")
        return gff_path

    def _download_cyanorak_gff(self) -> str:
        """Download Cyanorak GFF annotation file."""
        cyan_dir = os.path.join(self.data_dir, "cyanorak")
        os.makedirs(cyan_dir, exist_ok=True)
        gff_path = os.path.join(cyan_dir, f"{self.cyanorak_organism}.gff")
        if os.path.exists(gff_path):
            logger.info(f"Cyanorak GFF already exists at {gff_path}, skipping download")
            return gff_path
        url = f"https://cyanorak.sb-roscoff.fr/cyanorak/svc/export/organism/gff/{self.cyanorak_organism}"
        c = curl.Curl(url, silent=False)
        if c.result is None:
            raise ConnectionError(
                f"Failed to download Cyanorak GFF for {self.cyanorak_organism} from {url}"
            )
        with open(gff_path, 'w') as f:
            f.write(c.result)
        logger.info(f"Cyanorak GFF saved to {gff_path}")
        return gff_path

    def _download_cyanorak_gbk(self) -> str:
        """Download Cyanorak GenBank annotation file."""
        cyan_dir = os.path.join(self.data_dir, "cyanorak")
        os.makedirs(cyan_dir, exist_ok=True)
        gbk_path = os.path.join(cyan_dir, f"{self.cyanorak_organism}.gbk")
        if os.path.exists(gbk_path):
            logger.info(f"Cyanorak GBK already exists at {gbk_path}, skipping download")
            return gbk_path
        url = f"https://cyanorak.sb-roscoff.fr/cyanorak/svc/export/organism/gbk/{self.cyanorak_organism}"
        c = curl.Curl(url, silent=False)
        if c.result is None:
            raise ConnectionError(
                f"Failed to download Cyanorak GBK for {self.cyanorak_organism} from {url}"
            )
        with open(gbk_path, 'w') as f:
            f.write(c.result)
        logger.info(f"Cyanorak GBK saved to {gbk_path}")
        return gbk_path

    @validate_call
    def download_data(
        self,
        cache: bool = False,
        debug: bool = False,
        retries: int = 3,
    ) -> None:
        """Download genome files (if using accession-based config) and load GFF data.

        Args:
            cache: if True, uses pypath cached version; if False, forces download.
            debug: if True, turns on debug mode in pypath.
            retries: number of retries in case of download error.
        """
        # If accession-based, download files using pypath Curl
        if self.ncbi_accession:
            with ExitStack() as stack:
                if debug:
                    stack.enter_context(curl.debug_on())
                if not cache:
                    stack.enter_context(curl.cache_off())

                self.ncbi_gff_file = self._download_ncbi_genome()
                if self.cyanorak_organism:
                    self.cyan_gff_file = self._download_cyanorak_gff()
                    self.cyan_gbk_file = self._download_cyanorak_gbk()
                if self.ncbi_taxon_id:
                    self.taxonomy = _fetch_ncbi_taxonomy(
                        taxid=self.ncbi_taxon_id,
                        cache_dir=self.data_dir,
                    )

        if not self.ncbi_gff_file:
            raise ValueError(
                "NCBI GFF file is required. Provide ncbi_gff_file directly "
                "or set ncbi_accession."
            )

        # Branch: full merge (NCBI+Cyanorak) or NCBI-only
        has_cyanorak = all([self.cyan_gff_file, self.cyan_gbk_file])
        if has_cyanorak:
            self.data_df = self.load_gff_from_ncbi_and_cynorak(
                ncbi_gff_file=self.ncbi_gff_file,
                cyan_gff_file=self.cyan_gff_file,
                cyan_gbk_file=self.cyan_gbk_file,
            )
        else:
            logger.info("No Cyanorak data available; loading NCBI GFF only")
            self.data_df = self.load_gff_from_ncbi_only(
                ncbi_gff_file=self.ncbi_gff_file,
            )

        # Save mapping table to data_dir for use by the omics adapter
        mapping_dir = self.data_dir or os.path.dirname(self.ncbi_gff_file)
        os.makedirs(mapping_dir, exist_ok=True)
        mapping_path = os.path.join(mapping_dir, "gene_mapping.csv")
        self.data_df.to_csv(mapping_path, index=False)
        logger.info(f"Gene mapping table saved to {mapping_path}")

    def _get_gene_nodes(self) -> list[tuple]:
        """Generate gene nodes from the data.

        Returns:
            List of tuples: (node_id, label, properties_dict)
        """
        logger.info("Started writing gene nodes")

        node_list = []
        for _, row in self.data_df.iterrows():
            node_properties = {}
            for field in self.gene_node_fields:
                value = row.get(field)
                if not pd.isna(value):
                    if isinstance(value, (float, int)) and field in ('start', 'end', 'start_cyanorak', 'end_cyanorak'):
                        node_properties[field] = int(value)
                        continue
                    if isinstance(value, str):
                        value = unquote(value)
                    value = self.clean_text(value)
                    value = self._split_field(field, value)
                    node_properties[field] = value

            node_id = self.add_prefix_to_id(
                prefix="ncbigene",
                identifier=row.get("locus_tag"),
            )

            node_list.append((node_id, "gene", node_properties))

        logger.info(f"Finished writing {len(node_list)} gene nodes")
        return node_list

    def _get_cluster_nodes(self) -> list[tuple]:
        """Generate cluster nodes from unique cluster_number values in gene data.

        Returns:
            List of tuples: (node_id, label, properties_dict)
        """
        logger.info("Started writing cyanorak cluster nodes")

        node_list = []

        # Check if cluster_number column exists (not present in NCBI-only mode)
        if 'cluster_number' not in self.data_df.columns:
            logger.info("No cluster_number column in data - skipping cluster nodes")
            return node_list

        # Get unique cluster numbers, excluding NaN values
        cluster_numbers = self.data_df['cluster_number'].dropna().unique()

        for cluster_num in cluster_numbers:
            cluster_num_str = str(cluster_num).strip()
            if not cluster_num_str:
                continue

            node_id = self.add_prefix_to_id(
                prefix="cyanorak.cluster",
                identifier=cluster_num_str,
            )

            node_properties = {}
            for field in self.cluster_node_fields:
                if field == 'cluster_number':
                    node_properties[field] = cluster_num_str

            node_list.append((node_id, "cyanorak_cluster", node_properties))

        logger.info(f"Finished writing {len(node_list)} cyanorak cluster nodes")
        return node_list

    def _get_organism_node(self) -> list[tuple]:
        """Generate organism node for this assembly.

        Returns:
            List with a single tuple: (node_id, label, properties_dict)
        """
        if not self.ncbi_accession:
            return []

        node_id = self.add_prefix_to_id(
            prefix="insdc.gcf",
            identifier=self.ncbi_accession,
        )
        properties = {}
        if self.strain_name:
            properties['strain_name'] = self.strain_name
            properties['organism_name'] = self.strain_name
        if self.ncbi_taxon_id is not None:
            properties['ncbi_taxon_id'] = self.ncbi_taxon_id
        if self.clade:
            properties['clade'] = self.clade
        # Taxonomy ranks fetched from NCBI efetch (populated by download_data)
        for key in ('lineage', 'superkingdom', 'kingdom', 'phylum', 'tax_class',
                    'order', 'family', 'genus', 'species'):
            if self.taxonomy.get(key):
                properties[key] = self.taxonomy[key]

        logger.info(f"Created organism node {node_id} (strain: {self.strain_name})")
        return [(node_id, "organism", properties)]

    @validate_call
    def get_nodes(self) -> list[tuple]:
        """Generate all nodes (genes, clusters, and organism).

        Returns:
            List of tuples: (node_id, label, properties_dict)
        """
        node_list = []
        node_list.extend(self._get_organism_node())
        node_list.extend(self._get_gene_nodes())
        node_list.extend(self._get_cluster_nodes())
        return node_list

    @validate_call
    def get_edges(self) -> list[tuple]:
        """Generate gene → cluster edges.

        Returns:
            List of tuples: (edge_id, source_id, target_id, edge_type, properties_dict)
        """
        logger.info("Started writing gene edges")
        edge_list = []

        if GeneEdgeType.GENE_IN_CLUSTER not in self.edge_types:
            logger.info("GENE_IN_CLUSTER edge type not in edge_types, skipping gene-cluster edges")
            return edge_list

        for _, row in self.data_df.iterrows():
            cluster_num = row.get('cluster_number')
            locus_tag = row.get('locus_tag')

            # Skip if no cluster number or locus tag
            if pd.isna(cluster_num) or pd.isna(locus_tag):
                continue

            cluster_num_str = str(cluster_num).strip()
            if not cluster_num_str:
                continue

            gene_id = self.add_prefix_to_id(
                prefix="ncbigene",
                identifier=locus_tag,
            )
            cluster_id = self.add_prefix_to_id(
                prefix="cyanorak.cluster",
                identifier=cluster_num_str,
            )
            edge_id = f"{gene_id}_in_{cluster_id}"

            edge_list.append((
                edge_id,
                gene_id,
                cluster_id,
                "gene_in_cyanorak_cluster",
                {}  # no additional properties for now
            ))

        logger.info(f"Finished writing {len(edge_list)} gene-cluster edges")
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
            'tIGR_Role', 'tIGR_Role_description',
            'old_locus_tags',
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
        # Only keep columns that exist in the merged DataFrame (e.g. old_locus_tag_gene
        # may be absent for genomes that lack old_locus_tag in their GFF)
        available_cols = {k: v for k, v in ncbi_col_rename_map.items() if k in df.columns}
        df_filter = df[list(available_cols.keys())].rename(columns=available_cols)

        # Handle multiple old_locus_tag values (URL-encoded comma-separated in GFF)
        # e.g. "PMT0003%2CPMT_0003%2CRG24_RS00015" -> ["PMT0003", "PMT_0003", "RG24_RS00015"]
        if 'locus_tag' in df_filter.columns:
            df_filter['locus_tag'] = df_filter['locus_tag'].apply(
                lambda x: unquote(str(x)) if pd.notna(x) else x
            )
            # Store all old_locus_tags as a comma-separated string
            df_filter['old_locus_tags'] = df_filter['locus_tag']
        return df_filter

    def load_gff_from_ncbi_only(self, ncbi_gff_file: str) -> pd.DataFrame:
        """Load gene data from NCBI GFF3 file only (no Cyanorak data).

        Reuses load_gff() and ncbi_merge_cds_and_gene_entries(), then renames
        columns to match GeneNodeField values expected by get_nodes().
        """
        ncbi_df = self.load_gff(ncbi_gff_file)
        ncbi_df = self.ncbi_merge_cds_and_gene_entries(ncbi_df)

        # Same locus_tag priority as the merge path: prefer old_locus_tag,
        # fall back to locus_tag_ncbi (locus_tag may not exist if GFF lacks old_locus_tag)
        if 'locus_tag' in ncbi_df.columns:
            ncbi_df['locus_tag'] = ncbi_df['locus_tag'].fillna(ncbi_df['locus_tag_ncbi'])
        else:
            ncbi_df['locus_tag'] = ncbi_df['locus_tag_ncbi']

        # Rename to match GeneNodeField.GENE_NAMES
        ncbi_df = ncbi_df.rename(columns={'Name': 'gene_names'})

        # Drop rows without a usable locus_tag
        n_before = len(ncbi_df)
        ncbi_df = ncbi_df.dropna(subset=['locus_tag'])
        n_dropped = n_before - len(ncbi_df)
        if n_dropped > 0:
            logger.warning(f"Dropped {n_dropped} genes with no identifiable locus_tag")

        logger.info(
            f"NCBI-only gene load: {len(ncbi_df)} genes from {ncbi_gff_file}"
        )
        return ncbi_df

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
        # Drop Cyanorak entries without a locus_tag (no GBK mapping available)
        cyan_df = cyan_df.dropna(subset=['locus_tag'])

        # Handle multiple old_locus_tag values: explode to one row per tag for merge.
        # old_locus_tags column (comma-separated string) preserves all values.
        # If GFF lacks old_locus_tag, fall back to locus_tag_ncbi
        if 'locus_tag' not in ncbi_df.columns:
            ncbi_df['locus_tag'] = ncbi_df['locus_tag_ncbi']
        ncbi_df['locus_tag'] = ncbi_df['locus_tag'].str.split(',')
        ncbi_exploded = ncbi_df.explode('locus_tag')
        ncbi_exploded['locus_tag'] = ncbi_exploded['locus_tag'].str.strip()

        # Outer merge to include genes from either source, not just intersection
        merge_df = pd.merge(
            ncbi_exploded, cyan_df, on='locus_tag', how='outer',
            suffixes=['_ncbi', '_cyanorak'],
        )

        # --- Deduplication ---
        # NCBI genes may have multiple old_locus_tags producing multiple rows
        # after explode. Some tags may match Cyanorak and others may not.
        # Strategy: split into NCBI-sourced vs Cyanorak-only, dedup each
        # independently, then recombine.
        has_ncbi = merge_df['locus_tag_ncbi'].notna()
        ncbi_sourced = merge_df[has_ncbi].copy()
        cyan_only = merge_df[~has_ncbi].copy()

        # For NCBI genes: prefer rows that also matched Cyanorak (have ID).
        # Sort so matched rows (non-NaN Cyanorak ID) come first, then dedup
        # on locus_tag_ncbi to keep one row per NCBI gene.
        if not ncbi_sourced.empty:
            ncbi_sourced = ncbi_sourced.sort_values(
                by='ID', ascending=True, na_position='last',
            )
            dup_mask = ncbi_sourced.duplicated(subset=['locus_tag_ncbi'], keep=False)
            if dup_mask.any():
                dup_genes = ncbi_sourced.loc[dup_mask, 'locus_tag_ncbi'].unique()
                logger.warning(
                    f"{len(dup_genes)} NCBI genes have duplicate rows after merge. "
                    f"Keeping best match for each."
                )
            ncbi_sourced = ncbi_sourced.drop_duplicates(
                subset=['locus_tag_ncbi'], keep='first',
            )

        # For Cyanorak-only genes: remove any whose locus_tag was already
        # captured via an NCBI gene, then dedup on locus_tag.
        if not cyan_only.empty:
            matched_locus_tags = set(
                ncbi_sourced.loc[ncbi_sourced['ID'].notna(), 'locus_tag'].dropna()
            ) if not ncbi_sourced.empty else set()
            cyan_only = cyan_only[~cyan_only['locus_tag'].isin(matched_locus_tags)]
            cyan_only = cyan_only.drop_duplicates(subset=['locus_tag'], keep='first')

        # Compute counts before concat for logging
        n_ncbi_matched = int(ncbi_sourced['ID'].notna().sum()) if not ncbi_sourced.empty else 0
        n_ncbi_only = len(ncbi_sourced) - n_ncbi_matched
        n_cyan_only = len(cyan_only)

        merge_df = pd.concat([ncbi_sourced, cyan_only], ignore_index=True)

        # For NCBI genes without old_locus_tag, locus_tag is NaN — use
        # locus_tag_ncbi so they still get a valid node ID.
        merge_df['locus_tag'] = merge_df['locus_tag'].fillna(merge_df['locus_tag_ncbi'])

        # Drop rows where locus_tag is still NaN (shouldn't happen in practice)
        n_before = len(merge_df)
        merge_df = merge_df.dropna(subset=['locus_tag'])
        n_dropped = n_before - len(merge_df)
        if n_dropped > 0:
            logger.warning(f"Dropped {n_dropped} genes with no identifiable locus_tag")

        logger.info(
            f"Gene merge: {n_ncbi_matched} matched both sources, "
            f"{n_ncbi_only} NCBI-only, {n_cyan_only} Cyanorak-only, "
            f"{len(merge_df)} total"
        )

        merge_df = merge_df.rename(columns=self._get_final_merged_columns_map())
        return merge_df

    @validate_call
    def add_prefix_to_id(
        self, prefix: str = None, identifier: str = None, sep: str = ":"
    ) -> str:
        """
        Adds prefix to database id. Falls back to underscore-separated ID if
        prefix is not registered in bioregistry.
        """
        if self.add_prefix and identifier:
            curie = prefix + sep + identifier
            normalized = normalize_curie(curie)
            # Fall back to underscore-separated ID if normalize_curie returns None (unregistered prefix)
            if normalized is not None:
                return normalized
            return prefix + "_" + identifier

        return identifier

    def set_node_fields(self, gene_node_fields):
        if gene_node_fields:
            self.gene_node_fields = [field.value for field in gene_node_fields]
        else:
            self.gene_node_fields = [field.value for field in GeneNodeField]

    def set_cluster_node_fields(self, cluster_node_fields):
        if cluster_node_fields:
            self.cluster_node_fields = [field.value for field in cluster_node_fields]
        else:
            self.cluster_node_fields = [field.value for field in ClusterNodeField]

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

    def __init__(self, config_list_file: str, treatment_organisms_file: str = None, **kwargs):
        """
        Args:
            config_list_file: Path to a CSV file with columns:
                New format: ncbi_accession, cyanorak_organism, data_dir
                Legacy format: genome_dir, ncbi_gff, cyan_gff, cyan_gbk
                Lines starting with # are treated as comments and skipped.
            treatment_organisms_file: Optional path to a CSV file with columns:
                ncbi_taxon_id, organism_name
                Creates organism-only nodes (ncbitaxon:<taxid>) for treatment
                organisms that don't have loaded genomes.
            **kwargs: Additional arguments passed to each CyanorakNcbi.
        """
        self.adapters = []
        self.treatment_organisms_file = treatment_organisms_file
        with open(config_list_file, 'r') as f:
            # Filter out comment lines (starting with #) before parsing CSV
            lines = [line for line in f if not line.strip().startswith('#')]
            reader = csv.DictReader(lines)
            fieldnames = reader.fieldnames

            for row in reader:
                if 'ncbi_accession' in fieldnames:
                    # New accession-based format
                    # Treat empty cyanorak_organism as None (NCBI-only mode)
                    cyanorak_org = row.get('cyanorak_organism') or None
                    if cyanorak_org and cyanorak_org.strip() == '':
                        cyanorak_org = None
                    # Parse ncbi_taxon_id if present
                    taxon_id_str = row.get('ncbi_taxon_id')
                    ncbi_taxon_id = int(taxon_id_str) if taxon_id_str else None
                    clade = row.get('clade') or None
                    adapter = CyanorakNcbi(
                        ncbi_accession=row['ncbi_accession'],
                        cyanorak_organism=cyanorak_org,
                        data_dir=row.get('data_dir') or None,
                        strain_name=row.get('strain_name') or None,
                        ncbi_taxon_id=ncbi_taxon_id,
                        clade=clade,
                        **kwargs,
                    )
                else:
                    # Legacy path-based format
                    genome_dir = row['genome_dir']
                    cyan_gff = row.get('cyan_gff') or None
                    cyan_gbk = row.get('cyan_gbk') or None
                    adapter = CyanorakNcbi(
                        ncbi_gff_file=os.path.join(genome_dir, row['ncbi_gff']),
                        cyan_gff_file=os.path.join(genome_dir, cyan_gff) if cyan_gff else None,
                        cyan_gbk_file=os.path.join(genome_dir, cyan_gbk) if cyan_gbk else None,
                        **kwargs,
                    )
                self.adapters.append(adapter)
        logger.info(f"Loaded {len(self.adapters)} genome configs from {config_list_file}")

    def _get_treatment_organism_nodes(self) -> list[tuple]:
        """Create organism nodes from treatment_organisms_file.

        Returns nodes with ID ``ncbitaxon:<taxid>`` for genus-level or
        non-genome treatment organisms referenced in paperconfigs.
        """
        if not self.treatment_organisms_file:
            return []
        taxonomy_cache_dir = "cache/data/taxonomy"
        nodes = []
        with open(self.treatment_organisms_file, 'r') as f:
            lines = [line for line in f if not line.strip().startswith('#')]
            reader = csv.DictReader(lines)
            for row in reader:
                taxid = int(row['ncbi_taxon_id'])
                node_id = f"ncbitaxon:{taxid}"
                props = {
                    'organism_name': row.get('organism_name', ''),
                    'ncbi_taxon_id': taxid,
                }
                taxonomy = _fetch_ncbi_taxonomy(taxid=taxid, cache_dir=taxonomy_cache_dir)
                for key in ('lineage', 'superkingdom', 'kingdom', 'phylum', 'tax_class',
                            'order', 'family', 'genus', 'species'):
                    if taxonomy.get(key):
                        props[key] = taxonomy[key]
                nodes.append((node_id, "organism", props))
                logger.info(f"Created treatment organism node {node_id} ({row.get('organism_name', '')})")
        return nodes

    def download_data(self, **kwargs):
        for adapter in self.adapters:
            adapter.download_data(**kwargs)

    def get_nodes(self):
        nodes = []
        for adapter in self.adapters:
            nodes.extend(adapter.get_nodes())
        nodes.extend(self._get_treatment_organism_nodes())
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