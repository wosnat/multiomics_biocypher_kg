from __future__ import annotations
from pypath.share import curl

from bioregistry import normalize_curie
from biocypher._logger import logger

import collections
import csv
import json
import os
import xml.etree.ElementTree as ET

from pydantic import BaseModel, DirectoryPath, validate_call
from typing import Literal, Union

from enum import Enum, EnumMeta, auto

import pandas as pd
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
    CYANORAK_ROLE = 'cyanorak_Role'
    CYANORAK_ROLE_DESCRIPTION = 'cyanorak_Role_description'
    TIGR_ROLE = 'tIGR_Role'
    TIGR_ROLE_DESCRIPTION = 'tIGR_Role_description'
    CLUSTER_NUMBER = 'cluster_number'
    EC_NUMBERS = 'ec_numbers'
    FUNCTION_DESCRIPTION = 'function_description'
    # Gene naming
    GENE_NAME = 'gene_name'
    GENE_SYNONYMS = 'gene_synonyms'
    GENE_NAME_SOURCE = 'gene_name_source'
    PRODUCT_SOURCE = 'product_source'
    FUNCTION_DESCRIPTION_SOURCE = 'function_description_source'
    OLD_LOCUS_TAGS = 'old_locus_tags'
    # EggNOG / COG
    COG_CATEGORY = 'cog_category'
    EGGNOG_OGS = 'eggnog_ogs'
    EGGNOG_OG_DESCRIPTIONS = 'eggnog_og_descriptions'
    SEED_ORTHOLOG = 'seed_ortholog'
    MAX_ANNOT_LVL = 'max_annot_lvl'
    SEED_ORTHOLOG_EVALUE = 'seed_ortholog_evalue'
    PROTEIN_FAMILY = 'protein_family'
    # Gene Ontology
    GO_TERMS = 'go_terms'
    GO_TERM_DESCRIPTIONS = 'go_term_descriptions'
    # KEGG
    KEGG_KO = 'kegg_ko'
    KEGG_KO_DESCRIPTIONS = 'kegg_ko_descriptions'
    KEGG_PATHWAY = 'kegg_pathway'
    KEGG_MODULE = 'kegg_module'
    KEGG_REACTION = 'kegg_reaction'
    KEGG_BRITE = 'kegg_brite'
    # Pfam domains
    PFAM_IDS = 'pfam_ids'
    PFAM_DESCRIPTIONS = 'pfam_descriptions'
    # Specialized function
    CATALYTIC_ACTIVITY = 'catalytic_activity'
    TRANSMEMBRANE_REGIONS = 'transmembrane_regions'
    SIGNAL_PEPTIDE = 'signal_peptide'
    TRANSPORTER_CLASSIFICATION = 'transporter_classification'
    CAZY_IDS = 'cazy_ids'
    BIGG_REACTION = 'bigg_reaction'
    # Quality
    ANNOTATION_QUALITY = 'annotation_quality'


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
        ncbi_accession: str = None,
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

        self.ncbi_accession = ncbi_accession
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

    @validate_call
    def download_data(
        self,
        cache: bool = False,
        debug: bool = False,
        retries: int = 3,
    ) -> None:
        """Load gene annotations from pre-built gene_annotations_merged.json.

        Run 'bash scripts/prepare_data.sh' first to generate this file.
        """
        if not self.data_dir:
            raise ValueError("data_dir is required")
        merged_json_path = os.path.join(self.data_dir, "gene_annotations_merged.json")
        if not os.path.exists(merged_json_path):
            raise FileNotFoundError(
                f"gene_annotations_merged.json not found at {merged_json_path}. "
                f"Run 'bash scripts/prepare_data.sh' first."
            )
        with open(merged_json_path) as f:
            merged_data = json.load(f)
        self.data_df = pd.DataFrame.from_dict(merged_data, orient='index').reset_index(drop=True)
        logger.info(
            f"Loaded {len(self.data_df)} genes for {self.strain_name} "
            f"from gene_annotations_merged.json"
        )

        # Fetch taxonomy for organism node creation (cached to data_dir)
        if self.ncbi_taxon_id:
            self.taxonomy = _fetch_ncbi_taxonomy(
                taxid=self.ncbi_taxon_id,
                cache_dir=self.data_dir,
            )

    def _get_gene_nodes(self) -> list[tuple]:
        """Generate gene nodes from the data.

        Returns:
            List of tuples: (node_id, label, properties_dict)
        """
        logger.info("Started writing gene nodes")

        int_fields = {'start', 'end', 'start_cyanorak', 'end_cyanorak', 'annotation_quality'}
        float_fields = {'seed_ortholog_evalue'}
        node_list = []
        for _, row in self.data_df.iterrows():
            node_properties = {}
            for field in self.gene_node_fields:
                value = row.get(field)
                if pd.isna(value) if not isinstance(value, list) else False:
                    continue
                if isinstance(value, list):
                    # Already a list from merged JSON — clean text, no splitting needed
                    node_properties[field] = self.clean_text(value)
                elif field in int_fields:
                    node_properties[field] = int(value)
                elif field in float_fields:
                    node_properties[field] = float(value)
                elif isinstance(value, str):
                    value = unquote(value)
                    value = self.clean_text(value)
                    value = self._split_field(field, value)
                    node_properties[field] = value
                else:
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
        """Generate gene → cluster and gene → organism edges.

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

        # Gene → organism edges (one per gene, for all genes with a locus_tag)
        if self.ncbi_accession:
            organism_id = self.add_prefix_to_id(
                prefix="insdc.gcf",
                identifier=self.ncbi_accession,
            )
            for _, row in self.data_df.iterrows():
                locus_tag = row.get('locus_tag')
                if pd.isna(locus_tag):
                    continue
                gene_id = self.add_prefix_to_id(prefix="ncbigene", identifier=locus_tag)
                edge_list.append((
                    f"{gene_id}_belongs_to_{organism_id}",
                    gene_id,
                    organism_id,
                    "gene_belongs_to_organism",
                    {}
                ))
            logger.info(f"Finished writing gene-organism edges for {self.ncbi_accession}")

        return edge_list




    def _get_split_character(self, field: str) -> str:
        ''' get split character for multi-value string fields. Return None if not a multi-value field.
        Note: most fields come as Python lists from gene_annotations_merged.json and bypass this method.
        This is retained for any remaining comma-delimited string fields.
        '''
        comma_split_cols = [
            'cyanorak_Role', 'cyanorak_Role_description',
            'tIGR_Role', 'tIGR_Role_description',
            'old_locus_tags', 'ec_numbers',
        ]
        if field in comma_split_cols:
            return ','
        return None
    
    def _split_field(self, field: str, value: str) -> list[str] | str:
        ''' Clean up and split a text field based on its expected split character.'''
        split_char = self._get_split_character(field)
        if split_char is not None:
            repattern = re.escape(split_char) + r'(?! )'  # split_char not followed by space
            value = [v.strip() for v in re.split(repattern, str(value))]
        return value
        
    @validate_call
    def add_prefix_to_id(
        self, prefix: str | None = None, identifier: str | None = None, sep: str = ":"
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
    """Wrapper that reads a CSV file listing cyanobacteria genome configs
    and delegates to CyanorakNcbi instances."""

    def __init__(self, config_list_file: str, treatment_organisms_file: str = None, **kwargs):
        """
        Args:
            config_list_file: Path to a CSV file with columns:
                ncbi_accession, data_dir, strain_name, ncbi_taxon_id, clade
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

            for row in reader:
                taxon_id_str = row.get('ncbi_taxon_id')
                ncbi_taxon_id = int(taxon_id_str) if taxon_id_str else None
                clade = row.get('clade') or None
                adapter = CyanorakNcbi(
                    ncbi_accession=row['ncbi_accession'],
                    data_dir=row.get('data_dir') or None,
                    strain_name=row.get('strain_name') or None,
                    ncbi_taxon_id=ncbi_taxon_id,
                    clade=clade,
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