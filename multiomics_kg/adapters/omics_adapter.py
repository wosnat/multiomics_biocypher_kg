from __future__ import annotations

from biocypher._logger import logger
from pydantic import BaseModel, DirectoryPath, validate_call
from typing import Literal, Union, Optional
from enum import Enum, EnumMeta, auto
import math
import pandas as pd
import yaml
import os
from pathlib import Path
from bioregistry import normalize_curie

try:
    from multiomics_kg.adapters.pdf_publication_extraction import PDFPublicationExtractor
except ImportError:
    from pdf_publication_extraction import PDFPublicationExtractor


logger.debug(f"Loading module {__name__}.")


class OmicsEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()

class PublicationNodeType(Enum):
    """Define types of nodes the adapter can provide."""
    PUBLICATION = auto()
    TIME_SERIES_CLUSTER = auto()


class PublicationNodeField(Enum, metaclass=OmicsEnumMeta):
    """Fields for publication nodes."""
    PUBMED_ID = "pubmed_id"
    TITLE = "title"
    AUTHORS = "authors"
    JOURNAL = "journal"
    PUBLICATION_DATE = "publication_date"
    DOI = "doi"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class TimeSeriesNodeField(Enum, metaclass=OmicsEnumMeta):
    """Fields for timeseries cluster nodes."""
    NAME = 'name'
    METHOD = 'method'  # e.g., WGCNA, k-means
    DESCRIPTION = 'description'

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class OMICSEdgeType(Enum, metaclass=OmicsEnumMeta):
    """Edge types for omics data."""
    affects_expression_of = auto()
    cluster_in_publication = auto()
    molecular_in_cluster = auto()


class OMICSModel(BaseModel):
    """Pydantic model for OMICS adapter configuration."""
    edge_types: Union[list[OMICSEdgeType], None] = None
    test_mode: bool = False
    export_csv: bool = False
    output_dir: DirectoryPath | None = None
    add_prefix: bool = True
    organism: int | Literal["*"] | None = None


class OMICSAdapter:
    """
    Adapter for creating omics-related nodes and edges from publication config files.

    Parses YAML paperconfig files that specify statistical analyses (RNA-seq, proteomics,
    metabolomics) and creates publication, organism, and environmental condition nodes,
    plus affects_expression_of edges linking causes to genes.
    """

    def __init__(
        self,
        config_file: str = None,
        edge_types: Union[list[OMICSEdgeType], None] = None,
        test_mode: bool = False,
        export_csv: bool = False,
        output_dir: DirectoryPath | None = None,
        add_prefix: bool = True,
        organism: int | Literal["*"] | None = None,
    ):
        """
        Initialize the OMICS adapter.
        
        Args:
            config_file: Path to paperconfig.yaml file
            edge_types: List of edge types to include
            test_mode: If True, limit output for testing
            export_csv: If True, export nodes and edges as CSV
            output_dir: Directory for CSV output
            add_prefix: Whether to add prefix to identifiers
            organism: Filter by organism (not used for now)
        """
        model = OMICSModel(
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
        self.test_mode = model["test_mode"]

        # set edge types
        self.set_edge_types(edge_types=model["edge_types"])

        # set early_stopping, if test_mode true
        self.early_stopping = None
        if model["test_mode"]:
            self.early_stopping = 100

        # Load config file
        self.config_file = config_file
        self.config_data = {}
        if config_file:
            self.load_config(config_file)

        # pdf extractor for publication nodes
        self.pdf_extractor = PDFPublicationExtractor()

        self.data_source = "OMICS Adapter"
        self.data_version = "2026-01-01"
        self.data_licence = "CC BY 4.0"

    def set_edge_types(self, edge_types: Union[list[OMICSEdgeType], None]) -> None:
        """Set the edge types to include in the result."""
        if isinstance(edge_types, list):
            self.edge_types = edge_types
        else:
            self.edge_types = [
                OMICSEdgeType.affects_expression_of,
            ]

    def load_config(self, config_file: str) -> None:
        """
        Load and parse the paperconfig.yaml file.
        
        Args:
            config_file: Path to paperconfig.yaml file
        """
        if not os.path.exists(config_file):
            logger.error(f"Config file not found: {config_file}")
            return

        with open(config_file, 'r') as f:
            self.config_data = yaml.safe_load(f)
            logger.info(f"Loaded config from {config_file}")

    @validate_call
    def download_data(self, cache: bool = False) -> None:
        """Download necessary data files for the adapter.

        Args:
            cache (bool, optional): Whether to cache downloaded files. Defaults to False.
        """
        # Extract publication metadata
        publication = self.config_data.get('publication', {})
        if not isinstance(publication, dict):
            logger.warning(f"'publication' must be a dict, got {type(publication).__name__}. Skipping all tables.")
            return 

        # extract publication metadata from pdf
        pdf_path = publication.get('papermainpdf', None)
        if pdf_path and os.path.exists(pdf_path):
            self.extracted_data = self.pdf_extractor.extract_from_pdf(pdf_path)
        else:
            logger.warning(f"PDF path not found or not provided: {pdf_path}. Skipping publication nodes.")

    def get_publication_nodes(self) -> list[tuple]:
        """
        Generate publication nodes from cached data (includes study metadata).

        Returns:
            List of (node_id, label, properties) tuples
        """
        logger.info("Generating publication nodes from cached data")
        nodes = []

        data = self.extracted_data if hasattr(self, 'extracted_data') else {}

        if "publication" in data:
            pub = data["publication"]
            pub_id = self.get_publication_id()
            pub_id = self.add_prefix_to_id(prefix="doi", identifier=pub_id)
            pub_properties = {
                "title": pub.get("title"),
                "authors": pub.get("authors", []),
                "journal": pub.get("journal"),
                "publication_year": pub.get("publication_year"),
                "doi": pub.get("doi"),
                "pmid": pub.get("pmid"),
                "description": pub.get("description"),
                "abstract": pub.get("abstract"),
                "study_type": pub.get("study_type"),
                "organism": pub.get("organism", []),
            }
            pub_properties.update(self._get_default_properties())
            nodes.append((pub_id, "publication", pub_properties))

        return nodes


    def get_nodes(self) -> list[tuple]:
        """
        Generate publication, organism, and environmental condition nodes.

        Returns:
            List of (node_id, label, properties) tuples
        """
        logger.info("Generating nodes from omics config")
        node_list = []

        if not self.config_data:
            logger.warning("No config data loaded. Call load_config() first or provide config_file to __init__")
            return node_list

        if not hasattr(self, 'extracted_data'):
            self.download_data()

        publication = self.config_data.get('publication', {})
        if not isinstance(publication, dict):
            logger.warning(f"'publication' must be a dict, got {type(publication).__name__}. Skipping all tables.")
            return node_list

        # Create publication node
        node_list.extend(self.get_publication_nodes())

        # Collect unique treatment organisms from all analyses
        treatment_organisms = {}  # dict: taxid -> organism_name

        supp_materials = publication.get('supplementary_materials', {})
        if isinstance(supp_materials, dict):
            for table_key, table_data in supp_materials.items():
                if not isinstance(table_data, dict):
                    continue
                stat_analyses = table_data.get('statistical_analyses', [])
                if not isinstance(stat_analyses, list):
                    continue
                for stat_analysis in stat_analyses:
                    if not isinstance(stat_analysis, dict):
                        continue
                    treatment_org = stat_analysis.get('treatment_organism')
                    treatment_taxid = stat_analysis.get('treatment_taxid')
                    if treatment_org and treatment_taxid:
                        treatment_organisms[str(treatment_taxid)] = treatment_org

        # Create organism nodes for treatment organisms
        for taxid, organism_name in treatment_organisms.items():
            org_id = self.add_prefix_to_id(prefix="ncbitaxon", identifier=taxid)
            org_properties = self._get_default_properties()
            org_properties['organism_name'] = organism_name
            node_list.append((org_id, 'organism', org_properties))
            logger.info(f"Created organism node for treatment organism: {organism_name} (taxid: {taxid})")

        # Create environmental condition nodes from the config
        env_conditions = publication.get('environmental_conditions', {})
        if isinstance(env_conditions, dict):
            pub_id = self.get_publication_id()
            for env_id, env_data in env_conditions.items():
                if isinstance(env_data, dict):
                    unique_env_id = f"{pub_id}_{env_id}"
                    env_properties = self._get_default_properties()
                    for key, value in env_data.items():
                        env_properties[key] = value
                    env_properties['local_id'] = env_id
                    env_properties['publications'] = [pub_id]
                    node_list.append((unique_env_id, 'environmental_condition', env_properties))
                    logger.info(f"Created environmental condition node: {unique_env_id}")

        logger.info(f"Generated {len(node_list)} nodes (publication, organism, environmental condition)")
        return node_list

    def get_edges(self) -> list[tuple]:
        """
        Generate organism to gene expression association edges from the config file.
        
        Returns:
            List of (source_id, target_id, label, properties) tuples
        """
        logger.info("Generating organism to gene expression association edges from omics config")
        edge_list = []

        if not self.config_data:
            logger.warning("No config data loaded. Call load_config() first or provide config_file to __init__")
            return edge_list

        if not hasattr(self, 'extracted_data'):
            self.download_data()

        
        publication = self.config_data.get('publication', {})
        if not isinstance(publication, dict):
            logger.warning(f"'publication' must be a dict, got {type(publication).__name__}. Skipping all tables.")
            return edge_list

        supp_materials = publication.get('supplementary_materials', {})
        if not isinstance(supp_materials, dict):
            logger.warning(f"'supplementary_materials' must be a dict, got {type(supp_materials).__name__}. Skipping all tables.")
            return edge_list

        if not supp_materials:
            logger.warning("No supplementary materials found in config")
            return edge_list

        for table_key, table_data in supp_materials.items():
            if not isinstance(table_data, dict):
                logger.warning(f"Table '{table_key}' data must be a dict, got {type(table_data).__name__}. Skipping this table.")
                continue

            # Get file info
            filename = table_data.get('filename')
            if not filename:
                logger.warning(f"No 'filename' specified for table '{table_key}'. Skipping this table.")
                continue

            # Get statistical_analyses (required list format)
            stat_analyses = table_data.get('statistical_analyses', [])
            if not isinstance(stat_analyses, list):
                logger.warning(f"'statistical_analyses' in '{table_key}' must be a list, got {type(stat_analyses).__name__}. Skipping this table.")
                continue

            if not stat_analyses:
                logger.warning(f"No statistical analyses found in '{table_key}'")
                continue

            # Process each statistical analysis
            for idx, stat_analysis in enumerate(stat_analyses):
                if not isinstance(stat_analysis, dict):
                    logger.warning(f"Statistical analysis {idx} in '{table_key}' must be a dict, got {type(stat_analysis).__name__}. Skipping this analysis.")
                    continue

                # Load the data file
                edges_from_file = self._load_and_create_edges(
                    filename, stat_analysis
                )
                edge_list.extend(edges_from_file)

                if self.test_mode and len(edge_list) >= self.early_stopping:
                    break

        logger.info(f"Generated {len(edge_list)} organism to gene expression association edges")
        return edge_list

    def _load_and_create_edges(self, filename: str, analysis: dict) -> list[tuple]:
        """
        Load a data file and create expression association edges.

        The edge source is determined by the analysis config:
        - If environmental_treatment_condition_id is present, source is the environmental condition node
        - Otherwise, if treatment_organism/treatment_taxid are present, source is the organism node

        Args:
            filename: Path to the data file
            analysis: Dictionary with test metadata

        Returns:
            List of edge tuples
        """
        edges = []

        # Check if file exists
        if not os.path.exists(filename):
            logger.warning(f"Data file not found: {filename}")
            return edges

        try:
            # Determine edge source: environmental condition or organism
            env_condition_id = analysis.get('environmental_treatment_condition_id')

            if env_condition_id:
                # Use environmental condition as edge source
                pub_id = self.get_publication_id()
                source_id = f"{pub_id}_{env_condition_id}"
            else:
                # Fall back to organism as edge source
                treatment_organism = analysis.get('treatment_organism')
                treatment_taxid = analysis.get('treatment_taxid')

                if not treatment_organism:
                    logger.warning("No 'treatment_organism' or 'environmental_treatment_condition_id' field found in analysis. Skipping edge creation.")
                    return edges

                if not treatment_taxid:
                    logger.warning("No 'treatment_taxid' field found in analysis. Skipping edge creation.")
                    return edges

                source_id = self.add_prefix_to_id(prefix="ncbitaxon", identifier=str(treatment_taxid))

            # Load the data file, optionally skipping header rows
            skip_rows = analysis.get('skip_rows', 0)
            if skip_rows:
                df = pd.read_csv(filename, skiprows=skip_rows)
            else:
                df = pd.read_csv(filename)
            logger.info(f"Loaded {len(df)} rows from {filename}")

            # Get column mappings
            name_col = analysis.get('name_col', None)
            logfc_col = analysis.get('logfc_col', None)
            p_value_col = analysis.get('adjusted_p_value_col', None)
            pvalue_asterisk = analysis.get('pvalue_asterisk_in_logfc', False)

            # Validate columns exist
            missing_cols = []
            if name_col is None:
                logger.error(f"Required column 'name_col' not specified in analysis for file {filename}. Skipping this file.")
                return edges
            if logfc_col is None:
                logger.error(f"Required column 'logfc_col' not specified in analysis for file {filename}. Skipping this file.")
                return edges
            
            if name_col not in df.columns:
                missing_cols.append(name_col)
            if logfc_col not in df.columns:
                missing_cols.append(logfc_col)
            if missing_cols:
                logger.error(f"Required column(s) {missing_cols} not found in {filename}. Skipping this file.")
                return edges

            if not pvalue_asterisk:
                if p_value_col is None:
                    logger.warning(f"No 'adjusted_p_value_col' specified and 'pvalue_asterisk_in_logfc' is False in analysis for file {filename}. P-value will not be included in edges.")  
                elif p_value_col not in df.columns:
                    logger.warning(f"Column '{p_value_col}' not found in {filename}. P-value will not be included in edges.")
            


            # Track skipped rows for logging
            skipped_count = 0

            # Create edges for each row
            for idx, row in df.iterrows():
                if self.test_mode and idx >= self.early_stopping:
                    break

                # Get gene/protein identifier
                gene_id = row.get(name_col)
                if pd.isna(gene_id) or gene_id == '':
                    skipped_count += 1
                    continue

                # Create gene identifier with prefix
                gene_id = self.add_prefix_to_id(prefix="ncbigene", identifier=str(gene_id))

                # Skip rows with missing or non-numeric fold change values
                fc_float = None
                asterisk_significant = None
                if logfc_col in df.columns:
                    fc_val = row.get(logfc_col)
                    if pd.isna(fc_val):
                        skipped_count += 1
                        continue
                    fc_str = str(fc_val).strip()
                    if fc_str == '' or fc_str == 'NA':
                        skipped_count += 1
                        continue
                    # Handle asterisk-based significance markers in logFC column
                    if pvalue_asterisk and fc_str.endswith('*'):
                        asterisk_significant = True
                        fc_str = fc_str.rstrip('*')
                    elif pvalue_asterisk:
                        asterisk_significant = False
                    try:
                        fc_float = float(fc_str)
                    except (ValueError, TypeError):
                        skipped_count += 1
                        continue
                    if not math.isfinite(fc_float):
                        skipped_count += 1
                        continue

                # Extract edge properties
                edge_properties = self._get_default_properties()

                if fc_float is not None:
                    edge_properties['log2_fold_change'] = fc_float

                if asterisk_significant is not None:
                    # Asterisk indicates adjusted p-value < 0.1; use placeholder values
                    edge_properties['adjusted_p_value'] = 0.49 if asterisk_significant else 1.0
                elif p_value_col in df.columns and not pd.isna(row.get(p_value_col)):
                    try:
                        pval = float(row[p_value_col])
                        if math.isfinite(pval):
                            edge_properties['adjusted_p_value'] = pval
                    except (ValueError, TypeError):
                        pass

                # Determine expression direction from fold change
                if 'log2_fold_change' in edge_properties:
                    fc = edge_properties['log2_fold_change']
                    edge_properties['expression_direction'] = 'up' if fc > 0 else 'down'

                control_condition = analysis.get('control_condition')
                if control_condition:
                    edge_properties['control_condition'] = control_condition

                experimental_context = analysis.get('experimental_context')
                if experimental_context:
                    edge_properties['experimental_context'] = experimental_context

                timepoint = analysis.get('timepoint')
                if timepoint:
                    edge_properties['time_point'] = timepoint

                edge_properties['publications'] = [self.add_prefix_to_id(prefix="doi", identifier=self.get_publication_id())]

                # Create expression association edge
                # source: source_id (organism or env condition), target: gene_id
                edges.append((
                    None,
                    source_id,
                    gene_id,
                    'affects_expression_of',
                    edge_properties
                ))

            if skipped_count > 0:
                logger.info(f"Skipped {skipped_count} rows with empty or null identifiers in {filename}")

        except Exception as e:
            logger.error(f"Error processing file {filename}: {str(e)}")

        return edges


    def get_publication_id(self) -> str:
        """Get the PubMed ID from the config data, if available."""
        data = self.extracted_data 

        pub = data["publication"]
        # Use stored ID from cache (which uses DOI if available)
        pub_id = pub.get("publication_id") or pub.get("doi") or pub.get("pubmed_id") or f"pub_{pub.get('title', 'unknown')[:20]}"
        return str(pub_id)

    
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
        
    def _get_default_properties(self) -> dict:
        """
        Get default properties for nodes/edges
        """
                # generic properties for all edges for now
        properties = {
            "source": self.data_source,
            "licence": self.data_licence,
            "version": self.data_version,
        }
        return properties
        

class MultiOMICSAdapter:
    """Wrapper that reads a file listing paperconfig.yaml paths and delegates to OMICSAdapter instances."""

    def __init__(self, config_list_file: str, **kwargs):
        """
        Args:
            config_list_file: Path to a text file with one paperconfig.yaml path per line.
            **kwargs: Additional arguments passed to each OMICSAdapter (e.g. test_mode).
        """
        self.adapters = []
        with open(config_list_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    adapter = OMICSAdapter(config_file=line, **kwargs)
                    self.adapters.append(adapter)
        logger.info(f"Loaded {len(self.adapters)} config files from {config_list_file}")

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

    print("Running OMICSAdapter test")
    logger.info("Testing OMICSAdapter")

    config_dpath = 'data/Prochlorococcus/papers_and_supp/Aharonovich 2016/paperconfig.yaml'
    config_dpath = 'data/Prochlorococcus/papers_and_supp/bagby 2015/paperconfig.yaml'
    config_dpath = 'data/Prochlorococcus/papers_and_supp/biller 2016/paperconfig.yaml'
    print('pwd', os.getcwd())
    print("Config path:", config_dpath)
    print("Exists:", os.path.exists(config_dpath))

    adapter = OMICSAdapter(config_file=config_dpath,)

    adapter.download_data(cache=True)
    nodes = adapter.get_nodes()
    edges = adapter.get_edges()
    print(f"Generated {len(nodes)} nodes and {len(edges)} edges from omics config")

    import json
    with open('nodes.json', 'w') as nf:
        json.dump(nodes, nf, indent=2)

    with open('edges.json', 'w') as ef:
        json.dump(edges, ef, indent=2)

