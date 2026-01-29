from __future__ import annotations

from biocypher._logger import logger
from pydantic import BaseModel, DirectoryPath, validate_call
from typing import Literal, Union, Optional
from enum import Enum, EnumMeta, auto
import pandas as pd
import yaml
import os
from pathlib import Path
import uuid
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
    STUDY = auto()
    STATISTICAL_TEST = auto()
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


class StudyNodeField(Enum, metaclass=OmicsEnumMeta):
    """Fields for study nodes."""
    STUDY_ID = "study_id"
    TITLE = "title"
    DESCRIPTION = "description"
    ABSTRACT = "abstract"
    STUDY_TYPE = "study_type"
    ORGANISM = "organism"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class StatisticalTestNodeField(Enum, metaclass=OmicsEnumMeta):
    """Fields for statistical test nodes."""
    NAME = 'name'
    TEST_TYPE = 'test_type'
    CONTROL_CONDITION = 'control_condition'
    TREATMENT_CONDITION = 'treatment_condition'
    TIMEPOINT = 'timepoint'
    REFERENCE_TIMEPOINT = 'reference_timepoint'
    METHOD = 'method'
    ORGANISM = 'organism'

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
    study_published_in = auto()
    test_in_study = auto()
    study_has_time_series = auto()
    molecular_result_from_test = auto()
    cluster_in_study = auto()
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
    metabolomics) and creates corresponding statistical test nodes and molecular result
    edges in the knowledge graph.
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

        self._statistical_test_ids = dict()

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

        # pdf extractor for publication and study nodes
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
                OMICSEdgeType.test_in_study,
                OMICSEdgeType.molecular_result_from_test,
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

        # extract publication metadata from pdf and create publication and study nodes
        pdf_path = publication.get('papermainpdf', None)
        if pdf_path and os.path.exists(pdf_path):
            self.extracted_data = self.pdf_extractor.extract_from_pdf(pdf_path)
        else:
            logger.warning(f"PDF path not found or not provided: {pdf_path}. Skipping publication and study nodes.")
        pass

    def get_publication_nodes(self) -> list[tuple]:
        """
        Generate publication and study nodes from cached data.  

        Args:
            pdf_path: Path to the PDF file for which to generate nodes.
        Returns:
            List of (node_id, label, properties) tuples
        """       
        logger.info("Generating publication and study nodes from cached data")
        nodes = []

        data = self.extracted_data if hasattr(self, 'extracted_data') else {}

        if "publication" in data:
            pub = data["publication"]
            # Use stored ID from cache (which uses DOI if available)
            pub_id = self.get_publication_id()
            pub_id = self.add_prefix_to_id(prefix="doi", identifier=pub_id)
            pub_properties = {
                "title": pub.get("title"),
                "authors": pub.get("authors", []),
                "journal": pub.get("journal"),
                "publication_date": pub.get("publication_date"),
                "doi": pub.get("doi"),
                "pubmed_id": pub.get("pubmed_id"),
            }
            pub_properties.update(self._get_default_properties())
            nodes.append(
                (
                    pub_id,
                    "publication",
                    pub_properties
                )
            )

        if "study" in data:
            study = data["study"]
            # Use stored ID from cache - don't prefix with biolink
            study_id = study.get("study_id") 
            # No prefix for study IDs - use raw identifier
            study_properties = {
                "title": study.get("title"),
                "description": study.get("description"),
                "abstract": study.get("abstract"),
                "study_type": study.get("study_type"),
                "organism": study.get("organism", []),
            }
            study_properties.update(self._get_default_properties())
            nodes.append(

                (
                    study_id,
                    "study",
                    study_properties,
                )
            )

        return nodes

    def get_publication_study_edges(self) -> list[tuple]:
        """
        Generate study_published_in edges linking studies to publications.

        Returns:
            List of tuples (source_id, target_id, edge_type)
        """
        edges = []

        data = self.extracted_data if hasattr(self, 'extracted_data') else {}
        if "publication" in data and "study" in data:
            pub = data["publication"]
            study = data["study"]

            # Use stored IDs from cache (which use DOI if available)
            pub_id = self.get_publication_id()
            pub_id = self.add_prefix_to_id(prefix="doi", identifier=pub_id)
            study_id = study.get("study_id") 
            #study_id = self.add_prefix_to_id(prefix="biolink", identifier=study_id)
            edge_properties = self._get_default_properties()
            edges.append(
                (
                    None,
                    study_id,
                    pub_id,
                    "study_published_in",
                    edge_properties,
                )
            )

        return edges


    def get_nodes(self) -> list[tuple]:
        """
        Generate statistical test nodes from the config file.
        
        Returns:
            List of (node_id, label, properties) tuples
        """
        logger.info("Generating statistical test nodes from omics config")
        node_list = []


        if not self.config_data:
            logger.warning("No config data loaded. Call load_config() first or provide config_file to __init__")
            return node_list

        if not hasattr(self, 'extracted_data'):
            self.download_data()

        # Extract publication metadata
        publication = self.config_data.get('publication', {})
        if not isinstance(publication, dict):
            logger.warning(f"'publication' must be a dict, got {type(publication).__name__}. Skipping all tables.")
            return node_list

        # extract publication metadata from pdf and create publication and study nodes
        pub_study_nodes = self.get_publication_nodes()
        node_list.extend(pub_study_nodes)

        # Iterate through supplementary materials
        supp_materials = publication.get('supplementary_materials', {})
        if not isinstance(supp_materials, dict):
            logger.warning(f"'supplementary_materials' must be a dict, got {type(supp_materials).__name__}. Skipping all tables.")
            return node_list

        if not supp_materials:
            logger.warning("No supplementary materials found in config")
            return node_list

        for table_key, table_data in supp_materials.items():
            if not isinstance(table_data, dict):
                logger.warning(f"Table '{table_key}' data must be a dict, got {type(table_data).__name__}. Skipping this table.")
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

                node_id = self.get_statistical_test_id(stat_analysis)
                #node_id = self.add_prefix_to_id(prefix="biolink", identifier=node_id)
                properties = self._extract_test_properties(stat_analysis)
                properties.update(self._get_default_properties())
                node_list.append((node_id, 'test', properties))

                if self.test_mode and len(node_list) >= self.early_stopping:
                    break

        logger.info(f"Generated {len(node_list)} statistical test nodes")
        return node_list

    def get_edges(self) -> list[tuple]:
        """
        Generate molecular result edges from the config file.
        
        Returns:
            List of (source_id, target_id, label, properties) tuples
        """
        logger.info("Generating molecular result edges from omics config")
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

        edge_list.extend(self.get_publication_study_edges())

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

        logger.info(f"Generated {len(edge_list)} molecular result edges")
        return edge_list

    def _generate_test_id(self, analysis: dict) -> str:
        """
        Generate a unique ID for a statistical test.
        
        Args:
            analysis: Dictionary with test metadata
            
        Returns:
            A unique test ID
        """
        pub_id = self.get_publication_id()
        # Try to create a meaningful ID from test name
        test_id = analysis.get('id', '')
        if test_id:
            # Create a slug from the name
            slug = test_id.lower().replace(' ', '_')[:50]
            test_id = f"test_{slug}"
        else:
            logger.warning("No 'id' field found in analysis; generating random test ID")
            test_id = f"test_{uuid.uuid4().hex}"

        test_id = f"{pub_id}_{test_id}"
        return test_id

    def _extract_test_properties(self, analysis: dict) -> dict:
        """
        Extract properties for a statistical test node.
        
        Args:
            analysis: Dictionary with test metadata
            
        Returns:
            Dictionary of node properties
        """
        properties = {}

        # Map analysis fields to node properties
        field_mapping = {
            'name': StatisticalTestNodeField.NAME.value,
            'test_type': StatisticalTestNodeField.TEST_TYPE.value,
            'control_condition': StatisticalTestNodeField.CONTROL_CONDITION.value,
            'treatment_condition': StatisticalTestNodeField.TREATMENT_CONDITION.value,
            'timepoint': StatisticalTestNodeField.TIMEPOINT.value,
            'reference_timepoint': StatisticalTestNodeField.REFERENCE_TIMEPOINT.value,
            'method': StatisticalTestNodeField.METHOD.value,
            'organism': StatisticalTestNodeField.ORGANISM.value,
        }

        for key, prop_name in field_mapping.items():
            if key in analysis:
                value = analysis[key]
                if value is not None:
                    properties[prop_name] = value

        return properties

    def _load_and_create_edges(self, filename: str, analysis: dict) -> list[tuple]:
        """
        Load a data file and create edges for molecular results.
        
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

            # Generate test ID
            test_id = self.get_statistical_test_id(analysis)
            #test_id = self.add_prefix_to_id(prefix="biolink", identifier=test_id)

            # generate statistical test to study edge
            study = self.extracted_data.get('study', {})
            study_id = study.get('study_id')
            #study_id = self.add_prefix_to_id(prefix="biolink", identifier=study_id)
            if study_id:
                edge_properties = self._get_default_properties()
                edges.append((
                    None,
                    test_id,
                    study_id,
                    'test_in_study',
                    edge_properties,
                ))

            # Load the data file
            df = pd.read_csv(filename)
            logger.info(f"Loaded {len(df)} rows from {filename}")

            # Get column mappings
            name_col = analysis.get('name_col', 'Synonym')
            logfc_col = analysis.get('logfc_col', 'log2_fold_change')
            p_value_col = analysis.get('adjusted_p_value_col', 'adjusted_p_value')

            # Validate columns exist
            missing_cols = []
            if name_col not in df.columns:
                missing_cols.append(name_col)
            if logfc_col not in df.columns:
                logger.warning(f"Column '{logfc_col}' not found in {filename}. Fold change will not be included in edges.")
            if p_value_col not in df.columns:
                logger.warning(f"Column '{p_value_col}' not found in {filename}. P-value will not be included in edges.")
            
            if missing_cols:
                logger.error(f"Required column(s) {missing_cols} not found in {filename}. Skipping this file.")
                return edges


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

                # Extract edge properties
                edge_properties = self._get_default_properties()

                if logfc_col in df.columns and not pd.isna(row.get(logfc_col)):
                    edge_properties['log2_fold_change'] = float(row[logfc_col])

                if p_value_col in df.columns and not pd.isna(row.get(p_value_col)):
                    edge_properties['adjusted_p_value'] = float(row[p_value_col])

                # Determine direction from fold change
                if 'log2_fold_change' in edge_properties:
                    fc = edge_properties['log2_fold_change']
                    edge_properties['direction'] = 'up' if fc > 0 else 'down'
                    edge_properties['predicate'] = 'has_increased_amount' if fc > 0 else 'has_decreased_amount'

                edge_properties['publications'] = [self.add_prefix_to_id(prefix="doi", identifier=self.get_publication_id())]

                # Create edge
                edges.append((
                    None,
                    gene_id,
                    test_id,
                    'molecular_result_from_test',
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

    
    def get_statistical_test_id(self, analysis: dict) -> str:
        """Get the statistical test ID from the analysis data."""
        if analysis.get('id') and analysis['id'] in self._statistical_test_ids:
            return self._statistical_test_ids[analysis['id']]
        test_id = self._generate_test_id(analysis)
        self._statistical_test_ids[analysis['id']] = test_id
        return test_id
        return self._generate_test_id(analysis)


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
        

if __name__ == "__main__":

    print("Running OMICSAdapter test")
    logger.info("Testing OMICSAdapter")

    config_dpath = 'data/Prochlorococcus/papers_and_supp/Aharonovich 2016/paperconfig.yaml'
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

