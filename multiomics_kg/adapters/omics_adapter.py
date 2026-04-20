from __future__ import annotations

from biocypher._logger import logger
from pydantic import BaseModel, DirectoryPath, validate_call
from typing import Literal, Union, Optional
from enum import Enum, EnumMeta, auto
import math
import pandas as pd
import os
from pathlib import Path
from bioregistry import normalize_curie

from multiomics_kg.download.resolve_paper_ids import get_resolved_path
from multiomics_kg.utils.paperconfig_utils import (
    load_paperconfig, get_experiments, get_experiment_for_analysis, iter_analyses,
    parse_timepoint_hours,
)

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
    EXPERIMENT = auto()


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


class OMICSEdgeType(Enum, metaclass=OmicsEnumMeta):
    """Edge types for omics data."""
    changes_expression_of = auto()
    has_experiment = auto()
    tests_coculture_with = auto()


class OMICSModel(BaseModel):
    """Pydantic model for OMICS adapter configuration."""
    edge_types: Union[list[OMICSEdgeType], None] = None
    test_mode: bool = False
    export_csv: bool = False
    output_dir: DirectoryPath | None = None
    add_prefix: bool = True
    organism: int | Literal["*"] | None = None
    significance_mode: Literal["all", "significant_only"] = "all"
    default_pvalue_threshold: float | None = None
    default_logfc_threshold: float | None = None


def _convert_fold_change(fc_float: float, fold_change_type: str | None) -> float | None:
    """Convert fold-change value to log2 scale.

    Args:
        fc_float: The raw fold-change value from the CSV.
        fold_change_type: 'log2' (default), 'linear', or None (treated as log2).

    Returns:
        log2-scale fold-change, or None if conversion is invalid (e.g., linear FC <= 0).
    """
    if fold_change_type is None or fold_change_type == "log2":
        return fc_float
    if fold_change_type == "linear":
        if fc_float <= 0:
            return None
        return math.log2(fc_float)
    return fc_float


def _validate_fc_range(
    fc_values: list[float],
    fold_change_type: str | None,
    analysis_id: str,
    table_scope: str | None = None,
) -> None:
    """Warn if fold-change values look inconsistent with declared type.

    Heuristics:
    - linear FC must be all positive (negative linear FC is invalid)
    - log2 FC that is all positive AND all > 1.0 may actually be linear FC
      (unless table_scope indicates only upregulated genes)
    """
    if not fc_values:
        return
    fc_type = fold_change_type or "log2"
    has_negative = any(v < 0 for v in fc_values)
    all_positive = not has_negative
    min_val = min(fc_values)
    max_val = max(fc_values)

    if fc_type == "linear" and has_negative:
        logger.warning(
            "%s: fold_change_type is 'linear' but data contains negative values "
            "(min=%.3f). Negative values are invalid for linear FC — check if "
            "these are actually log2 values.",
            analysis_id, min_val,
        )
    elif fc_type == "log2" and all_positive and min_val > 1.0:
        # All values > 1.0 is suspicious for log2 — could be linear FC
        # But skip warning if table only has upregulated genes
        upregulated_scopes = {"significant_only"}
        if table_scope not in upregulated_scopes:
            logger.warning(
                "%s: fold_change_type is 'log2' but all values are > 1.0 "
                "(range %.3f–%.3f). These may be linear fold-changes. "
                "Consider setting fold_change_type: linear.",
                analysis_id, min_val, max_val,
            )


class OMICSAdapter:
    """
    Adapter for creating omics-related nodes and edges from publication config files.

    Parses YAML paperconfig files that specify statistical analyses (RNA-seq, proteomics,
    metabolomics) and creates publication, organism, and environmental condition nodes,
    plus expression edges linking causes to genes.
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
        significance_mode: Literal["all", "significant_only"] = "all",
        default_pvalue_threshold: float | None = 0.05,
        default_logfc_threshold: float | None = 1.0,
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
            significance_mode: "all" to create all edges with significant flag,
                "significant_only" to skip non-significant edges
            default_pvalue_threshold: Default adjusted p-value threshold for significance
            default_logfc_threshold: Default |log2FC| threshold for significance
        """
        model = OMICSModel(
            edge_types=edge_types,
            test_mode=test_mode,
            export_csv=export_csv,
            output_dir=output_dir,
            add_prefix=add_prefix,
            organism=organism,
            significance_mode=significance_mode,
            default_pvalue_threshold=default_pvalue_threshold,
            default_logfc_threshold=default_logfc_threshold,
        ).model_dump()

        self.export_csv = model["export_csv"]
        self.output_dir = model["output_dir"]
        self.add_prefix = model["add_prefix"]
        self.test_mode = model["test_mode"]
        self.significance_mode = model["significance_mode"]
        self.default_pvalue_threshold = model["default_pvalue_threshold"]
        self.default_logfc_threshold = model["default_logfc_threshold"]

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
                OMICSEdgeType.changes_expression_of,
                OMICSEdgeType.has_experiment,
                OMICSEdgeType.tests_coculture_with,
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

        self.config_data = load_paperconfig(Path(config_file))
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

        if not publication:
            logger.info("No publication block in config. Skipping publication nodes.")
            return

        # extract publication metadata from pdf
        pdf_path = publication.get('papermainpdf', None)
        if not pdf_path or not os.path.exists(pdf_path):
            raise FileNotFoundError(
                f"PDF path missing or not found for "
                f"{publication.get('papername', '?')}: {pdf_path!r}"
            )
        self.extracted_data = self.pdf_extractor.extract_from_pdf(pdf_path)
        if not self.extracted_data or "publication" not in self.extracted_data:
            raise RuntimeError(
                f"PDF extraction produced no publication block for {pdf_path}"
            )

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
                "title": self.clean_text(pub.get("title")),
                "authors": self.clean_text(pub.get("authors", [])),
                "journal": self.clean_text(pub.get("journal")),
                "publication_year": pub.get("publication_year"),
                "doi": self._get_override_doi() or pub.get("doi"),
                "pmid": pub.get("pmid"),
                "description": self.clean_text(pub.get("description")),
                "abstract": self.clean_text(pub.get("abstract")),
                "study_type": self.clean_text(pub.get("study_type")),
            }
            pub_properties.update(self._get_default_properties())
            nodes.append((pub_id, "publication", pub_properties))

        return nodes


    def get_nodes(self) -> list[tuple]:
        """
        Generate publication and experiment nodes.

        Organism nodes are created by the CyanorakNcbi adapter (single source
        of truth).  This adapter creates publication and experiment nodes.

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

        # Create experiment nodes from the experiments block
        experiments = get_experiments(self.config_data)
        pub_id = self.get_publication_id()
        for exp_key, exp in experiments.items():
            if not isinstance(exp, dict):
                continue
            experiment_id = f"{pub_id}_{exp_key}"

            # Count distinct timepoints referencing this experiment to determine is_time_course
            timepoints = set()
            for _, _, a in iter_analyses(self.config_data):
                if a.get("experiment") == exp_key and a.get("timepoint"):
                    timepoints.add(a.get("timepoint"))

            exp_props = {
                    "name": self.clean_text(exp.get("name", "")),
                    "organism_name": self.clean_text(exp.get("organism", "")),
                    "compartment": self.clean_text(exp.get("compartment", "whole_cell")),
                    "treatment_type": self._normalize_list_field(exp, "treatment_type"),
                    "treatment": self.clean_text(exp.get("treatment_condition", "")),
                    "control": self.clean_text(exp.get("control_condition", "")),
                    "experimental_context": self.clean_text(exp.get("experimental_context", "")),
                    "omics_type": self.clean_text(exp.get("omics_type", "")),
                    "statistical_test": self.clean_text(exp.get("test_type", "")),
                    "is_time_course": "true" if len(timepoints) > 1 else "false",
                    "medium": self.clean_text(exp.get("medium", "")),
                    "temperature": self.clean_text(exp.get("temperature", "")),
                    "light_condition": self.clean_text(exp.get("light_condition", "")),
                    "light_intensity": self.clean_text(exp.get("light_intensity", "")),
                    "table_scope": self.clean_text(exp.get("table_scope", "")),
                    "table_scope_detail": self.clean_text(exp.get("table_scope_detail", "")),
                    "background_factors": self._normalize_list_field(exp, "background_factors"),
            }
            partner = exp.get("treatment_organism", "")
            if partner:
                exp_props["coculture_partner"] = self.clean_text(partner)
            node_list.append((experiment_id, "experiment", exp_props))
            logger.info(f"Created experiment node: {experiment_id}")

        logger.info(f"Generated {len(node_list)} nodes (publication, experiment)")
        return node_list

    def get_edges(self) -> list[tuple]:
        """
        Generate experiment-to-gene expression edges, has_experiment edges,
        and tests_coculture_with edges from the config file.

        Returns:
            List of (edge_id, source_id, target_id, label, properties) tuples
        """
        logger.info("Generating expression edges from omics config")
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

        # Validate that all statistical analyses have a unique 'id' field within this publication
        seen_analysis_ids = set()
        for table_key, table_data in supp_materials.items():
            if not isinstance(table_data, dict):
                continue
            for idx, sa in enumerate(table_data.get('statistical_analyses', [])):
                if not isinstance(sa, dict):
                    continue
                analysis_id = sa.get('id')
                if not analysis_id:
                    raise ValueError(
                        f"Statistical analysis {idx} in '{table_key}' is missing required 'id' field."
                    )
                if analysis_id in seen_analysis_ids:
                    raise ValueError(
                        f"Duplicate statistical analysis id '{analysis_id}' in '{table_key}'. "
                        f"Each analysis must have a unique 'id' within the publication."
                    )
                seen_analysis_ids.add(analysis_id)

        pub_id_raw = self.get_publication_id()
        pub_id = self.add_prefix_to_id(prefix="doi", identifier=pub_id_raw)
        experiments = get_experiments(self.config_data)

        # ── Compute time_point_order per experiment ──
        # Collect analyses per experiment, sort by timepoint_hours (null last), then by id
        exp_analyses = {}  # exp_key -> [(analysis_dict, table_key, table_data)]
        for table_key, table_data in supp_materials.items():
            if not isinstance(table_data, dict):
                continue
            stat_analyses = table_data.get('statistical_analyses', [])
            if not isinstance(stat_analyses, list):
                continue
            for sa in stat_analyses:
                if not isinstance(sa, dict):
                    continue
                exp_key = sa.get('experiment')
                if exp_key:
                    exp_analyses.setdefault(exp_key, []).append((sa, table_key, table_data))

        # Build time_point_order mapping: analysis_id -> order
        analysis_tp_order = {}
        for exp_key, analyses_list in exp_analyses.items():
            # Sort by timepoint_hours (None sorts last), then by analysis id for stability
            def sort_key(item):
                sa = item[0]
                th = sa.get('timepoint_hours')
                if th is None:
                    # Try parsing from timepoint string
                    th = parse_timepoint_hours(sa.get('timepoint'))
                return (th if th is not None else float('inf'), sa.get('id', ''))
            sorted_analyses = sorted(analyses_list, key=sort_key)
            for order, (sa, _, _) in enumerate(sorted_analyses, start=1):
                analysis_tp_order[sa.get('id')] = order

        # ── Emit structural edges: has_experiment + tests_coculture_with ──
        for exp_key, exp in experiments.items():
            if not isinstance(exp, dict):
                continue
            experiment_id = f"{pub_id_raw}_{exp_key}"

            # has_experiment edge: Publication → Experiment
            edge_list.append((
                f"{pub_id_raw}_has_exp_{exp_key}",
                pub_id,
                experiment_id,
                'has_experiment',
                self._get_default_properties()
            ))

            # tests_coculture_with edge: Experiment → OrganismTaxon (coculture only)
            treatment_organism = exp.get('treatment_organism')
            if treatment_organism:
                treatment_accession = exp.get('treatment_assembly_accession')
                treatment_taxid = exp.get('treatment_taxid')
                if treatment_accession:
                    organism_id = self.add_prefix_to_id(
                        prefix="insdc.gcf", identifier=treatment_accession
                    )
                elif treatment_taxid:
                    organism_id = self.add_prefix_to_id(
                        prefix="ncbitaxon", identifier=str(treatment_taxid)
                    )
                else:
                    organism_id = None

                if organism_id:
                    edge_list.append((
                        f"{pub_id_raw}_coculture_{exp_key}",
                        experiment_id,
                        organism_id,
                        'tests_coculture_with',
                        self._get_default_properties()
                    ))

        # ── Emit expression edges: Experiment → Gene ──
        for table_key, table_data in supp_materials.items():
            if not isinstance(table_data, dict):
                logger.warning(f"Table '{table_key}' data must be a dict, got {type(table_data).__name__}. Skipping this table.")
                continue

            # Skip non-csv types (gene_clusters, id_translation, annotation_gff)
            # — these are handled by their own adapters/scripts.
            table_type = table_data.get("type", "csv")
            if table_type != "csv":
                continue

            filename = table_data.get('filename')
            if not filename:
                logger.warning(f"No 'filename' specified for table '{table_key}'. Skipping this table.")
                continue

            stat_analyses = table_data.get('statistical_analyses', [])
            if not isinstance(stat_analyses, list):
                logger.warning(f"'statistical_analyses' in '{table_key}' must be a list, got {type(stat_analyses).__name__}. Skipping this table.")
                continue

            if not stat_analyses:
                logger.warning(f"No statistical analyses found in '{table_key}'")
                continue

            # Get file-level settings
            skip_rows = table_data.get('skip_rows', 0)
            sep = table_data.get('sep', ',')

            for idx, stat_analysis in enumerate(stat_analyses):
                if not isinstance(stat_analysis, dict):
                    logger.warning(f"Statistical analysis {idx} in '{table_key}' must be a dict, got {type(stat_analysis).__name__}. Skipping this analysis.")
                    continue

                # Pass file-level settings into analysis if not already set
                if skip_rows and 'skip_rows' not in stat_analysis:
                    stat_analysis['skip_rows'] = skip_rows
                if sep and 'sep' not in stat_analysis:
                    stat_analysis['sep'] = sep
                # Pass table-level table_scope into analysis if not already set
                table_scope = table_data.get('table_scope')
                if table_scope and 'table_scope' not in stat_analysis:
                    stat_analysis['table_scope'] = table_scope

                # Determine experiment-level info
                exp_key = stat_analysis.get('experiment')
                if not exp_key or exp_key not in experiments:
                    logger.warning(f"Analysis '{stat_analysis.get('id')}' has no valid experiment reference. Skipping.")
                    continue

                experiment_id = f"{pub_id_raw}_{exp_key}"
                analysis_id = stat_analysis.get('id', '')
                tp_order = analysis_tp_order.get(analysis_id, 1)

                # Timepoint info
                timepoint = stat_analysis.get('timepoint')
                timepoint_hours = stat_analysis.get('timepoint_hours')
                if timepoint_hours is None:
                    timepoint_hours = parse_timepoint_hours(timepoint)
                growth_phase = stat_analysis.get('growth_phase')

                edges_from_file = self._load_and_create_edges(
                    filename, stat_analysis,
                    experiment_id=experiment_id,
                    time_point=timepoint,
                    time_point_order=tp_order,
                    time_point_hours=timepoint_hours,
                    growth_phase=growth_phase,
                )
                edge_list.extend(edges_from_file)

                if self.test_mode and len(edge_list) >= self.early_stopping:
                    break

        logger.info(f"Generated {len(edge_list)} edges (expression, has_experiment, tests_coculture_with)")
        return edge_list

    def _check_significance(
        self,
        fc_value: float | None,
        pvalue: float | None,
        asterisk_significant: bool | None,
        analysis: dict,
    ) -> str:
        """
        Determine whether a row is statistically significant.

        Priority order:
        1. prefiltered: true → always "significant"
        2. pvalue_asterisk_in_logfc → use asterisk_significant
        3. Thresholds (config overrides constructor defaults)
        4. "unknown" if no info available

        Returns:
            "significant", "not significant", or "unknown".
        """
        if analysis.get('prefiltered'):
            return 'significant'

        if asterisk_significant is not None:
            return 'significant' if asterisk_significant else 'not significant'

        pval_thresh = analysis.get('pvalue_threshold') or self.default_pvalue_threshold
        logfc_thresh = analysis.get('logfc_threshold') or self.default_logfc_threshold

        if pval_thresh is None and logfc_thresh is None:
            return 'unknown'

        sig = True
        if logfc_thresh is not None and fc_value is not None:
            sig = sig and (abs(fc_value) >= logfc_thresh)
        if pval_thresh is not None and pvalue is not None:
            sig = sig and (pvalue <= pval_thresh)

        return 'significant' if sig else 'not significant'

    def _load_and_create_edges(
        self,
        filename: str,
        analysis: dict,
        experiment_id: str = "",
        time_point: str | None = None,
        time_point_order: int = 1,
        time_point_hours: float | None = None,
        growth_phase: str | None = None,
    ) -> list[tuple]:
        """
        Load a data file and create expression association edges.

        Source is the experiment node; target is the gene.

        Args:
            filename: Path to the data file
            analysis: Dictionary with test metadata
            experiment_id: The experiment node ID (source of edges)
            time_point: Timepoint label (e.g., "3h", "day 18")
            time_point_order: 1-indexed order within the experiment
            time_point_hours: Numeric hours, or None if unparseable
            growth_phase: Growth phase label (e.g., "nutrient_limited"), or None

        Returns:
            List of edge tuples
        """
        edges = []

        # Check if file exists
        if not os.path.exists(filename):
            logger.warning(f"Data file not found: {filename}")
            return edges

        # Probe for pre-resolved CSV (written by resolve_paper_ids.py)
        _resolved_path = get_resolved_path(filename)
        _use_resolved = _resolved_path.exists()

        try:
            # Load the data file
            # _resolved.csv is always a clean comma-separated file with no extra header rows,
            # so use plain pd.read_csv() regardless of the original sep/skip_rows settings.
            if _use_resolved:
                df = pd.read_csv(str(_resolved_path))
            else:
                skip_rows = analysis.get('skip_rows', 0)
                sep = analysis.get('sep', ',')
                if skip_rows:
                    df = pd.read_csv(filename, sep=sep, skiprows=skip_rows)
                else:
                    df = pd.read_csv(filename, sep=sep)

            _resolved_col = 'resolved_locus_tag'
            use_locus_tag_col = _use_resolved and _resolved_col in df.columns
            if _use_resolved:
                if use_locus_tag_col:
                    logger.info(f"Using pre-resolved CSV: {_resolved_path.name} ({len(df)} rows)")
                else:
                    logger.warning(f"Pre-resolved CSV {_resolved_path.name} has no {_resolved_col} column; falling back to {filename}")
            else:
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
            raw_fc_values = []  # collect raw FC values for range validation
            pub_id_raw = self.get_publication_id()
            analysis_id = analysis.get('id', '')

            # Create edges for each row
            for idx, row in df.iterrows():
                if self.test_mode and idx >= self.early_stopping:
                    break

                # Get gene/protein identifier
                # Pre-resolved CSV: use the pre-computed locus_tag column;
                # rows with NaN locus_tag (unresolved) are skipped here.
                gene_id = row.get(_resolved_col) if use_locus_tag_col else row.get(name_col)
                if pd.isna(gene_id) or gene_id == '':
                    skipped_count += 1
                    continue

                # Strip whitespace and asterisks from gene IDs
                gene_id = str(gene_id).strip().strip('*').strip()
                if gene_id == '':
                    skipped_count += 1
                    continue

                # Create gene identifier with prefix
                gene_id = self.add_prefix_to_id(prefix="ncbigene", identifier=gene_id)

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
                    # Supports trailing "1.1 *", leading "* 1.1", or bare "1.1*"
                    if pvalue_asterisk:
                        has_asterisk = fc_str.startswith('*') or fc_str.endswith('*')
                        asterisk_significant = has_asterisk
                        fc_str = fc_str.strip('*').strip()
                    else:
                        # Always strip stray whitespace and stars from logFC values
                        fc_str = fc_str.strip('*').strip()
                    try:
                        fc_float = float(fc_str)
                    except (ValueError, TypeError):
                        skipped_count += 1
                        continue
                    if not math.isfinite(fc_float):
                        skipped_count += 1
                        continue
                    raw_fc_values.append(fc_float)

                # Extract edge properties — only per-gene/per-timepoint fields
                edge_properties = {}

                if fc_float is not None:
                    fold_change_type = analysis.get('fold_change_type', None)
                    converted = _convert_fold_change(fc_float, fold_change_type)
                    if converted is None:
                        logger.debug(
                            "Skipping row with invalid linear FC %.4f in %s",
                            fc_float, filename
                        )
                        skipped_count += 1
                        continue
                    edge_properties['log2_fold_change'] = converted

                pval = None
                if asterisk_significant is not None:
                    # Asterisk encodes significance at pvalue_threshold; use threshold as placeholder
                    pval_thresh = analysis.get('pvalue_threshold') or self.default_pvalue_threshold or 0.05
                    pval = pval_thresh if asterisk_significant else 1.0
                    edge_properties['adjusted_p_value'] = pval
                elif p_value_col and p_value_col in df.columns and not pd.isna(row.get(p_value_col)):
                    try:
                        pval = float(row[p_value_col])
                        if math.isfinite(pval):
                            edge_properties['adjusted_p_value'] = pval
                        else:
                            pval = None
                    except (ValueError, TypeError):
                        pass

                # Check significance and optionally filter
                significant = self._check_significance(fc_float, pval, asterisk_significant, analysis)
                edge_properties['significant'] = significant
                if self.significance_mode == "significant_only" and significant == 'not significant':
                    skipped_count += 1
                    continue

                # Determine expression direction from fold change
                if 'log2_fold_change' in edge_properties:
                    fc = edge_properties['log2_fold_change']
                    edge_properties['expression_direction'] = 'up' if fc > 0 else 'down'

                # Timepoint properties
                if time_point:
                    edge_properties['time_point'] = self.clean_text(time_point)
                edge_properties['time_point_order'] = time_point_order
                if time_point_hours is not None:
                    edge_properties['time_point_hours'] = time_point_hours
                if growth_phase:
                    edge_properties['growth_phase'] = self.clean_text(growth_phase)

                # Create expression association edge
                # source: experiment_id, target: gene_id
                edge_id = f"{pub_id_raw}_{analysis_id}_{gene_id}"
                edges.append((
                    edge_id,
                    experiment_id,
                    gene_id,
                    'changes_expression_of',
                    edge_properties
                ))

            if skipped_count > 0:
                src_label = _resolved_path.name if _use_resolved else filename
                logger.info(f"Skipped {skipped_count} rows with empty or null identifiers in {src_label}")

            # Validate that FC values are consistent with declared type
            _validate_fc_range(
                raw_fc_values,
                analysis.get('fold_change_type', None),
                analysis_id,
                table_scope=analysis.get('table_scope'),
            )

        except Exception as e:
            logger.error(f"Error processing file {filename}: {str(e)}")

        return edges


    def _get_override_doi(self) -> str | None:
        """Read optional doi from paperconfig publication block."""
        pub = self.config_data.get("publication", {}) or {}
        val = pub.get("doi")
        if not isinstance(val, str):
            return None
        val = val.strip()
        return val or None

    def get_publication_id(self) -> str:
        """Get the publication ID from config doi override, PDF extraction, or fallback."""
        override = self._get_override_doi()
        if hasattr(self, 'extracted_data') and "publication" in self.extracted_data:
            pub = self.extracted_data["publication"]
            extracted_doi = pub.get("doi")
            if override:
                if extracted_doi and extracted_doi != override:
                    logger.warning(
                        f"Config doi '{override}' disagrees with PDF-extracted doi "
                        f"'{extracted_doi}'; using config value."
                    )
                return override
            pub_id = pub.get("publication_id") or extracted_doi or pub.get("pubmed_id") or f"pub_{pub.get('title', 'unknown')[:20]}"
            return str(pub_id)
        if override:
            return override
        pub = self.config_data.get('publication', {})
        pub_id = pub.get("pubmed_id") or pub.get("papername") or "unknown"
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

    def _normalize_list_field(self, exp: dict, field: str) -> list[str]:
        """Normalize a paperconfig field to a list of cleaned strings."""
        val = exp.get(field, [])
        if isinstance(val, str):
            val = [val] if val else []
        return [self.clean_text(v) for v in val]

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

    def __init__(self, config_list_file: str | list[str], **kwargs):
        """
        Args:
            config_list_file: Path to a text file with one paperconfig.yaml path per line,
                or a list of such paths (for multi-organism support).
            **kwargs: Additional arguments passed to each OMICSAdapter (e.g. test_mode).
        """
        self.adapters = []
        list_files = config_list_file if isinstance(config_list_file, list) else [config_list_file]
        for lf in list_files:
            with open(lf, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        adapter = OMICSAdapter(config_file=line, **kwargs)
                        self.adapters.append(adapter)
            logger.info(f"Loaded configs from {lf}")
        logger.info(f"Total: {len(self.adapters)} paperconfig files from {len(list_files)} list file(s)")

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
    config_dpath = 'data/Prochlorococcus/papers_and_supp/biller 2018/paperconfig.yaml'
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

