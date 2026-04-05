"""
Gene Cluster Adapter

Reads gene_clusters entries from paperconfig.yaml files and emits:
- ClusteringAnalysis nodes (one per gene_clusters entry)
- GeneCluster nodes (one per unique cluster value in the CSV)
- Gene_in_gene_cluster membership edges
- Publication_has_clustering_analysis edges
- Clustering_analysis_has_gene_cluster edges
- Clusteringanalysis_belongs_to_organism edges
- Experiment_has_clustering_analysis edges
"""
import json
import logging
import re
from pathlib import Path

import pandas as pd

from multiomics_kg.download.resolve_paper_ids import get_resolved_path
from multiomics_kg.utils.paperconfig_utils import (
    load_paperconfig,
    load_all_paperconfigs,
    get_paper_name,
    iter_cluster_tables,
)

logger = logging.getLogger(__name__)

# Default cache path for PDF extraction (same as pdf_publication_extraction.py)
_DEFAULT_PDF_CACHE = Path(__file__).parent.parent.parent / "cache" / "pdf_extraction_cache.json"


def _clean_str(value) -> str:
    """Sanitize string for BioCypher CSV output."""
    if not isinstance(value, str):
        return str(value) if value is not None else ""
    return value.replace("'", "^").replace("|", ",")


def _make_cluster_id(doi: str, paper_name: str, entry_key: str, cluster_key: str) -> str:
    """Build cluster node ID: cluster:{doi_short}:{entry_key}:{cluster_key}.

    Uses DOI short form when available, falls back to paper_name slug.
    """
    if doi:
        doi_short = doi.rstrip("/").rsplit("/", 1)[-1]
    else:
        doi_short = re.sub(r"[^a-z0-9]+", "_", paper_name.lower()).strip("_")
    return f"cluster:{doi_short}:{entry_key}:{cluster_key}"


def _make_analysis_id(doi: str, paper_name: str, entry_key: str) -> str:
    """Build clustering analysis node ID."""
    if doi:
        doi_short = doi.rstrip("/").rsplit("/", 1)[-1]
    else:
        doi_short = re.sub(r"[^a-z0-9]+", "_", paper_name.lower()).strip("_")
    return f"clustering_analysis:{doi_short}:{entry_key}"


def _load_extraction_json(paperconfig_dir: Path, entry_key: str) -> dict:
    """Load extraction data. Tries new versioned structure first, falls back to legacy."""
    from multiomics_kg.extraction.cluster.run_manager import RunManager

    cache_dir = paperconfig_dir / ".extraction_cache"
    entry_dir = cache_dir / entry_key
    current_link = entry_dir / "current"
    if current_link.is_symlink() or current_link.is_dir():
        run_dir = current_link.resolve() if current_link.is_symlink() else current_link
        if run_dir.exists():
            result = {}
            stage_key_map = {1: "stage1_merged", 2: "stage2_results",
                             3: "stage3_validation", 4: "stage4_review"}
            for stage, filename in RunManager.STAGE_FILES.items():
                stage_path = run_dir / filename
                if stage_path.exists():
                    with open(stage_path) as f:
                        result[stage_key_map[stage]] = json.load(f)
            return result

    # Legacy fallback
    json_path = paperconfig_dir / f"cluster_extraction_{entry_key}.json"
    if not json_path.exists():
        return {}
    try:
        with open(json_path) as f:
            return json.load(f)
    except Exception:
        logger.warning("Failed to load extraction JSON: %s", json_path)
        return {}


def _get_extraction_cluster_data(extraction: dict, cluster_key: str) -> dict:
    """Get per-cluster data, respecting validation verdict and review status."""
    stage2 = extraction.get("stage2_results", {})
    stage3 = extraction.get("stage3_validation", {})
    stage4 = extraction.get("stage4_review", {})

    cluster_data = stage2.get(str(cluster_key), {})
    verdict = stage3.get(str(cluster_key), {}).get("verdict", "")
    review = stage4.get(str(cluster_key), {})
    review_status = review.get("status", "")

    # If review exists, it takes precedence over verdict
    if review_status:
        if review_status in ("approve", "edit"):
            result = dict(cluster_data)
            for field, value in review.get("edited_fields", {}).items():
                result[field] = value
            return result
        else:
            return {}  # reject, flag-issue, stale

    # No review — fall back to verdict check (legacy behavior)
    if verdict != "pass":
        return {}
    return cluster_data


def _resolve_csv_path(csv_path: str) -> tuple[Path, bool]:
    """Probe for pre-resolved CSV (written by resolve_paper_ids.py).

    Returns (path_to_use, is_resolved).
    """
    p = Path(csv_path)
    resolved = get_resolved_path(p)
    if resolved.exists():
        return resolved, True
    return p, False


def _load_pdf_cache(cache_path: Path = _DEFAULT_PDF_CACHE) -> dict:
    """Load the PDF extraction cache (DOI, title, etc.)."""
    if cache_path.exists():
        try:
            with open(cache_path) as f:
                return json.load(f)
        except Exception:
            pass
    return {}


def _cluster_val_to_str(val) -> str:
    """Convert a cluster column value to a clean string.

    Pandas reads integer columns with NaN as float (e.g. 5 → 5.0).
    This converts ``5.0`` → ``"5"`` while leaving genuine strings untouched.
    """
    if pd.isna(val):
        return ""
    if isinstance(val, float) and val == int(val):
        return str(int(val))
    return str(val).strip()


def _read_cluster_csv(csv_path: Path, table: dict) -> pd.DataFrame:
    """Read a cluster CSV, honouring optional ``skip_rows`` from paperconfig.

    ``skip_rows`` (int) tells pandas to skip that many leading rows before the
    real header row.  Useful when the CSV has multi-row headers (e.g. merged
    cells exported from Excel).
    """
    skip = table.get("skip_rows", 0)
    if skip:
        return pd.read_csv(csv_path, skiprows=skip)
    return pd.read_csv(csv_path)


class ClusterAdapter:
    """Adapter for one paperconfig's gene_clusters entries."""

    def __init__(self, config_file: str, test_mode: bool = False):
        self.config_file = config_file
        self.test_mode = test_mode
        self.config = load_paperconfig(Path(config_file))
        self.paper_name = get_paper_name(self.config, fallback_path=Path(config_file))
        self.doi = self._extract_doi()
        self._cluster_tables = list(iter_cluster_tables(self.config))
        # Populated by MultiClusterAdapter
        self._organism_lookup: dict[str, str] = {}
        # Directory containing paperconfig (for extraction JSON lookup)
        self._paperconfig_dir = Path(config_file).parent

    def _extract_doi(self) -> str:
        """Get DOI from paperconfig or PDF extraction cache."""
        pub = self.config.get("publication", {})
        # Direct doi field in config (rare but possible)
        doi = pub.get("doi", "")
        if doi:
            return doi
        # Look up from PDF extraction cache (same source as omics_adapter)
        pdf_path = pub.get("papermainpdf", "")
        if pdf_path:
            cache = _load_pdf_cache()
            cached = cache.get(pdf_path, {})
            doi = cached.get("publication", {}).get("doi", "")
        return doi or ""

    def get_nodes(self) -> list[tuple]:
        """Emit ClusteringAnalysis and GeneCluster nodes."""
        nodes = []
        for entry_key, table in self._cluster_tables:
            csv_path, _use_resolved = _resolve_csv_path(table["filename"])
            if not csv_path.exists():
                logger.warning(f"Cluster CSV not found: {csv_path}")
                continue

            # Resolved CSVs are clean single-header files; skip_rows doesn't apply
            if _use_resolved:
                df = pd.read_csv(str(csv_path))
            else:
                df = _read_cluster_csv(csv_path, table)
            cluster_col = table["cluster_col"]
            if cluster_col not in df.columns:
                logger.warning(
                    f"cluster_col '{cluster_col}' not found in {csv_path}; skipping"
                )
                continue

            # Clean cluster values (float 5.0 → str "5") and drop NaN
            df[cluster_col] = df[cluster_col].map(_cluster_val_to_str)
            df = df[df[cluster_col] != ""]

            organism = table.get("organism", "")

            # Derive cluster keys from CSV unique values
            unique_clusters = sorted(df[cluster_col].unique())

            # Load extraction JSON
            extraction = _load_extraction_json(self._paperconfig_dir, entry_key)

            # treatment_type may be a list or a string
            treatment_type = table.get("treatment_type", [])
            if isinstance(treatment_type, str):
                treatment_type = [treatment_type]

            # background_factors may be a list or a string
            background_factors = table.get("background_factors", [])
            if isinstance(background_factors, str):
                background_factors = [background_factors]

            # --- ClusteringAnalysis node ---
            analysis_id = _make_analysis_id(self.doi, self.paper_name, entry_key)
            analysis_props = {
                "name": _clean_str(table.get("name", f"Clustering {entry_key}")),
                "organism_name": _clean_str(organism),
                "cluster_method": _clean_str(table.get("cluster_method", "")),
                "cluster_count": len(unique_clusters),
                "total_gene_count": len(df),
                "omics_type": _clean_str(table.get("omics_type", "")),
                "treatment_type": treatment_type,
                "background_factors": background_factors,
                "treatment": _clean_str(table.get("treatment", "")),
                "light_condition": _clean_str(table.get("light_condition", "")),
                "cluster_type": _clean_str(table.get("cluster_type", "")),
                "experimental_context": _clean_str(table.get("experimental_context", "")),
            }
            nodes.append((analysis_id, "clustering_analysis", analysis_props))

            # --- GeneCluster nodes ---
            for cluster_key in unique_clusters:
                cluster_id = _make_cluster_id(
                    self.doi, self.paper_name, entry_key, cluster_key
                )

                mask = df[cluster_col] == cluster_key
                member_count = int(mask.sum())

                # Get extraction data for this cluster (empty if no JSON or verdict != pass)
                ext_data = _get_extraction_cluster_data(extraction, cluster_key)

                props = {
                    "id": _clean_str(ext_data.get("id", "")),
                    "name": _clean_str(ext_data.get("name", f"Cluster {cluster_key}")),
                    "organism_name": _clean_str(organism),
                    "member_count": member_count,
                    "functional_description": _clean_str(
                        ext_data.get("functional_description", "")
                    ),
                    "behavioral_description": _clean_str(
                        ext_data.get("behavioral_description", "")
                    ),
                    "peak_time_hours": ext_data.get("peak_time_hours"),
                    "period_hours": ext_data.get("period_hours"),
                }
                nodes.append((cluster_id, "gene_cluster", props))

        return nodes

    def get_edges(self) -> list[tuple]:
        """Emit edges:
        - gene_in_gene_cluster (GeneCluster → Gene)
        - publication_has_clustering_analysis (Publication → ClusteringAnalysis)
        - clustering_analysis_has_gene_cluster (ClusteringAnalysis → GeneCluster)
        - clusteringanalysis_belongs_to_organism (ClusteringAnalysis → OrganismTaxon)
        - experiment_has_clustering_analysis (Experiment → ClusteringAnalysis)
        """
        edges = []

        for entry_key, table in self._cluster_tables:
            csv_path, use_resolved = _resolve_csv_path(table["filename"])
            if not csv_path.exists():
                continue

            # When using pre-resolved CSV, read with plain pd.read_csv
            # (resolved CSVs are clean comma-separated with no extra header rows)
            if use_resolved:
                df = pd.read_csv(str(csv_path))
            else:
                df = _read_cluster_csv(csv_path, table)
            cluster_col = table["cluster_col"]
            if cluster_col not in df.columns:
                continue

            # Clean cluster values (float 5.0 → str "5") and drop NaN
            df[cluster_col] = df[cluster_col].map(_cluster_val_to_str)
            df = df[df[cluster_col] != ""]

            # Prefer pre-resolved resolved_locus_tag column; fall back to gene_id_col
            _resolved_col = "resolved_locus_tag"
            if use_resolved and _resolved_col in df.columns:
                gene_col = _resolved_col
                logger.info(f"Using pre-resolved CSV: {csv_path.name} ({len(df)} rows)")
            else:
                gene_col = table.get("gene_id_col", "")

            score_col = table.get("score_col")
            p_value_col = table.get("p_value_col")
            organism = table.get("organism", "")

            # Derive cluster keys from CSV
            unique_clusters = sorted(df[cluster_col].unique())

            analysis_id = _make_analysis_id(self.doi, self.paper_name, entry_key)

            # --- Publication → ClusteringAnalysis ---
            if self.doi:
                pub_id = f"doi:{self.doi}"
                pub_edge_id = f"pub_analysis__{analysis_id}"
                edges.append(
                    (pub_edge_id, pub_id, analysis_id, "publication_has_clustering_analysis", {})
                )

            # --- ClusteringAnalysis → OrganismTaxon ---
            if organism:
                org_id = self._organism_lookup.get(organism, "")
                if org_id:
                    org_edge_id = f"analysis_org__{analysis_id}"
                    edges.append(
                        (org_edge_id, analysis_id, org_id, "clusteringanalysis_belongs_to_organism", {})
                    )

            # --- Experiment → ClusteringAnalysis ---
            experiments_list = table.get("experiments", [])
            pub = self.config.get("publication", {})
            for exp_key in experiments_list:
                # Experiment IDs use raw DOI (no doi: prefix)
                if self.doi:
                    experiment_id = f"{self.doi}_{exp_key}"
                else:
                    experiment_id = f"{self.paper_name}_{exp_key}"
                exp_edge_id = f"exp_analysis__{analysis_id}__{exp_key}"
                edges.append(
                    (exp_edge_id, experiment_id, analysis_id, "experiment_has_clustering_analysis", {})
                )

            # --- Per-cluster edges ---
            rows_processed = 0
            for cluster_key in unique_clusters:
                cluster_id = _make_cluster_id(
                    self.doi, self.paper_name, entry_key, cluster_key
                )

                # ClusteringAnalysis → GeneCluster
                analysis_cluster_edge_id = f"analysis_cluster__{cluster_id}"
                edges.append(
                    (analysis_cluster_edge_id, analysis_id, cluster_id, "clustering_analysis_has_gene_cluster", {})
                )

                # Gene_in_gene_cluster edges
                mask = df[cluster_col] == cluster_key
                for _, row in df[mask].iterrows():
                    if self.test_mode and rows_processed >= 100:
                        break

                    gene_locus = row.get(gene_col, "")
                    if pd.isna(gene_locus) or not str(gene_locus).strip():
                        continue
                    gene_locus = str(gene_locus).strip()
                    gene_id = f"ncbigene:{gene_locus}"

                    edge_props = {}
                    if score_col and score_col in df.columns:
                        val = row.get(score_col)
                        if val is not None and pd.notna(val):
                            edge_props["membership_score"] = float(val)
                    if p_value_col and p_value_col in df.columns:
                        val = row.get(p_value_col)
                        if val is not None and pd.notna(val):
                            edge_props["p_value"] = float(val)

                    edge_id = f"{cluster_id}__{gene_locus}"
                    edges.append(
                        (edge_id, cluster_id, gene_id, "gene_in_gene_cluster", edge_props)
                    )
                    rows_processed += 1

        return edges

    def download_data(self, **kwargs):
        """No download needed."""
        pass


class MultiClusterAdapter:
    """Wrapper that reads paperconfig_files.txt and delegates to ClusterAdapter instances."""

    def __init__(
        self,
        config_list_file: str | list[str],
        genome_config_file: str = None,
        test_mode: bool = False,
        **kwargs,
    ):
        self._organism_lookup: dict[str, str] = {}
        if genome_config_file:
            self._organism_lookup = self._build_organism_lookup(genome_config_file)

        self.adapters: list[ClusterAdapter] = []
        list_files = config_list_file if isinstance(config_list_file, list) else [config_list_file]
        paperconfigs = load_all_paperconfigs([Path(lf) for lf in list_files])
        for pc_path, config in paperconfigs:
            supp = config.get("publication", {}).get("supplementary_materials", {})
            has_clusters = any(
                isinstance(v, dict) and v.get("type") == "gene_clusters"
                for v in supp.values()
            )
            if not has_clusters:
                continue
            adapter = ClusterAdapter(
                config_file=str(pc_path), test_mode=test_mode
            )
            adapter._organism_lookup = self._organism_lookup
            self.adapters.append(adapter)

    def _build_organism_lookup(self, genome_config_file: str) -> dict[str, str]:
        """Build {preferred_name: insdc.gcf:<accession>} lookup from genome CSV."""
        lookup: dict[str, str] = {}
        try:
            df = pd.read_csv(genome_config_file)
            for _, row in df.iterrows():
                name = row.get("preferred_name", "")
                accession = row.get("ncbi_accession", "")
                if name and accession:
                    lookup[str(name)] = f"insdc.gcf:{accession}"
        except Exception as e:
            logger.warning(f"Could not load genome config '{genome_config_file}': {e}")
        return lookup

    def download_data(self, **kwargs):
        """No download needed."""
        pass

    def get_nodes(self) -> list[tuple]:
        nodes = []
        for adapter in self.adapters:
            nodes.extend(adapter.get_nodes())
        return nodes

    def get_edges(self) -> list[tuple]:
        edges = []
        for adapter in self.adapters:
            edges.extend(adapter.get_edges())
        return edges
