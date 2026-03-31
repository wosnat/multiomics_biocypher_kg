"""
Gene Cluster Adapter

Reads gene_clusters entries from paperconfig.yaml files and emits
GeneCluster nodes with Gene_in_gene_cluster membership edges,
Publication_has_gene_cluster edges, and Genecluster_belongs_to_organism edges.
"""
import json
import logging
from pathlib import Path

import pandas as pd

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


def _make_cluster_id(doi: str, paper_name: str, organism: str, cluster_key: str) -> str:
    """Build cluster node ID: cluster:{doi_short}:{organism_suffix}:{cluster_key}.

    Uses DOI short form when available, falls back to paper_name slug.
    Adds organism suffix to disambiguate when a paper has clusters for
    multiple organisms (e.g., MED4 + MIT9313).
    """
    if doi:
        doi_short = doi.rsplit("/", 1)[-1] if "/" in doi else doi
    else:
        doi_short = paper_name.lower().replace(" ", "_")
    # Short organism suffix: "Prochlorococcus MED4" → "med4"
    org_suffix = organism.split()[-1].lower() if organism else "unknown"
    return f"cluster:{doi_short}:{org_suffix}:{cluster_key}"


def _resolve_csv_path(csv_path: str) -> Path:
    """Return the resolved CSV path if it exists, else the original path."""
    p = Path(csv_path)
    resolved = p.parent / f"{p.stem}_resolved{p.suffix}"
    if resolved.exists():
        return resolved
    return p


def _load_pdf_cache(cache_path: Path = _DEFAULT_PDF_CACHE) -> dict:
    """Load the PDF extraction cache (DOI, title, etc.)."""
    if cache_path.exists():
        try:
            with open(cache_path) as f:
                return json.load(f)
        except Exception:
            pass
    return {}


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
        """Emit GeneCluster nodes."""
        nodes = []
        for table_key, table in self._cluster_tables:
            csv_path = _resolve_csv_path(table["filename"])
            if not csv_path.exists():
                logger.warning(f"Cluster CSV not found: {csv_path}")
                continue

            df = pd.read_csv(csv_path)
            cluster_col = table["cluster_col"]
            if cluster_col not in df.columns:
                logger.warning(
                    f"cluster_col '{cluster_col}' not found in {csv_path}; skipping"
                )
                continue

            organism = table.get("organism", "")
            clusters_meta = table.get("clusters", {})
            for cluster_key, cluster_info in clusters_meta.items():
                # Use explicit id if provided, else fall back to the YAML key
                id_slug = cluster_info.get("id", cluster_key)
                cluster_id = _make_cluster_id(
                    self.doi, self.paper_name, organism, id_slug
                )

                mask = df[cluster_col].astype(str) == str(cluster_key)
                member_count = int(mask.sum())

                # treatment_type may be a list or a string
                treatment_type = table.get("treatment_type", [])
                if isinstance(treatment_type, str):
                    treatment_type = [treatment_type]

                props = {
                    "name": _clean_str(cluster_info.get("name", f"Cluster {cluster_key}")),
                    "source_paper": _clean_str(self.paper_name),
                    "organism_name": _clean_str(table.get("organism", "")),
                    "cluster_method": _clean_str(table.get("cluster_method", "")),
                    "cluster_type": _clean_str(cluster_info.get("cluster_type", "")),
                    "treatment_type": treatment_type,
                    "treatment": _clean_str(table.get("treatment", "")),
                    "omics_type": _clean_str(table.get("omics_type", "")),
                    "light_condition": _clean_str(table.get("light_condition", "")),
                    "member_count": member_count,
                    "functional_description": _clean_str(
                        cluster_info.get("functional_description", "")
                    ),
                    "behavioral_description": _clean_str(
                        cluster_info.get("behavioral_description", "")
                    ),
                    "peak_time_hours": cluster_info.get("peak_time_hours"),
                    "period_hours": cluster_info.get("period_hours"),
                    "experimental_context": _clean_str(
                        table.get("experimental_context", "")
                    ),
                }
                nodes.append((cluster_id, "gene_cluster", props))

        return nodes

    def get_edges(self) -> list[tuple]:
        """Emit Gene_in_gene_cluster, Publication_has_gene_cluster, and
        Genecluster_belongs_to_organism edges."""
        edges = []

        for table_key, table in self._cluster_tables:
            csv_path = _resolve_csv_path(table["filename"])
            if not csv_path.exists():
                continue

            df = pd.read_csv(csv_path)
            cluster_col = table["cluster_col"]
            if cluster_col not in df.columns:
                continue

            # Prefer pre-resolved locus_tag column; fall back to gene_id_col
            if "locus_tag" in df.columns:
                gene_col = "locus_tag"
            else:
                gene_col = table.get("gene_id_col", "")

            score_col = table.get("score_col")
            p_value_col = table.get("p_value_col")
            organism = table.get("organism", "")
            clusters_meta = table.get("clusters", {})

            rows_processed = 0
            for cluster_key, cluster_info in clusters_meta.items():
                id_slug = cluster_info.get("id", cluster_key)
                cluster_id = _make_cluster_id(
                    self.doi, self.paper_name, organism, id_slug
                )

                # Gene_in_gene_cluster edges
                mask = df[cluster_col].astype(str) == str(cluster_key)
                for _, row in df[mask].iterrows():
                    if self.test_mode and rows_processed >= 100:
                        break

                    gene_locus = row.get("locus_tag") if "locus_tag" in df.columns else row.get(gene_col, "")
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

                # Publication_has_gene_cluster edge
                if self.doi:
                    pub_id = f"doi:{self.doi}"
                    pub_edge_id = f"pub_cluster__{cluster_id}"
                    edges.append(
                        (
                            pub_edge_id,
                            pub_id,
                            cluster_id,
                            "publication_has_gene_cluster",
                            {},
                        )
                    )

                # Genecluster_belongs_to_organism edge
                if organism:
                    org_id = self._organism_lookup.get(organism, "")
                    if org_id:
                        org_edge_id = f"cluster_org__{cluster_id}"
                        edges.append(
                            (
                                org_edge_id,
                                cluster_id,
                                org_id,
                                "genecluster_belongs_to_organism",
                                {},
                            )
                        )

        return edges

    def download_data(self, **kwargs):
        """No download needed."""
        pass


class MultiClusterAdapter:
    """Wrapper that reads paperconfig_files.txt and delegates to ClusterAdapter instances."""

    def __init__(
        self,
        config_list_file: str,
        genome_config_file: str = None,
        test_mode: bool = False,
        **kwargs,
    ):
        self._organism_lookup: dict[str, str] = {}
        if genome_config_file:
            self._organism_lookup = self._build_organism_lookup(genome_config_file)

        self.adapters: list[ClusterAdapter] = []
        paperconfigs = load_all_paperconfigs(Path(config_list_file))
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
