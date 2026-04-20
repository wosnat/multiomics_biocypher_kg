"""
ObservationsAdapter — reads derived_metrics_table entries from paperconfig.yaml
files and emits:
  - DerivedMetric nodes (one per Experiment x metric_type)
  - publication_has_derived_metric edges (Publication -> DerivedMetric)
  - experiment_has_derived_metric edges (Experiment -> DerivedMetric)
  - derived_metric_belongs_to_organism edges (DerivedMetric -> OrganismTaxon)
  - derived_metric_quantifies_gene edges (DerivedMetric -> Gene, value_kind=numeric)
  - derived_metric_flags_gene edges (DerivedMetric -> Gene, value_kind=boolean)
  - derived_metric_classifies_gene edges (DerivedMetric -> Gene, value_kind=categorical)

Sibling of cluster_adapter.py — follows the same structural pattern.
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
    get_experiments,
    iter_derived_metrics_tables,
)
from multiomics_kg.vocab.non_de_evidence import DEFAULT_SKIP_TOKENS, VALID_BLANK_POLICIES

logger = logging.getLogger(__name__)

_DEFAULT_PDF_CACHE = Path(__file__).parent.parent.parent / "cache" / "pdf_extraction_cache.json"


def _clean_str(value) -> str:
    """Sanitize string for BioCypher CSV output (CLAUDE.md convention)."""
    if value is None:
        return ""
    if not isinstance(value, str):
        return str(value)
    return value.replace("'", "^").replace("|", ",")


def _make_derived_metric_id(
    doi: str, paper_name: str, entry_key: str, metric_type: str
) -> str:
    """DerivedMetric node ID: derived_metric:{doi_short}:{entry_key}:{metric_type}."""
    if doi:
        doi_short = doi.rstrip("/").rsplit("/", 1)[-1]
    else:
        doi_short = re.sub(r"[^a-z0-9]+", "_", paper_name.lower()).strip("_")
    return f"derived_metric:{doi_short}:{entry_key}:{metric_type}"


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
    if cache_path.exists():
        try:
            with open(cache_path) as f:
                return json.load(f)
        except Exception:
            pass
    return {}


def _parse_boolean_cell(
    value,
    true_tokens: list[str],
    false_tokens: list[str],
    skip_tokens: list[str],
    blank_policy: str,
) -> str | None:
    """Map a single CSV cell value to "true" / "false" / None (=skip).

    Hard-errors on unexpected tokens — per parent spec, the adapter must not
    silently coerce ambiguous values.
    """
    # NaN / None / empty string -> apply blank_policy
    if value is None or pd.isna(value):
        return _apply_blank_policy(blank_policy)
    s = str(value).strip()
    if s == "":
        return _apply_blank_policy(blank_policy)
    if s in true_tokens:
        return "true"
    if s in false_tokens:
        return "false"
    if s in skip_tokens:
        return None
    raise ValueError(
        f"Unexpected boolean token {s!r}: not in true_tokens={true_tokens}, "
        f"false_tokens={false_tokens}, or skip_tokens={skip_tokens}. "
        f"Paperconfig author must classify explicitly."
    )


def _apply_blank_policy(blank_policy: str) -> str | None:
    if blank_policy == "skip":
        return None
    if blank_policy == "true":
        return "true"
    if blank_policy == "false":
        return "false"
    raise ValueError(
        f"Invalid blank_policy {blank_policy!r}; must be one of {VALID_BLANK_POLICIES}"
    )


# --- Adapter classes (get_nodes / get_edges land in Tasks 5-11) ---


class ObservationsAdapter:
    """Adapter for one paperconfig's derived_metrics_table entries."""

    def __init__(self, config_file: str, test_mode: bool = False):
        self.config_file = config_file
        self.test_mode = test_mode
        self.config = load_paperconfig(Path(config_file))
        self.paper_name = get_paper_name(self.config, fallback_path=Path(config_file))
        self.doi = self._extract_doi()
        self._dm_entries = list(iter_derived_metrics_tables(self.config))
        self._organism_lookup: dict[str, str] = {}

    def _extract_doi(self) -> str:
        pub = self.config.get("publication", {})
        doi = pub.get("doi", "")
        if doi:
            return doi
        pdf_path = pub.get("papermainpdf", "")
        if pdf_path:
            cache = _load_pdf_cache()
            cached = cache.get(pdf_path, {})
            doi = cached.get("publication", {}).get("doi", "")
        return doi or ""

    def get_nodes(self) -> list[tuple]:
        return []  # filled in Tasks 5-7

    def get_edges(self) -> list[tuple]:
        return []  # filled in Tasks 8-11

    def download_data(self, **kwargs):
        pass


class MultiObservationsAdapter:
    """Wrapper that reads paperconfig_files.txt and delegates."""

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
        self.adapters: list[ObservationsAdapter] = []
        list_files = config_list_file if isinstance(config_list_file, list) else [config_list_file]
        paperconfigs = load_all_paperconfigs([Path(lf) for lf in list_files])
        for pc_path, config in paperconfigs:
            supp = config.get("publication", {}).get("supplementary_materials", {})
            has_dm = any(
                isinstance(v, dict) and v.get("type") == "derived_metrics_table"
                for v in supp.values()
            )
            if not has_dm:
                continue
            adapter = ObservationsAdapter(config_file=str(pc_path), test_mode=test_mode)
            adapter._organism_lookup = self._organism_lookup
            self.adapters.append(adapter)

    def _build_organism_lookup(self, genome_config_file: str) -> dict[str, str]:
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
