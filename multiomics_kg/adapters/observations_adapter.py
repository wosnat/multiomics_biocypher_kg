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

    def _denormalized_fields(self, experiment: dict) -> dict:
        """Derive the 9 parent-Experiment fields every DerivedMetric copies.

        Enforces spec invariant #9: DerivedMetric denormalized fields equal
        parent Experiment's values. Never re-read from the paperconfig entry.
        """
        def _as_list(v):
            if v is None:
                return []
            if isinstance(v, str):
                return [_clean_str(v)] if v else []
            return [_clean_str(x) for x in v]

        return {
            "organism_name": _clean_str(experiment.get("organism", "")),
            "publication_doi": _clean_str(self.doi),
            "compartment": _clean_str(experiment.get("compartment", "whole_cell")),
            "omics_type": _clean_str(experiment.get("omics_type", "")),
            "treatment_type": _as_list(experiment.get("treatment_type", [])),
            "background_factors": _as_list(experiment.get("background_factors", [])),
            "treatment": _clean_str(experiment.get("treatment_condition", "")),
            "light_condition": _clean_str(experiment.get("light_condition", "")),
            "experimental_context": _clean_str(experiment.get("experimental_context", "")),
        }

    def get_nodes(self) -> list[tuple]:
        """Emit DerivedMetric nodes — one per (entry × metric_type).

        Matches cluster_adapter's behaviour: if the source CSV is missing,
        skip the whole entry so we never emit dangling DM nodes with no
        measurement edges.
        """
        nodes = []
        experiments = get_experiments(self.config)
        for entry_key, entry in self._dm_entries:
            exp_key = entry.get("experiment")
            if not exp_key or exp_key not in experiments:
                logger.warning(
                    f"derived_metrics_table '{entry_key}' references unknown "
                    f"experiment '{exp_key}' — skipping"
                )
                continue
            csv_path, _ = _resolve_csv_path(entry["filename"])
            if not csv_path.exists():
                logger.warning(
                    f"derived_metrics_table CSV not found: {csv_path} — skipping"
                )
                continue
            exp = experiments[exp_key]
            denorm = self._denormalized_fields(exp)

            for metric in entry.get("metrics", []):
                metric_type = metric.get("metric_type", "")
                value_kind = metric.get("value_kind", "")
                if not metric_type or not value_kind:
                    logger.warning(
                        f"Metric in '{entry_key}' missing metric_type or value_kind — skipping"
                    )
                    continue

                dm_id = _make_derived_metric_id(
                    self.doi, self.paper_name, entry_key, metric_type
                )
                # experiment_id mirrors the omics_adapter's format: raw-doi + "_" + exp_key
                experiment_id = f"{self.doi}_{exp_key}" if self.doi else f"{self.paper_name}_{exp_key}"

                props = {
                    "name": _clean_str(metric.get("name", metric_type)),
                    "experiment_id": _clean_str(experiment_id),
                    "metric_type": _clean_str(metric_type),
                    "value_kind": _clean_str(value_kind),
                    "field_description": _clean_str(metric.get("field_description", "")),
                    **denorm,
                }

                if value_kind == "boolean":
                    props["rankable"] = "false"
                    props["has_p_value"] = "false"
                    props["unit"] = ""
                    props["allowed_categories"] = []
                elif value_kind == "categorical":
                    # Task 6 completes the categorical branch
                    props["rankable"] = "false"
                    props["has_p_value"] = "false"
                    props["unit"] = ""
                    ac = metric.get("allowed_categories", [])
                    if isinstance(ac, str):
                        ac = [ac] if ac else []
                    elif ac is None:
                        ac = []
                    props["allowed_categories"] = [_clean_str(c) for c in ac]
                elif value_kind == "numeric":
                    # Task 7 completes the numeric branch
                    props["rankable"] = _clean_str(metric.get("rankable", "false"))
                    props["has_p_value"] = _clean_str(metric.get("has_p_value", "false"))
                    props["unit"] = _clean_str(metric.get("unit", ""))
                    props["allowed_categories"] = []
                    pvt = metric.get("p_value_threshold")
                    if pvt is not None and pvt != "":
                        try:
                            props["p_value_threshold"] = float(pvt)
                        except (TypeError, ValueError):
                            logger.warning(
                                f"Invalid p_value_threshold {pvt!r} in '{entry_key}/{metric_type}' "
                                f"— expected float, skipping property"
                            )
                else:
                    logger.warning(
                        f"Unknown value_kind '{value_kind}' in '{entry_key}/{metric_type}' — skipping"
                    )
                    continue

                nodes.append((dm_id, "derived_metric", props))
        return nodes

    def get_edges(self) -> list[tuple]:
        """Emit binding + measurement edges for each DerivedMetric."""
        edges = []
        experiments = get_experiments(self.config)
        pub_id = f"doi:{self.doi}" if self.doi else None

        for entry_key, entry in self._dm_entries:
            exp_key = entry.get("experiment")
            if not exp_key or exp_key not in experiments:
                continue
            exp = experiments[exp_key]
            organism = exp.get("organism", "")

            csv_path, use_resolved = _resolve_csv_path(entry["filename"])
            if not csv_path.exists():
                logger.warning(f"derived_metrics_table CSV not found: {csv_path}")
                continue
            try:
                df = pd.read_csv(csv_path)
            except Exception as e:
                logger.warning(f"Failed to read derived_metrics_table CSV {csv_path}: {e}")
                continue

            # Column used to reach the Gene node
            if use_resolved and "resolved_locus_tag" in df.columns:
                gene_col = "resolved_locus_tag"
                logger.info(
                    f"Using pre-resolved CSV: {csv_path.name} ({len(df)} rows)"
                )
            else:
                gene_col = entry.get("name_col", "")

            experiment_id = (
                f"{self.doi}_{exp_key}" if self.doi else f"{self.paper_name}_{exp_key}"
            )

            for metric in entry.get("metrics", []):
                metric_type = metric.get("metric_type", "")
                value_kind = metric.get("value_kind", "")
                if not metric_type or not value_kind:
                    continue

                dm_id = _make_derived_metric_id(
                    self.doi, self.paper_name, entry_key, metric_type
                )

                # --- Binding edges ---
                if pub_id:
                    edges.append((
                        f"pub_dm__{dm_id}",
                        pub_id, dm_id, "publication_has_derived_metric", {},
                    ))
                edges.append((
                    f"exp_dm__{dm_id}__{exp_key}",
                    experiment_id, dm_id, "experiment_has_derived_metric", {},
                ))
                if organism and organism in self._organism_lookup:
                    edges.append((
                        f"dm_org__{dm_id}",
                        dm_id, self._organism_lookup[organism],
                        "derived_metric_belongs_to_organism", {},
                    ))

                # --- Measurement edges (Tasks 9–11) ---
                value_col = metric.get("value_col", "")
                if not value_col:
                    logger.warning(
                        f"'{entry_key}/{metric_type}' missing value_col — skipping measurement edges"
                    )
                    continue
                if value_col not in df.columns:
                    logger.warning(
                        f"value_col '{value_col}' not in {csv_path.name} for "
                        f"'{entry_key}/{metric_type}' — skipping measurement edges"
                    )
                    continue

                if value_kind == "boolean":
                    edges.extend(self._emit_boolean_edges(
                        df, gene_col, value_col, metric, dm_id, metric_type,
                    ))
                elif value_kind == "categorical":
                    edges.extend(self._emit_categorical_edges(
                        df, gene_col, value_col, metric, dm_id, metric_type,
                    ))
                elif value_kind == "numeric":
                    edges.extend(self._emit_numeric_edges(
                        df, gene_col, value_col, metric, dm_id, metric_type,
                    ))

        return edges

    def _emit_boolean_edges(self, df, gene_col, value_col, metric, dm_id, metric_type):
        """Emit derived_metric_flags_gene edges.

        Precondition: value_col is present in df.columns (guarded in get_edges).
        Task 9 implements this.
        """
        return []

    def _emit_categorical_edges(self, df, gene_col, value_col, metric, dm_id, metric_type):
        """Emit derived_metric_classifies_gene edges.

        Precondition: value_col is present in df.columns (guarded in get_edges).
        Task 10 implements this.
        """
        return []

    def _emit_numeric_edges(self, df, gene_col, value_col, metric, dm_id, metric_type):
        """Emit derived_metric_quantifies_gene edges.

        Precondition: value_col is present in df.columns (guarded in get_edges).
        Task 11 implements this.
        """
        return []

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
