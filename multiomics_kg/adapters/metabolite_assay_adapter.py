"""
MetaboliteAssayAdapter — reads metabolite_assays_table entries from paperconfig.yaml
files and emits:
  - MetaboliteAssay nodes (one per Experiment x metric_type)
  - publication_has_metabolite_assay edges (Publication -> MetaboliteAssay)
  - experiment_has_metabolite_assay edges (Experiment -> MetaboliteAssay)
  - metabolite_assay_belongs_to_organism edges (MetaboliteAssay -> OrganismTaxon)
  - assay_quantifies_metabolite edges (MetaboliteAssay -> Metabolite, value_kind=numeric)
  - assay_flags_metabolite edges (MetaboliteAssay -> Metabolite, value_kind=boolean)

Sibling of observations_adapter.py — follows the same structural pattern.

Phase 2 of the metabolite scaffold. Reads <stem>_resolved.csv (written by step 7,
resolve_paper_metabolites.py) to get pre-resolved metabolite_id per row. Falls back
to source CSV if no resolved file is present (in which case all rows resolve to
empty metabolite_id and emit no measurement edges).
"""
from __future__ import annotations

import logging
import re
import statistics
from pathlib import Path
from typing import Iterator

import pandas as pd
import yaml

from multiomics_kg.utils.paperconfig_utils import (
    load_paperconfig,
    load_all_paperconfigs,
    get_paper_name,
    get_experiments,
    iter_metabolite_assays_tables,
)

logger = logging.getLogger(__name__)

DEFAULT_NULL_VALUES = {"nd", "ND", "n.d.", "NA", "N/A"}
DEFAULT_MISSING_VALUES = {""}

_EMBEDDED_PATTERN = re.compile(
    r"^\s*([0-9.+\-eE]+)\s*\(\s*([0-9.+\-eEnNaA/]+)\s*\)\s*,\s*n\s*=\s*(\d+)\s*$"
)


def _clean_str(value) -> str:
    """Sanitize string for BioCypher CSV output (CLAUDE.md convention)."""
    if value is None:
        return ""
    if not isinstance(value, str):
        return str(value)
    return value.replace("'", "^").replace("|", ",")


def _make_metabolite_assay_id(
    doi: str, paper_name: str, entry_key: str, metric_type: str
) -> str:
    """MetaboliteAssay node ID: metabolite_assay:{doi_short}:{entry_key}:{metric_type}."""
    if doi:
        doi_short = doi.rstrip("/").rsplit("/", 1)[-1]
    else:
        doi_short = re.sub(r"[^a-z0-9]+", "_", paper_name.lower()).strip("_")
    return f"metabolite_assay:{doi_short}:{entry_key}:{metric_type}"


def _resolve_csv_path(csv_path: str) -> tuple[Path, bool]:
    """Probe for <stem>_resolved.csv (written by step 7).

    Returns (path_to_use, is_resolved). If no _resolved.csv exists, returns
    the source path with is_resolved=False (caller will see no metabolite_id
    column and emit no measurement edges).
    """
    p = Path(csv_path)
    resolved = p.with_name(p.stem + "_resolved.csv")
    if resolved.exists():
        return resolved, True
    return p, False


def _aggregate_replicates(
    raw_values,
    null_values: set[str],
    missing_values: set[str],
) -> tuple[float, float, int, int, list[float], str]:
    """Aggregate replicate cells into (mean, sd, n_replicates, n_non_zero, replicate_values, detection_status).

    null_values  → treated as 0.0 (not-detected, counted in n_replicates).
    missing_values → row excluded from aggregation entirely.
    Anything else → coerced to float; non-numeric cells are skipped.
    """
    parsed: list[float] = []
    for v in raw_values:
        if v is None or (isinstance(v, float) and pd.isna(v)):
            s = ""
        else:
            s = str(v).strip()
        if s in missing_values:
            continue
        if s in null_values:
            parsed.append(0.0)
            continue
        try:
            parsed.append(float(s))
        except (TypeError, ValueError):
            # Bad cell — treat as missing (skip)
            continue

    n_replicates = len(parsed)
    if n_replicates == 0:
        return 0.0, 0.0, 0, 0, [], "not_detected"
    n_non_zero = sum(1 for x in parsed if x != 0.0)
    mean = statistics.fmean(parsed)
    sd = statistics.stdev(parsed) if n_replicates >= 2 else 0.0

    if n_non_zero == 0:
        det = "not_detected"
    elif n_non_zero == n_replicates:
        det = "detected"
    else:
        det = "sporadic"

    return mean, sd, n_replicates, n_non_zero, parsed, det


def parse_embedded_mean_sd_n(cell) -> tuple[float, float, int] | None:
    """Parse cells like '0.00054 (8.8e-05), n=2' or 'nd'.

    Returns (mean, sd, n) or None for empty/unparseable cells. 'nd'/'ND' → (0,0,0).
    Non-numeric sd (e.g. 'NA') → 0.0.
    """
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return None
    s = str(cell).strip()
    if not s:
        return None
    if s.lower() in {"nd", "n.d."}:
        return 0.0, 0.0, 0
    m = _EMBEDDED_PATTERN.match(s)
    if not m:
        return None
    mean = float(m.group(1))
    try:
        sd = float(m.group(2))
    except ValueError:
        sd = 0.0
    n = int(m.group(3))
    return mean, sd, n


class MetaboliteAssayAdapter:
    """Adapter for one paperconfig's metabolite_assays_table entries."""

    def __init__(self, config_file: str, test_mode: bool = False):
        self.config_file = config_file
        self.test_mode = test_mode
        self.config = load_paperconfig(Path(config_file))
        self.paper_name = get_paper_name(self.config, fallback_path=Path(config_file))
        self.doi = (self.config.get("publication") or {}).get("doi", "") or ""
        self._entries = list(iter_metabolite_assays_tables(self.config))

    def _denormalized_fields(self, experiment: dict) -> dict:
        """Denormalized parent-Experiment block on every MetaboliteAssay node.

        Mirrors ObservationsAdapter._denormalized_fields. compartment defaults
        to whole_cell when the Experiment doesn't declare one.
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
            "compartment": _clean_str(experiment.get("compartment", "whole_cell") or "whole_cell"),
            "omics_type": _clean_str(experiment.get("omics_type", "METABOLOMICS")),
            "treatment_type": _as_list(experiment.get("treatment_type", [])),
            "background_factors": _as_list(experiment.get("background_factors", [])),
            "treatment": _clean_str(experiment.get("treatment_condition", "")),
            "light_condition": _clean_str(experiment.get("light_condition", "")),
            "experimental_context": _clean_str(experiment.get("experimental_context", "")),
        }

    def get_nodes(self) -> Iterator[tuple]:
        """Emit MetaboliteAssay nodes — one per (entry × metric_type)."""
        experiments = get_experiments(self.config)
        for entry_key, entry in self._entries:
            exp_key = entry.get("experiment")
            if not exp_key or exp_key not in experiments:
                logger.warning(
                    f"metabolite_assays_table '{entry_key}' references unknown "
                    f"experiment '{exp_key}' — skipping"
                )
                continue
            exp = experiments[exp_key]
            denorm = self._denormalized_fields(exp)
            organism = entry.get("organism") or exp.get("organism") or ""

            for assay in entry.get("assays", []):
                metric_type = assay.get("metric_type", "")
                value_kind = assay.get("value_kind", "")
                if not metric_type or not value_kind:
                    logger.warning(
                        f"assay in '{entry_key}' missing metric_type or value_kind — skipping"
                    )
                    continue
                node_id = _make_metabolite_assay_id(self.doi, self.paper_name, entry_key, metric_type)
                experiment_id = (
                    f"{self.doi}_{exp_key}" if self.doi else f"{self.paper_name}_{exp_key}"
                )
                props = {
                    "name": _clean_str(assay.get("name", metric_type)),
                    "experiment_id": _clean_str(experiment_id),
                    "metric_type": _clean_str(metric_type),
                    "value_kind": _clean_str(value_kind),
                    "unit": _clean_str(assay.get("unit", "")),
                    "rankable": _clean_str(assay.get("rankable", "false")),
                    "aggregation_method": _clean_str(
                        assay.get("aggregation_method")
                        or entry.get("aggregation_method", "mean_across_replicates")
                    ),
                    "field_description": _clean_str(assay.get("field_description", "")),
                    **denorm,
                }
                yield node_id, "metabolite_assay", props

    def get_edges(self) -> Iterator[tuple]:
        """Emit measurement + binding edges per (assay, metabolite, condition)."""
        experiments = get_experiments(self.config)
        for entry_key, entry in self._entries:
            exp_key = entry.get("experiment")
            if not exp_key or exp_key not in experiments:
                continue
            csv_str = entry.get("filename", "")
            if not csv_str:
                continue
            csv_path, is_resolved = _resolve_csv_path(csv_str)
            if not csv_path.exists():
                logger.warning(f"metabolite_assays_table CSV not found: {csv_path}")
                continue
            try:
                df = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
            except Exception as e:
                logger.warning(f"Cannot read {csv_path}: {e}")
                continue

            if not is_resolved:
                logger.warning(
                    f"metabolite_assays_table '{entry_key}': no _resolved.csv next to "
                    f"{csv_path} — run step 7 (prepare_data.sh --steps 7). "
                    f"Skipping measurement edges."
                )
                # Still emit binding edges so MetaboliteAssay nodes aren't orphans
                for assay in entry.get("assays", []):
                    metric_type = assay.get("metric_type", "")
                    if not metric_type:
                        continue
                    assay_node_id = _make_metabolite_assay_id(
                        self.doi, self.paper_name, entry_key, metric_type
                    )
                    yield from self._emit_binding_edges(assay_node_id, exp_key, entry)
                continue

            null_values = set(entry.get("null_values") or DEFAULT_NULL_VALUES)
            missing_values = set(entry.get("missing_values") or DEFAULT_MISSING_VALUES)

            for assay in entry.get("assays", []):
                metric_type = assay.get("metric_type", "")
                value_kind = assay.get("value_kind", "")
                if not metric_type or not value_kind:
                    continue
                assay_node_id = _make_metabolite_assay_id(
                    self.doi, self.paper_name, entry_key, metric_type
                )

                if value_kind == "numeric":
                    yield from self._emit_numeric_edges(
                        df, entry, assay, assay_node_id, null_values, missing_values,
                    )
                elif value_kind == "boolean":
                    yield from self._emit_boolean_edges(df, entry, assay, assay_node_id)

                yield from self._emit_binding_edges(assay_node_id, exp_key, entry)

    def _emit_numeric_edges(
        self, df, entry, assay, assay_node_id, null_values, missing_values,
    ):
        sample_cols_block = assay.get("sample_columns", [])
        drop_undetected = bool(entry.get("drop_undetected", False))
        for row_idx, row in df.iterrows():
            primary = str(row.get("metabolite_id", "") or "").strip()
            if not primary:
                continue
            for sc_idx, sc in enumerate(sample_cols_block):
                rcols = sc.get("replicate_columns") or []
                if not rcols:
                    continue
                raw_values = [row.get(c, "") for c in rcols]
                mean, sd, n_rep, n_nz, vals, det = _aggregate_replicates(
                    raw_values, null_values=null_values, missing_values=missing_values,
                )
                if n_rep == 0:
                    continue  # all-missing → no edge regardless of policy
                if drop_undetected and det == "not_detected":
                    continue
                cond_label = _clean_str(sc.get("condition_label", ""))
                edge_id = f"{assay_node_id}|{primary}|{sc_idx}|{row_idx}"
                props = {
                    "metric_type": _clean_str(assay.get("metric_type", "")),
                    "condition_label": cond_label,
                    "time_point": _clean_str(sc.get("time_point", "")),
                    "time_point_order": int(sc.get("time_point_order") or 0),
                    "time_point_hours": float(sc.get("time_point_hours", -1.0) or -1.0),
                    "value": float(mean),
                    "value_sd": float(sd),
                    "n_replicates": int(n_rep),
                    "n_non_zero": int(n_nz),
                    "replicate_values": vals,
                    "detection_status": det,
                }
                yield edge_id, assay_node_id, primary, "assay_quantifies_metabolite", props

    def _emit_boolean_edges(self, df, entry, assay, assay_node_id):
        sample_cols_block = assay.get("sample_columns", [])
        for row_idx, row in df.iterrows():
            primary = str(row.get("metabolite_id", "") or "").strip()
            if not primary:
                continue
            for sc_idx, sc in enumerate(sample_cols_block):
                flag_col = sc.get("flag_column")
                if not flag_col:
                    continue
                flag_true_value = str(sc.get("flag_true_value", "yes"))
                cell = str(row.get(flag_col, "") or "").strip()
                is_true = cell == flag_true_value
                cond_label = _clean_str(sc.get("condition_label", ""))
                edge_id = f"{assay_node_id}|{primary}|{sc_idx}|{row_idx}"
                props = {
                    "metric_type": _clean_str(assay.get("metric_type", "")),
                    "condition_label": cond_label,
                    "flag_value": "true" if is_true else "false",
                    "n_replicates": 1,
                    "n_positive": 1 if is_true else 0,
                }
                yield edge_id, assay_node_id, primary, "assay_flags_metabolite", props

    def _emit_binding_edges(self, assay_node_id, exp_key, entry):
        exp = (get_experiments(self.config) or {}).get(exp_key, {}) or {}
        organism = entry.get("organism") or exp.get("organism") or ""
        experiment_id = (
            f"{self.doi}_{exp_key}" if self.doi else f"{self.paper_name}_{exp_key}"
        )
        pub_id = self.doi or f"papername:{self.paper_name}"

        yield (
            f"pub__{assay_node_id}", pub_id, assay_node_id,
            "publication_has_metabolite_assay", {},
        )
        yield (
            f"exp__{assay_node_id}", experiment_id, assay_node_id,
            "experiment_has_metabolite_assay", {},
        )
        yield (
            f"org__{assay_node_id}", assay_node_id, organism,
            "metabolite_assay_belongs_to_organism", {},
        )


class MultiMetaboliteAssayAdapter:
    """Registry-driven wrapper. Reads paperconfig_files.txt, instantiates one
    MetaboliteAssayAdapter per paperconfig that has at least one
    metabolite_assays_table entry."""

    def __init__(self, paperconfig_paths=None, test_mode: bool = False):
        if paperconfig_paths is None:
            paperconfigs = load_all_paperconfigs()
            paperconfig_paths = [p for p, _ in paperconfigs]
        self.adapters: list[MetaboliteAssayAdapter] = []
        for p in paperconfig_paths:
            try:
                with open(p) as f:
                    cfg = yaml.safe_load(f) or {}
            except OSError:
                continue
            sm = (cfg.get("publication") or {}).get("supplementary_materials") or {}
            if any(
                isinstance(v, dict) and v.get("type") == "metabolite_assays_table"
                for v in sm.values()
            ):
                self.adapters.append(
                    MetaboliteAssayAdapter(config_file=str(p), test_mode=test_mode)
                )

    def get_nodes(self):
        for a in self.adapters:
            yield from a.get_nodes()

    def get_edges(self):
        for a in self.adapters:
            yield from a.get_edges()
