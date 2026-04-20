"""Controlled vocabularies for non-DE-evidence schema (DerivedMetric today;
AbundanceAnalysis in a future slice).

Design note: the registry is deliberately narrow. Only value_kind (which
drives adapter edge-type dispatch) is enforced centrally. Per-metric
metadata — rankable, has_p_value, p_value_threshold, unit,
allowed_categories, field_description — is declared inline on paperconfig
entries, so one paper's idiosyncratic metric doesn't force global vocab churn.
"""
from __future__ import annotations

COMPARTMENTS: frozenset[str] = frozenset()
EXTENDED_OMICS_TYPES: frozenset[str] = frozenset()
VALUE_KINDS: frozenset[str] = frozenset()
KNOWN_METRIC_TYPES: dict[str, str] = {}

DEFAULT_SKIP_TOKENS: tuple[str, ...] = ()
VALID_BLANK_POLICIES: tuple[str, ...] = ()

BUCKET_THRESHOLD_TOP_DECILE: int = 0
BUCKET_THRESHOLD_TOP_QUARTILE: int = 0
BUCKET_THRESHOLD_MID: int = 0
