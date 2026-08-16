"""Controlled-vocabulary contract — the single source of truth.

Adapters import their literals from here so an undeclared value cannot be
emitted; the ControlledVocabulary adapter publishes the same data into the
graph; kg-validity asserts the live graph agrees.

See docs/superpowers/specs/2026-08-16-vocabulary-contract-design.md §5.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

DEFAULT_PATH = Path("config/controlled_vocabularies.yaml")

# R5: a two-state fact is a meaningful categorical string. Native bool is not
# admissible — BioCypher does not round-trip it, and a bare `true` is
# unreadable in a result row.
VALUE_TYPES = {"string", "string_array", "float", "int", "bool_string"}

_REQUIRED = ("applies_to", "applies_to_kind", "property", "value_type",
             "closed", "description")


@dataclass(frozen=True)
class VocabEntry:
    applies_to: str
    applies_to_kind: str          # 'node' | 'edge'
    property: str
    value_type: str
    closed: bool
    description: str
    values: list[str] = field(default_factory=list)
    sparse: bool = False
    expected_empty: bool = False
    min_value: float | None = None
    max_value: float | None = None
    signal_count: int | None = None
    signals: list[str] = field(default_factory=list)

    @property
    def id(self) -> str:
        return f"{self.applies_to}.{self.property}"


def _validate(key: str, raw: dict[str, Any]) -> None:
    missing = [k for k in _REQUIRED if k not in raw]
    if missing:
        raise ValueError(
            f"controlled_vocabularies.yaml entry '{key}' is missing required "
            f"key(s): {', '.join(missing)}"
        )
    vt = raw["value_type"]
    if vt not in VALUE_TYPES:
        raise ValueError(
            f"controlled_vocabularies.yaml entry '{key}' declares "
            f"value_type='{vt}'. Allowed: {sorted(VALUE_TYPES)}. Note 'bool' "
            f"is forbidden by house rule R5 — use a meaningful two-state "
            f"string instead."
        )
    if raw["applies_to_kind"] not in ("node", "edge"):
        raise ValueError(
            f"controlled_vocabularies.yaml entry '{key}' has "
            f"applies_to_kind='{raw['applies_to_kind']}'; expected node|edge."
        )
    expected_id = f"{raw['applies_to']}.{raw['property']}"
    if key != expected_id:
        raise ValueError(
            f"controlled_vocabularies.yaml key '{key}' does not match its own "
            f"applies_to/property ('{expected_id}'). The key IS the node id."
        )


def load_vocabularies(path: str | Path = DEFAULT_PATH) -> dict[str, VocabEntry]:
    """Parse the YAML into VocabEntry objects keyed by '<applies_to>.<property>'."""
    with open(path) as f:
        raw_all = yaml.safe_load(f) or {}
    out: dict[str, VocabEntry] = {}
    for key, raw in raw_all.items():
        _validate(key, raw)
        out[key] = VocabEntry(
            applies_to=raw["applies_to"],
            applies_to_kind=raw["applies_to_kind"],
            property=raw["property"],
            value_type=raw["value_type"],
            closed=bool(raw["closed"]),
            description=" ".join(raw["description"].split()),
            values=list(raw.get("values") or []),
            sparse=bool(raw.get("sparse", False)),
            expected_empty=bool(raw.get("expected_empty", False)),
            min_value=raw.get("min_value"),
            max_value=raw.get("max_value"),
            signal_count=raw.get("signal_count"),
            signals=list(raw.get("signals") or []),
        )
    return out


def vocabularies_hash(entries: list[VocabEntry]) -> str:
    """Stable sha256 over the vocabulary set. Order-independent."""
    payload = sorted(
        json.dumps({
            "id": e.id, "value_type": e.value_type, "closed": e.closed,
            "values": sorted(e.values), "sparse": e.sparse,
            "expected_empty": e.expected_empty, "min_value": e.min_value,
            "max_value": e.max_value, "signal_count": e.signal_count,
            "signals": sorted(e.signals),
        }, sort_keys=True)
        for e in entries
    )
    digest = hashlib.sha256("\n".join(payload).encode()).hexdigest()
    return f"sha256:{digest}"


class _Vocab:
    """Lazy singleton so importing an adapter does not read the YAML."""

    def __init__(self, path: str | Path = DEFAULT_PATH):
        self._path = path
        self._entries: dict[str, VocabEntry] | None = None

    def _load(self) -> dict[str, VocabEntry]:
        if self._entries is None:
            self._entries = load_vocabularies(self._path)
        return self._entries

    def entries(self) -> list[VocabEntry]:
        return list(self._load().values())

    def get(self, applies_to: str, prop: str) -> VocabEntry | None:
        return self._load().get(f"{applies_to}.{prop}")

    def check(self, applies_to: str, prop: str, value: str) -> str:
        """Return *value*, raising if it is not declared for a closed vocabulary.

        Undeclared (applies_to, property) pairs pass through — the contract is
        seeded incrementally and this must not block unrelated adapters.
        """
        entry = self.get(applies_to, prop)
        if entry is None or not entry.closed:
            return value
        if value not in entry.values:
            raise ValueError(
                f"'{value}' is not a declared value of {entry.id}. "
                f"Declared: {sorted(entry.values)}. Add it to "
                f"config/controlled_vocabularies.yaml, or fix the emitter."
            )
        return value


VOCAB = _Vocab()
