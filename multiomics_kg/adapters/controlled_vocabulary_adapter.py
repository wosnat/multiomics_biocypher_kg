"""ControlledVocabulary adapter — publishes config/controlled_vocabularies.yaml
as queryable nodes so consumers stop hard-coding value sets.

Node-only, modelled on data_source_adapter.py. Booleans are emitted as
"true"/"false" strings per house rule R5.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterator

from multiomics_kg.utils.controlled_vocab import DEFAULT_PATH, load_vocabularies


def _clean_str(value: str) -> str:
    """Sanitize string for BioCypher CSV output (avoid ' and |)."""
    return value.replace("'", "^").replace("|", "")


def _b(value: bool) -> str:
    """R5: adapter-emitted booleans are meaningful strings, never native bool."""
    return "true" if value else "false"


class ControlledVocabularyAdapter:
    def __init__(self, config_path: str | Path = DEFAULT_PATH):
        self.config_path = Path(config_path)
        self._entries: list[Any] = []

    def download_data(self, cache: bool = True) -> None:
        self._entries = list(load_vocabularies(self.config_path).values())

    def entries(self) -> list[Any]:
        """Public accessor for the loaded VocabEntry objects (e.g. for hashing)."""
        if not self._entries:
            self.download_data()
        return self._entries

    def get_nodes(self) -> Iterator[tuple[str, str, dict[str, Any]]]:
        if not self._entries:
            self.download_data()
        for e in self._entries:
            props: dict[str, Any] = {
                "id": _clean_str(e.id),
                "applies_to": _clean_str(e.applies_to),
                "applies_to_kind": _clean_str(e.applies_to_kind),
                "property": _clean_str(e.property),
                "value_type": _clean_str(e.value_type),
                "closed": _b(e.closed),
                "sparse": _b(e.sparse),
                "expected_empty": _b(e.expected_empty),
                "exhaustive": _b(e.exhaustive),
                "values": [_clean_str(v) for v in e.values],
                "description": _clean_str(e.description),
            }
            # Sparse numeric/score metadata — omitted rather than nulled so the
            # absence is meaningful (a string vocabulary has no min_value).
            if e.min_value is not None:
                props["min_value"] = float(e.min_value)
            if e.max_value is not None:
                props["max_value"] = float(e.max_value)
            if e.signal_count is not None:
                props["signal_count"] = int(e.signal_count)
            if e.signals:
                props["signals"] = [_clean_str(s) for s in e.signals]
            yield e.id, "controlled_vocabulary", props

    def get_edges(self) -> Iterator[tuple]:
        return iter(())
