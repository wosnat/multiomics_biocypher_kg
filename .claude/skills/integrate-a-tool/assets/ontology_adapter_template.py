"""<Tool> ontology adapter — TEMPLATE (copy to multiomics_kg/adapters/<tool>_adapter.py).

Modeled verbatim on cazy_adapter.py — the canonical skeleton. Two classes:
- <Tool>Adapter        : per-strain, reads gene_annotations_merged.json, yields Gene_has_<x> edges.
- Multi<Tool>Adapter   : multi-strain orchestrator, OWNS the <NodeLabel> nodes (+ hierarchy edges).

Replace every <…> placeholder. Decide three things from references/decision-tree.md before editing:
  1. FLAT vs HIERARCHICAL  → _classify() + whether to emit <x>_is_a_<x> parent edges.
  2. SCORED vs BARE        → whether the Gene_has_<x> edge carries a properties dict.
  3. MULTI-CALL vs 1:1     → fan out one edge per list element, or emit ≤1 edge per gene.

The two immediate targets (psortb, signalp) are FLAT + SCORED + 1:1 + skip-the-sentinel — see the
"1:1 SCALAR variant" block in <Tool>Adapter.get_edges below. The default body shown is the general
MULTI-CALL + SCORED shape (closest to a scored cazy).
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Iterator

from multiomics_kg.utils.curie_utils import normalize_curie

logger = logging.getLogger(__name__)

# Controlled vocabulary → human-readable display name. Copy the distinct calls.json values here
# (skip no-signal sentinels like psortb "Unknown" / signalp "Other" — absence of an edge encodes them).
VOCAB: dict[str, str] = {
    # "<RAW_CALL>": "<Human readable name>",
    # "OuterMembrane": "Outer membrane",
}


def _clean_str(value: str | None) -> str:
    if value is None:
        return ""
    return value.replace("'", "^").replace("|", "")


def _gene_node_id(locus_tag: str) -> str:
    return normalize_curie(f"ncbigene:{locus_tag}") or f"ncbigene_{locus_tag}"


def _node_id(raw_id: str) -> str:
    # normalize_curie returns the colon CURIE only for a REGISTERED bioregistry prefix; otherwise it
    # returns None and the `or` fallback yields the UNDERSCORE form. An unregistered tool prefix gives
    # "psortb_OuterMembrane" (underscore), NOT "psortb:OuterMembrane". Register a bioregistry prefix to
    # earn the colon; otherwise write your unit-test assertions against the underscore id.
    return normalize_curie(f"<prefix>:{raw_id}") or f"<prefix>_{raw_id}"


def _classify(raw_id: str) -> tuple[int, str]:
    """Return (level, level_kind). FLAT ontology → always (0, "<tool>_class").

    HIERARCHICAL: derive depth from the id structure (see cazy_adapter._classify /
    tcdb level_kind from the pre-built cache). level 0 = broadest/root, mandatory.
    level_kind is OPTIONAL on a flat ontology — drop it from the node props if so.
    """
    return 0, "<tool>_class"


class <Tool>Adapter:
    """Per-strain adapter: yields Gene_has_<x> edges from gene_annotations_merged.json."""

    def __init__(self, genome_dir: Path, test_mode: bool = False) -> None:
        self.genome_dir = Path(genome_dir)
        self.test_mode = test_mode
        self._genes: dict = {}
        self._load()

    def _load(self) -> None:
        from multiomics_kg.utils.annotations_cache import load_merged_annotations
        self._genes = load_merged_annotations(self.genome_dir)

    def get_all_ids(self) -> set[str]:
        """Distinct raw calls observed in this strain (for the orchestrator's node set)."""
        ids: set[str] = set()
        for gene in self._genes.values():
            # ── MULTI-CALL (default): list-of-dicts merged field ──
            for rec in gene.get("<tool>_calls") or []:
                call = rec.get("call") if isinstance(rec, dict) else rec
                if call and call in VOCAB:
                    ids.add(call)
            # ── 1:1 SCALAR variant (psortb/signalp): replace the loop above with ──
            # call = gene.get("<tool>_call")          # scalar str field, e.g. psortb_localization
            # if call and call in VOCAB:
            #     ids.add(call)
        return ids

    def get_edges(self):
        count = 0
        for locus_tag, gene in self._genes.items():
            # ── MULTI-CALL + SCORED (default) ────────────────────────────────────
            for rec in gene.get("<tool>_calls") or []:
                call = rec.get("call") if isinstance(rec, dict) else rec
                if not call or call not in VOCAB:               # skips sentinels too
                    continue
                props: dict = {}
                score = rec.get("score") if isinstance(rec, dict) else None
                if score is not None:
                    props["score"] = float(score)               # numeric — NOT _clean_str'd
                yield (
                    f"{locus_tag}-has_<tool>-{call}",
                    _gene_node_id(locus_tag),
                    _node_id(call),
                    "gene_has_<x>",
                    props,                                       # {} if BARE (no score)
                )
                count += 1
                if self.test_mode and count >= 100:
                    return

            # ── 1:1 SCALAR variant (psortb / signalp): replace the loop above with ──
            # call = gene.get("<tool>_call")                     # scalar str field
            # if call and call in VOCAB:                         # skip Unknown/Other sentinel
            #     props = {"score": float(gene["<tool>_score"])} # or {"probability": ...}
            #     yield (f"{locus_tag}-has_<tool>-{call}", _gene_node_id(locus_tag),
            #            _node_id(call), "gene_has_<x>", props)
            #     count += 1
        logger.debug(f"<Tool>Adapter({self.genome_dir.name}): yielded {count} gene→<X> edges")


class Multi<Tool>Adapter:
    """Multi-strain orchestrator: owns <NodeLabel> nodes (+ parent edges if hierarchical)."""

    def __init__(self, genome_config_file: str, test_mode: bool = False) -> None:
        self.test_mode = test_mode
        self._strain_adapters: list[<Tool>Adapter] = []
        self._build_strain_adapters(genome_config_file)

    def _build_strain_adapters(self, genome_config_file: str) -> None:
        with open(genome_config_file, newline="", encoding="utf-8") as fh:
            lines = [line for line in fh if not line.lstrip().startswith("#")]
        for row in csv.DictReader(lines):
            data_dir = (row.get("data_dir") or "").strip()
            if not data_dir:
                continue
            self._strain_adapters.append(
                <Tool>Adapter(genome_dir=Path(data_dir), test_mode=self.test_mode)
            )
        logger.info(f"Multi<Tool>Adapter: loaded {len(self._strain_adapters)} strain adapters")

    def download_data(self, cache: bool = True) -> None:
        """No external resources — the vocabulary is in-process (VOCAB). Keep for adapter parity."""
        return

    def _all_observed_ids(self) -> set[str]:
        ids: set[str] = set()
        for adapter in self._strain_adapters:
            ids |= adapter.get_all_ids()
        return ids

    def get_nodes(self) -> Iterator[tuple[str, str, dict]]:
        # FLAT: emit one node per observed vocab member. HIERARCHICAL: expand to ancestors
        # first (see cazy_adapter._expand_to_class_family_subfamily).
        for raw_id in sorted(self._all_observed_ids()):
            level, level_kind = _classify(raw_id)
            props = {
                "name": _clean_str(VOCAB.get(raw_id, raw_id)),
                "<tool>_id": raw_id,
                "level": level,
                # "level_kind": level_kind,   # optional on a flat ontology — drop if unused
            }
            yield _node_id(raw_id), "<node label>", props        # label_in_input from schema_config
        logger.info(f"Multi<Tool>Adapter.get_nodes: {len(self._all_observed_ids())} <NodeLabel> nodes")

    def get_edges(self):
        # HIERARCHICAL only: emit <x>_is_a_<x> child→parent edges here (property-less {}),
        # exactly like cazy_adapter.get_edges step 1. FLAT: skip straight to the gene edges.
        gene_count = 0
        for adapter in self._strain_adapters:
            for edge in adapter.get_edges():
                yield edge
                gene_count += 1
        logger.info(f"Multi<Tool>Adapter.get_edges: {gene_count} gene edges")
