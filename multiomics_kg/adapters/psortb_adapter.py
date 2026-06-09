"""PSORTb subcellular-localization ontology adapter — flat + scored + 1:1.

Yields:
- SubcellularLocalization nodes (the 5 real PSORTb classes observed across
  configured strains; the "Unknown" sentinel is never emitted).
- Gene_has_subcellular_localization edges (Gene -> node), one per gene with a
  confident call, carrying the PSORTb confidence ``score`` as an edge property.

Structure mirrors cazy_adapter (the canonical skeleton) but is FLAT (no
hierarchy / no <x>_is_a_<x> edges), 1:1 (<=1 edge per gene), and SCORED (the
edge carries a ``score`` property — the first scored ontology edge, modeled on
Changes_expression_of). The merged field is the parallel scalar pair
(psortb_localization + psortb_score) written by build_gene_annotations.load_psortb.
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Iterator

from multiomics_kg.utils.curie_utils import normalize_curie
from multiomics_kg.utils.psortb import LOCALIZATION_VOCAB, display_name, is_kept

logger = logging.getLogger(__name__)


def _clean_str(value: str | None) -> str:
    if value is None:
        return ""
    return value.replace("'", "^").replace("|", "")


def _gene_node_id(locus_tag: str) -> str:
    return normalize_curie(f"ncbigene:{locus_tag}") or f"ncbigene_{locus_tag}"


def _node_id(raw_call: str) -> str:
    # "psortb" is not a bioregistry prefix, so normalize_curie returns None and
    # this falls back to "psortb_OuterMembrane" (underscore), consistent with the
    # cazy/tcdb id-fallback pattern.
    return normalize_curie(f"psortb:{raw_call}") or f"psortb_{raw_call}"


class SubcellularLocalizationAdapter:
    """Per-strain adapter: yields scored Gene_has_subcellular_localization edges."""

    def __init__(self, genome_dir: Path, test_mode: bool = False) -> None:
        self.genome_dir = Path(genome_dir)
        self.test_mode = test_mode
        self._genes: dict = {}
        self._load()

    def _load(self) -> None:
        from multiomics_kg.utils.annotations_cache import load_merged_annotations
        self._genes = load_merged_annotations(self.genome_dir)

    def get_all_localizations(self) -> set[str]:
        """Distinct real localization calls observed in this strain."""
        ids: set[str] = set()
        for gene in self._genes.values():
            call = gene.get("psortb_localization")
            if is_kept(call):
                ids.add(call)
        return ids

    def get_edges(self):
        count = 0
        for locus_tag, gene in self._genes.items():
            call = gene.get("psortb_localization")
            if not is_kept(call):          # skips None + "Unknown" sentinel
                continue
            props: dict = {}
            score = gene.get("psortb_score")
            if score is not None:
                props["score"] = float(score)
            yield (
                f"{locus_tag}-has_psortb-{call}",
                _gene_node_id(locus_tag),
                _node_id(call),
                "gene_has_subcellular_localization",
                props,
            )
            count += 1
            if self.test_mode and count >= 100:
                return
        logger.debug(
            f"SubcellularLocalizationAdapter({self.genome_dir.name}): "
            f"yielded {count} gene->localization edges"
        )


class MultiSubcellularLocalizationAdapter:
    """Multi-strain orchestrator: owns the SubcellularLocalization nodes (flat, level 0)."""

    def __init__(self, genome_config_file: str, test_mode: bool = False) -> None:
        self.test_mode = test_mode
        self._strain_adapters: list[SubcellularLocalizationAdapter] = []
        self._build_strain_adapters(genome_config_file)

    def _build_strain_adapters(self, genome_config_file: str) -> None:
        with open(genome_config_file, newline="", encoding="utf-8") as fh:
            lines = [line for line in fh if not line.lstrip().startswith("#")]
        for row in csv.DictReader(lines):
            data_dir = (row.get("data_dir") or "").strip()
            if not data_dir:
                continue
            self._strain_adapters.append(
                SubcellularLocalizationAdapter(
                    genome_dir=Path(data_dir), test_mode=self.test_mode
                )
            )
        logger.info(
            f"MultiSubcellularLocalizationAdapter: loaded "
            f"{len(self._strain_adapters)} strain adapters"
        )

    def download_data(self, cache: bool = True) -> None:
        """No external resources — the vocabulary is in-process (LOCALIZATION_VOCAB)."""
        return

    def _all_observed_localizations(self) -> set[str]:
        ids: set[str] = set()
        for adapter in self._strain_adapters:
            ids |= adapter.get_all_localizations()
        return ids

    def get_nodes(self) -> Iterator[tuple[str, str, dict]]:
        observed = self._all_observed_localizations()
        for raw_call in sorted(observed):
            props = {
                "name": _clean_str(display_name(raw_call)),
                "psortb_id": raw_call,
                "level": 0,        # flat ontology
            }
            yield _node_id(raw_call), "subcellular localization", props
        logger.info(
            f"MultiSubcellularLocalizationAdapter.get_nodes: "
            f"{len(observed)} SubcellularLocalization nodes"
        )

    def get_edges(self):
        # FLAT: no <x>_is_a_<x> hierarchy edges — go straight to the gene edges.
        gene_count = 0
        for adapter in self._strain_adapters:
            for edge in adapter.get_edges():
                yield edge
                gene_count += 1
        logger.info(
            f"MultiSubcellularLocalizationAdapter.get_edges: {gene_count} gene edges"
        )

    # Defensive: VOCAB import keeps the dependency explicit for readers / linters.
    _VOCAB = LOCALIZATION_VOCAB
