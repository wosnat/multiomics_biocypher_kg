"""SignalP signal-peptide-type ontology adapter — flat + scored + 1:1.

Yields:
- SignalPeptideType nodes (the real SignalP 6.0 types observed across configured
  strains; the "OTHER" no-signal sentinel is never emitted).
- Gene_has_signal_peptide_type edges (Gene -> node), one per gene with a confident
  signal-peptide call, carrying the winning-class ``probability`` plus the
  ``cleavage_site`` / ``cleavage_probability`` as edge properties.

Structure mirrors psortb_adapter (the canonical scored-flat-ontology skeleton):
FLAT (no hierarchy / no <x>_is_a_<x> edges), 1:1 (<=1 edge per gene), and SCORED
(the edge carries numeric properties — modeled on Changes_expression_of). The
merged fields are the parallel scalars (signalp_type / signalp_probability /
signalp_cleavage_site / signalp_cleavage_probability) written by
build_gene_annotations.load_signalp.
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Iterator

from multiomics_kg.utils.curie_utils import normalize_curie
from multiomics_kg.utils.signalp import SIGNALP_VOCAB, display_name, is_kept

logger = logging.getLogger(__name__)


def _clean_str(value: str | None) -> str:
    if value is None:
        return ""
    return value.replace("'", "^").replace("|", "")


def _gene_node_id(locus_tag: str) -> str:
    return normalize_curie(f"ncbigene:{locus_tag}") or f"ncbigene_{locus_tag}"


def _node_id(raw_call: str) -> str:
    # "signalp" is not a bioregistry prefix, so normalize_curie returns None and
    # this falls back to "signalp_SP" (underscore), consistent with the
    # psortb/cazy/tcdb id-fallback pattern.
    return normalize_curie(f"signalp:{raw_call}") or f"signalp_{raw_call}"


class SignalPeptideAdapter:
    """Per-strain adapter: yields scored Gene_has_signal_peptide_type edges."""

    def __init__(self, genome_dir: Path, test_mode: bool = False) -> None:
        self.genome_dir = Path(genome_dir)
        self.test_mode = test_mode
        self._genes: dict = {}
        self._load()

    def _load(self) -> None:
        from multiomics_kg.utils.annotations_cache import load_merged_annotations
        self._genes = load_merged_annotations(self.genome_dir)

    def get_all_types(self) -> set[str]:
        """Distinct real signal-peptide type calls observed in this strain."""
        ids: set[str] = set()
        for gene in self._genes.values():
            call = gene.get("signalp_type")
            if is_kept(call):
                ids.add(call)
        return ids

    def get_edges(self):
        count = 0
        for locus_tag, gene in self._genes.items():
            call = gene.get("signalp_type")
            if not is_kept(call):          # skips None + "OTHER" sentinel
                continue
            props: dict = {}
            prob = gene.get("signalp_probability")
            if prob is not None:
                props["probability"] = float(prob)
            cs_site = gene.get("signalp_cleavage_site")
            if cs_site is not None:
                props["cleavage_site"] = int(cs_site)
            cs_prob = gene.get("signalp_cleavage_probability")
            if cs_prob is not None:
                props["cleavage_probability"] = float(cs_prob)
            yield (
                f"{locus_tag}-has_signalp-{call}",
                _gene_node_id(locus_tag),
                _node_id(call),
                "gene_has_signal_peptide_type",
                props,
            )
            count += 1
            if self.test_mode and count >= 100:
                return
        logger.debug(
            f"SignalPeptideAdapter({self.genome_dir.name}): "
            f"yielded {count} gene->signal-peptide-type edges"
        )


class MultiSignalPeptideAdapter:
    """Multi-strain orchestrator: owns the SignalPeptideType nodes (flat, level 0)."""

    def __init__(self, genome_config_file: str, test_mode: bool = False) -> None:
        self.test_mode = test_mode
        self._strain_adapters: list[SignalPeptideAdapter] = []
        self._build_strain_adapters(genome_config_file)

    def _build_strain_adapters(self, genome_config_file: str) -> None:
        with open(genome_config_file, newline="", encoding="utf-8") as fh:
            lines = [line for line in fh if not line.lstrip().startswith("#")]
        for row in csv.DictReader(lines):
            data_dir = (row.get("data_dir") or "").strip()
            if not data_dir:
                continue
            self._strain_adapters.append(
                SignalPeptideAdapter(
                    genome_dir=Path(data_dir), test_mode=self.test_mode
                )
            )
        logger.info(
            f"MultiSignalPeptideAdapter: loaded "
            f"{len(self._strain_adapters)} strain adapters"
        )

    def download_data(self, cache: bool = True) -> None:
        """No external resources — the vocabulary is in-process (SIGNALP_VOCAB)."""
        return

    def _all_observed_types(self) -> set[str]:
        ids: set[str] = set()
        for adapter in self._strain_adapters:
            ids |= adapter.get_all_types()
        return ids

    def get_nodes(self) -> Iterator[tuple[str, str, dict]]:
        observed = self._all_observed_types()
        for raw_call in sorted(observed):
            props = {
                "name": _clean_str(display_name(raw_call)),
                "signalp_id": raw_call,
                "level": 0,        # flat ontology
            }
            yield _node_id(raw_call), "signal peptide type", props
        logger.info(
            f"MultiSignalPeptideAdapter.get_nodes: "
            f"{len(observed)} SignalPeptideType nodes"
        )

    def get_edges(self):
        # FLAT: no <x>_is_a_<x> hierarchy edges — go straight to the gene edges.
        gene_count = 0
        for adapter in self._strain_adapters:
            for edge in adapter.get_edges():
                yield edge
                gene_count += 1
        logger.info(
            f"MultiSignalPeptideAdapter.get_edges: {gene_count} gene edges"
        )

    # Defensive: VOCAB import keeps the dependency explicit for readers / linters.
    _VOCAB = SIGNALP_VOCAB
