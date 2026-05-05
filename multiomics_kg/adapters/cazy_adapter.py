"""CAZy ontology adapter — observed-only (no external download).

Yields:
- CazyFamily nodes (class + family + subfamily, only IDs observed in
  gene_annotations_merged.json across configured strains).
- Cazy_family_is_a_cazy_family parent edges (subfamily → family, family → class).
- Gene_has_cazy_family edges from per-strain merged JSON. Each gene attaches
  to the MOST SPECIFIC observed level (subfamily if present, else family).

Two-class shape mirrors functional_annotation_adapter.MultiPfamAnnotationAdapter
and adapters.tcdb_adapter (commit 5).
"""
from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Iterator

from multiomics_kg.utils.curie_utils import normalize_curie

from multiomics_kg.utils.cazy_utils import CAZY_CLASSES, parse_cazy_id

logger = logging.getLogger(__name__)


def _clean_str(value: str | None) -> str:
    if value is None:
        return ""
    return value.replace("'", "^").replace("|", "")


def _gene_node_id(locus_tag: str) -> str:
    return normalize_curie(f"ncbigene:{locus_tag}") or f"ncbigene_{locus_tag}"


def _cazy_node_id(cazy_id: str) -> str:
    return normalize_curie(f"cazy:{cazy_id}") or f"cazy_{cazy_id}"


def _most_specific_id(token: str) -> str | None:
    """Return the subfamily if present, else the family. None for malformed."""
    parsed = parse_cazy_id(token)
    if parsed is None:
        return None
    family, subfamily = parsed
    return subfamily or family


def _classify(cazy_id: str) -> tuple[int, str]:
    """Return (level, level_kind) for a class / family / subfamily ID."""
    if cazy_id in CAZY_CLASSES:
        return 0, "cazy_class"
    parsed = parse_cazy_id(cazy_id)
    if parsed is None:
        # Should be unreachable — the orchestrator only sees parse-valid IDs.
        return 1, "cazy_family"
    _family, subfamily = parsed
    if subfamily and cazy_id == subfamily:
        return 2, "cazy_subfamily"
    return 1, "cazy_family"


class CazyAnnotationAdapter:
    """Per-strain adapter: yields Gene_has_cazy_family edges from merged JSON."""

    def __init__(self, genome_dir: Path, test_mode: bool = False) -> None:
        self.genome_dir = Path(genome_dir)
        self.test_mode = test_mode
        self._genes: dict = {}
        self._load()

    def _load(self) -> None:
        json_path = self.genome_dir / "gene_annotations_merged.json"
        if not json_path.exists():
            logger.warning(f"gene_annotations_merged.json not found at {json_path}, skipping")
            return
        with open(json_path, encoding="utf-8") as fh:
            self._genes = json.load(fh)

    def get_all_cazy_ids(self) -> set[str]:
        ids: set[str] = set()
        for gene in self._genes.values():
            for token in gene.get("cazy_ids") or []:
                if not token:
                    continue
                parsed = parse_cazy_id(token)
                if parsed is None:
                    logger.debug(f"CAZy malformed token {token!r} (strain {self.genome_dir.name}), skipping")
                    continue
                family, subfamily = parsed
                ids.add(family)
                if subfamily:
                    ids.add(subfamily)
        return ids

    def get_edges(self):
        count = 0
        for locus_tag, gene in self._genes.items():
            for token in gene.get("cazy_ids") or []:
                target = _most_specific_id(token)
                if target is None:
                    logger.debug(
                        f"CAZy malformed token {token!r} on {locus_tag} "
                        f"(strain {self.genome_dir.name}), skipping"
                    )
                    continue
                yield (
                    f"{locus_tag}-has_cazy-{target}",
                    _gene_node_id(locus_tag),
                    _cazy_node_id(target),
                    "gene_has_cazy_family",
                    {},
                )
                count += 1
                if self.test_mode and count >= 100:
                    return
        logger.debug(
            f"CazyAnnotationAdapter({self.genome_dir.name}): yielded {count} gene→CAZy edges"
        )


class MultiCazyAnnotationAdapter:
    """Multi-strain orchestrator: owns CazyFamily nodes + parent edges."""

    def __init__(
        self,
        genome_config_file: str,
        test_mode: bool = False,
    ) -> None:
        self.test_mode = test_mode
        self._strain_adapters: list[CazyAnnotationAdapter] = []
        self._build_strain_adapters(genome_config_file)

    def _build_strain_adapters(self, genome_config_file: str) -> None:
        with open(genome_config_file, newline="", encoding="utf-8") as fh:
            lines = [line for line in fh if not line.lstrip().startswith("#")]
        reader = csv.DictReader(lines)
        for row in reader:
            data_dir = (row.get("data_dir") or "").strip()
            if not data_dir:
                continue
            self._strain_adapters.append(
                CazyAnnotationAdapter(genome_dir=Path(data_dir), test_mode=self.test_mode)
            )
        logger.info(
            f"MultiCazyAnnotationAdapter: loaded {len(self._strain_adapters)} strain adapters"
        )

    def download_data(self, cache: bool = True) -> None:
        """No external resources — CAZy hierarchy is in-process Python."""
        return

    def _all_observed_ids(self) -> set[str]:
        ids: set[str] = set()
        for adapter in self._strain_adapters:
            ids |= adapter.get_all_cazy_ids()
        return ids

    def _expand_to_class_family_subfamily(self, observed: set[str]) -> set[str]:
        """For every observed family/subfamily, add its class (and family if subfamily)."""
        out: set[str] = set()
        for tok in observed:
            parsed = parse_cazy_id(tok)
            if parsed is None:
                continue
            family, subfamily = parsed
            # Class code is the alphabetic prefix of the family
            for code in CAZY_CLASSES:
                if family.startswith(code) and family[len(code):].isdigit():
                    out.add(code)
                    break
            out.add(family)
            if subfamily:
                out.add(subfamily)
        return out

    def get_nodes(self) -> Iterator[tuple[str, str, dict]]:
        observed = self._all_observed_ids()
        all_ids = self._expand_to_class_family_subfamily(observed)
        for cazy_id in sorted(all_ids):
            level, level_kind = _classify(cazy_id)
            if cazy_id in CAZY_CLASSES:
                name = CAZY_CLASSES[cazy_id]
            else:
                name = cazy_id  # fallback
            props = {
                "name": _clean_str(name),
                "cazy_id": cazy_id,
                "level": level,
                "level_kind": level_kind,
            }
            yield _cazy_node_id(cazy_id), "cazy family", props
        logger.info(f"MultiCazyAnnotationAdapter.get_nodes: {len(all_ids)} CazyFamily nodes")

    def get_edges(self):
        observed = self._all_observed_ids()
        all_ids = self._expand_to_class_family_subfamily(observed)

        # 1. Parent edges: subfamily → family, family → class
        parent_count = 0
        for cazy_id in sorted(all_ids):
            if cazy_id in CAZY_CLASSES:
                continue
            parsed = parse_cazy_id(cazy_id)
            if parsed is None:
                continue
            family, subfamily = parsed
            cls_code = None
            for code in CAZY_CLASSES:
                if family.startswith(code) and family[len(code):].isdigit():
                    cls_code = code
                    break
            if cls_code is None:
                continue
            if subfamily and cazy_id == subfamily:
                parent = family
            else:
                parent = cls_code
            yield (
                f"{cazy_id}-parent-{parent}",
                _cazy_node_id(cazy_id),
                _cazy_node_id(parent),
                "cazy_family_is_a_cazy_family",
                {},
            )
            parent_count += 1

        # 2. Gene→CAZy edges via per-strain delegation. Every observed ID is a
        # node by construction (orchestrator's _expand_to_class_family_subfamily
        # covers all observed family + subfamily IDs and their classes), so no
        # filtering is needed — unlike the TCDB orchestrator, which has to drop
        # edges to unpruned IDs.
        gene_count = 0
        for adapter in self._strain_adapters:
            for edge in adapter.get_edges():
                yield edge
                gene_count += 1
        logger.info(
            f"MultiCazyAnnotationAdapter.get_edges: {parent_count} parent, {gene_count} gene edges"
        )
