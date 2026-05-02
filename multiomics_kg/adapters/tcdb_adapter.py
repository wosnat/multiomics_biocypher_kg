"""TCDB ontology adapter.

Yields:
- TcdbFamily nodes (only IDs in cache/data/tcdb/tcdb_pruned.json's kept_tcdb_ids).
- Tcdb_family_is_a_tcdb_family parent edges within the pruned hierarchy.
- Gene_has_tcdb_family edges from per-strain gene_annotations_merged.json
  (filtered to pruned IDs to avoid dangling edges).
- Tcdb_family_transports_metabolite edges on tc_specificity-level nodes only,
  using pre-resolved metabolite primary IDs from tcdb_pruned.json.

Two-class shape mirrors functional_annotation_adapter.MultiPfamAnnotationAdapter.
"""
from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Iterator

from bioregistry import normalize_curie

logger = logging.getLogger(__name__)

_TC_CLASS_NAMES = {
    "1": "Channels and Pores",
    "2": "Electrochemical Potential-driven Transporters",
    "3": "Primary Active Transporters",
    "4": "Group Translocators",
    "5": "Transmembrane Electron Carriers",
    "8": "Auxiliary Transport Proteins",
    "9": "Incompletely Characterized Transport Systems",
}


def _clean_str(value: str | None) -> str:
    if value is None:
        return ""
    return value.replace("'", "^").replace("|", "")


def _gene_node_id(locus_tag: str) -> str:
    """Normalize a locus_tag to the gene node ID format (matches Pfam adapter)."""
    return normalize_curie(f"ncbigene:{locus_tag}") or f"ncbigene_{locus_tag}"


def _tcdb_node_id(tcdb_id: str) -> str:
    return normalize_curie(f"tcdb:{tcdb_id}") or f"tcdb_{tcdb_id}"


class TcdbAnnotationAdapter:
    """Per-strain adapter: yields Gene_has_tcdb_family edges from gene_annotations_merged.json."""

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

    def get_all_tcdb_ids(self) -> set[str]:
        ids: set[str] = set()
        for gene in self._genes.values():
            for tc in gene.get("transporter_classification") or []:
                if tc:
                    ids.add(tc)
        return ids

    def get_edges(self):
        count = 0
        for locus_tag, gene in self._genes.items():
            for tc in gene.get("transporter_classification") or []:
                if not tc:
                    continue
                yield (
                    f"{locus_tag}-has_tcdb-{tc}",
                    _gene_node_id(locus_tag),
                    _tcdb_node_id(tc),
                    "gene_has_tcdb_family",
                    {},
                )
                count += 1
                if self.test_mode and count >= 100:
                    return
        logger.debug(
            f"TcdbAnnotationAdapter({self.genome_dir.name}): yielded {count} gene→TCDB edges"
        )


class MultiTcdbAnnotationAdapter:
    """Multi-strain orchestrator: owns TcdbFamily nodes + parent edges + substrate edges."""

    def __init__(
        self,
        genome_config_file: str,
        cache_root: Path,
        test_mode: bool = False,
        cache: bool = True,
    ) -> None:
        self.cache_root = Path(cache_root)
        self.test_mode = test_mode
        self.cache = cache
        self._strain_adapters: list[TcdbAnnotationAdapter] = []
        self._build_strain_adapters(genome_config_file)
        self._hierarchy: dict[str, dict] = {}
        self._kept_ids: set[str] = set()
        self._leaf_substrates: dict[str, list[str]] = {}

    def _build_strain_adapters(self, genome_config_file: str) -> None:
        with open(genome_config_file, newline="", encoding="utf-8") as fh:
            lines = [line for line in fh if not line.lstrip().startswith("#")]
        reader = csv.DictReader(lines)
        for row in reader:
            data_dir = (row.get("data_dir") or "").strip()
            if not data_dir:
                continue
            self._strain_adapters.append(
                TcdbAnnotationAdapter(genome_dir=Path(data_dir), test_mode=self.test_mode)
            )
        logger.info(
            f"MultiTcdbAnnotationAdapter: loaded {len(self._strain_adapters)} strain adapters"
        )

    def download_data(self, cache: bool = True) -> None:
        """Read tcdb_hierarchy.json + tcdb_pruned.json (both built by step 6)."""
        tcdb_dir = self.cache_root / "tcdb"
        hierarchy_path = tcdb_dir / "tcdb_hierarchy.json"
        pruned_path = tcdb_dir / "tcdb_pruned.json"
        if not hierarchy_path.exists() or not pruned_path.exists():
            raise FileNotFoundError(
                f"Missing TCDB cache file(s): {hierarchy_path} and/or {pruned_path}. "
                f"Run `bash scripts/prepare_data.sh --steps 6 --force` first."
            )
        self._hierarchy = json.loads(hierarchy_path.read_text())
        pruned = json.loads(pruned_path.read_text())
        self._kept_ids = set(pruned["kept_tcdb_ids"])
        self._leaf_substrates = pruned["leaf_substrates"]

    def get_nodes(self) -> Iterator[tuple[str, str, dict]]:
        if not self._kept_ids:
            self.download_data(cache=self.cache)

        emit_count = 0
        for tcdb_id in sorted(self._kept_ids):
            if self.test_mode and emit_count >= 100:
                break
            entry = self._hierarchy.get(tcdb_id, {})
            level = entry.get("level", 0)
            level_kind = entry.get("level_kind", "tc_class")
            raw_name = entry.get("name") or ""
            # Class fallback: source name often empty; pull from _TC_CLASS_NAMES
            if not raw_name and level_kind == "tc_class":
                raw_name = _TC_CLASS_NAMES.get(tcdb_id, "")
            # Other levels: fall back to the tcdb_id itself
            display_name = raw_name or tcdb_id

            props = {
                "name": _clean_str(display_name),
                "tcdb_id": tcdb_id,
                "level": level,
                "level_kind": level_kind,
            }
            if entry.get("superfamily"):
                props["superfamily"] = _clean_str(entry["superfamily"])
            yield _tcdb_node_id(tcdb_id), "tcdb family", props
            emit_count += 1
        logger.info(f"MultiTcdbAnnotationAdapter.get_nodes: {emit_count} TcdbFamily nodes")

    def get_edges(self):
        if not self._kept_ids:
            self.download_data(cache=self.cache)

        # 1. Parent edges within the pruned hierarchy
        parent_count = 0
        for tcdb_id in sorted(self._kept_ids):
            entry = self._hierarchy.get(tcdb_id, {})
            parent = entry.get("parent")
            if parent and parent in self._kept_ids:
                yield (
                    f"{tcdb_id}-parent-{parent}",
                    _tcdb_node_id(tcdb_id),
                    _tcdb_node_id(parent),
                    "tcdb_family_is_a_tcdb_family",
                    {},
                )
                parent_count += 1

        # 2. Gene→TcdbFamily edges (delegate to per-strain adapters; filter to kept IDs)
        gene_count = 0
        dropped = 0
        kept_node_ids = {_tcdb_node_id(tc) for tc in self._kept_ids}
        for adapter in self._strain_adapters:
            for edge in adapter.get_edges():
                if edge[2] not in kept_node_ids:
                    dropped += 1
                    continue
                yield edge
                gene_count += 1
        if dropped:
            logger.debug(
                f"MultiTcdbAnnotationAdapter: dropped {dropped} gene→TCDB edges to unpruned IDs"
            )

        # 3. Substrate edges (leaves only)
        sub_count = 0
        for leaf, primary_ids in sorted(self._leaf_substrates.items()):
            for primary in primary_ids:
                yield (
                    f"{leaf}-transports-{primary}",
                    _tcdb_node_id(leaf),
                    primary,
                    "tcdb_family_transports_metabolite",
                    {},
                )
                sub_count += 1

        logger.info(
            f"MultiTcdbAnnotationAdapter.get_edges: {parent_count} parent, "
            f"{gene_count} gene, {sub_count} substrate edges"
        )
