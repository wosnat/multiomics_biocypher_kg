"""TCDB ontology adapter.

Yields:
- TcdbFamily nodes (only IDs in cache/data/tcdb/tcdb_pruned.json's kept_tcdb_ids).
- Tcdb_family_is_a_tcdb_family parent edges within the pruned hierarchy.
- Gene_has_tcdb_family edges from per-strain gene_annotations_merged.json
  (filtered to pruned IDs to avoid dangling edges).
- Tcdb_family_transports_metabolite edges from tcdb_pruned.json["subtree_substrates"].
  Each kept TC ID maps to the union of substrate primary IDs from every
  tc_specificity descendant in the FULL TCDB hierarchy (rollup computed
  pre-pruning in step 6). The adapter emits one edge per (tc_id, primary_id)
  pair without traversing the hierarchy. Filter to source.level_kind =
  'tc_specificity' to recover leaf-only edges.

Two-class shape mirrors functional_annotation_adapter.MultiPfamAnnotationAdapter.
"""
from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Iterator

from multiomics_kg.utils.curie_utils import normalize_curie

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
        from multiomics_kg.utils.annotations_cache import load_merged_annotations
        self._genes = load_merged_annotations(self.genome_dir)

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
        self._subtree_substrates: dict[str, list[str]] = {}
        self._seed_aliases: dict[str, str] = {}

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
        # Pre-rolled-up substrate primary IDs per kept TC ID. Step 6 computes
        # the rollup over the FULL hierarchy (pre-pruning) so ancestor nodes
        # capture every TCDB descendant's substrates, not just gene-annotated
        # ones. Empty rollups are omitted from the file.
        self._subtree_substrates = pruned.get("subtree_substrates", {}) or {}
        # Remaps gene-annotated TCIDs not in the curated hierarchy to the
        # nearest curated ancestor (e.g. retired `3.A.1.35` → family `3.A.1`).
        # Older pruned files may not carry this field.
        self._seed_aliases = pruned.get("seed_aliases", {}) or {}

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

        # 2. Gene→TcdbFamily edges. Remap edges whose target isn't in the
        # curated hierarchy through `seed_aliases` (built by step 6) so retired
        # eggNOG TCIDs anchor onto the nearest curated ancestor instead of
        # being silently dropped.
        gene_count = 0
        remapped = 0
        dropped = 0
        kept_node_ids = {_tcdb_node_id(tc) for tc in self._kept_ids}
        alias_node_ids = {
            _tcdb_node_id(orig): _tcdb_node_id(anchor)
            for orig, anchor in self._seed_aliases.items()
        }
        for adapter in self._strain_adapters:
            for edge in adapter.get_edges():
                target = edge[2]
                if target in kept_node_ids:
                    yield edge
                    gene_count += 1
                    continue
                anchor_target = alias_node_ids.get(target)
                if anchor_target and anchor_target in kept_node_ids:
                    edge_id, source, _, label, props = edge
                    yield (edge_id, source, anchor_target, label, props)
                    gene_count += 1
                    remapped += 1
                else:
                    dropped += 1
        if remapped:
            logger.info(
                f"MultiTcdbAnnotationAdapter: re-anchored {remapped} gene→TCDB edges via seed_aliases"
            )
        if dropped:
            logger.warning(
                f"MultiTcdbAnnotationAdapter: dropped {dropped} gene→TCDB edges to unpruned IDs"
            )

        # 3. Substrate edges. Pre-rolled-up at step-6 time over the FULL
        # hierarchy, so ancestors carry substrates from every TCDB descendant
        # (not just gene-annotated ones). Recover leaf-only semantics by
        # filtering to source.level_kind = 'tc_specificity'.
        sub_count = 0
        for tcdb_id, primary_ids in sorted(self._subtree_substrates.items()):
            for primary in primary_ids:
                yield (
                    f"{tcdb_id}-transports-{primary}",
                    _tcdb_node_id(tcdb_id),
                    primary,
                    "tcdb_family_transports_metabolite",
                    {},
                )
                sub_count += 1

        logger.info(
            f"MultiTcdbAnnotationAdapter.get_edges: {parent_count} parent, "
            f"{gene_count} gene, {sub_count} substrate edges"
        )
