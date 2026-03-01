"""
Functional annotation adapter: gene → GO edges and GO hierarchy subset.

Creates:
- GO nodes (BiologicalProcess, CellularComponent, MolecularFunction) for only the
  ~5K GO terms referenced in gene annotations plus their ancestors (not all 30K).
- gene_involved_in_biological_process edges
- gene_located_in_cellular_component edges
- gene_enables_molecular_function edges
- GO-GO hierarchy edges (is_a, part_of, etc.) for the ancestor closure only.

This adapter replaces the heavyweight pypath-based go_adapter.py for the standard
build. The existing --go flag path is preserved for backward compatibility but should
not be run simultaneously (GO node IDs would conflict).
"""

import csv
import json
import logging
from pathlib import Path

from bioregistry import normalize_curie

from multiomics_kg.utils.go_utils import (
    NAMESPACE_TO_LABEL,
    compute_ancestry_closure,
    load_go_data,
    make_go_go_edge_label,
)

logger = logging.getLogger(__name__)

# Maps GO namespace string → gene→GO edge label_in_input value
_NS_TO_GENE_EDGE_LABEL: dict[str, str] = {
    "biological_process": "gene_involved_in_biological_process",
    "cellular_component": "gene_located_in_cellular_component",
    "molecular_function": "gene_enables_molecular_function",
}


def _go_node_id(go_id: str) -> str:
    """Normalize a GO term ID to the format BioCypher uses as node ID.

    e.g. "GO:0005737" → "go:0005737"  (matches go_adapter.py output)
    """
    return normalize_curie(f"go:{go_id}") or f"go_{go_id}"


def _gene_node_id(locus_tag: str) -> str:
    """Normalize a locus_tag to the gene node ID format."""
    return normalize_curie(f"ncbigene:{locus_tag}") or f"ncbigene_{locus_tag}"


class GoAnnotationAdapter:
    """
    Per-strain adapter: reads gene_annotations_merged.json and yields gene→GO edges.

    Args:
        genome_dir: directory containing ``gene_annotations_merged.json``
        go_data: full GO data dict from :func:`~multiomics_kg.utils.go_utils.load_go_data`
        test_mode: if True, stop after 100 genes
    """

    def __init__(self, genome_dir: Path, go_data: dict, test_mode: bool = False) -> None:
        self.genome_dir = Path(genome_dir)
        self.go_data = go_data
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
        logger.debug(f"Loaded {len(self._genes)} genes from {json_path}")

    def get_all_go_ids(self) -> set[str]:
        """Return all GO term IDs directly referenced by any gene in this strain."""
        ids: set[str] = set()
        for gene in self._genes.values():
            for go_id in gene.get("go_terms") or []:
                if go_id:
                    ids.add(go_id)
        return ids

    def get_edges(self):
        """
        Yield gene→GO edges as 5-tuples::

            (edge_id, gene_node_id, go_node_id, edge_label, properties)

        Only edges where the GO term is present in *go_data* and has a recognized
        namespace are emitted.  GO terms not found in *go_data* (obsolete/unknown) are
        silently skipped.
        """
        count = 0
        for locus_tag, gene in self._genes.items():
            go_terms = gene.get("go_terms") or []
            for go_id in go_terms:
                if not go_id:
                    continue
                entry = self.go_data.get(go_id)
                if entry is None:
                    continue
                edge_label = _NS_TO_GENE_EDGE_LABEL.get(entry["namespace"])
                if edge_label is None:
                    continue

                edge_id = f"{locus_tag}-go-{go_id}"
                yield (
                    edge_id,
                    _gene_node_id(locus_tag),
                    _go_node_id(go_id),
                    edge_label,
                    {},
                )
                count += 1

            if self.test_mode and count >= 100:
                break

        logger.debug(f"GoAnnotationAdapter({self.genome_dir.name}): yielded {count} gene→GO edges")


class MultiGoAnnotationAdapter:
    """
    Multi-strain adapter: aggregates GO nodes and gene→GO + GO-GO edges across all strains.

    Reads the same ``cyanobacteria_genomes.csv`` as :class:`MultiCyanorakNcbi`.
    Loads the GO data cache once and passes it to per-strain :class:`GoAnnotationAdapter`
    instances.

    Args:
        genome_config_file: path to ``cyanobacteria_genomes.csv``
        cache_root: project-level cache directory (e.g. ``Path("cache/data")``)
        test_mode: passed to each per-strain adapter (stop after 100 genes)
        cache: if False, force-rebuild the GO namespace cache (re-downloads OBO)
    """

    def __init__(
        self,
        genome_config_file: str,
        cache_root: Path,
        test_mode: bool = False,
        cache: bool = True,
    ) -> None:
        self.test_mode = test_mode
        self.go_data = load_go_data(Path(cache_root), force=not cache)
        self.adapters: list[GoAnnotationAdapter] = []
        self._build_adapters(genome_config_file)

    def _build_adapters(self, genome_config_file: str) -> None:
        with open(genome_config_file, newline="", encoding="utf-8") as fh:
            lines = [line for line in fh if not line.lstrip().startswith("#")]
        reader = csv.DictReader(lines)
        for row in reader:
            data_dir = row.get("data_dir", "").strip()
            if not data_dir:
                continue
            adapter = GoAnnotationAdapter(
                genome_dir=Path(data_dir),
                go_data=self.go_data,
                test_mode=self.test_mode,
            )
            self.adapters.append(adapter)
        logger.info(f"MultiGoAnnotationAdapter: loaded {len(self.adapters)} strain adapters")

    def _all_seed_go_ids(self) -> set[str]:
        """Collect all directly-referenced GO IDs across all strains."""
        seed: set[str] = set()
        for adapter in self.adapters:
            seed |= adapter.get_all_go_ids()
        return seed

    def get_nodes(self):
        """
        Yield GO nodes for all terms in the ancestry closure.

        1. Collect seed GO IDs from gene annotations across all strains.
        2. Compute closure (seed + all transitive ancestors in go_data).
        3. Yield ``(node_id, label, {"name": name})`` for each term in closure
           whose namespace is recognized.
        """
        seed = self._all_seed_go_ids()
        logger.info(f"Seed GO terms from gene annotations: {len(seed)}")
        closure = compute_ancestry_closure(seed, self.go_data)
        logger.info(f"GO ancestry closure size: {len(closure)} terms")

        count = 0
        for go_id in closure:
            entry = self.go_data.get(go_id)
            if entry is None:
                continue
            label = NAMESPACE_TO_LABEL.get(entry["namespace"])
            if label is None:
                continue
            yield (
                _go_node_id(go_id),
                label,
                {"name": entry["name"]},
            )
            count += 1

        logger.info(f"MultiGoAnnotationAdapter.get_nodes: yielded {count} GO nodes")

    def get_edges(self):
        """
        Yield both gene→GO edges and GO-GO hierarchy edges.

        GO-GO edges are emitted for each parent relationship where both child and
        parent are in the ancestry closure and the resulting edge label is in the
        allowed schema set.
        """
        # --- gene → GO edges ---
        gene_go_count = 0
        for adapter in self.adapters:
            for edge in adapter.get_edges():
                yield edge
                gene_go_count += 1

        logger.info(f"MultiGoAnnotationAdapter.get_edges: yielded {gene_go_count} gene→GO edges")

        # --- GO-GO hierarchy edges ---
        seed = self._all_seed_go_ids()
        closure = compute_ancestry_closure(seed, self.go_data)

        go_go_count = 0
        seen_edges: set[tuple] = set()
        for go_id in closure:
            entry = self.go_data.get(go_id)
            if entry is None:
                continue
            child_ns = entry.get("namespace", "")

            for parent_id, relation in entry.get("parents", []):
                if parent_id not in closure:
                    continue
                parent_entry = self.go_data.get(parent_id)
                if parent_entry is None:
                    continue
                parent_ns = parent_entry.get("namespace", "")

                edge_label = make_go_go_edge_label(child_ns, relation, parent_ns)
                if edge_label is None:
                    continue

                key = (go_id, parent_id, edge_label)
                if key in seen_edges:
                    continue
                seen_edges.add(key)

                yield (
                    None,
                    _go_node_id(go_id),
                    _go_node_id(parent_id),
                    edge_label,
                    {},
                )
                go_go_count += 1

        logger.info(f"MultiGoAnnotationAdapter.get_edges: yielded {go_go_count} GO-GO edges")
