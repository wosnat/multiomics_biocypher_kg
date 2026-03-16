"""OrthologGroup adapter.

Reads pre-computed ortholog_groups from gene_annotations_merged.json
(written by build_gene_annotations.py Phase 1) and yields:
  - OrthologGroup nodes (deduplicated across strains)
  - Gene_in_ortholog_group edges
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from biocypher._logger import logger


class OrthologGroupAdapter:
    """Per-strain: reads pre-computed ortholog_groups from gene_annotations_merged.json."""

    def __init__(self, genome_dir: Path, test_mode: bool = False):
        self.genome_dir = Path(genome_dir)
        self.test_mode = test_mode
        self._genes: dict = {}
        self._load()

    def _load(self):
        path = self.genome_dir / "gene_annotations_merged.json"
        if path.exists():
            with open(path) as fh:
                self._genes = json.load(fh)
        else:
            logger.warning(f"gene_annotations_merged.json not found at {path}")

    def get_og_memberships(self) -> list[tuple[str, dict]]:
        """Return (locus_tag, og_dict) pairs for all genes.

        Reads pre-computed ortholog_groups field from JSON
        (written by build_gene_annotations.py in Phase 1).
        """
        results = []
        for lt, gene in self._genes.items():
            for og in gene.get("ortholog_groups", []):
                results.append((lt, og))
            if self.test_mode and len(results) >= 100:
                break
        return results


class MultiOrthologGroupAdapter:
    """Multi-strain: yields OrthologGroup nodes + Gene_in_ortholog_group edges."""

    def __init__(self, genome_config_file: str, test_mode: bool = False):
        self.adapters: list[OrthologGroupAdapter] = []
        self._build_adapters(genome_config_file, test_mode)
        self.test_mode = test_mode

    def _build_adapters(self, genome_config_file: str, test_mode: bool):
        with open(genome_config_file, newline="", encoding="utf-8") as fh:
            lines = [line for line in fh if not line.lstrip().startswith("#")]
        reader = csv.DictReader(lines)
        for row in reader:
            data_dir = row.get("data_dir", "").strip()
            if not data_dir:
                continue
            self.adapters.append(
                OrthologGroupAdapter(genome_dir=Path(data_dir), test_mode=test_mode)
            )
        logger.info(f"OrthologGroupAdapter: loaded {len(self.adapters)} strains from {genome_config_file}")

    def download_data(self, **kwargs):
        """No-op: data already loaded in __init__."""
        pass

    def get_nodes(self):
        """Yield unique OrthologGroup nodes across all strains."""
        seen = set()
        node_list = []
        for adapter in self.adapters:
            for lt, og in adapter.get_og_memberships():
                og_id = og["og_id"]
                if og_id not in seen:
                    seen.add(og_id)
                    node_list.append((
                        og_id,
                        "ortholog_group",
                        {
                            "source": og["source"],
                            "taxonomic_level": og["taxonomic_level"],
                            "taxon_id": og["taxon_id"],
                        },
                    ))
        logger.info(f"OrthologGroupAdapter: {len(node_list)} unique OrthologGroup nodes")
        return node_list

    def get_edges(self):
        """Yield Gene_in_ortholog_group edges."""
        edge_list = []
        for adapter in self.adapters:
            for lt, og in adapter.get_og_memberships():
                gene_id = f"ncbigene:{lt}"
                edge_list.append((
                    f"{lt}-og-{og['og_id']}",
                    gene_id,
                    og["og_id"],
                    "gene_in_ortholog_group",
                    {},
                ))
        logger.info(f"OrthologGroupAdapter: {len(edge_list)} Gene_in_ortholog_group edges")
        return edge_list
