"""OrthologGroup adapter.

Reads pre-computed ortholog_groups from gene_annotations_merged.json
(written by build_gene_annotations.py Phase 1) and yields:
  - OrthologGroup nodes (deduplicated across strains)
  - Gene_in_ortholog_group edges
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from biocypher._logger import logger


def _clean_str(value: str) -> str:
    """Sanitize strings for BioCypher CSV export."""
    return value.replace("'", "^").replace("|", "")


def _consensus_value(values: list, exclude: set | None = None) -> str | None:
    """Return the most common non-null value, preferring values not in *exclude*.

    Falls back to the most common excluded value if nothing else is available.
    """
    non_null = [v for v in values if v]
    if not non_null:
        return None
    if exclude:
        preferred = [v for v in non_null if v not in exclude]
        if preferred:
            return Counter(preferred).most_common(1)[0][0]
    # Fall back to most common overall (including excluded)
    return Counter(non_null).most_common(1)[0][0]


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

    def get_og_memberships_with_gene_data(self) -> list[tuple[str, dict, dict]]:
        """Return (locus_tag, og_dict, gene_meta) triples.

        gene_meta contains product, gene_name, organism_strain for consensus computation.
        """
        results = []
        for lt, gene in self._genes.items():
            meta = {
                "product": gene.get("product"),
                "gene_name": gene.get("gene_name"),
                "organism_strain": gene.get("organism_strain"),
            }
            for og in gene.get("ortholog_groups", []):
                results.append((lt, og, meta))
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
        """Yield unique OrthologGroup nodes with consensus properties."""
        # First pass: collect members per OG
        og_info = {}  # og_id -> {"og": og_dict, "members": [gene_meta, ...]}
        for adapter in self.adapters:
            for lt, og, meta in adapter.get_og_memberships_with_gene_data():
                og_id = og["og_id"]
                if og_id not in og_info:
                    og_info[og_id] = {"og": og, "members": []}
                og_info[og_id]["members"].append(meta)

        # Second pass: compute consensus and emit nodes
        node_list = []
        for og_id, info in og_info.items():
            og = info["og"]
            members = info["members"]
            raw_name = og_id.split(":", 1)[1] if ":" in og_id else og_id

            # Consensus product: majority vote, preferring non-hypothetical
            consensus_product = _consensus_value(
                [m["product"] for m in members],
                exclude={"hypothetical protein", "conserved hypothetical protein"},
            )

            # Consensus gene name: most frequent non-null
            consensus_gene_name = _consensus_value(
                [m["gene_name"] for m in members],
            )

            # Organism stats
            org_strains = {m["organism_strain"] for m in members if m.get("organism_strain")}
            genera = sorted({s.split()[0] for s in org_strains if s})

            props = {
                "name": raw_name,
                "source": og["source"],
                "taxonomic_level": og["taxonomic_level"],
                "taxon_id": og["taxon_id"],
                "specificity_rank": og["specificity_rank"],
                "consensus_product": _clean_str(consensus_product) if consensus_product else None,
                "consensus_gene_name": _clean_str(consensus_gene_name) if consensus_gene_name else None,
                "member_count": len(members),
                "organism_count": len(org_strains),
                "genera": genera,
                "has_cross_genus_members": "cross_genus" if len(genera) > 1 else "single_genus",
            }
            node_list.append((og_id, "ortholog_group", props))

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
