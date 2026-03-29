"""OrthologGroup adapter.

Reads pre-computed ortholog_groups from gene_annotations_merged.json
(written by build_gene_annotations.py Phase 1) and yields:
  - OrthologGroup nodes (deduplicated across strains, with eggNOG descriptions
    and functional_description derived from member gene annotations)
  - Gene_in_ortholog_group edges
  - OG_has_cyanorak_role edges (majority-vote from member genes)
  - OG_in_cog_category edges (majority-vote from member genes)
"""

from __future__ import annotations

import csv
import json
import os
import sqlite3
from collections import Counter
from pathlib import Path

from biocypher._logger import logger

from multiomics_kg.adapters.functional_annotation_adapter import (
    COG_FUNCTIONAL_CATEGORIES,
    _cog_cat_node_id,
    _cyanorak_role_node_id,
)
from multiomics_kg.utils.cyanorak_role_utils import (
    full_role_description,
    parse_cyanorak_role_tree,
)

# Cyanorak role codes to exclude from functional_description (uninformative)
_UNINFORMATIVE_CYANORAK_CODES = frozenset({
    "R",      # "Other" (top-level)
    "R.2",    # "Conserved hypothetical proteins"
    "R.4",    # "Hypothetical proteins"
})

# COG categories to exclude from functional_description
_UNINFORMATIVE_COG_CATEGORIES = frozenset({
    "S",  # "Function unknown"
})


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


def _majority_codes(members: list[dict], field: str) -> list[str]:
    """Return codes that appear in >50% of members (majority vote).

    Each member's field is a list of codes. A gene with multiple codes
    contributes one vote per code. The threshold is >50% of member count.
    """
    n = len(members)
    if n == 0:
        return []
    counter: Counter = Counter()
    for m in members:
        for code in m.get(field, []):
            counter[code] += 1
    threshold = n / 2
    return sorted(code for code, count in counter.items() if count > threshold)


def _load_eggnog_descriptions(og_names: set[tuple[str, str]]) -> dict[tuple[str, str], str]:
    """Load descriptions for the given (og_name, level) pairs.

    Tries cache/data/eggnog/og_descriptions.json first (lightweight, Docker-safe),
    then falls back to querying eggnog.db directly.

    Returns a dict mapping (og_name, level) -> description string.
    """
    # 1. Try pre-computed cache file (written by build_og_descriptions.py / prepare_data step 5)
    project_root = Path(__file__).parent.parent.parent
    cache_path = project_root / "cache" / "data" / "eggnog" / "og_descriptions.json"
    if cache_path.exists():
        try:
            with open(cache_path) as f:
                raw = json.load(f)  # {"COG0592@2": "description", ...}
            descriptions: dict[tuple[str, str], str] = {}
            for key, desc in raw.items():
                if "@" in key:
                    name, level = key.rsplit("@", 1)
                    descriptions[(name, level)] = desc
            logger.info(f"Loaded {len(descriptions)} eggNOG OG descriptions from cache ({cache_path})")
            return descriptions
        except Exception as e:
            logger.warning(f"Error reading OG descriptions cache: {e}")

    # 2. Fallback: query eggnog.db directly
    eggnog_dir = os.environ.get("EGGNOG_DATA_DIR", "")
    if not eggnog_dir:
        env_path = project_root / ".env"
        if env_path.exists():
            for line in env_path.read_text().splitlines():
                line = line.strip()
                if line.startswith("EGGNOG_DATA_DIR="):
                    eggnog_dir = line.split("=", 1)[1].strip().strip('"').strip("'")
                    break

    if not eggnog_dir:
        logger.warning(
            "No OG descriptions cache and EGGNOG_DATA_DIR not set; "
            "run 'prepare_data.sh --steps 5' to build the cache"
        )
        return {}

    db_path = Path(eggnog_dir).expanduser() / "eggnog.db"
    if not db_path.exists():
        logger.warning(f"eggnog.db not found at {db_path}")
        return {}

    descriptions = {}
    try:
        conn = sqlite3.connect(str(db_path))
        unique_og_names = {name for name, _ in og_names}
        if not unique_og_names:
            conn.close()
            return {}

        batch_size = 500
        og_list = list(unique_og_names)
        for i in range(0, len(og_list), batch_size):
            batch = og_list[i:i + batch_size]
            placeholders = ",".join("?" * len(batch))
            cursor = conn.execute(
                f"SELECT og, level, description FROM og WHERE og IN ({placeholders})",
                batch,
            )
            for og, level, desc in cursor:
                if desc:
                    descriptions[(og, str(level))] = desc

        conn.close()
        logger.info(f"Loaded {len(descriptions)} eggNOG OG descriptions from {db_path}")
    except Exception as e:
        logger.warning(f"Error reading eggnog.db: {e}")

    return descriptions


def _parse_eggnog_og_id(og_id: str) -> tuple[str, str] | None:
    """Parse an og_id like 'eggnog:COG0592@2' into ('COG0592', '2').

    Returns None for non-eggnog IDs or unparseable formats.
    """
    if not og_id.startswith("eggnog:"):
        return None
    raw = og_id[len("eggnog:"):]  # e.g. "COG0592@2"
    if "@" not in raw:
        return None
    name, level = raw.rsplit("@", 1)
    return (name, level)


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

        gene_meta contains product, gene_name, organism_name,
        cog_category, and cyanorak_Role for consensus computation.
        """
        results = []
        for lt, gene in self._genes.items():
            meta = {
                "product": gene.get("product"),
                "gene_name": gene.get("gene_name"),
                "organism_name": gene.get("organism_name"),
                "cog_category": gene.get("cog_category", []),
                "cyanorak_Role": gene.get("cyanorak_Role", []),
            }
            for og in gene.get("ortholog_groups", []):
                results.append((lt, og, meta))
            if self.test_mode and len(results) >= 100:
                break
        return results


class MultiOrthologGroupAdapter:
    """Multi-strain: yields OrthologGroup nodes + Gene_in_ortholog_group edges
    + OG_has_cyanorak_role + OG_in_cog_category edges."""

    def __init__(
        self,
        genome_config_file: str,
        test_mode: bool = False,
        role_tree_file: str = "data/cyanorak_roles.csv",
    ):
        self.adapters: list[OrthologGroupAdapter] = []
        self._build_adapters(genome_config_file, test_mode)
        self.test_mode = test_mode

        # Load Cyanorak role tree for functional_description
        role_tree_path = Path(role_tree_file)
        if role_tree_path.exists():
            self.role_tree = parse_cyanorak_role_tree(role_tree_path)
        else:
            logger.warning(f"cyanorak_roles.csv not found at {role_tree_path}")
            self.role_tree = {}

        # Cache for OG info (populated by get_nodes, reused by get_edges)
        self._og_info: dict = {}

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

    def _build_functional_description(self, majority_cog: list[str], majority_cyanorak: list[str]) -> str | None:
        """Build a semicolon-separated functional description from majority roles/categories.

        Uses full hierarchical names for Cyanorak roles and COG category names.
        Filters out uninformative entries.
        """
        parts: list[str] = []

        # Cyanorak roles (hierarchical names)
        for code in majority_cyanorak:
            if code in _UNINFORMATIVE_CYANORAK_CODES:
                continue
            desc = full_role_description(code, self.role_tree)
            if desc:
                parts.append(desc)

        # COG categories (full names)
        for letter in majority_cog:
            if letter in _UNINFORMATIVE_COG_CATEGORIES:
                continue
            name = COG_FUNCTIONAL_CATEGORIES.get(letter)
            if name:
                parts.append(name)

        if not parts:
            return None
        return _clean_str("; ".join(parts))

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

        # Load eggNOG descriptions from local DB
        eggnog_keys: set[tuple[str, str]] = set()
        for og_id in og_info:
            parsed = _parse_eggnog_og_id(og_id)
            if parsed:
                eggnog_keys.add(parsed)
        eggnog_descs = _load_eggnog_descriptions(eggnog_keys)

        # Second pass: compute consensus and emit nodes
        node_list = []
        desc_count = 0
        func_desc_count = 0
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
            org_strains = {m["organism_name"] for m in members if m.get("organism_name")}
            genera = sorted({s.split()[0] for s in org_strains if s})

            # Description from eggNOG DB (null for Cyanorak groups)
            description = None
            parsed = _parse_eggnog_og_id(og_id)
            if parsed:
                description = eggnog_descs.get(parsed)
            if description:
                desc_count += 1

            # Majority-vote for functional annotations
            majority_cog = _majority_codes(members, "cog_category")
            majority_cyanorak = _majority_codes(members, "cyanorak_Role")

            # Store majority codes for edge generation
            info["majority_cog"] = majority_cog
            info["majority_cyanorak"] = majority_cyanorak

            # Build functional_description
            functional_description = self._build_functional_description(
                majority_cog, majority_cyanorak
            )
            if functional_description:
                func_desc_count += 1

            props = {
                "name": raw_name,
                "source": og["source"],
                "taxonomic_level": og["taxonomic_level"],
                "taxon_id": og["taxon_id"],
                "specificity_rank": og["specificity_rank"],
                "consensus_product": _clean_str(consensus_product) if consensus_product else None,
                "consensus_gene_name": _clean_str(consensus_gene_name) if consensus_gene_name else None,
                "description": _clean_str(description) if description else None,
                "functional_description": functional_description,
                "member_count": len(members),
                "organism_count": len(org_strains),
                "genera": genera,
                "has_cross_genus_members": "cross_genus" if len(genera) > 1 else "single_genus",
            }
            node_list.append((og_id, "ortholog_group", props))

        # Cache for get_edges
        self._og_info = og_info

        logger.info(
            f"OrthologGroupAdapter: {len(node_list)} unique OrthologGroup nodes, "
            f"{desc_count} with description, {func_desc_count} with functional_description"
        )
        return node_list

    def get_edges(self):
        """Yield Gene_in_ortholog_group + OG_has_cyanorak_role + OG_in_cog_category edges."""
        edge_list = []

        # 1. Gene_in_ortholog_group edges
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
        gene_og_count = len(edge_list)

        # 2. OG_has_cyanorak_role edges (from majority vote, computed in get_nodes)
        cyanorak_edge_count = 0
        for og_id, info in self._og_info.items():
            for code in info.get("majority_cyanorak", []):
                edge_list.append((
                    f"{og_id}-cyanorak_role-{code}",
                    og_id,
                    _cyanorak_role_node_id(code),
                    "og_has_cyanorak_role",
                    {},
                ))
                cyanorak_edge_count += 1

        # 3. OG_in_cog_category edges (from majority vote, computed in get_nodes)
        cog_edge_count = 0
        for og_id, info in self._og_info.items():
            for letter in info.get("majority_cog", []):
                edge_list.append((
                    f"{og_id}-cog_category-{letter}",
                    og_id,
                    _cog_cat_node_id(letter),
                    "og_in_cog_category",
                    {},
                ))
                cog_edge_count += 1

        logger.info(
            f"OrthologGroupAdapter: {gene_og_count} Gene_in_ortholog_group, "
            f"{cyanorak_edge_count} OG_has_cyanorak_role, "
            f"{cog_edge_count} OG_in_cog_category edges"
        )
        return edge_list
