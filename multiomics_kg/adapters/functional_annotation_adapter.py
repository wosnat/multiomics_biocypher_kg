"""
Functional annotation adapters: gene→GO, gene→EC, gene→KEGG, and gene→COG/role edges,
plus corresponding nodes.

GO section (MultiGoAnnotationAdapter):
- GO nodes (BiologicalProcess, CellularComponent, MolecularFunction) for only the
  ~5K GO terms referenced in gene annotations plus their ancestors (not all 30K).
- gene_involved_in_biological_process / gene_located_in_cellular_component /
  gene_enables_molecular_function edges
- GO-GO hierarchy edges (is_a, part_of, etc.) for the ancestor closure only.
The existing --go flag path is preserved for backward compatibility but should
not be run simultaneously (GO node IDs would conflict).

EC section (MultiEcAnnotationAdapter):
- EC number nodes (full Expasy hierarchy, with rich metadata) cached at
  cache/data/ec/ec_data.json.
- ec_number_is_a_ec_number hierarchy edges.
- gene_catalyzes_ec_number edges from gene annotations.
Replaces the standalone EC adapter in the main pipeline; EC adapter is kept for
backward compatibility.

KEGG section (MultiKeggAnnotationAdapter):
- 4-level KEGG hierarchy: KO → Pathway → Subcategory → Category.
  Data from KEGG REST API, cached at cache/data/kegg/kegg_data.json.
- gene_has_kegg_ko edges from gene annotations (kegg_ko field).
- ko_in_kegg_pathway, kegg_pathway_in_kegg_subcategory,
  kegg_subcategory_in_kegg_category hierarchy edges.

COG/Role section (MultiCogRoleAnnotationAdapter):
- CogFunctionalCategory nodes (25 standard letters, hardcoded).
- CyanorakRole nodes (full ~172-node tree from data/cyanorak_roles.csv).
- TigrRole nodes (only codes present in data, Pro/Syn strains only).
- gene_in_cog_category edges (all 13 strains, from cog_category field).
- gene_has_cyanorak_role edges (Pro/Syn 6 strains only).
- cyanorak_role_is_a_cyanorak_role hierarchy edges (full tree).
- gene_has_tigr_role edges (Pro/Syn 6 strains only).
"""

import csv
import json
import logging
from pathlib import Path

from bioregistry import normalize_curie

from multiomics_kg.adapters.ec_adapter import EC
from multiomics_kg.utils.kegg_utils import download_kegg_data
from multiomics_kg.utils.cyanorak_role_utils import parse_cyanorak_role_tree, full_role_description
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


def _clean_str(value: str) -> str:
    """Sanitize a string for Neo4j CSV import: replace single quotes and pipes."""
    return value.replace("'", "^").replace("|", "")


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
                {"name": _clean_str(entry["name"])},
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


# ---------------------------------------------------------------------------
# EC number annotation: gene → EC edges
# ---------------------------------------------------------------------------

def _ec_node_id(ec_number: str) -> str:
    """Normalize EC number to the eccode node ID format (matches ec_adapter.py)."""
    return normalize_curie(f"eccode:{ec_number}") or f"eccode_{ec_number}"


class EcAnnotationAdapter:
    """
    Per-strain adapter: reads gene_annotations_merged.json and yields gene→EC edges.

    Args:
        genome_dir: directory containing ``gene_annotations_merged.json``
        test_mode: if True, stop after 100 edges
    """

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
        logger.debug(f"Loaded {len(self._genes)} genes from {json_path}")

    def get_edges(self):
        """
        Yield gene→EC edges as 5-tuples::

            (edge_id, gene_node_id, ec_node_id, edge_label, properties)
        """
        count = 0
        for locus_tag, gene in self._genes.items():
            for ec_num in gene.get("ec_numbers") or []:
                if not ec_num:
                    continue
                yield (
                    f"{locus_tag}-ec-{ec_num}",
                    _gene_node_id(locus_tag),
                    _ec_node_id(ec_num),
                    "gene_catalyzes_ec_number",
                    {},
                )
                count += 1
                if self.test_mode and count >= 100:
                    logger.debug(
                        f"EcAnnotationAdapter({self.genome_dir.name}): test_mode stop at {count}"
                    )
                    return
        logger.debug(f"EcAnnotationAdapter({self.genome_dir.name}): yielded {count} gene→EC edges")


class MultiEcAnnotationAdapter:
    """
    Multi-strain adapter: owns all EC graph content (mirrors MultiGoAnnotationAdapter).

    Creates:
    - EC number nodes (full Expasy hierarchy with rich metadata), cached at
      ``cache_dir/ec_data.json``.
    - ``ec_number_is_a_ec_number`` hierarchy edges.
    - ``gene_catalyzes_ec_number`` edges from per-strain gene annotations.

    Args:
        genome_config_file: path to ``cyanobacteria_genomes.csv``
        cache_dir: directory for EC JSON cache (e.g. ``Path("cache/data/ec")``)
        test_mode: passed to EC adapter and per-strain adapters (stop after 100 items)
        cache: if False, force re-download even if cache exists
    """

    def __init__(
        self,
        genome_config_file: str,
        cache_dir: Path,
        test_mode: bool = False,
        cache: bool = True,
    ) -> None:
        self.test_mode = test_mode
        self._ec = EC(test_mode=test_mode, cache_dir=cache_dir)
        self._ec.download_ec_data(cache=cache)
        self._strain_adapters: list[EcAnnotationAdapter] = []
        self._build_strain_adapters(genome_config_file)

    def _build_strain_adapters(self, genome_config_file: str) -> None:
        with open(genome_config_file, newline="", encoding="utf-8") as fh:
            lines = [line for line in fh if not line.lstrip().startswith("#")]
        reader = csv.DictReader(lines)
        for row in reader:
            data_dir = row.get("data_dir", "").strip()
            if not data_dir:
                continue
            self._strain_adapters.append(
                EcAnnotationAdapter(genome_dir=Path(data_dir), test_mode=self.test_mode)
            )
        logger.info(f"MultiEcAnnotationAdapter: loaded {len(self._strain_adapters)} strain adapters")

    def get_nodes(self):
        """Yield EC number nodes (full Expasy hierarchy)."""
        yield from self._ec.get_nodes()

    def get_edges(self):
        """Yield EC hierarchy edges then gene→EC edges for all strains."""
        yield from self._ec.get_edges()
        count = 0
        for adapter in self._strain_adapters:
            for edge in adapter.get_edges():
                yield edge
                count += 1
        logger.info(f"MultiEcAnnotationAdapter.get_edges: yielded {count} gene→EC edges")


# ---------------------------------------------------------------------------
# KEGG annotation: gene→KO edges + 4-level KEGG hierarchy nodes/edges
# ---------------------------------------------------------------------------

def _ko_node_id(ko_id: str) -> str:
    """Normalize a KEGG Orthology ID to the node ID format."""
    return normalize_curie(f"kegg.orthology:{ko_id}") or f"kegg.orthology:{ko_id}"


def _pathway_node_id(pw_id: str) -> str:
    """Normalize a KEGG reference pathway ID to the node ID format."""
    return normalize_curie(f"kegg.pathway:{pw_id}") or f"kegg.pathway:{pw_id}"


def _subcat_node_id(code: str) -> str:
    """Node ID for a KEGG BRITE B-level subcategory."""
    return f"kegg.subcategory:{code}"


def _cat_node_id(code: str) -> str:
    """Node ID for a KEGG BRITE A-level category."""
    return f"kegg.category:{code}"


class KeggAnnotationAdapter:
    """
    Per-strain adapter: reads gene_annotations_merged.json and yields gene→KO edges.

    Args:
        genome_dir: directory containing ``gene_annotations_merged.json``
        test_mode: if True, stop after 100 edges
    """

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
        logger.debug(f"Loaded {len(self._genes)} genes from {json_path}")

    def get_all_ko_ids(self) -> set[str]:
        """Return all KO IDs (K#####) directly referenced by any gene in this strain."""
        ids: set[str] = set()
        for gene in self._genes.values():
            for ko_id in gene.get("kegg_ko") or []:
                if ko_id:
                    ids.add(ko_id)
        return ids

    def get_edges(self):
        """
        Yield gene→KO edges as 5-tuples::

            (edge_id, gene_node_id, ko_node_id, edge_label, properties)
        """
        count = 0
        for locus_tag, gene in self._genes.items():
            for ko_id in gene.get("kegg_ko") or []:
                if not ko_id:
                    continue
                yield (
                    f"{locus_tag}-kegg-{ko_id}",
                    _gene_node_id(locus_tag),
                    _ko_node_id(ko_id),
                    "gene_has_kegg_ko",
                    {},
                )
                count += 1
                if self.test_mode and count >= 100:
                    logger.debug(
                        f"KeggAnnotationAdapter({self.genome_dir.name}): test_mode stop at {count}"
                    )
                    return
        logger.debug(
            f"KeggAnnotationAdapter({self.genome_dir.name}): yielded {count} gene→KO edges"
        )


class MultiKeggAnnotationAdapter:
    """
    Multi-strain adapter: owns all KEGG graph content.

    Creates:
    - ``KeggOrthologousGroup`` nodes (K#####) from gene annotations + KEGG REST API names.
    - ``KeggPathway`` nodes (ko#####) from KEGG REST API.
    - ``KeggSubcategory`` nodes (B-level BRITE codes) from KEGG REST API.
    - ``KeggCategory`` nodes (A-level BRITE codes) from KEGG REST API.
    - ``gene_has_kegg_ko`` edges from per-strain gene annotations.
    - ``ko_in_kegg_pathway`` edges from KEGG REST API KO→Pathway links.
    - ``kegg_pathway_in_kegg_subcategory`` edges from BRITE hierarchy.
    - ``kegg_subcategory_in_kegg_category`` edges from BRITE hierarchy.

    Args:
        genome_config_file: path to ``cyanobacteria_genomes.csv``
        cache_root: project-level cache directory (e.g. ``Path("cache/data")``)
        test_mode: passed to per-strain adapters (stop after 100 gene→KO edges)
        cache: if False, force re-download KEGG data even if cache exists
    """

    def __init__(
        self,
        genome_config_file: str,
        cache_root: Path,
        test_mode: bool = False,
        cache: bool = True,
    ) -> None:
        self.test_mode = test_mode
        self.kegg_data = download_kegg_data(Path(cache_root), force=not cache)
        self._strain_adapters: list[KeggAnnotationAdapter] = []
        self._build_strain_adapters(genome_config_file)

    def _build_strain_adapters(self, genome_config_file: str) -> None:
        with open(genome_config_file, newline="", encoding="utf-8") as fh:
            lines = [line for line in fh if not line.lstrip().startswith("#")]
        reader = csv.DictReader(lines)
        for row in reader:
            data_dir = row.get("data_dir", "").strip()
            if not data_dir:
                continue
            self._strain_adapters.append(
                KeggAnnotationAdapter(genome_dir=Path(data_dir), test_mode=self.test_mode)
            )
        logger.info(f"MultiKeggAnnotationAdapter: loaded {len(self._strain_adapters)} strain adapters")

    def _all_ko_ids(self) -> set[str]:
        """Collect all KO IDs referenced across all strains."""
        ids: set[str] = set()
        for adapter in self._strain_adapters:
            ids |= adapter.get_all_ko_ids()
        return ids

    def get_nodes(self):
        """
        Yield KEGG hierarchy nodes (4 types), deduplicated.

        Order: KO → Pathway → Subcategory → Category (specific → general).
        Only nodes reachable from KOs referenced in the gene annotations are emitted.
        """
        ko_ids = self._all_ko_ids()
        logger.info(f"MultiKeggAnnotationAdapter: {len(ko_ids)} unique KO IDs across all strains")

        ko_names = self.kegg_data.get("ko_names", {})
        ko_to_pathways = self.kegg_data.get("ko_to_pathways", {})
        pathway_names = self.kegg_data.get("pathway_names", {})
        pathway_to_subcat = self.kegg_data.get("pathway_to_subcategory", {})
        subcat_names = self.kegg_data.get("subcategory_names", {})
        subcat_to_cat = self.kegg_data.get("subcategory_to_category", {})
        cat_names = self.kegg_data.get("category_names", {})

        # Collect pathway, subcategory, category IDs reachable from our KOs
        pw_ids: set[str] = set()
        for ko_id in ko_ids:
            pw_ids.update(ko_to_pathways.get(ko_id, []))

        subcat_ids: set[str] = set()
        for pw_id in pw_ids:
            sc = pathway_to_subcat.get(pw_id)
            if sc:
                subcat_ids.add(sc)

        cat_ids: set[str] = set()
        for sc in subcat_ids:
            cat = subcat_to_cat.get(sc)
            if cat:
                cat_ids.add(cat)

        # Yield KO nodes
        ko_count = 0
        for ko_id in sorted(ko_ids):
            name = _clean_str(ko_names.get(ko_id, ""))
            yield (_ko_node_id(ko_id), "kegg orthologous group", {"name": name})
            ko_count += 1
        logger.info(f"MultiKeggAnnotationAdapter.get_nodes: {ko_count} KO nodes")

        # Yield Pathway nodes
        pw_count = 0
        for pw_id in sorted(pw_ids):
            name = _clean_str(pathway_names.get(pw_id, ""))
            yield (_pathway_node_id(pw_id), "kegg pathway", {"name": name})
            pw_count += 1
        logger.info(f"MultiKeggAnnotationAdapter.get_nodes: {pw_count} pathway nodes")

        # Yield Subcategory nodes
        sc_count = 0
        for sc in sorted(subcat_ids):
            name = _clean_str(subcat_names.get(sc, ""))
            yield (_subcat_node_id(sc), "kegg subcategory", {"name": name})
            sc_count += 1
        logger.info(f"MultiKeggAnnotationAdapter.get_nodes: {sc_count} subcategory nodes")

        # Yield Category nodes
        cat_count = 0
        for cat in sorted(cat_ids):
            name = _clean_str(cat_names.get(cat, ""))
            yield (_cat_node_id(cat), "kegg category", {"name": name})
            cat_count += 1
        logger.info(f"MultiKeggAnnotationAdapter.get_nodes: {cat_count} category nodes")

    def get_edges(self):
        """
        Yield all KEGG edges in order:
        1. gene→KO (gene_has_kegg_ko) from per-strain adapters
        2. KO→Pathway (ko_in_kegg_pathway) from KEGG REST API
        3. Pathway→Subcategory (kegg_pathway_in_kegg_subcategory)
        4. Subcategory→Category (kegg_subcategory_in_kegg_category)
        """
        ko_to_pathways = self.kegg_data.get("ko_to_pathways", {})
        pathway_to_subcat = self.kegg_data.get("pathway_to_subcategory", {})
        subcat_to_cat = self.kegg_data.get("subcategory_to_category", {})

        # 1. gene → KO edges
        gene_ko_count = 0
        for adapter in self._strain_adapters:
            for edge in adapter.get_edges():
                yield edge
                gene_ko_count += 1
        logger.info(f"MultiKeggAnnotationAdapter.get_edges: {gene_ko_count} gene→KO edges")

        # Collect KOs and pathways actually emitted
        ko_ids = self._all_ko_ids()
        pw_ids: set[str] = set()
        for ko_id in ko_ids:
            pw_ids.update(ko_to_pathways.get(ko_id, []))

        # 2. KO → Pathway edges
        ko_pw_count = 0
        seen_ko_pw: set[tuple] = set()
        for ko_id in sorted(ko_ids):
            for pw_id in ko_to_pathways.get(ko_id, []):
                key = (ko_id, pw_id)
                if key in seen_ko_pw:
                    continue
                seen_ko_pw.add(key)
                yield (
                    f"{ko_id}-in-{pw_id}",
                    _ko_node_id(ko_id),
                    _pathway_node_id(pw_id),
                    "ko_in_kegg_pathway",
                    {},
                )
                ko_pw_count += 1
        logger.info(f"MultiKeggAnnotationAdapter.get_edges: {ko_pw_count} KO→Pathway edges")

        # 3. Pathway → Subcategory edges
        pw_sc_count = 0
        seen_pw_sc: set[tuple] = set()
        for pw_id in sorted(pw_ids):
            sc = pathway_to_subcat.get(pw_id)
            if not sc:
                continue
            key = (pw_id, sc)
            if key in seen_pw_sc:
                continue
            seen_pw_sc.add(key)
            yield (
                f"{pw_id}-in-{sc}",
                _pathway_node_id(pw_id),
                _subcat_node_id(sc),
                "kegg_pathway_in_kegg_subcategory",
                {},
            )
            pw_sc_count += 1
        logger.info(
            f"MultiKeggAnnotationAdapter.get_edges: {pw_sc_count} Pathway→Subcategory edges"
        )

        # 4. Subcategory → Category edges
        subcat_ids: set[str] = {pathway_to_subcat[pw] for pw in pw_ids if pw in pathway_to_subcat}
        sc_cat_count = 0
        for sc in sorted(subcat_ids):
            cat = subcat_to_cat.get(sc)
            if not cat:
                continue
            yield (
                f"{sc}-in-{cat}",
                _subcat_node_id(sc),
                _cat_node_id(cat),
                "kegg_subcategory_in_kegg_category",
                {},
            )
            sc_cat_count += 1
        logger.info(
            f"MultiKeggAnnotationAdapter.get_edges: {sc_cat_count} Subcategory→Category edges"
        )


# ---------------------------------------------------------------------------
# COG category + Cyanorak role + tIGR role annotation
# ---------------------------------------------------------------------------

# Standard COG functional category letter codes and names (NCBI/EggNOG).
COG_FUNCTIONAL_CATEGORIES: dict[str, str] = {
    "J": "Translation, ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "O": "Post-translational modification, protein turnover, chaperones",
    "X": "Mobilome: prophages, transposons",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
}


def _cog_cat_node_id(letter: str) -> str:
    """Node ID for a COG functional category letter."""
    return f"cog.category:{letter}"


def _cyanorak_role_node_id(code: str) -> str:
    """Node ID for a Cyanorak role code (unknown to bioregistry; raw string)."""
    return f"cyanorak.role:{code}"


def _tigr_role_node_id(code: str) -> str:
    """Node ID for a tIGR role code (unknown to bioregistry; raw string)."""
    return f"tigr.role:{code}"


class CogRoleAnnotationAdapter:
    """
    Per-strain adapter: reads gene_annotations_merged.json and yields
    gene→COG category, gene→CyanorakRole, and gene→TigrRole edges.

    Args:
        genome_dir: directory containing ``gene_annotations_merged.json``
        test_mode: if True, stop after 100 edges per edge type
    """

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
        logger.debug(f"Loaded {len(self._genes)} genes from {json_path}")

    def get_all_cyanorak_codes(self) -> set[tuple[str, str]]:
        """Return all (code, description) pairs from cyanorak_Role across all genes."""
        codes: set[tuple[str, str]] = set()
        for gene in self._genes.values():
            roles = gene.get("cyanorak_Role") or []
            descs = gene.get("cyanorak_Role_description") or []
            for code, desc in zip(roles, descs):
                if code:
                    codes.add((code, desc or ""))
        return codes

    def get_all_tigr_codes(self) -> set[tuple[str, str]]:
        """Return all (code, description) pairs from tIGR_Role across all genes."""
        codes: set[tuple[str, str]] = set()
        for gene in self._genes.values():
            roles = gene.get("tIGR_Role") or []
            descs = gene.get("tIGR_Role_description") or []
            for code, desc in zip(roles, descs):
                if code:
                    codes.add((code, desc or ""))
        return codes

    def get_edges(self):
        """
        Yield gene→COG category, gene→CyanorakRole, and gene→TigrRole edges.

        Edge 5-tuples: ``(edge_id, src_node_id, tgt_node_id, edge_label, properties)``
        Missing fields (e.g. Alteromonas has no cyanorak/tIGR roles) are silently skipped.
        """
        cog_count = cyr_count = tigr_count = 0
        for locus_tag, gene in self._genes.items():
            # gene → COG functional category
            for letter in gene.get("cog_category") or []:
                if not letter:
                    continue
                yield (
                    f"{locus_tag}-cogcat-{letter}",
                    _gene_node_id(locus_tag),
                    _cog_cat_node_id(letter),
                    "gene_in_cog_category",
                    {},
                )
                cog_count += 1
                if self.test_mode and cog_count >= 100:
                    logger.debug(
                        f"CogRoleAnnotationAdapter({self.genome_dir.name}): test_mode stop (COG)"
                    )
                    return

            # gene → Cyanorak role
            for code in gene.get("cyanorak_Role") or []:
                if not code:
                    continue
                yield (
                    f"{locus_tag}-cyrole-{code}",
                    _gene_node_id(locus_tag),
                    _cyanorak_role_node_id(code),
                    "gene_has_cyanorak_role",
                    {},
                )
                cyr_count += 1
                if self.test_mode and cyr_count >= 100:
                    logger.debug(
                        f"CogRoleAnnotationAdapter({self.genome_dir.name}): test_mode stop (CyanorakRole)"
                    )
                    return

            # gene → tIGR role
            for code in gene.get("tIGR_Role") or []:
                if not code:
                    continue
                yield (
                    f"{locus_tag}-tigrrole-{code}",
                    _gene_node_id(locus_tag),
                    _tigr_role_node_id(code),
                    "gene_has_tigr_role",
                    {},
                )
                tigr_count += 1
                if self.test_mode and tigr_count >= 100:
                    logger.debug(
                        f"CogRoleAnnotationAdapter({self.genome_dir.name}): test_mode stop (TigrRole)"
                    )
                    return

        logger.debug(
            f"CogRoleAnnotationAdapter({self.genome_dir.name}): "
            f"{cog_count} COG cat + {cyr_count} CyanorakRole + {tigr_count} TigrRole edges"
        )


class MultiCogRoleAnnotationAdapter:
    """
    Multi-strain adapter: owns all COG functional category, CyanorakRole, and
    TigrRole graph content.

    Creates:
    - ``CogFunctionalCategory`` nodes (25 standard letters, hardcoded).
    - ``CyanorakRole`` nodes (full tree from ``role_tree_file``; all ~172 nodes).
    - ``TigrRole`` nodes (only codes present in at least one strain's gene annotations).
    - ``gene_in_cog_category`` edges from per-strain adapters (all strains).
    - ``gene_has_cyanorak_role`` edges from per-strain adapters (Pro/Syn only).
    - ``cyanorak_role_is_a_cyanorak_role`` hierarchy edges (full tree).
    - ``gene_has_tigr_role`` edges from per-strain adapters (Pro/Syn only).

    Args:
        genome_config_file: path to ``cyanobacteria_genomes.csv``
        role_tree_file: path to ``data/cyanorak_roles.csv``
        test_mode: passed to per-strain adapters (stop after 100 edges per type)
    """

    def __init__(
        self,
        genome_config_file: str,
        role_tree_file: Path,
        test_mode: bool = False,
    ) -> None:
        self.test_mode = test_mode
        self.role_tree = parse_cyanorak_role_tree(Path(role_tree_file))
        logger.info(f"MultiCogRoleAnnotationAdapter: loaded {len(self.role_tree)} Cyanorak role nodes")
        self._strain_adapters: list[CogRoleAnnotationAdapter] = []
        self._build_strain_adapters(genome_config_file)

    def _build_strain_adapters(self, genome_config_file: str) -> None:
        with open(genome_config_file, newline="", encoding="utf-8") as fh:
            lines = [line for line in fh if not line.lstrip().startswith("#")]
        reader = csv.DictReader(lines)
        for row in reader:
            data_dir = row.get("data_dir", "").strip()
            if not data_dir:
                continue
            self._strain_adapters.append(
                CogRoleAnnotationAdapter(genome_dir=Path(data_dir), test_mode=self.test_mode)
            )
        logger.info(
            f"MultiCogRoleAnnotationAdapter: loaded {len(self._strain_adapters)} strain adapters"
        )

    def _all_tigr_codes(self) -> dict[str, str]:
        """Collect all unique tIGR (code → description) pairs across all strains."""
        codes: dict[str, str] = {}
        for adapter in self._strain_adapters:
            for code, desc in adapter.get_all_tigr_codes():
                if code not in codes:
                    codes[code] = desc
        return codes

    def get_nodes(self):
        """
        Yield COG functional category, CyanorakRole, and TigrRole nodes.

        1. 25 COG functional category nodes (always, hardcoded).
        2. All CyanorakRole nodes from the full parsed tree (always all ~172).
        3. TigrRole nodes: only codes present in at least one strain's gene annotations.
        """
        # 1. COG functional category nodes
        cog_count = 0
        for letter, name in sorted(COG_FUNCTIONAL_CATEGORIES.items()):
            yield (
                _cog_cat_node_id(letter),
                "cog functional category",
                {"code": letter, "name": _clean_str(name)},
            )
            cog_count += 1
        logger.info(f"MultiCogRoleAnnotationAdapter.get_nodes: {cog_count} COG category nodes")

        # 2. CyanorakRole nodes — full tree (all entries, not just leaf codes in data)
        cyr_count = 0
        for code, entry in sorted(self.role_tree.items()):
            yield (
                _cyanorak_role_node_id(code),
                "cyanorak role",
                {"code": code, "description": _clean_str(full_role_description(code, self.role_tree))},
            )
            cyr_count += 1
        logger.info(f"MultiCogRoleAnnotationAdapter.get_nodes: {cyr_count} CyanorakRole nodes")

        # 3. TigrRole nodes — only codes present in data
        tigr_codes = self._all_tigr_codes()
        tigr_count = 0
        for code, desc in sorted(tigr_codes.items()):
            yield (
                _tigr_role_node_id(code),
                "tigr role",
                {"code": code, "description": _clean_str(desc)},
            )
            tigr_count += 1
        logger.info(f"MultiCogRoleAnnotationAdapter.get_nodes: {tigr_count} TigrRole nodes")

    def get_edges(self):
        """
        Yield all COG/role edges:
        1. gene→COG category (all strains)
        2. gene→CyanorakRole (Pro/Syn strains; Alteromonas silently yields nothing)
        3. gene→TigrRole (Pro/Syn strains)
        4. CyanorakRole→parent (full hierarchy tree)
        """
        # 1–3. Per-strain gene→* edges
        cog_count = cyr_count = tigr_count = 0
        for adapter in self._strain_adapters:
            for edge in adapter.get_edges():
                label = edge[3]
                if label == "gene_in_cog_category":
                    cog_count += 1
                elif label == "gene_has_cyanorak_role":
                    cyr_count += 1
                elif label == "gene_has_tigr_role":
                    tigr_count += 1
                yield edge

        logger.info(
            f"MultiCogRoleAnnotationAdapter.get_edges: "
            f"{cog_count} gene→COG cat, {cyr_count} gene→CyanorakRole, "
            f"{tigr_count} gene→TigrRole edges"
        )

        # 4. CyanorakRole hierarchy: child → parent
        hier_count = 0
        for code, entry in sorted(self.role_tree.items()):
            parent = entry.get("parent")
            if not parent:
                continue
            yield (
                f"{code}-is_a-{parent}",
                _cyanorak_role_node_id(code),
                _cyanorak_role_node_id(parent),
                "cyanorak_role_is_a_cyanorak_role",
                {},
            )
            hier_count += 1

        logger.info(
            f"MultiCogRoleAnnotationAdapter.get_edges: {hier_count} CyanorakRole hierarchy edges"
        )
