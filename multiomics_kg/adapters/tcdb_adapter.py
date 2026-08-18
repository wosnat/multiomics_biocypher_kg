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
  pair without traversing the hierarchy.

  Note: since the 2026-08-06 pruning cleanup the hierarchy keeps only
  gene-annotated TC IDs + their ancestors, so filtering to
  source.level_kind = 'tc_specificity' selects the 246 specificity nodes genes
  actually annotate — not every specificity node in a reachable family. The
  ancestor rollup still carries the full substrate union either way.

Two-class shape mirrors functional_annotation_adapter.MultiPfamAnnotationAdapter.
"""
from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Iterator

from multiomics_kg.utils.controlled_vocab import VOCAB
from multiomics_kg.utils.curie_utils import normalize_curie

logger = logging.getLogger(__name__)

# NOTE: this module used to carry its own `_TC_CLASS_NAMES` copy as a fallback for
# `tc_class` nodes with an empty name. It was dead code: `build_tcdb_hierarchy`
# (prepare_data step 6) owns the authoritative table and fills every class name,
# so all 7 kept classes arrive named and the branch never fired. It could not even
# have helped the one unnamed class in the hierarchy — TC class `6` — because the
# adapter's copy had no entry for it. Class names live in step 6 alone now.


def _clean_str(value: str | None) -> str:
    if value is None:
        return ""
    return value.replace("'", "^").replace("|", "")


def _gene_node_id(locus_tag: str) -> str:
    """Normalize a locus_tag to the gene node ID format (matches Pfam adapter)."""
    return normalize_curie(f"ncbigene:{locus_tag}") or f"ncbigene_{locus_tag}"


def _tcdb_node_id(tcdb_id: str) -> str:
    return normalize_curie(f"tcdb:{tcdb_id}") or f"tcdb_{tcdb_id}"


def _pfam_node_id(pfam_accession: str) -> str:
    return normalize_curie(f"pfam:{pfam_accession}") or f"pfam_{pfam_accession}"


def _go_node_id(go_id: str) -> str:
    """"GO:0005737" → "go:0005737" (matches functional_annotation_adapter)."""
    return normalize_curie(f"go:{go_id}") or f"go_{go_id.replace(':', '_')}"


# GO namespace → the edge label TCDB-family bridge edges use for that namespace.
_GO_LABEL_TO_EDGE = {
    "biological process": "tcdb_family_involved_in_biological_process",
    "molecular function": "tcdb_family_enables_molecular_function",
    "cellular component": "tcdb_family_located_in_cellular_component",
}


class TcdbAnnotationAdapter:
    """Per-strain adapter: yields Gene_has_tcdb_family edges from gene_annotations_merged.json."""

    def __init__(self, genome_dir: Path, test_mode: bool = False) -> None:
        self.genome_dir = Path(genome_dir)
        self.test_mode = test_mode
        self._genes: dict = {}
        # {protein_id: {tcid: best candidate dict}} from the Phase-1 calls.json
        self._diamond: dict[str, dict[str, dict]] = {}
        self._load()

    def _load(self) -> None:
        from multiomics_kg.utils.annotations_cache import load_merged_annotations
        self._genes = load_merged_annotations(self.genome_dir)
        self._diamond = self._load_diamond_evidence()

    def _load_diamond_evidence(self) -> dict[str, dict[str, dict]]:
        """Read `<strain>.tcdb.calls.json` for per-call edge evidence.

        Mirrors interpro_adapter: the merge carries only the light id list, while
        the adapter reads the Phase-1 artifact directly for the rich per-call
        fields. Keyed by protein_id (WP_), then by tcid. When a protein has
        several candidates collapsing to the same tcid, the highest
        confidence_score wins (`calls` is already sorted descending).
        """
        path = self.genome_dir / "tcdb" / f"{self.genome_dir.name}.tcdb.calls.json"
        if not path.exists():
            return {}
        try:
            data = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            logger.warning(f"TcdbAnnotationAdapter({self.genome_dir.name}): unreadable {path}: {exc}")
            return {}
        out: dict[str, dict[str, dict]] = {}
        for wp, rec in data.items():
            if not isinstance(rec, dict):
                continue
            per_tc: dict[str, dict] = {}
            for cand in rec.get("calls", []):
                tcid = cand.get("tcid")
                if tcid and tcid not in per_tc:
                    per_tc[tcid] = cand
            if per_tc:
                out[str(wp).strip()] = per_tc
        return out

    def get_all_tcdb_ids(self) -> set[str]:
        ids: set[str] = set()
        for gene in self._genes.values():
            for tc in gene.get("transporter_classification") or []:
                if tc:
                    ids.add(tc)
        return ids

    def get_edges(self):
        """One Gene_has_tcdb_family edge per (gene, TC id), carrying source provenance.

        `transporter_classification` is the UNION of two independent sources, so a
        TC called by both yields ONE edge with sources=['eggnog','tcdb_diamond'] — that
        agreement is the cross-source corroboration signal, and materialising it as
        two parallel edges would destroy it.

        Diamond evidence (tier / confidence_score / identity / qcov / evalue /
        consensus_n) rides only on edges diamond actually called; eggNOG-only edges
        carry `sources` alone. Properties are sparse rather than null-filled,
        matching Gene_has_interpro_entry's treatment of nullable evalue.
        """
        count = 0
        for locus_tag, gene in self._genes.items():
            egn = set(gene.get("tcdb_eggnog_ids") or [])
            dia = set(gene.get("tcdb_diamond_ids") or [])
            protein_id = (gene.get("protein_id") or "").strip()
            evidence = self._diamond.get(protein_id, {})
            for tc in gene.get("transporter_classification") or []:
                if not tc:
                    continue
                # SORTED so `sources` has one canonical form. The seed_alias merge
                # in MultiTcdbAnnotationAdapter unions and sorts, so an unsorted
                # build here would leave two spellings of the same value
                # (['eggnog','tcdb_diamond'] vs ['tcdb_diamond','eggnog']) and any consumer
                # doing an equality check would silently miss one of them.
                sources = sorted(
                    ({"eggnog"} if tc in egn else set())
                    | ({"tcdb_diamond"} if tc in dia else set())
                )
                if not sources:
                    # Present in the union but attributable to neither per-source
                    # field — a stale merged file. Keep the edge (no annotation is
                    # lost) but leave provenance empty rather than guessing.
                    logger.debug(
                        f"TcdbAnnotationAdapter({self.genome_dir.name}): {locus_tag} "
                        f"TC {tc} has no source attribution; re-run prepare_data step 2"
                    )
                props: dict = {"sources": sources}
                cand = evidence.get(tc) if "tcdb_diamond" in sources else None
                if cand:
                    for key in ("tier", "consensus_n"):
                        if cand.get(key) is not None:
                            props[key] = int(cand[key])
                    for key in ("confidence_score", "identity", "qcov", "evalue"):
                        if cand.get(key) is not None:
                            props[key] = float(cand[key])
                yield (
                    f"{locus_tag}-has_tcdb-{tc}",
                    _gene_node_id(locus_tag),
                    _tcdb_node_id(tc),
                    "gene_has_tcdb_family",
                    props,
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
        pfam_node_ids: set[str] | None = None,
        go_terms: dict[str, str] | None = None,
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
        self._pfam_bridge: dict[str, dict[str, list[str]]] = {}
        self._go_bridge: dict[str, dict[str, list[str]]] = {}
        # None = endpoint node set not provided → emit NO bridge edges of that kind
        # (we cannot guarantee the endpoints exist). A provided set (even empty)
        # prunes to it. Mirrors MultiInterproAnnotationAdapter.pfam_node_ids.
        self.pfam_node_ids = (
            None if pfam_node_ids is None
            else {p.split(".")[0] for p in pfam_node_ids}
        )
        self.go_terms = go_terms

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
        # Cross-ontology bridges, pre-rolled onto kept nodes by step 6.
        # {tc_id: {xref_id: [curated 5-part TCIDs]}}
        self._pfam_bridge = pruned.get("pfam_bridge", {}) or {}
        self._go_bridge = pruned.get("go_bridge", {}) or {}

    def _compute_substrate_depth(self) -> dict[str, set[str]]:
        """{tc_id: {primary_id it is the MOST-SPECIFIC kept node for}}.

        A kept node is *most_specific* for a substrate when no kept child of it
        also carries that substrate. Checking DIRECT children is sufficient — never
        the whole subtree — because `_prune_tcdb` only ever walks *up* from a
        seed, so the kept set is ancestor-closed: any kept descendant of `t` has
        a kept child-of-`t` on its path to `t`, and the rollup is transitive, so
        that child necessarily carries the substrate too.
        """
        children: dict[str, list[str]] = {}
        for tc in self._kept_ids:
            parent = self._hierarchy.get(tc, {}).get("parent")
            if parent is not None and parent in self._kept_ids:
                children.setdefault(parent, []).append(tc)

        depth: dict[str, set[str]] = {}
        for tc, primaries in self._subtree_substrates.items():
            covered: set[str] = set()
            for child in children.get(tc, []):
                covered.update(self._subtree_substrates.get(child, ()))
            depth[tc] = set(primaries) - covered
        return depth

    def get_nodes(self) -> Iterator[tuple[str, str, dict]]:
        if not self._kept_ids:
            self.download_data(cache=self.cache)

        # NOT capped under test_mode, matching every sibling ontology adapter
        # (cazy / interpro / psortb / signalp all cap the per-gene EDGE loop and
        # emit the full ontology). The ontology is bounded reference data — the
        # expensive part is the 53K gene edges, which the per-strain adapter caps.
        # Capping nodes here left gene/parent/substrate/bridge edges pointing at
        # the ~1,400 unemitted nodes, and `skip_bad_relationships: true` dropped
        # them silently, so `--test` produced a quietly broken TCDB layer.
        emit_count = 0
        for tcdb_id in sorted(self._kept_ids):
            entry = self._hierarchy.get(tcdb_id, {})
            level = entry.get("level", 0)
            level_kind = entry.get("level_kind", "tc_class")
            # Falls back to the tcdb_id when the source has no name. TCDB ships
            # descriptions only for classes and families, so subclass /
            # subfamily / specificity nodes all render as their bare id.
            display_name = (entry.get("name") or "") or tcdb_id

            props = {
                "name": _clean_str(display_name),
                "tcdb_id": tcdb_id,
                "level": level,
                "level_kind": level_kind,
            }
            if entry.get("superfamily"):
                props["superfamily"] = _clean_str(entry["superfamily"])
            if entry.get("description"):
                props["description"] = _clean_str(entry["description"])
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
        # Remapping can make two distinct source TCIDs land on the SAME kept node:
        # e.g. eggNOG's retired `3.A.1.35` re-anchors to `3.A.1`, which diamond
        # already called directly. Emitting both would produce parallel edges whose
        # `sources` are split (['eggnog'] on one, ['tcdb_diamond'] on the other) —
        # destroying exactly the corroboration signal the single-edge model exists
        # to capture. So collapse per (gene, TC) after remapping: union `sources`
        # and keep the diamond evidence block from whichever edge carries it.
        merged_count = 0
        for adapter in self._strain_adapters:
            by_pair: dict[tuple[str, str], tuple[str, dict]] = {}
            order: list[tuple[str, str]] = []
            for edge in adapter.get_edges():
                edge_id, source, target, _label, props = edge
                if target not in kept_node_ids:
                    anchor_target = alias_node_ids.get(target)
                    if not (anchor_target and anchor_target in kept_node_ids):
                        dropped += 1
                        continue
                    target = anchor_target
                    remapped += 1
                key = (source, target)
                if key not in by_pair:
                    by_pair[key] = (edge_id, dict(props))
                    order.append(key)
                    continue
                merged_count += 1
                prev_id, prev_props = by_pair[key]
                sources = sorted(set(prev_props.get("sources", [])) | set(props.get("sources", [])))
                # Prefer whichever side carries diamond evidence; if both do, the
                # higher-confidence call wins.
                base = prev_props
                if "tier" in props and (
                    "tier" not in prev_props
                    or props.get("confidence_score", 0.0) > prev_props.get("confidence_score", 0.0)
                ):
                    base = props
                combined = dict(base)
                combined["sources"] = sources
                # Keep the lexicographically smaller id so the merge is deterministic.
                by_pair[key] = (min(prev_id, edge_id), combined)
            for key in order:
                edge_id, props = by_pair[key]
                yield (edge_id, key[0], key[1], "gene_has_tcdb_family", props)
                gene_count += 1
        if merged_count:
            logger.info(
                f"MultiTcdbAnnotationAdapter: merged {merged_count} gene→TCDB edge(s) "
                f"that collapsed onto the same (gene, TC) after seed_alias remapping"
            )
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
        # (not just gene-annotated ones).
        #
        # `substrate_depth` marks whether this node is the MOST-SPECIFIC kept
        # node carrying the substrate ('most_specific') or an ancestor of one
        # ('inherited'). Without it, "how many distinct transporter systems move
        # X" has no cheap answer: counting every level double-counts an ancestor
        # with its own descendant, and the old `level_kind = 'tc_specificity'`
        # filter selects only 466 of 11,263 edges — leaving 83% of transported
        # metabolites at transporter_count = 0 after the ancestor-only prune.
        #
        # NOT the same as "curated vs inherited": only tc_specificity nodes carry
        # their own `substrate_classes` and they are leaves, so curated-vs-inherited
        # is exactly `level_kind = 'tc_specificity'` and needs no new property.
        # Depth is the part that is not derivable without a hierarchy traversal.
        substrate_depth = self._compute_substrate_depth()
        sub_count = 0
        deepest_count = 0
        for tcdb_id, primary_ids in sorted(self._subtree_substrates.items()):
            for primary in primary_ids:
                is_deepest = primary in substrate_depth.get(tcdb_id, ())
                if is_deepest:
                    deepest_count += 1
                yield (
                    f"{tcdb_id}-transports-{primary}",
                    _tcdb_node_id(tcdb_id),
                    primary,
                    "tcdb_family_transports_metabolite",
                    # Categorical str, not bool: BioCypher mishandles boolean
                    # properties, so the KG uses string vocabularies throughout.
                    {"substrate_depth": VOCAB.check(
                        "Tcdb_family_transports_metabolite", "substrate_depth",
                        "most_specific" if is_deepest else "inherited")},
                )
                sub_count += 1

        # 4. Cross-ontology bridges (TcdbFamily → Pfam / GO).
        #
        # DIRECTION IS SEMANTIC, not incidental. These edges assert "TCDB's curated
        # reference proteins for this transport family carry this domain / GO
        # annotation" — a statement about the composition of the TC family.
        # Measured on 42 strains: read outward from a gene's known TC family, the
        # Pfam bridge corroborates 85% of calls and contradicts 2%. Read backwards
        # (xref → therefore a transporter) it is only ~31% precise, because these
        # domains also occur outside transport contexts. Hence TcdbFamily is the
        # SOURCE: the sound traversal is the natural one.
        #
        # `curated_tcids` preserves the original published 5-part TCIDs so the
        # roll-up onto surviving ancestors loses no precision.
        pfam_count = 0
        if self.pfam_node_ids is not None:
            for tcdb_id, xrefs in sorted(self._pfam_bridge.items()):
                if tcdb_id not in self._kept_ids:
                    continue
                for pfam_acc, curated in sorted(xrefs.items()):
                    if pfam_acc not in self.pfam_node_ids:
                        continue
                    yield (
                        f"{tcdb_id}-has_pfam-{pfam_acc}",
                        _tcdb_node_id(tcdb_id),
                        _pfam_node_id(pfam_acc),
                        "tcdb_family_has_pfam_domain",
                        {"curated_tcids": [_clean_str(c) for c in curated]},
                    )
                    pfam_count += 1

        go_count = 0
        if self.go_terms is not None:
            for tcdb_id, xrefs in sorted(self._go_bridge.items()):
                if tcdb_id not in self._kept_ids:
                    continue
                for go_id, curated in sorted(xrefs.items()):
                    label = self.go_terms.get(go_id)
                    edge_label = _GO_LABEL_TO_EDGE.get(label or "")
                    if edge_label is None:
                        continue
                    yield (
                        f"{tcdb_id}-has_go-{go_id}",
                        _tcdb_node_id(tcdb_id),
                        _go_node_id(go_id),
                        edge_label,
                        {"curated_tcids": [_clean_str(c) for c in curated]},
                    )
                    go_count += 1

        logger.info(
            f"MultiTcdbAnnotationAdapter.get_edges: {parent_count} parent, "
            f"{gene_count} gene, {sub_count} substrate "
            f"({deepest_count} most_specific / {sub_count - deepest_count} inherited), "
            f"{pfam_count} Pfam-bridge, {go_count} GO-bridge edges"
        )
