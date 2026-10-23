"""TCDB ontology adapter.

Yields (TCDB cleanup, spec docs/superpowers/specs/2026-10-22-tcdb-cleanup-design.md §5):
- TcdbFamily nodes for every id in tcdb_pruned.json's ``kept_tcdb_ids`` (gene-anchored
  seeds + ancestors, and range-only substrate-bearing systems under kept families),
  each carrying ``direction`` / ``direction_basis`` (TCDB-stated, ``tcdb_direction``).
- Tcdb_family_is_a_tcdb_family parent edges within the kept hierarchy.
- Gene_has_tcdb_family edges: ONE per (gene, TC node, ``tcdb_call``) — eggNOG ids,
  diamond tier calls and diamond best-hit systems are separate edges (no merge),
  each stamped with ``tcdb_gate`` (direct | annotated | below_gate).
- Tcdb_family_transports_metabolite / _chemical_role edges from ``system_substrates`` /
  ``system_roles``: tc_system sources only, no properties.
- Pfam / GO bridge edges on all kept nodes.

Two-class shape mirrors functional_annotation_adapter.MultiPfamAnnotationAdapter.
"""
from __future__ import annotations

import csv
import json
import logging
from collections import Counter
from pathlib import Path
from typing import Iterator

from multiomics_kg.utils.controlled_vocab import VOCAB
from multiomics_kg.utils.curie_utils import normalize_curie
from multiomics_kg.utils.tcdb_diamond import best_hit_tier
from multiomics_kg.utils.tcdb_utils import tcdb_direction

logger = logging.getLogger(__name__)

# tcdb_call → (sources, evidence rung). Shared evidence ladder:
# curated > signature > homology > family_inferred > domain_inferred.
# eggNOG is orthology transfer; both diamond kinds are direct blastp hits vs. the
# curated TCDB FASTA (strength within `homology` is read from `tier`).
_CALL_PROVENANCE = {
    "eggnog": (["eggnog"], "family_inferred"),
    "diamond_tier": (["tcdb_diamond"], "homology"),
    "diamond_best_hit": (["tcdb_diamond"], "homology"),
}
_TC_SYSTEM_LEVEL = 4


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


def tcdb_gate(tcdb_call: str, props: dict, target_level: int) -> str:
    """Which downstream rule a gene edge passes (spec §5.3), computed AFTER the
    seed_aliases remap so ``target_level`` is the level of the node it lands on.

    direct: rank-1 best hit with tier <= 2; eggNOG; tier-1 call — and in every arm
        the landed node must be a tc_system (direct reaches system substrates).
    annotated: a direct arm that landed on a non-system ancestor (alias remap);
        any other eggNOG edge; tier-2 call.
    below_gate: everything else.
    """
    tier = props.get("tier")
    if tcdb_call == "eggnog":
        passes, direct_arm = True, True
    elif tcdb_call == "diamond_tier":
        passes, direct_arm = tier in (1, 2), tier == 1
    elif tcdb_call == "diamond_best_hit":
        passes = props.get("rank") == 1 and tier is not None and tier <= 2
        direct_arm = passes
    else:
        passes = direct_arm = False
    if not passes:
        return "below_gate"
    return "direct" if direct_arm and target_level == _TC_SYSTEM_LEVEL else "annotated"


def _float_props(src: dict, keys: tuple[str, ...]) -> dict:
    return {k: float(src[k]) for k in keys if src.get(k) is not None}


class TcdbAnnotationAdapter:
    """Per-strain adapter: yields raw ``(locus_tag, tcdb_call, tc_id, props)`` claims.

    eggNOG claims come from the merged ``tcdb_eggnog_ids``; diamond claims are read
    from the Phase-1 ``<strain>.tcdb.calls.json`` directly (interpro_adapter pattern):
    one ``diamond_tier`` claim per entry of ``calls`` and one ``diamond_best_hit``
    claim per system in ``best_hit.hits`` — for EVERY protein with a best hit,
    including those with no tier call. A strain without calls.json yields eggNOG
    claims only. Alias remap, pruning, collapse and gating happen in the Multi wrapper.
    """

    def __init__(self, genome_dir: Path, test_mode: bool = False) -> None:
        self.genome_dir = Path(genome_dir)
        self.test_mode = test_mode
        self._genes: dict = {}
        # {protein_id: calls.json record}
        self._diamond: dict[str, dict] = {}
        self._load()

    def _load(self) -> None:
        from multiomics_kg.utils.annotations_cache import load_merged_annotations
        self._genes = load_merged_annotations(self.genome_dir)
        self._diamond = self._load_diamond_records()

    def _load_diamond_records(self) -> dict[str, dict]:
        path = self.genome_dir / "tcdb" / f"{self.genome_dir.name}.tcdb.calls.json"
        if not path.exists():
            return {}
        try:
            data = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            logger.warning(f"TcdbAnnotationAdapter({self.genome_dir.name}): unreadable {path}: {exc}")
            return {}
        out = {str(wp).strip(): rec for wp, rec in data.items() if isinstance(rec, dict)}
        no_hits = sum(1 for rec in out.values()
                      if isinstance(rec.get("best_hit"), dict) and rec["best_hit"].get("tcids")
                      and not rec["best_hit"].get("hits"))
        if no_hits:
            logger.warning(
                f"TcdbAnnotationAdapter({self.genome_dir.name}): {no_hits} best_hit block(s) "
                f"without per-system `hits` — regenerate with run_tcdb_diamond.py --from-tsv; "
                f"their best-hit edges are skipped")
        return out

    def _claims_for(self, gene: dict) -> Iterator[tuple[str, str, dict]]:
        for tc in gene.get("tcdb_eggnog_ids") or []:
            if tc:
                yield "eggnog", tc, {}
        rec = self._diamond.get((gene.get("protein_id") or "").strip())
        if not rec:
            return
        for cand in rec.get("calls") or []:
            tc = cand.get("tcid")
            if not tc:
                continue
            props: dict = {}
            if cand.get("tier") is not None:
                props["tier"] = int(cand["tier"])
            props.update(_float_props(cand, ("identity", "qcov", "scov", "evalue")))
            yield "diamond_tier", tc, props
        best = rec.get("best_hit")
        hits = (best or {}).get("hits") or [] if isinstance(best, dict) else []
        if not hits:
            return
        top = hits[0].get("bitscore")
        for i, hit in enumerate(hits):
            tc = hit.get("tcid")
            if not tc:
                continue
            props = {"rank": i + 1}
            if hit.get("identity") is not None:
                props["tier"] = best_hit_tier(float(hit["identity"]))
            props.update(_float_props(hit, ("identity", "qcov", "scov", "evalue", "bitscore")))
            if hit.get("bitscore") is not None and top:
                props["bitscore_ratio"] = (
                    1.0 if i == 0 else round(float(hit["bitscore"]) / float(top), 4))
            yield "diamond_best_hit", tc, props

    def get_edges(self) -> Iterator[tuple[str, str, str, dict]]:
        count = 0
        for locus_tag, gene in self._genes.items():
            for tcdb_call, tc, props in self._claims_for(gene):
                yield locus_tag, tcdb_call, tc, props
                count += 1
                if self.test_mode and count >= 100:
                    return
        logger.debug(f"TcdbAnnotationAdapter({self.genome_dir.name}): {count} TCDB claims")


class MultiTcdbAnnotationAdapter:
    """Multi-strain orchestrator: owns TcdbFamily nodes, hierarchy / gene / substrate /
    bridge edges."""

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
        self._system_substrates: dict[str, list[str]] = {}
        self._system_roles: dict[str, list[str]] = {}
        self._chemical_roles: dict[str, dict] = {}
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
        if "system_substrates" not in pruned:
            logger.warning(
                "tcdb_pruned.json has no system_substrates (pre-cleanup format) — no TCDB "
                "substrate edges; re-run prepare_data step 6")
        # {tc_system: [canonical id, ...]} — each system's OWN TCDB substrates.
        self._system_substrates = pruned.get("system_substrates", {}) or {}
        self._system_roles = pruned.get("system_roles", {}) or {}
        # Remaps gene-annotated TCIDs not in the curated hierarchy to the nearest
        # curated ancestor (e.g. retired `3.A.1.35` → family `3.A.1`); "" = dropped.
        self._seed_aliases = pruned.get("seed_aliases", {}) or {}
        # Cross-ontology bridges, pre-rolled onto kept nodes by step 6.
        # {tc_id: {xref_id: [curated 5-part TCIDs]}}
        self._pfam_bridge = pruned.get("pfam_bridge", {}) or {}
        self._go_bridge = pruned.get("go_bridge", {}) or {}
        self._chemical_roles = pruned.get("chemical_roles", {}) or {}

    def _level(self, tcdb_id: str) -> int:
        level = self._hierarchy.get(tcdb_id, {}).get("level")
        return int(level) if level is not None else tcdb_id.count(".")

    def _is_kept_system(self, tcdb_id: str) -> bool:
        return tcdb_id in self._kept_ids and self._level(tcdb_id) == _TC_SYSTEM_LEVEL

    def get_nodes(self) -> Iterator[tuple[str, str, dict]]:
        if not self._kept_ids:
            self.download_data(cache=self.cache)

        # NOT capped under test_mode, matching every sibling ontology adapter: the
        # ontology is bounded reference data, and capping nodes would leave gene /
        # parent / substrate / bridge edges dangling.
        emit_count = 0
        directions: Counter = Counter()
        for tcdb_id in sorted(self._kept_ids):
            entry = self._hierarchy.get(tcdb_id, {})
            level = self._level(tcdb_id)
            level_kind = entry.get("level_kind", "tc_class")
            # Falls back to the tcdb_id when TCDB has no name for the node.
            display_name = (entry.get("name") or "") or tcdb_id
            direction, basis = tcdb_direction(tcdb_id, self._hierarchy)
            directions[direction] += 1
            props = {
                "name": _clean_str(display_name),
                "tcdb_id": tcdb_id,
                "level": level,
                "level_kind": VOCAB.check("TcdbFamily", "level_kind", level_kind),
                "direction": VOCAB.check("TcdbFamily", "direction", direction),
            }
            if basis is not None:
                props["direction_basis"] = VOCAB.check("TcdbFamily", "direction_basis", basis)
            if entry.get("superfamily"):
                props["superfamily"] = _clean_str(entry["superfamily"])
            if entry.get("description"):
                props["description"] = _clean_str(entry["description"])
            yield _tcdb_node_id(tcdb_id), "tcdb family", props
            emit_count += 1
        logger.info(
            f"MultiTcdbAnnotationAdapter.get_nodes: {emit_count} TcdbFamily nodes "
            f"(direction {dict(sorted(directions.items()))})")
        for role_id, r in sorted(self._chemical_roles.items()):
            yield (role_id, "chemical role", {
                "name": _clean_str(r["name"]), "level": int(r["level"]),
                "level_kind": VOCAB.check("ChemicalRole", "level_kind", r["level_kind"]),
                "all_identifiers": [_clean_str(i) for i in r["all_identifiers"]]})

    def _gene_edges(self) -> Iterator[tuple]:
        """Gene_has_tcdb_family: alias remap → drop unkept → within-kind collapse
        (higher identity wins) → tcdb_gate on the landed node."""
        remapped = dropped = collapsed = 0
        by_call: Counter = Counter()
        by_gate: Counter = Counter()
        for adapter in self._strain_adapters:
            # {(locus_tag, target, tcdb_call): props}, insertion-ordered
            edges: dict[tuple[str, str, str], dict] = {}
            for locus_tag, tcdb_call, tc, props in adapter.get_edges():
                target = tc
                if target not in self._kept_ids:
                    anchor = self._seed_aliases.get(tc) or ""
                    if not anchor or anchor not in self._kept_ids:
                        dropped += 1
                        continue
                    target = anchor
                    remapped += 1
                key = (locus_tag, target, tcdb_call)
                prev = edges.get(key)
                if prev is not None:
                    collapsed += 1
                    if props.get("identity", -1.0) <= prev.get("identity", -1.0):
                        continue
                edges[key] = props
            for (locus_tag, target, tcdb_call), raw in edges.items():
                sources, evidence = _CALL_PROVENANCE[tcdb_call]
                gate = tcdb_gate(tcdb_call, raw, self._level(target))
                props = {
                    "sources": list(sources),
                    "evidence": evidence,
                    "tcdb_call": VOCAB.check("Gene_has_tcdb_family", "tcdb_call", tcdb_call),
                    "tcdb_gate": VOCAB.check("Gene_has_tcdb_family", "tcdb_gate", gate),
                    **raw,
                }
                by_call[tcdb_call] += 1
                by_gate[gate] += 1
                yield (
                    f"{locus_tag}-has_tcdb-{target}-{tcdb_call}",
                    _gene_node_id(locus_tag),
                    _tcdb_node_id(target),
                    "gene_has_tcdb_family",
                    props,
                )
        logger.info(
            f"MultiTcdbAnnotationAdapter: gene→TCDB edges by tcdb_call "
            f"{dict(sorted(by_call.items()))}, by tcdb_gate {dict(sorted(by_gate.items()))}; "
            f"{remapped} re-anchored via seed_aliases, {collapsed} within-kind duplicate(s) "
            f"collapsed, {dropped} dropped (target not kept)")
        self.last_gene_edge_counts = {"by_call": dict(by_call), "by_gate": dict(by_gate),
                                      "remapped": remapped, "collapsed": collapsed,
                                      "dropped": dropped}

    def get_edges(self):
        if not self._kept_ids:
            self.download_data(cache=self.cache)

        # 1. Parent edges within the kept hierarchy
        parent_count = 0
        for tcdb_id in sorted(self._kept_ids):
            parent = self._hierarchy.get(tcdb_id, {}).get("parent")
            if parent and parent in self._kept_ids:
                yield (
                    f"{tcdb_id}-parent-{parent}",
                    _tcdb_node_id(tcdb_id),
                    _tcdb_node_id(parent),
                    "tcdb_family_is_a_tcdb_family",
                    {},
                )
                parent_count += 1

        # 2. Gene→TcdbFamily edges, one per (gene, node, tcdb_call).
        yield from self._gene_edges()

        # 3. Substrate edges: tc_system sources only, each system's own TCDB
        # substrates (no roll-up; a family's range is the walk over its systems).
        sub_count = 0
        skipped = 0
        for tcdb_id, primary_ids in sorted(self._system_substrates.items()):
            if not self._is_kept_system(tcdb_id):
                skipped += 1
                continue
            for primary in primary_ids:
                yield (f"{tcdb_id}-transports-{primary}", _tcdb_node_id(tcdb_id), primary,
                       "tcdb_family_transports_metabolite", {})
                sub_count += 1

        # Chemical-role hierarchy + role-transport edges (same system-only rule).
        # Role ids are used verbatim: normalize_curie returns None for chebi.role:.
        for role_id, r in sorted(self._chemical_roles.items()):
            for parent in r["parents"]:
                yield (f"{role_id}-isa-{parent}", role_id, parent,
                       "chemical_role_is_a_chemical_role", {})
        role_count = 0
        for tcdb_id, roles in sorted(self._system_roles.items()):
            if not self._is_kept_system(tcdb_id):
                skipped += 1
                continue
            for role_id in roles:
                yield (f"t2r:{tcdb_id}:{role_id}", _tcdb_node_id(tcdb_id), role_id,
                       "tcdb_family_transports_chemical_role", {})
                role_count += 1
        if skipped:
            logger.warning(
                f"MultiTcdbAnnotationAdapter: skipped {skipped} system_substrates/system_roles "
                f"key(s) that are not kept tc_system nodes")

        # 4. Cross-ontology bridges (TcdbFamily → Pfam / GO), on ALL kept nodes —
        # gene-anchored and range-only alike (spec §5.5).
        #
        # DIRECTION IS SEMANTIC, not incidental. These edges assert "TCDB's curated
        # reference proteins for this transport family carry this domain / GO
        # annotation" — a statement about the composition of the TC family. Read
        # backwards (xref → therefore a transporter) it is only ~31% precise, because
        # these domains also occur outside transport contexts. Hence TcdbFamily is
        # the SOURCE. `curated_tcids` preserves the published 5-part TCIDs.
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
            f"{sub_count} substrate, {role_count} role, "
            f"{pfam_count} Pfam-bridge, {go_count} GO-bridge edges"
        )
