"""KEGG modules: the ``KeggModule`` ontology and per-organism module completeness.

Spec: ``docs/superpowers/specs/2026-11-18-kegg-modules-design.md``.

Reads the ``modules`` / ``module_ec_to_kos`` blocks step 6 writes into ``kegg_data.json`` and
each strain's ``gene_annotations_merged.json``, and emits:

* ``KeggModule`` nodes — modules (``level`` 3) and the ``ko00002`` class tree above them
  (type / class / subclass, ``level`` 0-2), one label, ``kegg.module:`` ids;
* ``Kegg_module_is_a_kegg_module`` (child → parent), ``Kegg_term_in_kegg_module`` (KO → module,
  ``steps`` + ``component``), ``Kegg_module_in_kegg_pathway`` (module → kept pathway map);
* ``Organism_has_kegg_module`` — one edge per (organism, module) with ≥ 1 step covered,
  scored by ``multiomics_kg.utils.kegg_modules.evaluate_organism``.

Dangling-proof: KO edges only for KOs the KEGG adapter emits (``ko_node_ids``), pathway edges
only for maps it emits (``pathway_node_ids``), and the module node set is fixed after the
evaluation, so every organism edge lands on a node.
"""
from __future__ import annotations

import csv
import json
import logging
import re
from collections import Counter
from pathlib import Path

from multiomics_kg.utils.curie_utils import normalize_curie
from multiomics_kg.utils.kegg_modules import ModuleResult, evaluate_organism, ko_roles, parse_steps
from multiomics_kg.utils.kegg_utils import load_kegg_data

logger = logging.getLogger(__name__)

_LEVEL_KINDS = ("kegg_module_type", "kegg_module_class", "kegg_module_subclass", "kegg_module")
_TYPE_OF_ROOT = {"Pathway modules": "pathway", "Signature modules": "signature"}


def _clean_str(value: str) -> str:
    return value.replace("'", "^").replace("|", "")


def _slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def _node_id(raw: str) -> str:
    return normalize_curie(f"kegg.module:{raw}") or f"kegg.module:{raw}"


def _class_id(path: list[str]) -> str:
    return _node_id(".".join(_slug(p) for p in path))


def _ko_node_id(ko: str) -> str:
    return normalize_curie(f"kegg.orthology:{ko}") or f"kegg.orthology:{ko}"


def _pathway_node_id(pw: str) -> str:
    return normalize_curie(f"kegg.pathway:{pw}") or f"kegg.pathway:{pw}"


def _organism_id(accession: str) -> str:
    return normalize_curie(f"insdc.gcf:{accession}") or f"insdc.gcf_{accession}"


def _is_full_ec(ec: str) -> bool:
    return bool(re.fullmatch(r"\d+\.\d+\.\d+\.n?\d+", ec))


class MultiKeggModuleAdapter:
    """KeggModule nodes + module edges + per-organism completeness for every registry strain."""

    def __init__(self, genome_config_file: str, cache_root: Path, ko_node_ids: set[str],
                 pathway_node_ids: set[str], test_mode: bool = False) -> None:
        kegg_data = load_kegg_data(Path(cache_root))
        self.modules: dict[str, dict] = kegg_data.get("modules", {})
        if not self.modules:
            logger.warning("MultiKeggModuleAdapter: kegg_data.json has no `modules` block "
                           "(re-run prepare_data step 6); no KeggModule content emitted")
        self.ec_to_kos: dict[str, list[str]] = kegg_data.get("module_ec_to_kos", {})
        self.ko_node_ids = set(ko_node_ids)
        self.pathway_node_ids = set(pathway_node_ids)
        self.test_mode = test_mode
        self.steps = {m: parse_steps(e["definition"]) for m, e in self.modules.items()}
        self.roles = {m: ko_roles(s) for m, s in self.steps.items()}
        self.results: dict[str, dict[str, ModuleResult]] = {}   # organism id → module → result
        self._evaluate(genome_config_file)
        self._module_ids = self._kept_modules()

    # ── evaluation ────────────────────────────────────────────────────────

    def _strains(self, genome_config_file: str):
        with open(genome_config_file, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(ln for ln in fh if not ln.lstrip().startswith("#")))
        seen = set()
        for row in rows:
            data_dir, acc = (row.get("data_dir") or "").strip(), (row.get("ncbi_accession") or "").strip()
            if not data_dir or not acc or data_dir in seen:
                continue
            seen.add(data_dir)
            yield acc, Path(data_dir)

    def _evaluate(self, genome_config_file: str) -> None:
        strains = list(self._strains(genome_config_file))
        if self.test_mode:
            strains = strains[:3]
        for acc, data_dir in strains:
            path = data_dir / "gene_annotations_merged.json"
            if not path.exists():
                logger.warning(f"MultiKeggModuleAdapter: {path} missing; no module edges for {acc}")
                continue
            genes = json.loads(path.read_text(encoding="utf-8"))
            kos: set[str] = set()
            ec_candidates: dict[str, set[str]] = {}
            for locus_tag, g in genes.items():
                live = {k.replace("ko:", "") for k in g.get("kegg_ko") or []} & self.ko_node_ids
                if live:
                    kos |= live
                    continue
                for ec in g.get("ec_numbers") or []:
                    if not _is_full_ec(ec):
                        continue
                    for ko in self.ec_to_kos.get(ec, ()):
                        ec_candidates.setdefault(ko, set()).add(locus_tag)
            res = evaluate_organism(self.steps, kos, ec_candidates)
            self.results[_organism_id(acc)] = {m: r for m, r in res.items() if r.steps_present > 0}

    def _kept_modules(self) -> set[str]:
        with_ko = {m for m, roles in self.roles.items() if set(roles) & self.ko_node_ids}
        with_edge = {m for res in self.results.values() for m in res}
        return with_ko | with_edge

    # ── output ────────────────────────────────────────────────────────────

    def _class_nodes(self) -> dict[str, tuple[list[str], str]]:
        out = {}
        for m in self._module_ids:
            path = self.modules[m].get("class_path") or []
            mtype = self.modules[m]["module_type"]
            for depth in range(1, len(path) + 1):
                out[_class_id(path[:depth])] = (path[:depth], _TYPE_OF_ROOT.get(path[0], mtype))
        return out

    def module_node_ids(self) -> set[str]:
        return {_node_id(m) for m in self._module_ids} | set(self._class_nodes())

    def get_nodes(self):
        counts = Counter()
        for cid, (path, mtype) in sorted(self._class_nodes().items()):
            level = len(path) - 1
            counts[_LEVEL_KINDS[level]] += 1
            yield (cid, "kegg module", {"name": _clean_str(path[-1]), "level": level,
                                        "level_kind": _LEVEL_KINDS[level], "module_type": mtype})
        for m in sorted(self._module_ids):
            e = self.modules[m]
            counts["kegg_module"] += 1
            yield (_node_id(m), "kegg module", {
                "name": _clean_str(e["name"]), "level": 3, "level_kind": "kegg_module",
                "module_type": e["module_type"], "module_id": m,
                "definition": _clean_str(e["definition"]), "step_count": len(self.steps[m])})
        logger.info(f"MultiKeggModuleAdapter.get_nodes: {dict(counts)}")

    def get_edges(self):
        counts = Counter()
        for cid, (path, _) in sorted(self._class_nodes().items()):
            if len(path) > 1:
                counts["is_a"] += 1
                parent = _class_id(path[:-1])
                yield (f"{cid}-is_a-{parent}", cid, parent, "kegg_module_is_a_kegg_module", {})
        for m in sorted(self._module_ids):
            mid = _node_id(m)
            path = self.modules[m].get("class_path") or []
            if path:
                counts["is_a"] += 1
                parent = _class_id(path)
                yield (f"{mid}-is_a-{parent}", mid, parent, "kegg_module_is_a_kegg_module", {})
            for ko, (steps, component) in sorted(self.roles[m].items()):
                if ko in self.ko_node_ids:
                    counts["ko"] += 1
                    yield (f"{ko}-in_module-{m}", _ko_node_id(ko), mid, "kegg_term_in_kegg_module",
                           {"steps": steps, "component": component})
            for pw in self.modules[m].get("pathways") or []:
                if pw in self.pathway_node_ids:
                    counts["pathway"] += 1
                    yield (f"{m}-in_pathway-{pw}", mid, _pathway_node_id(pw),
                           "kegg_module_in_kegg_pathway", {})
        ec_steps = 0
        for org in sorted(self.results):
            for m, r in sorted(self.results[org].items()):
                counts["organism"] += 1
                counts[r.status] += 1
                ec_steps += r.steps_present_via_ec
                props = {"status": r.status, "steps_total": r.steps_total,
                         "steps_present": r.steps_present,
                         "steps_present_via_ec": r.steps_present_via_ec,
                         "missing_step_indices": r.missing_step_indices,
                         "missing_step_definitions": [_clean_str(d) for d in r.missing_step_definitions],
                         "missing_kos": r.missing_kos}
                if r.ec_filled_step_indices:
                    props["ec_filled_step_indices"] = r.ec_filled_step_indices
                    props["ec_genes"] = [_clean_str(g) for g in r.ec_genes]
                yield (f"{org}-has_module-{m}", org, _node_id(m), "organism_has_kegg_module", props)
        logger.info(f"MultiKeggModuleAdapter.get_edges: {dict(counts)}; EC-filled steps {ec_steps}")
