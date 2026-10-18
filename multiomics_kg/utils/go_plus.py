"""go-plus parser: GO transport terms -> ChEBI substrates (RO:0004009 'has primary input')."""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import yaml

from multiomics_kg.utils.go_utils import read_go_release

_OBO = "http://purl.obolibrary.org/obo/"
_HAS_PRIMARY_INPUT = _OBO + "RO_0004009"
_TRANSPORT_ROOTS = {_OBO + "GO_0005215", _OBO + "GO_0006810"}  # transporter activity, transport
_VERSION_INFO = "http://www.w3.org/2002/07/owl#versionInfo"
_FILLER = re.compile(r"/(CHEBI|GOCHE)_(\d+)$")
_GO = re.compile(r"/GO_(\d+)$")

CACHE_NAME = "go_transport_substrates.json"


def _transport_closure(graph: dict) -> set[str]:
    children: dict[str, list[str]] = {}
    for e in graph.get("edges", []):
        if e.get("pred") == "is_a":
            children.setdefault(e["obj"], []).append(e["sub"])
    seen, stack = set(_TRANSPORT_ROOTS), list(_TRANSPORT_ROOTS)
    while stack:
        for c in children.get(stack.pop(), []):
            if c not in seen:
                seen.add(c)
                stack.append(c)
    return seen


def parse_transport_axioms(graph: dict) -> dict[str, list[str]]:
    transport = _transport_closure(graph)
    out: dict[str, set[str]] = {}
    for ax in graph.get("logicalDefinitionAxioms", []):
        m = _GO.search(ax.get("definedClassId", ""))
        if not m or not (set(ax.get("genusIds") or []) & transport):
            continue
        for r in ax.get("restrictions") or []:
            if r.get("propertyId") != _HAS_PRIMARY_INPUT:
                continue
            f = _FILLER.search(r.get("fillerId", ""))
            if f:
                out.setdefault(f"GO:{m.group(1)}", set()).add(f"CHEBI:{f.group(2)}")
    return {k: sorted(v) for k, v in sorted(out.items())}


def go_plus_release(graph: dict) -> str:
    for p in graph.get("meta", {}).get("basicPropertyValues", []):
        if p.get("pred") == _VERSION_INFO:
            return p["val"]
    raise ValueError("go-plus graph carries no owl:versionInfo")


def load_target_filter(path: Path) -> dict[str, str]:
    data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    return {t["chebi"]: t["name"] for t in data["targets"]}


def build_go_transport_cache(cache_root: Path) -> dict:
    cache_root = Path(cache_root)
    d = cache_root / "go_terms"
    raw, out = d / "raw" / "go-plus.json", d / CACHE_NAME
    ns_release = read_go_release(cache_root)
    if not raw.exists():
        if ns_release is None:
            from multiomics_kg.download.refresh_go_release import refresh_go_release
            refresh_go_release(cache_root)
            ns_release = read_go_release(cache_root)
        elif out.exists():
            committed = json.loads(out.read_text(encoding="utf-8"))
            if committed.get("release") == ns_release:
                return committed
        if not raw.exists():
            raise RuntimeError(
                f"go-plus raw file missing and committed {CACHE_NAME} does not match GO release "
                f"{ns_release}: run refresh_go_release")
    graph = json.loads(raw.read_text(encoding="utf-8"))["graphs"][0]
    rel = go_plus_release(graph)
    if rel != ns_release:
        raise RuntimeError(f"go-plus release {rel} != namespace cache release {ns_release}: "
                           "run refresh_go_release")
    cache = {"release": rel, "md5": hashlib.md5(raw.read_bytes()).hexdigest(),
             "axioms": parse_transport_axioms(graph)}
    out.write_text(json.dumps(cache, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return cache
