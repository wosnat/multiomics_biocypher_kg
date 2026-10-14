"""Stage M: KEGG compound BRITE → MetaboliteClass trees, pruned to written metabolites (spec §3.7, §6.2)."""
from __future__ import annotations

import html
import json
import re
from pathlib import Path

from multiomics_kg.utils import kegg_utils

COMPOUND_BRITE_TREES = {
    "br08001": "compounds_with_biological_roles", "br08002": "lipids",
    "br08003": "phytochemical_compounds", "br08005": "bioactive_peptides",
    "br08006": "endocrine_disrupting_compounds", "br08007": "pesticides",
    "br08009": "natural_toxins", "br08021": "glycosides"}
LEVEL_KINDS = ("kegg_compound_class", "kegg_compound_subclass", "kegg_compound_group", "kegg_compound_subgroup")
BR08001_LEVEL0 = ("Organic acids", "Lipids", "Carbohydrates", "Nucleic acids", "Peptides",
                  "Vitamins and cofactors", "Steroids", "Hormones and transmitters", "Antibiotics")
_LEAF = re.compile(r"^(C\d{5})\s")
_SUFFIX = re.compile(r"\s*\[(?:Fig|[A-Za-z]+:[^\]]*)\]")   # [Fig], [BR:…], [PATH:…], [KO:…], …
_TAG = re.compile(r"<[^>]+>")
_URL = "https://rest.kegg.jp/get/br:{code}/json"


def _clean(name: str) -> str:
    return _SUFFIX.sub("", _TAG.sub("", html.unescape(name))).strip()


def parse_compound_brite(tree_code, raw):
    nodes, members = {}, []

    def walk(children, path, parent):
        for i, c in enumerate(children, 1):
            m = _LEAF.match(c.get("name", ""))
            if m and not c.get("children"):
                if parent:
                    members.append((m.group(1), parent))
                continue
            if len(path) >= len(LEVEL_KINDS):
                raise ValueError(f"{tree_code}: class depth > {len(LEVEL_KINDS)} at {c.get('name')!r}")
            seg = path + [f"{'ABCD'[len(path)]}{i}"]
            cid = f"kegg.brite:{tree_code}." + ".".join(seg)
            nodes[cid] = {"parent": parent, "name": _clean(c["name"]), "level": len(path)}
            walk(c.get("children", []), seg, cid)

    walk(raw.get("children", []), [], None)
    if raw.get("children") and not members:
        raise ValueError(f"{tree_code}: parsed 0 compound members (KEGG BRITE format drift?)")
    return nodes, members


def prune_compound_brite(tree_code, nodes, members, *, canonical_of, written):
    kept_members = set()
    for c, cls in members:
        mid = canonical_of.get(f"kegg.compound:{c}", f"kegg.compound:{c}")
        if mid in written:
            kept_members.add((mid, cls))
    keep = set()
    for _, cls in kept_members:
        while cls and cls not in keep:
            keep.add(cls)
            cls = nodes[cls]["parent"]
    return {
        "tree_code": tree_code, "tree": COMPOUND_BRITE_TREES[tree_code],
        "nodes": [{"id": i, "parent": nodes[i]["parent"], "name": nodes[i]["name"], "level": nodes[i]["level"]}
                  for i in sorted(keep)],
        "members": [{"metabolite_id": m, "class_id": c} for m, c in sorted(kept_members)],
    }


def download_compound_brite(raw_dir: Path, force: bool = False) -> dict[str, dict]:
    raw_dir.mkdir(parents=True, exist_ok=True)
    out = {}
    for code in COMPOUND_BRITE_TREES:
        path = raw_dir / f"br_{code}.json"
        if path.exists() and not force:
            out[code] = json.loads(path.read_text())
            continue
        data = kegg_utils._fetch_json(_URL.format(code=code))
        path.write_text(json.dumps(data))
        out[code] = data
    return out
