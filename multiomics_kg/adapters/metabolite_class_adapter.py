"""MetaboliteClass nodes from KEGG compound BRITE (pruned in step 6; spec 2026-10-13 §3.7, §8.4).

Global data, so no Multi* wrapper. Membership edges are dangling-proofed against the metabolite
ids the metabolism adapter actually emitted (BRITE known_ko_ids precedent)."""
from __future__ import annotations

import json
import logging
from pathlib import Path

from multiomics_kg.utils.controlled_vocab import VOCAB
from multiomics_kg.utils.metabolite_identity.classes import COMPOUND_BRITE_TREES, LEVEL_KINDS, _clean as _clean_brite

logger = logging.getLogger(__name__)


def _clean_str(value: str) -> str:
    return value.replace("'", "^").replace("|", "") if isinstance(value, str) else value


class MetaboliteClassAdapter:
    def __init__(self, cache_root: Path, metabolite_node_ids: set[str], test_mode: bool = False):
        self.cache_root = Path(cache_root)
        self.metabolite_node_ids = set(metabolite_node_ids)
        self.test_mode = test_mode
        self._trees: dict[str, dict] = {}

    def download_data(self, cache: bool = True) -> None:
        for code in COMPOUND_BRITE_TREES:
            p = self.cache_root / "kegg" / f"brite_{code}.json"
            if p.exists():
                self._trees[code] = json.loads(p.read_text())
        if not self._trees:
            raise FileNotFoundError(f"no brite_br08xxx.json under {self.cache_root / 'kegg'}; run step 6")

    def get_nodes(self):
        for code, tree in sorted(self._trees.items()):
            for n in tree["nodes"]:
                yield (n["id"], "metabolite class", {
                    "name": _clean_str(_clean_brite(n["name"])), "tree": tree["tree"], "tree_code": code,
                    "level": int(n["level"]),
                    "level_kind": VOCAB.check("MetaboliteClass", "level_kind", LEVEL_KINDS[n["level"]])})

    def get_edges(self):
        skipped = 0
        for code, tree in sorted(self._trees.items()):
            for n in tree["nodes"]:
                if n["parent"]:
                    yield (f"{n['id']}-isa-{n['parent']}", n["id"], n["parent"],
                           "metabolite_class_is_a_metabolite_class", {})
            for m in tree["members"]:
                if m["metabolite_id"] not in self.metabolite_node_ids:
                    skipped += 1
                    continue
                yield (f"m2c:{m['metabolite_id']}:{m['class_id']}", m["metabolite_id"], m["class_id"],
                       "metabolite_in_metabolite_class", {"sources": ["kegg"], "evidence": "curated"})
        if skipped:
            logger.info(f"metabolite classes: skipped {skipped} memberships to unemitted metabolites")
