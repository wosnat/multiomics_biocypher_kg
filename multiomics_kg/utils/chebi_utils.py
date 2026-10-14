"""ChEBI ontology loader for step 6 (metabolite identity).

Reads chebi.obo.gz (release-tagged) into a term dict, answers the lookups the
identity pipeline needs, and builds the committed derived cache
cache/data/chebi/chebi_metabolite_xrefs.json (spec §6.4). Read by step 6 only.
ChEBI is CC BY 4.0.
"""
from __future__ import annotations

import gzip
from contextlib import contextmanager
import logging
import re
from pathlib import Path

import requests

log = logging.getLogger(__name__)

CHEBI_OBO_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz"
ROLE_ROOT = "CHEBI:50906"
ROLE_BRANCHES = frozenset({"CHEBI:24432", "CHEBI:51086", "CHEBI:33232"})
CARBOHYDRATE_ANCESTORS = frozenset({"CHEBI:16646", "CHEBI:63299", "CHEBI:24400", "CHEBI:78616"})
ISOTOPICALLY_MODIFIED = "CHEBI:139358"       # "isotopically modified compound"
ATOM = "CHEBI:33250"                          # "atom"
_MASS_NUMBER = re.compile(r"[A-Za-z]-\d+")   # "sulfur-32 atom", "iron-57 atom"

_REL = {"RO:0018033": "conjugates", "RO:0018034": "conjugates",
        "RO:0018036": "tautomers", "RO:0018039": "enantiomers"}
_PROP = re.compile(r'chemrof:(\S+) "(.*)" xsd:\w+')


def _new_term() -> dict:
    return {"name": None, "charge": None, "alt_ids": [], "kegg": [], "conjugates": [],
            "tautomers": [], "enantiomers": [], "is_a": [], "generalized_formula": None,
            "inchikey": None}


def parse_chebi_obo(path: Path) -> tuple[str, dict[str, dict]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    release, terms, cur = "", {}, None
    with opener(path, "rt", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("data-version:"):
                release = line.split(":", 1)[1].strip()
            if line == "[Term]":
                cur = _new_term()
                continue
            if line.startswith("["):
                cur = None
                continue
            if cur is None or ": " not in line:
                continue
            key, _, val = line.partition(": ")
            if key == "id":
                terms[val] = cur
            elif key == "name":
                cur["name"] = val
            elif key == "alt_id":
                cur["alt_ids"].append(val)
            elif key == "is_a":
                cur["is_a"].append(val.split(" ")[0])
            elif key == "is_obsolete" and val == "true":
                cur["obsolete"] = True
            elif key == "xref" and val.startswith("kegg.compound:"):
                cur["kegg"].append(val.split(" ")[0].split(":", 1)[1])
            elif key == "relationship":
                rel, tgt = val.split(" ")[:2]
                if rel in _REL:
                    cur[_REL[rel]].append(tgt)
            elif key == "property_value":
                m = _PROP.match(val)
                if m:
                    prop, pval = m.groups()
                    if prop == "generalized_empirical_formula":
                        cur["generalized_formula"] = pval
                    elif prop == "charge":
                        cur["charge"] = int(pval)
                    elif prop == "inchi_key_string":
                        cur["inchikey"] = pval
    terms = {k: v for k, v in terms.items() if not v.pop("obsolete", False)}
    for t in terms.values():
        for k in ("alt_ids", "kegg", "conjugates", "tautomers", "enantiomers", "is_a"):
            t[k] = sorted(set(t[k]))
    if not release:
        raise ValueError(f"{path}: no data-version header (ChEBI format drift?)")
    if not terms:
        raise ValueError(f"{path}: parsed 0 terms (ChEBI format drift?)")
    return release, terms


def _generic(formula: str | None) -> bool:
    # Local import keeps chebi_utils free of a package-level cycle.
    from multiomics_kg.utils.metabolite_identity.names import is_generic_formula, is_polymer_formula
    return formula is not None and (is_generic_formula(formula) or is_polymer_formula(formula))


class ChebiIndex:
    def __init__(self, terms: dict[str, dict], release: str):
        self.release = release
        self._terms = terms
        self._alt = {a: p for p, t in terms.items() for a in t.get("alt_ids", [])}
        self._kegg: dict[str, list[str]] = {}
        for cid, t in terms.items():
            for k in t.get("kegg", []):
                self._kegg.setdefault(k, []).append(cid)
        for v in self._kegg.values():
            v.sort(key=lambda c: int(c.split(":")[1]))
        self._anc: dict[str, frozenset[str]] = {}
        self.touched: set[str] = set()
        self.missing: set[str] = set()

    @staticmethod
    def _norm(cid: str) -> str:
        cid = cid.strip()
        local = cid.split(":", 1)[1] if ":" in cid else cid
        return f"CHEBI:{int(local)}" if local.isdigit() else f"CHEBI:{local}"

    def primary(self, cid: str) -> str | None:
        c = self._norm(cid)
        p = c if c in self._terms else self._alt.get(c)
        if p is None:
            self.missing.add(c)
        else:
            self.touched.add(p)
        return p

    def term(self, cid: str) -> dict | None:
        p = self.primary(cid)
        return self._terms.get(p) if p else None

    def kegg_terms(self, kegg_id: str) -> list[str]:
        out = self._kegg.get(kegg_id, [])
        self.touched.update(out)
        return list(out)

    def ancestors(self, cid: str) -> frozenset[str]:
        p = self.primary(cid)
        if p is None:
            return frozenset()
        if p not in self._anc:
            seen: set[str] = set()
            stack = list(self._terms[p].get("is_a", []))
            while stack:
                x = stack.pop()
                if x in seen:
                    continue
                seen.add(x)
                stack.extend(self._terms.get(x, {}).get("is_a", []))
            self._anc[p] = frozenset(seen)
        return self._anc[p]

    def is_role(self, cid: str) -> bool:
        p = self.primary(cid)
        if p is None:
            return False
        t = self._terms[p]
        if "is_role" in t:                      # loaded from the derived cache
            return t["is_role"]
        return p == ROLE_ROOT or ROLE_ROOT in self.ancestors(p)

    def is_class(self, cid: str) -> bool:
        p = self.primary(cid)
        if p is None:
            return False
        t = self._terms[p]
        if "is_class" in t:
            return t["is_class"]
        return t.get("inchikey") is None and (
            t.get("generalized_formula") is None or _generic(t.get("generalized_formula")))

    def is_isotope_variant(self, cid: str) -> bool:
        """A stable-/radio-isotope variant of a molecule or atom (e.g. the analytical standards
        "sulfur-32 atom", "((13)C)urea"): the term is under "isotopically modified compound",
        or it is an atom whose name carries a mass number. Structural, so it also catches
        neutral isotopes the isotopologue name pattern misses; "sulfur atom", "precorrin-4"
        and "plastoquinone-9" are not variants."""
        p = self.primary(cid)
        if p is None:
            return False
        anc = self.ancestors(p)
        if ISOTOPICALLY_MODIFIED in anc:
            return True
        return ATOM in anc and bool(_MASS_NUMBER.search(self._terms[p].get("name") or ""))

    def is_carbohydrate(self, cid: str) -> bool:
        p = self.primary(cid)
        return p is not None and bool((self.ancestors(p) | {p}) & CARBOHYDRATE_ANCESTORS)

    def secondary_ids(self, cid: str) -> list[str]:
        p = self.primary(cid)
        return list(self._terms[p].get("alt_ids", [])) if p else []

    def charge(self, cid: str) -> int | None:
        t = self.term(cid)
        return t.get("charge") if t else None

    def all_terms(self) -> dict[str, dict]:
        return self._terms

    @contextmanager
    def untracked(self):
        """Lookups inside the block do not count as touched / missing, so the committed ChEBI
        cache never depends on baseline-only lookups (previous build's ids)."""
        touched, missing = set(self.touched), set(self.missing)
        try:
            yield self
        finally:
            self.touched, self.missing = touched, missing


def build_chebi_cache(index: ChebiIndex) -> dict:
    terms = index.all_terms()
    keep = {c for c, t in terms.items() if t.get("kegg")} | set(index.touched)
    frontier = list(keep)
    while frontier:
        c = frontier.pop()
        t = terms.get(c)
        if t is None:
            continue
        for nxt in t["is_a"] + t["conjugates"] + t["tautomers"] + t["enantiomers"]:
            if nxt in terms and nxt not in keep:
                keep.add(nxt)
                frontier.append(nxt)
    out = {}
    for c in sorted(keep, key=lambda x: int(x.split(":")[1])):
        t = dict(terms[c])
        t["is_role"] = index.is_role(c)
        t["is_class"] = index.is_class(c)
        out[c] = t
    return {"release": index.release, "terms": out}


def index_from_cache(cache: dict) -> ChebiIndex:
    return ChebiIndex(cache["terms"], cache["release"])


def download_chebi_obo(raw_dir: Path, force: bool = False) -> Path:
    raw_dir.mkdir(parents=True, exist_ok=True)
    dest = raw_dir / "chebi.obo.gz"
    if dest.exists() and not force:
        return dest
    tmp = dest.with_suffix(".gz.part")
    with requests.get(CHEBI_OBO_URL, stream=True, timeout=300) as r:
        r.raise_for_status()
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(1 << 20):
                fh.write(chunk)
    tmp.replace(dest)
    return dest
