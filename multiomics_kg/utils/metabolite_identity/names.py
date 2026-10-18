"""Name keys, stereo descriptors and formula rules for metabolite identity (spec §3.3, §3.6)."""
from __future__ import annotations

import re

import chemparse

PERIODIC_TABLE = frozenset("""
H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se
Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb
Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm
Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og
""".split())

# A formula-like name: element symbols with optional counts, optional trailing charge.
_FORMULA_NAME = re.compile(r"^(?:[A-Z][a-z]?\d*)+(?:\d*[+-]|[+-]\d*)?$")
_SYMBOL = re.compile(r"[A-Z][a-z]?")
_WS = re.compile(r"\s+")
# Placeholder R: an R not starting one of Ra Rb Re Rf Rg Rh Ru. 'Rn' reads as R + n (spec §3.3).
_PLACEHOLDER_R = re.compile(r"R(?![abefghu])")
_POLYMER = re.compile(r"\)[nmx](?![a-z])")
_DESCRIPTOR = re.compile(
    r"(?:^|(?<=[-\s(,]))"
    r"(?:L-|D-|DL-|meso-"
    r"|\((?:\d*[RS])(?:,\d*[RS])*\)-"
    r"|\([+\-−]\)-"
    r"|(?:alpha|beta|α|β)-(?=[DL]-))")
_ANOMER = re.compile(r"(?:^|(?<=[-\s(,]))(?:alpha|beta|α|β)-", re.IGNORECASE)


def is_formula_like(name: str) -> bool:
    """Case-exact: every [A-Z][a-z]? token is an element symbol (so AMP and GABA are not formulas)."""
    s = name.strip()
    return bool(_FORMULA_NAME.match(s)) and all(sym in PERIODIC_TABLE for sym in _SYMBOL.findall(s))


def name_key(name: str) -> str:
    """Case-fold a name, except formulas (Co vs CO, NO vs No), where case is meaning."""
    s = _WS.sub(" ", name.strip())
    return s if is_formula_like(s) else s.casefold()


def is_stereo_qualified(name: str) -> bool:
    return bool(_DESCRIPTOR.search(name))


# Geometric (cis/trans, E/Z) descriptors: "cis-", "trans-", "all-trans-", "(E)-", "(2Z,4E)-".
# Also locant-prefixed / multiplied forms: "9,9'-dicis-", "9,15,9'-tricis-", "9,9',15-tri-cis-".
_GEOMETRIC = re.compile(
    r"(?:^|(?<=[-\s(,]))(?:(?:[\d',]+-)?(?:all-|(?:di|tri|tetra|poly)-?)?(?:cis|trans)-"
    r"|\((?:\d*[EZ])(?:,\d*[EZ])*\)-)", re.IGNORECASE)


def strip_geometric_tokens(name: str) -> str:
    """Case-folded name without geometric descriptors (spec §3.6: geometric isomers are
    separate compounds, never stereo forms of each other)."""
    return _GEOMETRIC.sub("", name).casefold().strip()


def strip_anomer_tokens(name: str) -> str:
    return _ANOMER.sub("", name).casefold().strip()


def is_polymer_formula(formula: str | None) -> bool:
    return bool(formula) and bool(_POLYMER.search(formula))


def is_generic_formula(formula: str | None) -> bool:
    if not formula:
        return False
    if "*" in formula or _PLACEHOLDER_R.search(formula):
        return True
    return any(sym not in PERIODIC_TABLE for sym in _SYMBOL.findall(formula))


def parse_elements(formula: str | None) -> list[str]:
    if not formula:
        return []
    cleaned = _PLACEHOLDER_R.sub("", formula).replace("*", "")
    cleaned = re.sub(r"\)[nmx]", ")", cleaned)
    cleaned = re.sub(r"[+-]\d*$|\d*[+-]$", "", cleaned)
    try:
        found = chemparse.parse_formula(cleaned).keys()
    except Exception:
        found = _SYMBOL.findall(cleaned)
    return sorted({e for e in found if e in PERIODIC_TABLE})


def normalize_hmdb(raw: str | None) -> str | None:
    if not raw:
        return None
    m = re.fullmatch(r"HMDB(\d+)", raw.strip())
    return f"HMDB{int(m.group(1)):07d}" if m else None
