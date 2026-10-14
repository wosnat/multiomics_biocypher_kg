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
_FORMULA_NAME = re.compile(r"^((?:[A-Za-z]{1,2}\d*)+)(?:\d*[+-]|[+-]\d*)?$")
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


_LOWER_SYMBOLS = frozenset(sym.lower() for sym in PERIODIC_TABLE)


def _segments_ok(letters: str) -> bool:
    """True if the letter run splits into element symbols (case-insensitive; Am+P reads AMP)."""
    ok = [False] * (len(letters) + 1)
    ok[0] = True
    for i in range(1, len(letters) + 1):
        for k in (1, 2):
            if i - k >= 0 and ok[i - k] and letters[i - k:i].lower() in _LOWER_SYMBOLS:
                ok[i] = True
    return ok[-1]


def is_formula_like(name: str) -> bool:
    s = name.strip()
    m = _FORMULA_NAME.match(s)
    if not m:
        return False
    return all(_segments_ok(run) for run in re.findall(r"[A-Za-z]+", m.group(1)))


_STRICT_FORMULA_NAME = re.compile(r"^(?:[A-Z][a-z]?\d*)+(?:\d*[+-]|[+-]\d*)?$")


def _is_case_exact_formula(s: str) -> bool:
    """Formula whose capitalisation is itself valid element notation (Co, CO, NH4+, HF)."""
    return bool(_STRICT_FORMULA_NAME.match(s)) and all(
        sym in PERIODIC_TABLE for sym in _SYMBOL.findall(s))


def name_key(name: str) -> str:
    """Case-fold a name, except case-exact formulas (Co vs CO, NO vs No), where case is meaning.

    The exemption is deliberately stricter than is_formula_like: AMP / GABA only read as formulas
    case-insensitively (Am+P, Ga+Ba), and exempting them would stop 'GABA' matching 'gaba'.
    """
    s = _WS.sub(" ", name.strip())
    return s if _is_case_exact_formula(s) else s.casefold()


def is_stereo_qualified(name: str) -> bool:
    return bool(_DESCRIPTOR.search(name))


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
    cleaned = re.sub(r"(?<=[0-9A-Za-z])n$", "", cleaned)
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
