"""Scrape TCDB node names (T6) into cache/data/tcdb/tcdb_names.json.

Two sources (verified live 2026-08-17 — see
docs/superpowers/specs/2026-08-12-tcdb-node-names-design.md §3):

  1. browse.php — one request; names every class, subclass and family
     (treeview markup: <div rel="<tc id>" class="entry"> +
     <div class="tcid name">).
  2. search/result.php?tc=<4-part> — one request PER FAMILY; the page covers
     the whole family: named-subfamily header rows
     (<strong><A id="X"></A>X:&nbsp;&nbsp;Name</strong> — absent when
     upstream has no name) and every 5-part system's curated Name cell.

Scope: the families carrying kept 4/5-part IDs in the committed
cache/data/tcdb/tcdb_pruned.json (~351 families), plus the full browse layer
(all classes/subclasses/families — decouples levels 0-2 from the gene set).

Hygiene (spec 5d): single-threaded, 2.5 s between requests, resumable page
cache under cache/data/tcdb/raw/pages/ (gitignored), exponential backoff,
abort after 3 consecutive family failures. Re-parsing never re-fetches.

Run:  uv run python -m multiomics_kg.download.scrape_tcdb_names
      [--force] [--delay 2.5] [--families 2.A.1 1.A.11]
"""
from __future__ import annotations

import argparse
import html as html_mod
import json
import logging
import re
import sys
import time
from datetime import date
from pathlib import Path

log = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).parent.parent.parent
TCDB_DIR = PROJECT_ROOT / "cache" / "data" / "tcdb"
PAGES_DIR = TCDB_DIR / "raw" / "pages"
NAMES_FILE = TCDB_DIR / "tcdb_names.json"
PRUNED_FILE = TCDB_DIR / "tcdb_pruned.json"

BROWSE_URL = "https://www.tcdb.org/browse.php"
RESULT_URL = "https://www.tcdb.org/search/result.php?tc={tc}"
USER_AGENT = (
    "multiomics-biocypher-kg/T6-names (academic KG build; "
    "contact: osnat.weissberg@gmail.com)"
)
DEFAULT_DELAY_S = 2.5
MAX_CONSECUTIVE_FAILURES = 3

NAME_CAP = 150
DESCRIPTION_CAP = 400


# ── text cleaning ────────────────────────────────────────────────────────────

def clean_html_fragment(s: str) -> str:
    """Tag-strip + entity-unescape + whitespace-collapse.

    Tolerates a malformed dangling tag at the end of the fragment (one live
    cell carries an unclosed "<br /").
    """
    t = re.sub(r"<[^>]+>", " ", s)
    t = re.sub(r"<[^>]*$", "", t)          # dangling unclosed tag at end
    t = html_mod.unescape(t)
    t = t.replace("\xa0", " ")
    return re.sub(r"\s+", " ", t).strip()


def strip_citations(s: str) -> str:
    """Remove (Author, 19xx/20xx ...) citation parentheticals, iteratively
    (innermost-first, so nested groups collapse over passes)."""
    prev = None
    t = s
    while prev != t:
        prev = t
        t = re.sub(r"\s*\(([^()]*\b(?:19|20)\d{2}[^()]*)\)", "", t)
    t = re.sub(r"\s+", " ", t).strip(" ;,")
    # a stripped trailing citation can leave "text ." — reattach the period
    t = re.sub(r"\s+\.", ".", t)
    return t


_ABBREV_AT_END = re.compile(
    r"(?:\b(?:et al|e\.g|i\.e|sp|spp|subsp|cf|ca|approx|St|var|str|vs)|"
    r"\b[A-Z])\.$"
)


def first_sentence(s: str) -> str:
    """First sentence, not splitting after common abbreviations or
    single-capital genus initials (E. coli)."""
    for match in re.finditer(r"\.\s+(?=[A-Z0-9])", s):
        prefix = s[: match.start() + 1]
        if _ABBREV_AT_END.search(prefix):
            continue
        return prefix
    return s


def cap_at_word_boundary(s: str, limit: int) -> str:
    if len(s) <= limit:
        return s
    cut = s[:limit].rsplit(" ", 1)[0].rstrip(" ,;:.")
    return cut + "…"


def make_name(full_text: str) -> str:
    return cap_at_word_boundary(first_sentence(strip_citations(full_text)), NAME_CAP)


def make_description(full_text: str) -> str:
    return cap_at_word_boundary(full_text, DESCRIPTION_CAP)


# ── parsers ──────────────────────────────────────────────────────────────────

_BROWSE_ENTRY = re.compile(
    r'<div rel="(\d[^"]*)" class="entry"[^>]*>\s*'
    r'<div class="tcid name"[^>]*>([^<]*)</div>',
)
_SUBFAMILY_HEADER = re.compile(
    r'<strong><A id="(\d\.[A-Z]\.\d+\.\d+)"></A>(.*?)</strong>', re.S
)
_SYSTEM_ROW = re.compile(
    r'HREF="/search/result\.php\?tc=(\d\.[A-Z]\.\d+\.\d+\.\d+)">\s*'
    r"[\d.A-Z]+\s*</A>.*?<td><div class='[0-9a-f]+'>(.*?)</div></td>",
    re.S,
)
_RESULT_CLUSTER_START = re.compile(r'<table id="result-cluster"')


def parse_browse(html: str) -> dict[str, str]:
    """browse.php → {tc_id: name} for classes, subclasses and families."""
    out: dict[str, str] = {}
    for tc_id, raw in _BROWSE_ENTRY.findall(html):
        name = clean_html_fragment(raw)
        # strip the "1.A.1:" prefix the cell repeats
        name = re.sub(rf"^{re.escape(tc_id)}\s*:\s*", "", name).strip()
        if name:
            out[tc_id] = name
    return out


def parse_family_page(
    html: str, family: str | None = None
) -> tuple[dict[str, str], dict[str, str]]:
    """result.php?tc=<4-part> → (subfamily_names, system_full_texts).

    Scoped from the result-cluster table start to end-of-page — the table
    embeds an inner </table> mid-way, so a non-greedy table match silently
    truncates (observed live: 2.A.1 lost subfamilies .71-.91). When `family`
    (3-part id) is given, matches are filtered to that family — the page
    embeds off-family 5-part IDs in a sidebar module. A subfamily with no
    <strong> header row has NO upstream name (spec 6d) — absence is the
    signal.
    """
    start = _RESULT_CLUSTER_START.search(html)
    scope = html[start.start():] if start else html
    prefix = f"{family}." if family else None

    subfams: dict[str, str] = {}
    for tc_id, raw in _SUBFAMILY_HEADER.findall(scope):
        if prefix and not tc_id.startswith(prefix):
            continue
        name = clean_html_fragment(raw)
        name = re.sub(rf"^{re.escape(tc_id)}\s*:\s*", "", name).strip()
        if name:
            subfams[tc_id] = name

    systems: dict[str, str] = {}
    for tc_id, raw in _SYSTEM_ROW.findall(scope):
        if prefix and not tc_id.startswith(prefix):
            continue
        text = clean_html_fragment(raw)
        if text:
            systems[tc_id] = text
    return subfams, systems
