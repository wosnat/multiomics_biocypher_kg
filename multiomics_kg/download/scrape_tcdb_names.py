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


# ── scope ────────────────────────────────────────────────────────────────────

def scoped_families(kept_ids: set[str]) -> dict[str, str]:
    """{family_3part: query_4part} for every family carrying a kept 4/5-part
    ID. Any 4-part in the family works as the result.php query key — the page
    returns the entire family (spec 3b)."""
    out: dict[str, str] = {}
    for tc in kept_ids:
        parts = tc.split(".")
        if len(parts) < 4:
            continue
        fam = ".".join(parts[:3])
        four = ".".join(parts[:4])
        if fam not in out or four < out[fam]:
            out[fam] = four
    return out


# ── fetch (resumable, slow, backoff) ─────────────────────────────────────────

def _http_get(url: str) -> str:
    import requests

    resp = requests.get(url, headers={"User-Agent": USER_AGENT}, timeout=120)
    resp.raise_for_status()
    return resp.text


def fetch_page(url: str, dest: Path, *, force: bool, delay: float,
               retries: int = 3) -> str:
    """Fetch url → dest. Resume = return the cached file when present (unless
    force). Sleeps `delay` BEFORE every network request; exponential backoff
    between retries. Raises after `retries` failures."""
    if dest.exists() and not force:
        return dest.read_text(encoding="utf-8", errors="replace")
    last_exc: Exception | None = None
    for attempt in range(retries):
        time.sleep(delay * (2 ** attempt) if attempt else delay)
        try:
            text = _http_get(url)
        except Exception as exc:  # noqa: BLE001 — retried, then re-raised below
            last_exc = exc
            log.warning(f"  fetch failed (attempt {attempt + 1}/{retries}): {url}: {exc}")
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(text, encoding="utf-8")
        return text
    raise RuntimeError(f"giving up on {url} after {retries} attempts") from last_exc


# ── artifact assembly ────────────────────────────────────────────────────────

def assemble_names(
    browse_names: dict[str, str],
    per_family_results: dict[str, tuple[dict[str, str], dict[str, str]]],
    kept_ids: set[str],
    scraped_families: list[str],
    failed_families: list[str],
) -> dict:
    """Build the artifact dict. Browse layer (levels 0-2) kept in FULL —
    removes gene-set coupling for those levels at ~150 KB. Deep layer (4/5
    part) trimmed to kept IDs (spec §4 scoped-scrape decision)."""
    names: dict[str, dict] = {}
    for tc_id, name in browse_names.items():
        names[tc_id] = {"name": name}

    for _fam, (subfams, systems) in per_family_results.items():
        for tc_id, name in subfams.items():
            if tc_id in kept_ids:
                names[tc_id] = {"name": name}
        for tc_id, full_text in systems.items():
            if tc_id in kept_ids:
                names[tc_id] = {
                    "name": make_name(full_text),
                    "description": make_description(full_text),
                }

    return {
        "meta": {
            "scraped_at": date.today().isoformat(),
            "source": "https://www.tcdb.org (browse.php + search/result.php per family)",
            "scraped_families": sorted(scraped_families),
            "failed_families": sorted(failed_families),
        },
        "names": dict(sorted(names.items())),
    }


# ── main ─────────────────────────────────────────────────────────────────────

def main(force: bool = False, delay: float = DEFAULT_DELAY_S,
         families: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    pruned = json.loads(PRUNED_FILE.read_text())
    kept_ids = set(pruned["kept_tcdb_ids"])
    fams = scoped_families(kept_ids)
    if families:
        fams = {f: q for f, q in fams.items() if f in set(families)}
    log.info(f"scope: {len(fams)} families (from {len(kept_ids)} kept TC ids)")

    browse_html = fetch_page(BROWSE_URL, PAGES_DIR / "browse.html",
                             force=force, delay=delay)
    browse_names = parse_browse(browse_html)
    log.info(f"browse.php: {len(browse_names)} class/subclass/family names")

    per_family: dict[str, tuple[dict[str, str], dict[str, str]]] = {}
    scraped: list[str] = []
    failed: list[str] = []
    consecutive = 0
    for i, (fam, query) in enumerate(sorted(fams.items()), 1):
        dest = PAGES_DIR / f"family_{fam}.html"
        try:
            page = fetch_page(RESULT_URL.format(tc=query), dest,
                              force=force, delay=delay)
        except Exception as exc:  # noqa: BLE001 — per-family isolation
            log.warning(f"[{i}/{len(fams)}] {fam}: FAILED ({exc})")
            failed.append(fam)
            consecutive += 1
            if consecutive >= MAX_CONSECUTIVE_FAILURES:
                log.error(
                    f"{consecutive} consecutive failures — tcdb.org is likely "
                    f"struggling. Aborting; re-run later to resume (cached "
                    f"pages are skipped)."
                )
                sys.exit(1)
            continue
        consecutive = 0
        subfams, systems = parse_family_page(page, family=fam)
        per_family[fam] = (subfams, systems)
        scraped.append(fam)
        if i % 25 == 0 or i == len(fams):
            log.info(f"[{i}/{len(fams)}] fetched+parsed (last: {fam}, "
                     f"{len(subfams)} subfam names, {len(systems)} systems)")

    artifact = assemble_names(browse_names, per_family, kept_ids, scraped, failed)
    NAMES_FILE.parent.mkdir(parents=True, exist_ok=True)
    NAMES_FILE.write_text(json.dumps(artifact, indent=1, ensure_ascii=False) + "\n",
                          encoding="utf-8")
    n = len(artifact["names"])
    log.info(f"wrote {NAMES_FILE} — {n} names, "
             f"{len(scraped)} families scraped, {len(failed)} failed")
    return n


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--force", action="store_true",
                        help="re-fetch pages even when cached")
    parser.add_argument("--delay", type=float, default=DEFAULT_DELAY_S,
                        help="seconds between requests (default 2.5)")
    parser.add_argument("--families", nargs="*", default=None,
                        help="restrict to these 3-part family ids (top-up mode)")
    args = parser.parse_args()
    main(force=args.force, delay=args.delay, families=args.families)
