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

Every run also refreshes the committed move map cache/data/tcdb/tcdb_moved.json
({old id: [targets]} from the "<old> has been moved to <new>" lines on retired family
pages), read from the page cache. Step 6 reads only that committed file, never the
gitignored pages. `--moved-only` rebuilds it from the cached pages without any network.

Each system entry also carries ``proteins`` (the row's protein cell: gene symbols,
accessions, source organism), which names a system whose name cell is empty.

Run:  uv run python -m multiomics_kg.download.scrape_tcdb_names
      [--force] [--delay 2.5] [--families 2.A.1 1.A.11] [--moved-only] [--reparse]
      (--reparse: offline rebuild of every entry from the cached pages)
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
from typing import Iterable

from multiomics_kg.download.utils.cli import load_genome_rows
from multiomics_kg.utils.gene_id_utils import load_gene_annotations

log = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).parent.parent.parent
TCDB_DIR = PROJECT_ROOT / "cache" / "data" / "tcdb"
PAGES_DIR = TCDB_DIR / "raw" / "pages"
NAMES_FILE = TCDB_DIR / "tcdb_names.json"
MOVED_FILE = TCDB_DIR / "tcdb_moved.json"
PRUNED_FILE = TCDB_DIR / "tcdb_pruned.json"
HIERARCHY_FILE = TCDB_DIR / "tcdb_hierarchy.json"

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

# Inline markup is removed without a space ("Na<sup>+</sup>" -> "Na+", "CO<sub>2</sub>"
# -> "CO2", "<i>P</i>yrimidine" -> "Pyrimidine"); every other tag (br, div, td, p, ...)
# is a word break.
_INLINE_TAG = re.compile(
    r"</?(?:sup|sub|i|em|b|strong|u|span|font|a|small|big|tt|var|sc)\b[^>]*>", re.IGNORECASE)


def clean_html_fragment(s: str) -> str:
    """Tag-strip + entity-unescape + whitespace-collapse + punctuation tidy.

    Tolerates a malformed dangling tag at the end of the fragment (one live
    cell carries an unclosed "<br /").
    """
    t = _INLINE_TAG.sub("", s)
    # A tag starts with a letter, '/' or '!' — a bare '<' is text ("(<1 µM)").
    t = re.sub(r"<[A-Za-z/!][^>]*>", " ", t)
    t = re.sub(r"<[A-Za-z/!][^>]*$", "", t)  # dangling unclosed tag at end
    t = html_mod.unescape(t)
    t = t.replace("\xa0", " ")
    t = re.sub(r"\s+", " ", t).strip()
    t = re.sub(r"([(\[])\s+", r"\1", t)          # "( βETE )" -> "(βETE)"
    return re.sub(r"\s+([)\],;.:])", r"\1", t)   # "text ." / "word ," -> attached


def strip_citations(s: str) -> str:
    """Remove (Author, 19xx/20xx ...) citation parentheticals, iteratively
    (innermost-first, so nested groups collapse over passes)."""
    prev = None
    t = s
    while prev != t:
        prev = t
        t = re.sub(r"\s*\(([^()]*\b(?:19|20)\d{2}\b[^()]*)\)", "", t)
    t = re.sub(r"\s+", " ", t).strip(" ;,")
    # a stripped trailing citation can leave "text ." — reattach the period
    t = re.sub(r"\s+\.", ".", t)
    return t


_ABBREV_AT_END = re.compile(
    r"(?:\b(?:et al|e\.g|i\.e|sp|spp|subsp|cf|ca|approx|St|var|str|vs)|"
    r"\b[A-Z])\.$"
)


def _paren_pairs(s: str) -> tuple[list[tuple[int, int]], list[int]]:
    """Matched ``(``/``)`` index pairs, and the indexes of ``(`` that never close."""
    pairs, opened = [], []
    for i, ch in enumerate(s):
        if ch == "(":
            opened.append(i)
        elif ch == ")" and opened:
            pairs.append((opened.pop(), i))
    return pairs, opened


def first_sentence(s: str) -> str:
    """First sentence, not splitting after common abbreviations or
    single-capital genus initials (E. coli), nor inside a parenthetical that
    closes later (a "(" TCDB never closes does not block the split)."""
    pairs, _ = _paren_pairs(s)
    for match in re.finditer(r"\.\s+(?=[A-Z0-9])", s):
        prefix = s[: match.start() + 1]
        if _ABBREV_AT_END.search(prefix):
            continue
        if any(o < match.start() < c for o, c in pairs):   # inside a parenthetical
            continue
        return prefix
    return s


def cap_at_word_boundary(s: str, limit: int) -> str:
    if len(s) <= limit:
        return s
    cut = s[:limit].rsplit(" ", 1)[0].rstrip(" ,;:.")
    return cut + "…"


_DANGLING_END = re.compile(
    r"\b(?:of|and|or|the|a|an|with|for|to|in|from|by|on|at|as|its|their)$", re.IGNORECASE)


def _drop_unclosed_parenthetical(s: str) -> str:
    """Cut at the first '(' that is never closed (the cap left it open), unless that
    would leave almost nothing or a phrase ending mid-way ("PetABCD of"); then only
    the '(' itself goes."""
    _, opened = _paren_pairs(s)
    if not opened:
        return s
    head = s[: opened[0]].rstrip(" ,;:")
    if len(head) >= 15 and not _DANGLING_END.search(head):
        return head
    return _remove_chars(s, opened)


def _remove_chars(s: str, idxs: list[int]) -> str:
    drop = set(idxs)
    out = "".join(ch for i, ch in enumerate(s) if i not in drop)
    return re.sub(r"\s+", " ", out).strip()


def _shorten(s: str, limit: int) -> str:
    """Drop whole parentheticals, longest first, until ``s`` fits ``limit``; the
    text after a removed parenthetical is kept ("X (long list) exporter, YejABEF")."""
    while len(s) > limit:
        pairs, _ = _paren_pairs(s)
        if not pairs:
            break
        o, c = max(pairs, key=lambda p: p[1] - p[0])
        s = re.sub(r"\s+", " ", (s[:o].rstrip() + " " + s[c + 1:].lstrip())).strip()
        s = re.sub(r"\s+([,;.:])", r"\1", s)
    return s


def make_name(full_text: str) -> str:
    s = first_sentence(strip_citations(full_text))
    _, never_closed = _paren_pairs(s)
    s = _remove_chars(s, never_closed)        # a "(" TCDB never closes is a typo
    s = _shorten(s, NAME_CAP)
    return _drop_unclosed_parenthetical(cap_at_word_boundary(s, NAME_CAP))


def make_description(full_text: str) -> str:
    return cap_at_word_boundary(full_text, DESCRIPTION_CAP)


# ── parsers ──────────────────────────────────────────────────────────────────

_BROWSE_ENTRY = re.compile(
    r'<div rel="(\d[^"]*)" class="entry"[^>]*>\s*'
    r'<div class="tcid name"[^>]*>(.*?)</div>', re.S
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
# The protein cell ("UrtA-E urea porter of <i>Anabaena</i> sp. PCC7120 <br /> UrtA (R) ...")
# is the last cell of the same row.
_PROTEIN_CELL = re.compile(
    r"""<td id="right-border"[^>]*>\s*<div class='[0-9a-f]+'>(.*?)</div>\s*</td>""", re.S)
_ORIGINALLY_ASSIGNED = re.compile(
    r"originally (?:assigned|classified) (?:to|as|under) (?:TC ?#?\s*)?"
    r"(\d\.[A-Za-z]\.\d+\.\d+\.\d+)")
PROTEINS_CAP = 300


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


def parse_protein_cells(html: str, family: str | None = None) -> dict[str, str]:
    """result.php?tc=<4-part> → ``{5-part id: protein cell text}``, lines joined by
    "; " ("UrtA-E urea porter of Anabaena sp. PCC7120; UrtA (R); UrtB (M); ...").
    The cell is the last one of the system's row (before the next row's link)."""
    start = _RESULT_CLUSTER_START.search(html)
    scope = html[start.start():] if start else html
    prefix = f"{family}." if family else None
    out: dict[str, str] = {}
    rows = list(_SYSTEM_ROW.finditer(scope))
    for i, m in enumerate(rows):
        tc_id = m.group(1)
        if prefix and not tc_id.startswith(prefix):
            continue
        end = rows[i + 1].start() if i + 1 < len(rows) else len(scope)
        row_rest = scope[m.end():end]
        row_rest = row_rest.split("</tr>", 1)[0] + "</tr>"
        cell = _PROTEIN_CELL.search(row_rest)
        if not cell:
            continue
        lines = [clean_html_fragment(x) for x in re.split(r"<br\s*/?>", cell.group(1))]
        text = "; ".join(x for x in lines if x)
        if text:
            out[tc_id] = cap_at_word_boundary(text, PROTEINS_CAP)
    return out


def parse_originally_assigned(html: str, family: str | None = None) -> dict[str, set[str]]:
    """``{old id: {row id}}`` from system rows whose text says "This system was
    originally assigned to TC# <old>" — the other wording TCDB uses for a move."""
    _, systems = parse_family_page(html, family=family)
    out: dict[str, set[str]] = {}
    for tc_id, text in systems.items():
        for old in _ORIGINALLY_ASSIGNED.findall(text):
            old = old.upper()
            if old != tc_id:
                out.setdefault(old, set()).add(tc_id)
    return out


# ── scope ────────────────────────────────────────────────────────────────────

def scoped_families(kept_ids: set[str],
                    hierarchy_ids: set[str] | None = None) -> dict[str, str]:
    """{family_3part: query_id} for every family carrying a kept 3/4/5-part
    ID. Any 4-part in the family works as the result.php query key — the page
    returns the entire family (spec 3b). A family kept only at 3-part level
    (no kept 4/5-part id) is queried by its smallest (numeric) 4-part child in
    `hierarchy_ids` (the full TCDB hierarchy keys); the live site returns an
    EMPTY page for a 3-part query. No such child (or no hierarchy) -> skipped."""
    out: dict[str, str] = {}
    family_only: set[str] = set()
    for tc in kept_ids:
        parts = tc.split(".")
        if len(parts) < 3:
            continue
        fam = ".".join(parts[:3])
        if len(parts) == 3:
            family_only.add(fam)
            continue
        four = ".".join(parts[:4])
        if fam not in out or four < out[fam]:
            out[fam] = four
    for fam in sorted(family_only - set(out)):
        kids = [t for t in (hierarchy_ids or ())
                if t.startswith(fam + ".") and t.count(".") == 3
                and t.rsplit(".", 1)[1].isdigit()]
        if not kids:
            log.warning(f"  family-only {fam}: no 4-part child in hierarchy; skipped")
            continue
        out[fam] = min(kids, key=lambda t: int(t.rsplit(".", 1)[1]))
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

def build_moved_map(pages_dir: Path, existing: dict | None = None) -> dict:
    """``{"meta": {...}, "moved": {old id: [targets]}}`` from every cached family page
    under ``pages_dir`` (no network). Merges with ``existing`` (target lists unioned), so a
    narrower scrape never drops a move. No date or page count is stored, so identical
    moves give a byte-identical file. Targets keep every conflicting move: step 6
    treats an old id with two targets as ambiguous."""
    from multiomics_kg.utils.tcdb_utils import parse_moved_tcids
    moved: dict[str, set[str]] = {
        old: set(t) for old, t in ((existing or {}).get("moved") or {}).items()}
    pages = sorted(Path(pages_dir).glob("family_*.html"))
    source_pages = set((existing or {}).get("meta", {}).get("source_pages") or [])
    for page in pages:
        text = page.read_text(encoding="utf-8", errors="replace")
        found = parse_moved_tcids([text])
        family = page.stem.removeprefix("family_")
        for old, targets in parse_originally_assigned(text, family=family).items():
            found.setdefault(old, set()).update(targets)
        if found:
            source_pages.add(page.name)
        for old, targets in found.items():
            moved.setdefault(old, set()).update(targets)
    return {
        "meta": {
            "source": "TCDB family pages (search/result.php), lines "
                      "<old> has been moved to <new>, and system rows saying "
                      "originally assigned to TC# <old>",
            "source_pages": sorted(source_pages),
        },
        "moved": {old: sorted(t) for old, t in sorted(moved.items())},
    }


def write_moved_map(pages_dir: Path = PAGES_DIR, out: Path = MOVED_FILE) -> dict:
    existing = json.loads(out.read_text()) if out.exists() else None
    art = build_moved_map(pages_dir, existing=existing)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(art, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    log.info(f"wrote {out} - {len(art['moved'])} moved old ids from "
             f"{len(art['meta']['source_pages'])} page(s)")
    return art

def collect_best_hit_tcids(genome_dirs: Iterable[Path]) -> set[str]:
    """Union of merged `tcdb_best_hit_ids` (best_hit_v1) over the given strain dirs."""
    out: set[str] = set()
    for d in genome_dirs:
        genes = load_gene_annotations(d) or {}
        for g in genes.values():
            out.update(t for t in g.get("tcdb_best_hit_ids") or [] if isinstance(t, str) and t)
    return out


def assemble_names(
    browse_names: dict[str, str],
    per_family_results: dict[str, tuple[dict[str, str], dict[str, str]]],
    kept_ids: set[str],
    scraped_families: list[str],
    failed_families: list[str],
    existing: dict | None = None,
    proteins: dict[str, str] | None = None,
    merge: bool = True,
) -> dict:
    """Build the artifact dict. Browse layer (levels 0-2) kept in FULL —
    removes gene-set coupling for those levels at ~150 KB. Deep layer (4/5
    part) trimmed to `kept_ids` (the wanted set: kept ∪ best-hit ids; spec §4
    scoped-scrape decision).

    With `existing` (the current artifact) the build is merge-preserving: every
    existing entry stays verbatim and only ids not yet named are added, so a
    re-scrape cannot rename an existing TcdbFamily; scraped/failed family lists
    are unions. Without it, a fresh build (``--reparse``: rebuild every entry from
    the cached pages, so parser fixes reach existing names).

    ``proteins`` (``{5-part id: protein cell}``) adds a ``proteins`` field to each
    system entry, and names a system whose name cell is empty after its protein
    cell's first line (e.g. "UrtA-E urea porter of Anabaena sp. PCC7120").
    ``meta.scraped_at`` only moves when the names change.
    """
    proteins = proteins or {}
    names: dict[str, dict] = dict((existing or {}).get("names", {})) if merge else {}
    for tc_id, name in browse_names.items():
        names.setdefault(tc_id, {"name": name})

    for _fam, (subfams, systems) in per_family_results.items():
        for tc_id, name in subfams.items():
            if tc_id in kept_ids:
                names.setdefault(tc_id, {"name": name})
        for tc_id, full_text in systems.items():
            if tc_id in kept_ids:
                names.setdefault(tc_id, {
                    "name": make_name(full_text),
                    "description": make_description(full_text),
                })
    for tc_id, cell in proteins.items():
        if tc_id not in kept_ids:
            continue
        entry = names.setdefault(tc_id, {"name": ""})
        entry.setdefault("proteins", cell)
        if not entry.get("name"):
            entry["name"] = make_name(cell.split("; ", 1)[0])

    prior = (existing or {}).get("meta", {})
    scraped = (set(prior.get("scraped_families", [])) if merge else set()) | set(scraped_families)
    # a family that has ever been scraped keeps valid names; "failed" = never scraped yet
    failed = ((set(prior.get("failed_families", [])) if merge else set())
              | set(failed_families)) - scraped
    names = dict(sorted(names.items()))
    unchanged = existing is not None and existing.get("names") == names
    return {
        "meta": {
            "scraped_at": prior.get("scraped_at") if unchanged and prior.get("scraped_at")
            else date.today().isoformat(),
            "source": "https://www.tcdb.org (browse.php + search/result.php per family)",
            "scraped_families": sorted(scraped),
            "failed_families": sorted(failed),
        },
        "names": names,
    }


# ── main ─────────────────────────────────────────────────────────────────────

def main(force: bool = False, delay: float = DEFAULT_DELAY_S,
         families: list[str] | None = None, reparse: bool = False) -> int:
    """``reparse``: no network — rebuild every entry from the cached pages (parser
    fixes reach existing names); an uncached page counts as failed."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    pruned = json.loads(PRUNED_FILE.read_text())
    kept_ids = set(pruned["kept_tcdb_ids"])
    # best_hit_v1: diamond best-hit 5-part systems are named too — since the TCDB
    # cleanup (2026-10-23) they are kept tc_system TcdbFamily nodes that the
    # diamond_best_hit Gene_has_tcdb_family edges target, and need a node name.
    best_hit_ids = collect_best_hit_tcids(PROJECT_ROOT / r["data_dir"] for r in load_genome_rows())
    wanted_ids = kept_ids | best_hit_ids
    hierarchy_ids = set(json.loads(HIERARCHY_FILE.read_text()))
    fams = scoped_families(wanted_ids, hierarchy_ids)
    if families:
        fams = {f: q for f, q in fams.items() if f in set(families)}
    log.info(f"scope: {len(fams)} families (from {len(kept_ids)} kept + "
             f"{len(best_hit_ids)} best-hit TC ids)")
    existing = json.loads(NAMES_FILE.read_text()) if NAMES_FILE.exists() else None

    def get(url: str, dest: Path) -> str:
        if reparse:
            if not dest.exists():
                raise FileNotFoundError(f"not cached: {dest.name}")
            return dest.read_text(encoding="utf-8", errors="replace")
        return fetch_page(url, dest, force=force, delay=delay)

    browse_html = get(BROWSE_URL, PAGES_DIR / "browse.html")
    browse_names = parse_browse(browse_html)
    log.info(f"browse.php: {len(browse_names)} class/subclass/family names")

    per_family: dict[str, tuple[dict[str, str], dict[str, str]]] = {}
    proteins: dict[str, str] = {}
    scraped: list[str] = []
    failed: list[str] = []
    consecutive = 0
    for i, (fam, query) in enumerate(sorted(fams.items()), 1):
        dest = PAGES_DIR / f"family_{fam}.html"
        try:
            page = get(RESULT_URL.format(tc=query), dest)
        except Exception as exc:  # noqa: BLE001 — per-family isolation
            log.warning(f"[{i}/{len(fams)}] {fam}: FAILED ({exc})")
            failed.append(fam)
            consecutive += 1
            if consecutive >= MAX_CONSECUTIVE_FAILURES and not reparse:
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
        proteins.update(parse_protein_cells(page, family=fam))
        scraped.append(fam)
        if i % 25 == 0 or i == len(fams):
            log.info(f"[{i}/{len(fams)}] fetched+parsed (last: {fam}, "
                     f"{len(subfams)} subfam names, {len(systems)} systems)")

    artifact = assemble_names(browse_names, per_family, wanted_ids, scraped, failed,
                              existing=existing, proteins=proteins, merge=not reparse)
    NAMES_FILE.parent.mkdir(parents=True, exist_ok=True)
    NAMES_FILE.write_text(json.dumps(artifact, indent=1, ensure_ascii=False) + "\n",
                          encoding="utf-8")
    n = len(artifact["names"])
    log.info(f"wrote {NAMES_FILE} — {n} names, "
             f"{len(scraped)} families scraped, {len(failed)} failed")
    write_moved_map()
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
    parser.add_argument("--reparse", action="store_true",
                        help="no network: rebuild every name from the cached pages "
                             "(parser fixes reach existing names)")
    parser.add_argument("--moved-only", action="store_true",
                        help="only rebuild tcdb_moved.json from the cached pages (no network)")
    args = parser.parse_args()
    if args.moved_only:
        logging.basicConfig(level=logging.INFO, format="%(message)s")
        write_moved_map()
    else:
        main(force=args.force, delay=args.delay, families=args.families,
             reparse=args.reparse)
