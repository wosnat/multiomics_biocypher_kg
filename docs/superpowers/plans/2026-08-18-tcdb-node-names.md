# TCDB Node Names (T6) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Name the 916 unnamed `TcdbFamily` nodes by scraping tcdb.org (browse.php + per-family result.php pages) into a committed `tcdb_names.json`, consumed by step 6's hierarchy build.

**Architecture:** A standalone, resumable, slow-cadence scraper (`multiomics_kg/download/scrape_tcdb_names.py`) fetches raw HTML into a gitignored page cache and writes a small committed names artifact. Step 6 (`build_kegg_metabolism_xrefs.py`) consumes the artifact when building `tcdb_hierarchy.json` (replacing the hardcoded `_TC_CLASS_NAMES`), and warns when kept TC IDs belong to families the scraper has not covered. The adapter gains a sparse `description` property on `tc_specificity` nodes; the `tcdbFamilyFullText` index gains `description`.

**Tech Stack:** Python 3 stdlib (`re`, `html`, `json`) + `requests` (already a dependency). No new dependencies.

**Spec:** `docs/superpowers/specs/2026-08-12-tcdb-node-names-design.md` (all design questions resolved; READY TO IMPLEMENT). Read it before starting — especially sections 3a/3b (verified live markup), 4 (scoped-scrape decision), 5b/5c/5d (name shape, fallback, hygiene).

## Global Constraints

- **Sanitise every string property before yielding from an adapter**: `.replace("'", "^").replace("|", "")`. Applies to `name` and `description`. HTML entities and tags must be unwrapped BEFORE sanitising. (Spec §7; CLAUDE.md.)
- **No native `bool` properties** (vocabulary contract R5). This plan adds none.
- **No new categorical property** → no `config/controlled_vocabularies.yaml` change. `name`/`description` are free text. The fallback marker is the predicate `t.name = t.tcdb_id` itself (spec 5c).
- `cache/data/tcdb/tcdb_hierarchy.json`, `tcdb_pruned.json`, and the new `cache/data/tcdb/tcdb_names.json` are **committed**; `cache/data/tcdb/raw/` (including the new `raw/pages/`) is **gitignored** (already covered by the `cache/data/tcdb/raw` pattern).
- Scrape cadence (spec 5d): single-threaded, **2.5 s between requests**, polite User-Agent with contact address, exponential backoff, abort after 3 consecutive family failures, resume = skip already-cached pages.
- Name shape (spec 5b): `name` = citation-strip → first sentence → 150-char word-boundary cap. `description` = cleaned full text (citations kept), 400-char word-boundary cap, `tc_specificity` only.
- Iteration loop: `bash scripts/prepare_data.sh --steps 6 --force` (~1 min, no network).
- `scripts/post-import.sh` and `scripts/post-import.cypher` must keep identical Cypher logic.
- Unit test suite: `pytest -m "not slow and not kg"`.

## File Structure

| File | Role |
|---|---|
| `multiomics_kg/download/scrape_tcdb_names.py` (create) | Fetchers + parsers + text cleaners + CLI. Importable; all pure functions unit-testable offline. |
| `tests/test_scrape_tcdb_names.py` (create) | Unit tests for cleaners, parsers, scoping, artifact assembly. |
| `multiomics_kg/download/build_kegg_metabolism_xrefs.py` (modify) | Consume names artifact in `build_tcdb_hierarchy`; delete `_TC_CLASS_NAMES`; unscraped-family guard; fix stale module docstring. |
| `tests/test_build_kegg_metabolism_xrefs.py` (modify) | Hierarchy-consumes-names tests. |
| `multiomics_kg/adapters/tcdb_adapter.py` (modify) | Emit sparse `description` on nodes. |
| `tests/test_tcdb_adapter.py` (modify) | Fixture + assertions for `description`; class-name fixture update. |
| `config/schema_config.yaml` (modify) | `description: str` on `tcdb family`. |
| `scripts/post-import.sh` + `scripts/post-import.cypher` (modify) | `tcdbFamilyFullText` drop + recreate with `description`. |
| `cache/data/tcdb/tcdb_names.json` (generated, committed) | The names artifact. |
| `.claude/skills/tcdb-diamond/SKILL.md` (modify) | "Refreshing TCDB names" section. |
| `CHANGELOG.md`, `CLAUDE.md`, `docs/kg-changes/tcdb-two-source-upgrade.md` (modify) | Documentation. |

### Names artifact format (`cache/data/tcdb/tcdb_names.json`)

```json
{
  "meta": {
    "scraped_at": "2026-08-18",
    "source": "https://www.tcdb.org (browse.php + search/result.php per family)",
    "scraped_families": ["1.A.11", "2.A.1"],
    "failed_families": []
  },
  "names": {
    "1": {"name": "Channels/Pores"},
    "1.A": {"name": "α-Type Channels"},
    "2.A.1": {"name": "The Major Facilitator Superfamily (MFS)"},
    "2.A.1.1": {"name": "The Sugar Porter (SP) Family"},
    "2.A.1.1.1": {
      "name": "Galactose:H+ symporter, GalP.",
      "description": "Galactose:H+ symporter, GalP. Also transports glucose, xylose, ... (Henderson and Giddens 1977; ...)."
    }
  }
}
```

- The **browse layer is kept in full** (all 8 classes + all ~56 subclasses + all ~2,127 families): ~150 KB extra removes gene-set coupling entirely for levels 0–2. Only the deep layer (4/5-part) is scoped to kept IDs. Total artifact ≈ 300 KB. (This deliberately exceeds the spec's ~100 KB estimate — decision recorded here: full browse layer costs nothing extra in requests and one fewer coupling.)
- `meta.scraped_families` is what the step-6 guard checks kept families against.

---

### Task 1: Text-cleaning helpers

**Files:**
- Create: `multiomics_kg/download/scrape_tcdb_names.py`
- Test: `tests/test_scrape_tcdb_names.py`

**Interfaces:**
- Produces (used by Tasks 2–3):
  - `clean_html_fragment(s: str) -> str` — tag-strip (tolerates dangling `<br /` at end), entity-unescape, whitespace-collapse.
  - `strip_citations(s: str) -> str` — removes `(... 19xx/20xx ...)` parentheticals, iteratively (handles nesting/multi-ref groups).
  - `first_sentence(s: str) -> str` — abbreviation-safe (`E. coli`, `et al.`, `e.g.`, single-capital genus initials).
  - `cap_at_word_boundary(s: str, limit: int) -> str` — ≤limit chars, cut at word boundary, `…` suffix, trailing punctuation stripped before the ellipsis.
  - `make_name(full_text: str) -> str` — `cap_at_word_boundary(first_sentence(strip_citations(full_text)), 150)`.
  - `make_description(full_text: str) -> str` — `cap_at_word_boundary(full_text, 400)`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_scrape_tcdb_names.py
"""Unit tests for the TCDB names scraper (T6)."""
from multiomics_kg.download import scrape_tcdb_names as m


class TestCleanHtmlFragment:
    def test_strips_tags_and_entities(self):
        s = "<p>Drug:H<sup>+</sup> Antiporter &alpha; &nbsp;porin</p>"
        assert m.clean_html_fragment(s) == "Drug:H + Antiporter α porin"

    def test_tolerates_dangling_unclosed_tag_at_end(self):
        # observed live: one cell ends with an unclosed "<br /"
        s = "Glycerol uptake permease Stl1.<br /"
        assert m.clean_html_fragment(s) == "Glycerol uptake permease Stl1."

    def test_collapses_whitespace(self):
        assert m.clean_html_fragment("a\n  b\t c") == "a b c"


class TestStripCitations:
    def test_removes_year_parenthetical(self):
        s = "Galactose:H+ symporter (Henderson and Giddens 1977). More text."
        assert m.strip_citations(s) == "Galactose:H+ symporter. More text."

    def test_keeps_non_citation_parenthetical(self):
        s = "The Sugar Porter (SP) Family"
        assert m.strip_citations(s) == "The Sugar Porter (SP) Family"

    def test_removes_multi_ref_group(self):
        s = "GalP (Henderson 1977; Patching et al., 2008). End."
        assert m.strip_citations(s) == "GalP. End."

    def test_handles_nested_parentheses(self):
        s = "Stl1 (similar to Stl1 of S. cerevisiae (Ferreira et al., 2005))."
        assert m.strip_citations(s) == "Stl1."


class TestFirstSentence:
    def test_plain_split(self):
        s = "MDR efflux pump, MdeA. Exports quaternary ammonium compounds."
        assert m.first_sentence(s) == "MDR efflux pump, MdeA."

    def test_does_not_split_on_genus_initial(self):
        s = "Transporter of E. coli origin. Second sentence."
        assert m.first_sentence(s) == "Transporter of E. coli origin."

    def test_does_not_split_on_et_al(self):
        s = "Characterized by Smith et al. Second clause continues here."
        assert m.first_sentence(s) == "Characterized by Smith et al. Second clause continues here."

    def test_no_terminal_period(self):
        s = "The glucose transport protein, GTP1"
        assert m.first_sentence(s) == "The glucose transport protein, GTP1"


class TestCaps:
    def test_short_string_unchanged(self):
        assert m.cap_at_word_boundary("short name", 150) == "short name"

    def test_long_string_capped_at_word_boundary(self):
        s = "word " * 60  # 300 chars
        out = m.cap_at_word_boundary(s.strip(), 150)
        assert len(out) <= 151  # 150 + ellipsis char
        assert out.endswith("…")
        assert not out[:-1].endswith(" ")

    def test_make_name_pipeline(self):
        s = ("Glycerol-P:Pi antiporter (Ambudkar et al. 1986). "
             "The 3-d structure is known.")
        assert m.make_name(s) == "Glycerol-P:Pi antiporter."
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_scrape_tcdb_names.py -v`
Expected: FAIL with `ModuleNotFoundError` / `ImportError` (module does not exist yet).

- [ ] **Step 3: Write the module skeleton with the cleaners**

```python
# multiomics_kg/download/scrape_tcdb_names.py
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
cache/data/tcdb/tcdb_pruned.json (~351 families), plus the full browse layer.

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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_scrape_tcdb_names.py -v`
Expected: all PASS. If `test_strips_tags_and_entities` disagrees on internal spacing (tag-strip inserts a space where `<sup>` sat), adjust the EXPECTED value to what the implementation produces **only if** the output is still readable ("Drug:H + Antiporter" is acceptable; "Drug:H+Antiporter" would also be acceptable via `re.sub(r"<[^>]+>", "", s)` — pick one convention and keep test + code consistent).

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/scrape_tcdb_names.py tests/test_scrape_tcdb_names.py
git commit -m "feat(tcdb-names): text-cleaning helpers for the T6 name scraper"
```

---

### Task 2: HTML parsers (browse.php + family result pages)

**Files:**
- Modify: `multiomics_kg/download/scrape_tcdb_names.py`
- Test: `tests/test_scrape_tcdb_names.py`

**Interfaces:**
- Consumes: `clean_html_fragment` (Task 1).
- Produces (used by Task 3):
  - `parse_browse(html: str) -> dict[str, str]` — `{tc_id: name}` for every 1/2/3-part entry on browse.php; names have the `"X:"` prefix stripped; zero-name entries omitted.
  - `parse_family_page(html: str) -> tuple[dict[str, str], dict[str, str]]` — `(subfamily_names, system_texts)`: 4-part → cleaned header name; 5-part → cleaned FULL text of the Name cell (uncapped — capping happens at artifact-assembly time). Both scoped to the `result-cluster` table (the page embeds off-family 5-part IDs in a sidebar module).

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_scrape_tcdb_names.py`:

```python
BROWSE_FIXTURE = """
<ul id="red" class="treeview-red">
  <li><span>
    <div rel="1" class="entry">
      <div class="tcid name">&nbsp;1:&nbsp;Channels/Pores</div>
    </div>
  </span>
  <ul><li><span>
    <div rel="1.A" class="entry">
      <div class="tcid name">&nbsp;1.A:&nbsp;&#945;-Type Channels</div>
    </div>
  </span>
  <ul><li><span>
    <div rel="1.A.1" class="entry" style="border-bottom:1px dotted #999;">
      <div class="tcid name" style="cursor:pointer;">&nbsp;1.A.1:&nbsp;The Voltage-gated Ion Channel (VIC) Superfamily </div>
    </div>
  </span></li></ul>
  </li></ul></li>
</ul>
"""

FAMILY_PAGE_FIXTURE = """
<A HREF="/search/result.php?tc=2.A.1">Read Family Description</A>
<table id="result-cluster" style="vertical-align:top">
  <tr><th>TCID</th><th>Name</th><th>Domain</th><th>Kingdom/Phylum</th><th class="right-border">Protein(s)</th></tr>
  <tr><td colspan="4" id="right-border"><strong><A id="2.A.1.1"></A>2.A.1.1:&nbsp;&nbsp;The Sugar Porter (SP) Family</strong></td></tr>
  <tr><td width="60" style="vertical-align:top">
    <A id="2.A.1.1.1"></A>
    <A valign="top" HREF="/search/result.php?tc=2.A.1.1.1">
2.A.1.1.1</A><br></td><td><div class='400aa753232406c9e9e7370875b4f60b'><p>Galactose:H<sup>+</sup> symporter, GalP (<a class="reflink" href="/search/result.php?tc=2.A.1#ref7167">Hern&aacute;ndez-Montalvo <em>et al.</em>, 2001</a>). Long tail text.</p></div></td>
    <td>Bacteria</td><td>Pseudomonadota</td>
    <td id="right-border" width="170"><div class='x'>GalP of <em>E. coli</em></div></td></tr>
  <tr><td colspan="4" id="right-border"><strong><A id="2.A.1.2"></A>2.A.1.2:&nbsp;&nbsp;The Drug:H<sup>+</sup> Antiporter-1 (12 Spanner) (DHA1) Family</strong></td></tr>
</table>
<div id="sidebar"><a HREF="/search/result.php?tc=9.Z.9.9.9">9.Z.9.9.9</a></div>
"""

UNNAMED_FAMILY_PAGE_FIXTURE = """
<table id="result-cluster" style="vertical-align:top">
  <tr><th>TCID</th><th>Name</th><th>Domain</th><th>Kingdom/Phylum</th><th class="right-border">Protein(s)</th></tr>
  <tr><td width="60" style="vertical-align:top">
    <A id="1.A.11.1"></A>
    <A id="1.A.11.1.1"></A>
    <A valign="top" HREF="/search/result.php?tc=1.A.11.1.1">
1.A.11.1.1</A><br></td><td><div class='cf8b031f73d53100bb24514ad08c11f0'><p>Ammonia transporter and regulatory sensor, AmtB.</p></div></td>
    <td>Bacteria</td><td>Pseudomonadota</td>
    <td id="right-border"><div class='y'>AmtB of <em>E. coli</em></div></td></tr>
</table>
"""


class TestParseBrowse:
    def test_parses_all_three_levels(self):
        names = m.parse_browse(BROWSE_FIXTURE)
        assert names["1"] == "Channels/Pores"
        assert names["1.A"] == "α-Type Channels"
        assert names["1.A.1"] == "The Voltage-gated Ion Channel (VIC) Superfamily"

    def test_no_4part_entries(self):
        assert all(k.count(".") <= 2 for k in m.parse_browse(BROWSE_FIXTURE))


class TestParseFamilyPage:
    def test_named_subfamily_headers(self):
        subfams, systems = m.parse_family_page(FAMILY_PAGE_FIXTURE)
        assert subfams["2.A.1.1"] == "The Sugar Porter (SP) Family"
        # inline markup inside the name must survive the capture (spec 6a):
        assert subfams["2.A.1.2"] == "The Drug:H + Antiporter-1 (12 Spanner) (DHA1) Family"

    def test_system_full_text(self):
        _, systems = m.parse_family_page(FAMILY_PAGE_FIXTURE)
        assert systems["2.A.1.1.1"].startswith("Galactose:H + symporter, GalP")
        assert "Long tail text." in systems["2.A.1.1.1"]
        assert "reflink" not in systems["2.A.1.1.1"]

    def test_sidebar_ids_excluded(self):
        _, systems = m.parse_family_page(FAMILY_PAGE_FIXTURE)
        assert "9.Z.9.9.9" not in systems

    def test_unnamed_family_has_no_headers_but_systems_parse(self):
        subfams, systems = m.parse_family_page(UNNAMED_FAMILY_PAGE_FIXTURE)
        assert subfams == {}
        assert systems["1.A.11.1.1"] == "Ammonia transporter and regulatory sensor, AmtB."
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_scrape_tcdb_names.py -k "ParseBrowse or ParseFamilyPage" -v`
Expected: FAIL with `AttributeError: ... has no attribute 'parse_browse'`.

- [ ] **Step 3: Implement the parsers**

Append to `scrape_tcdb_names.py` (after the cleaners):

```python
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
_RESULT_CLUSTER = re.compile(r'<table id="result-cluster".*?</table>', re.S)


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


def parse_family_page(html: str) -> tuple[dict[str, str], dict[str, str]]:
    """result.php?tc=<4-part> → (subfamily_names, system_full_texts).

    Both scoped to the result-cluster table. A subfamily with no <strong>
    header row has NO upstream name (spec 6d) — absence is the signal.
    """
    match = _RESULT_CLUSTER.search(html)
    if not match:
        return {}, {}
    tbl = match.group(0)

    subfams: dict[str, str] = {}
    for tc_id, raw in _SUBFAMILY_HEADER.findall(tbl):
        name = clean_html_fragment(raw)
        name = re.sub(rf"^{re.escape(tc_id)}\s*:\s*", "", name).strip()
        if name:
            subfams[tc_id] = name

    systems: dict[str, str] = {}
    for tc_id, raw in _SYSTEM_ROW.findall(tbl):
        text = clean_html_fragment(raw)
        if text:
            systems[tc_id] = text
    return subfams, systems
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_scrape_tcdb_names.py -v`
Expected: all PASS. The `"Drug:H + Antiporter-1"` / `"Galactose:H + symporter"` spacing must match whatever convention Task 1 locked in — keep fixture expectations consistent with `clean_html_fragment`.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/scrape_tcdb_names.py tests/test_scrape_tcdb_names.py
git commit -m "feat(tcdb-names): browse.php + family-page parsers"
```

---

### Task 3: Scope, fetch loop, artifact assembly, CLI

**Files:**
- Modify: `multiomics_kg/download/scrape_tcdb_names.py`
- Test: `tests/test_scrape_tcdb_names.py`

**Interfaces:**
- Consumes: `parse_browse`, `parse_family_page`, `make_name`, `make_description` (Tasks 1–2).
- Produces:
  - `scoped_families(kept_ids: set[str]) -> dict[str, str]` — `{family_3part: query_4part}`; query id is the lexicographically smallest kept 4-part prefix per family (any 4-part works as query key — the page covers the whole family).
  - `assemble_names(browse_names, per_family_results, kept_ids, scraped_families, failed_families) -> dict` — the full artifact dict (`meta` + `names`), deep layer trimmed to kept IDs.
  - `main(force: bool, delay: float, families: list[str] | None) -> int` — orchestrates; returns count of names written.
  - Page cache layout: `PAGES_DIR / "browse.html"`, `PAGES_DIR / f"family_{fam}.html"` (e.g. `family_2.A.1.html`).

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_scrape_tcdb_names.py`:

```python
class TestScopedFamilies:
    def test_derives_families_from_kept_deep_ids(self):
        kept = {"1", "1.A", "2.A.1", "2.A.1.3", "2.A.1.1.1", "3.A.1.5"}
        fams = m.scoped_families(kept)
        # 1-, 2-, 3-part kept ids do NOT create scrape work
        assert fams == {"2.A.1": "2.A.1.1", "3.A.1": "3.A.1.5"}

    def test_query_id_is_smallest_kept_4part(self):
        kept = {"2.A.1.7.2", "2.A.1.3"}
        assert m.scoped_families(kept) == {"2.A.1": "2.A.1.3"}


class TestAssembleNames:
    def test_deep_layer_trimmed_to_kept_ids(self):
        browse = {"1": "Channels/Pores", "2.A.1": "MFS"}
        per_family = {
            "2.A.1": (
                {"2.A.1.1": "The Sugar Porter (SP) Family",
                 "2.A.1.2": "DHA1"},
                {"2.A.1.1.1": "Galactose:H+ symporter, GalP (Henderson 1977). Tail.",
                 "2.A.1.1.2": "Not kept system."},
            )
        }
        kept = {"1", "2.A.1", "2.A.1.1", "2.A.1.1.1"}
        art = m.assemble_names(browse, per_family, kept,
                               scraped_families=["2.A.1"], failed_families=[])
        names = art["names"]
        assert names["1"] == {"name": "Channels/Pores"}          # full browse layer
        assert names["2.A.1.1"] == {"name": "The Sugar Porter (SP) Family"}
        assert "2.A.1.2" not in names                            # not kept
        assert "2.A.1.1.2" not in names                          # not kept
        assert names["2.A.1.1.1"]["name"] == "Galactose:H+ symporter, GalP."
        assert names["2.A.1.1.1"]["description"].startswith(
            "Galactose:H+ symporter, GalP (Henderson 1977). Tail.")
        assert art["meta"]["scraped_families"] == ["2.A.1"]

    def test_browse_layer_kept_even_when_not_in_kept_ids(self):
        art = m.assemble_names({"6": "Membrane Transporter Metabolons (MTM)"},
                               {}, set(), scraped_families=[], failed_families=[])
        assert art["names"]["6"] == {"name": "Membrane Transporter Metabolons (MTM)"}


class TestFetchResume:
    def test_skips_existing_page(self, tmp_path, monkeypatch):
        dest = tmp_path / "family_2.A.1.html"
        dest.write_text("cached")
        called = []
        monkeypatch.setattr(m, "_http_get", lambda url: called.append(url) or "fresh")
        out = m.fetch_page("http://x", dest, force=False, delay=0)
        assert out == "cached" and called == []

    def test_force_refetches(self, tmp_path, monkeypatch):
        dest = tmp_path / "family_2.A.1.html"
        dest.write_text("cached")
        monkeypatch.setattr(m, "_http_get", lambda url: "fresh")
        monkeypatch.setattr(m.time, "sleep", lambda s: None)
        out = m.fetch_page("http://x", dest, force=True, delay=0)
        assert out == "fresh" and dest.read_text() == "fresh"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_scrape_tcdb_names.py -k "ScopedFamilies or AssembleNames or FetchResume" -v`
Expected: FAIL with `AttributeError`.

- [ ] **Step 3: Implement scope + fetch + assembly + CLI**

Append to `scrape_tcdb_names.py`:

```python
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
    """Fetch url → dest. Resume = return cached file when present (unless
    force). Sleeps `delay` BEFORE every network request; exponential backoff
    between retries. Raises after `retries` failures."""
    if dest.exists() and not force:
        return dest.read_text(encoding="utf-8", errors="replace")
    last_exc: Exception | None = None
    for attempt in range(retries):
        time.sleep(delay * (2 ** attempt) if attempt else delay)
        try:
            text = _http_get(url)
        except Exception as exc:              # noqa: BLE001 — retried, then re-raised
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
        except Exception as exc:              # noqa: BLE001 — per-family isolation
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
        subfams, systems = parse_family_page(page)
        per_family[fam] = (subfams, systems)
        scraped.append(fam)
        if i % 25 == 0 or i == len(fams):
            log.info(f"[{i}/{len(fams)}] fetched+parsed (last: {fam}, "
                     f"{len(subfams)} subfam names, {len(systems)} systems)")

    artifact = assemble_names(browse_names, per_family, kept_ids, scraped, failed)
    NAMES_FILE.write_text(json.dumps(artifact, indent=1, sort_keys=False,
                                     ensure_ascii=False) + "\n", encoding="utf-8")
    n = len(artifact["names"])
    log.info(f"wrote {NAMES_FILE} — {n} names, "
             f"{len(scraped)} families scraped, {len(failed)} failed")
    return n


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--force", action="store_true",
                        help="re-fetch pages even when cached")
    parser.add_argument("--delay", type=float, default=DEFAULT_DELAY_S,
                        help="seconds between requests (default 2.5)")
    parser.add_argument("--families", nargs="*", default=None,
                        help="restrict to these 3-part family ids (top-up mode)")
    args = parser.parse_args()
    main(force=args.force, delay=args.delay, families=args.families)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_scrape_tcdb_names.py -v`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/scrape_tcdb_names.py tests/test_scrape_tcdb_names.py
git commit -m "feat(tcdb-names): resumable scoped fetch loop + artifact assembly + CLI"
```

---

### Task 4: Step-6 consumes the names artifact

**Files:**
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py`
- Test: `tests/test_build_kegg_metabolism_xrefs.py`

**Interfaces:**
- Consumes: the artifact file format (`{"meta": {...}, "names": {tc_id: {"name", ["description"]}}}`).
- Produces:
  - `build_tcdb_hierarchy(..., names_path: Path | None = None)` — new keyword-only-style trailing param; when the file exists, its names fill every level and `description` lands on `tc_specificity` entries.
  - `_TC_CLASS_NAMES` **deleted**.
  - `_warn_unscraped_families(kept_ids: set[str], names_path: Path) -> list[str]` — returns (and logs) sorted kept families absent from `meta.scraped_families`.

- [ ] **Step 1: Write the failing tests**

Look at the existing fixtures in `tests/test_build_kegg_metabolism_xrefs.py` around line 233 (`mod._build_tcdb_hierarchy(cache_root=cache_root)`) — reuse the same fixture-writing helpers the neighboring tests use for the 4 TSVs. Add:

```python
def _write_names_file(cache_root, names, scraped_families):
    import json
    tcdb_dir = cache_root / "tcdb"
    tcdb_dir.mkdir(parents=True, exist_ok=True)
    (tcdb_dir / "tcdb_names.json").write_text(json.dumps({
        "meta": {"scraped_at": "2026-08-18", "source": "test",
                 "scraped_families": scraped_families, "failed_families": []},
        "names": names,
    }))


def test_build_tcdb_hierarchy_consumes_names_artifact(tmp_path):
    # reuse the same raw-TSV fixture setup as
    # test_build_tcdb_hierarchy_seeds_5part_from_acc2tcid (copy its setup
    # lines verbatim), then:
    _write_names_file(cache_root, {
        "1": {"name": "Channels/Pores"},
        "1.A": {"name": "α-Type Channels"},
        "1.A.1.5": {"name": "The Epsilon Subfamily"},
        "1.A.1.5.2": {"name": "Short name.",
                      "description": "Short name. Much longer curated text."},
    }, scraped_families=["1.A.1"])
    mod._build_tcdb_hierarchy(cache_root=cache_root)
    h = json.loads((cache_root / "tcdb" / "tcdb_hierarchy.json").read_text())
    assert h["1"]["name"] == "Channels/Pores"           # replaces _TC_CLASS_NAMES
    assert h["1.A"]["name"] == "α-Type Channels"        # was always ""
    assert h["1.A.1.5"]["name"] == "The Epsilon Subfamily"
    assert h["1.A.1.5.2"]["name"] == "Short name."
    assert h["1.A.1.5.2"]["description"] == "Short name. Much longer curated text."


def test_build_tcdb_hierarchy_without_names_file_leaves_blank(tmp_path):
    # same TSV fixture setup, NO names file
    mod._build_tcdb_hierarchy(cache_root=cache_root)
    h = json.loads((cache_root / "tcdb" / "tcdb_hierarchy.json").read_text())
    assert h["1.A"]["name"] == ""
    assert "description" not in h["1.A.1.5.2"]
    # class names come ONLY from the artifact now — no hardcoded fallback
    assert h["1"]["name"] == ""


def test_warn_unscraped_families(caplog):
    import logging
    _write_names_file(cache_root, {}, scraped_families=["2.A.1"])
    with caplog.at_level(logging.WARNING):
        missing = mod._warn_unscraped_families(
            {"2.A.1.1", "3.A.1.5.9", "1.A"},  # 1.A too shallow to need scraping
            cache_root / "tcdb" / "tcdb_names.json")
    assert missing == ["3.A.1"]
    assert "3.A.1" in caplog.text and "scrape_tcdb_names" in caplog.text
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_build_kegg_metabolism_xrefs.py -k "names_artifact or without_names or unscraped" -v`
Expected: FAIL (`_TC_CLASS_NAMES` still active → class name assertion fails; `_warn_unscraped_families` missing).

- [ ] **Step 3: Implement**

In `build_kegg_metabolism_xrefs.py`:

1. **Delete** the `_TC_CLASS_NAMES` dict (lines ~89–96).
2. Add next to `_parse_tcdb_families`:

```python
def _load_tcdb_names(names_path: Path | None) -> tuple[dict[str, dict], list[str]]:
    """Read the committed tcdb_names.json (built by scrape_tcdb_names).

    Returns ({tc_id: {"name": str[, "description": str]}}, scraped_families).
    Missing file → empty maps (deep levels stay unnamed; the caller warns).
    """
    if names_path is None or not names_path.exists():
        return {}, []
    data = json.loads(names_path.read_text())
    return data.get("names", {}), data.get("meta", {}).get("scraped_families", [])
```

3. `build_tcdb_hierarchy(...)` gains a trailing parameter `names_path: Path | None = None`. At the top of the function body add `names, _scraped = _load_tcdb_names(names_path)`, then change each `ensure_*`:
   - `ensure_class`: `"name": names.get(cls, {}).get("name", "")`
   - `ensure_subclass`: `"name": names.get(subcls, {}).get("name", "")`
   - `ensure_family`: `"name": fam_descs.get(fam, "") or names.get(fam, {}).get("name", "")` (families.tsv stays primary; browse fills the 3 blanks)
   - `ensure_subfamily`: `"name": names.get(subfam, {}).get("name", "")`
   - `ensure_specificity`: `"name": names.get(tcid, {}).get("name", "")`, and after building `node`, add:
     ```python
     desc = names.get(tcid, {}).get("description", "")
     if desc:
         node["description"] = desc
     ```
4. `_build_tcdb_hierarchy` passes `names_path=tcdb_dir / "tcdb_names.json"`.
5. Add the guard function:

```python
def _warn_unscraped_families(kept_ids: set[str], names_path: Path) -> list[str]:
    """Spec §4 self-healing guard: kept 4/5-part TC IDs whose family the name
    scraper has not covered — loud staleness after a strain onboarding, never
    silent. Returns the sorted missing family list."""
    needed = {".".join(tc.split(".")[:3]) for tc in kept_ids
              if len(tc.split(".")) >= 4}
    if not names_path.exists():
        if needed:
            log.warning(
                f"  cache/data/tcdb/tcdb_names.json is MISSING — {len(needed)} "
                f"families' subfamily/specificity nodes will render as bare TC "
                f"ids. Run: uv run python -m multiomics_kg.download.scrape_tcdb_names"
            )
        return sorted(needed)
    _names, scraped = _load_tcdb_names(names_path)
    missing = sorted(needed - set(scraped))
    if missing:
        log.warning(
            f"  {len(missing)} kept TCDB families not covered by tcdb_names.json "
            f"(e.g. {missing[:5]}) — their deep nodes render as bare TC ids. "
            f"Top-up: uv run python -m multiomics_kg.download.scrape_tcdb_names "
            f"--families {' '.join(missing[:5])} ..."
        )
    return missing
```

6. In `main()`, right after the kept set is computed by `_prune_tcdb` (grep for the `_prune_tcdb(` call site), add:
   ```python
   _warn_unscraped_families(kept, cache_root / "tcdb" / "tcdb_names.json")
   ```
   (`kept` = whatever local name the call site binds the kept-ID set to.)
7. **Fix the stale module docstring** (~line 23): replace *"Step 6 also downloads the 3 TCDB reference TSVs (tc_classes, tc_subclasses, families) and writes the assembled hierarchy…"* with *"Step 6 also downloads the 4 TCDB reference TSVs (families, superfamilies, substrates, acc2tcid) and writes the assembled hierarchy… Node names come from the committed cache/data/tcdb/tcdb_names.json (built separately by `python -m multiomics_kg.download.scrape_tcdb_names` — see that module's docstring); step 6 itself never scrapes."*

- [ ] **Step 4: Run the tests**

Run: `pytest tests/test_build_kegg_metabolism_xrefs.py -v && pytest tests/test_tcdb_adapter.py -v`
Expected: the new tests PASS. `tests/test_tcdb_adapter.py::test_orchestrator_node_props` may FAIL if its fixture hierarchy relied on `_TC_CLASS_NAMES` — check how that test builds its fixture `tcdb_hierarchy.json`: if the fixture embeds `"name": "Channels and Pores"` directly it still passes (leave it); if it calls `build_tcdb_hierarchy` it now needs a names fixture. Fix forward in Task 5 if the failure is about `description`.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/build_kegg_metabolism_xrefs.py tests/test_build_kegg_metabolism_xrefs.py
git commit -m "feat(tcdb-names): step 6 consumes tcdb_names.json; delete _TC_CLASS_NAMES; unscraped-family guard"
```

---

### Task 5: Adapter `description` + schema + full-text index

**Files:**
- Modify: `multiomics_kg/adapters/tcdb_adapter.py` (in `get_nodes`, ~line 310)
- Modify: `config/schema_config.yaml` (`tcdb family` block, ~line 1260)
- Modify: `scripts/post-import.sh` (~line 92) and `scripts/post-import.cypher` (~line 73)
- Test: `tests/test_tcdb_adapter.py`

**Interfaces:**
- Consumes: hierarchy entries now carrying sparse `description` (Task 4).
- Produces: `TcdbFamily.description` (sparse, sanitized) node property; `tcdbFamilyFullText` covering it.

- [ ] **Step 1: Write the failing test**

In `tests/test_tcdb_adapter.py`, find the fixture that writes the test `tcdb_hierarchy.json` (used by `_make_orchestrator`) and add `"description": "Epsilon leaf curated text with 'quote' and |pipe|."` to the `1.A.1.5.2` entry. Then extend `test_orchestrator_node_props`:

```python
    # description: sparse, sanitized (quote → ^, pipe removed)
    nid, _label, props = nodes["tcdb:1.A.1.5.2"]
    assert props["description"] == "Epsilon leaf curated text with ^quote^ and pipe."
    # nodes without a description omit the key entirely
    nid, _label, props = nodes["tcdb:1.A.1.5"]
    assert "description" not in props
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_tcdb_adapter.py::test_orchestrator_node_props -v`
Expected: FAIL — `description` not emitted.

- [ ] **Step 3: Implement**

1. `tcdb_adapter.py` `get_nodes`, after the `superfamily` block:
```python
            if entry.get("description"):
                props["description"] = _clean_str(entry["description"])
```
2. `config/schema_config.yaml`, `tcdb family` properties — add after `superfamily`:
```yaml
    description: str        # sparse, tc_specificity only — scraped TCDB curated text, 400-char capped (T6)
```
   Also update the `name` comment to: `# TCDB entry name (scraped for subclass/subfamily/specificity since T6); falls back to tcdb_id when upstream has no name`
3. `scripts/post-import.sh` AND `scripts/post-import.cypher` — replace the `tcdbFamilyFullText` CREATE with the drop+recreate pattern already used by `interproEntryFullText` (full-text defs can't be ALTERed):
```cypher
// Full-text defs can't be ALTERed — drop + recreate so description is picked up
// even on reruns against an existing graph (T6 node names).
DROP INDEX tcdbFamilyFullText IF EXISTS;
CREATE FULLTEXT INDEX tcdbFamilyFullText IF NOT EXISTS
    FOR (t:TcdbFamily) ON EACH [t.name, t.tcdb_id, t.superfamily, t.description];
```
   Keep both files' Cypher identical.

- [ ] **Step 4: Run tests**

Run: `pytest tests/test_tcdb_adapter.py -v && pytest -m "not slow and not kg" -q`
Expected: all PASS (full unit suite guards against collateral damage).

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/tcdb_adapter.py config/schema_config.yaml scripts/post-import.sh scripts/post-import.cypher tests/test_tcdb_adapter.py
git commit -m "feat(tcdb-names): TcdbFamily.description property + fulltext index coverage"
```

---

### Task 6: Run the scrape (live, slow) and commit the artifact

**Files:**
- Generated: `cache/data/tcdb/tcdb_names.json` (committed), `cache/data/tcdb/raw/pages/*.html` (gitignored)

- [ ] **Step 1: Sanity-check scope offline**

```bash
uv run python -c "
import json
from multiomics_kg.download.scrape_tcdb_names import scoped_families, PRUNED_FILE
kept = set(json.loads(PRUNED_FILE.read_text())['kept_tcdb_ids'])
fams = scoped_families(kept)
print(len(fams), 'families to scrape')  # expect ~351
"
```

- [ ] **Step 2: Run the scraper in the background** (~20 min at 2.5 s/request)

```bash
uv run python -m multiomics_kg.download.scrape_tcdb_names 2>&1 | tee logs/scrape_tcdb_names.log
```
Run in background; monitor the log. If it aborts on consecutive failures (tcdb.org flaky), re-run later — cached pages resume.

- [ ] **Step 3: Inspect the artifact**

```bash
uv run python - <<'EOF'
import json
art = json.load(open('cache/data/tcdb/tcdb_names.json'))
names = art['names']
by_depth = {}
for k in names:
    by_depth.setdefault(k.count('.'), []).append(k)
for d in sorted(by_depth):
    print(f"{d+1}-part: {len(by_depth[d])}")
print("failed families:", art['meta']['failed_families'])
print("with description:", sum(1 for v in names.values() if 'description' in v))
# spot checks
for k in ('1', '1.A', '2.A.1', '3.F'):
    print(k, '->', names.get(k, {}).get('name'))
EOF
```
Expected: 8 classes, ~56 subclasses, ~2,127 families, ≤596 subfamilies (bimodal — expect a few hundred), ~286 specificities with descriptions; `failed_families` empty; `3.F` named. Eyeball ~10 random specificity names for garbage (unstripped HTML, citation debris).

- [ ] **Step 4: Commit**

```bash
git add cache/data/tcdb/tcdb_names.json
git commit -m "data(tcdb): scraped node names artifact (browse + 351 family pages, 2026-08-18)"
```

---

### Task 7: Rebuild step-6 artifacts

- [ ] **Step 1: Rebuild** (no network)

```bash
bash scripts/prepare_data.sh --steps 6 --force
tail -20 logs/prepare_data_step6.log
```
Expected: completes ~1 min; the unscraped-family warning does NOT fire (all kept families scraped); no empty-rebuild guard errors.

- [ ] **Step 2: Verify the diff is names-only**

```bash
git diff --stat cache/data/tcdb/ cache/data/kegg/
uv run python - <<'EOF'
import json
h = json.load(open('cache/data/tcdb/tcdb_hierarchy.json'))
unnamed = sum(1 for v in h.values() if not v.get('name'))
print('hierarchy entries:', len(h), '| unnamed:', unnamed)
print('with description:', sum(1 for v in h.values() if v.get('description')))
print("class 1:", h['1']['name'])   # expect "Channels/Pores"
print("class 6 present:", '6' in h, h.get('6', {}).get('name'))
EOF
```
Expected: `kegg_data.json` **unchanged** (names don't touch KEGG); `tcdb_pruned.json` unchanged (names don't affect pruning); `tcdb_hierarchy.json` gains names/descriptions only — entry COUNT identical to before. Class 6 is now named (it exists in the full hierarchy; still not kept).

- [ ] **Step 3: Run unit tests, commit**

```bash
pytest -m "not slow and not kg" -q
git add cache/data/tcdb/tcdb_hierarchy.json
git commit -m "data(tcdb): rebuild hierarchy with scraped names (step 6 --force)"
```

---

### Task 8: Docker rebuild + KG validation (Definition of Done)

- [ ] **Step 1: Rebuild the graph**

```bash
docker compose stop app deploy       # deploy holds the DB lock
docker compose up build import post-process   # full chain; ~long
docker compose up -d deploy app
```
Check `output/import.status` is 0 and `output/import.report` is empty (zero dangling — names add no edges, so any regression is a bug).

- [ ] **Step 2: DoD checks against the live graph**

```bash
uv run python - <<'EOF'
from neo4j import GraphDatabase
drv = GraphDatabase.driver("bolt://localhost:7687")
with drv.session() as s:
    # unnamed count — was 916
    rec = s.run("MATCH (t:TcdbFamily) WHERE t.name = t.tcdb_id "
                "RETURN t.level_kind AS k, count(*) AS n ORDER BY k").data()
    print("unnamed by level:", rec)
    # full-text search through a scraped name (DoD: verified via the INDEX)
    hits = s.run("CALL db.index.fulltext.queryNodes('tcdbFamilyFullText', "
                 "'sugar porter') YIELD node, score "
                 "RETURN node.tcdb_id, node.name LIMIT 5").data()
    print("fulltext 'sugar porter':", hits)
    hits2 = s.run("CALL db.index.fulltext.queryNodes('tcdbFamilyFullText', "
                  "'ammonia') YIELD node, score "
                  "RETURN node.tcdb_id, node.name LIMIT 5").data()
    print("fulltext 'ammonia':", hits2)
EOF
```
Expected: unnamed drops from 916 to only the genuinely-unnamed-upstream subfamilies (spec 6d bimodal — expect roughly 100–400, all `tc_subfamily`); subclass + specificity unnamed ≈ 0; both full-text queries return named nodes (e.g. `2.A.1.1` The Sugar Porter (SP) Family; `1.A.11.*` ammonia entries). Record the final unnamed number for the docs task.

- [ ] **Step 3: KG validity + snapshot**

```bash
pytest -m kg -v
```
If `tests/kg_validity/test_snapshot.py` fails on TcdbFamily `name` properties (snapshot predates the rename), regenerate and re-run:
```bash
uv run python tests/kg_validity/generate_snapshot.py
pytest -m kg -v
```

- [ ] **Step 4: Commit**

```bash
git add tests/kg_validity/snapshot_data.json
git commit -m "test(kg): regenerate snapshot after TCDB node naming"
```

---

### Task 9: Documentation

**Files:**
- Modify: `CHANGELOG.md` (`[Unreleased]`), `CLAUDE.md`, `docs/kg-changes/tcdb-two-source-upgrade.md`, `.claude/skills/tcdb-diamond/SKILL.md`, `docs/superpowers/specs/2026-08-12-tcdb-node-names-design.md`

- [ ] **Step 1: CHANGELOG.md** — under `[Unreleased]` → Added:

```markdown
- **TCDB node names (T6)**: 916 previously-unnamed `TcdbFamily` nodes (all
  subclasses, subfamilies, specificities + 3 families) now carry scraped
  upstream names; `tc_specificity` nodes gain a sparse 400-char `description`;
  `tcdbFamilyFullText` covers it. New committed artifact
  `cache/data/tcdb/tcdb_names.json` built by
  `python -m multiomics_kg.download.scrape_tcdb_names` (scoped per-family
  scrape, resumable, ~352 requests); step 6 consumes it and warns when a
  strain onboarding introduces families the scrape hasn't covered. The two
  wrong hardcoded class names are fixed (`Channels/Pores`, `Accessory Factors
  Involved in Transport`) and class 6 (MTM) is named. Unnamed remainder:
  <N — fill from Task 8> genuinely-unnamed-upstream subfamilies (render as
  bare TC id by design — spec 5c).
```

- [ ] **Step 2: CLAUDE.md** — in the TcdbFamily nodes bullet: update the properties list (`description` sparse), note that names are scraped (`tcdb_names.json`, scrape script, guard), and update the `tcdbFamilyFullText` field list in the post-import indexes bullet. Add `scrape_tcdb_names` to the step-6 section: names artifact is a step-6 *input* built separately (like `refresh_mnx.sh` for MNX).

- [ ] **Step 3: `docs/kg-changes/tcdb-two-source-upgrade.md`** — add a short section (after §6 "Node-set and count changes") titled "Node names (T6, 2026-08-18)": what changed for consumers (60% of ontology was reachable only by exact TC id; now searchable), the `description` property, the `t.name = t.tcdb_id` fallback-detection predicate, and the unnamed-remainder number from Task 8. Base wording on the CURRENT post-vocab-merge version of the doc (it renamed `substrate_depth` values etc.).

- [ ] **Step 4: `.claude/skills/tcdb-diamond/SKILL.md`** — add a "Refreshing TCDB names (T6)" section: when to run (new strain introduces unscraped families — step 6 warns; or a TCDB release), the command, `--families` top-up mode, the cadence/resume behavior, and that `tcdb_names.json` is committed while `raw/pages/` is not.

- [ ] **Step 5: Spec status** — update `docs/superpowers/specs/2026-08-12-tcdb-node-names-design.md` status line to **IMPLEMENTED (2026-08-18)** with the final unnamed count.

- [ ] **Step 6: Commit**

```bash
git add CHANGELOG.md CLAUDE.md docs/kg-changes/tcdb-two-source-upgrade.md .claude/skills/tcdb-diamond/SKILL.md docs/superpowers/specs/2026-08-12-tcdb-node-names-design.md
git commit -m "docs(tcdb-names): changelog + CLAUDE.md + kg-changes + skill docs for T6"
```

---

## Self-Review

- **Spec coverage:** §3a browse parsing → Task 2; §3b family-page parsing incl. 5-column change + `</strong>` capture + sidebar scoping → Task 2; §4 scoped scrape + full-browse-layer + names artifact + guard → Tasks 3–4; §5b name/description shape → Task 1 (caps 150/400); §5c bare-ID fallback → no code needed (adapter fallback already exists; verified by Task 8 unnamed count); §5d hygiene → Task 3 (delay/backoff/abort/resume) + Task 6 (execution); §7 sanitization → Task 5 adapter `_clean_str` on description (name already sanitized), vocabulary contract → no categorical props added; §7 stale docstring fix → Task 4 step 3.7; §8 DoD → Task 8 (unnamed drop, empty import.report, both test suites, full-text query); §9 docs → Task 9. Class-name fixes (§3a table) come free via artifact-over-hardcoded (Task 4). Gap check: `--refetch-raw` unification (spec 6f) is BLOCKED on the two broken substrate CGIs and intentionally out of scope here — the scraper is independent of them.
- **Placeholder scan:** one intentional fill-in: `<N — fill from Task 8>` in the CHANGELOG text — the number does not exist until the graph is rebuilt; Task 8 step 2 explicitly records it.
- **Type consistency:** `scoped_families(set[str]) -> dict[str,str]`; `parse_family_page -> tuple[dict,dict]` consumed as `(subfams, systems)` in Task 3's main and tests; `_load_tcdb_names -> tuple[dict, list]` used in Task 4 guard; artifact key paths (`meta.scraped_families`, `names.<id>.name/.description`) consistent across Tasks 3, 4, 6.
