# T6 — Naming the unnamed TcdbFamily nodes

**Status: READY TO IMPLEMENT — all design questions resolved. Live-site
verification complete (2026-08-17); scope revised to a SCOPED scrape (section
4, resolves 5a); 5b resolved empirically, 5c and 5d decided 2026-08-18.
Aligned 2026-08-18 with the controlled-vocabulary contract on main (section
7). Remaining external item: the two substrate CGIs were still HTTP 500 on
2026-08-17 (6f) — re-probe before a full `--refetch-raw`.**

`tcdb.org` served an SDSC "planned maintenance" stub for the whole site on
2026-08-12, so the original research was verified against **Wayback snapshots**
only. On **2026-08-17 the site was back online** and every check in section 6
was run against the live site; findings are recorded inline there. Headlines:

- **The scrape is far cheaper than designed**: `result.php?tc=<any 4-part>`
  returns the **entire family** (all subfamily headers + all 5-part system
  names) in one page, so the unit of work is the *family* (≤2,208 requests),
  not the subfamily (3,550) — see 6a/6c. The scoped-scrape decision (section
  4) then cuts it to **~352 requests** (the 351 families carrying kept deep
  nodes, plus `browse.php`).
- The 2015 header template survives; the column layout changed (6a).
- Subfamily naming upstream is **bimodal per family**: big superfamilies are
  100% named, many small families have no subfamily names at all — 5c's
  fallback question is confirmed real (6d).
- **New availability problem**: `getSubstrates.py` and `listSuperfamilies.py`
  (2 of the 4 TSV endpoints step 6 already depends on) return persistent
  HTTP 500 as of 2026-08-17 — a same-pass `--refetch-raw` (6f) is currently
  impossible for those two files, and `/public/` carries no static fallback
  for them.

Nothing is left to decide — implement per sections 4 and 5, following the
CLAUDE.md 5-phase workflow (section 7).

---

## 1. Problem, re-verified

916 of 1,515 `TcdbFamily` nodes render as their bare numeric ID because the node
has no name. Re-measured on the live graph on 2026-08-12 — identical to the
original report, so the graph has not drifted:

| `level_kind` | total | unnamed (`t.name = t.tcdb_id`) |
|---|---|---|
| `tc_class` | 7 | 0 |
| `tc_subclass` | 34 | **34** (all) |
| `tc_family` | 592 | 3 |
| `tc_subfamily` | 596 | **596** (all) |
| `tc_specificity` | 286 | **286** (all) |

The `tcdbFamilyFullText` index is built on `t.name`, so 60% of the ontology is
reachable only by exact TC ID.

Root cause is in `build_tcdb_hierarchy` in
`multiomics_kg/download/build_kegg_metabolism_xrefs.py`: `ensure_subclass`,
`ensure_subfamily` and `ensure_specificity` all hardcode `"name": ""`, and no
description source for those three levels is downloaded at all.
`tcdb_adapter.py:308` then falls back to `entry.get("name") or tcdb_id`.

### Scale in the full hierarchy vs the pruned graph

Names have to be attached in `tcdb_hierarchy.json`, which is the **full**
pre-prune reference. The pruned graph is a small slice of it:

| `level_kind` | full hierarchy | unnamed in full | kept (graph) | unnamed kept |
|---|---|---|---|---|
| `tc_class` | 8 | 1 (class `6`) | 7 | 0 |
| `tc_subclass` | 52 | 52 | 34 | 34 |
| `tc_family` | 2,208 | 29 | 592 | 3 |
| `tc_subfamily` | 3,550 | 3,550 | 596 | 596 |
| `tc_specificity` | 19,743 | 19,743 | 286 | 286 |

## 2. Feasibility verdict

**Feasible, but the only source is HTML scraping — TCDB publishes no bulk file
for these three levels.**

The complete published download set was enumerated from the archived
`tcdb.org/download.php` **and** the `/public/` directory index (both fetched):

| Endpoint | Content |
|---|---|
| `/public/tcdb` | all proteins, FASTA (13 M) |
| `/public/tcDoms.tar.gz` | tcDoms |
| `/cgi-bin/substrates/getSubstrates.py` | TC systems → substrates + ChEBI |
| `/cgi-bin/projectv/public/families.py` | **family definitions (3-part only)** |
| `/cgi-bin/substrates/listSuperfamilies.py` | systems/subfamilies/families → superfamilies |
| `/cgi-bin/projectv/public/acc2tcid.py` | UniProt/RefSeq acc → TC system |
| `/cgi-bin/projectv/public/refseq.py` | RefSeq → TCID |
| `/cgi-bin/projectv/public/go.py` | GO → TC system |
| `/cgi-bin/projectv/public/pdb.py` | PDB → TC system |
| `/cgi-bin/projectv/public/pfam.py` | Pfam → TC system |
| `/public/` extras | `betabarrel`, `human.csv`, `human_specific.csv`, `tcdb.dr`, `mim.tsv`, `pmid.tsv` |

Step 6 already consumes four of these. **None** carries a subclass, subfamily or
specificity name. TCDB's own FAQ states the definitions of classes, subclasses,
families and subfamilies are seen by *browsing* the TC system — consistent with
this.

Live re-check 2026-08-17: the CGI endpoints still exist but the **live
`/download.php` no longer links any of them** (it only links site pages); the
`/public/` index is the discoverable bulk inventory now and has grown static
TSV copies — `families.tsv`, `acc2tcid.tsv`, `refseq.tsv`, `pdb.tsv`,
`mim.tsv`, `pmid.tsv` (dated 2024-07-01) and `go.tsv`, `pfam.tsv` (2025-12-15),
plus the FASTA `tcdb` (2025-08-04). **No static substrates or superfamilies
file exists**, and both substrate CGIs are down — see 6f.

## 3. Sources found, with confidence

### 3a. `browse.php` — one request, VERIFIED LIVE 2026-08-17

Live page fetched 2026-08-17: 200, **686,668 bytes** (10× the 68,796-byte
2025-07-31 snapshot — the page now inlines the full class/subclass/family tree
as a treeview). Markup is uniform at all three levels:

```html
<div rel="1.A.1" class="entry" ...>
  <div class="tcid name">&nbsp;1.A.1:&nbsp;The Voltage-gated Ion Channel (VIC) Superfamily </div>
</div>
```

Parse `rel="<tc id>"` + the `tcid name` div text. A single page yields:

- **all 8 class names**, including class `6`
- **56 subclass names** (up from 54 in the 2025-07 snapshot) — covers
  **34/34 kept** and **52/52 in the full hierarchy**, `3.F` included (6b)
- **2,127 family names, zero blank** — fixes the 3 of 3 unnamed kept families

No 4-part IDs appear on the page, so it cannot help with subfamilies.

Class names as published (re-confirmed identical on the live 2026-08-17 page):

```
1  Channels/Pores
2  Electrochemical Potential-driven Transporters
3  Primary Active Transporters
4  Group Translocators
5  Transmembrane Electron Carriers
6  Membrane Transporter Metabolons (MTM)
8  Accessory Factors Involved in Transport
9  Incompletely Characterized Transport Systems
```

Subclass samples: `1.A` α-Type Channels · `1.B` β-Barrel Porins · `1.C`
Pore-Forming Toxins (Proteins and Peptides) · `1.D` Non-Ribosomally Synthesized
Channels · `1.E` Holins · `1.F` Vesicle Fusion Pores · `3.F` Artificial
(Unusual energy sources) Active Transporters.

**This obsoletes `_TC_CLASS_NAMES` rather than extending it.** Two of the seven
hardcoded names disagree with upstream:

| class | hardcoded | upstream |
|---|---|---|
| 1 | Channels **and** Pores | Channels**/**Pores |
| 8 | Auxiliary Transport Proteins | Accessory Factors Involved in Transport |
| 6 | *(missing)* | Membrane Transporter Metabolons (MTM) |

Class `6` is real, not a data artifact — it holds 9 nodes in our hierarchy
(`6.A.1` PTS-Glycolytic metabolon, `6.A.2` Plant Glycolysis-VDAC-Krebs Cycle
Metabolon, `6.A.3`, `6.A.4`, `6.A.5`, `6.B.1`). It is not currently
gene-annotated so it never reaches the graph.

### 3b. `search/result.php?tc=<4-part>` — the subfamily + specificity source, VERIFIED LIVE 2026-08-17

The 2015-template risk is **retired**. Live fetches (`tc=2.A.1.1`,
`tc=1.A.11.1`, `tc=9.B.1.1`, `tc=3.A.1.1`, `tc=2.A.6.1`, `tc=2.A.7.1`,
`tc=3.A.3.1`, `tc=2.A.3.1`) confirm the structure, with two changes vs 2015:

1. **A 4-part query returns the ENTIRE family**, not one subfamily:
   `tc=2.A.1.1` came back with all 91 subfamily anchors of `2.A.1` and every
   5-part system row (1.08 MB). So the unit of work is one request **per
   family** — any known 4-part child works as the query key. This cuts the
   full-hierarchy job from 3,550 to **≤2,208 requests**.
2. The table columns changed from 4 (`TCID | Name | Organismal Type | Example`)
   to **5** (`TCID | Name | Domain | Kingdom/Phylum | Protein(s)`). The
   `colspan="4"` on header rows survives unchanged (upstream never updated it —
   don't key the parser on it).

What survives from 2015, byte-for-byte:

- `<table id="result-cluster">`
- the named-subfamily section header row:
  ```html
  <tr><td colspan="4" id="right-border"><strong><A id="2.A.1.1"></A>2.A.1.1:&nbsp;&nbsp;The Sugar Porter (SP) Family</strong></td></tr>
  ```
  **Parser caveat**: subfamily names can contain inline markup
  (`The Drug:H<sup>+</sup> Antiporter-1 …`, `<sub>`, `<em>`), so capture up to
  `</strong>`, never `[^<]*`.
- every child system's curated `Name` cell, e.g. `2.A.1.1.1` → *"Galactose:H⁺
  symporter, GalP. Also transports glucose, xylose, …"* — with `<sup>`, `<em>`,
  `<i>`, HTML entities, and `<a class="reflink">` citation links to unwrap.

**Unnamed subfamilies emit no header row at all** (not an empty one): in
families without curated subfamily names the 4-part ID appears only as a bare
`<A id="…"></A>` anchor — or not at all (`3.A.3` shows one anchor for 17
subfamilies). Absence of a `<strong>` header == no upstream name. See 6d.

**Scoping caveat**: pages embed 5-part IDs from outside the queried family
(a sidebar module — `tc=2.A.1.1`'s page carries 1,041 distinct 5-part IDs).
Parse only inside `result-cluster`.

Revised cost at full-hierarchy scope: **≤2,208 requests** (one per family;
observed pages 22 KB–1.36 MB, ~1 s server time each), roughly **0.3–0.7 GB,
≈40 min at 1 req/s.**

### 3c. `tcfamilybrowse.php?tc=<3-part>` — specificity only, SUPERSEDED by 3b

Still live (`?tc=1.A.11` → 200, 82,527 bytes on 2026-08-17): the "Format for
printing" view with the family narrative plus a flat 5-part `Examples:` table.
But since a single `result.php` family fetch (3b) now yields subfamily *and*
specificity names at the same request count (one per family), this endpoint no
longer buys anything. Recorded for completeness only.

Class- and subclass-level variants are confirmed **unsupported** on the live
site: `?tc=1` returns a 43-byte stub and `?tc=1.A` a 606-byte references-only
fragment; `?tc=2.A.1.1` (4-part) also returns a 43-byte stub. There is no
8-request shortcut (6c).

### 3d. Local `tcdb.raw.faa` — free, but semantically wrong

`~/tools/TCDB/DB/tcdb.raw.faa` is already on disk for `/tcdb-diamond`.
24,281 headers, format `>gnl|TC-DB|<acc>|<5-part TCID> <description>`;
19,217 distinct 5-part TCIDs carry a description, covering **19,174/19,743
(97.1%)** of full specificity nodes and **285/286 kept**.

But the descriptions are raw UniProt FASTA text — per **protein**, not per
transport **system**, complete with `OS=` / `GN=` / `PE=` / `SV=` tokens, and
multi-protein systems yield several competing strings:

```
1.A.13.2.2 → Uncharacterized protein sll0103 OS=Synechocystis sp. (strain PCC 6803 / Kazusa) GN=sll0103 PE=4 SV=1
1.A.12.3.2 → Putative enzyme OS=Escherichia coli (strain 55989 / EAEC) GN=yfcF PE=4 SV=1
```

A `tc_specificity` node represents a transport system, so this is the wrong
noun. Keep it only as a documented last-resort fallback for specificity nodes
the scrape misses.

## 4. Direction agreed: SCOPED scrape (decided 2026-08-17)

`browse.php` **plus** a per-family `result.php` scrape **scoped to the kept
graph**. Two clearly separated fetchers, both writing cached files into the
gitignored `cache/data/tcdb/raw/`, both gated behind `--refetch-raw` so
`--force` stays a ~1 min no-network rebuild.

- **Classes / subclasses / families**: full coverage from the single
  `browse.php` fetch (it names every one, zero blanks — 3a). Scoping never
  applied to these levels. Delete `_TC_CLASS_NAMES`, fix class `6`, correct the
  two wrong class names.
- **Subfamilies / specificities**: scrape only the families that carry kept
  deep nodes. Measured against the current `tcdb_pruned.json`: the 596 kept
  subfamilies + 286 kept specificity nodes fall in **351 distinct families**,
  so the scrape is **~352 requests total (~6–10 min)** instead of 2,208
  (~40 min, 0.3–0.7 GB).
- **Artifact**: scoped names are ~900 strings (~100 KB), so the 5a size
  dilemma disappears. Recommended layout: one small committed
  `cache/data/tcdb/tcdb_names.json` holding all scraped names (all levels),
  consumed by `build_tcdb_hierarchy` at step-6 time — keeps
  `tcdb_hierarchy.json`'s existing shape untouched. Inlining into
  `tcdb_hierarchy.json` is acceptable too; minor implementation choice.

**The gene-set-coupling trade-off, accepted deliberately.** Scoping makes the
committed names gene-set-dependent: onboarding a new strain can introduce kept
TC nodes whose family was never scraped. This was originally grounds for
rejection (2026-08-12), reversed because:

1. **Precedent** — `interpro_reference.json`'s `description` field is already
   observed-only/corpus-coupled for the same size-gate reason, with the
   documented "re-run step 9 after a new strain lands" implication. Same
   pattern, same maintenance posture; add-a-strain already loops back through
   prepare_data.
2. **Self-healing at trivial cost** — step 6 already computes `kept_tcdb_ids`;
   extend the `dc76076c` guard to flag kept 4/5-part IDs whose family is
   absent from the names artifact (loud staleness, never silent), and the
   top-up fetch for one new strain is typically a handful of family pages.
3. **Failure mode is the status quo** — an unnamed new node renders its bare
   TC ID, which is exactly today's behaviour for all 916, until the next
   `--refetch-raw`.

What is given up: the "complete per upstream release" property of the names
artifact, and instant naming for new strains onboarded while tcdb.org is down
(demonstrated risk — names lag until it recovers).

## 5. Design questions — ALL RESOLVED

**5a. Committed-artifact size. — RESOLVED by the scoped scrape (section 4).**
The original dilemma was 19,743 specificity names (~3 MB) nearly doubling a
committed 4.5 MB `tcdb_hierarchy.json` for names of which only 286 are ever
used. Scoped to kept nodes, the artifact is ~900 names ≈ 100 KB; recommended
layout is a separate committed `tcdb_names.json` (see section 4), with
inline-in-hierarchy acceptable as a minor variant.

**5b. Shape of a `tc_specificity` name. — RESOLVED (2026-08-18, empirically):
citation-strip → first sentence → ~150-char word-boundary cap for `name`; full
cleaned text in `description`.** Tested on 2,230 real system names parsed from
six live family pages (2.A.1, 3.A.1, 1.A.11, 2.A.6, 3.A.3, 9.B.1):

- Raw names: median 150 chars, p90 747, max 27,415 (some cells are essays).
  Zero empty names.
- **First sentence** (after stripping `(Author, YYYY)` citation
  parentheticals, including nested/multi-ref groups): almost always a
  complete, self-standing label — *"Glycerol-P:Pi antiporter."*,
  *"α-Hemolysin exporter."*, *"Eye pigment precursor transporter, White."*
  Median 68 chars; 15% still exceed 150 (parenthetical-heavy entries), max
  705 — hence the cap.
- **First-N-words was tested and REJECTED**: at N=8 it chops mid-phrase and
  dangles open parentheses (*"…transporter of 560 aas and 12…"*) on anything
  longer than one sentence.
- Parser cautions: sentence-splitting must not break on abbreviations
  (*"E. coli"*, *"et al."*, *"e.g."*, *"2.0 Å"*, single-letter genus
  initials); tag-stripping must tolerate malformed fragments (one live cell
  carries a dangling unclosed `<br /`).
- `description` follows the InterproEntry precedent (400-char-capped there;
  pick the cap at implementation time — these are ≤27 KB pre-clean).

**5c. Fallback when upstream has no name. — RESOLVED (2026-08-18): keep the
bare TC ID.** The question is real (6d: naming is bimodal per family — big
superfamilies 100% named, many small families have **no** subfamily names at
all, so a fraction of the 596 kept subfamilies stays nameless after the
scrape). Decision: an honest bare ID beats a synthesized label like
`"2.A.1.1 (subfamily of The Major Facilitator Superfamily (MFS))"` — the
parent's name is one `Tcdb_family_is_a_tcdb_family` hop away, and synthesized
text would pollute the full-text index with the parent's tokens. No
`name_source` provenance property either: `t.name = t.tcdb_id` already *is*
the fallback marker (a property restating that predicate is exactly the
derivable-fact materialization rule R3 warns against), so no
`controlled_vocabularies.yaml` entry is needed.

**5d. Scrape hygiene. — RESOLVED (2026-08-18): standalone script, slow and
easy, resumable.** The scraper is a **separate script**, decoupled from step
6's fast no-network rebuild — code at
`multiomics_kg/download/scrape_tcdb_names.py` (repo convention: it's
upstream-reference fetching, sibling to `download_metabolism_reference.py`,
importable/unit-testable), with its invocation documented in the tcdb skill's
SKILL.md as the "refresh TCDB names" workflow. Step 6 only *consumes* the
committed names artifact.

Cadence — deliberately gentler than the server needs, because a rarely-run
offline script gains nothing from speed:

- Single-threaded, **one request per 2–3 s** (scoped 352 requests ≈ 15–20 min).
- Polite User-Agent carrying a contact address.
- Retry with exponential backoff on 5xx/timeouts; after a few consecutive
  failures **abort and resume later** — never hammer a struggling host (TCDB
  has demonstrated multi-day maintenance windows).
- **Fetch and parse are separate stages**: raw HTML pages cached under
  `cache/data/tcdb/raw/pages/` (gitignored); a re-run skips already-fetched
  families, and re-parsing never re-fetches. A page that 404s is recorded and
  skipped, not fatal.
- The `dc76076c` guard precedent extends to the names artifact: step 6 raises
  if a rebuild would come out empty beside a populated committed names file,
  and flags kept TC IDs whose family is missing from it (section 4).

## 6. Live-site verification — ALL CHECKS RUN 2026-08-17

**6a. Does the 2015 subfamily-page template still hold?** **Mostly yes — parser
design survives with two adjustments.** `result.php?tc=2.A.1.1` fetched live:
`table id="result-cluster"` and the `<strong><A id="…"></A>4-part:&nbsp;&nbsp;name</strong>`
header row are intact. Changes: (1) system rows now have **five** columns
(`TCID | Name | Domain | Kingdom/Phylum | Protein(s)`), not four; (2) subfamily
names can contain inline `<sup>`/`<sub>`/`<em>` markup, so capture to
`</strong>`. Bonus finding: the page covers the **whole family**, not just the
queried subfamily — full details in 3b.

**6b. Is `3.F` on the current `browse.php`?** **Yes** — `3.F: Artificial
(Unusual energy sources) Active Transporters`, with families `3.F.1` (LAMM),
`3.F.2` beneath it. Live browse.php carries 56 subclasses → **52/52 coverage**
of the full hierarchy, no permanent gap. (The page itself was rebuilt upstream
as a 686 KB treeview — see 3a for the new markup.)

**6c. Probe for a cheaper bulk path.** **Probed; the winner is 3b's own
whole-family behaviour (≤2,208 requests), and nothing cheaper exists.**
`tcfamilybrowse.php?tc=1` → 43-byte stub; `?tc=1.A` → 606-byte references-only
fragment; `?tc=<4-part>` → 43-byte stub (only 3-part works). `result.php` at
subclass/class level returns narrative-description pages with no child listing.
Live `/download.php` no longer lists CGI endpoints at all; `/public/` gained
static TSVs (`go.tsv` + `pfam.tsv` dated 2025-12-15, FASTA 2025-08-04, rest
2024-07-01) but nothing that names subfamilies or specificities.

**6d. Blank-name rate.** **Bimodal per family, not a uniform rate.** Sample of
8 families: the five large/curated ones are fully named — `2.A.1` 91/91,
`3.A.1` 98/98, `2.A.7` 34/34, `2.A.3` 15/15, `2.A.6` 9/9 — while three
smaller ones have **zero** named subfamilies (`1.A.11` 0/3, `9.B.1` 0/2,
`3.A.3` 0/17; P-type ATPase!). Unnamed subfamilies emit no header row at all.
Specificity (5-part) names were present for every system row in all 8 pages.
Answers 5c: a fallback story is genuinely needed.

**6e. Politeness.** No `robots.txt` (404), no stated rate limit
(Apache/2.4.52 Ubuntu). Observed ~1 s server time per family page, sizes
22 KB (`9.B.1`) – 1.36 MB (`3.A.1`). At 1 req/s a full-hierarchy scrape would
be ≈40 min and roughly 0.3–0.7 GB (replaces the ~400 MB / 1 h per-subfamily
estimate); the scoped scrape (section 4) is **~352 requests, ≈6–10 min**.
Given the site just came out of maintenance, keep 1 req/s + retries.

**6f. Re-download the four existing TSVs.** **Currently HALF-BLOCKED.** Live
status 2026-08-17: `families.py` 200 (138 KB) · `acc2tcid.py` 200 (483 KB) ·
`pfam.py` 200 (624 KB) · `go.py` 200 (2.5 MB) — but **`getSubstrates.py` and
`listSuperfamilies.py` both return persistent HTTP 500** (retried; stable
Apache error, `webdb-help@ucsd.edu` as contact). `/public/` has no static
substrates/superfamilies fallback. So a same-release refresh of all raw inputs
is impossible right now: either (a) wait for the substrate CGIs to recover
before running `--refetch-raw`, or (b) accept mixed vintages — the design
preference is (a); the empty-rebuild guard (`dc76076c`) protects against a 500
silently emptying `substrates.tsv` in the meantime. Re-check the two endpoints
before implementation.

## 7. Constraints carried from the task brief

- Follow the CLAUDE.md 5-phase workflow (Scope / Implement / Test / Review /
  Document); no phase skipped.
- **Sanitise every string property** — `.replace("'", "^").replace("|", "")`.
  TCDB names carry apostrophes and Greek letters (α, β, Δ, μ) as a matter of
  course, and `|` is the array delimiter; an unescaped one silently splits
  values. That exact bug was just fixed for UniProt in `c043241d` — do not
  reintroduce it. HTML entities and `<sup>`/`<em>` markup must be unwrapped
  before sanitising.
- **The controlled-vocabulary contract (landed on main 2026-08-18, after this
  spec was written) is binding** — see `docs/kg-changes/vocabulary-contract.md`
  and `docs/superpowers/specs/2026-08-16-vocabulary-contract-design.md`:
  - **R5** subsumes the old "BioCypher mishandles booleans" note: no native
    `bool` properties ever; a two-state fact is a meaningful categorical
    string (`bool_string`), and a `bool` in `schema_config.yaml` now *fails
    the vocabulary test suite* rather than just silently misbehaving.
  - Any **new categorical property** T6 introduces (e.g. a `name_source` /
    fallback marker out of 5c) must be declared in
    `config/controlled_vocabularies.yaml` — KG-minted values in lowercase
    `snake_case` (R1) — or the four-gate vocabulary test fails the build scan.
  - Scraped TCDB names themselves are free-text `name` properties, not
    vocabulary values — no declaration needed, R1's "external terms verbatim"
    spirit applies (keep upstream wording, post-sanitisation).
  - `TcdbFamily.level_kind` is already a declared closed vocabulary; T6 adds
    no new values to it.
- `cache/data/tcdb/tcdb_hierarchy.json` and `tcdb_pruned.json` are **committed**
  artifacts. The `dc76076c` empty-rebuild guard is intentional; extend it.
- Iteration loop: `bash scripts/prepare_data.sh --steps 6 --force` (~1 min, no
  network). `--refetch-raw` only when pulling upstream.
- **Fix the stale module docstring** (`build_kegg_metabolism_xrefs.py` ~line 22).
  It claims step 6 downloads *"the 3 TCDB reference TSVs (tc_classes,
  tc_subclasses, families)"*. None of those three names match what the code
  reads, there are four files, and no subclass file exists at all.

## 8. Definition of done

- Named nodes at subclass, subfamily and specificity level (source confirmed
  live, 6a — subfamily coverage is bounded by upstream's bimodal naming, 6d);
  `t.name = t.tcdb_id` count drops substantially from 916.
- The step-6 guard (extended `dc76076c` pattern) flags kept 4/5-part TC IDs
  whose family is absent from the committed names artifact, so gene-set
  staleness after a strain onboarding is loud, never silent.
- Full Docker rebuild with `output/import.report` still empty (zero dangling).
- `pytest -m "not slow and not kg"` and `pytest -m kg` both green.
- Names searchable through the `tcdbFamilyFullText` index, verified with a real
  full-text query — not by inspecting the property.

## 9. Documentation targets

The explorer has not consumed the current TCDB upgrade and it is unreleased, so
fold this into the **same** item rather than creating a new one:

- `docs/kg-changes/tcdb-two-source-upgrade.md` — see its §7 for the pattern.
  (Note: this doc was restructured by the 2026-08-18 vocabulary-contract merge
  — `is_promiscuous` deleted, `substrate_depth` values renamed to
  `most_specific`/`inherited`, scores rescaled — rebase T6's wording on the
  current version, not on pre-merge quotes.)
- `CHANGELOG.md` `[Unreleased]`
- the `TcdbFamily` bullet in `CLAUDE.md`
- `config/controlled_vocabularies.yaml` — only if T6 adds a categorical
  property (see sections 5c and 7); plain `name`/`description` strings need no
  entry

## 10. Related deferred work

**T4** — `TcdbFamily.gene_count` is a SUBTREE count (`*0..`) while
`InterproEntry.gene_count` is direct. A warning is in CLAUDE.md; not fixed. Out
of scope here, but both touch the same post-import block, so worth batching if
T6 lands first.
