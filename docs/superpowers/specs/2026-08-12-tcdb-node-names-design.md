# T6 — Naming the unnamed TcdbFamily nodes

**Status: BLOCKED on upstream availability. Research complete, design not settled.**

`tcdb.org` served an SDSC "planned maintenance" stub for the whole site on
2026-08-12 (nginx, `content-length: 1732`, body reads *"We'll be back soon.
Planned maintenance is in progress."*). That includes the four TSV endpoints
prepare_data step 6 already depends on. Nothing could be fetched live, so every
upstream claim below was verified against **Wayback snapshots**, not the live
site, and the snapshot age varies a lot per source. Section 6 lists what must be
re-checked against the live site before any code is written.

Resume by working section 6 top to bottom, then answering the four open design
questions in section 5.

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

## 3. Sources found, with confidence

### 3a. `browse.php` — one request, high confidence

Snapshot **2025-07-31** (`20250731190452`, 68,796 bytes). Recent, and the
structure matched the 2014 snapshot too, so the template is stable across a
decade. A single page yields:

- **all 8 class names**, including class `6`
- **54 subclass names** — covers **34/34 kept**, 51/52 in the full hierarchy
  (only `3.F` absent, present in our May-2026 TSVs; see 6b)
- **2,061 family names** — fixes 6 of the 29 unnamed families in the full
  hierarchy, **3 of 3 unnamed kept families**

No 4-part IDs appear on the page, so it cannot help with subfamilies.

Class names as published:

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
Channels · `1.E` Holins · `1.F` Vesicle Fusion Pores.

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

### 3b. `search/result.php?tc=<4-part>` — the only subfamily source, LOW confidence

Snapshot **2015-03-28** (`20150328052704`, 144,084 bytes) — **an 11-year-old
template, and this is the source carrying 882 of the 916 unnamed nodes.** It is
the single largest risk in this task.

One fetch names **both** remaining levels:

- the subfamily itself, as a section header row:
  ```html
  <table id="result-cluster">
  <tr><th>TCID</th><th>Name</th><th>Organismal Type</th><th class="right-border">Example</th></tr>
  <tr><td colspan="4" id="right-border"><strong><A id="2.A.1.1"></A>2.A.1.1:&nbsp;&nbsp;The Sugar Porter (SP) Family</strong></td></tr>
  ```
- every child system's curated `Name`, e.g. `2.A.1.1.1` → *"Galactose:H⁺
  symporter (also transports xylose) (Hernández-Montalvo et al., 2001). Relative
  substrate affinities of wild-type and mutant forms of the E. coli sugar
  transporter GalP have been determined by solid-state NMR (Patching et al.,
  2008)."*

Name cells contain `<sup>`, `<em>`, `<i>`, HTML entities, and
`<a class="reflink">` citation links that need unwrapping.

Cost at full-hierarchy scope: **3,550 requests, ~400 MB, ≈1 h at 1 req/s.**

### 3c. `tcfamilybrowse.php?tc=<3-part>` — specificity only, cannot replace 3b

Snapshot **2016-01-05** (`?tc=1.A.11`, 19,787 bytes). The "Format for printing"
view: family narrative plus a **flat** `Examples:` table (`TC# | Name |
Organismal Type | Example`) of 5-part systems with their curated names
("Ammonia transporter and regulatory sensor, AmtB (…)").

**No subfamily headings** — 4-part IDs appear only as prefixes of 5-part IDs.
So it gives specificity names in 2,208 requests over much smaller pages (1–40 KB,
~40 MB total) but leaves all 3,550 subfamilies unnamed. Useful only as a
cheaper specificity-only path, not as the primary source.

No archived snapshot exists for a subclass- or class-level variant
(`?tc=1.A`, `?tc=1`), suggesting those are unsupported. Worth one live probe
anyway (6c) — a class-level variant would collapse the whole job to 8 requests.

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

## 4. Direction agreed so far

`browse.php` **plus** the full per-subfamily scrape (option chosen 2026-08-12):
name every level across the whole hierarchy, delete `_TC_CLASS_NAMES`, fix class
`6`, and correct the two wrong class names. Two clearly separated fetchers, both
writing cached TSVs into the gitignored `cache/data/tcdb/raw/`, both gated behind
`--refetch-raw` so `--force` stays a ~1 min no-network rebuild.

Scoping the scrape to only the ~596 kept subfamilies was rejected: it would make
the committed `tcdb_hierarchy.json` gene-set-dependent, so adding a strain would
silently reintroduce unnamed nodes until someone re-scraped.

Nothing below this line is settled.

## 5. Open design questions

**5a. Committed-artifact size.** 19,743 specificity names at ~150 chars each adds
roughly 3 MB to a `tcdb_hierarchy.json` that is currently 4.5 MB — nearly
doubling a committed, indented, sort-keyed file, for names of which only 286 are
ever used. Options: inline them anyway; put names in a separate committed
`tcdb_names.json`; or carry names only for the levels above `tc_specificity` in
the hierarchy and resolve specificity names at prune time. Decide before writing
the parser, because it fixes the file layout.

**5b. Shape of a `tc_specificity` name.** The curated text is a full paragraph
with inline citations. A node `name` feeding a full-text index and an explorer
label probably wants something short. Options: whole string; first sentence;
truncate at N chars; or short synthesized `name` plus the full text in a
separate `description` property. Related: strip `(Author et al., YYYY)` citation
parentheticals, or keep them?

**5c. Fallback when upstream has no name.** Unknown how many of the 3,550
subfamilies are genuinely blank upstream (6d). Current behaviour — render the
bare TC ID — stays the default unless something better is chosen. A synthesized
label such as `"2.A.1.1 (subfamily of The Major Facilitator Superfamily (MFS))"`
is an option but arguably worse than an honest bare ID.

**5d. Scrape hygiene.** Rate limit, retry policy, resume-after-partial-failure,
and what happens when a single page 404s or the site goes back into maintenance
mid-run. Given TCDB's demonstrated fragility this needs to be non-fatal and
resumable, and the guard added in `dc76076c` (step 6 raises if a rebuild comes
out empty beside a populated committed artifact) is the right precedent to
extend to the name maps.

## 6. To re-check against the live site before writing code

**6a. Does the 2015 subfamily-page template still hold?** *(highest risk)* Fetch
`https://www.tcdb.org/search/result.php?tc=2.A.1.1` and confirm the
`table id="result-cluster"` structure, the `<strong>…4-part:&nbsp;&nbsp;name</strong>`
section-header row, and the four-column system rows. If the template changed,
3b's parser design is void and this section restarts.

**6b. Is `3.F` on the current `browse.php`?** It is in our May-2026 TSVs but
absent from the 2025-07 snapshot. Confirms whether browse.php gives 52/52
subclass coverage or 51/52 with one permanent gap.

**6c. Probe for a cheaper bulk path.** In order of value:
`tcfamilybrowse.php?tc=1` and `?tc=1.A` (a class-level printer view would reduce
the whole job to 8 requests); `tcfamilybrowse.php?tc=<4-part>`; a re-read of the
live `/public/` index for files newer than the 2024-08-27 batch; and the live
`/download.php` for endpoints added since the archived copy.

**6d. Blank-name rate.** Over a sample of subfamily pages, how many subfamilies
have no upstream name — this answers 5c.

**6e. Politeness.** Read `robots.txt` and any stated rate limit, and size the
delay for ~3,550 requests accordingly. Measure real page sizes to replace the
~400 MB estimate.

**6f. Re-download the four existing TSVs.** They are cached from 2026-05-05
(`families`, `superfamilies`, `substrates`, `acc2tcid`) and 2026-08-06
(`tcdb_go_map`, `tcdb_pfam_map`). Refreshing them in the same pass keeps the
hierarchy and the new name maps from the same upstream release.

## 7. Constraints carried from the task brief

- Follow the CLAUDE.md 5-phase workflow (Scope / Implement / Test / Review /
  Document); no phase skipped.
- **Sanitise every string property** — `.replace("'", "^").replace("|", "")`.
  TCDB names carry apostrophes and Greek letters (α, β, Δ, μ) as a matter of
  course, and `|` is the array delimiter; an unescaped one silently splits
  values. That exact bug was just fixed for UniProt in `c043241d` — do not
  reintroduce it. HTML entities and `<sup>`/`<em>` markup must be unwrapped
  before sanitising.
- BioCypher mishandles boolean properties; use categorical strings.
- `cache/data/tcdb/tcdb_hierarchy.json` and `tcdb_pruned.json` are **committed**
  artifacts. The `dc76076c` empty-rebuild guard is intentional; extend it.
- Iteration loop: `bash scripts/prepare_data.sh --steps 6 --force` (~1 min, no
  network). `--refetch-raw` only when pulling upstream.
- **Fix the stale module docstring** (`build_kegg_metabolism_xrefs.py` ~line 22).
  It claims step 6 downloads *"the 3 TCDB reference TSVs (tc_classes,
  tc_subclasses, families)"*. None of those three names match what the code
  reads, there are four files, and no subclass file exists at all.

## 8. Definition of done

- Named nodes at subclass and subfamily level, and at specificity level if 6a
  confirms the source; `t.name = t.tcdb_id` count drops substantially from 916.
- Full Docker rebuild with `output/import.report` still empty (zero dangling).
- `pytest -m "not slow and not kg"` and `pytest -m kg` both green.
- Names searchable through the `tcdbFamilyFullText` index, verified with a real
  full-text query — not by inspecting the property.

## 9. Documentation targets

The explorer has not consumed the current TCDB upgrade and it is unreleased, so
fold this into the **same** item rather than creating a new one:

- `docs/kg-changes/tcdb-two-source-upgrade.md` — see its §7 for the pattern
- `CHANGELOG.md` `[Unreleased]`
- the `TcdbFamily` bullet in `CLAUDE.md`

## 10. Related deferred work

**T4** — `TcdbFamily.gene_count` is a SUBTREE count (`*0..`) while
`InterproEntry.gene_count` is direct. A warning is in CLAUDE.md; not fixed. Out
of scope here, but both touch the same post-import block, so worth batching if
T6 lands first.
