# InterPro Multi-Ontology Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rebuild the InterProScan integration as a multi-ontology source: faceted per-member-DB calls.json, a new NCBIfam ontology, GO/Pfam/EC/CAZy feeds through the existing provenance machinery, and naming recovery — per the approved spec `docs/superpowers/specs/2026-08-17-interpro-multi-ontology-redesign-design.md`.

**Architecture:** Phase-1 re-parse (cached raw.json → faceted calls.json, no re-scan) → central reference builds (InterPro descriptions extension + new NCBIfam reference) → step-2 merge (declarative fields + `enrich_interpro_fields` rewrite) → KG adapters (rewritten `interpro_adapter`, new `ncbifam_adapter`) → post-import (rollups, 9th bucket, `has_any_edge` fix) → Docker rebuild + validation.

**Tech Stack:** Python 3 (uv), pytest, BioCypher, Neo4j (Docker), Cypher post-import, YAML-driven merge config.

## Global Constraints

(Every task's requirements implicitly include these — copied from the spec.)

- **String sanitization**: every string property yielded by an adapter goes through `_clean_str(value) -> value.replace("'", "^").replace("|", "")`. Curated descriptions WILL contain apostrophes.
- **No e-value cutoff anywhere** — member DBs pre-apply curated thresholds; evalue/score are evidence-only, never filters.
- **Count-don't-combine**: no cross-library `score` on the InterPro rollup or `Gene_has_interpro_entry` edge. Rollup carries `evalue` (min, nullable) + `evalue_library`. `Gene_has_ncbifam_family` keeps both `evalue` + `score` (single library, homogeneous HMMER bits).
- **NCBIfam node IDs are the underscore form** `ncbifam_TIGR00198` / `ncbifam_NF002735` (`ncbifam` is NOT a bioregistry prefix; `normalize_curie` falls back — psortb/signalp precedent).
- **Merged-JSON authority**: adapters emit edges ONLY for entries in the merged seed fields (`interpro_entries`, `ncbifam_ids`); calls.json only decorates evidence. Skew fails soft (missing decorations, never phantom edges).
- **Direct attribution**: gene Pfam/NCBIfam contributions come only from direct hits (`libraries.PFAM` / `libraries.NCBIFAM`); entry-mediated GO/EC/CAZy are gated (GO: FAMILY+DOMAIN; EC: single-EC FAMILY; CAZy: FAMILY+DOMAIN) and evidence-labeled.
- **Declarative-first**: YAML + transforms where the machinery reaches; bespoke Python only as named post-merge enrichment functions (`enrich_*` precedent).
- **References committed, raw gitignored**: `cache/data/ncbifam/raw/` and `cache/data/interpro/raw/` are gitignored; the reference JSONs are committed.
- **Description cap**: InterPro abstracts truncated to 400 chars, plain-texted.
- Test commands: fast suite `uv run pytest -m "not slow and not kg" -q`; KG suite `uv run pytest -m kg -v` (needs Docker graph).
- Commit messages end with `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.

## File Structure (locked decomposition)

| File | Responsibility |
|---|---|
| `multiomics_kg/utils/interproscan.py` | REWRITE — pure parser: raw JSON → faceted calls dict + QC summary. No filesystem. |
| `.claude/skills/interproscan-run/run_interproscan.py` | MODIFY — add `--normalize` (re-parse from raw.json, no Docker); stop writing entry_xrefs.json. |
| `multiomics_kg/utils/ncbifam.py` | CREATE — pure parser for `hmm_PGAP.tsv` rows → reference entries; family_type vocab constants. |
| `multiomics_kg/download/build_ncbifam_reference.py` | CREATE — download hmm_PGAP.tsv → `cache/data/ncbifam/ncbifam_reference.json`. |
| `multiomics_kg/utils/interpro_reference.py` + `multiomics_kg/download/build_interpro_reference.py` | MODIFY — add `description` extraction (abstract, 400-char cap). |
| `multiomics_kg/download/build_gene_annotations.py` | MODIFY — `load_interproscan()` for new format; rewrite `enrich_interpro_fields()` (gates + naming recovery + gene_name fill). |
| `config/gene_annotations_config.yaml` | MODIFY — `ncbifam_ids` field; DataSource `info_types`. |
| `multiomics_kg/adapters/interpro_adapter.py` | REWRITE internals — new calls format, node `description`, edge `evalue_library` (no `score`). |
| `multiomics_kg/adapters/ncbifam_adapter.py` | CREATE — NcbifamFamily nodes + Gene_has_ncbifam_family + Ncbifam_family_in_interpro_entry. |
| `config/schema_config.yaml` | MODIFY — NcbifamFamily node + 2 edges; Gene_has_interpro_entry property change. |
| `create_knowledge_graph.py` | MODIFY — wire ncbifam adapter + injections. |
| `scripts/post-import.cypher` + `scripts/post-import.sh` | MODIFY — indexes, rollups, 9th bucket, has_any_edge, F1.1 rules, ncbifam_family_count. Both files identical logic. |
| `config/uninformative_terms.yaml` | MODIFY — `interpro_entry` (name_patterns) + `ncbifam` (family_types — NEW rule kind) sections. |
| `scripts/prepare_data.sh` | MODIFY — reference steps run before step 2 (reorder default STEPS; step 9 becomes "central references", gains ncbifam). |
| Tests | `tests/test_interproscan.py` (rewrite), `tests/test_ncbifam.py` (new), `tests/test_build_gene_annotations.py` (extend), `tests/test_interpro_adapter.py` (rewrite), `tests/test_ncbifam_adapter.py` (new), `tests/test_interproscan_consistency.py` (new), `tests/kg_validity/` (extend). |
| Docs | `docs/kg-changes/interpro-multi-ontology.md` (new), superseded headers on 2 old docs, CLAUDE.md, interproscan-run SKILL.md. |

**Execution order note:** Tasks 1–4 (artifact layer) → 5–7 (references) → 8–11 (merge) → 12–15 (adapters/schema) → 16–17 (post-import) → 18–20 (pipeline runs + validation) → 21 (docs). Tasks 5–6 can run parallel to 1–4.

---

### Task 1: Rewrite the pure parser (`multiomics_kg/utils/interproscan.py`)

**Files:**
- Modify: `multiomics_kg/utils/interproscan.py` (full rewrite of parse path; keep module docstring updated)
- Test: `tests/test_interproscan.py` (rewrite)

**Interfaces:**
- Consumes: raw InterProScan `--formats JSON` document (`{"results": [...]}`).
- Produces (later tasks rely on these exact shapes):
  - `parse_interproscan_json(data: dict) -> dict[str, dict]` — WP_-keyed faceted records:
    ```
    {
      "md5": str, "match_count": int,
      "libraries": { "<LIB>": [ {"accession": str, "name": str|None, "ipr": str|None,
                                  "start": int|None, "end": int|None,
                                  "evalue": float|None, "score": float|None}, ... ] },
      "interpro_entries": { "<IPR>": {"type": str, "libraries": [str,...], "match_count": int,
                                       "start": int|None, "end": int|None,
                                       "evalue": float|None, "evalue_library": str|None} },
      "go_terms": { "<GO:NNNNNNN>": ["<IPR>", ...] }
    }
    ```
    Facets sparse (a `libraries` key exists only if that DB matched; `interpro_entries`/`go_terms` may be `{}`). NO `pathways` anywhere. Signature accessions version-stripped (`NF002735.2` → `NF002735`). Zero-match proteins kept with `match_count: 0`, `libraries: {}`.
  - `summarize(calls, *, strain, input_proteins, tool_version, applications, image_digest=None, wallclock_s=None, parse_failures=0, xrefs_requested=None) -> dict` — same QC fields as today EXCEPT: drop `proteins_with_pathways`, `distinct_pathways`, `pathway_databases`; `distribution` counts per-library match rows from the facets; `interpro_integrated_matches` = rows with non-null `ipr`.
  - DELETED: `parse_entry_xrefs`, `normalize_pathway_xref`, `_extract_matches`, `_aggregate` (helpers replaced).

- [ ] **Step 1: Write the failing tests** — replace `tests/test_interproscan.py` content with tests against a small inline raw-format fixture:

```python
"""Tests for the faceted InterProScan parser (multi-ontology redesign)."""
import pytest
from multiomics_kg.utils.interproscan import parse_interproscan_json, summarize


def _loc(start, end, evalue=None, score=None):
    return {"start": start, "end": end, "evalue": evalue, "score": score}


def _match(library, acc, desc, entry=None, locations=None):
    sig = {"accession": acc, "description": desc,
           "signatureLibraryRelease": {"library": library}}
    if entry:
        sig["entry"] = entry
    return {"signature": sig, "locations": locations or [_loc(1, 50)]}


ENTRY_FAM = {
    "accession": "IPR003686", "description": "Photosystem II PsbI", "type": "FAMILY",
    "goXRefs": [{"id": "GO:0015979"}, {"id": "GO:0009523"}],
    "pathwayXRefs": [{"databaseName": "MetaCyc", "id": "PWY-101"}],  # must be DROPPED
}
ENTRY_SF = {
    "accession": "IPR037271", "description": "PsbI superfamily",
    "type": "HOMOLOGOUS_SUPERFAMILY", "goXRefs": [{"id": "GO:0015979"}],
    "pathwayXRefs": [],
}

RAW = {"results": [{
    "md5": "abc", "xref": [{"id": "WP_000001.1"}],
    "matches": [
        _match("PFAM", "PF02532.18", "PSII PsbI", ENTRY_FAM,
               [_loc(1, 36, evalue=4.1e-18, score=76.3)]),
        _match("HAMAP", "MF_01316", "PSII reaction center I", ENTRY_FAM,
               [_loc(1, 36, score=17.4)]),
        _match("NCBIFAM", "NF002735.2", "photosystem II protein I", None,
               [_loc(1, 38, evalue=3.3e-23, score=92.7)]),
        _match("SUPERFAMILY", "SSF161041", "PsbI", ENTRY_SF, [_loc(1, 35)]),
    ],
}, {
    "md5": "def", "xref": [{"id": "WP_000002.1"}], "matches": [],
}]}


@pytest.fixture()
def calls():
    return parse_interproscan_json(RAW)


def test_libraries_facet_sparse_and_version_stripped(calls):
    rec = calls["WP_000001.1"]
    assert set(rec["libraries"]) == {"PFAM", "HAMAP", "NCBIFAM", "SUPERFAMILY"}
    pf = rec["libraries"]["PFAM"][0]
    assert pf["accession"] == "PF02532"          # version stripped
    assert pf["ipr"] == "IPR003686"
    assert pf["evalue"] == 4.1e-18 and pf["score"] == 76.3
    nf = rec["libraries"]["NCBIFAM"][0]
    assert nf["accession"] == "NF002735" and nf["ipr"] is None


def test_interpro_rollup_no_score_evalue_attributed(calls):
    ent = calls["WP_000001.1"]["interpro_entries"]["IPR003686"]
    assert ent["type"] == "FAMILY"
    assert ent["libraries"] == ["HAMAP", "PFAM"]
    assert ent["match_count"] == 2
    assert ent["evalue"] == 4.1e-18 and ent["evalue_library"] == "PFAM"
    assert "score" not in ent                     # count-don't-combine
    sf = calls["WP_000001.1"]["interpro_entries"]["IPR037271"]
    assert sf["evalue"] is None and sf["evalue_library"] is None


def test_go_terms_carry_entry_attribution(calls):
    go = calls["WP_000001.1"]["go_terms"]
    assert go["GO:0015979"] == ["IPR003686", "IPR037271"]
    assert go["GO:0009523"] == ["IPR003686"]


def test_no_pathways_anywhere(calls):
    rec = calls["WP_000001.1"]
    assert "pathways" not in rec
    assert all("pathways" not in e for e in rec["interpro_entries"].values())


def test_zero_match_protein_kept(calls):
    rec = calls["WP_000002.1"]
    assert rec["match_count"] == 0 and rec["libraries"] == {} \
        and rec["interpro_entries"] == {} and rec["go_terms"] == {}


def test_summarize_qc():
    calls = parse_interproscan_json(RAW)
    s = summarize(calls, strain="X", input_proteins=2,
                  tool_version="5.78-109.0", applications="ALL_DEFAULT")
    assert s["calls_made"] == 2 and s["proteins_no_match"] == 1
    assert s["total_matches"] == 4
    assert s["interpro_integrated_matches"] == 3
    assert s["distribution"] == {"HAMAP": 1, "NCBIFAM": 1, "PFAM": 1, "SUPERFAMILY": 1}
    assert s["proteins_with_go_terms"] == 1 and s["distinct_go_terms"] == 2
    assert "pathway_databases" not in s and "distinct_pathways" not in s
```

- [ ] **Step 2: Run tests, verify they fail** — `uv run pytest tests/test_interproscan.py -q` → FAIL (old shapes).

- [ ] **Step 3: Rewrite the parser.** Core implementation (module keeps its docstring style; delete `normalize_pathway_xref`, `parse_entry_xrefs`):

```python
def _strip_version(acc: str | None) -> str | None:
    return acc.split(".")[0] if acc else acc


def parse_interproscan_json(data: dict) -> dict[str, dict]:
    calls: dict[str, dict] = {}
    for result in data.get("results") or []:
        libraries: dict[str, list[dict]] = {}
        entries: dict[str, dict] = {}          # IPR -> rollup accumulator
        go_terms: dict[str, set[str]] = {}
        n_rows = 0
        for m in result.get("matches") or []:
            sig = m.get("signature") or {}
            library = (sig.get("signatureLibraryRelease") or {}).get("library")
            entry = sig.get("entry") or {}
            ipr = entry.get("accession")
            for loc in m.get("locations") or []:
                n_rows += 1
                row = {
                    "accession": _strip_version(sig.get("accession")),
                    "name": sig.get("description") or sig.get("name"),
                    "ipr": ipr,
                    "start": loc.get("start"), "end": loc.get("end"),
                    "evalue": loc.get("evalue"), "score": loc.get("score"),
                }
                if library:
                    libraries.setdefault(library, []).append(row)
                if ipr:
                    ent = entries.setdefault(ipr, {
                        "type": entry.get("type"), "libraries": set(),
                        "match_count": 0, "start": None, "end": None,
                        "evalue": None, "evalue_library": None,
                    })
                    ent["match_count"] += 1
                    if library:
                        ent["libraries"].add(library)
                    s, e = loc.get("start"), loc.get("end")
                    if s is not None and (ent["start"] is None or s < ent["start"]):
                        ent["start"] = s
                    if e is not None and (ent["end"] is None or e > ent["end"]):
                        ent["end"] = e
                    ev = loc.get("evalue")
                    if ev is not None and (ent["evalue"] is None or ev < ent["evalue"]):
                        ent["evalue"], ent["evalue_library"] = ev, library
                    for x in entry.get("goXRefs") or []:
                        if x.get("id"):
                            go_terms.setdefault(x["id"], set()).add(ipr)
        for lib in libraries:
            libraries[lib].sort(key=lambda r: (r["start"] or 0, r["accession"] or ""))
        record = {
            "md5": result.get("md5"),
            "match_count": n_rows,
            "libraries": libraries,
            "interpro_entries": {
                k: {**v, "libraries": sorted(v["libraries"])}
                for k, v in sorted(entries.items())
            },
            "go_terms": {k: sorted(v) for k, v in sorted(go_terms.items())},
        }
        for xref in result.get("xref") or []:
            if xref.get("id"):
                calls[xref["id"]] = dict(record)
    return calls
```

`summarize` update: iterate `libraries` facets for `distribution` / `total_matches`; `interpro_integrated_matches = sum(1 for lib rows if row["ipr"])`; GO counters from the `go_terms` facet keys; delete the pathway counters.

- [ ] **Step 4: Run tests** — `uv run pytest tests/test_interproscan.py -q` → PASS.
- [ ] **Step 5: Fix collateral fast-suite breakage in THIS module only** — `uv run pytest -m "not slow and not kg" -q 2>&1 | head -40`. Expected failures at this point: `tests/test_build_gene_annotations.py` (load_interproscan) and `tests/test_interpro_adapter.py` — these are Tasks 8/12; do NOT fix here. Confirm no OTHER module broke.
- [ ] **Step 6: Commit** — `git add multiomics_kg/utils/interproscan.py tests/test_interproscan.py && git commit -m "feat(interproscan): faceted pure parser — per-library facets, attributed GO, no pathways"`

---

### Task 2: `--normalize` mode in the runner + stop writing entry_xrefs

**Files:**
- Modify: `.claude/skills/interproscan-run/run_interproscan.py`

**Interfaces:**
- Consumes: Task 1's `parse_interproscan_json` + `summarize`.
- Produces: `--normalize` CLI flag — for each selected strain with `<strain>.interproscan.raw.json`, re-parse → overwrite `<strain>.interproscan.calls.json` + `<strain>.interproscan.skill_summary.json`; delete `<strain>.interproscan.entry_xrefs.json` if present; no Docker. Status row per strain (`NORMALIZED` / `NO_RAW` / `FAILED`).

- [ ] **Step 1: Add the flag + handler.** In `main()` add `p.add_argument("--normalize", action="store_true", help="Re-parse cached raw.json into calls.json (new faceted format); no Docker, no scan")`. Add a `normalize_strain(strain, data_dir) -> tuple[str, str]` function next to `run_strain`:

```python
def normalize_strain(strain: str, data_dir: Path) -> tuple[str, str]:
    out_dir = data_dir / "interproscan"
    raw_json = out_dir / f"{strain}.interproscan.raw.json"
    if not raw_json.exists():
        return "NO_RAW", "raw.json missing — re-run the scan or use committed calls.json"
    t0 = time.time()
    with raw_json.open() as f:
        data = json.load(f)
    calls = parse_interproscan_json(data)
    n_prot = len(calls)
    summary = summarize(
        calls, strain=strain, input_proteins=n_prot,
        tool_version=IPS_VERSION, applications="ALL_DEFAULT",
        wallclock_s=time.time() - t0,
    )
    (out_dir / f"{strain}.interproscan.calls.json").write_text(
        json.dumps(calls, indent=1, sort_keys=True) + "\n")
    (out_dir / f"{strain}.interproscan.skill_summary.json").write_text(
        json.dumps(summary, indent=1, sort_keys=True) + "\n")
    stale = out_dir / f"{strain}.interproscan.entry_xrefs.json"
    if stale.exists():
        stale.unlink()
    return "NORMALIZED", f"{n_prot} proteins"
```

Wire in `main()`: when `args.normalize`, loop `load_genome_rows(strains=args.strains)` calling `normalize_strain`, print the status table, and exit before any Docker path. Also delete the `entry_xrefs` write in `run_strain` (the block around line 234-239) so future SCANS also produce only the new artifacts.
- [ ] **Step 2: Smoke on MED4** — `uv run python .claude/skills/interproscan-run/run_interproscan.py --normalize --strains MED4`. Expected: `MED4 | NORMALIZED | 1872 proteins`; `MED4.interproscan.entry_xrefs.json` gone.
- [ ] **Step 3: Verify MED4 spot checks (spec §6 Layer 1)**:

```bash
jq '."WP_002805124.1".libraries.PFAM[0] | {accession, ipr}' cache/data/Prochlorococcus/genomes/MED4/interproscan/MED4.interproscan.calls.json
# {"accession": "PF02532", "ipr": "IPR003686"}
jq '."WP_002805124.1".libraries.NCBIFAM[0].ipr' …/MED4.interproscan.calls.json   # null
jq '."WP_002805124.1".go_terms."GO:0015979"' …/MED4.interproscan.calls.json      # ["IPR003686", ...]
jq '[.[] | .libraries.NCBIFAM[]?.accession] | index("TIGR00198")' …/MED4.interproscan.calls.json  # null (absent)
```

- [ ] **Step 4: Commit runner change only** (MED4 artifacts come with the full batch) — `git add .claude/skills/interproscan-run/run_interproscan.py && git commit -m "feat(interproscan-run): --normalize re-parse mode; drop entry_xrefs sidecar"` (checkout MED4 artifact changes first: `git checkout -- cache/`).

---

### Task 3: Normalize all 42 strains, commit artifacts

**Files:**
- Modify (generated): `cache/data/*/genomes/*/interproscan/<strain>.interproscan.{calls,skill_summary}.json` (42×2)
- Delete: `cache/data/*/genomes/*/interproscan/<strain>.interproscan.entry_xrefs.json` (42)

- [ ] **Step 1: Run the batch** — `uv run python .claude/skills/interproscan-run/run_interproscan.py --normalize 2>&1 | tail -50`. Expected: 42 `NORMALIZED` rows, zero `NO_RAW`/`FAILED`.
- [ ] **Step 2: Cross-strain spot checks** (spec §6): EZ55 `WP_156086936.1` + MIT1002 `WP_014977393.1` have `TIGR00198` under `libraries.NCBIFAM[].accession`; MED4 whole-file has none; MED4 `WP_011131653.1` has `TIGR00357`.
- [ ] **Step 3: Corpus sanity** — one-off python: count distinct NCBIFAM accessions across all new calls.json = **4,957**; total NCBIFAM rows ≈ 67,918; deduped (protein, accession) pairs ≈ 67,078.
- [ ] **Step 4: Verify entry_xrefs gone** — `ls cache/data/*/genomes/*/interproscan/*entry_xrefs* 2>/dev/null | wc -l` → 0.
- [ ] **Step 5: Commit** — `git add -A cache/data && git commit -m "data(interproscan): re-normalize all 42 strains to faceted calls.json; drop entry_xrefs sidecars"`

---

### Task 4: Update interproscan-run SKILL.md

**Files:**
- Modify: `.claude/skills/interproscan-run/SKILL.md`

- [ ] **Step 1:** Rewrite the **Output Schema** section to the faceted format (copy the record shape from Task 1's Interfaces block, field-by-field, with the sparse/`ipr: null`/version-stripped/no-pathways rules). Add `--normalize` to Quick Start (`uv run python .claude/skills/interproscan-run/run_interproscan.py --normalize [--strains X]`) with a note: normalization changes require raw.json; machines without it use committed calls.json. Replace the spot-check table with the four Layer-1 checks from Task 2 Step 3 + Task 3 Step 2 (exact jq one-liners + expected values). Note entry_xrefs.json no longer exists (central `interpro_reference.json` replaces it).
- [ ] **Step 2: Commit** — `git add .claude/skills/interproscan-run/SKILL.md && git commit -m "docs(interproscan-run): faceted output schema, --normalize workflow, refreshed spot checks"`

---

### Task 5: NCBIfam reference — parser util + builder + gitignore

**Files:**
- Create: `multiomics_kg/utils/ncbifam.py`
- Create: `multiomics_kg/download/build_ncbifam_reference.py`
- Modify: `.gitignore` (add `cache/data/ncbifam/raw/`)
- Test: `tests/test_ncbifam.py`

**Interfaces:**
- Produces:
  - `multiomics_kg/utils/ncbifam.py`:
    - `HYPOTH_FAMILY_TYPES = frozenset({"hypoth_equivalog", "hypoth_equivalog_domain"})`
    - `parse_hmm_pgap_rows(rows: Iterable[dict]) -> dict[str, dict]` — DictReader rows → `{unversioned_acc: {"name": product_name, "family_type": str, "gene_symbol": str|None, "description": str|None, "gene_synonyms": [str], "pmids": [str], "ec_numbers": [str], "go_terms": [str]}}`. Empty-string TSV cells → key omitted (sparse; lists omitted when empty). Column names (verified 2026-08-17): `#ncbi_accession, source_identifier, label, sequence_cutoff, domain_cutoff, hmm_length, family_type, for_structural_annotation, for_naming, for_AMRFinder, product_name, gene_symbol, gene_synonyms, ec_numbers, go_terms, pmids, taxonomic_range, taxonomic_range_name, taxonomic_rank_name, n_refseq_protein_hits, source, name_orig, hmm_name, comment`. `description` ← `comment`; multi-valued cells split on `,` (ec_numbers/go_terms/pmids) or space (gene_synonyms — verify empirically at implementation, adjust the split, and pin with a test).
  - `build_ncbifam_reference.py` — CLI (`--force`, `--refetch-raw`) mirroring `build_interpro_reference.py`: download `https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv` → `cache/data/ncbifam/raw/hmm_PGAP.tsv` (skip when present unless `--refetch-raw`), parse, write `cache/data/ncbifam/ncbifam_reference.json` (indented, sort_keys). Runnable: `uv run python -m multiomics_kg.download.build_ncbifam_reference`.

- [ ] **Step 1: Write failing tests** in `tests/test_ncbifam.py`:

```python
from multiomics_kg.utils.ncbifam import parse_hmm_pgap_rows, HYPOTH_FAMILY_TYPES


def _row(**kw):
    base = {"#ncbi_accession": "NF000001.1", "family_type": "equivalog",
            "product_name": "widget synthase", "gene_symbol": "", "comment": "",
            "gene_synonyms": "", "ec_numbers": "", "go_terms": "", "pmids": ""}
    base.update(kw)
    return base


def test_parse_basic_and_version_strip():
    ref = parse_hmm_pgap_rows([_row()])
    assert ref["NF000001"]["name"] == "widget synthase"
    assert ref["NF000001"]["family_type"] == "equivalog"
    assert "gene_symbol" not in ref["NF000001"]        # sparse
    assert "description" not in ref["NF000001"]


def test_parse_rich_row():
    ref = parse_hmm_pgap_rows([_row(**{"#ncbi_accession": "TIGR00198.1",
        "gene_symbol": "katG", "comment": "catalase-peroxidase HPI",
        "ec_numbers": "1.11.1.21", "pmids": "9006042,9871101"})])
    e = ref["TIGR00198"]
    assert e["gene_symbol"] == "katG"
    assert e["description"] == "catalase-peroxidase HPI"
    assert e["ec_numbers"] == ["1.11.1.21"] and e["pmids"] == ["9006042", "9871101"]


def test_hypoth_types_constant():
    assert "hypoth_equivalog" in HYPOTH_FAMILY_TYPES
    assert "hypoth_equivalog_domain" in HYPOTH_FAMILY_TYPES
```

- [ ] **Step 2: Run, verify FAIL** — `uv run pytest tests/test_ncbifam.py -q`.
- [ ] **Step 3: Implement** `multiomics_kg/utils/ncbifam.py` (pure) and `build_ncbifam_reference.py` (download + parse + write; copy the CLI/paths/logging skeleton from `build_interpro_reference.py`, swapping the fetch for a single-file download via `urllib.request.urlretrieve`).
- [ ] **Step 4: Run tests** → PASS.
- [ ] **Step 5: Build the real reference** — `uv run python -m multiomics_kg.download.build_ncbifam_reference`. Verify: `python3 -c "import json; d=json.load(open('cache/data/ncbifam/ncbifam_reference.json')); print(len(d)); print(d['TIGR00198']['gene_symbol'])"` → ~38,394 entries / `katG`.
- [ ] **Step 6: gitignore + commit** — add `cache/data/ncbifam/raw/` to `.gitignore`; `git add multiomics_kg/utils/ncbifam.py multiomics_kg/download/build_ncbifam_reference.py tests/test_ncbifam.py .gitignore cache/data/ncbifam/ncbifam_reference.json && git commit -m "feat(ncbifam): reference builder + committed ncbifam_reference.json from hmm_PGAP.tsv"`

---

### Task 6: InterPro reference gains `description`

**Files:**
- Modify: `multiomics_kg/utils/interpro_reference.py` (abstract extraction in the XML parse), `multiomics_kg/download/build_interpro_reference.py` (plumb through)
- Test: `tests/test_interpro_reference.py` (extend)

**Interfaces:**
- Produces: `interpro_reference.json` entries gain sparse `description: str` — first paragraph of the entry's `<abstract>`, HTML tags stripped, whitespace collapsed, truncated to 400 chars. Helper `clean_abstract(html_text: str, cap: int = 400) -> str` in `interpro_reference.py`.

- [ ] **Step 1: Failing test** (add to `tests/test_interpro_reference.py`):

```python
from multiomics_kg.utils.interpro_reference import clean_abstract

def test_clean_abstract_strips_tags_and_caps():
    html = "<p>This family represents <i>PsbI</i>, a small subunit.</p><p>Second para.</p>"
    out = clean_abstract(html)
    assert out == "This family represents PsbI, a small subunit."
    assert len(clean_abstract("<p>" + "x" * 1000 + "</p>")) == 400
```

- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement** — `clean_abstract`: take text up to the first `</p>` (or whole text if none), `re.sub(r"<[^>]+>", "", …)`, unescape entities (`html.unescape`), collapse whitespace, `[:400]`. In the streaming `interpro.xml.gz` pass (same place `db_xrefs` are parsed), capture the first `<abstract>` block per entry and store `description` when non-empty.
- [ ] **Step 4: Run tests → PASS.**
- [ ] **Step 5: Rebuild the reference** — raw XML must exist at `cache/data/interpro/raw/interpro.xml.gz` (it is gitignored; if missing run with `--refetch-raw`): `uv run python -m multiomics_kg.download.build_interpro_reference --force`. Verify size (`du -h cache/data/interpro/interpro_reference.json` — expect ≤ ~25 MB; if larger, apply the spec's documented fallback: keep descriptions only for entries observed in any committed calls.json + their ancestors) and content (`IPR003686` description non-empty).
- [ ] **Step 6: Commit** — code + tests + regenerated reference. `git commit -m "feat(interpro-reference): entry descriptions (400-char abstracts)"`

---

### Task 7: prepare_data ordering — references before merge

**Files:**
- Modify: `scripts/prepare_data.sh`

**Interfaces:** step 9 becomes "central reference builds" and runs BOTH `build_interpro_reference` and `build_ncbifam_reference`; the default STEPS order runs it before step 1/2. (Step keeps its number — flags/logs stay stable; ORDER is what the spec makes normative.)

- [ ] **Step 1:** Change `STEPS="0 1 2 3 4 5 6 7 8 9"` → `STEPS="0 9 1 2 3 4 5 6 7 8"` with a comment: `# 9 = central references (interpro + ncbifam); runs BEFORE the step-2 merge which consumes them (spec 2026-08-17 §4)`. In the `9)` case block, append a second `run_step 9` invocation (or extend the same log) for `uv run python -m multiomics_kg.download.build_ncbifam_reference` alongside the existing interpro one, honoring `--force`.
- [ ] **Step 2:** Sanity: `bash scripts/prepare_data.sh --steps 9` (no `--force` — both references exist, should skip/fast-noop). Check `logs/prepare_data_step9.log`.
- [ ] **Step 3: Commit** — `git commit -m "feat(prepare-data): step 9 = central references (interpro+ncbifam), ordered before merge"`

---

### Task 8: `load_interproscan()` for the faceted format

**Files:**
- Modify: `multiomics_kg/download/build_gene_annotations.py` (function `load_interproscan`, ~line 558)
- Test: `tests/test_build_gene_annotations.py` (extend)

**Interfaces:**
- Produces per-protein row for the merge (consumed by YAML fields + Task 9's enrich): `{"interpro_entries": [IPR,...], "pfam_signatures": [PF*,...], "ncbifam_ids": [acc,...], "hamap_descriptions": [str,...], "go_term_donors": {GO: [IPR,...]}}` — every key sparse (omitted when empty); row dropped when all empty. Accessions already version-stripped by the parser.

- [ ] **Step 1: Failing test** (add to `tests/test_build_gene_annotations.py`; write the new-format calls.json into `tmp_path`):

```python
def test_load_interproscan_faceted(tmp_path):
    from multiomics_kg.download.build_gene_annotations import load_interproscan
    d = tmp_path / "interproscan"; d.mkdir()
    calls = {"WP_1.1": {
        "md5": "x", "match_count": 3,
        "libraries": {
            "PFAM": [{"accession": "PF02532", "name": "PSII PsbI", "ipr": "IPR003686",
                      "start": 1, "end": 36, "evalue": 4.1e-18, "score": 76.3}],
            "NCBIFAM": [{"accession": "NF002735", "name": "psbI", "ipr": None,
                         "start": 1, "end": 38, "evalue": 3.3e-23, "score": 92.7}],
            "HAMAP": [{"accession": "MF_01316", "name": "Photosystem II reaction center protein I",
                       "ipr": "IPR003686", "start": 1, "end": 36, "evalue": None, "score": 17.4}],
        },
        "interpro_entries": {"IPR003686": {"type": "FAMILY", "libraries": ["HAMAP", "PFAM"],
                                           "match_count": 2, "start": 1, "end": 36,
                                           "evalue": 4.1e-18, "evalue_library": "PFAM"}},
        "go_terms": {"GO:0015979": ["IPR003686"]},
    }}
    (d / "S.interproscan.calls.json").write_text(json.dumps(calls))
    rows = load_interproscan(str(tmp_path), "S")
    r = rows["WP_1.1"]
    assert r["interpro_entries"] == ["IPR003686"]
    assert r["pfam_signatures"] == ["PF02532"]
    assert r["ncbifam_ids"] == ["NF002735"]
    assert r["hamap_descriptions"] == ["Photosystem II reaction center protein I"]
    assert r["go_term_donors"] == {"GO:0015979": ["IPR003686"]}
```

- [ ] **Step 2: Run → FAIL** (old format expectations).
- [ ] **Step 3: Implement** — rewrite the body: `entries = sorted(call.get("interpro_entries") or {})`; `libs = call.get("libraries") or {}`; `pfam_sigs = sorted({r["accession"] for r in libs.get("PFAM", []) if (r.get("accession") or "").startswith("PF")})`; `ncbifam = sorted({r["accession"] for r in libs.get("NCBIFAM", []) if r.get("accession")})`; `hamap = sorted({r["name"] for r in libs.get("HAMAP", []) if r.get("name")})`; `donors = call.get("go_terms") or {}`. Build the sparse row; skip when empty. Update the docstring (spec reference → 2026-08-17).
- [ ] **Step 4: Run → PASS.** Also fix any OLD-format tests for `load_interproscan` in this file (update fixtures to new format).
- [ ] **Step 5: Commit** — `git commit -m "feat(merge): load_interproscan surfaces faceted fields (ncbifam, hamap, GO donors)"`

---

### Task 9: Merge — YAML field + rewritten `enrich_interpro_fields` (gates, naming recovery, gene_name fill)

**Files:**
- Modify: `config/gene_annotations_config.yaml` (add `ncbifam_ids` field; update interproscan source comment + `info_types`)
- Modify: `multiomics_kg/download/build_gene_annotations.py` (`enrich_interpro_fields` rewrite; call-site passes `ncbifam_ref`; load ncbifam reference lazily like `interpro_ref`)
- Test: `tests/test_build_gene_annotations.py` (extend)

**Interfaces:**
- Consumes: Task 8 row shape; `interpro_reference.json` (types/xrefs/names); `ncbifam_reference.json` (product names, gene_symbols).
- Produces on the merged gene dict:
  - `ncbifam_ids: [acc,...]` (declarative passthrough from the interproscan source, `track_source: ncbifam_ids_source` — every token sourced `interproscan`)
  - `go_terms`/`ec_numbers`/`cazy_ids`/`pfam_ids` contributions exactly as the current Layer B, but GO gating now driven by `go_term_donors` attribution: a GO is added iff ≥1 donor entry's reference type ∈ {FAMILY, DOMAIN}; evidence `family` if any FAMILY donor else `domain`.
  - `alternate_functional_descriptions` += `[interpro] <entry name>` (FAMILY/DOMAIN, as today) + `[hamap] <desc>` + `[ncbifam] <product_name>` — each **skipped when it case-insensitively equals the gene's `product`** or is already present (dedup rule).
  - `gene_name` fill-if-empty from ncbifam `gene_symbol` (first symbol among the gene's `ncbifam_ids`, sorted for determinism); sets `gene_name_source = "ncbifam"`. Never overwrites.

- [ ] **Step 1: YAML** — in `gene_annotations_config.yaml` add after `interpro_entries` (~line 703):

```yaml
  ncbifam_ids:
    # Direct NCBIFAM HMM hits (curated prokaryotic families; TIGR* + NF*).
    # Surfaced from the faceted calls.json by load_interproscan(). Feeds the
    # NcbifamFamily ontology (adapter reads calls.json for edge evidence).
    type: union
    track_source: ncbifam_ids_source
    sources:
      - source: interproscan
        source_label: interproscan
        field: ncbifam_ids
```

Also update the interproscan `logical_sources` description/`info_types`-relevant comments (source block ~line 81) to list the new contributions (interpro_entries, ncbifam, pfam, go, ec, cazy, naming recovery). If `info_types` lives in this YAML for the DataSource node, update the list in place.
- [ ] **Step 2: Failing tests** (extend `tests/test_build_gene_annotations.py`; test `enrich_interpro_fields` directly):

```python
def _mk_gene(**kw):
    g = {"product": "photosystem II reaction center protein I",
         "interpro_entries": ["IPR003686"], "ncbifam_ids": ["NF002735"]}
    g.update(kw); return g

IPR_REF = {"IPR003686": {"name": "Photosystem II PsbI", "type": "FAMILY",
                         "go_terms": ["GO:0015979"], "ec_numbers": []},
           "IPR999999": {"name": "Some fold", "type": "HOMOLOGOUS_SUPERFAMILY",
                         "go_terms": ["GO:0000001"]}}
NCBIFAM_REF = {"NF002735": {"name": "photosystem II reaction center protein I",
                            "family_type": "equivalog", "gene_symbol": "psbI"}}

def test_go_gate_uses_donor_attribution():
    from multiomics_kg.download.build_gene_annotations import enrich_interpro_fields
    g = _mk_gene(interpro_entries=["IPR003686", "IPR999999"])
    row = {"go_term_donors": {"GO:0015979": ["IPR003686"], "GO:0000001": ["IPR999999"]}}
    enrich_interpro_fields(g, row, IPR_REF, NCBIFAM_REF)
    assert "GO:0015979" in g["go_terms"]
    assert "GO:0000001" not in g.get("go_terms", [])   # superfamily-only donor refused

def test_naming_recovery_dedup_against_product():
    from multiomics_kg.download.build_gene_annotations import enrich_interpro_fields
    g = _mk_gene()
    row = {"hamap_descriptions": ["Photosystem II reaction center protein I"],
           "ncbifam_ids": ["NF002735"]}
    enrich_interpro_fields(g, row, IPR_REF, NCBIFAM_REF)
    afd = g.get("alternate_functional_descriptions", [])
    # both tokens case-insensitively equal product -> both skipped
    assert not any(x.startswith("[hamap]") or x.startswith("[ncbifam]") for x in afd)

def test_gene_name_fill_if_empty_never_overwrites():
    from multiomics_kg.download.build_gene_annotations import enrich_interpro_fields
    g = _mk_gene()
    enrich_interpro_fields(g, {}, IPR_REF, NCBIFAM_REF)
    assert g["gene_name"] == "psbI" and g["gene_name_source"] == "ncbifam"
    g2 = _mk_gene(gene_name="psbI_existing", gene_name_source="uniprot")
    enrich_interpro_fields(g2, {}, IPR_REF, NCBIFAM_REF)
    assert g2["gene_name"] == "psbI_existing" and g2["gene_name_source"] == "uniprot"
```

- [ ] **Step 3: Run → FAIL** (new signature).
- [ ] **Step 4: Implement** — change signature to `enrich_interpro_fields(gene, ipr_row, interpro_ref, ncbifam_ref)`; GO block iterates `ipr_row.get("go_term_donors") or {}` (fallback: when the key is absent, derive from entries×reference exactly as today so old callers keep working); EC/CAZy/pfam blocks unchanged in logic; append the naming-recovery block:

```python
    product_lc = (gene.get("product") or "").strip().lower()
    def _maybe_desc(tag: str, text: str) -> None:
        text = (text or "").strip()
        if not text or text.lower() == product_lc:
            return
        desc_entries.append(f"[{tag}] {text}")

    for h in (ipr_row.get("hamap_descriptions") or []):
        _maybe_desc("hamap", h)
    symbols: list[str] = []
    for acc in (gene.get("ncbifam_ids") or []):
        meta = (ncbifam_ref or {}).get(acc) or {}
        _maybe_desc("ncbifam", meta.get("name") or "")
        if meta.get("gene_symbol"):
            symbols.append(meta["gene_symbol"])
    if symbols and not gene.get("gene_name"):
        gene["gene_name"] = sorted(symbols)[0]
        gene["gene_name_source"] = "ncbifam"
```

At the call site (~line 1297) load `ncbifam_ref` once (lazy, next to `interpro_ref`; from `cache/data/ncbifam/ncbifam_reference.json`, `{}` when missing) and pass it.
- [ ] **Step 5: Run tests → PASS**; run the full fast suite: `uv run pytest -m "not slow and not kg" -q` (only Task-12/13 adapter tests may still be red).
- [ ] **Step 6: Commit** — `git commit -m "feat(merge): ncbifam_ids field, donor-attributed GO gate, naming recovery + gene_name fill"`

---

### Task 10: Calls↔merge consistency test

**Files:**
- Create: `tests/test_interproscan_consistency.py`

**Interfaces:** Consumes committed calls.json (Task 3 format) + committed `gene_annotations_merged.json` (regenerated in Task 18). The test is the loud-failure guard for stale merges — it may legitimately FAIL between Task 3 and Task 18's merge re-run; mark accordingly.

- [ ] **Step 1: Write the test:**

```python
"""Calls.json <-> gene_annotations_merged.json consistency (spec 2026-08-17 §4.2/§6).

The two ungated passthrough fields must agree exactly per gene:
merged interpro_entries == calls interpro_entries keys, and
merged ncbifam_ids == calls libraries.NCBIFAM accessions.
A failure means calls.json was re-normalized without re-running
`prepare_data --steps 2` (or vice versa) — fix by re-running the merge.
"""
import csv, json
from pathlib import Path
import pytest

GENOMES_CSV = Path("data/Prochlorococcus/genomes/cyanobacteria_genomes.csv")


def _strains():
    with open(GENOMES_CSV, newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(l for l in fh if not l.lstrip().startswith("#"))]
    return [(r["strain_name"], Path(r["data_dir"])) for r in rows if r.get("data_dir")]


@pytest.mark.parametrize("strain,data_dir", _strains(), ids=lambda v: str(v))
def test_calls_merge_consistency(strain, data_dir):
    calls_p = data_dir / "interproscan" / f"{strain}.interproscan.calls.json"
    merged_p = data_dir / "gene_annotations_merged.json"
    if not calls_p.exists() or not merged_p.exists():
        pytest.skip("artifact missing")
    calls = json.loads(calls_p.read_text())
    merged = json.loads(merged_p.read_text())
    for locus_tag, gene in merged.items():
        wp = (gene.get("protein_id") or "").strip()
        call = calls.get(wp)
        if not call:
            assert not gene.get("interpro_entries"), f"{strain}:{locus_tag} merged has entries, calls lack protein"
            continue
        assert sorted(gene.get("interpro_entries") or []) == sorted(call.get("interpro_entries") or {}), \
            f"{strain}:{locus_tag} interpro_entries skew — re-run prepare_data --steps 2"
        calls_ncbifam = sorted({r["accession"] for r in (call.get("libraries") or {}).get("NCBIFAM", [])})
        assert sorted(gene.get("ncbifam_ids") or []) == calls_ncbifam, \
            f"{strain}:{locus_tag} ncbifam_ids skew — re-run prepare_data --steps 2"
```

(Adapt the genomes-CSV column names to the actual header — check `head -1` of the file; use the same columns `MultiInterproAnnotationAdapter._build_strain_adapters` reads: `data_dir`, and the strain name is `Path(data_dir).name`.)
- [ ] **Step 2: Run it** — `uv run pytest tests/test_interproscan_consistency.py -q`. At this point (before Task 18's merge) it is EXPECTED to fail/skew for strains whose merge predates the new fields; that's the test doing its job. If red: fine, note it; it must pass after Task 18.
- [ ] **Step 3: Commit** — `git commit -m "test(interproscan): calls<->merge consistency guard"`

---

### Task 11: uninformative_terms.yaml — new sections + third rule kind

**Files:**
- Modify: `config/uninformative_terms.yaml`

**Interfaces:** consumed by hand-written Cypher in Task 16 (the YAML documents the vocabulary; Cypher mirrors it — existing F1.1 pattern).

- [ ] **Step 1:** Append:

```yaml
interpro_entry:
  # Mirrored in post-import Cypher F1.1 (name-pattern rule, KEGG-KO precedent).
  name_patterns:
    - '^Protein of unknown function'
    - '^Domain of unknown function'
    - '^Uncharacterised protein family'
    - '^Uncharacterized protein family'

ncbifam:
  # THIRD RULE KIND (spec 2026-08-17 §5.5): property-valued — flags by
  # family_type, not id/name. hypoth_equivalog = "conserved family of unknown
  # function" (126 observed nodes, ~2.5%). Name-pattern fallback included.
  family_types:
    - hypoth_equivalog
    - hypoth_equivalog_domain
  name_patterns:
    - '(?i)hypothetical'
    - '(?i)uncharacterized'
    - '\bDUF\d'
```

- [ ] **Step 2: Commit** — `git commit -m "feat(uninformative-terms): interpro_entry + ncbifam sections (family_types rule kind)"`

---

### Task 12: Rewrite `interpro_adapter.py` for the faceted format

**Files:**
- Modify: `multiomics_kg/adapters/interpro_adapter.py`
- Test: `tests/test_interpro_adapter.py` (rewrite fixtures + assertions)

**Interfaces:**
- Consumes: faceted calls.json (Task 1 shape), merged JSON seeds, `interpro_reference.json` (now with `description`), injected `pfam_node_ids` / `ec_node_ids` (unchanged contract).
- Produces (changes only — everything not listed stays as-is, including hierarchy edges, Layer-A EC/CAZy routers, kept_ids pruning):
  - Nodes gain `description` (sparse: `_clean_str(ref.get("description"))` when non-empty).
  - `Gene_has_interpro_entry` edge props: `start`, `end`, `evalue` (nullable), **`evalue_library`** (str, present iff evalue is), `libraries` (str[]), `match_count`. **NO `score`.** Read DIRECTLY from the per-protein `interpro_entries` rollup (no more match grouping): the adapter iterates merged `interpro_entries` seeds (authority), looks up `calls[wp]["interpro_entries"][acc]` for evidence; missing rollup → edge with `{match_count: 0, libraries: []}` only (fail-soft skew).
  - `get_pfam_to_interpro()` now reads `calls[wp]["libraries"]["PFAM"]` rows (`accession`/`ipr`).
  - REMOVED: `aggregate_match_evidence` (rollup is precomputed at normalize time).

- [ ] **Step 1: Rewrite tests** — update `tests/test_interpro_adapter.py` fixtures to the faceted calls format; assert: edge props contain `evalue_library` and never `score`; node props contain `description` when the reference provides one; fail-soft edge for a merged seed missing from calls; pfam bridge from the PFAM facet; Layer-A behavior unchanged (keep those existing tests, fixtures updated).
- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement** per the Interfaces block. Edge emission core:

```python
    def get_edges(self):
        count = 0
        for locus_tag, gene in self._genes.items():
            protein_id = (gene.get("protein_id") or "").strip()
            entries = gene.get("interpro_entries") or []      # merged = authority
            if not protein_id or not entries:
                continue
            rollups = (self._calls.get(protein_id) or {}).get("interpro_entries") or {}
            for acc in entries:
                ent = rollups.get(acc)
                props: dict = {"match_count": 0, "libraries": []}
                if ent:
                    props = {"match_count": ent.get("match_count") or 0,
                             "libraries": [_clean_str(x) for x in ent.get("libraries") or []]}
                    for k in ("start", "end"):
                        if ent.get(k) is not None:
                            props[k] = ent[k]
                    if ent.get("evalue") is not None:
                        props["evalue"] = ent["evalue"]
                        props["evalue_library"] = _clean_str(ent.get("evalue_library"))
                yield (f"{locus_tag}-has_interpro-{acc}", _gene_node_id(locus_tag),
                       _interpro_node_id(acc), "gene_has_interpro_entry", props)
                count += 1
                if self.test_mode and count >= 100:
                    return
```

- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Commit** — `git commit -m "feat(interpro-adapter): faceted calls consumption, node descriptions, evalue_library (no score)"`

---

### Task 13: New `ncbifam_adapter.py`

**Files:**
- Create: `multiomics_kg/adapters/ncbifam_adapter.py`
- Test: `tests/test_ncbifam_adapter.py`

**Interfaces:**
- Consumes: merged JSON (`ncbifam_ids` seeds + `protein_id`), faceted calls.json (`libraries.NCBIFAM` evidence, `ipr` field for the bridge), `ncbifam_reference.json` (node metadata), injected `interpro_kept_ids: set[str] | None` (the InterPro adapter's kept node set — bridge target guarantee; `None` → no bridge edges).
- Produces:
  - `_ncbifam_node_id(acc) -> f"ncbifam_{acc}"` (underscore — Global Constraints).
  - `NcbifamAnnotationAdapter(genome_dir, test_mode)` — per-strain; `get_all_ncbifam_ids() -> set[str]`; `get_edges()` yields `(f"{locus_tag}-has_ncbifam-{acc}", gene_id, node_id, "gene_has_ncbifam_family", props)` — one edge per (gene, accession), evidence from the best NCBIFAM facet row for that accession (min evalue; carry that row's `score`, `start`, `end`; all sparse when null). Merged `ncbifam_ids` = authority; missing facet row → `{}` props (fail-soft).
  - `MultiNcbifamAdapter(genome_config_file, cache_root, interpro_kept_ids=None, test_mode=False)` — `download_data()` loads the reference; `get_nodes()` yields observed-only nodes `(node_id, "ncbifam family", {name, ncbifam_id, family_type, gene_symbol?, description?, level: 0})` (sparse props omitted; retired accessions absent from the reference → fallback `name` from any calls facet row's `name`, `family_type` omitted); `get_edges()` yields `Ncbifam_family_in_interpro_entry` bridge edges (`(f"{acc}-in_interpro-{ipr}", ncbifam_node_id, interpro_node_id, "ncbifam_family_in_interpro_entry", {})` for every observed (accession→ipr) pair from NCBIFAM facet rows with non-null `ipr`, pruned to `interpro_kept_ids`) + delegates per-strain gene edges.

- [ ] **Step 1: Write failing tests** — mirror `tests/test_interpro_adapter.py` structure: tmp_path genome dir with merged JSON + calls.json + a small reference; assert node shape (underscore id, family_type, description sparse), gene-edge evidence (evalue+score both present — single-library rule), one edge per (gene, acc) even with 2 facet rows, bridge pruning (`interpro_kept_ids=None` → no bridges; set → only members), retired-accession fallback node, `_clean_str` on name/description.
- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement** (copy the two-class skeleton from `interpro_adapter.py`; ~180 lines).
- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Commit** — `git commit -m "feat(ncbifam-adapter): NcbifamFamily nodes + gene edges + interpro bridge"`

---

### Task 14: Schema + wiring

**Files:**
- Modify: `config/schema_config.yaml`
- Modify: `create_knowledge_graph.py` (~line 312 block)

- [ ] **Step 1: Schema.** In `schema_config.yaml`:
  - Update `gene to interpro entry association` properties: delete `score: float`; add `evalue_library: str` (comment: attribution for the min-evalue; count-don't-combine — see spec §3.1).
  - Add to `interpro entry` properties: `description: str  # sparse — 400-char curated abstract`.
  - Add (next to the interpro block, following its comment style):

```yaml
ncbifam family:
  represented_as: node
  preferred_id: ncbifam
  label_in_input: ncbifam family
  is_a: biological entity
  properties:
    name: str               # curated product name (hmm_PGAP product_name)
    ncbifam_id: str         # e.g. "TIGR00198" / "NF002735" (unversioned)
    family_type: str        # equivalog | subfamily | domain | exception | hypoth_equivalog | ... (specificity LABEL, flat ontology)
    gene_symbol: str        # sparse
    description: str        # sparse — curated comment column
    level: int              # always 0 (flat)

gene to ncbifam family association:
  is_a: association
  represented_as: edge
  label_as_edge: Gene_has_ncbifam_family
  source: gene
  target: ncbifam family
  label_in_input: gene_has_ncbifam_family
  # Single-source scored edge (InterProScan NCBIFAM; homogeneous HMMER bits —
  # both evalue AND score are meaningful, unlike the cross-library interpro edge).
  properties:
    start: int
    end: int
    evalue: float
    score: float

ncbifam family in interpro entry:
  is_a: association
  represented_as: edge
  label_as_edge: Ncbifam_family_in_interpro_entry
  source: ncbifam family
  target: interpro entry
  label_in_input: ncbifam_family_in_interpro_entry
```

- [ ] **Step 2: Wiring.** In `create_knowledge_graph.py`, after the interpro adapter block: expose the interpro kept set (add a small method `kept_node_accessions() -> set[str]` on `MultiInterproAnnotationAdapter` returning `self._kept_ids(self._observed_ids())` — add it in this task with a 3-line unit test in `tests/test_interpro_adapter.py`), then:

```python
    from multiomics_kg.adapters.ncbifam_adapter import MultiNcbifamAdapter
    ncbifam_adapter = MultiNcbifamAdapter(
        genome_config_file=GENOME_CONFIG,
        cache_root="cache/data",
        interpro_kept_ids=interpro_adapter.kept_node_accessions(),
        test_mode=TEST_MODE,
    )
    ncbifam_adapter.download_data(cache=CACHE)
    bc.write_nodes(ncbifam_adapter.get_nodes())
    bc.write_edges(ncbifam_adapter.get_edges())
```

(Match the surrounding blocks' argument names — copy the interpro block's actual variable names for `GENOME_CONFIG`/`TEST_MODE`/`CACHE` from the file.)
- [ ] **Step 3: Smoke** — `uv run python create_knowledge_graph.py --test 2>&1 | tail -20` → completes; log lines show NcbifamFamily nodes + edges written.
- [ ] **Step 4: Commit** — `git commit -m "feat(schema+kg): NcbifamFamily node/edges wired; interpro edge evalue_library"`

---

### Task 15: `.gitignore` + limited-file rule check for interproscan (housekeeping)

**Files:**
- Modify: `.gitignore` (verify `cache/data/*/genomes/*/interproscan/*.limited_*` rule exists — add if missing; `cache/data/ncbifam/raw/` was Task 5)

- [ ] **Step 1:** `grep -n "interproscan" .gitignore` — confirm raw.json + limited rules; add whatever's missing.
- [ ] **Step 2: Commit** (or skip if no change).

---

### Task 16: Post-import Cypher — indexes, rollups, buckets, F1.1

**Files:**
- Modify: `scripts/post-import.cypher` AND `scripts/post-import.sh` (identical logic — edit both; the .sh groups statements into cypher-shell invocations)

Changes (each mirrored in both files):

- [ ] **Step 1: Indexes** (next to the InterproEntry block at ~line 78):

```cypher
// NCBIfam (flat curated-family ontology; family_type is the stratification key)
CREATE INDEX ncbifam_family_id_idx IF NOT EXISTS FOR (n:NcbifamFamily) ON (n.ncbifam_id);
CREATE INDEX ncbifam_family_type_idx IF NOT EXISTS FOR (n:NcbifamFamily) ON (n.family_type);
CREATE INDEX ncbifam_family_level_idx IF NOT EXISTS FOR (n:NcbifamFamily) ON (n.level);
CREATE FULLTEXT INDEX ncbifamFamilyFullText IF NOT EXISTS
    FOR (n:NcbifamFamily) ON EACH [n.name, n.gene_symbol, n.description];
```

Extend `interproEntryFullText` → `ON EACH [e.name, e.description]` (drop+recreate semantics: `DROP INDEX interproEntryFullText IF EXISTS;` before the CREATE — full-text defs can't be altered).
- [ ] **Step 2: NcbifamFamily rollups** (next to the InterproEntry computed block ~line 967):

```cypher
// ── NcbifamFamily computed properties (flat; direct counts) ───
MATCH (n:NcbifamFamily)
CALL (n) {
  OPTIONAL MATCH (n)<-[:Gene_has_ncbifam_family]-(g:Gene)
  RETURN count(DISTINCT g) AS gc, count(DISTINCT g.organism_name) AS oc
}
SET n.gene_count = gc, n.organism_count = oc;
```

(Copy the exact `CALL (n) { }` / `IN TRANSACTIONS` idiom used by the neighboring InterproEntry block — match it verbatim.)
- [ ] **Step 3: Gene routing** — next to `interpro_entry_count` (~line 1216) add `ncbifam_family_count` (same OPTIONAL MATCH/count pattern on `Gene_has_ncbifam_family`). In the `annotation_types` CASE list (~line 661) add `CASE WHEN EXISTS { (g)-[:Gene_has_ncbifam_family]->() } THEN ['ncbifam'] ELSE [] END`.
- [ ] **Step 4: F1.1 flags** (after the KEGG pattern rule ~line 564):

```cypher
// InterPro: unknown-function entries (name-pattern rule; uninformative_terms.yaml)
MATCH (t:InterproEntry)
WHERE t.name =~ '^Protein of unknown function.*'
   OR t.name =~ '^Domain of unknown function.*'
   OR t.name =~ '^Uncharacteri[sz]ed protein family.*'
SET t.is_uninformative = 'true';

// NCBIfam: typed rule (family_type — third rule kind) + name-pattern fallback
MATCH (t:NcbifamFamily)
WHERE t.family_type IN ['hypoth_equivalog', 'hypoth_equivalog_domain']
   OR t.name =~ '(?i).*hypothetical.*'
   OR t.name =~ '(?i).*uncharacterized.*'
   OR t.name =~ '.*DUF\\d.*'
SET t.is_uninformative = 'true';
```

- [ ] **Step 5: 9th bucket + has_any_edge** (the F1.2/F1.3 block ~lines 598–640):
  - Add to the WITH block: `EXISTS { (g)-[:Gene_has_ncbifam_family]->(t) WHERE t.is_uninformative IS NULL } AS has_ncbifam,`
  - Add `+ CASE WHEN has_ncbifam THEN 1 ELSE 0 END` to the informative_count sum.
  - Extend the `has_any_edge` relationship list with `|Gene_has_interpro_entry|Gene_has_ncbifam_family`.
  - Update the `SOURCE_BUCKETS` comment: `live (9): go, kegg, pfam, ec, role, reaction, transporter, cazy, ncbifam`.
  - In the F1.4 `informative_annotation_types` block add the parallel `ncbifam` source entry (same EXISTS-with-filter shape as the others).
- [ ] **Step 6: Diff-check the two files agree** — `bash -c "diff <(grep -o 'Gene_has_ncbifam[a-z_]*' scripts/post-import.cypher | sort -u) <(grep -o 'Gene_has_ncbifam[a-z_]*' scripts/post-import.sh | sort -u)"` → empty.
- [ ] **Step 7: Commit** — `git commit -m "feat(post-import): NcbifamFamily rollups/indexes, 9th bucket, has_any_edge fix, F1.1 rules"`

---

### Task 17: Bucket-count + fast-suite green gate

**Files:**
- Modify: whichever test pins the source-bucket count (find with `grep -rn "informative" tests/ | grep -i "bucket\|source"` — CLAUDE.md names a "bucket-count test") — update 8 → 9 and the bucket list.

- [ ] **Step 1:** Update the bucket-count test to 9 buckets including `ncbifam`.
- [ ] **Step 2:** Full fast suite green: `uv run pytest -m "not slow and not kg" -q` → all pass (consistency test may still skew until Task 18 — if red, confirm it's ONLY that file).
- [ ] **Step 3: Commit.**

---

### Task 18: Pipeline run — references + merge + consistency green

- [ ] **Step 1:** `bash scripts/prepare_data.sh --steps 9 --force` (both references rebuild; check `logs/prepare_data_step9.log`).
- [ ] **Step 2:** `bash scripts/prepare_data.sh --steps 2 --force` (merge all strains; ~check log for errors).
- [ ] **Step 3: Layer-2 spot checks** (spec §6): MED4 merged JSON —
  - ATP-synthase-c gene (protein `WP_002805169.1`): `ec_numbers` contains `7.1.2.2` and `ec_numbers_source["7.1.2.2"]` includes `interpro`.
  - Argininosuccinate-lyase gene (protein `WP_011131650.1`): none of IPR000362's 5 ECs is interpro-attributed.
  - MsrB gene (protein `WP_011131653.1`): `ncbifam_ids` contains `TIGR00357`.
  - Count skipped `[ncbifam]` echo tokens across MED4 (product == token): assert > 0 via a one-off script.
- [ ] **Step 4:** Consistency test green: `uv run pytest tests/test_interproscan_consistency.py -q` → PASS all strains.
- [ ] **Step 5:** Commit regenerated merged JSONs + resolved artifacts: `git add -A cache/data && git commit -m "data(merge): regenerate gene_annotations_merged with ncbifam/GO-donor/naming-recovery fields"`

---

### Task 19: Docker rebuild + KG validation

- [ ] **Step 1: Pre-rebuild snapshot** — `/omics-edge-snapshot` (record baseline).
- [ ] **Step 2: Rebuild** — stop `deploy`/`app`, run build → import → post-process → deploy (`docker compose up -d` per CLAUDE.md; re-running import needs deploy/app stopped first).
- [ ] **Step 3: Import gate** — `cat output/import.status` = 0; `import.report` has ZERO skipped relationships.
- [ ] **Step 4: Expected count movements** (spec §6 — judge against predictions): NcbifamFamily ~4,957 nodes; Gene_has_ncbifam_family ~67K; Gene_has_interpro_entry ~397K (~unchanged); ~185 bucket climbs; ~367–420 `no_evidence`→`catch_all_only` movers. Query each; deviations beyond ~5% → investigate before proceeding.
- [ ] **Step 5: `/omics-edge-snapshot` after** — expression edges byte-identical to Step 1.
- [ ] **Step 6: KG validity** — add Layer-3 spot-check assertions to `tests/kg_validity/` (new file `test_ncbifam.py`): katG edges exist for EZ55+MIT1002, zero for Prochlorococcus; `TIGR00198` node has `family_type='equivalog'` + non-null description; KT2440 `WP_010953880.1` GO edge to `GO:0022857` with `evidence='domain_inferred'` and `'interpro' IN sources`; ≥1 hypoth_equivalog node flagged `is_uninformative`; DUF-named InterproEntry flagged; NcbifamFamily count in [4900, 5050]. Run `uv run pytest -m kg -v` → green (except pre-existing known failures per CLAUDE.md orphan-proteins note).
- [ ] **Step 7: Snapshot regen** — `uv run python tests/kg_validity/generate_snapshot.py` (InterPro-sampled entries legitimately changed); `uv run pytest tests/kg_validity/test_snapshot.py -q` → PASS.
- [ ] **Step 8: Commit** tests + snapshot.

---

### Task 20: post-import validation dump comparison (regression guard)

- [ ] **Step 1:** `scripts/post-import-validate.sh > after.txt` against the rebuilt graph; diff against a pre-change baseline IF one was captured — this change is NOT a pure refactor, so differences are expected in interpro/quality sections; verify no UNEXPECTED sections changed (expression ranks, TCDB, metabolites untouched).
- [ ] **Step 2:** Note results in the plan/PR description.

---

### Task 21: Documentation + memory

**Files:**
- Create: `docs/kg-changes/interpro-multi-ontology.md`
- Modify: `docs/kg-changes/interproscan-extension.md`, `docs/kg-changes/interpro-two-layer.md` (superseded headers)
- Modify: `CLAUDE.md`
- Memory: update `MEMORY.md` pointer if needed

- [ ] **Step 1:** Write `docs/kg-changes/interpro-multi-ontology.md` — what shipped vs the spec: artifact format, NcbifamFamily counts (actual from Task 19), bucket resolution actuals, edge property changes (`evalue_library`, no `score`), naming-recovery stats (skipped-echo counts), links to the spec + plan.
- [ ] **Step 2:** Prepend to both old docs: `> **SUPERSEDED (2026-08-17)** by the multi-ontology redesign — see docs/kg-changes/interpro-multi-ontology.md. This document describes the PREVIOUS integration.`
- [ ] **Step 3:** CLAUDE.md updates: adapters list (`ncbifam_adapter.py`; interpro_adapter description), Neo4j labels (+`NcbifamFamily`, +`Gene_has_ncbifam_family`, +`Ncbifam_family_in_interpro_entry`), Key graph facts (new NcbifamFamily bullet; InterproEntry bullet updated: description, evalue_library/no-score, calls format; 9-bucket list + `ncbifam`; `has_any_edge` note; prepare_data step 9 = central references ordered before merge; step-2 mentions ncbifam reference consumption; Data Locations: `ncbifam_reference.json` + raw gitignored; interproscan calls.json format note).
- [ ] **Step 4:** Final full-suite run: `uv run pytest -m "not slow and not kg" -q` green.
- [ ] **Step 5: Commit** — `git commit -m "docs(interpro): multi-ontology redesign shipped — kg-changes doc, CLAUDE.md, superseded headers"`

---

## Self-Review (done at write time)

- **Spec coverage**: §3.1 → Tasks 1–3; §3.2 → Tasks 5–6 (+ no-pruned-artifact = no task, by design); §4 → Tasks 7–9, 18; §4.2 consistency → Task 10; §5.1–5.3 → Tasks 12–14; §5.4 → Task 16; §5.5 → Tasks 11+16; §5.6 → Tasks 9+16; §6 → Tasks 17–20; §7 → Tasks 3 (entry_xrefs deletion), 4, 21; §8 → out of scope, no tasks. Gap check: none found.
- **Type consistency**: `parse_interproscan_json` record shape (Task 1) = what Task 8 loads = what Task 12/13 adapters read; `ncbifam_ids` field name consistent across Tasks 8/9/10/13/16; node id `ncbifam_<acc>` consistent (Tasks 13/14/16); `evalue_library` consistent (Tasks 1/12/14).
- **Placeholders**: none; every code step has concrete content; two deliberate empirical checks are flagged as such (gene_synonyms delimiter in Task 5; genomes-CSV header in Task 10).
