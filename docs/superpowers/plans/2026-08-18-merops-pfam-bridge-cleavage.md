# MEROPS Pfam Bridge + Cleavage Specificity Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give the MEROPS layer a curated family→Pfam corroboration bridge
(`Merops_family_has_pfam_domain` + edge `pfam_support`), cleavage-specificity
properties on MeropsFamily nodes, contract registration for all merops
vocabularies, and two prepare_data improvements (step 9→10 move, `--rebuild`).

**Architecture:** Pure parsers in `utils/merops_diamond.py` feed
`build_merops_reference.py` (prepare_data step 10), which extends the committed
`merops_reference.json` with `pfam_bridge` + `cleavage` blocks; `merops_adapter`
emits the bridge edges (pruned by injected Pfam node ids, TCDB pattern) and the
node properties; post-import Cypher derives `pfam_support` per gene edge.

**Tech Stack:** Python 3.13 / pytest / BioCypher / Neo4j Cypher / bash.

**Spec:** `docs/superpowers/specs/2026-08-18-merops-pfam-bridge-cleavage-design.md`
(read it first — measured inputs, rejections, and vocabulary shapes live there).

## Global Constraints

- String properties never contain `'` or `|` (`_clean_str` / no-pipe-separator
  rule) — applies to computed strings like `cleavage_summary`.
- No native booleans anywhere (vocabulary contract R5): two-state facts are
  meaningful string pairs (`corroborated | uncorroborated`).
- `scripts/post-import.sh` and `scripts/post-import.cypher` must stay
  byte-identical in their Cypher logic (`diff` of merops-greps must be empty).
- Sparse properties are ABSENT, never empty strings/lists.
- Categorical literals emitted by adapters go through
  `VOCAB.check(applies_to, prop, value)` (`multiomics_kg/utils/controlled_vocab.py`).
- New-source args are keyword/optional so existing positional callers keep
  working.
- Every commit message ends with the Claude co-author line (repo convention).

---

### Task 1: Pure parsers — `parse_interpro_txt_stream`, `aggregate_cleavages`, `cleavage_properties`

**Files:**
- Modify: `multiomics_kg/utils/merops_diamond.py` (append after `type_example_names`)
- Test: `tests/test_merops_diamond.py` (append)

**Interfaces:**
- Produces: `parse_interpro_txt_stream(lines: Iterable[str]) -> dict[str, dict[str, int]]`
  — `{family: {pfam_acc: member_id_count}}`, member_id_count = distinct MEROPS
  identifiers backing the (family, Pfam) pair.
- Produces: `aggregate_cleavages(lines: Iterable[str]) -> dict[str, dict]`
  — `{family: {"p1": Counter, "physiological": int, "total": int}}`; P1 counts
  keep only `STANDARD_AA_3` residues.
- Produces: `cleavage_properties(agg: dict) -> dict` — the three node
  properties for ONE family's aggregate (empty dict when `total == 0`).
- Produces: `STANDARD_AA_3: frozenset[str]` (20 three-letter codes).

- [ ] **Step 1: Write the failing tests** (append to `tests/test_merops_diamond.py`)

```python
from multiomics_kg.utils.merops_diamond import (
    STANDARD_AA_3,
    aggregate_cleavages,
    cleavage_properties,
    parse_interpro_txt_stream,
)


def test_parse_interpro_txt_stream_counts_distinct_identifiers():
    lines = [
        '"A01A","A01.001","pepsin A","63-388","P0DJD8","PF00026/PF14543","IPR001461"',
        '"A01A","A01.001","pepsin A","63-388","B7Z719","PF00026/PF14543","IPR001461"',  # same id, dup accession
        '"A01A","A01.002","pepsin B","60-380","P27821","PF00026","IPR001461"',
        '"S01A","S01.001","chymotrypsin A","34-263","P00766","PF00089",""',
        '"M10B","M10.051","serralysin","1-470","P23694","",""',  # no Pfam -> ignored
    ]
    bridge = parse_interpro_txt_stream(lines)
    assert bridge["A01"]["PF00026"] == 2      # A01.001 + A01.002, dedup on identifier
    assert bridge["A01"]["PF14543"] == 1      # only A01.001
    assert bridge["S01"]["PF00089"] == 1
    assert "M10" not in bridge


def test_aggregate_cleavages_filters_nonstandard_p1():
    q = chr(39)  # single quote — the file quotes fields with it
    def row(ident, p1, kind):
        f = [q + "CLE1" + q, q + ident + q, "s", "s", "-", "-", "-",
             q + p1 + q, "-", "-", "-", "-", "ref", "NULL", "NULL", "NULL",
             "e", "NULL", "NULL", "NULL", "NULL", "NULL", q + kind + q]
        return "\t".join(f)
    agg = aggregate_cleavages([
        row("S01.001", "Lys", "physiological"),
        row("S01.001", "Arg", "synthetic"),
        row("S01.151", "TyI", "non-physiological"),   # modified residue: counted in total, not in p1
        row("A01.001", "Phe", "physiological"),
    ])
    assert agg["S01"]["total"] == 3
    assert agg["S01"]["physiological"] == 1
    assert agg["S01"]["p1"]["Lys"] == 1 and agg["S01"]["p1"]["Arg"] == 1
    assert "TyI" not in agg["S01"]["p1"]
    assert agg["A01"]["total"] == 1


def test_cleavage_properties_shapes():
    from collections import Counter
    props = cleavage_properties({
        "p1": Counter({"Lys": 36, "Arg": 34, "Glu": 11, "Ala": 9}),
        "physiological": 25, "total": 100,
    })
    assert props["cleavage_p1_residues"] == ["Lys", "Arg", "Glu"]  # >=10% share, max 3 (Ala at exactly 10% loses to top-3)
    assert props["known_cleavage_count"] == 100
    s = props["cleavage_summary"]
    # shares over the standard-P1 subtotal (90): 36/90=40%, 34/90=38%, 11/90=12%
    assert s == "cleaves after Lys (40%) / Arg (38%) / Glu (12%) - 100 known cleavages (25% physiological)"
    assert "'" not in s and "|" not in s
    # no standard-P1 data -> count-only summary, no residues key
    props2 = cleavage_properties({"p1": Counter(), "physiological": 3, "total": 10})
    assert "cleavage_p1_residues" not in props2
    assert props2["cleavage_summary"] == "10 known cleavages (30% physiological)"
    assert props2["known_cleavage_count"] == 10
    # no data at all -> empty (sparse discipline)
    assert cleavage_properties({"p1": Counter(), "physiological": 0, "total": 0}) == {}
```

Note the percentages: shares are computed over the **standard-residue P1
subtotal** (90 here), not `total` — 36/90 = 40%.

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_merops_diamond.py -q -k "interpro_txt or cleavage"`
Expected: FAIL with ImportError (names not defined).

- [ ] **Step 3: Implement in `multiomics_kg/utils/merops_diamond.py`**

Append (uses existing `id_family`, `re`, `collections` imports — add
`import csv` and `from collections import Counter` at the top if absent):

```python
# ============================================================================
# Phase-2 follow-up: Pfam bridge + cleavage specificity parsers
# (spec docs/superpowers/specs/2026-08-18-merops-pfam-bridge-cleavage-design.md)
# ============================================================================

STANDARD_AA_3 = frozenset({
    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val",
})


def parse_interpro_txt_stream(lines) -> dict[str, dict[str, int]]:
    """MEROPS interpro.txt → {family: {pfam_acc: member_id_count}}.

    One input row per curated UniProt member: subfam token, MEROPS identifier,
    name, peptidase-unit range, accession, slash-joined Pfam accessions, IPR.
    member_id_count = DISTINCT MEROPS identifiers backing the (family, Pfam)
    pair — accessions are heavily redundant (isoforms/orthologs), identifiers
    are the curated unit. Stream-safe: never holds the 182 MB file.
    """
    backing: dict[str, dict[str, set[str]]] = {}
    for row in csv.reader(lines):
        if len(row) < 6:
            continue
        merops_id = row[1].strip()
        family = id_family(merops_id)
        if family is None:
            continue
        for pf in row[5].split("/"):
            pf = pf.strip()
            if pf.startswith("PF"):
                backing.setdefault(family, {}).setdefault(pf, set()).add(merops_id)
    return {
        fam: {pf: len(ids) for pf, ids in sorted(pfams.items())}
        for fam, pfams in sorted(backing.items())
    }


def aggregate_cleavages(lines) -> dict[str, dict]:
    """MEROPS Substrate_search.txt → per-family cleavage aggregates.

    Tab-separated, single-quote-quoted; col 1 = MEROPS identifier, col 7 = P1
    residue, col 22 = record kind (physiological | non-physiological |
    synthetic | theoretical | pathological). Returns
    {family: {"p1": Counter, "physiological": int, "total": int}}. P1 counts
    keep only the 20 standard residues (STANDARD_AA_3) so the downstream
    vocabulary can be closed; modified residues (e.g. "TyI") still count
    toward "total". Caller decodes latin-1.
    """
    agg: dict[str, dict] = {}
    for line in lines:
        parts = [p.strip().strip("'") for p in line.split("\t")]
        if len(parts) < 23:
            continue
        family = id_family(parts[1])
        if family is None:
            continue
        rec = agg.setdefault(
            family, {"p1": Counter(), "physiological": 0, "total": 0}
        )
        rec["total"] += 1
        if parts[22] == "physiological":
            rec["physiological"] += 1
        if parts[7] in STANDARD_AA_3:
            rec["p1"][parts[7]] += 1
    return agg


def cleavage_properties(agg: dict) -> dict:
    """One family's aggregate → the three sparse MeropsFamily node properties.

    cleavage_p1_residues: P1 residues with >= 10% share of the standard-P1
    subtotal, max 3, frequency order. cleavage_summary: readable sentence
    (no ' or | — the clean-string rule applies to computed strings).
    known_cleavage_count: ALL curated records ("known cleavages" is MEROPS's
    own phrasing). Empty dict when the family has no records at all.
    """
    total = agg.get("total", 0)
    if total == 0:
        return {}
    phys_pct = round(100 * agg.get("physiological", 0) / total)
    tail = f"{total} known cleavages ({phys_pct}% physiological)"
    p1: Counter = agg.get("p1") or Counter()
    subtotal = sum(p1.values())
    top = [
        (res, round(100 * n / subtotal))
        for res, n in p1.most_common(3)
        if n / subtotal >= 0.10
    ] if subtotal else []
    props: dict = {"known_cleavage_count": total}
    if top:
        clause = " / ".join(f"{res} ({pct}%)" for res, pct in top)
        props["cleavage_p1_residues"] = [res for res, _ in top]
        props["cleavage_summary"] = f"cleaves after {clause} - {tail}"
    else:
        props["cleavage_summary"] = tail
    return props
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_merops_diamond.py -q`
Expected: all PASS (old Phase-1 tests included).

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/merops_diamond.py tests/test_merops_diamond.py
git commit -m "feat(merops): interpro.txt + Substrate_search.txt parsers (Pfam bridge, cleavage specificity)"
```

---

### Task 2: Builder — `pfam_bridge` + `cleavage` blocks in `merops_reference.json`

**Files:**
- Modify: `multiomics_kg/download/build_merops_reference.py`
- Regenerates: `cache/data/merops/merops_reference.json` (committed)

**Interfaces:**
- Consumes: Task 1's three parser functions.
- Produces: `merops_reference.json` gains top-level `"pfam_bridge":
  {family: {pfam_acc: member_id_count}}` and `"cleavage":
  {family: {cleavage_* props}}` — read by Task 6's adapter.

- [ ] **Step 1: Add the two raw files to `RAW_FILES`**

In `build_merops_reference.py`, extend the mapping (note: neither file exists
in the Phase-1 `$MEROPS_DATA_DIR/DB/` install, so `_locate_raw` will download
them into the gitignored `cache/data/merops/raw/` — 182 MB + 29 MB, one-time):

```python
RAW_FILES = {
    "family.txt": (f"{MEROPS_RELEASE_URL}/database_files/family.txt", "family.txt"),
    "clan.txt": (f"{MEROPS_RELEASE_URL}/database_files/clan.txt", "clan.txt"),
    # the Phase-1 runner keeps the pristine download as merops_scan.raw.lib
    "merops_scan.lib": (f"{MEROPS_RELEASE_URL}/merops_scan.lib", "merops_scan.raw.lib"),
    # Pfam bridge + cleavage specificity (2026-08-18 follow-up). Neither ships
    # with the Phase-1 tool install — they download straight into raw/.
    "interpro.txt": (f"{MEROPS_RELEASE_URL}/interpro.txt", "interpro.txt"),
    "Substrate_search.txt": (f"{MEROPS_RELEASE_URL}/Substrate_search.txt", "Substrate_search.txt"),
}
```

- [ ] **Step 2: Extend `build()`**

Import the new parsers at the top:

```python
from multiomics_kg.utils.merops_diamond import (
    aggregate_cleavages,
    cleavage_properties,
    parse_clan_txt,
    parse_family_txt,
    parse_interpro_txt_stream,
    type_example_names,
)
```

After the `example_names = ...` line, add (stream the big file — no
`read_text` on 182 MB):

```python
    interpro_path = _locate_raw("interpro.txt", refetch_raw)
    with open(interpro_path, encoding="latin-1", newline="") as fh:
        pfam_bridge = parse_interpro_txt_stream(fh)

    substrate_path = _locate_raw("Substrate_search.txt", refetch_raw)
    with open(substrate_path, encoding="latin-1") as fh:
        cleavage_agg = aggregate_cleavages(fh)
    cleavage = {
        fam: props
        for fam, agg in sorted(cleavage_agg.items())
        if (props := cleavage_properties(agg))
    }
```

Add both to the `ref` dict (after `"clans": clans,`):

```python
        "pfam_bridge": pfam_bridge,
        "cleavage": cleavage,
```

Extend the QC block (before the final write):

```python
    if not pfam_bridge:
        logger.warning("Parsed 0 pfam_bridge pairs — did interpro.txt column layout change?")
    if not cleavage:
        logger.warning("Parsed 0 cleavage profiles — did Substrate_search.txt column layout change?")
```

and the final log line gains
`"%d pfam-bridge families, %d cleavage families"` with
`len(pfam_bridge), len(cleavage)` appended to the existing args.

- [ ] **Step 3: Rebuild the committed reference**

Run: `uv run python -m multiomics_kg.download.build_merops_reference --force`
Expected: downloads interpro.txt (182 MB, minutes) + Substrate_search.txt
(29 MB); log line reports ~320 pfam-bridge families and ~250 cleavage
families; zero QC warnings.

- [ ] **Step 4: Spot-check the committed output**

```bash
jq '.pfam_bridge.S14' cache/data/merops/merops_reference.json     # {"PF00574": <n>}
jq '.cleavage.S01' cache/data/merops/merops_reference.json
```
Expected: S14 → exactly `PF00574`; S01 cleavage has `cleavage_p1_residues`
starting `["Lys", "Arg", ...]`, `known_cleavage_count` ≈ 39567, summary
containing "known cleavages".

- [ ] **Step 5: Run the fast merops tests (no regressions)**

Run: `uv run pytest tests/test_merops_diamond.py tests/test_merops_adapter.py -q`
Expected: PASS (adapter ignores unknown reference keys).

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/download/build_merops_reference.py cache/data/merops/merops_reference.json
git commit -m "feat(merops): reference gains pfam_bridge + cleavage blocks (interpro.txt, Substrate_search.txt)"
```

---

### Task 3: prepare_data — MEROPS reference → step 10, new `--rebuild` flag

**Files:**
- Modify: `scripts/prepare_data.sh`
- Modify: `multiomics_kg/adapters/merops_adapter.py` (error message only)

**Interfaces:**
- Produces: `bash scripts/prepare_data.sh --steps 10` builds the MEROPS
  reference; `--rebuild` ≡ `--steps 9 1 2 3 4 5 6 7 8 10 --force`.

- [ ] **Step 1: Move the merops sub-builder out of step 9 into a new `10)` case**

Delete the third `RUN_STEP_APPEND=1 run_step 9 ... build_merops_reference`
block from case `9)` and add after it:

```bash
        10)
            run_step 10 \
                "Build MEROPS reference cache (merops_reference.json + pfam_bridge + cleavage)" \
                "$LOG_DIR/prepare_data_step10.log" \
                uv run python -m multiomics_kg.download.build_merops_reference \
                    $FORCE \
                    $REFETCH_RAW
            ;;
```

Update the `*)` error line to `(valid: 0 1 2 3 4 5 6 7 8 9 10)`.

- [ ] **Step 2: Default steps, labels, usage**

- `STEPS="0 9 1 2 3 4 5 6 7 8"` → `STEPS="0 9 1 2 3 4 5 6 7 8 10"` and update
  the comment above it: step 9 = merge-consumed central references (InterPro +
  NCBIfam, BEFORE step 2); step 10 = MEROPS reference (KG-build-time only, no
  ordering constraint).
- The `echo "(step 1 = ..."` line: change the step-9 clause to
  `step 9 = central references: InterPro + NCBIfam reference caches, step 10 = MEROPS reference cache`.
- Usage comment block: add
  `#   ./scripts/prepare_data.sh --rebuild                          # all derived steps (1-10, 9 before 2) with --force; step-0 downloads excluded`.

- [ ] **Step 3: Add `--rebuild`**

In the arg parser, add a tracking var and the flag (order-independent conflict
detection):

```bash
USER_STEPS=0
REBUILD=0
```

(next to `FORCE=""`), then in the `case`:

```bash
        --rebuild)        REBUILD=1; shift ;;
```

change the `--steps` branch to also set `USER_STEPS=1`, and after the parse
loop:

```bash
if [[ $REBUILD -eq 1 ]]; then
    if [[ $USER_STEPS -eq 1 ]]; then
        echo "--rebuild and --steps are mutually exclusive" >&2; exit 1
    fi
    # All derived steps in DEPENDENCY order (9 before 2 — step 2 consumes the
    # step-9 reference caches). Step 0 (raw downloads) deliberately excluded.
    STEPS="9 1 2 3 4 5 6 7 8 10"
    FORCE="--force"
fi
```

- [ ] **Step 4: Update the adapter's fail-loudly message**

In `merops_adapter.py` `download_data()`, change
`--steps 9` → `--steps 10` in the `FileNotFoundError` text.

- [ ] **Step 5: Verify**

```bash
bash -n scripts/prepare_data.sh                                   # syntax
bash scripts/prepare_data.sh --rebuild --steps 2 ; echo "exit=$?" # expect error + exit 1
bash scripts/prepare_data.sh --steps 10                           # cache exists -> fast no-op message
uv run pytest tests/test_merops_adapter.py -q                     # error-message test still matches "prepare_data"
```
Expected: syntax OK; mutual-exclusion error; step 10 reports the existing
cache; tests PASS (the test regex matches `prepare_data`, not the step number).

- [ ] **Step 6: Commit**

```bash
git add scripts/prepare_data.sh multiomics_kg/adapters/merops_adapter.py
git commit -m "feat(prepare-data): merops reference moves to step 10; add --rebuild flag"
```

---

### Task 4: Controlled-vocabulary registration (10 entries)

**Files:**
- Modify: `config/controlled_vocabularies.yaml` (append a merops section)

**Interfaces:**
- Produces: `VOCAB.get("Gene_has_merops_family", "call_class")` etc. resolve;
  Task 6 emits literals through `VOCAB.check`.

- [ ] **Step 1: Append the merops section**

Append at the end of `config/controlled_vocabularies.yaml`:

```yaml
# ── MEROPS peptidase layer (shipped 2026-08-17; registered + extended
#    2026-08-18, spec docs/superpowers/specs/2026-08-18-merops-pfam-bridge-cleavage-design.md) ──

Gene_has_merops_family.call_class:
  applies_to: Gene_has_merops_family
  applies_to_kind: edge
  property: call_class
  value_type: string
  closed: true
  values: [peptidase, inhibitor, nonpeptidase_homolog]
  description: >
    Read-first verdict for a MEROPS diamond call. nonpeptidase_homolog = the
    best hit is a catalytically dead .9xx relative — fold evidence, NOT
    protease evidence. inhibitor = I-family (protease inhibitors are biology).
    Threshold-free: derived from best_hit_kind + the family type.

Gene_has_merops_family.tier:
  applies_to: Gene_has_merops_family
  applies_to_kind: edge
  property: tier
  value_type: int
  closed: false
  values: []
  min_value: 1
  max_value: 3
  description: >
    tcdb-diamond confidence band and truncation depth: 1 = MEROPS identifier
    (identity >= 70), 2 = subfamily (identity >= 40), 3 = family (floor only).
    NOT sparse — unlike Gene_has_tcdb_family.tier, every merops edge is
    diamond-sourced and carries one. ~92% tier 3 (expected remote homology).

Gene_has_merops_family.best_hit_kind:
  applies_to: Gene_has_merops_family
  applies_to_kind: edge
  property: best_hit_kind
  value_type: string
  closed: true
  values: [holotype, putative, nonpeptidase_homolog]
  description: >
    How well-characterized the matched MEROPS reference is: holotype =
    experimentally characterized peptidase (.001-.899), putative = predicted
    only (letter tails), nonpeptidase_homolog = proven dead relative (.9xx).

Gene_has_merops_family.pfam_support:
  applies_to: Gene_has_merops_family
  applies_to_kind: edge
  property: pfam_support
  value_type: string
  closed: true
  values: [corroborated, uncorroborated]
  description: >
    Whether a Pfam domain on this gene is curated into this MEROPS family via
    Merops_family_has_pfam_domain. Read in the sound direction only.

MeropsFamily.level_kind:
  applies_to: MeropsFamily
  applies_to_kind: node
  property: level_kind
  value_type: string
  closed: true
  values: [merops_clan, merops_family, merops_subfamily]
  description: MEROPS hierarchy depth name (level 0 / 1 / 2).

MeropsFamily.family_type:
  applies_to: MeropsFamily
  applies_to_kind: node
  property: family_type
  value_type: string
  closed: true
  values: [peptidase, inhibitor]
  description: >
    Whether this clan/family/subfamily holds peptidases or protease
    inhibitors (I-prefixed codes).

MeropsFamily.catalytic_type:
  applies_to: MeropsFamily
  applies_to_kind: node
  property: catalytic_type
  value_type: string
  closed: true
  values: [serine, cysteine, metallo, aspartic, threonine, glutamic,
           asparagine_lyase, mixed, unknown]
  sparse: true
  description: >
    Catalytic mechanism as a full word (MEROPS's one-letter codes spelled
    out). SPARSE — absent on inhibitor families, which have no catalytic type.

MeropsFamily.cleavage_p1_residues:
  applies_to: MeropsFamily
  applies_to_kind: node
  property: cleavage_p1_residues
  value_type: string_array
  closed: true
  values: [Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile,
           Leu, Lys, Met, Phe, Pro, Ser, Thr, Trp, Tyr, Val]
  sparse: true
  description: >
    Top P1 residues (Schechter-Berger: the residue N-terminal of the scissile
    bond) with >= 10% share of curated cleavages, max 3, frequency order.
    Filtered to the 20 standard residues. SPARSE — level-1 families with
    curated cleavage data only; read known_cleavage_count for support.

Merops_family_has_pfam_domain.member_id_count:
  applies_to: Merops_family_has_pfam_domain
  applies_to_kind: edge
  property: member_id_count
  value_type: int
  closed: false
  values: []
  min_value: 1
  description: >
    How many distinct curated MEROPS identifiers back this (family, Pfam)
    pair in MEROPS's interpro.txt. Support size, not a filter.

Gene.merops_classes:
  applies_to: Gene
  applies_to_kind: node
  property: merops_classes
  value_type: string_array
  closed: true
  values: [peptidase, inhibitor, nonpeptidase_homolog]
  description: >
    Sorted distinct call_class values across the gene's Gene_has_merops_family
    edges — the at-a-glance guard so a gene whose only call is a dead homolog
    never reads as "1 protease". Empty when the gene has no merops edge.
```

- [ ] **Step 2: Run the contract unit tests**

Run: `uv run pytest tests/test_controlled_vocab.py tests/test_controlled_vocabulary_adapter.py -q`
Expected: PASS (loader validates required keys/value_types; no count
assertion exists — verify with the output; if a count test DOES fail, bump it
by 10 and note it in the commit message).

- [ ] **Step 3: Commit**

```bash
git add config/controlled_vocabularies.yaml
git commit -m "feat(vocab): register merops vocabularies (7 shipped + 3 new entries)"
```

---

### Task 5: Schema — bridge edge + three node properties

**Files:**
- Modify: `config/schema_config.yaml`

**Interfaces:**
- Produces: BioCypher accepts `merops_family_has_pfam_domain` edge tuples and
  the three new node properties Task 6 emits.

- [ ] **Step 1: Add the three properties to the `merops family:` node block**

After `catalytic_type: str` add:

```yaml
    cleavage_p1_residues: str[]  # sparse — top P1 residues (>=10% share, max 3, 20 standard AAs)
    cleavage_summary: str        # sparse — readable specificity sentence ("cleaves after ... N known cleavages")
    known_cleavage_count: int    # sparse — curated cleavage records backing the profile
```

- [ ] **Step 2: Add the bridge edge after `merops family hierarchical association:`**

```yaml
# From MEROPS's published interpro.txt member map, aggregated to family level
# (183 pairs, median 1 Pfam/family). Semantics: "MEROPS's curated member
# proteins of this family carry this domain" — Tcdb_family_has_pfam_domain
# shape. Sound direction only: corroborate a gene's known family outward;
# never assign peptidase identity backward from a domain (shared folds —
# S09 alpha/beta-hydrolase, C26 GATase — are exactly the call_class trap).
merops family to pfam association:
  is_a: association
  represented_as: edge
  label_as_edge: Merops_family_has_pfam_domain
  source: merops family
  target: pfam
  label_in_input: merops_family_has_pfam_domain
  properties:
    member_id_count: int
```

- [ ] **Step 3: Sanity-check the yaml loads**

Run: `uv run python -c "import yaml; yaml.safe_load(open('config/schema_config.yaml')); print('OK')"`
Expected: OK

- [ ] **Step 4: Commit**

```bash
git add config/schema_config.yaml
git commit -m "feat(schema): Merops_family_has_pfam_domain edge + cleavage properties on merops family"
```

---

### Task 6: Adapter — bridge edges + cleavage properties + VOCAB literals

**Files:**
- Modify: `multiomics_kg/adapters/merops_adapter.py`
- Modify: `create_knowledge_graph.py` (one arg)
- Test: `tests/test_merops_adapter.py` (append)

**Interfaces:**
- Consumes: reference `pfam_bridge`/`cleavage` blocks (Task 2); vocabulary
  entries (Task 4); `pfam_adapter.all_pfam_ids() -> set[str]` (existing —
  returns emitted Pfam node ids, `pfam:PF*` CURIEs).
- Produces: `MultiMeropsAnnotationAdapter(..., pfam_node_ids: set[str] | None = None)`;
  bridge edge tuples `(f"{fam}-has_pfam-{pf}", <merops node id>, <pfam node id>,
  "merops_family_has_pfam_domain", {"member_id_count": int})`.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_merops_adapter.py`)

Extend `_write_reference` in place — add to its `ref` dict:

```python
        "pfam_bridge": {
            "S14": {"PF00574": 120},
            "S08": {"PF00082": 300},
            "M99": {"PF99999": 5},      # family never observed -> no edge
        },
        "cleavage": {
            "S14": {
                "cleavage_p1_residues": ["Met", "Leu"],
                "cleavage_summary": "cleaves after Met (40%) / Leu (20%) - 30 known cleavages (50% physiological)",
                "known_cleavage_count": 30,
            },
        },
```

then append tests:

```python
def _bridge_edges(m):
    return [e for e in m.get_edges() if e[3] == "merops_family_has_pfam_domain"]


def test_bridge_edges_emitted_for_kept_families(tmp_path):
    # rebuild the standard fixture but with pfam ids injected
    genes = {"LT001": {"protein_id": "WP_1.1", "merops_ids": ["S14"]}}
    calls = {"WP_1.1": {"calls": [_candidate()]}}
    genome_dir = _write_strain(tmp_path, genes, calls)
    config = _write_genome_config(tmp_path, [genome_dir])
    m = MultiMeropsAnnotationAdapter(
        genome_config_file=str(config),
        reference_path=_write_reference(tmp_path),
        pfam_node_ids={"pfam:PF00574", "pfam:PF00082"},
    )
    m.download_data()
    edges = _bridge_edges(m)
    assert edges == [(
        "S14-has_pfam-PF00574", "merops.family:S14", "pfam:PF00574",
        "merops_family_has_pfam_domain", {"member_id_count": 120},
    )]  # S08 not observed in this fixture; M99 never observed


def test_bridge_pruned_by_pfam_node_ids(tmp_path):
    genes = {"LT001": {"protein_id": "WP_1.1", "merops_ids": ["S14"]}}
    calls = {"WP_1.1": {"calls": [_candidate()]}}
    genome_dir = _write_strain(tmp_path, genes, calls)
    config = _write_genome_config(tmp_path, [genome_dir])
    m = MultiMeropsAnnotationAdapter(
        genome_config_file=str(config),
        reference_path=_write_reference(tmp_path),
        pfam_node_ids={"pfam:PF00082"},          # PF00574 absent from graph
    )
    m.download_data()
    assert _bridge_edges(m) == []


def test_bridge_absent_without_injection(multi):
    """pfam_node_ids=None (default) -> no bridge edges (dangling-proof)."""
    assert _bridge_edges(multi) == []


def test_cleavage_properties_on_family_nodes(multi):
    nodes = {nid: props for nid, _, props in multi.get_nodes()}
    s14 = nodes["merops.family:S14"]
    assert s14["cleavage_p1_residues"] == ["Met", "Leu"]
    assert s14["known_cleavage_count"] == 30
    assert "known cleavages" in s14["cleavage_summary"]
    # sparse discipline: families without data carry none of the three keys
    s08 = nodes["merops.family:S08"]
    assert "cleavage_p1_residues" not in s08 and "cleavage_summary" not in s08
    assert "known_cleavage_count" not in s08
    # clans/subfamilies never carry them
    assert "cleavage_summary" not in nodes["merops.clan:SK"]
```

(The existing `multi` fixture gains the new reference blocks automatically via
the `_write_reference` edit; its node/edge assertions must keep passing —
`test_node_properties` compares the full S14 dict by equality, so extend its
expected dict with the three fixture values:
`"cleavage_p1_residues": ["Met", "Leu"]`,
`"cleavage_summary": "cleaves after Met (40%) / Leu (20%) - 30 known cleavages (50% physiological)"`,
`"known_cleavage_count": 30`.)

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_merops_adapter.py -q`
Expected: new tests FAIL (`pfam_node_ids` unexpected kwarg / missing props);
`test_node_properties` FAILS until the implementation lands.

- [ ] **Step 3: Implement in `merops_adapter.py`**

1. Constructor: add `pfam_node_ids: set[str] | None = None` keyword arg;
   store `self.pfam_node_ids = pfam_node_ids`.
2. Add module helper (next to `_merops_node_id`):

```python
def _pfam_node_id(acc: str) -> str:
    return normalize_curie(f"pfam:{acc}") or f"pfam_{acc}"
```

3. In `_node_props`, level-1 branch, after the name/family_type lines:

```python
            for key, val in (ref.get("cleavage", {}).get(code) or {}).items():
                if key == "cleavage_summary":
                    props[key] = _clean_str(val)
                elif key == "cleavage_p1_residues":
                    props[key] = [_clean_str(v) for v in val]
                else:
                    props[key] = val
```

4. In `get_edges`, between the parent-edge loop and the per-strain delegation,
   add (kept level-1 families only; sound-direction bridge):

```python
        # 2. Family→Pfam bridge (MEROPS interpro.txt, family-level only).
        # Dangling-proof: emitted only for injected, existing Pfam node ids
        # (TCDB-bridge precedent) — pfam_node_ids=None -> no bridge edges.
        bridge_count = 0
        if self.pfam_node_ids is not None:
            for fam, pfams in sorted(self._ref().get("pfam_bridge", {}).items()):
                if fam not in all_codes:
                    continue
                for pf, n in sorted(pfams.items()):
                    pf_id = _pfam_node_id(pf)
                    if pf_id not in self.pfam_node_ids:
                        continue
                    yield (
                        f"{fam}-has_pfam-{pf}",
                        _merops_node_id(fam),
                        pf_id,
                        "merops_family_has_pfam_domain",
                        {"member_id_count": n},
                    )
                    bridge_count += 1
```

   (rename the local `all_codes = self._expand_with_ancestors(...)` variable if
   it is currently named differently — it is `all_codes` today) and extend the
   final log line with `{bridge_count} bridge`.
5. Route the two categorical edge literals through the contract in the
   per-strain `get_edges` (import `from multiomics_kg.utils.controlled_vocab
   import VOCAB` at the top):

```python
                props: dict = {
                    "call_class": VOCAB.check(
                        "Gene_has_merops_family", "call_class", call_class(cand)),
                    "best_hit_id": _clean_str(cand.get("best_hit_id")),
                    "best_hit_kind": VOCAB.check(
                        "Gene_has_merops_family", "best_hit_kind",
                        cand.get("best_hit_kind")),
                }
```

- [ ] **Step 4: Wire the injection in `create_knowledge_graph.py`**

In the merops block, add the arg (pfam_adapter is defined earlier in `main`):

```python
    merops_adapter = MultiMeropsAnnotationAdapter(
        genome_config_file='data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
        pfam_node_ids=pfam_adapter.all_pfam_ids(),
        test_mode=TEST_MODE,
    )
```

and update the block comment: bridge edges + cleavage properties, Pfam node
set injected so the bridge can never dangle.

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/test_merops_adapter.py tests/test_create_knowledge_graph.py -q`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/adapters/merops_adapter.py create_knowledge_graph.py tests/test_merops_adapter.py
git commit -m "feat(merops-adapter): Pfam bridge edges + cleavage node properties + VOCAB-checked literals"
```

---

### Task 7: Post-import — `pfam_support` + fulltext extension

**Files:**
- Modify: `scripts/post-import.sh`
- Modify: `scripts/post-import.cypher` (identical edits)

**Interfaces:**
- Consumes: `Merops_family_has_pfam_domain` edges (Task 6);
  `Merops_family_is_a_merops_family` hierarchy (existing).
- Produces: `Gene_has_merops_family.pfam_support` on every merops gene edge.

- [ ] **Step 1: Replace the merops fulltext index with the drop+recreate form**

In BOTH files, replace:

```
CREATE FULLTEXT INDEX meropsFamilyFullText IF NOT EXISTS
    FOR (m:MeropsFamily) ON EACH [m.name, m.merops_id, m.description];
```

with:

```
// Full-text defs can't be ALTERed — drop + recreate so cleavage_summary is
// picked up even on reruns against an existing graph.
DROP INDEX meropsFamilyFullText IF EXISTS;
CREATE FULLTEXT INDEX meropsFamilyFullText IF NOT EXISTS
    FOR (m:MeropsFamily) ON EACH [m.name, m.merops_id, m.description, m.cleavage_summary];
```

- [ ] **Step 2: Add the `pfam_support` block**

In BOTH files, directly after the MeropsFamily computed-properties block, add:

```
// Gene_has_merops_family.pfam_support: whether a Pfam on this gene is curated
// into this MEROPS family via Merops_family_has_pfam_domain (interpro.txt).
// Sound direction only — corroborates a known family; never assigns peptidase
// identity backward from a domain. R5 string pair, tcdb pfam_support pattern.
// The bridge attaches at family level (level 1); subfamily-attached gene edges
// check their parent family (single is-a hop).
CALL {
  MATCH (g:Gene)-[r:Gene_has_merops_family]->(m:MeropsFamily)
  WITH r, g, m,
       EXISTS {
         MATCH (g)-[:Gene_has_pfam]->(:Pfam)<-[:Merops_family_has_pfam_domain]-(fam:MeropsFamily)
         WHERE fam = m OR (m)-[:Merops_family_is_a_merops_family]->(fam)
       } AS pfam_ok
  SET r.pfam_support = CASE WHEN pfam_ok THEN 'corroborated' ELSE 'uncorroborated' END
} IN TRANSACTIONS OF 1000 ROWS;
```

- [ ] **Step 3: Verify the two files agree and gates hold**

```bash
diff <(grep -i "merops\|pfam_support" scripts/post-import.sh) <(grep -i "merops\|pfam_support" scripts/post-import.cypher) && echo IDENTICAL
uv run pytest tests/test_annotation_quality_buckets.py -q
```
Expected: IDENTICAL; 5 bucket tests PASS (merops still bucket-free).

- [ ] **Step 4: Commit**

```bash
git add scripts/post-import.sh scripts/post-import.cypher
git commit -m "feat(post-import): Gene_has_merops_family.pfam_support + cleavage_summary in meropsFamilyFullText"
```

---

### Task 8: kg-validity additions (run post-rebuild)

**Files:**
- Modify: `tests/kg_validity/test_merops.py` (append)

**Interfaces:**
- Consumes: the rebuilt graph (these tests auto-skip while Neo4j is
  unreachable and are expected to pass only after the Docker rebuild).

- [ ] **Step 1: Append the tests**

```python
# ── Pfam bridge + pfam_support + cleavage (2026-08-18 follow-up) ─────────────

def test_bridge_edge_count_in_range(run_query):
    n = run_query("MATCH ()-[r:Merops_family_has_pfam_domain]->() RETURN count(r) AS n")[0]["n"]
    assert 150 <= n <= 200, f"expected ~183 bridge edges, got {n}"


def test_bridge_no_dangling_ends(run_query):
    n = run_query("""
        MATCH (a)-[r:Merops_family_has_pfam_domain]->(b)
        WHERE NOT a:MeropsFamily OR NOT b:Pfam
        RETURN count(r) AS n
    """)[0]["n"]
    assert n == 0


def test_bridge_family_level_only_with_support(run_query):
    n = run_query("""
        MATCH (m:MeropsFamily)-[r:Merops_family_has_pfam_domain]->()
        WHERE m.level <> 1 OR coalesce(r.member_id_count, 0) < 1
        RETURN count(r) AS n
    """)[0]["n"]
    assert n == 0


def test_s14_bridges_to_clpp_domain(run_query):
    row = run_query("""
        MATCH (:MeropsFamily {merops_id: 'S14'})-[:Merops_family_has_pfam_domain]->(p:Pfam)
        RETURN p.id AS id
    """)
    assert [r["id"] for r in row] == ["pfam:PF00574"]


def test_pfam_support_vocabulary_and_consistency(run_query):
    rows = run_query("""
        MATCH ()-[r:Gene_has_merops_family]->()
        RETURN DISTINCT r.pfam_support AS v
    """)
    assert {r["v"] for r in rows} == {"corroborated", "uncorroborated"}
    n = run_query("""
        MATCH (g:Gene)-[r:Gene_has_merops_family]->(m:MeropsFamily)
        WHERE r.pfam_support = 'corroborated'
          AND NOT EXISTS {
            MATCH (g)-[:Gene_has_pfam]->(:Pfam)<-[:Merops_family_has_pfam_domain]-(fam:MeropsFamily)
            WHERE fam = m OR (m)-[:Merops_family_is_a_merops_family]->(fam)
          }
        RETURN count(r) AS n
    """)[0]["n"]
    assert n == 0


def test_cleavage_properties_sparse_and_sane(run_query):
    n = run_query("""
        MATCH (m:MeropsFamily)
        WHERE (m.cleavage_summary = '' OR m.cleavage_p1_residues = [])
           OR (m.cleavage_summary IS NOT NULL AND m.level <> 1)
           OR (m.cleavage_summary IS NOT NULL AND m.known_cleavage_count IS NULL)
        RETURN count(m) AS n
    """)[0]["n"]
    assert n == 0
    covered = run_query("""
        MATCH (m:MeropsFamily) WHERE m.cleavage_summary IS NOT NULL
        RETURN count(m) AS n
    """)[0]["n"]
    assert 60 <= covered <= 88


def test_s01_cleaves_after_basic_residues(run_query):
    row = run_query("""
        MATCH (m:MeropsFamily {merops_id: 'S01'})
        RETURN m.cleavage_p1_residues AS p1, m.cleavage_summary AS s,
               m.known_cleavage_count AS n
    """)
    assert row and "Lys" in row[0]["p1"] and "Arg" in row[0]["p1"]
    assert "known cleavages" in row[0]["s"]
    assert row[0]["n"] > 10000
```

- [ ] **Step 2: Confirm they skip cleanly pre-rebuild / fail honestly**

Run: `uv run pytest tests/kg_validity/test_merops.py -q`
Expected against the CURRENT (pre-rebuild) graph: the new tests FAIL (no
bridge edges yet) — do NOT fix; they gate the rebuild. If Neo4j is down they
skip. Either outcome is acceptable at this step; note which occurred.

- [ ] **Step 3: Commit**

```bash
git add tests/kg_validity/test_merops.py
git commit -m "test(kg-validity): merops Pfam bridge, pfam_support, cleavage assertions (post-rebuild gate)"
```

---

### Task 9: Docs — kg-changes, CLAUDE.md, CHANGELOG, backlog

**Files:**
- Modify: `docs/kg-changes/merops-extension.md`
- Modify: `CLAUDE.md`
- Modify: `CHANGELOG.md`
- Modify: `plans/backlog.md`

- [ ] **Step 1: `docs/kg-changes/merops-extension.md`** — append a dated
  section "2026-08-18 — Pfam bridge + pfam_support + cleavage specificity"
  covering: the new edge (count, `member_id_count`, verb rationale,
  read-direction warning), `pfam_support` (R5 pair, definition, post-import),
  the three cleavage properties (sparse, family-level, honesty caveats:
  support spread, member blending, non-independent events), vocabulary
  registration (10 entries incl. the 7-entry gap-fix), the measured GO
  rejection (338 genes vs all-kingdom noise), step 9→10 + `--rebuild`, and a
  placeholder line for the measured `pfam_support` split to fill in
  post-rebuild.
- [ ] **Step 2: `CLAUDE.md`** — update: MeropsFamily key-facts bullet (bridge
  edge + `pfam_support` + cleavage props + vocab registration), the
  relationships list (`Merops_family_has_pfam_domain` after
  `Merops_family_is_a_merops_family`), merops_adapter bullet, step 9 / new
  step 10 text + Data Locations (reference now includes `pfam_bridge` +
  `cleavage`; raw interpro.txt/Substrate_search.txt gitignored), the
  prepare_data usage block (`--rebuild`), `meropsFamilyFullText` field list.
- [ ] **Step 3: `CHANGELOG.md`** — extend the `[Unreleased]` MEROPS Added
  bullet with: Pfam bridge (~183 edges), `pfam_support`, cleavage
  specificity, vocabulary registration, step-10 move + `--rebuild`.
- [ ] **Step 4: `plans/backlog.md`** — delete the two bullets this work lands
  ("MEROPS cross-ontology bridges", "MEROPS cleavage specificity as node
  properties"); add one new bullet recording the measured GO-arm rejection
  (so the reasoning survives the bullet removal):

```markdown
- [ ] **MEROPS GO bridge — rejected on measurement (2026-08-18), revisit only
      with a concrete use case.** All-kingdom member rollup yields median 19 /
      max 389 GO terms per family incl. eukaryote-only terms; completeness win
      is 338 genes (~8%) vs the 1,311 that justified TCDB's GO bridge. Needs a
      filtering design (e.g. >= N supporting identifiers) before it can land.
      → `docs/superpowers/specs/2026-08-18-merops-pfam-bridge-cleavage-design.md`
```

- [ ] **Step 5: Commit**

```bash
git add docs/kg-changes/merops-extension.md CLAUDE.md CHANGELOG.md plans/backlog.md
git commit -m "docs(merops): Pfam bridge + cleavage follow-up — kg-changes, CLAUDE.md, CHANGELOG, backlog"
```

---

### Task 10: Validation gate + rebuild hand-off

**Files:** none new (verification only).

- [ ] **Step 1: Full fast suite**

Run: `uv run pytest -m "not slow and not kg" -q`
Expected: all PASS.

- [ ] **Step 2: Test-mode KG build**

Run: `uv run python create_knowledge_graph.py --test --output-dir /tmp/kg-test-bridge`
(BioCypher writes to `biocypher-out/<timestamp>/`.) Then verify:

```bash
ls -t biocypher-out | head -1                        # newest dir
grep -c "" biocypher-out/<newest>/Merops_family_has_pfam_domain-part000.csv   # > 0 rows
head -2 biocypher-out/<newest>/MeropsFamily-part000.csv                        # cleavage columns populated for a family row
```
Expected: bridge CSV present with >0 rows; MeropsFamily header includes the
three cleavage columns. (If GitHub rate-limits the Biolink fetch with
404/429, wait and retry — known environment flake, not a code failure.)

- [ ] **Step 3: Push**

```bash
git push git@github.com:wosnat/multiomics_biocypher_kg.git main
```

- [ ] **Step 4: Hand off the Docker rebuild** (two-clone workflow: user pulls
  and rebuilds in `multiomics_biocypher_kg`). Post-rebuild, run in order:

```bash
uv run python .claude/skills/omics-edge-snapshot/omics_edge_snapshot.py --compare pre_merops_integration   # or a fresh pre-bridge baseline captured before the rebuild
uv run pytest -m kg -q
uv run python tests/kg_validity/generate_snapshot.py && uv run pytest tests/kg_validity/test_snapshot.py -q
```

Expected: snapshot compare exit 0 (expression/DM/metabolism unchanged); kg
suite all PASS incl. Task 8's tests; snapshot regenerated. Then: fill the
measured `pfam_support` split into `docs/kg-changes/merops-extension.md`
(query: corroborated share by `call_class` — expect high for `peptidase`,
depressed for `nonpeptidase_homolog`), commit, push.
```
