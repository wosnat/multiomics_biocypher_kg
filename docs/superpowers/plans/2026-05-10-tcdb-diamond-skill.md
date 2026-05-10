# TCDB-Diamond Skill Implementation Plan (Phase 1)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a `/tcdb-diamond` skill that runs diamond blastp against the curated TCDB FASTA per strain, applies a tiered confidence policy + lightweight post-steps (consensus collapse, eggNOG agreement, class-9 tag), and writes per-strain `<strain>.tcdb.calls.json` artifacts. Phase 1 only — no merge into `gene_annotations_merged.json` or KG.

**Architecture:** Mirror `/eggnog-run`: opt-in heavy compute, per-strain artifacts in the strain cache, no KG-side coupling. Reuse Saier Lab's `extractTCDB.pl` (one-time `git clone` of TCDBtools to `~/tools/TCDB/TCDBtools/`) for the FASTA + diamond DB build. Tier policy and post-steps are pure Python in a separate utility module so they're unit-testable without filesystem or subprocess.

**Tech Stack:** Python 3.11 (project), `diamond` CLI (already in env), `perl` + `wget` (system tools, used by extractTCDB.pl), `python-dotenv` (TCDB_DATA_DIR override). Tests use `pytest` (project standard, marker `not slow and not kg`).

**Spec:** [docs/superpowers/specs/2026-05-10-tcdb-diamond-augmentation-design.md](../specs/2026-05-10-tcdb-diamond-augmentation-design.md)

---

## File Structure

| Path | Type | Responsibility |
|---|---|---|
| `multiomics_kg/utils/tcdb_diamond.py` | New module | Pure-Python tier policy, consensus collapse, eggNOG agreement, parsing helpers — testable without filesystem |
| `tests/test_tcdb_diamond.py` | New tests | Unit tests for `tcdb_diamond` module |
| `.claude/skills/tcdb-diamond/SKILL.md` | New | Skill manifest + workflow doc (mirrors `/eggnog-run`) |
| `.claude/skills/tcdb-diamond/run_tcdb_diamond.py` | New | Orchestrator: registry loop, DB-build, diamond invocation, post-step driver |
| `~/tools/TCDB/DB/` | New (created at runtime) | Skill-managed download cache: `tcdb.faa`, `tcdb.dmnd`, `tcdb_acc2tcid.tsv`, `refresh.log` |
| `cache/<org>/genomes/<strain>/tcdb/` | New (created at runtime) | Per-strain output: `<strain>.tcdb.tsv`, `<strain>.tcdb.calls.json`, `<strain>.tcdb.skill_summary.json` |

The split between `tcdb_diamond.py` (logic) and `run_tcdb_diamond.py` (filesystem + subprocess) keeps the testable surface large (logic) and the untestable surface small (orchestration).

`multiomics_kg/utils/tcdb_utils.py` (already exists) is **not modified** — it provides `tcdb_ancestors()` which Phase 2 will use; Phase 1 needs only simple TCID truncation, not full hierarchy walks.

---

## Task 1: Skill scaffold + empty utility module

**Files:**
- Create: `multiomics_kg/utils/tcdb_diamond.py`
- Create: `tests/test_tcdb_diamond.py`
- Create: `.claude/skills/tcdb-diamond/SKILL.md` (placeholder, filled in Task 14)
- Create: `.claude/skills/tcdb-diamond/run_tcdb_diamond.py` (placeholder, filled in Task 13)

- [ ] **Step 1: Create the empty utility module**

```python
# multiomics_kg/utils/tcdb_diamond.py
"""TCDB-vs-Diamond per-hit tier policy + per-protein post-steps.

Pure Python — no filesystem or subprocess. The orchestrator in
`.claude/skills/tcdb-diamond/run_tcdb_diamond.py` is responsible for I/O.
"""
from __future__ import annotations
```

- [ ] **Step 2: Create the empty test file**

```python
# tests/test_tcdb_diamond.py
"""Unit tests for multiomics_kg.utils.tcdb_diamond."""
from multiomics_kg.utils import tcdb_diamond  # noqa: F401
```

- [ ] **Step 3: Create skill placeholder files**

```bash
mkdir -p .claude/skills/tcdb-diamond
```

```markdown
<!-- .claude/skills/tcdb-diamond/SKILL.md -->
---
name: tcdb-diamond
description: Placeholder. Filled in Task 14.
user-invocable: true
---
```

```python
# .claude/skills/tcdb-diamond/run_tcdb_diamond.py
"""Placeholder. Filled in Task 13."""
```

- [ ] **Step 4: Verify import works**

Run: `uv run python -c "from multiomics_kg.utils import tcdb_diamond; import tests.test_tcdb_diamond"`
Expected: no output, no error.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/tcdb_diamond.py tests/test_tcdb_diamond.py .claude/skills/tcdb-diamond/
git commit -m "scaffold(tcdb-diamond): empty skill + utility module + test file"
```

---

## Task 2: TCID truncation helper

**Files:**
- Modify: `multiomics_kg/utils/tcdb_diamond.py`
- Modify: `tests/test_tcdb_diamond.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_tcdb_diamond.py — append
from multiomics_kg.utils.tcdb_diamond import truncate_tcid


def test_truncate_tcid_keeps_first_n_parts():
    assert truncate_tcid("1.A.11.1.5", 3) == "1.A.11"
    assert truncate_tcid("1.A.11.1.5", 4) == "1.A.11.1"
    assert truncate_tcid("1.A.11.1.5", 5) == "1.A.11.1.5"


def test_truncate_tcid_passthrough_when_already_short():
    assert truncate_tcid("1.A.11", 5) == "1.A.11"
    assert truncate_tcid("1.A", 3) == "1.A"


def test_truncate_tcid_empty_input_returns_empty():
    assert truncate_tcid("", 3) == ""
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_tcdb_diamond.py -v`
Expected: FAIL with `ImportError` or `AttributeError` on `truncate_tcid`.

- [ ] **Step 3: Implement**

```python
# multiomics_kg/utils/tcdb_diamond.py — append
def truncate_tcid(tcid: str, n_parts: int) -> str:
    """Truncate a TCID to its first n_parts dot-separated segments.

    A TCID has up to 5 parts (class.subclass.family.subfamily.specificity).
    Returns the input unchanged if it already has <= n_parts segments.
    """
    if not tcid:
        return ""
    parts = tcid.split(".")
    return ".".join(parts[:n_parts])
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_tcdb_diamond.py -v`
Expected: 3 passed.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/tcdb_diamond.py tests/test_tcdb_diamond.py
git commit -m "feat(tcdb-diamond): add truncate_tcid helper"
```

---

## Task 3: Per-hit tier classifier

The tier policy from spec §6.3 is encoded as a single function operating on parsed-row dicts. Returns the tier number (1, 2, 3) or `None` to drop the hit.

**Files:**
- Modify: `multiomics_kg/utils/tcdb_diamond.py`
- Modify: `tests/test_tcdb_diamond.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_tcdb_diamond.py — append
from multiomics_kg.utils.tcdb_diamond import classify_hit


def _hit(identity, qcov, scov, length=200, evalue=1e-10):
    """Build a hit dict with sensible defaults — only override what the test cares about."""
    return {
        "identity": identity, "qcov": qcov, "scov": scov,
        "length": length, "evalue": evalue,
    }


def test_classify_tier_1_high_identity_high_qcov():
    # 80%/85% identity over 80%/75% qcov — tier 1
    assert classify_hit(_hit(80.0, 85.0, 50.0)) == 1
    assert classify_hit(_hit(70.0, 70.0, 50.0)) == 1


def test_classify_tier_2_mid_identity_mid_qcov():
    # 50% identity / 65% qcov — passes tier 2 floor, fails tier 1
    assert classify_hit(_hit(50.0, 65.0, 30.0)) == 2
    assert classify_hit(_hit(40.0, 60.0, 30.0)) == 2


def test_classify_tier_3_gblast3_floor():
    # 35% identity (no tier 1/2), qcov ≥ 40 — tier 3
    assert classify_hit(_hit(35.0, 45.0, 30.0)) == 3
    # qcov < 40 BUT scov ≥ 40 — still tier 3 (gblast3 OR rule)
    assert classify_hit(_hit(35.0, 25.0, 50.0)) == 3
    # 25% identity is fine for tier 3 (no identity floor)
    assert classify_hit(_hit(25.0, 45.0, 30.0)) == 3


def test_classify_drops_hit_below_floor():
    # qcov AND scov both < 40 — fails gblast3 OR rule
    assert classify_hit(_hit(35.0, 30.0, 30.0)) is None
    # length < 50 — fails HSP-length floor
    assert classify_hit(_hit(80.0, 90.0, 90.0, length=49)) is None
    # e-value > 0.001 — fails e-value gate
    assert classify_hit(_hit(80.0, 90.0, 90.0, evalue=0.01)) is None


def test_classify_handles_boundary_values_inclusively():
    # All thresholds are >= (not >) — boundary inputs should pass
    assert classify_hit(_hit(70.0, 70.0, 0.0)) == 1
    assert classify_hit(_hit(40.0, 60.0, 0.0)) == 2
    assert classify_hit(_hit(0.0, 40.0, 0.0)) == 3
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_tcdb_diamond.py::test_classify_tier_1_high_identity_high_qcov -v`
Expected: FAIL with `ImportError` on `classify_hit`.

- [ ] **Step 3: Implement**

```python
# multiomics_kg/utils/tcdb_diamond.py — append
def classify_hit(hit: dict) -> int | None:
    """Assign a confidence tier (1/2/3) to a parsed diamond hit row.

    Returns None when the hit fails the tier-3 floor (drop it).

    Floor (gblast3-style): e-value <= 0.001, HSP length >= 50, AND
    (qcov >= 40 OR scov >= 40). Above the floor:
      - tier 1: identity >= 70 AND qcov >= 70
      - tier 2: identity >= 40 AND qcov >= 60
      - tier 3: floor only
    """
    if hit["evalue"] > 0.001:
        return None
    if hit["length"] < 50:
        return None
    if hit["qcov"] < 40.0 and hit["scov"] < 40.0:
        return None

    if hit["identity"] >= 70.0 and hit["qcov"] >= 70.0:
        return 1
    if hit["identity"] >= 40.0 and hit["qcov"] >= 60.0:
        return 2
    return 3
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_tcdb_diamond.py -v`
Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/tcdb_diamond.py tests/test_tcdb_diamond.py
git commit -m "feat(tcdb-diamond): add classify_hit tier policy"
```

---

## Task 4: Consensus collapse across top-N hits per protein

Spec §6.4-A: with up to 5 hits per query, the per-protein call collapses to the deepest level at which **all** top-N hits share the same prefix.

**Files:**
- Modify: `multiomics_kg/utils/tcdb_diamond.py`
- Modify: `tests/test_tcdb_diamond.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_tcdb_diamond.py — append
from multiomics_kg.utils.tcdb_diamond import consensus_collapse


def test_consensus_all_agree_at_5_part():
    hits = [{"tcid": "1.A.11.1.5"}, {"tcid": "1.A.11.1.5"}, {"tcid": "1.A.11.1.5"}]
    result = consensus_collapse(hits)
    assert result == {"tcid": "1.A.11.1.5", "agreement": "5_part", "n": 3}


def test_consensus_demote_to_4_part():
    # 4-part agreement, disagreement at position 5
    hits = [{"tcid": "1.A.11.1.5"}, {"tcid": "1.A.11.1.7"}]
    result = consensus_collapse(hits)
    assert result == {"tcid": "1.A.11.1", "agreement": "4_part", "n": 2}


def test_consensus_demote_to_3_part():
    # 3-part agreement only
    hits = [{"tcid": "1.A.11.1.5"}, {"tcid": "1.A.11.2.3"}]
    result = consensus_collapse(hits)
    assert result == {"tcid": "1.A.11", "agreement": "3_part", "n": 2}


def test_consensus_reject_below_3_part_agreement():
    # Different families
    hits = [{"tcid": "1.A.11.1.5"}, {"tcid": "2.A.7.4.3"}]
    result = consensus_collapse(hits)
    assert result is None


def test_consensus_single_hit_keeps_5_part():
    hits = [{"tcid": "1.A.11.1.5"}]
    result = consensus_collapse(hits)
    assert result == {"tcid": "1.A.11.1.5", "agreement": "5_part", "n": 1}


def test_consensus_empty_input_returns_none():
    assert consensus_collapse([]) is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_tcdb_diamond.py -v -k consensus`
Expected: FAIL with `ImportError`.

- [ ] **Step 3: Implement**

```python
# multiomics_kg/utils/tcdb_diamond.py — append
def consensus_collapse(hits: list[dict]) -> dict | None:
    """Collapse a list of hits (one query, top-N TCDB hits) into a per-protein call.

    Returns None when fewer than 3 leading parts agree across all hits.
    The returned dict has keys: tcid (str), agreement ("5_part" | "4_part" | "3_part"),
    n (number of hits considered).

    All hits are assumed to carry 5-part TCIDs (TCDB curates only at the
    tc_specificity leaves), so dot-prefix comparison is well-defined.
    """
    if not hits:
        return None
    parts_lists = [h["tcid"].split(".") for h in hits]

    for depth in (5, 4, 3):
        prefixes = {tuple(p[:depth]) for p in parts_lists}
        if len(prefixes) == 1:
            shared = parts_lists[0][:depth]
            return {
                "tcid": ".".join(shared),
                "agreement": f"{depth}_part",
                "n": len(hits),
            }
    return None
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_tcdb_diamond.py -v`
Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/tcdb_diamond.py tests/test_tcdb_diamond.py
git commit -m "feat(tcdb-diamond): add consensus_collapse for per-protein call"
```

---

## Task 5: eggNOG agreement tag

Spec §6.4-B: classify the diamond call's relationship to eggNOG's `KEGG_TC` value for the same protein.

**Files:**
- Modify: `multiomics_kg/utils/tcdb_diamond.py`
- Modify: `tests/test_tcdb_diamond.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_tcdb_diamond.py — append
from multiomics_kg.utils.tcdb_diamond import compute_egn_agreement


def test_egn_agreement_confirms_identical():
    assert compute_egn_agreement("1.A.11", "1.A.11") == "confirms"
    assert compute_egn_agreement("1.A.11.1.5", "1.A.11.1.5") == "confirms"


def test_egn_agreement_confirms_diamond_descendant():
    # eggNOG family-level, diamond a strict descendant — still confirms (consistent)
    assert compute_egn_agreement("1.A.11.1.5", "1.A.11") == "refines"
    assert compute_egn_agreement("1.A.11.1", "1.A.11") == "refines"


def test_egn_agreement_extends_when_eggnog_missing():
    assert compute_egn_agreement("1.A.11.1.5", None) == "extends"
    assert compute_egn_agreement("1.A.11.1.5", "") == "extends"


def test_egn_agreement_conflicts_different_family():
    assert compute_egn_agreement("1.A.11.1.5", "2.A.7") == "conflicts"
    assert compute_egn_agreement("1.A.11", "1.B.5") == "conflicts"


def test_egn_agreement_eggnog_descendant_of_diamond_is_conflict():
    # Diamond at 3-part, eggNOG at 5-part below it. Same family → not a conflict;
    # this case is rare (eggNOG rarely emits 5-part) but treat as confirms.
    assert compute_egn_agreement("1.A.11", "1.A.11.1.5") == "confirms"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_tcdb_diamond.py -v -k egn_agreement`
Expected: FAIL with `ImportError`.

- [ ] **Step 3: Implement**

```python
# multiomics_kg/utils/tcdb_diamond.py — append
def compute_egn_agreement(diamond_tcid: str, egn_tcid: str | None) -> str:
    """Tag the relationship between the diamond call and eggNOG's KEGG_TC.

    Returns one of: "confirms" | "refines" | "extends" | "conflicts".
    `egn_only` (eggNOG TC present, diamond absent) is not produced here —
    those proteins simply don't appear in the calls JSON in Phase 1.

    Rules:
      - confirms: identical TCIDs OR one is a strict prefix of the other
        AT family level or below (first 3 parts match)
      - refines: eggNOG TCID is a strict prefix of diamond TCID — same lineage,
        diamond went deeper. Reported separately from confirms because this
        is the headline specificity win.
      - extends: eggNOG had no TC; diamond produced one
      - conflicts: family-level (first 3 parts) disagrees
    """
    if not egn_tcid:
        return "extends"
    if diamond_tcid == egn_tcid:
        return "confirms"

    diamond_parts = diamond_tcid.split(".")
    egn_parts = egn_tcid.split(".")

    # Family-level disagreement → conflict
    if diamond_parts[:3] != egn_parts[:3]:
        return "conflicts"

    # Same family. Diamond strictly deeper than eggNOG → refines.
    # eggNOG strictly deeper than diamond → confirms (rare).
    if len(diamond_parts) > len(egn_parts) and diamond_parts[: len(egn_parts)] == egn_parts:
        return "refines"
    return "confirms"
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_tcdb_diamond.py -v`
Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/tcdb_diamond.py tests/test_tcdb_diamond.py
git commit -m "feat(tcdb-diamond): add compute_egn_agreement tag"
```

---

## Task 6: Class-9 tag

Spec §6.4-C: TCDB class `9.*` = "Incompletely Characterized Transport Systems". Tagged but not demoted.

**Files:**
- Modify: `multiomics_kg/utils/tcdb_diamond.py`
- Modify: `tests/test_tcdb_diamond.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_tcdb_diamond.py — append
from multiomics_kg.utils.tcdb_diamond import is_class_9


def test_is_class_9_matches_top_class():
    assert is_class_9("9.B.82.1.5") is True
    assert is_class_9("9.A.1") is True
    assert is_class_9("9") is True


def test_is_class_9_excludes_other_classes():
    assert is_class_9("1.A.11.1.5") is False
    assert is_class_9("8.A.1") is False  # auxiliary, not incompletely characterized
    assert is_class_9("19.A.1") is False  # not real but tests prefix matching


def test_is_class_9_handles_empty():
    assert is_class_9("") is False
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_tcdb_diamond.py -v -k class_9`
Expected: FAIL with `ImportError`.

- [ ] **Step 3: Implement**

```python
# multiomics_kg/utils/tcdb_diamond.py — append
def is_class_9(tcid: str) -> bool:
    """True iff TCID is in TCDB class 9 (Incompletely Characterized Transport Systems)."""
    if not tcid:
        return False
    return tcid.split(".", 1)[0] == "9"
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_tcdb_diamond.py -v`
Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/tcdb_diamond.py tests/test_tcdb_diamond.py
git commit -m "feat(tcdb-diamond): add is_class_9 tag"
```

---

## Task 7: Diamond TSV row parser

Diamond emits 8-column TSV rows. Parse one row into a typed dict.

**Files:**
- Modify: `multiomics_kg/utils/tcdb_diamond.py`
- Modify: `tests/test_tcdb_diamond.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_tcdb_diamond.py — append
from multiomics_kg.utils.tcdb_diamond import parse_diamond_row


def test_parse_diamond_row_extracts_typed_fields():
    line = "WP_011131900.1\tlcl|Q9I3F6-1.A.11.1.5\t87.4\t92.1\t89.7\t412\t1.2e-180\t650.5"
    row = parse_diamond_row(line)
    assert row["query_id"] == "WP_011131900.1"
    assert row["subject_id"] == "lcl|Q9I3F6-1.A.11.1.5"
    assert row["identity"] == 87.4
    assert row["qcov"] == 92.1
    assert row["scov"] == 89.7
    assert row["length"] == 412
    assert row["evalue"] == 1.2e-180
    assert row["bitscore"] == 650.5


def test_parse_diamond_row_returns_none_for_short_line():
    assert parse_diamond_row("only\ttwo\tcolumns") is None
    assert parse_diamond_row("") is None


def test_parse_diamond_row_returns_none_for_invalid_numeric():
    line = "WP_111.1\tlcl|X-1.A\tNOT_A_NUMBER\t90\t90\t100\t1e-10\t300"
    assert parse_diamond_row(line) is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_tcdb_diamond.py -v -k parse_diamond`
Expected: FAIL with `ImportError`.

- [ ] **Step 3: Implement**

```python
# multiomics_kg/utils/tcdb_diamond.py — append
def parse_diamond_row(line: str) -> dict | None:
    """Parse one diamond blastp output line (--outfmt 6, 8 columns) to a dict.

    Returns None when the row is malformed (wrong column count, non-numeric
    field). Caller is responsible for further filtering / classification.
    """
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 8:
        return None
    try:
        return {
            "query_id": parts[0],
            "subject_id": parts[1],
            "identity": float(parts[2]),
            "qcov": float(parts[3]),
            "scov": float(parts[4]),
            "length": int(parts[5]),
            "evalue": float(parts[6]),
            "bitscore": float(parts[7]),
        }
    except ValueError:
        return None
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_tcdb_diamond.py -v`
Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/tcdb_diamond.py tests/test_tcdb_diamond.py
git commit -m "feat(tcdb-diamond): add parse_diamond_row TSV parser"
```

---

## Task 8: TCDB FASTA header parser → acc2tcid map

Spec §6.1: parse the headers of `tcdb.faa` (downloaded by extractTCDB.pl) into an `(accession → tcid)` lookup.

The download URL embeds TCID in headers. From inspecting `extractTCDB.pl`, the format is `>lcl|<accession>-<5-part-tcid>` (e.g. `>lcl|Q9I3F6-1.A.11.1.5`). The diamond `subject_id` field is this whole `lcl|...` string, so we need to extract `(accession, tcid)` from it.

**Files:**
- Modify: `multiomics_kg/utils/tcdb_diamond.py`
- Modify: `tests/test_tcdb_diamond.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_tcdb_diamond.py — append
from multiomics_kg.utils.tcdb_diamond import parse_tcdb_subject_id


def test_parse_tcdb_subject_id_lcl_prefix():
    assert parse_tcdb_subject_id("lcl|Q9I3F6-1.A.11.1.5") == ("Q9I3F6", "1.A.11.1.5")


def test_parse_tcdb_subject_id_no_lcl_prefix():
    # If header was written without lcl|, accept it
    assert parse_tcdb_subject_id("Q9I3F6-1.A.11.1.5") == ("Q9I3F6", "1.A.11.1.5")


def test_parse_tcdb_subject_id_handles_dashes_in_accession():
    # UniProt accessions can contain dashes for isoforms — split on the LAST dash
    # before a dotted TCID
    result = parse_tcdb_subject_id("lcl|P12345-2-1.A.11.1.5")
    assert result == ("P12345-2", "1.A.11.1.5")


def test_parse_tcdb_subject_id_returns_none_when_no_tcid():
    assert parse_tcdb_subject_id("lcl|Q9I3F6") is None
    assert parse_tcdb_subject_id("") is None


def test_parse_tcdb_subject_id_validates_tcid_shape():
    # TCID must be at least 3 dotted parts to be plausible
    assert parse_tcdb_subject_id("lcl|Q9I3F6-1.A") is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_tcdb_diamond.py -v -k tcdb_subject`
Expected: FAIL with `ImportError`.

- [ ] **Step 3: Implement**

```python
# multiomics_kg/utils/tcdb_diamond.py — append
import re

_TCID_TAIL_RE = re.compile(r"-(\d+(?:\.[A-Za-z0-9]+){2,4})$")


def parse_tcdb_subject_id(subject_id: str) -> tuple[str, str] | None:
    """Extract (accession, tcid) from a TCDB FASTA-derived subject ID.

    Header format: ``[lcl|]<accession>-<TCID>`` where TCID is dot-separated
    with 3-5 parts (e.g. ``lcl|Q9I3F6-1.A.11.1.5``).

    Returns None when the subject ID does not contain a parseable TCID tail.
    Splits on the LAST dash followed by a dotted TCID — handles UniProt
    isoform accessions (e.g. ``P12345-2``) correctly.
    """
    if not subject_id:
        return None
    if subject_id.startswith("lcl|"):
        subject_id = subject_id[4:]
    match = _TCID_TAIL_RE.search(subject_id)
    if not match:
        return None
    tcid = match.group(1)
    accession = subject_id[: match.start()]
    if not accession:
        return None
    return accession, tcid
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_tcdb_diamond.py -v`
Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/tcdb_diamond.py tests/test_tcdb_diamond.py
git commit -m "feat(tcdb-diamond): add parse_tcdb_subject_id for FASTA headers"
```

---

## Task 9: eggNOG annotations reader

Read `<strain>.emapper.annotations` and return `{protein_id: kegg_tc}`. eggNOG's TSV uses `#` for comments and `#query` for the header; the `KEGG_TC` column is index 17.

**Files:**
- Modify: `multiomics_kg/utils/tcdb_diamond.py`
- Modify: `tests/test_tcdb_diamond.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_tcdb_diamond.py — append
import textwrap
from multiomics_kg.utils.tcdb_diamond import load_eggnog_kegg_tc


def test_load_eggnog_kegg_tc_extracts_per_protein_value(tmp_path):
    # Real eggNOG-mapper v2.1 column order, 21 columns
    content = textwrap.dedent("""\
        ## emapper-2.1.13
        #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
        WP_011131852.1\t59919.PMM0213\t2.76e-223\t617.0\tCOG3329@1|root\t1117\tS\tdesc\tsbtA\t-\t-\tko:K07086\t-\t-\t-\t-\tko00000\t1.A.11\t-\t-\tSbt_1
        WP_011131900.1\t59919.PMM0263\t0.0\t934.0\tCOG0004@1|root\t1117\tP\tAmm\tamtB\t-\t-\tko:K03320\t-\t-\t-\t-\tko00000,ko02000\t1.A.11\t-\t-\tAmmonium_transp
        WP_NOTC.1\t59919.PMM0001\t0.0\t100.0\tCOGZZZ\t1117\tS\tdesc\tx\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-
        ## end of file
        """)
    path = tmp_path / "MED4.emapper.annotations"
    path.write_text(content)

    result = load_eggnog_kegg_tc(path)
    assert result == {
        "WP_011131852.1": "1.A.11",
        "WP_011131900.1": "1.A.11",
        # WP_NOTC.1 has KEGG_TC = "-" → not in dict
    }


def test_load_eggnog_kegg_tc_handles_missing_file(tmp_path):
    assert load_eggnog_kegg_tc(tmp_path / "no_such_file") == {}


def test_load_eggnog_kegg_tc_skips_empty_or_dash_values(tmp_path):
    content = textwrap.dedent("""\
        ## emapper-2.1.13
        #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
        WP_A.1\tx\tx\tx\tx\tx\tS\tx\tx\t-\t-\tx\t-\t-\t-\t-\t-\t-\t-\t-\t-
        WP_B.1\tx\tx\tx\tx\tx\tS\tx\tx\t-\t-\tx\t-\t-\t-\t-\t-\t\t-\t-\t-
        WP_C.1\tx\tx\tx\tx\tx\tS\tx\tx\t-\t-\tx\t-\t-\t-\t-\t-\t1.A.11,3.A.1.27\t-\t-\t-
        """)
    path = tmp_path / "x.emapper.annotations"
    path.write_text(content)
    result = load_eggnog_kegg_tc(path)
    # Multi-value KEGG_TC: keep first one (rare; we don't try to merge)
    assert result == {"WP_C.1": "1.A.11"}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_tcdb_diamond.py -v -k eggnog_kegg_tc`
Expected: FAIL with `ImportError`.

- [ ] **Step 3: Implement**

```python
# multiomics_kg/utils/tcdb_diamond.py — append
from pathlib import Path


def load_eggnog_kegg_tc(annotations_path: Path) -> dict[str, str]:
    """Read an eggNOG-mapper .emapper.annotations file → {protein_id: kegg_tc}.

    Returns an empty dict when the file is absent. Skips rows whose KEGG_TC
    column is "-", empty, or missing. When KEGG_TC carries multiple values
    (rare; comma-separated), only the first is returned (we do not attempt
    to merge here — the merge happens in compute_egn_agreement at call time).
    """
    annotations_path = Path(annotations_path)
    if not annotations_path.exists():
        return {}

    result: dict[str, str] = {}
    KEGG_TC_COL = 17  # 0-indexed, post-emapper-v2.1 column order

    with open(annotations_path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                continue  # the #query header line
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= KEGG_TC_COL:
                continue
            protein_id = cols[0]
            kegg_tc_raw = cols[KEGG_TC_COL].strip()
            if not kegg_tc_raw or kegg_tc_raw == "-":
                continue
            # Multi-value: take first
            kegg_tc = kegg_tc_raw.split(",")[0].strip()
            if kegg_tc:
                result[protein_id] = kegg_tc
    return result
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_tcdb_diamond.py -v`
Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/tcdb_diamond.py tests/test_tcdb_diamond.py
git commit -m "feat(tcdb-diamond): add load_eggnog_kegg_tc reader"
```

---

## Task 10: Per-strain pipeline composer

Compose Tasks 3-9 into a single function that takes a diamond TSV path + an eggNOG annotations path and produces the full `calls.json` payload.

**Files:**
- Modify: `multiomics_kg/utils/tcdb_diamond.py`
- Modify: `tests/test_tcdb_diamond.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_tcdb_diamond.py — append
import textwrap
from multiomics_kg.utils.tcdb_diamond import build_strain_calls


def test_build_strain_calls_full_pipeline(tmp_path):
    tsv_content = (
        # 5 strong identical hits → tier 1 5-part call
        "WP_AAA.1\tlcl|Q1-1.A.11.1.5\t87.4\t92.1\t89.7\t412\t1e-180\t650\n"
        "WP_AAA.1\tlcl|Q2-1.A.11.1.5\t86.0\t91.0\t88.0\t410\t1e-179\t640\n"
        "WP_AAA.1\tlcl|Q3-1.A.11.1.5\t85.0\t90.0\t87.0\t408\t1e-178\t630\n"
        # Hits scattered across families → reject (consensus < 3-part)
        "WP_BBB.1\tlcl|Q4-1.A.11.1.5\t75.0\t75.0\t60.0\t300\t1e-100\t450\n"
        "WP_BBB.1\tlcl|Q5-2.A.7.4.3\t72.0\t73.0\t60.0\t300\t1e-95\t440\n"
        # Single hit, mid identity → tier 2 with 4-part TCID
        "WP_CCC.1\tlcl|Q6-1.A.11.1.5\t50.0\t65.0\t40.0\t250\t1e-50\t300\n"
        # Single hit, low identity passing tier 3 floor → tier 3 with 3-part
        "WP_DDD.1\tlcl|Q7-1.A.11.1.5\t30.0\t45.0\t30.0\t150\t1e-10\t120\n"
        # Class 9 hit → tagged
        "WP_EEE.1\tlcl|Q8-9.B.82.1.5\t80.0\t85.0\t75.0\t300\t1e-150\t500\n"
        # Hit failing the floor → dropped
        "WP_FFF.1\tlcl|Q9-1.A.11.1.5\t30.0\t30.0\t30.0\t150\t1e-10\t100\n"
    )
    tsv = tmp_path / "test.tcdb.tsv"
    tsv.write_text(tsv_content)

    egn_content = textwrap.dedent("""\
        ## emapper-2.1.13
        #query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
        WP_AAA.1\tx\tx\tx\tx\tx\tS\tx\tx\t-\t-\tx\t-\t-\t-\t-\t-\t1.A.11\t-\t-\t-
        WP_DDD.1\tx\tx\tx\tx\tx\tS\tx\tx\t-\t-\tx\t-\t-\t-\t-\t-\t9.A.1\t-\t-\t-
        """)
    egn = tmp_path / "test.emapper.annotations"
    egn.write_text(egn_content)

    calls, summary = build_strain_calls(tsv, egn)

    # WP_AAA.1: 3 strong consensus hits at 5-part → tier 1 specificity, refines eggNOG
    assert calls["WP_AAA.1"]["tcid"] == "1.A.11.1.5"
    assert calls["WP_AAA.1"]["tier"] == 1
    assert calls["WP_AAA.1"]["level_kind"] == "tc_specificity"
    assert calls["WP_AAA.1"]["consensus_agreement"] == "5_part"
    assert calls["WP_AAA.1"]["consensus_n"] == 3
    assert calls["WP_AAA.1"]["egn_agreement"] == "refines"
    assert calls["WP_AAA.1"]["egn_tcid"] == "1.A.11"
    assert calls["WP_AAA.1"]["incompletely_characterized"] is False

    # WP_BBB.1: scattered hits → rejected (not in calls)
    assert "WP_BBB.1" not in calls

    # WP_CCC.1: single tier-2 hit → 4-part subfamily
    assert calls["WP_CCC.1"]["tcid"] == "1.A.11.1"
    assert calls["WP_CCC.1"]["tier"] == 2
    assert calls["WP_CCC.1"]["level_kind"] == "tc_subfamily"
    assert calls["WP_CCC.1"]["egn_agreement"] == "extends"

    # WP_DDD.1: tier 3 → 3-part family; eggNOG conflict (different family)
    assert calls["WP_DDD.1"]["tcid"] == "1.A.11"
    assert calls["WP_DDD.1"]["tier"] == 3
    assert calls["WP_DDD.1"]["level_kind"] == "tc_family"
    assert calls["WP_DDD.1"]["egn_agreement"] == "conflicts"
    assert calls["WP_DDD.1"]["egn_tcid"] == "9.A.1"

    # WP_EEE.1: class 9 → tagged
    assert calls["WP_EEE.1"]["tcid"] == "9.B.82.1.5"
    assert calls["WP_EEE.1"]["incompletely_characterized"] is True
    assert calls["WP_EEE.1"]["egn_agreement"] == "extends"

    # WP_FFF.1: floor failure → not in calls
    assert "WP_FFF.1" not in calls

    # Summary
    assert summary["raw_hit_lines"] == 9
    assert summary["proteins_with_call"] == 4
    assert summary["proteins_rejected_by_consensus"] == 1
    assert summary["tier_distribution"] == {"1": 1, "2": 1, "3": 2}
    assert summary["agreement_distribution"]["refines"] == 1
    assert summary["agreement_distribution"]["extends"] == 2
    assert summary["agreement_distribution"]["conflicts"] == 1
    assert summary["agreement_distribution"]["confirms"] == 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_tcdb_diamond.py -v -k build_strain_calls`
Expected: FAIL with `ImportError`.

- [ ] **Step 3: Implement**

```python
# multiomics_kg/utils/tcdb_diamond.py — append
from collections import defaultdict


_TIER_TO_LEVEL_KIND = {
    1: "tc_specificity",
    2: "tc_subfamily",
    3: "tc_family",
}
_AGREEMENT_TO_PARTS = {"5_part": 5, "4_part": 4, "3_part": 3}


def build_strain_calls(
    tsv_path: Path,
    eggnog_annotations_path: Path,
) -> tuple[dict, dict]:
    """Run the full per-strain pipeline: parse TSV, classify, consensus, tag.

    Returns (calls, summary):
      calls: dict keyed by protein_id with the full §6.5 record shape
      summary: dict with raw counts (matches §6.6 stdout columns)
    """
    egn_lookup = load_eggnog_kegg_tc(Path(eggnog_annotations_path))

    # Group accepted (tier-classified) hits per query
    by_query: dict[str, list[dict]] = defaultdict(list)
    raw_lines = 0
    with open(tsv_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            raw_lines += 1
            row = parse_diamond_row(line)
            if row is None:
                continue
            tier = classify_hit(row)
            if tier is None:
                continue
            parsed = parse_tcdb_subject_id(row["subject_id"])
            if parsed is None:
                continue
            _, hit_tcid = parsed
            by_query[row["query_id"]].append({
                "tcid": hit_tcid,
                "tier": tier,
                "identity": row["identity"],
                "qcov": row["qcov"],
                "scov": row["scov"],
                "evalue": row["evalue"],
                "length": row["length"],
            })

    calls: dict[str, dict] = {}
    rejected = 0
    for query_id, hits in by_query.items():
        # Consensus collapse
        consensus = consensus_collapse(hits)
        if consensus is None:
            rejected += 1
            continue

        # Truncate the called TCID to the consensus depth, then look at the
        # MOST PER-HIT-CONSERVATIVE tier across the agreeing hits — i.e. the
        # worst tier of any hit. The level_kind is set from the consensus
        # depth (NOT the hit tier) since consensus may demote 5-part hits.
        n_parts = _AGREEMENT_TO_PARTS[consensus["agreement"]]
        called_tcid = truncate_tcid(consensus["tcid"], n_parts)
        worst_tier = max(h["tier"] for h in hits)
        # Effective tier: the deeper of (worst hit tier, depth-implied tier)
        depth_tier = {5: 1, 4: 2, 3: 3}[n_parts]
        effective_tier = max(worst_tier, depth_tier)

        # Best (highest-identity) hit drives the metadata fields
        best = max(hits, key=lambda h: h["identity"])

        egn_tcid = egn_lookup.get(query_id)
        agreement = compute_egn_agreement(called_tcid, egn_tcid)

        calls[query_id] = {
            "tcid": called_tcid,
            "level_kind": _TIER_TO_LEVEL_KIND[effective_tier],
            "tier": effective_tier,
            "identity": best["identity"],
            "qcov": best["qcov"],
            "scov": best["scov"],
            "evalue": best["evalue"],
            "length": best["length"],
            "consensus_n": consensus["n"],
            "consensus_agreement": consensus["agreement"],
            "egn_agreement": agreement,
            "egn_tcid": egn_tcid,
            "incompletely_characterized": is_class_9(called_tcid),
        }

    # Build summary
    tier_dist: dict[str, int] = defaultdict(int)
    agreement_dist: dict[str, int] = defaultdict(
        int, {"confirms": 0, "refines": 0, "extends": 0, "conflicts": 0}
    )
    for c in calls.values():
        tier_dist[str(c["tier"])] += 1
        agreement_dist[c["egn_agreement"]] += 1

    summary = {
        "raw_hit_lines": raw_lines,
        "proteins_with_hits": len(by_query),
        "proteins_with_call": len(calls),
        "proteins_rejected_by_consensus": rejected,
        "tier_distribution": dict(tier_dist),
        "agreement_distribution": dict(agreement_dist),
    }
    return calls, summary
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_tcdb_diamond.py -v`
Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/utils/tcdb_diamond.py tests/test_tcdb_diamond.py
git commit -m "feat(tcdb-diamond): add build_strain_calls pipeline composer"
```

---

## Task 11: TCDB DB build orchestrator

Wire up `extractTCDB.pl` invocation. Pure orchestration: given a target dir, ensure `tcdb.dmnd` is present, building it via Perl + diamond if not. Resolves the data dir from `TCDB_DATA_DIR` env var (default `~/tools/TCDB`).

**Files:**
- Modify: `.claude/skills/tcdb-diamond/run_tcdb_diamond.py`

- [ ] **Step 1: Add the imports + dotenv loader**

```python
# .claude/skills/tcdb-diamond/run_tcdb_diamond.py
"""Run diamond blastp vs. TCDB FASTA per strain. See spec
docs/superpowers/specs/2026-05-10-tcdb-diamond-augmentation-design.md.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import dotenv

REPO_ROOT = Path(__file__).resolve().parents[3]
GENOMES_CSV = REPO_ROOT / "data/Prochlorococcus/genomes/cyanobacteria_genomes.csv"
DEFAULT_TCDB_DATA_DIR = Path.home() / "tools" / "TCDB"

dotenv.load_dotenv(REPO_ROOT / ".env")
```

- [ ] **Step 2: Add the data-dir resolver**

```python
# Append to run_tcdb_diamond.py

def resolve_tcdb_data_dir() -> Path:
    """Resolve the TCDB data dir (env override or default)."""
    env = os.environ.get("TCDB_DATA_DIR")
    if env:
        return Path(env).expanduser()
    return DEFAULT_TCDB_DATA_DIR
```

- [ ] **Step 3: Add the DB-build function**

```python
# Append to run_tcdb_diamond.py

def build_tcdb_db(data_dir: Path, force: bool) -> Path:
    """Ensure ~/tools/TCDB/DB/tcdb.dmnd exists. Returns the path.

    Calls Saier Lab's extractTCDB.pl to download the FASTA and build the
    diamond DB. Skipped when the DB already exists and `force` is False.
    Requires:
      - extractTCDB.pl present at <data_dir>/TCDBtools/bin/
      - perl, wget, diamond on PATH
    """
    db_dir = data_dir / "DB"
    db_dir.mkdir(parents=True, exist_ok=True)
    dmnd = db_dir / "tcdb.dmnd"

    if dmnd.exists() and not force:
        return dmnd

    extract = data_dir / "TCDBtools" / "bin" / "extractTCDB.pl"
    if not extract.exists():
        print(
            f"ERROR: {extract} not found. Run:\n"
            f"  git clone https://github.com/SaierLaboratory/TCDBtools.git "
            f"{data_dir}/TCDBtools",
            file=sys.stderr,
        )
        sys.exit(1)

    for tool in ("perl", "wget", "diamond"):
        if not shutil.which(tool):
            print(f"ERROR: {tool} not on PATH (required by extractTCDB.pl)", file=sys.stderr)
            sys.exit(1)

    print(f"Building TCDB diamond DB in {db_dir} (this takes ~30s)...")
    for fmt in ("fasta", "diamond"):
        cmd = ["perl", str(extract), "-i", "tcdb", "-o", str(db_dir), "-f", fmt]
        result = subprocess.run(cmd, cwd=str(REPO_ROOT))
        if result.returncode != 0:
            print(f"ERROR: extractTCDB.pl -f {fmt} returned {result.returncode}", file=sys.stderr)
            sys.exit(1)

    if not dmnd.exists():
        print(f"ERROR: extractTCDB.pl completed but {dmnd} was not created", file=sys.stderr)
        sys.exit(1)

    return dmnd
```

- [ ] **Step 4: Verify the file parses**

Run: `uv run python -c "import importlib.util; spec=importlib.util.spec_from_file_location('rt','.claude/skills/tcdb-diamond/run_tcdb_diamond.py'); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); print(m.resolve_tcdb_data_dir())"`
Expected: prints `/home/<user>/tools/TCDB` (or whatever `TCDB_DATA_DIR` resolves to).

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/tcdb-diamond/run_tcdb_diamond.py
git commit -m "feat(tcdb-diamond): add TCDB DB build orchestrator"
```

---

## Task 12: Per-strain diamond runner (with log capture)

Wire up the `diamond blastp` invocation per spec §6.2. Diamond's stdout+stderr go to a per-strain log file under `logs/`, matching the existing convention (`*.log` is globally gitignored; `logs/prepare_data_step*.log`, `logs/eggnog_run_*.log` are the precedents).

**Files:**
- Modify: `.claude/skills/tcdb-diamond/run_tcdb_diamond.py`

- [ ] **Step 1: Add the LOGS_DIR constant near the other top-level constants**

```python
# Modify run_tcdb_diamond.py — add after GENOMES_CSV line
LOGS_DIR = REPO_ROOT / "logs"
```

- [ ] **Step 2: Add the run-diamond function (with log capture)**

```python
# Append to run_tcdb_diamond.py

def run_diamond(faa: Path, dmnd: Path, out_tsv: Path, threads: int, log_path: Path) -> bool:
    """Run diamond blastp for one strain. Returns True on success.

    Floor-only filtering at the diamond step: --evalue 0.001 only.
    Identity / coverage tiering happens in build_strain_calls (Python).

    Diamond's stdout + stderr are captured to `log_path` (overwritten on
    each run). The terminal only sees one progress line per strain.
    """
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "diamond", "blastp",
        "-q", str(faa),
        "-d", str(dmnd),
        "-o", str(out_tsv),
        "--outfmt", "6", "qseqid", "sseqid", "pident",
        "qcovhsp", "scovhsp", "length", "evalue", "bitscore",
        "--evalue", "0.001",
        "--max-target-seqs", "5",
        "--more-sensitive",
        "--threads", str(threads),
    ]
    print(f"\n>>> diamond blastp {faa.name} → {out_tsv} (log: {log_path.relative_to(REPO_ROOT)})")
    with open(log_path, "w") as logf:
        logf.write(f"$ {' '.join(cmd)}\n\n")
        logf.flush()
        result = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT)
    return result.returncode == 0


def truncate_faa(faa: Path, n_proteins: int, dest: Path) -> Path:
    """Copy the first N sequences of a FASTA to `dest`. Returns dest.

    Used by the --limit flag for fast end-to-end smoke tests against a small
    subset of one strain's proteome (~10-30s instead of ~5min).
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(faa) as src, open(dest, "w") as out:
        seen = 0
        for line in src:
            if line.startswith(">"):
                seen += 1
                if seen > n_proteins:
                    break
            out.write(line)
    return dest
```

- [ ] **Step 3: Verify the file still parses**

Run: `uv run python -c "import importlib.util; spec=importlib.util.spec_from_file_location('rt','.claude/skills/tcdb-diamond/run_tcdb_diamond.py'); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); print('ok')"`
Expected: prints `ok`.

- [ ] **Step 4: Commit**

```bash
git add .claude/skills/tcdb-diamond/run_tcdb_diamond.py
git commit -m "feat(tcdb-diamond): add per-strain diamond runner with log capture"
```

---

## Task 13: Genome-loop orchestrator + CLI

Tie everything together: load `cyanobacteria_genomes.csv`, iterate genome strains, run diamond + build_strain_calls, write outputs, print summary.

**Files:**
- Modify: `.claude/skills/tcdb-diamond/run_tcdb_diamond.py`

- [ ] **Step 1: Add the genome registry loader (mirror eggnog-run pattern)**

```python
# Append to run_tcdb_diamond.py

import csv

def load_genomes(strain_filter: str | None) -> list[dict]:
    """Parse cyanobacteria_genomes.csv, return rows for genome strains only.

    Skips reference_proteome_match and treatment organisms — only
    organism_type='genome_strain' rows have a `protein.faa`.
    """
    genomes: list[dict] = []
    with open(GENOMES_CSV) as f:
        reader = csv.DictReader(row for row in f if not row.strip().startswith("#"))
        for row in reader:
            if strain_filter and row["strain_name"] != strain_filter:
                continue
            if (row.get("organism_type") or "").strip() != "genome_strain":
                continue
            genomes.append(row)
    return genomes
```

- [ ] **Step 2: Add the per-strain processor**

```python
# Append to run_tcdb_diamond.py

from multiomics_kg.utils.tcdb_diamond import build_strain_calls


def process_strain(
    strain: str,
    data_dir_genome: Path,
    dmnd: Path,
    threads: int,
    force: bool,
    limit: int | None = None,
) -> tuple[str, str, dict | None]:
    """Run the per-strain pipeline. Returns (strain, status, summary_or_None).

    When `limit` is set, runs diamond against the first `limit` proteins from
    protein.faa (truncated copy in /tmp), and writes outputs to
    `<strain>.tcdb.limited_<N>.{tsv,calls.json,skill_summary.json}` to avoid
    overwriting full-run artifacts. Skip-if-exists logic also keys on the
    limit-suffixed name so re-running with a different limit is a fresh run.

    Status values:
      OK / SKIP_NO_FAA / SKIP_EXISTS / FAILED_DIAMOND / FAILED_NO_EGGNOG
    """
    faa = data_dir_genome / "protein.faa"
    if not faa.exists():
        return strain, "SKIP_NO_FAA", None

    out_dir = data_dir_genome / "tcdb"
    suffix = f".limited_{limit}" if limit else ""
    out_tsv = out_dir / f"{strain}.tcdb{suffix}.tsv"
    out_calls = out_dir / f"{strain}.tcdb{suffix}.calls.json"
    out_summary = out_dir / f"{strain}.tcdb{suffix}.skill_summary.json"

    if out_calls.exists() and not force:
        with open(out_summary) as f:
            return strain, "SKIP_EXISTS", json.load(f)

    if limit:
        truncated = Path("/tmp") / f"tcdb_diamond_{strain}_first{limit}.faa"
        faa = truncate_faa(faa, limit, truncated)

    log_path = LOGS_DIR / f"tcdb_diamond_{strain}{suffix}.log"
    if not run_diamond(faa, dmnd, out_tsv, threads, log_path):
        print(f"  see log: {log_path.relative_to(REPO_ROOT)}", file=sys.stderr)
        return strain, "FAILED_DIAMOND", None

    eggnog_path = data_dir_genome / "eggnog" / f"{strain}.emapper.annotations"
    if not eggnog_path.exists():
        # eggNOG missing → still proceed but agreement column will all be "extends"
        print(f"  WARN: {eggnog_path} missing — egn_agreement will all be 'extends'")

    calls, summary = build_strain_calls(out_tsv, eggnog_path)

    with open(out_calls, "w") as f:
        json.dump(calls, f, indent=2, sort_keys=True)
    with open(out_summary, "w") as f:
        json.dump(summary, f, indent=2)
    return strain, "OK", summary
```

- [ ] **Step 3: Add the main() entry point**

```python
# Append to run_tcdb_diamond.py

def main():
    parser = argparse.ArgumentParser(
        description="Run diamond vs. TCDB FASTA per strain (Phase 1)."
    )
    parser.add_argument("--strain", help="Run only this strain (e.g. MED4)")
    parser.add_argument("--force", action="store_true",
                        help="Re-run even if calls.json exists")
    parser.add_argument("--refresh-tcdb", action="store_true",
                        help="Re-download TCDB FASTA + diamond DB even if cached")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 4,
                        help="Diamond threads (default: os.cpu_count() or 4)")
    parser.add_argument("--limit", type=int, default=None,
                        help="Smoke test: run on first N proteins of each strain "
                             "only. Outputs go to <strain>.tcdb.limited_<N>.* "
                             "alongside (not replacing) full-run artifacts.")
    args = parser.parse_args()

    data_dir = resolve_tcdb_data_dir()
    print(f"TCDB data dir: {data_dir}")
    dmnd = build_tcdb_db(data_dir, force=args.refresh_tcdb)
    print(f"diamond DB: {dmnd}")

    genomes = load_genomes(args.strain)
    if not genomes:
        print(f"No genome_strain genomes found"
              f"{f' for strain {args.strain}' if args.strain else ''}.")
        sys.exit(1)

    results: list[tuple[str, str, dict | None]] = []
    for g in genomes:
        strain = g["strain_name"]
        data_dir_genome = REPO_ROOT / g["data_dir"].rstrip("/")
        results.append(
            process_strain(strain, data_dir_genome, dmnd, args.threads, args.force, args.limit)
        )

    # Status table
    print(f"\n{'='*92}")
    cols = ("Strain", "Status", "Hits", "T1", "T2", "T3",
            "confirms", "refines", "extends", "conflicts")
    print("{:<12} {:<14} {:>5} {:>4} {:>4} {:>4} {:>9} {:>8} {:>8} {:>10}".format(*cols))
    print("-" * 92)
    for strain, status, summary in results:
        if summary is None:
            print(f"{strain:<12} {status:<14}")
            continue
        td = summary["tier_distribution"]
        ad = summary["agreement_distribution"]
        print(
            "{:<12} {:<14} {:>5} {:>4} {:>4} {:>4} {:>9} {:>8} {:>8} {:>10}".format(
                strain, status, summary["proteins_with_call"],
                td.get("1", 0), td.get("2", 0), td.get("3", 0),
                ad.get("confirms", 0), ad.get("refines", 0),
                ad.get("extends", 0), ad.get("conflicts", 0),
            )
        )
    print("=" * 92)

    failed = [s for s, st, _ in results if st.startswith("FAILED")]
    if failed:
        print(f"\nFAILED strains: {', '.join(failed)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Verify the file still parses + main is callable**

Run: `uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --help`
Expected: argparse usage text including `--strain`, `--force`, `--refresh-tcdb`, `--threads`.

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/tcdb-diamond/run_tcdb_diamond.py
git commit -m "feat(tcdb-diamond): add genome-loop orchestrator + CLI"
```

---

## Task 14: SKILL.md write-up

Document the skill — frontmatter + workflow doc — mirroring the structure of `.claude/skills/eggnog-run/SKILL.md`.

**Files:**
- Modify: `.claude/skills/tcdb-diamond/SKILL.md`

- [ ] **Step 1: Replace the placeholder SKILL.md**

```markdown
---
name: tcdb-diamond
description: Run diamond blastp vs. the curated TCDB FASTA per strain to generate per-protein TCDB classifications, with tiered confidence (5-part / 4-part / 3-part) plus consensus + eggNOG-agreement tags. Phase 1 — produces inspectable `<strain>.tcdb.calls.json` artifacts; KG integration deferred to Phase 2.
argument-hint: [--strain <name> | --force | --refresh-tcdb | --threads <n>]
user-invocable: true
allowed-tools: Read, Bash(uv *), Bash(perl *), Bash(diamond *), Bash(git *)
---

# TCDB-Diamond Skill

Per-strain `diamond blastp` against the curated TCDB FASTA. Augments eggNOG's
family-level (3-part) `KEGG_TC` annotations with sequence-similarity-based calls
that can reach `tc_specificity` (5-part) when identity warrants it.

## Quick Start

```bash
# Run all genome strains (skip already-done)
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py

# Run a single strain
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --strain MED4

# Force re-run even if calls.json exists
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --force

# Re-download TCDB FASTA + diamond DB
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --refresh-tcdb

# Use more threads
uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --threads 8
```

## One-Time Setup

```bash
# Saier Lab's TCDB tools (used to download the FASTA + build the diamond DB)
git clone https://github.com/SaierLaboratory/TCDBtools.git ~/tools/TCDB/TCDBtools
```

Optional `.env` entry to relocate the TCDB data dir (default: `~/tools/TCDB`):

```
TCDB_DATA_DIR=/path/to/TCDB
```

System tools required: `perl`, `wget`, `diamond` (already required by other parts of the project).

## What It Does

1. Reads `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv`, filters to `organism_type='genome_strain'`.
2. Ensures `~/tools/TCDB/DB/tcdb.dmnd` exists, building it via `extractTCDB.pl -i tcdb -f diamond` if missing or `--refresh-tcdb`.
3. For each strain, runs `diamond blastp` of `<data_dir>/protein.faa` against the TCDB diamond DB, writing raw 8-column TSV to `<data_dir>/tcdb/<strain>.tcdb.tsv`.
4. Applies the per-hit tier policy + per-protein consensus collapse (Spec §6.3, §6.4-A).
5. Joins each protein's eggNOG `KEGG_TC` (from `<data_dir>/eggnog/<strain>.emapper.annotations`) and computes the agreement tag (Spec §6.4-B).
6. Tags class-9 calls (Spec §6.4-C).
7. Writes `<data_dir>/tcdb/<strain>.tcdb.calls.json` (per-protein records) and `<strain>.tcdb.skill_summary.json` (per-strain stats).
8. Prints a status table to stdout. Diamond's full per-strain stdout+stderr are captured to `logs/tcdb_diamond_<strain>.log` (auto-gitignored via `*.log`).

## Tier Policy

| Tier | Truncate to | Identity | Coverage rule | Notes |
|---|---|---|---|---|
| 1 | 5 parts (`tc_specificity`) | ≥ 70% | qcov ≥ 70% | High-confidence specificity call |
| 2 | 4 parts (`tc_subfamily`) | ≥ 40% | qcov ≥ 60% | Solid subfamily call |
| 3 | 3 parts (`tc_family`) | (no floor) | qcov ≥ 40% **OR** scov ≥ 40% | gblast3-style floor |

All tiers also require: e-value ≤ 0.001, HSP length ≥ 50 aa.

## Output Schema (`<strain>.tcdb.calls.json`)

Keyed by NCBI protein_id (WP_ accession):

```json
{
  "WP_011131900.1": {
    "tcid": "1.A.11.1.5",
    "level_kind": "tc_specificity",
    "tier": 1,
    "identity": 87.4,
    "qcov": 92.1,
    "scov": 89.7,
    "evalue": 1.2e-180,
    "length": 412,
    "consensus_n": 5,
    "consensus_agreement": "5_part",
    "egn_agreement": "refines",
    "egn_tcid": "1.A.11",
    "incompletely_characterized": false
  }
}
```

`egn_agreement` values: `confirms` (same TC) | `refines` (diamond deeper than eggNOG, same lineage) | `extends` (eggNOG had no TC) | `conflicts` (different family).

## Phase 2 (Future)

Phase 1 artifacts sit in the strain cache for inspection. Phase 2 (separate spec) will integrate them into `gene_annotations_merged.json` and the KG. See [docs/superpowers/specs/2026-05-10-tcdb-diamond-augmentation-design.md §7](../../../docs/superpowers/specs/2026-05-10-tcdb-diamond-augmentation-design.md).

## Workflow When Invoked

1. Verify one-time setup is done (`~/tools/TCDB/TCDBtools/bin/extractTCDB.pl` exists).
2. Run the skill: `uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py`.
3. Review the status table — note any FAILED strains.
4. Inspect `<data_dir>/tcdb/<strain>.tcdb.calls.json` for spot checks.
```

- [ ] **Step 2: Verify the markdown is well-formed**

Run: `head -10 .claude/skills/tcdb-diamond/SKILL.md`
Expected: shows the YAML frontmatter starting with `---\nname: tcdb-diamond`.

- [ ] **Step 3: Commit**

```bash
git add .claude/skills/tcdb-diamond/SKILL.md
git commit -m "docs(tcdb-diamond): add SKILL.md with workflow and tier policy"
```

---

## Task 15: Smoke run on ONE Alteromonas strain (MIT1002)

End-to-end verification on a **single strain**. Alteromonas is chosen because it stresses more code paths than a cyano: ~2× the proteome size (~4K vs ~2K proteins), ~3× the eggNOG-TC density (~437 vs ~120 TC genes), and very low UniProt curation coverage (~28% vs ~66%) — surfaces issues that MED4 wouldn't.

**Do not run all strains in this task.** That happens later, after the user has reviewed the MIT1002 output and confirmed the skill behaves as expected. This task is ONLY a verification of the skill on one strain.

**Files:**
- No code changes; verification only. May need to update SKILL.md or fix bugs found.

- [ ] **Step 1: Verify TCDBtools is cloned and TCDB_DATA_DIR is set**

Run: `ls ~/tools/TCDB/TCDBtools/bin/extractTCDB.pl && grep TCDB_DATA_DIR .env`
Expected: file exists; .env entry shows `TCDB_DATA_DIR=...`.

- [ ] **Step 2: Run the full unit-test suite first**

Run: `pytest tests/test_tcdb_diamond.py -v`
Expected: all tests pass (from Tasks 2-10).

- [ ] **Step 3a: Truncated smoke run (~10-30s, ~100 proteins)**

Run: `uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --strain MIT1002 --limit 100`
Expected:
- Prints "TCDB data dir: …" and "diamond DB: …".
- First run also triggers `extractTCDB.pl` to download the FASTA + build the diamond DB (~30s; reports "Building TCDB diamond DB in …"). Subsequent runs reuse it.
- One `>>> diamond blastp` line referencing `/tmp/tcdb_diamond_MIT1002_first100.faa`.
- Status table with `MIT1002 OK` plus per-tier and per-agreement counts. Counts will be small (~5-20 hits across 100 proteins) but non-zero.
- Outputs land at `cache/data/Alteromonas/genomes/MIT1002/tcdb/MIT1002.tcdb.limited_100.{tsv,calls.json,skill_summary.json}`.

- [ ] **Step 3b: Inspect the truncated calls.json**

Run:
```bash
ls cache/data/Alteromonas/genomes/MIT1002/tcdb/
python3 -c "
import json
with open('cache/data/Alteromonas/genomes/MIT1002/tcdb/MIT1002.tcdb.limited_100.calls.json') as f:
    calls = json.load(f)
print(f'protein count: {len(calls)}')
print(f'first 5 entries:')
for k in list(calls)[:5]:
    print(f'  {k}: {calls[k]}')
print(f'tier counts: ', {t: sum(1 for c in calls.values() if c[\"tier\"]==t) for t in (1,2,3)})
print(f'agreement counts: ', {a: sum(1 for c in calls.values() if c[\"egn_agreement\"]==a) for a in (\"confirms\",\"refines\",\"extends\",\"conflicts\")})
"
```
Expected:
- protein count is small but non-zero (typically 5-30 hits in the first 100 proteins of an Alteromonas FASTA).
- Each entry has all 13 keys from §6.5: `tcid`, `level_kind`, `tier`, `identity`, `qcov`, `scov`, `evalue`, `length`, `consensus_n`, `consensus_agreement`, `egn_agreement`, `egn_tcid`, `incompletely_characterized`.
- `tier` values are 1, 2, or 3.

- [ ] **Step 3c: Inspect the truncated diamond log**

Run: `tail -30 logs/tcdb_diamond_MIT1002.limited_100.log`
Expected: shows the diamond invocation command + diamond's own progress output (number of queries processed, timing).

- [ ] **Step 4: Pause for user review of truncated output**

**STOP HERE.** Show the inspection output (step 3b) to the user. They confirm the format and merge logic look correct. **Only after their OK** proceed to step 5 — the full-strain run.

If anything was wrong (parsing error, unexpected JSON shape, …) fix and re-run from step 3a.

- [ ] **Step 5: Full-strain run on MIT1002 (~2-5 min)**

Run: `uv run python .claude/skills/tcdb-diamond/run_tcdb_diamond.py --strain MIT1002`
Expected:
- One `>>> diamond blastp` line writing to `…/MIT1002.tcdb.tsv` (no `.limited_<N>` suffix).
- Status table with `MIT1002 OK` plus per-tier and per-agreement counts.
- Final outputs at `cache/data/Alteromonas/genomes/MIT1002/tcdb/MIT1002.tcdb.{tsv,calls.json,skill_summary.json}`.

- [ ] **Step 6: Inspect the full output**

Run:
```bash
python3 -c "
import json
with open('cache/data/Alteromonas/genomes/MIT1002/tcdb/MIT1002.tcdb.calls.json') as f:
    calls = json.load(f)
print(f'protein count: {len(calls)}')
print(f'tier counts: ', {t: sum(1 for c in calls.values() if c[\"tier\"]==t) for t in (1,2,3)})
print(f'agreement counts: ', {a: sum(1 for c in calls.values() if c[\"egn_agreement\"]==a) for a in (\"confirms\",\"refines\",\"extends\",\"conflicts\")})
print(f'incompletely_characterized count: ', sum(1 for c in calls.values() if c[\"incompletely_characterized\"]))
"
```
Expected:
- protein count is in the ~400-700 range (MIT1002 has ~427 eggNOG-TC genes; diamond may add `extends` for non-eggNOG-TC genes that hit TCDB directly).
- `egn_agreement` distribution is sensible: substantive `confirms`, some `refines` (the headline specificity win), notable `extends`, small minority `conflicts`.

- [ ] **Step 7: Confirm no regression in existing tests**

Run: `pytest -m "not slow and not kg" -q`
Expected: all pass — Phase 1 only adds new files, so existing tests should be unaffected.

- [ ] **Step 8: Pause for user review of full-strain output, hand back**

**STOP HERE.** Show the step 6 output to the user. They decide whether the full distributions look reasonable before approving an all-25-strain run. **The full run is a follow-up activity, not part of this plan.**

```bash
git status   # confirm only new files in cache/data/Alteromonas/genomes/MIT1002/tcdb/ + (if any) test/skill spot-fixes
git diff --stat
```

If anything was changed in response to issues found during inspection, commit it:

```bash
git add -A
git commit -m "verify(tcdb-diamond): smoke-test MIT1002 end-to-end"
```

Otherwise (no fixes needed), no commit is required — the cache outputs aren't source code. Hand back to the user for review.

---

## Self-Review

**Spec coverage check** — every section of [docs/superpowers/specs/2026-05-10-tcdb-diamond-augmentation-design.md](../specs/2026-05-10-tcdb-diamond-augmentation-design.md) maps to a task:

| Spec section | Task(s) |
|---|---|
| §3 Architecture (single component, no KG coupling) | All tasks (no `build_gene_annotations.py` changes) |
| §4 File layout (~/tools/TCDB/DB/, cache/.../tcdb/) | Task 11 (DB), Task 13 (per-strain output) |
| §5 One-time setup | Task 14 (documented in SKILL.md) |
| §6.1 TCDB DB build | Task 11 |
| §6.2 Per-strain diamond invocation | Task 12 |
| §6.3 Per-hit tier policy | Task 3 |
| §6.4-A Top-N consensus collapse | Task 4 |
| §6.4-B eggNOG agreement tag | Task 5 |
| §6.4-C Class-9 tag | Task 6 |
| §6.5 calls.json schema | Tasks 10, 13 (composer + write-out) |
| §6.6 Status table | Task 13 |
| §10 Acceptance criteria 1-5 | Task 15 (smoke run) |
| §10 Acceptance criterion 6 (no regression) | Task 15 step 5 |
| §10 Acceptance criterion 7 (existing tests pass) | Task 15 step 5 |

**Type/name consistency check** — verified across tasks:
- `parse_diamond_row` returns dict with keys `query_id`, `subject_id`, `identity`, `qcov`, `scov`, `length`, `evalue`, `bitscore` (Task 7) — used unchanged in Task 10
- `parse_tcdb_subject_id` returns `tuple[str, str] | None` (Task 8) — Task 10 destructures as `_, hit_tcid`
- `consensus_collapse` returns dict with keys `tcid`, `agreement`, `n` (Task 4) — Task 10 reads all three
- `compute_egn_agreement` returns one of `confirms | refines | extends | conflicts` (Task 5) — Task 10 stores it in `egn_agreement`
- `build_strain_calls` returns `(calls, summary)` (Task 10) — Task 13 unpacks both
- `process_strain` returns `(strain, status, summary | None)` (Task 13) — `main()` iterates this shape

**Placeholder scan** — no TBD/TODO/"add appropriate" — all code shown inline.

---

## Execution Handoff

**Plan complete and saved to [docs/superpowers/plans/2026-05-10-tcdb-diamond-skill.md](docs/superpowers/plans/2026-05-10-tcdb-diamond-skill.md). Two execution options:**

1. **Subagent-Driven (recommended)** — I dispatch a fresh subagent per task, review between tasks, fast iteration.
2. **Inline Execution** — Execute tasks in this session using executing-plans, batch execution with checkpoints.

**Which approach?**
