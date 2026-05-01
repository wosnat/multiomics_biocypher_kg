# Phase 1.2.1 — Unified pruned `kegg_data.json` cache

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor `kegg_data.json` from "full KEGG REST data, pruned at iteration time by each adapter" to "single gene-reachable pruned cache, both adapters read from it." Adds a new `Metabolite_in_pathway` edge type using compound→pathway data we already parse but don't currently use.

**Architecture:** Step 6 of `prepare_data.sh` becomes the single owner of the pruned cache. It walks every strain's `gene_annotations_merged.json` to identify gene-reachable {KOs, reactions, compounds, pathways}, prunes the raw KEGG data to that subset (with KO ∪ Rxn-reachable pathway semantic per Option B), enriches reactions/compounds with MNX cross-refs, and writes a single `cache/data/kegg/kegg_data.json` (~500 KB). The 5.4 MB full version currently committed in git is replaced. The separate `cache/data/kegg/kegg_metabolism_xrefs.json` is deleted (its content folds in). Both `kegg_anno_adapter` and `metabolism_adapter` become thin readers — no per-adapter pruning code.

**Tech Stack:** Python 3.10+, `sqlite3` (stdlib), `pytest` with `tmp_path` fixtures, the existing raw-cache infrastructure from commit `adfb153`.

**Spec basis:** Design conversation 2026-04-30 (memory: "questions - follow with me"). Locked decisions:
- Object-per-entity nested structure (kos / pathways / subcategories / categories / reactions / compounds)
- Strict mode: adapters error if `kegg_data.json` missing
- BRITE stays separate (12 raw `brite_ko*.json` files unchanged)
- `compound_to_pathways` flat dict dropped; folded into `compounds[c].pathways`
- New `Metabolite_in_pathway` edge with **Option B** semantic: only emit edges to pathways that are also KO- or Rxn-reachable (drops 83 compound-only pathways like Leishmaniasis, MAPK signaling, etc.)
- Replace existing committed `kegg_data.json` (5.4 MB → ~500 KB)

**Depends on:** Phase 1.2 complete (commits `be3f25a` → `91047f9`).

---

## Pruned `kegg_data.json` structure (target)

```json
{
  "kos": {
    "K02338": {
      "name": "DNA polymerase III subunit beta",
      "pathways": ["ko03030", "ko03430", "ko03440"]
    }
  },
  "pathways": {
    "ko03030": {
      "name": "DNA replication",
      "subcategory": "09124"
    }
  },
  "subcategories": {
    "09124": {"name": "Replication and repair", "category": "09120"}
  },
  "categories": {
    "09120": {"name": "Genetic Information Processing"}
  },
  "reactions": {
    "R00200": {
      "name": "ATP:pyruvate 2-O-phosphotransferase ...",
      "ec_numbers": ["2.7.1.40"],
      "pathways": ["ko00010", "ko00710"],
      "compounds": ["C00074", "C00008", "C00022", "C00002"],
      "mnxr_id": "MNXR101234",
      "rhea_ids": ["10828"],
      "mass_balance": "balanced",
      "reaction_class": "chemical"
    }
  },
  "compounds": {
    "C00031": {
      "name": "D-glucose",
      "formula": "C6H12O6",
      "mass": 180.063,
      "inchikey": "WQZGKKKJIJFFOK-GASJEMHNSA-N",
      "smiles": "OC[CH]1...",
      "mnxm_id": "MNXM41",
      "chebi_id": "17234",
      "hmdb_id": "HMDB0000122",
      "pathways": ["ko00010", "ko00500"]
    }
  }
}
```

Pathway-set rule: `pathways = (KOs reachable from genes) ∪ (Reactions reachable from genes) → all pathways linked from those`.

For **`compounds[c].pathways`**: filtered to that same set (Option B — drops compound-only pathways).

For **`reactions[r].pathways`**: same set (the `91047f9` filter from Phase 1.2 stays semantically — but is now intrinsic to the pruned cache rather than enforced separately).

---

## File structure

**Modify:**

| Path | Why |
|---|---|
| `multiomics_kg/utils/kegg_utils.py` | Remove `download_kegg_data()` (its role splits to step 6); keep raw fetchers + parsers; add `load_kegg_data(cache_root)` thin loader |
| `multiomics_kg/download/build_kegg_metabolism_xrefs.py` | Rewrite as unified pruned-cache builder; writes `kegg_data.json` directly. Optionally rename to `build_kegg_data.py` (deferred — see Task 2 note). |
| `multiomics_kg/adapters/kegg_annotation_adapter.py` | Iterate new structure (`kegg_data["kos"]`, `kegg_data["pathways"]`); drop iteration-time pruning |
| `multiomics_kg/adapters/metabolism_adapter.py` | Read from `kegg_data.json` instead of `kegg_metabolism_xrefs.json`; add `_metabolite_pathway_edges()` |
| `config/schema_config.yaml` | Add `metabolite to kegg pathway association` edge |
| `scripts/prepare_data.sh` | Step 6 description / log message updates |
| `tests/test_kegg_utils_metabolism.py` | Drop tests for removed `download_kegg_data()`; keep parser tests |
| `tests/test_build_kegg_metabolism_xrefs.py` | Rewrite for new pruned-cache shape |
| `tests/test_metabolism_adapter.py` | Update fixtures to new shape |
| `tests/test_kegg_annotation_adapter.py` | Update for new iteration shape |
| `tests/kg_validity/test_metabolism.py` | Add `Metabolite_in_pathway` count check |
| `cache/data/kegg/kegg_data.json` | Replace committed 5.4 MB full version with pruned ~500 KB version |
| `CLAUDE.md` | Document unified cache architecture + new edge type |

**Delete:**

| Path | Why |
|---|---|
| `cache/data/kegg/kegg_metabolism_xrefs.json` | Folded into `kegg_data.json` |

**Create:** none (everything is a refactor of existing files).

---

## Tasks

### Task 1: Schema — `metabolite_in_pathway` edge

**Files:**
- Modify: `config/schema_config.yaml`

Additive change. No code consumers yet (the adapter emission lands in Task 5).

- [ ] **Step 1: Add the edge entry**

After the existing `reaction to kegg pathway association:` block in `config/schema_config.yaml`, add:

```yaml
metabolite to kegg pathway association:
  is_a: association
  represented_as: edge
  label_as_edge: metabolite_in_pathway
  source: metabolite
  target: kegg term
  label_in_input: metabolite_in_pathway
```

- [ ] **Step 2: Verify YAML parses**

```bash
cd /home/osnat/github/multiomics_biocypher_kg && uv run python -c "import yaml; yaml.safe_load(open('config/schema_config.yaml'))"
```
Expected: no output, exit 0.

- [ ] **Step 3: Verify BioCypher accepts**

```bash
cd /home/osnat/github/multiomics_biocypher_kg && uv run python -c "from biocypher import BioCypher; bc = BioCypher(schema_config_path='config/schema_config.yaml', biocypher_config_path='config/biocypher_config.yaml'); print('OK')" 2>&1 | tail -3
```
Expected: prints `OK` with no traceback.

- [ ] **Step 4: Commit**

```bash
git add config/schema_config.yaml
git commit -m "metabolism: 1.2.1 — add metabolite_in_pathway edge to schema"
```

---

### Task 2: Step 6 rewrite — unified pruned `kegg_data.json` builder

**Files:**
- Modify: `multiomics_kg/download/build_kegg_metabolism_xrefs.py` (full rewrite)
- Modify: `tests/test_build_kegg_metabolism_xrefs.py` (full rewrite)
- Possibly modify: `multiomics_kg/utils/kegg_utils.py` (remove `download_kegg_data`, add `download_kegg_raw` if absent — see step 3 below)

This is the largest task. The module's role changes from "build a metabolism-only xrefs cache that consumes `kegg_data.json` written by `download_kegg_data()`" to "ensure raw KEGG cache is populated, then build a unified pruned `kegg_data.json` with new structure (folding in metabolism enrichment)."

**File rename note**: keeping the name `build_kegg_metabolism_xrefs.py` for now. Renaming to `build_kegg_data.py` would be cleaner but adds churn (path changes in `scripts/prepare_data.sh`, log file `prepare_data_step6.log`, tests). Defer the rename to a follow-up.

- [ ] **Step 1: Verify or add `download_kegg_raw()` in `kegg_utils.py`**

Read the current `kegg_utils.py`. Today `download_kegg_data()` does both download-and-parse-and-write. We need to ensure a function exists that just downloads raw to `cache/data/kegg/raw/<endpoint>.txt` (or `.json` for BRITE). The `_fetch_text_cached` and `_fetch_json_cached` helpers added in commit `adfb153` already do per-endpoint download-and-cache.

Add a new wrapper if not already present:

```python
def download_kegg_raw(cache_root: Path, force: bool = False) -> None:
    """Ensure all 9 KEGG raw endpoints are cached locally.

    Hits the network only for missing files (or when force=True).
    Output: cache/data/kegg/raw/<endpoint>.{txt,json}
    """
    raw_dir = Path(cache_root) / "kegg" / "raw"
    for url in [
        _KO_LIST_URL, _KO_PATHWAY_LINK_URL, _PATHWAY_KO_LIST_URL,
        _REACTION_LIST_URL, _COMPOUND_LIST_URL,
        _LINK_COMPOUND_REACTION_URL, _LINK_PATHWAY_REACTION_URL,
        _LINK_PATHWAY_COMPOUND_URL,
    ]:
        _fetch_text_cached(url, raw_dir=raw_dir, force=force)
    _fetch_json_cached(_BRITE_KO_URL, raw_dir=raw_dir, force=force)
```

(If `download_kegg_raw` already exists from `adfb153`, skip this step.)

Remove the old `download_kegg_data()` function — its role moves to step 6.

- [ ] **Step 2: Write the failing tests for the new build flow**

Replace the tests in `tests/test_build_kegg_metabolism_xrefs.py` with tests that exercise the new pruned-cache shape. Target structure:

```python
"""Unit tests for step 6: unified pruned kegg_data.json builder."""
from __future__ import annotations
import json
import sqlite3
from pathlib import Path

import pytest

from multiomics_kg.download import build_kegg_metabolism_xrefs as bx


def _write_strain_annotations(strain_dir, genes):
    strain_dir.mkdir(parents=True, exist_ok=True)
    (strain_dir / "gene_annotations_merged.json").write_text(json.dumps(genes))


def _make_resolver(tmp_path):
    """Minimal MNX resolver with synthetic compounds and reactions."""
    db = tmp_path / "resolver.db"
    conn = sqlite3.connect(db)
    conn.executescript("""
        CREATE TABLE compounds (mnxm_id TEXT PRIMARY KEY, name TEXT, formula TEXT,
            mass REAL, inchi TEXT, inchikey TEXT, smiles TEXT, charge INTEGER, reference TEXT);
        CREATE TABLE compound_aliases (source TEXT, value TEXT, mnxm_id TEXT,
            PRIMARY KEY (source, value, mnxm_id));
        CREATE INDEX idx_compound_aliases_mnxm ON compound_aliases(mnxm_id);
        CREATE TABLE reactions (mnxr_id TEXT PRIMARY KEY, mnx_equation TEXT,
            reference TEXT, classifs TEXT, is_balanced TEXT, is_transport TEXT);
        CREATE TABLE reaction_aliases (source TEXT, value TEXT, mnxr_id TEXT,
            PRIMARY KEY (source, value, mnxr_id));
        CREATE INDEX idx_reaction_aliases_mnxr ON reaction_aliases(mnxr_id);

        INSERT INTO compounds VALUES
            ('MNXM41', 'D-glucose', 'C6H12O6', 180.063, 'I=', 'WQZ-...', 'OC[CH]1...', 0, 'chebi:17234');
        INSERT INTO compound_aliases VALUES ('kegg.compound', 'C00031', 'MNXM41');
        INSERT INTO compound_aliases VALUES ('chebi', '17234', 'MNXM41');
        INSERT INTO reactions VALUES
            ('MNXR101234', '...', 'kegg.reaction:R00200', '2.7.1.40', 'B', NULL);
        INSERT INTO reaction_aliases VALUES ('kegg.reaction', 'R00200', 'MNXR101234');
    """)
    conn.commit()
    return conn


@pytest.fixture
def synthetic_kegg_raw(tmp_path):
    """Return a kegg_data dict matching the parser output, suitable as input to step 6."""
    return {
        "ko_names":               {"K02338": "DNA polymerase III subunit beta",
                                   "K01006": "pyruvate kinase"},
        "pathway_names":          {"ko00010": "Glycolysis", "ko00710": "Carbon fixation",
                                   "ko03030": "DNA replication",
                                   "ko00500": "Starch and sucrose metabolism"},
        "ko_to_pathways":         {"K02338": ["ko03030"], "K01006": ["ko00010"]},
        "pathway_to_subcategory": {"ko00010": "09101", "ko00710": "09102",
                                   "ko03030": "09124", "ko00500": "09101"},
        "subcategory_names":      {"09101": "Carbohydrate metabolism",
                                   "09102": "Energy metabolism",
                                   "09124": "Replication and repair"},
        "subcategory_to_category":{"09101": "09100", "09102": "09100", "09124": "09120"},
        "category_names":         {"09100": "Metabolism", "09120": "Genetic Information Processing"},
        "reaction_names":         {"R00200": "ATP:pyruvate ..."},
        "compound_names":         {"C00031": "D-glucose", "C00074": "PEP",
                                   "C99999": "compound-only-pathway-trigger"},
        "reaction_to_compounds":  {"R00200": ["C00031", "C00074"]},
        "reaction_to_pathways":   {"R00200": ["ko00010"]},
        "compound_to_pathways":   {"C00031": ["ko00010", "ko00500"],
                                   "C99999": ["ko99999"]},   # ko99999 is compound-only
    }


def test_build_pruned_kegg_data(tmp_path, monkeypatch, synthetic_kegg_raw):
    """End-to-end: 1 strain, 2 genes, ensure pruned structure is correct."""
    strain = tmp_path / "strain"
    _write_strain_annotations(strain, {
        "gene1": {"kegg_ko": ["K02338"], "kegg_reactions": []},
        "gene2": {"kegg_ko": ["K01006"], "kegg_reactions": ["R00200"]},
    })
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [{"data_dir": str(strain)}])

    conn = _make_resolver(tmp_path)
    out_path = tmp_path / "kegg_data.json"
    bx.build_pruned_kegg_data(synthetic_kegg_raw, conn, out_path)

    data = json.loads(out_path.read_text())

    # KOs: only gene-annotated KOs survive
    assert set(data["kos"].keys()) == {"K02338", "K01006"}
    assert data["kos"]["K02338"]["name"] == "DNA polymerase III subunit beta"
    assert data["kos"]["K02338"]["pathways"] == ["ko03030"]

    # Pathways: KO ∪ Rxn-reachable. K02338→ko03030, K01006→ko00010, R00200→ko00010
    # ko00710 not present — no KO/Rxn maps there.
    assert set(data["pathways"].keys()) == {"ko03030", "ko00010"}

    # Reactions: only gene-reachable
    assert set(data["reactions"].keys()) == {"R00200"}
    rxn = data["reactions"]["R00200"]
    assert rxn["name"].startswith("ATP:pyruvate")
    assert rxn["pathways"] == ["ko00010"]
    assert rxn["compounds"] == ["C00031", "C00074"]
    assert rxn["mnxr_id"] == "MNXR101234"
    assert rxn["ec_numbers"] == ["2.7.1.40"]
    assert rxn["mass_balance"] == "balanced"
    assert rxn["reaction_class"] == "chemical"

    # Compounds: only reaction-reachable
    assert set(data["compounds"].keys()) == {"C00031", "C00074"}
    glucose = data["compounds"]["C00031"]
    assert glucose["name"] == "D-glucose"
    assert glucose["formula"] == "C6H12O6"
    assert glucose["mnxm_id"] == "MNXM41"
    assert glucose["chebi_id"] == "17234"

    # Option B: compound.pathways pruned to KO ∪ Rxn-reachable set.
    # ko00500 is NOT in the set (no KO or Rxn maps there in the fixture).
    # ko00010 IS in the set.
    assert glucose["pathways"] == ["ko00010"]

    # Subcategories/categories present iff parent of an emitted pathway
    assert set(data["subcategories"].keys()) == {"09101", "09124"}
    assert data["subcategories"]["09101"]["category"] == "09100"
    assert set(data["categories"].keys()) == {"09100", "09120"}


def test_compound_only_pathways_dropped(tmp_path, monkeypatch, synthetic_kegg_raw):
    """The C99999→ko99999 compound-only pathway should not appear anywhere.

    Critical for Option B: even if a compound were in the pruned set, its pathways
    should be filtered to KO ∪ Rxn-reachable. C99999 isn't even in the pruned set
    (no reaction reaches it), so this is a double-defensive check.
    """
    strain = tmp_path / "strain"
    _write_strain_annotations(strain, {
        "gene1": {"kegg_ko": ["K02338"], "kegg_reactions": []},
    })
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [{"data_dir": str(strain)}])
    conn = _make_resolver(tmp_path)
    out_path = tmp_path / "kegg_data.json"
    bx.build_pruned_kegg_data(synthetic_kegg_raw, conn, out_path)
    data = json.loads(out_path.read_text())
    assert "C99999" not in data["compounds"]
    assert "ko99999" not in data["pathways"]


def test_compound_pathway_filter_drops_unevidenced(tmp_path, monkeypatch):
    """Compound pathways are filtered to KO ∪ Rxn-reachable.

    Force a case where a compound IS in the gene-reachable set (via R00200) but has
    a pathway (ko00500) that's not reachable from any gene-KO or gene-reaction.
    Expected: compound emitted, but ko00500 NOT in compound.pathways.
    """
    raw = {
        "ko_names":               {"K01006": "pyruvate kinase"},
        "pathway_names":          {"ko00010": "Glycolysis", "ko00500": "Starch and sucrose"},
        "ko_to_pathways":         {"K01006": ["ko00010"]},
        "pathway_to_subcategory": {"ko00010": "09101", "ko00500": "09101"},
        "subcategory_names":      {"09101": "Carbohydrate metabolism"},
        "subcategory_to_category":{"09101": "09100"},
        "category_names":         {"09100": "Metabolism"},
        "reaction_names":         {"R00200": "rxn"},
        "compound_names":         {"C00031": "D-glucose"},
        "reaction_to_compounds":  {"R00200": ["C00031"]},
        "reaction_to_pathways":   {"R00200": ["ko00010"]},
        "compound_to_pathways":   {"C00031": ["ko00010", "ko00500"]},  # ko00500 unevidenced
    }
    strain = tmp_path / "strain"
    _write_strain_annotations(strain, {
        "g1": {"kegg_ko": ["K01006"], "kegg_reactions": ["R00200"]},
    })
    monkeypatch.setattr(bx, "load_genome_rows", lambda: [{"data_dir": str(strain)}])
    conn = _make_resolver(tmp_path)
    out_path = tmp_path / "kegg_data.json"
    bx.build_pruned_kegg_data(raw, conn, out_path)
    data = json.loads(out_path.read_text())

    # ko00500 must NOT be in pathways top-level (Option B: only KO ∪ Rxn-reachable)
    assert "ko00500" not in data["pathways"]
    # And compound's pathway list filtered
    assert data["compounds"]["C00031"]["pathways"] == ["ko00010"]
```

Run: `cd /home/osnat/github/multiomics_biocypher_kg && uv run pytest tests/test_build_kegg_metabolism_xrefs.py -v`
Expected: 3 FAIL — `AttributeError: module 'build_kegg_metabolism_xrefs' has no attribute 'build_pruned_kegg_data'`.

- [ ] **Step 3: Rewrite the module**

Replace the existing `multiomics_kg/download/build_kegg_metabolism_xrefs.py` with:

```python
"""Step 6 — unified pruned KEGG data cache.

Walks every strain's gene_annotations_merged.json to identify gene-reachable
KEGG entities, prunes the raw KEGG dataset to that subset, enriches reactions
and compounds with MNX cross-refs, and writes a single
cache/data/kegg/kegg_data.json (~500 KB).

Both kegg_anno_adapter and metabolism_adapter read this file. No per-adapter
pruning at iteration time.

Pathway-set rule (Option B): pathway IDs in the cache are exactly those reachable
from gene-KOs or gene-reactions. compound→pathway lists are filtered to that set,
so compound-only pathways (e.g. ko05140 Leishmaniasis) are excluded.

Usage:
    uv run python -m multiomics_kg.download.build_kegg_metabolism_xrefs [--force]
"""
from __future__ import annotations

import argparse
import json
import logging
import sqlite3
from pathlib import Path

from multiomics_kg.download.utils.cli import load_genome_rows
from multiomics_kg.utils.gene_id_utils import load_gene_annotations
from multiomics_kg.utils import kegg_utils
from multiomics_kg.utils import metabolite_utils as mu

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent

KEGG_CACHE_DIR = PROJECT_ROOT / "cache" / "data" / "kegg"
KEGG_DATA_FILE = KEGG_CACHE_DIR / "kegg_data.json"
RESOLVER_DB = PROJECT_ROOT / "cache" / "data" / "mnx" / "metabolite_resolver.db"

log = logging.getLogger(__name__)


# ── Reachability ──────────────────────────────────────────────────────────────

def _gene_reachable_sets(raw: dict) -> dict[str, set[str]]:
    """Walk all strains to compute the gene-reachable {KOs, reactions, pathways, compounds}.

    Returns a dict with keys 'kos', 'rxns', 'cpds', 'pws'.
    Pathway set is KO ∪ Rxn-reachable (Option B).
    """
    kos: set[str] = set()
    rxns: set[str] = set()
    for row in load_genome_rows():
        genes = load_gene_annotations(row["data_dir"])
        if not genes:
            continue
        for g in genes.values():
            for ko in g.get("kegg_ko") or []:
                if isinstance(ko, str) and ko.startswith("K"):
                    kos.add(ko)
            for rxn in g.get("kegg_reactions") or []:
                if isinstance(rxn, str) and rxn.startswith("R"):
                    rxns.add(rxn)

    rxn_to_cpds = raw.get("reaction_to_compounds", {})
    cpds: set[str] = set()
    for r in rxns:
        cpds.update(rxn_to_cpds.get(r, []))

    ko_to_pw = raw.get("ko_to_pathways", {})
    rxn_to_pw = raw.get("reaction_to_pathways", {})
    pws: set[str] = set()
    for k in kos:
        pws.update(ko_to_pw.get(k, []))
    for r in rxns:
        pws.update(rxn_to_pw.get(r, []))

    log.info(
        f"Gene-reachable: {len(kos)} KOs, {len(rxns)} reactions, "
        f"{len(cpds)} compounds, {len(pws)} pathways (KO∪Rxn)"
    )
    return {"kos": kos, "rxns": rxns, "cpds": cpds, "pws": pws}


# ── Hierarchy parents ─────────────────────────────────────────────────────────

def _hierarchy_parents(raw: dict, pws: set[str]) -> tuple[set[str], set[str]]:
    """Compute the subcategory + category sets that need to be emitted."""
    pw_to_sub = raw.get("pathway_to_subcategory", {})
    sub_to_cat = raw.get("subcategory_to_category", {})
    subs = {pw_to_sub[p] for p in pws if p in pw_to_sub}
    cats = {sub_to_cat[s] for s in subs if s in sub_to_cat}
    return subs, cats


# ── MNX enrichment ────────────────────────────────────────────────────────────

def _enrich_reaction(rxn_id: str, raw: dict, conn: sqlite3.Connection,
                     allowed_pathways: set[str]) -> dict:
    out: dict = {
        "name": raw.get("reaction_names", {}).get(rxn_id, ""),
        "pathways": [p for p in raw.get("reaction_to_pathways", {}).get(rxn_id, []) if p in allowed_pathways],
        "compounds": list(raw.get("reaction_to_compounds", {}).get(rxn_id, [])),
        "ec_numbers": [],
        "mnxr_id": None,
        "rhea_ids": [],
        "mass_balance": "unbalanced",
        "reaction_class": "chemical",
    }
    cur = conn.cursor()
    cur.execute("SELECT mnxr_id FROM reaction_aliases WHERE source='kegg.reaction' AND value=? LIMIT 1", (rxn_id,))
    row = cur.fetchone()
    if row is None:
        return out
    out["mnxr_id"] = row[0]

    cur.execute("SELECT value FROM reaction_aliases WHERE mnxr_id=? AND source='rhea' ORDER BY value", (out["mnxr_id"],))
    out["rhea_ids"] = [r[0] for r in cur.fetchall()]

    cur.execute("SELECT classifs, is_balanced, is_transport FROM reactions WHERE mnxr_id=?", (out["mnxr_id"],))
    row = cur.fetchone()
    if row:
        classifs, is_balanced, is_transport = row
        if classifs:
            out["ec_numbers"] = [c.strip() for c in classifs.split(";") if c.strip()]
        out["mass_balance"] = "balanced" if (is_balanced or "").upper() == "B" else "unbalanced"
        out["reaction_class"] = "transport" if (is_transport or "").upper() == "T" else "chemical"
    return out


def _enrich_compound(cpd_id: str, raw: dict, conn: sqlite3.Connection,
                     allowed_pathways: set[str]) -> dict:
    out: dict = {
        "name": raw.get("compound_names", {}).get(cpd_id, ""),
        "formula": None, "mass": None, "inchikey": None, "smiles": None,
        "mnxm_id": None, "chebi_id": None, "hmdb_id": None,
        "pathways": [p for p in raw.get("compound_to_pathways", {}).get(cpd_id, []) if p in allowed_pathways],
    }
    cur = conn.cursor()
    cur.execute("SELECT mnxm_id FROM compound_aliases WHERE source='kegg.compound' AND value=? LIMIT 1", (cpd_id,))
    row = cur.fetchone()
    if row is None:
        return out
    out["mnxm_id"] = row[0]

    cur.execute("SELECT formula, mass, inchikey, smiles FROM compounds WHERE mnxm_id=?", (out["mnxm_id"],))
    row = cur.fetchone()
    if row:
        out["formula"] = row[0] or None
        out["mass"] = row[1] if row[1] is not None else None
        out["inchikey"] = row[2] or None
        out["smiles"] = row[3] or None

    cur.execute("SELECT value FROM compound_aliases WHERE mnxm_id=? AND source='chebi' ORDER BY value LIMIT 1", (out["mnxm_id"],))
    row = cur.fetchone()
    if row:
        out["chebi_id"] = row[0]
    cur.execute("SELECT value FROM compound_aliases WHERE mnxm_id=? AND source='hmdb' ORDER BY value LIMIT 1", (out["mnxm_id"],))
    row = cur.fetchone()
    if row:
        out["hmdb_id"] = row[0]
    return out


# ── Build ─────────────────────────────────────────────────────────────────────

def build_pruned_kegg_data(raw: dict, conn: sqlite3.Connection, out_path: Path) -> None:
    """Build the pruned kegg_data.json from raw KEGG data + MNX resolver."""
    sets = _gene_reachable_sets(raw)
    kos, rxns, cpds, pws = sets["kos"], sets["rxns"], sets["cpds"], sets["pws"]
    subs, cats = _hierarchy_parents(raw, pws)

    ko_to_pw = raw.get("ko_to_pathways", {})
    pw_to_sub = raw.get("pathway_to_subcategory", {})
    sub_to_cat = raw.get("subcategory_to_category", {})

    out = {
        "kos": {
            k: {
                "name": raw.get("ko_names", {}).get(k, ""),
                "pathways": [p for p in ko_to_pw.get(k, []) if p in pws],
            } for k in sorted(kos)
        },
        "pathways": {
            p: {
                "name": raw.get("pathway_names", {}).get(p, ""),
                **({"subcategory": pw_to_sub[p]} if p in pw_to_sub else {}),
            } for p in sorted(pws)
        },
        "subcategories": {
            s: {
                "name": raw.get("subcategory_names", {}).get(s, ""),
                **({"category": sub_to_cat[s]} if s in sub_to_cat else {}),
            } for s in sorted(subs)
        },
        "categories": {
            c: {"name": raw.get("category_names", {}).get(c, "")}
            for c in sorted(cats)
        },
        "reactions": {
            r: _enrich_reaction(r, raw, conn, pws) for r in sorted(rxns)
        },
        "compounds": {
            c: _enrich_compound(c, raw, conn, pws) for c in sorted(cpds)
        },
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, separators=(",", ":")))
    log.info(
        f"Wrote {out_path}: {len(out['kos'])} KOs, {len(out['pathways'])} pathways, "
        f"{len(out['reactions'])} reactions, {len(out['compounds'])} compounds"
    )


def _parse_raw_into_dict(cache_root: Path) -> dict:
    """Parse raw KEGG cache files into the same dict shape that download_kegg_data
    used to return (so build_pruned_kegg_data can stay agnostic of the on-disk layout).
    """
    raw_dir = cache_root / "kegg" / "raw"
    return {
        "ko_names":               kegg_utils._parse_ko_names(    (raw_dir / "list_ko.txt").read_text()),
        "pathway_names":          kegg_utils._parse_pathway_ko_names((raw_dir / "list_pathway_ko.txt").read_text()),
        "ko_to_pathways":         kegg_utils._parse_ko_to_pathways((raw_dir / "link_pathway_ko.txt").read_text()),
        "reaction_names":         kegg_utils._parse_reaction_names((raw_dir / "list_reaction.txt").read_text()),
        "compound_names":         kegg_utils._parse_compound_names((raw_dir / "list_compound.txt").read_text()),
        "reaction_to_compounds":  kegg_utils._parse_reaction_to_compounds((raw_dir / "link_compound_reaction.txt").read_text()),
        "reaction_to_pathways":   kegg_utils._parse_reaction_to_pathways((raw_dir / "link_pathway_reaction.txt").read_text()),
        "compound_to_pathways":   kegg_utils._parse_compound_to_pathways((raw_dir / "link_pathway_compound.txt").read_text()),
        # BRITE hierarchy supplements pathway names + builds the 4-level hierarchy
        **_parse_brite_supplements(raw_dir),
    }


def _parse_brite_supplements(raw_dir: Path) -> dict:
    """Parse br_ko00001.json into pathway_to_subcategory / subcategory_names /
    subcategory_to_category / category_names plus extra pathway names.
    """
    brite_json = json.loads((raw_dir / "br_ko00001.json").read_text())
    pw_to_sub, sub_names, sub_to_cat, cat_names, brite_pw_names = (
        kegg_utils._parse_brite_hierarchy(brite_json)
    )
    return {
        "pathway_to_subcategory": pw_to_sub,
        "subcategory_names":      sub_names,
        "subcategory_to_category": sub_to_cat,
        "category_names":         cat_names,
        # Note: pathway_names comes from list_pathway_ko.txt which is also parsed; BRITE
        # supplements those (but list_pathway_ko has more, so list wins on conflicts.)
    }


def main(force: bool = False) -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    cache_root = PROJECT_ROOT / "cache" / "data"

    if KEGG_DATA_FILE.exists() and not force:
        log.info(f"{KEGG_DATA_FILE} exists; use --force to rebuild.")
        return

    if not RESOLVER_DB.exists():
        raise FileNotFoundError(
            f"{RESOLVER_DB} missing — run prepare_data.sh step 0 sub-step 7 first."
        )

    log.info("Ensuring KEGG raw cache (downloads from KEGG REST if missing) ...")
    kegg_utils.download_kegg_raw(cache_root, force=force)

    log.info("Parsing raw KEGG into in-memory dict ...")
    raw = _parse_raw_into_dict(cache_root)

    log.info("Opening MNX resolver ...")
    conn = mu.open_resolver(RESOLVER_DB)

    log.info("Building pruned kegg_data.json ...")
    build_pruned_kegg_data(raw, conn, KEGG_DATA_FILE)

    conn.close()
    log.info("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Rebuild even if kegg_data.json exists.")
    args = parser.parse_args()
    main(force=args.force)
```

- [ ] **Step 4: Run unit tests**

```bash
cd /home/osnat/github/multiomics_biocypher_kg && uv run pytest tests/test_build_kegg_metabolism_xrefs.py -v
```
Expected: 3 PASS.

- [ ] **Step 5: Update `kegg_utils.py`**

Open `multiomics_kg/utils/kegg_utils.py`. Remove `download_kegg_data()` and `_RAW_FILES`/`_fetch_text_cached`/`_fetch_json_cached` infrastructure if they're now needed only by `download_kegg_raw()`. Add `download_kegg_raw()` (see Step 1) if absent.

Update tests in `tests/test_kegg_utils_metabolism.py`:
- Remove `test_download_kegg_data_includes_metabolism_keys` (covers removed function).
- Remove `test_download_kegg_data_uses_raw_cache_on_second_call` (same).
- Keep all `_parse_*` tests.
- Add a new test for `download_kegg_raw` that monkeypatches `_fetch_text` / `_fetch_json` and verifies the 9 raw files are written.

```python
def test_download_kegg_raw_populates_raw_dir(tmp_path, monkeypatch):
    """download_kegg_raw should populate cache/data/kegg/raw/ with 8 .txt + 1 .json."""
    text_calls = []
    json_calls = []
    monkeypatch.setattr(kegg_utils, "_fetch_text",
                        lambda url: (text_calls.append(url), "stub").__getitem__(1))
    monkeypatch.setattr(kegg_utils, "_fetch_json",
                        lambda url: (json_calls.append(url), {"children": []}).__getitem__(1))
    kegg_utils.download_kegg_raw(tmp_path, force=False)
    raw_dir = tmp_path / "kegg" / "raw"
    txt_files = sorted(p.name for p in raw_dir.glob("*.txt"))
    json_files = sorted(p.name for p in raw_dir.glob("*.json"))
    assert len(txt_files) == 8
    assert len(json_files) == 1
    assert "br_ko00001.json" in json_files
    # Second call without force: no new fetches
    text_calls.clear(); json_calls.clear()
    kegg_utils.download_kegg_raw(tmp_path, force=False)
    assert text_calls == [] and json_calls == []
    # Force: re-fetches
    kegg_utils.download_kegg_raw(tmp_path, force=True)
    assert len(text_calls) == 8 and len(json_calls) == 1
```

- [ ] **Step 6: Run kegg_utils tests**

```bash
cd /home/osnat/github/multiomics_biocypher_kg && uv run pytest tests/test_kegg_utils_metabolism.py tests/test_kegg_utils.py -v
```
Expected: all PASS.

- [ ] **Step 7: Commit**

```bash
git add multiomics_kg/download/build_kegg_metabolism_xrefs.py \
        multiomics_kg/utils/kegg_utils.py \
        tests/test_build_kegg_metabolism_xrefs.py \
        tests/test_kegg_utils_metabolism.py
git commit -m "metabolism: 1.2.1 — step 6 builds unified pruned kegg_data.json"
```

---

### Task 3: `kegg_annotation_adapter` refactor — read new structure

**Files:**
- Modify: `multiomics_kg/adapters/kegg_annotation_adapter.py`
- Modify: `tests/test_kegg_annotation_adapter.py`

The adapter has 4 emitter loops today (KO nodes, pathway nodes, subcategory nodes, category nodes + 3 hierarchy edge loops + gene→KO edge loop).

Today: each loop iterates a flat dict (e.g. `kegg_data["ko_names"]`) and applies per-strain pruning by checking gene KOs.

Tomorrow: the cache is already pruned. Each loop becomes a straight iteration of `kegg_data["kos"]` etc., with `name` / `pathways` / `subcategory` / `category` etc. read from the nested entry.

- [ ] **Step 1: Read the current adapter**

Read `multiomics_kg/adapters/kegg_annotation_adapter.py` carefully. Note its constructor signature (likely `MultiKeggAnnotationAdapter(genome_config_file, cache_root, test_mode, cache)`), its `download_data()` method (calls `download_kegg_data()`), and the iteration loops in `get_nodes()` / `get_edges()`.

Then read `tests/test_kegg_annotation_adapter.py` to see what fixtures look like.

- [ ] **Step 2: Rewrite for new structure**

Key changes in the adapter:
1. `download_data()` now ensures `kegg_data.json` exists; if not, raises a clear error pointing to `prepare_data.sh --steps 6`.
2. `_load()` returns the new pruned dict directly (using `kegg_utils.load_kegg_data(cache_root)`).
3. `get_nodes()` / `get_edges()` iterate the new structure:
   - `for ko_id, ko in data["kos"].items(): yield (..., ko["name"], ...)`
   - `for pw_id, pw in data["pathways"].items(): yield (..., pw["name"], ...)`
   - `for ko_id, ko in data["kos"].items(): for pw in ko["pathways"]: yield Kegg_term_is_a_kegg_term edge`
   - etc.
4. Remove per-strain pruning code — the cache is pre-pruned. Gene→KO edges still need per-strain iteration over `gene["kegg_ko"]`, but only emit if the KO is in `data["kos"]`.

Add `kegg_utils.load_kegg_data(cache_root) -> dict`:
```python
def load_kegg_data(cache_root: Path) -> dict:
    """Load the pruned kegg_data.json. Raises FileNotFoundError if missing."""
    p = Path(cache_root) / "kegg" / "kegg_data.json"
    if not p.exists():
        raise FileNotFoundError(
            f"{p} missing — run `bash scripts/prepare_data.sh --steps 6` first."
        )
    return json.loads(p.read_text())
```

- [ ] **Step 3: Update tests**

The test fixtures will need to switch from the flat dict shape to the nested object-per-entity shape. Most tests should still work with semantic-equivalent fixtures.

- [ ] **Step 4: Run tests**

```bash
cd /home/osnat/github/multiomics_biocypher_kg && uv run pytest tests/test_kegg_annotation_adapter.py -v
```
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/kegg_annotation_adapter.py \
        multiomics_kg/utils/kegg_utils.py \
        tests/test_kegg_annotation_adapter.py
git commit -m "metabolism: 1.2.1 — kegg_annotation_adapter reads new pruned structure"
```

---

### Task 4: `metabolism_adapter` refactor — read kegg_data.json + emit Metabolite_in_pathway

**Files:**
- Modify: `multiomics_kg/adapters/metabolism_adapter.py`
- Modify: `tests/test_metabolism_adapter.py`

- [ ] **Step 1: Update the adapter to read from `kegg_data.json`**

Today: `MetabolismAdapter.__init__` takes `xrefs_path`, defaults to `cache/data/kegg/kegg_metabolism_xrefs.json`. `_load()` loads that file.

Tomorrow: `MetabolismAdapter.__init__` takes `kegg_data_path`, defaults to `cache/data/kegg/kegg_data.json`. `_load()` loads that file.

Update the constructor signature, default path, and `_load()` body. The data structure changes from `{"reactions": {...}, "compounds": {...}}` to `{"kos": {...}, "pathways": {...}, ..., "reactions": {...}, "compounds": {...}}` — adapter only consumes `reactions` and `compounds` keys (kegg_anno_adapter handles the rest), so the iteration code stays largely the same.

Field naming changes inside reactions/compounds:
- Reactions: today `compound_ids` and `kegg_pathway_ids` → tomorrow `compounds` and `pathways`.
- Compounds: today no `pathways` field → tomorrow has `pathways` field for the new edge.

Update `_reaction_metabolite_edges` and `_reaction_pathway_edges` to read the new field names.

- [ ] **Step 2: Add `_metabolite_pathway_edges()` and emit it from `get_edges()`**

```python
def _metabolite_pathway_edges(self) -> Iterator[tuple[str, str, str, str, dict]]:
    data = self._load()
    for cpd_id, cpd in data.get("compounds", {}).items():
        for pw_id in cpd.get("pathways", []):
            edge_id = f"m2p:{cpd_id}:{pw_id}"
            yield (
                edge_id,
                f"kegg.compound:{cpd_id}",
                f"kegg.pathway:{pw_id}",
                "metabolite_in_pathway",
                {},
            )
```

In `MetabolismAdapter.get_edges()`:
```python
def get_edges(self):
    yield from self._reaction_metabolite_edges()
    yield from self._reaction_pathway_edges()
    yield from self._metabolite_pathway_edges()  # NEW
```

In `MultiMetabolismAdapter.get_edges()` add it after the other edge yields:
```python
def get_edges(self):
    yield from self._gene_reaction_edges()
    yield from self._reaction_metabolite_edges()
    yield from self._reaction_pathway_edges()
    yield from self._metabolite_pathway_edges()  # NEW
```

- [ ] **Step 3: Update tests**

Update fixtures in `tests/test_metabolism_adapter.py` to use the new structure (object-per-entity). Add a new test:

```python
def test_metabolite_pathway_edges(xrefs_file):
    """Each compound emits one edge per pathway in compound.pathways."""
    adapter = ma.MetabolismAdapter(kegg_data_path=xrefs_file)
    edges = list(adapter._metabolite_pathway_edges())
    pairs = {(e[1], e[2]) for e in edges}
    # Use the fixture's compound + pathways list
    assert ("kegg.compound:C00074", "kegg.pathway:ko00010") in pairs
    assert {e[3] for e in edges} == {"metabolite_in_pathway"}
```

Also update `test_get_edges_yields_all_three_types` → 4 types now:

```python
def test_get_edges_yields_all_four_types(...):
    ...
    labels = {e[3] for e in adapter.get_edges()}
    assert labels == {
        "gene_catalyzes_reaction",
        "reaction_has_metabolite",
        "reaction_in_kegg_pathway",
        "metabolite_in_pathway",
    }
```

- [ ] **Step 4: Run tests**

```bash
cd /home/osnat/github/multiomics_biocypher_kg && uv run pytest tests/test_metabolism_adapter.py -v
```
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/adapters/metabolism_adapter.py tests/test_metabolism_adapter.py
git commit -m "metabolism: 1.2.1 — metabolism_adapter reads kegg_data.json + emits metabolite_in_pathway"
```

---

### Task 5: Replace committed `kegg_data.json`, delete obsolete `kegg_metabolism_xrefs.json`

**Files:**
- Replace (in git): `cache/data/kegg/kegg_data.json` (5.4 MB full → ~500 KB pruned)
- Delete (from git): `cache/data/kegg/kegg_metabolism_xrefs.json`

- [ ] **Step 1: Run step 6 to regenerate the pruned cache**

```bash
cd /home/osnat/github/multiomics_biocypher_kg
rm -f cache/data/kegg/kegg_data.json cache/data/kegg/kegg_metabolism_xrefs.json
bash scripts/prepare_data.sh --steps 6 --force
ls -lh cache/data/kegg/kegg_data.json
```
Expected: produces a kegg_data.json ~500 KB (NOT 5.4 MB — verify size dropped).

- [ ] **Step 2: Verify the new structure with a quick grep**

```bash
uv run python -c "
import json
data = json.loads(open('cache/data/kegg/kegg_data.json').read())
print('Top-level keys:', sorted(data.keys()))
print('Sample KO:', list(data['kos'].items())[0])
print('Sample reaction:', list(data['reactions'].items())[0])
print('Sample compound:', list(data['compounds'].items())[0])
print('Counts:')
for k in ('kos','pathways','subcategories','categories','reactions','compounds'):
    print(f'  {k}: {len(data.get(k, {}))}')
"
```
Expected counts (from the empirical analysis we did): ~4534 KOs, ~377 pathways, ~25 subcategories, ~7 categories, ~2340 reactions, ~2184 compounds. Sample compound should have a `pathways` field.

- [ ] **Step 3: Run the full KG build (test mode for speed)**

```bash
uv run python create_knowledge_graph.py --test 2>&1 | tail -10
```
Expected: build completes without errors. Spot-check biocypher-out for Reaction-part000.csv, Metabolite-part000.csv, Metabolite_in_pathway-part000.csv.

- [ ] **Step 4: Stage the changes**

```bash
git add cache/data/kegg/kegg_data.json
git rm cache/data/kegg/kegg_metabolism_xrefs.json   # may not be in git yet — handle either way
git commit -m "metabolism: 1.2.1 — replace full kegg_data.json with pruned cache; drop kegg_metabolism_xrefs.json"
```

If `kegg_metabolism_xrefs.json` was never committed (was a gitignored cache), just `rm` it from disk:
```bash
rm -f cache/data/kegg/kegg_metabolism_xrefs.json
```

---

### Task 6: KG validity tests + step 6 wiring + CLAUDE.md

**Files:**
- Modify: `tests/kg_validity/test_metabolism.py`
- Modify: `scripts/prepare_data.sh`
- Modify: `CLAUDE.md`

- [ ] **Step 1: Add Metabolite_in_pathway count test**

In `tests/kg_validity/test_metabolism.py`, add:

```python
def test_metabolite_in_pathway_edges_present(neo4j_session):
    """Spec 1.2.1: ~5K-15K Metabolite_in_pathway edges (Option B-pruned)."""
    n = neo4j_session.run(
        "MATCH ()-[r:Metabolite_in_pathway]->() RETURN count(r) AS n"
    ).single()["n"]
    assert 5000 <= n <= 15000, f"{n} Metabolite_in_pathway edges (expected 5000-15000)"


def test_metabolite_pathways_only_kg_evidenced(neo4j_session):
    """Option B: Metabolite_in_pathway targets must also be reachable from genes."""
    # No metabolite-in-pathway edge should target a pathway with no Gene→KO edge to it
    # (since Option B prunes compound-only pathways at step 6).
    n_orphan = neo4j_session.run("""
        MATCH (m:Metabolite)-[:Metabolite_in_pathway]->(p:KeggTerm)
        WHERE NOT EXISTS {
            MATCH (g:Gene)-[:Gene_has_kegg_ko]->(:KeggTerm)-[:Kegg_term_is_a_kegg_term]->(p)
        }
        AND NOT EXISTS {
            MATCH (g:Gene)-[:Gene_catalyzes_reaction]->(:Reaction)-[:Reaction_in_kegg_pathway]->(p)
        }
        RETURN count(DISTINCT p) AS n
    """).single()["n"]
    assert n_orphan == 0, f"{n_orphan} pathways have only Metabolite-in-pathway edges (Option B violation)"
```

- [ ] **Step 2: Update step 6 description in `prepare_data.sh`**

In the header comment block, change "Build KEGG metabolism xrefs" → "Build pruned KEGG data cache (kegg_data.json) with metabolism enrichment" or similar. Update inline help text. Update the `run_step` description string.

- [ ] **Step 3: Update CLAUDE.md**

Edit `CLAUDE.md`:
- **Step 6** description: rewrite to reflect unified pruned cache role
- **Key Adapters** description for `metabolism_adapter.py` + brief note on `kegg_annotation_adapter` reading the pruned cache
- **Data Locations**: change `cache/data/kegg/kegg_metabolism_xrefs.json` reference → describe `kegg_data.json` as the pruned single source
- **Key graph facts**: add `Metabolite_in_pathway` edge with brief description
- **Actual Neo4j labels**: add `Metabolite_in_pathway` to relationships list

- [ ] **Step 4: Final test run**

```bash
cd /home/osnat/github/multiomics_biocypher_kg && uv run pytest -m "not slow and not kg" -q 2>&1 | tail -5
```
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add tests/kg_validity/test_metabolism.py scripts/prepare_data.sh CLAUDE.md
git commit -m "metabolism: 1.2.1 — KG validity tests + step 6 wiring + CLAUDE.md docs"
```

---

## Validation gate (manual, after all 6 tasks land)

```bash
# 1. Re-run step 6 (already done in Task 5, but force-rebuild to confirm idempotence)
bash scripts/prepare_data.sh --steps 6 --force

# 2. Full KG rebuild via Docker
docker compose down
docker compose up -d build
docker compose up -d import post-process deploy
docker compose logs -f post-process | grep -E "Done|ERROR"

# 3. Verify import.report has zero dangling edges
cat output/import.report   # should be empty or just headers

# 4. Run all KG validity tests
uv run pytest -m kg -v

# 5. Spot-check Metabolite_in_pathway edge counts via cypher-shell
docker exec deploy cypher-shell -a bolt://localhost:7687 \
  "MATCH ()-[r:Metabolite_in_pathway]->() RETURN count(r);"
```

Acceptance:
- `kegg_data.json` ≤ 600 KB
- Reaction count: ~2,340 (unchanged from Phase 1.2)
- Metabolite count: ~2,184 (unchanged)
- Pathway count (KeggTerm with `level_kind='pathway'`): ~377 (was 372 — adds 5 reaction-only)
- `Metabolite_in_pathway` edges: ~5K-15K
- import.report: zero dangling edges
- All KG validity tests pass

---

## Self-review checklist

1. **Spec coverage:** Each Phase 1.2.1 design decision in the header has a task that implements it. ✓
2. **No placeholders:** Every step has actual code or commands. ✓
3. **Type / name consistency:** `kegg_data.json` schema across Tasks 2/3/4 matches; field renames (`compound_ids` → `compounds`, `kegg_pathway_ids` → `pathways`) consistent. ✓
4. **Order discipline:** Schema first (additive), then step 6 owner change, then both adapters (in either order), then file replacement, then tests/docs. ✓
5. **Rollback safety:** Until Task 5 lands, the existing committed `kegg_data.json` (5.4 MB full) is still on disk and adapters work against the old structure. The new structure only goes live when Task 5 commits the replacement. Tasks 3 and 4 update adapters to read NEW structure — but won't run against actual data until Task 5. So the order **must be**: 1, 2, 3, 4, 5, 6 (no parallelism between Tasks 3-4 and Task 5).
