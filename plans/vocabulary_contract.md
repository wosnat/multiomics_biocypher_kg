# Controlled-Vocabulary Contract Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Publish the KG's controlled vocabularies as queryable `ControlledVocabulary` nodes, and align the InterPro/TCDB vocabularies to five house rules — so the explorer stops hard-coding value sets that drift silently into wrong answers.

**Architecture:** One YAML declares every vocabulary; a loader makes it importable by adapters (so an undeclared value cannot be emitted) and a node-only adapter emits it into the graph, following the existing `DataSourceAdapter` / `gene_annotations_config.yaml` pattern. Renames are applied at the emission site in each adapter and mirrored in the post-import Cypher. Four properties are deleted rather than renamed. Four test gates check observed-vs-declared at increasing cost.

**Tech Stack:** Python 3 · BioCypher · Neo4j 5.15 + cypher-shell · pytest · PyYAML · Docker Compose

**Spec:** `docs/superpowers/specs/2026-08-16-vocabulary-contract-design.md` (rev 5) — read §3 (house rules), §4 (rename table) and §5 (contract shape) before starting.

## Global Constraints

- **String sanitization is mandatory.** Every string property emitted from an adapter passes through a local `_clean_str(v) -> v.replace("'", "^").replace("|", "")`. `|` is the BioCypher array delimiter and `'` breaks `neo4j-admin import`.
- **`scripts/post-import.sh` and `scripts/post-import.cypher` must stay byte-identical in their Cypher logic.** Every change to one is applied to the other. The `.sh` is authoritative (Docker runs it); the `.cypher` is the non-Docker reference copy.
- **No native `bool` properties.** House rule R5. `value_type` does not admit `bool`; a two-state fact is a categorical string naming both states meaningfully.
- **Vocabulary values are lowercase `snake_case`.** House rule R1. Namespace-prefix a value only when the same property name carries values from different ontologies across labels (`level_kind` does; `interpro_type` does not).
- **Every `sources` value is the `id` of a `DataSource` node.** House rule R2. The eight ids are `ncbi, cyanorak, eggnog, uniprot, psortb, signalp, interproscan, tcdb_diamond`.
- **`evidence_score` is a float in `[0,1]`**, rounded to 3 decimals, on every annotation edge. House rule R4.
- **Do not change `annotation_quality`, `annotation_state`, `level_kind` values, or `DataSource` ids.** They are released and MCP-read. Out of scope.
- Run unit tests with `pytest -m "not slow and not kg"`. Graph tests need a running Docker Neo4j at `localhost:7687`.

---

## File Structure

| File | Responsibility | Task |
|---|---|---|
| `config/controlled_vocabularies.yaml` | **Create.** Single source of truth for every declared vocabulary | 1, 10 |
| `multiomics_kg/utils/controlled_vocab.py` | **Create.** Load the YAML; expose `VOCAB.check()` for adapters and `VOCAB.entries()` for the adapter/tests | 1 |
| `multiomics_kg/adapters/controlled_vocabulary_adapter.py` | **Create.** Node-only adapter emitting one node per vocabulary | 2 |
| `config/schema_config.yaml` | **Modify.** Add `controlled vocabulary` node type; apply renames; delete four properties | 2, 3, 5, 6, 7 |
| `create_knowledge_graph.py` | **Modify.** Register the adapter; write the vocabulary hash to the output dir | 2, 9 |
| `multiomics_kg/utils/interproscan.py` | **Modify.** Lowercase `interpro_type` + `libraries` at the parse site | 3 |
| `multiomics_kg/adapters/interpro_adapter.py` | **Modify.** Lowercase compare; delete Layer-A `ambiguous` / `source_db` | 3, 7 |
| `multiomics_kg/download/build_gene_annotations.py` | **Modify.** `interpro` → `interproscan` source label | 4 |
| `multiomics_kg/adapters/tcdb_adapter.py` | **Modify.** `diamond` → `tcdb_diamond`; `substrate_depth` values | 4, 6 |
| `multiomics_kg/utils/annotation_provenance.py` | **Modify.** Normalized float `evidence_score` | 5 |
| `scripts/post-import.{sh,cypher}` | **Modify.** Score normalization; bool→string; delete two flags, inline one threshold; hash stamp | 5, 6, 7, 9 |
| `tests/test_controlled_vocab.py` | **Create.** Loader + `check()` unit tests | 1 |
| `tests/test_controlled_vocabulary_adapter.py` | **Create.** Adapter emission tests | 2 |
| `tests/test_create_knowledge_graph.py` | **Modify.** Parameterize over build mode; add CSV vocabulary scan | 8 |
| `tests/kg_validity/test_controlled_vocabularies.py` | **Create.** Live-graph both-direction check + `DataSource` join + hash | 8 |
| `CLAUDE.md`, `CHANGELOG.md`, `docs/kg-changes/*.md` | **Modify.** House rules, renames, five chemistry definitions | 11 |

---

## Task 1: Vocabulary config + loader

**Files:**
- Create: `config/controlled_vocabularies.yaml`
- Create: `multiomics_kg/utils/controlled_vocab.py`
- Test: `tests/test_controlled_vocab.py`

**Interfaces:**
- Consumes: nothing
- Produces:
  - `VocabEntry` dataclass with fields `applies_to: str`, `applies_to_kind: str`, `property: str`, `value_type: str`, `closed: bool`, `values: list[str]`, `sparse: bool`, `expected_empty: bool`, `description: str`, `min_value: float | None`, `max_value: float | None`, `signal_count: int | None`, `signals: list[str]`, and a computed `id` property returning `f"{applies_to}.{property}"`
  - `load_vocabularies(path: str | Path = "config/controlled_vocabularies.yaml") -> dict[str, VocabEntry]` keyed by `id`
  - `VOCAB` module-level lazy singleton with `VOCAB.check(applies_to: str, property: str, value: str) -> str` (returns `value`, raises `ValueError` on undeclared) and `VOCAB.entries() -> list[VocabEntry]`
  - `vocabularies_hash(entries: list[VocabEntry]) -> str` returning `"sha256:<hex>"` over a canonical JSON dump

- [ ] **Step 1: Write the failing test**

```python
# tests/test_controlled_vocab.py
import pytest
from multiomics_kg.utils.controlled_vocab import (
    VOCAB, load_vocabularies, vocabularies_hash,
)


def test_loads_the_shipped_config():
    entries = load_vocabularies()
    assert "Gene_catalyzes_ec_number.evidence" in entries
    e = entries["Gene_catalyzes_ec_number.evidence"]
    assert e.applies_to_kind == "edge"
    assert e.value_type == "string"
    assert set(e.values) == {"curated", "family_inferred"}


def test_evidence_domain_is_per_edge_type():
    """domain_inferred is not a possible value on the EC edge (spec §5.2)."""
    entries = load_vocabularies()
    assert "domain_inferred" not in entries["Gene_catalyzes_ec_number.evidence"].values
    assert "domain_inferred" in entries[
        "Gene_involved_in_biological_process.evidence"].values


def test_check_returns_declared_value():
    assert VOCAB.check("Gene_catalyzes_ec_number", "evidence", "curated") == "curated"


def test_check_raises_on_undeclared_value():
    with pytest.raises(ValueError, match="domain_inferred"):
        VOCAB.check("Gene_catalyzes_ec_number", "evidence", "domain_inferred")


def test_check_passes_through_open_vocabularies():
    """closed: false means enumerate at runtime — check must not reject."""
    assert VOCAB.check("Gene", "gene_category", "anything at all") == "anything at all"


def test_bool_value_type_is_rejected():
    """R5: native bool is not an admissible value_type."""
    import yaml, tempfile, pathlib
    bad = {"Bad.prop": {"applies_to": "Bad", "applies_to_kind": "node",
                        "property": "prop", "value_type": "bool",
                        "closed": True, "values": [], "description": "x"}}
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as f:
        yaml.safe_dump(bad, f)
        path = f.name
    with pytest.raises(ValueError, match="bool"):
        load_vocabularies(path)
    pathlib.Path(path).unlink()


def test_hash_is_stable_and_order_independent():
    entries = list(load_vocabularies().values())
    h1 = vocabularies_hash(entries)
    h2 = vocabularies_hash(list(reversed(entries)))
    assert h1 == h2
    assert h1.startswith("sha256:")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_controlled_vocab.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'multiomics_kg.utils.controlled_vocab'`

- [ ] **Step 3: Write `config/controlled_vocabularies.yaml`**

Seed with the integration vocabularies only (spec §5.3 first group). Task 10 adds the second group. Keys are `<applies_to>.<property>`.

```yaml
# Controlled vocabularies published as ControlledVocabulary nodes.
#
# THE source of truth: adapters import their literals from here via
# multiomics_kg/utils/controlled_vocab.py, and tests assert the live graph
# matches. See docs/superpowers/specs/2026-08-16-vocabulary-contract-design.md.
#
# Required keys: applies_to, applies_to_kind, property, value_type, closed,
#                values, description
# Optional:      sparse, expected_empty, min_value, max_value,
#                signal_count, signals
#
# value_type ∈ {string, string_array, float, int, bool_string}
#   'bool' is NOT admissible — house rule R5.

# ── evidence: per edge type, domains genuinely differ (spec §5.2) ────────────
Gene_has_pfam.evidence:
  applies_to: Gene_has_pfam
  applies_to_kind: edge
  property: evidence
  value_type: string
  closed: true
  values: [curated, signature]
  description: >
    Inference strength for this annotation. 'signature' = a direct Pfam HMM
    hit. Note domain_inferred and family_inferred are NOT possible on this
    edge type.

Gene_catalyzes_ec_number.evidence:
  applies_to: Gene_catalyzes_ec_number
  applies_to_kind: edge
  property: evidence
  value_type: string
  closed: true
  values: [curated, family_inferred]
  description: >
    Inference strength. InterPro contributes EC only from FAMILY entries with
    a single EC, so domain_inferred is not possible on this edge type.

Gene_involved_in_biological_process.evidence: &go_evidence
  applies_to: Gene_involved_in_biological_process
  applies_to_kind: edge
  property: evidence
  value_type: string
  closed: true
  values: [curated, family_inferred, domain_inferred]
  description: >
    Inference strength. curated > signature > family_inferred >
    domain_inferred. InterPro contributes GO from FAMILY and DOMAIN entries.

Gene_enables_molecular_function.evidence:
  <<: *go_evidence
  applies_to: Gene_enables_molecular_function

Gene_located_in_cellular_component.evidence:
  <<: *go_evidence
  applies_to: Gene_located_in_cellular_component

Gene_has_cazy_family.evidence:
  applies_to: Gene_has_cazy_family
  applies_to_kind: edge
  property: evidence
  value_type: string
  closed: true
  values: [curated, family_inferred, domain_inferred]
  description: Inference strength for this CAZy assignment.

# ── sources: also per edge type (KG-IPT-011); values are DataSource ids ──────
Gene_has_pfam.sources:
  applies_to: Gene_has_pfam
  applies_to_kind: edge
  property: sources
  value_type: string_array
  closed: true
  values: [cyanorak, eggnog, interproscan, uniprot]
  description: >
    Who asserted this annotation. Every value is the id of a DataSource node.

Gene_catalyzes_ec_number.sources:
  applies_to: Gene_catalyzes_ec_number
  applies_to_kind: edge
  property: sources
  value_type: string_array
  closed: true
  values: [cyanorak, eggnog, interproscan, uniprot]
  description: Who asserted this annotation. Values are DataSource ids.

Gene_involved_in_biological_process.sources: &go_sources
  applies_to: Gene_involved_in_biological_process
  applies_to_kind: edge
  property: sources
  value_type: string_array
  closed: true
  values: [cyanorak, eggnog, interproscan, ncbi, uniprot]
  description: >
    Who asserted this annotation. Values are DataSource ids. ncbi contributes
    only to the GO edges.

Gene_enables_molecular_function.sources:
  <<: *go_sources
  applies_to: Gene_enables_molecular_function

Gene_located_in_cellular_component.sources:
  <<: *go_sources
  applies_to: Gene_located_in_cellular_component

Gene_has_cazy_family.sources:
  applies_to: Gene_has_cazy_family
  applies_to_kind: edge
  property: sources
  value_type: string_array
  closed: true
  values: [eggnog, interproscan]
  description: Who asserted this annotation. Values are DataSource ids.

Gene_has_tcdb_family.sources:
  applies_to: Gene_has_tcdb_family
  applies_to_kind: edge
  property: sources
  value_type: string_array
  closed: true
  values: [eggnog, tcdb_diamond]
  description: >
    Evidence sources for this transporter call. Values are DataSource ids.
    Both entries present = cross-source corroboration.

# ── evidence_score: normalized float, signals published (R4, KG-IPT-012) ─────
Gene_has_pfam.evidence_score: &score3
  applies_to: Gene_has_pfam
  applies_to_kind: edge
  property: evidence_score
  value_type: float
  closed: false
  values: []
  min_value: 0.0
  max_value: 1.0
  signal_count: 3
  signals: [multi_source, high_trust_assertion, not_domain_inferred]
  description: >
    Advisory ranking, never a filter. Fraction of independent corroborating
    signals that fired; round(score * signal_count) recovers the raw count.
    NOT calibrated against other edge types — compare within an edge type.
    Float: sort and threshold, never test equality.

Gene_catalyzes_ec_number.evidence_score:
  <<: *score3
  applies_to: Gene_catalyzes_ec_number
Gene_involved_in_biological_process.evidence_score:
  <<: *score3
  applies_to: Gene_involved_in_biological_process
Gene_enables_molecular_function.evidence_score:
  <<: *score3
  applies_to: Gene_enables_molecular_function
Gene_located_in_cellular_component.evidence_score:
  <<: *score3
  applies_to: Gene_located_in_cellular_component
Gene_has_cazy_family.evidence_score:
  <<: *score3
  applies_to: Gene_has_cazy_family

Gene_has_tcdb_family.evidence_score:
  applies_to: Gene_has_tcdb_family
  applies_to_kind: edge
  property: evidence_score
  value_type: float
  closed: false
  values: []
  min_value: 0.0
  max_value: 1.0
  signal_count: 5
  signals: [eggnog_called, agrees_across_sources, tier_le_2,
            pfam_corroborated, go_corroborated]
  description: >
    Advisory ranking, never a filter. Fraction of independent corroborating
    signals that fired; round(score * signal_count) recovers the raw count.
    NOT calibrated against other edge types. Float: sort and threshold, never
    test equality.

Gene.tcdb_evidence_score_max:
  applies_to: Gene
  applies_to_kind: node
  property: tcdb_evidence_score_max
  value_type: float
  closed: false
  values: []
  min_value: 0.0
  max_value: 1.0
  sparse: true
  description: >
    Strongest TCDB claim this gene has = max(evidence_score) over its
    Gene_has_tcdb_family edges. SPARSE BY DESIGN — absent when the gene has no
    TCDB edge, so "no transporter evidence" stays distinguishable from weak
    evidence. Use coalesce(x, -1.0) for a total order.

# ── R5 conversions: two-state facts as meaningful strings ────────────────────
Gene_has_tcdb_family.source_agreement:
  applies_to: Gene_has_tcdb_family
  applies_to_kind: edge
  property: source_agreement
  value_type: string
  closed: true
  values: [both_sources, single_source]
  description: >
    Whether eggNOG and diamond concur on this TC family. Agreement is
    HIERARCHICAL (ancestor or descendant), not exact-node — the two sources
    usually concur at different depths.

Gene_has_tcdb_family.pfam_support:
  applies_to: Gene_has_tcdb_family
  applies_to_kind: edge
  property: pfam_support
  value_type: string
  closed: true
  values: [corroborated, uncorroborated]
  description: >
    Whether a Pfam domain on this gene is curated into this TC family via
    Tcdb_family_has_pfam_domain. Read in the sound direction only.

Gene_has_tcdb_family.go_support:
  applies_to: Gene_has_tcdb_family
  applies_to_kind: edge
  property: go_support
  value_type: string
  closed: true
  values: [corroborated, uncorroborated]
  description: >
    Whether a GO term on this gene is curated onto this TC family. Not
    redundant with pfam_support — 11,565 edges carry exactly one of the two.

Tcdb_family_transports_metabolite.substrate_depth:
  applies_to: Tcdb_family_transports_metabolite
  applies_to_kind: edge
  property: substrate_depth
  value_type: string
  closed: true
  values: [most_specific, inherited]
  description: >
    A (node, substrate) fact, not a node fact — 2.A.1 can be most_specific for
    one substrate and inherited for another. most_specific = no kept child of
    this node carries the same substrate. Use this, NOT
    level_kind = 'tc_specificity', to recover leaf substrate semantics.

Gene.transport_substrate_resolution:
  applies_to: Gene
  applies_to_kind: node
  property: transport_substrate_resolution
  value_type: string
  closed: true
  sparse: true
  values: [resolved, family_inferred]
  description: >
    How specific this gene's transport substrate claim is. family_inferred =
    its deepest TC attachment lumps many substrates, so the count is
    reachability, not capability. SPARSE — absent when the gene has no TCDB
    edge, so "no evidence" stays distinguishable from a weak claim.

# ── InterPro ────────────────────────────────────────────────────────────────
InterproEntry.interpro_type:
  applies_to: InterproEntry
  applies_to_kind: node
  property: interpro_type
  value_type: string
  closed: true
  values: [family, domain, homologous_superfamily, repeat, conserved_site,
           active_site, binding_site, ptm]
  description: >
    InterPro entry class. PRIMARY stratification key for ORA over InterPro —
    breadth varies by type far more than by level.

InterproEntry.level_kind:
  applies_to: InterproEntry
  applies_to_kind: node
  property: level_kind
  value_type: string
  closed: true
  values: []
  expected_empty: true
  description: >
    InterPro depth tiers have no natural names. Emit null. Stratify ORA by
    (interpro_type, level) — interpro_type primary, level secondary.

Gene_has_interpro_entry.libraries:
  applies_to: Gene_has_interpro_entry
  applies_to_kind: edge
  property: libraries
  value_type: string_array
  closed: true
  values: [cdd, gene3d, hamap, ncbifam, panther, pfam, pirsf, prints,
           prosite_patterns, prosite_profiles, sfld, smart, superfamily]
  description: >
    InterProScan member databases that produced a match for this entry. The
    member-DB granularity of the signature-vs-inferred distinction — 'pfam'
    means a direct HMM hit.
```

- [ ] **Step 4: Write `multiomics_kg/utils/controlled_vocab.py`**

```python
"""Controlled-vocabulary contract — the single source of truth.

Adapters import their literals from here so an undeclared value cannot be
emitted; the ControlledVocabulary adapter publishes the same data into the
graph; kg-validity asserts the live graph agrees.

See docs/superpowers/specs/2026-08-16-vocabulary-contract-design.md §5.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

DEFAULT_PATH = Path("config/controlled_vocabularies.yaml")

# R5: a two-state fact is a meaningful categorical string. Native bool is not
# admissible — BioCypher does not round-trip it, and a bare `true` is
# unreadable in a result row.
VALUE_TYPES = {"string", "string_array", "float", "int", "bool_string"}

_REQUIRED = ("applies_to", "applies_to_kind", "property", "value_type",
             "closed", "description")


@dataclass(frozen=True)
class VocabEntry:
    applies_to: str
    applies_to_kind: str          # 'node' | 'edge'
    property: str
    value_type: str
    closed: bool
    description: str
    values: list[str] = field(default_factory=list)
    sparse: bool = False
    expected_empty: bool = False
    min_value: float | None = None
    max_value: float | None = None
    signal_count: int | None = None
    signals: list[str] = field(default_factory=list)

    @property
    def id(self) -> str:
        return f"{self.applies_to}.{self.property}"


def _validate(key: str, raw: dict[str, Any]) -> None:
    missing = [k for k in _REQUIRED if k not in raw]
    if missing:
        raise ValueError(
            f"controlled_vocabularies.yaml entry '{key}' is missing required "
            f"key(s): {', '.join(missing)}"
        )
    vt = raw["value_type"]
    if vt not in VALUE_TYPES:
        raise ValueError(
            f"controlled_vocabularies.yaml entry '{key}' declares "
            f"value_type='{vt}'. Allowed: {sorted(VALUE_TYPES)}. Note 'bool' "
            f"is forbidden by house rule R5 — use a meaningful two-state "
            f"string instead."
        )
    if raw["applies_to_kind"] not in ("node", "edge"):
        raise ValueError(
            f"controlled_vocabularies.yaml entry '{key}' has "
            f"applies_to_kind='{raw['applies_to_kind']}'; expected node|edge."
        )
    expected_id = f"{raw['applies_to']}.{raw['property']}"
    if key != expected_id:
        raise ValueError(
            f"controlled_vocabularies.yaml key '{key}' does not match its own "
            f"applies_to/property ('{expected_id}'). The key IS the node id."
        )


def load_vocabularies(path: str | Path = DEFAULT_PATH) -> dict[str, VocabEntry]:
    """Parse the YAML into VocabEntry objects keyed by '<applies_to>.<property>'."""
    with open(path) as f:
        raw_all = yaml.safe_load(f) or {}
    out: dict[str, VocabEntry] = {}
    for key, raw in raw_all.items():
        _validate(key, raw)
        out[key] = VocabEntry(
            applies_to=raw["applies_to"],
            applies_to_kind=raw["applies_to_kind"],
            property=raw["property"],
            value_type=raw["value_type"],
            closed=bool(raw["closed"]),
            description=" ".join(raw["description"].split()),
            values=list(raw.get("values") or []),
            sparse=bool(raw.get("sparse", False)),
            expected_empty=bool(raw.get("expected_empty", False)),
            min_value=raw.get("min_value"),
            max_value=raw.get("max_value"),
            signal_count=raw.get("signal_count"),
            signals=list(raw.get("signals") or []),
        )
    return out


def vocabularies_hash(entries: list[VocabEntry]) -> str:
    """Stable sha256 over the vocabulary set. Order-independent."""
    payload = sorted(
        {
            "id": e.id, "value_type": e.value_type, "closed": e.closed,
            "values": sorted(e.values), "sparse": e.sparse,
            "expected_empty": e.expected_empty, "min_value": e.min_value,
            "max_value": e.max_value, "signal_count": e.signal_count,
            "signals": sorted(e.signals),
        }
        for e in entries
    ) if False else sorted(
        (
            json.dumps({
                "id": e.id, "value_type": e.value_type, "closed": e.closed,
                "values": sorted(e.values), "sparse": e.sparse,
                "expected_empty": e.expected_empty, "min_value": e.min_value,
                "max_value": e.max_value, "signal_count": e.signal_count,
                "signals": sorted(e.signals),
            }, sort_keys=True)
            for e in entries
        )
    )
    digest = hashlib.sha256("\n".join(payload).encode()).hexdigest()
    return f"sha256:{digest}"


class _Vocab:
    """Lazy singleton so importing an adapter does not read the YAML."""

    def __init__(self, path: str | Path = DEFAULT_PATH):
        self._path = path
        self._entries: dict[str, VocabEntry] | None = None

    def _load(self) -> dict[str, VocabEntry]:
        if self._entries is None:
            self._entries = load_vocabularies(self._path)
        return self._entries

    def entries(self) -> list[VocabEntry]:
        return list(self._load().values())

    def get(self, applies_to: str, prop: str) -> VocabEntry | None:
        return self._load().get(f"{applies_to}.{prop}")

    def check(self, applies_to: str, prop: str, value: str) -> str:
        """Return *value*, raising if it is not declared for a closed vocabulary.

        Undeclared (applies_to, property) pairs pass through — the contract is
        seeded incrementally and this must not block unrelated adapters.
        """
        entry = self.get(applies_to, prop)
        if entry is None or not entry.closed:
            return value
        if value not in entry.values:
            raise ValueError(
                f"'{value}' is not a declared value of {entry.id}. "
                f"Declared: {sorted(entry.values)}. Add it to "
                f"config/controlled_vocabularies.yaml, or fix the emitter."
            )
        return value


VOCAB = _Vocab()
```

Simplify `vocabularies_hash` before committing — the `if False else` above is a drafting artifact. The body should be only the second branch:

```python
def vocabularies_hash(entries: list[VocabEntry]) -> str:
    """Stable sha256 over the vocabulary set. Order-independent."""
    payload = sorted(
        json.dumps({
            "id": e.id, "value_type": e.value_type, "closed": e.closed,
            "values": sorted(e.values), "sparse": e.sparse,
            "expected_empty": e.expected_empty, "min_value": e.min_value,
            "max_value": e.max_value, "signal_count": e.signal_count,
            "signals": sorted(e.signals),
        }, sort_keys=True)
        for e in entries
    )
    digest = hashlib.sha256("\n".join(payload).encode()).hexdigest()
    return f"sha256:{digest}"
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_controlled_vocab.py -v`
Expected: 7 passed. If `test_check_passes_through_open_vocabularies` fails, confirm `Gene.gene_category` is absent from the YAML at this stage — pass-through is the correct behaviour for an unseeded pair.

- [ ] **Step 6: Commit**

```bash
git add config/controlled_vocabularies.yaml multiomics_kg/utils/controlled_vocab.py tests/test_controlled_vocab.py
git commit -m "feat(vocab): controlled-vocabulary config + loader

Single source of truth for declared vocabularies. Adapters import their
literals from here so an undeclared value cannot be emitted. value_type
forbids native bool per house rule R5."
```

---

## Task 2: `ControlledVocabulary` adapter, schema, registration

**Files:**
- Create: `multiomics_kg/adapters/controlled_vocabulary_adapter.py`
- Modify: `config/schema_config.yaml` (add node type near the `data source` block, ~line 476)
- Modify: `create_knowledge_graph.py` (beside the `DataSourceAdapter` block, ~line 85)
- Test: `tests/test_controlled_vocabulary_adapter.py`

**Interfaces:**
- Consumes: `VOCAB.entries()`, `VocabEntry` from Task 1
- Produces: `ControlledVocabularyAdapter` with `download_data(cache: bool = True) -> None` and `get_nodes() -> Iterator[tuple[str, str, dict]]` yielding `(node_id, "controlled_vocabulary", properties)`; node ids are the bare `<applies_to>.<property>` key

- [ ] **Step 1: Write the failing test**

```python
# tests/test_controlled_vocabulary_adapter.py
from multiomics_kg.adapters.controlled_vocabulary_adapter import (
    ControlledVocabularyAdapter,
)


def _nodes():
    a = ControlledVocabularyAdapter()
    a.download_data()
    return {nid: props for nid, _label, props in a.get_nodes()}


def test_emits_one_node_per_vocabulary():
    nodes = _nodes()
    assert "Gene_catalyzes_ec_number.evidence" in nodes
    assert "InterproEntry.interpro_type" in nodes


def test_label_is_controlled_vocabulary():
    a = ControlledVocabularyAdapter()
    a.download_data()
    assert {label for _n, label, _p in a.get_nodes()} == {"controlled_vocabulary"}


def test_node_carries_the_declared_values():
    props = _nodes()["InterproEntry.interpro_type"]
    assert props["applies_to"] == "InterproEntry"
    assert props["applies_to_kind"] == "node"
    assert props["value_type"] == "string"
    assert "homologous_superfamily" in props["values"]
    assert props["closed"] == "true"          # R5: bool_string, not bool


def test_expected_empty_vocabulary_emits_empty_values():
    props = _nodes()["InterproEntry.level_kind"]
    assert props["values"] == []
    assert props["expected_empty"] == "true"


def test_score_vocabulary_publishes_its_signals():
    props = _nodes()["Gene_has_tcdb_family.evidence_score"]
    assert props["signal_count"] == 5
    assert "pfam_corroborated" in props["signals"]
    assert props["min_value"] == 0.0 and props["max_value"] == 1.0


def test_strings_are_sanitized():
    """No property value may contain ' or | (BioCypher CSV + array delimiter)."""
    for props in _nodes().values():
        for v in props.values():
            for s in (v if isinstance(v, list) else [v]):
                if isinstance(s, str):
                    assert "'" not in s and "|" not in s
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_controlled_vocabulary_adapter.py -v`
Expected: FAIL — `ModuleNotFoundError`

- [ ] **Step 3: Write the adapter**

```python
"""ControlledVocabulary adapter — publishes config/controlled_vocabularies.yaml
as queryable nodes so consumers stop hard-coding value sets.

Node-only, modelled on data_source_adapter.py. Booleans are emitted as
"true"/"false" strings per house rule R5.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterator

from multiomics_kg.utils.controlled_vocab import DEFAULT_PATH, load_vocabularies


def _clean_str(value: str) -> str:
    """Sanitize string for BioCypher CSV output (avoid ' and |)."""
    return value.replace("'", "^").replace("|", "")


def _b(value: bool) -> str:
    """R5: adapter-emitted booleans are meaningful strings, never native bool."""
    return "true" if value else "false"


class ControlledVocabularyAdapter:
    def __init__(self, config_path: str | Path = DEFAULT_PATH):
        self.config_path = Path(config_path)
        self._entries: list[Any] = []

    def download_data(self, cache: bool = True) -> None:
        self._entries = list(load_vocabularies(self.config_path).values())

    def get_nodes(self) -> Iterator[tuple[str, str, dict[str, Any]]]:
        if not self._entries:
            self.download_data()
        for e in self._entries:
            props: dict[str, Any] = {
                "id": _clean_str(e.id),
                "applies_to": _clean_str(e.applies_to),
                "applies_to_kind": _clean_str(e.applies_to_kind),
                "property": _clean_str(e.property),
                "value_type": _clean_str(e.value_type),
                "closed": _b(e.closed),
                "sparse": _b(e.sparse),
                "expected_empty": _b(e.expected_empty),
                "values": [_clean_str(v) for v in e.values],
                "description": _clean_str(e.description),
            }
            # Sparse numeric/score metadata — omitted rather than nulled so the
            # absence is meaningful (a string vocabulary has no min_value).
            if e.min_value is not None:
                props["min_value"] = float(e.min_value)
            if e.max_value is not None:
                props["max_value"] = float(e.max_value)
            if e.signal_count is not None:
                props["signal_count"] = int(e.signal_count)
            if e.signals:
                props["signals"] = [_clean_str(s) for s in e.signals]
            yield e.id, "controlled_vocabulary", props
```

- [ ] **Step 4: Add the schema node type**

In `config/schema_config.yaml`, immediately after the `data source:` block (ends ~line 489):

```yaml
# Controlled vocabularies published as data (spec 2026-08-16 §5). One node per
# (label, property) pair. `evidence` and `sources` get one node PER EDGE TYPE
# because their domains genuinely differ — a flat union would offer a filter
# value that can never match on that edge.
controlled vocabulary:
  is_a: information content entity
  represented_as: node
  preferred_id: controlled_vocabulary
  label_in_input: controlled_vocabulary
  properties:
    id: str                # '<applies_to>.<property>'
    applies_to: str        # node label or edge type
    applies_to_kind: str   # 'node' | 'edge'
    property: str
    value_type: str        # string | string_array | float | int | bool_string
    closed: str            # "true" | "false"  — false = enumerate at runtime
    sparse: str            # "true" | "false"  — property may be absent
    expected_empty: str    # "true" | "false"  — declared, no values by design
    values: str[]
    description: str
    min_value: float       # sparse: numeric vocabularies only
    max_value: float       # sparse
    signal_count: int      # sparse: evidence_score only
    signals: str[]         # sparse: evidence_score component names
```

- [ ] **Step 5: Register the adapter**

In `create_knowledge_graph.py`, after the `DataSourceAdapter` block (~line 89):

```python
from multiomics_kg.adapters.controlled_vocabulary_adapter import (
    ControlledVocabularyAdapter,
)

    # ControlledVocabulary nodes — the machine-readable vocabulary contract.
    # Consumers (MCP/explorer) read these instead of hard-coding value sets.
    controlled_vocab_adapter = ControlledVocabularyAdapter()
    controlled_vocab_adapter.download_data()
    bc.write_nodes(controlled_vocab_adapter.get_nodes())
```

- [ ] **Step 6: Run tests**

Run: `pytest tests/test_controlled_vocabulary_adapter.py -v && pytest -m "not slow and not kg" -q`
Expected: 6 passed in the new file; no regressions elsewhere.

- [ ] **Step 7: Commit**

```bash
git add multiomics_kg/adapters/controlled_vocabulary_adapter.py config/schema_config.yaml create_knowledge_graph.py tests/test_controlled_vocabulary_adapter.py
git commit -m "feat(vocab): ControlledVocabulary nodes

Publishes the vocabulary contract as queryable nodes, following the
DataSource adapter pattern. evidence and sources get one node per edge type
so a filter value can never be offered where it cannot match."
```

---

## Task 3: R1 — lowercase `interpro_type` and `libraries`

**Files:**
- Modify: `multiomics_kg/utils/interproscan.py:97` (`interpro_type`), `:124`/`:130` (`libraries`)
- Modify: `multiomics_kg/adapters/interpro_adapter.py:313` (emit), `:388`/`:404` (the `etype != "FAMILY"` guards)
- Test: `tests/test_interpro_adapter.py` (existing — check for the file first; if absent create it)

**Interfaces:**
- Consumes: `VOCAB.check` from Task 1
- Produces: `InterproEntry.interpro_type` and `Gene_has_interpro_entry.libraries` values in lowercase `snake_case`

**Critical:** [interpro_adapter.py:388](multiomics_kg/adapters/interpro_adapter.py#L388) and `:404` compare `etype != "FAMILY"`. Those comparisons must be lowercased **in the same change**, or the Layer-A gating inverts. Task 7 deletes the property those guards feed, so if Task 7 is done first this step is moot — do them in plan order and the guards disappear cleanly.

- [ ] **Step 1: Write the failing test**

```python
def test_interpro_type_is_lowercase_snake_case():
    from multiomics_kg.utils.interproscan import normalize_interpro_type
    assert normalize_interpro_type("HOMOLOGOUS_SUPERFAMILY") == "homologous_superfamily"
    assert normalize_interpro_type("Family") == "family"
    assert normalize_interpro_type("") == ""


def test_libraries_are_lowercase_snake_case():
    from multiomics_kg.utils.interproscan import normalize_library
    assert normalize_library("PROSITE_PATTERNS") == "prosite_patterns"
    assert normalize_library("GENE3D") == "gene3d"


def test_declared_vocabulary_matches_the_normalizer():
    """Every value the normalizer can produce must be declared."""
    from multiomics_kg.utils.controlled_vocab import VOCAB
    from multiomics_kg.utils.interproscan import normalize_interpro_type
    for raw in ["FAMILY", "DOMAIN", "HOMOLOGOUS_SUPERFAMILY", "REPEAT",
                "CONSERVED_SITE", "ACTIVE_SITE", "BINDING_SITE", "PTM"]:
        VOCAB.check("InterproEntry", "interpro_type", normalize_interpro_type(raw))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_interpro_adapter.py -k lowercase -v`
Expected: FAIL — `ImportError: cannot import name 'normalize_interpro_type'`

- [ ] **Step 3: Add the normalizers in `multiomics_kg/utils/interproscan.py`**

```python
def normalize_interpro_type(raw: str | None) -> str:
    """House rule R1 — lowercase snake_case. InterPro ships UPPERCASE."""
    return (raw or "").strip().lower()


def normalize_library(raw: str | None) -> str:
    """House rule R1 — lowercase snake_case member-DB name."""
    return (raw or "").strip().lower()
```

Apply at line 97: `"interpro_type": normalize_interpro_type(ipr_type),`
and at line 124: `libraries = sorted({normalize_library(m["library"]) for m in matches if m["library"]})`

- [ ] **Step 4: Update the adapter emission and guards**

`interpro_adapter.py:313` → `"interpro_type": _clean_str(normalize_interpro_type(ref.get("type"))),`
`interpro_adapter.py:388` and `:404` → the comparison becomes `etype != "family"` and `etype` is built with `normalize_interpro_type(ref.get("type"))` rather than `.upper()`.

- [ ] **Step 5: Run tests**

Run: `pytest tests/test_interpro_adapter.py -v && pytest -m "not slow and not kg" -q`
Expected: PASS. Grep for stragglers: `grep -rn '"FAMILY"\|PROSITE_PATTERNS' multiomics_kg/` should return nothing outside comments.

- [ ] **Step 6: Commit**

```bash
git add multiomics_kg/utils/interproscan.py multiomics_kg/adapters/interpro_adapter.py tests/test_interpro_adapter.py
git commit -m "refactor(interpro): lowercase interpro_type and libraries (R1)

Also lowercases the etype != FAMILY guards in the same change — leaving them
uppercase would invert Layer-A gating."
```

---

## Task 4: R2 — `sources` values become `DataSource` ids

**Files:**
- Modify: `multiomics_kg/download/build_gene_annotations.py:335` (`srcs.add("interpro")`)
- Modify: `multiomics_kg/adapters/tcdb_adapter.py:157`, `:168`, `:356` (the `"diamond"` literal)
- Modify: `multiomics_kg/utils/annotation_provenance.py:38` (the `{"eggnog", "interpro"}` dependent-pair collapse)
- Test: `tests/test_annotation_provenance.py` (existing), `tests/test_tcdb_adapter.py` (existing)

**Interfaces:**
- Produces: `sources` arrays whose every value is a `DataSource.id`

- [ ] **Step 1: Write the failing test**

```python
def test_source_values_are_all_data_source_ids():
    """R2: every declared sources value must be an id in gene_annotations_config."""
    import yaml
    from multiomics_kg.utils.controlled_vocab import load_vocabularies
    cfg = yaml.safe_load(open("config/gene_annotations_config.yaml"))
    ds_ids = {ls["id"]
              for src in cfg["sources"].values()
              for ls in src["logical_sources"]}
    for entry in load_vocabularies().values():
        if entry.property == "sources":
            undeclared = set(entry.values) - ds_ids
            assert not undeclared, (
                f"{entry.id} declares source value(s) {sorted(undeclared)} with "
                f"no matching DataSource id. Known ids: {sorted(ds_ids)}")


def test_interpro_source_label_is_interproscan():
    from multiomics_kg.utils.annotation_provenance import annotation_edge_props
    gene = {"go_terms_source": {"GO:1": ["interproscan"]},
            "go_terms_evidence": {"GO:1": "domain_inferred"}}
    props = annotation_edge_props(gene, "go_terms", "GO:1")
    assert props["sources"] == ["interproscan"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_annotation_provenance.py -k source -v`
Expected: FAIL — the config still emits `interpro`.

- [ ] **Step 3: Apply the renames**

- `build_gene_annotations.py:335` → `srcs.add("interproscan")`
- `annotation_provenance.py` → the dependent-pair collapse becomes `{"eggnog", "interproscan"}` and discards `"interproscan"`; update its docstring's source list to `ncbi|cyanorak|uniprot|eggnog|interproscan`
- `tcdb_adapter.py` → replace the three `"diamond"` literals with `"tcdb_diamond"`; update the docstring examples at `:132` and `:153`

- [ ] **Step 4: Run tests**

Run: `pytest -m "not slow and not kg" -q`
Expected: PASS. `grep -rn '"diamond"' multiomics_kg/` returns nothing.

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "refactor(provenance): sources values are DataSource ids (R2)

interpro -> interproscan, diamond -> tcdb_diamond. Provenance becomes
joinable against the DataSource nodes and the three-spellings problem goes
away."
```

---

## Task 5: R4 — normalize `evidence_score` to a float in [0,1]

**Files:**
- Modify: `multiomics_kg/utils/annotation_provenance.py` (the 3-signal score)
- Modify: `scripts/post-import.cypher:1188` and `:1212-1213`, plus the matching block in `scripts/post-import.sh`
- Modify: `config/schema_config.yaml:1302` (`tcdb_evidence_score: int` → `evidence_score: float`)
- Test: `tests/test_annotation_provenance.py`

**Interfaces:**
- Produces: `evidence_score: float` on all six gene→ontology edge types and on `Gene_has_tcdb_family`; `Gene.tcdb_evidence_score_max: float`

- [ ] **Step 1: Write the failing test**

```python
import pytest
from multiomics_kg.utils.annotation_provenance import annotation_edge_props


@pytest.mark.parametrize("sources,evidence,expected", [
    (["uniprot", "eggnog"], "curated",         1.0),    # 3/3
    (["eggnog"],            "curated",         0.667),  # 2/3
    (["eggnog"],            "family_inferred", 0.333),  # 1/3
    (["interproscan"],      "domain_inferred", 0.0),    # 0/3
])
def test_evidence_score_is_normalized_to_unit_interval(sources, evidence, expected):
    gene = {"go_terms_source": {"GO:1": sources},
            "go_terms_evidence": {"GO:1": evidence}}
    assert annotation_edge_props(gene, "go_terms", "GO:1")["evidence_score"] == expected


def test_score_is_rounded_to_three_decimals():
    gene = {"go_terms_source": {"GO:1": ["eggnog"]},
            "go_terms_evidence": {"GO:1": "family_inferred"}}
    score = annotation_edge_props(gene, "go_terms", "GO:1")["evidence_score"]
    assert isinstance(score, float)
    assert score == round(score, 3)


def test_round_recovers_the_raw_signal_count():
    """KG-IPT-012: round, never truncate — 0.333 * 3 = 0.999."""
    assert round(0.333 * 3) == 1
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_annotation_provenance.py -k normalized -v`
Expected: FAIL — score is currently `int` `0..3`.

- [ ] **Step 3: Normalize in `annotation_provenance.py`**

Replace the score assembly (keep the three `+= 1` signal blocks unchanged):

```python
    _SIGNAL_COUNT = 3   # module constant, mirrored in controlled_vocabularies.yaml

    # ... the three existing `score += 1` blocks are untouched ...

    props: dict = {
        "evidence": evidence,
        # R4: normalized to [0,1] so a Pfam score and a TCDB score are not
        # arithmetically comparable-but-wrong. round(x * signal_count)
        # recovers the raw count; NEVER truncate (0.333 * 3 = 0.999).
        "evidence_score": round(score / _SIGNAL_COUNT, 3),
    }
```

Update the module docstring: `evidence_score` is now "float in [0,1], advisory, never a filter; a ready sort key. Multiply by 3 and round to recover the signal count."

- [ ] **Step 4: Normalize the TCDB score in post-import**

In **both** `scripts/post-import.cypher` and `scripts/post-import.sh`, at the `r.tcdb_evidence_score =` assignment (~line 1188), rename the property and divide by 5:

```cypher
  SET r.source_agreement = ...,       // Task 6 sets these; keep the same block
      r.evidence_score = round(
        ( (CASE WHEN 'eggnog' IN r.sources THEN 1 ELSE 0 END)
        + (CASE WHEN agree   THEN 1 ELSE 0 END)
        + (CASE WHEN coalesce(r.tier, 99) <= 2 THEN 1 ELSE 0 END)
        + (CASE WHEN pfam_ok THEN 1 ELSE 0 END)
        + (CASE WHEN go_ok   THEN 1 ELSE 0 END) ) / 5.0, 3)
```

and at ~line 1212:

```cypher
  WITH g, max(r.evidence_score) AS best
  SET g.tcdb_evidence_score_max = best
```

Update the comment above it: the `coalesce(..., -1)` guidance becomes `coalesce(..., -1.0)` because `0.0` is now a legitimate value.

- [ ] **Step 5: Update the schema**

`config/schema_config.yaml:1302`: `tcdb_evidence_score: int` → `evidence_score: float`, and update the comment block above it from "(0-5)" to "(float 0-1; 5 signals, see controlled_vocabularies.yaml)".

- [ ] **Step 6: Run tests**

Run: `pytest -m "not slow and not kg" -q`
Expected: PASS. Then `diff <(grep -c . scripts/post-import.cypher) <(grep -c . scripts/post-import.sh)` is not a valid equality check (the `.sh` has bash wrapping) — instead visually confirm the two edited Cypher blocks are character-identical.

- [ ] **Step 7: Commit**

```bash
git add -u
git commit -m "refactor(provenance): evidence_score normalized to float [0,1] (R4)

One name, one scale across every annotation edge. Two integer scales under one
property name let an unread-contract consumer compare a Pfam 3 to a TCDB 3 and
be wrong by a factor; normalized they degrade to a comparable fraction. The
raw count is recoverable via round(score * signal_count), published in the
contract."
```

---

## Task 6: R5 — convert the remaining booleans to meaningful strings

**Files:**
- Modify: `scripts/post-import.{cypher,sh}` — `agrees_across_sources` / `pfam_corroborated` / `go_corroborated`
- Modify: `multiomics_kg/adapters/tcdb_adapter.py:442` — `substrate_depth` values
- Modify: `config/schema_config.yaml:1303-1305`, `:1329`
- Test: `tests/test_tcdb_adapter.py`

**Interfaces:**
- Produces: `source_agreement: both_sources|single_source`, `pfam_support`/`go_support: corroborated|uncorroborated`, `substrate_depth: most_specific|inherited`

- [ ] **Step 1: Write the failing test**

```python
def test_substrate_depth_uses_meaningful_values():
    from multiomics_kg.utils.controlled_vocab import VOCAB
    entry = VOCAB.get("Tcdb_family_transports_metabolite", "substrate_depth")
    assert set(entry.values) == {"most_specific", "inherited"}
    # and the adapter emits only declared values
    for v in ("most_specific", "inherited"):
        VOCAB.check("Tcdb_family_transports_metabolite", "substrate_depth", v)


def test_no_schema_property_is_declared_bool():
    """R5: native bool is forbidden graph-wide."""
    import re, pathlib
    text = pathlib.Path("config/schema_config.yaml").read_text()
    offenders = re.findall(r"^\s+([a-z_]+): bool\s*$", text, re.M)
    assert offenders == [], f"native bool properties remain: {offenders}"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_tcdb_adapter.py -k meaningful -v`
Expected: FAIL — five `: bool` properties remain in the schema.

- [ ] **Step 3: Convert `substrate_depth` in the adapter**

`tcdb_adapter.py:442`:

```python
                    {"substrate_depth": VOCAB.check(
                        "Tcdb_family_transports_metabolite", "substrate_depth",
                        "most_specific" if is_deepest else "inherited")},
```

Update the comment at `:415` and the `_compute_substrate_depth` docstring: "deepest" becomes "most_specific" throughout.

- [ ] **Step 4: Convert the three TCDB booleans in post-import**

In **both** scripts, the `SET` block at ~line 1185:

```cypher
  SET r.source_agreement = CASE WHEN agree   THEN 'both_sources' ELSE 'single_source' END,
      r.pfam_support     = CASE WHEN pfam_ok THEN 'corroborated' ELSE 'uncorroborated' END,
      r.go_support       = CASE WHEN go_ok   THEN 'corroborated' ELSE 'uncorroborated' END,
      r.evidence_score   = round( ... )   // from Task 5, unchanged here
```

The local variables `agree` / `pfam_ok` / `go_ok` stay native booleans — they are Cypher locals, never stored.

- [ ] **Step 5: Update the schema**

`config/schema_config.yaml:1303-1305`:

```yaml
    # R5: two-state facts are meaningful strings, never native bool.
    source_agreement: str    # both_sources | single_source (HIERARCHICAL agreement)
    pfam_support: str        # corroborated | uncorroborated
    go_support: str          # corroborated | uncorroborated
```

`:1329` keeps `substrate_depth: str`; update its comment to `most_specific | inherited`.

- [ ] **Step 6: Run tests**

Run: `pytest -m "not slow and not kg" -q`
Expected: PASS — including `test_no_schema_property_is_declared_bool`, because Task 7 has not run yet but the two `ambiguous: bool` entries are still present. **If that test fails on `ambiguous`, that is expected at this point** — either reorder to run Task 7 first, or mark this assertion `xfail` with a comment naming Task 7 and flip it there. Prefer running Task 7 immediately after.

- [ ] **Step 7: Commit**

```bash
git add -u
git commit -m "refactor(tcdb): two-state facts become meaningful strings (R5)

agrees_across_sources/pfam_corroborated/go_corroborated ->
source_agreement/pfam_support/go_support; substrate_depth deepest|ancestor ->
most_specific|inherited. ancestor described the node; the property describes a
(node, substrate) fact."
```

---

## Task 7: R3 + rev 4/5 — delete four properties

**Files:**
- Modify: `multiomics_kg/adapters/interpro_adapter.py:363-420` (Layer A)
- Modify: `scripts/post-import.{cypher,sh}` — lines ~922-953 (TCDB `is_promiscuous`), ~984-988 (InterPro `is_promiscuous`), ~1131 (the one read)
- Modify: `config/schema_config.yaml:1532`, `:1543`
- Test: `tests/test_interpro_adapter.py`, `tests/test_tcdb_adapter.py`

**Interfaces:**
- Produces: Layer-A router edges with **no properties**; `TcdbFamily` and `InterproEntry` without `is_promiscuous`

**Rationale (spec §3 R3 + §9.8):** all four restate something already in the graph — an adjacent node's property, a stored count, or a constant. Materialize traversals, not predicates.

- [ ] **Step 1: Write the failing test**

```python
def test_layer_a_edges_carry_no_properties():
    """rev 4/5: ambiguous is derivable and source_db is a constant."""
    import pathlib
    src = pathlib.Path("multiomics_kg/adapters/interpro_adapter.py").read_text()
    assert '"ambiguous"' not in src
    assert '"source_db"' not in src


def test_is_promiscuous_is_gone_from_post_import():
    import pathlib
    for p in ("scripts/post-import.cypher", "scripts/post-import.sh"):
        text = pathlib.Path(p).read_text()
        assert "is_promiscuous" not in text, f"{p} still sets or reads is_promiscuous"


def test_transport_substrate_resolution_threshold_is_inlined():
    import pathlib
    text = pathlib.Path("scripts/post-import.cypher").read_text()
    assert "metabolite_count, 0) >= 50" in text, (
        "the breadth threshold must be inlined into transport_substrate_resolution")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_interpro_adapter.py -k layer_a -v`
Expected: FAIL — both literals still present.

- [ ] **Step 3: Delete the Layer-A properties**

In `interpro_adapter.py`, both `yield` blocks (~line 399 and ~417) emit `{}` instead of the property dict, and the `amb = ...` lines are removed:

```python
            if ecs:
                for raw in ecs:
                    for norm in _normalize_ec_token(raw):
                        nid = _ec_node_id(norm)
                        if nid not in observed_ec:
                            continue
                        yield (
                            f"{acc}-related_ec-{norm}",
                            _interpro_node_id(acc),
                            nid,
                            "interpro_entry_related_to_ec_number",
                            {},
                        )
                        la_ec += 1
```

Apply the same shape to the CAZy block. Update the module docstring (~line 22): the `ambiguous` flag reference goes; state instead that the router carries no properties, and that a consumer derives the old flag from `count(r) > 1 OR n.interpro_type <> 'family'`.

`etype` becomes unused in this function — remove the local, and with it the `normalize_interpro_type` call added in Task 3 at that site (the node-emission call at `:313` stays).

- [ ] **Step 4: Delete the two schema property blocks**

`config/schema_config.yaml:1532` and `:1543` — remove `ambiguous: bool` and `source_db: str` from both `interpro entry related to ec number` and `interpro entry related to cazy family`, leaving each with no `properties:` block. Update the comment above (~line 1518) to say the routers carry no properties and why.

- [ ] **Step 5: Delete both `is_promiscuous` blocks and inline the threshold**

In **both** post-import scripts:
- Delete the comment + `SET t.is_promiscuous = ...` block (~lines 922-953) and the InterPro one (~984-988).
- At ~line 1131, replace the read:

```cypher
  WITH g,
       count(DISTINCT m_tr) AS tr_met_count,
       count(DISTINCT t) AS n_deepest,
       // Breadth threshold inlined (spec §3 R3): the old TcdbFamily
       // .is_promiscuous restated a predicate over metabolite_count, which the
       // node already publishes. Consumers apply their own cutoff to the count;
       // this is the KG's, and it lives in exactly one place.
       collect(DISTINCT (coalesce(t.level, 0) >= 2
                         AND coalesce(t.metabolite_count, 0) >= 50)) AS breadth
```

- [ ] **Step 6: Run tests**

Run: `pytest -m "not slow and not kg" -q`
Expected: PASS, including `test_no_schema_property_is_declared_bool` from Task 6 — the last two `: bool` entries are now gone.

- [ ] **Step 7: Commit**

```bash
git add -u
git commit -m "refactor: delete four derivable properties (R3, rev 4/5)

Layer-A ambiguous (both arms derivable; residue pre-pruning and arguably
wrong), Layer-A source_db (a hardcoded literal), and both is_promiscuous
flags (thresholds over counts the node publishes). Layer-A edges now carry no
properties, which is right for a pure router. The one internal read inlines
the threshold into transport_substrate_resolution.

Materialize traversals, not predicates."
```

---

## Task 8: Enforcement gates

**Files:**
- Modify: `tests/test_create_knowledge_graph.py` (parameterize over build mode)
- Create: `tests/kg_validity/test_controlled_vocabularies.py`

**Interfaces:**
- Consumes: `load_vocabularies` from Task 1; the `run_query` fixture already used across `tests/kg_validity/`

- [ ] **Step 1: Rewrite `tests/test_create_knowledge_graph.py`**

```python
"""Build-level gates for create_knowledge_graph.py.

Two modes share one helper: `--test` (100 items/adapter) runs by default so a
broken build surfaces in `pytest -m "not slow and not kg"`; the full build stays
behind @pytest.mark.slow. Both need a populated cache/ — neither is portable to
a fresh checkout.
"""

import subprocess
import sys
from pathlib import Path

import pytest

from multiomics_kg.utils.controlled_vocab import load_vocabularies

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _run_build(tmp_path: Path, test_mode: bool) -> tuple[Path, str]:
    """Run the build once; return (output_dir, combined stdout+stderr)."""
    out = tmp_path / "kg"
    cmd = [sys.executable, "create_knowledge_graph.py", "--output-dir", str(out)]
    if test_mode:
        cmd.append("--test")
    result = subprocess.run(
        cmd, cwd=str(PROJECT_ROOT), capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=3600,
    )
    combined = (result.stdout or "") + "\n" + (result.stderr or "")
    assert result.returncode == 0, (
        f"Build exited {result.returncode}.\nstderr:\n{(result.stderr or '')[-2000:]}")
    return out, combined


@pytest.fixture(scope="module")
def build_test_mode(tmp_path_factory):
    return _run_build(tmp_path_factory.mktemp("build_test"), test_mode=True)


@pytest.fixture(scope="module")
def build_full(tmp_path_factory):
    return _run_build(tmp_path_factory.mktemp("build_full"), test_mode=False)


def _assert_no_error_lines(combined: str) -> None:
    errors = [l for l in combined.splitlines() if l.strip().startswith("ERROR")]
    assert errors == [], f"Found {len(errors)} ERROR line(s):\n" + "\n".join(errors)


def _observed_values(out_dir: Path) -> dict[tuple[str, str], set[str]]:
    """Scan BioCypher CSVs -> {(label_or_edge_type, property): {values}}.

    BioCypher writes <Label>-part000.csv with a matching <Label>-header.csv.
    Array properties are pipe-delimited (biocypher_config.yaml).
    """
    import csv
    observed: dict[tuple[str, str], set[str]] = {}
    for header_file in out_dir.rglob("*-header.csv"):
        label = header_file.name.replace("-header.csv", "")
        cols = [c.split(":")[0] for c in
                next(csv.reader(open(header_file)))]
        for part in header_file.parent.glob(f"{label}-part*.csv"):
            with open(part) as f:
                for row in csv.reader(f, delimiter=";"):
                    for col, raw in zip(cols, row):
                        if not raw:
                            continue
                        for v in raw.split("|"):
                            observed.setdefault((label, col), set()).add(v.strip("'"))
    return observed


def _check_vocabularies(out_dir: Path, both_directions: bool) -> None:
    declared = load_vocabularies()
    observed = _observed_values(out_dir)
    problems = []
    for entry in declared.values():
        if not entry.closed:
            continue
        seen = observed.get((entry.applies_to, entry.property), set())
        undeclared = seen - set(entry.values)
        if undeclared:
            problems.append(
                f"{entry.id}: emitted undeclared value(s) {sorted(undeclared)}; "
                f"declared {sorted(entry.values)}")
        if both_directions and not entry.expected_empty:
            unseen = set(entry.values) - seen
            if unseen and seen:
                problems.append(
                    f"{entry.id}: declared but never emitted: {sorted(unseen)}")
    assert not problems, "Vocabulary drift:\n" + "\n".join(problems)


def test_build_no_errors_test_mode(build_test_mode):
    _assert_no_error_lines(build_test_mode[1])


def test_vocabularies_in_csvs_test_mode(build_test_mode):
    """observed <= declared only: 100 items/adapter cannot contain rare values."""
    _check_vocabularies(build_test_mode[0], both_directions=False)


@pytest.mark.slow
def test_build_no_errors_full(build_full):
    _assert_no_error_lines(build_full[1])


@pytest.mark.slow
def test_vocabularies_in_csvs_full(build_full):
    _check_vocabularies(build_full[0], both_directions=True)
```

- [x] **Step 2: Measure `--test` wall-clock (spec open item 1)**

Run: `time uv run python create_knowledge_graph.py --test --output-dir /tmp/kgtest`
Record the result **in this plan file** under Task 8. If it exceeds ~3 minutes, move the two test-mode tests behind a new `@pytest.mark.build` marker, register it in `pyproject.toml`, and change the documented default command to `pytest -m "not slow and not kg and not build"`.

**Measured:** `real 2m49.965s` (run against warm `cache/`); a repeat measurement via `pytest tests/test_create_knowledge_graph.py -k test_mode` (which also runs `--test` once) took `2m43.396s`. Both are under the ~3-minute threshold, so **no `@pytest.mark.build` marker was introduced** -- `test_build_no_errors_test_mode` and `test_vocabularies_in_csvs_test_mode` stay in the default `pytest -m "not slow and not kg"` run. This is close to the threshold (within ~10-15s of 3 minutes); worth revisiting if the build grows slower.

Also discovered during this measurement: `create_knowledge_graph.py --output-dir` does **not** control where BioCypher writes its main node/edge CSVs -- it only affects the optional `--go`/`--ec` exports. The real CSVs always land in `<repo>/biocypher-out/<timestamp>/` (BioCypher's own default `output_directory`, since neither `config/biocypher_config.yaml` nor the script overrides it). `tests/test_create_knowledge_graph.py::_run_build` discovers the real directory by diffing `biocypher-out/` before/after the subprocess run rather than trusting `--output-dir`.

Also discovered: `cache/data/*/genomes/*/gene_annotations_merged.json` was stale relative to Task 4's `interpro`->`interproscan` / `diamond`->`tcdb_diamond` renames (cache built before that commit). The CSV gate caught this immediately as real drift (`Gene_has_pfam.sources: undeclared ['interpro']`, etc.) -- proof the gate is not vacuous, independent of the deliberate-tamper test described in the Task 8 report. Fixed by rebuilding the caches: `bash scripts/prepare_data.sh --steps 2 --force` (~2m11s for all strains, no network I/O).

- [ ] **Step 3: Write the live-graph gate**

```python
# tests/kg_validity/test_controlled_vocabularies.py
"""Live-graph vocabulary contract checks.

The CSV gate is structurally blind to post-import-computed properties
(transport_substrate_resolution, tcdb_evidence_score_max, source_agreement,
pfam_support, go_support), so this suite is not redundant with it.
"""

import pytest

from multiomics_kg.utils.controlled_vocab import load_vocabularies, vocabularies_hash

pytestmark = pytest.mark.kg


@pytest.fixture(scope="module")
def declared():
    return load_vocabularies()


def test_every_vocabulary_has_a_node(run_query, declared):
    rows = run_query("MATCH (v:ControlledVocabulary) RETURN v.id AS id")
    assert {r["id"] for r in rows} == set(declared)


def test_observed_values_are_declared(run_query, declared):
    problems = []
    for e in declared.values():
        if not e.closed:
            continue
        if e.applies_to_kind == "node":
            q = (f"MATCH (n:{e.applies_to}) WHERE n.`{e.property}` IS NOT NULL "
                 f"UNWIND CASE WHEN n.`{e.property}` IS NULL THEN [] "
                 f"WHEN valueType(n.`{e.property}`) STARTS WITH 'LIST' "
                 f"THEN n.`{e.property}` ELSE [n.`{e.property}`] END AS v "
                 f"RETURN DISTINCT v AS v")
        else:
            q = (f"MATCH ()-[r:{e.applies_to}]->() WHERE r.`{e.property}` IS NOT NULL "
                 f"UNWIND CASE WHEN valueType(r.`{e.property}`) STARTS WITH 'LIST' "
                 f"THEN r.`{e.property}` ELSE [r.`{e.property}`] END AS v "
                 f"RETURN DISTINCT v AS v")
        seen = {r["v"] for r in run_query(q)}
        undeclared = seen - set(e.values)
        if undeclared:
            problems.append(f"{e.id}: undeclared {sorted(undeclared)}")
        if not e.expected_empty and seen:
            unseen = set(e.values) - seen
            if unseen:
                problems.append(f"{e.id}: declared but absent {sorted(unseen)}")
        if e.expected_empty:
            assert not seen, f"{e.id} is expected_empty but the graph has {seen}"
    assert not problems, "\n".join(problems)


def test_every_sources_value_joins_a_data_source(run_query):
    """R2: sources values are DataSource ids."""
    rows = run_query("""
        MATCH ()-[r]->() WHERE r.sources IS NOT NULL
        UNWIND r.sources AS s
        WITH DISTINCT s
        WHERE NOT EXISTS { MATCH (d:DataSource) WHERE d.id = s }
        RETURN collect(s) AS orphans
    """)
    assert rows[0]["orphans"] == [], (
        f"sources values with no DataSource node: {rows[0]['orphans']}")


def test_no_native_boolean_properties_remain(run_query):
    """R5. Checks the properties this change converted."""
    for label, prop in [("TcdbFamily", "is_promiscuous"),
                        ("InterproEntry", "is_promiscuous")]:
        rows = run_query(
            f"MATCH (n:{label}) WHERE n.`{prop}` IS NOT NULL RETURN count(n) AS n")
        assert rows[0]["n"] == 0, f"{label}.{prop} should have been deleted"


def test_layer_a_router_edges_carry_no_properties(run_query):
    for rel in ("Interpro_entry_related_to_ec_number",
                "Interpro_entry_related_to_cazy_family"):
        rows = run_query(
            f"MATCH ()-[r:{rel}]->() WHERE size(keys(r)) > 0 RETURN count(r) AS n")
        assert rows[0]["n"] == 0, f"{rel} edges still carry properties"


def test_schema_info_carries_the_vocabulary_hash(run_query, declared):
    rows = run_query(
        "MATCH (s:Schema_info {id:'schema_info'}) "
        "RETURN s.controlled_vocabularies_hash AS h")
    assert rows[0]["h"] == vocabularies_hash(list(declared.values())), (
        "Schema_info hash does not match the shipped config — the build and the "
        "checkout disagree about the vocabulary set.")
```

- [ ] **Step 4: Run the fast gates**

Run: `pytest -m "not slow and not kg" -q`
Expected: PASS. The `kg` suite cannot pass until Task 9 lands the hash and the graph is rebuilt (Task 12).

- [ ] **Step 5: Commit**

```bash
git add tests/test_create_knowledge_graph.py tests/kg_validity/test_controlled_vocabularies.py
git commit -m "test(vocab): four-gate vocabulary drift enforcement

Parameterizes the build test over --test and full modes, which also moves the
ERROR-line smoke check into the default run. The live-graph gate is not
redundant: the CSV scan is structurally blind to post-import-computed
properties."
```

---

## Task 9: Stamp `controlled_vocabularies_hash` on `Schema_info`

**Files:**
- Modify: `create_knowledge_graph.py` (write the hash next to the CSVs)
- Modify: `scripts/post-import.sh` Group 4 (~lines 1633-1692) and the matching block in `scripts/post-import.cypher`
- Modify: `tests/kg_validity/test_schema_info.py` (add the property to `STRING_PROPS`)

**Interfaces:**
- Consumes: `vocabularies_hash`, `VOCAB.entries()` from Task 1
- Produces: `Schema_info.controlled_vocabularies_hash: str`

**Why a file rather than computing it in the container:** the post-process container runs the neo4j image and drives `cypher-shell`; it has no Python. The build stage already writes to the shared output volume, so it writes the hash there and `post-import.sh` `cat`s it into a `-P` param — the same shape as `KG_RELEASE_HIGHLIGHTS`.

- [ ] **Step 1: Write the failing test**

Add to `tests/kg_validity/test_schema_info.py`:

```python
def test_controlled_vocabularies_hash_present_and_matches(schema_info):
    from multiomics_kg.utils.controlled_vocab import VOCAB, vocabularies_hash
    key = "controlled_vocabularies_hash"
    assert key in schema_info, (
        "Schema_info is missing the vocabulary hash. Post-import Group 4 did "
        "not read output/controlled_vocabularies.sha256 — check that the build "
        "stage wrote it.")
    assert schema_info[key] == vocabularies_hash(VOCAB.entries())
```

Also append `"controlled_vocabularies_hash"` to `STRING_PROPS`.

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/kg_validity/test_schema_info.py -k vocabular -v`
Expected: FAIL (or skip if Neo4j is down — start the graph first, or defer this run to Task 12).

- [ ] **Step 3: Write the hash from the build**

In `create_knowledge_graph.py`, immediately after `bc.write_nodes(controlled_vocab_adapter.get_nodes())`:

```python
    # Post-import stamps this onto Schema_info so kg_release_info can detect
    # vocabulary drift. Written to the output dir because the post-process
    # container runs the neo4j image and has no Python.
    from multiomics_kg.utils.controlled_vocab import vocabularies_hash
    hash_path = Path(output_dir) / "controlled_vocabularies.sha256"
    hash_path.write_text(vocabularies_hash(controlled_vocab_adapter._entries))
    logger.info("Wrote vocabulary hash to %s", hash_path)
```

Expose the entries properly rather than reaching into `_entries`: add to the adapter

```python
    def entries(self) -> list:
        if not self._entries:
            self.download_data()
        return list(self._entries)
```

and call `vocabularies_hash(controlled_vocab_adapter.entries())`.

- [ ] **Step 4: Stamp it in post-import Group 4**

In `scripts/post-import.sh`, before the `cypher-shell` invocation:

```bash
# Vocabulary hash written by the build stage. Empty when absent (a legacy
# build), which coalesces to '' and leaves the property null.
KG_VOCAB_HASH=""
if [ -r /output/controlled_vocabularies.sha256 ]; then
  KG_VOCAB_HASH=$(cat /output/controlled_vocabularies.sha256)
fi
```

Add the param `-P "vocab_hash => '${KG_VOCAB_HASH}'"` and, inside the heredoc, extend the first `SET`:

```cypher
    s.controlled_vocabularies_hash = CASE WHEN coalesce($vocab_hash, '') = ''
                                          THEN null ELSE $vocab_hash END,
```

Mirror the `SET` line into `scripts/post-import.cypher`'s Group 4 block.

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "feat(vocab): stamp controlled_vocabularies_hash on Schema_info

Lets kg_release_info detect vocabulary drift rather than the explorer
discovering it through a wrong answer. Written to the output volume by the
build because the post-process container has no Python."
```

---

## Task 10: Seed the released vocabularies (spec §5.3 second group)

**Files:**
- Modify: `config/controlled_vocabularies.yaml`
- Test: `tests/test_controlled_vocab.py`

**Interfaces:** no new code — data only.

**This task is a harvest, not a design.** Values are read out of the code and the live graph, not invented. Nothing here renames anything: these vocabularies are released and MCP-read.

- [ ] **Step 1: Harvest the value sets**

For each vocabulary below, get the authoritative list and record where it came from in the YAML `description`:

```bash
# From code/config:
grep -n "KNOWN_METRIC_TYPES" -A40 multiomics_kg/utils/paperconfig_utils.py   # metric_type
grep -n "omics_type\|compartment" config/schema_config.yaml | head -20
grep -n "BRITE_TREES\|tree_code" multiomics_kg/utils/brite_utils.py          # brite_tree

# From the live graph (start Docker first):
uv run python - <<'PY'
from neo4j import GraphDatabase
q = {
 "OrganismTaxon.organism_type": "MATCH (n:OrganismTaxon) RETURN DISTINCT n.organism_type AS v",
 "Experiment.omics_type":       "MATCH (n:Experiment) RETURN DISTINCT n.omics_type AS v",
 "Experiment.compartment":      "MATCH (n:Experiment) RETURN DISTINCT n.compartment AS v",
 "Experiment.table_scope":      "MATCH (n:Experiment) RETURN DISTINCT n.table_scope AS v",
 "Changes_expression_of.expression_status":
     "MATCH ()-[r:Changes_expression_of]->() RETURN DISTINCT r.expression_status AS v",
 "Metabolite.evidence_sources":
     "MATCH (n:Metabolite) UNWIND n.evidence_sources AS v RETURN DISTINCT v",
 "Gene.annotation_state":       "MATCH (n:Gene) RETURN DISTINCT n.annotation_state AS v",
}
d = GraphDatabase.driver("bolt://localhost:7687")
with d.session() as s:
    for k, cy in q.items():
        vals = sorted(x["v"] for x in s.run(cy) if x["v"] is not None)
        print(f"{k}: {vals}")
PY
```

Cross-check every list against the corresponding bullet in `CLAUDE.md` "Key graph facts". **A mismatch is a finding, not a nuisance** — record it and reconcile before declaring, because one of the two is wrong.

- [ ] **Step 2: Declare them**

Add entries following the Task 1 shape. Closed: `omics_type`, `compartment`, `table_scope`, `value_kind`, `metric_type`, `brite_tree`, `evidence_sources`, `expression_status`, `metric_bucket`, `detection_status`, `annotation_state`, `organism_type`, `treatment_type`, `background_factors`, and `level_kind` **one entry per label** (`TcdbFamily`, `CazyFamily`, `KeggTerm`, `BriteCategory`, `InterproEntry`-as-`expected_empty` from Task 1).

Open (`closed: false`, `values: []`): `Gene.gene_category` and `Experiment.growth_phases` — both are data-driven, so `list_filter_values` must enumerate them from the graph. Say so in the description.

Grandfathered `"true"`/`"false"` properties get `value_type: bool_string` and a description noting they predate R5 and are backlogged:

```yaml
DerivedMetric.rankable:
  applies_to: DerivedMetric
  applies_to_kind: node
  property: rankable
  value_type: bool_string
  closed: true
  values: ["true", "false"]
  description: >
    Whether this metric supports ranking. Stringified boolean predating house
    rule R5 — grandfathered because it is MCP-read; conversion to a meaningful
    pair is tracked in plans/backlog.md.
```

- [ ] **Step 3: Add the coverage test**

```python
def test_every_list_filter_values_filter_is_declared():
    """The 8 filters the MCP serves must all be declared, closed or open."""
    from multiomics_kg.utils.controlled_vocab import load_vocabularies
    declared = load_vocabularies()
    required = [
        "Gene.gene_category", "BriteCategory.tree", "Experiment.growth_phases",
        "DerivedMetric.metric_type", "DerivedMetric.value_kind",
        "Experiment.compartment", "Experiment.omics_type",
        "Metabolite.evidence_sources",
    ]
    missing = [k for k in required if k not in declared]
    assert not missing, f"list_filter_values filters not declared: {missing}"
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/test_controlled_vocab.py -v && pytest -m "not slow and not kg" -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add config/controlled_vocabularies.yaml tests/test_controlled_vocab.py
git commit -m "feat(vocab): seed the released vocabularies

Harvested from code and the live graph, cross-checked against CLAUDE.md. No
renames — these are released and MCP-read. gene_category and growth_phases are
declared open so list_filter_values enumerates them at runtime; the seven
stringified booleans are declared bool_string and grandfathered."
```

---

## Task 11: Documentation

**Files:**
- Modify: `CLAUDE.md`, `CHANGELOG.md`
- Modify: `docs/kg-changes/{interproscan-extension,interpro-two-layer,tcdb-two-source-upgrade}.md`
- Create: `docs/kg-changes/vocabulary-contract.md` (the explorer-facing what-changed doc)

- [ ] **Step 1: `CLAUDE.md`**

Add a "Controlled vocabularies" subsection under **Architecture → Adapter Pattern**, right after the string-sanitization block, stating the five house rules in one line each and pointing at `config/controlled_vocabularies.yaml` as the source of truth. Then update the affected "Key graph facts" bullets: the `TcdbFamily` bullet (drop `is_promiscuous`, note the documented threshold), the `InterproEntry` bullet (same, plus lowercase `interpro_type` values and lowercase `libraries`), the InterPro-two-layer bullet (`evidence_score` is a float; `sources` values are `interproscan`), the `Gene_has_tcdb_family` bullet (`tcdb_diamond`; `evidence_score`; the three renamed support properties), and the substrate-rollup bullet (`most_specific` / `inherited`). Add `ControlledVocabulary` to the "Actual Neo4j labels" node list.

- [ ] **Step 2: The five chemistry definitions (KG-IPT-004)**

Paste §7.1 of the spec verbatim into `docs/kg-changes/tcdb-two-source-upgrade.md`, beside the existing property tables. While in that file, fix the two stale items the explorer flagged: the `tcdb_evidence_score` distribution (now re-scaled to floats by R4 anyway) and the `coalesce(..., -1)` guidance, which becomes `-1.0`.

- [ ] **Step 3: `docs/kg-changes/vocabulary-contract.md`**

The consumer-facing doc. Sections: what `ControlledVocabulary` is and how to query it; the rename table (copy spec §4); the four deletions and their derivation recipes; the five house rules; and the migration notes for the explorer's entry-criteria queries (spec §0 rev 4/5 and §9.1).

- [ ] **Step 4: `CHANGELOG.md` under `## [Unreleased]`**

`### Breaking`: one bullet covering `TcdbFamily.is_promiscuous` deleted (with the derivation predicate), `evidence_score` integer → float on every annotation edge, `tcdb_evidence_score` → `evidence_score`, `tcdb_best_evidence_score` → `tcdb_evidence_score_max`, `substrate_depth` value rename, `interpro_type`/`libraries` lowercasing, and the `sources` value renames.
`### Added`: the `ControlledVocabulary` contract and the R5 boolean rule.
`### Highlights`: one plain-prose bullet — consumers can now ask the graph what values a filter accepts instead of guessing.

Also apply the KG-IPT-007 fix: in the existing biller 2016 `### Data` bullet, `control_condition` → `control`.

- [ ] **Step 5: Commit**

```bash
git add -u && git add docs/kg-changes/vocabulary-contract.md
git commit -m "docs: vocabulary contract, house rules, chemistry definitions

Adds the consumer-facing what-changed doc, the five verbatim chemistry-count
definitions (KG-IPT-004), and the CHANGELOG entries. Fixes KG-IPT-007:
CHANGELOG prose now uses the exact graph property name (control)."
```

---

## Task 12: Full validation gate

**Files:** none modified — this is the acceptance run.

**Do not skip and do not reorder.** The baseline must be captured *before* the rebuild, against the currently deployed graph.

- [ ] **Step 1: Capture the baseline (before rebuilding)**

```bash
bash scripts/post-import-validate.sh > /tmp/vocab-baseline.txt
uv run python - <<'PY' > /tmp/vocab-structural-before.txt
from neo4j import GraphDatabase
qs = {
 "InterproEntry": "MATCH (n:InterproEntry) RETURN count(n) AS n",
 "Gene_has_interpro_entry": "MATCH ()-[r:Gene_has_interpro_entry]->() RETURN count(r) AS n",
 "TcdbFamily": "MATCH (n:TcdbFamily) RETURN count(n) AS n",
 "Gene_has_tcdb_family": "MATCH ()-[r:Gene_has_tcdb_family]->() RETURN count(r) AS n",
 "Tcdb_family_transports_metabolite":
   "MATCH ()-[r:Tcdb_family_transports_metabolite]->() RETURN count(r) AS n",
 "GO_inferred": "MATCH ()-[r]->() WHERE type(r) STARTS WITH 'Gene_' AND "
                "r.evidence IN ['family_inferred','domain_inferred'] RETURN count(r) AS n",
}
d = GraphDatabase.driver("bolt://localhost:7687")
with d.session() as s:
    for k, q in qs.items():
        print(k, s.run(q).single()["n"])
PY
```

Expected (the explorer's step-1 acceptance figures): `InterproEntry` 12,999 · `Gene_has_interpro_entry` 397,342 · `TcdbFamily` 1,515 · `Gene_has_tcdb_family` 53,763 · `Tcdb_family_transports_metabolite` 11,263 · GO inferred 45,226.

- [ ] **Step 2: Unit tests**

Run: `pytest -m "not slow and not kg" -v`
Expected: all pass.

- [ ] **Step 3: Rebuild the graph**

```bash
docker compose down
docker compose up -d --build
docker compose logs -f post-process   # watch for the [timing] lines and Group 4
```

- [ ] **Step 4: Re-run the validate dump and diff**

```bash
bash scripts/post-import-validate.sh > /tmp/vocab-after.txt
diff /tmp/vocab-baseline.txt /tmp/vocab-after.txt
```

Expected: differences **only** in the renamed properties and the re-scaled score columns. Any changed count is a regression — stop and investigate rather than accepting it. Re-run the structural query and confirm the six counts are unchanged.

- [ ] **Step 5: Graph tests**

Run: `pytest -m kg -v`
Expected: all pass, including the new `tests/kg_validity/test_controlled_vocabularies.py` and the `Schema_info` hash check. The two known-failing orphan-protein tests are pre-existing (`CLAUDE.md` → Known Issues) and are not caused by this change.

- [ ] **Step 6: Edge snapshot + regression snapshot**

```bash
# expression edges must be untouched by this change
```
Run the `/omics-edge-snapshot` skill and confirm zero per-paper deltas. Then regenerate the node/edge snapshot, since renamed properties legitimately change it:
```bash
uv run python tests/kg_validity/generate_snapshot.py
git add tests/kg_validity/snapshot_data.json
```

- [ ] **Step 7: Run the explorer's entry criteria**

Paste and run the queries from spec §0 (rev 4/5) and §9.1. Confirm: `substrate_depth` is `most_specific` 4,381 / `inherited` 6,882 with no `deepest`/`ancestor`; `interpro_type` has 8 lowercase values; the derived predicates return 13 and 22; the derived `ambiguous` returns 4,865; `ControlledVocabulary` node count equals the declared count; `Schema_info.controlled_vocabularies_hash` is non-null.

- [ ] **Step 8: Commit**

```bash
git add tests/kg_validity/snapshot_data.json
git commit -m "test: regenerate snapshot after vocabulary renames

Full validation: post-import-validate diff shows only renamed properties and
re-scaled score columns; all six structural counts unchanged; pytest -m kg
green; omics-edge-snapshot shows zero per-paper deltas."
```

- [ ] **Step 9: Tell the explorer step 1 is done**

Per the agreed sequencing (spec §9.6), the explorer's W1 + W4 work is gated on this deploy. Send them: the redeployed build, `docs/kg-changes/vocabulary-contract.md`, and the confirmation that their entry criteria pass — noting the three queries that changed shape.

---

## Self-Review

**Spec coverage.** R1 → Task 3. R2 → Task 4. R3 → Task 7. R4 → Task 5. R5 → Tasks 6 + 7. §5 contract → Tasks 1, 2, 9, 10. §6 enforcement → Task 8. §7 answers → Task 11 (004, 007) and the reply doc; 002/005/006 need no code. §9 follow-ups: 009+013 → Task 7, 010 → Tasks 1+3, 011 → Task 1, 012 → Task 5. §8 validation → Task 12. §10 open items: 1 → Task 8 Step 2, 2 → Task 10, 3 and 4 and 5 stay backlogged.

**Ordering constraint.** Task 6's `test_no_schema_property_is_declared_bool` cannot pass until Task 7 deletes the two `ambiguous: bool` entries. Run 7 immediately after 6, or `xfail` that one assertion in 6 and flip it in 7. Called out in Task 6 Step 6.

**Type consistency.** `VOCAB.check(applies_to, property, value) -> str` is used identically in Tasks 3, 6. `VocabEntry.id` is `<applies_to>.<property>` and is the node id in Task 2, the YAML key in Task 1, and the join key in Task 8. `vocabularies_hash(entries: list[VocabEntry]) -> str` is called in Tasks 1, 8, 9 with `VOCAB.entries()` / `adapter.entries()`, both of which return `list[VocabEntry]`.

**Known drafting artifact.** Task 1 Step 4 shows `vocabularies_hash` with an `if False else` construction and then gives the correct body immediately below. Use the second version.
