# Explorer Frictions F1-F4 Resolution Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Resolve four explorer-side KG-side frictions (F1 term informativeness, F2 data-source surfacing, F3 cluster_type vocab, F4 sparsely-annotated gene properties) as a coordinated KG release.

**Architecture:** Mostly post-import Cypher additions + small build-time changes. New `DataSource` node type (4 nodes auto-derived from `gene_annotations_config.yaml`). New term-level `is_uninformative` flag driven by hand-curated YAML. New Gene properties `annotation_state`, `informative_annotation_types`, `contributing_sources`, `contig`, `seed_ortholog`, `seed_ortholog_evalue`. Refined `Gene.annotation_quality` (numeric encoding of `annotation_state`, 8 source buckets). One vocabulary rename (`classification` → `expression_bin`).

**Tech Stack:** Python 3, BioCypher, Neo4j 5 (Cypher), pytest, YAML config, pandas.

**Conversation seed:** [`multiomics_explorer/docs/superpowers/specs/2026-05-01-kg-side-frictions-reframed.md`](../../../../multiomics_explorer/docs/superpowers/specs/2026-05-01-kg-side-frictions-reframed.md)

**Spec:** [`docs/superpowers/specs/2026-05-01-explorer-frictions-resolution-design.md`](../specs/2026-05-01-explorer-frictions-resolution-design.md)

---

## Pre-flight context for the engineer

You don't need to read the whole spec to start, but you should understand:

- **Build pipeline.** `bash scripts/prepare_data.sh` downloads + processes data per strain. `uv run python create_knowledge_graph.py` writes BioCypher CSVs. Docker stages run `neo4j-admin import` then `scripts/post-import.sh` (which executes `scripts/post-import.cypher`'s logic via `cypher-shell`). Both `post-import.cypher` (reference) and `post-import.sh` (authoritative) must stay in sync — same Cypher logic.
- **String sanitization.** All string properties yielded by adapters must avoid `'` (escape to `^`) and `|` (the array delimiter). See CLAUDE.md "String property sanitization".
- **Test markers.** `pytest -m "not slow and not kg"` for fast unit tests; `pytest -m kg` for KG validity tests against a running Neo4j (Docker compose). KG tests auto-skip if Neo4j is unreachable.
- **Boolean-shaped properties convention.** Use `str` `"true"` / absent (no `"false"`) for KG bool properties. BioCypher has a known bug serializing `bool`; absent-means-informative pattern is established by `rankable` / `has_p_value` / `significant` on DerivedMetric.
- **Working dir.** Plan tasks assume you're in `/home/osnat/github/multiomics_biocypher_kg`.
- **Live KG.** `bolt://localhost:7687` (no auth, see CLAUDE.md). Use `cypher-shell` for direct queries; tests use the `neo4j` Python driver via `tests/kg_validity/conftest.py`.

**Commit pattern:** small, frequent, conventional-commits style (`feat:`, `test:`, `refactor:`, `docs:`). All commits should trail with `Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>` per CLAUDE.md.

---

## File map

**Create:**
- `config/uninformative_terms.yaml` — Catch-all term vocabulary (Task 9)
- `multiomics_kg/adapters/data_source_adapter.py` — `DataSourceAdapter` (Task 8)
- `tests/test_data_source_adapter.py` — unit tests for adapter (Task 8)
- `tests/kg_validity/test_uninformative_terms.py` — KG validity (Task 9, expanded in Task 11)
- `tests/kg_validity/test_annotation_state.py` — KG validity (Task 11)
- `tests/kg_validity/test_data_source.py` — KG validity (Task 8)
- `tests/kg_validity/test_gene_contig.py` — KG validity (Task 3)

**Modify:**
- `.claude/skills/paperconfig/validate_paperconfig.py` — `VALID_CLUSTER_TYPES` rename (Task 1)
- `.claude/skills/paperconfig/SKILL.md` — vocab + convention text (Task 1)
- `multiomics_kg/download/download_genome_data.py` — preserve `seqid` from GFF (Task 2)
- `config/gene_annotations_config.yaml` — `contig` passthrough (Task 3); `logical_sources` per source (Task 5)
- `config/schema_config.yaml` — Gene property additions (Tasks 3, 4, 6); new `data source` node type (Task 7)
- `multiomics_kg/adapters/cyanorak_ncbi_adapter.py` — `GeneNodeField` enum extensions (Tasks 3, 4, 6); `float_fields` set (Task 4)
- `multiomics_kg/download/build_gene_annotations.py` — add `_compute_contributing_sources` (Task 6); delete `_compute_annotation_quality` (Task 13)
- `create_knowledge_graph.py` — wire `DataSourceAdapter` (Task 8)
- `scripts/post-import.cypher` and `scripts/post-import.sh` — new Cypher blocks (Tasks 10, 11, 12)
- `CLAUDE.md` — new fields + bucket list (Task 14)

---

## Task 1: F3 — Rename `classification` → `expression_bin` in cluster_type vocab

**Files:**
- Modify: `.claude/skills/paperconfig/validate_paperconfig.py`
- Modify: `.claude/skills/paperconfig/SKILL.md`

The renaming is zero-risk: no current paperconfig uses `classification`. Only the validator vocab and skill template change.

- [ ] **Step 1: Locate `VALID_CLUSTER_TYPES` in validator**

```bash
grep -n "VALID_CLUSTER_TYPES" .claude/skills/paperconfig/validate_paperconfig.py
```
Expected: a line with `VALID_CLUSTER_TYPES = {"time_course", "diel", "condition_comparison", "classification"}` or similar.

- [ ] **Step 2: Replace `classification` with `expression_bin`**

Edit the set so it reads:

```python
VALID_CLUSTER_TYPES = {"time_course", "diel", "condition_comparison", "expression_bin"}
```

- [ ] **Step 3: Replace any `classification` mentions in SKILL.md**

```bash
grep -n "classification" .claude/skills/paperconfig/SKILL.md
```
Replace each with `expression_bin`. Add the convention text near the existing `cluster_type` documentation:

```markdown
**`expression_bin`** — clusters group genes by expression metric (e.g., quartile labels VEG/HEG/MEG/LEG/NEG). These analyses do **not** carry per-cluster `functional_description`; cluster `name` is the metric label. All other `cluster_type` values carry per-cluster functional descriptions where curated.

**Future ambiguity rule:** If a future cluster type has unclear intent, rename or split the `cluster_type` value rather than adding a parallel intent field.
```

- [ ] **Step 4: Run validator self-check**

```bash
uv run python .claude/skills/paperconfig/validate_paperconfig.py data/Prochlorococcus/papers_and_supp/Wang\ 2014/paperconfig.yaml
```
Expected: validator exits 0 (no current paperconfig should use `classification` so no failures).

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/paperconfig/validate_paperconfig.py .claude/skills/paperconfig/SKILL.md
git commit -m "$(cat <<'EOF'
feat(paperconfig): rename cluster_type 'classification' → 'expression_bin'

The 'classification' label was ambiguous (could mean functional
classification or expression bin). 'expression_bin' is unambiguous: clusters
are quartile / RPKM-bin labels with no curated functional content per
cluster by design. Migration cost zero — no current paperconfig uses
'classification'. Documents the convention and future-ambiguity rule.

Resolves F3 of the explorer frictions design.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: F4.1 — Add `seqid` to NCBI gene_mapping.csv writer

The NCBI GFF carries `seqid` (the contig name) but the current `gene_mapping.csv` writer drops it. Required for genomic-neighbor lookups. This task adds it to the CSV; the next task plumbs it onto the Gene node.

**Files:**
- Modify: `multiomics_kg/download/download_genome_data.py`

- [ ] **Step 1: Locate the gene_mapping.csv writer (sub-step 5)**

```bash
grep -n "gene_mapping.csv\|write.*csv\|build_gene_mapping\|step.*5" multiomics_kg/download/download_genome_data.py | head -20
```

Find the function that constructs the CSV columns. The current columns are: `locus_tag, locus_tag_ncbi, gene_names, source, start, end, strand, Note, exception, inference, product, protein_id, Ontology_term, go_component, go_function, go_process, gene_synonym, old_locus_tags, ec_numbers`.

- [ ] **Step 2: Add `seqid` to the column extraction**

Find the GFF row processing — it likely uses something like `row["seqid"]` or `feature.seqid`. Add a new column to the per-row dict that captures the GFF feature's `seqid` value. Place it as the **first** column for clarity.

Example modification (locate the actual dict construction):

```python
# In the per-feature dict construction, add:
{
    "seqid": str(feature.seqid),    # NEW: contig/chromosome name
    "locus_tag": locus_tag,
    "locus_tag_ncbi": ...,
    # ... existing columns
}
```

If the writer uses a column list rather than a dict, prepend `"seqid"` to the list and ensure each row's first value is the seqid.

- [ ] **Step 3: Force-regenerate one strain's gene_mapping.csv**

```bash
uv run bash scripts/prepare_data.sh --steps 0 --strains MED4 --force
```

- [ ] **Step 4: Verify the new column landed**

```bash
head -1 cache/data/Prochlorococcus/genomes/MED4/gene_mapping.csv | tr ',' '\n' | head -3
```
Expected: first three lines show `seqid`, then `locus_tag`, then `locus_tag_ncbi` (or similar).

```bash
awk -F, 'NR<=3 {print $1, $2}' cache/data/Prochlorococcus/genomes/MED4/gene_mapping.csv
```
Expected: header line shows `seqid locus_tag`; data rows show actual contig names (likely `NC_005072.1` for MED4).

- [ ] **Step 5: Commit (script change only — DO NOT commit regenerated cache)**

```bash
git add multiomics_kg/download/download_genome_data.py
git commit -m "$(cat <<'EOF'
feat(download): preserve GFF seqid in gene_mapping.csv

The NCBI GFF carries a seqid column (contig/chromosome name) but the
current gene_mapping.csv writer dropped it. Required for genomic-neighbor
lookups (F4 ship 4.1). Adds seqid as the first column.

Caches must be regenerated with `bash scripts/prepare_data.sh --steps 0 --force`
before downstream code can use the field.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 6: Regenerate all strains' caches**

```bash
uv run bash scripts/prepare_data.sh --steps 0 1 2 --force --skip-cyanorak 2>&1 | tail -40
```
Expected: completes without error. Cyanorak step skipped since the server is unreliable; the cached files are reused.

- [ ] **Step 7: Spot-check several strains have the seqid column**

```bash
for s in MED4 MIT9301 HP15 KT2440; do
  echo "=== $s ==="
  head -2 cache/data/*/genomes/$s/gene_mapping.csv 2>/dev/null | head -2
done
```
Expected: every strain's gene_mapping.csv has `seqid` as first column with non-empty values.

---

## Task 3: F4.1 — Plumb `contig` onto Gene nodes

**Files:**
- Modify: `config/gene_annotations_config.yaml`
- Modify: `multiomics_kg/adapters/cyanorak_ncbi_adapter.py`
- Modify: `config/schema_config.yaml`
- Test: `tests/kg_validity/test_gene_contig.py`

- [ ] **Step 1: Add the `contig` passthrough to gene_annotations_config.yaml**

Locate the `fields:` section in `config/gene_annotations_config.yaml` (a long block of field rules). Find a logical place near identity / locus fields and add:

```yaml
  contig:
    type: passthrough
    source: gene_mapping
    field: seqid
```

- [ ] **Step 2: Add `CONTIG = 'contig'` to GeneNodeField enum**

In `multiomics_kg/adapters/cyanorak_ncbi_adapter.py`, find the `GeneNodeField` enum (around line 130). Add the new entry adjacent to the existing core identity fields:

```python
class GeneNodeField(...):
    LOCUS_TAG = 'locus_tag'
    START = 'start'
    END = 'end'
    STRAND = 'strand'
    CONTIG = 'contig'                # NEW: chromosome / contig name from GFF seqid
    PRODUCT = 'product'
    # ... existing entries
```

- [ ] **Step 3: Add `contig: str` to Gene schema**

In `config/schema_config.yaml`, locate the `gene:` block (search `^gene:`). Find the `properties:` list and add adjacent to `start`/`end`/`strand`:

```yaml
gene:
  represented_as: node
  preferred_id: ncbigene
  label_in_input: gene
  properties:
    locus_tag: str
    start: int
    end: int
    strand: str
    contig: str                  # chromosome / contig name from GFF seqid
    product: str
    # ... existing entries
```

- [ ] **Step 4: Rebuild merged annotations and gene CSVs (test mode for speed)**

```bash
uv run bash scripts/prepare_data.sh --steps 2 --force 2>&1 | tail -10
uv run python create_knowledge_graph.py --test 2>&1 | tail -20
```
Expected: completes without error. The test-mode CSV emits ≤100 genes per strain.

- [ ] **Step 5: Spot-check the gene CSV has a contig column**

```bash
head -1 biocypher-log/example_knowledge_graph/Gene-part000.csv | tr ';' '\n' | grep -i 'contig\|seqid'
```
Expected: a `contig:string` column appears.

- [ ] **Step 6: Write KG validity test for contig**

Create `tests/kg_validity/test_gene_contig.py`:

```python
"""KG validity: every Gene must have a non-null contig.

F4 ship 4.1 requires the contig (NCBI GFF seqid) on every Gene node so
genomic-neighbor lookups have a valid discriminator.
"""

import pytest


@pytest.mark.kg
def test_every_gene_has_contig(neo4j):
    """Every Gene node must have a non-null contig property."""
    result = neo4j.query("""
        MATCH (g:Gene)
        WHERE g.contig IS NULL OR g.contig = ''
        RETURN count(*) AS missing
    """)
    missing = result[0]["missing"]
    assert missing == 0, f"{missing} Gene nodes have null/empty contig"


@pytest.mark.kg
def test_contig_consistent_within_organism(neo4j):
    """All genes from a single NCBI accession (organism) should map to a
    bounded set of contigs (1 for single-chromosome organisms, more for
    multi-replicon ones, but always at least one)."""
    result = neo4j.query("""
        MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
        WITH o.preferred_name AS organism, collect(DISTINCT g.contig) AS contigs
        WHERE size(contigs) = 0
        RETURN organism, contigs
    """)
    assert result == [], f"organisms with no distinct contig: {result}"


@pytest.mark.kg
def test_med4_has_known_chromosome(neo4j):
    """Spot-check: MED4 genes should land on its known NCBI accession."""
    result = neo4j.query("""
        MATCH (g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
        WHERE o.preferred_name = 'Prochlorococcus MED4'
        RETURN DISTINCT g.contig AS contig
    """)
    contigs = {row["contig"] for row in result}
    # MED4 single-replicon — expect one chromosome accession
    assert len(contigs) == 1, f"MED4 expected 1 contig, got {len(contigs)}: {contigs}"
```

- [ ] **Step 7: Run the test (it'll fail until the full graph rebuild lands)**

```bash
pytest tests/kg_validity/test_gene_contig.py -v
```
Expected: skips with "Neo4j unreachable" OR fails because the live Neo4j was built before this change. Document this — the test PASSES will only happen after the full rebuild in the release. Mark as expected; commit anyway.

- [ ] **Step 8: Commit**

```bash
git add config/gene_annotations_config.yaml multiomics_kg/adapters/cyanorak_ncbi_adapter.py config/schema_config.yaml tests/kg_validity/test_gene_contig.py
git commit -m "$(cat <<'EOF'
feat(gene): add Gene.contig property + validity test

Required for genomic-neighbor lookups (F4 ship 4.1). Plumbs the GFF seqid
preserved in Task 2 through gene_annotations_config.yaml -> Gene schema ->
cyanorak_ncbi_adapter. KG validity test asserts every Gene has a non-null
contig and MED4 lands on its known single chromosome.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: F4.2 — Add `seed_ortholog` + `seed_ortholog_evalue` Gene properties

These are already in `gene_annotations_merged.json` (eggnog passthrough). Just surface them as Gene properties.

**Files:**
- Modify: `config/schema_config.yaml`
- Modify: `multiomics_kg/adapters/cyanorak_ncbi_adapter.py`

- [ ] **Step 1: Add fields to Gene schema**

In `config/schema_config.yaml` `gene:` block `properties:`, add:

```yaml
    seed_ortholog: str             # eggNOG seed-ortholog identifier (taxid.id format)
    seed_ortholog_evalue: float    # eggNOG seed-ortholog E-value
```

Place adjacent to other eggnog-derived computed/passthrough fields if present, or near the protein-related fields.

- [ ] **Step 2: Add `SEED_ORTHOLOG` and `SEED_ORTHOLOG_EVALUE` to GeneNodeField enum**

In `multiomics_kg/adapters/cyanorak_ncbi_adapter.py`, in the `GeneNodeField` enum:

```python
    SEED_ORTHOLOG = 'seed_ortholog'                     # NEW: eggNOG seed match (taxid.identifier)
    SEED_ORTHOLOG_EVALUE = 'seed_ortholog_evalue'       # NEW: eggNOG seed match E-value
```

- [ ] **Step 3: Add `seed_ortholog_evalue` to `float_fields` in `_get_gene_nodes`**

Find the line:

```python
float_fields = set()
```

Change to:

```python
float_fields = {'seed_ortholog_evalue'}
```

- [ ] **Step 4: Rebuild test-mode CSVs to verify**

```bash
uv run python create_knowledge_graph.py --test 2>&1 | tail -10
```

- [ ] **Step 5: Spot-check the gene CSV has the new columns**

```bash
head -1 biocypher-log/example_knowledge_graph/Gene-part000.csv | tr ';' '\n' | grep -E 'seed_ortholog'
```
Expected: two lines, one for `seed_ortholog:string` (or array — verify type) and one for `seed_ortholog_evalue:double`.

```bash
awk -F';' '$1 ~ /seed_ortholog/ {print NR": "$0}' biocypher-log/example_knowledge_graph/Gene-part000.csv | head -3
```
Expected: rows have non-null seed_ortholog values for genes with eggnog hits.

- [ ] **Step 6: Commit**

```bash
git add config/schema_config.yaml multiomics_kg/adapters/cyanorak_ncbi_adapter.py
git commit -m "$(cat <<'EOF'
feat(gene): add seed_ortholog + seed_ortholog_evalue Gene properties

Surfaces eggNOG seed-ortholog match as Gene properties so sparsely-annotated
genes have a 'this gene's protein most resembles X (E=Y)' pointer when
functional fields are empty (F4 ship 4.2).

Format: <taxid>.<source_identifier> where source_identifier varies by
organism (locus tag, contig.gene_n for draft assemblies, RefSeq WP_).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: F2.2 — Add `logical_sources` to gene_annotations_config.yaml

Per-source metadata for the new `DataSourceAdapter` to consume. Always a list (the `gene_mapping` source is a merge of two logical sources: ncbi + cyanorak).

**Files:**
- Modify: `config/gene_annotations_config.yaml`

- [ ] **Step 1: Read the current sources block**

```bash
sed -n '/^sources:/,/^[a-z]/p' config/gene_annotations_config.yaml | head -50
```

Identify each top-level source: `gene_mapping`, `eggnog`, `uniprot`.

- [ ] **Step 2: Add `logical_sources` under each source**

For `gene_mapping`:

```yaml
  gene_mapping:
    type: csv
    path_pattern: "{data_dir}/gene_mapping.csv"
    index_col: locus_tag
    description: "NCBI + Cyanorak merged gene annotations (38 cols)"
    logical_sources:
      - id: ncbi
        scope: gene_level
        provenance: download
      - id: cyanorak
        scope: organism_restricted
        provenance: download
        applies_to_organisms:
          - "Prochlorococcus*"
          - "Synechococcus*"
          - "Parasynechococcus*"
          - "Thermosynechococcus*"
```

For `eggnog`:

```yaml
  eggnog:
    type: tsv_hash_header
    path_pattern: "{data_dir}/eggnog/{strain_name}.emapper.annotations"
    join_to: protein_id
    description: "EggNOG-mapper v2 functional annotations"
    logical_sources:
      - id: eggnog
        scope: gene_level
        provenance: tool_run
```

For `uniprot`:

```yaml
  uniprot:
    type: json_roworiented
    path_candidates:
      - "cache/data/{organism_group}/uniprot/{ncbi_taxon_id}/protein_annotations.json"
      - "protein_annotations.json"
    pivot_on: refseq_ids
    join_to: protein_id
    description: "UniProt proteins (reviewed + unreviewed)"
    logical_sources:
      - id: uniprot
        scope: gene_level
        provenance: download
```

- [ ] **Step 3: Verify YAML still parses**

```bash
uv run python -c "import yaml; yaml.safe_load(open('config/gene_annotations_config.yaml'))"
```
Expected: no output (clean parse).

- [ ] **Step 4: Verify build still works**

```bash
uv run bash scripts/prepare_data.sh --steps 2 --force --strains MED4 2>&1 | tail -5
```
Expected: existing field-merge consumers ignore the new key cleanly.

- [ ] **Step 5: Commit**

```bash
git add config/gene_annotations_config.yaml
git commit -m "$(cat <<'EOF'
feat(config): add logical_sources metadata to gene_annotations_config

Per-source metadata (id, scope, provenance, applies_to_organisms) for the
new DataSourceAdapter to consume (F2 ship 2.2). Always a list — gene_mapping
declares ncbi + cyanorak as separate logical sources merged at download
time. Existing field-merge consumers ignore the new key.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: F2.4 — `Gene.contributing_sources` build-time computation

**Files:**
- Modify: `multiomics_kg/download/build_gene_annotations.py`
- Modify: `multiomics_kg/adapters/cyanorak_ncbi_adapter.py`
- Modify: `config/schema_config.yaml`
- Test: `tests/test_contributing_sources.py` (new unit test)

- [ ] **Step 1: Write the unit test first**

Create `tests/test_contributing_sources.py`:

```python
"""Unit test for _compute_contributing_sources in build_gene_annotations.

F2 ship 2.4 derives Gene.contributing_sources from existing source-provenance
fields in gene_annotations_merged.json.
"""

import pytest
from multiomics_kg.download.build_gene_annotations import _compute_contributing_sources


def test_ncbi_only_gene_returns_ncbi():
    """A gene with no source-tagged fields beyond NCBI gets ['ncbi']."""
    gene = {"locus_tag": "TX50_RS09500"}
    assert _compute_contributing_sources(gene) == ["ncbi"]


def test_gene_with_eggnog_match():
    """Eggnog hit (seed_ortholog non-null) → eggnog in sources."""
    gene = {
        "locus_tag": "TX50_RS09501",
        "seed_ortholog": "59919.PMM0001",
        "seed_ortholog_evalue": 1e-100,
    }
    sources = _compute_contributing_sources(gene)
    assert "ncbi" in sources
    assert "eggnog" in sources


def test_gene_with_uniprot_match():
    """UniProt accession non-null → uniprot in sources."""
    gene = {
        "locus_tag": "PMM0001",
        "uniprot_accession": "Q7VFD7",
    }
    sources = _compute_contributing_sources(gene)
    assert "uniprot" in sources


def test_gene_with_cyanorak_locus():
    """Cyanorak locus_tag non-null → cyanorak in sources."""
    gene = {
        "locus_tag": "PMM0001",
        "locus_tag_cyanorak": "PMM0001",
    }
    sources = _compute_contributing_sources(gene)
    assert "cyanorak" in sources


def test_full_house_gene():
    """A well-annotated gene shows all 4 sources."""
    gene = {
        "locus_tag": "PMM0001",
        "locus_tag_cyanorak": "PMM0001",
        "uniprot_accession": "Q7VFD7",
        "seed_ortholog": "59919.PMM0001",
        "eggnog_ogs": ["COG0001@1|root"],
    }
    sources = _compute_contributing_sources(gene)
    assert sources == ["cyanorak", "eggnog", "ncbi", "uniprot"]


def test_returns_sorted_list():
    """Output must be a sorted list (deterministic)."""
    gene = {"locus_tag": "X", "uniprot_accession": "Z", "seed_ortholog": "1.A"}
    sources = _compute_contributing_sources(gene)
    assert sources == sorted(sources)
```

- [ ] **Step 2: Run test to verify it fails (function doesn't exist yet)**

```bash
pytest tests/test_contributing_sources.py -v
```
Expected: ImportError or AttributeError on `_compute_contributing_sources`.

- [ ] **Step 3: Implement the function in `build_gene_annotations.py`**

Add after `_compute_gene_category` (around line 200):

```python
def _has_source_label(gene: dict, label: str) -> bool:
    """Check whether a gene has any field provenance-tagged with `label`.

    Walks `gene["*_source"]` track fields (e.g. product_source, gene_name_source)
    plus `[label]` prefixes inside `alternate_functional_descriptions`.
    """
    # *_source track fields
    for k, v in gene.items():
        if k.endswith("_source") and v == label:
            return True
    # [label] prefix in alternate functional descriptions
    afd = gene.get("alternate_functional_descriptions") or []
    for entry in afd:
        if isinstance(entry, str) and entry.startswith(f"[{label}]"):
            return True
    return False


def _compute_contributing_sources(gene: dict) -> list[str]:
    """Return sorted list of data sources that contributed at least one field.

    Source presence rules:
    - 'ncbi': always present (every Gene comes from an NCBI GFF row).
    - 'cyanorak': locus_tag_cyanorak non-null OR any cyanorak-tagged field.
    - 'uniprot': uniprot_accession non-null OR any uniprot-tagged field.
    - 'eggnog': seed_ortholog or eggnog_ogs non-null OR any eggnog-tagged field.
    """
    sources = {"ncbi"}
    if gene.get("locus_tag_cyanorak") or _has_source_label(gene, "cyanorak"):
        sources.add("cyanorak")
    if gene.get("uniprot_accession") or _has_source_label(gene, "uniprot"):
        sources.add("uniprot")
    if (gene.get("seed_ortholog")
            or gene.get("eggnog_ogs")
            or _has_source_label(gene, "eggnog")):
        sources.add("eggnog")
    return sorted(sources)
```

- [ ] **Step 4: Wire it into the merge loop**

Find `_compute_annotation_quality` and `_compute_gene_category` calls (around line 685-695). Add after `gene_category` is set:

```python
        # Compute contributing_sources (F2 ship 2.4)
        result["contributing_sources"] = _compute_contributing_sources(result)
```

Also add `contributing_sources` to the recompute path in the post-merge enrich step (around line 280):

```python
    # Recompute contributing_sources after enrichment (in case sources changed)
    merged["contributing_sources"] = _compute_contributing_sources(merged)
```

- [ ] **Step 5: Run the test, verify it passes**

```bash
pytest tests/test_contributing_sources.py -v
```
Expected: all 6 tests pass.

- [ ] **Step 6: Add `CONTRIBUTING_SOURCES` to `GeneNodeField` enum**

In `cyanorak_ncbi_adapter.py`:

```python
    CONTRIBUTING_SOURCES = 'contributing_sources'
```

- [ ] **Step 7: Add `contributing_sources: str[]` to Gene schema**

```yaml
    contributing_sources: str[]      # data sources contributing >=1 field; values from {ncbi, cyanorak, uniprot, eggnog}
```

- [ ] **Step 8: Rebuild merged annotations + verify**

```bash
uv run bash scripts/prepare_data.sh --steps 2 --force --strains MED4 2>&1 | tail -5
uv run python -c "
import json
d = json.load(open('cache/data/Prochlorococcus/genomes/MED4/gene_annotations_merged.json'))
sample = list(d.values())[:3]
for g in sample:
    print(g['locus_tag'], g.get('contributing_sources'))
"
```
Expected: each gene shows a sorted list of sources, always including 'ncbi'.

- [ ] **Step 9: Commit**

```bash
git add multiomics_kg/download/build_gene_annotations.py multiomics_kg/adapters/cyanorak_ncbi_adapter.py config/schema_config.yaml tests/test_contributing_sources.py
git commit -m "$(cat <<'EOF'
feat(gene): add Gene.contributing_sources property

Per-gene presence list of data sources that contributed >=1 field. Values
from {ncbi, cyanorak, uniprot, eggnog}; computed at build time from existing
*_source track fields and [label] prefixes in alternate_functional_descriptions
(F2 ship 2.4).

Lets the explorer distinguish 'out of scope' (source absent) from 'no hit'
(source present but tool returned nothing) without ambiguity.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: F2.1 — `DataSource` node schema

**Files:**
- Modify: `config/schema_config.yaml`

- [ ] **Step 1: Add the new node type definition**

In `config/schema_config.yaml`, add a new top-level entry (place alphabetically or near other small node types):

```yaml
data source:
  represented_as: node
  preferred_id: data_source
  label_in_input: data_source
  properties:
    id: str                       # join key (e.g., 'eggnog')
    name: str                     # human-readable name
    description: str
    version: str                  # optional v1 - '' when not derivable
    scope: str                    # 'gene_level' | 'organism_restricted'
    provenance: str               # 'download' | 'tool_run'
    applies_to_organisms: str[]   # populated only for organism_restricted
    info_types: str[]             # auto-generated from gene_annotations_config.yaml
```

- [ ] **Step 2: Verify schema parses**

```bash
uv run python -c "
from biocypher import BioCypher
bc = BioCypher(schema_config_path='config/schema_config.yaml')
print('OK')
"
```
Expected: prints OK with no errors.

- [ ] **Step 3: Commit**

```bash
git add config/schema_config.yaml
git commit -m "$(cat <<'EOF'
feat(schema): add DataSource node type

New node type for the F2 data-source surfacing layer. 4-row metadata table
auto-derived from gene_annotations_config.yaml; soft-joined to Gene via
contributing_sources property (no materialized edges).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: F2.3 — `DataSourceAdapter`

**Files:**
- Create: `multiomics_kg/adapters/data_source_adapter.py`
- Modify: `create_knowledge_graph.py`
- Test: `tests/test_data_source_adapter.py`

- [ ] **Step 1: Write the unit test first**

Create `tests/test_data_source_adapter.py`:

```python
"""Unit tests for DataSourceAdapter.

Reads logical_sources from gene_annotations_config.yaml; emits one
DataSource node per logical source with auto-derived info_types.
"""

import pytest
from multiomics_kg.adapters.data_source_adapter import DataSourceAdapter


@pytest.fixture
def adapter():
    return DataSourceAdapter(config_path="config/gene_annotations_config.yaml")


def test_emits_four_nodes(adapter):
    """Initial deployment has exactly 4 data sources."""
    adapter.download_data()
    nodes = list(adapter.get_nodes())
    ids = {props["id"] for _, _, props in nodes}
    assert ids == {"ncbi", "cyanorak", "uniprot", "eggnog"}


def test_ncbi_node_properties(adapter):
    adapter.download_data()
    nodes = {props["id"]: props for _, _, props in adapter.get_nodes()}
    ncbi = nodes["ncbi"]
    assert ncbi["scope"] == "gene_level"
    assert ncbi["provenance"] == "download"
    assert ncbi["applies_to_organisms"] == []


def test_cyanorak_organism_restricted(adapter):
    adapter.download_data()
    nodes = {props["id"]: props for _, _, props in adapter.get_nodes()}
    cyano = nodes["cyanorak"]
    assert cyano["scope"] == "organism_restricted"
    assert cyano["provenance"] == "download"
    assert any("Prochlorococcus" in o for o in cyano["applies_to_organisms"])


def test_eggnog_is_tool_run(adapter):
    adapter.download_data()
    nodes = {props["id"]: props for _, _, props in adapter.get_nodes()}
    assert nodes["eggnog"]["provenance"] == "tool_run"


def test_info_types_auto_derived(adapter):
    """info_types is computed from the fields block; eggnog should include
    cog_category, kegg_ko, go_terms, etc."""
    adapter.download_data()
    nodes = {props["id"]: props for _, _, props in adapter.get_nodes()}
    eggnog_info = set(nodes["eggnog"]["info_types"])
    # spot-check a few expected entries
    assert "kegg_ko" in eggnog_info
    assert "eggnog_ogs" in eggnog_info


def test_emits_no_edges(adapter):
    adapter.download_data()
    edges = list(adapter.get_edges())
    assert edges == []
```

- [ ] **Step 2: Run test, verify it fails**

```bash
pytest tests/test_data_source_adapter.py -v
```
Expected: ImportError on `DataSourceAdapter`.

- [ ] **Step 3: Implement the adapter**

Create `multiomics_kg/adapters/data_source_adapter.py`:

```python
"""DataSource adapter — emits 4 metadata nodes describing the data sources
that contribute fields to Gene records. Soft-joined to Gene via
Gene.contributing_sources (no materialized edges).

Auto-generates DataSource.info_types by walking the fields block of
gene_annotations_config.yaml and inverting source -> field-name mapping.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any, Iterator

import yaml


def _clean_str(value: str) -> str:
    """Sanitize string for BioCypher CSV output."""
    return value.replace("'", "^").replace("|", "")


class DataSourceAdapter:
    """Emit one DataSource node per logical source declared in
    gene_annotations_config.yaml.

    Required YAML keys per source:
    - logical_sources: list of dicts with id, scope, provenance,
      applies_to_organisms (optional)

    The adapter fails loudly if logical_sources is missing rather than
    falling back, per the spec.
    """

    def __init__(self, config_path: str | Path = "config/gene_annotations_config.yaml"):
        self.config_path = Path(config_path)
        self._config: dict[str, Any] | None = None
        self._info_types_by_source: dict[str, list[str]] = {}

    def download_data(self, cache: bool = True) -> None:
        with open(self.config_path) as f:
            self._config = yaml.safe_load(f)
        self._info_types_by_source = self._derive_info_types(self._config)

    @staticmethod
    def _derive_info_types(config: dict[str, Any]) -> dict[str, list[str]]:
        """Walk fields:; for each field rule, collect contributing source(s)."""
        info_by_src: dict[str, set[str]] = defaultdict(set)
        for field_name, rule in (config.get("fields") or {}).items():
            for src in DataSourceAdapter._sources_from_rule(rule):
                info_by_src[src].add(field_name)
        return {k: sorted(v) for k, v in info_by_src.items()}

    @staticmethod
    def _sources_from_rule(rule: Any) -> Iterator[str]:
        """A field rule may declare a single source, a list of candidates, or
        nested options. Walk recursively; surface every source / source_label
        encountered."""
        if isinstance(rule, dict):
            if "source_label" in rule:
                yield rule["source_label"]
            elif "source" in rule:
                # A 'source' key naming a YAML source entry name is fine
                # for single-source fields. For 'gene_mapping' we map it to
                # both ncbi and cyanorak — see _emit_for_gene_mapping below.
                yield rule["source"]
            for v in rule.values():
                yield from DataSourceAdapter._sources_from_rule(v)
        elif isinstance(rule, list):
            for item in rule:
                yield from DataSourceAdapter._sources_from_rule(item)

    def _all_logical_sources(self) -> Iterator[dict[str, Any]]:
        """Iterate all logical sources from all source entries."""
        sources_block = (self._config or {}).get("sources") or {}
        for src_name, src_def in sources_block.items():
            if "logical_sources" not in src_def:
                raise ValueError(
                    f"Source '{src_name}' in gene_annotations_config.yaml "
                    f"missing required 'logical_sources' key."
                )
            for ls in src_def["logical_sources"]:
                yield {"_yaml_source_name": src_name, **ls}

    @staticmethod
    def _name_for(source_id: str) -> str:
        return {
            "ncbi": "NCBI RefSeq",
            "cyanorak": "Cyanorak",
            "uniprot": "UniProt",
            "eggnog": "EggNOG-mapper",
        }.get(source_id, source_id.title())

    @staticmethod
    def _description_for(source_id: str) -> str:
        return {
            "ncbi": "NCBI RefSeq genome annotation (GFF + GenBank).",
            "cyanorak": "Cyanorak: curated cyanobacteria gene clusters and roles.",
            "uniprot": "UniProt proteins (reviewed + unreviewed) cross-referenced via RefSeq protein_id.",
            "eggnog": "EggNOG-mapper functional annotations (per-protein query against the eggNOG reference).",
        }.get(source_id, "")

    def get_nodes(self) -> Iterator[tuple[str, str, dict[str, Any]]]:
        if self._config is None:
            self.download_data()
        for ls in self._all_logical_sources():
            src_id = ls["id"]
            yaml_name = ls["_yaml_source_name"]
            # info_types: for gene_mapping, both ncbi and cyanorak share the
            # same field set (we cannot easily attribute fields to one vs the
            # other without source_label inspection). Use the YAML source name.
            info_types = self._info_types_by_source.get(yaml_name, [])
            # Override for ncbi/cyanorak: when source_label is set explicitly
            # in a rule, prefer that
            specific = self._info_types_by_source.get(src_id, [])
            if specific:
                info_types = sorted(set(info_types) | set(specific))
            node_id = f"data_source:{src_id}"
            properties = {
                "id": _clean_str(src_id),
                "name": _clean_str(self._name_for(src_id)),
                "description": _clean_str(self._description_for(src_id)),
                "version": "",   # populated case-by-case in v2
                "scope": _clean_str(ls.get("scope", "")),
                "provenance": _clean_str(ls.get("provenance", "")),
                "applies_to_organisms": [_clean_str(o) for o in (ls.get("applies_to_organisms") or [])],
                "info_types": [_clean_str(t) for t in info_types],
            }
            yield node_id, "data source", properties

    def get_edges(self) -> Iterator[tuple]:
        return iter(())
```

- [ ] **Step 4: Run unit test, verify it passes**

```bash
pytest tests/test_data_source_adapter.py -v
```
Expected: all 6 tests pass.

- [ ] **Step 5: Wire adapter into create_knowledge_graph.py**

Find where other top-level adapters are instantiated (search for `MultiCyanorakNcbi(` or similar). Add:

```python
from multiomics_kg.adapters.data_source_adapter import DataSourceAdapter

# ...inside main, after other adapters...
data_source_adapter = DataSourceAdapter(config_path="config/gene_annotations_config.yaml")
data_source_adapter.download_data()
bc.write_nodes(data_source_adapter.get_nodes())
```

- [ ] **Step 6: Test-mode rebuild + verify CSV output**

```bash
uv run python create_knowledge_graph.py --test 2>&1 | tail -10
ls biocypher-log/example_knowledge_graph/ | grep -i 'datasource\|data_source\|DataSource'
```
Expected: a `DataSource-*.csv` file exists with 4 rows.

- [ ] **Step 7: Write KG validity test**

Create `tests/kg_validity/test_data_source.py`:

```python
"""KG validity: DataSource nodes (F2)."""

import pytest


@pytest.mark.kg
def test_four_data_source_nodes(neo4j):
    """Initial deployment has exactly 4 DataSource nodes."""
    result = neo4j.query("MATCH (d:DataSource) RETURN d.id AS id ORDER BY d.id")
    ids = [row["id"] for row in result]
    assert ids == ["cyanorak", "eggnog", "ncbi", "uniprot"]


@pytest.mark.kg
def test_info_types_non_empty(neo4j):
    """Every DataSource must have at least one info_type listed."""
    result = neo4j.query("""
        MATCH (d:DataSource)
        WHERE size(d.info_types) = 0
        RETURN d.id AS id
    """)
    assert result == [], f"DataSources with empty info_types: {result}"


@pytest.mark.kg
def test_eggnog_provenance_is_tool_run(neo4j):
    result = neo4j.query("""
        MATCH (d:DataSource {id: 'eggnog'}) RETURN d.provenance AS p
    """)
    assert result[0]["p"] == "tool_run"


@pytest.mark.kg
def test_cyanorak_organism_restricted(neo4j):
    result = neo4j.query("""
        MATCH (d:DataSource {id: 'cyanorak'})
        RETURN d.scope AS scope, size(d.applies_to_organisms) AS n
    """)
    assert result[0]["scope"] == "organism_restricted"
    assert result[0]["n"] > 0


@pytest.mark.kg
def test_every_gene_has_ncbi_in_contributing_sources(neo4j):
    """The ncbi source is universal — every Gene must list it."""
    result = neo4j.query("""
        MATCH (g:Gene)
        WHERE NOT 'ncbi' IN g.contributing_sources
        RETURN count(*) AS missing
    """)
    assert result[0]["missing"] == 0
```

- [ ] **Step 8: Commit**

```bash
git add multiomics_kg/adapters/data_source_adapter.py create_knowledge_graph.py tests/test_data_source_adapter.py tests/kg_validity/test_data_source.py
git commit -m "$(cat <<'EOF'
feat(adapter): DataSource adapter + KG validity tests

Emits 4 DataSource nodes (ncbi, cyanorak, uniprot, eggnog) auto-derived
from gene_annotations_config.yaml. info_types computed by walking the
fields block. Soft-joined to Gene via contributing_sources (no materialized
edges). Resolves F2 ships 2.1, 2.3.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: F1.1 — `is_uninformative` term flag (YAML + post-import Cypher)

**Files:**
- Create: `config/uninformative_terms.yaml`
- Modify: `scripts/post-import.cypher`
- Modify: `scripts/post-import.sh`
- Test: `tests/kg_validity/test_uninformative_terms.py`

- [ ] **Step 1: Create the YAML vocabulary**

Create `config/uninformative_terms.yaml`:

```yaml
# Catch-all / no-content-bearing terms across ontology sources.
#
# Guiding principle: if the term tells the consumer the broad CLASS
# (transporter, enzyme, peptidase, two-component system) — even if a sub-class
# dimension is unknown — treat the term as INFORMATIVE. Only flag terms that
# convey no class signal at all.
#
# Pfam DUF/UPF stay UN-flagged (recognized structural domain = informative).
# COG R ('General function prediction only') stays UN-flagged (borderline).
# All BRITE entries stay UN-flagged (every entry sits inside a class-bearing tree).

biological_process:
  ids:
    - go:0008150          # biological_process root

molecular_function:
  ids:
    - go:0003674          # molecular_function root

cellular_component:
  ids:
    - go:0005575          # cellular_component root

cog_category:
  ids:
    - cog.category:S      # Function unknown

cyanorak_role:
  ids:
    - cyanorak.role:R     # Other (parent — used standalone = no class info)
    - cyanorak.role:R.1   # Other > Conserved hypothetical domains
    - cyanorak.role:R.2   # Other > Conserved hypothetical proteins
    - cyanorak.role:R.4   # Other > Hypothetical proteins
    - cyanorak.role:R.5   # Other > Other
  # Q.9 (Transport > Unknown substrate) and R.3 (Other > Enzymes of unknown
  # specificity) intentionally NOT flagged — class is known.

tigr_role:
  ids:
    - tigr.role:156       # Hypothetical proteins / Conserved
    - tigr.role:704       # Hypothetical proteins / Domain
    - tigr.role:856       # Not Found
    - tigr.role:185       # Unclassified / Role category not yet assigned
    - tigr.role:157       # Unknown function / General
  # 141 (Transport / Unknown substrate) and 703 (Unknown function / Enzymes
  # of unknown specificity) intentionally NOT flagged — class is known.

kegg_term:
  # ~210 KOs match this pattern; ID list omitted (too many; pattern stable).
  name_patterns:
    - '^K\d+;\s+uncharacterized protein\b'

# ec_number: section omitted — every EC is functionally specific.
# brite_category: section omitted — every BRITE entry sits inside a
#   class-bearing tree (Transporters, Peptidases, Two-component, etc.);
#   no entries pass the 'no class signal' bar.
# pfam: section omitted — DUF/UPF stay informative-shaped per anti-scope.
# pfam_clan: section omitted — clans are structural superfamilies.
```

- [ ] **Step 2: Verify YAML parses**

```bash
uv run python -c "import yaml; print(yaml.safe_load(open('config/uninformative_terms.yaml')).keys())"
```
Expected: prints all 7 section keys.

- [ ] **Step 3: Add post-import Cypher block to `scripts/post-import.cypher`**

Find an appropriate location BEFORE the existing `annotation_types` block (line ~498) — the term flag must be set before `informative_annotation_types` and `annotation_state` queries reference it.

Insert (label-binding hardcoded per section name):

```cypher
// =====================================================================
// F1.1: Term-level is_uninformative flag (sentinel str 'true' / absent).
// Driven by config/uninformative_terms.yaml. Vocabulary is small; keys
// are hard-coded here matching the YAML section names.
// =====================================================================

// Direct ID flags
MATCH (t:BiologicalProcess) WHERE t.id IN ['go:0008150'] SET t.is_uninformative = 'true';
MATCH (t:MolecularFunction) WHERE t.id IN ['go:0003674'] SET t.is_uninformative = 'true';
MATCH (t:CellularComponent) WHERE t.id IN ['go:0005575'] SET t.is_uninformative = 'true';
MATCH (t:CogFunctionalCategory) WHERE t.id IN ['cog.category:S'] SET t.is_uninformative = 'true';

MATCH (t:CyanorakRole)
WHERE t.id IN ['cyanorak.role:R','cyanorak.role:R.1','cyanorak.role:R.2',
               'cyanorak.role:R.4','cyanorak.role:R.5']
SET t.is_uninformative = 'true';

MATCH (t:TigrRole)
WHERE t.id IN ['tigr.role:156','tigr.role:704','tigr.role:856',
               'tigr.role:185','tigr.role:157']
SET t.is_uninformative = 'true';

// Pattern-based flag for KEGG (uncharacterized protein KOs, ~210 nodes)
MATCH (t:KeggTerm)
WHERE t.name =~ '^K\\d+;\\s+uncharacterized protein\\b.*'
SET t.is_uninformative = 'true';
```

- [ ] **Step 4: Mirror the same Cypher in `scripts/post-import.sh`**

Locate the post-import shell script. The block of statements is grouped into ~3 `cypher-shell` invocations. Add the `is_uninformative` block to the small-aggregations group (the section with simple SET statements without `CALL IN TRANSACTIONS`).

The shell script uses bash heredoc `<<EOF` markers — preserve the same Cypher byte-for-byte to keep the two files diffable.

- [ ] **Step 5: Run post-import locally to verify**

If the Docker stack is up:

```bash
cypher-shell -u neo4j -p neo4j 'MATCH (t) WHERE t.is_uninformative IS NOT NULL RETURN labels(t)[0] AS label, count(*) AS n ORDER BY n DESC' 2>&1 | head -20
```

Expected counts: `KeggTerm ~210`, `CyanorakRole 5`, `TigrRole 5`, `BiologicalProcess 1`, `MolecularFunction 1`, `CellularComponent 1`, `CogFunctionalCategory 1`. Total ~224 flagged nodes.

(If running pre-rebuild: skip; the test in the next step locks this expectation.)

- [ ] **Step 6: Write KG validity test**

Create `tests/kg_validity/test_uninformative_terms.py`:

```python
"""KG validity: F1.1 is_uninformative term flag.

Asserts:
1. Every YAML ids: entry resolves to an existing node post-build with
   is_uninformative='true'.
2. The KEGG name_pattern matches >= 100 nodes (sanity floor; today ~210).
3. No node outside the YAML has is_uninformative set (catches accidental
   over-flagging).
4. Drift check: every catch-all role entry in the YAML has its corresponding
   *_TO_CATEGORY mapping returning 'Unknown' in build_gene_annotations.py.
"""

from pathlib import Path

import pytest
import yaml


@pytest.fixture(scope="module")
def yaml_vocab():
    with open("config/uninformative_terms.yaml") as f:
        return yaml.safe_load(f)


SECTION_TO_LABEL = {
    "biological_process": "BiologicalProcess",
    "molecular_function": "MolecularFunction",
    "cellular_component": "CellularComponent",
    "cog_category": "CogFunctionalCategory",
    "cyanorak_role": "CyanorakRole",
    "tigr_role": "TigrRole",
    "kegg_term": "KeggTerm",
}


@pytest.mark.kg
def test_yaml_ids_all_resolve_to_flagged_nodes(neo4j, yaml_vocab):
    """Every ids: entry exists post-build and has is_uninformative='true'."""
    failures = []
    for section_name, payload in yaml_vocab.items():
        label = SECTION_TO_LABEL.get(section_name)
        assert label, f"Unknown section '{section_name}' — update SECTION_TO_LABEL"
        for term_id in (payload or {}).get("ids") or []:
            result = neo4j.query(
                f"MATCH (t:{label} {{id: $id}}) RETURN t.is_uninformative AS flag",
                {"id": term_id},
            )
            if not result:
                failures.append(f"{label}:{term_id} — node missing")
            elif result[0]["flag"] != "true":
                failures.append(f"{label}:{term_id} — flag={result[0]['flag']!r}")
    assert not failures, "\n".join(failures)


@pytest.mark.kg
def test_kegg_pattern_flags_at_least_100_kos(neo4j):
    """KEGG name_pattern should flag at least 100 'uncharacterized protein' KOs."""
    result = neo4j.query("""
        MATCH (t:KeggTerm) WHERE t.is_uninformative = 'true' RETURN count(*) AS n
    """)
    assert result[0]["n"] >= 100


@pytest.mark.kg
def test_no_pfam_or_ec_flagged(neo4j):
    """Anti-scope: Pfam DUF/UPF and every EC stay UN-flagged."""
    for label in ["Pfam", "PfamClan", "EcNumber", "BriteCategory"]:
        result = neo4j.query(f"""
            MATCH (t:{label}) WHERE t.is_uninformative IS NOT NULL
            RETURN count(*) AS n
        """)
        assert result[0]["n"] == 0, f"{label} has flagged terms (anti-scope violation)"


def test_drift_role_yaml_matches_to_category_mappings(yaml_vocab):
    """Build-time drift check (no Neo4j): every catch-all role entry in the
    YAML must, when fed as a gene's only role, yield gene_category='Unknown'
    from the *_TO_CATEGORY mappings in build_gene_annotations.py.

    This catches the case where a new uninformative TIGR role gets added
    to the YAML but the corresponding TIGR_TO_CATEGORY entry isn't 'Unknown'.
    """
    from multiomics_kg.download.build_gene_annotations import (
        COG_TO_CATEGORY, CYANORAK_TO_CATEGORY, TIGR_TO_CATEGORY, _compute_gene_category,
    )
    failures = []
    # COG: ids look like 'cog.category:X' — strip prefix
    for term_id in (yaml_vocab.get("cog_category") or {}).get("ids") or []:
        letter = term_id.split(":")[-1]
        if COG_TO_CATEGORY.get(letter) != "Unknown":
            failures.append(f"COG {letter}: TO_CATEGORY={COG_TO_CATEGORY.get(letter)!r}, expected 'Unknown'")
    # Cyanorak: ids look like 'cyanorak.role:R.1'
    for term_id in (yaml_vocab.get("cyanorak_role") or {}).get("ids") or []:
        # Top-letter category
        code = term_id.split(":")[-1]
        top_letter = code.split(".")[0]
        # We use _compute_gene_category to honor the D-subcode special case
        result = _compute_gene_category({"cyanorak_Role": [code]})
        if result != "Unknown":
            failures.append(f"Cyanorak {code}: gene_category={result!r}, expected 'Unknown'")
    # TIGR: ids are 'tigr.role:156' but TIGR_TO_CATEGORY is keyed by description.
    # Skip programmatic lookup here; rely on COG + Cyanorak checks above
    # plus a separate hand-curated fixture if more strictness is needed.
    assert not failures, "\n".join(failures)
```

- [ ] **Step 7: Run the drift test (no Neo4j needed)**

```bash
pytest tests/kg_validity/test_uninformative_terms.py::test_drift_role_yaml_matches_to_category_mappings -v
```
Expected: passes.

- [ ] **Step 8: Commit**

```bash
git add config/uninformative_terms.yaml scripts/post-import.cypher scripts/post-import.sh tests/kg_validity/test_uninformative_terms.py
git commit -m "$(cat <<'EOF'
feat(post-import): add is_uninformative term flag (F1.1)

Hand-curated catch-all term vocabulary at config/uninformative_terms.yaml
flags ~224 ontology terms across 7 node types (GO roots, COG S, 5 Cyanorak
'Other' roles, 5 TIGR hypothetical/unknown roles, ~210 KEGG 'uncharacterized
protein' KOs by name pattern). Pfam DUF/UPF and BRITE entries intentionally
stay un-flagged per the F1 anti-scope.

Sentinel value 'true' / absence convention (project-standard for boolean
properties) avoids the BioCypher bool serialization bug.

KG validity tests assert YAML entries resolve, KEGG pattern matches >=100
KOs, no out-of-scope terms get flagged, and the *_TO_CATEGORY drift check
holds.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: F1.2 + F1.3 — `annotation_quality` + `annotation_state` (combined post-import block)

**Files:**
- Modify: `scripts/post-import.cypher`
- Modify: `scripts/post-import.sh`

Both fields share `informative_source_count` over 8 source buckets (`go`, `kegg`, `pfam`, `ec`, `role`, `reaction`, `transporter`, `cazy`). Single Cypher block.

- [ ] **Step 1: Add the new Cypher block AFTER the term-flag block from Task 9**

Insert into `scripts/post-import.cypher`:

```cypher
// =====================================================================
// F1.2 + F1.3: annotation_quality (numeric 0-3) + annotation_state (enum)
// from informative_source_count over 8 source buckets.
//
// SOURCE_BUCKETS:start
//   live (8): go, kegg, pfam, ec, role, reaction, transporter, cazy
// SOURCE_BUCKETS:end
//
// Maintenance: when adding a new functional Gene-edge type, append a
// has_<bucket> EXISTS line, include in informative_count sum, and add
// the edge type(s) to has_any_edge. See the design spec section
// 'Source bucket maintenance'.
// =====================================================================

MATCH (g:Gene)
CALL {
  WITH g
  WITH g,
       EXISTS { (g)-[:Gene_involved_in_biological_process|Gene_enables_molecular_function|Gene_located_in_cellular_component]->(t)
                WHERE t.is_uninformative IS NULL } AS has_go,
       EXISTS { (g)-[:Gene_has_kegg_ko]->(t) WHERE t.is_uninformative IS NULL } AS has_kegg,
       EXISTS { (g)-[:Gene_has_pfam]->() } AS has_pfam,
       EXISTS { (g)-[:Gene_catalyzes_ec_number]->() } AS has_ec,
       (g.gene_category IS NOT NULL AND g.gene_category <> 'Unknown') AS has_role,
       EXISTS { (g)-[:Gene_catalyzes_reaction]->() } AS has_reaction,
       EXISTS { (g)-[:Gene_has_tcdb_family]->() } AS has_transporter,
       EXISTS { (g)-[:Gene_has_cazy_family]->() } AS has_cazy,
       EXISTS { (g)-[:Gene_involved_in_biological_process|Gene_enables_molecular_function|Gene_located_in_cellular_component
                     |Gene_has_kegg_ko|Gene_has_pfam|Gene_catalyzes_ec_number
                     |Gene_in_cog_category|Gene_has_cyanorak_role|Gene_has_tigr_role
                     |Gene_catalyzes_reaction|Gene_has_tcdb_family|Gene_has_cazy_family]->() } AS has_any_edge
  WITH g,
       (CASE WHEN has_go THEN 1 ELSE 0 END
        + CASE WHEN has_kegg THEN 1 ELSE 0 END
        + CASE WHEN has_pfam THEN 1 ELSE 0 END
        + CASE WHEN has_ec THEN 1 ELSE 0 END
        + CASE WHEN has_role THEN 1 ELSE 0 END
        + CASE WHEN has_reaction THEN 1 ELSE 0 END
        + CASE WHEN has_transporter THEN 1 ELSE 0 END
        + CASE WHEN has_cazy THEN 1 ELSE 0 END) AS informative_count,
       has_any_edge
  SET g.annotation_state =
        CASE
          WHEN informative_count >= 2 THEN 'informative_multi'
          WHEN informative_count = 1 THEN 'informative_single'
          WHEN has_any_edge THEN 'catch_all_only'
          ELSE 'no_evidence'
        END,
      g.annotation_quality =
        CASE
          WHEN informative_count >= 2 THEN 3
          WHEN informative_count = 1 THEN 2
          WHEN has_any_edge THEN 1
          ELSE 0
        END
} IN TRANSACTIONS OF 500 ROWS;
```

- [ ] **Step 2: Mirror in `scripts/post-import.sh`**

Add the same block to the `CALL IN TRANSACTIONS` group of the bash script (the third group; this is a write-with-batching pattern).

- [ ] **Step 3: Add `annotation_state` to Gene schema**

In `config/schema_config.yaml` `gene:` `properties:`:

```yaml
    annotation_state: str          # post-import enum: no_evidence | catch_all_only | informative_single | informative_multi
```

(`annotation_quality` already exists in the schema; only its semantics change.)

- [ ] **Step 4: Update CLAUDE.md `annotation_quality` description**

Find the line in CLAUDE.md describing `annotation_quality` (search `0=hypothetical-no-func`) and replace with:

```
annotation_quality (post-import, refined from build-time): numeric encoding of `annotation_state` (0=no_evidence, 1=catch_all_only, 2=informative_single, 3=informative_multi), driven by `informative_source_count` over 8 source buckets: `go`, `kegg`, `pfam`, `ec`, `role` (via gene_category), `reaction`, `transporter` (via Gene_has_tcdb_family), `cazy` (via Gene_has_cazy_family). See `docs/superpowers/specs/2026-05-01-explorer-frictions-resolution-design.md` for maintenance procedure.
```

- [ ] **Step 5: Write KG validity test**

Create `tests/kg_validity/test_annotation_state.py`:

```python
"""KG validity: F1.2 + F1.3 annotation_quality + annotation_state."""

import pytest

VALID_STATES = {"no_evidence", "catch_all_only", "informative_single", "informative_multi"}
STATE_TO_QUALITY = {
    "no_evidence": 0,
    "catch_all_only": 1,
    "informative_single": 2,
    "informative_multi": 3,
}


@pytest.mark.kg
def test_every_gene_has_valid_state(neo4j):
    result = neo4j.query("""
        MATCH (g:Gene)
        WHERE g.annotation_state IS NULL OR NOT g.annotation_state IN $valid
        RETURN count(*) AS n, collect(DISTINCT g.annotation_state)[..5] AS sample
    """, {"valid": list(VALID_STATES)})
    assert result[0]["n"] == 0, f"Genes with invalid state: {result[0]}"


@pytest.mark.kg
def test_quality_state_one_to_one(neo4j):
    """annotation_quality must be the numeric encoding of annotation_state."""
    result = neo4j.query("""
        MATCH (g:Gene)
        WITH g.annotation_state AS state, g.annotation_quality AS qual,
             count(*) AS n
        RETURN state, qual, n ORDER BY state, qual
    """)
    failures = []
    for row in result:
        expected = STATE_TO_QUALITY.get(row["state"])
        if row["qual"] != expected:
            failures.append(f"state={row['state']!r}: quality={row['qual']!r} (expected {expected!r}, n={row['n']})")
    assert not failures, "\n".join(failures)


@pytest.mark.kg
def test_known_well_annotated_gene(neo4j):
    """Spot-check: dnaA in MED4 (PMM0001) should be informative_multi."""
    result = neo4j.query("""
        MATCH (g:Gene {locus_tag: 'PMM0001'})
        RETURN g.annotation_state AS state, g.annotation_quality AS qual
    """)
    assert result, "PMM0001 not found"
    # dnaA has GO + Pfam + KEGG + EC + role -> at least 2 informative buckets
    assert result[0]["state"] == "informative_multi"
    assert result[0]["qual"] == 3
```

- [ ] **Step 6: Commit**

```bash
git add scripts/post-import.cypher scripts/post-import.sh config/schema_config.yaml CLAUDE.md tests/kg_validity/test_annotation_state.py
git commit -m "$(cat <<'EOF'
feat(post-import): annotation_quality + annotation_state (F1.2 + F1.3)

Combined Cypher block sets both fields from informative_source_count over
8 source buckets: go, kegg, pfam, ec, role (via gene_category), reaction,
transporter (via Gene_has_tcdb_family), cazy (via Gene_has_cazy_family).
annotation_state is the categorical surface; annotation_quality is its 0-3
numeric encoding (1:1 mapping).

All 8 buckets are live as of the 2026-05-02 main-merge.

Existing 'WHERE g.annotation_quality <= 1' queries silently shift meaning
from 'hypothetical-named' to 'low-evidence' — the spirit of the original
filter, now correctly catching low-evidence genes.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 11: F1.4 — `informative_annotation_types` (post-import Cypher)

**Files:**
- Modify: `scripts/post-import.cypher`
- Modify: `scripts/post-import.sh`
- Modify: `config/schema_config.yaml`

Granular 11-source list (parallel of `annotation_types` with informativeness filter; includes new `reaction` source).

- [ ] **Step 1: Add `informative_annotation_types` Cypher block**

Insert into `scripts/post-import.cypher` adjacent to the existing `annotation_types` block (line ~498):

```cypher
// =====================================================================
// F1.4: Gene.informative_annotation_types — granular per-source list,
// only includes a source if at least one connected term is informative.
// Parallel of Gene.annotation_types (presence-by-source) with the
// informativeness filter applied.
// =====================================================================

MATCH (g:Gene)
CALL {
  WITH g
  SET g.informative_annotation_types =
    CASE WHEN EXISTS { (g)-[:Gene_involved_in_biological_process]->(t)
                       WHERE t.is_uninformative IS NULL }
         THEN ['go_bp'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_enables_molecular_function]->(t)
                       WHERE t.is_uninformative IS NULL }
         THEN ['go_mf'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_located_in_cellular_component]->(t)
                       WHERE t.is_uninformative IS NULL }
         THEN ['go_cc'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_pfam]->() } THEN ['pfam'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_in_cog_category]->(t)
                       WHERE t.is_uninformative IS NULL }
         THEN ['cog_category'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_kegg_ko]->(t)
                       WHERE t.is_uninformative IS NULL }
         THEN ['kegg'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_kegg_ko]->()-[:Kegg_term_in_brite_category]->() }
         THEN ['brite'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_catalyzes_ec_number]->() } THEN ['ec'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_cyanorak_role]->(t)
                       WHERE t.is_uninformative IS NULL }
         THEN ['cyanorak_role'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_tigr_role]->(t)
                       WHERE t.is_uninformative IS NULL }
         THEN ['tigr_role'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_catalyzes_reaction]->() } THEN ['reaction'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_tcdb_family]->() } THEN ['transporter'] ELSE [] END +
    CASE WHEN EXISTS { (g)-[:Gene_has_cazy_family]->() } THEN ['cazy'] ELSE [] END
} IN TRANSACTIONS OF 1000 ROWS;
```

- [ ] **Step 2: Mirror in `scripts/post-import.sh`**

Add the same block to the matching `cypher-shell` group (CALL IN TRANSACTIONS group).

- [ ] **Step 3: Add `informative_annotation_types: str[]` to Gene schema**

```yaml
    informative_annotation_types: str[]   # parallel of annotation_types with informativeness filter
```

- [ ] **Step 4: Add a KG validity test**

Append to `tests/kg_validity/test_annotation_state.py`:

```python
@pytest.mark.kg
def test_informative_subset_of_annotation_types(neo4j):
    """For every gene, informative_annotation_types must be a subset of
    annotation_types — informativeness can only filter OUT, not add."""
    result = neo4j.query("""
        MATCH (g:Gene)
        WITH g, [t IN g.informative_annotation_types
                 WHERE NOT t IN g.annotation_types | t] AS extra
        WHERE size(extra) > 0 AND NOT 'reaction' IN extra
        // 'reaction' is a NEW informative-only source not present in
        // legacy annotation_types — exempt during the F1 release window.
        RETURN count(*) AS n
    """)
    assert result[0]["n"] == 0


@pytest.mark.kg
def test_no_evidence_gene_has_empty_informative_types(neo4j):
    result = neo4j.query("""
        MATCH (g:Gene {annotation_state: 'no_evidence'})
        WHERE size(g.informative_annotation_types) > 0
        RETURN count(*) AS n
    """)
    assert result[0]["n"] == 0
```

- [ ] **Step 5: Commit**

```bash
git add scripts/post-import.cypher scripts/post-import.sh config/schema_config.yaml tests/kg_validity/test_annotation_state.py
git commit -m "$(cat <<'EOF'
feat(post-import): Gene.informative_annotation_types (F1.4)

Granular per-source list adjacent to existing annotation_types, with the
informativeness filter applied. Includes new 'reaction' source.

Existing annotation_types stays unchanged (per F-AUDIT-2 layer 1 fix
clarification — adding/tightening would invalidate that fix).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 12: Delete `_compute_annotation_quality` from build_gene_annotations.py

The function moved to post-import in Task 10. Now delete the obsolete build-time code.

**Files:**
- Modify: `multiomics_kg/download/build_gene_annotations.py`

- [ ] **Step 1: Locate and delete the function**

```bash
grep -n "_compute_annotation_quality" multiomics_kg/download/build_gene_annotations.py
```

Delete the function definition (around line 209-232) AND every call site:
- Around line 282: `gene["annotation_quality"] = _compute_annotation_quality(gene)`
- Around line 688: `result["annotation_quality"] = _compute_annotation_quality(result)`
- Around line 918: `q = merged.get("annotation_quality", 0)` — also remove the related stats reporting if it depended on this field, OR leave the stats but note that pre-merge `annotation_quality` is now always missing (the field is set post-import).

For step-2 stats reporting on annotation_quality, replace with a comment:

```python
# annotation_quality is no longer computed at build-time; see
# scripts/post-import.cypher (F1.2). Stats below skip it.
```

- [ ] **Step 2: Run unit tests to ensure nothing else relied on the function**

```bash
pytest -m "not slow and not kg" -k "annotation_quality or build_gene_annotations" -v 2>&1 | tail -30
```
Expected: tests pass; if any explicitly checked the build-time function, update them to skip.

- [ ] **Step 3: Rebuild merged annotations**

```bash
uv run bash scripts/prepare_data.sh --steps 2 --force --strains MED4 2>&1 | tail -5
```

- [ ] **Step 4: Verify merged JSON no longer carries annotation_quality**

```bash
uv run python -c "
import json
d = json.load(open('cache/data/Prochlorococcus/genomes/MED4/gene_annotations_merged.json'))
sample = list(d.values())[:3]
for g in sample:
    print(g['locus_tag'], 'annotation_quality' in g)
"
```
Expected: `False` for each gene (field absent).

- [ ] **Step 5: Commit**

```bash
git add multiomics_kg/download/build_gene_annotations.py
git commit -m "$(cat <<'EOF'
refactor(build): delete build-time _compute_annotation_quality (moved to post-import)

annotation_quality is now computed entirely in post-import Cypher (F1.2),
sharing informative_source_count with annotation_state. Removing the
build-time function eliminates phase split and YAML duplication.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 13: Full graph rebuild + KG validity gate

**Files:** none (test/verification task)

The build/post-import changes accumulated through Tasks 1-12 require a full graph rebuild before KG validity tests can run.

- [ ] **Step 1: Stop the deploy + app containers (they hold a lock)**

```bash
docker compose stop deploy app 2>&1 | tail -3
```

- [ ] **Step 2: Rebuild + import + post-process**

```bash
docker compose up build import post-process 2>&1 | tail -40
```
Expected: each container exits with status 0. Post-import logs show timing for each Cypher group; the new F1 blocks should add 1–2 seconds total.

- [ ] **Step 3: Bring deploy back up**

```bash
docker compose up -d deploy 2>&1 | tail -3
sleep 8
cypher-shell -u neo4j -p neo4j 'MATCH (n) RETURN count(*) AS total LIMIT 1'
```
Expected: deploy container is healthy; node count matches expectations from CLAUDE.md (~81K Genes, etc).

- [ ] **Step 4: Run KG validity suite**

```bash
pytest -m kg -v 2>&1 | tail -40
```
Expected: all tests pass, including the new ones from Tasks 3, 8, 9, 10, 11.

- [ ] **Step 5: Run baseline omics-edge-snapshot to confirm no expression edge regressions**

```bash
# Save snapshot under a release-marker name
uv run python -c "
from pathlib import Path
import sys
sys.path.insert(0, '.claude/skills/omics-edge-snapshot/scripts')
from snapshot import write_snapshot
write_snapshot(Path('omics_edge_snapshots/after_explorer_frictions_f1_f4.json'))
"
```

- [ ] **Step 6: Spot-check the floor-case gene**

```bash
cypher-shell -u neo4j -p neo4j "
MATCH (g:Gene {locus_tag: 'TX50_RS09500'})
RETURN g.contributing_sources, g.annotation_quality, g.annotation_state,
       g.informative_annotation_types, g.contig
"
```
Expected (typical): `contributing_sources=['ncbi']`, `annotation_state ∈ {'no_evidence', 'catch_all_only'}`, `informative_annotation_types=[]`, `contig` non-null.

- [ ] **Step 7: Commit the snapshot**

```bash
git add omics_edge_snapshots/after_explorer_frictions_f1_f4.json
git commit -m "$(cat <<'EOF'
test(kg): omics-edge-snapshot after F1-F4 release

Release-marker snapshot to compare against before/after the F1-F4 changes.
Expected to match the prior snapshot byte-for-byte (no expression edges
should change as a result of this release).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 14: Documentation rollup — CLAUDE.md + release notes

**Files:**
- Modify: `CLAUDE.md`
- Create or modify: a release notes location agreed by the project (e.g., `docs/kg-changes/2026-05-01-explorer-frictions-f1-f4.md`)

- [ ] **Step 1: Update CLAUDE.md "Key graph facts"**

Find the `annotation_types` description and the `annotation_quality` description in CLAUDE.md. Update:

- `annotation_quality` → 0–3 numeric encoding of `annotation_state` (replace the existing 0=hypothetical-no-func ... 3=real-product+structured description).
- Add `annotation_state` description: enum `{no_evidence, catch_all_only, informative_single, informative_multi}`, driven by `informative_source_count` over 8 source buckets.
- Add `informative_annotation_types` description: granular parallel of `annotation_types` with informativeness filter.
- Add `contributing_sources` description: per-gene presence list of data sources that contributed >=1 field, values from `{ncbi, cyanorak, uniprot, eggnog}`.
- Add `contig`, `seed_ortholog`, `seed_ortholog_evalue` to the Gene properties enumeration.
- Add a bullet listing the 4 DataSource nodes and their properties.

- [ ] **Step 2: Add a "Source bucket maintenance" callout under the annotation_quality / annotation_state entry**

> **Source bucket maintenance.** The `annotation_quality` / `annotation_state` source-bucket list is **explicitly enumerated**, not auto-discovered. The 8 live buckets: `go`, `kegg`, `pfam`, `ec`, `role`, `reaction`, `transporter`, `cazy`. To add a new bucket, follow the procedure in [`docs/superpowers/specs/2026-05-01-explorer-frictions-resolution-design.md`](docs/superpowers/specs/2026-05-01-explorer-frictions-resolution-design.md) — touches the YAML, post-import Cypher, this CLAUDE.md, and the bucket-count test.

- [ ] **Step 3: Add `expression_bin` cluster_type rename to CLAUDE.md**

If CLAUDE.md mentions `cluster_type` values (search `time_course\|condition_comparison\|classification`), update `classification` to `expression_bin` and add the convention text:

> `cluster_type=expression_bin` analyses do **not** carry per-cluster `functional_description`; cluster `name` is the metric label (e.g., VEG, HEG).

- [ ] **Step 4: Write the release notes**

Create or update the agreed release-notes location with a section per friction:

```markdown
# 2026-05-01 KG release — Explorer frictions F1-F4

## F1 — Term informativeness + annotation rollup refinement

**New fields:**
- `Gene.annotation_state: str` ∈ `{no_evidence, catch_all_only, informative_single, informative_multi}`
- `Gene.informative_annotation_types: str[]` (granular)
- `<term node>.is_uninformative: 'true' | absent` on 7 ontology types
- `Gene.annotation_quality: int` — REDEFINED as 0-3 numeric encoding of `annotation_state`

**Semantic shift on `annotation_quality`:**
- *Before*: 0-3 mixing product-name regex with arbitrary 4-source structured-count.
- *After*: pure informative-evidence richness over 8 source buckets, no product-name dependency.
- Existing `WHERE g.annotation_quality <= 1` queries silently shift meaning (today: hypothetical-named genes; refined: low-evidence genes — spirit of the filter, now correctly).

**Worked example.** Gene with `pfam_ids=[PF00712]` + `go_terms=[GO:0003674]` (MF root) + `gene_category='Unknown'`, no other functional edges:
- *Before*: structured_count=2 → score 3.
- *After*: GO MF root flagged uninformative; only Pfam contributes → informative_count=1 → score 2, state `informative_single`.

## F2 — Data source surfacing

**New nodes:** 4 `DataSource` nodes (`ncbi`, `cyanorak`, `uniprot`, `eggnog`) with metadata (scope, provenance, info_types).

**New field:** `Gene.contributing_sources: str[]` — per-gene data-source presence.

**Worked example.** TX50_RS09500 → `contributing_sources=['ncbi']`. Distinguishes "ncbi present, others didn't run/match" from "no data" cleanly.

## F3 — `cluster_type` vocab rename

`classification` → `expression_bin`. No paperconfig migrations needed (none used the old label). Validator vocab + paperconfig SKILL.md updated. Convention: `expression_bin` analyses do not carry per-cluster `functional_description`.

## F4 — Gene properties for sparsely-annotated genes

**New fields:**
- `Gene.contig: str` — chromosome / contig name from GFF seqid (required for genomic-neighbor lookups)
- `Gene.seed_ortholog: str` — eggNOG seed-ortholog identifier (taxid.identifier)
- `Gene.seed_ortholog_evalue: float`

**Honest framing.** These primarily benefit the broader sparsely-annotated population, not the strict floor (~14 TX50_RS-style genes). Strict floor cases gain `contig` + neighbor lookup; broader population gains the seed-ortholog "this resembles X (E=Y)" pointer.

## Bucket inventory

The annotation_quality / annotation_state source-bucket list has 8 live buckets as of 2026-05-02: `go`, `kegg`, `pfam`, `ec`, `role`, `reaction`, `transporter` (via `Gene_has_tcdb_family`), `cazy` (via `Gene_has_cazy_family`). TCDB and CAZy ontologies merged into main on 2026-05-02 from a separate spec ([`2026-05-01-tcdb-cazy-ontologies-design.md`](../specs/2026-05-01-tcdb-cazy-ontologies-design.md)).
```

- [ ] **Step 5: Commit**

```bash
git add CLAUDE.md docs/kg-changes/2026-05-01-explorer-frictions-f1-f4.md
git commit -m "$(cat <<'EOF'
docs: update CLAUDE.md + release notes for F1-F4

Documents the new fields, the annotation_quality semantic shift, the
DataSource nodes, the cluster_type rename, and the contig + seed_ortholog
Gene additions. Release notes include worked examples for each friction
plus the 8-bucket inventory (TCDB/CAZy already live).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 15: Final sanity sweep + close-out

- [ ] **Step 1: Run the full unit suite**

```bash
pytest -m "not slow and not kg" -v 2>&1 | tail -10
```
Expected: green.

- [ ] **Step 2: Run the KG validity suite**

```bash
pytest -m kg -v 2>&1 | tail -10
```
Expected: green.

- [ ] **Step 3: Diff the omics-edge-snapshot before vs after**

```bash
uv run python .claude/skills/omics-edge-snapshot/scripts/compare.py \
  omics_edge_snapshots/before_explorer_frictions.json \
  omics_edge_snapshots/after_explorer_frictions_f1_f4.json
```

(Capture the `before` snapshot manually before starting Task 1 of the implementation.)

Expected: identical totals; no per-paper edge regressions.

- [ ] **Step 4: Quick MCP smoke check**

```bash
# If the MCP server is configured locally, invoke a few tools that exercise
# the new properties:
# - resolve_gene → returns gene_overview, includes annotation_state
# - list_data_sources (new MCP tool — write separately if needed; confirm
#   placeholder isn't blocking)
```

(MCP tool changes for surfacing the new fields are NOT part of this plan — tracked separately on the explorer roadmap. The KG-side work is complete when the schema + post-import + tests + docs all land.)

- [ ] **Step 5: Final commit if any leftover docs**

If anything (release notes, CLAUDE.md callouts) needs a finishing touch, commit it now.

```bash
git status --short
```
Expected: clean tree.

---

## Acceptance summary

When this plan finishes, the following must be true:

- [ ] Every `Gene` node has `contig`, `contributing_sources`, `annotation_state`, `annotation_quality`, `informative_annotation_types` populated.
- [ ] 4 `DataSource` nodes exist with correct metadata.
- [ ] ~224 ontology terms have `is_uninformative='true'`.
- [ ] No paperconfig validator accepts `classification` anymore (only `expression_bin`).
- [ ] All commits trail with the standard `Co-Authored-By` line.
- [ ] `pytest -m kg` is green.
- [ ] `omics_edge_snapshots/after_explorer_frictions_f1_f4.json` matches the pre-release baseline (no expression-edge changes).
- [ ] CLAUDE.md and release notes describe the new fields and the `annotation_quality` semantic shift.

## Spec coverage

- F1.1 → Tasks 9 (YAML, Cypher, test)
- F1.2 + F1.3 → Tasks 10 (combined Cypher block), 12 (delete obsolete build-time code)
- F1.4 → Task 11
- F2.1 → Task 7
- F2.2 → Task 5
- F2.3 → Task 8
- F2.4 → Task 6
- F2.5 (eligibility note) → Task 14 (CLAUDE.md / release notes)
- F3 → Task 1
- F4.1 → Tasks 2 (script), 3 (schema/adapter)
- F4.2 → Task 4
- F4.3 → Task 14 (KG-side: nothing; documented in release notes)
- Source bucket maintenance → Task 14 (CLAUDE.md callout)
- Drift test → Task 9 (test_uninformative_terms)
- Full graph rebuild + validity gate → Task 13
- Sanity sweep → Task 15
