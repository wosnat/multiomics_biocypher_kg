# Plan: Protein Annotation Preprocessing Improvements

## Goal

Improve `protein_annotations.json` output quality by properly splitting concatenated
strings into arrays and cleaning up evidence-tag noise.
Most changes go in `config/protein_annotations_config.yaml` and
`multiomics_kg/download/utils/annotation_transforms.py`; small plumbing
goes in `build_protein_annotations.py` (or `annotation_helpers.py`).

---

## 1. New field type support: `passthrough_list` with `split_pattern`

Currently `_resolve_passthrough_list` splits on a fixed `delimiter` character.
Several UniProt fields are concatenated with a repeating keyword prefix (e.g.
`CATALYTIC ACTIVITY: ... CATALYTIC ACTIVITY: ...`) that cannot be split by a
single delimiter character.

**Change in `build_protein_annotations.py` → `_resolve_passthrough_list`:**

```python
# existing:  tokens = _coerce_to_tokens(raw, delimiter)
# new:
split_pattern = fconf.get("split_pattern")
if split_pattern:
    import re
    tokens = [t.strip() for t in re.split(split_pattern, str(raw)) if t.strip()]
else:
    tokens = _coerce_to_tokens(raw, delimiter)
```

This is a backward-compatible change; fields without `split_pattern` continue
to work exactly as before.

---

## 2. New transforms in `annotation_transforms.py`

Six transforms are needed. All are pure functions on a single string token.

| Transform name | Input | Output | Used for |
|---|---|---|---|
| `clean_function_description` | full cc_function string | clean str | `function_description` |
| `clean_catalytic_activity` | one `CATALYTIC ACTIVITY: …` chunk | clean str | `catalytic_activities` |
| `extract_cofactor_name` | one `COFACTOR: Name=X; …` chunk | name str | `cofactor_names` |
| `extract_pathway_name` | one `PATHWAY: …` chunk | clean str | `pathways` |
| `extract_tm_range` | one `TRANSMEM N..M; …` chunk | `'N..M'` str | `transmembrane_regions` |
| `extract_signal_range` | full `SIGNAL N..M; …` string | `'N..M'` str | `signal_peptide` |

### Implementation sketches

```python
_ECO_PATTERN = re.compile(r'\s*\{ECO:[^}]+\}[.,]?\s*')

def _tx_clean_function_description(value: str) -> str:
    """Strip 'FUNCTION:' prefix and inline ECO evidence tags."""
    s = re.sub(r'^FUNCTION:\s*', '', value.strip(), flags=re.IGNORECASE)
    s = _ECO_PATTERN.sub(' ', s).strip().rstrip('.')
    return s

def _tx_clean_catalytic_activity(value: str) -> str:
    """Strip 'CATALYTIC ACTIVITY:' prefix and ECO tags from one reaction chunk."""
    s = re.sub(r'^CATALYTIC ACTIVITY:\s*', '', value.strip(), flags=re.IGNORECASE)
    s = _ECO_PATTERN.sub(' ', s).strip().rstrip(';').strip()
    return s

def _tx_extract_cofactor_name(value: str) -> str:
    """'COFACTOR: Name=FMN; Xref=…' → 'FMN'"""
    m = re.match(r'COFACTOR:\s*Name=([^;]+)', value.strip(), re.IGNORECASE)
    return m.group(1).strip() if m else ""

def _tx_extract_pathway_name(value: str) -> str:
    """'PATHWAY: Energy metabolism; oxidative phosphorylation. {ECO:…}.' → clean string"""
    s = re.sub(r'^PATHWAY:\s*', '', value.strip(), flags=re.IGNORECASE)
    s = _ECO_PATTERN.sub(' ', s).strip().rstrip('.')
    return s

def _tx_extract_tm_range(value: str) -> str:
    """'TRANSMEM 32..50; /note="Helical"; …' → '32..50'"""
    m = re.search(r'TRANSMEM\s+(\d+\.\.\d+)', value)
    return m.group(1) if m else ""

def _tx_extract_signal_range(value: str) -> str:
    """'SIGNAL 1..26; /evidence="…"' → '1..26'"""
    m = re.search(r'SIGNAL\s+(\d+\.\.\d+)', value)
    return m.group(1) if m else ""
```

Register all six in `_TRANSFORMS`.

---

## 3. Config changes in `protein_annotations_config.yaml`

### 3a. `keywords` → list  (simple config-only change)

```yaml
# before:
keywords:
  type: passthrough
  field: keyword

# after:
keywords:
  type: passthrough_list
  field: keyword
  delimiter: ";"
```

### 3b. `function_description` — strip ECO tags

Replace single `strip_function_prefix` transform with the new combined one:

```yaml
function_description:
  type: passthrough
  field: cc_function
  transform: clean_function_description   # was: strip_function_prefix
```

### 3c. `catalytic_activity` → `catalytic_activities` (list)

```yaml
# remove old:
# catalytic_activity:
#   type: passthrough
#   field: cc_catalytic_activity

# add:
catalytic_activities:
  type: passthrough_list
  field: cc_catalytic_activity
  split_pattern: "(?=CATALYTIC ACTIVITY:)"
  transform: clean_catalytic_activity
```

### 3d. `cofactors` → `cofactor_names` (list)

```yaml
# remove old:
# cofactors:
#   type: passthrough
#   field: cc_cofactor

# add:
cofactor_names:
  type: passthrough_list
  field: cc_cofactor
  split_pattern: "(?=COFACTOR:)"
  transform: extract_cofactor_name
```

### 3e. `pathway_description` → `pathways` (list)

```yaml
# remove old:
# pathway_description:
#   type: passthrough
#   field: cc_pathway

# add:
pathways:
  type: passthrough_list
  field: cc_pathway
  split_pattern: "(?=PATHWAY:)"
  transform: extract_pathway_name
```

### 3f. `transmembrane_regions` → list of range strings

```yaml
# before:
transmembrane_regions:
  type: passthrough
  field: ft_transmem

# after:
transmembrane_regions:
  type: passthrough_list
  field: ft_transmem
  split_pattern: "(?=TRANSMEM)"
  transform: extract_tm_range
```

### 3g. `signal_peptide` → range string

```yaml
# before:
signal_peptide:
  type: passthrough
  field: ft_signal

# after:
signal_peptide:
  type: passthrough
  field: ft_signal
  transform: extract_signal_range
```

---

## 4. Downstream updates (ripple from renames)

Three fields are renamed: `catalytic_activity` → `catalytic_activities`,
`cofactors` → `cofactor_names`, `pathway_description` → `pathways`.

### 4a. `config/schema_config.yaml`

Update the Protein node property declarations:

| Old name | New name | Old type | New type |
|---|---|---|---|
| `catalytic_activity` | `catalytic_activities` | `str` | `str[]` |
| `cofactors` | `cofactor_names` | `str` | `str[]` |
| `pathway_description` | `pathways` | `str` | `str[]` |
| `keywords` | `keywords` | `str` | `str[]` |
| `transmembrane_regions` | `transmembrane_regions` | `str` | `str[]` |

Gene node also has `catalytic_activity: str` at line 221 and 249 — update to
`catalytic_activities: str[]` in both places.

### 4b. `multiomics_kg/adapters/uniprot_adapter.py`

Lines 164–168 and 179 reference old field names. Update the `entry.get(...)` keys
to match the new names:

```python
"catalytic_activities":  entry.get("catalytic_activities"),
"cofactor_names":        entry.get("cofactor_names"),
"pathways":              entry.get("pathways"),
"transmembrane_regions": entry.get("transmembrane_regions"),   # same name, now a list
"signal_peptide":        entry.get("signal_peptide"),          # same name, now '1..26'
"keywords":              entry.get("keywords"),                 # same name, now a list
```

### 4c. `config/gene_annotations_config.yaml`

Lines 434–452 pull protein-level fields onto gene nodes as passthrough from the
uniprot source. Update field names to match:

```yaml
# line ~434
catalytic_activities:
  type: passthrough
  source: uniprot
  field: catalytic_activities

# transmembrane_regions and signal_peptide keep same name,
# but are now list / short string respectively — no YAML change needed there.
```

---

## 5. Files touched (summary)

| File | Change |
|---|---|
| `config/protein_annotations_config.yaml` | 7 field edits (types, transforms, renames) |
| `multiomics_kg/download/utils/annotation_transforms.py` | 6 new transforms + registry entries |
| `multiomics_kg/download/build_protein_annotations.py` | `split_pattern` support in `_resolve_passthrough_list` (~5 lines) |
| `config/schema_config.yaml` | 5 property type updates across Gene + Protein nodes |
| `multiomics_kg/adapters/uniprot_adapter.py` | 3 key renames in `entry.get(...)` calls |
| `config/gene_annotations_config.yaml` | 1 field rename (`catalytic_activity` → `catalytic_activities`) |

---

## 6. Implementation order

1. `annotation_transforms.py` — add 6 transforms (self-contained, easy to unit-test)
2. `test_build_protein_annotations.py` — add test classes 7a–7c (run immediately to confirm transforms work)
3. `build_protein_annotations.py` — add `split_pattern` plumbing
4. `protein_annotations_config.yaml` — apply all 7 field changes
5. Run tests: `pytest tests/test_build_protein_annotations.py -v` (all new tests should pass)
6. Rebuild cache: `uv run python multiomics_kg/download/build_protein_annotations.py --force --strains EZ55`
7. Verify output shape with quick Python inspection
8. `schema_config.yaml` — update type declarations
9. `uniprot_adapter.py` — update field key names
10. `gene_annotations_config.yaml` — update `catalytic_activity` reference
11. Full test run: `pytest -m "not slow and not kg" -v`

---

## 7. Tests in `tests/test_build_protein_annotations.py`

The existing file follows a clear pattern: one `class Test<Thing>` per unit.
Add the classes below; do not modify existing classes (they still pass).

### 7a. New transform unit tests

Add one class per new transform, mirroring the existing style.

```python
class TestTxCleanFunctionDescription:
    def test_strips_function_prefix(self):
        s = "FUNCTION: Catalyzes the reduction. {ECO:0000256,PIRNR:PIRNR000001}."
        assert _tx_clean_function_description(s) == "Catalyzes the reduction"

    def test_strips_eco_tags_only_no_prefix(self):
        s = "Catalyzes the reduction. {ECO:0000256,PIRNR:PIRNR000001}."
        assert _tx_clean_function_description(s) == "Catalyzes the reduction"

    def test_no_eco_passthrough(self):
        assert _tx_clean_function_description("Catalyzes the reduction.") \
               == "Catalyzes the reduction"

    def test_empty_returns_empty(self):
        assert _tx_clean_function_description("") == ""


class TestTxCleanCatalyticActivity:
    def test_strips_catalytic_activity_prefix(self):
        s = "CATALYTIC ACTIVITY: Reaction=ATP + H2O = ADP + phosphate; EC=3.6.1.3; Evidence={ECO:0000256,ARBA:ARB00001}."
        result = _tx_clean_catalytic_activity(s)
        assert "CATALYTIC ACTIVITY:" not in result
        assert "Reaction=ATP" in result

    def test_strips_eco_evidence(self):
        s = "CATALYTIC ACTIVITY: Reaction=A = B; Evidence={ECO:0000256,ARBA:ARBA00001}."
        assert "{ECO:" not in _tx_clean_catalytic_activity(s)

    def test_empty_chunk_returns_empty(self):
        assert _tx_clean_catalytic_activity("") == ""


class TestTxExtractCofactorName:
    def test_extracts_simple_name(self):
        s = "COFACTOR: Name=FMN; Xref=ChEBI:CHEBI:58210; Evidence={ECO:0000256|HAMAP-Rule:MF_02225}"
        assert _tx_extract_cofactor_name(s) == "FMN"

    def test_extracts_metal_ion_name(self):
        s = "COFACTOR: Name=Mg(2+); Xref=ChEBI:CHEBI:18420; Evidence={ECO:0000256}"
        assert _tx_extract_cofactor_name(s) == "Mg(2+)"

    def test_no_match_returns_empty(self):
        assert _tx_extract_cofactor_name("unrelated text") == ""

    def test_empty_returns_empty(self):
        assert _tx_extract_cofactor_name("") == ""


class TestTxExtractPathwayName:
    def test_strips_prefix_and_eco(self):
        s = "PATHWAY: Energy metabolism; oxidative phosphorylation. {ECO:0000256,ARBA:ARBA00004673}."
        result = _tx_extract_pathway_name(s)
        assert result == "Energy metabolism; oxidative phosphorylation"

    def test_no_eco_still_strips_prefix_and_period(self):
        s = "PATHWAY: Carbohydrate biosynthesis."
        assert _tx_extract_pathway_name(s) == "Carbohydrate biosynthesis"

    def test_empty_returns_empty(self):
        assert _tx_extract_pathway_name("") == ""


class TestTxExtractTmRange:
    def test_extracts_single_range(self):
        s = 'TRANSMEM 32..50; /note="Helical"; /evidence="ECO:0000256|SAM:Phobius"'
        assert _tx_extract_tm_range(s) == "32..50"

    def test_extracts_range_ignoring_trailing_annotation(self):
        s = 'TRANSMEM 7..24; /note="Helical"; /evidence="ECO:0000256|SAM:Phobius"; TRANSMEM 60..82'
        # Only the first TRANSMEM match is extracted (single-token use)
        assert _tx_extract_tm_range(s) == "7..24"

    def test_no_transmem_returns_empty(self):
        assert _tx_extract_tm_range("no TM here") == ""


class TestTxExtractSignalRange:
    def test_extracts_range(self):
        s = 'SIGNAL 1..26; /evidence="ECO:0000256|HAMAP-Rule:MF_01961"'
        assert _tx_extract_signal_range(s) == "1..26"

    def test_no_signal_returns_empty(self):
        assert _tx_extract_signal_range("TRANSMEM 1..20") == ""

    def test_empty_returns_empty(self):
        assert _tx_extract_signal_range("") == ""
```

### 7b. `split_pattern` support in `_resolve_passthrough_list`

Add a test class to `TestProteinAnnotationBuilderBuildMerged` (or a separate
`TestPassthroughListSplitPattern`) covering the new plumbing:

```python
class TestPassthroughListSplitPattern:
    """Tests for the split_pattern parameter in passthrough_list fields."""

    def _make_builder(self, field_cfg: dict) -> ProteinAnnotationBuilder:
        return ProteinAnnotationBuilder({"fields": {"out": field_cfg}})

    def test_single_entry_no_split_needed(self):
        builder = self._make_builder({
            "type": "passthrough_list",
            "field": "data",
            "split_pattern": r"(?=CATALYTIC ACTIVITY:)",
            "transform": "clean_catalytic_activity",
        })
        row = {"data": "CATALYTIC ACTIVITY: Reaction=A = B; EC=1.1.1.1; Evidence={ECO:0000256}."}
        result = builder.build_merged("X", row)
        assert isinstance(result["out"], list)
        assert len(result["out"]) == 1
        assert "Reaction=A = B" in result["out"][0]

    def test_multi_entry_split_into_list(self):
        builder = self._make_builder({
            "type": "passthrough_list",
            "field": "data",
            "split_pattern": r"(?=CATALYTIC ACTIVITY:)",
        })
        row = {"data": "CATALYTIC ACTIVITY: Reaction=A = B; CATALYTIC ACTIVITY: Reaction=C = D;"}
        result = builder.build_merged("X", row)
        assert len(result["out"]) == 2

    def test_empty_chunks_after_split_are_dropped(self):
        builder = self._make_builder({
            "type": "passthrough_list",
            "field": "data",
            "split_pattern": r"(?=COFACTOR:)",
            "transform": "extract_cofactor_name",
        })
        row = {"data": "COFACTOR: Name=FMN; Xref=ChEBI:CHEBI:58210; COFACTOR: Name=Mg(2+);"}
        result = builder.build_merged("X", row)
        assert result["out"] == ["FMN", "Mg(2+)"]

    def test_split_pattern_takes_priority_over_delimiter(self):
        """When split_pattern is set, delimiter is ignored."""
        builder = self._make_builder({
            "type": "passthrough_list",
            "field": "data",
            "split_pattern": r"(?=TRANSMEM)",
            "transform": "extract_tm_range",
            "delimiter": ";",   # should be ignored
        })
        row = {"data": 'TRANSMEM 12..31; /note="Helical"; TRANSMEM 37..57; /note="Helical"'}
        result = builder.build_merged("X", row)
        assert result["out"] == ["12..31", "37..57"]
```

### 7c. Integration tests via `build_merged` for new field configs

Add a class that tests the full field config for each changed field using
realistic raw UniProt strings:

```python
# Realistic raw strings matching actual uniprot_raw_data.json format
RAW_MULTI_CATALYTIC = (
    "CATALYTIC ACTIVITY: Reaction=ATP + H2O = ADP + phosphate; EC=3.6.1.3;"
    " Evidence={ECO:0000256,ARBA:ARBA00048988};"
    " CATALYTIC ACTIVITY: Reaction=C + D = E; EC=1.2.3.4;"
    " Evidence={ECO:0000256,HAMAP-Rule:MF_01920}"
)
RAW_MULTI_COFACTOR = (
    "COFACTOR: Name=FMN; Xref=ChEBI:CHEBI:58210; Evidence={ECO:0000256|HAMAP};"
    " Note=Binds 1 FMN.; COFACTOR: Name=Mg(2+); Xref=ChEBI:CHEBI:18420;"
    " Evidence={ECO:0000256|HAMAP}"
)
RAW_MULTI_PATHWAY = (
    "PATHWAY: Cofactor biosynthesis; CoA: step 2/5. {ECO:0000256|HAMAP}.;"
    " PATHWAY: Cofactor biosynthesis; CoA: step 3/5. {ECO:0000256|HAMAP}."
)
RAW_MULTI_TM = (
    'TRANSMEM 12..31; /note="Helical"; /evidence="ECO:0000256|SAM:Phobius";'
    ' TRANSMEM 37..57; /note="Helical"; /evidence="ECO:0000256|SAM:Phobius"'
)
RAW_SIGNAL = 'SIGNAL 1..26; /evidence="ECO:0000256|HAMAP-Rule:MF_01961"'
RAW_KEYWORDS = "ATP-binding;DNA repair;Helicase;Hydrolase"

NEW_FIELDS_CONFIG = {
    "fields": {
        "catalytic_activities": {
            "type": "passthrough_list",
            "field": "cc_catalytic_activity",
            "split_pattern": r"(?=CATALYTIC ACTIVITY:)",
            "transform": "clean_catalytic_activity",
        },
        "cofactor_names": {
            "type": "passthrough_list",
            "field": "cc_cofactor",
            "split_pattern": r"(?=COFACTOR:)",
            "transform": "extract_cofactor_name",
        },
        "pathways": {
            "type": "passthrough_list",
            "field": "cc_pathway",
            "split_pattern": r"(?=PATHWAY:)",
            "transform": "extract_pathway_name",
        },
        "transmembrane_regions": {
            "type": "passthrough_list",
            "field": "ft_transmem",
            "split_pattern": r"(?=TRANSMEM)",
            "transform": "extract_tm_range",
        },
        "signal_peptide": {
            "type": "passthrough",
            "field": "ft_signal",
            "transform": "extract_signal_range",
        },
        "keywords": {
            "type": "passthrough_list",
            "field": "keyword",
            "delimiter": ";",
        },
    }
}


class TestNewFieldConfigs:
    def setup_method(self):
        self.builder = ProteinAnnotationBuilder(NEW_FIELDS_CONFIG)

    def test_catalytic_activities_is_list(self):
        result = self.builder.build_merged("X", {"cc_catalytic_activity": RAW_MULTI_CATALYTIC})
        assert isinstance(result["catalytic_activities"], list)
        assert len(result["catalytic_activities"]) == 2

    def test_catalytic_activities_clean_no_prefix_or_eco(self):
        result = self.builder.build_merged("X", {"cc_catalytic_activity": RAW_MULTI_CATALYTIC})
        for item in result["catalytic_activities"]:
            assert "CATALYTIC ACTIVITY:" not in item
            assert "{ECO:" not in item

    def test_cofactor_names_extracts_names(self):
        result = self.builder.build_merged("X", {"cc_cofactor": RAW_MULTI_COFACTOR})
        assert result["cofactor_names"] == ["FMN", "Mg(2+)"]

    def test_pathways_is_clean_list(self):
        result = self.builder.build_merged("X", {"cc_pathway": RAW_MULTI_PATHWAY})
        pathways = result["pathways"]
        assert isinstance(pathways, list)
        assert len(pathways) == 2
        assert all("{ECO:" not in p for p in pathways)
        assert all("PATHWAY:" not in p for p in pathways)

    def test_transmembrane_regions_is_list_of_ranges(self):
        result = self.builder.build_merged("X", {"ft_transmem": RAW_MULTI_TM})
        assert result["transmembrane_regions"] == ["12..31", "37..57"]

    def test_signal_peptide_is_range_string(self):
        result = self.builder.build_merged("X", {"ft_signal": RAW_SIGNAL})
        assert result["signal_peptide"] == "1..26"

    def test_keywords_is_list(self):
        result = self.builder.build_merged("X", {"keyword": RAW_KEYWORDS})
        assert result["keywords"] == ["ATP-binding", "DNA repair", "Helicase", "Hydrolase"]

    def test_missing_fields_omitted_not_none(self):
        result = self.builder.build_merged("X", {})
        assert "catalytic_activities" not in result
        assert "cofactor_names" not in result
        assert "transmembrane_regions" not in result
```

### 7d. Imports to add at top of test file

```python
from multiomics_kg.download.utils.annotation_transforms import (
    _tx_add_go_prefix,
    _tx_clean_catalytic_activity,      # new
    _tx_clean_function_description,    # new
    _tx_extract_cofactor_name,         # new
    _tx_extract_go_from_brackets,
    _tx_extract_pathway_name,          # new
    _tx_extract_signal_range,          # new
    _tx_extract_tm_range,              # new
    _tx_strip_function_prefix,
)
```

### 7e. What existing tests already cover (no changes needed)

- `TestNonempty`, `TestSplit`, `TestCoerceToTokens` — unchanged helpers
- `TestLoadUniprotColumnar` — loader is unchanged
- `TestProteinAnnotationBuilderBuildMerged` — existing `passthrough_list` tests
  still verify baseline behavior; `split_pattern` is additive, not breaking
- `TestProcessTaxid` — skip/force logic is unchanged

---

## 8. What is NOT changing

- GO terms remain as three separate fields (`go_cellular_components`,
  `go_biological_processes`, `go_molecular_functions`) — already proper lists.
- `pfam_ids` split into `pfam_ids` / `pfam_names` is deferred to a later pass
  (affects gene annotations more than protein annotations).
- `protein_synonyms`, `caution_notes`, `interaction_notes`, `domain_description`,
  `subcellular_location`, `functional_motifs` — left as-is; no clear improvement
  without deeper UniProt parsing.
- `proteome_ids` split deferred — the `'UPxxx: Chromosome'` format has low
  downstream value for the KG.
