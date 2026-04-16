# Reference Proteome Match Organisms — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `organism_type` property to all OrganismTaxon nodes and convert HP15 + AltMedDE from genome strains to reference proteome match organisms for Moreno 2023.

**Architecture:** Three new columns in `cyanobacteria_genomes.csv` (`organism_type`, `reference_database`, `reference_proteome`) flow through the CyanorakNcbi adapter into OrganismTaxon node properties. Treatment organisms get `organism_type: treatment` from the adapter. Moreno paperconfig organism names updated to match new `preferred_name` values. Schema config, tests, and docs updated.

**Tech Stack:** Python, YAML, CSV, Cypher, pytest

---

## File Map

| File | Action | Responsibility |
|---|---|---|
| `config/schema_config.yaml` | Modify | Add `organism_type`, `reference_database`, `reference_proteome` properties |
| `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv` | Modify | Add 3 columns; update HP15 + AltMedDE rows |
| `multiomics_kg/adapters/cyanorak_ncbi_adapter.py` | Modify | Read new CSV columns, set `organism_type` on genome + treatment nodes |
| `multiomics_kg/utils/gene_id_utils.py` | Modify | Update `ORGANISM_TO_GENOME_DIR` for renamed organisms |
| `data/Prochlorococcus/papers_and_supp/moreno 2023/paperconfig.yaml` | Modify | Update organism names + add assembly accessions to experiments |
| `scripts/validate_paperconfig.py` | Modify | Update `VALID_ORGANISMS` set with new names |
| `tests/test_cyanorak_ncbi_adapter.py` | Modify | Test `organism_type` on genome + treatment nodes |
| `tests/kg_validity/test_organism.py` | Modify | Update counts, names, add `organism_type` test |
| `docs/kg-changes/reference-proteome-match-organisms.md` | Create | Downstream what-changed doc |
| `docs/community_proteomics_marref_saga.md` | Modify | Reference spec, mark open questions answered |
| `scripts/post-import.sh` | Modify (if needed) | Index on `organism_type` |
| `scripts/post-import.cypher` | Modify (if needed) | Same |
| `cache/data/Alteromonas/genomes/Alt_MarRef/` | Delete + re-download | Clear old GCF_000020585.3 data, download GCA_003513035.1 |

---

### Task 1: Schema config — add organism_type properties

**Files:**
- Modify: `config/schema_config.yaml:311-335`

- [ ] **Step 1: Add properties to organism taxon schema**

In `config/schema_config.yaml`, add three properties to the `organism taxon` node after the existing `omics_types` property:

```yaml
organism taxon:
  represented_as: node
  preferred_id: insdc.gcf
  label_in_input: organism
  properties:
    organism_name: str
    preferred_name: str  # human-readable name, e.g. "Prochlorococcus MED4", "Alteromonas macleodii HOT1A3"
    strain_name: str
    ncbi_taxon_id: int
    clade: str           # research clade, e.g. HLI, HLII, LLIV (Prochlorococcus-specific)
    lineage: str         # full NCBI lineage string
    superkingdom: str    # e.g. Bacteria
    kingdom: str         # e.g. Bacillati
    phylum: str          # e.g. Cyanobacteriota
    tax_class: str       # e.g. Cyanophyceae
    order: str           # e.g. Synechococcales
    family: str          # e.g. Prochlorococcaceae
    genus: str           # e.g. Prochlorococcus
    species: str         # e.g. Prochlorococcus marinus
    organism_type: str   # genome_strain | treatment | reference_proteome_match
    reference_database: str   # e.g. "MarRef v6" (reference_proteome_match only)
    reference_proteome: str   # e.g. "UP000262181" (reference_proteome_match only)
    gene_count: int          # pre-computed: count of Gene_belongs_to_organism edges
    publication_count: int   # pre-computed: count of publications mentioning this organism
    experiment_count: int    # pre-computed: sum of experiment_count across matched publications
    treatment_types: str[]   # pre-computed: distinct treatment types across matched publications
    background_factors: str[] # pre-computed: distinct background_factors across matched publications
    omics_types: str[]       # pre-computed: distinct omics types across matched publications
```

- [ ] **Step 2: Commit**

```bash
git add config/schema_config.yaml
git commit -m "schema: add organism_type, reference_database, reference_proteome to OrganismTaxon"
```

---

### Task 2: Genomes CSV — add columns and update HP15/AltMedDE

**Files:**
- Modify: `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv`

- [ ] **Step 1: Add three new columns to the CSV header and all rows**

The header becomes:
```
ncbi_accession,cyanorak_organism,ncbi_taxon_id,strain_name,data_dir,clade,preferred_name,organism_type,reference_database,reference_proteome
```

All existing rows get `genome_strain,,` appended (organism_type=genome_strain, empty reference_database, empty reference_proteome).

- [ ] **Step 2: Update the HP15 row**

Change from:
```
GCF_000166295.1,,225937,HP15,cache/data/Marinobacter/genomes/HP15/,,Marinobacter adhaerens DSM 23420 / HP15
```

To:
```
GCF_000166295.1,,225937,HP15,cache/data/Marinobacter/genomes/HP15/,,Marinobacter (MarRef v6),reference_proteome_match,MarRef v6,GCF_000166295.1
```

- [ ] **Step 3: Update the AltMedDE row**

Change from:
```
GCF_000020585.3,,1774373,AltMedDE,cache/data/Alteromonas/genomes/AltMedDE/,,Alteromonas mediterranea DE
```

To:
```
GCA_003513035.1,,232,Alt_MarRef,cache/data/Alteromonas/genomes/Alt_MarRef/,,Alteromonas (MarRef v6),reference_proteome_match,MarRef v6,UP000262181
```

Note: three fields change — `ncbi_accession` (`GCF_000020585.3` → `GCA_003513035.1`), `ncbi_taxon_id` (`1774373` → `232`), `strain_name` (`AltMedDE` → `Alt_MarRef`), `data_dir` (`cache/data/Alteromonas/genomes/AltMedDE/` → `cache/data/Alteromonas/genomes/Alt_MarRef/`), and `preferred_name` (`Alteromonas mediterranea DE` → `Alteromonas (MarRef v6)`).

- [ ] **Step 4: Verify CSV parses correctly**

```bash
python -c "
import csv
with open('data/Prochlorococcus/genomes/cyanobacteria_genomes.csv') as f:
    lines = [l for l in f if not l.strip().startswith('#')]
    reader = csv.DictReader(lines)
    for row in reader:
        print(f\"{row['strain_name']:12s} {row.get('organism_type', 'MISSING'):30s} {row.get('preferred_name', '')}\")
"
```

Expected: all rows show `genome_strain` except HP15 and Alt_MarRef which show `reference_proteome_match`.

- [ ] **Step 5: Commit**

```bash
git add data/Prochlorococcus/genomes/cyanobacteria_genomes.csv
git commit -m "data: add organism_type columns to genomes CSV, convert HP15+AltMedDE to reference_proteome_match"
```

---

### Task 3: Clear AltMedDE cache and re-download with correct assembly

The AltMedDE cache directory currently has data from the wrong assembly (GCF_000020585.3 with MADE_RS locus tags). After changing the assembly to GCA_003513035.1 in the genomes CSV, clear the cache and re-download.

**Verified**: GCA_003513035.1 has PGAP annotations (4,305 protein-coding genes). The NCBI Datasets API and `_ncbi_download_genome` work with both GCF and GCA accessions — no code changes needed. The "also download GCA GFF" block in `download_genome_data.py` is harmlessly skipped when the primary accession is already GCA. All adapters use `insdc.gcf:` prefix consistently regardless of GCF/GCA.

**Files:**
- Delete: `cache/data/Alteromonas/genomes/AltMedDE/` (old folder, wrong assembly)

- [ ] **Step 1: Delete old cache data**

```bash
rm -rf cache/data/Alteromonas/genomes/AltMedDE/  # old folder from wrong assembly
```

- [ ] **Step 2: Re-run genome download for AltMedDE**

```bash
bash scripts/prepare_data.sh --strains Alt_MarRef --steps 0 --force
```

This downloads the GCA_003513035.1 assembly's GFF, protein FASTA, and GBFF files to `cache/data/Alteromonas/genomes/Alt_MarRef/`. Verify the locus tags are `DEH24_*`:

```bash
grep -m5 'locus_tag=' cache/data/Alteromonas/genomes/Alt_MarRef/genomic.gff | head -5
```

Expected: locus tags should show `DEH24_NNNNN` pattern.

- [ ] **Step 3: Run eggNOG for Alt_MarRef in the background**

`prepare_data.sh --steps 0` intentionally skips eggNOG (sub-step 4). Use the `/eggnog-run` skill to run it separately once the download (step 2 above) has completed and `protein.faa` exists:

```bash
# Run in background — takes ~15-30 min per strain
/eggnog-run Alt_MarRef
```

This is non-blocking — the rest of the plan can proceed without eggNOG annotations. The old AltMedDE eggNOG annotations were for the wrong assembly's proteins and must be regenerated.

**Important**: Steps 1 2 3 4 of prepare_data (gene_mapping, annotations, ID mapping, paper resolution) depend on the paperconfig and ORGANISM_TO_GENOME_DIR being updated first. Those are done in Tasks 5 and 6, then Task 7 runs the build steps.

- [ ] **Step 4: Commit**

No code changes to commit — this is a cache operation only.

---

### Task 4: CyanorakNcbi adapter — read and emit organism_type properties

**Files:**
- Modify: `multiomics_kg/adapters/cyanorak_ncbi_adapter.py:160-175` (CyanorakNcbi.__init__)
- Modify: `multiomics_kg/adapters/cyanorak_ncbi_adapter.py:290-325` (get_organism_node)
- Modify: `multiomics_kg/adapters/cyanorak_ncbi_adapter.py:460-477` (MultiCyanorakNcbi.__init__ CSV reader)
- Modify: `multiomics_kg/adapters/cyanorak_ncbi_adapter.py:494-508` (_get_treatment_organism_nodes)
- Test: `tests/test_cyanorak_ncbi_adapter.py`

- [ ] **Step 1: Write failing test for organism_type on genome organism**

Add to `tests/test_cyanorak_ncbi_adapter.py`:

```python
def test_organism_node_has_organism_type(med4_adapter):
    """Genome strain organisms should have organism_type='genome_strain'."""
    nodes = med4_adapter.get_nodes()
    org_nodes = [n for n in nodes if n[1] == "organism"]
    assert len(org_nodes) == 1
    props = org_nodes[0][2]
    assert props.get("organism_type") == "genome_strain"
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tests/test_cyanorak_ncbi_adapter.py::test_organism_node_has_organism_type -v
```

Expected: FAIL — `organism_type` not in properties.

- [ ] **Step 3: Add organism_type parameter to CyanorakNcbi.__init__**

In `cyanorak_ncbi_adapter.py`, add `organism_type`, `reference_database`, and `reference_proteome` to `__init__`:

```python
def __init__(
    self,
    ncbi_accession: str = None,
    data_dir: str = None,
    cyanorak_organism: str = None,
    strain_name: str = None,
    ncbi_taxon_id: int = None,
    clade: str = None,
    preferred_name: str = None,
    organism_type: str = "genome_strain",
    reference_database: str = None,
    reference_proteome: str = None,
    **kwargs,
):
```

Store them:
```python
self.organism_type = organism_type
self.reference_database = reference_database
self.reference_proteome = reference_proteome
```

- [ ] **Step 4: Emit organism_type in get_organism_node**

In the organism node properties dict (around line 300), add after the clade block:

```python
properties['organism_type'] = self.organism_type
if self.reference_database:
    properties['reference_database'] = self.reference_database
if self.reference_proteome:
    properties['reference_proteome'] = self.reference_proteome
```

- [ ] **Step 5: Pass new columns from MultiCyanorakNcbi CSV reader**

In `MultiCyanorakNcbi.__init__`, update the adapter construction (around line 468):

```python
adapter = CyanorakNcbi(
    ncbi_accession=row['ncbi_accession'],
    data_dir=row.get('data_dir') or None,
    cyanorak_organism=row.get('cyanorak_organism') or None,
    strain_name=row.get('strain_name') or None,
    ncbi_taxon_id=ncbi_taxon_id,
    clade=clade,
    preferred_name=row.get('preferred_name') or None,
    organism_type=row.get('organism_type') or 'genome_strain',
    reference_database=row.get('reference_database') or None,
    reference_proteome=row.get('reference_proteome') or None,
    **kwargs,
)
```

- [ ] **Step 6: Set organism_type on treatment organism nodes**

In `_get_treatment_organism_nodes`, add to the props dict:

```python
props = {
    'organism_name': organism_name,
    'preferred_name': organism_name,
    'ncbi_taxon_id': taxid,
    'organism_type': 'treatment',
}
```

- [ ] **Step 7: Run test to verify it passes**

```bash
pytest tests/test_cyanorak_ncbi_adapter.py::test_organism_node_has_organism_type -v
```

Expected: PASS

- [ ] **Step 8: Run full adapter test suite**

```bash
pytest tests/test_cyanorak_ncbi_adapter.py -v
```

Expected: all tests pass.

- [ ] **Step 9: Commit**

```bash
git add multiomics_kg/adapters/cyanorak_ncbi_adapter.py tests/test_cyanorak_ncbi_adapter.py
git commit -m "feat: emit organism_type on all OrganismTaxon nodes"
```

---

### Task 5: Update ORGANISM_TO_GENOME_DIR for renamed organisms

**Files:**
- Modify: `multiomics_kg/utils/gene_id_utils.py:63-72`

- [ ] **Step 1: Update the lookup entries for HP15 and AltMedDE**

In `gene_id_utils.py`, replace the HP15 and AltMedDE entries in `ORGANISM_TO_GENOME_DIR`:

```python
# CSV-ready papers batch (April 2026): Dominguez 2017, Fuszard 2012, Moreno 2023
"prochlorococcus marinus subsp. marinus ccmp1375 (ss120)": "cache/data/Prochlorococcus/genomes/SS120",
"prochlorococcus marinus ss120": "cache/data/Prochlorococcus/genomes/SS120",
"prochlorococcus ss120": "cache/data/Prochlorococcus/genomes/SS120",
"synechococcus sp. bl107": "cache/data/Synechococcus/genomes/BL107",
"synechococcus bl107": "cache/data/Synechococcus/genomes/BL107",
# Reference proteome match organisms (Moreno 2023 community fractions)
"marinobacter (marref v6)": "cache/data/Marinobacter/genomes/HP15",
"marinobacter adhaerens dsm 23420 / hp15": "cache/data/Marinobacter/genomes/HP15",
"marinobacter adhaerens hp15": "cache/data/Marinobacter/genomes/HP15",
"marinobacter hp15": "cache/data/Marinobacter/genomes/HP15",
"alteromonas (marref v6)": "cache/data/Alteromonas/genomes/Alt_MarRef",
"alteromonas mediterranea de": "cache/data/Alteromonas/genomes/Alt_MarRef",  # backward compat alias
```

Keep the old names as aliases — paperconfigs may still reference them until Task 5 updates them, and other code may use the old names.

- [ ] **Step 2: Commit**

```bash
git add multiomics_kg/utils/gene_id_utils.py
git commit -m "feat: add new organism names to ORGANISM_TO_GENOME_DIR for reference proteome match organisms"
```

---

### Task 6: Update Moreno 2023 paperconfig — organism names + assembly accessions

**Files:**
- Modify: `data/Prochlorococcus/papers_and_supp/moreno 2023/paperconfig.yaml`

- [ ] **Step 1: Replace all Marinobacter organism references**

Find all instances of:
```yaml
organism: "Marinobacter adhaerens DSM 23420 / HP15"
```

Replace with:
```yaml
organism: "Marinobacter (MarRef v6)"
```

There are ~19 occurrences across S4 experiments and supplementary table entries.

- [ ] **Step 2: Replace all Alteromonas organism references**

Find all instances of:
```yaml
organism: "Alteromonas mediterranea DE"
```

Replace with:
```yaml
organism: "Alteromonas (MarRef v6)"
```

There are ~20 occurrences across S3 experiments and supplementary table entries.

- [ ] **Step 3: Add assembly accession to Marinobacter and Alteromonas experiments**

For all Marinobacter S4 experiments, add `assembly_accession` to the experiment definition (not just the `treatment_assembly_accession` which is the coculture partner). This documents which reference proteome the data maps to:

```yaml
    light_low_glucose_marino_in_med4_proteomics:
      name: "Marinobacter in MED4 coculture: light +100 nM glucose vs light (control)"
      organism: "Marinobacter (MarRef v6)"
      assembly_accession: "GCF_000166295.1"
      ...
```

For all Alteromonas S3 experiments:

```yaml
    light_low_glucose_alt_in_med4_proteomics:
      name: "Alteromonas in MED4 coculture: light +100 nM glucose vs light (control)"
      organism: "Alteromonas (MarRef v6)"
      assembly_accession: "GCA_003513035.1"
      ...
```

Also add `assembly_accession` to the `csv` supplementary table entries for S3 and S4.

Note: `assembly_accession` is a new optional field. Verify that the omics adapter and paperconfig validator accept it without errors (it should be ignored by existing code since it's not a recognized field — check whether the validator rejects unknown fields).

- [ ] **Step 4: Verify YAML parses**

```bash
python -c "
import yaml
with open('data/Prochlorococcus/papers_and_supp/moreno 2023/paperconfig.yaml') as f:
    config = yaml.safe_load(f)
exps = config['publication']['experiments']
organisms = set(e.get('organism', '') for e in exps.values())
print('Organisms in experiments:', sorted(organisms))
"
```

Expected: should include `Marinobacter (MarRef v6)` and `Alteromonas (MarRef v6)` instead of old names.

- [ ] **Step 5: Commit**

```bash
git add "data/Prochlorococcus/papers_and_supp/moreno 2023/paperconfig.yaml"
git commit -m "data: update Moreno 2023 paperconfig organism names to MarRef v6 convention"
```

---

### Task 7: Update validate_paperconfig.py organism list

**Files:**
- Modify: `scripts/validate_paperconfig.py:134-135`

- [ ] **Step 1: Update VALID_ORGANISMS set**

Replace old names with new names:

```python
    "Alteromonas (MarRef v6)",
    "Marinobacter (MarRef v6)",
```

- [ ] **Step 2: Commit**

```bash
git add scripts/validate_paperconfig.py
git commit -m "fix: update VALID_ORGANISMS with renamed reference proteome match organisms"
```

---

### Task 8: Build gene_mapping + annotations + resolve paper IDs for Alt_MarRef

**Depends on**: Task 3 (download complete), Task 5 (ORGANISM_TO_GENOME_DIR updated), Task 6 (paperconfig updated)

Now that the new assembly is downloaded and the paperconfig + ORGANISM_TO_GENOME_DIR point to the right paths and names, build the annotation pipeline for Alt_MarRef.

- [ ] **Step 1: Run prepare_data steps 1-4 for Alt_MarRef**

```bash
bash scripts/prepare_data.sh --strains Alt_MarRef --steps 1 2 3 4 --force
```

- [ ] **Step 2: Verify gene_mapping has DEH24 locus tags**

```bash
head -5 cache/data/Alteromonas/genomes/Alt_MarRef/gene_mapping.csv
```

Expected: locus_tag column shows `DEH24_*` values.

- [ ] **Step 3: Verify Moreno S3 paper CSVs resolve**

```bash
ls cache/data/Alteromonas/genomes/Alt_MarRef/gene_id_mapping.json
# Check that resolve_paper_ids created _resolved.csv files for Moreno S3 tables
find "data/Prochlorococcus/papers_and_supp/moreno 2023/" -name "*s3*resolved*" | head -5
```

Expected: `gene_id_mapping.json` exists, and `_resolved.csv` files exist for Moreno S3 supplementary tables.

---

### Task 9: Update KG validity tests

**Files:**
- Modify: `tests/kg_validity/test_organism.py:50-88`

- [ ] **Step 1: Update organism count docstring and Alteromonas species test**

In `test_organism_count`, update the docstring: HP15 is now `Marinobacter (MarRef v6)`, AltMedDE is now `Alteromonas (MarRef v6)`. Total count stays 32.

```python
def test_organism_count(run_query):
    """32 OrganismTaxon nodes: 25 genome strains + 5 treatment organisms + 2 reference proteome match.

    Genome strains (25): 10 Pro (MED4, AS9601, MIT9301, MIT9312, MIT9313,
    MIT9303, NATL1A, NATL2A, RSP50, SS120), 6 Syn (CC9311, WH7803, WH8102,
    BL107, PCC7002, PCC7942, UTEX2973), 1 Thermosynechococcus (BP1),
    3 Alteromonas (MIT1002, EZ55, HOT1A3), 4 heterotrophs
    (W3-18-1, KT2440, DSS-3, MruberA).
    Reference proteome match (2): Marinobacter (MarRef v6), Alteromonas (MarRef v6).
    Treatment organisms (5): Phage, Alteromonas (genus), Vibrio
    parahaemolyticus, Meiothermus ruber, E. coli.
    """
    result = run_query("MATCH (o:OrganismTaxon) RETURN count(o) AS cnt")
    assert result[0]["cnt"] == 32, (
        f"Expected 32 OrganismTaxon nodes, got {result[0]['cnt']}"
    )
```

- [ ] **Step 2: Update Alteromonas species test**

The `test_alteromonas_strains_have_species` test expects 4 Alteromonas strains with `strain_name IS NOT NULL`. AltMedDE (now `Alteromonas (MarRef v6)`) still has `strain_name=AltMedDE`, so the count stays 4. But verify the species derived from `Alteromonas (MarRef v6)` — the species derivation logic in the adapter splits `preferred_name` into genus + lowercase second word. `"Alteromonas (MarRef v6)"` has `(MarRef` as second word which is not lowercase, so species won't be derived. We need to check if this test will fail.

The test asserts `species in {"Alteromonas macleodii", "Alteromonas mediterranea"}`. For `Alteromonas (MarRef v6)`, species derivation will fail (second word `(MarRef` starts with `(`), so `species` will be whatever NCBI taxonomy returns for taxid 232 (genus-level *Alteromonas*). Species will be null and the test will fail.

Update the test to handle reference proteome match organisms that may not have species:

```python
def test_alteromonas_strains_have_species(run_query):
    """Alteromonas genome_strain nodes must have a species property.

    Three strains (MIT1002, EZ55, HOT1A3) are A. macleodii.
    Reference proteome match organisms (Alteromonas MarRef v6) are excluded
    as their species may not be derivable from the preferred_name.
    """
    result = run_query("""
        MATCH (o:OrganismTaxon)
        WHERE o.genus = 'Alteromonas'
          AND o.strain_name IS NOT NULL
          AND o.organism_type = 'genome_strain'
        RETURN o.preferred_name AS name, o.species AS species
    """)
    assert len(result) == 3, f"Expected 3 Alteromonas genome strains, got {len(result)}"
    VALID_SPECIES = {"Alteromonas macleodii"}
    for r in result:
        assert r["species"] in VALID_SPECIES, (
            f"{r['name']} has species={r['species']!r}, expected one of {VALID_SPECIES}"
        )
```

- [ ] **Step 3: Add test for organism_type coverage**

```python
def test_organism_type_on_all_nodes(run_query):
    """Every OrganismTaxon node must have an organism_type property."""
    result = run_query("""
        MATCH (o:OrganismTaxon)
        WHERE o.organism_type IS NULL
        RETURN o.preferred_name AS name
    """)
    assert len(result) == 0, (
        f"OrganismTaxon nodes missing organism_type: {[r['name'] for r in result]}"
    )


def test_organism_type_values(run_query):
    """organism_type must be one of the three valid values."""
    result = run_query("""
        MATCH (o:OrganismTaxon)
        RETURN o.organism_type AS otype, count(o) AS cnt
        ORDER BY cnt DESC
    """)
    counts = {r["otype"]: r["cnt"] for r in result}
    assert set(counts.keys()) == {"genome_strain", "treatment", "reference_proteome_match"}, (
        f"Unexpected organism_type values: {counts}"
    )
    assert counts["genome_strain"] == 25
    assert counts["treatment"] == 5
    assert counts["reference_proteome_match"] == 2


def test_reference_proteome_match_properties(run_query):
    """Reference proteome match organisms must have reference_database and reference_proteome."""
    result = run_query("""
        MATCH (o:OrganismTaxon)
        WHERE o.organism_type = 'reference_proteome_match'
        RETURN o.preferred_name AS name,
               o.reference_database AS db,
               o.reference_proteome AS proteome
    """)
    assert len(result) == 2
    for r in result:
        assert r["db"] is not None, f"{r['name']} missing reference_database"
        assert r["proteome"] is not None, f"{r['name']} missing reference_proteome"
        assert "MarRef" in r["db"], f"{r['name']} unexpected reference_database: {r['db']}"
```

- [ ] **Step 4: Commit**

```bash
git add tests/kg_validity/test_organism.py
git commit -m "test: update KG validity tests for organism_type and renamed organisms"
```

---

### Task 10: Post-import index for organism_type (optional)

**Files:**
- Modify: `scripts/post-import.sh`
- Modify: `scripts/post-import.cypher`

- [ ] **Step 1: Add index on organism_type**

In both `scripts/post-import.sh` and `scripts/post-import.cypher`, add after the existing organism-related indexes:

```cypher
CREATE INDEX organism_type_idx IF NOT EXISTS FOR (o:OrganismTaxon) ON (o.organism_type);
```

- [ ] **Step 2: Commit**

```bash
git add scripts/post-import.sh scripts/post-import.cypher
git commit -m "post-import: add index on OrganismTaxon.organism_type"
```

---

### Task 11: Documentation

**Files:**
- Create: `docs/kg-changes/reference-proteome-match-organisms.md`
- Modify: `docs/community_proteomics_marref_saga.md`

- [ ] **Step 1: Create what-changed doc**

Create `docs/kg-changes/reference-proteome-match-organisms.md`:

```markdown
# Reference Proteome Match Organisms

**Date**: 2026-04-16
**Spec**: `docs/superpowers/specs/2026-04-16-reference-proteome-match-organisms-design.md`

## Changes

### New `organism_type` property on OrganismTaxon

All OrganismTaxon nodes now carry `organism_type`:
- `genome_strain` (25 organisms) — real genome assembly from a cultured isolate
- `treatment` (5 organisms) — non-genomic organisms used as coculture partners
- `reference_proteome_match` (2 organisms) — identified by matching experimental data against a multi-organism reference database

### New properties (reference_proteome_match only)

- `reference_database`: the database used for matching (e.g., "MarRef v6")
- `reference_proteome`: accession of the matched reference proteome

### Renamed organisms

| Old name | New name | organism_type |
|---|---|---|
| Marinobacter adhaerens DSM 23420 / HP15 | Marinobacter (MarRef v6) | reference_proteome_match |
| Alteromonas mediterranea DE | Alteromonas (MarRef v6) | reference_proteome_match |

### Assembly fix

Alteromonas entry changed from GCF_000020585.3 (wrong — A. mediterranea DE RefSeq) to GCA_003513035.1 (correct — the actual MarRef-matched reference proteome with DEH24_* locus tags).

### Query examples

```cypher
-- Filter to genome strains only (exclude community fractions)
MATCH (o:OrganismTaxon)
WHERE o.organism_type = 'genome_strain'
RETURN o.preferred_name

-- Find all reference proteome match organisms
MATCH (o:OrganismTaxon)
WHERE o.organism_type = 'reference_proteome_match'
RETURN o.preferred_name, o.reference_database, o.reference_proteome

-- Exclude community-fraction expression data
MATCH (e:Experiment)-[:Changes_expression_of]->(g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
WHERE o.organism_type <> 'reference_proteome_match'
RETURN e.name, g.locus_tag
```

### No changes to

- Edge types or edge properties
- Gene ID format (`ncbigene:<locus_tag>`)
- Experiment node properties
- Post-import computed properties (gene_count, etc.)
```

- [ ] **Step 2: Update community_proteomics_marref_saga.md**

Add at the top of the "Open work" section:

```markdown
> **Update 2026-04-16:** Questions 1-3 below have been answered by the reference proteome match organisms spec (`docs/superpowers/specs/2026-04-16-reference-proteome-match-organisms-design.md`). Summary: HP15 and AltMedDE converted to `organism_type: reference_proteome_match` with `reference_database: MarRef v6`. AltMedDE assembly fixed to GCA_003513035.1. See `docs/kg-changes/reference-proteome-match-organisms.md`.
```

- [ ] **Step 3: Update CLAUDE.md**

Update the following sections in `CLAUDE.md`:
- "Strains in graph" — change HP15 and AltMedDE descriptions to note they are reference proteome match organisms
- "Key graph facts" — update organism count breakdown: "25 genome strains + 5 treatment + 2 reference_proteome_match"
- OrganismTaxon properties list — add `organism_type`, `reference_database`, `reference_proteome`

- [ ] **Step 4: Commit**

```bash
git add docs/kg-changes/reference-proteome-match-organisms.md docs/community_proteomics_marref_saga.md CLAUDE.md
git commit -m "docs: add what-changed doc for reference proteome match organisms, update CLAUDE.md"
```

---

### Task 12: Run unit tests and verify

- [ ] **Step 1: Run full unit test suite**

```bash
pytest -m "not slow and not kg" -v
```

Expected: all tests pass. No regressions.

- [ ] **Step 2: Verify genomes CSV loads correctly end-to-end**

```bash
python -c "
from multiomics_kg.adapters.cyanorak_ncbi_adapter import MultiCyanorakNcbi
m = MultiCyanorakNcbi(
    'data/Prochlorococcus/genomes/cyanobacteria_genomes.csv',
    treatment_organisms_file='data/Prochlorococcus/treatment_organisms.csv'
)
print(f'Loaded {len(m.adapters)} genome adapters')
for a in m.adapters:
    print(f'  {a.strain_name:12s} organism_type={a.organism_type:30s} preferred_name={a.preferred_name}')
"
```

Expected: HP15 and Alt_MarRef show `organism_type=reference_proteome_match`, all others show `genome_strain`.
