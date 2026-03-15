# Plan: `gene_category` Property on Gene Nodes

## Goal

Add a `gene_category` string property to every Gene node, providing a high-level functional classification (26 controlled categories). This enables filtering/grouping genes by biological function in the explorer UI.

---

## Data Sources & Coverage

| Source | Field | Coverage | Available for |
|--------|-------|----------|---------------|
| Cyanorak Role | `cyanorak_Role` (list, dotted codes like `D.1.1`) | 57–76% | 10 cyanobacteria strains (NOT RSP50, NOT Alteromonas) |
| TIGR Role | `tIGR_Role` (list, numeric IDs) + `tIGR_Role_description` | 73–94% | same 10 strains |
| COG Category | `cog_category` (list, single letters) | 68–94% | **all 13 strains** |

**Resolution priority** (first non-null wins):
1. `cyanorak_Role[0]` → map top-level letter (or second-level for "D")
2. `tIGR_Role_description[0]` → map MainRole (text before " / ")
3. `cog_category[0]` → map single letter
4. Fallback: `"Unknown"`

This gives ~90%+ coverage for all organisms.

---

## Cyanorak Role Distribution (from live data)

Cyanorak roles are available for 10 of 13 strains (NOT RSP50, NOT Alteromonas). 14,616 / 35,226 genes (41.5%) have at least one role.

### Top-level letters

| Letter | Genes | Description |
|--------|------:|-------------|
| R | 3,811 | Conserved hypothetical / Unknown |
| D | 3,186 | Diverse cellular processes |
| B | 1,725 | Cofactors, vitamins, pigments |
| J | 1,561 | Photosynthesis and respiration |
| G | 1,251 | Energy / carbon metabolism |
| K | 1,191 | Protein synthesis |
| Q | 975 | Transport |
| A | 900 | Amino acid biosynthesis |
| L | 785 | Protein fate |
| E | 695 | Central intermediary metabolism (N, P, S metabolism) |
| F | 635 | DNA metabolism |
| C | 569 | Cell envelope |
| M | 521 | Nucleotide metabolism |
| H | 453 | Lipid metabolism |
| P | 412 | Transcription |
| N | 320 | Regulatory functions |
| O | 196 | Signal transduction |
| I | 28 | Mobile elements |

### Category D sub-role breakdown (3,186 genes)

D is a catch-all for "Diverse cellular processes". D.1 ("Adaptation/acclimation") dominates at 82.8%.

| Sub-role | Genes | Description |
|----------|------:|-------------|
| **D.1** | **2,637** | **Adaptation and acclimation** |
| D.1.9 | 579 | Other |
| D.1.7 | 404 | Trace metals |
| D.1.5 | 362 | Phosphorus |
| D.1.2 | 340 | Light acclimation |
| D.1.1 | 307 | Iron stress |
| D.1.4 | 249 | Oxidative stress |
| D.1.3 | 179 | Nitrogen |
| D.1 (no sub) | 75 | Unclassified adaptation |
| **D.2** | **147** | **Cell division** |
| **D.4** | **157** | **Chaperones and heat shock** |
| **D.3** | **109** | **Detoxification / toxin resistance** |
| D.7 | 36 | Circadian clock |
| D.6 | 23 | Programmed cell death |
| D.5 | 19 | Chemotaxis and motility |

### Cyanorak E sub-role breakdown (695 genes)

Cyanorak E is "Central intermediary metabolism" — nutrient/element metabolism, NOT energy production.

| Sub-role | Genes | Description |
|----------|------:|-------------|
| E.4 | 139 | Nitrogen metabolism |
| E.3 | 66 | Phosphorus metabolism |
| E.7 | 57 | Sulfur metabolism |
| E.1 | 34 | Amino sugars, glycogen |
| E.6 | 26 | Polysaccharides/glycoproteins |
| E.2 | 15 | One-carbon metabolism |
| E.5 | 8 | Polyamine biosynthesis |
| E.8 | 5 | Other |

---

## Normalized Categories (26 values)

These are the output values for `gene_category`. Based on COG functional categories with additions for photosynthesis, stress response, and central intermediary metabolism (critical for cyanobacteria research).

| # | Category | COG | Cyanorak | TIGR |
|---|----------|-----|----------|------|
| 1 | Amino acid metabolism | E | A | Amino acid biosynthesis |
| 2 | Carbohydrate metabolism | G | G (partial) | — |
| 3 | Cell cycle and division | D, Z | D.2 | — |
| 4 | Cell motility | N | D.5 | — |
| 5 | Cell wall and membrane | M, W | C | Cell envelope |
| 6 | Cellular processes | — | D.3, D.6, D.7 | Cellular processes |
| 7 | Central intermediary metabolism | — | E | Central intermediary metabolism |
| 8 | Coenzyme metabolism | H | B | Biosynthesis of cofactors... |
| 9 | Defense mechanisms | V | — | — |
| 10 | Energy production | C | G (partial) | Energy metabolism |
| 11 | Inorganic ion transport | P | — | — |
| 12 | Intracellular trafficking | U | — | — |
| 13 | Lipid metabolism | I | H | Fatty acid and phospholipid metabolism |
| 14 | Mobile elements | X | I | Mobile and extrachromosomal... |
| 15 | Nucleotide metabolism | F | M | Purines, pyrimidines... |
| 16 | Photosynthesis | — | J | — |
| 17 | Post-translational modification | O | L, D.4 | Protein fate |
| 18 | Regulatory functions | — | N | Regulatory functions |
| 19 | Replication and repair | L, B | F | DNA metabolism |
| 20 | Secondary metabolites | Q | — | — |
| 21 | Signal transduction | T | O | Signal transduction |
| 22 | Stress response and adaptation | — | D.1 | — |
| 23 | Transcription | K, A | P | Transcription |
| 24 | Translation | J | K | Protein synthesis |
| 25 | Transport | — | Q | Transport and binding proteins |
| 26 | Unknown | S, R, Y | R | Hypothetical; Unknown; Unclassified |

---

## Category coverage asymmetry across organism groups

**6 categories are only reachable via Cyanorak/TIGR and have no COG mapping.** This means they will be 0% for Alteromonas (3 strains) and RSP50 (no Cyanorak/TIGR data). Genes with these functions in those organisms will be assigned a different (typically broader) COG-derived category instead.

| Category | Source | Cyanobacteria-specific biology? | What happens for Alteromonas/RSP50 |
|----------|--------|---------------------------------|------------------------------------|
| **Photosynthesis** | Cyanorak J | **Yes** — Alteromonas is heterotrophic | N/A — no photosynthesis genes |
| **Stress response and adaptation** | Cyanorak D.1 | No — universal | Stress genes get COG-derived labels (often "Unknown" or "Inorganic ion transport" for nutrient stress) |
| **Central intermediary metabolism** | Cyanorak E, TIGR | No — universal | N/P/S metabolism genes land in COG E ("Amino acid metabolism") or P ("Inorganic ion transport") |
| **Transport** | Cyanorak Q, TIGR | No — universal | Transport genes land in COG E/P/G (which bundle transport + metabolism together) |
| **Regulatory functions** | Cyanorak N, TIGR | No — universal | Regulatory genes get COG K ("Transcription") or T ("Signal transduction") |
| **Cellular processes** | Cyanorak D.3/D.6/D.7, TIGR | No — universal | Small category; genes get various COG labels |

**Why this happens:** COG categories bundle transport + metabolism together (e.g., COG E = "Amino acid transport **and** metabolism", COG P = "Inorganic ion transport **and** metabolism"), while Cyanorak and TIGR separate them. There is no COG letter for generic "Transport" or "Regulatory functions".

**Impact:** For the 10 cyanobacteria strains with Cyanorak data, `gene_category` provides fine-grained functional classification. For Alteromonas and RSP50 (COG-only), the classification is coarser — functionally equivalent genes may get different category labels than their cyanobacteria orthologs. This is a data source limitation, not a design bug.

**Recommendation:** Accept this asymmetry. The explorer UI should note that category granularity varies by organism group. The raw annotation fields (`cyanorak_Role`, `tIGR_Role`, `cog_category`) remain available on gene nodes for users who need finer or more consistent classification.

---

## Mapping Tables

### COG letter → gene_category

NOTE: COG letter codes and Cyanorak letter codes use the SAME letters for DIFFERENT
functions. E.g., COG "E" = Amino acid metabolism, Cyanorak "E" = Central intermediary
metabolism (N/P/S). This is intentional — the two classification systems are independent.

```python
# COG single-letter functional categories (NCBI/eggNOG system).
# WARNING: COG letters ≠ Cyanorak letters! Same letter, different meaning.
# E.g., COG E = "Amino acid metabolism", Cyanorak E = "Central intermediary metabolism"
COG_TO_CATEGORY = {
    "A": "Transcription",           # RNA processing (rare, ~4 genes)
    "B": "Replication and repair",  # Chromatin (rare, ~8 genes)
    "C": "Energy production",
    "D": "Cell cycle and division",
    "E": "Amino acid metabolism",   # ≠ Cyanorak E (central intermediary metabolism)
    "F": "Nucleotide metabolism",
    "G": "Carbohydrate metabolism",
    "H": "Coenzyme metabolism",
    "I": "Lipid metabolism",
    "J": "Translation",
    "K": "Transcription",
    "L": "Replication and repair",
    "M": "Cell wall and membrane",
    "N": "Cell motility",
    "O": "Post-translational modification",
    "P": "Inorganic ion transport",
    "Q": "Secondary metabolites",
    "R": "Unknown",                 # General function prediction only
    "S": "Unknown",                 # Function unknown
    "T": "Signal transduction",
    "U": "Intracellular trafficking",
    "V": "Defense mechanisms",
    "W": "Cell wall and membrane",  # Extracellular structures (rare)
    "X": "Mobile elements",         # Mobilome
    "Y": "Unknown",                 # Nuclear structure (rare)
    "Z": "Cell cycle and division", # Cytoskeleton (rare)
}
```

### Cyanorak role → gene_category

Top-level letter mapping, except for "D" which uses second-level sub-codes:

```python
# Cyanorak top-level letter → gene_category.
# WARNING: Cyanorak letters ≠ COG letters! Same letter, different meaning.
# E.g., Cyanorak E = "Central intermediary metabolism", COG E = "Amino acid metabolism"
CYANORAK_TO_CATEGORY = {
    "A": "Amino acid metabolism",
    "B": "Coenzyme metabolism",
    "C": "Cell wall and membrane",
    # "D" handled by CYANORAK_D_SUBCODES below
    "E": "Central intermediary metabolism",  # ≠ COG E! N, P, S metabolism, amino sugars, polyamines
    "F": "Replication and repair",
    "G": "Carbohydrate metabolism", # Carbon metabolism incl. glycolysis, TCA, CO2 fixation.
                                    # TCA/glycolysis are arguably "Energy production" but we keep
                                    # all of Cyanorak G together since COG G maps here too.
    "H": "Lipid metabolism",
    "I": "Mobile elements",
    "J": "Photosynthesis",          # cyanobacteria-specific!
    "K": "Translation",
    "L": "Post-translational modification",
    "M": "Nucleotide metabolism",
    "N": "Regulatory functions",
    "O": "Signal transduction",
    "P": "Transcription",
    "Q": "Transport",
    "R": "Unknown",
}

# Cyanorak "D" sub-code mapping (second level)
CYANORAK_D_SUBCODES = {
    "D.1": "Stress response and adaptation",   # 2,637 genes — environmental acclimation
    "D.2": "Cell cycle and division",           # 147 genes
    "D.3": "Cellular processes",                # 109 genes — detoxification/toxin resistance
    "D.4": "Post-translational modification",   # 157 genes — chaperones/heat shock
    "D.5": "Cell motility",                     # 19 genes — chemotaxis
    "D.6": "Cellular processes",                # 23 genes — programmed cell death
    "D.7": "Cellular processes",                # 36 genes — circadian clock
}
# Fallback for bare "D" or unknown D.x: "Cellular processes"
```

### TIGR MainRole → gene_category

```python
TIGR_TO_CATEGORY = {
    "Amino acid biosynthesis": "Amino acid metabolism",
    "Biosynthesis of cofactors, prosthetic groups, and carriers": "Coenzyme metabolism",
    "Cell envelope": "Cell wall and membrane",
    "Cellular processes": "Cellular processes",
    "Central intermediary metabolism": "Central intermediary metabolism",
    "DNA metabolism": "Replication and repair",
    "Disrupted reading frame": "Unknown",
    "Energy metabolism": "Energy production",
    "Fatty acid and phospholipid metabolism": "Lipid metabolism",
    "Hypothetical proteins": "Unknown",
    "Mobile and extrachromosomal element functions": "Mobile elements",
    "Not Found": "Unknown",
    "Protein fate": "Post-translational modification",
    "Protein synthesis": "Translation",
    "Purines, pyrimidines, nucleosides, and nucleotides": "Nucleotide metabolism",
    "Regulatory functions": "Regulatory functions",
    "Signal transduction": "Signal transduction",
    "Transcription": "Transcription",
    "Transport and binding proteins": "Transport",
    "Unclassified": "Unknown",
    "Unknown function": "Unknown",
}
```

---

## Implementation

### Where: `build_gene_annotations.py` → `AnnotationBuilder.build_merged()`

Add `gene_category` computation after the existing `annotation_quality` block (~line 437), using the same pattern.

### Steps

1. **Define mapping dicts** as module-level constants in `build_gene_annotations.py` (or a new helper in `utils/annotation_helpers.py` if they're too large)

2. **Compute gene_category** in `build_merged()`:
```python
# gene_category — high-level functional classification
gene_category = None

# Priority 1: Cyanorak Role (cyanobacteria only)
cyanorak_roles = result.get("cyanorak_Role", [])
if cyanorak_roles:
    code = cyanorak_roles[0]
    top_letter = code.split(".")[0]
    if top_letter == "D":
        # Parse D sub-codes for finer classification
        parts = code.split(".")
        sub_key = ".".join(parts[:2]) if len(parts) >= 2 else "D"
        gene_category = CYANORAK_D_SUBCODES.get(sub_key, "Cellular processes")
    else:
        gene_category = CYANORAK_TO_CATEGORY.get(top_letter)

# Priority 2: TIGR Role
if not gene_category or gene_category == "Unknown":
    tigr_descs = result.get("tIGR_Role_description", [])
    if tigr_descs:
        main_role = tigr_descs[0].split(" / ")[0].strip()
        cat = TIGR_TO_CATEGORY.get(main_role)
        if cat and cat != "Unknown":
            gene_category = cat

# Priority 3: COG category
if not gene_category or gene_category == "Unknown":
    cog_cats = result.get("cog_category", [])
    if cog_cats:
        cat = COG_TO_CATEGORY.get(cog_cats[0])
        if cat:
            gene_category = cat

result["gene_category"] = gene_category or "Unknown"
```

3. **Add build-time assertion** — catch mapping bugs immediately:
```python
VALID_CATEGORIES = frozenset(COG_TO_CATEGORY.values()) | frozenset(CYANORAK_TO_CATEGORY.values()) \
    | frozenset(CYANORAK_D_SUBCODES.values()) | frozenset(TIGR_TO_CATEGORY.values()) | {"Unknown"}

# In build_merged(), after computing gene_category:
assert result["gene_category"] in VALID_CATEGORIES, \
    f"Invalid gene_category {result['gene_category']!r} for {result.get('locus_tag')}"
```

5. **Add to schema**: In `config/schema_config.yaml`, add `gene_category` to the Gene node properties (type: `str`).

6. **Add to stats**: Track category distribution in `process_strain()` stats for the coverage report.

7. **Rebuild**: `bash scripts/prepare_data.sh --steps 2 --force` to regenerate all `gene_annotations_merged.json`.

### What NOT to change
- No changes to the YAML config (`gene_annotations_config.yaml`) — the source fields are already defined
- No changes to adapters — `cyanorak_ncbi_adapter.py` already reads from `gene_annotations_merged.json`
- No index changes needed — `gene_category` is for filtering/grouping, not full-text search

---

## Edge Cases

| Case | Handling |
|------|----------|
| Gene has no role/COG data at all | `"Unknown"` |
| Cyanorak Role is "R" (hypothetical) but TIGR/COG has real category | TIGR/COG wins (priority fallthrough skips "Unknown"). ~20% of Cyanorak-R genes get rescued by COG; 80% stay "Unknown" (both sources agree). |
| Multi-valued COG (e.g., `["L", "U"]`) | Use first value only. `cog_category` is stored as a list of single letters (the `split_cog_category` transform in `gene_annotations_config.yaml` splits "LU" → `["L", "U"]`). 5.7% of genes (1,676/29,236) have multi-valued COG; taking the first is standard practice (COG lists primary category first). |
| Multi-valued Cyanorak/TIGR roles | Use first value only (same rationale). |
| Cyanorak "0" prefix (rare, 3 genes) | Not in mapping → falls through to TIGR/COG |
| COG "R" (General function prediction) | Maps to "Unknown" — same as "S" |

---

## Validation

### Build-time assertion

The `VALID_CATEGORIES` assertion (step 3) catches mapping bugs during build. No gene can get a category not in the controlled set.

### Unit tests (`tests/test_build_gene_annotations.py`)

Add/update tests in the existing `TestBuildMerged` class. The test fixtures `GM`, `EG`, `UP` need `cyanorak_Role`, `cyanorak_Role_description`, `tIGR_Role_description` fields added to `GM` (gene_mapping source) and `COG_category` already exists in `EG`.

**Test cases for gene_category:**

```python
# ── gene_category ────────────────────────────────────────────────────────

def test_category_from_cyanorak_role(self):
    # Cyanorak Role "F.1" → top letter "F" → "Replication and repair"
    gm = dict(GM, cyanorak_Role="F.1", cyanorak_Role_description="DNA replication, recombination, and repair")
    merged = self.builder.build_merged(gm, EG, UP)
    assert merged["gene_category"] == "Replication and repair"

def test_category_cyanorak_d1_stress(self):
    # Cyanorak D.1.3 → sub-code "D.1" → "Stress response and adaptation"
    gm = dict(GM, cyanorak_Role="D.1.3", cyanorak_Role_description="Nitrogen")
    merged = self.builder.build_merged(gm, EG, UP)
    assert merged["gene_category"] == "Stress response and adaptation"

def test_category_cyanorak_d2_cell_division(self):
    gm = dict(GM, cyanorak_Role="D.2", cyanorak_Role_description="Cell division")
    merged = self.builder.build_merged(gm, EG, UP)
    assert merged["gene_category"] == "Cell cycle and division"

def test_category_cyanorak_d4_chaperones(self):
    gm = dict(GM, cyanorak_Role="D.4", cyanorak_Role_description="Chaperones and heat shock")
    merged = self.builder.build_merged(gm, EG, UP)
    assert merged["gene_category"] == "Post-translational modification"

def test_category_cyanorak_d_bare_fallback(self):
    # Bare "D" with no sub-code → "Cellular processes"
    gm = dict(GM, cyanorak_Role="D", cyanorak_Role_description="Diverse")
    merged = self.builder.build_merged(gm, EG, UP)
    assert merged["gene_category"] == "Cellular processes"

def test_category_cyanorak_j_photosynthesis(self):
    gm = dict(GM, cyanorak_Role="J.6", cyanorak_Role_description="NADH dehydrogenase")
    merged = self.builder.build_merged(gm, EG, UP)
    assert merged["gene_category"] == "Photosynthesis"

def test_category_cyanorak_e_central_intermediary(self):
    # Cyanorak E → "Central intermediary metabolism" (NOT COG E which is amino acid)
    gm = dict(GM, cyanorak_Role="E.4", cyanorak_Role_description="Nitrogen metabolism")
    merged = self.builder.build_merged(gm, EG, UP)
    assert merged["gene_category"] == "Central intermediary metabolism"

def test_category_cyanorak_r_falls_through_to_cog(self):
    # Cyanorak R = unknown → should fall through to COG (EG has COG_category="L")
    gm = dict(GM, cyanorak_Role="R.2", cyanorak_Role_description="Conserved hypothetical")
    merged = self.builder.build_merged(gm, EG, UP)
    assert merged["gene_category"] == "Replication and repair"  # from COG L

def test_category_cog_only_no_cyanorak(self):
    # No cyanorak_Role → falls through to COG
    merged = self.builder.build_merged(GM, EG, UP)
    assert merged["gene_category"] == "Replication and repair"  # EG has COG_category="L"

def test_category_cog_s_unknown(self):
    # COG S → "Unknown"
    eg = dict(EG, COG_category="S")
    merged = self.builder.build_merged(GM, eg, UP)
    assert merged["gene_category"] == "Unknown"

def test_category_no_sources_unknown(self):
    # No Cyanorak, no TIGR, no COG → "Unknown"
    merged = self.builder.build_merged(GM, {}, {})
    assert merged["gene_category"] == "Unknown"

def test_category_tigr_fallback(self):
    # Cyanorak R (unknown) + TIGR "Energy metabolism / Photosynthesis" → "Energy production"
    gm = dict(GM, cyanorak_Role="R.2", cyanorak_Role_description="Conserved hypothetical",
              tIGR_Role_description="Energy metabolism / Photosynthesis")
    eg = dict(EG, COG_category="S")  # COG also unknown
    merged = self.builder.build_merged(gm, eg, UP)
    assert merged["gene_category"] == "Energy production"

def test_category_tigr_unknown_falls_through_to_cog(self):
    # Cyanorak R + TIGR "Hypothetical proteins / Conserved" → both unknown, falls to COG
    gm = dict(GM, cyanorak_Role="R.2", cyanorak_Role_description="Conserved hypothetical",
              tIGR_Role_description="Hypothetical proteins / Conserved")
    merged = self.builder.build_merged(gm, EG, UP)
    assert merged["gene_category"] == "Replication and repair"  # from COG L

def test_category_value_in_valid_set(self):
    # Ensure all possible outputs are in VALID_CATEGORIES
    merged = self.builder.build_merged(GM, EG, UP)
    assert merged["gene_category"] in VALID_CATEGORIES

def test_category_multi_valued_cog_uses_first(self):
    # COG "LU" split to ["L", "U"] → uses first = "L" → "Replication and repair"
    eg = dict(EG, COG_category="LU")
    merged = self.builder.build_merged(GM, eg, UP)
    assert merged["gene_category"] == "Replication and repair"
```

**Note on fixtures:** `cyanorak_Role` and `tIGR_Role_description` are gene_mapping fields (stored as comma-delimited strings in CSV, parsed to lists by `passthrough_list` type in config). The test fixtures pass them as strings; the `passthrough_list` resolver in `build_merged` splits them. If the config for these fields is not in the test's `MINIMAL_CONFIG`, these tests need to add the field definitions or use the full config.

### Post-rebuild checks

1. Coverage report in build output shows distribution per strain
2. Quick Cypher check: `MATCH (g:Gene) RETURN g.gene_category, count(*) ORDER BY count(*) DESC`

---

## Decisions Log

1. **Cyanorak "D" handling** → RESOLVED: parse sub-codes. D.1→"Stress response and adaptation", D.2→"Cell cycle and division", D.4→"Post-translational modification" (chaperones), D.5→"Cell motility", D.3/D.6/D.7→"Cellular processes".
2. **Cyanorak "E" mapping** → RESOLVED: Cyanorak E is NOT energy production. It covers nitrogen (E.4), phosphorus (E.3), sulfur (E.7), amino sugars, polyamines. Maps to **"Central intermediary metabolism"**. TIGR "Central intermediary metabolism" maps here too.
3. **Category coverage asymmetry** → ACCEPTED: 6 categories are cyanobacteria-only (no COG mapping). Only "Photosynthesis" is genuinely organism-specific. The other 5 represent universal functions that Alteromonas has but will get coarser COG-derived labels. This is a data source limitation. Document in explorer UI.
4. **"Unknown" handling** → use `"Unknown"` string (not null) for consistency and filterability.
