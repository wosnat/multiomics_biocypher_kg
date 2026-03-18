# Gene Node ID and Name Fields — Specification & Overlap Analysis

## Overview

Gene nodes carry **4 identifier/name properties** on the KG node: `locus_tag`, `gene_name`, `gene_name_synonyms`, and `all_identifiers`. Three formerly-redundant fields (`gene_synonyms`, `alternative_locus_tags`, `old_locus_tags`) were removed because their content was fully contained in `all_identifiers` and/or `gene_name_synonyms`. They remain in `gene_annotations_merged.json` for the gene ID mapping pipeline.

## Fields on Gene Nodes (in Neo4j)

### 1. `locus_tag` (scalar, str) — **Node ID basis**

| Attribute | Value |
|---|---|
| Coverage | 35,226/35,226 (100%) |
| Source | `gene_mapping.csv` → `locus_tag` column (merged from NCBI GFF + Cyanorak GFF) |
| Node ID | `ncbigene:<locus_tag>` |
| Format | Varies by strain (see table below) |

**Locus tag formats per strain:**

| Strain | `locus_tag` format | Example | Origin |
|---|---|---|---|
| MED4 | `PMM####` (4-digit) | `PMM0001` | Cyanorak original |
| MIT9313 | `PMT####` or `AKG35_RS#####` | `PMT0001`, `AKG35_RS12315` | Mixed: Cyanorak genes → PMT, NCBI-only → RS |
| MIT9312 | `PMT9312_####` or `BFV95_RS#####` | `PMT9312_0001` | Mixed |
| AS9601 | `A9601_####` or `GBY92_RS#####` | `A9601_00010` | Mixed |
| MIT9301 | `P9301_####` or `RPA94_RS#####` | `P9301_00010` | Mixed |
| NATL1A | `NATL1_####` or `AMD04_RS#####` | `NATL1_00010` | Mixed |
| NATL2A | `PMN2A_####` or `HWB67_RS#####` | `PMN2A_00010` | Mixed |
| RSP50 | `BSR22_RS#####` | `BSR22_RS00005` | NCBI only (no Cyanorak) |
| CC9311 | `sync_####` or `SYNC_RS#####` | `sync_0001` | Mixed |
| WH8102 | `SYNW####` or `ATH98_RS#####` | `SYNW0001` | Mixed |
| EZ55 | `ALTBGP6_RS#####` | `ALTBGP6_RS00350` | NCBI only |
| HOT1A3 | `ACZ81_#####` | `ACZ81_00005` | NCBI only |
| MIT1002 | `MIT1002_#####` or `ALT831_RS#####` | `MIT1002_00001` | Mixed |

**Rule**: When a gene exists in both Cyanorak and NCBI, `locus_tag` is the Cyanorak-style ID (short, original naming). When a gene is NCBI-only, `locus_tag` is the NCBI `_RS` format. The merge logic in `gene_mapping.csv` picks Cyanorak's locus_tag as primary when available.

---

### 2. `gene_name` (scalar, str)

| Attribute | Value |
|---|---|
| Coverage | 18,868/35,226 (53.6%) |
| Source priority | Cyanorak `gene_names_cyanorak` → UniProt `gene_symbol` → NCBI `gene_names` → EggNOG `Preferred_name` |
| Content | Biological gene symbol (e.g., `dnaN`, `rpoA`, `crtP`, `hli`) |
| Filter | Identifier-style strings (matching `^[A-Z][A-Z0-9]*\d{4,}` or `.*_.*`) are **rejected** and set to null |

Provenance is tracked in `gene_name_source` (in merged JSON, not on KG node).

---

### 3. `gene_synonyms` (array, str[])

| Attribute | Value |
|---|---|
| Coverage | 33,904/35,226 (96.2%) |
| Source | UNION of: Cyanorak `gene_names_cyanorak`, NCBI `gene_names`, NCBI `gene_synonym`, NCBI `old_locus_tags`, UniProt `gene_names` |
| Content | **Mixed bag** — contains both gene-name tokens (`dnaJ`, `pds`) AND locus-tag-style tokens (`PMM0017`, `TX50_RS00370`) |
| Post-processing | Canonical `gene_name` is **excluded** (removed in build_gene_annotations.py) |

**Key relationship**: `gene_synonyms ⊇ alternative_locus_tags ∪ gene_name_synonyms` (exact equality for 99.3% of MED4 genes; 13 edge cases from UniProt paralog merging, e.g. hli18/hli8 sharing a protein record).

---

### 4. `gene_name_synonyms` (array, str[])

| Attribute | Value |
|---|---|
| Coverage | 1,511/35,226 (4.3%) — very sparse |
| Source | Same sources as `gene_synonyms` but with **filter** `^([A-Z][A-Z0-9]*\d{4,}|.+_.+)$` applied as exclusion (rejects locus-tag patterns) |
| Content | Gene-name-like tokens only (e.g., `pds`, `dnaJ`, `hli18`) |
| Post-processing | Canonical `gene_name` excluded |

**Relationship to `gene_synonyms`**: strict subset. `gene_name_synonyms ⊂ gene_synonyms` always.

**Per-strain coverage** ranges from 0 (RSP50, no Cyanorak) to 196 (WH8102). This field is only populated when a gene has alternative biological names (not just locus tag aliases).

---

### 5. `alternative_locus_tags` (array, str[])

| Attribute | Value |
|---|---|
| Coverage | 33,904/35,226 (96.2%) — identical to `gene_synonyms` |
| Source | Same sources as `gene_synonyms` but with **filter** `^([A-Z][A-Z0-9]*\d{4,}|.+_.+)$` applied as inclusion (keeps only locus-tag patterns). Also includes `old_locus_tags`. UniProt `gene_names` deliberately excluded to avoid false Tier 1 conflicts from paralog merging. |
| Content | Locus-tag-like tokens only (e.g., `PMM0001`, `TX50_RS00020`, `RG24_RS00005`) |

**Key relationships**:
- `alternative_locus_tags ⊆ gene_synonyms` (always true)
- `old_locus_tags ⊆ alternative_locus_tags` (99.97%: 32,637/32,644 in KG)
- `locus_tag ∈ alternative_locus_tags` for 33,892/33,904 genes (96.2% of total)

**Redundancy**: For 27,675/33,904 genes (81.6%), `gene_synonyms == alternative_locus_tags` exactly (no gene-name tokens present).

---

### 6. `old_locus_tags` (array, str[])

| Attribute | Value |
|---|---|
| Coverage | 32,644/35,226 (92.7%) |
| Source | `gene_mapping.csv` → `old_locus_tags` column (from NCBI GFF `old_locus_tag` qualifier) |
| Content | Previous-generation locus tags from NCBI annotation updates |

**Common pattern**: For MED4, `old_locus_tags == [locus_tag]` for 94% of genes (1857/1976). This happens because MED4's primary locus_tag (`PMM####`) is Cyanorak-derived, and NCBI's `old_locus_tag` GFF field points back to the same `PMM####` ID as a historical tag. The current NCBI tag is `TX50_RS#####`.

**Key relationship**: `old_locus_tags ⊆ alternative_locus_tags` (by construction — `old_locus_tags` is one of the sources fed into `alternative_locus_tags`).

---

### 7. `all_identifiers` (array, str[])

| Attribute | Value |
|---|---|
| Coverage | 35,226/35,226 (100%) |
| Content | Union of identifiers for MCP `get_gene` lookup |
| Excludes | `locus_tag` and `gene_name` (they have their own scalar indexes) |

**Composition** (from `build_gene_annotations.py:740-753`):
```
all_identifiers = {locus_tag_ncbi, locus_tag_cyanorak, protein_id}
                  ∪ old_locus_tags
                  ∪ alternative_locus_tags
                  ∪ gene_name_synonyms
                  - {locus_tag, gene_name}
```

**Notable exclusions from `all_identifiers`**:
- `uniprot_accession` — **never included** (not in the union logic)
- `gene_name` — excluded by design (has scalar index)
- `locus_tag` — excluded by design (has scalar index)

**Notable inclusions**:
- `locus_tag_ncbi` (e.g., `TX50_RS00020`) — always present when it differs from `locus_tag`
- `locus_tag_cyanorak` (e.g., `CK_Pro_MED4_00001`) — always present for Cyanorak strains
- `protein_id` (e.g., `WP_011131639.1`) — RefSeq protein accession

---

## Fields in `gene_annotations_merged.json` Only (NOT on KG node)

These are used during build/mapping but not emitted by the adapter:

| Field | Coverage (MED4) | Content |
|---|---|---|
| `locus_tag_ncbi` | 95% | Current NCBI `_RS`-format locus tag (e.g., `TX50_RS00020`) |
| `locus_tag_cyanorak` | 99.3% (Pro/Syn only) | Cyanorak internal ID (e.g., `CK_Pro_MED4_00001`) |
| `protein_id` | 95% | RefSeq WP_ protein accession |
| `uniprot_accession` | 66.4% | UniProt accession (e.g., `Q7V3R7`) |

These 4 fields **do** flow into `all_identifiers` (except `uniprot_accession`), so they're searchable on the KG node through that array.

---

## Where Does Each ID Land? — Concrete Examples

### Example 1: PMM0004 / purF (MED4, typical gene with name, no name synonyms)

| ID / token | `locus_tag` | `gene_name` | `gene_synonyms` | `gene_name_synonyms` | `alternative_locus_tags` | `old_locus_tags` | `all_identifiers` |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| `PMM0004` (primary locus) | **YES** | | YES | | YES | YES | |
| `purF` (bio name) | | **YES** | | | | | |
| `TX50_RS00035` (NCBI RS) | | | | | | | YES |
| `CK_Pro_MED4_00004` (Cyanorak) | | | | | | | YES |
| `WP_011131642.1` (protein) | | | | | | | YES |
| `Q7V3R6` (UniProt) | | | | | | | **missing** (bug) |

**Note**: `PMM0004` appears in 4 places: `locus_tag`, `gene_synonyms`, `alternative_locus_tags`, `old_locus_tags`. It's excluded from `all_identifiers` (by design — `locus_tag` has its own index). `purF` appears only in `gene_name` (excluded from `all_identifiers` by design).

### Example 2: PMM0565 / dnaA (MED4, well-known gene, no name synonyms)

| ID / token | `locus_tag` | `gene_name` | `gene_synonyms` | `gene_name_synonyms` | `alternative_locus_tags` | `old_locus_tags` | `all_identifiers` |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| `PMM0565` | **YES** | | YES | | YES | YES | |
| `dnaA` | | **YES** | | | | | |
| `TX50_RS03015` | | | | | | | YES |
| `CK_Pro_MED4_00565` | | | | | | | YES |
| `WP_011132199.1` | | | | | | | YES |

Same pattern as PMM0004 — `dnaA` is only in `gene_name`. No synonyms because this gene has a single unambiguous name.

### Example 3: PMM0198 / hisC/cobC (MED4, gene with name synonyms)

| ID / token | `locus_tag` | `gene_name` | `gene_synonyms` | `gene_name_synonyms` | `alternative_locus_tags` | `old_locus_tags` | `all_identifiers` |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| `PMM0198` | **YES** | | YES | | YES | YES | |
| `hisC/cobC` (canonical) | | **YES** | | | | | |
| `cobC` (synonym) | | | YES | YES | | | YES |
| `hisC` (synonym) | | | YES | YES | | | YES |
| `TX50_RS01025` (NCBI RS) | | | YES | | YES | | YES |
| `CK_Pro_MED4_00198` | | | | | | | YES |
| `WP_011131837.1` | | | | | | | YES |

Here `gene_name_synonyms` adds value — `cobC` and `hisC` are real biological names discoverable through it. `TX50_RS01025` appears in `gene_synonyms` AND `alternative_locus_tags` AND `all_identifiers` (triple).

### Example 4: PMT0001 / dnaN (MIT9313, gene with multiple old locus tags)

| ID / token | `locus_tag` | `gene_name` | `gene_synonyms` | `gene_name_synonyms` | `alternative_locus_tags` | `old_locus_tags` | `all_identifiers` |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| `PMT0001` | **YES** | | YES | | YES | YES | |
| `dnaN` | | **YES** | | | | | |
| `PMT_0001` (underscore variant) | | | YES | | YES | YES | YES |
| `RG24_RS00005` (old NCBI RS) | | | YES | | YES | YES | YES |
| `AKG35_RS00005` (current NCBI RS) | | | | | | | YES |
| `CK_Pro_MIT9313_00001` | | | | | | | YES |
| `WP_011129380.1` | | | | | | | YES |

MIT9313 shows how `old_locus_tags` can have multiple entries (`PMT0001`, `PMT_0001`, `RG24_RS00005`). All appear in `gene_synonyms` and `alternative_locus_tags` too.

### Example 5: EZ55_01745 / rne (Alteromonas EZ55, no Cyanorak)

| ID / token | `locus_tag` | `gene_name` | `gene_synonyms` | `gene_name_synonyms` | `alternative_locus_tags` | `old_locus_tags` | `all_identifiers` |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| `EZ55_01745` | **YES** | | YES | | YES | YES | |
| `rne` | | **YES** | | | | | |
| `ALTBGP6_RS08545` (NCBI RS) | | | | | | | YES |
| `WP_156086307.1` | | | | | | | YES |

Alteromonas has no Cyanorak data. `locus_tag` = old-style from original annotation. `locus_tag_ncbi` (RS format) only surfaces in `all_identifiers`.

### Example 6: PMM0217 / mutS2 (MED4, gene with name synonyms from different sources)

| ID / token | `locus_tag` | `gene_name` | `gene_synonyms` | `gene_name_synonyms` | `alternative_locus_tags` | `old_locus_tags` | `all_identifiers` |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| `PMM0217` | **YES** | | YES | | YES | YES | |
| `mutS2` (canonical) | | **YES** | | | | | |
| `mutS` (synonym) | | | YES | YES | | | YES |
| `rqcU` (synonym) | | | YES | YES | | | YES |
| `TX50_RS01120` (NCBI RS) | | | YES | | YES | | YES |
| `CK_Pro_MED4_00217` | | | | | | | YES |
| `WP_011131856.1` | | | | | | | YES |

### Pattern Summary (after cleanup)

| ID type | Where it lands on KG node | Notes |
|---|---|---|
| **Primary locus tag** (e.g., `PMM0004`) | `locus_tag` (scalar) | Excluded from `all_identifiers` by design (has own index) |
| **Canonical gene name** (e.g., `dnaA`) | `gene_name` (scalar) | Excluded from `all_identifiers` by design (has own index) |
| **Gene name synonyms** (e.g., `cobC`, `mutS`) | `gene_name_synonyms` + `all_identifiers` | 2x (acceptable — enables "bio names only" queries) |
| **NCBI RS locus tag** (e.g., `TX50_RS01025`) | `all_identifiers` | 1x |
| **Old locus tags** (e.g., `PMT_0001`) | `all_identifiers` | 1x |
| **Cyanorak ID** (e.g., `CK_Pro_MED4_00004`) | `all_identifiers` | 1x |
| **NCBI RS (current)** via `locus_tag_ncbi` | `all_identifiers` | 1x |
| **Protein ID** (e.g., `WP_011131642.1`) | `all_identifiers` | 1x |
| **UniProt accession** (e.g., `Q7V3R6`) | `all_identifiers` | 1x (fixed — was missing) |

---

## Architecture: What's on the KG Node vs JSON Only

```
KG Gene node (4 ID/name fields):
  locus_tag           ← scalar, indexed (gene_locus_tag_idx)
  gene_name           ← scalar, indexed (gene_name_idx)
  gene_name_synonyms  ← array, bio names only (e.g. ["cobC", "hisC"])
  all_identifiers     ← array, all alt IDs (locus tags, protein IDs, Cyanorak IDs, name synonyms)

Full-text index (geneFullText):
  gene_summary + all_identifiers + gene_name_synonyms + alternate_functional_descriptions

gene_annotations_merged.json (build pipeline only, NOT on KG node):
  gene_synonyms           ← superset of gene_name_synonyms ∪ alternative_locus_tags
  alternative_locus_tags  ← locus-tag-pattern tokens (subset of all_identifiers)
  old_locus_tags          ← NCBI GFF old_locus_tag values (subset of alternative_locus_tags)
  locus_tag_ncbi          ← current NCBI RS-format tag
  locus_tag_cyanorak      ← Cyanorak internal ID
  protein_id              ← RefSeq WP_ accession
  uniprot_accession       ← UniProt accession
```

---

## Redundancy Analysis (historical — cleanup DONE)

The following fields were removed from the KG node because they were fully redundant:

### 1. `gene_synonyms` — REMOVED from KG

Was: union of all alternative names (mixed gene names + locus tags). For 81.6% of genes, `gene_synonyms == alternative_locus_tags` (no gene-name tokens). For the rest, `gene_synonyms = alternative_locus_tags ∪ gene_name_synonyms`. All tokens are now in `all_identifiers` and/or `gene_name_synonyms`.

### 2. `alternative_locus_tags` — REMOVED from KG

Was: locus-tag-pattern tokens only. Strict subset of `all_identifiers`.

### 3. `old_locus_tags` — REMOVED from KG

Was: historical NCBI locus tags. Strict subset of `alternative_locus_tags` ⊂ `all_identifiers`.

All three remain in `gene_annotations_merged.json` for the gene ID mapping pipeline (`gene_id_utils.py`, `build_gene_id_mapping.py`).

---

---

## Summary Table

| Field | On KG? | Coverage | Content | Indexed? |
|---|---|---|---|---|
| `locus_tag` | Yes (scalar) | 100% | Primary ID, used as node ID basis (`ncbigene:<locus_tag>`) | Scalar + full-text (via gene_summary) |
| `gene_name` | Yes (scalar) | 53.6% | Canonical biological name (e.g., `dnaA`, `crtP`). Priority: Cyanorak > UniProt > NCBI > EggNOG. Identifier-like values filtered out. | Scalar (gene_name_idx) |
| `gene_name_synonyms` | Yes (array) | 4.3% | Alternative biological names only — no locus tags. E.g., `["cobC", "hisC"]` for a gene whose canonical name is `hisC/cobC`. Canonical `gene_name` excluded. | Full-text (geneFullText) |
| `all_identifiers` | Yes (array) | 100% | Union of all alternative IDs for search: `locus_tag_ncbi`, `locus_tag_cyanorak`, `protein_id`, old locus tags, alternative locus tags, gene name synonyms. Excludes `locus_tag` and `gene_name` (they have scalar indexes). | Full-text (geneFullText) |
| `locus_tag_ncbi` | JSON only | ~95% | Current NCBI `_RS`-format locus tag | Via `all_identifiers` |
| `locus_tag_cyanorak` | JSON only | ~99% Pro/Syn | Cyanorak internal ID (e.g., `CK_Pro_MED4_00001`) | Via `all_identifiers` |
| `protein_id` | JSON only | ~95% | RefSeq WP_ protein accession | Via `all_identifiers` |
| `uniprot_accession` | JSON only | ~66% | UniProt accession | Via `all_identifiers` (fixed) |
| `gene_synonyms` | JSON only | ~96% | Mixed gene names + locus tags (= alt_locus_tags ∪ gene_name_synonyms) | Removed from KG |
| `alternative_locus_tags` | JSON only | ~96% | Locus-tag-pattern tokens | Removed from KG |
| `old_locus_tags` | JSON only | ~93% | Historical NCBI locus tags | Removed from KG |
