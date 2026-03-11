# Phase 1: Paperconfig Standardization

**Master plan:** `/home/osnat/github/multiomics_biocypher_kg/plans/schema_improvements_for_mcp.md`
**Status:** READY TO START
**Gate:** Agent G approval required before Phase 2 begins.

---

## Overview

Normalize all 24 active paperconfig files to use canonical vocabulary.
Update `validate_paperconfig.py` to enforce that vocabulary.
No adapter or schema code changes in this phase.

Active paperconfigs (23 papers + 1 strain resource file):
- `data/Prochlorococcus/papers_and_supp/MIT9313_resources/paperconfig.yaml` — no publication block; no changes needed
- `data/Prochlorococcus/papers_and_supp/Aharonovich 2016/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/Al-Hosani 2015/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/Anjur 2025/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/bagby 2015/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/barreto 2022/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/biller 2016/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/coe 2024/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/Fang 2019/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/he 2022/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/Hennon 2017/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/Lin 2015/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/lindell 2007/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/martiny 2006/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/steglich 2006/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/tetu 2019/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/thompson 2011/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/tolonen 2006/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/Read 2017/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/Thompson 2016/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/Weissberg_2025/paperconfig.yaml`
- `data/Prochlorococcus/papers_and_supp/Labban 2022/paperconfig.yaml` — commented out; skip

Note: `Labban 2022` is commented out in `paperconfig_files.txt` and is excluded from this phase.

---

## Canonical Vocabulary (authoritative reference)

### Canonical organism names (exhaustive allowed set)

For the `organism` field (the organism whose genes are measured):
```
Prochlorococcus MED4
Prochlorococcus AS9601
Prochlorococcus MIT9301
Prochlorococcus MIT9312
Prochlorococcus MIT9313
Prochlorococcus NATL1A
Prochlorococcus NATL2A
Prochlorococcus RSP50
Synechococcus WH8102
Synechococcus CC9311
Alteromonas macleodii HOT1A3
Alteromonas macleodii EZ55
Alteromonas macleodii MIT1002
```

Note: `Alteromonas macleodii MIT1002` (not `Alteromonas MIT1002`) is the canonical form. However, the existing paperconfigs already use `Alteromonas macleodii MIT1002` throughout biller 2016 and coe 2024 — verify during review.

For the `treatment_organism` field (the interacting organism), the same names apply plus treatment-only organisms from `treatment_organisms.csv`:
```
Phage
Marinobacter
Thalassospira
Pseudohoeflea
Alteromonas
```

### Canonical `condition_type` enum

```
growth_medium
nitrogen_stress
phosphorus_stress
iron_stress
salt_stress
carbon_stress
light_stress
darkness
plastic_stress
viral
coculture
growth_state
temperature_stress
```

### Canonical `test_type` enum

```
DESeq2
DESeq
edgeR
Rockhopper
microarray
microarray_Cyber-T
microarray_LPE
microarray_Goldenspike
```

---

## Per-File Change Inventory

This section lists every non-canonical value found by reading each file. Files with no changes needed are noted explicitly.

---

### MIT9313_resources/paperconfig.yaml
**No publication block. No `statistical_analyses`. No changes needed.**

---

### Aharonovich 2016/paperconfig.yaml

**Issues found:**

**1. organism field (line 57, 76):** `Alteromonas HOT1A3` → `Alteromonas macleodii HOT1A3`
The organism field in `supp_table_hot1a3_9313_vs_9313dil` and `supp_table_hot1a3_9313_vs_med4` reads:
```yaml
          organism: Alteromonas HOT1A3
```
Change to:
```yaml
          organism: Alteromonas macleodii HOT1A3
```

**2. treatment_organism field (lines 58, 78):** `Prochlorococcus marinus MIT9313` → `Prochlorococcus MIT9313`
In both `supp_table_hot1a3_9313_vs_9313dil` and `supp_table_hot1a3_9313_vs_med4`:
```yaml
          treatment_organism: Prochlorococcus marinus MIT9313
```
Change to:
```yaml
          treatment_organism: Prochlorococcus MIT9313
```

**3. treatment_condition in supp_table_2 (line 98) and supp_table_3 (line 123):** These read `Coculture with Alteromonas HOT1A3` — this is not an organism field, it is a free-text condition description string, not subject to the organism name enum. Leave as-is; it refers to the condition name, not the organism field. However, for consistency with canonical names in human-readable strings, consider updating to `Alteromonas macleodii HOT1A3`. Decision: **update treatment_condition free-text strings** in supp_table_2 and supp_table_3 to use the canonical organism name. These analyses belong to Prochlorococcus organisms so they use `treatment_organism: Alteromonas macleodii HOT1A3` (already correct on those rows). The treatment_condition text strings are free-form descriptions.

**Summary of changes for Aharonovich 2016:**

| Location | Field | Before | After |
|---|---|---|---|
| `supp_table_hot1a3_9313_vs_9313dil` > analysis | `organism` | `Alteromonas HOT1A3` | `Alteromonas macleodii HOT1A3` |
| `supp_table_hot1a3_9313_vs_med4` > analysis | `organism` | `Alteromonas HOT1A3` | `Alteromonas macleodii HOT1A3` |
| `supp_table_hot1a3_9313_vs_9313dil` > analysis | `treatment_organism` | `Prochlorococcus marinus MIT9313` | `Prochlorococcus MIT9313` |
| `supp_table_hot1a3_9313_vs_med4` > analysis | `treatment_organism` | `Prochlorococcus marinus MIT9313` | `Prochlorococcus MIT9313` |

**Exact search/replace strings for Agent A:**

Search 1 (appears twice — use replace_all):
```
          organism: Alteromonas HOT1A3
```
Replace with:
```
          organism: Alteromonas macleodii HOT1A3
```

Search 2 (appears twice — use replace_all):
```
          treatment_organism: Prochlorococcus marinus MIT9313
```
Replace with:
```
          treatment_organism: Prochlorococcus MIT9313
```

---

### Al-Hosani 2015/paperconfig.yaml

**Issues found:**

**1. condition_type `salt_stress` (line 16):** Already canonical. No change needed.

**No changes needed.**

Wait — let me verify: `alhosani_2015_seawater_control` has `condition_type: growth_medium` (canonical). `alhosani_2015_salt_acclimation` has `condition_type: salt_stress` (canonical). `test_type: DESeq` (canonical). `organism: Prochlorococcus AS9601` (canonical). All fields present.

**Result: NO CHANGES NEEDED.**

---

### Anjur 2025/paperconfig.yaml

**Issues found:**

All `condition_type` values are `growth_state` — canonical.
`test_type: edgeR` — canonical.
`organism: Prochlorococcus MIT9301` — canonical.
All required fields present (`id`, `type`, `treatment_condition`, `organism`, `name_col`, `logfc_col`, `test_type`, `control_condition`).

**Result: NO CHANGES NEEDED.**

---

### bagby 2015/paperconfig.yaml

**Issues found:**

**1. condition_type `gas_shock` (lines 7, 13, 20, 27):** Non-canonical. Must change to `carbon_stress`.

The file has four `environmental_conditions` entries all using `condition_type: gas_shock`:
- `bagby_2015_air`
- `bagby_2015_no_o2`
- `bagby_2015_low_co2`
- `bagby_2015_low_co2_no_o2`

**2. test_type `Affymetrix microarray` (lines 51, 77, 103, 128):** Non-canonical. Must change to `microarray`.

**Summary of changes for bagby 2015:**

| Location | Field | Before | After |
|---|---|---|---|
| All 4 env conditions | `condition_type` | `gas_shock` | `carbon_stress` |
| All 4 analyses | `test_type` | `Affymetrix microarray` | `microarray` |

**Exact search/replace strings for Agent A:**

Search 1 (appears 4 times — use replace_all):
```
      condition_type: gas_shock
```
Replace with:
```
      condition_type: carbon_stress
```

Search 2 (appears 4 times — use replace_all):
```
          test_type: Affymetrix microarray
```
Replace with:
```
          test_type: microarray
```

---

### barreto 2022/paperconfig.yaml

**Issues found:**

**1. condition_type `pco2` (lines 7, 13):** Non-canonical. Must change to `carbon_stress`.

Both `barreto_2022_ambient_pco2` and `barreto_2022_elevated_pco2` use `condition_type: pco2`.

**Summary of changes for barreto 2022:**

| Location | Field | Before | After |
|---|---|---|---|
| `barreto_2022_ambient_pco2` | `condition_type` | `pco2` | `carbon_stress` |
| `barreto_2022_elevated_pco2` | `condition_type` | `pco2` | `carbon_stress` |

All `organism` fields use `Alteromonas macleodii EZ55`, `Prochlorococcus MIT9312`, `Synechococcus CC9311`, `Synechococcus WH8102` — all canonical.
All `test_type: edgeR` — canonical.
All analyses have `treatment_condition` — present.

**Exact search/replace strings for Agent A:**

Search 1 (appears 2 times — use replace_all):
```
      condition_type: pco2
```
Replace with:
```
      condition_type: carbon_stress
```

---

### biller 2016/paperconfig.yaml

**Issues found:**

**1. condition_type `coculture` (lines 15, 22):** The `coculture` value IS in the canonical enum — canonical. No change needed for condition_type.

All `organism` fields: `Prochlorococcus NATL2A`, `Alteromonas macleodii MIT1002` — canonical.
All `test_type: DESeq2` — canonical.
All required fields present.

**Result: NO CHANGES NEEDED.**

---

### Biller 2018/paperconfig.yaml

**Issues found:**

**1. condition_type `light_stress` for extended darkness (line 15):** The `biller_2018_extended_darkness` condition uses `condition_type: light_stress`. Per the master plan, extended darkness maps to canonical value `darkness`. Must change.

The `biller_2018_diel_ld_control` condition uses `condition_type: growth_medium` — canonical.

All `organism` fields: `Prochlorococcus NATL2A`, `Alteromonas macleodii MIT1002` — canonical.
All `test_type: DESeq2` — canonical.
All required fields present.

**Summary of changes for Biller 2018:**

| Location | Field | Before | After |
|---|---|---|---|
| `biller_2018_extended_darkness` | `condition_type` | `light_stress` | `darkness` |

**Exact search/replace strings for Agent A:**

Search 1 (unique in file — verify before replace):
```
    biller_2018_extended_darkness:
      condition_type: light_stress
```
Replace with:
```
    biller_2018_extended_darkness:
      condition_type: darkness
```

---

### coe 2024/paperconfig.yaml

**Issues found:**

**1. condition_type `dark_tolerance` (line 15):** Non-canonical. Must change to `darkness`.

The `coe_2024_dark_tolerant_diel` condition uses `condition_type: dark_tolerance`.
The `coe_2024_parental_diel` condition uses `condition_type: growth_medium` — canonical.

All `organism` fields: `Prochlorococcus NATL2A`, `Alteromonas macleodii MIT1002` — canonical.
All `test_type: DESeq2` — canonical.
All required fields present.

**Summary of changes for coe 2024:**

| Location | Field | Before | After |
|---|---|---|---|
| `coe_2024_dark_tolerant_diel` | `condition_type` | `dark_tolerance` | `darkness` |

**Exact search/replace strings for Agent A:**

Search 1 (unique in file):
```
      condition_type: dark_tolerance
```
Replace with:
```
      condition_type: darkness
```

---

### Fang 2019/paperconfig.yaml

**Issues found:**

**1. condition_type `viral_lysis_products` (line 15):** Non-canonical. Must change to `viral`.

The `fang_2019_vdom_addition` condition uses `condition_type: viral_lysis_products`.
The `fang_2019_pro99_control` uses `condition_type: growth_medium` — canonical.

All `organism` fields: `Prochlorococcus MIT9313` — canonical.
All `test_type: DESeq2` — canonical.
All required fields present.

**Summary of changes for Fang 2019:**

| Location | Field | Before | After |
|---|---|---|---|
| `fang_2019_vdom_addition` | `condition_type` | `viral_lysis_products` | `viral` |

**Exact search/replace strings for Agent A:**

Search 1 (unique in file):
```
      condition_type: viral_lysis_products
```
Replace with:
```
      condition_type: viral
```

---

### he 2022/paperconfig.yaml

**Issues found:**

All `condition_type` values: `growth_medium`, `salt_stress` — canonical.
All `test_type: DESeq` — canonical.
All `organism` values: `Prochlorococcus NATL1A`, `Prochlorococcus MED4` — canonical.
All required fields present.

**Result: NO CHANGES NEEDED.**

---

### Hennon 2017/paperconfig.yaml

**Issues found:**

**1. condition_type `gas_shock` (lines 7, 15):** Non-canonical. Must change to `carbon_stress`.

Both `hennon_2017_ambient_co2` and `hennon_2017_elevated_co2` use `condition_type: gas_shock`.

**2. organism field `Alteromonas EZ55` (lines 95, 119):** Non-canonical. Must change to `Alteromonas macleodii EZ55`.

In `supp_table_3` and `supp_table_4` analyses:
```yaml
          organism: "Alteromonas EZ55"
```

**3. id_translation_ez55_author organism field (line 41):** `organism: "Alteromonas EZ55"` — this is in an `id_translation` entry, not a `statistical_analyses` entry. However, the organism field in `id_translation` blocks must also use canonical names (they determine which strain's gene ID mapping is used). Change to `Alteromonas macleodii EZ55`.

**Summary of changes for Hennon 2017:**

| Location | Field | Before | After |
|---|---|---|---|
| `hennon_2017_ambient_co2` | `condition_type` | `gas_shock` | `carbon_stress` |
| `hennon_2017_elevated_co2` | `condition_type` | `gas_shock` | `carbon_stress` |
| `id_translation_ez55_author` | `organism` | `"Alteromonas EZ55"` | `"Alteromonas macleodii EZ55"` |
| `supp_table_3` analysis | `organism` | `"Alteromonas EZ55"` | `"Alteromonas macleodii EZ55"` |
| `supp_table_4` analysis | `organism` | `"Alteromonas EZ55"` | `"Alteromonas macleodii EZ55"` |

All `test_type: edgeR` — canonical.
All required fields present (note: `supp_table_2` and `supp_table_3` have no `adjusted_p_value_col` but have `prefiltered: true` — acceptable).

**Exact search/replace strings for Agent A:**

Search 1 (appears 2 times — use replace_all):
```
      condition_type: gas_shock
```
Replace with:
```
      condition_type: carbon_stress
```

Search 2 (appears 3 times — use replace_all, organism field only):
```
      organism: "Alteromonas EZ55"
```
Replace with:
```
      organism: "Alteromonas macleodii EZ55"
```

Note: Search 2 with 6 leading spaces matches both the `id_translation` entry organism and the analysis organism fields. Verify all 3 occurrences after replace.

---

### Lin 2015/paperconfig.yaml

**Issues found:**

**1. condition_type `nutrient_stress` for phosphorus limitation (lines 11):** Non-canonical. In this paper the context is phosphorus limitation (`lin_2015_p_limited` description says "phosphorus limitation"). Must change to `phosphorus_stress`.

The `lin_2015_p_limited` condition uses `condition_type: nutrient_stress` with description "Prochlorococcus NATL2A under phosphorus limitation, P depleted from growth medium". Canonical form: `phosphorus_stress`.

The `lin_2015_p_replete_control` uses `condition_type: growth_medium` — canonical.

All `organism` fields: `Prochlorococcus NATL2A` — canonical.
All `test_type: DESeq2` — canonical.
All required fields present.

**Summary of changes for Lin 2015:**

| Location | Field | Before | After |
|---|---|---|---|
| `lin_2015_p_limited` | `condition_type` | `nutrient_stress` | `phosphorus_stress` |

**Exact search/replace strings for Agent A:**

Search 1 (unique in file — only one `nutrient_stress` occurrence):
```
      condition_type: nutrient_stress
```
Replace with:
```
      condition_type: phosphorus_stress
```

---

### lindell 2007/paperconfig.yaml

**Issues found:**

**1. test_type `Affymetrix microarray` (all 8 analyses):** Non-canonical. Must change to `microarray`.

The file has 8 analyses all with `test_type: "Affymetrix microarray"`.

All `organism` fields: `Prochlorococcus MED4` — canonical.
All required fields present.

**Summary of changes for lindell 2007:**

| Location | Field | Before | After |
|---|---|---|---|
| All 8 analyses | `test_type` | `"Affymetrix microarray"` | `microarray` |

**Exact search/replace strings for Agent A:**

Search 1 (appears 8 times — use replace_all):
```
          test_type: "Affymetrix microarray"
```
Replace with:
```
          test_type: microarray
```

---

### martiny 2006/paperconfig.yaml

**Issues found:**

**1. condition_type `nutrient_stress` for phosphorus starvation (line 15):** Non-canonical. The `martiny_2006_p_starvation` condition description says "P-depleted Pro99 medium". Canonical form: `phosphorus_stress`.

**2. test_type `Affymetrix microarray with Cyber-T Bayesian t-test` (all 9 analyses):** Non-canonical. Must change to `microarray_Cyber-T`.

All `organism` fields: `Prochlorococcus MED4`, `Prochlorococcus MIT9313` — canonical.
All required fields present.

**Summary of changes for martiny 2006:**

| Location | Field | Before | After |
|---|---|---|---|
| `martiny_2006_p_starvation` | `condition_type` | `nutrient_stress` | `phosphorus_stress` |
| All 9 analyses | `test_type` | `"Affymetrix microarray with Cyber-T Bayesian t-test"` | `microarray_Cyber-T` |

**Exact search/replace strings for Agent A:**

Search 1 (unique in file):
```
      condition_type: nutrient_stress
```
Replace with:
```
      condition_type: phosphorus_stress
```

Search 2 (appears 9 times — use replace_all):
```
          test_type: "Affymetrix microarray with Cyber-T Bayesian t-test"
```
Replace with:
```
          test_type: microarray_Cyber-T
```

---

### steglich 2006/paperconfig.yaml

**Issues found:**

All `condition_type` values: `light_stress` — canonical (this is actual light stress, not darkness).
All `test_type: "Affymetrix microarray with LPE test"` — non-canonical. Must change to `microarray_LPE`.

All `organism` fields: `Prochlorococcus MED4` — canonical.
All required fields present.

**Summary of changes for steglich 2006:**

| Location | Field | Before | After |
|---|---|---|---|
| All 6 analyses | `test_type` | `"Affymetrix microarray with LPE test"` | `microarray_LPE` |

**Exact search/replace strings for Agent A:**

Search 1 (appears 6 times — use replace_all):
```
          test_type: "Affymetrix microarray with LPE test"
```
Replace with:
```
          test_type: microarray_LPE
```

---

### tetu 2019/paperconfig.yaml

**Issues found:**

**1. condition_type `plastic_leachate_stress` (lines 15, 24):** Non-canonical. Must change to `plastic_stress`.

Both `tetu_2019_hdpe_leachate` and `tetu_2019_pvc_leachate` use `condition_type: plastic_leachate_stress`.

All `organism` fields: `Prochlorococcus MIT9312`, `Prochlorococcus NATL2A` — canonical.
All `test_type: DESeq2` — canonical.
All required fields present.

**Summary of changes for tetu 2019:**

| Location | Field | Before | After |
|---|---|---|---|
| `tetu_2019_hdpe_leachate` | `condition_type` | `plastic_leachate_stress` | `plastic_stress` |
| `tetu_2019_pvc_leachate` | `condition_type` | `plastic_leachate_stress` | `plastic_stress` |

**Exact search/replace strings for Agent A:**

Search 1 (appears 2 times — use replace_all):
```
      condition_type: plastic_leachate_stress
```
Replace with:
```
      condition_type: plastic_stress
```

---

### thompson 2011/paperconfig.yaml

**Issues found:**

**1. test_type `Affymetrix microarray with Bayesian t-test (Cyber-T)` (all 8 analyses):** Non-canonical. Must change to `microarray_Cyber-T`.

Note: This is a slightly different phrasing from martiny 2006 (`with Bayesian t-test (Cyber-T)` vs `with Cyber-T Bayesian t-test`). Both map to `microarray_Cyber-T`.

All `condition_type` values: `growth_medium`, `iron_stress` — canonical.
All `organism` fields: `Prochlorococcus MED4`, `Prochlorococcus MIT9313` — canonical.
All required fields present.

**Summary of changes for thompson 2011:**

| Location | Field | Before | After |
|---|---|---|---|
| All 8 analyses | `test_type` | `"Affymetrix microarray with Bayesian t-test (Cyber-T)"` | `microarray_Cyber-T` |

**Exact search/replace strings for Agent A:**

Search 1 (appears 8 times — use replace_all):
```
          test_type: "Affymetrix microarray with Bayesian t-test (Cyber-T)"
```
Replace with:
```
          test_type: microarray_Cyber-T
```

---

### tolonen 2006/paperconfig.yaml

**Issues found:**

**1. condition_type `nutrient_stress` for nitrogen deprivation (lines 16, 51):** Non-canonical. The nitrogen deprivation contexts (`tolonen_2006_n_deprivation_med4`, `tolonen_2006_n_deprivation_mit9313`) must change to `nitrogen_stress`.

**2. test_type `Affymetrix microarray with Goldenspike` (all 16 analyses):** Non-canonical. Must change to `microarray_Goldenspike`.

All `organism` fields: `Prochlorococcus MED4`, `Prochlorococcus MIT9313` — canonical.
All required fields present.

**Summary of changes for tolonen 2006:**

| Location | Field | Before | After |
|---|---|---|---|
| `tolonen_2006_n_deprivation_med4` | `condition_type` | `nutrient_stress` | `nitrogen_stress` |
| `tolonen_2006_n_deprivation_mit9313` | `condition_type` | `nutrient_stress` | `nitrogen_stress` |
| All 16 analyses | `test_type` | `"Affymetrix microarray with Goldenspike"` | `microarray_Goldenspike` |

**Exact search/replace strings for Agent A:**

Search 1 (appears 2 times — use replace_all):
```
      condition_type: nutrient_stress
```
Replace with:
```
      condition_type: nitrogen_stress
```

Search 2 (appears 16 times — use replace_all):
```
          test_type: "Affymetrix microarray with Goldenspike"
```
Replace with:
```
          test_type: microarray_Goldenspike
```

---

### Read 2017/paperconfig.yaml

**Issues found:**

**1. condition_type `nutrient_stress` for nitrogen depletion (line 16):** Non-canonical. The `read_2017_n_depleted` condition is nitrogen starvation — must change to `nitrogen_stress`.

**2. test_type `Rockhopper/DESeq` (all 3 analyses):** Non-canonical. Must change to `Rockhopper`.

All `organism` fields: `Prochlorococcus MED4` — canonical.
All required fields present.

**Summary of changes for Read 2017:**

| Location | Field | Before | After |
|---|---|---|---|
| `read_2017_n_depleted` | `condition_type` | `nutrient_stress` | `nitrogen_stress` |
| All 3 analyses | `test_type` | `"Rockhopper/DESeq"` | `Rockhopper` |

**Exact search/replace strings for Agent A:**

Search 1 (unique in file):
```
      condition_type: nutrient_stress
```
Replace with:
```
      condition_type: nitrogen_stress
```

Search 2 (appears 3 times — use replace_all):
```
          test_type: "Rockhopper/DESeq"
```
Replace with:
```
          test_type: Rockhopper
```

---

### Thompson 2016/paperconfig.yaml

**Issues found:**

**1. test_type `DESeq2/NOISeq` (all 13 analyses):** Non-canonical. Must change to `DESeq2`.

All `condition_type` values: `growth_medium`, `light_stress` — canonical.
All `organism` fields: `Prochlorococcus MED4` — canonical.
All required fields present (both `treatment_organism` analyses and `environmental_treatment_condition_id` analyses are valid).

**Summary of changes for Thompson 2016:**

| Location | Field | Before | After |
|---|---|---|---|
| All 13 analyses | `test_type` | `"DESeq2/NOISeq"` | `DESeq2` |

**Exact search/replace strings for Agent A:**

Search 1 (appears 13 times — use replace_all):
```
          test_type: "DESeq2/NOISeq"
```
Replace with:
```
          test_type: DESeq2
```

---

### Weissberg_2025/paperconfig.yaml

**Issues found:**

All `condition_type` values: `growth_state` — canonical.
All `test_type: DESeq2` — canonical.
All `organism` fields: `Alteromonas macleodii HOT1A3`, `Prochlorococcus MED4` — canonical.
All required fields present.

**Result: NO CHANGES NEEDED.**

---

## Files With No Changes Required

| File | Reason |
|---|---|
| `MIT9313_resources/paperconfig.yaml` | No publication block, no statistical_analyses |
| `Al-Hosani 2015/paperconfig.yaml` | All values already canonical |
| `Anjur 2025/paperconfig.yaml` | All values already canonical |
| `biller 2016/paperconfig.yaml` | All values already canonical |
| `he 2022/paperconfig.yaml` | All values already canonical |
| `Weissberg_2025/paperconfig.yaml` | All values already canonical |

---

## Complete Change Summary Table

| File | Category | Before | After | Count |
|---|---|---|---|---|
| Aharonovich 2016 | organism | `Alteromonas HOT1A3` | `Alteromonas macleodii HOT1A3` | 2 |
| Aharonovich 2016 | treatment_organism | `Prochlorococcus marinus MIT9313` | `Prochlorococcus MIT9313` | 2 |
| bagby 2015 | condition_type | `gas_shock` | `carbon_stress` | 4 |
| bagby 2015 | test_type | `Affymetrix microarray` | `microarray` | 4 |
| barreto 2022 | condition_type | `pco2` | `carbon_stress` | 2 |
| Biller 2018 | condition_type | `light_stress` | `darkness` | 1 |
| coe 2024 | condition_type | `dark_tolerance` | `darkness` | 1 |
| Fang 2019 | condition_type | `viral_lysis_products` | `viral` | 1 |
| Hennon 2017 | condition_type | `gas_shock` | `carbon_stress` | 2 |
| Hennon 2017 | organism (id_translation + analyses) | `Alteromonas EZ55` | `Alteromonas macleodii EZ55` | 3 |
| Lin 2015 | condition_type | `nutrient_stress` | `phosphorus_stress` | 1 |
| lindell 2007 | test_type | `"Affymetrix microarray"` | `microarray` | 8 |
| martiny 2006 | condition_type | `nutrient_stress` | `phosphorus_stress` | 1 |
| martiny 2006 | test_type | `"Affymetrix microarray with Cyber-T Bayesian t-test"` | `microarray_Cyber-T` | 9 |
| steglich 2006 | test_type | `"Affymetrix microarray with LPE test"` | `microarray_LPE` | 6 |
| tetu 2019 | condition_type | `plastic_leachate_stress` | `plastic_stress` | 2 |
| thompson 2011 | test_type | `"Affymetrix microarray with Bayesian t-test (Cyber-T)"` | `microarray_Cyber-T` | 8 |
| tolonen 2006 | condition_type | `nutrient_stress` | `nitrogen_stress` | 2 |
| tolonen 2006 | test_type | `"Affymetrix microarray with Goldenspike"` | `microarray_Goldenspike` | 16 |
| Read 2017 | condition_type | `nutrient_stress` | `nitrogen_stress` | 1 |
| Read 2017 | test_type | `"Rockhopper/DESeq"` | `Rockhopper` | 3 |
| Thompson 2016 | test_type | `"DESeq2/NOISeq"` | `DESeq2` | 13 |

**Total individual field changes: 91**

---

## Agent Ordering and Task Assignments

```
D (validator first)
  └─> A (paperconfig edits, after D done) + C (tests, can start with D)
        └─> G (code review)
              └─> F (validation run)
```

### Task D-1: Update validate_paperconfig.py (Agent D, FIRST)

**File:** `/home/osnat/github/multiomics_biocypher_kg/.claude/skills/paperconfig/validate_paperconfig.py`

**What to add** (insert as new constants near top of file, after the existing `VALID_TYPES` and `VALID_ID_TYPES` constants):

```python
CANONICAL_ORGANISM_NAMES = {
    "Prochlorococcus MED4",
    "Prochlorococcus AS9601",
    "Prochlorococcus MIT9301",
    "Prochlorococcus MIT9312",
    "Prochlorococcus MIT9313",
    "Prochlorococcus NATL1A",
    "Prochlorococcus NATL2A",
    "Prochlorococcus RSP50",
    "Synechococcus WH8102",
    "Synechococcus CC9311",
    "Alteromonas macleodii HOT1A3",
    "Alteromonas macleodii EZ55",
    "Alteromonas macleodii MIT1002",
    # Treatment-only organisms
    "Phage",
    "Marinobacter",
    "Thalassospira",
    "Pseudohoeflea",
    "Alteromonas",
}

CANONICAL_CONDITION_TYPES = {
    "growth_medium",
    "nitrogen_stress",
    "phosphorus_stress",
    "iron_stress",
    "salt_stress",
    "carbon_stress",
    "light_stress",
    "darkness",
    "plastic_stress",
    "viral",
    "coculture",
    "growth_state",
    "temperature_stress",
}

CANONICAL_TEST_TYPES = {
    "DESeq2",
    "DESeq",
    "edgeR",
    "Rockhopper",
    "microarray",
    "microarray_Cyber-T",
    "microarray_LPE",
    "microarray_Goldenspike",
}
```

**New validation logic to add in the `validate()` function:**

1. For each `environmental_conditions` entry, check `condition_type` against `CANONICAL_CONDITION_TYPES`. Error (not warning) if not in set. Print: `f"{env_key}: condition_type '{val}' not in canonical enum {sorted(CANONICAL_CONDITION_TYPES)}"`.

2. For each `statistical_analyses` entry:
   - Check `organism` value against `CANONICAL_ORGANISM_NAMES`. Error if not in set.
   - Check `treatment_organism` value (if present) against `CANONICAL_ORGANISM_NAMES`. Error if not in set.
   - Check `test_type` value against `CANONICAL_TEST_TYPES`. Error if not in set.
   - Check `treatment_condition` is present and non-empty. Error if absent or empty.

3. Add `"treatment_condition"` to `REQUIRED_ANALYSIS_FIELDS` list (currently at line 31).

**Insertion point for env condition check:** After the existing `env_conds` print statement (around line 134), add a loop:
```python
        for env_key, env_val in env_conds.items():
            ct = env_val.get("condition_type", "")
            if ct and ct not in CANONICAL_CONDITION_TYPES:
                errors.append(
                    f"environmental_conditions.{env_key}: condition_type '{ct}' "
                    f"not in canonical enum. Allowed: {sorted(CANONICAL_CONDITION_TYPES)}"
                )
```

**Insertion point for analysis organism check:** Inside the analysis loop (after the existing `atype` check, around line 238), add:
```python
            # Check organism name
            org = analysis.get("organism", "")
            if org and org not in CANONICAL_ORGANISM_NAMES:
                errors.append(
                    f"{aid}: organism '{org}' not in canonical organism names. "
                    f"Allowed: {sorted(CANONICAL_ORGANISM_NAMES)}"
                )

            # Check treatment_organism name if present
            treat_org = analysis.get("treatment_organism", "")
            if treat_org and treat_org not in CANONICAL_ORGANISM_NAMES:
                errors.append(
                    f"{aid}: treatment_organism '{treat_org}' not in canonical organism names. "
                    f"Allowed: {sorted(CANONICAL_ORGANISM_NAMES)}"
                )

            # Check test_type
            ttype = analysis.get("test_type", "")
            if ttype and ttype not in CANONICAL_TEST_TYPES:
                errors.append(
                    f"{aid}: test_type '{ttype}' not in canonical set. "
                    f"Allowed: {sorted(CANONICAL_TEST_TYPES)}"
                )

            # Check treatment_condition is non-empty
            tc = analysis.get("treatment_condition", "")
            if not tc:
                errors.append(f"{aid}: 'treatment_condition' is required and must be non-empty")
```

**Also add `"treatment_condition"` to `REQUIRED_ANALYSIS_FIELDS`:**

Current line 30-34:
```python
REQUIRED_ANALYSIS_FIELDS = [
    "type", "name", "id", "test_type",
    "control_condition", "treatment_condition",
    "organism", "name_col", "logfc_col",
]
```
`treatment_condition` is already in `REQUIRED_ANALYSIS_FIELDS` on line 32. Good — it is checked by the existing `for field in REQUIRED_ANALYSIS_FIELDS` loop. The new check above (checking non-empty) is in addition.

**Completion signal for D:** D is done when `uv run python .claude/skills/paperconfig/validate_paperconfig.py data/Prochlorococcus/papers_and_supp/Aharonovich\ 2016/paperconfig.yaml` reports the organism and condition_type errors for the pre-fix file (i.e., the validator actively catches them). A must verify this before starting.

---

### Task A-1: Standardize all paperconfig files (Agent A, after D-1)

**Order of files to process:** Process files with most changes first to build confidence; files with no changes last.

1. tolonen 2006 (18 changes — largest; good smoke test)
2. martiny 2006 (10 changes)
3. thompson 2011 (8 changes)
4. Thompson 2016 (13 changes)
5. lindell 2007 (8 changes)
6. steglich 2006 (6 changes)
7. bagby 2015 (8 changes)
8. barreto 2022 (2 changes)
9. Hennon 2017 (5 changes)
10. Aharonovich 2016 (4 changes)
11. Biller 2018 (1 change)
12. coe 2024 (1 change)
13. Fang 2019 (1 change)
14. Lin 2015 (1 change)
15. tetu 2019 (2 changes)
16. Read 2017 (4 changes)
17. Al-Hosani 2015 (0 changes — verify)
18. Anjur 2025 (0 changes — verify)
19. biller 2016 (0 changes — verify)
20. he 2022 (0 changes — verify)
21. Weissberg_2025 (0 changes — verify)

**After each file:** Run the validator to confirm:
```bash
uv run python .claude/skills/paperconfig/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/<Author Year>/paperconfig.yaml"
```
The validator must print `VALIDATION PASSED` with no ERRORS. Proceed to next file only after it passes.

**Validator run command template** (to be used after D finishes):
```bash
# Run validator on a single file
cd /home/osnat/github/multiomics_biocypher_kg
uv run python .claude/skills/paperconfig/validate_paperconfig.py "data/Prochlorococcus/papers_and_supp/tolonen 2006/paperconfig.yaml"
```

---

### Task C-1: Extend test_paperconfig_validation.py (Agent C, after D-1, parallel with A)

**File:** `/home/osnat/github/multiomics_biocypher_kg/tests/test_paperconfig_validation.py`

Read the current test file first. Then add or extend tests to assert:
1. All 23 active paperconfigs (excluding MIT9313_resources if it has no publication) pass the updated validator — this test likely already exists.
2. A synthetic paperconfig with `organism: "Alteromonas HOT1A3"` fails validation with an error mentioning the organism field.
3. A synthetic paperconfig with `condition_type: gas_shock` fails validation.
4. A synthetic paperconfig with `test_type: "Affymetrix microarray with Goldenspike"` fails validation.
5. A synthetic paperconfig missing `treatment_condition` in a statistical analysis fails validation.

C can write and run these tests immediately after D finishes — before A is done — because the tests use synthetic fixtures that do not require the actual paperconfig files to be already fixed.

---

## Test Commands and Expected Output

### Phase gate test suite

Run in order:

```bash
cd /home/osnat/github/multiomics_biocypher_kg

# 1. Validate each paperconfig individually (run after A finishes all files)
for f in \
  "data/Prochlorococcus/papers_and_supp/Aharonovich 2016/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/Al-Hosani 2015/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/Anjur 2025/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/bagby 2015/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/barreto 2022/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/biller 2016/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/Biller 2018/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/coe 2024/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/Fang 2019/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/he 2022/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/Hennon 2017/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/Lin 2015/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/lindell 2007/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/martiny 2006/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/steglich 2006/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/tetu 2019/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/thompson 2011/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/tolonen 2006/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/Read 2017/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/Thompson 2016/paperconfig.yaml" \
  "data/Prochlorococcus/papers_and_supp/Weissberg_2025/paperconfig.yaml"
do
  echo "=== $f ==="
  uv run python .claude/skills/paperconfig/validate_paperconfig.py "$f"
done
```
**Expected output for each file:** `VALIDATION PASSED` with 0 ERRORS.

```bash
# 2. Grep for old vocabulary (must return zero matches)
grep -r "gas_shock\|pco2\|dark_tolerance\|viral_lysis_products\|plastic_leachate_stress\|nutrient_stress\|salt_acclimation\|Affymetrix microarray\|Rockhopper/DESeq\|DESeq2/NOISeq" \
  data/Prochlorococcus/papers_and_supp/ \
  --include="paperconfig.yaml"
```
**Expected output:** No lines printed (grep returns exit code 1).

```bash
# 3. Grep for old organism names (must return zero matches)
grep -r "Alteromonas HOT1A3\|Alteromonas EZ55\|Prochlorococcus marinus" \
  data/Prochlorococcus/papers_and_supp/ \
  --include="paperconfig.yaml"
```
**Expected output:** No lines printed.

```bash
# 4. Run unit tests (excluding slow and kg)
pytest -m "not slow and not kg" -v
```
**Expected output:** All tests pass. Look specifically for `test_paperconfig_validation.py` — should show new vocabulary-specific tests passing.

```bash
# 5. Build smoke test
uv run python create_knowledge_graph.py --test
```
**Expected output:** Build completes without error. Expression edge counts should be unchanged from pre-Phase-1 baseline (no data loss from vocabulary normalization alone).

---

## Acceptance Criteria Checklist for Agent G

Agent G reviews after Agent A finishes all files and after `pytest -m "not slow and not kg"` passes.

### validate_paperconfig.py (Task D-1)
- [ ] `CANONICAL_ORGANISM_NAMES` constant defined with exactly the 18 names listed above
- [ ] `CANONICAL_CONDITION_TYPES` constant defined with exactly the 13 values listed above
- [ ] `CANONICAL_TEST_TYPES` constant defined with exactly the 8 values listed above
- [ ] Environmental condition `condition_type` checked against `CANONICAL_CONDITION_TYPES`; violation produces an ERROR (not a warning)
- [ ] Analysis `organism` field checked against `CANONICAL_ORGANISM_NAMES`; violation produces an ERROR
- [ ] Analysis `treatment_organism` field (if present) checked against `CANONICAL_ORGANISM_NAMES`; violation produces an ERROR
- [ ] Analysis `test_type` field checked against `CANONICAL_TEST_TYPES`; violation produces an ERROR
- [ ] Analysis `treatment_condition` checked for non-empty; missing or empty produces an ERROR
- [ ] Error messages include: field name, offending value, allowed values

### Paperconfig files (Task A-1)
- [ ] Zero occurrences of `gas_shock` in any paperconfig.yaml
- [ ] Zero occurrences of `pco2` as a condition_type in any paperconfig.yaml
- [ ] Zero occurrences of `dark_tolerance` in any paperconfig.yaml
- [ ] Zero occurrences of `viral_lysis_products` in any paperconfig.yaml
- [ ] Zero occurrences of `plastic_leachate_stress` in any paperconfig.yaml
- [ ] Zero occurrences of `nutrient_stress` in any paperconfig.yaml
- [ ] Zero occurrences of `salt_acclimation` in any paperconfig.yaml
- [ ] Zero occurrences of `Affymetrix microarray` in any paperconfig.yaml (with or without suffix)
- [ ] Zero occurrences of `Rockhopper/DESeq` in any paperconfig.yaml
- [ ] Zero occurrences of `DESeq2/NOISeq` in any paperconfig.yaml
- [ ] Zero occurrences of `Alteromonas HOT1A3` (without macleodii) in any paperconfig.yaml
- [ ] Zero occurrences of `Alteromonas EZ55` (without macleodii) in any paperconfig.yaml
- [ ] Zero occurrences of `Prochlorococcus marinus MIT9313` in any paperconfig.yaml
- [ ] All 23 paperconfigs pass `validate_paperconfig.py` with 0 ERRORS
- [ ] `Biller 2018/paperconfig.yaml` darkness condition has `condition_type: darkness` (not `light_stress`)
- [ ] `Hennon 2017/paperconfig.yaml` id_translation entry has `organism: "Alteromonas macleodii EZ55"`
- [ ] `Lin 2015/paperconfig.yaml` `lin_2015_p_limited` has `condition_type: phosphorus_stress`
- [ ] `tolonen 2006/paperconfig.yaml` nitrogen deprivation conditions have `condition_type: nitrogen_stress`
- [ ] `Read 2017/paperconfig.yaml` `read_2017_n_depleted` has `condition_type: nitrogen_stress`
- [ ] No other `condition_type` values changed (non-target files untouched)
- [ ] `biller 2016/paperconfig.yaml` `coculture` condition_type is still `coculture` (must NOT be changed)

### Tests (Task C-1)
- [ ] New tests in `test_paperconfig_validation.py` cover: bad organism name, bad condition_type, bad test_type, missing treatment_condition
- [ ] All new tests pass
- [ ] Existing tests still pass (no regressions)

### Build smoke test
- [ ] `uv run python create_knowledge_graph.py --test` completes without error
- [ ] Expression edge count in output CSVs is unchanged from pre-phase baseline

---

## Decisions Log for This Phase

**D-P1-1:** The `Biller 2018` extended darkness condition currently has `condition_type: light_stress`. The master plan mapping says `light_stress (extended darkness)` → `darkness`. However, the `steglich 2006` conditions also use `light_stress` for actual light stress experiments and those should stay as `light_stress`. The change to `Biller 2018` is surgical: only the `biller_2018_extended_darkness` key is changed, not the diel LD control (`biller_2018_diel_ld_control` stays as `growth_medium`).

**D-P1-2:** The `Hennon 2017` id_translation entry organism field (`Alteromonas EZ55`) must also be changed even though it is not inside `statistical_analyses`. The validator will check organism fields in `id_translation` blocks if/when D adds that check. For now, A must change it regardless to ensure consistency.

**D-P1-3:** Free-text strings in `treatment_condition`, `control_condition`, `name`, `description`, and `experimental_context` fields are NOT subject to the organism name enum — they are human-readable prose. Only the `organism` and `treatment_organism` YAML keys must use canonical names.

**D-P1-4:** `Alteromonas macleodii MIT1002` is used in biller 2016 and coe 2024. The canonical form is `Alteromonas macleodii MIT1002`. The agent spec for D lists it as `Alteromonas MIT1002` in the canonical set but the master plan and these files consistently use `Alteromonas macleodii MIT1002`. **Resolution: canonical form is `Alteromonas macleodii MIT1002`**. The validator `CANONICAL_ORGANISM_NAMES` set must use this form. Existing files already use it correctly.

**D-P1-5:** The `lindell 2007` paperconfig has no `environmental_conditions` block (Phase 2 will add one). For Phase 1 only the `test_type` values need fixing.

---

## Blockers

None at start of phase. If A discovers a paperconfig where the correct canonical mapping is ambiguous (e.g., a nutrient_stress condition that is both nitrogen and phosphorus), A must flag it to H before making the change. Do not guess — stop and report.

---

## Post-Phase Gate Record

*(To be filled in by Agent H after gate)*

- Gate result: PENDING
- Date: —
- Edge count before: —
- Edge count after: —
- Surprises: —
- Lessons learned: —
