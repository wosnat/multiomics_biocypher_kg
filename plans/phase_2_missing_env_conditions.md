# Phase 2: Add Missing Environmental Conditions

**Status:** Ready for execution
**Depends on:** Phase 1 (paperconfig standardization) complete

---

## Summary

Three coculture/phage papers need `environmental_conditions` blocks to describe baseline experimental conditions. These blocks create EnvironmentalCondition nodes in the graph, giving LLMs context about the experimental setup even though the expression edges route through OrganismTaxon (coculture partner / phage).

### Key Decision: No `environmental_treatment_condition_id` on analyses

**Finding from code review** (`omics_adapter.py`):

1. **EnvironmentalCondition nodes are created from the `environmental_conditions` block alone** — the adapter iterates all entries and creates nodes regardless of whether analyses reference them.
2. **`environmental_treatment_condition_id`** on analyses controls **edge routing** — if present, the edge source becomes the EnvironmentalCondition node instead of OrganismTaxon.
3. **`environmental_control_condition_id`** is NOT used by the adapter at all — it is purely paperconfig metadata.

**Therefore:** For these three papers, analyses must NOT get `environmental_treatment_condition_id` because:
- The treatment IS the coculture organism (or phage), not an environmental condition
- Adding `environmental_treatment_condition_id` would incorrectly reroute expression edges from OrganismTaxon to EnvironmentalCondition
- The env conditions describe the **baseline** (growth medium, temperature, light), not the experimental perturbation

**No changes to analyses are needed.** Adding the `environmental_conditions` block is sufficient.

---

## Current State

| Paper | `environmental_conditions` block | Action needed |
|-------|----------------------------------|---------------|
| Aharonovich 2016 | YES (`aharonovich_2016_pro99`, `condition_type: growth_medium`) | None — already complete |
| Biller 2016 | YES (3 conditions: axenic, coculture×2) | None — already complete |
| Lindell 2007 | **NO** | **Add `environmental_conditions` block** |

---

## Task T2.1: Add `environmental_conditions` to Lindell 2007

**File:** `data/Prochlorococcus/papers_and_supp/lindell 2007/paperconfig.yaml`

**Insert after** `papermainpdf:` line, before `supplementary_materials:`.

**YAML to add:**
```yaml
  environmental_conditions:
    lindell_2007_pro99:
      condition_type: growth_medium
      name: Pro99 medium phage infection baseline
      medium: Pro99 medium
      temperature: 24C
      light_condition: continuous_light
      light_intensity: 10 umol photons m-2 s-1
      description: Prochlorococcus MED4 grown in Pro99 medium at 24C under continuous light (10 umol photons m-2 s-1). Baseline conditions for phage infection time-course.
```

**Scientific basis:** Lindell et al. 2007 (Nature 449:83-86) Methods: MED4 grown in Pro99 medium under continuous cool white light at ~10 µmol photons m⁻² s⁻¹, 24°C — standard MED4 conditions matching Aharonovich 2016.

**What NOT to change:**
- Do NOT add `environmental_control_condition_id` or `environmental_treatment_condition_id` to any analysis
- All 8 analyses keep `treatment_organism: "Phage"` and `treatment_taxid: 10239`

---

## Test Commands

```bash
# 1. Paperconfig validation
pytest tests/test_paperconfig_validation.py -v

# 2. Unit tests (no regressions)
pytest -m "not slow and not kg" -v

# 3. Test mode build
uv run python create_knowledge_graph.py --test
```

**Expected:** New EnvironmentalCondition node `lindell_2007_pro99` appears in build output. Expression edge counts unchanged.

---

## Acceptance Criteria for Agent G

- [ ] Lindell 2007 has `environmental_conditions.lindell_2007_pro99` with `condition_type: growth_medium`
- [ ] No analyses in Lindell 2007 have `environmental_treatment_condition_id` or `environmental_control_condition_id`
- [ ] All 8 Lindell 2007 analyses still have `treatment_organism: "Phage"` and `treatment_taxid: 10239`
- [ ] Aharonovich 2016 and Biller 2016 are unchanged
- [ ] All tests pass
- [ ] Build succeeds with new EnvironmentalCondition node
