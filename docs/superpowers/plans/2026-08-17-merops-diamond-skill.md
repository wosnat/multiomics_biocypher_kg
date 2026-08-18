# MEROPS-Diamond Skill — Implementation Plan

Spec: `docs/superpowers/specs/2026-08-17-merops-diamond-design.md`
Pattern: tcdb-diamond clone via `/add-a-tool`. **Completed 2026-08-17.**

- [x] Step 0 — intent captured (spec §Step 0)
- [x] Step 1 — SKILL.md skeleton at `.claude/skills/merops-diamond/SKILL.md`
- [x] Step 1 — install path: `--refresh-merops` flag in runner (Flavor B; scan lib + family.txt/clan.txt → `~/tools/MEROPS/DB/`)
- [x] Step 2 — install run on this host; `~/tools/MEROPS/DB/merops.dmnd` exists (5,009 seqs, 594 families in reference); One-Time Setup filled
- [x] Step 3 — `multiomics_kg/utils/merops_diamond.py` (header/subject parsing, tier policy, consensus collapse, family reference parsing, consumption helpers)
- [x] Step 3 — `tests/test_merops_diamond.py` (33 tests, all pass)
- [x] Step 3 — `.claude/skills/merops-diamond/run_merops_diamond.py` (uses `load_genome_rows` + `tool_calls_io`; flags `--strains --force --limit --threads --refresh-merops`)
- [x] Step 3 — QC fields in skill_summary (input_proteins, parse_failures, tier/consensus/catalytic-type/entry-type/best-hit-kind distributions, wallclock_s, tool_version)
- [x] Step 3 — SKILL.md data-free sections (What It Does, Output Schema, QC field docs, spot-check table, Phase 2)
- [x] Step 4 — smoke test MED4 `--limit 100` via documented Quick Start (4 calls, schema verified); full MED4 spot checks pass (ClpP → S14/SK, FtsH1 → M41/MA)
- [x] Step 5 — full batch: 42 strains, zero FAILED, zero parse failures; 4,194 proteins / 4,254 candidates; cross-strain spot checks pass (WH7803 FtsH1 → M41/MA, MIT1002 S8 peptidase → S08/SB + SignalP SP)
- [x] Step 5 — SKILL.md finalized (Observed batch results, cross-strain narrative incl. signalp secretion cross-check, calibrated thresholds)
- [x] Step 5 — per-strain calls.json + skill_summary.json committed
- [x] `.gitignore` rule: `cache/data/*/genomes/*/merops/*.limited_*`
- [x] Step 6 — registered in `add-a-strain/SKILL.md` step 3 (Wave 1, single-strain + multi-strain blocks)
- [x] CLAUDE.md — skills list updated
- [x] `pytest -m "not slow and not kg"` passes

## Notes / decisions recorded

- Inhibitor families kept + flagged (`entry_type: "inhibitor"`), per approved spec.
- E-value floor 0.001 — exact tcdb-diamond parity; tier truncation is the quality gate.
- Ragged-hierarchy deviation: tier-2 claims in subfamily-less families land at
  `level_kind: "merops_family"` with `tier: 2` (confidence band ≠ claim depth).
- Tier-1 nearly absent in practice (6/4,254) — scan library has few marine
  bacterial holotypes; MEROPS-diamond is a family-level breadth source.
- Phase 2 deferred: `/integrate-a-tool` + fresh spec (MeropsFamily ontology,
  `Gene_has_merops_family`, MEROPS-published interpro.txt/GO_annotation.txt bridges).
