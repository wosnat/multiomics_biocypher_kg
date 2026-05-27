# SignalP KG integration — checkbox plan

Spec: `docs/superpowers/specs/2026-05-27-signalp-kg-integration-design.md`. Mirror PSORTb.

## Step 1 — normalize raw output → calls.json
- [ ] `multiomics_kg/utils/signalp.py`: `SIGNALP_VOCAB`, `OTHER_SENTINEL`, `is_kept`, `display_name`, `first_token`, `parse_cs_pos`, `parse_prediction_results`
- [ ] `tests/test_signalp.py`: parser + vocab unit tests (green)
- [ ] `signalp-run/run_signalp.py`: `--normalize` mode (reads `prediction_results.txt`, no Docker) → `<strain>.signalp.calls.json` + `.skill_summary.json`
- [ ] Run `--normalize` on all 40 strains; spot-check MED4 (`WP_011131644.1` → SP); commit calls.json
- [ ] `signalp-run/SKILL.md`: document `--normalize` + Output Schema (calls.json)

## Step 2 — gene-annotation merge
- [ ] `gene_annotations_config.yaml`: `signalp` source (logical_sources) + 4 field rules (signalp_type/probability/cleavage_site/cleavage_probability)
- [ ] `build_gene_annotations.py`: `load_signalp` + thread `sp` (optional-defaulted) through `_get_raw` + all six `_resolve_*` + `build_wide` + `build_merged` + `process_strain` row-join + `_compute_contributing_sources`
- [ ] `data_source_adapter.py`: `_name_for`/`_description_for` for signalp
- [ ] `tests/test_data_source_adapter.py`: five→six + signalp provenance
- [ ] Verify: `prepare_data.sh --steps 2 --strains MED4 --force`; jq field landed + contributing_sources; `pytest tests/test_build_gene_annotations.py tests/test_data_source_adapter.py` green

## Step 3 — 3A ontology track
- [ ] `schema_config.yaml`: `signal peptide type` node + `gene to signal peptide type association` edge (properties: probability/cleavage_site/cleavage_probability)
- [ ] `signalp_adapter.py`: `SignalPeptideAdapter` + `MultiSignalPeptideAdapter` (copy psortb_adapter; underscore id; score→edge props)
- [ ] wire `MultiSignalPeptideAdapter` into `create_knowledge_graph.py`
- [ ] `tests/test_signalp.py`: adapter node/edge tests

## Step 4 — post-import
- [ ] `post-import.sh` + `.cypher` (byte-identical): indexes (level, signalp_id, fulltext name); gene_count/organism_count; Gene.signal_peptide_type routing string; rank_by_probability on edges
- [ ] `post-import-validate.sh` baseline before rebuild

## Step 5 — live KG validation
- [ ] `/omics-edge-snapshot` before
- [ ] `docker compose up -d --build` + post-import
- [ ] `/omics-edge-snapshot` after (expression unchanged; new edges appear)
- [ ] `tests/kg_validity/test_signalp.py` + `test_data_source.py` count 5→6
- [ ] `pytest -m kg` green; regenerate `snapshot_data.json`

## Step 6 — docs
- [ ] `docs/kg-changes/signalp-extension.md`
- [ ] CLAUDE.md (Key graph facts + labels + indexes)

## Step 7 — hand-off
- [ ] confirm add-a-tool Phase-2 redirect + add-a-strain step-3/step-10 cover signalp
- [ ] `.gitignore` `*.limited_*` rule; MCP follow-up recorded
- [ ] `pytest -m "not slow and not kg"` green
