# InterProScan → KG integration plan (Phase 2)

Spec: `docs/superpowers/specs/2026-07-26-interproscan-kg-integration-design.md`
Track: 3A ontology (InterproEntry) · functional · hierarchical · scored edge.

## Step 0 — spec + plan
- [x] Inspect calls.json (join key WP_, GO/pathway empty, unintegrated %, entry-less gene count)
- [x] Route: ontology / functional / entries-only / hierarchy-in / GO-deferred (user-confirmed)
- [x] Write spec + this plan

## Step 1 — inspect runner & reuse parse utils
- [ ] Read `/interproscan-run` SKILL.md + `multiomics_kg/utils/interproscan.py` (reuse parsing)
- [ ] Confirm distinct interpro_type values + edge score/evalue nullability across a heterotroph

# ══ PHASE A — cache/data (all downloads via prepare_data.sh) ══

## Step 2a — InterPro reference download → prepare_data step 9
- [ ] `multiomics_kg/utils/interpro_reference.py` — pure parse of entry.list + ParentChildTreeFile → `{id:{name,type,parent,level}}`
- [ ] `multiomics_kg/download/build_interpro_reference.py` — download + I/O + cache `cache/data/interpro/interpro_reference.json` (committed); `--force`/`--refetch-raw`
- [ ] `scripts/prepare_data.sh` — add step 9 (case + default STEPS + header comment + usage) ; `tests/test_interpro_reference.py` (pure parse)

## Step 2b — gene-annotation merge (interproscan source)
- [ ] `config/gene_annotations_config.yaml`: `sources:` + `fields:` (`interpro_entries`) + `logical_sources:` for interproscan
- [ ] `build_gene_annotations.py`: `load_interproscan()` + thread defaulted `ip` arg through `_get_raw` → all six `_resolve_*` → `build_wide`/`build_merged` → `process_strain` row-level `protein_id` join + `_compute_contributing_sources`
- [ ] New source arg **optional-defaulted** (tests call positionally)
- [ ] DataSource: `data_source_adapter.py` `_name_for`/`_description_for` + bump node count in `test_data_source_adapter.py` (unit) **and** `tests/kg_validity/test_data_source.py` (live)

## Step 2c — regenerate caches + QA + COMMIT (gate)
- [ ] `prepare_data.sh --steps 9 --force` → reference JSON exists + sane counts
- [ ] `prepare_data.sh --steps 2 --strains MED4 --force` → `interpro_entries` lands + `contributing_sources` lists interproscan; DataSource node has real name
- [ ] `prepare_data.sh --steps 2 --force` (all 42 strains) — regenerate merged JSON caches
- [ ] `pytest tests/test_build_gene_annotations.py tests/test_data_source_adapter.py tests/test_interpro_reference.py -q` green
- [ ] **QA the cache diff, then COMMIT** (reference JSON + merged-JSON caches + config + build_gene_annotations + prepare_data + tests) — checkpoint before adapter work

# ══ PHASE B — adapter/schema/graph ══

## Step 3 — ontology track (3A)
- [ ] `config/schema_config.yaml` — `interpro entry` node (+ `level`, `interpro_type`) + `gene_has_interpro_entry` edge (properties block) + `interpro_entry_is_a_interpro_entry` + `pfam_in_interpro_entry` bridge
- [ ] `multiomics_kg/adapters/interpro_adapter.py` — `InterproAnnotationAdapter` (per-strain edges w/ aggregated evidence, reads calls.json) + `MultiInterproAnnotationAdapter` (nodes + hierarchy pruned to observed+ancestors + Pfam bridge pruned to injected Pfam node set, reads interpro_reference.json)
- [ ] Wire `MultiInterproAnnotationAdapter` into `create_knowledge_graph.py` (pass Pfam node set for bridge pruning, BRITE-`known_ko_ids` style)
- [ ] `tests/test_interpro_adapter.py` (pure parsing/aggregation/hierarchy)

## Step 4 — post-import rollups + indexes
- [ ] `post-import.sh` + `post-import.cypher` (byte-identical): InterproEntry `gene_count`/`organism_count`/`member_count` + `is_promiscuous` (ORA support, §8); indexes (`level`, `interpro_id`, `interpro_type`, full-text `name`); Gene `annotation_types` += interpro + `interpro_entry_count`
- [ ] `post-import-validate.sh` baseline/after = additions-only

## Step 5 — live KG validation gate
- [ ] `/omics-edge-snapshot` before → `docker compose up -d --build` → post-import → `/omics-edge-snapshot` after
- [ ] `tests/kg_validity/test_interpro.py`; `pytest -m kg`; regenerate `snapshot_data.json`

## Step 6 — docs (MCP-integration doc + release notes)
- [ ] `docs/kg-changes/interproscan-extension.md` — written **for the MCP/explorer to know how to integrate**: node/edge labels + full property list, id conventions, `interpro_type` vocab + `level` semantics, indexes/full-text, the `(interpro_type, level)` ORA stratification guidance + `is_promiscuous`, Pfam bridge, counts, validation results, out-of-scope
- [ ] CLAUDE.md — Key graph facts + Actual Neo4j labels + Data Locations (interpro_reference.json) + prepare_data step 9 + link the doc
- [ ] Update release notes (CHANGELOG.md entry for the InterPro extension)

## Step 7 — hand-off
- [ ] add-a-tool Phase-2 hand-off redirects here; add-a-strain step-10 covers merged field; MCP follow-up recorded
