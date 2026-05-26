# PSORTb subcellular-localization — KG integration (Phase 2)

Phase-2 partner to the Phase-1 artifact spec
[`2026-05-10-psortb-localization-design.md`](2026-05-10-psortb-localization-design.md).
Wires the committed `<strain>.psortb.calls.json` artifacts into the live KG via the
`integrate-a-tool` skill. **This doc supersedes §6 of the Phase-1 spec**, which proposed a
Gene-property approach; the `integrate-a-tool` decision tree routes PSORTb to the ONTOLOGY
track instead (closed shared vocabulary with cross-gene grouping semantics).

> Scope of THIS task (dogfooding run): Step 0–3A only. Step 4 (post-import), Step 5 (Docker
> rebuild + kg-validity), Step 6 (release notes), Step 7 (hand-off) are out of scope here and
> are listed in the plan for completeness only.

## Step 0 — intent (the four bullets the skill asks for)

1. **What the tool predicts + value space.** PSORTb v3.0.3 predicts bacterial subcellular
   localization per protein. Value space = closed 6-value vocabulary:
   `Cytoplasmic`, `CytoplasmicMembrane`, `Periplasmic`, `OuterMembrane`, `Extracellular`,
   `Unknown` (no-signal sentinel). Each call carries a confidence `score`.

2. **calls.json shape + join key + nulls/sentinels + scores.**
   - File: `cache/data/<org>/genomes/<strain>/psortb/<strain>.psortb.calls.json`.
   - Top-level dict **keyed by WP_ accession** (RefSeq `protein_id`) → join key is `protein_id`
     (same as eggNOG).
   - Per-record fields: `localization` (str), `score` (float), `secondary_localization`
     (always null in Phase 1), `secondary_score` (always null), `is_multi_localized`
     (always false → strictly 1:1), `is_unknown` (bool; true ⟺ `localization == "Unknown"`).
   - Observed scores: real (non-Unknown) calls 7.88–10.0; `Unknown` calls 2.0–6.58. Decision-tree
     reference quotes the full record range 2.0–10.0; the *kept* (non-Unknown) edge `score`
     is effectively ≥ 7.5. (NB: psortb-run SKILL.md says `score ∈ [0,10]` — the observed floor
     is 2.0, not 0.0; harmless for our `score: float` field but noted.)
   - Sentinel: skip `is_unknown` / `localization == "Unknown"` rows — absence of an edge encodes
     "no confident localization".

3. **Categorical-vs-property + functional-vs-structural verdicts.**
   - **Categorical → ONTOLOGY TRACK (3A).** Closed 5-real-value vocabulary shared across genes;
     a node de-duplicates the vocab and makes "all OuterMembrane genes in organism X" one hop.
   - **Scored** → the `score` rides as an **edge property** on `Gene_has_subcellular_localization`
     (the first scored ontology edge — borrow the `properties:` + rank pattern from
     `Changes_expression_of`).
   - **Structural** (where the gene *is*, not what it does) → do **NOT** fold into
     `Gene.annotation_types` / `informative_annotation_types` / `annotation_quality`. (Step 4.)

4. **Target labels + id prefix + cardinality.**
   - Node label `SubcellularLocalization` (BioCypher input phrase `subcellular localization`).
   - Node id CURIE prefix `psortb:` → `psortb:OuterMembrane` etc.
   - Flat ontology: `level = 0` on every node, no `level_kind`, no `<x>_is_a_<x>` hierarchy edge.
   - Edge `Gene_has_subcellular_localization` (input phrase `gene_has_subcellular_localization`),
     `properties: { score: float }`. 1:1 → ≤1 edge per gene.

## Checkbox plan

### Step 0 — spec + plan
- [x] This doc.

### Step 1 — inspect calls.json
- [x] `jq` the committed MED4 + KT2440 calls.json: 1872 MED4 records, WP_ keys, 6-value vocab,
      `is_multi_localized` always false, `secondary_*` always null, score 2.0–10.0.
- [x] Confirm no `multiomics_kg/utils/psortb.py` / `tool_calls_io.py` exists yet → write a pure
      vocab/parse helper.

### Step 2 — gene-annotation merge (front door)
- [ ] `config/gene_annotations_config.yaml`: `sources.psortb` + `fields.psortb_localization` +
      `fields.psortb_score` + `logical_sources: psortb`.
- [ ] `build_gene_annotations.py`: `load_psortb()` + thread `ps` through `_get_raw` /
      `_resolve_*` / `build_wide` / `build_merged` / `process_strain` + `_compute_contributing_sources`.
- [ ] Verify: `bash scripts/prepare_data.sh --steps 2 --strains MED4 --force`; field landed;
      `contributing_sources` lists `psortb`.

### Step 3A — ontology track
- [ ] `multiomics_kg/utils/psortb.py`: `LOCALIZATION_VOCAB` + pure helpers.
- [ ] `config/schema_config.yaml`: `subcellular localization` node (`level: int`) +
      `gene to subcellular localization association` edge with `properties: { score: float }`.
- [ ] `multiomics_kg/adapters/psortb_adapter.py`: `SubcellularLocalizationAdapter` (per-strain,
      1:1 scored edge) + `MultiSubcellularLocalizationAdapter` (owns 5 nodes, level 0, flat).
- [ ] Wire `MultiSubcellularLocalizationAdapter` into `create_knowledge_graph.py`.
- [ ] `tests/test_psortb.py`: pure vocab/parse + adapter edge/node unit tests.
- [ ] `pytest -m "not slow and not kg" tests/test_psortb.py -v`.

### Step 4–7 — OUT OF SCOPE for this run
- [ ] Step 4 post-import rollups + indexes (flat `gene_count`/`organism_count`, indexes, optional
      `subcellular_localization` Gene routing string + `rank_by_score`).
- [ ] Step 5 Docker rebuild + `/omics-edge-snapshot` + `pytest -m kg` + snapshot.
- [ ] Step 6 `docs/kg-changes/psortb-extension.md` + CLAUDE.md.
- [ ] Step 7 add-a-tool / add-a-strain hand-off + MCP follow-up.
