# SignalP → KG integration (Phase 2) — design

**Date:** 2026-05-27
**Skill:** `integrate-a-tool signalp`
**Phase-1 runner:** `.claude/skills/signalp-run/` (SignalP 6.0 via Docker — already run on 40 strains)
**Track:** 3A flat ontology + scored edge, 1:1, **structural** — an exact mirror of the PSORTb
integration (`docs/kg-changes/psortb-extension.md`).

This supersedes the deferred "Phase 2" sketch in the `signalp-run` SKILL.md and refines the
worked example in `integrate-a-tool/references/decision-tree.md` (two corrections noted below).

---

## Step-0 four bullets

1. **What SignalP predicts + value space.** One call per protein: a signal-peptide *type*
   ∈ `{SP, LIPO, TAT, TATLIPO, PILIN, OTHER}` plus six per-class probabilities (∈[0,1], sum≈1)
   and a cleavage site (when a peptide is present). `OTHER` = "no signal peptide" (the no-signal
   sentinel, ~95% of proteins). The KG keeps the **five real types** as ontology nodes and skips
   `OTHER` — absence of a `Gene_has_signal_peptide_type` edge encodes "no signal peptide".
   - `SP` = Sec/SPI · `LIPO` = Sec/SPII (lipoprotein) · `TAT` = Tat/SPI · `TATLIPO` = Tat/SPII
     (Tat lipoprotein) · `PILIN` = Sec/SPIII (pilin-like).

2. **calls.json shape + join key + nulls/sentinels.**
   - ⚠️ **Normalization prerequisite.** `signalp-run` predates the `add-a-tool` calls.json
     convention; it writes raw `signalp/output.json` + `signalp/prediction_results.txt`, **not**
     `<strain>.signalp.calls.json`. **Decision:** add a `--normalize` mode to `signalp-run` that
     parses the already-present `prediction_results.txt` into a standard calls.json (no Docker
     re-run; the raw output exists for all 40 strains).
   - **Correction vs decision-tree (a):** by normalizing *first*, the merge reads clean scalar
     fields via plain `passthrough`/`float`/`integer` — exactly like PSORTb — so **no new
     `annotation_transforms.py` transform is needed** (`_tx_parse_cs_pos` is obviated; CS parsing
     happens in the normalizer's pure parser, `multiomics_kg/utils/signalp.py`).
   - **calls.json:** dict keyed by **RefSeq WP_ accession** (first whitespace token of the FASTA
     header = `gene_mapping.protein_id`), value
     `{signalp_type, probability, cleavage_site, cleavage_probability}`. `cleavage_site` /
     `cleavage_probability` are `null` for `OTHER` and for kept types with no reported CS.
     `OTHER` rows are kept verbatim (mirrors PSORTb keeping `Unknown`).

3. **Categorical-vs-property + functional-vs-structural verdicts.**
   - **Categorical → ONTOLOGY track (3A).** A small, closed, shared vocabulary (5 types) with a
     numeric score → ontology node + `Gene_has_<x>` edge carrying the score as an **edge property**
     (the second scored ontology edge, after PSORTb; schema property-block modeled on
     `Changes_expression_of`).
   - **Structural.** Signal-peptide presence is *where/how the gene's product is handled* (a
     physical trait), not *what it does* → **stays out** of `Gene.annotation_types`,
     `informative_annotation_types`, and the `annotation_quality` 8-bucket count. Optional
     denormalized `Gene.signal_peptide_type` routing string + `Gene.signal_peptide_type_count`
     signal only.

4. **Target labels / ids / levels.**
   - **Node:** label `SignalPeptideType` (schema input label `signal peptide type`),
     id **`signalp_<TYPE>`** (e.g. `signalp_SP`).
     - **Correction vs decision-tree (b):** the worked example wrote `signalp:<TYPE>` (colon), but
       `normalize_curie('signalp:SP')` returns `None` — `signalp` is **not** a registered
       bioregistry prefix, so the canonical `normalize_curie(...) or 'signalp_<call>'` idiom falls
       back to **underscore** (`signalp_SP`), identical to `psortb_OuterMembrane`. Tests + adapter
       use the underscore form.
     - Props: `name` (str, display), `signalp_id` (str, raw code), `level` (int, `0` — flat),
       `gene_count` (post-import), `organism_count` (post-import).
   - **Edge:** `Gene_has_signal_peptide_type` (Gene → SignalPeptideType), properties:
     `probability: float` (winning class's likelihood, ∈[0,1]),
     `cleavage_site: int` (residue before the cut; omitted when null),
     `cleavage_probability: float` (omitted when null). The full 6-likelihood vector is
     **intentionally not stored** (Step-0 decision).
   - **No hierarchy** (`<x>_is_a_<x>`) — flat, level 0.

---

## Vocabulary (5 kept nodes, skip `OTHER`)

| code | node id | display name | SignalP class |
|---|---|---|---|
| `SP` | `signalp_SP` | Signal peptide (Sec/SPI) | Sec substrate, SPase I |
| `LIPO` | `signalp_LIPO` | Lipoprotein signal peptide (Sec/SPII) | Sec substrate, SPase II |
| `TAT` | `signalp_TAT` | TAT signal peptide (Tat/SPI) | Tat substrate, SPase I |
| `TATLIPO` | `signalp_TATLIPO` | TAT lipoprotein signal peptide (Tat/SPII) | Tat substrate, SPase II |
| `PILIN` | `signalp_PILIN` | Pilin-like signal peptide (Sec/SPIII) | Sec substrate, SPase III |

`OTHER` is the sentinel — kept verbatim in `gene_annotations_merged.json`, never a node/edge.

---

## Data path (mirror PSORTb end-to-end)

```
signalp/prediction_results.txt  (raw, present for 40 strains)
   │  signalp-run --normalize  (pure parser: multiomics_kg/utils/signalp.py)
   ▼
<strain>.signalp.calls.json     (WP_-keyed; committed under git)
   │  load_signalp() in build_gene_annotations.py  (mirror load_psortb)
   ▼
gene_annotations_merged.json[locus_tag].{signalp_type, signalp_probability,
                                          signalp_cleavage_site, signalp_cleavage_probability}
   │  signalp_adapter.py (mirror psortb_adapter.py): SignalPeptideAdapter + MultiSignalPeptideAdapter
   ▼
SignalPeptideType nodes (5, level 0) + Gene_has_signal_peptide_type edges (1:1, scored)
   │  post-import.{sh,cypher}
   ▼
gene_count / organism_count on node; indexes; Gene.signal_peptide_type routing string;
rank_by_probability on edges
```

---

## Files touched

**New:**
- `multiomics_kg/utils/signalp.py` — vocab + `is_kept` + `display_name` + pure parser
  (`parse_cs_pos`, `parse_prediction_results`, `first_token`)
- `tests/test_signalp.py` — unit tests for the parser + vocab
- `multiomics_kg/adapters/signalp_adapter.py` — `SignalPeptideAdapter` + `MultiSignalPeptideAdapter`
- `tests/kg_validity/test_signalp.py` — live-graph assertions (mirror `test_psortb.py`)
- `docs/kg-changes/signalp-extension.md` — release notes
- `cache/data/*/genomes/*/signalp/<strain>.signalp.calls.json` + `.skill_summary.json` (40 strains, committed)

**Edited:**
- `.claude/skills/signalp-run/{run_signalp.py, SKILL.md}` — `--normalize` mode + docs
- `config/gene_annotations_config.yaml` — `signalp` source + 4 field rules
- `multiomics_kg/download/build_gene_annotations.py` — `load_signalp` + thread `sp` through
  `_get_raw`, all six `_resolve_*`, `build_wide`, `build_merged`, `process_strain`,
  `_compute_contributing_sources`
- `multiomics_kg/adapters/data_source_adapter.py` — `_name_for`/`_description_for` for `signalp`
- `tests/test_data_source_adapter.py` — `test_emits_five_nodes` → six + signalp provenance
- `tests/kg_validity/test_data_source.py` — node count 5 → 6 + signalp assertion
- `config/schema_config.yaml` — `signal peptide type` node + `gene to signal peptide type association` edge
- `create_knowledge_graph.py` — wire `MultiSignalPeptideAdapter`
- `scripts/post-import.sh` + `scripts/post-import.cypher` — indexes + rollups + routing string + rank (byte-identical Cypher)
- `tests/kg_validity/snapshot_data.json` — regenerate
- `.gitignore` — `cache/data/*/genomes/*/signalp/*.limited_*`
- `CLAUDE.md` — Key graph facts + Actual Neo4j labels + indexes
- `.claude/skills/add-a-tool/SKILL.md` + `.claude/skills/add-a-strain/SKILL.md` — hand-off + step-10 already cover psortb; confirm/extend

**Threading landmine (from recipe):** new `sp` arg on `build_wide`/`build_merged`/`_resolve_union`
must be **optional-defaulted** (tests call them positionally). No field uses `extract_first_match`,
so `annotation_helpers.py` is **not** touched.

---

## Out of scope (follow-ups)

- MCP / explorer tool surface (no `gene_details` field add, no new MCP tool) — separate task.
- Re-running SignalP (normalization reads existing raw output).
- Migrating `signalp-run` to `tool_calls_io.py` (deferred, non-blocking).

## Validation gate (definition of done)

`/omics-edge-snapshot` before/after (expression edges unchanged, new `Gene_has_signal_peptide_type`
edges appear) → Docker rebuild → post-import → `pytest -m kg` green → `snapshot_data.json`
regenerated → `pytest -m "not slow and not kg"` green.
