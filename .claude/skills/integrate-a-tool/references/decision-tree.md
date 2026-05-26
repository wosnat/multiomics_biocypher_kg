# Decision tree — routing a tool to a track

The operational routing reference for Step 0 — two orthogonal decisions, both
made by inspecting the committed calls.json with `jq`.

## Decision 1 — categorical-vs-property (which track)

```
├─ Small controlled vocabulary (≲ a few hundred distinct values), shared across
│  genes, with hierarchy or cross-gene grouping semantics?
│        e.g. PSORTb localization {Cytoplasmic, OuterMembrane, …},
│             SignalP type {SP, LIPO, TAT, PILIN}, a TC/CAZy-style class tree
│   → ONTOLOGY TRACK (3A): new node type + Gene_has_<x> edges
│     (+ <x>_is_a_<x> hierarchy edges if there are levels). Copy cazy_adapter.
│
├─ …and the categorical call ALSO carries a numeric score / probability / tier?
│        e.g. PSORTb localization + score, SignalP type + probability,
│             tcdb-diamond TC family + identity/qcov + tier
│   → STILL ONTOLOGY TRACK (3A); the score rides as an EDGE PROPERTY on
│     Gene_has_<x>. The node de-duplicates the vocabulary; per-gene evidence
│     strength lives on the edge, exactly like Changes_expression_of carries
│     log2_fold_change + adjusted_p_value. No existing ontology edge carries
│     properties — a scored tool is the first, so borrow the schema +
│     post-import rank pattern from Changes_expression_of, NOT from an ontology.
│
└─ Per-gene numeric / scalar / continuous value, OR free-text, OR a single
   boolean flag, with NO shared cross-gene vocabulary?
        e.g. a standalone probability, a count, GC%, a yes/no with no taxonomy
   → GENE-PROPERTY TRACK (3B): new property on the Gene node in schema_config;
     cyanorak_ncbi_adapter reads it (after a GeneNodeField enum + int_fields/
     float_fields edit — the read is NOT automatic; see gene-annotation-merge-recipe.md).
```

**Multiple calls per gene is the ontology default.** A gene may carry many calls
(several TC families, several CAZy families, several localizations). The merged
field is a **list**; the adapter emits **one `Gene_has_<x>` edge per call**, fanning
one gene out to N edges pointing at N shared, de-duplicated nodes. 1:1 tools
(psortb, signalp) simply emit ≤1 edge per gene.

### Why a closed 5-value vocab is still an ontology, not a Gene property

Localization and SP-type are *shared, closed vocabularies*. A node makes
"all OuterMembrane genes in organism X" a one-hop `MATCH` and de-duplicates the
vocabulary across the whole graph. The 1:1 cardinality doesn't change that.
A Gene property (3B) is right only for a per-gene scalar with **no** shared
vocabulary (a GC%, a codon-usage score).

## Decision 2 — functional-vs-structural (which cross-node rollups)

This decision is independent of the track and drives whether Step 4 touches Gene
rollups at all (see post-import-patterns.md → "Cross-node rollups").

| | Functional | Structural |
|---|---|---|
| What it asserts | *what the gene does* — KO-/domain-/role-like evidence | *where the gene is / a physical trait* — localization, signal-peptide |
| Examples | a new EC-like or domain-like predictor | psortb localization, signalp signal-peptide |
| Fold into `Gene.annotation_types` + `informative_annotation_types`? | **yes** | **no** (would inflate functional-annotation coverage + skew `genes_by_function` `min_quality`) |
| Enter the `annotation_quality` 8-bucket count? | maybe — decide; follow the bucket-maintenance procedure | no |
| Optional `<tool>_count` Gene signal | optional | optional (surfaces via `gene_details` for free) |

**Both immediate targets (psortb, signalp) are STRUCTURAL** → they stay out of the
annotation arrays.

## Worked examples

Both are **3A flat ontology + edge score, 1:1, skip-the-no-signal-sentinel**.
Copy these when implementing either.

### PSORTb → `SubcellularLocalization` (no artifact prerequisite)

- **Phase-1 input:** `cache/data/<org>/genomes/<strain>/psortb/<strain>.psortb.calls.json`
  — dict **keyed by WP_ accession**; record
  `{localization, score (2.0–10.0), is_unknown, secondary_localization, secondary_score}`.
  `is_multi_localized` is always false → strictly 1:1.
- **Vocabulary (5 nodes, skip `Unknown`):** `Cytoplasmic`, `CytoplasmicMembrane`,
  `Periplasmic`, `OuterMembrane`, `Extracellular`.
- **Node:** label `SubcellularLocalization`, id `psortb:<Class>` (e.g.
  `psortb:OuterMembrane`); props `name`, `psortb_id`, `level = 0` (flat — omit
  `level_kind`).
- **Edge:** `Gene_has_subcellular_localization` (Gene → node),
  **`properties: { score: float }`** — the first scored ontology edge; model the
  schema `properties:` block on `Changes_expression_of`.
- **Merge (Step 2):** `join_to: protein_id` (WP_, exactly like eggNOG). Two scalar
  fields via `passthrough`: `psortb_localization: str`, `psortb_score: float`;
  skip rows where `is_unknown`.
- **Adapter (Step 3A):** `SubcellularLocalizationAdapter` (per-strain — reads the
  two fields, emits ≤1 scored edge/gene) + `MultiSubcellularLocalizationAdapter`
  (owns the 5 nodes at `level=0`). Copy cazy_adapter, then **delete the fan-out
  loop and the `*_is_a_*` hierarchy** and **add `score` to the edge tuple**.
- **Post-import (Step 4):** flat `gene_count`/`organism_count` (no `*0..`); scalar
  indexes on `level` + `psortb_id`, full-text on `name`; **optional** per-class
  `rank_by_score`. Do **NOT** fold into `annotation_types`/`annotation_quality`
  (structural); optionally denormalize a `subcellular_localization` routing string.
- **kg-validity (Step 5):** 5 nodes; no `Unknown` node; edge `score ∈ [2.0, 10.0]`;
  no orphan gene→node edges; rollup sanity.
- **Release notes:** `docs/kg-changes/psortb-extension.md`.

### SignalP → `SignalPeptideType`  ⚠️ calls.json prerequisite

- **Fix the artifact first:** `signalp-run` predates the calls.json convention and
  writes raw `signalp/output.json` (envelope + `.SEQUENCES` keyed by the **full
  FASTA header**), *not* a normalized `<strain>.signalp.calls.json`. Step-0 decision:
  either extend `signalp-run` to emit the standard calls.json, or have
  `load_signalp()` parse `output.json` directly.
- **Per-record:** `Prediction` (winning type), `Likelihood[6]` (∈[0,1], parallel to
  `Protein_types[6]`), `CS_pos` (e.g. `"Cleavage site between pos. 26 and 27.
  Probability 0.862107"`, `""` when none).
- **Key:** WP_ accession = first whitespace token of the header → reuse the existing
  `_tx_first_token_space` transform.
- **Vocabulary (5 nodes, skip `Other`):** `SP` (Sec/SPI), `LIPO` (Sec/SPII),
  `TAT` (Tat/SPI), `TATLIPO` (Tat/SPII), `PILIN` (Sec/SPIII).
- **Node:** label `SignalPeptideType`, id `signalp:<TYPE>` (e.g. `signalp:SP`);
  props `name`, `signalp_id`, `level = 0`.
- **Edge:** `Gene_has_signal_peptide_type`,
  `properties: { probability: float, cleavage_site: int|null, cleavage_probability: float|null }`
  (`probability` = the winning class's `Likelihood`). The full 6-vector is
  intentionally **not** stored — record that as a Step-0 decision.
- **Merge (Step 2):** `join_to: protein_id` via the header-split transform. Fields:
  `signalp_type: str`, `signalp_probability: float`, `signalp_cleavage_site: int`,
  `signalp_cleavage_probability: float`; skip `Prediction == Other`. New transform
  `_tx_parse_cs_pos` (parse `CS_pos` → position + prob), registered in `_TRANSFORMS`,
  unit-tested in `tests/test_annotation_transforms.py`.
- **Adapter / post-import / kg-validity:** identical shape to PSORTb. Validity:
  `probability ∈ [0,1]`; `cleavage_site` a positive int when present.
- **Release notes:** `docs/kg-changes/signalp-extension.md`.
