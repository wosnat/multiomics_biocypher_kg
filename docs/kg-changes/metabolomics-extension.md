# Metabolomics paper integration (Phase 2)

**Spec:** [docs/superpowers/specs/2026-05-03-metabolomics-paper-integration-design.md](../superpowers/specs/2026-05-03-metabolomics-paper-integration-design.md)
**Plan:** [docs/superpowers/plans/2026-05-03-metabolomics-paper-integration.md](../superpowers/plans/2026-05-03-metabolomics-paper-integration.md)
**Validation papers:** Capovilla 2023, Kujawinski 2023
**Lands as of 2026-05-04**

## What changed

### New node + edge types
- **`MetaboliteAssay`** — one per (Experiment × `value_kind`). Carries assay-level metabolomics semantics. Adapter: [metabolite_assay_adapter.py](../../multiomics_kg/adapters/metabolite_assay_adapter.py) (mirrors `observations_adapter.py`). Node IDs: `metabolite_assay:{doi_short}:{entry_key}:{metric_type}`.
- **`Assay_quantifies_metabolite`** — numeric edge. Replicate-aggregated. Carries `value`, `value_sd`, `n_replicates`, `n_non_zero`, `replicate_values`, adapter-set `detection_status ∈ {detected, sporadic, not_detected}`. Post-import: `rank_by_metric`, `metric_percentile`, `metric_bucket` (only when parent assay `rankable: "true"`).
- **`Assay_flags_metabolite`** — boolean edge. Carries `flag_value`, `n_replicates`, `n_positive`.
- Binding edges: `PublicationHasMetaboliteAssay`, `ExperimentHasMetaboliteAssay`, `MetaboliteAssayBelongsToOrganism` (PascalCase form mirrors existing DerivedMetric bindings).

### Property additions on existing types

| Type | Property | Source | Semantics |
|---|---|---|---|
| `Metabolite` | `evidence_sources` | adapter (extended) | gains `"metabolomics"` value when paper-measured |
| `Metabolite` | `measured_assay_count`, `measured_organisms` (str[]), `measured_compartments` (str[]), `measured_paper_count` | post-import | distinct counts / sorted distinct values across incoming Assay edges; `[]` on non-measured (KG-MET-016, 2026-05-05) |
| `TcdbFamily` | `is_promiscuous` (bool) | post-import | true when `metabolite_count >= 50 OR member_count >= 100`; flags ~30 of 12,883 families in the long tail; consumed by explorer family_inferred-dominance warnings (KG-MET-006, 2026-05-05) |
| `Metabolite` | `organism_count`, `organism_names` | post-import | UNION extended to include the measurement path |
| `Organism_has_metabolite` | `evidence_sources` (str[]) | post-import | values: `metabolism` \| `transport` \| `measured` |
| `Organism_has_metabolite` | `measured_assay_count`, `measured_compartments`, `measured_paper_count` | post-import | 0/[] on non-measured edges |
| `Experiment` | `metabolite_assay_count`, `metabolite_compartments`, `metabolite_count` | post-import | rollups |
| `Publication` | `metabolite_assay_count`, `metabolite_compartments`, `metabolite_count` | post-import | rollups |
| `OrganismTaxon` | `measured_metabolite_count` | post-import | rollup |

`Organism_has_metabolite` materialization extends to the measurement-only path: an organism gains a metabolite edge if any MetaboliteAssay anchored on that organism quantifies/flags the metabolite, even with no gene-catalysis or gene-transport path.

### Indexes added (post-import)

Scalar: `metabolite_assay_organism_idx`, `metabolite_assay_compartment_idx`, `metabolite_assay_metric_type_idx`, `metabolite_assay_value_kind_idx`, `metabolite_assay_experiment_idx`. Full-text: `metaboliteAssayFullText` on `name`, `field_description`, `treatment`, `experimental_context`.

### Conventions

- **`MetaboliteAssay.field_description` is the canonical normalisation-provenance field** for the metabolomics layer. Carries paper-specific text covering units, blank-correction, replicate-aggregation policy, and source-table provenance (e.g. *"Intracellular metabolite concentration in fg/cell, blank-corrected, replicate-aggregated; Capovilla 2023 Table sd03"*). Companion fields `value_kind`, `unit`, `metric_type`, `aggregation_method` round out the provenance picture. Indexed by `metaboliteAssayFullText`.
- **`Metabolite` nodes are compartment-agnostic.** Compartment lives on `MetaboliteAssay.compartment` and on `Assay_*_metabolite` edges' parent assay; the same `kegg.compound:` ID is reused across compartments. 92 metabolites are measured in two compartments (e.g. whole_cell + extracellular) against a single `Metabolite` node. Per-pair compartment rollup lives on `Organism_has_metabolite.measured_compartments`.

### paperconfig + pipeline

- New entry type `metabolite_assays_table` (mirrors `derived_metrics_table`). See [.claude/skills/paperconfig/SKILL.md](../../.claude/skills/paperconfig/SKILL.md).
- Per-paper `metabolite_aliases.yaml` for free-text-name → primary-ID overrides (winning over MNX hits).
- Step 6 ([build_kegg_metabolism_xrefs.py](../../multiomics_kg/download/build_kegg_metabolism_xrefs.py)) extended to harvest paper metabolites: extends `kegg_data.json` (with `evidence_sources` including `"metabolomics"`) and writes a new `cache/data/metabolomics/metabolite_id_mapping.json`.
- New step 7 ([resolve_paper_metabolites.py](../../multiomics_kg/download/resolve_paper_metabolites.py)): pure CSV rewriter — writes `<stem>_resolved.csv` for each `metabolite_assays_table` source with `metabolite_id` + `resolution_method` columns. Mirrors gene-side step 4. Does not load the MNX resolver.
- Step 6 perf split: `--force` rebuilds outputs from cached raw inputs (~1 min, no network); `--refetch-raw` also re-pulls raw KEGG REST + TCDB TSVs (use only on upstream releases). Frees the alias-iteration loop from network round-trips.

### Vocabulary

- Compartment vocab adds `extracellular` (was already in `METABOLITE_COMPARTMENT_VOCAB` for the metabolite_assays_table validator; now also in the central `COMPARTMENTS` frozenset used by Experiment validation, fixing a hole). Full set: `whole_cell, vesicle, exoproteome, extracellular, spent_medium, lysate`.
- `MIT0801` (Prochlorococcus, LLI ecotype) added to validator's `CANONICAL_GENOMIC_ORGANISMS` (was deployed in commit `9e6e183` but never registered with the validator).

## Validation results

| Paper | Source | Resolved / Total | Method mix |
|---|---|---|---|
| Capovilla 2023 | sd03 cellular concentrations | 16 / 16 (100%) | 10 `name_match` + 6 `alias_override` |
| Kujawinski 2023 | ChisholmPro KEGG export | 651 / 651 (100%, 7 entries × 93 rows) | `kegg_direct` (KEGG ID column) |
| Kujawinski 2023 | Table S2 presence flags | 93 / 93 (100%) | `name_match` + `alias_override` |

Both papers exceed spec §5.1 acceptance criterion A4 (`≥ 70%` Capovilla, `≥ 95%` Kujawinski).

## Known follow-ups

- biller 2022 metabolomics tables (S5–S7): partial. MS-feature-only Table S6 (4065 rows without compound IDs) is out of scope for v1 — would need a `MassFeature` node type. Embedded-uncertainty cells in S7 are in scope via `cell_format: embedded_mean_sd_n` parser; the entry shape is reserved.
- Lipid-class names (LOBSTAHS / LipidMaps notation in biller 2022 S1) often unresolved by MNX — future LipidMaps integration. `compound_class: lipid` hint reserved on the entry shape.
- Three-tier `metabolite_id_mapping.json` schema currently only fills `name_lookup`. Activates `specific_lookup` + `multi_lookup` when papers introduce ID columns (placeholder for cross-paper internal-ID harvesting; mirrors gene-side `gene_id_mapping.json`).
- Paired-comparison fold-change/p-value edges deferred — no current paper has them. If they appear, design a separate node/edge type rather than extending `MetaboliteAssay`.

## Out of scope (explicit)

- `MassFeature` node type for unidentified peaks
- Isotopologue tracing
- `Assay_classifies_metabolite` (no current paper needs categorical metabolite values)
- Auto-pivot of long-format CSVs (paperconfig is explicit about column → assay mapping)
- `UnresolvedMetabolite` node type — unresolved rows produce no edge
- Lipid-class hierarchy (LipidMaps integration)
