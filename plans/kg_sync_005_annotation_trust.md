# KG-SYNC-005 — annotation-trust surface (pre-alpha.7)

**Driver:** `multiomics_explorer/docs/kg-specs/2026-08-27-annotation-trust-kg-asks.md` (asks ONT-001…015,
KG verdicts §6, explorer acceptance §7, KG reply §8). All asks accepted; ships as one batch before the
0.1.0-alpha.7 cut. Status (2026-08-27): **all 6 steps done, all 7 acceptance criteria verified live**; uncommitted on `main` pending user review.

## Scope

### Changing
| Area | Change | Ask |
|---|---|---|
| adapters: kegg_annotation, cog/role (category_role), interpro, ncbifam, merops | constant `sources: str[]` on `Gene_has_kegg_ko` (`['eggnog']`), `Gene_in_cog_category` (`['eggnog']`), `Gene_has_cyanorak_role` + `Gene_has_tigr_role` (`['cyanorak']`), `Gene_has_interpro_entry` + `Gene_has_ncbifam_family` (`['interproscan']`), `Gene_has_merops_family` (`['merops_diamond']`) | ONT-007 |
| same adapters + tcdb_adapter | `evidence` on all 14 gene→ontology edge types. Ladder `curated > signature > homology > family_inferred > domain_inferred`. KO/COG/eggNOG-only TCDB → `family_inferred`; cyanorak/tigr roles → `curated`; interpro/ncbifam → `signature`; merops + diamond-sourced TCDB → `homology` (both-source TCDB → `homology`). Per-edge vocab entries keep their own subset. | ONT-008 |
| ncbifam_adapter | `score` → `bit_score`; retired accessions get `family_type = 'retired'` | ONT-011, ONT-004 |
| merops_adapter + schema + vocab | `MeropsFamily.family_type` → `family_class` | ONT-012 |
| post-import (.sh + .cypher) | `Gene_has_tcdb_family.attachment_depth: most_specific \| superseded`, materialized once; the 3 existing deepest-attachment sites read the property | ONT-010 |
| post-import | `Gene_has_merops_family.evidence_score` (2 signals: `tier <= 2`, `pfam_support = 'corroborated'`; `signal_count`, `signals`), `Gene.merops_evidence_score_max` (sparse) | ONT-006, ONT-013 |
| post-import | `MeropsFamily.peptidase_organism_count` | ONT-005 |
| post-import | `InterproEntry.gene_count`/`organism_count` → subtree; `direct_gene_count` on InterproEntry, TcdbFamily, CazyFamily, MeropsFamily | ONT-009 |
| post-import | subtree `gene_count`/`organism_count` + `direct_gene_count` on BiologicalProcess, MolecularFunction, CellularComponent (`is_a` ∪ `part_of`), EcNumber, KeggTerm (all levels), PfamClan, CyanorakRole, BriteCategory (direct only added); direct `gene_count`/`organism_count` on Pfam, TigrRole, CogFunctionalCategory | ONT-015 |
| controlled_vocabulary_adapter | vocab-derived "null-producing libraries" text for `Gene_has_interpro_entry.evalue` computed from calls.json at build (controlled_vocabulary_adapter or a builder hook) | ONT-002 |
| config/controlled_vocabularies.yaml | ~25 entries: `sources` ×7, `evidence` ×8 (incl. new TCDB), `Gene_has_merops_family.evidence_score`, `Gene.merops_evidence_score_max`, `Gene_has_tcdb_family.attachment_depth`, numeric entries (`Gene_has_interpro_entry.evalue`/`match_count`, `Gene_has_ncbifam_family.evalue`/`bit_score`/`start`/`end`, diamond props on tcdb + merops), `Gene_has_subcellular_localization.score`, `Gene_has_signal_peptide_type.probability`/`cleavage_site`/`cleavage_probability`, `NcbifamFamily.family_type` += `retired`, rename `MeropsFamily.family_type` → `family_class` | ONT-002/003/014 |
| docs | CLAUDE.md, `docs/kg-changes/interpro-multi-ontology.md` (strike `match_count` on ncbifam edge), `merops-extension.md`, `vocabulary-contract.md` (ladder + `homology`), CHANGELOG `### Breaking` (InterPro `gene_count` semantics; `score`→`bit_score`; `family_type`→`family_class`) + `### Added` | ONT-001 + all |

### Not changing
- No new node types, no new edge types. `Gene_has_interpro_entry` gets no `score` (count-don't-combine stands).
- `annotation_quality` buckets (9) and `informative_annotation_types` gating unchanged.
- `Gene_has_merops_family.confidence_score` keeps its name (KG-SYNC-003 stands).
- No `evidence` split of `homology` by tier.
- MCP/explorer code — separate repo, separate slice.

### Acceptance criteria
1. Every one of the 14 gene→ontology edge types carries `sources` and `evidence` on 100% of edges (live `count(r) = count(r.sources) = count(r.evidence)`).
2. `Gene_has_merops_family.evidence_score ∈ {0, 0.5, 1}` on all 4,257 edges; `round(score*2)` = fired signals.
3. MED4: 670 TCDB rows, 73 `attachment_depth = 'superseded'`; corpus `Gene.transported_metabolite_count` and `Metabolite.transporter_gene_count` byte-identical to pre-batch (`post-import-validate.sh` diff) — proves the materialized predicate equals the inline one.
4. `direct_gene_count` present on every hierarchical ontology label; `gene_count >= direct_gene_count` everywhere; InterproEntry `gene_count` change vs baseline reported (expected small).
5. `NcbifamFamily`: 0 nodes with null `family_type`; `Gene_has_ncbifam_family.bit_score` on 67,459 edges, `score` absent.
6. Vocab-drift test green; explorer §7.3 `(applies_to, property)` list fully covered by `ControlledVocabulary` nodes.
7. `pytest -m "not slow and not kg"` and `pytest -m kg` green; snapshot regenerated; `capture_annotation_state --compare` shows no bucket movement (this batch adds no bucket).

## Steps
1. Baselines: `post-import-validate.sh > baseline.txt`, `capture_annotation_state.py --save`, `/omics-edge-snapshot`.
2. Vocab yaml + schema_config (renames, new props).
3. Adapters (sources/evidence/bit_score/retired/family_class) + unit tests.
4. Post-import `.cypher` + `.sh` in lockstep (attachment_depth first, then counts, then merops score) + `tests/kg_validity` assertions for AC 1–5.
5. Docker rebuild; run AC checks; diff baselines.
6. Docs + CHANGELOG; notify explorer (§7.4: schema-baseline refresh).

## Progress log
- 2026-08-27: baselines captured (`scratch_baselines/`, gitignored scratch). Adapters, schema, vocab
  (110 entries), post-import (.cypher + .sh in lockstep), unit tests (2,413 pass), KG-validity tests
  (`tests/kg_validity/test_annotation_trust.py`), docs (CLAUDE.md, CHANGELOG, `annotation-trust-surface.md`,
  merops-extension / interpro-multi-ontology / vocabulary-contract addenda) all written. Deviations from §8
  of the asks doc: (a) `direct_gene_count` NOT emitted on PfamClan / BriteCategory (constant 0 — vacuous);
  (b) the InterPro null-e-value library list is measured by script and written into the vocab description,
  not recomputed at every build. Docker rebuild launched.
- 2026-08-27 (later): Docker rebuild OK (build/import/post-process exit 0; post-process group 3 55.7 s).
  AC1–AC7 verified: `pytest -m kg` 1,187 pass / 4 skip; `post-import-validate.sh` byte-identical to
  baseline; `capture_annotation_state --compare` unchanged everywhere; snapshot regenerated. Live: 110
  vocab nodes; merops score 151/3,768/338; attachment_depth 46,593/7,170; InterPro subtree≠direct on 115
  entries. §9 landed-note appended to the explorer asks doc.
