# Gene-ID mapping hygiene bundle (backlog #13–#16) — 2026-08-28

**Scope.** Four `build_gene_id_mapping` / `gene_id_utils` items, no schema change:
1. **GFF `Name` typing (root cause of #14 KT2440 + #16 MED4).** `extract_rows_from_annotation_gff`
   typed every `Name=` attribute as Tier-1 `locus_tag_ncbi`. On NCBI GFFs `Name` is the gene
   symbol whenever one exists (`Name=tnpB;gene=tnpB` — 1,338 of 5,587 KT2440 genes, 525 of
   1,920 MED4) and the protein accession on CDS rows. A shared symbol is gene-unique by
   declaration, so the 12 `tnpB` IS copies merged into `M8001_03425` (40 Tier-1 ids) and the
   5S rRNA `rrf` (RNA_41) merged into `frr` (PMM0521, whose synonym is also `rrf`).
   Fix: `Name == gene` → `gene_name` (Tier 3); `Name == protein_id` → `protein_id_refseq`
   (Tier 2); `Name == locus_tag/old_locus_tag` → skip (duplicate); anything else stays
   `locus_tag_ncbi`.
2. **Runaway guard (#13).** `GeneIdGraph.build_diagnostic_report` adds `tier1_count_median`,
   `tier1_count_max`, `runaway_genes` (Tier-1 count > 3 × median, ≥ 12) and a `[RUNAWAY]`
   warning. Report-only — the builder never fails a strain — but the warning is what makes
   the next merge visible instead of a manual re-check.
3. **`gene-` prefix heuristic (#15).** `_heuristic_candidates` also offers the value with a
   leading `gene-` / `cds-` / `rna-` GFF ID prefix stripped, so `gene-PMM0950` resolves through
   `PMM0950` when the prefixed form itself is absent or conflicted.
4. Re-run `prepare_data.sh --steps 3 4 8 --force`; compare per-strain Tier-1 max/median,
   conflicts and every `*_resolved_report.txt` against the baseline captured before the change.

**Not changing.** Paperconfigs; the GFF `gene` attribute handling; Tier-3 semantics; the
`id_translation` diamond generator.

**Acceptance.** Unit tests for the three changes; `pytest -m "not slow and not kg"` green;
after the rerun: KT2440 max Tier-1 ≈ median, `RNA_41` no longer maps to PMM0521, he 2022
MED4 unresolved < 94 with no new unresolved elsewhere, no strain loses resolved rows except
where a shared gene symbol was the only (wrong) bridge.

**Status — DONE 2026-08-28.** Also added on the way: Pass 2b (`heuristic_multi`, recovers kaur 2018's
unversioned accessions, −92 → 0) and Pass 3a (`multi_named`, recovers `atpB`/`atpC`-style rows).
Measured vs the pre-change baseline: supp tables 210,247 → 210,162 resolved rows (−85, all
multi-copy symbols; he 2022 +6), narrative mentions 1,373 → 1,315 (−58: multi-copy symbols such as
`hli16`/`psbA`/`pstS`, plus fadeev 2022's product strings that had case-insensitively matched a
Tier-1 token — junk resolutions, correctly gone); Tier-1 conflicts 2,300 → 114; max Tier-1 ids/gene
40 → 14; 0 runaways. `RNA_41` no longer resolves (rRNA gene, no Gene node — correct).
`pytest -m "not slow and not kg"` 2481 passed. Needs a Docker rebuild to reach the graph.
