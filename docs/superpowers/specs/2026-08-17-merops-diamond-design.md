# MEROPS-Diamond Phase 1 — Design

**Date:** 2026-08-17
**Status:** Phase 1 (per-strain artifacts only; KG integration deferred)
**Pattern:** tcdb-diamond clone (`docs/superpowers/specs/2026-05-10-tcdb-diamond-augmentation-design.md`)

## Step 0 — Intent

- **Predicts:** per-protein list-of-candidates — MEROPS peptidase/inhibitor
  classification (family / subfamily / identifier) via diamond blastp against
  the curated `merops_scan.lib`, with tiered confidence truncation following
  the tcdb-diamond policy. Extra per-candidate flags: non-peptidase homolog,
  putative entry, inhibitor family.
- **Input:** `protein.faa`
- **Install flavor:** B — `diamond` on PATH + downloadable reference
  (`merops_scan.lib` 1.9 MB + `family.txt`/`clan.txt`) → `~/tools/MEROPS/`
- **Triggers:** "run merops", "classify proteases / peptidases for strain X",
  "protease families", "peptidase annotation", "exoprotease classification",
  "secreted protease families"

## Why

The KG identifies proteases only through text matching of product/domain
annotations. MEROPS is the authoritative protease classification (clan →
family → subfamily → identifier); systematic assignment was Tier 1 in
`docs/paper/kg_improvements_exoenzyme_analysis.md` §8 and is the last Tier-1
item not yet integrated (SignalP + PSORTb are done). Biological payoff:
extracellular peptidase classification for the heterotroph DOM-degradation /
coculture story, joined with the existing SignalP + PSORTb layers.

## Reference data (verified 2026-08-17 against MEROPS release 12.5)

From `https://ftp.ebi.ac.uk/pub/databases/merops/current_release/`:

| File | Size | Use |
|---|---|---|
| `merops_scan.lib` | 1.9 MB, 5,009 seqs | The scan library — one representative peptidase-unit sequence per MEROPS identifier. Diamond DB source. |
| `database_files/family.txt` | small | family → clan, type-example name, peptidase/inhibitor flag |
| `database_files/clan.txt` | small | clan → description |

NOT used in Phase 1 but relevant for Phase 2 (all verified present):
`interpro.txt` (182 MB, per-UniProt-accession MEROPS↔Pfam↔InterPro map — the
structured corroboration bridge InterPro's own XML lacks: `db="MEROPS"` xrefs
= 0 in the 2026-08 InterPro release), `database_files/GO_annotation.txt`
(MEROPS GO map → potential `Merops_family_*` GO bridges, TCDB-parallel),
`Substrate_search.txt` (cleavage specificity → node properties).
`pepunit.lib` (457 MB) and `protease.lib` (867 MB) are deliberately not used —
the scan library is MEROPS's curated non-redundant search set.

### Scan-library header format (all 5,009 parse with one regex)

```
>MER0000002 - chymotrypsin A (cattle-type) (Bos taurus) [S01.001]#S01A#{peptidase unit: 34-263}~source XP_608091~
```

`MERNUM - name (organism) [merops_id]#subfamily_or_family#{unit}~source ...~`

- The `#...#` token is the subfamily (`S01A`) when the family has subfamilies,
  else the family code itself (2,212 of 5,009 headers).
- File is latin-1 encoded (0x96 en-dashes in names) — not UTF-8.

### MEROPS identifier tail shapes (observed distribution)

| Tail shape | Count | Meaning | `best_hit_kind` |
|---|---:|---|---|
| numeric 001–899 (`S01.001`) | 3,091 | characterized holotype | `holotype` |
| letter+num (`S01.A33`) | 1,463 | putative / unassigned peptidase | `putative` |
| `9xx` (`S01.971`) | 455 | **non-peptidase homolog** (catalytically dead) | `nonpeptidase_homolog` |

Catalytic types by family first letter: S 1,708 · M 1,089 · C 1,014 ·
**I 674 (inhibitors)** · A 319 · T 105 · P 53 (mixed) · N 24 · U 15 · G 8.
Inhibitor families are kept (protease inhibitors are biology) and flagged
`entry_type: "inhibitor"` with `catalytic_type: null`.

## Tier policy (tcdb-diamond parity)

Floor (identical): e-value ≤ 0.001, HSP length ≥ 50, qcov ≥ 40 OR scov ≥ 40.

| Tier | Identity | Coverage | Claim (level_kind) |
|---|---|---|---|
| 1 | ≥ 70% | qcov ≥ 70% | MEROPS identifier (`merops_id`) |
| 2 | ≥ 40% | qcov ≥ 60% | subfamily (`merops_subfamily`) |
| 3 | floor | floor | family (`merops_family`) |

Consensus collapse per **family** group (analog of TCDB's 3-part grouping):
agreement depth = `id` (one distinct identifier) > `subfamily` (one distinct
real subfamily) > `family`; weights 1.0 / 0.85 / 0.7.
`effective_tier = max(best_hit_tier, depth_tier)`; code truncated to the
depth the effective tier justifies.

**Ragged-hierarchy deviation from TCDB:** MEROPS families without subfamilies
(e.g., S33) make tier-2 claims land at family level — `tier: 2` with
`level_kind: "merops_family"`. TCDB could lock tier↔depth because every TCID
has 5 parts; MEROPS cannot. `tier` stays the confidence band; `level_kind`
describes the actual claim depth.

**Non-peptidase homologs:** a diamond hit cannot distinguish an active
peptidase from a dead homolog, but the subject's identifier tail can. Each
candidate carries `best_hit_kind` and `homolog_hit_fraction` (fraction of the
family group's hits that are `9xx` entries). No filtering — evidence rides on
the artifact, consumers threshold explicitly (tcdb-diamond principle).

## Pure sequence evidence

Same contract as tcdb-diamond post-2026-08-06: reads ONLY `protein.faa` + the
MEROPS diamond DB + the static family→clan reference. No eggNOG, no
`gene_annotations_merged.json` (avoids the step-2 cycle). There is no second
sequence source for MEROPS (eggNOG doesn't emit it; InterProScan has no MEROPS
member DB), so Phase 2 starts single-source; corroboration options are the
MEROPS-published `interpro.txt` Pfam bridge and/or a structure arm (Foldseek).

## Output

`cache/data/<org>/genomes/<strain>/merops/`:
- `<strain>.merops.tsv` — raw diamond 8-column output
- `<strain>.merops.calls.json` — keyed by WP_ accession, `{"calls": [candidate, ...]}` sorted by confidence_score desc
- `<strain>.merops.skill_summary.json` — per-strain QC

Candidate fields: `code` (called, truncated), `family`, `subfamily` (null when
not claimed / none exist), `clan` (null when family unassigned, clan codes
ending `-` are the unassigned sentinel), `catalytic_type` (single letter, null
for inhibitors), `entry_type` (`peptidase`|`inhibitor`), `level_kind`, `tier`,
`confidence_score`, `identity`, `qcov`, `scov`, `evalue`, `length`,
`consensus_n`, `consensus_agreement` (`id`|`subfamily`|`family`),
`best_hit_id`, `best_hit_mernum`, `best_hit_kind`, `homolog_hit_fraction`.

## Phase 2 sketch (deferred — separate spec via /integrate-a-tool)

Ontology track, CAZy-shaped: `MeropsFamily` nodes (clan → family → subfamily,
observed-only pruning; catalytic type as node **property**, not a level — clans
mix catalytic types, e.g. clan PA holds S1 + C3) + scored
`Gene_has_merops_family` edges. Merge front door: new `merops_diamond` logical
source. Bridges: `interpro.txt`-derived Pfam/InterPro links,
`GO_annotation.txt` GO links (both TCDB-bridge-shaped). No substrate→Metabolite
arm (peptidase substrates are proteins) — cleavage specificity becomes node
properties. Keep out of `annotation_quality` initially (single-source; tier-3
TCDB precedent).
