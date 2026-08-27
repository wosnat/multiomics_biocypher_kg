# Orphan Proteins Investigation Plan

## Problem

~46% of UniProt `Protein` nodes in the graph have no `Protein_belongs_to_organism` edge and no
`Gene_encodes_protein` edge (12,985 out of 28,392 proteins as of Feb 2026).

Two KG validity tests are failing because of this:
- `test_no_orphan_proteins` — requires 0 orphans, currently 12,985
- `test_no_orphan_proteins_without_gene` — threshold 15%, currently 45.7%

## Root Cause

Both edges are created in `multiomics_kg/adapters/uniprot_adapter.py` via a RefSeq WP_ join:

```python
# uniprot_adapter.py lines 205–218
for refseq in refseq_ids:
    for locus_tag, ncbi_acc in self._refseq_to_strains.get(refseq, []):
        yield None, protein_id, gene_id, "Gene_encodes_protein", props
        yield None, protein_id, org_id, "Protein_belongs_to_organism", props
```

A protein only gets organism/gene edges if:
1. UniProt reports at least one RefSeq WP_ ID for it, AND
2. That WP_ ID is present in `gene_mapping.csv` for one of our strains

Proteins that fail either condition are stranded nodes with no graph connections other
than GO/EC/pathway edges.

## Open Question: Regression or Pre-existing Gap?

The Feb 2026 UniProt adapter refactor (`fe5c2bb`, "refactor uniprot adapter and gene
mapping skills") made large changes to how proteins are fetched and cached. It is unclear
whether the 46% orphan rate existed before this refactor or was introduced by it.

## Investigation Steps

1. **Check which proteins are orphaned**
   ```cypher
   MATCH (p:Protein)
   WHERE NOT (p)-[:Protein_belongs_to_organism]->()
   RETURN p.id, p.protein_synonyms, p.refseq_ids
   LIMIT 20
   ```
   - Do the orphaned proteins have `refseq_ids` set? If yes → WP_ IDs are present in
     UniProt but not in gene_mapping.csv (NCBI mismatch).
   - If no refseq_ids → UniProt doesn't report WP_ cross-references for these proteins.

2. **Quantify: WP_ present but not in gene_mapping vs. no WP_ at all**
   ```cypher
   MATCH (p:Protein)
   WHERE NOT (p)-[:Protein_belongs_to_organism]->()
   RETURN
     count(CASE WHEN p.refseq_ids IS NOT NULL THEN 1 END) AS has_refseq,
     count(CASE WHEN p.refseq_ids IS NULL THEN 1 END) AS no_refseq
   ```

3. **Check gene_mapping.csv coverage** — compare WP_ IDs in UniProt protein cache vs.
   WP_ IDs in `gene_mapping.csv` files. If UniProt returns far more WP_ IDs than are
   in our gene maps, the NCBI data (gene_mapping.csv) is the limiting factor.

4. **Check if proteins are taxid-specific** — UniProt API is queried per taxid. Are
   orphans concentrated in specific strains? Could indicate a stale cache or failed
   download for some strains.

5. **Git bisect / diff** — compare `_load_gene_mapping()` and protein download logic
   before and after `fe5c2bb` to see if the WP_ join logic changed.

## Possible Fixes

### Option A: Accept the gap, relax test thresholds
If investigation shows this is a pre-existing, expected limitation (UniProt returns many
proteins without RefSeq WP_ IDs for cyanobacteria), update the tests:
- `test_no_orphan_proteins`: change `== 0` to `< 0.50` with a clear comment
- `test_no_orphan_proteins_without_gene`: raise threshold from 0.15 to ~0.55

### Option B: Add taxid-based fallback for Protein_belongs_to_organism
Even if a protein has no WP_ match, we know which taxid it was downloaded for.
Emit `Protein_belongs_to_organism` using the taxid → organism mapping directly,
without requiring the RefSeq join. This would fix orphan organism edges for all proteins
while still only creating `Gene_encodes_protein` where WP_ matches.

### Option C: Fix the RefSeq join upstream
Improve `gene_mapping.csv` to include more WP_ IDs (e.g. by pulling them from the
GenBank/GFF files more aggressively), or improve the UniProt download to request
RefSeq cross-references more completely.

## Investigation findings (2026-08-27, investigate-only — nothing changed)

Measured on the live graph (67,024 Protein nodes) and reproduced **exactly** from the
per-taxid `protein_annotations.json` + `gene_mapping.csv` caches, so the numbers below
are reproducible without Neo4j.

**25,441 orphans (38.0%)** = every protein whose UniProt `xref_refseq` does not hit a
`protein_id` in one of the taxid's `gene_mapping.csv` files. Both edges share that one
join, so "no organism edge" and "no gene edge" are the *same* 25,441 proteins.

### 1. Regression vs pre-existing — answered

- `Protein_belongs_to_organism`: **regression from `fe5c2bb`.** The pre-refactor
  `_create_edges_of_type` emitted the organism edge for *every* protein of the taxid
  (one edge per assembly, from `organism_id`, no RefSeq join). The refactor moved it
  inside the `refseq_to_strains` loop, so it now inherits the WP_ join's failure rate.
  Pre-refactor the organism-orphan count was 0 by construction.
- `Gene_encodes_protein`: **pre-existing.** Before and after the refactor it was the
  same RefSeq WP_ join (`_create_gene_to_protein_edges` → `get_edges`), so the ~38%
  gene-unlinked fraction predates `fe5c2bb`; the 15% threshold was never met.

### 2. Why the WP_ join fails — four populations

| Class | Proteins | What it is |
|---|---|---|
| `noWP_oln_rescuable` | **7,596** | UniProt entry has **no RefSeq xref at all** but its `gene_oln`/`gene_names` is a locus tag already in our `gene_mapping.csv` (e.g. MED4 `A8WI08` ↔ `PMM1805`). Uniform ~15% "tail" on every Pro/Syn strain; the largest single rescuable block. |
| `noWP_unrescuable` | 4,722 | No RefSeq xref and no locus-tag match. Concentrated in shared-taxid downloads (28108: 1,423; KT2440: 901; DSS-3: 623; PCC7942: 438). |
| `WP_foreign_proteome` | 7,496 | Has WP_, but the entry's `xref_proteomes` is a *different* UniProt proteome from the one our assembly matches — another isolate/assembly under the same taxid. 28108 (species-level *A. macleodii* taxid shared by MIT1002/EZ55/HOT1A3/BGP6) → `UP000509458` 3,433; DSS-3 → `UP000813672` 4,055. |
| `WP_not_in_our_assembly` | 5,627 | Has WP_, same proteome accession as our matched proteins, but the WP_ is **absent from our downloaded `protein.faa` too** (0/5,627 present, version-insensitive). `gene_mapping.csv` is NOT the limiter — its `protein_id` set equals the `protein.faa` id set for every strain checked. These are proteins the UniProt proteome carries that our RefSeq annotation build does not (different annotation release / GCA-vs-GCF build). 28108: 4,632; MruberA: 696; DSS-3: 171. |

Two structural drivers dominate the WP_ classes: (a) **species-level or shared
taxids** pull in every isolate UniProt knows (28108 alone is 9,637 of the 25,441
orphans = 80.6% orphan rate for that download), and (b) UniProt proteomes built on a
different annotation build than our `GCF_*` download.

Side finding: taxid **172827 resolves to *Meiothermus taiwanensis*** in UniProt
(`organism_name` on every entry), while the registry row is *Meiothermus ruber*
`GCF_000836395.1` — a registry/taxid mismatch worth a separate look.

### 3. Implications for the options above

- **Option B** (organism edge from taxid, no join) restores the pre-`fe5c2bb`
  behaviour and takes organism-orphans to 0 — but for shared taxids it would attach
  ~9.6K foreign-isolate proteins to all four *A. macleodii* strains. Needs a
  proteome-aware guard (`xref_proteomes` ∈ our assembly's proteome) to be correct,
  not just green.
- **Option C** is mostly *not* about `gene_mapping.csv` coverage; the only real
  upstream lever is the 7,596 `gene_oln` locus-tag join (secondary key when
  `xref_refseq` is empty), which would lift `Gene_encodes_protein` coverage from
  62.0% → ~73.3%.
- **Option A** thresholds: with B+OLN the honest targets are organism-orphans = 0
  and gene-unlinked ≈ 27% (foreign/other-build proteins are not linkable to genes we
  have).

## Implementation (2026-08-27) — code done, live verification pending

`multiomics_kg/adapters/uniprot_adapter.py`:

1. **Keep rule.** `_resolve(entry)` → (gene pairs, assemblies). A protein is emitted
   (node + edges) only when `assemblies` is non-empty: assemblies come from the
   RefSeq join, else the locus-tag join, plus any assembly whose *own* proteome the
   entry lists. Unlinkable entries are dropped instead of left as orphans.
2. **Own-proteome derivation** (`_derive_own_proteomes`): an assembly owns a UniProt
   proteome when its WP_-matched proteins cover ≥ `OWN_PROTEOME_MIN_FRACTION` (0.5)
   of that proteome's entries in the taxid download, ≥ `MIN_OWN_PROTEOME_SUPPORT`
   (20) of them, and no other assembly also covers a majority. Direction matters —
   measuring how many of *our* matched proteins list the proteome is fooled by
   conserved proteins (MIT1002 lists `UP000063991` on 59% of matches yet 80% of that
   proteome's WP_ are absent from its FASTA). Measured coverage: own proteomes
   71–100%; DSS-3's other build `UP000813672` 2%; every proteome under 28108 ≤ 27%
   → the 9,637 proteome-only *A. macleodii* proteins are unattributable and dropped.
3. **Organism edge** is emitted per assembly in `assemblies`, independent of the gene
   join — the pre-`fe5c2bb` semantics.
4. **Locus-tag fallback join** (`_locus_to_strains`, built by `_load_gene_mapping`
   from `locus_tag`/`locus_tag_ncbi`/`old_locus_tags`/`locus_tag_cyanorak`): used only
   when the RefSeq join has no hit, so a stale `gene_oln` cannot add a contradictory
   second gene edge.

Offline dry run over all 40 taxid caches (no Docker): **52,735 proteins kept
(−14,289), 51,927 gene edges, 55,667 organism edges, 100% with organism,
92.6% with gene (was 62.0%), 0 dangling.** Per-taxid breakdown is logged by
`download_data` (`… RefSeq-joined, … locus-tag-joined, … own-proteome only, … dropped`).

Tests: `tests/test_uniprot_adapter.py` (keep rule, proteome derivation incl. the
coverage-direction and shared-taxid cases, locus fallback, gene_mapping index) +
`tests/test_gene_protein_edges.py` updated for the tuple return. KG validity:
`test_no_orphan_proteins` → `== 0`, `test_no_orphan_proteins_without_gene` → `< 10%`.

Not done here: Docker rebuild, `pytest -m kg`, snapshot regeneration (the 10
sampled `uniprot:` ids in `snapshot_data.json` all survive the keep rule). The
Meiothermus 172827 naming drift is logged in `plans/backlog.md`, not changed.

## Status

- [x] Run Cypher investigation queries (step 1–4 above) — 2026-08-27
- [x] Determine whether gap is regression or pre-existing — organism edge: regression; gene edge: pre-existing
- [x] Choose fix (B + proteome gate + gene_oln join) and implement — 2026-08-27
- [x] Re-tighten kg tests (== 0 / < 10%) — verify on next rebuild
