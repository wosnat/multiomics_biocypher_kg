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

## Status

- [ ] Run Cypher investigation queries (step 1–4 above)
- [ ] Determine whether gap is regression or pre-existing
- [ ] Choose fix (A, B, or C) and implement
- [ ] Update or remove failing tests once fixed
