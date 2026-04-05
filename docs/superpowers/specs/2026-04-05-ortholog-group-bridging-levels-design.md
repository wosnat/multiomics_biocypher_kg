# Ortholog Group Bridging Levels

## Problem

The `ORGANISM_GROUP_LEVELS` whitelist in `ortholog_group_utils.py` only covers 3 organism groups (Prochlorococcus, Synechococcus, Alteromonas). Five recently added diverse organisms (Shewanella, Pseudomonas, Ruegeria, Meiothermus, Thermosynechococcus) only get Bacteria-level OGs (taxon_id=2). This means cross-organism comparisons between heterotrophs, and between cyanobacteria across genera, can only happen through very broad Bacteria-level groups (~8.7K nodes).

## Goal

Enable finer-grained cross-organism gene comparison for:
- **Heterotroph-to-heterotroph**: Alteromonas, Shewanella, Pseudomonas, Ruegeria
- **Cyanobacteria across genera**: Prochlorococcus, Synechococcus, Thermosynechococcus

Not in scope: cyanobacteria-to-heterotroph bridging (Bacteria-level COGs are sufficient), organism-specific family-level OGs for diverse organisms (no value for cross-comparison).

## Design

### Config change: `ORGANISM_GROUP_LEVELS`

Replace the `(target_tid, fallback_tid)` tuple with a list of `(taxon_id, specificity_rank)` pairs. Each level is extracted independently (no fallback logic).

```python
ORGANISM_GROUP_LEVELS: dict[str, list[tuple[int, int]]] = {
    "Prochlorococcus":     [(1212, 1), (1117, 2)],    # Prochloraceae, Cyanobacteria
    "Synechococcus":       [(1129, 1), (1117, 2)],    # Synechococcus, Cyanobacteria
    "Thermosynechococcus": [(1117, 2)],                # Cyanobacteria
    "Alteromonas":         [(72275, 1), (1224, 2)],    # Alteromonadaceae, Proteobacteria
    "Shewanella":          [(1224, 2)],                 # Proteobacteria
    "Pseudomonas":         [(1224, 2)],                 # Proteobacteria
    "Ruegeria":            [(1224, 2)],                 # Proteobacteria
    "Meiothermus":         [],                          # Bacteria-level only
}
```

### Bridging levels

| Level | Taxon ID | Rank | Organisms bridged |
|---|---|---|---|
| Cyanobacteria | 1117 | 2 | Prochlorococcus + Synechococcus + Thermosynechococcus |
| Proteobacteria | 1224 | 2 | Alteromonas + Shewanella + Pseudomonas + Ruegeria |

Meiothermus (Deinococcus-Thermus) has no shared intermediate level with any other organism group -- Bacteria-level OGs are sufficient.

### Behavior change

**Before:** Each gene gets at most 1 intermediate-level OG (target OR fallback).
- Pro/Syn genes with Prochloraceae OG do NOT get Cyanobacteria OG.
- Alteromonas genes with Alteromonadaceae OG do NOT get Gammaproteobacteria/Proteobacteria OG.

**After:** Each gene gets ALL configured levels independently.
- A Prochlorococcus gene gets Prochloraceae AND Cyanobacteria OGs.
- An Alteromonas gene gets Alteromonadaceae AND Proteobacteria OGs.
- A Shewanella gene gets Proteobacteria OG (plus Bacteria-level, as before).

Max OGs per gene: Cyanorak (0-1) + family (0-1) + intermediate (0-1) + Bacteria (0-1) = 2-4.

### Files changed

1. **`multiomics_kg/download/utils/ortholog_group_utils.py`**
   - Change `ORGANISM_GROUP_LEVELS` type from `dict[str, tuple[int, int]]` to `dict[str, list[tuple[int, int]]]`
   - Rewrite `extract_ortholog_groups()`: replace either/or fallback with loop over all configured levels
   - No change to `organism_group_from_path()` (path matching works for new groups since dir names match)

2. **No changes to `ortholog_group_adapter.py`** -- it reads whatever `ortholog_groups` list is in `gene_annotations_merged.json` and deduplicates by `og_id`.

3. **No schema changes** -- no new node or edge types.

### Rebuild steps

1. `bash scripts/prepare_data.sh --steps 2` -- rebuild `gene_annotations_merged.json` for all strains
2. Rebuild KG (`uv run python create_knowledge_graph.py`)
3. Verify OG counts increased as expected

### Expected impact

- ~2,500-4,000 new Proteobacteria-level OG nodes
- ~200-500 new Cyanobacteria-level OG nodes (most already exist from Pro/Syn fallback)
- ~15,000-25,000 new `Gene_in_ortholog_group` edges
- Pro/Syn genes gain explicit Cyanobacteria OG membership (previously only fallback)
- Alteromonas genes gain Proteobacteria OG membership (previously only Alteromonadaceae)

### Testing

1. **Unit tests** (`tests/test_ortholog_group_extraction.py`):
   - Test new `ORGANISM_GROUP_LEVELS` format: each organism group extracts all configured levels (not either/or)
   - Test that a Prochlorococcus gene gets BOTH Prochloraceae AND Cyanobacteria OGs
   - Test that an Alteromonas gene gets BOTH Alteromonadaceae AND Proteobacteria OGs
   - Test that a Shewanella gene gets Proteobacteria + Bacteria (no family level)
   - Test that a Meiothermus gene gets only Bacteria-level OG
   - Test deduplication: same `og_id` not added twice

2. **Integration verification** (after KG rebuild):
   - Query OG node counts by source/level — expect new Proteobacteria and Cyanobacteria rows
   - Verify Thermosynechococcus genes have Cyanobacteria-level OG membership
   - Verify Shewanella/Pseudomonas/Ruegeria genes have Proteobacteria-level OG membership
   - Verify cross-organism bridging: find shared Proteobacteria-level OGs between Alteromonas and Shewanella
   - Verify cross-organism bridging: find shared Cyanobacteria-level OGs between Prochlorococcus and Thermosynechococcus
   - Existing OG counts for Cyanorak and Bacteria levels should be unchanged
   - Run existing `pytest -m "not slow and not kg"` — all must pass

### Code review checklist

- `ORGANISM_GROUP_LEVELS` type annotation matches new structure
- `extract_ortholog_groups()` loop handles empty level list (Meiothermus)
- `organism_group_from_path()` matches all 8 organism groups in the dict
- No leftover fallback logic (target_tid/fallback_tid variables removed)
- Existing tests updated to match new config format
- No changes to `ortholog_group_adapter.py` (confirm it's truly untouched)
- String literals consistent: taxonomic level names come from eggnog data, not hardcoded

### Not changing

- Bacteria-level extraction (always applied regardless of whitelist)
- Cyanorak extraction (Pro/Syn only, from cluster_number)
- OrthologGroup node properties, consensus computation, majority-vote edges
- `ortholog_group_adapter.py` (reads from JSON, no logic changes)
- Schema config
- Gammaproteobacteria level (dropped -- Proteobacteria alone bridges all 4 heterotrophs)
- Family-level OGs for diverse organisms (no cross-comparison value)
