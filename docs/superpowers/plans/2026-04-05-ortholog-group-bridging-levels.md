# Ortholog Group Bridging Levels Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add intermediate taxonomic bridging levels (Cyanobacteria, Proteobacteria) to ortholog group extraction so diverse organisms can be compared cross-genus.

**Architecture:** Change `ORGANISM_GROUP_LEVELS` from a target+fallback tuple to a list of levels per organism group. Each gene extracts ALL configured levels independently (no either/or). The downstream adapter reads whatever is in `gene_annotations_merged.json` unchanged.

**Tech Stack:** Python, pytest

**Spec:** `docs/superpowers/specs/2026-04-05-ortholog-group-bridging-levels-design.md`

---

### Task 1: Update tests for new `ORGANISM_GROUP_LEVELS` structure

**Files:**
- Modify: `tests/test_ortholog_group_extraction.py`

- [ ] **Step 1: Update `TestOrganismGroupLevels` to expect 8 groups with list format**

Replace the existing `TestOrganismGroupLevels` class at the bottom of the file:

```python
class TestOrganismGroupLevels:
    def test_all_groups_defined(self):
        expected = {
            "Prochlorococcus", "Synechococcus", "Thermosynechococcus",
            "Alteromonas", "Shewanella", "Pseudomonas", "Ruegeria",
            "Meiothermus",
        }
        assert set(ORGANISM_GROUP_LEVELS.keys()) == expected

    def test_values_are_lists_of_tuples(self):
        for group, levels in ORGANISM_GROUP_LEVELS.items():
            assert isinstance(levels, list), f"{group} should be a list"
            for tid, rank in levels:
                assert isinstance(tid, int) and tid > 0, f"{group} taxon_id {tid}"
                assert isinstance(rank, int) and rank > 0, f"{group} rank {rank}"

    def test_meiothermus_empty(self):
        assert ORGANISM_GROUP_LEVELS["Meiothermus"] == []
```

- [ ] **Step 2: Add `TestOrganismGroupFromPath` cases for new organisms**

Add these test methods to the existing `TestOrganismGroupFromPath` class:

```python
    def test_thermosynechococcus(self):
        assert organism_group_from_path("cache/data/Thermosynechococcus/genomes/BP1/") == "Thermosynechococcus"

    def test_shewanella(self):
        assert organism_group_from_path("cache/data/Shewanella/genomes/W3-18-1/") == "Shewanella"

    def test_pseudomonas(self):
        assert organism_group_from_path("cache/data/Pseudomonas/genomes/KT2440/") == "Pseudomonas"

    def test_ruegeria(self):
        assert organism_group_from_path("cache/data/Ruegeria/genomes/DSS-3/") == "Ruegeria"

    def test_meiothermus(self):
        assert organism_group_from_path("cache/data/Meiothermus/genomes/MruberA/") == "Meiothermus"
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `pytest tests/test_ortholog_group_extraction.py -v`
Expected: `TestOrganismGroupLevels` tests FAIL (wrong type, missing keys). New `TestOrganismGroupFromPath` tests FAIL (return "unknown").

- [ ] **Step 4: Commit failing tests**

```bash
git add tests/test_ortholog_group_extraction.py
git commit -m "test: add failing tests for expanded ORGANISM_GROUP_LEVELS and new organism paths"
```

---

### Task 2: Update `ORGANISM_GROUP_LEVELS` and `extract_ortholog_groups()`

**Files:**
- Modify: `multiomics_kg/download/utils/ortholog_group_utils.py`

- [ ] **Step 1: Replace `ORGANISM_GROUP_LEVELS` with new format**

Replace lines 12-19 (the dict and its type annotation):

```python
# Whitelist: organism group → list of (target_taxon_id, specificity_rank).
# Each gene extracts ALL levels in the list independently.
# Bacteria-level COGs (taxon_id=2, rank=3) are always added separately.
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

- [ ] **Step 2: Rewrite `extract_ortholog_groups()` — replace fallback logic with loop**

Replace the section after `# 2b. Lowest-level OG` (lines 105-127) with:

```python
    # 2b. Configured intermediate levels — extract ALL independently
    levels = ORGANISM_GROUP_LEVELS.get(organism_group, [])
    for tid, rank in levels:
        if tid in parsed:
            og_part, level_name = parsed[tid]
            og_id = f"eggnog:{og_part}@{tid}"
            if og_id not in seen_ids:
                groups.append({
                    "og_id": og_id,
                    "source": "eggnog",
                    "taxonomic_level": level_name,
                    "taxon_id": tid,
                    "specificity_rank": rank,
                })
                seen_ids.add(og_id)

    return groups
```

- [ ] **Step 3: Run tests to verify they pass**

Run: `pytest tests/test_ortholog_group_extraction.py -v`
Expected: All `TestOrganismGroupLevels` and `TestOrganismGroupFromPath` tests PASS.

- [ ] **Step 4: Commit**

```bash
git add multiomics_kg/download/utils/ortholog_group_utils.py
git commit -m "feat: expand ORGANISM_GROUP_LEVELS to 8 organism groups with multi-level extraction"
```

---

### Task 3: Update existing tests for new multi-level behavior

**Files:**
- Modify: `tests/test_ortholog_group_extraction.py`

The behavior change is: Pro/Syn/Alt genes now get BOTH family AND intermediate levels, not either/or. Several existing tests need updating.

- [ ] **Step 1: Update `test_pro_gene_three_groups` → now expects 4 groups**

A Prochlorococcus gene with Prochloraceae AND Cyanobacteria in eggnog_ogs now gets both. Replace the test:

```python
    def test_pro_gene_all_levels(self):
        """Prochlorococcus gene with cluster + eggnog → 4 groups (cyanorak + bacteria + prochloraceae + cyanobacteria)."""
        gene = {
            "cluster_number": "CK_00000001",
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
                "1MKTR@1212|Prochloraceae",
                "3Q4X@1117|Cyanobacteria",
            ],
        }
        groups = extract_ortholog_groups(gene, "Prochlorococcus")
        assert len(groups) == 4

        # Cyanorak
        assert groups[0] == {
            "og_id": "cyanorak:CK_00000001",
            "source": "cyanorak",
            "taxonomic_level": "curated",
            "taxon_id": 0,
            "specificity_rank": 0,
        }
        # Bacteria COG
        assert groups[1] == {
            "og_id": "eggnog:COG0592@2",
            "source": "eggnog",
            "taxonomic_level": "Bacteria",
            "taxon_id": 2,
            "specificity_rank": 3,
        }
        # Prochloraceae (rank 1)
        assert groups[2] == {
            "og_id": "eggnog:1MKTR@1212",
            "source": "eggnog",
            "taxonomic_level": "Prochloraceae",
            "taxon_id": 1212,
            "specificity_rank": 1,
        }
        # Cyanobacteria (rank 2)
        assert groups[3] == {
            "og_id": "eggnog:3Q4X@1117",
            "source": "eggnog",
            "taxonomic_level": "Cyanobacteria",
            "taxon_id": 1117,
            "specificity_rank": 2,
        }
```

- [ ] **Step 2: Update `test_lowest_level_uses_whitelist_not_max_taxon_id`**

This test still validates cross-lineage OGs are excluded. But now the gene gets Prochloraceae AND Cyanobacteria. Replace:

```python
    def test_cross_lineage_og_excluded(self):
        """Cross-lineage OG (Pleurocapsales) excluded; both Prochloraceae and Cyanobacteria included."""
        gene = {
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
                "AAAA@52604|Pleurocapsales",     # cross-lineage, should NOT be picked
                "1MKTR@1212|Prochloraceae",       # rank 1
                "3Q4X@1117|Cyanobacteria",        # rank 2
            ],
        }
        groups = extract_ortholog_groups(gene, "Prochlorococcus")
        og_ids = [g["og_id"] for g in groups]
        assert "eggnog:1MKTR@1212" in og_ids
        assert "eggnog:3Q4X@1117" in og_ids
        assert "eggnog:AAAA@52604" not in og_ids
```

- [ ] **Step 3: Update `test_lowest_level_falls_back_to_cyanobacteria`**

The concept of "fallback" no longer applies. This test now just shows that a gene missing Prochloraceae still gets Cyanobacteria. Rename and simplify:

```python
    def test_missing_family_still_gets_intermediate(self):
        """Gene missing Prochloraceae@1212 still gets Cyanobacteria@1117."""
        gene = {
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
                "3Q4X@1117|Cyanobacteria",
                # No Prochloraceae@1212 entry
            ],
        }
        groups = extract_ortholog_groups(gene, "Prochlorococcus")
        assert len(groups) == 2
        assert groups[1]["og_id"] == "eggnog:3Q4X@1117"
        assert groups[1]["taxonomic_level"] == "Cyanobacteria"
        assert groups[1]["specificity_rank"] == 2
```

- [ ] **Step 4: Update `test_alt_gene_two_groups` → now expects 3 groups with Proteobacteria**

```python
    def test_alt_gene_with_proteobacteria(self):
        """Alteromonas gene gets bacteria + Alteromonadaceae + Proteobacteria."""
        gene = {
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
                "4648R@72275|Alteromonadaceae",
                "5ABCD@1224|Proteobacteria",
            ],
        }
        groups = extract_ortholog_groups(gene, "Alteromonas")
        assert len(groups) == 3
        assert groups[0]["og_id"] == "eggnog:COG0592@2"
        assert groups[1]["og_id"] == "eggnog:4648R@72275"
        assert groups[1]["taxonomic_level"] == "Alteromonadaceae"
        assert groups[2]["og_id"] == "eggnog:5ABCD@1224"
        assert groups[2]["taxonomic_level"] == "Proteobacteria"
```

- [ ] **Step 5: Run tests**

Run: `pytest tests/test_ortholog_group_extraction.py -v`
Expected: ALL pass.

- [ ] **Step 6: Commit**

```bash
git add tests/test_ortholog_group_extraction.py
git commit -m "test: update ortholog group tests for multi-level extraction behavior"
```

---

### Task 4: Add tests for new organism groups

**Files:**
- Modify: `tests/test_ortholog_group_extraction.py`

- [ ] **Step 1: Add tests for diverse organisms in `TestExtractOrthologGroups`**

Add these test methods:

```python
    def test_shewanella_gets_proteobacteria(self):
        """Shewanella gene gets bacteria + Proteobacteria (no family level)."""
        gene = {
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
                "2QCTB@267890|Shewanellaceae",      # NOT configured — skipped
                "1NI9C@1224|Proteobacteria",
            ],
        }
        groups = extract_ortholog_groups(gene, "Shewanella")
        assert len(groups) == 2
        og_ids = [g["og_id"] for g in groups]
        assert "eggnog:COG0592@2" in og_ids
        assert "eggnog:1NI9C@1224" in og_ids
        # Shewanellaceae should NOT be extracted
        assert "eggnog:2QCTB@267890" not in og_ids

    def test_ruegeria_gets_proteobacteria(self):
        """Ruegeria gene gets bacteria + Proteobacteria."""
        gene = {
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
                "XXXX@97050|Ruegeria",               # NOT configured — skipped
                "YYYY@1224|Proteobacteria",
            ],
        }
        groups = extract_ortholog_groups(gene, "Ruegeria")
        assert len(groups) == 2
        assert groups[1]["og_id"] == "eggnog:YYYY@1224"
        assert groups[1]["specificity_rank"] == 2

    def test_thermosynechococcus_gets_cyanobacteria(self):
        """Thermosynechococcus gene gets bacteria + Cyanobacteria."""
        gene = {
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
                "3Q4X@1117|Cyanobacteria",
            ],
        }
        groups = extract_ortholog_groups(gene, "Thermosynechococcus")
        assert len(groups) == 2
        assert groups[1]["og_id"] == "eggnog:3Q4X@1117"
        assert groups[1]["taxonomic_level"] == "Cyanobacteria"
        assert groups[1]["specificity_rank"] == 2

    def test_meiothermus_bacteria_only(self):
        """Meiothermus gene gets only bacteria-level (empty level list)."""
        gene = {
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
                "XXXX@1297|Deinococcus-Thermus",
            ],
        }
        groups = extract_ortholog_groups(gene, "Meiothermus")
        assert len(groups) == 1
        assert groups[0]["og_id"] == "eggnog:COG0592@2"

    def test_pseudomonas_gets_proteobacteria(self):
        """Pseudomonas gene gets bacteria + Proteobacteria."""
        gene = {
            "eggnog_ogs": [
                "COG0592@2|Bacteria",
                "ZZZZ@136845|Pseudomonas putida group",  # NOT configured — skipped
                "PPPP@1224|Proteobacteria",
            ],
        }
        groups = extract_ortholog_groups(gene, "Pseudomonas")
        assert len(groups) == 2
        assert groups[1]["og_id"] == "eggnog:PPPP@1224"
```

- [ ] **Step 2: Run tests**

Run: `pytest tests/test_ortholog_group_extraction.py -v`
Expected: ALL pass.

- [ ] **Step 3: Commit**

```bash
git add tests/test_ortholog_group_extraction.py
git commit -m "test: add ortholog group extraction tests for Shewanella, Ruegeria, Thermosynechococcus, Meiothermus, Pseudomonas"
```

---

### Task 5: Update stats logging in `build_gene_annotations.py`

**Files:**
- Modify: `multiomics_kg/download/build_gene_annotations.py:845-863`

The current stats distinguish `eggnog_bacteria` vs `eggnog_lowest`. With multiple intermediate levels, change to count by specificity_rank.

- [ ] **Step 1: Update the OG stats section**

Replace lines 845-863:

```python
    og_stats = {"has_og": 0, "cyanorak": 0, "eggnog_bacteria": 0, "eggnog_intermediate": 0, "eggnog_family": 0}
    for locus_tag, gene in merged_out.items():
        ogs = extract_ortholog_groups(gene, og_organism_group)
        if ogs:
            gene["ortholog_groups"] = ogs
            og_stats["has_og"] += 1
            for og in ogs:
                if og["source"] == "cyanorak":
                    og_stats["cyanorak"] += 1
                elif og["taxon_id"] == 2:
                    og_stats["eggnog_bacteria"] += 1
                elif og["specificity_rank"] == 2:
                    og_stats["eggnog_intermediate"] += 1
                else:
                    og_stats["eggnog_family"] += 1
    n = stats["total"] or 1
    print(f"\n  === {strain_name} ortholog groups ===")
    print(f"  Genes with OG:        {og_stats['has_og']} ({100 * og_stats['has_og'] // n}%)")
    print(f"  Cyanorak memberships: {og_stats['cyanorak']}")
    print(f"  EggNOG bacteria:     {og_stats['eggnog_bacteria']}")
    print(f"  EggNOG intermediate: {og_stats['eggnog_intermediate']}")
    print(f"  EggNOG family:       {og_stats['eggnog_family']}")
```

- [ ] **Step 2: Run full unit test suite**

Run: `pytest -m "not slow and not kg" -v`
Expected: ALL pass.

- [ ] **Step 3: Commit**

```bash
git add multiomics_kg/download/build_gene_annotations.py
git commit -m "feat: update OG stats logging for intermediate bridging levels"
```

---

### Task 6: Rebuild data and verify

- [ ] **Step 1: Rebuild gene annotations**

Run: `bash scripts/prepare_data.sh --steps 2`

This re-runs `build_gene_annotations.py` for all strains, populating the new intermediate-level OGs in `gene_annotations_merged.json`.

- [ ] **Step 2: Spot-check Shewanella annotations**

```bash
uv run python -c "
import json
with open('cache/data/Shewanella/genomes/W3-18-1/gene_annotations_merged.json') as f:
    genes = json.load(f)
with_proteo = sum(1 for g in genes.values()
    if any(og.get('taxon_id') == 1224 for og in g.get('ortholog_groups', [])))
print(f'Shewanella genes with Proteobacteria OG: {with_proteo}/{len(genes)}')
"
```

Expected: ~3,900 out of 4,058 genes have Proteobacteria-level OGs.

- [ ] **Step 3: Spot-check Thermosynechococcus annotations**

```bash
uv run python -c "
import json
with open('cache/data/Thermosynechococcus/genomes/BP1/gene_annotations_merged.json') as f:
    genes = json.load(f)
with_cyano = sum(1 for g in genes.values()
    if any(og.get('taxon_id') == 1117 for og in g.get('ortholog_groups', [])))
print(f'Thermosynechococcus genes with Cyanobacteria OG: {with_cyano}/{len(genes)}')
"
```

Expected: ~2,274 out of 2,449 genes have Cyanobacteria-level OGs.

- [ ] **Step 4: Spot-check Prochlorococcus now gets both levels**

```bash
uv run python -c "
import json
with open('cache/data/Prochlorococcus/genomes/MED4/gene_annotations_merged.json') as f:
    genes = json.load(f)
both = sum(1 for g in genes.values()
    if any(og.get('taxon_id') == 1212 for og in g.get('ortholog_groups', []))
    and any(og.get('taxon_id') == 1117 for og in g.get('ortholog_groups', [])))
print(f'MED4 genes with BOTH Prochloraceae AND Cyanobacteria OGs: {both}/{len(genes)}')
"
```

Expected: A significant number of MED4 genes now have both levels (previously either/or).

- [ ] **Step 5: Commit data verification notes**

No files to commit — this is a verification step only.

---

### Task 7: Update CLAUDE.md

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Update the ORGANISM_GROUP_LEVELS reference in CLAUDE.md**

Find the line in CLAUDE.md that says:

```
- OrthologGroup nodes: ~21K nodes
```

Update the count estimate and the description to mention 8 organism groups and the bridging levels. The exact counts will come from the rebuild, but the description should reflect that Proteobacteria and Cyanobacteria bridging levels are now extracted for diverse organisms.

Also update the bullet about `Gene_in_ortholog_group` edges to reflect the increased count (~100-110K expected).

- [ ] **Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: update CLAUDE.md with expanded ortholog group bridging levels"
```
