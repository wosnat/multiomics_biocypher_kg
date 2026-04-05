# OrthologGroup Node Enrichment Plan

## Goal

Add computed aggregate properties to OrthologGroup nodes:

| Property | Type | Description |
|---|---|---|
| `consensus_product` | str | Majority-vote product string from member genes |
| `consensus_gene_name` | str | Most frequent non-null `gene_name` among members |
| `member_count` | int | Total genes in group |
| `organism_count` | int | Distinct organisms represented |
| `genera` | str[] | List of distinct genera (e.g. `["Prochlorococcus", "Alteromonas"]`) |
| `has_cross_genus_members` | bool | Quick flag for groups bridging genera |

Also add: `specificity_rank` index in post-import.

---

## Decision: Adapter (Python) vs Post-Import (Cypher)

### Option A: Compute in adapter (Python) during build

**Pros:**
- Already have all member data in memory (`gene_annotations_merged.json` per strain)
- Full Python stdlib for aggregation (Counter, majority vote, genus extraction)
- Properties appear in CSVs → `neo4j-admin import` loads them natively, no extra step
- Easy to unit test with synthetic fixtures (existing test pattern)
- No dependency on running Neo4j instance
- Can handle complex logic cleanly (e.g., filtering "hypothetical protein" from consensus, tie-breaking)

**Cons:**
- Requires two-pass restructure of `get_nodes()`: first collect all members per OG, then emit nodes
- `get_nodes()` currently iterates adapters once; needs to iterate twice (once for aggregation, once for node emission) — or refactor into a single pass that collects members, then emits
- `genera` extraction requires parsing `organism_strain` strings (but this is trivial: first word)

### Option B: Compute in post-import.sh (Cypher)

**Pros:**
- No adapter code changes
- Declarative — easy to read as a single Cypher statement
- Could be added/removed without rebuilding CSVs

**Cons:**
- Cypher aggregation for majority-vote is awkward (nested `collect`, `unwind`, `count`, `ORDER BY`, `LIMIT 1` subqueries)
- `genera` as a list property requires string manipulation in Cypher (`split(g.organism_strain, ' ')[0]`)
- Post-import runs inside Docker via `cypher-shell` — harder to debug, no unit tests
- Adds runtime to post-import (currently just indexes, ~1s; aggregation over 84K edges + 21K nodes adds several seconds)
- Properties not in CSV schema → schema_config.yaml still needs updating, but BioCypher won't validate them at import time
- Two files to keep in sync (`post-import.sh` + `post-import.cypher`)
- `has_cross_genus_members` would need a second pass or subquery per OG

### Recommendation: **Option A (adapter)**

The data is already loaded in Python. The aggregation is simple (Counter + max). Unit testing is straightforward with the existing fixture pattern. The only cost is a minor refactor to two-pass node emission, which is clean and maintainable.

---

## Execution Plan

### Phase 1: Schema + Post-Import

**Files changed:**
- `config/schema_config.yaml` — add 6 new properties to `ortholog group`
- `scripts/post-import.sh` — add `specificity_rank` index
- `scripts/post-import.cypher` — same index (keep in sync)

**Acceptance criteria:**
- New properties declared in schema
- Index added to both post-import files

### Phase 2: Adapter Implementation

**Files changed:**
- `multiomics_kg/adapters/ortholog_group_adapter.py`

**Changes:**
1. Add `_clean_str()` helper (same as other adapters)
2. Modify `OrthologGroupAdapter` to also expose gene metadata:
   - New method `get_og_memberships_with_gene_data()` returning `(locus_tag, og_dict, gene_meta)` where `gene_meta = {product, gene_name, organism_strain}`
   - Keep existing `get_og_memberships()` unchanged (backward compat for edges)
3. Refactor `MultiOrthologGroupAdapter.get_nodes()`:
   - First pass: collect all members per OG ID into a dict `{og_id: {og_dict, members: [{product, gene_name, organism_strain}, ...]}}`
   - Second pass: for each OG, compute consensus properties and emit node
4. Consensus logic:
   - `consensus_product`: most common non-null product, excluding "hypothetical protein" and "conserved hypothetical protein" (fall back to those if nothing else). Use `_clean_str()`.
   - `consensus_gene_name`: most common non-null `gene_name`. Use `_clean_str()`.
   - `member_count`: `len(members)`
   - `organism_count`: `len(set(organism_strain for m in members))`
   - `genera`: `sorted(set(organism_strain.split()[0] for m in members if organism_strain))`
   - `has_cross_genus_members`: `len(genera) > 1`

**Acceptance criteria:**
- All 6 new properties populated on every OrthologGroup node
- `consensus_product` and `consensus_gene_name` can be `None` (if all members lack those fields)
- String properties sanitized (no `'` or `|`)
- `get_edges()` unchanged
- Existing `get_og_memberships()` API unchanged

### Phase 3: Unit Tests

**Files changed:**
- `tests/test_ortholog_group_adapter.py`

**New tests:**
1. `test_consensus_product_majority_vote` — 3 members with same product, 1 different → majority wins
2. `test_consensus_product_excludes_hypothetical` — members split between "hypothetical protein" and a real product → real product wins
3. `test_consensus_product_all_hypothetical` — falls back to "hypothetical protein" when nothing else
4. `test_consensus_gene_name` — most frequent non-null gene_name wins
5. `test_member_count` — matches expected count per OG
6. `test_organism_count` — correct distinct organism count
7. `test_genera_list` — correct genera extracted, sorted
8. `test_has_cross_genus_members_true` — OG spanning Pro + Alt → `True`
9. `test_has_cross_genus_members_false` — OG within single genus → `False`
10. `test_consensus_fields_none_when_all_null` — all members missing gene_name → `None`
11. `test_string_sanitization` — product with `'` or `|` → cleaned

**Fixture changes:**
- Add `product`, `gene_name`, `organism_strain` fields to `STRAIN1_DATA` and `STRAIN2_DATA`
- Add a third strain fixture with overlapping OGs for cross-genus testing

**Run:**
```bash
pytest tests/test_ortholog_group_adapter.py -v
```

**Acceptance criteria:**
- All new tests pass
- All existing tests still pass (no regressions)

### Phase 4: KG Validity Tests

**Files changed:**
- `tests/kg_validity/test_post_import.py`

**New tests (require running Neo4j):**
1. `test_ortholog_group_has_member_count` — all OGs have `member_count > 0`
2. `test_ortholog_group_has_organism_count` — all OGs have `organism_count > 0`
3. `test_ortholog_group_member_count_matches_edges` — `member_count` = actual edge count per OG
4. `test_ortholog_group_cross_genus_flag_consistent` — `has_cross_genus_members` matches actual genera count
5. `test_ortholog_group_consensus_product_not_all_null` — majority of OGs have a consensus_product
6. `test_specificity_rank_index_exists` — verify the new index was created

**Run:**
```bash
pytest tests/kg_validity/test_post_import.py -v
```

### Phase 5: Code Review Checklist

- [ ] Schema properties match adapter property keys exactly (string literals)
- [ ] `_clean_str()` applied to `consensus_product` and `consensus_gene_name`
- [ ] `genera` is a list (str[]), not a joined string — verify BioCypher handles array delimiter `|`
- [ ] No `'` or `|` in any new string values
- [ ] `post-import.sh` and `post-import.cypher` are identical in Cypher content
- [ ] Existing node properties (`name`, `source`, `taxonomic_level`, `taxon_id`, `specificity_rank`) unchanged
- [ ] Edge output unchanged (no new edge properties, same count)
- [ ] `get_og_memberships()` API backward-compatible
- [ ] Test fixtures updated consistently
- [ ] CLAUDE.md updated: key graph facts section mentions new properties

### Phase 6: Documentation

**Files changed:**
- `CLAUDE.md` — update "Key graph facts" with new OrthologGroup properties
- `docs/explorer_ortholog_migration.md` — mention enrichment properties

---

## Estimated Scope

- ~60 lines adapter code (aggregation logic + _clean_str)
- ~6 lines schema_config.yaml
- ~2 lines post-import (index)
- ~120 lines new unit tests
- ~60 lines new KG validity tests
