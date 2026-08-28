# Annotation-trust surface (KG-SYNC-005, 2026-08-27)

**Driver:** `multiomics_explorer/docs/kg-specs/2026-08-27-annotation-trust-kg-asks.md` (asks ONT-001…015;
KG verdicts §6, explorer acceptance §7, KG reply §8). **Plan:** `plans/kg_sync_005_annotation_trust.md`.
Ships before the 0.1.0-alpha.7 cut; the interpro / ncbifam / merops ontologies are unreleased, so the two
renames below are free.

## What changed

### 1. `sources` + `evidence` on every gene→ontology edge (ONT-007/008)

| edge type | `sources` | `evidence` | note |
|---|---|---|---|
| `Gene_involved_in_biological_process` / `_enables_molecular_function` / `_located_in_cellular_component` | per-token (unchanged) | curated / family_inferred / domain_inferred | unchanged |
| `Gene_catalyzes_ec_number`, `Gene_has_pfam`, `Gene_has_cazy_family` | per-token (unchanged) | unchanged | unchanged |
| `Gene_has_tcdb_family` | unchanged | **new**: `homology` (diamond called it) / `family_inferred` (eggNOG-only); both-source → `homology` | derived from `sources` in `tcdb_adapter._tcdb_evidence` |
| `Gene_has_kegg_ko`, `Gene_in_cog_category` | **new** `['eggnog']` | **new** `family_inferred` | constant — eggNOG-only fields |
| `Gene_has_cyanorak_role`, `Gene_has_tigr_role` | **new** `['cyanorak']` | **new** `curated` | constant |
| `Gene_has_interpro_entry`, `Gene_has_ncbifam_family` | **new** `['interproscan']` | **new** `signature` | constant; member DBs stay in `libraries` |
| `Gene_has_merops_family` | **new** `['merops_diamond']` | **new** `homology` | constant |

One ladder: **`curated > signature > homology > family_inferred > domain_inferred`**. `homology` is new —
a direct sequence hit against a curated reference DB; strength *within* it is `tier`. It is deliberately not
split by tier into more rungs, and it is orthogonal to `call_class` (a `nonpeptidase_homolog` MEROPS edge is
still `homology` evidence for the *placement*). Each edge type's vocab entry lists its own possible subset.

The explorer's ask assumed a per-token provenance map for KO/COG/roles; none exists (the merge tracks
`<field>_source` only for `ec_numbers`, `go_terms`, `pfam_ids`, `cazy_ids`, `ncbifam_ids`, `gene_name`,
`product`, `function_description`) and none is needed — those fields have exactly one contributing source.

### 2. MEROPS `evidence_score` (ONT-006/013)

`Gene_has_merops_family.evidence_score` (post-import, float `{0, 0.5, 1}`; R4 shape, `signal_count = 2`,
`signals = [tier_le_2, pfam_support]`). `call_class = 'peptidase'` is **not** a signal: it is the verdict
axis, and scoring it would rank a confidently placed dead homolog lower for being honest. Sparse
`Gene.merops_evidence_score_max` mirrors `tcdb_evidence_score_max` (`coalesce(x, -1.0)` for a total order).

### 3. `Gene_has_tcdb_family.attachment_depth` (ONT-010)

`most_specific | superseded`, materialized once in post-import; `superseded` = another edge of the SAME
gene lands on a descendant (`*1..4`) of this node — a less specific call, not a wrong one. The three
transport-arm rollups (`Gene.transported_metabolite_count`, `Metabolite.transporter_gene_count`,
`Organism_has_metabolite`) now read the property instead of re-deriving the predicate (verified
byte-identical counts against the pre-batch `post-import-validate.sh` dump). MED4: 670 edges, 73 superseded.

### 4. One `gene_count` semantics (ONT-009/015)

| label | `gene_count` / `organism_count` | `direct_gene_count` |
|---|---|---|
| BiologicalProcess, MolecularFunction, CellularComponent | **new**, subtree over `is_a ∪ part_of` | **new** |
| EcNumber, KeggTerm (all 4 levels), CyanorakRole | **new**, subtree | **new** |
| InterproEntry | **changed** DIRECT → subtree | **new** (= the old `gene_count`) |
| TcdbFamily, CazyFamily, MeropsFamily | subtree (unchanged) | **new** |
| Pfam, TigrRole, CogFunctionalCategory | **new**, direct | — (flat) |
| PfamClan | **new**, via member Pfams | — (would be constant 0) |
| BriteCategory, NcbifamFamily, SubcellularLocalization, SignalPeptideType | unchanged | — |

GO and KEGG are DAGs: a gene reachable through two parents is counted once per node, so sibling counts
do not sum to the parent's. Descendants are collected `DISTINCT` before genes are matched, so DAG path
multiplicity never multiplies rows. `MeropsFamily.peptidase_organism_count` added (ONT-005).

### 5. Hygiene (ONT-001/002/003/004/011/012/014)

- `Gene_has_ncbifam_family.match_count` **struck from the contract** — it was never emitted (adapter and
  schema only ever carried `start/end/evalue/score`), and is meaningful only on the cross-library InterPro edge.
- `Gene_has_ncbifam_family.score` → **`bit_score`** (`score` is PSORTb's confidence on
  `Gene_has_subcellular_localization`, a released and unrelated scale).
- `MeropsFamily.family_type` → **`family_class`** (R1b: `NcbifamFamily.family_type` holds an unrelated
  external vocabulary under the same name).
- 10 retired `NcbifamFamily` nodes: `family_type = 'retired'` (single KG-minted sentinel in an external set).
- ~50 new `ControlledVocabulary` entries: `sources`/`evidence` ×7+1, `evidence_score`/`attachment_depth`/
  `merops_evidence_score_max`, and the numeric edge props the explorer filters on (`evalue`, `match_count`,
  `bit_score`, `start`/`end`, the five diamond props on TCDB (sparse) and MEROPS (dense), PSORTb `score`,
  SignalP `probability`/`cleavage_site`/`cleavage_probability`). The InterPro `evalue` entry names the
  member DBs that never report one (COILS, HAMAP, MOBIDB_LITE, PANTHER, PRINTS, PROSITE_PATTERNS,
  PROSITE_PROFILES, SUPERFAMILY — measured over the committed calls.json, 2026-08-27; re-measure on an
  InterProScan release bump).

## Breaking

- `InterproEntry.gene_count` semantics: DIRECT → SUBTREE (use `direct_gene_count` for ORA).
- `Gene_has_ncbifam_family.score` → `bit_score`; `MeropsFamily.family_type` → `family_class`.

## Live numbers (2026-08-27 rebuild, 124,751 genes)

`ControlledVocabulary` 110 nodes. TCDB `evidence`: homology 40,598 / family_inferred 13,165;
`attachment_depth`: most_specific 46,593 / superseded 7,170. MEROPS `evidence_score`: 0 → 151,
0.5 → 3,768, 1.0 → 338. InterPro subtree `gene_count` differs from `direct_gene_count` on 115 of 12,999
entries (max +63). GO root `biological_process`: 66,484 subtree / 20,563 direct.

## Relationship-property indexes (2026-08-28, explorer HO-003 / R5)

The explorer's trust filters (`build_trust_filter_clause`: `r.evidence IN $evidence`,
`r.evidence_score >= $min_evidence_score`, `r.tier <= $max_tier`, `r.call_class IN
$call_class`, `any(s IN $sources WHERE s IN r.sources)`) run on all 14 gene→ontology
edge types. Only the types over ~100K edges are indexed — the rest are scanned inside
a term/gene-anchored match where an index does not pay:

| Index | Edge type | Property |
|---|---|---|
| `gene_go_bp_evidence_idx` / `gene_go_bp_evidence_score_idx` | `Gene_involved_in_biological_process` (539,873) | `evidence` / `evidence_score` |
| `gene_go_mf_evidence_idx` / `gene_go_mf_evidence_score_idx` | `Gene_enables_molecular_function` | `evidence` / `evidence_score` |
| `gene_go_cc_evidence_idx` / `gene_go_cc_evidence_score_idx` | `Gene_located_in_cellular_component` | `evidence` / `evidence_score` |
| `gene_pfam_evidence_idx` / `gene_pfam_evidence_score_idx` | `Gene_has_pfam` (177,453) | `evidence` / `evidence_score` |
| `gene_interpro_evidence_idx` | `Gene_has_interpro_entry` (404,191) | `evidence` only — constant-source edge, no `evidence_score` |

Not indexed: `sources` (list property; the `any(...)` predicate cannot use a range
index), `tier` / `call_class` (only on the 54K-edge TCDB and 4K-edge MEROPS types).
These are the KG's first relationship-property indexes.

## Verification

`tests/kg_validity/test_annotation_trust.py` (AC1–AC6 of the plan) + `pytest -m kg`; unit gates in
`tests/test_controlled_vocab.py`, `tests/test_tcdb_adapter.py`, `tests/test_cog_role_annotation_adapter.py`;
`scripts/post-import-validate.sh` byte-diff on the untouched sections; `capture_annotation_state.py --compare`
(no bucket movement expected — this batch adds no bucket).
