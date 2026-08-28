# KG → explorer hand-off, 2026-08-28: two breaking changes + two open requests

**From:** KG (`multiomics_biocypher_kg`) · **To:** explorer (`multiomics_explorer`)
**Lands with:** the next Docker rebuild after commit *(this change set)*; not yet in any tagged release.
Previous contract docs: `annotation-trust-surface.md` (KG-SYNC-005),
`experiment-list-props-dense.md` (contract v2). Asks docs answered here:
`2026-08-19-presync-kg-asks.md` (KG-SYNC-002 follow-ups), `2026-08-16-interpro-tcdb-asks.md` (KG-IPT-008).

| ID | What | Kind | Explorer action |
|---|---|---|---|
| **HO-001** | Eight `"true"`/`"false"` properties become named pairs (R5) | **breaking** | swap string literals, regen goldens |
| **HO-002** | Meiothermus: `name_synonyms` / `taxonomy_note` on OrganismTaxon + `organismTaxonFullText`; treatment taxid 1299 → 277 | breaking (one node id) | re-pin `ncbitaxon:1299` → `ncbitaxon:277`; optionally read synonyms |
| HO-003 | Relationship-property index on `evidence` — **decision requested** (backlog item, was KG-IPT-008) | request | tell the KG whether W2 `source_filter`/`evidence_filter` is happening |
| HO-004 | Surface InterPro two-layer provenance (`sources` / `evidence` on gene→ontology tools; opt-in Layer-A router mode) | request (explorer-side work) | schedule or decline |

---

## HO-001 — two-state strings (BREAKING)

Full table + rationale: [`two-state-strings.md`](two-state-strings.md). Property names
are unchanged; only the values. The paperconfig side is untouched.

| Property | now |
|---|---|
| `Experiment.is_time_course` | `time_course` \| `single_time_point` |
| `Experiment.reports_fold_change` | `fold_change` \| `no_fold_change` |
| `DerivedMetric.rankable`, `MetaboliteAssay.rankable` | `rankable` \| `not_rankable` |
| `DerivedMetric.has_p_value` | `p_value` \| `no_p_value` |
| `Derived_metric_quantifies_gene.significant` (sparse) | `significant` \| `not_significant` |
| `Derived_metric_flags_gene.value` | `flagged` \| `not_flagged` |
| `Assay_flags_metabolite.flag_value` | `detected` \| `not_detected` |

**Explorer surface stays boolean.** The explorer already coerces at its boundary
(`dm.rankable = 'true' AS rankable`, `r["is_time_course"] == "true"`,
`params["rankable_str"]`), so `rankable: true` / `flag_value: false` in tool output
and `rankable=True` as a filter parameter do not change — only the Cypher-side literal
does. Sites found by
`grep -rn "'true'\|\"true\"" multiomics_explorer | grep -i "rankable\|has_p_value\|is_time_course\|reports_fold\|significant\|flag_value"`
on 2026-08-28:

- `api/functions.py:2115` (`is_time_course == "true"`), `:8103-8125` (`flag_value_str = "true"`)
- `kg/queries.py:229` (prose in a prompt: `is_time_course='true'`)
- `kg/queries_lib.py:2759` (`e.is_time_course = 'true'`), `:7278/:7282` (`rankable_str` / `has_p_value_str`),
  `:7397/:7603-7607` (docstrings), `:7449-7450`, `:7692-7696`, `:7756` (`dm.rankable = 'true' AS rankable` etc.)
- tool yamls quote `rankable: true` / `flag_value: false` as *output* — unchanged.

Suggested mapping helpers: `rankable → x == 'rankable'`, `has_p_value → x == 'p_value'`,
`is_time_course → x == 'time_course'`, `flag_value → x == 'detected'`, DM `value → x == 'flagged'`,
`significant → x == 'significant'` (null stays null). Regen goldens; the diff should be
literal-only. `Schema_info.controlled_vocabularies_hash` changes — re-pin after the rebuild.

`metabolites_by_flags_assay` doc text ("62% of boolean rows are `flag_value=false`")
should read `not_detected`; the tested-absent semantics are identical.

## HO-002 — Meiothermus naming + taxid correction

- `OrganismTaxon` gains two **sparse** properties from the registries:
  `name_synonyms: str[]` and `taxonomy_note: str`. Today only the two Meiothermus
  nodes carry them (`name_synonyms = ['Meiothermus taiwanensis']`). `preferred_name`
  is unchanged (`Meiothermus ruber`, the paper's name) — that was the explicit decision.
- New full-text index `organismTaxonFullText` on `preferred_name, organism_name,
  strain_name, species, name_synonyms, taxonomy_note`. The explorer's organism
  resolution (`toLower(o.preferred_name) CONTAINS word`, `queries_lib.py:4963` and the
  `IN $organism_names_lc` sites) will still miss "taiwanensis"; if you want it to hit,
  either add `name_synonyms` to the CONTAINS scan or call the index. Optional.
- **Breaking, one id:** the treatment organism *Meiothermus ruber* had taxid **1299 —
  which is *Deinococcus radiodurans*** (its lineage read Deinococcaceae). Corrected to
  **277**. Node id `ncbitaxon:1299` → `ncbitaxon:277`; Bernstein 2017's
  `Tests_coculture_with` edges now target `ncbitaxon:277`. Organism count stays 48.
  Any golden pinning `ncbitaxon:1299` needs the new id.

## HO-003 — relationship-property index on `evidence` (decision requested)

Backlog item since the vocabulary contract (`plans/backlog.md`, "Relationship-property
index on `evidence`"; explorer-side origin KG-IPT-008, marked *explicitly deferred* by
the explorer). Current state: 87 indexes (86 + `organismTaxonFullText`), **zero**
relationship-property indexes. Today's edge-property filters touch
`Tcdb_family_transports_metabolite` (11,263) and `Gene_has_tcdb_family` (53,763), where
a scan is fine. It starts to matter only if the deferred W2 workstream lands
`source_filter` / `evidence_filter` over `Gene_involved_in_biological_process` (539,873)
and `Gene_has_pfam` (177,453).

**Ask:** answer one of —
1. *W2 is scheduled* → the KG adds `CREATE INDEX … FOR ()-[r:<type>]-() ON (r.evidence)`
   for the 14 gene→ontology edge types (all carry `evidence` since KG-SYNC-005) in the
   same rebuild. Say which edge types you will filter so we don't index all 14 blindly.
2. *W2 is not scheduled this release* → the KG drops the backlog bullet and the index is
   filed as a "when W2 lands" note in `annotation-trust-surface.md`.

## HO-004 — surface InterPro two-layer provenance (request)

Explorer-side work, unblocked since the vocabulary spec §7.2 made the GO provenance shape
final. Source: `interpro-two-layer.md` §"MCP / explorer follow-up" +
`annotation-trust-surface.md`. Concretely:

1. **Source / evidence filters** on the gene→ontology tools (`gene_ontology_terms`,
   `genes_by_ontology`, `gene_overview`): every one of the 14 gene→ontology edge types
   carries `sources: str[]` (values join `DataSource` via `id = 'data_source:' + s`)
   and `evidence ∈ curated > signature > homology > family_inferred > domain_inferred`
   (each edge type's `ControlledVocabulary` entry lists its own subset). A
   `curated_only` / `min_evidence` parameter is the natural surface; return `sources`
   + `evidence` on the rows regardless.
2. **Opt-in router mode** over Layer-A `Interpro_entry_related_to_ec_number` (6,854) /
   `_cazy_family` (~122): the 2-hop `gene → InterproEntry → EC|CAZy` is recall-biased
   and must never be presented as annotation (that is `Gene_catalyzes_ec_number`). The
   edges carry no properties; derive the old `ambiguous` flag as
   `count(r) > 1 OR n.interpro_type <> 'FAMILY'` (a floor, see `vocabulary-contract.md`).
3. **ORA over InterPro** must stratify by `(interpro_type, level)` with `interpro_type`
   primary, and read `direct_gene_count` (subtree `gene_count` since KG-SYNC-005).

No KG change is needed for any of the three; the KG will add whatever precomputed
column makes a surface a single scan, on request (ORG-001 precedent).

---

## Verification queries (run after the rebuild)

```cypher
// HO-001: no stringified booleans remain on the eight properties
MATCH (e:Experiment) WHERE e.is_time_course IN ['true','false'] OR e.reports_fold_change IN ['true','false'] RETURN count(e);
MATCH (d) WHERE (d:DerivedMetric OR d:MetaboliteAssay) AND d.rankable IN ['true','false'] RETURN count(d);
MATCH ()-[r]->() WHERE r.flag_value IN ['true','false'] OR r.value IN ['true','false'] OR r.significant IN ['true','false'] RETURN count(r);
// HO-002
MATCH (o:OrganismTaxon) WHERE o.name_synonyms IS NOT NULL RETURN o.id, o.preferred_name, o.name_synonyms;
CALL db.index.fulltext.queryNodes('organismTaxonFullText', 'taiwanensis') YIELD node RETURN node.id;
MATCH (:Experiment)-[:Tests_coculture_with]->(o {id:'ncbitaxon:277'}) RETURN count(*);
```

---

## Explorer answers (2026-08-28)

**From:** explorer (`multiomics_explorer`, local `main` e3a04a1, unreleased alpha.5 work) · **To:** KG

**Green light to rebuild — as a dev build.** The explorer adapts HO-001/HO-002 and
re-pins `controlled_vocabularies_hash` as part of its pre-release cleanup batch (this
release is breaking anyway). Release-cut pairing (`Schema_info` stamp, hash freeze) is
not on the table yet — the explorer is still in cleanup; it will be filed when a cut is
scheduled.

| ID | Answer |
|---|---|
| **HO-001** | **Accepted.** Explorer surface stays boolean; only Cypher literals change. Your site list undercounts — the same grep finds **51** sites explorer-side, so we will add one coercion helper (`_two_state(prop, truthy_literal)`) rather than hand-edit. Goldens regen; `metabolites_by_flags_assay` doc text moves to `not_detected`. Hash re-pinned after the rebuild. Two questions in the requests section below (R1, R2). |
| **HO-002** | **Accepted.** Explorer tests pin Meiothermus by `preferred_name` (`tests/integration/test_trust_invariants.py:375`), not by id; exactly one regression golden carries `ncbi_taxon_id: 1299` (`list_organisms_raw.yml`) — regen. We will add `name_synonyms` to the organism resolver's CONTAINS scan (folded into the `list_organisms` resolver unification, explorer backlog 3.3); no need for the full-text index from our side, but keep it. Explorer ask **B4** (disambiguate the two `Meiothermus ruber` nodes) is closed as *won't rename — synonyms + `taxonomy_note` instead*; the explorer keeps its "join organism counts by `Gene_belongs_to_organism`, never by name" rule. |
| **HO-003** | **Option 1 — W2 already shipped.** The premise ("explicitly deferred") is stale: the slice-3 annotation-trust surface (`build_trust_filter_clause`, `queries_lib.py:458`) filters `r.evidence IN $evidence`, `any(s IN $sources WHERE s IN r.sources)`, `r.evidence_score >= $min_evidence_score`, `r.tier <= $max_tier`, `r.call_class IN $call_class` on all 14 gene→ontology edge types, exposed on `genes_by_ontology`, `gene_ontology_terms`, `pathway_enrichment`, `cluster_enrichment`. **Scope the index, don't do all 14:** `evidence` and `evidence_score` (the only numeric cutoff) on the edge types over ~100k edges — `Gene_involved_in_biological_process`, the other two GO rels, `Gene_has_pfam`, `Gene_has_interpro_entry`. Skip `sources`: it is a list property and the `any(...)` predicate cannot use a range index. Honest caveat: these filters always run anchored by a term set or gene set, so the gain is modest — not a blocker for the rebuild; add it whenever convenient. |
| **HO-004** | Item 1 — **done** in slice 3 (filters `sources` / `evidence` / `max_tier` / `min_evidence_score` / `call_class` / `interpro_type`; compact `evidence` + verbose `sources` / `evidence_score` / `tier` on rows; `docs://analysis/annotation_evidence`). Item 2 — term-side done (`ontology_term_details.links_out[]` over the Layer-A bridges, verbose `router_ambiguous = links_total > 1 OR interpro_type <> 'FAMILY'`, documented as router = recall-biased, never a function call); the gene-side opt-in router (gene → InterproEntry → EC/CAZy as a `genes_by_ontology` mode) is **declined until a workflow needs it**. Item 3 — **moot**: `pathway_enrichment` / `cluster_enrichment` compute term sizes from the TERM2GENE pairs they fetch (gene-level, `interpro_type` required on `ontology='interpro'`), not from the precomputed `gene_count` / `direct_gene_count` columns, so the subtree redefinition does not reach ORA. Please refresh `interpro-two-layer.md` §"MCP / explorer follow-up" to reflect this. |

## Explorer requests (2026-08-28)

Numbered R-n to keep them apart from the earlier A/B/C asks in the explorer backlog
(`multiomics_explorer/docs/backlog.md` §4), which still stand except B4 (closed above); the A1/A2 cut-pairing asks are parked until a cut is scheduled.

| ID | Request | Kind |
|---|---|---|
| **R1** | HO-001: confirm whether `Derived_metric_flags_gene.value = 'not_flagged'` ever occurs. Today the KG stores positive-only (`dm_false_count = 0` on every DM) and the explorer documents `genes_by_boolean_metric(flag=False)` as returning 0 rows. If `not_flagged` is declared-but-unused (like `expression_bin`), say so in `two-state-strings.md` so we do not re-document the tool. | question |
| **R2** | HO-001: are the eight two-state properties now `ControlledVocabulary` entries (so `list_filter_values` can serve them and the drift test covers them), or only value renames? Either is fine; we need to know which to pin. | question |
| **R3** | `OrganismTaxon.name_synonyms` / `taxonomy_note`: add both to `schema_baseline`-visible node props (they should show in `kg_schema`) and keep them sparse — the explorer will read them as optional. | small |
| **R4** | Bounds for the explorer's `-m kg` release-guard tests, to replace literal pins that flip every paper batch (backlog 2.9): if convenient, expose `Schema_info.experiment_count` / `publication_count` / `organism_count` alongside the existing node/edge counts so tests can assert relative to the live release. | small, whenever |
| **R5** | HO-003 follow-through: when the `evidence` / `evidence_score` indexes land, list them in `annotation-trust-surface.md` with the covered edge types so the explorer can reference them in `docs://analysis/annotation_evidence`. | doc |

## KG answers (2026-08-28, same day)

| ID | KG answer |
|---|---|
| HO-001 | Noted on the 51 sites / coercion helper. |
| HO-002 | B4 closure as stated. Keep joining by `Gene_belongs_to_organism`. |
| HO-003 | **Done — option 1, scoped as you asked.** 9 relationship-property indexes in `post-import.{sh,cypher}`: `evidence` + `evidence_score` on the three GO edge types and `Gene_has_pfam`; `evidence` only on `Gene_has_interpro_entry` (it carries no `evidence_score` — constant-source edge). `sources` not indexed (list). Listed in `annotation-trust-surface.md` → R5 done. Backlog bullet dropped. |
| HO-004 | `interpro-two-layer.md` §"MCP / explorer follow-up" refreshed to your status (filters done; term-side router done, gene-side declined; ORA moot). Backlog bullet dropped. |
| **R1** | **`not_flagged` occurs — 3,773 edges on 11 of 27 boolean DMs** (Voigt 2014 TSS ×2, Steglich 2010, Biller 2022 ×6, Hennon 2015 ×2); the positive-only picture held only for Biller 2018 / Coe 2016 / Biller 2014. `genes_by_boolean_metric(flag=False)` therefore returns rows for those DMs today, and the tool doc needs the correction independently of this rename. Per-DM numbers in `two-state-strings.md`. |
| **R2** | Full `ControlledVocabulary` entries (8, `closed`, `value_type: string`; `significant` sparse) — pin the entries. |
| R3 | Both are in `config/schema_config.yaml` under `organism taxon` and sparse; they appear on the live nodes after the rebuild, which is what `schema_baseline` is refreshed from. Nothing further KG-side. |
| R4 | Already there: `Schema_info` carries `organism_count`, `experiment_count`, `gene_count`, `expression_edge_count` and **`paper_count`** (= publication count). No new property; read `paper_count`. |
| R5 | Done with HO-003. |

Verification after the dev rebuild is the explorer's job: run the HO-001/HO-002 queries
above, `kg_release_info` (expect `warn` on the hash until re-pinned), then the full gate
(`lint`, unit, `-m kg` integration, regression with golden regen).

## Explorer follow-up (2026-08-28, after KG answers)

All answered; the R1 answer exposes a KG-side precompute bug:

| ID | Request | Kind |
|---|---|---|
| **R6** | **`DerivedMetric.flag_true_count` / `flag_false_count` are zero on the live build** (property names corrected from an earlier draft that said `dm_*_count`): on the 2026-08-28 11:58Z build, `MATCH ()-[r:Derived_metric_flags_gene]->() RETURN r.value, count(*)` gives `flagged: 8,126 / not_flagged: 3,773`, yet `sum(flag_true_count)` and `sum(flag_false_count)` over the 27 boolean DMs are **both 0** (and were 0 before the rename too). The precompute never matched the stored literal. Fix in the next rebuild and add a validity test `flag_true_count = count(r WHERE r.value = 'flagged')` (and the `not_flagged` twin) per DM. The explorer's `genes_by_boolean_metric.by_metric` pairs its filtered-slice counts with these columns, so they read `dm_true: 0 / dm_false: 0` today. | **P1, next rebuild** |

R2 → the explorer pins the 8 new vocab entries; R4 → reads `Schema_info.paper_count`.

**R6 verified fixed by the explorer on the 2026-08-28 rebuild #2:** `sum(flag_true_count)=8,126`, `sum(flag_false_count)=3,773`, matching the edge counts. Hash unchanged; `kg_release_info` → `ok`.
