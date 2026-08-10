# InterPro as a gene-annotation merge source + union source attribution — design

**Date:** 2026-08-07
**Branch:** `interpro-go-pathways`
**Predecessors:**
- `docs/superpowers/specs/2026-07-22-interproscan-domains-design.md` (Phase 1 — calls.json)
- `docs/superpowers/specs/2026-07-26-interproscan-kg-integration-design.md` (Phase 2 — `InterproEntry` ontology)

**Status:** design proposed; not implemented

## 0. Intent

Promote `interproscan` from a **bookkeeping** source in the gene-annotation merge
to a **contributing annotator**, so InterPro's entry-level cross-references feed
the same ontology fields as UniProt / NCBI / Cyanorak / eggNOG.

This is deliberately *not* a new source. `gene_annotations_config.yaml` already
declares `interproscan` (§ `sources:`), contributing exactly one field —
`interpro_entries` — whose own comment describes it as feeding
`contributing_sources` + Gene routing. This spec widens what it contributes.

It also fixes a prerequisite defect: **union fields carry no provenance**
(§4). Because every field in scope here is a `union`, that fix is blocking, not
optional.

### Driving use cases (user, 2026-08-07)

| use case | what it needs | served by |
|---|---|---|
| "this hypothetical is upregulated — what is it" | any functional handle on uncharacterized ORFs | GO, Pfam, EC (§2); **sites + unintegrated signatures — deferred, §8.1** |
| "nutrient exchange between organisms" | transport, secretion, carbohydrate activity, heterotroph roles | CAZy (§2); NCBIfam GO/EC (§9.1d); existing TCDB/PSORTb/SignalP; **heterotroph roles partially — §9.1c ceiling** |
| "what does this mutation mean" | residue-level functional annotation | **not served by this spec — deferred, §8.1** |

### Two sources of entry xrefs — they complement, not duplicate

> **Superseded framing (2026-08-07).** An earlier draft argued reference-lookup
> *instead of* a `--goterms --pathways` re-scan. The re-scan has since happened
> (`b8172a29`), so both sources now exist and the question is what each is for.

**Per-strain scan artifacts** (`<strain>.interproscan.entry_xrefs.json`, from
`b8172a29`) carry GO + pathways for entries some protein actually matched:
69,759 proteins with GO, 76,377 with pathways across 42 strains.

**The reference cache** (`interpro_reference.json`, prepare_data step 9) carries
them for **every** InterPro entry — 14,803 with GO, 5,091 with pathways of
54,190. It remains necessary for two things the artifacts cannot supply:

1. **EC and CAZy.** InterProScan emits GO and pathway xrefs but **neither of
   these**. The reference is the only source — EC on 5,100 of 12,994 observed
   entries (39.2%), CAZy on 113. These are §2's largest and most
   carbon-relevant contributions respectively.
2. **is-a ancestor entries** that no gene matches directly, needed for
   ontology-level rollups.

Consequence for §2: the **GO** row is now largely redundant with the artifacts
for *gene→GO* purposes — prefer the artifacts there, since they are per-protein
and require no inference step. The reference's GO earns its place only at the
ontology level. **EC and CAZy are unaffected and remain reference-only.**

The xrefs are genuinely entry-level, which is what makes the two sources
equivalent where they overlap: over 100 MED4 proteins scanned with
`--goterms --pathways`, the reference lookup reproduced **100/100** identical GO
sets and **100/100** identical MetaCyc sets.

## 1. Agreed scope decisions

1. **Entry-level xrefs only.** These map onto existing node types and need no
   schema invention. Positional/site data and unintegrated signatures are a
   different modelling problem — separate spec (§8.1).
2. **Propagation must be type-filtered.** Entry-level xrefs inherit down to every
   gene matching the entry, so broad entries splatter. Filters are **per-field**,
   not global (§3) — the correct cut differs between EC and CAZy.
3. **Source attribution is required.** Implement `track_source` for `union`
   (§4). No field in §2 ships without it.
4. **Reference loading = lazy, inside step 2** (§5), mirroring `load_pfam_data`.
   No renumbering of `prepare_data.sh`.
5. **MetaCyc deferred** (§8.2) — no consumer, needs a new node type.
6. **NCBIfam gated** on the TIGRFAM→TigrRole question (§9.1).

## 2. Fields touched

All measured across 42 strains / 124,751 genes (2026-08-06). "ratio" =
net-new : corroborating — a low ratio means the source mostly *confirms* existing
calls, which is the signature of a safe merge.

| field | merge type | filter (§3) | net-new | corrob. | ratio | genes gaining first | verdict |
|---|---|---|---:|---:|---:|---:|---|
| `pfam_ids` | union | none needed | 13,977 | 133,249 | **0.10** | 3,438 | **ship** — signature-level, no inference |
| `go_terms` | union | none | 50,762 | 146,576 | 0.35 | 7,029 | **ship** |
| `ec_numbers` | union | FAMILY-only | 31,492 | 30,222 | 1.0 | 4,268 | **ship** |
| `cazy_ids` | union | exclude fold | 648 | 784 | 0.8 | 608 | **ship** |
| `alternate_functional_descriptions` | list append | — | — | — | — | — | ship (`[interpro]` prefix, mirrors `[pfam]`) |
| `go_terms` ← **NCBIfam** | union | — | 14,270 | 85,678 | **0.17** | 1,652 | **ship** — §9.1d |
| `ec_numbers` ← **NCBIfam** | union | — | 1,083 | 19,127 | **0.06** | 350 | **ship** — highest precision in spec, §9.1d |
| `ontology_terms` (raw NCBIfam accessions) | union | — | 67,459 | ~0 | — | 47,324 | **rejected** — no consumer, §9.1d |
| `interpro_pathways` | — | — | — | — | — | — | **deferred** — §8.2 |

Notes:

- **`pfam_ids` is the safest field here.** It is the only one sourced from a
  *direct HMM hit* (`matches[].signature_accession` where `library == "PFAM"`)
  rather than an entry-level inference, so it carries no propagation risk. A
  0.10 ratio means it is overwhelmingly corroborating eggNOG. `enrich_pfam_fields`
  already normalizes and filters to `PF*`, so it needs no new logic.
- **`go_terms` and `ec_numbers` / `cazy_ids` come from different places.** GO and
  pathways come from `interpro_reference.json`. EC and CAZy come from
  `interpro.xml.gz` db_xrefs and are **not currently in the reference cache** —
  step 9 must be extended to emit them (§6).
- `ontology_terms` currently has `filter_not: "^GO:"` and collects non-GO tokens
  from Cyanorak/NCBI. All 47,324 genes "gain their first" only because the field
  is near-empty today; there is no consumer for NCBIfam accessions. See §9.1.

## 3. Propagation filters

InterPro entry types encode a claim strength that the KG must respect:

- `FAMILY` — "this protein **is** an X" → licenses a functional assignment.
- `DOMAIN` — "this protein **contains** an X-like region".
- `HOMOLOGOUS_SUPERFAMILY` — **fold-level**; shared structure, unrelated function.

### 3.1 EC — `FAMILY` only

| filter | net-new | corrob. | ratio | first-EC genes |
|---|---:|---:|---:|---:|
| all entries | 117,542 | 34,487 | 3.4 | 10,346 |
| exclude fold | 85,945 | 34,288 | 2.5 | 8,057 |
| **FAMILY only** | **31,492** | **30,222** | **1.0** | **4,268** |
| FAMILY + single-EC | 10,750 | 22,148 | 0.5 | 3,398 |
| single-EC, any type | 31,911 | 25,084 | 1.3 | 8,017 |

`FAMILY only` is chosen because it is **semantic, not a tuned threshold**. Fold
exclusion alone is insufficient (ratio still 2.5) — the damage also comes from
multi-EC `DOMAIN` entries, e.g. `IPR006047` *Glycosyl hydrolase family 13,
catalytic domain* carrying 22 distinct ECs.

If recall is later preferred over precision, `single-EC any type` nearly doubles
first-EC genes at ratio 1.3 — but it is a cardinality heuristic, not a claim
about evidence strength.

### 3.2 EC specificity + a filter defect

InterPro EC xrefs are **76.2% complete 4-level**; 23.7% are 3-level. The partials
come in two shapes, and the existing `ec_numbers` filter treats them differently:

| form | count | vs `^\d+\.[\d\-]+\.[\d\-]+[\.\-]` |
|---|---:|---|
| `3.4.24.77` complete | 12,338 | passes |
| `2.1.1.-` dashed partial | 2,117 | passes |
| `3.4.21` **bare** partial | **1,731** | **rejected** |

The regex requires a 4th separator, so bare 3-level ECs are silently dropped
(~11% of InterPro's ECs). **Decision required:** normalize `3.4.21` → `3.4.21.-`
on ingest (recommended — a 3rd-level EC is still usable ORA signal), or accept
the loss. Either way the `normalize_ec` transform must be applied to the InterPro
source, as it is to every other EC source.

### 3.3 CAZy — exclude fold only, **keep DOMAIN**

| filter | net-new | corrob. | ratio | first-CAZy genes |
|---|---:|---:|---:|---:|
| all entries | 2,961 | 788 | 3.8 | 747 |
| **fold excluded** | **648** | **784** | **0.8** | **608** |

78% of the net-new traced to a single entry:

```
IPR017853 [HOMOLOGOUS_SUPERFAMILY] "Glycoside hydrolase superfamily"
  → CBM5, GH1, GH13, GH18, GH22   (5 unrelated families, ~457 genes each)
```

Excluding fold-level entries removes 78% of the noise, retains 81% of the
coverage gain, and loses **4** corroborating pairs.

**`DOMAIN` is deliberately kept for CAZy** — a CBM5 carbohydrate-binding module
is a genuine domain-level claim. This is why filters are per-field: the same
`DOMAIN` type that is noise for EC is signal for CAZy.

CAZy is small in absolute terms but is the largest *relative* gain in the spec:
1,308 genes carry CAZy today; +608 is a **+46% increase in CAZy gene coverage**,
directly serving carbon-metabolism questions.

## 4. Union source attribution (blocking prerequisite)

### 4.1 The defect

`track_source` is implemented **only** in `_resolve_single`.
`_resolve_union` never reads it. The `track_source: go_terms_source` declared on
`go_terms` in `gene_annotations_config.yaml` has therefore been a **silent
no-op**. Confirmed against MED4's merged JSON — the only `*_source` keys present
are `product_source`, `gene_name_source`, `function_description_source`.

Consequence today: `go_terms` is already a provenance-blind 3-source blend, and
`_compute_contributing_sources` (which walks `*_source` fields via
`_has_source_label`) **cannot see union contributions at all**.

Consequence if unfixed here: InterPro would be under-reported in
`contributing_sources`, and curated UniProt GO would be indistinguishable from
domain-**inferred** InterPro GO — a distinction that matters precisely because
these are inferred.

### 4.2 Design options

| | shape | assessment |
|---|---|---|
| **A** parallel list | `go_terms: [...]` + `go_terms_sources: [...]`, index-aligned | fragile — any consumer that filters or reorders one array breaks alignment |
| **B** per-token map | `go_terms_source: {"GO:0003677": ["uniprot","interpro"]}` | robust, self-describing, survives filtering; the only option that represents multi-source tokens honestly |
| **C** separate fields | `interpro_go_terms` alongside `go_terms` | zero risk to existing consumers, structural provenance; N fields per source, adapters must union |

**Recommendation: B as the end state, C as the shippable first slice.** They
compose — ship C, migrate to B when the union work lands. B matters because
multi-source tokens are the *common* case here, not the exception: 146,576 of the
InterPro GO pairs are already present from another source.

### 4.3 Scope of the union `track_source` work

- Implement per-token source capture in `_resolve_union`.
- Rewrite `_compute_contributing_sources` / `_has_source_label` to read the new
  structure (they currently assume scalar `*_source` fields).
- Decide compat for the 3 existing scalar `*_source` fields: keep as-is, or
  migrate to the same map shape.
- Decide whether **inferred vs curated** is a first-class attribute or implied by
  source identity.
- Size the merged-JSON delta.
- Blast radius: every union consumer — adapters, `annotation_quality`, the
  8-bucket source list, snapshot fixtures, KG-validity tests.

## 5. Ordering — lazy reference load inside step 2

The dependency is on the **global reference**, not per-strain calls.json (step 2
already requires the latter). So only the reference needs to be available.

**Precedent:** `load_pfam_data(cache_root)` is called directly inside
`build_gene_annotations.main()` — the Pfam reference is **not a numbered step**;
the merge fetches and caches its own global reference on demand.

**Decision: mirror it.** Step 2 calls `build_interpro_reference.build()`, which is
already idempotent (returns the cached dict unless `--force`). Step 9 remains the
explicit refresh / CI entry point. No renumbering, no `prepare_data.sh` churn,
and step 9 keeps its independence from steps 0–8.

Rejected: renumbering so the reference precedes step 2 (churns the script, all
docs, and muscle memory); adapter-only enrichment (bypasses
`annotation_quality` / `informative_annotation_types` and creates a second path
emitting into the same `Gene_involved_in_*` edge types).

## 6. Data flow / work items

1. **Extend step 9** to also parse EC + CAZy db_xrefs from `interpro.xml.gz`
   (GO + pathways already done). Same streaming parser, additional `include_dbs`.
   Emit sparse `ec_numbers` / `cazy_ids` keys per entry, alongside `go_terms`.
   - Store the entry `type` alongside so §3 filters are applicable downstream
     without a second lookup (already present as `type`).
2. **Implement union `track_source`** (§4.3) — blocking.
3. **Extend the `interproscan` source** in `gene_annotations_config.yaml` with the
   §2 field contributions + §3 filters. Filters are entry-type-dependent, so they
   need expression in the config (new `entry_type_filter`-style key) or a
   dedicated post-merge enrichment function.
   - **Recommendation: post-merge enrichment**, mirroring `enrich_pfam_fields`.
     The filter logic is type-aware and per-field; encoding it in the declarative
     config would strain the schema.
4. **Add an `hmm_PGAP.tsv` reference fetch** (§9.1d) — ~18 MB, same lazy-load
   pattern as §5, cached + committed like `interpro_reference.json`. Feeds
   `go_terms` / `ec_numbers` / descriptions from NCBIfam family records.
   Investigate the 27,361 unmatched NCBIfam hits (retired families?).
5. **Lazy-load the references in step 2** (§5).
6. **Rebuild** step 2 across 42 strains; rebuild KG; regenerate snapshot.

## 7. Validation gates

- Unit tests for the new parsers + filters (fixtures, not live downloads).
- `pytest -m "not slow and not kg"` clean — baseline is **2,127 passed**.
  *(Note: a fresh worktree without `.env` shows 28 failures / 43 errors in the
  omics adapter tests — unrelated, caused by a missing `OPENAI_API_KEY`.)*
- `/omics-edge-snapshot` before/after — expression edges must be unchanged.
- `pytest -m kg` after rebuild.
- Per-field assertion that measured deltas match §2 within tolerance; a large
  divergence means a filter regressed.
- Assert `contributing_sources` now includes `interproscan` for genes with
  InterPro-derived tokens (this is the observable proof §4 landed).

## 8. Deferred

### 8.1 Residue-level sites + unintegrated signatures — **separate spec**

Held off per user decision (2026-08-07). Recorded here because the source data is
**perishable**.

Measured by re-parsing the raw batch output (42 strains, 891,406 matches):

| | scale |
|---|---:|
| residue-level site records | **250,609** on **45,200 proteins (36% of proteome)** |
| ↳ on *unintegrated* matches (no IPR entry → currently dropped twice) | 184,281 (74%) |
| discontinuous domain locations (collapsed to start/end today) | 26,952 |

Site descriptions are exactly mutation-relevant: `active site` 11,886,
`NAD.` 7,326, `dimer interface` 6,702, `Proton acceptor.` 6,266,
`Substrate.` 5,303, `ATP binding site` 4,871, `Mg(2+)` 4,719 — with exact
residues (e.g. CDD `cd02204` on `WP_011131641.1`: ATP binding site at
E506, N510, G559, G560, N561).

Two reasons this needs its own spec rather than folding in here:
- **It does not fit the ontology pattern.** A site is positional and belongs to a
  (gene × signature) pair, not to a term. Needs a modelling decision — edge
  properties vs. a signature-level node type.
- **74% of it hangs off unintegrated matches**, for which no
  `Gene_has_interpro_entry` edge exists (Phase-2 §1 dropped them by design). That
  design decision would have to be revisited.

> **⚠️ Perishability.** The `*.interproscan.raw.json` files — now **17 GB** across
> 42 strains after the `--goterms --pathways` re-run — exist only in the working
> checkout. They are gitignored by design (see `.gitignore`), so a fresh clone has
> zero. Extraction is a cheap **re-parse** while they survive; once deleted,
> recovery costs the full ~27h re-scan. Extract before pruning the cache, even if
> modelling waits. This applies equally to §8.3.
>
> Measured on the current (post-re-run) artifacts: MED4 alone carries 3,788 site
> records. The 250,609 / 45,200-protein figures above are from the pre-re-run
> batch and are the right order of magnitude, but re-measure before building.

### 8.2 MetaCyc pathways

InterPro carries 79,683 MetaCyc xrefs, **100% `PWY-*` (pathways)**.

Not redundant with existing data, contrary to first assumption: MNX already
provides **20,242 `metacyc.reaction`** and **54,328 MetaCyc chem** xrefs — i.e.
reactions and compounds, *not* pathways. The two are complementary.

Deferred because there is no consumer and it needs a new node type — not because
it duplicates anything. Note the KG is closer to a MetaCyc layer than it appears:
MNX could supply reaction/compound identity if the pathway level is ever wired.

**There is no InterPro→KEGG mapping.** Verified against the full release:
`db="KEGG"` occurs **0** times. `KEGG` survives only as a legacy token in
`interpro.dtd`'s allowed-`db` list. KEGG pathways remain eggNOG-KO-derived
permanently.

### 8.3 PANTHER subfamily (match-level) GO — **now a re-parse, not a re-scan**

PANTHER `treeGrafter` attaches GO to the *subfamily* match (e.g. `PTHR30478:SF0`
BETA SLIDING CLAMP), which is more specific than the family-level entry. In the
100-protein smoke, **46 of 81** match-level GO terms were unreachable at entry
level.

**Status changed by `b8172a29`.** The re-run passed `--goterms`, so match-level
`goXRefs` are now present in the raw output (634 matches in MED4 alone). But
`parse_interproscan_json` still reads only `signature.entry.goXRefs` — match-level
`m["goXRefs"]` is dropped. So this is now **a parser change plus a re-parse of
existing raw JSON**, not a re-scan.

Directly serves *"this hypothetical is upregulated — what is it"*: a subfamily
call is a more specific functional statement than the family entry it rolls up to.

Same perishability caveat as §8.1 — it depends on the raw JSONs surviving.

## 9. Open questions

### 9.1 NCBIfam — **scoped 2026-08-07; verdict: ship, via a different route**

Scoping resolved this. Summary first, evidence below:

- **NCBIfam should ship**, but sourced from NCBI's `hmm_PGAP.tsv` (family-level
  GO / EC / product_name), **not** as raw accessions into `ontology_terms`. It is
  the highest-precision annotation source measured anywhere in this spec.
- **The TigrRole bridge is viable but capped**, and is now an *optional
  enhancement* rather than the gate on NCBIfam.

#### 9.1a Where `TigrRole` comes from

Cyanorak's `tIGR_Role` GFF attribute — **numeric TIGR role IDs** (`132`, `156`,
`125`) with `tIGR_Role_description` as `mainrole / subrole`. 1,765 populated genes
in MED4. `cyanorak` is `organism_restricted` → Pro/Syn only, confirming that
**Alteromonas, Shewanella, Pseudomonas, Ruegeria have no roles at all**.

#### 9.1b The authoritative role table is not readily available

- JCVI FTP (`ftp.jcvi.org/pub/data/TIGRFAMs/`) — **gone**.
- NCBI `hmm_PGAP.tsv` — 24 columns, **no role-link column**. NCBI dropped the
  TIGR role hierarchy when it absorbed TIGRFAMs.
- Wayback probe — inconclusive (rate-limited). Not chased further; sourcing
  `TIGRFAMS_ROLE_LINK` + `TIGR_ROLE_NAMES` from an archive remains an open lead.

#### 9.1c Empirical derivation works, but has a hard ceiling

Pro/Syn strains carry **both** signals, so accession→role can be learned and
transferred:

| | |
|---|---:|
| donor genes (Cyanorak role **and** TIGRFAM hit) | 13,932 across 21 strains |
| distinct TIGR accessions learned | 850 |
| **unambiguous (single role)** | **804 (94.6%)** |
| dominant role ≥80% | 24 |
| ambiguous | 22 |
| **role-less genes that would gain a role** | **13,854** |

94.6% single-role consistency confirms accession→role is a real functional
relation. Payoff lands exactly where the gap is: KT2440 757, HP15 702,
AltDE1 680, EZ55 661, BGP6 660, W3-18-1 657, BS11 650, HOT1A3 648 (14–28% of each
strain's genes).

**The ceiling:** of 2,172 distinct TIGR accessions on role-less genes, only
**798 (37%)** are learnable from Pro/Syn donors. **9,615 genes remain unmapped**,
needing 1,308 accessions that never appear in a donor — because they are
heterotroph-specific families absent from streamlined *Prochlorococcus*:
`TIGR02937` (ECF sigma factor, 242 genes), `TIGR02532` (236), `TIGR00229` (PAS
sensor, 205), `TIGR01782` (130).

Note the irony: those are environmental-sensing and signalling families —
precisely the ones most relevant to *"nutrient exchange between organisms"*. The
empirical bridge is systematically blind to the most interesting half. This is
why sourcing the real table is worth an archive hunt before building the
derivation.

#### 9.1d The better route — `hmm_PGAP.tsv` family-level GO/EC

`https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv` (~18 MB) indexes 68,695
family records keyed by both `#ncbi_accession` (`NF*`) and `source_identifier`
(`TIGR*`), carrying `go_terms`, `ec_numbers`, `product_name`, `gene_symbol`,
`family_type`, `taxonomic_range`.

Measured against actual hits (40,940 of 68,301 NCBIfam hits resolved to a family
record; 27,361 unmatched — likely retired/uncovered families, worth a look):

| via NCBIfam | net-new | corroborating | **ratio** | genes gaining first |
|---|---:|---:|---:|---:|
| GO | 14,270 | 85,678 | **0.17** | 1,652 |
| EC | 1,083 | 19,127 | **0.06** | 350 |

**An EC ratio of 0.06 is the cleanest signal in this entire spec** — better than
InterPro entry-level EC even after the FAMILY-only filter (ratio 1.0), and better
than Pfam (0.10). The reason is structural: NCBIfam families are curated at
*equivalog* level ("same function"), so their EC assignments are near-definitive
rather than inherited down a domain hierarchy.

Low yield (350 first-EC genes) but near-zero noise. This is a **precision**
source, and it also gives an authoritative `product_name` / `gene_symbol` —
making explicit the provenance that PGAP already baked into `Gene.product`.

**Revised recommendation:** add a `ncbifam` contribution feeding `go_terms`,
`ec_numbers` and `alternate_functional_descriptions` from `hmm_PGAP.tsv`. Treat
raw accessions into `ontology_terms` as still-unjustified (no consumer), and the
TigrRole bridge as a separate optional item pending the archive hunt.

**What it is:** NCBI's curated protein-family HMM collection, which **absorbed
TIGRFAMs** in 2021. Confirmed in the data: of 67,918 NCBIfam hits,
**40,834 are `TIGR*`** vs 27,084 `NF*` — 60% *are* TIGRFAMs.

**You already have this data and discard it.** MED4's Cyanorak `protein_domains`
column carries **3,242 `IPR*` and 709 `TIGR*` tokens**, and `enrich_pfam_fields`
drops both as "unresolved" (keeping only `PF*`). The route is also thin and
Pro/Syn-only — `cyanorak` is `organism_restricted`. InterProScan supplies 40,834
TIGR hits across **all 42 strains**.

**You also already have the evidence unlabelled.** NCBI PGAP uses NCBIfam to
assign product names, so `Gene.product` is *derived from* these models — the
accession is the missing provenance, not the information.

**The question that decides it:** `CyanorakRole` is Pro/Syn-restricted, so
**Alteromonas, Shewanella, Pseudomonas and Ruegeria have no functional roles at
all**. If TIGRFAM accessions can be mapped to TIGR roles (`TigrRole` nodes
already exist), NCBIfam stops being an inert accession dump and becomes **role
coverage for the heterotrophs** — the currently-dark half of "nutrient exchange
between organisms".

→ **Next action: scope the TIGRFAM→TigrRole mapping feasibility.** If it works,
`ontology_terms` (or a dedicated field) ships and likely earns a
`Gene_has_tigr_role` extension. If not, defer alongside MetaCyc.

### 9.2 `annotation_quality` semantics

7,029 genes gaining a first GO term shifts the `annotation_state` distribution.
Is domain-**inferred** GO "informative"? Defensible yes, but it is a semantic
call, and the 8-bucket source list is explicitly hand-maintained
(`docs/superpowers/specs/2026-05-01-explorer-frictions-resolution-design.md`).

### 9.3 Circularity with Pfam

eggNOG contributes Pfam; InterPro integrates Pfam signatures. Adding InterPro's
Pfam hits to `pfam_ids` is **corroboration, not independent evidence** — the 0.10
ratio reflects exactly that. Fine to ship, but any downstream confidence metric
that counts "number of independent sources agreeing" must not treat these two as
independent.
