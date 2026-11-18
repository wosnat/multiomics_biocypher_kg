# OrganismTaxon assembly quality — design (short)

**Date:** 2026-11-18 · **Size:** small change · **Status:** approved 2026-11-18 (revised after explorer review — see "Resolution" at the end)
**Backlog item:** `plans/backlog.md` → KG semantics → "Assembly level + CheckM completeness on OrganismTaxon"

## Intent

A comparative-genomics answer like "gene X is absent in strain Y" means something only if Y's
assembly is complete. Today `OrganismTaxon` has no assembly-quality field, so every such answer
needs an outside lookup (the spock_osnat `2026-11-04-aa_utilization_genomics` analysis had to
look assembly quality up by hand). This change adds NCBI's own assembly metadata to every
genome-backed organism node, so a consumer (the explorer, an LLM agent) can weigh an absence
call, or filter to complete genomes, in the same query.

## Goals / non-goals

- Goal: each genome-backed `OrganismTaxon` (116 `genome_strain` + 2 `reference_proteome_match`)
  carries assembly level, RefSeq status, contig stats and NCBI's CheckM result.
- Goal: values come straight from NCBI. The KG does not compute or re-run CheckM.
- Non-goal: `treatment` organisms (5). They have no assembly, so the properties are absent.
- Non-goal: switching any strain to a newer assembly version (see CB0101 below), and any
  per-gene "absence confidence" score. Consumers combine the fields themselves.
- Anti-goal: presenting CheckM completeness as the gene-absence signal. It under-scores
  streamlined genomes (MED4 is a closed genome at 99.45; NS01 is a closed genome at 87.31).
  `assembly_level` is the primary signal, and completeness grades the drafts.

## What NCBI has (measured 2026-11-18, Datasets v2 `genome/accession/{acc}/dataset_report`)

The data is not in our GFF or GBFF. The Datasets report is the source, from the same API that
step 0 already downloads genomes from.

| | count (of 118 registry accessions) |
|---|---|
| report returned | 118 (117 directly; CB0101 only with `filters.assembly_version=all_assemblies`, see below) |
| `assembly_level` | Complete Genome 71 · Contig 30 · Scaffold 11 · Chromosome 5 (CB0101 not counted) |
| `assembly_status` | current 107 · **suppressed 10** (CB0101 = previous) |
| CheckM `completeness` | 112. Missing on WH7803, HNLC1, HNLC2, RCC307 (all suppressed) and Alt_MarRef |
| CheckM `contamination` | 73 have a value, 39 omit it (see finding 3) |
| completeness, lowest | NS01 87.31 (Complete) · AltDE 96.52 (Complete) · PCC7001 97.03 (Scaffold) · median 99.53 |

### Three findings that shape the design

1. **10 of our assemblies are suppressed by RefSeq.** BS11, WH7803, HNLC1, HNLC2, A18-46.1,
   WH8103, A18-40, RCC307, RS9915 and BOUM118 are "removed by RefSeq staff". The reason is
   "unverified source organism" for 6 of them and "missing strain identifier" for WH7803 and
   HNLC2. The genomes and genes are still valid KG content: papers use WH7803 and BS11, and
   their annotation was downloaded before the suppression. But a consumer should be able to
   see the status, so the design adds `assembly_status` and the suppression reason.
2. **CB0101 is pinned to a superseded version.** The registry has `GCF_000179235.1`
   (contig-level, 2010). NCBI's current record is `GCF_000179235.2` (Complete Genome, 2019,
   same sequence per `revision_history`). A plain report request returns nothing for `.1`,
   so the fetch passes `filters.assembly_version=all_assemblies`. CB0101 then reads
   `assembly_level = contig`, `assembly_status = previous`, and
   `assembly_current_accession = GCF_000179235.2` so the node itself says a newer version
   exists. Moving the strain to `.2` is a separate strain-migration task, added to the
   backlog. It may also place some of CB0101's 17 unplaced Cyanorak genes.
3. **A missing `contamination` value means 0.** The 39 genomes without one are closed or
   clean genomes (MED4, MIT9312, AS9601, PCC7002, WH8102, …). The smallest value NCBI reports
   is 0.02, never 0.00, and proto3 JSON leaves out fields whose value is zero.
   **Design: when the report has a CheckM block but no `contamination`, store `0.0`.** When
   there is no CheckM block, both properties are absent. (The backlog said "absent ≠ 0".
   That was wrong for this source.)

## Approaches

- **A — committed per-strain report JSON, read by the organism adapter (chosen).** Step 0
  fetches each strain's report into `cache/data/<org>/genomes/<strain>/assembly_report.json`
  (committed, like the strain's `genomic.gff`). `cyanorak_ncbi_adapter._get_organism_node`
  reads it. The Docker build in the other clone needs no network access, and the values are
  reviewable in git diffs.
- B — fetch at KG build time, like `_fetch_ncbi_taxonomy`. Rejected: it adds a network
  dependency to the build, and values change without anyone seeing a diff.
- C — hand-maintained columns in `cyanobacteria_genomes.csv`. Rejected: 118 rows × 10 fields
  of upstream data typed by hand, and they go stale.

## Design

### New `OrganismTaxon` properties

All are sparse: absent on `treatment` organisms and wherever NCBI gives no value. The one
exception is `checkm_contamination` (finding 3).

| property | type | from report | notes |
|---|---|---|---|
| `assembly_level` | str | `assembly_info.assembly_level` | controlled vocab, lower snake case (R1): `complete_genome` \| `chromosome` \| `scaffold` \| `contig` |
| `assembly_status` | str | `assembly_info.assembly_status` | closed controlled vocab: `current` \| `suppressed` \| `previous`. **An unknown value fails the build** (NCBI's enum has more states than we observe) |
| `assembly_current_accession` | str | `current_accession` | only when it differs from the node's accession (today: CB0101 → `GCF_000179235.2`) |
| `assembly_report_date` | str | file's `fetched_at` | ISO date the report was fetched, so "suppressed as of when?" has an answer |
| `assembly_suppression_reason` | str | `assembly_info.suppression_reason` | free text, only when suppressed; the "contact info@…" boilerplate is stripped, leaving e.g. `unverified source organism` |
| `assembly_name` | str | `assembly_info.assembly_name` | e.g. `ASM1146v1` |
| `assembly_release_date` | str | `assembly_info.release_date` | ISO date |
| `genome_length` | int | `assembly_stats.total_sequence_length` | bp; NCBI sends a string, stored as int |
| `contig_count` | int | `assembly_stats.number_of_contigs` | |
| `contig_n50` | int | `assembly_stats.contig_n50` | bp |
| `gc_percent` | float | `assembly_stats.gc_percent` | NCBI's rounded value |
| `checkm_completeness` | float | `checkm_info.completeness` | 0–100; always cast to float (NCBI sends `100` as a JSON int, e.g. MIT9312) |
| `checkm_contamination` | float | `checkm_info.contamination` | 0–100. `0.0` = NCBI omitted the field (zero; finding 3) — the property description says so |
| `checkm_marker_set` | str | `checkm_info.checkm_marker_set` | verbatim taxon name (an external term, R1), e.g. `Prochlorococcus marinus` |
| `checkm_marker_set_rank` | str | `checkm_info.checkm_marker_set_rank` | verbatim, e.g. `species` \| `genus` \| `order` \| `domain`. A coarse marker set (e.g. `Bacteria`) gives a weaker completeness estimate |
| `checkm_version` | str | `checkm_info.checkm_version` | e.g. `v1.2.4` |

`ncbi_taxon_id` and the GCF accession (node id) already exist, so they are not duplicated.
`completeness_percentile` is left out: it is a percentile within NCBI's species set, absent on
most records, and not comparable across our genera.

**Reference-proteome-match organisms** (Marinobacter / Alteromonas MarRef v6): the values
describe the *reference* assembly the proteome match used, not the environmental sample. This
caveat goes in the **property** descriptions of the `assembly_*` / `checkm_*` entries, not in a
value description (it would otherwise read as applying to the 116 genome strains); consumers
tell the cases apart by `organism_type`.

### Fetch (step 0, `download_genome_data.py`)

- New `_ncbi_fetch_assembly_reports(genomes, force)` as part of sub-step 1. It runs for every
  strain, **including those whose genome files are already cached**, because those files skip
  the genome download today. It sends a GET to
  `…/genome/accession/{acc,acc,…}/dataset_report?filters.assembly_version=all_assemblies`
  in chunks of 40 (one 118-accession POST returned HTTP 500 on 2026-11-18) and writes each
  report to `<strain dir>/assembly_report.json` (indented, committed) as
  `{"fetched_at": "<ISO date>", "report": {…NCBI report verbatim…}}`. Without `--force` it
  only fetches strains that have no file yet, so `fetched_at` is what tells a reader how
  stale the RefSeq status may be.
- An accession with no report gets a warning and no file. The adapter then emits no
  assembly properties for it (absent, never `""`).

### Adapter, schema, vocabulary

- `cyanorak_ncbi_adapter`: a small `_assembly_properties(wrapped_report, accession) -> dict`
  helper (a pure function, unit-tested) merged into `_get_organism_node`; `_clean_str` is
  applied to string values; numeric fields are cast (`int` for the stats, `float` for
  CheckM and `gc_percent`). An `assembly_level` or `assembly_status` outside the declared
  vocabulary raises.
- `config/schema_config.yaml` → `organism taxon`: declare all 16 properties.
- `config/controlled_vocabularies.yaml`:
  - `OrganismTaxon.assembly_level` and `OrganismTaxon.assembly_status`: closed, sparse trust
    vocabularies with `value_descriptions` written for researchers. `complete_genome`:
    "the assembly is closed, so a missing gene is not an assembly gap — necessary for a
    confident absence call; still check the annotation". `suppressed`: "RefSeq removed the
    record (see `assembly_suppression_reason`); the genome and its genes remain valid KG
    content". `previous`: "a newer version exists (`assembly_current_accession`)".
  - Numeric entries, following `OrganismTaxon.peptidase_gene_count`: `checkm_completeness`
    and `checkm_contamination` (`float`, `min_value: 0`, `max_value: 100`, sparse),
    `genome_length`, `contig_count`, `contig_n50` (`int`, `min_value: 1`, sparse),
    `gc_percent` (`float`, 0–100, sparse).
  - `OrganismTaxon.checkm_marker_set_rank`: open (`closed: false`), sparse, values listed
    verbatim as NCBI reports them (R1 external term; observed 2026-11-18: `species`,
    `genus`, `family`, `order`, `domain`). The description says a coarser rank means a
    weaker completeness estimate. `checkm_marker_set` itself (a taxon name) gets no entry.
  - No entry for the free-text / identifier properties (`assembly_suppression_reason`,
    `assembly_name`, `assembly_release_date`, `assembly_current_accession`,
    `assembly_report_date`, `checkm_marker_set`, `checkm_version`) — the existing
    `OrganismTaxon` convention (only `organism_type` + the gene counts are declared).
  - This moves `Schema_info.controlled_vocabularies_hash`.

#### Vocabulary nodes (summary)

9 new `ControlledVocabulary` nodes, all `applies_to: OrganismTaxon` (153 → 162):

| kind | entries |
|---|---|
| closed categorical, with `value_descriptions` (trust vocabularies, +2) | `assembly_level`, `assembly_status` |
| open categorical | `checkm_marker_set_rank` |
| numeric float 0–100 | `checkm_completeness`, `checkm_contamination`, `gc_percent` |
| numeric int ≥ 1 | `genome_length`, `contig_count`, `contig_n50` |

The CLAUDE.md trust-vocabulary list and node count are updated with them.
- `config/gene_annotations_config.yaml`: add `assembly_quality` to the `ncbi` DataSource's
  `info_types`.
- `organismTaxonFullText` is unchanged: none of these fields is search text.

### Tests

- Unit: `_assembly_properties` on four fixture reports — MED4 (complete, contamination
  omitted → 0.0), MIT9312 (`completeness: 100` int → `100.0`), a suppressed draft with no
  CheckM block, and CB0101 (`previous`, `assembly_current_accession` set). Plus: an unknown
  `assembly_status` raises.
- Unit: the vocabulary loader accepts the new entries.
- KG (`tests/kg_validity/test_organism.py`): every `genome_strain` /
  `reference_proteome_match` node has `assembly_level`, `assembly_status` and
  `assembly_report_date`. No `treatment` node has any of them. Spot checks: MED4
  `complete_genome` / `current` / 99.45 / 0.0; WH7803 `suppressed`; CB0101 `contig` /
  `previous` / `GCF_000179235.2`.
- `test_schema_coverage.py` covers the declarations automatically.

### Docs

- CLAUDE.md: the OrganismTaxon key-graph-facts bullet and the Data Locations
  (`assembly_report.json`).
- CHANGELOG `[Unreleased]` → `### Added` + a `### Highlights` line ("organism nodes now say
  how complete their assembly is").
- Backlog: remove this item. Add "migrate CB0101 to GCF_000179235.2" and "review the 10
  RefSeq-suppressed assemblies".

## For the explorer to review

1. Do the property names fit how `list_organisms` / `kg_schema` present organisms? Should the
   `assembly_*` / `checkm_*` prefixes be kept, or would you rather have one nested-looking
   naming scheme?
2. Is `assembly_status` + `assembly_suppression_reason` enough for you to warn on suppressed
   genomes, or do you want a one-line caveat string?
3. Contamination: are you OK with storing `0.0` for the omitted value (finding 3)?
4. Which of these should `list_organisms` return by default? Our suggestion is
   `assembly_level`, `assembly_status` and `checkm_completeness`, with the rest available
   in the detail view. This is explorer-owned, so the KG side does not depend on the answer.

## Left to the builder

The chunk size and retry handling, the fixture trimming, and the exact wording of the
value descriptions.

## Explorer review (2026-11-18)

Overall sound. Approach A is right, the count tables reconcile (117 + CB0101 = 118;
completeness 112 = 118 − 5 − CB0101; 73 + 39 = 112), and none of the 14 property names
collide with existing `OrganismTaxon` keys in the live KG. Spot-checked against Datasets v2
for MED4 (`GCF_000011465.1`), MIT9312 (`GCF_000015965.1`), CC9311 (`GCF_000014585.1`) and
CB0101 (`GCF_000179235.1`).

### Fix before freezing

1. **Finding 3's percentile evidence is wrong.** MED4 omits `contamination` but has
   `completeness_percentile: 36.17`. CC9311 (`checkm_marker_set_rank: genus`) has
   `contamination: 0.12` but no percentile. The percentile's absence follows the marker-set
   rank (there is no species set to rank against), not the dropping of zero fields. Remove
   that sentence. The `0.0` conclusion still holds: the minimum reported value is 0.02, and
   proto3 JSON omits zero fields.
2. **Put the reference-proteome caveat at property level, not value level.** A caveat in the
   `complete_genome` value description would also apply to the 116 genome strains. Put it in
   the property `description` of the `assembly_*` / `checkm_*` entries. Consumers can tell
   the cases apart by `organism_type`.
3. **Soften "`complete_genome` is the gene-absence signal".** A closed assembly rules out
   assembly gaps, not annotation misses, and genes come from merged NCBI + Cyanorak calls.
   Suggested wording: "necessary for a confident absence call; still check annotation."
4. **`assembly_status` should fail loudly on unknown values.** NCBI's enum is larger than
   the three values observed (it has replaced / retired-style states). It is a closed
   vocabulary, so raise on an unknown value at build time rather than pass it through.
5. **Declare numeric vocabulary entries too.** Existing numeric organism properties (e.g.
   `OrganismTaxon.peptidase_gene_count`) have `ControlledVocabulary` entries with
   `min_value`. Add `checkm_completeness` / `checkm_contamination` (0–100) and the int
   stats for consistency.
6. **Cast `completeness` to float.** NCBI sends `100` as a JSON int (MIT9312).
7. **Record when the reports were fetched.** Without `--force`, files never refresh, and
   suppression status changes upstream. Add a `fetched_at` in each JSON file or a date on
   the `ncbi` DataSource `assembly_quality` entry, so "suppressed as of when?" has an answer.

### Answers to the four questions

1. **Naming:** keep it flat with the `assembly_*` / `checkm_*` prefixes. `list_organisms`
   returns flat columns, and the prefixes group well in `kg_schema`.
2. **Suppression:** `assembly_status` + `assembly_suppression_reason` is enough. The explorer
   writes its own warning, so the KG does not need a caveat string.
3. **Contamination `0.0`:** OK, provided the property description says "0.0 = NCBI omitted
   the field (zero)". Fix the evidence as in item 1.
4. **`list_organisms` defaults:** agreed. Compact view gets `assembly_level`,
   `assembly_status` and `checkm_completeness`; the rest go in verbose. The explorer will
   likely also add a `by_assembly_level` summary and possibly an `assembly_level` filter.
   That is explorer-scoped and does not block this change.

### Not blocking

- CB0101 reads `contig` / `previous` while `.2` is a complete genome, and nothing on the
  node hints at that. A short note in its property or vocabulary description would help
  until the migration happens.
- This change moves `controlled_vocabularies_hash`. The explorer's next absorb re-pins it
  and regenerates the goldens.

## Resolution (KG side, 2026-11-18)

All seven "fix before freezing" items accepted and folded into the design above:

1. Percentile sentence removed from finding 3; the `0.0` rule rests on the 0.02 minimum and
   proto3 zero-omission.
2. Reference-proteome caveat moved to the property descriptions.
3. `complete_genome` wording softened to "necessary for a confident absence call; still
   check the annotation".
4. `assembly_level` and `assembly_status` raise on unknown values at build time.
5. Numeric `ControlledVocabulary` entries added for the CheckM and stats properties.
6. `completeness` (and `contamination`, `gc_percent`) cast to float.
7. Each `assembly_report.json` wraps the report with `fetched_at`; surfaced on the node as
   `assembly_report_date`.

Non-blocking CB0101 note: handled with a sparse `assembly_current_accession` property
(KG call — data-driven from the report's `current_accession`, rather than a description note
that would go stale after the migration). The property count is now 16.
