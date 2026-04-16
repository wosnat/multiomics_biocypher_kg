# Community-fraction proteomics via MarRef — open issue

**Status**: open / deferred. Affects Moreno 2023 S3 + S4 today; will recur for any non-axenic community proteomics paper that uses a multi-organism reference DB (MarRef, MGnify, etc.).

**TL;DR**: Our current schema assigns each `Experiment` to a single `OrganismTaxon` (a strain), and every `Changes_expression_of` edge points at a `Gene` belonging to that strain. This breaks for proteomics studies that grow non-axenic cultures and assign peptides to whichever curated reference proteome best matches each spectrum — the "organism" in the data is not a cultured strain but a per-peptide reference-best-match. Until we design a schema-level solution, we proxy by treating the dominant matched reference proteome as the strain, with documented caveats on every affected experiment.

## Why it matters

Some proteomics papers in our corpus are non-axenic by design: the authors deliberately let environmental heterotrophs co-occur with their cyanobacterial cultures to mimic ocean conditions. The peptide spectra are searched against an **all-marine-organisms** reference database (e.g. MarRef v6 has ~4.5 million protein entries). Each spectrum is assigned to the single best-matching protein. The aggregated "Alteromonas proteome" or "Marinobacter proteome" reported in the paper is therefore:

- **Not a single cultured strain** — it's whichever environmental contaminants happened to grow.
- **Not a single reference proteome** — it's the union of best-match assignments across many MarRef entries, but in practice the best-match set funnels disproportionately to one or two well-curated reference proteomes per genus.
- **Identifier-coherent at the proteome level** — within one paper, all Alteromonas proteins typically end up labeled with the same locus-tag prefix because they all best-match the same MarRef reference.

Our schema's `Changes_expression_of (Experiment) → (Gene)` edges require a real gene node. A real gene node requires a real assembly. An assembly is a single strain. So the data has to be pinned to *some* strain or it can't be ingested.

## The Moreno 2023 case (concrete example)

Paper: Moreno-Cabezuelo et al., 2023, "Glucose Metabolism in Marine Picocyanobacteria", *Microbiology Spectrum*.

Methods (verbatim): *"to identify and categorize as many proteins as possible in the samples, we searched our data against an extensive database containing 4,565,004 entries from marine organisms (MarRef v6 database)"*.

Cultures were **non-axenic** (explicitly designed that way). Three perspectives were reported:

| Supp table | "Organism" in CSVs | Reality |
|---|---|---|
| S2 | One of the 5 cyano strains the lab grew (MED4, SS120, WH7803, WH8102, BL107) | Real cultured strain; clean schema fit |
| S3 | "Alteromonas sp." (genus taxid 232), GN=`DEH24_NNNNN` | MarRef best-matches to GenBank assembly **GCA_003513035.1** (Univ. of Queensland Alteromonas sp. scaffold WGS, UniProt proteome **UP000262181** — currently marked **Excluded — Redundant proteome** by UniProt) |
| S4 | "Marinobacter adhaerens (strain DSM 23420 / HP15)" (taxid 225937), GN=`HP15_NNNN` | MarRef best-matches to **GCF_000166295.1** (NCBI complete genome of HP15) |

The authors did not culture the heterotrophs — they reported what MarRef matched.

## What we tried for Alteromonas

1. **Initial guess** (wrong): I picked *Alteromonas mediterranea DE* (taxid 1774373, RefSeq `GCF_000020585.3`, locus prefix `MADE_RS`) on the basis of an unsourced hunch that `DEH24_` belonged to A. mediterranea DE. The paperconfig was written to use this strain.
2. **NCBI verification** (caught the error halfway): a verification agent confirmed that `DEH24_` is **not** an NCBI-registered locus prefix for `GCF_000020585.3` — that assembly uses `MADE_RS#####`. The agent flagged the discrepancy and suggested an `id_translation` bridge from `DEH24_*` → `MADE_RS_*` via Diamond protein matching (the same approach used for EZ55 and MIT1002).
3. **Bridge attempt** (suspended): I activated the S3 paperconfig blocks pointing at the canonical AltMedDE assembly with the bridge as deferred work (commit `895f1c9`). DEH24 IDs would have ~0% match rate until the bridge ran.
4. **Identifying the actual reference assembly** (correct): User pointed me at UniProt proteome `UP000262181` and NCBI assembly `GCA_003513035.1`. Confirmed via NCBI Datasets that this assembly's annotation uses `locus_tag=DEH24_NNNNN` for every CDS. This is what MarRef matched against.
5. **Realising the bigger problem**: The user pointed out that even with the right assembly, this isn't really "Alteromonas mediterranea DE" or any other single cultured strain — it's an environmental community fraction whose peptides happen to map best to one reference proteome. Tagging a single strain misrepresents the data.

## What we tried for Marinobacter

1. **HP15 chosen** because the S4 CSV `Description` column has explicit UniProt OS strings: `OS=Marinobacter adhaerens (strain DSM 23420 / HP15) OX=225937 GN=HP15_NNNN` for every protein. RefSeq `GCF_000166295.1` matches.
2. **Realised same caveat applies** — the HP15 strain assignment is the MarRef best-match, not what was in the cultures. HP15 is more directly evidenced than the Alteromonas case (every protein OS string says HP15) but the conceptual issue is identical.

## Decisions in the current PR (recorded for traceability)

- **Strain entries**: kept in `data/Prochlorococcus/genomes/cyanobacteria_genomes.csv`:
  - `AltMedDE` → `GCF_000020585.3` (still wrong assembly for the Moreno data — see "Open work" below).
  - `HP15` → `GCF_000166295.1` (correct reference proteome for Moreno's MarRef matches; community-fraction caveat applies).
- **Moreno S3 (Alteromonas)**: paperconfig blocks are present and active but match rates are low (DEH24 → MADE_RS bridge not built). 5 supp_table_s3 blocks in `data/Prochlorococcus/papers_and_supp/moreno 2023/paperconfig.yaml`.
- **Moreno S4 (Marinobacter)**: paperconfig blocks active, point at HP15. Match rates expected ~95%+.
- **EggNOG**: ran for HP15 + AltMedDE alongside SS120 + BL107 (added 2026-04-15).
- **Why we kept everything as-is**: this issue will recur for any future community-proteomics paper. Removing the strains and experiments now and re-adding them later is more disruptive than leaving them with documented caveats and fixing the schema once.

## Open work (out of scope for the current PR)

> **Update 2026-04-16:** Questions 1-3 below have been answered by the reference proteome match organisms spec (`docs/superpowers/specs/2026-04-16-reference-proteome-match-organisms-design.md`). Summary: HP15 and AltMedDE converted to `organism_type: reference_proteome_match` with `reference_database: MarRef v6`. AltMedDE assembly fixed to GCA_003513035.1, renamed to Alt_MarRef. See `docs/kg-changes/reference-proteome-match-organisms.md`.

The right fix needs a small schema design exercise. Key questions to answer:

1. **How should "community fraction matched against reference proteome X" be modeled?** Options to consider:
   - Genus-level `OrganismTaxon` nodes with an `Inferred_via_proteome` edge to the reference assembly's `OrganismTaxon`. DE edges point at the genus-level node and carry a `reference_proteome` property.
   - Keep the per-strain `OrganismTaxon` for the reference and add a flag on the experiment (`is_community_fraction: true`, `reference_proteome: GCA_003513035.1`).
   - Synthetic "community gene" nodes that aggregate per-locus_tag contributions across a reference proteome, decoupled from the strain `OrganismTaxon`.

2. **Should the reference assembly be "the strain we put in the KG", or should we also load the actual culture's intended target species/genus separately?** The latter is a 2-organism story per experiment.

3. **For Alteromonas DE specifically**: replace the current `GCF_000020585.3` (MADE_RS) entry with `GCA_003513035.1` (DEH24_) so Moreno S3 resolves natively, and rename the strain to reflect that it's a reference proteome, not a cultured strain. Likely worth doing as part of the schema work since the entry exists only to support Moreno.

4. **Generalize**: walk the rest of the corpus and identify other non-axenic / community proteomics papers that might be silently misrepresented by the current schema.

## Pointers

- Spec for the CSV-ready papers batch: [`docs/superpowers/specs/2026-04-14-csv-ready-papers-batch-design.md`](superpowers/specs/2026-04-14-csv-ready-papers-batch-design.md)
- Moreno paperconfig: [`data/Prochlorococcus/papers_and_supp/moreno 2023/paperconfig.yaml`](../data/Prochlorococcus/papers_and_supp/moreno%202023/paperconfig.yaml)
- Cyanobacteria genomes registry: [`data/Prochlorococcus/genomes/cyanobacteria_genomes.csv`](../data/Prochlorococcus/genomes/cyanobacteria_genomes.csv)
- Treatment organisms (genus-level entries): [`data/Prochlorococcus/treatment_organisms.csv`](../data/Prochlorococcus/treatment_organisms.csv)
- Related saga (single-strain reannotation, has different mechanics — not the same problem): EZ55 (`plans/ez55_deploy.md`) and MIT1002 (`plans/mit1002_deploy.md`)
