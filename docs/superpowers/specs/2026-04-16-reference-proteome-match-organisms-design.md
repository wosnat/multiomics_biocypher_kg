# Reference Proteome Match Organisms

**Date**: 2026-04-16
**Status**: Approved
**Scope**: Moreno 2023 S3 (Alteromonas) + S4 (Marinobacter) community-fraction proteomics

## Problem

Some proteomics papers grow non-axenic cultures and search peptide spectra against a multi-organism reference database (e.g., MarRef v6 with ~4.5M protein entries). Each spectrum is assigned to the single best-matching reference protein. The resulting "Alteromonas proteome" or "Marinobacter proteome" is not a single cultured strain -- it is the union of best-match assignments across many reference entries, though in practice the best-match set funnels to one well-curated reference proteome per genus.

The current schema assigns each Experiment to a single OrganismTaxon (a strain), and every `Changes_expression_of` edge points at a Gene belonging to that strain. This works mechanically but misrepresents the data -- it implies a cultured strain when the organism identity was inferred by database matching.

See `docs/community_proteomics_marref_saga.md` for the full history.

## Two distinct community-data concepts (only one in scope)

1. **Reference proteome match** (in scope): experimental data searched against a multi-organism reference DB. Gene IDs are inherited from the best-matching reference proteome entry. The organism identity is inferred. Example: Moreno 2023 MarRef matching.

2. **MAG (Metagenome-Assembled Genome)** (out of scope, deferred): de novo assembly from metagenome/metatranscriptome reads produces novel gene IDs that don't exist in any reference genome. Example: Ziegler 2025 heterotroph metaflye contigs. This requires a different approach (Diamond-based mapping, new gene nodes) and is deferred.

## Design

### New organism_type value

Add a new `organism_type` property to OrganismTaxon nodes. Values:

| Value | Meaning | Examples |
|---|---|---|
| `genome_strain` | Default. Real genome assembly from a cultured isolate. | MED4, MIT9313, EZ55 |
| `treatment` | Non-genomic organism used as an edge source (coculture partner, phage). No genes. | Phage, Vibrio parahaemolyticus |
| `reference_proteome_match` | Organism identified by matching experimental data against a multi-organism reference database. Genes are the reference proteome's genes. | Marinobacter (MarRef v6), Alteromonas (MarRef v6) |

The property is set by the CyanorakNcbi adapter based on the source data:

- **Genome organisms** (`cyanobacteria_genomes.csv`): new `organism_type` column. Default: `genome_strain`. Reference-match organisms set to `reference_proteome_match`.
- **Treatment organisms** (`treatment_organisms.csv`): adapter infers `organism_type: treatment` from the source file.

All OrganismTaxon nodes get this property explicitly -- no organism should have `organism_type` absent.

### Concrete changes for Moreno 2023

Two existing OrganismTaxon entries are converted:

| Current entry | New preferred_name | Assembly | Locus prefix | organism_type |
|---|---|---|---|---|
| `Marinobacter adhaerens HP15` (GCF_000166295.1) | `Marinobacter (MarRef v6)` | GCF_000166295.1 (unchanged) | HP15_NNNN | reference_proteome_match |
| `AltMedDE` (GCF_000020585.3, **wrong assembly**) | `Alteromonas (MarRef v6)` | GCA_003513035.1 (**fixed**) | DEH24_NNNNN | reference_proteome_match |

Both entries were added to the KG solely to support Moreno 2023 community-fraction data. Neither has any other data in the corpus.

### New properties on OrganismTaxon (reference_proteome_match only)

- `organism_type`: `"reference_proteome_match"`
- `reference_database`: `"MarRef v6"` (the database used for matching)
- `reference_proteome`: UniProt proteome accession or assembly accession (e.g., `"UP000262181"` for the Alteromonas entry, `"GCF_000166295.1"` for the Marinobacter entry)

These properties are null/absent on genome_strain and treatment organisms.

### Data source: cyanobacteria_genomes.csv

Three new columns added to the genomes CSV:

| Column | Default | Example (HP15 row) | Example (MED4 row) |
|---|---|---|---|
| `organism_type` | `genome_strain` | `reference_proteome_match` | `genome_strain` |
| `reference_database` | (empty) | `MarRef v6` | (empty) |
| `reference_proteome` | (empty) | `GCF_000166295.1` | (empty) |

The CyanorakNcbi adapter reads these columns and passes them through as OrganismTaxon node properties. For treatment organisms, `organism_type` is set to `"treatment"` by the adapter based on the source file (`treatment_organisms.csv`), with no additional columns needed there.

### Gene nodes

Reference proteome match organisms own their genes through the normal `Gene_belongs_to_organism` edges. Gene nodes go through the standard genome pipeline:

1. Download genome data (GFF, protein FASTA, GBFF)
2. Build gene_mapping.csv
3. Build gene_annotations_merged.json (eggNOG, UniProt, Cyanorak where applicable)
4. Build gene_id_mapping.json
5. OG membership, GO/KEGG/Pfam annotations -- all normal

Expression edges (`Changes_expression_of`) from Moreno experiments point at these genes.

### Gene IDs

Gene IDs use the standard `ncbigene:<locus_tag>` prefix, same as genome strains.

**Known compromise**: `ncbigene:` is semantically imprecise for reference-match organisms -- these genes are identified by proteome database matching, not authoritative NCBI gene records. However, changing the prefix for a subset of organisms would require touching every adapter that constructs gene IDs (CyanorakNcbi, omics_adapter, functional_annotation_adapter, ortholog_group_adapter, etc.). The cost is not justified for 2 organisms today.

**Future revisit trigger**: if the corpus grows more reference-match organisms, or an ID collision arises (e.g., someone publishes pure-culture HP15 RNA-seq and we need "real" HP15 Gene nodes alongside the MarRef-matched ones), introduce a `make_gene_id(locus_tag, organism_type)` utility and refactor all adapters to use it.

### Organism IDs

Follow existing conventions:

- Marinobacter: `insdc.gcf:GCF_000166295.1` (RefSeq assembly, unchanged)
- Alteromonas: `insdc.gcf:GCA_003513035.1` (GenBank assembly; no RefSeq version exists)

### Organism naming convention

`preferred_name` format: `"<Genus> (MarRef v6)"` -- genus-level, referencing the source database and version.

Names are paper-agnostic: if a future paper also matches Marinobacter peptides to HP15 via MarRef v6, it reuses `"Marinobacter (MarRef v6)"`. A different reference DB version or different best-match reference entry would get a new organism node.

### Reusability rules

- Same reference DB + same best-match reference proteome = shared organism node across papers
- Different reference DB version (e.g., MarRef v7) = new organism node
- Different best-match entry within the same DB = new organism node
- Same genus matched via a completely different DB (e.g., MGnify) = new organism node

### Schema and edge types

No new node types or edge types. All existing infrastructure works:

- `Gene_belongs_to_organism`: reference-match organism owns its genes
- `Changes_expression_of`: Experiment -> Gene (unchanged)
- `Has_experiment`: Publication -> Experiment (unchanged)
- `Gene_in_ortholog_group`: OG membership enables cross-organism lookup
- Post-import computed properties: gene_count, expression stats, etc. all work normally

The `organism_type` property is the only schema addition.

### What stays the same

- All existing adapters: no code changes needed for gene ID construction
- Experiment nodes: standard properties, no new flags
- Expression edges: unchanged format and semantics
- Post-import scripts: no changes
- MCP tools: queries work as-is; agents can use `organism_type` for provenance-aware answers

## Scope boundaries

- **In scope**: Moreno S3 + S4 reference-match organisms; `organism_type` property; assembly fix for Alteromonas
- **Out of scope**: MAG/de novo assembly organisms (Ziegler heterotrophs) -- different concept, deferred
- **Out of scope**: Gene ID prefix refactor (`ncbigene:` -> `refmatch:`) -- deferred until warranted
- **Out of scope**: Ziegler 2025 MED4 ID translation (solved problem, separate task)
- **Out of scope**: Synthetic community treatment representation (Ziegler `treatment_organism: "Synthetic heterotroph community"`)

## Downstream impact (what-changed doc)

After implementation, create `docs/kg-changes/reference-proteome-match-organisms.md` covering:

- **New `organism_type` property on all OrganismTaxon nodes**: `"genome_strain"` (25 organisms), `"treatment"` (5 organisms), `"reference_proteome_match"` (2 organisms). 32 total. Queryable for filtering/provenance.
- **New properties on reference_proteome_match organisms**: `reference_database`, `reference_proteome`. Absent on other organism types.
- **Renamed organisms**: `Marinobacter adhaerens HP15` → `Marinobacter (MarRef v6)`, `AltMedDE` → `Alteromonas (MarRef v6)`.
- **Fixed assembly**: Alteromonas entry changed from GCF_000020585.3 (wrong) to GCA_003513035.1 (correct MarRef-matched reference). Gene count and annotations will differ.
- **MCP / query impact**: agents filtering by organism name must use updated names. `organism_type` property available for provenance-aware queries (e.g., exclude community-fraction data when looking for pure-culture experiments).
- **No changes to**: edge types, gene ID format, experiment properties, expression edge properties, post-import scripts.

## Relationship to existing documentation

- `docs/community_proteomics_marref_saga.md`: update "Open work" section to reference this spec; mark questions 1-3 as answered
- `CLAUDE.md`: update organism counts and strain lists after implementation
- `cyanobacteria_genomes.csv`: modify HP15 and AltMedDE entries (preferred_name, assembly for AltMedDE, new columns)
