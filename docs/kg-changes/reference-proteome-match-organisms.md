# Reference Proteome Match Organisms

**Date**: 2026-04-16
**Spec**: `docs/superpowers/specs/2026-04-16-reference-proteome-match-organisms-design.md`

## Changes

### New `organism_type` property on OrganismTaxon

All OrganismTaxon nodes now carry `organism_type`:
- `genome_strain` (25 organisms) — real genome assembly from a cultured isolate
- `treatment` (5 organisms) — non-genomic organisms used as coculture partners
- `reference_proteome_match` (2 organisms) — identified by matching experimental data against a multi-organism reference database

### New properties (reference_proteome_match only)

- `reference_database`: the database used for matching (e.g., "MarRef v6")
- `reference_proteome`: accession of the matched reference proteome

### Renamed organisms

| Old name | New name | organism_type |
|---|---|---|
| Marinobacter adhaerens DSM 23420 / HP15 | Marinobacter (MarRef v6) | reference_proteome_match |
| Alteromonas mediterranea DE | Alteromonas (MarRef v6) | reference_proteome_match |

### Assembly fix

Alteromonas entry changed from GCF_000020585.3 (wrong -- A. mediterranea DE RefSeq) to GCA_003513035.1 (correct -- the actual MarRef-matched reference proteome with DEH24_* locus tags). Strain name changed from `AltMedDE` to `Alt_MarRef`. Taxid changed from 1774373 to 232 (genus-level Alteromonas sp.).

### Query examples

```cypher
-- Filter to genome strains only (exclude community fractions)
MATCH (o:OrganismTaxon)
WHERE o.organism_type = 'genome_strain'
RETURN o.preferred_name

-- Find all reference proteome match organisms
MATCH (o:OrganismTaxon)
WHERE o.organism_type = 'reference_proteome_match'
RETURN o.preferred_name, o.reference_database, o.reference_proteome

-- Exclude community-fraction expression data
MATCH (e:Experiment)-[:Changes_expression_of]->(g:Gene)-[:Gene_belongs_to_organism]->(o:OrganismTaxon)
WHERE o.organism_type <> 'reference_proteome_match'
RETURN e.name, g.locus_tag
```

### No changes to

- Edge types or edge properties
- Gene ID format (`ncbigene:<locus_tag>`)
- Experiment node properties
- Post-import computed properties (gene_count, etc.)
