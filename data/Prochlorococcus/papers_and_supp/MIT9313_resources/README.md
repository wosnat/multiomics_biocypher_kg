# MIT9313 Resources (not a paper)

**Citation:** n/a — this directory holds strain-level shared gene ID resources for *Prochlorococcus* MIT9313.
**DOI:** n/a
**Organism(s):** *Prochlorococcus* MIT9313
**Topic:** Registers a GenBank-derived locus-tag bridge table (`locus_tag` / `PMTid` / `P9313name` / `uniprot_id`) so that `build_gene_id_mapping.py` (prepare_data step 3) can resolve MIT9313 gene IDs used across multiple publications. No publication block; `omics_adapter` ignores it. Listed in `paperconfig_files.txt` so the mapping builder picks it up.

## Available data inventory

| File | Type | Content | KG status | Recommended action |
|------|------|---------|-----------|--------------------|
| `MIT9313_genbank.tsv` | TSV | GenBank-derived ID bridge for MIT9313: `locus_tag`, `PMTid` (old PMT-style), `P9313name` (P9313-style), `uniprot_id`, `product` | already in | — |
| `paperconfig.yaml` | YAML | Single `supplementary_materials` entry of type `id_translation` that maps the four ID columns onto the MIT9313 gene catalogue | already in | — |
| `paperconfig_orig.yaml` | YAML | Identical pre-resolution copy of paperconfig.yaml (no `_resolved` output because this is not a DE csv) | reference | — |

## Current paperconfig summary

- Experiments defined: 0 (no publication block)
- Statistical analyses (DE edges): 0
- Supplementary materials entry types: `id_translation` (1 entry: `mit9313_id_translation`)
- Organisms covered: *Prochlorococcus* MIT9313
- Table scope(s): n/a
- Non-DE evidence: none
- ID resolution: contributes `locus_tag` (NCBI RS-style), `PMTid`/`P9313name` as `old_locus_tag`, `uniprot_id` as `uniprot_accession` to the v2 three-tier gene ID mapping for MIT9313

## Recommended actions

1. **No action** — directory is already fully integrated as an ID-translation bridge; keep as-is.
2. **Skip** — no further data to add; the TSV is a pure identifier table, not per-gene evidence.

## Notes

- This is a *resources* directory, not a paper. Paperconfig deliberately omits a `publication:` block so `omics_adapter` treats it as inert.
- The bridge is critical for resolving MIT9313 gene IDs referenced by Thompson 2011 and Tolonen 2006 (both use `PMT0107`/`PMT_####` old locus tag style).
- See `docs/methods_gene_id_mapping.md` for how these bridge files feed the v2 mapping build.
