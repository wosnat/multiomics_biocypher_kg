# InterProScan per-strain domain annotation — design (Phase 1)

**Date:** 2026-07-22
**Skill:** `/interproscan-run`
**Status:** Phase 1 (per-strain inspectable artifacts; no KG coupling)

## Step 0 — intent

- **Predicts:** per-protein *list-of-matches*. InterProScan is an orchestrator
  that runs each protein against ~15 member databases (Pfam, NCBIfam/TIGRFAM,
  Hamap, PROSITE, SFLD, PANTHER, Gene3D, SUPERFAMILY, CDD, PRINTS, SMART,
  PIRSF, MobiDBLite, Coils) and maps every hit into a unified **InterPro entry**
  (`IPRxxxxxx`), attaching GO terms + pathway (Reactome/MetaCyc/KEGG) xrefs.
  One protein → many matches, each with member-DB signature accession, start/end
  coordinates, e-value/score, and (when integrated) an InterPro id + GO/pathway.
- **Input:** `protein.faa` (RefSeq WP_ accessions as FASTA ids).
- **Install flavor:** A (public Docker image `interpro/interproscan:5.78-109.0`)
  + a separately-downloaded data directory (~6.4 G compressed) under
  `~/tools/InterProScan/` mounted read-only at `/opt/interproscan/data`.
- **Triggers:** "run interproscan on X", "predict protein domains",
  "InterPro / domain annotation for all strains", "functional domains for the
  new strain".

## Why InterProScan, given we already have eggNOG/UniProt

- eggNOG gives Pfam/GO/KEGG at full coverage but *by orthology transfer*.
  InterProScan does **direct HMM/profile matching against curated models** — an
  independent method, so it is a **quality/confidence cross-check** on eggNOG,
  and adds per-domain **coordinates** + the unifying **InterPro entry id**.
- Genuinely new to the KG (not delivered by eggNOG/UniProt): **NCBIfam/TIGRFAM,
  Hamap, PROSITE, SFLD**, the InterPro-entry integration layer, and (with all
  apps) structural **Gene3D/SUPERFAMILY**.
- Biggest *quality* win is the **heterotrophs**, where UniProt coverage is thin
  and IPS matching is method-independent of that gap.

## Application-set policy

Default = **all default free apps** (omit `--applications`; IPS runs every
activated analysis; the licensed SignalP/TMHMM/Phobius are deactivated in the
free distribution). Rationale: Phase 1 produces inspectable artifacts only, the
redundancy-vs-eggNOG decision is a Phase-2 call made with real overlap numbers,
and the marginal cost is one-time compute (accepted) + already-committed disk.
`--applications APP1,APP2` overrides per-run (e.g. focused
`Pfam,NCBIFAM,Hamap,PROSITEPatterns,PROSITEProfiles,SFLD`).

## Output schema

`<strain>.interproscan.calls.json` — top-level dict keyed by WP_ accession.
Proteins that were processed but matched nothing are **included** with an empty
`matches` list (so absence of a key = "not processed", `match_count == 0` = "no
domain found"):

```json
{
  "WP_002805124.1": {
    "md5": "…",
    "match_count": 3,
    "interpro_entries": ["IPR000484", …],
    "go_terms": ["GO:0009772", …],
    "pathways": ["KEGG:00195", …],
    "libraries": ["PFAM", "NCBIFAM", …],
    "matches": [
      {
        "library": "PFAM",
        "signature_accession": "PF00124",
        "signature_description": "Photosynthetic reaction centre protein",
        "interpro_accession": "IPR000484",   // null if unintegrated
        "interpro_description": "…",           // null if unintegrated
        "interpro_type": "FAMILY",             // null if unintegrated
        "start": 1, "end": 300,
        "evalue": 1e-50,                       // null for pattern/profile hits
        "score": 120.0,                        // null when N/A
        "go_terms": ["GO:0009772"],
        "pathways": ["KEGG:00195"]
      }
    ]
  }
}
```

One entry per (match × location) so multi-region signatures keep their
coordinates. Sorted by `(start, evalue, signature_accession)`.

`<strain>.interproscan.skill_summary.json` — per-strain QC:

```json
{
  "strain": "MED4",
  "tool_version": "5.78-109.0",
  "image_digest": "sha256:…",
  "applications": "ALL_DEFAULT",
  "input_proteins": 1858,
  "calls_made": 1858,            // proteins that appear in IPS output
  "proteins_no_match": 16,       // match_count == 0
  "parse_failures": 0,
  "total_matches": 5321,
  "interpro_integrated_matches": 4800,
  "distribution": {"PFAM": 1600, "NCBIFAM": 900, …},  // matches per library
  "sentinel_rate": 0.0086,       // proteins_no_match / input_proteins
  "wallclock_s": 1234
}
```

## Phase 2 sketch (deferred, separate spec)

Two plausible surfaces, decided later with overlap data in hand:
1. **Merge into `gene_annotations_merged.json`** as an `interproscan` logical
   source — InterPro ids + coordinates + GO as fields, joined on `protein_id`.
2. **New KG ontology** — `InterProEntry` nodes + `Gene_has_interpro_entry`
   edges (carrying coordinates + e-value), an `Interpro_entry_is_a_interpro_entry`
   hierarchy, and `Interpro_entry_maps_to_go` cross-links. Route via
   `/integrate-a-tool`.

Phase-1 artifacts sit in the strain cache for inspection and are **not** wired
into `gene_annotations_merged.json` or any adapter yet.
