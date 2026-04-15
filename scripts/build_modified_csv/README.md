# `_modified.csv` builders

One script per omics paper that needs a derived view of its supplementary CSVs (sign-corrected log2FC, merged up/down halves, normalized column names, extracted locus-tag columns, etc.) before paperconfig ingestion.

## Convention

- One script per paper, named `build_<paperslug>_modified_csv.py`.
- Read the original supp CSV from `data/<Organism>/papers_and_supp/<paper>/`.
- **Never edit the original** — write `<original_name>_modified.csv` in the same directory.
- The paperconfig's `filename:` points at the `*_modified.csv`. Add a YAML comment near each `supp_table_*` block referencing this script for traceability.
- Re-running the script must be idempotent.

## Why a separate script (vs. inline)

The transformation is reproducible from a single script, so a fresh checkout of the repo can regenerate every `_modified.csv` deterministically. Keeps provenance: the original supp file matches what the paper deposited; the derivation lives in code.

## Existing builders

| Paper | Builder | Notes |
|---|---|---|
| Domínguez 2017 | `build_dominguez2017_modified_csv.py` | merges S3 up + S3 down → signed `log2_fold_change` |
| Fuszard 2012 | `build_fuszard2012_modified_csv.py` | adds `log2_fold_change = log2(Ratio of means)` per strain CSV |
| Moreno 2023 (cyano S2 + Marinobacter S4) | `build_moreno2023_modified_csv.py` | adds `log2FC.*` columns; per-strain override for SS120 `FC GmM/L` typo |
| Moreno 2023 (Alteromonas S3) | `build_moreno2023_s3_modified_csv.py` | extracts `DEH24_*` from `Description`, normalizes per-strain column quirks |

## Run

```bash
uv run python scripts/build_modified_csv/build_moreno2023_modified_csv.py
```

Pointers: [memory feedback note on the convention](../../../.claude/projects/-home-osnat-github-multiomics-biocypher-kg/memory/feedback_csv_modified_suffix.md).
