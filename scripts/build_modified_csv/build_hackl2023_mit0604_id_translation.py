#!/usr/bin/env python3
"""Generate the MIT0604 author-ID -> Cyanorak locus-tag bridge for Hackl 2023.

The GSE171435 DESeq2 tables use the authors' custom annotation ids
(``MIT0604_00001``), which are our Cyanorak ids minus the ``CK_Pro_``
prefix (``CK_Pro_MIT0604_00001``) — plan finding F4, measured 97.9%
tier-1. This script materializes that prefix relationship as a committed
``id_translation`` CSV: one row per ``CK_Pro_MIT0604_*`` key in the
strain's ``specific_lookup``, mapping the bare author id to the canonical
locus tag the key resolves to.

Deterministic — no diamond run, no network. Re-run after a MIT0604
mapping rebuild if the Cyanorak id set changes.

Output: data/Prochlorococcus/papers_and_supp/Hackl 2023/mit0604_id_translation.csv
    columns: mit0604_id, locus_tag
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

MAPPING = Path("cache/data/Prochlorococcus/genomes/MIT0604/gene_id_mapping.json")
OUT = Path("data/Prochlorococcus/papers_and_supp/Hackl 2023/mit0604_id_translation.csv")

PREFIX = "CK_Pro_"


def main() -> None:
    mapping = json.load(MAPPING.open())
    lookup = mapping["specific_lookup"]

    rows = sorted(
        (key[len(PREFIX):], locus_tag)
        for key, locus_tag in lookup.items()
        if key.startswith(PREFIX + "MIT0604_")
    )

    with OUT.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mit0604_id", "locus_tag"])
        w.writerows(rows)

    print(f"{len(rows)} bridge rows -> {OUT}")


if __name__ == "__main__":
    main()
