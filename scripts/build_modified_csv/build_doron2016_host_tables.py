#!/usr/bin/env python3
"""Extract host-gene DE tables from Doron 2016 GSE63690 processed files.

Each ``GSE63690_Processed_data_<host>.txt.gz`` mixes host probes with 232
``BSV9_gp*`` phage (Syn9) probes on the same array. Phage probes are
FILTERED OUT — the phage has no Gene nodes and every phage row would become
a dangling edge. Probe IDs are ``<locus_tag>_at`` (verified: every host
probe ends ``_at``; no duplicates after stripping), so the suffix strip
yields native locus tags directly (SynWH7803_* / SYNW* / Syncc8109_*).

Value columns are per-timepoint infected/control log-ratios (limma):
WH7803 in minutes (S0.C0 ... S240.C240), WH8102/WH8109 in hours
(SYN9.0.C0 / S0.C0 ... *4.C4). The single ``F``/``P.Value``/``adj.P.Val``
triplet is limma's OMNIBUS series-level test ("changed at some point"),
NOT a per-timepoint p-value — the paperconfig wires it onto the terminal
timepoint only (decision D7(a) in plans/geo_paperconfig_updates.md).

Outputs (committed, uncompressed):
    data/Synechococcus/papers_and_supp/Doron 2016/doron2016_<host>_host_de.csv
with an added ``locus_tag`` column and the original value columns unchanged.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Synechococcus/papers_and_supp/Doron 2016")

HOSTS = {
    # host key -> (source file suffix, probe id column)
    "wh7803": ("WH7803", "Probe.ID"),
    "wh8102": ("WH8102", "Probe_ID"),
    "wh8109": ("WH8109_Syn9", "Probe.ID"),
}


def main() -> None:
    for host, (suffix, id_col) in HOSTS.items():
        src = PAPER_DIR / f"GSE63690_Processed_data_{suffix}.txt.gz"
        df = pd.read_csv(src, sep="\t")

        phage = df[id_col].str.startswith("BSV9_")
        host_df = df[~phage].copy()

        stripped = host_df[id_col].str.replace(r"_at$", "", regex=True)
        assert stripped.str.len().lt(host_df[id_col].str.len()).all(), "non-_at probe"
        assert not stripped.duplicated().any(), "duplicate locus tags after strip"
        host_df.insert(0, "locus_tag", stripped)

        out = PAPER_DIR / f"doron2016_{host}_host_de.csv"
        host_df.to_csv(out, index=False)
        print(
            f"{host}: {len(df)} rows -> {len(host_df)} host rows "
            f"({int(phage.sum())} BSV9 phage probes dropped) -> {out.name}"
        )


if __name__ == "__main__":
    main()
