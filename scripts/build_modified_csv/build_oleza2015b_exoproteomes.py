#!/usr/bin/env python3
"""Convert Christie-Oleza et al. 2015 (Environ Microbiol 17:3781, "Functional
distinctness in the exoproteomes of marine Synechococcus") supplementary tables
to clean CSVs + per-strain protein FASTAs for KG integration.

Paper dir: data/Synechococcus/papers_and_supp/Oleza 2015b/

Inputs (never modified):
  emi12822-sup-0004-si.xls  Table S4 — exoproteome NSAF per strain (8 Synechococcus
                            strains), one row per homology cluster; each strain block
                            = ID / function / AA sequence / predicted secretion /
                            BLAST-P / replicate % NSAF / average. A protein with
                            average 0 is a homologue (BLASTp) of an exoprotein detected
                            in another strain that was NOT detected in this strain.
  emi12822-sup-0005-si.xls  Table S5 — WH8102 per-condition protein list (sequences)
  emi12822-sup-0006-si.xls  Table S6 — WH8102 comparative exoproteomics, 4 contrasts
                            (PatternLab TFold signed linear fold change + p-value)

Outputs:
  table_s4_<strain>_exoproteome_nsaf.csv  one per strain: strain, protein_id, function,
                                       avg_pct_nsaf (blank when not detected),
                                       detection_status, n_replicates, n_detected
  table_s6_<contrast>.csv              one per WH8102 contrast: protein_id, function,
                                       fold_change (paper, signed linear),
                                       log2_fold_change, p_value, section
  protein_seqs_<strain>.faa            per-strain FASTA of every protein accession the
                                       tables cite (source for the diamond_protein_match
                                       id_translation entries in the paperconfig)

Conventions:
  * Signed linear fold change (PatternLab TFold: +x = x-fold up in the treatment arm,
    -x = x-fold down; |FC| >= 1 always) -> log2FC = sign(FC) * log2(|FC|).
  * Detection status per strain from the biological-replicate NSAF values:
    detected = non-zero in every replicate, sporadic = non-zero in some,
    not_detected = zero in all (listed homologue only). BL107 was run twice per
    biological replicate (A1/A2, B1/B2 — the paper re-ran samples with 4x material);
    the two runs of a replicate are collapsed (replicate detected if either run is
    non-zero) so BL107 has 2 biological replicates like the other non-WH8102 strains.

Idempotent: rerun overwrites outputs.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd

PAPER_DIR = Path("data/Synechococcus/papers_and_supp/Oleza 2015b")
S4 = PAPER_DIR / "emi12822-sup-0004-si.xls"
S5 = PAPER_DIR / "emi12822-sup-0005-si.xls"
S6 = PAPER_DIR / "emi12822-sup-0006-si.xls"

SYN_STRAINS = ["BL107", "CC9311", "RS9916", "RS9917", "WH5701", "WH7803", "WH7805", "WH8102"]

S6_SHEETS = {
    "4 days vs 10 hours": "10h_vs_4d",
    "4 days vs 4 days co-culture": "coculture_vs_axenic_4d",
    "10 hours vs 10 hours dark": "dark_vs_light_10h",
    "10 hours vs 10 hours phage": "phage_vs_control_10h",
}


def clean_text(v) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    s = str(v).strip()
    return s.replace("|", ",").replace("'", "^")


def clean_seq(v) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    return re.sub(r"[^A-Za-z]", "", str(v)).upper()


def parse_s4(seqs: dict[str, dict[str, str]]) -> pd.DataFrame:
    raw = pd.read_excel(S4, header=None)
    strain_row = raw.iloc[2].tolist()
    hdr = raw.iloc[3].tolist()
    body = raw.iloc[4:]
    id_cols = [i for i, h in enumerate(hdr) if h == "ID"]
    records = []
    for k, i in enumerate(id_cols):
        strain = str(strain_row[i]).strip()
        if strain not in SYN_STRAINS:
            continue  # Prochlorococcus / SAR11 blocks: in-silico homologues only, no NSAF
        end = id_cols[k + 1] if k + 1 < len(id_cols) else len(hdr)
        block = {hdr[j]: j for j in range(i, end)}
        rep_cols = [j for j in range(i + 5, end) if hdr[j] not in ("average",) and pd.notna(hdr[j])]
        avg_col = block["average"]
        for _, row in body.iterrows():
            pid = clean_text(row.iloc[i])
            if not pid:
                continue
            seq = clean_seq(row.iloc[i + 2])
            if seq:
                seqs.setdefault(strain, {})[pid] = seq
            reps = {hdr[j]: pd.to_numeric(row.iloc[j], errors="coerce") for j in rep_cols}
            reps = {r: (0.0 if pd.isna(v) else float(v)) for r, v in reps.items()}
            if strain == "BL107":  # collapse technical re-runs A1/A2, B1/B2
                bio = {}
                for r, v in reps.items():
                    bio[r[0]] = max(bio.get(r[0], 0.0), v)
                reps = bio
            n_rep = len(reps)
            n_det = sum(1 for v in reps.values() if v > 0)
            avg = pd.to_numeric(row.iloc[avg_col], errors="coerce")
            avg = 0.0 if pd.isna(avg) else float(avg)
            if n_det == n_rep and n_rep > 0:
                status = "detected"
            elif n_det > 0:
                status = "sporadic"
            else:
                status = "not_detected"
            records.append({
                "strain": strain,
                "protein_id": pid,
                "annotated_function": clean_text(row.iloc[i + 1]),
                "avg_pct_nsaf": avg if avg > 0 else None,
                "detection_status": status,
                "n_replicates": n_rep,
                "n_replicates_detected": n_det,
            })
    df = pd.DataFrame(records)
    # A protein can sit in several homology-cluster rows; keep its strongest record.
    df["_rank"] = df["avg_pct_nsaf"].fillna(0)
    df = (df.sort_values(["strain", "protein_id", "_rank"], ascending=[True, True, False])
            .drop_duplicates(["strain", "protein_id"]))
    # Identical-sequence paralogs (e.g. RS9917's two cpcBA copies) cannot be told apart by
    # MS and collapse onto one current gene: keep only the most abundant copy per sequence.
    df["_seq"] = [seqs.get(s, {}).get(p, p) for s, p in zip(df["strain"], df["protein_id"])]
    before = len(df)
    df = (df.sort_values(["strain", "_seq", "_rank"], ascending=[True, True, False])
            .drop_duplicates(["strain", "_seq"]).drop(columns=["_rank", "_seq"]))
    print(f"  dropped {before - len(df)} identical-sequence duplicate rows")
    return df.sort_values(["strain", "avg_pct_nsaf"], ascending=[True, False]).reset_index(drop=True)


def parse_s5_seqs(seqs: dict[str, dict[str, str]]) -> None:
    df = pd.read_excel(S5, sheet_name="Protein list", header=2)
    for _, r in df.iterrows():
        pid = clean_text(r.get("NCBI reference"))
        seq = clean_seq(r.get("Protein sequence"))
        if pid and seq:
            seqs.setdefault("WH8102", {})[pid] = seq


def parse_s6(seqs: dict[str, dict[str, str]]) -> dict[str, pd.DataFrame]:
    out = {}
    for sheet, slug in S6_SHEETS.items():
        raw = pd.read_excel(S6, sheet_name=sheet, header=None)
        section = ""
        rows = []
        for idx in range(4, len(raw)):
            r = raw.iloc[idx]
            label = clean_text(r.iloc[0])
            # Section headers: "Significantly up/downregulated ..." (followed by a
            # "p-value < 0.05 and fold change > 2x" sub-line that belongs to it),
            # then "p-value < 0.05 but fold change < 2x", then "p-value > 0.05".
            if label and not label.lower().startswith("p-value < 0.05 and"):
                section = label
            pid = clean_text(r.iloc[1])
            if not pid:
                continue
            fc = pd.to_numeric(r.iloc[4], errors="coerce")
            p = pd.to_numeric(r.iloc[5], errors="coerce")
            seq = clean_seq(r.iloc[3])
            if seq:
                seqs.setdefault("WH8102", {})[pid] = seq
            if pd.isna(fc) or fc == 0:
                l2 = None
            else:
                l2 = math.copysign(math.log2(abs(fc)), fc)
            rows.append({
                "protein_id": pid,
                "annotated_function": clean_text(r.iloc[2]),
                "fold_change": fc,
                "log2_fold_change": l2,
                "p_value": p,
                "section": section,
            })
        out[slug] = pd.DataFrame(rows)
    return out


def write_fasta(strain: str, prot: dict[str, str]) -> Path:
    path = PAPER_DIR / f"protein_seqs_{strain}.faa"
    with open(path, "w") as f:
        for pid in sorted(prot):
            f.write(f">{pid}\n{prot[pid]}\n")
    return path


def main() -> None:
    seqs: dict[str, dict[str, str]] = {}
    s4 = parse_s4(seqs)
    for strain, g in s4.groupby("strain"):
        out = PAPER_DIR / f"table_s4_{strain}_exoproteome_nsaf.csv"
        g.to_csv(out, index=False)
        print(f"{S4.name} -> {out.name}: {len(g)} rows")
    print(s4.groupby(["strain", "detection_status"]).size().unstack(fill_value=0).to_string())

    parse_s5_seqs(seqs)
    for slug, df in parse_s6(seqs).items():
        p = PAPER_DIR / f"table_s6_{slug}.csv"
        df.to_csv(p, index=False)
        print(f"{S6.name} -> {p.name}: {len(df)} rows; sections {df['section'].value_counts().to_dict()}")

    for strain in SYN_STRAINS:
        fa = write_fasta(strain, seqs.get(strain, {}))
        print(f"  {fa.name}: {len(seqs.get(strain, {}))} sequences")


if __name__ == "__main__":
    main()
