"""Write barreto 2022 `<table>_modified.csv` with the `uniprot_acc` column blanked on
rows where it is not row-aligned (the source column is shifted by one row in
several blocks, e.g. PMT9312_1637/rpsH carried Q318J8 = rplF).

A row keeps its accession iff (a) it equals the KG's own UniProt join for that
gene, or (b) the gene has no KG accession and the row's `uniprot_def` shares
more informative words with the row's own `PRODUCT` than with any product on
the two rows either side (ties count as shifts). Rows whose texts are too
generic to score are kept. Re-run after a UniProt cache refresh.
"""
import json, re, sys, pathlib
import pandas as pd

STOP = {"protein", "putative", "probable", "family", "domain", "containing", "hypothetical",
        "uncharacterized", "unknown", "function", "the", "and", "type", "like", "subfamily",
        "conserved", "predicted", "possible"}
TABLES = [("pro_9312_anot.csv", "cache/data/Prochlorococcus/genomes/MIT9312"),
          ("syn_8102_anot.csv", "cache/data/Synechococcus/genomes/WH8102")]
ROOT = pathlib.Path(__file__).resolve().parent.parent
PAPER = ROOT / "data/Prochlorococcus/papers_and_supp/barreto 2022"


def toks(s):
    return {w for w in re.findall(r"[a-z0-9]+", str(s).lower()) if len(w) >= 2 and w not in STOP}


def score(a, b):
    ta, tb = toks(a), toks(b)
    return None if not ta or not tb else len(ta & tb)


def main():
    for src, genome_dir in TABLES:
        p = PAPER / src
        df = pd.read_csv(p, low_memory=False, dtype=str)
        merged = json.load(open(ROOT / genome_dir / "gene_annotations_merged.json"))
        lt_col = df.columns[0]
        c = {"kept_ref": 0, "blank_ref_mismatch": 0, "kept_text": 0, "blank_shift": 0,
             "blank_disagree": 0, "kept_unknown": 0}
        for i, r in df.iterrows():
            acc = r.get("uniprot_acc")
            if pd.isna(acc) or not str(acc).strip():
                continue
            g = merged.get(str(r[lt_col]).strip(), {})
            if g.get("uniprot_accession"):
                if g["uniprot_accession"] == acc:
                    c["kept_ref"] += 1
                else:
                    df.at[i, "uniprot_acc"] = ""; c["blank_ref_mismatch"] += 1
                continue
            own = score(r.get("PRODUCT"), r.get("uniprot_def"))
            nb = [s for s in (score(df.at[j, "PRODUCT"], r.get("uniprot_def"))
                              for j in (i - 2, i - 1, i + 1, i + 2) if 0 <= j < len(df)) if s is not None]
            if own is None:
                c["kept_unknown"] += 1
            elif nb and max(nb) >= own:          # a neighbour fits as well or better: shifted
                df.at[i, "uniprot_acc"] = ""; c["blank_shift"] += 1
            elif own >= 1:
                c["kept_text"] += 1
            else:
                df.at[i, "uniprot_acc"] = ""; c["blank_disagree"] += 1
        out = p.with_name(p.stem + "_modified.csv")
        df.to_csv(out, index=False)
        print(out.name, c)


if __name__ == "__main__":
    main()
