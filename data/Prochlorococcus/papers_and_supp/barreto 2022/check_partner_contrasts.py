"""Run:
    uv run --no-project --with rdata --with pandas --with openpyxl python -I check_partner_contrasts.py \
        Alteromonas_analysis/Alteromonas_analysis transcriptome_supplemental_tables/Supplemental_tables

Recompute Barreto 2022 EZ55 edgeR contrasts (approx., saturated cell-means model) and
match them against the published partner tables S9-S14.

Saturated design ~0+cyano+treatment+cyano:treatment => fitted cell means are (near) independent
of dispersion; logFC computed from TMM-normalised pooled cell means with edgeR-like prior count.
"""
import sys, warnings
import numpy as np, pandas as pd, rdata, openpyxl

warnings.simplefilter("ignore")
base, tabdir = sys.argv[1], sys.argv[2]
d = rdata.read_rda(f"{base}/ALL_EZ55_2.Rdata")
y = d["y"]
y = y if isinstance(y, dict) else y.value
c = y["counts"]  # xarray DataArray (genes x libraries)
counts = pd.DataFrame(np.asarray(c), index=list(c.coords[c.dims[0]].values),
                          columns=list(c.coords[c.dims[1]].values))
ws = openpyxl.load_workbook(f"{base}/targets.xlsx", read_only=True).worksheets[0]
rows = list(ws.iter_rows(values_only=True))
tg = pd.DataFrame(rows[1:], columns=rows[0])
assert [f.replace("_",".") for f in tg["file"]] == [c.replace("_",".") for c in counts.columns]
tg.index = counts.columns

drop = {"EZ55_00054","EZ55_00054.rRNA","EZ55_00055","EZ55_00055.tRNA","EZ55_00056","EZ55_00056.tRNA",
        "EZ55_00057","EZ55_00057.rRNA","EZ55_00058","EZ55_00058.rRNA","EZ55_00059","EZ55_00059.rRNA",
        "EZ55_00346","EZ55_00346.rRNA","EZ55_00347","EZ55_00347.rRNA","EZ55_00348","EZ55_00348.rRNA",
        "EZ55_02640","EZ55_03378","EZ55_03379","EZ55_03382","EZ55_03536","EZ55_03537","EZ55_03538",
        "EZ55_02636","EZ55_02637"}
counts = counts.loc[[g for g in counts.index if g not in drop]]
lib = counts.sum(0)
cpm = counts / lib * 1e6
counts = counts.loc[(cpm > 10 / lib.median() * 1e6).sum(1) >= 3]   # rough filterByExpr
lib = counts.sum(0)

def tmm(counts, lib):
    ref = counts.columns[np.argmin(abs((counts / lib).quantile(0.75) - (counts / lib).quantile(0.75).mean()))]
    f = {}
    for s in counts.columns:
        o, r = counts[s], counts[ref]
        ok = (o > 0) & (r > 0)
        po, pr = o[ok] / lib[s], r[ok] / lib[ref]
        M, A = np.log2(po / pr), 0.5 * np.log2(po * pr)
        keep = (M.rank(pct=True).between(0.3, 0.7)) & (A.rank(pct=True).between(0.05, 0.95))
        w = 1 / ((lib[s] - o[ok]) / lib[s] / o[ok] + (lib[ref] - r[ok]) / lib[ref] / r[ok])
        f[s] = 2 ** (np.sum(M[keep] * w[keep]) / np.sum(w[keep]))
    f = pd.Series(f); return f / np.exp(np.log(f).mean())

eff = lib * tmm(counts, lib)
cell = tg["Cyano"].astype(str) + "|" + tg["treatment"].astype(str)
print("cells:", cell.value_counts().to_dict())
pc = 0.125 * eff / eff.mean()
L = {}
for c in cell.unique():
    s = cell.index[cell == c]
    L[c] = np.log2((counts[s].add(pc[s], axis=1)).sum(1) / eff[s].sum())

cy = ["Pro", "Syn_CC9311", "Syn_W8102"]
neg, pos = "negative", "positive"
m_neg = sum(L[f"{k}|{neg}"] for k in cy) / 3
m_pos = sum(L[f"{k}|{pos}"] for k in cy) / 3
cand = {}
for k in cy:
    dev400 = L[f"{k}|{neg}"] - m_neg
    dev800 = L[f"{k}|{pos}"] - m_pos
    cand[f"{k}: dev at 400 (c3-5)"] = dev400
    cand[f"{k}: interaction (c7-9)"] = dev800 - dev400
    cand[f"{k}: dev at 800 (c3-5 + c7-9)"] = dev800
cand["co-culture vs axenic 400 (c2)"] = m_neg - L[f"None|{neg}"]
cand = pd.DataFrame(cand)

for t in ["S9", "S10", "S11", "S12", "S13", "S14"]:
    tb = pd.read_csv(f"{tabdir}/table_{t}.csv", encoding="utf-8-sig").set_index("symbol")
    j = tb.join(cand, how="inner")
    print(f"\n{t}: {len(tb)} rows, {len(j)} matched")
    for col in cand.columns:
        r = np.corrcoef(j["logFC"], j[col])[0, 1]
        mad = np.median(abs(j["logFC"] - j[col]))
        if abs(r) > 0.5:
            print(f"   r={r:+.3f}  median|diff|={mad:.3f}  {col}")
