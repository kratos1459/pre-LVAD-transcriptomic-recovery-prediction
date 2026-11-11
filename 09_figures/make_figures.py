#!/usr/bin/env python3
"""Stage 09: regenerate main figures from results/ (600 dpi PNG, white background).
Fig1 cohorts | Fig2 volcano | Fig4 top-30 | Fig5 LOCO | Fig6 novelty | Fig10 recovery CV+permutation | Fig11 held-out batch | Fig12 ROC | Fig13 tertiles"""
import sys
sys.path.insert(0, "../02_meta_lasso_shap")
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
from common import RES, COH

OUT = RES.parent / "figures"; OUT.mkdir(exist_ok=True)
B, G, O, R, K, GR = "#1D4ED8", "#059669", "#D97706", "#B91C1C", "#1F2937", "#9CA3AF"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.linestyle": "--", "grid.alpha": .35, "axes.titleweight": "bold", "axes.titlelocation": "left"})


def save(f, n):
    f.savefig(OUT / f"{n}.png", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(f)


m = pd.read_csv(COH / "all_meta.csv"); m = m[m.group.isin(["HF", "NF"])]
t = m.groupby(["cohort", "group"]).size().unstack(fill_value=0)
f, a = plt.subplots(figsize=(7, 4)); t[["HF", "NF"]].plot.barh(stacked=True, ax=a, color=[R, G], width=.7)
a.set_xlabel("Samples"); a.set_ylabel(""); a.set_title("Cohorts used for HF gene discovery"); a.legend(["Heart failure", "Non-failing"], frameon=False); save(f, "Fig1_cohorts")

d = pd.read_csv(RES / "meta_consensus.csv"); d["nl"] = -np.log10(d.p.clip(lower=1e-300))
sig = (d.fdr < 1e-6) & (d.g.abs() > .8)
f, a = plt.subplots(figsize=(6.5, 5)); a.scatter(d.g[~sig], d.nl[~sig], s=4, c=GR, alpha=.5)
a.scatter(d.g[sig & (d.g > 0)], d.nl[sig & (d.g > 0)], s=7, c=R); a.scatter(d.g[sig & (d.g < 0)], d.nl[sig & (d.g < 0)], s=7, c=B)
for _, r in d.dropna(subset=["consensus_rank_score"]).head(8).iterrows():
    a.annotate(r.gene, (r.g, r.nl), fontsize=8)
a.set_xlabel("Pooled Hedges g (HF minus non-failing)"); a.set_ylabel("-log10 p"); a.set_title("Random-effects meta-analysis"); save(f, "Fig2_volcano")

t30 = d.dropna(subset=["consensus_rank_score"]).head(30)[::-1]
f, a = plt.subplots(figsize=(6, 7))
sc = a.scatter(t30.lasso_freq, range(len(t30)), s=t30.mean_abs_shap * 900 + 20, c=t30.g, cmap="coolwarm", vmin=-2.5, vmax=2.5, edgecolor=K, lw=.4)
a.set_yticks(range(len(t30))); a.set_yticklabels(t30.gene, fontsize=8); a.set_xlabel("Bootstrap LASSO selection frequency")
plt.colorbar(sc, label="Hedges g"); a.set_title("Top 30 consensus genes"); save(f, "Fig4_top30")

l = pd.read_csv(RES / "loco_auc.csv"); l = l[l.variant == "standard"]
f, (a, b) = plt.subplots(1, 2, figsize=(10, 4))
for mod, c in zip(["L1", "L2", "RF", "XGB"], [B, G, O, R]):
    s = l[l.model == mod].groupby("N").auc.mean(); a.plot(s.index, s.values, "-o", c=c, label=mod)
a.set_xscale("log"); a.set_xlabel("Panel size (genes)"); a.set_ylabel("Mean held-out AUC"); a.legend(frameon=False); a.set_title("LOCO vs panel size")
s = l[(l.model == "L2") & (l.N == 20)]; b.bar(s.held_out, s.auc, color=B); b.set_ylim(.85, 1.01); b.tick_params(axis="x", rotation=60)
b.set_title("20-gene ridge AUC per cohort"); save(f, "Fig5_LOCO")

p = RES / "top150_annotated.csv"
if p.exists():
    n = pd.read_csv(p); f, a = plt.subplots(figsize=(6, 4.5))
    a.scatter(n.consensus_rank_score, np.log10(n.pubmed_hf_papers + 1), c=np.where(n.understudied, R, GR), s=14)
    a.set_xlabel("Consensus rank score"); a.set_ylabel("log10(PubMed HF papers + 1)"); a.set_title("Rank vs literature volume"); save(f, "Fig6_novelty")

cv = pd.read_csv(RES / "recovery_cv_auc.csv"); null = np.load(RES / "recovery_perm_null.npy")
f, (a, b) = plt.subplots(1, 2, figsize=(10, 4)); mm = cv.mean()
a.barh(mm.index[::-1], mm.values[::-1], xerr=cv.std().values[::-1], color=O); a.set_xlim(.5, .85); a.set_xlabel("Repeated 5-fold CV AUC"); a.set_title("Models")
b.hist(null, 30, color=GR); b.axvline(cv["enet|inner"].mean(), c=R, lw=2); b.set_xlabel("Permutation AUC"); b.set_title("500-label permutation null"); save(f, "Fig10_recovery_cv")

hb = pd.read_csv(RES / "recovery_heldout_batch.csv", index_col=0).iloc[:, 0]
hb.index = pd.MultiIndex.from_tuples([i.split("|") for i in hb.index]); hb = hb.unstack(0)
f, a = plt.subplots(figsize=(6, 4)); hb.plot.bar(ax=a, color=[B, O]); a.axhline(.5, c=K, ls=":"); a.set_ylabel("AUC on held-out batch"); a.set_title("Held-out sequencing batch"); save(f, "Fig11_heldout_batch")

o = pd.read_csv(RES / "recovery_oof.csv"); fp, tp, _ = roc_curve(o.y, o.oof_enet)
f, a = plt.subplots(figsize=(4.5, 4.5)); a.plot(fp, tp, c=B, lw=2); a.plot([0, 1], [0, 1], c=GR, ls="--")
a.set_xlabel("1 - specificity"); a.set_ylabel("Sensitivity"); a.set_title(f"OOF ROC (AUC {roc_auc_score(o.y, o.oof_enet):.2f})"); save(f, "Fig12_ROC")
o["tert"] = pd.qcut(o.oof_enet.rank(method="first"), 3, labels=["Low", "Mid", "High"]); r = o.groupby("tert", observed=True).y.mean() * 100
f, a = plt.subplots(figsize=(4.5, 4)); a.bar(r.index.astype(str), r.values, color=[B, GR, R]); a.set_ylabel("Responder rate (%)"); a.set_title("Risk tertiles"); save(f, "Fig13_tertiles")
print("figures:", sorted(p.name for p in OUT.glob("*.png")))
