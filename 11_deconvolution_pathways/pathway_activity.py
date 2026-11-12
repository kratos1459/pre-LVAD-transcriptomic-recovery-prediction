#!/usr/bin/env python3
"""Stage 11b: pathway- and transcription-factor-level recovery modelling.

Single-sample GSEA (gseapy) scores every pre-LVAD heart on MSigDB Hallmark gene sets and TRRUST
transcription-factor target sets. The same leakage-controlled pipeline as stage 05 is then run on
these low-dimensional, rank-based features, and the frozen model is applied once to PRJNA1312133.
Rationale: rank-based set scores are less sensitive to batch and platform than individual genes, so
if the gene-level signal is largely technical the set-level model should transfer better.

usage: python pathway_activity.py [../data/prjna]
"""
import sys, json, gzip, re
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "05_recovery"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "02_meta_lasso_shap"))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import gseapy
from scipy import stats
from sklearn.feature_selection import f_classif
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from pipeline import load_pre, adjust_batch, make_model
from common import COH, RES, SEED, bh

FIG = Path(__file__).resolve().parent.parent / "figures"
CACHE = Path(__file__).resolve().parent / "genesets.json"
plt.rcParams.update({"font.family": "serif", "font.size": 10, "axes.linewidth": 1.1})
REP, NSEL = 10, 40
LIBS = ["MSigDB_Hallmark_2020", "TRRUST_Transcription_Factors_2019"]

if CACHE.exists():
    gs = json.load(open(CACHE))
else:
    gs = {}
    for lib in LIBS:
        d = gseapy.get_library(name=lib, organism="Human")
        for k, v in d.items():
            if 10 <= len(v) <= 500:
                gs[f"{lib.split('_')[0]}|{k}"] = v
    json.dump(gs, open(CACHE, "w"))
print("gene sets:", len(gs))


def ssgsea(expr):
    r = gseapy.ssgsea(data=expr, gene_sets=gs, outdir=None, sample_norm_method="rank",
                      no_plot=True, threads=8, min_size=10, max_size=500)
    t = r.res2d.pivot(index="Term", columns="Name", values="NES").astype(float)
    return t[expr.columns]


X, y, b, m = load_pre()
samples = m["sample"].astype(str).values
expr = pd.read_csv(COH / "GSE328059_expr.csv.gz", index_col=0)[samples]
expr = expr.loc[expr.std(axis=1) > 0]
Sc = ssgsea(expr)
Sc.to_csv(RES / "pathway_scores_GSE328059.csv")
P = Sc.T.values
print("score matrix", P.shape)


def cv(Pv, y, b, repeats=REP, nsel=NSEL, adjust=True):
    aucs, oof_all = [], []
    for r in range(repeats):
        oof = np.zeros(len(y))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEED + r).split(Pv, y):
            Ptr, Pte = (adjust_batch(Pv[tr], y[tr], b[tr], Pv[te], b[te]) if adjust else (Pv[tr], Pv[te]))
            mu, sd = Ptr.mean(0), Ptr.std(0) + 1e-9
            Ztr, Zte = (Ptr - mu) / sd, (Pte - mu) / sd
            F_, _ = f_classif(Ztr, y[tr]); top = np.argsort(np.nan_to_num(F_))[::-1][:nsel]
            mo = make_model("enet", SEED + r).fit(Ztr[:, top], y[tr])
            oof[te] = mo.decision_function(Zte[:, top])
        aucs.append(roc_auc_score(y, oof)); oof_all.append(oof)
    return np.array(aucs), np.mean(oof_all, axis=0)


auc_adj, oof = cv(P, y, b)
auc_raw, _ = cv(P, y, b, adjust=False)
print(f"pathway-level CV AUC: batch-adjusted {auc_adj.mean():.3f} +/- {auc_adj.std():.3f}; unadjusted {auc_raw.mean():.3f}")

# held-out batch
hb = {}
for held in ["11943", "10931"]:
    te = b == held
    Ptr, Pte = adjust_batch(P[~te], y[~te], b[~te], P[te], b[te])
    mu, sd = Ptr.mean(0), Ptr.std(0) + 1e-9
    Ztr, Zte = (Ptr - mu) / sd, (Pte - mu) / sd
    F_, _ = f_classif(Ztr, y[~te]); top = np.argsort(np.nan_to_num(F_))[::-1][:NSEL]
    mo = make_model("enet", SEED).fit(Ztr[:, top], y[~te])
    hb[held] = float(roc_auc_score(y[te], mo.decision_function(Zte[:, top])))
print("held-out batch:", {k: round(v, 3) for k, v in hb.items()})

# univariate
t, p = stats.ttest_ind(P[y == 1], P[y == 0], equal_var=False)
uv = pd.DataFrame({"set": Sc.index, "t": t, "p": p, "fdr": bh(np.nan_to_num(p, nan=1))}).sort_values("p")
uv.to_csv(RES / "pathway_univariate.csv", index=False)
print(uv.head(8).round(4).to_string(index=False))

# frozen external test
ext = {}
if len(sys.argv) > 1 and Path(sys.argv[1]).exists():
    Q = Path(sys.argv[1]) / "quant"
    hdr = {}
    for l in gzip.open(Q.parent / "ref" / "gencode.v46.transcripts.fa.gz", "rt"):
        if l[0] == ">":
            f = l[1:].split("|"); hdr[f[0]] = f[5]
    rows = {}
    for d in sorted(Q.iterdir()):
        if re.search(r"_HF_(R|NR)\d+$", d.name):
            tt = pd.read_csv(d / "quant.sf", sep="\t", index_col=0); tt["g"] = tt.index.map(hdr)
            rows[d.name] = np.log2(tt.groupby("g").TPM.sum() + 1)
    E = pd.DataFrame(rows); E = E.loc[E.std(axis=1) > 0]
    yext = np.array([int(re.search(r"_HF_(R|NR)", c).group(1) == "R") for c in E.columns])
    Se = ssgsea(E)
    common = Sc.index.intersection(Se.index)
    Pd_, Pe = Sc.loc[common].T.values, Se.loc[common].T.values
    Pe = (Pe - Pe.mean(0)) / (Pe.std(0) + 1e-9) * (Pd_.std(0) + 1e-9) + Pd_.mean(0)
    Ptr, Pte = adjust_batch(Pd_, y, b, Pe, np.array(["other"] * len(Pe)))
    mu, sd = Ptr.mean(0), Ptr.std(0) + 1e-9
    Ztr, Zte = (Ptr - mu) / sd, (Pte - mu) / sd
    F_, _ = f_classif(Ztr, y); top = np.argsort(np.nan_to_num(F_))[::-1][:NSEL]
    mo = make_model("enet", SEED).fit(Ztr[:, top], y)
    sc = mo.decision_function(Zte[:, top])
    auc_e = roc_auc_score(yext, sc)
    rng = np.random.default_rng(SEED); boots = []
    while len(boots) < 5000:
        i = rng.integers(0, len(yext), len(yext))
        if 0 < yext[i].sum() < len(i): boots.append(roc_auc_score(yext[i], sc[i]))
    pe = stats.mannwhitneyu(sc[yext == 1], sc[yext == 0], alternative="two-sided", method="exact").pvalue
    ext = dict(auc=float(auc_e), ci=[float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))], exact_p=float(pe))
    pd.DataFrame({"sample": E.columns, "y": yext, "score": sc}).to_csv(RES / "pathway_external_scores.csv", index=False)
    print("external pathway-level:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in ext.items()})

summary = dict(n_sets=int(len(Sc)), cv_auc_batch_adjusted=[round(float(auc_adj.mean()), 3), round(float(auc_adj.std()), 3)],
               cv_auc_unadjusted=round(float(auc_raw.mean()), 3), held_out_batch={k: round(v, 3) for k, v in hb.items()},
               n_sets_fdr05=int((uv.fdr < 0.05).sum()), top_sets=uv.head(5)["set"].tolist(), external=ext)
json.dump(summary, open(RES / "pathway_summary.json", "w"), indent=1)
print(json.dumps(summary, indent=1))

fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.4))
top10 = uv.head(10).iloc[::-1]
ax[0].barh(range(10), -np.log10(top10.p), color=np.where(top10.t > 0, "#1f4e9e", "#b03a2e"), edgecolor="k", lw=.5)
ax[0].set_yticks(range(10)); ax[0].set_yticklabels([s.split("|")[-1][:34] for s in top10["set"]], fontsize=6.5)
ax[0].set_xlabel("$-\\log_{10}$ p (responder vs non-responder)"); ax[0].set_title("Most associated sets", fontsize=10)
bars = [auc_adj.mean(), hb["10931"], hb["11943"]] + ([ext["auc"]] if ext else [])
labs = ["cross-\nvalidated", "held-out\nbatch 10931", "held-out\nbatch 11943"] + (["external\n(n = 19)"] if ext else [])
ax[1].bar(range(len(bars)), bars, color=["#1f4e9e", "#2e86c1", "#2e86c1", "#e67e22"][:len(bars)], edgecolor="k", lw=.6)
ax[1].errorbar([0], [auc_adj.mean()], yerr=[auc_adj.std()], fmt="none", ecolor="k", capsize=4)
ax[1].axhline(0.5, color="k", ls="--", lw=1); ax[1].set_ylim(0.3, 0.9)
ax[1].set_xticks(range(len(bars))); ax[1].set_xticklabels(labs, fontsize=7.5); ax[1].set_ylabel("AUC")
ax[1].set_title("Pathway-level recovery model", fontsize=10)
ax[2].hist([oof[y == 0], oof[y == 1]], bins=12, color=["#b03a2e", "#1f4e9e"], edgecolor="k", lw=.4,
           label=["non-responder", "responder"])
ax[2].set_xlabel("Out-of-fold score"); ax[2].set_ylabel("Hearts"); ax[2].legend(fontsize=7); ax[2].set_title("Score distribution", fontsize=10)
for a, l in zip(ax, "abc"):
    a.text(-0.12, 1.05, f"({l})", transform=a.transAxes, fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG / "FigS18_pathway_activity.png", dpi=600)
fig.savefig(FIG / "FigS18_pathway_activity.tif", dpi=600, pil_kwargs={"compression": "tiff_lzw", "dpi": (600, 600)})
print("wrote", FIG / "FigS18_pathway_activity.tif")
