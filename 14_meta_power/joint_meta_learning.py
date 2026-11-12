#!/usr/bin/env python3
"""Stage 14: joint analysis of the two outcome-labelled cohorts, learning curve, seed stability,
and directional concordance of the pathway-level signal.

(1) Per-gene recovery effect sizes (Hedges g) in GSE328059 (n = 84) and PRJNA1312133 (n = 19),
    pooled by DerSimonian-Laird random effects with Benjamini-Hochberg correction.
(2) Learning curve: cross-validated AUC against training-set size, with extrapolation.
(3) Seed stability: the whole pipeline repeated under 20 master seeds.
(4) Directional concordance in the external cohort of the sets that reach FDR < 0.05 internally.

usage: python joint_meta_learning.py [../data/prjna]
"""
import sys, json, gzip, re
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "05_recovery"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "02_meta_lasso_shap"))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
from sklearn.feature_selection import f_classif
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split
from pipeline import load_pre, adjust_batch, make_model
from common import COH, RES, SEED, bh

FIG = Path(__file__).resolve().parent.parent / "figures"
plt.rcParams.update({"font.family": "serif", "font.size": 10, "axes.linewidth": 1.1})

X, y, b, m = load_pre(); Xv = X.values
full = pd.read_csv(COH / "GSE328059_expr.csv.gz", index_col=0)[m["sample"].astype(str)]


def hedges(a, c):
    n1, n2 = len(a), len(c)
    s = np.sqrt(((n1 - 1) * a.var(axis=0, ddof=1) + (n2 - 1) * c.var(axis=0, ddof=1)) / (n1 + n2 - 2)) + 1e-9
    d = (a.mean(0) - c.mean(0)) / s
    J = 1 - 3 / (4 * (n1 + n2) - 9)
    g = J * d
    se = np.sqrt((n1 + n2) / (n1 * n2) + g ** 2 / (2 * (n1 + n2 - 2)))
    return g, se


def dl_pool(g1, s1, g2, s2):
    G = np.vstack([g1, g2]); S = np.vstack([s1, s2]) ** 2
    w = 1 / S
    fixed = (w * G).sum(0) / w.sum(0)
    Q = (w * (G - fixed) ** 2).sum(0)
    c = w.sum(0) - (w ** 2).sum(0) / w.sum(0)
    tau2 = np.maximum(0, (Q - 1) / np.maximum(c, 1e-9))
    wr = 1 / (S + tau2)
    gp = (wr * G).sum(0) / wr.sum(0)
    sep = np.sqrt(1 / wr.sum(0))
    return gp, sep, tau2, Q


summary = {}
# ------------------------------------------------- (1) joint meta-analysis
if len(sys.argv) > 1 and Path(sys.argv[1]).exists():
    Q = Path(sys.argv[1]) / "quant"; hdr = {}
    for l in gzip.open(Q.parent / "ref" / "gencode.v46.transcripts.fa.gz", "rt"):
        if l[0] == ">":
            f = l[1:].split("|"); hdr[f[0]] = f[5]
    rw = {}
    for d in sorted(Q.iterdir()):
        if re.search(r"_HF_(R|NR)\d+$", d.name):
            t = pd.read_csv(d / "quant.sf", sep="\t", index_col=0); t["g"] = t.index.map(hdr)
            rw[d.name] = np.log2(t.groupby("g").TPM.sum() + 1)
    E = pd.DataFrame(rw); E = E.loc[E.std(axis=1) > 0]
    yext = np.array([int(re.search(r"_HF_(R|NR)", c).group(1) == "R") for c in E.columns])
    shared = full.index.intersection(E.index)
    A = full.loc[shared]; B = E.loc[shared]
    Az = A.sub(A.mean(1), axis=0).div(A.std(1) + 1e-9, axis=0).T.values
    Bz = B.sub(B.mean(1), axis=0).div(B.std(1) + 1e-9, axis=0).T.values
    g1, s1 = hedges(Az[y == 1], Az[y == 0])
    g2, s2 = hedges(Bz[yext == 1], Bz[yext == 0])
    gp, sep, tau2, Qs = dl_pool(g1, s1, g2, s2)
    zs = gp / sep; p = 2 * stats.norm.sf(abs(zs))
    M = pd.DataFrame({"gene": shared, "g_GSE328059": g1, "g_PRJNA": g2, "g_pooled": gp,
                      "se_pooled": sep, "z": zs, "p": p, "fdr": bh(p), "tau2": tau2,
                      "concordant": np.sign(g1) == np.sign(g2)}).sort_values("p")
    M.to_csv(RES / "joint_recovery_meta.csv", index=False)
    conc = float(M.concordant.mean())
    bp = stats.binomtest(int(M.concordant.sum()), len(M), 0.5).pvalue
    top = M.head(500)
    summary["joint_meta"] = dict(n_genes=int(len(M)), n_fdr05=int((M.fdr < 0.05).sum()),
                                 n_fdr10=int((M.fdr < 0.10).sum()), min_fdr=float(M.fdr.min()),
                                 sign_concordance=conc, concordance_binomial_p=float(bp),
                                 sign_concordance_top500=float(top.concordant.mean()),
                                 spearman_g=float(stats.spearmanr(g1, g2).statistic),
                                 top_genes=M.head(10)[["gene", "g_pooled", "fdr"]].to_dict("records"))
    print("joint meta:", json.dumps(summary["joint_meta"], indent=1, default=str)[:600])

# ------------------------------------------------- (2) learning curve
def cv_auc_sub(idx, seed, repeats=3):
    Xs, ys, bs = Xv[idx], y[idx], b[idx]
    a = []
    for r in range(repeats):
        oof = np.zeros(len(ys))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed + r).split(Xs, ys):
            Xtr, Xte = adjust_batch(Xs[tr], ys[tr], bs[tr], Xs[te], bs[te])
            mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
            Ztr, Zte = (Xtr - mu) / sd, (Xte - mu) / sd
            F_, _ = f_classif(Ztr, ys[tr]); top = np.argsort(np.nan_to_num(F_))[::-1][:200]
            oof[te] = make_model("enet", seed + r).fit(Ztr[:, top], ys[tr]).decision_function(Zte[:, top])
        a.append(roc_auc_score(ys, oof))
    return float(np.mean(a))


sizes = [30, 40, 50, 60, 70, 84]
rows = []
for n in sizes:
    vals = []
    for rep in range(10 if n < 84 else 1):
        if n == len(y): idx = np.arange(len(y))
        else:
            idx, _ = train_test_split(np.arange(len(y)), train_size=n, stratify=y, random_state=SEED + rep)
        vals.append(cv_auc_sub(idx, SEED + rep))
    rows.append(dict(n=n, auc_mean=float(np.mean(vals)), auc_sd=float(np.std(vals))))
LC = pd.DataFrame(rows); LC.to_csv(RES / "learning_curve.csv", index=False)
print(LC.round(3).to_string(index=False))


def inv(n, a, c):
    return a - c / np.sqrt(n)


try:
    popt, _ = curve_fit(inv, LC.n, LC.auc_mean, p0=[0.8, 1.0], maxfev=10000)
    ceiling = float(popt[0])
    need = {t: (float((popt[1] / max(popt[0] - t, 1e-6)) ** 2) if popt[0] > t else None) for t in [0.75, 0.80, 0.85]}
except Exception:
    ceiling, need = None, {}
summary["learning_curve"] = dict(points=LC.round(3).to_dict("records"), fitted_ceiling=ceiling,
                                 n_needed=need)
print("ceiling", ceiling, "n needed", need)

# ------------------------------------------------- (3) seed stability
seeds = [SEED + 7919 * i for i in range(20)]
st = [cv_auc_sub(np.arange(len(y)), s, repeats=2) for s in seeds]
summary["seed_stability"] = dict(n_seeds=len(seeds), mean=float(np.mean(st)), sd=float(np.std(st)),
                                 min=float(np.min(st)), max=float(np.max(st)))
print("seed stability", summary["seed_stability"])

# ------------------------------------------------- (4) pathway directional concordance
pw = RES / "pathway_univariate.csv"
if pw.exists() and "yext" in dir():
    import gseapy
    GSJ = json.load(open(Path(__file__).resolve().parent.parent / "11_deconvolution_pathways" / "genesets.json"))
    U = pd.read_csv(pw)
    sig = U[U.fdr < 0.05]
    r = gseapy.ssgsea(data=E, gene_sets=GSJ, outdir=None, sample_norm_method="rank", no_plot=True,
                      threads=8, min_size=10, max_size=500)
    Se = r.res2d.pivot(index="Term", columns="Name", values="NES").astype(float)[E.columns]
    common = [s for s in sig.set.values if s in Se.index]
    te, pe = stats.ttest_ind(Se.loc[common].T.values[yext == 1], Se.loc[common].T.values[yext == 0], equal_var=False)
    ti = sig.set_index("set").loc[common, "t"].values
    conc = float((np.sign(te) == np.sign(ti)).mean())
    bp = stats.binomtest(int((np.sign(te) == np.sign(ti)).sum()), len(common), 0.5).pvalue
    pd.DataFrame({"set": common, "t_internal": ti, "t_external": te, "p_external": pe,
                  "concordant": np.sign(te) == np.sign(ti)}).to_csv(RES / "pathway_external_concordance.csv", index=False)
    summary["pathway_concordance"] = dict(n_sets=len(common), concordance=conc, binomial_p=float(bp),
                                          n_external_p05=int((pe < 0.05).sum()))
    print("pathway concordance", summary["pathway_concordance"])

json.dump(summary, open(RES / "joint_meta_learning_summary.json", "w"), indent=1, default=str)

fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.3))
if "joint_meta" in summary:
    ax[0].scatter(g1, g2, s=3, color="#aec3df", alpha=.5)
    sel = M.head(50).gene.values; ii = [list(shared).index(g) for g in sel]
    ax[0].scatter(g1[ii], g2[ii], s=8, color="#b03a2e")
    ax[0].axhline(0, color="k", lw=.8); ax[0].axvline(0, color="k", lw=.8)
    ax[0].set_xlabel("Hedges g, GSE328059 (n = 84)"); ax[0].set_ylabel("Hedges g, PRJNA1312133 (n = 19)")
    ax[0].set_title(f"Per-gene recovery effects\nSpearman r = {stats.spearmanr(g1, g2).statistic:.2f}", fontsize=9)
ax[1].errorbar(LC.n, LC.auc_mean, yerr=LC.auc_sd, fmt="o-", color="#1f4e9e", capsize=3)
if ceiling:
    gx = np.linspace(25, 600, 200); ax[1].plot(gx, inv(gx, *popt), "--", color="#b03a2e", lw=1.2)
    ax[1].axhline(ceiling, color="grey", ls=":", lw=1)
    ax[1].set_xscale("log")
ax[1].axhline(0.5, color="k", ls="--", lw=1); ax[1].set_xlabel("Training-set size"); ax[1].set_ylabel("Cross-validated AUC")
ax[1].set_title("Learning curve", fontsize=10)
ax[2].hist(st, bins=10, color="#2e86c1", edgecolor="k", lw=.5)
ax[2].axvline(np.mean(st), color="k", ls="--", lw=1.2)
ax[2].set_xlabel("Cross-validated AUC"); ax[2].set_ylabel("Master seeds"); ax[2].set_title("Seed stability (20 runs)", fontsize=10)
for a, l in zip(ax, "abc"):
    a.text(-0.12, 1.05, f"({l})", transform=a.transAxes, fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG / "FigS20_joint_meta_learning.png", dpi=600)
fig.savefig(FIG / "FigS20_joint_meta_learning.tif", dpi=600, pil_kwargs={"compression": "tiff_lzw", "dpi": (600, 600)})
print("wrote", FIG / "FigS20_joint_meta_learning.tif")
