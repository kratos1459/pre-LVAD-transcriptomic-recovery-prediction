#!/usr/bin/env python3
"""Stage 05: pre-LVAD recovery prediction (84 labelled GSE328059 pre-LVAD hearts; 26 R / 58 NR)."""
import json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from scipy import stats
from sklearn.metrics import roc_auc_score
from pipeline import *
from common import bh

X, y, b, m = load_pre()
Xv = X.values
print("n", len(y), "R", y.sum(), "NR", (1 - y).sum(), "batches", pd.Series(b).value_counts().to_dict())
tab = pd.crosstab(b, y)
chi = stats.chi2_contingency(tab.loc[["10931", "11943"]])
print("batch x outcome\n", tab, "\nchi2 p (10931 vs 11943) =", round(chi[1], 4))
out = {"chi2_p_batch": chi[1]}

# 1. repeated CV across models / adjustment variants
variants = [(mod, "inner") for mod in ["enet", "ridge", "svm", "lasso", "rf", "xgb"]] + \
           [("lasso", "none"), ("enet", "none"), ("enet", "batch_only"), ("enet", "locscale")]
res = Parallel(n_jobs=10)(delayed(cv_auc)(Xv, y, b, mod, adj, 10, return_oof=True) for mod, adj in variants)
cv = {}
for (mod, adj), (a, oof) in zip(variants, res):
    cv[f"{mod}|{adj}"] = a
    print(f"{mod:6s} {adj:10s} AUC {a.mean():.3f} +/- {a.std():.3f}")
pd.DataFrame({k: v for k, v in cv.items()}).to_csv(RES / "recovery_cv_auc.csv", index=False)
oof_enet = res[0][1].mean(axis=0)
print("pooled repeat-averaged OOF enet AUC", round(roc_auc_score(y, oof_enet), 3))
q = pd.qcut(pd.Series(oof_enet).rank(method="first"), 3, labels=["low", "mid", "high"])
tert = pd.crosstab(q, y)
print("tertile responders", tert[1].to_dict(), "of", tert.sum(axis=1).to_dict())
lo, hi = tert.loc["low"], tert.loc["high"]
orr = (hi[1] / hi[0]) / (lo[1] / lo[0])
print("OR high vs low", round(orr, 1), "Fisher p", stats.fisher_exact([[hi[1], hi[0]], [lo[1], lo[0]]])[1])
pd.DataFrame({"sample": m["sample"].values, "y": y, "batch": b, "oof_enet": oof_enet}).to_csv(RES / "recovery_oof.csv", index=False)

# 2. 500 full-pipeline label permutations (enet, inner adjustment)
obs = cv["enet|inner"].mean()


def perm(i):
    rng = np.random.default_rng(SEED + 1000 + i)
    yp = rng.permutation(y)
    return cv_auc(Xv, yp, b, "enet", "inner", 1, seed=SEED + i).mean()


null = np.array(Parallel(n_jobs=14)(delayed(perm)(i) for i in range(500)))
p = (1 + (null >= obs).sum()) / 501
print(f"perm null mean {null.mean():.3f} p95 {np.percentile(null, 95):.3f}  observed {obs:.3f}  p = {p:.4f}")
np.save(RES / "recovery_perm_null.npy", null)

# 3. held-out batch
hb = {}
for held in ["11943", "10931"]:
    te = b == held
    for mod in ["enet", "ridge", "svm", "lasso", "rf", "xgb"]:
        s = run_fold(Xv[~te], y[~te], b[~te], Xv[te], b[te], mod, "inner")
        hb[(held, mod)] = roc_auc_score(y[te], s)
print("held-out batch AUC", {f"{k[0]}:{k[1]}": round(v, 3) for k, v in hb.items()})
pd.Series({f"{k[0]}|{k[1]}": v for k, v in hb.items()}).to_csv(RES / "recovery_heldout_batch.csv")

# 4. number of selected genes, immediate-early-gene exclusion
curve = {n: cv_auc(Xv, y, b, "lasso", "inner", 5, nsel=n).mean() for n in [3, 5, 10, 25, 50, 100, 200, 500]}
print("AUC vs nsel", {k: round(v, 3) for k, v in curve.items()})
keep = [i for i, g in enumerate(X.columns) if g not in IEG]
print("excluded IEG:", sum(g in IEG for g in X.columns), "AUC after exclusion", round(cv_auc(Xv[:, keep], y, b, "enet", "inner", 10).mean(), 3))

# 5. univariate Welch tests (batch-agnostic, as reported) + candidate genes
t, pw = stats.ttest_ind(Xv[y == 1], Xv[y == 0], equal_var=False)
fdr = bh(np.nan_to_num(pw, nan=1))
uv = pd.DataFrame({"gene": X.columns, "t": t, "p": pw, "fdr": fdr}).sort_values("p")
uv.to_csv(RES / "recovery_univariate.csv", index=False)
print("genes FDR<0.05:", int((fdr < 0.05).sum()), "top", uv.head(3)[["gene", "fdr"]].values.tolist())
for g in ["LRRN4CL", "ATF3"]:
    r = pd.read_csv(COH / "GSE328059_expr.csv.gz", index_col=0).loc[g, m["sample"].astype(str)].values if g in X.columns or True else None
    tt = stats.ttest_ind(r[y == 1], r[y == 0], equal_var=False)
    print(g, "Welch p", round(tt.pvalue, 4), "R-NR mean diff", round(r[y == 1].mean() - r[y == 0].mean(), 3))
json.dump({k: float(v) for k, v in out.items()}, open(RES / "recovery_summary.json", "w"))
