#!/usr/bin/env python3
"""Stage 02: random-effects meta-analysis, bootstrap LASSO stability, XGBoost SHAP, consensus rank."""
import numpy as np, pandas as pd, xgboost as xgb, shap
from joblib import Parallel, delayed
from sklearn.feature_selection import f_classif
from sklearn.linear_model import LogisticRegression
from common import *

N_BOOT = 500
coh = load_discovery()
print({c: (z.shape[1], int(y.sum())) for c, (z, y) in coh.items()}, "total", sum(z.shape[1] for z, _ in coh.values()))
G, V = {}, {}
for c, (z, y) in coh.items():
    G[c], V[c] = hedges_g(z, y)
G, V = pd.DataFrame(G), pd.DataFrame(V)
keep = G.notna().sum(axis=1) >= 6
meta = dersimonian_laird(G[keep], V[keep])
meta["fdr"] = bh(meta["p"].values)
meta["n_concordant"] = (np.sign(G[keep]).eq(np.sign(meta["g"]), axis=0) & G[keep].notna()).sum(axis=1)
print("genes tested", len(meta))

shared = G.dropna().index
print("shared by all cohorts", len(shared))
names = list(coh)
X = np.hstack([coh[c][0].loc[shared].values for c in names]).T
y = np.concatenate([coh[c][1] for c in names])
cid = np.concatenate([[i] * coh[c][0].shape[1] for i, c in enumerate(names)])


def one(b):
    rng = np.random.default_rng(SEED + b)
    idx = np.concatenate([rng.choice(np.where(cid == i)[0], (cid == i).sum(), replace=True) for i in range(len(names))])
    Xb, yb = X[idx], y[idx]
    F, _ = f_classif(Xb, yb)
    top = np.argsort(np.nan_to_num(F))[::-1][:500]
    m = LogisticRegression(penalty="l1", solver="liblinear", C=0.1, max_iter=500).fit(Xb[:, top], yb)
    sel = np.zeros(X.shape[1], bool)
    sel[top[m.coef_[0] != 0]] = True
    return sel


sel = np.array(Parallel(n_jobs=-1)(delayed(one)(b) for b in range(N_BOOT)))
lasso = pd.Series(sel.mean(axis=0), index=shared, name="lasso_freq")

# XGBoost + SHAP (restricted to the 2,000 genes with the largest |meta z| among the shared set for tractability)
cand = meta.loc[shared].reindex(shared)["z"].abs().sort_values(ascending=False).index[:2000]
ci = [list(shared).index(g) for g in cand]
dtr = xgb.XGBClassifier(n_estimators=300, max_depth=3, learning_rate=0.05, subsample=0.8, colsample_bytree=0.5, random_state=SEED, n_jobs=16, tree_method="hist")
dtr.fit(X[:, ci], y)
sv = shap.TreeExplainer(dtr).shap_values(X[:, ci])
shap_s = pd.Series(np.abs(sv).mean(axis=0), index=cand, name="mean_abs_shap")

res = meta.join(lasso).join(shap_s)
pr = lambda s: s.rank(pct=True)
res["consensus_rank_score"] = pd.concat([pr(res["z"].abs()), pr(res["lasso_freq"]), pr(res["mean_abs_shap"])], axis=1).mean(axis=1, skipna=False)  # NaN unless gene shared by all cohorts
res = res.sort_values("consensus_rank_score", ascending=False, na_position="last")
res.insert(0, "gene", res.index)
res.to_csv(RES / "meta_consensus.csv", index=False)
G.to_csv(RES / "hedges_g_per_cohort.csv")
sig = (res.fdr < 1e-6) & (res.g.abs() > 0.8)
print("consensus HF genes:", int(sig.sum()), "up", int((sig & (res.g > 0)).sum()), "down", int((sig & (res.g < 0)).sum()))
print(res.head(12)[["gene", "g", "fdr", "lasso_freq", "mean_abs_shap"]])
