#!/usr/bin/env python3
"""Stage 03: leave-one-cohort-out. Gene selection (meta-analysis, top N by |z|) uses TRAINING cohorts only."""
import sys
sys.path.insert(0, "../02_meta_lasso_shap")
import numpy as np, pandas as pd, xgboost as xgb
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from common import *

NS = [5, 10, 20, 50, 100, 200]
ERY = "HBB HBA1 HBA2 HBD HBG1 HBG2 HBE1 HBM HBQ1 ALAS2 SLC4A1 SLC25A37 SNCA BPGM AHSP GYPA GYPB GYPC EPB42 SPTA1 CA1".split()
coh = load_discovery()
names = list(coh)
GV = {c: hedges_g(*[coh[c][0], coh[c][1]]) for c in names}


def fit_score(name, Xtr, ytr, Xte):
    if name == "L1":
        m = LogisticRegression(penalty="l1", solver="liblinear", C=0.1, max_iter=500)
    elif name == "L2":
        m = LogisticRegression(C=0.1, max_iter=1000)
    elif name == "RF":
        m = RandomForestClassifier(300, random_state=SEED, n_jobs=1)
    else:
        m = xgb.XGBClassifier(n_estimators=200, max_depth=3, learning_rate=0.05, random_state=SEED, n_jobs=1, tree_method="hist")
    return m.fit(Xtr, ytr).predict_proba(Xte)[:, 1]


def fold(held, drop=None, strict=None):
    tr = [c for c in names if c != held]
    G = pd.DataFrame({c: GV[c][0] for c in tr})
    V = pd.DataFrame({c: GV[c][1] for c in tr})
    ok = G.notna().sum(axis=1) >= len(tr) - 1
    ok &= G.index.isin(coh[held][0].index)
    meta = dersimonian_laird(G[ok], V[ok])
    if drop:
        meta = meta.drop(index=[g for g in drop if g in meta.index])
    if strict is not None:  # restrict to under-studied genes inside the training fold
        top300 = meta["z"].abs().sort_values(ascending=False).index[:300]
        meta = meta.loc[[g for g in top300 if g in strict]]
    order = meta["z"].abs().sort_values(ascending=False).index
    rows = []
    for N in NS:
        genes = [g for g in order[:N] if all(g in coh[c][0].index for c in tr)]
        Xtr = np.hstack([coh[c][0].loc[genes].values for c in tr]).T
        ytr = np.concatenate([coh[c][1] for c in tr])
        Xte, yte = coh[held][0].loc[genes].values.T, coh[held][1]
        for mname in ["L1", "L2", "RF", "XGB"]:
            rows.append(dict(held_out=held, N=N, model=mname, auc=roc_auc_score(yte, fit_score(mname, Xtr, ytr, Xte)), n_genes=len(genes)))
    return rows


def run(tag, **kw):
    out = Parallel(n_jobs=min(len(names), 8))(delayed(fold)(h, **kw) for h in names)
    df = pd.DataFrame([r for o in out for r in o])
    df.insert(0, "variant", tag)
    return df


res = [run("standard"), run("no_erythrocyte", drop=ERY)]
allr = pd.concat(res)
allr.to_csv(RES / "loco_auc.csv", index=False)
s = allr[allr.variant == "standard"].groupby(["model", "N"]).auc.mean().unstack()
print(s.round(3))
print(allr[(allr.variant == "standard") & (allr.model == "L2") & (allr.N == 20)][["held_out", "auc"]].round(3).to_string(index=False))
print("no-ery L2 N=20 mean", allr[(allr.variant == "no_erythrocyte") & (allr.model == "L2") & (allr.N == 20)].auc.mean().round(3))
