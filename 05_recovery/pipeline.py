"""Leakage-controlled recovery pipeline: everything that touches outcome is fitted inside the training fold."""
import sys
sys.path.insert(0, "../02_meta_lasso_shap")
import numpy as np, pandas as pd, xgboost as xgb
from sklearn.feature_selection import f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from common import COH, RES, SEED

TOPVAR = 6000
IEG = ("ATF3 JUNB JUN JUND FOS FOSB FOSL1 FOSL2 EGR1 EGR2 EGR3 EGR4 KLF2 KLF4 KLF6 NR4A1 NR4A2 NR4A3 DUSP1 DUSP2 DUSP5 DUSP6 IER2 IER3 IER5 "
       "MYC SOCS3 BTG2 ZFP36 CYR61 CTGF NFKBIA HSPA1A HSPA1B HSPA6 GADD45B CEBPB").split()  # 37 stress-response genes


def load_pre():
    meta = pd.read_csv(COH / "all_meta.csv")
    m = meta[(meta.cohort == "GSE328059") & (meta.timepoint == "pre")].copy()
    x = pd.read_csv(COH / "GSE328059_expr.csv.gz", index_col=0)[m["sample"].astype(str)]
    y = (m["recovery"].values == "R").astype(int)
    b = m["batch"].astype(float).astype(int).astype(str).values
    b = np.where(np.isin(b, ["10931", "11943"]), b, "other")
    keep = x.var(axis=1).sort_values(ascending=False).index[:TOPVAR]
    return x.loc[keep].T, y, b, m


def design_batch(b):
    return np.column_stack([(b == "10931"), (b == "11943")]).astype(float)  # 'other' = reference


def adjust_batch(Xtr, ytr, btr, Xte, bte):
    """OLS of expression on batch dummies + outcome in TRAIN only; subtract batch terms from train and test."""
    Dtr = np.column_stack([np.ones(len(ytr)), design_batch(btr), ytr])
    beta, *_ = np.linalg.lstsq(Dtr, Xtr, rcond=None)
    bb = beta[1:3]
    return Xtr - design_batch(btr) @ bb, Xte - design_batch(bte) @ bb


def location_scale(Xtr, btr, Xte, bte):
    """ComBat-like per-batch centring/scaling using training statistics (outcome-agnostic)."""
    Xtr2, Xte2 = Xtr.copy(), Xte.copy()
    mu0, sd0 = Xtr.mean(0), Xtr.std(0) + 1e-9
    for k in np.unique(btr):
        i, j = btr == k, bte == k
        mu, sd = Xtr[i].mean(0), Xtr[i].std(0) + 1e-9
        Xtr2[i] = (Xtr[i] - mu) / sd * sd0 + mu0
        if j.any():
            Xte2[j] = (Xte[j] - mu) / sd * sd0 + mu0
    return Xtr2, Xte2


def make_model(name, seed):
    if name == "lasso":
        return LogisticRegression(penalty="l1", solver="liblinear", C=0.1, max_iter=1000)
    if name == "enet":
        return LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, C=0.1, max_iter=3000, random_state=seed)
    if name == "ridge":
        return LogisticRegression(C=0.005, max_iter=2000)
    if name == "svm":
        return LinearSVC(C=0.002, max_iter=20000)
    if name == "rf":
        return RandomForestClassifier(300, random_state=seed, n_jobs=1)
    return xgb.XGBClassifier(n_estimators=150, max_depth=2, learning_rate=0.05, subsample=0.8, colsample_bytree=0.5, random_state=seed, n_jobs=1, tree_method="hist")


def run_fold(Xtr, ytr, btr, Xte, bte, model, adjust="inner", nsel=200, seed=0):
    if adjust == "inner":
        Xtr, Xte = adjust_batch(Xtr, ytr, btr, Xte, bte)
    elif adjust == "locscale":
        Xtr, Xte = location_scale(Xtr, btr, Xte, bte)
    if adjust == "batch_only":
        m = LogisticRegression(C=1.0).fit(design_batch(btr), ytr)
        return m.predict_proba(design_batch(bte))[:, 1]
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
    Ztr, Zte = (Xtr - mu) / sd, (Xte - mu) / sd
    F, _ = f_classif(Ztr, ytr)
    top = np.argsort(np.nan_to_num(F))[::-1][:nsel]
    m = make_model(model, seed).fit(Ztr[:, top], ytr)
    return m.decision_function(Zte[:, top]) if hasattr(m, "decision_function") else m.predict_proba(Zte[:, top])[:, 1]


def cv_auc(X, y, b, model, adjust="inner", repeats=10, seed=SEED, nsel=200, return_oof=False):
    aucs, oofs = [], []
    for r in range(repeats):
        oof = np.zeros(len(y))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed + r).split(X, y):
            oof[te] = run_fold(X[tr], y[tr], b[tr], X[te], b[te], model, adjust, nsel, seed + r)
        aucs.append(roc_auc_score(y, oof))
        oofs.append(oof)
    return (np.array(aucs), np.array(oofs)) if return_oof else np.array(aucs)
