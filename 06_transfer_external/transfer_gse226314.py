#!/usr/bin/env python3
"""Stage 06d (UNTESTED until pseudobulk exists): cross-technology transfer + LVAD reversal in GSE226314.
Inputs from pseudobulk_gse226314.R: results/cohorts/GSE226314_pseudobulk_expr.csv, ..._by_celltype.csv, GSE226314_sample_condition.csv
(1) frozen GSE328059 elastic-net applied to within-cohort z-scored pre-LVAD pseudobulk (R vs NR), all nuclei + each cell type,
    exact two-sided permutation over all label assignments, BH across tests.
(2) reversal: paired post-minus-pre change of each consensus gene, signed by minus the HF effect."""
import sys, itertools, re
sys.path.insert(0, "../05_recovery"); sys.path.insert(0, "../02_meta_lasso_shap")
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from pipeline import *
from common import bh

cond = pd.read_csv(COH / "GSE226314_sample_condition.csv")
cond.columns = ["sample", "condition"]
c = cond.set_index("sample").condition.str.lower()
pre_R = [s for s in c.index if "responder pre" in c[s] and "non" not in c[s]]
pre_NR = [s for s in c.index if "non-responder pre" in c[s]]
print("pre-LVAD patients: R", len(pre_R), "NR", len(pre_NR))

Xd, yd, bd, _ = load_pre()


def frozen_score(P):  # P: genes x patients (log2 CPM)
    genes = [g for g in Xd.columns if g in P.index]
    Z = P.loc[genes].T.values
    Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-9)                      # within-cohort z-score
    Z = Z * Xd[genes].values.std(0) + Xd[genes].values.mean(0)   # derivation moments
    return run_fold(Xd[genes].values, yd, bd, Z, np.array(["other"] * len(Z)), "enet", "inner")


def exact_p(score, y):
    obs = roc_auc_score(y, score); k = int(y.sum()); ge = tot = 0
    for comb in itertools.combinations(range(len(y)), k):
        yp = np.zeros(len(y), int); yp[list(comb)] = 1
        ge_ = abs(roc_auc_score(yp, score) - .5) >= abs(obs - .5) - 1e-12; ge += ge_; tot += 1
    return obs, ge / tot


rows = []
all_pb = pd.read_csv(RES / "cohorts" / "GSE226314_pseudobulk_expr.csv", index_col=0)
ct = pd.read_csv(RES / "cohorts" / "GSE226314_pseudobulk_by_celltype.csv", index_col=0)
pats = pre_R + pre_NR; y = np.array([1] * len(pre_R) + [0] * len(pre_NR))
units = {"all_nuclei": all_pb[pats]}
for cell in sorted({re.sub(r".*__", "", x) for x in ct.columns}):
    cols = [f"{p}__{cell}" for p in pats]
    if all(x in ct.columns for x in cols):
        u = ct[cols]; u.columns = pats; units[cell] = u
for name, P in units.items():
    a, p = exact_p(frozen_score(P), y); rows.append(dict(unit=name, auc=a, p=p))
res = pd.DataFrame(rows); res["fdr"] = bh(res.p.values)
res.to_csv(RES / "transfer_gse226314.csv", index=False); print(res.round(3))
