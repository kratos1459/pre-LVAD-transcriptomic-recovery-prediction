#!/usr/bin/env python3
"""Stage 06b: freeze the GSE328059 elastic-net recovery model (in-fold-style batch adjustment fitted on all 84),
then score PRJNA1312133 pre-LVAD hearts ONCE: within-cohort standardise, rescale to derivation moments,
5,000-bootstrap AUC CI and exact permutation test over all C(19,10)=92,378 label assignments.
Usage: python frozen_external_test.py QUANT_DIR   (QUANT_DIR from quantify_prjna1312133.sh; tx2gene from GENCODE v46 gtf names)"""
import sys, itertools, gzip, re
from pathlib import Path
sys.path.insert(0, "../05_recovery"); sys.path.insert(0, "../02_meta_lasso_shap")
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from pipeline import *

Q = Path(sys.argv[1]) / "quant"
# transcript -> gene symbol from the GENCODE fasta headers (ENST|ENSG|...|SYMBOL|)
hdr = {}
for l in gzip.open(Q.parent / "ref" / "gencode.v46.transcripts.fa.gz", "rt"):
    if l[0] == ">":
        f = l[1:].split("|")
        hdr[f[0]] = f[5]
rows = {}
for d in Q.iterdir():
    if re.search(r"_HF_(R|NR)\d+$", d.name):
        t = pd.read_csv(d / "quant.sf", sep="\t", index_col=0)
        t["g"] = t.index.map(hdr)
        rows[d.name] = np.log2(t.groupby("g").TPM.sum() + 1)
E = pd.DataFrame(rows)
yext = np.array([int(re.search(r"_HF_(R|NR)", c).group(1) == "R") for c in E.columns])
print("external n", len(yext), "recovered", yext.sum())

Xd, yd, bd, _ = load_pre()
genes = [g for g in Xd.columns if g in E.index]
print(f"{len(genes)}/{Xd.shape[1]} model genes shared")
Xdv = Xd[genes].values
Xe = E.loc[genes].T.values
Xe = (Xe - Xe.mean(0)) / (Xe.std(0) + 1e-9)             # standardise within external cohort
Xe = Xe * (Xdv.std(0) + 1e-9) + Xdv.mean(0)             # rescale to derivation moments
score = run_fold(Xdv, yd, bd, Xe, np.array(["other"] * len(Xe)), "enet", "inner")  # frozen model, applied once
auc = roc_auc_score(yext, score)
rng = np.random.default_rng(SEED)
boots = []
while len(boots) < 5000:
    i = rng.integers(0, len(yext), len(yext))
    if 0 < yext[i].sum() < len(i):
        boots.append(roc_auc_score(yext[i], score[i]))
ge = tot = 0
for comb in itertools.combinations(range(len(yext)), int(yext.sum())):
    yp = np.zeros(len(yext), int); yp[list(comb)] = 1
    tot += 1; ge += roc_auc_score(yp, score) >= auc
print(f"AUC {auc:.3f}  95% CI {np.percentile(boots, 2.5):.2f}-{np.percentile(boots, 97.5):.2f}  exact permutation p = {ge / tot:.3f} ({tot} assignments)")
