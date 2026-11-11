#!/usr/bin/env python3
"""Stage 08 (UNTESTED TEMPLATE): secondary single-nucleus / perturbation cohorts with patient-level exact tests.
GSE285315 (paediatric DCM, 3 recovered vs 3 non-recovered), GSE301840 (6 paired pre/post LVAD), GSE253535 (1 paired patient),
GSE327946 (iPSC-CM LRRN4CL vs GFP). Fill CONFIG with the processed-file paths downloaded from GEO, then run.
Tests are patient-level (never cell-level): exact label permutation (unpaired) or exact sign-flip (paired)."""
import itertools
import numpy as np, pandas as pd, scanpy as sc

CONFIG = {  # accession: dict(path=<h5ad/10x mtx dir>, patient_col, group_col, time_col)
    "GSE285315": dict(path="", patient_col="patient", group_col="recovery", time_col=None),
    "GSE301840": dict(path="", patient_col="patient", group_col=None, time_col="timepoint"),
    "GSE253535": dict(path="", patient_col="patient", group_col=None, time_col="timepoint"),
    "GSE327946": dict(path="", patient_col="sample", group_col="treatment", time_col=None),
}


def pseudobulk(adata, patient_col):
    sc.pp.normalize_total(adata, target_sum=1e6); sc.pp.log1p(adata, base=2)
    X = pd.DataFrame(adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X, index=adata.obs_names, columns=adata.var_names)
    return X.groupby(adata.obs[patient_col].values).mean()


def perm_p(score, y):
    obs = score[y == 1].mean() - score[y == 0].mean(); k = int(y.sum()); ge = tot = 0
    for comb in itertools.combinations(range(len(y)), k):
        yp = np.zeros(len(y), int); yp[list(comb)] = 1
        ge += abs(score[yp == 1].mean() - score[yp == 0].mean()) >= abs(obs) - 1e-12; tot += 1
    return obs, ge / tot


def signflip_p(delta):
    obs = delta.mean(); n = len(delta); ge = 0
    for s in itertools.product([-1, 1], repeat=n):
        ge += abs((np.array(s) * delta).mean()) >= abs(obs) - 1e-12
    return obs, ge / 2 ** n


if __name__ == "__main__":
    cons = pd.read_csv("../results/meta_consensus.csv")
    cons = cons[(cons.fdr < 1e-6) & (cons.g.abs() > .8)]
    sign = cons.set_index("gene").g.pipe(np.sign)
    for acc, cfg in CONFIG.items():
        if not cfg["path"]:
            print(acc, "- set CONFIG path first"); continue
        ad = sc.read(cfg["path"]); pb = pseudobulk(ad, cfg["patient_col"])
        g = [x for x in sign.index if x in pb.columns]
        score = (pb[g].sub(pb[g].mean()).div(pb[g].std() + 1e-9) * sign[g]).mean(axis=1)  # consensus-signature score per patient
        print(acc, "patients", len(pb), "signature score computed for", len(g), "genes; add contrast (group/time) per cohort design")
