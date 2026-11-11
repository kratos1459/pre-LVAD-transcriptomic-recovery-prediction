#!/usr/bin/env python3
"""Stage 07: pre-LVAD molecular subtypes. All 218 GSE328059 pre-LVAD hearts (84 labelled + 134 unlabelled):
3,000 variable genes -> 30 PCs -> Harmony (batch) -> consensus K-means (K=2..6, 300 subsamples) -> K by minimum PAC."""
import sys
sys.path.insert(0, "../02_meta_lasso_shap")
import numpy as np, pandas as pd, harmonypy as hp
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from common import COH, RES, SEED

meta = pd.read_csv(COH / "all_meta.csv")
m = meta[(meta.cohort == "GSE328059") & meta.timepoint.isin(["pre", "unlabelled"])].reset_index(drop=True)
x = pd.read_csv(COH / "GSE328059_expr.csv.gz", index_col=0)[m["sample"].astype(str)]
print("pre-LVAD hearts", x.shape[1], "labelled", int(m.recovery.notna().sum()))
x = x.loc[x.var(axis=1).sort_values(ascending=False).index[:3000]]
# derivation and validation matrices were rlog-transformed separately -> centre/scale genes within each part first
part = m["cohort_part"].values
z = x.T.copy()
for p_ in np.unique(part):
    i = part == p_
    z.loc[i] = (z.loc[i] - z.loc[i].mean()) / z.loc[i].std()
z = z.values
pcs = PCA(30, random_state=SEED).fit_transform(z)
bt = m["batch"].astype(float).astype(int).astype(str)
hb_ = bt.where(bt.isin(["10931", "11943", "15792"]), "other")
hm = hp.run_harmony(pcs, pd.DataFrame({"batch": hb_.values}), "batch", random_state=SEED, verbose=False)
Z = np.asarray(hm.Z_corr)
Z = Z.T if Z.shape[0] == 30 else Z


def consensus(Z, k, n_iter=300, frac=0.8):
    n = len(Z)
    co, tot = np.zeros((n, n)), np.zeros((n, n))
    rng = np.random.default_rng(SEED + k)
    for _ in range(n_iter):
        idx = rng.choice(n, int(frac * n), replace=False)
        lab = KMeans(k, n_init=3, random_state=int(rng.integers(1e9))).fit_predict(Z[idx])
        same = (lab[:, None] == lab[None, :])
        co[np.ix_(idx, idx)] += same
        tot[np.ix_(idx, idx)] += 1
    return np.divide(co, tot, out=np.zeros_like(co), where=tot > 0)


rows, M = [], {}
for k in range(2, 7):
    C = consensus(Z, k)
    v = C[np.triu_indices_from(C, 1)]
    pac = ((v > 0.1) & (v < 0.9)).mean()
    M[k] = C
    rows.append(dict(K=k, PAC=pac))
pac = pd.DataFrame(rows)
print(pac.round(3).to_string(index=False))
K = int(pac.sort_values("PAC").K.iloc[0])
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
D = 1 - M[K]
np.fill_diagonal(D, 0)
lab = fcluster(linkage(squareform(D, checks=False), "average"), K, "maxclust")
m["subtype"] = ["S%d" % l for l in lab]
print("selected K =", K, m.subtype.value_counts().to_dict())
lb = m[m.recovery.notna()]
ct = pd.crosstab(lb.subtype, lb.recovery)
print(ct)
print("Fisher/chi2 p (subtype x recovery):", round(stats.chi2_contingency(ct)[1], 4) if ct.shape[0] > 1 else "n/a")
print("subtype x batch chi2 p:", stats.chi2_contingency(pd.crosstab(m.subtype, bt.where(bt.isin(["10931", "11943", "15792"]), "other")))[1])
m[["sample", "batch", "recovery", "subtype"]].to_csv(RES / "subtypes.csv", index=False)
pac.to_csv(RES / "subtype_pac.csv", index=False)
