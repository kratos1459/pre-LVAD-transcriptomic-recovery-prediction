"""Shared helpers: load harmonised cohorts, within-cohort z-scoring, Hedges g, DerSimonian-Laird."""
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
COH = ROOT / "results" / "cohorts"
RES = ROOT / "results"
SEED = 20261024
DISCOVERY = ["GSE57338", "GSE21610", "GSE46224", "GSE116250", "GSE141910", "GSE5406", "GSE328059", "GSE226314"]


def load_discovery(extra_ok=True):
    """Return dict cohort -> (genes x samples z-scored within cohort, y (1=HF,0=NF))."""
    meta = pd.read_csv(COH / "all_meta.csv")
    out = {}
    for c in DISCOVERY:
        f = COH / f"{c}_expr.csv.gz"
        if not f.exists():
            continue
        mm = meta[(meta.cohort == c) & meta.group.isin(["HF", "NF"])]
        x = pd.read_csv(f, index_col=0)[mm["sample"].astype(str).values]
        x = x.loc[x.notna().all(axis=1) & (x.std(axis=1) > 0)]
        z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1, ddof=1), axis=0)
        out[c] = (z, (mm["group"].values == "HF").astype(int))
    return out


def hedges_g(z, y):
    a, b = z.loc[:, y == 1], z.loc[:, y == 0]
    n1, n0 = a.shape[1], b.shape[1]
    sp = np.sqrt(((n1 - 1) * a.var(axis=1, ddof=1) + (n0 - 1) * b.var(axis=1, ddof=1)) / (n1 + n0 - 2))
    d = (a.mean(axis=1) - b.mean(axis=1)) / sp
    J = 1 - 3 / (4 * (n1 + n0) - 9)
    g = d * J
    v = (n1 + n0) / (n1 * n0) + g ** 2 / (2 * (n1 + n0))
    return g, v


def dersimonian_laird(G, V):
    """G,V: genes x cohorts DataFrames (NaN allowed). Returns DataFrame of pooled stats."""
    G, V = G.copy(), V.copy()
    w = 1 / V
    k = G.notna().sum(axis=1)
    sw = w.sum(axis=1)
    gf = (w * G).sum(axis=1) / sw
    Q = (w * (G.sub(gf, axis=0)) ** 2).sum(axis=1)
    C = sw - (w ** 2).sum(axis=1) / sw
    tau2 = ((Q - (k - 1)) / C).clip(lower=0)
    ws = 1 / (V.add(tau2, axis=0))
    sws = ws.sum(axis=1)
    g = (ws * G).sum(axis=1) / sws
    se = np.sqrt(1 / sws)
    z = g / se
    p = 2 * stats.norm.sf(np.abs(z))
    I2 = ((Q - (k - 1)) / Q).clip(lower=0) * 100
    return pd.DataFrame({"k": k, "g": g, "se": se, "z": z, "p": p, "tau2": tau2, "I2": I2})


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.minimum(q, 1)
