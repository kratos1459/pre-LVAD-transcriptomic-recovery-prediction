#!/usr/bin/env python3
"""Stage 13: head-to-head benchmark of published and candidate recovery predictors.

Every predictor goes through the identical leakage-controlled frame used for the main model:
 - fixed scores (single genes, signed gene-set scores) are evaluated directly, with and without
   batch centring, and once in the independent PRJNA1312133 cohort;
 - prior gene lists also restrict feature selection, so each list gets its own fitted elastic-net
   model under 10 x 5-fold cross-validation and on each held-out sequencing batch;
 - size-matched random gene sets give the null those fitted models must beat.

usage: python benchmark_signatures.py [../data/prjna]
"""
import sys, json, gzip, re
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "05_recovery"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "02_meta_lasso_shap"))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
from sklearn.feature_selection import f_classif
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from pipeline import load_pre, adjust_batch, make_model
from common import COH, RES, SEED, bh

FIG = Path(__file__).resolve().parent.parent / "figures"
GS = json.load(open(Path(__file__).resolve().parent.parent / "11_deconvolution_pathways" / "genesets.json"))
plt.rcParams.update({"font.family": "serif", "font.size": 10, "axes.linewidth": 1.1})
REP, NRAND = 10, 20

X, y, b, m = load_pre()
genes = list(X.columns); Xv = X.values
full = pd.read_csv(COH / "GSE328059_expr.csv.gz", index_col=0)[m["sample"].astype(str)]
cons = pd.read_csv(RES / "meta_consensus.csv")
top150 = pd.read_csv(RES / "top150_consensus.csv")
uv = pd.read_csv(RES / "recovery_univariate.csv")


def z(df):
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1) + 1e-9, axis=0)


Z = z(full)


def signed_score(up, dn, Zm=None):
    Zm = Z if Zm is None else Zm
    u = [g for g in up if g in Zm.index]; d = [g for g in dn if g in Zm.index]
    s = np.zeros(Zm.shape[1])
    if u: s += Zm.loc[u].mean(0).values
    if d: s -= Zm.loc[d].mean(0).values
    return s, len(u) + len(d)


def adj_batch_1d(v, bb):
    v = np.asarray(v, float).copy()
    for k in np.unique(bb): v[bb == k] -= v[bb == k].mean()
    return v


def auc2(s, yy):
    a = roc_auc_score(yy, s)
    return max(a, 1 - a)


hf_up = cons[(cons.fdr < 1e-6) & (cons.g > 0.8)].gene.dropna().tolist()
hf_dn = cons[(cons.fdr < 1e-6) & (cons.g < -0.8)].gene.dropna().tolist()
hall = {k.split("|")[-1]: v for k, v in GS.items() if k.startswith("MSigDB")}
trr = {k.split("|")[-1]: v for k, v in GS.items() if k.startswith("TRRUST")}
t150_up = top150[top150.g > 0].gene.dropna().tolist(); t150_dn = top150[top150.g < 0].gene.dropna().tolist()

fixed = {}
for g in ["LRRN4CL", "NPPB", "ATF3", "CAMK2D", "PCDH9", "MYL7", "NPPA", "MYH7"]:
    if g in Z.index: fixed[f"gene: {g}"] = (Z.loc[g].values, 1)
fixed["HF consensus signature (607)"] = signed_score(hf_up, hf_dn)
fixed["HF consensus top 150"] = signed_score(t150_up, t150_dn)
for s in ["Epithelial Mesenchymal Transition", "Myogenesis", "Inflammatory Response", "TGF-beta Signaling",
          "Oxidative Phosphorylation", "p53 Pathway", "Hypoxia", "Apoptosis"]:
    if s in hall: fixed[f"Hallmark: {s}"] = signed_score(hall[s], [])
for t in ["NFE2L2 human", "ATF3 human", "NFKB1 human", "TP53 human"]:
    if t in trr: fixed[f"TRRUST: {t}"] = signed_score(trr[t], [])

rows = []
for name, (s, n) in fixed.items():
    rows.append(dict(predictor=name, n_genes=n, auc_raw=auc2(s, y), auc_batch_centred=auc2(adj_batch_1d(s, b), y),
                     p_unadjusted=stats.mannwhitneyu(s[y == 1], s[y == 0], alternative="two-sided").pvalue))
FX = pd.DataFrame(rows).sort_values("auc_batch_centred", ascending=False)
FX["fdr_internal"] = bh(FX.p_unadjusted.values)
print(FX.round(3).to_string(index=False))


def cv_restricted(idx, repeats=REP, nsel=200):
    sub = Xv[:, idx]; aucs = []
    for r in range(repeats):
        oof = np.zeros(len(y))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEED + r).split(sub, y):
            Xtr, Xte = adjust_batch(sub[tr], y[tr], b[tr], sub[te], b[te])
            mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
            Ztr, Zte = (Xtr - mu) / sd, (Xte - mu) / sd
            k = min(nsel, Ztr.shape[1])
            F_, _ = f_classif(Ztr, y[tr]); top = np.argsort(np.nan_to_num(F_))[::-1][:k]
            oof[te] = make_model("enet", SEED + r).fit(Ztr[:, top], y[tr]).decision_function(Zte[:, top])
        aucs.append(roc_auc_score(y, oof))
    return np.array(aucs)


def heldout(idx, nsel=200):
    out = {}; sub = Xv[:, idx]
    for held in ["11943", "10931"]:
        te = b == held
        Xtr, Xte = adjust_batch(sub[~te], y[~te], b[~te], sub[te], b[te])
        mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
        Ztr, Zte = (Xtr - mu) / sd, (Xte - mu) / sd
        k = min(nsel, Ztr.shape[1])
        F_, _ = f_classif(Ztr, y[~te]); top = np.argsort(np.nan_to_num(F_))[::-1][:k]
        out[held] = float(roc_auc_score(y[te], make_model("enet", SEED).fit(Ztr[:, top], y[~te]).decision_function(Zte[:, top])))
    return out


lists = {"all 6,000 variable genes": genes,
         "HF consensus (607)": hf_up + hf_dn,
         "HF consensus top 150": top150.gene.dropna().tolist(),
         "Hallmark union": sorted({g for s in hall.values() for g in s}),
         "TRRUST union": sorted({g for s in trr.values() for g in s})}
gi = {g: i for i, g in enumerate(genes)}
rng = np.random.default_rng(SEED)
res = []
for name, gl in lists.items():
    idx = [gi[g] for g in dict.fromkeys(gl) if g in gi]
    if len(idx) < 10: continue
    a = cv_restricted(idx); h = heldout(idx)
    r = dict(gene_list=name, n_in_6000=len(idx), cv_auc=a.mean(), cv_sd=a.std(),
             heldout_11943=h["11943"], heldout_10931=h["10931"], random_mean=np.nan, random_p95=np.nan)
    if len(idx) < len(genes) * 0.9:
        draws = [cv_restricted(rng.choice(len(genes), len(idx), replace=False), repeats=1)[0] for _ in range(NRAND)]
        r["random_mean"] = float(np.mean(draws)); r["random_p95"] = float(np.percentile(draws, 95))
    res.append(r)
FT = pd.DataFrame(res)
FT.to_csv(RES / "benchmark_fitted_models.csv", index=False)
print(FT.round(3).to_string(index=False))

ext = {}; FXE = FX
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
    E = pd.DataFrame(rw); E = E.loc[E.std(axis=1) > 0]; Ze = z(E)
    yext = np.array([int(re.search(r"_HF_(R|NR)", c).group(1) == "R") for c in E.columns])
    erows = []
    for name in FX.predictor:
        if name.startswith("gene: "):
            g = name.split(": ")[1]
            if g not in Ze.index: continue
            se = Ze.loc[g].values
        elif "607" in name: se, _ = signed_score(hf_up, hf_dn, Ze)
        elif "top 150" in name: se, _ = signed_score(t150_up, t150_dn, Ze)
        elif name.startswith("Hallmark: "): se, _ = signed_score(hall[name.split(": ")[1]], [], Ze)
        else: se, _ = signed_score(trr[name.split(": ")[1]], [], Ze)
        erows.append(dict(predictor=name, auc_external=auc2(se, yext),
                          p_external=stats.mannwhitneyu(se[yext == 1], se[yext == 0], alternative="two-sided", method="exact").pvalue))
    EX = pd.DataFrame(erows); EX["fdr_external"] = bh(EX.p_external.values)
    FXE = FX.merge(EX, on="predictor", how="left")
    both = FX.merge(EX, on="predictor")
    ext = dict(n=int(len(yext)), n_tested=int(len(EX)), n_p05=int((EX.p_external < 0.05).sum()),
               n_fdr05_external=int((EX.fdr_external < 0.05).sum()),
               n_significant_in_both=int(((both.p_unadjusted < 0.05) & (both.p_external < 0.05)).sum()),
               best=EX.sort_values("auc_external").iloc[-1].to_dict())
    print(EX.sort_values("auc_external", ascending=False).round(3).to_string(index=False))
FXE.to_csv(RES / "benchmark_fixed_scores.csv", index=False)

summary = dict(n_fixed_predictors=int(len(FX)), n_fixed_fdr05_internal=int((FX.fdr_internal < 0.05).sum()),
               best_fixed_internal=FX.iloc[0][["predictor", "auc_batch_centred"]].to_dict(),
               n_fixed_p05_internal=int((FX.p_unadjusted < 0.05).sum()),
               fitted=FT.round(3).to_dict("records"), external=ext)
json.dump(summary, open(RES / "benchmark_summary.json", "w"), indent=1, default=str)
print(json.dumps(summary, indent=1, default=str))

fig, ax = plt.subplots(1, 2, figsize=(10.2, 4.4), gridspec_kw={"width_ratios": [1.25, 1]})
d = FXE.sort_values("auc_batch_centred")
ax[0].barh(range(len(d)), d.auc_batch_centred, color="#1f4e9e", edgecolor="k", lw=.5, label="84 labelled hearts")
if "auc_external" in d.columns:
    ax[0].plot(d.auc_external, range(len(d)), "o", color="#e67e22", ms=5, label="independent cohort (n = 19)")
ax[0].axvline(0.5, color="k", ls="--", lw=1); ax[0].legend(fontsize=7, loc="lower right")
ax[0].set_yticks(range(len(d))); ax[0].set_yticklabels(d.predictor, fontsize=6.5)
ax[0].set_xlabel("AUC for recovery"); ax[0].set_xlim(0.4, 0.95); ax[0].set_title("Fixed published and candidate scores", fontsize=10)
i = np.arange(len(FT)); w = 0.26
ax[1].bar(i - w, FT.cv_auc, w, yerr=FT.cv_sd, capsize=3, color="#1f4e9e", edgecolor="k", lw=.5, label="cross-validated")
ax[1].bar(i, FT.heldout_11943, w, color="#2e86c1", edgecolor="k", lw=.5, label="held-out 11943")
ax[1].bar(i + w, FT.heldout_10931, w, color="#7fb3d5", edgecolor="k", lw=.5, label="held-out 10931")
ax[1].plot(i, FT.random_p95, "k_", ms=14, label="random sets, 95th pct")
ax[1].axhline(0.5, color="k", ls="--", lw=1); ax[1].set_ylim(0.3, 0.95)
ax[1].set_xticks(i); ax[1].set_xticklabels([t.replace(" ", "\n", 1) for t in FT.gene_list], fontsize=6.5)
ax[1].set_ylabel("AUC"); ax[1].legend(fontsize=6.5); ax[1].set_title("Models restricted to prior gene lists", fontsize=10)
for a, l in zip(ax, "ab"):
    a.text(-0.1, 1.04, f"({l})", transform=a.transAxes, fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG / "FigS19_signature_benchmark.png", dpi=600)
fig.savefig(FIG / "FigS19_signature_benchmark.tif", dpi=600, pil_kwargs={"compression": "tiff_lzw", "dpi": (600, 600)})
print("wrote", FIG / "FigS19_signature_benchmark.tif")
