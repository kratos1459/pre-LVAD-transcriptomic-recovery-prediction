#!/usr/bin/env python3
"""Stage 15: robustness diagnostics.

(1) Heterogeneity of the consensus HF signature (I-squared, tau-squared), leave-one-cohort-out
    influence, and microarray versus RNA-seq concordance.
(2) How much of each principal component of the 84 pre-LVAD hearts is explained by sequencing
    batch and by outcome.
(3) Sex inferred from X/Y-linked transcripts: association with recovery and with batch, and the
    effect of adding sex to the model.
(4) Comparability of the 84 labelled hearts with the 134 withheld-label validation hearts.

usage: python signature_and_batch_diagnostics.py
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "05_recovery"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "02_meta_lasso_shap"))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.feature_selection import f_classif
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from pipeline import load_pre, adjust_batch, make_model
from common import COH, RES, SEED

FIG = Path(__file__).resolve().parent.parent / "figures"
plt.rcParams.update({"font.family": "serif", "font.size": 10, "axes.linewidth": 1.1})
ARRAY = ["GSE5406", "GSE21610", "GSE57338"]
S = {}

# ------------------------------------------------------------------ (1) signature heterogeneity
cons = pd.read_csv(RES / "meta_consensus.csv")
sig = cons[(cons.fdr < 1e-6) & (cons.g.abs() > 0.8)]
S["signature"] = dict(n=int(len(sig)), median_I2=float(sig.I2.median()),
                      pct_I2_gt50=float((sig.I2 > 50).mean() * 100), pct_I2_gt75=float((sig.I2 > 75).mean() * 100),
                      median_tau2=float(sig.tau2.median()), median_k=float(sig.k.median()))
per = pd.read_csv(RES / "hedges_g_per_cohort.csv", index_col=0)
per = per.loc[[g for g in sig.gene if g in per.index]]
cohorts = [c for c in per.columns]
loco = {}
for c in cohorts:
    sub = per.drop(columns=[c])
    mg = sub.mean(axis=1)
    keep = ((mg.abs() > 0.8) & (np.sign(sub).apply(lambda r: abs(r.dropna().sum()) == r.notna().sum(), axis=1)))
    loco[c] = int(keep.sum())
S["signature"]["loco_retained"] = loco
S["signature"]["loco_min_pct"] = round(100 * min(loco.values()) / len(per), 1)
arr = [c for c in cohorts if c in ARRAY]; rna = [c for c in cohorts if c not in ARRAY]
allg = pd.read_csv(RES / "hedges_g_per_cohort.csv", index_col=0)
ma, mr = allg[arr].mean(axis=1), allg[rna].mean(axis=1)
ok = ma.notna() & mr.notna()
S["signature"]["platform_spearman_all_genes"] = float(stats.spearmanr(ma[ok], mr[ok]).statistic)
s_ok = ok & ok.index.isin(sig.gene)
S["signature"]["platform_spearman_signature"] = float(stats.spearmanr(ma[s_ok], mr[s_ok]).statistic)
S["signature"]["platform_sign_agreement_signature"] = float((np.sign(ma[s_ok]) == np.sign(mr[s_ok])).mean())
print("signature:", json.dumps(S["signature"], indent=1))

# ------------------------------------------------------------------ (2) variance by batch and outcome
X, y, b, m = load_pre(); Xv = X.values
Z = (Xv - Xv.mean(0)) / (Xv.std(0) + 1e-9)
pc = PCA(10, random_state=SEED).fit(Z)
P = pc.transform(Z)
rows = []
for i in range(5):
    v = P[:, i]
    groups = [v[b == k] for k in np.unique(b) if (b == k).sum() > 1]
    fb = stats.f_oneway(*groups)
    ss_tot = ((v - v.mean()) ** 2).sum()
    ss_b = sum(len(g) * (g.mean() - v.mean()) ** 2 for g in groups)
    ss_y = sum((y == k).sum() * (v[y == k].mean() - v.mean()) ** 2 for k in [0, 1])
    rows.append(dict(PC=i + 1, variance_pct=100 * pc.explained_variance_ratio_[i],
                     batch_r2=100 * ss_b / ss_tot, outcome_r2=100 * ss_y / ss_tot,
                     batch_p=fb.pvalue, outcome_p=stats.ttest_ind(v[y == 1], v[y == 0], equal_var=False).pvalue))
V = pd.DataFrame(rows); V.to_csv(RES / "pc_variance_batch_outcome.csv", index=False)
S["pcs"] = V.round(3).to_dict("records")
print(V.round(2).to_string(index=False))

# ------------------------------------------------------------------ (3) inferred sex
full = pd.read_csv(COH / "GSE328059_expr.csv.gz", index_col=0)
lab = m["sample"].astype(str).values
ygenes = [g for g in ["RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY"] if g in full.index]
sexdf = None
if "XIST" in full.index and ygenes:
    xi = full.loc["XIST", lab].values; yy = full.loc[ygenes, lab].mean(0).values
    male = (yy - xi) > np.median(yy - xi)
    # k-means-free split: use the gap in the bimodal difference
    diff = yy - xi
    thr = (diff.max() + diff.min()) / 2 if diff.max() - diff.min() > 1 else np.median(diff)
    male = diff > thr
    sexdf = pd.DataFrame({"sample": lab, "XIST": xi, "Y_mean": yy, "inferred_male": male, "recovery": y, "batch": b})
    sexdf.to_csv(RES / "inferred_sex.csv", index=False)
    tab = pd.crosstab(male, y)
    fp = stats.fisher_exact(tab.values)[1] if tab.shape == (2, 2) else np.nan
    S["sex"] = dict(n_male=int(male.sum()), n_female=int((~male).sum()),
                    responders_male=int(y[male].sum()), responders_female=int(y[~male].sum()),
                    fisher_p_recovery=float(fp),
                    chi2_p_batch=float(stats.chi2_contingency(pd.crosstab(male, b).values)[1]))
    print("sex:", S["sex"])

# ------------------------------------------------------------------ model with and without sex
def cv(Xa, extra=None, repeats=10):
    a = []
    for r in range(repeats):
        oof = np.zeros(len(y))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEED + r).split(Xa, y):
            Xtr, Xte = adjust_batch(Xa[tr], y[tr], b[tr], Xa[te], b[te])
            if extra is not None:
                Xtr = np.column_stack([Xtr, extra[tr]]); Xte = np.column_stack([Xte, extra[te]])
            mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
            Ztr, Zte = (Xtr - mu) / sd, (Xte - mu) / sd
            F_, _ = f_classif(Ztr, y[tr]); top = np.argsort(np.nan_to_num(F_))[::-1][:200]
            oof[te] = make_model("enet", SEED + r).fit(Ztr[:, top], y[tr]).decision_function(Zte[:, top])
        a.append(roc_auc_score(y, oof))
    return float(np.mean(a)), float(np.std(a))


base = cv(Xv)
S["model_vs_sex"] = dict(base=base)
if sexdf is not None:
    S["model_vs_sex"]["with_sex"] = cv(Xv, sexdf.inferred_male.values.astype(float)[:, None])
    i_m = sexdf.inferred_male.values
    for lb, sel in [("male_only", i_m), ("female_only", ~i_m)]:
        if sel.sum() >= 25 and 3 <= y[sel].sum() <= sel.sum() - 3:
            Xs, ys, bs = Xv[sel], y[sel], b[sel]
            oof = np.zeros(len(ys))
            for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEED).split(Xs, ys):
                Xtr, Xte = adjust_batch(Xs[tr], ys[tr], bs[tr], Xs[te], bs[te])
                mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
                Ztr, Zte = (Xtr - mu) / sd, (Xte - mu) / sd
                F_, _ = f_classif(Ztr, ys[tr]); top = np.argsort(np.nan_to_num(F_))[::-1][:200]
                oof[te] = make_model("enet", SEED).fit(Ztr[:, top], ys[tr]).decision_function(Zte[:, top])
            S["model_vs_sex"][lb] = dict(n=int(sel.sum()), responders=int(ys.sum()), auc=float(roc_auc_score(ys, oof)))
print("model vs sex:", S["model_vs_sex"])

# ------------------------------------------------------------------ (4) derivation vs validation comparability
meta = pd.read_csv(COH / "all_meta.csv")
mm = meta[(meta.cohort == "GSE328059")]
unl = mm[mm.timepoint == "unlabelled"]["sample"].astype(str)
unl = [c for c in unl if c in full.columns]
common_genes = X.columns
A = full.loc[common_genes, lab].T.values
B = full.loc[common_genes, unl].T.values
pooled_sd = np.sqrt((A.var(0, ddof=1) + B.var(0, ddof=1)) / 2) + 1e-9
shift = (B.mean(0) - A.mean(0)) / pooled_sd
Zall = np.vstack([A, B]); Zall = (Zall - Zall.mean(0)) / (Zall.std(0) + 1e-9)
p2 = PCA(2, random_state=SEED).fit_transform(Zall)
lab_mask = np.r_[np.ones(len(A), bool), np.zeros(len(B), bool)]
S["derivation_vs_validation"] = dict(
    n_labelled=int(len(A)), n_unlabelled=int(len(B)),
    pct_genes_shift_gt1sd=float((np.abs(shift) > 1).mean() * 100),
    pct_genes_shift_gt2sd=float((np.abs(shift) > 2).mean() * 100),
    median_abs_shift=float(np.median(np.abs(shift))),
    pc1_separation_p=float(stats.ttest_ind(p2[lab_mask, 0], p2[~lab_mask, 0], equal_var=False).pvalue),
    pc1_auc=float(max(roc_auc_score(lab_mask, p2[:, 0]), 1 - roc_auc_score(lab_mask, p2[:, 0]))))
print("derivation vs validation:", S["derivation_vs_validation"])

json.dump(S, open(RES / "diagnostics_summary.json", "w"), indent=1, default=str)

fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.3))
ax[0].hist(sig.I2, bins=25, color="#2e86c1", edgecolor="k", lw=.4)
ax[0].axvline(sig.I2.median(), color="k", ls="--", lw=1.2)
ax[0].set_xlabel("$I^2$ (%)"); ax[0].set_ylabel("Consensus HF genes"); ax[0].set_title("Between-cohort heterogeneity", fontsize=10)
i = np.arange(len(V)); w = 0.38
ax[1].bar(i - w / 2, V.batch_r2, w, color="#b03a2e", edgecolor="k", lw=.5, label="sequencing batch")
ax[1].bar(i + w / 2, V.outcome_r2, w, color="#1f4e9e", edgecolor="k", lw=.5, label="recovery outcome")
ax[1].set_xticks(i); ax[1].set_xticklabels([f"PC{k}" for k in V.PC]); ax[1].set_ylabel("Variance explained (%)")
ax[1].legend(fontsize=7); ax[1].set_title("Structure of the 84 pre-LVAD hearts", fontsize=10)
ax[2].scatter(p2[lab_mask, 0], p2[lab_mask, 1], s=8, color="#1f4e9e", label="labelled (84)")
ax[2].scatter(p2[~lab_mask, 0], p2[~lab_mask, 1], s=8, color="#e67e22", label="withheld labels (134)")
ax[2].set_xlabel("PC1"); ax[2].set_ylabel("PC2"); ax[2].legend(fontsize=7)
ax[2].set_title("Derivation vs validation hearts", fontsize=10)
for a, l in zip(ax, "abc"):
    a.text(-0.12, 1.05, f"({l})", transform=a.transAxes, fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG / "FigS21_diagnostics.png", dpi=600)
fig.savefig(FIG / "FigS21_diagnostics.tif", dpi=600, pil_kwargs={"compression": "tiff_lzw", "dpi": (600, 600)})
print("wrote", FIG / "FigS21_diagnostics.tif")
