#!/usr/bin/env python3
"""Stage 11a: cell-type deconvolution of the bulk pre-LVAD hearts, and whether the recovery
signal survives adjustment for cellular composition.

Reference: GSE226314 single-nucleus pseudobulk by cell type (13 types). Markers are the genes most
specific to each type; fractions are estimated by non-negative least squares and normalised to 1.

Tests: (i) fractions vs recovery, (ii) fractions vs sequencing batch, (iii) fractions vs the S1/S2
consensus subtypes, (iv) recovery prediction from fractions alone, (v) the usual 6000-gene model
after in-fold removal of the composition component.

usage: python deconvolution.py
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "05_recovery"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "02_meta_lasso_shap"))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import nnls
from sklearn.feature_selection import f_classif
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from pipeline import load_pre, adjust_batch, make_model
from common import COH, RES, SEED, bh

FIG = Path(__file__).resolve().parent.parent / "figures"
plt.rcParams.update({"font.family": "serif", "font.size": 10, "axes.linewidth": 1.1})
NMARK, NSEL, REP = 100, 200, 10

pb = pd.read_csv(COH / "GSE226314_pseudobulk_by_celltype.csv", index_col=0)
nuc = pd.read_csv(COH / "GSE226314_nuclei_per_sample_celltype.csv")
ok = {f"{r.Var1}__{r.Var2}" for r in nuc.itertuples() if r.Freq >= 30}
pb = pb[[c for c in pb.columns if c in ok]]
ct = pd.Series([c.split("__")[-1] for c in pb.columns], index=pb.columns)
lin = np.power(2, pb) - 1
ref = lin.T.groupby(ct.values).mean().T
ref = ref[ref.sum(1) > 0]
print("reference", ref.shape, "cell types:", list(ref.columns))

X, y, b, m = load_pre()
expr = pd.read_csv(COH / "GSE328059_expr.csv.gz", index_col=0)
meta = pd.read_csv(COH / "all_meta.csv")
mm = meta[(meta.cohort == "GSE328059") & meta.timepoint.isin(["pre", "unlabelled"])]
bulk = expr[[c for c in mm["sample"].astype(str) if c in expr.columns]]
genes = ref.index.intersection(bulk.index)
ref, bulk = ref.loc[genes], bulk.loc[genes]
bulk_lin = np.power(2, bulk) - 1

spec = ref.div(ref.sum(1) + 1e-9, axis=0)
expressed = ref.max(1) > 5
marks = sorted({g for c in ref.columns for g in spec.loc[expressed, c].sort_values(ascending=False).head(NMARK).index})
R = ref.loc[marks].values
R = R / (R.sum(0, keepdims=True) + 1e-9)
print("markers used:", len(marks))

F = []
for s in bulk.columns:
    v = bulk_lin[s].loc[marks].values
    w, _ = nnls(R, v)
    F.append(w / (w.sum() + 1e-12))
F = pd.DataFrame(F, index=bulk.columns, columns=ref.columns)
F.to_csv(RES / "deconvolution_fractions.csv")
print(F.mean().round(3).sort_values(ascending=False).to_string())

lab = m["sample"].astype(str).values
Fl = F.loc[lab]

rows = []
subf = RES / "subtypes.csv"
sub = pd.read_csv(subf) if subf.exists() else None
sc = None
if sub is not None:
    sub.columns = [c.lower() for c in sub.columns]
    sc = pd.Series(sub["subtype"].astype(str).values, index=sub["sample"].astype(str)).reindex(F.index).values
for c in F.columns:
    v = Fl[c].values
    t_r = stats.ttest_ind(v[y == 1], v[y == 0], equal_var=False)
    i1, i2 = b == "10931", b == "11943"
    p_b = stats.ttest_ind(v[i1], v[i2], equal_var=False).pvalue if i1.sum() > 2 and i2.sum() > 2 else np.nan
    row = dict(cell_type=c, mean_fraction=v.mean(), R_mean=v[y == 1].mean(), NR_mean=v[y == 0].mean(),
               p_recovery=t_r.pvalue, p_batch=p_b)
    if sc is not None:
        s1, s2 = sc == "S1", sc == "S2"
        if True:
            if s1.sum() > 2 and s2.sum() > 2:
                row["S1_mean"] = F[c].values[s1].mean(); row["S2_mean"] = F[c].values[s2].mean()
                row["p_subtype"] = stats.ttest_ind(F[c].values[s1], F[c].values[s2], equal_var=False).pvalue
    rows.append(row)
A = pd.DataFrame(rows)
A["fdr_recovery"] = bh(A.p_recovery.fillna(1).values)
if "p_subtype" in A: A["fdr_subtype"] = bh(A.p_subtype.fillna(1).values)
A.sort_values("p_recovery").to_csv(RES / "deconvolution_associations.csv", index=False)
print(A.sort_values("p_recovery").round(4).to_string(index=False))


def cv_simple(Z, y, repeats=REP):
    a = []
    for r in range(repeats):
        oof = np.zeros(len(y))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEED + r).split(Z, y):
            mu, sd = Z[tr].mean(0), Z[tr].std(0) + 1e-9
            mo = make_model("enet", SEED + r).fit((Z[tr] - mu) / sd, y[tr])
            oof[te] = mo.decision_function((Z[te] - mu) / sd)
        a.append(roc_auc_score(y, oof))
    return np.array(a)


auc_frac = cv_simple(Fl.values, y)
print(f"recovery from cell fractions alone: AUC {auc_frac.mean():.3f} +/- {auc_frac.std():.3f}")


def remove_comp(Xtr, Ctr, Xte, Cte):
    Dtr = np.column_stack([np.ones(len(Ctr)), Ctr[:, :-1]])
    Dte = np.column_stack([np.ones(len(Cte)), Cte[:, :-1]])
    beta, *_ = np.linalg.lstsq(Dtr, Xtr, rcond=None)
    return Xtr - Dtr[:, 1:] @ beta[1:], Xte - Dte[:, 1:] @ beta[1:]


def cv_full(Xv, y, b, C=None, repeats=REP):
    a = []
    for r in range(repeats):
        oof = np.zeros(len(y))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEED + r).split(Xv, y):
            Xtr, Xte = adjust_batch(Xv[tr], y[tr], b[tr], Xv[te], b[te])
            if C is not None:
                Xtr, Xte = remove_comp(Xtr, C[tr], Xte, C[te])
            mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
            Ztr, Zte = (Xtr - mu) / sd, (Xte - mu) / sd
            F_, _ = f_classif(Ztr, y[tr]); top = np.argsort(np.nan_to_num(F_))[::-1][:NSEL]
            mo = make_model("enet", SEED + r).fit(Ztr[:, top], y[tr])
            oof[te] = mo.decision_function(Zte[:, top])
        a.append(roc_auc_score(y, oof))
    return np.array(a)


Xv = X.values
auc_base = cv_full(Xv, y, b)
auc_adj = cv_full(Xv, y, b, Fl.values)
print(f"transcriptome AUC: batch-adjusted {auc_base.mean():.3f}, + composition-adjusted {auc_adj.mean():.3f}")

summary = dict(n_cell_types=int(F.shape[1]), n_markers=len(marks),
               auc_fractions_only=[round(float(auc_frac.mean()), 3), round(float(auc_frac.std()), 3)],
               auc_transcriptome=[round(float(auc_base.mean()), 3), round(float(auc_base.std()), 3)],
               auc_transcriptome_composition_adjusted=[round(float(auc_adj.mean()), 3), round(float(auc_adj.std()), 3)],
               top_recovery_cell_type=str(A.sort_values("p_recovery").iloc[0].cell_type),
               min_p_recovery=float(A.p_recovery.min()), min_fdr_recovery=float(A.fdr_recovery.min()),
               min_fdr_subtype=float(A.fdr_subtype.min()) if "fdr_subtype" in A else None)
json.dump(summary, open(RES / "deconvolution_summary.json", "w"), indent=1)
print(json.dumps(summary, indent=1))

order = F.mean().sort_values(ascending=False).index[:8]
fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.4))
for k, c in enumerate(order):
    ax[0].boxplot([Fl[c].values[y == 1]], positions=[k - 0.18], widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor="#1f4e9e", alpha=.6), medianprops=dict(color="k"), showfliers=False)
    ax[0].boxplot([Fl[c].values[y == 0]], positions=[k + 0.18], widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor="#b03a2e", alpha=.6), medianprops=dict(color="k"), showfliers=False)
ax[0].set_xticks(np.arange(len(order))); ax[0].set_xticklabels(order, rotation=40, ha="right", fontsize=7)
ax[0].set_ylabel("Estimated fraction"); ax[0].set_title("Composition by outcome", fontsize=10)
ax[0].plot([], [], color="#1f4e9e", lw=6, alpha=.6, label="responder"); ax[0].plot([], [], color="#b03a2e", lw=6, alpha=.6, label="non-responder")
ax[0].legend(fontsize=7)
if "S1_mean" in A:
    w = 0.38; i = np.arange(len(order)); Ai = A.set_index("cell_type")
    ax[1].bar(i - w / 2, Ai.loc[order, "S1_mean"], w, color="#2e86c1", edgecolor="k", lw=.5, label="S1")
    ax[1].bar(i + w / 2, Ai.loc[order, "S2_mean"], w, color="#e67e22", edgecolor="k", lw=.5, label="S2")
    ax[1].set_xticks(i); ax[1].set_xticklabels(order, rotation=40, ha="right", fontsize=7)
    ax[1].set_ylabel("Estimated fraction"); ax[1].set_title("Composition by consensus subtype", fontsize=10); ax[1].legend(fontsize=7)
ax[2].bar([0, 1, 2], [auc_frac.mean(), auc_base.mean(), auc_adj.mean()],
          yerr=[auc_frac.std(), auc_base.std(), auc_adj.std()], capsize=4,
          color=["#7d3c98", "#1f4e9e", "#1e8449"], edgecolor="k", lw=.6)
ax[2].axhline(0.5, color="k", ls="--", lw=1)
ax[2].set_xticks([0, 1, 2]); ax[2].set_xticklabels(["fractions\nonly", "transcriptome", "transcriptome\nminus composition"], fontsize=8)
ax[2].set_ylabel("Cross-validated AUC"); ax[2].set_ylim(0.4, 0.9); ax[2].set_title("Recovery prediction", fontsize=10)
for a, l in zip(ax, "abc"):
    a.text(-0.12, 1.05, f"({l})", transform=a.transAxes, fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG / "FigS17_deconvolution.png", dpi=600)
fig.savefig(FIG / "FigS17_deconvolution.tif", dpi=600, pil_kwargs={"compression": "tiff_lzw", "dpi": (600, 600)})
print("wrote", FIG / "FigS17_deconvolution.tif")
