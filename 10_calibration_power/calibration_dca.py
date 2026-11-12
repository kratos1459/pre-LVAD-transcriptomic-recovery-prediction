#!/usr/bin/env python3
"""Stage 10a: calibration and decision-curve analysis for the pre-LVAD recovery model.

Internal: 10 x 5-fold out-of-fold predicted probabilities (elastic net, in-fold batch adjustment).
External: frozen model fitted on all 84 labelled hearts, applied once to PRJNA1312133 (n = 19).
Outputs CSVs, a JSON summary and 600-dpi TIFF / PNG figures in results/ and figures/.

usage: python calibration_dca.py [../data/prjna]
"""
import sys, json, gzip, re
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "05_recovery"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "02_meta_lasso_shap"))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
from scipy.special import logit, expit
from sklearn.feature_selection import f_classif
from sklearn.metrics import roc_auc_score, brier_score_loss
from sklearn.model_selection import StratifiedKFold
from statsmodels.api import GLM, families
from pipeline import load_pre, adjust_batch, make_model
from common import RES, SEED

FIG = Path(__file__).resolve().parent.parent / "figures"
OUT = RES
plt.rcParams.update({"font.family": "serif", "font.size": 10, "axes.linewidth": 1.1})
REP, NSEL = 10, 200


def fold_proba(Xtr, ytr, btr, Xte, bte, seed=0, nsel=NSEL):
    """Identical to pipeline.run_fold(enet, 'inner') but returns probabilities."""
    Xtr, Xte = adjust_batch(Xtr, ytr, btr, Xte, bte)
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
    Ztr, Zte = (Xtr - mu) / sd, (Xte - mu) / sd
    F, _ = f_classif(Ztr, ytr)
    top = np.argsort(np.nan_to_num(F))[::-1][:nsel]
    m = make_model("enet", seed).fit(Ztr[:, top], ytr)
    return m.predict_proba(Zte[:, top])[:, 1]


def cal_metrics(y, p, label):
    p = np.clip(np.asarray(p, float), 1e-6, 1 - 1e-6); y = np.asarray(y, int); lp = logit(p)
    citl = GLM(y, np.ones((len(y), 1)), family=families.Binomial(), offset=lp).fit()
    slope = GLM(y, np.column_stack([np.ones(len(y)), lp]), family=families.Binomial()).fit()
    brier = brier_score_loss(y, p); prev = y.mean()
    bs_ref = brier_score_loss(y, np.repeat(prev, len(y)))
    nb = min(4, max(2, len(y) // 20))
    q = pd.qcut(pd.Series(p).rank(method="first"), nb, labels=False).values
    ece = float(np.sum([abs(y[q == k].mean() - p[q == k].mean()) * (q == k).mean() for k in range(nb)]))
    z = np.sum((y - p) * (1 - 2 * p)) / np.sqrt(np.sum((1 - 2 * p) ** 2 * p * (1 - p)))
    return dict(set=label, n=int(len(y)), events=int(y.sum()), auc=float(roc_auc_score(y, p)),
                brier=float(brier), brier_skill=float(1 - brier / bs_ref),
                citl=float(citl.params[0]), citl_lo=float(citl.conf_int()[0][0]), citl_hi=float(citl.conf_int()[0][1]),
                slope=float(slope.params[1]), slope_lo=float(slope.conf_int()[1][0]), slope_hi=float(slope.conf_int()[1][1]),
                ece=ece, spiegelhalter_z=float(z), spiegelhalter_p=float(2 * stats.norm.sf(abs(z))),
                observed=float(y.mean()), expected=float(p.mean()), e_o=float(p.mean() / y.mean()))


def dca(y, p, pts=np.arange(0.02, 0.76, 0.01)):
    n = len(y); rows = []
    for pt in pts:
        w = pt / (1 - pt); pos = p >= pt
        tp, fp = float(np.sum(y[pos] == 1)), float(np.sum(y[pos] == 0))
        nb_model = tp / n - fp / n * w
        nb_all = y.sum() / n - (n - y.sum()) / n * w
        rows.append((pt, nb_model, nb_all, 0.0, (nb_model - nb_all) / w * 1000))
    return pd.DataFrame(rows, columns=["threshold", "nb_model", "nb_treat_all", "nb_treat_none", "net_reduction_per_1000"])


def cal_plot(ax, y, p, title, nb=4):
    y = np.asarray(y, int); p = np.asarray(p, float)
    q = pd.qcut(pd.Series(p).rank(method="first"), nb, labels=False).values
    xs = np.array([p[q == k].mean() for k in range(nb)]); ys = np.array([y[q == k].mean() for k in range(nb)])
    ns = np.array([(q == k).sum() for k in range(nb)]); se = np.sqrt(ys * (1 - ys) / ns)
    hi = max(0.9, float(p.max()) * 1.1, float(ys.max()) * 1.1); lim = [0, hi]
    ax.plot(lim, lim, "k--", lw=1, label="ideal")
    ax.errorbar(xs, ys, yerr=1.96 * se, fmt="o", color="#1f4e9e", ms=5, capsize=3, label="grouped observed")
    lp = logit(np.clip(p, 1e-6, 1 - 1e-6)); gr = np.linspace(lp.min(), lp.max(), 100)
    fit = GLM(y, np.column_stack([np.ones(len(y)), lp]), family=families.Binomial()).fit()
    ax.plot(expit(gr), expit(fit.params[0] + fit.params[1] * gr), color="#b03a2e", lw=1.6, label="logistic recalibration")
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xlabel("Predicted probability of recovery")
    ax.set_ylabel("Observed proportion"); ax.set_title(title, fontsize=10); ax.legend(fontsize=7, loc="upper left")


X, y, b, m = load_pre(); Xv = X.values
P = np.zeros((REP, len(y)))
for r in range(REP):
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEED + r).split(Xv, y):
        P[r, te] = fold_proba(Xv[tr], y[tr], b[tr], Xv[te], b[te], SEED + r)
oof_p = P.mean(0)
pd.DataFrame({"sample": m["sample"].values, "y": y, "batch": b, "oof_prob": oof_p}).to_csv(OUT / "calibration_oof_prob.csv", index=False)
res = [cal_metrics(y, oof_p, "internal_oof_84")]
dca_int = dca(y, oof_p); dca_int.insert(0, "set", "internal_oof_84")

dca_ext = None; p_ext = None; yext = None
if len(sys.argv) > 1 and Path(sys.argv[1]).exists():
    Q = Path(sys.argv[1]) / "quant"
    hdr = {}
    for l in gzip.open(Q.parent / "ref" / "gencode.v46.transcripts.fa.gz", "rt"):
        if l[0] == ">":
            f = l[1:].split("|"); hdr[f[0]] = f[5]
    rows = {}
    for d in sorted(Q.iterdir()):
        if re.search(r"_HF_(R|NR)\d+$", d.name):
            t = pd.read_csv(d / "quant.sf", sep="\t", index_col=0); t["g"] = t.index.map(hdr)
            rows[d.name] = np.log2(t.groupby("g").TPM.sum() + 1)
    E = pd.DataFrame(rows)
    yext = np.array([int(re.search(r"_HF_(R|NR)", c).group(1) == "R") for c in E.columns])
    genes = [g for g in X.columns if g in E.index]
    Xd = X[genes].values; Xe = E.loc[genes].T.values
    Xe = (Xe - Xe.mean(0)) / (Xe.std(0) + 1e-9) * (Xd.std(0) + 1e-9) + Xd.mean(0)
    p_ext = fold_proba(Xd, y, b, Xe, np.array(["other"] * len(Xe)), SEED)
    pd.DataFrame({"sample": E.columns, "y": yext, "prob": p_ext}).to_csv(OUT / "calibration_external_prob.csv", index=False)
    res.append(cal_metrics(yext, p_ext, "external_prjna1312133_19"))
    dca_ext = dca(yext, p_ext); dca_ext.insert(0, "set", "external_prjna1312133_19")

cal = pd.DataFrame(res); cal.to_csv(OUT / "calibration_metrics.csv", index=False)
pd.concat([d for d in [dca_int, dca_ext] if d is not None]).to_csv(OUT / "decision_curve.csv", index=False)
json.dump(res, open(OUT / "calibration_summary.json", "w"), indent=1)
print(cal.round(3).to_string(index=False))

nset = 2 if p_ext is not None else 1
fig, axes = plt.subplots(nset, 2, figsize=(7.2, 3.3 * nset), squeeze=False)
cal_plot(axes[0][0], y, oof_p, "Internal (84 hearts, out-of-fold)")
ax = axes[0][1]; ax.plot(dca_int.threshold, dca_int.nb_model, color="#1f4e9e", lw=1.8, label="recovery model")
ax.plot(dca_int.threshold, dca_int.nb_treat_all, color="grey", lw=1.2, ls="--", label="treat all as responders")
ax.axhline(0, color="k", lw=1, label="treat none")
ax.set_ylim(min(-0.05, float(dca_int.nb_model.min())), max(0.1, float(dca_int.nb_model.max()) * 1.3))
ax.set_xlabel("Threshold probability"); ax.set_ylabel("Net benefit"); ax.set_title("Internal decision curve", fontsize=10); ax.legend(fontsize=7)
if p_ext is not None:
    cal_plot(axes[1][0], yext, p_ext, "External (PRJNA1312133, 19 hearts)", nb=3)
    ax = axes[1][1]; ax.plot(dca_ext.threshold, dca_ext.nb_model, color="#1f4e9e", lw=1.8, label="recovery model")
    ax.plot(dca_ext.threshold, dca_ext.nb_treat_all, color="grey", lw=1.2, ls="--", label="treat all as responders")
    ax.axhline(0, color="k", lw=1, label="treat none")
    ax.set_ylim(min(-0.05, float(dca_ext.nb_model.min())), max(0.1, float(dca_ext.nb_model.max()) * 1.3))
    ax.set_xlabel("Threshold probability"); ax.set_ylabel("Net benefit"); ax.set_title("External decision curve", fontsize=10); ax.legend(fontsize=7)
for a, l in zip(np.ravel(axes), "abcd"):
    a.text(-0.14, 1.04, f"({l})", transform=a.transAxes, fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG / "FigS16_calibration_dca.png", dpi=600)
fig.savefig(FIG / "FigS16_calibration_dca.tif", dpi=600, pil_kwargs={"compression": "tiff_lzw", "dpi": (600, 600)})
print("wrote", FIG / "FigS16_calibration_dca.tif")
