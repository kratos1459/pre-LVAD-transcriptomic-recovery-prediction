#!/usr/bin/env python3
"""Stage 10c: sealed (pre-specified) predictions for the 134 GSE328059 validation hearts whose
recovery labels are withheld.

The elastic-net recovery model is fitted once on all 84 labelled pre-LVAD hearts exactly as in
stage 05 (in-fold-style batch adjustment, F-test selection of 200 genes), then applied once to the
134 unlabelled hearts. Two score versions are frozen:
  prob_raw       - model applied directly (the validation hearts are all from sequencing batches
                   outside the derivation batches, so no batch term applies to them);
  prob_rescaled  - same protocol as the frozen external test of stage 06b: standardise within the
                   validation cohort, rescale to the derivation moments, then score. PRIMARY score.
Predictions are hashed and deposited so the model can be tested prospectively when labels are released.

usage: python sealed_predictions.py
"""
import sys, json, hashlib, datetime
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "05_recovery"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "02_meta_lasso_shap"))
import numpy as np, pandas as pd
from sklearn.feature_selection import f_classif
from sklearn.model_selection import StratifiedKFold
from pipeline import load_pre, adjust_batch, make_model
from common import COH, RES, SEED


def fold_proba(Xtr, ytr, btr, Xte, bte, seed=0, nsel=200):
    Xtr, Xte = adjust_batch(Xtr, ytr, btr, Xte, bte)
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
    Ztr, Zte = (Xtr - mu) / sd, (Xte - mu) / sd
    F, _ = f_classif(Ztr, ytr)
    top = np.argsort(np.nan_to_num(F))[::-1][:nsel]
    return make_model("enet", seed).fit(Ztr[:, top], ytr).predict_proba(Zte[:, top])[:, 1]


X, y, b, m = load_pre()
meta = pd.read_csv(COH / "all_meta.csv")
u = meta[(meta.cohort == "GSE328059") & (meta.timepoint == "unlabelled")].copy()
expr = pd.read_csv(COH / "GSE328059_expr.csv.gz", index_col=0)
cols = [c for c in u["sample"].astype(str) if c in expr.columns]
if len(cols) < len(u):
    print(f"WARNING: {len(u) - len(cols)} of {len(u)} unlabelled samples absent from the expression matrix")
u = u[u["sample"].astype(str).isin(cols)]
Xd = X.values
Xu = expr.loc[X.columns, u["sample"].astype(str)].T.values
bu = u["batch"].astype(float).astype(int).astype(str).values
bu = np.where(np.isin(bu, ["10931", "11943"]), bu, "other")
Xu_rs = (Xu - Xu.mean(0)) / (Xu.std(0) + 1e-9) * (Xd.std(0) + 1e-9) + Xd.mean(0)

p_raw = fold_proba(Xd, y, b, Xu, bu, SEED)
p_rs = fold_proba(Xd, y, b, Xu_rs, np.array(["other"] * len(Xu)), SEED)

# honest decision threshold: prevalence-matched quantile of the OUT-OF-FOLD probabilities
oof = np.zeros(len(y))
for tr, te in StratifiedKFold(5, shuffle=True, random_state=SEED).split(Xd, y):
    oof[te] = fold_proba(Xd[tr], y[tr], b[tr], Xd[te], b[te], SEED)
thr = float(np.quantile(oof, 1 - y.mean()))

out = pd.DataFrame({"sample": u["sample"].values, "gsm": u["gsm"].values, "batch": bu,
                    "prob_rescaled_PRIMARY": np.round(p_rs, 6), "prob_raw": np.round(p_raw, 6),
                    "class_at_oof_threshold": np.where(p_rs >= thr, "R", "NR"),
                    "class_top31pct_within_cohort": np.where(
                        pd.Series(p_rs).rank(pct=True) > 1 - y.mean(), "R", "NR")})
out = out.sort_values("prob_rescaled_PRIMARY", ascending=False)
f = RES / "sealed_predictions_GSE328059_validation.csv"
out.to_csv(f, index=False)
h = hashlib.sha256(open(f, "rb").read()).hexdigest()

plan = {
    "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
    "file": f.name, "sha256": h, "n_predicted": int(len(out)),
    "derivation_prevalence": round(float(y.mean()), 3),
    "decision_threshold_from_oof": round(thr, 6),
    "predicted_responders_at_oof_threshold": int((out.class_at_oof_threshold == "R").sum()),
    "mean_predicted_probability_rescaled": round(float(p_rs.mean()), 3),
    "mean_predicted_probability_raw": round(float(p_raw.mean()), 3),
    "spearman_raw_vs_rescaled": round(float(pd.Series(p_raw).corr(pd.Series(p_rs), method="spearman")), 3),
    "model": "elastic-net logistic regression (l1_ratio 0.5, C 0.1), 6000 most variable genes, "
             "batch adjustment and F-test selection of 200 genes fitted on the 84 labelled hearts only",
    "distribution_shift_note": "All validation hearts come from sequencing runs outside the two derivation "
        "batches, so no fitted batch term applies to them and the raw probabilities are shifted upward "
        "relative to the derivation prevalence. The c-statistic depends only on ranks and is unaffected; "
        "calibration metrics are expected to be poor and are reported as such.",
    "prespecified_primary_analysis": "When the labels of the GSE328059 validation cohort are released, "
        "compute the c-statistic of prob_rescaled_PRIMARY against the released labels, with a "
        "5000-bootstrap 95% confidence interval and a two-sided permutation p-value. No refitting, "
        "rescaling or threshold change is permitted.",
    "prespecified_secondary_analyses": [
        "The same for prob_raw.",
        "Calibration-in-the-large, calibration slope, Brier score and Brier skill score.",
        "Decision-curve net benefit across threshold probabilities 0.02-0.75.",
        "Sensitivity: c-statistic within each sequencing batch separately."],
    "success_criterion": "A 95% confidence interval for the c-statistic excluding 0.5.",
    "power_note": "At the derivation prevalence (about 42 responders among 134), this cohort has about "
                  "97% power to detect a true c-statistic of 0.70 and about 80% power at 0.65.",
}
json.dump(plan, open(RES / "sealed_predictions_plan.json", "w"), indent=1)
print(json.dumps(plan, indent=1))
print("\nTop 3 / bottom 3:\n", pd.concat([out.head(3), out.tail(3)]).to_string(index=False))
