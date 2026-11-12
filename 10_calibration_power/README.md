# Stage 10 — calibration, decision curves, power/precision, sealed predictions

Additions to the leakage-controlled benchmark that close the "calibration and decision-curve analysis
were not performed" limitation and quantify what the external test could and could not have detected.

| script | what it does | outputs (in `../results`, figures in `../figures`) |
|---|---|---|
| `calibration_dca.py` | 10 x 5-fold out-of-fold predicted probabilities; calibration-in-the-large, calibration slope, Brier, Brier skill, ECE, Spiegelhalter z; decision curves. Same for the frozen external test (PRJNA1312133). | `calibration_metrics.csv`, `calibration_oof_prob.csv`, `calibration_external_prob.csv`, `decision_curve.csv`, `calibration_summary.json`, `FigS16_calibration_dca.{png,tif}` |
| `power_precision.py` | Hanley–McNeil analytic power, exact-permutation power by simulation, minimum detectable c-statistic at 80% power, sample size for a target confidence-interval width. | `power_by_cohort.csv`, `power_external_exact.csv`, `power_precision_summary.json` |
| `sealed_predictions.py` | Freezes predictions for the 134 withheld-label GSE328059 validation hearts with a SHA-256 hash and a pre-specified analysis plan. | `sealed_predictions_GSE328059_validation.csv`, `sealed_predictions_plan.json` |

Run after stage 05 (and stage 06 quantification for the external part):

```bash
cd 10_calibration_power
python calibration_dca.py ../data/prjna
python power_precision.py
python sealed_predictions.py
```

Environment: the `lvad` conda environment (scikit-learn, statsmodels, xgboost, matplotlib).
`SUPPLEMENT_TEXT.md` holds the ready-to-paste Methods, Results, table and figure legend.
