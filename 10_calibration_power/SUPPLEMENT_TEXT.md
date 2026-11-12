# Ready-to-paste text for the AHJ manuscript (stage 10)

## Methods — add to "Recovery modeling and batch audit"

Model performance was further assessed by calibration and clinical utility. Repeat-averaged
out-of-fold predicted probabilities were used to compute calibration-in-the-large (intercept of a
logistic model with the predicted logit as offset), the calibration slope, the Brier score, the
Brier skill score relative to the observed prevalence, the expected calibration error over quartiles
of predicted risk, and the Spiegelhalter z-test. Clinical utility was summarised by decision-curve
analysis across threshold probabilities of 0.02 to 0.75, comparing the model with treating all
patients as predicted responders and treating none. The same metrics were computed for the frozen
external test in PRJNA1312133.

The precision of the external test was quantified with the Hanley–McNeil variance of the
c-statistic. Because the permutation test applied to the c-statistic is equivalent to the exact
Wilcoxon rank-sum test, power was also obtained by simulation from a binormal model with that exact
test (20,000 replicates per scenario). Sample sizes required to estimate a c-statistic to a
pre-specified confidence-interval width were derived from the same variance.

Predictions for the 134 GSE328059 validation hearts with withheld labels were generated once with
the frozen model, hashed (SHA-256) and deposited with a pre-specified analysis plan, so that the
model can be tested prospectively when the labels are released.

## Results — add after "External replication fails..."

Calibration and clinical utility. Internally, out-of-fold predictions were well calibrated in the
large (calibration-in-the-large 0.11, 95% CI −0.38 to 0.59) with a calibration slope of 1.51 (95% CI
0.63 to 2.39), a Brier score of 0.178 (Brier skill score 0.17 relative to the observed prevalence),
an expected calibration error of 0.06 and a non-significant Spiegelhalter test (p = 0.64). Decision
curves showed net benefit above both default strategies for threshold probabilities above about 0.12,
with a net benefit of 0.122 at the prevalence threshold of 0.31, equivalent to 273 fewer patients per
1,000 incorrectly classified as responders than treating all as responders (Supplementary Figure S16a,b).
In the external cohort calibration was lost: the calibration slope fell to 0.07 (95% CI −1.45 to 1.58),
the Brier skill score was negative (−0.09, that is, worse than predicting the observed prevalence for
every patient), and the model gave no net benefit over treating all patients as responders below a
threshold probability of 0.45 (Supplementary Figure S16c,d).

What the external test could detect. With 10 recovered and 9 non-recovered hearts, the external cohort
had 44% power to detect a true c-statistic of 0.75 and 29% power at 0.70 (exact-permutation power 46%
and 30%, respectively); the minimum c-statistic detectable with 80% power was 0.85. The external result
is therefore compatible both with no discrimination and with moderate discrimination, and should be
read as inconclusive rather than as evidence of absence. To estimate a c-statistic of 0.75 with a 95%
confidence interval no wider than 0.10 would require approximately 484 patients (about 150 responders)
at the derivation prevalence, consistent with the rule of thumb of at least 100 events and 100
non-events for external validation of a binary prognostic model.

Sealed predictions. To allow prospective testing, probabilities for the 134 withheld-label GSE328059
validation hearts were frozen (SHA-256 b139bb3f…, Supplementary Table S5) with a pre-specified primary
analysis (c-statistic with a 5,000-bootstrap confidence interval and permutation p-value, no refitting).
That cohort would have approximately 97% power at a true c-statistic of 0.70 and 80% at 0.65.

## Discussion — optional sentence for Limitations

The external cohort was underpowered by design: it could not have detected anything below a
c-statistic of 0.85 with conventional power, so the failure to replicate constrains but does not
exclude a moderate effect.

## Supplementary Table S5 (new)

Frozen predictions for the 134 GSE328059 validation hearts (sample, GSM, sequencing batch, primary and
raw predicted probability, predicted class) with the pre-specified analysis plan and SHA-256 hash.
Source files: `sealed_predictions_GSE328059_validation.csv`, `sealed_predictions_plan.json`.

## Supplementary Figure S16 (new) — legend

Supplementary Figure S16. Calibration and clinical utility of the pre-LVAD recovery model.
(a) Calibration of repeat-averaged out-of-fold predictions in the 84 labelled hearts; points are
quartiles of predicted risk with 95% confidence intervals, the red line is the logistic recalibration
curve and the dashed line is perfect calibration. (b) Decision curve for the same predictions,
compared with treating all patients as predicted responders (grey dashed) and treating none (black).
(c, d) The same analyses for the frozen external test in PRJNA1312133 (19 hearts). Rendered at 600 dpi.

## Numbers to cite

| quantity | internal (n = 84) | external (n = 19) |
|---|---|---|
| c-statistic | 0.757 | 0.578 |
| Brier score | 0.178 | 0.273 |
| Brier skill score | 0.165 | −0.094 |
| calibration-in-the-large | 0.106 (−0.378 to 0.591) | 0.501 (−0.438 to 1.440) |
| calibration slope | 1.508 (0.631 to 2.385) | 0.067 (−1.447 to 1.581) |
| expected calibration error | 0.060 | 0.115 |
| Spiegelhalter p | 0.64 | 0.12 |
| net benefit exceeds treat-all from | threshold 0.12 | threshold 0.45 |

Power of the external cohort (10 responders, 9 non-responders): 10% at a true c-statistic of 0.60,
18% at 0.65, 29% at 0.70, 44% at 0.75, 63% at 0.80, 82% at 0.85. Minimum detectable c-statistic at
80% power: 0.845 (external), 0.689 (internal 84), 0.650 (withheld 134).
