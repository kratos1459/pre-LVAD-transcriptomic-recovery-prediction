#!/usr/bin/env python3
"""Stage 10b: how large an effect the external test could have detected, and how large a
validation cohort would need to be.

(1) Hanley-McNeil analytic power for a c-statistic vs the null c = 0.5.
(2) Exact-permutation power by simulation (the permutation test on AUC with n1/n0 samples is
    equivalent to the exact Wilcoxon rank-sum test, so power is simulated with that test).
(3) Minimum detectable c-statistic at 80% power.
(4) Sample size needed to estimate a c-statistic to a target confidence-interval width
    (precision criterion of Riley et al.), and the >=100 events / >=100 non-events rule of thumb.

usage: python power_precision.py
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from scipy.optimize import brentq

RES = Path(__file__).resolve().parent.parent / "results"
RNG = np.random.default_rng(20261024)
ALPHA = 0.05
NSIM = 20000


def se_auc(c, n1, n0):
    """Hanley & McNeil (1982) standard error of the c-statistic."""
    q1 = c / (2 - c); q2 = 2 * c ** 2 / (1 + c)
    return np.sqrt((c * (1 - c) + (n1 - 1) * (q1 - c ** 2) + (n0 - 1) * (q2 - c ** 2)) / (n1 * n0))


def power_analytic(c, n1, n0, alpha=ALPHA):
    z = stats.norm.ppf(1 - alpha / 2)
    return float(stats.norm.cdf((abs(c - 0.5) - z * se_auc(0.5, n1, n0)) / se_auc(c, n1, n0)))


def power_exact_perm(c, n1, n0, nsim=NSIM, alpha=ALPHA):
    """Binormal simulation; two-sided exact rank-sum test (= exact permutation test on AUC)."""
    d = np.sqrt(2) * stats.norm.ppf(c)  # binormal separation giving this AUC
    hits = 0
    for _ in range(nsim):
        a = RNG.normal(d, 1, n1); b = RNG.normal(0, 1, n0)
        if stats.mannwhitneyu(a, b, alternative="two-sided", method="exact").pvalue < alpha:
            hits += 1
    return hits / nsim


def min_detectable(n1, n0, target=0.80):
    return float(brentq(lambda c: power_analytic(c, n1, n0) - target, 0.501, 0.999))


def n_for_width(c, prev, width=0.10):
    """Total n so that the 95% CI of the c-statistic has the target width."""
    def f(n):
        n1 = max(n * prev, 2); n0 = max(n - n1, 2)
        return 2 * 1.96 * se_auc(c, n1, n0) - width
    return float(brentq(f, 20, 20000))


rows = []
# the two cohorts actually used
for lbl, n1, n0 in [("external PRJNA1312133", 10, 9), ("internal GSE328059 (84)", 26, 58),
                    ("withheld GSE328059 validation (134, assumed 31% responders)", 42, 92)]:
    for c in [0.60, 0.65, 0.70, 0.75, 0.80, 0.85]:
        rows.append(dict(cohort=lbl, n1=n1, n0=n0, true_c=c,
                         power_analytic=round(power_analytic(c, n1, n0), 3),
                         se_at_c=round(se_auc(c, n1, n0), 3)))
pw = pd.DataFrame(rows)
pw.to_csv(RES / "power_by_cohort.csv", index=False)

# exact-permutation power for the external cohort only (what the paper actually used)
ex = []
for c in [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90]:
    ex.append(dict(true_c=c, power_exact_perm=round(power_exact_perm(c, 10, 9), 3),
                   power_analytic=round(power_analytic(c, 10, 9), 3)))
exd = pd.DataFrame(ex); exd.to_csv(RES / "power_external_exact.csv", index=False)

# minimum detectable effect and required sample sizes
mde = {"external_19_mde_80pct": round(min_detectable(10, 9), 3),
       "internal_84_mde_80pct": round(min_detectable(26, 58), 3),
       "withheld_134_mde_80pct": round(min_detectable(42, 92), 3)}
need = {}
for c in [0.70, 0.75, 0.80]:
    for w in [0.10, 0.15]:
        n = n_for_width(c, 26 / 84, w)
        need[f"n_total_for_c{c}_CIwidth{w}"] = int(np.ceil(n))
        need[f"events_for_c{c}_CIwidth{w}"] = int(np.ceil(n * 26 / 84))
need["riley_rule_of_thumb"] = ">=100 events and >=100 non-events for external validation"
need["external_observed_c"] = 0.578
need["external_observed_ci"] = [0.27, 0.85]

summary = dict(alpha=ALPHA, nsim=NSIM, minimum_detectable_c_at_80pct_power=mde, precision_requirements=need)
json.dump(summary, open(RES / "power_precision_summary.json", "w"), indent=1)

print("Analytic power by cohort\n", pw.pivot(index="true_c", columns="cohort", values="power_analytic").to_string())
print("\nExternal cohort (10 R / 9 NR), exact permutation power\n", exd.to_string(index=False))
print("\nMinimum detectable c at 80% power:", mde)
print("\nSample size for a precise estimate:", {k: v for k, v in need.items() if k.startswith("n_total") or k.startswith("events")})
