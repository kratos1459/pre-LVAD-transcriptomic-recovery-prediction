# Stage 13 — head-to-head benchmark of recovery predictors

`benchmark_signatures.py` evaluates 22 fixed published and candidate scores (single genes, the
consensus HF signature, Hallmark and TRRUST set scores) in the 84 labelled hearts and once in
PRJNA1312133, and fits elastic-net models restricted to prior gene lists, compared with size-matched
random gene sets. Outputs: `benchmark_fixed_scores.csv`, `benchmark_fitted_models.csv`,
`benchmark_summary.json`, `FigS19_signature_benchmark.{png,tif}`.

Headline: no predictor is significant in both cohorts, and no prior gene list beats random sets of
the same size. `SUPPLEMENT_TEXT_stages_13_15.md` holds the manuscript text for stages 13-15.
