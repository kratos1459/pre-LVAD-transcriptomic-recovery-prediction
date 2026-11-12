# Stage 15 — robustness diagnostics

`signature_and_batch_diagnostics.py`: heterogeneity (I-squared, tau-squared), leave-one-cohort-out
retention and microarray/RNA-seq concordance of the consensus HF signature; variance of each
principal component explained by batch and by outcome; sex inferred from X/Y transcripts and its
relation to recovery; comparability of the 84 labelled and 134 withheld-label hearts. Outputs:
`pc_variance_batch_outcome.csv`, `inferred_sex.csv`, `diagnostics_summary.json`,
`FigS21_diagnostics.{png,tif}`.
