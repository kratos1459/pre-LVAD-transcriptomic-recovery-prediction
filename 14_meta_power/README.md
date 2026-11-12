# Stage 14 — joint analysis of both labelled cohorts, learning curve, seed stability

`joint_meta_learning.py`: DerSimonian-Laird pooling of per-gene recovery effects across GSE328059
and PRJNA1312133; learning curve by subsampling with an asymptotic fit; the pipeline repeated under
20 master seeds; directional concordance in the external cohort of the sets significant internally.
Outputs: `joint_recovery_meta.csv`, `learning_curve.csv`, `pathway_external_concordance.csv`,
`joint_meta_learning_summary.json`, `FigS20_joint_meta_learning.{png,tif}`.
