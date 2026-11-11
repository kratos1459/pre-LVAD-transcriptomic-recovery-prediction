#!/usr/bin/env bash
# Full pipeline. Stage 06 is heavy (hours, ~215 GB transfer); run separately.
set -e
(cd 01_data && bash download.sh && python3 prepare_cohorts.py)
(cd 02_meta_lasso_shap && python3 run_meta.py)
(cd 03_loco && python3 run_loco.py)
(cd 04_annotation && python3 run_annotation.py)
(cd 05_recovery && python3 run_recovery.py)
(cd 07_subtyping && python3 run_subtyping.py)
(cd 09_figures && python3 make_figures.py)
# Stage 06 (optional, heavy):
#   bash 06_transfer_external/quantify_prjna1312133.sh data/prjna 2 8
#   (cd 06_transfer_external && python3 frozen_external_test.py ../data/prjna)
#   Rscript 06_transfer_external/pseudobulk_gse226314.R GSE226314_global.rds results/cohorts
#   (cd 06_transfer_external && python3 transfer_gse226314.py)
