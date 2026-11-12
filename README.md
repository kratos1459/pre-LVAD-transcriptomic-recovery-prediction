# Cross-cohort HF signature and pre-LVAD recovery benchmark: analysis code

Code for: *Cross-Cohort Heart-Failure Signature Generalisation and the Limits of Pre-LVAD Transcriptomic Recovery Prediction: A Leakage-Controlled Multi-Cohort Benchmark* (Briefings in Bioinformatics submission).
**Code only.** All inputs are public (NCBI GEO / ENA) and are fetched by the scripts. Licence: MIT.

## Pipeline

| Stage | Folder | Entry point | Needs | Status |
|---|---|---|---|---|
| 01 | `01_data` | `download.sh`, `prepare_cohorts.py` | internet | run |
| 02 | `02_meta_lasso_shap` | `run_meta.py` | 01 | run |
| 03 | `03_loco` | `run_loco.py` | 01, 02 | run |
| 04 | `04_annotation` | `run_annotation.py` (PubMed, Enrichr, DGIdb) | 02, internet | run |
| 05 | `05_recovery` | `run_recovery.py` | 01 | run |
| 06a | `06_transfer_external` | `quantify_prjna1312133.sh` (Salmon, GENCODE v46) | ~215 GB transfer, salmon | see log |
| 06b | `06_transfer_external` | `frozen_external_test.py` | 06a | see log |
| 06c | `06_transfer_external` | `pseudobulk_gse226314.R` | R/Seurat, ~64 GB RAM, 17 GB GEO object | see log |
| 06d | `06_transfer_external` | `transfer_gse226314.py` | 06c | see log |
| 07 | `07_subtyping` | `run_subtyping.py` | 05 | run |
| 08 | `08_secondary` | `secondary_cohorts.py` | GEO processed files | template, untested |
| 09 | `09_figures` | `make_figures.py` | 02-05 results | run |
| 10 | `10_calibration_power` | `calibration_dca.py`, `power_precision.py`, `sealed_predictions.py` | 05 (06a for the external part) | run |
| 11 | `11_deconvolution_pathways` | `deconvolution.py`, `pathway_activity.py` | 01, 06c, 07 (06a for the external part) | run |

Quick start: `conda env create -f environment.yml && conda activate bib-hf && bash run_all.sh`. Seeds are fixed in `02_meta_lasso_shap/common.py`. Tested on Linux, Python 3.12, 16 cores, 64 GB RAM; stages 01-05, 07 and 09 take under an hour in total.

## Important notes
- **Re-implementation.** The original analysis scripts were not recovered; this code was re-implemented from the manuscript Methods (with AI assistance, see the manuscript disclosure). Results are close to, but not identical with, the manuscript. See `VALIDATION.md`.
- **Unspecified choices.** Classifier regularisation, XGBoost settings, the 37-gene immediate-early list, inclusion of HCM/PPCM samples in GSE141910 and the 21-gene erythrocyte list are reconstructions and are marked in the code.
- **GSE226314** (8th discovery cohort) enters the meta-analysis only after `pseudobulk_gse226314.R` has produced its pseudobulk; until then stages 02-03 use 7 cohorts.
- **Data not included.** Downloaded and derived data (`data/`, `results/`) are git-ignored. The one exception is `10_calibration_power/sealed/`, which holds the small summary tables and the frozen predictions for the 134 withheld-label GSE328059 validation hearts (with their SHA-256 hash and pre-specified analysis plan), deposited so that the model can be tested prospectively when those labels are released.

## Citation
See `CITATION.cff`. Archive a GitHub release on Zenodo and cite the resulting DOI.
