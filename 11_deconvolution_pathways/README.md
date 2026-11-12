# Stage 11 — cell-type deconvolution and pathway/TF-level modelling

| script | what it does | outputs |
|---|---|---|
| `deconvolution.py` | NNLS deconvolution of the 218 GSE328059 pre-LVAD hearts against the GSE226314 single-nucleus pseudobulk reference (13 cell types, 1,297 marker genes). Tests composition against recovery, sequencing batch and the consensus subtypes; predicts recovery from fractions alone; repeats the 6,000-gene model after in-fold removal of the composition component. | `deconvolution_fractions.csv`, `deconvolution_associations.csv`, `deconvolution_summary.json`, `FigS17_deconvolution.{png,tif}` |
| `pathway_activity.py` | Single-sample GSEA (gseapy) on MSigDB Hallmark and TRRUST TF-target sets (329 sets scored), then the stage-05 leakage-controlled pipeline on those scores, held-out-batch tests, univariate tests and a frozen external test on PRJNA1312133. | `pathway_scores_GSE328059.csv`, `pathway_univariate.csv`, `pathway_external_scores.csv`, `pathway_summary.json`, `FigS18_pathway_activity.{png,tif}`, `genesets.json` (cached libraries) |

```bash
cd 11_deconvolution_pathways
python deconvolution.py
python pathway_activity.py ../data/prjna
```

Needs stage 01 (cohorts), stage 06c (GSE226314 pseudobulk by cell type), stage 07 (subtypes) and,
for the external part, stage 06a. Internet is needed once to fetch the Enrichr libraries; they are
then cached in `genesets.json`.

Headline numbers are in `SUPPLEMENT_TEXT.md` together with ready-to-paste manuscript text.
