# Ready-to-paste text for the AHJ manuscript (stage 11)

## Methods — new paragraph

Cell-type composition and pathway-level modelling. Bulk pre-LVAD profiles were deconvolved against
a reference built from the GSE226314 single-nucleus pseudobulk (13 cell types, cell types with at
least 30 nuclei per sample). The 100 most specific genes per cell type (1,297 unique markers) formed
the signature matrix, and fractions were estimated by non-negative least squares and normalised to
sum to one. Fractions were compared between responders and non-responders, between sequencing
batches and between the consensus subtypes, used as the sole predictors of recovery, and removed
from the expression matrix inside each training fold (ordinary least squares on the fractions,
fitted in training data only) before repeating the main model. In parallel, every heart was scored
by single-sample GSEA on MSigDB Hallmark and TRRUST transcription-factor target sets (329 sets
passing size filters), and the same leakage-controlled pipeline, held-out-batch tests and frozen
external test were run on these rank-based set scores.

## Results — new paragraph (after the subtypes/reversal section)

Cellular composition explains part of the recovery signal. Deconvolution estimated a mean
cardiomyocyte fraction of 0.61, with fibroblasts 0.12 and smooth muscle 0.12. No cell type differed
between responders and non-responders after correction (smallest FDR 0.26; pericytes nominally
higher in responders, 0.034 vs 0.019, p = 0.02), but composition differed strongly by sequencing
batch (smooth muscle p < 0.001, pericytes p < 0.001, epicardium p < 0.001) and between the consensus
subtypes, which separated on fibroblast (0.14 vs 0.05) and cardiomyocyte (0.52 vs 0.64) content
(both FDR < 0.001) — confirming that the subtypes track regional sampling composition rather than a
recovery-related program. Cell fractions alone predicted recovery with a cross-validated AUC of 0.61
± 0.04, and removing the composition component inside each fold lowered the transcriptome model from
0.72 ± 0.05 to 0.68 ± 0.05, indicating that part, but not all, of the internal signal is cellular
composition (Supplementary Figure S17).

Pathway-level models are batch-robust but still do not replicate. Scoring the same hearts on 329
Hallmark and transcription-factor target sets gave a cross-validated AUC of 0.73 ± 0.05, essentially
unchanged without batch adjustment (0.73), and held-out-batch AUCs of 0.84 and 0.84 — substantially
higher than the gene-level model on the same held-out batches. Eighteen sets reached FDR < 0.05
(versus one gene at gene level), all lower in responders, led by NFE2L2 and ATF3 target sets, the p53
pathway and the UV-response-up set. Applied once to the independent cohort, however, the pathway-level
model again failed to replicate (AUC 0.54, 95% CI 0.25–0.82, exact p = 0.78; Supplementary Figure S18).
Robust, low-dimensional features therefore remove the batch sensitivity of the internal signal without
rescuing external replication, which points to cohort differences or limited power rather than feature
fragility as the reason replication fails.

## Discussion — sentences to add

Two further analyses sharpen the interpretation. First, composition is a real contributor: cell
fractions alone discriminate responders at an AUC of 0.61 and removing them costs the full model
about 0.04 AUC, so a transcriptomic classifier built on apical core tissue partly measures what was
sampled rather than how the myocardium will behave. Second, moving from genes to pathway and
transcription-factor target activities makes the internal signal markedly more portable across
sequencing batches (held-out-batch AUC 0.84 versus 0.67–0.75 at gene level) yet leaves external
replication unchanged, which argues that the barrier is the cohort and outcome definition rather
than the representation of the data.

## Supplementary Figure S17 — legend

Supplementary Figure S17. Cell-type composition of the pre-LVAD hearts. (a) Estimated fractions of
the eight most abundant cell types in responders and non-responders (boxes show median and
interquartile range, whiskers 1.5 x IQR, outliers omitted). (b) Composition of the two consensus
subtypes. (c) Cross-validated AUC for recovery prediction from cell fractions alone, from the
6,000-gene transcriptome, and from the transcriptome after in-fold removal of the composition
component (mean and SD over ten repeats). Rendered at 600 dpi.

## Supplementary Figure S18 — legend

Supplementary Figure S18. Pathway- and transcription-factor-level recovery modelling. (a) The ten
Hallmark or TRRUST sets most strongly associated with recovery; blue, higher in responders; red,
lower. (b) AUC of the set-level model under cross-validation, on each held-out sequencing batch, and
in the frozen external test. (c) Distribution of out-of-fold scores by outcome. Rendered at 600 dpi.

## Supplementary tables

- Supplementary Table S6. Estimated cell-type fractions for all 218 pre-LVAD hearts
  (`deconvolution_fractions.csv`) and their association with recovery, batch and subtype
  (`deconvolution_associations.csv`).
- Supplementary Table S7. Single-sample GSEA scores and univariate recovery associations for the 329
  Hallmark and TRRUST sets (`pathway_scores_GSE328059.csv`, `pathway_univariate.csv`).

## Numbers to cite

| quantity | value |
|---|---|
| mean cardiomyocyte / fibroblast / SMC fraction | 0.61 / 0.12 / 0.12 |
| best cell type for recovery (uncorrected) | pericytes, p = 0.020, FDR 0.26 |
| composition vs sequencing batch | SMC p < 0.0001, pericyte p = 0.0002, epicardium p = 0.0008 |
| composition vs consensus subtype | fibroblast 0.14 vs 0.05, cardiomyocyte 0.52 vs 0.64, both FDR < 0.001 |
| AUC from fractions alone | 0.612 ± 0.035 |
| AUC transcriptome / after composition removal | 0.717 ± 0.046 / 0.675 ± 0.051 |
| pathway-level CV AUC (adjusted / unadjusted) | 0.726 ± 0.048 / 0.733 |
| pathway-level held-out batch AUC | 0.844 (11943), 0.842 (10931) |
| sets at FDR < 0.05 | 18 of 329 |
| pathway-level external AUC | 0.544 (95% CI 0.25–0.82), exact p = 0.78 |

Note: AUC values here come from the public re-implementation, whose gene-level cross-validated AUC is
0.72 against 0.75 in the manuscript (see VALIDATION.md). Quote the pathway and composition numbers
relative to the 0.72 rebuild value, not to 0.75, when they are compared directly.
