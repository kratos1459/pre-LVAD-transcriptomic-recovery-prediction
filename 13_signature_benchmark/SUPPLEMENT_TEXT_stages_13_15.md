# Ready-to-paste text for stages 13-15

## Methods — new paragraph

Benchmark of published predictors, joint analysis and robustness. Published and candidate recovery
predictors were put through the same leakage-controlled frame: single genes (LRRN4CL, NPPB, ATF3,
CAMK2D, PCDH9, MYL7, NPPA, MYH7), signed scores for the consensus HF signature and for Hallmark and
TRRUST sets were evaluated as fixed scores in the 84 labelled hearts, with and without batch
centring, and once in PRJNA1312133; the same prior gene lists were also used to restrict feature
selection so that each list produced its own elastic-net model, which was compared with
size-matched random gene sets. Per-gene recovery effect sizes (Hedges g) were pooled across the two
outcome-labelled cohorts by DerSimonian-Laird random effects. A learning curve was obtained by
subsampling the 84 hearts (10 draws per size) and fitting AUC = a - c/sqrt(n); the whole pipeline was
also repeated under 20 master seeds. Heterogeneity of the consensus signature was summarised by
I-squared and tau-squared, by leave-one-cohort-out retention, and by microarray versus RNA-seq
concordance; Egger tests were not performed because at most eight cohorts contribute to any gene.
The variance of each principal component attributable to sequencing batch and to outcome was
computed by one-way analysis of variance. Sex was inferred from XIST against the mean of RPS4Y1,
DDX3Y, UTY, KDM5D and EIF1AY. Comparability of the 84 labelled and 134 withheld-label hearts was
assessed by standardised per-gene mean shifts and by principal components.

## Results — new paragraphs

No published or candidate predictor replicates. Of 22 fixed predictors, six were associated with
recovery in the 84 labelled hearts (ATF3 target set AUC 0.73, ATF3 itself 0.71, NFE2L2 targets 0.66,
LRRN4CL 0.64; five at FDR < 0.05) and two reached nominal significance in the independent cohort
(the consensus HF signature restricted to its top 150 genes, AUC 0.82, p = 0.017; Hallmark oxidative
phosphorylation, AUC 0.79, p = 0.035), but none survived correction externally and **no predictor was
significant in both cohorts**. When prior gene lists were instead used to restrict feature selection,
none outperformed size-matched random gene sets: the Hallmark union reached a cross-validated AUC of
0.72 against a random-set 95th percentile of 0.76, and the consensus HF signature 0.70 against 0.78
(Supplementary Figure S19, Supplementary Table S8).

The two labelled cohorts share no common gene-level signal. Pooling per-gene recovery effects across
GSE328059 and PRJNA1312133 by random effects left no gene at FDR < 0.05 (smallest FDR 0.056, CD1E).
Effect sizes were uncorrelated between the cohorts (Spearman r = 0.01) and agreed in direction for
50.3% of genes, indistinguishable from chance (binomial p = 0.43). The 18 gene sets that reached
FDR < 0.05 internally were concordant in direction in the external cohort for 7 of 18 (39%, p = 0.48),
and none reached p < 0.05 there.

How much data the problem needs. Cross-validated AUC rose from 0.60 at n = 30 to 0.72 at n = 84;
fitting AUC = a - c/sqrt(n) gave an asymptote of 0.84 and implied that about 190 labelled hearts
would be needed for a cross-validated AUC of 0.75 and about 1,000 for 0.80. Repeating the entire
pipeline under 20 master seeds gave 0.70 +/- 0.03 (range 0.65-0.76), so single-run estimates carry a
seed-level uncertainty of about +/- 0.05 (Supplementary Figure S20).

Robustness of the HF signature and structure of the cohort. Among the consensus HF genes the median
I-squared was 66%, with 78% above 50%, yet direction was unanimous between platforms (sign agreement
100%; Spearman r = 0.73 between microarray and RNA-seq pooled effects) and leave-one-cohort-out
removal retained at least 81.5% of the signature whichever cohort was dropped. In the 84 pre-LVAD
hearts, principal component 2 was 62.6% attributable to sequencing batch and 7.2% to outcome.
Inferred sex gave 69 male and 15 female hearts, with recovery in 26% of male and 53% of female
hearts (Fisher p = 0.06); adding sex to the model did not change performance, and the model was
weaker in men alone (AUC 0.63). The 134 withheld-label validation hearts were not exchangeable with
the 84 labelled hearts: 41% of genes shifted by more than one pooled standard deviation and the
first principal component separated the two groups completely (AUC 1.00), so the deposited sealed
predictions should be read with that distribution shift in mind (Supplementary Figure S21).

## Discussion — sentences to add

The benchmark is the strongest negative evidence in this study: across 22 published and candidate
predictors, prior gene lists that performed no better than random sets of the same size, and a joint
analysis of both labelled cohorts in which effect sizes were uncorrelated, nothing about the
pre-LVAD transcriptome transferred between cohorts. The learning curve suggests why this may be
partly a sample-size problem, with roughly 190 labelled hearts needed for a cross-validated AUC of
0.75, but the complete absence of cross-cohort agreement at gene and set level, and the fact that
the two labelled cohorts differ in platform and outcome ascertainment, point to cohort
heterogeneity as the more likely barrier. By contrast, the disease signature itself was unanimous in
direction across platforms, confirming that these data support robust descriptive biology and not,
on present evidence, individual-level prognosis.

## Figure legends

Supplementary Figure S19. Head-to-head benchmark of recovery predictors. (a) AUC of 22 fixed
published and candidate scores in the 84 labelled hearts after batch centring (bars) and in the
independent 19-heart cohort (points). (b) Elastic-net models restricted to prior gene lists under
cross-validation and on each held-out sequencing batch; dashes mark the 95th percentile of
size-matched random gene sets.

Supplementary Figure S20. Joint analysis and sample-size requirements. (a) Per-gene recovery effect
sizes in the two outcome-labelled cohorts; red points are the 50 genes with the smallest pooled
p-values. (b) Cross-validated AUC against training-set size with the fitted asymptotic curve.
(c) Distribution of cross-validated AUC across 20 master seeds.

Supplementary Figure S21. Robustness diagnostics. (a) Between-cohort heterogeneity of the consensus
HF signature. (b) Variance of each principal component of the 84 pre-LVAD hearts attributable to
sequencing batch and to outcome. (c) First two principal components of the 84 labelled and 134
withheld-label hearts.

## Supplementary tables

- Supplementary Table S8. Benchmark of fixed predictors (internal and external AUC, p and FDR) and of
  models restricted to prior gene lists against size-matched random sets.
- Supplementary Table S9. Joint random-effects meta-analysis of per-gene recovery effects across the
  two outcome-labelled cohorts.
- Supplementary Table S10. Learning-curve points, seed-stability runs, principal-component variance
  by batch and outcome, and inferred sex.

## Numbers to cite

| quantity | value |
|---|---|
| fixed predictors tested | 22 |
| significant internally (FDR < 0.05) | 5 |
| nominally significant externally | 2 of 22; none after FDR correction |
| significant in both cohorts | 0 |
| best prior gene list vs random sets | Hallmark union 0.72 vs random 95th pct 0.76 |
| joint meta-analysis, genes at FDR < 0.05 | 0 (smallest FDR 0.056, CD1E) |
| cross-cohort effect-size correlation | Spearman 0.01; sign agreement 50.3% (p = 0.43) |
| pathway direction concordance externally | 7 of 18 (39%), p = 0.48 |
| learning curve | 0.60 at n = 30, 0.72 at n = 84; asymptote 0.84; ~190 for 0.75, ~1,000 for 0.80 |
| seed stability (20 seeds) | 0.70 +/- 0.03 (0.65-0.76) |
| consensus signature heterogeneity | median I-squared 66%, 78% above 50% |
| leave-one-cohort-out retention | at least 81.5% |
| platform concordance (signature) | Spearman 0.73, sign agreement 100% |
| PC2 variance: batch vs outcome | 62.6% vs 7.2% |
| inferred sex | 69 male, 15 female; recovery 26% vs 53% (Fisher p = 0.06) |
| derivation vs validation shift | 41% of genes > 1 SD; PC1 separates completely (AUC 1.00) |
