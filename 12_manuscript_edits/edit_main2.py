import docx, copy
SRC="/home/claude/ahj_edit/AHJ_Main_Manuscript_revised.docx"
OUT="/home/claude/ahj_edit/AHJ_Main_Manuscript_revised_v2.docx"
d=docx.Document(SRC)

def find(sub):
    for p in d.paragraphs:
        if sub in p.text: return p
    raise KeyError(sub)

def after(anchor, text, heading=False):
    new=copy.deepcopy(anchor._p); anchor._p.addnext(new)
    p=docx.text.paragraph.Paragraph(new, anchor._parent)
    for r in list(p.runs): r._element.getparent().remove(r._element)
    r=p.add_run(text)
    if heading: r.bold=True; r.italic=True
    return p

# ---- Methods: new subsection after the composition/pathway methods paragraph
a=find("the same leakage-controlled pipeline, held-out-batch tests, and frozen external test were applied")
h=after(a, "Benchmark of published predictors, joint analysis, and robustness", heading=True)
after(h, "Published and candidate recovery predictors were put through the same leakage-controlled frame. Single genes (LRRN4CL, NPPB, ATF3, CAMK2D, PCDH9, MYL7, NPPA, MYH7) and signed scores for the consensus HF signature and for Hallmark and TRRUST sets were evaluated as fixed scores in the 84 labeled hearts, with and without batch centering, and once in PRJNA1312133; the same prior gene lists were also used to restrict feature selection, so that each list produced its own elastic-net model, which was compared with size-matched random gene sets. Per-gene recovery effect sizes were pooled across the two outcome-labeled cohorts by DerSimonian-Laird random effects. A learning curve was obtained by subsampling the 84 hearts (10 draws per size) and fitting AUC = a - c/sqrt(n), and the whole pipeline was repeated under 20 master seeds. Heterogeneity of the consensus signature was summarized by I2 and tau2, by leave-one-cohort-out retention, and by microarray versus RNA-seq concordance; Egger tests were not performed because at most eight cohorts contribute to any gene. The variance of each principal component attributable to batch and to outcome was computed by one-way analysis of variance. Sex was inferred from XIST against the mean of RPS4Y1, DDX3Y, UTY, KDM5D, and EIF1AY. Comparability of the 84 labeled and 134 withheld-label hearts was assessed by standardized per-gene mean shifts and by principal components.")

# ---- Results: new subsections after the pathway results paragraph
a=find("the set-level model again failed to replicate (AUC 0.54")
h=after(a, "No published predictor replicates, and prior gene lists do not beat chance", heading=True)
p1=after(h, "Of 22 fixed published and candidate predictors, six were associated with recovery in the 84 labeled hearts (ATF3 target set AUC 0.73, ATF3 itself 0.71, NFE2L2 targets 0.66, LRRN4CL 0.64; five at FDR < 0.05), and two reached nominal significance in the independent cohort (the consensus HF signature restricted to its top 150 genes, AUC 0.82, p = 0.017; Hallmark oxidative phosphorylation, AUC 0.79, p = 0.035). None survived correction externally, and no predictor was significant in both cohorts. When prior gene lists were instead used to restrict feature selection, none outperformed size-matched random gene sets: the Hallmark union reached a cross-validated AUC of 0.72 against a random-set 95th percentile of 0.76, and the consensus HF signature 0.70 against 0.78 (Supplementary Figure S19, Supplementary Table S8).")
p2=after(p1, "Pooling per-gene recovery effects across GSE328059 and PRJNA1312133 by random effects left no gene at FDR < 0.05 (smallest FDR 0.056, CD1E). Effect sizes were uncorrelated between the cohorts (Spearman r = 0.01) and agreed in direction for 50.3% of genes, indistinguishable from chance (binomial p = 0.43). The 18 sets that reached FDR < 0.05 internally were directionally concordant in the external cohort for 7 of 18 (39%, p = 0.48), and none reached p < 0.05 there (Supplementary Table S9).")
h2=after(p2, "Sample-size requirements and robustness", heading=True)
p3=after(h2, "Cross-validated AUC rose from 0.60 at n = 30 to 0.72 at n = 84; fitting AUC = a - c/sqrt(n) gave an asymptote of 0.84 and implied that about 190 labeled hearts would be needed for a cross-validated AUC of 0.75, and about 1,000 for 0.80. Repeating the entire pipeline under 20 master seeds gave 0.70 +/- 0.03 (range 0.65 to 0.76), so single-run estimates carry a seed-level uncertainty of about 0.05 AUC (Supplementary Figure S20).")
after(p3, "Among the consensus HF genes the median I2 was 66%, with 78% above 50%, yet direction was unanimous between platforms (sign agreement 100%; Spearman r = 0.73 between microarray and RNA-seq pooled effects), and leave-one-cohort-out removal retained at least 81.5% of the signature whichever cohort was dropped. In the 84 pre-LVAD hearts, principal component 2 was 62.6% attributable to sequencing batch and 7.2% to outcome. Inferred sex gave 69 male and 15 female hearts, with recovery in 26% of male and 53% of female hearts (Fisher p = 0.06); adding sex to the model did not change performance, and the model was weaker in men alone (AUC 0.63). The 134 withheld-label validation hearts were not exchangeable with the 84 labeled hearts: 41% of genes shifted by more than one pooled standard deviation, and the first principal component separated the two groups completely (Supplementary Figure S21).")

# ---- Discussion
a=find("Two further analyses sharpen this interpretation")
after(a, "The predictor benchmark is the strongest negative evidence here: across 22 published and candidate predictors, prior gene lists that performed no better than random sets of the same size, and a joint analysis in which per-gene effects were uncorrelated between the two labeled cohorts, nothing about the pre-LVAD transcriptome transferred between cohorts. The learning curve suggests that this is partly a sample-size problem, with roughly 190 labeled hearts needed for a cross-validated AUC of 0.75, but the complete absence of cross-cohort agreement at gene and set level, in cohorts that differ in platform and outcome ascertainment, points to cohort heterogeneity as the more likely barrier. The disease signature itself was unanimous in direction across platforms, so these data support robust descriptive biology and not, on present evidence, individual-level prognosis.")

# ---- Limitations addition
a=find("Deconvolution and pathway scoring are model-based")
for r in a.runs:
    if "Deconvolution and pathway scoring are model-based" in r.text:
        r.text=r.text.replace("Deconvolution and pathway scoring are model-based and inherit the assumptions of their reference data.",
        "Deconvolution and pathway scoring are model-based and inherit the assumptions of their reference data. The 134 withheld-label hearts are not exchangeable with the 84 labeled hearts, so the deposited predictions for them are a test of transportability as much as of the model. Sex was inferred from expression rather than reported, and the 20-seed spread of about 0.05 AUC should be borne in mind when comparing any two numbers in this study.")

# ---- Abstract
for p in d.paragraphs:
    if p.text.startswith("Results:"):
        p.runs[-1].text=p.runs[-1].text.rstrip()
        p.add_run(" In a head-to-head benchmark, none of 22 published or candidate predictors was significant in both labeled cohorts, prior gene lists did not beat size-matched random sets, and pooled per-gene effects were uncorrelated across cohorts (Spearman r = 0.01).")
        break
d.save(OUT); print("saved", OUT)
