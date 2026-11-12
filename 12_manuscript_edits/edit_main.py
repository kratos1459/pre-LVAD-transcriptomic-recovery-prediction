import docx, copy
from docx.enum.text import WD_ALIGN_PARAGRAPH
SRC="/mnt/user-data/uploads/022b555a-AHJ_Main_Manuscript.docx"
OUT="/home/claude/ahj_edit/AHJ_Main_Manuscript_revised.docx"
d=docx.Document(SRC)
P=d.paragraphs

def clone_after(anchor, text, bold_lead=None, heading=False):
    """Insert a new paragraph copying the anchor's formatting, after it."""
    new=copy.deepcopy(anchor._p)
    anchor._p.addnext(new)
    p=docx.text.paragraph.Paragraph(new, anchor._parent)
    for r in list(p.runs): r._element.getparent().remove(r._element)
    if bold_lead:
        r=p.add_run(bold_lead); r.bold=anchor.runs[0].bold if anchor.runs else False
    r=p.add_run(text)
    if heading: r.bold=True; r.italic=True
    return p

# ---- Methods: new subsections after para 19 (external validation methods)
a=P[19]
m2=clone_after(a, "Bulk pre-LVAD profiles were deconvolved against a reference built from the GSE226314 single-nucleus pseudobulk (13 cell types; cell types with at least 30 nuclei per sample). The 100 most specific genes per cell type (1,297 unique markers) formed the signature matrix, and fractions were estimated by non-negative least squares and normalized to sum to one. Fractions were compared between responders and non-responders, between sequencing batches, and between the consensus subtypes; used as the sole predictors of recovery; and removed from the expression matrix inside each training fold (ordinary least squares on the fractions, fitted in training data only) before the main model was refitted. In parallel, every heart was scored by single-sample gene set enrichment analysis [29] with GSEApy [30] on MSigDB Hallmark [31] and TRRUST transcription-factor target sets [32] (329 sets passing size filters), and the same leakage-controlled pipeline, held-out-batch tests, and frozen external test were applied to these rank-based set scores.")
h2=clone_after(a, "Cell-type composition and pathway-level modeling", heading=True)
m1=clone_after(a, "Repeat-averaged out-of-fold predicted probabilities were used to compute calibration-in-the-large (the intercept of a logistic model with the predicted logit as offset), the calibration slope, the Brier score, the Brier skill score relative to the observed prevalence, the expected calibration error across quartiles of predicted risk, and the Spiegelhalter z-test. Clinical utility was summarized by decision-curve analysis [26] across threshold probabilities of 0.02 to 0.75, comparing the model with treating all patients as predicted responders and treating none. The same metrics were computed for the frozen external test. The precision of that test was quantified with the Hanley-McNeil variance of the c-statistic [27]; because a permutation test on the c-statistic is equivalent to the exact Wilcoxon rank-sum test, power was also obtained by simulation from a binormal model with that exact test (20,000 replicates per scenario). Sample sizes required to estimate a c-statistic to a pre-specified confidence-interval width were derived from the same variance [28]. Predictions for the 134 validation hearts with withheld labels were generated once with the frozen model, hashed, and deposited with a pre-specified analysis plan (Supplementary Table S5).")
h1=clone_after(a, "Calibration, clinical utility, and precision", heading=True)

# ---- Results after para 28 (external replication)
a=P[28]
r2=clone_after(a, "With 10 recovered and 9 non-recovered hearts, the external cohort had 44% power to detect a true c-statistic of 0.75 and 29% power at 0.70 (exact-permutation power 46% and 30%); the smallest c-statistic detectable with 80% power was 0.85. The external result is therefore compatible both with no discrimination and with moderate discrimination. Estimating a c-statistic of 0.75 with a 95% confidence interval no wider than 0.10 would require approximately 484 patients (about 150 responders) at the derivation prevalence [28]. To allow prospective testing, probabilities for the 134 withheld-label validation hearts were frozen and deposited with their hash and analysis plan; that cohort would have approximately 97% power at a true c-statistic of 0.70 (Supplementary Table S5).")
r1=clone_after(a, "Internally, out-of-fold predictions were well calibrated in the large (calibration-in-the-large 0.11, 95% CI -0.38 to 0.59) with a calibration slope of 1.51 (95% CI 0.63 to 2.39), a Brier score of 0.178 (Brier skill score 0.17 relative to the observed prevalence), an expected calibration error of 0.06, and a non-significant Spiegelhalter test (p = 0.64). Decision curves showed net benefit above both default strategies for threshold probabilities above about 0.12, with a net benefit of 0.122 at the prevalence threshold of 0.31, equivalent to 273 fewer patients per 1,000 incorrectly classified as responders than treating all as responders. In the external cohort calibration was lost: the calibration slope fell to 0.07 (95% CI -1.45 to 1.58), the Brier skill score was negative (-0.09, worse than predicting the observed prevalence for every patient), and the model gave no net benefit over treating all patients as responders below a threshold probability of 0.45 (Supplementary Figure S16).")
h3=clone_after(a, "Calibration and clinical utility", heading=True)

# ---- Results after para 30 (subtypes/reversal)
a=P[30]
r4=clone_after(a, "Scoring the same hearts on 329 Hallmark and transcription-factor target sets gave a cross-validated AUC of 0.73 +/- 0.05, essentially unchanged without batch adjustment (0.73), and held-out-batch AUCs of 0.84 and 0.84 against 0.70 and 0.82 for the gene-level model in the same pipeline. Eighteen sets reached FDR < 0.05 (versus one gene at gene level), all lower in responders, led by NFE2L2 and ATF3 target sets, the p53 pathway, and the UV-response-up set. Applied once to the independent cohort, however, the set-level model again failed to replicate (AUC 0.54, 95% CI 0.25-0.82, exact p = 0.78; Supplementary Figure S18).")
r3=clone_after(a, "Deconvolution estimated a mean cardiomyocyte fraction of 0.61, with fibroblasts 0.12 and smooth muscle 0.12. No cell type differed between responders and non-responders after correction (smallest FDR 0.26; pericytes nominally higher in responders, 0.034 versus 0.019, p = 0.02), but composition differed by sequencing batch (smooth muscle p < 0.001, pericytes p < 0.001, epicardium p < 0.001) and between the consensus subtypes, which separated on fibroblast (0.14 versus 0.05) and cardiomyocyte (0.52 versus 0.64) content (both FDR < 0.001), confirming that the subtypes track regional sampling composition rather than a recovery-related program. Cell fractions alone predicted recovery with a cross-validated AUC of 0.61 +/- 0.04, and removing the composition component inside each fold reduced the transcriptome model by 0.04 AUC, indicating that part, but not all, of the internal signal reflects cellular composition (Supplementary Figure S17).")
h4=clone_after(a, "Cellular composition and pathway-level models", heading=True)

# ---- Discussion after para 33
a=P[33]
clone_after(a, "Two further analyses sharpen this interpretation. Composition is a real contributor: cell fractions alone discriminated responders at an AUC of 0.61, and removing them cost the full model about 0.04 AUC, so a classifier built on apical core tissue partly measures what was sampled rather than how the myocardium will behave. Moving from genes to pathway and transcription-factor target activities made the internal signal markedly more portable across sequencing batches (held-out-batch AUC 0.84 for both batches) yet left external replication unchanged, which argues that the barrier is the cohort and the outcome definition rather than the representation of the data. Internally the model was well calibrated and gave net benefit across clinically plausible thresholds; externally both calibration and net benefit were lost.")

# ---- Limitations sentence
for p in d.paragraphs:
    if "Calibration and decision-curve analysis were not performed." in p.text:
        for r in p.runs:
            if "Calibration and decision-curve analysis were not performed." in r.text:
                r.text=r.text.replace("Calibration and decision-curve analysis were not performed.",
                "The external cohort was also underpowered by design: it could not have detected anything below a c-statistic of 0.85 with conventional power, so the failure to replicate constrains but does not exclude a moderate effect. Deconvolution and pathway scoring are model-based and inherit the assumptions of their reference data.")
        break

# ---- Abstract results clause
for p in d.paragraphs:
    if p.text.startswith("Results:") and "607-gene" in p.text:
        p.runs[-1].text = p.runs[-1].text.rstrip()
        p.add_run(" Internally the model was well calibrated (slope 1.51) and gave net benefit, but calibration and net benefit were lost externally, where power to detect an AUC of 0.75 was 44%. Cell-type composition alone predicted recovery (AUC 0.61), and pathway-level models were batch-robust (held-out-batch AUC 0.84) yet still did not replicate externally (AUC 0.54).")
        break

# ---- Data and code availability
for p in d.paragraphs:
    if p.text.startswith("Data and code availability."):
        p.runs[-1].text=p.runs[-1].text.rstrip()
        p.add_run(" Frozen predictions for the 134 withheld-label validation hearts, with their SHA-256 hash and pre-specified analysis plan, and the code for the calibration, power, deconvolution, and pathway analyses are in the same repository.")
        break

# ---- References 26-29 appended after ref 25
ref25=None
for p in d.paragraphs:
    if p.text.strip().startswith("25. Al Arian T."): ref25=p
newrefs=[
 "26. Vickers AJ, Elkin EB. Decision curve analysis: a novel method for evaluating prediction models. Med Decis Making 2006;26:565-74.",
 "27. Hanley JA, McNeil BJ. The meaning and use of the area under a receiver operating characteristic (ROC) curve. Radiology 1982;143:29-36.",
 "28. Riley RD, Debray TPA, Collins GS, et al. Minimum sample size for external validation of a clinical prediction model with a binary outcome. Stat Med 2021;40:4230-51.",
 "29. Barbie DA, Tamayo P, Boehm JS, et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature 2009;462:108-12.",
 "30. Fang Z, Liu X, Peltz G. GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. Bioinformatics 2023;39:btac757.",
 "31. Liberzon A, Birger C, Thorvaldsdottir H, et al. The Molecular Signatures Database hallmark gene set collection. Cell Syst 2015;1:417-25.",
 "32. Han H, Cho JW, Lee S, et al. TRRUST v2: an expanded reference database of human and mouse transcriptional regulatory interactions. Nucleic Acids Res 2018;46:D380-6.",
]
# 23-25 inserted before 26 -> they must come after 22 but refs are numbered to 25 already; append 23..25 as 23/24/25? no: renumber-free approach
prev=ref25
for t in newrefs:
    prev=clone_after(prev, t)
d.save(OUT)
print("saved", OUT)
