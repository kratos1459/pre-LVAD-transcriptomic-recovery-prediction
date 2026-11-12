import docx, copy
from docx.shared import Inches
from docx.text.paragraph import Paragraph
SRC="/mnt/user-data/uploads/3e662ad9-AHJ_Supplementary_Material.docx"
OUT="/home/claude/ahj_edit/AHJ_Supplementary_Material_revised.docx"
IMG="/mnt/user-data/uploads/Samiran_PhD/_stage_lvad"
d=docx.Document(SRC)

def blank_from(template):
    new=copy.deepcopy(template._p)
    p=Paragraph(new, template._parent)
    for r in list(p.runs): r._element.getparent().remove(r._element)
    for dr in new.findall('.//{http://schemas.openxmlformats.org/drawingml/2006/main}graphic'):
        gp=dr.getparent()
        while gp is not None and gp.tag.split('}')[-1] not in ('r','p'): gp=gp.getparent()
        if gp is not None and gp.tag.split('}')[-1]=='r': gp.getparent().remove(gp)
    return p

def ins_after(anchor, template, text="", bold=False, italic=False):
    p=blank_from(template); anchor._p.addnext(p._p)
    if text:
        r=p.add_run(text); r.bold=bold; r.italic=italic
    return p

def ins_before(target, template, text="", bold=False, italic=False):
    p=blank_from(template); target._p.addprevious(p._p)
    if text:
        r=p.add_run(text); r.bold=bold; r.italic=italic
    return p

ps=d.paragraphs
body_tpl=[p for p in ps if p.text.strip().startswith("Harmonization.")][0]
cap_tpl=[p for p in ps if p.text.strip().startswith("Supplementary Figure S15.")][0]
figs_hdr=[p for p in ps if p.text.strip()=="Supplementary Figures"][0]
refs_hdr=[p for p in ps if p.text.strip()=="Supplementary References"][0]
tbl_tpl=[p for p in ps if p.text.strip().startswith("Supplementary Table S4.")][0]

# ---- Supplementary Methods (after the Harmonization paragraph block's section: append after last methods paragraph)
meth_last=[p for p in ps if p.text.strip().startswith("For GSE285315 snATAC-seq")][0]
m1=ins_after(meth_last, body_tpl, "Calibration, clinical utility, precision and sealed predictions. Out-of-fold probabilities were averaged over 10 repeats of 5-fold cross-validation. Calibration-in-the-large was the intercept of a logistic model of outcome on the predicted logit used as an offset, and the calibration slope its coefficient in an unconstrained logistic model; the expected calibration error was computed over quartiles of predicted risk. Decision curves report net benefit, TP/n - (FP/n)(pt/(1-pt)), across threshold probabilities pt of 0.02 to 0.75. Power for the external cohort used the Hanley-McNeil variance of the c-statistic and, independently, simulation from a binormal model scored with the exact Wilcoxon rank-sum test (20,000 replicates), which is equivalent to the permutation test used in the main analysis. Minimum sample sizes target a 95% confidence-interval width for the c-statistic. Predictions for the 134 withheld-label validation hearts were produced once with the model fitted on all 84 labeled hearts; because those hearts belong to sequencing runs outside both derivation batches, no fitted batch term applies to them, so the primary deposited score standardizes within the validation cohort and rescales to the derivation moments, as in the frozen external test (Spearman correlation with the unrescaled score 0.99).")
m2=ins_after(m1, body_tpl, "Cell-type deconvolution and pathway-level modeling. A reference signature matrix was built from the GSE226314 single-nucleus pseudobulk, averaging values back-transformed to linear scale across samples for each of 13 cell types, restricted to sample-cell-type combinations with at least 30 nuclei. The 100 most specific genes per cell type (1,297 unique markers, column-normalized) formed the signature matrix, and fractions for all 218 GSE328059 pre-LVAD hearts were estimated by non-negative least squares and scaled to sum to one. Composition was compared between responders and non-responders, between sequencing batches and between consensus subtypes with Welch tests and Benjamini-Hochberg correction. Recovery was then predicted from fractions alone, and the 6,000-gene model was refitted after removing the composition component inside each training fold (ordinary least squares of expression on fractions, fitted in training data only, subtracted from training and test data). Separately, single-sample gene set enrichment analysis was run with GSEApy (rank normalization) on the MSigDB Hallmark 2020 and TRRUST 2019 transcription-factor target libraries, keeping sets of 10 to 500 genes (329 sets scored). The stage-05 pipeline (in-fold batch adjustment, F-test selection of 40 sets, elastic net, 10 repeats of 5-fold cross-validation), the held-out-batch tests and the frozen external test were repeated on these scores.")

# ---- Table notes S5-S7 immediately before the Supplementary Figures heading
for t in ["Supplementary Table S5. Frozen predictions for the 134 GSE328059 validation hearts with withheld labels (sample, GSM, sequencing batch, primary and unrescaled predicted probability, predicted class), with the pre-specified analysis plan and the SHA-256 hash of the prediction file (b139bb3fb5cdcd7db703df2ffe1f7aa8634be61212d34341af4ace494685c5b4). The pre-specified primary analysis is the c-statistic against the released labels with a 5,000-bootstrap confidence interval and a two-sided permutation p-value, with no refitting, rescaling or threshold change.",
          "Supplementary Table S6. Estimated cell-type fractions for all 218 GSE328059 pre-LVAD hearts and their association with recovery, sequencing batch and consensus subtype. Provided as deconvolution_fractions.csv and deconvolution_associations.csv in the code repository.",
          "Supplementary Table S7. Single-sample GSEA scores for 329 Hallmark and TRRUST sets in the 84 labeled pre-LVAD hearts and their univariate association with recovery (Welch t, p, Benjamini-Hochberg FDR). Eighteen sets reach FDR < 0.05, all lower in responders. Provided as pathway_scores_GSE328059.csv and pathway_univariate.csv in the code repository."]:
    ins_before(figs_hdr, tbl_tpl, t)

# ---- Figures S16-S18 before Supplementary References
specs=[("FigS16_calibration_dca.png","Supplementary Figure S16. Calibration and clinical utility of the pre-LVAD recovery model. (a) Calibration of repeat-averaged out-of-fold predictions in the 84 labeled hearts; points are quartiles of predicted risk with 95% confidence intervals, the red line is the logistic recalibration curve, the dashed line is perfect calibration. (b) Decision curve for the same predictions against treating all patients as predicted responders (grey dashed) and treating none (black). (c, d) The same analyses for the frozen external test in PRJNA1312133 (19 hearts)."),
       ("FigS17_deconvolution.png","Supplementary Figure S17. Cell-type composition of the pre-LVAD hearts. (a) Estimated fractions of the eight most abundant cell types in responders and non-responders (box, median and interquartile range; whiskers, 1.5 x IQR; outliers omitted). (b) Composition of the two consensus subtypes. (c) Cross-validated AUC for recovery prediction from cell fractions alone, from the 6,000-gene transcriptome, and from the transcriptome after in-fold removal of the composition component (mean and SD over ten repeats)."),
       ("FigS18_pathway_activity.png","Supplementary Figure S18. Pathway- and transcription-factor-level recovery modeling. (a) The ten Hallmark or TRRUST sets most strongly associated with recovery; blue, higher in responders; red, lower. (b) AUC of the set-level model under cross-validation, on each held-out sequencing batch, and in the frozen external test. (c) Distribution of out-of-fold scores by outcome.")]
for fn, cap in specs:
    pic=ins_before(refs_hdr, cap_tpl)
    pic.add_run().add_picture(f"{IMG}/{fn}", width=Inches(6.0))
    ins_before(refs_hdr, cap_tpl, cap)
d.save(OUT)
print("saved", OUT)
