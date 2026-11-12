import docx, copy
from docx.shared import Inches
from docx.text.paragraph import Paragraph
SRC="/home/claude/ahj_edit/AHJ_Supplementary_Material_revised.docx"
OUT="/home/claude/ahj_edit/AHJ_Supplementary_Material_revised_v2.docx"
IMG="/mnt/user-data/uploads/Samiran_PhD/_stage_lvad"
d=docx.Document(SRC)

def blank_from(t):
    new=copy.deepcopy(t._p); p=Paragraph(new, t._parent)
    for r in list(p.runs): r._element.getparent().remove(r._element)
    for g in new.findall('.//{http://schemas.openxmlformats.org/drawingml/2006/main}graphic'):
        gp=g.getparent()
        while gp is not None and gp.tag.split('}')[-1] not in ('r','p'): gp=gp.getparent()
        if gp is not None and gp.tag.split('}')[-1]=='r': gp.getparent().remove(gp)
    return p

def after(anchor, tpl, text=""):
    p=blank_from(tpl); anchor._p.addnext(p._p)
    if text: p.add_run(text)
    return p

def before(target, tpl, text=""):
    p=blank_from(tpl); target._p.addprevious(p._p)
    if text: p.add_run(text)
    return p

def find(sub):
    for p in d.paragraphs:
        if sub in p.text: return p
    raise KeyError(sub)

body=find("Cell-type deconvolution and pathway-level modeling.")
cap=find("Supplementary Figure S18.")
tbl=find("Supplementary Table S7.")
refs=find("Supplementary References")

after(body, body, "Benchmark, joint analysis and robustness diagnostics. Fixed predictors (single genes, the consensus HF signature, Hallmark and TRRUST set scores) were scored as within-cohort z-scores, evaluated by AUC with and without per-batch centering, and tested by Mann-Whitney; the same prior gene lists were used to restrict feature selection for elastic-net models, each compared with 20 size-matched random gene sets drawn from the 6,000 variable genes. Per-gene recovery effects were computed as Hedges g with its standard error in each labeled cohort and pooled by DerSimonian-Laird random effects with Benjamini-Hochberg correction. Learning-curve points used stratified subsamples (10 draws per size, 3 cross-validation repeats each) with the asymptotic fit AUC = a - c/sqrt(n). Seed stability repeated the pipeline under 20 master seeds. Principal components were computed on the standardized 6,000-gene matrix, with the variance of each component attributed to batch and to outcome by one-way analysis of variance. Sex was inferred by thresholding the difference between the mean of RPS4Y1, DDX3Y, UTY, KDM5D and EIF1AY and XIST. Derivation-versus-validation comparability used standardized per-gene mean differences and the first two principal components of the combined 218 hearts.")

for t in ["Supplementary Table S8. Head-to-head benchmark of recovery predictors: AUC, p and FDR for 22 fixed scores in the 84 labeled hearts and in PRJNA1312133, and cross-validated and held-out-batch AUC for elastic-net models restricted to prior gene lists against size-matched random sets. Provided as benchmark_fixed_scores.csv and benchmark_fitted_models.csv.",
          "Supplementary Table S9. Joint random-effects meta-analysis of per-gene recovery effects across GSE328059 and PRJNA1312133 (Hedges g per cohort, pooled g, tau2, p, FDR, sign concordance), and directional concordance in the external cohort of the 18 gene sets significant internally. Provided as joint_recovery_meta.csv and pathway_external_concordance.csv.",
          "Supplementary Table S10. Learning-curve points, principal-component variance attributable to batch and outcome, and inferred sex with recovery status. Provided as learning_curve.csv, pc_variance_batch_outcome.csv and inferred_sex.csv."]:
    before(refs, tbl, t)

specs=[("FigS19_signature_benchmark.png","Supplementary Figure S19. Head-to-head benchmark of recovery predictors. (a) AUC of 22 fixed published and candidate scores in the 84 labeled hearts after batch centering (bars) and in the independent 19-heart cohort (points); the dashed line marks chance. (b) Elastic-net models restricted to prior gene lists under cross-validation and on each held-out sequencing batch; dashes mark the 95th percentile of size-matched random gene sets."),
       ("FigS20_joint_meta_learning.png","Supplementary Figure S20. Joint analysis and sample-size requirements. (a) Per-gene recovery effect sizes in the two outcome-labeled cohorts; red points are the 50 genes with the smallest pooled p-values. (b) Cross-validated AUC against training-set size with the fitted asymptotic curve; the dotted line is the fitted asymptote. (c) Distribution of cross-validated AUC across 20 master seeds."),
       ("FigS21_diagnostics.png","Supplementary Figure S21. Robustness diagnostics. (a) Between-cohort heterogeneity of the consensus HF signature; the dashed line is the median. (b) Variance of each principal component of the 84 pre-LVAD hearts attributable to sequencing batch and to outcome. (c) First two principal components of the 84 labeled and the 134 withheld-label hearts.")]
for fn, c in specs:
    pic=before(refs, cap)
    pic.add_run().add_picture(f"{IMG}/{fn}", width=Inches(6.0))
    before(refs, cap, c)
d.save(OUT); print("saved", OUT)
