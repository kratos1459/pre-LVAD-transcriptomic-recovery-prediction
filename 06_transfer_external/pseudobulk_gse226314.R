# Stage 06c: GSE226314 snRNA-seq -> patient-level and patient x cell-type pseudobulk log2 CPM.
# Input: GSE226314_global.rds (GEO supplementary, Seurat object). Needs ~64 GB RAM.
# Usage: Rscript pseudobulk_gse226314.R GSE226314_global.rds ../results/cohorts
suppressMessages({library(Seurat); library(Matrix)})
a <- commandArgs(TRUE); o <- readRDS(a[1]); out <- a[2]
md <- o@meta.data
sample_col <- intersect(c("orig.ident", "sample", "Sample"), colnames(md))[1]
ct_col <- intersect(c("celltype", "cell_type", "CellType", "seurat_annotations", "annotation"), colnames(md))[1]
cnt <- GetAssayData(o, assay = DefaultAssay(o), layer = "counts")
cpm_log <- function(m) { m <- as.matrix(m); log2(t(t(m) / colSums(m) * 1e6) + 1) }
pb <- function(group) { g <- factor(group); M <- sparse.model.matrix(~ 0 + g); colnames(M) <- levels(g); cnt %*% M }
all_pb <- cpm_log(pb(md[[sample_col]]))
write.csv(all_pb, file.path(out, "GSE226314_pseudobulk_expr.csv"))
ct_pb <- cpm_log(pb(paste(md[[sample_col]], md[[ct_col]], sep = "__")))
write.csv(ct_pb, file.path(out, "GSE226314_pseudobulk_by_celltype.csv"))
write.csv(unique(md[, c(sample_col, intersect("condition", colnames(md)))]), file.path(out, "GSE226314_sample_condition.csv"), row.names = FALSE)
# Condition labels (donor / responder pre-post / non-responder pre-post) come from GEO sample characteristics ("condition: ...").
# Gzip the pseudobulk and add rows to all_meta.csv (group HF = *pre, NF = Donor) so that 02_meta_lasso_shap picks up the 8th cohort.
