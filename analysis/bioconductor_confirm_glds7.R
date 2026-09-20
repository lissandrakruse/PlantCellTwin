suppressPackageStartupMessages({
  library(affy)
  library(limma)
  library(jsonlite)
  library(ath1121501cdf)
})

root <- "/workspace/scratch/1b7c2e5453f2"
site <- file.path(root, "orbital-cell-twin")
meta <- read.csv(file.path(root, "glds7_sample_metadata.csv"), stringsAsFactors=FALSE)
map <- read.csv(file.path(root, "glds7_probe_locus_map.csv"), stringsAsFactors=FALSE)
files <- list.files(file.path(root, "GSE56659_RAW"), pattern="[.]CEL[.]gz$", full.names=TRUE)
acc <- sub("_.*", "", basename(files))
files <- files[match(meta$accession, acc)]
stopifnot(!anyNA(files), all(file.exists(files)))

message("Reading 36 raw CEL files...")
raw <- ReadAffy(filenames=files, cdfname="ath1121501cdf")
message("Running canonical Bioconductor RMA...")
eset <- rma(raw, verbose=TRUE)
probe_expr <- exprs(eset)
colnames(probe_expr) <- sub("_.*", "", basename(sampleNames(eset)))
probe_expr <- probe_expr[, meta$accession, drop=FALSE]
write.csv(probe_expr, gzfile(file.path(root, "GSE56659_Bioconductor_RMA_expression.csv.gz")))

map <- map[match(rownames(probe_expr), map$probe),]
valid <- !is.na(map$locus) & nzchar(map$locus)
v <- apply(probe_expr[valid,,drop=FALSE], 1, var)
pick <- order(v, decreasing=TRUE)
candidate <- data.frame(probe=rownames(probe_expr)[valid][pick], locus=map$locus[valid][pick], variance=v[pick])
candidate <- candidate[!duplicated(candidate$locus),]
gene_expr <- probe_expr[candidate$probe,,drop=FALSE]
rownames(gene_expr) <- candidate$locus

all_results <- list()
summary_rows <- list()
for (tissue in c("whole plant", "shoot", "hypocotyl", "root")) {
  m <- meta[meta$tissue == tissue,]
  m$condition <- factor(m$condition, levels=c("ground", "flight"))
  m$plate <- factor(m$plate)
  design <- model.matrix(~ condition + plate, data=m)
  fit <- eBayes(lmFit(gene_expr[,m$accession,drop=FALSE], design))
  tt <- topTable(fit, coef="conditionflight", number=Inf, sort.by="none")
  out <- data.frame(tissue=tissue, locus=rownames(tt), log2fc=tt$logFC,
                    moderated_t=tt$t, p=tt$P.Value, fdr=tt$adj.P.Val)
  all_results[[tissue]] <- out
  summary_rows[[tissue]] <- data.frame(tissue=tissue, genes_tested=nrow(out),
    fdr_005=sum(out$fdr < 0.05), p001=sum(out$p <= 0.01))
}
results <- do.call(rbind, all_results)
rownames(results) <- NULL
write.csv(results, file.path(site, "bioconductor-limma-results.csv"), row.names=FALSE)

py <- read.csv(file.path(site, "rma-moderated-results.csv"), stringsAsFactors=FALSE)
concord <- list()
for (tissue in c("whole plant", "shoot", "hypocotyl", "root")) {
  a <- results[results$tissue == tissue,]
  b <- py[py$tissue == tissue,]
  x <- merge(a, b, by="locus", suffixes=c("_bioc", "_python"))
  ta <- head(x[order(x$p_bioc), "locus"], 300)
  tb <- head(x[order(x$p_python), "locus"], 300)
  common <- intersect(ta, tb)
  same <- if (length(common)) mean(sign(x$log2fc_bioc[match(common,x$locus)]) == sign(x$log2fc_python[match(common,x$locus)])) else NA
  concord[[tissue]] <- data.frame(tissue=tissue, shared_loci=nrow(x),
    effect_spearman=cor(x$log2fc_bioc, x$log2fc_python, method="spearman"),
    top300_overlap=length(common), top300_same_direction=same,
    bioc_fdr_005=sum(a$fdr < .05), python_fdr_005=sum(b$fdr < .05))
}
concord <- do.call(rbind, concord)
write.csv(concord, file.path(site, "bioconductor-python-concordance.csv"), row.names=FALSE)

payload <- list(
  engine=list(R=R.version.string, affy=as.character(packageVersion("affy")), limma=as.character(packageVersion("limma")), cdf=as.character(packageVersion("ath1121501cdf"))),
  workflow="ReadAffy -> affy::rma -> one representative probeset per Arabidopsis locus by label-independent overall variance -> limma tissue-specific model (~ condition + plate) -> eBayes",
  samples=ncol(probe_expr), probesets=nrow(probe_expr), loci=nrow(gene_expr),
  tissue_summary=do.call(rbind, summary_rows), concordance=concord,
  interpretation="Independent implementation check. It tests computational reproducibility, not biological validation."
)
write_json(payload, file.path(site, "bioconductor-summary.json"), pretty=TRUE, auto_unbox=TRUE, dataframe="rows", digits=8)
print(payload)
