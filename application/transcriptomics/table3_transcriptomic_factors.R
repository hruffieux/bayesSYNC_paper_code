## ---------------------------------------------------------------------------
## Summary table of the RNA-seq factors, for the main text: one row per factor,
## computed from the output_slim.RData of a analysis_transcriptomics.R run.
##
## Every column rests on the loadings or on patient-level summaries of the
## fitted factor trajectories, not on the shapes of the eigenfunctions, so the
## table does not depend on the edge behaviour of the trajectories.
##
##   share selected   proportion of genes with PPI > 0.5
##   key genes        largest positive and negative loadings
##   top-100 length   percentile of the median transcript length of the 100
##                    largest loadings among the analysed genes; extreme values
##                    point to a technical axis (long or short transcripts)
##   severity, hsCRP  rank correlation of the patient-level factor value with the
##                    severity class (B < C < D < E) and with the peak hsCRP, and
##                    the Kruskal-Wallis p-value for the class
##   batch            Kruskal-Wallis p-value for the extraction batch of the
##                    patient's first sample (as in section 12c of the script)
##   length bias      rank correlation with the patient-level length-bias index
##                    (as in section 12d)
##   enrichment       hallmark sets, by GSEA on the ranked loading magnitudes and
##                    by over-representation of the n_top_ora genes with the
##                    largest loadings, against the analysed genes
##
## The patient-level factor value is the mean of the fitted factor trajectory at
## the patient's own sampling times. It uses no extrapolated part of the
## trajectory and does not depend on the rotation within the factor made by the
## post-processing. Each factor is oriented so that this value correlates
## positively with the severity class, so that, for a factor associated with
## severity, genes with positive loadings are higher in the more severe patients.
##
## Writes table_rnaseq_factors.csv (all columns) and table_rnaseq_factors.tex (a
## booktabs fragment of the main columns) next to output_slim.RData. Labels are
## added by hand in labels_manual once the factors have been read.
##
## Usage, from the root of the repository:
##   Rscript application/transcriptomics/table3_transcriptomic_factors.R <run_dir>
## ---------------------------------------------------------------------------

rm(list = ls())


suppressPackageStartupMessages({
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(msigdbr)
})

# run directory, relative to main_dir, or given on the command line
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1, dir.exists(args[1]))
run_dir <- args[1]

# labels given by hand after reading the key genes and the enrichment, by factor
# number; NA leaves the label empty. Factor numbers are specific to a fit: these
# are the labels of Table 3, for the fit of the paper to the real data (N = 164,
# K = 4, CV0047 excluded). They are not applied to any other fit, in particular
# not to a fit to the simulated data.
labels_manual <- c(
  "1"  = "Neutrophil/lymphocyte",
  "3"  = "5S rRNA, pseudogenes",
  "6"  = "Interferon",
  "9"  = "Short transcripts, pseudogenes",
  "10" = "Erythroid")

thres_ppi <- 0.5
n_key <- 4        # key genes shown per sign
n_top_length <- 100
n_top_ora <- 200
min_pathway_size <- 5
max_pathway_size <- 200
fdr_show <- 0.05
seed <- 1


## ---------------------------------------------------------------------------
## Loading
## ---------------------------------------------------------------------------

e <- new.env()
load(file.path(run_dir, "output_slim.RData"), envir = e)
if (isTRUE(e$bool_simulated_data) || e$N != 164) {
  labels_manual <- character(0)
}
res <- e$res_slim
time_obs <- e$time_obs
df_samples <- e$df_comb_no_HC
subject_ids <- names(time_obs)
Q <- ncol(res$B_hat)
genes <- rownames(res$B_hat)
stopifnot(identical(names(res$list_h_hat), subject_ids) || all(subject_ids %in% names(res$list_h_hat)))


## ---------------------------------------------------------------------------
## Patient-level covariates
## ---------------------------------------------------------------------------

df_first <- df_samples[order(df_samples$days_from_sx_or_swab_imputed), ]
df_first <- df_first[!duplicated(df_first$subject_id), ]
severity <- as.character(df_first$severity[match(subject_ids, df_first$subject_id)])
sev_rank <- match(severity, c("B", "C", "D", "E"))
batch_first <- droplevels(as.factor(df_first$extractbatch[match(subject_ids, df_first$subject_id)]))
crp_peak <- tapply(df_samples$hsCRP_mg_l, df_samples$subject_id,
                   function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))[subject_ids]

# transcript lengths and per-sample length-bias index, as in section 12d
gene_length <- e$df_genes_sel[e$gene_id_of_label[genes], "Length"]
names(gene_length) <- genes
log_length <- log10(gene_length)
E_analysed <- t(as.matrix(df_samples[, genes, drop = FALSE]))
dev_mat <- E_analysed - rowMeans(E_analysed)
lb <- log_length - mean(log_length, na.rm = TRUE)
ok <- !is.na(lb)
length_bias <- colSums(dev_mat[ok, ] * lb[ok]) / sqrt(colSums(dev_mat[ok, ]^2) * sum(lb[ok]^2))
lb_patient <- tapply(length_bias, df_samples$subject_id, mean)[subject_ids]


## ---------------------------------------------------------------------------
## Patient-level factor values: mean of the fitted trajectory at the patient's
## own sampling times
## ---------------------------------------------------------------------------

H_patient <- sapply(seq_len(Q), function(q) {
  sapply(subject_ids, function(s) {
    mean(approx(res$time_g, res$list_h_hat[[s]][[q]], xout = time_obs[[s]], rule = 2)$y)
  })
})
colnames(H_patient) <- paste0("F", seq_len(Q))

# orientation: positive association with the severity class
sgn <- sign(sapply(seq_len(Q), function(q) cor(H_patient[, q], sev_rank, method = "spearman")))
sgn[sgn == 0] <- 1
H_patient <- sweep(H_patient, 2, sgn, "*")
B <- sweep(res$B_hat, 2, sgn, "*")


## ---------------------------------------------------------------------------
## Enrichment, with the de-duplication of section 14
## ---------------------------------------------------------------------------

df_map <- suppressMessages(bitr(e$selected_genes, fromType = "ENSEMBL",
                                toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db))
df_map <- df_map[!duplicated(df_map$ENSEMBL), , drop = FALSE]
df_map <- df_map[!is.na(df_map$ENTREZID) & !duplicated(df_map$ENTREZID), , drop = FALSE]
df_map$label <- e$gene_label_of_id[df_map$ENSEMBL]
df_map <- df_map[df_map$label %in% genes, , drop = FALSE]
rownames(df_map) <- df_map$label
cat("Genes mapped to Entrez identifiers: ", nrow(df_map), " / ", length(genes), "\n", sep = "")

df_t2g <- suppressWarnings(msigdbr(species = "Homo sapiens", category = "H"))
df_t2g <- data.frame(term = sub("^HALLMARK_", "", df_t2g$gs_name), gene = as.character(df_t2g$entrez_gene))

top_terms <- function(df, col_q = "p.adjust") {
  if (is.null(df) || nrow(df) == 0) return("")
  df <- df[df[[col_q]] < fdr_show, , drop = FALSE]
  if (nrow(df) == 0) return("none")
  df <- df[order(df[[col_q]]), , drop = FALSE]
  paste(sprintf("%s (%.1e)", tolower(df$ID[seq_len(min(2, nrow(df)))]), df[[col_q]][seq_len(min(2, nrow(df)))]),
        collapse = "; ")
}


## ---------------------------------------------------------------------------
## One row per factor
## ---------------------------------------------------------------------------

rows <- list()
for (q in seq_len(Q)) {

  b <- B[, q]
  sel <- res$ppi[, q] > thres_ppi
  o <- order(abs(b), decreasing = TRUE)
  top_len <- genes[o[seq_len(n_top_length)]]
  pct_len <- 100 * mean(gene_length < median(gene_length[top_len], na.rm = TRUE), na.rm = TRUE)
  # key genes: largest loadings of each sign among annotated genes, so that
  # Ensembl identifiers without a symbol do not fill the table
  b_annot <- b[!grepl("^ENSG", names(b))]
  key_pos <- names(sort(b_annot[b_annot > 0], decreasing = TRUE))[seq_len(min(n_key, sum(b_annot > 0)))]
  key_neg <- names(sort(b_annot[b_annot < 0]))[seq_len(min(n_key, sum(b_annot < 0)))]

  h <- H_patient[, q]
  ct_sev <- suppressWarnings(cor.test(h, sev_rank, method = "spearman", exact = FALSE))
  ct_crp <- suppressWarnings(cor.test(h, crp_peak, method = "spearman", exact = FALSE))
  ct_lb <- suppressWarnings(cor.test(h, lb_patient, method = "spearman", exact = FALSE))

  # GSEA on the loading magnitudes, as in section 14
  rank_list <- abs(b[rownames(df_map)])
  names(rank_list) <- df_map$ENTREZID
  rank_list <- sort(rank_list[!is.na(rank_list)], decreasing = TRUE)
  set.seed(seed)
  gsea <- tryCatch(as.data.frame(GSEA(geneList = rank_list, TERM2GENE = df_t2g, scoreType = "pos",
                                      minGSSize = min_pathway_size, maxGSSize = max_pathway_size,
                                      pvalueCutoff = 1, verbose = FALSE)),
                   error = function(err) NULL)

  # over-representation of the genes with the largest loadings
  target <- df_map[intersect(genes[o[seq_len(n_top_ora)]], rownames(df_map)), "ENTREZID"]
  ora <- tryCatch(as.data.frame(enricher(gene = target, universe = df_map$ENTREZID, TERM2GENE = df_t2g,
                                         minGSSize = min_pathway_size, maxGSSize = max_pathway_size,
                                         pvalueCutoff = 1)),
                  error = function(err) NULL)

  rows[[q]] <- data.frame(
    factor = q,
    label = if (as.character(q) %in% names(labels_manual)) labels_manual[[as.character(q)]] else NA,
    share_selected = mean(sel),
    key_genes_pos = paste(key_pos, collapse = ", "),
    key_genes_neg = paste(key_neg, collapse = ", "),
    top100_length_percentile = pct_len,
    rho_severity = unname(ct_sev$estimate), p_severity_kw = kruskal.test(h ~ factor(severity))$p.value,
    rho_hscrp = unname(ct_crp$estimate), p_hscrp = ct_crp$p.value,
    p_batch_kw = kruskal.test(h ~ batch_first)$p.value,
    rho_length_bias = unname(ct_lb$estimate), p_length_bias = ct_lb$p.value,
    gsea = top_terms(gsea), ora_top = top_terms(ora),
    stringsAsFactors = FALSE)
}
df_tab <- do.call(rbind, rows)

# A factor is flagged for a closer look as possibly technical when its top
# loadings sit at an extreme of the transcript-length distribution, or when it
# is very strongly associated with the extraction batch. The length-bias index is
# not used for the flag, since cell composition itself moves it (samples rich in
# neutrophils tend to have poorer RNA quality), so that biological factors can be
# correlated with it. The final call is made by hand, from the key genes.
df_tab$technical_flag <- with(df_tab, top100_length_percentile < 10 | top100_length_percentile > 90 |
                                p_batch_kw < 1e-6)

write.csv(df_tab, file.path(run_dir, "table_rnaseq_factors.csv"), row.names = FALSE)

old <- options(width = 200)
print(df_tab[, c("factor", "share_selected", "key_genes_pos", "key_genes_neg", "top100_length_percentile",
                 "rho_severity", "p_severity_kw", "rho_hscrp", "p_batch_kw", "rho_length_bias", "technical_flag")],
      digits = 2, row.names = FALSE)
cat("\nEnrichment (hallmark, FDR < ", fdr_show, "):\n", sep = "")
for (q in seq_len(Q)) cat(sprintf("F%d  GSEA: %s\n    ORA top %d: %s\n", q, df_tab$gsea[q], n_top_ora, df_tab$ora_top[q]))
options(old)


## ---------------------------------------------------------------------------
## LaTeX fragment of the main columns
## ---------------------------------------------------------------------------

fmt_p <- function(x) {
  sapply(x, function(v) {
    if (is.na(v)) return("--")
    if (v < 1e-3) {
      ex <- floor(log10(v))
      return(sprintf("$%.0f \\times 10^{%d}$", v / 10^ex, ex))
    }
    sprintf("%.3f", v)
  })
}
esc <- function(x) gsub("_", "\\\\_", x)

lines <- c("% Generated by table3_transcriptomic_factors.R. Do not edit by hand.",
           "\\begin{tabular}{llp{4.2cm}rrrl}",
           "\\toprule",
           "Factor & Label & Key genes (positive; negative) & Selected (\\%) & Severity $\\rho$ & hsCRP $\\rho$ & Enrichment (hallmark, GSEA) \\\\",
           "\\midrule")
for (q in seq_len(Q)) {
  r <- df_tab[q, ]
  lines <- c(lines, sprintf("%d & %s & %s; %s & %.0f & %.2f & %.2f & %s \\\\",
                            r$factor, ifelse(is.na(r$label), "", esc(r$label)),
                            esc(r$key_genes_pos), esc(r$key_genes_neg), 100 * r$share_selected,
                            r$rho_severity, r$rho_hscrp, esc(r$gsea)))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(run_dir, "table_rnaseq_factors.tex"))
cat("\nwritten: table_rnaseq_factors.csv and table_rnaseq_factors.tex in", run_dir, "\n")
