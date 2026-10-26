## ---------------------------------------------------------------------------
## Supplementary figure of Supplementary Material S9
## (figure_S_rnaseq_temporal.pdf) for the transcriptomic analysis, drawn from the output_slim.RData of a analysis_transcriptomics.R
## run.
##
## The figure shows the temporal structure of the two factors that recur across
## the sensitivity refits, interferon and erythroid.
##
##   (a) The estimated mean functions of the two factors, each averaged over the
##       50 genes with the largest loadings on that factor. The top-50 loadings
##       of each factor share a sign, so the average is a meaningful summary of
##       how those genes move over the window. Curves are centred on their mean
##       over the plotted window so that the two are comparable. These are
##       population quantities, so the whole window can be plotted.
##
##   (b) Expression of the canonical interferon marker genes,
##       mu_j(t) + b_jq h_iq(t) averaged over the markers, by severity class,
##       over days 10 to 30 and over the patients whose own sampling spans that
##       interval, so that no patient trajectory is extrapolated. Marker genes
##       load negatively on the factor, so the direction of the latent
##       trajectory alone says nothing about the direction of expression. The
##       erythroid factor is not shown here because its association with
##       severity is not consistent across intervals (Supplementary S9).
##
## All panels are on the scale of the standardised expression data. Samples run
## from day 1 to day 49, sparsely at both ends.
##
## Factors are matched to the marker sets by the number of markers among the 200
## largest loadings, so the script works on any run without editing factor
## numbers; the match is printed.
##
## Usage, from the root of the repository:
##   Rscript application/transcriptomics/plot_figure_S_transcriptomics_temporal.R <run_dir> [out_dir]
## The figure is written in out_dir, by default run_dir.
## ---------------------------------------------------------------------------

rm(list = ls())


args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1, dir.exists(args[1]))
run_dir <- args[1]
out_dir <- if (length(args) > 1) args[2] else run_dir

# canonical marker genes, used to identify the factors and for panel (b)
list_markers <- list(
  "Interferon" = c("STAT2", "RSAD2", "IFIH1", "IRF7", "GBP1", "HERC5", "CMPK2", "DDX58",
                   "ISG15", "MX1"),
  "Erythroid" = c("HBA1", "HBA2", "HBB", "HBD", "ALAS2", "SLC25A37", "AHSP", "GYPA",
                  "EPB42", "SLC4A1"))
n_top <- 50
day_min <- 0
day_max <- 49
day_from <- 10 # interval of panel (b); 37 patients' own sampling covers it in the
day_to <- 30   # N = 164 fit, against 23 for days 10 to 35

# grey ramp for the severity classes, as in the metabolomic figures
vec_col_sev <- c("B" = "grey72", "C" = "grey54", "D" = "grey32", "E" = "black")

e <- new.env()
load(file.path(run_dir, "output_slim.RData"), envir = e)
res <- e$res_slim

days <- res$time_g * (e$days_thres - 1)
ind_win <- which(days >= day_min & days <= day_max)
d_win <- days[ind_win]
ind_int <- which(days >= day_from & days <= day_to)
d_int <- days[ind_int]

B_hat <- res$B_hat
mat_mu <- do.call(cbind, res$list_mu_hat)
colnames(mat_mu) <- rownames(B_hat)

# factor matched to each marker set
vec_q <- sapply(list_markers, function(m) {
  hits <- sapply(seq_len(ncol(B_hat)), function(q)
    sum(m %in% rownames(B_hat)[order(abs(B_hat[, q]), decreasing = TRUE)[1:200]]))
  cat("markers in the top 200, by factor (", paste(m[1:2], collapse = ", "), "...): ",
      paste(hits, collapse = " "), "\n", sep = "")
  which.max(hits)
})

# severity class of each patient, from the first sample
df_first <- e$df_comb_no_HC[order(e$df_comb_no_HC$days_from_sx_or_swab_imputed), ]
df_first <- df_first[!duplicated(df_first$subject_id), ]
subj <- names(res$list_h_hat)
vec_sev <- as.character(df_first$severity[match(subj, df_first$subject_id)])

# patients whose own sampling covers the interval of panel (b)
covers <- sapply(subj, function(s) {
  d <- e$time_obs[[s]] * (e$days_thres - 1)
  min(d) <= day_from && max(d) >= day_to
})

# the n_top genes with the largest loadings on a factor, and their mean curve
get_top <- function(q) rownames(B_hat)[order(abs(B_hat[, q]), decreasing = TRUE)[seq_len(n_top)]]
mean_curve <- function(genes) rowMeans(mat_mu[, genes, drop = FALSE])

list_curve <- lapply(vec_q, function(q) {
  mu_q <- mean_curve(get_top(q))[ind_win]
  mu_q - mean(mu_q)   # centred so the factors are comparable
})
names(list_curve) <- names(vec_q)

# class means of the marker-gene expression over the interval
list_class <- lapply("Interferon", function(nm) {
  q <- vec_q[[nm]]
  m <- intersect(list_markers[[nm]], rownames(B_hat))
  mu_m <- rowMeans(mat_mu[, m, drop = FALSE])[ind_int]
  b_m <- mean(B_hat[m, q])
  x <- sapply(subj[covers], function(s) mu_m + b_m * res$list_h_hat[[s]][[q]][ind_int])
  sapply(names(vec_col_sev), function(k) rowMeans(x[, vec_sev[covers] == k, drop = FALSE]))
})
names(list_class) <- "Interferon"


pdf(file.path(out_dir, "figure_S_rnaseq_temporal.pdf"), width = 9, height = 3.8)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.6, 1.2), mgp = c(2.5, 0.8, 0))

## (a) mean functions of the two factors
ylim_a <- range(unlist(list_curve))
plot(NULL, xlim = c(day_min, day_max), ylim = ylim_a,
     xlab = "Days from symptom onset", ylab = "Mean function (SD units, centred)",
     main = "(a) Shared resolution of the response", cex.main = 1, font.main = 1)
abline(h = 0, col = "grey70", lty = 3)
for (k in seq_along(list_curve)) {
  lines(d_win, list_curve[[k]], lwd = 2.2, lty = k)
}
legend("topright", legend = names(list_curve), lty = seq_along(list_curve),
       lwd = 2.2, bty = "n", cex = 0.78, seg.len = 3, inset = c(0, -0.02))

## (b) interferon marker-gene expression by severity class
ylim_bc <- range(unlist(list_class))
for (k in seq_along(list_class)) {
  plot(NULL, xlim = c(day_from, day_to), ylim = ylim_bc,
       xlab = "Days from symptom onset", ylab = paste(names(list_class)[k], "markers (SD units)"),
       main = paste0("(", letters[k + 1], ") ", names(list_class)[k], " markers, by severity class"),
       cex.main = 1, font.main = 1)
  abline(h = 0, col = "grey70", lty = 3)
  for (cl in names(vec_col_sev)) {
    lines(d_int, list_class[[k]][, cl], lwd = 2.4, col = vec_col_sev[cl])
  }
  legend("topright", legend = paste0(names(vec_col_sev), " (", table(factor(vec_sev[covers],
         levels = names(vec_col_sev))), ")"), col = vec_col_sev, lwd = 2.4,
         bty = "n", cex = 0.78, ncol = 2, inset = c(0, -0.02))
}

invisible(dev.off())

cat("written:", file.path(out_dir, "figure_S_rnaseq_temporal.pdf"), "\n")
cat("factors:", paste0(names(vec_q), " F", vec_q, collapse = ", "), "\n")
cat("patients covering days", day_from, "to", day_to, ":", sum(covers), "\n")
cat("\nPanel (a), change over days", day_min, "to", day_max, "(SD units):\n")
for (k in seq_along(list_curve)) {
  cat(sprintf("  %-14s %+.2f\n", names(list_curve)[k],
              list_curve[[k]][length(ind_win)] - list_curve[[k]][1]))
}
cat("\nPanel (b), class means at day", day_from, "and day", day_to, ":\n")
for (k in seq_along(list_class)) {
  cat(names(list_class)[k], "\n")
  print(round(list_class[[k]][c(1, nrow(list_class[[k]])), ], 2))
}
