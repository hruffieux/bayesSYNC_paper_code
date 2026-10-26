## ---------------------------------------------------------------------------
## Diagnostics of the biomarker panel fit, computed from
## the output.RData written by analysis_biomarker_panel.R (save.image), without refitting.
##
##   - number of subjects and of observations per subject, by severity class
##     (Section 2);
##   - whether the subjects with extreme scores are those with few
##     observations. The "extreme subjects" figures show, for each factor and
##     component, the subjects with the smallest and the largest score, and
##     whether the association of the scores with severity survives adjustment
##     for the number of observations (Section 6.1);
##   - residual correlation against the time between two observations of the
##     same subject, the statistic reported in the misspecification study
##     (Section 5.5), so that the application can be read on the same scale
##     (Supplementary Material S6);
##   - sensitivity to the PPI threshold of the variables assigned to each
##     retained factor (Section 6.1).
##
## Everything is printed and also written to run_dir, next to output.RData.
##
## Usage, from the root of the repository:
##   Rscript application/biomarker_panel/diagnostics_biomarker_panel.R <run_dir>
## ---------------------------------------------------------------------------

rm(list = ls())

source("R/paths.R")

source("R/utils_misspec.R") # residual_correlogram()

# folder holding the output.RData of the run to examine
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1, dir.exists(args[1]))
run_dir <- args[1]

vec_thres_ppi <- c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99) # 0.5 is the threshold used in the paper
thres_ref <- 0.5


## ---------------------------------------------------------------------------
## Loading
##
## save.image stores the whole workspace of analysis_biomarker_panel.R, so the objects are
## taken from a separate environment rather than loaded over this one.
## ---------------------------------------------------------------------------

e <- new.env()
load(file.path(run_dir, "output.RData"), envir = e)
stopifnot(all(c("res", "time_obs", "df_info", "days_thres") %in% ls(e)))

res <- e$res
time_obs <- e$time_obs
days_thres <- e$days_thres
# prepare_data_irregular() rescales the times by days_thres - 1, so that the last
# admissible day (days_thres - 1) maps to 1
days_scale <- days_thres - 1
subject_ids <- names(time_obs)
n_obs <- sapply(time_obs, length)
severity <- as.character(e$df_info$severity[match(subject_ids, e$df_info$subject_id)])
stopifnot(!any(is.na(severity)))

retained <- which(res$factor_ppi > 0.5)
L_show <- 2 # leading components, as in the figures

sink(file.path(run_dir, "diagnostics.txt"), split = TRUE)


## ---------------------------------------------------------------------------
## Analysis set
## ---------------------------------------------------------------------------

cat("=== ANALYSIS SET ===\n")
cat("Subjects: ", length(subject_ids), "   total observations: ", sum(n_obs),
    "   variables: ", length(res$Y[[1]]), "\n", sep = "")
cat("Observations per subject: min ", min(n_obs), ", median ", median(n_obs),
    ", max ", max(n_obs), "\n", sep = "")
cat("\nNumber of observations per subject:\n")
print(table(n_obs))
cat("\nSeverity class by number of observations per subject:\n")
print(table(severity = severity, n_obs = n_obs))
cat("\nObservation times (days from symptom onset), quantiles:\n")
print(round(quantile(unlist(time_obs) * days_scale, c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)), 1))
cat("Retained factors: ", paste(retained, collapse = ", "), "\n\n", sep = "")


## ---------------------------------------------------------------------------
## Number of observations of the subjects with extreme scores
##
## A positive rank correlation between n_i and |score| means that subjects with
## more observations get larger scores in absolute value, as expected from the
## shrinkage of weakly observed subjects towards zero; a negative one means the
## extreme scores come disproportionately from subjects with few observations.
## ---------------------------------------------------------------------------

cat("=== EXTREME SCORES AND NUMBER OF OBSERVATIONS ===\n")
cat("(median number of observations over all subjects: ", median(n_obs), ")\n", sep = "")

df_extreme <- NULL
for (q in retained) {
  for (l in seq_len(L_show)) {
    sc <- res$list_Zeta_hat[[q]][, l]
    i_min <- which.min(sc)
    i_max <- which.max(sc)
    sp <- suppressWarnings(cor.test(n_obs, abs(sc), method = "spearman"))
    p_sev <- kruskal.test(sc ~ factor(severity))$p.value
    cat(sprintf("Factor %d, FPC %d: smallest score %s (n_i = %d, class %s), largest %s (n_i = %d, class %s); Spearman rho(n_i, |score|) = %.2f (p = %.2g); severity p = %.2g\n",
                q, l, subject_ids[i_min], n_obs[i_min], severity[i_min],
                subject_ids[i_max], n_obs[i_max], severity[i_max],
                sp$estimate, sp$p.value, p_sev))
    df_extreme <- rbind(df_extreme, data.frame(factor = q, fpc = l,
                                       subj_min = subject_ids[i_min], n_min = n_obs[i_min],
                                       subj_max = subject_ids[i_max], n_max = n_obs[i_max],
                                       rho_n_abs_score = unname(sp$estimate),
                                       p_rho = sp$p.value, p_severity = p_sev))
  }
}
write.csv(df_extreme, file.path(run_dir, "diagnostics_extreme_scores.csv"), row.names = FALSE)

cat("\nAbsolute score by number of observations (median over subjects):\n")
for (q in retained) {
  for (l in seq_len(L_show)) {
    cat(sprintf("Factor %d, FPC %d: ", q, l))
    print(round(tapply(abs(res$list_Zeta_hat[[q]][, l]), n_obs, median), 3))
  }
}
cat("\n")

# The number of observations differs between severity classes, and the scores
# of weakly observed subjects are shrunk towards zero, so an association between
# scores and severity could arise from the sampling alone. The partial rank
# correlation removes the number of observations from the ranks of both the
# score and the class (Section 6.1).
partial_rank_cor <- function(x, y, g) {
  cor(resid(lm(rank(x) ~ factor(g))), resid(lm(rank(y) ~ factor(g))))
}
sev_num <- match(severity, c("B", "C", "D", "E"))
cat("Mean number of observations by severity class:\n")
print(round(tapply(n_obs, severity, mean), 2))
cat("\nRank correlation of the scores with severity class:\n")
cat(sprintf("%-16s %8s %10s %14s %14s\n", "", "all", "adjusted", "n_i = 1", "n_i >= 2"))
for (q in retained) {
  for (l in seq_len(L_show)) {
    sc <- res$list_Zeta_hat[[q]][, l]
    rho_sub <- function(k) cor(sc[k], sev_num[k], method = "spearman")
    cat(sprintf("Factor %d, FPC %d %8.2f %10.2f %8.2f (%3d) %8.2f (%3d)\n", q, l,
                rho_sub(TRUE), partial_rank_cor(sc, sev_num, n_obs),
                rho_sub(n_obs == 1), sum(n_obs == 1), rho_sub(n_obs >= 2), sum(n_obs >= 2)))
  }
}
cat("\n")


## ---------------------------------------------------------------------------
## Residual correlogram (Supplementary Material S6)
##
## Bins of 0.1 on the rescaled time axis, about five days (0.1 x 49), since subjects have
## at most a handful of observations; the number of contributing pairs is shown
## so that sparsely populated bins are visible.
## ---------------------------------------------------------------------------

cat("=== RESIDUAL CORRELATION AGAINST THE TIME BETWEEN OBSERVATIONS ===\n")
df_corr <- residual_correlogram(res, time_obs, breaks = seq(0, 1, by = 0.1))
df_corr$days_low <- df_corr$dist_low * days_scale
df_corr$days_upp <- df_corr$dist_upp * days_scale
cat("Subjects with at least two observations: ", sum(n_obs >= 2), "\n", sep = "")
print(df_corr[, c("days_low", "days_upp", "n_pair", "corr")], row.names = FALSE, digits = 3)
write.csv(df_corr, file.path(run_dir, "diagnostics_residual_correlogram.csv"), row.names = FALSE)

pdf(file.path(run_dir, "residual_correlogram.pdf"), width = 5, height = 4.5)
ok <- df_corr$n_pair > 0
plot(df_corr$days_low[ok] + diff(df_corr$days_low)[1] / 2, df_corr$corr[ok], type = "o", pch = 19,
     xlab = "Days between the two observations", ylab = "Residual correlation",
     ylim = range(c(0, df_corr$corr[ok]), na.rm = TRUE))
abline(h = 0, lty = 2, col = "grey50")
dev.off()
cat("\n")


## ---------------------------------------------------------------------------
## Sensitivity of the selected variables to the PPI threshold
##
## For each retained factor, the number of variables with PPI above each
## threshold and the Jaccard overlap with the set selected at 0.5. The fitted
## trajectories and the factor structure do not depend on the threshold, which
## only decides which variables are reported as belonging to a factor.
## ---------------------------------------------------------------------------

cat("=== SENSITIVITY TO THE PPI THRESHOLD ===\n")
jaccard <- function(a, b) if (length(union(a, b)) == 0) NA else length(intersect(a, b)) / length(union(a, b))
df_thres <- NULL
for (q in retained) {
  ppi_q <- res$ppi[, q]
  sel_ref <- names(ppi_q)[ppi_q > thres_ref]
  for (thres in vec_thres_ppi) {
    sel <- names(ppi_q)[ppi_q > thres]
    df_thres <- rbind(df_thres, data.frame(factor = q, threshold = thres, n_selected = length(sel),
                                           jaccard_with_0.5 = jaccard(sel, sel_ref)))
  }
  cat(sprintf("Factor %d, PPI distribution: ", q))
  print(table(cut(ppi_q, c(-Inf, 0.1, 0.5, 0.9, 0.99, Inf))))
}
print(df_thres, row.names = FALSE, digits = 3)
write.csv(df_thres, file.path(run_dir, "diagnostics_ppi_threshold.csv"), row.names = FALSE)

sink()
