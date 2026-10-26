## ---------------------------------------------------------------------------
## Analysis of the simulated biomarker panel, compared with the values that
## generated it (Supplementary Material S10).
##
## Reads the output.RData of analysis_biomarker_panel.R run on the simulated
## data, and data/simulated_biomarker_panel_truth.rds written by
## simulate_cohort_data.R. Simulated and retained factors are matched one to
## one, by decreasing absolute correlation of their loadings; a simulated factor
## left without a retained factor is reported as not recovered. The signs of the estimated loadings, eigenfunctions and scores are set
## to agree with the truth before comparison (the model identifies them up to
## sign). Reported:
##   - the posterior inclusion probability of each candidate factor;
##   - for each matched pair, the correlation between true and estimated
##     loadings and scores, and the cosine between true and estimated
##     eigenfunctions, for the leading components;
##   - the rank correlation of the true and of the estimated leading scores
##     with the severity class.
##
## Writes, next to output.RData, compare_with_truth.txt, compare_with_truth.csv
## (one row per simulated factor) and
## figure_S_simulated_panel.pdf (loadings, leading eigenfunction and leading
## scores, estimated against true).
##
## Usage, from the root of the repository:
##   Rscript application/biomarker_panel/compare_with_truth_simulated.R <run_dir>
## ---------------------------------------------------------------------------

rm(list = ls())

source("R/paths.R")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1, dir.exists(args[1]))
run_dir <- args[1]

e <- new.env()
load(file.path(run_dir, "output.RData"), envir = e)
stopifnot(isTRUE(e$bool_simulated_data))
res <- e$res
truth <- readRDS(file.path(data_dir, "simulated_biomarker_panel_truth.rds"))
# the analysis of Supplementary Material S10 fits the simulated data with the
# seed that generated them, the default of analysis_biomarker_panel.R
if (!identical(as.integer(truth$seed), as.integer(e$seed))) {
  stop("the simulated data in ", data_dir, " were generated with seed ", truth$seed,
       " but the fit in ", run_dir, " used seed ", e$seed,
       "; run analysis_biomarker_panel.R without a seed argument")
}
cat("Simulated data and fit: seed", truth$seed, "\n")

subj <- names(e$time_obs)
stopifnot(all(subj %in% names(truth$list_zeta)))
sev_num <- match(truth$severity[match(subj, names(truth$list_zeta))], c("B", "C", "D", "E"))
stopifnot(identical(rownames(truth$B), rownames(res$B_hat)) || is.null(rownames(truth$B)))

days <- res$time_g * (e$days_thres - 1)
# truth and fit on the same grid, otherwise interpolate the truth
Phi_true_on <- function(k) {
  if (length(truth$time_g) == length(res$time_g) && max(abs(truth$time_g - res$time_g)) < 1e-8) {
    truth$list_Phi[[k]]
  } else {
    apply(truth$list_Phi[[k]], 2, function(x) approx(truth$time_g, x, xout = res$time_g)$y)
  }
}

# agreement of two eigenfunctions, the cosine of the angle between them in
# L2, approximated on the grid. Unlike a correlation it does not centre the
# curves, so it remains meaningful for a nearly constant eigenfunction.
cos_l2 <- function(x, y) sum(x * y) / sqrt(sum(x^2) * sum(y^2))

K_true <- ncol(truth$B)
retained <- which(res$factor_ppi > 0.5)
L_cmp <- 2 # leading components compared

out_file <- file.path(run_dir, "compare_with_truth.txt")
sink(out_file, split = TRUE)

cat("Patients:", length(subj), "  simulated factors:", K_true, "  retained factors:", length(retained), "\n")
cat("Factor PPI:", paste(sprintf("%.2f", res$factor_ppi), collapse = " "), "\n\n")

# one-to-one matching, by decreasing absolute loading correlation
mat_r <- matrix(NA, K_true, length(retained))
for (k in seq_len(K_true)) for (m in seq_along(retained)) mat_r[k, m] <- cor(truth$B[, k], res$B_hat[, retained[m]])
vec_q <- rep(NA, K_true)
a <- abs(mat_r)
while (any(!is.na(a)) && any(is.na(vec_q))) {
  ij <- which(a == max(a, na.rm = TRUE), arr.ind = TRUE)[1, ]
  vec_q[ij[1]] <- retained[ij[2]]
  a[ij[1], ] <- NA
  a[, ij[2]] <- NA
}

list_match <- list()
df_summary <- data.frame(factor = seq_len(K_true), n_retained = length(retained), recovered = !is.na(vec_q),
                         cor_loadings = NA, cos_phi_1 = NA, cor_scores_1 = NA,
                         rho_sev_true_1 = NA, rho_sev_est_1 = NA)
for (k in seq_len(K_true)) {
  if (is.na(vec_q[k])) {
    cat(sprintf("Simulated factor %d: not recovered%s\n", k,
                if (length(retained) > 0) sprintf(" (largest absolute loading correlation with a retained factor %.3f)",
                                                  max(abs(mat_r[k, ]))) else ""))
    next
  }
  q <- vec_q[k]
  r_kq <- mat_r[k, match(q, retained)]
  s_b <- sign(r_kq)
  Phi_true <- Phi_true_on(k)
  Z_true <- t(sapply(subj, function(s) truth$list_zeta[[s]][[k]]))
  Z_hat <- res$list_Zeta_hat[[q]][subj, , drop = FALSE]
  cat(sprintf("Simulated factor %d -> estimated factor %d: loading correlation %.3f\n", k, q, abs(r_kq)))
  for (l in seq_len(L_cmp)) {
    # sign of the component, from the eigenfunctions; scores follow, with the
    # sign of the loadings, since b h = (s b)(s h)
    c_phi <- cos_l2(Phi_true[, l], res$list_list_Phi_hat[[q]][, l])
    s_phi <- sign(c_phi)
    z_hat <- s_b * s_phi * Z_hat[, l]
    cat(sprintf("  component %d: eigenfunction cosine %.3f, score correlation %.3f,",
                l, abs(c_phi), cor(Z_true[, l], z_hat)))
    cat(sprintf(" severity rank correlation true %.2f, estimated %.2f\n",
                cor(Z_true[, l], sev_num, method = "spearman"), cor(z_hat, sev_num, method = "spearman")))
  }
  list_match[[length(list_match) + 1]] <- list(k = k, q = q, s_b = s_b, Phi_true = Phi_true, Z_true = Z_true)
  c_1 <- cos_l2(Phi_true[, 1], res$list_list_Phi_hat[[q]][, 1])
  z_1 <- s_b * sign(c_1) * Z_hat[, 1]
  df_summary[k, c("cor_loadings", "cos_phi_1", "cor_scores_1", "rho_sev_true_1", "rho_sev_est_1")] <-
    c(abs(r_kq), abs(c_1), cor(Z_true[, 1], z_1),
      cor(Z_true[, 1], sev_num, method = "spearman"), cor(z_1, sev_num, method = "spearman"))
}
write.csv(df_summary, file.path(run_dir, "compare_with_truth.csv"), row.names = FALSE)
sink()

# the simulated factors follow the columns of the calibration, factors 4
# (amino acids) and 5 (inflammatory) of the real fit; the figure shows the
# inflammatory factor first, as Figure 6 does
vec_label <- c("Amino-acid factor", "Inflammatory factor")
list_match <- list_match[order(match(sapply(list_match, `[[`, "k"), c(2, 1)))]

if (length(list_match) > 0) {
  vec_col_sev <- c("grey88", "grey62", "grey35", "grey5")
  pdf(file.path(run_dir, "figure_S_simulated_panel.pdf"), width = 9, height = 2.9 * length(list_match))
  par(mfrow = c(length(list_match), 3), mar = c(4, 4, 2.5, 1))
  for (mm in list_match) {
    k <- mm$k
    q <- mm$q
    phi_hat <- res$list_list_Phi_hat[[q]][, 1]
    s_phi <- sign(cos_l2(mm$Phi_true[, 1], phi_hat))
    plot(truth$B[, k], mm$s_b * res$B_hat[, q], pch = 20, xlab = "True loading",
         ylab = "Estimated loading", main = paste0(vec_label[k], ": loadings"), cex.main = 1)
    abline(0, 1, lty = 3)
    plot(days, mm$Phi_true[, 1], type = "l", lwd = 1.5, col = "grey60",
         ylim = range(mm$Phi_true[, 1], phi_hat), xlab = "Days from symptom onset",
         ylab = "", main = "Leading eigenfunction", cex.main = 1)
    lines(days, s_phi * phi_hat, lwd = 1.5)
    abline(h = 0, lty = 3)
    if (identical(mm, list_match[[1]])) legend("topright", bty = "n", lwd = 1.5, col = c("grey60", "black"),
                       legend = c("true", "estimated"), cex = 0.8)
    plot(mm$Z_true[, 1], mm$s_b * s_phi * res$list_Zeta_hat[[q]][subj, 1], pch = 21,
         bg = vec_col_sev[sev_num], col = "grey30", xlab = "True score", ylab = "Estimated score",
         main = "Leading scores", cex.main = 1)
    abline(0, 1, lty = 3)
    if (identical(mm, list_match[[1]])) legend("topleft", bty = "n", pch = 21, pt.bg = vec_col_sev, col = "grey30",
                       legend = c("B", "C", "D", "E"), cex = 0.8)
  }
  dev.off()
}
cat("written:", out_file, "\n")
