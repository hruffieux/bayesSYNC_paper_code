## ---------------------------------------------------------------------------
## Supplementary figure of Supplementary Material S9
## (figure_S_rnaseq_factors.pdf), the counterpart for the transcriptomic
## analysis of Figure 6, drawn from the output_slim.RData of an
## analysis_transcriptomics.R run.
##
## One row per interpretable factor (neutrophil/lymphocyte, interferon,
## erythroid):
##   (a) the 15 largest loadings, by absolute value;
##   (b), (c) the first two eigenfunctions, with their share of the variation of
##       the factor; days before day 5 and after day 45, where few patients
##       have a sample (the numbers are printed), are shaded;
##   (d) the scores on these two components, by severity class.
##
## Signs. A factor is identified up to the sign of its loadings and of its
## trajectories, and each component up to the joint sign of its eigenfunction
## and scores. Each factor is oriented so that most of its 50 largest loadings
## are positive, so that a higher trajectory means higher expression of those
## genes; each eigenfunction is then oriented to have a positive mean over days
## 5 to 45, its scores being flipped with it. B h_i is unchanged.
##
## The factor numbers are those of the N = 164 fit. On a fit to the simulated
## data, or another run, the three factors whose leading scores correlate most
## with the severity class are shown, without labels.
##
## Usage, from the root of the repository:
##   Rscript application/transcriptomics/plot_figure_S_transcriptomics_factors.R <run_dir> [out_dir]
## The figure is written in out_dir, by default run_dir.
## ---------------------------------------------------------------------------

rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 1:2, dir.exists(args[1]))
run_dir <- args[1]
out_dir <- if (length(args) == 2) args[2] else run_dir

e <- new.env()
load(file.path(run_dir, "output_slim.RData"), envir = e)
res <- e$res_slim
Q <- e$Q

vec_factors <- c("Neutrophil/lymphocyte" = 1, "Interferon" = 6, "Erythroid" = 10)

n_top_loadings <- 15
n_top_orient <- 50

days <- res$time_g * (e$days_thres - 1)
ind_5_45 <- which(days >= 5 & days <= 45)

subj <- names(res$list_h_hat)
df_first <- e$df_comb_no_HC[order(e$df_comb_no_HC$days_from_sx_or_swab_imputed), ]
df_first <- df_first[!duplicated(df_first$subject_id), ]
vec_sev <- as.character(df_first$severity[match(subj, df_first$subject_id)])
sev_num <- match(vec_sev, c("B", "C", "D", "E"))
stopifnot(!any(is.na(sev_num)))

if (isTRUE(e$bool_simulated_data) || e$N != 164) {
  rho <- sapply(1:Q, function(q) abs(cor(res$list_Zeta_hat[[q]][subj, 1], sev_num, method = "spearman")))
  vec_factors <- order(rho, decreasing = TRUE)[1:3]
  names(vec_factors) <- paste("Factor", vec_factors)
}
stopifnot(all(vec_factors <= Q))

# grey levels of Figure 6, spread slightly so that classes B and C are told apart
vec_col_sev <- c(B = "grey88", C = "grey62", D = "grey35", E = "grey5")
vec_lab_sev <- c(B = "B: Screen. symptomatic", C = "C: Hosp. no O2",
                 D = "D: Hosp. suppl. O2", E = "E: Hosp. assisted ventil.")

# patients with a sample in the shaded ends of the window
days_obs <- unlist(lapply(subj, function(s) e$time_obs[[s]] * (e$days_thres - 1)))
subj_obs <- rep(subj, sapply(e$time_obs[subj], length))
n_before_5 <- length(unique(subj_obs[days_obs < 5]))
n_after_45 <- length(unique(subj_obs[days_obs > 45]))

out_file <- file.path(out_dir, "figure_S_rnaseq_factors.pdf")
pdf(out_file, width = 11, height = 2.9 * length(vec_factors))
par(mfrow = c(length(vec_factors), 4), oma = c(0, 2, 0, 0))

for (k in seq_along(vec_factors)) {

  q <- vec_factors[k]
  b <- res$B_hat[, q]
  Phi <- res$list_list_Phi_hat[[q]][, 1:2]
  Zeta <- res$list_Zeta_hat[[q]][subj, 1:2]
  pve <- diff(c(0, res$list_cumulated_pve[[q]]))[1:2]

  # orientation, see the header
  s_q <- sign(sum(sign(b[order(abs(b), decreasing = TRUE)[1:n_top_orient]])))
  if (s_q == 0) s_q <- 1
  b <- s_q * b
  Phi <- s_q * Phi
  for (l in 1:2) {
    t_l <- sign(mean(Phi[ind_5_45, l]))
    if (t_l == 0) t_l <- 1
    Phi[, l] <- t_l * Phi[, l]
    Zeta[, l] <- t_l * Zeta[, l]
  }

  # (a) largest loadings
  par(mar = c(4, 7, 2.5, 1))
  b_top <- rev(b[order(abs(b), decreasing = TRUE)[1:n_top_loadings]])
  barplot(b_top, horiz = TRUE, las = 1, cex.names = 0.75, border = NA,
          col = ifelse(b_top > 0, "#ff7f6e", "skyblue1"), xlab = "Posterior mean of B")
  title(main = "Largest loadings", cex.main = 1)
  mtext(names(vec_factors)[k], side = 2, line = 6, cex = 0.8, font = 2)

  # (b), (c) eigenfunctions
  par(mar = c(4, 4, 2.5, 1))
  ylim <- range(Phi)
  for (l in 1:2) {
    plot(days, Phi[, l], type = "n", ylim = ylim, xlab = "Days from symptom onset",
         ylab = "", main = bquote(italic(l) == .(l)), cex.main = 1)
    rect(min(days), ylim[1] - 1, 5, ylim[2] + 1, col = "grey92", border = NA)
    rect(45, ylim[1] - 1, max(days), ylim[2] + 1, col = "grey92", border = NA)
    abline(h = 0, lty = 3)
    lines(days, Phi[, l], lwd = 1.5)
    legend("bottomright", bty = "n", legend = sprintf("PVE = %.1f%%", pve[l]), cex = 0.85)
    box()
  }

  # (d) scores by severity class
  plot(Zeta[, 1], Zeta[, 2], pch = 21, bg = vec_col_sev[vec_sev], col = "grey30", cex = 0.9,
       xlab = expression(italic(l) == 1), ylab = expression(italic(l) == 2),
       main = "Scores", cex.main = 1)
  if (k == 1) legend("topleft", bty = "n", pch = 21, pt.bg = vec_col_sev, col = "grey30",
                     legend = vec_lab_sev, cex = 0.7)
}

dev.off()
cat("Patients with a sample before day 5:", n_before_5, "; after day 45:", n_after_45, "\n")
cat("written:", out_file, "\n")
