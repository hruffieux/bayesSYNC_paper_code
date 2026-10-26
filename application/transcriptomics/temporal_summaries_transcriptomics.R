## ---------------------------------------------------------------------------
## Temporal summaries of the transcriptomic factors, for the text of Section 6.2
## and Supplementary Material S9 (and the last three columns of Table 3),
## computed from the output_slim.RData of an analysis_transcriptomics.R run.
##
## Four groups of quantities.
##
##   (a) For every factor, the shape of its leading eigenfunction: share of the
##       factor's variation it explains, number of sign changes between days 5
##       and 45 (the window outside which fewer than 15 patients have data),
##       the days at which it crosses zero, the number of components needed for
##       99% of the variation, and the rank correlation of the leading scores
##       with severity class. Shares are normalised within the L_max components
##       fitted, so reaching 100% at l = L_max is by construction.
##
##   (b) For the factors identified by canonical marker genes, the change in
##       the estimated mean functions averaged over the 50 genes with the
##       largest loadings on that factor. The mean functions carry the shared
##       time course; the latent trajectories are close to flat on average.
##
##   (c) For the same factors, the expression-scale contrast between day 10 and
##       day 35, mu_j(t) + b_jq h_iq(t) averaged over the marker genes, over the
##       patients whose own sampling covers days 10 to 35, so that no patient
##       trajectory is extrapolated. Marker genes load negatively on some
##       factors, so the direction of a latent trajectory says nothing about the
##       direction of expression; only this scale does.
##
##   (d) For every factor, the rank correlation of the leading scores with
##       severity class after adjusting for the number of samples, and for the
##       matched factors the share of the leading eigenfunction in the first
##       days of the window.
##
## Factors are matched to the marker sets by the number of markers among the
## 200 largest loadings, which is printed, so a weak match is visible. Factor
## numbers differ between runs.
##
## Usage, from the root of the repository:
##   Rscript application/transcriptomics/temporal_summaries_transcriptomics.R <run_dir>
## Writes rnaseq_temporal_summaries.txt into <run_dir>, and prints the same.
## ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1, dir.exists(args[1]))
run_dir <- args[1]

e <- new.env()
load(file.path(run_dir, "output_slim.RData"), envir = e)
res <- e$res_slim
Q <- e$Q
L <- e$L

out_file <- file.path(run_dir, "rnaseq_temporal_summaries.txt")
sink(out_file, split = TRUE)

days <- res$time_g * (e$days_thres - 1)
day_at <- function(d) which.min(abs(days - d))
ind_5_45 <- which(days >= 5 & days <= 45)

B_hat <- res$B_hat
genes_all <- rownames(B_hat)
mat_mu <- do.call(cbind, res$list_mu_hat)
colnames(mat_mu) <- genes_all

# severity class of each patient, from the first sample
df_first <- e$df_comb_no_HC[order(e$df_comb_no_HC$days_from_sx_or_swab_imputed), ]
df_first <- df_first[!duplicated(df_first$subject_id), ]
subj <- names(res$list_h_hat)
vec_sev <- as.character(df_first$severity[match(subj, df_first$subject_id)])
sev_num <- match(vec_sev, c("B", "C", "D", "E"))
stopifnot(!any(is.na(sev_num)))

cat("Run: ", basename(run_dir), "\n", sep = "")
cat("Patients: ", length(subj), ", classes ",
    paste(names(table(vec_sev)), table(vec_sev), sep = " ", collapse = ", "), "\n\n", sep = "")


## (a) leading eigenfunction of every factor ---------------------------------

sign_changes <- function(x) {
  s <- sign(x)
  s <- s[s != 0]
  sum(diff(s) != 0)
}
zero_days <- function(x, d) {
  i <- which(diff(sign(x)) != 0)
  round(d[i])
}

cat("(a) Leading eigenfunction of each factor\n")
cat(sprintf("%-7s %8s %10s %-18s %8s %10s %10s %10s\n", "factor", "PVE (%)", "sign ch.",
            "zero (days 5-45)", "for 99%", "|rho| sev", "p rho", "p KW"))
for (q in 1:Q) {
  phi_1 <- res$list_list_Phi_hat[[q]][, 1]
  cum_pve <- res$list_cumulated_pve[[q]]
  z_1 <- res$list_Zeta_hat[[q]][subj, 1]
  ct <- suppressWarnings(cor.test(z_1, sev_num, method = "spearman", exact = FALSE))
  p_kw <- kruskal.test(z_1, factor(vec_sev))$p.value
  cat(sprintf("%-7d %8.1f %10d %-18s %8d %10.2f %10.2g %10.2g\n", q, cum_pve[1],
              sign_changes(phi_1[ind_5_45]),
              paste(zero_days(phi_1[ind_5_45], days[ind_5_45]), collapse = ", "),
              which(cum_pve >= 99)[1], abs(ct$estimate), ct$p.value, p_kw))
}
cat("\n")


## factors matched to canonical marker sets ----------------------------------

list_markers <- list(
  "Interferon" = c("STAT2", "RSAD2", "IFIH1", "IRF7", "GBP1", "HERC5", "CMPK2", "DDX58",
                   "ISG15", "MX1"),
  "Erythroid" = c("HBA1", "HBA2", "HBB", "HBD", "ALAS2", "SLC25A37", "AHSP", "GYPA",
                  "EPB42", "SLC4A1"),
  "Proliferation" = c("PLK1", "MCM2", "RRM2", "TYMS", "NCAPD2", "TUBA1B", "TUBB", "GPI"))

top_200 <- lapply(1:Q, function(q) genes_all[order(abs(B_hat[, q]), decreasing = TRUE)[1:200]])

cat("Marker genes among the 200 largest loadings of each factor\n")
mat_hits <- sapply(list_markers, function(m) sapply(top_200, function(g) sum(m %in% g)))
rownames(mat_hits) <- paste0("F", 1:Q)
print(t(mat_hits))
vec_q <- apply(mat_hits, 2, which.max)
cat("\nMatched factors: ",
    paste0(names(vec_q), " F", vec_q, " (", apply(mat_hits, 2, max), " of ",
           lengths(list_markers), ")", collapse = ", "), "\n\n", sep = "")


## (b) mean functions over the top-50 genes of the matched factors ------------

n_top <- 50
cat("(b) Mean functions, averaged over the ", n_top, " largest loadings\n", sep = "")
for (k in seq_along(vec_q)) {
  q <- vec_q[k]
  g <- genes_all[order(abs(B_hat[, q]), decreasing = TRUE)[seq_len(n_top)]]
  mu <- rowMeans(mat_mu[, g, drop = FALSE])
  n_pos <- sum(B_hat[g, q] > 0)
  cat(sprintf("%-14s F%-3d change day 5 to 45 %+.2f, day 0 to 49 %+.2f (loadings: %d of %d positive)\n",
              names(vec_q)[k], q, mu[day_at(45)] - mu[day_at(5)],
              mu[length(mu)] - mu[1], n_pos, n_top))
}
cat("\n")


## (c) expression-scale contrast between day 10 and day 35 ---------------------

covers <- sapply(subj, function(s) {
  d <- e$time_obs[[s]] * (e$days_thres - 1)
  min(d) <= 10 && max(d) >= 35
})
cat("(c) Expression scale, marker genes, the ", sum(covers),
    " patients whose own sampling covers days 10 to 35\n", sep = "")
cat("    classes: ", paste(names(table(vec_sev[covers])), table(vec_sev[covers]),
                         sep = " ", collapse = ", "), "\n", sep = "")

for (k in seq_along(vec_q)) {
  q <- vec_q[k]
  m <- intersect(list_markers[[k]], genes_all)
  mu_m <- rowMeans(mat_mu[, m, drop = FALSE])
  b_m <- B_hat[m, q]
  # mean over markers of mu_j(t) + b_jq h_iq(t), for each patient
  x_at <- function(s, d) {
    i_d <- day_at(d)
    mu_m[i_d] + mean(b_m) * res$list_h_hat[[s]][[q]][i_d]
  }
  s_cov <- subj[covers]
  x10 <- sapply(s_cov, x_at, d = 10)
  x30 <- sapply(s_cov, x_at, d = 30)
  x35 <- sapply(s_cov, x_at, d = 35)
  ch <- x35 - x10
  ct <- suppressWarnings(cor.test(ch, sev_num[covers], method = "spearman", exact = FALSE))
  cls <- factor(vec_sev[covers], levels = c("B", "C", "D", "E"))
  cat(sprintf("%-14s F%-3d (%d markers)\n", names(vec_q)[k], q, length(m)))
  cat("   day 10 by class:      ", paste(sprintf("%+.2f", tapply(x10, cls, mean)), collapse = " / "), "\n")
  cat("   day 30 by class:      ", paste(sprintf("%+.2f", tapply(x30, cls, mean)), collapse = " / "), "\n")
  cat("   change 10 to 35:      ", paste(sprintf("%+.2f", tapply(ch, cls, mean)), collapse = " / "), "\n")
  cat(sprintf("   rank correlation of the change with class %+.2f (p = %.2g)\n", ct$estimate, ct$p.value))
}


## (d) scores against severity, allowing for the number of samples -----------
##
## The number of samples differs between severity classes, and the scores of
## patients with few samples are shrunk towards zero, so an association between
## scores and severity could arise from the sampling alone. The partial rank
## correlation removes the number of samples from the ranks of both the score
## and the class. Also printed, for the matched factors, the share of the
## squared leading eigenfunction before days 5 and 10, with the number of
## patients sampled in each interval, since weight placed where few patients
## have data rests on few observations.

partial_rank_cor <- function(x, y, g) {
  cor(resid(lm(rank(x) ~ factor(g))), resid(lm(rank(y) ~ factor(g))))
}
n_obs <- sapply(e$time_obs[subj], length)
cat("\n(d) Number of samples by severity class (mean):",
    paste(names(table(vec_sev)), sprintf("%.2f", tapply(n_obs, vec_sev, mean)), collapse = ", "), "\n")
cat(sprintf("%-7s %8s %10s %14s %14s\n", "factor", "rho sev", "adjusted", "n_i = 1", "n_i >= 2"))
for (q in 1:Q) {
  z_1 <- res$list_Zeta_hat[[q]][subj, 1]
  rho_sub <- function(k) cor(z_1[k], sev_num[k], method = "spearman")
  cat(sprintf("%-7d %8.2f %10.2f %8.2f (%3d) %8.2f (%3d)\n", q, rho_sub(TRUE),
              partial_rank_cor(z_1, sev_num, n_obs), rho_sub(n_obs == 1), sum(n_obs == 1),
              rho_sub(n_obs >= 2), sum(n_obs >= 2)))
}

days_obs <- unlist(lapply(subj, function(s) e$time_obs[[s]] * (e$days_thres - 1)))
subj_obs <- rep(subj, sapply(e$time_obs[subj], length))
cat("\nPatients with a sample before day 5:", length(unique(subj_obs[days_obs < 5])),
    "; between days 5 and 10:", length(unique(subj_obs[days_obs >= 5 & days_obs < 10])), "\n")
for (k in seq_along(vec_q)) {
  q <- vec_q[k]
  w <- res$list_list_Phi_hat[[q]][, 1]^2
  w <- w / sum(w)
  cat(sprintf("%-14s F%-3d share of the leading eigenfunction before day 5 %.2f, before day 10 %.2f\n",
              names(vec_q)[k], q, sum(w[days < 5]), sum(w[days < 10])))
}

sink()
cat("written:", out_file, "\n")
