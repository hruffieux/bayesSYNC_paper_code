## Redraws the figures of section 12 of analysis_transcriptomics.R from the
## output_slim.RData of a completed run, without the input data and without
## refitting. Written because the N = 164 run (deaths included) wrote every
## figure as a blank page although the saved fit is complete.
##
## Usage, from the root of the repository:
##   Rscript application/transcriptomics/replot_transcriptomics_from_slim.R <run_dir>
## The figures are written into <run_dir>, replacing those of the same name.
##
## Two differences from the figures of the main script, both due to what the
## slim file holds:
##   - the healthy-control interquartile band is absent from the trajectory
##     figures, since df_comb_no_HC has no control samples;
##   - the trajectory figures can only show the variables in vars_plotted, which
##     are the ones the main script plots anyway.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1, dir.exists(args[1]))
res_dir <- normalizePath(args[1])

output_file_slim <- file.path(res_dir, "output_slim.RData")
stopifnot(file.exists(output_file_slim))

load(output_file_slim)
res <- res_slim
source("R/plot_utils.R")

# The Y saved in the slim file is on the original log-CPM scale, whereas the
# fitted curves are on the standardised scale when bool_scale is TRUE (res$Y,
# which the main script plots, is dropped from res_slim). The transformation is
# the one the main script applies to df_comb after the fit.
if (bool_scale) {
  Y <- lapply(Y, function(y_i) {
    for (vv in names(y_i)) {
      y_i[[vv]] <- (y_i[[vv]] - res$mean_mean_across_subjects[vv]) / res$sd_mean_across_subjects[vv]
    }
    y_i
  })
}
res$Y <- Y

# plot_trajectory() indexes the fitted-curve lists by the position of the
# variable among all p variables, whereas the slim lists only hold vars_plotted.
# Padding them back to length p, with NULL elsewhere, keeps the indexing of the
# main script and therefore the same random variables and file name.
pad_to_all_var <- function(list_slim) {
  lapply(list_slim, function(x_i) {
    x_full <- vector("list", length(all_var))
    names(x_full) <- all_var
    x_full[names(x_i)] <- x_i
    x_full
  })
}
for (nm in c("list_Y_hat", "list_Y_low", "list_Y_upp")) {
  res[[nm]] <- pad_to_all_var(res[[nm]])
}

# Sanity check of the scale: the observations of a plotted variable should lie
# within the range of its fitted curves, not several units away from it.
vv_check <- vars_plotted[1]
obs_check <- unlist(lapply(res$Y, function(y_i) y_i[[vv_check]]))
fit_check <- unlist(lapply(res$list_Y_hat, function(y_i) y_i[[vv_check]]))
cat("Scale check on ", vv_check, ": observations ", paste(format(range(obs_check), digits = 3), collapse = " to "),
    ", fitted curves ", paste(format(range(fit_check), digits = 3), collapse = " to "), "\n", sep = "")

# no control samples in the slim file, see above
df_comb <- df_comb_no_HC

thres <- 0.5


## ---------------------------------------------------------------------------
## 12. Output (as in the main script, same random draws)
## ---------------------------------------------------------------------------

set.seed(seed)

vec_col_grey <- rep("white", N)
names(vec_col_grey) <- as.character(df_comb$severity[match(names(Y), df_comb$subject_id)])

vec_col_grey[names(vec_col_grey) %in% "A"] <- "grey90"
vec_col_grey[names(vec_col_grey) %in% "B"] <- "grey80"
vec_col_grey[names(vec_col_grey) %in% "C"] <- "grey65"
vec_col_grey[names(vec_col_grey) %in% "D"] <- "grey45"
vec_col_grey[names(vec_col_grey) %in% "E"] <- "grey5"

disp <- Q %/% 5 + 1 # integer division
pdf(paste0(res_dir, "/scores_1_vs_2.pdf"), height = 3*disp, width = 2.5*ceiling(Q/disp))
par(mfrow = c(disp, ceiling(Q/disp)), pty = "s")
null <- sapply(1:Q, function(q) plot(res$list_Zeta_hat[[q]][,1], res$list_Zeta_hat[[q]][,2],
                                     main = paste0("Scores factor q = ", q),
                                     xlab = "Scores FPC 1", ylab = "Scores FPC 2", pch = 20,
                                     col = vec_col_grey))
dev.off()

pdf(paste0(res_dir, "/mean_functions.pdf"), width = 8.5, height = 6)
par(mfrow = c(2, 3))
null <- sapply(sample(1:p, 6), function(j) plot_mean_functions(j, res$time_g, mu_g = NULL,
                                                               res$list_mu_hat))
dev.off()

pdf(paste0(res_dir, "/eigenfunctions.pdf"), width = 3*L, height = 3*Q)
cum_pve_thres <- 0.99
par(mfrow = c(Q, L))
null <- sapply(1:Q, function(q)
  plot_factor_eigenfunctions(q, res$time_g, Psi_g = NULL,
                             res$list_list_Phi_hat,
                             res$list_cumulated_pve,
                             cum_pve_thres = cum_pve_thres))
dev.off()

n_random <- 6
vec_random_i <- sample(1:N, n_random)
vec_random_j <- sample(1:p, n_random)
stopifnot(all(all_var[vec_random_j] %in% vars_plotted))
pdf(paste0(res_dir, "/traj_i_", paste0(vec_random_i, collapse = "-"), "_j_",
           paste0(vec_random_j, collapse = "-"), ".pdf"),
    width = 10, height = 6)
par(mfrow = c(2, 3))
for (id in 1:n_random) {
  plot_trajectory(vec_random_i[id], vec_random_j[id],
                  time_obs, res$time_g, res$Y,
                  res$list_Y_hat, res$list_Y_low, res$list_Y_upp, days_thres = days_thres,
                  df_comb = df_comb)
}
dev.off()

pdf(paste0(res_dir, "/h_i_extreme_subjects.pdf"), width = 9, height = 2.6*Q)
par(mfrow = c(Q, 4))
for (q in 1:Q) {
  for (l in 1:2) {
    plot_factor_trajectory(res$list_Zeta_hat, q, l,
                           res$time_g,
                           res$list_h_hat, res$list_h_low, res$list_h_upp,
                           days_thres = days_thres, df_comb = df_comb)
  }
}
dev.off()

n_var_display <- 4  # genes shown in the trajectory figure
for (q in 1:Q) {
  if (any(res$ppi[,q] > thres)) {
    top_B_hat_q <- res$B_hat[,q][res$ppi[,q] > thres]
    top_genes_q <- names(top_B_hat_q)[order(abs(top_B_hat_q), decreasing = TRUE)]
    n_show <- min(n_var_display, length(top_genes_q))
    pdf(paste0(res_dir, "/Y_i_extreme_subjects_variables_contributing_to_factor_q_", q,
               "_top_loading_magnitude.pdf"), width = 9, height = 2.2 * n_show)
    plot_factor_variables_trajectory(res$list_Zeta_hat, res$ppi, q, time_obs,
                                     res$time_g, res$Y, res$list_Y_hat,
                                     res$list_Y_low, res$list_Y_upp, thres_ppi = thres,
                                     df_comb = df_comb,
                                     vec_variables = head(top_genes_q, n_show),
                                     days_thres = days_thres, show_ppi = FALSE)
    dev.off()
  }
}


## ---------------------------------------------------------------------------
## 12a. Loadings and factor activation
## ---------------------------------------------------------------------------

n_top_loadings <- 20 # genes shown per factor

top_loadings_q <- function(res, q, thres, n_top) {
  sel <- res$ppi[, q] > thres
  if (!any(sel)) return(NULL)
  b <- res$B_hat[sel, q]
  b[order(abs(b), decreasing = TRUE)][1:min(n_top, sum(sel))]
}

pdf(paste0(res_dir, "/top_loadings_by_factor.pdf"),
    width = 4.2 * Q, height = 0.22 * n_top_loadings + 1.8)
par(mfrow = c(1, Q), mar = c(4.5, 8, 3, 1))
for (q in 1:Q) {
  b <- top_loadings_q(res, q, thres, n_top_loadings)
  if (is.null(b)) {
    plot.new()
    title(main = paste0("Factor ", q, "\n(no gene with PPI > ", thres, ")"))
    next
  }
  b <- rev(b) # largest at the top of the panel
  barplot(b, horiz = TRUE, las = 1, cex.names = 0.6,
          col = ifelse(b > 0, "grey30", "grey70"),
          border = NA, xlab = "Estimated loading",
          main = paste0("Factor ", q, " (PPI = ", format(res$factor_ppi[q], digits = 2), ")"))
  abline(v = 0, col = "black")
}
dev.off()

genes_union <- unique(unlist(lapply(1:Q, function(q) names(top_loadings_q(res, q, thres, n_top_loadings)))))

if (length(genes_union) > 0) {
  B_sub <- res$B_hat[genes_union, , drop = FALSE]
  B_sub[res$ppi[genes_union, , drop = FALSE] <= thres] <- NA
  ord <- order(apply(abs(B_sub), 1, which.max),
               apply(B_sub, 1, function(x) x[which.max(abs(x))]))
  B_sub <- B_sub[ord, , drop = FALSE]
  lim <- max(abs(B_sub), na.rm = TRUE)
  breaks <- seq(-lim, lim, length.out = 101)
  cols <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  pdf(paste0(res_dir, "/loadings_heatmap_top_genes.pdf"),
      width = 1.1 * Q + 3, height = 0.16 * nrow(B_sub) + 1.5)
  par(mfrow = c(1, 1), mar = c(4, 8, 3, 4))
  image(x = 1:Q, y = 1:nrow(B_sub), z = t(B_sub), col = cols, breaks = breaks,
        axes = FALSE, xlab = "Factor", ylab = "",
        main = paste0("Loadings of the top ", n_top_loadings, " genes per factor"))
  axis(1, at = 1:Q, labels = 1:Q)
  axis(2, at = 1:nrow(B_sub), labels = rownames(B_sub), las = 1, cex.axis = 0.5, tick = FALSE)
  box()
  dev.off()
}

pdf(paste0(res_dir, "/factor_ppi.pdf"), width = 1 * Q + 2, height = 3.5)
par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 1))
bp <- barplot(res$factor_ppi, names.arg = 1:Q, ylim = c(0, 1),
              col = ifelse(res$factor_ppi > 0.5, "grey30", "grey80"), border = NA,
              xlab = "Factor", ylab = "Posterior inclusion probability",
              main = "Factor activation")
abline(h = 0.5, lty = 2)
text(bp, pmin(res$factor_ppi + 0.05, 0.97), format(res$factor_ppi, digits = 2), cex = 0.7)
mtext(paste0("genes with PPI > ", thres, ": ",
             paste0(colSums(res$ppi > thres), collapse = ", ")), side = 1, line = 3.3, cex = 0.7)
dev.off()


## ---------------------------------------------------------------------------
## 12b. Scores against the number of observations
## ---------------------------------------------------------------------------

n_obs_Y <- sapply(time_obs, length)[names(Y)]

pdf(paste0(res_dir, "/abs_scores_vs_n_obs.pdf"), height = 3.2*Q, width = 6.5)
par(mfrow = c(Q, 2), pty = "s")
for (q in 1:Q) {
  for (l in 1:2) {
    abs_scores <- abs(res$list_Zeta_hat[[q]][names(Y), l])
    ct <- cor.test(n_obs_Y, abs_scores, method = "spearman", exact = FALSE)
    plot(jitter(n_obs_Y), abs_scores, pch = 20,
         xlab = "Number of observations", ylab = "|Score|",
         main = paste0("Factor q = ", q, ", FPC ", l),
         sub = paste0("Spearman rho = ", format(ct$estimate, digits = 2),
                      ", p = ", format(ct$p.value, digits = 2)))
  }
}
dev.off()


## ---------------------------------------------------------------------------
## 12c. Scores against the extraction batch
## ---------------------------------------------------------------------------

df_first <- df_comb_no_HC[order(df_comb_no_HC$days_from_sx_or_swab_imputed), ]
df_first <- df_first[!duplicated(df_first$subject_id), ]
batch_subj <- droplevels(as.factor(df_first$extractbatch[match(names(Y), df_first$subject_id)]))

pdf(paste0(res_dir, "/scores_vs_extractbatch.pdf"), height = 3.2*Q, width = 6.5)
par(mfrow = c(Q, 2))
for (q in 1:Q) {
  for (l in 1:2) {
    boxplot(res$list_Zeta_hat[[q]][names(Y), l] ~ batch_subj,
            xlab = "Extraction batch", ylab = "Scores",
            main = paste0("Factor q = ", q, ", FPC ", l))
  }
}
dev.off()


## ---------------------------------------------------------------------------
## 12d. Scores against the per-sample length-bias index
## ---------------------------------------------------------------------------

gene_length <- df_genes_sel[gene_id_of_label[rownames(res$ppi)], "Length"]
names(gene_length) <- rownames(res$ppi)
log_length <- log10(gene_length)

E_analysed <- t(as.matrix(df_comb_no_HC[, all_var, drop = FALSE]))
colnames(E_analysed) <- df_comb_no_HC$sample_id_v2

dev_mat <- E_analysed - rowMeans(E_analysed)
lb <- log_length[rownames(dev_mat)] - mean(log_length[rownames(dev_mat)], na.rm = TRUE)
ok <- !is.na(lb)

length_bias <- colSums(dev_mat[ok, ] * lb[ok]) /
  sqrt(colSums(dev_mat[ok, ]^2) * sum(lb[ok]^2))
lb_subj <- tapply(length_bias, df_comb_no_HC$subject_id, mean)[names(Y)]

pdf(paste0(res_dir, "/scores_vs_length_bias.pdf"), height = 3.2*Q, width = 6.5)
par(mfrow = c(Q, 2), pty = "s")
for (q in 1:Q) {
  for (l in 1:2) {
    plot(lb_subj, res$list_Zeta_hat[[q]][names(Y), l], pch = 20,
         xlab = "Length-bias index", ylab = "Scores",
         main = paste0("Factor q = ", q, ", FPC ", l))
  }
}
dev.off()


## ---------------------------------------------------------------------------
## 12g. Leading eigenfunctions
## ---------------------------------------------------------------------------

pdf(paste0(res_dir, "/leading_eigenfunctions.pdf"), height = 3 * Q, width = 6.5)
par(mfrow = c(Q, 2))
for (q in 1:Q) {
  for (l in 1:min(2, L)) {
    plot(res$time_g * days_thres, res$list_list_Phi_hat[[q]][, l], type = "l",
         xlab = "Days from symptom onset", ylab = "Eigenfunction",
         main = paste0("Factor ", q, ", l = ", l, " (",
                       format(res$list_cumulated_pve[[q]][l], digits = 3), "% cum PVE)"))
    abline(h = 0, lty = 2)
  }
}
dev.off()

cat("Figures written to ", res_dir, "\n", sep = "")
