## ---------------------------------------------------------------------------
## MEFISTO-base on the biomarker panel data (Supplementary Material S8).
##
## MEFISTO-base estimates one latent trajectory per factor, shared by all
## patients. The subject-specific version (patients as groups) is not run: in the
## simulations of Section 5.4 it did not complete within 36 hours at N = 40 in
## 11 of 25 runs, the application has N = 88, and the current implementation
## with patients as groups requires a common time grid or time warping, which
## the irregular sampling of the application does not provide.
##
## The input is the data of the fit reported in Section 6.1 (N = 88 patients, 57
## markers on the real data), taken from the output.RData of that fit, on the same standardised
## scale (res$Y) and with the same number of candidate factors (Q = 5), so that
## the two methods see identical data.
##
## Usage, from the root of the repository:
##   Rscript application/biomarker_panel/mefisto_base_biomarker_panel.R <run_dir>
## with the folder written by analysis_biomarker_panel.R.
##
## Outputs, written to mefisto_base_biomarker_panel under out_dir:
##   mefisto_base_covid.hdf5            the fitted MEFISTO-base model
##   comparison_mefisto_covid.txt       log: variance explained, smoothness
##                                      scales, factor matching, severity
##   comparison_mefisto_loadings.csv    loadings of both methods, matched factors
##   comparison_mefisto_loadings.pdf    side-by-side loadings of matched factors
##   mefisto_shared_trajectories.pdf    MEFISTO's shared factor trajectories
##
## Runs locally in a few minutes. MOFA2 needs its Python backend: with
## use_basilisk = TRUE it is installed by basilisk on first use, in the user's
## cache; set use_basilisk = FALSE if mofapy2 is already available to R.
## ---------------------------------------------------------------------------

rm(list = ls())

use_basilisk <- TRUE

source("R/paths.R")

source("R/utils_sim.R") # reshape_Y_to_df()
require(MOFA2)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1, dir.exists(args[1]))
fit_dir <- args[1]
res_dir_m <- file.path(out_dir, "mefisto_base_biomarker_panel")
dir.create(res_dir_m, showWarnings = FALSE)

seed <- 1
n_top <- 15 # markers per factor in the top-N overlap


## ---------------------------------------------------------------------------
## Data of the reported fit
## ---------------------------------------------------------------------------

e <- new.env()
load(file.path(fit_dir, "output.RData"), envir = e)
res <- e$res
time_obs <- e$time_obs
days_scale <- e$days_thres - 1 # times are rescaled by days_thres - 1 in prepare_data_irregular()
Q <- ncol(res$B_hat)
p <- nrow(res$B_hat)
subject_ids <- names(time_obs)
severity <- as.character(e$df_info$severity[match(subject_ids, e$df_info$subject_id)])
stopifnot(identical(names(res$Y[[1]]), rownames(res$B_hat)), !any(is.na(severity)))

# standardised data, as seen by our model
Y_std <- res$Y
names(Y_std) <- subject_ids
df_mofa <- reshape_Y_to_df(Y_std, time_obs)


## ---------------------------------------------------------------------------
## MEFISTO-base: no groups, time as the smooth covariate
## ---------------------------------------------------------------------------

outfile <- file.path(res_dir_m, "mefisto_base_covid.hdf5")
time_grid <- seq(0, 1, length.out = 200)

# the model is fitted only if it is not already on disk, so that the comparison
# and the figures can be redone without the Python backend
if (!file.exists(outfile)) {

sm_untrained <- create_mofa(data = df_mofa)
sm_untrained <- set_covariates(sm_untrained, covariates = "time")

data_opts <- get_default_data_options(sm_untrained)
model_opts <- get_default_model_options(sm_untrained)
model_opts$num_factors <- Q
train_opts <- get_default_training_options(sm_untrained)
train_opts$seed <- seed
mefisto_opts <- get_default_mefisto_options(sm_untrained)
mefisto_opts$new_values <- matrix(time_grid, nrow = 1)

sm_untrained <- prepare_mofa(sm_untrained,
                             model_options = model_opts,
                             mefisto_options = mefisto_opts,
                             training_options = train_opts,
                             data_options = data_opts)

runtime_mefisto <- system.time(sm <- run_mofa(sm_untrained, outfile = outfile,
                                              use_basilisk = use_basilisk))
} else {
  runtime_mefisto <- NA
}
# all factors kept, so that the matching sees every candidate
sm <- load_model(outfile, load_interpol_Z = TRUE, remove_inactive_factors = FALSE)

sink(file.path(res_dir_m, "comparison_mefisto_covid.txt"), split = TRUE)

cat("MEFISTO-base on the metabolomic data: N =", length(subject_ids), "patients,",
    length(unlist(time_obs)), "samples, p =", p, "markers, Q =", Q, "\n")
print(runtime_mefisto)

W <- sm@expectations$W[[1]]
W <- W[rownames(res$B_hat), , drop = FALSE]
colnames(W) <- paste0("M", seq_len(ncol(W)))

cat("\nVariance explained by each MEFISTO factor (R2, %):\n")
r2 <- tryCatch(get_variance_explained(sm)$r2_per_factor[[1]], error = function(err) NULL)
print(if (is.null(r2)) "not available" else round(r2, 2))

cat("\nSmoothness scale of each MEFISTO factor (1 = entirely smooth in time):\n")
sc <- tryCatch(get_scales(sm), error = function(err) NULL)
print(if (is.null(sc)) "not available" else round(sc, 3))


## ---------------------------------------------------------------------------
## Matching MEFISTO factors to the retained factors of our fit
##
## The first factor of the paper is the retained factor with the largest CRP
## loading in absolute value, the second is the other one.
## ---------------------------------------------------------------------------

retained <- which(res$factor_ppi > 0.5)
stopifnot(length(retained) == 2)
q_infl <- retained[which.max(abs(res$B_hat["log_CRP", retained]))]
q_amino <- setdiff(retained, q_infl)
ours <- c("first factor (inflammation)" = unname(q_infl), "second factor (amino acids)" = unname(q_amino))

C <- abs(cor(res$B_hat[, ours], W))
rownames(C) <- names(ours)
cat("\nAbsolute correlation between our loadings (rows) and MEFISTO's (columns):\n")
print(round(C, 3))

top_n <- function(v) names(sort(abs(v), decreasing = TRUE))[seq_len(n_top)]
match_m <- integer(0)
cat("\nBest MEFISTO match for each of our factors:\n")
for (k in seq_along(ours)) {
  o <- order(C[k, ], decreasing = TRUE)
  m <- o[1]
  match_m <- c(match_m, m)
  overlap <- length(intersect(top_n(res$B_hat[, ours[k]]), top_n(W[, m]))) / n_top
  cat(sprintf("  %s -> %s: |cor| %.3f (runner-up %.3f), top-%d overlap %.2f\n",
              names(ours)[k], colnames(W)[m], C[k, m], C[k, o[2]], n_top, overlap))
}
if (anyDuplicated(match_m)) {
  cat("NOTE: both of our factors match the same MEFISTO factor.\n")
}

# Display orientation of the paper: CRP loads positively on the first factor and
# the amino acids negatively on the second. MEFISTO's loadings are then oriented
# like ours. The same signs are applied to the factor values and trajectories.
amino <- intersect(c("Leucine", "Isoleucine", "Valine", "Methionine", "Tyrosine", "Lysine"), rownames(res$B_hat))
sgn_ours <- c(sign(res$B_hat["log_CRP", q_infl]), -sign(mean(res$B_hat[amino, q_amino])))
B_ours <- sweep(res$B_hat[, ours, drop = FALSE], 2, sgn_ours, "*")
sgn <- sign(sapply(seq_along(ours), function(k) cor(B_ours[, k], W[, match_m[k]])))
W_match <- sweep(W[, match_m, drop = FALSE], 2, sgn, "*")

df_load <- data.frame(marker = rownames(res$B_hat),
                      ours_first = B_ours[, 1], ppi_first = res$ppi[, q_infl],
                      mefisto_first = W_match[, 1],
                      ours_second = B_ours[, 2], ppi_second = res$ppi[, q_amino],
                      mefisto_second = W_match[, 2])
write.csv(df_load, file.path(res_dir_m, "comparison_mefisto_loadings.csv"), row.names = FALSE)


## ---------------------------------------------------------------------------
## Association with severity
##
## MEFISTO-base has one trajectory per factor shared by all patients, but its
## factor values at the samples can still differ between patients. They are
## averaged over each patient's samples and compared with the severity class,
## alongside our first-component scores.
## ---------------------------------------------------------------------------

Z <- sm@expectations$Z[[1]]
subj_idx <- as.integer(sub(".*_subject_", "", rownames(Z)))
Z_patient <- apply(Z[, match_m, drop = FALSE], 2, function(z) tapply(z, subj_idx, mean))
Z_patient <- Z_patient[as.character(seq_along(subject_ids)), , drop = FALSE]

cat("\nAssociation with severity class (Kruskal-Wallis p) and rank correlation with the class:\n")
for (k in seq_along(ours)) {
  z_m <- Z_patient[, k]
  z_o <- res$list_Zeta_hat[[ours[k]]][, 1]
  sev_rank <- match(severity, c("B", "C", "D", "E"))
  cat(sprintf("  %s: ours (FPC 1) p = %.2g, rho = %.2f | MEFISTO %s (patient mean) p = %.2g, rho = %.2f | rho(ours, MEFISTO) = %.2f\n",
              names(ours)[k],
              kruskal.test(z_o ~ factor(severity))$p.value, abs(cor(z_o, sev_rank, method = "spearman")),
              colnames(W)[match_m[k]],
              kruskal.test(z_m ~ factor(severity))$p.value, abs(cor(z_m, sev_rank, method = "spearman")),
              abs(cor(z_o, z_m, method = "spearman"))))
}
cat("(absolute rank correlations, since the sign of a factor is arbitrary)\n")

sink()


## ---------------------------------------------------------------------------
## Figures
## ---------------------------------------------------------------------------

# loadings of the matched factors, each column scaled to a maximum absolute value
# of one, since the two methods scale factors and loadings differently
L_show <- cbind(df_load$ours_first, df_load$mefisto_first, df_load$ours_second, df_load$mefisto_second)
L_show <- sweep(L_show, 2, apply(abs(L_show), 2, max), "/")
colnames(L_show) <- c("ours, 1", "MEFISTO, 1", "ours, 2", "MEFISTO, 2")
pal <- colorRampPalette(c("skyblue1", "white", "#ff7f6e"))(201)
pdf(file.path(res_dir_m, "comparison_mefisto_loadings.pdf"), width = 6, height = 11)
par(mar = c(6, 11, 2, 1))
image(seq_len(ncol(L_show)), seq_len(p), t(L_show[p:1, ]), col = pal, zlim = c(-1, 1),
      axes = FALSE, xlab = "", ylab = "")
axis(1, at = seq_len(ncol(L_show)), labels = colnames(L_show), las = 2)
axis(2, at = seq_len(p), labels = rev(gsub("log_", "", df_load$marker)), las = 1, cex.axis = 0.7)
abline(v = 2.5, lwd = 2)
box()
dev.off()

# MEFISTO's shared trajectories of the matched factors, with 95% bands
iz <- sm@interpolated_Z[[1]]
pdf(file.path(res_dir_m, "mefisto_shared_trajectories.pdf"), width = 8, height = 3.8)
par(mfrow = c(1, 2))
for (k in seq_along(ours)) {
  mu <- sgn[k] * iz$mean[match_m[k], ]
  s <- sqrt(iz$variance[match_m[k], ])
  days <- time_grid * days_scale
  plot(days, mu, type = "l", lwd = 2, ylim = range(c(mu - 1.96 * s, mu + 1.96 * s)),
       xlab = "Days from symptom onset", ylab = "Shared factor trajectory",
       main = paste0("MEFISTO-base, ", names(ours)[k]), cex.main = 0.9)
  lines(days, mu - 1.96 * s, lty = 2)
  lines(days, mu + 1.96 * s, lty = 2)
  abline(h = 0, col = "grey60", lty = 3)
}
dev.off()

cat("written to", res_dir_m, "\n")
