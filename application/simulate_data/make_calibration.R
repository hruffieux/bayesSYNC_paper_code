## ---------------------------------------------------------------------------
## Builds the two calibration files used by simulate_cohort_data.R from the
## fits of the paper to the real COVID-19 data. Included to document how the
## files were made; it cannot be run without the fits, which are not public.
##
## Only population-level quantities and design summaries are kept, namely
##
##   - the estimated loadings, mean functions and eigenfunctions of the
##     retained factors, and the standardisation constants of each variable;
##   - for each factor and component, the mean of the estimated scores within
##     each severity class and their pooled within-class standard deviation,
##     to which the average posterior variance of the scores is added (see
##     score_summaries below);
##   - the residual standard deviation of each variable;
##   - counts: patients per severity class, numbers of samples per patient and
##     sampling days by class, deaths within the window by class, samples by
##     extraction and sequencing batch;
##   - for the controls of the biomarker panel, the mean and standard deviation
##     of each marker, on the standardised scale of the fit; for the transcriptomic data, the mean and standard
##     deviation of log hsCRP by class;
##   - the names of the variables and, for the transcriptomic data, the public
##     gene annotation (Ensembl identifier, symbol, chromosome, length).
##
## No subject identifier and no subject-level value is kept.
##
## Usage, from the root of the repository:
##   Rscript application/simulate_data/make_calibration.R \
##     <biomarker panel output.RData> <transcriptomic run folder>
## with the output.RData of the N = 88 fit of analysis_biomarker_panel.R and
## the folder of the N = 164 fit of analysis_transcriptomics.R (holding
## output_slim.RData).
##
## Writes application/simulate_data/calibration/biomarker_panel_calibration.rds
## and application/simulate_data/calibration/transcriptomics_calibration.rds.
## ---------------------------------------------------------------------------

rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2, file.exists(args[1]), dir.exists(args[2]))

cal_dir <- file.path("application", "simulate_data", "calibration")
dir.create(cal_dir, showWarnings = FALSE, recursive = TRUE)

days_thres <- 50
classes <- c("B", "C", "D", "E")

# index of the grid point closest to each observation time
nearest_index <- function(time_g, t) {
  sapply(t, function(tt) which.min(abs(time_g - tt)))
}

# class means of the scores of each component and pooled within-class
# standard deviation, for one factor.
#
# The estimated scores are posterior means, shrunk towards zero, the more so for
# patients with few samples, so their spread understates that of the scores
# themselves. By the law of total variance, the variance of a score is the
# variance of its posterior mean plus the average of its posterior variance, so
# the average posterior variance (post_var, patients x components, NA where
# unavailable) is added to the within-class variance. Without it, a factor
# carried by weakly observed patients is simulated with too little signal.
score_summaries <- function(Zeta_hat, vec_class, post_var = NULL) {
  class_mean <- t(sapply(classes, function(cc) colMeans(Zeta_hat[vec_class == cc, , drop = FALSE])))
  resid <- Zeta_hat - class_mean[vec_class, , drop = FALSE]
  within_var <- colSums(resid^2) / (nrow(Zeta_hat) - length(classes))
  if (!is.null(post_var)) {
    mean_post_var <- colMeans(post_var)
    within_var <- within_var + ifelse(is.na(mean_post_var), 0, mean_post_var)
  }
  list(class_mean = class_mean, within_sd = sqrt(within_var))
}

# sampling design of the analysed patients, by class
design_summaries <- function(vec_class, time_obs) {
  n_obs <- sapply(time_obs, length)
  days <- lapply(time_obs, function(t) round(t * (days_thres - 1)))
  list(n_patients = table(factor(vec_class, levels = classes)),
       n_obs_by_class = table(factor(vec_class, levels = classes),
                              factor(n_obs, levels = 1:max(n_obs))),
       days_by_class = t(sapply(classes, function(cc) {
         table(factor(unlist(days[vec_class == cc]), levels = 0:(days_thres - 1)))
       })))
}


## ---------------------------------------------------------------------------
## Biomarker panel, N = 88, p = 57
## ---------------------------------------------------------------------------

m <- new.env()
load(args[1], envir = m)
res <- m$res

subj <- names(m$time_obs)
stopifnot(identical(subj, names(res$Y)))
df_subj <- m$df_info[match(subj, m$df_info$subject_id), ]
vec_class <- as.character(df_subj$severity)
stopifnot(all(vec_class %in% classes))

active <- which(res$factor_ppi > 0.5)
p <- length(m$all_var)

mu_g <- sapply(res$list_mu_hat, identity) # n_g x p, standardised scale
colnames(mu_g) <- m$all_var

B <- res$B_hat[, active, drop = FALSE]

# posterior variances of the scores. This fit predates the storage of the
# posterior covariances, but keeps the 95% credible ellipses of the first two
# components, built by ellipse::ellipse with radius sqrt(qchisq(0.95, 2)), so
# the half-range of an ellipse along a component is that radius times the
# posterior standard deviation. The other components are left without.
post_var_panel <- function(q) {
  L_q <- ncol(res$list_Zeta_hat[[q]])
  t(sapply(res$list_list_zeta_ellipse[[q]][subj], function(E) {
    v <- rep(NA, L_q)
    v[1:2] <- (apply(E, 2, function(x) diff(range(x))) / 2)^2 / qchisq(0.95, 2)
    v
  }))
}

# residual standard deviation of each marker, on the standardised scale, at
# the observation times
resid <- do.call(rbind, lapply(seq_along(subj), function(i) {
  idx <- nearest_index(res$time_g, m$time_obs[[i]])
  sapply(seq_len(p), function(j) res$Y[[i]][[j]] - res$list_Y_hat[[i]][[j]][idx])
}))
sd_eps <- sqrt(colMeans(resid^2))
names(sd_eps) <- m$all_var

# the workspace holds df_comb after the analysis script has standardised it
# (to draw the control bands on the scale of the fit), so these summaries of
# the controls are on the standardised scale
df_hc <- m$df_comb[m$df_comb$severity == "HC", m$all_var]

# the fit of the paper predates the vital-status code of the analysis script,
# so deaths within the window are recomputed from the clinical table as there
bool_dead <- df_subj$hospital_outcome %in% "dead" &
  !is.na(df_subj$hospital_outcome_days_from_start_sx_or_first_pos_swab) &
  df_subj$hospital_outcome_days_from_start_sx_or_first_pos_swab < days_thres
dead_in_window <- table(factor(vec_class, levels = classes),
                        factor(bool_dead, levels = c(FALSE, TRUE)))

cal_panel <- list(
  var_names = m$all_var,
  metabolite_names = names(m$df_ms),
  time_g = res$time_g,
  mu_g = round(mu_g, 5),
  B = round(B, 5),
  list_Phi = lapply(active, function(q) round(res$list_list_Phi_hat[[q]], 5)),
  list_scores = lapply(active, function(q) score_summaries(res$list_Zeta_hat[[q]][subj, , drop = FALSE],
                                                          vec_class, post_var_panel(q))),
  sd_eps = sd_eps,
  mean_var = res$mean_mean_across_subjects[m$all_var],
  sd_var = res$sd_mean_across_subjects[m$all_var],
  hc_mean_std = colMeans(df_hc, na.rm = TRUE),
  hc_sd_std = apply(df_hc, 2, sd, na.rm = TRUE),
  n_hc = nrow(df_hc),
  design = design_summaries(vec_class, m$time_obs),
  dead_in_window = dead_in_window,
  days_thres = days_thres)

saveRDS(cal_panel, file.path(cal_dir, "biomarker_panel_calibration.rds"), compress = "xz")
rm(m, res)


## ---------------------------------------------------------------------------
## Transcriptomic data, N = 164, p = 13,177
## ---------------------------------------------------------------------------

r <- new.env()
load(file.path(args[2], "output_slim.RData"), envir = r)
res <- r$res_slim

subj <- names(r$time_obs)
df_first <- r$df_comb_no_HC[order(r$df_comb_no_HC$days_from_sx_or_swab_imputed), ]
df_first <- df_first[match(subj, df_first$subject_id), ]
vec_class <- as.character(df_first$severity)
stopifnot(all(vec_class %in% classes))

genes <- r$all_var
p <- length(genes)
Q <- ncol(res$B_hat)

# the mean functions are smooth (K = 4 interior knots), so every fourth point
# of the grid of 200 is kept and the simulator interpolates linearly
ind_g <- unique(c(seq(1, length(res$time_g), by = 4), length(res$time_g)))
mu_g <- sapply(res$list_mu_hat[genes], function(mu) mu[ind_g])
colnames(mu_g) <- genes

B <- res$B_hat[genes, , drop = FALSE]

# residual standard deviation of each gene on the standardised scale; the
# saved Y is on the log-CPM scale and is standardised as in the fit
mean_var <- res$mean_mean_across_subjects[genes]
sd_var <- res$sd_mean_across_subjects[genes]
mu_full <- sapply(res$list_mu_hat[genes], identity)
ss <- rep(0, p)
n_tot <- 0
for (i in seq_along(subj)) {
  idx <- nearest_index(res$time_g, r$time_obs[[i]])
  Y_i <- do.call(cbind, r$Y[[i]][genes])
  Y_i <- sweep(sweep(matrix(Y_i, nrow = length(idx)), 2, mean_var), 2, sd_var, "/")
  H_i <- sapply(seq_len(Q), function(q) res$list_h_hat[[i]][[q]][idx])
  fitted <- mu_full[idx, , drop = FALSE] + matrix(H_i, nrow = length(idx)) %*% t(B)
  ss <- ss + colSums((Y_i - fitted)^2)
  n_tot <- n_tot + length(idx)
}
sd_eps <- sqrt(ss / n_tot)
names(sd_eps) <- genes

df_samples <- r$df_comb_no_HC[r$df_comb_no_HC$subject_id %in% subj, ]
log_crp <- tapply(log(df_samples$hsCRP_mg_l), as.character(df_samples$severity),
                  function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)))

dead_in_window <- table(factor(vec_class, levels = classes),
                        factor(r$df_status$dead_in_window[match(subj, r$df_status$subject_id)],
                               levels = c(FALSE, TRUE)))

df_genes <- r$df_genes_sel[match(genes, r$df_genes_sel$label),
                           c("GeneID", "Length", "Chr", "hgnc_symbol", "description")]
rownames(df_genes) <- NULL

cal_rnaseq <- list(
  genes = df_genes,
  time_g = res$time_g[ind_g],
  mu_g = round(mu_g, 4),
  B = round(B, 4),
  list_Phi = lapply(seq_len(Q), function(q) round(res$list_list_Phi_hat[[q]], 5)),
  time_g_Phi = res$time_g,
  list_scores = lapply(seq_len(Q), function(q) {
    post_var <- t(sapply(res$list_Cov_zeta_hat[[q]][subj], diag))
    score_summaries(res$list_Zeta_hat[[q]][subj, , drop = FALSE], vec_class, post_var)
  }),
  sd_eps = round(sd_eps, 4),
  mean_var = round(mean_var, 4),
  sd_var = round(sd_var, 4),
  log_crp = log_crp,
  batches = table(df_samples$extractbatch, df_samples$seqbatch),
  design = design_summaries(vec_class, r$time_obs),
  dead_in_window = dead_in_window,
  n_genes_before_filter = 22354,
  days_thres = days_thres)

saveRDS(cal_rnaseq, file.path(cal_dir, "transcriptomics_calibration.rds"), compress = "xz")

cat("Calibration files written to", cal_dir, "\n")
print(file.info(list.files(cal_dir, full.names = TRUE))[, "size", drop = FALSE])
