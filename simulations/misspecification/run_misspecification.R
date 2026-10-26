## ---------------------------------------------------------------------------
## Robustness of bayesSYNC under model misspecification, Section 5.5 of the
## paper (Table 2 and the residual correlogram reported in the text).
##
## Four departures from the model are studied: errors correlated along t,
## unmodelled smooth variation at the latent level, correlation between the
## latent processes, and informative drop-out.
##
## The first two mechanisms correspond to the two blocks of the marginal
## covariance implied by the model,
##
##   Cov(y_i^(j)(s), y_i^(j')(t))
##     = sum_q b_jq b_j'q sigma^(q)(s, t) + sigma^2_eps,j 1{j = j'} delta(s - t).
##
##   Scenario 1a, errors correlated along t. A variable-specific error process
##     with the same marginal variance as the white noise it replaces, so only
##     the temporal structure of the diagonal blocks changes. This is the
##     error term that the latent factor model does not include, and it leaves
##     the off-diagonal blocks, which identify the loadings and the latent
##     processes, untouched. The range rho is varied, rho = 0 giving the
##     well-specified baseline through the identical code path.
##
##   Scenario 1b, unmodelled smooth variation at the latent level. A smooth
##     perturbation of the latent trajectories, representing the tail of the
##     Karhunen-Loeve expansion discarded by the truncation at L. It reaches
##     every variable through the loadings and therefore perturbs the
##     off-diagonal blocks, where it competes with the factor structure. Its
##     variance is varied as a proportion kappa of the variance of the
##     structured latent part.
##
##   Scenario 1c, correlation between the latent processes. The scores of
##     different factors are made correlated while each factor keeps its own
##     eigenfunctions. Mixing the trajectories themselves would be absorbed
##     exactly, being the same model with different loadings; correlating the
##     scores of factors with differing eigenfunctions is the case the model
##     cannot represent whatever the number of factors, and is therefore the one
##     worth simulating.
##
##   Scenario 2, informative drop-out. Subjects with more severe latent
##     trajectories stop being observed earlier, so that truncation is
##     informative (missing not at random). A non-informative drop-out arm
##     calibrated to the same expected retention is included as a control, so
##     that the effect of informativeness is separated from the effect of
##     simply having fewer observations.
##
## The fitting call is identical in all settings: the model is always fitted
## under its own assumptions. Only the data-generating mechanism changes.
##
## Replicates are paired across settings: for a given replicate index the
## loadings, scores, eigenfunctions, mean functions and observation times are
## identical across all settings, and only the misspecification is switched on.
## Differences between settings are then paired, which reduces the Monte Carlo
## error of the comparison.
##
## One array task fits one (setting, replicate) pair: with 9 settings and 50
## replicates, job_id runs from 1 to 450, settings varying fastest (see below).
##
## Meant for a cluster, one array task per job
## (slurm/submit_misspecification.sh): each job was given 16 cores, 32 Gb of
## memory and up to 12 hours, over 450 tasks. Setting bool_quick_test checks
## that the pipeline goes through on a laptop, but its results are not
## meaningful (see the comment at bool_quick_test).
##
## Usage, from the root of the repository:
##   Rscript simulations/misspecification/run_misspecification.R <job_id>
##
## Writes, in its own folder under out_dir: summary_row.RData (the metrics,
## read by gather_misspecification.R), output.RData (the fit, without Y) and a
## few diagnostic plots.
## ---------------------------------------------------------------------------

rm(list = ls())

source("R/paths.R")
source("R/plot_utils.R")
source("R/utils_sim.R")
source("R/utils_misspec.R")

require(bayesSYNC)

job_id <- get_job_id(default = 1)


## ---------------------------------------------------------------------------
## Settings
##
## The number of variables is reduced relative to the main simulation study.
## Neither mechanism studied here involves p: the correlated error process acts
## along t within a variable, and the drop-out mechanism acts on the observation
## times. Lowering p buys the replicates and the settings needed to make the
## comparisons precise.
## ---------------------------------------------------------------------------

# small run to check that the pipeline goes through. At p = 100 each factor
# has only about 15 nonzero loadings and the model recovers one factor of three
# in the baseline, so a quick test says whether the run completes, not how the
# settings compare.
bool_quick_test <- F

p <- 1000
N <- 100
N_min <- 5
N_max <- 10

Q_true <- 3
Q <- 5

L_true <- 2
L <- 5

n_repl <- 50 # replicates per setting, spread over the job array

tol_abs <- 1e-3
tol_rel <- 1e-6
maxit <- 10000

vec_sd_zeta <- 1/(1:L_true)
vec_sd_eps <- rep(1, p)

a_om <- 1
b_om <- 10

K <- NULL

col <- "skyblue1"

bool_anneal <- T
if (bool_anneal) {
  anneal <- c(1, 1.9, 100)
} else {
  anneal <- NULL
}

bool_var_spec_prob <- FALSE

n_cpus <- 16 # must match --cpus-per-task in slurm/submit_misspecification.sh
bool_save <- T

if (bool_quick_test) {
  p <- 100
  N <- 30
  maxit <- 500
  n_cpus <- 4
  vec_sd_eps <- rep(1, p)
  # n_repl is deliberately left unchanged here: it only bounds the admissible
  # job identifiers and plays no part in the run itself, so lowering it would
  # restrict which settings can be reached locally without making the test faster
}


## ---------------------------------------------------------------------------
## Misspecification settings
##
## Row 1 is the common baseline: independent errors and complete follow-up. It
## is shared by both scenarios and is generated through the same code path as
## the misspecified settings.
##
## Every mechanism is run at two severities against the common baseline, which
## is enough to show both whether a departure matters and whether its effect
## grows with the departure, without multiplying the settings.
##
## Scenario 1a varies the range of the error correlation. With n_i between 5 and
## 10 observations on [0, 1] the typical gap between consecutive times is about
## 0.1 to 0.2, so rho = 0.15 is a moderate departure (correlation of about 0.5
## at a typical gap) and rho = 0.40 a severe one (about 0.8), where the error
## process starts to look like a subject-level smooth curve. A milder level was
## dropped as dominated: if the severe setting is benign, the mild one is too.
##
## Scenario 1b varies the variance of the unmodelled latent process, as a
## proportion kappa of the variance of the structured latent part, at a fixed
## range. The range is short relative to the leading eigenfunctions, which vary
## on the scale of the whole domain, since the discarded components of a
## Karhunen-Loeve expansion carry higher frequencies than the leading ones.
##
## Scenario 1c varies the correlation between the scores of any two factors at a
## given component, rho_fac = 0.3 being a moderate departure and 0.6 a strong
## one. The marginal variance of every score is unchanged, so only the
## dependence across factors differs from the baseline.
##
## Scenario 2 contrasts non-informative drop-out with strongly informative
## drop-out at a fixed expected retention, so that the effect of informativeness
## is separated from the effect of having fewer observations. The non-informative
## arm is the control and cannot be dropped; an intermediate informative level
## was removed, the contrast that carries the argument being control against
## strong dependence.
## ---------------------------------------------------------------------------

target_retention <- 0.6 # expected proportion of retained observations in both drop-out arms
q_dropout <- 1          # latent factor driving drop-out
min_obs <- 1            # minimum number of observations retained per subject;
# only used by the drop-out settings, so the baseline and
# the other settings keep the full schedule

rho_lat <- 0.15         # range of the unmodelled latent process, scenario 1b
corr_type_lat <- "matern32"

df_settings <- rbind(
  data.frame(scenario = "baseline",
             rho = 0, corr_type = "exponential", prop_corr = 1, kappa = 0,
             dropout = "none", beta_dropout = 0, rho_fac = 0,
             stringsAsFactors = FALSE),
  data.frame(scenario = "corr_error",
             rho = c(0.15, 0.40), corr_type = "exponential", prop_corr = 1,
             kappa = 0, dropout = "none", beta_dropout = 0, rho_fac = 0,
             stringsAsFactors = FALSE),
  data.frame(scenario = "latent_error",
             rho = 0, corr_type = "exponential", prop_corr = 1,
             kappa = c(0.1, 0.3), dropout = "none", beta_dropout = 0,
             rho_fac = 0,
             stringsAsFactors = FALSE),
  data.frame(scenario = "latent_corr",
             rho = 0, corr_type = "exponential", prop_corr = 1,
             kappa = 0, dropout = "none", beta_dropout = 0,
             rho_fac = c(0.3, 0.6),
             stringsAsFactors = FALSE),
  data.frame(scenario = "dropout",
             rho = 0, corr_type = "exponential", prop_corr = 1, kappa = 0,
             dropout = c("non_informative", "informative"),
             beta_dropout = c(0, 2), rho_fac = 0,
             stringsAsFactors = FALSE)
)

n_settings <- nrow(df_settings)

# the job array indexes settings within replicates, so that a partial array
# still yields complete replicates across all settings: job_id = 1, ...,
# n_settings covers every setting for replicate 1, job_id = n_settings + 1, ...,
# 2 * n_settings covers every setting for replicate 2, and so on. Locally, the
# whole set of settings is therefore covered by job_id = 1, ..., n_settings
setting_id <- ((job_id - 1) %% n_settings) + 1
repl_id <- ((job_id - 1) %/% n_settings) + 1

n_job <- n_settings * n_repl
if (is.na(job_id) || job_id < 1 || job_id > n_job) {
  stop(paste0("job_id must be an integer between 1 and n_settings * n_repl = ",
              n_job, "; the value received is ", job_id, "."))
}

setting <- df_settings[setting_id, , drop = FALSE]

# the data seed depends on the replicate only, hence the pairing across settings
seed <- repl_id
seed_noise <- 100000 + repl_id   # separate stream for the error process
seed_dropout <- 200000 + repl_id # separate stream for the drop-out draws
seed_latent <- 300000 + repl_id  # separate stream for the latent perturbation
seed_fac <- 400000 + repl_id     # separate stream for the factor correlation

set.seed(seed)

print(paste0("Setting ", setting_id, " of ", n_settings, ", replicate ", repl_id))
print(setting)


## ---------------------------------------------------------------------------
## Output directory
## ---------------------------------------------------------------------------

if (bool_save) {
  res_dir <- paste0(out_dir, "/R1_misspec_",
                    setting$scenario,
                    "_rho_", setting$rho, "_", setting$corr_type,
                    "_kappa_", setting$kappa, "_rho_lat_", rho_lat,
                    "_rhofac_", setting$rho_fac,
                    "_drop_", setting$dropout, "_beta_", setting$beta_dropout,
                    "_ret_", target_retention,
                    "_anneal_", paste0(anneal, collapse = "-"),
                    "_N_", N, "_p_", p, "_Q_true_", Q_true, "_Q_", Q,
                    "_L_true_", L_true, "_L_", L,
                    "_N_min_", N_min, "_N_max_", N_max,
                    "_seed_", seed, "/")
  dir.create(res_dir, recursive = TRUE)
  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"), type = "message")
} else {
  res_dir <- NULL
}


## ---------------------------------------------------------------------------
## Data generation
##
## The generator is called with a zero error standard deviation, so that $Y
## holds the noise-free signal. The unmodelled latent process is added first,
## since it belongs to the trajectories, and the observation error is added
## last. This keeps the signal identical across settings for a given replicate
## and leaves generate_bayesSYNC_data() unchanged.
## ---------------------------------------------------------------------------

# different mean functions across variables, to avoid irrelevant factors
# compensating for the mean structure when Q > Q_true
mu_func <- function(t, p, j = 1, alpha = 3) (-1)^j*alpha*sin((2*pi+j %% 5)*t + j / p)

n <- round(runif(N, N_min, N_max))

fpca_data <- generate_bayesSYNC_data(N = N, n = n, p = p, Q = Q_true, L = L_true,
                                     mu_func = mu_func,
                                     a_om = a_om, b_om = b_om,
                                     vec_sd_zeta = vec_sd_zeta,
                                     # the error process is added separately below;
                                     # if the generator rejects a zero standard
                                     # deviation, use rep(1e-8, p) instead, which
                                     # is numerically negligible here
                                     vec_sd_eps = rep(0, p),
                                     seed = seed)

latent_res <- add_latent_error(fpca_data, Y = fpca_data$Y,
                               kappa = setting$kappa,
                               rho_lat = rho_lat,
                               corr_type = corr_type_lat,
                               seed = seed_latent)

vec_tau_latent <- latent_res$vec_tau      # NULL when kappa = 0
vec_var_h_latent <- latent_res$vec_var_h  # NULL when kappa = 0

if (setting$kappa > 0) {
  print(paste0("Unmodelled latent process, standard deviation per factor: ",
               paste0(round(vec_tau_latent, 3), collapse = ", "),
               "; structured latent standard deviation per factor: ",
               paste0(round(sqrt(vec_var_h_latent), 3), collapse = ", ")))
}

corr_res <- add_latent_correlation(fpca_data, Y = latent_res$Y,
                                   rho_fac = setting$rho_fac,
                                   seed = seed_fac)

if (setting$rho_fac != 0) {
  print(paste0("Correlation between factor scores: ", setting$rho_fac,
               "; realised correlation between factors 1 and 2, component 1: ",
               round(cor(corr_res$Zeta_corr[[1]][, 1],
                         corr_res$Zeta_corr[[2]][, 1]), 3)))
}

Y <- add_correlated_noise(fpca_data,
                          vec_sd_eps = vec_sd_eps,
                          rho = setting$rho,
                          corr_type = setting$corr_type,
                          prop_corr = setting$prop_corr,
                          Y = corr_res$Y,
                          seed = seed_noise)

time_obs <- fpca_data$time_obs
# under scenario 1c the scores actually generating the data are the correlated
# ones, so they, and not fpca_data$Zeta, are what the estimates are compared to
Zeta <- corr_res$Zeta_corr
B <- fpca_data$B
pat <- fpca_data$gamma
time_g <- fpca_data$time_g
mu_g <- fpca_data$mu_g
Phi_g <- fpca_data$Phi_g


## ---------------------------------------------------------------------------
## Drop-out
##
## The hazard of missing visit k depends on the value of the latent factor at
## visit k, that is, on the value that would have been recorded. Values are
## standardised by the pooled standard deviation of the latent paths, so that
## beta_dropout is interpretable as a log-odds ratio per standard deviation of
## latent severity and is comparable across replicates.
## ---------------------------------------------------------------------------

if (setting$dropout == "none") {
  
  dropout_res <- NULL
  vec_dropout <- rep(FALSE, N)
  retention <- 1
  n_obs_used <- sapply(time_obs, length)
  
} else {
  
  list_h <- compute_factor_path_obs(fpca_data, q_dropout = q_dropout)
  sd_h <- sd(unlist(list_h))
  list_s <- lapply(list_h, function(h_i) h_i / sd_h)
  
  alpha_dropout <- calibrate_dropout_intercept(list_s,
                                               beta = setting$beta_dropout,
                                               target_retention = target_retention)
  
  print(paste0("Calibrated drop-out intercept: ", round(alpha_dropout, 3),
               " (expected retention ",
               round(expected_retention(list_s, alpha_dropout, setting$beta_dropout), 3), ")"))
  
  dropout_res <- apply_dropout(fpca_data, Y, list_s,
                               alpha = alpha_dropout,
                               beta = setting$beta_dropout,
                               min_obs = min_obs,
                               seed = seed_dropout)
  
  time_obs <- dropout_res$time_obs
  Y <- dropout_res$Y
  vec_dropout <- dropout_res$vec_dropout
  retention <- dropout_res$retention
  n_obs_used <- dropout_res$n_obs_new
  
  print(paste0("Realised retention: ", round(retention, 3),
               ", subjects with at least one lost visit: ", sum(vec_dropout),
               " out of ", N))
  print(summary(n_obs_used))
}


## ---------------------------------------------------------------------------
## Inference, identical to the well-specified runs
## ---------------------------------------------------------------------------

runtime_bayesSYNC <- system.time(fpca_res <- bayesSYNC(time_obs = time_obs,
                                                       Y = Y, L = L, Q = Q,
                                                       K = K, anneal = anneal,
                                                       tol_abs = tol_abs,
                                                       tol_rel = tol_rel,
                                                       maxit = maxit,
                                                       seed = seed,
                                                       bool_var_spec_prob = bool_var_spec_prob,
                                                       n_cpus = n_cpus,
                                                       bool_scale = F))

print(runtime_bayesSYNC)

ppi <- fpca_res$ppi
factor_ppi <- fpca_res$factor_ppi
B_hat <- fpca_res$B_hat
list_cumulated_pve <- fpca_res$list_cumulated_pve
elbo <- fpca_res$ELBO_iter

list_Y_hat <- fpca_res$list_Y_hat
list_Y_low <- fpca_res$list_Y_low
list_Y_upp <- fpca_res$list_Y_upp
list_mu_hat <- fpca_res$list_mu_hat
list_list_Phi_hat <- fpca_res$list_list_Phi_hat
list_Zeta_hat <- fpca_res$list_Zeta_hat
list_Cov_zeta_hat <- fpca_res$list_Cov_zeta_hat

# number of factors retained under the median probability model, computed
# before the permutation so that it refers to the full set of Q candidates
n_retained <- sum(factor_ppi > 0.5)
print(paste0("Number of retained factors: ", n_retained, " (true value ", Q_true, ")"))


## ---------------------------------------------------------------------------
## Alignment and metrics
##
## Same permutation, sign matching and rescaling as in the well-specified runs,
## so that all quantities are on the scale of the simulated ones.
## ---------------------------------------------------------------------------

perm_fpca_res <- match_factor_and_sign(B, B_hat, ppi, factor_ppi, Zeta,
                                       list_Zeta_hat, list_list_Phi_hat,
                                       list_cumulated_pve, list_Cov_zeta_hat)

B_hat <- perm_fpca_res$perm_B_hat
ppi <- perm_fpca_res$perm_ppi
factor_ppi <- perm_fpca_res$perm_factor_ppi

print(paste0("Factor PPI after permutation based on true loadings: ",
             paste0(format(factor_ppi, digits = 2), collapse = ", ")))

list_list_Phi_hat <- perm_fpca_res$perm_list_list_Phi_hat
list_Zeta_hat <- perm_fpca_res$perm_list_Zeta_hat
list_cumulated_pve <- perm_fpca_res$perm_list_cumulated_pve
list_Cov_zeta_hat <- perm_fpca_res$perm_list_Cov_zeta_hat

if (Q != Q_true) {
  list_Zeta_hat_untrimmed <- perm_fpca_res$perm_list_Zeta_hat_untrimmed
  list_list_Phi_hat_untrimmed <- perm_fpca_res$perm_list_list_Phi_hat_untrimmed
  B_hat_untrimmed <- perm_fpca_res$perm_B_hat_untrimmed
  ppi_untrimmed <- perm_fpca_res$perm_ppi_untrimmed
} else {
  list_Zeta_hat_untrimmed <- list_list_Phi_hat_untrimmed <- B_hat_untrimmed <- ppi_untrimmed <- NULL
}

# metrics already used in the well-specified simulations
Corr_zeta <- corr_zeta_vector(Zeta, list_Zeta_hat)
Corr_b <- corr_b_vector(B, B_hat)
Auc_b <- auc_roc_curve(B, ppi)
Err_mean <- error_mean_function(time_g, list_mu_hat, mu_g)
Err_eigen <- error_eigen_function(time_g, Phi_g, C_g = NULL, list_list_Phi_hat)

# metrics specific to the misspecification study: bias and attenuation of the
# scores, coverage of the score credible intervals, and the same quantities
# stratified by drop-out status
score_metrics <- score_bias_metrics(Zeta, list_Zeta_hat,
                                    list_Cov_zeta_hat = list_Cov_zeta_hat,
                                    group = vec_dropout)

rt <- runtime_bayesSYNC[["elapsed"]] / 60

# diagnostics on the fitted error process, computed from the unpermuted output
sigsq_eps_hat <- extract_sigsq_eps(fpca_res)
print(paste0("Estimated error variance, mean over variables: ",
             round(mean(sigsq_eps_hat), 3), " (true value ", vec_sd_eps[1]^2, ")"))

df_correlogram <- residual_correlogram(fpca_res, time_obs,
                                       breaks = seq(0, 1, by = 0.05),
                                       max_p = min(p, 1000), seed = seed)
print(head(df_correlogram, 5))

# a factor whose estimated loadings are essentially zero is left unrescaled by
# match_factor_and_sign(), which would put its scores and covariances on a
# different scale from the simulated ones; the smallest column norm flags it
min_norm_col_B_hat <- min(sqrt(colSums(B_hat[, 1:Q_true, drop = FALSE]^2)))

alignment_ok <- alignment_diagnostic(Zeta, list_Zeta_hat)
if (!alignment_ok) {
  warning(paste0("The component alignment is not the identity for at least one ",
                 "component; the score covariances are permuted by factor only, ",
                 "so the coverage metric should be interpreted with care for ",
                 "this run."))
}

summary_row <- build_summary_row(setting = setting, seed = seed, p = p,
                                 Corr_zeta = Corr_zeta, Corr_b = Corr_b,
                                 Auc_b = Auc_b, Err_mean = Err_mean,
                                 Err_eigen = Err_eigen, n_retained = n_retained,
                                 score_metrics = score_metrics,
                                 retention = retention,
                                 mean_n_obs = mean(n_obs_used),
                                 prop_dropout = mean(vec_dropout),
                                 rt = rt,
                                 sigsq_eps_hat = sigsq_eps_hat,
                                 df_correlogram = df_correlogram,
                                 alignment_ok = alignment_ok,
                                 min_norm_col_B_hat = min_norm_col_B_hat)

print(t(summary_row))


## ---------------------------------------------------------------------------
## Saving
##
## The one-row summary is saved separately from the full output, so that the
## across-job gathering step is a plain rbind of small files.
## ---------------------------------------------------------------------------

if (bool_save) {
  
  save(summary_row, file = file.path(res_dir, "summary_row.RData"))
  
  res <- create_named_list(setting, seed, Corr_zeta, Corr_b, Auc_b,
                           Err_mean, Err_eigen, score_metrics, n_retained,
                           factor_ppi, ppi, ppi_untrimmed, B, B_hat, B_hat_untrimmed,
                           pat, Zeta, list_Zeta_hat, list_Zeta_hat_untrimmed,
                           list_Cov_zeta_hat,
                           mu_g, list_mu_hat, Phi_g, list_list_Phi_hat,
                           list_list_Phi_hat_untrimmed,
                           list_cumulated_pve, elbo, rt,
                           sigsq_eps_hat, df_correlogram, alignment_ok,
                           min_norm_col_B_hat,
                           rho_lat, corr_type_lat,
                           vec_tau_latent, vec_var_h_latent,
                           time_obs, n_obs_used, vec_dropout, retention,
                           dropout_res)
  
  save(res, file = file.path(res_dir, "output.RData")) # Y is not stored, to keep the files small
}


## ---------------------------------------------------------------------------
## Diagnostic plots for a single run
##
## Kept minimal here; the across-setting summaries are produced by the
## gathering script.
## ---------------------------------------------------------------------------

if (bool_save) {
  
  pdf(paste0(res_dir, "/eigenfunctions.pdf"), width = 3*L, height = 3*Q_true)
  par(mfrow = c(Q_true, L))
  null <- sapply(1:Q_true, function(q)
    sapply(1:L, function(l) {
      if (l <= L_true) {
        res_Psi_g <- Phi_g
      } else {
        res_Psi_g <- NULL
      }
      plot_eigenfunctions(q, l, time_g, Psi_g = res_Psi_g, list_list_Phi_hat,
                          list_cumulated_pve, cum_pve_thres = 0.99,
                          col_sim = "black", col = col)
    })
  )
  dev.off()
  
  pdf(paste0(res_dir, "/scores_simulated_vs_estimated.pdf"), width = 4*L_true, height = 3*Q_true)
  par(mfrow = c(Q_true, L_true))
  for (q in 1:Q_true) {
    for (l in 1:L_true) {
      # droppers are highlighted, to show the direction of the bias induced by
      # informative truncation
      plot(Zeta[[q]][, l], list_Zeta_hat[[q]][, l],
           col = ifelse(vec_dropout, "red", "black"),
           pch = 19, cex = 0.7,
           xlab = "Simulated", ylab = "Estimated",
           main = paste0("Scores, q = ", q, ", l = ", l))
      abline(a = 0, b = 1, lty = 2)
      abline(lm(list_Zeta_hat[[q]][, l] ~ Zeta[[q]][, l]), col = "blue")
    }
  }
  dev.off()
  
  pdf(paste0(res_dir, "/residual_correlogram.pdf"), width = 5, height = 4.5)
  par(mfrow = c(1, 1))
  plot(df_correlogram$dist_mid, df_correlogram$corr, type = "o", pch = 19,
       xlab = "Lag between observation times",
       ylab = "Residual correlation",
       ylim = range(c(0, df_correlogram$corr), na.rm = TRUE),
       main = "Residual correlation along t")
  abline(h = 0, lty = 2)
  if (setting$rho > 0) {
    # correlation of the simulated error process, for reference; no reference
    # curve is drawn for the latent perturbation, whose imprint on the residuals
    # is filtered by the fit and is not the correlation of the process itself
    curve(exp(-x / setting$rho), add = TRUE, col = "red", lty = 3)
  }
  dev.off()
  
  if (!is.null(dropout_res)) {
    pdf(paste0(res_dir, "/dropout_diagnostics.pdf"), width = 9, height = 4.5)
    par(mfrow = c(1, 2))
    plot(Zeta[[q_dropout]][, 1], dropout_res$n_obs_new,
         xlab = paste0("Simulated score, factor ", q_dropout, ", component 1"),
         ylab = "Number of retained observations",
         pch = 19, cex = 0.7,
         main = "Drop-out against latent severity")
    hist(dropout_res$vec_time_dropout,
         xlab = "Time of the first lost visit", main = "Drop-out times")
    dev.off()
  }
}

if (bool_save) {
  sink()
  sink(type = "message")
}