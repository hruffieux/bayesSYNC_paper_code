## ---------------------------------------------------------------------------
## Comparison with MEFISTO, Section 5.4 of the paper (Figure 5, Table 1 and
## Supplementary Material S4 and S5).
##
## p = 100 variables, N in {20, 25, 30, 35, 40} subjects, Q = 2 latent factors
## with L = 3 FPCA components each, 2 to 10 observations per subject, zero mean
## functions and non-sparse loadings (every factor influences every variable
## with probability drawn from a Beta(1, 1) distribution). Each simulated data
## set is fitted with bayesSYNC and with one of the two MEFISTO settings, all
## with Q_max = 5 candidate factors:
##
##   MEFISTO        one latent trajectory per subject and factor: each subject
##                  is a MEFISTO group and model_groups = TRUE
##   MEFISTO-base   one latent trajectory per factor shared by all subjects:
##                  model_groups = FALSE
##
## One array task fits one (setting, N, replicate) triple:
##
##   job_id   1 to 130   MEFISTO,      N = 20, 25, 30, 35, 40 in blocks of 26
##   job_id 131 to 260   MEFISTO-base, N = 20, 25, 30, 35, 40 in blocks of 26
##
## so that seed = ((job_id - 1) %% 26) + 1. The paper uses 25 replicates per N;
## a 26th seed is included because seed 7 failed for MEFISTO at N = 25 and
## seed 26 replaced it. The seeds entering each summary are listed in
## summarise_mefisto_comparison.R. The seed is that of the data generation, so
## the two MEFISTO settings see the same simulated data for a given N and seed.
##
## MEFISTO runs within the 36-hour limit of the cluster only for a part of the
## replicates at N = 40 (14 of 25). All runtimes are serial (n_cpus = 1).
##
## MOFA2 calls its Python backend (mofapy2). With use_basilisk = TRUE it is
## installed by basilisk on first use; on the cluster use_basilisk = FALSE was
## used, with mofapy2 installed in the Python environment seen by R.
##
## Meant for a cluster, one array task per job
## (slurm/submit_mefisto_comparison.sh): each job was given one core, 64 Gb of
## memory and up to 36 hours. It is not meant to be run on a laptop.
##
## Usage, from the root of the repository:
##   Rscript simulations/mefisto_comparison/run_mefisto_comparison.R <job_id>
##
## Writes, in its own folder under out_dir: output.RData (fits, simulated data
## and metrics), the MEFISTO model as hdf5, and out.txt.
## ---------------------------------------------------------------------------

rm(list = ls())

source("R/paths.R")
source("R/plot_utils.R")
source("R/utils_sim.R")

require(bayesSYNC)
require(MOFA2)

job_id <- get_job_id(default = 1)

use_basilisk <- FALSE

# Both are FALSE for the runs of the paper. align_time = TRUE switches on time
# warping in MEFISTO and bool_speed_up = TRUE coarsens its optimisation; either
# changes the model that is fitted, and they were used only to test the code
# locally.
align_time <- FALSE
bool_speed_up <- FALSE

vec_N <- c(20, 25, 30, 35, 40)
n_seeds <- 26
n_job <- 2 * length(vec_N) * n_seeds
stopifnot(job_id >= 1, job_id <= n_job)

seed <- ((job_id - 1) %% n_seeds) + 1
N <- vec_N[((job_id - 1) %/% n_seeds) %% length(vec_N) + 1]

# TRUE: subject-specific trajectories (MEFISTO); FALSE: MEFISTO-base
mefisto_group_relationships <- (job_id <= n_job / 2)

set.seed(seed)

p <- 100

# all subjects share the same observation times when TRUE; FALSE in the paper
common_time_obs_for_all_subjects <- FALSE

if (common_time_obs_for_all_subjects) {
  N_subj <- 6
} else {
  N_min <- 2
  N_max <- 10
}

Q_true <- 2
Q <- 5
L_true <- 3
L <- 5

tol_abs <- 1e-3
tol_rel <- 1e-6
maxit <- 10000

mean_zeta <- 0
vec_mean_zeta <- rep(mean_zeta, L_true)
vec_sd_zeta <- 1/(1:L_true)
vec_sd_eps <- rep(1, p)

# for data generation: Beta(1, 1) probability of association, so that the
# loadings are not sparse (b_om = 10 in the other experiments)
a_om <- 1
b_om <- 1

# for inference: the default Beta(1, p) prior on the probability of association
c_0 <- 1
d_0 <- p

K <- NULL

n_repl <- 1

anneal <- c(1, 1.9, 100)

# grid on which the trajectories are evaluated and compared; coarser than the
# 1000 points of the other experiments, since MEFISTO interpolates on it
n_g <- 100
time_g <- seq(0, 1, length.out = n_g)

n_cpus <- 1 # serial runtimes, for comparability with MEFISTO
bool_save <- TRUE

if (bool_save) {
  res_dir <- paste0(out_dir, "/MEFISTO_corr_v2_2026_",
                    ifelse(mefisto_group_relationships, T, F),
                    ifelse(align_time, "align_time_T_TEST_ONLY_", ""),
                    "_speed_up_", ifelse(bool_speed_up, T, F),
                    "_mean_zeta_", mean_zeta,
                    "_anneal_", paste0(anneal, collapse = "-"),
                    "_N_", N, "_p_", p, "_Q_true_", Q_true,
                    "_Q_", Q, "_L_true_", L_true, "_L_", L,
                    ifelse(common_time_obs_for_all_subjects,
                           paste0("_N_subj_", N_subj),
                           paste0("_N_min_", N_min, "_N_max_", N_max)),
                    "_K_", K,
                    "_sd_zeta_", paste0(round(vec_sd_zeta, 2), collapse = "-"),
                    "_sd_eps_", paste0(round(unique(vec_sd_eps), 2), collapse = "-"),
                    "_a_om_", a_om, "_b_om_", b_om,
                    "_c_0_", c_0, "_d_0_", d_0, "_tol_abs_", tol_abs,
                    "_tol_rel_", tol_rel, "_maxit_", maxit,
                    "_n_repl_", n_repl, "_seed_", seed, "/")

  dir.create(res_dir, showWarnings = FALSE)

  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"), type = "message")
} else {
  res_dir <- NULL
}

# zero mean functions, since MEFISTO does not model variable-specific
# time-varying means
mu_func <- function(t, j = 1, p = 1) 0*t

list_res <- parallel::mclapply(1:n_repl, function(repl) {

  
  if (n_repl > 1) {
    repl_seed <- repl
  } else {
    repl_seed <- seed
  }
  set.seed(repl_seed) 
  
  print(paste0("=== SEED SET TO ", repl_seed, " ==="))
  
  if (common_time_obs_for_all_subjects) {
    
    n_i <- N_subj # number of obs per subject
    time_obs_i <- sort(runif(n_i, min = 0, max = 1))
    
    n <- rep(n_i, N)
    time_obs <- replicate(N, time_obs_i, simplify = FALSE)
    
    
    # variant of generate_bayesSYNC_data that also returns the noise-free
    # latent and observed trajectories on time_g, and allows a nonzero score
    # mean (zero here)
    fpca_data <- generate_bayesSYNC_data_mefisto(N = N, n = n, p = p, Q = Q_true, L = L_true,
                                                 mu_func = mu_func,
                                                 n_g = n_g,  
                                                 time_obs = time_obs,
                                                 a_om = a_om, b_om = b_om,
                                                 vec_mean_zeta = vec_mean_zeta, 
                                                 vec_sd_zeta = vec_sd_zeta,
                                                 vec_sd_eps = vec_sd_eps,
                                                 seed = repl_seed)
  } else {
    
    n <- round(runif(N, N_min, N_max))
    
    # variant of generate_bayesSYNC_data that also returns the noise-free
    # latent and observed trajectories on time_g, and allows a nonzero score
    # mean (zero here)
    fpca_data <- generate_bayesSYNC_data_mefisto(N = N, n = n, p = p, Q = Q_true, L = L_true,
                                                 mu_func = mu_func,
                                                 n_g = n_g,  
                                                 a_om = a_om, b_om = b_om,
                                                 vec_mean_zeta = vec_mean_zeta, 
                                                 vec_sd_zeta = vec_sd_zeta,
                                                 vec_sd_eps = vec_sd_eps,
                                                 seed = repl_seed)
  }


  time_obs <- fpca_data$time_obs
  Zeta <- fpca_data$Zeta
  B <- fpca_data$B
  pat <- fpca_data$gamma
  Y <- fpca_data$Y
  time_g <- fpca_data$time_g
  mu_g <- fpca_data$mu_g
  Phi_g <- fpca_data$Phi_g
  
  list_h_true <- fpca_data$list_h_true
  Y_true <- fpca_data$Y_true

  
  run_mefisto <- T
  if (run_mefisto) {
    data <- reshape_Y_to_df(Y, time_obs) # reformat Y for use by MOFA2
    
    indiv_as_groups <- T # yes, otherwise not displaying output for each subject
    if (indiv_as_groups) { # we will have Q factors *per* subject i, i.e., h_i^{(q)}(t) - same as in bayesSYNC
      names(data)[names(data) == "subject"] <- "group"
    }
    sm_untrained <- create_mofa(data = data)
    
    sm_untrained <- set_covariates(sm_untrained, covariates = "time")
    
    
    data_opts <- get_default_data_options(sm_untrained)
    data_opts$center_groups <- FALSE # In addition we do not center the data from individual groups as the data was already centered across all samples. This can be preferable if differences between groups should be maintained and the feature levels are comparable between groups.
    
    model_opts <- get_default_model_options(sm_untrained)
    model_opts$num_factors <- Q
    
    train_opts <- get_default_training_options(sm_untrained)
    train_opts$seed <- repl_seed
    
    mefisto_opts <- get_default_mefisto_options(sm_untrained)
    if (indiv_as_groups) {
      mefisto_opts$model_groups <- mefisto_group_relationships # if model_groups is set to T infer a group-group correlation matrix that indicates for each factor how the groups relate to one another. See also our evodevo tutorial for an example. Note that setting model_groups to True can be slow for large number of groups - in this case, we recommend setting it to False for initial analysis.
      
      if (align_time){
        cat("Time warping switched on: not the model fitted in the paper, for local tests only.\n")
        mefisto_opts$warping <- TRUE # whether to align the time or not
        mefisto_opts$warping_ref <- "subject_1"
        mefisto_opts$warping_freq <- 10
      } 
    }
    
    if (bool_speed_up){ # can speed up inference if needed, instead of using the default parameters
      train_opts$maxiter <- 100 # 1000 is the default
      mefisto_opts$n_grid <- 10
      mefisto_opts$start_opt <- 50 # must be multiple of warping_freq
      mefisto_opts$opt_freq <- 50
    }
    mefisto_opts$new_values <- matrix(time_g, nrow = 1) 
    
    sm_untrained <- prepare_mofa(sm_untrained, 
                                 model_options = model_opts,
                                 mefisto_options = mefisto_opts,
                                 training_options = train_opts,
                                 data_options = data_opts)
    
    outfile <- paste0(res_dir, "mefisto_model_repl_", repl_seed, ".hdf5")
    
    runtime_mefisto <- system.time(sm <- run_mofa(sm_untrained, 
                                                  outfile = outfile,
                                                  use_basilisk = use_basilisk))
    print(runtime_mefisto)
    
    sm <- load_model(outfile, load_interpol_Z = TRUE) #
    
    
    
    
    B_hat_mefisto <- sm@expectations$W$single_view # compare with B_hat
    nb_factors_retained_mefisto <- ncol(B_hat_mefisto) # keep here, before the next line that is going to trim further
    
    if (nb_factors_retained_mefisto < Q_true) {
      sm <- load_model(outfile, load_interpol_Z = TRUE, remove_inactive_factors = FALSE) # avoid removing the factors as fewer than Q_true
      B_hat_mefisto <- sm@expectations$W$single_view # compare with B_hat
    } 
    list_h_hat_mefisto <- lapply(sm@interpolated_Z, function(ll) t(ll$mean))
    list_h_hat_mefisto <- list_h_hat_mefisto[order(as.numeric(gsub("subject_", "", names(list_h_hat_mefisto))))]
    
    list_h_hat_variance_mefisto <- lapply(sm@interpolated_Z, function(ll) t(ll$variance))
    list_h_hat_variance_mefisto <- list_h_hat_variance_mefisto[order(as.numeric(gsub("subject_", "", names(list_h_hat_variance_mefisto))))]
  
    list_Y_hat_mefisto <- lapply(list_h_hat_mefisto, function(ll) ll %*% t(B_hat_mefisto))
  

    # needed to estimate low and upp
    list_h_hat_obs_mefisto <- sm@expectations$Z # only needed to estimate the error variance 
    list_h_hat_obs_mefisto <- list_h_hat_obs_mefisto[order(as.numeric(gsub("subject_", "", names(list_h_hat_obs_mefisto))))]
    
    perm_mefisto <- match_factor_and_sign_mefisto(B, B_hat_mefisto, list_h_hat_mefisto, list_h_hat_obs_mefisto, list_h_hat_variance_mefisto)
    B_hat_mefisto <- perm_mefisto$perm_B_hat
    list_h_hat_mefisto <- perm_mefisto$perm_list_h_hat
    list_h_hat_variance_mefisto <- perm_mefisto$perm_list_h_hat_variance
    list_h_hat_obs_mefisto <- perm_mefisto$perm_list_h_hat_obs
  
    
 
    list_Y_hat_obs_mefisto <- lapply(list_h_hat_obs_mefisto, function(ll) ll %*% t(B_hat_mefisto))
    
    # Compute the noise variance per feature for each group
    noise_var_from_residuals <- sapply(1:p, function(j) {
      var(unlist(sapply(1:N, function(i) Y[[i]][[j]] - list_Y_hat_obs_mefisto[[i]][,j])))
    })
    
    mat_sd_i_j <- sapply(1:p, function(j) sapply(1:N, function(i) sqrt(Reduce("+", lapply(1:Q_true, function(q) B_hat_mefisto[j, q]^2 * list_h_hat_variance_mefisto[[i]][[q]])) + noise_var_from_residuals[j])))# B assumed to be known exactly (variance = b^T Phi^T Cov_zeta Phi b) so sqrt(b^2) for sd
    
    list_Y_low_mefisto <- lapply(1:N, function(i) sapply(1:p, function(j) list_Y_hat_mefisto[[i]][,j] + qnorm(0.025) * mat_sd_i_j[i,j]))
    list_Y_upp_mefisto <- lapply(1:N, function(i) sapply(1:p, function(j) list_Y_hat_mefisto[[i]][,j] + qnorm(0.975) * mat_sd_i_j[i,j]))
    
    
    mat_coverage_ci_Y_mefisto <- sapply(1:p, function(j) sapply(1:N, function(i) {
      coverage_curve(Y[[i]][[j]], time_obs[[i]], list_Y_low_mefisto[[i]][,j], list_Y_upp_mefisto[[i]][,j], time_g)}))
    
    
    mat_mean_width_ci_Y_mefisto <- sapply(1:p, function(j) sapply(1:N, function(i) {
      credible_band_width(list_Y_low_mefisto[[i]][,j], list_Y_upp_mefisto[[i]][,j], time_g)$mean_width}))
    
    mat_median_width_ci_Y_mefisto <- sapply(1:p, function(j) sapply(1:N, function(i) {
      credible_band_width(list_Y_low_mefisto[[i]][,j], list_Y_upp_mefisto[[i]][,j], time_g)$median_width}))
    
    mat_integrated_width_ci_Y_mefisto <- sapply(1:p, function(j) sapply(1:N, function(i) {
      credible_band_width(list_Y_low_mefisto[[i]][,j], list_Y_upp_mefisto[[i]][,j], time_g)$integrated_width}))
    
    
    
  }

  list_hyper <- set_hyper(c_0 = c_0, d_0 = d_0)

  runtime_bayesSYNC <- system.time(fpca_res <- bayesSYNC(time_obs = time_obs,
                                                         Y = Y, L = L, Q = Q,
                                                         K = K,
                                                         n_g = n_g,
                                                         anneal = anneal,
                                                         list_hyper = list_hyper,
                                                         tol_abs = tol_abs,
                                                         tol_rel = tol_rel,
                                                         maxit = maxit,
                                                         seed = repl_seed,
                                                         n_cpus = 1, # serial, as MEFISTO
                                                         bool_scale = F)) # already roughly on the same scale

  print(runtime_bayesSYNC)
  

  ppi <- fpca_res$ppi
  factor_ppi <- fpca_res$factor_ppi
  
  factor_ppi_thres <- 0.5
  nb_factors_retained <- sum(factor_ppi > factor_ppi_thres)
  
  B_hat <- fpca_res$B_hat
  list_cumulated_pve <- fpca_res$list_cumulated_pve
  elbo <- fpca_res$ELBO_iter

  list_Y_hat <- fpca_res$list_Y_hat
  list_Y_low <- fpca_res$list_Y_low
  list_Y_upp <- fpca_res$list_Y_upp
  
  list_h_hat <- fpca_res$list_h_hat
  list_var_vec <- fpca_res$list_var_vec # to reconstruct  list_h_upp and  list_h_low after swapping the signs (using match_factor_and_sign below)
  
  list_mu_hat <- fpca_res$list_mu_hat
  list_list_Phi_hat <- fpca_res$list_list_Phi_hat
  list_Zeta_hat <- fpca_res$list_Zeta_hat
  list_Cov_zeta_hat <- fpca_res$list_Cov_zeta_hat
  


  perm_fpca_res <- match_factor_and_sign(B, B_hat, ppi, factor_ppi, Zeta,
                                         list_Zeta_hat, list_list_Phi_hat,
                                         list_cumulated_pve, list_Cov_zeta_hat,
                                         list_h_hat, list_var_vec)

  B_hat <- perm_fpca_res$perm_B_hat
  ppi <- perm_fpca_res$perm_ppi
  list_list_Phi_hat <- perm_fpca_res$perm_list_list_Phi_hat
  list_Zeta_hat <- perm_fpca_res$perm_list_Zeta_hat
  list_h_hat <- perm_fpca_res$perm_list_h_hat
  list_h_low <- perm_fpca_res$perm_list_h_low
  list_h_upp <- perm_fpca_res$perm_list_h_upp

  
  factor_ppi <- perm_fpca_res$perm_factor_ppi
  
  print(paste0("Factor PPI after permutation based on true loadings: ",
               paste0(format(factor_ppi, digits = 2), collapse = ", "))) # expect to see the largest values first
  
  
  if (Q!=Q_true) {
    list_Zeta_hat_untrimmed <- perm_fpca_res$perm_list_Zeta_hat_untrimmed
    list_list_Phi_hat_untrimmed <- perm_fpca_res$perm_list_list_Phi_hat_untrimmed
    B_hat_untrimmed <- perm_fpca_res$perm_B_hat_untrimmed
    ppi_untrimmed <- perm_fpca_res$perm_ppi_untrimmed
  } else {
    list_Zeta_hat_untrimmed <- list_list_Phi_hat_untrimmed <-  B_hat_untrimmed <- ppi_untrimmed <- NULL
  }
  
  
  list_cumulated_pve <- perm_fpca_res$perm_list_cumulated_pve
  list_Cov_zeta_hat <- perm_fpca_res$perm_list_Cov_zeta_hat
  
  list_list_zeta_ellipse <-  vector("list", length=Q_true)
  for (q in 1:Q_true) {
    list_list_zeta_ellipse[[q]] <-  vector("list", length=N)
    for(i in 1:N) {
      zeta_mean <- list_Zeta_hat[[q]][i,][1:2]
      
      zeta_ellipse <- ellipse::ellipse(
        list_Cov_zeta_hat[[q]][[i]][1:2, 1:2],
        centre = zeta_mean,
        level = 0.95
      )
      list_list_zeta_ellipse[[q]][[i]] <- zeta_ellipse
    }
  }
  
  mat_coverage_ci_Y <- sapply(1:p, function(j) sapply(1:N, function(i) {
    coverage_curve(Y[[i]][[j]], time_obs[[i]], list_Y_low[[i]][[j]], list_Y_upp[[i]][[j]], time_g)}))
  
  mat_mean_width_ci_Y <- sapply(1:p, function(j) sapply(1:N, function(i) {
    credible_band_width(list_Y_low[[i]][[j]], list_Y_upp[[i]][[j]], time_g)$mean_width}))
  
  mat_median_width_ci_Y <- sapply(1:p, function(j) sapply(1:N, function(i) {
    credible_band_width(list_Y_low[[i]][[j]], list_Y_upp[[i]][[j]], time_g)$median_width}))
  
  mat_integrated_width_ci_Y <- sapply(1:p, function(j) sapply(1:N, function(i) {
    credible_band_width(list_Y_low[[i]][[j]], list_Y_upp[[i]][[j]], time_g)$integrated_width}))

 
  vec_rmse_B_hat <- sapply(1:Q_true, function(q) sqrt(mean((B_hat[,q] - B[,q])^2)))
  vec_rmse_B_hat_mefisto <- sapply(1:Q_true, function(q) sqrt(mean((B_hat_mefisto[,q] - B[,q])^2)))
 
   
  # MEFISTO does not estimate mean functions, which is why the data are
  # generated with zero means and the latent trajectories can be compared directly
  mat_ise_list_h_hat <- sapply(1:Q_true, function(q) sapply(1:N, function(i) compute_ise(time_g, list_h_true[[i]][,q], list_h_hat[[i]][[q]])))
  mat_ise_list_h_hat_mefisto <- sapply(1:Q_true, function(q) sapply(1:N, function(i) compute_ise(time_g, list_h_true[[i]][,q], list_h_hat_mefisto[[i]][,q])))
  
  colnames(mat_ise_list_h_hat) <- colnames(mat_ise_list_h_hat_mefisto) <- paste0("factor_", 1:Q_true)  
  rownames(mat_ise_list_h_hat) <- rownames(mat_ise_list_h_hat_mefisto) <- paste0("subject_", 1:N)
  
  mat_ise_list_Y_hat <- sapply(1:p, function(j) sapply(1:N, function(i) compute_ise(time_g, Y_true[[i]][[j]], list_Y_hat[[i]][[j]])))
  mat_ise_list_Y_hat_mefisto <- sapply(1:p, function(j) sapply(1:N, function(i) compute_ise(time_g, Y_true[[i]][[j]], list_Y_hat_mefisto[[i]][,j])))
  
  colnames(mat_ise_list_Y_hat) <- colnames(mat_ise_list_Y_hat_mefisto) <- paste0("variable_", 1:p)  
  rownames(mat_ise_list_Y_hat) <- rownames(mat_ise_list_Y_hat_mefisto) <- paste0("subject_", 1:N)
  
  
  Corr_zeta <- corr_zeta_vector(Zeta, list_Zeta_hat)
  Corr_b <- corr_b_vector(B, B_hat)
  Corr_b_mefisto <- corr_b_vector(B, B_hat_mefisto)
  
  Auc_b <- auc_roc_curve(B, ppi)
  Auc_b_mefisto <- auc_roc_curve(B, abs(B_hat_mefisto))
  
  Err_mean <- error_mean_function(time_g, list_mu_hat, mu_g) / p # mean error of the mean function across all variables.
  Err_eigen <- error_eigen_function(time_g, Phi_g, C_g = NULL, list_list_Phi_hat)
  
  rt <- runtime_bayesSYNC[["elapsed"]] / 60
  rt_mefisto <- runtime_mefisto[["elapsed"]] / 60
  
  create_named_list(mat_ise_list_Y_hat, mat_ise_list_Y_hat_mefisto,
                    mat_ise_list_h_hat, mat_ise_list_h_hat_mefisto,
                    mat_coverage_ci_Y,
                    mat_mean_width_ci_Y,
                    mat_median_width_ci_Y,
                    mat_integrated_width_ci_Y,
                    mat_coverage_ci_Y_mefisto,
                    mat_mean_width_ci_Y_mefisto,
                    mat_median_width_ci_Y_mefisto,
                    mat_integrated_width_ci_Y_mefisto,
                    vec_rmse_B_hat, vec_rmse_B_hat_mefisto,
                    Corr_zeta, 
                    Corr_b, Corr_b_mefisto, 
                    Auc_b, Auc_b_mefisto, Err_mean, Err_eigen, 
                    rt, elbo,
                    time_obs, pat, B, Zeta, Y, mu_g, Phi_g,
                    ppi, ppi_untrimmed, factor_ppi, B_hat, B_hat_untrimmed,
                    list_mu_hat, list_list_Phi_hat, list_Zeta_hat, list_list_zeta_ellipse,
                    list_h_hat, list_h_low, list_h_upp,
                    list_Y_hat, list_Y_low, list_Y_upp,
                    list_h_true, list_h_hat_variance_mefisto,
                    list_cumulated_pve, list_Zeta_hat_untrimmed, list_list_Phi_hat_untrimmed,
                    rt_mefisto, B_hat_mefisto, list_h_hat_mefisto, list_Y_hat_mefisto,
                    nb_factors_retained, nb_factors_retained_mefisto)
  
}, mc.cores = n_cpus)

if (bool_save) {
  save(list_res, file = file.path(res_dir, "output.RData"))
}
