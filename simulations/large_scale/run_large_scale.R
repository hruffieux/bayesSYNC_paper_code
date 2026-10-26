## ---------------------------------------------------------------------------
## Large-scale simulation, Sections 5.2 and 5.3 of the paper (Figures 2 to 4).
##
## N = 100 subjects, p = 20,000 variables, Q = 3 latent factors with L = 2
## FPCA components each, 5 to 10 observations per subject, and a factor-specific
## probability of association drawn from a Beta(1, 10) distribution.
##
## One array task fits one (configuration, replicate) pair. The four
## configurations cross the number of candidate factors, Q_max = 3 (oracle) or 5,
## with vanilla or annealed variational inference:
##
##   job_id   1 to  50   Q_max = 3, vanilla VB
##   job_id  51 to 100   Q_max = 3, annealed VB
##   job_id 101 to 150   Q_max = 5, vanilla VB
##   job_id 151 to 200   Q_max = 5, annealed VB
##
## The replicate index is the seed of the data generation, so the four
## configurations of a replicate are fitted to the same simulated data set.
##
## Meant for a cluster, one array task per job (slurm/submit_large_scale.sh):
## each job was given 50 cores, 256 Gb of memory and up to 36 hours, and its
## output takes several tens of Gb on disk. It is not meant to be run on a
## laptop.
##
## Usage, from the root of the repository:
##   Rscript simulations/large_scale/run_large_scale.R <job_id>
##
## Writes, in its own folder under out_dir: output.RData (the fit, the
## simulated data and the metrics) and out.txt (the console output, from which
## plot_figure4_elbo.R reads the final ELBO and the number of iterations). At
## p = 20,000 each output.RData is several tens of Gb.
## ---------------------------------------------------------------------------

rm(list = ls())

source("R/paths.R")
source("R/plot_utils.R")
source("R/utils_sim.R")

require(bayesSYNC)

job_id <- get_job_id(default = 151)

n_repl_total <- 50 # replicates per configuration

df_config <- expand.grid(bool_anneal = c(FALSE, TRUE), Q = c(3, 5))
stopifnot(job_id >= 1, job_id <= n_repl_total * nrow(df_config))

config_id <- ((job_id - 1) %/% n_repl_total) + 1
seed <- ((job_id - 1) %% n_repl_total) + 1
set.seed(seed)

p <- 20000
N <- 100
N_min <- 5
N_max <- 10

Q_true <- 3
Q <- df_config$Q[config_id]

L_true <- 2
L <- 5

tol_abs <- 1e-3
tol_rel <- 1e-6
maxit <- 10000

vec_sd_zeta <- 1/(1:L_true)
vec_sd_eps <- rep(1, p)

a_om <- 1
b_om <- 10

K <- NULL

n_repl <- 1 # one replicate per array task

if (df_config$bool_anneal[config_id]) {
  anneal <- c(1, 1.9, 100)
} else {
  anneal <- NULL
}

bool_var_spec_prob <- FALSE

n_cpus <- 50 # used inside bayesSYNC; must match --cpus-per-task in slurm/submit_large_scale.sh
bool_save <- TRUE

if (bool_save) {
  res_dir <- paste0(out_dir, "/BS_r_anneal_", paste0(anneal, collapse = "-"),
                    "_N_", N, "_p_", p, "_Q_true_", Q_true,
                    "_Q_", Q, "_L_true_", L_true, "_L_", L,
                    "_N_min_", N_min, "_N_max_", N_max,
                    "_K_", K,
                    "_sd_zeta_", paste0(round(vec_sd_zeta, 2), collapse = "-"),
                    "_sd_eps_", paste0(round(unique(vec_sd_eps), 2), collapse = "-"),
                    ifelse(bool_var_spec_prob, "_omega_j_q", ""),
                    "_a_om_", a_om, "_b_om_", b_om, "_tol_abs_", tol_abs,
                    "_tol_rel_", tol_rel, "_maxit_", maxit, "_n_cpus_", max(floor(n_cpus/n_repl), 1),
                    "_n_repl_", n_repl, "_seed_", seed, "/")
  dir.create(res_dir, showWarnings = FALSE)
  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"), type = "message")
} else {
  res_dir <- NULL
}

# a different mean function for each variable, so that factors in excess of
# Q_true cannot be used to compensate for a mean structure shared across variables
mu_func <- function(t, p, j = 1, alpha = 3) (-1)^j*alpha*sin((2*pi+j %% 5)*t + j / p)


list_res <- parallel::mclapply(1:n_repl, function(repl) {

  n <- round(runif(N, N_min, N_max))

  fpca_data <- generate_bayesSYNC_data(N = N, n = n, p = p, Q = Q_true, L = L_true,
                                       mu_func = mu_func,
                                       a_om = a_om, b_om = b_om,
                                       vec_sd_zeta = vec_sd_zeta,
                                       vec_sd_eps = vec_sd_eps,
                                       seed = seed)

  time_obs <- fpca_data$time_obs
  Zeta <- fpca_data$Zeta
  B <- fpca_data$B
  pat <- fpca_data$gamma
  Y <- fpca_data$Y
  time_g <- fpca_data$time_g
  mu_g <- fpca_data$mu_g
  Phi_g <- fpca_data$Phi_g

  runtime_bayesSYNC <- system.time(fpca_res <- bayesSYNC(time_obs = time_obs,
                                                         Y = Y, L = L, Q = Q,
                                                         K = K, anneal = anneal,
                                                         tol_abs = tol_abs,
                                                         tol_rel = tol_rel,
                                                         maxit = maxit,
                                                         seed = seed,
                                                         bool_var_spec_prob = bool_var_spec_prob,
                                                         n_cpus = max(floor(n_cpus/n_repl), 1),
                                                         bool_scale = F)) # already roughly on the same scale

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

  # permutation and sign matching to the simulated factors, used only to
  # compute the performance metrics (Section 5.1)
  perm_fpca_res <- match_factor_and_sign(B, B_hat, ppi, factor_ppi, Zeta,
                                         list_Zeta_hat, list_list_Phi_hat,
                                         list_cumulated_pve, list_Cov_zeta_hat)

  B_hat <- perm_fpca_res$perm_B_hat
  ppi <- perm_fpca_res$perm_ppi
  factor_ppi <- perm_fpca_res$perm_factor_ppi

  print(paste0("Factor PPI after permutation based on true loadings: ",
               paste0(format(factor_ppi, digits = 2), collapse = ", "))) # expect to see the largest values first

  list_list_Phi_hat <- perm_fpca_res$perm_list_list_Phi_hat
  list_Zeta_hat <- perm_fpca_res$perm_list_Zeta_hat

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

  # 95% credible ellipses and posterior standard deviations of the scores,
  # used by the score panels of Figure 2
  list_list_zeta_ellipse <-  list_sd_hat <- vector("list", length=Q_true)
  for (q in 1:Q_true) {
    list_list_zeta_ellipse[[q]] <- vector("list", length=N)
    sd_mat <- matrix(NA, nrow = N, ncol = L)
    for(i in 1:N) {
      zeta_mean <- list_Zeta_hat[[q]][i,][1:2]

      zeta_ellipse <- ellipse::ellipse(
        list_Cov_zeta_hat[[q]][[i]][1:2, 1:2],
        centre = zeta_mean,
        level = 0.95
      )
      list_list_zeta_ellipse[[q]][[i]] <- zeta_ellipse
      sd_mat[i, ] <- sqrt(diag(list_Cov_zeta_hat[[q]][[i]]))
    }
    colnames(sd_mat) <- paste0("FPC_", 1:L)
    rownames(sd_mat) <- paste0("subject_", 1:N)

    list_sd_hat[[q]] <- sd_mat

  }

  Corr_zeta <- corr_zeta_vector(Zeta, list_Zeta_hat)
  Corr_b <- corr_b_vector(B, B_hat)

  Auc_b <- auc_roc_curve(B, ppi)

  Err_mean <- error_mean_function(time_g, list_mu_hat, mu_g)
  Err_eigen <- error_eigen_function(time_g, Phi_g, C_g = NULL, list_list_Phi_hat)

  rt <- runtime_bayesSYNC[["elapsed"]] / 60

  create_named_list(Corr_zeta, Corr_b, Auc_b, Err_mean, Err_eigen, rt, elbo,
                    time_obs, pat, B, Zeta, Y, mu_g, Phi_g,
                    ppi, ppi_untrimmed, factor_ppi, B_hat, B_hat_untrimmed,
                    list_mu_hat, list_list_Phi_hat, list_Zeta_hat, list_list_zeta_ellipse,
                    list_sd_hat,
                    list_Y_hat, list_Y_low, list_Y_upp,
                    list_cumulated_pve, list_Zeta_hat_untrimmed, list_list_Phi_hat_untrimmed)

}, mc.cores = n_repl)

if (bool_save) {
  save(list_res, file = file.path(res_dir, "output.RData"))
}
