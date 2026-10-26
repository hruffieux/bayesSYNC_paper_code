## ---------------------------------------------------------------------------
## Gathers the performance metrics of the 50 replicates of each configuration
## of run_large_scale.R into one file per configuration, which
## plot_figure3_annealing.R then reads. The per-replicate output.RData files are
## large, so only the metrics are kept.
##
## Usage, from the root of the repository:
##   Rscript simulations/large_scale/gather_large_scale.R
##
## Writes, for each configuration, a folder BS_r_summary_... under out_dir
## holding output.RData. Stops if the output of any replicate is missing, since
## the summaries are computed over all replicates.
## ---------------------------------------------------------------------------

rm(list = ls())

source("R/paths.R")
source("R/utils_sim.R")

p <- 20000
N <- 100
N_min <- 5
N_max <- 10

Q_true <- 3
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

bool_var_spec_prob <- FALSE

n_cpus <- 50 # only enters the folder names written by run_large_scale.R

n_repl <- 1
sum_n_repl <- 50

list_anneal <- list(NULL, c(1, 1.9, 100))

for (Q in c(3, 5)) {

  for (anneal in list_anneal) {

    sum_res_dir <- paste0(out_dir, "/BS_r_summary_anneal_", paste0(anneal, collapse = "-"),
                          "_summary_mu_diff_N_", N, "_p_", p, "_Q_true_", Q_true,
                          "_Q_", Q, "_L_true_", L_true, "_L_", L,
                          "_N_min_", N_min, "_N_max_", N_max,
                          "_K_", K,
                          "_sd_zeta_", paste0(round(vec_sd_zeta, 2), collapse = "-"),
                          "_sd_eps_", paste0(round(unique(vec_sd_eps), 2), collapse = "-"),
                          ifelse(bool_var_spec_prob, "_omega_j_q", ""),
                          "_a_om_", a_om, "_b_om_", b_om, "_tol_abs_", tol_abs,
                          "_tol_rel_", tol_rel, "_maxit_", maxit, "_n_cpus_", max(floor(n_cpus/n_repl), 1),
                          "_n_repl_", sum_n_repl, "/")
    dir.create(sum_res_dir, showWarnings = FALSE)

    list_Corr_zeta <- list_Corr_b <- list_Auc_b <- list_Err_mean <- list_Err_eigen <- list_rt <- list_pat <- list_ppi <- list_factor_ppi <- list_elbo <- NULL

    for (seed in 1:sum_n_repl) {

      res_seed_dir <- paste0(out_dir, "/BS_r_anneal_", paste0(anneal, collapse = "-"),
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

      out_seed_file <- paste0(res_seed_dir, "output.RData")

      if (file.exists(out_seed_file)) {
        load(out_seed_file)
      } else {
        stop(paste0("No output for seed = ", seed, ", Q = ", Q,
                    ", anneal = ", paste0(anneal, collapse = "-")))
      }

      list_Corr_zeta <- c(list_Corr_zeta, lapply(list_res, "[[", "Corr_zeta"))
      list_Corr_b <- c(list_Corr_b, lapply(list_res, "[[", "Corr_b"))
      list_Auc_b <- c(list_Auc_b, lapply(list_res, "[[", "Auc_b"))
      list_Err_mean <- c(list_Err_mean, lapply(list_res, "[[", "Err_mean"))
      list_Err_eigen <- c(list_Err_eigen, lapply(list_res, "[[", "Err_eigen"))

      list_rt <- c(list_rt, lapply(list_res, "[[", "rt"))

      list_pat <- c(list_pat, lapply(list_res, "[[", "pat"))
      list_ppi <- c(list_ppi, lapply(list_res, "[[", "ppi"))
      list_factor_ppi <- c(list_factor_ppi, lapply(list_res, "[[", "factor_ppi"))

      list_elbo <- c(list_elbo, lapply(list_res, "[[", "elbo"))

      rm(list_res)

    }

    save(list_Corr_zeta, list_Corr_b, list_Auc_b, list_Err_mean, list_Err_eigen,
         list_rt, list_pat, list_ppi, list_factor_ppi, list_elbo,
         file = file.path(sum_res_dir, "output.RData"))

  }
}
