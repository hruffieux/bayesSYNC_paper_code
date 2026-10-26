## ---------------------------------------------------------------------------
## Panels of Figure 3: effect of annealing on estimation accuracy and factor
## selection, over the 50 replicates of the large-scale simulation.
##
## For vanilla VB (panel a) and annealed VB (panel b), compares Q_max = 3
## (oracle) with Q_max = 5 on the posterior factor inclusion probabilities,
## the AUC for the variable inclusion of the loadings, the ISE of the mean
## functions and the ISE of the eigenfunctions. The score correlations quoted
## in the caption are printed, as means (standard errors) over replicates.
## Error bars are one standard error over replicates.
##
## Reads the files written by gather_large_scale.R.
##
## Usage, from the root of the repository:
##   Rscript simulations/large_scale/plot_figure3_annealing.R
##
## Writes, for each of the two inference schemes, a folder BS_r_subset_summary_...
## under out_dir holding the panels as separate pdf files.
## ---------------------------------------------------------------------------

rm(list = ls())

source("R/paths.R")
source("R/plot_utils.R")
source("R/utils_sim.R") # also loads pROC

require(matrixStats)

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

n_repl <- 1
sum_n_repl <- 50

bool_var_spec_prob <- FALSE

n_cpus <- 50 # only enters the folder names written by run_large_scale.R
bool_save <- TRUE

list_param <- list_Q <- list("Q = 3" = 3, "Q = 5" = 5)

vec_col <- c("skyblue1", "honeydew3")

for (bool_anneal in c(FALSE, TRUE)) {

  if (bool_anneal) {
    anneal <- c(1, 1.9, 100)
    alpha <- 1
  } else {
    anneal <- NULL
    alpha <- 0.5
  }

  part_res_dir <- paste0("anneal_", paste0(anneal, collapse = "-"), "_list_Q_", paste0(list_Q, collapse = "-"))

  res_list <- lapply(seq_along(list_param), function(id) {

    Q <- list_param[[id]]

    r_file <- paste0(out_dir, "/BS_r_summary_anneal_", paste0(anneal, collapse = "-"),
                     "_summary_mu_diff_N_", N, "_p_", p, "_Q_true_", Q_true,
                     "_Q_", Q, "_L_true_", L_true, "_L_", L,
                     "_N_min_", N_min, "_N_max_", N_max,
                     "_K_", K,
                     "_sd_zeta_", paste0(round(vec_sd_zeta, 2), collapse = "-"),
                     "_sd_eps_", paste0(round(unique(vec_sd_eps), 2), collapse = "-"),
                     ifelse(bool_var_spec_prob, "_omega_j_q", ""),
                     "_a_om_", a_om, "_b_om_", b_om, "_tol_abs_", tol_abs,
                     "_tol_rel_", tol_rel, "_maxit_", maxit, "_n_cpus_", max(floor(n_cpus/n_repl), 1),
                     "_n_repl_", sum_n_repl, "/output.RData")
    load(r_file)

    create_named_list(list_Corr_zeta, list_Corr_b, list_Auc_b, list_Err_mean, list_Err_eigen,
                      list_rt, list_pat, list_ppi, list_factor_ppi)
  })

  names(res_list) <- names(list_param)

  if (bool_save) {
    sum_res_dir <- paste0(out_dir, "/BS_r_subset_summary_", part_res_dir,
                          "_N_", N, "_p_", p, "_Q_true_", Q_true,
                          "_L_true_", L_true, "_L_", L,
                          "_N_min_", N_min, "_N_max_", N_max,
                          "_K_", K,
                          "_sd_zeta_", paste0(round(vec_sd_zeta, 2), collapse = "-"),
                          "_sd_eps_", paste0(round(unique(vec_sd_eps), 2), collapse = "-"),
                          ifelse(bool_var_spec_prob, "_omega_j_q", ""),
                          "_a_om_", a_om, "_b_om_", b_om, "_tol_abs_", tol_abs,
                          "_tol_rel_", tol_rel, "_maxit_", maxit, "_n_cpus_", max(floor(n_cpus/n_repl), 1),
                          "_n_repl_", sum_n_repl, "/")
    dir.create(sum_res_dir, showWarnings = FALSE)
  } else {
    sum_res_dir <- NULL
  }

  ## Variable inclusion: AUC over all factors, from the loading PPIs
  if (bool_save) {
    pdf(paste0(sum_res_dir, "/boxplots_AUC_variable_contributions_all_factors.pdf"), width = 3.5, height = 6)
  }

  mat_auc <- NULL
  for (ii in seq_along(list_param)) {
    list_roc_obj <- lapply(seq_along(res_list[[ii]]$list_pat), function(i) roc(as.vector(res_list[[ii]]$list_pat[[i]]), as.vector(res_list[[ii]]$list_ppi[[i]])))
    mat_auc <- cbind(mat_auc, sapply(list_roc_obj, function(ll) auc(ll)))
  }
  colnames(mat_auc) <- names(list_param)

  par(mfrow = c(1, 1))
  plot_bar_with_error(mat_auc, ylab = "AUC", main_title = "Variable inclusion", bar_colors = adjustcolor(vec_col, alpha.f = alpha))

  if (bool_save) {
    dev.off()
  }

  ## Score correlations, by component
  if (bool_save) {
    pdf(paste0(sum_res_dir, "/boxplots_cor_zeta.pdf"), width = 7, height = 6)
  }

  par(mfrow = c(ceiling(L_true / 2), 2))

  for (l in 1:L_true) {
    mat_zeta_l <- NULL

    for (ii in seq_along(list_param)) {
      mat_zeta_l <- cbind(mat_zeta_l, unlist(lapply(res_list[[ii]]$list_Corr_zeta, function(ll) ll[, l])))
    }

    colnames(mat_zeta_l) <- names(list_param)

    plot_bar_with_error(
      mat_zeta_l,
      ylab = "Correlation",
      main_title = paste0("Scores, l = ", l),
      bar_colors = adjustcolor(vec_col, alpha.f = alpha),
      ylim_max = 1.2*max(unlist(lapply(seq_along(list_param),
                                       function(id) mean(unlist(res_list[[id]]$list_Corr_zeta), na.rm = T) +
                                         sd(unlist(res_list[[id]]$list_Corr_zeta), na.rm = T)/sqrt(sum(!is.na(unlist(res_list[[id]]$list_Corr_zeta)))))))

    )
  }

  ## ISE of the mean functions (averaged over the p variables) and of the
  ## eigenfunctions
  if (bool_save) {
    dev.off()
    pdf(paste0(sum_res_dir, "/boxplots_ISE_functions.pdf"), width = 7.2, height = 3.65)
  }

  par(mfrow = c(ceiling(L_true / 2), 3))

  mat_ise_mean <- NULL
  for (ii in seq_along(list_param)) {
    mat_ise_mean <- cbind(mat_ise_mean, unlist(lapply(res_list[[ii]]$list_Err_mean, function(ll) ll / p)))
  }

  colnames(mat_ise_mean) <- names(list_param)

  ylim_max_ise <- 1.2*max(c(unlist(lapply(seq_along(list_param), function(id) mean(unlist(res_list[[id]]$list_Err_eigen)) + sd(unlist(res_list[[id]]$list_Err_eigen))/sqrt(sum_n_repl)))))

  plot_bar_with_error(
    mat_ise_mean,
    ylab = "Integrated squared error",
    main_title = paste0("Mean function"),
    bar_colors = adjustcolor(vec_col, alpha.f = alpha),
    error_bar_length = 0.15,
    ylim_max = 0.9
  )

  for (l in 1:L_true) {
    mat_ise_eigen_l <- NULL

    for (ii in seq_along(list_param)) {
      mat_ise_eigen_l <- cbind(mat_ise_eigen_l, unlist(lapply(res_list[[ii]]$list_Err_eigen, function(ll) ll[, l])))
    }

    colnames(mat_ise_eigen_l) <- names(list_param)

    plot_bar_with_error(
      mat_ise_eigen_l,
      ylab = "Integrated squared error",
      main_title = paste0("Eigenfunction, l = ", l, ""),
      bar_colors = adjustcolor(vec_col, alpha.f = alpha),
      error_bar_length = 0.15,
      ylim_max = max(ylim_max_ise, 0.9)  )
  }

  if (bool_save) {
    dev.off()
  }

  ## Posterior factor inclusion probabilities, one file per value of Q_max
  par(mfrow = c(1,1))

  for (id in 1:2) {

    mat_factor_ppi <- do.call(rbind, res_list[[id]]$list_factor_ppi)
    colnames(mat_factor_ppi) <- paste0("factor ", 1:ncol(mat_factor_ppi))
    avg_factor_ppi <- colMeans(mat_factor_ppi)
    se_factor_ppi <- colSds(mat_factor_ppi) / sqrt(nrow(mat_factor_ppi))
    se_factor_ppi[se_factor_ppi == 0] <- 0.001 # so that a zero standard error still shows as a bar

    if (bool_save) {
      pdf(paste0(sum_res_dir, "/average_factor_ppi_Q_", list_param[[id]], ".pdf"), width = 4.7, height = 5.4)
    }
    ze_barplot <- barplot(
      avg_factor_ppi,
      beside = TRUE,
      col = adjustcolor(vec_col[id], alpha.f = alpha),
      ylim = c(0, 1.1),
      space = 0.7,
      names.arg = colnames(mat_factor_ppi),
      ylab = "Posterior probability of inclusion",
      main = "Factor relevance"
    )

    error.bar <- function(x, y, upper, lower = upper, length = 0.1, lwd = 2, ...) {
      upper <- pmin(upper, 1)
      arrows(x, y + upper, x, y - lower, angle = 90, code = 3, length = length, lwd = lwd, ...)
    }

    error.bar(ze_barplot,avg_factor_ppi,se_factor_ppi, length = 0.1, lwd = 2)
    if (bool_save) {
      dev.off()
    }

  }

  ## Score correlations quoted in the caption: mean (standard error) over
  ## replicates, for Q_max = 3 and Q_max = 5
  cat("\nScore correlations, anneal = ", paste0(anneal, collapse = "-"), "\n", sep = "")
  for (l in 1:L_true) {
    mat_zeta_l <- NULL

    for (ii in seq_along(list_param)) {
      mat_zeta_l <- cbind(mat_zeta_l, unlist(lapply(res_list[[ii]]$list_Corr_zeta, function(ll) ll[, l])))
    }

    means <- colMeans(mat_zeta_l, na.rm = T)
    se <- apply(mat_zeta_l, 2, function(x) sd(x, na.rm = T) / sqrt(sum(!is.na(x))))
    print(paste0("l = ", l, ", ", names(list_param), ": ",
                 format(means, digits = 3), " (", format(se, digits = 2), ")"))

  }
}
