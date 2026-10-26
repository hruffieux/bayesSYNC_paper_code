## ---------------------------------------------------------------------------
## Summaries of the MEFISTO comparison for one value of N and one MEFISTO
## setting, over the replicates listed in vec_seed below. Read by
## plot_figure5_table1.R, which combines the values of N and the two settings.
##
## The replicates are seeds 1 to 25 except where a MEFISTO run failed, as
## follows. MEFISTO, N = 25: seed 7 failed (MEFISTO removed too many factors)
## and seed 26 replaced it. MEFISTO, N = 40: 14 of the 25 runs completed within
## the 36-hour limit. MEFISTO-base, N = 20: seeds 15 and 16 failed in MOFA2.
## MEFISTO-base, N = 25: seed 7 failed, as for MEFISTO, and was not replaced.
## MEFISTO-base, N = 35: seeds 6 and 7 failed. The bayesSYNC fits are those of
## the same replicates.
##
## Also writes, for every replicate, the estimated and simulated latent
## trajectories of each subject (Supplementary Material S4 shows those of
## seed 15 at N = 30), the ROC curves of Figure 5(a) (at N = 30) and paired
## boxplots of every metric.
##
## Usage, from the root of the repository:
##   Rscript simulations/mefisto_comparison/summarise_mefisto_comparison.R <N> <setting>
## with <setting> either "mefisto" or "mefisto_base".
##
## Writes runtime_summary.RData and metrics_summary.RData, with the figures, in
## a folder summarise_MEFISTO_corr_v2_2026_... under out_dir.
## ---------------------------------------------------------------------------

rm(list = ls())

source("R/paths.R")
source("R/plot_utils.R")
source("R/utils_sim.R")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2, args[2] %in% c("mefisto", "mefisto_base"))

set.seed(1)

p <- 100
N <- as.numeric(args[1])
stopifnot(N %in% c(20, 25, 30, 35, 40))

bool_speed_up <- FALSE

common_time_obs_for_all_subjects <- FALSE

if (common_time_obs_for_all_subjects) {
  N_subj <- 6
} else {
  N_min <- 2
  N_max <- 10
}

# TRUE: subject-specific trajectories (MEFISTO); FALSE: MEFISTO-base
mefisto_group_relationships <- (args[2] == "mefisto")

if (mefisto_group_relationships) {
  if (N %in% c(20, 30, 35)) {
    vec_seed <- 1:25
  } else if (N == 25) {
    vec_seed <- c(1:6, 8:26) # seed 7 failed, MEFISTO removing too many factors
  } else if (N == 40) {
    vec_seed <- c(1, 4, 7, 9, 12, 13, 15, 16, 17, 18, 20, 21, 22, 25) # runs completed within 36 hours
  }
  col_mefisto <- "#ff7f6e"
} else {
  if (N == 20) {
    vec_seed <- c(1:14, 17:25) # MOFA2 failed for seeds 15 and 16
  } else if (N == 25) {
    vec_seed <- c(1:6, 8:25) # seed 7 failed, MEFISTO removing too many factors
  } else if (N == 35) {
    vec_seed <- c(1:5, 8:25)
  } else if (N %in% c(30, 40)) {
    vec_seed <- 1:25
  }
  col_mefisto <- adjustcolor("#ff7f6e", alpha.f = 0.35)
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

a_om <- 1
b_om <- 1

c_0 <- 1
d_0 <- p

K <- NULL

n_repl <- 1
col <- "skyblue1"

anneal <- c(1, 1.9, 100)

n_g <- 100
time_g <- seq(0, 1, length.out = n_g)

bool_save <- TRUE

if (bool_save) {
  res_dir <- paste0(out_dir, "/summarise_MEFISTO_corr_v2_2026_",
                    ifelse(mefisto_group_relationships, T, F),
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
                    "_n_repl_", n_repl, "_nb_seeds_", length(vec_seed), "/")

  dir.create(res_dir, showWarnings = FALSE)

  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
} else {
  res_dir <- NULL
}


list_ppi <- list_pat <- list_B_hat_mefisto <- vec_nb_factors_retained <- vec_nb_factors_retained_mefisto <- NULL

## --- Trajectory-level errors (ISE) ---
list_mat_ise_list_Y_hat         <- NULL
list_mat_ise_list_Y_hat_mefisto <- NULL

## --- Latent factor curve errors (ISE) ---
list_mat_ise_list_h_hat         <- NULL
list_mat_ise_list_h_hat_mefisto <- NULL

## --- CI coverage and widths for Y ---
list_mat_coverage_ci_Y         <- NULL
list_mat_mean_width_ci_Y       <- NULL
list_mat_median_width_ci_Y     <- NULL
list_mat_integrated_width_ci_Y <- NULL

list_mat_coverage_ci_Y_mefisto         <- NULL
list_mat_mean_width_ci_Y_mefisto       <- NULL
list_mat_median_width_ci_Y_mefisto     <- NULL
list_mat_integrated_width_ci_Y_mefisto <- NULL

## --- Loadings recovery ---
list_vec_rmse_B_hat         <- NULL
list_vec_rmse_B_hat_mefisto <- NULL

## --- Runtime ---
list_rt         <- NULL
list_rt_mefisto <- NULL

list_Corr_b <- list_Corr_b_mefisto <- list_Auc_b <- list_Auc_b_mefisto <- NULL

list_factor_ppi <- vec_nb_factors_retained <- vec_nb_factors_retained_mefisto <- NULL

list_mat_coverage_h <- list_mat_coverage_h_mefisto <- list_mat_width_h <- list_mat_width_h_mefisto <- NULL

for (seed in vec_seed) {
  
  print("SEED")
  print(seed)
  
  repl_file <- paste0(out_dir, "/MEFISTO_corr_v2_2026_", 
                    ifelse(mefisto_group_relationships, T, F), 
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
                    "_tol_rel_", tol_rel, "_maxit_", maxit, # "_n_cpus_", max(floor(n_cpus/n_repl), 1),
                    "_n_repl_", n_repl, "_seed_", seed, "/output.RData")
  
  if (file.exists(repl_file)) {
    
    load(repl_file)
    
  } else {
    stop("Not all output files exist.")
  }
  
  list_pat <- c(list_pat, lapply(list_res, "[[", "pat"))
  list_ppi <- c(list_ppi, lapply(list_res, "[[", "ppi"))
  list_B_hat_mefisto <- c(list_B_hat_mefisto, lapply(list_res, "[[", "B_hat_mefisto"))
  
  list_h_hat              <- lapply(list_res, "[[", "list_h_hat")[[1]]
  list_h_low              <- lapply(list_res, "[[", "list_h_low")[[1]]
  list_h_upp              <- lapply(list_res, "[[", "list_h_upp")[[1]]
  list_h_true             <- lapply(list_res, "[[", "list_h_true")[[1]]
  
  list_h_hat_mefisto      <- lapply(list_res, "[[", "list_h_hat_mefisto")[[1]]
  list_h_hat_variance_mefisto <- lapply(list_res, "[[", "list_h_hat_variance_mefisto")[[1]]
  
  for (q in 1:Q_true) {
    
    if (bool_save) {
      pdf(paste0(res_dir, "/interpol_factors_facet_subject_bayesSYNC_factor_", q, "_seed_", seed, ".pdf"), width = 0.75*N, height = N/2)
    }
    print(plot_interpolated_factors_bayesSYNC(list_h_hat, list_h_low, list_h_upp, time_g, factors = q, list_h_true = list_h_true))
    if (bool_save) {
      dev.off()
    }
    
  }
  
  for (q in 1:Q_true) {
    
    if (bool_save) {
      pdf(paste0(res_dir, "/interpol_factors_facet_subject_mefisto_factor_", q, "_seed_", seed, ".pdf"), width = 0.75*N, height = N/2)
    }
    print(plot_interpolated_factors_mefisto(list_h_hat_mefisto, list_h_hat_variance_mefisto, time_g, factors = q, list_h_true = list_h_true))
    if (bool_save) {
      dev.off()
    }
    
  }
  
  
  list_mat_coverage_h <- c(list_mat_coverage_h, 
                           list(coverage_latent(list_h_true, list_h_low, list_h_upp)))

  
  list_h_low_mefisto <- lapply(seq_along(list_h_hat_mefisto), function(i) {
    lapply(1:ncol(list_h_hat_mefisto[[i]]), function(q) {
      mean_q <- list_h_hat_mefisto[[i]][,q]
      sd_q   <- sqrt(list_h_hat_variance_mefisto[[i]][,q])
      mean_q - 1.96 * sd_q
    })
  })
  
  list_h_upp_mefisto <- lapply(seq_along(list_h_hat_mefisto), function(i) {
    lapply(1:ncol(list_h_hat_mefisto[[i]]), function(q) {
      mean_q <- list_h_hat_mefisto[[i]][,q]
      sd_q   <- sqrt(list_h_hat_variance_mefisto[[i]][,q])
      mean_q + 1.96 * sd_q
    })
  })
  
  
  list_mat_coverage_h_mefisto <- c(list_mat_coverage_h_mefisto, 
                           list(coverage_latent(list_h_true, list_h_low_mefisto, list_h_upp_mefisto)))
  
  
  list_mat_width_h <- c(
    list_mat_width_h,
    list(width_latent(list_h_low, list_h_upp))
  )
  
  list_mat_width_h_mefisto <- c(
    list_mat_width_h_mefisto,
    list(width_latent(list_h_low_mefisto, list_h_upp_mefisto))
  )
  
  
  ## --- Trajectory-level errors (ISE) ---
  list_mat_ise_list_Y_hat <- c(
    list_mat_ise_list_Y_hat,
    lapply(list_res, "[[", "mat_ise_list_Y_hat")
  )
  
  list_mat_ise_list_Y_hat_mefisto <- c(
    list_mat_ise_list_Y_hat_mefisto,
    lapply(list_res, "[[", "mat_ise_list_Y_hat_mefisto")
  )
  
  ## --- Latent factor curve errors (ISE) ---
  list_mat_ise_list_h_hat <- c(
    list_mat_ise_list_h_hat,
    lapply(list_res, "[[", "mat_ise_list_h_hat")
  )
  
  list_mat_ise_list_h_hat_mefisto <- c(
    list_mat_ise_list_h_hat_mefisto,
    lapply(list_res, "[[", "mat_ise_list_h_hat_mefisto")
  )
  
  ## --- CI coverage and widths for Y ---
  list_mat_coverage_ci_Y <- c(
    list_mat_coverage_ci_Y,
    lapply(list_res, "[[", "mat_coverage_ci_Y")
  )
  
  list_mat_mean_width_ci_Y <- c(
    list_mat_mean_width_ci_Y,
    lapply(list_res, "[[", "mat_mean_width_ci_Y")
  )
  
  list_mat_median_width_ci_Y <- c(
    list_mat_median_width_ci_Y,
    lapply(list_res, "[[", "mat_median_width_ci_Y")
  )
  
  list_mat_integrated_width_ci_Y <- c(
    list_mat_integrated_width_ci_Y,
    lapply(list_res, "[[", "mat_integrated_width_ci_Y")
  )
  
  list_mat_coverage_ci_Y_mefisto <- c(
    list_mat_coverage_ci_Y_mefisto,
    lapply(list_res, "[[", "mat_coverage_ci_Y_mefisto")
  )
  
  list_mat_mean_width_ci_Y_mefisto <- c(
    list_mat_mean_width_ci_Y_mefisto,
    lapply(list_res, "[[", "mat_mean_width_ci_Y_mefisto")
  )
  
  list_mat_median_width_ci_Y_mefisto <- c(
    list_mat_median_width_ci_Y_mefisto,
    lapply(list_res, "[[", "mat_median_width_ci_Y_mefisto")
  )
  
  list_mat_integrated_width_ci_Y_mefisto <- c(
    list_mat_integrated_width_ci_Y_mefisto,
    lapply(list_res, "[[", "mat_integrated_width_ci_Y_mefisto")
  )
  
  ## --- Loadings recovery ---
  list_vec_rmse_B_hat <- c(
    list_vec_rmse_B_hat,
    lapply(list_res, "[[", "vec_rmse_B_hat")
  )
  
  list_vec_rmse_B_hat_mefisto <- c(
    list_vec_rmse_B_hat_mefisto,
    lapply(list_res, "[[", "vec_rmse_B_hat_mefisto")
  )

  list_Auc_b <- c(
    list_Auc_b,
    lapply(list_res, "[[", "Auc_b")
  )
  
  list_Auc_b_mefisto <- c(
    list_Auc_b_mefisto,
    lapply(list_res, "[[", "Auc_b_mefisto")
  )
  
  list_Corr_b <- c(
    list_Corr_b,
    lapply(list_res, "[[", "Corr_b")
  )
  
  list_Corr_b_mefisto <- c(
    list_Corr_b_mefisto,
    lapply(list_res, "[[", "Corr_b_mefisto")
  )
  
  ## --- Runtime ---
  list_rt <- c(
    list_rt,
    lapply(list_res, "[[", "rt")
  )
  
  list_rt_mefisto <- c(
    list_rt_mefisto,
    lapply(list_res, "[[", "rt_mefisto")
  )

  vec_nb_factors_retained <- c(
    vec_nb_factors_retained,
    lapply(list_res, "[[", "nb_factors_retained")
  )  
  
  vec_nb_factors_retained_mefisto <- c(
    vec_nb_factors_retained_mefisto,
    lapply(list_res, "[[", "nb_factors_retained_mefisto")
  )
  
  list_factor_ppi <- c(list_factor_ppi, lapply(list_res, "[[", "factor_ppi"))
  
}


ise_Y_bayesSYNC <- collapse_mat(list_mat_ise_list_Y_hat)
ise_Y_mefisto   <- collapse_mat(list_mat_ise_list_Y_hat_mefisto)

summary_ise_Y <- rbind(
  bayesSYNC = mean_se(ise_Y_bayesSYNC),
  MEFISTO   = mean_se(ise_Y_mefisto)
)
print(summary_ise_Y)


## --- Trajectory-level ISE ---
if (bool_save) {
  pdf(paste0(res_dir, "/boxplot_ISE_Y.pdf"), width = 4, height = 4)
}
paired_boxplot(
  ise_Y_bayesSYNC, ise_Y_mefisto,
  main = "ISE(Y)", ylab = "ISE"
)
if (bool_save) {
  dev.off()
}


## --- Coverage ---
cov_Y_bayesSYNC <- collapse_mat(list_mat_coverage_ci_Y)
cov_Y_mefisto   <- collapse_mat(list_mat_coverage_ci_Y_mefisto)

## --- CI width ---
width_Y_bayesSYNC <- collapse_mat(list_mat_integrated_width_ci_Y)
width_Y_mefisto   <- collapse_mat(list_mat_integrated_width_ci_Y_mefisto)

summary_ci_Y <- rbind(
  bayesSYNC_coverage = mean_se(cov_Y_bayesSYNC),
  MEFISTO_coverage   = mean_se(cov_Y_mefisto),
  bayesSYNC_width   = mean_se(width_Y_bayesSYNC),
  MEFISTO_width     = mean_se(width_Y_mefisto)
)
print(summary_ci_Y)

if (bool_save) {
  pdf(paste0(res_dir, "/boxplot_Coverage_Y.pdf"), width = 4, height = 4)
}
paired_boxplot(
  cov_Y_bayesSYNC, cov_Y_mefisto,
  main = "Coverage(Y)", ylab = "Coverage"
)
if (bool_save) {
  dev.off()
}


if (bool_save) {
  pdf(paste0(res_dir, "/boxplot_Width_Y.pdf"), width = 4, height = 4)
}
paired_boxplot(
  width_Y_bayesSYNC, width_Y_mefisto,
  main = "CI width(Y)", ylab = "Width"
)
if (bool_save) {
  dev.off()
}


## --- Latent factor ISE ---
ise_h_bayesSYNC <- collapse_mat(list_mat_ise_list_h_hat)
ise_h_mefisto   <- collapse_mat(list_mat_ise_list_h_hat_mefisto)

summary_ise_h <- rbind(
  bayesSYNC = mean_se(ise_h_bayesSYNC),
  MEFISTO   = mean_se(ise_h_mefisto)
)
print(summary_ise_h)

if (bool_save) {
  pdf(paste0(res_dir, "/boxplot_ISE_h.pdf"), width = 4, height = 4)
}
paired_boxplot(
  ise_h_bayesSYNC, ise_h_mefisto,
  main = "ISE(h)", ylab = "ISE"
)
if (bool_save) {
  dev.off()
}

## --- Loadings RMSE ---

rmse_B_bayesSYNC <- unlist(list_vec_rmse_B_hat)
rmse_B_mefisto   <- unlist(list_vec_rmse_B_hat_mefisto)

summary_rmse_B <- rbind(
  bayesSYNC = mean_se(rmse_B_bayesSYNC),
  MEFISTO   = mean_se(rmse_B_mefisto)
)
print(summary_rmse_B)

if (bool_save) {
  pdf(paste0(res_dir, "/boxplot_RMSE_B.pdf"), width = 4, height = 4)
}
paired_boxplot(
  rmse_B_bayesSYNC, rmse_B_mefisto,
  main = "RMSE(B)", ylab = "RMSE"
)
if (bool_save) {
  dev.off()
}

## --- Loadings correlation ---

corr_B_bayesSYNC <- unlist(list_Corr_b)
corr_B_mefisto   <- unlist(list_Corr_b_mefisto)

## --- Loadings AUC ---
auc_B_bayesSYNC <- unlist(list_Auc_b)
auc_B_mefisto   <- unlist(list_Auc_b_mefisto)

summary_B_recovery <- rbind(
  bayesSYNC_corr = mean_se(corr_B_bayesSYNC),
  MEFISTO_corr   = mean_se(corr_B_mefisto),
  bayesSYNC_auc  = mean_se(auc_B_bayesSYNC),
  MEFISTO_auc    = mean_se(auc_B_mefisto)
)
print(summary_B_recovery)

if (bool_save) {
  pdf(paste0(res_dir, "/boxplot_Corr_B.pdf"), width = 4, height = 4)
}
paired_boxplot(
  corr_B_bayesSYNC, corr_B_mefisto,
  main = "Corr(B)", ylab = "Correlation"
)
if (bool_save) {
  dev.off()
}


if (bool_save) {
  pdf(paste0(res_dir, "/boxplot_Corr_B.pdf"), width = 4, height = 4)
}
paired_boxplot(
  corr_B_bayesSYNC, corr_B_mefisto,
  main = "Corr(B)", ylab = "Correlation"
)
if (bool_save) {
  dev.off()
}


## --- Runtime ---

rt_bayesSYNC <- unlist(list_rt)
rt_mefisto   <- unlist(list_rt_mefisto)

summary_rt <- rbind(
  bayesSYNC = mean_se(rt_bayesSYNC),
  MEFISTO   = mean_se(rt_mefisto)
)
print(summary_rt)


if (bool_save) {
  pdf(paste0(res_dir, "/boxplot_Runtime.pdf"), width = 4, height = 4)
}
paired_boxplot(
  rt_bayesSYNC, rt_mefisto,
  main = "Runtime", ylab = "Seconds"
)
if (bool_save) {
  dev.off()
}


## --- Number of factors retained ---
nb_factors_bayesSYNC <- unlist(vec_nb_factors_retained, use.names = FALSE)
nb_factors_mefisto   <- unlist(vec_nb_factors_retained_mefisto, use.names = FALSE)

summary_nb_factors <- rbind(
  bayesSYNC = mean_se(nb_factors_bayesSYNC),
  MEFISTO   = mean_se(nb_factors_mefisto)
)

print(summary_nb_factors)

if (bool_save) {
  pdf(paste0(res_dir, "/boxplot_Nb_factors_retained.pdf"), width = 4, height = 4)
}
paired_boxplot(
  nb_factors_bayesSYNC,
  nb_factors_mefisto,
  main = "Number of factors retained",
  ylab = "Count"
)
if (bool_save) {
  dev.off()
}

cov_h_bayesSYNC <- collapse_mat(list_mat_coverage_h)
cov_h_mefisto   <- collapse_mat(list_mat_coverage_h_mefisto)

width_h_bayesSYNC <- collapse_mat(list_mat_width_h)
width_h_mefisto   <- collapse_mat(list_mat_width_h_mefisto)

if (bool_save) {
  pdf(paste0(res_dir, "/boxplot_Coverage_h.pdf"), width = 4, height = 4)
}

paired_boxplot(
  cov_h_bayesSYNC,
  cov_h_mefisto,
  main = "Latent factor coverage",
  ylab = "Pointwise coverage"
)

if (bool_save) {
  dev.off()
}


if (bool_save) {
  pdf(paste0(res_dir, "/boxplot_Width_h.pdf"), width = 4, height = 4)
}

paired_boxplot(
  width_h_bayesSYNC,
  width_h_mefisto,
  main = "Latent factor CI width",
  ylab = "Mean width"
)

if (bool_save) {
  dev.off()
}

# Convert list to matrix
prob_mat <- do.call(rbind, list_factor_ppi)

rownames(prob_mat) <- paste0(seq_len(nrow(prob_mat)))
colnames(prob_mat) <- paste0(seq_len(ncol(prob_mat)))

library(reshape2)

df_long <- melt(prob_mat)
colnames(df_long) <- c("Replicate", "Index", "Probability")

library(ggplot2)

if (bool_save) {
  pdf(paste0(res_dir, "/factor_ppi_across_replicates.pdf"), width = 4.5*Q_true, height = 4.5)
}

pp <- ggplot(df_long, aes(x = Index, y = Replicate, fill = Probability)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient(
    low = "grey95",
    high = "black",
    limits = c(0, 1)
  ) +
  scale_y_discrete(limits = rev(levels(df_long$Replicate))) +
  labs(
    title = "Factor posterior probabilities",
    x = "Factor",
    y = "Replicate",
    fill = "Probability"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8)
  )
print(pp)
if (bool_save) {
  dev.off()
}

if (bool_save) {
  pdf(paste0(res_dir, "/histograms_nb_factors.pdf"), width = 4.5*Q_true, height = 4.5)
}

# Common breaks and limits
breaks <- seq(min(c(nb_factors_bayesSYNC, nb_factors_mefisto))-0.5, 
              max(c(nb_factors_bayesSYNC, nb_factors_mefisto))+0.5, by = 1)

h1 <- hist(nb_factors_bayesSYNC, breaks = breaks, plot = FALSE)
h2 <- hist(nb_factors_mefisto,   breaks = breaks, plot = FALSE)

ylim <- c(0, max(c(h1$counts, h2$counts)))

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plot_hist <- function(h, subtitle) {
  plot(h,
       col = "white",
       border = "black",
       xlim = c(-0.5, 4.5),
       ylim = ylim,
       lwd = 1.5,                      # thicker bin contours
       main = paste("Number of factors retained\n", subtitle),
       xlab = "Number of factors",
       ylab = "Count")
  
  # Highlight bin x = 2  (interval [1.5, 2.5])
  idx <- which(h$mids == 2)
  rect(1.5, 0, 2.5, h$counts[idx],
       col = col_mefisto,
       border = NA)
  
  # Draw bins manually
  rect(h$breaks[-length(h$breaks)],
       0,
       h$breaks[-1],
       h$counts,
       col = adjustcolor("white", alpha.f = 0),
       border = "black",
       lwd = 1.5)   # <-- bin contour thickness

}

plot_hist(h1, "BayesSYNC")
plot_hist(h2, "MEFISTO")
if (bool_save) {
  dev.off()
}


if (bool_save) {
  pdf(paste0(res_dir, "/average_ROC_curves.pdf"), width = 4.5*Q_true, height = 4.5)
}
par(mfrow = c(1, Q_true))
for (q in 1:Q_true) {
  print(plot_roc_list(list_ppi = lapply(list_ppi, function(ll) ll[,q]),
                      list_pat = lapply(list_pat, function(ll) ll[,q]),
                      col = c("skyblue2", col_mefisto), main = paste0("Loadings for factor q = ", q),
                      leg = c("bayesSYNC", "MEFISTO"), bool_leg = T))
  print(plot_roc_list(list_ppi = lapply(list_B_hat_mefisto, function(ll) abs(ll[,q])),
                      list_pat = lapply(list_pat, function(ll) ll[,q]),
                      col = col_mefisto, main = paste0("Loadings for factor q = ", q),
                      add = T,  bool_leg = F))
}
if (bool_save) {
  dev.off()
}



par(mfrow = c(1, 1))

library(knitr)

summary_table <- rbind(
  ISE_Y_bayesSYNC          = mean_se(ise_Y_bayesSYNC),
  ISE_Y_MEFISTO            = mean_se(ise_Y_mefisto),
  Coverage_Y_bayesSYNC     = mean_se(cov_Y_bayesSYNC),
  Coverage_Y_MEFISTO       = mean_se(cov_Y_mefisto),
  Width_Y_bayesSYNC        = mean_se(width_Y_bayesSYNC),
  Width_Y_MEFISTO          = mean_se(width_Y_mefisto),
  ISE_h_bayesSYNC          = mean_se(ise_h_bayesSYNC),
  ISE_h_MEFISTO            = mean_se(ise_h_mefisto),
  Coverage_h_bayesSYNC     = mean_se(cov_h_bayesSYNC),
  Coverage_h_MEFISTO       = mean_se(cov_h_mefisto),
  Width_h_bayesSYNC        = mean_se(width_h_bayesSYNC),
  Width_h_MEFISTO          = mean_se(width_h_mefisto),
  RMSE_B_bayesSYNC         = mean_se(rmse_B_bayesSYNC),
  RMSE_B_MEFISTO           = mean_se(rmse_B_mefisto),
  AUC_B_bayesSYNC         = mean_se(auc_B_bayesSYNC),
  AUC_B_MEFISTO           = mean_se(auc_B_mefisto),
  NbFactors_bayesSYNC      = mean_se(nb_factors_bayesSYNC),
  NbFactors_MEFISTO        = mean_se(nb_factors_mefisto),
  Runtime_bayesSYNC        = mean_se(rt_bayesSYNC),
  Runtime_MEFISTO          = mean_se(rt_mefisto)
)

latex_table <- data.frame(
  Metric = rownames(summary_table),
  Value  = sprintf("%.3f (%.3f)", summary_table[, "mean"], summary_table[, "se"])
)

kable(latex_table, format = "latex", booktabs = TRUE)


rt_bayesSYNC <- unlist(list_rt)
rt_mefisto   <- unlist(list_rt_mefisto)

runtime_df <- data.frame(
  seed      = vec_seed,
  bayesSYNC = rt_bayesSYNC,
  MEFISTO   = rt_mefisto,
  N         = N
)

if (bool_save) {
  save(runtime_df, file = file.path(res_dir, "runtime_summary.RData"))
}


metrics_summary <- data.frame(
  Metric = c(
    "AUC loadings",
    "RMSE loadings",
    "ISE(Y)",
    "Coverage(Y)",
    "Width(Y)",
    "Nb factors"
  ),
  bayesSYNC_mean = c(
    mean_se(auc_B_bayesSYNC)["mean"],
    mean_se(rmse_B_bayesSYNC)["mean"],
    mean_se(ise_Y_bayesSYNC)["mean"],
    mean_se(cov_Y_bayesSYNC)["mean"],
    mean_se(width_Y_bayesSYNC)["mean"],
    mean_se(nb_factors_bayesSYNC)["mean"]
  ),
  bayesSYNC_se = c(
    mean_se(auc_B_bayesSYNC)["se"],
    mean_se(rmse_B_bayesSYNC)["se"],
    mean_se(ise_Y_bayesSYNC)["se"],
    mean_se(cov_Y_bayesSYNC)["se"],
    mean_se(width_Y_bayesSYNC)["se"],
    mean_se(nb_factors_bayesSYNC)["se"]
  ),
  MEFISTO_mean = c(
    mean_se(auc_B_mefisto)["mean"],
    mean_se(rmse_B_mefisto)["mean"],
    mean_se(ise_Y_mefisto)["mean"],
    mean_se(cov_Y_mefisto)["mean"],
    mean_se(width_Y_mefisto)["mean"],
    mean_se(nb_factors_mefisto)["mean"]
  ),
  MEFISTO_se = c(
    mean_se(auc_B_mefisto)["se"],
    mean_se(rmse_B_mefisto)["se"],
    mean_se(ise_Y_mefisto)["se"],
    mean_se(cov_Y_mefisto)["se"],
    mean_se(width_Y_mefisto)["se"],
    mean_se(nb_factors_mefisto)["se"]
  ),
  N = N
)

if (bool_save) {
  save(metrics_summary,
       file = file.path(res_dir, "metrics_summary.RData"))
}


# [here produce code to report: boxplot pairs, one for bayesSYNC one for MEFISTO, reporting:
# - Runtime
# - Metrics on the final trajectories: (i) ISE (ii) CI widths and (iii) coverage (averaged across all subjects and all variables)
# - Errors on the latent factor curves 1 and 2: ISE (averaged across all subjects)
# - Errors oadings 1 and 2: RMSE
# also produce code to put these information in a table, each time reporting the mean across replicates along with the se in parentheses. table that has a latex format
