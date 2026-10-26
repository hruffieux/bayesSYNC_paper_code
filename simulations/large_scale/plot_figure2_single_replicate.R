## ---------------------------------------------------------------------------
## Panels of Figure 2: qualitative recovery on a single replicate of the
## large-scale simulation (N = 100, p = 20,000, Q = 3, L = 2), fitted with
## Q_max = 5 and annealed VB. This is replicate 1, job_id = 151 of
## run_large_scale.R.
##
##   (a) observed and fitted trajectories with 95% pointwise prediction bands,
##       and the noise-free simulated trajectories, for subsets of subjects and
##       variables drawn at random (one file per random draw, 25 draws)
##   (b) estimated loadings of the first 50 variables, with the signs of the
##       nonzero simulated loadings marked
##   (c) simulated and estimated eigenfunctions
##   (d) simulated and estimated scores
##
## Supplementary Figure S1 (simulated and estimated loadings side by side) is
## drawn from the same fit by plot_figure_S_loadings.R.
##
## Usage, from the root of the repository:
##   Rscript simulations/large_scale/plot_figure2_single_replicate.R
##
## Writes the panels as pdf files in a folder R1_plot_2026_single_repl_... under
## out_dir.
## ---------------------------------------------------------------------------

rm(list = ls())

source("R/paths.R")
source("R/plot_utils.R")
source("R/utils_sim.R")

require(bayesSYNC)
require(gplots)

seed <- 1
set.seed(seed)

p <- 20000
N <- 100
N_min <- 5
N_max <- 10

Q_true <- 3
Q <- 5
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

anneal <- c(1, 1.9, 100)

n_repl <- 1

bool_var_spec_prob <- FALSE

n_cpus <- 50 # only enters the folder names written by run_large_scale.R
bool_save <- TRUE

if (bool_save) {
  res_dir <- paste0(out_dir, "/R1_plot_2026_single_repl_mu_diff_N_", N, "_p_", p, "_Q_true_", Q_true,
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

repl_file <- paste0(out_dir, "/BS_r_anneal_", paste0(anneal, collapse = "-"),
                    "_N_", N, "_p_", p, "_Q_true_", Q_true,
                    "_Q_", Q, "_L_true_", L_true, "_L_", L,
                    "_N_min_", N_min, "_N_max_", N_max,
                    "_K_", K,
                    "_sd_zeta_", paste0(round(vec_sd_zeta, 2), collapse = "-"),
                    "_sd_eps_", paste0(round(unique(vec_sd_eps), 2), collapse = "-"),
                    ifelse(bool_var_spec_prob, "_omega_j_q", ""),
                    "_a_om_", a_om, "_b_om_", b_om, "_tol_abs_", tol_abs,
                    "_tol_rel_", tol_rel, "_maxit_", maxit, "_n_cpus_", max(floor(n_cpus/n_repl), 1),
                    "_n_repl_", n_repl, "_seed_", seed, "/output.RData")

if (file.exists(repl_file)) {
  load(repl_file)
} else {
  stop(paste0("Output of the single replicate not found: ", repl_file))
}

n_g <- 1000
time_g <- seq(0, 1, length.out = n_g)


## ---------------------------------------------------------------------------
## Noise-free simulated trajectories. The generator returns the mean functions
## and the eigenfunctions evaluated on its own grid of n_g points, which is the
## grid used here, so
##
##   y_i^(j)(t) = mu^(j)(t) + sum_q b_jq sum_l zeta_il^(q) psi_l^(q)(t)
##
## reproduces the signal that the observations were drawn around, before the
## observation error was added. Held as a list of lists matching list_Y_hat.
## ---------------------------------------------------------------------------

compute_true_trajectories <- function(res, time_g, vec_i, vec_j) {

  stopifnot(!is.null(res$B), !is.null(res$Zeta), !is.null(res$mu_g), !is.null(res$Phi_g))
  stopifnot(length(res$mu_g[[1]]) == length(time_g))

  Q <- length(res$Zeta)

  # latent trajectories on the grid, one N x n_g matrix per factor
  list_h_g <- lapply(seq_len(Q), function(q) res$Zeta[[q]] %*% t(res$Phi_g[[q]]))

  list_Y_true <- vector("list", length(res$Y))
  names(list_Y_true) <- names(res$Y)

  for (i in vec_i) {
    Y_true_i <- vector("list", length(res$Y[[i]]))
    names(Y_true_i) <- names(res$Y[[i]])
    for (j in vec_j) {
      shared <- Reduce(`+`, lapply(seq_len(Q), function(q) res$B[j, q] * list_h_g[[q]][i, ]))
      Y_true_i[[j]] <- res$mu_g[[j]] + shared
    }
    list_Y_true[[i]] <- Y_true_i
  }

  list_Y_true
}


## (a) trajectories, for 25 random draws of five subjects and five variables;
## the panel of the paper shows one of them

for (seed in 1:25) {

  set.seed(seed)

  n_random <- 5
  vec_random_i <- sort(sample(1:N, n_random))
  vec_random_j <- sort(sample(1:p, n_random))

  # only the subjects and variables actually plotted are reconstructed, since
  # doing so for all p = 20000 variables on a grid of 1000 points is needless
  list_Y_true <- compute_true_trajectories(list_res[[1]], time_g,
                                           vec_random_i, vec_random_j)
  if (bool_save) {
    pdf(paste0(res_dir, "/data_i_", paste0(vec_random_i, collapse = "-"), "_j_",
               paste0(vec_random_j, collapse = "-"), "_repl.pdf"),
        width = 12.5, height = 2.65)
  }
  par(mfrow = c(1, 5))
  for (id in 1:n_random) {
    random_i <- vec_random_i[id]
    random_j <- vec_random_j[id]
    plot_trajectory(random_i, random_j,
                    list_res[[1]]$time_obs, time_g, list_res[[1]]$Y,
                    list_res[[1]]$list_Y_hat, list_res[[1]]$list_Y_low, list_res[[1]]$list_Y_upp,
                    bool_disp_curves = F)
  }
  if (bool_save) {
    dev.off()
    pdf(paste0(res_dir, "/traj_i_", paste0(vec_random_i, collapse = "-"), "_j_",
               paste0(vec_random_j, collapse = "-"), "_repl.pdf"),
        width = 12.5, height = 2.65)
  }
  par(mfrow = c(1, 5))
  for (id in 1:n_random) {
    random_i <- vec_random_i[id]
    random_j <- vec_random_j[id]
    plot_trajectory(random_i, random_j,
                    list_res[[1]]$time_obs, time_g, list_res[[1]]$Y,
                    list_res[[1]]$list_Y_hat, list_res[[1]]$list_Y_low, list_res[[1]]$list_Y_upp,
                    Y_true = list_Y_true)
  }
  if (bool_save) {
    dev.off()
  }
}

## (d) scores
if (bool_save) {
  pdf(paste0(res_dir, "/scores_simulated_vs_estimated.pdf"), width = 4*L_true, height = 3*Q_true)
}
par(mfrow = c(Q_true, L_true))
null <- sapply(1:Q_true, function(q)
  sapply(1:L_true, function(l) plot_scores(q, l,
                                           list_res[[1]]$Zeta,
                                           list_res[[1]]$list_Zeta_hat,
                                           list_sd_hat = list_res[[1]]$list_sd_hat)))

## (c) eigenfunctions, shown up to 99% of the variance explained
if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "/eigenfunctions.pdf"), width = 3*L, height = 3*Q_true)
}
par(mfrow = c(Q_true, L))
null <- sapply(1:Q_true, function(q)
  sapply(1:L, function(l) {
    if (l <= L_true) {
      res_repl_Psi_g <- list_res[[1]]$Phi_g
    } else {
      res_repl_Psi_g <- NULL
    }
    plot_eigenfunctions(q, l, time_g, Psi_g = res_repl_Psi_g, list_res[[1]]$list_list_Phi_hat,
                        list_res[[1]]$list_cumulated_pve,
                        cum_pve_thres = 0.99, col_sim = "black", col = "skyblue1")
  })
)
if (bool_save) {
  dev.off()
}


## (b) estimated loadings of the first 50 variables, on a scale with a white
## band around zero, and the signs of the nonzero simulated loadings

trunc <- 50
pal_length <-  trunc

n_white <- round(0.25 * pal_length)
n_side  <- (pal_length - n_white) / 2

my_palette <- c(
  colorRampPalette(c("skyblue1", "white"))(n_side),
  rep("white", n_white),
  colorRampPalette(c("white", "#ff7f6e"))(n_side)
)

true_B_trunc <- list_res[[1]]$B[1:min(p, trunc), ]
B_hat_trunc <- list_res[[1]]$B_hat[1:min(p, trunc), ]

true_sign_trunc <- matrix("", nrow = nrow(true_B_trunc), ncol = ncol(true_B_trunc))
true_sign_trunc[true_B_trunc > 0] <- "+"
true_sign_trunc[true_B_trunc < 0] <- "-"

d <- B_hat_trunc
rownames(d) <- paste0(1:nrow(d))
colnames(d) <- paste0(1:ncol(d))

# width of the white band around zero. Taken over the nonzero entries only:
# the loadings are sparse, so the 10th percentile of all absolute values is
# exactly zero as soon as a tenth of them are, which would leave duplicate
# break points and heatmap.2 rejects those. Floored away from zero for the same
# reason.
vec_abs_nonzero <- abs(d)[!is.na(d) & d != 0]
eps <- quantile(vec_abs_nonzero, 0.1, na.rm = TRUE)
eps <- max(eps, 1e-6 * max(abs(d), na.rm = TRUE))

breaks <- c(
  seq(min(d, na.rm = TRUE), -eps, length.out = n_side + 1),
  seq(-eps,  eps, length.out = n_white + 1)[-1],
  seq(eps, max(d, na.rm = TRUE), length.out = n_side + 1)[-1]
)

if (bool_save) {
  pdf(paste0(res_dir, "/estimated_B", ifelse(p>trunc, "_truncated", ""), ".pdf"),
      width = 1.4*Q_true, height = 20, paper='special')
}

par(mar = c(5, 4, 4, 8))

heatmap.2(d,
          breaks = breaks,
          dendrogram='none',
          Rowv=FALSE,
          Colv=FALSE,
          trace='none',
          main = "Est. factor loadings",
          ylab = "Variable",
          xlab = "Factor",
          margin = c(8, 16),
          cellnote = true_sign_trunc,
          notecex=1.0,
          notecol="black",
          col=my_palette,
          key=F,
          keysize = 1.0,
          cexRow = 1.0,
          cexCol = 1.0,
          scale="none")

par(mar = c(5, 4, 4, 2) + 0.1)

if (bool_save) {
  dev.off()
}
