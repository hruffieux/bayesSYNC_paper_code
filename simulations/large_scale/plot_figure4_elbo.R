## ---------------------------------------------------------------------------
## Panels of Figure 4: optimisation behaviour with and without annealing, over
## the 50 replicates of the large-scale simulation.
##
##   left    difference in final ELBO, annealed VB minus vanilla VB, for
##           Q_max = 3 and Q_max = 5, paired by replicate
##   right   total number of iterations to convergence, including the
##           annealing iterations, for the four configurations
##
## The final ELBO and the number of iterations are read from the out.txt
## written by each run of run_large_scale.R, so that the large output.RData
## files need not be loaded.
##
## Usage, from the root of the repository:
##   Rscript simulations/large_scale/plot_figure4_elbo.R
##
## Writes boxplots_difference_elbo.pdf and boxplots_number_of_iterations.pdf in
## a folder BS_r_elbo_... under out_dir.
## ---------------------------------------------------------------------------

rm(list = ls())

source("R/paths.R")

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

vec_Q <- c(3, 5)
list_anneal <- list("no_anneal" = NULL,
                    "anneal_strong" = c(1, 1.9, 100))

bool_var_spec_prob <- FALSE
n_cpus <- 50 # only enters the folder names written by run_large_scale.R

sum_n_repl <- 50
bool_save <- TRUE
if (bool_save) {
  sum_res_dir <- paste0(out_dir, "/BS_r_elbo_",
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

# one column per configuration, in the order (Q = 3, vanilla), (Q = 3,
# annealed), (Q = 5, vanilla), (Q = 5, annealed)
elbo_values <- rt_values <- nb_it_values <- matrix(NA, nrow = sum_n_repl, ncol = length(vec_Q) * length(list_anneal))

i <- 1

for (Q in vec_Q) {

  for (anneal_id in seq_along(list_anneal)) {

    anneal <- list_anneal[[anneal_id]]

    for (SS in 1:sum_n_repl) {

      base_dir <- paste0(out_dir, "/BS_r_anneal_", paste0(anneal, collapse = "-"),
                         "_N_", N, "_p_", p, "_Q_true_", Q_true,
                         "_Q_", Q, "_L_true_", L_true, "_L_", L,
                         "_N_min_", N_min, "_N_max_", N_max,
                         "_K_", K,
                         "_sd_zeta_", paste0(round(vec_sd_zeta, 2), collapse = "-"),
                         "_sd_eps_", paste0(round(unique(vec_sd_eps), 2), collapse = "-"),
                         ifelse(bool_var_spec_prob, "_omega_j_q", ""),
                         "_a_om_", a_om, "_b_om_", b_om, "_tol_abs_", tol_abs,
                         "_tol_rel_", tol_rel, "_maxit_", maxit, "_n_cpus_", max(floor(n_cpus/n_repl), 1),
                         "_n_repl_", n_repl, "_seed_")

      file_path <- file.path(paste0(base_dir, SS), "out.txt")

      if (file.exists(file_path)) {

        lines <- readLines(file_path)

        # runtime: the line below the header of the printed system.time()
        elapsed_index <- grep("elapsed", lines)
        if (length(elapsed_index) > 0 && elapsed_index < length(lines)) {
          elapsed_line <- lines[elapsed_index + 1]
          elapsed_time <- as.numeric(strsplit(elapsed_line, "\\s+")[[1]][3])
          rt_values[SS, i] <- elapsed_time / 60
        } else {
          warning(paste("Elapsed time not found in", file_path))
        }

        # total number of iterations, including the annealing iterations
        nb_it_line <- grep("Convergence obtained after", lines, value = TRUE)
        if (length(nb_it_line) > 0) {
          nb_it_values[SS, i] <- as.numeric(gsub("\\D", "", nb_it_line))
        } else {
          warning(paste("Number of iterations not found in", file_path))
        }

        elbo_line <- grep("Optimal marginal log-likelihood variational lower bound \\(ELBO\\)", lines, value = TRUE)
        if (length(elbo_line) > 0) {
          elbo_values[SS, i] <- as.numeric(sub(".*=\\s*([0-9eE.-]+).*", "\\1", elbo_line))
        } else {
          warning(paste("Final ELBO not found in", file_path))
        }

      } else {
        warning(paste("File not found:", file_path))
      }
    }
    i <- i + 1

  }

}

colnames(rt_values) <- colnames(nb_it_values) <- colnames(elbo_values) <- paste0("Q = ", rep(vec_Q, each = 2), ", ", rep(c("VB", "annealed VB"), length(vec_Q)))

colors <- c(adjustcolor("skyblue1", alpha.f = 0.5),
            "skyblue1",
            adjustcolor("honeydew3", alpha.f = 0.5),
            "honeydew3")

if (bool_save) {
  pdf(paste0(sum_res_dir, "/boxplots_number_of_iterations.pdf"), width = 6, height = 6)
}
boxplot(nb_it_values, col = colors, main = "Total number of iterations", lwd = 1.5, lty = 1, ylim = c(0, 650))
if (bool_save) {
  dev.off()
}

# paired differences, annealed minus vanilla, within each value of Q_max
diff_elbo_Q_3 <- elbo_values[,2] - elbo_values[,1]
diff_elbo_Q_5 <- elbo_values[,4] - elbo_values[,3]

cat("Proportion of replicates with a higher final ELBO under annealing: Q = 3, ",
    format(mean(diff_elbo_Q_3 > 0, na.rm = TRUE), digits = 2), "; Q = 5, ",
    format(mean(diff_elbo_Q_5 > 0, na.rm = TRUE), digits = 2), "\n", sep = "")

if (bool_save) {
  pdf(paste0(sum_res_dir, "/boxplots_difference_elbo.pdf"), width = 3.7, height = 6)
}
boxplot(cbind(diff_elbo_Q_3, diff_elbo_Q_5), col = c(colors[2], colors[4]),
        main = "Difference final ELBO \n(annealed VB - vanilla VB)", lwd = 1.5, lty = 1)
abline(h = 0, lty = 3)
if (bool_save) {
  dev.off()
}
