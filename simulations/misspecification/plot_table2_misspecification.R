## ---------------------------------------------------------------------------
## Manuscript table and figure for the misspecification study.
##
## Reads df_all.RData written by gather_misspecification.R and produces the
## two objects the paper needs, rather than the diagnostic grid that script
## draws for checking runs.
##
## The reporting is deliberately split in two, because averaging accuracy over
## replicates that recovered every factor together with replicates that missed
## one mixes outcomes that are not comparable. The mixture favours whichever
## setting over-retains, since the matching step then has more candidates to
## align against the truth: on a plain average the correlated-error settings
## appear to recover the eigenfunctions twenty times better than the
## well-specified baseline, which is an artefact of retention and not a result.
##
##   Panel and block 1: how many factors were retained, and how often all of
##     them were found. A result in its own right, and the headline for the
##     correlated-error scenario.
##   Panel and block 2: accuracy, restricted to the replicates that recovered
##     every true factor, so that the comparison is like for like.
##
## Usage, from the root of the repository:
##   Rscript simulations/misspecification/plot_table2_misspecification.R
##
## Writes table_misspecification.tex (Table 2), figure_misspecification.pdf and
## figure_residual_correlogram.pdf in R1_misspec_summary under out_dir.
## ---------------------------------------------------------------------------

rm(list = ls())

source("R/paths.R")

sum_res_dir <- file.path(out_dir, "R1_misspec_summary/")

load(file.path(sum_res_dir, "df_all.RData"))


## ---------------------------------------------------------------------------
## Tidying
##
## The summary rows are assembled column by column, so numeric fields can arrive
## as lists. Coerce once, here, rather than at every use.
## ---------------------------------------------------------------------------

as_num <- function(x) suppressWarnings(as.numeric(unlist(x)))

vec_num <- c("auc_b", "corr_zeta_avg", "ise_eigen_avg", "ise_mean_per_variable",
             "coverage_zeta_avg", "slope_zeta_avg", "n_retained", "seed",
             "sigsq_eps_mean", "min_norm_col_B_hat", "retention", "mean_n_obs",
             "corr_resid_bin1", "corr_resid_bin2", "corr_resid_bin3")
for (nm in intersect(vec_num, names(df_all))) df_all[[nm]] <- as_num(df_all[[nm]])

vec_corr_cols <- grep("^corr_zeta_q[0-9]+_l[0-9]+$", names(df_all), value = TRUE)
Q_true <- max(as.integer(sub("^corr_zeta_q([0-9]+)_l[0-9]+$", "\\1", vec_corr_cols)))

# a true factor is undetected when its matched column of B_hat is essentially
# zero, which is what min_norm_col_B_hat records
df_all$all_detected <- df_all$min_norm_col_B_hat >= 1e-6

# ordering and display names for the paper
df_order <- data.frame(
  label = c("baseline", "rho = 0.15", "rho = 0.4", "kappa = 0.1", "kappa = 0.3",
            "factor corr = 0.3", "factor corr = 0.6",
            "non_informative, beta = 0", "informative, beta = 2"),
  group = c("", rep("Correlated errors", 2), rep("Unmodelled latent variation", 2),
            rep("Correlated factors", 2), rep("Drop-out", 2)),
  pretty = c("-",
             "$\\rho = 0.15$", "$\\rho = 0.40$",
             "$\\kappa = 0.1$", "$\\kappa = 0.3$",
             "$\\varrho = 0.3$", "$\\varrho = 0.6$",
             "Non-informative", "Informative"),
  plain = c("Baseline", "rho = 0.15", "rho = 0.40", "kappa = 0.1", "kappa = 0.3",
            "fac. corr = 0.3", "fac. corr = 0.6", "Non-inform.", "Informative"),
  stringsAsFactors = FALSE)
df_order <- df_order[df_order$label %in% unique(df_all$label), ]

df_all$label <- factor(df_all$label, levels = df_order$label)


## ---------------------------------------------------------------------------
## Summaries
## ---------------------------------------------------------------------------

mean_se <- function(x) {
  x <- x[!is.na(x)]
  c(mean = mean(x), se = sd(x) / sqrt(length(x)), n = length(x))
}

summarise_by <- function(metric, subset = rep(TRUE, nrow(df_all))) {
  do.call(rbind, lapply(df_order$label, function(lab) {
    v <- df_all[[metric]][df_all$label == lab & subset]
    as.data.frame(t(mean_se(v)))
  }))
}

# block 1, every run
df_ret <- do.call(rbind, lapply(df_order$label, function(lab) {
  d <- df_all[df_all$label == lab, ]
  data.frame(n_ret = mean(d$n_retained, na.rm = TRUE),
             n_ret_se = sd(d$n_retained, na.rm = TRUE) / sqrt(sum(!is.na(d$n_retained))),
             prop_all = mean(d$all_detected, na.rm = TRUE),
             n_runs = nrow(d))
}))

# block 2, replicates that recovered every true factor
ok <- df_all$all_detected %in% TRUE
vec_acc <- c("auc_b", "ise_eigen_avg", "corr_zeta_avg", "coverage_zeta_avg",
             "sigsq_eps_mean", "corr_resid_bin1")
list_acc <- lapply(vec_acc, summarise_by, subset = ok)
names(list_acc) <- vec_acc


## ---------------------------------------------------------------------------
## Table fragment for the manuscript
##
## Written as a fragment so that main.tex keeps the caption and the label and
## the numbers are never retyped. Booktabs is already loaded there.
## ---------------------------------------------------------------------------

fmt <- function(m, se, digits = 3) {
  ifelse(is.na(m), "---",
         sprintf(paste0("%.", digits, "f (%.", digits, "f)"), m, se))
}

# mechanism names are split over the two rows of each block, so that the first
# column stays narrow enough for the table to fit the text width at \scriptsize,
# which the manuscript sets around the fragment
list_group_lines <- list("Correlated errors" = c("Correlated", "errors"),
                         "Unmodelled latent variation" = c("Unmodelled", "latent variation"),
                         "Correlated factors" = c("Correlated", "factors"),
                         "Drop-out" = c("Drop-out", ""))

vec_lines <- character(0)
add <- function(...) vec_lines <<- c(vec_lines, paste0(...))

add("% Generated by plot_table2_misspecification.R. Do not edit by hand.")
add("\\begin{tabular}{llrrrrr}")
add("\\toprule")
add(" & & \\multicolumn{2}{c}{Factors} & \\multicolumn{3}{c}{Accuracy, runs recovering all factors} \\\\")
add("\\cmidrule(lr){3-4} \\cmidrule(lr){5-7}")
add("Mechanism & Setting & Retained & All found & AUC & Eigenfn.\\ ISE & Score corr. \\\\")
add("\\midrule")

group_prev <- NA
k_in_group <- 0
for (k in seq_len(nrow(df_order))) {
  g <- df_order$group[k]
  if (!identical(g, group_prev)) {
    if (k > 1) add("\\addlinespace")
    k_in_group <- 0
  }
  k_in_group <- k_in_group + 1
  group_prev <- g
  if (!nzchar(g)) {
    cell_g <- "Baseline" # the baseline row carries no setting of its own
  } else {
    lines_g <- list_group_lines[[g]]
    if (is.null(lines_g)) lines_g <- g # a mechanism without a split keeps its full name
    cell_g <- if (k_in_group <= length(lines_g)) lines_g[k_in_group] else ""
  }
  add(cell_g, " & ", df_order$pretty[k],
      " & ", sprintf("%.2f", df_ret$n_ret[k]),
      " & ", sprintf("%.0f\\%%", 100 * df_ret$prop_all[k]),
      " & ", fmt(list_acc$auc_b$mean[k], list_acc$auc_b$se[k]),
      " & ", fmt(list_acc$ise_eigen_avg$mean[k], list_acc$ise_eigen_avg$se[k]),
      " & ", fmt(list_acc$corr_zeta_avg$mean[k], list_acc$corr_zeta_avg$se[k]),
      " \\\\")
}
add("\\bottomrule")
add("\\end{tabular}")

writeLines(vec_lines, file.path(sum_res_dir, "table_misspecification.tex"))
cat("written: table_misspecification.tex\n")


## ---------------------------------------------------------------------------
## Figure
##
## Four panels, one per part of the argument, rather than one panel per metric.
## Points are means with Monte Carlo standard errors over replicates; the
## baseline is drawn as a horizontal reference line so that each departure is
## read against it.
## ---------------------------------------------------------------------------

col_group <- c("grey35", "#d95f02", "#1b9e77", "#7570b3", "#e7298a")
names(col_group) <- c("", "Correlated errors", "Unmodelled latent variation",
                      "Correlated factors", "Drop-out")
vec_col <- vapply(df_order$group,
                  function(g) if (!nzchar(g)) "grey35" else unname(col_group[[g]]),
                  character(1))

panel <- function(m, se, ylab, main, baseline_ref = TRUE, ylim = NULL) {
  x <- seq_along(m)
  if (is.null(ylim)) {
    ylim <- range(c(m - 2 * se, m + 2 * se), na.rm = TRUE)
    ylim <- ylim + c(-1, 1) * 0.08 * diff(ylim)
  }
  plot(x, m, type = "n", xaxt = "n", xlab = "", ylab = ylab, main = main,
       ylim = ylim, cex.main = 1.0)
  if (baseline_ref) abline(h = m[1], col = "grey70", lty = 2)
  ind <- which(se > 0)
  if (length(ind) > 0) {
    arrows(x[ind], (m - se)[ind], x[ind], (m + se)[ind], angle = 90, code = 3,
           length = 0.03, col = vec_col[ind], lwd = 1.5)
  }
  points(x, m, pch = 19, col = vec_col, cex = 1.1)
  axis(1, at = x, labels = FALSE)
  text(x, par("usr")[3] - 0.06 * diff(par("usr")[3:4]),
       labels = df_order$plain, srt = 40, adj = 1, xpd = TRUE, cex = 0.8)
}

pdf(file.path(sum_res_dir, "figure_misspecification.pdf"), width = 9, height = 7.5)
par(mfrow = c(2, 2), mar = c(6.5, 4.2, 3, 1))

panel(df_ret$n_ret, df_ret$n_ret_se, "Factors retained",
      "(a) Number of factors retained", baseline_ref = FALSE,
      ylim = c(2, 5.3))
abline(h = Q_true, col = "grey70", lty = 3)

panel(list_acc$ise_eigen_avg$mean, list_acc$ise_eigen_avg$se,
      "Integrated squared error", "(b) Eigenfunction error")

panel(list_acc$corr_zeta_avg$mean, list_acc$corr_zeta_avg$se,
      "Correlation with the truth", "(c) Subject-level scores")

panel(list_acc$coverage_zeta_avg$mean, list_acc$coverage_zeta_avg$se,
      "Empirical coverage", "(d) Coverage of 95% score intervals")

invisible(dev.off())
cat("written: figure_misspecification.pdf\n")


## ---------------------------------------------------------------------------
## Residual correlogram across settings
##
## The diagnostic the application is asked to report, shown here on data where
## the answer is known: flat at zero when the errors really are uncorrelated,
## positive and decaying with the lag when they are not.
## ---------------------------------------------------------------------------

mat_resid <- sapply(c("corr_resid_bin1", "corr_resid_bin2", "corr_resid_bin3"),
                    function(cn) tapply(df_all[[cn]], df_all$label, mean, na.rm = TRUE))
mat_resid <- mat_resid[df_order$label, , drop = FALSE]

pdf(file.path(sum_res_dir, "figure_residual_correlogram.pdf"), width = 5.5, height = 4.5)
par(mar = c(4.2, 4.2, 2.5, 1))
matplot(1:3, t(mat_resid), type = "b", pch = 19, lty = 1, col = vec_col,
        xaxt = "n", ylim = range(c(0, mat_resid), na.rm = TRUE) * 1.1,
        xlab = "Lag between observation times", ylab = "Residual correlation",
        main = "Residual correlation by lag")
axis(1, at = 1:3, labels = c("short", "medium", "long"))
abline(h = 0, col = "grey70", lty = 2)
legend("topright", legend = df_order$plain,
       col = vec_col, pch = 19, bty = "n", cex = 0.7)
invisible(dev.off())
cat("written: figure_residual_correlogram.pdf\n")
