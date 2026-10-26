## ---------------------------------------------------------------------------
## Gather the misspecification runs across settings and replicates.
##
## Each job of the array writes a one-row summary, so the gathering step is a
## plain rbind. The paired structure of the design is exploited twice: the
## across-replicate averages are reported with their Monte Carlo standard
## errors, and the differences with respect to the baseline setting are
## computed within replicate, which is the comparison with the smallest Monte
## Carlo error.
##
## Usage, from the root of the repository:
##   Rscript simulations/misspecification/gather_misspecification.R
##
## Writes df_all.RData, summary_metrics.csv, paired_differences.csv and a
## diagnostic grid of boxplots in R1_misspec_summary under out_dir.
## ---------------------------------------------------------------------------

rm(list = ls())

source("R/paths.R")

sum_res_dir <- file.path(out_dir, "R1_misspec_summary/")
dir.create(sum_res_dir, recursive = TRUE, showWarnings = FALSE)


## ---------------------------------------------------------------------------
## Collect the one-row summaries
## ---------------------------------------------------------------------------

vec_dirs <- list.dirs(out_dir, recursive = FALSE)
vec_dirs <- vec_dirs[grepl("^R1_misspec_", basename(vec_dirs))]
vec_dirs <- vec_dirs[!grepl("_summary", basename(vec_dirs))]

list_rows <- lapply(vec_dirs, function(dir_run) {
  file_run <- file.path(dir_run, "summary_row.RData")
  if (!file.exists(file_run)) {
    return(NULL) # job not finished or failed
  }
  load(file_run) # loads summary_row
  summary_row
})

vec_missing <- vec_dirs[sapply(list_rows, is.null)]
if (length(vec_missing) > 0) {
  print(paste0("Runs without a summary: ", length(vec_missing)))
  print(basename(vec_missing))
}

df_all <- do.call(rbind, list_rows[!sapply(list_rows, is.null)])

# label identifying a setting, used for grouping
df_all$label <- with(df_all,
                     ifelse(scenario == "corr_error",
                            paste0("rho = ", rho),
                            ifelse(scenario == "latent_error",
                                   paste0("kappa = ", kappa),
                                   ifelse(scenario == "latent_corr",
                                          paste0("factor corr = ", rho_fac),
                                          ifelse(scenario == "dropout",
                                                 paste0(dropout, ", beta = ", beta_dropout),
                                                 "baseline")))))

# every scenario must map to a label of its own: a scenario missing from the
# chain above would fall through to "baseline", which would both merge it into
# the baseline group and corrupt the paired differences below, since those match
# on the seed within df_all$label == "baseline"
stopifnot(!any(df_all$scenario != "baseline" & df_all$label == "baseline"))

save(df_all, file = file.path(sum_res_dir, "df_all.RData"))

print(table(df_all$label))

# runs where the component alignment was not the identity, for which the score
# coverage should be interpreted with care
if (any(df_all$alignment_ok == 0, na.rm = TRUE)) {
  print("Runs with a non-identity component alignment, by setting:")
  print(tapply(df_all$alignment_ok == 0, df_all$label, sum, na.rm = TRUE))
}

# runs where a factor was left unrescaled because its estimated loadings were
# essentially zero, i.e. where a true factor was not detected
if (any(df_all$min_norm_col_B_hat < 1e-6, na.rm = TRUE)) {
  print("Runs with an undetected factor, by setting:")
  print(tapply(df_all$min_norm_col_B_hat < 1e-6, df_all$label, sum, na.rm = TRUE))
}


## ---------------------------------------------------------------------------
## Across-replicate averages with Monte Carlo standard errors
## ---------------------------------------------------------------------------

vec_metrics <- c("auc_b", "ise_mean_per_variable", "ise_eigen_avg",
                 "corr_zeta_avg", "slope_zeta_avg", "coverage_zeta_avg",
                 "n_retained", "sigsq_eps_mean",
                 "corr_resid_bin1", "corr_resid_bin2", "corr_resid_bin3",
                 "retention", "mean_n_obs", "runtime_min")

# Q_true and L_true are not stored explicitly; they are recoverable from the
# per-factor, per-component columns written by flatten() in build_summary_row()
vec_corr_cols <- grep("^corr_zeta_q[0-9]+_l[0-9]+$", names(df_all), value = TRUE)
Q_true <- max(as.integer(sub("^corr_zeta_q([0-9]+)_l[0-9]+$", "\\1", vec_corr_cols)))
L_true <- max(as.integer(sub("^corr_zeta_q[0-9]+_l([0-9]+)$", "\\1", vec_corr_cols)))
print(paste0("Q_true = ", Q_true, ", L_true = ", L_true))


## ---------------------------------------------------------------------------
## Factor recovery, reported separately from accuracy
##
## Averaging accuracy over replicates that found every factor and replicates
## that missed one mixes two qualitatively different outcomes, and the mixture
## favours whichever setting over-retains, since the matching step then has more
## candidates to align against the truth. Correlated errors push the retained
## count to the maximum, so on a plain average they appear to recover the
## eigenfunctions better than the well-specified baseline, which is an artefact.
## Retention is therefore reported on its own, and the accuracy metrics are
## reported again below restricted to the runs that recovered every factor.
## ---------------------------------------------------------------------------

df_all$all_detected <- df_all$min_norm_col_B_hat >= 1e-6
df_all$exact_Q <- df_all$n_retained == Q_true

df_recovery <- do.call(rbind, lapply(split(df_all, df_all$label), function(df_lab) {
  data.frame(label = df_lab$label[1],
             n_runs = nrow(df_lab),
             mean_n_retained = mean(df_lab$n_retained, na.rm = TRUE),
             median_n_retained = median(df_lab$n_retained, na.rm = TRUE),
             prop_exact_Q = mean(df_lab$exact_Q, na.rm = TRUE),
             prop_under = mean(df_lab$n_retained < Q_true, na.rm = TRUE),
             prop_over = mean(df_lab$n_retained > Q_true, na.rm = TRUE),
             prop_all_detected = mean(df_lab$all_detected, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
rownames(df_recovery) <- NULL
write.csv(df_recovery, file.path(sum_res_dir, "factor_recovery.csv"), row.names = FALSE)
print(df_recovery)

summarise_metric <- function(metric) {
  agg_mean <- tapply(df_all[[metric]], df_all$label, mean, na.rm = TRUE)
  agg_se <- tapply(df_all[[metric]], df_all$label,
                   function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  agg_med <- tapply(df_all[[metric]], df_all$label, median, na.rm = TRUE)
  data.frame(label = names(agg_mean),
             metric = metric,
             mean = as.vector(agg_mean),
             median = as.vector(agg_med),
             mc_se = as.vector(agg_se),
             stringsAsFactors = FALSE)
}

df_summary <- do.call(rbind, lapply(vec_metrics, summarise_metric))
write.csv(df_summary, file.path(sum_res_dir, "summary_metrics.csv"), row.names = FALSE)

print(df_summary)

# the same summaries restricted to the runs that recovered every true factor, so
# that the accuracy comparison is not driven by how many factors were retained
df_full <- df_all[df_all$all_detected %in% TRUE, ]
print(paste0("Runs with every factor detected: ", nrow(df_full), " of ", nrow(df_all)))

summarise_metric_cond <- function(metric) {
  agg_mean <- tapply(df_full[[metric]], df_full$label, mean, na.rm = TRUE)
  agg_med <- tapply(df_full[[metric]], df_full$label, median, na.rm = TRUE)
  agg_se <- tapply(df_full[[metric]], df_full$label,
                   function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  agg_n <- tapply(df_full[[metric]], df_full$label, function(x) sum(!is.na(x)))
  data.frame(label = names(agg_mean), metric = metric,
             n = as.vector(agg_n), mean = as.vector(agg_mean),
             median = as.vector(agg_med), mc_se = as.vector(agg_se),
             stringsAsFactors = FALSE)
}

df_summary_cond <- do.call(rbind, lapply(vec_metrics, summarise_metric_cond))
write.csv(df_summary_cond,
          file.path(sum_res_dir, "summary_metrics_full_recovery.csv"), row.names = FALSE)
print(df_summary_cond)


## ---------------------------------------------------------------------------
## Paired differences with respect to the baseline
##
## For each replicate the misspecified run and the baseline run share the same
## simulated loadings, scores, eigenfunctions and observation times, so the
## within-replicate difference isolates the effect of the misspecification.
## ---------------------------------------------------------------------------

df_base <- df_all[df_all$label == "baseline", ]

paired_difference <- function(metric) {
  
  vec_labels <- setdiff(unique(df_all$label), "baseline")
  
  do.call(rbind, lapply(vec_labels, function(lab) {
    
    df_lab <- df_all[df_all$label == lab, ]
    ind_match <- match(df_lab$seed, df_base$seed)
    diff_lab <- df_lab[[metric]] - df_base[[metric]][ind_match]
    diff_lab <- diff_lab[!is.na(diff_lab)]
    
    data.frame(label = lab,
               metric = metric,
               n_paired = length(diff_lab),
               mean_diff = mean(diff_lab),
               mc_se = sd(diff_lab) / sqrt(length(diff_lab)),
               stringsAsFactors = FALSE)
  }))
}

df_paired <- do.call(rbind, lapply(vec_metrics, paired_difference))
write.csv(df_paired, file.path(sum_res_dir, "paired_differences.csv"), row.names = FALSE)

print(df_paired)


## ---------------------------------------------------------------------------
## Informative against non-informative drop-out
##
## The comparison with the baseline confounds two things, the loss of
## observations and the dependence of that loss on the latent trajectory. Both
## drop-out arms are calibrated to the same expected retention, so the
## within-replicate difference between them isolates informativeness alone. This
## is the contrast that answers the comment on subjects whose follow-up is
## truncated by death, and it is not obtainable from the differences above.
## ---------------------------------------------------------------------------

lab_inf <- grep("^informative", unique(df_all$label), value = TRUE)
lab_non <- grep("^non_informative", unique(df_all$label), value = TRUE)

if (length(lab_inf) > 0 && length(lab_non) > 0) {

  df_non <- df_all[df_all$label == lab_non[1], ]

  dropout_contrast <- function(metric) {
    do.call(rbind, lapply(lab_inf, function(lab) {
      df_lab <- df_all[df_all$label == lab, ]
      ind_match <- match(df_lab$seed, df_non$seed)
      diff_lab <- df_lab[[metric]] - df_non[[metric]][ind_match]
      diff_lab <- diff_lab[!is.na(diff_lab)]
      data.frame(label = paste0(lab, " vs ", lab_non[1]),
                 metric = metric,
                 n_paired = length(diff_lab),
                 mean_diff = mean(diff_lab),
                 mc_se = sd(diff_lab) / sqrt(length(diff_lab)),
                 stringsAsFactors = FALSE)
    }))
  }

  df_dropout <- do.call(rbind, lapply(vec_metrics, dropout_contrast))
  write.csv(df_dropout, file.path(sum_res_dir, "dropout_contrast.csv"), row.names = FALSE)
  print(df_dropout)

} else {
  print("No drop-out arms found; the informativeness contrast is skipped.")
}


## ---------------------------------------------------------------------------
## Boxplots across settings
## ---------------------------------------------------------------------------

vec_order <- c("baseline",
               sort(grep("^rho = ", unique(df_all$label), value = TRUE)),
               sort(grep("^kappa = ", unique(df_all$label), value = TRUE)),
               sort(grep("^factor corr = ", unique(df_all$label), value = TRUE)),
               sort(grep("beta", unique(df_all$label), value = TRUE)))
vec_order <- vec_order[vec_order %in% unique(df_all$label)]

df_all$label <- factor(df_all$label, levels = vec_order)

pdf(file.path(sum_res_dir, "boxplots_misspecification.pdf"),
    width = 4 * 3, height = 4 * ceiling(length(vec_metrics) / 3))
par(mfrow = c(ceiling(length(vec_metrics) / 3), 3))
for (metric in vec_metrics) {
  boxplot(df_all[[metric]] ~ df_all$label,
          xlab = "", ylab = metric, main = metric, las = 2, col = "white")
}
dev.off()