## ---------------------------------------------------------------------------
## Biomarker panel analysis of the COVID-19 cohort, Section 6.1 of the paper
## (Figures 6 and 7).
##
## 57 circulating biomarkers (CRP, cytokines, complement components,
## glycoproteins, immune cell populations and metabolites) in symptomatic
## PCR-positive patients of severity classes B to E, over the first seven weeks
## after symptom onset (days 0 to 49). The model is fitted with Q_max = 5
## candidate factors and L_max = 5 components, with annealed variational
## inference, and the variables are standardised before the fit.
##
## The real data cannot be made public. By default the script runs on the
## simulated data written by application/simulate_data/simulate_cohort_data.R,
## which have the same format; set bool_simulated_data to FALSE and
## COVID_DATA_DIR to the folder holding the preprocessed data to run it on the
## real data.
##
## Usage, from the root of the repository:
##   Rscript application/biomarker_panel/analysis_biomarker_panel.R [seed] [excl_dead]
##
## The analysis of the paper is the default, with the patients who died within
## the window retained (N = 88 on the real data, seed 1). Giving excl_dead as
## second argument excludes them, which is the sensitivity analysis of
## Supplementary Material S7 (N = 83 on the real data, fitted with seed 20).
##
## Writes, in its own folder under out_dir: output.RData (the whole workspace,
## read by the other scripts of this folder), out.txt, and the panels of
## Figures 6 and 7:
##   B.pdf                          loadings, Figure 6(a)
##   eigenfunctions.pdf             Figure 6(b)
##   scores_1_vs_2.pdf, h_i_extreme_subjects.pdf     Figure 6(c)
##   Y_i_extreme_subjects_..._top_loading_magnitude.pdf, Figure 7
##   boxplots_Factor_*_FPC_*_conditional_on_severity_*.pdf  scores by class,
##                                  with the overall ANOVA p-value
## The gene set enrichment and over-representation analyses of the loadings on
## KEGG pathways are printed to out.txt.
## ---------------------------------------------------------------------------

rm(list = ls())

source("R/paths.R")
source("R/fun_application.R")
source("R/plot_utils.R")

require(tidyr)
require(dplyr)
require(bayesSYNC)

args <- commandArgs(trailingOnly = TRUE)
seed <- if (length(args) >= 1) as.integer(args[1]) else 1
set.seed(seed)

# TRUE: drop the patients who died within the window (sensitivity analysis);
# FALSE: keep them (analysis reported in the paper)
excl_dead <- length(args) >= 2 && args[2] == "excl_dead"

bool_simulated_data <- TRUE

type <- "metabo_cyto_cplt_glyco_main_ct" # CRP, cytokines, complement, glycoproteins, main cell types and metabolites

if (bool_simulated_data) {

  # df_comb (one row per sample: identifiers, clinical variables and markers),
  # df_info (clinical variables), all_var (the 57 markers) and
  # metabolite_names (the 36 metabolites, used for the enrichment analysis)
  load(file.path(data_dir, "simulated_biomarker_panel.RData"))
  type <- paste0("simulated_", type)

} else {

  covid_data_dir <- Sys.getenv("COVID_DATA_DIR")

  version_ms <- "20210630"
  version_nmr <- "201120"
  version_glyc <- "010721"
  version_lipid <- "20210614"
  version_ques <- "20210505"
  version_info <- "20210427_fu"
  version_flow <- "HC2_210429"
  version_cytof <- "201005"

  keep_replicates <- FALSE

  load(file.path(covid_data_dir,
                 paste0("preprocessed_data/metabolites_ms_v", version_ms,
                        "_nmr_v", version_nmr,
                        "_glyc_v", version_glyc,
                        "_lipids_v", version_lipid,
                        "_clin_v", version_info,
                        "_ques_v", version_ques,
                        ifelse(keep_replicates, "_with_replicates", ""),
                        ".RData")))

  load(file.path(covid_data_dir,
                 paste0("preprocessed_data/cell_types_flow_v", version_flow,
                        "_cytof_v", version_cytof,
                        "_clin_v", version_info, ".RData")))

  df_comb <- merge(df_ms_comb, df_flow_main, by="row.names")

  all_var <- c("log_CRP", "IL10", "IL1B", "IL6", "C3a", "C3c", "TCC", "GlycA", "GlycB", "SPC",
               sort(names(df_flow_main)), sort(names(df_ms)))

  metabolite_names <- names(df_ms)
}


days_thres <- 50
tnew <- 0:(days_thres-1) # timepoints to evaluate the prediction at

min_timepoints <- 1

selected_severity_groups <- c("B", "C", "D", "E") # symptomatic PCR-positive patients

df_comb_no_HC <- choose_severity_groups(df_comb, selected_severity_groups)
df_comb_no_HC <- df_comb_no_HC[df_comb_no_HC$days_from_sx_or_swab_imputed < days_thres, ]


# ---------------------------------------------------------------------------
# Vital status: sensitivity analysis excluding subjects who died
# (trajectories of subjects who died are informatively truncated at the time
# of death)
# ---------------------------------------------------------------------------

# excl_dead is set from the command line above

# "in_window": drop only the subjects who died within the days_thres window, i.e.
#              those whose observed trajectory is actually truncated by death.
#              These are the subjects affected by the truncation; the default.
# "any":       drop every subject with hospital outcome "dead", whenever death
#              occurred, including after the end of the window.
dead_definition <- c("in_window", "any")[1]

# Which values of hospital_outcome count as a death. The recorded categories are
# dead, discharged, hospital, other_facility, other_hospital, palliative, and a
# sizeable number of missing values.
#   - "dead" is unambiguous.
#   - "palliative" is not a recorded death but is not a survival either: the
#     subject may have died without the field being updated. Add it here to test
#     the sensitivity of the conclusions to that choice.
#   - a missing outcome is an unknown status, not a survival. Such subjects are
#     counted separately below and the death counts are lower bounds.
dead_labels <- c("dead")

# Subject-level vital status table built from the clinical data.
# Note: the plotting functions (plot_trajectory, plot_factor_trajectory, ...)
# force hospital_outcome to "alive" for severity classes HC, A and B. That
# recoding is not applied here, so that the counts below reflect the recorded
# outcome for every subject, including those in class B.
get_vital_status <- function(df_info, days_thres, dead_labels) {

  stopifnot(all(c("subject_id", "severity", "hospital_outcome",
                  "hospital_outcome_days_from_start_sx_or_first_pos_swab",
                  "last_observation_days_from_start_sx_or_first_pos_swab") %in% names(df_info)))

  # one row per subject (df_info has one row per sample)
  df_status <- df_info[!duplicated(df_info$subject_id),
                       c("subject_id", "severity", "hospital_outcome",
                         "hospital_outcome_days_from_start_sx_or_first_pos_swab",
                         "last_observation_days_from_start_sx_or_first_pos_swab"),
                       drop = FALSE]

  df_status$dead <- df_status$hospital_outcome %in% dead_labels
  df_status$unknown_status <- is.na(df_status$hospital_outcome) # unknown, not a survival

  # day of death, counted from symptom onset or first positive swab
  df_status$time_death <- rep(NA_real_, nrow(df_status))
  df_status$time_death[df_status$dead] <-
    df_status$hospital_outcome_days_from_start_sx_or_first_pos_swab[df_status$dead]

  # time of the event: death for those who died, last observation otherwise.
  # Used by the plotting functions to cut the fitted curve at the time of death.
  df_status$time_event <- df_status$last_observation_days_from_start_sx_or_first_pos_swab
  df_status$time_event[df_status$dead] <- df_status$time_death[df_status$dead]

  # death occurring inside the analysis window, i.e. trajectory truncated by death.
  # A subject who died after the window has a fully observed trajectory over the
  # window and is therefore not concerned by the informative truncation issue.
  df_status$dead_in_window <- df_status$dead & !is.na(df_status$time_death) &
    df_status$time_death < days_thres

  # a subject recorded as dead but with no date cannot be classified
  if (any(df_status$dead & is.na(df_status$time_death))) {
    warning(paste0(sum(df_status$dead & is.na(df_status$time_death)),
                   " subject(s) recorded as dead but with no date of death; ",
                   "counted as dying outside the window.\n"))
  }

  df_status

}

df_status <- get_vital_status(df_info, days_thres, dead_labels)

if (excl_dead) {

  stopifnot(dead_definition %in% c("in_window", "any"))
  subj_to_drop <- df_status$subject_id[if (dead_definition == "any") df_status$dead else df_status$dead_in_window]

  n_before <- length(unique(df_comb_no_HC$subject_id))
  df_comb_no_HC <- df_comb_no_HC[!(df_comb_no_HC$subject_id %in% subj_to_drop), , drop = FALSE]
  n_after <- length(unique(df_comb_no_HC$subject_id))

  cat(paste0("Excluding subjects who died (definition: ", dead_definition, "): ",
             n_before - n_after, " subjects removed before data preparation.\n"))

}


bool_rescale_time <- TRUE
data <- prepare_data_irregular(df_comb_no_HC, all_var, min_nb_timepoints = min_timepoints,
                               days_thres = days_thres, bool_rescale_time = bool_rescale_time,
                               Z = NULL)

list_time_obs_Y <- enforce_common_grid(data$time_obs, data$Y, min_nb_timepoints = min_timepoints)

time_obs <- list_time_obs_Y$time_obs
Y <- list_time_obs_Y$Y
N <- length(Y)
p <- data$p
N
p


# ---------------------------------------------------------------------------
# How many of the N subjects actually entering the model died?
# Y is built after all filtering steps (severity classes, time window, missing
# values, common grid, minimum number of time points), so names(Y) is the
# definitive list of analysed subjects.
# ---------------------------------------------------------------------------

df_status_Y <- df_status[match(names(Y), df_status$subject_id), , drop = FALSE]
stopifnot(nrow(df_status_Y) == N)

if (any(is.na(df_status_Y$subject_id))) {
  warning(paste0(sum(is.na(df_status_Y$subject_id)),
                 " analysed subject(s) not found in df_info; their vital status is unknown.\n"))
}

cat("\n=== VITAL STATUS OF THE ", N, " ANALYSED SUBJECTS ===\n", sep = "")
cat("Death defined as hospital_outcome in: ", paste0(dead_labels, collapse = ", "), "\n", sep = "")
cat("Deaths (any time): ", sum(df_status_Y$dead, na.rm = TRUE), " / ", N, "\n", sep = "")
cat("Deaths within the ", days_thres, "-day window: ",
    sum(df_status_Y$dead_in_window, na.rm = TRUE), " / ", N, "\n", sep = "")
cat("Vital status unknown (missing hospital outcome): ",
    sum(df_status_Y$unknown_status, na.rm = TRUE), " / ", N,
    ". The two counts above are lower bounds.\n", sep = "")

cat("\nRecorded hospital outcome:\n")
print(table(df_status_Y$hospital_outcome, useNA = "ifany"))

cat("\nBreakdown by severity class (rows) and recorded outcome (columns):\n")
print(table("severity" = droplevels(factor(df_status_Y$severity)),
            "outcome" = df_status_Y$hospital_outcome, useNA = "ifany"))

cat("\nNumber of observations per subject, by vital status:\n")
print(tapply(sapply(time_obs, length), df_status_Y$dead, summary))

if (any(df_status_Y$dead_in_window, na.rm = TRUE)) {
  cat("\nDay of death (from symptom onset or first positive swab) for deaths within the window:\n")
  print(sort(df_status_Y$time_death[which(df_status_Y$dead_in_window)]))
}
if (any(df_status_Y$dead & !df_status_Y$dead_in_window, na.rm = TRUE)) {
  cat("\nDay of death for deaths occurring after the window (trajectories not truncated):\n")
  print(sort(df_status_Y$time_death[which(df_status_Y$dead & !df_status_Y$dead_in_window)]))
}
cat("=====================================\n\n")


subject_ids <- names(Y)


n_obs <- sapply(time_obs, length)
N_min <- min(n_obs)
N_max <- max(n_obs)
N_min
N_max

L <- 5
K <- NULL

Q <- 5 # over-specified; the factors are pruned by their posterior inclusion probabilities

anneal <- c(1, 1.9, 100)

tol_abs <- 1e-3
tol_rel <- 1e-6
maxit <- 10000

n_cpus <- 1 # used inside bayesSYNC

bool_var_spec_prob <- F
bool_scale <- T # whether to put all the variables on a similar scale

bool_save <- T

if (bool_save) {
  res_dir <- paste0(out_dir, "/R1_covid_",
                    ifelse(excl_dead, paste0("excl_dead_", dead_definition, "_"), ""), # sensitivity analysis excluding subjects who died
                    "anneal_", paste0(anneal, collapse = "-"), "_",
                    "min_timepts_", min_timepoints, "_",
                    ifelse(bool_var_spec_prob, "omega_j_q_", ""), type,
                    "_N_", N, "_p_", p,
                    "_Q_", Q, "_L_", L, "_N_min_", N_min, "_N_max_",
                    N_max, "_K_", K,
                    "_tol_abs_", tol_abs, "_tol_rel_",
                    tol_rel, "_maxit_", maxit, "_scale_", bool_scale,
                    "_n_cpus_", n_cpus, "_seed_", seed, "/")
  dir.create(res_dir, showWarnings = FALSE)
  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"), type = "message")
  output_file <- file.path(res_dir, "output.RData")
} else {
  res_dir <- output_file <- NULL
}


if (!is.null(output_file) && file.exists(output_file)) {
  # a completed fit with the same settings is reused rather than refitted;
  # delete its folder to refit, for instance after regenerating the data
  # printed with cat, since messages are diverted to err.txt
  cat("Reusing the existing fit in ", output_file, " instead of refitting. ",
      "Delete this folder to refit.\n", sep = "")
  load(output_file)
  source("R/fun_application.R")
  source("R/plot_utils.R")
} else {

  show_progress <- F
  runtime_bayesSYNC <- system.time(res <- bayesSYNC(time_obs = time_obs,
                                                    Y = Y, L = L, Q = Q,
                                                    K = K,
                                                    anneal = anneal,
                                                    tol_abs = tol_abs,
                                                    tol_rel = tol_rel,
                                                    maxit = maxit,
                                                    seed = seed,
                                                    n_cpus = n_cpus,
                                                    bool_scale = bool_scale,
                                                    bool_var_spec_prob = bool_var_spec_prob,
                                                    show_factor_ppi_progress = show_progress))

  res_original <- res

  # The sign of a factor is not determined by the data. On the real data
  # (N = 88, seed 1) the two retained factors are factors 4 and 5, and their
  # signs are set so that the inflammatory markers load positively on the
  # first factor of the paper (factor 5) and the amino acids negatively on the
  # second (factor 4); the latent trajectories are reversed with them, through
  # the eigenfunctions for factor 5 and the scores for factor 4. The factor
  # numbers belong to that fit, so nothing is reversed on the simulated data.
  if (!bool_simulated_data && !excl_dead) {

    vec_id <- c(5)
    for (id in vec_id) {
      res$list_list_Phi_hat[[id]] <- -1*res$list_list_Phi_hat[[id]]
      res$B_hat[,id] <- -1*res$B_hat[,id]

      for (i in seq_along(res$list_h_hat)) {

        res$list_h_low[[i]][[id]] <- -1*res$list_h_low[[i]][[id]]
        res$list_h_hat[[i]][[id]] <- -1*res$list_h_hat[[i]][[id]]
        res$list_h_upp[[i]][[id]] <- -1*res$list_h_upp[[i]][[id]]

      }
    }

    vec_id_b <- c(4)
    for (id_b in vec_id_b) {

      res$B_hat[,id_b] <- -1*res$B_hat[,id_b]
      res$list_Zeta_hat[[id_b]] <- -1*res$list_Zeta_hat[[id_b]]

      for (i in seq_along(res$list_h_hat)) {

        res$list_h_low[[i]][[id_b]] <- -1*res$list_h_low[[i]][[id_b]]
        res$list_h_hat[[i]][[id_b]] <- -1*res$list_h_hat[[i]][[id_b]]
        res$list_h_upp[[i]][[id_b]] <- -1*res$list_h_upp[[i]][[id_b]]

      }
    }
  }

  if (bool_scale) { # in order to be able to display the HC bands on the correct scale.
    for (vv in all_var) {
      df_comb[, vv] <-  (df_comb[, vv] - res$mean_mean_across_subjects[vv]) / res$sd_mean_across_subjects[vv]
    }
  }

}

print(runtime_bayesSYNC)

# saved as soon as the fit is available, so that a failure in the figures below
# does not lose it; saved again at the end with the enrichment results
if (bool_save) {
  save.image(output_file)
}

print(paste0("Factor PPIs: ",
             paste0(format(res$factor_ppi, digits = 2), collapse = ", "))) # expect to see the largest values first

set.seed(seed)

disp <-  Q %/% 5 + 1 # integer division
if (bool_save) {
  pdf(paste0(res_dir, "/scores_1_vs_2.pdf"), height = 3*disp, width = 2.5*ceiling(Q/disp))
}

vec_col_grey <- rep("white", N)
names(vec_col_grey) <- as.character(df_info$severity[match(names(Y), df_info$subject_id)])

vec_col_grey[names(vec_col_grey) %in% "A"] <- "grey90"
vec_col_grey[names(vec_col_grey) %in% "B"] <- "grey80"
vec_col_grey[names(vec_col_grey) %in% "C"] <- "grey65"
vec_col_grey[names(vec_col_grey) %in% "D"] <- "grey45"
vec_col_grey[names(vec_col_grey) %in% "E"] <- "grey5"

par(mfrow = c(disp, ceiling(Q/disp)), pty = "s")
null <- sapply(1:Q, function(q) plot(res$list_Zeta_hat[[q]][,1], res$list_Zeta_hat[[q]][,2],
                                     main = paste0("Scores factor q = ", q),
                                     xlab = "Scores FPC 1", ylab = "Scores FPC 2", pch = 20,
                                     col = vec_col_grey))
if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "/mean_functions.pdf"), width = 8.5, height = 6)
}
par(mfrow = c(2, 3))
sapply(sample(1:p, 6), function(j) plot_mean_functions(j, res$time_g, mu_g = NULL,
                                                       res$list_mu_hat))
if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "/eigenfunctions.pdf"), width = 3*L, height = 3*Q)
}
cum_pve_thres <- 0.99
par(mfrow = c(Q, L)) # only showing functions up to 100 % of variability explained reached.
null <- sapply(1:Q, function(q)
  plot_factor_eigenfunctions(q, res$time_g, Psi_g = NULL,
                             res$list_list_Phi_hat,
                             res$list_cumulated_pve,
                             cum_pve_thres = cum_pve_thres))
if (bool_save) {
  dev.off()
}

n_random <- 6
vec_random_i <- sample(1:N, n_random)
vec_random_j <- sample(1:p, n_random)
if (bool_save) {
  pdf(paste0(res_dir, "/traj_i_", paste0(vec_random_i, collapse = "-"), "_j_",
             paste0(vec_random_j, collapse = "-"), ".pdf"),
      width = 10, height = 6)
}
par(mfrow = c(2, 3))
for (id in 1:n_random) {
  random_i <- vec_random_i[id]
  random_j <- vec_random_j[id]
  plot_trajectory(random_i, random_j,
                  time_obs, res$time_g, res$Y,
                  res$list_Y_hat, res$list_Y_low, res$list_Y_upp, days_thres = days_thres,
                  df_comb = df_comb)
}
if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "/h_i_extreme_subjects.pdf"),
      width = 9, height = 2.6*Q)
}

par(mfrow = c(Q, 4))
for (q in 1:Q) {
  for (l in 1:2) {
    plot_factor_trajectory(res$list_Zeta_hat, q, l,
                           res$time_g,
                           res$list_h_hat, res$list_h_low, res$list_h_upp,
                           days_thres = days_thres, df_comb = df_comb)
  }

}
if (bool_save) {
  dev.off()
}
for (q in 1:Q) {
  thres_ppi <- 0.5
  if (bool_save) {
    pdf(paste0(res_dir, "/Y_i_extreme_subjects_ppi_thres_", thres_ppi, "_variables_contributing_to_factor_q_", q, ".pdf"),
        width = 9, height = 13)
  }

  plot_factor_variables_trajectory(res$list_Zeta_hat, res$ppi, q, time_obs,
                                   res$time_g, res$Y, res$list_Y_hat,
                                   res$list_Y_low, res$list_Y_upp, thres_ppi = thres_ppi,
                                   df_comb = df_comb,
                                   days_thres = days_thres, show_ppi = FALSE)

  if (bool_save) {
    dev.off()
  }
}


# heatmap.2 only works for matrices with at least two columns and rows, and
# needs distinct breaks, so not when every loading is zero (all factors pruned)
if (Q > 1 && max(abs(res$B_hat)) > 0) {

  require(gplots)

  trunc <- 800
  pal_length <-  trunc
  my_palette <-  colorRampPalette(c("skyblue1", "white", "#ff7f6e"))(n = pal_length)

  if (bool_save) {
    pdf(paste0(res_dir, "/B", ifelse(p>200, "_truncated", ""), ".pdf"),
        width = 2*Q, height = 12, paper='special')
  }
  d <- res$B_hat[1:min(p, trunc), ]
  m <- max(abs(d))
  breaks <- seq(-m, m, length.out = length(my_palette) + 1)

  heatmap.2(d,
            breaks = breaks,
            dendrogram='none',
            Rowv=FALSE,
            Colv=FALSE,
            trace='none',
            density.info="none",
            main = "Factor loadings",
            ylab = "",
            xlab = "Factor q",
            col=my_palette,  # Set the custom color palette
            key=T,       # Show the color key
            keysize = 1.0,  # Adjust the size of the color key
            key.title = "Posterior mean of B",
            key.xlab = "",
            margins = c(8, 16),
            cexRow = 1,   # Adjust the size of row labels
            cexCol = 1.0,   # Adjust the size of column labels
            labRow = gsub("log_", "", rownames(d)),
            labCol = paste0("PPI: ", format(res$factor_ppi, digits = 2)),
            scale="none")   # Turn off scaling

  if (bool_save) {
    dev.off()
  }
}

thres <- 0.5
for (q in 1:Q) {
  if (any(res$ppi[,q]>thres)) {
    top_ppi_q <- res$ppi[,q][res$ppi[,q]>thres]
    top_B_hat_q <- res$B_hat[,q][res$ppi[,q]>thres]
    o <- order(top_ppi_q, decreasing = TRUE)
    print(paste0("Top PPI multiplied by sign B_hat for factor q = ", q, ": "))
    print(top_ppi_q[o]*sign(top_B_hat_q[o]))

    o_magnitude <- order(abs(top_B_hat_q), decreasing = TRUE)
    print(top_B_hat_q[o_magnitude])

    if (bool_save) {
      pdf(paste0(res_dir, "/Y_i_extreme_subjects_ppi_thres_", thres_ppi, "_variables_contributing_to_factor_q_", q, "_top_loading_magnitude.pdf"),
          width = 9, height = 13)
    }
    vec_variables <- names(top_B_hat_q[o_magnitude])[1:6]
    plot_factor_variables_trajectory(res$list_Zeta_hat, res$ppi, q, time_obs,
                                     res$time_g, res$Y, res$list_Y_hat,
                                     res$list_Y_low, res$list_Y_upp, thres_ppi = thres_ppi,
                                     df_comb = df_comb,vec_variables = vec_variables,
                                     days_thres = days_thres, show_ppi = FALSE)

    if (bool_save) {
      dev.off()
    }
  }
}

print("omega_hat")
print(res$omega_hat)

if (bool_save) {
  save.image(output_file)
}


print(res$list_cumulated_pve)

# FPC scores by severity class, with the overall ANOVA p-value (Section 6.1)
method <- "t.test"
bool_log10_display <- F

scores <- NULL
for (q in 1:Q) {
  scores <- cbind(scores, res$list_Zeta_hat[[q]][,1:2])
  colnames(scores)[(ncol(scores)-1):ncol(scores)] <- paste0("Factor_", q, "_FPC_", 1:2)
}
vec_var <- colnames(scores)

df_scores <- data.frame(scores, "severity" = names(vec_col_grey))
res_tests <- plot_boxplots(df_scores,
                           vec_var = vec_var,
                           group_var = "severity",
                           vec_var_to_log = NULL,
                           method = method,
                           main = "Scores FPC by severity class",
                           xlab = "Severity class",
                           ylab = "Scores",
                           vec_col = vec_col_grey,
                           save_path = res_dir,
                           bool_log10_disp = bool_log10_display,
                           bool_overall = T,
                           width = 4.6,
                           height = 3.7)


bool_enrichment <- TRUE
min_pathway_size <- 5
max_pathway_size <- 200 # Limiting the maximum size of pathways is typically recommended to remove overly broad or generic pathways that may be statistically significant simply because they contain many genes, not because they represent a biologically coherent process.

if (bool_enrichment) {

  # The KEGG pathway library is retrieved through the KEGG REST API on first
  # use and cached in data_dir. Some computing clusters block the API (HTTP
  # error 403), in which case run this step once on a machine with access.
  file_kegg <- file.path(data_dir, "kegg_pathway_library.RData")

  if (file.exists(file_kegg)) {
    load(file_kegg)
  } else {
    df_library <- retrieve_pathway_library()
    save(df_library, file = file_kegg)
  }

  cat("=== ENRICHMENT ANALYSIS FOR METABOLITES CONTRIBUTING TO FACTORS ===\n",
      "=== with pathway size filter:",  min_pathway_size, "===\n")

  for (factor_id in c(1:Q)) {

    cat(paste0("----- FACTOR: ", factor_id, " -----\n",
               "--- factor PPI: ", format(res$factor_ppi[factor_id], digits = 3), " ---\n",
               "---------------------\n"))

    for (useGSEA in c(T, F)) {

      if (useGSEA) {

        cat("*GSEA analysis (ranking based on estimated loadings)*\n")
        target_list <- metabolite_names

        bool_ppi_based <- F
        if (bool_ppi_based) {
          rank_list <- res$ppi[metabolite_names, factor_id]
          scoreType <- "pos"
        } else {
          rank_list <- abs(res$B_hat[metabolite_names, factor_id])
          scoreType <-  "pos" # "std"
        }

        background_list <- NULL

        no_results <- (length(unique(rank_list)) == 1)

      } else {

        # default now manual for both
        cat("\n*Hypergeometric test with background (selection based on PPI > 0.5)*\n")

        background_list <- metabolite_names
        thres <- 0.5
        target_list <- intersect(background_list,
                                 rownames(res$ppi)[res$ppi[, factor_id]>thres]) # intersect to use only metabolites, no CRP, cell types etc

        cat("Target list (PPI > 0.5): \n")
        print(target_list)

        rank_list <- NULL
        scoreType <- NULL

        no_results <- (length(target_list) == 0)
      }

      if (no_results) {
        cat("Inactive factor: no enrichement analysis performed.\n")
      } else {

        set.seed(seed)

        df_enrichment <- do_enrichment_analysis(target_list,
                                                background_list = background_list,
                                                rank_list = rank_list,
                                                pathway_library = df_library$pathway_library,
                                                pathway_names = df_library$pathway_names,
                                                term2gene = df_library$term2gene,
                                                useGSEA = useGSEA,
                                                min_pathway_size = min_pathway_size,
                                                max_pathway_size = max_pathway_size,
                                                scoreType = scoreType,
                                                kegg_dir = data_dir,
                                                bool_strict_metabolic_pathway_library = T)

        if (!is.null(df_enrichment)) {
          print(df_enrichment[df_enrichment$pvalue<0.05,])
        }

      }

    }

  }

}


if (bool_save) {
  save.image(output_file)
}

