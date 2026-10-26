## ---------------------------------------------------------------------------
## Simulated data emulating the two COVID-19 data sets of Section 6, which
## cannot be made public:
##
##   - a biomarker panel of 57 circulating markers on 88 patients, one to four
##     samples each (Section 6.1);
##   - whole-blood transcriptomic data on 22,354 genes, of which 13,177 are
##     expressed in blood, on 164 patients, one to four samples each
##     (Section 6.2), with 60 healthy-control samples.
##
## The two data sets share 83 patients, as the real ones do. Both are written in
## the format read by the analysis scripts, so that every step of the
## application can be run on them unchanged.
##
## The data are generated from the model of the paper, with the parameters
## estimated on the real data, which are read from the two files in
## calibration/ (see make_calibration.R): for each variable, the estimated mean
## function, loadings on the retained factors and residual standard deviation;
## for each factor, the estimated eigenfunctions; and, for each factor and
## component, the class means and the within-class standard deviation of the
## estimated scores, so that the latent processes are associated with severity
## roughly as in the real data. The numbers of patients per severity class, of
## samples per patient and the distribution of sampling days within each class
## follow the real design.
##
## What the simulation does not reproduce:
##   - only the analysed design is simulated, so every marker is observed at
##     every sample of the panel (the real panel has missing values, and only
##     the 88 patients with complete data enter the analysis);
##   - deaths are drawn with the class-specific frequencies of the real data
##     and are not informative: a simulated death occurs after the last sample
##     and does not depend on the latent trajectory;
##   - the two data sets of a shared patient are related only through the
##     severity class;
##   - age, sex and body mass index are drawn independently of everything else,
##     and there is no batch effect, although the batch labels are drawn with
##     the real proportions;
##   - the genes not expressed in blood are drawn as noise around a common
##     detection floor, and the healthy controls of both data sets are drawn
##     independently for each variable;
##   - the voom precision weights are all set to one.
##
## Results obtained on these data illustrate the use of the code and are not
## expected to reproduce the numbers of the paper.
##
## Usage, from the root of the repository:
##   Rscript application/simulate_data/simulate_cohort_data.R
## The seed is fixed at 8, which gives the data sets described in the README
## and the biomarker panel analysed in Supplementary Material S10.
##
## Writes, in data_dir (data/ by default):
##   simulated_biomarker_panel.RData   df_comb, df_info, all_var, metabolite_names
##   simulated_biomarker_panel_truth.rds  generating loadings, eigenfunctions and
##                                     scores of the panel (Supplementary S10)
##   simulated_transcriptomics.rds     list with E, weights, targets, genes and
##                                     clinical, as the real EList (about 100 Mb)
## ---------------------------------------------------------------------------

rm(list = ls())

source("R/paths.R")

seed <- 8 # the data set of Supplementary Material S10
set.seed(seed)

cal_dir <- file.path("application", "simulate_data", "calibration")
cal_panel <- readRDS(file.path(cal_dir, "biomarker_panel_calibration.rds"))
cal_rnaseq <- readRDS(file.path(cal_dir, "transcriptomics_calibration.rds"))

classes <- c("B", "C", "D", "E")
days_thres <- cal_panel$days_thres
stopifnot(days_thres == cal_rnaseq$days_thres)

n_common <- 83 # patients in both data sets
n_hc_rnaseq <- 60 # healthy-control samples of the transcriptomic data, 28 of which carry no identifier
n_hc_rnaseq_no_id <- 28


## ---------------------------------------------------------------------------
## Cohort: patients, severity classes and vital status
##
## The shared patients are allocated to the classes in proportion to the
## smaller of the two class sizes, so that both data sets keep their class
## sizes exactly.
## ---------------------------------------------------------------------------

n_panel <- as.vector(cal_panel$design$n_patients[classes])
n_rnaseq <- as.vector(cal_rnaseq$design$n_patients[classes])
names(n_panel) <- names(n_rnaseq) <- classes

n_min <- pmin(n_panel, n_rnaseq)
o_c <- floor(n_min * n_common / sum(n_min))
while (sum(o_c) < n_common) {
  frac <- n_min * n_common / sum(n_min) - o_c
  frac[o_c >= n_min] <- -Inf
  o_c[which.max(frac)] <- o_c[which.max(frac)] + 1
}

df_cohort <- do.call(rbind, lapply(classes, function(cc) {
  data.frame(severity = cc,
             in_panel = c(rep(TRUE, o_c[cc]), rep(TRUE, n_panel[cc] - o_c[cc]), rep(FALSE, n_rnaseq[cc] - o_c[cc])),
             in_rnaseq = c(rep(TRUE, o_c[cc]), rep(FALSE, n_panel[cc] - o_c[cc]), rep(TRUE, n_rnaseq[cc] - o_c[cc])),
             stringsAsFactors = FALSE)
}))
df_cohort <- df_cohort[sample(nrow(df_cohort)), ]
df_cohort$subject_id <- sprintf("SIM%04d", seq_len(nrow(df_cohort)))
rownames(df_cohort) <- df_cohort$subject_id

# class-specific frequency of death within the window, pooled over the two
# real data sets
dead_panel <- cal_panel$dead_in_window[classes, "TRUE"]
dead_rnaseq <- cal_rnaseq$dead_in_window[classes, "TRUE"]
prob_dead <- (dead_panel + dead_rnaseq) / (n_panel + n_rnaseq)
df_cohort$dead <- rbinom(nrow(df_cohort), 1, prob_dead[df_cohort$severity]) == 1

df_cohort$gender <- sample(c("F", "M"), nrow(df_cohort), replace = TRUE, prob = c(0.45, 0.55))
df_cohort$age <- round(pmin(pmax(rnorm(nrow(df_cohort), 60, 15), 18), 95))
df_cohort$bmi <- round(pmax(rnorm(nrow(df_cohort), 28, 5), 16), 1)


## ---------------------------------------------------------------------------
## Sampling days, by data set, from the real class-specific distributions of
## the number of samples per patient and of the sampling days. A small mass is
## added to every day so that a patient can always be given distinct days.
## ---------------------------------------------------------------------------

draw_days <- function(cc, design) {
  n_i <- sample(seq_len(ncol(design$n_obs_by_class)), 1, prob = design$n_obs_by_class[cc, ])
  prob_day <- design$days_by_class[cc, ] + 0.1
  sort(sample(0:(days_thres - 1), n_i, prob = prob_day))
}

list_days_panel <- lapply(df_cohort$subject_id[df_cohort$in_panel], function(s) {
  draw_days(df_cohort[s, "severity"], cal_panel$design)
})
names(list_days_panel) <- df_cohort$subject_id[df_cohort$in_panel]

list_days_rnaseq <- lapply(df_cohort$subject_id[df_cohort$in_rnaseq], function(s) {
  draw_days(df_cohort[s, "severity"], cal_rnaseq$design)
})
names(list_days_rnaseq) <- df_cohort$subject_id[df_cohort$in_rnaseq]

# day of death after the last sample of the patient, within the window; the
# end of follow-up of the others after the end of the window
last_day <- sapply(df_cohort$subject_id, function(s) {
  max(c(list_days_panel[[s]], list_days_rnaseq[[s]]))
})
df_cohort$time_death <- ifelse(df_cohort$dead,
                               pmin(last_day + sample(1:7, nrow(df_cohort), replace = TRUE), days_thres - 1),
                               NA)
df_cohort$hospital_outcome <- ifelse(df_cohort$dead, "dead", "discharged")
df_cohort$hospital_outcome_days_from_start_sx_or_first_pos_swab <-
  ifelse(df_cohort$dead, df_cohort$time_death, last_day + sample(1:20, nrow(df_cohort), replace = TRUE))
df_cohort$last_observation_days_from_start_sx_or_first_pos_swab <-
  ifelse(df_cohort$dead, df_cohort$time_death, days_thres + sample(10:100, nrow(df_cohort), replace = TRUE))


## ---------------------------------------------------------------------------
## Observations from the model
##
##   y_ij(t) = m_j + s_j [ mu_j(t) + sum_q b_jq h_iq(t) + eps_ij(t) ],
##   h_iq(t) = sum_l zeta_iql psi_ql(t),
##
## on the standardised scale of the fit, mapped back to the original scale of
## each variable with its mean m_j and standard deviation s_j. The scores are
## Gaussian with the class mean and the within-class standard deviation of the
## estimated scores.
## ---------------------------------------------------------------------------

# linear interpolation weights of the points t on the grid time_g
interp_weights <- function(time_g, t) {
  k <- findInterval(t, time_g, all.inside = TRUE)
  w <- (t - time_g[k]) / (time_g[k + 1] - time_g[k])
  list(k = k, w = w)
}

interp_rows <- function(mat_g, iw) {
  (1 - iw$w) * mat_g[iw$k, , drop = FALSE] + iw$w * mat_g[iw$k + 1, , drop = FALSE]
}

simulate_values <- function(cal, list_days, time_g_Phi = cal$time_g) {

  p <- nrow(cal$B)
  Q <- ncol(cal$B)

  list_Y <- lapply(names(list_days), function(s) {

    cc <- df_cohort[s, "severity"]
    t <- list_days[[s]] / (days_thres - 1)

    list_zeta <- vector("list", Q)
    H <- sapply(seq_len(Q), function(q) {
      sc <- cal$list_scores[[q]]
      zeta <- sc$class_mean[cc, ] + sc$within_sd * rnorm(length(sc$within_sd))
      list_zeta[[q]] <<- zeta
      as.vector(interp_rows(cal$list_Phi[[q]], interp_weights(time_g_Phi, t)) %*% zeta)
    })
    H <- matrix(H, nrow = length(t))

    Y_std <- interp_rows(cal$mu_g, interp_weights(cal$time_g, t)) + H %*% t(cal$B) +
      sweep(matrix(rnorm(length(t) * p), nrow = length(t)), 2, cal$sd_eps, "*")

    Y_i <- sweep(sweep(Y_std, 2, cal$sd_var, "*"), 2, cal$mean_var, "+")
    attr(Y_i, "zeta") <- list_zeta
    Y_i
  })

  # simulated scores, one list per patient with one vector per factor, kept so
  # that estimates can be compared with the truth
  list_zeta <- lapply(list_Y, attr, "zeta")
  names(list_zeta) <- names(list_days)
  Y <- do.call(rbind, lapply(list_Y, function(Y_i) { attr(Y_i, "zeta") <- NULL; Y_i }))
  attr(Y, "zeta") <- list_zeta
  Y
}


## ---------------------------------------------------------------------------
## Biomarker panel
## ---------------------------------------------------------------------------

all_var <- cal_panel$var_names
metabolite_names <- cal_panel$metabolite_names

Y_panel <- simulate_values(cal_panel, list_days_panel)
list_zeta_panel <- attr(Y_panel, "zeta")
attr(Y_panel, "zeta") <- NULL
colnames(Y_panel) <- all_var

df_samples_panel <- do.call(rbind, lapply(names(list_days_panel), function(s) {
  data.frame(sample_id_v2 = paste0(s, "_D", list_days_panel[[s]]),
             subject_id = s,
             days_from_sx_or_swab_imputed = list_days_panel[[s]],
             stringsAsFactors = FALSE)
}))

# healthy controls, from the mean and standard deviation of each marker in the
# real controls, independently across markers
n_hc <- cal_panel$n_hc
Y_hc <- sapply(all_var, function(vv) {
  cal_panel$mean_var[vv] + cal_panel$sd_var[vv] *
    (cal_panel$hc_mean_std[vv] + cal_panel$hc_sd_std[vv] * rnorm(n_hc))
})
Y_hc <- matrix(Y_hc, nrow = n_hc, dimnames = list(NULL, all_var))
df_samples_hc <- data.frame(sample_id_v2 = sprintf("SIMHC%02d", seq_len(n_hc)),
                            subject_id = sprintf("SIMHC%02d", seq_len(n_hc)),
                            days_from_sx_or_swab_imputed = NA,
                            stringsAsFactors = FALSE)

cols_clin <- c("severity", "gender", "age", "hospital_outcome",
               "hospital_outcome_days_from_start_sx_or_first_pos_swab",
               "last_observation_days_from_start_sx_or_first_pos_swab")

df_info <- rbind(cbind(df_samples_panel, df_cohort[df_samples_panel$subject_id, cols_clin]),
                 cbind(df_samples_hc,
                       data.frame(severity = "HC",
                                  gender = sample(c("F", "M"), n_hc, replace = TRUE),
                                  age = round(runif(n_hc, 25, 70)),
                                  hospital_outcome = NA,
                                  hospital_outcome_days_from_start_sx_or_first_pos_swab = NA,
                                  last_observation_days_from_start_sx_or_first_pos_swab = NA)))
rownames(df_info) <- df_info$sample_id_v2

df_comb <- cbind(df_info, rbind(Y_panel, Y_hc))
rownames(df_comb) <- df_comb$sample_id_v2

save(df_comb, df_info, all_var, metabolite_names,
     file = file.path(data_dir, "simulated_biomarker_panel.RData"))

# generating values, for the comparison of Supplementary Material S10
saveRDS(list(B = cal_panel$B, list_Phi = cal_panel$list_Phi, time_g = cal_panel$time_g,
             list_zeta = list_zeta_panel, severity = df_cohort[names(list_zeta_panel), "severity"],
             seed = seed),
        file.path(data_dir, "simulated_biomarker_panel_truth.rds"))
cat("Biomarker panel: ", sum(df_cohort$in_panel), " patients, ", nrow(df_samples_panel),
    " samples, ", sum(df_cohort$dead[df_cohort$in_panel]), " deaths in the window, ",
    n_hc, " control samples\n", sep = "")


## ---------------------------------------------------------------------------
## Transcriptomic data
## ---------------------------------------------------------------------------

genes <- cal_rnaseq$genes

E_expr <- simulate_values(cal_rnaseq, list_days_rnaseq, time_g_Phi = cal_rnaseq$time_g_Phi)
attr(E_expr, "zeta") <- NULL

df_samples_rnaseq <- do.call(rbind, lapply(names(list_days_rnaseq), function(s) {
  data.frame(sampleID = paste0(s, "_D", list_days_rnaseq[[s]]),
             recordID = s,
             days_since = list_days_rnaseq[[s]],
             stringsAsFactors = FALSE)
}))

# controls, independently across genes around the mean and standard deviation
# of the patients
E_hc <- sapply(seq_len(nrow(genes)), function(j) {
  cal_rnaseq$mean_var[j] + cal_rnaseq$sd_var[j] * rnorm(n_hc_rnaseq)
})
df_hc_rnaseq <- data.frame(sampleID = c(sprintf("SIMHC%02d", seq_len(n_hc_rnaseq - n_hc_rnaseq_no_id)),
                                        as.character(900000 + seq_len(n_hc_rnaseq_no_id))),
                           recordID = c(sprintf("SIMHC%02d", seq_len(n_hc_rnaseq - n_hc_rnaseq_no_id)),
                                        rep(NA, n_hc_rnaseq_no_id)),
                           days_since = NA,
                           stringsAsFactors = FALSE)

# genes not expressed in blood: noise around a detection floor below the
# lowest mean of the expressed genes
n_undetected <- cal_rnaseq$n_genes_before_filter - nrow(genes)
n_samples <- nrow(df_samples_rnaseq) + n_hc_rnaseq
detection_floor <- min(cal_rnaseq$mean_var) - 2
E_undetected <- matrix(rnorm(n_undetected * n_samples, detection_floor, 0.3), nrow = n_undetected)

E <- rbind(t(rbind(E_expr, E_hc)), E_undetected)

genes_undetected <- data.frame(GeneID = sprintf("SIMUND%05d", seq_len(n_undetected)),
                               Length = sample(genes$Length, n_undetected, replace = TRUE),
                               Chr = NA, hgnc_symbol = "",
                               description = "simulated transcript not expressed in blood",
                               stringsAsFactors = FALSE)
df_genes <- rbind(genes, genes_undetected)

targets <- rbind(df_samples_rnaseq, df_hc_rnaseq)
bool_patient <- !is.na(targets$days_since)
idx_cohort <- match(targets$recordID, df_cohort$subject_id) # NA for the controls
targets$severity_group <- ifelse(bool_patient, df_cohort$severity[idx_cohort], "HC")
targets$gender <- ifelse(bool_patient, df_cohort$gender[idx_cohort],
                         sample(c("F", "M"), nrow(targets), replace = TRUE))
targets$age <- ifelse(bool_patient, df_cohort$age[idx_cohort], round(runif(nrow(targets), 25, 70)))
targets$bmi <- ifelse(bool_patient, df_cohort$bmi[idx_cohort], NA)
targets$hosp_group <- ifelse(bool_patient,
                             c(B = "mild", C = "moderate", D = "severe", E = "severe")[targets$severity_group], NA)
targets$resp <- NA

# batch labels with the real proportions; the sequencing batch is determined by
# the extraction batch, as in the real data
tab_batch <- cal_rnaseq$batches
targets$extractbatch <- sample(rownames(tab_batch), nrow(targets), replace = TRUE, prob = rowSums(tab_batch))
targets$seqbatch <- colnames(tab_batch)[apply(tab_batch, 1, which.max)][match(targets$extractbatch, rownames(tab_batch))]

# hsCRP drawn log-normally with the class-specific mean and standard deviation
# of the real log values
targets$hsCRP_mg_l <- NA
for (cc in classes) {
  idx <- which(targets$severity_group == cc)
  targets$hsCRP_mg_l[idx] <- round(exp(rnorm(length(idx), cal_rnaseq$log_crp[[cc]]["mean"],
                                              cal_rnaseq$log_crp[[cc]]["sd"])), 2)
}

dimnames(E) <- list(df_genes$GeneID, targets$sampleID)

clinical <- df_cohort[df_cohort$in_rnaseq,
                      c("subject_id", "hospital_outcome",
                        "hospital_outcome_days_from_start_sx_or_first_pos_swab",
                        "last_observation_days_from_start_sx_or_first_pos_swab")]
rownames(clinical) <- NULL

elist <- list(E = E,
              weights = matrix(1, nrow = nrow(E), ncol = ncol(E)),
              targets = targets,
              genes = df_genes,
              clinical = clinical)

saveRDS(elist, file.path(data_dir, "simulated_transcriptomics.rds"))
cat("Transcriptomic data: ", sum(df_cohort$in_rnaseq), " patients, ", nrow(df_samples_rnaseq),
    " samples, ", sum(df_cohort$dead[df_cohort$in_rnaseq]), " deaths in the window, ",
    n_hc_rnaseq, " control samples, ", nrow(E), " genes\n", sep = "")
cat("Patients in both data sets: ", sum(df_cohort$in_panel & df_cohort$in_rnaseq), "\n", sep = "")
