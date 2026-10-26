## ---------------------------------------------------------------------------
## Agreement between the severity axes of the two COVID-19 analyses, over the
## patients common to both: the first-component scores of the inflammatory
## factor of the biomarker panel analysis (Figure 6(c)) and the patient-level
## value of the neutrophil and lymphocyte factor of the transcriptomic analysis
## (mean of the estimated factor trajectory at the patient's own sampling
## times, as in table3_transcriptomic_factors.R).
##
## Both quantities are associated with severity class, so their raw correlation
## is expected to go partly through it. The partial rank correlation given the
## class, and given the class and the mean sampling day in each analysis, says
## how much of the agreement remains beyond severity. Samples of the two
## analyses are not taken on the same days.
##
## On the real data the factors are factor 5 of the N = 88 panel fit (CRP loads
## on it) and factor 1 of the N = 164 transcriptomic fit. On a fit to the
## simulated data the factor numbers differ, and the factors are chosen as the
## retained panel factor with the largest CRP loading and the transcriptomic
## factor whose patient-level value is most associated with severity class.
## Each quantity is oriented so that it increases with severity, since the
## signs are not determined by the data.
##
## Usage, from the root of the repository:
##   Rscript application/transcriptomics/cross_analysis_severity.R <panel output.RData> <transcriptomic run_dir>
## ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2, file.exists(args[1]), dir.exists(args[2]))

m <- new.env()
load(args[1], envir = m)
r <- new.env()
load(file.path(args[2], "output_slim.RData"), envir = r)

bool_simulated <- isTRUE(m$bool_simulated_data) || isTRUE(r$bool_simulated_data)
if (!bool_simulated) {
  q_panel <- 5
  q_rnaseq <- 1
}

# mean of the estimated factor trajectory at the patient's own sampling times
level_at_own <- function(res, time_obs, q) {
  sapply(names(res$list_h_hat), function(s) {
    h <- res$list_h_hat[[s]][[q]]
    mean(sapply(time_obs[[s]], function(t) h[which.min(abs(res$time_g - t))]))
  })
}

# severity class from the first sample, coded 1 to 4 for classes B to E
first_sev <- function(df, ids) {
  df <- df[order(df$days_from_sx_or_swab_imputed), ]
  df <- df[!duplicated(df$subject_id), ]
  match(as.character(df$severity[match(ids, df$subject_id)]), c("B", "C", "D", "E"))
}

mean_day <- function(time_obs, days_thres) sapply(time_obs, function(t) mean(t) * (days_thres - 1))

orient <- function(x, sev) x * sign(cor(x, sev, method = "spearman"))

if (bool_simulated) {
  retained_m <- which(m$res$factor_ppi > 0.5)
  q_panel <- retained_m[which.max(abs(m$res$B_hat["log_CRP", retained_m]))]
  ids_tmp <- names(r$res_slim$list_h_hat)
  sev_tmp <- first_sev(r$df_comb_no_HC, ids_tmp)
  rho_tmp <- sapply(seq_len(ncol(r$res_slim$B_hat)), function(q) {
    abs(cor(level_at_own(r$res_slim, r$time_obs[ids_tmp], q)[ids_tmp], sev_tmp, method = "spearman"))
  })
  q_rnaseq <- which.max(rho_tmp)
  cat("Factors compared: panel ", q_panel, ", transcriptomic ", q_rnaseq, "\n", sep = "")
}

ids_m <- rownames(m$res$list_Zeta_hat[[q_panel]])
sev_m <- first_sev(m$df_comb_no_HC, ids_m)
z_m <- orient(m$res$list_Zeta_hat[[q_panel]][, 1], sev_m)
names(z_m) <- ids_m

ids_r <- names(r$res_slim$list_h_hat)
sev_r <- first_sev(r$df_comb_no_HC, ids_r)
v_r <- orient(level_at_own(r$res_slim, r$time_obs[ids_r], q_rnaseq)[ids_r], sev_r)

common <- intersect(ids_m, ids_r)
sev_c <- sev_m[match(common, ids_m)]
stopifnot(all(sev_c == sev_r[match(common, ids_r)]))
x <- z_m[common]
y <- v_r[common]

cat("Patients: panel ", length(ids_m), ", transcriptomic ", length(ids_r),
    ", common ", length(common), "\n", sep = "")
print(table(factor(c("B", "C", "D", "E")[sev_c], levels = c("B", "C", "D", "E"))))

# partial rank correlation: ranks regressed on the covariates, residuals correlated
partial <- function(a, b, covs) {
  ra <- residuals(lm(rank(a) ~ ., data = covs))
  rb <- residuals(lm(rank(b) ~ ., data = covs))
  ct <- cor.test(ra, rb)
  c(ct$estimate, p = ct$p.value)
}

ct <- cor.test(x, y, method = "spearman", exact = FALSE)
cv <- data.frame(sev = factor(sev_c))
md <- data.frame(day_panel = mean_day(m$time_obs[common], m$days_thres),
                 day_rnaseq = mean_day(r$time_obs[common], r$days_thres))
p1 <- partial(x, y, cv)
p2 <- partial(x, y, cbind(cv, md))

cat(sprintf("Rank correlation %.2f (p = %.2g)\n", ct$estimate, ct$p.value))
cat(sprintf("Within severity class %.2f (p = %.2g)\n", p1[1], p1[2]))
cat(sprintf("Within class, adjusting for mean sampling days %.2f (p = %.2g)\n", p2[1], p2[2]))
cat(sprintf("Association with class: panel %.2f, transcriptomic %.2f\n",
            cor(x, sev_c, method = "spearman"), cor(y, sev_c, method = "spearman")))
cat("Median absolute difference in mean sampling day:", median(abs(md$day_panel - md$day_rnaseq)), "days\n")
for (k in sort(unique(sev_c))) {
  ok <- sev_c == k
  ctk <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  cat(sprintf("  class %s (n = %d): %.2f (p = %.2g)\n", c("B", "C", "D", "E")[k], sum(ok),
              ctk$estimate, ctk$p.value))
}
