## ---------------------------------------------------------------------------
## Joint association of the transcriptomic factors with severity, for
## Supplementary Material S9, computed from the output_slim.RData of an
## analysis_transcriptomics.R run.
##
## Several factors are associated with the severity class one at a time
## (Table 3). To see whether they carry distinct information, the severity class
## (B < C < D < E) is regressed on the leading scores of all factors together,
## in a proportional-odds model, since the class is ordinal. The leading score
## of each factor is used, fixed in advance, rather than a selection of
## components. Scores are standardised, so each odds ratio is per standard
## deviation of the score, and reported in the direction of increasing severity
## since the sign of a score is arbitrary. The scores are estimates, shrunk more
## for patients with fewer samples, and their uncertainty is not carried into
## the model, so the odds ratios are best read as descriptive.
##
## Uses MASS (polr), a recommended package shipped with R.
##
## Usage, from the root of the repository:
##   Rscript application/transcriptomics/severity_regression_transcriptomics.R <run_dir>
## Writes severity_regression.txt into <run_dir>, and prints the same.
## ---------------------------------------------------------------------------

rm(list = ls())

suppressPackageStartupMessages(library(MASS))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1, dir.exists(args[1]))
run_dir <- args[1]

e <- new.env()
load(file.path(run_dir, "output_slim.RData"), envir = e)
res <- e$res_slim
Q <- e$Q

subj <- names(res$list_h_hat)
df_first <- e$df_comb_no_HC[order(e$df_comb_no_HC$days_from_sx_or_swab_imputed), ]
df_first <- df_first[!duplicated(df_first$subject_id), ]
vec_sev <- as.character(df_first$severity[match(subj, df_first$subject_id)])
stopifnot(!any(is.na(vec_sev)))

Z <- sapply(1:Q, function(q) res$list_Zeta_hat[[q]][subj, 1])
colnames(Z) <- paste0("F", 1:Q)
Z <- scale(Z)
df <- data.frame(sev = factor(vec_sev, levels = c("B", "C", "D", "E"), ordered = TRUE), Z)

out_file <- file.path(run_dir, "severity_regression.txt")
sink(out_file, split = TRUE)

cat("Patients:", nrow(df), "\n")
C <- cor(Z)
cat("Largest absolute correlation between the leading scores of two factors:",
    round(max(abs(C[upper.tri(C)])), 2), "\n\n")

m_all <- polr(as.formula(paste("sev ~", paste(colnames(Z), collapse = " + "))), df, Hess = TRUE)
s <- coef(summary(m_all))[colnames(Z), ]
p_wald <- 2 * pnorm(-abs(s[, "t value"]))
tab <- data.frame(factor = colnames(Z),
                  OR_per_SD = exp(abs(s[, "Value"])),
                  lower = exp(abs(s[, "Value"]) - 1.96 * s[, "Std. Error"]),
                  upper = exp(abs(s[, "Value"]) + 1.96 * s[, "Std. Error"]),
                  p = p_wald,
                  bonferroni = p_wald < 0.05 / Q)
cat("Proportional-odds model on the leading scores of all factors\n")
cat("(odds ratio per SD, in the direction of increasing severity, with 95% Wald interval)\n")
print(format(tab, digits = 2), row.names = FALSE)

# do the other factors add to the one most associated with severity?
q_top <- which.max(abs(s[, "Value"]))
m_top <- polr(as.formula(paste("sev ~", colnames(Z)[q_top])), df, Hess = TRUE)
lr <- as.numeric(2 * (logLik(m_all) - logLik(m_top)))
cat(sprintf("\nLikelihood-ratio test of the other %d factors added to %s: chi2 = %.1f, p = %.2g\n",
            Q - 1, colnames(Z)[q_top], lr, pchisq(lr, Q - 1, lower.tail = FALSE)))

sink()
cat("written:", out_file, "\n")
