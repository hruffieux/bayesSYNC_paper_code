## ---------------------------------------------------------------------------
## Summary over the replicates of the simulated biomarker panel run by
## run_simulated_replicates.sh (Supplementary Material S10): how often each
## simulated factor is recovered, and the accuracy of the recovered factors.
##
## Usage, from the root of the repository:
##   Rscript application/biomarker_panel/summarise_simulated_replicates.R
## ---------------------------------------------------------------------------

source("R/paths.R")

root <- file.path(out_dir, "simulated_replicates")
files <- Sys.glob(file.path(root, "seed_*", "output", "*", "compare_with_truth.csv"))
stopifnot(length(files) > 0)

df <- do.call(rbind, lapply(files, function(f) {
  d <- read.csv(f)
  d$seed <- as.integer(sub(".*seed_([0-9]+)/output.*", "\\1", f))
  d
}))
df <- df[order(df$seed, df$factor), ]

cat("Replicates:", length(unique(df$seed)), "\n")
cat("Number of retained factors:\n")
print(table(df$n_retained[!duplicated(df$seed)]))
cat("\nBy simulated factor, over the replicates in which it is recovered (median, range):\n")
for (k in sort(unique(df$factor))) {
  d <- df[df$factor == k, ]
  cat(sprintf("Factor %d: recovered in %d of %d\n", k, sum(d$recovered), nrow(d)))
  d <- d[d$recovered, ]
  for (v in c("cor_loadings", "cos_phi_1", "cor_scores_1", "rho_sev_true_1", "rho_sev_est_1")) {
    if (nrow(d) > 0) cat(sprintf("  %-15s %.2f (%.2f to %.2f)\n", v, median(d[[v]]), min(d[[v]]), max(d[[v]])))
  }
}
write.csv(df, file.path(root, "summary_simulated_replicates.csv"), row.names = FALSE)
