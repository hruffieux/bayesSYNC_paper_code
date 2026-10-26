## ---------------------------------------------------------------------------
## Match the factors of two fitted runs and report how well they correspond.
##
## Written for the mortality sensitivity analysis, where the same data are fitted
## with and without the subjects whose follow-up is truncated by death, but it
## applies to any pair of runs on the same variables, for instance two values of
## K or two seeds.
##
## Factor labels are not comparable across runs, so the factors have to be
## matched before anything can be said about consistency. Two statistics are
## computed on the common set of variables.
##
##   Loading correlation. The absolute Pearson correlation between the loading
##     columns, which uses every variable and is invariant to the sign flip that
##     the model leaves undetermined. This is the primary criterion.
##   Restricted correlation. The same correlation computed only over the union
##     of the two top-N sets. The loadings are sparse against a dense background
##     of near-zero values, so a correlation over all p variables is diluted by
##     the variables that carry no signal in either factor, and understates a
##     correspondence that is in fact clear. This is usually the most readable
##     of the three.
##   Top-N overlap. The overlap between the N variables with the largest
##     absolute loadings in each factor. This is the complement that matters
##     here, because the posterior inclusion probabilities select most of the
##     panel and a comparison based on the selected sets is dominated by their
##     size rather than by which variables carry the signal.
##
## Matching is greedy on the correlation matrix, which is adequate when the
## correspondence is clear. The full matrix is printed and plotted so that an
## ambiguous case is visible rather than silently resolved, and the margin to
## the runner-up is reported for every matched pair.
##
## Usage, from the root of the repository:
##   Rscript application/biomarker_panel/compare_runs.R <run_dir_a> <run_dir_b> [label_a] [label_b]
## The correspondence figure is written in run_dir_b.
## ---------------------------------------------------------------------------

rm(list = ls())

## ---------------------------------------------------------------------------
## The two runs to compare. Give the directories written by the application
## script; each must hold output_slim.RData or output.RData.
## ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2)
dir_a <- args[1]   # for example the run including the subjects who died
dir_b <- args[2]   # and the run excluding them
lab_a <- if (length(args) >= 3) args[3] else "run A"
lab_b <- if (length(args) >= 4) args[4] else "run B"

n_top <- 200       # variables per factor in the top-N comparison
thres_ppi <- 0.5   # only used for the reported selection counts
# A match is flagged as weak on the restricted correlation rather than the
# correlation over all variables, which is diluted by the dense near-zero
# background and can sit near 0.4 for a correspondence that is not in doubt.
min_cor <- 0.5     # applied to cor_top, below
min_margin <- 0.1  # and a match this close to the runner-up is flagged as ambiguous

res_dir_cmp <- dir_b # where the correspondence figure is written


## ---------------------------------------------------------------------------
## Loading
## ---------------------------------------------------------------------------

load_fit <- function(dir_run) {

  stopifnot(nzchar(dir_run))
  file_slim <- file.path(dir_run, "output_slim.RData")
  file_full <- file.path(dir_run, "output.RData")
  file_run <- if (file.exists(file_slim)) file_slim else file_full
  if (!file.exists(file_run)) {
    stop("No output_slim.RData or output.RData in ", dir_run)
  }

  e <- new.env()
  load(file_run, envir = e)

  # the slim file holds res_slim, the full one res; both keep B_hat and ppi
  nm <- intersect(c("res_slim", "res"), ls(e))
  if (length(nm) == 0) stop("Neither res_slim nor res found in ", file_run)
  fit <- get(nm[1], envir = e)

  stopifnot(!is.null(fit$B_hat), !is.null(rownames(fit$B_hat)))
  list(B_hat = fit$B_hat, ppi = fit$ppi, factor_ppi = fit$factor_ppi,
       file = file_run)
}

fit_a <- load_fit(dir_a)
fit_b <- load_fit(dir_b)

vec_common <- intersect(rownames(fit_a$B_hat), rownames(fit_b$B_hat))
cat("Variables: ", nrow(fit_a$B_hat), " in ", lab_a, ", ", nrow(fit_b$B_hat),
    " in ", lab_b, ", ", length(vec_common), " in common\n", sep = "")
stopifnot(length(vec_common) > 0)

B_a <- fit_a$B_hat[vec_common, , drop = FALSE]
B_b <- fit_b$B_hat[vec_common, , drop = FALSE]
Q_a <- ncol(B_a)
Q_b <- ncol(B_b)


## ---------------------------------------------------------------------------
## Correspondence
## ---------------------------------------------------------------------------

# a factor whose loadings are constant has no correlation with anything; keep
# the entry as NA rather than letting cor() warn for every pair
safe_cor <- function(x, y) {
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(x, y)
}

mat_cor <- outer(seq_len(Q_a), seq_len(Q_b),
                 Vectorize(function(i, j) safe_cor(B_a[, i], B_b[, j])))
dimnames(mat_cor) <- list(paste0(lab_a, " q", seq_len(Q_a)),
                          paste0(lab_b, " q", seq_len(Q_b)))

top_set <- function(B, q, n) {
  v <- abs(B[, q])
  v <- v[v > 0]
  names(sort(v, decreasing = TRUE))[seq_len(min(n, length(v)))]
}

list_top_a <- lapply(seq_len(Q_a), top_set, B = B_a, n = n_top)
list_top_b <- lapply(seq_len(Q_b), top_set, B = B_b, n = n_top)

# correlation restricted to the union of the two top-N sets
mat_cor_top <- outer(seq_len(Q_a), seq_len(Q_b),
                     Vectorize(function(i, j) {
                       ind <- union(top_set(B_a, i, n_top), top_set(B_b, j, n_top))
                       if (length(ind) < 3) return(NA_real_)
                       safe_cor(B_a[ind, i], B_b[ind, j])
                     }))
dimnames(mat_cor_top) <- list(paste0(lab_a, " q", seq_len(Q_a)),
                              paste0(lab_b, " q", seq_len(Q_b)))

mat_top <- outer(seq_len(Q_a), seq_len(Q_b),
                 Vectorize(function(i, j) {
                   a <- list_top_a[[i]]; b <- list_top_b[[j]]
                   if (length(a) == 0 || length(b) == 0) return(NA_real_)
                   length(intersect(a, b)) / min(length(a), length(b))
                 }))
dimnames(mat_top) <- dimnames(mat_cor)

cat("\n=== absolute loading correlation, ", lab_a, " (rows) against ", lab_b,
    " (columns) ===\n", sep = "")
print(round(abs(mat_cor), 2))

cat("\n=== absolute correlation over the union of the top-", n_top,
    " loadings ===\n", sep = "")
print(round(abs(mat_cor_top), 2))

cat("\n=== top-", n_top, " overlap ===\n", sep = "")
print(round(mat_top, 2))


## ---------------------------------------------------------------------------
## Greedy matching on the loading correlation
## ---------------------------------------------------------------------------

mat_work <- abs(mat_cor)
mat_work[is.na(mat_work)] <- -Inf

df_match <- NULL
for (k in seq_len(min(Q_a, Q_b))) {

  ind <- which(mat_work == max(mat_work), arr.ind = TRUE)[1, ]
  i <- ind[1]; j <- ind[2]
  if (!is.finite(mat_work[i, j])) break

  # margin to the best alternative for this row, before it is struck out
  vec_row <- mat_work[i, ]
  margin <- max(vec_row) - max(vec_row[-j], -Inf)

  df_match <- rbind(df_match, data.frame(
    q_a = i, q_b = j,
    abs_cor = abs(mat_cor[i, j]),
    cor_top = abs(mat_cor_top[i, j]),
    sign = ifelse(mat_cor[i, j] >= 0, "+", "-"),
    top_overlap = mat_top[i, j],
    margin = margin,
    n_sel_a = if (is.null(fit_a$ppi)) NA else sum(fit_a$ppi[, i] > thres_ppi),
    n_sel_b = if (is.null(fit_b$ppi)) NA else sum(fit_b$ppi[, j] > thres_ppi),
    stringsAsFactors = FALSE))

  mat_work[i, ] <- -Inf
  mat_work[, j] <- -Inf
}

df_match$weak <- df_match$cor_top < min_cor | is.na(df_match$cor_top)
df_match$ambiguous <- df_match$margin < min_margin

cat("\n=== matched factors ===\n")
print(df_match, row.names = FALSE, digits = 3)

cat("\nUnmatched in ", lab_a, ": ",
    paste(setdiff(seq_len(Q_a), df_match$q_a), collapse = ", "), "\n", sep = "")
cat("Unmatched in ", lab_b, ": ",
    paste(setdiff(seq_len(Q_b), df_match$q_b), collapse = ", "), "\n", sep = "")

cat("\nMatches with |correlation| below ", min_cor, ": ", sum(df_match$weak), "\n", sep = "")
cat("Matches within ", min_margin, " of their runner-up: ", sum(df_match$ambiguous),
    "\n", sep = "")
cat("Median absolute correlation over matched pairs, all variables: ",
    round(median(df_match$abs_cor), 3), ", restricted to the top ", n_top, ": ",
    round(median(df_match$cor_top, na.rm = TRUE), 3), "\n", sep = "")


## ---------------------------------------------------------------------------
## Correspondence figure
## ---------------------------------------------------------------------------

pdf(file.path(res_dir_cmp, "factor_correspondence.pdf"), width = 5.5, height = 5)
par(mar = c(5, 5, 3, 2))
m <- abs(mat_cor)
image(seq_len(Q_a), seq_len(Q_b), m, col = grey(seq(1, 0, length.out = 64)),
      zlim = c(0, 1), xlab = paste(lab_a, "factor"), ylab = paste(lab_b, "factor"),
      main = "Absolute loading correlation", axes = FALSE)
axis(1, at = seq_len(Q_a)); axis(2, at = seq_len(Q_b)); box()
points(df_match$q_a, df_match$q_b, pch = 1, col = "#d95f02", cex = 2.4, lwd = 2)
invisible(dev.off())
cat("\nwritten: ", file.path(res_dir_cmp, "factor_correspondence.pdf"), "\n", sep = "")
