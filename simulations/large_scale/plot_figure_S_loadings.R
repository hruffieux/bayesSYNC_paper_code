## ---------------------------------------------------------------------------
## Supplementary Figure S1 (Supplementary Material S3): simulated and estimated
## loadings of the single replicate of Section 5.2, drawn in the style of
## Figure 2(b), with one framed panel per matrix, variables numbered 1 to 50 on
## the side, a shared colour scale with its legend, and the signs of the
## nonzero simulated loadings marked on the estimated panel.
##
## Reads the output.RData of the run used for Figure 2 (replicate 1, Q_max = 5,
## annealed VB; job_id = 151 of run_large_scale.R), which holds
## list_res[[1]]$B (simulated) and list_res[[1]]$B_hat (estimated, already
## matched to the simulated factors).
##
## Usage, from the root of the repository:
##   Rscript simulations/large_scale/plot_figure_S_loadings.R [output.RData] [out.pdf]
## Without arguments, reads the fit from out_dir and writes
## figure_S_loadings_single_repl.pdf in out_dir.
## ---------------------------------------------------------------------------

source("R/paths.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  args <- c(file.path(out_dir, paste0("BS_r_anneal_1-1.9-100_N_100_p_20000_Q_true_3_Q_5_L_true_2_L_5",
                                      "_N_min_5_N_max_10_K__sd_zeta_1-0.5_sd_eps_1_a_om_1_b_om_10",
                                      "_tol_abs_0.001_tol_rel_1e-06_maxit_10000_n_cpus_50_n_repl_1_seed_1"),
                      "output.RData"),
            file.path(out_dir, "figure_S_loadings_single_repl.pdf"))
}
stopifnot(length(args) == 2, file.exists(args[1]))

e <- new.env()
load(args[1], envir = e)
B <- e$list_res[[1]]$B
B_hat <- e$list_res[[1]]$B_hat

trunc <- 50 # variables shown, as in Figure 2(b)
B <- B[seq_len(trunc), , drop = FALSE]
B_hat <- B_hat[seq_len(trunc), , drop = FALSE]
stopifnot(ncol(B) == ncol(B_hat))
Q <- ncol(B)


## colour scale shared by the two panels, as in Figure 2(b): blue for negative
## loadings, red for positive ones and a white band around zero, so that the
## many exact zeros of the simulated matrix and the near-zero estimates both
## read as white

n_col <- 50
n_white <- round(0.25 * n_col)
n_side <- (n_col - n_white) / 2
pal <- c(colorRampPalette(c("skyblue1", "white"))(n_side),
         rep("white", n_white),
         colorRampPalette(c("white", "#ff7f6e"))(n_side))

all_val <- c(B, B_hat)
lim_lo <- min(all_val)
lim_hi <- max(all_val)
# width of the white band, taken over the nonzero entries because most
# loadings are exactly zero
eps <- quantile(abs(all_val)[all_val != 0], 0.1)
breaks <- c(seq(lim_lo, -eps, length.out = n_side + 1),
            seq(-eps, eps, length.out = n_white + 1)[-1],
            seq(eps, lim_hi, length.out = n_side + 1)[-1])


## one framed panel; variable 1 at the top, as in Figure 2(b)

draw_panel <- function(M, title, signs = NULL) {
  image(x = seq_len(Q), y = seq_len(trunc), z = t(M[trunc:1, , drop = FALSE]),
        col = pal, breaks = breaks, axes = FALSE, xlab = "", ylab = "")
  if (!is.null(signs)) {
    ind <- which(signs != 0, arr.ind = TRUE)
    text(ind[, 2], trunc + 1 - ind[, 1], ifelse(signs[ind] > 0, "+", "\u2212"), cex = 0.9)
  }
  box(lwd = 1.2)
  axis(1, at = seq_len(Q), labels = seq_len(Q), tick = FALSE, line = -0.6, cex.axis = 0.85)
  mtext("factors", side = 1, line = 1.6, cex = 0.75)
  axis(4, at = trunc:1, labels = seq_len(trunc), las = 1, tick = FALSE, line = -0.7,
       cex.axis = 0.55)
  mtext(title, side = 3, line = 0.5, cex = 0.85, font = 2)
}


cairo_pdf(args[2], width = 5.6, height = 7.2)
layout(matrix(1:3, nrow = 1), widths = c(1, 1, 0.75))

par(mar = c(3, 0.8, 2.5, 2.2), oma = c(0, 2.2, 0, 0))
draw_panel(B, "Simulated")
mtext("variables (first 50 out of 20,000)", side = 2, line = 0.8, cex = 0.7)
draw_panel(B_hat, "Estimated", signs = sign(B))

## legends: colour scale and signs, as in Figure 2(b)
par(mar = c(3, 1.2, 2.5, 3.2))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
bar_lo <- 0.08; bar_hi <- 0.48
y_br <- seq(bar_lo, bar_hi, length.out = length(pal) + 1)
rect(0.05, y_br[-length(y_br)], 0.3, y_br[-1], col = pal, border = NA)
rect(0.05, bar_lo, 0.3, bar_hi, border = "black", lwd = 0.8)
ticks <- pretty(c(lim_lo, lim_hi), n = 5)
ticks <- ticks[ticks >= lim_lo & ticks <= lim_hi]
# position of each tick on the bar, through the same breaks as the colours
pos <- approx(breaks, seq(bar_lo, bar_hi, length.out = length(breaks)), xout = ticks)$y
segments(0.3, pos, 0.36, pos, lwd = 0.8)
text(0.4, pos, sub("-", "\u2212", format(ticks)), adj = 0, cex = 0.75)
text(0.02, bar_hi + 0.06, "Loadings", adj = 0, cex = 0.8, font = 2)

text(0.02, 0.78, "Sign of\nsimulated\nloadings", adj = c(0, 0.5), cex = 0.8, font = 2)
text(0.12, 0.70, "+", cex = 0.9)
text(0.12, 0.67, "\u2212", cex = 0.9)

invisible(dev.off())
cat("written:", args[2], "\n")
