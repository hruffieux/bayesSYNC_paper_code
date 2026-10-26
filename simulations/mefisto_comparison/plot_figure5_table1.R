## ---------------------------------------------------------------------------
## Figure 5(b), Supplementary Figure S5 (runtime on the log scale) and the
## numbers of Table 1, from the summaries written by
## summarise_mefisto_comparison.R for every N and both MEFISTO settings.
##
## The bayesSYNC row of the table and the bayesSYNC runtimes are those of the
## replicates paired with subject-specific MEFISTO, except at N = 40, where
## only 14 MEFISTO runs completed and the 25 replicates paired with
## MEFISTO-base are used instead. Figure 5(a) is the average_ROC_curves.pdf
## written by summarise_mefisto_comparison.R at N = 30 for each setting.
##
## Usage, from the root of the repository:
##   Rscript simulations/mefisto_comparison/plot_figure5_table1.R
##
## Writes the runtime plots and metrics_across_N.tex (with a version marking
## the better of the two MEFISTO settings in bold) in
## results_across_N_MEFISTO_with_no_group under out_dir.
## ---------------------------------------------------------------------------

rm(list = ls())

source("R/paths.R")

library(dplyr)
library(ggplot2)
library(tidyr)
library(knitr)
library(kableExtra)

bool_save <- TRUE

############################################################
# COLORS
############################################################

col_bayesSYNC <- "skyblue2"
col_mefisto <- "#ff7f6e"
col_mefisto_base <- adjustcolor("#ff7f6e", alpha.f = 0.35)

plot_dir <- file.path(out_dir, "results_across_N_MEFISTO_with_no_group")

if (bool_save) {
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
}

Ns <- c(20,25,30,35,40)

############################################################
# DIRECTORY DEFINITIONS
############################################################

base_dir_grp_1 <- "summarise_MEFISTO_corr_v2_2026_TRUE_speed_up_FALSE_mean_zeta_0_anneal_1-1.9-100_N_"
base_dir_grp_2 <- "_p_100_Q_true_2_Q_5_L_true_3_L_5_N_min_2_N_max_10_K__sd_zeta_1-0.5-0.33_sd_eps_1_a_om_1_b_om_1_c_0_1_d_0_100_tol_abs_0.001_tol_rel_1e-06_maxit_10000_n_repl_1_nb_seeds_"

base_dir_base_1 <- "summarise_MEFISTO_corr_v2_2026_FALSE_speed_up_FALSE_mean_zeta_0_anneal_1-1.9-100_N_"
base_dir_base_2 <- base_dir_grp_2

n_seeds_grp <- c(25,25,25,25,14)
n_seeds_base <- c(23,24,25,23,25)

res_dirs_grp <- paste0(out_dir, "/",
                       base_dir_grp_1, Ns,
                       base_dir_grp_2, n_seeds_grp)

res_dirs_base <- paste0(out_dir, "/",
                        base_dir_base_1, Ns,
                        base_dir_base_2, n_seeds_base)

############################################################
# LOAD RUNTIME RESULTS
############################################################

runtime_grp <- lapply(res_dirs_grp, function(d){
  load(file.path(d,"runtime_summary.RData"))
  runtime_df
}) %>% bind_rows()

runtime_grp_long <- runtime_grp %>%
  pivot_longer(cols=c("bayesSYNC","MEFISTO"),
               names_to="Method",
               values_to="Runtime")

runtime_base <- lapply(res_dirs_base, function(d){
  load(file.path(d,"runtime_summary.RData"))
  runtime_df
}) %>% bind_rows()

runtime_base_long <- runtime_base %>%
  transmute(
    N,
    Method="MEFISTO-base",
    Runtime=MEFISTO
  )

runtime_long <- bind_rows(runtime_grp_long, runtime_base_long)

############################################################
# Replace BayesSYNC runtime for N=40
############################################################

bayes_base_40 <- runtime_base %>%
  filter(N==40) %>%
  transmute(
    N=40,
    Method="bayesSYNC",
    Runtime=bayesSYNC
  )

runtime_long <- runtime_long %>%
  filter(!(Method=="bayesSYNC" & N==40))

runtime_long <- bind_rows(runtime_long, bayes_base_40)

############################################################
# RUNTIME SUMMARY
############################################################

runtime_summary <- runtime_long %>%
  group_by(N,Method) %>%
  summarise(
    mean=mean(Runtime),
    sd=sd(Runtime),
    .groups="drop"
  )

############################################################
# RUNTIME PLOTS
############################################################

method_colors <- c(
  bayesSYNC=col_bayesSYNC,
  MEFISTO=col_mefisto,
  `MEFISTO-base`=col_mefisto_base
)

p_runtime_linear <- ggplot(runtime_summary,
                           aes(x=N,y=mean,color=Method,group=Method))+
  geom_line(linewidth=1)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),
                width=1,linewidth=0.7)+
  scale_color_manual(values=method_colors)+
  labs(title="Runtime vs number of subjects",
       x="Number of subjects (N)",
       y="Minutes")+
  theme_minimal()+
  theme(panel.grid=element_blank())

if(bool_save){
  ggsave(file.path(plot_dir,"runtime_vs_N_linear.pdf"),
         p_runtime_linear,width=5,height=4)
}

p_runtime_log <- p_runtime_linear +
  scale_y_log10() +
  labs(title="Runtime vs number of subjects (log scale)",
       y="Minutes (log scale)")

if(bool_save){
  ggsave(file.path(plot_dir,"runtime_vs_N_log.pdf"),
         p_runtime_log,width=5,height=4)
}



############################################################
# Base R runtime plots (linear + log, no scientific labels)
############################################################

rt_bayes <- runtime_summary %>% filter(Method == "bayesSYNC")
rt_mefi  <- runtime_summary %>% filter(Method == "MEFISTO")
rt_mefi_b <- runtime_summary %>% filter(Method == "MEFISTO-base")

# Disable scientific notation globally in base R plots
old_scipen <- getOption("scipen")
options(scipen=10)

if (bool_save) {
  pdf(file.path(plot_dir, "runtime_vs_N_linear_base.pdf"), width = 5.3, height = 4.65)
}

ylim_range <- range(c(
  rt_bayes$mean - rt_bayes$sd, rt_bayes$mean + rt_bayes$sd,
  rt_mefi$mean  - rt_mefi$sd, rt_mefi$mean  + rt_mefi$sd,
  rt_mefi_b$mean - rt_mefi_b$sd, rt_mefi_b$mean + rt_mefi_b$sd
))

plot(rt_bayes$N, rt_bayes$mean,
     type = "b", pch = 16, lwd = 2,
     col = col_bayesSYNC,
     ylim = ylim_range,
     xlab = "Number of subjects (N)",
     ylab = "Minutes",
     main = "Runtime vs number of subjects",
     axes = FALSE)

axis(1)
axis(2) # will use plain numbers
box()

arrows(rt_bayes$N,
       rt_bayes$mean - rt_bayes$sd,
       rt_bayes$N,
       rt_bayes$mean + rt_bayes$sd,
       angle = 90, code = 3, length = 0.05,
       col = col_bayesSYNC)

lines(rt_mefi$N, rt_mefi$mean,
      type = "b", pch = 16, lwd = 2,
      col = col_mefisto)

arrows(rt_mefi$N,
       rt_mefi$mean - rt_mefi$sd,
       rt_mefi$N,
       rt_mefi$mean + rt_mefi$sd,
       angle = 90, code = 3, length = 0.05,
       col = col_mefisto)

lines(rt_mefi_b$N, rt_mefi_b$mean,
      type = "b", pch = 16, lwd = 2,
      col = col_mefisto_base)

arrows(rt_mefi_b$N,
       rt_mefi_b$mean - rt_mefi_b$sd,
       rt_mefi_b$N,
       rt_mefi_b$mean + rt_mefi_b$sd,
       angle = 90, code = 3, length = 0.05,
       col = col_mefisto_base)

legend("topleft",
       legend = c("bayesSYNC","MEFISTO","MEFISTO-base"),
       col = c(col_bayesSYNC,col_mefisto,col_mefisto_base),
       lwd = 2, pch = 16, bty = "n")

if (bool_save) dev.off()

############################################################
# Log scale (plain numeric labels)
############################################################

if (bool_save) {
  pdf(file.path(plot_dir, "runtime_vs_N_log_base.pdf"), width = 5.3, height = 4.65)
}

ymin <- min(c(
  rt_bayes$mean - rt_bayes$sd,
  rt_mefi$mean  - rt_mefi$sd,
  rt_mefi_b$mean - rt_mefi_b$sd
))
ymin <- ifelse(ymin <= 0,
               min(c(rt_bayes$mean[rt_bayes$mean>0],
                     rt_mefi$mean[rt_mefi$mean>0],
                     rt_mefi_b$mean[rt_mefi_b$mean>0]))/2,
               ymin)

ymax <- max(c(
  rt_bayes$mean + rt_bayes$sd,
  rt_mefi$mean  + rt_mefi$sd,
  rt_mefi_b$mean + rt_mefi_b$sd
))

plot(rt_bayes$N, rt_bayes$mean,
     type = "b", pch = 16, lwd = 2,
     col = col_bayesSYNC,
     log = "y",
     ylim = c(ymin, ymax),
     xlab = "Number of subjects (N)",
     ylab = "Minutes (log scale)       ",
     main = "Runtime vs number of subjects (log scale)",
     axes = FALSE)

axis(1)
axis(2, las=1) # plain numbers along y-axis
box()

arrows(rt_bayes$N,
       rt_bayes$mean - rt_bayes$sd,
       rt_bayes$N,
       rt_bayes$mean + rt_bayes$sd,
       angle = 90, code = 3, length = 0.05,
       col = col_bayesSYNC)

lines(rt_mefi$N, rt_mefi$mean,
      type = "b", pch = 16, lwd = 2,
      col = col_mefisto)

arrows(rt_mefi$N,
       rt_mefi$mean - rt_mefi$sd,
       rt_mefi$N,
       rt_mefi$mean + rt_mefi$sd,
       angle = 90, code = 3, length = 0.05,
       col = col_mefisto)

lines(rt_mefi_b$N, rt_mefi_b$mean,
      type = "b", pch = 16, lwd = 2,
      col = col_mefisto_base)

arrows(rt_mefi_b$N,
       rt_mefi_b$mean - rt_mefi_b$sd,
       rt_mefi_b$N,
       rt_mefi_b$mean + rt_mefi_b$sd,
       angle = 90, code = 3, length = 0.05,
       col = col_mefisto_base)

legend("bottomright",
       legend = c("bayesSYNC","MEFISTO","MEFISTO-base"),
       col = c(col_bayesSYNC,col_mefisto,col_mefisto_base),
       lwd = 2, pch = 16, bty = "n")

if (bool_save) dev.off()

# Restore original scipen
options(scipen=old_scipen)




############################################################
# LOAD METRICS
############################################################

metrics_grp <- lapply(res_dirs_grp,function(d){
  load(file.path(d,"metrics_summary.RData"))
  metrics_summary
}) %>% bind_rows()

metrics_base <- lapply(res_dirs_base,function(d){
  load(file.path(d,"metrics_summary.RData"))
  metrics_summary
}) %>% bind_rows()

############################################################
# Replace BayesSYNC metrics for N=40 ONLY
############################################################

bayes_metrics_40 <- metrics_base %>%
  filter(N==40) %>%
  select(Metric,N,
         bayesSYNC_mean,
         bayesSYNC_se)

metrics_grp <- metrics_grp %>%
  left_join(bayes_metrics_40,
            by=c("Metric","N"),
            suffix=c("", "_base"))

metrics_grp <- metrics_grp %>%
  mutate(
    bayesSYNC_mean = ifelse(!is.na(bayesSYNC_mean_base),
                            bayesSYNC_mean_base,
                            bayesSYNC_mean),
    bayesSYNC_se   = ifelse(!is.na(bayesSYNC_se_base),
                            bayesSYNC_se_base,
                            bayesSYNC_se)
  ) %>%
  select(-bayesSYNC_mean_base,-bayesSYNC_se_base)

############################################################
# Add MEFISTO-base metrics
############################################################

metrics_base_mefisto <- metrics_base %>%
  mutate(
    MEFISTO_base_mean=MEFISTO_mean,
    MEFISTO_base_se=MEFISTO_se
  ) %>%
  select(Metric,N,MEFISTO_base_mean,MEFISTO_base_se)

metrics_all <- metrics_grp %>%
  left_join(metrics_base_mefisto,by=c("Metric","N"))

############################################################
# LONG FORMAT
############################################################

metrics_long <- metrics_all %>%
  pivot_longer(
    cols=c(
      bayesSYNC_mean,bayesSYNC_se,
      MEFISTO_mean,MEFISTO_se,
      MEFISTO_base_mean,MEFISTO_base_se
    ),
    names_to=c("Method",".value"),
    names_pattern="(MEFISTO_base|bayesSYNC|MEFISTO)_(.*)"
  ) %>%
  mutate(Method=recode(Method,
                       bayesSYNC="bayesSYNC",
                       MEFISTO="MEFISTO",
                       MEFISTO_base="MEFISTO-base"))

metrics_long <- metrics_long %>%
  mutate(Value=sprintf("%.3f (%.3f)",mean,se))

############################################################
# METRICS TABLE
############################################################

table_df <- metrics_long %>%
  select(Metric,Method,N,Value) %>%
  pivot_wider(names_from=N,values_from=Value)

latex_metrics <- kable(
  table_df,
  format="latex",
  booktabs=TRUE,
  caption="Comparison of BayesSYNC, MEFISTO and MEFISTO-base across numbers of subjects. Values are mean (SE)."
)

writeLines(
  paste0(
    "\\documentclass{article}
\\usepackage{booktabs}
\\usepackage[margin=1in]{geometry}
\\begin{document}
\\begin{center}
\\section*{Comparison Across Metrics}
",
    latex_metrics,
    "
\\end{center}
\\end{document}"
  ),
  file.path(plot_dir,"metrics_across_N.tex")
)

############################################################
# BOLD TABLE
############################################################

metrics_bold <- metrics_long %>%
  group_by(Metric,N) %>%
  mutate(best = ifelse(Method!="bayesSYNC",
                       mean==min(mean),
                       FALSE)) %>%
  mutate(Value = ifelse(best,
                        paste0("\\textbf{",Value,"}"),
                        Value))

table_df_bold <- metrics_bold %>%
  select(Metric,Method,N,Value) %>%
  pivot_wider(names_from=N,values_from=Value)

latex_metrics_bold <- kable(
  table_df_bold,
  format="latex",
  escape=FALSE,
  booktabs=TRUE,
  caption="Best result between MEFISTO and MEFISTO-base in bold."
)

writeLines(
  paste0(
    "\\documentclass{article}
\\usepackage{booktabs}
\\usepackage[margin=1in]{geometry}
\\begin{document}
\\begin{center}
\\section*{Comparison Across Metrics (Best in Bold)}
",
    latex_metrics_bold,
    "
\\end{center}
\\end{document}"
  ),
  file.path(plot_dir,"metrics_across_N_bold.tex")
)