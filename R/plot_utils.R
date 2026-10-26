plot_scores <- function(q,
                        l,
                        Zeta,
                        list_Zeta_hat,
                        list_sd_hat = NULL  # optional: posterior SDs of estimated scores
) {
  
  q_res <- q
  l_res_orth <- l
  
  N <- nrow(Zeta[[1]])
  
  ## True scores
  Zeta_ql <- sapply(
    1:N,
    function(i) Zeta[[q_res]][i, l_res_orth]
  )
  
  ## Estimated scores
  Zeta_hat_ql <- sapply(
    1:N,
    function(i) list_Zeta_hat[[q_res]][i, l_res_orth]
  )
  
  par(pty = "s")
  
  plot(
    Zeta_ql, Zeta_hat_ql,
    xlab = "Simulated",
    ylab = "Estimated",
    pch = 20,
    main = paste0(
      "Scores\nFactor q = ", q,
      ", FPC l = ", l
    )
  )
  
  ## Add 95% credible intervals if SDs are provided
  if (!is.null(list_sd_hat)) {
    
    sd_hat_ql <- sapply(
      1:N,
      function(i) list_sd_hat[[q_res]][i, l_res_orth]
    )
    
    y_low <- Zeta_hat_ql - 1.96 * sd_hat_ql
    y_upp <- Zeta_hat_ql + 1.96 * sd_hat_ql
    
    segments(
      x0 = Zeta_ql,
      y0 = y_low,
      x1 = Zeta_ql,
      y1 = y_upp,
      lwd = 1
    )
  }
  
  ## Optional reference
  # abline(a = 0, b = 1, lty = 2)
}



library(ggplot2)

plot_bar_with_error_ggplot <- function(mat, l_values = NULL, ylab = "", main_title = "", bar_colors = c("skyblue1", "honeydew3")) {
  # Convert the matrix to a long-format data frame
  mat_long <- as.data.frame(mat)
  mat_long$l <- factor(rep(l_values, each = nrow(mat)))  # Add l-values if provided
  mat_long$Sample <- factor(rep(1:nrow(mat), times = ncol(mat)))  # Sample identifiers
  mat_long <- reshape2::melt(mat_long, id.vars = c("Sample", "l"), variable.name = "Condition", value.name = "Value")
  
  # Calculate means and standard errors by Condition and l
  summary_df <- mat_long %>%
    group_by(Condition, l) %>%
    summarise(
      Mean = mean(Value),
      SE = sd(Value) / sqrt(n()),
      .groups = "drop"
    )
  
  # Create the ggplot
  p <- ggplot(summary_df, aes(x = Condition, y = Mean, fill = Condition)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "black") +
    geom_errorbar(
      aes(ymin = Mean - SE, ymax = Mean + SE),
      width = 0.2,
      position = position_dodge(0.9)
    ) +
    scale_fill_manual(values = bar_colors) +
    facet_wrap(~l, scales = "free_y") +
    labs(title = main_title, y = ylab, x = "Condition") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}


# Function to create a barplot with error bars
plot_bar_with_error <- function(mat, 
                                bar_colors = c("skyblue1", "honeydew3"), 
                                ylab = "", 
                                main_title = "", 
                                error_bar_length = 0.2, 
                                error_bar_lwd = 2,
                                ylim_max = NULL,
                                log_scale = FALSE) {
  # Calculate means and standard errors
  means <- colMeans(mat, na.rm = T)
  se <- apply(mat, 2, function(x) sd(x, na.rm = T) / sqrt(sum(!is.na(x))))
  
  # Define plot limits
  # Define plot limits
  if (is.null(ylim_max)) {
    ylim_max <- if (log_scale) {
      1.2 * max(means * exp(se)) # Use exp(se) for log-transformed upper bounds
    } else {
      1.2 * max(means + se)
    }
  }
  
  # If log_scale is TRUE, ensure all values are >0 (log scale cannot handle zeros or negatives)
  if (log_scale && any(means <= 0)) {
    stop("Means must be greater than zero for log scale plotting.")
  }
  
  # Plot the barplot
  ze_barplot <- barplot(
    means, 
    beside = TRUE, 
    col = bar_colors, 
    log = if (log_scale) "y" else "", 
    ylim = c(if (log_scale) 0.1 else 0, ylim_max), 
    names.arg = colnames(mat), 
    ylab = ylab, 
    main = main_title
  )
  
  # Function to add error bars
  error.bar <- function(x, y, upper, lower = upper, length = 0.1, lwd = 2, ...) {
    if (log_scale) {
      arrows(x, y * exp(upper), x, y / exp(lower), angle = 90, code = 3, length = length, lwd = lwd, ...)
    } else {
      arrows(x, y + upper, x, y - lower, angle = 90, code = 3, length = length, lwd = lwd, ...)
    }
  }
  
  # Add the error bars
  error.bar(ze_barplot, means, se, length = error_bar_length, lwd = error_bar_lwd)
}


plot_factor_trajectory <- function(list_Zeta_hat, q, l, time_g, h_hat, h_low, h_upp, days_thres, df_comb, offset = 0.1,
                                   bool_rescale_time = T) {
  
  fac <- ifelse(bool_rescale_time, days_thres-1, 1) # scaling factor back to the original time window
  xlim_max <- ifelse(bool_rescale_time, days_thres, 1)
  
  scores_q <- list_Zeta_hat[[q]]
  vec_extreme_scores_q <- c(min(scores_q[,l]), max(scores_q[,l]))
  vec_subj_extreme_q <- c(rownames(scores_q)[which.min(scores_q[,l])],
                          rownames(scores_q)[which.max(scores_q[,l])]) # extreme subjects
  
  list_ylim <- c(min(c(h_low[[vec_subj_extreme_q[1]]][[q]], h_low[[vec_subj_extreme_q[2]]][[q]],
                       h_upp[[vec_subj_extreme_q[1]]][[q]], h_upp[[vec_subj_extreme_q[2]]][[q]]))-offset,
            max(c(h_low[[vec_subj_extreme_q[1]]][[q]], h_low[[vec_subj_extreme_q[2]]][[q]],
                  h_upp[[vec_subj_extreme_q[1]]][[q]], h_upp[[vec_subj_extreme_q[2]]][[q]]))+offset)
  
  for (id in 1:2) {
    
    score <- vec_extreme_scores_q[id]
    subj <- vec_subj_extreme_q[id]
    
    df_comb$time <- df_comb$hospital_outcome_days_from_start_sx_or_first_pos_swab
    df_comb$time[!(df_comb$hospital_outcome %in% "dead")] <- df_comb$last_observation_days_from_start_sx_or_first_pos_swab[!(df_comb$hospital_outcome %in% "dead")]
    df_comb$hospital_outcome[df_comb$severity %in% c("HC", "A", "B")] <- "alive"
    df_comb$status <- ifelse(df_comb$hospital_outcome == "dead", 2, 1) # 1 = cencored, 2 = dead
    
    status <- unique(df_comb$status[df_comb$subject_id %in% subj])
    time <- unique(df_comb$time[df_comb$subject_id %in% subj])
    time <- time[!is.na(time)]
    
    stopifnot(length(status) == 1 & length(time) == 1)
    
    time_grid_subj_id <- time_g
    if (bool_rescale_time) {
      time_grid_subj_id <- time_grid_subj_id * fac
    }
    
    
    if (status == 2 & time < days_thres) { # cut prediction if patient died within the analysis window
      
      plot(time_grid_subj_id[time_grid_subj_id<=time], h_hat[[subj]][[q]][time_grid_subj_id<=time], 
           type = "l", xlab = "Days from symptom onset", ylab=paste0("q = ", q), 
           main = paste0(subj, ", FPC ", l, " = ", format(score, digits = 2), "\n",  
                         "(", unique(df_comb$severity[df_comb$subject_id %in% subj]), 
                         ", gender: ", unique(df_comb$gender[df_comb$subject_id %in% subj]), 
                         ", age: ", unique(df_comb$age[df_comb$subject_id %in% subj]), 
                         ")"), ylim = list_ylim, 
           xlim = c(0, xlim_max),
            col="black",lwd=2)
      points(time_grid_subj_id[which.min(abs(time_grid_subj_id-time))], 
             h_hat[[subj]][[q]][which.min(abs(time_grid_subj_id-time))], 
             pch = 4, cex = 1.25, lwd = 1.5)
      if (!is.null(h_upp)){
        lines(time_grid_subj_id[time_grid_subj_id<=time], 
              h_upp[[subj]][[q]][time_grid_subj_id<=time],
              col="black",lwd=2,lty=2)
      }
      if (!is.null(h_low)){
        lines(time_grid_subj_id[time_grid_subj_id<=time], 
              h_low[[subj]][[q]][time_grid_subj_id<=time],
              col="black",lwd=2,lty=2)
      }
    } else {
      plot(time_grid_subj_id, h_hat[[subj]][[q]], 
           type = "l", xlab = "Days from symptom onset", ylab=paste0("q = ", q), 
           main = paste0(subj, ", FPC ", l, " = ", format(score, digits = 2), "\n",  
                         "(", unique(df_comb$severity[df_comb$subject_id %in% subj]), 
                         ", gender: ", unique(df_comb$gender[df_comb$subject_id %in% subj]), 
                         ", age: ", unique(df_comb$age[df_comb$subject_id %in% subj]), 
                         ")"), ylim = list_ylim, 
           xlim = c(0, xlim_max),
           col="black",lwd=2)
      if (!is.null(h_upp)){
        lines(time_grid_subj_id, 
              h_upp[[subj]][[q]],
              col="black",lwd=2,lty=2)
      }
      if (!is.null(h_low)){
        lines(time_grid_subj_id, 
              h_low[[subj]][[q]],
              col="black",lwd=2,lty=2)
      }
    }
    abline(h = 0, lty = 3)
  }
  
}


plot_factor_variables_trajectory <- function(list_Zeta_hat, ppi, q, time_obs, 
                                             time_g, Y, Y_hat, Y_low, Y_upp, 
                                             days_thres, df_comb, 
                                             vec_variables = NULL, thres_ppi = 0.5, offset = 0.1,
                                             show_ppi = T) {
  
  scores_q <- res$list_Zeta_hat[[q]]
  
  if (is.null(vec_variables)) {
    vec_variables <- rownames(ppi)[ppi[,q]>thres_ppi]
    if (length(vec_variables)==0) {
      vec_variables <- rownames(ppi)[1:min(4, nrow(ppi))]
    }
  }
  
  fac <- ifelse(!is.null(days_thres), days_thres-1, 1) # scaling factor back to the original time window
  xlim_max <- ifelse(!is.null(days_thres), days_thres, 1)
  
  par(mfrow = c(min(length(vec_variables), 5), 4))
  for (vv in vec_variables) {
    
    for (l in 1:2) {
    
      vec_extreme_scores_q <- c(min(scores_q[,l]), max(scores_q[,l]))
      vec_subj_extreme_q <- c(rownames(scores_q)[which.min(scores_q[,l])],
                              rownames(scores_q)[which.max(scores_q[,l])]) # extreme subjects
      
      list_ylim <- c(min(c(Y[[vec_subj_extreme_q[1]]][[vv]], Y[[vec_subj_extreme_q[1]]][[vv]],
                           Y_low[[vec_subj_extreme_q[1]]][[vv]], Y_low[[vec_subj_extreme_q[2]]][[vv]]))-offset,
                     max(c(Y[[vec_subj_extreme_q[1]]][[vv]], Y[[vec_subj_extreme_q[1]]][[vv]],
                           Y_upp[[vec_subj_extreme_q[1]]][[vv]], Y_upp[[vec_subj_extreme_q[2]]][[vv]]))+offset)
      
      if (any(df_comb$severity == "HC")) {
        lo <- quantile(df_comb[df_comb$severity == "HC", vv], probs = 0.25, na.rm = TRUE)
        up <- quantile(df_comb[df_comb$severity == "HC", vv], probs = 0.75, na.rm = TRUE)
        
        list_ylim <- c(min(list_ylim[1], lo), max(list_ylim[2], up))
      }
      
      for (id in 1:2) {
        
        score <- vec_extreme_scores_q[id]
        subj <- vec_subj_extreme_q[id]
        
        df_comb$time <- df_comb$hospital_outcome_days_from_start_sx_or_first_pos_swab
        df_comb$time[!(df_comb$hospital_outcome %in% "dead")] <- df_comb$last_observation_days_from_start_sx_or_first_pos_swab[!(df_comb$hospital_outcome %in% "dead")]
        df_comb$hospital_outcome[df_comb$severity %in% c("HC", "A", "B")] <- "alive"
        df_comb$status <- ifelse(df_comb$hospital_outcome == "dead", 2, 1) # 1 = cencored, 2 = dead
        
        status <- unique(df_comb$status[df_comb$subject_id %in% subj])
        time <- unique(df_comb$time[df_comb$subject_id %in% subj])
        time <- time[!is.na(time)]
        
        stopifnot(length(status) == 1 & length(time) == 1)
        
        time_grid_subj_id <- time_g
        if (!is.null(days_thres)) {
          time_grid_subj_id <- time_grid_subj_id * fac
        }
        
        
        plot(fac*time_obs[[subj]], Y[[subj]][[vv]], type = "p", pch = 20, xlab = "Days from symptom onset", 
             main = paste0(subj, ", FPC ", l, " = ", format(score, digits = 2), "\n",
                      "(", unique(df_comb$severity[df_comb$subject_id %in% subj]), 
                      ", gender: ", unique(df_comb$gender[df_comb$subject_id %in% subj]), 
                      ", age: ", unique(df_comb$age[df_comb$subject_id %in% subj]), 
                      ")"),
             ylim = list_ylim, 
             xlim = c(0, xlim_max), 
             # xlim = c(min(time_g), max(time_g)),
             ylab = paste0(ifelse(l == 1 & id == 1, gsub("log_", "", vv), ""), 
                           ifelse(show_ppi, paste0(" (PPI = ", format(ppi[vv,q], digits = 2), ")"), "")),
             col = "grey40")
        
        if (status == 2 & time < days_thres) { # cut prediction if patient died within the analysis window
          lines(time_grid_subj_id[time_grid_subj_id<=time],
                Y_hat[[subj]][[vv]][time_grid_subj_id<=time],
                col="black",lwd=2)
          points(time_grid_subj_id[which.min(abs(time_grid_subj_id-time))], 
                 Y_hat[[subj]][[vv]][which.min(abs(time_grid_subj_id-time))], 
                 pch = 4, cex = 1.25, lwd = 1.5)
          if (!is.null(Y_upp)){
            lines(time_grid_subj_id[time_grid_subj_id<=time], 
                  Y_upp[[subj]][[vv]][time_grid_subj_id<=time],
                  col="black",lwd=2,lty=2)
          }
          if (!is.null(Y_low)){
            lines(time_grid_subj_id[time_grid_subj_id<=time], 
                  Y_low[[subj]][[vv]][time_grid_subj_id<=time],
                  col="black",lwd=2,lty=2)
          }
        } else {
          lines(time_grid_subj_id,
                Y_hat[[subj]][[vv]],
                col="black",lwd=2)
          if (!is.null(Y_upp)){
            lines(time_grid_subj_id, 
                  Y_upp[[subj]][[vv]],
                  col="black",lwd=2,lty=2)
          }
          if (!is.null(Y_low)){
            lines(time_grid_subj_id, 
                  Y_low[[subj]][[vv]],
                  col="black",lwd=2,lty=2)
          }
        }
        
        if (any(df_comb$severity == "HC")) {
          rect(-5, lo, days_thres+5, up, col = adjustcolor("gray", alpha.f=0.5), border = "gray")
        }
        
      }  
    }
    
  }

}


trunc_perf <-  function(perf, fpr.stop) { # perf = list with one entry per replication
  
  for (iperf in seq_along(perf@x.values)){
    ind = which(perf@x.values[[iperf]] <= fpr.stop)
    perf@y.values[[iperf]] = perf@y.values[[iperf]][ind]
    perf@x.values[[iperf]] = perf@x.values[[iperf]][ind]
  }
  
  perf
}

compute_roc_list <- function(list_ppi, list_pat, fpr_stop = 1) {
  list_vec_rank <- lapply(list_ppi, function(ppi) as.numeric(as.factor(ppi)))
  list_vec_pat <- lapply(list_pat, function(pat) as.numeric(pat))
  
  require(ROCR)
  pred <- ROCR::prediction(list_vec_rank, list_vec_pat)
  perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
  
  perf_pauc <- performance(pred, measure = "auc", fpr.stop = fpr_stop)
  perf_pauc_norm <- perf_pauc
  
  if (fpr_stop < 1) {
    perf <- trunc_perf(perf, fpr_stop)
    perf_pauc_norm@y.values = lapply(perf_pauc@y.values, function(x) x / fpr_stop)
  } 
  
  create_named_list(perf, perf_pauc_norm)
}

plot_roc_list <- function(list_ppi, list_pat, leg = NULL, bool_leg = T, col = "black", main = "Average ROC curve", fpr_stop = 1, add = F) {
  
  par(pty = "s")
  
  list_roc_auc <- compute_roc_list(list_ppi, list_pat, fpr_stop = fpr_stop)
  
  output_roc <- list_roc_auc$perf
  output_pauc_norm <- unlist(list_roc_auc$perf_pauc_norm@y.values)

  if (length(list_ppi)>1) {
    if (is.null(leg) & bool_leg) {
      leg <- paste0(ifelse(fpr_stop < 1, "pAUC = ", "AUC = "), format(mean(output_pauc_norm), digits = 2), " (", format(round(sd(output_pauc_norm) / sqrt(length(list_ppi)), 2), nsmall = 2), ")")
    }
    ROCR::plot(output_roc, col = col[1], type = "l", main = main,
    lty = 1, lwd = 2, avg="vertical", spread.estimate="stderror", 
    spread.scale = 1.96, add = add, ylim = c(0, 1))
  } else {
    if (is.null(leg) & bool_leg) {
      leg <- paste0(ifelse(fpr_stop < 1, "pAUC = ", "AUC = "), format(output_pauc_norm, digits = 2))
    }
    ROCR::plot(output_roc, main = main, col = col, type = "l", lty = 1, lwd = 2, add = add, ylim = c(0, 1))
  }
  
  if (fpr_stop == 1) {
    abline(0,1)
  }
  if (bool_leg) {
    legend("bottomright", leg, col = col, lty = 1, lwd = 2, bty = "n")  
  }
  
}




plot_mean_functions <- function(j, time_g, mu_g, list_mu_hat, col_sim = "grey", col = "black") {
  
  if (!is.null(mu_g)) {
    plot(time_g, mu_g[[j]], type= "l",xlab = " ", ylab= " ", 
         main=paste("Mean function, j = ", j), col = col_sim)
    lines(time_g, list_mu_hat[[j]], type="l", col= col, lwd = 2)
  } else {
    plot(time_g, list_mu_hat[[j]], type= "l",xlab = " ", ylab= " ", 
         main=paste("Mean function, j = ", j), col= col, lwd = 2)
  }

}

plot_factor_eigenfunctions <- function(q, time_g, Psi_g, list_list_Phi_hat, list_cumulated_pve = NULL, cum_pve_thres = 0.99) {
  
  mar = c(3.7, 4, 2.1, 1.1) 
  par(mar = mar)
  
  L_max <- ncol(list_list_Phi_hat[[q]])
  
  ind_pve <- which(list_cumulated_pve[[q]] < cum_pve_thres*100)
  L_cum_max <- min(ind_pve[length(ind_pve)] + 1, L_max)

  if (!is.null(Psi_g)) {
    ylim <- c(min(c(as.vector(Psi_g[[q]][,1:L_cum_max, drop = F]), as.vector(list_list_Phi_hat[[q]][,1:L_cum_max, drop = F]))),
              max(c(as.vector(Psi_g[[q]][,1:L_cum_max, drop = F]), as.vector(list_list_Phi_hat[[q]][,1:L_cum_max, drop = F]))))
  } else {
    ylim <- c(min(as.vector(list_list_Phi_hat[[q]][,1:L_cum_max, drop = F])),
              max(as.vector(list_list_Phi_hat[[q]][,1:L_cum_max, drop = F])))
  }

  for (l in 1:L_max) {
    if (is.nan(list_cumulated_pve[[q]][l])) { # not capturing any variation - factor PPI likely close to zero?
      
      plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=c(0,1))
      
    } else {
      
      if (!is.null(list_cumulated_pve) && l>1 && list_cumulated_pve[[q]][l-1] > cum_pve_thres*100) { # this eigenfunction captures less than 1% of the variation
        plot(1, type="n", xlab="", ylab= ifelse(l == 1, paste("q = ", q), ""), xlim=c(0, 1), ylim=ylim)
      } else {
        if (!is.null(Psi_g)) {
          plot(time_g, Psi_g[[q]][,l],type= "l",xlab = paste("l = ", l), ylab= ifelse(l == 1, paste("q = ", q), ""), ylim = ylim,
               main=paste0(ifelse(is.null(list_cumulated_pve), "", 
                                 paste0("cum PVE = ", format(list_cumulated_pve[[q]][l], digits = 2), "%")))) # [[q]][, l]
          lines(time_g, list_list_Phi_hat[[q]][,l], type="l", col= "black", lwd = 2)
        } else {
          plot(time_g, list_list_Phi_hat[[q]][,l], type= "l",xlab = paste("l = ", l), ylab= ifelse(l == 1, paste("q = ", q), ""),  ylim = ylim,
               main=paste0(ifelse(is.null(list_cumulated_pve), "", paste0("cum PVE = ", 
                                                                         format(list_cumulated_pve[[q]][l], digits = 2), "%"))), # [[q]][, l]
               col= "black", lwd = 2) # [[q]][, l]
        }
        abline(h = 0, lty = 2)
      } 
    }
  }
  
}



plot_eigenfunctions <- function(q,l, time_g, Psi_g, list_list_Phi_hat, list_cumulated_pve = NULL, cum_pve_thres = 0.99,
                                col_sim = "grey", col = "black") {
  
  if (is.nan(list_cumulated_pve[[q]][l])) { # not capturing any variation - factor PPI likely close to zero?
    
    plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 1))
    
  } else {
   
    if (!is.null(list_cumulated_pve) && l>1 && list_cumulated_pve[[q]][l-1] > cum_pve_thres*100) { # this eigenfunction captures less than 1% of the variation
      plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 1))
    } else {
      if (!is.null(Psi_g)) {
        plot(time_g, Psi_g[[q]][,l],type= "l",xlab = " ", ylab= " ", 
             main=paste("Factor q = ", q, ", FPC l = ", l, "\n",
                        ifelse(is.null(list_cumulated_pve), "", 
                               paste0("cumulated PVE = ", format(list_cumulated_pve[[q]][l], digits = 2), "%"))),
             col = col_sim, lwd = 2 ) # [[q]][, l]
        lines(time_g, list_list_Phi_hat[[q]][,l], type="l", col= col, lwd = 2)
      } else {
        plot(time_g, list_list_Phi_hat[[q]][,l], type= "l",xlab = " ", ylab= " ", 
             main=paste("Factor q = ", q, ", FPC l = ", l, "\n",
                        ifelse(is.null(list_cumulated_pve), "", paste0("cumulated PVE = ", 
                                                                       format(list_cumulated_pve[[q]][l], digits = 2), "%"))), # [[q]][, l]
             col= col, lwd = 2) # [[q]][, l]
      }
    } 
  }
}


plot_trajectory <- function(i, j, time_obs, time_g, Y, Y_hat, Y_low, Y_upp, days_thres = NULL, 
                            df_comb = NULL, offset = 0.1, bool_disp_curves = T,
                            Y_true = NULL, col_true = "grey60") {

  # Y_true holds the noise-free simulated trajectories on time_g, with the same
  # nested structure as Y_hat. It is only available in simulations, where the
  # truth is known, and is left NULL for the application.

  fac <- ifelse(!is.null(days_thres), days_thres-1, 1) # scaling factor back to the original time window
  xlim_max <- ifelse(!is.null(days_thres), days_thres, 1)

  
  if (is.null(Y_hat)) {
    list_ylim <- c(min(Y[[i]][[j]])-offset,
                   max(Y[[i]][[j]])+offset)
  } else {
    list_ylim <- c(min(Y[[i]][[j]],
                       Y_low[[i]][[j]])-offset,
                   max(Y[[i]][[j]],
                       Y_upp[[i]][[j]])+offset)
  }

  if (!is.null(Y_true)) {
    list_ylim <- c(min(list_ylim[1], min(Y_true[[i]][[j]]) - offset),
                   max(list_ylim[2], max(Y_true[[i]][[j]]) + offset))
  }
  
  if (!is.null(df_comb) && any(df_comb$severity == "HC")) {
    lo <- quantile(df_comb[df_comb$severity == "HC", names(Y[[i]])[j]], probs = 0.25, na.rm = TRUE)
    up <- quantile(df_comb[df_comb$severity == "HC", names(Y[[i]])[j]], probs = 0.75, na.rm = TRUE)
    
    list_ylim <- c(min(list_ylim[1], lo), max(list_ylim[2], up))
  }
  
  subj_id <- names(Y)[i]
  
  if (!is.null(df_comb)) {
    df_comb$time <- df_comb$hospital_outcome_days_from_start_sx_or_first_pos_swab
    df_comb$time[!(df_comb$hospital_outcome %in% "dead")] <- df_comb$last_observation_days_from_start_sx_or_first_pos_swab[!(df_comb$hospital_outcome %in% "dead")]
    df_comb$hospital_outcome[df_comb$severity %in% c("HC", "A", "B")] <- "alive"
    df_comb$status <- ifelse(df_comb$hospital_outcome == "dead", 2, 1) # 1 = cencored, 2 = dead
    
    status <- unique(df_comb$status[df_comb$subject_id %in% subj_id])
    time <- unique(df_comb$time[df_comb$subject_id %in% subj_id])
    time <- time[!is.na(time)]
    
    stopifnot(length(status) == 1 & length(time) == 1)
  }
    
  time_grid_subj_id <- time_g
  if (!is.null(days_thres)) {
    time_grid_subj_id <- time_grid_subj_id * fac
  }

  plot(fac*time_obs[[i]], Y[[i]][[j]], type = "p", pch = 19, xlab = "time", ylab="fit", 
        main = paste0("Subject ", i, ", variable ",j, 
                      ifelse(!is.null(df_comb), paste0("\n",
                      "(", unique(df_comb$severity[df_comb$subject_id %in% subj_id]), 
                      ", gender: ", unique(df_comb$gender[df_comb$subject_id %in% subj_id]), 
                      ", age: ", unique(df_comb$age[df_comb$subject_id %in% subj_id]), 
                      ")"), "")),
       ylim = list_ylim, 
       xlim = c(0, xlim_max),
       col = "grey40")
  
  # drawn before the fitted curve so that the fit remains legible on top of it
  if (!is.null(Y_true)) {
    # drawn thicker than the fitted curve so that it stays visible where the two
    # coincide, which is most of the domain when the fit is good
    lines(time_grid_subj_id, Y_true[[i]][[j]], col = col_true, lwd = 4)
  }

  if (bool_disp_curves) {
    if (!is.null(df_comb) && (status == 2 & time < days_thres)) { # cut prediction if patient died within the analysis window
      lines(time_grid_subj_id[time_grid_subj_id<=time],
            Y_hat[[i]][[j]][time_grid_subj_id<=time],
            col="black",lwd=2)
      points(time_grid_subj_id[which.min(abs(time_grid_subj_id-time))], 
             Y_hat[[i]][[j]][which.min(abs(time_grid_subj_id-time))], 
             pch = 4, cex = 1.25, lwd = 1.5)
      if (!is.null(Y_upp)){
        lines(time_grid_subj_id[time_grid_subj_id<=time], 
              Y_upp[[i]][[j]][time_grid_subj_id<=time],
              col="black",lwd=2,lty=2)
      }
      if (!is.null(Y_low)){
        lines(time_grid_subj_id[time_grid_subj_id<=time], 
              Y_low[[i]][[j]][time_grid_subj_id<=time],
              col="black",lwd=2,lty=2)
      }
    } else {
      lines(time_grid_subj_id,
            Y_hat[[i]][[j]],
            col="black",lwd=2)
      if (!is.null(Y_upp)){
        lines(time_grid_subj_id, 
              Y_upp[[i]][[j]],
              col="black",lwd=2,lty=2)
      }
      if (!is.null(Y_low)){
        lines(time_grid_subj_id, 
              Y_low[[i]][[j]],
              col="black",lwd=2,lty=2)
      }
      
    }
  }
  
  if (!is.null(df_comb) && any(df_comb$severity == "HC")) {
    rect(-5, lo, days_thres+5, up, col = adjustcolor("gray", alpha.f=0.5), border = "gray")
  }
  
}





plot_interpolated_factors_mefisto <- function(list_h_hat,
                                              list_h_hat_variance,
                                              time_g,
                                              factors = "all",
                                              list_h_true = NULL) {
  
  library(ggplot2)
  library(dplyr)
  
  N <- length(list_h_hat)
  Q <- ncol(list_h_hat[[1]])
  
  if (length(list_h_hat_variance) != N)
    stop("list_h_hat_variance must have the same structure as list_h_hat")
  
  df <- do.call(rbind, lapply(seq_len(N), function(i) {
    do.call(rbind, lapply(seq_len(Q), function(q) {
      out <- data.frame(
        group  = factor(paste0("subject_", i),
                        levels = paste0("subject_", seq_len(N))),
        factor = factor(q),
        time   = time_g,
        mean   = list_h_hat[[i]][, q],
        variance = list_h_hat_variance[[i]][, q]
      )
      if (!is.null(list_h_true)) {
        out$true <- list_h_true[[i]][, q]
      }
      out
    }))
  }))
  
  df <- df %>%
    mutate(
      sd   = sqrt(variance),
      ymin = mean - 1.96 * sd,
      ymax = mean + 1.96 * sd
    )
  
  if (factors != "all") {
    df <- df %>% filter(factor %in% factors)
  }
  
  gg_interpol <- ggplot(df, aes(x = time, colour = group))
  
  if (!is.null(list_h_true)) {
    gg_interpol <- gg_interpol +
      geom_line(aes(y = true),
                colour = "black",
                linewidth = 0.5)
  }
  
  gg_interpol +
    ## mean curve
    geom_line(aes(y = mean),
              linewidth = 0.7) +
    
    ## credible interval bounds (dashed)
    geom_line(aes(y = ymin),
              linetype = "dashed",
              linewidth = 0.7) +
    geom_line(aes(y = ymax),
              linetype = "dashed",
              linewidth = 0.7) +
    
    ## optional ribbon (kept subtle)
    geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = group),
                alpha = 0.15,
                colour = NA) +
    
    facet_wrap(
      ~ group,
      labeller = as_labeller(function(x)
        gsub("subject_", "Subject ", x))
    ) +
    
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    ) +
    labs(
      x = "Time",
      y = "Factor value"
    )
}




plot_interpolated_factors_bayesSYNC <- function(list_h_hat,
                                                list_h_low,
                                                list_h_upp,
                                                time_g,
                                                factors = "all",
                                                list_h_true = NULL) {
  
  library(ggplot2)
  library(dplyr)
  
  N <- length(list_h_hat)
  Q <- length(list_h_hat[[1]])
  
  if (length(list_h_low) != N || length(list_h_upp) != N)
    stop("list_h_low and list_h_upp must have the same structure as list_h_hat")
  
  df <- do.call(rbind, lapply(seq_len(N), function(i) {
    do.call(rbind, lapply(seq_len(Q), function(q) {
      out <- data.frame(
        group  = factor(paste0("subject_", i),
                        levels = paste0("subject_", seq_len(N))),
        factor = factor(q),
        time   = time_g,
        mean   = list_h_hat[[i]][[q]],
        ymin   = list_h_low[[i]][[q]],
        ymax   = list_h_upp[[i]][[q]]
      )
      if (!is.null(list_h_true)) {
        out$true <- list_h_true[[i]][, q]
      }
      out
    }))
  }))
  
  if (factors != "all") {
    df <- df %>% filter(factor %in% factors)
  }
  
  gg_interpol <- ggplot(df, aes(x = time, colour = group))
  
  if (!is.null(list_h_true)) {
    gg_interpol <- gg_interpol +
      geom_line(aes(y = true),
                colour = "black",
                linewidth = 0.5)
  }
  
  gg_interpol +
    geom_line(aes(y = mean), linewidth = 0.7) +
    geom_line(aes(y = ymin), linetype = "dashed", linewidth = 0.7) +
    geom_line(aes(y = ymax), linetype = "dashed", linewidth = 0.7) +
    geom_ribbon(
      aes(ymin = ymin, ymax = ymax, fill = group, group = group),
      alpha = 0.15,
      colour = NA
    ) +
    facet_wrap(
      ~ group,
      labeller = as_labeller(function(x)
        gsub("subject_", "Subject ", x))
    ) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    ) +
    labs(x = "Time", y = "Factor value")
}


mean_se <- function(x) {
  x <- as.numeric(x)
  c(mean = mean(x, na.rm = TRUE),
    se   = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
}

collapse_mat <- function(mat_list) {
  unlist(lapply(mat_list, function(m) as.numeric(m)))
}

paired_boxplot <- function(x1, x2, labels = c("bayesSYNC", "MEFISTO"),
                           main = "", ylab = "") {
  
  stopifnot(length(x1) == length(x2))
  
  boxplot(
    list(x1, x2),
    names = labels,
    main  = main,
    ylab  = ylab,
    outline = FALSE
  )
  
}

