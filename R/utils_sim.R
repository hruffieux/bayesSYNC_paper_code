library(pROC)

create_named_list <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))
}


check_structure <- function(x, struct, type, size = NULL,
                            null_ok = FALSE,  inf_ok = FALSE, na_ok = FALSE) {
  if (type == "double") {
    bool_type <-  is.double(x)
    type_mess <- "a double-precision "
  } else if (type == "integer") {
    bool_type <- is.integer(x)
    type_mess <- "an integer "
  } else if (type == "numeric") {
    bool_type <- is.numeric(x)
    type_mess <- "a numeric "
  } else if (type == "logical") {
    bool_type <- is.logical(x)
    type_mess <- "a boolean "
  } else if (type == "string") {
    bool_type <- is.character(x)
    type_mess <- "string "
  }
  
  bool_size <- TRUE # for case size = NULL (no assertion on the size/dimension)
  size_mess <- ""
  if (struct == "vector") {
    bool_struct <- is.vector(x) & (length(x) > 0) # not an empty vector
    if (!is.null(size)) {
      bool_size <- length(x) %in% size
      size_mess <- paste0(" of length ", paste0(size, collapse=" or "))
    }
  } else if (struct == "matrix") {
    bool_struct <- is.matrix(x) & (length(x) > 0) # not an empty matrix
    if (!is.null(size)) {
      bool_size <- all(dim(x) == size)
      size_mess <- paste0(" of dimension ", size[1], " x ", size[2])
    }
  }
  
  correct_obj <- bool_struct & bool_type & bool_size
  
  bool_null <- is.null(x)
  
  if (!is.list(x) & type != "string") {
    na_mess <- ""
    if (!na_ok) {
      if (!bool_null) correct_obj <- correct_obj & !any(is.na(x))
      na_mess <- " without missing value"
    }
    
    inf_mess <- ""
    if (!inf_ok) {
      if (!bool_null) correct_obj <- correct_obj & all(is.finite(x[!is.na(x)]))
      inf_mess <- ", finite"
    }
  } else {
    na_mess <- ""
    inf_mess <- ""
  }
  
  null_mess <- ""
  if (null_ok) {
    correct_obj <- correct_obj | bool_null
    null_mess <- " or must be NULL"
  }
  
  if(!(correct_obj)) {
    stop(paste0(deparse(substitute(x)), " must be a non-empty ", type_mess, struct,
                size_mess, inf_mess, na_mess, null_mess, "."))
  }
}



factor_recoveries <-function(true_B, B_hat, mu_q_gamma,list_Zeta_hat, list_list_Phi_hat){
  p <- dim(true_B)[1]
  Q <- dim(true_B)[2]
  true_factors <- lapply(1:Q, function(q) which(true_B[, q] != 0))
  rev_factors <- lapply(1:Q, function(q) which(mu_q_gamma[, q]>0.5))
  perm_factors <- rep (0, Q)
  prop_factors <- rep (0,Q)
  perm_sign <- rep(0,Q)
  perm_B_hat <- matrix (0, nrow = p, ncol=Q)
  perm_mu_q_gamma <- matrix (0, nrow = p, ncol=Q)
  perm_list_Zeta_hat <- vector("list", length = Q)
  perm_list_list_Phi_hat <- vector("list", length = Q)
  for (q in 1:Q){
    for (q_tilde in 1:Q){
      if (all(true_factors[[q]] %in% rev_factors[[q_tilde]])){
        perm_factors[q] <- q_tilde
        prop_factors[q] <- 1
      }
    }
    if (perm_factors[q]==0){ #The factor is not totally recovered. 
      shared_ele <-sapply(1:Q, function(q_tilde)
          sum(sapply(true_factors[q], function(e) e %in% rev_factors[[q_tilde]]))
      )
      perm_factors[q] <- which.max(shared_ele)
      prop_factors[q] <- shared_ele[which.max(shared_ele)]/ length(true_factors[[q]])
    }
    
    perm_sign[q] <- ifelse(sign(true_B[true_factors[[q]][1],q]) == sign(B_hat[true_factors[[q]][1],perm_factors[q]]), 1, -1)
    perm_B_hat[,q]<- perm_sign[q]*B_hat[,perm_factors[q]]
    perm_mu_q_gamma[,q]<- mu_q_gamma[,perm_factors[q]]
    perm_list_Zeta_hat[[q]] <- perm_sign[q]*list_Zeta_hat[[perm_factors[q]]]
    perm_list_list_Phi_hat[[q]]<- list_list_Phi_hat[[perm_factors[q]]]
  }
  res<- create_named_list(perm_factors, prop_factors, perm_sign,perm_B_hat, 
                          perm_mu_q_gamma, perm_list_Zeta_hat,perm_list_list_Phi_hat)
}

auc_roc_curve <- function(B, mu_q_gamma, plot = F) {
  
  # Create a new plot for each factor loading
  Q_true <- ncol(B)
  par(mfrow=c(1, Q_true), mar=c(5, 5, 4, 2))
  
  auc_score <- rep(1,Q_true)
  
  for (q in 1:Q_true) {
    true_pattern <- ifelse(B[,q]!=0, 1, 0)
    if (length(unique(true_pattern)) != 1){
      # alg_factor <- ifelse(mu_q_gamma[, q] > 0.5, 1, 0)
      # roc_obj <- roc(true_pattern, alg_factor)
      roc_obj <- roc(true_pattern, mu_q_gamma[, q])
      auc_score[q]<- auc(roc_obj)
      if (plot == T) {
        plot(roc_obj, main=paste("ROC Curve - Factor", q), col="dodgerblue3")
        legend("bottomright", legend=paste("AUC =", round(auc(roc_obj), 2)), bty = "n", col="dodgerblue3", lty=1, lwd = 2, cex=0.8)
      }
    }
  }
  return (auc_score)
}


sign_recoveries <- function(list_Zeta_hat, Zeta,list_list_Phi_hat){
  Q <- dim(Zeta[[1]])[2]
  L <- dim(Zeta[[1]])[1]
  
  perm_list_Zeta_hat <- list_Zeta_hat
  perm_list_list_Phi_hat <- list_list_Phi_hat
  
  perm_matrix <- matrix (0, nrow= Q, ncol = L)
  perm_sign_fpca <- matrix(0,nrow = Q, ncol = L)
  for (q in 1:Q){
    for (l in 1:L){
      Zeta_true <- sapply(Zeta, function(Zeta_i) Zeta_i[l,q])
      #Zeta_hat_orth <- sapply(1:nrow(perm_list_Zeta_hat[[q]]),function(i) perm_list_Zeta_hat[[q]][i,l])
      corr_zeta <- sapply (1:L, function(l_tilde) 
                           sapply (1:Q, function(q_tilde) 
                             cor(Zeta_true, sapply(1:nrow(list_Zeta_hat[[q_tilde]]),function(i) 
                               list_Zeta_hat[[q_tilde]][i,l_tilde]))))
      corr_zeta[is.na(corr_zeta)] <- 0
      q_tilde <- which(abs(corr_zeta) == max(abs(corr_zeta)), arr.ind = TRUE)[1,1]
      l_tilde <- which(abs(corr_zeta) == max(abs(corr_zeta)), arr.ind = TRUE)[1,2]
      perm_sign_fpca[q,l] <- sign(corr_zeta[q_tilde,l_tilde])
      perm_list_Zeta_hat[[q]][,l] <- sign(corr_zeta[q_tilde,l_tilde])*list_Zeta_hat[[q_tilde]][,l_tilde]
      perm_list_list_Phi_hat[[q]][,l]<- sign(corr_zeta[q_tilde,l_tilde])*list_list_Phi_hat [[q_tilde]][,l_tilde]
    }
  }
  res<- create_named_list(perm_sign_fpca, perm_list_Zeta_hat, perm_list_list_Phi_hat)
}

corr_zeta_vector <- function(Zeta, list_Zeta_hat){

  Q <- length(Zeta)
  N <- nrow(Zeta[[1]])
  L <- ncol(Zeta[[1]])
  
  corr_zeta <- matrix (0, nrow = Q, ncol = L)
  for (q in 1:Q){
    for (l in 1:L){
      zeta_ql <- sapply(1:N, function(i) Zeta[[q]][i,l])
      zeta_hat_orth_ql <- sapply(1:N, function(i) list_Zeta_hat[[q]][i,l])
      corr_zeta[q,l] <- cor(zeta_ql, zeta_hat_orth_ql)
    }
  }
  return(corr_zeta)
}

corr_b_vector <-function(B, B_hat){
  Q <- ncol(B)
  corr_b <- sapply(1:Q, function(q) cor(B[,q], B_hat[,q]))
  return(corr_b)
}

compute_ise <- function(time_g, fct_g, fct_est) {
  
  trapint(time_g, (fct_est - fct_g)^2) 
  
}

error_mean_function <- function(time_g, list_mu_hat,mu_g ){
  p <- length (list_mu_hat)
  ise_mu_tot <- 0
  ise_mu_vec <- rep (0, p)
  for (j in 1:p){
    ise_mu_vec[j] <- trapint(time_g, (list_mu_hat[[j]] - mu_g[[j]])^2) 
    ise_mu_tot <- ise_mu_tot+ise_mu_vec[j]
  }
  return(ise_mu_tot)
}

error_eigen_function <- function(time_g, Psi_g, C_g, list_list_Phi_hat){
  
  Q_true <- length(Psi_g)
  L_true <- ncol(Psi_g[[1]])
  ise_mfpca_vec <- matrix(0, nrow = Q_true, ncol = L_true)
  for (q in 1:Q_true){
    for (l in 1:L_true){
      ise_mfpca_vec[q,l] <- trapint(time_g, (list_list_Phi_hat[[q]][,l] - Psi_g[[q]][,l])^2)
    }
  }
  return(ise_mfpca_vec)
}

plot_correlation_b <- function (q, Corr_b){
 corr_b_q <- sapply(Corr_b, function(cor) cor[q])
 boxplot(corr_b_q, 
         names = c(q),
         xlab = "",
         ylab = "Correlation",
         ylim = c(0, 1),
         las = 2,
         main = "Correlation factors")
 return(corr_b_q)
}

plot_auc_score <- function(q, Auc_b){
  auc_b_q <- sapply(Auc_b, function(auc) auc[q])
  boxplot(auc_b_q, 
          names = c(q),
          xlab = "",
          ylab = "AUC",
          las = 2,
          ylim = c(0, 1),
          main = paste0("AUC factor q = ", q))
  return(auc_b_q)
}

plot_correlation_zeta <- function(q, l, Corr_zeta){
  corr_zeta_ql <- sapply(Corr_zeta, function(cor) cor[q,l])
  boxplot(corr_zeta_ql, 
          names = c(q),
          xlab = "",
          ylab = "Correlation",
          las = 2,
          ylim = c(0, 1),
          main = paste0("Correlation scores, q = ", q, ", l = ", l))
  return(corr_zeta_ql)
}

plot_err_eigen <- function(q, l, Err_eigen, ylim){
  eigen_ql <- sapply(Err_eigen, function(err) err[q,l])
  boxplot(eigen_ql, 
          names = c(q),
          xlab = "",
          ylab = "ISE",
          las = 2,
          ylim = ylim,
          main = paste0("Errors eigenfunction, q = ", q, ", l = ", l))
  return(eigen_ql)
}

plot_err_mean <- function(j, Err_mean){
  err_mean_j <- sapply(Err_mean, function(err) err[j])
  boxplot(err_mean_j, 
          names = c(j),
          xlab = "",
          ylab = "ISE",
          las = 2,
          main = paste0("Errors mean function, j = ", j)
  )
  return(err_mean_j)
}
  

generate_bayesSYNC_data_mefisto <- function (N, n, p, Q, L, mu_func, 
                                     time_obs = NULL,
                                     n_g = 1000, 
                                     a_om = 1, b_om = p, 
                                     vec_sd_zeta = NULL, 
                                     vec_sd_eps = NULL,
                                     seed = NULL,
                                     vec_mean_zeta = NULL){
  
  check_structure(seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  check_structure(n, "vector", "numeric", N)
  
  if (is.null(time_obs)) {
    time_obs <- lapply(n, runif)
    time_obs <- lapply(time_obs, sort)
  }
  
  if (is.null(vec_sd_zeta)) {
    vec_sd_zeta <- 1/(1:L)
  } else {
    stopifnot(length(vec_sd_zeta) == L)
  }
  
  if (!is.null(vec_sd_zeta)) {
    stopifnot(length(vec_mean_zeta) == L)
  }
  
  if(is.null(vec_sd_eps)){
    vec_sd_eps <- rep(1, p) 
  } else if (length(vec_sd_eps) == 1) {
    vec_sd_eps <- rep(vec_sd_eps, p)
  } else {
    stopifnot(length(vec_sd_eps) == p)
  }
  
  omega <- rbeta(Q, shape1 = a_om, shape2 = b_om)
  gamma <- matrix(rbinom(p*Q, size = 1, prob = omega), nrow = p, ncol = Q)
  
  # whether all factors should be "active", 
  # i.e., represent the contribution of at least one variable j = 1, ..., p
  bool_active <- T  
  if (bool_active) {
    if(any(colSums(gamma) == 0)) {
      for (q in which(colSums(gamma) == 0)) {
        gamma[sample(1:p, 1), q] <- 1
      }
    }
  }
  
  B_normal <- matrix(rnorm(p*Q, mean = 0, sd = 1), nrow = p, ncol = Q)
  
  B <- B_normal * gamma
  
  if (is.null(vec_mean_zeta)) {
    Zeta <- lapply(1:Q, function(q) { # list of length Q, containing matrices of size N x L 
      matrix(rnorm(N * L, 0, rep(vec_sd_zeta, each = N)), nrow = N, ncol = L)
    })
  } else {
    Zeta <- lapply(1:Q, function(q) { # list of length Q, containing matrices of size N x L 
      matrix(rnorm(N * L, rep(vec_mean_zeta, each = N), rep(vec_sd_zeta, each = N)), nrow = N, ncol = L)
    })
  }

  
  mu <- lapply(1:N, function(i) {
    lapply(1:p, function(j) {
      mu_func(time_obs[[i]], j = j, p = p)
    })
  })
  
  res_splines <- construct_splines(Q = Q, L = L)
  splines <- res_splines$splines
  vec_K <- res_splines$vec_K
  
  list_ind_spline_functions <- lapply(1:Q, function(q) {
    sample(1:(ncol(Splinets::evspline(splines[[1]]$os))-1), L) # remove 1 as the first column is the vector of timepoints
  })
  
  Phi <- lapply(1:Q, function(q) {
    lapply(1:N, function(i) {
      eval_qi <- Splinets::evspline(splines[[q]]$os, x = time_obs[[i]])[,-1,drop =FALSE]
      eval_qi[, list_ind_spline_functions[[q]], drop = F]
      # eval_qi[, 1:L, drop = F]
    })
  })
  
  Y <- lapply(1:N, function(i) { 
    Y_i <- lapply(1:p, function(j) {
      
      tmp <- sapply(1:Q, function(q) B[j,q] * Phi[[q]][[i]] %*% Zeta[[q]][i, ])
      if (is.vector(tmp)) { # to deal with cases with a single observation
        tmp <- t(as.matrix(tmp))
      }
      mu[[i]][[j]] + 
        rowSums(tmp) + 
        rnorm(length(time_obs[[i]]), mean = 0, sd = vec_sd_eps[j])
    })
    names(Y_i) <- paste0("variable_", 1:p)
    Y_i
  })
  
  names(Y) <- paste0("subject_", 1:N)
  
  time_g <- seq(0, 1, length.out = n_g)
  
  mu_g <- lapply(1:p, function(j) mu_func(time_g, j = j, p = p))
  
  Phi_g <- lapply(1:Q, function(q) {
    eval_q <- Splinets::evspline(splines[[q]]$os, x = time_g)[,-1,drop =FALSE]
    eval_q[, list_ind_spline_functions[[q]], drop = F]
    # eval_qi[, 1:L, drop = F]
  })
  
  
  list_h_true <- lapply(1:N, function(i) { 
    tmp <- sapply(1:Q, function(q) Phi_g[[q]] %*% Zeta[[q]][i, ])
    if (is.vector(tmp)) { # to deal with cases with a single observation
      tmp <- t(as.matrix(tmp))
    }
    colnames(tmp) <- paste0("Factor", 1:Q)
    tmp
  })
  names(list_h_true) <- paste0("subject_", 1:N)
  
  Y_true <- lapply(1:N, function(i) { 
    Y_i <- lapply(1:p, function(j) {
      
      tmp <- sapply(1:Q, function(q) B[j,q] * Phi_g[[q]] %*% Zeta[[q]][i, ])
      if (is.vector(tmp)) { # to deal with cases with a single observation
        tmp <- t(as.matrix(tmp))
      }
      mu_g[[j]] + rowSums(tmp)
    })
    names(Y_i) <- paste0("variable_", 1:p)
    Y_i
  })
  names(Y_true) <- paste0("subject_", 1:N)
  
  create_named_list(time_obs, Y, B, omega, gamma, mu, Phi, Zeta,
                    vec_sd_zeta, vec_sd_eps, 
                    time_g, mu_g, Phi_g,  Y_true, list_h_true, vec_K)
}



match_factor_and_sign_mefisto <- function(B, B_hat, list_h_hat, list_h_hat_obs, list_h_hat_variance = NULL) {
  
  perm_factor <- match_factors_mefisto(B, B_hat, list_h_hat, list_h_hat_obs, list_h_hat_variance)
  
  best_perm <- perm_factor$best_perm
  perm_sign_factor <- perm_factor$perm_sign_factor
  
  perm_B_hat <- perm_factor$perm_B_hat
  perm_B_hat_untrimmed <- perm_factor$perm_B_hat_untrimmed
  
  perm_list_h_hat <- perm_factor$perm_list_h_hat
  perm_list_h_hat_untrimmed <- perm_factor$perm_list_h_hat_untrimmed
  
  perm_list_h_hat_obs <- perm_factor$perm_list_h_hat_obs
  perm_list_h_hat_obs_untrimmed <- perm_factor$perm_list_h_hat_obs_untrimmed
  
  perm_list_h_hat_variance <- perm_factor$perm_list_h_hat_variance
  perm_list_h_hat_variance_untrimmed <- perm_factor$perm_list_h_hat_variance_untrimmed
  
  
  perm_sign_fpca <- perm_factor$perm_sign_fpca
  
  bool_rescale_loadings_and_scores <- T # so they match the simulated loadings and scores, respectively
  if (bool_rescale_loadings_and_scores) {
    # B is identifiable up to multiplicative sign on its columns
    norm_col_B <- sqrt(colSums(B^2)) # this is unknown in practice.
    norm_col_B_hat <- sqrt(colSums(perm_B_hat^2))
    
    Q_true <- ncol(B)
    for (q in 1:Q_true) {
      perm_B_hat[,q] <- perm_B_hat_untrimmed[,q] <- perm_B_hat[,q] * norm_col_B[q] / norm_col_B_hat[q]
    }
    
    for(i in 1:length(list_h_hat)){
      for(q in 1:Q_true) {
        perm_list_h_hat[[i]][,q] <- perm_list_h_hat_untrimmed[[i]][,q] <- perm_list_h_hat[[i]][,q] * norm_col_B_hat[q] / norm_col_B[q] 
        perm_list_h_hat_obs[[i]][,q] <- perm_list_h_hat_obs_untrimmed[[i]][,q] <- perm_list_h_hat_obs[[i]][,q] * norm_col_B_hat[q] / norm_col_B[q] 
        
        if (!is.null(list_h_hat_variance)) {
          perm_list_h_hat_variance[[i]][,q] <- perm_list_h_hat_variance_untrimmed[[i]][,q] <- perm_list_h_hat_variance[[i]][,q] * norm_col_B_hat[q]^2 / norm_col_B[q]^2 
        }
      }
    }
    
  }
  
  create_named_list(best_perm, perm_sign_factor, 
                    perm_B_hat, perm_list_h_hat, perm_list_h_hat_obs,
                    perm_list_h_hat_variance, perm_list_h_hat_variance_untrimmed,
                    perm_B_hat_untrimmed, perm_list_h_hat_untrimmed, perm_list_h_hat_obs_untrimmed)
  
}


# Example code to reshape Y into a data.frame
reshape_Y_to_df <- function(Y, time_obs) {
  # Initialize an empty data frame to store the results
  result <- data.frame(sample = integer(),
                       subject = character(),
                       time = numeric(),
                       feature = character(),
                       value = numeric(),
                       stringsAsFactors = FALSE)
  
  # Loop over the groups
  for (group_idx in seq_along(Y)) {
    group_name <- names(Y)[group_idx]
    group_data <- Y[[group_idx]]
    time_obs_i <- time_obs[[group_idx]]
    
    # Loop over the features within each group
    for (feature_idx in seq_along(group_data)) {
      feature_name <- names(group_data)[feature_idx]
      feature_values <- group_data[[feature_idx]]
      
      # Create a data frame for the current group and feature
      df <- data.frame(
        sample = paste0("sample_", seq_along(feature_values), "_subject_", group_idx),
        subject = group_name,
        time = time_obs_i,
        feature = feature_name,
        value = feature_values
      )
      
      # Append the data frame to the result
      result <- rbind(result, df)
    }
  }
  
  return(result)
}

match_factors_mefisto <- function(B, B_hat, list_h_hat, list_h_hat_obs, list_h_hat_variance = NULL) {
  
  Q_true <- ncol(B)
  Q <- ncol(B_hat)
  
  if (Q_true > Q) {
    stop("The number of estimated factors, Q, must be larger than the number of true factors.")
  } else if (Q_true < Q) {
    warning("Dropping superfluous factors based on the true loadings.")
  }
  true_pat <- 1*(abs(B) > 0)
  
  # Generate all possible permutations of column indices
  perms <- gtools::permutations(Q, Q_true)
  
  # Initialize minimum distance and best permutation
  max_abs_correlation <- -Inf
  best_perm <- NULL
  
  # Iterate over all permutations to find the best one
  for (i in 1:nrow(perms)) {
    permuted_B_hat <- B_hat[, perms[i, ], drop = F]
    abs_correlation <- sum(sapply(1:Q_true, function(q) abs(cor(B[,q], permuted_B_hat[,q]))), na.rm = T) 
    sign_correlation <- sapply(1:Q_true, function(q) sign(cor(B[,q], permuted_B_hat[,q])))
    
    if (abs_correlation > max_abs_correlation) {
      max_abs_correlation <- abs_correlation
      perm_sign_factor <- sign_correlation
      best_perm <- perms[i, ]
    }
  }
  
  # Reorder B with the best permutation found
  perm_list_h_hat <- lapply(list_h_hat, function(ll) sweep(ll[, best_perm, drop = F], 2, perm_sign_factor, FUN = "*"))
  perm_list_h_hat_obs <- lapply(list_h_hat_obs, function(ll) sweep(ll[, best_perm, drop = F], 2, perm_sign_factor, FUN = "*"))
  
  if (!is.null(list_h_hat_variance)) {
    perm_list_h_hat_variance <- lapply(list_h_hat_variance, function(ll) ll[, best_perm, drop = F])
    perm_list_h_hat_variance_untrimmed <- lapply(seq_along(list_h_hat_variance), function(i) cbind(perm_list_h_hat_variance[[i]], 
                                                                                          list_h_hat_variance[[i]][, -best_perm, drop = F]))
  } else {
    perm_list_h_hat_variance <- perm_list_h_hat_variance_untrimmed <- NULL
  }
  perm_B_hat <- sweep(B_hat[, best_perm, drop = F], 2, perm_sign_factor, FUN = "*")
  
  perm_B_hat_untrimmed <- cbind(perm_B_hat, B_hat[, -best_perm, drop = F])
  perm_list_h_hat_untrimmed <- lapply(seq_along(list_h_hat), function(i) cbind(perm_list_h_hat[[i]], 
                                                                                       list_h_hat[[i]][, -best_perm, drop = F]))
  perm_list_h_hat_obs_untrimmed <- lapply(seq_along(list_h_hat_obs), function(i) cbind(perm_list_h_hat_obs[[i]], 
                                                                               list_h_hat_obs[[i]][, -best_perm, drop = F]))
  
  create_named_list(best_perm, perm_sign_factor, 
                    perm_B_hat_untrimmed, perm_list_h_hat_untrimmed, perm_list_h_hat_obs_untrimmed,
                    perm_list_h_hat_variance, perm_list_h_hat_variance_untrimmed,
                    perm_B_hat, perm_list_h_hat, perm_list_h_hat_obs)
}



coverage_curve <- function(y_obs, time_obs, y_hat_low, y_hat_upp, time_g) {
  # Interpolate credible interval bounds at observed time points
  y_low_interp <- approx(time_g, y_hat_low, xout = time_obs)$y
  y_upp_interp <- approx(time_g, y_hat_upp, xout = time_obs)$y
  
  # Count how many observed values fall within the interpolated interval
  in_interval <- (y_obs >= y_low_interp) & (y_obs <= y_upp_interp)
  
  # Compute coverage proportion
  coverage <- mean(in_interval, na.rm = TRUE)  # Ignore NA values from interpolation
  
  return(coverage)
}


coverage_latent <- function(list_h_true, list_h_low, list_h_upp) {
  
  N <- length(list_h_true)
  Q <- ncol(list_h_true[[1]])
  
  cov_mat <- matrix(NA, nrow = N, ncol = Q)
  
  for (i in 1:N) {
    for (q in 1:Q) {
      
      h_true <- list_h_true[[i]][, q]
      h_low  <- list_h_low[[i]][[q]]
      h_upp  <- list_h_upp[[i]][[q]]
      
      in_interval <- (h_true >= h_low) & (h_true <= h_upp)
      cov_mat[i, q] <- mean(in_interval, na.rm = TRUE)
    }
  }
  
  cov_mat
}

width_latent <- function(list_h_low, list_h_upp) {
  
  N <- length(list_h_low)
  Q <- length(list_h_low[[1]])
  
  width_mat <- matrix(NA, nrow = N, ncol = Q)
  
  for (i in 1:N) {
    for (q in 1:Q) {
      width_mat[i, q] <- mean(
        list_h_upp[[i]][[q]] - list_h_low[[i]][[q]],
        na.rm = TRUE
      )
    }
  }
  
  width_mat
}



credible_band_width <- function(y_hat_low, y_hat_upp, time_g) {
  # Compute the pointwise credible band width
  band_width <- y_hat_upp - y_hat_low
  
  # Summary statistics
  mean_width <- mean(band_width, na.rm = TRUE)  # Average width
  median_width <- median(band_width, na.rm = TRUE)  # Median width
  max_width <- max(band_width, na.rm = TRUE)  # Maximum width
  min_width <- min(band_width, na.rm = TRUE)  # Minimum width
  
  # Integrated Band Width (IBW) using trapezoidal integration
  integrated_width <- sum(diff(time_g) * (head(band_width, -1) + tail(band_width, -1)) / 2)
  
  # Store results
  list(
    mean_width = mean_width,
    median_width = median_width,
    max_width = max_width,
    min_width = min_width,
    integrated_width = integrated_width,
    pointwise_width = band_width  # Optional: return full vector for visualization
  )
}





  