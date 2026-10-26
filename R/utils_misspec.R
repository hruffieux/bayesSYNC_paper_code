## ---------------------------------------------------------------------------
## Helper functions for the robustness-under-misspecification simulations.
##
## Two misspecification mechanisms are implemented, both as post-processing
## steps applied to data produced by generate_bayesSYNC_data(). The generator
## itself is left untouched, and so is the bayesSYNC fitting call: the model is
## always fitted under its own (white-noise, complete-follow-up) assumptions,
## which is exactly what makes the setting a misspecification study.
##
##   Scenario 1, correlated errors along t:
##     the data are generated with vec_sd_eps = 0 (noise-free signal) and the
##     error process is added afterwards with a user-specified temporal
##     correlation. Setting rho = 0 reproduces the white-noise case through the
##     identical code path, so it serves as the internal baseline.
##
##   Scenario 2, informative drop-out:
##     observations are removed from the first visit at which a subject drops
##     out, with a drop-out hazard depending on the current (unrecorded) value
##     of a latent factor process. The mechanism is therefore missing not at
##     random. A non-informative arm (beta_dropout = 0) calibrated to the same
##     expected retention is used as control.
## ---------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## Scenario 1: correlated errors along t
## ---------------------------------------------------------------------------

#' Correlation matrix of the error process at a set of observation times
#'
#' @param t_vec numeric vector of observation times for one subject.
#' @param rho range parameter, on the scale of T = [0, 1]. rho = 0 gives
#'   independent errors.
#' @param corr_type "exponential" (continuous-time analogue of AR(1),
#'   equivalently Matern with smoothness 1/2, non-differentiable paths) or
#'   "matern32" (Matern with smoothness 3/2, once-differentiable paths, hence a
#'   harder case since the error then looks like signal).
#'
#' @return an n_i x n_i correlation matrix.
build_corr_matrix <- function(t_vec, rho, corr_type = c("exponential", "matern32")) {
  
  corr_type <- match.arg(corr_type)
  n_i <- length(t_vec)
  
  if (n_i == 1 || rho <= 0) {
    return(diag(n_i))
  }
  
  d <- abs(outer(t_vec, t_vec, "-"))
  
  R <- switch(corr_type,
              exponential = exp(-d / rho),
              matern32 = {
                u <- sqrt(3) * d / rho
                (1 + u) * exp(-u)
              })
  
  R
}


#' Add an error process with temporal correlation to a noise-free dataset
#'
#' The marginal error variance is vec_sd_eps^2 whatever the value of rho, so
#' that the signal-to-noise ratio is held fixed across settings and only the
#' correlation structure changes. All p variables of a given subject share the
#' same observation times, hence a single Cholesky factor per subject.
#'
#' @param fpca_data output of generate_bayesSYNC_data() called with
#'   vec_sd_eps = rep(0, p), i.e. containing the noise-free signal in $Y.
#' @param vec_sd_eps vector of length p of error standard deviations.
#' @param rho range parameter of the error correlation, see build_corr_matrix().
#' @param corr_type correlation family, see build_corr_matrix().
#' @param prop_corr proportion of the error variance carried by the correlated
#'   component; the remaining 1 - prop_corr is a white-noise nugget. Default 1,
#'   i.e. the whole error process is correlated.
#' @param Y responses to which the error process is added. Defaults to the
#'   noise-free signal held in fpca_data$Y, and is supplied explicitly when a
#'   latent-level perturbation has already been added by add_latent_error(), so
#'   that the observation noise is always the last term added.
#' @param seed seed used for the error draws only, so that the signal and the
#'   noise come from separate, reproducible streams.
#'
#' @return the list Y with the error process added, same structure as fpca_data$Y.
add_correlated_noise <- function(fpca_data, vec_sd_eps, rho,
                                 corr_type = "exponential",
                                 prop_corr = 1, Y = NULL, seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  stopifnot(prop_corr >= 0, prop_corr <= 1)
  
  if (is.null(Y)) {
    Y <- fpca_data$Y
  }
  N <- length(Y)
  p <- length(Y[[1]])
  
  stopifnot(length(vec_sd_eps) == p)
  
  for (i in seq_len(N)) {
    
    t_i <- fpca_data$time_obs[[i]]
    n_i <- length(t_i)
    
    R <- build_corr_matrix(t_i, rho, corr_type)
    
    # mixture of a correlated component and a white-noise nugget; the diagonal
    # stays at one, so the marginal variance is unchanged
    R <- prop_corr * R + (1 - prop_corr) * diag(n_i)
    
    # small jitter for numerical stability when times are nearly tied
    C <- chol(R + 1e-8 * diag(n_i)) # upper triangular, R = t(C) %*% C
    
    # rows of E are independent draws with covariance R, then scaled variable-wise
    E <- matrix(rnorm(p * n_i), nrow = p, ncol = n_i) %*% C
    E <- E * vec_sd_eps # recycling is by row, i.e. variable j is scaled by vec_sd_eps[j]
    
    Y_i <- Y[[i]] # local copy to avoid repeated copying of the nested list
    for (j in seq_len(p)) {
      Y_i[[j]] <- Y_i[[j]] + E[j, ]
    }
    Y[[i]] <- Y_i
  }
  
  Y
}


## ---------------------------------------------------------------------------
## Scenario 1b: unmodelled smooth variation at the latent level
##
## The counterpart of add_correlated_noise() for the off-diagonal blocks of the
## marginal covariance. The implied covariance of the fitted model is
##
##   Cov(y_i^(j)(s), y_i^(j')(t))
##     = sum_q b_jq b_j'q sigma^(q)(s, t) + sigma^2_eps,j 1{j = j'} delta(s - t),
##
## so a variable-specific error process perturbs the diagonal blocks only, and
## leaves the loadings and the latent processes identified by the off-diagonal
## blocks. A perturbation of the latent trajectories themselves is transmitted
## by the loading matrix to every variable loading on the factor, and therefore
## enters the off-diagonal blocks, where it competes with the factor structure.
##
## The perturbation represents the part of the Karhunen-Loeve expansion left out
## by the truncation at L, that is, the tail sum_{l > L} zeta_il^(q) psi_l^(q)(t).
## It is generated as a smooth process rather than as explicit extra components,
## so that the target of the metrics stays exactly what it is in the baseline,
## namely mu, B, the first L_true eigenfunctions and their scores, with the
## perturbation acting purely as unmodelled nuisance.
## ---------------------------------------------------------------------------

#' Add an unmodelled smooth process to the latent factor trajectories
#'
#' The perturbation delta_i^(q) is drawn independently across subjects and
#' factors from a mean-zero process with correlation R and variance
#' kappa * Var(h^(q)), where the variance of the structured latent part is
#' estimated empirically over subjects and observation times. The parameter
#' kappa is therefore the proportion of latent variation that the fitted
#' expansion cannot represent, and is comparable across factors and replicates.
#'
#' A Matern 3/2 correlation with a short range is used by default: the discarded
#' components of a Karhunen-Loeve expansion carry higher frequencies than the
#' leading ones, while remaining smooth.
#'
#' @param fpca_data output of generate_bayesSYNC_data().
#' @param Y responses to perturb, typically the noise-free signal in
#'   fpca_data$Y, before the observation noise is added.
#' @param kappa variance of the unmodelled process, as a proportion of the
#'   variance of the structured latent part. kappa = 0 leaves Y unchanged.
#' @param rho_lat range parameter of the unmodelled process.
#' @param corr_type correlation family, see build_corr_matrix().
#' @param seed seed used for the perturbation draws only.
#'
#' @return a list with the perturbed Y, the per-factor standard deviations of
#'   the perturbation and the empirical variances of the structured part.
add_latent_error <- function(fpca_data, Y, kappa, rho_lat,
                             corr_type = "matern32", seed = NULL) {
  
  if (kappa <= 0) {
    return(list(Y = Y, vec_tau = NULL, vec_var_h = NULL))
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  B <- fpca_data$B
  Q <- ncol(B)
  N <- length(Y)
  p <- length(Y[[1]])
  
  # empirical variance of the structured latent part of each factor, over
  # subjects and observation times; the process has mean zero by construction
  vec_var_h <- sapply(seq_len(Q), function(q) {
    h_q <- unlist(compute_factor_path_obs(fpca_data, q_dropout = q))
    mean(h_q^2)
  })
  vec_tau <- sqrt(kappa * vec_var_h)
  
  for (i in seq_len(N)) {
    
    t_i <- fpca_data$time_obs[[i]]
    n_i <- length(t_i)
    
    R <- build_corr_matrix(t_i, rho_lat, corr_type)
    C <- chol(R + 1e-8 * diag(n_i)) # upper triangular, R = t(C) %*% C
    
    # rows of Delta are the perturbations of the Q latent trajectories at the
    # observation times of subject i
    Delta <- matrix(rnorm(Q * n_i), nrow = Q, ncol = n_i) %*% C
    Delta <- Delta * vec_tau # recycling is by row, i.e. factor q is scaled by vec_tau[q]
    
    # the perturbation reaches variable j through its loadings, exactly as the
    # structured latent trajectories do
    E <- B %*% Delta # p x n_i
    
    Y_i <- Y[[i]] # local copy to avoid repeated copying of the nested list
    for (j in seq_len(p)) {
      Y_i[[j]] <- Y_i[[j]] + E[j, ]
    }
    Y[[i]] <- Y_i
  }
  
  create_named_list(Y, vec_tau, vec_var_h)
}


## ---------------------------------------------------------------------------
## Scenario 1c: correlation between the latent factor processes
##
## The model assumes Cov(h^(q)(s), h^(q')(t)) = 0 for q != q'. Whether that
## assumption restricts what the model can fit depends on how the correlation
## arises. Mixing the trajectories themselves, h -> M h, is absorbed exactly:
## it is the same model with loadings B M, so nothing is lost and the departure
## would be invisible. The binding case is correlation between the scores of
## factors whose eigenfunctions differ, since the score of factor q' is then
## paired with the eigenfunction of factor q and the resulting covariance is no
## longer a sum of terms each combining one loading pattern with one covariance
## function, which is the only form the model can represent. Increasing Q does
## not recover it. This scenario generates that case.
## ---------------------------------------------------------------------------

#' Eigenfunction evaluations of one latent factor at the observation times
#'
#' Companion to compute_factor_path_obs(), which returns the trajectories; this
#' returns the basis evaluations themselves, needed when the scores are altered.
#'
#' @param fpca_data output of generate_bayesSYNC_data().
#' @param q index of the latent factor.
#'
#' @return a list of length N of n_i x L matrices of eigenfunction evaluations.
eval_eigenfunctions_obs <- function(fpca_data, q) {

  N <- length(fpca_data$time_obs)

  if (!is.null(fpca_data$Phi)) {

    lapply(seq_len(N), function(i) as.matrix(fpca_data$Phi[[q]][[i]]))

  } else {

    time_g <- fpca_data$time_g
    Phi_g_q <- fpca_data$Phi_g[[q]]
    L <- ncol(Phi_g_q)

    lapply(seq_len(N), function(i) {
      t_i <- fpca_data$time_obs[[i]]
      eval_i <- sapply(seq_len(L), function(l) {
        approx(time_g, Phi_g_q[, l], xout = t_i, rule = 2)$y
      })
      if (is.vector(eval_i)) { # single observation time
        eval_i <- t(as.matrix(eval_i))
      }
      eval_i
    })
  }
}


#' Correlate the latent factor processes through their scores
#'
#' For each FPCA component l, the vector of scores across the Q factors is
#' replaced by a correlated version with the same marginal variances, using the
#' Cholesky factor of an equicorrelation matrix. Components are treated
#' independently of one another, so only the dependence across factors changes.
#' Each factor keeps its own eigenfunctions, which is what makes the departure
#' one the model cannot represent.
#'
#' The perturbation is applied as a difference and propagated through the
#' loadings, so that everything else in the replicate, including the loadings,
#' the eigenfunctions, the mean functions and the observation times, is left
#' exactly as it was.
#'
#' @param fpca_data output of generate_bayesSYNC_data().
#' @param Y responses to perturb, before the observation noise is added.
#' @param rho_fac correlation between the scores of any two factors, at a given
#'   component. rho_fac = 0 leaves Y unchanged. Must leave the equicorrelation
#'   matrix positive definite, that is, rho_fac > -1 / (Q - 1).
#' @param seed seed used for the mixing draws only.
#'
#' @return a list with the perturbed Y and the correlated scores, the latter
#'   being the quantity the fitted model should be compared against.
add_latent_correlation <- function(fpca_data, Y, rho_fac, seed = NULL) {

  Zeta <- fpca_data$Zeta

  if (rho_fac == 0) {
    return(list(Y = Y, Zeta_corr = Zeta))
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  B <- fpca_data$B
  Q <- length(Zeta)
  L <- ncol(Zeta[[1]])
  N <- length(Y)
  p <- length(Y[[1]])

  stopifnot(Q > 1, rho_fac > -1 / (Q - 1), rho_fac < 1)

  R <- matrix(rho_fac, nrow = Q, ncol = Q)
  diag(R) <- 1
  C <- chol(R) # upper triangular, so that Z %*% C has row covariance R

  # scores of component l across factors, subjects in rows
  Zeta_corr <- Zeta
  for (l in seq_len(L)) {
    Z_l <- sapply(seq_len(Q), function(q) Zeta[[q]][, l]) # N x Q
    if (is.vector(Z_l)) { # single subject
      Z_l <- t(as.matrix(Z_l))
    }
    Z_l_corr <- Z_l %*% C
    for (q in seq_len(Q)) {
      Zeta_corr[[q]][, l] <- Z_l_corr[, q]
    }
  }

  # eigenfunction evaluations, computed once per factor
  list_Psi <- lapply(seq_len(Q), function(q) eval_eigenfunctions_obs(fpca_data, q))

  for (i in seq_len(N)) {

    n_i <- length(fpca_data$time_obs[[i]])

    # change in the latent trajectories of subject i, Q x n_i
    Delta <- matrix(0, nrow = Q, ncol = n_i)
    for (q in seq_len(Q)) {
      d_q <- Zeta_corr[[q]][i, ] - Zeta[[q]][i, ]
      Delta[q, ] <- as.vector(list_Psi[[q]][[i]] %*% d_q)
    }

    E <- B %*% Delta # p x n_i, the change reaching each variable

    Y_i <- Y[[i]]
    for (j in seq_len(p)) {
      Y_i[[j]] <- Y_i[[j]] + E[j, ]
    }
    Y[[i]] <- Y_i
  }

  create_named_list(Y, Zeta_corr)
}


## ---------------------------------------------------------------------------
## Scenario 2: informative drop-out
## ---------------------------------------------------------------------------

#' Latent factor process evaluated at the observation times
#'
#' Used as the driver of the drop-out hazard. If the generator returned the
#' subject-specific eigenfunction evaluations ($Phi), they are used directly;
#' otherwise the grid evaluations ($Phi_g) are interpolated, which is accurate
#' to numerical precision for smooth splines on a grid of 1000 points.
#'
#' @param fpca_data output of generate_bayesSYNC_data().
#' @param q_dropout index of the latent factor driving drop-out.
#'
#' @return a list of length N of numeric vectors h_i^{(q)}(t_i1), ..., h_i^{(q)}(t_ini).
compute_factor_path_obs <- function(fpca_data, q_dropout = 1) {
  
  N <- length(fpca_data$time_obs)
  Zeta_q <- fpca_data$Zeta[[q_dropout]]
  
  if (!is.null(fpca_data$Phi)) {
    
    lapply(seq_len(N), function(i) {
      as.vector(fpca_data$Phi[[q_dropout]][[i]] %*% Zeta_q[i, ])
    })
    
  } else {
    
    time_g <- fpca_data$time_g
    Phi_g_q <- fpca_data$Phi_g[[q_dropout]]
    L <- ncol(Phi_g_q)
    
    lapply(seq_len(N), function(i) {
      t_i <- fpca_data$time_obs[[i]]
      eval_i <- sapply(seq_len(L), function(l) {
        approx(time_g, Phi_g_q[, l], xout = t_i, rule = 2)$y
      })
      if (is.vector(eval_i)) { # single observation time
        eval_i <- t(as.matrix(eval_i))
      }
      as.vector(eval_i %*% Zeta_q[i, ])
    })
  }
}


#' Expected proportion of retained observations under a given drop-out hazard
#'
#' The drop-out hazard at visit k is expit(alpha + beta * s_ik), with the
#' baseline visit always retained. Once a subject drops out, all subsequent
#' visits are lost, so the probability of observing visit k is the product of
#' the survival probabilities up to k. The expectation below is exact given the
#' realised latent paths, which allows the intercept to be calibrated without
#' Monte Carlo error.
#'
#' @param list_s list of standardised latent values at the observation times.
#' @param alpha intercept of the drop-out hazard.
#' @param beta slope of the drop-out hazard, per standard deviation of the
#'   latent factor. beta = 0 gives non-informative (MCAR) drop-out.
#'
#' @return the expected proportion of retained observations.
expected_retention <- function(list_s, alpha, beta) {
  
  num <- 0
  den <- 0
  
  for (i in seq_along(list_s)) {
    n_i <- length(list_s[[i]])
    haz <- plogis(alpha + beta * list_s[[i]])
    haz[1] <- 0 # baseline visit always retained
    num <- num + sum(cumprod(1 - haz))
    den <- den + n_i
  }
  
  num / den
}


#' Calibrate the drop-out intercept to a target retention
#'
#' Called separately for each drop-out arm, so that the informative and the
#' non-informative arms discard the same expected amount of data and any
#' difference in performance is attributable to the informativeness of the
#' mechanism rather than to the sample size.
#'
#' @param list_s list of standardised latent values at the observation times.
#' @param beta slope of the drop-out hazard.
#' @param target_retention target expected proportion of retained observations.
#' @param interval search interval for the intercept.
#'
#' @return the calibrated intercept.
calibrate_dropout_intercept <- function(list_s, beta, target_retention,
                                        interval = c(-20, 20)) {
  
  f <- function(a) expected_retention(list_s, a, beta) - target_retention
  
  # expected_retention is decreasing in alpha
  if (f(interval[1]) < 0) {
    warning("Target retention is above the attainable range; using the lower bound.")
    return(interval[1])
  }
  if (f(interval[2]) > 0) {
    warning("Target retention is below the attainable range; using the upper bound.")
    return(interval[2])
  }
  
  uniroot(f, interval = interval)$root
}


#' Apply drop-out to a simulated dataset
#'
#' The hazard at visit k depends on the value of the latent factor at visit k,
#' that is, on the value that would have been recorded had the subject not
#' dropped out. The mechanism is therefore missing not at random. Truncation is
#' monotone: all visits from the drop-out time onwards are removed.
#'
#' @param fpca_data output of generate_bayesSYNC_data().
#' @param Y list of responses to be truncated, possibly already carrying the
#'   error process added by add_correlated_noise().
#' @param list_s list of standardised latent values at the observation times.
#' @param alpha intercept of the drop-out hazard.
#' @param beta slope of the drop-out hazard.
#' @param min_obs minimum number of observations retained per subject. Subjects
#'   drawn to fall below this are truncated at min_obs instead, which keeps the
#'   fit well posed at the price of a slightly higher realised retention. The
#'   default of 1 leaves the mechanism undistorted, since bayesSYNC borrows
#'   information across subjects and accepts single-observation subjects.
#' @param seed seed used for the drop-out draws only.
#'
#' @return a list with the truncated time_obs and Y, the drop-out indicators and
#'   summaries of the realised truncation.
apply_dropout <- function(fpca_data, Y, list_s, alpha, beta,
                          min_obs = 1, seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  time_obs <- fpca_data$time_obs
  N <- length(time_obs)
  p <- length(Y[[1]])
  
  n_obs_full <- sapply(time_obs, length)
  n_obs_new <- rep(NA, N)
  vec_time_dropout <- rep(NA, N) # time of the first missed visit, NA for completers
  
  for (i in seq_len(N)) {
    
    t_i <- time_obs[[i]]
    n_i <- length(t_i)
    
    haz <- plogis(alpha + beta * list_s[[i]])
    haz[1] <- 0 # baseline visit always retained
    
    drawn <- runif(n_i) < haz
    k_drop <- which(drawn)[1] # index of the first missed visit, NA if none
    
    if (is.na(k_drop)) {
      n_keep <- n_i
    } else {
      # the minimum number of observations overrides the drawn drop-out time,
      # which can only happen for subjects with very few scheduled visits
      n_keep <- max(k_drop - 1, min(min_obs, n_i))
      if (n_keep < n_i) {
        vec_time_dropout[i] <- t_i[n_keep + 1] # first time point actually lost
      }
    }
    
    if (n_keep < n_i) {
      ind_keep <- seq_len(n_keep)
      time_obs[[i]] <- t_i[ind_keep]
      Y_i <- Y[[i]]
      for (j in seq_len(p)) {
        Y_i[[j]] <- Y_i[[j]][ind_keep]
      }
      Y[[i]] <- Y_i
    }
    
    n_obs_new[i] <- n_keep
  }
  
  vec_dropout <- n_obs_new < n_obs_full # subjects with at least one lost visit
  retention <- sum(n_obs_new) / sum(n_obs_full)
  
  create_named_list(time_obs, Y, n_obs_full, n_obs_new, vec_dropout,
                    vec_time_dropout, retention, alpha, beta)
}


## ---------------------------------------------------------------------------
## Diagnostics on the fitted error process
## ---------------------------------------------------------------------------

#' Variable-specific error variances of a fitted model
#'
#' bayesSYNC_core() computes sigsq_eps as the posterior mean of the
#' inverse-gamma variational posterior, but the entry is currently commented out
#' of the returned list. It is recovered here from quantities that are returned,
#' using the construction of the prediction bands in orthonormalise():
#'
#'   list_Y_upp[[i]][[j]] - list_Y_hat[[i]][[j]]
#'     = qnorm(0.975) * sqrt(sum_q B_hat[j,q]^2 * list_var_vec[[q]][[i]] + sigsq_eps[j]),
#'
#' where B_hat is mu_q_b and list_var_vec is returned as well. The identity holds
#' at every grid point, so the median across the grid is taken for stability.
#' If sigsq_eps is added back to the returned list, it is used directly.
#'
#' The quantities are on the scale of the data actually fitted, that is, the
#' standardised scale when bool_scale = TRUE.
#'
#' @param fpca_res output of bayesSYNC().
#' @param i_ref index of the subject used for the reconstruction; any subject
#'   gives the same answer up to numerical error.
#'
#' @return a vector of length p of error variances.
extract_sigsq_eps <- function(fpca_res, i_ref = 1) {
  
  if (!is.null(fpca_res$sigsq_eps)) {
    return(fpca_res$sigsq_eps)
  }
  
  B_hat <- fpca_res$B_hat # unpermuted, matching list_var_vec
  Q <- ncol(B_hat)
  p <- nrow(B_hat)
  
  var_latent_i <- lapply(seq_len(Q), function(q) fpca_res$list_var_vec[[q]][[i_ref]])
  
  sapply(seq_len(p), function(j) {
    sd_i_j <- (fpca_res$list_Y_upp[[i_ref]][[j]] - fpca_res$list_Y_hat[[i_ref]][[j]]) / qnorm(0.975)
    var_signal_j <- Reduce("+", lapply(seq_len(Q), function(q) B_hat[j, q]^2 * var_latent_i[[q]]))
    median(sd_i_j^2 - var_signal_j)
  })
}


#' Empirical correlogram of the residuals along t
#'
#' Direct check of the white-noise assumption: under a correctly specified error
#' process the residual correlation is zero at every non-zero lag, whereas
#' correlated errors leave a positive and decaying correlation that the fit
#' cannot absorb. The same function can be applied to a fit on the real data, so
#' the assumption can be assessed empirically rather than only argued.
#'
#' Fitted values are returned on the dense grid; each observation time is mapped
#' to its nearest grid point, whose spacing is 1e-3 for the default n_g = 1000,
#' so the induced error is negligible for smooth fits. All variables of a
#' subject share one time grid, so the sums of residual cross-products over
#' variables are obtained by a single crossproduct per subject.
#'
#' @param fpca_res output of bayesSYNC(). Its $Y entry is used, so that the
#'   residuals are on the same scale as the fitted values whatever bool_scale.
#' @param time_obs list of observation times, as passed to bayesSYNC().
#' @param breaks bin edges for the lag between two observation times.
#' @param max_p optional number of variables to subsample, for speed.
#' @param seed seed for the subsampling of variables.
#'
#' @return a data frame with one row per lag bin, giving the number of
#'   contributing pairs and the estimated residual correlation.
residual_correlogram <- function(fpca_res, time_obs,
                                 breaks = seq(0, 1, by = 0.05),
                                 max_p = NULL, seed = NULL) {
  
  Y <- fpca_res$Y # possibly standardised, consistently with list_Y_hat
  list_Y_hat <- fpca_res$list_Y_hat
  time_g <- fpca_res$time_g
  n_g <- length(time_g)
  delta_g <- time_g[2] - time_g[1]
  
  N <- length(Y)
  p_all <- length(Y[[1]])
  
  if (!is.null(max_p) && max_p < p_all) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    ind_var <- sort(sample(seq_len(p_all), max_p))
  } else {
    ind_var <- seq_len(p_all)
  }
  p <- length(ind_var)
  
  n_bin <- length(breaks) - 1
  sum_prod <- rep(0, n_bin)
  n_pair <- rep(0, n_bin)
  sum_sq <- 0
  n_sq <- 0
  
  for (i in seq_len(N)) {
    
    t_i <- time_obs[[i]]
    n_i <- length(t_i)
    if (n_i < 2) {
      next # no pair contributes, but the variance term is still informative
    }
    
    ind_g <- pmin(pmax(round((t_i - time_g[1]) / delta_g) + 1, 1), n_g)
    
    R_i <- matrix(unlist(lapply(ind_var, function(j) {
      Y[[i]][[j]] - list_Y_hat[[i]][[j]][ind_g]
    })), nrow = p, ncol = n_i, byrow = TRUE)
    
    S_i <- crossprod(R_i) # entry (k, l) is the sum over variables of r_k r_l
    
    sum_sq <- sum_sq + sum(diag(S_i))
    n_sq <- n_sq + p * n_i
    
    d_i <- abs(outer(t_i, t_i, "-"))
    ind_upper <- upper.tri(d_i)
    bin_i <- cut(d_i[ind_upper], breaks = breaks, include.lowest = TRUE, labels = FALSE)
    val_i <- S_i[ind_upper]
    
    for (b in seq_len(n_bin)) {
      sel <- which(bin_i == b)
      if (length(sel) > 0) {
        sum_prod[b] <- sum_prod[b] + sum(val_i[sel])
        n_pair[b] <- n_pair[b] + length(sel) * p
      }
    }
  }
  
  var_resid <- sum_sq / n_sq
  
  data.frame(dist_low = breaks[-length(breaks)],
             dist_upp = breaks[-1],
             dist_mid = (breaks[-length(breaks)] + breaks[-1]) / 2,
             n_pair = n_pair,
             corr = ifelse(n_pair > 0, (sum_prod / pmax(n_pair, 1)) / var_resid, NA),
             var_resid = var_resid,
             stringsAsFactors = FALSE)
}


#' Check that the component alignment behaved as a permutation to the identity
#'
#' match_sign_components() matches each simulated component (q, l) to the
#' estimated component maximising the absolute score correlation over all pairs
#' (q_tilde, l_tilde), so it can in principle pull a component from a different
#' factor. The posterior covariances are permuted by factor only, not by
#' component, so a cross-component match would leave the coverage metric indexed
#' against the wrong variance. In the well-specified runs the match is the
#' identity and the question does not arise, but it is worth verifying under
#' misspecification.
#'
#' The check below is a necessary condition rather than a proof: it verifies
#' that, after alignment, each simulated component is best matched by the
#' estimated component carrying the same indices.
#'
#' @param Zeta simulated scores.
#' @param list_Zeta_hat aligned estimated scores.
#'
#' @return TRUE when the alignment is consistent with the identity for every
#'   component, FALSE otherwise.
alignment_diagnostic <- function(Zeta, list_Zeta_hat) {
  
  Q_true <- length(Zeta)
  L_true <- ncol(Zeta[[1]])
  
  ok <- matrix(NA, nrow = Q_true, ncol = L_true)
  
  for (q in seq_len(Q_true)) {
    for (l in seq_len(L_true)) {
      
      corr_all <- sapply(seq_len(L_true), function(l_tilde) {
        sapply(seq_len(Q_true), function(q_tilde) {
          abs(cor(Zeta[[q]][, l], list_Zeta_hat[[q_tilde]][, l_tilde]))
        })
      })
      corr_all <- matrix(corr_all, nrow = Q_true, ncol = L_true)
      corr_all[is.na(corr_all)] <- 0
      
      ind_max <- which(corr_all == max(corr_all), arr.ind = TRUE)[1, ]
      ok[q, l] <- (ind_max[1] == q) && (ind_max[2] == l)
    }
  }
  
  all(ok)
}


## ---------------------------------------------------------------------------
## Metrics specific to the misspecification study
## ---------------------------------------------------------------------------

#' Bias, attenuation and coverage of the subject-level scores
#'
#' Complements the score correlations already reported elsewhere. Correlation is
#' invariant to shrinkage, so it can stay high while the scores are attenuated
#' towards zero, which is precisely the failure mode expected under informative
#' drop-out. The regression slope of the estimated on the true score captures
#' that attenuation, and the coverage of the marginal 95% credible intervals
#' captures whether the variational posterior remains honest.
#'
#' Scores and covariances are assumed to have been aligned and rescaled by
#' match_factor_and_sign(), so that they are on the scale of the simulated ones.
#'
#' @param Zeta simulated scores, list of length Q_true of N x L_true matrices.
#' @param list_Zeta_hat estimated scores, same structure.
#' @param list_Cov_zeta_hat list of length Q of lists of length N of L x L
#'   posterior covariance matrices. Coverage is skipped if NULL.
#' @param group optional logical vector of length N defining two subgroups
#'   (used here for droppers versus completers), for stratified summaries.
#'
#' @return a list of Q_true x L_true matrices of bias, slope and coverage, plus
#'   the stratified biases when group is supplied.
score_bias_metrics <- function(Zeta, list_Zeta_hat, list_Cov_zeta_hat = NULL,
                               group = NULL) {
  
  Q_true <- length(Zeta)
  L_true <- ncol(Zeta[[1]])
  N <- nrow(Zeta[[1]])
  
  bias <- slope <- coverage <- matrix(NA, nrow = Q_true, ncol = L_true)
  bias_group_1 <- bias_group_0 <- matrix(NA, nrow = Q_true, ncol = L_true)
  
  # the stratified summaries are only defined when both groups are non-empty,
  # which is not the case in the settings without drop-out
  if (!is.null(group) && (sum(group) == 0 || sum(!group) == 0)) {
    group <- NULL
  }
  
  for (q in seq_len(Q_true)) {
    for (l in seq_len(L_true)) {
      
      zeta_true <- Zeta[[q]][, l]
      zeta_hat <- list_Zeta_hat[[q]][, l]
      
      bias[q, l] <- mean(zeta_hat - zeta_true)
      slope[q, l] <- unname(coef(lm(zeta_hat ~ zeta_true))[2])
      
      if (!is.null(list_Cov_zeta_hat)) {
        sd_ql <- sapply(seq_len(N), function(i) sqrt(list_Cov_zeta_hat[[q]][[i]][l, l]))
        coverage[q, l] <- mean(abs(zeta_hat - zeta_true) <= qnorm(0.975) * sd_ql)
      }
      
      if (!is.null(group)) {
        bias_group_1[q, l] <- mean((zeta_hat - zeta_true)[group])
        bias_group_0[q, l] <- mean((zeta_hat - zeta_true)[!group])
      }
    }
  }
  
  create_named_list(bias, slope, coverage, bias_group_1, bias_group_0)
}


#' Flatten the metrics of one run into a single-row data frame
#'
#' Keeps the gathering step across cluster jobs to a plain rbind.
#'
#' @return a one-row data frame of scalar summaries.
build_summary_row <- function(setting, seed, p, Corr_zeta, Corr_b, Auc_b,
                              Err_mean, Err_eigen, n_retained, score_metrics,
                              retention, mean_n_obs, prop_dropout, rt,
                              sigsq_eps_hat = NA, df_correlogram = NULL,
                              alignment_ok = NA, min_norm_col_B_hat = NA) {
  
  Q_true <- nrow(Corr_zeta)
  L_true <- ncol(Corr_zeta)
  
  flatten <- function(mat, prefix) {
    out <- as.vector(t(mat))
    names(out) <- paste0(prefix, "_q", rep(seq_len(nrow(mat)), each = ncol(mat)),
                         "_l", rep(seq_len(ncol(mat)), times = nrow(mat)))
    out
  }
  
  # residual correlation in the first three non-zero lag bins, which is where a
  # correlated error process is most visible given the observation density
  if (!is.null(df_correlogram)) {
    corr_resid <- df_correlogram$corr[1:min(3, nrow(df_correlogram))]
    names(corr_resid) <- paste0("corr_resid_bin", seq_along(corr_resid))
  } else {
    corr_resid <- c(corr_resid_bin1 = NA, corr_resid_bin2 = NA, corr_resid_bin3 = NA)
  }
  
  vec <- c(auc_b = mean(Auc_b),
           corr_b = mean(Corr_b),
           ise_mean_per_variable = Err_mean / p, # error_mean_function() returns a total over the p variables
           ise_eigen_avg = mean(Err_eigen),
           corr_zeta_avg = mean(Corr_zeta),
           n_retained = n_retained,
           coverage_zeta_avg = mean(score_metrics$coverage, na.rm = TRUE),
           slope_zeta_avg = mean(score_metrics$slope),
           bias_zeta_avg = mean(score_metrics$bias),
           sigsq_eps_mean = mean(sigsq_eps_hat),
           sigsq_eps_median = median(sigsq_eps_hat),
           alignment_ok = as.numeric(alignment_ok),
           min_norm_col_B_hat = min_norm_col_B_hat,
           retention = retention,
           mean_n_obs = mean_n_obs,
           prop_dropout = prop_dropout,
           runtime_min = rt,
           corr_resid,
           flatten(Corr_zeta, "corr_zeta"),
           flatten(Err_eigen, "ise_eigen"),
           flatten(score_metrics$slope, "slope_zeta"),
           flatten(score_metrics$coverage, "coverage_zeta"),
           flatten(score_metrics$bias_group_1, "bias_zeta_droppers"),
           flatten(score_metrics$bias_group_0, "bias_zeta_completers"))
  
  setting_row <- setting
  rownames(setting_row) <- NULL
  
  cbind(setting_row, seed = seed, as.data.frame(t(vec)))
}