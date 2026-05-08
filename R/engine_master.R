
###########################
### Helper functions.   ###
###########################

# Helper: polychoric-like Pearson correlation of two ordinal variables
# created by thresholding standard normals with correlation r.
ordinal_corr_from_latent <- function(r, thrA, thrB, scoresA = NULL, scoresB = NULL) {
  KA <- length(thrA) + 1L
  KB <- length(thrB) + 1L
  if (is.null(scoresA)) scoresA <- seq_len(KA)
  if (is.null(scoresB)) scoresB <- seq_len(KB)
  cutsA <- c(-Inf, thrA, Inf)
  cutsB <- c(-Inf, thrB, Inf)
  Sigma <- matrix(c(1, r, r, 1), 2, 2)

  P <- matrix(0, KA, KB)
  for (i in 1:KA) for (j in 1:KB) {
    lower <- c(cutsA[i], cutsB[j]); upper <- c(cutsA[i+1], cutsB[j+1])
    P[i, j] <- mvtnorm::pmvnorm(lower = lower, upper = upper, mean = c(0,0), sigma = Sigma)
  }
  P <- P / sum(P)
  pA <- rowSums(P); pB <- colSums(P)

  EA <- sum(pA * scoresA); EB <- sum(pB * scoresB)
  VA <- sum(pA * (scoresA - EA)^2); VB <- sum(pB * (scoresB - EB)^2)
  EAB <- sum(P * outer(scoresA, scoresB))
  (EAB - EA * EB) / sqrt(VA * VB)
}

# Compute latent *item-level* correlation between Z_Ai and Z_Bj given the measurement model
latent_item_corr_AB <- function(rho, la, lb, ca, cb, theta_ai, theta_bj) {
  # Item models (standard in your setup):
  # Z_Ai = la * A + ca * B + e_ai
  # Z_Bj = cb * A + lb * B + e_bj
  # Var(A)=Var(B)=1, Cov(A,B)=rho, e's independent of A,B and each other
  cov_ij <- la * cb + ca * lb + rho * (la * lb + ca * cb)
  var_i  <- la^2 + ca^2 + 2 * rho * la * ca + theta_ai
  var_j  <- lb^2 + cb^2 + 2 * rho * lb * cb + theta_bj
  cov_ij / sqrt(var_i * var_j)
}

# Compute latent item-item correlation within A (or within B)
latent_item_corr_within <- function(rho, l1, l2, c1, c2, theta1, theta2) {
  # Z1 = l1*A + c1*B + e1, Z2 = l2*A + c2*B + e2
  cov_12 <- l1 * l2 + c1 * c2 + rho * (l1 * c2 + c1 * l2)
  var1   <- l1^2 + c1^2 + 2 * rho * l1 * c1 + theta1
  var2   <- l2^2 + c2^2 + 2 * rho * l2 * c2 + theta2
  cov_12 / sqrt(var1 * var2)
}


### Analytic confidence intervals ###
compute_ci_analytic <- function(
    r_obs,
    r_SE = NULL,
    r_ci = NULL,
    conf_level = 0.95,
    correct_fun
) {
  if (is.null(r_SE) && is.null(r_ci)) {
    stop("For analytic CI, provide either r_SE or r_ci.")
  }

  if (!is.null(r_ci)) {
    if (!is.numeric(r_ci) || length(r_ci) != 2L || any(!is.finite(r_ci))) {
      stop("r_ci must be a numeric vector of length 2: c(lower, upper).")
    }
    if (any(abs(r_ci) >= 1)) {
      stop("Both values in r_ci must be strictly between -1 and 1.")
    }

    r_low <- min(r_ci)
    r_up  <- max(r_ci)

  } else {
    alpha <- 1 - conf_level

    z_obs <- atanh(r_obs)
    z_se  <- r_SE / (1 - r_obs^2)

    z_low <- z_obs + qnorm(alpha / 2) * z_se
    z_up  <- z_obs + qnorm(1 - alpha / 2) * z_se

    r_low <- tanh(z_low)
    r_up  <- tanh(z_up)
  }

  rho_low <- correct_fun(r_low)
  rho_up  <- correct_fun(r_up)

  CI_lower <- min(rho_low, rho_up, na.rm = TRUE)
  CI_upper <- max(rho_low, rho_up, na.rm = TRUE)

  list(
    CI_lower = CI_lower,
    CI_upper = CI_upper,
    r_CI_lower_raw = r_low,
    r_CI_upper_raw = r_up,
    method = if (!is.null(r_ci)) {
      "analytic_observed_r_ci_endpoint"
    } else {
      "analytic_fisher_z_endpoint"
    }
  )
}

########################
### Continuous items ###
########################

r_pred_general_cont_items <- function(
    rho,
    lambdaA, lambdaB,          # vectors (preferred) or scalars (replicated)
    thetaA = NULL, thetaB = NULL,   # now robust to NULL / length 0
    cA_vec = NULL, cB_vec = NULL,
    cA = 0, cB = 0, overlapA = NULL, overlapB = NULL, k = 0
){
  # -----------------------------------------
  # Resolve lambda lengths
  # -----------------------------------------
  nA <- length(lambdaA)
  nB <- length(lambdaB)

  # -----------------------------------------
  # Automatic residual variances if missing/empty
  # -----------------------------------------
  if (is.null(thetaA) || length(thetaA) == 0L) {
    thetaA <- 1 - lambdaA^2
  } else if (length(thetaA) == 1L) {
    thetaA <- rep(thetaA, nA)
  } else if (length(thetaA) != nA) {
    stop("thetaA must be NULL/length 0 (auto), length 1, or same length as lambdaA.")
  }

  if (is.null(thetaB) || length(thetaB) == 0L) {
    thetaB <- 1 - lambdaB^2
  } else if (length(thetaB) == 1L) {
    thetaB <- rep(thetaB, nB)
  } else if (length(thetaB) != nB) {
    stop("thetaB must be NULL/length 0 (auto), length 1, or same length as lambdaB.")
  }

  # -----------------------------------------
  # Cross-loading vectors (robust to NULL / length 0 / all NA)
  # -----------------------------------------
  if (is.null(cA_vec) || length(cA_vec) == 0L || all(is.na(cA_vec))) {
    cA_vec <- rep(0, nA)
    if (!is.null(overlapA) && length(overlapA) > 0L) {
      cA_vec[overlapA] <- cA
    } else if (!is.null(k) && k > 0L) {
      cA_vec[(nA - k + 1L):nA] <- cA
    }
  } else if (length(cA_vec) != nA) {
    stop("cA_vec must be NULL/length 0 (auto), or the same length as lambdaA.")
  }

  if (is.null(cB_vec) || length(cB_vec) == 0L || all(is.na(cB_vec))) {
    cB_vec <- rep(0, nB)
    if (!is.null(overlapB) && length(overlapB) > 0L) {
      cB_vec[overlapB] <- cB
    } else if (!is.null(k) && k > 0L) {
      cB_vec[(nB - k + 1L):nB] <- cB
    }
  } else if (length(cB_vec) != nB) {
    stop("cB_vec must be NULL/length 0 (auto), or the same length as lambdaB.")
  }

  # -----------------------------------------
  # Cross-scale covariance
  # -----------------------------------------
  covAB <- 0
  for (i in 1:nA) for (j in 1:nB) {
    covAB <- covAB + (
      lambdaA[i] * lambdaB[j] * rho +
        lambdaA[i] * cB_vec[j] +
        cA_vec[i] * lambdaB[j] +
        cA_vec[i] * cB_vec[j] * rho
    )
  }

  # -----------------------------------------
  # Within-scale variances
  # -----------------------------------------
  varA_sum <- 0
  for (i in 1:nA) for (j in 1:nA) {
    if (i == j) {
      varA_sum <- varA_sum + (lambdaA[i]^2 + cA_vec[i]^2 +
                                2 * rho * lambdaA[i] * cA_vec[i] + thetaA[i])
    } else {
      varA_sum <- varA_sum + (
        lambdaA[i] * lambdaA[j] +
          cA_vec[i] * cA_vec[j] +
          rho * (lambdaA[i] * cA_vec[j] + cA_vec[i] * lambdaA[j])
      )
    }
  }

  varB_sum <- 0
  for (i in 1:nB) for (j in 1:nB) {
    if (i == j) {
      varB_sum <- varB_sum + (lambdaB[i]^2 + cB_vec[i]^2 +
                                2 * rho * lambdaB[i] * cB_vec[i] + thetaB[i])
    } else {
      varB_sum <- varB_sum + (
        lambdaB[i] * lambdaB[j] +
          cB_vec[i] * cB_vec[j] +
          rho * (lambdaB[i] * cB_vec[j] + cB_vec[i] * lambdaB[j])
      )
    }
  }

  covAB / sqrt(varA_sum * varB_sum)
}


## ---------- Numeric inversion for continuous model ----------
rho_correction_continuous <- function(
    r_obs,
    lambdaA, lambdaB,
    thetaA = NULL, thetaB = NULL,
    cA_vec = NULL, cB_vec = NULL,
    cA = 0, cB = 0,
    overlapA = NULL, overlapB = NULL,
    k = 0,
    lower = -0.999, upper = 0.999,
    tol = 1e-8, maxiter = 100
){

  pred_r <- function(rho) {
    r_pred_general_cont_items(
      rho = rho,
      lambdaA = lambdaA,
      lambdaB = lambdaB,
      thetaA  = thetaA,
      thetaB  = thetaB,
      cA_vec = cA_vec,
      cB_vec = cB_vec,
      cA = cA,
      cB = cB,
      overlapA = overlapA,
      overlapB = overlapB,
      k = k
    )
  }

  make_result <- function(rho_hat,
                          boundary = FALSE,
                          boundary_side = NA_character_,
                          boundary_reason = NA_character_,
                          feasible_r_min = NA_real_,
                          feasible_r_max = NA_real_) {
    list(
      rho_hat = rho_hat,
      boundary = boundary,
      boundary_side = boundary_side,
      boundary_reason = boundary_reason,
      feasible_r_min = feasible_r_min,
      feasible_r_max = feasible_r_max
    )
  }

  f <- function(rho) pred_r(rho) - r_obs

  ## Feasible observed range implied by the model
  r_min <- pred_r(-1)
  r_max <- pred_r(1)

  if (!is.finite(r_min) || !is.finite(r_max)) {
    feasible_r_min <- NA_real_
    feasible_r_max <- NA_real_
  } else {
    feasible_r_min <- min(r_min, r_max)
    feasible_r_max <- max(r_min, r_max)
  }

  ## Immediate boundary handling
  if (is.finite(feasible_r_min) && r_obs < feasible_r_min) {
    return(make_result(
      rho_hat = -1,
      boundary = TRUE,
      boundary_side = "lower",
      boundary_reason = "r_obs_below_feasible_min",
      feasible_r_min = feasible_r_min,
      feasible_r_max = feasible_r_max
    ))
  }

  if (is.finite(feasible_r_max) && r_obs > feasible_r_max) {
    return(make_result(
      rho_hat = 1,
      boundary = TRUE,
      boundary_side = "upper",
      boundary_reason = "r_obs_above_feasible_max",
      feasible_r_min = feasible_r_min,
      feasible_r_max = feasible_r_max
    ))
  }

  ## Try requested bracket first
  fL <- f(lower)
  fU <- f(upper)

  if (!is.finite(fL) || !is.finite(fU) || fL * fU > 0) {
    grid <- seq(-0.999, 0.999, length.out = 201)
    vals <- sapply(grid, f)

    ok <- is.finite(vals[-length(vals)]) & is.finite(vals[-1L])
    idx <- which(ok & (sign(vals[-length(vals)]) != sign(vals[-1L])))

    if (!length(idx)) {
      ## In principle this should now be rare because we already checked feasible range.
      ## Still keep a safe fallback.
      if (f(1) < 0) {
        return(make_result(
          rho_hat = 1,
          boundary = TRUE,
          boundary_side = "upper",
          boundary_reason = "no_root_found_clamped_upper",
          feasible_r_min = feasible_r_min,
          feasible_r_max = feasible_r_max
        ))
      } else {
        return(make_result(
          rho_hat = -1,
          boundary = TRUE,
          boundary_side = "lower",
          boundary_reason = "no_root_found_clamped_lower",
          feasible_r_min = feasible_r_min,
          feasible_r_max = feasible_r_max
        ))
      }
    }

    lower <- grid[idx[1]]
    upper <- grid[idx[1] + 1]
  }

  ## Root finding
  root_val <- tryCatch(
    uniroot(f, lower = lower, upper = upper, tol = tol, maxiter = maxiter)$root,
    error = function(e) NA_real_
  )

  if (is.finite(root_val)) {
    return(make_result(
      rho_hat = root_val,
      boundary = FALSE,
      boundary_side = NA_character_,
      boundary_reason = NA_character_,
      feasible_r_min = feasible_r_min,
      feasible_r_max = feasible_r_max
    ))
  }

  ## Final fallback
  if (f(1) < 0) {
    return(make_result(
      rho_hat = 1,
      boundary = TRUE,
      boundary_side = "upper",
      boundary_reason = "uniroot_failed_clamped_upper",
      feasible_r_min = feasible_r_min,
      feasible_r_max = feasible_r_max
    ))
  } else {
    return(make_result(
      rho_hat = -1,
      boundary = TRUE,
      boundary_side = "lower",
      boundary_reason = "uniroot_failed_clamped_lower",
      feasible_r_min = feasible_r_min,
      feasible_r_max = feasible_r_max
    ))
  }
}


###########################
### Ordinal item scores ###
###########################

r_pred_general_ordinal_items <- function(
    rho,
    lambdaA, lambdaB,          # vectors (preferred) or scalars (replicated)
    thetaA = NULL,  thetaB = NULL,   # now robust to NULL / length 0
    cA_vec = NULL, cB_vec = NULL,
    cA = 0, cB = 0, overlapA = NULL, overlapB = NULL, k = 0,
    thrA, thrB,                # thresholds supplied by user
    threshold_scale = c("z","raw"),
    scoresA = NULL, scoresB = NULL
){
  threshold_scale <- match.arg(threshold_scale)

  # Resolve lambda lengths
  nA <- length(lambdaA)
  nB <- length(lambdaB)
  if (nA < 1L || nB < 1L) stop("lambdaA and lambdaB must have length >= 1.")

  # Residual variances: allow NULL / length 0 -> 1 - lambda^2
  if (is.null(thetaA) || length(thetaA) == 0L) {
    thetaA <- 1 - lambdaA^2
  } else if (length(thetaA) == 1L) {
    thetaA <- rep(thetaA, nA)
  } else if (length(thetaA) != nA) {
    stop("thetaA must be NULL/length 0 (auto), length 1, or same length as lambdaA.")
  }

  if (is.null(thetaB) || length(thetaB) == 0L) {
    thetaB <- 1 - lambdaB^2
  } else if (length(thetaB) == 1L) {
    thetaB <- rep(thetaB, nB)
  } else if (length(thetaB) != nB) {
    stop("thetaB must be NULL/length 0 (auto), length 1, or same length as lambdaB.")
  }

  # Cross-loading vectors
  if (is.null(cA_vec) || length(cA_vec) == 0L || all(is.na(cA_vec))) {
    cA_vec <- rep(0, nA)
    if (!is.null(overlapA) && length(overlapA) > 0L) {
      cA_vec[overlapA] <- cA
    } else if (!is.null(k) && k > 0L) {
      cA_vec[(nA - k + 1L):nA] <- cA
    }
  } else if (length(cA_vec) != nA) {
    stop("cA_vec must be NULL/length 0 (auto), or same length as lambdaA.")
  }

  if (is.null(cB_vec) || length(cB_vec) == 0L || all(is.na(cB_vec))) {
    cB_vec <- rep(0, nB)
    if (!is.null(overlapB) && length(overlapB) > 0L) {
      cB_vec[overlapB] <- cB
    } else if (!is.null(k) && k > 0L) {
      cB_vec[(nB - k + 1L):nB] <- cB
    }
  } else if (length(cB_vec) != nB) {
    stop("cB_vec must be NULL/length 0 (auto), or same length as lambdaB.")
  }

  # Per-item SDs
  sdA <- sqrt(lambdaA^2 + cA_vec^2 + 2 * rho * lambdaA * cA_vec + thetaA)
  sdB <- sqrt(lambdaB^2 + cB_vec^2 + 2 * rho * lambdaB * cB_vec + thetaB)

  # Build per-item z-thresholds
  if (threshold_scale == "z") {
    thrA_z_list <- replicate(nA, thrA, simplify = FALSE)
    thrB_z_list <- replicate(nB, thrB, simplify = FALSE)
  } else { # "raw"
    thrA_z_list <- lapply(sdA, function(s) thrA / s)
    thrB_z_list <- lapply(sdB, function(s) thrB / s)
  }

  # Precompute marginal variances of ordinal scores
  VA_i <- numeric(nA)
  for (i in 1:nA) {
    thr_i <- thrA_z_list[[i]]
    KA <- length(thr_i) + 1L
    scA <- if (is.null(scoresA)) 1:KA else scoresA
    pA  <- diff(pnorm(c(-Inf, thr_i, Inf)))
    EA  <- sum(pA * scA)
    VA_i[i] <- sum(pA * (scA - EA)^2)
  }
  VB_j <- numeric(nB)
  for (j in 1:nB) {
    thr_j <- thrB_z_list[[j]]
    KB <- length(thr_j) + 1L
    scB <- if (is.null(scoresB)) 1:KB else scoresB
    pB  <- diff(pnorm(c(-Inf, thr_j, Inf)))
    EB  <- sum(pB * scB)
    VB_j[j] <- sum(pB * (scB - EB)^2)
  }

  # Between-scale covariance of ordinal sums
  covAB <- 0
  for (i in 1:nA) for (j in 1:nB) {
    r_lat <- latent_item_corr_AB(rho, lambdaA[i], lambdaB[j],
                                 cA_vec[i], cB_vec[j],
                                 thetaA[i], thetaB[j])
    r_ord <- ordinal_corr_from_latent(r_lat, thrA_z_list[[i]], thrB_z_list[[j]],
                                      scoresA, scoresB)
    covAB <- covAB + r_ord * sqrt(VA_i[i] * VB_j[j])
  }

  # Within-scale variances
  varA_sum <- sum(VA_i)  # diag
  if (nA > 1) {
    for (i in 1:(nA-1)) for (j in (i+1):nA) {
      r_lat <- latent_item_corr_within(rho, lambdaA[i], lambdaA[j],
                                       cA_vec[i], cA_vec[j],
                                       thetaA[i], thetaA[j])
      r_ord <- ordinal_corr_from_latent(r_lat, thrA_z_list[[i]], thrA_z_list[[j]],
                                        scoresA, scoresA)
      varA_sum <- varA_sum + 2 * r_ord * sqrt(VA_i[i] * VA_i[j])
    }
  }

  varB_sum <- sum(VB_j)
  if (nB > 1) {
    for (i in 1:(nB-1)) for (j in (i+1):nB) {
      r_lat <- latent_item_corr_within(rho, lambdaB[i], lambdaB[j],
                                       cB_vec[i], cB_vec[j],
                                       thetaB[i], thetaB[j])
      r_ord <- ordinal_corr_from_latent(r_lat, thrB_z_list[[i]], thrB_z_list[[j]],
                                        scoresB, scoresB)
      varB_sum <- varB_sum + 2 * r_ord * sqrt(VB_j[i] * VB_j[j])
    }
  }

  covAB / sqrt(varA_sum * varB_sum)
}



## ---------- Numeric inversion for ordinal model ----------
rho_correction_ordinal <- function(
    r_obs,
    # measurement model parameters
    lambdaA, lambdaB,
    thetaA,  thetaB,
    cA_vec = NULL, cB_vec = NULL,
    cA = 0, cB = 0,
    overlapA = NULL, overlapB = NULL,
    k = 0,
    thrA, thrB,
    threshold_scale = c("z", "raw"),
    scoresA = NULL, scoresB = NULL,
    # numeric inversion controls
    lower = -0.999, upper = 0.999,
    tol = 1e-6, maxiter = 100
){
  threshold_scale <- match.arg(threshold_scale)

  pred_r <- function(rho) {
    r_pred_general_ordinal_items(
      rho             = rho,
      lambdaA         = lambdaA,
      lambdaB         = lambdaB,
      thetaA          = thetaA,
      thetaB          = thetaB,
      cA_vec          = cA_vec,
      cB_vec          = cB_vec,
      cA              = cA,
      cB              = cB,
      overlapA        = overlapA,
      overlapB        = overlapB,
      k               = k,
      thrA            = thrA,
      thrB            = thrB,
      threshold_scale = threshold_scale,
      scoresA         = scoresA,
      scoresB         = scoresB
    )
  }

  make_result <- function(rho_hat,
                          boundary = FALSE,
                          boundary_side = NA_character_,
                          boundary_reason = NA_character_,
                          feasible_r_min = NA_real_,
                          feasible_r_max = NA_real_) {
    list(
      rho_hat = rho_hat,
      boundary = boundary,
      boundary_side = boundary_side,
      boundary_reason = boundary_reason,
      feasible_r_min = feasible_r_min,
      feasible_r_max = feasible_r_max
    )
  }

  f <- function(rho) pred_r(rho) - r_obs

  ## Feasible observed range implied by the ordinal measurement model
  r_min <- pred_r(-1)
  r_max <- pred_r(1)

  if (!is.finite(r_min) || !is.finite(r_max)) {
    feasible_r_min <- NA_real_
    feasible_r_max <- NA_real_
  } else {
    feasible_r_min <- min(r_min, r_max)
    feasible_r_max <- max(r_min, r_max)
  }

  ## Immediate boundary handling
  if (is.finite(feasible_r_min) && r_obs < feasible_r_min) {
    return(make_result(
      rho_hat = -1,
      boundary = TRUE,
      boundary_side = "lower",
      boundary_reason = "r_obs_below_feasible_min",
      feasible_r_min = feasible_r_min,
      feasible_r_max = feasible_r_max
    ))
  }

  if (is.finite(feasible_r_max) && r_obs > feasible_r_max) {
    return(make_result(
      rho_hat = 1,
      boundary = TRUE,
      boundary_side = "upper",
      boundary_reason = "r_obs_above_feasible_max",
      feasible_r_min = feasible_r_min,
      feasible_r_max = feasible_r_max
    ))
  }

  ## Try requested bracket first
  fL <- f(lower)
  fU <- f(upper)

  ## If no sign change, search on a grid
  if (!is.finite(fL) || !is.finite(fU) || fL * fU > 0) {
    grid <- seq(-0.999, 0.999, length.out = 201)
    vals <- sapply(grid, f)

    ok <- is.finite(vals[-length(vals)]) & is.finite(vals[-1L])
    idx <- which(ok & (sign(vals[-length(vals)]) != sign(vals[-1L])))

    if (!length(idx)) {
      ## Should now be rare because feasible-range check was already done
      if (f(1) < 0) {
        return(make_result(
          rho_hat = 1,
          boundary = TRUE,
          boundary_side = "upper",
          boundary_reason = "no_root_found_clamped_upper",
          feasible_r_min = feasible_r_min,
          feasible_r_max = feasible_r_max
        ))
      } else {
        return(make_result(
          rho_hat = -1,
          boundary = TRUE,
          boundary_side = "lower",
          boundary_reason = "no_root_found_clamped_lower",
          feasible_r_min = feasible_r_min,
          feasible_r_max = feasible_r_max
        ))
      }
    }

    if (idx[1] == length(grid)) {
      lower <- grid[length(grid) - 1]
      upper <- grid[length(grid)]
    } else {
      lower <- grid[idx[1]]
      upper <- grid[idx[1] + 1]
    }
  }

  ## Root finding
  root_val <- tryCatch(
    uniroot(f, lower = lower, upper = upper, tol = tol, maxiter = maxiter)$root,
    error = function(e) NA_real_
  )

  if (is.finite(root_val)) {
    return(make_result(
      rho_hat = root_val,
      boundary = FALSE,
      boundary_side = NA_character_,
      boundary_reason = NA_character_,
      feasible_r_min = feasible_r_min,
      feasible_r_max = feasible_r_max
    ))
  }

  ## Final fallback
  if (f(1) < 0) {
    return(make_result(
      rho_hat = 1,
      boundary = TRUE,
      boundary_side = "upper",
      boundary_reason = "uniroot_failed_clamped_upper",
      feasible_r_min = feasible_r_min,
      feasible_r_max = feasible_r_max
    ))
  } else {
    return(make_result(
      rho_hat = -1,
      boundary = TRUE,
      boundary_side = "lower",
      boundary_reason = "uniroot_failed_clamped_lower",
      feasible_r_min = feasible_r_min,
      feasible_r_max = feasible_r_max
    ))
  }
}


# Coerce model function #
.coerce_model <- function(model = c("congeneric","tau","parallel"),
                          lambdaA, lambdaB, thetaA, thetaB,
                          cA_vec = NULL, cB_vec = NULL,
                          overlapA = NULL, overlapB = NULL, k = 0L) {

  model <- match.arg(model)

  # infer nA/nB
  nA <- length(lambdaA)
  nB <- length(lambdaB)

  # recover lengths if lambdas are scalar but other inputs carry item-length
  if (nA == 1L) {
    if (!is.null(thetaA) && length(thetaA) > 1L) nA <- length(thetaA)
    if (!is.null(cA_vec) && length(cA_vec) > 1L) nA <- length(cA_vec)
    if (!is.null(overlapA) && length(overlapA) > 0L) nA <- max(nA, max(overlapA))
    if (!is.null(k) && k > 0L) nA <- max(nA, k)
    lambdaA <- rep(lambdaA, nA)
  }
  if (nB == 1L) {
    if (!is.null(thetaB) && length(thetaB) > 1L) nB <- length(thetaB)
    if (!is.null(cB_vec) && length(cB_vec) > 1L) nB <- length(cB_vec)
    if (!is.null(overlapB) && length(overlapB) > 0L) nB <- max(nB, max(overlapB))
    if (!is.null(k) && k > 0L) nB <- max(nB, k)
    lambdaB <- rep(lambdaB, nB)
  }

  # theta scalars -> vectors
  if (!is.null(thetaA) && length(thetaA) == 1L) thetaA <- rep(thetaA, nA)
  if (!is.null(thetaB) && length(thetaB) == 1L) thetaB <- rep(thetaB, nB)

  # enforce restrictions
  if (model %in% c("tau", "parallel")) {
    lambdaA <- rep(mean(lambdaA), nA)
    lambdaB <- rep(mean(lambdaB), nB)
  }
  if (model == "parallel") {
    if (!is.null(thetaA)) thetaA <- rep(mean(thetaA), nA)
    if (!is.null(thetaB)) thetaB <- rep(mean(thetaB), nB)
  }

  list(lambdaA = lambdaA, lambdaB = lambdaB, thetaA = thetaA, thetaB = thetaB)
}


#######################
### Master function ###
#######################
rho_correction_master <- function(
    r_obs,
    type = c("continuous", "ordinal"),
    spec = NULL,
    # continuous/ordinal shared parameters
    lambdaA, lambdaB,
    thetaA, thetaB,
    cA = 0, cB = 0,
    k = 0,
    cA_vec = NULL, cB_vec = NULL,
    overlapA = NULL, overlapB = NULL,
    # ordinal-only extras
    thrA = NULL, thrB = NULL,
    threshold_scale = c("z", "raw"),
    scoresA = NULL, scoresB = NULL,
    # CI / bootstrap
    r_SE = NULL,
    r_ci = NULL,
    ci = c("analytic", "bootstrap", "none"),
    n_boot = 0,
    conf_level = 0.95,
    return_samples = FALSE,
    # numeric inversion controls
    lower = -0.999, upper = 0.999,
    tol = 1e-6, maxiter = 100,
    # progress callback
    progress_callback = NULL,
    warn_boundary = TRUE
) {
  type <- match.arg(type)
  threshold_scale <- match.arg(threshold_scale)
  ci <- match.arg(ci)

  ## ---- 0. Basic input validation ----
  if (!is.numeric(r_obs) || length(r_obs) != 1L || !is.finite(r_obs)) {
    stop("r_obs must be a single finite numeric value.")
  }
  if (abs(r_obs) >= 1) {
    stop("r_obs must be strictly between -1 and 1.")
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be a single number strictly between 0 and 1.")
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1L || !is.finite(n_boot) ||
      n_boot < 0 || n_boot != as.integer(n_boot)) {
    stop("n_boot must be a single non-negative integer.")
  }
  if (!is.null(r_SE)) {
    if (!is.numeric(r_SE) || length(r_SE) != 1L || !is.finite(r_SE) || r_SE < 0) {
      stop("r_SE must be NULL or a single finite non-negative numeric value.")
    }
  }
  if (!is.logical(return_samples) || length(return_samples) != 1L || is.na(return_samples)) {
    stop("return_samples must be TRUE or FALSE.")
  }

  if (ci == "analytic" && is.null(r_SE) && is.null(r_ci)) {
    stop("For ci = 'analytic', provide either r_SE or r_ci.")
  }

  if (ci == "bootstrap" && is.null(r_SE)) {
    stop("For ci = 'bootstrap', provide r_SE.")
  }

  if (ci == "bootstrap" && n_boot <= 0) {
    stop("n_boot must be > 0 when ci = 'bootstrap'.")
  }

  if (return_samples && ci != "bootstrap") {
    warning("return_samples is only used when ci = 'bootstrap'.")
  }

   if (ci == "analytic" && is.null(r_SE) && is.null(r_ci)) {
    stop("For ci = 'analytic', provide either r_SE or r_ci.")
  }

  if (!is.null(r_ci)) {
    if (!is.numeric(r_ci) || length(r_ci) != 2L || any(!is.finite(r_ci))) {
      stop("r_ci must be NULL or a numeric vector of length 2: c(lower, upper).")
    }
    if (any(abs(r_ci) >= 1)) {
      stop("Both values in r_ci must be strictly between -1 and 1.")
    }
    r_ci <- sort(r_ci)
  }

  if (ci == "analytic" && is.null(r_SE) && is.null(r_ci)) {
    stop("For ci = 'analytic', provide either r_SE or r_ci.")
  }

  # --- infer measurement model from spec (if provided) ---
  model <- "congeneric"
  if (!is.null(spec)) {
    if (!is.null(spec$model)) {
      model <- tolower(spec$model)
    } else {
      cls <- tolower(class(spec))
      if (any(grepl("parallel", cls))) model <- "parallel"
      else if (any(grepl("tau", cls))) model <- "tau"
      else if (any(grepl("congeneric", cls))) model <- "congeneric"
    }
  }

  # --- coerce measurement model once ---
  tmp_par <- .coerce_model(
    model    = model,
    lambdaA  = lambdaA,  lambdaB  = lambdaB,
    thetaA   = thetaA,   thetaB   = thetaB,
    cA_vec   = cA_vec,   cB_vec   = cB_vec,
    overlapA = overlapA, overlapB = overlapB,
    k        = k
  )

  lambdaA <- tmp_par$lambdaA
  lambdaB <- tmp_par$lambdaB
  thetaA  <- tmp_par$thetaA
  thetaB  <- tmp_par$thetaB

  ## ---- 1. Single correction for a given r value ----
  correct_one <- function(r_val) {
    if (type == "continuous") {
      rho_correction_continuous(
        r_obs    = r_val,
        lambdaA  = lambdaA,
        lambdaB  = lambdaB,
        thetaA   = thetaA,
        thetaB   = thetaB,
        cA_vec   = cA_vec,
        cB_vec   = cB_vec,
        cA       = cA,
        cB       = cB,
        overlapA = overlapA,
        overlapB = overlapB,
        k        = k,
        lower    = lower,
        upper    = upper,
        tol      = tol,
        maxiter  = maxiter
      )
    } else {
      rho_correction_ordinal(
        r_obs    = r_val,
        lambdaA  = lambdaA,
        lambdaB  = lambdaB,
        thetaA   = thetaA,
        thetaB   = thetaB,
        cA_vec   = cA_vec,
        cB_vec   = cB_vec,
        cA       = cA,
        cB       = cB,
        overlapA = overlapA,
        overlapB = overlapB,
        k        = k,
        thrA     = thrA,
        thrB     = thrB,
        threshold_scale = threshold_scale,
        scoresA  = scoresA,
        scoresB  = scoresB,
        lower    = lower,
        upper    = upper,
        tol      = tol,
        maxiter  = maxiter
      )
    }
  }

  ## ---- 2. Point estimate ----
  fit0 <- correct_one(r_obs)

  rho_hat         <- fit0$rho_hat
  boundary        <- isTRUE(fit0$boundary)
  boundary_side   <- fit0$boundary_side
  boundary_reason <- fit0$boundary_reason
  feasible_r_min  <- fit0$feasible_r_min
  feasible_r_max  <- fit0$feasible_r_max

  if (isTRUE(boundary) && isTRUE(warn_boundary)) {
    warning(boundary_reason, call. = FALSE)
  }

  ## ---- 3. No CI requested ----
  if (ci == "none") {
    out <- list(
      rho_hat = rho_hat,
      boundary = boundary,
      boundary_side = boundary_side,
      boundary_reason = boundary_reason,
      feasible_r_min = feasible_r_min,
      feasible_r_max = feasible_r_max,
      method = "point_estimate_only",
      ci_method_used = "none",
      conf_level = conf_level,
      n_boot = n_boot,
      type = type,
      model = model,
      ci_computed = FALSE
    )
    class(out) <- "latentoverlap_result"
    return(out)
  }

  ## ---- 3b. Degenerate CI if r_SE == 0 ----
  if (!is.null(r_SE) && r_SE == 0 && is.null(r_ci)) {
    out <- list(
      rho_hat = rho_hat,
      CI_lower = rho_hat,
      CI_upper = rho_hat,
      SE_corrected = 0,
      n_boundary_boot = if (ci == "bootstrap") 0 else NA_integer_,
      prop_boundary_boot = if (ci == "bootstrap") 0 else NA_real_,
      ci_boundary = boundary,
      boundary = boundary,
      boundary_side = boundary_side,
      boundary_reason = boundary_reason,
      feasible_r_min = feasible_r_min,
      feasible_r_max = feasible_r_max,
      method = if (ci == "bootstrap") "parametric_bootstrap_fisher_z" else "analytic_fisher_z_endpoint",
      ci_method_used = ci,
      conf_level = conf_level,
      n_boot = n_boot,
      type = type,
      model = model,
      ci_computed = TRUE
    )

    if (return_samples && ci == "bootstrap") {
      out$rho_samples <- rep(rho_hat, n_boot)
      out$boundary_samples <- rep(boundary, n_boot)
    }

    class(out) <- "latentoverlap_result"
    return(out)
  }

  ## ---- 4. Common Fisher-z quantities ----
  alpha <- 1 - conf_level

  if (!is.null(r_SE)) {
    z_obs <- atanh(r_obs)
    z_se  <- r_SE / (1 - r_obs^2)
  }

  ## ---- 5. Analytic CI ----
  if (ci == "analytic") {

    if (!is.null(r_ci)) {
      r_low_raw <- r_ci[1]
      r_up_raw  <- r_ci[2]
    } else {
      z_low <- z_obs + stats::qnorm(alpha / 2) * z_se
      z_up  <- z_obs + stats::qnorm(1 - alpha / 2) * z_se

      r_low_raw <- tanh(z_low)
      r_up_raw  <- tanh(z_up)
    }

    r_low <- r_low_raw
    r_up  <- r_up_raw

    ci_boundary <- FALSE

    if (!is.null(feasible_r_min) && is.finite(feasible_r_min)) {
      if (r_low < feasible_r_min) {
        r_low <- feasible_r_min
        ci_boundary <- TRUE
      }
      if (r_up < feasible_r_min) {
        r_up <- feasible_r_min
        ci_boundary <- TRUE
      }
    }

    if (!is.null(feasible_r_max) && is.finite(feasible_r_max)) {
      if (r_low > feasible_r_max) {
        r_low <- feasible_r_max
        ci_boundary <- TRUE
      }
      if (r_up > feasible_r_max) {
        r_up <- feasible_r_max
        ci_boundary <- TRUE
      }
    }

    fit_low <- correct_one(r_low)
    fit_up  <- correct_one(r_up)

    rho_low <- fit_low$rho_hat
    rho_up  <- fit_up$rho_hat

    CI_lower <- min(rho_low, rho_up, na.rm = TRUE)
    CI_upper <- max(rho_low, rho_up, na.rm = TRUE)

    ci_boundary <- isTRUE(ci_boundary) ||
      isTRUE(fit_low$boundary) ||
      isTRUE(fit_up$boundary)

    if (isTRUE(ci_boundary) && isTRUE(warn_boundary)) {
      warning(
        "Analytic CI touches the feasible observed-r boundary; consider ci = 'bootstrap' as a sensitivity check.",
        call. = FALSE
      )
    }

    out <- list(
      rho_hat = rho_hat,
      CI_lower = CI_lower,
      CI_upper = CI_upper,
      SE_corrected = (CI_upper - CI_lower) / (2 * stats::qnorm(1 - alpha / 2)),
      ci_boundary = ci_boundary,
      boundary = boundary,
      boundary_side = boundary_side,
      boundary_reason = boundary_reason,
      feasible_r_min = feasible_r_min,
      feasible_r_max = feasible_r_max,
      r_CI_lower_raw = r_low_raw,
      r_CI_upper_raw = r_up_raw,
      r_CI_lower_used = r_low,
      r_CI_upper_used = r_up,
      method = if (!is.null(r_ci)) "analytic_observed_r_ci_endpoint" else "analytic_fisher_z_endpoint",
      ci_method_used = "analytic",
      conf_level = conf_level,
      n_boot = n_boot,
      type = type,
      model = model,
      ci_computed = TRUE
    )

    class(out) <- "latentoverlap_result"
    return(out)
  }

  ## ---- 6. Bootstrap CI (parametric, via Fisher z) ----
  z_samples <- stats::rnorm(n_boot, mean = z_obs, sd = z_se)
  r_samples <- tanh(z_samples)

  rho_samples      <- numeric(n_boot)
  boundary_samples <- logical(n_boot)

  for (i in seq_len(n_boot)) {
    if (!is.null(progress_callback)) progress_callback(i)

    tmp <- correct_one(r_samples[i])
    rho_samples[i]      <- tmp$rho_hat
    boundary_samples[i] <- isTRUE(tmp$boundary)
  }

  CI <- stats::quantile(
    rho_samples,
    probs = c(alpha / 2, 1 - alpha / 2),
    na.rm = TRUE,
    names = FALSE
  )

  out <- list(
    rho_hat = rho_hat,
    CI_lower = CI[1],
    CI_upper = CI[2],
    SE_corrected = stats::sd(rho_samples, na.rm = TRUE),
    n_boundary_boot = sum(boundary_samples, na.rm = TRUE),
    prop_boundary_boot = mean(boundary_samples, na.rm = TRUE),
    ci_boundary = any(boundary_samples, na.rm = TRUE),
    boundary = boundary,
    boundary_side = boundary_side,
    boundary_reason = boundary_reason,
    feasible_r_min = feasible_r_min,
    feasible_r_max = feasible_r_max,
    method = "parametric_bootstrap_fisher_z",
    ci_method_used = "bootstrap",
    conf_level = conf_level,
    n_boot = n_boot,
    type = type,
    model = model,
    ci_computed = TRUE
  )

  if (return_samples) {
    out$rho_samples <- rho_samples
    out$boundary_samples <- boundary_samples
  }

  class(out) <- "latentoverlap_result"
  out
}
