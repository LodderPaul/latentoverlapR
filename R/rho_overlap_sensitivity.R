# Internal helper -----------------------------------------------------------

.make_overlap_inputs <- function(
    vary = c("multiplier", "common", "item"),
    value,
    cA = 0, cB = 0,
    cA_vec = NULL, cB_vec = NULL,
    overlapA = NULL, overlapB = NULL,
    focal_scale = c("A", "B"),
    focal_item = NULL
) {
  vary <- match.arg(vary)
  focal_scale <- match.arg(focal_scale)

  out <- list(
    cA = cA,
    cB = cB,
    cA_vec = cA_vec,
    cB_vec = cB_vec,
    overlapA = overlapA,
    overlapB = overlapB
  )

  if (vary == "multiplier") {
    out$cA <- value * cA
    out$cB <- value * cB
    if (!is.null(cA_vec)) out$cA_vec <- value * cA_vec
    if (!is.null(cB_vec)) out$cB_vec <- value * cB_vec
    return(out)
  }

  if (vary == "common") {
    if (!is.null(cA_vec)) {
      out$cA_vec <- cA_vec
      idxA <- which(abs(cA_vec) > 0)
      if (length(idxA)) out$cA_vec[idxA] <- value
    } else {
      out$cA <- value
    }

    if (!is.null(cB_vec)) {
      out$cB_vec <- cB_vec
      idxB <- which(abs(cB_vec) > 0)
      if (length(idxB)) out$cB_vec[idxB] <- value
    } else {
      out$cB <- value
    }

    return(out)
  }

  if (vary == "item") {
    if (is.null(focal_item) || length(focal_item) != 1L || !is.finite(focal_item)) {
      stop("For vary = 'item', focal_item must be a single valid item index.")
    }

    focal_item <- as.integer(focal_item)

    if (focal_scale == "A") {
      if (is.null(cA_vec)) stop("For focal_scale = 'A', cA_vec must be provided.")
      if (focal_item < 1L || focal_item > length(cA_vec)) {
        stop("focal_item is out of bounds for cA_vec.")
      }
      out$cA_vec <- cA_vec
      out$cA_vec[focal_item] <- value
    }

    if (focal_scale == "B") {
      if (is.null(cB_vec)) stop("For focal_scale = 'B', cB_vec must be provided.")
      if (focal_item < 1L || focal_item > length(cB_vec)) {
        stop("focal_item is out of bounds for cB_vec.")
      }
      out$cB_vec <- cB_vec
      out$cB_vec[focal_item] <- value
    }

    return(out)
  }

  out
}


# Public function ----------------------------------------------------------
#' Sensitivity of corrected association to assumed measurement overlap
#'
#' Varies the assumed measurement-overlap intensity and recomputes the corrected
#' latent association at each value of a user-specified grid. This is useful when
#' the degree of overlap is uncertain and researchers want to evaluate how robust
#' their corrected association is to alternative overlap assumptions.
#'
#' @param r_obs Observed total-score correlation. Must be strictly between -1 and 1.
#' @param vary Character string indicating how overlap is varied. `"multiplier"`
#'   multiplies the supplied overlap values by each grid value. `"common"`
#'   replaces all non-zero overlap values by the grid value. `"item"` varies one
#'   focal item-specific cross-loading.
#' @param grid Numeric vector of overlap values or multipliers to evaluate.
#' @param focal_scale Character string. Either `"A"` or `"B"`. Only used when
#'   `vary = "item"` and indicates whether an item-specific cross-loading in
#'   `cA_vec` or `cB_vec` is varied.
#' @param focal_item Integer. Item index to vary when `vary = "item"`.
#' @param type Character string. Either `"continuous"` for continuous item scores
#'   or `"ordinal"` for ordinal item scores.
#' @param measurement_model A measurement model specification created with
#'   `reliability_model()`, `tau_model()`, `parallel_model()`, or
#'   `congeneric_model()`.
#' @param cA,cB Scalar cross-loading values for overlapping items from scale A on
#'   construct B and from scale B on construct A. Used when item-specific vectors
#'   are not supplied.
#' @param k_overlap Number of overlapping items per scale when using scalar
#'   cross-loading values.
#' @param overlapA,overlapB Optional integer vectors identifying overlapping items
#'   in scale A and scale B.
#' @param cA_vec,cB_vec Optional item-specific cross-loading vectors for scale A
#'   and scale B.
#' @param thrA,thrB Threshold vectors for ordinal items. Required when
#'   `type = "ordinal"`.
#' @param threshold_scale Scale of the thresholds. Either `"z"` for thresholds on
#'   the standard-normal latent-response scale or `"raw"` for raw thresholds.
#' @param scoresA,scoresB Optional score values assigned to ordinal response
#'   categories.
#' @param ci Confidence interval method. `"analytic"` applies the correction to
#'   the confidence limits of the observed correlation. `"bootstrap"` uses a
#'   parametric Fisher-z bootstrap. `"none"` returns only point estimates.
#' @param r_SE Standard error of the observed total-score correlation. Required
#'   for `ci = "bootstrap"` and used for `ci = "analytic"` when `r_ci` is not
#'   supplied.
#' @param r_ci Optional numeric vector of length 2 giving the lower and upper
#'   confidence limits of the observed total-score correlation. Can be used with
#'   `ci = "analytic"` instead of `r_SE`.
#' @param n_boot Number of bootstrap draws. Only used when `ci = "bootstrap"`.
#' @param conf_level Confidence level for analytic or bootstrap confidence
#'   intervals.
#' @param return_fits Logical. If `TRUE`, returns the full correction object for
#'   each grid value.
#' @param warn_boundary Logical. If `TRUE`, boundary-related information is
#'   retained. Warnings from individual grid evaluations are suppressed to avoid
#'   repeated messages.
#' @param show_progress Logical. If `TRUE`, shows a text progress bar.
#' @param progress_callback Optional function called at each grid value as
#'   `progress_callback(i, n)`.
#' @param ... Additional arguments passed to [rho_correction()].
#'
#' @return An object of class `"rho_overlap_sensitivity"` with the following
#'   elements:
#' \describe{
#'   \item{critical_overlap}{The first grid value at which the confidence
#'   interval includes zero, if available.}
#'   \item{observed_overlap}{The grid value closest to 1, representing the
#'   originally assumed overlap under `vary = "multiplier"`.}
#'   \item{baseline_overlap_value}{The non-zero baseline overlap value when a
#'   single common value can be inferred.}
#'   \item{ci_includes_zero_at_observed}{Logical indicating whether the
#'   confidence interval includes zero at the observed overlap value.}
#'   \item{results}{A data frame with one row per grid value, including the
#'   corrected estimate, confidence interval limits, boundary information, and
#'   zero-inclusion indicators.}
#'   \item{fits}{If `return_fits = TRUE`, a list of full correction results.}
#' }
#'
#' @details
#' The function repeatedly calls [rho_correction()] while changing the assumed
#' overlap intensity. With `vary = "multiplier"`, the grid is interpreted as a
#' multiplier applied to the originally supplied overlap values. For example,
#' when `cA = cB = .30`, `grid = 1` corresponds to cross-loadings of `.30`,
#' whereas `grid = 2` corresponds to cross-loadings of `.60`.
#'
#' With `vary = "multiplier"`, the grid values multiply the originally supplied
#' overlap values. For example, if `cA = cB = .30`, then `grid = 1` evaluates the
#' original overlap and `grid = 2` evaluates cross-loadings of `.60`.
#'
#' With `vary = "common"`, all non-zero overlap values are replaced by each grid
#' value. This is useful when researchers want to ask how large the overlap would
#' have to be, in absolute cross-loading units, for conclusions to change.
#'
#' With `vary = "item"`, only one item-specific cross-loading is varied. This is
#' useful when uncertainty concerns a particular item rather than the full overlap
#' structure.
#'
#' When confidence intervals are requested, the output identifies the first grid
#' value at which the corrected confidence interval includes zero. This can help
#' evaluate how strong the assumed overlap would need to be for the corrected
#' association to no longer be distinguishable from zero.
#'
#' @section Methods:
#' Objects of class `"rho_overlap_sensitivity"` have the following methods:
#' \describe{
#'   \item{\code{print(x)}}{Displays a summary including the observed overlap,
#'   critical overlap (if available), and whether the confidence interval
#'   includes zero.}
#'   \item{\code{plot(x)}}{Plots the corrected association as a function of the
#'   overlap value, optionally including confidence interval bands and markers
#'   for observed and critical overlap.}
#' }
#'
#' @examples
#' mod <- reliability_model(
#'   relA = .80,
#'   relB = .85,
#'   kA = 6,
#'   kB = 6
#' )
#'
#' ## Sensitivity to multiplying the assumed overlap
#' sens <- rho_overlap_sensitivity(
#'   r_obs = .60,
#'   type = "continuous",
#'   measurement_model = mod,
#'   cA = .30,
#'   cB = .30,
#'   k_overlap = 1,
#'   ci = "analytic",
#'   r_SE = .10,
#'   grid = seq(0, 2, by = .25),
#'   show_progress = FALSE
#' )
#'
#' head(sens$results)
#'
#' ## Plot the sensitivity curve
#' plot(sens)
#'
#' ## Use reported confidence limits for the observed correlation
#' sens_rci <- rho_overlap_sensitivity(
#'   r_obs = .60,
#'   type = "continuous",
#'   measurement_model = mod,
#'   cA = .30,
#'   cB = .30,
#'   k_overlap = 1,
#'   ci = "analytic",
#'   r_ci = c(.45, .72),
#'   grid = seq(0, 2, by = .25),
#'   show_progress = FALSE
#' )
#'
#' head(sens_rci$results)
#'
#'#' ## Sensitivity by replacing all non-zero overlap values
#' sens_common <- rho_overlap_sensitivity(
#'   r_obs = .60,
#'   type = "continuous",
#'   measurement_model = mod,
#'   cA_vec = c(0, 0, 0, 0, .30, .20),
#'   cB_vec = c(.25, .15, 0, 0, 0, 0),
#'   ci = "analytic",
#'   r_SE = .10,
#'   vary = "common",
#'   grid = seq(0, .60, by = .10),
#'   show_progress = FALSE
#' )
#'
#' head(sens_common$results)
#'
#' ## Sensitivity for one item-specific overlap parameter
#' sens_item <- rho_overlap_sensitivity(
#'   r_obs = .60,
#'   type = "continuous",
#'   measurement_model = mod,
#'   cA_vec = c(0, 0, 0, 0, .30, .20),
#'   cB_vec = c(.25, .15, 0, 0, 0, 0),
#'   ci = "analytic",
#'   r_SE = .10,
#'   vary = "item",
#'   focal_scale = "A",
#'   focal_item = 5,
#'   grid = seq(0, .60, by = .10),
#'   show_progress = FALSE
#' )
#'
#' head(sens_item$results)
#'
#' plot(sens_item)
#'
#'
#' @export
rho_overlap_sensitivity <- function(
    r_obs,
    vary = c("multiplier", "common", "item"),
    focal_scale = c("A", "B"),
    focal_item = NULL,
    grid = seq(0, 2, by = 0.01),
    type = c("continuous", "ordinal"),
    measurement_model,
    cA = 0, cB = 0,
    k_overlap = 0,
    overlapA = NULL, overlapB = NULL,
    cA_vec = NULL, cB_vec = NULL,
    thrA = NULL, thrB = NULL,
    threshold_scale = c("z", "raw"),
    scoresA = NULL, scoresB = NULL,
    ci = c("analytic", "bootstrap", "none"),
    r_SE = NULL,
    r_ci = NULL,
    n_boot = 0, conf_level = 0.95,
    return_fits = FALSE,
    warn_boundary = TRUE,
    show_progress = interactive(),
    progress_callback = NULL,
    ...
) {
  type <- match.arg(type)
  threshold_scale <- match.arg(threshold_scale)
  ci <- match.arg(ci)
  vary <- match.arg(vary)
  focal_scale <- match.arg(focal_scale)

  if (!inherits(measurement_model, "latentoverlap_spec")) {
    stop("measurement_model must be a latentoverlap_spec.")
  }

  if (!is.numeric(r_obs) || length(r_obs) != 1L || !is.finite(r_obs) || abs(r_obs) >= 1) {
    stop("r_obs must be a single finite numeric value strictly between -1 and 1.")
  }

  if (!is.numeric(grid) || length(grid) < 2L || any(!is.finite(grid))) {
    stop("grid must be a numeric vector with at least two finite values.")
  }

  if (!is.logical(show_progress) || length(show_progress) != 1L || is.na(show_progress)) {
    stop("show_progress must be TRUE or FALSE.")
  }

  grid <- sort(unique(grid))

  no_overlap_info <- is.null(cA_vec) &&
    is.null(cB_vec) &&
    identical(cA, 0) &&
    identical(cB, 0)

  if (no_overlap_info) {
    stop("No overlap parameters supplied to vary.")
  }

  if (type == "ordinal" && (is.null(thrA) || is.null(thrB))) {
    stop("For type = 'ordinal', provide thrA and thrB.")
  }

  if (ci == "bootstrap") {
    if (is.null(r_SE) || !is.finite(r_SE) || r_SE < 0) {
      stop("For ci = 'bootstrap', provide a non-negative r_SE.")
    }
    if (!is.numeric(n_boot) || length(n_boot) != 1L || !is.finite(n_boot) || n_boot <= 0) {
      stop("For ci = 'bootstrap', n_boot must be > 0.")
    }
  }

  if (ci == "analytic" && is.null(r_SE) && is.null(r_ci)) {
    stop("For ci = 'analytic', provide either r_SE or r_ci.")
  }

  if (ci == "bootstrap") {
    if (is.null(r_SE) || !is.numeric(r_SE) || length(r_SE) != 1L ||
        !is.finite(r_SE) || r_SE < 0) {
      stop("For ci = 'bootstrap', provide a non-negative r_SE.")
    }
    if (!is.numeric(n_boot) || length(n_boot) != 1L ||
        !is.finite(n_boot) || n_boot <= 0) {
      stop("For ci = 'bootstrap', n_boot must be > 0.")
    }
  }

  n <- length(grid)

  rho_hat <- rep(NA_real_, n)
  CI_lower <- rep(NA_real_, n)
  CI_upper <- rep(NA_real_, n)
  SE_corrected <- rep(NA_real_, n)

  boundary <- rep(FALSE, n)
  boundary_side <- rep(NA_character_, n)
  boundary_reason <- rep(NA_character_, n)
  feasible_r_min <- rep(NA_real_, n)
  feasible_r_max <- rep(NA_real_, n)

  n_boundary_boot <- rep(NA_integer_, n)
  prop_boundary_boot <- rep(NA_real_, n)

  fits <- if (return_fits) vector("list", n) else NULL

  pb <- NULL
  if (isTRUE(show_progress)) {
    pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (i in seq_len(n)) {
    s <- grid[i]

    if (!is.null(progress_callback)) {
      progress_callback(i, n)
    }
    if (!is.null(pb)) utils::setTxtProgressBar(pb, i)

    ov <- .make_overlap_inputs(
      vary = vary,
      value = s,
      cA = cA,
      cB = cB,
      cA_vec = cA_vec,
      cB_vec = cB_vec,
      overlapA = overlapA,
      overlapB = overlapB,
      focal_scale = focal_scale,
      focal_item = focal_item
    )

    fit <- rho_correction(
      r_obs = r_obs,
      type = type,
      measurement_model = measurement_model,
      cA = ov$cA,
      cB = ov$cB,
      k_overlap = k_overlap,
      overlapA = ov$overlapA,
      overlapB = ov$overlapB,
      cA_vec = ov$cA_vec,
      cB_vec = ov$cB_vec,
      thrA = thrA,
      thrB = thrB,
      threshold_scale = threshold_scale,
      scoresA = scoresA,
      scoresB = scoresB,
      ci = ci,
      r_SE = r_SE,
      r_ci = r_ci,
      n_boot = if (ci == "bootstrap") n_boot else 0,
      conf_level = conf_level,
      warn_boundary = FALSE,
      ...
    )

    rho_hat[i] <- fit$rho_hat
    if (!is.null(fit$CI_lower)) CI_lower[i] <- fit$CI_lower
    if (!is.null(fit$CI_upper)) CI_upper[i] <- fit$CI_upper
    if (!is.null(fit$SE_corrected)) SE_corrected[i] <- fit$SE_corrected

    boundary[i] <- isTRUE(fit$boundary)
    if (!is.null(fit$boundary_side)) boundary_side[i] <- fit$boundary_side
    if (!is.null(fit$boundary_reason)) boundary_reason[i] <- fit$boundary_reason
    if (!is.null(fit$feasible_r_min)) feasible_r_min[i] <- fit$feasible_r_min
    if (!is.null(fit$feasible_r_max)) feasible_r_max[i] <- fit$feasible_r_max

    if (!is.null(fit$n_boundary_boot)) n_boundary_boot[i] <- fit$n_boundary_boot
    if (!is.null(fit$prop_boundary_boot)) prop_boundary_boot[i] <- fit$prop_boundary_boot

    if (return_fits) fits[[i]] <- fit
  }

  results <- data.frame(
    overlap_value = grid,
    rho_hat = rho_hat,
    CI_lower = CI_lower,
    CI_upper = CI_upper,
    SE_corrected = SE_corrected,
    boundary = boundary,
    boundary_side = boundary_side,
    boundary_reason = boundary_reason,
    feasible_r_min = feasible_r_min,
    feasible_r_max = feasible_r_max,
    n_boundary_boot = n_boundary_boot,
    prop_boundary_boot = prop_boundary_boot,
    stringsAsFactors = FALSE
  )

  if (ci != "none") {
    results$ci_includes_zero <-
      !is.na(results$CI_lower) &
      !is.na(results$CI_upper) &
      results$CI_lower <= 0 &
      results$CI_upper >= 0
  } else {
    results$ci_includes_zero <- NA
  }

  results$rho_crosses_zero <- FALSE
  sign_ok <- !is.na(results$rho_hat)
  results$rho_crosses_zero[sign_ok] <- abs(results$rho_hat[sign_ok]) < .Machine$double.eps^0.5

  if (vary == "multiplier") {
    idx_obs <- which.min(abs(results$overlap_value - 1))
    observed_overlap <- results$overlap_value[idx_obs]
  } else {

    baseline_overlap_value <- NA_real_

    nonzero_overlap_values <- c()

    if (!is.null(cA_vec)) {
      nonzero_overlap_values <- c(
        nonzero_overlap_values,
        cA_vec[abs(cA_vec) > 0]
      )
    }

    if (!is.null(cB_vec)) {
      nonzero_overlap_values <- c(
        nonzero_overlap_values,
        cB_vec[abs(cB_vec) > 0]
      )
    }

    if (length(nonzero_overlap_values) == 0L) {
      if (!identical(cA, 0)) nonzero_overlap_values <- c(nonzero_overlap_values, cA)
      if (!identical(cB, 0)) nonzero_overlap_values <- c(nonzero_overlap_values, cB)
    }

    if (length(nonzero_overlap_values) > 0L) {
      observed_overlap <- nonzero_overlap_values[1]
    } else {
      observed_overlap <- NA_real_
    }

    idx_obs <- which.min(abs(results$overlap_value - observed_overlap))
  }

  if (ci != "none") {
    idx_crit <- which(results$ci_includes_zero)[1]
    critical_overlap <- if (length(idx_crit)) results$overlap_value[idx_crit] else NA_real_
    ci_includes_zero_at_observed <- results$ci_includes_zero[idx_obs]
  } else {
    critical_overlap <- NA_real_
    ci_includes_zero_at_observed <- NA
  }

  nonzero_overlap_values <- c()

  if (!is.null(cA_vec)) nonzero_overlap_values <- c(nonzero_overlap_values, cA_vec[abs(cA_vec) > 0])
  if (!is.null(cB_vec)) nonzero_overlap_values <- c(nonzero_overlap_values, cB_vec[abs(cB_vec) > 0])

  if (length(nonzero_overlap_values) == 0L) {
    if (!identical(cA, 0)) nonzero_overlap_values <- c(nonzero_overlap_values, cA)
    if (!identical(cB, 0)) nonzero_overlap_values <- c(nonzero_overlap_values, cB)
  }

  baseline_overlap_value <- NA_real_
  if (length(nonzero_overlap_values) > 0L) {
    unique_vals <- unique(round(nonzero_overlap_values, 10))
    if (length(unique_vals) == 1L) {
      baseline_overlap_value <- unique_vals[1]
    }
  }

  out <- list(
    critical_overlap = critical_overlap,
    observed_overlap = observed_overlap,
    baseline_overlap_value = baseline_overlap_value,
    ci_includes_zero_at_observed = ci_includes_zero_at_observed,
    vary = vary,
    focal_scale = if (vary == "item") focal_scale else NULL,
    focal_item = if (vary == "item") focal_item else NULL,
    ci = ci,
    conf_level = conf_level,
    type = type,
    grid = grid,
    results = results
  )

  if (return_fits) out$fits <- fits

  class(out) <- "rho_overlap_sensitivity"
  out
}
