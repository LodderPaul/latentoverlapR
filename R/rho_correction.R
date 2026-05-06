#' Correct an observed correlation for measurement overlap
#'
#' Re-expresses an observed total-score correlation as the latent construct-level
#' correlation implied by a specified measurement model. The correction accounts
#' for attenuation due to measurement error and for inflation due to overlapping
#' measurement indicators.
#'
#' @param r_obs Observed correlation between two total scores. Must be strictly
#'   between -1 and 1.
#' @param type Character string. Either `"continuous"` for continuous item scores
#'   or `"ordinal"` for ordinal item scores.
#' @param measurement_model A measurement model specification created with
#'   [reliability_model()], [parallel_model()], [tau_model()], or [congeneric_model()].
#' @param cA,cB Cross-loading size for overlapping items from scale A on construct
#'   B and from scale B on construct A. Used when item-specific vectors are not
#'   provided.
#' @param k_overlap Number of overlapping items per scale when using scalar
#'   cross-loading values.
#' @param overlapA,overlapB Optional integer vectors indicating which items in
#'   scale A and scale B overlap.
#' @param cA_vec,cB_vec Optional item-specific cross-loading vectors for scale A
#'   and scale B.
#' @param thrA,thrB Thresholds for ordinal items. Required when
#'   `type = "ordinal"`.
#' @param threshold_scale Scale of the thresholds. Either `"z"` for thresholds on
#'   the standard-normal latent-response scale or `"raw"` for raw thresholds.
#' @param scoresA,scoresB Optional score values for ordinal response categories.
#' @param ci Confidence interval method. `"analytic"` applies the correction to
#'   the confidence limits of the observed correlation and is the recommended
#'   default. `"bootstrap"` uses a parametric Fisher-z bootstrap. `"none"` returns
#'   only the point estimate.
#' @param r_SE Standard error of the observed correlation. Required for
#'   `ci = "bootstrap"` and used for `ci = "analytic"` when `r_ci` is not supplied.
#' @param r_ci Optional numeric vector of length 2 giving the lower and upper
#'   confidence limits of the observed correlation. Can be used with
#'   `ci = "analytic"` instead of `r_SE`.
#' @param n_boot Number of parametric bootstrap draws. Only used when
#'   `ci = "bootstrap"`.
#' @param conf_level Confidence level for the interval. Defaults to .95.
#' @param warn_boundary Logical. If `TRUE`, warns when the inversion reaches the
#'   feasible boundary implied by the measurement model.
#' @param ... Additional arguments passed to the internal correction routines.
#'
#' @return An object of class `"latentoverlap_result"` containing at least:
#' \describe{
#'   \item{rho_hat}{Corrected latent correlation.}
#'   \item{CI_lower, CI_upper}{Confidence interval limits, if requested.}
#'   \item{SE_corrected}{Approximate standard error of the corrected estimate, if available.}
#'   \item{boundary}{Logical indicating whether the solution hit the feasible boundary.}
#'   \item{feasible_r_min, feasible_r_max}{Feasible observed correlation range implied by the model.}
#'   \item{method}{Method used for the correction or confidence interval.}
#' }
#'
#' @details
#' The correction is conditional on the supplied measurement model. Therefore,
#' the corrected correlation should be interpreted as the latent correlation
#' implied by the assumed loadings, reliabilities, thresholds, and overlap
#' structure.
#'
#' For confidence intervals, the recommended default is `ci = "analytic"`.
#' If `r_SE` is supplied, a confidence interval for the observed correlation is
#' first constructed using Fisher's z-transformation and then the correction is
#' applied to its endpoints. Alternatively, users may supply the confidence
#' limits of the observed correlation directly via `r_ci`.
#'
#' The bootstrap option uses the same Fisher-z uncertainty model but repeatedly
#' draws observed correlations and corrects each draw. This can be slow,
#' especially for ordinal models.
#'
#' @examples
#' ## Example 1: Continuous items, reliability model, point estimate only
#' rho_correction(
#'   r_obs = .60,
#'   type = "continuous",
#'   measurement_model = reliability_model(
#'     relA = .80,
#'     relB = .85,
#'     kA = 6,
#'     kB = 6
#'   ),
#'   cA = .30,
#'   cB = .30,
#'   k_overlap = 1,
#'   ci = "none"
#' )
#'
#' ## Example 2: Continuous items with analytic CI from r_SE
#' rho_correction(
#'   r_obs = .60,
#'   type = "continuous",
#'   measurement_model = reliability_model(
#'     relA = .80,
#'     relB = .85,
#'     kA = 6,
#'     kB = 6
#'   ),
#'   cA = .30,
#'   cB = .30,
#'   k_overlap = 1,
#'   ci = "analytic",
#'   r_SE = .10
#' )
#'
#' ## Example 3: Continuous items with analytic CI from reported r confidence limits
#' rho_correction(
#'   r_obs = .60,
#'   type = "continuous",
#'   measurement_model = reliability_model(
#'     relA = .80,
#'     relB = .85,
#'     kA = 6,
#'     kB = 6
#'   ),
#'   cA = .30,
#'   cB = .30,
#'   k_overlap = 1,
#'   ci = "analytic",
#'   r_ci = c(.45, .72)
#' )
#'
#' ## Example 4: Ordinal items with symmetric five-category thresholds
#' rho_correction(
#'   r_obs = .60,
#'   type = "ordinal",
#'   measurement_model = reliability_model(
#'     relA = .80,
#'     relB = .85,
#'     kA = 6,
#'     kB = 6
#'   ),
#'   cA = .30,
#'   cB = .30,
#'   k_overlap = 1,
#'   thrA = c(-1.645, -0.643, 0.643, 1.645),
#'   thrB = c(-1.645, -0.643, 0.643, 1.645),
#'   ci = "analytic",
#'   r_SE = .10
#' )
#'
#' ## Example 5: Ordinal items with a reported CI for the observed correlation
#' rho_correction(
#'   r_obs = .60,
#'   type = "ordinal",
#'   measurement_model = reliability_model(
#'     relA = .80,
#'     relB = .85,
#'     kA = 6,
#'     kB = 6
#'   ),
#'   cA = .30,
#'   cB = .30,
#'   k_overlap = 1,
#'   thrA = c(-1.645, -0.643, 0.643, 1.645),
#'   thrB = c(-1.645, -0.643, 0.643, 1.645),
#'   ci = "analytic",
#'   r_ci = c(.45, .72)
#' )
#'
#' ## Example 6: Item-specific overlap using cross-loading vectors
#' rho_correction(
#'   r_obs = .60,
#'   type = "continuous",
#'   measurement_model = tau_model(
#'     lambdaA = rep(.80, 6),
#'     lambdaB = rep(.80, 6),
#'     thetaA = rep(.36, 6),
#'     thetaB = rep(.36, 6)
#'   ),
#'   cA_vec = c(0, 0, 0, 0, .30, .20),
#'   cB_vec = c(.25, .15, 0, 0, 0, 0),
#'   ci = "analytic",
#'   r_SE = .10
#' )
#'
#' ## Example 7: Parametric bootstrap CI
#' ## Kept small here for speed; use a larger n_boot in applied work.
#' rho_correction(
#'   r_obs = .60,
#'   type = "continuous",
#'   measurement_model = reliability_model(
#'     relA = .80,
#'     relB = .85,
#'     kA = 6,
#'     kB = 6
#'   ),
#'   cA = .30,
#'   cB = .30,
#'   k_overlap = 1,
#'   ci = "bootstrap",
#'   r_SE = .10,
#'   n_boot = 50
#' )
#'
#' @export



#' @export
rho_correction <- function(
    r_obs,
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
    n_boot = 0,
    conf_level = 0.95,
    warn_boundary = TRUE,
    ...
) {
  type <- match.arg(type)
  threshold_scale <- match.arg(threshold_scale)
  ci <- match.arg(ci)

  if (!inherits(measurement_model, "latentoverlap_spec")) {
    stop("measurement_model must be a latentoverlap_spec (use tau_model(), reliability_model(), or spec_congeneric()).")
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

  if (ci == "bootstrap") {
    if (is.null(r_SE) || !is.numeric(r_SE) || length(r_SE) != 1L ||
        !is.finite(r_SE) || r_SE < 0) {
      stop("For ci = 'bootstrap', provide a non-negative r_SE.")
    }
    if (n_boot <= 0) {
      stop("For ci = 'bootstrap', n_boot must be > 0.")
    }
  }

  if (!is.null(r_SE)) {
    if (!is.numeric(r_SE) || length(r_SE) != 1L ||
        !is.finite(r_SE) || r_SE < 0) {
      stop("r_SE must be NULL or a single finite non-negative numeric value.")
    }
  }

  if (type == "ordinal" && (is.null(thrA) || is.null(thrB))) {
    stop("For type = 'ordinal', provide thrA and thrB.")
  }

  model <- "congeneric"
  if (!is.null(measurement_model$model)) {
    model <- tolower(measurement_model$model)
  } else {
    cls <- tolower(class(measurement_model))
    if (any(grepl("parallel", cls))) model <- "parallel"
    else if (any(grepl("tau", cls))) model <- "tau"
    else if (any(grepl("congeneric", cls))) model <- "congeneric"
  }

  lambdaA <- measurement_model$lambdaA
  lambdaB <- measurement_model$lambdaB
  thetaA  <- measurement_model$thetaA
  thetaB  <- measurement_model$thetaB

  tmp <- .coerce_model(
    model    = model,
    lambdaA  = lambdaA,
    lambdaB  = lambdaB,
    thetaA   = thetaA,
    thetaB   = thetaB,
    cA_vec   = cA_vec,
    cB_vec   = cB_vec,
    overlapA = overlapA,
    overlapB = overlapB,
    k        = k_overlap
  )

  rho_correction_master(
    r_obs = r_obs,
    type = type,
    spec = measurement_model,
    lambdaA = tmp$lambdaA,
    lambdaB = tmp$lambdaB,
    thetaA  = tmp$thetaA,
    thetaB  = tmp$thetaB,
    cA = cA,
    cB = cB,
    k = k_overlap,
    cA_vec = cA_vec,
    cB_vec = cB_vec,
    overlapA = overlapA,
    overlapB = overlapB,
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
    warn_boundary = warn_boundary,
    ...
  )
}
