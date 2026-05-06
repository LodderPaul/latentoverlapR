#' Specify a tau-equivalent measurement model
#'
#' Creates a measurement-model specification in which items within each scale
#' have equal factor loadings. Residual variances may differ across items if
#' supplied as vectors.
#'
#' @param lambdaA,lambdaB Factor loading(s) for scale A and scale B. Each may be
#'   a scalar or a vector of length `kA` / `kB`. If unequal vectors are supplied,
#'   their mean is used to enforce tau-equivalence.
#' @param kA,kB Number of items in scale A and scale B. Required when the
#'   corresponding loading and residual variance are supplied as scalars.
#' @param thetaA,thetaB Optional residual variance(s) for scale A and scale B.
#'   May be `NULL`, scalar, or a vector of length `kA` / `kB`.
#'
#' @return An object of class `"latentoverlap_spec"` for use in [rho_correction()].
#'
#' @examples
#' tau_model(
#'   lambdaA = .80,
#'   lambdaB = .75,
#'   kA = 6,
#'   kB = 5,
#'   thetaA = .36,
#'   thetaB = .44
#' )
#'
#' @export
tau_model <- function(lambdaA, lambdaB, kA = NULL, kB = NULL, thetaA = NULL, thetaB = NULL) {

  # infer number of items if not supplied
  if (is.null(kA)) {
    if (length(lambdaA) > 1) kA <- length(lambdaA)
    else if (!is.null(thetaA) && length(thetaA) > 1) kA <- length(thetaA)
    else stop("tau_model: argument 'kA' is missing. Provide kA when lambdaA/thetaA are scalars.")
  }
  if (is.null(kB)) {
    if (length(lambdaB) > 1) kB <- length(lambdaB)
    else if (!is.null(thetaB) && length(thetaB) > 1) kB <- length(thetaB)
    else stop("tau_model: argument 'kB' is missing. Provide kB when lambdaB/thetaB are scalars.")
  }

  # lambdaA
  if (length(lambdaA) == 1L) {
    la <- rep(lambdaA, kA)
  } else if (length(lambdaA) == kA) {
    if (length(unique(lambdaA)) > 1L) {
      warning("tau_model: lambdaA contains unequal values; using their mean to enforce tau-equivalence.")
      la <- rep(mean(lambdaA), kA)
    } else {
      la <- lambdaA
    }
  } else {
    stop("tau_model: lambdaA must have length 1 or kA.")
  }

  # lambdaB
  if (length(lambdaB) == 1L) {
    lb <- rep(lambdaB, kB)
  } else if (length(lambdaB) == kB) {
    if (length(unique(lambdaB)) > 1L) {
      warning("tau_model: lambdaB contains unequal values; using their mean to enforce tau-equivalence.")
      lb <- rep(mean(lambdaB), kB)
    } else {
      lb <- lambdaB
    }
  } else {
    stop("tau_model: lambdaB must have length 1 or kB.")
  }

  # thetaA: tau-equivalent does not require equal residuals, so keep vector if provided
  if (is.null(thetaA)) {
    ta <- NULL
  } else if (length(thetaA) == 1L) {
    ta <- rep(thetaA, kA)
  } else if (length(thetaA) == kA) {
    ta <- thetaA
  } else {
    stop("tau_model: thetaA must be NULL, length 1, or kA.")
  }

  # thetaB
  if (is.null(thetaB)) {
    tb <- NULL
  } else if (length(thetaB) == 1L) {
    tb <- rep(thetaB, kB)
  } else if (length(thetaB) == kB) {
    tb <- thetaB
  } else {
    stop("tau_model: thetaB must be NULL, length 1, or kB.")
  }

  structure(
    list(model = "tau", lambdaA = la, lambdaB = lb, thetaA = ta, thetaB = tb),
    class = "latentoverlap_spec"
  )
}

#' Specify a congeneric measurement model
#'
#' Creates a measurement-model specification in which item loadings may differ
#' across items within each scale. This is the most flexible measurement-model
#' specification currently supported by the correction function.
#'
#' @param lambdaA,lambdaB Numeric vectors of factor loadings for scale A and
#'   scale B.
#' @param thetaA,thetaB Optional numeric vectors of residual variances for scale A
#'   and scale B. If supplied, each must have the same length as the corresponding
#'   loading vector.
#'
#' @return An object of class `"latentoverlap_spec"` for use in [rho_correction()].
#'
#' @examples
#' congeneric_model(
#'   lambdaA = c(.70, .75, .80, .85, .90),
#'   lambdaB = c(.65, .70, .75, .80),
#'   thetaA = 1 - c(.70, .75, .80, .85, .90)^2,
#'   thetaB = 1 - c(.65, .70, .75, .80)^2
#' )
#'
#' @export
congeneric_model <- function(lambdaA, lambdaB, thetaA = NULL, thetaB = NULL) {
  stopifnot(length(lambdaA) >= 1, length(lambdaB) >= 1)

  if (!is.null(thetaA) && length(thetaA) != length(lambdaA)) {
    stop("thetaA must be NULL or same length as lambdaA.")
  }
  if (!is.null(thetaB) && length(thetaB) != length(lambdaB)) {
    stop("thetaB must be NULL or same length as lambdaB.")
  }

  structure(
    list(model = "congeneric", lambdaA = lambdaA, lambdaB = lambdaB, thetaA = thetaA, thetaB = thetaB),
    class = "latentoverlap_spec"
  )
}

#' Specify a parallel measurement model
#'
#' Creates a measurement-model specification in which items within each scale
#' have equal factor loadings and equal residual variances.
#'
#' @param lambdaA,lambdaB Factor loading(s) for scale A and scale B. Each may be
#'   a scalar or a vector of length `kA` / `kB`. If unequal vectors are supplied,
#'   their mean is used to enforce the parallel model.
#' @param kA,kB Number of items in scale A and scale B. Required when the
#'   corresponding loading and residual variance are supplied as scalars.
#' @param thetaA,thetaB Optional residual variance(s) for scale A and scale B.
#'   If `NULL`, residual variances are set to `1 - lambda^2`. If unequal vectors
#'   are supplied, their mean is used to enforce the parallel model.
#'
#' @return An object of class `"latentoverlap_spec"` for use in [rho_correction()].
#'
#' @examples
#' parallel_model(
#'   lambdaA = .80,
#'   lambdaB = .75,
#'   kA = 6,
#'   kB = 5
#' )
#'
#' @export
parallel_model <- function(lambdaA, lambdaB, kA = NULL, kB = NULL,
                          thetaA = NULL, thetaB = NULL) {

  # infer number of items if not supplied
  if (is.null(kA)) {
    if (length(lambdaA) > 1) kA <- length(lambdaA)
    else if (!is.null(thetaA) && length(thetaA) > 1) kA <- length(thetaA)
    else stop("parallel_model: argument 'kA' is missing. Provide kA when lambdaA/thetaA are scalars.")
  }
  if (is.null(kB)) {
    if (length(lambdaB) > 1) kB <- length(lambdaB)
    else if (!is.null(thetaB) && length(thetaB) > 1) kB <- length(thetaB)
    else stop("parallel_model: argument 'kB' is missing. Provide kB when lambdaB/thetaB are scalars.")
  }

  # lambdaA: must be equal within scale
  if (length(lambdaA) == 1L) {
    la <- rep(lambdaA, kA)
  } else if (length(lambdaA) == kA) {
    if (length(unique(lambdaA)) > 1L) {
      warning("parallel_model: lambdaA contains unequal values; using their mean to enforce parallel model.")
      la <- rep(mean(lambdaA), kA)
    } else {
      la <- lambdaA
    }
  } else {
    stop("parallel_model: lambdaA must have length 1 or kA.")
  }

  # lambdaB
  if (length(lambdaB) == 1L) {
    lb <- rep(lambdaB, kB)
  } else if (length(lambdaB) == kB) {
    if (length(unique(lambdaB)) > 1L) {
      warning("parallel_model: lambdaB contains unequal values; using their mean to enforce parallel model.")
      lb <- rep(mean(lambdaB), kB)
    } else {
      lb <- lambdaB
    }
  } else {
    stop("parallel_model: lambdaB must have length 1 or kB.")
  }

  # thetaA: must be equal within scale
  if (is.null(thetaA)) {
    ta <- rep(1 - mean(la)^2, kA)
  } else if (length(thetaA) == 1L) {
    ta <- rep(thetaA, kA)
  } else if (length(thetaA) == kA) {
    if (length(unique(thetaA)) > 1L) {
      warning("parallel_model: thetaA contains unequal values; using their mean to enforce parallel model.")
      ta <- rep(mean(thetaA), kA)
    } else {
      ta <- thetaA
    }
  } else {
    stop("parallel_model: thetaA must be NULL, length 1, or kA.")
  }

  # thetaB
  if (is.null(thetaB)) {
    tb <- rep(1 - mean(lb)^2, kB)
  } else if (length(thetaB) == 1L) {
    tb <- rep(thetaB, kB)
  } else if (length(thetaB) == kB) {
    if (length(unique(thetaB)) > 1L) {
      warning("parallel_model: thetaB contains unequal values; using their mean to enforce parallel model.")
      tb <- rep(mean(thetaB), kB)
    } else {
      tb <- thetaB
    }
  } else {
    stop("parallel_model: thetaB must be NULL, length 1, or kB.")
  }

  structure(
    list(model = "parallel", lambdaA = la, lambdaB = lb, thetaA = ta, thetaB = tb),
    class = "latentoverlap_spec"
  )
}

#' Specify a parallel model from scale reliabilities
#'
#' Converts total-score reliabilities and scale lengths into a parallel
#' measurement-model specification. This is useful when only scale reliabilities
#' are available rather than item-level factor loadings and residual variances.
#'
#' @param relA,relB Total-score reliability estimates for scale A and scale B.
#'   Values must be strictly between 0 and 1.
#' @param kA,kB Number of items in scale A and scale B. Each must be at least 2.
#'
#' @return An object of class `"latentoverlap_spec"` for use in [rho_correction()].
#'
#' @details
#' Under a parallel measurement model, the total-score reliability is used to
#' infer a common item loading and residual variance for each scale. This provides
#' a convenient approximation when published studies report only Cronbach's alpha,
#' omega, or another total-score reliability estimate.
#'
#' @examples
#' reliability_model(
#'   relA = .80,
#'   relB = .85,
#'   kA = 6,
#'   kB = 6
#' )
#'
#' @export
reliability_model <- function(relA, relB, kA, kB) {
  stopifnot(
    is.numeric(relA), length(relA) == 1L, relA > 0, relA < 1,
    is.numeric(relB), length(relB) == 1L, relB > 0, relB < 1,
    is.numeric(kA), length(kA) == 1L, kA >= 2,
    is.numeric(kB), length(kB) == 1L, kB >= 2
  )

  # Convert total-score reliability to item-level lambda^2
  lambdaA_sq <- relA / (kA - relA * (kA - 1))
  lambdaB_sq <- relB / (kB - relB * (kB - 1))

  # Numerical safeguard
  if (lambdaA_sq <= 0 || lambdaA_sq >= 1) {
    stop("Implied lambdaA^2 is outside (0,1); check relA and kA.")
  }
  if (lambdaB_sq <= 0 || lambdaB_sq >= 1) {
    stop("Implied lambdaB^2 is outside (0,1); check relB and kB.")
  }

  lambdaA <- sqrt(lambdaA_sq)
  lambdaB <- sqrt(lambdaB_sq)

  thetaA <- 1 - lambdaA_sq
  thetaB <- 1 - lambdaB_sq

  parallel_model(
    lambdaA = lambdaA,
    lambdaB = lambdaB,
    thetaA = thetaA,
    thetaB = thetaB,
    kA = kA,
    kB = kB
  )
}
