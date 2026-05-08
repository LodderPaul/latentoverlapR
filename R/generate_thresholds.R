#' Generate Thresholds for Ordinal Item Scores
#'
#' Generates threshold parameters for ordinal item responses under
#' symmetric, positively skewed, or negatively skewed response distributions.
#' The returned thresholds are defined on the standard normal latent-response
#' scale and can be used in \code{\link{rho_correction}} when
#' \code{type = "ordinal"}.
#'
#' Thresholds are currently implemented for ordinal items with
#' 2 to 7 response categories.
#'
#' @param n_categories Integer specifying the number of ordinal response
#' categories. Must be between 2 and 7.
#'
#' @param shape Character string specifying the response distribution shape.
#' One of:
#' \itemize{
#'   \item \code{"Symmetric"} for approximately symmetric category probabilities.
#'   \item \code{"Positive_skewness"} for positively skewed response distributions
#'   with higher prevalence of lower response categories.
#'   \item \code{"Negative_skewness"} for negatively skewed response distributions
#'   with higher prevalence of higher response categories.
#' }
#'
#' @return
#' A numeric vector containing threshold values on the latent-response scale.
#'
#' @details
#' The thresholds are predefined values chosen to produce plausible ordinal
#' response distributions for simulation studies and sensitivity analyses.
#' They can be supplied directly to the \code{thrA} and \code{thrB} arguments
#' of \code{\link{rho_correction}}.
#'
#' @examples
#' # Symmetric thresholds for a 5-category item
#' generate_thresholds(
#'   n_categories = 5,
#'   shape = "Symmetric"
#' )
#'
#' # Positively skewed thresholds for a 4-category item
#' generate_thresholds(
#'   n_categories = 4,
#'   shape = "Positive_skewness"
#' )
#'
#' @export
generate_thresholds <- function(n_categories,
                                shape = c("Symmetric", "Positive_skewness", "Negative_skewness"),
                                min_prob = 0.01,
                                skew_strength = 0.55,
                                return_probs = FALSE) {
  shape <- match.arg(shape)
  K <- as.integer(n_categories)

  if (K < 2 || K > 7) {
    stop("generate_thresholds: implemented only for 2-7 categories.")
  }

  if (min_prob <= 0 || min_prob * K >= 1) {
    stop("min_prob must be > 0 and min_prob * n_categories must be < 1.")
  }

  # Helper: convert category probabilities to normal thresholds
  probs_to_thresholds <- function(probs) {
    if (any(probs <= 0)) stop("All category probabilities must be > 0.")
    probs <- probs / sum(probs)
    qnorm(cumsum(probs)[-length(probs)])
  }

  # Symmetric target category probabilities
  sym_probs <- switch(
    as.character(K),
    "2" = c(.50, .50),
    "3" = c(.20, .60, .20),
    "4" = c(.10, .40, .40, .10),
    "5" = c(.05, .20, .50, .20, .05),
    "6" = c(.05, .15, .30, .30, .15, .05),
    "7" = c(.03, .09, .18, .40, .18, .09, .03)
  )

  if (shape == "Symmetric") {
    probs <- sym_probs
  } else {
    # Monotone skewed target probabilities.
    # Positive_skewness: many responses in low categories, few in high categories.
    raw <- exp(-skew_strength * (0:(K - 1)))
    probs_pos <- raw / sum(raw)

    # Enforce minimum category probability
    probs_pos <- pmax(probs_pos, min_prob)
    probs_pos <- probs_pos / sum(probs_pos)

    if (shape == "Positive_skewness") {
      probs <- probs_pos
    } else {
      probs <- rev(probs_pos)
    }
  }

  # Final safeguard
  probs <- pmax(probs, min_prob)
  probs <- probs / sum(probs)

  thresholds <- probs_to_thresholds(probs)

  if (return_probs) {
    return(list(
      thresholds = thresholds,
      probabilities = probs
    ))
  }

  thresholds
}



