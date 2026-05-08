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
#' @param min_prob Minimum allowable probability for any response category.
#' Used to avoid extremely sparse categories. Must be greater than 0 and small
#' enough that \code{min_prob * n_categories < 1}.
#'
#' @param skew_strength Numeric value controlling the severity of skewness for
#' asymmetric response distributions. Larger values produce stronger skewness.
#' Values close to 0 produce relatively mild skewness, whereas larger values
#' increasingly concentrate responses in the extreme categories. Typical values
#' range from approximately 0.2 (mild skewness) to 1.0 (strong skewness), with
#' the default value of 0.55 producing moderate skewness. Only used when
#' \code{shape} is \code{"Positive_skewness"} or
#' \code{"Negative_skewness"}.
#'
#' @param return_probs Logical. If \code{FALSE} (default), only the thresholds
#' are returned. If \code{TRUE}, a list containing both thresholds and category
#' probabilities is returned.
#'
#' @return
#' If \code{return_probs = FALSE}, a numeric vector containing threshold values
#' on the latent-response scale.
#'
#' If \code{return_probs = TRUE}, a list with:
#' \itemize{
#'   \item \code{thresholds}: threshold values on the latent-response scale.
#'   \item \code{probabilities}: category probabilities used to generate the thresholds.
#' }
#'
#' @details
#' The thresholds are generated from target category probabilities using the
#' inverse cumulative standard normal distribution. Symmetric thresholds are
#' based on predefined approximately symmetric response distributions, whereas
#' skewed thresholds are generated from exponentially decreasing category
#' probabilities controlled by \code{skew_strength}.
#'
#' The resulting thresholds can be supplied directly to the \code{thrA} and
#' \code{thrB} arguments of \code{\link{rho_correction}}.
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
#' # Return both thresholds and category probabilities
#' generate_thresholds(
#'   n_categories = 5,
#'   shape = "Negative_skewness",
#'   return_probs = TRUE
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



