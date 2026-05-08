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
                                shape = c("Symmetric", "Positive_skewness", "Negative_skewness")) {
  shape <- match.arg(shape)
  K <- as.integer(n_categories)

  thr_table <- list(
    `2` = list(
      Symmetric           = c( 0.000),
      Positive_skewness   = c( 1.052),
      Negative_skewness   = c(-1.052)
    ),
    `3` = list(
      Symmetric           = c(-0.842,  0.842),
      Positive_skewness   = c( 0.994,  2.721),
      Negative_skewness   = c(-2.721, -0.994)
    ),
    `4` = list(
      Symmetric           = c(-1.282,  0.000,  1.282),
      Positive_skewness   = c( 0.925,  2.376,  3.693),
      Negative_skewness   = c(-3.693, -2.376, -0.925)
    ),
    `5` = list(
      Symmetric           = c(-1.645, -0.643, 0.643,   1.645),
      Positive_skewness   = c( 0.822,  1.944,  3.248,  4.312),
      Negative_skewness   = c(-4.312, -3.248, -1.944, -0.822)
    ),
    `6` = list(
      Symmetric           = c(-1.645, -0.842,  0.000,  0.842,  1.645),
      Positive_skewness   = c( 0.750,  1.817,  2.776,  3.690,  4.525),
      Negative_skewness   = c(-4.525, -3.690, -2.776, -1.817, -0.750)
    ),
    `7` = list(
      Symmetric           = c(-1.881, -1.175, -0.524,  0.524, 1.175, 1.881),
      Positive_skewness   = c( 0.623,  1.501,  2.250,  3.066,  3.795,  4.430),
      Negative_skewness   = c(-4.430, -3.795, -3.066, -2.250, -1.501, -0.623)
    )
  )

  key <- as.character(K)
  if (!key %in% names(thr_table)) {
    stop("generate_thresholds: implemented only for 2-7 categories.")
  }

  thr_table[[key]][[shape]]
}
