#' Plot sensitivity analysis results
#'
#' Visualizes how the corrected association changes as a function of the assumed
#' measurement overlap.
#'
#' @param x A `"rho_overlap_sensitivity"` object.
#' @param show_ci Logical. If `TRUE`, displays confidence interval bands when available.
#' @param show_critical Logical. If `TRUE`, displays the critical overlap where the
#'   confidence interval includes zero.
#' @param show_observed Logical. If `TRUE`, displays the originally assumed overlap.
#' @param ... Further arguments (ignored).
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' \dontrun{
#' sens <- rho_overlap_sensitivity(...)
#' plot(sens)
#' }
#'
#' @export
plot.rho_overlap_sensitivity <- function(
    x,
    show_ci = TRUE,
    show_critical = TRUE,
    show_observed = TRUE,
    ...
) {
  if (!inherits(x, "rho_overlap_sensitivity")) {
    stop("x must be a rho_overlap_sensitivity object.")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }

  df <- x$results

  if (!all(c("overlap_value", "rho_hat") %in% names(df))) {
    stop("x$results must contain at least 'overlap_value' and 'rho_hat'.")
  }

  vary <- if (!is.null(x$vary)) x$vary else "multiplier"

  baseline_overlap_value <- if (!is.null(x$baseline_overlap_value)) {
    x$baseline_overlap_value
  } else {
    NA_real_
  }

  observed_overlap <- if (!is.null(x$observed_overlap)) x$observed_overlap else NA_real_
  critical_overlap <- if (!is.null(x$critical_overlap)) x$critical_overlap else NA_real_

  fmt <- function(z) {
    ifelse(is.finite(z), format(round(z, 3), nsmall = 0), NA_character_)
  }

  if (vary == "multiplier") {

    use_crossloading_labels <- is.finite(baseline_overlap_value)

    observed_label <- if (use_crossloading_labels) {
      fmt(observed_overlap * baseline_overlap_value)
    } else {
      fmt(observed_overlap)
    }

    critical_label <- if (is.finite(critical_overlap)) {
      if (use_crossloading_labels) {
        fmt(critical_overlap * baseline_overlap_value)
      } else {
        fmt(critical_overlap)
      }
    } else {
      NA_character_
    }

    x_axis_label <- "Overlap multiplier"

    subtitle_text <- if (use_crossloading_labels) {
      paste0(
        "Assumed overlap = ", observed_label,
        " (multiplier = ", fmt(observed_overlap), ")",
        if (!is.na(critical_label)) {
          paste0(" | Critical overlap = ", critical_label)
        } else {
          ""
        }
      )
    } else {
      paste0(
        "Assumed overlap multiplier = ", fmt(observed_overlap),
        if (!is.na(critical_label)) {
          paste0(" | Critical overlap multiplier = ", critical_label)
        } else {
          ""
        }
      )
    }

    title_text <- "Sensitivity to multiplying the assumed overlap"

  } else if (vary == "common") {

    observed_label <- fmt(observed_overlap)
    critical_label <- if (is.finite(critical_overlap)) fmt(critical_overlap) else NA_character_

    x_axis_label <- "Common overlap value"

    subtitle_text <- paste0(
      "Assumed common overlap value = ", observed_label,
      if (!is.na(critical_label)) {
        paste0(" | Critical overlap value = ", critical_label)
      } else {
        ""
      }
    )

    title_text <- "Sensitivity to a common overlap value"

  } else if (vary == "item") {

    observed_label <- fmt(observed_overlap)
    critical_label <- if (is.finite(critical_overlap)) fmt(critical_overlap) else NA_character_

    focal_text <- ""
    if (!is.null(x$focal_scale) && !is.null(x$focal_item)) {
      focal_text <- paste0(" for item ", x$focal_item, " in scale ", x$focal_scale)
    }

    x_axis_label <- paste0("Item-specific overlap value", focal_text)

    subtitle_text <- paste0(
      "Assumed focal overlap value = ", observed_label,
      if (!is.na(critical_label)) {
        paste0(" | Critical focal overlap value = ", critical_label)
      } else {
        ""
      }
    )

    title_text <- paste0("Sensitivity to one item-specific overlap value", focal_text)

  } else {

    x_axis_label <- "Overlap value"
    subtitle_text <- NULL
    title_text <- "Sensitivity of corrected association to assumed overlap"

  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$overlap_value, y = .data$rho_hat)
  )

  if (isTRUE(show_ci) &&
      all(c("CI_lower", "CI_upper") %in% names(df)) &&
      any(!is.na(df$CI_lower)) &&
      any(!is.na(df$CI_upper))) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$CI_lower, ymax = .data$CI_upper),
        alpha = 0.2
      )
  }

  p <- p +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::labs(
      x = x_axis_label,
      y = expression(hat(rho)),
      title = title_text,
      subtitle = subtitle_text
    ) +
    ggplot2::theme_minimal(base_size = 13)

  if (isTRUE(show_observed) &&
      is.finite(observed_overlap)) {

    p <- p +
      ggplot2::geom_vline(
        xintercept = observed_overlap,
        linetype = 3
      )

    obs_idx <- which.min(abs(df$overlap_value - observed_overlap))
    obs_df <- df[obs_idx, , drop = FALSE]

    p <- p +
      ggplot2::geom_point(
        data = obs_df,
        ggplot2::aes(x = .data$overlap_value, y = .data$rho_hat),
        inherit.aes = FALSE,
        size = 2.5
      )
  }

  if (isTRUE(show_critical) &&
      is.finite(critical_overlap)) {

    p <- p +
      ggplot2::geom_vline(
        xintercept = critical_overlap,
        linetype = 2
      )

    crit_idx <- which.min(abs(df$overlap_value - critical_overlap))
    crit_df <- df[crit_idx, , drop = FALSE]

    p <- p +
      ggplot2::geom_point(
        data = crit_df,
        ggplot2::aes(x = .data$overlap_value, y = .data$rho_hat),
        inherit.aes = FALSE,
        shape = 21,
        size = 2.5,
        stroke = 1
      )
  }

  p
}
