#' @export
print.latentoverlap_result <- function(x, digits = 3, ...) {
  cat("\n")
  cat("Latent overlap correction result\n")
  cat("--------------------------------\n")

  if (!is.null(x$type)) {
    cat("Type:   ", x$type, "\n", sep = "")
  }
  if (!is.null(x$model)) {
    cat("Model:  ", x$model, "\n", sep = "")
  }
  if (!is.null(x$method)) {
    cat("Method: ", x$method, "\n", sep = "")
  }

  cat("\n")
  cat("Corrected association (rho_hat): ",
      formatC(x$rho_hat, digits = digits, format = "f"), "\n", sep = "")

  if (!is.null(x$CI_lower) && !is.null(x$CI_upper)) {
    cl_pct <- if (!is.null(x$conf_level)) round(100 * x$conf_level, 1) else 95
    cat(cl_pct, "% CI: [",
        formatC(x$CI_lower, digits = digits, format = "f"), ", ",
        formatC(x$CI_upper, digits = digits, format = "f"), "]\n", sep = "")
  }

  if (!is.null(x$SE_corrected)) {
    cat("SE(corrected): ",
        formatC(x$SE_corrected, digits = digits, format = "f"), "\n", sep = "")
  }

  if (!is.null(x$boundary) && isTRUE(x$boundary)) {
    cat("\n")
    cat("Point estimate hit boundary")
    if (!is.null(x$boundary_side) && !is.na(x$boundary_side)) {
      cat(" (", x$boundary_side, ")", sep = "")
    }
    cat(".\n")
  }

  if (!is.null(x$prop_boundary_boot) && is.finite(x$prop_boundary_boot) &&
      x$prop_boundary_boot > 0) {
    cat("Bootstrap draws hitting boundary: ",
        formatC(100 * x$prop_boundary_boot, digits = digits, format = "f"),
        "%", sep = "")
    if (!is.null(x$n_boundary_boot) && !is.null(x$n_boot)) {
      cat(" (", x$n_boundary_boot, "/", x$n_boot, ")", sep = "")
    }
    cat("\n")
  }

  if (!is.null(x$feasible_r_min) && !is.null(x$feasible_r_max) &&
      is.finite(x$feasible_r_min) && is.finite(x$feasible_r_max)) {
    cat("Feasible observed r range under model: [",
        formatC(x$feasible_r_min, digits = digits, format = "f"), ", ",
        formatC(x$feasible_r_max, digits = digits, format = "f"), "]\n", sep = "")
  }

  cat("\n")
  invisible(x)
}
