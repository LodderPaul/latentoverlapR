#' @export
print.latentoverlap_fit <- function(x, ...) {
  cls <- class(x)
  on.exit(class(x) <- cls, add = TRUE)
  class(x) <- "list"

  # Put rho_hat first; keep the rest in original order
  nm <- names(x)
  if (!is.null(nm) && "rho_hat" %in% nm) {
    nm2 <- c("rho_hat", setdiff(nm, "rho_hat"))
    x <- x[nm2]
  }

  print(x, ...)
  invisible(x)
}
