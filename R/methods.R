#' Print a fitted penalized principal machine
#'
#' @param x An object of class `"ppm"` returned by [ppm()].
#' @param ... Currently ignored.
#' @return The input object `x`, invisibly.
#' @seealso [ppm()], [summary.ppm()]
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200 * 6), 200, 6)
#' y <- x[, 1] / (0.5 + (x[, 2] + 1)^2) + 0.2 * rnorm(200)
#' fit <- ppm(x, y, loss = "lssvm", lambda = 0.01)
#' print(fit)
#' @export
print.ppm <- function(x, ...) {
  cat("Penalized Principal Machine (P2M) for Sufficient Dimension Reduction\n")
  cat("\nCall:\n  ")
  print(x$call)
  cat(sprintf("\nLoss: %s   Penalty: %s   lambda = %g   gamma = %g   C = %g   H = %d\n",
              x$loss, x$penalty, x$lambda, x$gamma, x$C, x$H))
  cat(sprintf("Data: n = %d, p = %d (%s response)\n", x$n, x$p, x$ytype))

  k  <- min(5L, length(x$evalues))
  ev <- round(x$evalues[seq_len(k)], 8L)
  cat("Leading eigenvalues of the working matrix M:\n  ",
      paste(ev, collapse = ", "), "\n", sep = "")
  invisible(x)
}


#' Summarize a fitted penalized principal machine
#'
#' Reports the estimated basis of the central subspace for a given working
#' dimension `d`, together with the predictors selected by the row-group
#' penalty (those whose loadings are non-zero across the leading `d`
#' directions).
#'
#' @param object An object of class `"ppm"` returned by [ppm()].
#' @param d Working structural dimension; the number of leading eigenvectors
#'   to report. Default `2`.
#' @param tol Tolerance below which a row L2-norm is treated as zero when
#'   determining the selected variables. Default `1e-6`.
#' @param ... Currently ignored.
#' @return Invisibly, a list with elements `d`, `basis` (the `p` by `d`
#'   estimated basis), `selected` (indices of the selected predictors) and
#'   `n.selected`.
#' @seealso [ppm()], [print.ppm()]
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200 * 6), 200, 6)
#' y <- x[, 1] / (0.5 + (x[, 2] + 1)^2) + 0.2 * rnorm(200)
#' fit <- ppm(x, y, loss = "lssvm", lambda = 0.01)
#' summary(fit, d = 2)
#' @export
summary.ppm <- function(object, d = 2, tol = 1e-6, ...) {
  if (!is.numeric(d) || length(d) != 1L || d < 1 || d != as.integer(d))
    stop("'d' must be a single positive integer.", call. = FALSE)
  d <- min(as.integer(d), object$p)

  basis <- object$evectors[, seq_len(d), drop = FALSE]
  rn <- if (!is.null(colnames(object$x))) {
    colnames(object$x)
  } else {
    paste0("x", seq_len(object$p))
  }
  rownames(basis) <- rn
  colnames(basis) <- paste0("Dir", seq_len(d))

  row_norm <- sqrt(rowSums(basis^2))
  selected <- which(row_norm > tol)

  cat("Penalized Principal Machine (P2M) summary\n")
  cat(sprintf("Loss: %s   Penalty: %s   lambda = %g\n",
              object$loss, object$penalty, object$lambda))
  cat(sprintf("Working dimension d = %d\n\n", d))
  cat("Estimated basis of the central subspace:\n")
  print(round(basis, 4L))
  cat(sprintf("\nSelected variables (%d of %d): %s\n",
              length(selected), object$p,
              if (length(selected)) paste(rn[selected], collapse = ", ") else "none"))

  invisible(list(d = d, basis = basis,
                 selected = selected, n.selected = length(selected)))
}


#' Print a penalized principal machine cross-validation result
#'
#' @param x An object of class `"ppm_tune"` returned by [ppm_tune()].
#' @param ... Currently ignored.
#' @return The input object `x`, invisibly.
#' @seealso [ppm_tune()]
#' @export
print.ppm_tune <- function(x, ...) {
  cat("Cross-validation for the penalized principal machine\n")
  cat(sprintf("Loss: %s   Penalty: %s   d = %d   folds = %d\n",
              x$loss, x$penalty, x$d, x$n.fold))
  cat(sprintf("Grid: %d candidates in [%.4g, %.4g]\n",
              length(x$lambda), min(x$lambda), max(x$lambda)))
  cat(sprintf("Selected lambda = %.5g  (mean dCor = %.4f)\n",
              x$opt.lambda, max(x$dcor, na.rm = TRUE)))
  invisible(x)
}
