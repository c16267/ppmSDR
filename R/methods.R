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

  ## algorithm / convergence line -------------------------------------------
  alg <- if (is.null(x$algorithm)) "GCD" else x$algorithm
  if (alg == "GCD") {
    cat("Algorithm: GCD (exact; no line search)\n")
  } else {
    ls_txt <- if (isTRUE(x$line.search)) {
      nh <- x$n.halving
      nh <- nh[!is.na(nh)]
      sprintf("step-halving line search on (max.halving = %d; %d of %d updates damped)",
              x$max.halving, sum(nh > 0), length(x$n.halving))
    } else {
      "line search off"
    }
    cat(sprintf("Algorithm: %s, %s\n", alg, ls_txt))
    if (!is.null(x$iter)) {
      st <- if (is.null(x$status)) "" else x$status
      if (identical(st, "line.search.halted"))
        st <- sprintf("line.search.halted: no descent step at iteration %s, last accepted iterate returned", x$iter)
      cat(sprintf("Iterations: %s (%s)\n",
                  if (is.na(x$iter)) "NA" else as.character(x$iter), st))
    }
  }

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
  if (!is.null(object$algorithm) && object$algorithm != "GCD")
    cat(sprintf("Algorithm: %s (line search %s), %s iterations, %s\n",
                object$algorithm,
                if (isTRUE(object$line.search)) "on" else "off",
                if (is.null(object$iter) || is.na(object$iter)) "NA" else object$iter,
                if (is.null(object$status)) "" else object$status))
  cat(sprintf("Working dimension d = %d\n\n", d))
  cat("Estimated basis of the central subspace:\n")
  print(round(basis, 4L))
  sel_txt <- sprintf("Selected variables (%d of %d): %s",
                     length(selected), object$p,
                     if (length(selected)) paste(rn[selected], collapse = ", ") else "none")
  cat("\n", paste(strwrap(sel_txt, width = getOption("width"), exdent = 2),
                  collapse = "\n"), "\n", sep = "")

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
  if (!is.null(x$line.search))
    cat(sprintf("Line search: %s\n", if (isTRUE(x$line.search)) "on" else "off"))
  cat(sprintf("Grid: %d candidates in [%.4g, %.4g]\n",
              length(x$lambda), min(x$lambda), max(x$lambda)))
  cat(sprintf("Selected lambda = %.5g  (mean dCor = %.4f)\n",
              x$opt.lambda, max(x$dcor, na.rm = TRUE)))
  invisible(x)
}
