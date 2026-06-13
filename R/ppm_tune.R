#' Cross-Validation for the Penalized Principal Machine Tuning Parameter
#'
#' Selects the sparsity parameter `lambda` of a penalized principal machine by
#' \eqn{K}-fold cross-validation. For each candidate `lambda`, the model is
#' fitted on the training folds and the held-out distance correlation (dCor;
#' Szekely, Rizzo and Bakirov, 2007) between the response and the predictors
#' projected onto the estimated `d`-dimensional central subspace is recorded.
#' The `lambda` maximizing the average held-out dCor is returned, following the
#' criterion of Shin, Wu, Zhang and Liu (2017).
#'
#' @inheritParams ppm
#' @param d Working structural dimension used to form the held-out sufficient
#'   predictors. Default `2`.
#' @param n.fold Number of cross-validation folds. Default `5`. For
#'   reproducible folds, call [set.seed()] before `ppm_tune()`.
#' @param lambda Optional numeric vector of candidate values. If `NULL`
#'   (default), a log-spaced grid of length `nlambda` is built from
#'   `lambda.max` down to `lambda.max * lambda.min.ratio`.
#' @param nlambda Length of the automatically generated grid. Default `20`.
#' @param lambda.max Largest candidate value in the generated grid. Default `1`.
#' @param lambda.min.ratio Ratio of the smallest to the largest candidate in
#'   the generated grid. Default `1e-3`.
#' @param verbose Logical; if `TRUE`, progress is reported via [message()].
#'   Default `FALSE`.
#'
#' @return An object of S3 class `"ppm_tune"`, a list with elements
#' \describe{
#'   \item{opt.lambda}{the selected value of `lambda`.}
#'   \item{lambda}{the candidate grid (decreasing).}
#'   \item{dcor}{the average held-out distance correlation at each candidate.}
#'   \item{fit}{a [ppm()] object refitted on all data at `opt.lambda`.}
#'   \item{d, loss, penalty, n.fold, call}{the settings and the matched call.}
#' }
#'
#' @references
#' Shin, S. J., Wu, Y., Zhang, H. H. and Liu, Y. (2017) Principal weighted
#' support vector machines for sufficient dimension reduction in binary
#' classification. *Biometrika*, 104(1), 67--81. \doi{10.1093/biomet/asw057}
#'
#' Szekely, G. J., Rizzo, M. L. and Bakirov, N. K. (2007) Measuring and testing
#' dependence by correlation of distances. *The Annals of Statistics*, 35(6),
#' 2769--2794. \doi{10.1214/009053607000000505}
#'
#' @seealso [ppm()]
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 1000; p <- 10
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] / (0.5 + (x[, 2] + 1)^2) + 0.2 * rnorm(n)
#' cv <- ppm_tune(x, y, loss = "lssvm", d = 2, n.fold = 5)
#' cv$opt.lambda
#' summary(cv$fit, d = 2)
#' }
#'
#' @export
ppm_tune <- function(x, y, loss = "lssvm", d = 2, H = 10, C = 1, gamma = 3.7,
                     penalty = c("grSCAD", "grLasso", "grMCP"),
                     max.iter = 100, n.fold = 5, lambda = NULL,
                     nlambda = 20, lambda.max = 0.1, lambda.min.ratio = 1e-3,
                     verbose = FALSE, ...) {

  cl <- match.call()

  loss_map <- c(lssvm = "pplssvm", wlssvm = "ppwlssvm",
                svm = "ppsvm",    wsvm  = "ppwsvm",
                l2svm = "ppl2svm", wl2svm = "ppwl2svm",
                logit = "pplr",   wlogit = "ppwlr",
                asls = "ppasls",  qr     = "ppqr")
  loss <- tolower(as.character(loss)[1L])
  if (!loss %in% names(loss_map))
    stop("Unknown loss '", loss, "'. Choose one of: ",
         paste(names(loss_map), collapse = ", "), ".", call. = FALSE)
  penalty <- match.arg(penalty)

  if (!is.numeric(d) || length(d) != 1L || d < 1 || d != as.integer(d))
    stop("'d' must be a single positive integer.", call. = FALSE)
  if (!is.numeric(n.fold) || length(n.fold) != 1L || n.fold < 2 ||
      n.fold != as.integer(n.fold))
    stop("'n.fold' must be a single integer greater than or equal to 2.",
         call. = FALSE)

  chk <- .validate_ppm_input(x, y, loss)
  x <- chk$x; y <- chk$y; n <- chk$n
  d <- min(as.integer(d), chk$p)

  ## candidate grid (decreasing) --------------------------------------------
  if (is.null(lambda)) {
    if (lambda.max <= 0 || lambda.min.ratio <= 0 || lambda.min.ratio >= 1)
      stop("Require lambda.max > 0 and 0 < lambda.min.ratio < 1.",
           call. = FALSE)
    grid <- exp(seq(log(lambda.max), log(lambda.max * lambda.min.ratio),
                    length.out = nlambda))
  } else {
    if (!is.numeric(lambda) || any(lambda <= 0))
      stop("'lambda' must be a vector of positive values.", call. = FALSE)
    grid <- sort(unique(lambda), decreasing = TRUE)
  }
  L <- length(grid)

  estimator <- get(loss_map[[loss]], mode = "function",
                   envir = environment(ppm_tune))
  fmls  <- names(formals(estimator))
  extra <- list(...)
  extra <- extra[names(extra) %in% fmls]

  ## random folds (uses the caller's RNG state; no seed set inside) ----------
  folds <- sample(rep_len(seq_len(n.fold), n))
  value <- matrix(NA_real_, n.fold, L)

  for (ll in seq_len(L)) {
    for (jj in seq_len(n.fold)) {
      ts <- which(folds == jj)
      tr <- setdiff(seq_len(n), ts)
      args <- list(x = x[tr, , drop = FALSE], y = y[tr], H = as.integer(H),
                   C = C, lambda = grid[ll], gamma = gamma, penalty = penalty,
                   max.iter = as.integer(max.iter))
      args <- c(args[names(args) %in% fmls], extra)
      B <- do.call(estimator, args)$evectors[, seq_len(d), drop = FALSE]
      scores <- x[ts, , drop = FALSE] %*% B
      value[jj, ll] <- .ppm_dcor(y[ts], scores)
    }
    if (verbose)
      message(sprintf("lambda = %-10.5g  mean dCor = %.4f",
                      grid[ll], mean(value[, ll], na.rm = TRUE)))
  }

  m <- colMeans(value, na.rm = TRUE)
  names(m) <- signif(grid, 4L)
  sel <- which.max(m)
  opt.lambda <- grid[sel]

  fit <- ppm(x, y, loss = loss, H = H, C = C, lambda = opt.lambda,
             gamma = gamma, penalty = penalty, max.iter = max.iter, ...)

  out <- list(opt.lambda = opt.lambda, lambda = grid, dcor = m, d = d,
              loss = loss, penalty = penalty, n.fold = as.integer(n.fold),
              fit = fit, call = cl)
  class(out) <- "ppm_tune"
  out
}


## Held-out dependence between the response and the sufficient predictors.
## Returns 0 for a degenerate (near-constant) projection so that over-shrunk
## lambdas are never selected.
.ppm_dcor <- function(yv, scores) {
  if (!is.matrix(scores)) scores <- as.matrix(scores)
  if (all(abs(scores) < .Machine$double.eps^0.5)) return(0)
  if (any(apply(scores, 2L, stats::sd) < .Machine$double.eps^0.5)) {
    scores <- scores[, apply(scores, 2L, stats::sd) >=
                       .Machine$double.eps^0.5, drop = FALSE]
    if (ncol(scores) == 0L) return(0)
  }
  val <- energy::dcor(yv, scores)
  if (!is.finite(val)) 0 else val
}
