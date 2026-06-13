#' Penalized Principal Machine for Sufficient Dimension Reduction
#'
#' Fits a penalized principal machine (P2M), a sparse sufficient dimension
#' reduction (SDR) estimator, through a single group coordinate descent (GCD)
#' engine. A principal machine (PM) estimates the basis of the central subspace
#' by solving a family of convex-loss problems over several cutoffs (slices);
#' the penalized version adds a row-group sparsity penalty so that dimension
#' reduction and variable selection are performed simultaneously.
#'
#' @details
#' `ppm()` is a single front-end that dispatches to the loss-specific solver
#' selected by `loss`. Two families are supported, following the taxonomy of
#' Shin and Shin (2024):
#'
#' * **Response-based PM (RPM)** for a continuous response, where the loss is
#'   fixed and the pseudo-response varies across slices:
#'   `"lssvm"` (least squares, P2LSM), `"l2svm"` (L2-hinge, P2L2M),
#'   `"svm"` (hinge, P2SVM), `"logit"` (logistic, P2LR),
#'   `"asls"` (asymmetric least squares, P2AR) and `"qr"` (quantile, P2QR).
#' * **Loss-based PM (LPM)** for a binary response coded internally as
#'   \eqn{\{-1, +1\}}, where the loss varies across slices:
#'   `"wlssvm"` (P2WLSM), `"wl2svm"` (P2WL2M), `"wsvm"` (P2WSVM) and
#'   `"wlogit"` (P2WLR).
#'
#' Acronyms used above: SDR (sufficient dimension reduction), PM (principal
#' machine), P2M (penalized principal machine), GCD (group coordinate descent),
#' SVM (support vector machine). The penalty `penalty` is one of the group
#' LASSO (least absolute shrinkage and selection operator), the group SCAD
#' (smoothly clipped absolute deviation) or the group MCP (minimax concave
#' penalty), passed as `"grLasso"`, `"grSCAD"` or `"grMCP"`.
#'
#' The basis of the central subspace is estimated by the leading eigenvectors
#' of the working matrix \eqn{M = \sum_{k=1}^{H} \beta_k \beta_k^\top}, where
#' \eqn{\beta_k} is the slope estimated at the \eqn{k}-th cutoff.
#'
#' @param x A numeric matrix or data frame of predictors, of dimension
#'   `n` (observations) by `p` (variables).
#' @param y A response vector of length `n`. For continuous-type losses a
#'   numeric vector is expected; for weighted losses (those whose name starts
#'   with `"w"`) a two-class response is expected and is recoded internally to
#'   \eqn{\{-1, +1\}}.
#' @param loss Character string selecting the loss function. One of
#'   `"lssvm"`, `"wlssvm"`, `"svm"`, `"wsvm"`, `"l2svm"`, `"wl2svm"`,
#'   `"logit"`, `"wlogit"`, `"asls"`, `"qr"`. Default `"lssvm"`.
#' @param H Number of cutoffs (slices); a single integer `>= 2`. Default `10`.
#' @param C Positive cost parameter that balances the loss against the
#'   covariance term. Default `1`.
#' @param lambda Positive regularization parameter controlling sparsity.
#'   Default `0.1`. In practice `lambda` should be selected by cross-validation.
#' @param gamma Concavity parameter of the SCAD/MCP penalty; must exceed `2`.
#'   Default `3.7`.
#' @param penalty Penalty type: `"grSCAD"` (default), `"grLasso"` or `"grMCP"`.
#' @param max.iter Maximum number of GCD iterations. Default `100`.
#' @param ... Additional arguments passed to the underlying solver. The most
#'   useful is `ridge`, a small non-negative ridge constant added for numerical
#'   stability, which is accepted by the iterative solvers
#'   (`logit`, `wlogit`, `svm`, `wsvm`, `qr`, `lssvm`, `wlssvm`, `wl2svm`).
#'
#' @return An object of S3 class `"ppm"`, a list containing:
#' \describe{
#'   \item{M}{the estimated working matrix (a `p` by `p` symmetric matrix).}
#'   \item{evalues, evectors}{the eigenvalues and eigenvectors of `M`; the
#'     leading `d` eigenvectors estimate the basis of the central subspace.}
#'   \item{x, y}{the (validated) input data.}
#'   \item{loss, penalty, lambda, gamma, C, H, max.iter}{the fitting settings.}
#'   \item{ytype}{`"continuous"` or `"binary"`.}
#'   \item{n, p}{the sample size and the number of predictors.}
#'   \item{call}{the matched call.}
#' }
#'
#' @references
#' Li, B., Artemiou, A. and Li, L. (2011)
#' Principal support vector machines for linear and nonlinear sufficient
#' dimension reduction. *The Annals of Statistics*, 39(6), 3182--3210.
#' \doi{10.1214/11-AOS932}
#'
#' Shin, S. J. and Artemiou, A. (2017)
#' Penalized principal logistic regression for sparse sufficient dimension
#' reduction. *Computational Statistics & Data Analysis*, 111, 48--58.
#' \doi{10.1016/j.csda.2016.12.003}
#'
#' Breheny, P. and Huang, J. (2015)
#' Group descent algorithms for nonconvex penalized linear and logistic
#' regression models with grouped predictors. *Statistics and Computing*,
#' 25, 173--187. \doi{10.1007/s11222-013-9424-2}
#'
#' @seealso [print.ppm()], [summary.ppm()]
#'
#' @examples
#' set.seed(1)
#' n <- 1000; p <- 10
#' B <- matrix(0, p, 2); B[1, 1] <- B[2, 2] <- 1
#' x <- matrix(rnorm(n * p), n, p)
#' y <- (x %*% B[, 1]) / (0.5 + (x %*% B[, 2] + 1)^2) + 0.2 * rnorm(n)
#'
#' ## penalized principal least-squares SVM (P2LSM) with the group SCAD penalty
#' fit <- ppm(x, y, loss = "lssvm", penalty = "grSCAD", lambda = 0.01)
#' round(fit$evectors[, 1:2], 3)
#' print(fit)
#' summary(fit)
#'
#' \donttest{
#' ## binary response: penalized principal asymmetric least squares (P2AR)
#' yb <- sign(y)
#' fitw <- ppm(x, yb, loss = "asls", penalty = "grSCAD", lambda = 0.04)
#' round(fitw$evectors[, 1:2], 3)
#' print(fitw)
#' summary(fitw)
#' }
#' \donttest{
#' data(boston)
#' xb <- scale(as.matrix(boston[, setdiff(names(boston), "medv")]))
#' yb <- as.numeric(scale(boston$medv))
#' fit_b <- ppm(xb, yb, loss = "lssvm", penalty = "grSCAD", lambda = 8e-3)
#' summary(fit_b, d = 2)
#' }
#'
#' @export
ppm <- function(x, y, loss = "lssvm", H = 10, C = 1, lambda = 0.01,
                gamma = 3.7, penalty = c("grSCAD", "grLasso", "grMCP"),
                max.iter = 100, ...) {

  cl <- match.call()

  loss_map <- c(lssvm = "pplssvm", wlssvm = "ppwlssvm",
                svm = "ppsvm",    wsvm  = "ppwsvm",
                l2svm = "ppl2svm", wl2svm = "ppwl2svm",
                logit = "pplr",   wlogit = "ppwlr",
                asls = "ppasls",  qr     = "ppqr")

  loss <- tolower(as.character(loss)[1L])
  if (!loss %in% names(loss_map)) {
    stop("Unknown loss '", loss, "'. Choose one of: ",
         paste(names(loss_map), collapse = ", "), ".", call. = FALSE)
  }
  penalty <- match.arg(penalty)

  ## ---- scalar hyper-parameter checks ------------------------------------
  .pos_scalar <- function(v, nm) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0)
      stop("'", nm, "' must be a single positive number.", call. = FALSE)
  }
  .pos_scalar(C, "C")
  .pos_scalar(lambda, "lambda")
  .pos_scalar(gamma, "gamma")
  if (gamma <= 2)
    stop("'gamma' must be greater than 2 for the SCAD/MCP penalties.",
         call. = FALSE)
  if (!is.numeric(H) || length(H) != 1L || H < 2 || H != as.integer(H))
    stop("'H' must be a single integer greater than or equal to 2.",
         call. = FALSE)
  if (!is.numeric(max.iter) || length(max.iter) != 1L ||
      max.iter < 1 || max.iter != as.integer(max.iter))
    stop("'max.iter' must be a single positive integer.", call. = FALSE)

  ## ---- data validation / response coding --------------------------------
  chk <- .validate_ppm_input(x, y, loss)

  ## ---- dispatch ---------------------------------------------------------
  estimator <- get(loss_map[[loss]], mode = "function",
                   envir = environment(ppm))
  args <- list(x = chk$x, y = chk$y, H = as.integer(H), C = C,
               lambda = lambda, gamma = gamma, penalty = penalty,
               max.iter = as.integer(max.iter))
  extra <- list(...)
  fmls  <- names(formals(estimator))
  args  <- c(args[names(args) %in% fmls],
             extra[names(extra) %in% fmls])

  fit <- do.call(estimator, args)

  out <- list(
    M = fit$M, evalues = fit$evalues, evectors = fit$evectors,
    x = chk$x, y = chk$y,
    loss = loss, penalty = penalty, lambda = lambda, gamma = gamma,
    C = C, H = as.integer(H), max.iter = as.integer(max.iter),
    ytype = chk$ytype, n = chk$n, p = chk$p, call = cl
  )
  class(out) <- "ppm"
  out
}


## Internal input validator (not exported).
.validate_ppm_input <- function(x, y, loss) {
  if (missing(x) || missing(y))
    stop("Both 'x' and 'y' must be provided.", call. = FALSE)

  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x))     x <- as.matrix(x)
  if (!is.numeric(x))
    stop("'x' must be a numeric matrix or data frame.", call. = FALSE)
  if (anyNA(x))
    stop("'x' contains missing values; remove or impute them first.",
         call. = FALSE)
  if (any(!is.finite(x)))
    stop("'x' contains non-finite values (Inf or NaN).", call. = FALSE)

  n <- nrow(x); p <- ncol(x)
  if (n < 2L) stop("'x' must have at least two rows.", call. = FALSE)
  if (p < 1L) stop("'x' must have at least one column.", call. = FALSE)

  if (length(y) != n)
    stop(sprintf("length(y) = %d does not match nrow(x) = %d.",
                 length(y), n), call. = FALSE)
  if (anyNA(y))
    stop("'y' contains missing values.", call. = FALSE)

  weighted <- grepl("^w", loss)
  uy <- unique(y)

  if (weighted) {
    if (length(uy) != 2L)
      stop("Weighted losses (here '", loss,
           "') require a binary response with exactly two classes.",
           call. = FALSE)
    lev  <- sort(uy)
    ymap <- ifelse(y == lev[2L], 1, -1)
    if (!isTRUE(all.equal(as.numeric(lev), c(-1, 1))))
      message("Note: binary response recoded as '", lev[1L],
              "' -> -1 and '", lev[2L], "' -> +1.")
    y <- ymap
    ytype <- "binary"
  } else {
    if (is.factor(y) || is.character(y))
      stop("Loss '", loss, "' expects a numeric response. ",
           "For binary classification use a weighted loss ",
           "(e.g. 'wsvm', 'wlssvm', 'wlogit', 'wl2svm').", call. = FALSE)
    y <- as.numeric(y)
    ytype <- if (length(uy) <= 2L) "binary" else "continuous"
    if (ytype == "binary"){}
  }

  list(x = x, y = y, n = n, p = p, ytype = ytype, weighted = weighted)
}
