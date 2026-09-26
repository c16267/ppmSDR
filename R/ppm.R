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
#' @section Algorithms and the step-halving line search:
#' The solvers fall into three classes, and the class determines the default
#' of `line.search`:
#'
#' * **GCD** (`"lssvm"`, `"wlssvm"`): the squared loss gives an exact group
#'   penalized least squares problem, solved by a single GCD run. There is no
#'   outer iteration and `line.search` does not apply (it is ignored with a
#'   warning if set to `TRUE`).
#' * **Iterative GCD** (`"logit"`, `"wlogit"`, `"asls"`, `"l2svm"`,
#'   `"wl2svm"`): the loss is replaced at every outer iteration by a local
#'   quadratic approximation (IRLS for the logistic loss; residual signs or
#'   the active set frozen for the asymmetric squared and squared hinge
#'   losses). Such an approximation is not a majorizer, so a plain GCD update
#'   does not by itself guarantee descent of the penalized objective
#'   \eqn{Q(\theta)}. For `"logit"`, `"wlogit"`, `"l2svm"` and `"wl2svm"`
#'   the default `line.search = TRUE` therefore pairs every update with a
#'   step-halving safeguard. For `"asls"` (P2AR) the default is
#'   `line.search = FALSE`, i.e. the plain updates used in the numerical
#'   studies of the accompanying article; the safeguard can be switched on,
#'   and it is advisable to do so when a large cost parameter is used (the
#'   plain sweep is stable for moderate values, roughly `C <= 5`; for larger
#'   `C` it can diverge and the solver then stops with an informative error).
#' * **MM-GCD** (`"svm"`, `"wsvm"`, `"qr"`): the hinge and check losses are
#'   replaced by quadratic majorizers, so descent follows from the
#'   majorization-minimization argument and the default is
#'   `line.search = FALSE`. The safeguard can still be switched on. The
#'   majorizers use a floor `eps` that safeguards a zero margin (hinge) or a
#'   zero residual (check loss). For `"svm"` and `"wsvm"` the margin is
#'   scale-free and `eps = 1e-6`. For `"qr"` the residual has the scale of
#'   the response and, since the quantile-regression solution interpolates
#'   `p + 1` observations exactly, a tiny absolute floor makes the MM weights
#'   `1 / (4 eps)` explode as the iterations approach the solution; the
#'   default is therefore relative, `eps = 0.05 * sd(y)`, which keeps the
#'   iterations stable and can be overridden through `...`.
#'
#' With `line.search = TRUE`, if \eqn{\hat\theta^{(t)}} denotes the candidate
#' returned by the GCD step at iteration \eqn{t}, the accepted update is
#' \deqn{\theta^{(t+1)} = \theta^{(t)} + 2^{-m_t}\{\hat\theta^{(t)} - \theta^{(t)}\},}
#' where \eqn{m_t \ge 0} is the smallest integer for which
#' \eqn{Q(\theta^{(t+1)}) \le Q(\theta^{(t)})}, and \eqn{Q} is the original
#' penalized P2M objective
#' \deqn{Q(\theta) = \sum_{k=1}^{h}\Big[\theta_k^\top\tilde\Sigma\theta_k
#'   + \frac{c}{n}\sum_{i=1}^{n} w_{ik} L_k(\tilde y_{ik}, \theta_k^\top\tilde x_i)\Big]
#'   + \sum_{j=1}^{p} p_\lambda(\|\theta_{(j)}\|),}
#' evaluated on the objective itself and not on the quadratic surrogate. If no
#' acceptable step is found within `max.halving` halvings, the current iterate
#' is retained and the algorithm stops. The returned `status` is then
#' `"converged"` when the rejected candidate step already satisfied the
#' stopping rule below (the plain iteration would have stopped there as
#' well), and `"line.search.halted"` otherwise, meaning that the GCD
#' direction failed to decrease \eqn{Q} at a non-converged iterate and the
#' last accepted iterate is returned. No warning is issued; inspect `status`,
#' `iter` and `n.halving` (also shown by [print.ppm()]) to see whether and
#' where the safeguard intervened.
#' The safeguard guarantees monotone non-increase of \eqn{Q} along the
#' accepted iterates; it coincides with the plain update (\eqn{m_t = 0}) at
#' every iteration at which the plain update does not increase \eqn{Q}, so
#' `line.search = FALSE` reproduces the plain iterative GCD exactly.
#'
#' All iterative solvers start at \eqn{\theta^{(0)} = 0} and stop once
#' \eqn{\|\theta^{(t+1)} - \theta^{(t)}\|_\infty / (1 + \|\theta^{(t)}\|_\infty) < 10^{-5}}
#' or after `max.iter` outer iterations.
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
#'   Default `0.01`. In practice `lambda` should be selected by cross-validation.
#' @param gamma Concavity parameter of the SCAD/MCP penalty; must exceed `2`.
#'   Default `3.7`.
#' @param penalty Penalty type: `"grSCAD"` (default), `"grLasso"` or `"grMCP"`.
#' @param max.iter Maximum number of outer (iterative GCD / MM-GCD) iterations.
#'   Default `100`.
#' @param line.search Logical or `NULL`. Whether to apply the step-halving
#'   line search on the original penalized objective described in the
#'   *Algorithms* section. The default `NULL` selects `TRUE` for the
#'   iterative GCD losses `"logit"`, `"wlogit"`, `"l2svm"` and `"wl2svm"`,
#'   and `FALSE` for all other losses (`"asls"`, the MM-GCD losses `"svm"`,
#'   `"wsvm"`, `"qr"`, and the exact GCD losses). Setting `FALSE` for an
#'   iterative GCD loss gives the plain (undamped) updates used in the
#'   numerical studies of the accompanying article; setting `TRUE` for
#'   `"asls"` or for an MM-GCD loss adds the safeguard on top of the
#'   corresponding quadratic step. The argument is ignored, with a warning,
#'   for the exact GCD losses `"lssvm"` and `"wlssvm"`.
#' @param max.halving Maximum number of step halvings \eqn{m_{\max}} tried at
#'   each iteration when `line.search = TRUE`. Default `10`.
#' @param ... Additional arguments passed to the underlying solver. The most
#'   useful are `ridge`, a small non-negative ridge constant added for
#'   numerical stability, which is accepted by the iterative solvers
#'   (`logit`, `wlogit`, `svm`, `wsvm`, `qr`, `lssvm`, `wlssvm`, `wl2svm`),
#'   and `eps`, the positive floor of the quadratic majorizer for the MM-GCD
#'   losses (`svm`, `wsvm`: default `1e-6`; `qr`: default `0.05 * sd(y)`, see
#'   the *Algorithms* section).
#'
#' @return An object of S3 class `"ppm"`, a list containing:
#' \describe{
#'   \item{M}{the estimated working matrix (a `p` by `p` symmetric matrix).}
#'   \item{evalues, evectors}{the eigenvalues and eigenvectors of `M`; the
#'     leading `d` eigenvectors estimate the basis of the central subspace.}
#'   \item{theta}{the estimated coefficient matrix \eqn{\Theta} of dimension
#'     `(p + 1)` by `H - 1`, whose `k`-th column holds the intercept (first
#'     row) and the slope \eqn{\beta_k} of the `k`-th machine.}
#'   \item{x, y}{the (validated) input data.}
#'   \item{loss, penalty, lambda, gamma, C, H, max.iter}{the fitting settings.}
#'   \item{algorithm}{`"GCD"`, `"iterative GCD"` or `"MM-GCD"`.}
#'   \item{line.search, max.halving}{the line-search settings actually used.}
#'   \item{iter}{the number of outer iterations performed.}
#'   \item{status}{`"converged"`, `"max.iter"` or `"line.search.halted"`.}
#'   \item{objective}{when `line.search = TRUE`, the value of the penalized
#'     objective \eqn{Q(\theta^{(t)})} after each accepted iteration
#'     (non-increasing by construction); otherwise `NULL`.}
#'   \item{n.halving}{when `line.search = TRUE`, the number of halvings
#'     \eqn{m_t} at each iteration (`NA` at a halted iteration);
#'     otherwise `NULL`.}
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
#' ## penalized principal logistic regression (P2LR): iterative GCD, so the
#' ## step-halving line search is on by default
#' fit_lr <- ppm(x, y, loss = "logit", penalty = "grSCAD", lambda = 0.01)
#' fit_lr$status
#' fit_lr$n.halving            # m_t at each iteration (0 = plain update)
#' all(diff(fit_lr$objective) <= 0)   # Q is non-increasing
#'
#' ## the plain (undamped) updates used in the article's numerical studies
#' fit_lr0 <- ppm(x, y, loss = "logit", penalty = "grSCAD", lambda = 0.01,
#'                line.search = FALSE)
#'
#' \donttest{
#' ## binary response with a two-dimensional central subspace spanned by
#' ## (x1, x2): penalized principal asymmetric least squares (P2AR), plain
#' ## updates by default, and penalized principal weighted logistic regression
#' yb <- sign(x[, 1] + x[, 2]^3 / 3 + 0.2 * rnorm(n))
#' fitw <- ppm(x, yb, loss = "asls", penalty = "grSCAD", lambda = 0.2)
#' round(fitw$evectors[, 1:2], 3)
#' print(fitw)
#' summary(fitw)
#' fit_wlr <- ppm(x, yb, loss = "wlogit", penalty = "grSCAD", lambda = 0.005)
#' summary(fit_wlr)
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
                max.iter = 100, line.search = NULL, max.halving = 10, ...) {

  cl <- match.call()

  loss_map <- .ppm_loss_map()

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

  ## ---- line-search settings ----------------------------------------------
  algorithm   <- .ppm_algorithm(loss)
  line.search <- .ppm_resolve_line_search(loss, line.search)
  max.halving <- .ppm_check_max_halving(max.halving)

  ## ---- data validation / response coding --------------------------------
  chk <- .validate_ppm_input(x, y, loss)

  ## ---- dispatch ---------------------------------------------------------
  estimator <- get(loss_map[[loss]], mode = "function",
                   envir = environment(ppm))
  args <- list(x = chk$x, y = chk$y, H = as.integer(H), C = C,
               lambda = lambda, gamma = gamma, penalty = penalty,
               max.iter = as.integer(max.iter),
               line.search = line.search, max.halving = max.halving)
  extra <- list(...)
  fmls  <- names(formals(estimator))
  args  <- c(args[names(args) %in% fmls],
             extra[names(extra) %in% fmls])

  fit <- do.call(estimator, args)

  out <- list(
    M = fit$M, evalues = fit$evalues, evectors = fit$evectors,
    theta = fit$theta,
    x = chk$x, y = chk$y,
    loss = loss, penalty = penalty, lambda = lambda, gamma = gamma,
    C = C, H = as.integer(H), max.iter = as.integer(max.iter),
    algorithm = algorithm,
    line.search = isTRUE(fit$line.search), max.halving = max.halving,
    iter = fit$iter, status = fit$status,
    objective = fit$objective, n.halving = fit$n.halving,
    ytype = chk$ytype, n = chk$n, p = chk$p, call = cl
  )
  class(out) <- "ppm"
  out
}


## ---------------------------------------------------------------------------
## Internal helpers shared by ppm() and ppm_tune() (not exported).
## ---------------------------------------------------------------------------

## loss name -> internal solver
.ppm_loss_map <- function() {
  c(lssvm = "pplssvm", wlssvm = "ppwlssvm",
    svm = "ppsvm",    wsvm  = "ppwsvm",
    l2svm = "ppl2svm", wl2svm = "ppwl2svm",
    logit = "pplr",   wlogit = "ppwlr",
    asls = "ppasls",  qr     = "ppqr")
}

## Algorithm class of a loss, following Table 1 of the article.
.ppm_algorithm <- function(loss) {
  if (loss %in% c("lssvm", "wlssvm")) return("GCD")
  if (loss %in% c("logit", "wlogit", "asls", "l2svm", "wl2svm"))
    return("iterative GCD")
  if (loss %in% c("svm", "wsvm", "qr")) return("MM-GCD")
  stop("Unknown loss '", loss, "'.", call. = FALSE)
}

## Default of line.search when the user leaves it NULL: the step-halving
## safeguard is on for the iterative GCD losses that use a local quadratic
## approximation of a smooth or saturating loss (logit, wlogit, l2svm,
## wl2svm); it is off for asls (plain updates, as in the article's numerical
## studies), for the MM-GCD losses (descent follows from majorization) and
## for the exact GCD losses (no outer iteration).
.ppm_line_search_default <- function(loss) {
  loss %in% c("logit", "wlogit", "l2svm", "wl2svm")
}

## Resolve the line.search argument:
##   NULL  -> loss-dependent default (see .ppm_line_search_default);
##   TRUE  -> allowed for iterative GCD and MM-GCD; ignored (warning) for GCD.
.ppm_resolve_line_search <- function(loss, line.search) {
  alg <- .ppm_algorithm(loss)
  if (is.null(line.search)) return(.ppm_line_search_default(loss))
  if (!is.logical(line.search) || length(line.search) != 1L ||
      is.na(line.search))
    stop("'line.search' must be NULL, TRUE or FALSE.", call. = FALSE)
  if (line.search && alg == "GCD") {
    warning("'line.search' does not apply to loss '", loss,
            "' (exact GCD on a single quadratic problem); it is ignored.",
            call. = FALSE)
    return(FALSE)
  }
  line.search
}

.ppm_check_max_halving <- function(max.halving) {
  if (!is.numeric(max.halving) || length(max.halving) != 1L ||
      !is.finite(max.halving) || max.halving < 0 ||
      max.halving != as.integer(max.halving))
    stop("'max.halving' must be a single non-negative integer.", call. = FALSE)
  as.integer(max.halving)
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
  }

  list(x = x, y = y, n = n, p = p, ytype = ytype, weighted = weighted)
}
