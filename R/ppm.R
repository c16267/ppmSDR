#' Penalized Principal Machine Wrapper Function
#'
#' Provides a unified interface for fitting various penalized principal machine (P$^2$PM) estimators for sufficient dimension reduction (SDR) in regression and classification. The \code{ppm()} function calls the appropriate penalized SDR estimator based on the specified loss function, allowing users to select among squared loss, logistic, hinge, quantile, and their weighted versions for classification and regression.
#'
#' @param x Numeric matrix. The input data matrix (observations in rows, variables in columns).
#' @param y Numeric or factor vector. The response variable.
#' @param H Integer. Number of quantile slices or thresholds to use (default: 10).
#' @param C Numeric. Regularization parameter for the loss (default: 1).
#' @param loss Character. Specifies the loss type. One of \code{"lssvm"}, \code{"wlssvm"}, \code{"svm"}, \code{"wsvm"}, \code{"l2svm"}, \code{"wl2svm"}, \code{"logit"}, \code{"wlogit"}, \code{"asls"}, or \code{"qr"}.
#' @param penalty Character. Penalty type for group selection. One of \code{"grSCAD"}, \code{"grLasso"}, or \code{"grMCP"} (default: "grSCAD").
#' @param lambda Numeric. Regularization parameter for the penalty (default: 0.001).
#' @param gamma Numeric. Penalty shape parameter for SCAD/MCP (default: 3.7).
#' @param max.iter Integer. Maximum number of iterations (default: 100).
#' @param ... Additional arguments passed to the underlying method.
#'
#' @details
#' This wrapper selects among ten state-of-the-art penalized SDR estimators, including penalized principal least squares SVM (P$^2$LSM), penalized principal SVM (P$^2$SVM), penalized principal L2-hinge SVM (P$^2$L2M), penalized principal quantile regression (P$^2$QR), penalized principal asymmetric least squares (P$^2$AR), penalized principal logistic regression (P$^2$LR), and their weighted counterparts for classification.
#'
#' Each method is implemented in a dedicated function and uses an efficient optimization algorithm such as group coordinate descent (GCD) or majorization-minimization (MM) for computation-friendly and scalable estimation. The argument \code{loss} selects the corresponding method, and all relevant parameters are passed automatically.
#'
#' @return A list containing:
#'   \item{evectors}{Estimated basis directions (principal SDR vectors) as columns.}
#'   \item{evalues}{Associated eigenvalues.}
#'   \item{Mn}{Estimated working matrix (if returned by the underlying method).}
#'   \item{x}{Original input matrix.}
#'
#' @references
#'
#' Artemiou, A. and Dong, Y. (2016). Sufficient dimension reduction via principal lq support vector machine, \emph{Electronic Journal of Statistics}, 10: 783–805.
#'
#' Artemiou, A., Dong, Y. and Shin, S. J. (2021). Real-time sufficient dimension reduction through principal least squares support vector machines, \emph{Pattern Recognition}, 112: 107768.
#'
#' Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized likelihood and its oracle properties, \emph{Journal of the American Statistical Association}, 96: 1348–1360.
#'
#' Hunter, D. R. and Lange, K. (2004). A tutorial on MM algorithms, \emph{The American Statistician}, 58(1): 30–37.
#'
#' Jang, H. J., Shin, S. J. and Artemiou, A. (2023). Principal weighted least square support vector machine: An online dimension-reduction tool for binary classification, \emph{Computational Statistics & Data Analysis}, 187: 107818.
#'
#' Kim, B. and Shin, S. J. (2019). Principal weighted logistic regression for sufficient dimension reduction in binary classification, \emph{Journal of the Korean Statistical Society}, 48(2): 194–206.
#'
#' Li, B., Artemiou, A. and Li, L. (2011). Principal support vector machines for linear and nonlinear sufficient dimension reduction, \emph{Annals of Statistics}, 39(6): 3182–3210.
#'
#' Shin, J. and Shin, S. J. (2024). A concise overview of principal support vector machines and its generalization, \emph{Communications for Statistical Applications and Methods}, 31(2): 235–246.
#'
#' Shin, J., Shin, S. J. and Artemiou, A. (2024). The R package psvmsdr: A unified algorithm for sufficient dimension reduction via principal machines, \emph{arXiv preprint arXiv:2409.01547}.
#'
#' Shin, S. J. and Artemiou, A. (2017). Penalized principal logistic regression for sparse sufficient dimension reduction, \emph{Computational Statistics & Data Analysis}, 111: 48–58.
#'
#'
#'
#' @examples
#' \donttest{
#' library(MASS)
#' set.seed(1)
#' n <- 1000
#' p   <- 10
#' B <- matrix(0, p, 2)
#' B[1,1] <- B[2,2] <- 1
#' x <- mvrnorm(n, rep(0, p), diag(1,p))
#' y <- (x %*% B[,1]/(0.5 + (x %*% B[,2] + 1)^2)) + 0.2*rnorm(n, 0, 1)
#' y.binary <- sign(y)
#'
#' # SDR in Regression: lssvm
#' result1 <- ppm(x, y, H=10, C=1, loss="lssvm", penalty = 'grSCAD', lambda=0.003)
#' round(result1$evectors[,1:2], 5)
#'
#' # SDR in Regression: svm
#' result2 <- ppm(x, y, H=10, C=1, loss="svm", penalty = 'grSCAD', lambda=0.0001)
#' round(result2$evectors[,1:2], 5)

#' # SDR in Classification: wlssvm
#' result3 <- ppm(x, y.binary, H=10, C=1, loss="wlssvm", penalty = 'grSCAD', lambda=0.0005)
#' round(result3$evectors[,1:2], 5)
#'
#' # SDR in Classification: wlogit
#' result4 <- ppm(x, y.binary, H=10, C=20, loss="wlogit", penalty = 'grSCAD', lambda=0.07)
#' round(result4$evectors[,1:2], 5)
#'}
#'
#' @export
#' @import Matrix
#' @import grpreg
#' @import blockmatrix
#' @import Rfast
#' @importFrom stats cov glm model.matrix na.omit quantile relevel residuals rnorm
#' @importFrom Matrix bdiag
#' @importFrom methods as

ppm <- function(x, y, H = 10, C = 1, loss = "lssvm", penalty = "grSCAD",
                lambda = 0.001, gamma = 3.7, max.iter = 100, ...) {



  # Check for required arguments
  if (missing(x) || missing(y)) {
    stop("Both 'x' and 'y' must be provided.")
  }

  # Define loss_map: maps loss names to corresponding functions
  loss_map <- list(
    lssvm  = pplssvm,
    wlssvm = ppwlssvm,
    svm    = ppsvm,
    wsvm   = ppwsvm,
    l2svm  = ppl2svm,
    wl2svm = ppwl2svm,
    logit  = pplr,
    wlogit = ppwlr,
    asls   = ppasls,
    qr     = ppqr
  )
  loss <- tolower(loss)

  if (!(loss %in% names(loss_map))) {
    stop("Unknown loss function. Choose from: ", paste(names(loss_map), collapse = ", "))
  }

  y_vec <- as.vector(y)
  y_unique <- unique(y_vec)

  if (grepl("^w", loss)) {
    if (length(y_unique) != 2 || !all(sort(y_unique) == c(-1, 1))) {
      stop("For weighted methods (loss starting with 'w'), the response 'y' must be strictly binary with values +1 and -1.")
    }
  } else {
    if ((length(y_unique) == 2) &&
        (is.numeric(y_vec) || is.logical(y_vec) || is.integer(y_vec) || is.factor(y_vec))) {
      message("Note: For continuous loss functions, you provided a binary response y (two unique values: ",
              paste(y_unique, collapse=", "), "). Please check if this is intended.")
    }
  }

  call_args <- c(
    list(
      x = x, y = y, H = H, C = C,
      penalty = penalty, lambda = lambda,
      gamma = gamma, max.iter = max.iter
    ),
    list(...)
  )

  do.call(loss_map[[loss]], call_args)
}
