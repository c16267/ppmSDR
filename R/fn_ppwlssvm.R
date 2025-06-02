#' Penalized Principal Weighted Least Squares SVM (P\eqn{^2}WLSM) for Sparse Sufficient Dimension Reduction
#'
#' Implements the penalized principal weighted least squares support vector machine (P\eqn{^2}WLSM) for sparse sufficient dimension reduction (SDR)
#'
#' @param x Numeric predictor matrix (n x p), where n is the sample size and p is the number of predictors.
#' @param y Binary response vector of length n, coded as \code{-1, 1}.
#' @param H Number of quantile slices for principal machine construction (default: 10).
#' @param C Regularization parameter.
#' @param lambda Penalty parameter for sparsity. If NULL, cross-validation or an information criterion will be used for selection.
#' @param gamma Regularization parameter for SCAD/MCP penalty (default: 3.7).
#' @param penalty Penalty type: \code{"grSCAD"}, \code{"grLasso"}, or \code{"grMCP"}.
#' @param max.iter Maximum number of iterations (default: 100).
#'
#' @details
#' This function fits a penalized principal machine using a weighted least squares SVM loss, suitable for binary responsess. It supports group penalties for structured sparsity and variable selection. Only use with binary responses (\code{-1, 1}).
#'
#' @return
#' A list with components:
#'   \item{evalues}{Eigenvalues of the estimated sufficient dimension reduction matrix.}
#'   \item{evectors}{Eigenvectors spanning the SDR subspace.}
#'   \item{x}{(Optional) Original predictor matrix.}
#'
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   n <- 1000; p <- 10
#'   x <- MASS::mvrnorm(n, rep(0, p), diag(1, p))
#'   y <- (x %*% B[,1] / (0.5 + (x %*% B[,2] + 1)^2)) + 0.2 * rnorm(n)
#'   y.binary <- sign(y)
#'   fit <- ppwlssvm(x, y.binary, H = 10, C = 1, lambda = 0.0003, gamma = 3.7,
#'                   penalty = "grSCAD", max.iter = 100)
#'   fit$evectors[, 1:2]
#' }
#'
#' @export
#' @import Matrix
#' @import MASS
#' @import grpreg

ppwlssvm <- function(x, y, H=10, C=1, lambda, gamma=3.7, penalty="grSCAD", max.iter=100){

  n <- nrow(x)
  p <- ncol(x)

  #step1
  # x centering
  bar.x <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))

  #step2
  qprob <- seq(1/H, 1-1/H, by = 1/H)
  h <- length(qprob)
  tmp.y <- rep(y, times=h)
  pi.grid <- rep(qprob, each=n)


  # S matrix
  Sigma <- cov(x)  #time consuming
  Sigma.tilde <- cbind(rep(0, p+1), rbind(rep(0, p),Sigma))

  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- Sigma.tilde
  S <- Matrix::bdiag(temp)

  # weighted X.tilde
  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- x.tilde
  W <- X.tilde <- Matrix::bdiag(temp)


  weight <- (1-pi.grid)*as.numeric(tmp.y==1)+(pi.grid)*(as.numeric(tmp.y==-1)) #s,s
  sqrt_weight <- sqrt(weight)


  #multiply weight
  W <- X.tilde <- as.matrix(sqrt_weight*W)
  u.tilde <- sqrt(weight)*tmp.y

  #big X and Y
  a.eig <- eigen((n/C)*S+t(W)%*%W )  ##time!!##
  G.tilde.big <- (C)*a.eig$vectors %*% (sqrt(a.eig$values) * solve(a.eig$vectors))
  xi <- (Matrix::chol2inv(G.tilde.big)%*%t(W))%*%u.tilde
  xi.tilde <- xi <- (as.vector(xi))

  #restart the original process
  group <- colnames(G.tilde.big) <- rep(c(1:(p+1)), times=h)

  #initialization
  J <- max(as.integer(factor(group)))
  init.theta <- rnorm(h*(p+1), 0, 0.1)
  res <- xi.tilde - G.tilde.big %*% init.theta
  res_vec <- c()
  set.seed(1)
  obj_grpreg <- grpreg::grpreg(X=G.tilde.big, y=xi.tilde, group=group, penalty, family="gaussian",
                         lambda=lambda, alpha=1, eps=1e-5, max.iter=max.iter, dfmax=p,
                         gmax=length(unique(group)), gamma=gamma)


    #extract betas
    theta <- obj_grpreg$beta
    theta  <- theta[-1]

    #drop intercept
    intercept_ind <- which(((1:(h*(p+1))) %% (p+1) == 1))
    beta_mat <- matrix(theta, ncol = h)[-intercept_ind, ]

    #Working matrix (Mn)
    Mn <- matrix(0, p, p)
    for (h in 1:length(qprob)) Mn <- Mn + beta_mat[,h, drop = F] %*% t(beta_mat[,h, drop = F])
    rslt <- list("evalues" = eigen(Mn)$values, "evectors"=eigen(Mn)$vectors, "x"=x)
  return(rslt)
}


