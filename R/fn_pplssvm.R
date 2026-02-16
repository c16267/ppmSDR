#' Penalized Principal Least Squares Support Vector Machine (P\eqn{^2}LSSVM) for Sparse Sufficient Dimension Reduction
#'
#' This function implements the penalized principal least squares support vector machine (P\eqn{^2}LSSVM) approach
#' for sparse sufficient dimension reduction (SDR).
#'
#' @param x Numeric predictor matrix (n x p), where n is the sample size and p is the number of predictors.
#' @param y Numeric response vector of length n.
#' @param H The number of slicing point (default: 10).
#' @param C Regularization parameter (default: 1).
#' @param lambda Penalty parameter for sparsity.
#' @param gamma Regularization parameter for SCAD/MCP penalty (default: 3.7).
#' @param penalty Penalty type: `"grSCAD"` (default), `"grLasso"`, or `"grMCP"`.
#' @param max.iter Maximum number of iterations for the algorithm (default: 100).
#'
#' @details
#' This method estimates a sparse basis for the central subspace by minimizing a penalized principal least squares SVM objective.
#' It supports group penalties to induce row-wise sparsity and efficiently computes solutions even for high-dimensional data.
#'
#' @return
#' A list containing:
#'   \item{evalues}{Eigenvalues of the estimated central subspace matrix.}
#'   \item{evectors}{Eigenvectors (columns) corresponding to the estimated sufficient directions.}
#'   \item{x}{The original predictor matrix.}
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   n <- 100; p <- 10
#'   B <- matrix(0, p, 2)
#'   B[1,1] <- B[2,2] <- 1
#'   x <- MASS::mvrnorm(n, rep(0, p), diag(1, p))
#'   y <- (x %*% B[,1] / (0.5 + (x %*% B[,2] + 1)^2)) + 0.2 * rnorm(n)
#'   fit <- pplssvm(x, y, H = 10, C = 1, lambda = 0.03, gamma = 3.7,
#'                penalty = "grSCAD", max.iter = 100)
#'   fit$evectors[, 1:2]
#' }
#'
#' @export
#' @import Matrix
#' @importFrom stats cov glm model.matrix na.omit quantile relevel residuals rnorm
#' @importFrom Matrix bdiag

pplssvm <- function(x, y, H=10, C=1, lambda, gamma=3.7, penalty="grSCAD", max.iter=100)
{

  n <- nrow(x)
  p <- ncol(x)

  #step1
  #Y.tilde : dichotomize
  qprob <- seq(1/H, 1-1/H, by = 1/H)
  h <- length(qprob)
  qy <- quantile(y, qprob)
  U <- sapply(qy, function(x) 2*(x>y) - 1)
  y.tilde <- c(U)  # hxn

  #step2
  # x centering
  bar.x <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))


  #S : expanded covariance matrix # h(p+1) x h(p+1)
  Sigma.hat <- cov(x)
  Sigma.hat.star <- cbind(rep(0,p+1),rbind(rep(0,p), Sigma.hat))
  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- Sigma.hat.star
  S <- Matrix::bdiag(temp)  # h(p+1) x h(p+1)


  #W : expanded data matrix
  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- x.tilde
  W <- as.matrix(Matrix::bdiag(temp)) # nh x h(p+1)

  #G.tilde  #(p+1)h x (p+1)h
  a.eig <- eigen((n/C)*S+t(W)%*%W)
  G.tilde.big <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
  #G.tilde.big <- sqrtm( (n/C)*S+t(W)%*%W )

  #xi # (p+1)h
  #cholesk decompositon
  xi <- C*(Matrix::chol2inv(G.tilde.big)%*%t(W))%*%y.tilde
  #xi <- sqrtm(solve((n/C)*S+t(W)%*%W))%*%t(W)%*%y.tilde
  xi.tilde <- xi <- (as.vector(xi))

  #restart the original process


  #grouping
  group <- colnames(G.tilde.big) <- rep(c(1:(p+1)), times=h)


  #initialization
  J <- max(as.integer(factor(group)))
  init.theta <- rnorm(h*(p+1), 0, 1)
  res <- xi.tilde - G.tilde.big %*% init.theta
  res_vec <- c()

  if(penalty == "grSCAD"){
    for(iter in 1:max.iter){
      for(j in seq_along(integer(J))){
        ind <- which(colnames(G.tilde.big)==j)
        z_j <- (C/n)*(t(G.tilde.big[,ind,drop=F])%*%res + init.theta[ind]); #z_j <- as.vector(z_j)
        z_norm <- norm(z_j, "2")
        theta_j_new <- scad_firm_thresholding(z_norm, lambda1=lambda, lambda2=2, gamma)*z_j/z_norm
        resid_j_new <- res - (G.tilde.big[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        init.theta[ind] <- theta_j_new
        res <- c(resid_j_new)
      }
    }

  }else if(penalty == "grLasso")
    for(iter in 1:max.iter){
      for(j in seq_along(integer(J))){
        ind <- which(colnames(G.tilde.big)==j)
        z_j <- (C/n)*(t(G.tilde.big[,ind,drop=F])%*%res + init.theta[ind]); #z_j <- as.vector(z_j)
        z_norm <- norm(z_j, "2")
        theta_j_new <- soft_thresholding(z_norm, lambda=lambda)/(1 + 2)*z_j/z_norm
        resid_j_new <- res - (G.tilde.big[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        init.theta[ind] <- theta_j_new
        res <- c(resid_j_new)
      }
    }else if(penalty == "grMCP"){
      for(iter in 1:max.iter){
        for(j in seq_along(integer(J))){
          ind <- which(colnames(G.tilde.big)==j)
          z_j <- (C/n)*(t(G.tilde.big[,ind,drop=F])%*%res + init.theta[ind]); #z_j <- as.vector(z_j)
          z_norm <- norm(z_j, "2")
          theta_j_new <- firm_thresholding(z_norm, lambda1=lambda, lambda2=2, gamma)*z_j/z_norm
          resid_j_new <- res - (G.tilde.big[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
          init.theta[ind] <- theta_j_new
          res <- c(resid_j_new)
        }
      }
    }


  #Working matrix (Mn)
  theta  <- init.theta

  intercept_ind <- which(((1:(h*(p+1))) %% (p+1) == 1))
  beta_mat <- matrix(theta, ncol = h)[-intercept_ind, ]

  Mn <- matrix(0, p, p)
  for (h in 1:length(qprob)) Mn <- Mn + beta_mat[,h, drop = F] %*% t(beta_mat[,h, drop = F])

  rslt <- list("evalues" = eigen(Mn)$values, "evectors"=eigen(Mn)$vectors, "x"=x)
  return(rslt)
}

