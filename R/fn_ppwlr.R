#' Penalized Principal Weighted Logistic Regression (P\eqn{^2}WLR) for Sparse Sufficient Dimension Reduction
#'
#' Implements the penalized principal weighted logistic regression (P\eqn{^2}WLR) method for sparse sufficient dimension reduction (SDR), typically for binary responses.
#'
#' @param x Numeric predictor matrix (n x p), where n is the sample size and p is the number of predictors.
#' @param y Binary response vector of length n, coded as \code{-1, 1}.
#' @param H Number of quantile slices for principal machine construction.
#' @param C Regularization parameter.
#' @param lambda Penalty parameter for sparsity.
#' @param gamma Regularization parameter for SCAD/MCP penalty.
#' @param penalty Penalty type: \code{"grSCAD"}, \code{"grLasso"}, or \code{"grMCP"}.
#' @param max.iter Maximum number of iterations (default: 100).
#'
#' @details
#' This function fits a penalized principal machine using weighted logistic regression loss, which is well-suited for a binary classification in sufficient dimension reduction. Group penalties (SCAD, Lasso, MCP) are available for structured sparsity and variable selection. Only use with binary responses (\code{-1, 1}).
#'
#' @return
#' A list with components:
#'   \item{Mn}{Estimated sufficient dimension reduction matrix.}
#'   \item{evalues}{Eigenvalues of \code{Mn}.}
#'   \item{evectors}{Eigenvectors spanning the SDR subspace.}
#'   \item{x}{Original predictor matrix.}
#'
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   n <- 1000; p <- 10
#'   x <- MASS::mvrnorm(n, rep(0, p), diag(1, p))
#'   y <- (x %*% B[,1] / (0.5 + (x %*% B[,2] + 1)^2)) + 0.2 * rnorm(n)
#'   y.binary <- sign(y)
#'   fit <- ppwlr(x, y.binary, H = 10, C = 10, lambda = 0.02, gamma = 3.7,
#'                penalty = "grSCAD", max.iter = 100)
#'   fit$evectors[, 1:2]
#' }
#'
#' @export
#' @import MASS
#' @import Matrix


ppwlr <- function(x, y, H = 10, C = 1, lambda = 0.001, gamma = 3.7, penalty = 'grSCAD', max.iter=100){
  n <- nrow(x)
  p <- ncol(x)

  # centering
  bar.x <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))

  # u vector Y nhx1
  qprob <- seq(1/H, 1-1/H, by = 1/H)
  step <- 1/H; h <- length(qprob)
  pi.grid <- seq(step, 1-step, by = step)
  u <- rep(y, h)

  #Omega matrix
  Omega_mat <- matrix( ((u+1)/2)*rep(pi.grid, each=length(y)) + ((-u+1)/2)*rep(1-pi.grid, each=length(y)), ncol=h)
  Omega_vec <- as.vector(Omega_mat)

  # S matrix
  Sigma <- cov(x)  #time consuming
  Sigma.tilde <- cbind(rep(0, p+1), rbind(rep(0, p),Sigma))

  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- Sigma.tilde
  S <- Matrix::bdiag(temp)

  # X.tilde
  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- x.tilde
  W <- X.tilde <- Matrix::bdiag(temp)

  init.theta <- theta0 <- rnorm((h*(p+1)), 0, 1)
  tol <- 10^-5

  if(penalty=="grSCAD"){
    for(iter in 1:max.iter){
      # initialization
      pi <- pi0 <- as.vector(1/(1 + exp(c(u) * (X.tilde %*% init.theta))))
      tmp_pi <- pi * (1- pi)*Omega_vec;
      sqrt_tmp_pi <- sqrt(tmp_pi)
      G <- G0 <- Matrix::Diagonal(x=tmp_pi)
      #u.tilde <- u0.tilde <- as.vector(W %*% init.theta + (1/tmp_pi * (u * pi)*Omega_vec))
      u.tilde <- u0.tilde <- as.vector(W %*% init.theta*Omega_vec + (1/tmp_pi * (u * pi)))

      #multiply G:Hessian matrix
      u.tilde_new <- as.vector(sqrt_tmp_pi*u.tilde)
      W_new <- as.matrix(sqrt_tmp_pi*W)

      a.eig <- eigen(S+t(W_new)%*%W_new )  ##time!!##
      big_X <- (a.eig$vectors %*% (sqrt(a.eig$values) * solve(a.eig$vectors)))


      #cholesky decompositon
      big_Y <- as.vector(Matrix::chol2inv(big_X)%*%(t(W_new)%*%u.tilde_new))

      group <- colnames(big_X) <- rep(c(1:(p+1)), times=h)
      J <- max(as.integer(factor(group)))

      res <- big_Y - big_X %*% init.theta
      old.theta <- init.theta
      res.old <- res
      for(j in 1:J){
        ind <- which(colnames(big_X)==j)
        z_j <- (C/(n))*(t(big_X[,ind,drop=F])%*%res)+ init.theta[ind];
        z_norm <- norm(z_j, "2")
        theta_j_new <- scad_firm_thresholding(z_norm, lambda1=lambda, lambda2=1, gamma=gamma)*z_j/z_norm
        resid_j_new <- res - (big_X[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        init.theta[ind] <- theta_j_new
        res <- c(resid_j_new)
      }
      res.new <- res
      new.theta <- init.theta
      delta <- mean(abs(abs(res.new)-abs(res.old)))
      #print(paste("iteration :", iter))
      if (delta < tol) break
    }
  }else if(penalty == "grLasso"){
    for(iter in 1:max.iter){
      # initialization
      pi <- pi0 <- as.vector(1/(1 + exp(c(u) * (X.tilde %*% init.theta))))
      tmp_pi <- pi * (1- pi)*Omega_vec;
      sqrt_tmp_pi <- sqrt(tmp_pi)
      G <- G0 <- Matrix::Diagonal(x=tmp_pi)
      u.tilde <- u0.tilde <- as.vector(W %*% init.theta + (1/tmp_pi * (u * pi)*Omega_vec))

      #multiply G:Hessian matrix
      u.tilde_new <- as.vector(sqrt_tmp_pi*u.tilde)
      W_new <- as.matrix(sqrt_tmp_pi*W)

      ####how to###
      a.eig <- eigen(S+t(W_new)%*%W_new )  ##time!!##
      big_X <- a.eig$vectors %*% (sqrt(a.eig$values) * solve(a.eig$vectors))


      #cholesky decompositon
      big_Y <- as.vector(Matrix::chol2inv(big_X)%*%(t(W_new)%*%u.tilde_new))

      group <- colnames(big_X) <- rep(c(1:(p+1)), times=h)
      J <- max(as.integer(factor(group)))

      ind_intercept <- which(colnames(big_X)==1)

      res <- big_Y - big_X %*% init.theta
      old.theta <- init.theta
      res.old <- res
      for(j in 1:J){
        ind <- which(colnames(big_X)==j)
        z_j <- (C/(n))*(t(big_X[,ind,drop=F])%*%res)+ init.theta[ind];
        z_norm <- norm(z_j, "2")
        theta_j_new <- soft_thresholding(z_norm, lambda=lambda)/(1 + 2)*z_j/z_norm
        resid_j_new <- res - (big_X[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        init.theta[ind] <- theta_j_new
        res <- c(resid_j_new)
      }
      res.new <- res
      new.theta <- init.theta
      delta <- mean(abs(abs(res.new)-abs(res.old)))
      print(paste("iteration :", iter))
      if (delta < tol) break
    }
  }else if(penalty == "grMCP"){
    for(iter in 1:max.iter){
      # initialization
      pi <- pi0 <- as.vector(1/(1 + exp(c(u) * (X.tilde %*% init.theta))))
      tmp_pi <- pi * (1- pi)*Omega_vec;
      sqrt_tmp_pi <- sqrt(tmp_pi)
      G <- G0 <- Matrix::Diagonal(x=tmp_pi)
      u.tilde <- u0.tilde <- as.vector(W %*% init.theta + (1/tmp_pi * (u * pi)*Omega_vec))

      #multiply G:Hessian matrix
      u.tilde_new <- as.vector(sqrt_tmp_pi*u.tilde)
      W_new <- as.matrix(sqrt_tmp_pi*W)

      ####how to###
      a.eig <- eigen(S+t(W_new)%*%W_new )  ##time!!##
      big_X <- a.eig$vectors %*% (sqrt(a.eig$values) * solve(a.eig$vectors))


      #cholesky decompositon
      big_Y <- as.vector(Matrix::chol2inv(big_X)%*%(t(W_new)%*%u.tilde_new))

      group <- colnames(big_X) <- rep(c(1:(p+1)), times=h)
      J <- max(as.integer(factor(group)))

      ind_intercept <- which(colnames(big_X)==1)

      res <- big_Y - big_X %*% init.theta
      old.theta <- init.theta
      res.old <- res
      for(j in 1:J){
        ind <- which(colnames(big_X)==j)
        z_j <- (C/(n))*(t(big_X[,ind,drop=F])%*%res)+ init.theta[ind];
        z_norm <- norm(z_j, "2")
        theta_j_new <- firm_thresholding(z_norm, lambda1=lambda, lambda2=1, gamma=gamma)*z_j/z_norm
        resid_j_new <- res - (big_X[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        init.theta[ind] <- theta_j_new
        res <- c(resid_j_new)
      }
      res.new <- res
      new.theta <- init.theta
      delta <- mean(abs(abs(res.new)-abs(res.old)))
      print(paste("iteration :", iter))
      if (delta < tol) break
    }
  }
  #Working matrix (Mn)
  theta  <- init.theta

  intercept_ind <- which(((1:(h*(p+1))) %% (p+1) == 1))
  beta_mat <- matrix(theta, ncol = h)[-intercept_ind, ]

  Mn <- matrix(0, p, p)
  for (h in 1:length(qprob)) Mn <- Mn + beta_mat[,h, drop = F] %*% t(beta_mat[,h, drop = F])

  rslt <- list("Mn"=Mn, "evalues" = eigen(Mn)$values, "evectors"=eigen(Mn)$vectors, 'x'=x)
  return(rslt)
}

