#' Penalized Principal Support Vector Machine (P\eqn{^2}SVM, MM-GCD-based) for Sparse Sufficient Dimension Reduction
#'
#' This function implements the penalized principal support vector machine (P\eqn{^2}SVM) approach for sparse sufficient dimension reduction (SDR), using the MM-GCD algorithm for efficient computation.
#'
#' @param x Numeric predictor matrix (n x p), where n is the sample size and p is the number of predictors.
#' @param y Numeric response vector of length n.
#' @param H The number of quantile slices (default: 10).
#' @param C Regularization parameter (default: 1).
#' @param lambda Penalty parameter for sparsity.
#' @param gamma Regularization parameter for SCAD/MCP penalty (default: 3.7).
#' @param penalty Penalty type: \code{"grSCAD"} (default), \code{"grLasso"}, or \code{"grMCP"}.
#' @param max.iter Maximum number of iterations for the algorithm (default: 100).
#'
#' @details
#' This function estimates a sparse basis for the central subspace by solving a penalized principal support vector machine objective with an MM-GCD optimization algorithm, enabling scalable and accurate sparse SDR.
#'
#' @return
#' A list containing:
#'   \item{Mn}{Estimated central subspace matrix.}
#'   \item{evalues}{Eigenvalues of the estimated matrix.}
#'   \item{evectors}{Eigenvectors (columns) corresponding to the estimated sufficient directions.}
#'   \item{x}{The original predictor matrix.}
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   n <- 300; p <- 10
#'   B <- matrix(0, p, 2)
#'   B[1,1] <- B[2,2] <- 1
#'   x <- MASS::mvrnorm(n, rep(0, p), diag(1, p))
#'   y <- (x %*% B[,1] / (0.5 + (x %*% B[,2] + 1)^2)) + 0.2 * rnorm(n)
#'   fit <- ppsvm(x, y, H = 10, C = 100, lambda = 0.0001, gamma = 3.7,
#'          penalty = "grSCAD", max.iter = 100)
#'   fit$evectors[, 1:2]
#' }
#'
#' @export
#' @import MASS
#' @import Matrix


ppsvm <- function(x, y, H=10, C=1,  lambda, gamma=3.7, penalty="grSCAD", max.iter=100)
{
  n <- nrow(x)
  p <- ncol(x)

  # centering
  bar.x <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))

  # u vector Y nhx1
  qprob <- seq(1/H, 1-1/H, by = 1/H)
  h <- length(qprob)
  qy <- quantile(y, qprob)
  U <- sapply(qy, function(x) 2*(x>y) - 1)
  u <- c(U)

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

  init.theta <- theta0 <- rep(0.1, (h*(p+1)))
  tol = 10^-5

  if(penalty=="grSCAD"){
    for(iter in 1:max.iter){
      m_t <- u*(W%*%init.theta)/n;
      c_t <- as.vector(abs(1-m_t))
      omega_t <- 1/(4*c_t)
      sqrt_omega_t <- sqrt(omega_t)
      u.tilde <- (1+c_t)*u

      #multiply omega:Hessian matrix
      u.tilde <- as.vector(sqrt_omega_t*u.tilde)
      W <- X.tilde <- as.matrix(sqrt_omega_t*W)  #time!!##

      #a.eig <- eigen( (ncol(W)/C)*S+t(W)%*%W )  ##time!!##
      a.eig <- eigen(S+t(W)%*%W )  ##time!!##

      big_X <- a.eig$vectors %*% (sqrt(a.eig$values) * solve(a.eig$vectors))


      #cholesky decompositon
      big_Y <- as.vector(Matrix::chol2inv(big_X)%*%(t(W)%*%u.tilde))

      group <- colnames(big_X) <- rep(c(1:(p+1)), times=h)
      J <- max(as.integer(factor(group)))

      res <- big_Y - big_X %*% init.theta
      old.theta <- init.theta
      res.old <- res

      for(j in 1:J){
        ind <- which(colnames(big_X)==j)
        z_j <- (C/n)*(t(big_X[,ind,drop=F])%*%res) + init.theta[ind];
        z_norm <- norm(z_j, "2")
        theta_j_new <- scad_firm_thresholding(z_norm, lambda1=lambda, lambda2=2, gamma=gamma)*(z_j/z_norm)
        resid_j_new <- res - (big_X[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        init.theta[ind] <- theta_j_new
        res <- c(resid_j_new)
        #print(paste("z_norm:", z_norm))
      }

      res.new <- res
      new.theta <- init.theta

      delta <- mean(abs(abs(res.new)-abs(res.old)))
      #delta <- mean( abs((new.theta)-(old.theta))/abs(old.theta) )
      #print(paste("iteration :", iter, "delta:", delta))
      if (delta < tol) break
    }
  }else if(penalty == "grLasso"){
    for(iter in 1:max.iter){
      #iterative weight c_0 for response
      m_t <- u*(W%*%init.theta)/n;
      c_t <- as.vector(abs(1-m_t))
      omega_t <- 1/(4*c_t)
      sqrt_omega_t <- sqrt(omega_t)
      u.tilde <- (1+c_t)*u

      #multiply omega:Hessian matrix
      u.tilde <- as.vector(sqrt_omega_t*u.tilde)
      W <- X.tilde <- as.matrix(sqrt_omega_t*W)  #time!!##

      a.eig <- eigen(S+t(W)%*%W )  ##time!!##
      big_X <- a.eig$vectors %*% (sqrt(a.eig$values) * solve(a.eig$vectors))


      #cholesky decompositon
      big_Y <- as.vector(Matrix::chol2inv(big_X)%*%(t(W)%*%u.tilde))

      group <- colnames(big_X) <- rep(c(1:(p+1)), times=h)
      J <- max(as.integer(factor(group)))

      res <- big_Y - big_X %*% init.theta
      old.theta <- init.theta
      res.old <- res

      for(j in 1:J){
        ind <- which(colnames(big_X)==j)
        z_j <- (C/n)*(t(big_X[,ind,drop=F])%*%res) + init.theta[ind];
        z_norm <- norm(z_j, "2")
        theta_j_new <- soft_thresholding(z_norm, lambda=lambda)/(1 + 2)*z_j/z_norm
        resid_j_new <- res - (big_X[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        init.theta[ind] <- theta_j_new
        res <- c(resid_j_new)
      }

      res.new <- res
      new.theta <- init.theta

      delta <- mean(abs(abs(res.new)-abs(res.old)))
      #print(paste("iteration :", iter, "delta:", delta))
      if (delta < tol) break
    }
  }else if(penalty == "grMCP"){
    for(iter in 1:max.iter){
      #iterative weight c_0 for response
      m_t <- u*(W%*%init.theta)/n;
      c_t <- as.vector(abs(1-m_t))
      omega_t <- 1/(4*c_t)
      sqrt_omega_t <- sqrt(omega_t)
      u.tilde <- (1+c_t)*u

      #multiply omega:Hessian matrix
      u.tilde <- as.vector(sqrt_omega_t*u.tilde)
      W <- X.tilde <- as.matrix(sqrt_omega_t*W)  #time!!##

      a.eig <- eigen((n/C)*S+t(W)%*%W )  ##time!!##
      big_X <- a.eig$vectors %*% (sqrt(a.eig$values) * solve(a.eig$vectors))


      #cholesky decompositon
      big_Y <- as.vector(Matrix::chol2inv(big_X)%*%(t(W)%*%u.tilde))

      group <- colnames(big_X) <- rep(c(1:(p+1)), times=h)
      J <- max(as.integer(factor(group)))

      res <- big_Y - big_X %*% init.theta
      old.theta <- init.theta
      res.old <- res

      for(j in 1:J){
        ind <- which(colnames(big_X)==j)
        z_j <- (C/n)*(t(big_X[,ind,drop=F])%*%res) + init.theta[ind];
        z_norm <- norm(z_j, "2")
        theta_j_new <- firm_thresholding(z_norm, lambda1=lambda, lambda2=2, gamma=gamma)*z_j/z_norm
        resid_j_new <- res - (big_X[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        init.theta[ind] <- theta_j_new
        res <- c(resid_j_new)
      }

      res.new <- res
      new.theta <- init.theta

      delta <- mean(abs(abs(res.new)-abs(res.old)))
      #print(paste("iteration :", iter, "delta:", delta))
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



