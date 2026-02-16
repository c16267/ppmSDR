#' Penalized Principal Logistic Regression for Sparse Sufficient Dimension Reduction
#'
#' This function implements the penalized principal logistic regression (P\eqn{^2}LR) approach
#' for sparse sufficient dimension reduction (SDR), as proposed in Shin et al. (2024).
#'
#' @param x Numeric predictor matrix (n x p), where n is the sample size and p is the number of predictors.
#' @param y Numeric response vector of length n.
#' @param H The number of slicing quantiles (default: 10).
#' @param C Regularization parameter (default: 1).
#' @param lambda Penalty parameter for sparsity (default: 0.001).
#' @param gamma Regularization parameter for SCAD/MCP penalty (default: 3.7).
#' @param penalty Penalty type: `"grSCAD"` (default), `"grLasso"`, or `"grMCP"`.
#' @param max.iter Maximum number of iterations for the algorithm (default: 100).
#'
#' @details
#' This method estimates a sparse basis for the central subspace by minimizing a penalized principal logistic regression objective.
#' Supports group penalties to induce row-wise sparsity and efficiently computes solutions even for high-dimensional data.
#'
#' @return
#' A list containing:
#'   \item{Mn}{The estimated kernel matrix for the central subspace.}
#'   \item{evalues}{Eigenvalues of the estimated central subspace matrix.}
#'   \item{evectors}{Eigenvectors (columns) corresponding to the estimated sufficient directions.}
#'   \item{x}{The original predictor matrix.}
#'
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   n <- 100; p <- 10
#'   B <- matrix(0, p, 2)
#'   B[1,1] <- B[2,2] <- 1
#'   x <- MASS::mvrnorm(n, rep(0, p), diag(1, p))
#'   y <- (x %*% B[,1] / (0.5 + (x %*% B[,2] + 1)^2)) + 0.2 * rnorm(n)
#'   fit <- pplr(x, y, H = 10, C = 1, lambda = 0.00005, gamma = 3.7,
#'                penalty = "grSCAD", max.iter = 100)
#'   fit$evectors[, 1:2]
#' }
#'
#' @export
#' @import Matrix
#' @importFrom stats cov glm model.matrix na.omit quantile relevel residuals rnorm
#' @importFrom Matrix bdiag
#'
pplr <- function(x, y, H = 10, C = 1,  lambda, gamma=3.7, penalty = "grSCAD", max.iter = 100){
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

  init.theta <- theta0 <- rep(0.2, (h*(p+1)))
  tol <-  10^-5
  if(penalty=="grSCAD"){
    for(iter in 1:max.iter){
      # initialization
      pi <- pi0 <- as.vector(1/(1 + exp(c(u) * (X.tilde %*% init.theta))))
      tmp_pi <- pi * (1- pi);
      sqrt_tmp_pi <- sqrt(tmp_pi)
      G <- G0 <- Matrix::Diagonal(x=tmp_pi)
      u.tilde <- u0.tilde <- as.vector(W %*% init.theta + (1/tmp_pi * (u * pi)))

      #multiply G:Hessian matrix
      u.tilde <- as.vector(sqrt_tmp_pi*u.tilde)
      W <- X.tilde <- as.matrix(sqrt_tmp_pi*W)  #time!!##

      ####how to###
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
        z_j <- (C/(n))*(t(big_X[,ind,drop=F])%*%res)+ init.theta[ind];
        z_norm <- norm(z_j, "2")
        theta_j_new <- scad_firm_thresholding(z_norm, lambda1=lambda, lambda2=2, gamma=gamma)*z_j/z_norm
        resid_j_new <- res - (big_X[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        init.theta[ind] <- theta_j_new
        res <- c(resid_j_new)
      }

      res.new <- res
      new.theta <- init.theta
      #print(round(res,5))

      delta <- mean(abs(abs(res.new)-abs(res.old)))
      #print(paste("iteration :", iter))
      if (delta < tol) break
    }
  }else if(penalty == "grLasso"){
      for(iter in 1:max.iter){
        # initialization
        pi <- pi0 <- as.vector(1/(1 + exp(c(u) * (X.tilde %*% init.theta))))
        tmp_pi <- pi * (1- pi);
        sqrt_tmp_pi <- sqrt(tmp_pi)
        G <- G0 <- Matrix::Diagonal(x=tmp_pi)
        u.tilde <- u0.tilde <- as.vector(W %*% init.theta + (1/tmp_pi * (u * pi)))

        #multiply G:Hessian matrix
        u.tilde <- as.vector(sqrt_tmp_pi*u.tilde)
        W <- X.tilde <- as.matrix(sqrt_tmp_pi*W)  #time!!##

        ####how to###
        a.eig <- eigen(S+t(W)%*%W )  ##time!!##
        big_X <- a.eig$vectors %*% (sqrt(a.eig$values) * solve(a.eig$vectors))


        #cholesky decomposition
        big_Y <- as.vector(Matrix::chol2inv(big_X)%*%(t(W)%*%u.tilde))

        group <- colnames(big_X) <- rep(c(1:(p+1)), times=h)
        J <- max(as.integer(factor(group)))

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
        print(round(res,5))

        delta <- mean(abs(abs(res.new)-abs(res.old)))
        print(paste("iteration :", iter))
        if (delta < tol) break
      }
    }else if(penalty == "grMCP"){
          for(iter in 1:max.iter){
            # initialization
            pi <- pi0 <- as.vector(1/(1 + exp(c(u) * (X.tilde %*% init.theta))))
            tmp_pi <- pi * (1- pi);
            sqrt_tmp_pi <- sqrt(tmp_pi)
            G <- G0 <- Matrix::Diagonal(x=tmp_pi)
            u.tilde <- u0.tilde <- as.vector(W %*% init.theta + (1/tmp_pi * (u * pi)))

            #multiply G:Hessian matrix
            u.tilde <- as.vector(sqrt_tmp_pi*u.tilde)
            W <- X.tilde <- as.matrix(sqrt_tmp_pi*W)  #time!!##

            ####how to###
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
              z_j <- (C/(n))*(t(big_X[,ind,drop=F])%*%res)+ init.theta[ind];
              z_norm <- norm(z_j, "2")
              theta_j_new <- firm_thresholding(z_norm, lambda1=lambda, lambda2=2, gamma=gamma)*z_j/z_norm
              resid_j_new <- res - (big_X[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
              init.theta[ind] <- theta_j_new
              res <- c(resid_j_new)
            }

            res.new <- res
            new.theta <- init.theta
            print(round(res,5))

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


