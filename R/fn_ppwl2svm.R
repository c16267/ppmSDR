#' Penalized Principal Weighted L2-Hinge SVM (P\eqn{^2}WL2M) for Sparse Sufficient Dimension Reduction
#'
#' This function implements the penalized principal weighted L2-hinge support vector machine (P\eqn{^2}WL2M) method for sparse sufficient dimension reduction (SDR), typically for binary responses. Group penalties (SCAD, Lasso, MCP) are supported for variable selection.
#'
#' @param x Numeric predictor matrix (n x p), where n is the sample size and p is the number of predictors.
#' @param y Binary response vector of length n, coded as \code{-1, 1}.
#' @param H Number of quantile slices for principal machine construction (default: 10).
#' @param C Regularization parameter (default: 1).
#' @param lambda Penalty parameter for sparsity.
#' @param gamma Regularization parameter for SCAD/MCP penalty (default: 3.7).
#' @param penalty Penalty type: \code{"grSCAD"} (default), \code{"grLasso"}, or \code{"grMCP"}.
#' @param max.iter Maximum number of iterations (default: 100).
#'
#' @details
#' This function solves a group-penalized principal machine with weighted L2 hinge loss, suitable for classification tasks. Efficient group coordinate descent (GCD) is used for parameter estimation. Use only with binary responses (\code{-1, 1}).
#'
#' @return
#' A list with components:
#'   \item{evalues}{Eigenvalues of the estimated subspace matrix.}
#'   \item{evectors}{Eigenvectors spanning the sufficient dimension reduction subspace.}
#'   \item{x}{The original predictor matrix.}
#'
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   n <- 1000; p <- 10
#'   x <- MASS::mvrnorm(n, rep(0, p), diag(1, p))
#'   y <- (x %*% B[,1] / (0.5 + (x %*% B[,2] + 1)^2)) + 0.2 * rnorm(n)
#'   y.binary <- sign(y)
#'   fit <- ppwl2svm(x, y.binary, H = 10, C = 5, lambda = 0.08, gamma = 3.7,
#'             penalty = "grSCAD", max.iter = 100)
#'   fit$evectors[, 1:2]
#' }
#'
#' @export
#' @import MASS
#' @import Matrix


ppwl2svm <- function(x, y, H = 10, C = 1, lambda, gamma = 3.7, penalty = "grSCAD", max.iter=100)
{

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
  S <- bdiag(temp)

  # weighted X.tilde
  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- x.tilde
  W <- X.tilde <- Matrix::bdiag(temp)

  # applying square root of the class weight onto both X and Y
  weight <- (1-pi.grid)*as.numeric(tmp.y==1)+(pi.grid)*(as.numeric(tmp.y==-1)) #s,s
  sqrt_weight <- sqrt(weight)

  W <- as.matrix(sqrt_weight*W)
  tmp.y <- sqrt(weight)*tmp.y

  #initialization
  set.seed(1)
  init.theta <- rnorm(h*(p+1), 0, 0.2)
  tau.vec <- rep(0, each=n* h)


  if(penalty == "grSCAD"){

    for(iter in 1:max.iter){

      u <- (1-tau.vec)*as.matrix(as.integer(tmp.y - as.matrix(W)%*%init.theta <= 0), ncol=length(y)) #I_{-} at t-th iteration
      u.c <- (tau.vec)*as.matrix(as.integer(tmp.y - as.matrix(W)%*%init.theta > 0), ncol=length(y)) #I_{+} at t-th iteration

      y.tilde <- tmp.y*sqrt(u)
      y.tilde.c <- tmp.y*sqrt(u.c)

      u.mat <- matrix(rep(sqrt(u),ncol(W)),ncol=ncol(W))
      u.c.mat <- matrix(rep(sqrt(u.c),ncol(W)),ncol=ncol(W))

      w.tilde <- u.mat*W
      w.tilde.c <-u.c.mat*W

      a.eig <- eigen(0.5*S+t(w.tilde)%*%w.tilde)
      G.tilde.big <- (a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)) #big X^{t}

      xi <- (solve(G.tilde.big)%*%t(w.tilde)%*%y.tilde) #big Y^{t}
      xi.tilde <- xi <- (as.vector(xi));

      res <- xi.tilde - G.tilde.big %*% init.theta
      n.new <- length(xi.tilde)

      #grouping
      group <- colnames(G.tilde.big) <- rep(c(1:(p+1)), times=h)
      J <- max(as.integer(factor(group)))
      old.theta <- init.theta

      for(j in seq_along(integer(J))){
        ind <- which(colnames(G.tilde.big)==j)
        z_j <- (C/(n))*((t(G.tilde.big[,ind,drop=F])%*%res)) + init.theta[ind]; #z_j <- as.vector(z_j)
        z_norm <- norm(z_j, "2")
        theta_j_new <- scad_firm_thresholding(z_norm, lambda1=lambda, lambda2=1, gamma=gamma)*z_j/z_norm

        resid_j_new_pos <- res - (G.tilde.big[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        init.theta[ind] <- theta_j_new
        res <- c(resid_j_new_pos)
      }
      new.theta <- init.theta
      delta <- max(na.omit(abs(new.theta - old.theta)/abs(old.theta)))
      if (delta < 1.0e-5) break

    }

  }else if(penalty == "grLasso"){
    for(iter in 1:max.iter){

      u <- (1-tau.vec)*as.matrix(as.integer(tmp.y - as.matrix(W)%*%init.theta <= 0), ncol=length(y))
      u.c <- (tau.vec)*as.matrix(as.integer(tmp.y - as.matrix(W)%*%init.theta > 0), ncol=length(y))

      y.tilde <- tmp.y*sqrt(u)
      y.tilde.c <- tmp.y*sqrt(u.c)

      u.mat <- matrix(rep(sqrt(u),ncol(W)),ncol=ncol(W))
      u.c.mat <- matrix(rep(sqrt(u.c),ncol(W)),ncol=ncol(W))

      w.tilde <- u.mat*W

      a.eig <- eigen(0.5*S+t(w.tilde)%*%w.tilde)
      G.tilde.big <- (a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors))

      xi <- (solve(G.tilde.big)%*%t(w.tilde)%*%y.tilde)
      xi.tilde <- xi <- (as.vector(xi));

      res <- xi.tilde - G.tilde.big %*% init.theta
      n.new <- length(xi.tilde)

      #grouping
      group <- colnames(G.tilde.big) <- rep(c(1:(p+1)), times=h)
      J <- max(as.integer(factor(group)))
      old.theta <- init.theta

      for(j in seq_along(integer(J))){
        ind <- which(colnames(G.tilde.big)==j)
        z_j <- (C/(n))*((t(G.tilde.big[,ind,drop=F])%*%res)) + init.theta[ind]; #z_j <- as.vector(z_j)
        z_norm <- norm(z_j, "2")
        theta_j_new <- soft_thresholding(z_norm, lambda=lambda)/(1 + 2)*z_j/z_norm
        resid_j_new_pos <- res - (G.tilde.big[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        init.theta[ind] <- theta_j_new
        res <- c(resid_j_new_pos)
      }
      new.theta <- init.theta
      delta <- max(na.omit(abs(new.theta - old.theta)/abs(old.theta)))
      if (delta < 1.0e-5) break
    }
  }else if(penalty == "grMCP"){
    for(iter in 1:max.iter){

      u <- (1-tau.vec)*as.matrix(as.integer(tmp.y - as.matrix(W)%*%init.theta <= 0), ncol=length(y))
      u.c <- (tau.vec)*as.matrix(as.integer(tmp.y - as.matrix(W)%*%init.theta > 0), ncol=length(y))

      y.tilde <- tmp.y*sqrt(u)
      y.tilde.c <- tmp.y*sqrt(u.c)

      u.mat <- matrix(rep(sqrt(u),ncol(W)),ncol=ncol(W))
      u.c.mat <- matrix(rep(sqrt(u.c),ncol(W)),ncol=ncol(W))

      w.tilde <- u.mat*W
      w.tilde.c <-u.c.mat*W

      a.eig <- eigen(0.5*S+t(w.tilde)%*%w.tilde)
      G.tilde.big <- (a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors))

      xi <- (solve(G.tilde.big)%*%t(w.tilde)%*%y.tilde)
      xi.tilde <- xi <- (as.vector(xi));

      res <- xi.tilde - G.tilde.big %*% init.theta
      n.new <- length(xi.tilde)

      #grouping
      group <- colnames(G.tilde.big) <- rep(c(1:(p+1)), times=h)
      J <- max(as.integer(factor(group)))
      old.theta <- init.theta

      for(j in seq_along(integer(J))){
        ind <- which(colnames(G.tilde.big)==j)
        z_j <- (C/(n))*((t(G.tilde.big[,ind,drop=F])%*%res)) + init.theta[ind]; #z_j <- as.vector(z_j)
        z_norm <- norm(z_j, "2")
        theta_j_new <- firm_thresholding(z_norm, lambda1=lambda, lambda2=1, gamma=gamma)*z_j/z_norm

        resid_j_new_pos <- res - (G.tilde.big[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        init.theta[ind] <- theta_j_new
        res <- c(resid_j_new_pos)
      }
      new.theta <- init.theta
      delta <- max(na.omit(abs(new.theta - old.theta)/abs(old.theta)))
      if (delta < 1.0e-5) break
    }
  }

  #Working matrix (Mn)
  theta  <- init.theta

  intercept_ind <- which(((1:(h*(p+1))) %% (p+1) == 1))
  beta_mat <- matrix(theta, ncol = h)[-intercept_ind, ]

  Mn <- matrix(0, p, p)
  for (h in 1:length(qprob)) Mn <- Mn + beta_mat[,h, drop = F] %*% t(beta_mat[,h, drop = F])

  rslt <- list("evalues" = eigen(Mn)$values, "evectors"=eigen(Mn)$vectors, 'x'=x)
  return(rslt)
}



