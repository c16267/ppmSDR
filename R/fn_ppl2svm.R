#' Penalized Principal L2-Hinge Support Vector Machine for Sparse Sufficient Dimension Reduction
#'
#' This function implements the penalized principal L2-hinge support vector machine (P\eqn{^2}L2M) approach
#' for sparse sufficient dimension reduction (SDR)
#'
#' @param x A numeric matrix of predictors (n x p), where n is the number of samples and p is the number of predictors.
#' @param y A numeric response vector of length n.
#' @param H The number of slicing or quantiles (default: 10).
#' @param C Regularization parameter (default: 1).
#' @param lambda Penalty parameter controlling sparsity (default: 0.001).
#' @param gamma Regularization parameter for SCAD/MCP penalty (default: 3.7).
#' @param penalty Type of penalty: `"grSCAD"` (default), `"grLasso"`, or `"grMCP"`.
#' @param max.iter Maximum number of iterations for the algorithm (default: 100).
#'
#' @details
#' The method estimates the sparse basis of the central subspace by minimizing a penalized principal L2-hinge SVM objective.
#' The algorithm supports group penalties to encourage row-wise sparsity, and efficiently computes the solution even in high-dimensional settings.
#'
#' @return
#' A list with the following components:
#'   \item{evalues}{Eigenvalues of the estimated matrix spanning the central subspace.}
#'   \item{evectors}{Corresponding eigenvectors (columns) for the estimated sufficient directions.}
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
#'   fit <- ppl2svm(x, y, H = 10, C = 1, lambda = 0.2, gamma = 3.7,
#'                   penalty = "grSCAD", max.iter = 100)
#'   fit$evectors[, 1:2]
#' }
#'
#' @export
#' @import Matrix
#' @importFrom stats cov glm model.matrix na.omit quantile relevel residuals rnorm
#' @importFrom Matrix bdiag

ppl2svm <- function(x, y, H = 10, C = 1, lambda = 0.001, gamma = 3.7, penalty = "grSCAD", max.iter = 100)
{


  n <- nrow(x)
  p <- ncol(x)

  #step1 :Build a Design matrix
  bar.x <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))


  #step2
  #Y.tilde : dichotomize
  qprob <- seq(1/H, 1-1/H, by = 1/H)
  h <- length(qprob)
  qy <- quantile(y, qprob)
  U <- sapply(qy, function(x) 2*(x>y) - 1)
  tmp.y <- c(U)  # hxn

  #S : expanded covariance matrix # h(p+1) x h(p+1)
  Sigma.hat <- cov(x)
  Sigma.hat.star <- cbind(rep(0,p+1),rbind(rep(0,p), Sigma.hat))
  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- Sigma.hat.star
  S <- as.matrix(Matrix::bdiag(temp))  # h(p+1) x h(p+1)

  #W : expanded data matrix
  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- x.tilde
  W <- Matrix::bdiag(temp) # nh x h(p+1)
  W <- as.matrix(W)

  #initialization
  init.theta <- rnorm(h*(p+1), 0, 1)
  tau.vec <- rep(0, each=n* h)


  if(penalty == "grSCAD"){

    for(iter in 1:max.iter){

      u <- (1-tau.vec)*as.matrix(as.integer(tmp.y - as.matrix(W)%*%init.theta <= 0), ncol=length(y)) #I_{-}
      u.c <- (tau.vec)*as.matrix(as.integer(tmp.y - as.matrix(W)%*%init.theta > 0), ncol=length(y)) #I_{+}

      y.tilde <- tmp.y*sqrt(u)
      y.tilde.c <- tmp.y*sqrt(u.c)

      u.mat <- matrix(rep(sqrt(u),ncol(W)),ncol=ncol(W))
      u.c.mat <- matrix(rep(sqrt(u.c),ncol(W)),ncol=ncol(W))

      w.tilde <- u.mat*W
      w.tilde.c <-u.c.mat*W

      a.eig <- eigen(0.5*S+t(w.tilde)%*%w.tilde)
      G.tilde.big <- (a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors))

      b.eig <- eigen(0.5*S+t(w.tilde.c)%*%w.tilde.c)
      G.tilde.big.c <- (b.eig$vectors %*% diag(sqrt(b.eig$values)) %*% solve(b.eig$vectors))

      xi <- (solve(G.tilde.big)%*%t(w.tilde)%*%y.tilde)
      #xi.c <-(solve(G.tilde.big.c)%*%t(w.tilde.c)%*%y.tilde.c)

      xi.tilde <- xi <- (as.vector(xi));
      #xi.tilde.c <- xi.c <- (as.vector(xi.c))

      res <- xi.tilde - G.tilde.big %*% init.theta
      #neg.resid <- xi.tilde.c - G.tilde.big.c %*% init.theta
      #res <- pos.resid + neg.resid

      n.new <- length(xi.tilde)

      #grouping
      group <- colnames(G.tilde.big) <- rep(c(1:(p+1)), times=h)
      J <- max(as.integer(factor(group)))
      old.theta <- init.theta

      for(j in seq_along(integer(J))){
        ind <- which(colnames(G.tilde.big)==j)
        z_j <- (C/(n))*((t(G.tilde.big[,ind,drop=F])%*%res)) + init.theta[ind]; #z_j <- as.vector(z_j)
        z_norm <- norm(z_j, "2")
        theta_j_new <- scad_firm_thresholding(z_norm, lambda1=lambda, lambda2=2, gamma=gamma)*z_j/z_norm

        resid_j_new_pos <- res - (G.tilde.big[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        #resid_j_new_neg <- neg.resid - (G.tilde.big.c[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        init.theta[ind] <- theta_j_new
        res <- c(resid_j_new_pos)
        #neg.resid <- c(resid_j_new_neg)
      }
      new.theta <- init.theta
      delta <- max(na.omit(abs(new.theta - old.theta)/abs(old.theta)))
      #print(paste("iter:", iter, "delta=",delta))
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
        w.tilde.c <-u.c.mat*W

        a.eig <- eigen(0.5*S+t(w.tilde)%*%w.tilde)
        G.tilde.big <- (a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors))

        b.eig <- eigen(0.5*S+t(w.tilde.c)%*%w.tilde.c)
        G.tilde.big.c <- (b.eig$vectors %*% diag(sqrt(b.eig$values)) %*% solve(b.eig$vectors))

        xi <- (solve(G.tilde.big)%*%t(w.tilde)%*%y.tilde)
        #xi.c <-(solve(G.tilde.big.c)%*%t(w.tilde.c)%*%y.tilde.c)

        xi.tilde <- xi <- (as.vector(xi));
        #xi.tilde.c <- xi.c <- (as.vector(xi.c))

        res <- xi.tilde - G.tilde.big %*% init.theta
        #neg.resid <- xi.tilde.c - G.tilde.big.c %*% init.theta
        #res <- pos.resid + neg.resid

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
          #resid_j_new_neg <- neg.resid - (G.tilde.big.c[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
          init.theta[ind] <- theta_j_new
          res <- c(resid_j_new_pos)
          #neg.resid <- c(resid_j_new_neg)
        }
        new.theta <- init.theta
        delta <- max(na.omit(abs(new.theta - old.theta)/abs(old.theta)))
        #print(paste("iter:", iter, "delta=",delta))
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

    b.eig <- eigen(0.5*S+t(w.tilde.c)%*%w.tilde.c)
    G.tilde.big.c <- (b.eig$vectors %*% diag(sqrt(b.eig$values)) %*% solve(b.eig$vectors))

    xi <- (solve(G.tilde.big)%*%t(w.tilde)%*%y.tilde)
    #xi.c <-(solve(G.tilde.big.c)%*%t(w.tilde.c)%*%y.tilde.c)

    xi.tilde <- xi <- (as.vector(xi));
    #xi.tilde.c <- xi.c <- (as.vector(xi.c))

    res <- xi.tilde - G.tilde.big %*% init.theta
    #neg.resid <- xi.tilde.c - G.tilde.big.c %*% init.theta
    #res <- pos.resid + neg.resid

    n.new <- length(xi.tilde)

    #grouping
    group <- colnames(G.tilde.big) <- rep(c(1:(p+1)), times=h)
    J <- max(as.integer(factor(group)))
    old.theta <- init.theta

    for(j in seq_along(integer(J))){
      ind <- which(colnames(G.tilde.big)==j)
      z_j <- (C/(n))*((t(G.tilde.big[,ind,drop=F])%*%res)) + init.theta[ind]; #z_j <- as.vector(z_j)
      z_norm <- norm(z_j, "2")
      theta_j_new <- firm_thresholding(z_norm, lambda1=lambda, lambda2=2, gamma=gamma)*z_j/z_norm

      resid_j_new_pos <- res - (G.tilde.big[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
      #resid_j_new_neg <- neg.resid - (G.tilde.big.c[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
      init.theta[ind] <- theta_j_new
      res <- c(resid_j_new_pos)
      #neg.resid <- c(resid_j_new_neg)
    }
    new.theta <- init.theta
    delta <- max(na.omit(abs(new.theta - old.theta)/abs(old.theta)))
    #print(paste("iter:", iter, "delta=",delta))
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

