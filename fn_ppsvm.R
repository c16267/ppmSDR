# ---- penalized principal SVM (MM-GCD based) ---- #
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
  S <- bdiag(temp)
  
  # X.tilde
  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- x.tilde
  W <- X.tilde <- bdiag(temp)
  
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



