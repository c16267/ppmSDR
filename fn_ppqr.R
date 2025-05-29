# ---- penalized principal quantile regression (MM-GCD based) ---- #
ppqr <- function(x, y, H=10, C=1, lambda, gamma=3.7, penalty="grSCAD", max.iter=100)
{
  n <- nrow(x)
  p <- ncol(x)

  # Step 1: generate block matrix
  # centering
  bar.x <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))

  #step2 : we don't need to dichotomize Y
  qprob <- seq(1/H, 1-1/H, by = 1/H)
  h <- length(qprob)
  tmp.y <- rep(y, times=h)


  # S matrix
  Sigma <- cov(x)  #time consuming
  Sigma.tilde <- cbind(rep(0, p+1), rbind(rep(0, p),Sigma))

  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- Sigma.tilde
  S <- bdiag(temp)

  # weighted X.tilde
  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- x.tilde
  W <- X.tilde <- bdiag(temp)

  init.theta <- theta0 <- rep(0.1, (h*(p+1)))
  tau_vec <- rep(qprob, each=n) # nhx1
  epsilon <- rep(1+0.0032*(qprob), each=n)   #propositon#3 Hunter&Lange (2000)
  tol <- 10^-5
  
  if(penalty=="grSCAD"){
    for(iter in 1:max.iter){
      r_t <- as.vector(tmp.y - (W%*%init.theta))
      #print(mean(r_t))
      K_t <- (r_t^2 + 4*((r_t)^2)*(tau_vec*(1-tau_vec)))/(r_t^2)
      omega_t <- K_t / (4*(abs(r_t)+epsilon))
      #omega_t <- omega_t/n
      sqrt_omega_t <- sqrt(omega_t)
      u.tilde <- tmp.y + abs(r_t)*(2*tau_vec -1)
      
      u.tilde <- as.vector(sqrt_omega_t*u.tilde)
      W <- X.tilde <- as.matrix(sqrt_omega_t*W)

      a.eig <- eigen(S+t(W)%*%W )  ##time!!##
      big_X <- a.eig$vectors %*% (sqrt(a.eig$values) * solve(a.eig$vectors))
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
        #print(paste("z_norm:", round(z_norm,3)))
        #print(paste("z_j:", round(z_j,3)))
      }
      res.new <- res
      new.theta <- init.theta

      delta <- mean(abs(abs(res.new)-abs(res.old)))
      if (delta < tol) break
    }
  }else if(penalty == "grLasso"){
    for(iter in 1:max.iter){
      r_t <- as.vector(tmp.y - (W%*%init.theta))
      K_t <- (r_t^2 + 4*((r_t)^2)*(tau_vec*(1-tau_vec)))/(r_t^2)
      omega_t <- K_t / (4*(abs(r_t)+epsilon))
      #omega_t <- omega_t/n
      sqrt_omega_t <- sqrt(omega_t)
      u.tilde <- tmp.y + abs(r_t)*(2*tau_vec -1)
      
      u.tilde <- as.vector(sqrt_omega_t*u.tilde)
      W <- X.tilde <- as.matrix(sqrt_omega_t*W)
      
      a.eig <- eigen(S+t(W)%*%W )  ##time!!##
      big_X <- a.eig$vectors %*% (sqrt(a.eig$values) * solve(a.eig$vectors))
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
      r_t <- as.vector(tmp.y - (W%*%init.theta))
      
      K_t <- (r_t^2 + 4*((r_t)^2)*(tau_vec*(1-tau_vec)))/(r_t^2)
      omega_t <- K_t / (4*(abs(r_t)+epsilon))
      #omega_t <- omega_t/n
      sqrt_omega_t <- sqrt(omega_t)
      u.tilde <- tmp.y + abs(r_t)*(2*tau_vec -1)
      
      u.tilde <- as.vector(sqrt_omega_t*u.tilde)
      W <- X.tilde <- as.matrix(sqrt_omega_t*W)
      
      a.eig <- eigen(S+t(W)%*%W )  ##time!!##
      big_X <- a.eig$vectors %*% (sqrt(a.eig$values) * solve(a.eig$vectors))
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
  beta_mat <- matrix(theta[-intercept_ind ], ncol = h)
  
  Mn <- matrix(0, p, p)
  for (h in 1:length(qprob)) Mn <- Mn + beta_mat[,h, drop = F] %*% t(beta_mat[,h, drop = F])
  

  rslt <- list("Mn"=Mn, "evalues" = eigen(Mn)$values, "evectors"=eigen(Mn)$vectors, 'x'=x)
  return(rslt)
}
  



