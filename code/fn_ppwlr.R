# ---- penalized principal weighted logistic regression ---- #
ppwlr <- function(x, y, H, C, lambda, gamma, penalty, max.iter=300){
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
  S <- bdiag(temp)
  
  # X.tilde
  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- x.tilde
  W <- X.tilde <- bdiag(temp)
  
  init.theta <- theta0 <- rnorm((h*(p+1)), 0, .1)
  tol <- 10^-4
  
  if(penalty=="grSCAD"){
    for(iter in 1:max.iter){
      # initialization
      pi <- pi0 <- as.vector(1/(1 + exp(c(u) * (X.tilde %*% init.theta))))
      tmp_pi <- pi * (1- pi)*Omega_vec;
      sqrt_tmp_pi <- sqrt(tmp_pi)
      G <- G0 <- Diagonal(x=tmp_pi)  
      #u.tilde <- u0.tilde <- as.vector(W %*% init.theta + (1/tmp_pi * (u * pi)*Omega_vec))  
      u.tilde <- u0.tilde <- as.vector(W %*% init.theta*Omega_vec + (1/tmp_pi * (u * pi)))
      
      #multiply G:Hessian matrix
      u.tilde_new <- as.vector(sqrt_tmp_pi*u.tilde)
      W_new <- as.matrix(sqrt_tmp_pi*W) 
      
      a.eig <- eigen((1/(C/n))*S+t(W_new)%*%W_new )  ##time!!##
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
        #z_j <- (C/10)*(t(big_X[,ind,drop=F])%*%res)+ init.theta[ind];
        z_norm <- norm(z_j, "2")
        theta_j_new <- scad_firm_thresholding(z_norm, lambda1=lambda, lambda2=2, gamma=gamma)*z_j/z_norm
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
      G <- G0 <- Diagonal(x=tmp_pi)  
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
      G <- G0 <- Diagonal(x=tmp_pi)  
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

