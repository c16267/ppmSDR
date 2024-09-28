##function arguments##
# x:data matrix
# y:response (vector), should be continuous
# H:number of k for loss function
# C:cost parameter, usually set to 1
# lambda: nonconvex penalty parameter, 
#         if not specified, the function automatically choose a proper lambda via grid search
# gamma: nonconvex shape parameter
# penalty: type of penalty function: 'grSCAD', 'grLasso', 'grMCP'. 
# max.iter: maximum number of iteration
# tol: thresholding for stopping the iteration

pplr <- function(x, y, H, C,  lambda, gamma, penalty, max.iter=100, tol=1.0e-4){
  require(Matrix)
  require(expm)
  require(blockmatrix)
  require(Rfast)
  
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
  
  init.theta <- theta0 <- rep(0.2, (h*(p+1)))
  
  if(penalty=="grSCAD"){
    for(iter in 1:max.iter){
      # initialization
      pi <- pi0 <- as.vector(1/(1 + exp(c(u) * (X.tilde %*% init.theta))))
      tmp_pi <- pi * (1- pi);
      sqrt_tmp_pi <- sqrt(tmp_pi)
      G <- G0 <- Diagonal(x=tmp_pi)  
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
        G <- G0 <- Diagonal(x=tmp_pi)  
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
            G <- G0 <- Diagonal(x=tmp_pi)  
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
    
    rslt <- list("Mn"=Mn, values = eigen(Mn)$values, "vectors"=eigen(Mn)$vectors, 'x'=x)
    return(rslt)
}


