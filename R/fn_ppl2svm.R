##function arguments##
# x:data matrix
# y:response (vector), can be either continuous or binary (+1,-1) type
# H:number of k for loss function
# C:cost parameter, usually set to 1
# lambda: nonconvex penalty parameter
# gamma: nonconvex shape parameter
# penalty: type of penalty function: 'grSCAD', 'grLasso', 'grMCP'. 
# max.iter: maximum number of iteration
  
ppl2svm <- function(x, y, H, C, lambda, gamma, penalty, max.iter=100)
{
  
  require(Matrix)
  require(expm)
  require(blockmatrix)
  require(Rfast)
  
  n <- nrow(x)
  p <- ncol(x)
  
  #step1 :Build a Design matrix
  bar.x <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))
  
  
  #step2 : we don't need to dichotomize Y
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
  S <- as.matrix(bdiag(temp))  # h(p+1) x h(p+1)
  
  #W : expanded data matrix 
  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- x.tilde
  W <- bdiag(temp) # nh x h(p+1)
  W <- as.matrix(W) 
  
  #initialization
  init.theta <- rnorm(h*(p+1), 0, 1) 
  tau.vec <- rep(0, each=n* h)
  
  
  if(penalty == "grSCAD"){
    
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
        theta_j_new <- scad_firm_thresholding(z_norm, lambda1=lambda, lambda2=2, gamma=gamma)*z_j/z_norm
        
        resid_j_new_pos <- res - (G.tilde.big[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        #resid_j_new_neg <- neg.resid - (G.tilde.big.c[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        init.theta[ind] <- theta_j_new
        res <- c(resid_j_new_pos)
        #neg.resid <- c(resid_j_new_neg)
      }
      new.theta <- init.theta
      delta <- max(na.omit(abs(new.theta - old.theta)/abs(old.theta)))
      print(paste("iter:", iter, "delta=",delta))
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
        print(paste("iter:", iter, "delta=",delta))
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
    print(paste("iter:", iter, "delta=",delta))
    if (delta < 1.0e-5) break
   }
  }
  
  #Working matrix (Mn)
  theta  <- init.theta
  
  intercept_ind <- which(((1:(h*(p+1))) %% (p+1) == 1))
  beta_mat <- matrix(theta, ncol = h)[-intercept_ind, ]
  
  Mn <- matrix(0, p, p)
  for (h in 1:length(qprob)) Mn <- Mn + beta_mat[,h, drop = F] %*% t(beta_mat[,h, drop = F])
  
  rslt <- list("values" = eigen(Mn)$values, "vectors"=eigen(Mn)$vectors, 'x'=x)
  return(rslt)
}

