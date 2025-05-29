sparseSIR_MAX_SCAD <- function(x, y, h, lambda, tol = 1.0e-4, maxiter=200)
{
  n <- length(y)
  p <- ncol(x)
  fk <- trans(y, h)
  
  # initialization
  beta0 <- matrix(0, p, h)
  for (j in 1:h) beta0[,j] <- solve(t(x) %*% x) %*% t(x) %*% fk[,j]
  beta <- beta0
  
  X <- matrix(0, n*h, (h+1)*p)
  for (j in 1:h) {
    row.id <- ((j-1)*n+1):(j*n)
    col.id <- ((j-1)*p+1):(j*p)
    X[row.id, col.id] <- x
  }
  
  # start updating
  for (iter in 1:maxiter) { 
    
    xi <- apply(abs(beta), 1, max) #sup-norm
    
    
    # quadratic term
    DMat <- t(X) %*% (X) + diag(rep(1.0e-10, p*(h+1)))
    
    # linear
    fkVec <- c(fk) 
    
    cVec <- c(rep(0, p*h), lambda - p2(xi, lambda))
    dVec <- (fkVec %*% X) - cVec/2
    
    # constraint
    aMat1 <- diag(rep(1, p*h))
    aMat2 <- matrix(0, p, p*h)
    
    for (j in 1:h) 
    {
      col.id <- ((j-1)*p+1):(j*p)
      diag(aMat2[1:p, col.id]) <- 1
    }
    
    aMat <- cbind(rbind(aMat1, aMat2), rbind(-aMat1, aMat2))
    bVec <- rep(0, 2*(p*h))
    
    obj <- solve.QP(DMat,dVec,aMat,bVec)
    
    new.beta <- matrix(obj$solution[1:(p*h)], p, h)
    new.xi <- obj$solution[-(1:(p*h))]
    
    delta <- mean(abs(new.beta - beta)/abs(beta)) 
    
    if (delta < tol) break
    
    beta <- new.beta
  }
  
  est <- eigen(beta %*% t(beta)/h)$vectors
  return(est)
}

L2sparseSIR <- function(x, y, h, lambda, tol = 1.0e-4, maxiter=200)
{
  n <- length(y)
  p <- ncol(x)
  fk <- trans(y, h)
  
  # initialization
  beta0 <- matrix(0, p, h)
  for (j in 1:h) beta0[,j] <- solve(t(x) %*% x) %*% t(x) %*% fk[,j]
  beta <- beta0
  
  X <- matrix(0, n*h, (h+1)*p)
  for (j in 1:h) {
    row.id <- ((j-1)*n+1):(j*n)
    col.id <- ((j-1)*p+1):(j*p)
    X[row.id, col.id] <- x
  }
  xi <- c()
  # start updating
  for (iter in 1:maxiter) {   
    
    for(r in 1:nrow(beta)){
      xi[r] <-  norm(beta[r,], '2') # group-wise L2 penalty  
    }
    
    #xi <- apply(abs(beta), 1, max) sup-norm
    
    # quadratic term
    DMat <- t(X) %*% (X) + diag(rep(1.0e-10, p*(h+1)))
    
    # linear
    fkVec <- c(fk) 
    
    cVec <- c(rep(0, p*h), lambda - p2(xi, lambda))
    dVec <- (fkVec %*% X) - cVec/2
    
    # constraint
    aMat1 <- diag(rep(1, p*h))
    aMat2 <- matrix(0, p, p*h)
    
    for (j in 1:h) 
    {
      col.id <- ((j-1)*p+1):(j*p)
      diag(aMat2[1:p, col.id]) <- 1
    }
    
    aMat <- cbind(rbind(aMat1, aMat2), rbind(-aMat1, aMat2))
    bVec <- rep(0, 2*(p*h))
    
    obj <- solve.QP(DMat,dVec,aMat,bVec)
    
    new.beta <- matrix(obj$solution[1:(p*h)], p, h)
    new.xi <- obj$solution[-(1:(p*h))]
    
    delta <- mean(abs(new.beta - beta)/abs(beta)) 
    
    if (delta < tol) break
    
    beta <- new.beta
  }
  
  est <- eigen(beta %*% t(beta)/h)$vectors
  
  return(est)
}


trans <- 
  function(y, h) 
  {
    n <- length(y)
    prob <- seq(1/h, 1-1/h, by = 1/h)
    temp <- quantile(y, prob)
    obj <- matrix(0, n, h-1)
    a <- sapply(temp, function(x) as.numeric(y < x))
    b <- sapply(temp, function(x) as.numeric(y >= x))
    
    a <- cbind(a, rep(1, n))
    b <- cbind(rep(1, n), b)
    
    v <- a*b
    
    return(v)
  }
