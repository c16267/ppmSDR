##################################################################
#L2-SCAD penalized sparse SDR with logit via D.C. algorithm
##################################################################

GLSCADPlogit <- function(x, y, H, C, lambda, tol=1.0e-4, maxiter)
{
  require(kernlab)
  n <- nrow(x)
  p <- ncol(x)
  
  # centering
  bar.x <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))
  
  # u vector Y nx1
  qprob <- seq(1/H, 1-1/H, by = 1/H)
  h <- length(qprob)
  qy <- quantile(y, qprob)
  U <- sapply(qy, function(x) 2*(x>y) - 1)
  u <- c(U)
  
  # S matrix
  Sigma <- cov(x)
  Sigma.tilde <- cbind(rep(0, p+1), rbind(rep(0, p),Sigma))
  
  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- Sigma.tilde
  S <- bdiag(temp)
  
  # X.tilde
  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- x.tilde
  X.tilde <- bdiag(temp)
  
  # initialization
  theta <- theta0 <- rep(0, h*(p+1))
  pi <- pi0 <- as.vector(1/(1 + exp(c(u) * (X.tilde %*% theta))))
  W <- W0 <- Diagonal(x=pi * (1- pi))  #time consuming
  u.tilde <- u0.tilde <- as.vector(X.tilde %*% theta + Diagonal(x=1/(pi * (1- pi))) %*% (u * pi)) #time consuming

  
  Beta <- Beta0 <- matrix(theta, p+1, h)[-1,]
  #xi <- xi0 <- apply(abs(Beta), 1, max)
  xi <- c()
  for(r in 1:nrow(Beta)){
    xi[r] <-  norm(Beta[r,], '2') # group-wise L2 penalty  
  }
  
  p2 <- function(theta, lambda)
  {
    value <- rep(0, length(theta))
    if (lambda != 0){
      temp1 <- (3.7 * lambda - theta)
      temp2 <- temp1 * (temp1 > 0)
      
      denom <- ((3.7-1)*lambda)
      value <- lambda * (1 - temp2/denom) * (theta > lambda)
    }
    value
  }
  
  for (iter in 1:maxiter)
  {
    # Dmat
    temp <- (t(X.tilde) %*% W) %*% X.tilde  #matrix multiplication time consuming
    temp.D <- C/n * temp + S   #c=1
    Dmat <- matrix(0, h*(p+1) + p, h*(p+1) + p)
    Dmat[1:(h*(p+1)), 1:(h*(p+1))] <- as.matrix(temp.D)
    
    Dmat <- Dmat + diag(rep(1.0e-10, h*(p+1) + p))
    
    # dvec
    dvec1 <- C/n * as.vector(t(u.tilde) %*% W %*% X.tilde) #c=1
    dvec2 <- -(lambda - p2(xi, lambda))
    dvec <- c(dvec1, dvec2)
    
    # Amat
    A1 <- diag(rep(c(0, rep(-1,p)), h))
    A2 <- diag(rep(c(0, rep(+1,p)), h))
    del <- 1 + (0:(h-1))*(p+1)
    A1 <- A1[,-del]
    A2 <- A2[,-del]
    A.up <- cbind(A1, A2)
    
    A3 <- diag(rep(1, p))
    A.lo <- NULL
    for (jj in 1:(2*h)) A.lo <- cbind(A3, A.lo)
    
    Amat <- rbind(A.up, A.lo)
    
    # bvec
    bvec <- rep(0, ncol(A.up))
    obj <- solve.QP(Dmat, dvec, Amat, bvec)
    Theta <- matrix(obj$solution[1:(h*(p+1))], p+1, h)
    
    new.theta <- c(Theta)
    new.xi <- obj$solution[-(1:(h*(p+1)))]
    
    #delta <- mean(abs(new.theta - theta), na.rm=T)
    delta <- mean(abs((new.theta - theta)/theta), na.rm = T)
    if (delta < tol) break
    
    theta <- new.theta
    pi <- as.vector(1/(1 + exp(c(u) * (X.tilde %*% theta))))
    W <- Diagonal(x=pi * (1 - pi))
    u.tilde <- as.vector(X.tilde %*% theta + Diagonal(x=1/(pi * (1- pi))) %*% (u * pi))
    
    xi <- new.xi
  }
  
  beta <- round(matrix(theta, ncol = h)[-1,], 8)
  
  result <- eigen(beta %*% t(beta))
  return(result)
}

