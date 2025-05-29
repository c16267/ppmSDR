# ---- penalized principal weighted least square SVM ---- #
ppwlssvm <- function(x, y, H=10, C=1, lambda=NULL, gamma=3.7, penalty="grSCAD", max.iter=100){

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
  W <- X.tilde <- bdiag(temp)

  
  weight <- (1-pi.grid)*as.numeric(tmp.y==1)+(pi.grid)*(as.numeric(tmp.y==-1)) #s,s
  sqrt_weight <- sqrt(weight)
  
  
  #multiply weight 
  W <- X.tilde <- as.matrix(sqrt_weight*W)
  u.tilde <- sqrt(weight)*tmp.y
  
  #big X and Y
  a.eig <- eigen(0.5*S+t(W)%*%W )  ##time!!##
  G.tilde.big <- (C)*a.eig$vectors %*% (sqrt(a.eig$values) * solve(a.eig$vectors))
  xi <- (Matrix::chol2inv(G.tilde.big)%*%t(W))%*%u.tilde
  xi.tilde <- xi <- (as.vector(xi))
  
  #restart the original process
  group <- colnames(G.tilde.big) <- rep(c(1:(p+1)), times=h)

  #initialization
  J <- max(as.integer(factor(group)))
  
  #setup Lambda
  if(is.null(lambda)){
    yy <- newY(y, family='gaussian') # yy -> centered yy, m = dimension of y
    XG <- newXG(X=x.tilde, g=c(1:J), m=rep(1,p+1), ncolY=attr(yy, 'm'), bilevel=FALSE)
    
    lambdas <-  setupLambda(XG$X, yy, XG$g, family='gaussian', penalty, alpha=1, lambda.min=0.0001, log.lambda=FALSE, nlambda=15, XG$m)
    orders <- lambdas < 1
    lambdas <- lambdas[orders]
    
    dcor_vec <- c()

    for(l in 1:length(lambdas)){  
      obj_grpreg <- grpreg(X=G.tilde.big, y=xi.tilde, group=group, penalty, family="gaussian",
                lambda=lambdas[l], alpha=1, eps=1e-5, max.iter=max.iter, dfmax=p,
               gmax=length(unique(group)), gamma=gamma)
      
      
      #extract betas
      theta <- obj_grpreg$beta
      theta  <- theta[-1]
      
      #drop intercept
      intercept_ind <- which(((1:(h*(p+1))) %% (p+1) == 1))
      beta_mat <- matrix(theta, ncol = h)[-intercept_ind, ]
      
      #Working matrix (Mn)
      Mn <- matrix(0, p, p)
      for (h in 1:length(qprob)) Mn <- Mn + beta_mat[,h, drop = F] %*% t(beta_mat[,h, drop = F])
      tmp <- list("values" = eigen(Mn)$values, "vectors"=eigen(Mn)$vectors)
      obj <- tmp$vectors[,1:q,drop = F]
      
      if(sum(tmp$vectors[,1]) != 1){
        #value <- d2(tmp$vectors[,1:q], B[,1:q])  
        value <- dcor(y, x %*% obj)$dcor 
        dcor_vec[l] <- value
      }else{
        dcor_vec[l] <- 0
      }
    }
    
    #best.lambda <- lambdas[which.min(f_norm)]
    best.lambda <- lambdas[which.max(dcor_vec)]
    
    
    obj_grpreg <- grpreg(X=G.tilde.big, y=xi.tilde, group=group, penalty, family="gaussian",
                         lambda=best.lambda, alpha=1, eps=1e-5, max.iter=max.iter, dfmax=p,
                         gmax=length(unique(group)), gamma=gamma)
    
    
    #extract betas
    theta <- obj_grpreg$beta
    theta  <- theta[-1]
    
    #drop intercept
    intercept_ind <- which(((1:(h*(p+1))) %% (p+1) == 1))
    beta_mat <- matrix(theta, ncol = h)[-intercept_ind, ]
    
    #Working matrix (Mn)
    Mn <- matrix(0, p, p)
    for (h in 1:length(qprob)) Mn <- Mn + beta_mat[,h, drop = F] %*% t(beta_mat[,h, drop = F])
    rslt <- list("lambda"=lambdas, "best.lambda"=best.lambda,"evalues" = eigen(Mn)$values, "evectors"=eigen(Mn)$vectors, "object"=obj_grpreg)
  }else{
    obj_grpreg <- grpreg(X=G.tilde.big, y=xi.tilde, group=group, penalty, family="gaussian",
                         lambda=lambda, alpha=1, eps=1e-5, max.iter=max.iter, dfmax=p,
                         gmax=length(unique(group)), gamma=gamma)
    
    
    #extract betas
    theta <- obj_grpreg$beta
    theta  <- theta[-1]
    
    #drop intercept
    intercept_ind <- which(((1:(h*(p+1))) %% (p+1) == 1))
    beta_mat <- matrix(theta, ncol = h)[-intercept_ind, ]
    
    #Working matrix (Mn)
    Mn <- matrix(0, p, p)
    for (h in 1:length(qprob)) Mn <- Mn + beta_mat[,h, drop = F] %*% t(beta_mat[,h, drop = F])
    rslt <- list("best.lambda"=lambda,"evalues" = eigen(Mn)$values, "evectors"=eigen(Mn)$vectors, "object"=obj_grpreg, "x"=x)
  }
  return(rslt)
}


