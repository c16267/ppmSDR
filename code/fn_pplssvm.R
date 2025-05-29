# ---- penalized principal least square SVM ---- #
######################################################
# Least squares PSDR + grLASSO, grMCP, grSCAD Penalty
######################################################

pplssvm <- function(x, y, H=10, C=1, lambda, gamma=3.7, penalty="grSCAD", max.iter=100)
{
  
  n <- nrow(x)
  p <- ncol(x)
  
  #step1
  #Y.tilde : dichotomize
  qprob <- seq(1/H, 1-1/H, by = 1/H)
  h <- length(qprob)
  qy <- quantile(y, qprob)
  U <- sapply(qy, function(x) 2*(x>y) - 1)
  y.tilde <- c(U)  # hxn
  
  #step2
  # x centering
  bar.x <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))
  
  
  #S : expanded covariance matrix # h(p+1) x h(p+1)
  Sigma.hat <- cov(x)
  Sigma.hat.star <- cbind(rep(0,p+1),rbind(rep(0,p), Sigma.hat))
  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- Sigma.hat.star
  S <- bdiag(temp)  # h(p+1) x h(p+1)
  
  
  #W : expanded data matrix 
  temp <- as.list(1:h)
  for (k in 1:h)  temp[[k]] <- x.tilde
  W <- bdiag(temp) # nh x h(p+1)
  
  #G.tilde  #(p+1)h x (p+1)h
  #a.eig <- eigen(0.5*S+t(W)%*%W) 
  #G.tilde.big <- C*a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors) 
  G.tilde.big <- sqrtm( (n/C)*S+t(W)%*%W )
  
  #xi # (p+1)h
  #cholesk decompositon
  #xi <- C*(Matrix::chol2inv(G.tilde.big)%*%t(W))%*%y.tilde
  xi <- sqrtm(solve((n/C)*S+t(W)%*%W))%*%t(W)%*%y.tilde
  xi.tilde <- xi <- (as.vector(xi))
  
  #restart the original process
  
  
  #grouping
  group <- colnames(G.tilde.big) <- rep(c(1:(p+1)), times=h)
  
  
  #initialization
  J <- max(as.integer(factor(group)))
  init.theta <- rnorm(h*(p+1), 0, 1)
  res <- xi.tilde - G.tilde.big %*% init.theta 
  res_vec <- c()
  
  if(penalty == "grSCAD"){
    for(iter in 1:max.iter){
      for(j in seq_along(integer(J))){
        ind <- which(colnames(G.tilde.big)==j)
        z_j <- (C/n)*(t(G.tilde.big[,ind,drop=F])%*%res + init.theta[ind]); #z_j <- as.vector(z_j)
        z_norm <- norm(z_j, "2")
        theta_j_new <- scad_firm_thresholding(z_norm, lambda1=lambda, lambda2=2, gamma)*z_j/z_norm
        resid_j_new <- res - (G.tilde.big[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        init.theta[ind] <- theta_j_new
        res <- c(resid_j_new)
      }
    }
    
  }else if(penalty == "grLasso")
    for(iter in 1:max.iter){
      for(j in seq_along(integer(J))){
        ind <- which(colnames(G.tilde.big)==j)
        z_j <- (C/n)*(t(G.tilde.big[,ind,drop=F])%*%res + init.theta[ind]); #z_j <- as.vector(z_j)
        z_norm <- norm(z_j, "2")
        theta_j_new <- soft_thresholding(z_norm, lambda=lambda)/(1 + 2)*z_j/z_norm
        resid_j_new <- res - (G.tilde.big[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
        init.theta[ind] <- theta_j_new
        res <- c(resid_j_new)
      }
    }else if(penalty == "grMCP"){
      for(iter in 1:max.iter){
        for(j in seq_along(integer(J))){
          ind <- which(colnames(G.tilde.big)==j)
          z_j <- (C/n)*(t(G.tilde.big[,ind,drop=F])%*%res + init.theta[ind]); #z_j <- as.vector(z_j)
          z_norm <- norm(z_j, "2")
          theta_j_new <- firm_thresholding(z_norm, lambda1=lambda, lambda2=2, gamma)*z_j/z_norm
          resid_j_new <- res - (G.tilde.big[,ind,drop=F])%*%(theta_j_new - init.theta[ind])
          init.theta[ind] <- theta_j_new
          res <- c(resid_j_new)
        }
      } 
    }
  
  
  #Working matrix (Mn)
  theta  <- init.theta
  
  intercept_ind <- which(((1:(h*(p+1))) %% (p+1) == 1))
  beta_mat <- matrix(theta, ncol = h)[-intercept_ind, ]
  
  Mn <- matrix(0, p, p)
  for (h in 1:length(qprob)) Mn <- Mn + beta_mat[,h, drop = F] %*% t(beta_mat[,h, drop = F])
  
  rslt <- list("evalues" = eigen(Mn)$values, "evectors"=eigen(Mn)$vectors, "x"=x)
  return(rslt)
}


# 
# pplssvm <- function(x, y, H, C, lambda=NULL, gamma, penalty, max.iter){
#   require(grpreg)
#   require(Matrix)
#   require(energy)
#   n <- nrow(x)
#   p <- ncol(x)
#   
#   #step1
#   #Y.tilde : dichotomize
#   qprob <- seq(1/H, 1-1/H, by = 1/H)
#   h <- length(qprob)
#   qy <- quantile(y, qprob)
#   U <- sapply(qy, function(x) 2*(x>y) - 1)
#   y.tilde <- c(U)  # hxn
# 
#   #step2
#   # x centering
#   bar.x <- apply(x, 2, mean)
#   x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))
#   
#   
#   #S : expanded covariance matrix # h(p+1) x h(p+1)
#   Sigma.hat <- cov(x)
#   Sigma.hat.star <- cbind(rep(0,p+1),rbind(rep(0,p), Sigma.hat))
#   temp <- as.list(1:h)
#   for (k in 1:h)  temp[[k]] <- Sigma.hat.star
#   S <- bdiag(temp)  # h(p+1) x h(p+1)
#   
#   
#   #W : expanded data matrix 
#   temp <- as.list(1:h)
#   for (k in 1:h)  temp[[k]] <- x.tilde
#   W <- bdiag(temp) # nh x h(p+1)
# 
#   
#   #G.tilde  #(p+1)h x (p+1)h
#   a.eig <- eigen(0.5*S+t(W)%*%W) 
#   G.tilde.big <- (C)*a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors) 
#   #G.tilde.big <- C*sqrtm(S+t(W)%*%W)
# 
#   #xi # (p+1)h
#   #solve
#   #xi <- C*(solve(G.tilde.big)%*%t(W))%*%y.tilde
#   #cholesky decompositon
#   xi <- (C)*(Matrix::chol2inv(G.tilde.big)%*%t(W))%*%y.tilde
#   xi.tilde <- xi <- (as.vector(xi))
#   
#   #restart the original process
#   group <- colnames(G.tilde.big) <- rep(c(1:(p+1)), times=h)
# 
#   #initialization
#   J <- max(as.integer(factor(group)))
#   
#   #setup Lambda
#   if(is.null(lambda)){
#     yy <- newY(y, family='gaussian') # yy -> centered yy, m = dimension of y
#     XG <- newXG(X=x.tilde, g=c(1:J), m=rep(1,p+1), ncolY=attr(yy, 'm'), bilevel=FALSE)
#     
#     lambdas <-  setupLambda(XG$X, yy, XG$g, family='gaussian', penalty, alpha=1, lambda.min=0.0001, log.lambda=FALSE, nlambda=20, XG$m)
#     orders <- lambdas < 1
#     lambdas <- lambdas[orders]
#     
#     dcor_vec <- c()
# 
#     for(l in 1:length(lambdas)){  
#       obj_grpreg <- grpreg(X=G.tilde.big, y=xi.tilde, group=group, penalty, family="gaussian",
#                 lambda=lambdas[l], alpha=1, eps=1e-5, max.iter=max.iter, dfmax=p,
#                gmax=length(unique(group)), gamma=gamma)
#       
#       
#       #extract betas
#       theta <- obj_grpreg$beta
#       theta  <- theta[-1]
#       
#       #drop intercept
#       intercept_ind <- which(((1:(h*(p+1))) %% (p+1) == 1))
#       beta_mat <- matrix(theta, ncol = h)[-intercept_ind, ]
#       
#       #Working matrix (Mn)
#       Mn <- matrix(0, p, p)
#       for (h in 1:length(qprob)) Mn <- Mn + beta_mat[,h, drop = F] %*% t(beta_mat[,h, drop = F])
#       tmp <- list("values" = eigen(Mn)$values, "vectors"=eigen(Mn)$vectors)
#       obj <- tmp$vectors[,1:d,drop = F]
#       
#       if(sum(tmp$vectors[,1]) != 1){
#         #value <- d2(tmp$vectors[,1:q], B[,1:q])  
#         value <- dcor(y, x %*% obj)
#         dcor_vec[l] <- value
#       }else{
#         dcor_vec[l] <- 0
#       }
#     }
#     
#     #best.lambda <- lambdas[which.min(f_norm)]
#     best.lambda <- lambdas[which.max(dcor_vec)]
#     
#     
#     obj_grpreg <- grpreg(X=G.tilde.big, y=xi.tilde, group=group, penalty, family="gaussian",
#                          lambda=best.lambda, alpha=1, eps=1e-5, max.iter=max.iter, dfmax=p,
#                          gmax=length(unique(group)), gamma=gamma)
#     
#     
#     #extract betas
#     theta <- obj_grpreg$beta
#     theta  <- theta[-1]
#     
#     #drop intercept
#     intercept_ind <- which(((1:(h*(p+1))) %% (p+1) == 1))
#     beta_mat <- matrix(theta, ncol = h)[-intercept_ind, ]
#     
#     #Working matrix (Mn)
#     Mn <- matrix(0, p, p)
#     for (h in 1:length(qprob)) Mn <- Mn + beta_mat[,h, drop = F] %*% t(beta_mat[,h, drop = F])
#     rslt <- list("lambda"=lambdas, "best.lambda"=best.lambda,"evalues" = eigen(Mn)$values, "evectors"=eigen(Mn)$vectors, "object"=obj_grpreg)
#   }else{
#     obj_grpreg <- grpreg(X=G.tilde.big, y=xi.tilde, group=group, penalty, family="gaussian",
#                          lambda=lambda, alpha=1, eps=1e-5, max.iter=max.iter, dfmax=p,
#                          gmax=length(unique(group)), gamma=gamma)
#     
#     
#     #extract betas
#     theta <- obj_grpreg$beta
#     theta  <- theta[-1]
#     
#     #drop intercept
#     intercept_ind <- which(((1:(h*(p+1))) %% (p+1) == 1))
#     beta_mat <- matrix(theta, ncol = h)[-intercept_ind, ]
#     
#     #Working matrix (Mn)
#     Mn <- matrix(0, p, p)
#     for (h in 1:length(qprob)) Mn <- Mn + beta_mat[,h, drop = F] %*% t(beta_mat[,h, drop = F])
#     rslt <- list("best.lambda"=lambda,"evalues" = eigen(Mn)$values, "evectors"=eigen(Mn)$vectors, "object"=obj_grpreg, "x"=x)
#   }
#   return(rslt)
# }

