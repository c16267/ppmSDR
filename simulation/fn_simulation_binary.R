###############################
#### This code is for Table 2.
###############################
##################################################################################
##lambda for the penalized methods are needed to be tuned by "tune." functions.
##below codes are not tuned yet.
##################################################################################
rm(list = ls())
require(dr)
require(kernlab)
require(Matrix)
require(quadprog)
require(energy)
require(spls)
library(grpreg)
library(Matrix)
library(multcomp)
library(expm)
library(blockmatrix)
library(ggplot2)
library(Rfast)
library(psvmSDR)
library(common)

#change the working directory to match your PC environment.

setwd("/Users/shin/Dropbox/shared_JMShin/(ver. 240524) [JCGS] Computationally efficient algorithm for a sparse dimension reduction using a squared loss/penalized-PSDR/R/")
files.sources = list.files()
sapply(files.sources, source)


case <- 3 #case is either 1: dense or 3: sparse
n.sim <- 100

Fnorm.list <- as.list(1:12)


temp <- matrix(0,4*3, 3)
temp[,1] <- c(rep(1,3),rep(2,3), rep(3,3), rep(4,3))
temp[,2] <- rep(case, 12)
temp[,3] <- c(rep(c(10,20,30),4))
colnames(temp) <- c("model","case","p")

iter.fnorm <- matrix(rep(0, 5*n.sim),ncol=5)
colnames(iter.fnorm) <- c("PWLSSVM","PWLR","PALSR","P" %p% supsc("2") %p% "WLR", "P" %p% supsc("2") %p% "AR")



for (a in 1:nrow(temp)) {
  
  n <- 1000 # sample size
  model <- temp[a,1]
  case <- temp[a,2]
  p <- temp[a,3]
  max.iter <- maxiter <- 100
  H <- 5
  h <- 1.0e-5; eps <- 1.0e-5
  C <- 1
  delta <- 0.05
  tol <- 1.0e-5
  
  # true cs
  B <- matrix(0, p, 2)
  add_term <- rnorm(p*2, 0, 1)*10^-3
  if (case == 1) {B[,c(1:2)] <- (1/sqrt(p))+add_term; q = 2}
  if (case == 3) {B[1,1] <- B[2,2] <- 1; q = 2}

  
  for (ii in 1:n.sim) {
    
    set.seed(1234 + ii)
    Sigma <- diag(1,p)
    x <- mvrnorm(n, rep(0, p), Sigma)
    eps <- rnorm(n, 0, 1)
    
    x1 <- x %*% B[,1]
    x2 <- x %*% B[,2]
    
    # responses
    {
      y <- if (model == 1) (x1/(0.5 + (x2 + 1)^2))+0.2*eps
      else if (model == 2) (x1+0.5)*((x2-0.5)^2)+0.2*eps
      else if (model == 3) sin(x1)/exp(x2) + 0.2*eps
      else if (model == 4) cos(x1+1)/exp(x2) + 0.2*eps
    }    
    
    y <- sign(y)

    # # 1. PWLS-SVM
    B1 <- psdr(x, y, h=H, lambda=C, eta=delta, max.iter=maxiter, loss='wlssvm')$evectors[,1:q,drop = F]
    round(B1, 2)
    
    # 2. PWLR
    B2 <- psdr(x, y, h=H, lambda=C, eta=delta, max.iter=maxiter, loss='wlogit')$evectors[,1:q,drop = F]
    round(B2, 2)
    
    
    # 3. PALSR
    B3 <- psdr(x, y, h=H, lambda=C,  eta=delta, max.iter=maxiter, loss="asls")$evectors[,1:q,drop = F];
    round(B3,2)

    
    #4. PPWLR
    ########################
    # for automatic tuning #
    ########################
    
    if(a %in% 1:10){
      lambda_grid_pwlogit <- seq(0.001, 0.0035, length=20)
    }else{
      lambda_grid_pwlogit <- seq(0.005, 0.01, length=20)
    }
    best_lambda_pwlogit <- tryCatch(tune_wpplr(x=x, y=y, d=q, H=H, C=5, lambda.grid=lambda_grid_pwlogit, gamma=3.7, penalty='grSCAD', n.fold=3, max.iter=100)$opt.lambda, error = function(e) best_lambda_pwlogit)
    
    #best_lambda_pwlogit <- best_lambda
    B4 <- tryCatch(wpplr(x, y, H=H, C=5, lambda=best_lambda_pwlogit, gamma=3.7, penalty='grSCAD', max.iter=100, tol=1.0e-4)$vectors[,1:q,drop = F], error=function(e) B4)
    round(B4,2)
    d2(B4[,1:q],B[,1:q])
    
    
    
    # 5. PPAR
    ########################
    # for automatic tuning #
    ########################
    if(a %in% 1:10){
      lambda_grid_ppals <- seq(0.01, 0.05, length=30)
    }else{
      lambda_grid_ppals <- seq(0.006, 0.009, length=30)
    }
    best_lambda_ppals <- tryCatch(tune_ppalssvm(x, y, d=q, H=H, C=5, lambda.grid=lambda_grid_ppals, gamma=3.7, penalty='grSCAD', max.iter, n.fold=3)$opt.lambda, error = function(e) best_lambda_ppals)
    
    B5 <- tryCatch(ppalssvm(x, y, H=H, C=5, lambda=best_lambda_ppals, gamma=3.7, penalty='grSCAD', max.iter=maxiter)$vectors[,1:q,drop = F], error=function(e) B5)
    round(B5, 3)
    d2(B5[,1:q],B[,1:q])

    
    #Frobenius norm
    iter.fnorm[ii,] <- tryCatch(round(c(d2(B1[,1:q],B[,1:q]), d2(B2[,1:q],B[,1:q]), d2(B3[,1:q],B[,1:q]), d2(B4[,1:q],B[,1:q]), d2(B5[,1:q],B[,1:q])), 3), error=function(e) iter.fnorm[ii-1,])
    print(paste("model:",a, "iteration: ",ii))
  }
  Fnorm.list[[a]] <- iter.fnorm 
}
