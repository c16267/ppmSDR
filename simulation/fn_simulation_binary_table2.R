
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
#setwd("/Users/shin/Desktop/Penalize PM/")


source('fn_minor_Sparse_SDR.R')
source('fn_pplssvm.R')
source('fn_ppwlssvm.R')
source('fn_pplr.R')
source('fn_ppwlr.R')
source('fn_ppasls.R')
source('fn_ppl2svm.R')
source('fn_ppwl2svm.R')
source('fn_ppsvm.R')
source('fn_ppwsvm.R')
source('fn_ppqr.R')


source('fn_tune_pplssvm.R')
source('fn_tune_ppwlssvm.R')
source('fn_tune_pplr.R')
source('fn_tune_ppwlr.R')
source('fn_tune_ppasls.R')
source('fn_tune_ppl2svm.R')
source('fn_tune_ppwl2svm.R')
source('fn_tune_ppsvm.R')
source('fn_tune_ppwsvm.R')
source('fn_tune_ppqr.R')


#####
#case
#####
case <- 3
n.sim <- 100

Fnorm.list <- list()

temp <- matrix(0,2*3, 3)
temp[,1] <- c(rep(1,3),rep(2,3))
temp[,2] <- rep(case, 6)
temp[,3] <- c(rep(c(10,20,30),2))
colnames(temp) <- c("model","case","p")
fnorm_list <- vector("list", 7)
fnorm_list_q1 <- vector("list", 7)

a=1;ii=1

for (a in 1:nrow(temp)) {
  
  n <- 1000 # sample size
  model <- temp[a,1]
  case <- temp[a,2]
  p <- temp[a,3]
  max.iter <- maxiter <- 30
  H <- 10
  h <- 1.0e-5; eps <- 1.0e-5
  C <- 1
  delta <- 0.05
  tol <- 1.0e-5
  
  # true cs
  B <- matrix(0, p, 2)
  if (case == 3) {B[1,1] <- B[2,2] <- 1; q = 2}
  
  iter.fnorm <- matrix(rep(0, 13*n.sim),ncol=13)
  colnames(iter.fnorm) <- c("save", "pwlssvm", "pwlr", "pasls", "pwl2svm", "pwsvm", "pqr",  "ppwlssvm", "ppwlr", "ppasls", "ppwl2svm", "ppwsvm", "ppqr")
  
  for (ii in 1:n.sim) {
    
    set.seed(2025 + ii)
    Sigma <- diag(1,p)
    x <- mvrnorm(n, rep(0, p), Sigma)
    eps <- rnorm(n, 0, 1)
    q = 2
    
    x1 <- x %*% B[,1]
    x2 <- x %*% B[,2]
    
    # responses
    {
      y <- if (model == 1) (x1/(0.5 + (x2 + 1)^2))+0.2*eps
      else if (model == 2) (x1+0.5)*((x2-0.5)^2)+0.2*eps
    }    
    
    y <- sign(y)
    
    # 1. SAVE
    B1 <- tryCatch(dr(y ~ x, method="save")$evectors[,1:q,drop = F], error=function(e) B1)
    round(B1,2)
    d2(B1[,1:q],B[,1:q])
    
    #1A - PWLSM   
    B1A <- tryCatch(psdr(x, y, lambda=3, h=H, loss="wlssvm")$evectors[,1:q,drop = F], error=function(e) B1A)
    round(B1A ,2)
    d2(B1A[,1:q],B[,1:q])
    
    
    #1B - PWLR   
    B1B <- tryCatch(psdr(x, y, h=H, loss="wlogit")$evectors[,1:q,drop = F], error=function(e)  B1B)
    round(B1B ,2)
    d2(B1B[,1:q],B[,1:q])
    
    #1C - PAR  
    B1C <- tryCatch(psdr(x, y, h=H, loss="asls")$evectors[,1:q,drop = F], error=function(e)  B1C)
    round(B1C ,2)
    d2(B1C[,1:q],B[,1:q])
    
    
    #1D - PWL2M
    B1D <- tryCatch(psdr(x, y, h=H, lambda = 1, loss="wl2svm")$evectors[,1:q,drop = F], error=function(e)  B1D)
    round(B1D ,2)
    d2(B1D[,1:q],B[,1:q])
    
    #1E - PWSVM
    B1E <- tryCatch(psdr(x, y, h=H, lambda = 1, loss="wsvm")$evectors[,1:q,drop = F], error=function(e)  B1E)
    round(B1E ,2)
    d2(B1E[,1:q],B[,1:q])
    
    
    #1F - PQR 
    B1F <- tryCatch(psdr(x, y, h=H, lambda = 1, loss="qr")$evectors[,1:q,drop = F], error=function(e)  B1F)
    round(B1F ,2)
    d2(B1F[,1:q],B[,1:q])
    
    
    #2. PPWLSSVM
    if(p == 30){
      C = 10
    }else{
      C = 3
    }
    lambda_grid_ppwlssvm <- seq(0.0001, 0.0001, length=20)
    #best_lambda_ppwlssvm <- tryCatch(tune.ppwlssvm(x=x, y=y, d=q, H=H, C=C, lambda.grid=lambda_grid_ppwlssvm, gamma=3.7, penalty='grSCAD', n.fold=3, max.iter)$opt.lambda, error = function(e) best_lambda_ppwlssvm)
    
    B2 <- tryCatch(ppwlssvm(x, y, H=H, C=C, lambda=0.0001, gamma=3.7, penalty='grSCAD', max.iter=max.iter)$evectors[,1:q,drop = F], error=function(e) B2)
    round(B2,2)
    d2(B2[,1:q],B[,1:q])
    
    
    #3. PPWLR
    lambda_grid_ppwlr <- seq(0.0001, 0.01, length=20)
    best_lambda_ppwlr <- tryCatch(tune.ppwlr(x=x, y=y, d=q, H=H, C=3, lambda.grid=lambda_grid_ppwlr, gamma=3.7, penalty='grSCAD', n.fold=3, max.iter)$opt.lambda, error = function(e) best_lambda_ppwlr)
    B3 <- tryCatch(ppwlr(x, y, H=H, C=3, lambda=best_lambda_ppwlr, gamma=3.7, penalty='grSCAD', max.iter=max.iter)$evectors[,1:q,drop = F], error=function(e) B3)
    round(B3,2)
    d2(B3[,1:q],B[,1:q])
    
    
    #4. PPASLS
    lambda_grid_ppasls <- seq(0.0001, 0.01, length=20)
    best_lambda_ppasls <- tryCatch(tune.ppasls(x=x, y=y, d=q, H=H, C=3, lambda.grid=lambda_grid_ppasls, gamma=3.7, penalty='grSCAD', n.fold=3, max.iter)$opt.lambda, error = function(e) best_lambda_ppasls)
    B4 <- tryCatch(ppasls(x, y, H=H, C=3, lambda=best_lambda_ppasls, gamma=3.7, penalty='grSCAD', max.iter=max.iter)$evectors[,1:q,drop = F], error=function(e) B4)
    round(B4,2)
    d2(B4[,1:q],B[,1:q])
    
    #5. PPWL2M
    lambda_grid_ppwl2svm <- seq(0.0001, 0.01, length=20)
    best_lambda_ppwl2svm <- tryCatch(tune.ppwl2svm(x=x, y=y, d=q, H=H, C=3, lambda.grid=lambda_grid_ppwl2svm, gamma=3.7, penalty='grSCAD', n.fold=3, max.iter)$opt.lambda, error = function(e) best_lambda_ppwl2svm)
    B5 <- tryCatch(ppwl2svm(x, y, H=H, C=3, lambda=best_lambda_ppwl2svm, gamma=3.7, penalty='grSCAD', max.iter=max.iter)$evectors[,1:q,drop = F], error=function(e) B5)
    round(B5,2)
    d2(B5[,1:q],B[,1:q])
    
    #6. PPWSVM
    lambda_grid_ppwsvm <- seq(0.000000001, 0.0000025, length=20)
    best_lambda_ppwsvm <- tryCatch(tune.ppwsvm(x=x, y=y, d=q, H=H, C=3, lambda.grid=lambda_grid_ppwsvm, gamma=3.7, penalty='grSCAD', n.fold=3, max.iter)$opt.lambda, error = function(e) lambda_grid_ppwsvm)
    B6 <- tryCatch(ppwsvm(x, y, H=H, C=3, lambda=best_lambda_ppwsvm, gamma=3.7, penalty='grSCAD', max.iter=max.iter)$evectors[,1:q,drop = F], error=function(e) B6)
    round(B6,2)
    d2(B6[,1:q],B[,1:q])
    
    
    #6. PPQR
    lambda_grid_ppqr <- seq(0.000000001, 0.000004, length=20)
    best_lambda_ppqr <- tryCatch(tune.ppqr(x=x, y=y, d=q, H=H, C=10, lambda.grid=lambda_grid_ppqr, gamma=3.7, penalty='grSCAD', n.fold=3, max.iter)$opt.lambda, error = function(e) lambda_grid_ppqr)
    B7 <- tryCatch(ppqr(x, y, H=H, C=10, lambda=0.000003, gamma=3.7, penalty='grSCAD', max.iter=max.iter)$evectors[,1:q,drop = F], error=function(e) B7)
    round(B7,2)
    d2(B7[,1:q],B[,1:q])
    
    
    #Frobenius norm
    iter.fnorm[ii,] <- tryCatch(round(c(d2(B1[,1:q],B[,1:q]),d2(B1A[,1:q],B[,1:q]), d2(B1B[,1:q],B[,1:q]), d2(B1C[,1:q],B[,1:q]), d2(B1D[,1:q],B[,1:q]), d2(B1E[,1:q],B[,1:q]), d2(B1F[,1:q],B[,1:q]),
                                        d2(B2[,1:q],B[,1:q]), d2(B3[,1:q],B[,1:q]), d2(B4[,1:q],B[,1:q]), d2(B5[,1:q],B[,1:q]),
                                        d2(B6[,1:q],B[,1:q]), d2(B7[,1:q],B[,1:q])), 3), error=function(e) iter.fnorm[ii-1,])
    
    
    print(paste("model:", model, "dimension:", p, "iteration: ",ii))
  }
  Fnorm.list[[a]] <- iter.fnorm 
}


