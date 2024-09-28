###############################
#### This code is for Table 1.
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

#change the working directory to match your PC environment.

setwd("/Users/shin/Dropbox/shared_JMShin/(ver. 240524) [JCGS] Computationally efficient algorithm for a sparse dimension reduction using a squared loss/penalized-PSDR/R/")
files.sources = list.files()
sapply(files.sources, source)

n.sim <- 100
case <- 3 #case is either 1: dense or 3: sparse

Fnorm.temp <- as.list(1:4)
Fnorm.temp.sd <- as.list(1:4)

for(k in 1:4){
  Fnorm.temp[[k]] <- matrix(rep(0, 4*8), ncol=8)
  Fnorm.temp.sd[[k]] <- matrix(rep(0, 4*8), ncol=8)
}


temp <- matrix(0,3*4, 3)
temp[,1] <- c(rep(1,3), rep(2,3),rep(3,3),rep(4,3))
temp[,2] <- rep(case, 12)
temp[,3] <- (1:3)*10
colnames(temp) <- c("model","case","p")


iter.fnorm <- matrix(rep(0, 6*n.sim),ncol=6)
colnames(iter.fnorm) <- c("SIR","PSIR","PPLR", "PPLSVM", "PPL2M", "PPAR")
  
  
for (a in 1:nrow(temp)) {
    
    n <- 300 # sample size
    model <- temp[a,1]
    p    <- temp[a,3]
    max.iter <- maxiter <- 100
    H <- 10
    h <- 1.0e-5; eps <- 1.0e-5
    C <- 0.1
    delta <- 0.05
    # true cs
    B <- matrix(0, p, 2)
    add_term <- rnorm(p*2, 0, 1)*10^-3
    if (case == 1) {B[,c(1:2)] <- (1/sqrt(p))+add_term; q = 2}
    if (case == 3) {B[1,1] <- B[2,2] <- 1; q = 2}
    
    # structure dimensionality

    for (ii in 1:n.sim) {
      
      set.seed(ii)
      Sigma <- diag(1,p)
      x <- mvrnorm(n, rep(0, p), Sigma)
      eps <- rnorm(n, 0, 1)
      
      x1 <- x %*% B[,1]
      x2 <- x %*% B[,2]
      
      # responses
      {
        y <- if (model == 1) (x1/(0.5 + (x2 + 1)^2))+0.2*eps
        else if (model == 2) (x1+0.5)*(x2-0.5)^2+0.2*eps
        else if (model == 3) sin(x1)/exp(x2) + 0.2*eps
        else if (model == 4) sin(x2) + (x1+1)^2 + 0.2*eps
        else if (model == 5) cos(x1+1)/exp(x2) + 0.2*eps
      }    
      
      
      # 1. SIR
      B1 <- dr(y~x, method='sir')$evectors[,1:q,drop = F]
      round(B1, 2)
      
      
      # 2.Penalized-SIR
      # lambda tuning
      #L2 SCAD
      if(p == 10){
        if (case == 3) {nlambda <- 10; lambda.max <- 50; n.fold = 3; PSIR_grid <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)), 1.0e-4)
        }else{nlambda <- 10; lambda.max <- 50; n.fold = 3; PSIR_grid <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)), 1.0e-4)}
      }else if(p == 20){
        if (case == 3) {nlambda <- 10; lambda.max <- 50; n.fold = 3; PSIR_grid <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)), 1.0e-4)
        }else{nlambda <- 10; lambda.max <- 50; n.fold = 3; PSIR_grid <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)), 1.0e-4)}
      }else if(p == 30){
        if (case == 3) {nlambda <- 10; lambda.max <- 50; n.fold = 3; PSIR_grid <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)), 1.0e-4)
        }else{nlambda <- 10; lambda.max <- 50; n.fold = 3; PSIR_grid <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)), 1.0e-4)}
      }

      
      best_lambda_psir <- tune.sSIR.L2(x, y, d=q, H, lambda.grid=PSIR_grid, n.fold)$opt.lambda
      B2 <- L2sparseSIR(x, y, H, lambda=best_lambda_psir)[,1:q,drop = F];
      round(B2, 3)
      
      
      # 3. PPLR
      lambda_grid_plogit <- seq(0.00001, 0.005, length=20);n.fold=3
      best_lambda_plogit <- tryCatch(tune_pplr(x=x, y=y, d=q, H=H, C=1, lambda.grid=lambda_grid_plogit, penalty="grSCAD", gamma=3.7, n.fold=3, max.iter=maxiter)$opt.lambda, error = function(e) best_lambda_plogit)
      B3 <- pplr(x, y, H=H, C=1, lambda=best_lambda_plogit, gamma=3.7, penalty="grSCAD", max.iter=100, tol=1.0e-4)$vectors[,1:q,drop = F]
      round(B3,5)
      
      
      
      # 4. PPLSM
      lambda_grid_ppls <- seq(0.000001, 0.1, length=20);n.fold=3
      best_lambda_ppls_nopack <- tryCatch(tune_pplssvm(x=x, y=y, d=q, H, C=1, lambda.grid = lambda_grid_ppls, gamma=3.7, penalty="grSCAD", max.iter=maxiter, n.fold)$opt.lambda, error = function(e) best_lambda_ppls_nopack)
      B4 <- pplssvm(x, y, H, C=1, lambda=best_lambda_ppls_nopack, gamma=3.7, penalty='grSCAD', max.iter)$vectors[,1:q,drop = F]
      round(B4, 5)
      
      
      # 5. PPL2M
      lambda_grid_ppl2svm <- seq(0.000001, 0.1, length=20);n.fold=3
      best_lambda_ppl2svm <- tryCatch(tune_ppl2svm(x=x, y=y, d=q, H, C=1, lambda.grid = lambda_grid_ppl2svm, gamma=3.7, penalty="grSCAD", max.iter=maxiter, n.fold)$opt.lambda, error = function(e) best_lambda_ppls_nopack)
      B5 <- ppl2svm(x, y, H, C=1, lambda=best_lambda_ppl2svm, gamma=3.7, penalty='grSCAD', max.iter)$vectors[,1:q,drop = F]
      round(B5, 5)
      
      
      
      # 6. PPALS-SVM
      lambda_grid_ppals <- seq(0.0001, 0.1, length=20);n.fold=3
      best_lambda_ppals <- tryCatch(tune_ppalssvm(x, y, d=q, H, C=.1, lambda.grid=lambda_grid_ppals, gamma=3.7, penalty='grSCAD', max.iter, n.fold=3)$opt.lambda, error = function(e) best_lambda_ppals)
      B6 <- ppalssvm(x, y, H, C=.1, lambda=best_lambda_ppals, gamma=3.7, penalty='grSCAD', max.iter=maxiter)$vectors[,1:q,drop = F]
      round(B6, 5)

      
      #distance correlation
      iter.fnorm[ii,] <- round(c(d2(B1[,1:q],B[,1:q]), d2(B2[,1:q],B[,1:q]), d2(B3[,1:q],B[,1:q]), d2(B4[,1:q],B[,1:q]),
                                 d2(B5[,1:q],B[,1:q]), d2(B6[,1:q],B[,1:q])), 5)
      
      print(paste("case: ",case, "model and dim: ",a, "iter: ", ii))
    }
}

  