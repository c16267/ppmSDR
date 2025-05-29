###############################
#### This code is for Table 1.
###############################
##################################################################################
##lambda for the penalized methods are needed to be tuned by "tune." functions.
##below codes are not tuned yet.
##################################################################################
rm(list=ls())
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
source('fn_sparse_SIR.R')

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
source('fn_tune_sparse_SIR.R')

n.sim <- 100
case <- 3  #case is either 1: dense or 3: sparse

Fnorm.temp <- as.list(1:6)

temp <- matrix(0,3*2, 3)
temp[,1] <- c(rep(1,3), rep(2,3))
temp[,2] <- rep(case, 6)
temp[,3] <- (1:3)*10
colnames(temp) <- c("model","case","p")

iter.fnorm <- matrix(rep(0, 8*n.sim), ncol=8)
colnames(iter.fnorm) <- c("SIR","PSIR", "PPLSVM", "PPLR", "PPAR", "PPL2M", "PPSVM", "PPQR")


for (a in 1:nrow(temp)) {
  
  n <- 300 # sample size
  model <- temp[a,1]
  p    <- temp[a,3]
  max.iter <- maxiter <- 100
  H <- 10
  h <- 1.0e-5; eps <- 1.0e-5
  C <- 1
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
    
    
    # 3. PPLSM
    lambda_grid_ppls <- seq(0.00001, 0.5, length=20);n.fold=3
    best_lambda_ppls_nopack <- tryCatch(tune.pplssvm(x=x, y=y, d=q, H, C=1, lambda.grid = lambda_grid_ppls, gamma=3.7, penalty="grSCAD", max.iter=maxiter, n.fold)$opt.lambda, error = function(e) best_lambda_ppls_nopack)
    B3 <- pplssvm(x, y, H, C=1, lambda=best_lambda_ppls_nopack, gamma=3.7, penalty='grSCAD', max.iter)$evectors[,1:q,drop = F]
    round(B3, 5)
    
    
    # 4. PPLR
    lambda_grid_plogit <- seq(0.00001, 0.005, length=20);n.fold=3
    best_lambda_plogit <- tryCatch(tune.pplr(x=x, y=y, d=q, H=H, C=1, lambda.grid=lambda_grid_plogit, penalty="grSCAD", gamma=3.7, n.fold=3, max.iter=maxiter)$opt.lambda, error = function(e) best_lambda_plogit)
    B4 <- pplr(x, y, H=H, C=1, lambda=best_lambda_plogit, gamma=3.7, penalty="grSCAD", max.iter=max.iter, tol=1.0e-4)$evectors[,1:q,drop = F]
    round(B4,5)
    
    
    
    # 5. PPAR
    lambda_grid_ppals <- seq(0.0001, 0.01, length=20);n.fold=3
    best_lambda_ppals <- tryCatch(tune.ppasls(x, y, d=q, H, C=.1, lambda.grid=lambda_grid_ppals, gamma=3.7, penalty='grSCAD', max.iter, n.fold=3)$opt.lambda, error = function(e) best_lambda_ppals)
    B5 <- ppasls(x, y, H, C=0.1, lambda=best_lambda_ppals, gamma=3.7, penalty='grSCAD', max.iter=maxiter)$evectors[,1:q,drop = F]
    round(B5, 5)

    
    # 6. PPL2M
    lambda_grid_ppl2svm <- seq(0.000001, 0.1, length=20);n.fold=3
    best_lambda_ppl2svm <- tryCatch(tune.ppl2svm(x=x, y=y, d=q, H, C=1, lambda.grid = lambda_grid_ppl2svm, gamma=3.7, penalty="grSCAD", max.iter=maxiter, n.fold)$opt.lambda, error = function(e) best_lambda_ppl2svm)
    B6 <- ppl2svm(x, y, H, C=1, lambda=best_lambda_ppl2svm, gamma=3.7, penalty='grSCAD', max.iter)$evectors[,1:q,drop = F]
    round(B6, 5)
    
    # 7. PPSVM
    lambda_grid_psvm <- seq(0.0000001, 0.0001, length=20);n.fold=3
    best_lambda_ppsvm <- tryCatch(tune.ppsvm(x=x, y=y, d=q, H, C=100, lambda.grid = lambda_grid_psvm, gamma=3.7, penalty="grSCAD", max.iter=maxiter, n.fold)$opt.lambda, error = function(e) best_lambda_ppsvm)
    B7 <- ppsvm(x, y, H, C=100, lambda=best_lambda_ppsvm, gamma=3.7, penalty='grSCAD', max.iter)$evectors[,1:q,drop = F]
    round(B7, 5)
    
    # 8. PPQR
    lambda_grid_ppqr <- seq(0.000001, 0.001, length=20);n.fold=3
    best_lambda_ppqr <- tryCatch(tune.ppqr(x=x, y=y, d=q, H, C=0.05, lambda.grid = lambda_grid_ppqr, gamma=3.7, penalty="grSCAD", max.iter=maxiter, n.fold)$opt.lambda, error = function(e) best_lambda_ppqr)
    B8 <- ppqr(x, y, H, C=0.05, lambda=0.00001, gamma=3.7, penalty='grSCAD', max.iter)$evectors[,1:q,drop = F]
    round(B8, 5)
    
    
    #distance correlation
    iter.fnorm[ii,] <- round(c(d2(B1[,1:q],B[,1:q]), d2(B2[,1:q],B[,1:q]), d2(B3[,1:q],B[,1:q]), d2(B4[,1:q],B[,1:q]),
                               d2(B5[,1:q],B[,1:q]), d2(B6[,1:q],B[,1:q]), d2(B7[,1:q],B[,1:q]), d2(B8[,1:q],B[,1:q])), 10)
    
    print(paste("model:", model, "dimension:", p, "iteration: ",ii))
  }
  Fnorm.temp[[a]] <- iter.fnorm
}

