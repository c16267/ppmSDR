#################################################################
#### This code is for real data analysis for boston housing data.
#################################################################
##################################################################################
##lambda for the penalized methods are needed to be tuned by "tune." functions.
##below codes are not tuned yet.
##################################################################################

rm(list=ls())
library(dr)
library(MASS)
library(kernlab)
library(Matrix)
library(quadprog)
library(spls)
library(grpreg)
library(Matrix)
library(multcomp)
library(expm)
library(blockmatrix)
library(ggplot2)
library(energy)
library(dplyr)
library(common)
library(reporter)
library(magrittr)
library(tidyr)
library(mlbench)

#change the working directory to match your PC environment.
#set you PC working environment
#setwd("/Users/shin/Desktop/Penalize PM/")

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}


#for non-penalized sufficient dimension reduction
source('fn_minor_Sparse_SDR.R')
source('fn_pplssvm.R')
source('fn_pplssvm_nopack.R')
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

data("BostonHousing")
BostonHousing <- BostonHousing[, -c(4)]

X <- BostonHousing[,-13]
X$zn <- X$zn+0.1
Y <- BostonHousing[,"medv"]
x <- scale(X)
y <- scale(Y)


max.iter <- maxiter <- 100
H <- 10
h <- 3
C <- 1
lambda <- 0.3
gamma <- 3.7; penalty <- 'grSCAD'; log.lambda <- FALSE;
nfolds <- 3
n.sim <- 100

test_dcor_vec<- matrix(NA, nrow=nfolds, ncol=15)
test_dcor_mat <- matrix(NA, nrow=n.sim, ncol=15)
colnames(test_dcor_mat) <- c("No SDR","SIR","PLR","PLS_SVM","PALS_SVM", "PL2M", "PSVM", "PQR","PSIR","PPLR_SVM","PPLS_SVM","PPALS_SVM", "PPL2SVM", "PPSVM", "PPQR")
colnames(test_dcor_vec) <- c("No SDR","SIR","PLR","PLS_SVM","PALS_SVM", "PL2M", "PSVM", "PQR","PSIR","PPLR_SVM","PPLS_SVM","PPALS_SVM", "PPL2SVM", "PPSVM", "PPQR")

final_dcor_mat <- matrix(NA, nrow=dim(x)[2], ncol=15)
colnames(final_dcor_mat) <- c("No SDR","SIR","PLR","PLS_SVM","PALS_SVM", "PL2M", "PSVM", "PQR","PSIR","PPLR_SVM","PPLS_SVM","PPALS_SVM", "PPL2SVM", "PPSVM", "PPQR")

beta_list_ppls_d2 <- list(length=n.sim)
beta_list_ppals_d2 <- list(length=n.sim)
beta_list_pplr_d2 <- list(length=n.sim)

beta_list_ppals_d3 <- list(length=n.sim)
beta_list_ppls_d3 <- list(length=n.sim)
beta_list_pplr_d3 <- list(length=n.sim)

for(d in 1:dim(x)[2]){
  
  for(ii in 1:n.sim){
    
    for(ff in 1:nfolds){  
      
      set.seed(ii+ff+2025)
      
      #shuffle
      shuff_ind <- sample(nrow(x))
      x_shuff <- as.matrix(x[shuff_ind,])
      y_shuff <- as.matrix(y[shuff_ind])
      
      #dividen points
      folds <- cut(seq(1,nrow(x_shuff)), breaks=nfolds, labels=FALSE)
      test_error_linear <- rep(NA, nfolds)
      test_error_loess <- rep(NA, nfolds)
      test_dcor <- rep(NA, nfolds)
      
      #step1.
      ind_te <- which(folds==ff, arr.ind=TRUE)
      x_te <- x_shuff[ind_te, ] ; y_te <- y_shuff[ind_te]
      x_tr <- x_shuff[-ind_te, ] ; y_tr <- y_shuff[-ind_te]
      p <- ncol(x_tr);n<- nrow(x_tr)
      
      init.theta <- rnorm(mean=0, sd=1, n=p)
      
      
      # 1. SIR
      B1 <- dr(y_tr~x_tr, nslices=h, method='sir')$evectors
      round(B1, 2)
      
      
      # 2.PLR
      B2 <- psdr(x=x_tr, y=y_tr, h=h, lambda=0.5, max.iter=20,loss='logit')$evectors
      round(B2, 2)

      
      # 3. PLS-SVM
      B3 <- psdr(x_tr, y_tr, h=h, lambda=0.5,max.iter=20, loss='lssvm')$evectors
      round(B3, 2)
      
      # 4. PALS-SVM
      B4 <- psdr(x=x_tr, y=y_tr, h=h, lambda=0.5, max.iter=20, loss="asls")$evectors
      round(B4,2)
      
      # 5. PL2-SVM
      B5 <- psdr(x=x_tr, y=y_tr, h=h,lambda=0.3, max.iter=20, loss="l2svm")$evectors
      round(B5,2)
      
      # 6. PSVM
      B6 <- psdr(x=x_tr, y=y_tr, h=h, lambda=0.5, max.iter=20, loss="svm")$evectors
      round(B6,2)
      
      
      # 7. PQR
      B7 <- psdr(x=x_tr, y=y_tr, h=h, lambda=0.3, max.iter=50, loss="qr")$evectors
      round(B7,2)
      
      
      
      # 5.Penalized-SIR
      
      #8. PSIR
      nlambda <- 10; lambda.max <- 1; n.fold = 3;
      PSIR_grid <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)), 1.0e-4)
      best_lambda_psir <- tune.sSIR.L2(x=x_tr, y=y_tr, d=d, H, lambda.grid=PSIR_grid, n.fold=3)$opt.lambda
      B8 <- L2sparseSIR(x=x_tr, y=y_tr, H, lambda=best_lambda_psir)
      round(B8, 3)
      
      
      # 9. PPLR-SVM
      nlambda <- 20; lambda.max <- 0.000001; n.fold = 2; 
      lambda_grid_plogit <- seq(0.000015, 0.000055, length = 20)
      best_lambda_plogit <- tryCatch(tune.pplr(x=x_tr, y=y_tr, d=d, H=H, C=40, lambda.grid=lambda_grid_plogit, n.fold=n.fold, penalty = "grSCAD", max.iter=maxiter, gamma=3.7)$opt.lambda, error = function(e) best_lambda_plogit)
      B9 <- tryCatch(pplr(x=x_tr, y=y_tr, H=H, C=40, lambda=best_lambda_plogit, penalty="grSCAD", tol = 1.0e-5, max.iter=maxiter)$evectors, error = function(e) B9)
      round(B9,5)
      
      
      # 10. PPLS-SVM
      lambda_grid_ppls <- seq(0.01, 0.6, length = 10)
      best_lambda_ppls <- tryCatch(tune.pplssvm(x=x_tr, y=y_tr, d=d, H=H, C=1, lambda.grid = lambda_grid_ppls, gamma=3.7, penalty="grSCAD", max.iter=maxiter, n.fold=2)$opt.lambda, error = function(e) best_lambda_ppls)
      B10 <- tryCatch(pplssvm_nopack(x=x_tr, y=y_tr, H, C=2, lambda=best_lambda_ppls, gamma=3.7, penalty='grSCAD', max.iter)$evectors, error = function(e) B10)
      round(B10, 2)
      

      # 11. PPALS-SVM
      nlambda <- 20; lambda.max <- 0.0001; n.fold =3;
      lambda_grid_ppals <- seq(0.05, 0.2, length = 10)
      best_lambda_ppals <- tryCatch(tune.ppasls(x=x_tr, y=y_tr, d=d, H, C=1, lambda.grid=lambda_grid_ppals, gamma=3.7, penalty='grSCAD', max.iter, n.fold=2)$opt.lambda, error = function(e) best_lambda_ppals)
      B11 <- tryCatch(ppasls(x=x_tr, y=y_tr, H, C=1, lambda=best_lambda_ppals, gamma=3.7, penalty='grSCAD', max.iter=maxiter)$evectors, error = function(e) B11)
      round(B11, 5)


      # 12. PPL2-SVM
      nlambda <- 20; lambda.max <- 0.0001; n.fold =3;
      lambda_grid_ppl2svm <- seq(0.1, 0.3, length = 10)
      best_lambda_ppl2svm <- tryCatch(tune.ppl2svm(x=x_tr, y=y_tr, d=d, H, C=1, lambda.grid=lambda_grid_ppl2svm, gamma=3.7, penalty='grSCAD', max.iter, n.fold=2)$opt.lambda, error = function(e) best_lambda_ppl2svm)
      B12 <- tryCatch(ppl2svm(x=x_tr, y=y_tr, H, C=1, lambda=best_lambda_ppl2svm, gamma=3.7, penalty='grSCAD', max.iter=maxiter)$evectors, error = function(e) B12)
      round(B12, 5)
      
      
      # 13. PPSVM
      nlambda <- 20; lambda.max <- 0.000001; n.fold =3;
      lambda_grid_ppsvm <- seq(0.000015, 0.000025, length = 10)
      best_lambda_ppsvm <- tryCatch(tune.ppsvm(x=x_tr, y=y_tr, d=d, H, C=10, lambda.grid=lambda_grid_ppsvm, gamma=3.7, penalty='grSCAD', max.iter, n.fold=2)$opt.lambda, error = function(e) best_lambda_ppsvm)
      B13 <- tryCatch(ppsvm(x=x_tr, y=y_tr, H, C=3, lambda=best_lambda_ppsvm, gamma=3.7, penalty='grSCAD', max.iter=maxiter)$evectors, error = function(e) B13)
      round(B13, 5)
      

      # 14. PPQR
      nlambda <- 20; lambda.max <- 0.0000001; n.fold =3;
      lambda_grid_ppqr <- seq(0.000000015, 0.000000025, length = 20)
      best_lambda_ppqr <- tryCatch(tune.ppqr(x=x_tr, y=y_tr, d=d, H, C=10, lambda.grid=lambda_grid_ppsvm, gamma=3.7, penalty='grSCAD', max.iter, n.fold=2, tol=1.0e-5)$opt.lambda, error = function(e) best_lambda_ppqr)
      B14 <- tryCatch(ppqr(x=x_tr, y=y_tr, H, C=1, lambda=best_lambda_ppqr, gamma=3.7, penalty='grSCAD', max.iter=maxiter)$evectors, error = function(e) B14)
      round(B14, 5)
      
      
      #1.sir
      sir_evec <- round(B1, 3)[,1:d]
      x.sir<- x_te %*% sir_evec
      
      #2.plr
      plr_svm_evec <- round(B2, 3)[,1:d]
      x.plr_svm<- x_te %*% plr_svm_evec
      
      #3.pls
      pls_svm_evec <- round(B3,3)[,1:d]
      x.pls_svm<- x_te %*% pls_svm_evec
      
      #4. pals
      pals_svm_evec <- round(B4,3)[,1:d]
      x.pals_svm<- x_te %*% pals_svm_evec
      
      #5. pl2svm
      pl2svm_evec <- round(B5,3)[,1:d]
      x.pl2svm <- x_te %*% pl2svm_evec
      
      #6. psvm
      psvm_evec <- round(B6,3)[,1:d]
      x.psvm <- x_te %*% psvm_evec
      
      
      #7. pqr
      pqr_evec <- round(B7,3)[,1:d]
      x.pqr <- x_te %*% pqr_evec
      
      
      #8.psir
      psir_evec <- round(B8,3)[,1:d]
      x.psir <- x_te %*% psir_evec
      
      
      #9.pplr
      pplr_svm_evec <- round(B9,3)[,1:d]
      x.pplr_svm <- x_te %*% pplr_svm_evec
      
      
      #10.ppls
      ppls_svm_evec <- round(B10, 5)[,1:d]
      x.ppls_svm<- x_te %*% ppls_svm_evec
      
      
      #11.ppals
      ppals_svm_evec <- round(B11,3)[,1:d]
      x.ppals_svm <- x_te %*% ppals_svm_evec
      
      #12. ppl2svm
      ppl2svm_svm_evec <- round(B12,3)[,1:d]
      x.ppl2svm <- x_te %*% ppl2svm_svm_evec
      
      #13. ppsvm
      ppsvm_evec <- round(B13,3)[,1:d]
      x.ppsvm <- x_te %*% ppsvm_evec
      
      #14. ppqr
      ppqr_evec <- round(B14,3)[,1:d]
      x.ppqr <- x_te %*% ppqr_evec
      
      
      #projection
      
      test_dcor_vec[ff,] <- c(energy::dcor(y_te, x_te[,1:12]), energy::dcor(y_te, x.sir), energy::dcor(y_te, x.plr_svm),energy::dcor(y_te, x.pls_svm),energy::dcor(y_te, x.pals_svm),
                              energy::dcor(y_te, x.pl2svm), energy::dcor(y_te, x.psvm), energy::dcor(y_te, x.pqr),
                              energy::dcor(y_te, x.psir), energy::dcor(y_te, x.pplr_svm),energy::dcor(y_te, x.ppls_svm), energy::dcor(y_te, x.ppals_svm),energy::dcor(y_te, x.ppl2svm),
                              energy::dcor(y_te, x.ppsvm), energy::dcor(y_te, x.ppqr))
    

    print(paste("nfold :", ff))

    }
    
    test_dcor_mat[ii,] <- round(apply(test_dcor_vec,2,'mean'),3)
    
    
    print(paste("dimension:",d,"iteration :", ii))
    if(d == 2){
      beta_list_ppls_d2[[ii]] <- ppls_svm_evec
      beta_list_ppals_d2[[ii]] <- ppals_svm_evec
      beta_list_pplr_d2[[ii]] <- pplr_svm_evec
    }
    if(d == 3){
      beta_list_ppls_d3[[ii]] <- ppls_svm_evec
      beta_list_ppals_d3[[ii]] <- ppals_svm_evec
      beta_list_pplr_d3[[ii]] <- pplr_svm_evec
    }
  }
  
  final_dcor_mat[d,] <- apply(test_dcor_mat, 2, 'mean')
}