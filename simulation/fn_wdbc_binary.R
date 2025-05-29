#################################################################################
#### This code is for real data analysis for Dignostic Wisconsin Breast Cancer Data.
#################################################################################
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
library(mlbench)
library(RCurl)
library(mlbench)
library(class)
library(scatterplot3d)
library(psvmSDR)
library(common)
library(reporter)
library(magrittr)
library(dplyr)

#change the working directory to match your PC environment.
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

######################
### load data
######################
UCI_data_URL <- getURL('https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data')
names <- c('id_number', 'diagnosis', 'radius_mean', 
           'texture_mean', 'perimeter_mean', 'area_mean', 
           'smoothness_mean', 'compactness_mean', 
           'concavity_mean','concave_points_mean', 
           'symmetry_mean', 'fractal_dimension_mean',
           'radius_se', 'texture_se', 'perimeter_se', 
           'area_se', 'smoothness_se', 'compactness_se', 
           'concavity_se', 'concave_points_se', 
           'symmetry_se', 'fractal_dimension_se', 
           'radius_worst', 'texture_worst', 
           'perimeter_worst', 'area_worst', 
           'smoothness_worst', 'compactness_worst', 
           'concavity_worst', 'concave_points_worst', 
           'symmetry_worst', 'fractal_dimension_worst')
breast_cancer <- read.table(textConnection(UCI_data_URL), sep = ',', col.names = names)

breast_cancer$id_number <- NULL
wisc <- breast_cancer
x.wisc <- matrix(unlist(wisc[,-c(1)]), ncol = 30)
x.wisc <- scale(x.wisc)
y.wisc <- y.wisc.binary <- 2*as.numeric(as.factor(unlist(wisc[,1]))) - 3 #m=+1, b=-1
x <- x.wisc; y <- y.wisc

n.sim <- 100

maxiter <- 100
H <- 10; h=1.0e-6; eta=delta=0.01;eps=1e-05;C<-.1;
nfolds <- 3

test_error_logit_vec_nosdr <- rep(NA, nfolds)
test_error_logit_vec_save <- rep(NA, nfolds)
test_error_logit_vec_pHd<- rep(NA, nfolds)

test_error_logit_vec_wsvm <- rep(NA, nfolds)
test_error_logit_vec_wplr <- rep(NA, nfolds)
test_error_logit_vec_pals <- rep(NA, nfolds)
test_error_logit_vec_pqr <- rep(NA, nfolds)

test_error_logit_vec_ppwsvm <- rep(NA, nfolds)
test_error_logit_vec_ppwlr <- rep(NA, nfolds)
test_error_logit_vec_ppals <- rep(NA, nfolds)
test_error_logit_vec_ppqr <- rep(NA, nfolds)

test_error_list <- list(length=n.sim)
test_error_mat <- matrix(NA, nrow=n.sim, ncol=11)
lambda_list <- list(length=n.sim)
lambda_mat <- matrix(NA, nrow=n.sim, ncol=nfolds)

colnames(test_error_mat) <- c("No_SDR", "SAVE","pHd", "pwsvm","wplr", "pqr",  "par","ppwsvm", "ppwlr", "ppqr", "ppar")


for(d in 1:10){
  
for(ii in 1:n.sim){
  
  for(ff in 1:nfolds){  
    
    set.seed(ii+ff+2025)
    #shuffle
    shuff_ind <- sample(nrow(x))
    x_shuff <- as.matrix(x[shuff_ind,])
    y_shuff <- as.matrix(y[shuff_ind])
    
    #divide points
    folds <- cut(seq(1,nrow(x_shuff)), breaks=nfolds, labels=FALSE)
    test_error_knn <- rep(NA, nfolds)
    test_error_logit <- rep(NA, nfolds)
    
    #step1.
    ind_te <- which(folds==ff, arr.ind=TRUE)
    x_te <- x_shuff[ind_te, ] ; y_te <- y_shuff[ind_te]
    x_tr <- x_shuff[-ind_te, ] ; y_tr <- y_shuff[-ind_te]
    p <- ncol(x_tr);n<- nrow(x_tr)

    #color, pch setting
    colors <- rep(NA,length(y_te))
    pchs <- rep(NA,length(y_te))
    
    for(i in 1:length(y_te)){
      if(as.numeric(y_te[i] == 1)){
        colors[i] <- "red"
      }else{
        colors[i] <- "blue"
      }
    }
    for(j in 1:length(y_te)){
      if(as.numeric(y_te[j] == 1)){
        pchs[j] <- 1
      }else{
        pchs[j] <- 3
      }
    }

  init.theta <- rnorm(ncol(x.wisc), 0, 0.1)

  # 1. No sdr
  # 0:benign, 1:malignant
  fit.knn <- knn(train=x_tr, test=x_te, cl=y_tr, k=7)
  test_error_logit_vec_nosdr[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
  
  #2. SAVE
  wisc.save <- dr(y_tr~x_tr, method='save')
  save <- round(wisc.save$evectors[,1:d],5)
  x.save_te <- x_te %*% save
  x.save_tr <- x_tr %*% save
  
  fit.knn <- knn(train=x.save_tr, test=x.save_te, cl=y_tr, k=7)
  test_error_logit_vec_save[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
  
  # 3. pHd
  wisc.pHd <- dr(y_tr~x_tr, method='phdy')
  pHd <- round(wisc.pHd$evectors[,1:d],5)
  x.pHd_te <- x_te %*% pHd
  x.pHd_tr <- x_tr %*% pHd
  
  fit.knn <- knn(train=x.pHd_tr, test=x.pHd_te, cl=y_tr, k=7)
  test_error_logit_vec_pHd[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
  
  
  #4. wsvm
  wisc.wsvm <- psdr(x_tr, y_tr, init.theta, h=10, lambda=1, eps=1e-05, eta = 0.1, max.iter=20, loss="wsvm")
  wsvm <- round(wisc.wsvm$evectors[,1:d],5)
  x.wsvm_te <- x_te %*% wsvm
  x.wsvm_tr <- x_tr %*% wsvm
  
  fit.knn <- knn(train=x.wsvm_tr, test=x.wsvm_te, cl=y_tr, k=7)
  test_error_logit_vec_wsvm[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)

  
  ##5.wlogit
  
  wisc.obj <- psdr(x_tr, y_tr, init.theta, h=10, lambda=1, eps=1e-05, eta = 0.05, max.iter=20, loss="wlogit")
  
  wlogit <- round(wisc.obj$evectors[,1:d],5)
  x.wlogit_te <- x_te %*% wlogit
  x.wlogit_tr <- x_tr %*% wlogit
  
  fit.knn <- knn(train=x.wlogit_tr, test=x.wlogit_te, cl=y_tr, k=7)
  test_error_logit_vec_wplr[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
  
  
  #6. PQR
  wisc.obj <- psdr(x_tr, y_tr, init.theta, h=10, lambda=0.1, eta=0.03, eps=1e-05, max.iter=maxiter, loss="qr")
  wpqr <- round(wisc.obj$evectors[,1:d],5)
  x.wpqr_te <- x_te %*% wpqr
  x.wpqr_tr <- x_tr %*% wpqr
  
  fit.knn <- knn(train=x.wpqr_tr, test=x.wpqr_te, cl=y_tr, k=7)
  test_error_logit_vec_pqr[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
  
  
  # 7. pals
  wisc.obj_pals <- psdr(x_tr, y_tr, init.theta, h=10, lambda=1, eps=1e-03, eta=0.1, max.iter=maxiter, loss="asls")
  
  wpals <- round(wisc.obj_pals$evectors,5)[,1:d]
  x.pals_te <- x_te %*% wpals
  x.pals_tr <- x_tr %*% wpals
  
  fit.knn <- knn(train=x.pals_tr, test=x.pals_te, cl=y_tr, k=7)
  test_error_logit_vec_pals[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)

  
  #8. ppwsvm
  grid_ppwsvm <- seq(0.00000001, 0.0000008, length=20)
  best_lambda_ppwsvm <- tryCatch(tune.ppwsvm(x=x_tr, y=y_tr, d=d, H=5, C=20, lambda.grid=grid_ppwsvm, gamma=3.7, penalty='grSCAD', max.iter=100, n.fold=2)$opt.lambda, error = function(e) best_lambda_ppwsvm)
  wisc.obj_ppwsvm <- ppwsvm(x_tr, y_tr, H=10, C=20, lambda=best_lambda_ppwsvm, gamma=3.7, penalty='grSCAD', max.iter=100)
  evec_ppwsvm <- round(wisc.obj_ppwsvm$evectors, 6)[,1:d]
  x.ppwsvm_te <- x_te %*% evec_ppwsvm
  x.ppwsvm_tr <- x_tr %*% evec_ppwsvm

  fit.knn <- knn(train=x.ppwsvm_tr, test=x.ppwsvm_te, cl=y_tr, k=7)
  test_error_logit_vec_ppwsvm[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
  
  
  #9. wPPLR
  grid_ppwlr <- seq(0.02, 0.035, length=20)
  best_lambda_ppwlr <- tryCatch(tune.ppwlr(x=x_tr, y=y_tr, d=d, H=5, C=1, lambda.grid=grid_ppwlr, gamma=3.7, penalty='grSCAD', max.iter=100, n.fold=2)$opt.lambda, error = function(e) best_lambda_ppwlr)
  wisc.obj_wpplr <- ppwlr(x_tr, y_tr, H=5, C=1, lambda=best_lambda_ppwlr, gamma=3.7, penalty='grSCAD', max.iter=maxiter, tol=1.0e-4)
  evec_ppwlr <- round(wisc.obj_wpplr$evectors, 6)[,1:d]
  x.pplr_te <- x_te %*% evec_ppwlr
  x.pplr_tr <- x_tr %*% evec_ppwlr
  
  fit.knn <- knn(train=x.pplr_tr, test=x.pplr_te, cl=y_tr, k=7)
  test_error_logit_vec_ppwlr[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
  
  
  #10. ppqr
  grid_ppqr <- seq(0.02, 0.035, length=20)
  best_lambda_ppqr <- tryCatch(tune.ppqr(x=x_tr, y=y_tr, d=d, H=5, C=1, lambda.grid=grid_ppqr, gamma=3.7, penalty='grSCAD', max.iter=100, n.fold=2)$opt.lambda, error = function(e) best_lambda_ppqr)
  wisc.obj_ppqr <- ppwlr(x_tr, y_tr, H=5, C=1, lambda=best_lambda_ppqr, gamma=3.7, penalty='grSCAD', max.iter=maxiter, tol=1.0e-4)
  evec_ppqr <- round(wisc.obj_ppqr$evectors, 6)[,1:d]
  x.ppqr_te <- x_te %*% evec_ppqr
  x.ppqr_tr <- x_tr %*% evec_ppqr
  
  fit.knn <- knn(train=x.ppqr_tr, test=x.ppqr_te, cl=y_tr, k=7)
  test_error_logit_vec_ppqr[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
  
  
  #11. ppasls
  
  #tuning lambda#
  grid_ppals <-  seq(0.01, 0.025, length=20)
  best_lambda_ppals <- tryCatch(tune.ppasls(x=x_tr, y=y_tr, d=d, H=5, C=1, lambda.grid=grid_ppals, gamma=3.7, penalty='grSCAD', max.iter=100, n.fold=2)$opt.lambda, error = function(e) best_lambda_ppals)
  wisc.obj_ppals <- ppasls(x_tr, y_tr, H=5, C=1, lambda=best_lambda_ppals, gamma=3.7, penalty='grSCAD', max.iter=maxiter)
  evec_ppals <- round(wisc.obj_ppals$evectors, 6)[,1:d]
  x.ppals_te <- x_te %*% evec_ppals
  x.ppals_tr <- x_tr %*% evec_ppals

  fit.knn <- knn(train=x.ppals_tr, test=x.ppals_te, cl=y_tr, k=7)
  test_error_logit_vec_ppals[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
  
  
  print(paste("dimension:",d, ", iteration :", ii, ", fold :", ff))
  
  }

  
  test_error_mat[ii,] <- round(c( mean(test_error_logit_vec_nosdr),mean(test_error_logit_vec_save),mean(test_error_logit_vec_pHd),
                                  mean(test_error_logit_vec_wsvm), mean(test_error_logit_vec_wplr),mean(test_error_logit_vec_pals),mean(test_error_logit_vec_pqr),
                                  mean(test_error_logit_vec_ppwsvm),mean(test_error_logit_vec_ppwlr), mean(test_error_logit_vec_ppals),mean(test_error_logit_vec_ppqr)),4)
  
}
  test_error_list[[d]] <- test_error_mat
}


#plotting
error_mat_overall <- matrix(NA, ncol=11, nrow=10);
colnames(error_mat_overall) <- c("No_SDR", "SAVE","pHd", "pwsvm","wplr", "pqr",  "par","ppwsvm", "ppwlr", "ppqr", "ppar")

for(l in 1:10){
  error_mat_overall[l,] <- apply(test_error_list[[l]], 2, 'mean')
}

method_order <- c("No_SDR", "SAVE", "pHd",
                  "PWSVM", "PWLR", "PQR", "PAR",       # SDR
                  "P.WSVM", "P.WLR", "P.QR", "P.AR")   # sparse SDR
method_labels <- c(
  "No_SDR", 
  "SAVE", 
  "pHd", 
  "PWSVM", 
  "PWLR", 
  "PQR", 
  "PAR",
  expression(P^2*WSVM), 
  expression(P^2*WLR), 
  expression(P^2*QR), 
  expression(P^2*AR)
)
my_colors <- c(
  "No_SDR"="red",         # baseline
  "SAVE"="chartreuse4",           
  "pHd"="#377EB8",            
  "PWSVM"="#4DAF4A",  "P.WSVM"="#4DAF4A",  
  "PWLR"="#FF7F00",   "P.WLR"="#FF7F00",   
  "PQR"="#984EA3",    "P.QR"="#984EA3",    
  "PAR"="#00CED1",    "P.AR"="#00CED1"     
)
my_linetypes <- c(
  "No_SDR"="solid",
  "SAVE"="solid",
  "pHd"="solid",
  "PWSVM"="solid", "P.WSVM"="dashed",
  "PWLR"="solid",  "P.WLR"="dashed",
  "PQR"="solid",   "P.QR"="dashed",
  "PAR"="solid",   "P.AR"="dashed"
)
my_shapes <- c(
  "No_SDR"=1, "SAVE"=2, "pHd"=3,
  "PWSVM"=4, "P.WSVM"=8,
  "PWLR"=5, "P.WLR"=9,
  "PQR"=6, "P.QR"=10,
  "PAR"=7, "P.AR"=11
)

ggplot(time_mat_long_v2, aes(x=d, y=test_error, colour=method, group=method)) + 
  geom_line(aes(linetype=method, color=method), size=1) +
  geom_point(aes(color=method, shape=method), size=3) +
  scale_color_manual(
    values = my_colors,
    labels = method_labels
  ) +
  scale_linetype_manual(
    values = my_linetypes,
    labels = method_labels
  ) +
  scale_shape_manual(
    values = my_shapes,
    labels = method_labels
  ) +
  ggtitle("7-NN test error rate") +
  ylab("Error rate") + 
  xlab("Sufficient predictor (d)") +
  scale_x_continuous(breaks=c(1,3,5,7,9,10)) +
  theme_bw(base_size = 14) +
  theme(
    legend.key.size = unit(1.2, 'cm'),
    legend.title = element_text(size=15),
    legend.text = element_text(size=15), 
    text = element_text(size = 16),
    axis.text = element_text(size=16),
    axis.title = element_text(size=16),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.minor = element_line(color = "grey90", linetype = "dashed"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  labs(
    colour = "Method",
    linetype = "Method",
    shape = "Method"
  )
