##################################################################################
## This code belongs to real data analysis for Figure 3.
##################################################################################
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
library(mlbench)
library(common)
library(reporter)
library(magrittr)
library(tidyverse)
library(psvmSDR)

#change the working directory to match your PC environment.

setwd("/Users/shin/Dropbox/shared_JMShin/(ver. 240524) [JCGS] Computationally efficient algorithm for a sparse dimension reduction using a squared loss/penalized-PSDR/R/")
files.sources = list.files()
sapply(files.sources, source)



data("BostonHousing")
BostonHousing <- BostonHousing[, -c(4)]

X <- BostonHousing[,-13]
X$zn <- X$zn+0.1
Y <- BostonHousing[,"medv"]
x <- scale(X)
y <- scale(Y)


max.iter <- maxiter <- 200
H <- 5
h <- 1.0e-3; eps <- 1.0e-3
C <- 10^-1
delta <- 0.0005
lambda <- 0.01
gamma <- 3.7; penalty <- 'grSCAD'; log.lambda <- FALSE;
nfolds <- 2
n.sim <- 1

test_dcor_vec<- matrix(NA, nrow=nfolds, ncol=10)
test_dcor_mat <- matrix(NA, nrow=n.sim, ncol=10)

colnames(test_dcor_mat) <- c("No SDR","SIR","PLR","PLS_SVM","PALS_SVM","PSIR","PPLR_SVM","PPLS_SVM","PPALS_SVM", "PPL2SVM")
colnames(test_dcor_vec) <- c("No SDR","SIR","PLR","PLS_SVM","PALS_SVM","PSIR","PPLR_SVM","PPLS_SVM","PPALS_SVM", "PPL2SVM")

final_dcor_mat <- matrix(NA, nrow=dim(x)[2], ncol=(dim(x)[2]-7))

beta_list_ppls_d2 <- list(length=n.sim)
beta_list_ppals_d2 <- list(length=n.sim)
beta_list_ppls_d3 <- list(length=n.sim)
beta_list_ppals_d3 <- list(length=n.sim)



for(d in 1:(dim(x)[2]-7)){
  for(ii in 1:n.sim){
    for(ff in 1:nfolds){  
      
      set.seed(ii+ff+2024)
      
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
      
      
      #fitting
      
      # 1. SIR
      B1 <- dr(y_tr~x_tr, nslices=20, method='sir')$evectors


      # 3. PLR
      B2 <- psdr(x=x_tr, y=y_tr, h=20, lambda=lambda, eta=delta, max.iter=30, loss='logit')$evectors
      round(B2, 2)


      # 3. PLSSVM
      B3 <- psdr(x_tr, y_tr, h=20, lambda=lambda, eta=0.15, max.iter=30, loss='lssvm')$evectors
      round(B3, 2)


      # 4. PALSR
      B4 <- psdr(x=x_tr, y=y_tr, h=20, lambda=lambda, eta=0.3, max.iter=50, loss="asls")$evectors
      round(B4,2)

      ##################################################################################
      ##lambda for the penalized methods are needed to be tuned by "tune." functions.
      ##below codes are not tuned yet.
      ##################################################################################
      
      # 5.PSIR
      #tuning lambda
      nlambda <- 10; lambda.max <- 3; n.fold = 3; PSIR_grid <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)), 1.0e-4)
      best_lambda_psir <- tune.sSIR.L2(x=x_tr, y=y_tr, d=d, H, lambda.grid=PSIR_grid, n.fold=3)$opt.lambda

      B5 <- L2sparseSIR(x=x_tr, y=y_tr, H, lambda=best_lambda_psir)
      round(B5, 3)


      # 6. PPLR
      nlambda <- 2; lambda.max <- 0.01; n.fold = 3; lambda_grid_plogit <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)))
      best_lambda_plogit <- tryCatch(tune_pplr(x=x_tr, y=y_tr, d=d, H=H, C=1, lambda.grid=lambda_grid_plogit, n.fold=n.fold, penalty, max.iter=maxiter, gamma=3.7)$opt.lambda, error = function(e) best_lambda_plogit)
      best_lambda_plogit <- 0.0001
      
      B6 <- tryCatch(pplr(x=x_tr, y=y_tr, H=H, C=.1, lambda=best_lambda_plogit, gamma=3.7, penalty, max.iter=maxiter)$vectors, error = function(e) B6)
      round(B6,5)

      

      #7. PPLSM
      lambda_grid_ppls <- seq(0.0001, 0.001, length=20)
      best_lambda_ppls_nopack <- tryCatch(tune_pplssvm(x, y, d, H, C, lambda.grid=lambda_grid_ppls, gamma=3.7, penalty, max.iter=100, n.fold=3)$opt.lambda, error = function(e) best_lambda_ppls_nopack)
      B7 <- tryCatch(pplssvm(x, y, H, C, lambda=best_lambda_ppls_nopack, gamma=3.7, penalty, max.iter=100)$vectors, error = function(e) B7)
      round(B7, 2)
      

      # 8. PPAR
      lambda_grid_ppals <- seq(0.001, 0.1, length=5)
      best_lambda_ppals <- tryCatch(tune_ppalssvm(x=x_tr, y=y_tr, d=d, H, C=1, lambda.grid=lambda_grid_ppals, gamma=3.7, penalty='grSCAD', max.iter, n.fold=3)$opt.lambda, error = function(e) best_lambda_ppals)
      B8 <- ppalssvm(x=x_tr, y=y_tr, H, C=1, lambda=best_lambda_ppals, gamma=3.7, penalty='grSCAD', max.iter=maxiter)$vectors
      round(B8, 5)
      
      # 9. PPL2M
      lambda_grid_ppl2svm <- seq(0.001, 0.3, length=3)
      best_lambda_ppl2svm <- tryCatch(tune_ppl2svm(x=x_tr, y=y_tr, d=d, H, C=1, lambda.grid=lambda_grid_ppl2svm, gamma=3.7, penalty='grSCAD', max.iter, n.fold=3)$opt.lambda, error = function(e) best_lambda_ppl2svm)
      
      B9 <- ppl2svm(x=x_tr, y=y_tr, H, C=1, lambda=best_lambda_ppl2svm, gamma=3.7, penalty='grSCAD', max.iter=maxiter)$vectors
      round(B9, 5)
    
      #SIR
      sir_evec <- round(B1, 3)[,1:d]
      x.sir<- x_te %*% sir_evec

      #PLR
      plr_svm_evec <- round(B2, 3)[,1:d]
      x.plr_svm<- x_te %*% plr_svm_evec

      #PLSSVM
      pls_svm_evec <- round(B3,3)[,1:d]
      x.pls_svm<- x_te %*% pls_svm_evec


      #PALSR
      pals_svm_evec <- round(B4,3)[,1:d]
      x.pals_svm<- x_te %*% pals_svm_evec

      #PSIR
      psir_evec <- round(B5,3)[,1:d]
      x.psir <- x_te %*% psir_evec


      #PPLR
      pplr_svm_evec <- round(B6,3)[,1:d]
      x.pplr_svm <- x_te %*% pplr_svm_evec


      #PPLSM
      ppls_svm_evec <- round(B7, 5)[,1:d]
      x.ppls_svm<- x_te %*% ppls_svm_evec


      #PPAR
      ppals_svm_evec <- round(B8,3)[,1:d]
      x.ppals_svm <- x_te %*% ppals_svm_evec
      
      #PPL2M
      ppl2svm_svm_evec <- round(B9,3)[,1:d]
      x.ppl2svm <- x_te %*% ppl2svm_svm_evec
      
      #projection
      test_dcor_vec[ff,] <- abs(round(c(dcor(y_te, x_te[,1:12])$dcor, dcor(y_te, x.sir)$dcor, dcor(y_te, x.plr_svm)$dcor,dcor(y_te, x.pls_svm)$dcor,dcor(y_te, x.pals_svm)$dcor,
                                        dcor(y_te, x.psir)$dcor, dcor(y_te, x.pplr_svm)$dcor,dcor(y_te, x.ppls_svm)$dcor, dcor(y_te, x.ppals_svm)$dcor,dcor(y_te, x.ppl2svm)$dcor),3))

      
      print(paste("nfold :", ff))
      
    }
    
    test_dcor_mat[ii,] <- round(apply(test_dcor_vec,2,'mean'),3)
    print(paste("dimension:",d,"iteration :", ii))
    
    #save variable selection results.
    if(d == 3){
      beta_list_ppls_d3[[ii]] <- ppls_svm_evec
      beta_list_ppals_d3[[ii]] <- ppals_svm_evec
    }
  }
  final_dcor_mat[d,] <- apply(test_dcor_mat, 2, 'mean')
}






#################
#variable selection.
##################
#PPLSM
##################
beta_list_ppls_d3
beta_list_ppals_d3

select_ppls_beta1<- matrix(0, nrow=n.sim, ncol=12)
select_ppls_beta2<- matrix(0, nrow=n.sim, ncol=12)
select_ppls_beta3<- matrix(0, nrow=n.sim, ncol=12)

for(l in 1:n.sim){
  for(k in 1:12){
    if(abs(beta_list_ppls_d3[[l]][k,1]) > 0.0001){select_ppls_beta1[l,k] <- 1}
  }
}

for(l in 1:n.sim){
  for(k in 1:12){
    if(abs(beta_list_ppls_d3[[l]][k,2]) > 0.0001){select_ppls_beta2[l,k] <- 1}
  }
}

for(l in 1:n.sim){
  for(k in 1:12){
    if(abs(beta_list_ppls_d3[[l]][k,3]) > 0.0001){select_ppls_beta3[l,k] <- 1}
  }
}

beta1_rslt_ppls <- apply(select_ppls_beta1,2,'sum')
beta2_rslt_ppls <- apply(select_ppls_beta1,2,'sum')
beta3_rslt_ppls <- apply(select_ppls_beta1,2,'sum')

df_long_ppls <- data.frame(variable=rep(c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10",
                                          "x11","x12"),3), count=c(beta1_rslt_ppls,beta2_rslt_ppls,beta3_rslt_ppls),
                           beta=c(rep("beta1",12),rep("beta2",12),rep("beta3",12)), method=rep("PPLS-SVM", 12*3))


#########
#PPAR
#########
select_ppals_beta1<- matrix(0, nrow=n.sim, ncol=12)
select_ppals_beta2<- matrix(0, nrow=n.sim, ncol=12)
select_ppals_beta3<- matrix(0, nrow=n.sim, ncol=12)


for(l in 1:n.sim){
  for(k in 1:12){
    if(abs(beta_list_ppals_d3[[l]][k,1]) > 0.0001){select_ppals_beta1[l,k] <- 1}
  }
}

for(l in 1:n.sim){
  for(k in 1:12){
    if(abs(beta_list_ppals_d3[[l]][k,2]) > 0.0001){select_ppals_beta2[l,k] <- 1}
  }
}

for(l in 1:n.sim){
  for(k in 1:12){
    if(abs(beta_list_ppals_d3[[l]][k,3]) > 0.0001){select_ppals_beta3[l,k] <- 1}
  }
}

beta1_rslt_ppals <- apply(select_ppals_beta1, 2, 'sum')
beta2_rslt_ppals <- apply(select_ppals_beta2, 2, 'sum')
beta3_rslt_ppals <- apply(select_ppals_beta2, 2, 'sum')


df_long_ppals <- data.frame(variable=rep(c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10",
                                           "x11","x12"),3), count=c(beta1_rslt_ppals,beta2_rslt_ppals,beta3_rslt_ppals),
                            beta=c(rep("beta1",12),rep("beta2",12),rep("beta3",12)), method=rep("PPALS-SVM", 12*3))

#merge two results
df_long <- as.data.frame(rbind(df_long_ppls, df_long_ppals))
df_long$variable <- factor(df_long$variable,levels = c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10",
                                                       "x11","x12"))
df_long$method2 <- rep(c("P" %p% supsc("2") %p% "LSM" , "P" %p% supsc("2") %p% "AR"), each=36)



ggplot(data=df_long, aes(variable, count, fill=beta))+
  geom_bar(position="stack", stat="identity")+coord_flip() + 
  facet_wrap(~factor(method2,c("P" %p% supsc("2") %p% "LSM" , "P" %p% supsc("2") %p% "AR")))+
  xlab("variable")+ylab("Total number of selected")+
  labs(fill = expression(hat(bold(beta))))+
  scale_fill_discrete(labels=c(expression(hat(bold(beta))[1],hat(beta)[2],hat(beta)[3])))+
  theme(legend.key.size = unit(1.2, 'cm'), legend.title = element_text(size=15), legend.text = element_text(size=15), 
        text = element_text(size = 16), axis.text = element_text(size=16), axis.title = element_text(size=18),
        strip.text = element_text(size = 20),
        , panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.major = element_line(colour="gray", size=0.5, linetype="dashed"), panel.grid.minor = element_line(colour="gray", linetype="dashed"))




