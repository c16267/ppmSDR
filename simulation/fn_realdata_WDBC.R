##################################################################################
## This code belongs to real data analysis for Table 3 and Figure 2.
##################################################################################
##################################################################################
##lambda for the penalized methods are needed to be tuned by "tune." functions.
##below codes are not tuned yet.
##################################################################################
rm(list = ls())
library(dr)
library(kernlab)
library(Matrix)
library(quadprog)
library(energy)
library(spls)
library(grpreg)
library(Matrix)
library(multcomp)
library(expm)
library(blockmatrix)
library(ggplot2)
library(Rfast)
library(mlbench)
library(RCurl)
library(mlbench)
library(class)
library(scatterplot3d)
library(psvmSDR)
library(common)
library(Rmisc)
library(dplyr)
library(tidyr)
library(reporter)
library(magrittr)
library(tidyverse)

#change the working directory to match your PC environment.

setwd("/Users/shin/Dropbox/shared_JMShin/(ver. 240524) [JCGS] Computationally efficient algorithm for a sparse dimension reduction using a squared loss/penalized-PSDR/R/")
files.sources = list.files()
sapply(files.sources, source)


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
nfolds <- 3
maxiter <- 100
H <- 10; h=1.0e-6; delta=0.01; eps=10^-4;C<-.1


test_error_logit_vec_save <- rep(NA, nfolds)
test_error_logit_vec_pHd<- rep(NA, nfolds)
test_error_logit_vec_wls <- rep(NA, nfolds)
test_error_logit_vec_wplr <- rep(NA, nfolds)
test_error_logit_vec_pals <- rep(NA, nfolds)
test_error_logit_vec_ppals <- rep(NA, nfolds)
test_error_logit_vec_nosdr <- rep(NA, nfolds)
test_error_logit_vec_wpplr <- rep(NA, nfolds)

test_error_list <- list(length=n.sim)
test_error_mat <- matrix(NA, nrow=n.sim, ncol=8)
lambda_list <- list(length=n.sim)
lambda_mat <- matrix(NA, nrow=n.sim, ncol=nfolds)

colnames(test_error_mat) <- c("No_SDR", "SAVE","pHd", "PWLSSVM","PWLR","PALSR","P" %p% supsc("2") %p% "AR","P" %p% supsc("2") %p% "WLR")   


for(d in 1:10){
  for(ii in 1:n.sim){
    for(ff in 1:nfolds){  
      
      set.seed(ii+ff)
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
      
      init.theta <- rnorm(ncol(x.wisc),0,1)
      
      # 0. No sdr
      # 0:benign, 1:malignant
      fit.knn <- knn(train=x_tr, test=x_te, cl=y_tr, k=7)
      test_error_logit_vec_nosdr[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
      
      #1. SAVE
      wisc.save <- dr(y_tr~x_tr, method='save')
      save <- round(wisc.save$evectors[,1:d],5)
      x.save_te <- x_te %*% save
      x.save_tr <- x_tr %*% save
      
      fit.knn <- knn(train=x.save_tr, test=x.save_te, cl=y_tr, k=7)
      test_error_logit_vec_save[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
      
      #2. pHd
      wisc.pHd <- dr(y_tr~x_tr, method='phdy')
      pHd <- round(wisc.pHd$evectors[,1:d],5)
      x.pHd_te <- x_te %*% pHd
      x.pHd_tr <- x_tr %*% pHd
      
      fit.knn <- knn(train=x.pHd_tr, test=x.pHd_te, cl=y_tr, k=7)
      test_error_logit_vec_pHd[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
      
      
      #3. PWLSSVM
      wisc.wls <- psdr(x_tr, y_tr, h=10, lambda=.1, max.iter=maxiter, loss="wlssvm")
      wls <- round(wisc.wls$evectors[,1:d],5)
      x.wls_te <- x_te %*% wls
      x.wls_tr <- x_tr %*% wls
      
      fit.knn <- knn(train=x.wls_tr, test=x.wls_te, cl=y_tr, k=7)
      test_error_logit_vec_wls[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
      
      #4.PWLR
      wisc.obj <- psdr(x_tr, y_tr, h=10, lambda=1, eta=0.1,max.iter=maxiter, loss="wlogit")
      wlogit <- round(wisc.obj$evectors[,1:d],5)
      x.wlogit_te <- x_te %*% wlogit
      x.wlogit_tr <- x_tr %*% wlogit
      
      fit.knn <- knn(train=x.wlogit_tr, test=x.wlogit_te, cl=y_tr, k=7)
      test_error_logit_vec_wplr[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
      
      
      #5. PALSR
      wisc.obj_pals <- psdr(x_tr, y_tr, h=10, lambda=.1, eta=0.1, max.iter=maxiter, loss="asls")
      wpals <- round(wisc.obj_pals$evectors,5)[,1:d]
      x.pals_te <- x_te %*% wpals
      x.pals_tr <- x_tr %*% wpals
      
      fit.knn <- knn(train=x.pals_tr, test=x.pals_te, cl=y_tr, k=7)
      test_error_logit_vec_pals[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
      
      
      #6. P2AR
      #tuning lambda#
      grid_ppals <- seq(0.0001, 0.1, length=20)
      best_lambda_ppals <- tryCatch(tune_ppalssvm(x_tr, y_tr, d=d, H=10, C=1, lambda.grid=grid_ppals, gamma=3.7, penalty='grSCAD', max.iter=maxiter, n.fold=3)$opt.lambda, error = function(e) best_lambda_ppals)
      
      
      wisc.obj_ppals <- ppalssvm(x_tr, y_tr, H=10, C=1, lambda=best_lambda_ppals, gamma=3.7, penalty='grSCAD', max.iter=maxiter)
      ppals <- round(wisc.obj_ppals$vectors, 6)[,1:d]
      x.ppals_te <- x_te %*% ppals
      x.ppals_tr <- x_tr %*% ppals
      
      fit.knn <- knn(train=x.ppals_tr, test=x.ppals_te, cl=y_tr, k=7)
      test_error_logit_vec_ppals[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
      
      
      
      #7. P2WLR
      #tuning lambda#
      grid_ppwlr <- seq(0.0001, 0.1, length=20)
      best_lambda_ppwlr <- tryCatch(tune_wpplr(x_tr, y_tr, d=d, H=10, C=1, lambda.grid=grid_ppwlr, gamma=3.7, penalty='grSCAD', max.iter=maxiter, n.fold=3)$opt.lambda, error = function(e) best_lambda_ppals)
      
      wisc.obj_wpplr <- wpplr(x_tr, y_tr, H=10, C=1, lambda=best_lambda_ppwlr, gamma=3.7, penalty='grSCAD', max.iter=maxiter, tol=1.0e-4)
      wpplr.rslt <- round(wisc.obj_wpplr$vectors, 6)[,1:d]
      x.pplr_te <- x_te %*% wpplr.rslt
      x.pplr_tr <- x_tr %*% wpplr.rslt
      
      fit.knn <- knn(train=x.pplr_tr, test=x.pplr_te, cl=y_tr, k=7)
      test_error_logit_vec_wpplr[ff] <- round(sum(y_te != fit.knn)/nrow(x_te), 5)
      
      
      
      print(paste("dimension:",d, "iteration :", ii, "fold :", ff))
      
    }
    
    test_error_mat[ii,] <- round(c( mean(test_error_logit_vec_nosdr),mean(test_error_logit_vec_save),mean(test_error_logit_vec_pHd),
                                    mean(test_error_logit_vec_wls), mean(test_error_logit_vec_wplr),
                                    mean(test_error_logit_vec_pals),mean(test_error_logit_vec_ppals), mean(test_error_logit_vec_wpplr)),4)
    
    print(paste("iteration :", ii))
  }
  test_error_list[[d]] <- test_error_mat
}

##############
#for plotting#
##############

error_mat_overall <- matrix(NA, ncol=8, nrow=n.sim);
colnames(error_mat_overall) <- c("No_SDR", "SAVE","pHd", "PWLSSVM","PWLR","PALSR","P" %p% supsc("2") %p% "AR","P" %p% supsc("2") %p% "WLR")   
for(l in 1:10){
  error_mat_overall[l,] <- apply(test_error_list[[l]], 2, 'mean')
}

error_mat_overall <- cbind(error_mat_overall, 1:10)
colnames(error_mat_overall) <- c("No_SDR", "SAVE","pHd", "PWLSSVM","PWLR","PALSR","P" %p% supsc("2") %p% "AR","P" %p% supsc("2") %p% "WLR","d")   

df.time.wide_v2 <- data.frame(error_mat_overall)
time_mat_long_v2 <- gather(df.time.wide_v2, key="method", value="test_error", 1:8, factor_key=TRUE)

ggplot(time_mat_long_v2, aes(x=d, y=test_error, colour=method, group=method)) + 
  geom_line(aes(linetype=method, color=method), size=1)+
  geom_point(aes(color=method, shape=method), size=3)+
  scale_shape_manual(values=seq(1,8))+
  ggtitle("7-NN test error rate")+ylab("Error rate")+xlab(expression(d))+
  scale_x_continuous(breaks=c(1,5,10))+
  theme(legend.key.size = unit(1.2, 'cm'), legend.title = element_text(size=16), legend.text = element_text(size=16), 
        text = element_text(size = 16), axis.text = element_text(size=16), axis.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.major = element_line(colour="gray", size=0.5, linetype="dashed"), panel.grid.minor = element_line(colour="gray", linetype="dashed"))



