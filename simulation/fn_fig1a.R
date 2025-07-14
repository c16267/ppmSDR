##################################################################
#### This code is for computation time for DC vs. GCD for figure 1(a)
##################################################################
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
library(Rfast)
library(dplyr)
library(common)
library(reporter)
library(magrittr)
library(tidyr)

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
source('fn_pplr.R')
source('fn_pplr_DC.R')


n.grid <- seq(1000, 10000, by=1000)
temp <- matrix(0,length(n.grid), 3)
temp[,2] <- c(rep(3,length(n.grid)))
temp[,3] <- n.grid

n.sim <- 100
time.list <- as.list(1:length(n.grid))

time_mat <- matrix(0, nrow=n.sim, ncol=2)
iter.fnorm <- matrix(0, nrow=n.sim, ncol=2)
colnames(time_mat) <- c("PPLR-DC","PPLR-GCD")
colnames(iter.fnorm) <- c("PPLR-DC","PPLR-GCD")


for (a in 1:nrow(temp)) {

  n <- temp[a,3] # sample size
  model <- 1
  case <- temp[a,2]
  p   <- 10
  max.iter <- maxiter <- 30
  H <- 10
  h <- 1.0e-5; eps <- 1.0e-5
  C <- 1
  delta <- 0.1
  
  # true cs
  B <- matrix(0, p, 2)
  add_term <- rnorm(p*2, 0, 1)*10^-3
  if (case == 1) {B[,c(1:2)] <- (1/sqrt(p))+add_term; q = 2}
  if (case == 2) {B[1:(p/2),1] <- B[-(1:(p/2)),2] <- 1/sqrt(p/2); q = 2}
  if (case == 3) {B[1,1] <- B[2,2] <- 1; q = 2}
  
  # initialize theta for psdr
  init.theta <- rnorm(mean=0, sd=1, n=p)
  
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
    
    
    #1. PPLR-DC
    start_time_pplr_dc <- Sys.time()
    B1 <- tryCatch(GLSCADPlogit(x, y, H, C, lambda=0.01, tol = 1.0e-4, maxiter=30)$vectors[,1:q,drop = F], error = function(e) B1)
    round(B1,2)
    end_time_pplr_dc <- Sys.time()
    time_mat[ii,1] <-as.numeric(end_time_pplr_dc - start_time_pplr_dc,units="secs")

    
    # 2. PPLR-GCD
    start_time_pplr_gcd <- Sys.time()
    B2 <- tryCatch(pplr(x, y, H, C=1, lambda=0.0001, gamma=3.7, penalty='grSCAD', max.iter, tol = 1.0e-4)$evectors[,1:q,drop = F], error = function(e) B2)
    round(B2,2)
    end_time_pplr_gcd <- Sys.time()
    time_mat[ii,2] <- as.numeric(end_time_pplr_gcd - start_time_pplr_gcd,units="secs")

    print(paste("sample size :", n.grid[a], "iteration : ", ii))
  }
  
  time.list[[a]] <- time_mat
  
}
