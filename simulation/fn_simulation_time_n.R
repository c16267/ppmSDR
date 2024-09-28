###############################
#### This code is for Figure 1.
###############################
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
library(dplyr)
library(tidyr)
library(Rmisc)
library(common)
library(reporter)
library(magrittr)
library(tidyverse)

#change the working directory to match your PC environment.
setwd("/Users/shin/Dropbox/shared_JMShin/(ver. 240524) [JCGS] Computationally efficient algorithm for a sparse dimension reduction using a squared loss/penalized-PSDR/R/")
files.sources = list.files()
sapply(files.sources, source)


n.sim <- 100
n.grid <- c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
temp <- matrix(0,length(n.grid), 3)
temp[,2] <- c(rep(3,length(n.grid)))
temp[,3] <- n.grid
time.list <- as.list(1:length(n.grid))

time_mat <- matrix(0, nrow=n.sim, ncol=4)
iter.fnorm <- matrix(0, nrow=n.sim, ncol=4)
colnames(time_mat) <- c("P" %p% supsc("2") %p% "LR-DC", "P" %p% supsc("2") %p% "LR-GCD", "P" %p% supsc("2") %p% "LSM", "P" %p% supsc("2") %p% "AR")
colnames(iter.fnorm) <- c("P" %p% supsc("2") %p% "LR-DC", "P" %p% supsc("2") %p% "LR-GCD", "P" %p% supsc("2") %p% "LSM", "P" %p% supsc("2") %p% "AR")



for (a in 1:nrow(temp)) {
  n <- temp[a,3] # sample size
  model <- 1
  case <- temp[a,2]
  p   <- 10
  max.iter <- maxiter <- 30
  H <- 3
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
      else if (model == 3) 0.5*x1^3 + (2*x2+3)*eps #asymmetric Dong #3
    }   
    
    
    # 1. PPLR-DC
    start_time_pplr_dc <- Sys.time()
    B1 <- tryCatch(SCADPlogit(x, y, H, C, lambda=0.01, tol = 1.0e-4, maxiter=max.iter)$vectors[,1:q,drop = F], error = function(e) B1)
    round(B1,2)
    end_time_pplr_dc <- Sys.time()
    time_mat[ii,1] <-as.numeric(end_time_pplr_dc - start_time_pplr_dc,units="secs")
    
    
    # 2. PPLR-GCD
    start_time_pplr_gcd <- Sys.time()
    B2 <- tryCatch(pplr(x, y, H, C, lambda=0.0001, gamma=3.7, penalty='grSCAD', max.iter=max.iter, tol = 1.0e-4)$vectors[,1:q,drop = F], error = function(e) B2)
    round(B2,2)
    end_time_pplr_gcd <- Sys.time()
    time_mat[ii,2] <- as.numeric(end_time_pplr_gcd - start_time_pplr_gcd,units="secs")
    
    # 3. PPLSM
    start_time_ppls <- Sys.time()
    B3 <- tryCatch(pplssvm(x, y, H, C=C, lambda=0.005, gamma=3.7, penalty='grSCAD', max.iter=maxiter)$vectors[,1:q,drop = F], error = function(e) B3)
    end_time_ppls <- Sys.time()
    time_mat[ii, 3] <- as.numeric(end_time_ppls - start_time_ppls, units='secs')
    round(B3, 2)
    
    
    # 4. PPAR
    start_time_ppals <- Sys.time()
    B4 <- tryCatch(ppalssvm(x, y, H, C, lambda=0.005446, gamma=3.7, penalty='grSCAD', max.iter=20)$vectors[,1:q,drop = F], error = function(e) B4)
    end_time_ppals <- Sys.time()
    time_mat[ii, 4] <- as.numeric(end_time_ppals - start_time_ppals, units='secs')
    round(B4, 2)
    
    print(paste("sample size :", n.grid[a], "iteration : ", ii))
  }
  time.list[[a]] <- time_mat
}

########
#plotting
########
time.wide <- rbind(time.list[[1]], time.list[[2]],time.list[[3]], time.list[[4]], time.list[[5]], time.list[[6]], time.list[[7]], time.list[[8]], time.list[[9]],
                   time.list[[10]])

mean.table <- matrix(, ncol=4, nrow=10)
var.table <- matrix(, ncol=4, nrow=10)

n.sim <- nrow(time.list[[1]])

for(i in 1:10){
  mean.table[i,] <- colMeans(time.wide[((i-1)*n.sim+1):(i*n.sim),])
  var.table[i,] <- apply(time.wide[((i-1)*n.sim+1):(i*n.sim),], 2, 'sd')/sqrt(n.sim)
}

mean.table <- data.frame(mean.table)
var.table <- data.frame(var.table)

time.wide <- cbind(mean.table, c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))

colnames(time.wide) <- c("P2LR_DC", "P2LR_GCD", "P2LSM", "P2AR", "n")
colnames(var.table) <- c("P2LR_DC", "P2LR_GCD", "P2LSM", "P2AR")

df.time.wide <- data.frame(time.wide)

par(mar=c(5,5,5,5),oma=c(1,1,1,1))
plot(df.time.wide$n, df.time.wide$P2LR_DC, col='blue', lty=2, type='b',
     xlab ="sample size", ylab = "Execution time in second", lwd=2.2, xaxt='n', cex.lab=1.5, pch=19)
axis(1, at=df.time.wide$n, las=1, cex=1.5, labels=c("1K", "2K", "3K", "4K", "5K", "6K","7K","8K","9K","10K"))
lines(df.time.wide$n, df.time.wide$P2LR_GCD, col='red', lty=3, type='b', lwd=2.2, pch=6)
lines(df.time.wide$n, df.time.wide$P2LSM, col='forestgreen', lty=4, type='b', lwd=2.2, pch=4)
lines(df.time.wide$n, df.time.wide$P2AR, col='violet', lty=5, type='b', lwd=2.2, pch=5)

lines(df.time.wide$n, df.time.wide$P2LR_DC + 2*var.table$P2LR_DC, col='blue3', lty=5, lwd=1.2)
lines(df.time.wide$n, df.time.wide$P2LR_DC - 2*var.table$P2LR_DC, col='blue3', lty=5, lwd=1.2)

lines(df.time.wide$n, df.time.wide$P2LR_GCD + 2*var.table$P2LR_GCD, col='red', lty=5, lwd=1.2)
lines(df.time.wide$n, df.time.wide$P2LR_GCD - 2*var.table$P2LR_GCD, col='red', lty=5, lwd=1.2)

lines(df.time.wide$n, df.time.wide$P2LSM + 2*var.table$P2LSM, col='forestgreen', lty=5, lwd=1.2)
lines(df.time.wide$n, df.time.wide$P2LSM - 2*var.table$P2LSM, col='forestgreen', lty=5, lwd=1.2)

lines(df.time.wide$n, df.time.wide$P2AR + 2*var.table$P2AR, col='violet', lty=5, lwd=1.2)
lines(df.time.wide$n, df.time.wide$P2AR - 2*var.table$P2AR, col='violet', lty=5, lwd=1.2)

grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 1)      # Grid line width
legend("topleft", legend=c("P" %p% supsc("2") %p% "LR-DC", "P" %p% supsc("2") %p% "LR-GCD", "P" %p% supsc("2") %p% "LSM", "P" %p% supsc("2") %p% "AR"),
       col=c("blue", "violet", "red", "forestgreen"), lty=c(2,5,3,4), cex=1.4, lwd=2, pch=c(19,5,6,4))



