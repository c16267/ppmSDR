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


n.grid <- seq(10000, 100000, by=10000)
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

time.wide <- rbind(time.list[[1]], time.list[[2]],time.list[[3]], time.list[[4]], time.list[[5]], time.list[[6]], time.list[[7]], time.list[[8]], time.list[[9]],
                   time.list[[10]])
colnames(time.wide) <- c("P" %p% supsc("2") %p% "LR" %p% "-DC" ,"P" %p% supsc("2") %p% "LR" %p% "-GCD")   
method_name <- c("P" %p% supsc("2") %p% "LR" %p% "-DC" ,"P" %p% supsc("2") %p% "LR" %p% "-GCD")   
mean.table <- matrix(, ncol=2, nrow=length(n.grid))
var.table <- matrix(, ncol=2, nrow=length(n.grid))
n.sim <- nrow(time.list[[1]])

for(i in 1:length(n.grid)){
  mean.table[i,] <- colMeans(time.wide[((i-1)*n.sim+1):(i*n.sim),])
  var.table[i,] <- apply(time.wide[((i-1)*n.sim+1):(i*n.sim),], 2, 'sd')/sqrt(n.sim)
}


colnames(mean.table) <- method_name
colnames(var.table) <- method_name
var.table <- data.frame(var.table)
colnames(var.table) <- method_name

time.wide <- cbind(mean.table, n.grid <- seq(10000, 100000, by=10000))
colnames(time.wide) <- c(method_name, "n")

df.time.wide <- data.frame(time.wide)
colnames(df.time.wide) <- c(method_name, "n")

time_mat_long <- gather(df.time.wide, key="method", value="time", 1:2, factor_key=TRUE)



################################################################################################
par(mar=c(5,5,5,5),oma=c(1,1,1,1))
plot(df.time.wide$n, df.time.wide$`P²LR-DC`, col='blue', lty=2, type='b',
     xlab ="sample size", ylab = "Execution time in second", lwd=2.2, xaxt='n', cex.lab=1.5, pch=19, ylim=c(0, max(df.time.wide$`P²LR-DC`)))
axis(1, at=df.time.wide$n, las=1, cex=1.5, labels=c("10K", "20K", "30K", "40K", "50K", "60K","70K","80K","90K","100K"))
lines(df.time.wide$n, df.time.wide$`P²LR-GCD`, col='red', lty=3, type='b', lwd=2.2, pch=6)

lines(df.time.wide$n, df.time.wide$`P²LR-DC` + 2*var.table$`P²LR-DC`, col='blue3', lty=5, lwd=1.2)
lines(df.time.wide$n, df.time.wide$`P²LR-DC` - 2*var.table$`P²LR-DC`, col='blue3', lty=5, lwd=1.2)

lines(df.time.wide$n, df.time.wide$`P²LR-GCD` + 2*var.table$`P²LR-GCD`, col='red', lty=5, lwd=1.2)
lines(df.time.wide$n, df.time.wide$`P²LR-GCD` - 2*var.table$`P²LR-GCD`, col='red', lty=5, lwd=1.2)


grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 1)      # Grid line width
legend("topleft", legend=c(expression("P"^2~"LR-DC"), expression("P"^2~"LR-GCD")),
       col=c("blue",  "red"), lty=c(2,3), cex=1.3, lwd=2, pch=c(19,6))


################################################################################################