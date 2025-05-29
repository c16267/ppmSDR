##################################################################
#### This code is for computation time for LPM-type P2Ms 
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


n.grid <- round(exp(seq(log(10000), log(1000000), length.out = 8)))
temp <- matrix(0,length(n.grid), 3)
temp[,2] <- c(rep(3,length(n.grid)))
temp[,3] <- n.grid

n.sim <- 100
time.list <- as.list(1:length(n.grid))

time_mat <- matrix(0, nrow=n.sim, ncol=6)
iter.fnorm <- matrix(0, nrow=n.sim, ncol=6)
colnames(time_mat) <-  c("ppwlssvm", "ppwlr", "ppasls", "ppwl2svm", "ppwsvm", "ppqr")
colnames(iter.fnorm) <- c("ppwlssvm", "ppwlr", "ppasls", "ppwl2svm", "ppwsvm", "ppqr")


for (a in 1:nrow(temp)) {
  n <- temp[a,3] # sample size
  model <- 1
  case <- temp[a,2]
  p   <- 10
  max.iter <- maxiter <- 20
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
    
    y <- sign(y)
    
    #1. SAVE
    # start_time_save <- Sys.time()
    # B1 <- tryCatch(dr(y ~ x, method ="save")$evectors[,1:q,drop = F], error = function(e) B1)
    # round(B1,2)
    # end_time_save <- Sys.time()
    # time_mat[ii,1] <- as.numeric(end_time_save - start_time_save,units="secs")
    
    # 2. ppwlssvm
    start_time_ppwlssvm <- Sys.time()
    B2 <- tryCatch(ppwlssvm(x, y, H, C, lambda=0.002, gamma=3.7, penalty='grSCAD', max.iter)$evectors[,1:q,drop = F], error = function(e) B2)
    round(B2,2)
    end_time_ppwlssvm <- Sys.time()
    time_mat[ii,1] <- as.numeric(end_time_ppwlssvm - start_time_ppwlssvm,units="secs")
    
    # 3. ppwlr
    start_time_ppwlr <- Sys.time()
    B3 <- tryCatch(ppwlr(x, y, H, C=5, lambda=0.0002, gamma=3.7, penalty='grSCAD', max.iter)$evectors[,1:q,drop = F], error = function(e) B3)
    end_time_ppwlr <- Sys.time()
    time_mat[ii, 2] <- as.numeric(end_time_ppwlr - start_time_ppwlr, units='secs')
    round(B3, 2)
    
    
    # 4. ppasls
    start_time_ppals <- Sys.time()
    B4 <- tryCatch(ppasls(x, y, H, C=5, lambda=0.03, gamma=3.7, penalty='grSCAD', max.iter)$evectors[,1:q,drop = F], error = function(e) B4)
    end_time_ppals <- Sys.time()
    time_mat[ii, 3] <- as.numeric(end_time_ppals - start_time_ppals, units='secs')
    round(B4, 2)
    
    # 5. ppwl2svm
    start_time_ppwl2svm <- Sys.time()
    B5 <- tryCatch(ppwl2svm(x, y, H, C=5, lambda=0.005, gamma=3.7, penalty='grSCAD', max.iter)$evectors[,1:q,drop = F], error = function(e) B5)
    end_time_ppwl2svm <- Sys.time()
    time_mat[ii, 4] <- as.numeric(end_time_ppwl2svm - start_time_ppwl2svm, units='secs')
    round(B5, 2)
    
    
    # 6. ppwsvm
    start_time_ppwsvm <- Sys.time()
    B6 <- tryCatch(ppwsvm(x, y, H, C=10,  lambda = 0.000001, gamma=3.7, penalty="grSCAD", max.iter, tol=10^-6)$evectors[,1:q,drop = F], error = function(e) B6)
    end_time_ppwsvm <- Sys.time()
    time_mat[ii, 5] <- as.numeric(end_time_ppwsvm - start_time_ppwsvm, units='secs')
    round(B6, 2)
    
    
    # 7. ppqr
    start_time_ppqr <- Sys.time()
    B7 <- tryCatch(ppqr(x, y, H, C=10,  lambda = 0.000008, gamma=3.7, penalty="grSCAD", max.iter, tol=10^-6)$evectors[,1:q,drop = F], error = function(e) B7)
    end_time_ppqr <- Sys.time()
    time_mat[ii, 6] <- as.numeric(end_time_ppqr - start_time_ppqr, units='secs')
    round(B7, 2)
    
    print(paste("sample size :", n.grid[a], "iteration : ", ii))
  }
  
  time.list[[a]] <- time_mat
  
}

###############
# --- plotting
###############
n.grid <- round(exp(seq(log(10000), log(1000000), length.out = 8)))
time.wide <- rbind(time.list[[1]], time.list[[2]],time.list[[3]], time.list[[4]], time.list[[5]], time.list[[6]], time.list[[7]], time.list[[8]])
colnames(time.wide) <- c("P" %p% supsc("2") %p% "WLSM", "P" %p% supsc("2") %p% "WLR","P" %p% supsc("2") %p% "AR", "P" %p% supsc("2") %p% "WL2M", "P" %p% supsc("2") %p% "WSVM", "P" %p% supsc("2") %p% "QR")  
method_name <- c("P" %p% supsc("2") %p% "WLSM", "P" %p% supsc("2") %p% "WLR","P" %p% supsc("2") %p% "AR", "P" %p% supsc("2") %p% "WL2M", "P" %p% supsc("2") %p% "WSVM", "P" %p% supsc("2") %p% "QR")
mean.table <- matrix(, ncol=6, nrow=length(n.grid))
var.table <- matrix(, ncol=6, nrow=length(n.grid))
n.sim <- nrow(time.list[[1]])

for(i in 1:length(n.grid)){
  mean.table[i,] <- colMeans(time.wide[((i-1)*n.sim+1):(i*n.sim),])
  var.table[i,] <- apply(time.wide[((i-1)*n.sim+1):(i*n.sim),], 2, 'sd')/sqrt(n.sim)
}


colnames(mean.table) <- method_name
colnames(var.table) <- method_name
var.table <- data.frame(var.table)
colnames(var.table) <- method_name

time.wide <- cbind(mean.table, exp(seq(log(10000), log(1000000), length.out = 8)))
colnames(time.wide) <- c(method_name, "n")
df.time.wide <- data.frame(time.wide)
colnames(df.time.wide) <- c(method_name, "n")
time_mat_long <- gather(df.time.wide, key="method", value="time", 1:6, factor_key=TRUE)

###############
# --- plotting
###############
par(mar=c(5,5,5,5),oma=c(1,1,1,1))
method_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")
#method_colors <- c("#e7298a", "#7570b3", "#66a61e", "#e6ab02", "#1b9e77", "#d95f02")
plot(df.time.wide$n, df.time.wide$`P²WLSM`, col=method_colors[2], lty=2, type='b',
     xlab = expression("Sample size"), ylab = "Execution time in second", lwd=2.2, xaxt='n', cex.lab=1.5, pch=2, ylim=c(0, max(df.time.wide$`P²WL2M`)))
custom_ticks <- c(10000, 100000, 200000, 400000, 600000, 800000, 1000000)
custom_labels <- c("10K", "100K", "200K", "400K", "600K", "800K", "1M")
axis(1, at=custom_ticks, labels=custom_labels, las=1, cex=1.5)

lines(df.time.wide$n, df.time.wide$`P²WLR`, col=method_colors[1], lty=1, type='b', lwd=2.2, pch=1)
lines(df.time.wide$n, df.time.wide$`P²AR`, col=method_colors[3], lty=3, type='b', lwd=2.2, pch=3)
lines(df.time.wide$n, df.time.wide$`P²WL2M`, col=method_colors[4], lty=4, type='b', lwd=2.2, pch=4)
lines(df.time.wide$n, df.time.wide$`P²WSVM`, col=method_colors[5], lty=5, type='b', lwd=2.2, pch=5)
lines(df.time.wide$n, df.time.wide$`P²QR`, col=method_colors[6], lty=6, type='b', lwd=2.2, pch=6)


lines(df.time.wide$n, df.time.wide$`P²WLSM` + 2*var.table$`P²WLSM`, col=method_colors[2], lty=2, lwd=1.2)
lines(df.time.wide$n, df.time.wide$`P²WLSM` - 2*var.table$`P²WLSM`, col=method_colors[2], lty=2, lwd=1.2)

lines(df.time.wide$n, df.time.wide$`P²WLR` + 2*var.table$`P²WLR`, col=method_colors[1], lty=1, lwd=1.2)
lines(df.time.wide$n, df.time.wide$`P²WLR` - 2*var.table$`P²WLR`, col=method_colors[1], lty=1, lwd=1.2)

lines(df.time.wide$n, df.time.wide$`P²AR` + 2*var.table$`P²AR`, col=method_colors[3], lty=3, lwd=1.2)
lines(df.time.wide$n, df.time.wide$`P²AR` - 2*var.table$`P²AR`, col=method_colors[3], lty=3, lwd=1.2)

lines(df.time.wide$n, df.time.wide$`P²WL2M` + 2*var.table$`P²WL2M`, col=method_colors[4], lty=4, lwd=1.2)
lines(df.time.wide$n, df.time.wide$`P²WL2M` - 2*var.table$`P²WL2M`, col=method_colors[4], lty=4, lwd=1.2)


lines(df.time.wide$n, df.time.wide$`P²WSVM` + 2*var.table$`P²WSVM`, col=method_colors[5], lty=5, lwd=1.2)
lines(df.time.wide$n, df.time.wide$`P²WSVM` - 2*var.table$`P²WSVM`, col=method_colors[5], lty=5, lwd=1.2)

lines(df.time.wide$n, df.time.wide$`P²QR` + 2*var.table$`P²QR`, col=method_colors[6], lty=6, lwd=1.2)
lines(df.time.wide$n, df.time.wide$`P²QR` - 2*var.table$`P²QR`, col=method_colors[6], lty=6, lwd=1.2)

grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 1)      # Grid line width
legend("topleft", legend=c(expression("P"^2~"WL2M"),  expression("P"^2~"AR"), expression("P"^2~"WSVM"), expression("P"^2~"QR"), expression("P"^2~"WLR"), expression("P"^2~"WLSM")),
       col=c("#e7298a", "#7570b3", "#66a61e", "#e6ab02", "#1b9e77", "#d95f02"), lty=c(4,3,5,6,1,2), cex=1.3, lwd=2, pch=c(4,3,5,6,1,2))

