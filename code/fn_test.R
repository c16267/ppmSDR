rm(list=ls())

setwd("/Users/shin/Dropbox/shared_JMShin/(ver. 240524) [JCGS] Computationally efficient algorithm for a sparse dimension reduction using a squared loss/penalized-PSDR/R/")
files.sources = list.files()
sapply(files.sources, source)

library(MASS)

n <- 1000 # sample size
p   <- 5
max.iter <- maxiter <- 100
H <- 10
h <- 1.0e-5; eps <- 1.0e-5
C <- 1
delta <- 0.1

# true cs
B <- matrix(0, p, 2)
B[1,1] <- B[2,2] <- 1; d=q = 2
Sigma <- diag(1,p)
x <- mvrnorm(n, rep(0, p), Sigma)
eps <- rnorm(n, 0, 1)
x1 <- x %*% B[,1]
x2 <- x %*% B[,2]
model=1
# responses
{
  y <- if (model == 1) (x1/(0.5 + (x2 + 1)^2))+0.2*eps
  else if (model == 2) (x1+0.5)*(x2-0.5)^2+0.2*eps
  else if (model == 3) 0.5*x1^3 + (2*x2+3)*eps #asymmetric Dong #3
}   
y.binary <- sign(y)


#method <- c("pplsvm", "ppals", "pplr", "wpplr", "ppl2svm")
###########
##1.P2LSM##
###########
rslt_pplssvm <- spsvmSDR(x, y, H, C, method="pplsvm", lambda=0.01, gamma=3.7, penalty="grSCAD", max.iter=100)
round(rslt_pplssvm$vectors, 3)

###########
##2.P2AR ##
###########
rslt_ppals <- spsvmSDR(x, y, H=5, C, method="ppals", lambda=0.05, gamma=3.7, penalty="grSCAD", max.iter=100)
round(rslt_ppals$vectors, 3)

#for binary
rslt_ppals_binary <- spsvmSDR(x, y=y.binary, H, C=5, method="ppals", lambda=0.041, gamma=3.7, penalty="grSCAD", max.iter=100)
round(rslt_ppals_binary$vectors, 3)

###########
##3.P2LR ##
###########
rslt_pplr <- spsvmSDR(x, y, H, C, method="pplr", lambda=0.0002, gamma=3.7, penalty="grSCAD", max.iter=100)
round(rslt_pplr$vectors, 3)

###########
##4.P2L2M #
###########
rslt_ppl2svm <- spsvmSDR(x, y, H, C, method="ppl2svm", lambda=0.05, gamma=3.7, penalty="grSCAD", max.iter=100)
round(rslt_ppl2svm$vectors, 3)

###########
##5.P2WLR #
###########
rslt_wpplr <- spsvmSDR(x, y=y.binary, H, C=5, method="wpplr", lambda=0.002, gamma=3.7, penalty="grSCAD", max.iter=100)
round(rslt_wpplr$vectors, 3)
