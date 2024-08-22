rm(list=ls())

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
source('fn_minor_pPSDR.R')
source('fn_pplssvm.R')
source('fn_pplr.R')
source('fn_wpplr.R')
source('fn_ppl2svm.R')
source('fn_ppals.R')
source('fn_spsdr.R')


n <- 300 # sample size
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
rslt_pplssvm <- spsvmSDR(x, y, H, C, method="pplsvm", lambda=0.01, gamma=3.7, penalty="grSCAD", max.iter=100)
rslt_pplssvm$vectors

rslt_ppals <- spsvmSDR(x, y, H, C, method="ppals", lambda=0.05, gamma=3.7, penalty="grSCAD", max.iter=100)
rslt_ppals$vectors

rslt_pplr <- spsvmSDR(x, y, H, C, method="pplr", lambda=0.0002, gamma=3.7, penalty="grSCAD", max.iter=100)
rslt_pplr$vectors

rslt_ppl2svm <- spsvmSDR(x, y, H, C, method="ppl2svm", lambda=0.05, gamma=3.7, penalty="grSCAD", max.iter=100)
rslt_ppl2svm$vectors


rslt_wpplr <- spsvmSDR(x, y=y.binary, H, C, method="wpplr", lambda=0.001, gamma=3.7, penalty="grSCAD", max.iter=100)
rslt_wpplr$vectors
