# ppmSDR: Penalized Principal Machine for Sparse Sufficient Dimension Reduction
A unified and computationally efficient R package for sparse sufficient dimension reduction (SDR) using Penalized Principal Machines (P²M) and Group Coordinate Descent (GCD).

## Overview

The **ppmSDR** package provides a unified interface and efficient algorithms for sparse sufficient dimension reduction (SDR) in regression and classification. It implements the **Penalized Principal Machine (P²M)** family, generalizing principal support vector machine (PSVM) by allowing a wide range of convex loss functions and modern sparsity-inducing penalties.  
Efficient computation is achieved via the Group Coordinate Descent (GCD) and MM-GCD algorithms, making the package scalable to large and high-dimensional data.

### Key Features
Unified interface for various penalized principal machine estimators (P²M), supporting both regression and classification

Implements state-of-the-art sparse SDR methods:
P²LSM / P²WLSM (Least Squares SVM, weighted)
P²LR / P²WLR (Logistic, weighted)
P²L2M / P²WL2M (L2-Hinge, weighted)
P²SVM / P²WSVM (Hinge, weighted)
P²QR (Quantile regression)
P²AR (Asymmetric least squares)

Group SCAD, MCP, and Lasso penalties for row-wise sparsity and variable selection
Fast, scalable optimization (GCD, MM-GCD)
Cross-validation functions for penalty parameter tuning
Functions for simulation and real data analysis (Boston Housing, Wisconsin Breast Cancer)

## Installation

```r
# Install from GitHub (requires devtools)
devtools::install_github("c16267/ppmSDR")
```


### Repository Structure

1. R/ : Core functions for penalized PM estimation and tuning
 - fn_pplssvm.R : Penalized principal least squares SVM (pplssvm)
 - fn_ppwlssvm.R : Penalized principal weighted least squares SVM (ppwlssvm)
 - fn_ppasls.R : Penalized principal asymmetric least squares (ppasls)
 - fn_ppqr.R : Penalized principal quantile regression (ppqr)
 - fn_ppl2svm.R : Penalized principal L2-hinge SVM (ppl2svm)
 - fn_ppwl2svm.R : Penalized principal weighted L2-hinge SVM (ppwl2svm)
 - fn_pplr.R : Penalized principal logistic regression (pplr)
 - fn_ppwlr.R : Penalized principal weighted logistic regression (ppwlr)
 - fn_ppsvm.R : Penalized principal support vector machine (ppsvm)
 - fn_ppwsvm.R : Penalized principal weighted support vector machine (ppwsvm)
 - fn_minor_pPSDR.R : Auxiliary functions (thresholding, etc.)

2. data/ : Example datasets
Boston Housing
Breast Cancer

3. simulation/ : Scripts to reproduce simulation studies
fn_simulation_continuous.R
fn_simulation_binary.R
fn_simulation_time_n.R

## Main Functions

| Function   | Description                                              | 
| ---------- | -------------------------------------------------------- |
| `pplssvm`  | Penalized principal least squares SVM (P²LSM)            |
| `ppasls`   | Penalized principal asymmetric least squares (P²AR)      |
| `ppl2svm`  | Penalized principal L2-hinge SVM (P²L2M)                 |
| `pplr`     | Penalized principal logistic regression (P²LR)           |
| `ppsvm`    | Penalized principal SVM (P²SVM, MM-GCD)                  |
| `ppqr`     | Penalized principal quantile regression (P²QR)           |
| `ppwlssvm` | Penalized principal weighted least squares SVM (P²WLSM)  |
| `ppwlr`    | Penalized principal weighted logistic regression (P²WLR) |
| `ppwl2svm` | Penalized principal weighted L2-hinge SVM (P²WL2M)       |
| `ppwsvm`   | Penalized principal weighted SVM (P²WSVM, MM-GCD)        |

### Unified Wrapper
ppm: A unified wrapper function to fit any penalized PM estimator with a single interface. Selects loss, penalty, and method automatically via arguments.

## Example Usage

```r
library(ppmSDR)

# Generate data
set.seed(1)
n <- 200; p <- 10
B <- matrix(0, p, 2); B[1,1] <- B[2,2] <- 1
x <- MASS::mvrnorm(n, rep(0, p), diag(1, p))
y <- (x %*% B[,1]/(0.5 + (x %*% B[,2] + 1)^2)) + 0.2*rnorm(n)

# Fit penalized principal least squares SVM (P2LSM)
fit <- pplssvm(x, y, H = 10, C = 1, lambda = 0.01, gamma = 3.7, penalty = "grSCAD", max.iter = 100)
fit$evectors[,1:2]

# Unified wrapper (any method)
fit2 <- ppm(x, y, H = 10, C = 1, loss = "lssvm", penalty = "grSCAD", lambda = 0.01)
fit2$evectors[,1:2]
```

## References

- Artemiou, A. and Dong, Y. (2016). Sufficient dimension reduction via principal lq support vector machine, *Electronic Journal of Statistics*, 10: 783–805.
- Artemiou, A., Dong, Y. and Shin, S. J. (2021). Real-time sufficient dimension reduction through principal least squares support vector machines, *Pattern Recognition*, 112: 107768.
- Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized likelihood and its oracle properties, *Journal of the American Statistical Association*, 96: 1348–1360.
- Hunter, D. R. and Lange, K. (2004). A tutorial on MM algorithms, *The American Statistician*, 58(1): 30–37.
- Jang, H. J., Shin, S. J. and Artemiou, A. (2023). Principal weighted least square support vector machine: An online dimension-reduction tool for binary classification, *Computational Statistics & Data Analysis*, 187: 107818.
- Kim, B. and Shin, S. J. (2019). Principal weighted logistic regression for sufficient dimension reduction in binary classification, *Journal of the Korean Statistical Society*, 48(2): 194–206.
- Li, B., Artemiou, A. and Li, L. (2011). Principal support vector machines for linear and nonlinear sufficient dimension reduction, *Annals of Statistics*, 39(6): 3182–3210.
- Shin, J. and Shin, S. J. (2024). A concise overview of principal support vector machines and its generalization, *Communications for Statistical Applications and Methods*, 31(2): 235–246.
- Shin, J., Shin, S. J. and Artemiou, A. (2024). The R package psvmsdr: A unified algorithm for sufficient dimension reduction via principal machines, *arXiv preprint arXiv:2409.01547*.
- Shin, S. J. and Artemiou, A. (2017). Penalized principal logistic regression for sparse sufficient dimension reduction, *Computational Statistics & Data Analysis*, 111: 48–58.



