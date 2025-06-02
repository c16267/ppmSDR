# ppmSDR: Penalized Principal Machine for Sparse Sufficient Dimension Reduction
A unified and computationally efficient R package for sparse sufficient dimension reduction (SDR) using Penalized Principal Machines (P$^2$M) and Group Coordinate Descent (GCD).

## Overview

The **ppmSDR** package provides a unified interface and efficient algorithms for sparse sufficient dimension reduction (SDR) in regression and classification. It implements the **Penalized Principal Machine (P$^2$M)** family, generalizing principal support vector machine (PSVM) by allowing a wide range of convex loss functions and modern sparsity-inducing penalties.  
Efficient computation is achieved via the Group Coordinate Descent (GCD) and MM-GCD algorithms, making the package scalable to large and high-dimensional data.

### Key Features

- Unified interface for various penalized principal machine estimators (P$^2$M), supporting both regression and classification
- Implements state-of-the-art sparse SDR methods:
    - P$^2$LSM / P$^2$WLSM (Least Squares SVM, weighted)
    - P$^2$LR / P$^2$WLR (Logistic, weighted)
    - P$^2$L2M / P$^2$WL2M (L2-Hinge, weighted)
    - P$^2$SVM / P$^2$WSVM (Hinge, weighted)
    - P$^2$QR (Quantile regression)
    - P$^2$AR (Asymmetric least squares)
- Group SCAD, MCP, and Lasso penalties for row-wise sparsity and variable selection
- Fast, scalable optimization (GCD, MM-GCD)
- Cross-validation functions for penalty parameter tuning
- Functions for simulation and real data analysis (Boston Housing, Wisconsin Breast Cancer)
- Detailed documentation and examples

---

## Installation

```r
# Install from GitHub (requires devtools)
devtools::install_github("c16267/ppmSDR")
```
