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


\section*{Repository Structure}

\begin{itemize}
    \item \textbf{R/} : Core functions for penalized PM estimation and tuning
    \begin{itemize}
        \item \texttt{fn\_pplssvm.R} : Penalized principal least squares SVM (\texttt{pplssvm})
        \item \texttt{fn\_ppals.R} : Penalized principal asymmetric least squares (\texttt{ppalssvm})
        \item \texttt{fn\_ppl2svm.R} : Penalized principal L2-hinge SVM (\texttt{ppl2svm})
        \item \texttt{fn\_pplr.R} : Penalized principal logistic regression (\texttt{pplr})
        \item \texttt{fn\_wpplr.R} : Weighted penalized principal logistic regression (\texttt{wpplr})
        \item \texttt{fn\_spsdr.R} : Unified wrapper for penalized PMs (\texttt{spsvmSDR})
        \item \texttt{fn\_test.R} : Example usage of \texttt{spsvmSDR}
        \item \texttt{fn\_tune\_*.R} : Cross-validation for optimal lambda selection
        \item \texttt{fn\_minor\_pPSDR.R} : Auxiliary functions (thresholding, etc.)
        \item \texttt{fn\_penalized\_logit\_dc.R}, \texttt{fn\_sparse\_SIR.R}, \texttt{fn\_tune\_sparse\_SIR.R} : Other SDR competitors
    \end{itemize}
    \item \textbf{data/} : Example datasets (Boston Housing, Breast Cancer)
    \item \textbf{simulation/} : Scripts to reproduce simulation studies
    \begin{itemize}
        \item \texttt{fn\_simulation\_continuous.R}
        \item \texttt{fn\_simulation\_binary.R}
        \item \texttt{fn\_simulation\_time\_n.R}
    \end{itemize}
\end{itemize}

\section*{Main Functions}

\begin{table}[ht]
    \centering
    \begin{tabular}{lll}
        \textbf{Function} & \textbf{Description} & \textbf{Penalty Options} \\
        \hline
        \texttt{pplssvm} & Penalized principal least squares SVM (P$^2$LSM) & SCAD, Lasso, MCP \\
        \texttt{ppalssvm} & Penalized principal asymmetric least squares (P$^2$AR) & SCAD, Lasso, MCP \\
        \texttt{ppl2svm} & Penalized principal L2-hinge SVM (P$^2$L2M) & SCAD, Lasso, MCP \\
        \texttt{pplr} & Penalized principal logistic regression (P$^2$LR) & SCAD, Lasso, MCP \\
        \texttt{wpplr} & Weighted penalized principal logistic regression (P$^2$WLR) & SCAD, Lasso, MCP \\
        \texttt{ppsvm} & Penalized principal SVM (P$^2$SVM, MM-GCD) & SCAD, Lasso, MCP \\
        \texttt{ppqr} & Penalized principal quantile regression (P$^2$QR) & SCAD, Lasso, MCP \\
        \texttt{ppwlssvm} & Penalized principal weighted least squares SVM (P$^2$WLSM) & SCAD, Lasso, MCP \\
        \texttt{ppwlr} & Penalized principal weighted logistic regression (P$^2$WLR) & SCAD, Lasso, MCP \\
        \texttt{ppwl2svm} & Penalized principal weighted L2-hinge SVM (P$^2$WL2M) & SCAD, Lasso, MCP \\
        \texttt{ppwsvm} & Penalized principal weighted SVM (P$^2$WSVM, MM-GCD) & SCAD, Lasso, MCP \\
    \end{tabular}
    \caption{Summary of main functions in \texttt{ppmSDR} and their penalty options.}
\end{table}

\section*{Unified Wrapper}

\begin{itemize}
    \item \texttt{ppm}: A unified wrapper function to fit any penalized PM estimator with a single interface. Selects loss, penalty, and method automatically via arguments.
    \item Tuning functions (\texttt{tune\_pplssvm}, \texttt{tune\_ppalssvm}, etc.): Cross-validation for sparsity parameter (\texttt{lambda}) selection.
\end{itemize}

\section*{Example Usage}

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



