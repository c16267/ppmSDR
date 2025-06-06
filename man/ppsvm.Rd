% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fn_ppsvm.R
\name{ppsvm}
\alias{ppsvm}
\title{Penalized Principal Support Vector Machine (P\eqn{^2}SVM, MM-GCD-based) for Sparse Sufficient Dimension Reduction}
\usage{
ppsvm(
  x,
  y,
  H = 10,
  C = 1,
  lambda,
  gamma = 3.7,
  penalty = "grSCAD",
  max.iter = 100
)
}
\arguments{
\item{x}{Numeric predictor matrix (n x p), where n is the sample size and p is the number of predictors.}

\item{y}{Numeric response vector of length n.}

\item{H}{The number of quantile slices (default: 10).}

\item{C}{Regularization parameter (default: 1).}

\item{lambda}{Penalty parameter for sparsity.}

\item{gamma}{Regularization parameter for SCAD/MCP penalty (default: 3.7).}

\item{penalty}{Penalty type: \code{"grSCAD"} (default), \code{"grLasso"}, or \code{"grMCP"}.}

\item{max.iter}{Maximum number of iterations for the algorithm (default: 100).}
}
\value{
A list containing:
\item{Mn}{Estimated central subspace matrix.}
\item{evalues}{Eigenvalues of the estimated matrix.}
\item{evectors}{Eigenvectors (columns) corresponding to the estimated sufficient directions.}
\item{x}{The original predictor matrix.}
}
\description{
This function implements the penalized principal support vector machine (P\eqn{^2}SVM) approach for sparse sufficient dimension reduction (SDR), using the MM-GCD algorithm for efficient computation.
}
\details{
This function estimates a sparse basis for the central subspace by solving a penalized principal support vector machine objective with an MM-GCD optimization algorithm, enabling scalable and accurate sparse SDR.
}
\examples{
\dontrun{
  set.seed(1)
  n <- 300; p <- 10
  B <- matrix(0, p, 2)
  B[1,1] <- B[2,2] <- 1
  x <- MASS::mvrnorm(n, rep(0, p), diag(1, p))
  y <- (x \%*\% B[,1] / (0.5 + (x \%*\% B[,2] + 1)^2)) + 0.2 * rnorm(n)
  fit <- ppsvm(x, y, H = 10, C = 100, lambda = 0.0001, gamma = 3.7,
         penalty = "grSCAD", max.iter = 100)
  fit$evectors[, 1:2]
}

}
