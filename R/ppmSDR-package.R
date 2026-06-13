#' ppmSDR: Penalized Principal Machines for Sufficient Dimension Reduction
#'
#' A unified, computation-friendly framework for penalized principal machines
#' (P2M), a class of sparse sufficient dimension reduction (SDR) estimators for
#' regression and binary classification. All estimators are fitted by a single
#' group coordinate descent (GCD) algorithm that accommodates differentiable
#' and non-differentiable losses through quadratic approximation, together with
#' the group LASSO, SCAD and MCP penalties.
#'
#' The main entry point is [ppm()], which dispatches to the loss-specific
#' solver selected by the `loss` argument and returns an object of class
#' `"ppm"`. See [print.ppm()] and [summary.ppm()] for inspection methods.
#'
#' @references
#' Li, B., Artemiou, A. and Li, L. (2011) Principal support vector machines for
#' linear and nonlinear sufficient dimension reduction. *The Annals of
#' Statistics*, 39(6), 3182--3210. \doi{10.1214/11-AOS932}
#'
#' Shin, S. J. and Artemiou, A. (2017) Penalized principal logistic regression
#' for sparse sufficient dimension reduction. *Computational Statistics & Data
#' Analysis*, 111, 48--58. \doi{10.1016/j.csda.2016.12.003}
#'
#' Breheny, P. and Huang, J. (2015) Group descent algorithms for nonconvex
#' penalized linear and logistic regression models with grouped predictors.
#' *Statistics and Computing*, 25, 173--187. \doi{10.1007/s11222-013-9424-2}
#'
#' @keywords internal
#' @aliases ppmSDR-package
#' @importFrom stats cov quantile
#' @importFrom Matrix bdiag
#' @importFrom energy dcor
#' @importFrom stats sd
#' @importFrom grpreg grpreg
"_PACKAGE"
