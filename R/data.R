#' Boston Housing Data
#'
#' The Boston housing data of Harrison and Rubinfeld (1978), as distributed in
#' \pkg{MASS}. Following the real-data analysis in the package's accompanying
#' article, the binary Charles-River dummy variable (`chas`) is removed,
#' leaving twelve continuous predictors and the response `medv`. Values are on
#' their natural scale; standardize them before fitting (see Examples).
#'
#' @format A data frame with 506 rows and 13 variables:
#' \describe{
#'   \item{crim}{per-capita crime rate by town.}
#'   \item{zn}{proportion of residential land zoned for lots over 25,000 sq ft.}
#'   \item{indus}{proportion of non-retail business acres per town.}
#'   \item{nox}{nitrogen-oxides concentration (parts per 10 million).}
#'   \item{rm}{average number of rooms per dwelling.}
#'   \item{age}{proportion of owner-occupied units built prior to 1940.}
#'   \item{dis}{weighted mean of distances to five Boston employment centres.}
#'   \item{rad}{index of accessibility to radial highways.}
#'   \item{tax}{full-value property-tax rate per $10,000.}
#'   \item{ptratio}{pupil-teacher ratio by town.}
#'   \item{b}{\eqn{1000 (Bk - 0.63)^2}, where \eqn{Bk} is the proportion of the
#'     Black population by town. See the note below.}
#'   \item{lstat}{lower-status percentage of the population.}
#'   \item{medv}{median value of owner-occupied homes in $1000s (the response).}
#' }
#'
#' @note The variable `b` embeds a racist premise from the original 1978 study
#'   and is retained only to reproduce the published analysis. It should not be
#'   used uncritically; see the discussion in the documentation of
#'   `MASS::Boston`.
#'
#' @source Harrison, D. and Rubinfeld, D. L. (1978) Hedonic prices and the
#'   demand for clean air. *Journal of Environmental Economics and Management*,
#'   5, 81--102. Distributed in the \pkg{MASS} package.
#'
#' @examples
#' data(boston)
#' x <- scale(as.matrix(boston[, setdiff(names(boston), "medv")]))
#' y <- as.numeric(scale(boston$medv))
#' fit <- ppm(x, y, loss = "lssvm", penalty = "grSCAD", lambda = 8e-3)
#' summary(fit, d = 2)
"boston"


#' Wisconsin Diagnostic Breast Cancer (WDBC) Data
#'
#' Diagnostic measurements for 569 breast masses from the Wisconsin Diagnostic
#' Breast Cancer study (Street, Wolberg and Mangasarian; 1993). Each record has
#' thirty real-valued features computed from a digitized image of a fine-needle
#' aspirate, namely the mean, standard error (`_se`) and worst (largest) value
#' of ten cell-nucleus characteristics, together with a benign/malignant
#' diagnosis. Features are on their natural scale; standardize them before
#' fitting (see Examples).
#'
#' @format A data frame with 569 rows and 31 variables: a factor `diagnosis`
#'   with levels `"B"` (benign, 357 cases) and `"M"` (malignant, 212 cases),
#'   followed by thirty numeric features. The features are the `_mean`, `_se`
#'   and `_worst` summaries of: `radius`, `texture`, `perimeter`, `area`,
#'   `smoothness`, `compactness`, `concavity`, `concave_points`, `symmetry`
#'   and `fractal_dimension` (for example `radius_mean`, `radius_se`,
#'   `radius_worst`).
#'
#' @details For the weighted-loss principal machines (`"wsvm"`, `"wlssvm"`,
#'   `"wlogit"`, `"wl2svm"`) the response is coded as malignant `= +1` and
#'   benign `= -1`, as shown in the Examples.
#'
#' @source Street, W. N., Wolberg, W. H. and Mangasarian, O. L. (1993) Nuclear
#'   feature extraction for breast tumor diagnosis. *Biomedical Image
#'   Processing and Biomedical Visualization*, 1905, 861--870. UCI Machine
#'   Learning Repository, Breast Cancer Wisconsin (Diagnostic) Data Set.
#'
#' @examples
#' data(wdbc)
#' x <- scale(as.matrix(wdbc[, -1]))
#' y <- ifelse(wdbc$diagnosis == "M", 1, -1)
#' fit <- ppm(x, y, loss = "wl2svm", penalty = "grSCAD", lambda = 0.3)
#' summary(fit, d = 2)
#' B      <- fit$evectors[, 1:2]
#' scores <- x %*% B
#' plot(scores[, 1], scores[, 2], col  = ifelse(y == 1, "red", "blue"),
#'  pch  = ifelse(y == 1, 17, 1),
#'  xlab = "1st SDR direction", ylab = "2nd SDR direction")
#' legend("topright", legend = c("malignant (+1)", "benign (-1)"),
#'  col = c("red", "blue"), pch = c(17, 1))
"wdbc"
