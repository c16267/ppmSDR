################################################
# Wrapper function for Penalize Principal Machine 
################################################
setwd("/Users/shin/Desktop/Penalize PM")
source("fn_minor_Sparse_SDR.R")
source("fn_ppasls.R")
source("fn_ppl2svm.R")
source("fn_pplr.R")
source("fn_pplssvm.R")
source("fn_ppqr.R")
source("fn_ppsvm.R")
source("fn_ppwlssvm.R")
source("fn_ppwlr.R")
source("fn_ppwl2svm.R")
source("fn_ppwsvm.R")


ppm <- function(x, y, H = 10, C = 1, loss = "lssvm", penalty = "grSCAD",
                lambda = 0.001, gamma = 3.7, max.iter = 100, ...) {
  
  required_pkgs <- c("MASS", "Matrix", "expm", "grpreg")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Check for required arguments
  if (missing(x) || missing(y)) {
    stop("Both 'x' and 'y' must be provided.")
  }
  
  # Define loss_map: maps loss names to corresponding functions
  loss_map <- list(
    lssvm  = pplssvm,
    wlssvm = ppwlssvm,
    svm    = ppsvm,
    wsvm   = ppwsvm,
    l2svm  = ppl2svm,
    wl2svm = ppwl2svm,
    logit  = pplr,
    wlogit = ppwlr,
    asls   = ppasls,
    qr     = ppqr
  )
  loss <- tolower(loss)
  
  if (!(loss %in% names(loss_map))) {
    stop("Unknown loss function. Choose from: ", paste(names(loss_map), collapse = ", "))
  }
  
  y_vec <- as.vector(y)
  y_unique <- unique(y_vec)
  
  if (grepl("^w", loss)) {
    if (length(y_unique) != 2 || !all(sort(y_unique) == c(-1, 1))) {
      stop("For weighted methods (loss starting with 'w'), the response 'y' must be strictly binary with values +1 and -1.")
    }
  } else {
    if ((length(y_unique) == 2) &&
        (is.numeric(y_vec) || is.logical(y_vec) || is.integer(y_vec) || is.factor(y_vec))) {
      message("Note: For continuous loss functions, you provided a binary response y (two unique values: ",
              paste(y_unique, collapse=", "), "). Please check if this is intended.")
    }
  }
  
  call_args <- c(
    list(
      x = x, y = y, H = H, C = C,
      penalty = penalty, lambda = lambda,
      gamma = gamma, max.iter = max.iter
    ),
    list(...)
  )
  
  do.call(loss_map[[loss]], call_args)
}


##
#example
##
set.seed(1)
n <- 1000
p   <- 10
B <- matrix(0, p, 2)
B[1,1] <- B[2,2] <- 1
x <- mvrnorm(n, rep(0, p), diag(1,p))
y <- (x %*% B[,1]/(0.5 + (x %*% B[,2] + 1)^2)) + 0.2*rnorm(n, 0, 1)
y.binary <- sign(y)

#############################
# -- SDR in Regression  -- #
#############################

# -- PLSSVM vs. P2LSM
result <- ppm(x, y, H=10, C=1, loss="lssvm", penalty = 'grSCAD', lambda=0.1)
round(result$evectors[,1:2], 5)

result0 <- psdr(x, y, h=10, loss="lssvm")
round(result0$evectors[,1:2], 5)


# -- PSVM vs. P2SVM
result11 <- ppm(x, y, H=10, C=1, loss="svm", penalty = 'grSCAD', lambda=0.0001)
round(result11$evectors[,1:2], 5)

result10 <- psdr(x, y, h=10, loss="svm")
round(result10$evectors[,1:2], 5)

######################################
# -- SDR in Binary Classification -- #
######################################

# -- PWLSSVM vs. P2WLSM
result3 <- ppm(x, y.binary, H=10, C=1, loss="wlssvm", penalty = 'grSCAD', lambda=0.002)
round(result3$evectors[,1:2], 5)

result30 <- psdr(x, y.binary, h=10, loss="wlssvm")
round(result30$evectors[,1:2], 5)


# -- PWLR vs. P2WLR
result4 <- ppm(x, y.binary, H=10, C=15, loss="wlogit", penalty = 'grSCAD', lambda=0.02)
round(result4$evectors[,1:2], 5)

result5 <- psdr(x, y.binary, h=10, loss="wlogit")
round(result5$evectors[,1:2], 5)





