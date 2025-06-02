library(testthat)
library(ppmSDR)

test_that("ppm returns a list with expected components", {
  set.seed(1)
  n <- 300; p <- 10
  B <- matrix(0, p, 2); B[1,1] <- B[2,2] <- 1
  x <- MASS::mvrnorm(n, rep(0, p), diag(1, p))
  y <- (x %*% B[,1]/(0.5 + (x %*% B[,2] + 1)^2)) + 0.2*rnorm(n)

  res <- ppm(x, y, H = 10, C = 1, loss = "lssvm", lambda = 0.03)
  expect_type(res, "list")
  expect_true("evectors" %in% names(res))
  expect_true("evalues" %in% names(res))
})

