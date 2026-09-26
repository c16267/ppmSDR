# Internal penalized-principal-machine (P2M) estimators.
# These are NOT exported; the user-facing entry point is ppm().
# Each solver returns: M (working matrix), evalues, evectors, x, theta (the
# (p+1) x h coefficient matrix Theta, intercepts in the first row) and the
# convergence diagnostics line.search, iter, status, objective, n.halving.
# Helper operators (thresholding, block-diagonal algebra, the objective
# Q(theta) and the step-halving line search) live in utils-internal.R.
#
# Algorithm classes (Table 1 of the article):
#   GCD           : pplssvm, ppwlssvm            (exact; no outer iteration)
#   iterative GCD : pplr, ppwlr, ppasls, ppl2svm, ppwl2svm
#                   (local quadratic approximation; step-halving line search
#                    on the original objective Q, default line.search = TRUE
#                    except for ppasls, which uses the plain updates unless
#                    line.search = TRUE is requested)
#   MM-GCD        : ppsvm, ppwsvm, ppqr
#                   (quadratic majorizer; descent follows from the MM
#                    argument, so the line search is optional and off by
#                    default, line.search = FALSE)

# Common epilogue: working matrix, eigen-decomposition and diagnostics.
.ppm_finish <- function(theta, h_eff, x, diag) {
  Theta    <- matrix(theta, ncol = h_eff)                   # (p+1) x h
  beta_mat <- Theta[-1L, , drop = FALSE]                    # drop intercepts
  p  <- nrow(beta_mat)
  Mn <- matrix(0, p, p)
  for (k in 1:h_eff) Mn <- Mn + beta_mat[, k, drop = FALSE] %*% t(beta_mat[, k, drop = FALSE])
  eg <- eigen(Mn)
  c(list(M = Mn, evalues = eg$values, evectors = eg$vectors, x = x,
         theta = Theta), diag)
}


# ============================================================================
# 1. ppasls  (P2AR; iterative GCD, plain updates by default)
# ============================================================================
ppasls <- function(x, y, H, C, lambda, gamma = 3.7, penalty = "grSCAD", max.iter = 100,
                   line.search = FALSE, max.halving = 10L) {
  n <- nrow(x)
  p <- ncol(x)

  if (!penalty %in% c("grSCAD", "grLasso", "grMCP")) stop("Invalid penalty")

  bar.x   <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))

  qprob   <- seq(1/H, 1 - 1/H, by = 1/H)
  h_eff   <- length(qprob)
  qy      <- quantile(y, qprob)
  tmp.y   <- rep(y, times = h_eff)

  Sigma.hat       <- cov(x)
  Sigma.hat.star  <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma.hat))

  init.theta <- rep(0, h_eff * (p + 1))
  tol        <- 1e-5
  tau.vec    <- rep(qprob, each = n)

  ## ---- line-search bookkeeping -------------------------------------------
  Q.path <- numeric(max.iter); m.path <- integer(max.iter); status <- "max.iter"
  if (line.search) {
    ctx   <- .ppm_make_ctx("asls", x.tilde, Sigma.hat.star,
                           Y = matrix(y, n, h_eff), tau = qprob,
                           C = C, lambda = lambda, gamma = gamma, penalty = penalty)
    Q.old <- .ppm_objective(init.theta, ctx)
  }

  for (iter in 1:max.iter) {
    old.theta <- init.theta
    Wtheta <- numeric(n * h_eff)
    for (k in 1:h_eff) {
      theta_k <- init.theta[((k - 1) * (p + 1) + 1):(k * (p + 1))]
      Wtheta[((k - 1) * n + 1):(k * n)] <- as.vector(x.tilde %*% theta_k)
    }

    u   <- (1 - tau.vec) * as.integer(tmp.y - Wtheta <= 0)
    u.c <- (tau.vec)     * as.integer(tmp.y - Wtheta >  0)

    y.tilde   <- tmp.y * sqrt(u)
    y.tilde.c <- tmp.y * sqrt(u.c)

    A_blocks <- vector("list", h_eff)
    B_blocks <- vector("list", h_eff)

    for (k in 1:h_eff) {
      idx_k <- ((k - 1) * n + 1):(k * n)
      wTw_k   <- crossprod(x.tilde * sqrt(u[idx_k]))
      wTw_c_k <- crossprod(x.tilde * sqrt(u.c[idx_k]))

      # Each of the two transformed Gram matrices carries n S / (2c), so that
      # their sum carries n S / c as in (11); see Supplement S.2.3.
      A_blocks[[k]] <- (n / (2 * C)) * Sigma.hat.star + wTw_k
      B_blocks[[k]] <- (n / (2 * C)) * Sigma.hat.star + wTw_c_k
    }

    A_sqrt_blocks <- sqrt_sym_blockdiag(A_blocks)
    B_sqrt_blocks <- sqrt_sym_blockdiag(B_blocks)

    rhs   <- numeric(h_eff * (p + 1))
    rhs.c <- numeric(h_eff * (p + 1))
    for (k in 1:h_eff) {
      idx_k <- ((k - 1) * n + 1):(k * n)
      rhs[((k - 1) * (p + 1) + 1):(k * (p + 1))]   <- t(x.tilde * sqrt(u[idx_k])) %*% y.tilde[idx_k]
      rhs.c[((k - 1) * (p + 1) + 1):(k * (p + 1))] <- t(x.tilde * sqrt(u.c[idx_k])) %*% y.tilde.c[idx_k]
    }

    xi   <- solve_blockdiag(A_sqrt_blocks, rhs)
    xi.c <- solve_blockdiag(B_sqrt_blocks, rhs.c)

    pos.resid <- as.vector(xi)   - multiply_blockdiag(A_sqrt_blocks, init.theta)
    neg.resid <- as.vector(xi.c) - multiply_blockdiag(B_sqrt_blocks, init.theta)
    res       <- pos.resid + neg.resid

    group <- rep(1:(p + 1), times = h_eff)
    J     <- p + 1

    for (j in 1:J) {
      ind <- which(group == j)
      G_cols   <- get_var_columns(A_sqrt_blocks, j)
      G_c_cols <- get_var_columns(B_sqrt_blocks, j)

      z_j <- (C / n) * (t(G_cols) %*% pos.resid + t(G_c_cols) %*% neg.resid) + 2 * init.theta[ind]
      z_norm <- norm(z_j, "2")

      if (z_norm < 1e-12) {
        theta_j_new <- rep(0, length(ind))
      } else {
        if (penalty == "grSCAD") {
          theta_j_new <- scad_firm_thresholding(z_norm, lambda1=lambda, lambda2=1, gamma=gamma) * z_j / z_norm
        } else if (penalty == "grLasso") {
          theta_j_new <- soft_thresholding(z_norm, lambda=lambda)/(1 + 2) * z_j / z_norm
        } else if (penalty == "grMCP") {
          theta_j_new <- firm_thresholding(z_norm, lambda1=lambda, lambda2=2, gamma=gamma) * z_j / z_norm
        }
        theta_j_new <- as.vector(theta_j_new)
      }

      pos.resid <- pos.resid - G_cols   %*% (theta_j_new - init.theta[ind])
      neg.resid <- neg.resid - G_c_cols %*% (theta_j_new - init.theta[ind])
      init.theta[ind] <- theta_j_new

      pos.resid <- as.vector(pos.resid)
      neg.resid <- as.vector(neg.resid)
    }

    .ppm_check_finite(init.theta, iter)

    ## ---- step-halving safeguard on the original objective Q ---------------
    if (line.search) {
      ls <- .ppm_step_halving(old.theta, init.theta, Q.old, ctx, max.halving, tol)
      init.theta <- ls$theta; Q.path[iter] <- ls$Q; m.path[iter] <- ls$m
      if (!ls$accepted) {
        status <- if (ls$halted) "line.search.halted" else "converged"; break
      }
      Q.old <- ls$Q
    }

    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))

    if (delta < tol) { status <- "converged"; break }
  }

  .ppm_finish(init.theta, h_eff, x,
              .ppm_diag(line.search, iter, status, Q.path, m.path))
}


# ============================================================================
# 2. ppl2svm  (P2L2M; iterative GCD)
# ============================================================================
ppl2svm <- function(x, y, H, C, lambda, gamma = 3.7, penalty = "grSCAD", max.iter = 100,
                    line.search = TRUE, max.halving = 10L) {
  n <- nrow(x)
  p <- ncol(x)

  if (!penalty %in% c("grSCAD", "grLasso", "grMCP")) stop("Invalid penalty")

  bar.x   <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))

  qprob <- seq(1/H, 1 - 1/H, by = 1/H)
  h_eff <- length(qprob)
  qy    <- quantile(y, qprob)
  U     <- sapply(qy, function(z) 2 * (z > y) - 1)
  tmp.y <- c(U)                                   # stacked y_tilde_{ik}, machine by machine

  Sigma.hat      <- cov(x)
  Sigma.hat.star <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma.hat))

  init.theta <- rep(0, h_eff * (p + 1))
  tol        <- 1e-5

  ## ---- line-search bookkeeping -------------------------------------------
  Q.path <- numeric(max.iter); m.path <- integer(max.iter); status <- "max.iter"
  if (line.search) {
    ctx   <- .ppm_make_ctx("l2", x.tilde, Sigma.hat.star, Y = matrix(tmp.y, n, h_eff),
                           C = C, lambda = lambda, gamma = gamma, penalty = penalty)
    Q.old <- .ppm_objective(init.theta, ctx)
  }

  for (iter in 1:max.iter) {
    old.theta <- init.theta

    # f_{ik}^(t) = theta_k^(t)' x_tilde_i  (unweighted fitted values)
    Wtheta <- numeric(n * h_eff)
    for (k in 1:h_eff) {
      theta_k <- init.theta[((k - 1) * (p + 1) + 1):(k * (p + 1))]
      Wtheta[((k - 1) * n + 1):(k * n)] <- as.vector(x.tilde %*% theta_k)
    }

    # ---- active set: 1 - y_tilde * f > 0 (margin violators) -----------------
    u <- as.integer(1 - tmp.y * Wtheta > 0)
    y.tilde <- tmp.y * sqrt(u)

    A_blocks <- vector("list", h_eff)
    for (k in 1:h_eff) {
      idx_k <- ((k - 1) * n + 1):(k * n)
      wTw_k <- crossprod(x.tilde * sqrt(u[idx_k]))
      A_blocks[[k]] <- 0.5 * Sigma.hat.star + wTw_k
    }

    A_sqrt_blocks <- sqrt_sym_blockdiag(A_blocks)

    rhs <- numeric(h_eff * (p + 1))
    for (k in 1:h_eff) {
      idx_k <- ((k - 1) * n + 1):(k * n)
      rhs[((k - 1) * (p + 1) + 1):(k * (p + 1))] <- t(x.tilde * sqrt(u[idx_k])) %*% y.tilde[idx_k]
    }

    xi <- solve_blockdiag(A_sqrt_blocks, rhs)
    res <- as.vector(xi) - multiply_blockdiag(A_sqrt_blocks, init.theta)

    group <- rep(1:(p + 1), times = h_eff)
    J     <- p + 1

    for (j in 1:J) {
      ind <- which(group == j)
      G_cols <- get_var_columns(A_sqrt_blocks, j)

      z_j <- (C / n) * (t(G_cols) %*% res) + init.theta[ind]
      z_norm <- norm(z_j, "2")

      if (z_norm < 1e-12) {
        theta_j_new <- rep(0, length(ind))
      } else {
        if (penalty == "grSCAD") {
          theta_j_new <- scad_firm_thresholding(z_norm, lambda1 = lambda, lambda2 = 2, gamma = gamma) * z_j / z_norm
        } else if (penalty == "grLasso") {
          theta_j_new <- soft_thresholding(z_norm, lambda = lambda) / (1 + 2) * z_j / z_norm
        } else if (penalty == "grMCP") {
          theta_j_new <- firm_thresholding(z_norm, lambda1 = lambda, lambda2 = 2, gamma = gamma) * z_j / z_norm
        }
        theta_j_new <- as.vector(theta_j_new)
      }

      res <- res - G_cols %*% (theta_j_new - init.theta[ind])
      init.theta[ind] <- theta_j_new
      res <- as.vector(res)
    }

    .ppm_check_finite(init.theta, iter)

    ## ---- step-halving safeguard on the original objective Q ---------------
    if (line.search) {
      ls <- .ppm_step_halving(old.theta, init.theta, Q.old, ctx, max.halving, tol)
      init.theta <- ls$theta; Q.path[iter] <- ls$Q; m.path[iter] <- ls$m
      if (!ls$accepted) {
        status <- if (ls$halted) "line.search.halted" else "converged"; break
      }
      Q.old <- ls$Q
    }

    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))
    if (delta < tol) { status <- "converged"; break }
  }

  .ppm_finish(init.theta, h_eff, x,
              .ppm_diag(line.search, iter, status, Q.path, m.path))
}

# ============================================================================
# 3. pplr  (P2LR; iterative GCD)
# ============================================================================
pplr <- function(x, y, H, C, lambda, gamma = 3.7, penalty = "grSCAD", max.iter = 100, ridge = 1e-10,
                 line.search = TRUE, max.halving = 10L) {
  n <- nrow(x)
  p <- ncol(x)

  if (!penalty %in% c("grSCAD", "grLasso", "grMCP")) stop("Invalid penalty")

  bar.x   <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))

  qprob <- seq(1/H, 1 - 1/H, by = 1/H)
  h_eff <- length(qprob)
  qy    <- quantile(y, qprob)
  U     <- sapply(qy, function(z) 2 * (z > y) - 1)
  u     <- c(U)                                   # stacked y_tilde_{ik}, machine by machine

  Sigma             <- cov(x)
  Sigma.tilde       <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma))
  Sigma.tilde.ridge <- Sigma.tilde + diag(ridge, p + 1)

  init.theta <- rep(0, h_eff * (p + 1))
  tol <- 1e-5

  ## ---- line-search bookkeeping -------------------------------------------
  Q.path <- numeric(max.iter); m.path <- integer(max.iter); status <- "max.iter"
  if (line.search) {
    ctx   <- .ppm_make_ctx("logit", x.tilde, Sigma.tilde, Y = matrix(u, n, h_eff),
                           C = C, lambda = lambda, gamma = gamma, penalty = penalty)
    Q.old <- .ppm_objective(init.theta, ctx)
  }

  for (iter in 1:max.iter) {
    old.theta <- init.theta

    # f_{ik}^(t) = theta_k^(t)' x_tilde_i  (unweighted fitted values)
    Wtheta <- numeric(n * h_eff)
    for (k in 1:h_eff) {
      theta_k <- init.theta[((k - 1) * (p + 1) + 1):(k * (p + 1))]
      idx_k   <- ((k - 1) * n + 1):(k * n)
      Wtheta[idx_k] <- as.vector(x.tilde %*% theta_k)
    }

    # ---- IRLS quantities, recomputed from the current fit ------------------
    q_vec   <- 1 / (1 + exp(u * Wtheta))          # q = 1 / (1 + exp(y f))
    h_vec   <- q_vec * (1 - q_vec)                # H = q (1 - q)
    sqrt_h  <- sqrt(h_vec)

    u.tilde     <- Wtheta + u / (1 - q_vec)       # u = f + y / (1 - q)
    u.tilde_new <- sqrt_h * u.tilde

    A_blocks <- vector("list", h_eff)
    for (k in 1:h_eff) {
      idx_k       <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * sqrt_h[idx_k]
      A_blocks[[k]] <- Sigma.tilde.ridge + crossprod(Wk_weighted)
    }

    A_sqrt_blocks <- sqrt_sym_blockdiag(A_blocks)

    rhs <- numeric(h_eff * (p + 1))
    for (k in 1:h_eff) {
      idx_k       <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * sqrt_h[idx_k]
      rhs[((k - 1) * (p + 1) + 1):(k * (p + 1))] <- t(Wk_weighted) %*% u.tilde_new[idx_k]
    }
    big_Y <- solve_blockdiag(A_blocks, rhs)
    res <- big_Y - multiply_blockdiag(A_sqrt_blocks, init.theta)

    group <- rep(1:(p + 1), times = h_eff)
    J     <- p + 1

    for (j in 1:J) {
      ind <- which(group == j)
      G_cols <- get_var_columns(A_sqrt_blocks, j)

      z_j <- (C / n) * (t(G_cols) %*% res) + init.theta[ind]
      z_norm <- norm(z_j, "2")

      if (z_norm < 1e-12) {
        theta_j_new <- rep(0, length(ind))
      } else {
        if (penalty == "grSCAD") {
          theta_j_new <- scad_firm_thresholding(z_norm, lambda1 = lambda, lambda2 = 2, gamma = gamma) * z_j / z_norm
        } else if (penalty == "grLasso") {
          theta_j_new <- soft_thresholding(z_norm, lambda = lambda) / (1 + 2) * z_j / z_norm
        } else if (penalty == "grMCP") {
          theta_j_new <- firm_thresholding(z_norm, lambda1 = lambda, lambda2 = 2, gamma = gamma) * z_j / z_norm
        }
        theta_j_new <- as.vector(theta_j_new)
      }

      res <- res - G_cols %*% (theta_j_new - init.theta[ind])
      init.theta[ind] <- theta_j_new
      res <- as.vector(res)
    }

    .ppm_check_finite(init.theta, iter)

    ## ---- step-halving safeguard on the original objective Q ---------------
    if (line.search) {
      ls <- .ppm_step_halving(old.theta, init.theta, Q.old, ctx, max.halving, tol)
      init.theta <- ls$theta; Q.path[iter] <- ls$Q; m.path[iter] <- ls$m
      if (!ls$accepted) {
        status <- if (ls$halted) "line.search.halted" else "converged"; break
      }
      Q.old <- ls$Q
    }

    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))
    if (delta < tol) { status <- "converged"; break }
  }

  .ppm_finish(init.theta, h_eff, x,
              .ppm_diag(line.search, iter, status, Q.path, m.path))
}


# ============================================================================
# 4. ppsvm  (P2SVM; MM-GCD, line search optional)
# ============================================================================
ppsvm <- function(x, y, H = 10, C = 1, lambda, gamma = 3.7, penalty = "grSCAD",
                  max.iter = 100, ridge = 1e-10, eps = 1e-6,
                  line.search = FALSE, max.halving = 10L) {
  n <- nrow(x)
  p <- ncol(x)

  if (!penalty %in% c("grSCAD", "grLasso", "grMCP")) stop("Invalid penalty")
  if (eps <= 0) stop("eps must be positive")

  bar.x   <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))

  qprob <- seq(1/H, 1 - 1/H, by = 1/H)
  h_eff <- length(qprob)
  qy    <- quantile(y, qprob)
  u     <- sapply(qy, function(z) 2 * (z > y) - 1)
  u     <- c(u)                                   # stacked y_tilde_{ik}, machine by machine

  Sigma       <- cov(x)
  Sigma.tilde <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma))

  init.theta <- rep(0, h_eff * (p + 1))
  tol <- 1e-5

  ## ---- line-search bookkeeping -------------------------------------------
  Q.path <- numeric(max.iter); m.path <- integer(max.iter); status <- "max.iter"
  if (line.search) {
    ctx   <- .ppm_make_ctx("hinge", x.tilde, Sigma.tilde, Y = matrix(u, n, h_eff),
                           C = C, lambda = lambda, gamma = gamma, penalty = penalty)
    Q.old <- .ppm_objective(init.theta, ctx)
  }

  for (iter in 1:max.iter) {
    old.theta <- init.theta

    # f_{ik}^(t) = theta_k^(t)' x_tilde_i  (unweighted, no division by n)
    Wtheta <- numeric(n * h_eff)
    for (k in 1:h_eff) {
      theta_k <- init.theta[((k - 1) * (p + 1) + 1):(k * (p + 1))]
      idx_k   <- ((k - 1) * n + 1):(k * n)
      Wtheta[idx_k] <- as.vector(x.tilde %*% theta_k)
    }

    # ---- Safeguarded MM weights (recomputed each iteration) ----------------
    a_t          <- pmax(abs(1 - u * Wtheta), eps)   # a = max(|1 - y f|, eps)
    omega_t      <- 1 / (4 * a_t)                    # Omega = 1 / (4 a)
    sqrt_omega_t <- sqrt(omega_t)

    u.tilde     <- (1 + a_t) * u                     # working response (1 + a) y
    u.tilde_new <- sqrt_omega_t * u.tilde

    A_blocks <- vector("list", h_eff)
    for (k in 1:h_eff) {
      idx_k       <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * sqrt_omega_t[idx_k]
      A_blocks[[k]] <- Sigma.tilde + crossprod(Wk_weighted)
    }

    A_sqrt_blocks <- sqrt_sym_blockdiag(A_blocks, ridge = ridge)

    rhs <- numeric(h_eff * (p + 1))
    for (k in 1:h_eff) {
      idx_k       <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * sqrt_omega_t[idx_k]
      rhs[((k - 1) * (p + 1) + 1):(k * (p + 1))] <- t(Wk_weighted) %*% u.tilde_new[idx_k]
    }
    big_Y <- solve_blockdiag(A_blocks, rhs, ridge = ridge)
    res <- big_Y - multiply_blockdiag(A_sqrt_blocks, init.theta)

    group <- rep(1:(p + 1), times = h_eff)
    J     <- p + 1

    for (j in 1:J) {
      ind <- which(group == j)
      G_cols <- get_var_columns(A_sqrt_blocks, j)

      z_j <- (C / n) * (t(G_cols) %*% res) + init.theta[ind]
      z_norm <- norm(z_j, "2")

      if (z_norm < 1e-12) {
        theta_j_new <- rep(0, length(ind))
      } else {
        if (penalty == "grSCAD") {
          theta_j_new <- scad_firm_thresholding(z_norm, lambda1 = lambda, lambda2 = 2, gamma = gamma) * z_j / z_norm
        } else if (penalty == "grLasso") {
          theta_j_new <- soft_thresholding(z_norm, lambda = lambda) / (1 + 2) * z_j / z_norm
        } else if (penalty == "grMCP") {
          theta_j_new <- firm_thresholding(z_norm, lambda1 = lambda, lambda2 = 2, gamma = gamma) * z_j / z_norm
        }
        theta_j_new <- as.vector(theta_j_new)
      }

      res <- res - G_cols %*% (theta_j_new - init.theta[ind])
      init.theta[ind] <- theta_j_new
      res <- as.vector(res)
    }

    .ppm_check_finite(init.theta, iter)

    ## ---- optional step-halving safeguard on the original objective Q ------
    if (line.search) {
      ls <- .ppm_step_halving(old.theta, init.theta, Q.old, ctx, max.halving, tol)
      init.theta <- ls$theta; Q.path[iter] <- ls$Q; m.path[iter] <- ls$m
      if (!ls$accepted) {
        status <- if (ls$halted) "line.search.halted" else "converged"; break
      }
      Q.old <- ls$Q
    }

    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))
    if (delta < tol) { status <- "converged"; break }
  }

  .ppm_finish(init.theta, h_eff, x,
              .ppm_diag(line.search, iter, status, Q.path, m.path))
}


# ============================================================================
# 5. ppqr  (P2QR; MM-GCD, line search optional)
#
# The floor eps of the check-loss majorizer, a_ik = max(|r_ik|, eps), is by
# default set relative to the response scale, eps = 0.05 * sd(y). A very small
# absolute floor (e.g. 1e-6) makes the MM weights 1/(4 a_ik) explode near the
# quantile-regression solution, where p + 1 residuals vanish exactly, and the
# fixed-step group updates then fail to converge or diverge. With the relative
# floor the majorizer coincides with that of Hunter and Lange (2000) for every
# residual with |r_ik| >= eps and the gap at the current iterate is at most
# eps / 4 per observation.
# ============================================================================
ppqr <- function(x, y, H = 10, C = 1, lambda, gamma = 3.7, penalty = "grSCAD",
                 max.iter = 100, ridge = 1e-10, eps = NULL,
                 line.search = FALSE, max.halving = 10L) {
  n <- nrow(x)
  p <- ncol(x)

  if (!penalty %in% c("grSCAD", "grLasso", "grMCP")) stop("Invalid penalty")
  if (is.null(eps)) eps <- 0.05 * stats::sd(y)
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0)
    stop("eps must be a single positive number (for a constant response sd(y) = 0; supply 'eps' explicitly)")

  bar.x   <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))

  qprob <- seq(1/H, 1 - 1/H, by = 1/H)
  h_eff <- length(qprob)
  tmp.y <- rep(y, times = h_eff)                  # y stacked machine by machine

  Sigma             <- cov(x)
  Sigma.tilde       <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma))
  Sigma.tilde.ridge <- Sigma.tilde + diag(ridge, p + 1)

  init.theta <- rep(0, h_eff * (p + 1))
  tau_vec    <- rep(qprob, each = n)              # tau_k for each (i, k)
  tol        <- 1e-5

  ## ---- line-search bookkeeping -------------------------------------------
  Q.path <- numeric(max.iter); m.path <- integer(max.iter); status <- "max.iter"
  if (line.search) {
    ctx   <- .ppm_make_ctx("qr", x.tilde, Sigma.tilde, Y = matrix(y, n, h_eff), tau = qprob,
                           C = C, lambda = lambda, gamma = gamma, penalty = penalty)
    Q.old <- .ppm_objective(init.theta, ctx)
  }

  for (iter in 1:max.iter) {
    old.theta <- init.theta

    # f_{ik}^(t) = theta_k^(t)' x_tilde_i  (unweighted fitted values)
    Wtheta <- numeric(n * h_eff)
    for (k in 1:h_eff) {
      theta_k <- init.theta[((k - 1) * (p + 1) + 1):(k * (p + 1))]
      idx_k   <- ((k - 1) * n + 1):(k * n)
      Wtheta[idx_k] <- as.vector(x.tilde %*% theta_k)
    }

    # ---- Safeguarded MM quantities (recomputed each iteration) -------------
    r_t          <- as.vector(tmp.y - Wtheta)     # r = y - f
    a_t          <- pmax(abs(r_t), eps)           # a = max(|r|, eps)
    omega_t      <- 1 / (4 * a_t)                 # Q = 1 / (4 a)
    sqrt_omega_t <- sqrt(omega_t)

    u.tilde     <- tmp.y + (2 * tau_vec - 1) * a_t   # u = y + (2 tau - 1) a
    u.tilde_new <- sqrt_omega_t * u.tilde

    A_blocks <- vector("list", h_eff)
    for (k in 1:h_eff) {
      idx_k       <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * sqrt_omega_t[idx_k]
      A_blocks[[k]] <- Sigma.tilde.ridge + crossprod(Wk_weighted)
    }

    A_sqrt_blocks <- sqrt_sym_blockdiag(A_blocks)

    rhs <- numeric(h_eff * (p + 1))
    for (k in 1:h_eff) {
      idx_k       <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * sqrt_omega_t[idx_k]
      rhs[((k - 1) * (p + 1) + 1):(k * (p + 1))] <- t(Wk_weighted) %*% u.tilde_new[idx_k]
    }
    big_Y <- solve_blockdiag(A_blocks, rhs)
    res <- big_Y - multiply_blockdiag(A_sqrt_blocks, init.theta)

    group <- rep(1:(p + 1), times = h_eff)
    J     <- p + 1

    for (j in 1:J) {
      ind <- which(group == j)
      G_cols <- get_var_columns(A_sqrt_blocks, j)

      z_j <- (C / n) * (t(G_cols) %*% res) + init.theta[ind]
      z_norm <- norm(z_j, "2")

      if (z_norm < 1e-12) {
        theta_j_new <- rep(0, length(ind))
      } else {
        if (penalty == "grSCAD") {
          theta_j_new <- scad_firm_thresholding(z_norm, lambda1 = lambda, lambda2 = 2, gamma = gamma) * z_j / z_norm
        } else if (penalty == "grLasso") {
          theta_j_new <- soft_thresholding(z_norm, lambda = lambda) / (1 + 2) * z_j / z_norm
        } else if (penalty == "grMCP") {
          theta_j_new <- firm_thresholding(z_norm, lambda1 = lambda, lambda2 = 2, gamma = gamma) * z_j / z_norm
        }
        theta_j_new <- as.vector(theta_j_new)
      }

      res <- res - G_cols %*% (theta_j_new - init.theta[ind])
      init.theta[ind] <- theta_j_new
      res <- as.vector(res)
    }

    .ppm_check_finite(init.theta, iter)

    ## ---- optional step-halving safeguard on the original objective Q ------
    if (line.search) {
      ls <- .ppm_step_halving(old.theta, init.theta, Q.old, ctx, max.halving, tol)
      init.theta <- ls$theta; Q.path[iter] <- ls$Q; m.path[iter] <- ls$m
      if (!ls$accepted) {
        status <- if (ls$halted) "line.search.halted" else "converged"; break
      }
      Q.old <- ls$Q
    }

    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))
    if (delta < tol) { status <- "converged"; break }
  }

  .ppm_finish(init.theta, h_eff, x,
              .ppm_diag(line.search, iter, status, Q.path, m.path))
}

# ============================================================================
# 6. pplssvm  (P2LSM; exact GCD on a single quadratic problem, no line search)
# ============================================================================
pplssvm <- function(x, y, H = 10, C = 1, lambda, gamma = 3.7, penalty = "grSCAD", max.iter = 100, ridge = 1e-10) {
  n <- nrow(x)
  p <- ncol(x)

  if (!penalty %in% c("grSCAD", "grLasso", "grMCP")) stop("Invalid penalty")

  bar.x   <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))

  qprob <- seq(1/H, 1 - 1/H, by = 1/H)
  h_eff <- length(qprob)
  qy    <- quantile(y, qprob)
  U     <- sapply(qy, function(z) 2 * (z > y) - 1)
  y.tilde <- c(U)

  Sigma.hat      <- cov(x)
  Sigma.hat.star <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma.hat))

  # For LS-SVM, W is constant, so A_blocks is calculated exactly once
  A_blocks <- vector("list", h_eff)
  for (k in 1:h_eff) {
    A_blocks[[k]] <- (n/C) * Sigma.hat.star + crossprod(x.tilde) + diag(ridge, p + 1)
  }

  A_sqrt_blocks <- sqrt_sym_blockdiag(A_blocks)

  rhs <- numeric(h_eff * (p + 1))
  for (k in 1:h_eff) {
    idx_k <- ((k - 1) * n + 1):(k * n)
    rhs[((k - 1) * (p + 1) + 1):(k * (p + 1))] <- C * t(x.tilde) %*% y.tilde[idx_k]
  }

  xi <- solve_blockdiag(A_blocks, rhs)
  xi.tilde <- as.vector(xi)

  group <- rep(1:(p + 1), times = h_eff)
  J     <- p + 1

  init.theta <- rep(0, h_eff * (p + 1))
  tol        <- 1e-5
  res <- xi.tilde - multiply_blockdiag(A_sqrt_blocks, init.theta)
  status <- "max.iter"

  for (iter in 1:max.iter) {
    old.theta <- init.theta
    for (j in 1:J) {
      ind <- which(group == j)
      G_cols <- get_var_columns(A_sqrt_blocks, j)

      z_j <- (C / n) * (t(G_cols) %*% res) + init.theta[ind]
      z_norm <- norm(z_j, "2")

      if (z_norm < 1e-12) {
        theta_j_new <- rep(0, length(ind))
      } else {
        if (penalty == "grSCAD") {
          theta_j_new <- scad_firm_thresholding(z_norm, lambda1 = lambda, lambda2 = 2, gamma = gamma) * z_j / z_norm
        } else if (penalty == "grLasso") {
          theta_j_new <- soft_thresholding(z_norm, lambda = lambda) / (1 + 2) * z_j / z_norm
        } else if (penalty == "grMCP") {
          theta_j_new <- firm_thresholding(z_norm, lambda1 = lambda, lambda2 = 2, gamma = gamma) * z_j / z_norm
        }
        theta_j_new <- as.vector(theta_j_new)
      }

      res <- res - G_cols %*% (theta_j_new - init.theta[ind])
      init.theta[ind] <- theta_j_new
      res <- as.vector(res)
    }

    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))
    if (delta < tol) { status <- "converged"; break }
  }

  .ppm_finish(init.theta, h_eff, x,
              .ppm_diag(FALSE, iter, status, NULL, NULL))
}


# ============================================================================
# 7. ppwlssvm  (P2WLSM; exact GCD via grpreg, no line search)
# ============================================================================
ppwlssvm <- function(x, y, H = 10, C = 1, lambda, gamma = 3.7, penalty = "grSCAD", max.iter = 100, ridge = 1e-10) {
  n <- nrow(x); p <- ncol(x)
  bar.x <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))

  qprob <- seq(1/H, 1-1/H, by = 1/H)
  h_eff <- length(qprob)
  tmp.y <- rep(y, times = h_eff)
  pi.grid <- rep(qprob, each = n)

  weight <- (1 - pi.grid) * as.numeric(tmp.y == 1) + pi.grid * as.numeric(tmp.y == -1)
  sqrt_weight <- sqrt(weight)
  u.tilde <- sqrt_weight * tmp.y

  Sigma.hat <- cov(x)
  Sigma.hat.star <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma.hat))

  A_blocks <- vector("list", h_eff)
  for (k in 1:h_eff) {
    idx_k <- ((k - 1) * n + 1):(k * n)
    Wk_weighted <- x.tilde * sqrt_weight[idx_k]
    A_blocks[[k]] <- (n / C) * Sigma.hat.star + crossprod(Wk_weighted) + diag(ridge, p + 1)
  }

  A_sqrt_blocks <- sqrt_sym_blockdiag(A_blocks)
  G.tilde.big <- C * as.matrix(Matrix::bdiag(A_sqrt_blocks))

  rhs <- numeric(h_eff * (p + 1))
  for (k in 1:h_eff) {
    idx_k <- ((k - 1) * n + 1):(k * n)
    Wk_weighted <- x.tilde * sqrt_weight[idx_k]
    rhs[((k - 1) * (p + 1) + 1):(k * (p + 1))] <- t(Wk_weighted) %*% u.tilde[idx_k]
  }

  xi <- solve_blockdiag(A_sqrt_blocks, rhs)
  xi.tilde <- as.vector(xi) / C

  group <- rep(1:(p + 1), times = h_eff)
  colnames(G.tilde.big) <- group

  obj_grpreg <- grpreg::grpreg(X = G.tilde.big, y = xi.tilde, group = group, penalty = penalty, family = "gaussian",
                       lambda = lambda, alpha = 1, eps = 1e-5, max.iter = max.iter, dfmax = p,
                       gmax = length(unique(group)), gamma = gamma)

  theta <- obj_grpreg$beta[-1]
  iter  <- if (!is.null(obj_grpreg$iter)) as.integer(obj_grpreg$iter[1L]) else NA_integer_

  .ppm_finish(theta, h_eff, x,
              .ppm_diag(FALSE, iter, "converged", NULL, NULL))
}

# ============================================================================
# 8. ppwlr  (P2WLR; iterative GCD)
# ============================================================================
ppwlr <- function(x, y, H, C, lambda, gamma = 3.7,
                  penalty = "grSCAD", max.iter = 100, ridge = 1e-10,
                  line.search = TRUE, max.halving = 10L) {

  n <- nrow(x); p <- ncol(x)
  if (!penalty %in% c("grSCAD", "grLasso", "grMCP")) stop("Invalid penalty")

  bar.x   <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))

  qprob   <- seq(1/H, 1 - 1/H, by = 1/H)
  h_eff   <- length(qprob)
  pi.grid <- rep(qprob, each = n)

  u <- rep(y, h_eff)                              # y in {-1, 1}, stacked machine by machine

  # class weights: w_ik = 1 - pi_k (y = 1),  pi_k (y = -1)
  Omega_vec <- (1 - pi.grid) * as.numeric(u == 1) + pi.grid * as.numeric(u == -1)

  Sigma       <- cov(x)
  Sigma.tilde <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma))

  init.theta <- rep(0, h_eff * (p + 1))
  tol <- 1e-5

  # Penalty-dependent S scaling factor (kept as in the previous version)
  S_scale <- if (penalty == "grSCAD") n / C else 1

  ## ---- line-search bookkeeping -------------------------------------------
  Q.path <- numeric(max.iter); m.path <- integer(max.iter); status <- "max.iter"
  if (line.search) {
    ctx   <- .ppm_make_ctx("logit", x.tilde, Sigma.tilde, Y = matrix(u, n, h_eff),
                           W = matrix(Omega_vec, n, h_eff),
                           C = C, lambda = lambda, gamma = gamma, penalty = penalty)
    Q.old <- .ppm_objective(init.theta, ctx)
  }

  for (iter in 1:max.iter) {
    old.theta <- init.theta

    # f_{ik}^(t) = theta_k^(t)' x_tilde_i  (unweighted fitted values)
    Wtheta <- numeric(n * h_eff)
    for (k in 1:h_eff) {
      theta_k <- init.theta[((k - 1) * (p + 1) + 1):(k * (p + 1))]
      idx_k   <- ((k - 1) * n + 1):(k * n)
      Wtheta[idx_k] <- as.vector(x.tilde %*% theta_k)
    }

    # ---- IRLS quantities (class weight enters the curvature only) ----------
    q_vec   <- 1 / (1 + exp(u * Wtheta))          # q = 1 / (1 + exp(y f))
    h_vec   <- Omega_vec * q_vec * (1 - q_vec)    # H = w q (1 - q)
    sqrt_h  <- sqrt(h_vec)

    u.tilde     <- Wtheta + u / (1 - q_vec)       # u = f + y / (1 - q), for every penalty
    u.tilde_new <- sqrt_h * u.tilde

    A_blocks <- vector("list", h_eff)
    for (k in 1:h_eff) {
      idx_k       <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * sqrt_h[idx_k]
      A_blocks[[k]] <- S_scale * Sigma.tilde + crossprod(Wk_weighted)
    }

    A_sqrt_blocks <- sqrt_sym_blockdiag(A_blocks, ridge = ridge)

    rhs <- numeric(h_eff * (p + 1))
    for (k in 1:h_eff) {
      idx_k       <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * sqrt_h[idx_k]
      rhs[((k - 1) * (p + 1) + 1):(k * (p + 1))] <- t(Wk_weighted) %*% u.tilde_new[idx_k]
    }

    big_Y <- solve_blockdiag(A_blocks, rhs, ridge = ridge)
    res   <- big_Y - multiply_blockdiag(A_sqrt_blocks, init.theta)

    group <- rep(1:(p + 1), times = h_eff)
    J     <- p + 1

    for (j in 1:J) {
      ind    <- which(group == j)
      G_cols <- get_var_columns(A_sqrt_blocks, j)

      z_j    <- (C / n) * (t(G_cols) %*% res) + init.theta[ind]
      z_norm <- norm(z_j, "2")

      if (z_norm < 1e-12) {
        theta_j_new <- rep(0, length(ind))
      } else {
        if (penalty == "grSCAD") {
          theta_j_new <- scad_firm_thresholding(z_norm, lambda1 = lambda, lambda2 = 2, gamma = gamma) * z_j / z_norm
        } else if (penalty == "grLasso") {
          theta_j_new <- soft_thresholding(z_norm, lambda = lambda) / (1 + 2) * z_j / z_norm
        } else if (penalty == "grMCP") {
          # kept as in the previous version (lambda2 = 1)
          theta_j_new <- firm_thresholding(z_norm, lambda1 = lambda, lambda2 = 1, gamma = gamma) * z_j / z_norm
        }
        theta_j_new <- as.vector(theta_j_new)
      }

      res <- res - G_cols %*% (theta_j_new - init.theta[ind])
      init.theta[ind] <- theta_j_new
      res <- as.vector(res)
    }

    .ppm_check_finite(init.theta, iter)

    ## ---- step-halving safeguard on the original objective Q ---------------
    if (line.search) {
      ls <- .ppm_step_halving(old.theta, init.theta, Q.old, ctx, max.halving, tol)
      init.theta <- ls$theta; Q.path[iter] <- ls$Q; m.path[iter] <- ls$m
      if (!ls$accepted) {
        status <- if (ls$halted) "line.search.halted" else "converged"; break
      }
      Q.old <- ls$Q
    }

    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))
    if (delta < tol) { status <- "converged"; break }
  }

  .ppm_finish(init.theta, h_eff, x,
              .ppm_diag(line.search, iter, status, Q.path, m.path))
}

# ============================================================================
# 9. ppwl2svm  (P2WL2M; iterative GCD)
# ============================================================================
ppwl2svm <- function(x, y, H = 10, C = 1, lambda, gamma = 3.7, penalty = "grSCAD", max.iter = 100, ridge = 1e-4,
                     line.search = TRUE, max.halving = 10L) {
  n <- nrow(x); p <- ncol(x)
  bar.x   <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))

  qprob <- seq(1/H, 1 - 1/H, by = 1/H)
  h_eff <- length(qprob)
  tmp.y <- rep(y, times = h_eff)                  # y in {-1, 1}, stacked machine by machine
  pi.grid <- rep(qprob, each = n)

  # class weights: w_ik = 1 - pi_k (y = 1),  pi_k (y = -1)
  weight <- (1 - pi.grid) * as.numeric(tmp.y == 1) + pi.grid * as.numeric(tmp.y == -1)
  sqrt_weight <- sqrt(weight)
  weighted_y <- sqrt_weight * tmp.y

  Sigma <- cov(x)
  Sigma.hat.star <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma))
  init.theta <- rep(0, h_eff * (p + 1))
  tol        <- 1e-5

  ## ---- line-search bookkeeping -------------------------------------------
  Q.path <- numeric(max.iter); m.path <- integer(max.iter); status <- "max.iter"
  if (line.search) {
    ctx   <- .ppm_make_ctx("l2", x.tilde, Sigma.hat.star, Y = matrix(tmp.y, n, h_eff),
                           W = matrix(weight, n, h_eff),
                           C = C, lambda = lambda, gamma = gamma, penalty = penalty)
    Q.old <- .ppm_objective(init.theta, ctx)
  }

  for (iter in 1:max.iter) {
    old.theta <- init.theta

    # f_{ik}^(t) = theta_k^(t)' x_tilde_i  (UNWEIGHTED fitted values;
    # the margin must not be scaled by sqrt(w_ik))
    f_t <- numeric(n * h_eff)
    for (k in 1:h_eff) {
      theta_k <- init.theta[((k - 1) * (p + 1) + 1):(k * (p + 1))]
      idx_k   <- ((k - 1) * n + 1):(k * n)
      f_t[idx_k] <- as.vector(x.tilde %*% theta_k)
    }

    # ---- active set: 1 - y * f > 0 (margin violators) ----------------------
    u_mask <- as.integer(1 - tmp.y * f_t > 0)
    y.tilde <- weighted_y * sqrt(u_mask)          # sqrt(w_ik) * 1(A_k) * y_ik

    A_blocks <- vector("list", h_eff)
    for (k in 1:h_eff) {
      idx_k <- ((k - 1) * n + 1):(k * n)
      wTw_k <- crossprod(x.tilde * (sqrt_weight[idx_k] * sqrt(u_mask[idx_k])))
      A_blocks[[k]] <- 0.5 * Sigma.hat.star + wTw_k
    }

    A_sqrt_blocks <- sqrt_sym_blockdiag(A_blocks, ridge = ridge)
    rhs <- numeric(h_eff * (p + 1))
    for (k in 1:h_eff) {
      idx_k <- ((k - 1) * n + 1):(k * n)
      rhs[((k - 1) * (p + 1) + 1):(k * (p + 1))] <- t(x.tilde * (sqrt_weight[idx_k] * sqrt(u_mask[idx_k]))) %*% y.tilde[idx_k]
    }

    xi <- solve_blockdiag(A_sqrt_blocks, rhs, ridge = ridge)
    res <- as.vector(xi) - multiply_blockdiag(A_sqrt_blocks, init.theta)

    group <- rep(1:(p + 1), times = h_eff)
    J     <- p + 1

    for (j in 1:J) {
      ind <- which(group == j)
      G_cols <- get_var_columns(A_sqrt_blocks, j)
      z_j <- (C / n) * (t(G_cols) %*% res) + init.theta[ind]
      z_norm <- norm(z_j, "2")

      if (z_norm < 1e-12) theta_j_new <- rep(0, length(ind))
      else {
        if (penalty == "grSCAD") theta_j_new <- scad_firm_thresholding(z_norm, lambda1 = lambda, lambda2 = 2, gamma = gamma) * z_j / z_norm
        else if (penalty == "grLasso") theta_j_new <- soft_thresholding(z_norm, lambda = lambda) / (1 + 2) * z_j / z_norm
        else if (penalty == "grMCP") theta_j_new <- firm_thresholding(z_norm, lambda1 = lambda, lambda2 = 2, gamma = gamma) * z_j / z_norm
        theta_j_new <- as.vector(theta_j_new)
      }
      res <- res - G_cols %*% (theta_j_new - init.theta[ind])
      init.theta[ind] <- theta_j_new
      res <- as.vector(res)
    }

    .ppm_check_finite(init.theta, iter)

    ## ---- step-halving safeguard on the original objective Q ---------------
    if (line.search) {
      ls <- .ppm_step_halving(old.theta, init.theta, Q.old, ctx, max.halving, tol)
      init.theta <- ls$theta; Q.path[iter] <- ls$Q; m.path[iter] <- ls$m
      if (!ls$accepted) {
        status <- if (ls$halted) "line.search.halted" else "converged"; break
      }
      Q.old <- ls$Q
    }

    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))
    if (delta < tol) { status <- "converged"; break }
  }

  .ppm_finish(init.theta, h_eff, x,
              .ppm_diag(line.search, iter, status, Q.path, m.path))
}

# ============================================================================
# 10. ppwsvm  (P2WSVM; MM-GCD, line search optional)
# ============================================================================
ppwsvm <- function(x, y, H = 10, C = 1, lambda, gamma = 3.7, penalty = "grSCAD",
                   max.iter = 100, ridge = 1e-10, eps = 1e-6,
                   line.search = FALSE, max.halving = 10L) {
  n <- nrow(x); p <- ncol(x)
  if (!penalty %in% c("grSCAD", "grLasso", "grMCP")) stop("Invalid penalty")
  if (eps <= 0) stop("eps must be positive")

  bar.x   <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))

  qprob <- seq(1/H, 1 - 1/H, by = 1/H)
  h_eff <- length(qprob)
  tmp.y <- rep(y, times = h_eff)                  # y in {-1, 1}, stacked machine by machine
  pi.grid <- rep(qprob, each = n)

  Sigma <- cov(x)
  Sigma.tilde <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma))
  init.theta <- rep(0, h_eff * (p + 1))
  tol <- 1e-5

  # class weights: w_ik = 1 - pi_k (y = 1),  pi_k (y = -1)
  weight_t_base <- (1 - pi.grid) * as.numeric(tmp.y == 1) + pi.grid * as.numeric(tmp.y == -1)

  ## ---- line-search bookkeeping -------------------------------------------
  Q.path <- numeric(max.iter); m.path <- integer(max.iter); status <- "max.iter"
  if (line.search) {
    ctx   <- .ppm_make_ctx("hinge", x.tilde, Sigma.tilde, Y = matrix(tmp.y, n, h_eff),
                           W = matrix(weight_t_base, n, h_eff),
                           C = C, lambda = lambda, gamma = gamma, penalty = penalty)
    Q.old <- .ppm_objective(init.theta, ctx)
  }

  for (iter in 1:max.iter) {
    old.theta <- init.theta

    # f_{ik}^(t) = theta_k^(t)' x_tilde_i  (unweighted, no division by n)
    Wtheta <- numeric(n * h_eff)
    for (k in 1:h_eff) {
      theta_k <- init.theta[((k - 1) * (p + 1) + 1):(k * (p + 1))]
      idx_k   <- ((k - 1) * n + 1):(k * n)
      Wtheta[idx_k] <- as.vector(x.tilde %*% theta_k)
    }

    # ---- Safeguarded MM weights (recomputed each iteration) ----------------
    a_t          <- pmax(abs(1 - tmp.y * Wtheta), eps)   # a = max(|1 - y f|, eps)
    omega_t      <- weight_t_base / (4 * a_t)            # Omega = w / (4 a)
    sqrt_omega_t <- sqrt(omega_t)

    u.tilde     <- (1 + a_t) * tmp.y                     # working response (1 + a) y
    u.tilde_new <- sqrt_omega_t * u.tilde

    A_blocks <- vector("list", h_eff)
    for (k in 1:h_eff) {
      idx_k <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * sqrt_omega_t[idx_k]
      A_blocks[[k]] <- Sigma.tilde + crossprod(Wk_weighted)
    }

    A_sqrt_blocks <- sqrt_sym_blockdiag(A_blocks, ridge = ridge)
    rhs <- numeric(h_eff * (p + 1))
    for (k in 1:h_eff) {
      idx_k <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * sqrt_omega_t[idx_k]
      rhs[((k - 1) * (p + 1) + 1):(k * (p + 1))] <- t(Wk_weighted) %*% u.tilde_new[idx_k]
    }

    big_Y <- solve_blockdiag(A_blocks, rhs, ridge = ridge)
    res <- big_Y - multiply_blockdiag(A_sqrt_blocks, init.theta)

    group <- rep(1:(p + 1), times = h_eff)
    J     <- p + 1

    for (j in 1:J) {
      ind <- which(group == j)
      G_cols <- get_var_columns(A_sqrt_blocks, j)
      z_j <- (C / n) * (t(G_cols) %*% res) + init.theta[ind]
      z_norm <- norm(z_j, "2")

      if (z_norm < 1e-12) theta_j_new <- rep(0, length(ind))
      else {
        if (penalty == "grSCAD") theta_j_new <- scad_firm_thresholding(z_norm, lambda1 = lambda, lambda2 = 2, gamma = gamma) * z_j / z_norm
        else if (penalty == "grLasso") theta_j_new <- soft_thresholding(z_norm, lambda = lambda) / (1 + 2) * z_j / z_norm
        else if (penalty == "grMCP") theta_j_new <- firm_thresholding(z_norm, lambda1 = lambda, lambda2 = 2, gamma = gamma) * z_j / z_norm
        theta_j_new <- as.vector(theta_j_new)
      }
      res <- res - G_cols %*% (theta_j_new - init.theta[ind])
      init.theta[ind] <- theta_j_new
      res <- as.vector(res)
    }

    .ppm_check_finite(init.theta, iter)

    ## ---- optional step-halving safeguard on the original objective Q ------
    if (line.search) {
      ls <- .ppm_step_halving(old.theta, init.theta, Q.old, ctx, max.halving, tol)
      init.theta <- ls$theta; Q.path[iter] <- ls$Q; m.path[iter] <- ls$m
      if (!ls$accepted) {
        status <- if (ls$halted) "line.search.halted" else "converged"; break
      }
      Q.old <- ls$Q
    }

    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))
    if (delta < tol) { status <- "converged"; break }
  }

  .ppm_finish(init.theta, h_eff, x,
              .ppm_diag(line.search, iter, status, Q.path, m.path))
}
