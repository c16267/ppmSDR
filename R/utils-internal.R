# ============================================================================
# Internal helper operators for the group coordinate descent (GCD) engine.
# None of these are exported. They are shared by all P2M solvers in
# estimators.R. Kept as a single source to avoid duplicate definitions.
# ============================================================================

# Soft-thresholding operator (group LASSO update)
soft_thresholding <- function(z, lambda) {
  if (z > lambda) return(z - lambda)
  if (z < -lambda) return(z + lambda)
  return(0)
}

# Firm-thresholding operator (group MCP update)
firm_thresholding <- function(z, lambda1, lambda2, gamma) {
  s <- 0
  if (z > 0) s <- 1
  else if (z < 0) s <- -1
  if (abs(z) <= lambda1) return(0)
  else if (abs(z) <= gamma * lambda1 * (1 + lambda2)) {
    return(s * (abs(z) - lambda1) / (1 + lambda2 - 1 / gamma))
  } else {
    return(z / (1 + lambda2))
  }
}

# SCAD-modified firm-thresholding operator (group SCAD update)
scad_firm_thresholding <- function(z, lambda1, lambda2, gamma) {
  s <- 0
  if (z > 0) {
    s <- 1
  } else if (z < 0) {
    s <- -1
  } else {
    s <- 0
  }

  if (abs(z) <= lambda1) {
    return(0)
  } else if (abs(z) <= (lambda1 * (1 + lambda2) + lambda1)) {
    return(s * (abs(z) - lambda1) / (1 + lambda2))
  } else if (abs(z) <= gamma * lambda1 * (1 + lambda2)) {
    return(s * (abs(z) - gamma * lambda1 / (gamma - 1)) /
             (1 - 1 / (gamma - 1) + lambda2))
  } else {
    return(z / (1 + lambda2))
  }
}

# Per-slice symmetric square root via block-diagonal structure
sqrt_sym_blockdiag <- function(block_list, ridge = 0) {
  h <- length(block_list)
  out <- vector("list", h)
  for (k in 1:h) {
    A_k <- block_list[[k]]
    if (ridge > 0) A_k <- A_k + diag(ridge, nrow(A_k))
    eig <- eigen(A_k, symmetric = TRUE)
    out[[k]] <- eig$vectors %*%
      diag(sqrt(pmax(eig$values, 0))) %*%
      t(eig$vectors)
  }
  out
}

# Solve (G x = rhs) where G is block-diagonal (list of blocks)
solve_blockdiag <- function(block_list, rhs, ridge = 0) {
  h <- length(block_list)
  m <- nrow(block_list[[1]])
  stopifnot(length(rhs) == h * m)
  out <- numeric(h * m)
  for (k in 1:h) {
    idx <- ((k - 1) * m + 1):(k * m)
    A_k <- block_list[[k]]
    if (ridge > 0) A_k <- A_k + diag(ridge, nrow(A_k))
    out[idx] <- tryCatch(
      solve(A_k, rhs[idx]),
      error = function(e) stop(
        "The quadratic subproblem of machine ", k, " is singular (",
        conditionMessage(e), "). This typically happens when the GCD sweep ",
        "has diverged, e.g. for a large cost parameter 'C'; reduce 'C' or set ",
        "line.search = TRUE.", call. = FALSE))
  }
  out
}

# Divergence guard called after every GCD sweep of the iterative solvers.
.ppm_check_finite <- function(theta, iter) {
  if (!all(is.finite(theta)))
    stop("The GCD sweep diverged at iteration ", iter,
         " (non-finite coefficients); reduce 'C' or set line.search = TRUE.",
         call. = FALSE)
  invisible(TRUE)
}

# Compute (G x) where G is block-diagonal
multiply_blockdiag <- function(block_list, x) {
  h <- length(block_list)
  m <- nrow(block_list[[1]])
  stopifnot(length(x) == h * m)
  out <- numeric(h * m)
  for (k in 1:h) {
    idx <- ((k - 1) * m + 1):(k * m)
    out[idx] <- as.vector(block_list[[k]] %*% x[idx])
  }
  out
}

# Extract a column-subset across all blocks for the GCD inner loop
get_var_columns <- function(block_list, var_idx) {
  h <- length(block_list)
  m <- nrow(block_list[[1]])
  out <- matrix(0, h * m, h)
  for (k in 1:h) {
    row_idx <- ((k - 1) * m + 1):(k * m)
    out[row_idx, k] <- block_list[[k]][, var_idx]
  }
  out
}


# ============================================================================
# Objective function Q(theta) of the penalized principal machine (P2M) and
# the step-halving line search (Section 4 and Remark 1 of the article).
#
#   Q(theta) = sum_k theta_k' Sigma~ theta_k
#            + (c/n) sum_k sum_i w_ik L_k(y~_ik, theta_k' x~_i)
#            + sum_{j=1}^{p} p_lambda(||theta_(j)||),
#
# where theta = vec(Theta) is stacked machine by machine (column-wise), the
# intercept row of Theta is left unpenalized, and theta_(j) is the j-th
# predictor row of Theta across the h machines.
#
# The safeguarded update is
#   theta^{(t+1)} = theta^{(t)} + 2^{-m_t} (theta_hat^{(t)} - theta^{(t)}),
# where m_t >= 0 is the smallest integer with Q(theta^{(t+1)}) <= Q(theta^{(t)}),
# and the algorithm stops at theta^{(t)} if no such m_t <= max.halving exists.
# ============================================================================

# Group penalty p_lambda(t) for a vector of non-negative row norms t.
.ppm_penalty <- function(t, lambda, gamma, penalty) {
  t <- abs(t)
  switch(penalty,
    grLasso = lambda * t,
    grSCAD  = ifelse(t <= lambda, lambda * t,
                ifelse(t <= gamma * lambda,
                       (2 * gamma * lambda * t - t^2 - lambda^2) /
                         (2 * (gamma - 1)),
                       (gamma + 1) * lambda^2 / 2)),
    grMCP   = ifelse(t <= gamma * lambda,
                     lambda * t - t^2 / (2 * gamma),
                     gamma * lambda^2 / 2),
    stop("Invalid penalty")
  )
}

# Element-wise loss L(y, f) for the loss family used by a solver.
# Y and F are n x h matrices; tau is a length-h vector (asls/qr only).
.ppm_loss_matrix <- function(family, Y, F, tau = NULL) {
  switch(family,
    ls    = (1 - Y * F)^2,
    logit = {                              # log(1 + exp(-m)), overflow-safe
      m <- Y * F
      pmax(-m, 0) + log1p(exp(-abs(m)))
    },
    l2    = pmax(0, 1 - Y * F)^2,
    hinge = pmax(0, 1 - Y * F),
    asls  = {                              # tau (r >= 0), (1 - tau) (r < 0)
      r  <- Y - F
      tt <- matrix(tau, nrow(r), ncol(r), byrow = TRUE)
      r^2 * ifelse(r >= 0, tt, 1 - tt)
    },
    qr    = {                              # check loss r (tau - 1{r < 0})
      r  <- Y - F
      tt <- matrix(tau, nrow(r), ncol(r), byrow = TRUE)
      r * (tt - (r < 0))
    },
    stop("Unknown loss family '", family, "'")
  )
}

# Context needed to evaluate Q(theta); built once per solver call.
#   x.tilde    : n x (p+1) design with intercept column
#   Sigma.star : (p+1) x (p+1), diag(0, Sigma_hat)
#   Y          : n x h matrix of (pseudo-)responses, one column per machine
#   W          : n x h matrix of class weights w_ik, or NULL for w_ik = 1
#   tau        : length-h asymmetry/quantile levels, or NULL
.ppm_make_ctx <- function(family, x.tilde, Sigma.star, Y, W = NULL, tau = NULL,
                          C, lambda, gamma, penalty) {
  list(family = family, x.tilde = x.tilde, Sigma.star = Sigma.star,
       Y = Y, W = W, tau = tau, n = nrow(x.tilde), h = ncol(Y),
       C = C, lambda = lambda, gamma = gamma, penalty = penalty)
}

# Q(theta) for the stacked coefficient vector theta = vec(Theta).
.ppm_objective <- function(theta, ctx) {
  Theta <- matrix(theta, ncol = ctx$h)                 # (p+1) x h
  F     <- ctx$x.tilde %*% Theta                        # n x h fitted values
  cov.term  <- sum(Theta * (ctx$Sigma.star %*% Theta))  # sum_k theta_k' S~ theta_k
  L <- .ppm_loss_matrix(ctx$family, ctx$Y, F, ctx$tau)
  if (!is.null(ctx$W)) L <- ctx$W * L
  loss.term <- (ctx$C / ctx$n) * sum(L)
  rn <- sqrt(rowSums(Theta[-1L, , drop = FALSE]^2))     # predictor rows only
  pen.term  <- sum(.ppm_penalty(rn, ctx$lambda, ctx$gamma, ctx$penalty))
  cov.term + loss.term + pen.term
}

# Step-halving line search on the original objective Q.
# Returns the accepted iterate, its objective value, the number of halvings
# m_t and a flag telling whether an acceptable step was found. When no step
# is acceptable, `halted` is TRUE only if the rejected plain step would not
# have satisfied the stopping rule, i.e. its relative sup-norm is at least
# `tol`; the rejection of a smaller step means the iteration has converged
# and is reported as "converged".
.ppm_step_halving <- function(theta.old, theta.hat, Q.old, ctx,
                              max.halving = 10L, tol = 1e-5, ls.tol = 1e-10) {
  d   <- theta.hat - theta.old
  bar <- Q.old + ls.tol * (1 + abs(Q.old))   # numerical slack for "<="
  for (m in 0:max.halving) {
    theta.try <- theta.old + 2^(-m) * d
    Q.try     <- .ppm_objective(theta.try, ctx)
    if (is.finite(Q.try) && Q.try <= bar)
      return(list(theta = theta.try, Q = Q.try, m = m, accepted = TRUE,
                  halted = FALSE))
  }
  delta.hat <- max(abs(d)) / (1 + max(abs(theta.old)))
  list(theta = theta.old, Q = Q.old, m = NA_integer_, accepted = FALSE,
       halted = delta.hat >= tol)
}

# Assemble the convergence diagnostics returned by every iterative solver.
.ppm_diag <- function(line.search, iter, status, Q.path, m.path) {
  list(line.search = line.search, iter = iter, status = status,
       objective = if (line.search) Q.path[seq_len(iter)] else NULL,
       n.halving = if (line.search) m.path[seq_len(iter)] else NULL)
}
