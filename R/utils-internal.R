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
    out[idx] <- solve(A_k, rhs[idx])
  }
  out
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
