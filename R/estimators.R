# Internal penalized-principal-machine (P2M) estimators.
# These are NOT exported; the user-facing entry point is ppm().
# Each solver returns: M (working matrix), evalues, evectors, x.
# Helper operators (thresholding, block-diagonal algebra) live in utils-internal.R.

# ============================================================================
# 1. ppasls
# ============================================================================
ppasls <- function(x, y, H, C, lambda, gamma = 3.7, penalty = "grSCAD", max.iter = 100) {
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
      
      A_blocks[[k]] <- 0.5 * Sigma.hat.star + wTw_k
      B_blocks[[k]] <- 0.5 * Sigma.hat.star + wTw_c_k
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
    
    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))
    
    if (delta < tol) break
  }
  
  intercept_ind <- which(((1:(h_eff * (p + 1))) %% (p + 1) == 1))
  beta_mat <- matrix(init.theta, ncol = h_eff)[-intercept_ind, ]
  
  Mn <- matrix(0, p, p)
  for (k in 1:h_eff) Mn <- Mn + beta_mat[, k, drop = FALSE] %*% t(beta_mat[, k, drop = FALSE])
  
  list(M = Mn, evalues = eigen(Mn)$values, evectors = eigen(Mn)$vectors, x = x)
}


# ============================================================================
# 2. ppl2svm
# ============================================================================
ppl2svm <- function(x, y, H, C, lambda, gamma = 3.7, penalty = "grSCAD", max.iter = 100) {
  n <- nrow(x)
  p <- ncol(x)
  
  if (!penalty %in% c("grSCAD", "grLasso", "grMCP")) stop("Invalid penalty")
  
  bar.x   <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))
  
  qprob <- seq(1/H, 1 - 1/H, by = 1/H)
  h_eff <- length(qprob)
  qy    <- quantile(y, qprob)
  U     <- sapply(qy, function(z) 2 * (z > y) - 1)
  tmp.y <- c(U)
  
  Sigma.hat      <- cov(x)
  Sigma.hat.star <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma.hat))
  
  init.theta <- rep(0, h_eff * (p + 1))
  tol        <- 1e-5
  tau.vec    <- rep(0, n * h_eff)
  
  for (iter in 1:max.iter) {
    old.theta <- init.theta
    Wtheta <- numeric(n * h_eff)
    for (k in 1:h_eff) {
      theta_k <- init.theta[((k - 1) * (p + 1) + 1):(k * (p + 1))]
      Wtheta[((k - 1) * n + 1):(k * n)] <- as.vector(x.tilde %*% theta_k)
    }
    
    u <- (1 - tau.vec) * as.integer(tmp.y - Wtheta <= 0)
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
    
    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))
    if (delta < tol) break
  }
  
  intercept_ind <- which(((1:(h_eff * (p + 1))) %% (p + 1) == 1))
  beta_mat <- matrix(init.theta, ncol = h_eff)[-intercept_ind, ]
  
  Mn <- matrix(0, p, p)
  for (k in 1:h_eff) Mn <- Mn + beta_mat[, k, drop = FALSE] %*% t(beta_mat[, k, drop = FALSE])
  
  list(M = Mn, evalues = eigen(Mn)$values, evectors = eigen(Mn)$vectors, x = x)
}


# ============================================================================
# 3. pplr
# ============================================================================
pplr <- function(x, y, H, C, lambda, gamma = 3.7, penalty = "grSCAD", max.iter = 100, ridge = 1e-10) {
  n <- nrow(x)
  p <- ncol(x)
  
  if (!penalty %in% c("grSCAD", "grLasso", "grMCP")) stop("Invalid penalty")
  
  bar.x   <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))
  
  qprob <- seq(1/H, 1 - 1/H, by = 1/H)
  h_eff <- length(qprob)
  qy    <- quantile(y, qprob)
  U     <- sapply(qy, function(z) 2 * (z > y) - 1)
  u     <- c(U)
  
  Sigma             <- cov(x)
  Sigma.tilde       <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma))
  Sigma.tilde.ridge <- Sigma.tilde + diag(ridge, p + 1)
  
  init.theta <- rep(0, h_eff * (p + 1))
  tol <- 1e-5
  omega_cum <- rep(1, n * h_eff)
  
  for (iter in 1:max.iter) {
    old.theta <- init.theta
    Wtheta <- numeric(n * h_eff)
    for (k in 1:h_eff) {
      theta_k <- init.theta[((k - 1) * (p + 1) + 1):(k * (p + 1))]
      idx_k   <- ((k - 1) * n + 1):(k * n)
      Wtheta[idx_k] <- as.vector((x.tilde * omega_cum[idx_k]) %*% theta_k)
    }
    
    pi_vec       <- 1 / (1 + exp(u * Wtheta))
    tmp_pi       <- pi_vec * (1 - pi_vec)
    sqrt_tmp_pi  <- sqrt(tmp_pi)
    
    u.tilde   <- Wtheta + (1 / tmp_pi) * (u * pi_vec)
    omega_cum <- omega_cum * sqrt_tmp_pi
    u.tilde   <- sqrt_tmp_pi * u.tilde
    
    A_blocks <- vector("list", h_eff)
    for (k in 1:h_eff) {
      idx_k       <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * omega_cum[idx_k]
      A_blocks[[k]] <- Sigma.tilde.ridge + crossprod(Wk_weighted)
    }
    
    A_sqrt_blocks <- sqrt_sym_blockdiag(A_blocks)
    
    rhs <- numeric(h_eff * (p + 1))
    for (k in 1:h_eff) {
      idx_k       <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * omega_cum[idx_k]
      rhs[((k - 1) * (p + 1) + 1):(k * (p + 1))] <- t(Wk_weighted) %*% u.tilde[idx_k]
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
    
    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))
    if (delta < tol) break
  }
  
  intercept_ind <- which(((1:(h_eff * (p + 1))) %% (p + 1) == 1))
  beta_mat <- matrix(init.theta, ncol = h_eff)[-intercept_ind, ]
  
  Mn <- matrix(0, p, p)
  for (k in 1:h_eff) Mn <- Mn + beta_mat[, k, drop = FALSE] %*% t(beta_mat[, k, drop = FALSE])
  
  list(M = Mn, evalues = eigen(Mn)$values, evectors = eigen(Mn)$vectors, x = x)
}


# ============================================================================
# 4. ppsvm
# ============================================================================
ppsvm <- function(x, y, H = 10, C = 1, lambda, gamma = 3.7, penalty = "grSCAD", max.iter = 100, ridge = 1e-10) {
  n <- nrow(x)
  p <- ncol(x)
  
  if (!penalty %in% c("grSCAD", "grLasso", "grMCP")) stop("Invalid penalty")
  
  bar.x   <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))
  
  qprob <- seq(1/H, 1 - 1/H, by = 1/H)
  h_eff <- length(qprob)
  qy    <- quantile(y, qprob)
  u     <- sapply(qy, function(z) 2 * (z > y) - 1)
  u     <- c(u)
  
  Sigma       <- cov(x)
  Sigma.tilde <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma))
  
  init.theta <- rep(0, h_eff * (p + 1))
  tol <- 1e-5
  omega_cum <- rep(1, n * h_eff)
  
  for (iter in 1:max.iter) {
    old.theta <- init.theta
    Wtheta <- numeric(n * h_eff)
    for (k in 1:h_eff) {
      theta_k <- init.theta[((k - 1) * (p + 1) + 1):(k * (p + 1))]
      idx_k   <- ((k - 1) * n + 1):(k * n)
      Wtheta[idx_k] <- as.vector((x.tilde * omega_cum[idx_k]) %*% theta_k)
    }
    
    m_t          <- u * Wtheta / n
    c_t          <- as.vector(abs(1 - m_t))
    omega_t      <- 1 / (4 * c_t)
    sqrt_omega_t <- sqrt(omega_t)
    
    u.tilde   <- (1 + c_t) * u
    omega_cum <- omega_cum * sqrt_omega_t
    u.tilde   <- sqrt_omega_t * u.tilde
    
    A_blocks <- vector("list", h_eff)
    for (k in 1:h_eff) {
      idx_k       <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * omega_cum[idx_k]
      A_blocks[[k]] <- Sigma.tilde + crossprod(Wk_weighted)
    }
    
    A_sqrt_blocks <- sqrt_sym_blockdiag(A_blocks, ridge = 1e-10)
    
    rhs <- numeric(h_eff * (p + 1))
    for (k in 1:h_eff) {
      idx_k       <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * omega_cum[idx_k]
      rhs[((k - 1) * (p + 1) + 1):(k * (p + 1))] <- t(Wk_weighted) %*% u.tilde[idx_k]
    }
    big_Y <- solve_blockdiag(A_blocks, rhs, ridge = 1e-10)
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
    
    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))
    if (delta < tol) break
  }
  
  intercept_ind <- which(((1:(h_eff * (p + 1))) %% (p + 1) == 1))
  beta_mat <- matrix(init.theta, ncol = h_eff)[-intercept_ind, ]
  
  Mn <- matrix(0, p, p)
  for (k in 1:h_eff) Mn <- Mn + beta_mat[, k, drop = FALSE] %*% t(beta_mat[, k, drop = FALSE])
  
  list(M = Mn, evalues = eigen(Mn)$values, evectors = eigen(Mn)$vectors, x = x)
}


ppqr <- function(x, y, H = 10, C = 1, lambda, gamma = 3.7, penalty = "grSCAD", max.iter = 100, ridge = 1e-10) {
  n <- nrow(x)
  p <- ncol(x)
  
  if (!penalty %in% c("grSCAD", "grLasso", "grMCP")) stop("Invalid penalty")
  
  bar.x   <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))
  
  qprob <- seq(1/H, 1 - 1/H, by = 1/H)
  h_eff <- length(qprob)
  tmp.y <- rep(y, times = h_eff)
  
  Sigma             <- cov(x)
  Sigma.tilde       <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma))
  Sigma.tilde.ridge <- Sigma.tilde + diag(ridge, p + 1)
  
  init.theta <- rep(0, h_eff * (p + 1))
  tau_vec    <- rep(qprob, each = n)
  epsilon    <- rep(1 + 0.0032 * qprob, each = n)
  tol        <- 1e-5
  omega_cum  <- rep(1, n * h_eff)
  
  # 수치적 에러 방지: r_t가 0일 때의 NaN을 막기 위해 수학적으로 약분된 상수 형태로 사용
  K_t_const <- 1 + 4 * tau_vec * (1 - tau_vec)
  
  for (iter in 1:max.iter) {
    old.theta <- init.theta
    Wtheta <- numeric(n * h_eff)
    for (k in 1:h_eff) {
      theta_k <- init.theta[((k - 1) * (p + 1) + 1):(k * (p + 1))]
      idx_k   <- ((k - 1) * n + 1):(k * n)
      Wtheta[idx_k] <- as.vector((x.tilde * omega_cum[idx_k]) %*% theta_k)
    }
    
    r_t          <- as.vector(tmp.y - Wtheta)
    
    # K_t를 r_t^2로 나누지 않음으로써 0/0 (NaN) 에러 완벽 차단
    omega_t      <- K_t_const / (4 * (abs(r_t) + epsilon))
    sqrt_omega_t <- sqrt(omega_t)
    
    u.tilde   <- tmp.y + abs(r_t) * (2 * tau_vec - 1)
    omega_cum <- omega_cum * sqrt_omega_t
    u.tilde   <- sqrt_omega_t * u.tilde
    
    A_blocks <- vector("list", h_eff)
    for (k in 1:h_eff) {
      idx_k       <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * omega_cum[idx_k]
      A_blocks[[k]] <- Sigma.tilde.ridge + crossprod(Wk_weighted)
    }
    
    A_sqrt_blocks <- sqrt_sym_blockdiag(A_blocks)
    
    rhs <- numeric(h_eff * (p + 1))
    for (k in 1:h_eff) {
      idx_k       <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * omega_cum[idx_k]
      rhs[((k - 1) * (p + 1) + 1):(k * (p + 1))] <- t(Wk_weighted) %*% u.tilde[idx_k]
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
    
    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))
    if (delta < tol) break
  }
  
  # 인덱싱 에러 수정: 행렬 변환 후에는 절편이 무조건 1번째 행에 위치합니다.
  beta_mat <- matrix(init.theta, ncol = h_eff)[-1, , drop = FALSE]
  
  Mn <- matrix(0, p, p)
  for (k in 1:h_eff) Mn <- Mn + beta_mat[, k, drop = FALSE] %*% t(beta_mat[, k, drop = FALSE])
  
  list(M = Mn, evalues = eigen(Mn)$values, evectors = eigen(Mn)$vectors, x = x)
}
# ============================================================================
# 6. pplssvm (Newly Created Block-Diagonal Accelerated + 3 Penalties)
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
    if (delta < tol) break
  }
  
  intercept_ind <- which(((1:(h_eff * (p + 1))) %% (p + 1) == 1))
  beta_mat <- matrix(init.theta, ncol = h_eff)[-intercept_ind, ]
  
  Mn <- matrix(0, p, p)
  for (k in 1:h_eff) Mn <- Mn + beta_mat[, k, drop = FALSE] %*% t(beta_mat[, k, drop = FALSE])
  
  list(M = Mn, evalues = eigen(Mn)$values, evectors = eigen(Mn)$vectors, x = x)
}


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
  
  # v1의 치명적 버그(chol2inv)를 수정하고 수학적으로 올바른 G^{-1} W^T u 도출
  xi <- solve_blockdiag(A_sqrt_blocks, rhs)
  xi.tilde <- as.vector(xi) / C
  
  group <- rep(1:(p + 1), times = h_eff)
  colnames(G.tilde.big) <- group
  
  obj_grpreg <- grpreg::grpreg(X = G.tilde.big, y = xi.tilde, group = group, penalty = penalty, family = "gaussian",
                       lambda = lambda, alpha = 1, eps = 1e-5, max.iter = max.iter, dfmax = p,
                       gmax = length(unique(group)), gamma = gamma)
  
  theta <- obj_grpreg$beta[-1]
  beta_mat <- matrix(theta, ncol = h_eff)[-1, , drop = FALSE] # 절편 제거 완벽 대응
  
  Mn <- matrix(0, p, p)
  for (k in 1:h_eff) Mn <- Mn + beta_mat[, k, drop = FALSE] %*% t(beta_mat[, k, drop = FALSE])
  
  list(M = Mn, evalues = eigen(Mn)$values, evectors = eigen(Mn)$vectors, x = x)
}

# ============================================================================
# 2. ppwlr (Penalized Principal Weighted Logistic Regression) - FINAL FIX
# ============================================================================
ppwlr <- function(x, y, H, C, lambda, gamma = 3.7,
                     penalty = "grSCAD", max.iter = 100, ridge = 1e-10) {
  
  n <- nrow(x); p <- ncol(x)
  bar.x   <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))
  
  qprob   <- seq(1/H, 1 - 1/H, by = 1/H)
  h_eff   <- length(qprob)
  pi.grid <- qprob
  
  u <- rep(y, h_eff)
  
  Omega_mat <- matrix(
    ((u + 1) / 2) * rep(pi.grid, each = n) +
      ((-u + 1) / 2) * rep(1 - pi.grid, each = n),
    ncol = h_eff
  )
  Omega_vec <- as.vector(Omega_mat)
  
  Sigma       <- cov(x)
  Sigma.tilde <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma))
  
  # init at 0 (lambda_max warm-start point)
  init.theta <- rep(0, h_eff * (p + 1))
  tol <- 1e-5
  
  # Penalty-dependent S scaling factor (see header note 3)
  S_scale <- if (penalty == "grSCAD") n / C else 1
  
  for (iter in 1:max.iter) {
    old.theta <- init.theta
    
    # W * theta (W is unweighted block-diagonal of x.tilde)
    Wtheta <- numeric(n * h_eff)
    for (k in 1:h_eff) {
      theta_k <- init.theta[((k - 1) * (p + 1) + 1):(k * (p + 1))]
      idx_k   <- ((k - 1) * n + 1):(k * n)
      Wtheta[idx_k] <- as.vector(x.tilde %*% theta_k)
    }
    
    # IRLS weights
    pi_vec      <- 1 / (1 + exp(u * Wtheta))
    tmp_pi      <- pi_vec * (1 - pi_vec) * Omega_vec
    sqrt_tmp_pi <- sqrt(tmp_pi)
    
    # FIX 2: penalty-dependent u.tilde formula
    if (penalty == "grSCAD") {
      # grSCAD: u.tilde = W*theta * Omega_vec + (1/tmp_pi) * (u * pi)
      u.tilde <- Wtheta * Omega_vec + (1 / tmp_pi) * (u * pi_vec)
    } else {
      # grLasso, grMCP: u.tilde = W*theta + (1/tmp_pi) * (u * pi * Omega_vec)
      u.tilde <- Wtheta + (1 / tmp_pi) * (u * pi_vec * Omega_vec)
    }
    u.tilde_new <- sqrt_tmp_pi * u.tilde
    
    # Per-slice A blocks with FIX 3: penalty-dependent S scaling
    A_blocks <- vector("list", h_eff)
    for (k in 1:h_eff) {
      idx_k       <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * sqrt_tmp_pi[idx_k]
      A_blocks[[k]] <- S_scale * Sigma.tilde + crossprod(Wk_weighted)
    }
    
    A_sqrt_blocks <- sqrt_sym_blockdiag(A_blocks, ridge = ridge)
    
    rhs <- numeric(h_eff * (p + 1))
    for (k in 1:h_eff) {
      idx_k       <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * sqrt_tmp_pi[idx_k]
      rhs[((k - 1) * (p + 1) + 1):(k * (p + 1))] <-
        t(Wk_weighted) %*% u.tilde_new[idx_k]
    }
    
    big_Y <- solve_blockdiag(A_blocks, rhs, ridge = ridge)
    res   <- big_Y - multiply_blockdiag(A_sqrt_blocks, init.theta)
    
    group   <- rep(1:(p + 1), times = h_eff)
    J       <- p + 1
    
    for (j in 1:J) {
      ind    <- which(group == j)
      G_cols <- get_var_columns(A_sqrt_blocks, j)
      
      z_j    <- (C / n) * (t(G_cols) %*% res) + init.theta[ind]
      z_norm <- norm(z_j, "2")
      
      if (z_norm < 1e-12) {
        theta_j_new <- rep(0, length(ind))
      } else {
        # FIX 4: penalty-specific thresholding matching original exactly
        if (penalty == "grSCAD") {
          theta_j_new <- scad_firm_thresholding(z_norm,
                                                lambda1 = lambda,
                                                lambda2 = 2,
                                                gamma = gamma) * z_j / z_norm
        } else if (penalty == "grLasso") {
          theta_j_new <- soft_thresholding(z_norm, lambda = lambda) /
            (1 + 2) * z_j / z_norm
        } else if (penalty == "grMCP") {
          # Note: original uses lambda2=1 here (NOT 2)
          theta_j_new <- firm_thresholding(z_norm, lambda1 = lambda,
                                           lambda2 = 1, gamma = gamma) *
            z_j / z_norm
        } else {
          stop(sprintf("Unknown penalty: %s", penalty))
        }
        theta_j_new <- as.vector(theta_j_new)
      }
      
      res <- res - G_cols %*% (theta_j_new - init.theta[ind])
      init.theta[ind] <- theta_j_new
      res <- as.vector(res)
    }
    
    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))
    if (delta < tol) break
  }
  
  intercept_ind <- which(((1:(h_eff * (p + 1))) %% (p + 1) == 1))
  beta_mat <- matrix(init.theta, ncol = h_eff)[-intercept_ind, ]
  
  Mn <- matrix(0, p, p)
  for (k in 1:h_eff) {
    Mn <- Mn + beta_mat[, k, drop = FALSE] %*% t(beta_mat[, k, drop = FALSE])
  }
  
  list(M = Mn, evalues = eigen(Mn)$values, evectors = eigen(Mn)$vectors, x = x)
}

# ============================================================================
# 3. ppwl2svm (Penalized Principal Weighted L2 Hinge SVM)
# ============================================================================
ppwl2svm <- function(x, y, H = 10, C = 1, lambda, gamma = 3.7, penalty = "grSCAD", max.iter = 100, ridge = 1e-4) {
  n <- nrow(x); p <- ncol(x)
  bar.x   <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))
  
  qprob <- seq(1/H, 1 - 1/H, by = 1/H)
  h_eff <- length(qprob)
  tmp.y <- rep(y, times = h_eff)
  pi.grid <- rep(qprob, each = n)
  
  weight <- (1 - pi.grid) * as.numeric(tmp.y == 1) + pi.grid * as.numeric(tmp.y == -1)
  sqrt_weight <- sqrt(weight)
  weighted_y <- sqrt_weight * tmp.y
  
  Sigma <- cov(x)
  Sigma.hat.star <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma))
  init.theta <- rep(0, h_eff * (p + 1))
  tol        <- 1e-5
  
  for (iter in 1:max.iter) {
    old.theta <- init.theta
    Wtheta <- numeric(n * h_eff)
    for (k in 1:h_eff) {
      theta_k <- init.theta[((k - 1) * (p + 1) + 1):(k * (p + 1))]
      idx_k   <- ((k - 1) * n + 1):(k * n)
      Wtheta[idx_k] <- as.vector((x.tilde * sqrt_weight[idx_k]) %*% theta_k)
    }
    
    u_mask <- as.integer(weighted_y - Wtheta <= 0)
    y.tilde <- weighted_y * sqrt(u_mask)
    
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
    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))
    if (delta < tol) break
  }
  
  beta_mat <- matrix(init.theta, ncol = h_eff)[-1, , drop = FALSE]
  Mn <- matrix(0, p, p)
  for (k in 1:h_eff) Mn <- Mn + beta_mat[, k, drop = FALSE] %*% t(beta_mat[, k, drop = FALSE])
  list(M = Mn, evalues = eigen(Mn)$values, evectors = eigen(Mn)$vectors, x = x)
}

# ============================================================================
# 4. ppwsvm (Penalized Principal Weighted Hinge SVM)
# ============================================================================
ppwsvm <- function(x, y, H = 10, C = 1, lambda, gamma = 3.7, penalty = "grSCAD", max.iter = 100, ridge = 1e-10) {
  n <- nrow(x); p <- ncol(x)
  bar.x   <- apply(x, 2, mean)
  x.tilde <- cbind(rep(1, n), t(t(x) - bar.x))
  
  qprob <- seq(1/H, 1 - 1/H, by = 1/H)
  h_eff <- length(qprob)
  tmp.y <- rep(y, times = h_eff)
  pi.grid <- rep(qprob, each = n)
  
  Sigma <- cov(x)
  Sigma.tilde <- cbind(rep(0, p + 1), rbind(rep(0, p), Sigma))
  init.theta <- rep(0, h_eff * (p + 1))
  tol <- 1e-5
  
  weight_t_base <- (1 - pi.grid) * as.numeric(tmp.y == 1) + pi.grid * as.numeric(tmp.y == -1)
  omega_cum <- rep(1, n * h_eff)
  
  for (iter in 1:max.iter) {
    old.theta <- init.theta
    Wtheta <- numeric(n * h_eff)
    for (k in 1:h_eff) {
      theta_k <- init.theta[((k - 1) * (p + 1) + 1):(k * (p + 1))]
      idx_k   <- ((k - 1) * n + 1):(k * n)
      Wtheta[idx_k] <- as.vector((x.tilde * omega_cum[idx_k]) %*% theta_k)
    }
    
    m_t <- as.vector(tmp.y * Wtheta) / n
    c_t <- abs(1 - m_t)
    omega_t <- weight_t_base / (4 * c_t)
    sqrt_omega_t <- sqrt(omega_t)
    
    u.tilde <- (1 + c_t) * tmp.y
    u.tilde_new <- sqrt_omega_t * u.tilde
    omega_cum <- omega_cum * sqrt_omega_t
    
    A_blocks <- vector("list", h_eff)
    for (k in 1:h_eff) {
      idx_k <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * omega_cum[idx_k]
      A_blocks[[k]] <- Sigma.tilde + crossprod(Wk_weighted)
    }
    
    A_sqrt_blocks <- sqrt_sym_blockdiag(A_blocks, ridge = ridge)
    rhs <- numeric(h_eff * (p + 1))
    for (k in 1:h_eff) {
      idx_k <- ((k - 1) * n + 1):(k * n)
      Wk_weighted <- x.tilde * omega_cum[idx_k]
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
    delta <- max(abs(init.theta - old.theta)) / (1 + max(abs(old.theta)))
    if (delta < tol) break
  }
  
  beta_mat <- matrix(init.theta, ncol = h_eff)[-1, , drop = FALSE]
  Mn <- matrix(0, p, p)
  for (k in 1:h_eff) Mn <- Mn + beta_mat[, k, drop = FALSE] %*% t(beta_mat[, k, drop = FALSE])
  list(M = Mn, evalues = eigen(Mn)$values, evectors = eigen(Mn)$vectors, x = x)
}
