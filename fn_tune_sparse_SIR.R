tune.sSIR <- function(x, y, d, H, lambda.grid, n.fold)
{
  require(energy)
  n <- length(y)
  p <- ncol(x)
  cat(n.fold, "-fold CV.\n", sep = "")
  L <- length(lambda.grid)
  value <- matrix(0, n.fold, L)
  
  for (ll in 1:L){
    for (jj in 1:n.fold){
      new.index <- sample(n)
      step <- round(n/n.fold)
      
      ts.id <- ((jj-1)*step+1):min(n,((jj)*step))
      tr.id <- unlist(setdiff(1:n, ts.id))
      
      obj <- sparseSIR_MAX_SCAD(x[tr.id,], y[tr.id], H, lambda=lambda.grid[ll])[,1:d,drop = F]

      value[jj,ll] <- energy::dcor(y[ts.id,drop = F], x[ts.id,,drop = F] %*% obj)
    }
  }
  tmp <- apply(value, 2, mean)
  names(tmp) <- lambda.grid
  
  sel <- which.max(tmp)

  list(opt.lambda = lambda.grid[sel], dcor = round(tmp,3))
}

tune.sSIR.L2 <- function(x, y, d, H, lambda.grid, n.fold)
{
  require(energy)
  n <- length(y)
  p <- ncol(x)
  cat(n.fold, "-fold CV.\n", sep = "")
  L <- length(lambda.grid)
  value <- matrix(0, n.fold, L)
  
  for (ll in 1:L){
    for (jj in 1:n.fold){
      new.index <- sample(n)
      step <- round(n/n.fold)
      
      ts.id <- ((jj-1)*step+1):min(n,((jj)*step))
      tr.id <- unlist(setdiff(1:n, ts.id))
      
      obj <- L2sparseSIR(x[tr.id,], y[tr.id], H, lambda=lambda.grid[ll])[,1:d,drop = F]
      
      value[jj,ll] <- energy::dcor(y[ts.id,drop = F], x[ts.id,,drop = F] %*% obj)

    }
  }
  tmp <- apply(value, 2, mean)
  names(tmp) <- lambda.grid
  
  sel <- which.max(tmp)

  list(opt.lambda = lambda.grid[sel], dcor = round(tmp,3))
}

