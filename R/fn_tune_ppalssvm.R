##cross-validataion for Group coordinate descent : Sparse PALS GCD
##function arguments##
# x:data matrix
# y:response (vector), can be either continuous or binary (+1,-1) type
# H:number of k for loss function
# d:true structural dimension
# C:cost parameter, usually set to 1
# lambda grid: if not specified, the function automatically choose a proper lambda via grid search
# gamma: nonconvex shape parameter
# penalty: type of penalty function: 'grSCAD', 'grLasso', 'grMCP'. 
# max.iter: maximum number of iteration
# n.fold: number of cross validation
# it returns optimal 'lambda' which maximizes distance correlation


tune_ppalssvm <- function(x, y, d, H, C, lambda.grid, gamma, penalty, max.iter, n.fold)
{
  require(grpreg)
  require(Matrix)
  require(expm)
  require(blockmatrix)
  require(Rfast)
  n <- length(y)
  
  cat(n.fold, "-fold CV.\n", sep = "")
  L <- length(lambda.grid)
  value <- matrix(0, n.fold, L)
  
  for (ll in 1:L){
    for (jj in 1:n.fold){
      new.index <- sample(n)
      step <- round(n/n.fold)
      
      ts.id <- ((jj-1)*step+1):min(((jj)*step), n)
      tr.id <- unlist(setdiff(1:n, ts.id))
      
      #G.C.D.
      obj <- GCD_SPALS(x=x[tr.id,,drop = F], y=y[tr.id,drop = F], H=H, C=C, lambda=lambda.grid[ll], gamma, penalty, max.iter)$vectors[,1:d,drop = F]
      
      value[jj,ll] <- dcor(y[ts.id,drop = F], x[ts.id,,drop = F] %*% obj)$dcor
    }
  }
  tmp <- apply(value, 2, mean)
  names(tmp) <- lambda.grid
  
  sel <- which.max(tmp)

  list(opt.lambda = lambda.grid[sel], dcor = round(tmp,3))
}



