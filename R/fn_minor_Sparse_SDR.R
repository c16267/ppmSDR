# Soft-thresholding operator
soft_thresholding <- function(z, lambda) {
  if (z > lambda) return(z - lambda)
  if (z < -lambda) return(z + lambda)
  return(0)
}

# Firm-thresholding operator
firm_thresholding <- function(z, lambda1, lambda2, gamma) {
  s <- 0
  if (z > 0) s <- 1
  else if (z < 0) s <- -1
  if (abs(z) <= lambda1) return(0)
  else if (abs(z) <= gamma*lambda1*(1+lambda2)) return(s*(abs(z)-lambda1)/(1+lambda2-1/gamma))
  else return(z/(1+lambda2))
}


# SCAD-modified firm-thresholding operator
scad_firm_thresholding <- function(z, lambda1, lambda2, gamma) {
  s <- 0
  if(z > 0){
    s <- 1
  }else if(z < 0){
    s <- -1
  }else{
    s <- 0
  }

  if (abs(z) <= lambda1){
    return(0)
  }else if (abs(z) <= (lambda1*(1+lambda2)+lambda1)){
    return(s*(abs(z)-lambda1)/(1+lambda2))
  }else if (abs(z) <= gamma*lambda1*(1+lambda2)){
    return(s*(abs(z)-gamma*lambda1/(gamma-1))/(1-1/(gamma-1)+lambda2))
  }else{
    return(z/(1+lambda2))
  }
}


norm <- function(x, p) {
  return(sqrt(sum(x^2)))
}


fn_orthonormalization <- function(X, y, g=NULL, m=NULL, bilevel=FALSE, family="gaussian") {

  # 1. Validate and coerce X into a numeric matrix
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (storage.mode(X) == "integer") storage.mode(X) <- "double"

  # 2. Check for missing values in the predictor matrix
  if (any(is.na(X))) {
    stop("Missing data (NA's) detected in X. You must eliminate missing data before passing X to ppmSDR.", call.=FALSE)
  }

  # 3. Handle undefined group assignments (default: each feature is treated as its own independent group)
  if (is.null(g)) {
    g <- 1:ncol(X)
  }

  # 4. Verify that the length of the group vector matches the number of columns in X
  if (length(g) != ncol(X)) {
    stop("Dimensions of group is not compatible with X", call.=FALSE)
  }

  # Extract column names, or assign default names (V1, V2, ...) if they are missing
  xnames <- if (is.null(colnames(X))) paste("V", 1:ncol(X), sep="") else colnames(X)

  # Convert group assignments into contiguous integers
  gf <- factor(g)
  g <- as.integer(gf)
  J <- max(g)
  lev <- levels(gf)

  # 5. Initialize group multipliers (m) with NA if not explicitly provided
  if (is.null(m)) {
    m <- rep(NA, length(lev))
  }

  # Create a structure containing group information, factor levels, and multipliers
  G <- structure(g, levels=lev, m=m)

  # 6. Feature-level standardization (Calls the C function for computational efficiency)
  std <- .Call("standardize", X, PACKAGE = "ppmSDR")
  XX <- std[[1]]      # Standardized matrix
  center <- std[[2]]  # Column means
  scale <- std[[3]]   # Column standard deviations
  nz <- which(scale > 1e-6) # Indices of non-constant columns (removes zero-variance features)

  # 7. Reorder columns so that members of the same group are adjacent to each other
  G <- reorderG(G, attr(G, 'm'), bilevel)
  if (attr(G, 'reorder')) XX <- XX[, attr(G, 'ord')]

  # 8. Group-level standardization (Orthonormalization)
  if (!bilevel) {
    XX <- orthogonalize(XX, G)
    g <- attr(XX, "group") # Update group assignments after orthonormalization
  }

  # 9. Set default group multipliers (e.g., square root of the group size) if they were initially missing
  m_attr <- attr(G, 'm')
  if (all(is.na(m_attr))) {
    m <- if (bilevel) rep(1, max(g)) else sqrt(table(g[g!=0]))
  } else {
    m <- m_attr
  }

  # 10. Center the response vector (y) if the family is gaussian
  if (family == "gaussian") {
    meanY <- mean(y)
    y <- y - meanY
    attr(y, "mean") <- meanY
  }

  attr(y, "m") <- 1

  # Return a comprehensive list containing all processed elements required for downstream modeling
  return(list(X=XX, g=g, m=m, reorder=attr(G, 'reorder'), ord.inv=attr(G, 'ord.inv'),
              names=xnames, center=center, scale=scale, nz=nz, "y.std"=y))
}


# fn_orthonormalization <- function(X, y, g=NULL, m, ncolY=dim(y)[2], bilevel=FALSE, family="gaussian") {
#
#   if (!inherits(X, "matrix")) {
#     tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
#     if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
#   }
#   if (storage.mode(X)=="integer") storage.mode(X) <- "double"
#   if (any(is.na(X))) stop("Missing data (NA's) detected in X.  You must eliminate missing data (e.g., by removing cases, removing features, or imputation) before passing X to grpreg", call.=FALSE)
#   if (length(g) != ncol(X)) stop ("Dimensions of group is not compatible with X", call.=FALSE)
#   xnames <- if (is.null(colnames(X))) paste("V", 1:ncol(X), sep="") else colnames(X)
#
#
#   gf <- factor(g)
#   g <- as.integer(gf)
#   J <- max(g)
#   lev <- levels(gf)
#   m <- rep(NA, length(lev))
#   G <- structure(g, levels=lev, m=m)
#
#   # Feature-level standardization
#   std <- .Call("standardize", X, PACKAGE = "ppmSDR")
#   XX <- std[[1]]
#   center <- std[[2]]  #colmean
#   scale <- std[[3]]   #col standard error
#   nz <- which(scale > 1e-6)                # non-constant columns
#
#   # Reorder groups, if necessary
#   G <- reorderG(G, attr(G, 'm'), bilevel)
#   if (attr(G, 'reorder')) XX <- XX[, attr(G, 'ord')]
#
#
#   # Group-level standardization
#   if (!bilevel) {
#     XX <- orthogonalize(XX, G)
#     g <- attr(XX, "group")
#   }
#
#   # Set group multiplier if missing
#   m <- attr(G, 'm')
#   if (all(is.na(m))) {
#     m <- if (bilevel) rep(1, max(g)) else sqrt(table(g[g!=0]))
#   }
#
#   # centering y
#   if (family=="gaussian") {
#     meanY <- mean(y)
#     y <- y - meanY
#     attr(y, "mean") <- meanY
#   }
#
#   attr(y, "m") <- 1
#   #rslt <- list("orthogonalized.x"=X.tilde, "QL"=T, "group"=g, "center"=center, "scale"=scale, "nz"=nz, "y.std"=y)
#   return(list(X=XX, g=g, m=m, reorder=attr(G, 'reorder'), ord.inv=attr(G, 'ord.inv'), names=xnames,
#               center=center, scale=scale, nz=nz, "y.std"=y))
# }
#



orthogonalize <- function(X, group) {
  n <- nrow(X)
  J <- max(group)
  T <- vector("list", J)
  XX <- matrix(0, nrow=nrow(X), ncol=ncol(X))
  XX[, which(group==0)] <- X[, which(group==0)]
  for (j in seq_along(integer(J))) {
    ind <- which(group==j)
    if (length(ind)==0) next
    SVD <- svd(X[, ind, drop=FALSE], nu=0)
    r <- which(SVD$d > 1e-10)
    T[[j]] <- sweep(SVD$v[, r, drop=FALSE], 2, sqrt(n)/SVD$d[r], "*")
    XX[, ind[r]] <- X[, ind] %*% T[[j]]
  }
  nz <- !apply(XX==0, 2, all)
  XX <- XX[, nz, drop=FALSE]
  attr(XX, "orthogonalized.x") <- XX
  attr(XX, "T") <- T
  attr(XX, "group") <- group[nz]
  XX
}

unorthogonalize <- function(b, XX, group, intercept=TRUE) {
  ind <- !sapply(attr(XX, "T"), is.null)
  T <- Matrix::bdiag(attr(XX, "T")[ind])

  if (intercept) {
    ind0 <- c(1, 1+which(group==0))
    val <- as.matrix(rbind(b[ind0, , drop=FALSE], T %*% b[-ind0, , drop=FALSE]))
  } else if (sum(group==0)) {
    ind0 <- which(group==0)
    val <- as.matrix(rbind(b[ind0, , drop=FALSE], T %*% b[-ind0, , drop=FALSE]))
  } else {
    val <- as.matrix(T %*% b)
  }
  return(val)
}


my.unorthogonalize <- function(b, XX, group, intercept=TRUE) {
  ind <- !sapply(attr(XX, "T"), is.null)
  T <- Matrix::bdiag(attr(XX, "T")[ind])

  if (intercept) {
    ind0 <- c(1, 1+which(group==0))
    val <- as.matrix(rbind(b[ind0, , drop=FALSE], T %*% b[-ind0, , drop=FALSE]))
  } else if (sum(group==0)) {
    ind0 <- which(group==0)
    val <- as.matrix(rbind(b[ind0, , drop=FALSE], T %*% b[-ind0, , drop=FALSE]))
  } else {
    val <- as.matrix(T %*% b)
  }
  return(val)
}

my.unstandardize <- function(b, XG) {
  beta <- rep(0, length(XG$scale))
  beta[XG$nz] <- b / XG$scale[XG$nz]
  beta
}

setupG <- function(group, m, bilevel) {
  gf <- factor(group)
  if (any(levels(gf)=='0')) {
    g <- as.integer(gf) - 1
    lev <- levels(gf)[levels(gf)!='0']
  } else {
    g <- as.integer(gf)
    lev <- levels(gf)
  }
  if (is.numeric(group) | is.integer(group)) {
    lev <- paste0("G", lev)
  }
  if (missing(m)) {
    m <- rep(NA, length(lev))
    names(m) <- lev
  } else {
    #if (all.equal(sort(names(m)), sort(group)))
    TRY <- try(as.integer(group)==g)
    if (inherits(TRY, 'try-error') || any(!TRY)) stop('Attempting to set group.multiplier is ambiguous if group is not a factor', call.=FALSE)
    if (length(m) != length(lev)) stop("Length of group.multiplier must equal number of penalized groups", call.=FALSE)
    if (storage.mode(m) != "double") storage.mode(m) <- "double"
    if (any(m < 0)) stop('group.multiplier cannot be negative', call.=FALSE)
  }
  structure(g, levels=lev, m=m)
}

subsetG <- function(g, nz) {
  lev <- attr(g, 'levels')
  m <- attr(g, 'm')
  new <- g[nz]
  dropped <- setdiff(g, new)
  if (length(dropped)) {
    lev <- lev[-dropped]
    m <- m[-dropped]
    gf <- factor(new)
    new <- as.integer(gf) - 1*any(levels(gf)=='0')
  }
  structure(new, levels=lev, m=m)
}

reorderG <- function(g, m, bilevel) {
  og <- g
  lev <- attr(g, 'levels')
  m <- attr(g, 'm')
  if (any(g==0)) {
    g <- as.integer(relevel(factor(g), "0"))-1
  }
  if (any(order(g) != 1:length(g))) {
    reorder <- TRUE
    gf <- factor(g)
    if (any(levels(gf)=="0")) {
      gf <- relevel(gf, "0")
      g <- as.integer(gf) - 1
    } else {
      g <- as.integer(gf)
    }
    ord <- order(g)
    ord.inv <- match(1:length(g), ord)
    g <- g[ord]
  } else {
    reorder <- FALSE
    ord <- ord.inv <- NULL
  }
  structure(g, levels=lev, m=m, ord=ord, ord.inv=ord.inv, reorder=reorder)
}

multiX <- function(X, m) {
  p <- ncol(X)
  n <- nrow(X)
  A <- matrix(0, m*n, m*p)
  for (i in 1:m) {
    A[m*(1:n)-i+1, m*(1:p)-i+1] <- X
  }
  cbind(matrix(as.double(diag(m)), m*n, m, byrow=TRUE)[,2:m], A)
}

multiG <- function(g, ncolY) {
  structure(c(rep(0, ncolY-1), rep(g, each=ncolY)),
            levels=attr(g, 'levels'),
            m=attr(g, 'm'))
}



newXG <- function(X, g, m, ncolY, bilevel) {
  # Coerce X to matrix
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (any(is.na(X))) stop("Missing data (NA's) detected in X.  You must eliminate missing data (e.g., by removing cases, removing features, or imputation) before passing X to grpreg", call.=FALSE)
  if (length(g) != ncol(X)) stop ("Dimensions of group is not compatible with X", call.=FALSE)
  xnames <- if (is.null(colnames(X))) paste("V", 1:ncol(X), sep="") else colnames(X)

  # Setup group
  G <- setupG(g, m, bilevel)

  # Reconfigure for multiple outcomes, if necessary
  if (ncolY > 1) {
    X <- multiX(X, ncolY)
    G <- multiG(G, ncolY)
  }

  # Feature-level standardization
  std <- .Call("standardize", X, PACKAGE = "ppmSDR")
  XX <- std[[1]]
  center <- std[[2]]  #colmean
  scale <- std[[3]]   #col standard error
  nz <- which(scale > 1e-6)                # non-constant columns
  if (length(nz) != ncol(X)) {
    XX <- XX[, nz, drop=FALSE]
    G <- subsetG(G, nz)
  }

  # Reorder groups, if necessary
  G <- reorderG(G, attr(G, 'm'), bilevel)
  if (attr(G, 'reorder')) XX <- XX[, attr(G, 'ord')]

  # Group-level standardization
  if (!bilevel) {
    XX <- orthogonalize(XX, G)
    g <- attr(XX, "group")
  } else {
    g <- as.integer(G)
  }

  # Set group multiplier if missing
  m <- attr(G, 'm')
  if (all(is.na(m))) {
    m <- if (bilevel) rep(1, max(g)) else sqrt(table(g[g!=0]))
  }

  # Return
  return(list(X=XX, g=g, m=m, reorder=attr(G, 'reorder'), ord.inv=attr(G, 'ord.inv'), names=xnames,
              center=center, scale=scale, nz=nz))
}

# centering y
newY <- function(y, family) {
  if (is.data.frame(y)) y <- as.matrix(y)
  if (is.matrix(y)) {
    d <- dim(y)
    y <- t(y)
  } else {
    d <- c(length(y), 1)
  }

  # Convert fuzzy binomial data
  if (family=="binomial" && typeof(y) != "logical") {
    tab <- table(y)
    if (length(tab) > 2) stop("Attemping to use family='binomial' with non-binary data", call.=FALSE)
    if (!identical(names(tab), c("0", "1"))) {
      message(paste0("Logistic regression modeling Pr(y=", names(tab)[2], ")"))
      y <- as.double(as.character(y) == names(tab)[2])
      if (d[2] > 1) attr(y, "dim") <- d
    }
  }

  # Convert to double, if necessary
  if (typeof(y) != "double") {
    tryCatch(storage.mode(y) <- "double", warning=function(w) {stop("y must be numeric or able to be coerced to numeric", call.=FALSE)})
  }
  if (any(is.na(y))) stop("Missing data (NA's) detected in outcome y.  You must eliminate missing data (e.g., by removing cases or imputation) before passing y to grpreg", call.=FALSE)

  # Handle multi
  if (is.matrix(y)) {
    if (ncol(y) > 1) {
      if (is.null(colnames(y))) paste("Y", 1:ncol(y), sep="")
    }
    attributes(y) <- NULL
  }

  if (family=="gaussian") {
    meanY <- mean(y)
    y <- y - meanY
    attr(y, "mean") <- meanY
  }
  attr(y, "m") <- d[2]
  y
}

setupLambda <- function(X, y, group, family, penalty, alpha, lambda.min, log.lambda, nlambda, group.multiplier) {

  # Fit to unpenalized covariates
  n <- length(y)
  K <- table(group)
  K1 <- if (min(group)==0) cumsum(K) else c(0, cumsum(K))
  storage.mode(K1) <- "integer"
  if (K1[1]!=0) {
    fit <- glm(y~X[, group==0], family=family)
  } else {
    fit <- glm(y~1, family=family)
  }

  ## Determine lambda.max
  if (family=="gaussian") {
    r <- fit$residuals
  } else {
    w <- fit$weights
    if (max(w) < 1e-4) stop("Unpenalized portion of model is already saturated; exiting...", call.=FALSE)
    r <- residuals(fit, "working")*w
  }
  if (strtrim(penalty, 2) == "gr") {
    zmax <- .Call("maxgrad", X, r, K1, as.double(group.multiplier), PACKAGE = "ppmSDR") / n
  } else {
    zmax <- .Call("maxprod", X, r, K1, as.double(group.multiplier), PACKAGE = "ppmSDR") / n
  }
  #lambda.max <- zmax/alpha
  lambda.max <- 3

  # if (log.lambda) { # lambda sequence on log-scale
  #   if (lambda.min==0) {
  #     lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)), 0)
  #   } else {
  #     lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), length=nlambda))
  #   }
  # } else { # lambda sequence on linear-scale
  #   if (lambda.min==0) {
  #     lambda <- c(seq(lambda.max, 0.001*lambda.max, length = nlambda-1), 0)
  #   } else {
  #     lambda <- seq(lambda.max, lambda.min*lambda.max, length = nlambda)
  #   }
  # }
  if (lambda.min==0) {
   lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)), 0)
  } else {
   lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), length=nlambda))
  }
  lambda
}

setupLambda.gBridge <- function(X, y, group, family, alpha, lambda.min, lambda.max, nlambda, gamma, group.multiplier) {
  ## Fit to unpenalized covariates
  n <- length(y)
  ind <- which(group!=0)
  if (length(ind)!=length(group)) {
    fit <- glm(y~X[, group==0], family=family)
  } else {
    fit <- glm(y~1, family=family)
  }

  ## Guess lambda.max
  if (missing(lambda.max)) {
    if (family=="gaussian") {
      z <- crossprod(X[, ind], fit$residuals) / n
      a <- .35
    } else {
      z <- crossprod(X[, ind], fit$weights * residuals(fit, "working")) / n
      a <- .2
    }
    lambda.max <- max(abs(z)/group.multiplier[group])*a^(1-gamma)/(gamma*alpha)
  }
  if (lambda.min==0) {
    lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), len=nlambda-1)), 0)
  } else {
    lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), len=nlambda))
  }
  return(rev(lambda))
}

unstandardize <- function(b, XG) {
  beta <- matrix(0, nrow=1+length(XG$scale), ncol=ncol(b))
  beta[1 + XG$nz,] <- b[-1,] / XG$scale[XG$nz]
  beta[1,] <- b[1,] - crossprod(XG$center, beta[-1, , drop=FALSE])
  beta
}



d2 <- function(Bhat, B) {
  # This function reports performance in terms of the Frobenius norm of the projection matrix difference

  return(norm(B%*%solve(t(B)%*%B)%*%t(B)-Bhat%*%solve(t(Bhat)%*%Bhat)%*%t(Bhat),"f"))
}

p2 <- function(theta, lambda)
{
  value <- rep(0, length(theta))
  if (lambda != 0){
    temp1 <- (3.7 * lambda - theta)
    temp2 <- temp1 * (temp1 > 0)

    denom <- ((3.7-1)*lambda)
    value <- lambda * (1 - temp2/denom) * (theta > lambda)
  }
  value
}

