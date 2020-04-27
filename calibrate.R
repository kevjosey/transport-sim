##########################################################################
### PURPOSE: Functions for Extending Inference using Entropy Balancing ###
### BY: Kevin Josey                                                    ###
##########################################################################

### Start auxillary functions ###

lagrange <- function(coefs, A, b) {
  
  temp <- sum(exp(-A %*% coefs))
  out <- temp + sum(b * coefs)
  return(out)
  
}

# Estimating equation for the SATE
esteq_sate <- function(X, Y, Z, weights, target, tau) {
  
  eq1 <- Z*weights*X - target
  eq2 <- (1 - Z)*weights*X - target
  eq3 <- weights*(Z*(Y - tau) - (1 - Z)*Y)
  
  eq <- c(eq1, eq2, eq3) 
  return(eq)
  
}

# Estimating equation for the PATE
esteq_pate <- function(S, X, Y, Z, weights, base_weights, target, tau) {
  
  eq1 <- S*(Z*weights*X - target)
  eq2 <- S*((1 - Z)*weights*X - target)
  eq3 <- (1 - S)*(base_weights*X - target)
  eq4 <- S*weights*(Z*(Y - tau) - (1 - Z)*Y)
  
  eq <- c(eq1, eq2, eq3, eq4) 
  return(eq)
  
}

esteq_sate_mom <- function(X, Y, Z, weights, target, tau) {
  
  eq1 <- weights*X - target
  eq2 <- weights*(Z*(Y - tau) - (1 - Z)*Y)
  
  eq <- c(eq1, eq2) 
  return(eq)
  
}

esteq_pate_mom <- function(S, X, Y, Z, weights, base_weights, target, tau) {
  
  eq1 <- S*(weights*X - target)
  eq2 <- (1 - S)*(base_weights*X - target)
  eq3 <- S*weights*(Z*(Y - tau) - (1 - Z)*Y)
  
  eq <- c(eq1, eq2, eq3) 
  return(eq)
  
}

### End auxillary functions ###

# mom == method of moments
calibrate <- function(X, Z, target, mom = FALSE, optim_ctrl = list(maxit = 500, reltol = 1e-20), ...) {
  
  if (!is.matrix(X))
    stop("X must be a matrix")
  
  if (!is.vector(target))
    stop("Z must be a vector")
  
  if (length(Z) != nrow(X))
    stop("length(Z) != nrow(X)")
  
  if (!is.vector(target))
    stop("target must be a vector")
  
  if (length(target) != ncol(X))
    stop("length(target) != ncol(X)")
  
  extraArgs <- list(...)
  
  if (length(extraArgs)) {
    
    arg <- names(formals(stats::optim))
    indx <- match(names(extraArgs), arg, nomatch = 0)
    if (any(indx == 0)) 
      stop(paste("Argument", names(extraArgs)[indx == 0], "not matched"))
    
  }
  
  fn <- match.fun(lagrange)
  n_1 <- length(Z)
  
  if (mom) {
    
    A <- cbind(X, (2*Z - 1))
    b <- c(n_1*target, 0)
    
  } else {
    
    A <- cbind(Z*X, (1 - Z)*X)
    b <- c(n_1*target, n_1*target)
    
  }
  
  # initialize coefs
  coefs_init <- rep(0, times = ncol(A))
  opt <- stats::optim(coefs_init, fn, method = "BFGS",
                      A = A, b = b, control = optim_ctrl)
  
  converged <- ifelse(opt$convergence == 0, TRUE, FALSE)
  coefs <- opt$par
  
  if (converged)
    weights <- c( exp(-A %*% coefs) )
  else
    weights <- NA
  
  out <- list(weights = weights,
              coefs = coefs,
              converged = converged,
              X = X, Z = Z, target = target,
              coefs_init = coefs_init,
              optim_ctrl = optim_ctrl,
              mom = mom)
  
  class(out) <- "calibrate"
  return(out)
  
}

# Estimation with target sample moments. Y is from the trial sample only.
estimate_sate <- function(obj, Y, ...) {
  
  if (!inherits(obj, "calibrate"))
    stop("obj must be of class \"calibrate\"")
  
  weights <- obj$weights
  coefs <- obj$coefs
  X <- obj$X
  Z <- obj$Z
  n <- length(Z)
  m <- ncol(X)
  target <- obj$target
  
  tau <- sum(weights*(2*Z - 1)*Y)/sum(Z*weights)
  conv <- TRUE
  
  if(obj$mom) {
    
    U <- matrix(0, ncol = m, nrow = m)
    v <- rep(0, times = m + 1)
    meat <- matrix(0, ncol = m + 1, nrow = m + 1)
      
    for (i in 1:n) {
      
      U[1:m,1:m] <- U[1:m,1:m] - weights[i] * (X[i,] %*% t(X[i,]))
      v[1:m] <- v[1:m] - (2*Z[i] - 1) * weights[i] * (Y[i] - Z[i]*tau) * X[i,]
      v[m + 1] <- v[m + 1] - weights[i]*Z[i]
      s <- esteq_sate_mom(X = X[i,], Y = Y[i], Z = Z[i], weights = weights[i], target = target, tau = tau)
      meat <- meat + s %*% t(s)
      
    }
    
    invbread <- matrix(0, nrow = m + 1, ncol = m + 1)
    invbread[1:m,1:m] <- U
    invbread[m + 1, ] <- v
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error")) {
      
      sandie_conv <- FALSE
      warning("Sandwich estimator is singular.")
      
    } else {
      
      sandwich <- bread %*% meat %*% t(bread)
      variance <- sandwich[m + 1, m + 1]
      
    }
    
  } else {
    
    U <- matrix(0, ncol = 2*m, nrow = 2*m)
    v <- rep(0, times = 2*m + 1)
    meat <- matrix(0, ncol = 2*m + 1, nrow = 2*m + 1)
    
    for (i in 1:n) {
      
      U[1:m,1:m] <- U[1:m,1:m] - Z[i] * weights[i] * (X[i,] %*% t(X[i,]))
      U[(m + 1):(2*m),(m + 1):(2*m)] <- U[(m + 1):(2*m),(m + 1):(2*m)] - (1 - Z[i]) * weights[i] * (X[i,] %*% t(X[i,]))
      v[1:m] <- v[1:m] - Z[i] * weights[i] * (Y[i] - tau) * X[i,]
      v[(m + 1):(2*m)] <- v[(m + 1):(2*m)] + (1 - Z[i]) * weights[i] * Y[i] * X[i,]
      v[2*m + 1] <- v[2*m + 1] - weights[i]*Z[i]
      s <- esteq_sate(X = X[i,], Y = Y[i], Z = Z[i], weights = weights[i], target = target, tau = tau)
      meat <- meat + s %*% t(s)
      
    }
    
    invbread <- matrix(0, nrow = 2*m + 1, ncol = 2*m + 1)
    invbread[1:(2*m),1:(2*m)] <- U
    invbread[2*m + 1, ] <- v
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error")) {
      
      sandie_conv <- FALSE
      warning("Sandwich estimator is singular.")
      
    } else {
      
      sandwich <- bread %*% meat %*% t(bread)
      variance <- sandwich[2*m + 1, 2*m + 1]
      
    }
    
  }
  
  out <- list(estimate = tau, variance = variance)
  return(out)
  
}

# Estimation with individual-level data. Y and Z are arbitrary for the target sample, 
# but must be included. Set them to 0.
estimate_pate <- function(obj, S, X, Y, Z, base_weights = NULL, ...) {
  
  if (!inherits(obj, "calibrate"))
    stop("obj must be of class \"calibrate\"")
  
  weights <- rep(1, times = length(S))
  weights[S == 1] <- obj$weights
  coefs <- obj$coefs
  n_0 <- sum(1 - S)
  n_1 <- sum(S)
  n <- n_1 + n_0
  m <- ncol(X)
  target <- obj$target
  
  if (is.null(base_weights))
    base_weights <- rep(1, times = length(S))
  
  if (length(base_weights) != length(S))
    stop("base_weights must have the same length as S")
  
  tau <- sum(S*(weights*(2*Z - 1)*Y)/sum(S*Z*weights))
  
  
  if (obj$mom) {
    
    U <- matrix(0, ncol = 2*m, nrow = 2*m)
    v <- rep(0, times = 3*m + 1)
    meat <- matrix(0, ncol = 2*m + 1, nrow = 2*m + 1)
    
    for (i in 1:n) {
      
      U[1:m,1:m] <- U[1:m,1:m] - S[i] * weights[i] * X[i,] %*% t(X[i,])
      U[1:m, (m + 1):(2*m)] <- U[1:m, (2*m + 1):(3*m)] - diag(S[i], m, m)
      U[(m + 1):(2*m),(m + 1):(2*m)] <- U[(2*m + 1):(3*m),(2*m + 1):(3*m)] - diag((1 - S[i]), m, m)
      
      v[1:m] <- v[1:m] - S[i]* weights[i] * (Y[i] - Z[i]*tau) * X[i,]
      v[3*m + 1] <- v[3*m + 1] - S[i]*weights[i]*Z[i]
      s <- esteq_pate_mom(X = X[i,], Y = Y[i], Z = Z[i], S = S[i], weights = weights[i], 
                             base_weights = base_weights[i], target = target, tau = tau)
      meat <- meat + s %*% t(s)
      
    }
    
    invbread <- matrix(0, nrow = 2*m + 1, ncol = 2*m + 1)
    invbread[1:(2*m),1:(2*m)] <- U
    invbread[2*m + 1, ] <- v
    invbread[(m + 1):(2*m),(m + 1):(2*m)] <- (n_1/n_0)*invbread[(m + 1):(2*m),(m + 1):(2*m)]
    meat[(m + 1):(2*m),(m + 1):(2*m)] <- (n_1/n_0)*meat[(m + 1):(2*m),(m + 1):(2*m)]
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error")) {
      
      sandwich <- NA
      variance <- NA
      
    } else {
      
      sandwich <- bread %*% meat %*% t(bread)
      variance <- sandwich[2*m + 1, 2*m + 1]
      
    }
    
  } else {
  
    U <- matrix(0, ncol = 3*m, nrow = 3*m)
    v <- rep(0, times = 3*m + 1)
    meat <- matrix(0, ncol = 3*m + 1, nrow = 3*m + 1)
    
    for (i in 1:n) {
      
      U[1:m,1:m] <- U[1:m,1:m] - S[i]*Z[i] * weights[i] * X[i,] %*% t(X[i,])
      U[(m + 1):(2*m),(m + 1):(2*m)] <- U[(m + 1):(2*m),(m + 1):(2*m)] - S[i]*(1 - Z[i]) * weights[i] * X[i,] %*% t(X[i,])
      U[1:m, (2*m + 1):(3*m)] <- U[1:m, (2*m + 1):(3*m)] - diag(S[i], m, m)
      U[(m + 1):(2*m),(2*m + 1):(3*m)] <- U[(m + 1):(2*m),(2*m + 1):(3*m)] - diag(S[i], m, m)
      U[(2*m + 1):(3*m),(2*m + 1):(3*m)] <- U[(2*m + 1):(3*m),(2*m + 1):(3*m)] - diag((1 - S[i]), m, m)
      
      v[1:m] <- v[1:m] - S[i]*Z[i] * weights[i] * (Y[i] - tau) * X[i,]
      v[(m + 1):(2*m)] <- v[(m + 1):(2*m)] + S[i]*(1 - Z[i]) * weights[i] * Y[i] * X[i,]
      v[3*m + 1] <- v[3*m + 1] - S[i]*weights[i]*Z[i]
      s <- esteq_pate(X = X[i,], Y = Y[i], Z = Z[i], S = S[i], weights = weights[i], 
                      base_weights = base_weights[i], target = target, tau = tau)
      meat <- meat + s %*% t(s)
      
    }
    
    invbread <- matrix(0, nrow = 3*m + 1, ncol = 3*m + 1)
    invbread[1:(3*m),1:(3*m)] <- U
    invbread[3*m + 1, ] <- v
    invbread[(2*m + 1):(3*m),(2*m + 1):(3*m)] <- (n_1/n_0)*invbread[(2*m + 1):(3*m),(2*m + 1):(3*m)]
    meat[(2*m + 1):(3*m),(2*m + 1):(3*m)] <- (n_1/n_0)*meat[(2*m + 1):(3*m),(2*m + 1):(3*m)]
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error")) {
      
      sandwich <- NA
      variance <- NA
      
    } else {
      
      sandwich <- bread %*% meat %*% t(bread)
      variance <- sandwich[3*m + 1, 3*m + 1]
      
    }
    
  }
  
  out <- list(estimate = tau, variance = variance)
  return(out)
  
}
