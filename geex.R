
OM_EE <- function(data){
  Z <- data$Z
  S <- data$S
  Y <- data$Y
  X <- as.matrix(data[,2:(ncol(data) - 2)])
  
  function(theta){ 
    #outcome model 
    beta <- theta[1:ncol(X)]
    alpha <- theta[(ncol(X)+1):(2*ncol(X))]
    mu1 <- theta[2*ncol(X) + 1]
    mu0 <- theta[2*ncol(X) + 2]
    muate <- theta[2*ncol(X) + 3]
    m_Z1 <- X %*% beta
    m_Z0 <- X %*% alpha
    ols_Z1 <- crossprod(X, (S*Z)*(Y - m_Z1))
    ols_Z0 <- crossprod(X, (S*(1-Z))*(Y - m_Z0))
    #estimates
    mean1 <- (1-S)*(m_Z1-mu1) 
    mean0 <- (1-S)*(m_Z0-mu0) 
    ate <- (1-S)*(m_Z1-m_Z0-muate) 
    c(ols_Z1,ols_Z0,mean1, mean0, ate)
  }
}

IOW_EE <- function(data) { 
  
  Z <-data$Z
  S <- data$S
  Y <- data$Y  
  X <- as.matrix(data[,2:(ncol(data) - 2)])
  
  function(theta){
    #participation model
    lp <- X %*% theta[1:ncol(X)]
    ps <- plogis(lp)
    score_eqns <- crossprod(X, S-ps)
    #treatment model
    lp2 <- theta[ncol(X) + 1] 
    pa <- plogis(lp2)
    score_eqns2<-crossprod(1, S*(Z - pa) )
    w = (Z * S*(1-ps))/(ps*pa) + ((1 - Z)*S*(1-ps))/(ps*(1-pa)) 
    #outcome model
    m_Z1 <- theta[ncol(X) + 2] 
    m_Z0 <- theta[ncol(X) + 3] 
    linear_eqns1 <- crossprod(1, (S*Z*w)*(Y -  m_Z1) )
    linear_eqns0 <- crossprod(1, (S*(1-Z)*w)*(Y - m_Z0) )
    mu1 <- theta[ncol(X) + 4]
    mu0 <- theta[ncol(X) + 5]
    muate <- theta[ncol(X) + 6]
    #estimates
    mean1 <- (1-S)*(m_Z1 - mu1) 
    mean0 <- (1-S)*(m_Z0 - mu0)
    ate <- (1-S)*(m_Z1 - m_Z0 - muate)
    c(score_eqns, score_eqns2, linear_eqns1, linear_eqns0, mean1, mean0, ate)
  }
  
}

DR_EE <- function(data) {
  
  Z <- data$Z
  S <- data$S
  Y <- data$Y
  X <- as.matrix(data[,2:(ncol(data) - 2)])
  
  function(theta){
    #participation model
    lp  <- X %*% theta[1:ncol(X)]
    ps <- plogis(lp)
    score_eqns <- crossprod(X, S-ps)
    #treatment model
    lp2  <- theta[ncol(X) + 1] 
    pa <- plogis(lp2)
    score_eqns2 <- crossprod(rep(1, times = nrow(X)), S*(Z - pa) )
    w <- (Z * S*(1-ps))/(ps*pa) + ((1 - Z)*S*(1-ps))/(ps*(1-pa)) 
    #outcome model
    beta <- theta[(ncol(X) + 2):(2*ncol(X) + 1)]
    alpha <- theta[(2*ncol(X) + 2):(3*ncol(X) + 1)]
    mu1 <- theta[3*ncol(X) + 2]
    mu0 <- theta[3*ncol(X) + 3]
    mu <- theta[3*ncol(X) + 4]
    m_Z1 <- X %*% beta
    m_Z0 <- X %*% alpha
    ols_Z1 <- crossprod(X, (S*Z)*(Y - m_Z1))
    ols_Z0 <- crossprod(X, (S*(1-Z))*(Y - m_Z0))
    ey1 <- w*S*Z*(Y-m_Z1) + (1-S)*m_Z1
    ey0 <- w*S*(1-Z)*(Y-m_Z0) + (1-S)*m_Z0
    #estimates
    mean1 <- ey1-(1-S)*mu1
    mean0 <- ey0-(1-S)*mu0
    ate <- ey1-ey0- (1-S)*mu
    c(score_eqns, score_eqns2, ols_Z1, ols_Z0, mean1, mean0, ate)   
  }
}
