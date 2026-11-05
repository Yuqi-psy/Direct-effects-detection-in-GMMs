# compute the pi_ik when including covariate that predicts the class membership 
update_pi_ik <- function(
  X, # n*p_x  matrix for covariates without intercept
  params # parameters set including logit intercept alpha and slope beta 
) {
  # extract parameters (first class as reference class)
  alpha <- as.numeric(params$alpha) # K-1 vector
  Beta <- as.matrix(params$Beta) # (k-1)*p_x matrix 
  
  # 1.compute Eta: Eta[i, k-1] = alpha_k + x_i^T beta_k,  k = 2..K
  eta <- sweep((as.matrix(X) %*% t(Beta)), 2, alpha, "+")
  Eta <- cbind(0, eta) # add the reference group
  
  # 2. compute Pi_ik: Pi_ik = exp(Eta)/sum^K_k=1 exp(Eta) 
  m     <- apply(Eta, 1, max)     # length n
  numer <- exp(Eta - m)
  denom <- rowSums(numer) 
  Pi_ik <- numer/denom # n*k. rowsums equal to 1
  
  Pi_ik
}

# this is used to compute the logit coefficients of covs in prediting class membership
update_gate <- function(
    X_class = NULL, # design matrix including covs
    pi_ik   = NULL, # mixing proportions of each i and each k
    params  = NULL, # current gate parameters 
    tau     = NULL,
    dcrit = 1e-6,
    maxiter = 991119,
    verbose = FALSE,
    
    Debug = FALSE
){
  if(Debug) browser()
  
  # 1. prepare
  X_class <- cbind(1,X_class) # n x q0
  n <- nrow(X_class)
  q0 <- ncol(X_class)
  K   <- ncol(pi_ik)
  K_nr  <- K - 1
  iXXX <- solve(t(X_class) %*% X_class) %*% t(X_class)
  
  
  alpha <- as.numeric(params$gate$alpha)       # K-1 vector
  Beta  <- as.matrix(params$gate$Beta)         # (k-1) x (q0-1) matrix 
  Gamma_nr <- cbind(alpha, Beta)               # (K-1) x q0
  tau_nr   <- tau[, -1, drop = FALSE] 
  tau_r    <- tau[, 1, drop = TRUE]
  pi_nr <- pi_ik[, -1, drop = FALSE] 
  
  eta_nr <- X_class %*% t(Gamma_nr)                  # n x (K-1)
  exp_eta <- exp(eta_nr)                 # n x (K-1)
  denom   <- 1 + rowSums(exp_eta)        # length n
  pi_r    <- 1 / denom                   # length n
  pi_nr   <- exp_eta / denom             # n x (K-1)
  dev_old <- -2 * (sum(tau_r * log(pi_r))
                   +sum(tau_nr * log(pi_nr)))
  
  iter <- 0
  dif  <- Inf
  while(iter < maxiter){
    iter = iter + 1
    
    Z_nr <- eta_nr + 4 * (tau_nr - pi_nr) ## ! might need to change the constrant
    
    # MM least squares update:
    Gamma_t_new  <- iXXX %*% Z_nr                    # q0 x (K-1)
    Gamma_nr_new <- t(Gamma_t_new)                   # (K-1) x q0
    
    # update eta
    eta_nr <- X_class %*% Gamma_t_new                  # n x (K-1)
    
    # recompute pi (inline baseline-softmax) 
  
    exp_eta <- exp(eta_nr)                 # n x (K-1)
    denom   <- 1 + rowSums(exp_eta)        # length n
    pi_r    <- 1 / denom                   # length n
    pi_nr   <- exp_eta / denom             # n x (K-1)
    
    # compute the deviance
    dev_new <- -2 * (sum(tau_r * log(pi_r))
                     +sum(tau_nr * log(pi_nr)))
    dif <- 2 * (dev_old - dev_new) / (dev_old + dev_new)
    
    if (verbose) cat(iter, dev_old, dev_new, dif, "\n")
    
    Gamma_nr <- Gamma_nr_new
    if (!is.finite(dif) || dif < dcrit) break
    
    dev_old <- dev_new
  }
  
  eta_nr <- X_class %*% t(Gamma_nr)
  exp_eta <- exp(eta_nr)                     
  denom   <- 1 + rowSums(exp_eta)            
  pi_r    <- 1 / denom                       
  pi_nr   <- exp_eta / denom
  dev_new <- -2 * (sum(tau_r * log(pi_r))
                   +sum(tau_nr * log(pi_nr)))

  
  params$gate$alpha <- Gamma_nr[,1]
  params$gate$Beta  <- Gamma_nr[,-1, drop = FALSE]
  
  list(
    params  = params,
    deviance = dev_new,
    iter = iter,
    pi_ik = cbind(pi_r, pi_nr)
  )
}
