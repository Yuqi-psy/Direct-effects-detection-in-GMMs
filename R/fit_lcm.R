# ---------- Main: weighted LCM estimation ----------
fit_lcm <- function(
    y,        # character vector of repeated-measure column names (length T)
    ID,       # character: ID column name
    data,     # data.frame
    growth_covariates = NULL,           # covs names
    cov_pred = c("none","intercept","slope","both"),
    weights  = NULL,        # if NULL, then weights are rep(1,n)
    active   = NULL,        # estimate fixed effects set 1, not estimate set 0.
    lcm_design = NULL,
    
    max_iter = 1000,
    tol = 1e-6, active_tol = 1e-4,
    
    penalty = FALSE,
    lambda = c(0,0),
    
    
    verbose = FALSE,
    Debug = FALSE,
    freeze_beta   = FALSE, 
    freeze_theta  = FALSE,   
    freeze_sigma2 = FALSE,   
    beta_init     = NULL,   
    theta_init    = NULL,
    sigma2_init   = NULL
    
){
  if (Debug) browser()
  gmm_setup() 
  cov_pred <- match.arg(cov_pred)
  penalty_k <- NULL
  active_out <- NULL
  
  # ---- basic prep ----
  if (is.null(lcm_design)) {
    lcm_design <- prep_lcm_design(
      y = y,
      data = data,
      growth_covariates = growth_covariates,
      cov_pred = cov_pred,
      penalty = penalty
    )
  }
  stopifnot(lcm_design$n == nrow(data), lcm_design$Tobs == length(y))
  
  n          <- lcm_design$n
  Tobs       <- lcm_design$Tobs
  y_vec      <- lcm_design$y_vec
  t_long     <- lcm_design$t_long
  x          <- lcm_design$x
  z_i        <- lcm_design$z_i
  X_0        <- lcm_design$X_0
  X_1        <- lcm_design$X_1

  
  if (is.null(weights)) weights <- rep(1, n)
  w <- pmax(as.numeric(weights), 1e-12)

  
  # ---- initialization ----

  if (is.null(beta_init)) {
    beta   <- solve(t(x) %*% x, t(x) %*% y_vec)     # OLS
  } else {
    beta   <- as.matrix(beta_init)
  }
  # checking: when penalization applied, their must be covariates to predict the growth factors
  if (penalty) {
    if (is.null(beta_init)) {
      stop("When penalty = TRUE, beta_init must be provided.")
    }
    if (is.null(growth_covariates) || length(growth_covariates) == 0L || cov_pred == "none") {
      stop("Penalty is currently only supported when growth_covariates are provided and cov_pred is not 'none'.")
    }
  
  }
  
  res0 <- y_vec - x %*% beta
  if (is.null(sigma2_init)) {
    sigma2 <- as.numeric(var(res0))
  } else {
    sigma2 <- as.numeric(sigma2_init)
  }
  
  if (is.null(theta_init)) {
    theta  <- diag(2) * 0.1                         
  } else {
    theta  <- as.matrix(theta_init)
  }
  
  # histories
  cd_loglik_history <- numeric(max_iter)
  beta_history     <- vector("list", max_iter)
  theta_history    <- vector("list", max_iter)
  sigma2_history   <- numeric(max_iter)
  
  # ---- MM/ECM-style loop ----
  for (iter in 1:max_iter) {
    
    # 1) Build Ω from CURRENT (theta, sigma2)
    omega <- update_omega(z = z_i, theta = theta, sigma2 = sigma2, n = n)
    
    # 2) Update beta (weighted GLS)  
    if (!freeze_beta) {
      if (penalty){


        beta_out <- update_beta_penalty(
          X_0 = X_0, X_1 = X_1,
          beta_ini = beta,
          omega = omega, y = y_vec, n = n, weights = w,
          lambda = lambda, trace = F, Debug = Debug)
        
        beta <- beta_out$beta_std
        penalty_k <- beta_out$penalty
        
      } else if (!is.null(active)) {
        
       beta_out <- update_beta_active(
         X_0 = X_0,
         X_1 = X_1, 
         omega = omega,        
         y = y_vec,            # nT x 1
         n = n,
         active = active,       # when refit, for estimated parameters (1), not estimate parameters (0)
         weights = w,
         Debug = Debug
       )
       beta <- beta_out$beta_std
       if(!is.null(active)) active_out <- beta_out$active
      }
        else {
        beta <- update_beta(
          x       = x,
          omega   = omega,
          y       = y_vec,
          n       = n,
          weights = w,
          Debug   = Debug
        ) 
      }
    } else if (!is.null(beta_init)) {
      beta <- as.matrix(beta_init)                          
    } 
    
    # 3) Update sigma^2 (weighted MM) 
    if (!freeze_sigma2) {                                  
      sigma2 <- update_sigma2(
        y       = y_vec,
        x       = x,
        beta    = beta,
        omega   = omega, 
        sigma2  = sigma2,
        n       = n,
        weights = w,
        Debug   = Debug
      )
    } else if (!is.null(sigma2_init)) {
      sigma2 <- as.numeric(sigma2_init)                  
    }
    
    # 4) Update theta (weighted MM) 
    if (!freeze_theta) {                                   
      th_out <- update_theta2(
        sigma2  = sigma2,
        omega   = omega,  
        y       = y_vec,
        x       = x,
        z       = z_i,
        beta    = beta,
        theta   = theta,
        n       = n,
        weights = w,
        Debug   = Debug
      )
      theta <- th_out$theta
    } else if (!is.null(theta_init)) {
      theta <- as.matrix(theta_init)                   
    }
    
   
    omega_q <- update_omega(z = z_i, theta = theta, sigma2 = sigma2, n = n) 
    
    beta_history[[iter]]    <- beta
    theta_history[[iter]]   <- theta
    sigma2_history[iter]    <- sigma2
    


    qparts <- update_loglike(
      omega         = omega_q,      
      y             = y_vec,
      x             = x,
      beta          = beta,
      n             = n,
      weights       = w,
      include_const = TRUE,
      Debug         = Debug
    )
    
    # Public penalty is computed on standardized DE coefficients.
    # if (penalty) {
    #   bd <- drop(beta)[-(1:2)]
    #   penalty_k <- list(pen1 = lambda[1] * sum(abs(bd)),
    #                     pen2 = lambda[2] * sum(bd^2))
    # }
    # 
    # penalized CD log likelihood values
    Qk_pen <- qparts$Qk_no_pi
    if (penalty && !is.null(penalty_k)) {
      Qk_pen <- Qk_pen - (penalty_k$pen1 + penalty_k$pen2)
    }
    
    cd_loglik_history[iter] <- if (penalty) Qk_pen else qparts$Qk_no_pi
    
    # 6) convergence check on log-likelihood (monotone ascent expected)
    if (iter > 1) {
      diff_ll <- cd_loglik_history[iter] - cd_loglik_history[iter - 1]
      if (verbose)
        cat(sprintf("[LCM iter %d] ΔQk = %.3e\n", iter, diff_ll))
      if (diff_ll < tol) {
        if (verbose) cat(sprintf("LCM Converged at iter %d: Δll = %.3e\n", iter, diff_ll))
        break
      }
    }
  }
  
  # extract active set index
  if (penalty) {
    beta_pen <- drop(beta)[-(1:2)] # STANDARDIZED coefficients
    active_out <- as.integer(abs(beta_pen) > active_tol)
    
  } else if (!is.null(active)) {
    active_out <- as.integer(active != 0)
    
  } else {
    active_out <- NULL
  }
  
  
  iter <- min(iter, max_iter)
  out <- list(beta = beta, theta = theta, sigma2 = sigma2, 
              iter = iter, penalty_k = penalty_k, active = active_out, 
              cd_loglik = cd_loglik_history[seq_len(iter)], beta_values = beta_history[seq_len(iter)], 
              theta_values = theta_history[seq_len(iter)], sigma_values = sigma2_history[seq_len(iter)], 
              Qk_no_pi = qparts$Qk_no_pi, Qk_pen = Qk_pen, sum_w_logdet = qparts$sum_w_logdet, 
              quad_wsum = qparts$quad_wsum, const_term = qparts$const_term, 
              Nk = qparts$Nk, parameter_scale = "standardized", active_scale = "standardized")
  return(out)
}
