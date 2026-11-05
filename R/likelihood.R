update_tau <- function(
    Y,             # n x T (row i is y_i^T)
    Xi,            # T x 2 (time basis: intercept + time)
    Zi,            # T x p_z
    params,        # list: K, pi, beta[[k]], Psi[[k]], sigma2[k]
    pi_ik = NULL,  # n x K
    C = NULL,      # [ADD] n x q growth covariates (time-invariant)
    cov_pred = c("none","intercept","slope","both"),  # [ADD]
    eps = 1e-12,
    jitter = 1e-8,
    include_const = TRUE,
    protection = FALSE
){
  cov_pred <- match.arg(cov_pred)
  
  n    <- nrow(Y)
  Tobs <- ncol(Y)
  K    <- params$K
  
  # --- gate part: log pi (unchanged) ---
  if (!is.null(pi_ik)) {
    pi_ik <- as.matrix(pi_ik)
    if (protection) {
      pi_ik <- pmax(pi_ik, eps)
      pi_ik <- pi_ik / rowSums(pi_ik)
    }
    logpi <- log(pi_ik)    # n x K
  } else {
    if (protection) pi <- pmax(as.numeric(params$pi), eps)
    else            pi <- as.numeric(params$pi)
    pi <- pi / sum(pi)
    logpi <- matrix(rep(log(pi), each = n), nrow = n, ncol = K)
  }
  
  # --- small input checks (fail fast, avoids silent mis-spec) ---
  if (cov_pred != "none") {
    if (is.null(C)) stop("cov_pred != 'none' but C is NULL. Pass C = data[, growth_covariates].")
    C <- as.matrix(C)
    q <- ncol(C)
  } else {
    q <- 0L
  }
  
  tvec <- as.numeric(Xi[, 2])  # time values (length T), works for your Xi=cbind(1,0:(T-1))
  
  logP <- matrix(NA_real_, n, K)
  
  for (k in 1:K) {
    beta_k   <- as.numeric(params$beta[[k]])
    Psi_k    <- as.matrix(params$Psi[[k]])
    sigma2_k <- as.numeric(params$sigma2[k])
    
    # --- [MOD] build mu_ik (n x T) ---
    # baseline: Xa_k
    a_k <- beta_k[1:2]
    mu_time <- as.numeric(Xi %*% a_k)               # length T
    mu_ik   <- matrix(rep(mu_time, each = n), nrow = n, ncol = Tobs)  # n x T
    
    if (cov_pred != "none") {
      # expected beta length
      p_need <- 2L + if (cov_pred == "both") 2L*q else q
      if (length(beta_k) != p_need) {
        stop(sprintf("Class %d: length(beta_k)=%d but expected %d given cov_pred='%s' and q=%d.",
                     k, length(beta_k), p_need, cov_pred, q))
      }
      
      # intercept covariate effects: add constant shift across time
      if (cov_pred %in% c("intercept","both")) {
        gammaI <- beta_k[3:(2+q)]
        cI <- as.numeric(C %*% gammaI)              # n
        mu_ik <- mu_ik + cI                         # recycled across columns
      }
      
      # slope covariate effects: add time-varying shift proportional to t
      if (cov_pred %in% c("slope","both")) {
        if (cov_pred == "both") {
          gammaS <- beta_k[(3+q):(2+2*q)]
        } else { # slope only
          gammaS <- beta_k[3:(2+q)]
        }
        cS <- as.numeric(C %*% gammaS)              # n
        mu_ik <- mu_ik + outer(cS, tvec, "*")       # n x T
      }
    }
    
    # --- Omega_k (unchanged) ---
    Omega_k <- Zi %*% Psi_k %*% t(Zi) + sigma2_k * diag(Tobs)
    Omega_k <- (Omega_k + t(Omega_k)) / 2
    
    Rk <- chol(Omega_k)
    logdet <- 2 * sum(log(diag(Rk)))
    
    # residuals now use mu_ik (n x T)
    E_k <- Y - mu_ik                                # n x T
    Yt  <- forwardsolve(t(Rk), t(E_k))              # T x n
    qf  <- colSums(Yt^2)
    
    logP[, k] <- logpi[, k] - 0.5 * (logdet + qf)
  }
  
  m     <- apply(logP, 1, max)
  Zexp  <- exp(logP - m)
  denom <- rowSums(Zexp)
  
  Tau <- Zexp / denom
  
  if (protection) loglik <- sum(m + log(pmax(denom, .Machine$double.eps)))
  else            loglik <- sum(m + log(denom))
  if (include_const) loglik <- loglik - 0.5 * n * Tobs * log(2*base::pi)
  
  Nk <- colSums(Tau)
  
  list(
    Tau    = Tau,
    Nk     = Nk,
    loglik = as.numeric(loglik),
    logP   = logP
  )
}


# Observed log-likelihood for GMM
# If logP is given: use it directly (fast path, no extra checks)
# Else: recompute from (Y, Xi, Zi, params) (slow path)
loglik_GMM <- function(
    Y = NULL,          # n x T matrix (rows are subjects); needed if recomputing
    Xi = NULL,         # T x 2 (time basis: intercept + time)
    Zi = NULL,         # T x p_z
    params = NULL,     # list: K, pi, beta[[k]], Psi[[k]], sigma2[k]
    pi_ik = NULL,      # n x K individual mixing proportion
    C = NULL,          # [ADD] n x q growth covariates (time-invariant)
    cov_pred = c("none","intercept","slope","both"),  # [ADD]
    logP = NULL,       # n x K matrix: log pi_k + log N(y_i | mu_ik, Omega_k), no constant
    include_const = TRUE,
    eps = 1e-12,
    protection = FALSE
){
  cov_pred <- match.arg(cov_pred)
  
  # -------- Case A: fast path using logP --------
  if (!is.null(logP)) {
    m  <- apply(logP, 1, max)
    ll <- sum(m + log(rowSums(exp(logP - m))))
    if (include_const && !is.null(Y)) {
      ll <- ll - 0.5 * nrow(Y) * ncol(Y) * log(2*base::pi)
    }
    return(as.numeric(ll))
  }
  
  # -------- Case B: recompute logP from (Y, Xi, Zi, params) --------
  n    <- nrow(Y); Tobs <- ncol(Y); K <- params$K
  
  # Input Pi_ik or Pi_k ------------------------------------
  if (!is.null(pi_ik)) {
    pi_ik <- as.matrix(pi_ik)
    if (protection) {
      pi_ik <- pmax(pi_ik, eps)
      pi_ik <- pi_ik / rowSums(pi_ik)
    }
    logpi <- log(pi_ik)    # n x K
  } else {
    if (protection) pi <- pmax(as.numeric(params$pi), eps)
    else            pi <- as.numeric(params$pi)
    pi <- pi / sum(pi)
    logpi <- matrix(rep(log(pi), each = n), nrow = n, ncol = K)
  }
  
  # [ADD] checks & helpers
  if (cov_pred != "none") {
    if (is.null(C)) stop("cov_pred != 'none' but C is NULL. Pass C = data[, growth_covariates].")
    C <- as.matrix(C)
    q <- ncol(C)
    if (ncol(Xi) < 2) stop("Xi must have at least 2 columns (intercept + time) when cov_pred != 'none'.")
  } else {
    q <- 0L
  }
  tvec <- if (ncol(Xi) >= 2) as.numeric(Xi[, 2]) else 0:(Tobs-1)  # safe fallback
  
  logP <- matrix(NA_real_, n, K)
  
  for (k in 1:K) {
    beta_k   <- as.numeric(params$beta[[k]])
    Psi_k    <- as.matrix(params$Psi[[k]])
    sigma2_k <- as.numeric(params$sigma2[k])
    
    # ========= [MOD] build mu_ik (n x T) =========
    a_k <- beta_k[1:2]
    mu_time <- as.numeric(Xi %*% a_k)                               # length T
    mu_ik   <- matrix(rep(mu_time, each = n), nrow = n, ncol = Tobs) # n x T
    
    if (cov_pred != "none") {
      p_need <- 2L + if (cov_pred == "both") 2L*q else q
      if (length(beta_k) != p_need) {
        stop(sprintf("Class %d: length(beta_k)=%d but expected %d given cov_pred='%s' and q=%d.",
                     k, length(beta_k), p_need, cov_pred, q))
      }
      
      if (cov_pred %in% c("intercept","both")) {
        gammaI <- beta_k[3:(2+q)]
        cI <- as.numeric(C %*% gammaI)          # n
        mu_ik <- mu_ik + cI                      # add to all time points
      }
      
      if (cov_pred %in% c("slope","both")) {
        if (cov_pred == "both") {
          gammaS <- beta_k[(3+q):(2+2*q)]
        } else {
          gammaS <- beta_k[3:(2+q)]
        }
        cS <- as.numeric(C %*% gammaS)          # n
        mu_ik <- mu_ik + outer(cS, tvec, "*")   # n x T
      }
    }
    
    # ========= Omega_k (same as update_tau for stability) =========
    Omega_k <- Zi %*% Psi_k %*% t(Zi) + sigma2_k * diag(Tobs)
    Omega_k <- (Omega_k + t(Omega_k)) / 2   # [MOD] numerical symmetry
    
    # Cholesky
    Rk <- chol(Omega_k)
    logdet <- 2 * sum(log(diag(Rk)))
    
    # residuals using mu_ik
    E_k <- Y - mu_ik                          # n x T
    Yt  <- forwardsolve(t(Rk), t(E_k))        # T x n
    qf  <- colSums(Yt^2)
    
    logP[, k] <- logpi[, k] - 0.5 * (logdet + qf)
  }
  
  m  <- apply(logP, 1, max)
  ll <- sum(m + log(rowSums(exp(logP - m))))
  if (include_const) ll <- ll - 0.5 * n * Tobs * log(2*base::pi)
  return(as.numeric(ll))
}


#  Weighted within-class log-likelihood 
update_loglike <- function(
    omega, # individual omega
    y,
    x,
    beta,
    n,
    include_const = TRUE,  # include sum(w)*T*log(2π) if TRUE
    weights = rep(1,n),
    Debug = FALSE
) {
  if(Debug)browser()
  ll <- 0
  Tobs <- nrow(omega)
  
  log_det <- as.numeric(determinant( omega, logarithm = T)$modulus)
  res <- y- (x %*% beta)
  #quad_form <- t(res) %*% kronecker(diag(weights, n, n),  pracma::pinv(omega)) %*% res # too expensive?
  E <- matrix(res, nrow = Tobs, ncol = n, byrow = FALSE)
  P <- pracma::pinv(omega)
  quad_form <- sum(weights * colSums(E * (P %*% E)))
  
  const_term <- if (include_const) sum(weights) * Tobs * log(2*base::pi) else 0
  
  ll <- -0.5 * (sum(weights)*log_det + quad_form + const_term)
  
  return(list(
    Qk_no_pi     = as.numeric(ll),
    sum_w_logdet = as.numeric(sum(weights) * log_det),
    quad_wsum    = as.numeric(quad_form ),
    const_term   = as.numeric(const_term),
    Nk           = as.numeric(sum(weights))
  ))
  
  

}
