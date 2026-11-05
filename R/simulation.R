sim_gmm_cov <- function(
    n,
    time,                  # scalar Tobs, or numeric time vector
    K,
    beta_list,             # length K, each length q
    theta_list,            # length K, each q x q
    sigma2_vec,            # length K
    
    X = NULL,              # n x p covariates; used for BOTH gate and growth factors
    
    # only used when X is NULL and you want the function to generate X
    x_type = NULL,         # length p: "cont" or "bin"
    x_mean = NULL,         # length p, only for cont
    x_sd   = NULL,         # length p, only for cont
    x_prob = NULL,         # length p, only for bin
    
    alpha = NULL,          # length K-1
    Beta  = NULL,          # (K-1) x p
    
    Gamma_list = NULL,     # length K list, each q x p
    pi = NULL,             # length K, used only when p = 0
    
    Lambda = NULL,
    Debug = FALSE,
    seed  = NULL
){
  if (!is.null(seed)) set.seed(seed)
  if (Debug) browser()
  
  n <- as.integer(n)
  
  # -------------------------
  # 1. time / Lambda
  # -------------------------
  if (length(time) == 1L) {
    Tobs  <- as.integer(time)
    t_vec <- 0:(Tobs - 1)
  } else {
    t_vec <- as.numeric(time)
    Tobs  <- length(t_vec)
  }
  
  q <- length(beta_list[[1]])
  
  if (is.null(Lambda)) {
    if (q != 2L) {
      stop("If Lambda is NULL, q must be 2. Otherwise please provide Lambda.")
    }
    Lambda <- cbind(1, t_vec)   # intercept + linear slope
  } else {
    Lambda <- as.matrix(Lambda)
    stopifnot(nrow(Lambda) == Tobs, ncol(Lambda) == q)
  }
  
  # -------------------------
  # 2. generate X if needed
  # -------------------------
  if (is.null(X)) {
    stopifnot(!is.null(x_type))
    
    p <- length(x_type)
    X <- matrix(NA_real_, nrow = n, ncol = p)
    
    if (is.null(x_mean)) x_mean <- rep(0, p)
    if (is.null(x_sd))   x_sd   <- rep(1, p)
    if (is.null(x_prob)) x_prob <- rep(0.5, p)
    
    for (j in seq_len(p)) {
      if (x_type[j] %in% c("cont", "continuous")) {
        X[, j] <- rnorm(n, mean = x_mean[j], sd = x_sd[j])
      } else if (x_type[j] %in% c("bin", "binary")) {
        X[, j] <- rbinom(n, size = 1, prob = x_prob[j])
      } else {
        stop("x_type must be 'cont' or 'bin'.")
      }
    }
    
    colnames(X) <- paste0("x", seq_len(p))
    
  } else {
    X <- as.matrix(X)
    stopifnot(nrow(X) == n)
    if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  }
  
  p <- ncol(X)
  
  # -------------------------
  # 3. basic checks
  # -------------------------
  stopifnot(length(beta_list)  == K,
            length(theta_list) == K,
            length(sigma2_vec) == K)
  
  for (k in seq_len(K)) {
    stopifnot(length(beta_list[[k]]) == q)
    stopifnot(all(dim(theta_list[[k]]) == c(q, q)))
    stopifnot(sigma2_vec[k] > 0)
  }
  
  # -------------------------
  # 4. Gamma_list
  # -------------------------
  if (is.null(Gamma_list)) {
    Gamma_list <- replicate(K, matrix(0, q, p), simplify = FALSE)
  } else {
    stopifnot(is.list(Gamma_list), length(Gamma_list) == K)
    for (k in seq_len(K)) {
      Gamma_list[[k]] <- as.matrix(Gamma_list[[k]])
      stopifnot(nrow(Gamma_list[[k]]) == q, ncol(Gamma_list[[k]]) == p)
    }
  }
  
  # -------------------------
  # 5. class probabilities
  # -------------------------
  if (p == 0) {
    stopifnot(!is.null(pi))
    pi <- as.numeric(pi)
    stopifnot(length(pi) == K, all(pi >= 0))
    pi <- pi / sum(pi)
    
    pi_ik <- matrix(rep(pi, each = n), nrow = n, ncol = K, byrow = FALSE)
    colnames(pi_ik) <- paste0("class", seq_len(K))
    
  } else {
    stopifnot(!is.null(alpha), !is.null(Beta))
    stopifnot(length(alpha) == K - 1)
    
    Beta <- as.matrix(Beta)
    stopifnot(all(dim(Beta) == c(K - 1, p)))
    
    eta_nr   <- sweep(X %*% t(Beta), 2, alpha, "+")   # n x (K-1)
    eta_full <- cbind(0, eta_nr)                      # baseline: class 1
    
    row_max <- apply(eta_full, 1, max)
    exp_eta <- exp(eta_full - row_max)
    pi_ik   <- exp_eta / rowSums(exp_eta)
    
    colnames(pi_ik) <- paste0("class", seq_len(K))
  }
  
  # -------------------------
  # 6. draw class
  # -------------------------
  class <- apply(pi_ik, 1, function(prob) sample.int(K, 1, prob = prob))
  
  # -------------------------
  # 7. simulate eta and Y
  # -------------------------
  eta <- matrix(NA_real_, nrow = n, ncol = q)
  Y   <- matrix(NA_real_, nrow = n, ncol = Tobs)
  
  for (k in seq_len(K)) {
    idx <- which(class == k)
    nk  <- length(idx)
    if (nk == 0L) next
    
    beta_k   <- as.numeric(beta_list[[k]])
    theta_k  <- as.matrix(theta_list[[k]])
    sigma2_k <- sigma2_vec[k]
    Gk       <- Gamma_list[[k]]   # q x p
    
    mu_eta <- matrix(beta_k, nrow = nk, ncol = q, byrow = TRUE)
    if (p > 0) {
      mu_eta <- mu_eta + X[idx, , drop = FALSE] %*% t(Gk)
    }
    
    u_k <- MASS::mvrnorm(n = nk, mu = rep(0, q), Sigma = theta_k)
    if (nk == 1L) u_k <- matrix(u_k, nrow = 1L)
    
    eta_k <- mu_eta + u_k
    eps_k <- matrix(rnorm(nk * Tobs, mean = 0, sd = sqrt(sigma2_k)),
                    nrow = nk, ncol = Tobs)
    
    Y[idx, ]   <- eta_k %*% t(Lambda) + eps_k
    eta[idx, ] <- eta_k
  }
  
  # -------------------------
  # 8. output
  # -------------------------
  eta_df <- as.data.frame(eta)
  colnames(eta_df) <- paste0("eta", seq_len(q))
  
  Y_df <- as.data.frame(Y)
  colnames(Y_df) <- paste0("y.", 0:(Tobs - 1))
  
  dat <- data.frame(
    id    = seq_len(n),
    class = factor(class, levels = seq_len(K)),
    X,
    eta_df,
    Y_df,
    check.names = FALSE
  )
  
  list(
    data       = dat,
    X          = X,
    pi_ik      = pi_ik,
    class      = class,
    alpha      = alpha,
    Beta       = Beta,
    Gamma_list = Gamma_list,
    beta_list  = beta_list,
    theta_list = theta_list,
    sigma2_vec = sigma2_vec,
    Lambda     = Lambda,
    time       = t_vec
  )
}
