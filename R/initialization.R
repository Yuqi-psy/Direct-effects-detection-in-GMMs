kmeans_init <- function(
    y_names, ID, data, K, 
    growth_covariates = NULL,           
    cov_pred = c("none","intercept","slope","both"),
    pool_fit = NULL,
    use_residuals = TRUE,
    refine_classes = TRUE,
    refine_max_iter = 10,
    refine_tol = 1e-3,
    lcm_design = NULL,
    min_class_n = NULL,
    nstart = 10,
    jitter = 0.02,
    seed = NULL,
    Debug = TRUE
){
  cov_pred <- match.arg(cov_pred)
  
  if (!is.null(seed)) set.seed(seed)
  if (Debug) browser()
  if (is.null(lcm_design)) 
  lcm_design <- prep_lcm_design(y_names, data, growth_covariates, cov_pred)
  data <- lcm_design$data_std
  
  # -------------------------
  # Basic shapes
  # -------------------------
  n    <- nrow(data)
  if (is.null(min_class_n)) {
    min_class_n <- max(5L, ceiling(0.01 * n))
  }
  Tobs <- length(y_names)
  Y    <- as.matrix(data[, y_names, drop = FALSE])   # n x T
  Zi   <- cbind(1, 0:(Tobs - 1))                     # T x 2
  
  has_direct <- !is.null(growth_covariates) && length(growth_covariates) > 0 && cov_pred != "none"
  
  # -------------------------
  # Pooled LCM if not provided
  # -------------------------
  if (is.null(pool_fit)) {
    pool_fit <- fit_lcm(
      y = y_names, ID = ID, data = data, growth_covariates = growth_covariates,
      cov_pred = cov_pred, weights = rep(1, n), lcm_design = lcm_design,
      max_iter = 2000, tol = 1e-6, penalty = FALSE, Debug = Debug)
  }
  if (!is.null(pool_fit$beta)) {
    beta_full_pool <- as.numeric(pool_fit$beta)
    Psi_pool <- as.matrix(pool_fit$theta)
    sigma2_pool <- as.numeric(pool_fit$sigma2)
  } else {
    beta_full_pool <- as.numeric(pool_fit$`Fixed effects (Beta)`)
    Psi_pool <- as.matrix(pool_fit$`Random effects (Theta)`)
    sigma2_pool <- as.numeric(pool_fit$`Residual variance (sigma2)`)
  }
  # A supplied pool_fit must already use this STD block design.
  if (length(beta_full_pool) != ncol(lcm_design$x)) stop("pool_fit beta order/dimensions must match the block LCM design.")
  # -------------------------
  # Extract pooled parameters
  # -------------------------

  beta0_pool     <- beta_full_pool[1:2]      # baseline intercept + slope
  gamma_pool     <- beta_full_pool[-(1:2)]   # direct effects part (possibly numeric(0))

  
  # -------------------------
  # Whitening matrix
  # -------------------------
  Omega_pool <- Zi %*% Psi_pool %*% t(Zi) + sigma2_pool * diag(Tobs)
  R          <- chol(Omega_pool)
  
  # whiten Y
  Yw <- t(backsolve(R, t(Y), transpose = FALSE))   # n x T
  
  # whiten baseline design only
  X0w <- backsolve(R, Zi, transpose = FALSE)       # T x 2
  
  # -------------------------
  # Build pooled whitened mean
  # -------------------------
  if (!has_direct) {
    # no direct effects: common pooled mean for everyone
    mu_pool_w <- as.numeric(X0w %*% beta0_pool)    # length T
    
    Fmat <- if (use_residuals) {
      sweep(Yw, 2, mu_pool_w, "-")
    } else {
      Yw
    }
    
  } else {
    # direct effects present: pooled mean is subject-specific
    X_cov <- as.matrix(data[, growth_covariates, drop = FALSE])   # n x p
    p_cov <- ncol(X_cov)
    
    mu_pool_w <- matrix(NA_real_, nrow = n, ncol = Tobs)
    
    for (i in 1:n) {
      Xi_i <- Zi   # start with baseline block
        # append x_j * intercept-loading
        # append x_j * slope-loading
        # append x_j * Zi = [x_j*intercept, x_j*slope]
      if (cov_pred %in% c("intercept", "both"))
        Xi_i <- cbind(Xi_i, matrix(rep(X_cov[i, ], each = Tobs), Tobs))
      if (cov_pred %in% c("slope", "both"))
        Xi_i <- cbind(Xi_i, matrix(rep(X_cov[i, ], each = Tobs), Tobs) * Zi[, 2])
      
      Xiw_i <- backsolve(R, Xi_i, transpose = FALSE)
      mu_pool_w[i, ] <- as.numeric(Xiw_i %*% beta_full_pool)
    }
    
    Fmat <- if (use_residuals) {
      Yw - mu_pool_w
    } else {
      Yw
    }
  }
  
  # -------------------------
  # K-means in whitened space
  # -------------------------
  km <- stats::kmeans(Fmat, centers = K, nstart = nstart, iter.max = 100)
  clusters <- km$cluster
  centers  <- km$centers
  
  # -------------------------
  # Hard responsibilities
  # -------------------------
  Tau0 <- matrix(0, n, K)
  Tau0[cbind(1:n, clusters)] <- 1
  
  # -------------------------
  # Back-map initial parameters per class
  # Only update baseline intercept/slope;
  # keep gamma at pooled value
  # -------------------------
  S00 <- crossprod(X0w)   # 2 x 2
  
  beta_list   <- vector("list", K)
  Psi_list    <- vector("list", K)
  sigma2_vec  <- numeric(K)
  pi_vec      <- numeric(K)
  refine_iter <- integer(K)
  
  for (k in seq_len(K)) {
    idx <- which(clusters == k)
    n_k <- length(idx)
    pi_vec[k] <- ifelse(n_k > 0, n_k / n, 1e-12)
    
    # -------------------------
    # build a class-specific beta seed from kmeans center
    # -------------------------
    if (n_k > 0) {
      rbar_w_k  <- colMeans(Fmat[idx, , drop = FALSE])
      delta_b_k <- solve(S00, crossprod(X0w, rbar_w_k))
      beta0_k   <- beta0_pool + as.numeric(delta_b_k)
    } else {
      beta0_k <- beta0_pool
    }
    
    beta0_k <- beta0_k +
      rnorm(length(beta0_k), 0, jitter * pmax(1e-6, abs(beta0_k) + 1))
    
    if (!has_direct) {
      beta_seed_k <- beta0_k
    } else {
      # gamma_pool is only used as an initial seed here
      beta_seed_k <- c(beta0_k, gamma_pool)
    }
    
    # -------------------------
    # fallback: empty or tiny cluster
    # -------------------------
    if (!refine_classes || n_k < min_class_n) {
      beta_list[[k]]  <- beta_seed_k
      Psi_list[[k]]   <- Psi_pool
      sigma2_vec[k]   <- sigma2_pool
      refine_iter[k]  <- 0L
      next
    }
    
    # -------------------------
    # short hard-weighted class-specific LCM refinement
    # -------------------------
    fit_k <- tryCatch(
      fit_lcm(
        y = y_names,
        ID = ID,
        data = data,
        growth_covariates = growth_covariates,
        cov_pred = cov_pred,
        weights = Tau0[, k],          # hard class weights from kmeans
        lcm_design = lcm_design,
        max_iter = refine_max_iter,
        tol = refine_tol,
        penalty = FALSE,
        verbose = FALSE,
        Debug = FALSE,
        beta_init = matrix(beta_seed_k, ncol = 1),
        theta_init = Psi_pool,
        sigma2_init = sigma2_pool
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit_k)) {
      beta_list[[k]]  <- beta_seed_k
      Psi_list[[k]]   <- Psi_pool
      sigma2_vec[k]   <- sigma2_pool
      refine_iter[k]  <- 0L
    } else {
      beta_list[[k]]  <- as.numeric(fit_k$beta)
      Psi_list[[k]]   <- as.matrix(fit_k$theta)
      sigma2_vec[k]   <- as.numeric(fit_k$sigma2)
      refine_iter[k]  <- fit_k$iter
    }
  }
  
  # -------------------------
  # Normalize pi
  # -------------------------
  pi_vec <- pmax(pi_vec, 1e-12)
  pi_vec <- pi_vec / sum(pi_vec)
  

  
  params0 <- list(
    K = K,
    pi = pi_vec,
    beta = beta_list,
    Psi = Psi_list,
    sigma2 = sigma2_vec
  )
  
  list(
    params     = params0,
    Tau0       = Tau0,
    clusters   = clusters,
    centers    = centers,
    Zi         = Zi,
    R          = R,
    Omega_pool = Omega_pool,
    refine_iter = refine_iter,
    pool_fit   = list(
      beta   = beta_full_pool,
      theta  = Psi_pool,
      sigma2 = sigma2_pool
    )
  )
}


LG_fit <- function(
    param_path,
    tau_path = NULL,
    K,
    p_class,
    p_growth = p_class,
    ref_ext = 1L,
    active = NULL
){
  # -----------------------------
  # Part 1: read parameter txt
  # -----------------------------
  vals <- scan(param_path, what = numeric(), sep = "\n", dec = ".", quiet = TRUE)
  
  # Active positions follow the internal beta order, excluding the two means.
  if (is.null(active)) {
    active <- replicate(K, rep(1L, 2L * p_growth), simplify = FALSE)
  }
  if (is.list(active)) {
    if (length(active) != K || any(lengths(active) != 2L * p_growth)) {
      stop("active must contain K vectors, each of length 2 * p_growth.")
    }
    active_mat <- matrix(unlist(active), nrow = K, ncol = 2L * p_growth,
                         byrow = TRUE)
  } else {
    active_mat <- as.matrix(active)
  }
  if (!identical(dim(active_mat), c(as.integer(K), 2L * as.integer(p_growth))) ||
      anyNA(active_mat) || any(!active_mat %in% c(0, 1))) {
    stop("active must be a K x (2 * p_growth) matrix of zeros and ones.")
  }
  
  active_i <- active_mat[, seq_len(p_growth), drop = FALSE]
  active_s <- active_mat[, p_growth + seq_len(p_growth), drop = FALSE]
  
  expected_n <- (K - 1) +                # gate intercepts
    (K - 1) * p_class +                  # gate slopes
    K +                                 # latent intercept means
    sum(active_i) +                      # active DEs on latent intercept
    K +                                 # latent slope means
    sum(active_s) +                      # active DEs on latent slope
    K * 4                               # sigma2 + 3 Psi params per class
  
  if (length(vals) != expected_n) {
    stop("The txt parameter count does not match active: ",
         length(vals), " values found; ", expected_n, " expected.")
  }
  
  if (length(ref_ext) != 1L || is.na(ref_ext) ||
      !ref_ext %in% seq_len(K)) {
    stop("Please provide a valid ref_ext between 1 and K.")
  }
  
  i <- 0L
  take <- function(n) {
    idx <- i + seq_len(n)
    i <<- i + n
    vals[idx]
  }
  
  # 1) gate intercepts from external software
  alpha_raw <- take(K - 1)
  
  # 2) gate slopes from external software
  Beta_raw <- matrix(
    take((K - 1) * p_class),
    nrow = K - 1,
    ncol = p_class,
    byrow = FALSE
  )
  
  # -------------------------------------------------
  # convert external gate to internal reference = 1
  # external rows are assumed ordered by class labels
  # excluding ref_ext, in ascending order
  # -------------------------------------------------
  ext_classes <- setdiff(seq_len(K), ref_ext)
  
  a_all <- numeric(K)
  B_all <- matrix(0, nrow = K, ncol = p_class)
  
  a_all[ext_classes] <- alpha_raw
  B_all[ext_classes, ] <- Beta_raw
  
  # internal gate: class 2..K relative to class 1
  nonref <- seq_len(K)[-1L]
  alpha <- a_all[nonref] - a_all[1]
  Beta_gate <- B_all[nonref, , drop = FALSE] -
    matrix(B_all[1, ], nrow = K - 1, ncol = p_class, byrow = TRUE)
  
  # 3) latent intercept means
  mu_i <- take(K)
  
  # 4) DEs on latent intercept
  # Logical matrix indexing follows column order: covariate, then class.
  de_i <- matrix(0, nrow = K, ncol = p_growth)
  de_i[active_i == 1L] <- take(sum(active_i))
  
  # 5) latent slope means
  mu_s <- take(K)
  
  # 6) DEs on latent slope
  de_s <- matrix(0, nrow = K, ncol = p_growth)
  de_s[active_s == 1L] <- take(sum(active_s))
  
  # 7) sigma2 and Psi
  sigma2 <- numeric(K)
  Psi <- vector("list", K)
  
  for (k in seq_len(K)) {
    sigma2[k] <- take(1)
    v <- take(3)  # var(intercept), cov, var(slope)
    
    Psi[[k]] <- matrix(
      c(v[1], v[2],
        v[2], v[3]),
      nrow = 2,
      byrow = TRUE
    )
  }
  
  # 8) beta
  # internal order required by prep_lcm_design / update_tau:
  # [mu_i, mu_s, DEs on intercept, DEs on slope]
  beta <- lapply(seq_len(K), function(k) {
    c(mu_i[k], mu_s[k], de_i[k, ], de_s[k, ])
  })
  
  # 9) gate probabilities at raw covariates = 0
  # internal reference = class 1; these are not sample-average proportions.
  eta0 <- c(0, alpha)
  eta0 <- eta0 - max(eta0)
  pi_at_zero <- exp(eta0) / sum(exp(eta0))
  
  params <- list(
    K = K,
    pi = pi_at_zero,
    beta = beta,
    Psi = Psi,
    sigma2 = sigma2,
    gate = list(
      alpha = alpha,
      Beta = Beta_gate
    )
  )
  
  # -----------------------------
  # Part 2: read posterior probs
  # -----------------------------
  Tau <- NULL
  ID_vec <- NULL
  pi_source <- "gate_at_raw_zero"
  
  if (!is.null(tau_path)) {
    tau_dat <- read.csv(tau_path, check.names = FALSE, stringsAsFactors = FALSE)
    names(tau_dat) <- trimws(names(tau_dat))
    
    class_cols <- paste0("Class#", seq_len(K))
    
    if (!all(class_cols %in% names(tau_dat))) {
      stop(
        "Some posterior probability columns are missing. Needed: ",
        paste(class_cols, collapse = ", ")
      )
    }
    
    Tau <- as.matrix(tau_dat[, class_cols, drop = FALSE])
    storage.mode(Tau) <- "double"
    
    rs <- rowSums(Tau)
    if (any(!is.finite(rs))) {
      stop("Tau contains non-finite values.")
    }
    if (any(abs(rs - 1) > 1e-4)) {
      warning("Some rows of Tau do not sum to 1 within tolerance.")
    }
    
    id_idx <- which(names(tau_dat) == "ID")
    if (length(id_idx) >= 1L) {
      ID_vec <- tau_dat[[id_idx[1L]]]
    }
    
    # Posterior mean proportions, using the same LG class labels as params.
    params$pi <- colMeans(Tau)
    pi_source <- "mean_posterior"
  }
  
  # -----------------------------
  # Output
  # -----------------------------
  LG_fit <- list(
    params = params,
    Tau = Tau,
    active = lapply(seq_len(K), function(k) as.integer(active_mat[k, ])),
    param_scale = "raw",
    pi_at_zero = pi_at_zero,
    pi_source = pi_source
  )
  
  if (!is.null(ID_vec)) {
    LG_fit$ID <- ID_vec
  }
  
  LG_fit
}

