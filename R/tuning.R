# this file contains the process of turning the penalty parameter lambda
tune_psr <- function(
    data,
    y,
    ID,
    K,
    class_covariates = NULL,
    growth_covariates = NULL,
    cov_pred = "both",
    prefit = NULL,                 # must be a fitted object, not params only
    init_ps = NULL,
    nlambda = 20,
    lambda_min_ratio = 1e-3,
    lambda_expand = 1.5,
    a = 0.5, # weights of penalty
    pre_lambdas = NULL,
    warm_start = TRUE,
    seed = NULL,
    lcm_design = NULL, active_tol = 1e-04
){
  
  # -----------------------------------
  # step 1: prepare and get baseline fit
  # -----------------------------------
  # design objects
  gmm_setup()
  data_design <- lcm_design
  if (is.null(data_design)) 
    data_design <- prep_lcm_design(
      y = y, 
      data = data, 
      growth_covariates = growth_covariates, 
      cov_pred = cov_pred, 
      penalty = TRUE, 
      class_covariates = class_covariates)
  n      <- data_design$n
  X0     <- data_design$X_0
  X1_pen <- data_design$X_1_pen
  y_vec  <- data_design$y_vec
  Zi     <- data_design$z_i
  
  if (ncol(X1_pen) == 0) {
    stop("X_1_pen has zero columns, no penalized direct effects to tune.")
  }
  
  
  if (is.null(prefit)) {
    #active_zero <- replicate(K, rep(0L, ncol(X1_pen)), simplify = FALSE)
    fit_fulldie <- fit_gmm(
      y = y, ID = ID, data = data, K = K, max_iter = 1000, tol = 1e-06, 
      nstart = 10, seed = seed, verbose = T, Debug = F, a = a, class_covariates = class_covariates, 
      growth_covariates = growth_covariates, cov_pred = cov_pred, free_mode = "pi_beta_sigma2_theta", 
      penalty = FALSE, init_params = init_ps, active = NULL, Lambda = NULL, lcm_design = data_design)
    
  } else {
    fit_fulldie <- prefit
  }
  

  # extract fitted objects
  tau_unpen <- fit_fulldie$Tau
  psi_unpen <- fit_fulldie$params$Psi
  sigma2_unpen <- fit_fulldie$params$sigma2

  
  # -----------------------------------
  # step 2: compute lambda max and lambda sequence
  # -----------------------------------
  
  lam_info <- list(lambda_max_global = NA_real_, lambda_max_k = rep(NA_real_, K))
  if (!is.null(pre_lambdas)) {
    lambda_seq <- pre_lambdas
  }
  else {
    lam_info <- lambda_max(
      y = y_vec, 
      X0 = X0, 
      X1 = X1_pen, 
      tau = tau_unpen, 
      psi = psi_unpen, 
      sigma2 = sigma2_unpen, 
      Zi = Zi, 
      a = a)
    lam_max <- lam_info$lambda_max_global
    
    if (!is.finite(lam_max) || lam_max <= 0) 
      stop("Computed lambda_max_global is not a positive finite number.")
    lambda_seq <- exp(seq(
      from = log(lam_max * lambda_min_ratio), 
      to = log(lam_max * lambda_expand), 
      length.out = nlambda))
  }
  
  
  # -----------------------------------
  # step 3: run lambda path
  # keep the smallest lambda within each consecutive active-set run
  # -----------------------------------
  
  path_summary <- vector("list", length(lambda_seq))
  refit_candidates <- list()
  
  pen_start <- if (is.null(init_ps)) 
    fit_fulldie$params
  else init_ps
  prev_sig <- NULL
  current_candidate <- NULL
  last_j <- 0L
  
  make_sig <- function(active) {
    paste(
      vapply(active, function(x) paste(as.integer(x), collapse = ""), character(1)),
      collapse = "|"
    )
  }
  
  for (j in seq_along(lambda_seq)) {
    
    lam_j <- lambda_seq[j]
    cat(sprintf("\n[path %02d/%02d] lambda = %.6f\n", j, length(lambda_seq), lam_j))
    
    fit_pen <- fit_gmm(
      y = y,
      ID = ID,
      data = data,
      K = K,
      max_iter = 1000,
      tol = 1e-6,
      nstart = 1,
      seed = seed,
      verbose = T,
      Debug = FALSE,
      a = a,
      class_covariates = class_covariates,
      growth_covariates = growth_covariates,
      cov_pred = cov_pred,
      free_mode = "pi_beta_sigma2_theta",
      penalty = TRUE,
      init_params =  if (warm_start) pen_start else fit_fulldie$params,
      active = NULL,
      Lambda = rep(lam_j, K), 
      active_tol = active_tol, 
      lcm_design = data_design
    )
    
    active_j <- fit_pen$active
    print(active_j)
    sig_j <- make_sig(active_j)
    zero_check <- all(as.integer(unlist(active_j)) == 0L)
    
    # light summary for every lambda
    path_summary[[j]] <- list(
      index = j,
      lambda = lam_j,
      signature = sig_j,
      active = active_j
    )
    last_j <- j
    
    # always update penalized warm start
    pen_start <- fit_pen$params
    
    # online update of refit candidate
    if (is.null(prev_sig)) {
      
      current_candidate <- list(
        index = j,
        lambda = lam_j,
        signature = sig_j,
        active = active_j,
        params = fit_pen$params
      )
      
    } else if (identical(sig_j, prev_sig)) {

      
    } else {
      
      # close previous run
      refit_candidates[[length(refit_candidates) + 1]] <- current_candidate
      
      # start new run
      current_candidate <- list(
        index = j,
        lambda = lam_j,
        signature = sig_j,
        active = active_j,
        params = fit_pen$params
      )
    }
    
    prev_sig <- sig_j
    # early stop rule:
    # once active set is all zero, stop the path
    if (zero_check) {
      cat("All active sets are zero. Stop lambda path early.\n")
      break
    }
  }
  # save the last candidate
  if (!is.null(current_candidate)) {
    refit_candidates[[length(refit_candidates) + 1]] <- current_candidate
  }
  if (last_j < length(path_summary)) {
    path_summary <- path_summary[seq_len(last_j)]
  }
  
  
  # -----------------------------------
  # step 4: refit + compute BIC
  # -----------------------------------
  refit_fits <- vector("list", length(refit_candidates))
  ll_vec <- numeric(length(refit_candidates))
  p_active_vec <- numeric(length(refit_candidates))
  bic_vec <- numeric(length(refit_candidates))
  
  for (m in seq_along(refit_candidates)) {
    
    cand <- refit_candidates[[m]]
    cat(sprintf("\n[refit %02d/%02d] lambda = %.6f\n",
                m, length(refit_candidates), cand$lambda))
    
    fit_refit <- fit_gmm(
      y = y, ID = ID, data = data, K = K, max_iter = 1000, tol = 1e-06, 
      nstart = 1, seed = NULL, verbose = FALSE, Debug = FALSE, a = a, class_covariates = class_covariates, 
      growth_covariates = growth_covariates, cov_pred = cov_pred, free_mode = "pi_beta_sigma2_theta", 
      penalty = FALSE, init_params = cand$params, active = cand$active, Lambda = NULL, lcm_design = data_design)
    
    
    ll_refit <- tail(fit_refit$loglik_hist, 1)
    p_active <- sum(unlist(cand$active))
    BIC_pen  <- -2 * ll_refit + log(n) * p_active
    
    fit_refit$BIC <- list(
      ll = ll_refit,
      n = n,
      p_active = p_active,
      BIC_pen = BIC_pen
    )
    
    refit_fits[[m]] <- fit_refit
    ll_vec[m] <- ll_refit
    p_active_vec[m] <- p_active
    bic_vec[m] <- BIC_pen
  }
  
  best_index <- which.min(bic_vec)
  best_fit <- refit_fits[[best_index]]
  best_lambda <- refit_candidates[[best_index]]$lambda
  
  refit_summary <- data.frame(
    index = seq_along(refit_candidates),
    lambda = sapply(refit_candidates, `[[`, "lambda"),
    signature = sapply(refit_candidates, `[[`, "signature"),
    p_active = p_active_vec,
    loglik = ll_vec,
    BIC_pen = bic_vec
  )
  
  
  out <- list(
    covariate_scaling = data_design$covariate_scaling, parameter_scale = "standardized", 
    fit_fulldie = fit_fulldie, prep = data_design, lambda_max_global = lam_info$lambda_max_global, 
    lambda_max_k = lam_info$lambda_max_k, lambda_seq = lambda_seq, path_summary = path_summary, refit_candidates = refit_candidates, 
    refit_fits = refit_fits, refit_summary = refit_summary, best_index = best_index, best_lambda = best_lambda, 
    best_fit = best_fit)

  
  return(out)
}


tune_emm <- function(
    data,
    y,
    ID,
    K,
    class_covariates = NULL,
    growth_covariates = NULL,
    cov_pred = "both",
    prefit = NULL,                 # must be a fitted object, not params only
    nlambda = 20,
    lambda_min_ratio = 1e-3,
    lambda_expand = 1.5,
    a = 0.5, # weights of penalty
    warm_start = TRUE,
    seed = NULL,
    pre_lambdas = NULL,
    Debug = FALSE,
    lcm_design = NULL, active_tol = 1e-04
){
  if(Debug) browser()
  # -----------------------------------
  # step 1: prepare and get baseline fit
  # -----------------------------------
  # design objects
  gmm_setup()
  data_design <- lcm_design
  if (is.null(data_design)) 
    data_design <- prep_lcm_design(
      y = y, data = data, growth_covariates = growth_covariates, 
      cov_pred = cov_pred, penalty = TRUE, class_covariates = class_covariates)
  n      <- data_design$n
  X0     <- data_design$X_0
  X1_pen <- data_design$X_1_pen
  y_vec  <- data_design$y_vec
  Zi     <- data_design$z_i
  
  if (ncol(X1_pen) == 0) {
    stop("X_1_pen has zero columns: no penalized direct effects to tune.")
  }
  
  
  if (is.null(prefit)) {
    active_zero <- replicate(K, rep(0L, ncol(X1_pen)), simplify = FALSE)
    fit_nodie <- fit_emm(
      y = y,
      ID = ID,
      data = data,
      K = K,
      max_iter = 1000, 
      tol = 1e-6,
      nstart = 10,
      seed = seed,
      verbose = 1,
      Debug = FALSE,
      a = a,
      class_covariates = class_covariates,
      growth_covariates = growth_covariates,
      cov_pred = cov_pred,
      free_mode = "pi_beta_sigma2_theta",
      penalty = FALSE,
      init_params = NULL,
      active = NULL,
      
      emm = FALSE,
      emm_beta_only = FALSE,
      emm_m2_max_iter = 500,
      emm_m2_tol = 1e-6,
      
      Lambda = NULL,
      lcm_design = data_design
    )
  } else {
    fit_nodie <- prefit
  }
  
  
  # extract fitted objects
  tau_unpen <- fit_nodie$Tau
  psi_unpen <- fit_nodie$params$Psi
  sigma2_unpen <- fit_nodie$params$sigma2
  
  
  # -----------------------------------
  # step 2: compute lambda max and lambda sequence
  # -----------------------------------
  lam_info <- list(lambda_max_global = NA_real_, lambda_max_k = rep(NA_real_, K))
  
  if (!is.null(pre_lambdas)) {
    lambda_seq <- pre_lambdas
  }
  else {
    lam_info <- lambda_max(
      y = y_vec, X0 = X0, 
      X1 = X1_pen, tau = tau_unpen, 
      psi = psi_unpen, sigma2 = sigma2_unpen, 
      Zi = Zi, a = a)
    lam_max <- lam_info$lambda_max_global
    
    if (!is.finite(lam_max) || lam_max <= 0) 
      stop("Computed lambda_max_global is not a positive finite number.")
    
    lambda_seq <- exp(seq(
      from = log(lam_max * lambda_min_ratio), 
      to = log(lam_max * lambda_expand), 
      length.out = nlambda))
  }
  
  # -----------------------------------
  # step 3: run EMM lambda path
  # keep ALL lambda fits, no refit
  # -----------------------------------
  
  path_summary <- vector("list", length(lambda_seq))
  path_fits    <- vector("list", length(lambda_seq))
  
  pen_start <- fit_nodie$params
  
  make_sig <- function(active) {
    paste(vapply(active, function(x) paste(as.integer(x), collapse = ""), character(1)), collapse = "|")
  }
  
  # optional: constant part of parameter count for a fuller BIC
  r <- ncol(Zi)
  p_theta <- r * (r + 1) / 2
  
  p_gate <- if (is.null(class_covariates)) {
    K - 1
  } else {
    (K - 1) * (1 + length(class_covariates))
  }
  
  p_const <- K * (ncol(X0) + p_theta + 1) + p_gate
  
  for (j in seq_along(lambda_seq)) {
    
    lam_j <- lambda_seq[j]
    cat(sprintf("\n[path %02d/%02d] lambda = %.6f\n", j, length(lambda_seq), lam_j))
    
    fit_pen <- fit_emm(
      y = y,
      ID = ID,
      data = data,
      K = K,
      max_iter = 1000,
      tol = 1e-6,
      nstart = 1,
      seed = seed,
      verbose = T,
      Debug = FALSE,
      a = a,
      class_covariates = class_covariates,
      growth_covariates = growth_covariates,
      cov_pred = cov_pred,
      free_mode = "pi_beta_sigma2_theta",
      penalty = TRUE,
      init_params = if (warm_start) pen_start else fit_nodie$params,
      active = NULL,
      emm = TRUE,
      emm_beta_only = FALSE,
      emm_m2_max_iter = 500,
      emm_m2_tol = 1e-6,
      Lambda = rep(lam_j, K),
      active_tol = active_tol, lcm_design = data_design
    )
    
    # extract active set / signature
    active_j <- fit_pen$active
    sig_j <- make_sig(active_j)
    zero_check <- all(as.integer(unlist(active_j)) == 0L)
    
    # observed log-likelihood from this EMM fit
    ll_j <- tail(fit_pen$loglik_hist, 1)
    
    # number of selected direct effects
    p_active_j <- sum(unlist(active_j))
    
    # two BIC versions
    BIC_active_only_j <- -2 * ll_j + log(n) * p_active_j
    BIC_full_j        <- -2 * ll_j + log(n) * (p_const + p_active_j)
    
    # store BIC inside fit object as well
    fit_pen$BIC <- list(
      ll = ll_j,
      n = n,
      p_active = p_active_j,
      p_const = p_const,
      BIC_active_only = BIC_active_only_j,
      BIC_full = BIC_full_j
    )
    
    # light summary for every lambda
    path_summary[[j]] <- data.frame(
      index = j,
      lambda = lam_j,
      signature = sig_j,
      logLik = ll_j,
      p_active = p_active_j,
      p_const = p_const,
      BIC_active_only = BIC_active_only_j,
      BIC_full = BIC_full_j,
      stringsAsFactors = FALSE
    )
    
    # store full fit
    path_fits[[j]] <- fit_pen
    
    # warm start for next lambda
    pen_start <- fit_pen$params
    
    if (zero_check) {
      cat("All active sets are zero. Stop lambda path early.\n")
      break
    }
  }
  
  path_summary <- do.call(rbind, path_summary)
  
  # best model by BIC_full
  best_idx <- which.min(path_summary$BIC_full)
  
  out <- list(
    covariate_scaling = data_design$covariate_scaling, 
    parameter_scale = "standardized", 
    fit_nodie = fit_nodie, 
    prep = data_design, 
    lambda_max_global = lam_info$lambda_max_global, 
    lambda_max_k = lam_info$lambda_max_k, 
    lambda_seq = lambda_seq, 
    path_summary = path_summary, 
    path_fits = path_fits, 
    best_index = best_idx, 
    best_lambda = path_summary$lambda[best_idx], 
    best_fit = path_fits[[best_idx]])
  
  return(out)
}


lambda_max <- function(
    y = NULL,         # nt x 1 stacked outcome
    X0 = NULL,        # nt x 2 unpenalized design matrix
    X1 = NULL,        # nt x p penalized design matrix
    Zi = NULL,        # T x 2 growth loading matrix
    tau = NULL,       # n x K posterior matrix
    psi = NULL,       # list of K q x q random-effect covariance matrices
    sigma2 = NULL,    # length-K residual variance vector
    a = NULL
){

  stopifnot(!is.null(y), !is.null(X1), !is.null(Zi),
            !is.null(tau), !is.null(psi), !is.null(sigma2))
  
  y <- as.matrix(y)
  K_class <- ncol(tau)
  n <- nrow(tau)
  Tobs <- nrow(Zi)
  
  # construct X0 = stacked Zi for all subjects
  #X0 <- do.call(rbind, replicate(n, Zi, simplify = FALSE))  # nt x 2
  
  # dimension checks
  stopifnot(nrow(X0) == n * Tobs)
  stopifnot(nrow(X1) == n * Tobs)
  stopifnot(nrow(y)  == n * Tobs)
  
  lambda_k <- numeric(K_class)
  
  for (k in seq_len(K_class)) {
    # class-specific marginal covariance
    omega_k <- Zi %*% psi[[k]] %*% t(Zi) + as.numeric(sigma2[k]) * diag(Tobs)
    omega_k_w <- kronecker(diag(tau[, k]), pracma::pinv(omega_k))
    
    # reduced fit with all direct effects fixed to 0
    XtWX <- crossprod(X0, omega_k_w %*% X0)
    XtWy <- crossprod(X0, omega_k_w %*% y)
    beta_k <- MASS::ginv(XtWX) %*% XtWy
    
    # reduced residual
    r_k <- y - X0 %*% beta_k
    
    # score vector for penalized block
    score_k <- as.vector(crossprod(X1, omega_k_w %*% r_k))
    
    # class-specific lambda max
    lambda_k[k] <- max(abs(score_k))
  }
  
  lambda_max_global <- max(lambda_k)
  
  return(list(
    lambda_max_global = lambda_max_global / a,
    lambda_max_k = lambda_k / a
  ))
}
