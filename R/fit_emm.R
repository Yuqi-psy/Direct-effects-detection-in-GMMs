# freeze parameters
# EMM- to reduce the penalty bias, we maximise two M step at each EM iteration.
# We refer the first M steo as the M_pem step, in which we specifiy all possible DEs
#and impose the penalty term on the direct effects
# We refer the second M step as the M_nonpen step, we specify the active set of DEs retrived from the M_pem step
# and maximise the Q function without impose penalty.



# freeze parameters
fit_emm <- function(
    y, 
    ID, 
    data, 
    K,
    penalty = FALSE,
    Lambda = NULL,
    a = 0.5,
    active = NULL,
    
    class_covariates = NULL,   # specify covs name
    growth_covariates = NULL,
    cov_pred = NULL,
    max_iter = 100, 
    tol = 1e-5,
    nstart = 10, 
    seed = NULL,
    init_params = NULL,
    
    verbose = TRUE,   # 0/1/2 
    Debug = FALSE,
    
    emm = FALSE,
    emm_beta_only = FALSE,
    emm_m2_max_iter = 500,
    emm_m2_tol = 1e-6,
    
    free_mode = "all",  # "pi_only","pi_beta","pi_sigma2","pi_theta",
    # "pi_beta_sigma2","pi_beta_theta","pi_sigma2_theta",
    # "pi_beta_sigma2_theta"(=all),"all"
    lcm_design = NULL, active_tol = 1e-04
){
  if (!is.null(seed)) set.seed(seed)
  if(Debug) browser( )
  gmm_setup() 
  
  verbose_level <- if (isTRUE(verbose)) 1L else if (isFALSE(verbose)) 0L else as.integer(verbose)
  if (is.na(verbose_level)) verbose_level <- 1L
  logf <- function(level, fmt, ...) {
    if (verbose_level >= level) cat(sprintf(fmt, ...))
  }
  sigma_flag   <- 1e-6
  ll_jump_flag <- 100
  lcm_iter_cap <- 3000L   # keep consistent with LCM_estimation_weight_freeze(max_iter = 3000)
  
  boundary_flag_hist   <- logical(max_iter)
  boundary_class_hist  <- vector("list", max_iter)
  boundary_reason_hist <- vector("list", max_iter)
  boundary_info        <- NULL # history
  
  # -- basic shapes and design matrices --
  n    <- nrow(data)
  Tobs <- length(y)
  Y    <- as.matrix(data[, y, drop = FALSE])            # n x T
  Xi   <- cbind(1, 0:(Tobs - 1))                        # T x p_x
  Zi   <- Xi                                            # T x p_z 
  # ----  growth-covariate design for within-class mean ----
  if (is.null(cov_pred)) cov_pred <- "none"
  cov_pred <- match.arg(cov_pred, c("none","intercept","slope","both"))
  # init_params and all internal/output coefficients are STANDARDIZED.
  if (is.null(lcm_design)) {
    lcm_design <- prep_lcm_design(
      y = y,
      data = data,
      growth_covariates = growth_covariates,
      cov_pred = cov_pred,
      penalty = penalty,
      class_covariates = class_covariates
    )
  }
  data <- lcm_design$data_std
  
  C <- if (!is.null(growth_covariates) && length(growth_covariates) > 0L && cov_pred != "none") {
    lcm_design$C   # n x q
  } else {
    NULL
  }
  
  if (is.null(class_covariates) || length(class_covariates) == 0L) {
    use_covariates <- FALSE
    X_class <- NULL
  } else {
    use_covariates <- TRUE
    X_class <- lcm_design$X_class # n x q_x
  } # design matrix with class-membership covariates ---
  
  
  # -- step 0: pooled fit -> whitening K-means init --
  
  
  # input initial valus
  if (!is.null(init_params)){
    init <- init_params
    params <- init
  } else {
    init <- kmeans_init(y_names = y, ID = ID, data = data, 
                        growth_covariates = growth_covariates,           # covs names
                        cov_pred = cov_pred,
                        K = K, pool_fit = NULL,
                        refine_classes = TRUE, refine_max_iter = 10,
                        refine_tol = 1e-3,
                        use_residuals = TRUE, nstart = nstart,
                        jitter = 0.02, seed = seed, Debug = Debug, 
                        lcm_design = lcm_design)
    params <- init$params
  }
  
  if (penalty) {
    params$Lambda <- Lambda
  } else {
    params$Lambda <- NULL
  }
  
  # initial beta to match LCM design when growth_covariates are used
  if (!is.null(growth_covariates) && length(growth_covariates) > 0L && cov_pred != "none") {
    q <- ncol(C)
    p_target <- 2L + if (cov_pred == "both") 2L*q else q   # 2 + q or 2 + 2q
    
    params$beta <- lapply(params$beta, function(b){
      b <- as.numeric(b)
      length(b) <- p_target   
      b[is.na(b)] <- 0       
      b
    })
  }
  
  
  if (use_covariates&& (is.null(params$gate))) {
    q <- ncol(X_class)
    
    pi0 <- pmax(as.numeric(params$pi), 1e-8)
    pi0 <- pi0/sum(pi0)
    
    params$gate <- list(
      alpha = log(pi0[2:K] / pi0[1]),                    # logit intercept α_2,...,α_K
      Beta  = matrix(0, nrow = K - 1, ncol = q) # logit slope β_2,...,β_K
    )
  } 
  
  if (is.null(init_params) && !is.null(init$Tau0) && use_covariates) {
    pi_ik0 <- update_pi_ik(X_class, params$gate)
    
    gate0 <- update_gate(
      X_class = X_class,
      pi_ik   = pi_ik0,
      params  = params,
      tau     = init$Tau0,
      verbose = FALSE,
      Debug   = Debug
    )
    
    params$gate <- gate0$params$gate
  }
  
  
  # final initial conditional mixing proportions
  pi_ik0 <- if (use_covariates) {
    update_pi_ik(X_class, params$gate)
  } else {
    NULL
  } # initialize mixing proportions and gate parameters
  
  
  
  # -- optional: compute initial loglik --
  ll0 <- loglik_GMM(Y = Y, Xi = Xi, Zi = Zi, 
                    C = C, cov_pred = cov_pred, 
                    params = params, include_const = TRUE,
                    pi_ik = pi_ik0)            # NEW input pi_ik
  logf(1, "[Init] loglik = %.6f\n", ll0)        # [KEEP]
  
  # ---- histories ----
  loglik_hist   <- numeric(max_iter)            # [KEEP] M-step (observed) loglik history
  # ===== keep M1 / M2 iteration histories separately =====
  lcm_iter_hist_m1 <- matrix(NA_integer_, nrow = max_iter, ncol = K)
  lcm_iter_hist_m2 <- matrix(NA_integer_, nrow = max_iter, ncol = K)
  
  # backward-compatible combined history (optional)
  lcm_iter_hist <- matrix(NA_integer_, nrow = max_iter, ncol = K)
  
  if (penalty) {
    penalty_total_hist <- numeric(max_iter)
    loglik_pen_hist    <- numeric(max_iter)
  }
  
  
  #////////////////  free_mode  //////////////////////
  valid_modes <- c(
    "pi_only",
    "pi_beta", "pi_sigma2", "pi_theta",
    "pi_beta_sigma2", "pi_beta_theta", "pi_sigma2_theta",
    "pi_beta_sigma2_theta",  # (= all)
    "all"                    # alias
  )
  if (!free_mode %in% valid_modes) {
    stop(sprintf("free_mode must be one of: %s", paste(valid_modes, collapse = ", ")))
  }
  
  if (free_mode == "all") free_mode <- "pi_beta_sigma2_theta"
  
  update_pi <- TRUE
  allow_beta   <- grepl("beta",   free_mode, fixed = TRUE)
  allow_sigma2 <- grepl("sigma2", free_mode, fixed = TRUE)
  allow_theta  <- grepl("theta",  free_mode, fixed = TRUE)
  if (free_mode == "pi_only") {
    allow_beta <- allow_sigma2 <- allow_theta <- FALSE
  }
  
  logf(1,
       "[Mode] %s | update_pi=TRUE | allow_beta=%s allow_sigma2=%s allow_theta=%s\n",
       free_mode, allow_beta, allow_sigma2, allow_theta)
  logf(1, "[EMM] emm=%s | emm_beta_only=%s | m2_max_iter=%d | m2_tol=%.1e\n",
       emm, emm_beta_only, emm_m2_max_iter, emm_m2_tol)
  #////////////////////// end free_mode //////////////////////
  
  # -- EM loop --
  ll_prev  <- ll0                
  obj_prev <- ll0
  active_hist <- vector("list", max_iter)
  for (iter in 1:max_iter) {
    
    # ========= E-step =========
    params_old <- params
    
    pi_ik <- if (use_covariates) {
      update_pi_ik(X_class, params_old$gate)  # n x K
    } else {
      NULL
    }
    
    eout <- update_tau(
      Y = Y, Xi = Xi, Zi = Zi,
      params = params_old,
      pi_ik  = pi_ik,
      C = C, cov_pred = cov_pred,   
      include_const = TRUE,
      protection = FALSE
    )
    
    Tau_old <- eout$Tau
    Nk_old  <- eout$Nk
    
    
    # ----  class-mass monitoring (useful in low separation) ----
    minNk <- min(Nk_old)
    maxNk <- max(Nk_old)
    
    # ========= M-step: mixing proportions =========
    if (use_covariates) {
      
      gout <- update_gate(
        X_class = X_class,
        pi_ik   = pi_ik,
        params  = params_old,
        tau     = Tau_old,
        verbose = FALSE,   # [MOD] suppress inner printing; we log at outer level
        Debug   = Debug
      ) # update gate parameter if covs
      params$gate <- gout$params$gate      
      pi_ik       <- gout$pi_ik           
      params$pi   <- colMeans(pi_ik)       
    } else {
      params$pi <- as.numeric(Nk_old / n)
      pi_ik <- matrix(rep(params$pi, each = n), nrow = n, ncol = K)
    }
    
    
    if (free_mode == "pi_only") {  
      ll_new <- loglik_GMM( Y = Y, Xi = Xi, Zi = Zi, 
                            C = C, cov_pred = cov_pred, 
                            params = params,  pi_ik = pi_ik,
                            include_const = TRUE)
      loglik_hist[iter] <- ll_new
      
      # ---- [KEEP] concise per-iteration report: observed ll + Δll + Nk info ----
      ll_prev_print <- if (iter == 1L) ll0 else ll_prev
      improv_print  <- ll_new - ll_prev_print
      logf(1, "[EM %03d] ll=%.6f  Δll=%+.3e  minNk=%.1f maxNk=%.1f\n",
           iter, ll_new, improv_print, minNk, maxNk)
      
      if (iter > 1) {
        improv <- ll_new - ll_prev
        if (improv < -1e-8)
          warning(sprintf("[Iter %d] Observed log-likelihood decreased by %.3e.", iter, improv))
        if (improv >= 0 && improv < tol) {
          logf(1, "Converged at iter %d: Δℓ = %.3e\n", iter, improv)  # [KEEP]
          Tau_final <- update_tau(Y = Y, Xi = Xi, Zi = Zi, 
                                  pi_ik = pi_ik, C = C, cov_pred = cov_pred, 
                                  params = params, include_const = TRUE)$Tau
          
          # trim & return
          loglik_hist   <- loglik_hist[1:iter]
          lcm_iter_hist <- lcm_iter_hist[1:iter, , drop = FALSE]
          return(list(
            covariate_scaling = lcm_design$covariate_scaling,
            parameter_scale = "standardized", active_scale = "standardized",
            active = NULL,
            params = params,
            Tau    = Tau_final,
            loglik_hist   = loglik_hist,
            lcm_iter_hist = lcm_iter_hist,
            penalty_total_hist = if (penalty) penalty_total_hist[1:iter] else NULL,
            loglik_pen_hist    = if (penalty) loglik_pen_hist[1:iter] else NULL,
            init   = init,
            iter   = iter
          ))
        }
      }
      ll_prev  <- ll_new           
      obj_prev <- ll_new           
      next
    }
    
    # ===== M1 step -- penalized fit =====
    iter_k1 <- integer(K)
    iter_k2 <- integer(K)
    penalty_total <- 0
    
    active_list <- vector("list", K)
    
    # store M1 results 
    beta_m1   <- vector("list", K)
    Psi_m1    <- vector("list", K)
    sigma2_m1 <- numeric(K)
    
    for (k in 1:K) {
      w_k <- as.numeric(Tau_old[, k])
      active_k <- if (!is.null(active)) active[[k]] else NULL
      
      if (penalty) {
        lambda_k <- as.numeric(Lambda[k])
        lambda_k <- c(
          lambda_k * a,
          lambda_k * (1 - a)
        )
      } else {
        lambda_k <- c(0, 0)
      }
      
      fit1_k <- fit_lcm(
        y = y,
        ID = ID,
        data = data,
        growth_covariates = growth_covariates,
        cov_pred = cov_pred,
        lcm_design = lcm_design,
        weights = w_k,
        active = active_k,
        active_tol = active_tol,
        max_iter = 3000,
        tol = 1e-7,
        
        penalty = penalty,
        lambda = lambda_k,
        
        Debug = Debug,
        verbose = FALSE,
        
        freeze_beta   = !allow_beta,
        freeze_theta  = !allow_theta,
        freeze_sigma2 = !allow_sigma2,
        
        beta_init   = params$beta[[k]],
        theta_init  = params$Psi[[k]],
        sigma2_init = params$sigma2[k]
      )
      
      iter_k1[k] <- if (!is.null(fit1_k$iter)) as.integer(fit1_k$iter) else NA_integer_
      
      beta_m1[[k]]   <- as.numeric(fit1_k$beta)
      Psi_m1[[k]]    <- as.matrix(fit1_k$theta)
      sigma2_m1[k]   <- as.numeric(fit1_k$sigma2)
      
      if (penalty && !is.null(fit1_k$penalty_k)) {
        penalty_total <- penalty_total + fit1_k$penalty_k$pen1 + fit1_k$penalty_k$pen2
      }
      
      if (!is.null(fit1_k$active)) {
        active_list[[k]] <- as.integer(fit1_k$active)
      } else {
        active_list[[k]] <- NULL
      }
    }
    
    if (!emm) {
      for (k in 1:K) {
        if (allow_beta)   params$beta[[k]] <- beta_m1[[k]]
        if (allow_theta)  params$Psi[[k]]  <- Psi_m1[[k]]
        if (allow_sigma2) params$sigma2[k] <- sigma2_m1[k]
      }
    }
    
    
    # ---- EMM pre-check after M1: decide whether M2 is allowed ----
    m1_hard_class <- logical(K)
    m1_warn_class <- logical(K)
    m1_hard_reason <- vector("list", K)
    m1_warn_reason <- vector("list", K)
    
    for (k in seq_len(K)) {
      if (sigma2_m1[k] < sigma_flag) {
        m1_hard_class[k] <- TRUE
        m1_hard_reason[[k]] <- c(m1_hard_reason[[k]], "sigma2_near_zero_m1")
      }
      if (!is.na(iter_k1[k]) && iter_k1[k] >= lcm_iter_cap) {
        m1_warn_class[k] <- TRUE
        m1_warn_reason[[k]] <- c(m1_warn_reason[[k]], "LCM_hit_max_iter_m1")
      }
    }
    
    m1_hard_stop <- emm && any(m1_hard_class)
    if (m1_hard_stop) {
      
      m1_msg <- paste(
        sprintf(
          "class %d: %s",
          which(m1_hard_class),
          vapply(m1_hard_reason[m1_hard_class], paste, character(1), collapse = ", ")
        ),
        collapse = " | "
      )
      
      if (any(m1_warn_class)) {
        warn_msg <- paste(
          sprintf(
            "class %d: %s",
            which(m1_warn_class),
            vapply(m1_warn_reason[m1_warn_class], paste, character(1), collapse = ", ")
          ),
          collapse = " | "
        )
        m1_msg <- paste(c(m1_msg, warn_msg), collapse = " | ")
      }
      
      logf(1, "[EM %03d][M1] Hard failure detected. Skip M2 and stop fitting.\n", iter)
      if (nzchar(m1_msg)) logf(1, "         Reason: %s\n", m1_msg)
      logf(2, "         LCM iters M1: %s\n", paste(iter_k1, collapse = ", "))
      
      boundary_flag_hist[iter] <- TRUE
      boundary_class_hist[[iter]] <- m1_hard_class | m1_warn_class
      boundary_reason_hist[[iter]] <- Map(c, m1_hard_reason, m1_warn_reason)
      
      boundary_info <- list(
        stop_iter      = iter,
        stage          = "M1",
        iter_k_m1      = iter_k1,
        sigma2_m1      = sigma2_m1,
        hard_class_m1  = m1_hard_class,
        hard_reason_m1 = m1_hard_reason,
        warn_class_m1  = m1_warn_class,
        warn_reason_m1 = m1_warn_reason,
        act_problem    = active_list,
        params_m1_last = list(
          beta   = beta_m1,
          Psi    = Psi_m1,
          sigma2 = sigma2_m1,
          active = active_list
        )
      )
      
      last_valid_iter <- iter - 1L
      
      return(list(
        covariate_scaling = lcm_design$covariate_scaling,
        parameter_scale = "standardized", active_scale = "standardized",
        active = NULL,
        params = params_old,
        Tau    = Tau_old,
        act = active_list,
        act_flag = "problem_iter_m1",
        params_flag = "previous_valid_iter",
        
        boundary_flag_hist   = boundary_flag_hist[1:iter],
        boundary_class_hist  = boundary_class_hist[1:iter],
        boundary_reason_hist = boundary_reason_hist[1:iter],
        boundary_info = boundary_info,
        
        loglik_hist = if (last_valid_iter > 0L) loglik_hist[1:last_valid_iter] else numeric(0),
        lcm_iter_hist = if (last_valid_iter > 0L) lcm_iter_hist[1:last_valid_iter, , drop = FALSE] else lcm_iter_hist[0, , drop = FALSE],
        penalty_total_hist = if (penalty && last_valid_iter > 0L) penalty_total_hist[1:last_valid_iter] else NULL,
        loglik_pen_hist    = if (penalty && last_valid_iter > 0L) loglik_pen_hist[1:last_valid_iter] else NULL,
        
        status = "boundary_stop",
        init   = init,
        iter   = last_valid_iter,
        emm = emm,
        emm_beta_only = emm_beta_only,
        lcm_iter_hist_m1 = lcm_iter_hist_m1[1:iter, , drop = FALSE],
        lcm_iter_hist_m2 = lcm_iter_hist_m2[1:max(0, iter - 1), , drop = FALSE],
        params_m1_last = list(
          beta   = beta_m1,
          Psi    = Psi_m1,
          sigma2 = sigma2_m1,
          active = active_list
        )
      ))
    }
    
    
    # ===== M2 step -- support unpenalized fit =====
    if (emm) {
      
      # keep last M1 params for later comparison
      params_m1_last <- list(
        beta = beta_m1,
        Psi = Psi_m1,
        sigma2 = sigma2_m1,
        active = active_list
      )
      
      for (k in 1:K) {
        w_k <- as.numeric(Tau_old[, k])
        
        fit2_k <- fit_lcm(
          y = y,
          ID = ID,
          data = data,
          growth_covariates = growth_covariates,
          cov_pred = cov_pred,
          lcm_design = lcm_design,
          weights = w_k,
          active = active_list[[k]],     # fixed support from M1
          active_tol = active_tol,
          max_iter = emm_m2_max_iter,
          tol = emm_m2_tol,
          
          penalty = FALSE,               # no penalty in M2
          lambda = c(0, 0),
          
          Debug = Debug,
          verbose = FALSE,
          
          
          freeze_beta   = !allow_beta,
          freeze_theta  = if (emm_beta_only) TRUE else !allow_theta,
          freeze_sigma2 = if (emm_beta_only) TRUE else !allow_sigma2,
          
          # initialize from M1 results
          beta_init   = beta_m1[[k]],
          theta_init  = Psi_m1[[k]],
          sigma2_init = sigma2_m1[k]
        )
        
        iter_k2[k] <- if (!is.null(fit2_k$iter)) as.integer(fit2_k$iter) else NA_integer_
        
        if (allow_beta)   params$beta[[k]] <- as.numeric(fit2_k$beta)
        if (!emm_beta_only && allow_theta)  params$Psi[[k]]  <- as.matrix(fit2_k$theta)
        if (!emm_beta_only && allow_sigma2) params$sigma2[k] <- as.numeric(fit2_k$sigma2)
        
        # if beta-only debiasing, keep theta/sigma2 from M1
        if (emm_beta_only) {
          if (allow_theta)  params$Psi[[k]]  <- Psi_m1[[k]]
          if (allow_sigma2) params$sigma2[k] <- sigma2_m1[k]
        }
      }
    }
    
    
    active_hist[[iter]] <- lapply(active_list, as.integer)
    
    make_sig <- function(active) {
      paste(
        vapply(active, function(x) paste(as.integer(x), collapse = ""), character(1)),
        collapse = "|"
      )
    }
    
    cat(sprintf("[EM %03d][M1] active = %s\n", iter, make_sig(active_list)))
    
    lcm_iter_hist_m1[iter, ] <- iter_k1
    lcm_iter_hist_m2[iter, ] <- if (emm) iter_k2 else rep(NA_integer_, K)
    
    # combined history for backward compatibility
    lcm_iter_hist[iter, ] <- if (emm) pmax(iter_k1, iter_k2, na.rm = TRUE) else iter_k1
    
    # ========= Evaluate current iteration (fixed τ_t) =========
    ll_new <- loglik_GMM(Y = Y, Xi = Xi, Zi = Zi, C = C, cov_pred = cov_pred, params = params,  pi_ik = pi_ik,
                         include_const = TRUE)
    ll_pen <- ll_new
    if (penalty) ll_pen <- ll_new - penalty_total
    # if EMM check convergence based on observed LL instead of penalized LL
    if (emm) {
      obj_new <- ll_new
    } else {
      obj_new <- if (penalty) ll_pen else ll_new
    }
    
  
    
    # ---- final check after current iteration ----
    ll_prev_check <- if (iter == 1L) ll0 else ll_prev
    ll_improv_check <- ll_new - ll_prev_check
    
    hard_class <- logical(K)
    warn_class <- logical(K)
    hard_reason <- vector("list", K)
    warn_reason <- vector("list", K)
    
    iter_k_check <- if (emm) pmax(iter_k1, iter_k2, na.rm = TRUE) else iter_k1
    
    for (k in seq_len(K)) {
      # 1) residual variance nearly collapsed -> hard stop
      if (params$sigma2[k] < sigma_flag) {
        hard_class[k] <- TRUE
        hard_reason[[k]] <- c(hard_reason[[k]], "sigma2_near_zero")
      }
      
      # 2) inner LCM hit iteration cap -> warning only
      if (!is.na(iter_k_check[k]) && iter_k_check[k] >= lcm_iter_cap) {
        warn_class[k] <- TRUE
        warn_reason[[k]] <- c(warn_reason[[k]], "LCM_hit_max_iter")
      }
    }
    
    # 3) abnormal observed loglik jump -> hard stop
    jump_flag <- (iter >= 5L && ll_improv_check > ll_jump_flag)
    hard_stop_flag <- any(hard_class) || jump_flag
    
    # combine hard + warning for history (keep old interface)
    boundary_solution_flag <- hard_stop_flag || any(warn_class)
    boundary_class <- hard_class | warn_class
    boundary_reason <- Map(c, hard_reason, warn_reason)
    
    boundary_flag_hist[iter] <- boundary_solution_flag
    boundary_class_hist[[iter]] <- boundary_class
    boundary_reason_hist[[iter]] <- boundary_reason
    
    # warning only: continue fitting
    if (any(warn_class)) {
      warn_msg <- paste(
        sprintf(
          "class %d: %s",
          which(warn_class),
          vapply(warn_reason[warn_class], paste, character(1), collapse = ", ")
        ),
        collapse = " | "
      )
      
      logf(1, "[EM %03d] Warning detected, continue fitting.\n", iter)
      logf(1, "         Warning: %s\n", warn_msg)
      if (emm) {
        logf(2, "         LCM iters M1: %s\n", paste(iter_k1, collapse = ", "))
        logf(2, "         LCM iters M2: %s\n", paste(iter_k2, collapse = ", "))
      } else {
        logf(2, "         LCM iters: %s\n", paste(iter_k1, collapse = ", "))
      }
    }
    
    # hard stop: return previous valid iterate, but keep current active set
    if (hard_stop_flag) {
      
      msg_parts <- c(
        if (any(hard_class)) {
          paste(
            sprintf(
              "class %d: %s",
              which(hard_class),
              vapply(hard_reason[hard_class], paste, character(1), collapse = ", ")
            ),
            collapse = " | "
          )
        },
        if (any(warn_class)) {
          paste(
            sprintf(
              "class %d: %s",
              which(warn_class),
              vapply(warn_reason[warn_class], paste, character(1), collapse = ", ")
            ),
            collapse = " | "
          )
        },
        if (jump_flag) paste0("abnormal_ll_jump=", format(ll_improv_check, scientific = TRUE, digits = 3))
      )
      
      msg_all <- paste(msg_parts[nzchar(msg_parts)], collapse = " | ")
      
      logf(1, "[EM %03d] Hard failure detected. Stop fitting.\n", iter)
      if (nzchar(msg_all)) logf(1, "         Reason: %s\n", msg_all)
      logf(1, "         Current ll=%.6f  Δll=%+.3e  minNk=%.1f maxNk=%.1f\n",
           ll_new, ll_improv_check, minNk, maxNk)
      
      if (emm) {
        logf(2, "         LCM iters M1: %s\n", paste(iter_k1, collapse = ", "))
        logf(2, "         LCM iters M2: %s\n", paste(iter_k2, collapse = ", "))
      } else {
        logf(2, "         LCM iters: %s\n", paste(iter_k1, collapse = ", "))
      }
      
      # Diagnostic only: current responsibilities may be unavailable at a boundary.
      Tau_new <- tryCatch(
        update_tau(Y = Y, Xi = Xi, Zi = Zi, params = params,
                   pi_ik = pi_ik, C = C, cov_pred = cov_pred,
                   include_const = TRUE)$Tau,
        error = function(e) NULL
      )
      boundary_info <- list(
        stop_iter      = iter,
        stage          = "final",
        ll_current     = ll_new,
        ll_prev        = ll_prev_check,
        delta_ll       = ll_improv_check,
        minNk          = minNk,
        maxNk          = maxNk,
        iter_k_m1      = iter_k1,
        iter_k_m2      = if (emm) iter_k2 else NULL,
        sigma2         = params$sigma2,
        class_flag     = boundary_class,
        class_reason   = boundary_reason,
        hard_class     = hard_class,
        hard_reason    = hard_reason,
        warn_class     = warn_class,
        warn_reason    = warn_reason,
        jump_flag      = jump_flag,
        act_problem    = active_list,
        params_problem = params,
        Tau_problem    = Tau_new
      )
      
      last_valid_iter <- iter - 1L
      
      return(list(
        covariate_scaling = lcm_design$covariate_scaling,
        parameter_scale = "standardized", active_scale = "standardized",
        active = NULL,
        params = params_old,         # previous valid iterate
        Tau    = Tau_old,            # previous valid responsibilities
        
        # keep name unchanged for downstream PSR / EMM extraction
        act = active_list,
        act_flag = "problem_iter",
        params_flag = "previous_valid_iter",
        
        # diagnostic only
        params_problem = params,
        Tau_problem    = Tau_new,
        
        loglik_hist = if (last_valid_iter > 0L) loglik_hist[1:last_valid_iter] else numeric(0),
        lcm_iter_hist = if (last_valid_iter > 0L) lcm_iter_hist[1:last_valid_iter, , drop = FALSE] else lcm_iter_hist[0, , drop = FALSE],
        penalty_total_hist = if (penalty && last_valid_iter > 0L) penalty_total_hist[1:last_valid_iter] else NULL,
        loglik_pen_hist    = if (penalty && last_valid_iter > 0L) loglik_pen_hist[1:last_valid_iter] else NULL,
        
        boundary_flag_hist   = boundary_flag_hist[1:iter],
        boundary_class_hist  = boundary_class_hist[1:iter],
        boundary_reason_hist = boundary_reason_hist[1:iter],
        boundary_info = boundary_info,
        
        status = "boundary_stop",
        init   = init,
        iter   = last_valid_iter,
        emm = emm,
        emm_beta_only = emm_beta_only,
        lcm_iter_hist_m1 = lcm_iter_hist_m1[1:iter, , drop = FALSE],
        lcm_iter_hist_m2 = lcm_iter_hist_m2[1:iter, , drop = FALSE],
        params_m1_last = if (emm) params_m1_last else NULL
      ))
    }
    #-------
    
    
    loglik_hist[iter] <- ll_new
    if (penalty && !emm) {
      penalty_total_hist[iter] <- penalty_total
      loglik_pen_hist[iter]    <- ll_pen
    }
    
    # ----  concise per-iteration report: observed ll + Δll + Nk info ----
    if (penalty) {
      obj_prev_print <- if (iter == 1L) ll0 else obj_prev
      improv_print   <- obj_new - obj_prev_print
      logf(1, "[EM %03d] ll=%.6f  obj=%.6f  Δobj=%+.3e  minNk=%.1f maxNk=%.1f\n",
           iter, ll_new, obj_new, improv_print, minNk, maxNk)
    } else {
      ll_prev_print <- if (iter == 1L) ll0 else ll_prev
      improv_print  <- ll_new - ll_prev_print
      logf(1, "[EM %03d] ll=%.6f  Δll=%+.3e  minNk=%.1f maxNk=%.1f\n",
           iter, ll_new, improv_print, minNk, maxNk)
    }
    if (emm) {
      logf(2, "         LCM iters M1: %s\n", paste(iter_k1, collapse = ", "))
      logf(2, "         LCM iters M2: %s\n", paste(iter_k2, collapse = ", "))
    } else {
      logf(2, "         LCM iters: %s\n", paste(iter_k1, collapse = ", "))
    }
    
    
    # ========= Convergence =========
    
    if (iter > 1) {
      improv <- obj_new - obj_prev
      rel_change <- 2*(obj_new - obj_prev)/(abs(obj_new) + abs(obj_prev)) 
      
      if (improv < -1e-8) {
        if (emm) {
          warning(sprintf("[Iter %d] Observed log-likelihood decreased by %.3e.", iter, improv))
        } else if (penalty) {
          warning(sprintf("[Iter %d] Penalized objective decreased by %.3e.", iter, improv))
        } else {
          warning(sprintf("[Iter %d] Observed log-likelihood decreased by %.3e.", iter, improv))
        }
      }
      
      if (improv >= 0 && rel_change < tol) {
        if (emm) {
          logf(1, "Converged at iter %d: Δℓ = %.3e\n", iter, improv)
        } else if (penalty) {
          logf(1, "Converged at iter %d: Δobj = %.3e\n", iter, improv)
        } else {
          logf(1, "Converged at iter %d: Δℓ = %.3e\n", iter, improv)
        }
      
        
        pi_ik_final <- if (use_covariates) {
          update_pi_ik(X_class, params$gate)
        } else {
          matrix(rep(params$pi, each = n), nrow = n, ncol = K)
        }
        
        Tau_final <- update_tau(
          Y = Y, Xi = Xi, Zi = Zi, params = params,
          pi_ik = pi_ik_final,
          C = C, cov_pred = cov_pred,
          include_const = TRUE
        )$Tau
        
        if(use_covariates){pi_hat_tau <- colMeans(Tau_final)} else {pi_hat_tau <- NULL}
        
        # trim & return
        loglik_hist   <- loglik_hist[1:iter]
        lcm_iter_hist <- lcm_iter_hist[1:iter, , drop = FALSE]
        return(list(
          covariate_scaling = lcm_design$covariate_scaling,
          parameter_scale = "standardized", active_scale = "standardized",
          params = params,
          Tau    = Tau_final,
          pi_hat_tau = pi_hat_tau,
          loglik_hist   = loglik_hist,
          lcm_iter_hist = lcm_iter_hist,
          penalty_total_hist = if (penalty&& !emm) penalty_total_hist[1:iter] else NULL,
          loglik_pen_hist    = if (penalty&& !emm) loglik_pen_hist[1:iter] else NULL,
          emm = emm,
          emm_beta_only = emm_beta_only,
          lcm_iter_hist_m1 = lcm_iter_hist_m1[1:iter, , drop = FALSE],
          lcm_iter_hist_m2 = lcm_iter_hist_m2[1:iter, , drop = FALSE],
          params_m1_last = if (emm) params_m1_last else NULL,
          active             = active_list,
          init   = init,
          iter   = iter
        ))
      }
    }
    ll_prev  <- ll_new          
    obj_prev <- obj_new    
    }
      
  
  
  logf(1, "Reached max_iter without meeting tol.\n")  
  pi_ik_final <- if (use_covariates) {
    update_pi_ik(X_class, params$gate)
  } else {
    matrix(rep(params$pi, each = n), nrow = n, ncol = K)
  }
  
  Tau_final <- update_tau(
    Y = Y, Xi = Xi, Zi = Zi, params = params,
    pi_ik = pi_ik_final,
    C = C, cov_pred = cov_pred,
    include_const = TRUE
  )$Tau
  
  if(use_covariates){pi_hat_tau <- colMeans(Tau_final)} else {pi_hat_tau <- NULL}
  
  # trim histories (full length if no convergence)
  loglik_hist   <- loglik_hist[1:max_iter]
  lcm_iter_hist <- lcm_iter_hist[1:max_iter, , drop = FALSE]
  if (penalty) {
    penalty_total_hist <- penalty_total_hist[1:max_iter] 
    loglik_pen_hist    <- loglik_pen_hist[1:max_iter]   
  }
  
  list(
    covariate_scaling = lcm_design$covariate_scaling,
    parameter_scale = "standardized", active_scale = "standardized",
    params = params,
    Tau    = Tau_final,
    pi_hat_tau = pi_hat_tau,
    loglik_hist   = loglik_hist,
    lcm_iter_hist = lcm_iter_hist,
    penalty_total_hist = if (penalty&& !emm) penalty_total_hist else NULL,   
    loglik_pen_hist    = if (penalty&& !emm) loglik_pen_hist else NULL, 
    emm = emm,
    emm_beta_only = emm_beta_only,
    lcm_iter_hist_m1 = lcm_iter_hist_m1[1:iter, , drop = FALSE],
    lcm_iter_hist_m2 = lcm_iter_hist_m2[1:iter, , drop = FALSE],
    params_m1_last = if (emm) params_m1_last else NULL,
    active  = active_list,
    init   = init,
    iter   = max_iter
  )
}
