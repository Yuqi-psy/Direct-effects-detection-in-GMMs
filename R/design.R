prep_lcm_design <- function (
    y, 
    data, 
    growth_covariates = NULL, 
    cov_pred = c("none", "intercept", "slope", "both"), 
    penalty = FALSE, 
    class_covariates = NULL
    ){
  
    cov_pred <- match.arg(cov_pred)
    
    cc <- if (is.null(class_covariates))
      character()
    else class_covariates
    
    gc <- if (is.null(growth_covariates) || cov_pred == "none") 
      character()
    else growth_covariates
    covars <- unique(c(cc, gc))
    data_std <- data
    
    # standardizes covariates
    center_vec <- scale_vec <- setNames(numeric(), character())
    if (length(covars)) {
      if (!all(covars %in% names(data)) || !all(vapply(data[, covars, drop = FALSE], is.numeric, logical(1)))) 
        stop("Covariates must be existing numeric columns.")
      C_raw <- as.matrix(data[, covars, drop = FALSE])
      if (any(!is.finite(C_raw))) 
        stop("Covariates must contain finite values.")
      center_vec <- colMeans(C_raw)
      scale_vec <- apply(C_raw, 2, sd)
      if (any(!is.finite(scale_vec) | scale_vec <= 0)) 
        stop("Each covariate must have a positive finite SD.")
      data_std[, covars] <- sweep(sweep(C_raw, 2, center_vec, "-"), 2, scale_vec, "/")
    }
    
    
    # prepare design matrix
    n <- nrow(data)
    Tobs <- length(y)
    y_vec <- matrix(t(as.matrix(data[, y, drop = FALSE])), nrow = n * Tobs, ncol = 1)
    t_long <- rep(0:(Tobs - 1), times = n)
    X_0 <- cbind(1, time = t_long)
    X_1 <- matrix(0, nrow = n * Tobs, ncol = 0)
    C <- if (length(gc)) 
      as.matrix(data_std[, gc, drop = FALSE])
    else NULL
    if (length(gc)) {
      C_long <- C[rep(seq_len(n), each = Tobs), , drop = FALSE]
      if (cov_pred == "intercept") 
        X_1 <- C_long
      if (cov_pred == "slope") 
        X_1 <- C_long * t_long
      if (cov_pred == "both") 
        X_1 <- cbind(C_long, C_long * t_long)
    }
    X_class <- if (length(cc)) 
      as.matrix(data_std[, cc, drop = FALSE])
    else NULL
    
    
    scaling <- list(center = center_vec, scale = scale_vec, class_covariates = cc, growth_covariates = gc, 
                    cov_pred = cov_pred)
    
    
    list(n = n, Tobs = Tobs, y_vec = y_vec, t_long = t_long, data_std = data_std, C = C, X_class = X_class, 
         x = cbind(X_0, X_1), z_i = cbind(1, 0:(Tobs - 1)), X_0 = X_0, X_1 = X_1, X_1_pen = X_1, center_vec = center_vec, 
         scale_vec = scale_vec, covariate_scaling = scaling, 
         meta = list(y_names = y, growth_covariates = if (is.null(growth_covariates)) character() else growth_covariates, 
                                                                         cov_pred = cov_pred, penalty_prepped = penalty, covariate_scale = "standardized"))
  }


gmm_convert_params <- function (
    params, 
    covariate_scaling, # extract from 'prep_lcm_design( )' output
    from = c("standardized", "raw"), 
    to = c("raw", "standardized")
    ){
  
    from <- match.arg(from)
    to <- match.arg(to)
    
    sc <- covariate_scaling
    gc <- sc$growth_covariates
    cc <- sc$class_covariates
    mu <- sc$center[gc]
    ss <- sc$scale[gc]
    q <- length(gc)
    p_de <- if (sc$cov_pred == "both")   2 * q
    else if (sc$cov_pred == "none")  0
    else q
    ii <- if (sc$cov_pred %in% c("intercept", "both")) 2 + seq_len(q)
    else integer()
    si <- if (sc$cov_pred == "both") 2 + q + seq_len(q)
    else if (sc$cov_pred == "slope") 2 + seq_len(q)
    else integer()
    
    # convert within class parameters
    params$beta <- lapply(params$beta, function(b) {
      b <- as.numeric(b)
      if (length(b) == 2 && p_de > 0) 
        b <- c(b, rep(0, p_de))
      if (length(b) != 2 + p_de) 
        stop("beta length must match the growth design.")
      for (base_idx in 1:2) {
        idx <- list(ii, si)[[base_idx]]
        if (!length(idx) || from == to) 
          next
        if (from == "raw") {
          b[base_idx] <- b[base_idx] + sum(b[idx] * mu)
          b[idx] <- b[idx] * ss
        }
        else {
          b[idx] <- b[idx]/ss
          b[base_idx] <- b[base_idx] - sum(b[idx] * mu)
        }
      }
      b
    })
    
    # convert gate parameters
    if (!is.null(params$gate)) {
      B <- as.matrix(params$gate$Beta)
      aa <- as.numeric(params$gate$alpha)
      if (ncol(B) != length(cc) || nrow(B) != length(aa)) 
        stop("gate dimensions must match class_covariates.")
      if (from != to) {
        if (from == "raw") {
          aa <- aa + as.numeric(B %*% sc$center[cc])
          B <- sweep(B, 2, sc$scale[cc], "*")
        }
        else {
          B <- sweep(B, 2, sc$scale[cc], "/")
          aa <- aa - as.numeric(B %*% sc$center[cc])
        }
      }
      colnames(B) <- cc
      params$gate$alpha <- aa
      params$gate$Beta <- B
    }
    params
  }