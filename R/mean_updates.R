
# Beta: full unpenalized update.
update_beta <- function(
    x,
    omega,
    y,
    n,
    weights = rep(1, n), # length n vector of subject weights
    Debug = FALSE
) {
  if (Debug) browser()
  
  # Compute t(x) %*% omega.w using the same pseudo-inverse and subject blocks.
  P <- pracma::pinv(omega)
  Tobs <- nrow(omega)
  w_long <- rep(weights, each = Tobs)
  
  x_blocks <- matrix(x, nrow = Tobs, ncol = n * ncol(x))
  Z <- matrix(t(P) %*% x_blocks, nrow = n * Tobs, ncol = ncol(x))
  A <- t(Z * w_long)
  rownames(A) <- colnames(x)
  
  # Preserve the original solve/inverse expression in this experiment.
  beta <- solve(A %*% x) %*% (A %*% y)
  return(beta)
}



update_beta_penalty <- function(
    X_0,
    X_1, 
    omega,       # t*t
    y,           # nt*1
    n,
    beta_ini = NULL,
    weights = rep(1, n),
    lambda,
    trace = FALSE,
    Debug = FALSE
){
  
  if (Debug) browser()
  eps <- 1e-10
  
  # -------------------------------------------------
  # outer penalty: used for reporting Q_pen = Q - P
  # inner penalty: used for penalized beta optimization
  # because max(Q - P) is min(q + 2P)
  # -------------------------------------------------
  lambda_outer <- as.numeric(lambda)
  lambda_inner <- 2 * lambda_outer
  
  # prepare
  X <- cbind(X_0, X_1)
  p0 <- ncol(X_0); p1 <- ncol(X_1)
  if (!p1) stop("No candidate DE columns were supplied.")
  
  b <- as.matrix(beta_ini)
  if (nrow(b) != p0 + p1 || ncol(b) != 1L) stop("beta_ini dimensions do not match the design.")
  beta_0 <- b[seq_len(p0), , drop = FALSE]
  beta_1 <- b[p0 + seq_len(p1), , drop = FALSE]
  
 # omega.w <- kronecker(diag(weights, n, n), pracma::pinv(omega))
 # A_0 <- t(X_0) %*% omega.w
 # A_1 <- t(X_1_pen) %*% omega.w
  
  P <- pracma::pinv(omega)
  Tobs <- nrow(omega)
  w_long <- rep(weights, each = Tobs)
  make_A <- function(X_part) {
    X_blocks <- matrix(X_part, nrow = Tobs, ncol = n * ncol(X_part))
    Z <- matrix(t(P) %*% X_blocks, nrow = n * Tobs, ncol = ncol(X_part))
    A <- t(Z * w_long)
    rownames(A) <- colnames(X_part)
    A
  }
  quad_loss <- function(r) {
    E <- matrix(r, nrow = Tobs, ncol = n, byrow = FALSE)
    sum(weights * colSums(E * (P %*% E)))
  }
  A_0 <- make_A(X_0); A_1 <- make_A(X_1); Ip <- diag(p1)
  penalty <- lambda_inner[1L] * sum(abs(beta_1)) + lambda_inner[2L] * sum(beta_1^2)
  loss.old <- quad_loss(y - X %*% b) + penalty
  # initial inner loss = q + 2P

  
  G00_inv <- solve(A_0 %*% X_0)
  G11 <- A_1 %*% X_1
  maxiter <- 1000L
  for (iter in seq_len(maxiter)) {
    beta_0 <- G00_inv %*% (A_0 %*% (y - X_1 %*% beta_1))
    y_1 <- y - X_0 %*% beta_0
    bold <- beta_1
    
    # Explicit nrow preserves the intended diagonal when p1 = 1.
    # Exact-zero locking is retained from the uploaded MM implementation.
   D1 <- if (lambda_inner[1L] == 0) matrix(0, p1, p1) else
     lambda_inner[1L]/2 * diag(1/abs(drop(bold)), nrow = p1)
    # D1 <- if (lambda_inner[1L] == 0) {
    #   matrix(0, p1, p1)
    # } else {
    #   lambda_inner[1L] / 2 *
    #     diag(1 / pmax(abs(drop(bold)), eps), nrow = p1)
    # }
    # 
    D2 <- lambda_inner[2L] * Ip
    Dtot <- D1 + D2
    
    dd <- diag(Dtot); idx <- which(dd > 1/eps)
    
    if (length(idx)) { dd[idx] <- 1/eps; diag(Dtot) <- dd }
    beta_1 <- solve(G11 + Dtot, A_1 %*% y_1)
    
    if (length(idx)) beta_1[idx] <- bold[idx]
    penalty <- lambda_inner[1L] * sum(abs(beta_1)) + lambda_inner[2L] * sum(beta_1^2)
    b <- rbind(beta_0, beta_1)
    loss.new <- quad_loss(y - X %*% b) + penalty
    delta_loss <- loss.new - loss.old
    stop_rule <- 2 * (loss.old - loss.new)/(abs(loss.old) + abs(loss.new))
    
    if (trace) cat(sprintf("[beta iter %03d] loss.old=%.10f loss.new=%.10f delta=%+.3e relchg=%+.3e\n",
                           iter, loss.old, loss.new, delta_loss, stop_rule))
    if (delta_loss > 1e-8) cat(sprintf("[beta iter %03d] WARNING: beta loss increased by %.3e\n", iter, delta_loss))
    if (delta_loss <= 0 && is.finite(stop_rule) && abs(stop_rule) < 1e-6) break
    if (iter == maxiter) warning("Maximum number of iterations reached")
    loss.old <- loss.new
  }
  list(y = y, X = X, beta = b, beta_std = b, loss = loss.new,
       penalty = list(pen1 = lambda_outer[1L] * sum(abs(beta_1)),
                      pen2 = lambda_outer[2L] * sum(beta_1^2)),
       parameter_scale = "standardized")
}


update_beta_active <- function(
    X_0,
    X_1,          # raw design matrix for penalized fixed effects
    omega,        # T x T
    y,            # nT x 1
    n,
    active,       # when refit, for estimated parameters (1), not estimate parameters (0)
    weights = rep(1, n),
    Debug = FALSE
) {
  # IMPORTANT:
  # active must correspond to the column order of X_1,
  # and beta_raw is returned as c(beta_0, beta_1)
  
  if (Debug) browser()
  active <- as.logical(active)
  p0 <- ncol(X_0)
  p1 <- ncol(X_1)
  stopifnot(length(active) == p1)
  
  # Weighted precision: apply blocks without constructing the nT x nT matrix.
  P <- pracma::pinv(omega)
  Tobs <- nrow(omega)
  w_long <- rep(weights, each = Tobs)
  
  # Return the same quantity as t(X_part) %*% omega.w.
  make_A <- function(X_part) {
    X_blocks <- matrix(X_part, nrow = Tobs, ncol = n * ncol(X_part))
    Z <- matrix(t(P) %*% X_blocks,
                nrow = n * Tobs, ncol = ncol(X_part))
    A <- t(Z * w_long)
    rownames(A) <- colnames(X_part)
    A
  }
  
  # Initialize full beta_1 with zeros.
  beta_1 <- matrix(0, nrow = p1, ncol = 1)
  
  if (!any(active)) {
    # Only unpenalized block remains.
    A0 <- make_A(X_0)
    beta_0 <- solve(A0 %*% X_0, A0 %*% y)
  } else {
    X_1_act <- X_1[, active, drop = FALSE]
    X_refit <- cbind(X_0, X_1_act)
    
    A <- make_A(X_refit)
    beta_hat <- solve(A %*% X_refit, A %*% y)
    
    beta_0 <- beta_hat[1:p0, , drop = FALSE]
    beta_1[active, 1] <- beta_hat[(p0 + 1):nrow(beta_hat), 1]
  }
  
  beta_std <- rbind(beta_0, beta_1)
  
  X_full <- cbind(X_0, X_1)
  r <- y - X_full %*% beta_std
  E <- matrix(r, nrow = Tobs, ncol = n, byrow = FALSE)
  loss <- as.numeric(sum(weights * colSums(E * (P %*% E))))
  output <- list(y = y, X = X_full, beta_0 = beta_0, beta_1 = beta_1, 
                 beta_std = beta_std, active = active, loss = loss)
  return(output)
}

















