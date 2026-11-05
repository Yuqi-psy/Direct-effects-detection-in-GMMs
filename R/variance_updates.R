# Omega

update_omega <- function(
    z,
    theta,
    sigma2,
    n
){
  #browser()
  
  omega <- z %*% theta %*% t(z) + as.numeric(sigma2) * diag(nrow(z))
  
  return(omega)
  
}


# Sigma2 (weighted)
update_sigma2 <- function(
    y,
    x,
    beta,
    omega,
    sigma2,
    n,
    weights =  rep(1, n),
    sigma_floor = 1e-8,
    Debug = FALSE
){
  if(Debug) browser()
  
  # preparing
  #omega_inv <- solve(omega)
  omega_inv <- pracma::pinv(omega) # general methods take inverse
  omega_inv2 <- omega_inv %*% omega_inv
  res <- y - x %*% beta
  design_res <- matrix(res, nrow = nrow(omega), ncol = n)
  
  #denom <- sum(weights*(diag(t(design_res) %*% omega_inv2 %*% design_res)))
  denom <- sum(
    weights * colSums(design_res * (omega_inv2 %*% design_res))
  )
  numer <- sum(weights)* sum(diag(omega_inv))
  
  ratio <- as.numeric(denom / numer)
  if (!is.finite(ratio) || ratio <= 0) ratio <- 1
  
  sigma2_new <- as.numeric(sigma2) * sqrt(ratio)
  sigma2_new <- max(as.numeric(sigma2_new), sigma_floor)
  
  return(sigma2_new)
}


# update theta for real data probelm
update_theta2 <- function(
    sigma2,
    omega,
    y,
    x,
    z,
    beta,
    theta,
    n,
    weights = rep(1, n),  
    
    freeze_A = FALSE, 
    freeze_B = FALSE, 
    A_fixed = NULL, 
    B_fixed = NULL,
    
    theta_floor = 1e-8,
    Debug = FALSE
){
  
  if(Debug) browser()
  
  is.positive.definite <- function(mat) {
    mat <- Re((mat + t(mat)) / 2)
    eigenvalues <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
    all(eigenvalues > theta_floor)
  } 
  
  omega <- Re((omega + t(omega)) / 2) # force symmetry
  omega_inv <- pracma::pinv(omega) # persdo inverse
  
  res <- y - x %*% beta
  design_res <- matrix(res, nrow = nrow(omega), ncol = n, byrow = FALSE)
  
  R <- theta %*% t(z) %*% omega_inv %*% design_res
  
  if (freeze_A) {
    A <- A_fixed
  } else {
    # A = R diag(weights) R' is PSD for nonnegative weights.
    A <- tcrossprod(sweep(R, 2L, weights, FUN = "*"), R)
    A <- Re((A + t(A)) / 2)
    
    eigA <- base::eigen(A, symmetric = TRUE)
    roundoff_A <- 64 * .Machine$double.eps *
      max(1, max(abs(eigA$values)))
    
    min_A <- min(eigA$values)
    
    # Fall back to nearPD only for a materially negative eigenvalue.
    if (min_A < -roundoff_A) {
      warning(sprintf(
        "A has minimum eigenvalue %.3e below -%.3e; applying nearPD.",
        min_A, roundoff_A
      ), call. = FALSE)
      A <- as.matrix(Matrix::nearPD(A)$mat)
      A <- Re((A + t(A)) / 2)
      
      # Clip roundoff-sized negatives only; retain ALL positive eigenvalues.
    } else if (min_A < 0) {
      A <- tcrossprod(sweep(
        eigA$vectors, 2L, sqrt(pmax(eigA$values, 0)), FUN = "*"
      ))
      A <- Re((A + t(A)) / 2)
    }
  }
  
  if (freeze_B) {
    B <- B_fixed
  } else {
    B <- sum(weights) * (t(z) %*% omega_inv %*% z)
    B <- Re((B + t(B)) / 2) # protection
    if (!is.positive.definite(B)) {
      B <- B + 1e-8 * diag(nrow(B))
      B <- Re((B + t(B)) / 2)
      if (!is.positive.definite(B)) {
        B <- as.matrix(Matrix::nearPD(B)$mat)
      }
    }
  }
  
  if (!freeze_A || !freeze_B) {
    L <- t(Matrix::chol(B))
    L_inv <- pracma::pinv(L)
    
    C <- t(L) %*% A %*% L
    C <- Re((C + t(C)) / 2)
    
    eigC <- eigen(C, symmetric = TRUE)
    eigC$values[eigC$values < theta_floor] <- theta_floor
    D <- eigC$vectors %*% diag(sqrt(eigC$values), length(eigC$values)) %*% t(eigC$vectors)
    D <- Re((D + t(D)) / 2)
    
    theta_new <- t(L_inv) %*% D %*% L_inv
    theta_new <- Re((theta_new + t(theta_new)) / 2)
    
    eigT <- eigen(theta_new, symmetric = TRUE)
    eigT$values[eigT$values < theta_floor] <- theta_floor
    theta_new <- eigT$vectors %*% diag(eigT$values, length(eigT$values)) %*% t(eigT$vectors)
    theta_new <- Re((theta_new + t(theta_new)) / 2)
    
  } else {
    theta_new <- theta
  }
  
  return(list(theta = theta_new, A = A, B = B))
}