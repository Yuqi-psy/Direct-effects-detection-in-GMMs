source("load_functions.R")

# Run from the repository root after loading the original GMM functions.
#project_dir <- getwd()
project_dir <- "/Users/liuyuqi/Library/CloudStorage/OneDrive-UniversiteitLeiden/Project_2/P2B/project2B/Project2B"

# Input files and LG starting values.
c19prc_data_path <- file.path(
  project_dir, "real_data", "C19PRC", "LG_syntax", "C19PRC_GAD_no_allzero.csv"
)
charls_data_path <- file.path(project_dir, "real_data", "CHARLS", "CHARLS_cog.csv")
c19prc_lg_dir <- dirname(c19prc_data_path)
charls_lg_dir <- file.path(project_dir, "real_data", "CHARLS", "LG_syntax")
sim_lg_dir <- file.path(project_dir, "sim_data")

# 1. C19PRC
C19PRC_nonzero <- read.csv(c19prc_data_path)

ID <- "ID"
y <- paste0("y.", 1:5)
covars <- c("age", "loneliness", "death_anxiety", "oci")
design <- prep_lcm_design(
  y = y, data = C19PRC_nonzero, class_covariates = covars, growth_covariates = covars,
  cov_pred = "both", penalty = TRUE
)

# Preliminary R fit.
Raw_fit <- fit_gmm(
  data = C19PRC_nonzero, y = y, ID = ID, K = 2, class_covariates = covars,
  growth_covariates = covars, cov_pred = "both", penalty = FALSE, init_params = NULL,
  max_iter = 1000, tol = 1e-6, nstart = 5, seed = 123, verbose = TRUE, lcm_design = design
)

# LG starting values on the raw scale.
LG_C19PRC_2C <- LG_fit(
  param_path = file.path(c19prc_lg_dir, "C19PRC_pars_2c.txt"),
  tau_path = file.path(c19prc_lg_dir, "C19PRC_classification_2c.csv"), K = 2, p_class = 4,
  p_growth = 4, ref_ext = 1
)
LG_std <- LG_C19PRC_2C
LG_std$params <- gmm_convert_params(
  LG_C19PRC_2C$params, design$covariate_scaling, from = "raw", to = "standardized"
)

PSR <- tune_psr(
  data = C19PRC_nonzero, y = y, ID = ID, K = 2, nlambda = 40, class_covariates = covars,
  growth_covariates = covars, cov_pred = "both", prefit = LG_std, lcm_design = design
)

EMM <- tune_emm(
  data = C19PRC_nonzero, y = y, ID = ID, K = 2, nlambda = 40, class_covariates = covars,
  growth_covariates = covars, cov_pred = "both", prefit = LG_std, pre_lambdas = PSR$lambda_seq,
  lcm_design = design
)

# Refine the C19PRC lambda grid.
lambda_check <- sort(unique(c(
  PSR$lambda_seq,
  seq(23.8059401, 60.7962123, by = 0.5),
  23.8059401, 60.7962123
)))

PSR_lambdacheck <- tune_psr(
  data = C19PRC_nonzero, y = y, ID = ID, K = 2, nlambda = 40, class_covariates = covars,
  growth_covariates = covars, cov_pred = "both", prefit = LG_std, pre_lambdas = lambda_check,
  lcm_design = design
)

EMM_lambdacheck <- tune_emm(
  data = C19PRC_nonzero, y = y, ID = ID, K = 2, nlambda = 40, class_covariates = covars,
  growth_covariates = covars, cov_pred = "both", prefit = LG_std, pre_lambdas = lambda_check,
  lcm_design = design
)

# 2. CHARLS
gmm_setup()

CHARLS_cog <- read.csv(charls_data_path, check.names = FALSE)
ID_charls <- "id"
y_charls <- paste0("y.", 1:4)
covs_charls <- c("x1", "x4", "x5")

design_charls <- prep_lcm_design(
  y = y_charls, data = CHARLS_cog, class_covariates = covs_charls,
  growth_covariates = covs_charls, cov_pred = "both", penalty = TRUE
)

LG_charls <- LG_fit(
  param_path = file.path(charls_lg_dir, "CHARLS_K2_initial_parameters.txt"),
  tau_path = file.path(charls_lg_dir, "CHARLS_K2_initial_classification.csv"), K = 2, p_class = 3,
  p_growth = 3, ref_ext = 1
)

initial_classification <- read.csv(
  file.path(charls_lg_dir, "CHARLS_K2_initial_classification.csv"), check.names = FALSE
)
i_initial <- match(CHARLS_cog[[ID_charls]], initial_classification[[1]])
LG_charls$Tau <- LG_charls$Tau[i_initial, , drop = FALSE]
LG_charls$ID <- CHARLS_cog[[ID_charls]]

# LG order: gate, two means per factor, then class-specific variances.
base_values <- scan(file.path(charls_lg_dir, "CHARLS_K2_base_parameters.txt"),
                    what = numeric(), quiet = TRUE)

base_classification <- read.csv(
  file.path(charls_lg_dir, "CHARLS_K2_base_classification.csv"), check.names = FALSE
)
base_tau <- as.matrix(base_classification[, c("Class#1", "Class#2")])
LG_charls_base <- list(
  params = list(
    K = 2, pi = colMeans(base_tau),
    beta = list(c(base_values[5], base_values[7]),
                c(base_values[6], base_values[8])),
    Psi = list(matrix(base_values[c(10, 11, 11, 12)], 2, 2),
               matrix(base_values[c(14, 15, 15, 16)], 2, 2)),
    sigma2 = base_values[c(9, 13)],
    gate = list(alpha = base_values[1],
                Beta = matrix(base_values[2:4], 1, 3,
                              dimnames = list(NULL, covs_charls)))
  ),
  Tau = base_tau, ID = base_classification[[1]]
)

LG_charls_std <- LG_charls
LG_charls_std$params <- gmm_convert_params(
  LG_charls$params, design_charls$covariate_scaling, from = "raw", to = "standardized"
)
LG_charls_base_std <- LG_charls_base
LG_charls_base_std$params <- gmm_convert_params(
  LG_charls_base$params, design_charls$covariate_scaling, from = "raw", to = "standardized"
)

prefit_charls <- LG_charls_std
LG_full_min_eigen <- numeric(2)
for (k in 1:2) {
  LG_full_min_eigen[k] <- min(eigen(
    LG_charls_std$params$Psi[[k]], symmetric = TRUE, only.values = TRUE
  )$values)
}

# R prefit when the LG covariance is inadmissible.
if (any(LG_full_min_eigen <= 0) || any(LG_charls_std$params$sigma2 <= 0)) {
  init_charls <- LG_charls_base_std$params
  prefit_charls <- fit_gmm(
    data = CHARLS_cog, y = y_charls, ID = ID_charls, K = 2, class_covariates = covs_charls,
    growth_covariates = covs_charls, cov_pred = "both", penalty = FALSE, active = NULL,
    init_params = init_charls, lcm_design = design_charls, max_iter = 1000, tol = 1e-6,
    nstart = 1, seed = 123, verbose = TRUE
  )
}

PSR_charls <- tune_psr(
  data = CHARLS_cog, y = y_charls, ID = ID_charls, K = 2, class_covariates = covs_charls,
  growth_covariates = covs_charls, cov_pred = "both", prefit = prefit_charls, nlambda = 40,
  lambda_min_ratio = 1e-2, lambda_expand = 1.2, lcm_design = design_charls
)

EMM_charls <- tune_emm(
  data = CHARLS_cog, y = y_charls, ID = ID_charls, K = 2, class_covariates = covs_charls,
  growth_covariates = covs_charls, cov_pred = "both", prefit = prefit_charls,
  pre_lambdas = PSR_charls$lambda_seq, lambda_min_ratio = 1e-2, lambda_expand = 1.2,
  lcm_design = design_charls, Debug = FALSE
)

# Refine the CHARLS lambda grid.
lambda_check <- sort(unique(c(
  PSR_charls$lambda_seq,
  seq(54.548518, 164.662237, by = 0.5),
  54.548518, 164.662237
)))

PSR_charls_lambda_check <- tune_psr(
  data = CHARLS_cog, y = y_charls, ID = ID_charls, K = 2, class_covariates = covs_charls,
  growth_covariates = covs_charls, cov_pred = "both", prefit = prefit_charls,
  pre_lambdas = lambda_check, lcm_design = design_charls, seed = 123
)

EMM_charls_lambda_check <- tune_emm(
  data = CHARLS_cog, y = y_charls, ID = ID_charls, K = 2, class_covariates = covs_charls,
  growth_covariates = covs_charls, cov_pred = "both", prefit = prefit_charls,
  pre_lambdas = lambda_check, lcm_design = design_charls, Debug = FALSE
)

# 3. Simulated data
K <- 3

beta_list <- list(
  c(0.0,  0.50),
  c(3.0, -0.50),
  c(6.0,  0.20)
)

theta_list <- list(
  matrix(c(0.40, 0.10,
           0.10, 0.10), 2, 2, byrow = TRUE),
  matrix(c(0.35, 0.11,
           0.11, 0.15), 2, 2, byrow = TRUE),
  matrix(c(0.45, 0.10,
           0.10, 0.20), 2, 2, byrow = TRUE)
)

sigma2_vec <- c(0.50, 0.55, 0.40)

alpha <- c(-0.4, 0.3)

Beta <- rbind(
  c( 0.6, -0.3,  0.5),
  c(-0.2,  0.4, -0.3)
)

Gamma_list <- list(
  matrix(c(
    0.5, 0.0, -0.3,
    0.0, 0.2, 0.2
  ), nrow = 2, byrow = TRUE),

  matrix(c(
    0.3, 0.0, 0.3,
    0.2, 0.0, 0.0
  ), nrow = 2, byrow = TRUE),

  matrix(c(
    -0.3, 0.0, 0.2,
    0.0, 0.0, -0.3
  ), nrow = 2, byrow = TRUE)
)

sim1 <- sim_gmm_cov(
  n = 400, time = c(0, 1, 2, 3, 4, 5), K = K, beta_list = beta_list, theta_list = theta_list,
  sigma2_vec = sigma2_vec, x_type = c("cont", "cont", "cont"), x_mean = c(0, 0, 0),
  x_sd = c(1, 1, 1), x_prob = c(0.5, 0.5, 0.4), alpha = alpha, Beta = Beta,
  Gamma_list = Gamma_list, seed = 123
)

design_sim1 <- prep_lcm_design(
  y = paste0("y.", 0:5), data = sim1$data, class_covariates = paste0("x", 1:3),
  growth_covariates = paste0("x", 1:3), cov_pred = "both", penalty = TRUE
)

# Unpenalized baseline and full-DE R fits.
R_noDE <- fit_gmm(
  y = paste0("y.", 0:5), ID = "id", data = sim1$data, K = 3, max_iter = 1000, tol = 1e-6,
  nstart = 5, seed = NULL, verbose = TRUE, Debug = FALSE, class_covariates = paste0("x", 1:3),
  growth_covariates = paste0("x", 1:3), cov_pred = "both", free_mode = "pi_beta_sigma2_theta",
  penalty = FALSE, init_params = NULL, active = replicate(3, rep(0L, 6), simplify = FALSE),
  Lambda = NULL
)

R_fullDE <- fit_gmm(
  y = paste0("y.", 0:5), ID = "id", data = sim1$data, K = 3, max_iter = 1000, tol = 1e-6,
  nstart = 5, seed = NULL, verbose = TRUE, Debug = FALSE, class_covariates = paste0("x", 1:3),
  growth_covariates = paste0("x", 1:3), cov_pred = "both", free_mode = "pi_beta_sigma2_theta",
  penalty = FALSE, init_params = NULL, active = NULL, Lambda = NULL
)

LG_sim1_3k <- LG_fit(
  param_path = file.path(sim_lg_dir, "sim1_K3_initial_parameters.txt"),
  tau_path = file.path(sim_lg_dir, "sim1_k3_initial_classification.csv"), K = 3, p_class = 3,
  p_growth = 3, ref_ext = 1
)
LG_sim1_std <- LG_sim1_3k
LG_sim1_std$params <- gmm_convert_params(
  LG_sim1_3k$params, design_sim1$covariate_scaling, from = "raw", to = "standardized"
)

PSR_sim1 <- tune_psr(
  data = sim1$data, y = paste0("y.", 0:5), ID = "id", K = 3, nlambda = 40,
  class_covariates = paste0("x", 1:3), growth_covariates = paste0("x", 1:3), cov_pred = "both",
  prefit = LG_sim1_std, lcm_design = design_sim1
)

EMM_sim1 <- tune_emm(
  data = sim1$data, y = paste0("y.", 0:5), ID = "id", K = 3, nlambda = 40,
  class_covariates = paste0("x", 1:3), growth_covariates = paste0("x", 1:3), cov_pred = "both",
  prefit = LG_sim1_std, pre_lambdas = PSR_sim1$lambda_seq, lcm_design = design_sim1
)

# Refine the simulated-data lambda grid.
lambda_check3 <- sort(unique(c(
  PSR_sim1$lambda_seq,
  seq(37.3613567, 115.0938354, by = 0.5),
  37.3613567, 115.0938354
)))
PSR_sim1_refine <- tune_psr(
  data = sim1$data, y = paste0("y.", 0:5), ID = "id", K = 3, pre_lambdas = lambda_check3,
  class_covariates = paste0("x", 1:3), growth_covariates = paste0("x", 1:3), cov_pred = "both",
  prefit = LG_sim1_std, lcm_design = design_sim1
)
EMM_sim1_refine <- tune_emm(
  data = sim1$data, y = paste0("y.", 0:5), ID = "id", K = 3, pre_lambdas = lambda_check3,
  class_covariates = paste0("x", 1:3), growth_covariates = paste0("x", 1:3), cov_pred = "both",
  prefit = LG_sim1_std, lcm_design = design_sim1
)
