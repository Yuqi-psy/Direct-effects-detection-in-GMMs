# this file use for set up the environment for estimating GMM

gmm_setup <- local({
  initialized <- FALSE
  
  function() {
    if (initialized) return(invisible(NULL))
    
    functions_path <- c(
      "R/design.R",
      "R/initialization.R",
      "R/likelihood.R",
      "R/class_membership.R",
      "R/mean_updates.R",
      "R/variance_updates.R",
      "R/fit_lcm.R",
      "R/fit_gmm.R",
      "R/fit_emm.R",
      "R/tuning.R",
      "R/simulation.R"
    )
    invisible(lapply(functions_path, source))
    
    packages <- c("tidyr", "dplyr", "knitr", "Matrix", "MASS", "expm", "pracma")
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is required but not installed.", pkg))
      }
      # pracma is used via ::; keep expm::sqrtm as the attached sqrtm.
      if (pkg != "pracma" && !paste0("package:", pkg) %in% search()) {
        suppressPackageStartupMessages(
          library(pkg, character.only = TRUE)
        )
      }
    }
    
    initialized <<- TRUE
    invisible(NULL)
  }
})
