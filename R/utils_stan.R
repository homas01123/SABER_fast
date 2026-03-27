#' Compare BayesianTools vs Stan/HMC Performance
#'
#' Runs both inversion methods on the same data and compares:
#' - Runtime
#' - Convergence (Rhat, ESS)
#' - Parameter estimates
#' - Posterior uncertainty
#'
#' @param rrs Data frame with wavelength and rrs_0m
#' @param par_inversed Vector of parameter names
#' @param par_fixed List of fixed parameters
#' @param bt_iterations BayesianTools iterations (default = 30000)
#' @param stan_iterations Stan post-warmup iterations (default = 2000)
#' @param verbose Print detailed output?
#'
#' @return List with comparison results
#' @export
compare_mcmc_methods <- function(
    rrs,
    par_inversed,
    par_fixed = NULL,
    bt_iterations = 30000,
    stan_iterations = 2000,
    verbose = TRUE
) {
  
  if (verbose) {
    message("\n", strrep("=", 70))
    message("MCMC METHOD COMPARISON: BayesianTools vs Stan/HMC")
    message(strrep("=", 70), "\n")
  }
  
  # Prepare lower/upper bounds for BayesianTools
  lower <- sapply(par_inversed, function(p) {
    switch(p,
           "chl" = 0.1,
           "a_g_440" = 0.01,
           "a_nap_440" = 0.001,
           "bb_p_550" = 0.0001,
           "a_g_slope" = 0.010,
           "bb_p_gamma" = 0.2,
           "h_w" = 0.5,
           0.0  # default for benthic fractions
    )
  })
  
  upper <- sapply(par_inversed, function(p) {
    switch(p,
           "chl" = 50,
           "a_g_440" = 3,
           "a_nap_440" = 1,
           "bb_p_550" = 0.1,
           "a_g_slope" = 0.025,
           "bb_p_gamma" = 2.0,
           "h_w" = 50,
           1.0  # default for benthic fractions
    )
  })
  
  # ========================================================================
  # METHOD 1: BayesianTools (DEzs)
  # ========================================================================
  
  if (verbose) message("Running BayesianTools (DEzs)...\n")
  
  start_time_bt <- Sys.time()
  
  tryCatch({
    result_bt <- inverse_mcmc(
      rrs = rrs,
      forward_model = "am03",
      par_inversed = par_inversed,
      lower = lower,
      upper = upper,
      par_fixed = par_fixed,
      iterations = bt_iterations,
      burnin = bt_iterations / 3,
      sampler = "DEzs"
    )
    
    runtime_bt <- as.numeric(difftime(Sys.time(), start_time_bt, units = "secs"))
    
    if (verbose) {
      message(sprintf("  Runtime: %.1f seconds (%.1f minutes)\n", 
                      runtime_bt, runtime_bt / 60))
    }
    
    bt_success <- TRUE
    
  }, error = function(e) {
    message("  ERROR: BayesianTools failed - ", conditionMessage(e), "\n")
    result_bt <- NULL
    runtime_bt <- NA
    bt_success <- FALSE
  })
  
  # ========================================================================
  # METHOD 2: Stan/HMC (NUTS)
  # ========================================================================
  
  if (verbose) message("Running Stan/HMC (NUTS)...\n")
  
  start_time_stan <- Sys.time()
  
  tryCatch({
    result_stan <- inverse_mcmc_stan(
      rrs = rrs,
      forward_model = "am03",
      par_inversed = par_inversed,
      par_fixed = par_fixed,
      iterations = stan_iterations,
      warmup = stan_iterations / 2,
      chains = 4,
      adapt_delta = 0.9,
      return_fit = TRUE,
      verbose = FALSE  # Suppress Stan output for comparison
    )
    
    runtime_stan <- as.numeric(difftime(Sys.time(), start_time_stan, units = "secs"))
    
    if (verbose) {
      message(sprintf("  Runtime: %.1f seconds (%.1f minutes)\n", 
                      runtime_stan, runtime_stan / 60))
    }
    
    stan_success <- TRUE
    
  }, error = function(e) {
    message("  ERROR: Stan/HMC failed - ", conditionMessage(e), "\n")
    result_stan <- NULL
    runtime_stan <- NA
    stan_success <- FALSE
  })
  
  # ========================================================================
  # COMPARISON
  # ========================================================================
  
  if (verbose) {
    message("\n", strrep("=", 70))
    message("RESULTS COMPARISON")
    message(strrep("=", 70), "\n")
    
    # Runtime comparison
    message("RUNTIME:")
    if (bt_success) {
      message(sprintf("  BayesianTools: %.1f sec (%.1f min)", 
                      runtime_bt, runtime_bt / 60))
    }
    if (stan_success) {
      message(sprintf("  Stan/HMC:      %.1f sec (%.1f min)", 
                      runtime_stan, runtime_stan / 60))
    }
    if (bt_success && stan_success) {
      speedup <- runtime_bt / runtime_stan
      message(sprintf("  Speedup:       %.1fx faster with Stan\n", speedup))
    }
    
    # Parameter estimates comparison
    message("PARAMETER ESTIMATES:")
    message(sprintf("  %-15s %15s %15s %15s", 
                    "Parameter", "BayesianTools", "Stan/HMC", "Difference"))
    message(strrep("-", 70))
    
    for (param in par_inversed) {
      bt_val <- if (bt_success) result_bt[[param]] else NA
      stan_val <- if (stan_success) result_stan$estimates[[param]] else NA
      diff <- if (!is.na(bt_val) && !is.na(stan_val)) {
        abs(bt_val - stan_val) / bt_val * 100
      } else NA
      
      message(sprintf("  %-15s %15.4f %15.4f %14.1f%%", 
                      param, bt_val, stan_val, diff))
    }
    
    # Convergence diagnostics (Stan only)
    if (stan_success) {
      message("\nSTAN CONVERGENCE DIAGNOSTICS:")
      diag <- result_stan$diagnostics
      message(sprintf("  Divergent transitions: %d", 
                      sum(diag$num_divergent)))
      message(sprintf("  Max treedepth exceeded: %d", 
                      sum(diag$num_max_treedepth)))
      
      # Check Rhat and ESS from summary
      summary <- result_stan$summary
      rhat_max <- max(summary$rhat, na.rm = TRUE)
      ess_min <- min(summary$ess_bulk, na.rm = TRUE)
      
      message(sprintf("  Max Rhat: %.4f (should be < 1.01)", rhat_max))
      message(sprintf("  Min ESS: %.0f (should be > 400)", ess_min))
      
      if (rhat_max > 1.01) {
        message("  WARNING: Rhat > 1.01, chains may not have converged!")
      }
      if (ess_min < 400) {
        message("  WARNING: ESS < 400, increase iterations!")
      }
    }
    
    message("\n", strrep("=", 70), "\n")
  }
  
  # Return structured comparison
  comparison <- list(
    bayesiantools = if (bt_success) {
      list(
        estimates = result_bt,
        runtime = runtime_bt,
        success = TRUE
      )
    } else {
      list(success = FALSE)
    },
    
    stan = if (stan_success) {
      list(
        estimates = result_stan$estimates,
        fit = result_stan$fit,
        diagnostics = result_stan$diagnostics,
        summary = result_stan$summary,
        runtime = runtime_stan,
        success = TRUE
      )
    } else {
      list(success = FALSE)
    },
    
    speedup = if (bt_success && stan_success) runtime_bt / runtime_stan else NA
  )
  
  return(invisible(comparison))
}


#' Validate Stan Model Compilation
#'
#' Checks if Stan is properly installed and can compile models.
#'
#' @return Logical, TRUE if Stan is ready
#' @export
validate_stan_installation <- function() {
  
  message("\n", strrep("=", 70))
  message("STAN INSTALLATION VALIDATION")
  message(strrep("=", 70), "\n")
  
  # Check cmdstanr package
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    message("✗ cmdstanr package not installed")
    message("\nInstall with:")
    message('  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))')
    return(FALSE)
  } else {
    message("✓ cmdstanr package installed")
  }
  
  # Check CmdStan
  stan_version <- tryCatch({
    cmdstanr::cmdstan_version()
  }, error = function(e) {
    NULL
  })
  
  if (is.null(stan_version)) {
    message("✗ CmdStan not found")
    message("\nInstall with:")
    message('  cmdstanr::install_cmdstan()')
    return(FALSE)
  } else {
    message(sprintf("✓ CmdStan version %s installed", stan_version))
  }
  
  # Check Stan model file
  stan_file <- system.file("stan", "am03_shallow.stan", package = "SABER")
  
  if (!file.exists(stan_file) || stan_file == "") {
    message("✗ Stan model file not found")
    message("  Expected: inst/stan/am03_shallow.stan")
    return(FALSE)
  } else {
    message("✓ Stan model file found")
  }
  
  # Try compiling model
  message("\nTesting model compilation...")
  
  compile_success <- tryCatch({
    mod <- cmdstanr::cmdstan_model(stan_file, quiet = TRUE)
    message("✓ Model compiled successfully!")
    TRUE
  }, error = function(e) {
    message("✗ Model compilation failed:")
    message("  ", conditionMessage(e))
    FALSE
  })
  
  if (compile_success) {
    message("\n", strrep("=", 70))
    message("STAN IS READY TO USE!")
    message(strrep("=", 70), "\n")
    return(TRUE)
  } else {
    message("\n", strrep("=", 70))
    message("STAN SETUP INCOMPLETE - See messages above")
    message(strrep("=", 70), "\n")
    return(FALSE)
  }
}
