#' Enhanced Inverse MCMC with Quantum Yield Retrieval and Optical Properties
#'
#' Extended version of inverse_mcmc that retrieves quantum yield (phi_f) along
#' with standard OACs, and returns additional optical properties (Ed, E0, K, PAR).
#'
#' @param rrs Data frame with wavelength and rrs_0m columns
#' @param par_inversed Vector of parameter names to inverse (must include "phi_f")
#' @param prior Prior function (see make_prior_bundle) or NULL for uniform
#' @param lower Numeric vector of lower bounds matching par_inversed
#' @param upper Numeric vector of upper bounds matching par_inversed
#' @param best Numeric vector of best guess values (optional)
#' @param par_fixed Named vector of fixed parameters
#' @param iterations Number of MCMC iterations (default = 30000)
#' @param burnin Number of burnin iterations (default = 5000)
#' @param sampler MCMC sampler ("DEzs", "DREAMzs", etc.)
#' @param lat Latitude for Ed calculation (default = 49)
#' @param lon Longitude for Ed calculation (default = -68)
#' @param date_time POSIXct datetime in UTC (default = Sys.time())
#' @param return_full_output Logical, return optical properties? (default = TRUE)
#' @param verbose Logical, print progress messages (default = TRUE)
#'
#' @return If return_full_output = FALSE: named vector of MAP estimates and uncertainties
#'         If return_full_output = TRUE: list with:
#'           - par_estimates: MAP estimates with uncertainties
#'           - rrs_modeled: Modeled Rrs using MAP estimates
#'           - optical_properties: List with Ed_0m, E0_0m, Kd, Ku_W, Ku_B, PAR
#'           - mcmc_output: Full BayesianTools output object
#'
#' @export
inverse_mcmc_sicf <- function(rrs,
                               par_inversed,
                               prior = NULL,
                               lower = NULL,
                               upper = NULL,
                               best = NULL,
                               par_fixed = NULL,
                               iterations = 30000,
                               burnin = 5000,
                               sampler = "DEzs",
                               lat = 49,
                               lon = -68,
                               date_time = Sys.time(),
                               return_full_output = TRUE,
                               verbose = TRUE) {
  
  # Validate that phi_f is in par_inversed
  if (!"phi_f" %in% par_inversed) {
    warning("'phi_f' not in par_inversed. Adding it automatically.")
    par_inversed <- c(par_inversed, "phi_f")
    lower <- c(lower, 0.005)
    upper <- c(upper, 0.03)
    if (!is.null(best)) best <- c(best, 0.02)
  }
  
  # Add lat, lon, date_time to par_fixed if not already present
  if (is.null(par_fixed)) par_fixed <- c()
  
  if (!"lat" %in% names(par_fixed)) par_fixed["lat"] <- lat
  if (!"lon" %in% names(par_fixed)) par_fixed["lon"] <- lon
  if (!"date_time" %in% names(par_fixed)) {
    par_fixed["date_time"] <- as.numeric(date_time)
  }
  
  # Create likelihood function using am03_sicf model
  likelihood <- objective_factory(
    model = "am03_sicf",
    objective = "log-ll",
    rrs_observed = rrs,
    par_inversed = par_inversed,
    par_fixed = par_fixed
  )
  
  # Create Bayesian setup
  setup <- BayesianTools::createBayesianSetup(
    prior = prior,
    likelihood = likelihood,
    lower = lower,
    best = best,
    upper = upper,
    names = par_inversed,
    parallel = FALSE
  )
  
  # Check setup
  BayesianTools::checkBayesianSetup(setup)
  
  # Run MCMC
  if (verbose) {
    message("\n========================================")
    message("Running MCMC with SICF model...")
    message("Parameters to retrieve: ", paste(par_inversed, collapse = ", "))
    message("Iterations: ", iterations, " | Burnin: ", burnin)
    message("========================================\n")
  }
  
  out <- BayesianTools::runMCMC(
    bayesianSetup = setup,
    settings = list(
      iterations = iterations,
      burnin = burnin,
      message = verbose
    ),
    sampler = sampler
  )
  
  # Calculate uncertainties (standard deviation of chains)
  estimates_sd <- purrr::map_df(
    .x = out[["chain"]],
    ~ apply(.x[, 1:(ncol(.x) - 3)], 2, sd)
  )
  estimates_sd <- colMeans(estimates_sd)
  
  # Get MAP estimates
  map_values <- BayesianTools::MAP(out)[[1]]
  
  # Combine estimates with uncertainties
  par_estimates <- stats::setNames(
    c(map_values, estimates_sd),
    c(names(map_values), paste0(names(map_values), "_sd"))
  )
  
  if (verbose) {
    message("\n========================================")
    message("MCMC Complete!")
    message("========================================")
    message("MAP Estimates:")
    for (i in seq_along(map_values)) {
      message(sprintf("  %s = %.6f ± %.6f", 
                      names(map_values)[i], 
                      map_values[i], 
                      estimates_sd[i]))
    }
    message("========================================\n")
  }
  
  # If return_full_output, compute optical properties
  if (return_full_output) {
    
    # Combine MAP estimates with fixed parameters
    par_complete <- c(map_values, par_fixed)
    
    # Prepare inputs
    inputs <- input_am03_sicf(par_complete, rrs)
    inputs$return_components <- TRUE  # Request all components
    
    # Run forward model with return_components = TRUE
    forward_result <- forward_am03_sicf(
      wavelength = inputs$wavelength,
      iop = inputs$iop,
      water_type = inputs$water_type,
      theta_view = inputs$theta_view,
      theta_sun = inputs$theta_sun,
      h_w = inputs$h_w,
      r_b = inputs$r_b,
      chl = inputs$chl,
      a_dg_443 = inputs$a_dg_443,
      phi_f = inputs$phi_f,
      include_sicf = inputs$include_sicf,
      lat = inputs$lat,
      lon = inputs$lon,
      date_time = inputs$date_time,
      return_components = TRUE
    )
    
    # Extract optical properties
    optical_properties <- list(
      Ed_0m = forward_result$Ed_0m,
      E0_0m = forward_result$E0_0m,
      Kd = forward_result$Kd,
      Ku_W = forward_result$Ku_W,
      Ku_B = forward_result$Ku_B,
      PAR = forward_result$PAR,
      wavelength = forward_result$wavelength
    )
    
    # Return full output
    return(list(
      par_estimates = par_estimates,
      rrs_modeled = forward_result$rrs_total,
      rrs_elastic = forward_result$rrs_elastic,
      rrs_sicf = forward_result$rrs_sicf,
      optical_properties = optical_properties,
      mcmc_output = out
    ))
    
  } else {
    # Return only parameter estimates
    return(par_estimates)
  }
}


#' Enhanced Inverse Gradient with Quantum Yield Retrieval and Optical Properties
#'
#' Extended version of inverse_gradient that retrieves quantum yield (phi_f) along
#' with standard OACs, and returns additional optical properties (Ed, E0, K, PAR).
#'
#' @param rrs Data frame with wavelength and rrs_0m columns
#' @param objective_fct Objective function name ("log-ll", "rss", etc.)
#' @param optim_mtd Optimization method ("L-BFGS-B", "levenberg-marqardt", etc.)
#' @param par_inversed Vector of parameter names to inverse (must include "phi_f")
#' @param par_fixed Named vector of fixed parameters
#' @param lower_b Numeric vector of lower bounds
#' @param upper_b Numeric vector of upper bounds
#' @param init_val Numeric vector of initial values
#' @param lat Latitude for Ed calculation (default = 49)
#' @param lon Longitude for Ed calculation (default = -68)
#' @param date_time POSIXct datetime in UTC (default = Sys.time())
#' @param return_full_output Logical, return optical properties? (default = TRUE)
#' @param verbose Logical, print progress messages
#'
#' @return If return_full_output = FALSE: named vector of estimates and uncertainties
#'         If return_full_output = TRUE: list with:
#'           - par_estimates: Optimal estimates with uncertainties
#'           - rrs_modeled: Modeled Rrs using optimal estimates
#'           - optical_properties: List with Ed_0m, E0_0m, Kd, Ku_W, Ku_B, PAR
#'           - optim_output: Full optimization output
#'
#' @export
inverse_gradient_sicf <- function(rrs,
                                   objective_fct = "log-ll",
                                   optim_mtd = "L-BFGS-B",
                                   par_inversed,
                                   par_fixed = NULL,
                                   lower_b = NULL,
                                   upper_b = NULL,
                                   init_val = NULL,
                                   lat = 49,
                                   lon = -68,
                                   date_time = Sys.time(),
                                   return_full_output = TRUE,
                                   verbose = TRUE) {
  
  # Validate that phi_f is in par_inversed
  if (!"phi_f" %in% par_inversed) {
    warning("'phi_f' not in par_inversed. Adding it automatically.")
    par_inversed <- c(par_inversed, "phi_f")
    lower_b <- c(lower_b, 0.005)
    upper_b <- c(upper_b, 0.03)
    if (!is.null(init_val)) init_val <- c(init_val, 0.02)
  }
  
  # Add lat, lon, date_time to par_fixed
  if (is.null(par_fixed)) par_fixed <- c()
  if (!"lat" %in% names(par_fixed)) par_fixed["lat"] <- lat
  if (!"lon" %in% names(par_fixed)) par_fixed["lon"] <- lon
  if (!"date_time" %in% names(par_fixed)) {
    par_fixed["date_time"] <- as.numeric(date_time)
  }
  
  # Use the standard inverse_gradient function
  result <- inverse_gradient(
    rrs = rrs,
    forward_model = "am03_sicf",
    objective_fct = objective_fct,
    optim_mtd = optim_mtd,
    par_inversed = par_inversed,
    par_fixed = par_fixed,
    lower_b = lower_b,
    upper_b = upper_b,
    init_val = init_val,
    verbose = verbose
  )
  
  # If return_full_output, compute optical properties
  if (return_full_output) {
    
    # Extract parameter estimates (without _sd suffix)
    par_names <- par_inversed
    par_values <- result[par_names]
    
    # Combine with fixed parameters
    par_complete <- c(par_values, par_fixed)
    
    # Prepare inputs
    inputs <- input_am03_sicf(par_complete, rrs)
    inputs$return_components <- TRUE
    
    # Run forward model
    forward_result <- forward_am03_sicf(
      wavelength = inputs$wavelength,
      iop = inputs$iop,
      water_type = inputs$water_type,
      theta_view = inputs$theta_view,
      theta_sun = inputs$theta_sun,
      h_w = inputs$h_w,
      r_b = inputs$r_b,
      chl = inputs$chl,
      a_dg_443 = inputs$a_dg_443,
      phi_f = inputs$phi_f,
      include_sicf = inputs$include_sicf,
      lat = inputs$lat,
      lon = inputs$lon,
      date_time = inputs$date_time,
      return_components = TRUE
    )
    
    # Extract optical properties
    optical_properties <- list(
      Ed_0m = forward_result$Ed_0m,
      E0_0m = forward_result$E0_0m,
      Kd = forward_result$Kd,
      Ku_W = forward_result$Ku_W,
      Ku_B = forward_result$Ku_B,
      PAR = forward_result$PAR,
      wavelength = forward_result$wavelength
    )
    
    # Return full output
    return(list(
      par_estimates = result,
      rrs_modeled = forward_result$rrs_total,
      rrs_elastic = forward_result$rrs_elastic,
      rrs_sicf = forward_result$rrs_sicf,
      optical_properties = optical_properties
    ))
    
  } else {
    return(result)
  }
}
