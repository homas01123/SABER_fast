make_prior <- function(sample_data, dist = "weibull", name = "param") {
  fit <- fitdistrplus::fitdist(sample_data, dist)
  list(
    fit = fit,
    density = function(x) {
      do.call(paste0("d", dist), list(x = x, shape = fit$estimate[1], scale = fit$estimate[2], log = TRUE))
    },
    sampler = function(n = 1) {
      do.call(paste0("r", dist), list(n = n, shape = fit$estimate[1], scale = fit$estimate[2]))
    },
    name = name
  )
}


make_prior_bundle <- function(priors, lower, upper, best_guess = NULL) {
  density_fn <- function(par) {
    sum(mapply(function(p, val) p$density(val), priors, par))
  }

  sampler_fn <- function(n = 1) {
    matrix(unlist(lapply(priors, function(p) p$sampler(n))), ncol = length(priors))
  }

  list(
    prior = BayesianTools::createPrior(
      density = density_fn,
      sampler = sampler_fn,
      lower = lower,
      upper = upper,
      best = best_guess
    ),
    names = sapply(priors, function(p) p$name)
  )
}



#' SABER inverse model using MCMC sampling
#'
#' Retrieve optically deep water OSCs and benthic variables from input Rrs and wavelength
#'
#' @param rrs data-frame of wavelengths [nm] and sub-surface Rrs (must be named as rrs_0m) [1/sr]
#' @param forward_model the SA forward model to be used (e.g. "am03")
#' @param par_inversed vector of parameter names to be inversed (e.g. c("chl", "a_g_440", "bb_p_550"))
#' @param prior prior function (see make_prior_bundle) or can set as NULL for uniform prior
#' @param lower numeric vector containing lower bounds for each parameter in par_inversed. The order must match par_inversed.
#' @param upper numeric vector containing upper bounds for each parameter in par_inversed. The order must match par_inversed.
#' @param best numeric vector containing best guess values for each parameter in par_inversed. The order must match par_inversed.
#' @param par_fixed named list of fixed parameters to be passed to the forward model (e.g. c(water_type = 2, theta_sun = 30, theta_view = 0...))
#' @param iterations number of MCMC iterations (default = 10000, recommended > 15000)
#' @param burnin number of burnin iterations (default = 2000)
#' @param sampler MCMC sampler to be used (default = "DEzs", see BayesianTools documentation for other options)
#' @param return_full_output logical, return optical properties when forward_model = "am03_sicf"? (default = FALSE)

#' @return numeric vector of parameter estimates and their standard deviations (e.g. c(chl = 2.5, chl_sd = 0.3, a_g_440 = 0.1, a_g_440_sd = 0.02, bb_p_550 = 0.01, bb_p_550_sd = 0.003))
#'         If forward_model = "am03_sicf" and return_full_output = TRUE, returns list with par_estimates, rrs_modeled, rrs_elastic, rrs_sicf, optical_properties, mcmc_output
#'
#' @references Mukherjee, S., Mabit, R. and Bélanger, S. (2025), A semi-analytical Bayesian estimate retrieval algorithm for the inversion of remote-sensing reflectance in optically deep and shallow waters. Limnol Oceanogr Methods. https://doi.org/10.1002/lom3.70004
#'
#' @export
inverse_mcmc <- function(
    rrs,
    forward_model,
    par_inversed,
    prior = NULL,
    lower = NULL,
    best = NULL,
    upper = NULL,
    par_fixed = NULL,
    iterations = 10000,
    burnin = 2000,
    sampler = "DEzs",
    return_full_output = FALSE) {

  # Auto-add phi_f if using am03_sicf model and phi_f not in par_inversed
  if (forward_model == "am03_sicf" && !"phi_f" %in% par_inversed) {
    warning("'phi_f' not in par_inversed when using am03_sicf model. Adding it automatically.")
    par_inversed <- c(par_inversed, "phi_f")
    lower <- c(lower, 0.005)
    upper <- c(upper, 0.03)
    if (!is.null(best)) best <- c(best, 0.02)
  }

  # Separate numeric and non-numeric parameters from par_fixed
  par_meta <- NULL
  if (!is.null(par_fixed)) {
    # Identify non-numeric metadata parameters
    meta_params <- c("sicf_model", "depth_integration")
    meta_names <- intersect(names(par_fixed), meta_params)
    
    if (length(meta_names) > 0) {
      par_meta <- par_fixed[meta_names]
      par_fixed <- par_fixed[!names(par_fixed) %in% meta_names]
    }
  }

  likelihood <- objective_factory(
    model = forward_model,
    objective = "log-ll",
    rrs_observed = rrs,
    par_inversed = par_inversed,
    par_fixed = par_fixed,
    par_meta = par_meta
  )

  setup <- BayesianTools::createBayesianSetup(
    prior = prior,
    likelihood = likelihood,
    lower = lower,
    best = best,
    upper = upper,
    names = par_inversed,
    parallel = FALSE
  )

  BayesianTools::checkBayesianSetup(setup)

  out <- BayesianTools::runMCMC(
    bayesianSetup = setup,
    settings = list(
      iterations = iterations,
      burnin = burnin,
      message = TRUE
    ),
    sampler = sampler
  )

  estimates_sd <- purrr::map_df(
    .x = out[["chain"]],
    ~ apply(.x[, 1:(ncol(.x) - 3)], 2, sd)
  )

  estimates_sd <- colMeans(estimates_sd)
  
  map_values <- BayesianTools::MAP(out)[[1]]

  par_estimates <- stats::setNames(
    c(map_values, estimates_sd),
    c(names(map_values), paste0(names(map_values), "_sd"))
  )

  # If using am03_sicf model with full output, compute optical properties
  if (forward_model == "am03_sicf" && return_full_output) {
    
    # Combine MAP estimates with fixed parameters
    # Convert par_fixed from list to numeric vector if needed
    if (!is.null(par_fixed)) {
      if (is.list(par_fixed)) {
        par_fixed_vec <- unlist(par_fixed)
        # Preserve names from list
        if (is.null(names(par_fixed_vec))) {
          names(par_fixed_vec) <- names(par_fixed)
        }
      } else {
        par_fixed_vec <- par_fixed
      }
      par_complete <- c(map_values, par_fixed_vec)
      par_complete <- par_complete[order(names(par_complete))]
    } else {
      par_complete <- map_values
    }
    
    # Prepare inputs
    inputs <- input_am03_sicf(par_complete, rrs, par_meta)
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
      sicf_model = inputs$sicf_model,
      depth_integration = inputs$depth_integration,
      return_components = TRUE
    )
    
    # Extract optical properties
    optical_properties <- list(
      Ed_0m = forward_result$Ed_0m,
      E0_0m = forward_result$E0_0m,
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
  }

  return(par_estimates)
}
