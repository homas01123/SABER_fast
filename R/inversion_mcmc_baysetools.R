#' Build a lognormal BayesianTools prior for deep-water inversion parameters
#'
#' Every parameter in par_inversed is given an independent lognormal prior whose
#' log-space mean is the geometric midpoint of [lower_b, upper_b] and whose
#' log-space sd spans half the log range.
#'
#' @param par_inversed character vector of parameter names to invert
#' @param lower_b named numeric vector of lower bounds (must cover all par_inversed)
#' @param upper_b named numeric vector of upper bounds (must cover all par_inversed)
#' @return a BayesianTools prior object
#' @export
make_bt_prior_deep <- function(par_inversed, lower_b, upper_b) {
  lower   <- lower_b[par_inversed]
  upper   <- upper_b[par_inversed]
  log_mid <- 0.5 * (log(lower) + log(upper))
  log_sd  <- 0.5 * (log(upper) - log(lower))
  BayesianTools::createPrior(
    density = function(par) {
      sum(dnorm(log(par), mean = log_mid, sd = log_sd, log = TRUE) - log(par))
    },
    sampler = function(n = 1) {
      mat <- matrix(NA_real_, nrow = n, ncol = length(par_inversed))
      for (i in seq_len(n))
        mat[i, ] <- exp(rnorm(length(par_inversed), mean = log_mid, sd = log_sd))
      mat
    },
    lower = lower,
    upper = upper
  )
}


#' Build a BayesianTools prior for 2-class shallow water (mix_sand parametrisation)
#'
#' OAC parameters receive independent lognormal priors.
#' mix_sand in [0, 1] receives a Beta(2, 2) prior (broad, symmetric, zero-avoiding).
#'
#' @param par_inversed character vector of parameter names to invert
#' @param lower_b named numeric vector of lower bounds
#' @param upper_b named numeric vector of upper bounds
#' @return a BayesianTools prior object
#' @export
make_bt_prior_shallow_2class <- function(par_inversed, lower_b, upper_b) {
  lower   <- lower_b[par_inversed]
  upper   <- upper_b[par_inversed]
  log_mid <- ifelse(par_inversed == "mix_sand", NA_real_, 0.5 * (log(lower) + log(upper)))
  log_sd  <- ifelse(par_inversed == "mix_sand", NA_real_, 0.5 * (log(upper) - log(lower)))
  BayesianTools::createPrior(
    density = function(par) {
      total <- 0
      for (j in seq_along(par)) {
        if (par_inversed[j] == "mix_sand") {
          total <- total + dbeta(par[j], shape1 = 2, shape2 = 2, log = TRUE)
        } else {
          total <- total + dnorm(log(par[j]), mean = log_mid[j], sd = log_sd[j], log = TRUE) - log(par[j])
        }
      }
      total
    },
    sampler = function(n = 1) {
      mat <- matrix(NA_real_, nrow = n, ncol = length(par_inversed))
      for (i in seq_len(n)) {
        for (j in seq_along(par_inversed)) {
          if (par_inversed[j] == "mix_sand") {
            mat[i, j] <- rbeta(1, 2, 2)
          } else {
            mat[i, j] <- exp(rnorm(1, mean = log_mid[j], sd = log_sd[j]))
          }
        }
      }
      mat
    },
    lower = lower,
    upper = upper
  )
}


#' Build a BayesianTools prior for N-class shallow water (soft Dirichlet + lognormal OACs)
#'
#' Parameters whose names start with r_rs_b_ are treated as benthic fractions and
#' receive a soft Dirichlet prior: (alpha-1)*sum(log(b)) - 50*(sum(b)-1)^2.
#' All other parameters receive independent lognormal priors.
#'
#' @param par_inversed character vector of parameter names to invert
#' @param lower_b named numeric vector of lower bounds
#' @param upper_b named numeric vector of upper bounds
#' @param alpha Dirichlet concentration parameter (default 2)
#' @return a BayesianTools prior object
#' @export
make_bt_prior_shallow_nclass <- function(par_inversed, lower_b, upper_b, alpha = 2) {
  is_benthic <- grepl("^r_rs_b_", par_inversed)
  lower   <- lower_b[par_inversed]
  upper   <- upper_b[par_inversed]
  log_mid <- ifelse(is_benthic, NA_real_, 0.5 * (log(lower) + log(upper)))
  log_sd  <- ifelse(is_benthic, NA_real_, 0.5 * (log(upper) - log(lower)))
  BayesianTools::createPrior(
    density = function(par) {
      b_vals <- par[is_benthic]
      oac_lp <- sum(
        dnorm(log(par[!is_benthic]), mean = log_mid[!is_benthic], sd = log_sd[!is_benthic], log = TRUE) -
          log(par[!is_benthic])
      )
      dir_lp <- (alpha - 1) * sum(log(b_vals + 1e-9)) - 50 * (sum(b_vals) - 1)^2
      oac_lp + dir_lp
    },
    sampler = function(n = 1) {
      nb  <- sum(is_benthic)
      mat <- matrix(NA_real_, nrow = n, ncol = length(par_inversed))
      for (i in seq_len(n)) {
        for (j in which(!is_benthic))
          mat[i, j] <- exp(rnorm(1, mean = log_mid[j], sd = log_sd[j]))
        g <- rgamma(nb, shape = alpha, rate = 1)
        g <- g / sum(g)
        mat[i, which(is_benthic)] <- g
      }
      mat
    },
    lower = lower,
    upper = upper
  )
}



#' Build an adaptive lognormal BayesianTools prior centred on per-observation spectral estimates
#'
#' Like \code{make_bt_prior_deep} but the log-space mean for each parameter is set
#' from \code{modes} (e.g. OC3 / band-ratio estimates from \code{estimate_prior_modes})
#' rather than the geometric midpoint of the bounds.
#'
#' @param par_inversed character vector of parameter names to invert
#' @param lower_b named numeric vector of lower bounds (must cover all par_inversed)
#' @param upper_b named numeric vector of upper bounds (must cover all par_inversed)
#' @param modes named numeric vector of prior mode estimates (must cover all par_inversed)
#' @return a BayesianTools prior object
#' @export
make_bt_prior_adaptive <- function(par_inversed, lower_b, upper_b, modes) {
  lower   <- lower_b[par_inversed]
  upper   <- upper_b[par_inversed]
  modes_c <- pmax(lower * 1.01, pmin(upper * 0.99, modes[par_inversed]))
  log_mid <- log(modes_c)
  log_sd  <- pmax(0.3,
    pmin(abs(log(upper) - log_mid), abs(log_mid - log(lower))) * 0.5
  )
  BayesianTools::createPrior(
    density = function(par) {
      sum(dnorm(log(par), mean = log_mid, sd = log_sd, log = TRUE) - log(par))
    },
    sampler = function(n = 1) {
      mat <- matrix(NA_real_, nrow = n, ncol = length(par_inversed))
      for (i in seq_len(n))
        for (j in seq_along(par_inversed)) {
          repeat {
            v <- exp(rnorm(1, log_mid[j], log_sd[j]))
            if (v >= lower[j] && v <= upper[j]) { mat[i, j] <- v; break }
          }
        }
      mat
    },
    lower = lower,
    upper = upper
  )
}


#' Build an adaptive 2-class shallow-water BayesianTools prior
#'
#' OAC parameters receive independent lognormal priors centred on \code{modes}.
#' \code{mix_sand} in [0, 1] receives a Beta(2, 2) prior.
#' A \code{repeat\{\}} bounded sampler prevents out-of-bounds start values.
#'
#' @param par_inversed character vector of parameter names to invert
#' @param lower_b named numeric vector of lower bounds
#' @param upper_b named numeric vector of upper bounds
#' @param modes named numeric vector of prior mode estimates (OAC parameters only)
#' @return a BayesianTools prior object
#' @export
make_bt_prior_shallow_2class_adaptive <- function(par_inversed, lower_b, upper_b, modes) {
  mix_idx <- which(par_inversed == "mix_sand")
  oac_idx <- which(par_inversed != "mix_sand")
  lower   <- lower_b[par_inversed]
  upper   <- upper_b[par_inversed]
  modes_c <- pmax(lower[oac_idx] * 1.01, pmin(upper[oac_idx] * 0.99, modes[par_inversed[oac_idx]]))
  log_mid <- log(modes_c)
  log_sd  <- pmax(0.3,
    pmin(abs(log(upper[oac_idx]) - log_mid), abs(log_mid - log(lower[oac_idx]))) * 0.5
  )
  BayesianTools::createPrior(
    density = function(par) {
      oac_lp <- sum(dnorm(log(par[oac_idx]), mean = log_mid, sd = log_sd, log = TRUE) -
                    log(par[oac_idx]))
      mix_lp <- dbeta(par[mix_idx], 2, 2, log = TRUE)
      oac_lp + mix_lp
    },
    sampler = function(n = 1) {
      mat <- matrix(NA_real_, nrow = n, ncol = length(par_inversed))
      for (i in seq_len(n)) {
        for (k in seq_along(oac_idx)) {
          j <- oac_idx[k]
          repeat {
            v <- exp(rnorm(1, log_mid[k], log_sd[k]))
            if (v >= lower[j] && v <= upper[j]) { mat[i, j] <- v; break }
          }
        }
        mat[i, mix_idx] <- max(0.001, min(0.999, rbeta(1, 2, 2)))
      }
      mat
    },
    lower = lower,
    upper = upper
  )
}


#' SABER inverse model using MCMC sampling
#'
#' Retrieve optically deep water OSCs and benthic variables from input Rrs and wavelength
#'
#' @param rrs data-frame of wavelengths [nm] and sub-surface Rrs (must be named as rrs_0m) [1/sr]
#' @param forward_model the SA forward model to be used (e.g. "am03")
#' @param par_inversed vector of parameter names to be inversed (e.g. c("chl", "a_dg_440", "bb_p_550"))
#' @param prior prior function (see make_prior_bundle) or can set as NULL for uniform prior
#' @param lower numeric vector containing lower bounds for each parameter in par_inversed. The order must match par_inversed.
#' @param upper numeric vector containing upper bounds for each parameter in par_inversed. The order must match par_inversed.
#' @param best numeric vector containing best guess values for each parameter in par_inversed. The order must match par_inversed.
#' @param par_fixed named list of fixed parameters to be passed to the forward model (e.g. c(water_type = 2, theta_sun = 30, theta_view = 0...))
#' @param iterations number of MCMC iterations (default = 10000, recommended > 15000)
#' @param burnin number of burnin iterations (default = 2000)
#' @param sampler MCMC sampler to be used (default = "DEzs", see BayesianTools documentation for other options)
#' @param return_full_output logical, return optical properties when forward_model = "am03_sicf"? (default = FALSE)

#' @return numeric vector of parameter estimates and their standard deviations (e.g. c(chl = 2.5, chl_sd = 0.3, a_dg_440 = 0.1, a_dg_440_sd = 0.02, bb_p_550 = 0.01, bb_p_550_sd = 0.003))
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
    return_full_output = FALSE,
    spectral_weights = NULL) {

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
    model            = forward_model,
    objective        = "log-ll",
    rrs_observed     = rrs,
    par_inversed     = par_inversed,
    par_fixed        = par_fixed,
    par_meta         = par_meta,
    spectral_weights = spectral_weights
  )

  setup <- BayesianTools::createBayesianSetup(
    prior      = prior,
    likelihood = likelihood,
    lower      = if (is.null(prior)) lower else NULL,
    best       = if (is.null(prior)) best  else NULL,
    upper      = if (is.null(prior)) upper else NULL,
    names      = par_inversed,
    parallel   = FALSE
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
