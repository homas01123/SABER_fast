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

#' @return numeric vector of parameter estimates and their standard deviations (e.g. c(chl = 2.5, chl_sd = 0.3, a_g_440 = 0.1, a_g_440_sd = 0.02, bb_p_550 = 0.01, bb_p_550_sd = 0.003))
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
    sampler = "DEzs") {


  likelihood <- objective_factory(
    model = forward_model,
    objective = "log-ll",
    rrs_observed = rrs,
    par_inversed = par_inversed,
    par_fixed = par_fixed
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
    ~ apply(.x[, ncol(.x) - 3:ncol(.x)], 2, sd)
  )

  estimates_sd <- colMeans(estimates_sd)

  par_estimates <- stats::setNames(
    c(BayesianTools::MAP(out)[[1]], estimates_sd),
    c(names(BayesianTools::MAP(out)[[1]]), paste0(names(BayesianTools::MAP(out)[[1]]), "_sd"))
  )

  return(par_estimates)
}
