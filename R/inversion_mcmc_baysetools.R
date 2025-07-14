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
