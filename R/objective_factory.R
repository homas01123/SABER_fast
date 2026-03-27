#' Make objective function from forward model and objective
#'
#' @param model Name of the forward model
#' @param objective Name of the objective function
#' @param rrs_observed Observed Rrs data
#' @param par_fixed Named list or data frame of fixed parameters
#' @param minimize Logical, if TRUE returns negative objective for minimization
#' @param log_prior_fn Optional function `function(par)` returning the log-prior
#'   density for the **inverted** parameters (as a named numeric vector).
#'   When supplied and `minimize = TRUE`, the factory returns
#'   `-(log_ll + log_prior)`, i.e. MAP optimisation.
#'   When `minimize = FALSE` (MCMC likelihood path) the prior is ignored here
#'   because BayesianTools handles it via `createBayesianSetup(prior = ...)`.
#'
#' @return A function returning the model objective/likelihood to be used in an optimization process
#' @export

objective_factory <- function(model, objective, rrs_observed, par_inversed, par_fixed = NULL, par_meta = NULL, minimize = FALSE, log_prior_fn = NULL, spectral_weights = NULL) {
  prepare_input <- get_input_preparer(paste0("input_", model))
  forward_model <- get_forward_model(model)
  objective_fn_raw <- get_objective_function(objective)

  if (is.null(prepare_input)) stop(paste("Unknown prepare input for: ", paste0("input_", model)))
  if (is.null(forward_model)) stop(paste("Unknown forward model: ", model))
  if (is.null(objective_fn_raw)) stop(paste("Unknown objective function: ", objective))

  # When spectral_weights are supplied, wrap the objective to inject them.
  # The weights vector is captured in the closure and never travels through par.
  objective_function <- if (!is.null(spectral_weights)) {
    force(spectral_weights)
    function(modelled, observed, par)
      objective_fn_raw(modelled, observed, par, weights = spectral_weights)
  } else {
    objective_fn_raw
  }

  complete_par <- function(par) {
    # par is numeric vector from optimizer
    names(par) <- par_inversed
    
    # Add numeric fixed parameters
    if (!is.null(par_fixed)) {
      # If par_fixed is a list, convert numeric portion to named vector
      if (is.list(par_fixed)) {
        par_fixed_vec <- unlist(par_fixed)
        # Ensure names are preserved from the list
        if (is.null(names(par_fixed_vec))) {
          names(par_fixed_vec) <- names(par_fixed)
        }
      } else {
        par_fixed_vec <- par_fixed
      }
      # Combine as numeric vector with names preserved
      par <- c(par, par_fixed_vec)
    }
    
    par[order(names(par))]
  }

  function(par) {
    par_complete <- complete_par(par)
    
    # Pass par_meta only for models that support it (am03_sicf)
    if (model == "am03_sicf" && !is.null(par_meta)) {
      inputs <- prepare_input(par_complete, rrs_observed, par_meta)
    } else {
      inputs <- prepare_input(par_complete, rrs_observed)
    }
    
    rrs_modeled <- forward_model(inputs)

    result <- objective_function(
      modelled = rrs_modeled,
      observed = rrs_observed$rrs_0m,
      par = par_complete
    )
    
    # For minimization (gradient methods): return -(log_ll + log_prior)
    # This is MAP estimation — the prior is a regularisation penalty.
    # When minimize = FALSE (MCMC likelihood path) the prior is intentionally
    # excluded here; BayesianTools combines likelihood + prior internally via
    # createBayesianSetup(prior = ...).  Do NOT add log_prior_fn there.
    if (minimize && objective == "log-ll") {
      log_prior_val <- if (!is.null(log_prior_fn)) {
        # par_complete is named and sorted; extract the inversed subset only
        par_for_prior <- par_complete[par_inversed]
        lp <- tryCatch(log_prior_fn(par_for_prior), error = function(e) -Inf)
        if (is.na(lp) || !is.finite(lp)) -Inf else lp
      } else 0.0
      return(-(result + log_prior_val))
    } else {
      return(result)
    }
  }
}
