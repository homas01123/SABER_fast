#' Build an objective / likelihood function from a forward model
#'
#' Combines a registered forward model with an objective function into a single
#' callable \code{f(par)}.  Used internally by [inverse_mcmc()] and
#' [inverse_gradient()]; exposed for advanced workflows.
#'
#' @param model       Forward model name: \code{"am03"}, \code{"am03_sicf"}, or
#'   \code{"lee98"}.
#' @param objective   Objective: \code{"log-ll"} (Gaussian log-likelihood) or
#'   \code{"lee99"} (spectral error index).
#' @param rrs_observed  Data frame with columns \code{wavelength} and
#'   \code{rrs_0m}.
#' @param par_inversed  Character vector of inverted parameter names.
#' @param par_fixed   Named list of fixed parameters passed to the forward
#'   model.
#' @param par_meta    Named list of non-numeric metadata
#'   (\code{sicf_model}, \code{depth_integration}).  Used only for
#'   \code{"am03_sicf"}.
#' @param minimize    Logical.  Return \code{-objective} for minimisation?
#'   Default \code{FALSE} (returns log-likelihood for MCMC).
#' @param log_prior_fn  Optional \code{function(par)} returning the log-prior.
#'   Active only when \code{minimize = TRUE}; ignored on the MCMC path.
#' @param spectral_weights  Optional named numeric per-band weight vector
#'   (see \code{\link{create_spectral_weights}}).
#'
#' @return A function \code{f(par)} → scalar objective value.
#'
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

    # Guard: return a penalty when parameters are physically invalid.
    # This is critical for bound-free optimisers (e.g. Levenberg-Marquardt) that
    # can walk into negative IOP space.  Rather than crashing inside the forward
    # model we return a large cost (minimisation) or -Inf (likelihood), which
    # pushes the optimiser back towards the feasible region.
    inputs <- tryCatch({
      # Pass par_meta only for models that support it (am03_sicf)
      if (model == "am03_sicf" && !is.null(par_meta)) {
        prepare_input(par_complete, rrs_observed, par_meta)
      } else {
        prepare_input(par_complete, rrs_observed)
      }
    }, error = function(e) NULL)

    if (is.null(inputs)) {
      return(if (minimize) .Machine$double.xmax / 2 else -Inf)
    }

    rrs_modeled <- tryCatch(
      forward_model(inputs),
      error = function(e) NULL
    )

    if (is.null(rrs_modeled) || !all(is.finite(rrs_modeled))) {
      return(if (minimize) .Machine$double.xmax / 2 else -Inf)
    }

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
