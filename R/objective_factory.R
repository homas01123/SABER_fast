#' Make objective function from forward model and objective
#'
#' @param model Name of the forward model
#' @param objective Name of the objective function
#' @param rrs_observed Observed Rrs data
#' @param par_fixed Named list or data frame of fixed parameters
#' @param minimize Logical, if TRUE returns negative objective for minimization
#'
#' @return A function returning the model objective/likelihood to be used in an optimization process
#' @export

objective_factory <- function(model, objective, rrs_observed, par_inversed, par_fixed = NULL, minimize = FALSE) {
  prepare_input <- get_input_preparer(paste0("input_", model))
  forward_model <- get_forward_model(model)
  objective_function <- get_objective_function(objective)

  if (is.null(prepare_input)) stop(paste("Unknown prepare input for: ", paste0("input_", model)))
  if (is.null(forward_model)) stop(paste("Unknown forward model: ", model))
  if (is.null(objective_function)) stop(paste("Unknown objective function: ", objective))

  complete_par <- function(par) {
    if (!is.null(par_fixed)) {
      # if (any(par_fixed$name %in% names(par))) {
      #   stop("Parameters cannot be defined in both `par_init` and `par_fixed`.")
      # }

      par <- c(par, par_fixed)
    }

    par[order(names(par))]
  }

  function(par) {
    # Restore names destroyed by mcmc sampling
    names(par) <- par_inversed
    par <- complete_par(par)
    inputs <- prepare_input(par, rrs_observed)
    rrs_modeled <- forward_model(inputs)

    result <- objective_function(
      modelled = rrs_modeled,
      observed = rrs_observed$rrs_0m,
      par = par
    )
    
    # For minimization (gradient methods), negate log-likelihood
    if (minimize && objective == "log-ll") {
      return(-result)
    } else {
      return(result)
    }
  }
}
