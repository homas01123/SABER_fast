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

objective_factory <- function(model, objective, rrs_observed, par_inversed, par_fixed = NULL, par_meta = NULL, minimize = FALSE) {
  prepare_input <- get_input_preparer(paste0("input_", model))
  forward_model <- get_forward_model(model)
  objective_function <- get_objective_function(objective)

  if (is.null(prepare_input)) stop(paste("Unknown prepare input for: ", paste0("input_", model)))
  if (is.null(forward_model)) stop(paste("Unknown forward model: ", model))
  if (is.null(objective_function)) stop(paste("Unknown objective function: ", objective))

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
    
    # For minimization (gradient methods), negate log-likelihood
    if (minimize && objective == "log-ll") {
      return(-result)
    } else {
      return(result)
    }
  }
}
