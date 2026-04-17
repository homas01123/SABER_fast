#' @title Forward model input preparation Registry
#' @description List of forward model input preparers available to SABER
.input_preparer_registry <- new.env(parent = emptyenv())

#' @title Forward Model Registry
#' @description List of forward models available to SABER
.forward_model_registry <- new.env(parent = emptyenv())

#' @title Objective Function Registry
#' @description List of objective functions available to SABER
.objective_function_registry <- new.env(parent = emptyenv())

#' Register a new input preparer
#' @param name Name of the preparer; must follow the pattern \code{input_<modelName>}.
#' @param fn   Function \code{f(par, rrs, par_meta)} returning a named list of
#'   forward-model inputs.
register_input_preparer <- function(name, fn) {
  .input_preparer_registry[[name]] <- fn
}

#' Register a new forward model
#' @param name Name of the model (e.g. \code{"am03"}).
#' @param fn   Function \code{f(inputs)} returning a numeric Rrs vector.
register_forward_model <- function(name, fn) {
  .forward_model_registry[[name]] <- fn
}

#' Register a new objective function
#' @param name Name of the objective (e.g. \code{"log-ll"}, \code{"lee99"}).
#' @param fn   Function \code{f(modelled, observed, par, ...)} returning a scalar.
register_objective_function <- function(name, fn) {
  .objective_function_registry[[name]] <- fn
}

#' Get registered model or error function
get_input_preparer <- function(name) {
  get0(name, envir = .input_preparer_registry)
}

get_forward_model <- function(name) {
  get0(name, envir = .forward_model_registry)
}

get_objective_function <- function(name) {
  get0(name, envir = .objective_function_registry)
}

#' List registered input preparers
#' @return Character vector of registered input preparer names.
#' @export
list_input_preparer <- function() {
  ls(.input_preparer_registry)
}

#' List registered forward models
#' @return Character vector of registered forward model names.
#' @export
list_forward_models <- function() {
  ls(.forward_model_registry)
}

#' List registered objective functions
#' @return Character vector of registered objective function names.
#' @export
list_objective_functions <- function() {
  ls(.objective_function_registry)
}
