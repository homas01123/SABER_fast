#' Perform deterministic inversion based on Gradient-based/Newtown-based optimization methods.
#'
#' @author Soham Mukherjee, Raphael Mabit
#'
#' @param rrs A tibble with the remote sensing reflectance data.
#' @param forward_model c("am03", "lee98")
#' @param objective_fct c("log-ll", "SSR", "lee99") SSR is not implemented in the code
#' @param optim_mtd c("nelder-mead", "BFGS", "CG", "L-BFGS-B", "sann", "brent","levenberg-marqardt", "auglag")
#' @param par_inversed a tibble with the initial parameter to be inverted.
#' @param par_fixed a tibble with parameter considered known, hence not retrieved during
#'  optimization. If defined here they must not be defined in `par_inversed`.
#' @param lower_b lower boundary for `optim_mtd = "L-BFGS-B"`.
#' If not provided will be calculated in `parse_inverse_parameter`.
#' @param init_val Initial values for parameters to be inverted.
#' @param upper_b same as `lower_b` but for maximum possible values..
#' @param verbose Boolean, if TRUE prints additional information
#' @return A named vector with the maximum likelihood estimates and their standard deviations.
#'
#' @export
inverse_gradient <- function(
    rrs,
    forward_model,
    objective_fct,
    optim_mtd,
    par_inversed,
    par_fixed = NULL,
    lower_b = NULL,
    init_val = NULL,
    upper_b = NULL,
    verbose = F) {
  rlang::inform(paste0("\033[0;33m", "###################################################################", "\033[0m", "\n"))
  rlang::inform(paste0("\033[0;39m", "########### ALL GOOD THINGS ARE WILD & FREE, LET'S RUN FREE #######", "\033[0m", "\n"))
  rlang::inform(paste0("\033[0;32m", "###################################################################", "\033[0m", "\n"))

  minimization_fct <- objective_factory(
    model = forward_model,
    objective = objective_fct,
    rrs_observed = rrs,
    par_fixed = par_fixed,
    par_inversed = par_inversed,
    minimize = TRUE
  )

  # Instantiate initial values
  params <- parse_inverse_parameter(
    par_df = data.frame("name" = par_inversed, "value" = init_val),
    optim_mtd=optim_mtd,
    lower_b = lower_b,
    upper_b = upper_b,
    verbose = verbose
  )

  par <- params$par
  names(par) <- params$names
  lower_b <- params$lower
  upper_b <- params$upper
  parscale <- params$parscale

  # Optimization ------------------------------------------------------------

  start.time <- Sys.time()

  if (optim_mtd == "L-BFGS-B") {
    optim_result <- optim(
      par = par,
      fn = minimization_fct,
      method = optim_mtd,
      lower = lower_b,
      upper = upper_b,
      control = list(
        parscale = parscale,
        fnscale = 1,  # Since we're minimizing
        maxit = 1000
      )
    )
  }

  if (optim_mtd == "Nelder-Mead" |
    optim_mtd == "SANN" |
    optim_mtd == "Brent") {
    optim_result <- optim(
      par = par,
      fn = minimization_fct,
      method = optim_mtd,
      control = list(parscale = parscale),
      hessian = FALSE
    )
  }

  if (optim_mtd == "levenberg-marqardt") {
    # Note: Does not support bounds directly

    lm_result <- marqLevAlg::marqLevAlg(
      b = par,
      fn = minimization_fct,
      minimize = TRUE,  # Explicitly set minimize flag
      print.info = F
    )

    n_params <- length(par)
    # Fill in the lower triangular part from the vector
    varcov_matrix[lower.tri(varcov_matrix, diag = TRUE)] <- lm_result$v

    # Make it symmetric by copying lower triangle to upper triangle
    varcov_matrix[upper.tri(varcov_matrix)] <- t(varcov_matrix)[upper.tri(varcov_matrix)]

    # Create optim_result structure compatible with other methods
    optim_result <- list(
      par = lm_result$b,
      value = lm_result$fn.value,
      convergence = ifelse(lm_result$istop == 1, 0, 1),  # 0 = success, 1 = failure
      message = lm_result$message,

      # Store variance-covariance matrix for uncertainty calculation
      varcov = varcov_matrix
    )
  }

  if (optim_mtd == "auglag") {
    print("Augmented Lagriangian with equality constraints will be used for inversion")

    get_fraction_indices <- function(par_names, prefix = "rb_") {
      which(grepl(paste0("^", prefix), par_names))
    }

    # Identify indices dynamically
    fraction_indices <- get_fraction_indices(names(par))

    # bounds for areal fractions
    fheq <- function(pars) {
      sum(pars[fraction_indices]) - 1 # must sum to 1
    }

    fhin <- function(pars) {
      pars[fraction_indices] # All should be ≥ 0
    }

    optim_result <- alabama::auglag(
      fn = minimization_fct,
      par = par,
      heq = fheq,
      hin = fhin,
      control.outer = list(trace = F, method = "nlminb")
    )
  }

  # Calculate uncertainty ---------------------------------------------------

  # Calculate hessian matrix for var-covar matrix
  if (optim_mtd == "auglag") {
    hessian_inverse <- optim_result$hessian
  } else if (optim_mtd == "levenberg-marqardt") {
    hessian_inverse <- optim_result$varcov
  } else {
    # For standard errors, we need Hessian of POSITIVE log-likelihood
    # Since minimization_fct returns negative log-likelihood, negate the Hessian
    hessian_neg_ll <- numDeriv::hessian(
      x = optim_result$par,
      func = minimization_fct
    )
    hessian_inverse <- hessian_neg_ll  # Negate to get Hessian of positive log-likelihood
  }

  if (verbose) {
    rownames(hessian_inverse) <- par_inversed
    colnames(hessian_inverse) <- par_inversed
    rlang::inform(paste0("\033[0;32m", "#################### VAR-COV HESSIAN MATRIX #########################", "\033[0m", "\n"))
    prmatrix(hessian_inverse)
  }

  param_estimate <- optim_result$par

  param_sd <- tryCatch({
    # Check for singular matrix
    if (abs(det(hessian_inverse)) < 1e-5) {
      warning("Hessian is nearly singular - using pseudoinverse")
      varcov <- MASS::ginv(hessian_inverse)  # Generalized inverse
    } else {
      varcov <- solve(hessian_inverse)
    }
    sqrt(abs(diag(varcov)))  # Use abs() to handle numerical errors
  },
  error = myFun
  )

  # param_sd <- tryCatch(
  #   {
  #     sqrt(diag(solve(hessian_inverse)))
  #   }, # solve for diagonal elements to get sd
  #   objective = NA
  # )

  end.time <- Sys.time()

  if (!is.numeric(param_sd)) {
    rlang::warn(
      paste0("\033[0;31m", "Failed to calculate diagonal of hessian from
             high degree of correlation, coerce to NA", "\033[0m", "\n")
    )
    param_sd <- rep(NA, length(param_estimate))
  }

  # Maximum Likelihood Estimates - Convert to named vector format like inverse_mcmc
  param_names <- c(par_inversed, paste0(par_inversed, "_sd"))

  par_estimates <- stats::setNames(
    c(param_estimate, param_sd),
    param_names
  )

  # # Maximum Likelihood Estimates
  # mle <- tibble(
  #   "name" = par_inversed,
  #   "estimate" = param_estimate,
  #   "sd" = param_sd
  # )

  if (verbose) {
    if (optim_result$convergence == 0) {
      # convergence <- "TRUE"
      rlang::inform(paste0("\033[0;32m", "CONVERGENCE: GLOBAL", "\033[0m", "\n"))
    } else {
      # convergence = "FALSE"
      rlang::inform(paste0("\033[0;34m", "CONVERGENCE: LOCAL", "\033[0m", "\n"))
    }

    time_taken <- end.time - start.time
    rlang::inform(glue::glue("time.elapsed: ", time_taken))
    # return(list(mle, "convergence"= convergence))
  }

  return(par_estimates)
}


#' @export
myFun <- function(x) {
  NA
}

