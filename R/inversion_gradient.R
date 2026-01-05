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
#' @param return_full_output logical, return optical properties when forward_model = "am03_sicf"? (default = FALSE)
#' @return A named vector with the maximum likelihood estimates and their standard deviations.
#'         If forward_model = "am03_sicf" and return_full_output = TRUE, returns list with par_estimates, rrs_modeled, rrs_elastic, rrs_sicf, optical_properties
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
    verbose = F,
    return_full_output = FALSE) {
      
  rlang::inform(paste0("\033[0;33m", "###################################################################", "\033[0m", "\n"))
  rlang::inform(paste0("\033[0;39m", "########### ALL GOOD THINGS ARE WILD & FREE, LET'S RUN FREE #######", "\033[0m", "\n"))
  rlang::inform(paste0("\033[0;32m", "###################################################################", "\033[0m", "\n"))

  # Auto-add phi_f if using am03_sicf model and phi_f not in par_inversed
  if (forward_model == "am03_sicf" && !"phi_f" %in% par_inversed) {
    warning("'phi_f' not in par_inversed when using am03_sicf model. Adding it automatically.")
    par_inversed <- c(par_inversed, "phi_f")
    lower_b <- c(lower_b, 0.005)
    upper_b <- c(upper_b, 0.03)
    if (!is.null(init_val)) init_val <- c(init_val, 0.02)
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

  minimization_fct <- objective_factory(
    model = forward_model,
    objective = objective_fct,
    rrs_observed = rrs,
    par_fixed = par_fixed,
    par_meta = par_meta,
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
    if (verbose) {
      message("Using Augmented Lagrangian with equality constraints for inversion")
    }

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
    message("\n#################### VAR-COV HESSIAN MATRIX #########################\n")
    prmatrix(hessian_inverse)
    message(paste0("Absolute determinant of Hessian: ", abs(det(hessian_inverse))))
  }

  param_estimate <- optim_result$par

  param_sd <- tryCatch({
    # Check for singular matrix
    if (abs(det(hessian_inverse)) < 1e-5 | abs(det(hessian_inverse)) > 1e10) {
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

  # If using am03_sicf model with full output, compute optical properties
  if (forward_model == "am03_sicf" && return_full_output) {
    
    # Extract parameter estimates (without _sd suffix)
    par_names <- par_inversed
    par_values <- par_estimates[par_names]
    
    # Combine with fixed parameters
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
      par_complete <- c(par_values, par_fixed_vec)
      par_complete <- par_complete[order(names(par_complete))]
    } else {
      par_complete <- par_values
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
      optical_properties = optical_properties
    ))
  }

  return(par_estimates)
}


#' @export
myFun <- function(x) {
  NA
}

