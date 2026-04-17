#' Retrieve bio-optical parameters via deterministic optimisation
#'
#' Inverts sub-surface \eqn{R_{rs}} against a semi-analytical forward model
#' using gradient-based or derivative-free optimisers.  Returns MAP estimates
#' and approximate standard deviations derived from the Hessian (or the
#' var-cov matrix returned by the Levenberg-Marquardt solver).
#'
#' @param rrs          Data frame with columns \code{wavelength} \[nm\] and
#'   \code{rrs_0m} \[sr\eqn{^{-1}}\].
#' @param forward_model  One of \code{"am03"}, \code{"am03_sicf"}, or
#'   \code{"lee98"}.
#' @param objective_fct  Objective: \code{"log-ll"} (Gaussian log-likelihood),
#'   \code{"rss"} (residual sum of squares), or \code{"lee99"} (spectral error
#'   index).
#' @param optim_mtd  Optimisation algorithm: \code{"L-BFGS-B"},
#'   \code{"Nelder-Mead"}, \code{"levenberg-marqardt"}, or \code{"auglag"}.
#' @param par_inversed  Character vector of parameters to retrieve
#'   (e.g. \code{c("chl", "a_dg_440", "bb_p_550")}).
#' @param par_fixed  Named list of parameters held fixed.  Must not overlap
#'   with \code{par_inversed}.
#' @param lower_b  Named numeric lower bounds.  Required for \code{"L-BFGS-B"};
#'   auto-filled from \code{par_inversed} values otherwise.
#' @param upper_b  Named numeric upper bounds.  Same semantics as \code{lower_b}.
#' @param init_val  Numeric initial values, in the same order as
#'   \code{par_inversed}.
#' @param verbose  Logical.  Print optimiser diagnostics?  Default \code{FALSE}.
#' @param log_prior_fn  Optional \code{function(par)} returning the log-prior
#'   for the inverted parameters.  When supplied the optimiser minimises
#'   \eqn{-(log\_ll + log\_prior)}, i.e. MAP estimation.  \code{NULL} (default)
#'   gives standard MLE.
#' @param spectral_weights  Optional named numeric vector of per-band weights
#'   (see \code{\link{create_spectral_weights}}).
#' @param lat,lon,date_time  Latitude [°N], longitude [°E], and UTC
#'   \code{POSIXct} timestamp.  Required when
#'   \code{forward_model = "am03_sicf"}.
#'
#' @return Named numeric vector of MAP estimates followed by their SDs:
#'   \code{c(chl = 2.1, a_dg_440 = 0.08, ..., chl_sd = 0.4, ...)}.
#'
#' @examples
#' \dontrun{
#' cfg <- make_inversion_params(depth_m = NA, rrs_df = rrs_obs)
#' est <- inverse_gradient(
#'   rrs           = rrs_obs,
#'   forward_model = cfg$fwd_model,
#'   objective_fct = "log-ll",
#'   optim_mtd     = "L-BFGS-B",
#'   par_inversed  = cfg$par_inv,
#'   par_fixed     = cfg$par_fixed,
#'   lower_b       = cfg$lower_b,
#'   upper_b       = cfg$upper_b,
#'   init_val      = cfg$init_val,
#'   log_prior_fn  = cfg$log_prior_fn
#' )
#' }
#'
#' @author Soham Mukherjee, Raphael Mabit
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
    log_prior_fn = NULL,
    spectral_weights = NULL,
    lat = NULL,
    lon = NULL,
    date_time = NULL) {
      
  rlang::inform(paste0("\033[0;33m", "###################################################################", "\033[0m", "\n"))
  rlang::inform(paste0("\033[0;39m", "########### ALL GOOD THINGS ARE WILD & FREE #######", "\033[0m", "\n"))
  rlang::inform(paste0("\033[0;32m", "###################################################################", "\033[0m", "\n"))

  # Auto-add phi_f if using am03_sicf model and phi_f not in par_inversed
  if (forward_model == "am03_sicf" && !"phi_f" %in% par_inversed) {
    warning("'phi_f' not in par_inversed when using am03_sicf model. Adding it automatically.")
    par_inversed <- c(par_inversed, "phi_f")
    lower_b <- c(lower_b, 0.005)
    upper_b <- c(upper_b, 0.03)
    if (!is.null(init_val)) init_val <- c(init_val, 0.02)
  }

  # Inject lat/lon/date_time into par_fixed for SICF forward model
  if (forward_model == "am03_sicf") {
    if (is.null(par_fixed)) par_fixed <- list()
    if (!is.null(lat)       && !"lat"       %in% names(par_fixed)) par_fixed[["lat"]]       <- lat
    if (!is.null(lon)       && !"lon"       %in% names(par_fixed)) par_fixed[["lon"]]       <- lon
    if (!is.null(date_time) && !"date_time" %in% names(par_fixed)) par_fixed[["date_time"]] <- as.numeric(date_time)
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
    model            = forward_model,
    objective        = objective_fct,
    rrs_observed     = rrs,
    par_fixed        = par_fixed,
    par_meta         = par_meta,
    par_inversed     = par_inversed,
    minimize         = TRUE,
    log_prior_fn     = log_prior_fn,
    spectral_weights = spectral_weights
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
    # marqLevAlg has no bounds support. Running in raw parameter space causes:
    #   1. Negative probes → log(negative) = NaN in log_prior_fn / log_ll
    #   2. Parameters spanning 4+ decades → near-singular Hessian → bad SEs
    #
    # Fix: log-reparametrize all positive-definite parameters (OAC + h_w + sd).
    # Benthic fractions (r_rs_b_*, mix_sand) live in [0,1] and stay in natural
    # space — they rarely go negative from a good starting point.
    #
    # Delta method converts log-space SEs back to natural-space SEs:
    #   se_natural[i] = par_natural[i] * se_log[i]   (for log-transformed params)

    is_log_par <- !grepl("^r_rs_b_", par_inversed) & par_inversed != "mix_sand"
    log_idx    <- which(is_log_par)

    # Transform init values: log for positive-definite params
    par_lm          <- par
    par_lm[log_idx] <- log(par[log_idx])

    # Objective wrapper: accepts log-space inputs, calls raw objective with
    # back-transformed (natural-space) values
    minimization_fct_log <- local({
      fn      <- minimization_fct
      log_idx <- log_idx
      function(p) {
        p_nat          <- p
        p_nat[log_idx] <- exp(p[log_idx])
        fn(p_nat)
      }
    })

    lm_result <- marqLevAlg::marqLevAlg(
      b          = par_lm,
      fn         = minimization_fct_log,
      minimize   = TRUE,
      print.info = FALSE
    )

    # Back-transform to natural space
    par_opt          <- lm_result$b
    par_opt[log_idx] <- exp(lm_result$b[log_idx])

    # Reconstruct symmetric var-cov from lower-triangular vector (log-space)
    n_params     <- length(par_lm)
    varcov_log   <- matrix(0, nrow = n_params, ncol = n_params)
    varcov_log[lower.tri(varcov_log, diag = TRUE)] <- lm_result$v
    varcov_log[upper.tri(varcov_log)] <- t(varcov_log)[upper.tri(varcov_log)]

    lm_varcov_is_direct <- lm_result$istop %in% c(1L, 3L)

    optim_result <- list(
      par              = par_opt,
      value            = lm_result$fn.value,
      convergence      = ifelse(lm_result$istop == 1L, 0, 1),
      message          = lm_result$message,
      varcov           = varcov_log,
      varcov_is_direct = lm_varcov_is_direct,
      # stored for delta-method SE conversion in the uncertainty block below
      is_log_par       = is_log_par,
      par_natural      = par_opt
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
    if (optim_mtd == "levenberg-marqardt") {
      # hessian_inverse is the log-space var-cov (or Hessian) from marqLevAlg.
      # Two cases depending on istop:
      #   varcov_is_direct = TRUE  (istop 1): v is already the inverted Hessian
      #   varcov_is_direct = FALSE (istop 2/3/4): v is raw Hessian, must solve()
      # In BOTH cases we must apply the delta method because the optimisation
      # ran in log-space for positive-definite parameters.
      if (isTRUE(optim_result$varcov_is_direct)) {
        varcov_log <- hessian_inverse                       # already var-cov
      } else {
        if (abs(det(hessian_inverse)) < 1e-10) {
          warning("L-M log-space Hessian is singular - using pseudoinverse for SEs")
          varcov_log <- MASS::ginv(hessian_inverse)
        } else {
          varcov_log <- solve(hessian_inverse)
        }
      }
      # Delta method: se_natural[i] = par_natural[i] * se_log[i]
      se_log                          <- sqrt(abs(diag(varcov_log)))
      se_nat                          <- se_log
      se_nat[optim_result$is_log_par] <- optim_result$par_natural[optim_result$is_log_par] *
                                           se_log[optim_result$is_log_par]
      se_nat
    } else {
      # L-BFGS-B / Nelder-Mead: hessian_inverse is the natural-space Hessian of
      # -log L computed by numDeriv::hessian().  Invert to get var-cov.
      if (abs(det(hessian_inverse)) < 1e-5 | abs(det(hessian_inverse)) > 1e10) {
        warning("Hessian is nearly singular - using pseudoinverse")
        varcov <- MASS::ginv(hessian_inverse)
      } else {
        varcov <- solve(hessian_inverse)
      }
      sqrt(abs(diag(varcov)))
    }
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

