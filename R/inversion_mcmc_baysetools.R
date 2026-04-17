#' Build a lognormal BayesianTools prior for deep-water inversion parameters
#'
#' Every parameter in par_inversed is given an independent lognormal prior whose
#' log-space mean is the geometric midpoint of [lower_b, upper_b] and whose
#' log-space sd spans half the log range.
#'
#' @param par_inversed character vector of parameter names to invert
#' @param lower_b named numeric vector of lower bounds (must cover all par_inversed)
#' @param upper_b named numeric vector of upper bounds (must cover all par_inversed)
#' @return a BayesianTools prior object
#' @export
make_bt_prior_deep <- function(par_inversed, lower_b, upper_b) {
  lower   <- lower_b[par_inversed]
  upper   <- upper_b[par_inversed]
  log_mid <- 0.5 * (log(lower) + log(upper))
  log_sd  <- 0.5 * (log(upper) - log(lower))
  BayesianTools::createPrior(
    density = function(par) {
      sum(dnorm(log(par), mean = log_mid, sd = log_sd, log = TRUE) - log(par))
    },
    sampler = function(n = 1) {
      mat <- matrix(NA_real_, nrow = n, ncol = length(par_inversed))
      for (i in seq_len(n))
        mat[i, ] <- exp(rnorm(length(par_inversed), mean = log_mid, sd = log_sd))
      mat
    },
    lower = lower,
    upper = upper
  )
}


#' Build a BayesianTools prior for 2-class shallow water (mix_sand parametrisation)
#'
#' OAC parameters receive independent lognormal priors.
#' mix_sand in [0, 1] receives a Beta(2, 2) prior (broad, symmetric, zero-avoiding).
#'
#' @param par_inversed character vector of parameter names to invert
#' @param lower_b named numeric vector of lower bounds
#' @param upper_b named numeric vector of upper bounds
#' @return a BayesianTools prior object
#' @export
make_bt_prior_shallow_2class <- function(par_inversed, lower_b, upper_b) {
  lower   <- lower_b[par_inversed]
  upper   <- upper_b[par_inversed]
  log_mid <- ifelse(par_inversed == "mix_sand", NA_real_, 0.5 * (log(lower) + log(upper)))
  log_sd  <- ifelse(par_inversed == "mix_sand", NA_real_, 0.5 * (log(upper) - log(lower)))
  BayesianTools::createPrior(
    density = function(par) {
      total <- 0
      for (j in seq_along(par)) {
        if (par_inversed[j] == "mix_sand") {
          total <- total + dbeta(par[j], shape1 = 2, shape2 = 2, log = TRUE)
        } else {
          total <- total + dnorm(log(par[j]), mean = log_mid[j], sd = log_sd[j], log = TRUE) - log(par[j])
        }
      }
      total
    },
    sampler = function(n = 1) {
      mat <- matrix(NA_real_, nrow = n, ncol = length(par_inversed))
      for (i in seq_len(n)) {
        for (j in seq_along(par_inversed)) {
          if (par_inversed[j] == "mix_sand") {
            mat[i, j] <- rbeta(1, 2, 2)
          } else {
            mat[i, j] <- exp(rnorm(1, mean = log_mid[j], sd = log_sd[j]))
          }
        }
      }
      mat
    },
    lower = lower,
    upper = upper
  )
}


#' Build a BayesianTools prior for N-class shallow water (soft Dirichlet + lognormal OACs)
#'
#' Parameters whose names start with r_rs_b_ are treated as benthic fractions and
#' receive a soft Dirichlet prior: (alpha-1)*sum(log(b)) - 50*(sum(b)-1)^2.
#' All other parameters receive independent lognormal priors.
#'
#' @param par_inversed character vector of parameter names to invert
#' @param lower_b named numeric vector of lower bounds
#' @param upper_b named numeric vector of upper bounds
#' @param alpha Dirichlet concentration parameter (default 2)
#' @return a BayesianTools prior object
#' @export
make_bt_prior_shallow_nclass <- function(par_inversed, lower_b, upper_b, alpha = 2) {
  is_benthic <- grepl("^r_rs_b_", par_inversed)
  lower   <- lower_b[par_inversed]
  upper   <- upper_b[par_inversed]
  log_mid <- ifelse(is_benthic, NA_real_, 0.5 * (log(lower) + log(upper)))
  log_sd  <- ifelse(is_benthic, NA_real_, 0.5 * (log(upper) - log(lower)))
  BayesianTools::createPrior(
    density = function(par) {
      b_vals <- par[is_benthic]
      oac_lp <- sum(
        dnorm(log(par[!is_benthic]), mean = log_mid[!is_benthic], sd = log_sd[!is_benthic], log = TRUE) -
          log(par[!is_benthic])
      )
      dir_lp <- (alpha - 1) * sum(log(b_vals + 1e-9)) - 50 * (sum(b_vals) - 1)^2
      oac_lp + dir_lp
    },
    sampler = function(n = 1) {
      nb  <- sum(is_benthic)
      mat <- matrix(NA_real_, nrow = n, ncol = length(par_inversed))
      for (i in seq_len(n)) {
        for (j in which(!is_benthic)) {
          repeat {
            v <- exp(rnorm(1, mean = log_mid[j], sd = log_sd[j]))
            if (v >= lower[j] && v <= upper[j]) { mat[i, j] <- v; break }
          }
        }
        g <- rgamma(nb, shape = alpha, rate = 1)
        g <- g / sum(g)
        mat[i, which(is_benthic)] <- g
      }
      mat
    },
    lower = lower,
    upper = upper
  )
}



#' Build an adaptive lognormal BayesianTools prior centred on per-observation spectral estimates
#'
#' Like \code{make_bt_prior_deep} but the log-space mean for each parameter is set
#' from \code{modes} (e.g. OC3 / band-ratio estimates from \code{estimate_prior_modes})
#' rather than the geometric midpoint of the bounds.
#'
#' @param par_inversed character vector of parameter names to invert
#' @param lower_b named numeric vector of lower bounds (must cover all par_inversed)
#' @param upper_b named numeric vector of upper bounds (must cover all par_inversed)
#' @param modes named numeric vector of prior mode estimates (must cover all par_inversed)
#' @return a BayesianTools prior object
#' @export
make_bt_prior_adaptive <- function(par_inversed, lower_b, upper_b, modes) {
  lower   <- lower_b[par_inversed]
  upper   <- upper_b[par_inversed]
  modes_c <- pmax(lower * 1.01, pmin(upper * 0.99, modes[par_inversed]))
  log_mid <- log(modes_c)
  log_sd  <- pmax(0.3,
    pmin(abs(log(upper) - log_mid), abs(log_mid - log(lower))) * 0.5
  )
  BayesianTools::createPrior(
    density = function(par) {
      sum(dnorm(log(par), mean = log_mid, sd = log_sd, log = TRUE) - log(par))
    },
    sampler = function(n = 1) {
      mat <- matrix(NA_real_, nrow = n, ncol = length(par_inversed))
      for (i in seq_len(n))
        for (j in seq_along(par_inversed)) {
          repeat {
            v <- exp(rnorm(1, log_mid[j], log_sd[j]))
            if (v >= lower[j] && v <= upper[j]) { mat[i, j] <- v; break }
          }
        }
      mat
    },
    lower = lower,
    upper = upper
  )
}


#' Build an adaptive 2-class shallow-water BayesianTools prior
#'
#' OAC parameters receive independent lognormal priors centred on \code{modes}.
#' \code{mix_sand} in [0, 1] receives a Beta(2, 2) prior.
#' A \code{repeat\{\}} bounded sampler prevents out-of-bounds start values.
#'
#' @param par_inversed character vector of parameter names to invert
#' @param lower_b named numeric vector of lower bounds
#' @param upper_b named numeric vector of upper bounds
#' @param modes named numeric vector of prior mode estimates (OAC parameters only)
#' @return a BayesianTools prior object
#' @export
make_bt_prior_shallow_2class_adaptive <- function(par_inversed, lower_b, upper_b, modes) {
  mix_idx <- which(par_inversed == "mix_sand")
  oac_idx <- which(par_inversed != "mix_sand")
  lower   <- lower_b[par_inversed]
  upper   <- upper_b[par_inversed]
  modes_c <- pmax(lower[oac_idx] * 1.01, pmin(upper[oac_idx] * 0.99, modes[par_inversed[oac_idx]]))
  log_mid <- log(modes_c)
  log_sd  <- pmax(0.3,
    pmin(abs(log(upper[oac_idx]) - log_mid), abs(log_mid - log(lower[oac_idx]))) * 0.5
  )
  BayesianTools::createPrior(
    density = function(par) {
      oac_lp <- sum(dnorm(log(par[oac_idx]), mean = log_mid, sd = log_sd, log = TRUE) -
                    log(par[oac_idx]))
      mix_lp <- dbeta(par[mix_idx], 2, 2, log = TRUE)
      oac_lp + mix_lp
    },
    sampler = function(n = 1) {
      mat <- matrix(NA_real_, nrow = n, ncol = length(par_inversed))
      for (i in seq_len(n)) {
        for (k in seq_along(oac_idx)) {
          j <- oac_idx[k]
          repeat {
            v <- exp(rnorm(1, log_mid[k], log_sd[k]))
            if (v >= lower[j] && v <= upper[j]) { mat[i, j] <- v; break }
          }
        }
        mat[i, mix_idx] <- max(0.001, min(0.999, rbeta(1, 2, 2)))
      }
      mat
    },
    lower = lower,
    upper = upper
  )
}


#' Retrieve bio-optical parameters via MCMC sampling
#'
#' Inverts sub-surface \eqn{R_{rs}} against a semi-analytical forward model
#' using Bayesian MCMC (via \pkg{BayesianTools}).  Returns MAP estimates and
#' posterior standard deviations.
#'
#' @param rrs  Data frame with columns \code{wavelength} \[nm\] and
#'   \code{rrs_0m} \[sr\eqn{^{-1}}\].
#' @param forward_model  One of \code{"am03"}, \code{"am03_sicf"}, or
#'   \code{"lee98"}.
#' @param par_inversed  Character vector of parameters to retrieve
#'   (e.g. \code{c("chl", "a_dg_440", "bb_p_550")}).
#' @param prior  A BayesianTools prior object (e.g. from
#'   \code{\link{make_bt_prior_deep}}) or \code{NULL} for a uniform prior
#'   defined by \code{lower} / \code{upper}.
#' @param lower  Named numeric lower bounds.  Used only when
#'   \code{prior = NULL}.
#' @param upper  Named numeric upper bounds.  Used only when
#'   \code{prior = NULL}.
#' @param best   Named numeric initial values.  Used only when
#'   \code{prior = NULL}.
#' @param par_fixed  Named list of parameters held fixed during inversion.
#' @param iterations  Total MCMC iterations.  Default \code{10000}
#'   (recommended \eqn{\ge 15000}).
#' @param burnin  Burn-in iterations discarded before posterior summary.
#'   Default \code{2000}.
#' @param sampler  BayesianTools sampler name.  Default \code{"DEzs"}.
#' @param spectral_weights  Optional named numeric vector of per-band weights
#'   (see \code{\link{create_spectral_weights}}).
#' @param lat,lon,date_time  Latitude [°N], longitude [°E], and UTC
#'   \code{POSIXct} timestamp.  Required when
#'   \code{forward_model = "am03_sicf"}.
#'
#' @return Named numeric vector of MAP estimates followed by posterior SDs:
#'   \code{c(chl = 2.5, a_dg_440 = 0.10, ..., chl_sd = 0.3, ...)}.
#'
#' @examples
#' \dontrun{
#' cfg <- make_inversion_params(depth_m = NA, rrs_df = rrs_obs)
#' est <- inverse_mcmc(
#'   rrs           = rrs_obs,
#'   forward_model = cfg$fwd_model,
#'   par_inversed  = cfg$par_inv,
#'   prior         = cfg$bt_prior,
#'   par_fixed     = cfg$par_fixed,
#'   iterations    = 15000,
#'   burnin        = 3000
#' )
#' }
#'
#' @references Mukherjee, S., Mabit, R. and Bélanger, S. (2025).
#'   Limnol. Oceanogr. Methods. \doi{10.1002/lom3.70004}
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
    sampler = "DEzs",
    spectral_weights = NULL,
    lat = NULL,
    lon = NULL,
    date_time = NULL) {

  # Auto-add phi_f if using am03_sicf model and phi_f not in par_inversed
  if (forward_model == "am03_sicf" && !"phi_f" %in% par_inversed) {
    warning("'phi_f' not in par_inversed when using am03_sicf model. Adding it automatically.")
    par_inversed <- c(par_inversed, "phi_f")
    lower <- c(lower, 0.005)
    upper <- c(upper, 0.03)
    if (!is.null(best)) best <- c(best, 0.02)
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

  likelihood <- objective_factory(
    model            = forward_model,
    objective        = "log-ll",
    rrs_observed     = rrs,
    par_inversed     = par_inversed,
    par_fixed        = par_fixed,
    par_meta         = par_meta,
    spectral_weights = spectral_weights
  )

  setup <- BayesianTools::createBayesianSetup(
    prior      = prior,
    likelihood = likelihood,
    lower      = if (is.null(prior)) lower else NULL,
    best       = if (is.null(prior)) best  else NULL,
    upper      = if (is.null(prior)) upper else NULL,
    names      = par_inversed,
    parallel   = FALSE
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
    ~ apply(.x[, 1:(ncol(.x) - 3)], 2, sd)
  )

  estimates_sd <- colMeans(estimates_sd)

  map_values <- BayesianTools::MAP(out)[[1]]

  par_estimates <- stats::setNames(
    c(map_values, estimates_sd),
    c(names(map_values), paste0(names(map_values), "_sd"))
  )

  return(par_estimates)
}
