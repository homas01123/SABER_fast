#' SABER inverse model using Stan/HMC sampling
#'
#' Uses Hamiltonian Monte Carlo (HMC) with No-U-Turn Sampler (NUTS) via Stan.
#' Much more efficient than random-walk MCMC for complex posteriors, especially
#' in shallow dark waters with low SNR.
#'
#' @param rrs data-frame with columns: wavelength [nm] and rrs_0m [1/sr]
#' @param forward_model Forward model name (default = "am03")
#' @param par_inversed Vector of parameter names to retrieve
#' @param prior List of prior specifications (optional, uses weakly informative defaults)
#' @param lower Numeric vector of lower bounds
#' @param upper Numeric vector of upper bounds
#' @param par_fixed Named list of fixed parameters (water_type, theta_sun, etc.)
#' @param iterations Number of post-warmup iterations (default = 2000)
#' @param warmup Number of warmup/adaptation iterations (default = 1000)
#' @param chains Number of MCMC chains (default = 4)
#' @param adapt_delta Target acceptance rate (0.8-0.99, higher = more accurate, slower)
#' @param max_treedepth Maximum tree depth for NUTS (default = 12)
#' @param backend Which Stan backend to use: "pure_stan" (default), "cpp_optimized", 
#'        or "hybrid". See Details.
#' @param return_fit Return full Stan fit object? (default = FALSE)
#' @param verbose Print Stan sampling progress? (default = TRUE)
#'
#' @return Named vector of MAP estimates and standard deviations.
#'         If return_fit=TRUE, returns list with estimates and fit object.
#'
#' @details
#' **Why HMC/NUTS is better than DEzs/DREAMzs:**
#' 
#' 1. **Gradient-guided exploration**: Uses log-posterior gradients to propose
#'    efficient moves, avoiding random walk behavior.
#'    
#' 2. **Handles correlations**: Explores along posterior contours, not trapped
#'    by parameter correlations (e.g., Chl vs a_g).
#'    
#' 3. **Fewer samples needed**: ~2000-5000 effective samples vs 30,000+ for DEzs.
#' 
#' 4. **Better for shallow waters**: Multi-modal posteriors (different benthic
#'    mixtures) are explored more thoroughly.
#'    
#' 5. **Automatic tuning**: NUTS adapts step size and trajectory length.
#'
#' **When to increase adapt_delta:**
#' - If you see "divergent transitions" warnings
#' - Low SNR / dark waters: use 0.95-0.99
#' - Start with 0.8, increase if needed
#'
#' **Stan Backend Options:**
#' 
#' - `"pure_stan"` (default): All code in Stan language, fully vectorized.
#'   Easy to modify, good performance (~10-20 min for 71 bands).
#'   
#' - `"cpp_optimized"`: Uses external C++ functions (colleague's approach).
#'   Fastest (~3-5 min), requires C++ header compilation.
#'   
#' - `"hybrid"`: C++ for IOP calculation (bottleneck), Stan for forward model.
#'   Balanced speed (~7-12 min) and flexibility. RECOMMENDED for development.
#'
#' @references
#' Hoffman, M. D., & Gelman, A. (2014). The No-U-Turn sampler: adaptively 
#' setting path lengths in Hamiltonian Monte Carlo. JMLR, 15(1), 1593-1623.
#'
#' @export
inverse_mcmc_stan <- function(
    rrs,
    forward_model = "am03",
    par_inversed,
    prior = NULL,
    lower = NULL,
    upper = NULL,
    par_fixed = NULL,
    iterations = 2000,
    warmup = 1000,
    chains = 4,
    adapt_delta = 0.8,
    max_treedepth = 12,
    backend = c("pure_stan", "hybrid", "cpp_optimized"),
    return_fit = FALSE,
    verbose = TRUE
) {
  
  # Match backend argument
  backend <- match.arg(backend)
  
  # Check Stan installation
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "Package 'cmdstanr' required for Stan-based inversion.\n",
      "Install with: install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos')))\n",
      "Then install Stan: cmdstanr::install_cmdstan()"
    )
  }
  
  # Check if CmdStan is installed
  cmdstan_ver <- tryCatch(
    cmdstanr::cmdstan_version(error_on_NA = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(cmdstan_ver) || is.na(cmdstan_ver)) {
    stop(
      "CmdStan not found. Install with: cmdstanr::install_cmdstan()\n",
      "See: https://mc-stan.org/cmdstanr/articles/cmdstanr.html"
    )
  }
  
  # Validate inputs
  if (!all(c("wavelength", "rrs_0m") %in% names(rrs))) {
    stop("rrs must contain columns 'wavelength' and 'rrs_0m'")
  }
  
  # Get wavelength grid
  wavelength <- rrs$wavelength
  
  # Get all cached/interpolated data from C (uses existing cache system)
  # Cache is automatically built if wavelength grid changes
  cached <- .Call("c_get_stan_data", wavelength)
  
  a_w <- cached$a_w
  bb_w <- cached$bb_w
  a0 <- cached$a0
  a1 <- cached$a1
  r_b_matrix <- cached$r_b_matrix
  r_b_names <- cached$r_b_colnames
  
  # Check for shallow vs deep water mode
  has_benthic <- any(grepl("^r_rs_b_|^r_b_", par_inversed)) || 
                 any(grepl("^h_w$", par_inversed))
  
  # Prepare Stan data
  stan_data <- list(
    n_wl = length(wavelength),
    wavelength = wavelength,
    rrs_obs = rrs$rrs_0m,
    sigma = rep(max(mean(rrs$rrs_0m) * 0.03, 1e-5), length(wavelength)),
    
    # Fixed parameters (with defaults)
    water_type = if (is.null(par_fixed$water_type)) 2 else par_fixed$water_type,
    theta_sun_deg = if (is.null(par_fixed$theta_sun)) 30 else par_fixed$theta_sun,
    theta_view_deg = if (is.null(par_fixed$theta_view)) 0 else par_fixed$theta_view,
    
    # Pre-computed tables
    a_w = a_w,
    bb_w = bb_w,
    a0 = a0,
    a1 = a1
  )
  
  # Add shallow water specific data if needed
  if (has_benthic) {
    stan_data$shallow <- 1L
    stan_data$r_b_matrix <- r_b_matrix
    stan_data$n_benthic <- ncol(r_b_matrix)
  }
  
  # Select Stan model based on backend and water depth
  if (has_benthic) {
    stan_filename <- switch(
      backend,
      pure_stan = "am03_shallow.stan",
      cpp_optimized = "am03_cpp_optimized.stan",
      hybrid = "am03_hybrid.stan"
    )
  } else {
    # Deep water - only pure_stan supported currently
    if (backend != "pure_stan") {
      warning("Deep water mode only supports 'pure_stan' backend. Switching to pure_stan.")
      backend <- "pure_stan"
    }
    stan_filename <- "am03_deep.stan"
  }
  
  # Locate Stan model file
  stan_file <- system.file("stan", stan_filename, package = "SABER")
  
  if (!file.exists(stan_file) || stan_file == "") {
    stop(
      "Stan model file not found: ", stan_filename, "\n",
      "Available models: am03, lee98\n",
      "Contact package maintainers if this model is not yet implemented in Stan."
    )
  }
  
  if (verbose) {
    message("Compiling Stan model (", backend, " backend)...")
  }
  
  # Configure compilation based on backend
  if (backend %in% c("cpp_optimized", "hybrid")) {
    # C++ external functions require header include path
    include_path <- system.file("stan/include", package = "SABER")
    
    if (!dir.exists(include_path)) {
      stop(
        "C++ headers not found. This should not happen - please reinstall SABER.\n",
        "Expected: ", include_path
      )
    }
    
    if (verbose) {
      message("  Using C++ external functions from: ", include_path)
    }
    
    # Path to the user header file (required when using --allow-undefined)
    user_header <- file.path(include_path, "rtm_stan_funcs.hpp")
    
    # Compile with C++ headers and allow undefined functions (they're in the header)
    mod <- cmdstanr::cmdstan_model(
      stan_file, 
      include_paths = include_path,
      stanc_options = list("allow-undefined" = TRUE),
      cpp_options = list(
        stan_threads = FALSE,
        USER_HEADER = user_header
      ),
      quiet = !verbose
    )
  } else {
    # Pure Stan - no external dependencies
    mod <- cmdstanr::cmdstan_model(stan_file, quiet = !verbose)
  }
  
  if (verbose) {
    message("\nRunning HMC/NUTS sampling...")
    message(sprintf("  Chains: %d (parallel)", chains))
    message(sprintf("  Warmup: %d iterations", warmup))
    message(sprintf("  Sampling: %d iterations", iterations))
    message(sprintf("  adapt_delta: %.3f", adapt_delta))
    message("")
  }
  
  # Run Stan sampling
  fit <- mod$sample(
    data = stan_data,
    iter_warmup = warmup,
    iter_sampling = iterations,
    chains = chains,
    parallel_chains = chains,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    refresh = if(verbose) 500 else 0,
    show_messages = verbose
  )
  
  # Check for sampling issues
  if (verbose) {
    message("\n", strrep("=", 70))
    message("SAMPLING DIAGNOSTICS")
    message(strrep("=", 70))
    
    fit$diagnostic_summary()
    
    n_divergent <- sum(fit$diagnostic_summary()$num_divergent)
    if (n_divergent > 0) {
      warning(
        sprintf("\n%d divergent transitions detected!", n_divergent),
        "\nIncrease adapt_delta (e.g., 0.95 or 0.99) for better accuracy.",
        "\nSee: https://mc-stan.org/misc/warnings.html#divergent-transitions"
      )
    }
    
    max_treedepth_exceeded <- sum(fit$diagnostic_summary()$num_max_treedepth)
    if (max_treedepth_exceeded > 0) {
      warning(
        sprintf("\n%d iterations exceeded max_treedepth!", max_treedepth_exceeded),
        "\nIncrease max_treedepth (e.g., 15) if you see this warning."
      )
    }
  }
  
  # Extract posterior summary
  draws_summary <- fit$summary(variables = par_inversed)
  
  # Extract MAP estimate (mode of posterior)
  # For Stan, use posterior mean as point estimate (more stable than mode)
  par_mean <- draws_summary$mean
  par_sd <- draws_summary$sd
  
  names(par_mean) <- par_inversed
  names(par_sd) <- paste0(par_inversed, "_sd")
  
  estimates <- c(par_mean, par_sd)
  
  if (verbose) {
    message("\n" , strrep("=", 70))
    message("PARAMETER ESTIMATES (posterior mean ± sd)")
    message(strrep("=", 70))
    for (i in seq_along(par_inversed)) {
      message(sprintf("  %25s: %8.4f ± %8.4f", 
                      par_inversed[i], par_mean[i], par_sd[i]))
    }
    message(strrep("=", 70), "\n")
  }
  
  # Return results
  if (return_fit) {
    return(list(
      estimates = estimates,
      fit = fit,
      summary = draws_summary,
      diagnostics = fit$diagnostic_summary()
    ))
  } else {
    return(estimates)
  }
}
