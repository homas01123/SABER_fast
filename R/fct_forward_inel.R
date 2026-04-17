#'
#' This module provides two approaches for modeling chlorophyll fluorescence:
#'
#' 1. Semi-Analytical Model (Gilerson et al. 2007)
#'    - Fast parameterized approach based on empirical relationships
#'    - Uses Gaussian spectral shapes for emission
#'    - Suitable for operational satellite processing
#'
#' 2. Analytical Model (Full Radiative Transfer with WRF)
#'    - Rigorous wavelength redistribution function (WRF) approach
#'    - Numerical integration over excitation-emission matrix
#'    - Based on Ocean Optics Web Book formulations
#'    - More computationally intensive but physically complete
#'
#' Both models exclude fDOM (CDOM fluorescence) and focus on phytoplankton
#' chlorophyll-a fluorescence. The primary emission is centered at 685 nm
#' with a secondary peak near 730 nm.
#'
#' References:
#' - Gilerson et al. (2007) Opt. Express 15, 15702-15721
#' - Ocean Optics Web Book: Inelastic Scattering
#' - Gregg & Carder (1990) for Ed calculations
#'
#' ============================================================================

# ============================================================================
# Ed Caching System
# ============================================================================
# Create package-level cache for Gregg & Carder Ed calculations
# This prevents redundant computation during inversion iterations where
# Ed parameters (wavelength, sunzen, lat, lon, date_time) remain constant
.Ed_cache <- new.env(hash = TRUE, parent = emptyenv())

#' Compute Ed with Caching
#'
#' Wrapper around Gregg & Carder Ed calculation that caches results.
#' Cache key is based on all input parameters.
#'
#' @param wavelength Vector of wavelengths (nm)
#' @param sunzen_deg Solar zenith angle (degrees)
#' @param lat Latitude (decimal degrees)
#' @param lon Longitude (decimal degrees)
#' @param date_time POSIXct datetime in UTC
#' @param verbose Print cache hit/miss messages
#'
#' @return List with Ed_0m (subsurface Ed) and E0_0m (subsurface scalar irradiance)
#' @keywords internal
.compute_Ed_cached <- function(wavelength, sunzen_deg, lat, lon, date_time, verbose = FALSE) {

  # ---------------------------------------------------------------------------
  # Fast-path cache check: build a lightweight key WITHOUT parsing date or
  # concatenating all wavelengths (both are expensive at 22k MCMC calls/pixel).
  # Wavelength fingerprint = first + last + count: sufficient to distinguish
  # all grids encountered in a single inversion run.
  # ---------------------------------------------------------------------------
  nw     <- length(wavelength)
  dt_num <- as.numeric(date_time)          # cheap — just reads the numeric attribute
  wv_key <- sprintf("%.2f_%.2f_%d", wavelength[1L], wavelength[nw], nw)

  if (sunzen_deg >= 0) {
    cache_key <- paste(wv_key, round(sunzen_deg, 4L), round(lat, 4L),
                       round(lon, 4L), dt_num, sep = "|")
    if (exists(cache_key, envir = .Ed_cache)) {
      if (verbose) message("Ed cache HIT - using cached Gregg & Carder result")
      return(get(cache_key, envir = .Ed_cache))
    }
  }

  # Cache miss (or sunzen was the -99 sentinel) — now do the heavier date parse
  jday_no    <- lubridate::yday(date_time)
  time_no    <- format(date_time, "%T")
  time_parts <- as.numeric(strsplit(time_no, ":")[[1L]])
  time_dec   <- time_parts[1L] + time_parts[2L] / 60

  # Resolve sunzen from lat/lon/time when the sentinel was supplied
  if (sunzen_deg < 0) {
    sunzen_deg <- .gc_sunang(iday = jday_no, hr = time_dec, xlon = lon, ylat = lat)
    cache_key  <- paste(wv_key, round(sunzen_deg, 4L), round(lat, 4L),
                        round(lon, 4L), dt_num, sep = "|")
    if (exists(cache_key, envir = .Ed_cache)) {
      if (verbose) message("Ed cache HIT - using cached Gregg & Carder result")
      return(get(cache_key, envir = .Ed_cache))
    }
  }

  if (verbose) message("Ed cache MISS - computing Gregg & Carder (will be cached)")

  # Calculate Ed at 0+ using internal Gregg & Carder model
  Ed_gc <- .gc_irradiance(
    the     = sunzen_deg,
    lam.sel = wavelength,
    hr      = time_dec,
    jday    = jday_no,
    rlon    = lon,
    rlat    = lat
  )

  # Calculate Fresnel reflectance
  rhoF <- .gc_sfcrfl(theta = sunzen_deg, ws = 5)

  # Convert Ed from 0+ to 0- (subsurface)
  Ed_0m <- (Ed_gc$Edir * (1 - rhoF$rod)) + (Ed_gc$Edif * (1 - rhoF$ros))

  # E0 ≈ Ed × (1 + 1/μ_d);  μ_d = 0.85 for clear sky
  E0_0m <- Ed_0m * (1 + 1 / 0.85)

  result <- list(
    Ed_0m      = Ed_0m,
    E0_0m      = E0_0m,
    sunzen_deg = sunzen_deg
  )

  assign(cache_key, result, envir = .Ed_cache)

  return(result)
}

#' Clear Ed Cache
#'
#' Remove all cached Ed calculations. Useful for testing or memory management.
#' @export
clear_Ed_cache <- function() {
  rm(list = ls(envir = .Ed_cache), envir = .Ed_cache)
  message("Ed cache cleared")
}

#' Get Ed Cache Info
#'
#' Return information about the current cache state.
#' @return List with number of cached entries
#' @export
get_Ed_cache_info <- function() {
  n_entries <- length(ls(envir = .Ed_cache))
  list(
    n_cached = n_entries,
    cache_exists = exists(".Ed_cache", envir = parent.env(environment()))
  )
}

# ============================================================================
# WRF Caching System
# ============================================================================
# Cache pre-computed WRF matrices for different phi_f values
# Allows fast interpolation during MCMC inversion where phi_f varies
.WRF_cache <- new.env(hash = TRUE, parent = emptyenv())

#' Build WRF Cache
#'
#' Pre-computes a single base WRF matrix (phi_f = 1) for the given wavelength
#' grid. The actual quantum yield is applied at compute time inside
#' \code{sicf_analytical()}, exploiting the linearity
#' WRF(i,j; phi_f) = phi_f * WRF_base(i,j).
#'
#' @param wavelength Vector of wavelengths (nm)
#' @param verbose Print progress messages
#'
#' @return Invisible NULL (cache is stored in .WRF_cache environment)
#' @export
build_WRF_cache <- function(wavelength, verbose = TRUE) {

  nwv    <- length(wavelength)
  wv_key <- sprintf("%.4f_%.4f_%d", wavelength[1L], wavelength[nwv], nwv)

  if (verbose) {
    message(sprintf("Building WRF base matrix for %d wavelengths...", nwv))
  }

  start_time <- Sys.time()

  # Build one base matrix with phi_f = 1; actual phi_f is passed at compute time
  WRF_base <- .Call("c_discretize_wrf", as.double(wavelength), 1.0)

  assign("wavelength_key", wv_key,     envir = .WRF_cache)
  assign("wavelength",     wavelength, envir = .WRF_cache)
  assign("WRF_base",       WRF_base,   envir = .WRF_cache)

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  if (verbose) {
    message(sprintf("WRF base matrix built in %.2f seconds (%d x %d)",
                    elapsed, nwv, nwv))
  }

  invisible(NULL)
}

#' Ensure WRF Cache is Initialized
#'
#' Automatically initialize WRF cache if not already built for given wavelengths.
#' Called internally by sicf_analytical() - users don't need to call this directly.
#'
#' @param wavelength Vector of wavelengths (nm)
#' @param verbose Print initialization messages
#'
#' @return Invisible NULL
#' @keywords internal
.ensure_WRF_cache <- function(wavelength, verbose = FALSE) {

  # Check if cache exists
  if (!exists("wavelength_key", envir = .WRF_cache)) {
    if (verbose) {
      message("WRF cache not found - initializing automatically...")
    }
    build_WRF_cache(wavelength, verbose = verbose)
    return(invisible(NULL))
  }

  # Check if wavelength matches
  nwv    <- length(wavelength)
  wv_key <- sprintf("%.4f_%.4f_%d", wavelength[1L], wavelength[nwv], nwv)
  cached_wv_key <- get("wavelength_key", envir = .WRF_cache)

  if (wv_key != cached_wv_key) {
    if (verbose) {
      message("WRF cache wavelength mismatch - rebuilding cache...")
    }
    build_WRF_cache(wavelength, verbose = verbose)
  }

  invisible(NULL)
}

#' Get Cached WRF Base Matrix
#'
#' Returns the pre-computed base WRF matrix (phi_f = 1).  The caller is
#' responsible for scaling by the actual phi_f value.
#'
#' @param wavelength Vector of wavelengths (nm)
#'
#' @return WRF base matrix (N x N, phi_f = 1)
#' @keywords internal
.get_WRF_cached <- function(wavelength) {

  # Ensure cache is initialized (auto-build if needed)
  .ensure_WRF_cache(wavelength, verbose = FALSE)

  get("WRF_base", envir = .WRF_cache)
}

#' Clear WRF Cache
#'
#' Remove all cached WRF matrices. Useful for testing or memory management.
#' @export
clear_WRF_cache <- function() {
  rm(list = ls(envir = .WRF_cache), envir = .WRF_cache)
  message("WRF cache cleared")
}

#' Get WRF Cache Info
#'
#' Return information about the current WRF cache state.
#' @return List with cache details
#' @export
get_WRF_cache_info <- function() {
  if (!exists("wavelength_key", envir = .WRF_cache)) {
    return(list(
      initialized = FALSE,
      n_wavelengths = 0,
      n_phi_f_values = 0,
      phi_f_range = c(NA, NA)
    ))
  }

  wavelength <- get("wavelength", envir = .WRF_cache)

  list(
    initialized = TRUE,
    n_wavelengths = length(wavelength),
    wavelength_range = range(wavelength)
  )
}

# ============================================================================

#' Semi-Analytical SICF Model (Gilerson et al. 2007)
#'
#' Calculates Sun-Induced Chlorophyll Fluorescence using the parameterized
#' approach of Gilerson et al. (2007). This model uses empirical relationships
#' for fluorescence line height at 685 nm combined with Gaussian spectral shapes.
#'
#' @param c_chl Chlorophyll-a concentration (mg m^-3)
#' @param a_dg_443 CDOM+NAP absorption coefficient at 443 nm (m^-1)
#' @param wavelength Vector of wavelengths (nm)
#' @param phi_f Fluorescence quantum yield (dimensionless, typically 0.01-0.02)
#' @param coeff_x Vector of 3 empirical coefficients (default: c(0.0992, 0.40, 0.078))
#' @param Ed_source Either "gregg_carder" for analytical Ed or file path to Ed CSV
#' @param use_analytic_Ed Logical, use Gregg & Carder model? (TRUE/FALSE)
#' @param sunzen_deg Solar zenith angle at surface (degrees). If < 0, calculated from geometry
#' @param lat Latitude (decimal degrees)
#' @param lon Longitude (decimal degrees)
#' @param date_time POSIXct datetime object in UTC
#' @param return_radiance Logical, return Lf instead of Rrs_sicf? (default FALSE)
#'
#' @return Data frame with columns: wavelength, Rrs_sicf (or Lf if return_radiance=TRUE)
#'
#' @references
#' Gilerson, A., Zhou, J., Hlaing, S., Ioannou, I., Schalles, J., Gross, B.,
#' Moshary, F., and Ahmed, S. (2007). "Fluorescence component in the reflectance
#' spectra from coastal waters. Dependence on water composition,"
#' Opt. Express 15, 15702-15721.
#'
#' @export
sicf_semi_analytical <- function(c_chl,
                                  a_dg_443,
                                  wavelength = seq(400, 800, 10),
                                  phi_f = 0.02,
                                  coeff_x = c(0.0992, 0.40, 0.078),
                                  Ed_source = "gregg_carder",
                                  use_analytic_Ed = TRUE,
                                  sunzen_deg = 30,
                                  lat = 49,
                                  lon = -68,
                                  date_time = as.POSIXct("2019-08-18 20:50:00", tz = "UTC"),
                                  return_radiance = FALSE,
                                  verbose = FALSE) {

  # Input validation
  if (c_chl <= 0) stop("c_chl must be positive")
  if (a_dg_443 < 0) stop("a_dg_443 must be non-negative")
  if (phi_f < 0 || phi_f > 0.1) warning("phi_f outside typical range (0.002-0.02)")

  # Calculate Ed at 0- (subsurface)
  if (use_analytic_Ed && Ed_source == "gregg_carder") {
    if (verbose) message("SICF: Using Gregg & Carder (1990) to estimate Ed")

    # Use cached Ed computation
    tryCatch({
      Ed_result <- .compute_Ed_cached(
        wavelength = wavelength,
        sunzen_deg = sunzen_deg,
        lat = lat,
        lon = lon,
        date_time = date_time,
        verbose = verbose
      )

      # Extract results
      Ed_0m_sicf <- Ed_result$Ed_0m
      E0_0m <- Ed_result$E0_0m
      sunzen_deg <- Ed_result$sunzen_deg  # Update if it was calculated

    }, error = function(e) {
      stop(paste("Error in Gregg & Carder Ed calculation:", conditionMessage(e)))
    })

  } else {
    # Load Ed from file
    if (verbose) message("SICF: Loading Ed from file - results may be less accurate")

    if (!file.exists(Ed_source)) {
      stop(paste("Ed file not found:", Ed_source))
    }

    Ed_sim <- read.csv(file = Ed_source, header = TRUE, skip = 9)
    Ed_sim_interp <- Hmisc::approxExtrap(
      x = Ed_sim$Wavelength,
      y = log(Ed_sim$Ed_total.W.m.2.nm.),
      xout = wavelength,
      method = "linear"
    )$y

    Ed_sim_interp <- exp(Ed_sim_interp)  # Extrapolate (NO × 100)
    Ed_sim_interp_0m <- 0.96 * Ed_sim_interp   # Convert to underwater E0-
    Ed_0m_sicf <- Ed_sim_interp_0m
    E0_0m      <- Ed_0m_sicf * (1 + 1 / 0.85)  # approximate scalar irradiance
  }

  # -------------------------------------------------------------------------
  # Core SICF math: use C implementation on the hot inversion path
  # (return_radiance = FALSE) for maximum speed; fall back to R when the
  # caller needs Lf (diagnostic / visualisation use).
  # -------------------------------------------------------------------------
  if (!return_radiance) {

    Rrs_sicf <- .Call(
      "c_sicf_rrs_semi_analytical",
      as.double(wavelength),
      as.double(Ed_0m_sicf),
      as.double(c_chl),
      as.double(a_dg_443),
      as.double(phi_f),
      as.double(coeff_x),   # length-3 vector [c1, c2, c3]
      13.5                  # Lf_685 scale factor
    )

    result <- list(
      Rrs_sicf = Rrs_sicf,
      Ed       = Ed_0m_sicf,
      E0       = E0_0m
    )
    if (verbose) message("Subsurface (0-) Rrs equivalent to SICF calculated")

  } else {

    # R path: kept for diagnostic / return_radiance = TRUE callers
    Lf_685_height <- phi_f * coeff_x[1] * c_chl /
                     (1 + coeff_x[2] * a_dg_443 + coeff_x[3] * c_chl)
    Lf_685 <- Lf_685_height * Ed_0m_sicf[which.min(abs(wavelength - 685))] / 13.5

    d1 <- wavelength - 685
    d2 <- wavelength - 730
    ln2x4 <- 4 * log(2)
    shape <- exp(-ln2x4 * (d1/25)^2) + 0.3 * exp(-ln2x4 * (d2/50)^2)
    Lf <- Lf_685 * shape
    Rrs_sicf <- Lf / Ed_0m_sicf

    result <- list(
      wavelength = wavelength,
      Lf         = Lf,
      Ed         = Ed_0m_sicf,
      E0         = E0_0m,
      Rrs_sicf   = Rrs_sicf
    )
    if (verbose) message("Subsurface (0-) fluorescence radiance (Lf) calculated")

  }

  return(result)
}


#' Wavelength Redistribution Function for Chlorophyll Fluorescence
#'
#' Computes the WRF for chlorophyll-a fluorescence following Ocean Optics
#' Web Book formulations. The WRF describes the probability of a photon
#' absorbed at wavelength λ_ex being re-emitted at wavelength λ_em.
#'
#' @param wave_ex Excitation wavelength (nm)
#' @param wave_em Emission wavelength (nm)
#' @param quant_phi Fluorescence quantum yield (default 0.02)
#'
#' @return WRF value (dimensionless)
#'
#' @details
#' The WRF is defined as:
#' WRF(λ_ex, λ_em) = φ_f × g(λ_ex) × h(λ_em) × λ_ex/λ_em
#'
#' where:
#' - g(λ_ex) = 1 for 370 nm ≤ λ_ex ≤ 690 nm, 0 otherwise (excitation range)
#' - h(λ_em) = Dual Gaussian emission shape:
#'     Primary peak at 685 nm (FWHM = 25 nm)
#'     Secondary peak at 730 nm (FWHM = 50 nm, 30% amplitude)
#' - φ_f = fluorescence quantum yield
#' - λ_ex/λ_em = energy conservation factor
#'
#' @keywords internal
.wrf_chl <- function(wave_ex, wave_em, quant_phi = 0.02) {

  # Primary peak at 685 nm
  wavec0_primary <- 685.0   # Center wavelength of primary emission (nm)
  fwhm_primary <- 25.0      # Full width at half maximum (nm)
  sigmac_primary <- fwhm_primary / (2 * sqrt(2 * log(2)))  # Convert FWHM to σ

  # Secondary peak at 730 nm
  wavec0_secondary <- 730.0  # Center wavelength of secondary emission (nm)
  fwhm_secondary <- 50.0     # Full width at half maximum (nm)
  sigmac_secondary <- fwhm_secondary / (2 * sqrt(2 * log(2)))  # Convert FWHM to σ
  amplitude_secondary <- 0.3  # Secondary peak is 30% of primary

  PhiChl <- quant_phi

  # g(λ_ex): Excitation range function
  # g = 1 for 370 ≤ λ_ex ≤ 690 nm, 0 otherwise
  gchl <- ifelse(wave_ex < 370.0 | wave_ex > 690.0, 0.0, 1.0)

  # h(λ_em): Dual Gaussian emission spectral shape
  # Primary peak (normalized Gaussian)
  factor1_primary <- 1.0 / (sigmac_primary * sqrt(2.0 * pi))
  factor2_primary <- 0.5 / (sigmac_primary * sigmac_primary)
  hchl_primary <- factor1_primary * exp(-factor2_primary * (wave_em - wavec0_primary)^2)

  # Secondary peak (normalized Gaussian, scaled by amplitude)
  factor1_secondary <- 1.0 / (sigmac_secondary * sqrt(2.0 * pi))
  factor2_secondary <- 0.5 / (sigmac_secondary * sigmac_secondary)
  hchl_secondary <- factor1_secondary * exp(-factor2_secondary * (wave_em - wavec0_secondary)^2)

  # Combined emission shape
  hchl <- hchl_primary + amplitude_secondary * hchl_secondary

  # Complete WRF with energy conservation
  wrfchl <- PhiChl * gchl * hchl * wave_ex / wave_em

  return(wrfchl)
}


#' Discretize Wavelength Redistribution Function for SICF
#'
#' Performs numerical integration of the WRF over discretized wavelength bands
#' to create a matrix suitable for computing fluorescence contributions.
#'
#' @param wavelength Vector of wavelength band centers (nm)
#' @param quant_phi Fluorescence quantum yield
#'
#' @return Matrix (Nwave × Nwave) where element [i,j] contains the integrated
#'         WRF for excitation band i contributing to emission band j
#'
#' @details
#' The discretization performs a double integration:
#' WRF_discrete[i,j] = ∫∫ WRF(λ_ex, λ_em) dλ_ex dλ_em
#'
#' over the wavelength ranges defined by bands i and j, with 1 nm step size.
#' Only wavelengths where λ_ex < λ_em contribute (Stokes shift requirement).
#'
#' @keywords internal
.discretize_wrf_chl <- function(wavelength, quant_phi = 0.02) {
  # Fast path: C implementation (50-100x speedup vs nested R loops)
  return(.Call("c_discretize_wrf", as.double(wavelength), as.double(quant_phi)))

  deltaw <- 1  # Integration step size (nm)  # nolint (dead code — R fallback preserved for reference)
  waveb <- wavelength
  Nwave <- length(waveb)

  # Initialize WRF matrix
  WRF_Chl <- matrix(0, nrow = Nwave, ncol = Nwave)

  # Pre-compute band widths (vectorized)
  band_widths <- c(deltaw, diff(waveb))

  # Pre-compute number of integration steps
  nsteps <- pmax(1, as.integer(band_widths / deltaw))

  # Pre-compute band centers
  wave_centers <- waveb + 0.5 * deltaw

  # Pre-compute factors
  factors <- deltaw * deltaw / band_widths

  # Loop over emission bands (j) - cannot fully vectorize outer loop
  for (j in 2:Nwave) {
    nj <- nsteps[j]

    # Create emission wavelength grid for this band (vectorized)
    wavej_grid <- wave_centers[j] + seq(0, nj - 1) * deltaw

    # Loop over excitation bands (i < j, Stokes shift)
    # Pre-compute for all i at once
    i_range <- 1:(j - 1)
    ni_vec <- nsteps[i_range]

    for (idx in seq_along(i_range)) {
      i <- i_range[idx]
      ni <- ni_vec[idx]

      # Create excitation wavelength grid for this band (vectorized)
      wavei_grid <- wave_centers[i] + seq(0, ni - 1) * deltaw

      # Vectorized double integration using outer product
      # Create all combinations of wavei and wavej
      wavei_matrix <- matrix(rep(wavei_grid, each = nj), nrow = ni * nj)
      wavej_matrix <- matrix(rep(wavej_grid, times = ni), nrow = ni * nj)

      # Compute WRF for all combinations at once (fully vectorized!)
      wrf_values <- .wrf_chl(wavei_matrix, wavej_matrix, quant_phi = quant_phi)

      # Sum over all integration points
      sumchl <- sum(wrf_values)

      WRF_Chl[i, j] <- factors[j] * sumchl
    }
  }

  return(WRF_Chl)
}

#' Discretize WRF - Fully Vectorized Alternative (Advanced)
#'
#' This is an experimental fully-vectorized version that pre-computes
#' ALL wavelength combinations at once. Uses more memory but much faster.
#'
#' @keywords internal
.discretize_wrf_chl_vectorized <- function(wavelength, quant_phi = 0.02) {

  deltaw <- 1
  waveb <- wavelength
  Nwave <- length(waveb)

  # Pre-compute all band properties (vectorized)
  band_widths <- c(deltaw, diff(waveb))
  nsteps <- pmax(1, as.integer(band_widths / deltaw))
  wave_centers <- waveb + 0.5 * deltaw
  factors <- deltaw * deltaw / band_widths

  # Pre-compute maximum grid size needed
  max_ni <- max(nsteps)
  max_nj <- max(nsteps)

  # Create master excitation and emission grids
  # Use expand.grid-like approach but more memory efficient
  excitation_grids <- lapply(seq_along(waveb), function(i) {
    wave_centers[i] + seq(0, nsteps[i] - 1) * deltaw
  })

  emission_grids <- lapply(seq_along(waveb), function(j) {
    wave_centers[j] + seq(0, nsteps[j] - 1) * deltaw
  })

  # Initialize WRF matrix
  WRF_Chl <- matrix(0, nrow = Nwave, ncol = Nwave)

  # Compute WRF for all valid (i < j) pairs
  for (j in 2:Nwave) {
    wavej_vec <- emission_grids[[j]]
    nj <- length(wavej_vec)

    # Vectorize over all i < j at once
    for (i in 1:(j - 1)) {
      wavei_vec <- excitation_grids[[i]]
      ni <- length(wavei_vec)

      # Use outer product to create all combinations efficiently
      # This is KEY: outer() is implemented in C and very fast
      wrf_matrix <- outer(wavei_vec, wavej_vec, function(ex, em) {
        .wrf_chl(ex, em, quant_phi)
      })

      # Sum all contributions
      WRF_Chl[i, j] <- factors[j] * sum(wrf_matrix)
    }
  }

  return(WRF_Chl)
}


#' Calculate Phytoplankton Absorption Spectrum
#'
#' Estimates phytoplankton absorption coefficient from chlorophyll concentration
#' using specific absorption coefficient and package effect.
#'
#' @param c_chl Chlorophyll-a concentration (mg m^-3)
#' @param wavelength Vector of wavelengths (nm)
#'
#' @return Vector of phytoplankton absorption coefficients (m^-1)
#'
#' @details
#' Uses the relationship: a_phy(λ) = a*_phy(λ) × Chl
#' where a*_phy is the chlorophyll-specific absorption coefficient.
#'
#' For simplicity, uses typical values from NOMAD or similar dataset.
#' For more accuracy, load actual a*_phy(λ) data.
#'
#' @keywords internal
.calculate_a_phy <- function(c_chl, wavelength) {

  # Simplified phytoplankton-specific absorption
  # These are approximate values - for production use, load actual a*_phy data
  # from NOMAD or similar datasets

  # Simple parameterization based on Bricaud et al. (1995, 1998, 2004)
  # a_phy(λ) = A(λ) × Chl^B(λ)

  # Approximate coefficients at key wavelengths
  # For full implementation, use lookup table with interpolation

  # Very simplified for demonstration - peak absorption near 440 and 675 nm
  a_phy <- rep(0, length(wavelength))

  for (i in seq_along(wavelength)) {
    wv <- wavelength[i]

    # Rough spectral shape - peaks at blue (440) and red (675)
    if (wv < 500) {
      # Blue peak
      a_star <- 0.05 * exp(-((wv - 440)/30)^2) + 0.01
    } else if (wv > 600) {
      # Red peak
      a_star <- 0.03 * exp(-((wv - 675)/25)^2) + 0.005
    } else {
      # Green minimum
      a_star <- 0.008
    }

    # Package effect: a_phy = a*_phy × Chl^0.65
    a_phy[i] <- a_star * c_chl^0.65
  }

  return(a_phy)
}

# ============================================================================
# Internal Function: Gaussian Vertical Chlorophyll Profile
# ============================================================================
#' @keywords internal
.gaussian_chl_profile <- function(chl_surface, depths, z_max = 15, sigma = 8) {
  # Realistic subsurface chlorophyll maximum (SCM)
  # Surface value = chl_surface, peak at z_max with ~2x enhancement
  peak_enhancement <- 2.0
  chl_z <- chl_surface + (chl_surface * (peak_enhancement - 1)) * exp(-0.5 * ((depths - z_max) / sigma)^2)
  return(chl_z)
}

#' Plot Depth Profiles of E0 and Fluorescence Attenuation
#'
#' Creates diagnostic plots showing how scalar irradiance and fluorescence
#' radiance decay with depth for different wavelengths.
#'
#' @param wavelength Vector of wavelengths (nm)
#' @param E0_surface Scalar irradiance at surface (W m^-2 nm^-1)
#' @param Kd Diffuse attenuation coefficients (m^-1)
#' @param a_phy Phytoplankton absorption (m^-1)
#' @param WRF_matrix Wavelength redistribution function matrix
#' @param phi_f Fluorescence quantum yield
#' @param c_chl Chlorophyll concentration for plot title
#' @param a_dg_443 CDOM absorption for plot title
#' @param bb_p_550 Particulate backscattering for plot title
#'
#' @return List with E0_z and Lf_z matrices, plus ggplot objects
#'
#' @keywords internal
.plot_depth_profiles_sicf <- function(wavelength, E0_surface, Kd, a_phy,
                                      WRF_matrix, phi_f, c_chl, a_dg_443, bb_p_550) {

  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("scales", quietly = TRUE)) {
    warning("ggplot2 and scales required for depth profiles")
    return(NULL)
  }

  library(ggplot2)
  library(scales)

  # Depth grid (0 to 25 m, 0.5 m steps - reduced to avoid extremely small Lf values)
  depths <- seq(0, 25, by = 0.5)
  n_depths <- length(depths)
  n_wave <- length(wavelength)

  # Generate Gaussian chlorophyll profile (SCM at 15m with 2x surface enhancement)
  chl_z <- .gaussian_chl_profile(c_chl, depths, z_max = 15, sigma = 8)

  # Calculate E0(λ, z) = E0(λ, 0-) × exp(-Kd × z)
  E0_z <- matrix(0, nrow = n_depths, ncol = n_wave)
  for (i in 1:n_wave) {
    E0_z[, i] <- E0_surface[i] * exp(-Kd[i] * depths)
  }

  # Calculate Lf(λ_em, z) for each depth
  # Note: Use surface a_phy (from input parameter) for all depths to match original behavior
  # The depth variation comes from E0 attenuation, not from changing absorption
  Lf_z <- matrix(0, nrow = n_depths, ncol = n_wave)

  for (z_idx in 1:n_depths) {
    for (j in 1:n_wave) {
      sum_fluor <- 0
      for (i in 1:n_wave) {
        if (WRF_matrix[i, j] > 0) {
          flux_absorbed <- E0_z[z_idx, i] * a_phy[i]  # Use surface a_phy, not depth-varying
          flux_emitted <- flux_absorbed * WRF_matrix[i, j]
          sum_fluor <- sum_fluor + flux_emitted
        }
      }
      Lf_z[z_idx, j] <- sum_fluor / (4 * pi)
    }
  }

  # Check for Lf values
  Lf_max <- max(Lf_z, na.rm = TRUE)
  if (!is.finite(Lf_max) || Lf_max <= 0 || Lf_max < 1e-20) {
    warning("All Lf values are zero or extremely small (<1e-20)! Cannot create meaningful Lf plot.")

    # Return E0 plot only, with NULL for Lf plot
    # Still create E0 plot...
  }

  # Select representative wavelengths for plotting
  wv_plot <- c(440, 550, 685, 730)
  wv_indices <- sapply(wv_plot, function(w) which.min(abs(wavelength - w)))
  wv_colors <- c("440" = "#0000FF", "550" = "#00AA00",
                 "685" = "#FF0000", "730" = "#AA0000")

  # Calculate euphotic depths (90% light attenuation = 2.3/Kd)
  euphotic_depths <- 2.3 / Kd[wv_indices]

  # Create data frame for E0 profiles
  E0_df <- data.frame()
  for (i in seq_along(wv_plot)) {
    E0_df <- rbind(E0_df, data.frame(
      depth = depths,
      E0 = E0_z[, wv_indices[i]],
      chl = chl_z,
      wavelength = factor(wv_plot[i]),
      euphotic = euphotic_depths[i]
    ))
  }

  # Create data frame for Lf profiles (NO SCALING - use actual values)
  Lf_max_val <- max(Lf_z, na.rm = TRUE)
  Lf_min_val <- min(Lf_z[Lf_z > 0], na.rm = TRUE)

  Lf_df <- data.frame()
  for (i in seq_along(wv_plot)) {
    Lf_df <- rbind(Lf_df, data.frame(
      depth = depths,
      Lf = Lf_z[, wv_indices[i]],  # Use actual Lf values, no scaling
      chl = chl_z,
      wavelength = factor(wv_plot[i]),
      euphotic = euphotic_depths[i]
    ))
  }

  # Plot E0(z) with chlorophyll on top axis
  # Create an ULTRA-HIGH-RESOLUTION dataframe for filled curve appearance
  depths_fine <- seq(0, 25, by = 0.001)  # 0.001m = 1mm resolution for dense filled curve
  chl_z_fine <- .gaussian_chl_profile(c_chl, depths_fine, z_max = 15, sigma = 8)
  chl_df <- data.frame(depth = depths_fine, chl = chl_z_fine)

  chl_range <- range(chl_z_fine)
  E0_range_log <- range(log10(E0_df$E0[E0_df$E0 > 0]), na.rm = TRUE)

  # Check if chlorophyll varies enough to plot secondary axis
  has_chl_variation <- diff(chl_range) > 0.01 * chl_range[1]

  if (has_chl_variation) {
    # Transform Chl to E0 scale for ultra-dense overlay line
    chl_to_E0_slope <- diff(E0_range_log) / diff(chl_range)
    chl_to_E0_intercept <- E0_range_log[1] - chl_to_E0_slope * chl_range[1]
    chl_df$E0_position <- 10^(chl_to_E0_intercept + chl_to_E0_slope * chl_df$chl)

    p_E0 <- ggplot(E0_df, aes(y = depth)) +
      geom_line(aes(x = E0, color = wavelength), size = 1.3) +
      geom_line(data = chl_df, aes(x = E0_position, y = depth),
                color = "grey50", size = 0.8, alpha = 0.6) +
      geom_hline(data = E0_df %>% dplyr::distinct(wavelength, euphotic),
                 aes(yintercept = euphotic, color = wavelength),
                 linetype = "dashed", size = 1.3, alpha = 0.6) +
      scale_y_reverse(limits = c(25, 0), breaks = seq(0, 25, 5)) +
      scale_x_log10(name = expression(paste(E[0], "(", lambda, ", z) [W m"^{-2}, " nm"^{-1}, "]")),
                    labels = trans_format("log10", math_format(10^.x)),
                    sec.axis = sec_axis(~ chl_range[1] + (log10(.) - E0_range_log[1]) /
                                          diff(E0_range_log) * diff(chl_range),
                                        name = expression(paste("Chl [mg m"^{-3}, "]")))) +
      scale_color_manual(name = expression(paste(lambda, " [nm]")), values = wv_colors) +
      theme_bw() +
      theme(axis.text.x.bottom = element_text(size = 20, color = 'black'),
            axis.text.x.top = element_text(size = 18, color = 'grey30'),
            axis.text.y = element_text(size = 20, color = 'black'),
            axis.title.x.bottom = element_text(size = 25),
            axis.title.x.top = element_text(size = 22, color = 'grey30'),
            axis.title.y = element_text(size = 25),
            axis.ticks.length = unit(.25, "cm"),
            legend.position = c(0.15, 0.25),
            legend.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size = 18),
            legend.background = element_rect(fill = NA, size = 0.5, linetype = "solid", colour = 0),
            legend.key = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_line(colour = "black", size = 0.5, linetype = "dotted"),
            panel.grid.minor = element_line(colour = "grey80", linewidth = 0.2, linetype = "solid"),
            plot.margin = unit(c(0.5, 1.0, 0.5, 0.5), "cm"),
            panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) +
      labs(y = "Depth [m]")
  } else {
    # No Chl variation - skip secondary axis
    p_E0 <- ggplot(E0_df, aes(y = depth)) +
      geom_line(aes(x = E0, color = wavelength), size = 1.3) +
      geom_hline(data = E0_df %>% dplyr::distinct(wavelength, euphotic),
                 aes(yintercept = euphotic, color = wavelength),
                 linetype = "dashed", size = 1.3, alpha = 0.6) +
      scale_y_reverse(limits = c(25, 0), breaks = seq(0, 25, 5)) +
      scale_x_log10(name = expression(paste(E[0], "(", lambda, ", z) [W m"^{-2}, " nm"^{-1}, "]")),
                    labels = trans_format("log10", math_format(10^.x))) +
      scale_color_manual(name = expression(paste(lambda, " [nm]")), values = wv_colors) +
      theme_bw() +
      theme(axis.text.x.bottom = element_text(size = 20, color = 'black'),
            axis.text.y = element_text(size = 20, color = 'black'),
            axis.title.x.bottom = element_text(size = 25),
            axis.title.y = element_text(size = 25),
            axis.ticks.length = unit(.25, "cm"),
            legend.position = c(0.15, 0.25),
            legend.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size = 18),
            legend.background = element_rect(fill = NA, size = 0.5, linetype = "solid", colour = 0),
            legend.key = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_line(colour = "black", size = 0.5, linetype = "dotted"),
            panel.grid.minor = element_line(colour = "grey80", linewidth = 0.2, linetype = "solid"),
            plot.margin = unit(c(0.5, 1.0, 0.5, 0.5), "cm"),
            panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) +
      labs(y = "Depth [m]")
  }

  # Plot Lf(z) with chlorophyll on top axis - use actual values
  Lf_label <- expression(paste(L[f], "(", lambda, ", z) [W m"^{-2}, " sr"^{-1}, " nm"^{-1}, "]"))

  Lf_range_log <- suppressWarnings(range(log10(Lf_df$Lf[Lf_df$Lf > 0]), na.rm = TRUE))

  # Try log scale first (should work with realistic Lf values)
  if (is.finite(Lf_range_log[1]) && is.finite(Lf_range_log[2]) && diff(Lf_range_log) > 0.5) {
    # Log scale works - use it with chlorophyll secondary axis

    # Check if chlorophyll varies enough for secondary axis
    if (diff(chl_range) > 0.01 * chl_range[1]) {
      # Transform Chl to Lf scale for ultra-dense overlay line
      chl_to_Lf_slope <- diff(Lf_range_log) / diff(chl_range)
      chl_to_Lf_intercept <- Lf_range_log[1] - chl_to_Lf_slope * chl_range[1]
      chl_df$Lf_position <- 10^(chl_to_Lf_intercept + chl_to_Lf_slope * chl_df$chl)

      p_Lf <- ggplot(Lf_df, aes(y = depth)) +
        geom_line(aes(x = Lf, color = wavelength), size = 1.3) +
        geom_line(data = chl_df, aes(x = Lf_position, y = depth),
                  color = "grey50", size = 0.8, alpha = 0.6) +
        geom_hline(data = Lf_df %>% dplyr::distinct(wavelength, euphotic),
                   aes(yintercept = euphotic, color = wavelength),
                   linetype = "dashed", size = 1.3, alpha = 0.6) +
        scale_y_reverse(limits = c(25, 0), breaks = seq(0, 25, 5)) +
        scale_x_log10(name = Lf_label,
                      labels = trans_format("log10", math_format(10^.x)),
                      sec.axis = sec_axis(~ chl_range[1] + (log10(.) - Lf_range_log[1]) /
                                            diff(Lf_range_log) * diff(chl_range),
                                          name = expression(paste("Chl [mg m"^{-3}, "]")))) +
        scale_color_manual(name = expression(paste(lambda, " [nm]")), values = wv_colors) +
        theme_bw() +
        theme(axis.text.x.bottom = element_text(size = 20, color = 'black'),
              axis.text.x.top = element_text(size = 18, color = 'grey30'),
              axis.text.y = element_text(size = 20, color = 'black'),
              axis.title.x.bottom = element_text(size = 25),
              axis.title.x.top = element_text(size = 22, color = 'grey30'),
              axis.title.y = element_text(size = 25),
              axis.ticks.length = unit(.25, "cm"),
              legend.position = c(0.15, 0.25),
              legend.title = element_text(size = 20, face = "bold"),
              legend.text = element_text(size = 18),
              legend.background = element_rect(fill = NA, size = 0.5, linetype = "solid", colour = 0),
              legend.key = element_blank(),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = "black", size = 0.5, linetype = "dotted"),
              panel.grid.minor = element_line(colour = "grey80", linewidth = 0.2, linetype = "solid"),
              plot.margin = unit(c(0.5, 1.0, 0.5, 0.5), "cm"),
              panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) +
        labs(y = "Depth [m]")
    } else {
      # No Chl variation - skip secondary axis
      p_Lf <- ggplot(Lf_df, aes(y = depth)) +
        geom_line(aes(x = Lf, color = wavelength), size = 1.3) +
        geom_hline(data = Lf_df %>% dplyr::distinct(wavelength, euphotic),
                   aes(yintercept = euphotic, color = wavelength),
                   linetype = "dashed", size = 1.3, alpha = 0.6) +
        scale_y_reverse(limits = c(25, 0), breaks = seq(0, 25, 5)) +
        scale_x_log10(name = Lf_label,
                      labels = trans_format("log10", math_format(10^.x))) +
        scale_color_manual(name = expression(paste(lambda, " [nm]")), values = wv_colors) +
        theme_bw() +
        theme(axis.text.x.bottom = element_text(size = 20, color = 'black'),
              axis.text.y = element_text(size = 20, color = 'black'),
              axis.title.x.bottom = element_text(size = 25),
              axis.title.y = element_text(size = 25),
              axis.ticks.length = unit(.25, "cm"),
              legend.position = c(0.15, 0.25),
              legend.title = element_text(size = 20, face = "bold"),
              legend.text = element_text(size = 18),
              legend.background = element_rect(fill = NA, size = 0.5, linetype = "solid", colour = 0),
              legend.key = element_blank(),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = "black", size = 0.5, linetype = "dotted"),
              panel.grid.minor = element_line(colour = "grey80", linewidth = 0.2, linetype = "solid"),
              plot.margin = unit(c(0.5, 1.0, 0.5, 0.5), "cm"),
              panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) +
        labs(y = "Depth [m]")
    }
  } else {
    # Fallback to linear scale
    warning("Using linear scale for Lf plot (log scale not feasible).")
    Lf_range <- range(Lf_df$Lf, na.rm = TRUE)

    p_Lf <- ggplot(Lf_df, aes(y = depth)) +
      geom_line(aes(x = Lf, color = wavelength), size = 1.3) +
      geom_hline(data = Lf_df %>% dplyr::distinct(wavelength, euphotic),
                 aes(yintercept = euphotic, color = wavelength),
                 linetype = "dashed", size = 1.3, alpha = 0.6) +
      scale_y_reverse(limits = c(25, 0), breaks = seq(0, 25, 5)) +
      scale_x_continuous(name = Lf_label,
                         labels = function(x) format(x, scientific = TRUE, digits = 3)) +
      scale_color_manual(name = expression(paste(lambda, " [nm]")), values = wv_colors) +
      theme_bw() +
      theme(axis.text.x.bottom = element_text(size = 20, color = 'black'),
            axis.text.y = element_text(size = 20, color = 'black'),
            axis.title.x.bottom = element_text(size = 25),
            axis.title.y = element_text(size = 25),
            axis.ticks.length = unit(.25, "cm"),
            legend.position = c(0.15, 0.25),
            legend.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size = 18),
            legend.background = element_rect(fill = NA, size = 0.5, linetype = "solid", colour = 0),
            legend.key = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_line(colour = "black", size = 0.5, linetype = "dotted"),
            panel.grid.minor = element_line(colour = "grey80", linewidth = 0.2, linetype = "solid"),
            plot.margin = unit(c(0.5, 1.0, 0.5, 0.5), "cm"),
            panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) +
      labs(y = "Depth [m]")
  }
  return(list(
    E0_z = E0_z,
    Lf_z = Lf_z,
    depths = depths,
    chl_z = chl_z,
    plot_E0 = p_E0,
    plot_Lf = p_Lf,
    euphotic_depths = euphotic_depths
  ))
}


#' Analytical SICF Model (Full Radiative Transfer)
#'
#' Calculates SICF using complete wavelength redistribution approach with
#' numerical integration over excitation-emission matrix. This is the most
#' rigorous model based on fundamental radiative transfer theory.
#'
#' @param c_chl Chlorophyll-a concentration (mg m^-3)
#' @param a_phy Phytoplankton absorption spectrum (m^-1). If NULL, calculated from c_chl
#' @param wavelength Vector of wavelengths (nm)
#' @param phi_f Fluorescence quantum yield (default 0.02)
#' @param Ed_source Either "gregg_carder" for analytical Ed or file path to Ed CSV
#' @param use_analytic_Ed Use Gregg & Carder model? (TRUE/FALSE)
#' @param sunzen_deg Solar zenith angle (degrees)
#' @param lat Latitude (decimal degrees)
#' @param lon Longitude (decimal degrees)
#' @param date_time POSIXct datetime in UTC
#' @param scalar_irradiance Use scalar irradiance E0 instead of Ed? (default TRUE)
#' @param depth_resolved Calculate depth-integrated fluorescence with attenuation? (default FALSE)
#' @param a_dg_443 CDOM+NAP absorption at 443 nm (m^-1), required if depth_resolved=TRUE
#' @param bb_p_550 Particulate backscattering at 550 nm (m^-1), required if depth_resolved=TRUE
#' @param plot_depth_profiles Create diagnostic plots of E0(z) and Lf(z)? (default FALSE)
#'
#' @return Data frame with columns: wavelength, Rrs_sicf. If plot_depth_profiles=TRUE,
#'         also returns depth_profiles list with E0_z and Lf_z matrices.
#'
#' @details
#' **Surface-only mode (depth_resolved = FALSE)**:
#' Lf(λ_em) = (1/4π) × Σ[E0(λ_ex) × a_phy(λ_ex) × WRF(λ_ex, λ_em)]
#'
#' Assumes all fluorescence occurs at surface (z=0-).
#'
#' **Depth-integrated mode (depth_resolved = TRUE)**:
#' Lf(λ_em, z=0-) = (1/4π) × Σ[E0(λ_ex, 0-) × a_phy(λ_ex) × WRF(λ_ex, λ_em) / (Kd_ex + Kd_em)]
#'
#' where Kd = (a_total + bb_total) / μ_d with μ_d ≈ 0.8.
#' The factor 1/(Kd_ex + Kd_em) represents depth integration of exponentially
#' attenuated excitation and emission: ∫[0→∞] exp(-(Kd_ex + Kd_em)×z) dz
#'
#' Requires a_dg_443 and bb_p_550 to calculate total IOPs using iop_from_oac().
#'
#' @export
sicf_analytical <- function(c_chl,
                             a_phy = NULL,
                             wavelength = seq(400, 800, 10),
                             phi_f = 0.02,
                             Ed_source = "gregg_carder",
                             use_analytic_Ed = TRUE,
                             sunzen_deg = 30,
                             lat = 49,
                             lon = -68,
                             date_time = as.POSIXct("2019-08-18 20:50:00", tz = "UTC"),
                             scalar_irradiance = TRUE,
                             depth_resolved = FALSE,
                             a_dg_443 = NULL,
                             bb_p_550 = NULL,
                             plot_depth_profiles = FALSE,
                             verbose = FALSE) {

  # Input validation
  if (c_chl <= 0) stop("c_chl must be positive")
  if (phi_f < 0 || phi_f > 0.1) warning("phi_f outside typical range (0.002-0.02)")

  # Check depth_resolved requirements
  if (depth_resolved) {
    if (is.null(a_dg_443)) stop("a_dg_443 required when depth_resolved=TRUE")
    if (is.null(bb_p_550)) stop("bb_p_550 required when depth_resolved=TRUE")
    if (verbose) message("SICF Analytical: DEPTH-RESOLVED mode - calculating Kd and depth integration")
  } else {
    if (verbose) message("SICF Analytical: SURFACE-ONLY mode - computing full wavelength redistribution model")
  }

  # a_phy: validate if supplied externally; otherwise defer to iop_from_oac
  # (Bricaud model) below — used for both surface-only and depth-resolved paths.
  if (!is.null(a_phy) && length(a_phy) != length(wavelength)) {
    stop("a_phy must have same length as wavelength")
  }

  # Calculate scalar irradiance E0 at 0-
  if (use_analytic_Ed && Ed_source == "gregg_carder") {
    if (verbose) message("SICF Analytical: Using Gregg & Carder (1990) for E0 calculation")

    # Use cached Ed computation
    Ed_result <- .compute_Ed_cached(
      wavelength = wavelength,
      sunzen_deg = sunzen_deg,
      lat = lat,
      lon = lon,
      date_time = date_time,
      verbose = verbose
    )

    # Extract results
    Ed_0m <- Ed_result$Ed_0m
    E0_0m <- Ed_result$E0_0m
    sunzen_deg <- Ed_result$sunzen_deg  # Update if it was calculated

    # Store Ed for Rrs calculation (same as semi-analytical)
    Ed_for_rrs <- Ed_0m

    # Convert to scalar irradiance if requested
    if (scalar_irradiance) {
      # E0 = Ed / cos(θ_w) where θ_w is refracted sun angle
      sun_view <- snell_law(theta_view = 0, theta_sun = sunzen_deg)
      E0_0m_corrected <- Ed_0m / cos(sun_view$sun_w)
      irrad <- E0_0m_corrected  # NO × 100
      if (verbose) message("SICF Analytical: Using scalar irradiance E0")
    } else {
      irrad <- Ed_0m  # NO × 100
      if (verbose) message("SICF Analytical: Using downwelling irradiance Ed")
    }

  } else {
    # Load from file
    if (verbose) message("SICF Analytical: Loading Ed from file")

    if (!file.exists(Ed_source)) {
      stop(paste("Ed file not found:", Ed_source))
    }

    Ed_sim <- read.csv(file = Ed_source, header = TRUE, skip = 9)
    Ed_sim_interp <- Hmisc::approxExtrap(
      x = Ed_sim$Wavelength,
      y = log(Ed_sim$Ed_total.W.m.2.nm.),
      xout = wavelength,
      method = "linear"
    )$y

    irrad <- exp(Ed_sim_interp)  # NO × 100

    if (scalar_irradiance) {
      # Approximate conversion Ed -> E0
      # E0 ≈ Ed / (2 × cos(θ))
      # Use simplified conversion
      irrad <- irrad / (2 * cos(sunzen_deg * pi/180))
      if (verbose) message("SICF Analytical: Converted Ed to approximate E0")
    }

    irrad <- 0.96 * irrad  # Surface transmission
    Ed_for_rrs <- irrad  # Store for Rrs calculation
  }

  # Discretize WRF into matrix (with automatic caching for performance)
  if (verbose) message("SICF Analytical: Retrieving wavelength redistribution function...")

  # Use cached WRF base matrix (phi_f = 1); actual phi_f applied in C
  WRF_matrix <- .get_WRF_cached(wavelength)

  # Calculate Kd for depth-resolved mode
  Kd_wavelength <- NULL
  if (depth_resolved) {
    if (verbose) {
      message("SICF Analytical: Calculating diffuse attenuation coefficients (Kd)...")
    }

    # Build parameter vector for iop_from_oac
    par <- c(
      "chl" = c_chl,
      "a_dg_440" = a_dg_443,  # Will be adjusted to 440 nm inside iop_from_oac
      "bb_p_550" = bb_p_550
    )

    # Get total IOPs (a_w + a_phy + a_g, bb_w + bb_p) and extract a_phy
    iop_total <- iop_from_oac(wavelength = wavelength, par = par)

    # Reuse Bricaud a_phy from iop_from_oac — consistent with Kd, no extra loop
    a_phy <- iop_total$a_phy

    # Calculate Kd = (a + bb) / μ_d
    # μ_d ≈ 0.8 is average cosine for diffuse light field
    mu_d <- 0.8
    Kd_wavelength <- (iop_total$a + iop_total$bb) / mu_d

    if (verbose) {
      message(sprintf("SICF Analytical: Kd range = %.4f to %.4f m⁻¹",
                      min(Kd_wavelength), max(Kd_wavelength)))
    }
  } else if (is.null(a_phy)) {
    # Surface-only mode: get Bricaud a_phy from iop_from_oac (chl only)
    if (verbose) message("SICF Analytical: No a_phy input. Calculating from from Bricaud model (1998)")
    a_phy <- iop_from_oac(wavelength = wavelength,
                          par = c("chl" = c_chl))$a_phy
  }

  # Calculate fluorescence radiance by integration
  # Lf(λ_em) = (1/4π) × Σ[E0(λ_ex) × a_phy(λ_ex) × WRF(λ_ex, λ_em)]
  # If depth_resolved: multiply by 1/(Kd_ex + Kd_em) for depth integration
  if (verbose) message("SICF Analytical: Computing fluorescence radiance via WRF integration")

  Nwave <- length(wavelength)

  if (depth_resolved) {
    # C path: fused double loop — no N×N temporary matrix allocations
    # Lf[j] = (1/4π) Σ_i { E0[i] × a_phy[i] × WRF[i,j] / (Kd[i] + Kd[j]) }
    Lf <- .Call("c_sicf_depth_integrated",
                as.double(Kd_wavelength),
                as.double(a_phy),
                as.double(irrad),
                WRF_matrix,
                as.double(phi_f))

    if (verbose) message("SICF Analytical: Depth integration applied with Kd attenuation")

  } else {
    # VECTORIZED: Surface-only fluorescence
    # Pre-compute flux absorbed for all wavelengths
    flux_absorbed <- irrad * a_phy  # Vector operation

    # Compute flux emitted matrix
    # flux_emitted[i,j] = flux_absorbed[i] × WRF[i,j]
    flux_emitted_matrix <- WRF_matrix * flux_absorbed  # Broadcasting

    # Sum over excitation wavelengths to get emission at each wavelength
    Lf <- colSums(flux_emitted_matrix) / (4 * pi)

    if (verbose) message("SICF Analytical: Surface-only fluorescence")
  }

  # Convert to Rrs using the SAME Ed as semi-analytical model
  # Rrs = Lf / Ed (both models now use consistent Ed)
  Rrs_sicf <- Lf / Ed_for_rrs

  result <- list(
    wavelength = wavelength,
    E0       = E0_0m,
    Ed       = Ed_0m,
    Rrs_sicf = Rrs_sicf
  )

  # Generate depth profiles if requested
  if (plot_depth_profiles && depth_resolved) {
    if (verbose) message("SICF Analytical: Generating depth profile diagnostics")

    depth_data <- .plot_depth_profiles_sicf(
      wavelength = wavelength,
      E0_surface = irrad,
      Kd = Kd_wavelength,
      a_phy = a_phy,
      WRF_matrix = WRF_matrix,
      phi_f = phi_f,
      c_chl = c_chl,
      a_dg_443 = a_dg_443,
      bb_p_550 = bb_p_550
    )

    result <- list(
      data = result,
      depth_profiles = depth_data
    )
  }

  if (verbose) message("SICF Analytical: Completed successfully")

  return(result)
}


#' Compare Semi-Analytical and Analytical SICF Models
#'
#' Runs both SICF models across a range of parameter values and returns
#' combined results for comparison and sensitivity analysis.
#'
#' @param c_chl_range Vector of Chl concentrations to test (mg m^-3)
#' @param a_dg_443_range Vector of CDOM+NAP absorption values to test (m^-1)
#' @param bb_p_550_range Vector of particulate backscattering values to test (m^-1)
#' @param wavelength Vector of wavelengths (nm)
#' @param phi_f Fluorescence quantum yield
#' @param Ed_source Ed source ("gregg_carder" or file path)
#' @param sunzen_deg Solar zenith angle (degrees)
#' @param lat Latitude
#' @param lon Longitude
#' @param date_time POSIXct datetime
#' @param depth_resolved Use depth-integrated fluorescence? (default FALSE)
#' @param plot_depth_profiles Generate depth profile plots? (default FALSE, only if depth_resolved=TRUE)
#'
#' @return List with:
#'   - results: Data frame with all model runs
#'   - parameters: Input parameter grid
#'   - depth_profiles: (if requested) List of depth profile data
#'
#' @export
compare_sicf_models <- function(c_chl_range,
                                 a_dg_443_range,
                                 bb_p_550_range = NULL,
                                 wavelength = seq(400, 800, 10),
                                 phi_f = 0.02,
                                 Ed_source = "gregg_carder",
                                 sunzen_deg = 30,
                                 lat = 49,
                                 lon = -68,
                                 date_time = as.POSIXct("2019-08-18 20:50:00", tz = "UTC"),
                                 depth_resolved = FALSE,
                                 plot_depth_profiles = FALSE) {

  message("\n========================================")
  message("SICF MODEL COMPARISON")
  if (depth_resolved) {
    message("MODE: DEPTH-RESOLVED with attenuation")
  } else {
    message("MODE: SURFACE-ONLY (no depth integration)")
  }
  message("========================================\n")

  # Validate inputs for depth_resolved mode
  if (depth_resolved && is.null(bb_p_550_range)) {
    stop("bb_p_550_range required when depth_resolved=TRUE")
  }

  # Set default bb_p_550 for surface-only mode
  if (!depth_resolved) {
    bb_p_550_range <- 0.001  # Dummy value, not used
  }

  # Create parameter grid
  param_grid <- expand.grid(
    c_chl = c_chl_range,
    a_dg_443 = a_dg_443_range,
    bb_p_550 = bb_p_550_range,
    stringsAsFactors = FALSE
  )

  n_scenarios <- nrow(param_grid)
  message(sprintf("Running %d scenarios...\n", n_scenarios))

  # Storage for results
  all_results <- list()
  depth_profile_list <- list()

  # Loop over parameter combinations
  for (i in 1:n_scenarios) {
    chl_i <- param_grid$c_chl[i]
    adg_i <- param_grid$a_dg_443[i]
    bbp_i <- param_grid$bb_p_550[i]

    if (depth_resolved) {
      message(sprintf("Scenario %d/%d: Chl=%.2f mg/m³, a_dg(443)=%.3f m⁻¹, bb_p(550)=%.4f m⁻¹",
                      i, n_scenarios, chl_i, adg_i, bbp_i))
    } else {
      message(sprintf("Scenario %d/%d: Chl=%.2f mg/m³, a_dg(443)=%.3f m⁻¹",
                      i, n_scenarios, chl_i, adg_i))
    }

    # Semi-Analytical Model
    message("  Running Semi-Analytical model...")
    sa_result <- sicf_semi_analytical(
      c_chl = chl_i,
      a_dg_443 = adg_i,
      wavelength = wavelength,
      phi_f = phi_f,
      Ed_source = Ed_source,
      use_analytic_Ed = TRUE,
      sunzen_deg = sunzen_deg,
      lat = lat,
      lon = lon,
      date_time = date_time
    )
    sa_result$model <- "Semi-Analytical"
    sa_result$c_chl <- chl_i
    sa_result$a_dg_443 <- adg_i
    sa_result$bb_p_550 <- bbp_i

    # Analytical Model
    message("  Running Analytical model...")

    # Generate depth profiles only for first scenario (to avoid too many plots)
    plot_profiles_now <- plot_depth_profiles && (i == 1)

    analytical_result <- sicf_analytical(
      c_chl = chl_i,
      a_phy = NULL,  # Will be calculated from chl
      wavelength = wavelength,
      phi_f = phi_f,
      Ed_source = Ed_source,
      use_analytic_Ed = TRUE,
      sunzen_deg = sunzen_deg,
      lat = lat,
      lon = lon,
      date_time = date_time,
      scalar_irradiance = TRUE,
      depth_resolved = depth_resolved,
      a_dg_443 = adg_i,
      bb_p_550 = bbp_i,
      plot_depth_profiles = plot_profiles_now
    )

    # Extract data if depth profiles were generated
    if (is.list(analytical_result) && !is.null(analytical_result$depth_profiles)) {
      depth_profile_list[[i]] <- analytical_result$depth_profiles
      analytical_result <- analytical_result$data
    }

    analytical_result$model <- "Analytical"
    analytical_result$c_chl <- chl_i
    analytical_result$a_dg_443 <- adg_i
    analytical_result$bb_p_550 <- bbp_i

    # Combine
    all_results[[2*i - 1]] <- sa_result
    all_results[[2*i]] <- analytical_result

    message("")
  }

  # Combine all results
  combined_results <- do.call(rbind, all_results)

  message("========================================")
  message("COMPARISON COMPLETE")
  message("========================================\n")

  output <- list(
    results = combined_results,
    parameters = param_grid
  )

  # Add depth profiles if generated
  if (length(depth_profile_list) > 0) {
    output$depth_profiles <- depth_profile_list
  }

  return(output)
}


#' Plot SICF Model Comparison
#'
#' Creates visualization comparing semi-analytical and analytical SICF models.
#'
#' @param comparison_results Output from compare_sicf_models()
#' @param plot_type Type of plot: "spectra", "peak", or "integrated"
#'
#' @return ggplot object
#'
#' @export
plot_sicf_comparison <- function(comparison_results, plot_type = "spectra") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 required for plotting")
  }

  library(ggplot2)

  df <- comparison_results$results

  if (plot_type == "spectra") {
    # Full spectral comparison
    p <- ggplot(df, aes(x = wavelength, y = Rrs_sicf,
                        color = model, linetype = model)) +
      geom_line(linewidth = 0.8) +
      facet_grid(a_dg_443 ~ c_chl,
                 labeller = label_both) +
      scale_color_manual(values = c("Semi-Analytical" = "#E41A1C",
                                     "Analytical" = "#377EB8")) +
      labs(
        title = "SICF Model Comparison: Spectral Profiles",
        x = "Wavelength (nm)",
        y = expression(paste(R[rs]^{SICF}, " (sr"^-1, ")")),
        color = "Model",
        linetype = "Model"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      )

  } else if (plot_type == "peak") {
    # Peak fluorescence comparison
    df_peak <- df %>%
      dplyr::group_by(model, c_chl, a_dg_443) %>%
      dplyr::summarize(
        peak_Rrs = max(Rrs_sicf, na.rm = TRUE),
        .groups = "drop"
      )

    p <- ggplot(df_peak, aes(x = c_chl, y = peak_Rrs,
                              color = model, shape = factor(a_dg_443))) +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(model, a_dg_443))) +
      scale_x_log10() +
      scale_color_manual(values = c("Semi-Analytical" = "#E41A1C",
                                     "Analytical" = "#377EB8")) +
      labs(
        title = "SICF Peak Fluorescence Comparison",
        x = "Chlorophyll-a (mg/m³, log scale)",
        y = expression(paste("Peak ", R[rs]^{SICF}, " (sr"^-1, ")")),
        color = "Model",
        shape = "a_dg(443)"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "right")

  } else if (plot_type == "integrated") {
    # Integrated fluorescence signal
    df_int <- df %>%
      dplyr::group_by(model, c_chl, a_dg_443) %>%
      dplyr::summarize(
        integrated_Rrs = sum(Rrs_sicf, na.rm = TRUE),
        .groups = "drop"
      )

    p <- ggplot(df_int, aes(x = c_chl, y = integrated_Rrs,
                             color = model, shape = factor(a_dg_443))) +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(model, a_dg_443))) +
      scale_x_log10() +
      scale_color_manual(values = c("Semi-Analytical" = "#E41A1C",
                                     "Analytical" = "#377EB8")) +
      labs(
        title = "SICF Integrated Signal Comparison",
        x = "Chlorophyll-a (mg/m³, log scale)",
        y = expression(paste("Integrated ", R[rs]^{SICF})),
        color = "Model",
        shape = "a_dg(443)"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "right")

  } else {
    stop("plot_type must be 'spectra', 'peak', or 'integrated'")
  }

  return(p)
}
