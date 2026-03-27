#' Combined Forward Model: Elastic + Inelastic (SICF) Scattering
#'
#' Computes total Rrs including both elastic scattering (AM03) and
#' Sun-Induced Chlorophyll Fluorescence (SICF). Also returns intermediate
#' optical properties: Ed, E0, and PAR.
#'
#' @param wavelength Vector of wavelengths [nm]
#' @param iop List with elements `a` and `bb`, same length as wavelength
#' @param water_type Either 1 or 2 (default = 2)
#' @param theta_sun Sun zenith angle [degrees]
#' @param theta_view Sensor view angle [degrees]
#' @param h_w Optional: water depth [m] (enables shallow mode)
#' @param r_b Optional: bottom reflectance vector (same length as wavelength)
#' @param chl Chlorophyll-a concentration [mg/m³] (required for SICF)
#' @param a_dg_443 CDOM+NAP absorption at 443nm [1/m] (required for SICF)
#' @param phi_f Fluorescence quantum yield (default = 0.02)
#' @param include_sicf Logical, include SICF component? (default = TRUE)
#' @param sicf_model Character, SICF model to use: "semi_analytical" or "analytical" (default = "semi_analytical")
#' @param depth_integration Logical, use depth-integrated analytical model? Only applies when sicf_model = "analytical" (default = FALSE)
#' @param lat Latitude [decimal degrees] (default = 49)
#' @param lon Longitude [decimal degrees] (default = -68)
#' @param date_time POSIXct datetime in UTC (default = current time)
#' @param return_components Logical, return all components separately? (default = FALSE)
#' @param verbose Logical, print diagnostic messages? (default = FALSE)
#'
#' @return If return_components = FALSE: numeric vector of total Rrs
#'         If return_components = TRUE: list with:
#'           - rrs_total: Total Rrs (elastic + SICF)
#'           - rrs_elastic: Elastic component only
#'           - rrs_sicf: SICF component only
#'           - Ed_0m: Subsurface downwelling irradiance [W/m²/nm]
#'           - E0_0m: Scalar irradiance at surface [W/m²/nm]
#'           - PAR: Photosynthetically Available Radiation [μmol photons/m²/s]
#'
#' @export
forward_am03_sicf <- function(wavelength, iop, water_type = 2,
                               theta_sun, theta_view,
                               h_w = NULL, r_b = NULL,
                               chl, a_dg_443,
                               phi_f = 0.02,
                               include_sicf = TRUE,
                               sicf_model = "analytical",
                               depth_integration = FALSE,
                               lat = 49, lon = -68,
                               date_time = Sys.time(),
                               return_components = FALSE,
                               verbose = FALSE) {

  # ========================================================================
  # PART 0: Robust Input Validation
  # ========================================================================

  # Validate wavelength
  if (!is.numeric(wavelength) || length(wavelength) == 0) {
    stop("'wavelength' must be a non-empty numeric vector")
  }
  if (any(is.na(wavelength)) || any(wavelength <= 0)) {
    stop("'wavelength' must contain positive values without NAs")
  }

  # Validate IOP
  if (!is.list(iop) || !all(c("a", "bb") %in% names(iop))) {
    stop("'iop' must be a list with elements 'a' and 'bb'")
  }
  if (!is.numeric(iop$a) || !is.numeric(iop$bb)) {
    stop("'iop$a' and 'iop$bb' must be numeric vectors")
  }
  if (length(iop$a) != length(wavelength) || length(iop$bb) != length(wavelength)) {
    stop("'iop$a' and 'iop$bb' must have the same length as 'wavelength'")
  }
  if (any(iop$a < 0, na.rm = TRUE) || any(iop$bb < 0, na.rm = TRUE)) {
    stop("'iop$a' and 'iop$bb' must be non-negative")
  }

  # Validate water_type
  if (!water_type %in% c(1, 2)) {
    stop("'water_type' must be either 1 or 2")
  }

  # Validate angles
  if (!is.numeric(theta_sun) || length(theta_sun) != 1 ||
      theta_sun < 0 || theta_sun > 90) {
    stop("'theta_sun' must be a single numeric value between 0 and 90 degrees")
  }
  if (!is.numeric(theta_view) || length(theta_view) != 1 ||
      theta_view < 0 || theta_view > 90) {
    stop("'theta_view' must be a single numeric value between 0 and 90 degrees")
  }

  # Validate shallow water parameters
  if (!is.null(h_w)) {
    if (!is.numeric(h_w) || length(h_w) != 1 || h_w <= 0) {
      stop("'h_w' must be a positive numeric value or NULL")
    }
    if (is.null(r_b)) {
      warning("'h_w' is specified but 'r_b' is NULL - shallow water mode may not work correctly")
    }
  }
  if (!is.null(r_b)) {
    if (!is.numeric(r_b) || length(r_b) != length(wavelength)) {
      stop("'r_b' must be a numeric vector with the same length as 'wavelength' or NULL")
    }
  }

  # Validate SICF parameters
  if (!is.logical(include_sicf) || length(include_sicf) != 1) {
    stop("'include_sicf' must be a single logical value (TRUE or FALSE)")
  }

  if (include_sicf) {
    if (missing(chl) || !is.numeric(chl) || length(chl) != 1 || chl < 0) {
      stop("'chl' must be a non-negative numeric value when include_sicf = TRUE")
    }
    if (missing(a_dg_443) || !is.numeric(a_dg_443) || length(a_dg_443) != 1 || a_dg_443 < 0) {
      stop("'a_dg_443' must be a non-negative numeric value when include_sicf = TRUE")
    }
    if (!is.numeric(phi_f) || length(phi_f) != 1 || phi_f < 0 || phi_f > 0.1) {
      stop("'phi_f' must be a numeric value between 0 and 0.1")
    }

    # Validate sicf_model
    if (!sicf_model %in% c("semi_analytical", "analytical")) {
      stop("'sicf_model' must be either 'semi_analytical' or 'analytical'")
    }

    # Validate depth_integration
    if (!is.logical(depth_integration) || length(depth_integration) != 1) {
      stop("'depth_integration' must be a single logical value (TRUE or FALSE)")
    }
    if (depth_integration && sicf_model != "analytical") {
      warning("'depth_integration' is only used when sicf_model = 'analytical'. Ignoring.")
    }
  }

  # Validate geographic parameters
  if (!is.numeric(lat) || length(lat) != 1 || lat < -90 || lat > 90) {
    stop("'lat' must be a numeric value between -90 and 90 degrees")
  }
  if (!is.numeric(lon) || length(lon) != 1 || lon < -180 || lon > 180) {
    stop("'lon' must be a numeric value between -180 and 180 degrees")
  }

  # Validate date_time
  if (!inherits(date_time, "POSIXct") && !inherits(date_time, "POSIXt")) {
    tryCatch({
      date_time <- as.POSIXct(date_time, tz = "UTC")
    }, error = function(e) {
      stop("'date_time' must be a POSIXct object or convertible to one")
    })
  }

  # Validate return_components
  if (!is.logical(return_components) || length(return_components) != 1) {
    stop("'return_components' must be a single logical value (TRUE or FALSE)")
  }

  if (include_sicf) {
    if (missing(chl) || missing(a_dg_443)) {
      stop("chl and a_dg_443 are required when include_sicf = TRUE")
    }
  }

  # ========================================================================
  # PART 1: Compute Elastic Scattering Component (AM03)
  # ========================================================================

  rrs_elastic <- forward_am03(
    wavelength = wavelength,
    iop = iop,
    water_type = water_type,
    theta_sun = theta_sun,
    theta_view = theta_view,
    h_w = h_w,
    r_b = r_b
  )

  # # ========================================================================
  # # PART 2: Calculate Ed and E0 using Gregg & Carder
  # # ========================================================================
  #
  # # Load Gregg & Carder data
  # Cops::GreggCarder.data()
  #
  # # Extract date/time components
  # jday_no <- lubridate::yday(date_time)
  # time_no <- format(date_time, "%T")
  # time_dec <- sapply(strsplit(as.character(time_no), ":"), function(x) {
  #   x <- as.numeric(x)
  #   x[1] + x[2]/60
  # })
  #
  # # Calculate Ed at 0+ using Gregg & Carder
  # Ed_gc <- tryCatch({
  #   Cops::GreggCarder.f(
  #     the = theta_sun,
  #     lam.sel = wavelength,
  #     hr = time_dec,
  #     jday = jday_no,
  #     rlon = lon,
  #     rlat = lat,
  #     debug = FALSE
  #   )
  # }, error = function(e) {
  #   warning("Gregg & Carder Ed calculation failed: ", e$message)
  #   return(NULL)
  # })
  #
  # if (!is.null(Ed_gc) && !all(is.na(Ed_gc))) {
  #   Ed0_0p <- Ed_gc$Ed      # Total Ed at 0+
  #   Ed0_dir_0p <- Ed_gc$Edir  # Direct component
  #   Ed0_dif_0p <- Ed_gc$Edif  # Diffuse component
  #
  #   # Calculate Fresnel reflectance
  #   rhoF <- Cops::GreggCarder.sfcrfl(rad = 180/pi, theta = theta_sun, ws = 5)
  #
  #   # Convert Ed from 0+ to 0- (subsurface)
  #   Ed_0m <- (Ed0_dir_0p * (1 - rhoF$rod)) + (Ed0_dif_0p * (1 - rhoF$ros))
  #
  #   # Calculate E0 (scalar irradiance) at surface
  #   # E0 ≈ Ed × (1 + 1/μ_d) where μ_d is average cosine for diffuse light
  #   mu_d <- 0.85  # Typical value for clear sky
  #   E0_0m <- Ed_0m * (1 + 1/mu_d)
  #
  # } else {
  #   warning("Ed calculation failed, returning NA for Ed, E0, and PAR")
  #   Ed_0m <- rep(NA_real_, length(wavelength))
  #   E0_0m <- rep(NA_real_, length(wavelength))
  # }


  # ========================================================================
  # PART 3: Calculate SICF Component
  # ========================================================================

  if (include_sicf) {

    # Switch between SICF models
    if (sicf_model == "semi_analytical") {

      # Use semi-analytical SICF model (Gilerson et al. 2007)
      sicf_result <- sicf_semi_analytical(
        c_chl = chl,
        a_dg_443 = a_dg_443,
        wavelength = wavelength,
        phi_f = phi_f,
        Ed_source = "gregg_carder",
        use_analytic_Ed = TRUE,
        sunzen_deg = theta_sun,
        lat = lat,
        lon = lon,
        date_time = date_time,
        return_radiance = FALSE,
        verbose = verbose
      )

      rrs_sicf <- sicf_result$Rrs_sicf

    } else if (sicf_model == "analytical") {

      # Use analytical SICF model (WRF-based)
      if (depth_integration) {
        # Depth-resolved analytical model (requires bb_p_550)
        bb_p_550 <- if (!is.null(iop$bb) && length(iop$bb) > 0) {
          # Extract bb_p_550 from IOP (approximate at 550nm)
          idx_550 <- which.min(abs(wavelength - 550))
          iop$bb[idx_550]
        } else {
          0.005  # Default value if not available
        }

        sicf_result <- sicf_analytical(
          c_chl = chl,
          a_dg_443 = a_dg_443,
          bb_p_550 = bb_p_550,
          wavelength = wavelength,
          phi_f = phi_f,
          Ed_source = "gregg_carder",
          use_analytic_Ed = TRUE,
          sunzen_deg = theta_sun,
          lat = lat,
          lon = lon,
          date_time = date_time,
          depth_resolved = TRUE,
          verbose = verbose
        )
      } else {
        # Surface analytical model
        sicf_result <- sicf_analytical(
          c_chl = chl,
          a_dg_443 = a_dg_443,
          wavelength = wavelength,
          phi_f = phi_f,
          Ed_source = "gregg_carder",
          use_analytic_Ed = TRUE,
          sunzen_deg = theta_sun,
          lat = lat,
          lon = lon,
          date_time = date_time,
          depth_resolved = FALSE,
          verbose = verbose
        )
      }

      rrs_sicf <- sicf_result$Rrs_sicf
    }

    # Total Rrs = Elastic + SICF
    rrs_total <- rrs_elastic + rrs_sicf
    Ed_0m <- sicf_result$Ed
    E0_0m <- sicf_result$E0

  } else {
    rrs_sicf <- rep(0, length(wavelength))
    rrs_total <- rrs_elastic
    Ed_0m <- NULL
    E0_0m <- NULL
  }


  # ========================================================================
  # PART 4: Calculate PAR (Photosynthetically Available Radiation)
  # ========================================================================

  # PAR is integrated over 400-700 nm
  # Convert from W/m²/nm to μmol photons/m²/s
  # E [μmol/m²/s] = ∫ Ed(λ) × λ/(h×c×Na) dλ
  # Where h = Planck's constant, c = speed of light, Na = Avogadro's number

  if (!all(is.na(Ed_0m))) {
    # Constants
    h <- 6.62607015e-34  # Planck's constant [J·s]
    c <- 2.99792458e8    # Speed of light [m/s]
    Na <- 6.02214076e23  # Avogadro's number [1/mol]

    # Find PAR wavelengths (400-700 nm)
    par_idx <- which(wavelength >= 400 & wavelength <= 700)

    if (length(par_idx) > 0) {
      wl_par <- wavelength[par_idx]
      Ed_par <- Ed_0m[par_idx]

      # Convert wavelength to meters
      wl_m <- wl_par * 1e-9

      # Energy per photon: E_photon = h × c / λ
      E_photon <- h * c / wl_m

      # Convert Ed [W/m²/nm] to photon flux [photons/m²/s/nm]
      # Ed [W/m²/nm] / E_photon [J/photon] = [photons/m²/s/nm]
      photon_flux <- Ed_par / E_photon

      # Integrate over wavelength (trapezoidal rule)
      # Result in photons/m²/s
      PAR_photons <- pracma::trapz(wl_par, photon_flux)

      # Convert to μmol photons/m²/s
      PAR <- PAR_photons / Na * 1e6

    } else {
      PAR <- NA_real_
    }
  } else {
    PAR <- NA_real_
  }

  # ========================================================================
  # PART 5: Return Results
  # ========================================================================

  if (return_components) {
    return(list(
      rrs_total = rrs_total,
      rrs_elastic = rrs_elastic,
      rrs_sicf = rrs_sicf,
      Ed_0m = Ed_0m,
      E0_0m = E0_0m,
      PAR = PAR,
      wavelength = wavelength
    ))
  } else {
    return(rrs_total)
  }
}


#' Input Preparer for Combined Forward Model (AM03 + SICF)
#'
#' Prepares inputs for the combined forward model that includes both
#' elastic scattering and SICF components.
#'
#' @param par Named vector of parameters
#' @param rrs Data frame with wavelength column
#'
#' @return List of inputs for forward_am03_sicf
#'
#' @keywords internal
input_am03_sicf <- function(par, rrs, par_meta = NULL) {
  # Extract metadata from par_meta (separate vector)
  sicf_model <- if (!is.null(par_meta) && "sicf_model" %in% names(par_meta)) {
    as.character(par_meta["sicf_model"])
  } else {
    "semi_analytical"
  }

  depth_integration <- if (!is.null(par_meta) && "depth_integration" %in% names(par_meta)) {
    as.logical(par_meta["depth_integration"])
  } else {
    FALSE
  }

  # Auto-initialize WRF cache for analytical model (happens once per wavelength set)
  if (sicf_model == "analytical") {
    .ensure_WRF_cache(rrs$wavelength, verbose = TRUE)
  }

  # par is already a pure numeric vector - use it directly
  # Calculate IOPs from OACs (requires numeric named vector)
  iop <- iop_from_oac(rrs$wavelength, par)

  # Extract benthic reflectance fractions
  r_b_fraction_vec <- par[grep("^r_rs_b", names(par))]

  # Check if all r_b fractions are NA or NULL
  if (all(is.na(r_b_fraction_vec)) || length(r_b_fraction_vec) == 0) {
    r_b <- NULL
  } else {
    r_b <- compute_r_rs_b_lmm(fractions = r_b_fraction_vec)
  }

  # Handle h_w
  h_w <- if (is.null(par["h_w"]) || is.na(par["h_w"])) NULL else par["h_w"]

  # Extract phi_f (quantum yield)
  phi_f <- if ("phi_f" %in% names(par)) par["phi_f"] else 0.02

  # Extract chl and a_dg for SICF
  chl <- par["chl"]

  # Calculate a_dg at 443nm
  # a_dg = a_g + a_nap at 443nm (approximately 440nm)
  a_dg_440 <- if ("a_dg_440" %in% names(par)) par["a_dg_440"] else 0
  a_nap_440 <- if ("a_nap_440" %in% names(par)) par["a_nap_440"] else 0
  a_dg_443 <- a_dg_440 + a_nap_440

  # Extract geometry info
  lat <- if ("lat" %in% names(par)) par["lat"] else 49
  lon <- if ("lon" %in% names(par)) par["lon"] else -68

  # Extract or create date_time
  if ("date_time" %in% names(par)) {
    date_time <- as.POSIXct(par["date_time"], origin = "1970-01-01", tz = "UTC")
  } else {
    date_time <- Sys.time()
  }

  # sicf_model and depth_integration already extracted at the beginning

  # Ensure water_type exists and is valid
  if (!("water_type" %in% names(par)) || is.na(par["water_type"])) {
    stop("'water_type' parameter is missing or NA. Available parameters: ",
         paste(names(par), collapse = ", "))
  }

  list(
    wavelength = rrs$wavelength,
    iop = iop,
    water_type = as.numeric(par["water_type"]),
    theta_view = as.numeric(par["theta_view"]),
    theta_sun = as.numeric(par["theta_sun"]),
    h_w = h_w,
    r_b = r_b,
    chl = chl,
    a_dg_443 = a_dg_443,
    phi_f = phi_f,
    include_sicf = TRUE,
    sicf_model = sicf_model,
    depth_integration = depth_integration,
    lat = lat,
    lon = lon,
    date_time = date_time,
    return_components = FALSE
  )
}


