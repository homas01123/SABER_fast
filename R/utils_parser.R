#' Interpolate non-water IOPs to a target wavelength grid
#'
#' Linearly interpolates absorption and fits a power-law to backscattering,
#' both resampled to \code{wavelength}.
#'
#' @param a           Data frame with columns \code{wavelength} and \code{a}
#'   [m\eqn{^{-1}}].
#' @param bb          Data frame with columns \code{wavelength} and \code{bb}
#'   [m\eqn{^{-1}}].
#' @param wavelength  Numeric vector of target wavelengths [nm].
#' @param verbose     Logical.  Plot interpolated vs. observed IOPs?
#'   Default \code{FALSE}.
#'
#' @return Tibble with columns \code{wavelength}, \code{a}, \code{bb}.
#'
#' @author Soham Mukherjee, Raphael Mabit
#' @export

parse_iop <- function(
    a,
    bb,
    wavelength,
    verbose = F) {
  a_non_water <- approx(
    x = a$wavelength, y = a$a,
    xout = wavelength, method = "linear"
  )$y

  # For bb we create a power law model and use it for interpolation

  model <- nls(
    bb ~ a * wavelength^b,
    start = list(a = bb$bb[3], b = 1),
    data = bb,
    control = list(maxiter = 100, warnOnly = T)
  )

  # bbp_hs[j,i] = refbbp[j,1]*((waveletngth_hs[i,1]/555)^-(refexponent[j,1]))
  bb_non_water <- coef(model)[1] * wavelength^coef(model)[2]

  # bb_non_water <- approx(x= iop$wavelength, y = iop$bb,
  #                        xout = wavelength, method = "linear")$y

  iop <- tibble(
    wavelength,
    "a" = a_non_water,
    "bb" = bb_non_water
  )

  if (verbose) {
    # Plot interpolated non water absorption
    plot(a$wavelength, a$a,
      xlab = "wavelength",
      ylab = "non-water absorption [m^-1]"
    )
    lines(wavelength, a_non_water, col = "red", lwd = 3)

    # Plot power law fitted backscatering
    plot(bb$wavelength, bb$bb,
      xlab = "wavelength",
      ylab = "non-water backscatter [m^-1]"
    )
    lines(wavelength, bb_non_water, col = "red", lwd = 3)
  }

  return(iop)
}

#' parse_inverse_name
#'
#' a function to parse and prepare the init name to be passed to optimization function

parse_inverse_parameter <- function(
    par_df,
    optim_mtd,
    lower_b = NULL,
    upper_b = NULL,
    verbose = F) {
  if (optim_mtd == "L-BFGS-B" & is.null(lower_b)) {
    lower_b <- dplyr::case_when(
      par_df$name %in% c("chl", "a_g_440", "bb_p_550") ~ par_df$value - 0.8 * par_df$value,
      par_df$name == "h_w" ~ 0.5,
      stringr::str_detect(par_df$name, "^rb_") ~ 0.01,
      par_df$name == "sd" ~ 1e-5,
      TRUE ~ NA_real_
    )
  }

  if (optim_mtd == "L-BFGS-B" & is.null(upper_b)) {
    upper_b <- dplyr::case_when(
      par_df$name %in% c("chl", "a_g_440", "bb_p_550") ~ par_df$value + 5 * par_df$value,
      par_df$name == "h_w" ~ 10,
      stringr::str_detect(par_df$name, "^rb_") ~ 1,
      par_df$name == "sd" ~ 1,
      TRUE ~ NA_real_
    )
  }

  if (optim_mtd == "L-BFGS-B" & verbose) {
    rlang::inform(
      glue::glue("{par_df$name} lower boundary set at: {lower_b} upper boundary set at: {upper_b}")
    )
  }

  list(
    par = par_df$value,
    names = par_df$name,
    lower = lower_b,
    upper = upper_b,
    parscale = abs(par_df$value)
  )
}

# Inversion configuration constants (package defaults)

# Default set of inverted OAC parameters
.PAR_INVERSED_DEEP <- c("chl", "a_dg_440", "bb_p_550")

# Universe of default bounds covering all parameters that can be promoted from
# par_fixed to par_inv. Extend this table to add new promotable parameters.
.ALL_LOWER <- c(
  chl        = 0.03,   a_dg_440   = 0.002,  bb_p_550   = 0.0004,
  a_dg_s     = 0.003,  bb_p_gamma = 0.1,    phi_f      = 0.005,
  h_w        = 0.1
)
.ALL_UPPER <- c(
  chl        = 1200,   a_dg_440   = 6.0,    bb_p_550   = 0.1,
  a_dg_s     = 0.030,  bb_p_gamma = 2.0,    phi_f      = 0.06,
  h_w        = 30
)

# Parameters that are structural / geometry inputs and can never be inverted
.FIXED_ONLY <- c("lat", "lon", "water_type", "theta_view")

# SD bounds are on the raw Rrs scale [sr^-1].
# Typical Rrs range 1e-4 – 5e-3 sr^-1; noise is 0.5–50% of signal.
.SD_LOWER <- 1e-5
.SD_UPPER <- 5e-3

# Internal spectral helpers
# KD490 empirical algorithm (Mueller 2000 / Lee et al.)
.estimate_kd490 <- function(rrs_df) {
  r <- approx(rrs_df$wavelength, rrs_df$rrs_0m, xout = c(443, 490, 560), rule = 2)$y
  if (any(!is.finite(r)) || any(r <= 0)) return(NA_real_)
  X  <- log10(max(r[1], r[2]) / r[3])
  kd <- 10^(-0.8515 - 1.8263*X + 0.4683*X^2 - 0.8185*X^3 - 0.4278*X^4)
  max(0.005, kd)
}

# Multi-criterion optical shallow classification.
# A pixel is treated as optically shallow only when ALL three conditions hold:
#   (1) physical depth  <= depth_threshold
#   (2) 2 * Kd(490) * depth < optical_depth_threshold  (bottom visible)
#   (3) Rrs(665) / Rrs(560) >= red_green_min  (benthic red signal present)
.is_optically_shallow <- function(rrs_df, depth_m,
                                   depth_threshold         = 20,
                                   optical_depth_threshold = 4.0,
                                   red_green_min           = 0.10) {
  if (is.na(depth_m) || depth_m > depth_threshold) return(FALSE)
  r665 <- approx(rrs_df$wavelength, rrs_df$rrs_0m, xout = 665, rule = 2)$y
  r560 <- approx(rrs_df$wavelength, rrs_df$rrs_0m, xout = 560, rule = 2)$y
  if (!is.finite(r665) || !is.finite(r560) || r560 <= 0) return(FALSE)
  if (r665 / r560 < red_green_min) return(FALSE)
  kd <- .estimate_kd490(rrs_df)
  if (is.na(kd)) return(FALSE)
  2 * kd * abs(depth_m) < optical_depth_threshold
}

# Per-observation empirical prior modes for OAC parameters.
# chl     : OC3 polynomial (O'Reilly et al. 1998)
# a_dg_440: blue-ratio absorption proxy
# bb_p_550: green Rrs proxy (Gordon et al. 1988)
.estimate_prior_modes <- function(rrs_df) {
  defaults <- c(chl = 1.0, a_dg_440 = 0.10, bb_p_550 = 0.005)
  r <- tryCatch(
    setNames(
      approx(rrs_df$wavelength, rrs_df$rrs_0m, xout = c(412, 443, 490, 560), rule = 2)$y,
      c("r412", "r443", "r490", "r560")
    ),
    error = function(e) NULL
  )
  if (is.null(r) || any(!is.finite(r)) || any(r <= 0)) return(defaults)

  X        <- log10(max(r["r443"], r["r490"]) / r["r560"])
  chl      <- max(0.02, min(80, 10^(0.3272 - 2.994*X + 2.722*X^2 - 1.226*X^3 - 0.597*X^4)))
  a_dg_440 <- max(0.002, min(2.5, 10^(-1.3 - 1.8 * log10(r["r412"] / r["r443"]))))
  bb_p_550 <- max(0.001, min(0.20, r["r560"] * 4.0))

  c(chl = unname(chl), a_dg_440 = unname(a_dg_440), bb_p_550 = unname(bb_p_550))
}

#' Build a per-pixel inversion configuration
#'
#' Constructs bounds, initial values, and a BayesianTools prior for a single
#' pixel.  Automatically selects optically-shallow treatment and appends
#' benthic fraction parameters when depth and spectral criteria are met.
#' Benthic reflectance classes must be loaded beforehand with
#' \code{\link{select_benthic_classes}}.
#'
#' @param depth_m   Bathymetric depth \[m, positive down\], or \code{NULL} /
#'   \code{NA} to force optically-deep treatment.
#' @param rrs_df    Data frame with columns \code{wavelength} and
#'   \code{rrs_0m} \[sr\eqn{^{-1}}\].
#' @param par_inv   Character vector of OAC parameters to invert.  Default
#'   \code{c("chl", "a_dg_440", "bb_p_550")}.  Additional invertible OAC
#'   parameters: \code{a_dg_s}, \code{bb_p_gamma}, \code{phi_f}.
#'   \strong{Do not} include shallow-branch parameters (\code{h_w},
#'   \code{mix_sand}, \code{r_rs_b_*}) here — they are appended automatically.
#'   Use \code{fixed} to pin any of them to a constant value.
#' @param lower     Optional named numeric lower-bound overrides (applies to
#'   both OAC and auto-appended shallow params).
#' @param upper     Optional named numeric upper-bound overrides (applies to
#'   both OAC and auto-appended shallow params).
#' @param fixed     Optional named list of parameters to hold constant (remove
#'   from the inverted set and pass to \code{par_fixed}).  Can target any
#'   parameter, e.g. \code{list(h_w = 5, a_dg_s = 0.015)}.  Parameters listed
#'   here must not also appear in \code{par_inv}.
#' @param sicf      Logical.  Use the SICF forward model (\code{"am03_sicf"})
#'   and automatically promote \code{phi_f} into \code{par_inv}?
#'   Default \code{FALSE}.
#' @param sicf_model  SICF variant: \code{"analytical"} (default) or
#'   \code{"semi_analytical"}.
#' @param depth_integration  Logical.  Depth-integrated SICF?  Default
#'   \code{TRUE}.
#' @param lat,lon,date_time  Latitude [°N], longitude [°E], and UTC
#'   \code{POSIXct} timestamp.  Required when \code{sicf = TRUE}.
#' @param theta_sun  Solar zenith angle [°].  Default \code{30}.
#' @param depth_threshold          Depth [m] threshold for shallow detection.
#'   Default \code{20}.
#' @param optical_depth_threshold  \eqn{2 K_d z} threshold for shallow
#'   detection.  Default \code{4.0}.
#' @param red_green_min  Minimum \eqn{R_{rs}(665)/R_{rs}(560)} ratio for
#'   benthic-signal detection.  Default \code{0.10}.
#' @param water_type  Integer coastal-water type for the forward model.
#'   Default \code{2}.
#' @param force_shallow  Logical.  If \code{TRUE}, skip the spectral/depth
#'   shallow-detection test and always use the shallow forward model.
#'   Default \code{FALSE}.
#'
#' @return Named list: \code{shallow}, \code{h_w}, \code{fwd_model},
#'   \code{par_inv}, \code{lower_b}, \code{upper_b}, \code{init_val},
#'   \code{bt_prior}, \code{log_prior_fn}, \code{par_fixed}.
#'
#' @examples
#' \dontrun{
#' cfg <- make_inversion_params(depth_m = NA, rrs_df = rrs_df)
#' est <- inverse_mcmc(rrs_df, cfg$fwd_model, cfg$par_inv,
#'   prior = cfg$bt_prior, par_fixed = cfg$par_fixed)
#' }
#'
#' @export
make_inversion_params <- function(
    depth_m,
    rrs_df,
    par_inv                 = .PAR_INVERSED_DEEP,
    lower                   = NULL,
    upper                   = NULL,
    fixed                   = NULL,
    sicf                    = FALSE,
    sicf_model              = "analytical",
    depth_integration       = TRUE,
    water_type              = 2,
    lat                     = NULL,
    lon                     = NULL,
    date_time               = NULL,
    theta_sun               = 30,
    depth_threshold         = 20,
    optical_depth_threshold = 4.0,
    red_green_min           = 0.10,
    force_shallow           = FALSE) {

  # ---- Validate: structural params can never be inverted ------------------
  bad <- intersect(par_inv, .FIXED_ONLY)
  if (length(bad) > 0)
    stop("Parameters cannot be inverted: ", paste(bad, collapse = ", "))

  # ---- Validate: shallow-specific params must not appear in par_inv -------
  # h_w, mix_sand, r_rs_b_* are auto-appended in the shallow branch.
  # To hold any of them constant, use: fixed = list(h_w = 5)
  auto_shallow <- par_inv[par_inv %in% c("h_w", "mix_sand") |
                            grepl("^r_rs_b_", par_inv)]
  if (length(auto_shallow) > 0)
    stop("Shallow-branch parameters are auto-appended and cannot be in par_inv: ",
         paste(auto_shallow, collapse = ", "),
         ".\nTo pin a value use: fixed = list(", auto_shallow[1], " = <value>).")

  # ---- Validate: fixed must be a named list, no overlap with par_inv ------
  if (!is.null(fixed)) {
    if (!is.list(fixed) || is.null(names(fixed)))
      stop("`fixed` must be a named list, e.g. list(h_w = 5, a_dg_440 = 0.01).")
    overlap <- intersect(par_inv, names(fixed))
    if (length(overlap) > 0)
      stop("Parameters appear in both par_inv and fixed: ",
           paste(overlap, collapse = ", "))
  }

  # ---- Validate: sicf requires geo metadata --------------------------------
  if (sicf) {
    missing_geo <- c(
      if (is.null(lat))       "lat",
      if (is.null(lon))       "lon",
      if (is.null(date_time)) "date_time"
    )
    if (length(missing_geo) > 0)
      stop("sicf = TRUE requires: ", paste(missing_geo, collapse = ", "))
  }

  # ---- Step 1: OAC parameter set ------------------------------------------
  # par_inv contains only OAC-type parameters (chl, a_dg_440, bb_p_550,
  # a_dg_s, bb_p_gamma, phi_f).  Shallow params (h_w, r_rs_b_*, mix_sand)
  # are appended automatically in the shallow branch below.
  par_oac <- if (sicf && !"phi_f" %in% par_inv) c(par_inv, "phi_f") else par_inv
  if (!is.null(fixed))
    par_oac <- par_oac[!par_oac %in% names(fixed)]

  # ---- Step 2: Bounds for OAC params --------------------------------------
  lo_oac <- .ALL_LOWER[par_oac]
  hi_oac <- .ALL_UPPER[par_oac]
  unknown <- par_oac[is.na(lo_oac)]
  if (length(unknown) > 0)
    stop("No default bounds for: ", paste(unknown, collapse = ", "),
         ". Supply lower/upper.")
  if (!is.null(lower)) {
    nm <- intersect(names(lower), names(lo_oac))
    for (n in nm) lo_oac[n] <- lower[[n]]
  }
  if (!is.null(upper)) {
    nm <- intersect(names(upper), names(hi_oac))
    for (n in nm) hi_oac[n] <- upper[[n]]
  }

  # ---- Step 3: OAC initial values -----------------------------------------
  emp      <- .estimate_prior_modes(rrs_df)
  init_oac <- exp(0.5 * (log(lo_oac) + log(hi_oac)))   # geometric midpoint
  known    <- names(emp)[names(emp) %in% par_oac]
  init_oac[known] <- emp[known]                          # overwrite with empirical
  init_oac <- pmax(lo_oac * 1.01, pmin(hi_oac * 0.99, init_oac))

  # ---- Step 4: par_fixed --------------------------------------------------
  par_fixed <- list(
    water_type        = water_type,
    theta_sun         = theta_sun,
    theta_view        = 0,
    a_dg_s            = 0.017,
    bb_p_gamma        = 0.5,
    lat               = lat,
    lon               = lon,
    date_time         = if (!is.null(date_time)) as.numeric(date_time) else NULL,
    sicf_model        = sicf_model,
    depth_integration = depth_integration
  )
  par_fixed[intersect(names(par_fixed), par_oac)] <- NULL
  if (!is.null(fixed)) par_fixed[names(fixed)] <- fixed

  # ---- Step 5: Shallow detection ------------------------------------------
  is_shallow <- isTRUE(force_shallow) ||
    (!is.na(depth_m) &&
     .is_optically_shallow(rrs_df, depth_m,
                           depth_threshold, optical_depth_threshold, red_green_min))

  # ==========================================================================
  # DEEP BRANCH
  # ==========================================================================
  if (!is_shallow) {
    p  <- c(par_oac, "sd")
    lb <- setNames(c(lo_oac,   .SD_LOWER), p)
    ub <- setNames(c(hi_oac,   .SD_UPPER), p)
    iv <- setNames(c(init_oac, 5e-4),      p)

    bt_prior <- make_bt_prior_adaptive(p, lb, ub, iv)
    log_prior_fn <- local({
      lm <- log(init_oac)
      ls <- pmax(0.3, pmin(abs(log(ub[par_oac]) - lm),
                            abs(lm - log(lb[par_oac]))) * 0.5)
      pn <- par_oac
      function(par) sum(dnorm(log(par[pn]), mean = lm, sd = ls, log = TRUE) - log(par[pn]))
    })

    return(list(
      shallow      = FALSE,
      h_w          = NULL,
      fwd_model    = if (sicf) "am03_sicf" else "am03",
      par_inv      = p,
      lower_b      = lb,
      upper_b      = ub,
      init_val     = iv,
      bt_prior     = bt_prior,
      log_prior_fn = log_prior_fn,
      par_fixed    = par_fixed
    ))
  }

  # ==========================================================================
  # SHALLOW BRANCH — auto-append h_w then benthic fractions
  # ==========================================================================
  sel_cls <- getOption("SABER.selected_classes", default = character(0))
  n_cl    <- length(sel_cls)
  if (n_cl == 0)
    stop("Optically shallow pixel but no benthic classes loaded. ",
         "Call select_benthic_classes() first.")

  fwd_model <- if (sicf) "am03_sicf" else if (n_cl == 2) "am03_2b" else "am03"

  # Working vectors: start from OAC and extend one param at a time.
  # Using explicit c() extension ensures positional alignment: names(lb) == p,
  # so make_bt_prior_shallow_*() can subscript by position without NA risk.
  lo   <- lo_oac
  hi   <- hi_oac
  iv_v <- init_oac
  pars <- par_oac

  # -- h_w: append unless pinned via fixed -----------------------------------
  if (!"h_w" %in% names(par_fixed)) {
    hw_lo   <- if (!is.null(lower) && "h_w" %in% names(lower)) lower[["h_w"]] else .ALL_LOWER[["h_w"]]
    hw_hi   <- if (!is.null(upper) && "h_w" %in% names(upper)) upper[["h_w"]] else .ALL_UPPER[["h_w"]]
    hw_init <- if (!is.na(depth_m))
                 max(hw_lo * 1.01, min(hw_hi * 0.99, depth_m))
               else
                 exp(0.5 * (log(hw_lo) + log(hw_hi)))
    pars <- c(pars, "h_w")
    lo   <- c(lo,   h_w = hw_lo)
    hi   <- c(hi,   h_w = hw_hi)
    iv_v <- c(iv_v, h_w = hw_init)
  }

  # Capture OAC + h_w names/bounds/inits NOW, before benthic params are
  # appended.  These are the continuous positive parameters that get lognormal
  # priors in the gradient log_prior_fn.  Without this regularisation the
  # gradient optimiser has no constraint on OAC params and falls into
  # equifinality traps (e.g. chl→0, a_dg→∞ trades off perfectly at high SNR).
  oac_hw_pars <- pars
  oac_hw_lo   <- lo
  oac_hw_hi   <- hi
  oac_hw_iv   <- iv_v

  # -- Benthic fractions: append unless pinned via fixed --------------------
  if (n_cl == 2) {
    if (!"mix_sand" %in% names(par_fixed)) {
      mx_lo <- if (!is.null(lower) && "mix_sand" %in% names(lower)) lower[["mix_sand"]] else 0.001
      mx_hi <- if (!is.null(upper) && "mix_sand" %in% names(upper)) upper[["mix_sand"]] else 0.999
      pars  <- c(pars, "mix_sand")
      lo    <- c(lo,   mix_sand = mx_lo)
      hi    <- c(hi,   mix_sand = mx_hi)
      iv_v  <- c(iv_v, mix_sand = 0.5)
    }
    p  <- c(pars, "sd")
    lb <- setNames(c(lo,   .SD_LOWER), p)
    ub <- setNames(c(hi,   .SD_UPPER), p)
    iv <- setNames(c(iv_v, 5e-4),      p)

    bt_prior     <- make_bt_prior_shallow_2class_adaptive(p, lb, ub, iv)
    log_prior_fn <- local({
      pn   <- oac_hw_pars
      lm_v <- log(oac_hw_iv)
      ls_v <- pmax(0.3, pmin(abs(log(oac_hw_hi) - lm_v),
                              abs(lm_v - log(oac_hw_lo))) * 0.5)
      function(par) {
        oac_lp <- sum(dnorm(log(par[pn]), mean = lm_v, sd = ls_v, log = TRUE) -
                        log(par[pn]))
        mix_lp <- if ("mix_sand" %in% names(par))
                    dbeta(par[["mix_sand"]], 2, 2, log = TRUE) else 0
        oac_lp + mix_lp
      }
    })

  } else {
    free_cls <- sel_cls[!sel_cls %in% names(par_fixed)]
    for (nm in free_cls) {
      nm_lo <- if (!is.null(lower) && nm %in% names(lower)) lower[[nm]] else 0.001
      nm_hi <- if (!is.null(upper) && nm %in% names(upper)) upper[[nm]] else 1.0
      pars  <- c(pars, nm)
      lo    <- c(lo,   setNames(nm_lo,       nm))
      hi    <- c(hi,   setNames(nm_hi,       nm))
      iv_v  <- c(iv_v, setNames(1 / n_cl,   nm))
    }
    p  <- c(pars, "sd")
    lb <- setNames(c(lo,   .SD_LOWER), p)
    ub <- setNames(c(hi,   .SD_UPPER), p)
    iv <- setNames(c(iv_v, 5e-4),      p)

    bt_prior     <- make_bt_prior_shallow_nclass(p, lb, ub)
    log_prior_fn <- local({
      pn   <- oac_hw_pars
      lm_v <- log(oac_hw_iv)
      ls_v <- pmax(0.3, pmin(abs(log(oac_hw_hi) - lm_v),
                              abs(lm_v - log(oac_hw_lo))) * 0.5)
      function(par) {
        oac_lp <- sum(dnorm(log(par[pn]), mean = lm_v, sd = ls_v, log = TRUE) -
                        log(par[pn]))
        b      <- par[grepl("^r_rs_b_", names(par))]
        dir_lp <- if (length(b) > 0)
                    (2 - 1) * sum(log(b + 1e-9)) - 50 * (sum(b) - 1)^2 else 0
        oac_lp + dir_lp
      }
    })
  }

  list(
    shallow      = TRUE,
    h_w          = depth_m,
    fwd_model    = fwd_model,
    par_inv      = p,
    lower_b      = lb,
    upper_b      = ub,
    init_val     = iv,
    bt_prior     = bt_prior,
    log_prior_fn = log_prior_fn,
    par_fixed    = par_fixed
  )
}

# OAC parameter names recognised by iop_from_oac (C implementation)
.IOP_OAC_NAMES <- c("chl", "a_dg_440", "a_nap_440", "bb_p_550",
                     "a_dg_s", "a_nap_s", "bb_p_gamma")

# ============================================================================
#' Compute L2 optical products from inversion results
#'
#' Reconstruct Rrs and compute derived L2 products from inversion results
#'
#' Re-runs the forward model at the retrieved parameters, computes
#' goodness-of-fit metrics (RMSE, MAPE, R\eqn{^2}), and optionally derives
#' IOPs, \eqn{K_d}, PAR, and diagnostic plots.
#'
#' @param rrs_df       Data frame with columns \code{wavelength} \[nm\] and
#'   \code{rrs_0m} \[sr\eqn{^{-1}}\].
#' @param cfg          List returned by [make_inversion_params()].
#' @param result_mcmc  Named numeric vector from [inverse_mcmc()], or
#'   \code{NULL}.
#' @param result_grad  Named numeric vector from [inverse_gradient()], or
#'   \code{NULL}.
#' @param compute_iop  Return full IOP spectra (\code{a}, \code{bb},
#'   \code{a_phy}, \code{a_dg}, \code{bb_p})?  Default \code{TRUE}.
#' @param compute_kd   Return diffuse attenuation
#'   \eqn{K_d = (a + b_b)/0.8}?  Default \code{TRUE}.
#' @param compute_par  Compute PAR via the Gregg-Carder irradiance model?
#'   Requires \code{lat}, \code{lon}, \code{date_time} in
#'   \code{cfg$par_fixed}.  Default \code{FALSE}.
#' @param plot_rrs  Build a ggplot2 Rrs reconstruction figure?  Default
#'   \code{TRUE}.
#' @param plot_iop  Build a 4-panel IOP + \eqn{K_d} figure?  Default
#'   \code{TRUE}.
#'
#' @return Named list with elements \code{wavelength}, \code{rrs_obs},
#'   \code{rrs_recon}, \code{gof}, \code{par_retrieved}; optionally
#'   \code{iop}, \code{Kd}, \code{par}, and \code{plots}.
#'
#' @export
compute_l2_products <- function(
    rrs_df,
    cfg,
    result_mcmc  = NULL,
    result_grad  = NULL,
    compute_iop  = TRUE,
    compute_kd   = TRUE,
    compute_par  = FALSE,
    plot_rrs     = TRUE,
    plot_iop     = TRUE
) {

  # --- Validate ---------------------------------------------------------------
  if (is.null(result_mcmc) && is.null(result_grad))
    stop("At least one of result_mcmc or result_grad must be non-NULL.")

  if (compute_par) {
    need <- c("lat", "lon", "date_time")
    missing_geo <- need[!need %in% names(cfg$par_fixed) |
                          vapply(cfg$par_fixed[need], is.null, logical(1))]
    if (length(missing_geo) > 0)
      stop("compute_par = TRUE requires cfg$par_fixed to contain: ",
           paste(missing_geo, collapse = ", "))
  }

  # --- Observational wavelength grid -----------------------------------------
  wl      <- sort(rrs_df$wavelength)
  rrs_obs <- rrs_df$rrs_0m[order(rrs_df$wavelength)]

  # --- Process each result ----------------------------------------------------
  .proc <- function(result, label) {
    if (is.null(result)) return(NULL)
    tryCatch(
      .compute_l2_one(result, wl, rrs_obs, cfg, compute_par),
      error = function(e) {
        warning("compute_l2_products: failed for '", label, "': ", conditionMessage(e))
        NULL
      })
  }

  out_mc <- .proc(result_mcmc, "mcmc")
  out_gr <- .proc(result_grad, "grad")

  # --- Assemble output --------------------------------------------------------
  .pluck <- function(field) {
    list(mcmc = if (!is.null(out_mc)) out_mc[[field]] else NULL,
         grad = if (!is.null(out_gr)) out_gr[[field]] else NULL)
  }

  res <- list(
    wavelength    = wl,
    rrs_obs       = rrs_obs,
    rrs_recon     = .pluck("rrs_recon"),
    gof           = .pluck("gof"),
    par_retrieved = .pluck("par_retrieved")
  )

  if (compute_iop) res$iop <- .pluck("iop")
  if (compute_kd)  res$Kd  <- .pluck("Kd")
  if (compute_par) res$par <- .pluck("par")

  # --- Build ggplot2 figures ------------------------------------------------
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plots <- list()
    if (plot_rrs)
      plots$rrs_recon <- .make_rrs_plot(
        wl          = wl,
        rrs_obs     = rrs_obs,
        rrs_mc      = if (!is.null(out_mc)) out_mc$rrs_recon else NULL,
        rrs_gr      = if (!is.null(out_gr)) out_gr$rrs_recon else NULL,
        result_mcmc = result_mcmc,
        result_grad = result_grad,
        gof_mc      = if (!is.null(out_mc)) out_mc$gof else NULL,
        gof_gr      = if (!is.null(out_gr)) out_gr$gof else NULL
      )
    if (plot_iop)
      plots$iop <- .make_iop_plot(
        wl     = wl,
        iop_mc = if (!is.null(out_mc)) out_mc$iop else NULL,
        iop_gr = if (!is.null(out_gr)) out_gr$iop else NULL,
        kd_mc  = if (!is.null(out_mc)) out_mc$Kd  else NULL,
        kd_gr  = if (!is.null(out_gr)) out_gr$Kd  else NULL
      )
    if (length(plots) > 0) res$plots <- plots
  }

  res
}

# Internal: compute products for a single result vector ---------------------
.compute_l2_one <- function(result, wl, rrs_obs, cfg, compute_par) {

  pf <- cfg$par_fixed

  # --- 1. Strip posterior SD entries ----------------------------------------
  nms  <- names(result)
  keep <- !grepl("_sd$", nms) & nms != "sd"
  est  <- result[keep]

  # --- 2. Build OAC parameter vector ----------------------------------------
  # Merge: retrieved values take precedence; fill gaps from par_fixed.
  # Start as a list to avoid type coercion when par_fixed contains strings.
  oac_pool <- as.list(est)
  for (nm in names(pf)) {
    if (!nm %in% names(oac_pool) && !is.null(pf[[nm]]))
      oac_pool[[nm]] <- pf[[nm]]
  }
  oac_nms <- intersect(names(oac_pool), .IOP_OAC_NAMES)
  oac_par <- setNames(as.numeric(unlist(oac_pool[oac_nms])), oac_nms)
  oac_par <- oac_par[!is.na(oac_par)]

  # --- 3. Compute IOPs -------------------------------------------------------
  iop <- iop_from_oac(wl, oac_par)

  # --- 4. Benthic reflectance (shallow only) ---------------------------------
  r_b <- NULL
  h_w <- NULL
  if (isTRUE(cfg$shallow)) {
    sel_cls <- getOption("SABER.selected_classes", default = character(0))
    r_b     <- .resolve_benthic_fractions(est, sel_cls)
    h_w     <- if ("h_w" %in% names(est)) as.numeric(est[["h_w"]]) else cfg$h_w
  }

  # --- 5. Forward model parameters from par_fixed ---------------------------
  water_type <- as.integer(pf[["water_type"]] %||% 2L)
  theta_sun  <- as.numeric(pf[["theta_sun"]]  %||% 30)
  theta_view <- as.numeric(pf[["theta_view"]] %||% 0)

  # --- 6. Reconstruct Rrs ---------------------------------------------------
  rrs_recon <- if (cfg$fwd_model == "am03_sicf") {

    chl_val    <- as.numeric(oac_pool[["chl"]]      %||% 1)
    adg443_val <- as.numeric(oac_pool[["a_dg_440"]] %||% 0.05)  # 440 ≈ 443
    phi_f_val  <- as.numeric(
      if ("phi_f" %in% names(est)) est[["phi_f"]]
      else pf[["phi_f"]] %||% 0.02)

    forward_am03_sicf(
      wavelength       = wl,
      iop              = iop,
      water_type       = water_type,
      theta_sun        = theta_sun,
      theta_view       = theta_view,
      h_w              = h_w,
      r_b              = r_b,
      chl              = chl_val,
      a_dg_443         = adg443_val,
      phi_f            = phi_f_val,
      include_sicf     = TRUE,
      sicf_model       = pf[["sicf_model"]]       %||% "analytical",
      depth_integration = isTRUE(pf[["depth_integration"]]),
      lat              = as.numeric(pf[["lat"]]       %||% 49),
      lon              = as.numeric(pf[["lon"]]       %||% -68),
      date_time        = if (!is.null(pf[["date_time"]]))
                           as.POSIXct(pf[["date_time"]], origin = "1970-01-01", tz = "UTC")
                         else Sys.time(),
      return_components = FALSE
    )

  } else {
    forward_am03(
      wavelength = wl,
      iop        = iop,
      water_type = water_type,
      theta_sun  = theta_sun,
      theta_view = theta_view,
      h_w        = h_w,
      r_b        = r_b
    )
  }

  # --- 7. Goodness-of-fit ---------------------------------------------------
  .gof <- function(mod) {
    valid <- is.finite(rrs_obs) & is.finite(mod) & rrs_obs > 0
    if (sum(valid) < 2) return(list(rmse = NA_real_, mape = NA_real_, r2 = NA_real_))
    o <- rrs_obs[valid]; m <- mod[valid]
    r  <- o - m
    list(
      rmse = sqrt(mean(r^2)),
      mape = 100 * mean(abs(r) / o),
      r2   = 1 - sum(r^2) / sum((o - mean(o))^2)
    )
  }
  gof <- .gof(rrs_recon)

  # --- 8. Derived products (always computed; caller controls public exposure)
  w    <- pure_water_iop(wl)
  bb_p <- pmax(0, iop$bb - w$bb_w)
  out <- list(
    rrs_recon     = rrs_recon,
    gof           = gof,
    par_retrieved = oac_par,
    iop = list(
      a     = iop$a,
      bb    = iop$bb,
      a_phy = iop$a_phy,
      a_dg  = iop$a_dg,
      bb_p  = bb_p
    ),
    Kd = (iop$a + iop$bb) / 0.8
  )

  if (compute_par) {
    comp <- forward_am03_sicf(
      wavelength       = wl,
      iop              = iop,
      water_type       = water_type,
      theta_sun        = theta_sun,
      theta_view       = theta_view,
      h_w              = h_w,
      r_b              = r_b,
      chl              = as.numeric(oac_pool[["chl"]]      %||% 1),
      a_dg_443         = as.numeric(oac_pool[["a_dg_440"]] %||% 0.05),
      phi_f            = as.numeric(
        if ("phi_f" %in% names(est)) est[["phi_f"]]
        else pf[["phi_f"]] %||% 0.02),
      include_sicf     = isTRUE(cfg$fwd_model == "am03_sicf"),
      sicf_model       = pf[["sicf_model"]]       %||% "analytical",
      depth_integration = isTRUE(pf[["depth_integration"]]),
      lat              = as.numeric(pf[["lat"]]       %||% 49),
      lon              = as.numeric(pf[["lon"]]       %||% -68),
      date_time        = if (!is.null(pf[["date_time"]]))
                           as.POSIXct(pf[["date_time"]], origin = "1970-01-01", tz = "UTC")
                         else Sys.time(),
      return_components = TRUE
    )
    out$par <- comp$PAR
  }

  out
}

# NULL-coalescing operator (internal, avoids rlang dependency)
`%||%` <- function(x, y) if (!is.null(x)) x else y

# Internal: "name = value ± sd" annotation string for each retrieved OAC param
.format_retrieved_params <- function(result, label = NULL) {
  if (is.null(result)) return("")
  nms <- names(result)
  oac <- nms[!grepl("_sd$", nms) & nms != "sd" &
               !grepl("^r_rs_b_|^mix_sand$|^h_w$", nms)]
  if (length(oac) == 0) return("")
  lines <- vapply(oac, function(p) {
    val  <- result[[p]]
    sdnm <- paste0(p, "_sd")
    if (sdnm %in% nms)
      sprintf("  %s = %.3g \u00b1 %.2g", p, val, result[[sdnm]])
    else
      sprintf("  %s = %.3g", p, val)
  }, character(1))
  prefix <- if (!is.null(label)) paste0(label, "\n") else ""
  paste0(prefix, paste(lines, collapse = "\n"))
}

# Internal: forward-reconstruction ggplot (Rrs observed vs MCMC vs Gradient) -
.make_rrs_plot <- function(wl, rrs_obs, rrs_mc, rrs_gr,
                            result_mcmc, result_grad,
                            gof_mc, gof_gr) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)

  sources  <- "Observed"
  rrs_list <- list(rrs_obs)
  if (!is.null(rrs_mc)) {
    sources  <- c(sources, "MCMC")
    rrs_list <- c(rrs_list, list(rrs_mc))
  }
  if (!is.null(rrs_gr)) {
    sources  <- c(sources, "Gradient")
    rrs_list <- c(rrs_list, list(rrs_gr))
  }

  gof_df <- data.frame(
    wavelength = rep(wl, length(sources)),
    rrs        = unlist(rrs_list),
    source     = rep(sources, each = length(wl)),
    stringsAsFactors = FALSE
  )
  ymax <- max(gof_df$rrs, na.rm = TRUE) * 1.15

  # GoF statistics block (top-left)
  gof_parts <- character(0)
  if (!is.null(gof_mc))
    gof_parts <- c(gof_parts, sprintf(
      "MCMC  RMSE=%.2e  MAPE=%.1f%%  R\u00b2=%.3f",
      gof_mc$rmse, gof_mc$mape, gof_mc$r2))
  if (!is.null(gof_gr))
    gof_parts <- c(gof_parts, sprintf(
      "Grad  RMSE=%.2e  MAPE=%.1f%%  R\u00b2=%.3f",
      gof_gr$rmse, gof_gr$mape, gof_gr$r2))
  stats_label <- paste(gof_parts, collapse = "\n")

  # Retrieved parameter block (top-right, with ±SD)
  par_parts <- c(
    .format_retrieved_params(result_mcmc, "MCMC"),
    .format_retrieved_params(result_grad,  "Grad")
  )
  par_label <- paste(par_parts[nchar(par_parts) > 0], collapse = "\n\n")

  col_map <- c("Observed" = "black", "MCMC" = "steelblue", "Gradient" = "firebrick")
  lty_map <- c("Observed" = "dashed", "MCMC" = "solid",    "Gradient" = "solid")

  ggplot2::ggplot(gof_df,
                  ggplot2::aes(x = wavelength, y = rrs,
                               color = source, linetype = source)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.3) +
    ggplot2::scale_color_manual(values = col_map) +
    ggplot2::scale_linetype_manual(values = lty_map) +
    ggplot2::annotate("text",
                      x = min(wl) + 5, y = ymax * 0.98,
                      label = stats_label,
                      hjust = 0, vjust = 1,
                      size = 3.5, family = "mono", color = "grey20") +
    ggplot2::annotate("text",
                      x = max(wl) - 5, y = ymax * 0.98,
                      label = par_label,
                      hjust = 1, vjust = 1,
                      size = 3.0, family = "mono", color = "grey30") +
    ggplot2::labs(
      x        = expression(paste("Wavelength (", lambda, ") [nm]")),
      y        = expression(paste(italic(R)[rs], " [sr"^{-1}, "]")),
      color    = NULL,
      linetype = NULL) +
    ggplot2::coord_cartesian(xlim = range(wl), ylim = c(0, ymax)) +
    ggplot2::scale_x_continuous(breaks = seq(400, 700, 50)) +
    ggplot2::scale_y_continuous(
      labels = function(x) format(x, scientific = FALSE, digits = 3)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text         = ggplot2::element_text(size = 13, color = "black"),
      axis.title        = ggplot2::element_text(size = 15),
      axis.ticks.length = ggplot2::unit(0.2, "cm"),
      legend.position   = c(0.15, 0.7),
      legend.title      = ggplot2::element_blank(),
      legend.text       = ggplot2::element_text(size = 13),
      legend.background = ggplot2::element_rect(fill = NA),
      panel.grid.major  = ggplot2::element_line(colour = "grey70", linewidth = 0.4,
                                                 linetype = "dotted"),
      panel.grid.minor  = ggplot2::element_line(colour = "grey85", linewidth = 0.2),
      panel.border      = ggplot2::element_rect(colour = "black", fill = NA,
                                                 linewidth = 1.2),
      plot.margin       = ggplot2::unit(c(0.5, 0.8, 0.5, 0.5), "cm"))
}

# Internal: 4-panel component IOP+ ggplot (a_phy / a_dg / bb_p / Kd) -------------
.make_iop_plot <- function(wl, iop_mc, iop_gr, kd_mc = NULL, kd_gr = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  
  .rows <- function(iop_list, kd_vec, label) {
    if (is.null(iop_list) && is.null(kd_vec)) return(NULL)
    components <- c("a_phy", "a_dg", "bb_p")
    values     <- c(iop_list$a_phy, iop_list$a_dg, iop_list$bb_p)
    if (!is.null(kd_vec)) {
      components <- c(components, "Kd")
      values     <- c(values, kd_vec)
    }
    data.frame(
      wavelength = rep(wl, length(components)),
      value      = values,
      component  = rep(components, each = length(wl)),
      source     = label,
      stringsAsFactors = FALSE
    )
  }
  df <- rbind(.rows(iop_mc, kd_mc, "MCMC"), .rows(iop_gr, kd_gr, "Gradient"))
  if (is.null(df)) return(NULL)

  comp_labels <- c(
    "a_phy" = "a_phy  [m\u207b\u00b9]",
    "a_dg"  = "a_dg   [m\u207b\u00b9]",
    "bb_p"  = "bb_p   [m\u207b\u00b9]",
    "Kd"    = "Kd     [m\u207b\u00b9]"
  )
  df$component <- factor(df$component, levels = c("a_phy", "a_dg", "bb_p", "Kd"))

  col_map <- c("MCMC" = "steelblue", "Gradient" = "firebrick")

  ggplot2::ggplot(df, ggplot2::aes(x = wavelength, y = value,
                                   color = source, linetype = source)) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::geom_point(size = 1.3) +
    ggplot2::facet_wrap(~ component, scales = "free_y", ncol = 1,
                        labeller = ggplot2::as_labeller(comp_labels)) +
    ggplot2::scale_color_manual(values = col_map) +
    ggplot2::scale_linetype_manual(
      values = c("MCMC" = "solid", "Gradient" = "solid")) +
    ggplot2::labs(
      x        = expression(paste("Wavelength (", lambda, ") [nm]")),
      y        = NULL,
      color    = NULL,
      linetype = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(400, 700, 50)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text         = ggplot2::element_text(size = 11, color = "black"),
      axis.title.x      = ggplot2::element_text(size = 13),
      strip.text        = ggplot2::element_text(size = 12, face = "bold"),
      legend.position   = "bottom",
      legend.title      = ggplot2::element_blank(),
      legend.text       = ggplot2::element_text(size = 12),
      panel.grid.major  = ggplot2::element_line(colour = "grey70", linewidth = 0.4,
                                                 linetype = "dotted"),
      panel.grid.minor  = ggplot2::element_line(colour = "grey85", linewidth = 0.2),
      panel.border      = ggplot2::element_rect(colour = "black", fill = NA,
                                                 linewidth = 1.2),
      plot.margin       = ggplot2::unit(c(0.5, 0.8, 0.3, 0.5), "cm"))
}

# ============================================================================
#' Create per-band spectral weights for inversion
#'
#' Builds a named numeric weight vector aligned to a wavelength grid.  Higher
#' weights tighten the effective noise on that band:
#' \deqn{\sigma_{\mathrm{eff},\lambda} = \sigma / \sqrt{w_\lambda}}
#' Three spectral regions (blue, green, red) receive flat weights with cosine
#' tapers at region boundaries.  Set \code{weight_sicf} to add a two-Gaussian
#' bump (peaks 685 nm and 730 nm) over the SICF window.  Supply
#' \code{custom_weights} to bypass the region logic entirely.
#'
#' @param wavelength  Numeric vector of wavelengths [nm].
#' @param weight_blue,weight_green,weight_red  Weight scalars for the blue,
#'   green, and red regions.  Default \code{1}.
#' @param weight_sicf  Either \code{FALSE} (off, default) or a numeric peak
#'   weight for the SICF window (670–750 nm, peaks at 685 and 730 nm).  Must
#'   be \eqn{\ge} \code{weight_red}.
#' @param blue_range,green_range,red_range  Two-element [nm] region boundaries.
#'   Defaults: \code{c(400, 500)}, \code{c(500, 600)}, \code{c(600, 700)}.
#' @param taper_nm    Cosine taper half-width [nm] at region boundaries.
#'   Default \code{10}.  Set to \code{0} for step-function weights.
#' @param custom_exclude_ranges  List of \code{c(lo, hi)} [nm] vectors.
#'   Matched bands are set to \code{1e-4} (effectively excluded).  Applied
#'   after all other rules.
#' @param custom_weights  Named numeric vector (wavelength as character name,
#'   one decimal place, e.g. \code{"412.0"}).  Mutually exclusive with the
#'   colour weight arguments.
#'
#' @return Named numeric vector of length \code{length(wavelength)}.  Pass
#'   directly to \code{spectral_weights} in [inverse_mcmc()] or
#'   [inverse_gradient()].
#'
#' @examples
#' wl <- seq(400, 750, by = 0.5)
#' sw <- create_spectral_weights(wl,
#'   weight_blue = 5, weight_green = 2, weight_red = 3,
#'   weight_sicf = 8, taper_nm = 10)
#'
#' @export
create_spectral_weights <- function(
    wavelength,
    weight_blue   = 1,
    weight_green  = 1,
    weight_red    = 1,
    weight_sicf   = FALSE,
    blue_range    = c(400, 500),
    green_range   = c(500, 600),
    red_range     = c(600, 700),
    taper_nm      = 10,
    custom_exclude_ranges = NULL,
    custom_weights        = NULL
) {
  stopifnot(is.numeric(wavelength), length(wavelength) >= 1)
  stopifnot(is.numeric(taper_nm), length(taper_nm) == 1, taper_nm >= 0)

  wl <- wavelength

  # ---------- Guard: custom_weights is mutually exclusive with color args --
  if (!is.null(custom_weights)) {
    if (!missing(weight_blue) || !missing(weight_green) ||
        !missing(weight_red)  || !missing(weight_sicf))
      stop(paste(
        "custom_weights is mutually exclusive with weight_blue, weight_green,",
        "weight_red, and weight_sicf. Set only one or the other."))
    if (!is.numeric(custom_weights) || is.null(names(custom_weights)))
      stop("custom_weights must be a named numeric vector.")

    # Build output from custom_weights; default unmatched wavelengths to 1.
    wl_nms <- format(round(wl, 1), nsmall = 1)
    w <- setNames(rep(1, length(wl)), wl_nms)
    matched <- intersect(wl_nms, names(custom_weights))
    w[matched] <- custom_weights[matched]

  } else {

    # ------ Validate color weight scalars ------------------------------------
    stopifnot(is.numeric(weight_blue),  length(weight_blue)  == 1, weight_blue  >= 0)
    stopifnot(is.numeric(weight_green), length(weight_green) == 1, weight_green >= 0)
    stopifnot(is.numeric(weight_red),   length(weight_red)   == 1, weight_red   >= 0)
    if (!isFALSE(weight_sicf)) {
      if (!is.numeric(weight_sicf) || length(weight_sicf) != 1 || weight_sicf < 0)
        stop("weight_sicf must be FALSE or a non-negative numeric scalar.")
      if (weight_sicf < weight_red)
        stop("weight_sicf must be >= weight_red (the SICF peaks sit above the red baseline).")
    }

    # ------ Step 1: flat per-region weights ----------------------------------
    # Bands outside all three regions keep the neutral weight of 1.
    w <- rep(1, length(wl))

    in_blue  <- wl >= blue_range[1]  & wl <= blue_range[2]
    in_green <- wl >= green_range[1] & wl <= green_range[2]
    in_red   <- wl >= red_range[1]   & wl <= red_range[2]

    w[in_blue]  <- weight_blue
    w[in_green] <- weight_green
    w[in_red]   <- weight_red

    # ------ Step 2: cosine taper at all finite region boundaries -------------
    if (taper_nm > 0) {
      boundaries <- unique(c(
        if (is.finite(blue_range[2]))  blue_range[2],
        if (is.finite(green_range[2])) green_range[2],
        if (is.finite(red_range[2]))   red_range[2]
      ))

      for (b in boundaries) {
        left_w  <- w[which.min(abs(wl - (b - taper_nm)))]
        right_w <- w[which.min(abs(wl - (b + taper_nm)))]
        idx_tap <- which(wl >= (b - taper_nm) & wl <= (b + taper_nm))
        if (length(idx_tap) == 0) next
        t     <- (wl[idx_tap] - (b - taper_nm)) / (2 * taper_nm)
        alpha <- (1 - cos(pi * t)) / 2   # 0 → 1
        w[idx_tap] <- left_w * (1 - alpha) + right_w * alpha
      }
    }

    # ------ Step 3: SICF two-Gaussian bump (670–750 nm) ----------------------
    # Profile = G(685, σ=7) + 0.65 · G(730, σ=12), then normalised to [0,1].
    # Final weight in the SICF window: weight_red + (weight_sicf - weight_red)
    # × normalised_profile — floor is weight_red, peak is weight_sicf.
    if (!isFALSE(weight_sicf) && is.numeric(weight_sicf)) {
      sicf_lo <- 670;  sicf_hi <- 750
      idx_s   <- which(wl >= sicf_lo & wl <= sicf_hi)
      if (length(idx_s) > 0) {
        g1    <- exp(-0.5 * ((wl[idx_s] - 685) / 7 ) ^ 2)
        g2    <- 0.65 * exp(-0.5 * ((wl[idx_s] - 730) / 12) ^ 2)
        shape <- g1 + g2
        shape <- shape / max(shape)     # normalise peak to 1
        w[idx_s] <- weight_red + (weight_sicf - weight_red) * shape
      }
    }
  }

  # ---------- Step 4: apply exclusion ranges (near-zero epsilon) -----------
  EXCLUDE_EPS <- 1e-4
  if (!is.null(custom_exclude_ranges)) {
    if (!is.list(custom_exclude_ranges))
      stop("custom_exclude_ranges must be a list of c(lo, hi) numeric vectors.")
    for (rng in custom_exclude_ranges) {
      if (length(rng) != 2 || !is.numeric(rng))
        stop("Each custom_exclude_ranges entry must be a length-2 numeric vector c(lo, hi).")
      w[wl >= rng[1] & wl <= rng[2]] <- EXCLUDE_EPS
    }
  }

  names(w) <- format(round(wl, 1), nsmall = 1)
  w
}
