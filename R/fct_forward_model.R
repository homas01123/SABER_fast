#' Albert & Mobley (2003) Forward Model
#'
#' Computes Rrs below surface using analytical model from Albert & Mobley (2003)
#'
#' @param wavelength vector of wavelengths \[nm\]
#' @param iop list with elements `a` and `bb`, same length as wavelength
#' @param water_type either 1 or 2 (default = 2)
#' @param theta_sun sun zenith angle \[degrees\]
#' @param theta_view sensor view angle \[degrees\]
#' @param h_w optional: water depth \[m\] (enables shallow mode)
#' @param r_b optional: bottom reflectance vector (same length as wavelength)
#'
#' @return numeric vector of subsurface Rrs
#'
#' @references Albert, A. and Mobley, C.D. (2003) ‘An analytical model for subsurface irradiance and remote sensing reflectance in deep and shallow case-2 waters’, Optics Express, 11(22), pp. 2873–2890. Available at: https://doi.org/10.1364/OE.11.002873.
#'
#' @export
forward_am03 <- function(wavelength, iop, water_type = 2,
                         theta_sun, theta_view,
                         h_w = NULL, r_b = NULL) {
  stopifnot(is.numeric(wavelength))
  stopifnot(is.list(iop), all(c("a", "bb") %in% names(iop)))
  stopifnot(is.numeric(iop$a), is.numeric(iop$bb))
  stopifnot(
    length(iop$a) == length(wavelength),
    length(iop$bb) == length(wavelength)
  )

  .Call(
    "c_forward_am03",
    as.numeric(wavelength),
    list(a = iop$a, bb = iop$bb),
    as.integer(water_type),
    as.numeric(theta_sun),
    as.numeric(theta_view),
    h_w, r_b
  )
}

#' Function to retrieve Rrs_b
#' `Rrs_Bottom <- (obs_rrs - Rrs_below_deep * (1 - Ars1 * exp(-zB * (Kd + kuW)))) / (Ars2 * exp(-zB * (Kd + kuB)))`
#'
#' @export
retrieve_r_rs_b_am03 <- function(
    wavelength,
    iop,
    r_rs,
    water_type = 2,
    theta_sun,
    theta_view,
    h_w
    ) {
  stopifnot(is.numeric(wavelength))
  stopifnot(is.list(iop), all(c("a", "bb") %in% names(iop)))
  stopifnot(is.numeric(iop$a), is.numeric(iop$bb))
  stopifnot(
    length(iop$a) == length(wavelength),
    length(iop$bb) == length(wavelength)
  )

  .Call(
    "c_retrieve_r_rs_b_am03",
    as.numeric(wavelength),
    list(a = iop$a, bb = iop$bb),
    as.numeric(r_rs),
    as.integer(water_type),
    as.numeric(theta_sun),
    as.numeric(theta_view),
    h_w
  )
}

# forward_am03 <- function(
#     wavelength,
#     iop,
#     water_type = 2,
#     theta_view,
#     theta_sun,
#     optically_shallow = F,
#     h_w = NULL,
#     r_b = NULL,
#     verbose = F) {
#   # Define IOPs  ------------------------------------------------------------
#
#   # extinction coeff. [1/m]
#   ext <- iop$a + iop$bb
#
#   # single back scattering albedo
#   omega_b <- iop$bb / ext
#
#   # RT model ----------------------------------------------------------------
#
#   geometry <- snell_law(theta_view = theta_view, theta_sun = theta_sun)
#   sun_w <- geometry$sun_w
#   view_w <- geometry$view_w
#
#   ## Remote Sensing Reflectance below the water surface
#   if (water_type == "1") {
#     f_rs <- 0.095 # [1/sr]
#   } else if (water_type == "2") {
#     f_rs <- 0.0512 * (1 + (4.6659 * omega_b) +
#       (-7.8387 * (omega_b^2)) +
#       (5.4571 * (omega_b^3))
#     ) *
#       (1 + (0.1098 / cos(sun_w))) *
#       (1 + (0.4021 / cos(view_w))) # [1/sr]
#   }
#
#   rrs_0m_deep <- f_rs * omega_b # [1/sr]
#
#   if (optically_shallow) {
#     # Attenuation Coefficients
#     if (water_type == 1) {
#       k0 <- 1.0395 # case 1
#     } else if (water_type == 2) {
#       k0 <- 1.0546 # case 2
#     }
#
#     Kd <- k0 * (ext / cos(sun_w))
#     kuW <- (ext / cos(view_w)) * ((1 + omega_b)^3.5421) * (1 - (0.2786 / cos(sun_w)))
#     kuB <- (ext / cos(view_w)) * ((1 + omega_b)^2.2658) * (1 - (0.0577 / cos(sun_w)))
#
#     # Final calculation for shallow Rrs
#     Ars1 <- 1.1576
#     Ars2 <- 1.0389 # Parametric coeffs for shallow water
#
#     rrs_0m_shallow <- rrs_0m_deep *
#       (1 - (Ars1 * exp(-h_w * (Kd + kuW)))) + Ars2 *
#         r_b * exp(-h_w * (Kd + kuB))
#
#     rrs_0m <- rrs_0m_shallow
#   } else {
#     # Optically deep water
#     rrs_0m <- rrs_0m_deep
#   }
#
#   # rrs_0m <- tibble(
#   #   wavelength,
#   #   rrs_0m
#   # )
#
#   return(
#     rrs_0m
#   )
# }

#' Lee 1998 forward model
#' Compute Rrs given IOP
#'
#' @author Soham Mukherjee
#'
#' @param iop a tibble with column {wavelength, a, bb}, non-water component of absorption and backscattering.
#' @param optically_shallow bolean, is the water optically shallow ? If TRUE, provide `h_w`, `rb_fraction`, `rb`.
#' @param h_w water column height [m].
#' @param rrs_botom remote sensing bottom reflectance. see function `rrs_bottom_lmm`
#' @param wavelength Requested wavelength for the simulation
#'
#' @references Lee, z. et al. (1998) ‘Hyperspectral remote sensing for shallow waters. I. A semianalytical model’, Applied Optics, 37(27), pp. 6329–6338. Available at: https://doi.org/10.1364/AO.37.006329.
#'

forward_lee98 <- function(
    wavelength,
    iop,
    theta_view,
    theta_sun,
    optically_shallow,
    h_w = NULL,
    r_b = NULL,
    verbose = F) {
  # Define IOPs  ------------------------------------------------------------

  # extinction coeff. [1/m]
  ext <- iop$a + iop$bb

  # single back scattering albedo
  omega_b <- iop$bb / ext

  geometry <- snell_law(theta_view = theta_view, theta_sun = theta_sun)
  sun_w <- geometry$sun_w
  view_w <- geometry$view_w

  ## Remote Sensing Reflectance 0m
  p1 <- 0.084
  p2 <- 0.17
  k1w <- 1.03
  k2w <- 2.04
  k1b <- 1.04
  k2b <- 5.04
  q1 <- 1 # (for viewing angle=0)

  rrs_0m_deep <- q1 * (p1 + p2 * omega_b) * omega_b

  if (optically_shallow) {
    # Attenuation Coefficients
    mu_s <- cos(sun_w)
    mu_v <- cos(view_w)
    du_w <- k1w * sqrt(1 + k2w * omega_b)
    du_b <- k1b * sqrt(1 + k2b * omega_b)

    rrs_0m_shallow <- rrs_0m_deep *
      (1 - exp(-(1 / mu_s + du_w / mu_v) * ext * h_w)) +
      r_b$r_b * exp(-(1 / mu_s + du_b / mu_v) * ext * h_w)

    rrs_0m <- rrs_0m_shallow
  } else {
    # Optically deep water
    rrs_0m <- rrs_0m_deep
  }

  if (verbose) {
    plot(wavelength, rrs_0m,
      xlab = "wavelength",
      ylab = "Rrs 0m [m^-1]"
    )
  }

  return(tibble(
    wavelength,
    rrs_0m
  ))
}
