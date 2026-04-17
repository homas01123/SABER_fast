#' Gregg & Carder (1990) spectral solar irradiance model
#'
#' Reference spectral data used by the Gregg & Carder atmospheric irradiance
#' model. Columns: \code{lam} (wavelength, nm), \code{Fobar} (extraterrestrial
#' solar irradiance, W/m2/nm), \code{oza} (ozone absorption coefficient),
#' \code{ag} (oxygen/gas absorption), \code{aw} (water vapour absorption).
#'
#' @format A data frame with 581 rows and 5 columns.
#' @source Gregg, W.W. & Carder, K.L. (1990). A simple spectral solar irradiance
#'   model for cloudless maritime atmospheres. Limnology and Oceanography, 35(8),
#'   1657-1675.
"gc_data"


# ---------------------------------------------------------------------------
# Internal Gregg & Carder implementation — no external package dependency
# ---------------------------------------------------------------------------

#' Compute solar zenith angle
#'
#' @param iday Julian day (integer)
#' @param hr Time in decimal hours GMT
#' @param xlon Longitude, W negative (decimal degrees)
#' @param ylat Latitude, N positive (decimal degrees)
#' @return Solar zenith angle in degrees
#' @keywords internal
.gc_sunang <- function(iday, hr, xlon, ylat) {
  rad   <- 180 / pi
  thez  <- 360 * (iday - 1) / 365
  rthez <- thez / rad
  sdec  <- 0.396372 - 22.91327 * cos(rthez) + 4.02543 * sin(rthez) -
            0.387205 * cos(2 * rthez) + 0.051967 * sin(2 * rthez) -
            0.154527 * cos(3 * rthez) + 0.084798 * sin(3 * rthez)
  rsdec <- sdec / rad
  tc  <- 0.004297 + 0.107029 * cos(rthez) - 1.837877 * sin(rthez) -
         0.837378 * cos(2 * rthez) - 2.342824 * sin(2 * rthez)
  xha <- (hr - 12) * 15 + xlon + tc
  if (xha >  180) xha <- xha - 360
  if (xha < -180) xha <- xha + 360
  rlat <- ylat / rad
  rha  <- xha  / rad
  costmp <- sin(rlat) * sin(rsdec) + cos(rlat) * cos(rsdec) * cos(rha)
  acos(costmp) * rad
}


#' Compute Navy marine aerosol parameters
#'
#' @param rh  Relative humidity (%)
#' @param am  Air mass type (1-10)
#' @param wsm Mean wind speed (m/s)
#' @param ws  Current wind speed (m/s)
#' @param vis Visibility (km)
#' @return List: beta, alpha, wa (single-scattering albedo), asymp
#' @keywords internal
.gc_navaer <- function(rh, am, wsm, ws, vis) {
  ro   <- c(0.03, 0.24, 2.0)
  r    <- c(0.1,  1.0,  10.0)
  rlam <- 0.55
  if (rh >= 100) rh <- 99.9
  frh  <- ((2 - rh / 100) / (6 * (1 - rh / 100)))^0.333

  a    <- numeric(3)
  a[1] <- 2000 * am * am
  a[2] <- max(0.5, 5.866 * (wsm - 2.2))
  a[3] <- max(1.4e-5, 0.01527 * (ws - 2.2) * 0.05)

  # Size distribution at three radii (fully vectorised over i)
  dndr <- vapply(seq_len(3), function(n) {
    arg  <- (log(r[n] / (frh * ro)))^2
    sum(a * exp(-arg) / frh)
  }, numeric(1))

  # Least-squares fit in log-log space
  lr   <- log10(r)
  ld   <- log10(dndr)
  gama <- sum(lr * ld) / sum(lr^2)
  alpha <- -(gama + 3)

  beta  <- (3.91 / vis) * rlam^alpha
  asymp <- if (alpha > 1.2) 0.65 else if (alpha < 0) 0.82 else -0.14167 * alpha + 0.82
  w0    <- (-0.0032 * am + 0.972) * exp(3.06e-4 * rh)

  list(beta = beta, alpha = alpha, wa = w0, asymp = asymp)
}


#' Compute Fresnel surface reflectance
#'
#' @param theta Solar zenith angle (degrees)
#' @param ws    Wind speed (m/s)
#' @return List: rod (direct reflectance), ros (diffuse reflectance)
#' @keywords internal
.gc_sfcrfl <- function(theta, ws) {
  rad  <- 180 / pi
  rn   <- 1.341
  rof  <- 0; rosps <- 0.066
  if (ws > 4) {
    rof  <- if (ws <= 7) (6.2e-4 + 1.56e-3 / ws) * 1.2e3 * 2.2e-5 * ws^2 - 4e-4 else
                         (0.49e-3 + 0.065e-3 * ws) * 1.2e3 * 4.5e-5 * ws^2 - 4e-5 * ws^2
    rosps <- 0.057
  }
  if (theta < 50 || ws < 2) {
    rospd <- if (theta == 0) 0.0211 else {
      rtheta  <- theta / rad
      sintr   <- sin(rtheta) / rn
      rthetar <- asin(sintr)
      sinp    <- (sin(rtheta - rthetar) / sin(rtheta + rthetar))^2
      tanp    <- (tan(rtheta - rthetar) / tan(rtheta + rthetar))^2
      0.5 * (sinp + tanp)
    }
  } else {
    a <- 5.25e-4 * ws + 0.065
    b <- -1.67e-3 * ws + 0.074
    rospd <- a * exp(b * (theta - 60))
  }
  list(rod = rospd + rof, ros = rosps + rof)
}


#' Compute spectral atmospheric transmittance and downwelling irradiance
#'
#' @keywords internal
.gc_atmodd <- function(lam, theta, oza, ag, aw, sco3, p,
                       wv, rh, am, wsm, ws, vis, Fo) {
  rad    <- 180 / pi
  p0     <- 1013.25
  cosunz <- cos(theta / rad)
  rex    <- -1.253
  rm     <- 1 / (cosunz + 0.15 * (93.885 - theta)^rex)
  rmp    <- p / p0 * rm
  otmp   <- sqrt(cosunz^2 + 44 / 6370)
  rmo    <- (1 + 22 / 6370) / otmp

  nav  <- .gc_navaer(rh, am, wsm, ws, vis)
  eta  <- -nav$alpha
  alg  <- log(1 - nav$asymp)
  afs  <- alg * (1.459  + alg * (0.1595 + alg * 0.4129))
  bfs  <- alg * (0.0783 + alg * (-0.3824 - alg * 0.5874))
  Fa   <- 1 - 0.5 * exp((afs + bfs * cosunz) * cosunz)

  rfl  <- .gc_sfcrfl(theta, ws)

  rlam <- lam * 1e-3
  tr   <- 1 / (115.6406 * rlam^4 - 1.335 * rlam^2)
  rtra <- exp(-tr * rmp)
  to   <- oza * sco3
  otra <- exp(-to * rmo)
  ta   <- nav$beta * rlam^eta
  atra <- exp(-ta * rm)
  taa  <- exp(-(1 - nav$wa)  * ta * rm)
  tas  <- exp(-nav$wa * ta * rm)

  gtmp  <- (1 + 118.3 * ag  * rmp)^0.45
  gtra  <- exp(-1.41 * ag  * rmp / gtmp)
  wtmp  <- (1 + 20.07 * aw  * wv * rm)^0.45
  wtra  <- exp(-0.2385 * aw * wv * rm / wtmp)

  Edir <- Fo * cosunz * rtra * otra * atra * gtra * wtra
  dray <- Fo * cosunz * gtra * wtra * otra * taa * 0.5 * (1 - rtra^0.95)
  daer <- Fo * cosunz * gtra * wtra * otra * rtra^1.5 * taa * Fa * (1 - tas)
  Edif <- dray + daer

  list(Edir = Edir, Edif = Edif, Ed = Edir + Edif)
}


#' Gregg & Carder spectral solar irradiance at the sea surface (0+)
#'
#' Computes direct and diffuse downwelling spectral irradiance just above
#' the sea surface (0+) using the Gregg & Carder (1990) maritime atmosphere
#' model.  Fresnel surface reflectance is NOT applied here — call
#' \code{.gc_sfcrfl} separately and apply to convert to subsurface (0-).
#' Uses the internal \code{gc_data} spectral reference table — no external
#' package required.
#'
#' @param jday   Julian day (integer, 1-365)
#' @param rlon   Longitude, W negative (decimal degrees)
#' @param rlat   Latitude,  N positive (decimal degrees)
#' @param lam.sel Wavelength vector (nm, default 350:700)
#' @param the    Solar zenith angle (degrees); \code{-99} to compute from
#'   \code{hr}/\code{rlon}/\code{rlat}
#' @param hr     Time in decimal hours GMT; required when \code{the < -90}
#' @param Vi     Visibility (km, default 15)
#' @param am     Air mass type 1-10 (default 1)
#' @param wsm    Mean wind speed m/s (default 4)
#' @param ws     Current wind speed m/s (default 6)
#' @param pres   Atmospheric pressure mb (default 1013.25)
#' @param rh     Relative humidity % (default 80)
#' @param wv     Water vapour cm (default 1.5)
#'
#' @return List: \code{lam}, \code{Edir} (direct), \code{Edif} (diffuse),
#'   \code{Ed} = Edir + Edif (all in W/m2/nm at 0+), plus
#'   \code{sunzen_deg} (resolved solar zenith angle in degrees).
#' @keywords internal
.gc_irradiance <- function(jday, rlon, rlat, lam.sel = 350:700,
                            the = -99, hr = -99,
                            Vi = 15, am = 1, wsm = 4, ws = 6,
                            pres = 1013.25, rh = 80, wv = 1.5) {
  pi2 <- 2 * pi
  to3 <- 235 + (150 + 40 * sin(0.9865 * (jday - 30)) +
                20 * sin(3 * rlon)) * sin(1.28 * rlat)^2
  sco3 <- to3 * 1e-3

  if (the < -90) {
    if (hr < -90) stop("Specify either solar zenith angle (the) or time (hr).")
    theta <- .gc_sunang(jday, hr, rlon, rlat)
  } else {
    theta <- the
  }

  if (theta >= 90) stop("Sun is below horizon for given geometry.")

  # Interpolate spectral reference data onto requested wavelengths
  d    <- gc_data
  Fobar <- approx(d$lam, d$Fobar, lam.sel)$y
  oza   <- approx(d$lam, d$oza,   lam.sel)$y
  ag    <- approx(d$lam, d$ag,    lam.sel)$y
  aw    <- approx(d$lam, d$aw,    lam.sel)$y

  # Earth-Sun distance correction
  Fo <- Fobar * (1 + 1.67e-2 * cos(pi2 * (jday - 3) / 365))^2

  res <- .gc_atmodd(lam.sel, theta, oza, ag, aw, sco3, pres,
                    wv, rh, am, wsm, ws, Vi, Fo)

  list(lam = lam.sel, Edir = res$Edir, Edif = res$Edif,
       Ed = res$Ed, sunzen_deg = theta)
}
