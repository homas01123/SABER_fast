#' parse_iop
#'
#' take tibbles for non water absorption and backscattering
#'
#' @author Raphael Mabit, Soham Mukherjee
#'
#' @param a tibble with columns {a, wavelength}
#' @param bb tibble with columns {bb, wavelength}
#' @param wavelength a vector with the requested interpolation wavelength
#' @param verbose Boolean, produce input and interpolated IOPs plot
#'
#' @returns A tibble with columns {wavelength, a, bb}
#'
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
