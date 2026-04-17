#' iop_from_oac
#' compute Inherent Optical Proprieties from Optically Active constituents.
#'
#' @param wavelength Numeric vector of wavelengths \[nm\]
#' @param par a vector with name {chl, a_dg_440, a_nap_440, bb_p_550}. Optionally provide a_g, bb_p slope with names {a_dg_s, a_nap_s, bb_p_gamma}.
#' Total IOP will be computed with the combination of oac you provide.
#' If you don't provide any, the returned value will be equal to the pure water IOP.
#'
#' @export

iop_from_oac <- function(wavelength, par) {
  stopifnot(is.numeric(wavelength))
  stopifnot(is.numeric(par))
  stopifnot(!is.null(names(par)))

  .Call("c_iop_from_oac", as.numeric(wavelength), par)
}


# iop_from_oac <- function(
#     wavelength,
#     oac,
#     rrs = NULL,
#     optically_shallow = F
#     ) {
#
#   if (any(oac < 0)) return(NA)
#
#   # Phyto concentration to absorption  --------------------------------------
#   if ("chl" %in% names(oac)) {
#
#     a0 <- memoised_a0_a1_phyto(wavelength)$a0
#     a1 <- memoised_a0_a1_phyto(wavelength)$a1
#
#     # Plankton absorption (1/m)
#     aph_440 <- 0.06 * (oac["chl"])^0.65 # [mg/m^3] #Prieur & Satyendranath (1981)
#
#     a_phy <- (a0 + a1 * log(aph_440)) * aph_440
#
#     # if (any(a_phy < 0)) {
#     #   # rlang::warn("Some a_phy inferiro to 0")
#     #   a_phy[a_phy < 0] <- 0
#     # }
#   } else {
#     a_phy <- 0
#   }
#
#   # CDOM absorption ---------------------------------------------------------
#   if ("a_g_440" %in% names(oac)) {
#     ## CDOM+NAP absorption coefficient [1/m]
#
#     Ga_CDOM <- 1
#     Oa_CDOM <- 0
#
#     abs_CDM_440 <- (Ga_CDOM * oac["a_g_440"]) + Oa_CDOM # [1/m], CDOM abs. coeff. at 440 [nm]
#
#     # 3 possibilities for CDOM slope: manual, parametric, default
#
#     if (all(c("a_g_s_g", "a_g_s_d") %in% names(oac))) {
#       cdom_slope <- oac["a_g_s_g"] + oac["a_g_s_d"]
#     } else if (optically_shallow && tibble::is_tibble(rrs)) {
#       cdom_slope <- 0.015 + (
#         0.002 / (0.6 + (rrs$rrs[which.min(abs(rrs$wavelength - 443))] /
#                           rrs$rrs[which.min(abs(rrs$wavelength - 555))])))
#     } else {
#       cdom_slope <- 0.017 # Model Default
#     }
#
#     a_g <- abs_CDM_440 * exp(-cdom_slope * (wavelength - 440))
#
#   } else {
#     a_g <- 0
#   }
#
#   # Non algal particle absorption (a_nap)  ----------------------------------
#   if ("a_nap_440" %in% names(oac)) {
#     if ("a_nap_s_d" %in% names(oac)) {
#       s_nap = as.numeric(oac["a_nap_s_d"])
#     } else {
#       s_nap <- 0.01160 # Model Default
#     }
#
#     # TODO: why those coeff ?
#     Ga_nap <- 1# [m^2/mg]
#     Oa_nap <- 0
#
#     a_nap_440 <-  (Ga_nap * oac["a_nap_440"]) + Oa_nap# [1/m], SPM abs. coeff. at 440 [nm]
#
#     a_nap <- a_nap_440 * exp(-s_nap * (wavelength - 440))
#
#     # abs_CDM = abs_CDOM + abs_X
#   } else {
#     a_nap <- 0
#   }
#
#   # Particulate backscattering ----------------------------------------------
#   if ("bb_p_550" %in% names(oac)) {
#
#     # 3 possibilities for bb slope: manual, parametric, default
#     if ("bb_p_gamma" %in% names(oac)) {
#       bb_gamma <- oac["bb_p_gamma"]
#     } else if (optically_shallow && tibble::is_tibble(rrs)) {
#       bb_gamma <- 2 * (1 - (1.2 * exp(
#         -0.9 * (rrs$rrs[which.min(abs(rrs$wavelength - 443))] /
#                   rrs$rrs[which.min(abs(rrs$wavelength - 555))])
#       )))
#     } else {
#       bb_gamma <- 0.46 # Model Default
#     }
#
#     bb_p <- oac["bb_p_550"] * ((wavelength / 550) ^ -bb_gamma)
#
#   } else {
#     bb_p <- 0
#   }
#
#   # Pure water IOP ----------------------------------------------------------
#   iop_w <- pure_water_iop(wavelength)
#
#   out <- list("a" = iop_w$a_w + a_phy + a_g, "bb" = iop_w$bb_w + bb_p)
#
#   # out <- tibble(
#   #   wavelength,
#   #   "a_phy" = if ("chl" %in% names(oac)) a_phy else NULL,
#   #   "a_nap" = if ("a_nap_440" %in% names(oac)) a_nap else NULL,
#   #   "a_g" = if ("a_g_440" %in% names(oac)) a_g else NULL,
#   #   "bb_p" = if ("bb_p_550" %in% names(oac)) bb_p else NULL,
#   #   "a" = iop_w$a_w + a_phy + a_nap + a_g,
#   #   "bb" = iop_w$bb_w + bb_p
#   # )
#
#   # if (verbose) {
#   #   # Assume `out` is a tibble with at least one column 'wavelength' and others are spectra
#   #   cols_to_plot <- setdiff(names(out), "wavelength")
#   #
#   #   ply <- purrr::reduce(
#   #     cols_to_plot,
#   #     .init = plot_ly(out, x = ~wavelength),
#   #     .f = function(p, colname) {
#   #       add_lines(p, y = out[[colname]], name = colname)
#   #     }
#   #   )
#   #
#   #   # Display the plot
#   #   print(ply)
#   # }
#
#   return(out)
# }
