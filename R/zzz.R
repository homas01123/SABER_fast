.onLoad <- function(libname, pkgname) {
  # C data ------------------------------------------------------------------
  data("a_w", package = pkgname, envir = environment())
  .Call("c_load_pure_water", a_w$wavelength, a_w$a_w)

  data("a0_a1_phyto", package = pkgname, envir = environment())
  .Call("c_load_a0_a1", a0_a1_phyto$wavelength, a0_a1_phyto$a0, a0_a1_phyto$a1)

  data("r_rs_b_gamache", package = pkgname, envir = environment())
  r_rs_b <- r_rs_b_gamache %>%
    dplyr::select(class, wavelength, r_rs_b_mean) %>%
    tidyr::pivot_wider(
      names_from = "class",
      values_from = "r_rs_b_mean",
      names_prefix = "r_rs_b_"
    )
  .Call("c_load_r_rs_b", r_rs_b[[1]], as.matrix(r_rs_b[,-1]))

  # Registry ----------------------------------------------------------------
  register_input_preparer("input_am03", function(par, rrs) {
    input_am03(par, rrs)
  })

  register_forward_model("am03", function(inputs) {
    forward_am03(
      wavelength = inputs$wavelength,
      iop = inputs$iop,
      water_type = inputs$water_type,
      theta_view = inputs$theta_view,
      theta_sun = inputs$theta_sun,
      h_w = inputs$h_w,
      r_b = inputs$r_b
    )
  })

  register_forward_model("lee98", function(inputs) {
    forward_lee98(
      wavelength = inputs$wavelength,
      iop = inputs$iop,
      theta_view = inputs$theta_view,
      theta_sun = inputs$theta_sun,
      optically_shallow = inputs$optically_shallow,
      h_w = inputs$h_w,
      r_b = inputs$r_b
    )
  })

  register_objective_function("log-ll", function(modelled, observed, par) {
    log_ll(modelled = modelled, observed = observed, sd = par[["sd"]])
  })

  register_objective_function("rss", function(modelled, observed, par) {
    rss(modelled = modelled, observed = observed)
  })

  register_objective_function("lee99", function(modelled, observed, par) {
    lee99(modelled = modelled, observed = observed, wavelength = par[["wavelength"]])
  })
}

.onUnload <- function(libpath)
{
  # Tell the DLL to free every malloc() it owns
  invisible(.Call("c_saber_reset_tables"))
  library.dynam.unload("SABER", libpath)
}
