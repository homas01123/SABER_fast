.onLoad <- function(libname, pkgname) {
  # C data ------------------------------------------------------------------
  data("a_w", package = pkgname, envir = environment())
  .Call("c_load_pure_water", a_w$wavelength, a_w$a_w)

  data("a0_a1_phyto", package = pkgname, envir = environment())
  .Call("c_load_a0_a1", a0_a1_phyto$wavelength, a0_a1_phyto$a0, a0_a1_phyto$a1)

  # data("r_rs_b_gamache", package = pkgname, envir = environment()) # to use Baie-Gamache Rb
  # r_rs_b <- r_rs_b_gamache %>%
  #   dplyr::select(class, wavelength, r_rs_b_mean) %>%
  #   tidyr::pivot_wider(
  #     names_from = "class",
  #     values_from = "r_rs_b_mean",
  #     names_prefix = "r_rs_b_"
  #   )
  # .Call("c_load_r_rs_b", r_rs_b[[1]], as.matrix(r_rs_b[,-1]))

  tryCatch({
    data("r_rs_b_egsl", package = pkgname, envir = environment()) # to use combined EGSL Rb
    r_rs_b <- r_rs_b_egsl %>%
      dplyr::select(class, wavelength, r_rs_b_mean) %>%
      tidyr::pivot_wider(
        names_from = "class",
        values_from = "r_rs_b_mean",
        names_prefix = "r_rs_b_"
      )

    # Ensure wavelength is numeric
    wavelength_vec <- as.numeric(r_rs_b[[1]])
    r_rs_b_matrix <- as.matrix(r_rs_b[,-1])

    # CRITICAL: Ensure column names are set (C code requires them)
    colnames(r_rs_b_matrix) <- names(r_rs_b)[-1]

    # Validate before calling C
    if (!is.numeric(wavelength_vec)) {
      stop("Wavelength must be numeric")
    }
    if (!is.matrix(r_rs_b_matrix) || !is.numeric(r_rs_b_matrix)) {
      stop("Bottom reflectance data must be a numeric matrix")
    }
    if (is.null(colnames(r_rs_b_matrix))) {
      stop("Matrix column names are NULL")
    }

    .Call("c_load_r_rs_b", wavelength_vec, r_rs_b_matrix)
    options(SABER.available_classes = names(r_rs_b)[-1])
    #message("Successfully loaded ", length(names(r_rs_b)[-1]), " benthic classes")
  }, error = function(e) {
    warning("Failed to load benthic reflectance data: ", e$message, immediate. = TRUE)
  })

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
  
  register_input_preparer("input_am03_sicf", function(par, rrs, par_meta = NULL) {
    input_am03_sicf(par, rrs, par_meta)
  })
  
  register_forward_model("am03_sicf", function(inputs) {
    forward_am03_sicf(
      wavelength = inputs$wavelength,
      iop = inputs$iop,
      water_type = inputs$water_type,
      theta_view = inputs$theta_view,
      theta_sun = inputs$theta_sun,
      h_w = inputs$h_w,
      r_b = inputs$r_b,
      chl = inputs$chl,
      a_dg_443 = inputs$a_dg_443,
      phi_f = inputs$phi_f,
      include_sicf = inputs$include_sicf,
      sicf_model = inputs$sicf_model,
      depth_integration = inputs$depth_integration,
      lat = inputs$lat,
      lon = inputs$lon,
      date_time = inputs$date_time,
      return_components = inputs$return_components
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
