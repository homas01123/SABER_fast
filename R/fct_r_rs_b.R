#' Load bottom reflectance data
#' @param wavelength numeric vector
#' @param rrs_b matrix (columns = class names)
#' @export
load_r_rs_b <- function(wavelength, rrs_b) {
  stopifnot(is.numeric(wavelength), is.matrix(rrs_b))
  .Call("c_load_r_rs_b", as.numeric(wavelength), rrs_b)
}

#' Compute interpolated Rrs_b mixture
#' @param fractions named numeric vector (e.g., c(sand = 0.6, seagrass = 0.4))
#' @export
compute_r_rs_b_lmm <- function(fractions) {
  .Call("c_compute_r_rs_b_lmm", fractions)
}

#' Select and load specific benthic classes
#'
#' @param classes Character vector of class names to use (e.g., c("Mud_2019", "Eelgrass_2019", "Sand_2019"))
#' @export
select_benthic_classes <- function(classes) {
  # Get full dataset
  data("r_rs_b_egsl", package = "SABER", envir = environment())

  # Available classes
  available <- unique(r_rs_b_egsl$class)

  # Validate user selection
  if (!all(classes %in% available)) {
    missing <- classes[!classes %in% available]
    stop(paste("Classes not found:", paste(missing, collapse = ", "),
               "\nAvailable classes:", paste(available, collapse = ", ")))
  }

  # Filter selected classes
  r_rs_b_selected <- r_rs_b_egsl %>%
    dplyr::filter(class %in% classes) %>%
    dplyr::select(class, wavelength, r_rs_b_mean) %>%
    tidyr::pivot_wider(
      names_from = "class",
      values_from = "r_rs_b_mean",
      names_prefix = "r_rs_b_"
    )

  # Reload into C cache with proper conversion
  wavelength_vec <- as.numeric(r_rs_b_selected[[1]])
  # Convert to matrix while preserving column names
  r_rs_b_matrix <- as.matrix(r_rs_b_selected[,-1])
  # Ensure column names are set
  colnames(r_rs_b_matrix) <- names(r_rs_b_selected)[-1]
  
  # Debug: verify structure before calling C
  stopifnot(
    is.numeric(wavelength_vec),
    is.matrix(r_rs_b_matrix),
    !is.null(colnames(r_rs_b_matrix))
  )
  
  load_r_rs_b(wavelength_vec, r_rs_b_matrix)
  # Reset the interpolated cache; it will be rebuilt at the correct wavelength
  # grid on the next call to iop_from_oac() or forward_am03().
  .Call("c_reset_cache")

  # Store selected classes in options
  options(SABER.selected_classes = paste0("r_rs_b_", classes))

  message("Loaded benthic classes: ", paste(classes, collapse = ", "))
  invisible(classes)
}

#' List available benthic classes
#' @export
list_benthic_classes <- function() {
  data("r_rs_b_egsl", package = "SABER", envir = environment())
  unique(r_rs_b_egsl$class)
}

# Resolve benthic reflectance vector from a retrieved parameter vector.
# Handles two cases:
#   1. mix_sand  – 2-class linear mixing; selected_classes must have length 2.
#   2. r_rs_b_*  – N-class fractions already in par_vec.
# Returns a spectral r_b vector (same wavelength grid as the class data), or NULL.
.resolve_benthic_fractions <- function(par_vec, selected_classes) {
  nms     <- names(par_vec)
  if (is.null(nms)) return(NULL)

  # Case 1: 2-class mixing parameter
  idx_mix <- which(nms == "mix_sand")
  if (length(idx_mix) > 0) {
    mix <- as.numeric(par_vec[[idx_mix[1]]])
    if (is.na(mix)) return(NULL)
    fracs <- setNames(c(mix, 1 - mix),
                      c(selected_classes[1], selected_classes[2]))
    return(compute_r_rs_b_lmm(fracs))
  }

  # Case 2: explicit per-class fractions
  keep  <- grepl("^r_rs_b_", nms) & !grepl("_sd$", nms)
  if (!any(keep)) return(NULL)
  fracs <- setNames(as.numeric(unlist(par_vec[keep])), nms[keep])
  if (all(is.na(fracs))) return(NULL)
  compute_r_rs_b_lmm(fracs)
}
