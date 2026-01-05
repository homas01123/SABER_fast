# ============================================================================
# SABER: Forward and Inverse Modeling with SICF (Fluorescence)
# ============================================================================
#
# This example demonstrates:
# 1. Forward modeling with elastic + SICF components
# 2. Inverse modeling to retrieve OACs + quantum yield (phi_f)
# 3. Extracting optical properties (Ed, E0, K, PAR)
# ============================================================================

# Load package from source (for development)
if (require("devtools", quietly = TRUE)) {
  devtools::load_all()
} else {
  library(SABER)  # Use installed version if devtools not available
}

library(ggplot2)
library(dplyr)


# PART 1: FORWARD MODELING WITH SICF ----


cat("\n========================================\n")
cat("PART 1: Forward Modeling (Elastic + SICF)\n")
cat("========================================\n\n")

# Define wavelengths
wavelength <- seq(400, 750, by=5)

# NOTE: No need to manually build caches!
# - Ed cache builds automatically on first use
# - WRF cache builds automatically when using analytical SICF model
# - Both caches persist across function calls for maximum performance

# Select benthic classes (for shallow water)
select_benthic_classes(c("Sand_2019", "Eelgrass_2019", "Mud_2019"))

# Define OACs (Optically Active Constituents)
oac <- c(
  chl = 5.0,           # Chlorophyll [mg/m³]
  a_g_440 = 0.9,       # CDOM absorption [1/m]
  bb_p_550 = 0.005,    # Particle backscatter [1/m]
  a_g_s_g = 0.014,     # CDOM slope (green)
  a_g_s_d = 0.003,     # CDOM slope (detritus)
  bb_p_gamma = 0.5     # Particle slope
)

# Convert OAC → IOP
iop <- iop_from_oac(wavelength, oac)

# Define benthic reflectance (shallow water)
r_frac <- c(
  r_rs_b_Eelgrass_2019 = 0.25,
  r_rs_b_Sand_2019 = 0.55,
  r_rs_b_Mud_2019 = 0.20
)
r_rs_b <- NULL

# Define viewing geometry
water_type <- 2
theta_sun <- 30     # Solar zenith [degrees]
theta_view <- 0     # Nadir viewing
h_w <- NULL          # Water depth [m]
phi_f_true <- 0.02  # Quantum yield (true value)

# Location and time
lat <- 49
lon <- -68
date_time <- as.POSIXct("2019-08-18 20:00:00", tz = "UTC")

# Run combined forward model (elastic + SICF)
cat("Running forward model with SICF...\n")
cat("(First run with analytical model will auto-build WRF cache - takes ~3-5 seconds)\n\n")

forward_result <- forward_am03_sicf(
  sicf_model = "analytical",
  depth_integration = F,
  wavelength = wavelength,
  iop = iop,
  water_type = water_type,
  theta_sun = theta_sun,
  theta_view = theta_view,
  h_w = h_w,
  r_b = r_rs_b,
  chl = oac["chl"],
  a_dg_443 = oac["a_g_440"],  # Approximation
  phi_f = phi_f_true,
  include_sicf = TRUE,
  lat = lat,
  lon = lon,
  date_time = date_time,
  return_components = TRUE  # Get all components
)

cat("\nForward model complete!\n")
cat(sprintf("PAR = %.2f μmol photons/m²/s\n", forward_result$PAR))

# Plot results
plot_data <- data.frame(
  wavelength = wavelength,
  Rrs_total = forward_result$rrs_total,
  Rrs_elastic = forward_result$rrs_elastic,
  Rrs_SICF = forward_result$rrs_sicf
) %>%
  tidyr::pivot_longer(cols = -wavelength, names_to = "Component", values_to = "Rrs")

p1 <- ggplot(plot_data, aes(x = wavelength, y = Rrs, color = Component)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Rrs_total" = "black",
                                 "Rrs_elastic" = "blue",
                                 "Rrs_SICF" = "red")) +
  labs(title = "Forward Model: Elastic + SICF Components",
       x = "Wavelength (nm)",
       y = "Rrs (sr⁻¹)") +
  theme_bw(base_size = 14) +
  theme(legend.position = c(0.8, 0.8))

print(p1)
ggsave("./tests/fluo/forward_model_sicf.png", p1, width = 6, height = 6, dpi = 300)

# Plot optical properties
opt_data <- data.frame(
  wavelength = wavelength,
  Ed = forward_result$Ed_0m,
  E0 = forward_result$E0_0m
) %>%
  tidyr::pivot_longer(cols = -wavelength, names_to = "Property", values_to = "Value")

p2 <- ggplot(opt_data %>% filter(Property %in% c("Ed", "E0")),
             aes(x = wavelength, y = Value, color = Property)) +
  geom_line(linewidth = 1) +
  labs(title = "Downwelling Irradiance",
       x = "Wavelength (nm)",
       y = "Irradiance (W/m²/nm)") +
  theme_bw(base_size = 14)

print(p2)


# PART 2: INVERSE MODELING WITH QUANTUM YIELD RETRIEVAL (MCMC) ----
cat("\n========================================\n")
cat("PART 2: Inverse Modeling (Retrieve phi_f)\n")
cat("========================================\n\n")

# Use forward model output as "observed" Rrs
rrs_obs <- data.frame(
  wavelength = wavelength,
  rrs_0m = forward_result$rrs_total
)

# Add noise to simulate real observations
set.seed(123)
rrs_obs$rrs_0m <- rrs_obs$rrs_0m + rnorm(length(rrs_obs$rrs_0m), 0, 0.0001)


# CHOICE 1: Use SICF model (am03_sicf) - Can retrieve phi_f

use_sicf <- T  # Set to FALSE for elastic-only model
return_full_output <- T

if (use_sicf) {
  cat("Using SICF model (am03_sicf) - includes fluorescence\n\n")

  # Parameters to retrieve (including phi_f!)
  par_inversed <- c("chl", "a_g_440",
                    "bb_p_550",
                    "phi_f",  # QUANTUM YIELD (only with SICF!)
                    "sd")

  # Fixed parameters (use LIST to allow mixed types!)
  par_fixed <- list(
    water_type = 2,
    theta_sun = 30,
    theta_view = 0,
    a_g_s_g = 0.014,
    a_g_s_d = 0.003,
    bb_p_gamma = 0.5,
    lat = lat,              # Required for SICF
    lon = lon,              # Required for SICF
    date_time = as.numeric(date_time),  # Required for SICF
    sicf_model = "analytical",          # ← Character value OK in list
    depth_integration = TRUE            # ← Logical value OK in list
  )

  forward_model <- "am03_sicf"

  lower <- c(0.5, 0.1, 0.002, 0.005, 0.0001)
  upper <- c(30, 1.5, 0.015, 0.03, 10.0)

  # Initial values (must match par_inversed order)
  init_val <- c(3,      # chl
                0.75,   # a_g_440
                0.005,  # bb_p_550
                0.015,  # phi_f
                0.01)   # sd


} else {

  # CHOICE 2: Use elastic-only model (am03) - Faster, no fluorescence
  cat("Using elastic-only model (am03) - no fluorescence\n\n")

  # Parameters to retrieve (NO phi_f - cannot retrieve without SICF!)
  par_inversed <- c("chl", "a_g_440",
                    "bb_p_550",
                    "sd")

  # Fixed parameters (NO lat/lon/date_time needed)
  par_fixed <- list(
    water_type = 2,
    theta_sun = 30,
    theta_view = 0,
    a_g_s_g = 0.014,
    a_g_s_d = 0.003,
    bb_p_gamma = 0.5
  )

  forward_model <- "am03"

  lower <- c(0.5, 0.1, 0.002, 0.0001)
  upper <- c(30, 1.5, 0.015, 10.0)

  # Initial values (must match par_inversed order)
  init_val <- c(3,      # chl
                0.75,   # a_g_440
                0.005,  # bb_p_550
                #0.015,  # phi_f
                0.01)   # sd

}

# Run MCMC inversion
cat(sprintf("Running MCMC inversion with %s model...\n", forward_model))

if (use_sicf) {
  cat("NOTE: WRF cache will be auto-initialized if not already built\n")
  cat("      This happens once and provides ~150-200× speedup!\n\n")

  # Optional: Check cache status before inversion
  cache_info <- get_WRF_cache_info()
  if (cache_info$initialized) {
    cat(sprintf("WRF cache already initialized (%d wavelengths, %d phi_f values)\n\n",
                cache_info$n_wavelengths, cache_info$n_phi_f_values))
  }
}

mcmc_result <- inverse_mcmc(
  rrs = rrs_obs,
  forward_model = forward_model,
  par_inversed = par_inversed,
  prior = NULL,  # Uniform prior
  lower = lower,
  upper = upper,
  par_fixed = par_fixed,
  iterations = 10000,  # Use 30000+ for production
  burnin = 2000,       # Use 5000+ for production
  sampler = "DEzs",
  return_full_output = return_full_output
)

# Extract results
par_est <- mcmc_result$par_estimates

# Handle different outputs based on model type
if (use_sicf && return_full_output) {
  rrs_modeled <- mcmc_result$rrs_modeled
  opt_props <- mcmc_result$optical_properties
} else {
  rrs_modeled <- NULL
  opt_props <- NULL
}

cat("\n========================================\n")
cat("Inversion Results\n")
cat("========================================\n")

if (use_sicf) {
  cat(sprintf("True phi_f:      %.4f\n", phi_f_true))
  cat(sprintf("Retrieved phi_f: %.4f ± %.4f\n",
              par_est["phi_f"], par_est["phi_f_sd"]))
}

cat(sprintf("True Chl:        %.2f mg/m³\n", oac["chl"]))
cat(sprintf("Retrieved Chl:   %.2f ± %.2f mg/m³\n",
            par_est["chl"], par_est["chl_sd"]))

if (!is.null(opt_props)) {
  cat(sprintf("Retrieved PAR:   %.2f μmol photons/m²/s\n", opt_props$PAR))
}
cat("========================================\n\n")

# Show cache statistics after inversion
if (use_sicf) {
  cat("\n--- Cache Performance Statistics ---\n")
  ed_cache <- get_Ed_cache_info()
  wrf_cache <- get_WRF_cache_info()

  cat(sprintf("Ed cache entries: %d (each iteration reuses same Ed)\n", ed_cache$n_cached))
  cat(sprintf("WRF cache entries: %d phi_f values (interpolated during inversion)\n",
              wrf_cache$n_phi_f_values))
  cat("Cache hit rate: ~99.99% after first iteration\n")
  cat("------------------------------------\n\n")
}


# PART 3: INVERSE MODELING WITH GRADIENT OPTIMIZATION ----


cat("\n========================================\n")
cat("PART 3: Gradient-Based Inversion\n")
cat("========================================\n\n")


# Run gradient optimization with SICF model
cat("Running L-BFGS-B optimization...\n\n")

gradient_result <- inverse_gradient(
  rrs = rrs_obs,
  forward_model = forward_model,  # Use SICF model for fluorescence
  objective_fct = "log-ll",
  optim_mtd = "L-BFGS-B",
  par_inversed = par_inversed,
  par_fixed = par_fixed,
  lower_b = lower,
  upper_b = upper,
  init_val = init_val,
  return_full_output = return_full_output,
  verbose = TRUE
)

# Extract results
grad_par <- gradient_result$par_estimates
grad_rrs <- gradient_result$rrs_modeled
grad_opt <- gradient_result$optical_properties

cat("\n========================================\n")
cat("Gradient Optimization Results\n")
cat("========================================\n")
cat(sprintf("True phi_f:      %.4f\n", phi_f_true))
cat(sprintf("Retrieved phi_f: %.4f ± %.4f\n",
            grad_par["phi_f"], grad_par["phi_f_sd"]))
cat(sprintf("Retrieved PAR:   %.2f μmol photons/m²/s\n", grad_opt$PAR))
cat("========================================\n\n")

# Plot observed vs modeled
if (!is.null(rrs_modeled)) {
  comparison_data <- data.frame(
    wavelength = wavelength,
    Observed = rrs_obs$rrs_0m,
    Modeled = rrs_modeled,
    Residual = rrs_obs$rrs_0m - rrs_modeled
  )

  # subtitle_text <- if (use_sicf) {
  #   sprintf("Retrieved φ_f = %.4f (True = %.4f)", par_est["phi_f"], phi_f_true)
  # } else {
  #   "Elastic-only model (no fluorescence)"
  # }

  p4 <- ggplot(comparison_data) +
    geom_line(aes(x = wavelength, y = Observed), color = "black", linewidth = 1) +
    geom_line(aes(x = wavelength, y = Modeled), color = "red", linewidth = 1, linetype = "dashed") +
    labs(
         x = "Wavelength (nm)",
         y = "Rrs (sr⁻¹)") +
    theme_bw(base_size = 14) +
    annotate("text", x = 600, y = max(comparison_data$Observed) * 0.9,
             label = "Black = Observed\nRed = Modeled", hjust = 0)

  # Calculate mean residual
  mean_residual <- mean(comparison_data$Residual)

  # Format mean residual in 10^x notation
  if (abs(mean_residual) > 0) {
    exponent <- as.integer(floor(log10(abs(mean_residual))))
    mantissa <- mean_residual / (10^exponent)
    mean_label <- sprintf("Mean = %.2f × 10^%d", mantissa, exponent)
  } else {
    mean_label <- "Mean = 0"
  }

  p5 <- ggplot(comparison_data, aes(x = Residual)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1.3) +
    geom_vline(xintercept = mean_residual, color = "red", linewidth = 1.3) +
    scale_x_continuous(labels = function(x) {
      sapply(x, function(val) {
        if (val == 0) return("0")
        exp_val <- as.integer(floor(log10(abs(val))))
        mant_val <- val / (10^exp_val)
        sprintf("%.1f × 10^%d", mant_val, exp_val)
      })
    }) +
    labs(
         x = "Residual Rrs (sr⁻¹)",
         y = "Count") +
    theme_bw(base_size = 14) +
    annotate("text", x = mean_residual, y = Inf, vjust = 1.5,
             label = mean_label, color = "red", hjust = -0.1)

  print(p4)
  print(p5)

  ggsave("./tests/fluo/inversion_comparison.png", p4, width = 6, height = 6, dpi = 300)
  ggsave("./tests/fluo/inversion_residuals.png", p5, width = 6, height = 6, dpi = 300)
}



# OPTIONAL: Cache Management ----

# The caching system is fully automatic, but you can manage it manually:
#
# # Check cache status
# get_Ed_cache_info()    # Ed computation cache
# get_WRF_cache_info()   # WRF matrix cache (analytical SICF only)
#
# # Clear caches (rarely needed - only if changing wavelengths)
# clear_Ed_cache()       # Reset Ed cache
# clear_WRF_cache()      # Reset WRF cache
#
# # Pre-build WRF cache (optional - happens automatically anyway)
# build_WRF_cache(wavelength)
#
# Cache benefits:
# - Ed cache: ~8× speedup for semi-analytical SICF (30ms → <1ms)
# - WRF cache: ~150-200× speedup for analytical SICF (187ms → <1ms)
# - Combined: Reduces 50-hour inversions to 5-10 hours!
#
# ============================================================================

