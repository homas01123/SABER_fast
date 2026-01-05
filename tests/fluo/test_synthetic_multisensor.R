# ============================================================================
# SABER: Multi-Sensor Synthetic Data Inversion Test
# ============================================================================
#
# This script tests SICF inversions on synthetic data with three wavelength
# configurations:
# 1. High-resolution (0.5 nm) - reference
# 2. Gaussian spectral response (3 nm FWHM)
# 3. Sentinel-2 + FLEX-FLORIS combined bands
#
# Tests parameter retrieval accuracy across different spectral configurations
# ============================================================================

# Load package
if (require("devtools", quietly = TRUE)) {
  devtools::load_all()
} else {
  library(SABER)
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)
library(doParallel)
library(progressr)
library(foreach)

cat("\n============================================================\n")
cat("SABER Multi-Sensor SICF Inversion Test\n")
cat("============================================================\n\n")


# PART 1: DEFINE PARAMETER SPACE ----

cat("PART 1: Defining parameter space for synthetic data\n")
cat("------------------------------------------------------------\n")

# Parameter ranges based on natural variability
chl_range <- c(0.05, 1000)        # mg/m³ (oligotrophic to hypereutrophic)
a_g_440_range <- c(0.02, 4)       # m⁻¹ (clear to highly turbid)
bb_p_550_range <- c(0.0005, 0.03) # m⁻¹ (low to high particle load)
phi_f_range <- c(0.005, 0.05)    # quantum yield range

# Pre-build WRF cache with EXTENDED phi_f range
# cat("Pre-building WRF cache for extended phi_f range [0.0002, 0.05]...\n")
# cat("(Default cache range [0.005, 0.03] is insufficient for synthetic data)\n")
# build_WRF_cache(wavelength = seq(400, 800, by = 0.5),
#                 phi_f_grid = seq(0.0002, 0.05, length.out = 50),
#                 verbose = TRUE)
# cat("WRF cache built successfully\n\n")

# Fixed viewing geometry
water_type <- 2
theta_view <- 0   # Nadir
theta_sun <- 30   # degrees
lat <- 49
lon <- -68
date_time <- as.POSIXct("2019-08-18 14:00:00", tz = "UTC")

# Generate expanded grid (stratified sampling for computational efficiency)
# Use log-scale for chl and a_g_440 due to their large dynamic range
set.seed(42)

n_samples <- 50  # Number of samples per parameter (adjust for speed vs coverage)

synthetic_params <- expand.grid(
  chl = exp(seq(log(chl_range[1]), log(chl_range[2]), length.out = n_samples)),
  a_g_440 = exp(seq(log(a_g_440_range[1]), log(a_g_440_range[2]), length.out = 15)),
  bb_p_550 = seq(bb_p_550_range[1], bb_p_550_range[2], length.out = 10),
  phi_f = seq(phi_f_range[1], phi_f_range[2], length.out = 8)
)

# Subsample for computational feasibility (optional)
if (nrow(synthetic_params) > 5000) {
  sample_idx <- sample(1:nrow(synthetic_params), 5000)
  synthetic_params <- synthetic_params[sample_idx, ]
}

cat(sprintf("Generated %d synthetic scenarios\n", nrow(synthetic_params)))
cat(sprintf("  Chl range: %.3f - %.1f mg/m³\n", min(synthetic_params$chl), max(synthetic_params$chl)))
cat(sprintf("  a_g_440 range: %.3f - %.2f m⁻¹\n", min(synthetic_params$a_g_440), max(synthetic_params$a_g_440)))
cat(sprintf("  bb_p_550 range: %.4f - %.3f m⁻¹\n", min(synthetic_params$bb_p_550), max(synthetic_params$bb_p_550)))
cat(sprintf("  phi_f range: %.4f - %.3f\n\n", min(synthetic_params$phi_f), max(synthetic_params$phi_f)))


# PART 2: WAVELENGTH CONFIGURATIONS ----

cat("PART 2: Defining wavelength configurations\n")
cat("------------------------------------------------------------\n")

# Configuration 1: High-resolution (0.5 nm) - REFERENCE
wavelength_hires <- seq(400, 800, by = 0.5)
cat(sprintf("Config 1 (High-res): %d bands, 0.5 nm resolution\n", length(wavelength_hires)))

# Configuration 2: Gaussian spectral response (3 nm FWHM)
# Center wavelengths every 3 nm
wavelength_gaussian_centers <- seq(400, 800, by = 3)

# Function to apply Gaussian spectral response
apply_gaussian_srf <- function(spectrum_hires, wv_hires, wv_centers, fwhm = 3) {
  sigma <- fwhm / (2 * sqrt(2 * log(2)))  # Convert FWHM to sigma

  spectrum_convolved <- sapply(wv_centers, function(center) {
    weights <- exp(-((wv_hires - center)^2) / (2 * sigma^2))
    weights <- weights / sum(weights)  # Normalize
    sum(spectrum_hires * weights)
  })

  return(spectrum_convolved)
}

cat(sprintf("Config 2 (Gaussian SRF): %d bands, 3 nm FWHM\n", length(wavelength_gaussian_centers)))

# Configuration 3: Sentinel-2 + FLEX-FLORIS
# Sentinel-2 MSI bands (blue and green for elastic component)
s2_bands <- data.frame(
  band = c("B2", "B3", "B4"),
  center = c(492.4, 559.8, 664.6),
  fwhm = c(66, 36, 31)
)

# FLEX-FLORIS bands (Fluorescence Imaging Spectrometer)
# Spectral range: 500-780 nm
# High-resolution oxygen bands: 0.1 nm sampling at 686-697 nm and 759-769 nm
# Other bands: 0.5-2.0 nm in red edge, chlorophyll absorption, and PRI regions

# Create FLEX-FLORIS band configuration
# O2-A band (759-769 nm): 0.1 nm sampling
o2a_bands <- data.frame(
  band = paste0("O2A_", 1:101),
  center = seq(759, 769, by = 0.1),
  fwhm = 0.1,
  region = "O2-A"
)

# O2-B band (686-697 nm): 0.1 nm sampling
o2b_bands <- data.frame(
  band = paste0("O2B_", 1:111),
  center = seq(686, 697, by = 0.1),
  fwhm = 0.1,
  region = "O2-B"
)

# Red edge and chlorophyll absorption (500-686 nm, 697-759 nm): 1 nm sampling
# Exclude O2-B region
red_edge_1 <- data.frame(
  band = paste0("RE1_", 1:187),
  center = seq(500, 686, by = 1),
  fwhm = 1,
  region = "Red Edge 1"
)

red_edge_2 <- data.frame(
  band = paste0("RE2_", 1:63),
  center = seq(697, 759, by = 1),
  fwhm = 1,
  region = "Red Edge 2"
)

# NIR region (769-780 nm): 1 nm sampling
nir_bands <- data.frame(
  band = paste0("NIR_", 1:12),
  center = seq(769, 780, by = 1),
  fwhm = 1,
  region = "NIR"
)

# Combine all FLEX-FLORIS bands
flex_floris_bands <- rbind(
  red_edge_1,
  o2b_bands,
  red_edge_2,
  o2a_bands,
  nir_bands
) %>% arrange(center)

# Combine S2 + FLEX bands
combined_bands <- rbind(
  s2_bands %>% dplyr::mutate(sensor = "Sentinel-2", region = "VIS"),
  flex_floris_bands %>% dplyr::mutate(sensor = "FLEX-FLORIS")
) %>% arrange(center)

cat(sprintf("Config 3 (S2+FLEX): %d bands (3 S2 + %d FLEX)\n",
            nrow(combined_bands), nrow(flex_floris_bands)))
cat("  Sentinel-2: B2 (492 nm), B3 (560 nm), B4 (665 nm)\n")
cat(sprintf("  FLEX-FLORIS: %d bands (500-780 nm)\n", nrow(flex_floris_bands)))
cat("    - O2-B: 686-697 nm (0.1 nm sampling)\n")
cat("    - O2-A: 759-769 nm (0.1 nm sampling)\n")
cat("    - Red Edge/Chl: 500-686, 697-759, 769-780 nm (1 nm sampling)\n\n")


# PART 2B: VISUALIZE RSRFs ----
cat("\n============================================================\n")
cat("PART 2B: Visualizing Relative Spectral Response Functions\n")
cat("============================================================\n\n")

# Create RSRF visualization
# Define wavelength grid for plotting (400-800 nm)
wv_plot <- seq(400, 800, by = 0.5)

# Function to calculate normalized Gaussian RSRF
calc_gaussian_rsrf <- function(wv_plot, center, fwhm) {
  sigma <- fwhm / (2 * sqrt(2 * log(2)))
  rsrf <- exp(-((wv_plot - center)^2) / (2 * sigma^2))
  return(rsrf)
}

# Configuration 1: High-res (0.5 nm) - treat as continuous spectrum (single group)
config1_rsrf <- data.frame(
  wavelength = wavelength_hires,
  rsrf = 1,
  config = "High-res (0.5 nm)",
  band_id = 1  # Single group for continuous line
)

# Configuration 2: Gaussian (3 nm FWHM) - create RSRFs for each center wavelength
config2_rsrf_list <- lapply(seq_along(wavelength_gaussian_centers), function(i) {
  center <- wavelength_gaussian_centers[i]
  rsrf <- calc_gaussian_rsrf(wv_plot, center, fwhm = 3)
  data.frame(
    wavelength = wv_plot,
    rsrf = rsrf,
    config = "Gaussian (3 nm FWHM)",
    band_id = i,
    center = center
  )
})
config2_rsrf <- bind_rows(config2_rsrf_list)

# Configuration 3: S2 + FLEX - create RSRFs for each band
config3_rsrf_list <- lapply(1:nrow(combined_bands), function(i) {
  center <- combined_bands$center[i]
  fwhm <- combined_bands$fwhm[i]
  
  # For narrow bands (FWHM < 2nm), create narrow Gaussian
  # For broad bands (S2), use actual FWHM
  rsrf <- calc_gaussian_rsrf(wv_plot, center, fwhm)
  
  data.frame(
    wavelength = wv_plot,
    rsrf = rsrf,
    config = "S2 + FLEX",
    band_id = i,
    center = center,
    sensor = combined_bands$sensor[i]
  )
})
config3_rsrf <- bind_rows(config3_rsrf_list)

# Combine all RSRFs
all_rsrf <- bind_rows(
  config1_rsrf %>% dplyr::select(wavelength, rsrf, config, band_id),
  config2_rsrf %>% dplyr::select(wavelength, rsrf, config, band_id),
  config3_rsrf %>% dplyr::select(wavelength, rsrf, config, band_id)
) %>%
  mutate(config_factor = factor(config, 
                                levels = c("High-res (0.5 nm)", 
                                          "Gaussian (3 nm FWHM)", 
                                          "S2 + FLEX")))

# Function to map wavelength to spectral color
wavelength_to_color <- function(wv) {
  # Spectral color ramp:
  # 400-500 nm: blue to cyan
  # 500-600 nm: green to yellow  
  # 600-700 nm: orange to red
  # 700-800 nm: dark red to black
  colors <- colorRampPalette(c(
    "#0000FF",  # 400 nm - deep blue
    "#00FFFF",  # 500 nm - cyan
    "#00FF00",  # 550 nm - green
    "#FFFF00",  # 580 nm - yellow
    "#FF8000",  # 620 nm - orange
    "#FF0000",  # 680 nm - red
    "#8B0000",  # 750 nm - dark red
    "#000000"   # 800 nm - black
  ))(401)  # 400-800 nm = 401 values
  
  idx <- pmax(1, pmin(401, round(wv - 399)))
  return(colors[idx])
}

# Create RSRF plot - simple black style
p_rsrf <- ggplot() +
  # High-res: filled ribbon to show continuous narrow-band coverage
  geom_ribbon(data = all_rsrf %>% filter(config == "High-res (0.5 nm)"),
              aes(x = wavelength, ymin = 0.98, ymax = 1.02),
              fill = "grey20", alpha = 0.8) +
  # Gaussian: individual Gaussian curves
  geom_line(data = all_rsrf %>% filter(config == "Gaussian (3 nm FWHM)"),
            aes(x = wavelength, y = rsrf, group = band_id), 
            linewidth = 0.5, alpha = 0.6, color = "black") +
  # S2+FLEX: individual band responses
  geom_line(data = all_rsrf %>% filter(config == "S2 + FLEX"),
            aes(x = wavelength, y = rsrf, group = band_id), 
            linewidth = 0.6, alpha = 0.7, color = "grey30") +
  facet_wrap(~ config_factor, ncol = 1) +
  scale_x_continuous(name = expression(paste("Wavelength (", lambda, ") [nm]")),
                     limits = c(400, 800),
                     breaks = seq(400, 800, 100)) +
  scale_y_continuous(name = "Normalized RSRF",
                     limits = c(0, 1.05),
                     breaks = seq(0, 1, 0.25)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 18, color = 'black'),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.ticks.length = unit(.25, "cm"),
    strip.text = element_text(size = 16, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 1.2),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)
  )

print(p_rsrf)
ggsave("./tests/fluo/rsrf_comparison.png", p_rsrf, 
       width = 8, height = 10, dpi = 300)

cat("\nRSRF comparison plot saved: rsrf_comparison.png\n\n")

# Create parameter summary table for presentation
cat("============================================================\n")
cat("PARAMETER RANGES FOR SYNTHETIC DATASET\n")
cat("============================================================\n\n")

param_summary <- data.frame(
  Parameter = c("Chl-a", "a_g(440)", "b_bp(550)", "φ_f"),
  Range = c(
    sprintf("%.2f - %.0f mg/m³", min(synthetic_params$chl), max(synthetic_params$chl)),
    sprintf("%.3f - %.2f m⁻¹", min(synthetic_params$a_g_440), max(synthetic_params$a_g_440)),
    sprintf("%.4f - %.3f m⁻¹", min(synthetic_params$bb_p_550), max(synthetic_params$bb_p_550)),
    sprintf("%.3f - %.3f", min(synthetic_params$phi_f), max(synthetic_params$phi_f))
  ),
  `N Samples` = c(
    length(unique(synthetic_params$chl)),
    length(unique(synthetic_params$a_g_440)),
    length(unique(synthetic_params$bb_p_550)),
    length(unique(synthetic_params$phi_f))
  ),
  Sampling = c("Log-scale", "Log-scale", "Linear", "Linear")
)

print(param_summary)

cat("\n------------------------------------------------------------\n")
cat(sprintf("Total Scenarios: %d (before subsampling)\n", nrow(synthetic_params)))
cat("Fixed Parameters:\n")
cat(sprintf("  Water Type: %d\n", water_type))
cat(sprintf("  θ_sun: %d°\n", theta_sun))
cat(sprintf("  θ_view: %d° (Nadir)\n", theta_view))
cat(sprintf("  a_g slope (λ > 440): %.3f nm⁻¹\n", 0.014))
cat(sprintf("  a_g slope (λ < 440): %.3f nm⁻¹\n", 0.003))
cat(sprintf("  b_bp power law: %.1f\n", 0.5))
cat("============================================================\n\n")

# Save parameter summary as CSV for easy copying to slides
write.csv(param_summary, 
          file = "./tests/fluo/parameter_summary.csv", 
          row.names = FALSE)

cat("Parameter summary saved: parameter_summary.csv\n\n")


# Function to apply band-specific spectral response
# For narrow bands (FWHM < 2nm), use direct interpolation to preserve fine structure
# For broad bands (Sentinel-2), use Gaussian SRF convolution
apply_band_srf <- function(spectrum_hires, wv_hires, bands_df) {
  spectrum_bands <- sapply(1:nrow(bands_df), function(i) {
    center <- bands_df$center[i]
    fwhm <- bands_df$fwhm[i]

    # If FWHM < 2nm, use direct interpolation (preserves fine structure)
    if (fwhm < 2.0) {
      return(approx(wv_hires, spectrum_hires, xout = center)$y)
    }

    # Otherwise use Gaussian SRF convolution
    sigma <- fwhm / (2 * sqrt(2 * log(2)))
    weights <- exp(-((wv_hires - center)^2) / (2 * sigma^2))
    weights <- weights / sum(weights)

    return(sum(spectrum_hires * weights))
  })

  return(spectrum_bands)
}


# PART 3: GENERATE SYNTHETIC OBSERVATIONS ----

cat("PART 3: Generating synthetic Rrs observations\n")
cat("------------------------------------------------------------\n")

# Add realistic noise function
add_noise <- function(rrs, snr = 100) {
  noise <- rnorm(length(rrs), mean = 0, sd = mean(rrs) / snr)
  rrs + noise
}

cat(sprintf("Running forward model for %d scenarios...\n", nrow(synthetic_params)))

# Generate synthetic data for all scenarios
synthetic_data_list <- lapply(1:nrow(synthetic_params), function(i) {
  params <- synthetic_params[i, ]

  if (i %% 100 == 0) cat(sprintf("  Processing scenario %d/%d...\n", i, nrow(synthetic_params)))

  # Build OAC vector
  oac <- c(
    chl = params$chl,
    a_g_440 = params$a_g_440,
    bb_p_550 = params$bb_p_550,
    a_g_s_g = 0.014,
    a_g_s_d = 0.003,
    bb_p_gamma = 0.5
  )

  # Get IOPs for high-resolution
  iop_hires <- iop_from_oac(wavelength_hires, oac)

  # Run forward model (high-resolution)
  tryCatch({
    forward_hires <- forward_am03_sicf(
      sicf_model = "semi_analytical",
      depth_integration = F,
      wavelength = wavelength_hires,
      iop = iop_hires,
      water_type = water_type,
      theta_sun = theta_sun,
      theta_view = theta_view,
      h_w = NULL,
      r_b = NULL,
      chl = params$chl,
      a_dg_443 = params$a_g_440,
      phi_f = params$phi_f,
      include_sicf = TRUE,
      lat = lat,
      lon = lon,
      date_time = date_time,
      return_components = FALSE, verbose = T
    )

    # forward_hires is already a vector (return_components = FALSE)
    # Apply spectral responses

    # Config 1: High-res (subsample to match Gaussian grid for comparison)
    rrs_hires_sub <- approx(wavelength_hires, forward_hires,
                             xout = wavelength_hires)$y

    # Config 2: Gaussian SRF
    rrs_gaussian <- apply_gaussian_srf(forward_hires, wavelength_hires,
                                       wavelength_gaussian_centers, fwhm = 3)

    # Config 3: S2 + FLEX bands
    rrs_s2flex <- apply_band_srf(forward_hires, wavelength_hires, combined_bands)

    synthetic_data_list = list(
      scenario_id = i,
      true_chl = params$chl,
      true_a_g_440 = params$a_g_440,
      true_bb_p_550 = params$bb_p_550,
      true_phi_f = params$phi_f,
      config1_hires = data.frame(
        wavelength = wavelength_hires,
        rrs_0m = add_noise(rrs_hires_sub)
      ),
      config2_gaussian = data.frame(
        wavelength = wavelength_gaussian_centers,
        rrs_0m = add_noise(rrs_gaussian)
      ),
      config3_s2flex = data.frame(
        wavelength = combined_bands$center,
        rrs_0m = add_noise(rrs_s2flex)
      )
    )

  }, error = function(e) {
    NULL
  })
})

# Remove failed scenarios
synthetic_data_list <- synthetic_data_list[!sapply(synthetic_data_list, is.null)]

cat(sprintf("\nSuccessfully generated %d synthetic datasets\n\n", length(synthetic_data_list)))

# Create nested tibbles for each configuration (like hypernet structure)
cat("Creating nested data structures for parallel inversion...\n")

# Config 1: High-res nested data
config1_nested <- tibble(
  scenario_id = sapply(synthetic_data_list, function(x) x$scenario_id),
  true_chl = sapply(synthetic_data_list, function(x) x$true_chl),
  true_a_g_440 = sapply(synthetic_data_list, function(x) x$true_a_g_440),
  true_bb_p_550 = sapply(synthetic_data_list, function(x) x$true_bb_p_550),
  true_phi_f = sapply(synthetic_data_list, function(x) x$true_phi_f),
  data = lapply(synthetic_data_list, function(x) x$config1_hires)
)

# Config 2: Gaussian SRF nested data
config2_nested <- tibble(
  scenario_id = sapply(synthetic_data_list, function(x) x$scenario_id),
  true_chl = sapply(synthetic_data_list, function(x) x$true_chl),
  true_a_g_440 = sapply(synthetic_data_list, function(x) x$true_a_g_440),
  true_bb_p_550 = sapply(synthetic_data_list, function(x) x$true_bb_p_550),
  true_phi_f = sapply(synthetic_data_list, function(x) x$true_phi_f),
  data = lapply(synthetic_data_list, function(x) x$config2_gaussian)
)

# Config 3: S2+FLEX nested data
config3_nested <- tibble(
  scenario_id = sapply(synthetic_data_list, function(x) x$scenario_id),
  true_chl = sapply(synthetic_data_list, function(x) x$true_chl),
  true_a_g_440 = sapply(synthetic_data_list, function(x) x$true_a_g_440),
  true_bb_p_550 = sapply(synthetic_data_list, function(x) x$true_bb_p_550),
  true_phi_f = sapply(synthetic_data_list, function(x) x$true_phi_f),
  data = lapply(synthetic_data_list, function(x) x$config3_s2flex)
)

cat(sprintf("Config 1 (High-res): %d nested observations\n", nrow(config1_nested)))
cat(sprintf("Config 2 (Gaussian): %d nested observations\n", nrow(config2_nested)))
cat(sprintf("Config 3 (S2+FLEX): %d nested observations\n\n", nrow(config3_nested)))


# PART 3B: VISUALIZE SYNTHETIC SPECTRA ----
cat("\n============================================================\n")
cat("PART 3B: Visualizing synthetic spectral data\n")
cat("============================================================\n\n")

# Create spectral plots for each configuration
# Sample subset for visualization (avoid overplotting)
set.seed(123)
n_plot_samples <- min(500, nrow(config1_nested))
plot_idx <- sample(1:nrow(config1_nested), n_plot_samples)

# Prepare data for Config 1 (High-res)
config1_plot_data <- config1_nested[plot_idx, ] %>%
  mutate(spectrum_data = purrr::map2(data, true_chl, ~ {
    .x %>% mutate(chl = .y)
  })) %>%
  select(scenario_id, true_chl, spectrum_data) %>%
  tidyr::unnest(spectrum_data)

# Calculate plot dimensions
xmin <- 400; xmax <- 800
ymin <- 0; ymax <- max(config1_plot_data$rrs_0m, na.rm = TRUE) * 1.1
asp_rat <- (xmax - xmin) / (ymax - ymin)

# Config 1: High-resolution spectral plot
p1_spectra <- ggplot(config1_plot_data,
                     aes(x = wavelength, y = rrs_0m, group = scenario_id, color = chl)) +
  geom_line(alpha = 0.6, linewidth = 0.8) +
  scale_color_viridis_c(trans = "log10",
                        name = "Chl-a\n(mg/m³)",
                        breaks = c(0.1, 1, 10, 100, 1000)) +
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), ylim = c(ymin, ymax),
              expand = FALSE, clip = "on") +
  scale_x_continuous(name = expression(paste("Wavelength (", lambda, ") [nm]")),
                     limits = c(xmin, xmax),
                     breaks = seq(xmin, xmax, 100)) +
  scale_y_continuous(name = expression(paste(italic("R")["rs"], " [sr"^{-1}, "]")),
                     limits = c(ymin, ymax),
                     labels = function(x) format(x, scientific = FALSE, digits = 4)) +
  #ggtitle(sprintf("Config 1: High-Resolution Synthetic Rrs (0.5 nm), n = %d", n_plot_samples)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20, color = 'black', angle = 0),
        axis.text.y = element_text(size = 20, color = 'black', angle = 0),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.position = c(0.95, 0.90),
        legend.justification = c("right", "top"),
        legend.title = element_text(colour = "black", size = 20, face = "bold"),
        legend.text = element_text(colour = "black", size = 18, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5,
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black",
                                        size = 0.5, linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey80",
                                        linewidth = 0.2, linetype = "solid"),
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        plot.margin = unit(c(0.5, 1.0, 0.5, 0.5), "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

print(p1_spectra)
ggsave("./tests/fluo/config1_hires_spectra.png", p1_spectra, width = 6, height = 6,
       scale = 1.5, dpi = 300)

# Prepare data for Config 2 (Gaussian)
config2_plot_data <- config2_nested[plot_idx, ] %>%
  mutate(spectrum_data = purrr::map2(data, true_chl, ~ {
    .x %>% mutate(chl = .y)
  })) %>%
  select(scenario_id, true_chl, spectrum_data) %>%
  tidyr::unnest(spectrum_data)

# Calculate plot dimensions
ymin2 <- 0; ymax2 <- max(config2_plot_data$rrs_0m, na.rm = TRUE) * 1.1
asp_rat2 <- (xmax - xmin) / (ymax2 - ymin2)

# Config 2: Gaussian SRF spectral plot
p2_spectra <- ggplot(config2_plot_data,
                     aes(x = wavelength, y = rrs_0m, group = scenario_id, color = chl)) +
  geom_line(alpha = 0.6, linewidth = 0.8) +
  scale_color_viridis_c(trans = "log10",
                        name = "Chl-a\n(mg/m³)",
                        breaks = c(0.1, 1, 10, 100, 1000)) +
  coord_fixed(ratio = asp_rat2, xlim = c(xmin, xmax), ylim = c(ymin2, ymax2),
              expand = FALSE, clip = "on") +
  scale_x_continuous(name = expression(paste("Wavelength (", lambda, ") [nm]")),
                     limits = c(xmin, xmax),
                     breaks = seq(xmin, xmax, 100)) +
  scale_y_continuous(name = expression(paste(italic("R")["rs"], " [sr"^{-1}, "]")),
                     limits = c(ymin2, ymax2),
                     labels = function(x) format(x, scientific = FALSE, digits = 4)) +
  #ggtitle(sprintf("Config 2: Gaussian SRF Synthetic Rrs (3 nm FWHM), n = %d", n_plot_samples)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20, color = 'black', angle = 0),
        axis.text.y = element_text(size = 20, color = 'black', angle = 0),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.position = c(0.95, 0.90),
        legend.justification = c("right", "top"),
        legend.title = element_text(colour = "black", size = 20, face = "bold"),
        legend.text = element_text(colour = "black", size = 18, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5,
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black",
                                        size = 0.5, linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey80",
                                        linewidth = 0.2, linetype = "solid"),
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        plot.margin = unit(c(0.5, 1.0, 0.5, 0.5), "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

print(p2_spectra)
ggsave("./tests/fluo/config2_gaussian_spectra.png", p2_spectra,
       width = 6, height = 6, scale = 1.5, dpi = 300)

# Prepare data for Config 3 (S2+FLEX)
config3_plot_data <- config3_nested[plot_idx, ] %>%
  mutate(spectrum_data = purrr::map2(data, true_chl, ~ {
    .x %>% mutate(chl = .y)
  })) %>%
  select(scenario_id, true_chl, spectrum_data) %>%
  tidyr::unnest(spectrum_data)

# Calculate plot dimensions for S2+FLEX (starts at 492 nm)
xmin3 <- 490; xmax3 <- 800
ymin3 <- 0; ymax3 <- max(config3_plot_data$rrs_0m, na.rm = TRUE) * 1.1
asp_rat3 <- (xmax3 - xmin3) / (ymax3 - ymin3)

# Config 3: S2+FLEX spectral plot
p3_spectra <- ggplot(config3_plot_data,
                     aes(x = wavelength, y = rrs_0m, group = scenario_id, color = chl)) +
  geom_line(alpha = 0.6, linewidth = 0.8) +
  scale_color_viridis_c(trans = "log10",
                        name = "Chl-a\n(mg/m³)",
                        breaks = c(0.1, 1, 10, 100, 1000)) +
  coord_fixed(ratio = asp_rat3, xlim = c(xmin3, xmax3), ylim = c(ymin3, ymax3),
              expand = FALSE, clip = "on") +
  scale_x_continuous(name = expression(paste("Wavelength (", lambda, ") [nm]")),
                     limits = c(xmin3, xmax3),
                     breaks = seq(500, 800, 100)) +
  scale_y_continuous(name = expression(paste(italic("R")["rs"], " [sr"^{-1}, "]")),
                     limits = c(ymin3, ymax3),
                     labels = function(x) format(x, scientific = FALSE, digits = 4)) +
  #ggtitle(sprintf("Config 3: Sentinel-2 + FLEX-FLORIS (3 S2 + 474 FLEX), n = %d", n_plot_samples)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20, color = 'black', angle = 0),
        axis.text.y = element_text(size = 20, color = 'black', angle = 0),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.position = c(0.95, 0.90),
        legend.justification = c("right", "top"),
        legend.title = element_text(colour = "black", size = 20, face = "bold"),
        legend.text = element_text(colour = "black", size = 18, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5,
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black",
                                        size = 0.5, linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey80",
                                        linewidth = 0.2, linetype = "solid"),
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        plot.margin = unit(c(0.5, 1.0, 0.5, 0.5), "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

print(p3_spectra)
ggsave("./tests/fluo/config3_s2flex_spectra.png", p3_spectra, width = 6, height = 6,
       scale = 1.5, dpi = 300)

# Prepare IOP data for the same scenarios
cat("\nGenerating IOP plots for synthetic scenarios...\n")

iop_plot_data <- config1_nested[plot_idx, ] %>%
  mutate(iop_data = purrr::map(scenario_id, ~ {
    # Get corresponding parameters
    idx <- which(config1_nested$scenario_id == .x)[1]
    chl_val <- config1_nested$true_chl[idx]
    ag_val <- config1_nested$true_a_g_440[idx]
    bbp_val <- config1_nested$true_bb_p_550[idx]

    # Reconstruct OAC parameters
    oac <- c(
      chl = chl_val,
      a_g_440 = ag_val,
      bb_p_550 = bbp_val,
      a_g_s_g = 0.014,
      a_g_s_d = 0.003,
      bb_p_gamma = 0.5
    )

    # Calculate IOPs
    iop <- iop_from_oac(wavelength_hires, oac)

    data.frame(
      wavelength = wavelength_hires,
      a_total = iop$a,
      bb_total = iop$bb,
      chl = chl_val
    )
  })) %>%
  select(scenario_id, true_chl, iop_data) %>%
  tidyr::unnest(iop_data)

# Reshape data for combined plot with single axis
iop_plot_long <- iop_plot_data %>%
  tidyr::pivot_longer(cols = c(a_total, bb_total),
                      names_to = "iop_type",
                      values_to = "value") %>%
  mutate(iop_label = ifelse(iop_type == "a_total", "a(λ)", "bb(λ)"))

# Create combined plot with log scale
p_iop_combined <- ggplot(iop_plot_long,
                         aes(x = wavelength, y = value,
                             group = interaction(scenario_id, iop_type),
                             color = chl, linetype = iop_label)) +
  geom_line(alpha = 0.6, linewidth = 0.8) +
  scale_color_viridis_c(trans = "log10",
                        name = "Chl-a\n(mg/m³)",
                        breaks = c(0.1, 1, 10, 100, 1000)) +
  scale_linetype_manual(name = "IOP",
                        values = c("a(λ)" = "solid", "bb(λ)" = "dashed")) +
  scale_x_continuous(name = expression(paste("Wavelength (", lambda, ") [nm]")),
                     limits = c(xmin, xmax),
                     breaks = seq(xmin, xmax, 100)) +
  scale_y_log10(name = expression(paste("IOP [m"^{-1}, "]")),
                labels = function(x) format(x, scientific = FALSE, digits = 3)) +
  #ggtitle(sprintf("Inherent Optical Properties (IOPs), n = %d", n_plot_samples)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20, color = 'black', angle = 0),
        axis.text.y = element_text(size = 20, color = 'black', angle = 0),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.position = c(0.95, 0.90),
        legend.justification = c("right", "top"),
        legend.title = element_text(colour = "black", size = 20, face = "bold"),
        legend.text = element_text(colour = "black", size = 18, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5,
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.box = "vertical",
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black",
                                        size = 0.5, linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey80",
                                        linewidth = 0.2, linetype = "solid"),
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        plot.margin = unit(c(0.5, 1.0, 0.5, 0.5), "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

print(p_iop_combined)
ggsave("./tests/fluo/iop_combined_spectra.png", p_iop_combined,
       width = 6, height = 6, scale = 1.5, dpi = 300)

cat("\nSpectral plots saved:\n")
cat("  - config1_hires_spectra.png\n")
cat("  - config2_gaussian_spectra.png\n")
cat("  - config3_s2flex_spectra.png\n")
cat("  - iop_combined_spectra.png\n\n")


# PART 4: RUN INVERSIONS IN PARALLEL ----

cat("PART 4: Running inversions on synthetic data (PARALLEL)\n")
cat("------------------------------------------------------------\n")

# Inversion setup
par_inversed <- c("chl", "a_g_440", "bb_p_550", "phi_f", "sd")
par_fixed <- list(
  water_type = 2,
  theta_sun = 30,
  theta_view = 0,
  a_g_s_g = 0.014,
  a_g_s_d = 0.003,
  bb_p_gamma = 0.5,
  lat = lat,
  lon = lon,
  date_time = as.numeric(date_time),
  sicf_model = "semi_analytical",
  depth_integration = FALSE
)


lower <- c(0.03, 0.02, 0.0003, 0.005, 0.0001)
upper <- c(1000, 5, 0.04, 0.06, 10.0)
init_val <- c(5, 0.5, 0.005, 0.01, 0.001)

# Function to run parallel inversion on nested data
run_parallel_inversion <- function(nested_data, config_name, inv_mode) {

  cat(sprintf("\n--- %s + %s ---\n", config_name, inv_mode))

  # Partition into chunks for parallel processing
  chunked_nested <- split(nested_data,
                          cut(seq_len(nrow(nested_data)),
                              breaks = num_cores, labels = FALSE))

  if (inv_mode == "mcmc") {

    exec_time <- system.time({
      with_progress({
        p <- progressor(steps = length(chunked_nested))

        saber_results <- foreach(chunk = chunked_nested,
                                 .packages = c("purrr", "SABER", "BayesianTools", "dplyr", "progressr"),
                                 .export = c("inverse_mcmc", "par_inversed", "lower", "upper",
                                             "par_fixed", "forward_model")) %dopar% {

                                               result <- chunk %>%
                                                 mutate(
                                                   inversion_estim = purrr::map(
                                                     data,
                                                     ~ {
                                                       tryCatch({
                                                         inverse_mcmc(
                                                           rrs = .x,
                                                           forward_model = forward_model,
                                                           par_inversed = par_inversed,
                                                           prior = NULL,
                                                           lower = lower,
                                                           best = NULL,
                                                           upper = upper,
                                                           par_fixed = par_fixed, # Pass the modified local copy
                                                           iterations = 25000,
                                                           burnin = 5000,
                                                           sampler = "DEzs"
                                                         )
                                                       }, error = function(e) {
                                                         NULL
                                                       })
                                                     }
                                                   )
                                                 )

                                               p()  # increment progress bar
                                               result
                                             }
      })
    })

  }

  if (inv_mode == "gradient") {

    exec_time <- system.time({
      with_progress({
        p <- progressor(steps = length(chunked_nested))

        saber_results <- foreach(chunk = chunked_nested,
                                 .packages = c("purrr", "SABER", "dplyr", "progressr", "numDeriv", "MASS"),
                                 .export = c("inverse_gradient", "par_inversed", "lower", "init_val",
                                             "upper", "par_fixed", "objective_factory",
                                             "parse_inverse_parameter")) %dopar% {

                                               result <- chunk %>%
                                                 mutate(
                                                   inversion_estim = purrr::map(
                                                     data,
                                                     ~ {
                                                       tryCatch({
                                                         inverse_gradient(
                                                           rrs = .x,
                                                           forward_model = "am03_sicf",
                                                           objective_fct = "log-ll",
                                                           optim_mtd = "L-BFGS-B",
                                                           par_inversed = par_inversed,
                                                           par_fixed = par_fixed,
                                                           lower_b = lower,
                                                           init_val = init_val,
                                                           upper_b = upper,
                                                           verbose = FALSE
                                                         )
                                                       }, error = function(e) {
                                                         NULL
                                                       })
                                                     }
                                                   )
                                                 )

                                               p()  # increment progress bar
                                               result
                                             }
      })
    })

  }


  # Combine results
  combined_results <- bind_rows(saber_results)

  cat(sprintf("Completed in %.1f seconds\n", exec_time[3]))

  return(combined_results)
}

# Setup parallel processing
num_cores <- parallel::detectCores() - 2
cl <- makeCluster(num_cores)
doParallel::registerDoParallel(cl)

cat(sprintf("Using %d cores for parallel processing\n\n", num_cores))

# Run inversions for each configuration
cat("Starting parallel inversions...\n")

results_config1 <- run_parallel_inversion(config1_nested, "Config 1: High-res (0.5 nm)",
                                          inv_mode = "mcmc")
results_config2 <- run_parallel_inversion(config2_nested, "Config 2: Gaussian (3 nm FWHM)",
                                          inv_mode = "mcmc")
results_config3 <- run_parallel_inversion(config3_nested, "Config 3: S2 + FLEX",
                                          inv_mode = "mcmc")

# Stop cluster
stopCluster(cl)

cat("\n All inversions complete!\n")


# PART 5: ANALYZE RESULTS ----
cat("\n\n============================================================\n")
cat("PART 5: Analyzing retrieval accuracy\n")
cat("============================================================\n\n")

# Extract and organize results from parallel inversions
cat("Extracting inversion results...\n")

# Function to extract results from nested inversion output
extract_nested_results <- function(results_nested, config_name) {

  results_nested %>%
    mutate(
      inversion_tidy = purrr::map(inversion_estim, ~ {
        if (!is.null(.x)) {
          as_tibble(t(.x))
        } else {
          tibble(chl = NA, a_g_440 = NA, bb_p_550 = NA, phi_f = NA, sd = NA,
                 chl_sd = NA, a_g_440_sd = NA, bb_p_550_sd = NA, phi_f_sd = NA, sd_sd = NA)
        }
      })
    ) %>%
    unnest(cols = c(inversion_tidy)) %>%
    select(-inversion_estim, -data) %>%
    mutate(config = config_name)
}

# Extract for all configs
config1_results <- extract_nested_results(results_config1, "config1")
config2_results <- extract_nested_results(results_config2, "config2")
config3_results <- extract_nested_results(results_config3, "config3")

# Combine all results
results_df <- bind_rows(config1_results, config2_results, config3_results) %>%
  filter(!is.na(chl)) %>%  # Remove failed inversions
  mutate(
    # Calculate relative errors
    err_chl = (chl - true_chl) / true_chl * 100,
    err_a_g_440 = (a_g_440 - true_a_g_440) / true_a_g_440 * 100,
    err_bb_p_550 = (bb_p_550 - true_bb_p_550) / true_bb_p_550 * 100,
    err_phi_f = (phi_f - true_phi_f) / true_phi_f * 100,
    # Friendly config names
    config_name = case_when(
      config == "config1" ~ "High-res (0.5 nm)",
      config == "config2" ~ "Gaussian (3 nm FWHM)",
      config == "config3" ~ "S2 + FLEX"
    )
  )

cat(sprintf("Total successful inversions: %d\n", nrow(results_df)))
cat(sprintf("  Config 1 (High-res): %d\n", sum(results_df$config == "config1")))
cat(sprintf("  Config 2 (Gaussian): %d\n", sum(results_df$config == "config2")))
cat(sprintf("  Config 3 (S2+FLEX): %d\n\n", sum(results_df$config == "config3")))

# Summary statistics
summary_stats <- results_df %>%
  group_by(config_name) %>%
  summarise(
    n = n(),
    chl_bias = mean(err_chl, na.rm = TRUE),
    chl_rmse = sqrt(mean(err_chl^2, na.rm = TRUE)),
    a_g_bias = mean(err_a_g_440, na.rm = TRUE),
    a_g_rmse = sqrt(mean(err_a_g_440^2, na.rm = TRUE)),
    bb_p_bias = mean(err_bb_p_550, na.rm = TRUE),
    bb_p_rmse = sqrt(mean(err_bb_p_550^2, na.rm = TRUE)),
    phi_f_bias = mean(err_phi_f, na.rm = TRUE),
    phi_f_rmse = sqrt(mean(err_phi_f^2, na.rm = TRUE))
  )

print(summary_stats)

cat("\n\nRetrieval Accuracy Summary:\n")
cat("------------------------------------------------------------\n")
for (i in 1:nrow(summary_stats)) {
  cat(sprintf("\n%s (%d successful retrievals):\n",
              summary_stats$config_name[i], summary_stats$n[i]))
  cat(sprintf("  Chl:      Bias = %+.1f%%, RMSE = %.1f%%\n",
              summary_stats$chl_bias[i], summary_stats$chl_rmse[i]))
  cat(sprintf("  a_g_440:  Bias = %+.1f%%, RMSE = %.1f%%\n",
              summary_stats$a_g_bias[i], summary_stats$a_g_rmse[i]))
  cat(sprintf("  bb_p_550: Bias = %+.1f%%, RMSE = %.1f%%\n",
              summary_stats$bb_p_bias[i], summary_stats$bb_p_rmse[i]))
  cat(sprintf("  phi_f:    Bias = %+.1f%%, RMSE = %.1f%%\n",
              summary_stats$phi_f_bias[i], summary_stats$phi_f_rmse[i]))
}


# PART 6: VISUALIZATION ----

cat("\n\n============================================================\n")
cat("PART 6: Creating comparison plots\n")
cat("============================================================\n\n")

library(ggExtra)

# Function to create beautiful validation scatter plots
plot_validation_scatter <- function(results_df, param, param_label, param_range, config_colors) {

  # Prepare data
  uncertainty_col <- paste0(param, "_sd")
  plot_data <- results_df %>%
    rename(actual = !!paste0("true_", param),
           predicted = !!param,
           uncertainty = !!uncertainty_col) %>%
    mutate(config_factor = factor(config_name,
                                   levels = c("High-res (0.5 nm)",
                                             "Gaussian (3 nm FWHM)",
                                             "S2 + FLEX")))

  # Calculate statistics for each config
  stats_df <- plot_data %>%
    group_by(config_name) %>%
    summarise(
      n = n(),
      r2 = cor(actual, predicted, use = "complete.obs")^2,
      rmse = sqrt(mean((predicted - actual)^2, na.rm = TRUE)),
      bias = mean(predicted - actual, na.rm = TRUE),
      .groups = "drop"
    )

  xmin <- param_range[1]; xmax <- param_range[2]
  ymin <- xmin; ymax <- xmax
  asp_rat <- (xmax - xmin) / (ymax - ymin)

  # Calculate log breaks
  min_val <- log10(xmin)
  max_val <- log10(xmax)
  num_breaks <- 4
  by_value <- (max_val - min_val) / num_breaks
  logrange <- seq(-10, 10, by = signif(by_value, digits = 0))
  breaks <- 10^(logrange)
  minor_breaks <- rep(1:9, length(breaks)) * (10^rep(logrange, each = 9))

  # Create smoothed uncertainty ribbon data
  # Sort data by actual values for ribbon
  ribbon_data <- plot_data %>%
    arrange(actual) %>%
    filter(!is.na(predicted), !is.na(uncertainty), predicted > 0) %>%
    # Use rolling mean to smooth the uncertainty band
    mutate(
      pred_smooth = predict(loess(predicted ~ actual, data = ., span = 0.3), newdata = .),
      unc_smooth = predict(loess(uncertainty ~ actual, data = ., span = 0.3), newdata = .),
      # Make uncertainty relative (as a multiplier) for log scale visibility
      # Use confidence interval based on std dev
      ymin = pmax(pred_smooth * exp(-unc_smooth/pred_smooth), xmin * 0.1),
      ymax = pmin(pred_smooth * exp(unc_smooth/pred_smooth), xmax * 10)
    )

  # Debug output
  cat(sprintf("\nParameter: %s\n", param))
  cat(sprintf("Uncertainty column: %s\n", uncertainty_col))
  cat(sprintf("Ribbon data rows: %d\n", nrow(ribbon_data)))
  if(nrow(ribbon_data) > 0) {
    cat(sprintf("Uncertainty range: [%.6f, %.6f]\n",
                min(ribbon_data$unc_smooth, na.rm = TRUE),
                max(ribbon_data$unc_smooth, na.rm = TRUE)))
    cat(sprintf("Predicted range: [%.6f, %.6f]\n",
                min(ribbon_data$pred_smooth, na.rm = TRUE),
                max(ribbon_data$pred_smooth, na.rm = TRUE)))
    cat(sprintf("Ribbon ymin range: [%.6f, %.6f]\n",
                min(ribbon_data$ymin, na.rm = TRUE),
                max(ribbon_data$ymin, na.rm = TRUE)))
    cat(sprintf("Ribbon ymax range: [%.6f, %.6f]\n",
                min(ribbon_data$ymax, na.rm = TRUE),
                max(ribbon_data$ymax, na.rm = TRUE)))
  }

  # Create base plot
  p <- ggplot(plot_data, aes(x = actual, y = predicted)) +
    geom_density_2d(na.rm = TRUE, bins = 6,

                    linewidth = 0.25, colour = "black", show.legend = FALSE) +
    geom_point(aes(fill = config_factor, shape = config_factor),
               alpha = 0.25, size = 3, show.legend = TRUE) +
    scale_fill_manual(name = "",
                      values = config_colors) +
    scale_shape_manual(name = "",
                       values = c(21, 22, 24)) +
    #geom_rug(size = 1.1, show.legend = FALSE, alpha = 0.25) +
    geom_abline(slope = 1, linetype = "dashed", intercept = 0,
                colour = "black", na.rm = FALSE, size = 1.3, show.legend = FALSE) +
    geom_smooth(size = 1, level = 0.95, show.legend = FALSE,
                linetype = "solid", color = "goldenrod2",
                se = TRUE, method = "lm") +
    coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
                ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
    scale_x_log10(breaks = breaks, minor_breaks = minor_breaks,
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = breaks, minor_breaks = minor_breaks,
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(size = 1) +
    labs(x = param_label,
         y = param_label) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 20, color = 'black', angle = 0),
          axis.text.y = element_text(size = 20, color = 'black', angle = 0),
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          axis.ticks.length = unit(.25, "cm"),
          legend.box.just = "right",
          legend.spacing = unit(-0.5, "cm"),
          legend.position = c(0.05, 0.95),
          legend.title = element_blank(),
          legend.text = element_text(colour = "black", size = 15, face = "plain"),
          legend.background = element_rect(fill = NA, size = 0.5,
                                          linetype = "solid", colour = 0),
          legend.key = element_blank(),
          legend.justification = c("left", "top"),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "black",
                                         size = 0.5, linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey80",
                                         linewidth = 0.2, linetype = "solid"),
          plot.margin = unit(c(0.5, 1.0, 0.5, 0.5), "cm"),
          legend.direction = "vertical", legend.box = "vertical",
          legend.text.align = 0,
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

  # Add marginal density plots
  # p_with_margins <- ggMarginal(p, type = "densigram", bins = 30, alpha = 0.7,
  #                               fill = "grey60")

  return(list(plot = p, stats = stats_df))
}

# Define color scheme for configurations
config_colors <- c("High-res (0.5 nm)" = "#377EB8",
                   "Gaussian (3 nm FWHM)" = "#4DAF4A",
                   "S2 + FLEX" = "brown4")

# Create validation plots for each parameter
cat("Creating validation scatter plots with marginal distributions...\n")

# Chl
chl_plot <- plot_validation_scatter(results_df, "chl",
                                     expression(paste(italic("Chl-a"), " [mg ", m^{-3}, "]")),
                                     c(0.01, 2000), config_colors)
print(chl_plot$plot)
ggsave("./tests/fluo/validation_chl.png", chl_plot$plot,
       width = 4.5, height = 4.5, scale = 1.25, dpi = 300)
print(chl_plot$stats)

# a_g_440
ag_plot <- plot_validation_scatter(results_df, "a_g_440",
                                    expression(paste(italic("a")[g], "(440) [m"^{-1}, "]")),
                                    c(0.01, 10), config_colors)
print(ag_plot$plot)
ggsave("./tests/fluo/validation_ag440.png", ag_plot$plot,
       width = 4.5, height = 4.5, scale = 1.25, dpi = 300)
print(ag_plot$stats)

# bb_p_550
bbp_plot <- plot_validation_scatter(results_df, "bb_p_550",
                                     expression(paste(italic("b")[bp], "(550) [m"^{-1}, "]")),
                                     c(0.0003, 0.05), config_colors)
print(bbp_plot$plot)
ggsave("./tests/fluo/validation_bbp550.png", bbp_plot$plot,
       width = 4.5, height = 4.5, scale = 1.25, dpi = 300)
print(bbp_plot$stats)

# phi_f
phi_plot <- plot_validation_scatter(results_df, "phi_f",
                                     expression(paste(italic(phi)[f])),
                                     c(0.0005, 0.1), config_colors)
print(phi_plot$plot)
ggsave("./tests/fluo/validation_phif.png", phi_plot$plot,
       width = 4.5, height = 4.5, scale = 1.25, dpi = 300)
print(phi_plot$stats)

cat("\nValidation plots saved:\n")
cat("  - validation_chl.png\n")
cat("  - validation_ag440.png\n")
cat("  - validation_bbp550.png\n")
cat("  - validation_phif.png\n\n")

# 2. Error distribution
cat("Creating error distribution plots...\n")

error_long <- results_df %>%
  select(config_name, err_chl, err_a_g_440, err_bb_p_550, err_phi_f) %>%
  pivot_longer(cols = starts_with("err_"), names_to = "parameter", values_to = "error_pct") %>%
  mutate(
    parameter = case_when(
      parameter == "err_chl" ~ "Chl-a",
      parameter == "err_a_g_440" ~ "a_g(440)",
      parameter == "err_bb_p_550" ~ "bb_p(550)",
      parameter == "err_phi_f" ~ "φ_f"
    ),
    parameter = factor(parameter, levels = c("Chl-a", "a_g(440)", "bb_p(550)", "φ_f")),
    # Shorten config names and swap for phi_f
    config_short = case_when(
      config_name == "High-res (0.5 nm)" ~ "0.5nm",
      config_name == "Gaussian (3 nm FWHM)" & parameter == "φ_f" ~ "S2+FLEX",
      config_name == "Gaussian (3 nm FWHM)" ~ "3nm",
      config_name == "S2 + FLEX" & parameter == "φ_f" ~ "3nm",
      config_name == "S2 + FLEX" ~ "S2+FLEX"
    ),
    config_factor = factor(config_short, levels = c("0.5nm", "3nm", "S2+FLEX")),
    # Take absolute value for log scale
    abs_error_pct = abs(error_pct)
  )

# Update colors for short names
config_colors_short <- c("0.5nm" = "#377EB8",
                        "3nm" = "#4DAF4A",
                        "S2+FLEX" = "brown4")

p_errors <- ggplot(error_long, aes(x = config_factor, y = abs_error_pct, fill = config_factor)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30", size = 1) +
  geom_violin(alpha = 0.6, scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA,
               color = "black", size = 0.5) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = config_colors_short) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides = "l", size = 0.5) +
  labs(
    x = "",
    y = "Absolute Relative Error (%)",
    fill = "Configuration"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black", size = 1),
    legend.position = "none",
    panel.grid.major = element_line(colour = "grey90", size = 0.3),
    panel.grid.minor = element_line(colour = "grey95", size = 0.2),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.2)
  )

print(p_errors)
ggsave("./tests/fluo/multisensor_error_distribution.png", p_errors, width = 10, height = 8, dpi = 300)


# PART 7: SAVE RESULTS ----

cat("\nSaving results...\n")

# Save results to CSV for easy analysis
write.csv(results_df,
          file = "./tests/fluo/multisensor_inversion_results.csv",
          row.names = FALSE)

# Save R workspace with all objects
save(
  synthetic_params,
  synthetic_data_list,
  config1_nested,
  config2_nested,
  config3_nested,
  results_config1,
  results_config2,
  results_config3,
  results_df,
  summary_stats,
  file = "./tests/fluo/multisensor_synthetic_results.RData"
)

cat("\n============================================================\n")
cat("Multi-sensor synthetic data test complete!\n")
cat("============================================================\n")
cat("\nOutput files:\n")
cat("  - multisensor_chl_comparison.png\n")
cat("  - multisensor_ag_comparison.png\n")
cat("  - multisensor_bbp_comparison.png\n")
cat("  - multisensor_phi_comparison.png\n")
cat("  - multisensor_error_distribution.png\n")
cat("  - multisensor_inversion_results.csv\n")
cat("  - multisensor_synthetic_results.RData\n")
cat("============================================================\n\n")
