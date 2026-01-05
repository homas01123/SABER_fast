library(terra)
library(ncdf4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(SABER)
library(cowplot)
library(ggspatial)  # For scale bar and north arrow
library(ggpubr)     # For extracting legends

# ============================================================================
# BENTHIC REFLECTANCE VISUALIZATION
# ============================================================================
# This script visualizes benthic reflectance from shallow water inversion in two ways:
# 1. FRACTIONAL RGB: Color mixing based on normalized benthic class fractions
# 2. SPECTRAL RGB: RGB composite from mixed benthic reflectance spectra
# ============================================================================

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input files
INPUT_NC <- "./tests/sat/s2_l2a_mp/S2A_MSI_2019_08_28_15_39_27_T19UEQ_L2W.nc"
OUTPUT_NC <- "./tests/sat/s2_l2b_mp/S2A_MSI_2019_08_28_15_39_27_T19UEQ_inversion_shallow_fixed.nc"

# Output directory
OUTPUT_DIR <- "./tests/sat/s2_l2b_mp"

# Study area extent (Lat/Lon WGS84)
study_area_coords <- matrix(c(
  -68.107338, 49.191128,
  -68.529625, 49.193147,
  -68.532028, 49.01423,
  -68.106308, 49.015356,
  -68.107338, 49.191128
), ncol = 2, byrow = TRUE)

# ============================================================================
# APPROACH SELECTION
# ============================================================================
# Set which approach to run:
# 1 = Fractional RGB (color mixing of inverted fractions)
# 2 = Spectral RGB (linear mixing of benthic endmember spectra)
# 3 = Algebraic R_B (direct solution from forward model)
# "all" = Run all three approaches
# ============================================================================

RUN_APPROACH <- 3  # Change this to select approach

# ============================================================================
# APPROACH 3: IOP SOURCE SELECTION
# ============================================================================
# For Approach 3 (Algebraic R_B), select IOP source:
# "bio_optical" = Use bio-optical relationships (OSCs → IOPs via iop_from_oac)
# "in_situ"     = Use field-measured IOPs from nearest station (constant across image)
# ============================================================================

IOP_SOURCE <- "in_situ"  # Options: "bio_optical" or "in_situ"

# In-situ IOP configuration (only used if IOP_SOURCE = "in_situ")
IOP_DATA_DIR <- "Y:/soham_data/IOP/surface_iops/"
IOP_STATION_DB <- "Y:/soham_data/IOP/biogeochemistry_wiseman.csv"
IOP_QC_KILDIR <- "Y:/soham_data/IOP/Kildir_IOP.Process_Log.csv"
IOP_QC_SAUCIER <- "Y:/soham_data/IOP/saucierQC.csv"
MANUAL_IOP_SELECTION <- FALSE
MANUAL_STATION_NAME <- "OUT-F18"

# Benthic class colors (distinct colors for RGB mixing)
BENTHIC_COLORS <- list(
  Eelgrass = c(R = 0.1, G = 0.9, B = 0.2),               # Bright green
  Sand = c(R = 0.5, G = 0.3, B = 0.15),                  # Dark brown
  Mud = c(R = 0.85, G = 0.75, B = 0.0)                   # Saturated dark yellow (high contrast)
)

cat("\n========================================\n")
cat("BENTHIC REFLECTANCE VISUALIZATION\n")

cat(sprintf("Selected approach: %s\n", 
            if(is.numeric(RUN_APPROACH)) paste0("Approach ", RUN_APPROACH) else "All approaches"))
if (RUN_APPROACH == 3 || RUN_APPROACH == "all") {
  cat(sprintf("IOP Source: %s\n", IOP_SOURCE))
}


# ============================================================================
# HELPER FUNCTIONS FOR IN-SITU IOP LOADING
# ============================================================================

# Local implementation of Snell's law (needed for both in-situ and bio-optical IOPs)
snell_law_local <- function(theta_view_deg, theta_sun_deg) {
  # Index of refractions
  n_air <- 1.0    # Air
  n_w <- 1.33     # Water
  
  # Convert degrees to radians
  theta_view <- theta_view_deg * (pi / 180)
  theta_sun <- theta_sun_deg * (pi / 180)
  
  # Angles inside the water (radians)
  view_w <- asin((n_air / n_w) * sin(theta_view))
  sun_w <- asin((n_air / n_w) * sin(theta_sun))
  
  # Fresnel reflectance
  rho_L <- 0.5 * abs(
    ((sin(theta_view - view_w)^2) / (sin(theta_view + view_w)^2)) + 
    ((tan(theta_view - view_w)^2) / (tan(theta_view + view_w)^2))
  )
  
  return(list(view_w = view_w, sun_w = sun_w, rho_L = rho_L))
}

# Function to get nearest field station from image center
get_nearest_station <- function(image_extent, iop_station_db, iop_qc_kildir, iop_qc_saucier, 
                                manual_selection = FALSE, manual_station = "OUT-F18") {
  
  if (manual_selection) {
    cat(sprintf("  Using manually selected station: %s\n", manual_station))
    return(manual_station)
  }
  
  # Get image center in UTM
  center_x <- mean(c(image_extent[1], image_extent[2]))
  center_y <- mean(c(image_extent[3], image_extent[4]))
  
  # Convert to lat/lon (assuming UTM zone 19N, WGS84)
  library(terra)
  center_utm <- matrix(c(center_x, center_y), ncol = 2)
  v <- vect(center_utm, crs = "+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs")
  center_wgs <- project(v, "+proj=longlat +datum=WGS84")
  image_latlon <- crds(center_wgs)
  
  cat(sprintf("  Image center: Lon=%.4f, Lat=%.4f\n", image_latlon[1], image_latlon[2]))
  
  # Load QC-passed stations
  kildir_qc <- read.csv(iop_qc_kildir, header = TRUE)
  saucier_qc <- read.csv(iop_qc_saucier, header = TRUE)
  
  kildir_stations <- kildir_qc$StationID[kildir_qc$ASPH == "Y" & kildir_qc$HS6 == "Y"]
  saucier_stations <- saucier_qc$Station[saucier_qc$station.keep == "TRUE"]
  
  qc_stations <- c(as.character(kildir_stations), as.character(saucier_stations))
  
  # Load station locations
  station_db <- read.csv(iop_station_db, header = TRUE)
  station_db <- station_db[station_db$depth == 0 & station_db$boat != "zodiac", ]
  
  # Filter to QC-passed stations
  station_db <- station_db[station_db$station %in% qc_stations, ]
  
  if (nrow(station_db) == 0) {
    stop("No QC-passed stations found in database")
  }
  
  # Calculate distances
  library(geosphere)
  distances <- distHaversine(
    cbind(station_db$lon, station_db$lat),
    matrix(rep(image_latlon, nrow(station_db)), ncol = 2, byrow = TRUE)
  )
  
  nearest_idx <- which.min(distances)
  nearest_station <- station_db$station[nearest_idx]
  nearest_dist_km <- distances[nearest_idx] / 1000
  
  cat(sprintf("  Nearest station: %s (%.2f km away)\n", nearest_station, nearest_dist_km))
  
  return(as.character(nearest_station))
}

# Function to load in-situ IOPs and calculate attenuation coefficients
load_insitu_iops <- function(station_name, iop_data_dir, wavelengths, 
                             theta_sun, theta_view = 0, water_type = 2) {
  
  cat(sprintf("\nLoading in-situ IOPs for station: %s\n", station_name))
  
  # Find IOP files
  iop_files <- list.files(iop_data_dir, full.names = TRUE)
  
  a_file <- grep(paste0("abs_surf_", station_name, ".csv$"), iop_files, value = TRUE)
  bb_file <- grep(paste0("bb_surf_", station_name, ".csv$"), iop_files, value = TRUE)
  
  if (length(a_file) == 0 || length(bb_file) == 0) {
    stop(sprintf("IOP files not found for station %s in %s", station_name, iop_data_dir))
  }
  
  cat(sprintf("  Absorption file: %s\n", basename(a_file)))
  cat(sprintf("  Backscatter file: %s\n", basename(bb_file)))
  
  # Load IOP data
  a_data <- read.csv(a_file, header = TRUE)
  bb_data <- read.csv(bb_file, header = TRUE)
  
  # Extract wavelength and values (adjust column names as needed)
  # Assuming columns: wave, at_w (for absorption) and bbp (for backscatter)
  a_wave <- a_data$wave
  a_vals <- a_data$at_w
  
  bb_wave <- bb_data$wave
  bb_vals <- bb_data$bbp
  
  # Get pure water IOPs
  a_w <- pure_water_iop(wavelength = wavelengths)$a
  
  # Pure water backscattering (from WISE code)
  if (water_type == 1) {
    b1 <- 0.00144  # Case 1
  } else {
    b1 <- 0.00111  # Case 2
  }
  lambda1 <- 500
  bb_w <- b1 * (wavelengths / lambda1)^(-4.32)
  
  # Interpolate measured non-water absorption to S2 wavelengths
  a_non_water <- approx(x = a_wave, y = a_vals, xout = wavelengths, method = "linear", rule = 2)$y
  a_non_water[is.na(a_non_water)] <- 0
  
  # Calculate total absorption
  a_total <- a_w + a_non_water
  
  # Fit power-law model to backscatter data (bb = bb_550 * (lambda/550)^-gamma)
  # Detect sensor type
  if (length(bb_wave) == 6) {
    sensor_wl <- c(394, 420, 470, 532, 620, 700)
    cat("  Sensor: HS-6 VSF\n")
  } else if (length(bb_wave) == 9) {
    sensor_wl <- c(412, 440, 488, 510, 532, 595, 650, 676, 715)
    cat("  Sensor: BB-9 VSF\n")
  } else {
    stop("Unknown backscatter sensor configuration")
  }
  
  # Fit power-law: bb = b * (lambda/555)^z
  x <- 555 / sensor_wl
  y <- bb_vals
  y[y < 0 | y > 0.1] <- NA
  
  if (all(is.na(y))) {
    stop("All backscatter values are NA or out of range")
  }
  
  # Remove NAs for fitting
  valid <- !is.na(y)
  x_fit <- x[valid]
  y_fit <- y[valid]
  
  # Fit power-law model
  tryCatch({
    model <- nls(y_fit ~ b * x_fit^z, 
                 start = list(b = y_fit[which.min(abs(sensor_wl[valid] - 555))], z = 1),
                 control = list(maxiter = 100, warnOnly = TRUE))
    
    bb_555 <- coef(model)[1]
    gamma <- coef(model)[2]
    
    cat(sprintf("  Power-law fit: bb_555 = %.6f, gamma = %.4f\n", bb_555, gamma))
    
  }, error = function(e) {
    cat("  Warning: Power-law fit failed, using mean values\n")
    bb_555 <- mean(y_fit, na.rm = TRUE)
    gamma <- 1.0
  })
  
  # Calculate bb at S2 wavelengths using power-law
  bb_non_water <- bb_555 * ((wavelengths / 555)^(-gamma))
  bb_total <- bb_w + bb_non_water
  
  # Apply Snell's law for underwater angles
  snell_result <- snell_law_local(theta_view, theta_sun)
  view_w <- snell_result$view_w
  sun_w <- snell_result$sun_w
  
  # Calculate attenuation coefficients
  ext <- a_total + bb_total
  omega_b <- bb_total / (a_total + bb_total + 1e-10)
  
  if (water_type == 1) {
    k0 <- 1.0395
  } else {
    k0 <- 1.0546
  }
  
  Kd <- k0 * (ext / cos(sun_w))
  KuW <- (ext / cos(view_w)) * ((1 + omega_b)^3.5421) * (1 - (0.2786 / cos(sun_w)))
  KuB <- (ext / cos(view_w)) * ((1 + omega_b)^2.2658) * (1 - (0.0577 / cos(sun_w)))
  
  # Calculate deep water rrs using f_rs formula
  if (water_type == 1) {
    f_rs <- 0.095
  } else {
    f_rs <- 0.0512 * (1 + (4.6659 * omega_b) + (-7.8387 * (omega_b^2)) + (5.4571 * (omega_b^3))) *
      (1 + (0.1098 / cos(sun_w))) * (1 + (0.4021 / cos(view_w)))
  }
  
  rrs_deep <- f_rs * omega_b
  
  cat(sprintf("  ✓ Calculated attenuation coefficients for %d wavelengths\n", length(wavelengths)))
  
  return(list(
    a = a_total,
    bb = bb_total,
    Kd = Kd,
    KuW = KuW,
    KuB = KuB,
    rrs_deep = rrs_deep,
    station = station_name
  ))
}

# ============================================================================
# 1. LOAD INPUT DATA (S2 L2W - for wavelengths)
# ============================================================================

cat("\nLoading input S2 L2W data...\n")
cat(sprintf("  File: %s\n", INPUT_NC))

input_rast <- rast(INPUT_NC)

# Extract wavelengths from Rrs band names
rrs_bands <- grep("^Rrs_", names(input_rast), value = TRUE)
wavelengths <- as.numeric(gsub("^Rrs_", "", rrs_bands))

cat(sprintf("  Found %d Rrs bands\n", length(rrs_bands)))
cat("  Wavelengths: ", paste(wavelengths, collapse = ", "), " nm\n")

# ============================================================================
# 2. LOAD INVERSION OUTPUT DATA
# ============================================================================

cat("\nLoading inversion output data...\n")
cat(sprintf("  File: %s\n", OUTPUT_NC))

output_rast <- rast(OUTPUT_NC)

cat(sprintf("  Raster dimensions: %d rows x %d cols\n", 
            nrow(output_rast), ncol(output_rast)))
cat(sprintf("  Total layers: %d\n", nlyr(output_rast)))

# Extract normalized benthic class fractions
benthic_normalized_layers <- grep("_normalized$", names(output_rast), value = TRUE)

cat("\nNormalized benthic layers found:\n")
for (layer in benthic_normalized_layers) {
  cat(sprintf("  - %s\n", layer))
}

# Extract specific layers
eelgrass_norm <- output_rast[["r_rs_b_Eelgrass_2019_normalized"]]
sand_norm <- output_rast[["r_rs_b_Sand_2019_normalized"]]
mud_norm <- output_rast[["r_rs_b_Mud_2019_normalized"]]

# Verify normalization (should sum to ~1)
benthic_sum <- eelgrass_norm + sand_norm + mud_norm
sum_values <- values(benthic_sum, na.rm = TRUE)
cat(sprintf("\nNormalization check: Mean sum = %.6f (should be ~1.0)\n", 
            mean(sum_values, na.rm = TRUE)))

# ============================================================================
# 3. APPROACH 1: FRACTIONAL RGB COLOR MIXING
# ============================================================================

if (RUN_APPROACH == 1 || RUN_APPROACH == "all") {
  
cat("\n========================================\n")
cat("APPROACH 1: FRACTIONAL RGB COLOR MIXING\n")

cat("Creating color-mixed benthic classification image...\n")

# Create RGB raster by mixing colors based on fractions
# Step 1: Mix colors (composition)
R_mixed <- (eelgrass_norm * BENTHIC_COLORS$Eelgrass["R"] + 
            sand_norm * BENTHIC_COLORS$Sand["R"] + 
            mud_norm * BENTHIC_COLORS$Mud["R"])

G_mixed <- (eelgrass_norm * BENTHIC_COLORS$Eelgrass["G"] + 
            sand_norm * BENTHIC_COLORS$Sand["G"] + 
            mud_norm * BENTHIC_COLORS$Mud["G"])

B_mixed <- (eelgrass_norm * BENTHIC_COLORS$Eelgrass["B"] + 
            sand_norm * BENTHIC_COLORS$Sand["B"] + 
            mud_norm * BENTHIC_COLORS$Mud["B"])

# Step 2: Apply gamma correction for better contrast
# Gamma < 1 brightens mid-tones, Gamma > 1 darkens them
# Here we use the max fraction to modulate intensity
max_fraction <- max(c(eelgrass_norm, sand_norm, mud_norm))

cat("  Applying intensity scaling with gamma correction...\n")

# Use a power function to scale intensity more naturally
# This prevents over-brightening of light colors and over-darkening of dark colors
gamma <- 0.5  # Adjust this value: <1 = brighter, >1 = darker
intensity_scale <- max_fraction ^ gamma

# Step 3: Scale RGB by intensity
R_channel <- R_mixed * intensity_scale
G_channel <- G_mixed * intensity_scale
B_channel <- B_mixed * intensity_scale

# Stack RGB channels
rgb_fractional <- c(R_channel, G_channel, B_channel)
names(rgb_fractional) <- c("R", "G", "B")

cat("RGB fractional raster created\n")

# Convert to dataframe for ggplot
df_fractional <- as.data.frame(rgb_fractional, xy = TRUE, na.rm = TRUE)

cat(sprintf("  Valid pixels: %d\n", nrow(df_fractional)))

# Clip values to 0-1 range
df_fractional$R <- pmin(pmax(df_fractional$R, 0), 1)
df_fractional$G <- pmin(pmax(df_fractional$G, 0), 1)
df_fractional$B <- pmin(pmax(df_fractional$B, 0), 1)

# Calculate brightness (intensity) for diagnostics
df_fractional$brightness <- (df_fractional$R + df_fractional$G + df_fractional$B) / 3

cat(sprintf("  Brightness range: %.3f - %.3f (mean: %.3f)\n",
            min(df_fractional$brightness, na.rm = TRUE),
            max(df_fractional$brightness, na.rm = TRUE),
            mean(df_fractional$brightness, na.rm = TRUE)))

# Create RGB plot WITHOUT equalization (preserve true colors)
cat("Creating fractional RGB plot with scale bar and north arrow...\n")

p_fractional <- ggplot(df_fractional, aes(x = x, y = y)) +
  geom_raster(aes(fill = rgb(R, G, B))) +
  scale_fill_identity() +
  scale_y_reverse() +  # Flip Y-axis to match raster orientation
  coord_fixed() +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12, color = 'black', angle = 15, hjust = 1),
    axis.text.y = element_text(size = 12, color = 'black'),
    axis.title.x = element_text(size = 20, margin = margin(t = 15)),
    axis.title.y = element_text(size = 20, margin = margin(t = 15)),
    axis.ticks.length = unit(.25, "cm"),
    panel.grid.major = element_line(color = "grey50", linewidth = 0.5, linetype = "dotted"),
    panel.grid.minor = element_line(color = "grey80", linewidth = 0.2),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering)

# Create legend data frame
legend_df <- data.frame(
  Class = c("Eelgrass", "Sand", "Mud"),
  R = c(BENTHIC_COLORS$Eelgrass["R"], BENTHIC_COLORS$Sand["R"], BENTHIC_COLORS$Mud["R"]),
  G = c(BENTHIC_COLORS$Eelgrass["G"], BENTHIC_COLORS$Sand["G"], BENTHIC_COLORS$Mud["G"]),
  B = c(BENTHIC_COLORS$Eelgrass["B"], BENTHIC_COLORS$Sand["B"], BENTHIC_COLORS$Mud["B"]),
  x = 1,
  y = 1:3
)

# # Create color legend plot
# p_legend <- ggplot(legend_df, aes(x = x, y = y)) +
#   geom_tile(aes(fill = rgb(R, G, B)), width = 0.8, height = 0.8) +
#   geom_text(aes(label = Class), x = 1.5, hjust = 0, size = 5) +
#   scale_fill_identity() +
#   xlim(0.5, 3) +
#   ylim(0.5, 3.5) +
#   labs(title = "Benthic Classes") +
#   theme_void() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#     plot.margin = margin(10, 10, 10, 10)
#   )

# Combine plot with legend
# combined_fractional <- plot_grid(p_fractional, p_legend, ncol = 2, 
#                                  rel_widths = c(0.75, 0.25))

# Save plot
output_fractional <- file.path(OUTPUT_DIR, "benthic_fractional_rgb.png")
ggsave(output_fractional, p_fractional, width = 12, height = 12, dpi = 300, units = "in")

cat(sprintf("✓ Saved fractional RGB plot: %s\n", basename(output_fractional)))

}  # End of Approach 1

# ============================================================================
# 4. APPROACH 2: MIXED BENTHIC REFLECTANCE SPECTRA
# ============================================================================

if (RUN_APPROACH == 2 || RUN_APPROACH == "all") {

cat("\n========================================\n")
cat("APPROACH 2: MIXED BENTHIC REFLECTANCE SPECTRA\n")

cat("Computing mixed benthic reflectance for each pixel...\n")

# Get the benthic reflectance data directly from SABER package
data("r_rs_b_egsl", package = "SABER", envir = environment())

# Extract benthic fraction data as dataframe
df_benthic <- as.data.frame(c(eelgrass_norm, sand_norm, mud_norm), 
                            xy = TRUE, na.rm = TRUE)
names(df_benthic)[3:5] <- c("r_rs_b_Eelgrass_2019", "r_rs_b_Sand_2019", "r_rs_b_Mud_2019")

cat(sprintf("  Processing %d pixels...\n", nrow(df_benthic)))

cat("  Computing mixed benthic reflectance (vectorized)...\n")

# Filter for our three classes and prepare spectra
sand_data <- r_rs_b_egsl %>% 
  filter(class == "Sand_2019") %>%
  filter(wavelength >= 400 & wavelength <= 760) %>%
  arrange(wavelength)

eelgrass_data <- r_rs_b_egsl %>% 
  filter(class == "Eelgrass_2019") %>%
  filter(wavelength >= 400 & wavelength <= 760) %>%
  arrange(wavelength)

mud_data <- r_rs_b_egsl %>% 
  filter(class == "Mud_2019") %>%
  filter(wavelength >= 400 & wavelength <= 760) %>%
  arrange(wavelength)

benthic_wavelengths <- sand_data$wavelength
cat(sprintf("  Benthic data has %d wavelengths\n", length(benthic_wavelengths)))

# Get S2 wavelengths in valid range
wavelengths_filtered <- wavelengths[wavelengths >= 400 & wavelengths <= 760]
cat(sprintf("  Interpolating to %d S2 wavelengths\n", length(wavelengths_filtered)))

# Interpolate each benthic class spectrum to S2 wavelengths
sand_spectrum <- approx(sand_data$wavelength, sand_data$r_rs_b_mean, 
                       xout = wavelengths_filtered, rule = 2)$y
eelgrass_spectrum <- approx(eelgrass_data$wavelength, eelgrass_data$r_rs_b_mean, 
                            xout = wavelengths_filtered, rule = 2)$y
mud_spectrum <- approx(mud_data$wavelength, mud_data$r_rs_b_mean, 
                      xout = wavelengths_filtered, rule = 2)$y

# Vectorized computation: linear mixing
# r_rs_b_mixed = fraction_sand * sand_spectrum + fraction_eelgrass * eelgrass_spectrum + fraction_mud * mud_spectrum
r_rs_b_mixed <- (
  outer(df_benthic$r_rs_b_Sand_2019, mud_spectrum) +
  outer(df_benthic$r_rs_b_Eelgrass_2019, eelgrass_spectrum) +
  outer(df_benthic$r_rs_b_Mud_2019, sand_spectrum )
)

cat("  ✓ Mixed reflectance computed (vectorized)\n")

# Add spectral data to dataframe
colnames(r_rs_b_mixed) <- paste0("Rrs_", wavelengths_filtered)
df_benthic_spectral <- cbind(df_benthic[, c("x", "y")], r_rs_b_mixed)

cat(sprintf("  Valid spectral pixels: %d\n", nrow(df_benthic_spectral)))

# ============================================================================
# 5. CREATE RGB COMPOSITE FROM BENTHIC REFLECTANCE SPECTRA
# ============================================================================

cat("\nCreating RGB composite from benthic reflectance spectra...\n")

# Select wavelengths closest to RGB (Red ~665, Green ~560, Blue ~490)
# Find closest available wavelengths in filtered range
rgb_target_wl <- c(665, 560, 490)
rgb_actual_wl <- sapply(rgb_target_wl, function(target) {
  wavelengths_filtered[which.min(abs(wavelengths_filtered - target))]
})

cat(sprintf("  Using wavelengths: R=%d, G=%d, B=%d nm\n", 
            rgb_actual_wl[1], rgb_actual_wl[2], rgb_actual_wl[3]))

# Extract RGB bands
R_band <- df_benthic_spectral[[paste0("Rrs_", rgb_actual_wl[1])]]
G_band <- df_benthic_spectral[[paste0("Rrs_", rgb_actual_wl[2])]]
B_band <- df_benthic_spectral[[paste0("Rrs_", rgb_actual_wl[3])]]

cat(sprintf("  R-band range: %.4f - %.4f\n", min(R_band, na.rm=TRUE), max(R_band, na.rm=TRUE)))
cat(sprintf("  G-band range: %.4f - %.4f\n", min(G_band, na.rm=TRUE), max(G_band, na.rm=TRUE)))
cat(sprintf("  B-band range: %.4f - %.4f\n", min(B_band, na.rm=TRUE), max(B_band, na.rm=TRUE)))

# Simple linear scaling to 0-1 based on typical benthic reflectance range
# Benthic Rrs_b typically ranges from 0 to ~0.2 sr^-1
# Use a fixed upper limit for consistent visualization
max_rrs_b <- 0.15  # Typical maximum benthic reflectance

# Apply gamma correction for better contrast (brightens the image)
gamma_spectral <- 16  # <1 = brighter, >1 = darker

df_benthic_spectral$R_norm <- pmin((R_band / max_rrs_b) ^ gamma_spectral, 1.0)
df_benthic_spectral$G_norm <- pmin((G_band / max_rrs_b) ^ gamma_spectral, 1.0)
df_benthic_spectral$B_norm <- pmin((B_band / max_rrs_b) ^ gamma_spectral, 1.0)

cat(sprintf("  Applied gamma correction (γ=%.2f) for better contrast\n", gamma_spectral))

cat(sprintf("  Normalized RGB ranges: R[%.3f-%.3f], G[%.3f-%.3f], B[%.3f-%.3f]\n",
            min(df_benthic_spectral$R_norm, na.rm=TRUE), max(df_benthic_spectral$R_norm, na.rm=TRUE),
            min(df_benthic_spectral$G_norm, na.rm=TRUE), max(df_benthic_spectral$G_norm, na.rm=TRUE),
            min(df_benthic_spectral$B_norm, na.rm=TRUE), max(df_benthic_spectral$B_norm, na.rm=TRUE)))

# Create RGB plot
cat("Creating spectral RGB plot with three intensity scales...\n")

# Main RGB composite plot
p_spectral_rgb <- ggplot(df_benthic_spectral, aes(x = x, y = y)) +
  geom_raster(aes(fill = rgb(R_norm, G_norm, B_norm))) +
  scale_fill_identity() +
  scale_y_reverse() +  # Flip Y-axis to match raster orientation
  coord_fixed() +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12, color = 'black', angle = 15, hjust = 1),
    axis.text.y = element_text(size = 12, color = 'black'),
    axis.title.x = element_text(size = 20, margin = margin(t = 15)),
    axis.title.y = element_text(size = 20, margin = margin(t = 15)),
    axis.ticks.length = unit(.25, "cm"),
    panel.grid.major = element_line(color = "grey50", linewidth = 0.5, linetype = "dotted"),
    panel.grid.minor = element_line(color = "grey80", linewidth = 0.2),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering)

# Create three individual band plots for legend
# R band (665 nm) - Red scale
p_R_band <- ggplot(df_benthic_spectral, aes(x = x, y = y, fill = R_band)) +
  geom_raster() +
  scale_fill_gradient(
    low = "black", 
    high = "red",
    limits = c(0, max(R_band, na.rm = TRUE)),
    name = expression(R[B](665~nm)),
    breaks = seq(0, max(R_band, na.rm = TRUE), length.out = 5),
    labels = function(x) sprintf("%.3f", x)
  ) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.5, "cm")
  )

# G band (560 nm) - Green scale
p_G_band <- ggplot(df_benthic_spectral, aes(x = x, y = y, fill = G_band)) +
  geom_raster() +
  scale_fill_gradient(
    low = "black", 
    high = "green",
    limits = c(0, max(G_band, na.rm = TRUE)),
    name = expression(R[B](560~nm)),
    breaks = seq(0, max(G_band, na.rm = TRUE), length.out = 5),
    labels = function(x) sprintf("%.3f", x)
  ) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.5, "cm")
  )

# B band (490 nm) - Blue scale
p_B_band <- ggplot(df_benthic_spectral, aes(x = x, y = y, fill = B_band)) +
  geom_raster() +
  scale_fill_gradient(
    low = "black", 
    high = "blue",
    limits = c(0, max(B_band, na.rm = TRUE)),
    name = expression(R[B](490~nm)),
    breaks = seq(0, max(B_band, na.rm = TRUE), length.out = 5),
    labels = function(x) sprintf("%.3f", x)
  ) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.5, "cm")
  )

# Extract legends
library(ggpubr)
legend_R <- get_legend(p_R_band)
legend_G <- get_legend(p_G_band)
legend_B <- get_legend(p_B_band)

# Combine legends vertically (stacked on right side)
combined_legends <- plot_grid(legend_R, legend_G, legend_B, ncol = 1, align = "v")

# Combine RGB plot with legends (legends on right side)
p_spectral <- plot_grid(
  p_spectral_rgb, 
  combined_legends, 
  ncol = 2, 
  rel_widths = c(0.80, 0.20)
)

# Save spectral RGB plot
output_spectral <- file.path(OUTPUT_DIR, "benthic_spectral_rgb.png")
ggsave(output_spectral, p_spectral, width = 14, height = 12, dpi = 300, units = "in")

cat(sprintf("✓ Saved spectral RGB plot: %s\n", basename(output_spectral)))

}  # End of Approach 2

# ============================================================================
# 5. APPROACH 3: ALGEBRAIC R_B RETRIEVAL FROM FORWARD MODEL
# ============================================================================

if (RUN_APPROACH == 3 || RUN_APPROACH == "all") {

cat("\n========================================\n")
cat("APPROACH 3: ALGEBRAIC R_B RETRIEVAL\n")

cat("Computing benthic reflectance by algebraically solving forward model...\n")

# Load inverted IOP parameters from output NetCDF
cat("\nLoading inverted IOP parameters...\n")

chl_raster <- output_rast[["chl"]]
a_g_440_raster <- output_rast[["a_g_440"]]
a_g_s_raster <- output_rast[["a_g_s"]]
bb_p_550_raster <- output_rast[["bb_p_550"]]
bb_p_gamma_raster <- output_rast[["bb_p_gamma"]]

# ============================================================================
# STEP 1: DEFINE STUDY AREA AND CROP ALL INPUTS
# ============================================================================
cat("\nDefining study area polygon...\n")

# Study area coordinates (provided by user)
study_area_coords <- matrix(c(
  -68.107338, 49.191128,  # NE
  -68.529625, 49.193147,  # NW
  -68.532028, 49.01423,   # SW
  -68.106308, 49.015356,  # SE
  -68.107338, 49.191128   # Close polygon
), ncol = 2, byrow = TRUE)

study_area <- vect(study_area_coords, type = "polygons", crs = "EPSG:4326")

# Convert to UTM 19N
s2_crs <- crs(input_rast)
study_area_utm <- project(study_area, s2_crs)

cat(sprintf("Study area extent (UTM): [%.1f, %.1f] x [%.1f, %.1f]\n",
            ext(study_area_utm)[1], ext(study_area_utm)[2],
            ext(study_area_utm)[3], ext(study_area_utm)[4]))

# ============================================================================
# STEP 2: CROP S2 L2W TO STUDY AREA
# ============================================================================
cat("\nCropping S2 L2W to study area...\n")

# Crop S2 to study area extent
s2_cropped <- crop(input_rast, study_area_utm)

cat(sprintf("S2 cropped: %d rows x %d cols, Res=[%.3f, %.3f]\n",
            nrow(s2_cropped), ncol(s2_cropped),
            res(s2_cropped)[1], res(s2_cropped)[2]))

# Use cropped S2 as reference grid
reference_grid <- s2_cropped[[1]]

# Extract Rrs for all bands (400-860 nm range)
rrs_bands_all <- grep("^Rrs_", names(s2_cropped), value = TRUE)
wavelengths_all <- as.numeric(gsub("^Rrs_", "", rrs_bands_all))
wavelengths_range <- wavelengths_all[wavelengths_all >= 400 & wavelengths_all <= 860]
rrs_bands_range <- paste0("Rrs_", wavelengths_range)

cat(sprintf("Processing %d wavelengths: %s nm\n", 
            length(wavelengths_range), 
            paste(wavelengths_range, collapse = ", ")))

rrs_aligned <- s2_cropped[[rrs_bands_range]]

# ============================================================================
# STEP 3: GEOREFERENCE INVERSION OUTPUT USING CROPPED S2
# ============================================================================
cat("\nGeoreferencing inversion output...\n")

# Check if inversion has CRS
inv_crs <- crs(output_rast)
if (inv_crs == "") {
  cat("Inversion has no CRS - assigning from cropped S2...\n")
  
  # Assign CRS to inversion rasters
  crs(chl_raster) <- s2_crs
  crs(a_g_440_raster) <- s2_crs
  crs(a_g_s_raster) <- s2_crs
  crs(bb_p_550_raster) <- s2_crs
  crs(bb_p_gamma_raster) <- s2_crs
  crs(output_rast) <- s2_crs
}

# Resample inversion IOPs to match cropped S2 grid
cat("Resampling inversion IOPs to reference grid...\n")
chl_aligned <- resample(chl_raster, reference_grid, method = "bilinear")
a_g_440_aligned <- resample(a_g_440_raster, reference_grid, method = "bilinear")
a_g_s_aligned <- resample(a_g_s_raster, reference_grid, method = "bilinear")
bb_p_550_aligned <- resample(bb_p_550_raster, reference_grid, method = "bilinear")
bb_p_gamma_aligned <- resample(bb_p_gamma_raster, reference_grid, method = "bilinear")

cat(sprintf("IOPs aligned: %d rows x %d cols\n", nrow(chl_aligned), ncol(chl_aligned)))

# ============================================================================
# STEP 4: CROP AND RESAMPLE BATHYMETRY TO STUDY AREA
# ============================================================================
cat("\nProcessing bathymetry...\n")

# Load bathymetry
bathy_file <- "c:/R/SABER_fast/tests/sat/bathymetry_mp/Asli_bhara_hua_NONNA_bath.tiff"
bathy_full <- rast(bathy_file)

# Crop bathymetry to study area (in lat/lon)
bathy_cropped <- crop(bathy_full, study_area)

# Reproject and resample to reference grid
cat("Reprojecting bathymetry to UTM and resampling...\n")
bathy_aligned <- project(bathy_cropped, reference_grid, method = "bilinear")

cat(sprintf("Bathymetry aligned: %d rows x %d cols\n", nrow(bathy_aligned), ncol(bathy_aligned)))

# ============================================================================
# STEP 5: PROCESS SOLAR ZENITH ANGLE
# ============================================================================
cat("\nProcessing solar zenith angle...\n")

# Load solar zenith angle
sza_file <- "./tests/sat/s2_l2a_mp/S2A_2019_08_28_sza.tif"
if (file.exists(sza_file)) {
  sza_raster <- rast(sza_file)
  
  # Crop and resample to reference grid
  sza_cropped <- crop(sza_raster, study_area_utm)
  sza_aligned <- resample(sza_cropped, reference_grid, method = "bilinear")
  
  cat(sprintf("SZA aligned: %d rows x %d cols\n", nrow(sza_aligned), ncol(sza_aligned)))
} else {
  cat("SZA file not found, using constant: 41.13 deg\n")
  sza_aligned <- reference_grid * 0 + 41.13
}

# ============================================================================
# STEP 6: VERIFY ALL GRIDS MATCH
# ============================================================================
cat("\nVerifying grid alignment...\n")

all_match <- TRUE

# Check IOPs
if (!all(ext(chl_aligned) == ext(reference_grid))) {
  cat("ERROR: IOP extent mismatch!\n")
  all_match <- FALSE
}

# Check Rrs
if (!all(ext(rrs_aligned) == ext(reference_grid))) {
  cat("ERROR: Rrs extent mismatch!\n")
  all_match <- FALSE
}

# Check bathymetry
if (!all(ext(bathy_aligned) == ext(reference_grid))) {
  cat("ERROR: Bathymetry extent mismatch!\n")
  all_match <- FALSE
}

# Check SZA
if (!all(ext(sza_aligned) == ext(reference_grid))) {
  cat("ERROR: SZA extent mismatch!\n")
  all_match <- FALSE
}

if (all_match) {
  cat("All grids aligned successfully!\n")
  cat(sprintf("Common extent: [%.1f, %.1f] x [%.1f, %.1f]\n",
              ext(reference_grid)[1], ext(reference_grid)[2], 
              ext(reference_grid)[3], ext(reference_grid)[4]))
  cat(sprintf("Common resolution: [%.3f, %.3f] m\n", 
              res(reference_grid)[1], res(reference_grid)[2]))
  cat(sprintf("Dimensions: %d rows x %d cols\n",
              nrow(reference_grid), ncol(reference_grid)))
} else {
  stop("Grid alignment failed!")
}

# ============================================================================
# STEP 7: FILTER TO SHALLOW WATER (0.5-10m)
# ============================================================================
cat("\nFiltering to shallow water (0.5-10m)...\n")

# Create masks
depth_mask <- bathy_aligned >= 0.5 & bathy_aligned <= 10.0
iop_mask <- !is.na(chl_aligned) & !is.na(a_g_440_aligned) & !is.na(bb_p_550_aligned)

# Combined mask
valid_mask <- depth_mask & iop_mask

cat(sprintf("Shallow water pixels: %d (%.1f%%)\n",
            sum(values(depth_mask), na.rm = TRUE),
            100 * sum(values(depth_mask), na.rm = TRUE) / length(values(depth_mask))))

cat(sprintf("Pixels with valid IOPs: %d (%.1f%%)\n",
            sum(values(iop_mask), na.rm = TRUE),
            100 * sum(values(iop_mask), na.rm = TRUE) / length(values(iop_mask))))

cat(sprintf("Valid pixels (shallow + IOPs): %d (%.1f%%)\n",
            sum(values(valid_mask), na.rm = TRUE),
            100 * sum(values(valid_mask), na.rm = TRUE) / length(values(valid_mask))))

# Extract valid pixel indices
pixels_valid <- which(values(valid_mask))

if (length(pixels_valid) == 0) {
  stop("No pixels found with shallow water + valid IOPs!")
}

cat(sprintf("Processing %d valid pixels...\n", length(pixels_valid)))

# Extract coordinates
coords <- xyFromCell(reference_grid, pixels_valid)

# Extract all variables for valid pixels
df_algebraic <- data.frame(
  x = coords[, 1],
  y = coords[, 2],
  h_w = values(bathy_aligned)[pixels_valid],
  chl = values(chl_aligned)[pixels_valid],
  a_g_440 = values(a_g_440_aligned)[pixels_valid],
  a_g_s = values(a_g_s_aligned)[pixels_valid],
  bb_p_550 = values(bb_p_550_aligned)[pixels_valid],
  bb_p_gamma = values(bb_p_gamma_aligned)[pixels_valid],
  theta_sun = values(sza_aligned)[pixels_valid]  # Solar zenith angle
)

# Extract Rrs for all wavelengths
for (i in seq_along(rrs_bands_range)) {
  band_name <- rrs_bands_range[i]
  df_algebraic[[band_name]] <- values(rrs_aligned[[i]])[pixels_valid]
}

cat(sprintf("  ✓ Extracted data for %d pixels\n", nrow(df_algebraic)))

# ============================================================================
# LOAD IN-SITU IOPs (if IOP_SOURCE = "in_situ")
# ============================================================================

INSITU_IOPS <- NULL  # Will hold pre-computed IOPs if using in-situ data

if (IOP_SOURCE == "in_situ") {
  cat("\n========================================\n")
  cat("LOADING IN-SITU FIELD IOPs\n")
  
  
  # Get image extent for station finding
  image_ext <- ext(input_rast)
  
  # Find nearest station
  station_name <- get_nearest_station(
    image_extent = as.vector(image_ext),
    iop_station_db = IOP_STATION_DB,
    iop_qc_kildir = IOP_QC_KILDIR,
    iop_qc_saucier = IOP_QC_SAUCIER,
    manual_selection = MANUAL_IOP_SELECTION,
    manual_station = MANUAL_STATION_NAME
  )
  
  # Load IOPs and calculate attenuation coefficients
  # Use average theta_sun from data
  avg_theta_sun <- mean(df_algebraic$theta_sun, na.rm = TRUE)
  cat(sprintf("  Using average solar zenith angle: %.2f deg\n", avg_theta_sun))
  
  INSITU_IOPS <- load_insitu_iops(
    station_name = station_name,
    iop_data_dir = IOP_DATA_DIR,
    wavelengths = wavelengths_range,
    theta_sun = avg_theta_sun,
    theta_view = 0,
    water_type = 2
  )
  
  cat("\n  ✓ In-situ IOPs loaded and ready for algebraic R_B retrieval\n")
  cat(sprintf("  Station: %s (constant IOPs across entire image)\n", INSITU_IOPS$station))
  
} else {
  cat("\n========================================\n")
  cat("USING BIO-OPTICAL MODEL IOPs\n")
  
  cat("  IOPs will be calculated from OSCs (chl, a_g, bb_p) using iop_from_oac\n")
  cat("  IOPs vary spatially across the image\n")
}

# Algebraic solution for each pixel
cat("\nComputing algebraic R_B for each pixel...\n")
cat("This may take a few minutes...\n")

# Function to compute R_B algebraically for one pixel
compute_rb_algebraic <- function(pixel_data, wavelengths, insitu_iops = NULL) {
  
  # Get Rrs(0+) observations (above-water) for all wavelengths
  rrs_0p_obs <- as.numeric(pixel_data[paste0("Rrs_", wavelengths)])
  
  # CRITICAL: Convert above-water Rrs(0+) to subsurface rrs(0-)
  rrs_obs <- rrs_0p_to_0m(rrs_0p_obs)
  
  # Extract parameters
  h_w <- pixel_data$h_w
  theta_sun <- pixel_data$theta_sun
  theta_view <- 0  # Nadir viewing
  water_type <- 2  # Case 2 water
  
  # ========================================================================
  # OPTION 1: Use in-situ field IOPs (constant across image)
  # ========================================================================
  if (!is.null(insitu_iops)) {
    # Use pre-computed IOPs from nearest field station
    a <- insitu_iops$a
    bb <- insitu_iops$bb
    Kd <- insitu_iops$Kd
    KuW <- insitu_iops$KuW
    KuB <- insitu_iops$KuB
    rrs_dp <- insitu_iops$rrs_deep
    
  # ========================================================================
  # OPTION 2: Use bio-optical relationships (OSCs → IOPs)
  # ========================================================================
  } else {
    # Extract OSC parameters
    chl <- pixel_data$chl
    a_g_440 <- pixel_data$a_g_440
    a_g_s <- pixel_data$a_g_s
    bb_p_550 <- pixel_data$bb_p_550
    bb_p_gamma <- pixel_data$bb_p_gamma
    
    # Prepare parameter vector for IOP calculation
    par <- c(
      chl = chl,
      a_g_440 = a_g_440,
      a_g_s_g = a_g_s,
      bb_p_550 = bb_p_550,
      bb_p_gamma = bb_p_gamma
    )
    
    # Calculate IOPs using iop_from_oac
    iops <- tryCatch({
      iop_from_oac(wavelength = wavelengths, par = par)
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(iops) || !is.list(iops)) {
      return(rep(NA_real_, length(wavelengths)))
    }
    
    # Calculate deep water rrs
    rrs_dp <- tryCatch({
      forward_am03(
        wavelength = wavelengths,
        iop = iops,
        water_type = water_type,
        theta_sun = theta_sun,
        theta_view = theta_view,
        h_w = NULL,  # Deep water (no bottom)
        r_b = NULL
      )
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(rrs_dp) || any(is.na(rrs_dp))) {
      return(rep(NA_real_, length(wavelengths)))
    }
    
    # Apply Snell's law for underwater angles
    snell_result <- snell_law_local(theta_view, theta_sun)
    view_w <- snell_result$view_w
    sun_w <- snell_result$sun_w
    
    # Extract IOPs
    a <- iops$a
    bb <- iops$bb
    
    # Calculate single scattering albedo and extinction
    omega_b <- bb / (a + bb + 1e-10)
    ext <- a + bb
    
    # Calculate attenuation coefficients
    if (water_type == 1) {
      k0 <- 1.0395  # Case 1
    } else {
      k0 <- 1.0546  # Case 2
    }
    
    Kd <- k0 * (ext / cos(sun_w))
    KuW <- (ext / cos(view_w)) * ((1 + omega_b)^3.5421) * (1 - (0.2786 / cos(sun_w)))
    KuB <- (ext / cos(view_w)) * ((1 + omega_b)^2.2658) * (1 - (0.0577 / cos(sun_w)))
  }
  
  # ========================================================================
  # ALGEBRAIC SOLUTION (same for both IOP sources)
  # ========================================================================
  # WISE-style formula with Ars coefficients
  Ars1 <- 1.1576  # Parametric coefficient for shallow water
  Ars2 <- 1.0389  # Parametric coefficient for shallow water
  
  # Calculate R_B using WISE formula
  # R_B = [rrs_obs - rrs_dp * (1 - Ars1 * exp(-h*(Kd+KuW)))] / [Ars2 * exp(-h*(Kd+KuB))]
  exp_term_W <- exp(-h_w * (Kd + KuW))
  exp_term_B <- exp(-h_w * (Kd + KuB))
  
  R_B <- (rrs_obs - rrs_dp * (1 - Ars1 * exp_term_W)) / (Ars2 * exp_term_B + 1e-10)
  
  # NO CONSTRAINTS - Return all R_B values including negative and >1
  # This allows for diagnostic analysis of model performance
  # Only mask NaN and Inf
  R_B[is.nan(R_B) | is.infinite(R_B)] <- NA
  
  # Return R_B with additional diagnostic info
  result <- list(
    R_B = R_B,
    rrs_obs = rrs_obs,
    rrs_dp = rrs_dp,
    Kd = Kd,
    KuW = KuW,
    KuB = KuB
  )
  
  return(result)
}

# Apply to all pixels - vectorized processing
cat("  Processing pixels (vectorized)...\n")

# Use lapply for vectorized processing (can be parallelized if needed)
rb_results <- lapply(1:nrow(df_algebraic), function(i) {
  if (i %% 1000 == 0) {
    cat(sprintf("  Processing pixel %d/%d (%.1f%%)...\n", 
                i, nrow(df_algebraic), 100 * i / nrow(df_algebraic)))
  }
  compute_rb_algebraic(df_algebraic[i, ], wavelengths_range, insitu_iops = INSITU_IOPS)
})

cat("  ✓ Algebraic R_B computation complete\n")

# Extract R_B values and diagnostic info
rb_matrix <- do.call(rbind, lapply(rb_results, function(x) x$R_B))
colnames(rb_matrix) <- paste0("R_B_", wavelengths_range)

# Store diagnostic spectral data for sample pixels
rrs_obs_matrix <- do.call(rbind, lapply(rb_results, function(x) x$rrs_obs))
rrs_dp_matrix <- do.call(rbind, lapply(rb_results, function(x) x$rrs_dp))
colnames(rrs_obs_matrix) <- paste0("rrs_obs_", wavelengths_range)
colnames(rrs_dp_matrix) <- paste0("rrs_dp_", wavelengths_range)

df_algebraic_rb <- cbind(df_algebraic[, c("x", "y", "h_w", "chl", "a_g_440")], 
                         rb_matrix, rrs_obs_matrix, rrs_dp_matrix)

# Diagnostic: Check how many pixels have negative/invalid R_B
cat("\nDiagnostic - R_B validity by category:\n")
for (wl in wavelengths_range) {
  col_name <- paste0("R_B_", wl)
  total <- nrow(df_algebraic_rb)
  valid <- sum(!is.na(df_algebraic_rb[[col_name]]))
  negative_count <- sum(df_algebraic_rb[[col_name]] < 0, na.rm = TRUE)
  above_one <- sum(df_algebraic_rb[[col_name]] > 1, na.rm = TRUE)
  between_0_1 <- sum(df_algebraic_rb[[col_name]] >= 0 & df_algebraic_rb[[col_name]] <= 1, na.rm = TRUE)
  cat(sprintf("  %s nm: Total=%d, Valid=%d (%.1f%%), Negative=%d (%.1f%%), >1=%d (%.1f%%), 0-1=%d (%.1f%%)\n", 
              wl, total, valid, 100*valid/total, 
              negative_count, 100*negative_count/total,
              above_one, 100*above_one/total,
              between_0_1, 100*between_0_1/total))
}

# ============================================================================
# SPECTRAL DIAGNOSTIC PLOTS - Sample pixels from each R_B category
# ============================================================================

cat("\n========================================\n")
cat("SPECTRAL DIAGNOSTIC ANALYSIS\n")

cat("Creating spectral comparison plots for different R_B ranges...\n")

# Categorize pixels based on R_B at 560 nm (green band - typically most stable)
ref_wl <- 560
rb_ref_col <- paste0("R_B_", ref_wl)

df_algebraic_rb$rb_category <- NA
df_algebraic_rb$rb_category[df_algebraic_rb[[rb_ref_col]] < 0] <- "negative"
df_algebraic_rb$rb_category[df_algebraic_rb[[rb_ref_col]] >= 0 & df_algebraic_rb[[rb_ref_col]] <= 1] <- "valid"
df_algebraic_rb$rb_category[df_algebraic_rb[[rb_ref_col]] > 1] <- "above_one"

cat(sprintf("  Reference wavelength: %d nm\n", ref_wl))
cat(sprintf("  Negative R_B: %d pixels\n", sum(df_algebraic_rb$rb_category == "negative", na.rm=TRUE)))
cat(sprintf("  Valid R_B (0-1): %d pixels\n", sum(df_algebraic_rb$rb_category == "valid", na.rm=TRUE)))
cat(sprintf("  R_B > 1: %d pixels\n", sum(df_algebraic_rb$rb_category == "above_one", na.rm=TRUE)))

# Sample pixels from each category
set.seed(42)
n_samples <- 5  # Number of samples per category

sample_negative <- df_algebraic_rb %>% 
  filter(rb_category == "negative") %>% 
  sample_n(min(n_samples, n()))

sample_valid <- df_algebraic_rb %>% 
  filter(rb_category == "valid") %>% 
  sample_n(min(n_samples, n()))

sample_above <- df_algebraic_rb %>% 
  filter(rb_category == "above_one") %>% 
  sample_n(min(n_samples, n()))

cat(sprintf("\nSampled %d pixels per category for spectral plots\n", n_samples))

# Function to create spectral comparison plot for one pixel
create_spectral_plot <- function(pixel_row, wavelengths, title_suffix = "") {
  # Extract rrs_obs and rrs_dp
  rrs_obs <- as.numeric(pixel_row[paste0("rrs_obs_", wavelengths)])
  rrs_dp <- as.numeric(pixel_row[paste0("rrs_dp_", wavelengths)])
  rb_vals <- as.numeric(pixel_row[paste0("R_B_", wavelengths)])
  
  # Filter to 400-750nm range for better visualization
  wl_filter <- wavelengths >= 400 & wavelengths <= 750
  wavelengths_plot <- wavelengths[wl_filter]
  rrs_obs_plot <- rrs_obs[wl_filter]
  rrs_dp_plot <- rrs_dp[wl_filter]
  rb_vals_plot <- rb_vals[wl_filter]
  
  # Create dataframe for plotting
  df_spec <- data.frame(
    wavelength = rep(wavelengths_plot, 3),
    rrs = c(rrs_obs_plot, rrs_dp_plot, rb_vals_plot),
    type = factor(rep(c("Observed (Rrs 0m)", "Deep water (no benthic)", "Retrieved R_B"), 
                      each = length(wavelengths_plot)),
                  levels = c("Observed (Rrs 0m)", "Deep water (no benthic)", "Retrieved R_B"))
  )
  
  # Get depth and other info
  h_w <- pixel_row$h_w
  rb_560 <- pixel_row[[paste0("R_B_", ref_wl)]]
  
  p <- ggplot(df_spec, aes(x = wavelength, y = rrs, color = type, linetype = type)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Observed (Rrs 0m)" = "black", 
                                   "Deep water (no benthic)" = "blue", 
                                   "Retrieved R_B" = "red")) +
    scale_linetype_manual(values = c("Observed (Rrs 0m)" = "solid", 
                                      "Deep water (no benthic)" = "dashed", 
                                      "Retrieved R_B" = "solid")) +
    scale_y_log10(labels = scales::label_number(accuracy = 0.0001)) +  # LOG SCALE for better visualization
    annotation_logticks(sides = "l") +  # Add log tick marks
    labs(
      title = sprintf("Spectral Comparison %s", title_suffix),
      subtitle = sprintf("Depth=%.2fm, R_B@%dnm=%.3f, Category=%s", 
                        h_w, ref_wl, rb_560, pixel_row$rb_category),
      x = "Wavelength (nm)",
      y = "Reflectance (sr^-1) [LOG SCALE]",
      color = "Spectrum",
      linetype = "Spectrum"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom",
      panel.grid.minor = element_line(color = "gray90"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    geom_hline(yintercept = 1e-10, linetype = "dotted", color = "gray50")
  
  return(p)
}

# Create plots for each category
plots_negative <- list()
plots_valid <- list()
plots_above <- list()

cat("\nGenerating spectral plots...\n")

# Only create plots if samples exist
if (nrow(sample_negative) > 0) {
  for (i in 1:nrow(sample_negative)) {
    plots_negative[[i]] <- create_spectral_plot(sample_negative[i, ], wavelengths_range, 
                                                 sprintf("- Negative R_B (#%d)", i))
  }
  cat(sprintf("  Created %d plots for negative R_B\n", length(plots_negative)))
}

if (nrow(sample_valid) > 0) {
  for (i in 1:nrow(sample_valid)) {
    plots_valid[[i]] <- create_spectral_plot(sample_valid[i, ], wavelengths_range, 
                                             sprintf("- Valid R_B (#%d)", i))
  }
  cat(sprintf("  Created %d plots for valid R_B\n", length(plots_valid)))
}

if (nrow(sample_above) > 0) {
  for (i in 1:nrow(sample_above)) {
    plots_above[[i]] <- create_spectral_plot(sample_above[i, ], wavelengths_range, 
                                             sprintf("- R_B > 1 (#%d)", i))
  }
  cat(sprintf("  Created %d plots for R_B > 1\n", length(plots_above)))
}

# Combine plots into grids
if (length(plots_negative) > 0) {
  cat(sprintf("  Combining %d negative R_B plots...\n", length(plots_negative)))
  combined_negative <- plot_grid(plotlist = plots_negative, ncol = 2)
  output_negative <- file.path(OUTPUT_DIR, "spectral_diagnostic_negative.png")
  ggsave(output_negative, combined_negative, width = 16, height = 4*ceiling(length(plots_negative)/2), 
         dpi = 300, units = "in")
  cat(sprintf("✓ Saved negative R_B spectral plots: %s\n", basename(output_negative)))
} else {
  cat("  No negative R_B samples to plot\n")
}

if (length(plots_valid) > 0) {
  cat(sprintf("  Combining %d valid R_B plots...\n", length(plots_valid)))
  combined_valid <- plot_grid(plotlist = plots_valid, ncol = 2)
  output_valid <- file.path(OUTPUT_DIR, "spectral_diagnostic_valid.png")
  ggsave(output_valid, combined_valid, width = 16, height = 4*ceiling(length(plots_valid)/2), 
         dpi = 300, units = "in")
  cat(sprintf("✓ Saved valid R_B spectral plots: %s\n", basename(output_valid)))
} else {
  cat("  No valid R_B samples to plot\n")
}

if (length(plots_above) > 0) {
  cat(sprintf("  Combining %d R_B>1 plots...\n", length(plots_above)))
  combined_above <- plot_grid(plotlist = plots_above, ncol = 2)
  output_above <- file.path(OUTPUT_DIR, "spectral_diagnostic_above1.png")
  ggsave(output_above, combined_above, width = 16, height = 4*ceiling(length(plots_above)/2), 
         dpi = 300, units = "in")
  cat(sprintf("✓ Saved R_B>1 spectral plots: %s\n", basename(output_above)))
} else {
  cat("  No R_B>1 samples to plot\n")
}

cat(sprintf("  Valid pixels (at least one band): %d\n", sum(!is.na(rb_matrix[, 1]))))

# ============================================================================
# CREATE INPUT Rrs RGB COMPOSITE (BEFORE ALGEBRAIC INVERSION)
# ============================================================================
cat("\n========================================\n")
cat("INPUT Rrs RGB COMPOSITE\n")

cat("Creating RGB composite from input Rrs (0m corrected)...\n")

# Select wavelengths closest to RGB
rgb_target_wl <- c(665, 560, 490)
rgb_actual_wl_rrs <- sapply(rgb_target_wl, function(target) {
  wavelengths_range[which.min(abs(wavelengths_range - target))]
})

cat(sprintf("  Using wavelengths: R=%d, G=%d, B=%d nm\n", 
            rgb_actual_wl_rrs[1], rgb_actual_wl_rrs[2], rgb_actual_wl_rrs[3]))

# Extract Rrs RGB bands
rrs_r_col <- paste0("Rrs_", rgb_actual_wl_rrs[1])
rrs_g_col <- paste0("Rrs_", rgb_actual_wl_rrs[2])
rrs_b_col <- paste0("Rrs_", rgb_actual_wl_rrs[3])

# Get Rrs values (subsurface corrected)
R_band_rrs <- df_algebraic_rb[[paste0("rrs_obs_", rgb_actual_wl_rrs[1])]]
G_band_rrs <- df_algebraic_rb[[paste0("rrs_obs_", rgb_actual_wl_rrs[2])]]
B_band_rrs <- df_algebraic_rb[[paste0("rrs_obs_", rgb_actual_wl_rrs[3])]]

cat(sprintf("  R-band range: %.6f - %.6f sr^-1\n", min(R_band_rrs, na.rm=TRUE), max(R_band_rrs, na.rm=TRUE)))
cat(sprintf("  G-band range: %.6f - %.6f sr^-1\n", min(G_band_rrs, na.rm=TRUE), max(G_band_rrs, na.rm=TRUE)))
cat(sprintf("  B-band range: %.6f - %.6f sr^-1\n", min(B_band_rrs, na.rm=TRUE), max(B_band_rrs, na.rm=TRUE)))

# Use percentile-based normalization
max_R_rrs <- quantile(R_band_rrs, 0.98, na.rm = TRUE)
max_G_rrs <- quantile(G_band_rrs, 0.98, na.rm = TRUE)
max_B_rrs <- quantile(B_band_rrs, 0.98, na.rm = TRUE)

cat(sprintf("  Using 98th percentile for normalization: R=%.6f, G=%.6f, B=%.6f\n", 
            max_R_rrs, max_G_rrs, max_B_rrs))

# Normalize to [0,1] and apply gamma correction for better visualization
df_rrs_rgb <- df_algebraic_rb %>%
  mutate(
    R_norm = pmin(pmax((R_band_rrs / max_R_rrs)^0.5, 0), 1),  # Gamma 0.5 for brightness
    G_norm = pmin(pmax((G_band_rrs / max_G_rrs)^0.5, 0), 1),
    B_norm = pmin(pmax((B_band_rrs / max_B_rrs)^0.5, 0), 1)
  ) %>%
  filter(!is.na(R_norm), !is.na(G_norm), !is.na(B_norm),
         is.finite(R_norm), is.finite(G_norm), is.finite(B_norm))

cat(sprintf("  Valid RGB pixels: %d\n", nrow(df_rrs_rgb)))

# Create RGB plot with reference style
p_rrs_rgb <- ggplot(df_rrs_rgb, aes(x = x, y = y)) +
  geom_raster(aes(fill = rgb(R_norm, G_norm, B_norm))) +
  scale_fill_identity() +
  coord_equal() +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  labs(title = "Input Rrs (Subsurface Corrected)",
       x = "Easting (m)", y = "Northing (m)") +
  theme_minimal()

output_rrs_rgb <- file.path(OUTPUT_DIR, "input_rrs_rgb.png")
ggsave(output_rrs_rgb, p_rrs_rgb, width = 12, height = 10, dpi = 300, units = "in", bg = "white")
cat(sprintf("Saved input Rrs RGB: %s\n", basename(output_rrs_rgb)))

# ============================================================================
# CREATE ALGEBRAIC R_B RGB COMPOSITE (AFTER INVERSION)
# ============================================================================
cat("\n========================================\n")
cat("ALGEBRAIC R_B RGB COMPOSITE\n")

cat("Creating RGB composite from algebraic R_B...\n")

# Select wavelengths closest to RGB
rgb_actual_wl_alg <- sapply(rgb_target_wl, function(target) {
  wavelengths_range[which.min(abs(wavelengths_range - target))]
})

cat(sprintf("  Using wavelengths: R=%d, G=%d, B=%d nm\n", 
            rgb_actual_wl_alg[1], rgb_actual_wl_alg[2], rgb_actual_wl_alg[3]))

R_band_alg <- df_algebraic_rb[[paste0("R_B_", rgb_actual_wl_alg[1])]]
G_band_alg <- df_algebraic_rb[[paste0("R_B_", rgb_actual_wl_alg[2])]]
B_band_alg <- df_algebraic_rb[[paste0("R_B_", rgb_actual_wl_alg[3])]]

cat(sprintf("  R-band range: %.4f - %.4f\n", min(R_band_alg, na.rm=TRUE), max(R_band_alg, na.rm=TRUE)))
cat(sprintf("  G-band range: %.4f - %.4f\n", min(G_band_alg, na.rm=TRUE), max(G_band_alg, na.rm=TRUE)))
cat(sprintf("  B-band range: %.4f - %.4f\n", min(B_band_alg, na.rm=TRUE), max(B_band_alg, na.rm=TRUE)))

# Print percentiles to understand distribution
cat(sprintf("  R-band percentiles: 10%%=%.4f, 50%%=%.4f, 90%%=%.4f\n", 
            quantile(R_band_alg, 0.10, na.rm=TRUE), 
            quantile(R_band_alg, 0.50, na.rm=TRUE),
            quantile(R_band_alg, 0.90, na.rm=TRUE)))
cat(sprintf("  G-band percentiles: 10%%=%.4f, 50%%=%.4f, 90%%=%.4f\n", 
            quantile(G_band_alg, 0.10, na.rm=TRUE), 
            quantile(G_band_alg, 0.50, na.rm=TRUE),
            quantile(G_band_alg, 0.90, na.rm=TRUE)))
cat(sprintf("  B-band percentiles: 10%%=%.4f, 50%%=%.4f, 90%%=%.4f\n", 
            quantile(B_band_alg, 0.10, na.rm=TRUE), 
            quantile(B_band_alg, 0.50, na.rm=TRUE),
            quantile(B_band_alg, 0.90, na.rm=TRUE)))

# Use percentile-based normalization for better dynamic range
# Use 95th percentile as max to avoid outliers
max_R <- quantile(R_band_alg, 0.95, na.rm = TRUE)
max_G <- quantile(G_band_alg, 0.95, na.rm = TRUE)
max_B <- quantile(B_band_alg, 0.95, na.rm = TRUE)

cat(sprintf("  Using 95th percentile for normalization: R=%.4f, G=%.4f, B=%.4f\n", 
            max_R, max_G, max_B))

# Clip to [0, max_percentile] range and normalize
df_algebraic_rgb <- df_algebraic_rb %>%
  mutate(
    R_norm = pmin(pmax(R_band_alg / max_R, 0), 1),
    G_norm = pmin(pmax(G_band_alg / max_G, 0), 1),
    B_norm = pmin(pmax(B_band_alg / max_B, 0), 1)
  ) %>%
  filter(!is.na(R_norm), !is.na(G_norm), !is.na(B_norm),
         is.finite(R_norm), is.finite(G_norm), is.finite(B_norm))

cat(sprintf("  Valid RGB pixels after filtering: %d\n", nrow(df_algebraic_rgb)))

# Check coordinate ranges for orientation verification
cat(sprintf("  X range: %.2f to %.2f\n", min(df_algebraic_rgb$x), max(df_algebraic_rgb$x)))
cat(sprintf("  Y range: %.2f to %.2f\n", min(df_algebraic_rgb$y), max(df_algebraic_rgb$y)))

# Create plot with reference style
p_algebraic_rgb <- ggplot(df_algebraic_rgb, aes(x = x, y = y)) +
  geom_raster(aes(fill = rgb(R_norm, G_norm, B_norm))) +
  scale_fill_identity() +
  coord_equal() +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  labs(title = "Algebraic R_B (Benthic Reflectance)",
       x = "Easting (m)", y = "Northing (m)") +
  theme_minimal()

# Save plot
output_algebraic <- file.path(OUTPUT_DIR, "benthic_algebraic_rgb.png")
ggsave(output_algebraic, p_algebraic_rgb, width = 12, height = 12, dpi = 300, units = "in")

cat(sprintf("Saved algebraic RGB: %s\n", basename(output_algebraic)))

# Save algebraic R_B to NetCDF with proper georeferencing
cat("\nSaving algebraic R_B to NetCDF...\n")

# Convert to rasters (use reference grid)
r_template_alg <- reference_grid  # Match inversion output extent
rb_rasters_alg <- list()

cell_indices_alg <- cellFromXY(r_template_alg, 
                               cbind(df_algebraic_rb$x, df_algebraic_rb$y))

for (wl in wavelengths_range) {
  col_name <- paste0("R_B_", wl)
  if (col_name %in% names(df_algebraic_rb)) {
    r_temp <- rast(r_template_alg)
    vals <- rep(NA_real_, ncell(r_temp))
    
    valid_cells <- !is.na(cell_indices_alg) & cell_indices_alg > 0 & cell_indices_alg <= ncell(r_temp)
    vals[cell_indices_alg[valid_cells]] <- df_algebraic_rb[[col_name]][valid_cells]
    
    values(r_temp) <- vals
    names(r_temp) <- col_name
    rb_rasters_alg[[col_name]] <- r_temp
  }
}

if (length(rb_rasters_alg) > 0) {
  r_rb_stack_alg <- rast(rb_rasters_alg)
  
  # Write to NetCDF using ncdf4 with PROPER GEOREFERENCING (matching input L2W structure)
  nx_alg <- ncol(r_rb_stack_alg)
  ny_alg <- nrow(r_rb_stack_alg)
  ext_vals_alg <- ext(r_rb_stack_alg)
  
  # Create coordinate arrays (cell centers) - MATCHING INPUT FILE
  xvals_alg <- seq(ext_vals_alg[1], ext_vals_alg[2], length.out = nx_alg)
  yvals_alg <- seq(ext_vals_alg[3], ext_vals_alg[4], length.out = ny_alg)
  
  # Define dimensions with proper CF-compliant names (EXACT MATCH to input)
  xdim_alg <- ncdim_def("x", "m", xvals_alg, 
                        longname = "x coordinate of projection",
                        create_dimvar = TRUE)
  ydim_alg <- ncdim_def("y", "m", yvals_alg, 
                        longname = "y coordinate of projection",
                        create_dimvar = TRUE)
  
  # Define CRS/grid_mapping variable for GDAL/QGIS compatibility
  crs_var <- ncvar_def(
    name = "transverse_mercator",
    units = "",
    dim = list(),  # Scalar variable
    missval = NULL,
    longname = "CRS definition",
    prec = "integer"
  )
  
  # Define lat/lon coordinate variables (2D grids)
  lon_var <- ncvar_def(
    name = "lon",
    units = "degrees_east",
    dim = list(xdim_alg, ydim_alg),
    missval = -9999,
    longname = "longitude",
    prec = "double"
  )
  
  lat_var <- ncvar_def(
    name = "lat",
    units = "degrees_north",
    dim = list(xdim_alg, ydim_alg),
    missval = -9999,
    longname = "latitude",
    prec = "double"
  )
  
  # Define R_B variables
  var_list_alg <- list(crs_var, lon_var, lat_var)  # Start with CRS and coordinate variables
  
  for (i in 1:nlyr(r_rb_stack_alg)) {
    var_name <- names(r_rb_stack_alg)[i]
    wl_num <- gsub("R_B_", "", var_name)
    var_list_alg[[i+3]] <- ncvar_def(
      name = var_name,
      units = "1",
      dim = list(xdim_alg, ydim_alg),
      missval = -9999,
      longname = sprintf("Benthic reflectance at %s nm", wl_num),
      prec = "float"
    )
  }
  
  output_algebraic_nc <- file.path(OUTPUT_DIR, "benthic_algebraic_Rb.nc")
  ncout_alg <- nc_create(output_algebraic_nc, var_list_alg, force_v4 = TRUE)
  
  # Write CRS variable (scalar, value = 1)
  ncvar_put(ncout_alg, "transverse_mercator", 1)
  
  # Calculate and write lat/lon grids
  cat("Computing lat/lon coordinate grids...\n")
  xy_grid <- expand.grid(x = xvals_alg, y = yvals_alg)
  lonlat <- project(as.matrix(xy_grid), from = crs(r_rb_stack_alg), to = "EPSG:4326")
  lon_matrix <- matrix(lonlat[, 1], nrow = nx_alg, ncol = ny_alg)
  lat_matrix <- matrix(lonlat[, 2], nrow = nx_alg, ncol = ny_alg)
  
  ncvar_put(ncout_alg, "lon", lon_matrix)
  ncvar_put(ncout_alg, "lat", lat_matrix)
  
  cat(sprintf("Lon range: %.4f to %.4f\n", min(lon_matrix), max(lon_matrix)))
  cat(sprintf("Lat range: %.4f to %.4f\n", min(lat_matrix), max(lat_matrix)))
  
  # Add CF-1.7 grid_mapping attributes to CRS variable
  ncatt_put(ncout_alg, "transverse_mercator", "grid_mapping_name", "transverse_mercator")
  ncatt_put(ncout_alg, "transverse_mercator", "longitude_of_central_meridian", -69.0)
  ncatt_put(ncout_alg, "transverse_mercator", "latitude_of_projection_origin", 0.0)
  ncatt_put(ncout_alg, "transverse_mercator", "scale_factor_at_central_meridian", 0.9996)
  ncatt_put(ncout_alg, "transverse_mercator", "false_easting", 500000.0)
  ncatt_put(ncout_alg, "transverse_mercator", "false_northing", 0.0)
  ncatt_put(ncout_alg, "transverse_mercator", "semi_major_axis", 6378137.0)
  ncatt_put(ncout_alg, "transverse_mercator", "inverse_flattening", 298.257223563)
  ncatt_put(ncout_alg, "transverse_mercator", "spatial_ref", crs(r_rb_stack_alg, proj=TRUE))
  ncatt_put(ncout_alg, "transverse_mercator", "GeoTransform", 
            paste(ext_vals_alg[1], 10, 0, ext_vals_alg[4], 0, -10))
  
  # Write R_B data
  for (i in 1:nlyr(r_rb_stack_alg)) {
    vals <- values(r_rb_stack_alg[[i]], mat = TRUE)
    vals[is.na(vals)] <- -9999
    ncvar_put(ncout_alg, var_list_alg[[i+3]], vals)
  }
  
  # ========================================================================
  # ADD COMPREHENSIVE GEOREFERENCING ATTRIBUTES (MATCHING INPUT L2W FORMAT)
  # ========================================================================
  
  # Get CRS from template
  crs_proj4 <- crs(r_rb_stack_alg, proj=TRUE)
  
  # Global attributes - MATCHING ACOLITE L2W STRUCTURE
  ncatt_put(ncout_alg, 0, "Conventions", "CF-1.7")
  ncatt_put(ncout_alg, 0, "title", "SABER Algebraic Benthic Reflectance (R_B)")
  ncatt_put(ncout_alg, 0, "source", "Sentinel-2 MSI L2W - Algebraic retrieval from forward model")
  ncatt_put(ncout_alg, 0, "institution", "SABER Project")
  ncatt_put(ncout_alg, 0, "date_created", as.character(Sys.time()))
  ncatt_put(ncout_alg, 0, "product_type", "NetCDF")
  ncatt_put(ncout_alg, 0, "sensor", "S2A_MSI")
  
  # Projection information (MATCHING INPUT)
  ncatt_put(ncout_alg, 0, "proj4_string", crs_proj4)
  ncatt_put(ncout_alg, 0, "projection_key", "transverse_mercator")
  
  # Spatial extent (MATCHING INPUT FORMAT)
  ncatt_put(ncout_alg, 0, "xrange", c(ext_vals_alg[1], ext_vals_alg[2]))
  ncatt_put(ncout_alg, 0, "yrange", c(ext_vals_alg[3], ext_vals_alg[4]))
  ncatt_put(ncout_alg, 0, "pixel_size", c(10, -10))  # Match S2 10m resolution
  
  # Data dimensions (MATCHING INPUT)
  ncatt_put(ncout_alg, 0, "data_dimensions", c(ny_alg, nx_alg))
  ncatt_put(ncout_alg, 0, "data_elements", ny_alg * nx_alg)
  
  # NetCDF projection flag (MATCHING INPUT) - use integer instead of logical
  ncatt_put(ncout_alg, 0, "netcdf_projection", as.integer(1))
  
  # Add coordinate variable attributes (CF-1.7 compliant)
  ncatt_put(ncout_alg, "x", "standard_name", "projection_x_coordinate")
  ncatt_put(ncout_alg, "x", "long_name", "x coordinate of projection")
  ncatt_put(ncout_alg, "x", "units", "m")
  ncatt_put(ncout_alg, "x", "axis", "X")
  
  ncatt_put(ncout_alg, "y", "standard_name", "projection_y_coordinate")
  ncatt_put(ncout_alg, "y", "long_name", "y coordinate of projection")
  ncatt_put(ncout_alg, "y", "units", "m")
  ncatt_put(ncout_alg, "y", "axis", "Y")
  
  # Add attributes to R_B variables
  for (i in 1:nlyr(r_rb_stack_alg)) {
    var_name <- names(r_rb_stack_alg)[i]
    wl_num <- gsub("R_B_", "", var_name)
    ncatt_put(ncout_alg, var_name, "grid_mapping", "transverse_mercator")  # CRITICAL for GDAL/QGIS
    ncatt_put(ncout_alg, var_name, "standard_name", "benthic_reflectance")
    ncatt_put(ncout_alg, var_name, "wavelength", as.numeric(wl_num))
    ncatt_put(ncout_alg, var_name, "wavelength_units", "nm")
    ncatt_put(ncout_alg, var_name, "valid_range", c(-10, 10))  # Allow all values for diagnostics
    ncatt_put(ncout_alg, var_name, "comment", "No constraints applied - includes negative and >1 values for diagnostics")
  }
  
  # Close the file
  nc_close(ncout_alg)
  
  cat(sprintf("✓ Saved algebraic R_B NetCDF: %s\n", basename(output_algebraic_nc)))
  cat(sprintf("  Bands: %d wavelengths\n", nlyr(r_rb_stack_alg)))
  cat("  Georeferencing: COMPLETE (CF-1.7 compliant, matching input L2W structure)\n")
  cat(sprintf("  Extent: X[%.2f, %.2f], Y[%.2f, %.2f]\n", 
              ext_vals_alg[1], ext_vals_alg[2], ext_vals_alg[3], ext_vals_alg[4]))
  cat(sprintf("  Dimensions: %d rows x %d cols\n", ny_alg, nx_alg))
  cat(sprintf("  Projection: %s\n", crs_proj4))
}

}  # End of Approach 3

# ============================================================================
# 6. CREATE COMPARISON PLOT
# ============================================================================

if (RUN_APPROACH == "all") {

cat("\n========================================\n")
cat("CREATING COMPARISON PLOT\n")


# Load original S2 RGB for comparison
s2_rgb_file <- file.path(OUTPUT_DIR, "inv_shallow_RGB_input.png")

if (file.exists(s2_rgb_file)) {
  # Create combined comparison (would need to load the image)
  cat("  Note: For full comparison, manually arrange plots\n")
} else {
  cat("  Original RGB not found, creating standalone comparison\n")
}

# Create three-way comparison of all approaches
combined_comparison <- plot_grid(
  p_fractional + theme(legend.position = "none") + labs(title = "A) Fractional RGB"),
  p_spectral_rgb + labs(title = "B) Spectral RGB"),
  p_algebraic_rgb + labs(title = "C) Algebraic R_B"),
  ncol = 3,
  label_size = 14
)

output_comparison <- file.path(OUTPUT_DIR, "benthic_comparison_all.png")
ggsave(output_comparison, combined_comparison, width = 24, height = 8, dpi = 300, units = "in")

cat(sprintf("✓ Saved comparison plot: %s\n", basename(output_comparison)))

}  # End of Comparison plot

# ============================================================================
# 7. SAVE BENTHIC SPECTRAL RASTERS (OPTIONAL)
# ============================================================================

if (RUN_APPROACH == 2 || RUN_APPROACH == "all") {

cat("\nSaving benthic spectral rasters to NetCDF...\n")

# Convert spectral dataframe back to rasters (VECTORIZED)
spectral_rasters <- list()

# Create template raster
r_template <- rast(eelgrass_norm)

# Get all cell indices for the xy coordinates (vectorized)
cell_indices <- cellFromXY(r_template, 
                          cbind(df_benthic_spectral$x, df_benthic_spectral$y))

# Process all wavelengths
for (wl in wavelengths_filtered) {
  col_name <- paste0("Rrs_", wl)
  if (col_name %in% names(df_benthic_spectral)) {
    # Create empty raster
    r_temp <- rast(r_template)
    
    # Initialize with NA
    vals <- rep(NA_real_, ncell(r_temp))
    
    # Vectorized assignment: directly assign values to cells
    valid_cells <- !is.na(cell_indices) & cell_indices > 0 & cell_indices <= ncell(r_temp)
    vals[cell_indices[valid_cells]] <- df_benthic_spectral[[col_name]][valid_cells]
    
    values(r_temp) <- vals
    names(r_temp) <- paste0("r_rs_b_mixed_", wl)
    spectral_rasters[[paste0("r_rs_b_mixed_", wl)]] <- r_temp
  }
}

if (length(spectral_rasters) > 0) {
  # Stack all spectral bands
  r_spectral_stack <- rast(spectral_rasters)
  
  cat(sprintf("Creating NetCDF with %d bands using ncdf4...\n", length(spectral_rasters)))
  
  # Get dimensions from the raster
  nx <- ncol(r_spectral_stack)
  ny <- nrow(r_spectral_stack)
  ext_vals <- ext(r_spectral_stack)
  
  # Create dimension variables
  xvals <- seq(ext_vals[1], ext_vals[2], length.out = nx)
  yvals <- seq(ext_vals[3], ext_vals[4], length.out = ny)
  
  # Define dimensions with standard names for SEADAS compatibility
  xdim <- ncdim_def("x", "meters", xvals)
  ydim <- ncdim_def("y", "meters", yvals)
  
  # Create wavelength dimension for band information
  wl_vals <- as.numeric(gsub("r_rs_b_mixed_", "", names(r_spectral_stack)))
  wl_dim <- ncdim_def("wavelength", "nanometers", wl_vals)
  
  # Create a list to store all variable definitions
  var_list <- list()
  
  cat("  Creating variable definitions:\n")
  for (i in 1:nlyr(r_spectral_stack)) {
    var_name <- names(r_spectral_stack)[i]
    cat(sprintf("    %d. %s\n", i, var_name))
    
    var_list[[i]] <- ncvar_def(
      name = var_name,
      units = "sr^-1",
      dim = list(xdim, ydim),
      missval = -9999,
      longname = paste("Mixed benthic reflectance at", gsub("r_rs_b_mixed_", "", var_name), "nm"),
      prec = "float"
    )
  }
  
  # Create the NetCDF file
  output_spectral_nc <- file.path(OUTPUT_DIR, "benthic_mixed_spectra.nc")
  ncout <- nc_create(output_spectral_nc, var_list, force_v4 = TRUE)
  
  # Write data for each variable
  cat("  Writing data to NetCDF:\n")
  for (i in 1:nlyr(r_spectral_stack)) {
    var_name <- names(r_spectral_stack)[i]
    cat(sprintf("    Writing %s...\n", var_name))
    
    # Extract values and convert to matrix
    vals <- values(r_spectral_stack[[i]], mat = TRUE)
    vals[is.na(vals)] <- -9999
    
    # Transpose matrix (terra uses row-major, ncdf4 uses column-major)
    vals <- t(vals)
    
    # Write to NetCDF
    ncvar_put(ncout, var_list[[i]], vals)
  }
  
  # Add global attributes
  ncatt_put(ncout, 0, "title", "SABER Mixed Benthic Reflectance Spectra")
  ncatt_put(ncout, 0, "source", "Sentinel-2 shallow water inversion")
  ncatt_put(ncout, 0, "date_created", as.character(Sys.time()))
  ncatt_put(ncout, 0, "crs", as.character(crs(r_spectral_stack)))
  ncatt_put(ncout, 0, "wavelengths", paste(wavelengths_filtered, collapse = ", "))
  ncatt_put(ncout, 0, "Conventions", "CF-1.6")
  ncatt_put(ncout, 0, "instrument", "MSI")
  ncatt_put(ncout, 0, "platform", "Sentinel-2A")
  ncatt_put(ncout, 0, "product_name", "benthic_mixed_spectra")
  
  # Add coordinate system attributes
  ncatt_put(ncout, "x", "standard_name", "projection_x_coordinate")
  ncatt_put(ncout, "x", "long_name", "x coordinate of projection")
  ncatt_put(ncout, "x", "axis", "X")
  
  ncatt_put(ncout, "y", "standard_name", "projection_y_coordinate")
  ncatt_put(ncout, "y", "long_name", "y coordinate of projection")
  ncatt_put(ncout, "y", "axis", "Y")
  
  # Close the file
  nc_close(ncout)
  
  cat(sprintf("✓ Saved spectral rasters: %s\n", basename(output_spectral_nc)))
  cat(sprintf("  Bands: %d wavelengths\n", nlyr(r_spectral_stack)))
} else {
  cat("  Warning: No spectral rasters created\n")
}

}  # End of Spectral rasters save

# ============================================================================
# 8. SUMMARY STATISTICS
# ============================================================================

cat("\n========================================\n")
cat("BENTHIC REFLECTANCE SUMMARY\n")


# Fractional composition statistics (only if Approach 1 or 2 ran)
if (exists("df_benthic") && nrow(df_benthic) > 0) {
  cat("\nFractional Composition (normalized):\n")
  frac_stats <- df_benthic %>%
    summarize(
      Eelgrass_mean = mean(r_rs_b_Eelgrass_2019, na.rm = TRUE),
      Eelgrass_sd = sd(r_rs_b_Eelgrass_2019, na.rm = TRUE),
      Sand_mean = mean(r_rs_b_Sand_2019, na.rm = TRUE),
      Sand_sd = sd(r_rs_b_Sand_2019, na.rm = TRUE),
      Mud_mean = mean(r_rs_b_Mud_2019, na.rm = TRUE),
      Mud_sd = sd(r_rs_b_Mud_2019, na.rm = TRUE)
    )
  print(frac_stats)

# Dominant class statistics
df_benthic_dominant <- df_benthic %>%
  mutate(
    dominant_class = case_when(
      r_rs_b_Eelgrass_2019 > r_rs_b_Sand_2019 & r_rs_b_Eelgrass_2019 > r_rs_b_Mud_2019 ~ "Eelgrass",
      r_rs_b_Sand_2019 > r_rs_b_Eelgrass_2019 & r_rs_b_Sand_2019 > r_rs_b_Mud_2019 ~ "Sand",
      TRUE ~ "Mud"
    )
  )

  dominant_counts <- table(df_benthic_dominant$dominant_class)
  cat("\nDominant Benthic Class Distribution:\n")
  print(dominant_counts)
  cat(sprintf("  Eelgrass: %.1f%%\n", 100 * dominant_counts["Eelgrass"] / sum(dominant_counts)))
  cat(sprintf("  Sand: %.1f%%\n", 100 * dominant_counts["Sand"] / sum(dominant_counts)))
  cat(sprintf("  Mud: %.1f%%\n", 100 * dominant_counts["Mud"] / sum(dominant_counts)))
}

# Spectral reflectance statistics (only if Approach 2 ran)
if (exists("df_benthic_spectral") && nrow(df_benthic_spectral) > 0) {
  cat("\nMixed Benthic Reflectance Statistics (Approach 2):\n")
  
  # Statistics for each wavelength
  spectral_stats <- df_benthic_spectral %>%
    select(starts_with("Rrs_")) %>%
    summarize(across(everything(), 
                    list(mean = ~mean(., na.rm = TRUE),
                         sd = ~sd(., na.rm = TRUE),
                         min = ~min(., na.rm = TRUE),
                         max = ~max(., na.rm = TRUE))))
  
  cat(sprintf("  Wavelengths: %d bands\n", length(wavelengths_filtered)))
  cat(sprintf("  Mean reflectance range: %.4f - %.4f\n",
              min(as.numeric(spectral_stats[, grep("_mean$", names(spectral_stats))]), na.rm = TRUE),
              max(as.numeric(spectral_stats[, grep("_mean$", names(spectral_stats))]), na.rm = TRUE)))
}

# Algebraic R_B statistics (only if Approach 3 ran)
if (exists("df_algebraic_rb") && nrow(df_algebraic_rb) > 0) {
  cat("\nAlgebraic R_B Statistics (Approach 3):\n")
  cat(sprintf("  Shallow water pixels (0.5-10m): %d\n", 
              ifelse(exists("df_algebraic"), nrow(df_algebraic), 0)))
  cat(sprintf("  Valid R_B pixels after computation: %d\n", 
              ifelse(exists("rb_matrix"), sum(!is.na(rb_matrix[, 1])), 0)))
  cat(sprintf("  Valid RGB pixels after filtering: %d\n", nrow(df_algebraic_rb)))
  
  # R_B range for RGB bands
  cat("\n  R_B range for RGB bands:\n")
  if (exists("wavelengths_range")) {
    for (wl in wavelengths_range[c(4, 3, 2)]) {  # R, G, B
      col_name <- paste0("R_B_", wl)
      if (col_name %in% names(df_algebraic_rb)) {
        rb_vals <- df_algebraic_rb[[col_name]]
        rb_vals <- rb_vals[!is.na(rb_vals)]
        if (length(rb_vals) > 0) {
          cat(sprintf("    %d nm: %.4f - %.4f (mean: %.4f)\n", 
                      wl, min(rb_vals), max(rb_vals), mean(rb_vals)))
        }
      }
    }
  }
}

cat("\n========================================\n")
cat("VISUALIZATION COMPLETE\n")

cat("Generated files:\n")
if (exists("output_fractional") && file.exists(output_fractional)) {
  cat(sprintf("  1. %s\n", basename(output_fractional)))
}
if (exists("output_spectral") && file.exists(output_spectral)) {
  cat(sprintf("  2. %s\n", basename(output_spectral)))
}
if (exists("output_algebraic") && file.exists(output_algebraic)) {
  cat(sprintf("  3. %s\n", basename(output_algebraic)))
}
if (exists("output_comparison") && file.exists(output_comparison)) {
  cat(sprintf("  4. %s\n", basename(output_comparison)))
}
if (exists("output_spectral_nc") && file.exists(output_spectral_nc)) {
  cat(sprintf("  5. %s\n", basename(output_spectral_nc)))
}
if (exists("output_algebraic_nc") && file.exists(output_algebraic_nc)) {
  cat(sprintf("  6. %s\n", basename(output_algebraic_nc)))
}
cat("\nLocation:", OUTPUT_DIR, "\n")
cat("\nApproaches implemented:\n")
cat("  - Fractional RGB: Color mixing of inverted benthic fractions\n")
cat("  - Spectral RGB: Linear mixing of benthic endmember spectra\n")
cat("  - Algebraic R_B: Direct solution from forward model inversion\n")
cat("========================================\n\n")
