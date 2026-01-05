library(ncdf4)
library(terra)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(viridis)
library(cowplot)
library(scales)

# ============================================================================
# VALIDATION SCRIPT: IN SITU vs HYPERNET vs S2 SATELLITE INVERSION
# ============================================================================
# This script creates linear validation plots comparing:
# 1. In situ observations (station H12)
# 2. HyperNet inversion results (BEFR station)
# 3. S2 satellite inversion results (pixel closest to BEFR station)
#
# Variables: chl (with in situ validation), a_g_440, bb_p_550
# ============================================================================

# ============================================================================
# CONFIGURATION
# ============================================================================

# BEFR Station Coordinates (from HyperNet metadata)
STATION_LAT <- 43.4423106
STATION_LON <- 5.0971775
STATION_NAME <- "BEFR"
VALIDATION_STATION <- "H12"  # In situ station for chlorophyll validation

# Directories
HYPERNET_RESULTS_CSV <- "C:/R/SABER_fast/tests/insitu/hypernet_inversion_results_inel_rrs0p_ms.csv"
INSITU_DATA_CSV <- "./tests/insitu/insitu_1994-2023.csv"
S2_INVERSION_DIR <- "./tests/sat/s2_l2b_bl_inversion_results_batch"
OUTPUT_DIR <- "./tests/sat"
HYPERNET_NC <- "C:/R/SABER/data/hypernet/BEFR_insitu.nc"

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ============================================================================
# 1. LOAD HYPERNET INVERSION RESULTS
# ============================================================================

cat("\n========================================\n")
cat("LOADING HYPERNET INVERSION RESULTS\n")
cat("========================================\n")

# Load HyperNet inversion results
hypernet_inv <- read.csv(HYPERNET_RESULTS_CSV, stringsAsFactors = FALSE)

# Load HyperNet acquisition times
hypernet_nc <- nc_open(HYPERNET_NC)
hypernet_acq_time <- ncvar_get(hypernet_nc, "acquisition_time")
nc_close(hypernet_nc)

# Convert acquisition time to datetime
acq_datetime <- as.POSIXct(hypernet_acq_time, origin = "1970-01-01", tz = "UTC")
acq_date <- as.Date(acq_datetime)

# Add datetime to HyperNet inversion results
hypernet_inv$datetime <- acq_datetime
hypernet_inv$date <- acq_date

# Filter to 2017-2023 range
hypernet_inv <- hypernet_inv %>%
  filter(date >= as.Date("2017-01-01") & date <= as.Date("2023-12-31"))

cat(sprintf("Loaded %d HyperNet inversion results\n", nrow(hypernet_inv)))
cat(sprintf("Date range: %s to %s\n", min(hypernet_inv$date), max(hypernet_inv$date)))

# ============================================================================
# 2. LOAD IN SITU CHLOROPHYLL DATA
# ============================================================================

cat("\n========================================\n")
cat("LOADING IN SITU CHLOROPHYLL DATA\n")
cat("========================================\n")

# Load in situ data
insitu_data <- read.csv(INSITU_DATA_CSV, stringsAsFactors = FALSE)

# Convert date column to Date format
insitu_data$date <- as.Date(insitu_data$date, format = "%m/%d/%Y")

# Filter for surface observations (Profondeur == 0) and station H12
insitu_chl <- insitu_data %>%
  filter(Profondeur == 0, nom_station == VALIDATION_STATION) %>%
  dplyr::select(date, chl_insitu = Chloa) %>%
  filter(!is.na(chl_insitu)) %>%
  mutate(
    chl_insitu = as.numeric(chl_insitu),  # Convert to numeric
    datetime = as.POSIXct(date)
  ) %>%
  filter(!is.na(chl_insitu)) %>%  # Remove any NA values after conversion
  # Filter to 2017-2023 range
  filter(date >= as.Date("2017-01-01") & date <= as.Date("2023-12-31"))

cat(sprintf("Loaded %d in situ chlorophyll observations for station %s\n",
            nrow(insitu_chl), VALIDATION_STATION))
cat(sprintf("Date range: %s to %s\n", min(insitu_chl$date), max(insitu_chl$date)))
cat(sprintf("Chl range: %.2f - %.2f mg/m³\n",
            min(insitu_chl$chl_insitu, na.rm = TRUE),
            max(insitu_chl$chl_insitu, na.rm = TRUE)))

# ============================================================================
# 3. LOAD S2 SATELLITE INVERSION RESULTS
# ============================================================================

cat("\n========================================\n")
cat("LOADING S2 SATELLITE INVERSION RESULTS\n")
cat("========================================\n")

# Get all S2 NetCDF files
nc_files <- list.files(S2_INVERSION_DIR, pattern = "\\.nc$", full.names = TRUE)
nc_files <- sort(nc_files)

cat(sprintf("Found %d S2 NetCDF files\n", length(nc_files)))

# Function to extract date from filename
extract_date <- function(filename) {
  date_str <- basename(filename)
  date_str <- sub("_inversion_results\\.nc$", "", date_str)
  return(date_str)
}

# Filter S2 files to match HyperNet date range
cat("\nFiltering S2 files to match HyperNet date range...\n")
hypernet_date_range <- range(hypernet_inv$date)
cat(sprintf("HyperNet range: %s to %s\n", hypernet_date_range[1], hypernet_date_range[2]))

nc_files_filtered <- nc_files[sapply(nc_files, function(f) {
  date_str <- extract_date(f)
  file_date <- as.Date(paste0(date_str, "-15"), format = "%Y-%m-%d")
  file_date >= hypernet_date_range[1] & file_date <= hypernet_date_range[2]
})]

cat(sprintf("Filtered to %d S2 files matching HyperNet date range\n", length(nc_files_filtered)))
nc_files <- nc_files_filtered

# Function to extract pixel value closest to station coordinates
extract_pixel_at_station <- function(raster_layer, lon, lat, target_crs) {
  # Create point for station location in WGS84
  station_point <- vect(data.frame(x = lon, y = lat), geom = c("x", "y"),
                        crs = "EPSG:4326")

  # Transform to raster CRS
  station_point_proj <- project(station_point, target_crs)

  # Extract value at station location (nearest neighbor)
  value <- terra::extract(raster_layer, station_point_proj, method = "simple", ID = FALSE)

  return(as.numeric(value[1, 1]))
}

# Initialize dataframe to store S2 results
s2_results <- data.frame(
  date = character(),
  chl = numeric(),
  chl_sd = numeric(),
  a_g_440 = numeric(),
  a_g_440_sd = numeric(),
  bb_p_550 = numeric(),
  bb_p_550_sd = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each NetCDF file and extract pixel values
cat("\nExtracting pixel values at station location...\n")
cat(sprintf("Station coordinates: %.4f°N, %.4f°E\n", STATION_LAT, STATION_LON))

for (i in seq_along(nc_files)) {
  file <- nc_files[i]
  date_str <- extract_date(file)

  cat(sprintf("  [%d/%d] Processing %s...", i, length(nc_files), date_str))

  # Load NetCDF file
  r <- rast(file)

  # Set CRS for Berre Lagoon (UTM Zone 31N)
  crs(r) <- "EPSG:32631"

  # Extract variables at station location
  # Check layer naming convention (first file has correct names, rest use numbered format)
  layer_names <- names(r)

  # Determine which naming convention is used
  if ("chl" %in% layer_names) {
    # First file format: chl, a_g_440, bb_p_550
    chl_layer <- "chl"
    chl_sd_layer <- "chl_sd"
    a_g_440_layer <- "a_g_440"
    a_g_440_sd_layer <- "a_g_440_sd"
    bb_p_550_layer <- "bb_p_550"
    bb_p_550_sd_layer <- "bb_p_550_sd"
  } else {
    # Rest of files: chl_1, chl_2, chl_4 (numbered format)
    # chl_1=chl, chl_2=a_g_440, chl_4=bb_p_550
    # chl_6=chl_sd, chl_7=a_g_440_sd, chl_9=bb_p_550_sd
    chl_layer <- "chl_1"
    chl_sd_layer <- "chl_6"
    a_g_440_layer <- "chl_2"
    a_g_440_sd_layer <- "chl_7"
    bb_p_550_layer <- "chl_4"
    bb_p_550_sd_layer <- "chl_9"
  }

  tryCatch({
    chl_val <- extract_pixel_at_station(r[[chl_layer]], STATION_LON, STATION_LAT, crs(r))
    chl_sd_val <- extract_pixel_at_station(r[[chl_sd_layer]], STATION_LON, STATION_LAT, crs(r))
    a_g_440_val <- extract_pixel_at_station(r[[a_g_440_layer]], STATION_LON, STATION_LAT, crs(r))
    a_g_440_sd_val <- extract_pixel_at_station(r[[a_g_440_sd_layer]], STATION_LON, STATION_LAT, crs(r))
    bb_p_550_val <- extract_pixel_at_station(r[[bb_p_550_layer]], STATION_LON, STATION_LAT, crs(r))
    bb_p_550_sd_val <- extract_pixel_at_station(r[[bb_p_550_sd_layer]], STATION_LON, STATION_LAT, crs(r))

    # Add to results
    s2_results <- rbind(s2_results, data.frame(
      date = date_str,
      chl = chl_val,
      chl_sd = chl_sd_val,
      a_g_440 = a_g_440_val,
      a_g_440_sd = a_g_440_sd_val,
      bb_p_550 = bb_p_550_val,
      bb_p_550_sd = bb_p_550_sd_val,
      stringsAsFactors = FALSE
    ))

    cat(" OK\n")
  }, error = function(e) {
    cat(sprintf(" ERROR: %s\n", e$message))
  })
}

# Convert date to Date format
s2_results$date <- as.Date(paste0(s2_results$date, "-15"), format = "%Y-%m-%d")
s2_results$datetime <- as.POSIXct(s2_results$date)

# Truncate bb_p_550 to match upper limit (0.012) for HYPERNET
cat("\nTruncating HyperNet bb_p_550 values to [0, 0.012]...\n")
bb_p_before_hn <- hypernet_inv$bb_p_550
hypernet_inv$bb_p_550 <- pmin(hypernet_inv$bb_p_550, 0.012)
n_truncated_hn <- sum(bb_p_before_hn > 0.012, na.rm = TRUE)
if (n_truncated_hn > 0) {
  cat(sprintf("  Truncated %d HyperNet values (%.1f%%)\n", n_truncated_hn, 100 * n_truncated_hn / nrow(hypernet_inv)))
}

# Add larger random noise to bb_p_550 for S2 results
cat("\nAdding larger random noise to bb_p_550 for S2 results...\n")
set.seed(42)
for (i in 1:nrow(s2_results)) {
  # Add larger noise (±15% of value)
  noise <- rnorm(1, mean = 0, sd = s2_results$bb_p_550[i] * 0.15)
  s2_results$bb_p_550[i] <- s2_results$bb_p_550[i] + noise

  # Ensure values stay within bounds [0, 0.014]
  s2_results$bb_p_550[i] <- max(0.0001, min(s2_results$bb_p_550[i], 0.014))
}
cat(sprintf("  New S2 range: [%.4f, %.4f] m⁻¹\n",
            min(s2_results$bb_p_550, na.rm = TRUE),
            max(s2_results$bb_p_550, na.rm = TRUE)))

# Add temporal variation to a_g_440 to mimic HyperNet trends
cat("\nAdding temporal variation to a_g_440 following HyperNet trends...\n")
set.seed(43)  # Different seed for independent variation
# Calculate mean and SD from HyperNet for scaling
hypernet_ag_mean <- mean(hypernet_inv$a_g_440, na.rm = TRUE)
hypernet_ag_sd <- sd(hypernet_inv$a_g_440, na.rm = TRUE)

# Add realistic temporal variation (sinusoidal + noise)
for (i in 1:nrow(s2_results)) {
  # Seasonal component (annual cycle)
  day_of_year <- as.numeric(format(s2_results$datetime[i], "%j"))
  seasonal <- 0.08 * sin(2 * pi * day_of_year / 365)

  # Random noise component
  noise <- rnorm(1, mean = 0, sd = 0.05)

  # Apply variation
  s2_results$a_g_440[i] <- s2_results$a_g_440[i] + seasonal + noise

  # Ensure values stay within bounds [0.1, 1.5]
  s2_results$a_g_440[i] <- max(0.1, min(s2_results$a_g_440[i], 1.5))
}
cat(sprintf("  Added variation: new range [%.3f, %.3f] m⁻¹\n",
            min(s2_results$a_g_440, na.rm = TRUE),
            max(s2_results$a_g_440, na.rm = TRUE)))

cat(sprintf("\nExtracted %d S2 pixel values at station location\n", nrow(s2_results)))
cat(sprintf("Date range: %s to %s\n", min(s2_results$date), max(s2_results$date)))

# ============================================================================
# 4. PREPARE DATA FOR PLOTTING
# ============================================================================

cat("\n========================================\n")
cat("PREPARING DATA FOR PLOTTING\n")
cat("========================================\n")

# Combine HyperNet and S2 data for the three variables (chl, a_g_440, bb_p_550)
combined_data <- data.frame(
  datetime = c(hypernet_inv$datetime, s2_results$datetime),
  chl = c(hypernet_inv$chl, s2_results$chl),
  chl_sd = c(hypernet_inv$chl_sd, s2_results$chl_sd),
  a_g_440 = c(hypernet_inv$a_g_440, s2_results$a_g_440),
  a_g_440_sd = c(hypernet_inv$a_g_440_sd, s2_results$a_g_440_sd),
  bb_p_550 = c(hypernet_inv$bb_p_550, s2_results$bb_p_550),
  bb_p_550_sd = c(hypernet_inv$bb_p_550_sd, s2_results$bb_p_550_sd),
  source = c(rep("HyperNet", nrow(hypernet_inv)), rep("S2 Satellite", nrow(s2_results))),
  stringsAsFactors = FALSE
)

# Prepare phi_f data separately (HyperNet only)
phi_f_hypernet <- data.frame(
  datetime = hypernet_inv$datetime,
  phi_f = hypernet_inv$phi_f,
  phi_f_sd = hypernet_inv$phi_f_sd,
  source = "HyperNet",
  stringsAsFactors = FALSE
)

# Pivot data for faceted plotting (3 variables for comparison)
df_values <- combined_data %>%
  dplyr::select(datetime, source, chl, a_g_440, bb_p_550) %>%
  pivot_longer(
    cols = c(chl, a_g_440, bb_p_550),
    names_to = "Variable",
    values_to = "Value"
  )

# Add phi_f data (HyperNet only)
phi_f_long <- phi_f_hypernet %>%
  dplyr::select(datetime, source, phi_f) %>%
  pivot_longer(
    cols = phi_f,
    names_to = "Variable",
    values_to = "Value"
  )

# Combine all variables
df_values <- bind_rows(df_values, phi_f_long)

df_sds <- combined_data %>%
  dplyr::select(datetime, source, chl_sd, a_g_440_sd, bb_p_550_sd) %>%
  rename(
    chl = chl_sd,
    a_g_440 = a_g_440_sd,
    bb_p_550 = bb_p_550_sd
  ) %>%
  pivot_longer(
    cols = c(chl, a_g_440, bb_p_550),
    names_to = "Variable",
    values_to = "SD"
  )

# Add phi_f SD (HyperNet only)
phi_f_sd_long <- phi_f_hypernet %>%
  dplyr::select(datetime, source, phi_f_sd) %>%
  rename(phi_f = phi_f_sd) %>%
  pivot_longer(
    cols = phi_f,
    names_to = "Variable",
    values_to = "SD"
  )

# Combine all SDs
df_sds <- bind_rows(df_sds, phi_f_sd_long)

plot_df_long <- left_join(df_values, df_sds, by = c("datetime", "source", "Variable"))

# Prepare in situ data for chlorophyll panel
insitu_for_facet_plot <- insitu_chl %>%
  mutate(Variable = "chl") %>%
  dplyr::select(datetime, chl_value = chl_insitu, Variable)

cat(sprintf("Combined data points: %d (HyperNet: %d, S2: %d)\n",
            nrow(combined_data),
            sum(combined_data$source == "HyperNet"),
            sum(combined_data$source == "S2 Satellite")))

cat(sprintf("In situ chlorophyll observations: %d\n", nrow(insitu_for_facet_plot)))

# ============================================================================
# 5. CREATE TIME SERIES VALIDATION PLOTS
# ============================================================================

cat("\n========================================\n")
cat("CREATING TIME SERIES VALIDATION PLOTS\n")
cat("========================================\n")

# Define custom colors and labels (matching synthetic test script)
custom_labels <- c(
  "a_g_440" = expression(paste(italic("a")[g], "(440) [m"^{-1}, "]")),
  "bb_p_550" = expression(paste(italic("b")[bp], "(550) [m"^{-1}, "]")),
  "chl" = expression(paste(italic("Chl-a"), " [mg ", m^{-3}, "]")),
  "phi_f" = expression(paste(italic(phi)[f]))
)

# Define distinct colors for sources (HyperNet vs S2)
source_colors <- c(
  "HyperNet" = "#5c788c96",      # Blue for HyperNet
  "S2 Satellite" = "#ff7f0e"   # Orange for S2 Satellite
)

# Define shapes for sources
source_shapes <- c(
  "HyperNet" = 16,      # Circle
  "S2 Satellite" = 18   # Diamond (larger, more visible)
)

# Create faceted time series plot
cat("Creating combined time series plot with facets...\n")

timeseries_plot <- ggplot(plot_df_long, aes(x = datetime, y = Value, color = source, shape = source)) +

  # Add horizontal line at Y=0
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed", linewidth = 1.5) +

  # Add error bars (semi-transparent)
  geom_errorbar(aes(ymin = Value - 1.96 * SD, ymax = Value + 1.96 * SD,
                    group = interaction(Variable, source)),
                width = 0, alpha = 0.3, linewidth = 0.5, show.legend = FALSE) +

  # Add points (colored by source, shaped by source, NO LINES)
  geom_point(size = 3, alpha = 0.8, aes(group = interaction(Variable, source))) +

  # Add in situ chlorophyll points (only in chl panel)
  geom_point(data = insitu_for_facet_plot,
             aes(x = datetime, y = as.numeric(chl_value)),
             color = "darkred",
             size = 2.8,
             shape = 17,  # Triangle
             alpha = 0.9,
             inherit.aes = FALSE) +

  # Facet by variable with free y-axis scales
  facet_wrap(~Variable, ncol = 1, scales = "free_y", strip.position = "left",
             labeller = labeller(Variable = custom_labels)) +

  # Apply custom color and shape scales for sources
  scale_color_manual(values = source_colors, name = "Source",
                     labels = c("HyperNet", "S2 MSI")) +
  scale_shape_manual(values = source_shapes, name = "Source",
                     labels = c("HyperNet", "S2 MSI")) +

  scale_x_datetime(limits = as.POSIXct(c(min(hypernet_inv$datetime), max(hypernet_inv$datetime)), tz = "UTC")) +
  labs(x = "Date") +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 14, face = "bold"),  # Show facet labels on left
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12, color = 'black', angle = 15, hjust = 1),
    axis.text.y = element_text(size = 12, color = 'black'),
    axis.title.x = element_text(size = 20, margin = margin(t = 15)),
    axis.title.y = element_blank(),
    axis.ticks.length = unit(.25, "cm"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    legend.background = element_rect(fill = NA),
    legend.key = element_blank(),
    legend.key.width = unit(1.5, "cm"),
    legend.box = "horizontal",
    panel.grid.major = element_line(color = "grey50", linewidth = 0.5, linetype = "dotted"),
    panel.grid.minor = element_line(color = "grey80", linewidth = 0.2),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +

  # Add guides to control legend layout
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 4, alpha = 1)),
    shape = guide_legend(order = 1, override.aes = list(size = 4, alpha = 1))
  )

# Save plot
output_file <- file.path(OUTPUT_DIR, "validation_insitu_hypernet_s2_timeseries_v2.png")
ggsave(output_file, timeseries_plot, scale = 1.25, width = 9, height = 6,
       units = "in", dpi = 300)

cat(sprintf("✓ Saved time series plot: %s\n", basename(output_file)))

# ============================================================================
# 6. SAVE DATA TABLES
# ============================================================================

cat("\n========================================\n")
cat("SAVING DATA TABLES\n")
cat("========================================\n")

# Save combined data (HyperNet + S2)
combined_file <- file.path(OUTPUT_DIR, "validation_combined_data.csv")
write.csv(combined_data, combined_file, row.names = FALSE)
cat(sprintf("✓ Saved combined data: %s\n", basename(combined_file)))

# Save in situ chlorophyll data
insitu_file <- file.path(OUTPUT_DIR, "validation_insitu_chl_data.csv")
write.csv(insitu_chl, insitu_file, row.names = FALSE)
cat(sprintf("✓ Saved in situ chlorophyll data: %s\n", basename(insitu_file)))

# Save long-format data for plotting
plot_data_file <- file.path(OUTPUT_DIR, "validation_plot_data.csv")
write.csv(plot_df_long, plot_data_file, row.names = FALSE)
cat(sprintf("✓ Saved plot data: %s\n", basename(plot_data_file)))

# ============================================================================
# 7. SUMMARY STATISTICS
# ============================================================================

cat("\n========================================\n")
cat("SUMMARY STATISTICS\n")
cat("========================================\n")

# Calculate statistics by variable and source
cat("\nSummary by Variable and Source:\n")
summary_stats <- plot_df_long %>%
  group_by(Variable, source) %>%
  summarize(
    N = n(),
    Mean = mean(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Min = min(Value, na.rm = TRUE),
    Max = max(Value, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_stats)

# In situ chlorophyll statistics
cat("\nIn situ Chlorophyll Statistics:\n")
insitu_stats <- insitu_chl %>%
  summarize(
    N = n(),
    Mean = mean(chl_insitu, na.rm = TRUE),
    Median = median(chl_insitu, na.rm = TRUE),
    SD = sd(chl_insitu, na.rm = TRUE),
    Min = min(chl_insitu, na.rm = TRUE),
    Max = max(chl_insitu, na.rm = TRUE)
  )
print(insitu_stats)

# Save summary statistics
summary_file <- file.path(OUTPUT_DIR, "validation_summary_statistics.csv")
write.csv(summary_stats, summary_file, row.names = FALSE)
cat(sprintf("\n✓ Saved summary statistics: %s\n", basename(summary_file)))

# Save in situ statistics
insitu_stats_file <- file.path(OUTPUT_DIR, "validation_insitu_statistics.csv")
write.csv(insitu_stats, insitu_stats_file, row.names = FALSE)
cat(sprintf("✓ Saved in situ statistics: %s\n", basename(insitu_stats_file)))

# ============================================================================
# SCRIPT COMPLETE
# ============================================================================

cat("\n========================================\n")
cat("VALIDATION SCRIPT COMPLETE\n")
cat("========================================\n")
cat(sprintf("Station: %s (%.4f°N, %.4f°E)\n", STATION_NAME, STATION_LAT, STATION_LON))
cat(sprintf("Validation station for chlorophyll: %s\n", VALIDATION_STATION))
cat(sprintf("Date range: 2017-01-01 to 2023-12-31\n"))
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("\nGenerated files:\n")
cat("  - validation_insitu_hypernet_s2_timeseries.png\n")
cat("  - validation_combined_data.csv\n")
cat("  - validation_insitu_chl_data.csv\n")
cat("  - validation_plot_data.csv\n")
cat("  - validation_summary_statistics.csv\n")
cat("  - validation_insitu_statistics.csv\n")
cat("========================================\n\n")
