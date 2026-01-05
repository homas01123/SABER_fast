library(terra)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tidyverse)
library(scales)
library(SABER)
library(parallel)
library(doParallel)

library(raster)
library(ggplot2)
library(rasterVis)
library(ggspatial)
library(sf)
library(cowplot)
library(progressr)
library(xml2)
library(suncalc)
library(lubridate)
#gc()

# S2 SABER INVERSION FOR OPTICALLY SHALLOW WATERS ----

# ============================================================================
# FUNCTION: EXTRACT SOLAR AND VIEWING ANGLES FROM L1C XML METADATA
# ============================================================================
extract_angles_from_l1c <- function(l2w_file, l1c_base_dir = "./tests/sat/s2_l1c_mp") {
  # Parse L2W filename to extract date/time info
  # Format: S2X_MSI_YYYY_MM_DD_HH_MM_SS_TXXXXX_L2W.nc
  l2w_basename <- basename(l2w_file)

  # Extract components
  parts <- strsplit(l2w_basename, "_")[[1]]
  satellite <- parts[1]  # S2A or S2B
  year <- parts[3]
  month <- parts[4]
  day <- parts[5]
  hour <- parts[6]
  minute <- parts[7]
  second <- parts[8]
  tile <- parts[9]

  # Construct datetime for suncalc fallback
  datetime_str <- sprintf("%s-%s-%s %s:%s:%s", year, month, day, hour, minute, second)
  image_datetime <- as.POSIXct(datetime_str, tz = "UTC")

  # Construct L1C search pattern
  # L1C format: S2X_MSIL1C_YYYYMMDDTHHMMSS_NXXXX_RXXX_TXXXXX_YYYYMMDDTXXXXXX.SAFE
  date_str <- sprintf("%s%s%s", year, month, day)
  time_str <- sprintf("%s%s%s", hour, minute, second)

  cat(sprintf("\nSearching for L1C file matching: %s %s %s %s\n",
              satellite, date_str, time_str, tile))

  # Search for matching L1C .SAFE directory
  l1c_pattern <- sprintf("%s_MSIL1C_%sT%s_*_%s_*.SAFE",
                         satellite, date_str, time_str, tile)

  l1c_dirs <- list.dirs(l1c_base_dir, full.names = TRUE, recursive = FALSE)
  l1c_match <- l1c_dirs[grep(l1c_pattern, basename(l1c_dirs))]

  if (length(l1c_match) == 0) {
    cat(sprintf("WARNING: No matching L1C directory found for pattern: %s\n", l1c_pattern))
    cat("Using suncalc package to estimate sun angles from image location and time\n")

    # Calculate center from study area coordinates (will be passed as parameter)
    # This is a placeholder - actual calculation done in main script
    # For now, return NULL and handle in main script
    return(list(
      mean_sza = NA,
      mean_saa = NA,
      mean_vza = 0.0,
      mean_vaa = 0.0,
      sza_raster = NULL,
      saa_raster = NULL,
      vza_raster = NULL,
      vaa_raster = NULL,
      image_datetime = image_datetime,
      needs_suncalc = TRUE
    ))
  }

  l1c_dir <- l1c_match[1]
  cat(sprintf("Found L1C: %s\n", basename(l1c_dir)))

  # Find MTD_TL.xml inside GRANULE subdirectory
  granule_dir <- file.path(l1c_dir, "GRANULE")
  if (!dir.exists(granule_dir)) {
    cat("WARNING: GRANULE directory not found in L1C .SAFE\n")
    cat("Falling back to suncalc estimation\n")

    return(list(
      mean_sza = NA,
      mean_saa = NA,
      mean_vza = 0.0,
      mean_vaa = 0.0,
      sza_raster = NULL,
      saa_raster = NULL,
      vza_raster = NULL,
      vaa_raster = NULL,
      image_datetime = image_datetime,
      needs_suncalc = TRUE
    ))
  }

  granule_subdirs <- list.dirs(granule_dir, full.names = TRUE, recursive = FALSE)
  if (length(granule_subdirs) == 0) {
    cat("WARNING: No subdirectories in GRANULE\n")

    return(list(
      mean_sza = NA,
      mean_saa = NA,
      mean_vza = 0.0,
      mean_vaa = 0.0,
      sza_raster = NULL,
      saa_raster = NULL,
      vza_raster = NULL,
      vaa_raster = NULL,
      image_datetime = image_datetime,
      needs_suncalc = TRUE
    ))
  }

  xml_file <- file.path(granule_subdirs[1], "MTD_TL.xml")
  if (!file.exists(xml_file)) {
    cat(sprintf("WARNING: MTD_TL.xml not found at: %s\n", xml_file))

    return(list(
      mean_sza = NA,
      mean_saa = NA,
      mean_vza = 0.0,
      mean_vaa = 0.0,
      sza_raster = NULL,
      saa_raster = NULL,
      vza_raster = NULL,
      vaa_raster = NULL,
      image_datetime = image_datetime,
      needs_suncalc = TRUE
    ))
  }

  cat(sprintf("Reading XML: %s\n", basename(xml_file)))

  # Read XML
  xml_data <- read_xml(xml_file)

  # Extract mean sun angles
  mean_sza <- xml_find_first(xml_data, ".//Mean_Sun_Angle/ZENITH_ANGLE") %>%
    xml_text() %>% as.numeric()
  mean_saa <- xml_find_first(xml_data, ".//Mean_Sun_Angle/AZIMUTH_ANGLE") %>%
    xml_text() %>% as.numeric()

  cat(sprintf("Mean Solar Zenith Angle: %.2f deg\n", mean_sza))
  cat(sprintf("Mean Solar Azimuth Angle: %.2f deg\n", mean_saa))

  # Extract mean viewing angles (across all bands)
  mean_vza_nodes <- xml_find_all(xml_data, ".//Mean_Viewing_Incidence_Angle")

  vza_values <- numeric(0)
  vaa_values <- numeric(0)

  for (node in mean_vza_nodes) {
    vza <- xml_find_first(node, ".//ZENITH_ANGLE") %>% xml_text() %>% as.numeric()
    vaa <- xml_find_first(node, ".//AZIMUTH_ANGLE") %>% xml_text() %>% as.numeric()
    vza_values <- c(vza_values, vza)
    vaa_values <- c(vaa_values, vaa)
  }

  mean_vza <- mean(vza_values, na.rm = TRUE)
  mean_vaa <- mean(vaa_values, na.rm = TRUE)

  cat(sprintf("Mean Viewing Zenith Angle: %.2f deg\n", mean_vza))
  cat(sprintf("Mean Viewing Azimuth Angle: %.2f deg\n", mean_vaa))

  # Extract sun angle grids
  sun_zenith_node <- xml_find_first(xml_data, ".//Sun_Angles_Grid/Zenith")
  sza_col_step <- xml_find_first(sun_zenith_node, ".//COL_STEP") %>%
    xml_text() %>% as.numeric()
  sza_row_step <- xml_find_first(sun_zenith_node, ".//ROW_STEP") %>%
    xml_text() %>% as.numeric()
  sza_values_nodes <- xml_find_all(sun_zenith_node, ".//VALUES")
  sza_values_text <- xml_text(sza_values_nodes)

  # Parse grid values
  sza_grid <- lapply(sza_values_text, function(row) {
    vals <- strsplit(row, " +")[[1]]
    vals <- vals[vals != ""]
    as.numeric(vals)
  })
  sza_matrix <- do.call(rbind, sza_grid)

  # Get tile geocoding info
  tile_ulx <- xml_find_first(xml_data, ".//Geoposition[@resolution='10']/ULX") %>%
    xml_text() %>% as.numeric()
  tile_uly <- xml_find_first(xml_data, ".//Geoposition[@resolution='10']/ULY") %>%
    xml_text() %>% as.numeric()
  tile_crs <- xml_find_first(xml_data, ".//HORIZONTAL_CS_CODE") %>% xml_text()

  # Calculate extent for the angle grid
  grid_ncols <- ncol(sza_matrix)
  grid_nrows <- nrow(sza_matrix)

  angle_xmin <- tile_ulx
  angle_xmax <- tile_ulx + (grid_ncols) * sza_col_step
  angle_ymax <- tile_uly
  angle_ymin <- tile_uly - (grid_nrows) * sza_row_step

  # Create rasters for sun angles
  r_sza <- rast(sza_matrix,
                extent = c(angle_xmin, angle_xmax, angle_ymin, angle_ymax),
                crs = tile_crs)
  names(r_sza) <- "sza"

  # Sun Azimuth Grid
  sun_azimuth_node <- xml_find_first(xml_data, ".//Sun_Angles_Grid/Azimuth")
  saa_values_nodes <- xml_find_all(sun_azimuth_node, ".//VALUES")
  saa_values_text <- xml_text(saa_values_nodes)

  saa_grid <- lapply(saa_values_text, function(row) {
    vals <- strsplit(row, " +")[[1]]
    vals <- vals[vals != ""]
    as.numeric(vals)
  })
  saa_matrix <- do.call(rbind, saa_grid)

  r_saa <- rast(saa_matrix,
                extent = c(angle_xmin, angle_xmax, angle_ymin, angle_ymax),
                crs = tile_crs)
  names(r_saa) <- "saa"

  cat(sprintf("SZA grid: %d x %d, range: %.2f - %.2f deg\n",
              nrow(sza_matrix), ncol(sza_matrix),
              min(sza_matrix, na.rm=TRUE), max(sza_matrix, na.rm=TRUE)))

  # Extract viewing angle grids (Band 0 as representative)
  # Viewing angles are stored per band and detector - use Band 0 (B1) as representative
  cat("\nExtracting viewing angle grids...\n")

  # Try to get viewing angles for first band
  viewing_angle_nodes <- xml_find_all(xml_data, ".//Viewing_Incidence_Angles_Grids[@bandId='0']")

  if (length(viewing_angle_nodes) > 0) {
    # Get zenith angles for all detectors of band 0
    vza_all_detectors <- list()

    for (va_node in viewing_angle_nodes) {
      vza_node <- xml_find_first(va_node, ".//Zenith")
      vza_values_nodes <- xml_find_all(vza_node, ".//VALUES")
      vza_values_text <- xml_text(vza_values_nodes)

      vza_grid <- lapply(vza_values_text, function(row) {
        vals <- strsplit(row, " +")[[1]]
        vals <- vals[vals != ""]
        as.numeric(vals)
      })
      vza_matrix_detector <- do.call(rbind, vza_grid)
      vza_all_detectors[[length(vza_all_detectors) + 1]] <- vza_matrix_detector
    }

    # Average across detectors
    vza_matrix <- Reduce("+", vza_all_detectors) / length(vza_all_detectors)

    # Create VZA raster
    r_vza <- rast(vza_matrix,
                  extent = c(angle_xmin, angle_xmax, angle_ymin, angle_ymax),
                  crs = tile_crs)
    names(r_vza) <- "vza"

    cat(sprintf("VZA grid: %d x %d, range: %.2f - %.2f deg\n",
                nrow(vza_matrix), ncol(vza_matrix),
                min(vza_matrix, na.rm=TRUE), max(vza_matrix, na.rm=TRUE)))

    # Get azimuth angles
    vaa_all_detectors <- list()

    for (va_node in viewing_angle_nodes) {
      vaa_node <- xml_find_first(va_node, ".//Azimuth")
      vaa_values_nodes <- xml_find_all(vaa_node, ".//VALUES")
      vaa_values_text <- xml_text(vaa_values_nodes)

      vaa_grid <- lapply(vaa_values_text, function(row) {
        vals <- strsplit(row, " +")[[1]]
        vals <- vals[vals != ""]
        as.numeric(vals)
      })
      vaa_matrix_detector <- do.call(rbind, vaa_grid)
      vaa_all_detectors[[length(vaa_all_detectors) + 1]] <- vaa_matrix_detector
    }

    vaa_matrix <- Reduce("+", vaa_all_detectors) / length(vaa_all_detectors)

    r_vaa <- rast(vaa_matrix,
                  extent = c(angle_xmin, angle_xmax, angle_ymin, angle_ymax),
                  crs = tile_crs)
    names(r_vaa) <- "vaa"

  } else {
    cat("WARNING: Could not extract viewing angle grids, using mean values\n")
    r_vza <- NULL
    r_vaa <- NULL
  }

  return(list(
    mean_sza = mean_sza,
    mean_saa = mean_saa,
    mean_vza = mean_vza,
    mean_vaa = mean_vaa,
    sza_raster = r_sza,
    saa_raster = r_saa,
    vza_raster = r_vza,
    vaa_raster = r_vaa
  ))
}

# This script processes Sentinel-2 L2W data for shallow water conditions
# with benthic reflectance contribution (Eelgrass, Sand, Mud)


# Function :: Convert raster to dataframe for ggplot ----
raster_to_df <- function(rast) {
  df <- as.data.frame(rast, xy = TRUE)
  names(df)[3] <- "value"
  return(df)
}


## DEFINE STUDY AREA BOUNDARY ----

# Create polygon from coordinates
study_area_coords <- matrix(c(
  -68.107338, 49.191128,
  -68.529625, 49.193147,
  -68.532028, 49.01423,
  -68.106308, 49.015356,
  -68.107338, 49.191128
), ncol = 2, byrow = TRUE)

study_area_poly <- st_polygon(list(study_area_coords))
study_area_sf <- st_sfc(study_area_poly, crs = 4326)  # WGS84

## DEPTH THRESHOLDS FOR WATER TYPE CLASSIFICATION ----

# Minimum depth for inversion (pixels below this are skipped)
DEPTH_MIN <- 0.25  # meters

# Maximum depth for shallow water classification (optically shallow with benthic influence)
# Pixels > DEPTH_SHALLOW_MAX are classified as deep water (no benthic influence)
DEPTH_SHALLOW_MAX <- 3.0  # meters

cat("\nStudy area boundary defined\n")
cat(sprintf("Bounding box: [%.6f, %.6f] to [%.6f, %.6f]\n",
            min(study_area_coords[,1]), min(study_area_coords[,2]),
            max(study_area_coords[,1]), max(study_area_coords[,2])))


## LOAD S2 L2W DATA ----

s2_l2w_file <- "./tests/sat/s2_l2a_mp/S2B_MSI_2019_07_07_15_49_31_T19UEQ_L2W.nc"

s2_img_name <- basename(s2_l2w_file)

#unlist to remove the .nc part
s2_img_name <- unlist(strsplit(s2_img_name, "\\."))[1]

cat("\nLoading S2 L2W data from:", s2_l2w_file, "\n")

# Load the L2W NetCDF file
s2_l2w <- rast(s2_l2w_file)

cat("Available bands:\n")
print(names(s2_l2w))

# Extract Rrs bands (look for Rrs or rhos_* patterns in L2W)
# ACOLITE L2W typically has Rrs_* for water-leaving reflectance
rrs_bands <- grep("^Rrs_", names(s2_l2w), value = TRUE)

if (length(rrs_bands) == 0) {
  cat("No Rrs_* bands found. Checking for rhos_* bands...\n")
  rrs_bands <- grep("^rhos_", names(s2_l2w), value = TRUE)
}

cat("\nUsing reflectance bands:\n")
print(rrs_bands)

# Select the Rrs spectral bands
spectral_raster_full <- s2_l2w[[rrs_bands[1:6]]]

cat(sprintf("\nFull image dimensions: %d rows x %d cols x %d bands\n",
            nrow(spectral_raster_full), ncol(spectral_raster_full), nlyr(spectral_raster_full)))


## CLIP TO STUDY AREA ----

cat("\nClipping raster to study area...\n")

# Convert study area to raster CRS
study_area_vect <- vect(study_area_sf)
study_area_transformed <- project(study_area_vect, crs(spectral_raster_full))

# Crop and mask the raster
spectral_raster <- crop(spectral_raster_full, study_area_transformed)
spectral_raster <- mask(spectral_raster, study_area_transformed)

cat(sprintf("Clipped image dimensions: %d rows x %d cols x %d bands\n",
            nrow(spectral_raster), ncol(spectral_raster), nlyr(spectral_raster)))

# Calculate size reduction
original_pixels <- ncell(spectral_raster_full)
clipped_pixels <- sum(values(!is.na(spectral_raster[[1]])))
cat(sprintf("Pixels reduced from %d to %d (%.1f%% reduction)\n",
            original_pixels, clipped_pixels,
            (1 - clipped_pixels/original_pixels) * 100))


## LOAD AND RESAMPLE BATHYMETRY DATA ----

cat("\nLoading bathymetry data...\n")
bathy_file <- "c:/R/SABER_fast/tests/sat/bathymetry_mp/Asli_bhara_hua_NONNA_bath.tiff"

if (file.exists(bathy_file)) {
  # Load bathymetry
  bathy_full <- rast(bathy_file)

  cat(sprintf("Original bathymetry: %d rows x %d cols\n",
              nrow(bathy_full), ncol(bathy_full)))
  cat(sprintf("Bathymetry extent: [%.4f, %.4f, %.4f, %.4f]\n",
              ext(bathy_full)[1], ext(bathy_full)[2],
              ext(bathy_full)[3], ext(bathy_full)[4]))

  # Resample bathymetry to match spectral raster grid
  cat("Resampling bathymetry to S2 grid (this may take a moment)...\n")

  # Reproject bathymetry to match spectral raster CRS if needed
  if (crs(bathy_full) != crs(spectral_raster)) {
    cat("Reprojecting bathymetry to match S2 CRS...\n")
    bathy_full <- project(bathy_full, crs(spectral_raster))
  }

  bathy_resampled <- resample(bathy_full, spectral_raster, method = "bilinear")

  # Mask to study area
  bathy_resampled <- mask(bathy_resampled, study_area_transformed)

  cat(sprintf("Resampled bathymetry: %d rows x %d cols\n",
              nrow(bathy_resampled), ncol(bathy_resampled)))

  # Get bathymetry statistics
  bathy_vals <- values(bathy_resampled)
  bathy_valid <- bathy_vals[!is.na(bathy_vals)]

  if (length(bathy_valid) > 0) {
    cat(sprintf("Bathymetry depth range: %.2f - %.2f m\n",
                min(bathy_valid), max(bathy_valid)))
    cat(sprintf("Mean depth: %.2f m, Median: %.2f m\n",
                mean(bathy_valid), median(bathy_valid)))

    # Flag for using bathymetry constraints
    use_bathy_constraint <- TRUE
  } else {
    cat("Warning: No valid bathymetry data in study area\n")
    use_bathy_constraint <- FALSE
    bathy_resampled <- NULL
  }
} else {
  cat("Warning: Bathymetry file not found, proceeding without depth constraint\n")
  use_bathy_constraint <- FALSE
  bathy_resampled <- NULL
}

# Plot to visualize
cat("\nPlotting spectral data and bathymetry...\n")
plot(spectral_raster[[1:4]])

# Plot bathymetry if available
if (use_bathy_constraint && !is.null(bathy_resampled)) {
  cat("Plotting bathymetry...\n")
  plot(bathy_resampled, main = "Bathymetry (m)")
}

dims <- dim(spectral_raster)
nrows <- dims[1]
ncols <- dims[2]
nbands <- dims[3]

cat(sprintf("\nImage dimensions: %d rows x %d cols x %d bands\n", nrows, ncols, nbands))


## PREPARE DATA FOR INVERSION ----

cat("\nPreparing data for inversion...\n")

# Convert raster to dataframe
df <- as.data.frame(spectral_raster, xy = TRUE, cells = TRUE, na.rm = TRUE) %>%
  mutate(
    row = rowFromCell(spectral_raster, cell),
    col = colFromCell(spectral_raster, cell),
    ensemble = paste0("id_", row_number())
  )

# Add bathymetry data if available
if (use_bathy_constraint && !is.null(bathy_resampled)) {
  bathy_df <- as.data.frame(bathy_resampled, xy = TRUE, cells = TRUE, na.rm = FALSE)
  names(bathy_df)[4] <- "bathy_depth"

  df <- df %>%
    left_join(bathy_df %>% dplyr::select(cell, bathy_depth), by = "cell")

  # Check coverage alignment between S2 water pixels and bathymetry
  cat("\n=== BATHYMETRY COVERAGE ANALYSIS ===\n")
  cat(sprintf("Total S2 water pixels: %d\n", nrow(df)))
  cat(sprintf("Pixels with valid bathymetry: %d\n", sum(!is.na(df$bathy_depth))))
  cat(sprintf("Pixels with NA bathymetry: %d\n", sum(is.na(df$bathy_depth))))
  cat(sprintf("Coverage: %.1f%%\n", 100 * sum(!is.na(df$bathy_depth)) / nrow(df)))

  if (sum(is.na(df$bathy_depth)) > 0) {
    cat("\nWARNING: Bathymetry has gaps in S2 water pixel coverage!\n")
    cat("Option 1: Filter to only pixels with valid bathymetry (RECOMMENDED)\n")
    cat("Option 2: Use unconstrained approach for all pixels\n")

    # FILTER: Keep only pixels with valid bathymetry
    cat("\n==> Filtering to pixels with valid bathymetry only\n")
    df_original_count <- nrow(df)
    df <- df %>% filter(!is.na(bathy_depth))
    cat(sprintf("Filtered from %d to %d pixels (%.1f%% retained)\n",
                df_original_count, nrow(df),
                100 * nrow(df) / df_original_count))
  } else {
    cat("\n==> Perfect alignment! All S2 water pixels have bathymetry data\n")
  }

  # Get bathymetry statistics for valid pixels
  bathy_valid <- df$bathy_depth[!is.na(df$bathy_depth)]
  if (length(bathy_valid) > 0) {
    cat(sprintf("\nBathymetry statistics for valid pixels:\n"))
    cat(sprintf("  Range: %.2f - %.2f m\n", min(bathy_valid), max(bathy_valid)))
    cat(sprintf("  Mean: %.2f m, Median: %.2f m\n", mean(bathy_valid), median(bathy_valid)))
    cat(sprintf("  Quantiles: %.2f (25%%), %.2f (75%%)\n",
                quantile(bathy_valid, 0.25), quantile(bathy_valid, 0.75)))

    # Classify pixels by depth for hybrid inversion approach
    too_shallow <- sum(bathy_valid < DEPTH_MIN)
    shallow_water <- sum(bathy_valid >= DEPTH_MIN & bathy_valid <= DEPTH_SHALLOW_MAX)
    deep_water <- sum(bathy_valid > DEPTH_SHALLOW_MAX)

    cat(sprintf("\nHybrid inversion classification:\n"))
    if (too_shallow > 0) cat(sprintf("  < %.2fm: %d pixels (%.1f%%) - SKIPPED (too shallow)\n",
                                     DEPTH_MIN, too_shallow, 100 * too_shallow / length(bathy_valid)))
    cat(sprintf("  %.2f-%.1fm: %d pixels (%.1f%%) - OPTICALLY SHALLOW (with benthic)\n",
                DEPTH_MIN, DEPTH_SHALLOW_MAX, shallow_water, 100 * shallow_water / length(bathy_valid)))
    cat(sprintf("  > %.1fm: %d pixels (%.1f%%) - OPTICALLY DEEP (no benthic)\n",
                DEPTH_SHALLOW_MAX, deep_water, 100 * deep_water / length(bathy_valid)))

    # Filter out pixels < DEPTH_MIN (too shallow for inversion)
    if (too_shallow > 0) {
      cat(sprintf("\n==> Removing %d pixels with depth < %.2fm\n", too_shallow, DEPTH_MIN))
      df_before_filter <- nrow(df)
      df <- df %>% filter(bathy_depth >= DEPTH_MIN)
      cat(sprintf("Filtered from %d to %d pixels\n", df_before_filter, nrow(df)))
    }

    # Add water type classification
    df <- df %>%
      mutate(water_type_class = ifelse(bathy_depth > DEPTH_SHALLOW_MAX, "deep", "shallow"))

    cat(sprintf("\nFinal pixel counts for inversion:\n"))
    cat(sprintf("  Shallow water (%.2f-%.1fm): %d pixels\n",
                DEPTH_MIN, DEPTH_SHALLOW_MAX, sum(df$water_type_class == "shallow")))
    cat(sprintf("  Deep water (>%.1fm): %d pixels\n",
                DEPTH_SHALLOW_MAX, sum(df$water_type_class == "deep")))
  }
  cat("====================================\n\n")
}


## EXTRACT SOLAR AND VIEWING ZENITH ANGLES DYNAMICALLY FROM L1C XML ----

cat("\nExtracting solar and viewing zenith angles from L1C metadata...\n")

# Extract angles dynamically for this specific image
angle_data <- extract_angles_from_l1c(s2_l2w_file)

# Check if suncalc fallback is needed
if (!is.null(angle_data$needs_suncalc) && angle_data$needs_suncalc) {
  cat("Calculating sun position from image extent...\n")

  # Get extent from loaded raster (in projected CRS)
  extent_utm <- ext(spectral_raster_full)

  # Convert extent corners to lat/lon
  corners_utm <- matrix(c(
    extent_utm[1], extent_utm[3],  # xmin, ymin (SW corner)
    extent_utm[2], extent_utm[4]   # xmax, ymax (NE corner)
  ), ncol = 2, byrow = TRUE)

  corners_latlon <- project(corners_utm,
                            from = crs(spectral_raster_full),
                            to = "EPSG:4326")

  # Calculate center
  center_lon <- mean(corners_latlon[, 1])
  center_lat <- mean(corners_latlon[, 2])

  cat(sprintf("Image center: %.4f°E, %.4f°N\n", center_lon, center_lat))

  # Calculate sun position
  sun_pos <- getSunlightPosition(date = angle_data$image_datetime,
                                  lat = center_lat,
                                  lon = center_lon)

  mean_sza <- 90 - (sun_pos$altitude * 180 / pi)
  mean_saa <- sun_pos$azimuth * 180 / pi

  cat(sprintf("Calculated SZA: %.2f deg, SAA: %.2f deg\n", mean_sza, mean_saa))

  # Update angle_data
  angle_data$mean_sza <- mean_sza
  angle_data$mean_saa <- mean_saa
}

mean_sza <- angle_data$mean_sza
mean_saa <- angle_data$mean_saa
mean_vza <- angle_data$mean_vza
mean_vaa <- angle_data$mean_vaa

# Initialize flags for pixel-wise angle usage
use_pixel_sza <- FALSE
use_pixel_vza <- FALSE
sza_raster <- NULL
vza_raster <- NULL

# Threshold for coefficient of variation (CV = sd/mean)
# If CV > 0.05 (5%), use pixel-wise angles
CV_THRESHOLD <- 0.05

if (!is.null(angle_data$sza_raster)) {
  # Resample SZA to L2W grid
  cat("\nProcessing solar zenith angles...\n")
  l2w_raster_template <- spectral_raster_full[[1]]
  sza_raster_full <- resample(angle_data$sza_raster, l2w_raster_template, method = "bilinear")

  # Crop and mask to study area
  sza_raster <- crop(sza_raster_full, study_area_transformed)
  sza_raster <- mask(sza_raster, study_area_transformed)

  # Calculate statistics over study area
  sza_values <- values(sza_raster)
  sza_values <- sza_values[!is.na(sza_values)]

  if (length(sza_values) > 0) {
    mean_sza_area <- mean(sza_values)
    sd_sza_area <- sd(sza_values)
    cv_sza <- sd_sza_area / mean_sza_area

    cat(sprintf("SZA statistics over study area:\n"))
    cat(sprintf("  Range: %.2f - %.2f deg\n", min(sza_values), max(sza_values)))
    cat(sprintf("  Mean: %.2f deg\n", mean_sza_area))
    cat(sprintf("  SD: %.2f deg\n", sd_sza_area))
    cat(sprintf("  CV: %.4f (%.2f%%)\n", cv_sza, cv_sza * 100))

    # Decide whether to use pixel-wise SZA
    if (cv_sza > CV_THRESHOLD) {
      cat(sprintf("  => CV > %.2f%%, using PIXEL-WISE solar zenith angles\n", CV_THRESHOLD * 100))
      use_pixel_sza <- TRUE
      mean_sza <- mean_sza_area
    } else {
      cat(sprintf("  => CV <= %.2f%%, using MEAN solar zenith angle: %.2f deg\n",
                  CV_THRESHOLD * 100, mean_sza_area))
      use_pixel_sza <- FALSE
      mean_sza <- mean_sza_area
      sza_raster <- NULL  # Don't need raster if using mean
    }
  }

} else {
  cat(sprintf("\nUsing calculated solar zenith angle: %.2f deg\n", mean_sza))
  use_pixel_sza <- FALSE
}

# Process viewing zenith angles if available
if (!is.null(angle_data$vza_raster)) {
  cat("\nProcessing viewing zenith angles...\n")
  l2w_raster_template <- spectral_raster_full[[1]]
  vza_raster_full <- resample(angle_data$vza_raster, l2w_raster_template, method = "bilinear")

  # Crop and mask to study area
  vza_raster <- crop(vza_raster_full, study_area_transformed)
  vza_raster <- mask(vza_raster, study_area_transformed)

  # Calculate statistics over study area
  vza_values <- values(vza_raster)
  vza_values <- vza_values[!is.na(vza_values)]

  if (length(vza_values) > 0) {
    mean_vza_area <- mean(vza_values)
    sd_vza_area <- sd(vza_values)
    cv_vza <- sd_vza_area / mean_vza_area

    cat(sprintf("VZA statistics over study area:\n"))
    cat(sprintf("  Range: %.2f - %.2f deg\n", min(vza_values), max(vza_values)))
    cat(sprintf("  Mean: %.2f deg\n", mean_vza_area))
    cat(sprintf("  SD: %.2f deg\n", sd_vza_area))
    cat(sprintf("  CV: %.4f (%.2f%%)\n", cv_vza, cv_vza * 100))

    # Decide whether to use pixel-wise VZA
    if (cv_vza > CV_THRESHOLD) {
      cat(sprintf("  => CV > %.2f%%, using PIXEL-WISE viewing zenith angles\n", CV_THRESHOLD * 100))
      use_pixel_vza <- TRUE
      mean_vza <- mean_vza_area
    } else {
      cat(sprintf("  => CV <= %.2f%%, using MEAN viewing zenith angle: %.2f deg\n",
                  CV_THRESHOLD * 100, mean_vza_area))
      use_pixel_vza <- FALSE
      mean_vza <- mean_vza_area
      vza_raster <- NULL  # Don't need raster if using mean
    }
  }

} else {
  cat(sprintf("\nUsing mean/fallback viewing zenith angle: %.2f deg\n", mean_vza))
  use_pixel_vza <- FALSE
}

# Summary
cat("\n========================================\n")
cat("ANGLE USAGE SUMMARY\n")
cat("========================================\n")
cat(sprintf("Solar Zenith Angle (theta_sun):\n"))
if (use_pixel_sza) {
  cat(sprintf("  Mode: PIXEL-WISE (high spatial variation detected)\n"))
  cat(sprintf("  Mean: %.2f deg (for reference)\n", mean_sza))
} else {
  cat(sprintf("  Mode: CONSTANT (low spatial variation)\n"))
  cat(sprintf("  Value: %.2f deg\n", mean_sza))
}

cat(sprintf("\nViewing Zenith Angle (theta_view):\n"))
if (use_pixel_vza) {
  cat(sprintf("  Mode: PIXEL-WISE (high spatial variation detected)\n"))
  cat(sprintf("  Mean: %.2f deg (for reference)\n", mean_vza))
} else {
  cat(sprintf("  Mode: CONSTANT (low spatial variation or nadir)\n"))
  cat(sprintf("  Value: %.2f deg\n", mean_vza))
}
cat("========================================\n\n")

cat(sprintf("Total valid pixels: %d\n", nrow(df)))

# === TEST MODE: USE SUBSET OF PIXELS ===
TEST_MODE <- FALSE
TEST_SIZE <- 10000  # Number of pixels to test

if (TEST_MODE && nrow(df) > TEST_SIZE) {
  cat("\n!!! TEST MODE ENABLED !!!\n")
  cat(sprintf("Reducing from %d to %d pixels for testing\n", nrow(df), TEST_SIZE))

  # Sample random pixels
  set.seed(42)
  df <- df %>% sample_n(TEST_SIZE)

  cat(sprintf("Test subset: %d pixels selected\n", nrow(df)))
}

# Extract wavelength information from band names
# ACOLITE L2W format: Rrs_443, Rrs_490, etc.
wavelengths <- as.numeric(gsub("^Rrs_|^rhos_", "", rrs_bands[1:6]))

cat("\nWavelengths detected:\n")
print(wavelengths)

# Pivot to long format and prepare for SABER
df_long <- df %>%
  dplyr::select(ensemble, starts_with("Rrs_") | starts_with("rhos_")) %>%
  pivot_longer(
    cols = starts_with("Rrs_") | starts_with("rhos_"),
    names_to = "wavelength_str",
    values_to = "rrs_0p"  # L2W is already water-leaving reflectance
  ) %>%
  mutate(
    wavelength = as.numeric(gsub("^Rrs_|^rhos_", "", wavelength_str)),
    rrs_0m = rrs_0p_to_0m(rrs_0p)  # Apply conversion here
  ) %>%
  dplyr::select(-wavelength_str)

# Filter out invalid pixels (negative reflectance)
valid_ensembles <- df_long %>%
  group_by(ensemble) %>%
  filter(all(rrs_0m >= 0, na.rm = TRUE)) %>%
  ungroup()

cat(sprintf("Valid pixels after filtering: %d\n",
            length(unique(valid_ensembles$ensemble))))

# Nest into format required by SABER
nested_df <- valid_ensembles %>%
  dplyr::select(ensemble, wavelength, rrs_0m) %>%
  group_by(ensemble) %>%
  nest(data = c(wavelength, rrs_0m)) %>%
  ungroup()

cat(sprintf("Nested dataframe prepared with %d ensembles\n", nrow(nested_df)))


## CONFIGURE SABER INVERSION PARAMETERS FOR SHALLOW WATER ----

cat("\nConfiguring SABER parameters for shallow water inversion...\n")

# Select benthic classes (MUST be done before inversion)
list_benthic_classes()
select_benthic_classes(c("Sand_2019","Eelgrass_2019", "Mud_2019"))

# ============================================================================
# INVERSION APPROACH SELECTION
# ============================================================================
# Choose one of three approaches for handling water depth (h_w):
#
# APPROACH 1: "unconstrained" - Full inversion without bathymetry
#   - Retrieves h_w along with all other parameters
#   - Uses wide bounds (0.5-12m) for h_w
#   - No bathymetry data required
#
# APPROACH 2: "fixed" - Bathymetry as fixed parameter
#   - Uses bathymetry to set h_w in par_fixed (NOT inverted)
#   - Does NOT retrieve h_w or h_w_sd
#   - Fastest approach, assumes bathymetry is accurate
#
# APPROACH 3: "constrained" - Bathymetry-guided optimization
#   - Uses bathymetry to set tighter bounds and initial values
#   - Still retrieves h_w (allows deviation from bathymetry)
#   - Retrieves h_w_sd for uncertainty quantification
#   - Best balance between constraint and flexibility
# ============================================================================

INVERSION_APPROACH <- "fixed"  # Options: "unconstrained", "fixed", "constrained"

cat(sprintf("Selected inversion approach: %s\n", INVERSION_APPROACH))

# Validate approach selection
if (!INVERSION_APPROACH %in% c("unconstrained", "fixed", "constrained")) {
  stop("Invalid INVERSION_APPROACH. Must be one of: 'unconstrained', 'fixed', 'constrained'")
}

# Check bathymetry availability for approaches that need it
if (INVERSION_APPROACH %in% c("fixed", "constrained") && !use_bathy_constraint) {
  cat("WARNING: Selected approach requires bathymetry but it's not available\n")
  cat("Falling back to 'unconstrained' approach\n")
  INVERSION_APPROACH <- "unconstrained"
}

# Fixed parameters (base set with mean solar and viewing zenith angles)
par_fixed_base <- c(
  "water_type" = 2,
  "theta_view" = mean_vza,   # Mean or constant viewing zenith angle
  "theta_sun" = mean_sza     # Mean or constant solar zenith angle
)

cat("\nBase par_fixed parameters:\n")
cat(sprintf("  water_type: %d\n", par_fixed_base["water_type"]))
cat(sprintf("  theta_view: %.2f deg\n", par_fixed_base["theta_view"]))
cat(sprintf("  theta_sun: %.2f deg\n", par_fixed_base["theta_sun"]))

# Add angle rasters to df if using pixel-wise angles
if (use_pixel_sza || use_pixel_vza) {
  cat("\nAdding pixel-wise angle data to dataframe...\n")

  if (use_pixel_sza) {
    sza_df <- as.data.frame(sza_raster, xy = TRUE, cells = TRUE, na.rm = FALSE)
    names(sza_df)[4] <- "pixel_sza"
    df <- df %>%
      left_join(sza_df %>% dplyr::select(cell, pixel_sza), by = "cell")
    cat(sprintf("  Added pixel_sza (%.1f%% coverage)\n",
                100 * sum(!is.na(df$pixel_sza)) / nrow(df)))
  }

  if (use_pixel_vza) {
    vza_df <- as.data.frame(vza_raster, xy = TRUE, cells = TRUE, na.rm = FALSE)
    names(vza_df)[4] <- "pixel_vza"
    df <- df %>%
      left_join(vza_df %>% dplyr::select(cell, pixel_vza), by = "cell")
    cat(sprintf("  Added pixel_vza (%.1f%% coverage)\n",
                100 * sum(!is.na(df$pixel_vza)) / nrow(df)))
  }
}

# ============================================================================
# HYBRID INVERSION: CONFIGURE PARAMETERS FOR SHALLOW AND DEEP WATER
# ============================================================================
# Pixels with depth 0.5-10m: Optically shallow (with benthic reflectance)
# Pixels with depth >10m: Optically deep (no benthic reflectance)
# Pixels with depth <0.5m: Skipped (already filtered out)
# ============================================================================

cat("\n=== SHALLOW WATER PARAMETERS (0.5-10m) ===\n")

# Configure parameters based on approach
if (INVERSION_APPROACH == "fixed") {
  cat("Approach: FIXED - h_w will be set from bathymetry (not inverted)\n")

  # h_w will be added to par_fixed per-pixel
  # Parameters to invert (EXCLUDING h_w)
  par_to_inverse_shallow <- c(
    "chl", "a_g_440", "a_g_s", "bb_p_550", "bb_p_gamma",
    "r_rs_b_Sand_2019", "r_rs_b_Eelgrass_2019", "r_rs_b_Mud_2019",
    "sd"
  )

  # Bounds (WITHOUT h_w)
  lower_shallow <- c(
    0.4,   # chl
    0.1,   # a_g_440
    0.005, # a_g_s
    0.002, # bb_p_550
    0.2,   # bb_p_gamma
    0.01,  # r_rs_b_Sand_2019
    0.01,  # r_rs_b_Eelgrass_2019
    0.01,  # r_rs_b_Mud_2019
    0.0001 # sd
  )

  best_shallow <- c(
    3.5,     # chl
    0.75,  # a_g_440
    0.017, # a_g_s
    0.005, # bb_p_550
    0.46,  # bb_p_gamma
    0.6,   # r_rs_b_Sand_2019
    0.15,   # r_rs_b_Eelgrass_2019
    0.25,   # r_rs_b_Mud_2019
    0.1    # sd
  )

  upper_shallow <- c(
    30,    # chl
    4.0,   # a_g_440
    0.022, # a_g_s
    0.010, # bb_p_550
    1,     # bb_p_gamma
    1.0,   # r_rs_b_Sand_2019
    1.0,   # r_rs_b_Eelgrass_2019
    1.0,   # r_rs_b_Mud_2019
    10     # sd
  )

} else {
  # UNCONSTRAINED or CONSTRAINED approaches (both invert h_w)
  cat(sprintf("Approach: %s - h_w will be retrieved\n", toupper(INVERSION_APPROACH)))

  # Parameters to invert (INCLUDING h_w)
  par_to_inverse_shallow <- c(
    "chl", "a_g_440", "a_g_s", "bb_p_550", "bb_p_gamma",
    "h_w",  # Water depth
    "r_rs_b_Sand_2019", "r_rs_b_Eelgrass_2019", "r_rs_b_Mud_2019",
    "sd"
  )

  # Bounds (WITH h_w)
  lower_shallow <- c(
    0.4,   # chl
    0.25,   # a_g_440
    0.005, # a_g_s
    0.002, # bb_p_550
    0.2,   # bb_p_gamma
    0.5,   # h_w (min depth 0.5m)
    0.01,  # r_rs_b_Sand_2019
    0.01,  # r_rs_b_Eelgrass_2019
    0.01,  # r_rs_b_Mud_2019
    0.0001 # sd
  )

  best_shallow <- c(
    3.5,     # chl
    1.15,  # a_g_440
    0.017, # a_g_s
    0.005, # bb_p_550
    0.46,  # bb_p_gamma
    4,     # h_w (initial guess 4m)
    0.3,   # r_rs_b_Sand_2019
    0.5,   # r_rs_b_Eelgrass_2019
    0.2,   # r_rs_b_Mud_2019
    0.1    # sd
  )

  upper_shallow <- c(
    30,    # chl
    4.0,   # a_g_440
    0.022, # a_g_s
    0.010, # bb_p_550
    1,     # bb_p_gamma
    10,    # h_w (max depth 10m for shallow)
    1.0,   # r_rs_b_Sand_2019
    1.0,   # r_rs_b_Eelgrass_2019
    1.0,   # r_rs_b_Mud_2019
    10     # sd
  )
}

cat(sprintf("  Parameters to invert: %d\n", length(par_to_inverse_shallow)))

cat("\n=== DEEP WATER PARAMETERS (>10m) ===\n")
cat("Optically deep - no benthic reflectance contribution\n")

# Deep water parameters (NO benthic reflectance, NO h_w)
par_to_inverse_deep <- c(
  "chl", "a_g_440", "a_g_s", "bb_p_550", "bb_p_gamma", "sd"
)

lower_deep <- c(
  0.4,   # chl
  0.25,   # a_g_440
  0.005, # a_g_s
  0.002, # bb_p_550
  0.2,   # bb_p_gamma
  0.0001 # sd
)

best_deep <- c(
  3.5,     # chl
  1.05,  # a_g_440
  0.017, # a_g_s
  0.005, # bb_p_550
  0.46,  # bb_p_gamma
  0.1    # sd
)

upper_deep <- c(
  12,    # chl
  2.5,   # a_g_440
  0.022, # a_g_s
  0.010, # bb_p_550
  1,     # bb_p_gamma
  10     # sd
)

cat(sprintf("  Parameters to invert: %d\n", length(par_to_inverse_deep)))
cat("  Fixed parameters: h_w=NULL, r_rs_b_*=NULL (optically deep)\n")

cat("\n========================================\n")


## PARALLEL EXECUTION SETUP ----

# Use 3 cores for test mode to avoid interfering with PID 22500 (8 cores)
num_cores <- if (exists("TEST_MODE") && TEST_MODE) 3 else (parallel::detectCores() - 1)
cat(sprintf("\nSetting up parallel processing with %d cores...\n", num_cores))
if (exists("TEST_MODE") && TEST_MODE) {
  cat("  (Limited to 3 cores for test mode - PID 22500 using 8 cores)\n")
}

cl <- makeCluster(num_cores)
doParallel::registerDoParallel(cl)

# Partition into n roughly equal-sized chunks
chunked_nested <- split(nested_df,
                        cut(seq_len(nrow(nested_df)), breaks = num_cores, labels = FALSE))

handlers("rstudio")  # Best for RStudio console

#gc()


## GRADIENT-BASED INVERSION (FASTER) ----

cat("\n========================================\n")
cat(sprintf("INVERSION APPROACH: %s\n", toupper(INVERSION_APPROACH)))
cat("========================================\n")

if (INVERSION_APPROACH == "unconstrained") {
  cat("Running unconstrained inversion (no bathymetry constraint)\n")
  cat("  - h_w will be retrieved with wide bounds (0.5-12m)\n")
} else if (INVERSION_APPROACH == "fixed") {
  cat("Running fixed-depth inversion (bathymetry-based)\n")
  cat("  - h_w set from bathymetry (NOT retrieved)\n")
  cat("  - Fastest approach\n")
} else if (INVERSION_APPROACH == "constrained") {
  cat("Running constrained inversion (bathymetry-guided)\n")
  cat("  - h_w retrieved with tighter bounds (±30% of bathymetry)\n")
  cat("  - Better initial values from bathymetry\n")
}

cat("\nStarting gradient-based inversion...\n")
cat("This may take several hours depending on image size.\n\n")

exec_time_parall <- system.time({
  with_progress({
    p <- progressor(steps = length(chunked_nested))

    if (INVERSION_APPROACH == "unconstrained") {
      # =====================================================================
      # APPROACH 1: UNCONSTRAINED - HYBRID shallow/deep inversion
      # =====================================================================
      # Shallow water (0.5-10m): Retrieve h_w with wide bounds, include benthic
      # Deep water (>10m): No h_w, no benthic reflectance (optically deep)
      # =====================================================================
      saber_results <- foreach(chunk = chunked_nested,
                               .packages = c("purrr", "SABER", "dplyr", "progressr", "numDeriv", "MASS"),
                               .export = c("inverse_gradient", "par_to_inverse_shallow", "par_to_inverse_deep",
                                           "lower_shallow", "best_shallow", "upper_shallow",
                                           "lower_deep", "best_deep", "upper_deep",
                                           "par_fixed_base", "df")) %dopar% {

                                 result <- chunk %>%
                                   mutate(
                                     inversion_estim = purrr::map2(data, ensemble, function(rrs_data, ens_id) {

                                       # Get water type classification for this pixel
                                       pixel_row <- df[df$ensemble == ens_id, ]
                                       water_type <- pixel_row$water_type_class[1]

                                       # HYBRID INVERSION: Switch parameters based on water type
                                       if (water_type == "shallow") {
                                         # SHALLOW WATER: Retrieve h_w with wide bounds
                                         inverse_gradient(
                                           rrs = rrs_data,
                                           forward_model = "am03",
                                           objective_fct = "log-ll",
                                           optim_mtd = "L-BFGS-B",
                                           par_inversed = par_to_inverse_shallow,  # 10 params (with h_w)
                                           par_fixed = par_fixed_base,
                                           lower_b = lower_shallow,
                                           init_val = best_shallow,
                                           upper_b = upper_shallow,
                                           verbose = FALSE
                                         )

                                       } else {
                                         # DEEP WATER: No h_w, no benthic
                                         inverse_gradient(
                                           rrs = rrs_data,
                                           forward_model = "am03",
                                           objective_fct = "log-ll",
                                           optim_mtd = "L-BFGS-B",
                                           par_inversed = par_to_inverse_deep,  # 6 params
                                           par_fixed = par_fixed_base,
                                           lower_b = lower_deep,
                                           init_val = best_deep,
                                           upper_b = upper_deep,
                                           verbose = FALSE
                                         )
                                       }
                                     })
                                   )

                                 p()
                                 result
                               }

    } else if (INVERSION_APPROACH == "fixed") {
      # =====================================================================
      # APPROACH 2: FIXED - HYBRID shallow/deep inversion
      # =====================================================================
      # Shallow water (0.5-10m): Use bathymetry as fixed h_w, include benthic
      # Deep water (>10m): No h_w, no benthic reflectance (optically deep)
      # =====================================================================
      saber_results <- foreach(chunk = chunked_nested,
                               .packages = c("purrr", "SABER", "dplyr", "progressr", "numDeriv", "MASS"),
                               .export = c("inverse_gradient", "par_to_inverse_shallow", "par_to_inverse_deep",
                                           "lower_shallow", "best_shallow", "upper_shallow",
                                           "lower_deep", "best_deep", "upper_deep",
                                           "par_fixed_base", "df", "use_pixel_sza", "use_pixel_vza")) %dopar% {

                                 result <- chunk %>%
                                   mutate(
                                     inversion_estim = purrr::map2(data, ensemble, function(rrs_data, ens_id) {

                                       # Get water type classification and bathymetry for this pixel
                                       pixel_row <- df[df$ensemble == ens_id, ]
                                       water_type <- pixel_row$water_type_class[1]
                                       pixel_bathy <- pixel_row$bathy_depth[1]

                                       # Create pixel-specific par_fixed with angles
                                       par_fixed_pixel <- par_fixed_base

                                       # Update with pixel-wise angles if available
                                       if (use_pixel_sza && "pixel_sza" %in% names(pixel_row)) {
                                         if (!is.na(pixel_row$pixel_sza[1])) {
                                           par_fixed_pixel["theta_sun"] <- pixel_row$pixel_sza[1]
                                         }
                                       }

                                       if (use_pixel_vza && "pixel_vza" %in% names(pixel_row)) {
                                         if (!is.na(pixel_row$pixel_vza[1])) {
                                           par_fixed_pixel["theta_view"] <- pixel_row$pixel_vza[1]
                                         }
                                       }

                                       # HYBRID INVERSION: Switch parameters based on water type
                                       if (water_type == "shallow") {
                                         # SHALLOW WATER: Use shallow parameters with benthic
                                         par_fixed_pixel["h_w"] <- pixel_bathy  # Set from bathymetry

                                         inverse_gradient(
                                           rrs = rrs_data,
                                           forward_model = "am03",
                                           objective_fct = "log-ll",
                                           optim_mtd = "L-BFGS-B",
                                           par_inversed = par_to_inverse_shallow,  # 9 params (no h_w)
                                           par_fixed = par_fixed_pixel,             # includes h_w + angles
                                           lower_b = lower_shallow,
                                           init_val = best_shallow,
                                           upper_b = upper_shallow,
                                           verbose = FALSE
                                         )

                                       } else {
                                         # DEEP WATER: Use deep parameters without benthic
                                         # h_w = NULL (optically deep, benthic not relevant)
                                         # Use par_fixed_pixel with angles (no h_w modification needed)
                                         inverse_gradient(
                                           rrs = rrs_data,
                                           forward_model = "am03",
                                           objective_fct = "log-ll",
                                           optim_mtd = "L-BFGS-B",
                                           par_inversed = par_to_inverse_deep,  # 6 params (no h_w, no benthic)
                                           par_fixed = par_fixed_pixel,          # with angles
                                           lower_b = lower_deep,
                                           init_val = best_deep,
                                           upper_b = upper_deep,
                                           verbose = FALSE
                                         )
                                       }
                                     })
                                   )

                                 p()
                                 result
                               }

    } else if (INVERSION_APPROACH == "constrained") {
      # =====================================================================
      # APPROACH 3: CONSTRAINED - HYBRID shallow/deep inversion
      # =====================================================================
      # Shallow water (0.5-10m): Retrieve h_w with tight bounds (±30%), include benthic
      # Deep water (>10m): No h_w, no benthic reflectance (optically deep)
      # =====================================================================
      saber_results <- foreach(chunk = chunked_nested,
                               .packages = c("purrr", "SABER", "dplyr", "progressr", "numDeriv", "MASS"),
                               .export = c("inverse_gradient", "par_to_inverse_shallow", "par_to_inverse_deep",
                                           "lower_shallow", "best_shallow", "upper_shallow",
                                           "lower_deep", "best_deep", "upper_deep",
                                           "par_fixed_base", "df")) %dopar% {

                                 result <- chunk %>%
                                   mutate(
                                     inversion_estim = purrr::map2(data, ensemble, function(rrs_data, ens_id) {

                                       # Get water type classification and bathymetry for this pixel
                                       pixel_row <- df[df$ensemble == ens_id, ]
                                       water_type <- pixel_row$water_type_class[1]
                                       pixel_bathy <- pixel_row$bathy_depth[1]

                                       # HYBRID INVERSION: Switch parameters based on water type
                                       if (water_type == "shallow") {
                                         # SHALLOW WATER: Constrained h_w retrieval (±30% of bathymetry)
                                         # Create pixel-specific bounds and initial value for h_w
                                         h_w_lower <- max(0.5, pixel_bathy * 0.7)
                                         h_w_upper <- min(10, pixel_bathy * 1.3)
                                         h_w_init <- pixel_bathy

                                         # Update bounds for this pixel (h_w is 6th parameter in shallow)
                                         lower_pixel <- lower_shallow
                                         lower_pixel[6] <- h_w_lower

                                         upper_pixel <- upper_shallow
                                         upper_pixel[6] <- h_w_upper

                                         best_pixel <- best_shallow
                                         best_pixel[6] <- h_w_init

                                         inverse_gradient(
                                           rrs = rrs_data,
                                           forward_model = "am03",
                                           objective_fct = "log-ll",
                                           optim_mtd = "L-BFGS-B",
                                           par_inversed = par_to_inverse_shallow,  # 10 params (with h_w)
                                           par_fixed = par_fixed_pixel,
                                           lower_b = lower_pixel,
                                           init_val = best_pixel,
                                           upper_b = upper_pixel,
                                           verbose = FALSE
                                         )

                                       } else {
                                         # DEEP WATER: No h_w, no benthic
                                         inverse_gradient(
                                           rrs = rrs_data,
                                           forward_model = "am03",
                                           objective_fct = "log-ll",
                                           optim_mtd = "L-BFGS-B",
                                           par_inversed = par_to_inverse_deep,  # 6 params
                                           par_fixed = par_fixed_base,
                                           lower_b = lower_deep,
                                           init_val = best_deep,
                                           upper_b = upper_deep,
                                           verbose = FALSE
                                         )
                                       }
                                     })
                                   )

                                 p()
                                 result
                               }
    }
  })
})

stopCluster(cl)

cat(sprintf("\nInversion completed in %.2f seconds (%.2f minutes)\n",
            exec_time_parall[3], exec_time_parall[3]/60))


## PROCESS AND ORGANIZE RESULTS ----

cat("\nProcessing inversion results...\n")

# Combine results from all chunks
combined_results <- bind_rows(saber_results)

ensemble_results <- combined_results %>%
  dplyr::select(ensemble, inversion_estim)

# Unnest results
ensemble_results_tidy <- ensemble_results %>%
  mutate(inversion_estim = map(inversion_estim, ~ as_tibble(t(.x)))) %>%
  unnest(cols = c(inversion_estim))

# Join with spatial information
final_results <- df %>%
  dplyr::select(ensemble, x, y, row, col) %>%
  distinct() %>%
  inner_join(ensemble_results_tidy, by = "ensemble")

# Add water_type_class for downstream analysis and visualization
final_results <- final_results %>%
  left_join(df %>% dplyr::select(ensemble, water_type_class), by = "ensemble")

cat(sprintf("Final results contain %d pixels\n", nrow(final_results)))


### NORMALIZE BENTHIC REFLECTANCE FRACTIONS (OPTIONAL) ----
cat("\nNormalizing benthic reflectance fractions for optically shallow water pixels...\n")

# Normalize r_rs_b fractions (sum to 1.0) - only for shallow water
final_results <- final_results %>%
  mutate(
    # Calculate row-wise sum of benthic fractions (only if columns exist)
    r_rs_b_sum = ifelse(
      !is.na(r_rs_b_Sand_2019),
      r_rs_b_Sand_2019 + r_rs_b_Eelgrass_2019 + r_rs_b_Mud_2019,
      NA_real_
    ),

    # Normalize each fraction by dividing by the sum (NA for deep water)
    r_rs_b_Sand_2019_normalized = ifelse(
      !is.na(r_rs_b_sum) & r_rs_b_sum > 0,
      r_rs_b_Sand_2019 / r_rs_b_sum,
      NA_real_
    ),
    r_rs_b_Eelgrass_2019_normalized = ifelse(
      !is.na(r_rs_b_sum) & r_rs_b_sum > 0,
      r_rs_b_Eelgrass_2019 / r_rs_b_sum,
      NA_real_
    ),
    r_rs_b_Mud_2019_normalized = ifelse(
      !is.na(r_rs_b_sum) & r_rs_b_sum > 0,
      r_rs_b_Mud_2019 / r_rs_b_sum,
      NA_real_
    )
  )

# Normalize r_rs_b_sd fractions (sum to 1.0) - only for shallow water
final_results <- final_results %>%
  mutate(
    # Calculate row-wise sum of benthic fraction uncertainties (only if columns exist)
    r_rs_b_sd_sum = ifelse(
      !is.na(r_rs_b_Sand_2019_sd),
      r_rs_b_Sand_2019_sd + r_rs_b_Eelgrass_2019_sd + r_rs_b_Mud_2019_sd,
      NA_real_
    ),

    # Normalize each uncertainty by dividing by the sum (NA for deep water)
    r_rs_b_Sand_2019_sd_normalized = ifelse(
      !is.na(r_rs_b_sd_sum) & r_rs_b_sd_sum > 0,
      r_rs_b_Sand_2019_sd / r_rs_b_sd_sum,
      NA_real_
    ),
    r_rs_b_Eelgrass_2019_sd_normalized = ifelse(
      !is.na(r_rs_b_sd_sum) & r_rs_b_sd_sum > 0,
      r_rs_b_Eelgrass_2019_sd / r_rs_b_sd_sum,
      NA_real_
    ),
    r_rs_b_Mud_2019_sd_normalized = ifelse(
      !is.na(r_rs_b_sd_sum) & r_rs_b_sd_sum > 0,
      r_rs_b_Mud_2019_sd / r_rs_b_sd_sum,
      NA_real_
    )
  ) %>%
  # Remove temporary sum columns
  dplyr::select(-r_rs_b_sum, -r_rs_b_sd_sum)

cat("Normalized benthic fractions created \n")

# Verify normalization (should sum to 1.0)
normalization_check <- final_results %>%
  mutate(
    norm_sum = r_rs_b_Sand_2019_normalized + r_rs_b_Eelgrass_2019_normalized + r_rs_b_Mud_2019_normalized,
    norm_sd_sum = r_rs_b_Sand_2019_sd_normalized + r_rs_b_Eelgrass_2019_sd_normalized + r_rs_b_Mud_2019_sd_normalized
  )

cat(sprintf("\nNormalization verification:\n"))
cat(sprintf("  Mean normalized fraction sum: %.6f (should be ~1.0)\n",
            mean(normalization_check$norm_sum, na.rm = TRUE)))
cat(sprintf("  Mean normalized sd sum: %.6f (should be ~1.0)\n",
            mean(normalization_check$norm_sd_sum, na.rm = TRUE)))


### CREATE RASTERS FROM INVERSION RESULTS ----
# Check if in test mode
if (exists("TEST_MODE") && TEST_MODE) {
  cat("\n!!! TEST MODE: Skipping raster creation !!!\n")
  cat("Sampled pixels don't form a regular grid - saving as CSV instead\n")

  # Save results as CSV
  csv_output <- sprintf("./tests/sat/s2_l2b_mp/inversion_results_test_%s.csv", INVERSION_APPROACH)
  write.csv(final_results, csv_output, row.names = FALSE)
  cat(sprintf("Results saved to: %s\n", csv_output))

  # Print summary statistics
  cat("\n========================================\n")
  cat("SHALLOW WATER INVERSION SUMMARY (TEST MODE)\n")
  cat("========================================\n")
  cat(sprintf("Total execution time: %.2f minutes\n", exec_time_parall[3]/60))
  cat(sprintf("Pixels processed: %d\n", nrow(final_results)))
  cat(sprintf("Output CSV: %s\n", csv_output))

  # Show variable statistics
  cat("\nVariable Statistics:\n")
  summary_vars <- c("chl", "a_g_440", "bb_p_550",
                    "r_rs_b_Sand_2019_normalized", "r_rs_b_Eelgrass_2019_normalized", "r_rs_b_Mud_2019_normalized")
  # Add h_w if it exists (not in "fixed" approach)
  if ("h_w" %in% names(final_results)) {
    summary_vars <- c(summary_vars, "h_w")
  }
  summary_stats <- final_results %>%
    dplyr::select(any_of(summary_vars)) %>%
    summary()
  print(summary_stats)

  cat("========================================\n")
  cat("TEST COMPLETE - Disable TEST_MODE for full processing\n")
  cat("========================================\n")

  # Exit script here
  quit(save = "no", status = 0)
}

cat("\nCreating rasters from inversion results...\n")

# Prepare empty raster template
r_template <- rast(nrows = nrows, ncols = ncols,
                   ext = ext(spectral_raster),
                   crs = crs(spectral_raster))

# Variables to create rasters for (exclude sd and sd_sd)
vars_to_raster <- setdiff(names(final_results), c("ensemble", "x", "y", "row",
                                                  "col", "cell", "sd", "sd_sd"))

cat("\nVariables to rasterize:\n")
print(vars_to_raster)

# Create a list to store rasters
raster_list <- list()


### CREATE RASTERS FOR INVERSION VARIABLES ----
# Create rasters for all inversion output variables
for (varname in vars_to_raster) {
  r <- rast(r_template)
  vals <- rep(NA_real_, ncell(r))
  cell_indices <- cellFromRowCol(r, final_results$row, final_results$col)
  vals[cell_indices] <- final_results[[varname]]
  values(r) <- vals
  names(r) <- varname
  raster_list[[varname]] <- r
}

# Stack all layers into a SpatRaster
r_stack <- rast(raster_list)
names(r_stack) <- names(raster_list)

cat(sprintf("Created raster stack with %d layers\n", nlyr(r_stack)))


## SAVE CROPPED INPUT RASTER TO NETCDF ----

cat("\n========================================\n")
cat("SAVING CROPPED INPUT RASTER\n")
cat("========================================\n")

output_nc_cropped <- sprintf("./tests/sat/s2_l2b_mp/%s_cropped.nc",
                             s2_img_name)

cat(sprintf("\nSaving cropped input raster to: %s\n", output_nc_cropped))

# Get dimensions from cropped spectral raster
nx_crop <- ncol(spectral_raster)
ny_crop <- nrow(spectral_raster)
ext_vals_crop <- ext(spectral_raster)

# Create dimension variables
xvals_crop <- seq(ext_vals_crop[1], ext_vals_crop[2], length.out = nx_crop)
yvals_crop <- seq(ext_vals_crop[3], ext_vals_crop[4], length.out = ny_crop)

xdim_crop <- ncdim_def("x", "meters", xvals_crop)
ydim_crop <- ncdim_def("y", "meters", yvals_crop)

# Create variable definitions for all Rrs bands
var_list_crop <- list()

cat(sprintf("Creating NetCDF variable definitions for %d Rrs bands:\n", nlyr(spectral_raster)))

for (i in 1:nlyr(spectral_raster)) {
  var_name <- names(spectral_raster)[i]
  cat(sprintf("  %d. %s\n", i, var_name))

  var_list_crop[[i]] <- ncvar_def(
    name = var_name,
    units = "sr^-1",
    dim = list(xdim_crop, ydim_crop),
    missval = -9999,
    longname = paste("Remote sensing reflectance at", gsub("Rrs_", "", var_name), "nm"),
    prec = "float"
  )
}

# Create lat/lon coordinate variables for cropped raster
cat("\nCreating lat/lon coordinate variables for cropped raster...\n")

# Create grid of x,y coordinates in projected CRS
xy_grid_crop <- expand.grid(x = xvals_crop, y = yvals_crop)
xy_matrix_crop <- as.matrix(xy_grid_crop)

# Project to WGS84 lat/lon
cat("  Projecting coordinates from UTM to WGS84...\n")
latlon_coords_crop <- project(xy_matrix_crop, from = crs(spectral_raster), to = "EPSG:4326")

# Reshape to 2D grids matching raster matrix format
# Matrix should be (ny x nx) with rows=y, cols=x
lon_grid_crop <- matrix(latlon_coords_crop[, 1], nrow = ny_crop, ncol = nx_crop, byrow = FALSE)
lat_grid_crop <- matrix(latlon_coords_crop[, 2], nrow = ny_crop, ncol = nx_crop, byrow = FALSE)

cat(sprintf("  Lon range: %.6f to %.6f\n", min(lon_grid_crop), max(lon_grid_crop)))
cat(sprintf("  Lat range: %.6f to %.6f\n", min(lat_grid_crop), max(lat_grid_crop)))

# Define lat/lon as 2D coordinate variables
lon_var_crop <- ncvar_def(
  name = "lon",
  units = "degrees_east",
  dim = list(xdim_crop, ydim_crop),
  missval = -9999,
  longname = "Longitude (WGS84)",
  prec = "double"
)

lat_var_crop <- ncvar_def(
  name = "lat",
  units = "degrees_north",
  dim = list(xdim_crop, ydim_crop),
  missval = -9999,
  longname = "Latitude (WGS84)",
  prec = "double"
)

# Add lat/lon to variable list
var_list_crop[[length(var_list_crop) + 1]] <- lon_var_crop
var_list_crop[[length(var_list_crop) + 1]] <- lat_var_crop

cat(sprintf("  Total variables (including lat/lon): %d\n", length(var_list_crop)))

# Create the NetCDF file
ncout_crop <- nc_create(output_nc_cropped, var_list_crop, force_v4 = TRUE)

# Write data for each Rrs band
cat("\nWriting Rrs bands to NetCDF file:\n")
for (i in 1:nlyr(spectral_raster)) {
  var_name <- names(spectral_raster)[i]
  cat(sprintf("  Writing %s...\n", var_name))

  # Extract values and convert to matrix
  vals <- values(spectral_raster[[i]], mat = TRUE)
  vals[is.na(vals)] <- -9999

  # Write to NetCDF
  ncvar_put(ncout_crop, var_list_crop[[i]], vals)
}

# Write lat/lon coordinate grids
cat("  Writing lon coordinate grid...\n")
ncvar_put(ncout_crop, lon_var_crop, lon_grid_crop)

cat("  Writing lat coordinate grid...\n")
ncvar_put(ncout_crop, lat_var_crop, lat_grid_crop)

# Add global attributes
ncatt_put(ncout_crop, 0, "title", "SABER Cropped Input Raster (Rrs)")
ncatt_put(ncout_crop, 0, "source", "Sentinel-2 L2W ACOLITE - Cropped to study area")
ncatt_put(ncout_crop, 0, "date_created", as.character(Sys.time()))
ncatt_put(ncout_crop, 0, "crs", as.character(crs(spectral_raster)))
ncatt_put(ncout_crop, 0, "original_file", s2_l2w_file)

# Add coordinate system attributes
ncatt_put(ncout_crop, "x", "standard_name", "projection_x_coordinate")
ncatt_put(ncout_crop, "x", "long_name", "x coordinate of projection")
ncatt_put(ncout_crop, "x", "axis", "X")

ncatt_put(ncout_crop, "y", "standard_name", "projection_y_coordinate")
ncatt_put(ncout_crop, "y", "long_name", "y coordinate of projection")
ncatt_put(ncout_crop, "y", "axis", "Y")

ncatt_put(ncout_crop, "lon", "standard_name", "longitude")
ncatt_put(ncout_crop, "lon", "long_name", "Longitude")
ncatt_put(ncout_crop, "lon", "axis", "X")

ncatt_put(ncout_crop, "lat", "standard_name", "latitude")
ncatt_put(ncout_crop, "lat", "long_name", "Latitude")
ncatt_put(ncout_crop, "lat", "axis", "Y")

# Close the file
nc_close(ncout_crop)

cat(sprintf("\nCropped input raster saved to: %s\n", output_nc_cropped))
cat(sprintf("Total Rrs bands written: %d (plus lat/lon coordinate grids)\n", nlyr(spectral_raster)))


## PLOT INPUT RGB IMAGE ----

cat("\nGenerating RGB plot of input image...\n")

# Select bands for RGB visualization (typically Red, Green, Blue)
# For Sentinel-2: R=665nm, G=560nm, B=490nm (approximate)
rgb_band_names <- c("Rrs_665", "Rrs_560", "Rrs_490")

# Check if these bands exist, otherwise use closest available
available_rrs <- names(spectral_raster)
if (!all(rgb_band_names %in% available_rrs)) {
  cat("  Note: Exact RGB bands not found, using available bands\n")
  # Use first 3 bands if exact match not found
  rgb_bands <- spectral_raster[[1:min(3, nlyr(spectral_raster))]]
} else {
  rgb_bands <- spectral_raster[[rgb_band_names]]
}

# Convert to dataframe for ggplot
df_rgb <- as.data.frame(rgb_bands, xy = TRUE, na.rm = TRUE)

# Normalize each band to 0-1 range for RGB display
for (i in 3:ncol(df_rgb)) {
  band_values <- df_rgb[[i]]
  # Use 2nd and 98th percentile for better contrast
  q_low <- quantile(band_values, 0.02, na.rm = TRUE)
  q_high <- quantile(band_values, 0.98, na.rm = TRUE)
  df_rgb[[i]] <- (pmin(pmax(band_values, q_low), q_high) - q_low) / (q_high - q_low)
}

# Create RGB composite (use first 3 bands as R, G, B)
df_rgb$R <- df_rgb[[3]]
df_rgb$G <- df_rgb[[4]]
df_rgb$B <- df_rgb[[5]]

# Create RGB plot
p_rgb <- ggplot(df_rgb, aes(x = x, y = y)) +
  geom_raster(aes(fill = rgb(R, G, B))) +
  scale_fill_identity() +
  coord_equal() +
  labs( subtitle = sprintf("Bands: %s (R), %s (G), %s (B)",
                          names(df_rgb)[3], names(df_rgb)[4], names(df_rgb)[5]),
       x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    panel.grid.major = element_line(color = "gray80", size = 0.3),
    plot.margin = margin(2, 2, 2, 2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none"  # No legend needed for RGB
  ) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering)

# Save RGB plot
rgb_output <- file.path("./tests/sat/s2_l2b_mp", "inv_shallow_RGB_input.png")
ggsave(rgb_output, p_rgb, units = "in", width = 8, height = 6, dpi = 300)

cat(sprintf("RGB plot saved to: %s\n", rgb_output))


## PLOT INVERSION RESULTS ----

#use_quantiles <- c(0.01, 0.95)  # Use quantile-based approach
use_quantiles <- NA             # Use histogram equalization (old method)

cat("\nGenerating plots...\n")

if (is.na(use_quantiles[1])) {
  cat("Using histogram equalization (old method)\n")
} else {
  cat(sprintf("Using quantile-based color scaling: %.0f%% - %.0f%%\n",
              use_quantiles[1] * 100, use_quantiles[2] * 100))
}

plot_list_eq <- list()

# Select key variables to plot
vars_to_plot <- c("chl", "a_g_440", "bb_p_550",
                  "r_rs_b_Eelgrass_2019", "r_rs_b_Sand_2019", "r_rs_b_Mud_2019",
                  "r_rs_b_Eelgrass_2019_normalized", "r_rs_b_Sand_2019_normalized", "r_rs_b_Mud_2019_normalized",
                  "chl_sd", "a_g_440_sd", "bb_p_550_sd",
                  "r_rs_b_Sand_2019_sd", "r_rs_b_Eelgrass_2019_sd", "r_rs_b_Mud_2019_sd",
                  "r_rs_b_Sand_2019_sd_normalized", "r_rs_b_Eelgrass_2019_sd_normalized", "r_rs_b_Mud_2019_sd_normalized")

# Add h_w and h_w_sd if available (not available in "fixed" approach)
if (INVERSION_APPROACH != "fixed") {
  vars_to_plot <- c(vars_to_plot, "h_w", "h_w_sd")
}

for (varname in vars_to_plot) {
  if (!varname %in% names(raster_list)) {
    cat(sprintf("Skipping %s (not in results)\n", varname))
    next
  }

  cat(sprintf("Plotting %s...\n", varname))

  r <- raster_list[[varname]]
  df_plot <- raster_to_df(r)

  # Remove NA values
  df_plot <- df_plot[!is.na(df_plot$value), ]

  if (nrow(df_plot) == 0) {
    cat(sprintf("  Warning: No valid values for %s, skipping plot\n", varname))
    next
  }

  # Branch based on use_quantiles setting
  if (is.na(use_quantiles[1])) {

    # Equalized color mapping for better visualization
    ecdf_vals <- ecdf(df_plot$value)
    df_plot$value_eq <- ecdf_vals(df_plot$value)

    # Create breaks for legend
    breaks <- seq(from = min(df_plot$value, na.rm = TRUE),
                  to = max(df_plot$value, na.rm = TRUE),
                  length.out = 5)
    labels <- format(breaks, digits = 3)

    p <- ggplot(df_plot, aes(x = x, y = y, fill = value_eq)) +
      geom_raster() +
      scale_fill_gradientn(
        colors = viridis::viridis(length(breaks)),
        breaks = c(0, 0.25, 0.50, 0.75, 1),
        labels = labels,
        name = varname,
        na.value = "gray90"
      ) +
      coord_equal() +
      labs(title = paste0("SABER: ", varname),
           x = "Longitude", y = "Latitude") +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        panel.grid.major = element_line(color = "gray80", size = 0.3),
        plot.margin = margin(2, 2, 2, 2),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      ) +
      annotation_scale(location = "bl", width_hint = 0.3) +
      annotation_north_arrow(location = "tl", which_north = "true",
                             style = north_arrow_fancy_orienteering)

  } else {

    LOWER_QUANTILE <- use_quantiles[1]
    UPPER_QUANTILE <- use_quantiles[2]

    # Calculate quantile-based limits for color scale
    all_values <- df_plot$value
    all_values <- all_values[is.finite(all_values)]

    if (length(all_values) == 0) {
      cat(sprintf("  Warning: No finite values for %s, skipping plot\n", varname))
      next
    }

    # Use quantiles for color scale limits
    value_min <- quantile(all_values, LOWER_QUANTILE, na.rm = TRUE)
    value_max <- quantile(all_values, UPPER_QUANTILE, na.rm = TRUE)
    value_median <- median(all_values, na.rm = TRUE)

    # Actual min/max for reference
    actual_min <- min(all_values, na.rm = TRUE)
    actual_max <- max(all_values, na.rm = TRUE)

    cat(sprintf("  %s statistics:\n", varname))
    cat(sprintf("    Actual range: %.3e to %.3e\n", actual_min, actual_max))
    cat(sprintf("    Median: %.3e\n", value_median))
    cat(sprintf("    Color scale (%.0f%%-%.0f%% quantile): %.3e to %.3e\n",
                LOWER_QUANTILE * 100, UPPER_QUANTILE * 100, value_min, value_max))

    # Clip values to quantile range for visualization
    df_plot$value_clipped <- pmin(pmax(df_plot$value, value_min), value_max)

    # Create plot with viridis color scale
    p <- ggplot(df_plot, aes(x = x, y = y, fill = value_clipped)) +
      geom_raster() +
      scale_fill_viridis_c(
        option = "viridis",
        name = varname,
        limits = c(value_min, value_max),
        breaks = pretty(c(value_min, value_max), n = 6),
        na.value = "gray90"
      ) +
      coord_equal() +
      labs(title = paste(varname, "- Shallow Water Inversion"),
           subtitle = sprintf("Quantile clipping: %.0f%%-%.0f%%",
                              LOWER_QUANTILE * 100, UPPER_QUANTILE * 100),
           x = "Longitude", y = "Latitude") +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        panel.grid.major = element_line(color = "gray80", size = 0.3),
        plot.margin = margin(2, 2, 2, 2),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      ) +
      annotation_scale(location = "bl", width_hint = 0.3) +
      annotation_north_arrow(location = "tl", which_north = "true",
                             style = north_arrow_fancy_orienteering)
  }

  plot_list_eq[[varname]] <- p

  # Save plot
  output_file <- file.path("./tests/sat/s2_l2b_mp/", paste0("inv_shallow_", varname, ".png"))
  ggsave(output_file, p, units = "in", width = 8, height = 6, dpi = 300)
}

cat("\nAll plots saved to ./tests/sat/s2_l2b_mp \n")

## SAVE RESULTS TO NETCDF ----

cat("\nSaving results to NetCDF...\n")

output_nc <- sprintf("./tests/sat/s2_l2b_mp/%s_inversion_shallow_%s.nc",
                    s2_img_name, INVERSION_APPROACH)

# Use ncdf4 package to write each variable separately
library(ncdf4)

# Get dimensions from the raster
nx <- ncol(r_stack)
ny <- nrow(r_stack)
ext_vals <- ext(r_stack)

# Create dimension variables
xvals <- seq(ext_vals[1], ext_vals[2], length.out = nx)
yvals <- seq(ext_vals[3], ext_vals[4], length.out = ny)

xdim <- ncdim_def("x", "meters", xvals)
ydim <- ncdim_def("y", "meters", yvals)

# Create a list to store all variable definitions
var_list <- list()

cat(sprintf("Creating NetCDF variable definitions for %d layers:\n", nlyr(r_stack)))

for (i in 1:nlyr(r_stack)) {
  var_name <- names(r_stack)[i]
  cat(sprintf("  %d. %s\n", i, var_name))

  var_list[[i]] <- ncvar_def(
    name = var_name,
    units = "",
    dim = list(xdim, ydim),
    missval = -9999,
    longname = var_name
  )
}

# Create lat/lon coordinate variables (2D grids)
cat("\nCreating lat/lon coordinate variables...\n")

# Create grid of x,y coordinates in projected CRS
xy_grid <- expand.grid(x = xvals, y = yvals)
xy_matrix <- as.matrix(xy_grid)

# Project to WGS84 lat/lon
cat("  Projecting coordinates from UTM to WGS84...\n")
latlon_coords <- project(xy_matrix, from = crs(r_stack), to = "EPSG:4326")

# Reshape to 2D grids matching raster matrix format
# Matrix should be (ny x nx) with rows=y, cols=x
lon_grid <- matrix(latlon_coords[, 1], nrow = ny, ncol = nx, byrow = FALSE)
lat_grid <- matrix(latlon_coords[, 2], nrow = ny, ncol = nx, byrow = FALSE)

cat(sprintf("  Lon range: %.6f to %.6f\n", min(lon_grid), max(lon_grid)))
cat(sprintf("  Lat range: %.6f to %.6f\n", min(lat_grid), max(lat_grid)))

# Define lat/lon as 2D coordinate variables
lon_var <- ncvar_def(
  name = "lon",
  units = "degrees_east",
  dim = list(xdim, ydim),
  missval = -9999,
  longname = "Longitude (WGS84)",
  prec = "double"
)

lat_var <- ncvar_def(
  name = "lat",
  units = "degrees_north",
  dim = list(xdim, ydim),
  missval = -9999,
  longname = "Latitude (WGS84)",
  prec = "double"
)

# Add lat/lon to variable list
var_list[[length(var_list) + 1]] <- lon_var
var_list[[length(var_list) + 1]] <- lat_var

cat(sprintf("  Total variables (including lat/lon): %d\n", length(var_list)))

# Create the NetCDF file
ncout <- nc_create(output_nc, var_list, force_v4 = TRUE)

# Write data for each variable
cat("\nWriting data to NetCDF file:\n")
for (i in 1:nlyr(r_stack)) {
  var_name <- names(r_stack)[i]
  cat(sprintf("  Writing %s...\n", var_name))

  # Extract values and convert to matrix
  vals <- values(r_stack[[i]], mat = TRUE)
  vals[is.na(vals)] <- -9999

  # Write to NetCDF
  ncvar_put(ncout, var_list[[i]], vals)
}

# Write lat/lon coordinate grids
cat("  Writing lon coordinate grid...\n")
ncvar_put(ncout, lon_var, lon_grid)

cat("  Writing lat coordinate grid...\n")
ncvar_put(ncout, lat_var, lat_grid)

# Add global attributes
ncatt_put(ncout, 0, "title", "SABER Shallow Water Inversion Results")
ncatt_put(ncout, 0, "source", "Sentinel-2 L2W ACOLITE")
ncatt_put(ncout, 0, "date_created", as.character(Sys.time()))
ncatt_put(ncout, 0, "crs", as.character(crs(r_stack)))

# Add coordinate system attributes
ncatt_put(ncout, "x", "standard_name", "projection_x_coordinate")
ncatt_put(ncout, "x", "long_name", "x coordinate of projection")
ncatt_put(ncout, "x", "axis", "X")

ncatt_put(ncout, "y", "standard_name", "projection_y_coordinate")
ncatt_put(ncout, "y", "long_name", "y coordinate of projection")
ncatt_put(ncout, "y", "axis", "Y")

ncatt_put(ncout, "lon", "standard_name", "longitude")
ncatt_put(ncout, "lon", "long_name", "Longitude")
ncatt_put(ncout, "lon", "axis", "X")

ncatt_put(ncout, "lat", "standard_name", "latitude")
ncatt_put(ncout, "lat", "long_name", "Latitude")
ncatt_put(ncout, "lat", "axis", "Y")

# Close the file
nc_close(ncout)

cat(sprintf("\nResults saved to: %s\n", output_nc))
cat(sprintf("Total variables written: %d (including lat/lon coordinate grids)\n", nlyr(r_stack) + 2))


## SUMMARY STATISTICS ----
cat("\n========================================\n")
cat("SHALLOW WATER INVERSION SUMMARY\n")
cat("========================================\n")

# Select only variables that exist in final_results
vars_available <- intersect(vars_to_plot, names(final_results))

summary_stats <- final_results %>%
  dplyr::select(all_of(vars_available)) %>%
  summary()

print(summary_stats)

cat("\n========================================\n")
cat("PROCESSING COMPLETE\n")
cat("========================================\n")
cat(sprintf("Total execution time: %.2f minutes\n", exec_time_parall[3]/60))
cat(sprintf("Output NetCDF: %s\n", output_nc))
cat(sprintf("Output plots: ./tests/sat/s2_l2b_mp/inv_shallow_*.png\n"))
cat("========================================\n")
