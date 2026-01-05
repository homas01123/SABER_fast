# Extract Solar and Viewing Angles from Sentinel-2 L1C XML Metadata
# These angles will be used in SABER inversion as par_fixed values

library(xml2)
library(terra)
library(dplyr)

# Paths
xml_file <- "./tests/sat/s2_l1c_mp/S2A_MSIL1C_20190828T153601_N0500_R111_T19UEQ_20230705T010212.SAFE/S2A_MSIL1C_20190828T153601_N0500_R111_T19UEQ_20230705T010212.SAFE/GRANULE/L1C_T19UEQ_A021844_20190828T153558/MTD_TL.xml"
l2w_file <- "./tests/sat/s2_l2a_mp/S2A_MSI_2019_08_28_15_39_27_T19UEQ_L2W.nc"

cat("\n=== EXTRACTING ANGLES FROM S2 L1C XML ===\n")

# Read XML
xml_data <- read_xml(xml_file)

# Extract mean sun angles (backup if grid extraction fails)
mean_sza <- xml_find_first(xml_data, ".//Mean_Sun_Angle/ZENITH_ANGLE") %>% xml_text() %>% as.numeric()
mean_saa <- xml_find_first(xml_data, ".//Mean_Sun_Angle/AZIMUTH_ANGLE") %>% xml_text() %>% as.numeric()

cat(sprintf("Mean Solar Zenith Angle: %.2f deg\n", mean_sza))
cat(sprintf("Mean Solar Azimuth Angle: %.2f deg\n", mean_saa))

# Extract sun angle grids
cat("\nExtracting Sun Angle Grids...\n")

# Sun Zenith Grid
sun_zenith_node <- xml_find_first(xml_data, ".//Sun_Angles_Grid/Zenith")
sza_col_step <- xml_find_first(sun_zenith_node, ".//COL_STEP") %>% xml_text() %>% as.numeric()
sza_row_step <- xml_find_first(sun_zenith_node, ".//ROW_STEP") %>% xml_text() %>% as.numeric()
sza_values_nodes <- xml_find_all(sun_zenith_node, ".//VALUES")
sza_values_text <- xml_text(sza_values_nodes)

# Parse grid values
sza_grid <- lapply(sza_values_text, function(row) {
  vals <- strsplit(row, " +")[[1]]
  vals <- vals[vals != ""]
  as.numeric(vals)
})
sza_matrix <- do.call(rbind, sza_grid)

cat(sprintf("  SZA grid dimensions: %d x %d\n", nrow(sza_matrix), ncol(sza_matrix)))
cat(sprintf("  SZA grid resolution: %d m x %d m\n", sza_col_step, sza_row_step))
cat(sprintf("  SZA range: %.2f - %.2f deg\n", min(sza_matrix, na.rm=TRUE), max(sza_matrix, na.rm=TRUE)))

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

cat(sprintf("  SAA grid dimensions: %d x %d\n", nrow(saa_matrix), ncol(saa_matrix)))
cat(sprintf("  SAA range: %.2f - %.2f deg\n", min(saa_matrix, na.rm=TRUE), max(saa_matrix, na.rm=TRUE)))

# Get tile geocoding info
tile_ulx <- xml_find_first(xml_data, ".//Geoposition[@resolution='10']/ULX") %>% xml_text() %>% as.numeric()
tile_uly <- xml_find_first(xml_data, ".//Geoposition[@resolution='10']/ULY") %>% xml_text() %>% as.numeric()
tile_xdim <- xml_find_first(xml_data, ".//Geoposition[@resolution='10']/XDIM") %>% xml_text() %>% as.numeric()
tile_ydim <- xml_find_first(xml_data, ".//Geoposition[@resolution='10']/YDIM") %>% xml_text() %>% as.numeric()
tile_nrows <- xml_find_first(xml_data, ".//Size[@resolution='10']/NROWS") %>% xml_text() %>% as.numeric()
tile_ncols <- xml_find_first(xml_data, ".//Size[@resolution='10']/NCOLS") %>% xml_text() %>% as.numeric()
tile_crs <- xml_find_first(xml_data, ".//HORIZONTAL_CS_CODE") %>% xml_text()

cat(sprintf("\nTile Info:\n"))
cat(sprintf("  CRS: %s\n", tile_crs))
cat(sprintf("  Upper-left corner: (%d, %d)\n", tile_ulx, tile_uly))
cat(sprintf("  Dimensions: %d x %d at 10m resolution\n", tile_nrows, tile_ncols))

# Create raster from sun angle grids
# Grid step is 5000m, starting from ULX, ULY
grid_ncols <- ncol(sza_matrix)
grid_nrows <- nrow(sza_matrix)

# Calculate extent for the angle grid
# Grid starts at tile UL corner and extends by grid_step intervals
# Note: Y-axis goes DOWN (tile_ydim is negative), so ymin < ymax
angle_xmin <- tile_ulx
angle_xmax <- tile_ulx + (grid_ncols) * sza_col_step
angle_ymax <- tile_uly
angle_ymin <- tile_uly - (grid_nrows) * sza_row_step  # Subtract because moving down

cat(sprintf("\nAngle grid extent: [%.0f, %.0f, %.0f, %.0f]\n", 
            angle_xmin, angle_xmax, angle_ymin, angle_ymax))

# Create rasters for sun angles
r_sza <- rast(sza_matrix, 
              extent = c(angle_xmin, angle_xmax, angle_ymin, angle_ymax),
              crs = tile_crs)
names(r_sza) <- "sza"

r_saa <- rast(saa_matrix,
              extent = c(angle_xmin, angle_xmax, angle_ymin, angle_ymax),
              crs = tile_crs)
names(r_saa) <- "saa"

cat("\nSun angle rasters created\n")

# Load L2W image to get target grid
cat("\nLoading L2W image for target grid...\n")
l2w_raster <- rast(l2w_file, subds = "Rrs_443")

cat(sprintf("L2W dimensions: %d x %d\n", nrow(l2w_raster), ncol(l2w_raster)))
cat(sprintf("L2W extent: [%.0f, %.0f, %.0f, %.0f]\n", 
            ext(l2w_raster)[1], ext(l2w_raster)[2], 
            ext(l2w_raster)[3], ext(l2w_raster)[4]))

# Resample sun angles to L2W grid (bilinear interpolation)
cat("\nResampling sun angles to L2W grid...\n")
r_sza_resampled <- resample(r_sza, l2w_raster, method = "bilinear")
r_saa_resampled <- resample(r_saa, l2w_raster, method = "bilinear")

cat(sprintf("Resampled SZA range: %.2f - %.2f deg\n", 
            min(values(r_sza_resampled), na.rm=TRUE), 
            max(values(r_sza_resampled), na.rm=TRUE)))

# Save as GeoTIFF for use in inversion
output_sza <- "./tests/sat/s2_l2a_mp/S2A_2019_08_28_sza.tif"
output_saa <- "./tests/sat/s2_l2a_mp/S2A_2019_08_28_saa.tif"

writeRaster(r_sza_resampled, output_sza, overwrite = TRUE)
writeRaster(r_saa_resampled, output_saa, overwrite = TRUE)

cat(sprintf("\nSaved:\n"))
cat(sprintf("  Solar Zenith: %s\n", output_sza))
cat(sprintf("  Solar Azimuth: %s\n", output_saa))

# Extract viewing angles (mean per band)
cat("\n=== EXTRACTING VIEWING ANGLES ===\n")
cat("Note: Viewing angles vary by detector and band\n")
cat("Using mean viewing incidence angles for each band:\n\n")

# Get mean viewing angles
mean_vza_nodes <- xml_find_all(xml_data, ".//Mean_Viewing_Incidence_Angle")

vza_by_band <- data.frame(
  bandId = character(),
  vza = numeric(),
  vaa = numeric(),
  stringsAsFactors = FALSE
)

for (node in mean_vza_nodes) {
  band_id <- xml_attr(node, "bandId")
  vza <- xml_find_first(node, ".//ZENITH_ANGLE") %>% xml_text() %>% as.numeric()
  vaa <- xml_find_first(node, ".//AZIMUTH_ANGLE") %>% xml_text() %>% as.numeric()
  
  vza_by_band <- rbind(vza_by_band, data.frame(
    bandId = band_id,
    vza = vza,
    vaa = vaa
  ))
}

print(vza_by_band)

# Calculate overall mean VZA (across all bands)
mean_vza <- mean(vza_by_band$vza, na.rm = TRUE)
mean_vaa <- mean(vza_by_band$vaa, na.rm = TRUE)

cat(sprintf("\nOverall Mean Viewing Zenith Angle: %.2f deg\n", mean_vza))
cat(sprintf("Overall Mean Viewing Azimuth Angle: %.2f deg\n", mean_vaa))

cat("\n=== SUMMARY ===\n")
cat(sprintf("For SABER inversion par_fixed:\n"))
cat(sprintf("  theta_sun: Use pixel-wise SZA from %s\n", basename(output_sza)))
cat(sprintf("  theta_view: Use mean VZA = %.2f deg (or 0 for nadir approximation)\n", mean_vza))
cat(sprintf("\nNote: S2 viewing angles are near-nadir (%.2f deg), so theta_view=0 is reasonable approximation\n", mean_vza))

cat("\n=== DONE ===\n")
