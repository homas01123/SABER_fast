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

# ============================================================================
# BATCH PROCESSING SCRIPT FOR MULTIPLE SENTINEL-2 IMAGES
# ============================================================================
# This script processes multiple S2 images across a year range using the same
# inversion workflow as s2_run_saber.R, but automated for batch processing.
#
# INVERSION METHOD: Uses gradient-based optimization (L-BFGS-B) by default
#   - MUCH FASTER than MCMC (seconds vs minutes per image)
#   - Still provides parameter estimates and standard errors
#   - Recommended for batch processing
#
# To use MCMC instead, set INVERSION_METHOD <- "mcmc" in configuration below
# ============================================================================

# Function :: Convert raster to dataframe for ggplot ----
raster_to_df <- function(rast) {
  df <- as.data.frame(rast, xy = TRUE)
  names(df)[3] <- "value"
  return(df)
}

# Function :: Load the Sentinel L2A Rrs raster ----
get_s2_files_by_month <- function(base_path,
                                  months_of_interest,
                                  file_pattern = "\\.nc$",
                                  mid_day = "mid") {
  # 1. List all relevant files
  all_nc_files <- list.files(
    path = base_path,
    pattern = file_pattern,
    recursive = TRUE,
    full.names = TRUE
  )

  # 2. Split into main and ancillary
  main_nc_files <- all_nc_files[!grepl("_anc\\.nc$", all_nc_files)]
  anc_nc_files  <- all_nc_files[grepl("_anc\\.nc$", all_nc_files)]

  # 3. Output list
  selected_files <- list()

  # 4. Loop over requested months
  for (month_path in months_of_interest) {

    # Get all main files for this month
    month_main_files <- main_nc_files[grepl(month_path, main_nc_files)]

    if (length(month_main_files) == 0) {
      selected_files[[month_path]] <- list(ref = NA, anc = NA)
      next
    }

    # Extract acquisition dates (e.g., 20170306T)
    date_tags <- regmatches(
      month_main_files,
      regexpr("\\d{8}T", month_main_files)
    )

    days <- as.integer(substr(date_tags, 7, 8))  # Extract day

    # Sort by day
    ordered_idx <- order(days)
    month_main_files <- month_main_files[ordered_idx]
    days <- days[ordered_idx]

    # Decide selection strategy
    select_idx <- switch(mid_day,
                         "start" = 1,
                         "end"   = length(days),
                         "mid"   = which.min(abs(days - 15)),
                         as.integer(mid_day))  # or a specific numeric day

    # Handle edge cases
    if (select_idx < 1 || select_idx > length(month_main_files)) {
      selected_files[[month_path]] <- list(ref = NA, anc = NA)
      next
    }

    selected_main <- month_main_files[select_idx]

    # Try finding matching ancillary file
    anc_tag <- sub("\\.nc$", "_anc.nc", basename(selected_main))
    matching_anc <- anc_nc_files[basename(anc_nc_files) == anc_tag]

    selected_files[[month_path]] <- list(
      ref = selected_main,
      anc = ifelse(length(matching_anc) > 0, matching_anc, NA)
    )
  }

  return(selected_files)
}

# Function :: Get S2 images within a year range ----
get_s2_files_by_year_range <- function(base_path,
                                       year_start,
                                       year_end,
                                       file_pattern = "\\.nc$",
                                       mid_day = "mid") {
  # Generate all year/month combinations in the range
  years <- year_start:year_end
  months <- sprintf("%02d", 1:12)
  
  # Create all possible year/month paths
  year_month_paths <- expand.grid(year = years, month = months, stringsAsFactors = FALSE) %>%
    mutate(path = paste(year, month, sep = "/")) %>%
    pull(path)
  
  # Get files for all these paths
  selected_files <- get_s2_files_by_month(
    base_path = base_path,
    months_of_interest = year_month_paths,
    file_pattern = file_pattern,
    mid_day = mid_day
  )
  
  # Filter out months with no data (NA entries)
  selected_files <- selected_files[!sapply(selected_files, function(x) is.na(x$ref))]
  
  return(selected_files)
}

# Function :: Process a single S2 image ----
process_single_s2_image <- function(s2_boa_ref,
                                    s2_boa_anc,
                                    image_id,
                                    par_to_inverse,
                                    lower,
                                    best,
                                    upper,
                                    par_fixed_base,
                                    inversion_method = "gradient",
                                    iterations = 22000,
                                    burnin = 5000,
                                    sampler = "DEzs",
                                    output_dir = "./",
                                    save_plots = TRUE,
                                    save_netcdf = TRUE,
                                    cl = NULL,
                                    num_cores = NULL) {
  
  cat("\n========================================\n")
  cat("Processing image:", image_id, "\n")
  cat("========================================\n")
  
  # Load raster data
  s2_reflectance <- rast(s2_boa_ref)
  s2_anc <- rast(s2_boa_anc)
  
  # Get dimensions
  dims <- dim(s2_reflectance)
  nrows <- dims[1]
  ncols <- dims[2]
  nbands <- dims[3]
  
  # Extract spectral bands
  selected_band_names <- names(s2_reflectance)[1:6]
  spectral_raster <- s2_reflectance[[selected_band_names]]
  
  # Convert to dataframe
  df <- as.data.frame(spectral_raster, xy = TRUE, cells = TRUE, na.rm = TRUE) %>%
    mutate(
      row = rowFromCell(spectral_raster, cell),
      col = colFromCell(spectral_raster, cell),
      ensemble = paste0("id_", row_number())
    )
  
  # Pivot to long format and convert rrs_0p to rrs_0m
  df_long <- df %>%
    dplyr::select(ensemble, starts_with("Rrs_wl=")) %>%
    pivot_longer(
      cols = starts_with("Rrs_wl="),
      names_to = "wavelength",
      names_prefix = "Rrs_wl=",
      values_to = "rrs_0p"
    ) %>%
    mutate(
      wavelength = as.numeric(wavelength),
      rrs_0m = rrs_0p_to_0m(rrs_0p)
    )
  
  # Filter valid pixels
  valid_ensembles <- df_long %>%
    group_by(ensemble) %>%
    filter(all(rrs_0m >= 0, na.rm = TRUE)) %>%
    ungroup()
  
  # Nest data
  nested_df <- valid_ensembles %>%
    dplyr::select(ensemble, wavelength, rrs_0m) %>%
    group_by(ensemble) %>%
    nest(data = c(wavelength, rrs_0m)) %>%
    ungroup()
  
  cat("Valid pixels:", nrow(nested_df), "\n")
  
  # Update par_fixed with image-specific angles
  par_fixed <- par_fixed_base
  par_fixed["theta_view"] <- mean(values(s2_reflectance$vza), na.rm = TRUE)
  par_fixed["theta_sun"] <- mean(values(s2_reflectance$sza), na.rm = TRUE)
  
  # Set up parallel cluster if not provided
  cleanup_cluster <- FALSE
  if (is.null(cl)) {
    if (is.null(num_cores)) {
      num_cores <- parallel::detectCores() - 4
    }
    cat("Setting up parallel cluster with", num_cores, "cores...\n")
    cl <- makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
    cleanup_cluster <- TRUE
  }
  
  # Partition into chunks
  chunked_nested <- split(nested_df,
                          cut(seq_len(nrow(nested_df)), breaks = num_cores, labels = FALSE))
  
  # Progress handlers
  handlers("rstudio")
  
  gc()
  
  # Run inversion based on method
  if (inversion_method == "mcmc") {
    cat("Running MCMC inversion...\n")
    exec_time <- system.time({
      with_progress({
        p <- progressor(steps = length(chunked_nested))
        
        saber_results <- foreach(chunk = chunked_nested,
                                 .packages = c("purrr", "SABER", "BayesianTools", "dplyr", "progressr"),
                                 .export = c("inverse_mcmc", "par_to_inverse", "lower", "upper", "par_fixed")) %dopar% {
                                   
                                   result <- chunk %>%
                                     mutate(
                                       inversion_estim = purrr::map(data, ~ inverse_mcmc(
                                         rrs = .x,
                                         forward_model = "am03",
                                         par_inversed = par_to_inverse,
                                         prior = NULL,
                                         lower = lower,
                                         best = NULL,
                                         upper = upper,
                                         par_fixed = par_fixed,
                                         iterations = iterations,
                                         burnin = burnin,
                                         sampler = sampler
                                       ))
                                     )
                                   
                                   p()
                                   result
                                 }
      })
    })
  } else {  # gradient
    cat("Running gradient-based inversion...\n")
    exec_time <- system.time({
      with_progress({
        p <- progressor(steps = length(chunked_nested))
        
        saber_results <- foreach(chunk = chunked_nested,
                                 .packages = c("purrr", "SABER", "dplyr", "progressr", "numDeriv", "MASS"),
                                 .export = c("inverse_gradient", "par_to_inverse", "lower", "best",
                                             "upper", "par_fixed")) %dopar% {
                                               
                                               result <- chunk %>%
                                                 mutate(
                                                   inversion_estim = purrr::map(data, ~ inverse_gradient(
                                                     rrs = .x,
                                                     forward_model = "am03",
                                                     objective_fct = "log-ll",
                                                     optim_mtd = "L-BFGS-B",
                                                     par_inversed = par_to_inverse,
                                                     par_fixed = par_fixed,
                                                     lower_b = lower,
                                                     init_val = best,
                                                     upper_b = upper,
                                                     verbose = FALSE
                                                   ))
                                                 )
                                               
                                               p()
                                               result
                                             }
      })
    })
  }
  
  cat("Inversion completed in", exec_time[3], "seconds\n")
  
  # Clean up cluster if we created it
  if (cleanup_cluster) {
    stopCluster(cl)
  }
  
  # Combine results
  combined_results <- bind_rows(saber_results)
  
  ensemble_results <- combined_results %>%
    dplyr::select(ensemble, inversion_estim)
  
  ensemble_results_tidy <- ensemble_results %>%
    mutate(inversion_estim = map(inversion_estim, ~ as_tibble(t(.x)))) %>%
    unnest(cols = c(inversion_estim))
  
  final_results <- df %>%
    dplyr::select(ensemble, x, y, row, col) %>%
    distinct() %>%
    inner_join(ensemble_results_tidy, by = "ensemble")
  
  # Create rasters
  r_template <- rast(nrows = nrows, ncols = ncols,
                     ext = ext(spectral_raster),
                     crs = crs(spectral_raster))
  
  vars_to_raster <- setdiff(names(final_results), c("ensemble", "x", "y", "row",
                                                     "col", "cell", "sd", "sd_sd"))
  
  raster_list <- list()
  
  for (varname in vars_to_raster) {
    r <- rast(r_template)
    vals <- rep(NA_real_, ncell(r))
    cell_indices <- cellFromRowCol(r, final_results$row, final_results$col)
    vals[cell_indices] <- final_results[[varname]]
    values(r) <- vals
    names(r) <- varname
    raster_list[[varname]] <- r
  }
  
  r_stack <- rast(raster_list)
  names(r_stack) <- names(raster_list)
  
  # Save plots if requested
  if (save_plots && length(vars_to_raster) > 0) {
    cat("Generating plots...\n")
    plot_vars <- vars_to_raster[!vars_to_raster %in% c("a_g_s", "bb_p_gamma")]
    
    for (varname in plot_vars) {
      r <- raster_list[[varname]]
      df_plot <- raster_to_df(r)
      
      if (nrow(df_plot) > 0 && !all(is.na(df_plot$value))) {
        ecdf_vals <- ecdf(df_plot$value)
        df_plot$value_eq <- ecdf_vals(df_plot$value)
        
        breaks <- seq(from = min(df_plot$value, na.rm = TRUE), 
                      to = max(df_plot$value, na.rm = TRUE), length.out = 5)
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
          coord_fixed() +
          labs(x = "Longitude", y = "Latitude") +
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
        
        ggsave(paste0(output_dir, "/inv_plot_", gsub("/", "-", image_id), "_", varname, ".png"), p,
               units = "in", width = 6, height = 4.5, dpi = 300)
      }
    }
  }
  
  # Save NetCDF if requested
  if (save_netcdf) {
    nc_filename <- paste0(output_dir, "/", gsub("/", "-", image_id), "_inversion_results.nc")
    cat("\nSaving results to NetCDF...\n")
    
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
    ncout <- nc_create(nc_filename, var_list, force_v4 = TRUE)
    
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
    ncatt_put(ncout, 0, "title", "SABER Deep Water Inversion Results")
    ncatt_put(ncout, 0, "source", "Sentinel-2 L2A ACOLITE")
    ncatt_put(ncout, 0, "date_created", as.character(Sys.time()))
    ncatt_put(ncout, 0, "crs", as.character(crs(r_stack)))
    ncatt_put(ncout, 0, "image_id", image_id)
    
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
    
    cat(sprintf("\nResults saved to: %s\n", nc_filename))
    cat(sprintf("Total variables written: %d (including lat/lon coordinate grids)\n", nlyr(r_stack) + 2))
  }
  
  cat("Image processing complete!\n")
  
  return(list(
    image_id = image_id,
    results = final_results,
    rasters = r_stack,
    exec_time = exec_time[3],
    n_pixels = nrow(nested_df)
  ))
}

# Function :: Run inversion on multiple images ----
run_multi_image_inversion <- function(base_path,
                                      year_start,
                                      year_end,
                                      par_to_inverse,
                                      lower,
                                      best,
                                      upper,
                                      par_fixed_base,
                                      inversion_method = "gradient",
                                      iterations = 22000,
                                      burnin = 5000,
                                      sampler = "DEzs",
                                      output_dir = "./inversion_results",
                                      save_plots = TRUE,
                                      save_netcdf = TRUE,
                                      num_cores = NULL) {
  
  # Get all images in the year range
  cat("========================================\n")
  cat("BATCH PROCESSING: MULTIPLE S2 IMAGES\n")
  cat("========================================\n")
  cat("Searching for S2 images between", year_start, "and", year_end, "...\n")
  
  selected_files <- get_s2_files_by_year_range(
    base_path = base_path,
    year_start = year_start,
    year_end = year_end,
    mid_day = "mid"
  )
  
  cat("Found", length(selected_files), "images to process\n\n")
  
  if (length(selected_files) == 0) {
    stop("No images found in the specified year range")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n\n")
  }
  
  # Set up parallel cluster (persistent across all images)
  if (is.null(num_cores)) {
    num_cores <- parallel::detectCores() - 3
  }
  cat("Setting up persistent parallel cluster with", num_cores, "cores...\n")
  cl <- makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  
  # Process each image
  all_results <- list()
  successful_count <- 0
  failed_count <- 0
  
  for (i in seq_along(selected_files)) {
    image_id <- names(selected_files)[i]
    image_files <- selected_files[[i]]
    
    cat(sprintf("\n[%d/%d] Processing: %s\n", i, length(selected_files), image_id))
    
    tryCatch({
      result <- process_single_s2_image(
        s2_boa_ref = image_files$ref,
        s2_boa_anc = image_files$anc,
        image_id = image_id,
        par_to_inverse = par_to_inverse,
        lower = lower,
        best = best,
        upper = upper,
        par_fixed_base = par_fixed_base,
        inversion_method = inversion_method,
        iterations = iterations,
        burnin = burnin,
        sampler = sampler,
        output_dir = output_dir,
        save_plots = save_plots,
        save_netcdf = save_netcdf,
        cl = cl,  # Pass existing cluster
        num_cores = num_cores
      )
      
      all_results[[image_id]] <- result
      successful_count <- successful_count + 1
      
    }, error = function(e) {
      cat("ERROR processing image", image_id, ":", e$message, "\n")
      failed_count <- failed_count + 1
    })
    
    gc()  # Garbage collection between images
  }
  
  # Stop cluster
  stopCluster(cl)
  
  # Summary
  cat("\n========================================\n")
  cat("BATCH PROCESSING COMPLETE\n")
  cat("========================================\n")
  cat("Total images found    :", length(selected_files), "\n")
  cat("Successfully processed:", successful_count, "\n")
  cat("Failed                :", failed_count, "\n")
  cat("Output directory      :", output_dir, "\n")
  cat("========================================\n\n")
  
  if (successful_count > 0) {
    cat("Processing times per image:\n")
    for (image_id in names(all_results)) {
      cat(sprintf("  %-20s : %6.2f sec (%d pixels)\n", 
                  image_id, 
                  all_results[[image_id]]$exec_time,
                  all_results[[image_id]]$n_pixels))
    }
    cat("\n")
  }
  
  return(all_results)
}

# ============================================================================
# MAIN EXECUTION: Configure and run batch processing
# ============================================================================

# Configuration ----
BASE_PATH <- "./tests/sat/s2_l2a_berre_lagoon"
YEAR_START <- 2021
YEAR_END <- 2023  # Change to process multiple years (e.g., 2019)
OUTPUT_DIR <- "./tests/sat/s2_l2b_bl_inversion_results_batch"

# Inversion method (GRADIENT IS MUCH FASTER!)
INVERSION_METHOD <- "gradient"  # "gradient" (FAST - seconds per image) or "mcmc" (SLOW - minutes per image)

# MCMC parameters (only used if INVERSION_METHOD == "mcmc")
MCMC_ITERATIONS <- 22000
MCMC_BURNIN <- 5000
MCMC_SAMPLER <- "DEzs"

# Define fixed parameters (angles will be updated per image)
par_fixed_base <- c(
  "water_type" = 2,
  "theta_view" = NULL,  # Will be set per image
  "theta_sun" = NULL,   # Will be set per image
  "h_w" = NULL,
  "r_rs_b_Sand_2019" = NULL,
  "r_rs_b_Eelgrass_2019" = NULL,
  "r_rs_b_Mud_2019" = NULL
)

# Define parameters to inverse
par_to_inverse <- c(
  "chl", "a_g_440", "a_g_s", "bb_p_550", "bb_p_gamma", "sd"
)

# Parameter bounds
lower <- c(0.5, 0.1, 0.001, 0.003, 0.2, 0.0001)
best <- c(5, 0.5, 0.017, 0.005, 0.5, 0.01)
upper <- c(50, 2, 0.025, 0.015, 1, 10)

# Run batch processing ----
cat("Starting batch inversion processing...\n")
cat("Inversion method:", INVERSION_METHOD, "\n\n")

all_results <- run_multi_image_inversion(
  base_path = BASE_PATH,
  year_start = YEAR_START,
  year_end = YEAR_END,
  par_to_inverse = par_to_inverse,
  lower = lower,
  best = best,
  upper = upper,
  par_fixed_base = par_fixed_base,
  inversion_method = INVERSION_METHOD,
  iterations = MCMC_ITERATIONS,
  burnin = MCMC_BURNIN,
  sampler = MCMC_SAMPLER,
  output_dir = OUTPUT_DIR,
  save_plots = TRUE,
  save_netcdf = TRUE,
  num_cores = NULL  # NULL = auto-detect (num_cores - 4)
)

# Optional: Save all results to RData file ----
save(all_results, file = paste0(OUTPUT_DIR, "/batch_results_", YEAR_START, "-", YEAR_END, ".RData"))
cat("All results saved to:", paste0(OUTPUT_DIR, "/batch_results_", YEAR_START, "-", YEAR_END, ".RData\n"))

cat("\nBatch processing script complete!\n")
