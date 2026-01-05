library(terra)
library(ggplot2)
library(viridis)
library(scales)
library(ggspatial)
library(dplyr)
library(tidyr)

# ============================================================================
# TIME SERIES VISUALIZATION OF CHLOROPHYLL FROM INVERSION RESULTS
# ============================================================================

# ============================================================================
# USER-DEFINED PARAMETERS
# ============================================================================
LOWER_QUANTILE <- 0.01  # Lower quantile clip (e.g., 0.01 = 1st percentile)
UPPER_QUANTILE <- 0.99   # Upper quantile clip (e.g., 0.95 = 95th percentile)

# Variables to plot
# NetCDF layer naming can be either:
# 1. Named format: "chl", "a_g_440", "bb_p_550", "chl_sd", "a_g_440_sd", "bb_p_550_sd"
# 2. Numbered format: "chl_1", "chl_2", "chl_4", "chl_6", "chl_7", "chl_9"
#    where chl_1=chl, chl_2=a_g_440, chl_4=bb_p_550
#    and   chl_6=chl_sd, chl_7=a_g_440_sd, chl_9=bb_p_550_sd
VARIABLES <- list(
  list(name = "chl", named_layer = "chl", numbered_layer = "chl_1", 
       sd_named = "chl_sd", sd_numbered = "chl_6", 
       label = "Chl-a", units = "mg/m³"),
  list(name = "a_g_440", named_layer = "a_g_440", numbered_layer = "chl_2",
       sd_named = "a_g_440_sd", sd_numbered = "chl_7",
       label = "a_g(440)", units = "m⁻¹"),
  list(name = "bb_p_550", named_layer = "bb_p_550", numbered_layer = "chl_4",
       sd_named = "bb_p_550_sd", sd_numbered = "chl_9",
       label = "bb_p(550)", units = "m⁻¹")
)
# ============================================================================

# Directory containing NetCDF files
nc_dir <- "./tests/sat/s2_l2b_bl_inversion_results_batch"

# Get all NetCDF files
nc_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
nc_files <- sort(nc_files)

cat("Found", length(nc_files), "NetCDF files\n")

# Function to extract date from filename
extract_date <- function(filename) {
  # Extract YYYY-MM from filename like "2017-01_inversion_results.nc"
  date_str <- basename(filename)
  date_str <- sub("_inversion_results\\.nc$", "", date_str)
  return(date_str)
}

# Function to convert raster to dataframe
raster_to_df <- function(rast) {
  df <- as.data.frame(rast, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "value"
  return(df)
}

# Load all variable rasters
cat("\nLoading data from NetCDF files...\n")

# Initialize lists for each variable (include containers for sd_*)
data_lists <- list()
for (var in VARIABLES) {
  data_lists[[var$name]] <- list()
  # also create container for sd_* layers
  data_lists[[paste0(var$name, "_sd")]] <- list()
}
dates <- character()

for (i in seq_along(nc_files)) {
  file <- nc_files[i]
  date <- extract_date(file)
  
  cat(sprintf("  [%d/%d] Loading %s...\n", i, length(nc_files), date))
  
  # Load the NetCDF file
  r <- rast(file)
  
  # Check layer naming convention (first file may have named layers, rest numbered)
  layer_names <- names(r)

  # Extract each variable using adaptive layer detection
  for (var in VARIABLES) {
    # Try named layer first, fall back to numbered
    if (var$named_layer %in% layer_names) {
      var_layer <- var$named_layer
      sd_layer <- var$sd_named
      convention <- "named"
    } else if (var$numbered_layer %in% layer_names) {
      var_layer <- var$numbered_layer
      sd_layer <- var$sd_numbered
      convention <- "numbered"
    } else {
      cat(sprintf("    Warning: %s layer not found (tried %s and %s)\n", 
                  var$label, var$named_layer, var$numbered_layer))
      next
    }
    
    # Extract variable raster
    if (var_layer %in% layer_names) {
      var_raster <- r[[var_layer]]
      data_lists[[var$name]][[date]] <- var_raster
      cat(sprintf("    %s: Extracted from %s (%s format)\n", var$label, var_layer, convention))
    }
    
    # Extract SD layer if available
    if (sd_layer %in% layer_names) {
      sd_raster <- r[[sd_layer]]
      data_lists[[paste0(var$name, "_sd")]][[date]] <- sd_raster
      cat(sprintf("    %s SD: Extracted from %s\n", var$label, sd_layer))
    }
  }
  
  if (i == 1) {
    dates <- date
  } else {
    dates <- c(dates, date)
  }
}

cat("\nSuccessfully loaded data for", length(dates), "time points\n")

# Calculate global statistics with user-defined quantile clipping for each variable
global_limits <- list()

for (var in VARIABLES) {
  # Collect values across time points; some variables may be missing
  vals_list <- data_lists[[var$name]]
  if (length(vals_list) == 0 || all(sapply(vals_list, is.null))) {
    cat(sprintf("\nNo data found for %s - skipping global statistic computation\n", var$label))
    global_limits[[var$name]] <- list(min = NA, max = NA, median = NA, actual_min = NA, actual_max = NA)
    next
  }

  all_values <- unlist(lapply(vals_list, function(r) values(r, na.rm = TRUE)))
  all_values <- all_values[is.finite(all_values) & all_values > 0]
  
  if (length(all_values) == 0) {
    cat(sprintf("\nNo finite positive values found for %s - skipping\n", var$label))
    global_limits[[var$name]] <- list(min = NA, max = NA, median = NA, actual_min = NA, actual_max = NA)
    next
  }

  # Use quantiles for color scale limits
  global_min <- quantile(all_values, LOWER_QUANTILE, na.rm = TRUE)
  global_max <- quantile(all_values, UPPER_QUANTILE, na.rm = TRUE)
  global_median <- median(all_values, na.rm = TRUE)
  
  # Actual min/max for reference
  actual_min <- min(all_values, na.rm = TRUE)
  actual_max <- max(all_values, na.rm = TRUE)
  
  global_limits[[var$name]] <- list(
    min = global_min,
    max = global_max,
    median = global_median,
    actual_min = actual_min,
    actual_max = actual_max
  )
  
  cat(sprintf("\n%s statistics (%s):\n", var$label, var$units))
  cat(sprintf("  Actual Min: %.3e\n", actual_min))
  cat(sprintf("  Actual Max: %.3e\n", actual_max))
  cat(sprintf("  Median: %.3e\n", global_median))
  cat(sprintf("  Color scale (%.0f%%-%.0f%% quantile): %.3e - %.3e\n", 
              LOWER_QUANTILE * 100, UPPER_QUANTILE * 100, global_min, global_max))
}

# Create plots for each variable
cat("\nGenerating time series plots...\n")

library(cowplot)

for (var in VARIABLES) {
  cat(sprintf("\n=== Processing %s ===\n", var$label))
  
  plot_list <- list()
  var_data <- data_lists[[var$name]]
  var_limits <- global_limits[[var$name]]
  
  for (i in seq_along(dates)) {
    date <- dates[i]
    r <- var_data[[date]]
    
    if (is.null(r)) {
      cat(sprintf("  [%d/%d] Skipping %s (no data)...\n", i, length(dates), date))
      next
    }
    
    cat(sprintf("  [%d/%d] Creating plot for %s...\n", i, length(dates), date))
    
    # Convert to dataframe
    df <- raster_to_df(r)
    
    # Remove zeros and negative values
    df <- df[df$value > 0, ]
    
    # Clip values to quantile range
    df$value_clipped <- pmin(pmax(df$value, var_limits$min), var_limits$max)
    
    # Format date as Month/Year (e.g., "January 2017" or "Jan 2017")
    date_formatted <- format(as.Date(paste0(date, "-15")), "%B %Y")
    
    # Create plot with viridis color scale
    p <- ggplot(df, aes(x = x, y = y, fill = value_clipped)) +
      geom_raster() +
      scale_fill_viridis_c(
        option = "viridis",
        name = paste0(var$label, "\n(", var$units, ")"),
        limits = c(var_limits$min, var_limits$max),
        breaks = pretty(c(var_limits$min, var_limits$max), n = 6),
        na.value = "gray90",
        guide = "none"  # Remove legend
      ) +
      coord_fixed(ratio = 1) +  # Fixed aspect ratio instead of coord_equal
      labs(
        title = date_formatted,
        x = NULL,
        y = NULL
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "none",  # Remove legend
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(10, 10, 10, 10),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    plot_list[[date]] <- p
    
    # Save individual plot with consistent dimensions
    output_file <- file.path(nc_dir, paste0(var$name, "_timeseries_", date, ".png"))
    ggsave(output_file, p, width = 6, height = 6, dpi = 300, units = "in")
  }
  
  

  # If sd_* data are present, create SD plots (per-time and combined)
  sd_key <- paste0(var$name, "_sd")
  sd_data <- data_lists[[sd_key]]
  if (!is.null(sd_data) && length(sd_data) > 0 && !all(sapply(sd_data, is.null))) {
    cat(sprintf("  %s: found sd_* data, creating sd plots...\n", var$label))
    sd_plot_list <- list()
    sd_limits_vals <- unlist(lapply(sd_data, function(r) if (!is.null(r)) values(r, na.rm = TRUE) else NULL))
    sd_limits_vals <- sd_limits_vals[is.finite(sd_limits_vals)]
    if (length(sd_limits_vals) > 0) {
      sd_min <- quantile(sd_limits_vals, LOWER_QUANTILE, na.rm = TRUE)
      sd_max <- quantile(sd_limits_vals, UPPER_QUANTILE, na.rm = TRUE)
    } else {
      sd_min <- NA; sd_max <- NA
    }

    for (i in seq_along(dates)) {
      date <- dates[i]
      r_sd <- sd_data[[date]]
      if (is.null(r_sd)) next
      df_sd <- raster_to_df(r_sd)
      df_sd <- df_sd[is.finite(df_sd$value), ]
      df_sd$value_clipped <- if (!is.na(sd_min)) pmin(pmax(df_sd$value, sd_min), sd_max) else df_sd$value

      date_formatted <- format(as.Date(paste0(date, "-15")), "%B %Y")

      p_sd <- ggplot(df_sd, aes(x = x, y = y, fill = value_clipped)) +
        geom_raster() +
        scale_fill_viridis_c(option = "viridis", name = paste0("sd_", var$label, "\n(", var$units, ")"), 
                            limits = c(sd_min, sd_max), na.value = "gray90", guide = "none") +
        coord_fixed(ratio = 1) + 
        labs(title = paste0(date_formatted, " (SD)"), x = NULL, y = NULL) + 
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(10, 10, 10, 10),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
        )

      sd_plot_list[[date]] <- p_sd
      out_sd <- file.path(nc_dir, paste0(var$name, "_timeseries_", date, "_sd.png"))
      ggsave(out_sd, p_sd, width = 6, height = 6, dpi = 300, units = "in")
    }

    if (length(sd_plot_list) > 0) {
      combined_sd_file <- file.path(nc_dir, paste0(var$name, "_timeseries_combined_sd.png"))
      combined_sd_plot <- plot_grid(plotlist = sd_plot_list, ncol = 3, nrow = 4, align = "hv", axis = "tb")
      ggsave(combined_sd_file, combined_sd_plot, width = 18, height = 20, dpi = 300, units = "in")
      cat(sprintf("  %s sd plots saved!\n", var$label))
    }
  }
}

cat("\n========================================\n")
cat(sprintf("TIME SERIES VISUALIZATION COMPLETE\n"))
cat(sprintf("Quantile clipping: %.0f%% - %.0f%%\n", LOWER_QUANTILE * 100, UPPER_QUANTILE * 100))
cat("========================================\n")
cat("Individual plots saved as: VARIABLE_timeseries_YYYY-MM.png\n")
cat("Combined figures saved as: VARIABLE_timeseries_combined.png\n")
cat("Variables:", paste(sapply(VARIABLES, function(v) v$name), collapse = ", "), "\n")
cat("Location:", nc_dir, "\n")
cat("========================================\n\n")

# Create summary statistics tables for each variable
cat("Creating summary statistics tables...\n")

for (var in VARIABLES) {
  var_data <- data_lists[[var$name]]
  
  summary_df <- data.frame(
    Date = dates,
    Min = sapply(dates, function(d) {
      if (!is.null(var_data[[d]])) min(values(var_data[[d]], na.rm = TRUE)) else NA
    }),
    Mean = sapply(dates, function(d) {
      if (!is.null(var_data[[d]])) mean(values(var_data[[d]], na.rm = TRUE)) else NA
    }),
    Median = sapply(dates, function(d) {
      if (!is.null(var_data[[d]])) median(values(var_data[[d]], na.rm = TRUE)) else NA
    }),
    Max = sapply(dates, function(d) {
      if (!is.null(var_data[[d]])) max(values(var_data[[d]], na.rm = TRUE)) else NA
    }),
    SD = sapply(dates, function(d) {
      if (!is.null(var_data[[d]])) sd(values(var_data[[d]], na.rm = TRUE)) else NA
    })
  )
  
  cat(sprintf("\n%s statistics:\n", var$label))
  print(summary_df)
  
  # Save summary table
  csv_file <- file.path(nc_dir, paste0(var$name, "_timeseries_statistics.csv"))
  write.csv(summary_df, csv_file, row.names = FALSE)
  cat(sprintf("Saved to: %s\n", basename(csv_file)))
}

cat("\n*** To adjust color scale, edit LOWER_QUANTILE and UPPER_QUANTILE at the top of the script ***\n")
