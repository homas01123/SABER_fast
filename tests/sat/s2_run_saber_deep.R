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

# Load S2 data into work space ----
image_idx = "2017/04"#, "2017/08", "2017/12"
s2_selected <- get_s2_files_by_month(
  base_path = "./sat/s2_l2a_berre_lagoon",
  months_of_interest = image_idx,
  file_pattern = "\\.nc$",
  mid_day = "mid"  # Can be "start", "end", "mid", or specific number
)

s2_boa_ref = s2_selected[[1]]$ref
s2_boa_anc = s2_selected[[1]]$anc

s2_reflectance <- rast(s2_boa_ref)  # spectral data
s2_anc <- rast(s2_boa_anc)          # ancillary data

plot(s2_reflectance)
plot(s2_anc)

dims <- dim(s2_reflectance)
nrows <- dims[1]
ncols <- dims[2]
nbands <- dims[3]

selected_band_names <- names(s2_reflectance)[1:6]
spectral_raster <- s2_reflectance[[selected_band_names]]


df <- as.data.frame(spectral_raster, xy = TRUE, cells = TRUE, na.rm = TRUE) %>%
  mutate(
    row = rowFromCell(spectral_raster, cell),
    col = colFromCell(spectral_raster, cell),
    ensemble = paste0("id_", row_number())
  )


# Pivot all spectral bands to long format
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
    rrs_0m = rrs_0p_to_0m(rrs_0p)  # Apply conversion here
  )

# Filter out invalid pixels (e.g., all bands < 0)
valid_ensembles <- df_long %>%
  group_by(ensemble) %>%
  filter(all(rrs_0m >= 0, na.rm = TRUE)) %>%
  ungroup()

# Now nest into the desired format
nested_df <- valid_ensembles %>%
  dplyr::select(ensemble, wavelength, rrs_0m) %>%
  group_by(ensemble) %>%
  nest(data = c(wavelength, rrs_0m)) %>%
  ungroup()


# # Subset the nested dataframe to a smaller size for testing
# nested_df_subset = nested_df[1:22,]

# Invert S2 image reflectance using SABER ----

## Prepare SABER parameters ----

### Read the in situ data for Berre lagoon ----

analyze_variable <- function(variable, date_range) {
  # Read the CSV file
  insitu_data <- read.csv("./sat/BDD_HYDRO_GIPREB_compil_1994-2023.csv", stringsAsFactors = FALSE)

  # Convert date column to Date format
  insitu_data$date <- as.Date(insitu_data$date, format = "%m/%d/%Y")

  # Convert date_range to Date format
  start_date <- as.Date(date_range[1], format = "%m/%d/%Y")
  end_date <- as.Date(date_range[2], format = "%m/%d/%Y")

  # Filter for surface observations only
  surface_data <- insitu_data[insitu_data$Profondeur_m == "S", ]

  # Filter for date range
  filtered_data <- surface_data[surface_data$date >= start_date & surface_data$date <= end_date, ]

  # Extract the variable of interest
  var_data <- filtered_data[[variable]]

  # Remove missing values
  var_data <- var_data[!is.na(var_data) & var_data != 999999]

  # Calculate summary statistics
  summary_stats <- list(
    variable = variable,
    date_range = paste(date_range[1], "to", date_range[2]),
    count = length(var_data),
    mean = mean(var_data, na.rm = TRUE),
    median = median(var_data, na.rm = TRUE),
    min = min(var_data, na.rm = TRUE),
    max = max(var_data, na.rm = TRUE),
    sd = sd(var_data, na.rm = TRUE)
  )

  return(summary_stats)
}

par_fixed <- c(
  "water_type" = 2,
  "theta_view" = mean(values(s2_reflectance$vza), na.rm = TRUE),
  "theta_sun" = mean(values(s2_reflectance$sza), na.rm = TRUE)
  ,"h_w" = NULL,
  "r_rs_b_saccharina" = NULL,
  "r_rs_b_cca" = NULL,
  "r_rs_b_gravel" = NULL
  #, sd = 0.5
)


par_to_inverse = c(#"theta_sun",
  "chl", "a_g_440", "a_g_s", "bb_p_550", "bb_p_gamma",

  # "h_w",  "r_rs_b_saccharina" , "r_rs_b_cca" , "r_rs_b_gravel" ,

  "sd"
)

lower = c(0.5, 0.1, 0.001, 0.003, 0.2,
          # 0.5, 0.1, 0.1, 0.1,
          0.0001
)

best = c(10, 0.75, 0.017, 0.005, 0.46,
         #  4, 0.5, 0.1, 0.25,
         0.1
)

upper = c(30, 2, 0.025, 0.01, 1,
          #  12, 1, 1, 1,
          10
)

## Parallel execution ----
num_cores <- parallel::detectCores() - 2
cl <- makeCluster(num_cores)
doParallel::registerDoParallel(cl)

# Partition into n roughly equal-sized chunks
chunked_nested <- split(nested_df,
                        cut(seq_len(nrow(nested_df)), breaks = num_cores, labels = FALSE))



handlers("rstudio")  # Best for RStudio console

gc()

### MCMC Inversion ----
exec_time_parall <- system.time({
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
                                     iterations = 22000,
                                     burnin = 5000,
                                     sampler = "DEzs"
                                   ))
                                 )

                               p()  # increment progress bar
                               result
                             }
  })
})

#### Graident based inversion ----
exec_time_parall <- system.time({
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
                                     objective_fct = "log-ll",  # or "SSR", "lee99"
                                     optim_mtd = "L-BFGS-B",    # or "levenberg-marqardt", "auglag", etc.
                                     par_inversed = par_to_inverse,
                                     par_fixed = par_fixed,
                                     lower_b = lower,
                                     init_val = best,
                                     upper_b = upper,
                                     verbose = FALSE
                                   ))
                                 )

                               p()  # increment progress bar
                               result
                             }
  })
})

stopCluster(cl)

# Save the inversion results into memory ----

combined_results <- bind_rows(saber_results)


ensemble_results <- combined_results %>%
  dplyr::select(ensemble, inversion_estim)


ensemble_results_tidy <- ensemble_results %>%
  mutate(inversion_estim = map(inversion_estim, ~ as_tibble(t(.x)))) %>%  # transpose vector to a row tibble
  unnest(cols = c(inversion_estim))



final_results <- df %>%
  dplyr::select(ensemble, x, y, row, col) %>%
  distinct() %>%
  inner_join(ensemble_results_tidy, by = "ensemble")


# Prepare empty raster template with same dims and extent as original
r_template <- rast(nrows = nrows, ncols = ncols,
                   ext = ext(spectral_raster),
                   crs = crs(spectral_raster))

# List of variables to create rasters for (exclude sd and sd_sd)
vars_to_raster <- setdiff(names(final_results), c("ensemble", "x", "y", "row",
                                                  "col", "cell", "sd", "sd_sd"))

# Create a list to store rasters
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

#Test plotting the rasters (No visual enhancement)
for (varname in vars_to_raster) {
  print(paste("Plotting", varname))
  plot(raster_list[[varname]], main = varname, col = terrain.colors(50))
}


# Stack all layers into a SpatRaster
r_stack <- rast(raster_list)
names(r_stack) <- names(raster_list)

## Plot the SABER retrieved results as spatial rasters ----
plot_list_eq <- list()

for (varname in vars_to_raster[-c(3,5,8,10)]) {
  cat("Plotting (contrast-enhanced) ", varname, "\n")

  r <- raster_list[[varname]]
  df <- raster_to_df(r)

  # Equalized color mapping
  ecdf_vals <- ecdf(df$value)
  df$value_eq <- ecdf_vals(df$value)  # this spreads out color contrast

  # Keep track of original value range for legend
  breaks <- seq(from = min(df$value), to = max(df$value), length.out = 5)
  labels <- format(breaks, digits = 3)

  # Build color palette mapped to equalized values, but labeled with actual ones
  p <- ggplot(df, aes(x = x, y = y, fill = value_eq)) +
    geom_raster() +
    scale_fill_gradientn(
      colors = viridis::viridis(length(breaks)),
      breaks = #ecdf_vals(breaks),
        c(0,0.25,0.50,0.75,1),       # break positions in equalized space
      labels = labels,                  # labels in original value space
      name = varname,
      na.value = "gray90"
    ) +
    coord_equal() +
    labs(title = paste(varname), x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major = element_line(color = "gray80", size = 0.3),
      plot.margin = margin(2, 2, 2, 2),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )+
    annotation_scale(location = "bl", width_hint = 0.3) +
    annotation_north_arrow(location = "tl", which_north = "true",
                           style = north_arrow_fancy_orienteering)

  plot_list_eq[[varname]] <- p
  print(p)
  ggsave(paste0("./inv_plot_", gsub("/", "-", image_idx),"_",varname,".png"), p,
         units = "in",
         width = 6, height = 4.5, dpi = 300)
}


# Write the inversion results to  NetCDF file ----
writeCDF(
  r_stack,
  filename = "C:/Users/muks0001/Downloads/s2_l2a_berre_lagoon/2017/08/18/S2A_MSIL2A_20170818_inversion_results.nc",
  varname = names(r_stack),
  overwrite = TRUE
)






# inverse_mcmc(
#   rrs = nested_df_subset[[2]][[1]],
#   forward_model = "am03",
#   par_inversed = par_to_inverse,
#   prior = NULL,
#   lower = lower,
#   best = NULL,
#   upper = upper,
#   par_fixed = par_fixed,
#   iterations = 25000,
#   burnin = 5000,
#   sampler = "DEzs"
# )


# exec_time_parall <- system.time({
#   saber_results <- foreach(chunk = chunked_nested,
#                      .packages = c("purrr", "SABER", "BayesianTools", "dplyr"),
#                      .export = c("inverse_mcmc", "par_to_inverse", "lower", "upper", "par_fixed")) %dopar% {
#                        chunk %>%
#                          mutate(
#                            inversion_estim = purrr::map(data, ~ inverse_mcmc(
#                              rrs = .x,
#                              forward_model = "am03",
#                              par_inversed = par_to_inverse,
#                              prior = NULL,
#                              lower = lower,
#                              best = NULL,
#                              upper = upper,
#                              par_fixed = par_fixed,
#                              iterations = 25000,
#                              burnin = 5000,
#                              sampler = "DEzs"
#                            ))
#                          )
#                      }
# })
