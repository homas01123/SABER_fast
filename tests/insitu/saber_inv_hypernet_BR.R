library(ncdf4)
library(terra)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tidyverse)
library(scales)
#library(SABER)
library(parallel)
library(doParallel)
library(viridis)
library(ggh4x)
library(BayesianTools)

library(raster)
library(ggplot2)
library(rasterVis)
library(ggspatial)
library(sf)
library(cowplot)
library(progressr)

# Function :: Interpolate the input matrix to a given wavelength ----
interp_vectorized <- function(row, wavelength_obs, wavelength_out){

  obsdata_interp = approx(x= wavelength_obs, y = row, xout = wavelength_out, method = "linear")$y

  return(obsdata_interp)
}

# Function :: Read the in situ data for Berre lagoon ----
analyze_hypernet_insitu_data <- function(variable, date_range) {
  # Read the CSV file
  insitu_data <- read.csv("./tests/insitu/insitu_1994-2023.csv",
                          stringsAsFactors = FALSE, skip = 0)


  # Convert date column to Date format
  insitu_data$date <- as.Date(insitu_data$date, format = "%m/%d/%Y")

  # Convert date_range to Date format
  start_date <- as.Date(date_range[1], format = "%m/%d/%Y")
  end_date <- as.Date(date_range[2], format = "%m/%d/%Y")

  # Filter for surface observations only
  surface_data <- insitu_data[insitu_data$Profondeur == 0, ]

  # Filter for date range
  filtered_data <- surface_data[surface_data$date >= start_date & surface_data$date <= end_date, ]

  # Extract the variable of interest
  var_data <- filtered_data[[variable]]

  # Remove missing values
  var_data <- as.numeric(var_data[!is.na(var_data) & var_data != 999999])

  # Calculate summary statistics
  summary_stats <- list(
    filtered_data = filtered_data,
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

# Load Hypernet Rrs and ancillary data ----
hypernet = nc_open("C:/R/SABER/data/hypernet/BEFR_insitu.nc")
hypernet_db = nc_open("C:/R/SABER/data/hypernet/BEFR_insitu.nc")
insitu_data <- read.csv("./tests/insitu/insitu_1994-2023.csv",
                        stringsAsFactors = FALSE)

print(hypernet)
names(hypernet$var)
names(hypernet$dim)

hypernet_Rrs = ncvar_get(hypernet, "Rrs")

hypernet_sunzen = ncvar_get(hypernet, "solar_zenith_angle")
hypernet_vzen = ncvar_get(hypernet, "viewing_zenith_angle")

hypernet_acq_time = ncvar_get(hypernet, "acquisition_time")

acq_datetime <- as.POSIXct(hypernet_acq_time, origin = "1970-01-01", tz = "UTC") # Convert to POSIXct (seconds since 1970-01-01 UTC)
acq_date <- format(acq_datetime, "%Y-%m-%d") # Format in same style as in situ (YYYY-mm-dd)

hypernet_time_df <- tibble(
  acquisition_time = hypernet_acq_time,
  datetime = acq_datetime,
  date = as.Date(acq_date)
)

hypernet_wl = hypernet[["dim"]][["wl"]][["vals"]]
hypernet_flags = ncvar_get(hypernet, "quality_flag")

rownames(hypernet_Rrs) = hypernet_wl
colnames(hypernet_Rrs) = paste0("id_", seq(1, ncol(hypernet_Rrs), 1))
hypernet_Rrs = as.data.frame(t(hypernet_Rrs))  #Transpose to have wavelength in columns

matplot((hypernet_Rrs), col = viridis::viridis(500), type = "l", lwd=3)


wave_interp = c(412, 443, 465, 489, 510, 532, 555, 567, 589, 625, 667, 683,
                695, 710, 730, 747)

#wave_interp = seq(400, 750, 3.5) #hyperspectral wavelengths



hypernet_Rrs_ms = t(apply(hypernet_Rrs, 1, interp_vectorized, wavelength_obs = hypernet_wl,
                          wavelength_out = wave_interp))
colnames(hypernet_Rrs_ms) = wave_interp

hypernet_Rrs_ms = as.data.frame(hypernet_Rrs_ms)

# Convert rownames to a proper column
hypernet_long <- hypernet_Rrs_ms %>%
  rownames_to_column("ensemble") %>%
  pivot_longer(
    cols = -ensemble,
    names_to = "wavelength",
    values_to = "rrs_0p"
  ) %>%
  mutate(wavelength = as.numeric(wavelength),
        # rrs_0m = rrs_0p_to_0m(rrs_0p) # Apply conversion here
         rrs_0m = rrs_0p
  )

# Nest into SABER structure
hypernet_Rrs_ms_nested <- hypernet_long %>%
  dplyr::select(ensemble, wavelength, rrs_0m) %>%
  group_by(ensemble) %>%
  nest(data = c(wavelength, rrs_0m)) %>%
  ungroup()

hypernet_Rrs_ms_nested <- hypernet_Rrs_ms_nested %>%
  mutate(
    theta_sun = hypernet_sunzen,  # Add sun zenith angle
    datetime = hypernet_time_df$datetime  # Add acquisition datetime
  )

## Plot the hyperspectral HYPERNET Rrs ----
plot_data <- unnest(hypernet_Rrs_ms_nested, cols = c(data))
ggplot(plot_data, aes(x = wavelength, y = rrs_0m, color = theta_sun, group = ensemble)) +
  geom_line(alpha = 0.8) + # Use geom_line to connect the points for each spectrum
  scale_color_viridis_c(direction = -1) + # A nice color scale for continuous data
  labs(
    title = "Hyperspectral Rrs Colored by Solar Zenith Angle",
    x = "Wavelength (nm)",
    y = expression(paste("Remote Sensing Reflectance (sr"^{-1}, ")")), # Formats a nice y-axis label
    color = "Sun Zenith\nAngle (°)"
  ) +
  theme_minimal()

# Perform the inversion on hypernet Rrs----
## Get prior ranges from in situ data ----
chl_summary = analyze_hypernet_insitu_data(variable = names(insitu_data)[grep("chl", names(insitu_data),
                                                                              ignore.case = TRUE)],

                                           date_range = c(min(hypernet_time_df$date),
                                                          max(hypernet_time_df$date)))
## Prepare SABER inversion parameters ----

# Parameter bounds

use_sicf <- T  # Set to FALSE for elastic-only model
return_full_output <- F

site_lat <- 43.4423106
site_lon <- 5.0971775
date_time <- as.POSIXct("2020-06-15 10:30:00", tz = "UTC")

if (use_sicf) {
  cat("Using SICF model (am03_sicf) - includes fluorescence\n\n")

  # Parameters to retrieve (including phi_f!)
  par_inversed <- c(#"theta_sun",
    "chl", "a_g_440","a_g_s",
    #"a_g_s_d",
    "bb_p_550", "bb_p_gamma",

    # "a_nap_440",  "a_nap_sd",

    "phi_f",

    # "h_w",  "r_rs_b_saccharina" , "r_rs_b_cca" , "r_rs_b_gravel" ,

    "sd"
  )

  # Fixed parameters (use LIST to allow mixed types!)
  par_fixed <- list(
    water_type = 2,
    theta_sun = 30,
    theta_view = 40.008,
    #a_g_s_g = 0.014,
    #a_g_s_d = 0.003,
    #bb_p_gamma = 0.5,
    lat = site_lat,              # Required for SICF
    lon = site_lon,              # Required for SICF
    date_time = as.numeric(date_time),  # Required for SICF
    sicf_model = "semi_analytical",          # ← Character value OK in list
    depth_integration = FALSE            # ← Logical value OK in list
  )

  forward_model <- "am03_sicf"

  lower <- c(0.5, 0.1, 0.001, 0.003, 0.2, 0.01, 0.0001)
  init_val <- c(5, 0.5, 0.017, 0.005, 0.5, 0.02, 0.01)
  upper <- c(50, 2, 0.025, 0.015, 1, 0.03, 10)


} else {

  # CHOICE 2: Use elastic-only model (am03) - Faster, no fluorescence
  cat("Using elastic-only model (am03) - no fluorescence\n\n")

  # Parameters to retrieve (NO phi_f - cannot retrieve without SICF!)
  par_inversed <- c(#"theta_sun",
    "chl", "a_g_440","a_g_s",
    #"a_g_s_d",
    "bb_p_550", "bb_p_gamma",

    # "a_nap_440",  "a_nap_sd",

    # "h_w",  "r_rs_b_saccharina" , "r_rs_b_cca" , "r_rs_b_gravel" ,

    "sd"
  )

  # Fixed parameters (NO lat/lon/date_time needed)
  par_fixed <- list(
    water_type = 2,
    theta_sun = 30,
    theta_view = 40.008
    #, a_g_s_g = 0.014,
    #a_g_s_d = 0.003,
    #bb_p_gamma = 0.5
  )

  forward_model <- "am03"

  lower <- c(0.5, 0.1, 0.001, 0.003, 0.2, 0.0001)
  init_val <- c(5, 0.5, 0.017, 0.005, 0.5,0.01)
  upper <- c(50, 2, 0.025, 0.015, 1, 10)

}

## SABER inversion test with single spectra inversion ----
test_inverse_mcmc = inverse_mcmc(
  rrs = hypernet_Rrs_ms_nested$data[[10]],
  forward_model = forward_model,
  par_inversed = par_inversed,
  prior = NULL,
  lower = lower,
  best = NULL,
  upper = upper,
  par_fixed = par_fixed,
  iterations = 30000,
  burnin = 5000,
  sampler = "DEzs"
)
print(test_inverse_mcmc)


test_inverse_grad = inverse_gradient(
  rrs = hypernet_Rrs_ms_nested$data[[10]],
  forward_model = forward_model,
  objective_fct = "log-ll",
  #optim_mtd = "levenberg-marqardt",
  optim_mtd = "L-BFGS-B",
  par_inversed = par_inversed,
  par_fixed = par_fixed,
  lower_b  = lower,
  init_val = init_val,
  upper_b  = upper,
  verbose = T
)
print(test_inverse_grad)

## SABER inversion with parallel execution ----
gc()

num_cores <- parallel::detectCores() - 2
cl <- makeCluster(num_cores)
doParallel::registerDoParallel(cl)

# Partition into n roughly equal-sized chunks
chunked_nested <- split(hypernet_Rrs_ms_nested,
                        cut(seq_len(nrow(hypernet_Rrs_ms_nested)), breaks = num_cores, labels = FALSE))

# chunked_nested_test <- split(hypernet_Rrs_ms_nested[1:20,],
#                         cut(seq_len(nrow(hypernet_Rrs_ms_nested[1:20,])),
#                             breaks = num_cores, labels = FALSE))

### MCMC based inversion ----
exec_time_parall <- system.time({
  with_progress({
    p <- progressor(steps = length(chunked_nested))

    saber_results <- foreach(chunk = chunked_nested,
                             .packages = c("purrr", "SABER", "BayesianTools", "dplyr", "progressr"),
                             .export = c("inverse_mcmc", "par_inversed", "lower", "upper",
                                         "par_fixed", "forward_model")) %dopar% {

                             result <- chunk %>%
                               mutate(
                                 inversion_estim = purrr::pmap(
                                   list(data, theta_sun, date_time),
                                   ~ {
                                     # Create a local copy and update with observed values
                                     par_fixed_local <- par_fixed
                                     par_fixed_local["theta_sun"] <- ..2
                                     par_fixed_local["date_time"] <- as.numeric(..3)

                                     inverse_mcmc(
                                       rrs = ..1,
                                       forward_model = forward_model,
                                       par_inversed = par_inversed,
                                       prior = NULL,
                                       lower = lower,
                                       best = NULL,
                                       upper = upper,
                                       par_fixed = par_fixed_local,
                                       iterations = 22000,
                                       burnin = 5000,
                                       sampler = "DEzs"
                                     )
                                   }
                                 )
                               )
                             p()  # increment progress bar
                               result
                             }
  })
})

### Gradient based inversion ----
# Setup for inverse_gradient parallel execution
exec_time_parall <- system.time({
  with_progress({
    p <- progressor(steps = length(chunked_nested))

    saber_results <- foreach(chunk = chunked_nested,
                             .packages = c("purrr", "SABER", "dplyr", "progressr", "numDeriv", "MASS"),
                             .export = c("inverse_gradient", "par_inversed", "lower_b", "init_val",
                                         "upper_b", "par_fixed", "objective_factory", "forward_model",
                                         "parse_inverse_parameter", "myFun")) %dopar% {

                                         result <- chunk %>%
                                           mutate(
                                             inversion_estim = purrr::pmap(
                                               list(data, theta_sun, datetime),
                                               ~ {
                                                 # Create a local copy and update with observed values
                                                 par_fixed_local <- par_fixed
                                                 par_fixed_local["theta_sun"] <- ..2
                                                 par_fixed_local["date_time"] <- as.numeric(..3)

                                                 inverse_gradient(
                                                   rrs = ..1,
                                                   forward_model = forward_model,
                                                   objective_fct = "log-ll",
                                                   optim_mtd = "L-BFGS-B",
                                                   par_inversed = par_inversed,
                                                   par_fixed = par_fixed_local,
                                                   lower_b = lower_b,
                                                   init_val = init_val,
                                                   upper_b = upper_b,
                                                   verbose = FALSE
                                                 )
                                               }
                                             )
                                           )
                                          p()  # increment progress bar
                                           result
                                         }
  })
})

stopCluster(cl)

# Save the inversion results into memory and disc----
combined_results <- bind_rows(saber_results)
ensemble_results <- combined_results %>%
  dplyr::select(ensemble, inversion_estim)
ensemble_results_tidy <- ensemble_results %>%
  mutate(inversion_estim = map(inversion_estim, ~ as_tibble(t(.x)))) %>%  # transpose vector to a row tibble
  unnest(cols = c(inversion_estim))

write.csv(x = ensemble_results_tidy, file = "./sat/hypernet_inversion_results_rrs0p_ms.csv",
          quote = F, col.names = T, sep = ",")

# Load the OG SABER inversion on hypernetdata Rrs
fit_results_hs <- read.csv("C:/R/SABER/outputs/hypernet_inversion_results_old.csv", header = TRUE)
fit_results_hs$id <- seq(1, nrow(fit_results_hs), 1)

# Visualize the inversion retrieved parameters ----
plot_df <- bind_cols(ensemble_results_tidy, hypernet_time_df) #current inversion
#====================== OR ============================#
plot_df <- bind_cols(fit_results_hs, hypernet_time_df) #OG inversion

df_values <- plot_df %>%
  dplyr::select(datetime, chl, a_g_440, bb_p_550, phi_f) %>%
  pivot_longer(
    cols = -datetime,
    names_to = "Variable",
    values_to = "Value"
  )

df_sds <- plot_df %>%
  dplyr::select(datetime, chl_sd, a_g_440_sd, bb_p_550_sd, phi_f_sd) %>%
  rename( # Rename columns to remove "_sd" for easier pivoting
    chl = chl_sd,
    a_g_440 = a_g_440_sd,
    bb_p_550 = bb_p_550_sd,
    phi_f = phi_f_sd
  ) %>%
  pivot_longer(
    cols = -datetime,
    names_to = "Variable",
    values_to = "SD"
  )

plot_df_long <- left_join(df_values, df_sds, by = c("datetime", "Variable"))

## Timeseries plot of retrieved OSCs with [chl] validation ----
custom_labels <- c(
  "a_g_440" = expression(paste(italic(a)[dg](443), " [m"^{-1},"]")),
  "bb_p_550" = expression(paste(italic(b)[bp](555), " [m"^{-1},"]")),
  "chl" = expression(paste("[",italic("chl"),"]", " [mg m"^{-3},"]")),
  "phi_f" = expression(paste(italic(phi)[f]))
)

custom_colors <- c(
  "a_g_440" = "goldenrod",
  "bb_p_550" = "#365C8DFF",
  "chl" = "#4AC16DFF",
  "phi_f" = "purple"
)

insitu_chl_data_prepared <- chl_summary$filtered_data %>%
  #filter for the specific station
  dplyr::filter(nom_station == "H12") %>%

  dplyr::select(
    datetime = date,
    chl_value = `Chloa`
  ) %>%
  dplyr::mutate(

    datetime = as.POSIXct(datetime)
  ) %>%

  tidyr::drop_na(chl_value)

insitu_for_facet_plot <- insitu_chl_data_prepared %>%
  mutate(Variable = "chl")

timeseries_plot <- ggplot(plot_df_long, aes(x = datetime, y = Value, color = Variable)) +

  geom_ribbon(aes(ymin = Value - 1.96 * SD, ymax = Value + 1.96 * SD, fill = Variable),
              alpha = 0.3, linetype = "blank", show.legend = TRUE) +

  geom_line(linewidth = 0.5) +

  geom_point(data = insitu_for_facet_plot,
             aes(y = as.numeric(chl_value)),
             color = "darkred",
             size = 3,
             shape = 19) +

facet_wrap(~Variable, ncol = 1, scales = "free_y", strip.position = "left") +

  scale_color_manual(values = custom_colors, labels = custom_labels) +
  scale_fill_manual(values = custom_colors, labels = custom_labels) +

  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_blank(), # Hides the facet titles on top
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12, color = 'black', angle = 15, hjust = 1),
    axis.text.y = element_text(size = 12, color = 'black'),
    axis.title.x = element_text(size = 20, margin = margin(t = 15)),
    axis.title.y = element_blank(),
    axis.ticks.length = unit(.25, "cm"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.background = element_rect(fill = NA),
    legend.key = element_blank(),
    legend.key.width = unit(1.5, "cm"),
    panel.grid.major = element_line(color = "grey50", linewidth = 0.5, linetype = "dotted"),
    panel.grid.minor = element_line(color = "grey80", linewidth = 0.2),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  labs(x = "Date") # Add x-axis label


print(timeseries_plot)

ggsave(paste0("./hypernet_inv_param_inel_mcmc.png"),
       plot = timeseries_plot,
       scale = 1.25, width = 9, height = 6,
       units = "in",dpi = 300)


## Scatter plots of Rrs vs retrieved parameters ----
hypernet_Rrs_ms = as.data.frame(hypernet_Rrs_ms)
hypernet_Rrs_ms$id = rownames(hypernet_Rrs_ms)

fit_results_hs$id = hypernet_Rrs_ms$id
Rrs_fitparam_comb <- inner_join(fit_results_hs, hypernet_Rrs_ms, by = "id")


plot_rrs_vs_param <- function(data, xvar, yvar, xlim, ylim,
                              xlab_expr, ylab_expr,
                              point_color = "darkgreen",
                              density_bins = 6) {

  asp_rat <- (xlim[2] - xlim[1]) / (ylim[2] - ylim[1])
  xstp <- (xlim[2] - xlim[1]) / 5
  ystp <- (ylim[2] - ylim[1]) / 5

  sd_col <- paste0("sd_", yvar)

  g <- ggplot(data = data, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(color = point_color, alpha = 0.8, size = 3) +

    # Conditional ribbon for uncertainty
    {
      if (sd_col %in% names(data)) {
        geom_ribbon(aes(ymin = .data[[yvar]] - .data[[sd_col]],
                        ymax = .data[[yvar]] + .data[[sd_col]]),
                    fill = "black", alpha = 0.66, colour = NA)
      }
    } +

    geom_density_2d(bins = density_bins, linewidth = 0.25, size = 1.1, color = "black") +

    geom_smooth(size=1,level = 0.66,show.legend = F,linetype = "solid",
                color="grey",
                se= T, method = "loess")+

    geom_hline(yintercept = 0, color = "darkred", size = 1.3, linetype = "dashed") +

    coord_fixed(ratio = asp_rat, xlim = xlim, ylim = ylim, expand = FALSE, clip = "on") +

    scale_x_continuous(name = xlab_expr, limits = xlim, breaks = seq(xlim[1], xlim[2], xstp)) +
    scale_y_continuous(name = ylab_expr, limits = ylim, breaks = seq(ylim[1], ylim[2], ystp)) +

    theme_bw() +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 20, color = 'black'),
      axis.text.y = element_text(size = 20, color = 'black'),
      axis.title.x = element_text(size = 25),
      axis.title.y = element_text(size = 25),
      axis.ticks.length = unit(0.25, "cm"),
      legend.position = "none",
      panel.grid.major = element_line(colour = "black", size = 0.5, linetype = "dotted"),
      panel.grid.minor = element_line(colour = "grey80", linewidth = 0.2),
      panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
      plot.margin = unit(c(0.5, 1.0, 0.5, 0.5), "cm")
    )

  # Add marginal histogram
  g_out <- ggMarginal(p = g, type = "densigram", bins = 45, color = "grey", groupFill = FALSE)

  return(g_out)
}


plot_chl = plot_rrs_vs_param(
  data = Rrs_fitparam_comb,
  xvar = "695",
  yvar = "chl",
  xlim = c(0, 0.005),
  ylim = c(-5, 20),
  xlab_expr = expression(paste(italic("R")["rs"]("0"^"-", 683), " in situ [sr"^{-1}, "]")),
  ylab_expr = expression(paste("[", italic("chl"), "]", italic(" predicted"), " [mg m"^{-3}, "]")),
  point_color = "darkgreen",
  density_bins = 10
)

ggsave(paste0("./outputs/hypernet_rrs683_chl.png"), plot = plot_chl,
       scale = 1.25, width = 4.5, height = 4.5, units = "in",dpi = 300)

plot_adg443 = plot_rrs_vs_param(
  data = Rrs_fitparam_comb,
  xvar = "443",
  yvar = "adg443",
  xlim = c(0, 0.010),
  ylim = c(-0.5, 2),
  xlab_expr = expression(paste(italic("R")["rs"]("0"^"-", 443), " in situ [sr"^{-1}, "]")),
  ylab_expr = expression(paste(italic("a")["dg"](443),italic("predicted"), " [", "m"^-1, "]")),
  point_color = "goldenrod",
  density_bins = 10
)

ggsave(paste0("./outputs/hypernet_rrs443_adg443.png"), plot = plot_adg443,
       scale = 1.25, width = 4.5, height = 4.5, units = "in",dpi = 300)

plot_bbp555 = plot_rrs_vs_param(
  data = Rrs_fitparam_comb,
  xvar = "555",
  yvar = "bbp555",
  xlim = c(0, 0.020),
  ylim = c(0, 0.012),
  xlab_expr = expression(paste(italic("R")["rs"]("0"^"-", 555), " in situ [sr"^{-1}, "]")),
  ylab_expr = expression(paste(italic("b")["bp"](555),italic("predicted"), " [", "m"^-1, "]")),
  point_color = "navyblue",
  density_bins = 10
)

ggsave(paste0("./outputs/hypernet_rrs555_bbp555.png"), plot = plot_bbp555,
       scale = 1.25, width = 4.5, height = 4.5, units = "in",dpi = 300)



plot_dual_yaxis_fixed <- function(data,
                                  id_col = "id_num",        # numeric ID column
                                  param_col = "chl",        # e.g., chl, adg443, bbp555
                                  sd_col = "sd_chl",
                                  rrs_col = "683",          # e.g., Rrs wavelength
                                  y1_limits = c(-5, 20),
                                  y2_limits = c(0, 0.005),
                                  ylab_primary = expression(paste("[", italic("chl"), "]", italic(" predicted"), " [mg m"^{-3}, "]")),
                                  ylab_secondary = expression(paste(italic("R")["rs"]("0"^"-", 683), " in situ [sr"^{-1}, "]")),
                                  x_limits = c(0, 2000),
                                  x_breaks = 500,
                                  y1_color = "darkgreen",
                                  y2_color = "steelblue",
                                  show_primary_smooth = TRUE,
                                  show_secondary_smooth = TRUE) {

  # Scaling factor to align Rrs values with primary Y-axis
  scale_factor <- diff(y1_limits) / diff(y2_limits)

  # Create new variable in data for scaled Rrs
  data$rrs_scaled <- data[[rrs_col]] * scale_factor

  p <- ggplot(data, aes(x = .data[[id_col]])) +

    # Primary parameter uncertainty ribbon
    # geom_ribbon(
    #   aes(
    #     ymin = .data[[param_col]] - .data[[sd_col]],
    #     ymax = .data[[param_col]] + .data[[sd_col]]
    #   ),
    #   fill = y1_color, alpha = 0.25
    # ) +

    # Primary parameter line
    geom_point(aes(y = (.data[[param_col]])), color = y1_color, linewidth = 1.1, alpha = 0.5) +

    # Optional trendline for primary Y
    {if (show_primary_smooth) geom_smooth(aes(y = .data[[param_col]]),
                                          color = y1_color, linetype = "dashed", se = T, size=2.5)} +

    # Rrs line, scaled to primary Y-axis
    geom_point(aes(y = (rrs_scaled)), color = y2_color, linewidth = 1.1, alpha = 0.5) +

    # Optional trendline for secondary Y
    {if (show_secondary_smooth) geom_smooth(aes(y = rrs_scaled),
                                            color = y2_color, linetype = "dashed", se = T, size=2.5)} +

    scale_x_continuous(
      name = "Observation ID",
      limits = x_limits,
      breaks = seq(x_limits[1], x_limits[2], by = x_breaks)
    ) +

    # Primary Y axis + secondary Y axis
    scale_y_continuous(
      name = ylab_primary,
      limits = y1_limits,
      sec.axis = sec_axis(~ . / scale_factor,
                          name = ylab_secondary)
    ) +

    theme_bw() +
    theme(
      axis.title.y.left = element_text(color = y1_color, size = 14),
      axis.text.y.left = element_text(color = y1_color, size = 12),
      axis.title.y.right = element_text(color = y2_color, size = 14),
      axis.text.y.right = element_text(color = y2_color, size = 12),
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey70", linetype = "dotted"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
      plot.margin = unit(c(1, 1, 1, 1), "lines")
    )

  return(p)
}

chl_series = plot_dual_yaxis_fixed(
  data = Rrs_fitparam_comb,
  id_col = "X",
  param_col = "chl",
  sd_col = "sd_chl",
  rrs_col = "695",
  y1_limits = c(0, 20),
  y2_limits = c(0, 0.005), y2_color = "black",
  x_limits = c(0, 2000),
  x_breaks = 500
)

ggsave(paste0("./outputs/hypernet_rrs683_chl_series.png"), plot = chl_series,
       scale = 1.25, width = 6, height = 4.5, units = "in",dpi = 300)

adg443_series = plot_dual_yaxis_fixed(
  data = Rrs_fitparam_comb,
  id_col = "X",
  param_col = "adg443",
  sd_col = "sd_adg443",
  rrs_col = "443",
  y1_limits = c(0, 2),
  y2_limits = c(0, 0.01),
  x_limits = c(0, 2000),
  x_breaks = 500, y1_color = "goldenrod", y2_color = "black",
  ylab_primary = expression(paste(italic("a")["dg"](443),italic("predicted"), " [", "m"^-1, "]")),
  ylab_secondary = expression(paste(italic("R")["rs"]("0"^"-", 443), " in situ [sr"^{-1}, "]"))
)

ggsave(paste0("./outputs/hypernet_rrs443_adg443_series.png"), plot = adg443_series,
       scale = 1.25, width = 6, height = 4.5, units = "in",dpi = 300)

bbp555_series = plot_dual_yaxis_fixed(
  data = Rrs_fitparam_comb,
  id_col = "X",
  param_col = "bbp555",
  sd_col = "sd_bbp555",
  rrs_col = "555",
  y1_limits = c(0, 0.012),
  y2_limits = c(0, 0.020),
  x_limits = c(0, 2000),
  x_breaks = 500, y1_color = "darkblue", y2_color = "black",
  ylab_primary = expression(paste(italic("b")["bp"](555),italic("predicted"), " [", "m"^-1, "]")),
  ylab_secondary = expression(paste(italic("R")["rs"]("0"^"-", 555), " in situ [sr"^{-1}, "]"))
)

ggsave(paste0("./outputs/hypernet_rrs555_bbp555_series.png"), plot = bbp555_series,
       scale = 1.25, width = 6, height = 4.5, units = "in",dpi = 300)
