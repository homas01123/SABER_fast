## code to prepare `r_rs_b_egsl` dataset goes here
library(readr)
library(tidyr)
library(dplyr)

r_rs_b_egsl <- read_csv(fs::path_package("SABER", "extdata", "r_rs_b_egsl.csv")) %>%
  pivot_longer(-wavelength, names_to = "class", values_to = "r_rs_b_mean")

usethis::use_data(r_rs_b_egsl, overwrite = TRUE)
