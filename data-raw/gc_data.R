## code to prepare `gc_data` dataset goes here
## Source: Gregg, W.W. & Carder, K.L. (1990). A simple spectral solar irradiance
##   model for cloudless maritime atmospheres. Limnology and Oceanography, 35(8), 1657-1675.
library(readr)

gc_data <- read_csv(fs::path_package("SABER", "extdata", "gc_data.csv"))

usethis::use_data(gc_data, overwrite = TRUE)
