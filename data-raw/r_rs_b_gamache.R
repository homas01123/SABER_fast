## code to prepare `r_rs_b_gamache` dataset goes here
library(readr)

r_rs_b_gamache <- read_csv(fs::path_package("SABER", "extdata", "r_rs_b_gamache.csv"))

usethis::use_data(r_rs_b_gamache, overwrite = TRUE)
