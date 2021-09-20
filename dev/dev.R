

usethis::use_data(abiotic_stresses_data, version = 3, overwrite = T)



devtools::document()
devtools::install()

devtools::build_vignettes()
